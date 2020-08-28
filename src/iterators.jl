
"""
    AxisIterator(axis, window_size[, first_pad=nothing, last_pad=nothing,
                 stride=nothing, dilation=nothing])

Creates an iterator for indexing ranges of elements within `axis`.
"""
struct AxisIterator{B<:AbstractRange{<:Integer},W<:AbstractRange{<:Integer}}
    bounds::B
    window::W

    AxisIterator{B,W}(bounds::B, window::W) where {B,W} = new{B,W}(bounds, window)
end

function AxisIterator(axis, sz; first_pad=nothing, last_pad=nothing, stride=nothing, dilation=nothing)
    return _iter(
        eachindex(axis),
        _k2i(axis, sz),
        _k2i(axis, first_pad),
        _k2i(axis, last_pad),
        _k2i(axis, stride),
        _k2i(axis, dilation)
    )
end

Base.length(w::AxisIterator) = length(getfield(w, :bounds))

@inline function Base.iterate(w::AxisIterator)
    itr = iterate(getfield(w, :bounds))
    if itr === nothing
        return nothing
    else
        return _iterate(first(itr), getfield(w, :window)), last(itr)
    end
end

@inline function Base.iterate(w::AxisIterator, state)
    itr = iterate(getfield(w, :bounds), state)
    if itr === nothing
        return nothing
    else
        return _iterate(first(itr), getfield(w, :window)), last(itr)
    end
end

_iterate(itr, w::AbstractUnitRange) = (first(w) + itr):(last(w) + itr)
_iterate(itr, w::OrdinalRange) = (first(w) + itr):step(w):(last(w) + itr)

@inline function Base.first(itr::AxisIterator)
    return _iterate(first(getfield(itr, :bounds)), getfield(itr, :window))
end

@inline function Base.last(itr::AxisIterator)
    return _iterate(last(getfield(itr, :bounds)), getfield(itr, :window))
end

Base.show(io::IO, ::MIME"text/plain", itr::AxisIterator) = print_axis_iterator(io, itr)

function print_axis_iterator(io, itr::AxisIterator)
    print(io, "AxisIterator(")
    print(io, "($(first(itr))):" * "$(step(getfield(itr, :bounds)))" *":($(last(itr)))")
    print(io, ")")
end

"""
    AxesIterator

N-dimensional iterator of `AxisIterator`s.
"""
struct AxesIterator{I<:Tuple}
    iterators::I
end

function AxesIterator(
    axs::NTuple{N,Any},
    window_size::NTuple{N,Any};
    first_pad=ntuple(_ -> nothing, Val(N)),
    last_pad=ntuple(_ -> nothing, Val(N)),
    stride=ntuple(_ -> nothing, Val(N)),
    dilation=ntuple(_ -> nothing, Val(N))
) where {N}

    itrs = map(_iter, axs, map(_to_size, axs, sz), first_pad, last_pad, stride, dilation)
    return AxesIterator{typeof(itrs)}(itrs)
end

function AxesIterator(A::AbstractArray, sz; kwargs...)
    return iterate_indices(A; window_size=sz, kwargs...)
end

Base.length(itr::AxesIterator) = prod(map(length, getfield(itr, :iterators)))

inc(::Tuple{}, ::Tuple{}, ::Tuple{}) = ()
@inline function inc(itr::Tuple{Any}, olditr::Tuple{Any}, state::Tuple{Any})
    newitr = iterate(first(itr), first(state))
    if newitr === nothing
        return nothing
    else
        return (first(newitr),), (last(newitr),)
    end
end

@inline function inc(itr::Tuple{Any,Vararg{Any}}, olditr::Tuple{Any,Vararg{Any}}, state::Tuple{Any,Vararg{Any}})
    subitr = first(itr)
    substate = first(state)
    nextsubitr = iterate(subitr, substate)
    if nextsubitr === nothing
        subitr, substate = iterate(subitr)
        nextitr = inc(tail(itr), tail(olditr), tail(state))
        if nextitr === nothing
            return nothing
        else
            return (subitr, first(nextitr)...), (substate, last(nextitr)...)
        end
    else
        return (first(nextsubitr), tail(olditr)...), (last(nextsubitr), tail(state)...)
    end
end

@inline firstinc(itr::Tuple{Any}) = iterate(first(itr))

@inline function firstinc(itr::Tuple{Any,Any})
    itr_i = iterate(first(itr))
    if itr_i === nothing
        return nothing
    else
        nextitr = firstinc(tail(itr))
        if nextitr === nothing
            return nothing
        else
            return (first(itr_i), first(nextitr)), (last(itr_i), last(nextitr))
        end
    end
end

@inline function firstinc(itr::Tuple{Any,Vararg{Any}})
    itr_i = iterate(first(itr))
    if itr_i === nothing
        return nothing
    else
        nextitr = firstinc(tail(itr))
        if nextitr === nothing
            return nothing
        else
            return (first(itr_i), first(nextitr)...), (last(itr_i), last(nextitr)...)
        end
    end
end

function Base.iterate(itr::AxesIterator)
    newitr = firstinc(itr.iterators)
    if newitr === nothing
        return nothing
    else
        return Iterators.ProductIterator(first(newitr)), newitr
    end
end

function Base.iterate(itr::AxesIterator, state)
    if state === nothing
        return nothing
    else
        newitrs = inc(itr.iterators, first(state), last(state))
        if newitrs === nothing
            return nothing
        else
            return (Iterators.ProductIterator(first(newitrs)), newitrs)
        end
    end
end

_first(itr) = map(first, getfield(itr, :iterators))
_last(itr) = map(last, getfield(itr, :iterators))

Base.first(itr::AxesIterator) = Iterators.ProductIterator(_first(itr))
Base.last(itr::AxesIterator) = Iterators.ProductIterator(_last(itr))

function Base.show(io::IO, ::MIME"text/plain", itr::AxesIterator)
    print(io, "AxesIterator:\n")
    itrs = getfield(itr, :iterators)
    N = length(itrs)
    for i in 1:N
        print(io, " • ")
        print_axis_iterator(io, itrs[i])
        if i != N
            print(io, "\n")
        end
    end
end

###
###
###

### ensure that x is a proper integer type instead of a key
_k2i(axis, ::Nothing) = nothing
function _k2i(axis, x)
    if is_key(x)
        if keys_type(axis) <: AbstractUnitRange{<:Integer} 
            return Integer(x)
        else
            return Integer(div(x, step(keys(axis))))
        end
    else
        return Integer(x)
    end
end

@inline function _iter(
    inds::AbstractUnitRange{<:Integer},
    ws::Union{Nothing,Integer},
    fp::Union{Nothing,Integer},
    lp::Union{Nothing,Integer},
    s::Union{Nothing,Integer},
    d::Union{Nothing,Integer}
)
    if ws === nothing
        if s !== nothing
            @warn "Ignoring `stride` argument because cannot have a strides" *
                   " between windows without windows argument."
        end
        stop = lp === nothing ? last(inds) : last(inds) - lp
        if d === nothing
            if fp === nothing
                if known_first(inds) === oneunit(eltype(inds))
                    return Base.OneTo(stop)
                else
                    return first(inds):stop
                end
            else
                return (first(inds) + fp):stop
            end
        else
            return (fp === nothing ? first(inds) : first(inds) + fp):d:stop
        end
    else
        if d === nothing
            window = Base.OneTo(ws)
        else
            fi = first(inds)
            window = fi:d:(fi + ws - 1)
        end

        start = fp === nothing ? (first(inds) - 1) : (first(inds) - 1 + fp)
        stop = lp === nothing ? (last(inds) - ws) : (last(inds) - ws - lp)
        #=
        if stride === nothing
            bounds = range(start, step=sz, stop=stop)
        else
            bounds = range(start, step=_to_size(axis, stride) + sz, stop = stop)
        end
        =#

        bounds = range(start, step = (s === nothing ? ws : (ws + s)), stop=stop)

        return AxisIterator{typeof(bounds),typeof(window)}(bounds, window)
    end
end

function _iter(inds::Tuple, ws::Tuple, fp::Tuple, lp::Tuple, s::Tuple, d::Tuple)
    itrs = map(
        _iter,
        inds,
        map(_k2i, inds, ws),
        map(_k2i, inds, fp),
        map(_k2i, inds, lp),
        map(_k2i, inds, s),
        map(_k2i, inds, d)
    )

    return AxesIterator{typeof(itrs)}(itrs)
end

# TODO document this
"""
    iterate_indices(axis[, window_size, first_pad=nothing, last_pad=nothing,
                    stride=nothing, dilation=nothing])

Produces an iterator of windows that can be used to iterate along axes.

* `window_size`: the size of each window. Defaults to the entirety length of the axis.
* `first_pad`: number of elements to skip along the beginning of an axis. defaults to `0`
* `last_pad`: number of elements to skip at the end of an axis. defaults to `0`
* `stride`: strides between windows. defaults to `0`
* `dilation`: steps between each element within a given window. defaults to `1`

## Examples

```jldoctest axis_iterator_examples
julia> using AxisIndices

julia> axis = Axis(range(2.0, step=3.0, length=20))
Axis(2.0:3.0:59.0 => Base.OneTo(20))

julia> iterate_indices(axis, window_size=3)
AxisIterator((1:3):3:(16:18))
```
The final print out indicates that the first window is `1:3` and all subsequent
iterations move by `3` until reaching `16:18`.


The size of the window may be determined by providing an explicit size or the size
in terms of the keys of an axis.
```jldoctest axis_iterator_examples
julia> collect(AxisIterator(axis, 3))
6-element Array{Any,1}:
 1:3
 4:6
 7:9
 10:12
 13:15
 16:18

julia> collect(iterate_indices(axis, window_size=9.0))
6-element Array{Any,1}:
 1:3
 4:6
 7:9
 10:12
 13:15
 16:18

```

The iterator may start with padding from the beginning..
```jldoctest axis_iterator_examples
julia> collect(iterate_indices(axis, window_size=3, first_pad=1))
6-element Array{Any,1}:
 2:4
 5:7
 8:10
 11:13
 14:16
 17:19

julia> collect(iterate_indices(axis, window_size=9.0, first_pad=3.0))
6-element Array{Any,1}:
 2:4
 5:7
 8:10
 11:13
 14:16
 17:19

```
...and the end.
```jldoctest axis_iterator_examples
julia> collect(iterate_indices(axis, window_size=3, first_pad=1, last_pad=2))
5-element Array{Any,1}:
 2:4
 5:7
 8:10
 11:13
 14:16

julia> collect(iterate_indices(axis, window_size=9.0, first_pad=3.0, last_pad=6.0))
5-element Array{Any,1}:
 2:4
 5:7
 8:10
 11:13
 14:16

```

The window can be dilated so that a regular but non-continuous range of elements
are indexed.
```jldoctest axis_iterator_examples
julia> collect(iterate_indices(axis, window_size=3, first_pad=1, last_pad=2, dilation=2))
5-element Array{Any,1}:
 2:2:4
 5:2:7
 8:2:10
 11:2:13
 14:2:16

julia> collect(iterate_indices(axis, window_si=9.0, first_pad=3.0, last_pad=6.0, dilation=6.0))
5-element Array{Any,1}:
 2:2:4
 5:2:7
 8:2:10
 11:2:13
 14:2:16

```

Regular strides can be placed between each iteration.
```jldoctest axis_iterator_examples
julia> collect(iterate_indices(axis, window_size=3, first_pad=1, last_pad=2, stride=2))
3-element Array{Any,1}:
 2:4
 7:9
 12:14

julia> collect(iterate_indices(axis, window_size=9.0, first_pad=3.0, last_pad=6.0, stride=6.0))
3-element Array{Any,1}:
 2:4
 7:9
 12:14

```
"""
function iterate_indices(
    axis::AbstractUnitRange{<:Integer};
    window_size=nothing,
    first_pad=nothing,
    last_pad=nothing,
    stride=nothing,
    dilation=nothing
)

    return _iter(
        eachindex(axis),
        _k2i(axis, window_size),
        _k2i(axis, first_pad),
        _k2i(axis, last_pad),
        _k2i(axis, stride),
        _k2i(axis, dilation)
    )
end


"""

    iterate_indices(A[, window_size, first_pad=nothing, last_pad=nothing,
                    stride=nothing, dilation=nothing])

A

## Examples
```jldoctest
julia> using AxisIndices

julia> AxesIterator(CartesianAxes((20, 20, 20)), (3,3,3))
AxesIterator:
 • AxisIterator((1:3):3:(16:18))
 • AxisIterator((1:3):3:(16:18))
 • AxisIterator((1:3):3:(16:18))

```
"""
function iterate_indices(
    A;
    window_size=nothing,
    first_pad=nothing,
    last_pad=nothing,
    stride=nothing,
    dilation=nothing
)

    return _iter(
        indices(A),
        _to_tuple(window_size, Val(ndims(A))),
        _to_tuple(first_pad, Val(ndims(A))),
        _to_tuple(last_pad, Val(ndims(A))),
        _to_tuple(stride, Val(ndims(A))),
        _to_tuple(dilation, Val(ndims(A)))
    )
end

_to_tuple(x::NTuple{N,<:Any}, ::Val{N}) where {N} = x
_to_tuple(x, nd::Val{N}) where {N} = ntuple(_ -> x, nd)
function _to_tuple(x::NTuple{N,<:Any}, ::Val{M}) where {N,M}
    @assert N < M "received arguments for $N dimensions but array only has $M."
    return (x..., ntuple(_ -> x, Val(M - N))...)
end


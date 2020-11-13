
# TODO it would probably be better to return a tuple of indices instead of a product
# iterator for NDSlidingWindow

### ensure that x is a proper integer type instead of a key
_k2i(axis, ::Nothing, d) = d
_k2i(axis, x::Integer, d) = Int(x)
_k2i(axis, x::StaticInt, d) = x
@inline _k2i(axis, x, d) = Int(div(x, step(keys(axis))))

_to_tuple(x::NTuple{N,<:Any}, ::Val{N}) where {N} = x
_to_tuple(x, nd::Val{N}) where {N} = ntuple(_ -> x, nd)
function _to_tuple(x::NTuple{N,<:Any}, ::Val{M}) where {N,M}
    @assert N < M "received arguments for $N dimensions but array only has $M."
    return (x..., ntuple(_ -> x, Val(M - N))...)
end

"""
    SlidingWindow{B,W}

An iterators for a sliding window along one dimension of an axis.
"""
struct SlidingWindow{B<:AbstractRange,W<:AbstractRange}
    bounds::B
    window::W
end

"""
    NDSlidingWindow{I}

N-dimensional iterator of `SlidingWindow`s.
"""
struct NDSlidingWindow{I<:Tuple}
    sliding_windows::I
end

Base.ndims(x::NDSlidingWindow) = ndims(typeof(x))
Base.ndims(::Type{T}) where {N,I<:Tuple{Vararg{Any,N}},T<:NDSlidingWindow{I}} = N
Base.size(x::NDSlidingWindow) = map(length, x.sliding_windows)
Base.IteratorSize(::Type{T}) where {T<:NDSlidingWindow} = Base.HasShape{ndims(T)}()
@generated function Base.eltype(::Type{T}) where {I,T<:NDSlidingWindow{I}}
    return Iterators.ProductIterator{Tuple{map(eltype, I.parameters)...}}
end

# used for `collapse`
@inline function _all_firsts(itr::SlidingWindow)
    wf = static_first(itr.window)
    return (static_first(itr.bounds) + wf):static_step(itr.bounds):(static_last(itr.bounds) + wf)
end
_all_firsts(itr::AbstractRange) = itr

@inline function sliding_window(
    x;
    window_size=nothing,
    first_pad=nothing,
    last_pad=nothing,
    stride=nothing,
    dilation=nothing
)

    if known_step(x) === 1
        return _sliding_window(x, window_size, first_pad, last_pad, stride, dilation)
    else
        return NDSlidingWindow(_sliding_window(axes(x), window_size, first_pad, last_pad, stride, dilation))
    end
end

@inline function _sliding_window(x::AbstractRange, ws, fp, lp, s, d)
    return _sliding_window(
        eachindex(x),
        _k2i(x, ws, static_length(x)),
        _k2i(x, fp, Zero()),
        _k2i(x, lp, Zero()),
        _k2i(x, s, Zero()),
        _k2i(x, d, One())
    )
end
@inline function _sliding_window(inds::AbstractRange, ws::Integer, fp::Integer, lp::Integer, s::Integer, d::Integer)
    fi = static_first(inds)
    w = fi:d:(fi + ws - oneunit(ws))
    b = (static_first(inds) - oneunit(fp) + fp):(ws + s):(static_last(inds) - ws - lp)
    return SlidingWindow{typeof(b),typeof(w)}(b, w)
end
@inline function _sliding_window(axs::Tuple{Vararg{<:Any,N}}, ws, fp, lp, s, d) where {N}
    D = Val(N)
    return _sliding_window(
        axs,
        _to_tuple(ws, D),
        _to_tuple(fp, D),
        _to_tuple(lp, D),
        _to_tuple(s, D),
        _to_tuple(d, D)
    )
end
@inline function _sliding_window(axs::Tuple, ws::Tuple, fp::Tuple, lp::Tuple, s::Tuple, d::Tuple)
    return map(_sliding_window, axs, ws, fp, lp, s, d)
end

@inline function Base.iterate(w::SlidingWindow)
    itr = _iterate_bounds(getfield(w, :bounds))
    if itr === nothing
        return nothing
    else
        return _iterate(first(itr), getfield(w, :window)), last(itr)
    end
end
@inline function Base.iterate(w::SlidingWindow, state)
    itr = _iterate_bounds(getfield(w, :bounds), state)
    if itr === nothing
        return nothing
    else
        return _iterate(first(itr), getfield(w, :window)), last(itr)
    end
end

@inline function _iterate_bounds(b::AbstractRange)
    if isempty(b)
        return nothing
    elseif known_step(b) === nothing
        # we make this `Int` b/c subsequent iterations won't have a static first value
        # so this makes the type consistent.
        next = Int(first(b))
        return (next, next)
    else
        return (first(b), first(b))
    end
end
@inline function _iterate_bounds(b::AbstractRange, i)
    if i == last(b)
        return nothing
    else
        next = i + static_step(b)
        return (next, next)
    end
end

@inline function _iterate(itr, w::AbstractRange)
    if known_step(w) === 1
        return OptionallyStaticUnitRange(static_first(w) + itr, static_last(w) + itr)
    else
        return OptionallyStaticStepRange(static_first(w) + itr, static_step(w), static_last(w) + itr)
    end
end

@inline function Base.eltype(::Type{T}) where {B,W,T<:SlidingWindow{B,W}}
    if known_step(B) === nothing
        F = Int
        L = Int
    else
        F = known_first(W) === nothing ? Int : StaticInt
        L = known_last(W) === nothing ? Int : StaticInt
    end
    if known_step(W) === 1
        return ArrayInterface.OptionallyStaticUnitRange{F,L}
    else
        if known_step(W) === nothing
            return ArrayInterface.OptionallyStaticStepRange{F,Int,L}
        else
            return ArrayInterface.OptionallyStaticStepRange{F,StaticInt{known_step(W)},L}
        end
    end
end

@inline function Base.first(itr::SlidingWindow)
    return _iterate(static_first(getfield(itr, :bounds)), getfield(itr, :window))
end

@inline function Base.last(itr::SlidingWindow)
    return _iterate(static_last(getfield(itr, :bounds)), getfield(itr, :window))
end

Base.length(itr::SlidingWindow) = static_length(getfield(itr, :bounds))
Base.length(itr::NDSlidingWindow) = prod(map(length, getfield(itr, :sliding_windows)))

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

function Base.iterate(itr::NDSlidingWindow)
    newitr = firstinc(itr.sliding_windows)
    if newitr === nothing
        return nothing
    else
        return Iterators.ProductIterator(first(newitr)), newitr
    end
end

function Base.iterate(itr::NDSlidingWindow, state)
    if state === nothing
        return nothing
    else
        newitrs = inc(itr.sliding_windows, first(state), last(state))
        if newitrs === nothing
            return nothing
        else
            return (Iterators.ProductIterator(first(newitrs)), newitrs)
        end
    end
end

_first(itr) = map(first, getfield(itr, :sliding_windows))
_last(itr) = map(last, getfield(itr, :sliding_windows))

Base.first(itr::NDSlidingWindow) = Iterators.ProductIterator(_first(itr))
Base.last(itr::NDSlidingWindow) = Iterators.ProductIterator(_last(itr))

function Base.show(io::IO, ::MIME"text/plain", itr::NDSlidingWindow)
    print(io, "NDSlidingWindow:\n")
    itrs = getfield(itr, :sliding_windows)
    N = length(itrs)
    for i in 1:N
        print(io, " • ")
        print(io, itrs[i])
        if i != N
            print(io, "\n")
        end
    end
end

"""
    lag(A::AbstractArray, n::Integer)

Shift the elements of `A` along the the axis of the dimension `dim` by `nshift`
elements later. If `dim` is not specified then the dimension returned by
`timedim` is used. If `A` does not have a time dimension then the last dimension
is assumed to be the time dimension.

## Examples
```jldoctest
julia> using TimeAxes

julia> using Unitful: s

julia> A = NamedAxisArray{(:time,)}(collect(1:5), (1:5)s)
5-element NamedAxisArray{Int64,1}
 • time - 1 s:1 s:5 s

  1 s   1
  2 s   2
  3 s   3
  4 s   4
  5 s   5

julia> lag(A, 1)
4-element NamedAxisArray{Int64,1}
 • time - 2 s:1 s:5 s

  2 s   1
  3 s   2
  4 s   3
  5 s   4

julia> lag([1 2 3; 4 5 6; 7 8 9], 1, 1)
2×3 Array{Int64,2}:
 1  2  3
 4  5  6

julia> lag([1 2 3; 4 5 6; 7 8 9], 1, 2)
3×2 Array{Int64,2}:
 1  2
 4  5
 7  8

```
"""
@inline lag(x, nshift; dims) = _lag(x, nshift, to_dims(x, dims))
@inline function _lag(x, nshift::Integer, dim)
    inds = _lag_indices(Val(N), axes(x), nshift, dim)
    return unsafe_reconstruct(
        x,
        @inbounds(x[inds...]);  # TODO ensure that we are really producing inbounds indices
        axes=to_axes(x, inds)
    )
end
@inline function lag(x, nshift::Int)
    if has_timedim(A)
        return _lag(x, nshift, timedim(x))
    else
        return _lag(x, nshift, N)
    end
end

@inline function _lag_indices(::Val{N}, inds::Tuple, nshift::Integer, dim::Integer) where {N}
    ntuple(Val(N)) do i
        if i === dim
            index = getfield(inds, i)
            firstindex(index):(lastindex(index) - nshift)
        else
            Colon()
        end
    end
end

"""
    lead(A::AbstractArray, nshift::Integer[, dim::Integer])

Shift the elements of `A` along the the axis of the dimension `dim` by `nshift`
elements earlier. If `dim` is not specified then the dimension returned by
`timedim` is used. If `A` does not have a time dimension then the last dimension
is assumed to be the time dimension.

## Examples
```jldoctest
julia> using TimeAxes

julia> using Unitful: s

julia> A = NamedAxisArray{(:time,)}(collect(1:5), (1:5)s)
5-element NamedAxisArray{Int64,1}
 • time - 1 s:1 s:5 s

  1 s   1
  2 s   2
  3 s   3
  4 s   4
  5 s   5

julia> lead(A, 1)
4-element NamedAxisArray{Int64,1}
 • time - 1 s:1 s:4 s

  1 s   2
  2 s   3
  3 s   4
  4 s   5

julia> lead([1 2 3; 4 5 6; 7 8 9], 1, 1)
2×3 Array{Int64,2}:
 4  5  6
 7  8  9

julia> lead([1 2 3; 4 5 6; 7 8 9], 1, 2)
3×2 Array{Int64,2}:
 2  3
 5  6
 8  9

```
"""
@inline lead(x, nshift::Integer; dims) = _lead(x, nshift, to_dims(x, dims))
@inline function _lead(x, nshift::Integer, dims)
    inds = _lead_indices(Val(N), axes(x), nshift, dim)
    return unsafe_reconstruct(
        x,
        @inbounds(x[inds...]);  # TODO ensure that we are really producing inbounds indices
        axes=to_axes(x, inds)
    )
end
@inline function lead(x::AbstractArray{T,N}, nshift::Int) where {T,N}
    if has_timedim(x)
        return _lead(x, nshift, timedim(x))
    else
        return _lead(x, nshift, N)
    end
end

@inline function _lead_indices(::Val{N}, inds::Tuple, nshift::Integer, dim::Integer) where {N}
    ntuple(Val(N)) do i
        index = getfield(inds, i)
        if i === dim
            index = getfield(inds, i)
            (firstindex(index) + nshift):lastindex(index)
        else
            Colon()
        end
    end
end

"""
    collapse(f, x; dims, window_size=nothing, first_pad=nothing, last_pad=nothing,
             stride=nothing, dilation=nothing)

collapse and array `x` along a dimension or tuple of dimensions `dims`. The indices along
collapsed dimensions are specified by additional kwargs. If no other kwargs are specified
the entire dimension is collapsed to a single value.

```jldoctest
julia> x = NamedAxisArray(reshape(1:40, 2, 10, 2); x = [:one, :two], y = 1.0:10.0, z = [:three, :four]);

julia> collapse(sum, x; dims=:y, window_size=2)
2×5×2 NamedDimsArray(AxisArray(::Array{Int64,3}
  • axes:
     x = [:one, :two]
     y = 1.0:2.0:9.0
     z = [:three, :four]
))
[:, :, three] =
        1.0  3.0   5.0   7.0   9.0
  :one  4    12    20    28    36
  :two  6    14    22    30    38

[:, :, four] =
        1.0   3.0   5.0   7.0   9.0
  :one  44    52    60    68    76
  :two  46    54    62    70    78


```
"""
@inline function collapse(
    f,
    x;
    dims,
    window_size=nothing,
    first_pad=nothing,
    last_pad=nothing,
    stride=nothing,
    dilation=nothing
)

    return apply_collapse(
        f,
        x,
        NDSlidingWindow(_collapse_iterator(x, dims, window_size, first_pad, last_pad, stride, dilation))
    )
end

# produce the iterator that collapses the array
@inline function _collapse_iterator(x, dims, ws, fp, lp, s, d)
    return _collapse_iterator(x, (dims,), ws, fp, lp, s, d)
end
@inline function _collapse_iterator(x, dims::Tuple{Vararg{<:Any,N}}, ws, fp, lp, s, d) where {N}
    D = Val(N)
    return __collapse_iterator(
        axes(x),
        NamedDims.dim(dimnames(x), dims),  # just in case we have named dimensions
        _to_tuple(ws, D),
        _to_tuple(fp, D),
        _to_tuple(lp, D),
        _to_tuple(s, D),
        _to_tuple(d, D)
    )

end
@inline function __collapse_iterator(axs::Tuple{Vararg{<:Any,D}}, dims::Tuple{Vararg{Int,N}}, ws, fp, lp, s, d) where {D,N}
    ntuple(Val(D)) do i
        dims_index = 0
        @inbounds for dim_i in OneTo(N)
            if dims[dim_i] === i
                dims_index = dim_i
                break
            end
        end
        if dims_index === 0
            return eachindex(@inbounds(axs[i]))
        else
            return _sliding_window(
                @inbounds(axs[i]),
                @inbounds(ws[dims_index]),
                @inbounds(fp[dims_index]),
                @inbounds(lp[dims_index]),
                @inbounds(s[dims_index]),
                @inbounds(d[dims_index])
            )
        end
    end
end

# collapse_axes
collapse_axis(axis, itr) = ArrayInterface.to_axis(axis, _all_firsts(itr))

@inline function apply_collapse(f, x, itr)
    if has_namedims(x)
        return unsafe_reconstruct(
            x,
            [f(@inbounds(view(x, i.iterators...))) for i in itr];
            axis = map(collapse_axis, axes(x), itr.sliding_windows),
            dimnames = Val(dimnames(x))
        )
    else
        return unsafe_reconstruct(
            x,
            [f(@inbounds(view(x, i.iterators...))) for i in itr];
            axis = map(collapse_axis, axes(x), itr.sliding_windows)
        )
    end
end

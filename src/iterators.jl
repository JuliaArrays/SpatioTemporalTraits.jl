# TODO it would probably be better to return a tuple of indices instead of a product
# iterator for NDAxisIterator

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
    AxisIterator{B,W}
"""
struct AxisIterator{B<:AbstractRange,W<:AbstractRange}
    bounds::B
    window::W
end

"""
    NDAxisIterator{I}

N-dimensional iterator of `AxisIterator`s.
"""
struct NDAxisIterator{I<:Tuple}
    axis_iterators::I
end

Base.ndims(x::NDAxisIterator) = ndims(typeof(x))
Base.ndims(::Type{T}) where {N,I<:Tuple{Vararg{Any,N}},T<:NDAxisIterator{I}} = N
Base.size(x::NDAxisIterator) = map(length, x.axis_iterators)
Base.IteratorSize(::Type{T}) where {T<:NDAxisIterator} = Base.HasShape{ndims(T)}()
@generated function Base.eltype(::Type{T}) where {I,T<:NDAxisIterator{I}}
    return Iterators.ProductIterator{Tuple{map(eltype, I.parameters)...}}
end

# used for `collapse`
@inline function _all_firsts(itr::AxisIterator)
    wf = static_first(itr.window)
    return (static_first(itr.bounds) + wf):static_step(itr.bounds):(static_last(itr.bounds) + wf)
end
_all_firsts(itr::AbstractRange) = itr

@inline function axis_iterator(
    x;
    window_size=nothing,
    first_pad=nothing,
    last_pad=nothing,
    stride=nothing,
    dilation=nothing
)

    if known_step(x) === 1
        return _axis_iterator(x, window_size, first_pad, last_pad, stride, dilation)
    else
        return NDAxisIterator(_axis_iterator(axes(x), window_size, first_pad, last_pad, stride, dilation))
    end
end

@inline function _axis_iterator(x::AbstractRange, ws, fp, lp, s, d)
    return _axis_iterator(
        eachindex(x),
        _k2i(x, ws, static_length(x)),
        _k2i(x, fp, Zero()),
        _k2i(x, lp, Zero()),
        _k2i(x, s, Zero()),
        _k2i(x, d, One())
    )
end
@inline function _axis_iterator(inds::AbstractRange, ws::Integer, fp::Integer, lp::Integer, s::Integer, d::Integer)
    fi = static_first(inds)
    w = fi:d:(fi + ws - oneunit(ws))
    b = (static_first(inds) - oneunit(fp) + fp):(ws + s):(static_last(inds) - ws - lp)
    return AxisIterator{typeof(b),typeof(w)}(b, w)
end
@inline function _axis_iterator(axs::Tuple{Vararg{<:Any,N}}, ws, fp, lp, s, d) where {N}
    D = Val(N)
    return _axis_iterator(
        axs,
        _to_tuple(ws, D),
        _to_tuple(fp, D),
        _to_tuple(lp, D),
        _to_tuple(s, D),
        _to_tuple(d, D)
    )
end
@inline function _axis_iterator(axs::Tuple, ws::Tuple, fp::Tuple, lp::Tuple, s::Tuple, d::Tuple)
    return map(_axis_iterator, axs, ws, fp, lp, s, d)
end

@inline function Base.iterate(w::AxisIterator)
    itr = _iterate_bounds(getfield(w, :bounds))
    if itr === nothing
        return nothing
    else
        return _iterate(first(itr), getfield(w, :window)), last(itr)
    end
end
@inline function Base.iterate(w::AxisIterator, state)
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

@inline function Base.eltype(::Type{T}) where {B,W,T<:AxisIterator{B,W}}
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

@inline function Base.first(itr::AxisIterator)
    return _iterate(static_first(getfield(itr, :bounds)), getfield(itr, :window))
end

@inline function Base.last(itr::AxisIterator)
    return _iterate(static_last(getfield(itr, :bounds)), getfield(itr, :window))
end

Base.length(itr::AxisIterator) = static_length(getfield(itr, :bounds))
Base.length(itr::NDAxisIterator) = prod(map(length, getfield(itr, :axis_iterators)))

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

function Base.iterate(itr::NDAxisIterator)
    newitr = firstinc(itr.axis_iterators)
    if newitr === nothing
        return nothing
    else
        return Iterators.ProductIterator(first(newitr)), newitr
    end
end

function Base.iterate(itr::NDAxisIterator, state)
    if state === nothing
        return nothing
    else
        newitrs = inc(itr.axis_iterators, first(state), last(state))
        if newitrs === nothing
            return nothing
        else
            return (Iterators.ProductIterator(first(newitrs)), newitrs)
        end
    end
end

_first(itr) = map(first, getfield(itr, :axis_iterators))
_last(itr) = map(last, getfield(itr, :axis_iterators))

Base.first(itr::NDAxisIterator) = Iterators.ProductIterator(_first(itr))
Base.last(itr::NDAxisIterator) = Iterators.ProductIterator(_last(itr))

function Base.show(io::IO, ::MIME"text/plain", itr::NDAxisIterator)
    print(io, "NDAxisIterator:\n")
    itrs = getfield(itr, :axis_iterators)
    N = length(itrs)
    for i in 1:N
        print(io, " • ")
        print(io, itrs[i])
        if i != N
            print(io, "\n")
        end
    end
end


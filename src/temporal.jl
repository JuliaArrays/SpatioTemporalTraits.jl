
"""
    is_temporal(x) -> StaticBool

Returns `True` if `x` refers to time.
"""
is_temporal(x::Symbol) = is_temporal(static(x))
is_temporal(x) = is_temporal(typeof(x))
is_temporal(::Type{T}) where {T} = static(false)
is_temporal(::Type{<:Union{StaticSymbol{:time},StaticSymbol{:Time}}}) = static(true)
is_temporal(::Type{<:Dates.AbstractTime}) = static(true)

function ArrayInterface.to_dims(::Type{T}, ::typeof(is_temporal)) where {T}
    d = _find_timedim(static(ndims(T)), dimnames(T))
    if d === nothing
        throw(ArgumentError("$T does not have a temporal dimension."))
    else
        return d
    end
end

_find_timedim(d::StaticInt{0}, n::Tuple) where {D} = nothing
@inline function _find_timedim(d::StaticInt{D}, n::Tuple) where {D}
    _find_timedim(is_temporal(@inbounds(getfield(n, D))), d, n)
end
_find_timedim(::True, d::StaticInt{D}, n::Tuple) where {D} = d
_find_timedim(::False, d::StaticInt{D}, n::Tuple) where {D} = _find_timedim(d - static(1), n)

"""
    has_timedim(x) -> Bool

Returns 'true' if `x` has a time dimension.
"""
has_timedim(x) = _find_timedim(static(ndims(x)), dimnames(x)) !== nothing

"""
    timedim(::Type{T}) -> Integer

Returns the dimension associated with time.
"""
timedim(x) = timedim(typeof(x))
timedim(::Type{T}) where {T} = ArrayInterface.to_dims(T, is_temporal)

"""
    ntimes(x) -> Integer

Returns the number of elements along the time dimension.
"""
ntimes(x) = size(x, timedim(x))

"""
    times(x)

Returns a collection of times associated with `x`. If `x` has a specific dimension
associated with time, this returns the keys along that dimension.
"""
@inline times(x) = keys(ArrayInterface.axes(x, is_temporal))

"""
    time_last(x)

Last time point along the time axis.
"""
time_last(x) = last(times(x))

"""
    time_first(x)

First time point along the time axis.
"""
time_first(x) = first(times(x))

"""
    time_step(x)

The time step/interval between each element.
"""
time_step(x) = step(times(x))

"""
    duration(x)

Duration of the event along the time axis.
"""
duration(x) = time_last(x) - time_first(x) + time_step(x)

"""
    sampling_rate(x)

Number of samples per second.
"""
sampling_rate(x) = inv(time_step(x))

"""
    assert_timedim_last(x)

Throw an error if the `x` has a time dimension that is not the last dimension.
"""
@inline function assert_timedim_last(x)
    if _find_timedim(static(ndims(x)), dimnames(x)) != ndims(x)
        throw(ArgumentError("time dimension is not last"))
    else
        return nothing
    end
end



"""
    is_time(x) -> StaticBool

Returns `True` if `x` refers to time.
"""
is_time(x::Symbol) = is_time(static(x))
is_time(x) = is_time(typeof(x))
is_time(::Type{T}) where {T} = False()
is_time(::Type{StaticSymbol{:time}}) = True()
is_time(::Type{StaticSymbol{:Time}}) = True()

@inline function ArrayInterface.to_dims(::Type{T}, ::typeof(is_time)) where {T}
    out = _to_timedim(dimnames(T), nstatic(Val(ndims(T))))
    if out === nothing
        ArrayInterface.no_dimname_error(T, :time)
    else
        return out
    end
end
_to_timedim(::Tuple{}, ::Tuple{}) = nothing
function _to_timedim(dn::Tuple{Vararg{Any}}, dims::Tuple{Vararg{Any}})
    if is_spatial(first(dn))
        return first(dims)
    else
        return _to_spatialdims(tail(dn), tail(dims))
    end
end

"""
    timedim(::Type{T}) -> Integer

Returns the dimension associated with time.
"""
timedim(x) = timedim(typeof(x))
timedim(::Type{T}) where {T} = ArrayInterface.to_dims(T, is_time)

"""
    ntimes(x) -> Integer

Returns the number of elements along the time dimension.
"""
ntimes(x) = size(x, timedim(x))

"""
    time_axis(x)

Returns the axis associated with time.
"""
time_axis(x) = axes(x, timedim(x))

"""
    times(x)

Returns a collection of times associated with `x`. If `x` has a specific dimension
associated with time, this returns the keys along that dimension.
"""
times(x) = keys(time_axis(x))

"""
    time_end(x)

Last time point along the time axis.
"""
time_end(x) = last(times(x))

"""
    onset(x)

First time point along the time axis.
"""
onset(x) = first(times(x))

"""
    time_step(x)

The time step/interval between each element.
"""
time_step(x) = step(times(x))

"""
    duration(x)

Duration of the event along the time axis.
"""
duration(x) = time_end(x) - onset(x) + time_step(x)

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
    if has_timedim(x)
        if timedim(x) === ndims(x)
            return nothing
        else
            error("time dimension is not last")
        end
    else
        return nothing
    end
end


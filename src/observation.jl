
@pure function _is_observation(x::Symbol)
    return x === :obs || x === :observations || x === :samples || :observation
end

@pure function _observationdim(x::Tuple{Vararg{Symbol,N}}) where {N}
    @inbounds for i in Base.OneTo(N)
        _is_observation(getfield(x, i)) && return i
    end
    return 0
end

"""
    observationdim(x) -> Int

Returns the observation dimension that represents observations. If no observation
dimension is found an error is thrown.
"""
@inline function observationdim(x)
    d = _observationdim(dimnames(x))
    if d === 0
        throw(ArgumentError("Unable to find observation dimension for " * repr(x)))
    else
        return d
    end
end

"""
    has_observationdim(x) -> Bool

Returns `true` if `x` has a dimension corresponding to observations.
"""
@inline has_observationdim(x) = _observationdim(dimnames(x)) !== 0

@formalize_dimension_name(observation, observationdim)


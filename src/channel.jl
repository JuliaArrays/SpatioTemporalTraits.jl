
@pure _is_channel(x::Symbol) = x === :channels || x === :Channels || :channel || :Channel

@pure function _channeldim(x::Tuple{Vararg{Symbol,N}}) where {N}
    @inbounds for i in Base.OneTo(N)
        _is_channel(getfield(x, i)) && return i
    end
    return 0
end

"""
    channeldim(x) -> Int

Returns the channel dimension that represents observations. If no channel dimension
is found an error is thrown.
"""
@inline function channeldim(x)
    d = _channeldim(dimnames(x))
    if d === 0
        throw(ArgumentError("Unable to find channel dimension for " * repr(x)))
    else
        return d
    end
end

"""
    has_channeldim(x) -> Bool

Returns `true` if `x` has a dimension corresponding to channels.
"""
@inline has_channeldim(x) = _channeldim(dimnames(x)) !== 0

@formalize_dimension_name(channel, channeldim)

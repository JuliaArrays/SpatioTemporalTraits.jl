
#= TODO we shouldn't need these
@generated function _spatialdims(::Val{L}) where {L}
    out = Expr(:tuple)
    for i in Base.OneTo(length(L))
        if _is_spatial(L[i])
            push!(out.args, i)
        end
    end
    quote
        return $out
    end
end

spatialdims(img::AbstractMappedArray) = spatialdims(parent(img))
function spatialdims(img::AbstractMultiMappedArray)
    ps = traititer(spatialdims, parent(img)...)
    checksame(ps)
end
@inline function spatialdims(img::SubArray)
    return _subarray_offset(0, spatialdims(parent(img)), img.indices...)
end
@inline _subarray_offset(off, x, i::Real, inds...) =
    _subarray_offset(off-1, tail(x), inds...)
@inline _subarray_offset(off, x, i, inds...) =
    (x[1]+off, _subarray_offset(off, tail(x), inds...)...)
_subarray_offset(off, x::Tuple{}) = ()
@inline function spatialdims(img::Base.PermutedDimsArray{T,N,perm,iperm}) where {T,N,perm,iperm}
    return _getindex_tuple(spatialdims(parent(img)), iperm)
end
=#

#=
spatial_size(img::AbstractMappedArray) = size_spatial(parent(img))
function spatial_size(img::AbstractMultiMappedArray)
    ps = traititer(spatial_size, parent(img)...)
    checksame(ps)
end

@inline function spatial_size(img::SubArray)
    return _subarray_filter(spatial_size(parent(img)), img.indices...)
end
@inline function spatial_size(img::Base.PermutedDimsArrays.PermutedDimsArray{T,N,perm,iperm}) where {T,N,perm,iperm}
    return _getindex_tuple(spatial_size(parent(img)), iperm)
end
=#

"""
    affine_map(x) -> AffineMap

Returns and affine map. By default using `spatial_directions` and `pixel_spacing`
are used to constuct the mapping.
"""
function affine_map(x)
    return AffineMap(_spatial_directions_to_rotation(RotMatrix, spatial_directions(x)),
                     _pixelspacing_to_linearmap(pixel_spacing(x)))
end

function _pixelspacing_to_linearmap(ps::NTuple{2,T}) where {T}
    return @inbounds LinearMap(SVector(Float64(ps[1]), Float64(ps[2]), 0.0))
end

function _pixelspacing_to_linearmap(ps::NTuple{3,T}) where {T}
    return @inbounds LinearMap(SVector(Float64(ps[1]), Float64(ps[2]), Float64(ps[3])))
end

function _spatial_directions_to_rotation(::Type{R}, sd::NTuple{2,NTuple{2,T}}) where {R,T}
    return @inbounds R(SMatrix{3,3,Float64,9}(
        sd[1][1], sd[2][1], 0,
        sd[1][2], sd[2][2], 0,
               0,        0, 1)
    )
end

function _spatial_directions_to_rotation(::Type{R}, sd::NTuple{3,NTuple{3,T}}) where {R,T}
    return @inbounds R(SMatrix{3,3,Float64,9}(
        sd[1][1], sd[2][1], sd[3][1],
        sd[1][2], sd[2][2], sd[3][2],
        sd[1][3], sd[2][3], sd[3][3])
    )
end

# FIXME
spatial_directions(img::AbstractMappedArray) = spatial_directions(parent(img))
function spatial_directions(img::AbstractMultiMappedArray)
    ps = traititer(spatial_directions, parent(img)...)
    checksame(ps)
end
@inline function spatial_directions(img::SubArray)
    return _subarray_filter(spatial_directions(parent(img)), getfield(img, :indices)...)
end

@inline function spatial_directions(img::Base.PermutedDimsArray{T,N,perm}) where {T,N,perm}
    return _getindex_tuple(spatial_directions(parent(img)), perm)
end

function _spatial_directions(ps::NTuple{N,Any}) where N
    return ntuple(i->ntuple(d->d==i ? ps[d] : zero(ps[d]), Val(N)), Val(N))
end

@inline _subarray_filter(x, i::Real, inds...) = _subarray_filter(tail(x), inds...)
@inline _subarray_filter(x, i, inds...) = (x[1], _subarray_filter(tail(x), inds...)...)
_subarray_filter(x::Tuple{}) = ()

@inline _getindex_tuple(t::Tuple, inds::Tuple) =
    (t[first(inds)], _getindex_tuple(t, tail(inds))...)
_getindex_tuple(t::Tuple, ::Tuple{}) = ()

@inline traititer(f, A, rest...) = (f(A), traititer(f, rest...)...)
traititer(f) = ()

function checksame(t::Tuple)
    val1 = t[1]
    @assert all(p -> p == val1, t)
    return val1
end

# TODO decide what to do about MappedArrays

using MappedArrays
using MappedArrays: AbstractMultiMappedArray, MultiMappedArray, ReadonlyMultiMappedArray

function has_dimnames(::Type{ReadonlyMultiMappedArray{T,N,AAs,F}}) where {T,N,AAs,F}
    return _multi_array_has_dimnames(AAs)
end
function has_dimnames(::Type{MultiMappedArray{T,N,AAs,F,Finv}}) where {T,N,AAs,F,Finv}
    return _multi_array_has_dimnames(AAs)
end

# FIXME this doesn't account for when there are incompatable names from multiple arrays
@inline function _multi_array_has_dimnames(::Type{T}) where {T}
    for T_i in T.parameters
        has_dimnames(T_i) && return true
    end
    return false
end

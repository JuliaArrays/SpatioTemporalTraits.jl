
using Test

using SpatioTemporalTraits
using ArrayInterface
using ArrayInterface: dimnames, to_dims
using Metadata
using Static
using Unitful
using Unitful: m, mm, ft, s

struct KeyedAxis{K,P<:AbstractUnitRange{Int}} <: AbstractUnitRange{Int}
    keys::K
    parent::P
end
KeyedAxis(k) = KeyedAxis(k, eachindex(k))

Base.keys(x::KeyedAxis) = getfield(x, :keys)

Base.parent(x::KeyedAxis) = getfield(x, :parent)

Base.first(x::KeyedAxis) = first(parent(x))
Base.last(x::KeyedAxis) = last(parent(x))

ArrayInterface.parent_type(::Type{<:KeyedAxis{K,P}}) where {K,P} = P

Base.@propagate_inbounds Base.getindex(x::KeyedAxis, i::Integer) = parent(x)[i]
Base.@propagate_inbounds Base.getindex(x::KeyedAxis, i::AbstractRange) = parent(x)[i]

struct NamedIndices{dnames,T,N,P<:Union{CartesianIndices{N},LinearIndices{N}}} <: AbstractArray{T,N}
    parent::P

    function NamedIndices{D}(x::Union{LinearIndices,CartesianIndices}) where {D}
        N = ndims(x)
        return new{D::NTuple{N,Symbol},eltype(x),N,typeof(x)}(x)
    end
end
NamedIndices(x::Union{LinearIndices,CartesianIndices}) = NamedIndices{ntuple(_->:_, Val(ndims(x)))}(x)

ArrayInterface.dimnames(::Type{<:NamedIndices{dnames}}) where {dnames} = static(dnames)
function ArrayInterface.dimnames(::Type{<:NamedIndices{dnames}}, dim) where {dnames}
    if dim > length(dnames)
        return static(:_)
    else
        return static(getfield(dnames, Int(dim)))
    end
end

ArrayInterface.parent_type(::Type{<:NamedIndices{<:Any,<:Any,<:Any,P}}) where {P} = P

Base.parent(x::NamedIndices) = getfield(x, :parent)

Base.getindex(x::NamedIndices, i::Vararg) = parent(x)[i...]

Base.size(x::NamedIndices) = size(parent(x))
#Base.size(x::NamedIndices, dim) = Int(ArrayInterface.size(x, to_dims(x, dim)))

Base.axes(x::NamedIndices) = ArrayInterface.axes(parent(x))
Base.axes(x::NamedIndices, dim::Integer) = _axes(axes(x), dim)
function _axes(axs::Tuple{Vararg{Any,N}}, dim) where {N}
    if dim > N
        return static(1):static(1)
    else
        return @inbounds(axs[dim])
    end
end


x = NamedIndices{(:x, :y, :z, :time)}(
    CartesianIndices((
        KeyedAxis((1:2:9)m),
        KeyedAxis((1.0:2:9)ft),
        KeyedAxis((2:2:10)mm),
        KeyedAxis((1:4)s)
    )));

no_time = view(x, :, :, 1:2, 1);

## No time dimension
@test_throws ArgumentError timedim(no_time)
@test !has_timedim(no_time)
#@test ntimes(no_time) == 1
@test @inferred(pixel_spacing(no_time)) === (2m, 2.0ft, 2mm)
@test @inferred(spatial_directions(no_time)) === ((2m, 0.0ft, 0mm), (0m, 2.0ft, 0mm), (0m, 0.0ft, 2mm))
@test @inferred(spatialdims(no_time)) === (static(1), static(2), static(3))
@test @inferred(spatial_order(no_time)) === (static(:x), static(:y), static(:z))
@test @inferred(spatial_indices(no_time)) == ((1:2:9)m, (1.0:2.0:9.0)ft, (2:2:4)mm)
@test @inferred(origin(no_time)) == (1m, 1.0ft, 2mm)
@test @inferred(spatial_last(no_time)) == (9m, 9.0ft, 4mm)

mx = attach_metadata(no_time, Dict{Symbol,Any}());
@test origin(mx) == (1m, 1.0ft, 2mm)
metadata!(mx, :origin, (0,0,0))
@test origin(mx) == (0, 0, 0)
@test spatial_directions(mx) === ((2m, 0.0ft, 0mm), (0m, 2.0ft, 0mm), (0m, 0.0ft, 2mm))
metadata!(mx, :spatial_directions, ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
@test spatial_directions(mx) === ((1, 0, 0), (0, 1, 0), (0, 0, 1))

## time dimension
assert_timedim_last(x)
@test @inferred(times(x)) == (1:4)s
@test @inferred(has_timedim(x))
@test @inferred(timedim(x)) === static(4)
@test @inferred(ntimes(x)) == 4



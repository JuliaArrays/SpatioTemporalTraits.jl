
"""
    is_spatial(x::StaticSymbol) -> StaticBool

Returns `static(true)` if `x` refers to a spatial dimension. By default, any dimension that
is not one of the following is considered spatial:
* time, Time
* channels, Channels
* observations, Observations, obs
"""
is_spatial(x::Symbol) = is_spatial(static(x))
is_spatial(x) = is_spatial(typeof(x))
is_spatial(::Type{T}) where {T} = static(false)
is_spatial(::Type{StaticSymbol{sym}}) where {sym} = static(true)
const NotSpatial = Union{StaticSymbol{:time},StaticSymbol{:Time}, StaticSymbol{:Channels},
    StaticSymbol{:channels},StaticSymbol{:color}, StaticSymbol{:Color},
    StaticSymbol{:observations},StaticSymbol{:Observations},StaticSymbol{:obs}}
is_spatial(::Type{NotSpatial}) = static(false)


@inline function ArrayInterface.to_dims(::Type{T}, ::typeof(is_spatial)) where {T}
    _to_sdims(static(0), dimnames(T), Static.nstatic(Val(ndims(T))))
end

_to_sdims(::StaticInt, ::Tuple{}, ::Tuple{}) = ()
@inline function _to_sdims(
    ::StaticInt{N},
    dn::Tuple{StaticSymbol{name},Vararg{Any}},
    dims::Tuple{StaticInt{D},Vararg{Any}}
) where {N,name,D}

    if known(is_spatial(static(name)))
        if N === 2
            return (static(D),)
        else
            return (static(D), _to_sdims(static(N + 1), tail(dn), tail(dims))...)
        end
    else
        return _to_sdims(static(N + 1), tail(dn), tail(dims))
    end
end

"""
    spatialdims(x) -> Tuple{Vararg{StaticInt}}

Return a tuple of the spatial dimensions of `x`. 
"""
spatialdims(@nospecialize(x)) = spatialdims(typeof(x))
@inline spatialdims(::Type{T}) where {T} = _spatialdims(has_dimnames(T), T)
@inline _spatialdims(::True, ::Type{T}) where {T} = to_dims(T, is_spatial)
@inline _spatialdims(::False, ::Type{T}) where {T} = ntuple(identity, Val(min(ndims(x), 3)))

@inline spatial_axes(x) = _spatial_axes(x, spatialdims(x))
_spatial_axes(x, ::Tuple{}) = ()
@inline function _spatial_axes(x, d::Tuple{StaticInt{N},Vararg{Any}}) where {N}
    (ArrayInterface.axes(x, static(N)), _spatial_axes(x, tail(d))...)
end

"""
    spatial_order(x) -> Tuple

Returns the `dimnames` of `x` that correspond to spatial dimensions, as determined by
[`is_spatial`](@ref). At most three dimensions will be considered spatial.
"""
spatial_order(@nospecialize(x)) = spatial_order(typeof(x))
@inline spatial_order(::Type{T}) where {T} = Static.permute(dimnames(T), spatialdims(T))

"""
    sdims(x) -> Int

Return the number of spatial dimensions that `x` has.

See also: [`spatialdims`](@ref), [`is_spatial`](@ref)
"""
@inline sdims(x) = length(spatialdims(x))

"""
    spatial_indices(x)

Returns the keys along each spatial dimension.
"""
@inline spatial_indices(x) = map(keys, spatial_axes(x))
spatial_indices(x::MetaArray) = spatial_indices(parent(x))
# first get the spatial_indices of the parent array and then index into the ones that aren't
# dropped
function spatial_indices(x::SubArray)
    map(getindex,
        Static.permute(spatial_indices(parent(x)), to_spatialdims(typeof(x))),
        map(subdim -> getfield(x, :indices)[subdim], map(d -> to_parent_dims(x, d), spatialdims(x)))
    )
end
@inline function spatial_indices(x::Union{PermutedDimsArray,Adjoint,Transpose})
    Static.permute(spatial_indices(parent(x)), to_spatialdims(typeof(x)))
end

"""
    pixel_spacing(x)

Return a tuple representing the separation between adjacent pixels along each axis
of the image. Derived from the step size of each element of `spatial_indices`.
"""
@inline pixel_spacing(x) = map(_axis_pixel_spacing, spatial_indices(x))
@inline function pixel_spacing(x::SubArray)
    Static.permute(pixel_spacing(parent(x)), to_spatialdims(typeof(x)))
end
_axis_pixel_spacing(x::AbstractRange) = step(x)
_axis_pixel_spacing(x) = 1

"""
    origin(x) -> Tuple

Returns the spatial origin of `x`.
"""
origin(x) = spatial_first(x)
@inline origin(x::MetaArray) = getmeta(spatial_first, x, :origin)
@inline function origin(x::Union{PermutedDimsArray,SubArray,Adjoint,Transpose})
    Static.permute(origin(parent(x)), to_spatialdims(typeof(x)))
end

"""
    spatial_directions(img) -> (axis1, axis2, ...)

Return a tuple-of-tuples, each `axis[i]` representing the displacement
vector between adjacent pixels along spatial axis `i` of the image
array, relative to some external coordinate system ("physical
coordinates").

By default this is computed from `pixel_spacing`, but you can set this
manually using ImageMeta.
"""
@inline spatial_directions(x::X) where {X} = _spatial_directions(pixel_spacing(x))
function _spatial_directions(ps::NTuple{N,Any}) where {N}
    return ntuple(Val(N)) do i
        ntuple(Val(N)) do d
            if d === i
                getfield(ps, d)
            else
                zero(getfield(ps, d))
            end
        end
    end
end
@inline function spatial_directions(x::MetaArray)
    if has_metadata(x, :spatial_directions)
        return metadata(x, :spatial_directions)
    else
        return spatial_directions(parent(x))
    end
end
@inline function spatial_direction(x::Union{PermutedDimsArray,SubArray,Adjoint,Transpose})
    Static.permute(spatial_directions(parent(x)), to_spatialdims(typeof(x)))
end

"""
    spatial_first(x)

The first position along each spatial indices (i.e., where each spatial axis starts).
"""
spatial_first(x) = map(first, spatial_indices(x))

"""
    spatial_last(x)

The last position along each spatial indices (i.e., where each spatial axis stops).
"""
spatial_last(x) = map(last, spatial_indices(x))

function to_spatialdims(::Type{T}) where {T}
    _to_spatialdims(spatialdims(parent_type(T)), to_parent_dims(T))
end
@generated function _to_spatialdims(::SD, ::TD) where {SD,TD}
    out = Expr(:tuple)
    for p in TD.parameters
        i = findfirst(==(known(p)), known(SD))
        if i !== nothing
            push!(out.args, :(static($i)))
        end
    end
    Expr(:block, Expr(:meta, :inline), out)
end

"""
    width(x)

Returns the size of the dimension corresponding to width.
"""
width(x) = _width(spatialdims(x), x)
_width(::Union{Tuple{},Tuple{Any}}, x) = 1
_width(dims::Union{Tuple{Any,Any},Tuple{Any,Any,Any}}, x) = size(x, getfield(dims, 2))

"""
    height(x)

Returns the size of the dimension corresponding to height.
"""
height(x) = _height(spatialdims(x), x)
_height(::Tuple{}, x) = 1
_height(dims::Tuple, x) = size(x, getfield(dims, 1))

"""
    depth(x)

Returns the size of the dimension corresponding to depth.
"""
depth(x) = _depth(spatialdims(x), x)
_depth(::Union{Tuple{},Tuple{Any},Tuple{Any,Any}}, x) = 1
_depth(dims::Tuple{Any,Any,Any}, x) = size(x, getfield(dims, 3))


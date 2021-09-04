
"""
    is_spatial(x::StaticSymbol) -> StaticBool

Returns `static(true)` if `x` refers to a spatial dimension. By default, any dimension that
isn't `:time` or `:Time` is considered spatial.
"""
is_spatial(x::Symbol) = is_spatial(static(x))
is_spatial(x) = is_spatial(typeof(x))
is_spatial(::Type{T}) where {T} = static(false)
is_spatial(::Type{StaticSymbol{sym}}) where {sym} = static(true)
is_spatial(::Type{StaticSymbol{:time}}) = static(false)
is_spatial(::Type{StaticSymbol{:Time}}) = static(false)

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

#=
function _to_spatialdims(dn::Tuple{Vararg{Any}}, dims::Tuple{Vararg{Any}})
    return _to_spatialdims(is_spatial(first(dn)), dn, dims)
end
function _to_spatialdims(::True, dn::Tuple{Vararg{Any}}, dims::Tuple{Vararg{Any}})
    return (first(dims), _to_spatialdims(tail(dn), tail(dims))...)
end
function _to_spatialdims(::False, dn::Tuple{Vararg{Any}}, dims::Tuple{Vararg{Any}})
    return _to_spatialdims(tail(dn), tail(dims))
end
_to_spatialdims(dn::Tuple{}, dims::Tuple{}) = ()

_max3(::Tuple{}) = ()
_max3(x::Tuple{Any}) = x
_max3(x::Tuple{Any,Any}) = x
_max3(x::Tuple{Any,Any,Any}) = x
function _max3(x::Tuple{Any,Any,Any,Vararg{Any}})
    (@inbounds(getfield(x, 1)), @inbounds(getfield(x, 2)), @inbounds(getfield(x, 3)))
end

_spatialdims()
=#


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
sdims(@nospecialize(x)) = sdims(typeof(x))
sdims(::Type{T}) where {T} = _sdims(has_dimnames(T), T)
_sdims(::True, ::Type{T}) where {T} = length(spatialdims(x))
_sdims(::False, ::Type{T}) where {T} = min(ndims(x), 3)

# FIXME we should use something like `is_regularly_sampled` to ensure this makes sense
"""
    spatial_size(x) -> Tuple{Vararg{Int}}

Return a tuple listing the sizes of the spatial dimensions of the image.
"""
@inline spatial_size(x) =  map(_spatial_length, spatial_indices(x))
_spatial_length(s) = last(s) - first(s) + step(s)

"""
    spatial_indices(x)

Returns the keys along each spatial dimension.
"""
@inline spatial_indices(x) = map(keys, spatial_axes(x))
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
@inline function origin(x::MetaArray)
    if has_metadata(x, :origin)
        return metadata(x, :origin)
    else
        return spatial_first(x)
    end
end
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
@inline function spatial_directions(x::X) where {X}
    if has_metadata(x, :spatial_directions)
        return spatial_directions(parent(x))
    else
        return _spatial_directions(pixel_spacing(x))
    end
end
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


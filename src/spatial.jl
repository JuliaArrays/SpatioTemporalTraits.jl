
"""
    is_spatial(x) -> StaticBool

Returns `True` if `x` refers to a spatial dimension.
"""
is_spatial(x::Symbol) = is_spatial(static(x))
is_spatial(x) = is_spatial(typeof(x))
is_spatial(::Type{T}) where {T} = False()
is_spatial(::Type{StaticSymbol{sym}}) where {sym} = True()
is_spatial(::Type{StaticSymbol{:time}}) = False()
is_spatial(::Type{StaticSymbol{:Time}}) = False()

@inline function ArrayInterface.to_dims(::Type{T}, ::typeof(is_spatial)) where {T}
    return _max3(_to_spatialdims(dimnames(T), nstatic(Val(ndims(T)))))
end

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
    return (@inbounds(x[1]), @inbounds(x[2]), @inbounds(x[3]))
end

"""
    spatialdims(x) -> Tuple{Vararg{Int}}

Return a tuple listing the spatial dimensions of `img`.
Note that a better strategy may be to use ImagesAxes and take slices along the time axis.
"""
spatialdims(x) = spatialdims(typeof(x))
function spatialdims(::Type{T}) where {T}
    if has_dimnames(T)
        return to_dims(T, is_spatial)
    else
        return ntuple(identity, Val(min(ndims(x), 3)))
    end
end

"""
    spatial_order(x) -> Tuple{Vararg{Symbol}}

Returns the `dimnames` of `x` that correspond to spatial dimensions.
"""
spatial_order(x) = spatial_order(typeof(x))
@inline function spatial_order(::Type{T}) where {T}
    if has_dimnames(T)
        return eachop(getindex, dimnames(T), spatialdims(T))
    else
        throw(MethodError(spatial_order, (T,)))
    end
end

"""
    spatial_axes(x) -> Tuple

Returns a tuple of each axis corresponding to a spatial dimensions.
"""
@inline function spatial_axes(x)
    if has_dimnames(x)
        return eachop(axes, x, spatialdims(x))
    else
        return ntuple(i -> axes(x, i), Val(min(ndims(x), 3)))
    end
end

"""
    spatial_size(x) -> Tuple{Vararg{Int}}

Return a tuple listing the sizes of the spatial dimensions of the image.
"""
@inline spatial_size(x) =  eachop(size, x, spatialdims(x))

"""
    spatial_keys(x)

Returns the keys along each spatial dimension.
"""
@inline spatial_keys(x) = map(keys, spatial_axes(x))

"""
    pixel_spacing(x)

Return a tuple representing the separation between adjacent pixels along each axis
of the image. Derived from the step size of each element of `spatial_keys`.
"""
@inline pixel_spacing(x) = map(axis_pixel_spacing, spatial_axes(x))
axis_pixel_spacing(x) = _axis_pixel_spacing(keys(x))
_axis_pixel_spacing(x::AbstractRange) = step(x)
_axis_pixel_spacing(x) = 1

"""
    spatial_directions(img) -> (axis1, axis2, ...)

Return a tuple-of-tuples, each `axis[i]` representing the displacement
vector between adjacent pixels along spatial axis `i` of the image
array, relative to some external coordinate system ("physical
coordinates").

By default this is computed from `pixel_spacing`, but you can set this
manually using ImageMeta.
"""
spatial_directions(img::AbstractArray) = _spatial_directions(pixel_spacing(img))
function _spatial_directions(ps::NTuple{N,Any}) where {N}
    return ntuple(Val(N)) do i
        ntuple(Val(N)) do d
            if d === i
                ps[d]
            else
                zero(ps[d])
            end
        end
    end
end

"""
    spatial_offset(x)

The offset of each dimension (i.e., where each spatial axis starts).
"""
spatial_offset(x) = map(_axis_spatial_offset, spatial_axes(x))
_axis_spatial_offset(x) = first(keys(x))

# FIXME what qualifies as a spatial dimension needs to be documented somewhere
"""
    sdims(x)

Return the number of spatial dimensions in the image. Defaults to the same as
`ndims`, but with `NamedDimsArray` you can specify that some dimensions correspond
to other quantities (e.g., time) and thus not included by `sdims`.
"""
sdims(x) = sdims(typeof(x))
function sdims(::Type{T}) where {T}
    if has_dimnames(T)
        return length(spatialdims(x))
    else
        return min(ndims(x), 3)
    end
end


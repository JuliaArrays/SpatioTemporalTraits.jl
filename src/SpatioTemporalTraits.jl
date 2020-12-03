module SpatioTemporalTraits

using ArrayInterface
using ArrayInterface: dimnames, has_dimnames, known_first, known_last, known_step,
    static_first, static_last, static_length, static_step, StaticInt, to_dims,
    One, OptionallyStaticUnitRange, OptionallyStaticStepRange, Zero
using Base: tail, OneTo
using Base.Cartesian: @nif

export
    assert_timedim_last,
    duration,
    each_time,
    has_timedim,
    ntimes,
    onset,
    pixel_spacing,
    sampling_rate,
    select_time,
    spatialdims,
    spatial_order,
    spatial_axes,
    spatial_size,
    spatial_keys,
    spatial_offset,
    spatial_directions,
    timedim,
    times,
    time_end,
    time_step

struct NotSpatial{DimName} end
_not_spatial(::NotSpatial{DimName}) where {DimName} = false
@inline not_spatial(x::Symbol) = _not_spatial(NotSpatial{x}())

"""
    @not_spatial(x::Symbol)

Declares dimensions with the name `x` to not be spatial dimensions.
"""
macro not_spatial(x)
    esc(quote
        SpatioTemporalTraits._not_spatial(::SpatioTemporalTraits.NotSpatial{$x}) = true
        nothing
    end)
end

"""
    spatialdims(x) -> Tuple{Vararg{Int}}

Return a tuple listing the spatial dimensions of `img`.
Note that a better strategy may be to use ImagesAxes and take slices along the time axis.
"""
@inline function spatialdims(x)
    if has_dimnames(x)
        return _spatialdims(Val(dimnames(x)))
    else
        return ntuple(identity, Val(min(ndims(x), 3)))
    end
end

macro ifelse(condition, trueout, falseout)
    esc(Expr(:if, condition, trueout, falseout))
end

function _generate_spatialdims_body(N::Int)
    quote
        @nif(
            $N,
            d1 -> !not_spatial(getfield(L, d1)),
            d1 -> @nif(
                $N - d1,
                d2 -> !not_spatial(getfield(L, d2 + d1)),
                d2 -> @nif(
                    $N - (d2 + d1),
                    d3 -> !not_spatial(getfield(L, d3 + d2 + d1)),
                    d3 -> (d1, d2 + d1, d3 + d2 + d1),
                    d3 -> @ifelse(
                        !not_spatial(getfield(L, d3 + d2 + d1)),
                        (d1, d2 + d1, d3 + d2 + d1),
                        (d1, d2 + d1)
                    )
                ),
                d2 -> @ifelse(
                    !not_spatial(getfield(L, d2 + d1)),
                    (d1, d2 + d1),
                    (d1,)
                )
            ),
            d1 -> @ifelse(
                !not_spatial(getfield(L, d1)),
                (d1,),
                ()
            )
        )
    end
end

@generated function _spatialdims(::Val{L}) where {L}
    _generate_spatialdims_body(length(L))
end

@generated function _filter_by_spatialdims(::Val{D}, x::Tuple) where {D}
    out = Expr(:tuple)
    for i in D
        push!(out.args, :(getfield(x, $i)))
    end
    return Expr(:block, Expr(:meta, :inline), out)
end

function _generate_namedim_body(N::Int, f::Function)
    quote
        Base.Cartesian.@nif $N d->$(f)(getfield(L, d)) d->d d->0
    end
end

"""
    spatial_order(x) -> Tuple{Vararg{Symbol}}

Returns the `dimnames` of `x` that correspond to spatial dimensions.
"""
@inline function spatial_order(x)
    if has_dimnames(x)
        return _filter_by_spatialdims(spatialdims(x), dimnames(x))
    else
        throw(MethodError(spatial_order, (x,)))
    end
end

"""
    spatial_axes(x) -> Tuple

Returns a tuple of each axis corresponding to a spatial dimensions.
"""
@inline function spatial_axes(x)
    if has_dimnames(x)
        return _filter_by_spatialdims(spatialdims(x), axes(x))
    else
        return ntuple(i -> axes(x, i), Val(min(ndims(x), 3)))
    end
end

"""
    spatial_size(x) -> Tuple{Vararg{Int}}

Return a tuple listing the sizes of the spatial dimensions of the image.
"""
@inline spatial_size(x) =  map(i -> size(x, i), spatialdims(x))

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
    spatial_offset(x)

The offset of each dimension (i.e., where each spatial axis starts).
"""
spatial_offset(x) = map(axis_spatial_offset, spatial_axes(x))
axis_spatial_offset(x) = first(keys(x))

# FIXME what qualifies as a spatial dimension needs to be documented somewhere
"""
    sdims(x)

Return the number of spatial dimensions in the image. Defaults to the same as
`ndims`, but with `NamedDimsArray` you can specify that some dimensions correspond
to other quantities (e.g., time) and thus not included by `sdims`.
"""
@inline function sdims(x)
    if has_dimnames(x)
        return length(spatialdims(x))
    else
        return min(ndims(x), 3)
    end
end

"""
    spatial_directions(img) -> (axis1, axis2, ...)

Return a tuple-of-tuples, each `axis[i]` representing the displacement
vector between adjacent pixels along spatial axis `i` of the image
array, relative to some external coordinate system ("physical
coordinates").

By default this is computed from `pixel_spacing`, but you can set this
manually using ImagesMeta.
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

abstract type AbstractNamedDim{name} <: Function end

ArrayInterface.to_dims(x, d::AbstractNamedDim) = d(x)

Base.show(io::IO, p::AbstractNamedDim) = _show_nameddim(io, p)
Base.show(io::IO, ::MIME"text/plain", p::AbstractNamedDim) = _show_nameddim(io, p)
function _show_nameddim(io, nd::AbstractNamedDim{name}) where {name}
    fxnname = Symbol(name, :dim)
    nmethods = length(methods(nd))
    if nmethods == 1
        print(io, "$fxnname (generic function with 1 method)")
    else
        print(io, "$fxnname (generic function with $nmethods methods)")
    end
end

"""
    @defdim(name, dim_not_found[, is_a_spatialdim::Bool=false])

Defines internal code for giving assigning multiple `Symbol` names to a single meaning
for a dimension. For example, the following code formalizes internal code to identify and
access the time dimension (as specified by the first argument `time`) and the second
argument is an anonymous function that acts on anything that doesn't have a time dimension.
In this case, if a time dimension isn't identified we assume it is 1 more than the last
dimension.

```julia
@defdim time (x -> ndims(x) + 1)

@is_time :time

@is_time :Time
```

The last two arguments use the `@is_time` macro produced by `@defdim` to assert that any
dimensions with named `:time` or `:Time` should be considered time dimensions.

"""
macro defdim(name, dim_not_found, is_a_spatialdim::Bool=false)
    IsName = Symbol(:Is, uppercasefirst(string(name)))

    _is_name = Symbol(:_is_, name)

    is_name = Symbol(:is_, name)

    _namedim = Symbol(:_, name, :dim)

    namedim_no_error = Symbol(name, :dim_no_error)

    _namedim = Symbol(:_, name, :dim)

    NameDim = Symbol(uppercasefirst(string(name)), :Dim)
    namedim = Symbol(name, :dim)

    if is_a_spatialdim
        macro_is_name = Expr(:macro,
            Expr(:call, is_name, :DimName),
            Expr(:block, Expr(:quote, Expr(:block,
                esc(Expr(:function,
                    Expr(:call, _is_name, Expr(:(::), Expr(:curly, IsName, Expr(:$, :DimName)))),
                    true
                )),
                nothing
        ))))
    else
        macro_is_name = Expr(:macro,
            Expr(:call, is_name, :DimName),
            Expr(:block, Expr(:quote, Expr(:block,
                esc(Expr(:function,
                    Expr(:call, _is_name, Expr(:(::), Expr(:curly, IsName, Expr(:$, :DimName)))),
                    true
                )),
                Expr(:macrocall,
                    Expr(:., :SpatioTemporalTraits, QuoteNode(Symbol("@not_spatial"))),
                    LineNumberNode(@__LINE__),
                    Expr(:$, :DimName)
                ),
                nothing
        ))))
    end

    namedim_doc = """
        $(namedim)(x) -> Int

    Returns the dimension that represents $(name).
    """

    macro_doc = """
        @$(is_name)(x::Symbol)

    Declares any dimension with the name `x` to is a $(name) dimension. This also declares
    any dimension with the name `x` is _not_ a spatial dimension.
    """

    has_namedim = Symbol(:has_, name, :dim)
    has_namedim_doc = """
        $(has_namedim)(x) -> Bool

    Returns `true` if `x` has a dimension corresponding to $(name).
    """
 
    nnames = Symbol(:n, name, :s)
    nnames_doc = """
        $nnames(x) -> Int

    Returns the number of unique $name values. If `x` has a dimension corresponding to $name
    then the length of the corresponding axis is returned.
    """

    names = Symbol(name, :s)
    names_doc = """
        $names(x)

    Returns the keys corresponding to the $name axis.
    """

    name_axis = Symbol(name, :_axis)
    name_axis_doc = """
        $name_axis(x)

    Returns the axis corresponding to the $name dimension.
    """

    select_name = Symbol(:select_, name)
    select_name_doc = """
        $select_name(x, i)

    Return a view of all the data of `x` where the index for the $name dimension equals `i`.
    """

    each_name = Symbol(:each_, name)
    each_name_doc = """
        $each_name(x)

    Create a generator that iterates over the $name dimensions `x`, returning views that select
    all the data from the other dimensions in `x`.
    """

    esc(quote
        struct $(NameDim) <: SpatioTemporalTraits.AbstractNamedDim{$(name)} end
        const $namedim = $(NameDim)()

        struct $(IsName){DimName} end

        $(_is_name)(::$(IsName){DimName}) where {DimName} = false

        @inline $is_name(x::Symbol) = $_is_name($IsName{x}())

        @generated function $_namedim(::Val{L}) where {L}
            SpatioTemporalTraits._generate_namedim_body(length(L))
        end

        @doc $namedim_doc
        @inline function $namedim(x)
            if has_dimnames(x)
                d = $_namedim(Val(dimnames(x)))
                if d === 0
                    return $dim_not_found(x)
                else
                    return d
                end
            else
                return $dim_not_found(x)
            end
        end

        @doc $has_namedim_doc
        @inline function $has_namedim(x)
            return has_dimnames(x) && ($_namedim(Val(dimnames(x))) !== 0)
        end

        $macro_is_name

        @doc $macro_doc $(Symbol("@", is_name))

        @doc $nnames_doc
        @inline $nnames(x) = Base.size(x, $namedim(x))

        @doc $names_doc
        @inline $names(x) = keys(axes(x, $namedim(x)))

        @doc $name_axis_doc
        @inline $name_axis(x) = eachindex(axes(x, $namedim(x)))

        @doc $select_name_doc
        @inline $select_name(x, i) = selectdim(x, $namedim(x), i)

        @doc $each_name_doc
        @inline $each_name(x) = eachslice(x, dims=$namedim(x))

         nothing
    end)
end

@defdim time (x -> ndims(x) + 1)

@is_time :time

@is_time :Time

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

end # module


"""
    has_dimnames(x) -> Bool

Returns `true` if `x` has names for each dimension.
"""
has_dimnames(x) = has_dimnames(typeof(x))
function has_dimnames(::Type{T}) where {T}
    if parent_type(T) <: T
        return false
    else
        return has_dimnames(parent_type(T))
    end
end

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


# TODO port docs to SpatioTemporalTraits
"""
    named_axes(A) -> NamedTuple{names}(axes)

Returns a `NamedTuple` where the names are the dimension names and each indice
is the corresponding dimensions's axis. If dimnesion names are not defined for `x`
default names are returned. `x` should have an `axes` method.

```julia
julia> using AxisIndices

julia> A = reshape(1:24, 2,3,4);

julia> named_axes(A)
(dim_1 = Base.OneTo(2), dim_2 = Base.OneTo(3), dim_3 = Base.OneTo(4))

julia> named_axes(NamedAxisArray{(:a, :b, :c)}(A))
(a = SimpleAxis(Base.OneTo(2)), b = SimpleAxis(Base.OneTo(3)), c = SimpleAxis(Base.OneTo(4)))
```
"""
function named_axes(x)
    if has_dimnames(x)
        return NamedTuple{dimnames(x)}(axes(x))
    else
        return NamedTuple{default_names(Val(ndims(x)))}(axes(x))
    end
end

@generated default_names(::Val{N}) where {N} = :($(ntuple(i -> Symbol(:dim_, i), N)))

macro formalize_dimension_name(name, name_dim, dim_noerror_name)
    nname = Symbol(:n, name)
    nname_doc = """
        $nname(x) -> Int

    Returns the size along the dimension corresponding to the $name.
    """

    name_keys = Symbol(name, :_keys)
    name_keys_doc = """
        $name_keys(x)

    Returns the keys corresponding to the $name axis
    """

    name_indices = Symbol(name, :_indices)
    name_indices_doc = """
        $name_indices(x)

    Returns the indices corresponding to the $name axis
    """
 
    has_namedim = Symbol(:has_, name, :dim)
    has_namedim_doc = """
        $has_namedim(x) -> Bool

    Returns `true` if `x` has a dimension corresponding to $name.
    """

    name_axis = Symbol(name, :_axis)
    name_axis_doc = """
        $name_axis(x)

    Returns the axis corresponding to the $name dimension.
    """

    name_axis_itr = """
        $name_axis(x, window_size=nothing[; first_pad=nothing, last_pad=nothing, stride=nothing, dilation=nothing])

    Returns an `AxisIterator` along the $name axis.
    """


    name_select = Symbol(:select_, name)
    name_select_doc = """
        $name_select(x, i)

    Return a view of all the data of `x` where the index for the $name dimension equals `i`.
    """

    each_name = Symbol(:each_, name)
    each_name_doc = """
        $each_name(x)

    Create a generator that iterates over the $name dimensions `A`, returning views that select
    all the data from the other dimensions in `A`.
    """
 
    esc(quote

        @doc $nname_doc
        @inline $nname(x) = Base.size(x, $name_dim(x))

        @doc $name_keys_doc
        @inline $name_keys(x) = keys(axes(x, $name_dim(x)))

        @doc $name_indices_doc
        @inline $name_indices(x) = ArrayInterface.indices(axes(x, $name_dim(x)))

        @doc $has_namedim_doc
        @inline $has_namedim(x) = !($dim_noerror_name(dimnames(x)) === 0)

        @doc $name_axis_doc
        @inline $name_axis(x) = axes(x, $name_dim(x))

        @doc $name_axis_itr
        @inline function $name_axis(x, sz; kwargs...)
            return SpatioTemporalTraits.iterate_indices(axes(x, $name_dim(x)); window_size=sz, kwargs...)
        end

        @doc $name_select_doc
        @inline $name_select(x, i) = selectdim(x, $name_dim(x), i)

        @doc $each_name_doc
        @inline $each_name(x) = eachslice(x, dims=$name_dim(x))
 
       nothing
    end)
end


# we have to abuse `@pure` a bit with names to ensure constant propagation
@pure _is_observation(x::Symbol) = (x === :obs) | (x === :observations) | (x === :samples)
@pure function _is_channel(x::Symbol)
    return x === :channels || x === :Channels || x === :Color || x === :color
end

@pure function _find_observation_dim(x::Tuple{Vararg{Symbol,N}}) where {N}
    for i in Base.OneTo(N)
        _is_observation(getfield(x, i, false)) && return i
    end
    return 0
end

@pure function _find_channel_dim(x::Tuple{Vararg{Symbol,N}}) where {N}
    for i in Base.OneTo(N)
        _is_channel(getfield(x, i, false)) && return i
    end
    return 0
end

"""
    observationdim(x) -> Int

Returns the observation dimension that represents observations. If no observation
dimension is found an error is thrown.
"""
@inline function observationdim(x)
    d = _find_observation_dim(dimnames(x))
    if d === 0
        throw(ArgumentError("Unable to find observation dimension for " * repr(x)))
    else
        return d
    end
end

@formalize_dimension_name(observation, observationdim, _find_observation_dim)

"""
    channeldim(x) -> Int

Returns the channel dimension that represents observations. If no channel dimension
is found an error is thrown.
"""
@inline function channeldim(x)
    d = _find_channel_dim(dimnames(x))
    if d === 0
        throw(ArgumentError("Unable to find channel dimension for " * repr(x)))
    else
        return d
    end
end

@formalize_dimension_name(channel, channeldim, _find_channel_dim)

@pure _is_time(x::Symbol) = x === :time || x === :Time

@pure function _find_time_dim(x::Tuple{Vararg{Symbol,N}}) where {N}
    for i in Base.OneTo(N)
        _is_time(getfield(x, i, false)) && return i
    end
    return 0
end

"""
    timedim(x) -> Int

Returns the dimension that represents time. If no time dimension is found it is
assumed that this is is a single time point and the time dimension is `ndims(x) + 1`.
"""
@inline function timedim(x)
    d = _find_time_dim(dimnames(x))
    if d === 0
        return ndims(x) + 1
    else
        return d
    end
end


@formalize_dimension_name(time, timedim, _find_time_dim)


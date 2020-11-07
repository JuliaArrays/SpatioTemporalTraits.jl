
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

macro formalize_dimension_name(name, name_dim)
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

    name_indices = Symbol(name, :_indices)
    name_indices_doc = """
        $name_indices(x)

    Returns the indices corresponding to the $name axis.
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

    iterate_name = Symbol(:iterate_, name)
    iterate_name_doc = """
        $iterate_name(x; kwargs...)

    """

    esc(quote
        @doc $nnames_doc
        @inline $nnames(x) = Base.size(x, $name_dim(x))

        @doc $names_doc
        @inline $names(x) = keys(axes(x, $name_dim(x)))

        @doc $name_indices_doc
        @inline $name_indices(x) = eachindex(axes(x, $name_dim(x)))

        @doc $select_name_doc
        @inline $select_name(x, i) = selectdim(x, $name_dim(x), i)

        @doc $each_name_doc
        @inline $each_name(x) = eachslice(x, dims=$name_dim(x))

        @doc $iterate_name_doc
        @inline $iterate_name(x; kwargs...) = axis_iterator(axes(x, $name_dim(x)); kwargs...)

        nothing
    end)
end

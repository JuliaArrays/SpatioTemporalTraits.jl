
@pure _is_time(x::Symbol) = x === :time || x === :Time
@pure function _time_dim(x::Tuple{Vararg{Symbol,N}}) where {N}
    @inbounds for i in Base.OneTo(N)
        _is_time(getfield(x, i)) && return i
    end
    return 0
end

"""
    timedim(x) -> Int

Returns the dimension that represents time. If no time dimension is found it is
assumed that this is is a single time point and the time dimension is `ndims(x) + 1`.
"""
@inline function timedim(x)
    d = _timedim(dimnames(x))
    if d === 0
        return ndims(x) + 1
    else
        return d
    end
end

"""
    has_timedim(x) -> Bool

Returns `true` if `x` has a dimension corresponding to time.
"""
@inline has_timedim(x) = _time_dim(dimnames(x)) !== 0

@formalize_dimension_name(time, timedim)

"""
    time_end(x)

Last time point along the time axis.
"""
time_end(x) = last(time_keys(x))

"""
    onset(x)

First time point along the time axis.
"""
onset(x) = first(time_keys(x))

"""
    time_step(x)

The time step/interval between each element.
"""
time_step(x) = step(time_keys(x))

"""
    duration(x)

Duration of the event along the time axis.
"""
duration(x) = time_end(x) - onset(x) + time_step(x)

"""
    sampling_rate(x)

Number of samples per second.
"""
sampling_rate(x) = 1 / time_step(x)

"""
    assert_timedim_last(x)

Throw an error if the `x` has a time dimension that is not the last dimension.
"""
@inline function assert_timedim_last(x::T) where {T}
    if has_timedim(x)
        if timedim(x) === ndims(T)
            return nothing
        else
            error("time dimension is not last")
        end

    else
        return nothing
    end
end


"""
    lag(A::AbstractArray, n::Integer)

Shift the elements of `A` along the the axis of the dimension `dim` by `nshift`
elements later. If `dim` is not specified then the dimension returned by
`timedim` is used. If `A` does not have a time dimension then the last dimension
is assumed to be the time dimension.

## Examples
```jldoctest
julia> using TimeAxes

julia> using Unitful: s

julia> A = NamedAxisArray{(:time,)}(collect(1:5), (1:5)s)
5-element NamedAxisArray{Int64,1}
 • time - 1 s:1 s:5 s

  1 s   1
  2 s   2
  3 s   3
  4 s   4
  5 s   5

julia> lag(A, 1)
4-element NamedAxisArray{Int64,1}
 • time - 2 s:1 s:5 s

  2 s   1
  3 s   2
  4 s   3
  5 s   4

julia> lag([1 2 3; 4 5 6; 7 8 9], 1, 1)
2×3 Array{Int64,2}:
 1  2  3
 4  5  6

julia> lag([1 2 3; 4 5 6; 7 8 9], 1, 2)
3×2 Array{Int64,2}:
 1  2
 4  5
 7  8

```
"""

@inline function lag(A::AbstractArray{T,N}, nshift::Integer, dim::Integer) where {T,N}
    return A[_lag_indices(Val(N), axes(A), nshift, dim)...]
end

@inline function lag(A::AxisArray{T,N}, nshift::Integer, dim::Int) where {T,N}
    indexing_indices = _lag_indices(Val(N), axes(A), nshift, dim)
    p = parent(A)[indexing_indices...]
    axs = _shift_axes(axes(A), indexing_indices, axes(p), dim)
    return unsafe_reconstruct(A, p, axs)
end

function lag(A::NamedAxisArray, dim::Int, n::Int)
    return NamedDimsArray{dimnames(A)}(lag(parent(A), dim, n))
end

@inline function lag(A::AbstractArray{T,N}, nshift::Int) where {T,N}
    if has_timedim(A)
        return lag(A, nshift, timedim(A))
    else
        return lag(A, nshift, N)
    end
end

@inline function _lag_indices(::Val{N}, inds::Tuple, nshift::Integer, dim::Integer) where {N}
    ntuple(Val(N)) do i
        if i === dim
            index = getfield(inds, i)
            firstindex(index):(lastindex(index) - nshift)
        else
            Colon()
        end
    end
end

"""
    lead(A::AbstractArray, nshift::Integer[, dim::Integer])

Shift the elements of `A` along the the axis of the dimension `dim` by `nshift`
elements earlier. If `dim` is not specified then the dimension returned by
`timedim` is used. If `A` does not have a time dimension then the last dimension
is assumed to be the time dimension.

## Examples
```jldoctest
julia> using TimeAxes

julia> using Unitful: s

julia> A = NamedAxisArray{(:time,)}(collect(1:5), (1:5)s)
5-element NamedAxisArray{Int64,1}
 • time - 1 s:1 s:5 s

  1 s   1
  2 s   2
  3 s   3
  4 s   4
  5 s   5

julia> lead(A, 1)
4-element NamedAxisArray{Int64,1}
 • time - 1 s:1 s:4 s

  1 s   2
  2 s   3
  3 s   4
  4 s   5

julia> lead([1 2 3; 4 5 6; 7 8 9], 1, 1)
2×3 Array{Int64,2}:
 4  5  6
 7  8  9

julia> lead([1 2 3; 4 5 6; 7 8 9], 1, 2)
3×2 Array{Int64,2}:
 2  3
 5  6
 8  9

```
"""
@inline function lead(A::AbstractArray{T,N}, nshift::Integer, dim::Integer) where {T,N}
    return A[_lead_indices(Val(N), axes(A), nshift, dim)...]
end

@inline function lead(A::AxisArray{T,N}, nshift::Integer, dim::Int) where {T,N}
    indexing_indices = _lead_indices(Val(N), axes(A), nshift, dim)
    p = parent(A)[indexing_indices...]
    axs = _shift_axes(axes(A), indexing_indices, axes(p), dim)
    return unsafe_reconstruct(A, p, axs)
end

function lead(A::NamedAxisArray, dim::Int, n::Int)
    return NamedDimsArray{dimnames(A)}(lead(parent(A), dim, n))
end

@inline function lead(A::AbstractArray{T,N}, nshift::Int) where {T,N}
    if has_timedim(A)
        return lead(A, nshift, timedim(A))
    else
        return lead(A, nshift, N)
    end
end

@inline function _lead_indices(::Val{N}, inds::Tuple, nshift::Integer, dim::Integer) where {N}
    ntuple(Val(N)) do i
        index = getfield(inds, i)
        if i === dim
            index = getfield(inds, i)
            (firstindex(index) + nshift):lastindex(index)
        else
            Colon()
        end
    end
end

@inline function _noshift_axis(axis::A, inds::I) where {A,I}
    if is_indices_axis(axis)
        return to_axis(axis, nothing, inds, false)
    else
        return to_axis(axis, keys(axis), inds, false)
    end
end

@inline function _shift_axis(
    axis::AbstractAxis,
    indexing_index::AbstractUnitRange,
    parent_index::AbstractUnitRange
)

    if is_indices_axis(axis)
        return to_axis(
            axis,
            nothing,
            parent_index,
            false,
        )
    else
        return to_axis(
            axis,
            @inbounds(keys(axis)[indexing_index]),
            parent_index,
            false,
        )
    end
end

_shift_axes(::Tuple{}, ::Tuple{}, ::Tuple{}, dim::Int) = ()
@inline function _shift_axes(
    old_axes::Tuple,
    indexing_indices::Tuple,
    parent_indices::Tuple,
    dim::Int
)

    if dim === 1
        return (
            _shift_axis(
                first(old_axes),
                first(indexing_indices),
                first(parent_indices),
            ),
            map(_noshift_axis, tail(old_axes), tail(parent_indices))...
        )
    else
        return (
            _noshift_axis(
                first(old_axes),
                first(parent_indices)
            ),
            _shift_axes(
                tail(old_axes),
                tail(indexing_indices),
                tail(parent_indices),
                dim - 1
            )...
        )
    end
end

"""
    collapse(f, x; dims, window_size=nothing, first_pad=nothing, last_pad=nothing,
             stride=nothing, dilation=nothing)

collapse and array `x` along a dimension or tuple of dimensions `dims`. The indices along
collapsed dimensions are specified by additional kwargs. If no other kwargs are specified
the entire dimension is collapsed to a single value.

```jldoctest
julia> x = NamedAxisArray(reshape(1:40, 2, 10, 2); x = [:one, :two], y = 1.0:10.0, z = [:three, :four]);

julia> collapse(sum, x; dims=:y, window_size=2)
julia> collapse(sum, x; dims=:y, window_size=2)
2×5×2 NamedDimsArray(AxisArray(::Array{Int64,3}
  • axes:
     x = [:one, :two]
     y = 1.0:2.0:9.0
     z = [:three, :four]
))
[:, :, three] =
        1.0  3.0   5.0   7.0   9.0
  :one  4    12    20    28    36
  :two  6    14    22    30    38

[:, :, four] =
        1.0   3.0   5.0   7.0   9.0
  :one  44    52    60    68    76
  :two  46    54    62    70    78


```
"""
@inline function collapse(
    f,
    x;
    dims,
    window_size=nothing,
    first_pad=nothing,
    last_pad=nothing,
    stride=nothing,
    dilation=nothing
)

    return apply_collapse(
        f,
        x,
        NDAxisIterator(_collapse_iterator(x, dims, window_size, first_pad, last_pad, stride, dilation))
    )
end

# produce the iterator that collapses the array
@inline function _collapse_iterator(x, dims, ws, fp, lp, s, d)
    return _collapse_iterator(x, (dims,), ws, fp, lp, s, d)
end
@inline function _collapse_iterator(x, dims::Tuple{Vararg{<:Any,N}}, ws, fp, lp, s, d) where {N}
    D = Val(N)
    return __collapse_iterator(
        axes(x),
        NamedDims.dim(dimnames(x), dims),  # just in case we have named dimensions
        _to_tuple(ws, D),
        _to_tuple(fp, D),
        _to_tuple(lp, D),
        _to_tuple(s, D),
        _to_tuple(d, D)
    )

end
@inline function __collapse_iterator(axs::Tuple{Vararg{<:Any,D}}, dims::Tuple{Vararg{Int,N}}, ws, fp, lp, s, d) where {D,N}
    ntuple(Val(D)) do i
        dims_index = 0
        @inbounds for dim_i in OneTo(N)
            if dims[dim_i] === i
                dims_index = dim_i
                break
            end
        end
        if dims_index === 0
            return eachindex(@inbounds(axs[i]))
        else
            return _axis_iterator(
                @inbounds(axs[i]),
                @inbounds(ws[dims_index]),
                @inbounds(fp[dims_index]),
                @inbounds(lp[dims_index]),
                @inbounds(s[dims_index]),
                @inbounds(d[dims_index])
            )
        end
    end
end

# collapse_axes
collapse_axis(axis, itr) = ArrayInterface.to_axis(axis, _all_firsts(itr))

# all arguments should be ready at this point for applying the collapse
@inline function apply_collapse(f, x::NamedDimsArray{L}, itr) where {L}
    return NamedDimsArray{L}(apply_collapse(f, parent(x), itr))
end
@inline function apply_collapse(f, x::AxisArray, itr)
    ArrayInterface.unsafe_reconstruct(
        x,
        apply_collapse(f, parent(x), itr);
        axes=map(collapse_axis, axes(x), itr.axis_iterators)
    )
end
@inline apply_collapse(f, x, itr) = [f(@inbounds(x[i.iterators...])) for i in itr]


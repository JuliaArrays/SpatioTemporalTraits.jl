# SpatioTemporalTraits

# Warning

At the time of writing this, `SpatioTemporalTraits` is not tested on any specific array implementation.
I fully intend to implement formal testing for this in the AxisIndices.jl package and make myself available for questions/feedback to interested parties maintaining other array types that wish to incorporate these traits.
This does not mean the examples or code shouldn't work.
Rather, it just means the examples and code aren't necessarily useful to end users in the interim.
They may still be informative for those developing code that wants to generically provide support for spatiotemporal data.

# Introduction

`SpatioTemporalTraits` provides spatial and temporal traits on top of [ArrayInterface.jl](https://github.com/SciML/ArrayInterface.jl).
Types that use ArrayInterface.jl have inherent compatibility with methods that utilize `SpatioTemporalTraits`.
For example, the following method returns `x` if it's last dimension is a time dimension and otherwise permutes the last dimensions to be the time dimension.

```julia
using SpatioTemporalTraits

function ensure_timedim_last(x)
    d = timedim(x)  # returns the integer that corresponds to the time dimension
    N = ndims(x)
    if N >= d
        return x
    else
        perms = ntuple(i -> ifelse(N >= i, i + 1, i), Val(N))
        return permutedims(x, perms)
    end
end
```

By default the only dimensions that are considered "time dimensions" are those that return `:time` or `:Time` when `ArrayInterface.dimnames` is called (e.g., `ArrayInterface.dimnames(x::TimeMatrix) = (:observations, :time)` would return `2` when `timedim(x::TimeMatrix)` is called).
Conversely, any all dimensions are considered spatial dimensions except for those with the names `:time` or `:Time`.
This can be changed by adding methods to `SpatioTemporalTraits.is_time` and `SpatioTemporalTraits.is_spatial`
For example, the "frequency" can be made a formal time dimension by doing the following:

```julia
using SpatioTemporalTraits
using Static

SpatioTemporalTraits.is_time(::Type{typeof(static(:frequency))}) = static(true)
SpatioTemporalTraits.is_spatial(::Type{typeof(static(:frequency))}) = static(false)
```

## Non-dependency compatibility with SpatioTemporalTraits

The following method may be called as `ensure_dim_last(x, timedim)` to produce the same result as the previous example, but it doesn't depend directly on SpatioTemporalTraits.
```julia
using ArrayInterface

function ensure_dim_last(x, dim)
    d = ArrayInterface.to_dims(x, dim)  # returns the integer that corresponds to the `dim`
    N = ndims(x)
    if N >= d
        return x
    else
        perms = ntuple(i -> ifelse(N >= i, i + 1, i), Val(N))
        return permutedims(x, perms)
    end
end

```


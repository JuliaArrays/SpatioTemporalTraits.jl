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
There may be instances you want to specify a different name refers to time, which is what `SpatioTemporalTraits.@is_time` is for.
For example, `SpatioTemporalTraits.@is_time :hammer_time` could be used for specifying a dimension of times that you "can't touch this".

## Defining new dimensions

We can create new meaningful dimensions using `@defdim`.
The first argument passed to `@defdim` is the suffix used for methods created (see the docstring for `@defdim` to get a full list of methods created).
The second argument is a method called when no dimension is found corresponding to the newly defined dimension.
In this instance an error is thrown when `x` doesn't have a dimension named `:observations` or `:obs`.
```julia
@defdim(
    observations,
    (x -> throw(ArgumentError("$x does not have a dimension corresponding to 'observations'")))
)

@is_observations :observations

@is_observations :obs
```
`@defdim` creates the `@is_observation` macro for formalizing which `Symbol` names should be considered observation dimensions.
This will internally assert that anything that refers to an observation is not spatial.
If we want to define a dimension that does refer to a spatial dimension we can do the following:
```julia
@defdim(
    horizontal,
    (x -> throw(ArgumentError("$x does not have a dimension corresponding to 'horizontal'"))),
    true
)

@is_horizontal :horizontal
```


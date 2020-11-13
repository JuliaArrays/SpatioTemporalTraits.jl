module SpatioTemporalTraits

using ArrayInterface
using ArrayInterface: known_first, known_last, known_step, static_first, static_last,
    static_length, static_step, StaticInt, One, OptionallyStaticUnitRange,
    OptionallyStaticStepRange, Zero
using Base: tail, OneTo
using Base.Cartesian: @nif
using Metadata
using NamedDims

export
    assert_timedim_last,
    collapse,
    duration,
    each_time,
    has_timedim,
    lag,
    lead,
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

# TODO remove these
has_dimnames(x) = has_dimnames(typeof(x))
has_dimnames(::Type{T}) where {T} = false
has_dimnames(::Type{T}) where {T<:NamedDimsArray} = true

#include("sliding_window.jl")
#include("spatial.jl")
#include("names.jl")
#include("time.jl")
include("names.jl")
#include("methods.jl")

end # module

module SpatioTemporalTraits

using ArrayInterface
using ArrayInterface: False, True, StaticInt, StaticSymbol
using ArrayInterface: dimnames, has_dimnames, eachop, to_dims, static, nstatic
using Base: tail

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
    time_axis,
    time_end,
    time_step

include("spatial.jl")
include("temporal.jl")

end # module

module SpatioTemporalTraits

using ArrayInterface
using ArrayInterface: dimnames, has_dimnames, to_dims, parent_type, to_parent_dims
using Dates
using Metadata
using Metadata: MetaArray, NoMetadata
using Static
using Base: tail

using LinearAlgebra

export
    assert_timedim_last,
    duration,
    has_timedim,
    ntimes,
    origin,
    pixel_spacing,
    sampling_rate,
    sdims,
    spatialdims,
    spatial_order,
    spatial_size,
    spatial_indices,
    spatial_first,
    spatial_last,
    spatial_directions,
    timedim,
    times,
    time_first,
    time_last,
    time_step,
    width,
    depth,
    height

include("spatial.jl")
include("temporal.jl")

end # module

module SpatioTemporalTraits

using ArrayInterface
using ArrayInterface: known_first, known_last, known_step, static_first, static_last,
    static_length, static_step, StaticInt, One, OptionallyStaticUnitRange,
    OptionallyStaticStepRange, Zero
using AxisIndices
using Base: @pure, tail, OneTo
using CoordinateTransformations
using MappedArrays
using MappedArrays: AbstractMultiMappedArray, MultiMappedArray, ReadonlyMultiMappedArray
using Metadata
using NamedDims
using Rotations

export 
    assert_timedim_last,
    collapse,
    duration,
    each_channel,
    each_observation,
    each_time,
    has_channeldim,
    has_observationdim,
    has_timedim,
    iterate_channel,
    iterate_observation,
    iterate_time,
    lag,
    lead,
    nchannels,
    nobservations,
    ntimes,
    observationdim,
    observations,
    observation_indices,
    onset,
    sampling_rate,
    select_channel,
    select_observation,
    select_time,
    timedim,
    times,
    time_end,
    time_indices,
    time_step

include("iterators.jl")
include("names.jl")
include("observation.jl")
include("channel.jl")
include("time.jl")
include("spatial.jl")

end # module

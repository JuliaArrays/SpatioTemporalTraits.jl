module SpatioTemporalTraits

using ArrayInterface
using ArrayInterface: indices

using AxisIndices

using CoordinateTransformations

using MappedArrays
using MappedArrays: AbstractMultiMappedArray, MultiMappedArray, ReadonlyMultiMappedArray

using Metadata
using NamedDims
using Rotations


using Base: @pure

include("iterators.jl")
include("names.jl")
include("time.jl")
include("spatial.jl")

end # module

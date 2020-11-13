
using AxisIndices
using ArrayInterface
using SpatioTemporalTraits
using Test

val_spatialdims(x) = Val(spatialdims(x))

nda = NamedDimsArray{(:w,:x,:y,:z)}(rand(2, 2,2,2));

@test @inferred(val_spatialdims(nda)) === Val((1, 2, 3))

SpatioTemporalTraits.@not_spatial :y

@test @inferred(val_spatialdims(nda)) === Val((1, 2, 4))

SpatioTemporalTraits.@not_spatial :w

@test @inferred(val_spatialdims(nda)) === Val((2,4))

SpatioTemporalTraits.@not_spatial :x
@test @inferred(val_spatialdims(nda)) === Val((4,))

SpatioTemporalTraits.@not_spatial :z
@test @inferred(val_spatialdims(nda)) === Val(())

include("sliding_window.jl")
# TODO include("names.jl")
include("time.jl")

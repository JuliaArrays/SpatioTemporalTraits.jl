# SpatioTemporalTraits


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaarrays.github.io/SpatioTemporalTraits.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaarrays.github.io/SpatioTemporalTraits.jl/dev)
[![CI](https://github.com/JuliaArrays/SpatioTemporalTraits.jl/workflows/CI/badge.svg)](https://github.com/JuliaArrays/SpatioTemporalTraits.jl/actions?query=workflow%3ACI)
[![CI (Julia nightly)](https://github.com/JuliaArrays/SpatioTemporalTraits.jl/workflows/CI%20(Julia%20nightly)/badge.svg)](https://github.com/JuliaArrays/SpatioTemporalTraits.jl/actions?query=workflow%3A%22CI+%28Julia+nightly%29%22)
[![Build status](https://badge.buildkite.com/a2db252d92478e1d7196ee7454004efdfb6ab59496cbac91a2.svg?branch=master)](https://buildkite.com/julialang/SpatioTemporalTraits-dot-jl)
[![codecov](https://codecov.io/gh/JuliaArrays/SpatioTemporalTraits.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaArrays/SpatioTemporalTraits.jl)


# Introduction

`SpatioTemporalTraits` serves as a relatively low-level source of spatiotemporal traits, allowing other packages the opportunity to use a common interface for their unique types.
Many of these traits take advantage of [ArrayInterface.jl](https://github.com/JuliaArrays/ArrayInterface.jl) and have sensible defaults for subtypes of `AbstractArray` or any other type that uses `ArrayInterface`.


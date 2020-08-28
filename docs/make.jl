
using Documenter
using SpatioTemporalTraits

makedocs(;
    modules=[SpatioTemporalTraits],
    format=Documenter.HTML(),
    pages=["index.md"],
    repo="https://github.com/JuliaArrays/SpatioTemporalTraits.jl/blob/{commit}{path}#L{line}",
    sitename="SpatioTemporalTraits.jl",
)

deploydocs(
    repo = "github.com/JuliaArrays/SpatioTemporalTraits.jl.git",
)



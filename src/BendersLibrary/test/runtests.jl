using Test
using BendersLibrary

@testset "BendersLibrary.jl" begin

    include("test_types.jl")
    include("test_autodecomposeDisjunctive.jl")

end
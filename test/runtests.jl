using NonLinearReactionAdvectionDiffusionWithFrontData
using Test
using Documenter
using ForwardDiff
using LinearAlgebra

@testset "DocTests" begin
    #    doctest(NonLinearReactionAdvectionDiffusionWithFrontData)
end

include("utils.jl")

# TODO: they are broken due to changing signatures
#include("direct_check.jl")
#include("adjoint_check.jl")

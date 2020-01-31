using NonLinearReactionAdvectionDiffusionWithFrontData
using Test
using Documenter
using ForwardDiff
using LinearAlgebra

@testset "DocTests" begin
    #    doctest(NonLinearReactionAdvectionDiffusionWithFrontData)
end

# workaround GR warnings
ENV["GKSwstype"] = "100"

include("utils.jl")

include("direct_check.jl")
# TODO: they are broken due to changing signatures
#include("adjoint_check.jl")

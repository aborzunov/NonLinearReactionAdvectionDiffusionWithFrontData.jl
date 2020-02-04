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
include("adjoint_check.jl")

using NonLinearReactionAdvectionDiffusionWithFrontData
using Test
using Documenter
using ForwardDiff
using LinearAlgebra

@testset "DocTests" begin
    #    doctest(NonLinearReactionAdvectionDiffusionWithFrontData)
end

include("utils.jl")

include("direct_check.jl")
include("adjoint_check.jl")
include("gradient_check.jl")
include("degenerated_check.jl")

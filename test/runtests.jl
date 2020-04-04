using NonLinearReactionAdvectionDiffusionWithFrontData;
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using NonLinearReactionAdvectionDiffusionWithFrontData: apply_on_dynamic_mesh, Φ;

using Test
using Documenter
using ForwardDiff
using LinearAlgebra

@testset "DocTests" begin
    #    doctest(NonLinearReactionAdvectionDiffusionWithFrontData)
end

include("utils.jl")

include("direct_check.jl")
include("direct_nonuniform_check.jl")

#include("adjoint_check.jl")
#include("gradient_check.jl")
#include("degenerated_check.jl")

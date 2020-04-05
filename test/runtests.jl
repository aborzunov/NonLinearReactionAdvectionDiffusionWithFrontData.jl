using NonLinearReactionAdvectionDiffusionWithFrontData;
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using NonLinearReactionAdvectionDiffusionWithFrontData: apply_on_dynamic_mesh, Φ;

using Test
using Documenter
using ForwardDiff
using LinearAlgebra

include("utils.jl")

@time @testset "Якобиан прямой задачи на статической сетке.      " begin include("direct_jacobian.jl") end
@time @testset "Якобиан прямой задачи на динамической сетке.     " begin include("direct_jacobian_nonuniform.jl") end
@time @testset "Модельная невзяка для прямой задачи              " begin include("direct_check.jl") end

@time @testset "Якобиан сопряженной задачи на статической сетке. " begin include("adjoint_jacobian.jl") end
#include("adjoint_check.jl")

@time @testset "Градиент на точных данных                        " begin include("gradient_check.jl") end
#include("degenerated_check.jl")

using NonLinearReactionAdvectionDiffusionWithFrontData;

using Test
using Documenter
using ForwardDiff
using LinearAlgebra

include("utils.jl")

@time @testset "Якобиан прямой задачи на статической сетке.      " begin include("direct_jacobian.jl") end
@time @testset "Якобиан прямой задачи на динамической сетке.     " begin include("direct_jacobian_nonuniform.jl") end
@time @testset "Модельная невзяка для прямой задачи              " begin include("direct_check.jl") end

@time @testset "Якобиан сопряженной задачи на статической сетке. " begin include("adjoint_jacobian.jl") end
@time @testset "Модельная невязка для сопряженной задачи         " begin include("adjoint_check.jl") end

@time @testset "Неоднородность на точных данных                  " begin include("heterogeneity_check.jl") end
@time @testset "Вырожденые корни                                 " begin include("degenerated_check.jl") end
#
@time @testset "Градиент, функционал                             " begin include("functional.jl") end

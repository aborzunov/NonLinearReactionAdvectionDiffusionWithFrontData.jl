using NonLinearReactionAdvectionDiffusionWithFrontData;

using Test
using Documenter
using SafeTestsets

@time @safetestset "Вспомогательные                                  " begin include("utils.jl") end

@time @safetestset "Якобиан прямой задачи на статической сетке.      " begin include("direct_jacobian.jl") end
@time @safetestset "Якобиан прямой задачи на динамической сетке.     " begin include("direct_jacobian_nonuniform.jl") end
@time @safetestset "Модельная невзяка для прямой задачи              " begin include("direct_check.jl") end

@time @safetestset "Якобиан сопряженной задачи на статической сетке. " begin include("adjoint_jacobian.jl") end
@time @safetestset "Модельная невязка для сопряженной задачи         " begin include("adjoint_check.jl") end

@time @safetestset "Неоднородность на точных данных                  " begin include("heterogeneity_check.jl") end
@time @safetestset "Новая аппроксимация дельта-функции               " begin include("delta_check.jl") end
@time @safetestset "Вырожденые корни                                 " begin include("degenerated_check.jl") end

@time @safetestset "Градиент, функционал                             " begin include("functional.jl") end

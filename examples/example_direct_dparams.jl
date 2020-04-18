# ## Решение прямой задачи
#
# Этот файл будет использоваться в тестировании.
#
# Задаем параметры с помощью стандартной функции [`NonLinearReactionAdvectionDiffusionWithFrontData.dparams`](@ref)
using NonLinearReactionAdvectionDiffusionWithFrontData

a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = NonLinearReactionAdvectionDiffusionWithFrontData.dparams();
#
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
nothing #hide

# ## Генерация априорной информации
#
# Генерируем экспериментальные данные функцией [`NonLinearReactionAdvectionDiffusionWithFrontData.generate_obs_data`](@ref).
ϕl, ϕr, ϕ, f1_data, f2_data = NonLinearReactionAdvectionDiffusionWithFrontData.generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
nothing # hide

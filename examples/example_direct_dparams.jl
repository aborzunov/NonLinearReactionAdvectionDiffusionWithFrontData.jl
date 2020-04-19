# ## [Пример №2, краткий, стандартные параметры](@id de_2)
#
# В этом примере:
# - получаем стандартные входные параметры,
# - решается прямая задача,
# - генерируются экспериментальные данные.
# Только в отличии от [Подробного примера №1](@ref de_1)
# будем пользоваться функциями-сокращениями.
#
# Этот файл будет использоваться в тестировании.
#
# Задаем параметры с помощью стандартной функции [`NonLinearReactionAdvectionDiffusionWithFrontData.dparams`](@ref)
using NonLinearReactionAdvectionDiffusionWithFrontData

a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = NonLinearReactionAdvectionDiffusionWithFrontData.dparams();
#
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
nothing #hide

# Генерируем экспериментальные данные функцией [`NonLinearReactionAdvectionDiffusionWithFrontData.generate_obs_data`](@ref).
ϕl, ϕr, ϕ, f1_data, f2_data = NonLinearReactionAdvectionDiffusionWithFrontData.generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
nothing # hide

# # Эксперимент №2

# Здесь, выберем эмпирический параметр ``w`` в аппроксимации
# дельта-функции так, чтобы на каждом временном шаге в неоднородности были
# ненулевые элементы.

# ## Набор параметров
α       = 0.005;        # Параметр регуляризации
w       = 0.005;        # Эмпирический параметр регуляризации
S       = 5000;         # Количество итераций
β       = 0.001;        # Шаг минимизации
#
x_tp    = 0.1;          # Стартовое местоположение фронта
T_end   = 0.42;         # Регулируем конечное местоположение фронта
ε       = 0.03;         # Крутизна фронта

in("Travis", keys(ENV)) && S = 300

# -----------------------------------------------------------------------------
using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: heterogeneity_map;
using Plots; gr();

a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp=0.1,
                                                        ε = 0.03,
                                                        Nx = 200,
                                                        Mt = 250,
                                                        T_end = 0.42);
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения")

# Найдем начальное приближение
q_guess = initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ, α);
plot(Xₙ, qₙ, label="Истинное");
plot!(Xₙ, q_guess, label="Найденное")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# В следующей строке, мы записываем все параметры в latex формате,
# не забываем экранировать все спецсимволы.
# Дальше, мы отобразим все эти параметры на графике.
using Printf;
params = "\$w = $(w), \\varepsilon = $(ε), "*
"\\beta = $(β), f_1 \\in $(@sprintf("[%.2f, %.2f]", extrema(f1_data)...)), "*
"\\alpha = $(α)\$"
nothing #hide
# -----------------------------------------------------------------------------




# ## Выбор эмпирического параметра в аппроксимации дельта-функции
# -----------------------------------------------------------------------------
# Ненулевых элементов должно быть мало
hmap = heterogeneity_map(XX, N, Tₘ, M, u, f1_data, f2_data, w);
@info "Ненулевых элементов $(length( findall( x -> ! isapprox(x, 0), hmap)))"

heterogeneityP = heatmap(hmap, yflip=true)
# -----------------------------------------------------------------------------



# ## Старт с нулевого начального приближения
# -----------------------------------------------------------------------------
q₀ = zero(q_guess);
@time qs, Js, Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data, S = S, β = β, w = w)

# эскиз
noguessP = minimization_draft(qₙ, Qs, Xₙ, N, Js, zoom = true, annotate_string = params)
# -----------------------------------------------------------------------------




# ## Старт с найденного приближения
# -----------------------------------------------------------------------------
q₀ = q_guess;
@time qs, Js, Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data, S = S, β = β, w = w)

# эскиз
withguessP = minimization_draft(qₙ, Qs, Xₙ, N, Js,
                        zoom = true, annotate_string = params)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
savefig(directP, "direct2.png")
savefig(heterogeneityP, "heterogeneity2.png")
savefig(noguessP, "noguess2.png")
savefig(withguessP, "withguess2.png")
# -----------------------------------------------------------------------------

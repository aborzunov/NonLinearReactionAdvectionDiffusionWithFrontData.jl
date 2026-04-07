# # Эксперимент №2. Удачный пример.

# ## Набор параметров
α       = 0.004;        # Параметр регуляризации
w       = 0.0001;       # Эмпирический параметр регуляризации
S       = 150000;       # Количество итераций (верхняя граница; Armijo+tol_J обычно останавливают раньше)
β       = 1.0;          # Начальный шаг; line search Armijo подбирает фактический шаг
#
x_tp    = 0.05;         # Стартовое местоположение фронта
T_end   = 0.47;         # Регулируем конечное местоположение фронта
ε       = 10^(-1.5);    # Крутизна фронта
Nx      = 15000;          # Число интервалов по ``X``
Mt      = 30000;         # Число интервалов по ``T``

# -----------------------------------------------------------------------------
using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: heterogeneity_map;
using Serialization;
using Plots; gr();
using Dierckx;
# -----------------------------------------------------------------------------

# ## Решение на точных данных

# ### Решение прямой задачи для генерирования экспериментальной информации
# -----------------------------------------------------------------------------
@info "perfect_example: dparams() (крупная сетка Nx=$(Nx), Mt=$(Mt))..."
a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
                                                         ε = ε,
                                                         Nx = Nx,
                                                         Mt = Mt,
                                                         T_end = T_end);
@info "perfect_example: solve() #1 (Nx=$(Nx), Mt=$(Mt))..."
@time u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ; showProgress = true);
@info "perfect_example: generate_obs_data() #1..."
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ; showProgress = true);
directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения")
savefig(directP, "perfect_direct.png");
@info "perfect_example: сохранён perfect_direct.png"
# -----------------------------------------------------------------------------

# ### Найдем начальное приближение
# -----------------------------------------------------------------------------
@info "perfect_example: initial_guess()..."
q_guess = initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ, 0.005);
initial_guessP = plot(Xₙ, qₙ, xlabel = "X", ylabel = "q(x)", label="Истинное");
initial_guessP = plot!(Xₙ, q_guess, label="Найденное")
#savefig("initial_guess.svg");
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# В следующей строке, мы записываем все параметры в latex формате,
# не забываем экранировать все спецсимволы.
# Дальше, мы отобразим все эти параметры на графике.
using Printf;
params = "\$w = $(w), \\varepsilon = $(@sprintf("%.4f", ε)), "*
"\\beta = $(β), f_1 \\in $(@sprintf("[%.2f, %.2f]", extrema(f1_data)...)), "*
"\\alpha = $(α)\$"
nothing #hide
# -----------------------------------------------------------------------------


# ### Старт с найденного приближения на точных данных
# -----------------------------------------------------------------------------
@info "perfect_example: переходим на рабочую сетку Nx=500, Mt=1000..."
spl = Spline1D(Xₙ, q_guess);
Nx      = 500;          # Число интервалов по ``X``
Mt      = 1000;         # Число интервалов по ``T``
a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
                                                         ε = ε,
                                                         Nx = Nx,
                                                         Mt = Mt,
                                                         T_end = T_end);
@info "perfect_example: solve() #2 (Nx=$(Nx), Mt=$(Mt))..."
@time u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ; showProgress = true);
@info "perfect_example: generate_obs_data() #2..."
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ; showProgress = true);
q₀ = spl(Xₙ);

# -----------------------------------------------------------------------------
@info "perfect_example: minimize() с начальным приближением (S=$(S))..."
@time qs, Js, Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data,
                            S = S, β = β, w = w, showProgress = true,
                            linesearch = true, tol_J = 1e-10, tol_grad = 1e-8)
@info "perfect_example: minimize() с начальным приближением завершён (итераций: $(length(Js)), J_final = $(Js[end]))"
serialize("perfect_qs_withguess.jld", qs);
serialize("perfect_Js_withguess.jld", Js);
serialize("perfect_Qs_withguess.jld", Qs);
@info "perfect_example: данные сериализованы → perfect_{qs,Js,Qs}_withguess.jld"
# -----------------------------------------------------------------------------

# ### Эскиз процесса минимизации
# -----------------------------------------------------------------------------
@info "perfect_example: отрисовка результатов (withguess)..."
a, b, c = minimization_draft(qₙ, Qs, Xₙ, N, Js, zoom = true,
                              annotate_string = params, zoom_chunk = 17//18)
plot(a, size = (800, 800))
plot(b)
plot(c)
withguessP = plot(a, b, c, layout = (1,3), size = (2400, 800));
savefig(withguessP, "perfect_withguess.png")
@info "perfect_example: сохранён perfect_withguess.png"
nothing; #hide
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# ------------------ Старт с нулевого приближения -----------------------------
# -----------------------------------------------------------------------------
@info "perfect_example: minimize() с нулевого приближения (S=$(S))..."
q₀ = zero(Xₙ)
# -----------------------------------------------------------------------------
@time qs, Js, Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data,
                            S = S, β = β, w = w, showProgress = true,
                            linesearch = true, tol_J = 1e-10, tol_grad = 1e-8)
@info "perfect_example: minimize() с нулевого приближения завершён (итераций: $(length(Js)), J_final = $(Js[end]))"
serialize("perfect_qs_noguess.jld", qs);
serialize("perfect_Js_noguess.jld", Js);
serialize("perfect_Qs_noguess.jld", Qs);
@info "perfect_example: данные сериализованы → perfect_{qs,Js,Qs}_noguess.jld"
# -----------------------------------------------------------------------------

# ### Эскиз процесса минимизации
# -----------------------------------------------------------------------------
@info "perfect_example: отрисовка результатов (noguess)..."
a, b, c = minimization_draft(qₙ, Qs, Xₙ, N, Js, zoom = true,
                              annotate_string = params, zoom_chunk = 17//18)
plot(a, size = (800, 800))
plot(b)
plot(c)
noguessP = plot(a, b, c, layout = (1,3), size = (2400, 800));
savefig(noguessP, "perfect_noguess.png")
@info "perfect_example: сохранён perfect_noguess.png"
nothing; #hide
# -----------------------------------------------------------------------------

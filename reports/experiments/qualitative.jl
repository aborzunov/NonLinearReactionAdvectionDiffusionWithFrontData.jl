
# # Эксперимент №3 Зашумленные


# ## Набор параметров
α       = 0.004;        # Параметр регуляризации
w       = 0.0001;       # Эмпирический параметр регуляризации
S       = 3000;        # Количество итераций (верхняя граница; Armijo+tol_J обычно останавливают раньше)
β       = 1.0;          # Начальный шаг; line search Armijo подбирает фактический шаг
#
x_tp    = 0.05;         # Стартовое местоположение фронта
T_end   = 0.47;         # Регулируем конечное местоположение фронта
ε       = 10^(-1.5);    # Крутизна фронта
Nx      = 500;          # Число интервалов по ``X``
Mt      = 1000;         # Число интервалов по ``T``
δ       = 0.00008;      # Модуль величины помех

using Dates;
timestamp           = String(Dates.format(now(), "yyyy-mm-dd_HH-MM"));
wdir                = String("Expirement3___" * timestamp);
prefix              = String("Exp3_" * timestamp * "_");
n_noised_data       = prefix * "noised_data.png";
n_noised_velocity   = prefix * "noised_velocity.png";
n_Js_min_data       = prefix * "Js_noguess.jld";
n_Qs_min_data       = prefix * "Qs_noguess.jld";

using Markdown
info = """Вы зпустили эксперимент №3 (qualitative). Сценарий скрипта:
  - Создадим каталог `$(wdir)`.
  - Перейдем в него и все результаты поместим там.
  - Все файлы будем сохранять с префиксом `$(prefix)`
  - Сгенерируем экспериментальные данные на сетке `$(Nx+1) × $(Mt+1)`.
  - Зашумим экспериментальные данные f₁, f₂ гауссовским шумом с модулем `$(δ)`
    * Сохраним график зашумленных экспериментальных данных `$(n_noised_data)`
    * Сохраним график скорости фронта из зашумленных данных `$(n_noised_velocity)`
  - Запустим процесс минимизации с нулевого начального приближения на ЧИСТЫХ данных
    * Сериализуем данные и сохраним `$(n_Js_min_data)`, `$(n_Qs_min_data)`
"""
display(Markdown.parse(info))

@info "Переходим в $(wdir)"
run(`mkdir $(wdir)`)
run(`cd $(wdir)`)

# -----------------------------------------------------------------------------
using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: heterogeneity_map;
using Serialization;
using Plots; gr();
using Dierckx;
# -----------------------------------------------------------------------------


# ## Решение на точных данных

# ### Решение прямой задачи для генерирования экспериментальной информации ----------------------------------------------------------------------------- a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
# -----------------------------------------------------------------------------
a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
                                                         ε = ε,
                                                         Nx = Nx,
                                                         Mt = Mt,
                                                         T_end = T_end);
@info "solve() (Nx=$(Nx), Mt=$(Mt))..."
@time u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ; showProgress = true);
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения")
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

# ## Слабо зашумленные данные
# -----------------------------------------------------------------------------
using Random
using NonLinearReactionAdvectionDiffusionWithFrontData: front_velocity;
rng = MersenneTwister(1234)

# Создадим мелкий шум.
f1_data_noised = f1_data + randn(rng, length(f1_data)) * δ;
f2_data_noised = f2_data + randn(rng, length(f2_data)) * δ;
h = Xₙ[2] - Xₙ[1];
δₙ = sqrt(sum( [ (f1_data[n] - f1_data_noised[n])^2 * h for n in 1:N+1]));
@info δₙ

plot(f1_data_noised[1:5:end], Tₘ[1:5:end], ylims = (0, T_end * 1.1), label = "", title="Зашумленные \$f_1, δₙ = $(δₙ)\$", xlim = (0,1))
plot!(f1_data[1:5:end], Tₘ[1:5:end], label = "")

#plot!(Tₘ, f2_data_noised, label="Зашумленные \$f_2\$")
savefig(n_noised_data);


# Убедимся, что численное дифференцирование зашумленных функций --- это плохо.
v_f1 = front_velocity(f1_data_noised, Tₘ, M);
plot(Tₘ[5:end-5], v_f1[5:end-5], label="Скорость зашумленного фронта \$\\dfrac{df_1}{dt}\$")
savefig(n_noised_velocity);

# -----------------------------------------------------------------------------
# ------------------ Старт с нулевого приближения -----------------------------
# -----------------------------------------------------------------------------
q₀ = zero(Xₙ)
# -----------------------------------------------------------------------------
@time qs, noised_noguess_Js, noised_noguess_Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data,
                            S = S, β = β, w = w, showProgress = true,
                            linesearch = true, tol_J = 1e-10, tol_grad = 1e-8)
serialize(n_Js_min_data, noised_noguess_Js);
serialize(n_Qs_min_data, noised_noguess_Qs);
# -----------------------------------------------------------------------------

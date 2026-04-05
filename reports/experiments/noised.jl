# # Эксперимент №2 Зашумленные


# ## Набор параметров
α       = 0.004;        # Параметр регуляризации
w       = 0.0001;       # Эмпирический параметр регуляризации
S       = 50000;        # Количество итераций
β       = 0.001;        # Шаг минимизации
#
x_tp    = 0.05;         # Стартовое местоположение фронта
T_end   = 0.47;         # Регулируем конечное местоположение фронта
ε       = 10^(-1.5);    # Крутизна фронта
Nx      = 500;          # Число интервалов по ``X``
Mt      = 1000;         # Число интервалов по ``T``
δ       = 0.00008;      # Модуль величины помех

using Dates;
timestamp           = String(Dates.format(now(), "yyyy-mm-dd_HH-MM"));
wdir                = String("Expirement2___" * timestamp);
prefix              = String("Exp2_" * timestamp * "_");
n_direct_draft      = prefix * "direct_draft.png";
n_direct_gif        = prefix * "direct.gif";
n_noised_data       = prefix * "noised_data.png";
n_noised_velocity   = prefix * "noised_velocity.png";
n_initial_guess     = prefix * "initial_guess.png";
n_Js_min_data_guess = prefix * "Js_withguess.jld";
n_Qs_min_data_guess = prefix * "Qs_withguess.jld";
n_Js_min_data       = prefix * "Js_noguess.jld";
n_Qs_min_data       = prefix * "Qs_noguess.jld";
n_guess_draft       = prefix * "minimization_draft_withguess.png"
n_noguess_draft     = prefix * "minimization_draft_noguess.png"
n_minim_gif_guess   = prefix * "minimization_withguess.gif";
n_minim_gif_noguess = prefix * "minimization_noguess.gif";
n_minim_avi_guess   = prefix * "minimization_withguess.avi";
n_minim_avi_noguess = prefix * "minimization_noguess.avi";
n_minim_mp4_guess   = prefix * "minimization_withguess.mp4";
n_minim_mp4_noguess = prefix * "minimization_noguess.mp4";

using Markdown
info = """Вы зпустили эксперимент №2. Сценарий скрипта:
  - Создадим каталог `$(wdir)`.
  - Перейдем в него и все результаты поместим там.
  - Все файлы будем сохранять с префиксом `$(prefix)`
  - Сгенерируем экспериментальные данные на сетке `$(Nx+1) × $(Mt+1)`
    * Сохраним эскиз прямого решения `$(n_direct_draft)`
    * Сохраним анимацию решения прямой задачи, с помощью которой генерировали
        экспериментальные даные `$(n_direct_gif)`.
  - Зашумим экспериментальные данные f₁, f₂ равномерным распределение с модулем `$(δ)`
    * Сохраним график зашумленных экспериментальных данных `$(n_noised_data)`
    * Сохраним график скорости фронта из зашумленных данных `$(n_noised_velocity)`
  - Найдем начальное приближение с параметром регуляризации `$(α)`
    * Сохраним начальное приближение `$(n_initial_guess)`
  - Начнем процесс минимизации используя найденное приближение
    * Сериализуем данные и сохраним `$(n_Js_min_data_guess)`, `$(n_Qs_min_data_guess)`
    * Сохраним Черновик процесса минимизации `$(n_guess_draft)`
    * Сохраним Анимацию процесса минимизации `$(n_minim_gif_guess)`
    * Сохраним mp4 Анимации процесса минимизации `$(n_minim_mp4_guess)`
    * Сохраним mp4 Анимации процесса минимизации `$(n_minim_avi_guess)`
  - Начнем процесс минимизации с нулевого начального приближения
    * Сериализуем данные и сохраним `$(n_Js_min_data)`, `$(n_Qs_min_data)`
    * Сохраним Черновик процесса минимизации `$(n_noguess_draft)`
    * Сохраним Анимацию процесса минимизации `$(n_minim_gif_noguess)`
    * Сохраним mp4 Анимации процесса минимизации `$(n_minim_mp4_noguess)`
    * Сохраним mp4 Анимации процесса минимизации `$(n_minim_avi_noguess)`
"""
display(Markdown.parse(info))

@info "Переходим в $(wdir)"
run(`mkdir $(wdir)`)
run(`cd $(wdir)`)

# -----------------------------------------------------------------------------
using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: heterogeneity_map;
using Serialization;
using Plots; pyplot();
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
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения")
make_gif(u, Xₙ, Tₘ, ϕl, ϕr, ϕ, f1_data, f2_data, name = n_direct_gif);
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
h = Xₙ[2] - Xₙ[1]
δₙ = sqrt(sum( [ (f1_data[n] - f1_data_noised[n])^2 * h for n in 1:N+1]));
@info δₙ
plot(f1_data_noised[1:5:end], Tₘ[1:5:end], ylims = (0, T_end * 1.1), label = "", title="Зашумленные \$f_1, δₙ = $(δₙ)\$", xlim = (0,1))
plot!(f1_data[1:5:end], Tₘ[1:5:end], label = "")
#plot!(Tₘ, f2_data_noised, label="Зашумленные \$f_2\$")

savefig(n_noised_data);


# Убедимся, что численное дифференцирование зашумленных функций --- это плохо.
v_f1 = front_velocity(f1_data_noised, Tₘ, M);
plot(Tₘ, v_f1[5:end-5], label="Скорость зашумленного фронта \$\\dfrac{df_1}{dt}\$")
savefig(n_noised_velocity);

# ### Найдем начальное приближение
# -----------------------------------------------------------------------------
q_guess = initial_guess(f1_data_noised, Xₙ, N, Tₘ, M, ulₘ, urₘ, 0.005);
initial_guessP = plot(Xₙ, qₙ, xlabel = "X", ylabel = "q(x)", label="Истинное");
initial_guessP = plot!(Xₙ, q_guess, label="Найденное")
savefig(n_initial_guess);
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
spl = Spline1D(Xₙ, q_guess);
Nx      = 500;          # Число интервалов по ``X``
Mt      = 1000;         # Число интервалов по ``T``
a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
                                                         ε = ε,
                                                         Nx = Nx,
                                                         Mt = Mt,
                                                         T_end = T_end);
q₀ = spl(Xₙ);

# ### Старт с найденного приближения на точных данных
# -----------------------------------------------------------------------------
@time qs, noised_guess_Js, noised_guess_Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data_noised, f2_data_noised,
                            S = S, β = β, w = w, showProgress = true)
serialize(n_Js_min_data_guess, noised_guess_Js);
serialize(n_Qs_min_data_guess, noised_guess_Qs);
# -----------------------------------------------------------------------------

# ### Отрисовка результатов минимизации
# -----------------------------------------------------------------------------
a, b, c = minimization_draft(qₙ, noised_guess_Qs, Xₙ, N, noised_guess_Js, zoom = true,
                              annotate_string = params, zoom_chunk = 17//18)
plot(a, size = (800, 800))
plot(b)
plot(c)
withguessP = plot(a, b, c, layout = (3,1), size = (1000, 2400));
savefig(withguessP, n_guess_draft)
nothing; #hide
# -----------------------------------------------------------------------------
frames = collect(1:div(S, 500):S);
make_minimzation_gif(
    noised_guess_Js, noised_guess_Qs, qₙ, Xₙ,
    name = n_minim_gif_guess,
    frames_to_write = frames,
    β = β
)
run(`gif2mp4 $(n_minim_gif_guess) $(n_minim_mp4_guess)`)
run(`ffmpeg -i  $(n_minim_gif_guess) $(n_minim_avi_guess) -vcodec copy -acodec copy`)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# ------------------ Старт с нулевого приближения -----------------------------
# -----------------------------------------------------------------------------
q₀ = zero(Xₙ)
# -----------------------------------------------------------------------------
@time qs, noised_noguess_Js, noised_noguess_Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data,
                            S = S, β = β, w = w, showProgress = true)
serialize(n_Js_min_data, noised_noguess_Js);
serialize(n_Qs_min_data, noised_noguess_Qs);
# -----------------------------------------------------------------------------

# ### Эскиз процесса минимизации
# -----------------------------------------------------------------------------
a, b, c = minimization_draft(qₙ, noised_noguess_Qs, Xₙ, N, noised_noguess_Js, zoom = true,
                              annotate_string = params, zoom_chunk = 17//18)
plot(a, size = (800, 800))
plot(b)
plot(c)
withguessP = plot(a, b, c, layout = (3,1), size = (1000, 2400));
savefig(withguessP, n_noguess_draft)
nothing; #hide
frames = collect(1:div(S, 500):S);
make_minimzation_gif(
    noised_noguess_Js, noised_noguess_Qs, qₙ, Xₙ,
    name = n_minim_gif_noguess,
    frames_to_write = frames,
    β = β
)
run(`gif2mp4 $(n_minim_gif_noguess) $(n_minim_mp4_noguess)`)
run(`ffmpeg -i  $(n_minim_gif_noguess) $(n_minim_avi_noguess) -vcodec copy -acodec copy`)
# -----------------------------------------------------------------------------

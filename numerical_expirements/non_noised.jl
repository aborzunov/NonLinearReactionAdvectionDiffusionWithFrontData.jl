# # Эксперимент №1 Не зашумленные

# ## Набор параметров
α       = 0.004;        # Параметр регуляризации
w       = 0.0001;       # Эмпирический параметр регуляризации
S       =  5000;        # Количество итераций
β       = 0.001;        # Шаг минимизации
#
x_tp    = 0.05;         # Стартовое местоположение фронта
T_end   = 0.47;         # Регулируем конечное местоположение фронта
ε       = 10^(-1.5);    # Крутизна фронта
Nx      = 15000;          # Число интервалов по ``X``
Mt      = 30000;         # Число интервалов по ``T``

using Dates;
timestamp           = String(Dates.format(now(), "yyyy-mm-dd_HH-MM"));
wdir                = String("Expirement1___" * timestamp);
prefix              = String("Exp1_" * timestamp * "_");
n_direct_draft      = prefix * "direct_draft.png";
n_direct_gif        = prefix * "direct.gif";
n_data              = prefix * "data.png";
n_velocity          = prefix * "velocity.png";
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
info = """Вы зпустили эксперимент №1. Сценарий скрипта:
  - Создадим каталог `$(wdir)`.
  - Перейдем в него и все результаты поместим там.
  - Все файлы будем сохранять с префиксом `$(prefix)`
  - Сгенерируем экспериментальные данные на сетке `$(Nx+1) × $(Mt+1)`
    * Сохраним эскиз прямого решения `$(n_direct_draft)`
    * Сохраним анимацию решения прямой задачи, с помощью которой генерировали
        экспериментальные даные `$(n_direct_gif)`.
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

# ### Решение прямой задачи для генерирования экспериментальной информации
@info "Решение прямой задачи для генерирования экспериментальной информации"
# -----------------------------------------------------------------------------
a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
                                                         ε = ε,
                                                         Nx = Nx,
                                                         Mt = Mt,
                                                         T_end = T_end);
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения")
savefig(directP, n_direct_draft);
# -----------------------------------------------------------------------------

# ### Найдем начальное приближение
@info "Найдем начальное приближение"
# -----------------------------------------------------------------------------
q_guess = initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ, 0.005);
initial_guessP = plot(Xₙ, qₙ, xlabel = "X", ylabel = "q(x)", label="Истинное");
initial_guessP = plot!(Xₙ, q_guess, label="Найденное")
savefig(n_initial_guess);
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


# -----------------------------------------------------------------------------
spl = Spline1D(Xₙ, q_guess);
Nx      = 500;          # Число интервалов по ``X``
Mt      = 1000;         # Число интервалов по ``T``
a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
                                                         ε = ε,
                                                         Nx = Nx,
                                                         Mt = Mt,
                                                         T_end = T_end);
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
q₀ = spl(Xₙ);

# ### Старт с найденного приближения на точных данных
# -----------------------------------------------------------------------------
@time qs, guess_Js, guess_Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data,
                            S = S, β = β, w = w, showProgress = true)
serialize(n_Js_min_data_guess, guess_Js);
serialize(n_Qs_min_data_guess, guess_Qs);
# -----------------------------------------------------------------------------

# ### Отрисовка результатов минимизации
# -----------------------------------------------------------------------------
a, b, c = minimization_draft(qₙ, guess_Qs, Xₙ, N, guess_Js, zoom = true,
                              annotate_string = params, zoom_chunk = 17//18)
plot(a, size = (800, 800))
plot(b)
plot(c)
withguessP = plot(a, b, c, layout = (3,1), size = (1000, 2400));
savefig(withguessP, n_guess_draft)
nothing; #hide
# -----------------------------------------------------------------------------
frames = collect(1:div(S, 50):S);
make_minimzation_gif(
    guess_Js, guess_Qs, qₙ, Xₙ,
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
@time qs, noguess_Js, noguess_Qs = minimize(q₀, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data,
                            S = S, β = β, w = w, showProgress = true)

serialize(n_Js_min_data, noguess_Js);
serialize(n_Qs_min_data, noguess_Qs);
# -----------------------------------------------------------------------------

# ### Эскиз процесса минимизации
# -----------------------------------------------------------------------------
a, b, c = minimization_draft(qₙ, noguess_Qs, Xₙ, N, noguess_Js, zoom = true,
                              annotate_string = params, zoom_chunk = 17//18)
plot(a, size = (800, 800))
plot(b)
plot(c)
withguessP = plot(a, b, c, layout = (3,1), size = (1000, 2400));
savefig(withguessP, n_noguess_draft)
nothing; #hide
frames = collect(1:div(S, 50):S);
make_minimzation_gif(
    noguess_Js, noguess_Qs, qₙ, Xₙ,
    name = n_minim_gif_noguess,
    frames_to_write = frames,
    β = β
)
run(`gif2mp4 $(n_minim_gif_noguess) $(n_minim_mp4_noguess)`)
run(`ffmpeg -i  $(n_minim_gif_noguess) $(n_minim_avi_noguess) -vcodec copy -acodec copy`)
# -----------------------------------------------------------------------------

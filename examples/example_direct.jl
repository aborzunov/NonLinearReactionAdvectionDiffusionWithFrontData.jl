#' ## Решение прямой задачи
using NonLinearReactionAdvectionDiffusionWithFrontData

u_l(t) = -8 #+ cos(2*π * t);
u_r(t) =  4 #+ (1 + sin(2*π * t));
q(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
                                # задаче необходимо определить, но при генерации априорной
                                # информации мы задаем некоторый коэффициент, который,
                                # собственно, после имея априорную информацию и будем определять.
ε = 0.2;                        # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 0.28;                # Область по T
N, M = 50, 80;                  # Кол-во разбиений по X, T
h = (b-a)/N;                    # шаг по X
τ = (T-t₀)/M;                   # шаг по T
Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
qₙ =      q.(Xₙ);               # Сеточные значения коэффициента линейного усиления
ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
y₀ = u_init.(Xₙ);               # Начальные условия
nothing #hide

#' ## Нахождение решения
#' Все массивы передаются внутрь функции `solve` полностью, вместе с граничными точками.
#' Внутри, они локально модифицируются, и на вход [`directRP`](@ref), [`∂directRP_∂y`](@ref)
#' подаются без крайних точек.
u = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
nothing #hide

#' Запись gif только решения
make_gif(u, Xₙ, Tₘ; frame_skip = div(M,50), frames_to_write=M, name="example_direct.gif")
nothing #hide

#' ![solution gif](example_direct.gif)

#' ## Генерация априорной информации
# Вырожденные корни
# TODO: FIX phidetermination
y = u[:,1] # TMP FIXURE
ϕl = phidetermination(qₙ, y, ulₘ, Xₙ, N::Int);
# Разворачиваем сетку по Х, а после — решение
ϕr = phidetermination(qₙ, y, urₘ, Xₙ[end:-1:1], N::Int);
ϕr = ϕr[end:-1:1];
# Полуразность вырожденных корней
ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N);
# Положение переходного слоя
f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);
# Значение функции на переходном слое
f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);

# Можно нарисовать пятый шаг по времени
make_plot(u, Xₙ, Tₘ, 5, ϕl, ϕr, f1, f2)

#' Запись **только** mp4 вместе с вырожденными корнями
make_gif(u, Xₙ, Tₘ, ϕl, ϕr, f1, f2; frame_skip = div(M,50), frames_to_write=M,
         convert2mp4 = true, name="example_direct_with_f1_f2.gif")
nothing #hide

#' ![solution mp4 with degenerated](example_direct_with_f1_f2.mp4)

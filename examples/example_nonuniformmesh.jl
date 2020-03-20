#' ## Решение прямой задачи
using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: meshformation;

u_l(t) = -8;                    # ГУ
u_r(t) =  4;                    #
q(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
                                # задаче необходимо определить, но при генерации априорной
                                # информации мы задаем некоторый коэффициент, который,
                                # собственно, после имея априорную информацию и будем определять.
ε = 0.01;                        # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 0.28;                # Область по T
x_tp = 0.25;
N, M = 20, 400;                  # Кол-во разбиений по X, T
h = (b-a)/N;                    # шаг по X
τ = (T-t₀)/M;                   # шаг по T
Xₙ = meshformation(a, b, x_tp, ε, N, 1.0, 8.0);   # Сетка по Х
N = 11N;
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
qₙ =      q.(Xₙ);               # Сеточные значения коэффициента линейного усиления
ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
y₀ = u_init.(Xₙ, ε=ε, x_tp = x_tp);               # Начальные условия
nothing #hide

#' Все массивы передаются внутрь функции `solve` полностью, вместе с граничными точками.
#'
#' Внутри, они локально модифицируются, и на вход
#' [`NonLinearReactionAdvectionDiffusionWithFrontData.directRP`](@ref),
#' [`NonLinearReactionAdvectionDiffusionWithFrontData.∂directRP_∂y`](@ref)
#' подаются без крайних точек.
u = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
nothing #hide


#' ## Генерация априорной информации
ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
ϕr = phidetermination(qₙ, urₘ, reverse(Xₙ), N, Tₘ, M);                      # Нужно подать инвертированную сетку
ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N, M);       # Серединный корень
f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);   # Положение переходного слоя
f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);  # Значение функции на переходном слое
nothing # hide

#' ## Визуализация

#' Нарисовать пятый шаг по времени
m = 2;
make_plot(u, Xₙ, Tₘ, m, y_lim=extrema(u[:, m]*1.05))

#' Запись gif одного только решения
@info "$( splitdir(@__FILE__)[2] ) Рисует решение прямой задачи."
make_gif(u, Xₙ, Tₘ; name="example_direct.gif", frames_to_write=collect( 1:20 ), convert2mp4 = true)

#' Запись **только** mp4 вместе с вырожденными корнями
@info "$( splitdir(@__FILE__)[2] ) Рисует решение с вырожденными корнями и информации о переходном слое."
make_gif(u, Xₙ, Tₘ, ϕl, ϕr; convert2mp4 = true, name="example_direct_with_f1_f2.gif")

#' ## Решение прямой задачи
using NonLinearReactionAdvectionDiffusionWithFrontData

u_l(t) = -8 + 4*sin(2*π / T * t);# Прямая задача может быть с неоднородными ГУ
u_r(t) =  4 + sin(-2*π / T * t);# Но в дальнейшем, будем использовать только однородные.
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

#' ## Визуализация

#' Нарисовать пятый шаг по времени
make_plot(u, Xₙ, Tₘ, 5, ϕl, ϕr, f1, f2)

#' Запись gif одного только решения
@info "$( splitdir(@__FILE__)[2] ) Рисует решение прямой задачи."
make_gif(u, Xₙ, Tₘ; name="example_direct.gif")

#' Запись **только** mp4 вместе с вырожденными корнями
@info "$( splitdir(@__FILE__)[2] ) Рисует решение с вырожденными корнями и информации о переходном слое."
make_gif(u, Xₙ, Tₘ, ϕl, ϕr, f1, f2; convert2mp4 = true, name="example_direct_with_f1_f2.gif")

#' ## Решение прямой задачи
#' Сопряженная задача определяется для сеточного решения некоторой прямой задачи,
#' и некоторых априорных данных.
using NonLinearReactionAdvectionDiffusionWithFrontData

#' Сначала, сгенирируем априорные данные, на увеличенном количестве узлов.
u_l(t) = -8 #+ cos(2*π * t);
u_r(t) =  4 #+ (1 + sin(2*π * t));
q(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
                                # задаче необходимо определить, но при генерации априорной
                                # информации мы задаем некоторый коэффициент, который,
                                # собственно, после имея априорную информацию и будем определять.
ε = 0.2;                        # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 0.28;                # Область по T
N, M = 250, 400;                # Увеличенное Кол-во разбиений по X, T
h = (b-a)/N;                    # шаг по X
τ = (T-t₀)/M;                   # шаг по T
Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
qₙ =      q.(Xₙ);               # Сеточные значения коэффициента линейного усиления
ulₙ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₙ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
y₀ = u_init.(Xₙ);               # Начальные условия
nothing #hide

u = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₙ, urₙ, qₙ);

# Вырожденные корни
# TODO: FIX PHIDETERMINATION
y = u[:,1]
ϕl = phidetermination(qₙ, y, ulₙ, Xₙ, N::Int);
# Разворачиваем сетку по Х, а после — решение
ϕr = phidetermination(qₙ, y, urₙ, Xₙ[end:-1:1], N::Int);
ϕr = ϕr[end:-1:1];
# Полуразность вырожденных корней
ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N);
# Положение переходного слоя
f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);
# Значение функции на переходном слое
f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);

#' ## Решение сопряженной задачи

Uₙₘ = u; # Сохраним старую матрицу
y₀ = [0.0 for i in 1:N+1];
ψl = [0.0 for i in 1:M+1];
ψr = [0.0 for i in 1:M+1];

# tmp vars
# X = Xₙ[2:N];
# q_ = qₙ[2:N];
# y = y₀[2:N];
# T = Tₘ[end:-1:1];
# U = Uₙₘ[2:N, :];
# NonLinearReactionAdvectionDiffusionWithFrontData.adjointRP(y, 1, X, N, T, M, ε, ψl, ψr, q_, U, f1, f2)
# NonLinearReactionAdvectionDiffusionWithFrontData.∂adjointRP_∂y(y, 1, X, N, T, M, ε, ψl, ψr, q_, U, f1, f2)

ψ = solve_adjoint(y₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, Uₙₘ, f1, f2)

# На отрисовку, решение сопряженной задачи передадим в инвертированном времени.
make_gif(ψ[:,end:-1:1], Xₙ, Tₘ[end:-1:1]; frame_skip = div(M, 350), frames_to_write=81, name="adjoint.gif", convert2mp4 = true)

#' ![](adjoint.mp4)

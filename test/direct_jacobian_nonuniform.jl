using NonLinearReactionAdvectionDiffusionWithFrontData;
using NonLinearReactionAdvectionDiffusionWithFrontData: shishkin_mesh;
using LinearAlgebra;

u_l(t) = -8;                    # ГУ
u_r(t) =  4;                    #
qf(x) = 4*sin(3 * π * x);       # Коэффициент линейного усиления, который в обратной
# задаче необходимо определить, но при генерации априорной
# информации мы задаем некоторый коэффициент, который,
# собственно, после имея априорную информацию и будем определять.
ε = 0.01;                       # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 0.50;                # Область по T
x_tp = 0.12;                    # Положение переходного слоя
M = 500;                        # Кол-во разбиений по X, T
τ = (T-t₀)/M;                   # шаг по T
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
nothing #hide

#' Новое поведение будет реализоваться с помощью передачи функции создания сетки внутрь `solve` в качестве аргумента.
#' Эта функция должна принимать только один аргумент — положение переходного слоя и формировать соответствующую сетку.
#' В пакете есть формирование кусочно-равномерной сетки со сгущениями на границе и на переходном слое.
#' Создадим замыкание этой функции, которое будет иметь нужную сигнатуру.
mshfrm(x_tp) = shishkin_mesh(a, b, x_tp, ε, 40, 0.5, 1.0, 1.0, 0.25);
Xₙ = mshfrm(x_tp); ;                            # Сетка по Х
N = length(Xₙ) - 1                              # Примем за N длину сетки, что получилась в итоге.
qₙ =    qf.(Xₙ);                                # Сеточные значения коэффициента линейного усиления
y₀ = u_init.(Xₙ, ε=ε, x_tp = x_tp);             # Начальные условия

# Подготовим массивы, выбросив граничные точки
# ведь тестируемая функция — для внутреннего использования
qq = qₙ[2:N];
y = y₀[2:N];

dl, d, du = NonLinearReactionAdvectionDiffusionWithFrontData.DRP_y(y, 1, Xₙ, N, ε, ulₘ, urₘ, qq)

@test length(dl) == N - 2
@test length(d)  == N - 1
@test length(du) == N - 2

@test count(x -> x == 0, dl) == 0
@test count(x -> x == 0, d)  == 0
@test count(x -> x == 0, du) == 0

tdjac = NonLinearReactionAdvectionDiffusionWithFrontData.∂DRP_∂y(y, 1, Xₙ, N, ε, ulₘ, urₘ, qq)
@test typeof(tdjac) <: Tridiagonal
@test tdjac == Tridiagonal( dl, d, du )

# Сравним с якобианом автодифференцирования
jac = NonLinearReactionAdvectionDiffusionWithFrontData.∂directRP_∂y(y, 1, Xₙ, N, ε, ulₘ, urₘ, qq)
@test isapprox(tdjac, jac)

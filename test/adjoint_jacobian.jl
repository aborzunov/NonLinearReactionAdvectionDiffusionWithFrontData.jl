using NonLinearReactionAdvectionDiffusionWithFrontData;
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using NonLinearReactionAdvectionDiffusionWithFrontData: ARP_y, ∂ARP_∂y;
using LinearAlgebra;


#' Сначала, сгенирируем априорные данные, на увеличенном количестве узлов.
u_l(t) = -8 #+ cos(2*π * t);
u_r(t) =  4 #+ (1 + sin(2*π * t));
qf(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
                                # задаче необходимо определить, но при генерации априорной
                                # информации мы задаем некоторый коэффициент, который,
                                # собственно, после имея априорную информацию и будем определять.
ε = 0.2;                        # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0.0, 0.28;              # Область по T
N, M = 50, 80;                  # Увеличенное Кол-во разбиений по X, T
h = (b-a)/N;                    # шаг по X
τ = (T-t₀)/M;                   # шаг по T
Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
qₙ =      qf.(Xₙ);               # Сеточные значения коэффициента линейного усиления
ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
u₀ = u_init.(Xₙ);               # Начальные условия
#
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
nothing #hide

#' ## Генерация априорной информации
ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
ϕr = phidetermination(reverse(qₙ), urₘ, reverse(Xₙ), N, Tₘ, M);             # Нужно подать инвертированную сетку
ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
ϕ = Φ(ϕl, ϕr, N, M);                                                        # Серединный корень
f1_data = f1(ϕ, u, Xₙ, N, M);                                               # Положение переходного слоя
f2_data = f2(f1_data, u, Xₙ, N, M);                                         # Значение функции на переходном слое

# Подготовим массивы, выбросив граничные точки
# ведь тестируемая функция — для внутреннего использования
# И развернем массивы по времени
qq = qₙ[2:N];
y = u₀[2:N];
U = u[2:N, :];
U = reverse(U, dims=2);
f1_data = reverse(f1_data);
f2_data = reverse(f2_data);
Tₘ = reverse(Tₘ);

dl, d, du = NonLinearReactionAdvectionDiffusionWithFrontData.ARP_y(y, 1, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qq, U, f1_data, f2_data)

@test length(dl) == N - 2
@test length(d)  == N - 1
@test length(du) == N - 2

@test count(x -> x == 0, dl) == 0
@test count(x -> x == 0, d)  == 0
@test count(x -> x == 0, du) == 0

tdjac = NonLinearReactionAdvectionDiffusionWithFrontData.∂ARP_∂y(y, 1, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qq, U, f1_data, f2_data)
@test typeof(tdjac) <: Tridiagonal
@test tdjac == Tridiagonal( dl, d, du )

# Сравним с якобианом автодифференцирования
# Передадим **любой** априорный параметр для аппроксимации дельта-функции
# он никак не участвует в вычислении якобиана,
# но при использовании якобиана через автодифференцирование, требуется, чтобы он вызывался ТОЧНО так же
# как и функция, якобиан от которой вычисляется.
jac = NonLinearReactionAdvectionDiffusionWithFrontData.∂adjointRP_∂y(y, 1, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qq, U, f1_data, f2_data, 0.0001)
@test isapprox(tdjac, jac)

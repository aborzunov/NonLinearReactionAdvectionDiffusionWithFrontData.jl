
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


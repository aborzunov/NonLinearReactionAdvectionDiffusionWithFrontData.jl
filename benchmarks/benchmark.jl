using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: DRP_y, ∂DRP_∂y, directRP, ∂directRP_∂y;
using NonLinearReactionAdvectionDiffusionWithFrontData: ARP_y, ∂ARP_∂y, adjointRP, ∂adjointRP_∂y;
using BenchmarkTools
using LinearAlgebra


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
u₀ = u_init.(Xₙ);               # Начальные условия
nothing #hide

# Подготовим массивы, выбросив граничные точки
# ведь тестируемая функция — для внутреннего использования
X = Xₙ[2:N];
qq = qₙ[2:N];
y = u₀[2:N];

dl, d, du = DRP_y(y, 1, X, N, ε, ulₘ, urₘ, qq)

b = @benchmark DRP_y($y, 1, $X, $N, $ε, $ulₘ, $urₘ, $qq);
@info "Сформировать 3 вектора диагоналей якобиана прямой задачи:" b

b = @benchmark ∂DRP_∂y($y, 1, $X, $N, $ε, $ulₘ, $urₘ, $qq);
@info "Трехдиагональная матрица `∂DRP_∂y`:" b

b = @benchmark ∂directRP_∂y($y, 1, $X, $N, $ε, $ulₘ, $urₘ, $qq);
@info "Якобиан прямой задачи с помощью автоматического дифференциарования:" b

b = @benchmark solve($u₀, $Xₙ, $N, $Tₘ, $M, $ε, $ulₘ, $urₘ, $qₙ, $directRP, $∂DRP_∂y);
@info "Решение прямой задачи с трехдиагональным якобианом" b

b = @benchmark solve($u₀, $Xₙ, $N, $Tₘ, $M, $ε, $ulₘ, $urₘ, $qₙ, $directRP, $∂directRP_∂y);
@info "Решение прямой задачи с якобианом афтодифференцирования" b

# Готовимся к решению сопряженной задачи
u = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
#' ## Генерация априорной информации
ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
ϕr = phidetermination(qₙ, urₘ, reverse(Xₙ), N, Tₘ, M);                      # Нужно подать инвертированную сетку
ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N, M);       # Серединный корень
f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);   # Положение переходного слоя
f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);  # Значение функции на переходном слое

# Подготовим массивы, выбросив граничные точки
# ведь тестируемая функция — для внутреннего использования
U = u[2:N, :];

dl, d, du = ARP_y(y, 1, X, N, Tₘ, M, ε, ulₘ, urₘ, qq, U, f1, f2)

b = @benchmark ARP_y($y, 1, $X, $N, $Tₘ, $M, $ε, $ulₘ, $urₘ, $qq, $U, $f1, $f2);
@info "Сформировать 3 вектора диагоналей якобиана сопряженной задачи:" b

b = @benchmark ∂ARP_∂y($y, 1, $X, $N, $Tₘ, $M, $ε, $ulₘ, $urₘ, $qq, $U, $f1, $f2);
@info "Трехдиагональная матрица `∂ARP_∂y`:" b

# Здесь значения якобиана вычисляются неправильно, ведь нам нужно еще
# развернуть вектора f1, f2 и U(по второй координате).
# Помните — Функция для внутреннего использования, поэтому сейчас мы
# должны сделать это сами.
b = @benchmark ∂adjointRP_∂y($y, 1, $X, $N, $reverse(Tₘ), $M, $ε, $ulₘ, $urₘ, $qq, $U, $f1, $f2);
@info "Якобиан сопряженной задачи с помощью автоматического дифференциарования:" b

b = @benchmark solve_adjoint($u₀, $Xₙ, $N, $Tₘ, $M, $ε, $ulₘ, $urₘ, $qₙ, $u, $f1, $f2, $adjointRP, $∂ARP_∂y);
@info "Решение сопряженной задачи с трехдиагональным якобианом" b

b = @benchmark solve_adjoint($u₀, $Xₙ, $N, $Tₘ, $M, $ε, $ulₘ, $urₘ, $qₙ, $u, $f1, $f2, $adjointRP, $∂adjointRP_∂y);
@info "Решение сопряженной задачи с якобианом афтодифференцирования" b

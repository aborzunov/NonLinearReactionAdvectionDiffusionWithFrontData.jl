#' ## Генерация априорной информации

using NonLinearReactionAdvectionDiffusionWithFrontData
using LaTeXStrings
using Plots; gr()
using ProgressMeter
#
u_l(t) = -8 #+ cos(2*π * t);
u_r(t) =  4 #+ (1 + sin(2*π * t));
q(x) = sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
                                # задаче необходимо определить, но при генерации априорной
                                # информации мы задаем некоторый коэффициент, который,
                                # собственно, после имея априорную информацию и будем определять.
ε = 0.2;                        # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 0.36;                # Область по T
N, M = 100, 80;                 # Кол-во разбиений по X, T
h = (b-a)/N;                    # шаг по X
τ = (T-t₀)/M;                   # шаг по T
Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
qₙ =      q.(Xₙ);               # Сеточные значения коэффициента линейного усиления
ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
u₀ = u_init.(Xₙ);               # Начальные условия
nothing #hide

u = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
ϕr = phidetermination(qₙ, urₘ, reverse(Xₙ), N, Tₘ, M);                      # Нужно подать инвертированную сетку
ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N, M);       # Серединный корень
f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);   # Положение переходного слоя
f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);  # Значение функции на переходном слое
nothing #hide
#########################################################################################

q₀ = [ 0 for x in qₙ];
ψ₀ = zeros(N+1);
ψl = zeros(M+1);
ψr = zeros(M+1);
S = 100;
β = 0.1;
qf, Js, Qs = minimize(q₀, u₀, ulₘ, urₘ, ψ₀, ψl, ψr, Xₙ, N, Tₘ, M, ε, f1, f2, S = S, β = β)
frames = [1:20; 21:div(S,50):S];    # Первые двадцать без пропусков
frames = collect(1:div(S, 50):S);   # С пропусками, чтобы всего было 50 кадров
make_minimzation_gif(Js, Qs, qₙ, Xₙ, name = "Minimization.gif", frames_to_write = frames, β = β)

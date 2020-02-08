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

u = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
ϕr = phidetermination(qₙ, urₘ, reverse(Xₙ), N, Tₘ, M);                      # Нужно подать инвертированную сетку
ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N, M);       # Серединный корень
f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);   # Положение переходного слоя
f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);  # Значение функции на переходном слое
nothing #hide
#########################################################################################

q₀ = [ 0 for i in 1:N+1];
ψ₀ = zeros(N+1);
ψl = zeros(M+1);
ψr = zeros(M+1);
S = 500;
β = 0.00001;
qf, Js, Qs = minimize(q₀, u₀, ulₘ, urₘ, ψ₀, ψl, ψr, Xₙ, N, Tₘ, M, ε, f1, f2, S = S, β = β)

frames_to_write = collect(1:div(S, 50):S);
a = Animation()
@showprogress "Composing mp4.." for s in frames_to_write
    pQs = plot(xlabel = "x", ylabel="q(x)", size = (800, 800) )
    pQs = plot!(Xₙ, qₙ, label="q(x)")
    pQs = scatter!(Xₙ, Qs[:,s], title="Искомая qˢ(x) при s = $(s)", label=L"q^s(x)")
    pJs = plot(title="Значение функционала на шаге s = $(s)", size = (800, 800))
    pJs = plot!(1:s, Js[1:s], xlims = (1,S), xlabel = "s", ylabel="J(q)", ylimits = (0, maximum(Js)))
    pJs = annotate!(div(S, 10), maximum(Js)*0.1, "β = $(β)")
    pJs = annotate!(div(S, 10), maximum(Js)*0.05, "S = $(S)")
    p = plot(pQs, pJs, size = (1600, 800) )
    frame(a);
end
g = mp4(a, "Minimization.mp4")


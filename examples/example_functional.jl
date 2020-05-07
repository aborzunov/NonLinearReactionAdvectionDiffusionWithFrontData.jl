#' ## Генерация априорной информации

using NonLinearReactionAdvectionDiffusionWithFrontData;
using NonLinearReactionAdvectionDiffusionWithFrontData: phidetermination, Φ;
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;

u_l(t) = -8 #+ cos(2*π * t);
u_r(t) =  4 #+ (1 + sin(2*π * t));
qf(x) = sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
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
qₙ =      qf.(Xₙ);               # Сеточные значения коэффициента линейного усиления
ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
u₀ = u_init.(Xₙ);               # Начальные условия
nothing #hide

u, XX, TP   = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);                              # Нам вернут решение, сетку, f_1^s
#
ϕl          = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
ϕr          = phidetermination(reverse(qₙ), urₘ, reverse(Xₙ), N, Tₘ, M);             # Нужно подать инвертированную сетку
ϕr          = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
ϕ           = Φ(ϕl, ϕr, N, M);                                                       # Серединный корень
f1_data     = f1(ϕ, u, Xₙ, N, M);                                                    # Положение переходного слоя
f2_data     = f2(f1_data, u, Xₙ, N, M);                                              # Значение функции на переходном слое
nothing #hide
#########################################################################################

q₀ = [ 0.60 * x for x in qₙ];
ψ₀ = zeros(N+1);
ψl = zeros(M+1);
ψr = zeros(M+1);
S = 100;
β = 0.05;
q_final, Js, Qs = minimize(q₀, u₀, ulₘ, urₘ, ψ₀, ψl, ψr, Xₙ, N, Tₘ, M, ε, f1_data, f2_data, S = S, β = β, w = 0.0003);

frames = [1:20; 21:div(S,50):S];    # Первые двадцать без пропусков
frames = collect(1:div(S, 50):S);   # С пропусками, чтобы всего было 50 кадров
make_minimzation_gif(Js, Qs, qₙ, Xₙ, name = "Minimization_uniform.gif", frames_to_write = frames, β = β)

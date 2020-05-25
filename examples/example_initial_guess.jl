# TODO: переделать файл в пример использования готового кода
# TODO: все подробности --- в src/doc/asymptotic/asymptotic_example.jl
using NonLinearReactionAdvectionDiffusionWithFrontData
using Dierckx
using LinearAlgebra;
using Plots; pyplot();

u_l(t) = -8;
u_r(t) =  4;
qf(x) = 4*sin(3 * π * x);       # Коэффициент линейного усиления,
ε = 0.15;                        # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 0.55;                # Область по T
N, M = 50, 80;                  # Кол-во разбиений по X, T
h = (b-a)/N;                    # шаг по X
τ = (T-t₀)/M;                   # шаг по T
Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
qₙ =     qf.(Xₙ);               # Сеточные значения коэффициента лин. усиления
ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
u₀ = u_init.(Xₙ, ε=ε, x_tp = 0.1);               # Начальные условия
nothing #hide
##
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
nothing # hide

# Набросаем эскиз
psol = plot(title="Эскиз решения");
plot!(Xₙ, u[:, 1], label="u(x, 0)");
plot!(Xₙ, u[:, div(end,2)], label="u(x, 0.5)");
plot!(Xₙ, u[:, end], label="u(x, 1)");
plot!(Xₙ, ϕl[:,1], label="");
plot!(Xₙ, ϕr[:,1], label="");
plot!(Xₙ, ϕ[:,1], label="")

# Вычислим ``\frac{d f_1}{dt}``.
df1dt = zero(f1_data);
df1dt[1] = (-3 * f1_data[1] + 2 * f1_data[2] - f1_data[3]) / (2 * τ);
for m in 2:M-2
    df1dt[m] = (f1_data[m+1] - f1_data[m-1]) / (2 * τ);
end
df1dt[M+1] = (f1_data[M+1 - 2] - 2 * f1_data[M+1 - 1] + 3 * f1_data[M+1]) / (2 * τ);
pvel = plot(title="Скорость движения переходного слоя");
plot!(Tₘ, f1_data, label="f_1");
plot!(Tₘ, df1dt, label="v_m")
#
# Последние точки некрасивые, обрежем их
x = f1_data[2:end-1];
y = df1dt[2:end-1];
tm = Tₘ[2:end-1];
#
pvel = plot(title="Скорость движения переходного слоя");
plot!(tm, x, label="f_1");
plot!(tm, y, label="v_m")

# Найдем минимальный и максимальный номер узлов сетки Xₙ, которые пересекает
# переходный слой
# .     Узлы сетки
# x     f1(t) положение переходного слоя
#                     k                                   K
# .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
#                   x x x x x x x x x x x x x x x x x x x x x
k = findfirst(z -> z > minimum(x), Xₙ)
K = findfirst(z -> z > maximum(x), Xₙ) - 1
@info "Выбор точек, которые прошел фронт" minimum(x) Xₙ[k] maximum(x) Xₙ[K]

vspl = Spline1D(x, y);
Vn = zero(Xₙ);
# Если мы запрашиваем значения функции для аргумента, который не входит в
# f1_data, вернем 0
for n in k:K
        Vn[n] = vspl(Xₙ[n])
end
scatter!(psol, Xₙ, Vn, label="v_n")


K = K - 3; # Подгонка, потому что скорость на правом крае резко изменяется
Al = LowerTriangular(h * ones(K-k+1, K-k+1));
Ar = UpperTriangular(h * ones(K-k+1, K-k+1));
Al[:, 1]    .= k*h;
Ar[:, end]  .= (N+1-K)*h;

A = -(Al + Ar) / 2.0;
# B составим из ненулевых v_n, без последней точки
B = Vn[k:K];
B = B .- ( -(ulₘ[1] + ulₘ[1])  / 2.0);
qx = A \ B
l = [qx[1] for i in 1:k-1];     # Дополняем начальное приближение константами
r = [qx[end] for i in K+1:N+1]; # Дополняем Справа
q_guess = [l; qx; r];
plot(Xₙ, q_guess, label="Начальное приближение")
plot!(Xₙ, qₙ, label="Истинное")

savefig("guess.png")

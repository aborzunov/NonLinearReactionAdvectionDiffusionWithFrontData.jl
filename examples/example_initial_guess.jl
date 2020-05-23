using NonLinearReactionAdvectionDiffusionWithFrontData
using Dierckx
using LinearAlgebra;

a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = NonLinearReactionAdvectionDiffusionWithFrontData.dparams();
##
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
nothing #hide

# Генерируем экспериментальные данные функцией
# [`NonLinearReactionAdvectionDiffusionWithFrontData.generate_obs_data`](@ref).
ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
nothing # hide

# Вычислим ``\frac{d f_1}{dt}``. ``f_1`` у нас определен на сетке ``T_M``,
# будем использовать формулу правой разности

# Для однообразности, создадим массив такого же размера
# Для всех точек, кроме последней, найдем правую разность
# Для последней -- левую.
df1dt = similar(f1_data);
for k in 1:M
    df1dt[k] = (f1_data[k+1] - f1_data[k]) / (Tₘ[k+1] - Tₘ[k]);
end
df1dt[end] = (f1_data[end-1] - f1_data[end]) / (Tₘ[end-1] - Tₘ[end]);

vspl = Spline1D(f1_data, df1dt);
# Наивный вариант аппроксимации не очень удобен
# Здесь Dierckx создаст интерполяционный объект, который определит значения
# функции вне `f1_data` константами.
Vn = vspl(Xₙ);

# Если мы запрашиваем значения функции для аргумента, который не входит в
# f1_data, вернем 0
for n in 1:N+1
    if Xₙ[n] < f1_data[1]
        Vn[n] = 0.0
    elseif Xₙ[n] > f1_data[end]
        Vn[n] = 0.0
    else
        Vn[n] = vspl(Xₙ[n])
    end
end
plot(Vn)

# Найдем минимальный и максимальный номер узлов сетки Xₙ, которые пересекает
# переходный слой
k = findfirst(x -> x > minimum(f1_data), Xₙ)
K = findfirst(x -> x > maximum(f1_data), Xₙ)

h = Xₙ[2] - Xₙ[1];
Al = LowerTriangular(h * ones(K-k, K-k));
Ar = UpperTriangular(h * ones(K-k, K-k));
Al[:, 1]    .= k*h;
Ar[:, end]  .= K*h;

A = Al + Ar;
B = filter(!iszero, Vn) .- (ulₘ[1] + urₘ[1]);
qx = A \ B
l = [qx[1] for i in 1:k-1];
r = [qx[end] for i in K+1:N+1];
q_guess = [l; qx; r]
plot(Xₙ, q_guess, label="Начальное приближение")
plot!(Xₙ, qₙ, label="Истинное")

savefig("guess.png")

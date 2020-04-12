# Функционал, его градиент и минимизация

@doc raw"""
    J_q(uˢ::Matrix, ψˢ::Matrix,
        Xₙ::Array, N::Int,
        Tₘ::Vector, M::Int) -> Vector

Вычисляет градиент ``J_q(x) = \int_0^T u^s(x, t) \psi^s(x, t) dt ``.
"""
function J_q(uˢ::Matrix, ψˢ::Matrix,
             Xₙ::Matrix, N::Int,
             Tₘ::Vector, M::Int,
             X_q::Vector = [NaN, NaN])

    size(Xₙ) == (N+1, M+1) ||
    throw(ArgumentError("size(Xₙ) == $(size(Xₙ)), N == $(N), M == $(M) " *
                        "Массив Xₙ должен иметь размерность (N+1, M+1)."))

    isDynamicMesh = ! isapprox(Xₙ[:, 1], Xₙ[:, end])
    # Когда мы понимаем, что у нас динамическая сетка,
    # мы должны убедиться, что нам передали необходимую вспомогательную сетку.
    if isDynamicMesh
        if all(isnan.(X_q))
            throw(ArgumentError("Xₙ передан в виде матрицы, значит подразумевается," *
                                "\nчто используется динамическая сетка. Это требует" *
                                "\nпередачи вспомогательной сетки `X_q`, для ``q``"))
        end
    end
    @assert length(Tₘ) == M+1


    # заполним рабочие массивы специальным образом
    # Если сетка динамическая,
    # Нам нужна интерполяция решений на вспомогательную сетку `X_q`
    if isDynamicMesh
        J_q = zero(X_q);                # Градиент будет определен на вспомогательной сетке
        u = zeros(length(X_q), M+1);
        ψ = zeros(length(X_q), M+1);
        for m in 1:M+1
            X = Xₙ[:, m];
            uspl = Spline1D(X, uˢ[:, m]);
            ψspl = Spline1D(X, ψˢ[:, m]);
            u[:, m] = uspl(X_q);
            ψ[:, m] = ψspl(X_q);
        end
    else                        # Если все вычислялось на статических сетках
        u = uˢ;                 # То будем работать с исходными массивами
        ψ = ψˢ;
        J_q = zeros(N+1);       # Градиент определяется на сетке Xₙ
    end

    #  Здесь уже не важно, на какой сетке определены `u`, `ψ`
    #  Вычислим интеграл по формуле трапеций.
    #  ``\int_0^T u^s(x, t) \psi^s(x, t) dt ``
    for m in 1:M
            τ = Tₘ[m+1] - Tₘ[m]
            J_q  += (( u[:, m] .* ψ[:, m] ) + ( u[:, m+1] .* ψ[:, m+1] )).* τ / 2
    end

    return -J_q
end

@doc raw"""
    J(uˢ::Matrix, Xₙ::Array, N::Int,
      Tₘ::Vector, M::Int,
      f1::Vector, f2::Vector,
      qₙˢ::Vector, α::Real = 0.0) -> Real

Вычисляет функционал `` J(\mathbf{x}) = \int_0^T \left( u(f_1(t), t; q^s) - f_2(t) \right)^2 +
\alpha \int_0^1 q^2(x) dx ``.

Использует [`f2`](@ref) для вычисления ``u(f_1(t), t)`` и после по формуле трапеций.
"""
function J(uˢ::Matrix, Xₙ::Array, N::Int,
           Tₘ::Vector, M::Int,
           f1::Vector, f2::Vector,
           qₙˢ::Vector, α::Real = 0.0)
    J = 0.0
    ξ = zeros(M+1);

    ξ = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, uˢ, Xₙ, N, M);

    for m in 1:M
            τ = Tₘ[m+1] - Tₘ[m];
            J += ( (ξ[m] - f2[m])^2 + (ξ[m+1] - f2[m+1])^2) * τ / 2
    end

    return J
end

@doc raw"""
    minimize(q₀::Vector, u₀::Vector,
             ulₘ::Vector, urₘ::Vector,
             ψ₀::Vector,
             ψl::Vector, ψr::Vector,
             Xₙ, N,
             Tₘ::Vector, M,
             ε,
             f1_data::Vector, f2_data::Vector;
             S::Int = 10,
             β::Real = 0.01,
             w::Real = 0.0001,
             create_mesh::Function = x -> [NaN, NaN])  -> Vector, Vector, Matrix


Вернет `qˢ` на оригинальной сетке `Xₙ[:,1]`.
"""
function minimize(q₀::Vector, u₀::Vector,
                  ulₘ::Vector, urₘ::Vector,
                  ψ₀::Vector,
                  ψl::Vector, ψr::Vector,
                  Xₙ, N,
                  Tₘ::Vector, M,
                  ε,
                  f1_data::Vector, f2_data::Vector;
                  S::Int = 10,
                  β::Real = 0.01,
                  w::Real = 0.0001,
                  create_mesh::Function = x -> [NaN, NaN])

    isDynamicMesh = all(isnan.(create_mesh(Xₙ[end - div(end,2)])))
    @info isDynamicMesh ? "Используем оригинальную сетку Xₙ для q₀" : "Используем отдельную сетку для q₀"

    # Создадим интерполяционный объект для q
    # Он остается постоянным, вдоль всего цикла минимизации
    # Сетка Xₙ — статическая, определенная только для начального
    # приближения, то, что она двигается — дело сугубо внутрннее
    # для прямой и сопряженной задачи.
    qspl = Spline1D(Xₙ, q₀)
    k = 100;                                 # Кол-во интервалов в вспомогательной сетке
    X_q = [ first(Xₙ) + n * (last(Xₙ) - first(Xₙ))/Float64(k) for n in 0:k]
    # вектор сеточных значений q на вспомогательной сетке.
    q_aux = qspl(X_q);

    # Проверим, что нам передали q₀ именно на сетке Xₙ
    if ! isapprox(qspl.(Xₙ), q₀)
        throw(ArgumentError("Сеточные значения q₀ определены на какой-то другой сетке, а не Xₙ"))
    end

    #' ## Подготовка к итерационному процессу
    qˢ = q₀;                    # Приближение вектора `q` на текущей итерации.
    uˢ = zeros(N+1, M+1);       # Решение прямой задачи на текущем шаге итерации.
    ψˢ = zeros(N+1, M+1);       # Решение сопряженной задачи на текущем шаге итерации.
    J_values = zeros(S);        # История значений градиента.
    Q_values = zeros(N+1, S);   # История вектора `q`.

    @showprogress "Iterating minimization loop S=$(S)... " for s in 1:S

        Q_values[:, s] .= qˢ; # Сохраним вектор q на текущей итерации

        # Нахождения решения прямой задачи внутри итерационного числа
        # в аргументах отличается лишь вектором `qˢ`.
        # Решение `u` зависит от него как от параметра.
        uˢ, XXˢ, TPˢ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qˢ, create_mesh = create_mesh);

        # Решаем сопряженную задачу, на вход которой подадим
        # `uˢ`      — Параметр, изменяющийся внутри цикла
        # `f1, f2`  — Параметры, сформированные единожды до цикла

        ψˢ = solve_adjoint(ψ₀, XXˢ, N, Tₘ, M, ε, ψl, ψr, qˢ, uˢ, f1_data, f2_data, w=w);

        J_values[s] = J(uˢ, XXˢ, N, Tₘ, M, f1_data, f2_data, qˢ);

        ∇J = J_q(uˢ, ψˢ[:,end:-1:1], XXˢ, N, Tₘ, M, X_q);

        # Делаем шаг минимизации
        q_aux = q_aux - β * ∇J

        qs_spl = Spline1D(X_q, q_aux)
        qˢ = qs_spl(Xₙ);

    end

    return (qˢ, J_values, Q_values)
end

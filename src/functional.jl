# Функционал, его градиент и минимизация

@doc raw"""
    J_q(uˢ::Matrix, ψˢ::Matrix,
             Xₙ::Vector, N::Int,
             Tₘ::Vector, M::Int) -> Vector

"""
function J_q(uˢ::Matrix, ψˢ::Matrix,
             Xₙ::Vector, N::Int,
             Tₘ::Vector, M::Int)

    @assert length(Xₙ) == N+1
    @assert length(Tₘ) == M+1

    J_q = zeros(N+1);

    for m in 1:M+1
        if m != M+1
            τ = Tₘ[m+1] - Tₘ[m]
            J_q  += (( uˢ[:, m] .* ψˢ[:, m] ) + ( uˢ[:, m+1] .* ψˢ[:, m+1] )).* τ / 2
        else
            τ = -Tₘ[m-1] + Tₘ[m]
            J_q  += (( uˢ[:, m] .* ψˢ[:, m] ) + ( uˢ[:, m-1] .* ψˢ[:, m-1] )).* τ / 2
        end
    end

    return -J_q
end

@doc raw"""
    J(uˢ::Matrix, Xₙ::Vector, N::Int,
           Tₘ::Vector, M::Int,
           f1::Vector, f2::Vector,
           qₙˢ::Vector, α::Real = 0.0) -> Real

"""
function J(uˢ::Matrix, Xₙ::Vector, N::Int,
           Tₘ::Vector, M::Int,
           f1::Vector, f2::Vector,
           qₙˢ::Vector, α::Real = 0.0)
    J = 0.0
    ξ = zeros(M+1);

    for m in 1:M+1
        # ищем номер узла в сетке Xₙ в котором координата `x` больше координаты положения фронта `f1(t)`
        k = 0;
        for n in 1:N+1
            if Xₙ[n] > f1[m]
                k = n;
                break;
            end
        end

        # Выбираем 8 элементов значений аппроксимируемой функции и сетки по X
        y =  uˢ[k-3:k+4, m]
        x = Xₙ[k-3:k+4]

        spl = Spline1D(x, y)

        # Вычисляем ``uˢ(f1(t),t;qˢ)``
        ξ[m] = spl(f1[m])
    end

    for m in 1:M+1
        if m != M+1
            τ = Tₘ[m+1] - Tₘ[m];
            J += ( (ξ[m] - f2[m])^2 + (ξ[m+1] - f2[m+1])^2) * τ / 2
        else
            τ = -Tₘ[m-1] + Tₘ[m];
            J += ( (ξ[m] - f2[m])^2 + (ξ[m-1] - f2[m-1])^2) * τ / 2
        end
    end

    return J
end

function minimize(q₀, u₀, ul, ur,
                  ψ₀, ψl, ψr,
                  Xₙ, N,
                  Tₘ, M,
                  ε,
                  f1, f2;
                  S = 10,
                  β = 0.01)

    #' ## Подготовка к итерационному процессу
    #
    qˢ = q₀;              # Приближение вектора `q` на текущей итерации.
    uˢ = zeros(N+1, M+1);                        # Решение прямой задачи на текущем шаге итерации.
    ψˢ = zeros(N+1, M+1);                        # Решение сопряженной задачи на текущем шаге итерации.
    J_values = zeros(S);                    # История значений градиента.
    Q_values = zeros(N+1, S);               # История вектора `q`.

    @showprogress "Iterating minimization loop S=$(S)... " for s in 1:S

        Q_values[:, s] .= qˢ; # Сохраним вектор q на текущей итерации

        # Нахождения решения прямой задачи внутри итерационного числа
        # в аргументах отличается лишь вектором `qˢ`.
        # Решение `u` зависит от него как от параметра.
        uˢ = solve(u₀, Xₙ, N, Tₘ, M, ε, ul, ur, qˢ);

        # Решаем сопряженную задачу, на вход которой подадим
        # `uˢ`      — Параметр, изменяющийся внутри цикла
        # `f1, f2`  — Параметры, сформированные единожды до цикла

        ψˢ = solve_adjoint(ψ₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, qˢ, uˢ, f1, f2);

        J_values[s] = J(uˢ, Xₙ, N, Tₘ, M, f1, f2, qˢ);

        ∇J = J_q(uˢ, ψˢ[:,end:-1:1], Xₙ, N, Tₘ, M);

        qˢ = qˢ - β * ∇J;

    end

    return (qˢ, J_values, Q_values)
end

# Все, что касается сопряженной задачи

@doc raw"""
function adjointRP(y, m::Int, t, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2)

!!! note
    `u`, `Xₙ` передаются как есть, вместе с граничными точками, внутри функции они локально модифицируются,
    для сохранения индексации
"""
function adjointRP(y, m::Int,
                   Xₙ::Vector, N::Int,
                   Tₘ::Vector, M::Int,
                   ε::Real, qₙ::Vector,
                   u::Matrix, f1::Vector, f2::Vector)

    # ВНИМАНИЕ! Выбрасываем первый и последний элемент из сетки Xₙ, для сохранения одинаковой индексации
    Xₙ = Xₙ[2:end-1]
    # Внимание! Выбрасываем первый и последний элемент на каждом временном слое из матрицы `u` для сохранения индексации
    u = u[2:end-1, :]

    @assert size(u)     == (N-1, M+1)

    @assert length(Xₙ)  == N-1   # Модифицированная Сетка по `x` размера N-1
    @assert length(y)   == N-1   # Проверка длины вектор-функции решения на текущем временном шаге.
    @assert length(qₙ)  == N-1   # Проверка длины вектора "неоднородности"

    @assert length(f1)  == M+1
    @assert length(f2)  == M+1

    RP = similar(y)             # Создаем вектор того же типа и размера
    h = Xₙ[2] - Xₙ[1];          # Нижестоящие формулы приведены для равномерной сетки. Вычислим её шаг.

    RP[1]        = - ε * (y[2]   - 2*y[1] )/(h^2) + u[1, m]   * (y[2] )/(2*h) + y[1]   * qₙ[1] - 2 * delta( Xₙ[1], Xₙ, f1[m]) * ( u[1, m] - f2[m] )
    for n in 2:N-2
        RP[n]    = - ε * (y[n+1] - 2*y[n] + y[n-1])/(h^2) + u[n, m]   * (y[n+1] - y[n-1])/(2*h) + y[n]   * qₙ[n] - 2 * delta( Xₙ[n], Xₙ, f1[m]) * ( u[n, m] - f2[m] )
    end
    RP[N-1]      = - ε * (2*y[N-1] + y[N-2])/(h^2) + u[N-1, m] * ( -y[N-2] )/(2*h) + y[N-1] * qₙ[N-1] - 2 * delta( Xₙ[N-1], Xₙ, f1[m]) * ( u[N-1, m] - f2[m] )

    return RP;
end

@doc raw"""

!!! warning
    Массивы `Xₙ`, `Tₘ`, `u`, `f1`, `f2` Передаются **как есть**!
    Они переварачиваются внутри функции локально.
"""
function solve_adjoint(y₀::Vector, Xₙ::Vector, N::Int,
                       Tₘ::Vector, M::Int,
                       ε::Real, qₙ::Vector,
                       u::Matrix, f1::Vector, f2::Vector;
                       α::Complex = complex(0.5, 0.5))

    if Tₘ[end] - Tₘ[1] < 0
        throw(ArgumentError("Time mesh should be passed in default direction [0, T], not in reversed"))
    end

    Tₘ = reverse(Tₘ);
    f1 = reverse(f1);
    f2 = reverse(f2);
    u = u[:,end:-1:1]

    # `u` – матрица содержащая искомую функцию на каждом временном шаге
    u = zeros(N+1, M+1);
    # Запишем граничные условия в матрицу `u` для нулевого шага по времени.
    u[1, :]   .= 0.0;
    u[N+1, :] .= 0.0;
    # Запишем искомый вектор, здесь он соответствует начальным условиям переданным внутрь функции
    u[2:N, M+1] = y₀;

    y = y₀;
    # Создаем якобиан, через замыкание функции
    j(y, m, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2) = ForwardDiff.jacobian( z -> adjointRP(z, m, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2), y)

    for m in 1:M

        W = (I - α * (Tₘ[m+1] - Tₘ[m]) * j(y, m, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2)) \ adjointRP(y, m, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2)
        y = y .+ (Tₘ[m+1] - Tₘ[m])  * real(W);

        # Запишем найденный вектор.
        # Запишем граничные условия в матрицу `u` для нулевого шага по времени.
        # Т.к. `u` имеет размеры (N+1, M+1), то как и для Tₘ не забудем сместить нумерацию на 1.
        u[2:N, m+1] = y

    end

    return u
end

module NonLinearReactionAdvectionDiffusionWithFrontData

using ForwardDiff
using LinearAlgebra
using Dierckx
using LaTeXStrings
using Plots
using Printf
using Missings

export u_init, f, solve!;
export make_plot, make_gif;
export phidetermination, Φ, f1, f2;
export delta, adjointRP, solve_adjoint;

@doc raw"""
    u_init(x::Real; ε = 0.2) -> Real

Начальные условия в виде $(x^2 - x -2) -6 \tanh( -3 \xi)$, где $\xi = \frac{x - 0.25}{ε}$.
!!! warning
    Граничные условия для этих начальных условий должны быть заданы как `(-8, 4)`.

!!! note
    Вы вольны устанавливать начальные условия на свое усмотрения, эта функция внесена в модуль для удобства повествования.

# Example
```jldoctest
julia> NonLinearReactionAdvectionDiffusionWithFrontData.u_init.(0:0.1:1)
11-element Array{Float64,1}:
 -7.993366656356917
 -7.958156688432881
 -5.970893714323723
  1.600893714323723
  3.628156688432881
  3.7433666563569172
  3.759669571706625
  3.7899835485135176
  3.8399991809276504
  3.9099999592207864
  3.9999999979697227
```
"""
function u_init(x::Real; ε::Real = 0.2)
    ξ = (x - 0.25) / ε;
    return ((x^2 - x - 2) - 6 * tanh(-3 * ξ))
end

@doc raw"""
    f(y::Array{<:Real, 1}, t::Real, Xₙ::Array{<:Real, 1},
      N::Int, ε::Real, u_l::Function, u_r::Function, q::Function)

Задает `qₙ = [ q(x) for x in Xₙ[2:N-1] ]` и вызывает [`f`](@ref).
"""
function f(y::Array{<:Real, 1}, t::Real, Xₙ::Array{<:Real, 1}, N::Int, ε::Real, u_l::Function, u_r::Function, q::Function)
    @assert length(Xₙ) == N+1   # Сетка по `x` размера N+1
    @assert length(y)  == N-1   # Проверка длины вектор-функции решения на текущем временном шаге.

    qₙ = [ q(x) for x in Xₙ[2:N] ]

    f(y, t, Xₙ, N, ε, u_l::Function, u_r::Function, qₙ)
end

@doc raw"""
    f(y::Vector, 1}, t::Real, Xₙ::Vector, 1},
      N::Int, ε::Real, u_l::Function, u_r::Function, qₙ::Vector) -> ::Vector

Функция вычисляет вектор правой части с помощью конечно-разностной аппроксимации пространственных производных.

Если записать систему в следующем виде:
```math
    \left\{
    \begin{aligned}
        &\dfrac{d \mathbf{\textbf{y}}}{d t} = \mathbf{\textbf{f}} \, (\mathbf{\textbf{y}},t), \quad t \in (t_0,T],\\
        &\mathbf{\textbf{y}}(t_0) = \mathbf{\textbf{y}}_{init},
    \end{aligned}
    \right.
```
где
```math
    \begin{aligned}
        &\mathbf{\textbf{y}} = \big(u_1 \; u_2 \;  \ldots \; u_{N - 1} \big)^T, \\
        &\mathbf{\textbf{f}} = \big(f_1 \; f_2 \; \ldots \; f_{N - 1}\big)^T, \\
        &\mathbf{\textbf{y}}_{init} = \big(u_{init} (x_1) \; u_{init} (x_2) \; \ldots \; u_{init} (x_{N - 1}) \big)^T.
    \end{aligned}
```
$u_init(x_n)$ вычисляется с помощью [`u_init(x)`](@ref).


То текущая функция определяет вектор-функцию $\mathbf{\textbf{f}}$ следующим образом:
```math
    \begin{aligned}
        &f_1 =       \varepsilon \dfrac{y_{2}        - 2y_1       + u_{left}(t)}{h^2} + y_1       \dfrac{y_{2}        - u_{left}(t)}{2h} - q(x_1) y_1, \\
        &f_n =       \varepsilon \dfrac{y_{n + 1}    - 2y_n       + y_{n - 1}}{h^2}   + y_n       \dfrac{y_{n + 1}    - y_{n - 1}}{2h}   - q(x_n) u_n, \quad n=\overline{2, N-2}, \\
        &f_{N - 1} = \varepsilon \dfrac{u_{right}(t) - 2y_{N - 1} + y_{N - 2}}{h^2}   + y_{N - 1} \dfrac{u_{right}(t) - y_{N - 2}}{2h}   - q(x_{N-1}) y_{N-1}.
    \end{aligned}
```

# Arguments
- `y::Array{<:Real, 1}`:  Вектор размера `N-1` решения системы в текущий момент времени
- `t::Real`:  Текущий момент времени.
- `Xₙ::Array{<:Real, 1}`: Пространственная сетка по `x`.
- `N::Int`: Число -интервалов- сеткию
- `ε::Real`: Малый параметр при старшей производной.
- `u_l::Function`: Функция левого ГУ.
- `u_r::Function`: Функция правого ГУ.
- `q::Vector`: Вектор размера `N-1` представляющий "неоднородность", см. постановку задачи.

!!! note
    Длина вектора `length(Xₙ)` равняется `N+1`, сетка передается полностью, вместе с граничными точками.
!!! warning
    Функция работает по формулам для *равномерной* сетки!.
"""
function f(y::Array{<:Real, 1}, t::Real, Xₙ::Array{<:Real, 1}, N::Int, ε::Real, u_l::Function, u_r::Function, qₙ::Vector)
    @assert length(Xₙ) == N+1   # Сетка по `x` размера N+1
    @assert length(y)  == N-1   # Проверка длины вектор-функции решения на текущем временном шаге.
    @assert length(qₙ) == N-1   # Проверка длины вектора "неоднородности"
    RP = similar(y)              # Создаем вектор того же типа и размера
    h = Xₙ[2] - Xₙ[1];          # Нижестоящие формулы приведены для равномерной сетки. Вычислим её шаг.

    RP[1]        = ε * (y[2]   - 2*y[1]   + u_l(t))/(h^2) + y[1]   * (y[2]   - u_l(t))/(2*h) - y[1]   * qₙ[1]
    for n in 2:N-2
        RP[n]    = ε * (y[n+1] - 2*y[n]   + y[n-1])/(h^2) + y[n]   * (y[n+1] - y[n-1])/(2*h) - y[n]   * qₙ[n]
    end
    RP[N-1]      = ε * (u_r(t) - 2*y[N-1] + y[N-2])/(h^2) + y[N-1] * (u_r(t) - y[N-2])/(2*h) - y[N-1] * qₙ[N-1]

    return RP
end

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

@doc raw"""
    solve!(y::Array{<:Real, 1}, Xₙ::Array{<:Real, 1}, Tₘ::Array{<:Real, 1}, N::Int,
                M::Int, ε::Real, u_l::Function, u_r::Function, qₙ::Vector; α::Complex = complex(0.5, 0.5))

Alias для вызова `solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, qₙ, f, j); `,
где `f = NonLinearReactionAdvectionDiffusionWithFrontData.f`,
`j(y, t, Xₙ, N, ε, u_l, u_r, q) = ForwardDiff.jacobian( z -> f(z, t, Xₙ, N, ε, u_l, u_r, qₙ), y)`.
"""
function solve!(y::Array{<:Real, 1}, Xₙ::Array{<:Real, 1}, Tₘ::Array{<:Real, 1}, N::Int,
                M::Int, ε::Real, u_l::Function, u_r::Function, qₙ::Vector; α::Complex = complex(0.5, 0.5))

    j(y, t, Xₙ, N, ε, u_l, u_r, q) = ForwardDiff.jacobian( z -> f(z, t, Xₙ, N, ε, u_l, u_r, qₙ), y)
    solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, qₙ, f, j);
end

@doc raw"""
    solve!(y::Array{<:Real, 1}, Xₙ::Array{<:Real, 1}, Tₘ::Array{<:Real, 1}, N::Int,
                M::Int, ε::Real, u_l::Function, u_r::Function, qₙ::Vector,
                RP::Function, jac::Function ; α::Complex = complex(0.5, 0.5)) -> Matrix

Функция, которая находит решение с помощью одностадийной схемы Розенброка с комплексным коэффициентом.

На каждом временном шаге, решение находится как:
```math
    \begin{aligned}
        &\mathbf{\textbf{y}}(t_{m + 1}) = \mathbf{\textbf{y}}(t_m) + (t_{m + 1} - t_m) \, \mathrm{Re} \, \mathbf{\textbf{w}}_1,\\
    \end{aligned}
```
где ``W_1`` находится из
```math
\begin{aligned}
    &\left[\mathbf{\textbf{E}} - \alpha \, (t_{m + 1} - t_m) \, \mathbf{\textbf{f}}_\mathbf{\textbf{y}}\Big(\mathbf{\textbf{y}}(t_m),t_m\Big)\right] \, \mathbf{\textbf{w}}_1 = \\
    &\qquad\qquad\qquad\qquad\quad = \mathbf{\textbf{f}} \, \Big(\mathbf{\textbf{y}}(t_m),\frac{t_{m + 1} + t_m}{2}\Big).
\end{aligned}
```

``\mathbf{f}_\mathbf{y}(\mathbf{y}(t_m), t_m)`` — якобиан функции [`f`](@ref) по вектору ``y`` (в момент времени ``t_m``) в момент времени ``t\_m``.
Эта матрица Якоби имеет следущие ненулевые элементы.
```math
\begin{aligned}
    & \left(f_y\right)_{1,1}  & \equiv & \frac{\partial f_1}{\partial y_1} & = & \varepsilon\dfrac{-2}{h^2} - \dfrac{y_{2} - u_{left}(t)}{2h} + q(x_1), \\

     & \left(f_y\right)_{n,n - 1}  & \equiv & \frac{\partial f_n}{\partial y_{n - 1}} & = & \varepsilon \dfrac{1}{h^2} + \dfrac{y_{n}}{2h}, \quad n=\overline{2, N-1},\\

     & \left(f_y\right)_{n,n}  & \equiv & \frac{\partial f_n}{\partial y_{n}} & = &  -\varepsilon \dfrac{2}{h^2} - \dfrac{y_{n+1} - y_{n-1}}{2h} + q(x_n), \quad n=\overline{2, N-2},\\

     & \left(f_y\right)_{n,n + 1}  & \equiv & \frac{\partial f_n}{\partial y_{n + 1}} & = & \varepsilon \dfrac{1}{h^2} - \dfrac{y_{n}}{2h}, \quad n=\overline{1, N-2},\\

     & \left(f_y\right)_{N - 1,N - 1}  & \equiv & \frac{\partial f_{N - 1}}{\partial y_{N - 1}} & = &  \varepsilon \dfrac{-2}{h^2} - \dfrac{u_{right}(t) - y_{N - 2}}{2h} + q(x_N).
\end{aligned}
```


# Arguments
- `y::Array{<:Real, 1}`: Вектор решения системы в текущий момент времени.
- `Xₙ::Array{<:Real, 1}`: Пространственная сетка по `x`.
- `Tₘ::Array{<:Real, 1}`: Пространственная сетка по `x`.
- `N::Int`: Число -интервалов- сетки.
- `M::Int`: Число -интервалов- сетки.
- `ε::Real`: Малый параметр при старшей производной.
- `u_l::Function`: Функция левого ГУ.
- `u_r::Function`: Функция правого ГУ.
- `qₙ::Vector`: Вектор размера `N-1` представляющий "неоднородность", см. постановку задачи.
- `RP::Function`: Функция вычисления вектора ``f`` правой части.
- `jac::Function`: Якобиан функции вычисления правой части по вектору `y` — ``f_y``.
- `α::Complex`: Коэффициент схемы. При `α = 0` — схема Эйлера, при `α = complex(0.5, 0.5)` — схема Розенброка с комплексным коэффициентом.

# Return
Матрицу `N+1, M+1` с искомой функцией на каждом временном шаге.

!!! note

    Длина вектора `length(Xₙ)` равняется `N+1`, сетка передается полностью, вместе с граничными точками.

    Длина вектора `length(Tₘ)` равняется `M+1`, сетка передается полностью, вместе с граничными точками.
"""
function solve!(y::Vector, Xₙ::Vector, Tₘ::Vector, N::Int,
                M::Int, ε::Real, u_l::Function, u_r::Function, qₙ::Vector,
                RP::Function, jac::Function ; α::Complex = complex(0.5, 0.5))

    # `u` – матрица содержащая искомую функцию на каждом временном шаге
    u = zeros(N+1, M+1);
    # Запишем граничные условия в матрицу `u` для нулевого шага по времени.
    u[1, 1]   = u_l(Tₘ[1]); u[N+1, 1] = u_r(Tₘ[1]);
    # Запишем искомый вектор, здесь он соответствует начальным условиям переданным внутрь функции
    u[2:N, 1] = y;

    for m in 1:M

        W = (I - α * (Tₘ[m+1] - Tₘ[m]) * jac(y, Tₘ[m], Xₙ, N, ε, u_l, u_r, qₙ)) \ RP(y, (Tₘ[m+1] + Tₘ[m])/2, Xₙ, N, ε, u_l, u_r, qₙ)
        y = y .+ (Tₘ[m+1] - Tₘ[m])  * real(W);

        # Запишем найденный вектор.
        # Запишем граничные условия в матрицу `u` для нулевого шага по времени.
        # Т.к. `u` имеет размеры (N+1, M+1), то как и для Tₘ не забудем сместить нумерацию на 1.
        u[1, m+1]   = u_l(Tₘ[m+1]);
        u[N+1, m+1] = u_r(Tₘ[m+1]);
        u[2:N, m+1] = y

    end

    return u;
end

@doc raw"""
    delta(x, Xₙ, x₀ = 0) -> ::Number

Вычисляет конечно разностную аппроксимацию дельта функции ``\delta(x; x₀)`` на сетке `Xₙ` исходя из ``\int\limits_a^b \delta(x; x₀) dx = 1``.
"""
function delta(x::Real, Xₙ::Vector, x₀::Real = 0)
    @assert Xₙ[1] <= x <= Xₙ[end]
    @assert Xₙ[1] <= x₀ <= Xₙ[end]

    N = length(Xₙ) - 1
    out = 0.0;

    for n in 1:N
        if ( x >= Xₙ[n]) && ( x < Xₙ[n+1] )
            if ( x₀ >= Xₙ[n] ) && ( x₀ < Xₙ[n+1])
                out = 1 / ( Xₙ[n+1] - Xₙ[n] )
            end
        end
    end
    return out;
end

@doc raw"""
    phidetermination(q, u, u_border, Xₙ, N::Int)

Решает ОДУ для нахождения вырожденных корней.
"""
function phidetermination(q, u, u_border, Xₙ, N::Int)
    @assert length(Xₙ) == N + 1

    phi = similar(Xₙ);
    phi[1] = u_border;

    for n = 1:N
        phi[n+1] = phi[n] + (Xₙ[n+1] - Xₙ[n]) * ( q[n+1] + q[n] ) / 2;
    end

    return phi
end

@doc raw"""
    f1(ϕ, u, Xₙ, N, M)

Находит значение искомой функции на переходном слое ``f_1(t) = u(x_{t.p.}, t)`` путем поиска точки пересечения ``u(x,t)`` и ``ϕ(x)``.

Точка пересечения находится путем интерполяции функции ``u(x,t) - ϕ(x) = 0``.
"""
function f1(ϕ, u, Xₙ, N, M)
    @assert length(Xₙ) == N+1
    @assert size(u) == (N+1, M+1)


    f1 = zeros(M+1);

    # Необходимо найти абсциссу пересечения Φ и u(x, t), т.е. аргумент х на каждом временном шаге
    for m in 1:M+1
        # Начальные условия при разработки программы заданы таким образом,
        # что в начале массива `V` образуются отрицательные элементы,
        # в конце — положительные
        V = u[:,m] - ϕ
        @assert V[1] < 0
        @assert V[end] > 0

        k = 1;
        # Найдем номер элемента, после которого будут стоять положительный элемент
        for n in 1:N
            if V[n]*V[n+1] < 0
                k = n;
                break;
            end
        end

        # try-catch Заглушка, на случай если разбиений по Х мало,
        # а переходный слой подошел слишком близко к границе
        # тогда k + 3 > N
        try
            # Выберем 8 ближайших элементов, по К.Р.А. от нуля функции
            x = Xₙ[k-3:k+4]
            y = V[k-3:k+4]

            # Теперь построим функцию ``x = v(y)``, она задана сеточными значениями выше
            # С помощью интерполяции найдем ``v(0)``
            # WARN интерполяция производится на неравномерной сетке: см массив `y`

            spl = Spline1D(y, x)

            f1[m] =spl(0);

        catch e
            isa(e, BoundsError) || rethrow(e)
            @warn "Bound Error" m k
            break
        end

    end

    return f1
end

function f2(f1, u, Xₙ, N, M)
    @assert size(u) == (N+1, M+1)
    @assert length(f1) == M+1;

    f2 = similar(f1);

    # находим ординату на каждом временном шаге
    for m in 1:M+1
        k = 1;

        # ищем место в массивах
        for n in 1:N+1
            if Xₙ[n] > f1[m]
                k = n
                break
            end
        end

        # try-catch Заглушка, на случай если разбиений по Х мало,
        # а переходный слой подошел слишком близко к границе
        # тогда k + 3 > N
        try
        x = Xₙ[k-3:k+4]
        y = u[k-3:k+4,m]

        spl = Spline1D(x, y)

        f2[m] = spl(f1[m])

        catch e
            isa(e, BoundsError) || rethrow(e)
            @warn "Bound Error" m k
            break
        end
    end

    return f2
end

@doc raw"""
    Φ(ϕ_l::Vector, ϕ_r::Vector, N::Int) -> ::Vector

Вычисляет значение функции на переходном слое ``(\phi_l - \phi_r)/2`` с помощью левого ``\phi_l`` и правого ``\phi_r`` вырожденного корня.
"""
function Φ(ϕ_l::Vector, ϕ_r::Vector, N::Int)
    Φ = similar(ϕ_l)
    Φ = abs.(ϕ_l - ϕ_r)/2 + ϕ_l

    return Φ;
end

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

    J_q = similar(Xₙ);

    for m in 1:M+1
        τ = Tₘ[m+1] - Tₘ[m]
        J_q  += uˢ[:, m] * ψˢ[:, m] .* τ
    end

    return J_q
end

# {{{
@doc raw"""

!!! warning
    Оно не работает!
"""
function jacobian(y, t, Xₙ, ε, u_l, u_r)
    @assert length(Xₙ) == N+1   # Сетка по `x` размера N+1
    @assert length(y)  == N-1   # Проверка длины вектор-функции решения на текущем временном шаге.
    h = Xₙ[2] - Xₙ[1];          # Нижестоящие формулы приведены для равномерной сетки. Вычислим её шаг.


    # Якобиан в первой и последней строке имеет 2 ненулевых элемента
    # Во всех остальных по 3
    I = zeros((N-1 - 2)*3 + 6);
    J = zeros((N-1 - 2)*3 + 6);
    V = zeros((N-1 - 2)*3 + 6);

    # Номера строк первого и второго элемента
    I[1] = 1;
    I[2] = 1;

    # Номера столбцов первого и второго элемента
    J[1] = 1;
    J[2] = 2;

    out = zeros(N-1,N-1)
    out[1,1] = -2*ε/h^2 - (y[2] - u_l)/(2*h) + q(Xₙ[1+1]);
    out[1,2] =    ε/h^2 - y[1]/(2*h);

    for n in 2:N-2
        out[n,n-1] =    ε/h^2 + y[n]/(2*h);
        out[n,n]   = -2*ε/h^2 - (y[n+1] - y[n-1])/(2*h) + q(Xₙ[n+1]);
        out[n,n+1] =    ε/h^2 - y[n]/(2*h)
    end
    out[N-1,N-2] = ε/h^2 + y[N-1]/(2*h)
    out[N-1,N-1] = -2 * ε/h^2 - (p.ur(t) - y[N-2])/(2*h) + q(Xₙ[N-1+1]);

    return out
end
# }}}

# {{{ Plotting
@doc raw"""
    make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector, analitical=nothing;
            frames_to_write::Int = -1, frame_skip::Int=-1, name="solution.gif")

Рисует gif анимацию решения каждый `frame_skip` кадр, вплоть по `frames_to_write`-ый кадр, сохраняет под именем `name`.

Так же рисует аналитическое решение `analitic(x,t)`, если таково передано.

TODO: Fix doc
"""
function make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector,
                  ϕ_l::Vector = missings(2), ϕ_r::Vector = missings(2),
                  f1::Vector = missings(2), f2::Vector = missings(2),
                  analitical = nothing;
                  frames_to_write::Int = -1, frame_skip::Int=-1,
                  name = "solution.gif", convert2mp4 = false)
    N,M = size(u) .-1
    if frames_to_write < 0
        frames_to_write = M;
    end
    if frame_skip < 0
        frame_skip = div(frames_to_write, 40) + 1
    end

    a = Animation()
    for m in 1:frame_skip:frames_to_write
        make_plot(u, Xₙ, Tₘ, m, ϕ_l, ϕ_r, f1, f2)
        frame(a)
    end

    g = gif(a, name, show_msg=false)
    convert2mp4  && run(`gif2mp4 $(name) $(replace(name, "gif" => "mp4")) \&`)

    return g
end

function make_plot(u::Matrix, Xₙ::Vector, Tₘ::Vector, m::Int,
                   ϕ_l::Vector = missings(2), ϕ_r::Vector = missings(2),
                   f1::Vector = missings(2), f2::Vector = missings(2),
                   analitical = nothing)

    yl = extrema(u[:,:]).*1.05;

    # График, оси, подписи осей и пр.
    pl = plot(size=(800, 600), xlabel="x", ylabel="u(x)", ylims=yl)

    # Начальное условие
    plot!(Xₙ, u[:,1], line=:dash, label="u_inital")

    # Найденное решение
    plot!(Xₙ, u[:,m], label="u(x,t)", color=:blue)

    # Точки сетки на найденной функции и их проекция на ось Х
    scatter!(Xₙ, u[:,m], color=:blue, label="", markersize=3)
    scatter!(Xₙ, [0 for x in Xₙ], color=:black, label="", markersize=2)

    # Надпись слева внизу с текущим временем
    annotate!(0.0, 0.9*first(yl), Plots.text(@sprintf("t = %.2f",Tₘ[m]), 16, :left ))

    check(x::Vector) = !any(ismissing.(x))
    if ( check(ϕ_l) ) && ( check(ϕ_r)) && ( check(f1) ) && ( check(f2) )
        ϕ = Φ(ϕ_l, ϕ_r, length(Xₙ)-1);
        plot!(Xₙ, ϕ_l, label=L"\phi_l", color=:darkgoldenrod)
        plot!(Xₙ, ϕ_r, label=L"\phi_r", color=:darkgoldenrod)
        plot!(Xₙ, ϕ, label=L"\widetilde{\Phi}", color=:gold)

        # Плавающая вслед за пунктиром надпись
        annotate!(0.0, f2[m] + 0.5, Plots.text(@sprintf("f2(t) = %.2f",f2[m]), 16, :left ))
        plot!([0, f1[m]], [f2[m], f2[m]], line=:dash, color=:black, label="")
        # красная точка слева, около подписи f2
        scatter!( [0], [f2[m]], color=:red, label="")

        # Плавающая вслед за пунктиром надпись
        annotate!(f1[m] + 0.01, 0.95 * first(yl), Plots.text(@sprintf("f1(t) = %.2f",f1[m]), 16, :left ))
        plot!([f1[m], f1[m]], [yl[1], f2[m]], line=:dash, color=:black, label="")
        # красная точка, нормировочный коэффициент, чтобы влезло в кадр
        scatter!( [f1[m]], [yl[1]*0.99], color=:red, label="")

        # красная точка на пересечении пунктиров априорной информации
        scatter!( [f1[m]], [f2[m]], color=:red, label="")
    end

    # Аналитическое решение
    if analitical != nothing && typeof( analitical ) <: Function
        plot!(Xₙ, analitical.(Xₙ, Tₘ[m]), label="analitical(x,t)", color=:green, linewidth = 5, alpha=0.3)
    end

    return pl
end

# }}}

end # module

# Файл содержит все, что касается решение прямой задачи

@doc raw"""
    f(y::Array{<:Real, 1}, t::Real, Xₙ::Array{<:Real, 1},
      N::Int, ε::Real, u_l::Function, u_r::Function, q::Function)

Задает `qₙ = [ q(x) for x in Xₙ[2:N-1] ]` и вызывает [`f`](@ref).
"""
function f(y::Vector, t::Real, Xₙ::Vector, N::Int, ε::Real, u_l::Function, u_r::Function, q::Function)
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

TODO: fix docs
# Arguments
- `y::Array{<:Real, 1}`:  Вектор размера `N-1` решения системы в текущий момент времени
- `t::Real`:  Текущий момент времени.
- `Xₙ::Array{<:Real, 1}`: Пространственная сетка по `x`.
- `N::Int`: Число -интервалов- сеткию
- `ε::Real`: Малый параметр при старшей производной.
- `ulₙ::Function`: Функция левого ГУ.
- `urₙ::Function`: Функция правого ГУ.
- `q::Vector`: Вектор размера `N-1` представляющий "неоднородность", см. постановку задачи.

!!! warning
    Все вектора передаются вместе с граничными точками.

    Длина должна равняться `N+1`.

!!! note
    Функция работает по формулам для **равномерной** сетки!.
"""
function directRP(y::Vector, m::Int,
                  Xₙ::Vector, N::Int,
                  ε::Real,
                  ulₙ::Vector, urₙ::Vector,
                  qₙ::Vector)

    if length(Xₙ) != N+1   # Сетка по `x` размера N+1
        throw(ArgumentError("""length(Xₙ) == $(length(Xₙ)), N == $(N)
        Массив Xₙ должен передаваться как есть!
        он локально модифицируется внутри функции."""))
    end
    if length(qₙ) != N+1
        throw(ArgumentError("""length(qₙ) == $(length(qₙ)), N == $(N)
        Массив qₙ должен передаваться как есть!
        он локально модифицируется внутри функции."""))
    end
    if length(y) != N-1
        throw(ArgumentError("""length(y) == $(length(y)), N == $(N)
        Массив `y` должен иметь размерность N-1. См. Описание решения."""))
    end

    # Векторы, которые подлежат локальному изменению
    # Дальше, в вычислениях должны использоваться только их локальные копии
    X = strip_borderPoints(Xₙ, N);
    q = strip_borderPoints(qₙ, N);

    RP = similar(y)             # Создаем вектор того же типа и размера
    h = X[2] - X[1];            # Нижестоящие формулы приведены для равномерной сетки. Вычислим её шаг.

    RP[1]        = ε * (y[2]   - 2*y[1]   + ulₙ[m])/(h^2) + y[1]   * (y[2]   - ulₙ[m])/(2*h) - y[1]   * q[1]
    for n in 2:N-2
        RP[n]    = ε * (y[n+1] - 2*y[n]   + y[n-1])/(h^2) + y[n]   * (y[n+1] - y[n-1])/(2*h) - y[n]   * q[n]
    end
    RP[N-1]      = ε * (urₙ[m] - 2*y[N-1] + y[N-2])/(h^2) + y[N-1] * (urₙ[m] - y[N-2])/(2*h) - y[N-1] * q[N-1]

    return RP
end

@doc raw"""
    ∂directRP_∂y(y::Vector, m::Int,
                  Xₙ::Vector, N::Int,
                  ε::Real,
                  ulₙ::Vector, urₙ::Vector,
                  qₙ::Vector)

Функция якобиана для `adjointRP`.
"""
function ∂directRP_∂y(y::Vector, m::Int,
                  Xₙ::Vector, N::Int,
                  ε::Real,
                  ulₙ::Vector, urₙ::Vector,
                  qₙ::Vector)
    return ForwardDiff.jacobian( z -> directRP(z, m, Xₙ, N, ε, ulₙ, urₙ, qₙ), y)
end

@doc raw"""
    solve(y₀::Vector, Xₙ::Vector, N::Int,
               Tₘ::Vector, M::Int,
               ε::Real, ulₙ::Vector, urₙ::Vector,
               qₙ::Vector,
               RP::Function = directRP,
               jac::Function = ∂directRP_∂y;
               α::Complex = complex(0.5, 0.5))

TODO: fix docs
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
function solve(y₀::Vector, Xₙ::Vector, N::Int,
               Tₘ::Vector, M::Int,
               ε::Real, ulₙ::Vector, urₙ::Vector,
               qₙ::Vector,
               RP::Function = directRP,
               jac::Function = ∂directRP_∂y;
               α::Complex = complex(0.5, 0.5))

    # `u` – матрица содержащая искомую функцию на каждом временном шаге
    u = zeros(N+1, M+1);
    # Запишем граничные условия в матрицу `u` для всех шагов по времени.
    u[1  , :]   = ulₙ;
    u[N+1, :] = urₙ;
    # Запишем искомый вектор, здесь он соответствует начальным условиям переданным внутрь функции
    u[:, 1] = y₀;
    # Создадим временный вектор
    y = y₀[2:N];

    # m=1 соответствует начальным условиям
    for m in 2:M+1
        τ = (m != M+1 ?
             (Tₘ[m+1] - Tₘ[m]) :
            -(Tₘ[m-1] - Tₘ[m]));

        j = jac(y, m, Xₙ, N, ε, ulₙ, urₙ, qₙ);
        rp = RP(y, m, Xₙ, N, ε, ulₙ, urₙ, qₙ);

        W = (I - α * τ * j) \ rp;
        y = y .+ τ * real(W);

        # Запишем найденный вектор.
        # Т.к. `u` имеет размеры (N+1, M+1), то как и для Tₘ не забудем сместить нумерацию на 1.
        u[2:N, m] = y

    end

    return u;
end

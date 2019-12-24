module NonLinearReactionAdvectionDiffusionWithFrontData

using ForwardDiff
using LinearAlgebra
using Plots
using Printf

export u_init, f, solve!, make_gif;

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
function solve!(y::Array{<:Real, 1}, Xₙ::Array{<:Real, 1}, Tₘ::Array{<:Real, 1}, N::Int,
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

@doc raw"""
    make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector, analitical=nothing;
            frames_to_write::Int = -1, frame_skip::Int=-1, name="solution.gif")

Рисует gif анимацию решения каждый `frame_skip` кадр, вплоть по `frames_to_write`-ый кадр, сохраняет под именем `name`.

Так же рисует аналитическое решение `analitic(x,t)`, если таково передано.
"""
function make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector, analitical=nothing; frames_to_write::Int = -1, frame_skip::Int=-1, name="solution.gif")
    N,M = size(u)
    if frames_to_write < 0
        frames_to_write = M;
    end
    if frame_skip < 0
        frame_skip = div(frames_to_write, 80) + 1
    end

    yl = extrema(u[:,1:frames_to_write]).*1.05;

    a = Animation()
    for m in 1:frame_skip:frames_to_write
        # График, оси, подписи осей и пр.
        plot(size=(800, 600), xlabel="x", ylabel="u(x)", ylims=yl)
        # Найденное решение
        plot!(Xₙ, u[:,m], label="u(x,t)", color=:blue)
        # Начальное условие
        plot!(Xₙ, u[:,1], line=:dash, label="u_inital")
        # Точки сетки на найденной функции и их проекция на ось Х
        scatter!(Xₙ, u[:,m], color=:blue, label="", markersize=3)
        scatter!(Xₙ, [0 for i in 1:N+1], color=:black, label="", markersize=2)
        # Надпись слева внизу с текущим временем
        annotate!(0.0, 0.9*first(yl), Plots.text(@sprintf("t = %.2f",Tₘ[m]), 16, :left ))
        # Аналитическое решение
        if analitical != nothing
            plot!(Xₙ, analitical.(Xₙ, Tₘ[m]), label="analitical(x,t)", color=:green, linewidth = 5, alpha=0.3)
        end
        frame(a)
    end
    gif(a, name);
end

end # module

module NonLinearReactionAdvectionDiffusionWithFrontData

using ForwardDiff
using LinearAlgebra
using Plots

@doc raw"""
    u_init(x::Real; ε = 0.2) -> Real

Начальные условия в виде $(x^2 - x -2) -6 \tanh( -3 \xi)$, где $\xi = \frac{x - 0.25}{ε}$.
Граничные условия для этих начальных условий должны быть заданы как `(-8, 4)`.

# Example
```jldoctest
julia> u_init.(0:0.1:1)
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
      N::Int, ε::Real, u_l::Function, u_r::Function, q::Function) -> ::Array{<:Real, 1}

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
\mathbf{\textbf{y}} = \big(u_1 \; u_2 \;  \ldots \; u_{N - 1} \big)^T$, $\mathbf{\textbf{f}}
= \big(f_1 \; f_2 \; \ldots \; f_{N - 1}\big)^T$ и $\mathbf{\textbf{y}}_{init} = \big(u_{init}
(x_1) \; u_{init} (x_2) \; \ldots \; u_{init} (x_{N - 1}) \big)^T.
```


То текущая функция определяет вектор-функцию $\mathbf{\textbf{f}}$ следующим образом:
```math
    f_1 = \varepsilon \dfrac{y_{2} - 2y_1 + u_{left}(t)}{h^2} + y_1 \dfrac{y_{2} - u_{left}(t)}{2h} + q(x_1) y_1,
    f_n = \varepsilon \dfrac{y_{n + 1} - 2y_n + y_{n - 1}}{h^2} + y_n \dfrac{y_{n + 1} - y_{n - 1}}{2h} + q(x_n) u_n, \quad \text{для $n=\overline{2, N-2}$},
    f_{N - 1} = \varepsilon \dfrac{u_{right}(t) - 2y_{N - 1} + y_{N - 2}}{h^2} + y_{N - 1} \dfrac{u_{right}(t) - y_{N - 2}}{2h} + q(x_{N-1}) y_{N-1}.
```

# Arguments
- `y::Array{<:Real, 1}`: Вектор решения системы в текущий момент времени.
- `t::Real`: Текущий момент времени.
- `Xₙ::Array{<:Real, 1}`: Пространственная сетка по `x`.
- `N::Int`: Число -интервалов- сеткию
- `ε::Real`: Малый параметр при старшей производной.
- `u_l::Function`: Функция левого ГУ.
- `u_r::Function`: Функция правого ГУ.
- `q::Function`: Функция "неоднородности", см. постановку задачи.

!!! Длина вектора `length(Xₙ)` равняется `N+1`, сетка передается полностью, вместе с граничными точками.
"""
function f(y::Array{<:Real, 1}, t::Real, Xₙ::Array{<:Real, 1}, N::Int, ε::Real, u_l::Function, u_r::Function, q::Function)
    @assert length(Xₙ) == N+1   # Сетка по `x` размера N+1
    @assert length(y)  == N-1   # Проверка длины вектор-функции решения на текущем временном шаге.
    f = similar(y)              # Создаем вектор того же типа и размера
    h = Xₙ[2] - Xₙ[1];          # Нижестоящие формулы приведены для равномерной сетки. Вычислим её шаг.

    f[1]        = ε * (y[2]   - 2*y[1]   + u_l(t))/(h^2) + y[1]   * (y[2]   - u_l(t))/(2*h) - y[1]   * q( Xₙ[1+1] )
    for n in 2:N-2
        f[n]    = ε * (y[n+1] - 2*y[n]   + y[n-1])/(h^2) + y[n]   * (y[n+1] - y[n-1])/(2*h) - y[n]   * q( Xₙ[n+1] )
    end
    f[N-1]      = ε * (u_r(t) - 2*y[N-1] + y[N-2])/(h^2) + y[N-1] * (u_r(t) - y[N-2])/(2*h) - y[N-1] * q( Xₙ[N-1 + 1] )

end

@doc raw"""
    solve!(y::Array{<:Real, 1}, Xₙ::Array{<:Real, 1}, Tₘ::Array{<:Real, 1}, N::Int,
            M::Int, ε::Real, u_l::Function, u_r::Function, q::Function)

# Arguments
- `y::Array{<:Real, 1}`: Вектор решения системы в текущий момент времени.
- `Xₙ::Array{<:Real, 1}`: Пространственная сетка по `x`.
- `Tₘ::Array{<:Real, 1}`: Пространственная сетка по `x`.
- `N::Int`: Число -интервалов- сетки.
- `M::Int`: Число -интервалов- сетки.
- `ε::Real`: Малый параметр при старшей производной.
- `u_l::Function`: Функция левого ГУ.
- `u_r::Function`: Функция правого ГУ.
- `q::Function`: Функция "неоднородности", см. постановку задачи.

!!! Длина вектора `length(Xₙ)` равняется `N+1`, сетка передается полностью, вместе с граничными точками.
"""
function solve!(y::Array{<:Real, 1}, Xₙ::Array{<:Real, 1}, Tₘ::Array{<:Real, 1}, N::Int,
                M::Int, ε::Real, u_l::Function, u_r::Function, q::Function; α = complex(0.5, 0.5))

    # `u` – матрица содержащая искомую функцию на каждом временном шаге
    u = zeros(M+1, N+1);
    # Запишем граничные условия в матрицу `u` на шаге `m`
    u[1,1] = u_l(t);
    u[N+1,1] = u_r(t);
    # Запишем найденный вектор
    u[2:N, 1] = y

    for m in 1:M
        W = (I - α * (Tₘ[m+1] - Tₘ[m]) * 1) \ f(y, (Tₘ[m+1] + Tₘ[m])/2, Xₙ, N, ε, u_l, u_r, q)
        y = y .+ (Tₘ[m+1] - Tₘ[m])  * real(W);
        # Запишем граничные условия в матрицу `u` на шаге `m`
        u[1, m] = u_l(t);
        u[N+1, m] = u_r(t);
        # Запишем найденный вектор
        u[2:N, m] = y
    end

    return u;
end



# {{{
function jacobian(y, t, ε, h, u_l, u_r)


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

function make_gif(u::Matrix, Xn; frame_skip=1)
    N,M = size(u)
    if M > 60
        frame_skip = div(M, 60)
    end

    a = Animation()
    for m in 1:frame_skip:M
        plot(Xn, u[m,:])
        frame(a)
    end
    gif(a, "solution.gif");
end
# }}}

end # module

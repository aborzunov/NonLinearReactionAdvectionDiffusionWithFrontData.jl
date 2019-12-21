module NonLinearReactionAdvectionDiffusionWithFrontData

using ForwardDiff
using LinearAlgebra
using Plots

struct DirectProblem
    a::Real
    b::Real
    t₀::Real
    T::Real
    ε::Real
    N::Int
    M::Int
    ul::Function
    ur::Function
    q::Function
    Xₙ
    Tₘ
    h::Real
    τ::Real
end

function DirectProblem()
    l(t) = -1* cos(π*t/2) - 1; # Граничные условия на левой  границе
    r(t) = 1 * sin(π*t/2) + 1;   # Граничные условия на правой границе
    a, b = 0, 1
    t, T = 0, 1
    N, M = 20, 30
    Xn, h = create_mesh(a, b, N)
    Tm, tau = create_mesh(t, T, M)
    DirectProblem(a, b, t, T, 0.1, N, M, l, r, x->sin(3π * x), Xn, Tm, h, tau )
end

function Base.show(io::IO, p::DirectProblem)
    println(io, "Space domain               [a,b] = [$(p.a), $(p.b)]")
    println(io, "Time domain                [t₀, T] = [$(p.t₀), $(p.T)]")
    println(io, "Small parameter            ε = $(p.ε)")
    println(io, "Left Border Condition      uₗ(0) = $(p.ul(0))")
    println(io, "Right Border Condition     uᵣ(0) = $(p.ur(0))")
    println(io, "Time step                  τ = $(p.τ)")
    println(io, "Space step                 h = $(p.h)")
    #show(io, "Funtion q at space domain  q(Xₙ) = $(p.q.(p.Xₙ))")
    #show(io, p.Xₙ)
    #show(io, p.Tₘ)
end


@doc raw"""
    u_initial(x, p::DirectProblem) -> Real

Возвращает значение функции при данного `x` в начальный момент времени для задачи `p`.

Начальные условия заданы в виде $ \tanh( h \frac{x - x\_{tp}}{ε}) + g $,  где $h,g$ — нормировочные
константы для согласования с граничными условиями. Они вычисляются автоматически, исходя из переданного
агрумента `dp`.
"""
function u_initial(x; p::DirectProblem=DirectProblem(), tp = 0.5)
    l = p.ul(0); # creating shortcut for border condition at `t=0`
    r = p.ur(0); # creating shortcut for border condition at `t=0`
    a = p.a; b = p.b;
    # Коэффициенты, для согласования начального и граничных условий
    h = (l-r)/(tanh((p.a - tp)/p.ε) - tanh( (b-tp)/p.ε))
    g = r - h * tanh( (p.b - tp)/p.ε)
    return h*tanh((x-tp) / p.ε) + g;
end
function u_initial(p::DirectProblem=DirectProblem(), tp = 0.5)
    l = p.ul(0); # creating shortcut for border condition at `t=0`
    r = p.ur(0); # creating shortcut for border condition at `t=0`
    a = p.a; b = p.b;
    # Коэффициенты, для согласования начального и граничных условий
    h = (l-r)/(tanh((p.a - tp)/p.ε) - tanh( (p.b-tp)/p.ε))
    g = r - h * tanh( (p.b - tp)/p.ε)
    return [h*tanh((x-tp) / p.ε) + g for x in p.Xₙ[2:end-1]];
end

@doc raw"""
    u_init(x::Real; ε = 0.2) -> Real

Начальные условия в виде $(x^2 - x -2) -6 \tanh( -3 \xi)$, где $\xi = \frac{x - 0.25}{ε}$.
Граничные условия для этих начальных условий должны быть заданы как `(-8, 4)`.

#Example
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

function create_mesh(l::Real, r::Real, N::Int, type::Symbol= :uniform)
    if type === :uniform
        h = (r - l)/N
        out = [l + i*h for i in 0:N]
    else
        thow(ArgumentError("There is unknown type $type"))
    end
    return out, h
end

@doc raw"""
    f(y::Array{<:Real, 1}, t::Real, Xₙ::Array{<:Real, 1},
      N::Int, ε::Real, u_l::Function, u_r::Function, q::Function) -> ::Array{<:Real, 1}

Функция вычисляет вектор правой части с помощью конечно-разностной аппроксимации пространственных производных.

Если записать систему в следующем виде:
```math
\begin{equation}
    \label{Ex_2_System_of_ODEs}
    \left\{
    \begin{aligned}
        &\dfrac{d \mathbf{\textbf{y}}}{d t} = \mathbf{\textbf{f}} \, (\mathbf{\textbf{y}},t), \quad t \in (t_0,T],\\
        &\mathbf{\textbf{y}}(t_0) = \mathbf{\textbf{y}}_{init},
    \end{aligned}
    \right.
\end{equation}
где $\mathbf{\textbf{y}} = \big(u_1 \; u_2 \;  \ldots \; u_{N - 1} \big)^T$, $\mathbf{\textbf{f}}
= \big(f_1 \; f_2 \; \ldots \; f_{N - 1}\big)^T$ и $\mathbf{\textbf{y}}_{init} = \big(u_{init}
(x_1) \; u_{init} (x_2) \; \ldots \; u_{init} (x_{N - 1}) \big)^T$.
```
То текущая функция определяет вектор-функцию $\mathbf{\textbf{f}}$ следующим образом:
```math
\begin{align*}
    &f_1 = \varepsilon \dfrac{y_{2} - 2y_1 + u_{left}(t)}{h^2} + y_1 \dfrac{y_{2} - u_{left}(t)}{2h} + q(x_1) y_1,
\end{align*}
\begin{align*}
    &f_n = \varepsilon \dfrac{y_{n + 1} - 2y_n + y_{n - 1}}{h^2} + y_n \dfrac{y_{n + 1} - y_{n - 1}}{2h} + q(x_n) u_n, \quad \text{для $n=\overline{2, N-2}$},
\end{align*}
\begin{align*}
    &f_{N - 1} = \varepsilon \dfrac{u_{right}(t) - 2y_{N - 1} + y_{N - 2}}{h^2} + y_{N - 1} \dfrac{u_{right}(t) - y_{N - 2}}{2h} + q(x_{N-1}) y_{N-1}.
\end{align*}
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
    @assert length(y) == N-1    # Проверка верности длины вектор
    @assert length(Xₙ) == N+1   # Сетка по `x` размера N+1
    f = similar(y)              # Создаем вектор того же типа и размера
    h = Xₙ[5] - Xₙ[4];

    f[1] = ε * (y[2] - 2*y[1] + u_l(t))/(h^2) + y[1] * (y[2] - u_l(t))/(2*h) - y[1] * q( Xₙ[1+1] )
    for n in 2:N-2
        f[n] = ε * (y[n+1] - 2*y[n] + y[n-1])/(h^2) + y[n] * (y[n+1] - y[n-1])/(2*h) - y[n] * q( Xₙ[n+1] )
    end
    f[N-1]= ε * (u_r(t) - 2*y[N-1] + y[N-2])/(h^2) + y[N-1] * (u_r(t) - y[N-2])/(2*h) - y[N-1] * q( Xₙ[N-1 + 1] )

end

@doc raw"""
    solve!(y::Array{<:Real, 1}, Xₙ::Array{<:Real, 1}, Tₘ::Array{<:Real, 1}, N::Int,
            M::Int, ε::Real, u_l::Function, u_r::Function, q::Function)

"""
function solve!(y::Array{<:Real, 1}, Xₙ::Array{<:Real, 1}, Tₘ::Array{<:Real, 1}, N::Int,
                M::Int, ε::Real, u_l::Function, u_r::Function, q::Function)

    # `u` – матрица содержащая искомую функцию на каждом временном шаге
    u = zeros(M+1, N+1);
    # Запишем граничные условия в матрицу `u` на шаге `m`
    u[1,1] = u_l(t);
    u[N+1,1] = u_r(t);
    # Запишем найденный вектор
    u[2:N, 1] = y

    for m in 1:M
        W = (I - 0 * complex(0.5, 0.5) * (Tₘ[m+1] - Tₘ[m]) * 1) \ f(y, (Tₘ[m+1] + Tₘ[m])/2, Xₙ, N, ε, u_l, u_r, q)
        y = y .+ (Tₘ[m+1] - Tₘ[m])  * real(W);
        # Запишем граничные условия в матрицу `u` на шаге `m`
        u[1, m] = u_l(t);
        u[N+1, m] = u_r(t);
        # Запишем найденный вектор
        u[2:N, m] = y
    end

    return u;
end

function make_gif(u::Matrix, Xn; frame_skip=1)
    M,N = size(u)

    a = Animation()
    for m in 1:frame_skip:M
        plot(Xn, u[m,:])
        frame(a)
    end
    gif(a, "solution.gif");
end


# {{{
function ∂u∂x(y,i,t, p)
    out = 0.0;
    l = p.ul;
    r = p.ur;
    h = p.h;
    if i == 1
        out = (y[i+1] - l(t))/(2*h);
    elseif i == p.N-1
        out = (r(t) - y[N-2])/(2*h)
    else
        out = (y[i+1] - y[i-1])/(2*h)
    end
    return out;
end

@doc raw"""
    ∂²u∂x²(u,i,t)

`u` — сеточные значения фукции.
`i` — номер узла сетки, для которого записывается аппроксимация.
`t` — время, от него зависят граничные условия.

Конечно разнастная аппроксимация $\dfrac{y_{2} - 2y_1 + u_{left}(t)}{h^2}%
"""
function ∂²u∂x²(u,i,t,p)
    out = 0.0;
    ul = p.ul;
    ur = p.ur;
    h = p.h;
    if i == 1
        out = (u[i+1] - 2*u[i] + ul(t))/(h^2)
    elseif i == p.N-1
        out = (ur(t)  - 2*u[i] + u[i-1])/(h^2)
    else
        out = (u[i+1] - 2*u[i] + u[i-1])/(h^2)
    end
    return out;
end
# }}}

@deprecate rp(y, t, p::DirectProblem) f(y::Array{<:Real, 1}, t::Real, Xₙ::Array{<:Real, 1}, N::Int, ε::Real, u_l::Function, u_r::Function, q::Function)
function rp(y, t, p::DirectProblem)
    @assert length(y) == p.N-1
    N = p.N;
    out = similar(y);
    h = p.h;

    out[1] = p.ε * ( y[2] - 2*y[1] + p.ul(t) )/( h^2 ) + y[1] * ( y[2] - p.ul(t))/(2h) - p.q(p.Xₙ[2]) * y[1];
    for n in 2:N-2
        out[n] = p.ε * ( y[n+1] - 2*y[n] + y[n] )/h^2 + y[n] * ( y[n+1] - y[n-1])/(2/h) - p.q(p.Xₙ[n+1]) * y[n];
    end
    out[N-1] = p.ε * ( p.ur(t) - 2*y[N-1] + y[N-2] )/h^2 + y[N-1] * ( p.ur(t) - y[N-2])/(2/h) - p.q(p.Xₙ[N]) * y[N-1];

    #for n in 1:N-1
        #out[n] += (p.ε*π^2*(1 - 2t) + 2t + π * (1 - 2t)^2 * cos(π * p.Xₙ[n+1]) - p.q(p.Xₙ[n+1]) * (1 - 2t) ) * sin(π * p.Xₙ[n+1])
    #end
    return out;
end

function jacobian(y, t, p)
    N = length(y) + 1
    e = p.ε;
    h = p.h;

    out = zeros(N-1,N-1)
    out[1,1] = -2e/h^2 - (y[2] - p.ul(t))/(2*h) + p.q(p.Xₙ[1]);
    out[1,2] = e/h^2 - y[1]/(2*h);

    for n in 2:N-2
        out[n,n-1] = e/h^2 + y[n]/(2*h);
        out[n,n] = -2*e/h^2 - (y[n+1] - y[n-1])/(2*h) + p.q(p.Xₙ[n]);
        out[n,n+1] = e/h^2 - y[n]/(2*h)
    end
    out[N-1,N-2] = e/h^2 + y[N-1]/(2*h)
    out[N-1,N-1] = -2 * e/h^2 - (p.ur(t) - y[N-2])/(2*h) + p.q(p.Xₙ[N-1]);

    return out
end

function solve(p)
    N = p.N;
    M = p.M;
    T = p.Tₘ;
    out = zeros(N+1, M+1);
    u = u_initial(p)

    ∂rp∂y(y,t) = ForwardDiff.jacobian( z-> rp(z, t, p), y)
    for m in 2:M
        out[2:end-1,m]=u;
        out[1,m] = p.ul(T[m]);
        out[end,m] = p.ur(T[m]);
        #z = plot(u,label="u")
        #z = plot!(rp(u,(T[m+1] + T[m])/2,p),label="rp")
        #display(z)
        #sleep(0.5)
        W = (I - 0 * complex(0.5, 0.5) * p.τ * jacobian(u, T[m], p)) \ rp(u,(T[m+1] + T[m])/2,p);
        u = real(W) .* p.τ;

    end

    return u, out
    #jacobian(u,p.τ,p)
    #∂rp∂y(u,p.τ)
end

end # module

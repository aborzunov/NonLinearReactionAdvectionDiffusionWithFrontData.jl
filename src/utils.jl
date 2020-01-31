# Всякие вспомогательные функции

@doc raw"""
    strip_borderPoints(a::Vector, N) -> Vector

Функция для внутренного использования.

Входящий массив должен быть размера `N+1`.
Обрезает граничные точки слева и справа.
Возвращает массив размера `N-1`.

# Example
```jldoctest
julia> N = 10; a = collect(1:N+1);

julia> NonLinearReactionAdvectionDiffusionWithFrontData.strip_borderPoints(a, N)
9-element Array{Int64,1}:
  2
  3
  4
  5
  6
  7
  8
  9
 10
```
"""
function strip_borderPoints(a::Vector, N)
    @assert length(a) == N+1
    return a[2:N]
end

@doc raw"""
    u_init(x::Real; ε = 0.2) -> Real

Начальные условия в виде $(x^2 - x -2) -6 \tanh( -3 \xi)$, где $\xi = \frac{x - 0.25}{ε}$.

!!! note
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
    delta(x, Xₙ, x₀ = 0) -> ::Real

Вычисляет конечно разностную аппроксимацию дельта функции ``\delta(x; x₀)`` на сетке `Xₙ` исходя из ``\int\limits_a^b \delta(x; x₀) dx = 1``.
"""
function delta(x::Real, Xₙ::Vector, x₀::Real = 0)
    if ! (Xₙ[1] <= x <= Xₙ[end])
        throw(DomainError("`x` находится вне `Xₙ`"))
    end
    if ! (Xₙ[1] <= x₀ <= Xₙ[end])
        throw(DomainError("`x₀` находится вне `Xₙ`"))
    end

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

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
    @assert length(a) == N+1 "Размерность входного вектора отличается от $(N+1)"
    return a[2:N]
end

function strip_borderPoints(a::Matrix, N)
    @assert size(a,1) == N+1 "Размерность входного вектора отличается от $(N+1)"
    return a[2:N,:]
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
 -7.868156688432881
 -5.900893714323723
  1.650893714323725
  3.6581566884328813
  3.753366656356917
  3.749669571706625
  3.7599835485135173
  3.7899991809276505
  3.839999959220786
  3.909999997969723
  3.999999999898918
```
"""
function u_init(x::Real; ε::Real = 0.2)
    ξ = (x - 0.15) / ε;
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
    if ! (Xₙ[1] <= x₀ < Xₙ[end])
        throw(DomainError("`x₀` находится вне полуоткрытого отрезка ``[Xₙ[1], Xₙ[end])`` [$(Xₙ[1]), $(Xₙ[end]))"))
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

@doc raw"""
    deltaw(n::Int, x₀::Real, Xₙ::Vector, N::Int)

Вычисляет ``\delta(x - x₀)``, где `x = Xₙ[n]`.
`n` — номер узла в сетке `Xₙ`(без граничных точек).
Использует более точную конечно разностную аппроксимацию дельта функции [`δw`](@ref).
"""
function deltaw(n::Int, x₀::Real, Xₙ::Vector, N::Int)

    @assert length(Xₙ) == N-1
    @assert 1 <= n <= N-1 "`n` вне корректного отрезка"

    out = δw(Xₙ[n] - x₀)

    return out;
end

@doc raw"""
    δw(x::Real; ω = 0.001)

Линейная аппроксимация ``δ(x)``, `ω` — эвристический парамет.
Подбирается так, чтобы невязка сопряженной задачи
``2 \delta( x - f_1(t)) ( u(x,t) - f_2(t))`` при подстановке
в неё `u(x,t;q)` истинного `q` обнулялась почти везде, **но не всюду**.
"""
function δw(x::Real; ω = 0.001)
    if abs(x/ω) <= 1
        return ( 1 - abs(x/ω))/ω
    else
        return 0.0
    end
end

# Всякие вспомогательные функции

@doc raw"""
    shishkin_mesh(a::Real, b::Real,
                  x_tp::Real, ε::Real,
                  N::Int = 50,
                  C_i::Real = 1.0, K_i::Real = 1.0,
                  C_b::Real = 1.0, K_b::Real = 1.0)

Возвращает кусочнораномерную сетку со сгущением на границах и на переходном слое.

# Arguments
- `a::Int`      Левая граница.
- `b::Int`      Правая граница.
- `x_tp::Real`  Положение переходного слоя.
- `ε::Real`     Малый параметр при старшей производной.
- `C_i::Real`   Масштабирующий коэффициент для ширины внутреннего сгущения.
- `K_i::int`    Масштабирующий коэффициент количества интервалов внутреннего сгущения.
- `C_b::Real`   Масштабирующий коэффициент для ширин пограничных сгущений.
- `K_b::Int`    Масштабирующий коэффициент количеств интервалов пограничных сгущений.
- `K::Int`      Кол-во интервалов вне всех сгущений.

# Return
Кусчноравномерную сетк, на эскизе чтоками изображены узлы сетки.
```md
.....  .  .  .  .  .  .  .  .  .  .......  .  .  .  .  .  .  .  .  .  .  ......
```
"""
function shishkin_mesh(a::Real, b::Real,
                       x_tp::Real, ε::Real,
                       N::Int = 50,
                       C_i::Real = 1.0, K_i::Real = 1.0,
                       C_b::Real = 1.0, K_b::Real = 1.0)

    if ! ( a < x_tp < b )
        throw(ArgumentError("``a, x_tp, b`` должна образовывать строго упорядоченную тройку: $((a, x_tp, b))"))
    end

    N_i = Int(K_i * N);
    N_b = Int(K_b * N);
    x = zeros(N + N_i + 2*N_b + 1);

    if ! all( x -> x>10, [N_i, N_b, N] )
        throw(ArgumentError("Количество точек в сгущениях или вне их меньше 10. Используйте другую комбинацию `N, K_i, K_b` " *
                            "N, N_b, N_i = $((N, N_b, N_i))"));
    end

        # Вычисляем долю узлов от N, которые будут лежать слева от внутреннего слоя
        # Проверяем, не пересекаются ли "зоны влияния" внутреннего и пограничного слоёв
        if ((x_tp - abs(C_i*ε*log(ε))) - (a + abs(C_b*ε*log(ε)))) >= 0

            if ((b - abs(C_b*ε*log(ε))) - (x_tp + abs(C_i*ε*log(ε)))) >= 0
                k_left = ((x_tp - abs(C_i*ε*log(ε))) - (a + abs(C_b*ε*log(ε))))/(b - a - 2*abs(C_i*ε*log(ε)) - 2*abs(C_b*ε*log(ε)));
                @show k_left
            else
                k_left = 1;
            end
        else
            k_left = 0;
        end
        # Вычисляем число узлов, которые буду лежать слева и справа от внутреннего слоя
        N_left = Int(round(k_left*N));
        N_right = N - N_left;

        # Вычисляем соответствующие шаги

        # Проверяем, не пересекаются ли "зоны влияния" внутреннего и пограничного слоёв
        # В разных случаях действуем по-разному.
        if ((x_tp - abs(C_i*ε*log(ε))) - (a + abs(C_b*ε*log(ε)))) >= 0
            if ((b - abs(C_b*ε*log(ε))) - (x_tp + abs(C_i*ε*log(ε)))) >= 0
                h_int = 2*abs(C_i*ε*log(ε))/(N_i);
                h_left = ((x_tp - abs(C_i*ε*log(ε))) - (a + abs(C_b*ε*log(ε))))/N_left;
                h_right = ((b - abs(C_b*ε*log(ε))) - (x_tp + abs(C_i*ε*log(ε))))/N_right;
                h_bound = abs(C_b*ε*log(ε))/(N_b);
            else
                h_bound = abs(C_b*ε*log(ε))/(N_b);
                h_left = ((x_tp - abs(C_i*ε*log(ε))) - (a + abs(C_b*ε*log(ε))))/N_left;
                h_int = (b - (x_tp - abs(C_i*ε*log(ε))))/(N_i + N_b);
            end
        else
            h_int = ((x_tp + abs(C_i*ε*log(ε))) - a)/(N_b + N_i);
            h_right = ((b - abs(C_b*ε*log(ε))) - (x_tp + abs(C_i*ε*log(ε))))/N_right;
            h_bound = abs(C_b*ε*log(ε))/(N_b);
        end

        x[1] = a;

        # Проверяем, не пересекаются ли "зоны влияния" внутреннего и пограничного слоёв
        # В разных случаях действуем по-разному.

        if ((x_tp - abs(C_i*ε*log(ε))) - (a + abs(C_b*ε*log(ε)))) >= 0
            if ((b - abs(C_b*ε*log(ε))) - (x_tp + abs(C_i*ε*log(ε)))) >= 0
                for n = 1:(N_b)
                    x[n + 1] = x[n] + h_bound;
                end
                for n = (N_b + 1):(N_b + N_left)
                    x[n + 1] = x[n] + h_left;
                end
                for n = (N_b + N_left + 1):(N_b + N_left + N_i)
                    x[n + 1] = x[n] + h_int;
                end
                for n = (N_b + N_left + N_i + 1):(N_b + N + N_i)
                    x[n + 1] = x[n] + h_right;
                end
                for n = (N_b + N + N_i + 1):(N + N_i + 2*N_b)
                    x[n + 1] = x[n] + h_bound;
                end
            else
                for n = 1:(N_b)
                    x[n + 1] = x[n] + h_bound;
                end
                for n = (N_b + 1):(N_b + N_left)
                    x[n + 1] = x[n] + h_left;
                end
                for n = (N_b + N_left + 1):(N + N_i + 2*N_b)
                    x[n + 1] = x[n] + h_int;
                end
            end
        else
            for n = 1:(N_b + N_left + N_i)
                x[n + 1] = x[n] + h_int;
            end
            for n = (N_b + N_left + N_i + 1):(N_b + N + N_i)
                x[n + 1] = x[n] + h_right;
            end
            for n = (N_b + N + N_i + 1):(N + N_i + 2*N_b)
                x[n + 1] = x[n] + h_bound;
            end
        end


    return x
end

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
    u_init(x::Real; ε::Real = 0.2, x_tp = 0.15) -> Real

Начальные условия в виде $(x^2 - x -2) -6 \tanh( -3 \xi)$, где $\xi = \frac{x - x_tp}{ε}$.
`x_tp` — положение внутрннего переходного слоя.

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
function u_init(x::Real; ε::Real = 0.2, x_tp = 0.15)
    ξ = (x - x_tp) / ε;
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

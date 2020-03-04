# Содержание:
# Функция вычисления            правой части сопряженной задачи         |adjointRP|
# Функция вычисления якобиана   правой части сопряженной задачи         |∂adjointRP_∂y|
# Функция нахождения решения сопряженной задачи                         |solve_adjoint|
#

@doc raw"""
    heterogenety(n::Int, m::Int,
                Xₙ::Vector, N::Int,
                Tₘ::Vector, M::Int,
                Uₙₘ::Matrix,
                f1::Vector, f2::Vector)

Неоднородность выражающая невязку текущего решения с искомым.
`` 2 \delta( x - f_1(t) ) ( u^s(x,t) - f_2(t) ) ``.
"""
function heterogenety(n::Int, m::Int,
                      Xₙ::Vector, N::Int,
                      Tₘ::Vector, M::Int,
                      Uₙₘ::Matrix,
                      f1::Vector, f2::Vector)

    @assert length(Xₙ) == N-1

    @assert length(Tₘ) == M+1
    @assert length(f1) == M+1
    @assert length(f2) == M+1

    @assert size(Uₙₘ) == (N-1, M+1)

    return 2 * deltaw( n, f1[m], Xₙ, N) * ( Uₙₘ[n, m] - f2[m] )
end

@doc raw"""
function adjointRP(y, m::Int, t, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2)

!!! note
    `u`, `Xₙ` передаются как есть, вместе с граничными точками, внутри функции они локально модифицируются,
    для сохранения индексации
"""
function adjointRP(y::Vector, m::Int,
                   Xₙ::Vector, N::Int,
                   Tₘ::Vector, M::Int,
                   ε::Real,
                   ulₘ::Vector, urₘ::Vector,
                   qₙ::Vector,
                   Uₙₘ::Matrix, f1::Vector, f2::Vector)


    @assert length(Xₙ) == N+1
    @assert length(qₙ) == N-1
    @assert length(y)  == N-1

    @assert (Tₘ[end] - Tₘ[1]) < 0
    @assert length(Tₘ) == M+1
    @assert length(ulₘ) == M+1
    @assert length(urₘ) == M+1
    @assert length(f1) == M+1
    @assert length(f2) == M+1

    @assert size(Uₙₘ) == (N-1, M+1)

    @assert m < length(ulₘ)

    RP = zero(y)                # Создаем нулевой вектор того же типа и размера
    h = Xₙ[2] - Xₙ[1];          # Нижестоящие формулы приведены для равномерной сетки. Вычислим её шаг.
    X = Xₙ[2:N];

    RP[1]        = - ε * (y[2]   - 2*y[1] + ulₘ[m])/(h^2)   + Uₙₘ[1, m]   * (y[2] - ulₘ[m])/(2*h) + y[1]        * qₙ[1]     - heterogenety(1, m, X, N, Tₘ, M, Uₙₘ, f1, f2)
    for n in 2:N-2
        RP[n]    = - ε * (y[n+1] - 2*y[n] + y[n-1])/(h^2)   + Uₙₘ[n, m]   * (y[n+1] - y[n-1])/(2*h) + y[n]      * qₙ[n]     - heterogenety(n, m, X, N, Tₘ, M, Uₙₘ, f1, f2)
    end
    RP[N-1]      = - ε * (urₘ[m] - 2*y[N-1] + y[N-2])/(h^2) + Uₙₘ[N-1, m] * (urₘ[m] - y[N-2] )/(2*h) + y[N-1]   * qₙ[N-1]   - heterogenety(N-1, m, X, N, Tₘ, M, Uₙₘ, f1, f2)

    return RP;
end

@doc raw""""
    ARP_y(y::Vector, m::Int,
          Xₙ::Vector, N::Int,
          Tₘ::Vector, M::Int,
          ε::Real,
          ulₘ::Vector, urₘ::Vector,
          qₙ::Vector,
          Uₙₘ::Matrix, f1::Vector, f2::Vector) ->
          (::Vector, ::Vector, ::Vector)

Возращает три диагонали якоибана правой части сопряженной задачи.

# Return
Возвращает три вектора `dl`, `d`, `dl` элементов якоибана.
Нижняя и верхняя диагонали длины ``N-2``, главная — ``N-1``.
```
    [ d[1]  du[1]                         ]
    [ dl[1] d[2]  du[2]                   ]
    [ 0     dl[2] d[3] du[3]              ]
    [           ...  ...  ...             ]
    [                ...  ...     du[N-2] ]
    [                     dl[N-2] d[N-1]  ]
```

!!! warning
    Работает по формулам равномерной сетки.
"""
function ARP_y(y::Vector, m::Int,
               Xₙ::Vector, N::Int,
               Tₘ::Vector, M::Int,
               ε::Real,
               ulₘ::Vector, urₘ::Vector,
               qₙ::Vector,
               Uₙₘ::Matrix, f1::Vector, f2::Vector)

    @assert length(Xₙ) == N+1
    @assert length(qₙ) == N-1
    @assert length(y) == N-1
    @assert m < length(ulₘ)

    d = zeros(N-1);
    dl = zeros(N-2);
    du = zeros(N-2);

    # Вычислим шаг равномерной сетки
    h = Xₙ[2] - Xₙ[1];


    d[1]  = 2 * ε / h^2 + qₙ[1]
    du[1] = -ε / h^2    + Uₙₘ[1, m] / (2*h)

    # Сдвиг индексов для `du` — см # Return документации функции
    for n in 2:N-2
        dl[n - 1]   = -ε / h^2      - Uₙₘ[n,m] / (2*h)
        d[n]        = 2 * ε /h^2    + qₙ[n]
        du[n]       = -ε / h^2      + Uₙₘ[n,m] / (2*h)
    end

    dl[N-2] = -ε / h^2      - Uₙₘ[N-1,m] / (2*h)
    d[N-1]  =  2 * ε / h^2 + qₙ[N-1]

    return dl, d, du
end

@doc raw"""
    ∂ARP_∂y(y::Vector, m::Int,
          Xₙ::Vector, N::Int,
          Tₘ::Vector, M::Int,
          ε::Real,
          ulₘ::Vector, urₘ::Vector,
          qₙ::Vector,
          Uₙₘ::Matrix, f1::Vector, f2::Vector) -> Tridiagonal

Обертка фнукции [`ARP_y`](@ref), которая возвращает `Tridiagonal( ARP_y(...))`
трехдиагональную матрицы из векторов, которые возвращает `ARP_y`.
"""
function ∂ARP_∂y(y::Vector, m::Int,
               Xₙ::Vector, N::Int,
               Tₘ::Vector, M::Int,
               ε::Real,
               ulₘ::Vector, urₘ::Vector,
               qₙ::Vector,
               Uₙₘ::Matrix, f1::Vector, f2::Vector)
    return Tridiagonal( (ARP_y(y, m, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, Uₙₘ, f1, f2))... )
end

@doc raw"""
    ∂adjointRP_∂y(y::Vector, m::Int,
                   Xₙ::Vector, N::Int,
                   Tₘ::Vector, M::Int,
                   ε::Real,
                   ulₘ::Vector, urₘ::Vector,
                   qₙ::Vector,
                   Uₙₘ::Matrix, f1::Vector, f2::Vector)

Функция якобиана для `adjointRP`.
Размерности входных векторов такие же, как и у [`adjointRP`](@ref).

!!! warn
    Функция не входит в публичный API, поэтому размерность
    входных векторов отличается от аналогичных в [`solve_adjoint`](@ref).
     - `Xₙ`, `qₙ`, `y` размера `N-1`!
"""
function ∂adjointRP_∂y(y::Vector, m::Int,
                   Xₙ::Vector, N::Int,
                   Tₘ::Vector, M::Int,
                   ε::Real,
                   ulₘ::Vector, urₘ::Vector,
                   qₙ::Vector,
                   Uₙₘ::Matrix, f1::Vector, f2::Vector)
    return ForwardDiff.jacobian( z -> adjointRP(z, m, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, Uₙₘ, f1, f2), y)
end

@doc raw"""

!!! warning
    Массивы `Xₙ`, `Tₘ`, `u`, `f1`, `f2` Передаются **как есть**!
    Они переварачиваются внутри функции локально.
"""
function solve_adjoint(y₀::Vector, Xₙ::Vector, N::Int,
                       Tₘ::Vector, M::Int,
                       ε::Real, ulₘ::Vector, urₘ::Vector,
                       qₙ::Vector,
                       Uₙₘ::Matrix, f1::Vector, f2::Vector,
                       RP::Function = adjointRP,
                       jac::Function = ∂ARP_∂y;
                       α::Complex = complex(0.5, 0.5))
    # Checking lengths
    # {{{
    if length(Xₙ) != N+1
        throw(ArgumentError("""length(Xₙ) == $(length(Xₙ)), N == $(N)
        Массив Xₙ должен иметь размерность N+1."""))
    end
    if length(qₙ) != N+1
        throw(ArgumentError("""length(qₙ) == $(length(qₙ)), N == $(N)
        Массив qₙ должен иметь размерность N+1."""))
    end
    if length(y₀) != N+1
        throw(ArgumentError("""length(y) == $(length(y)), N == $(N)
        Массив `y` должен иметь размерность N+1."""))
    end

    if length(Tₘ) != M+1
        throw(ArgumentError("""length(Tₘ) == $(length(Tₘ)), M == $(M)
        Массив `Tₘ` должен иметь размерность M+1"""))
    end
    if length(ulₘ) != M+1
        throw(ArgumentError("""length(ulₘ) == $(length(ulₘ)), M == $(M)
        Массив `ulₘ` должен иметь размерность M+1"""))
    end
    if length(urₘ) != M+1
        throw(ArgumentError("""length(y) == $(length(y)), M == $(M)
        Массив `ulₘ` должен иметь размерность M+1"""))
    end
    if length(f1) != M+1
        throw(ArgumentError("""length(f1) == $(length(f1)), M == $(M)
        Массив `f1` должен иметь размерность M+1."""))
    end
    if length(f2) != M+1
        throw(ArgumentError("""length(y) == $(length(f2)), M == $(M)
        Массив `f2` должен иметь размерность M+1."""))
    end

    if size(Uₙₘ) != (N+1, M+1)
        throw(ArgumentError("size(Uₙₘ) == $(size(Uₙₘ)), != $((N+1,M+1))
                            Матрица решения сопряженной задачи должа быть размера N+1, M+1."))
    end
    # }}}

    if Tₘ[end] - Tₘ[1] < 0
        throw(ArgumentError("Сетка по времени должна передаваться в стандартном
                            виде, а не развернутом, то же касается остальных массивов."))
    end

    Tₘ  = reverse(Tₘ);          # Инвертируем сетку по времени
    ulₘ = reverse(ulₘ);         # И все векторы сеточных значений по времени
    ulₘ = reverse(urₘ);
    f1  = reverse(f1);
    f2  = reverse(f2);
    Uₙₘ   = Uₙₘ[:,end:-1:1]

    u = zeros(N+1, M+1);        # Решение
    u[1  , :]   = ulₘ;          # Запись левого  ГУ на каждом шаге
    u[N+1, :]   = urₘ;          # Запись правого ГУ на каждом шаге
    u[:, 1] = y₀;               # Запись начального условия

    # Векторы, которые подлежат локальному изменению
    # Дальше, в вычислениях должны использоваться только их локальные копии
    y = strip_borderPoints(y₀,  N);      # Вектор содержащий решение на текущем шаге
    q = strip_borderPoints(qₙ,  N);      # Сетка без граничных точек
    U = strip_borderPoints(Uₙₘ, N)

    for m in 1:M
        τ = (Tₘ[m+1] - Tₘ[m]);

        rp = RP(y, m, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, q, U, f1, f2)
        j = jac(y, m, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, q, U, f1, f2)

        W = (I - α * τ * j) \ rp;
        y = y .+ τ * real(W);

        u[2:N, m+1] = y       # Сохраним вектор решения на следующем шаге

    end
    return u
end

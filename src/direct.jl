# Содержание
# Функция вычисления правой части прямой задачи     |directRP|
# Функция вычисления якобиана правой части          |∂directRP_∂y|
# Функция нахождения решения прямой задачи          |solve|

@doc raw"""
    directRP(y::Vector, m::Int,
             Xₙ::Vector, N::Int,
             ε::Real,
             ulₘ::Vector, urₘ::Vector,
             qₙ::Vector)

Вычисляет вектор правой части с помощью конечно-разностной
аппроксимации пространственных производных.

# Arguments
- `y::Vector`:          Вектор размера `N-1`, решение в текущий момент без ГТ.
- `m::Int`:             Номер шага в сетке по времени.
- `Xₙ::Vector`:         Вектор размера `N-1` сетки по ``X``, без ГТ.
- `N::Int`:             Число **узлов** в полной сетке по ``X``.
- `ε::Real`:            Малый параметр при старшей производной.
- `ulₘ::Vector`:        Вектор сеточных значений левого ГУ.
- `urₘ::Vector`:        Вектор сеточных значений правого ГУ.
- `qₙ::Vector`:         Сеточные значения неоднородности без ГТ.

!!! warn
    Функция не входит в публичный API, поэтому размерность
    входных векторов отличается от аналогичных в [`solve`](@ref).
     - `Xₙ`, `qₙ`, `y` размера `N-1`!

!!! note
    Функция работает по формулам для **равномерной** сетки!.
"""
function directRP(y::Vector, m::Int,
                  Xₙ::Vector, N::Int,
                  ε::Real,
                  ulₘ::Vector, urₘ::Vector,
                  qₙ::Vector)

    @assert length(Xₙ) == N-1
    @assert length(qₙ) == N-1
    @assert length(y)  == N-1
    @assert m < length(ulₘ)


    RP = zero(y)                    # Создаем нулевой вектор того же типа и размера
    h = Xₙ[2] - Xₙ[1];              # Нижестоящие формулы приведены для равномерной сетки. Вычислим её шаг.

    RP[1]        = ε * (y[2]   - 2*y[1]   + ulₘ[m])/(h^2) + y[1]   * (y[2]   - ulₘ[m])/(2*h) - y[1]   * qₙ[1]
    for n in 2:N-2
        RP[n]    = ε * (y[n+1] - 2*y[n]   + y[n-1])/(h^2) + y[n]   * (y[n+1] - y[n-1])/(2*h) - y[n]   * qₙ[n]
    end
    RP[N-1]      = ε * (urₘ[m] - 2*y[N-1] + y[N-2])/(h^2) + y[N-1] * (urₘ[m] - y[N-2])/(2*h) - y[N-1] * qₙ[N-1]

    return RP
end

@doc raw""""
    f_y(y::Int, m::Int,
        Xₙ::Vector, N::Int,
        ε::Real,
        ulₘ::Vector, urₘ::Vector,
        qₙ::Vector) -> (::Vector, ::Vector, ::Vector)

Возращает три диагонали якоибана правой части прямой задачи.

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
function f_y(y::Vector, m::Int,
             Xₙ::Vector, N::Int,
             ε::Real,
             ulₘ::Vector, urₘ::Vector,
             qₙ::Vector)

    @assert length(Xₙ) == N-1
    @assert length(qₙ) == N-1
    @assert length(y) == N-1
    @assert m < length(ulₘ)

    d = zeros(N-1);
    dl = zeros(N-2);
    du = zeros(N-2);

    # Вычислим шаг равномерной сетки
    h = Xₙ[2] - Xₙ[1];


    d[1] = - 2 * ε / h^2 + (y[2] - ulₘ[m])/(2*h) - qₙ[1]
    du[1] = ε / h^2      + y[1] / (2*h)

    # Сдвиг индексов для `du` — см # Return документации функции
    for n in 2:N-2
        dl[n - 1]   = ε / h^2      - y[n] / (2*h)
        d[n]        = -2 * ε /h^2 + (y[n+1] - y[n-1]) / (2h) - qₙ[n]
        du[n]   = ε / h^2      + y[n] / (2*h)
    end

    d[N-1] = - 2 * ε / h^2 + (urₘ[m] - y[N-2])/(2*h) - qₙ[N-1]
    dl[N-2] = ε / h^2      - y[N-1] / (2*h)

    return dl, d, du
end

@doc raw"""
∂f_∂y(y::Vector, m::Int,
               Xₙ::Vector, N::Int,
               ε::Real,
               ulₘ::Vector, urₘ::Vector,
               qₙ::Vector) -> Tridiagonal

Обертка фнукции [`f_y`](@ref), которая возвращает `Tridiagonal( f_y(...)`
трехдиагональную матрицы из векторов, которые возвращает `f_y`.
"""
function ∂f_∂y(y::Vector, m::Int,
               Xₙ::Vector, N::Int,
               ε::Real,
               ulₘ::Vector, urₘ::Vector,
               qₙ::Vector)
    return Tridiagonal( (f_y(y, m, Xₙ, N, ε, ulₘ, urₘ, qₙ))... )
end

@doc raw"""
    ∂directRP_∂y(y::Vector, m::Int,
                  Xₙ::Vector, N::Int,
                  ε::Real,
                  ulₘ::Vector, urₘ::Vector,
                  qₙ::Vector)

Функция якобиана для `adjointRP`.
Размерности входных векторов такие же, как и у [`directRP`](@ref).

!!! warn
    Функция не входит в публичный API, поэтому размерность
    входных векторов отличается от аналогичных в [`solve`](@ref).
     - `Xₙ`, `qₙ`, `y` размера `N-1`!
"""
function ∂directRP_∂y(y::Vector, m::Int,
                  Xₙ::Vector, N::Int,
                  ε::Real,
                  ulₘ::Vector, urₘ::Vector,
                  qₙ::Vector)
    return ForwardDiff.jacobian( z -> directRP(z, m, Xₙ, N, ε, ulₘ, urₘ, qₙ), y)
end

@doc raw"""
    solve(y₀::Vector, Xₙ::Vector, N::Int,
               Tₘ::Vector, M::Int,
               ε::Real, ulₘ::Vector, urₘ::Vector,
               qₙ::Vector,
               RP::Function = directRP,
               jac::Function = ∂directRP_∂y;
               α::Complex = complex(0.5, 0.5)) -> Matrix —

Функция, которая находит решение с помощью одностадийной схемы Розенброка с комплексным коэффициентом.

# Arguments
- `y₀::Vector`:         Сеточные значения начального условия.
- `Xₙ::Vector`:         Пространственная сетка по ``X``.
- `N::Int`:             Число **узлов** сетки по ``X``.
- `Tₘ::Vector`:         Пространственная сетка по ``t``.
- `M::Int`:             Число **узлов** сетки по ``t``.
- `ε::Real`:            Малый параметр при старшей производной.
- `ulₘ::Function`:      Сеточные значения левого ГУ.
- `urₘ::Function`:      Сеточные значения правого ГУ.
- `qₙ::Vector`:         Сеточные значения "неоднородности", см. постановку задачи.
- `RP::Function`:       Функция вычисления вектора правой части.
- `jac::Function`:      Якобиан правой части по вектору `y` — ``f_y``.
- `α::Complex`:         Коэффициент схемы. При `α = 0` — схема Эйлера, при `α = complex(0.5, 0.5)` — схема Розенброка с комплексным коэффициентом.

# Return
Матрицу размера ``(N+1, M+1)``, содержащую значения искомой функции на сетках ``X\_n, T\_m``,.

!!! info
     - `Xₙ`, `qₙ`, `y₀`     векторы размера ``N+1``!
     - `Tₘ`, `ulₘ`, `urₘ`   векторы размера ``M+1``!
"""
function solve(y₀::Vector, Xₙ::Vector, N::Int,
               Tₘ::Vector, M::Int,
               ε::Real, ulₘ::Vector, urₘ::Vector,
               qₙ::Vector,
               RP::Function = directRP,
               jac::Function = ∂f_∂y;
               α::Complex = complex(0.5, 0.5))

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

    u = zeros(N+1, M+1);        # Решение
    u[1  , :]   = ulₘ;          # Запись левого  ГУ на каждом шаге
    u[N+1, :]   = urₘ;          # Запись правого ГУ на каждом шаге
    u[:, 1] = y₀;               # Запись начального условия

    # Векторы, которые подлежат локальному изменению
    # Дальше, в вычислениях должны использоваться только их локальные копии
    y = strip_borderPoints(y₀, N);      # Вектор содержащий решение на текущем шаге
    X = strip_borderPoints(Xₙ, N);      # Сетка без граничных точек
    q = strip_borderPoints(qₙ, N);      # Сетка без граничных точек

    for m in 1:M
        τ = (Tₘ[m+1] - Tₘ[m]);

        j = jac(y, m, X, N, ε, ulₘ, urₘ, q);
        rp = RP(y, m, X, N, ε, ulₘ, urₘ, q);

        W = (I - α * τ * j) \ rp;
        y = y .+ τ * real(W);

        u[2:N, m+1] = y       # Сохраним вектор решения на следующем шаге
    end

    return u;
end

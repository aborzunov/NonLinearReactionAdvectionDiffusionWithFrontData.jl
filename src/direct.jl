# Содержание
# Функция вычисления правой части прямой задачи     :   |directRP|
# Функция вычисления трехдиагонального
#                       якобиана правой части       :   |DRP_y|
# Обертка `DRP_y`                                   :   |∂DRP_∂y|
# Функция вычисления якобиана правой части          :   |∂directRP_∂y|
# Функция нахождения решения прямой задачи          :   |solve|

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
- `Xₙ::Vector`:         Вектор размера `N+1` сетки по ``X``, вместе с ГТ.
- `N::Int`:             Число **интервалов** в полной сетке по ``X``.
- `ε::Real`:            Малый параметр при старшей производной.
- `ulₘ::Vector`:        Вектор сеточных значений левого ГУ.
- `urₘ::Vector`:        Вектор сеточных значений правого ГУ.
- `qₙ::Vector`:         Вектор размера `N-1` сеточные значений  неоднородности без ГТ.

!!! warn
    Функция не входит в публичный API, поэтому размерность
    входных векторов отличается от аналогичных в [`solve`](@ref).
     - `qₙ`, `y` размера `N-1`!
     - `Xₙ` размера `N+1`.

"""
function directRP(y::Vector, m::Int,
                  Xₙ::Vector, N::Int,
                  ε::Real,
                  ulₘ::Vector, urₘ::Vector,
                  qₙ::Vector)

    @assert length(Xₙ) == N+1
    @assert length(qₙ) == N-1
    @assert length(y)  == N-1
    @assert m < length(ulₘ)

    RP = zero(y)                    # Создаем нулевой вектор того же типа и размера

    # Здесь нужно применить сдвиг индексов в `Xₙ` на +1!
    # Xₙ[1] соответствует первому узлу нашей сетки — ``x_0``
    # Xₙ[2] — ``x_1``
    # Xₙ[3] — ``x_2``
    #
    # Для удобства, формулы записаны в следующем виде
    # d[1] = ... Xₙ[1+1+1] ...,
    # первая единица соответствует единице в d[1],
    # вторая — берем следующий элемент,
    # третья — совершаем сдвиг индексов для Xₙ
    RP[1]     = 2ε / (Xₙ[1+1+1] - Xₙ[1-1+1]) * ( (y[2] - y[1]) /  (Xₙ[1+1+1] - Xₙ[1+1]) - (y[1] - ulₘ[m])/(Xₙ[1+1] - Xₙ[1-1+1])) +
    + y[1]   * (y[2]   - ulₘ[m])/(Xₙ[1+1+1] - Xₙ[1-1+1]) - y[1] * qₙ[1]

    # Xₙ[n  ] — соответствует ``x_{n-1}``
    # Xₙ[n+1] — соответствует ``x_{n  }``
    # Xₙ[n+2] — соответствует ``x_{n+1}``
    #
    # Для удобства, формулы записаны в следующем виде
    # d[n] = ... Xₙ[n-1+1],
    # первое слагаемое соответствует d[n],
    # второе — берем предыдущий элемент,
    # третье — совершаем сдвиг индексов для Xₙ
    for n in 2:N-2
        RP[n] = 2ε / (Xₙ[n+1+1] - Xₙ[n-1+1]) * ( (y[n+1] - y[n]) /  (Xₙ[n+1+1] - Xₙ[n+1]) - (y[n] - y[n-1])/(Xₙ[n+1] - Xₙ[n-1+1])) +
        + y[n]   * (y[n+1] - y[n-1])/(Xₙ[n+1+1] - Xₙ[n-1+1]) - y[n]   * qₙ[n]
    end

    RP[N-1]   =  2ε / (Xₙ[N+1] - Xₙ[N-2+1]) * ( (urₘ[m] - y[N-1]) /  (Xₙ[N+1] - Xₙ[N-1+1]) - (y[N-1] - y[N-2])/(Xₙ[N-1+1] - Xₙ[N-2+1])) +
    + y[N-1] * (urₘ[m] - y[N-2])/(Xₙ[N+1] - Xₙ[N-2+1]) - y[N-1] * qₙ[N-1]

    return RP
end

@doc raw""""
    DRP(y::Int, m::Int,
        Xₙ::Vector, N::Int,
        ε::Real,
        ulₘ::Vector, urₘ::Vector,
        qₙ::Vector) -> (::Vector, ::Vector, ::Vector)

Возращает три диагонали якоибана правой части прямой задачи,
чтобы после сформировать трехдигональную матрицу `Tridiagonal`.


!!! warn
    Функция не входит в публичный API, поэтому размерность
    входных векторов отличается от аналогичных в [`solve`](@ref).
     - `qₙ`, `y` размера `N-1`!
     - `Xₙ` размера `N+1`.

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

# Example
```
julia> # Данные необохдимо подготовить, например, как в examples/example_direct.jl

julia> # Подготовим массивы, выбросив граничные точки

julia> # ведь тестируемая функция — для внутреннего использования

julia> qq = qₙ[2:N];

julia> y = y₀[2:N];

julia> dl, d, du = NonLinearReactionAdvectionDiffusionWithFrontData.DRP_y(y, 1, Xₙ, N, ε, ulₘ, urₘ, qq);

julia> Tridiagonal( dl, d, du )
49×49 Tridiagonal{Float64,Array{Float64,1}}:
 -991.038   305.462      ⋅         ⋅         ⋅         ⋅         ⋅     …       ⋅         ⋅         ⋅         ⋅         ⋅        ⋅
[...]
```
"""
function DRP_y(y::Vector, m::Int,
             Xₙ::Vector, N::Int,
             ε::Real,
             ulₘ::Vector, urₘ::Vector,
             qₙ::Vector)

    @assert length(Xₙ) == N+1
    @assert length(qₙ) == N-1
    @assert length(y) == N-1
    @assert m < length(ulₘ)

    dl = zeros(N-2);        # Поддиагональные элементы
    d = zeros(N-1);         #    Диагональные элементы
    du = zeros(N-2);        # Наддиагональные элементы

    # Здесь нужно применить сдвиг индексов на +1 для Xₙ!
    # Xₙ[1] соответствует ПЕРВОМУ узлу нашей сетки — ``x_0``
    # Xₙ[2] — ``x_1``
    # Xₙ[3] — ``x_2``
    #
    # Для удобства, формулы записаны в следующем виде
    # d[1] = ... Xₙ[1+1+1] ...,
    # первая единица соответствует единице в d[1],
    # вторая — берем следующий элемент,
    # третья — совершаем сдвиг индексов для Xₙ
    d[1] = 2ε / ( Xₙ[1+1+1] - Xₙ[1-1+1]) * ( -1 / (Xₙ[1+1+1] - Xₙ[1+1]) - 1 / (Xₙ[1+1] - Xₙ[1-1+1])) +
            + (y[2] - ulₘ[m]) / ( Xₙ[1+1+1] - Xₙ[1-1+1])  - qₙ[1]

    du[1] = 2ε / ( Xₙ[1+1+1] - Xₙ[1-1+1]) * ( 1/ (Xₙ[1+1+1] - Xₙ[1+1]) ) +
            + y[1] / (Xₙ[1+1+1] - Xₙ[1-1+1])

    # Xₙ[n  ] — соответствует ``x_{n-1}``
    # Xₙ[n+1] — соответствует ``x_{n  }``
    # Xₙ[n+2] — соответствует ``x_{n+1}``
    #
    # Для удобства, формулы записаны в следующем виде
    # d[n] = ... Xₙ[n-1+1],
    # первое слагаемое соответствует d[n],
    # второе — берем предыдущий элемент,
    # третье — совершаем сдвиг индексов для Xₙ
    for n in 2:N-2
        # Сдвиг индексов для поддиагональных элементов `dl`
        # см # Return документации функции
        dl[n - 1] = 2ε / ( Xₙ[n+1+1] - Xₙ[n-1+1]) * ( 1 / (Xₙ[n+1] - Xₙ[n-1+1])) -
                    y[n] / ( Xₙ[n+1+1] - Xₙ[n-1+1])

        d[n]      = 2ε / ( Xₙ[n+1+1] - Xₙ[n-1+1]) * ( -1 / (Xₙ[n+1+1] - Xₙ[n+1]) - 1 / (Xₙ[n+1] - Xₙ[n-1+1])) +
                + (y[n+1] - y[n-1]) / ( Xₙ[n+1+1] - Xₙ[n-1+1])  - qₙ[n]

        du[n]     = 2ε / ( Xₙ[n+1+1] - Xₙ[n-1+1]) * ( 1/ (Xₙ[n+1+1] - Xₙ[n+1]) ) +
                + y[n] / (Xₙ[n+1+1] - Xₙ[n-1+1])
    end

    dl[N-2] = 2ε / ( Xₙ[N+1] - Xₙ[N-2+1]) * ( 1 / (Xₙ[N-1+1] - Xₙ[N-2+1])) -
                 y[N-1] / ( Xₙ[N+1] - Xₙ[N-2+1])
    d[N-1] =  2ε / ( Xₙ[N+1] - Xₙ[N-2+1]) * ( -1 / (Xₙ[N+1] - Xₙ[N-1+1]) - 1 / (Xₙ[N-1+1] - Xₙ[N-2+1]))             +
                + (urₘ[m] - y[N-2]) / ( Xₙ[N+1] - Xₙ[N-2+1])  - qₙ[N-1]

    return dl, d, du
end

@doc raw"""
    ∂DRP_∂y(y::Vector, m::Int,
            Xₙ::Vector, N::Int,
            ε::Real,
            ulₘ::Vector, urₘ::Vector,
            qₙ::Vector) -> Tridiagonal

Обертка фнукции [`DRP_y`](@ref), которая возвращает `Tridiagonal( DRP_y(...))`
трехдиагональную матрицу из векторов, которые возвращает `DRP_y(...)`.

!!! warn
    Функция не входит в публичный API, поэтому размерность
    входных векторов отличается от аналогичных в [`solve`](@ref).
     - `qₙ`, `y` размера `N-1`!
     - `Xₙ` размера `N+1`.
"""
function ∂DRP_∂y(y::Vector, m::Int,
               Xₙ::Vector, N::Int,
               ε::Real,
               ulₘ::Vector, urₘ::Vector,
               qₙ::Vector)
    return Tridiagonal( (DRP_y(y, m, Xₙ, N, ε, ulₘ, urₘ, qₙ))... )
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
     - `qₙ`, `y` размера `N-1`!
     - `Xₙ` размера `N+1`.
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
          α::Complex = complex(0.5, 0.5)) -> Matrix

Функция, которая находит решение прямой задачи с помощью
одностадийной схемы Розенброка с комплексным коэффициентом.

# Arguments
- `y₀::Vector`:         Сеточные значения начального условия.
- `Xₙ::Vector`:         Пространственная сетка по ``X``.
- `N::Int`:             Число **интервалов** сетки по ``X``.
- `Tₘ::Vector`:         Пространственная сетка по ``t``.
- `M::Int`:             Число **интервалов** сетки по ``t``.
- `ε::Real`:            Малый параметр при старшей производной.
- `ulₘ::Function`:      Сеточные значения левого ГУ.
- `urₘ::Function`:      Сеточные значения правого ГУ.
- `qₙ::Vector`:         Сеточные значения "неоднородности", см. постановку задачи.
- `RP::Function`:       Функция вычисления вектора правой части — `directRP`.
- `jac::Function`:      Якобиан правой части по вектору `y` — `∂DRP_∂y`.
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
               jac::Function = ∂DRP_∂y;
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
    q = strip_borderPoints(qₙ, N);      # Сетка без граничных точек

    for m in 1:M
        τ = (Tₘ[m+1] - Tₘ[m]);

        j = jac(y, m, Xₙ, N, ε, ulₘ, urₘ, q);
        rp = RP(y, m, Xₙ, N, ε, ulₘ, urₘ, q);

        W = (I - α * τ * j) \ rp;
        y = y .+ τ * real(W);

        u[2:N, m+1] = y       # Сохраним вектор решения на следующем шаге
    end

    return u;
end

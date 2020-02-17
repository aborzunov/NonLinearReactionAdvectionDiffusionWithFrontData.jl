@testset "Якобиан сопряженной задачи в трехдиагональном виде." begin

    #' Сначала, сгенирируем априорные данные, на увеличенном количестве узлов.
    u_l(t) = -8 #+ cos(2*π * t);
    u_r(t) =  4 #+ (1 + sin(2*π * t));
    q(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
                                    # задаче необходимо определить, но при генерации априорной
                                    # информации мы задаем некоторый коэффициент, который,
                                    # собственно, после имея априорную информацию и будем определять.
    ε = 0.2;                        # Малый параметр при старшей производной
    a, b = 0, 1;                    # Область по X
    t₀, T = 0.0, 1.0;               # Область по T
    N, M = 250, 500;                 # Увеличенное Кол-во разбиений по X, T
    h = (b-a)/N;                    # шаг по X
    τ = (T-t₀)/M;                   # шаг по T
    Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
    Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
    qₙ =      q.(Xₙ);               # Сеточные значения коэффициента линейного усиления
    ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
    urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
    y₀ = u_init.(Xₙ);               # Начальные условия
    nothing #hide

    u = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);

    #' ## Генерация априорной информации
    ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
    ϕr = phidetermination(qₙ, urₘ, reverse(Xₙ), N, Tₘ, M);                      # Нужно подать инвертированную сетку
    ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
    ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N, M);       # Серединный корень
    f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);   # Положение переходного слоя
    f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);  # Значение функции на переходном слое

    # Подготовим массивы, выбросив граничные точки
    # ведь тестируемая функция — для внутреннего использования
    # И развернем массивы по времени
    X = Xₙ[2:N];
    qq = qₙ[2:N];
    y = y₀[2:N];
    U = u[2:N, :];
    U = reverse(U, dims=2);
    f1 = reverse(f1);
    f2 = reverse(f2);
    Tₘ = reverse(Tₘ);

    dl, d, du = NonLinearReactionAdvectionDiffusionWithFrontData.ARP_y(y, 1, X, N, Tₘ, M, ε, ulₘ, urₘ, qq, U, f1, f2)

    @test length(dl) == N - 2
    @test length(d)  == N - 1
    @test length(du) == N - 2

    @test count(x -> x == 0, dl) == 0
    @test count(x -> x == 0, d)  == 0
    @test count(x -> x == 0, du) == 0

    tdjac = NonLinearReactionAdvectionDiffusionWithFrontData.∂ARP_∂y(y, 1, X, N, Tₘ, M, ε, ulₘ, urₘ, qq, U, f1, f2)
    @test typeof(tdjac) <: Tridiagonal
    @test tdjac == Tridiagonal( dl, d, du )

    # Сравним с якобианом автодифференцирования
    jac = NonLinearReactionAdvectionDiffusionWithFrontData.∂adjointRP_∂y(y, 1, X, N, Tₘ, M, ε, ulₘ, urₘ, qq, U, f1, f2)
    @test isapprox(tdjac, jac)

end

# Юнит тест проверяет корректность решения прямой задачи.
# Алгоритм описан в /docs/src/direct/adjoint_check.md
# Возвращает решение, аналитическое решение, сетку по X, сетку по T.
# ДОЛЖЕН быть в конце файла
@testset "Модельная невзяка для сопряженной задачи." begin

    using NonLinearReactionAdvectionDiffusionWithFrontData: heterogenety, adjointRP, ∂ARP_∂y;

    #' Сначала, сгенирируем априорные данные, на увеличенном количестве узлов.
    u_l(t) = -8 #+ cos(2*π * t);
    u_r(t) =  4 #+ (1 + sin(2*π * t));
    q(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
                                    # задаче необходимо определить, но при генерации априорной
                                    # информации мы задаем некоторый коэффициент, который,
                                    # собственно, после имея априорную информацию и будем определять.
    ε = 0.2;                        # Малый параметр при старшей производной
    a, b = 0, 1;                    # Область по X
    t₀, T = 0.0, 1.0;               # Область по T
    N, M = 250, 500;                 # Увеличенное Кол-во разбиений по X, T
    h = (b-a)/N;                    # шаг по X
    τ = (T-t₀)/M;                   # шаг по T
    Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
    Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
    qₙ =      q.(Xₙ);               # Сеточные значения коэффициента линейного усиления
    ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
    urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
    y₀ = u_init.(Xₙ);               # Начальные условия
    nothing #hide

    u = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);

    #' ## Генерация априорной информации
    ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
    ϕr = phidetermination(qₙ, urₘ, reverse(Xₙ), N, Tₘ, M);                      # Нужно подать инвертированную сетку
    ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
    ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N, M);       # Серединный корень
    f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);   # Положение переходного слоя
    f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);  # Значение функции на переходном слое

    #' ## Решение сопряженной задачи


    #######################################################################################
    # Теперь переходим непосредственно к проверки корректности решения сопряженной задачи #
    #######################################################################################
    #
    ################################ Модельная функция ####################################
    raw"""
        g(n::Int, X = Xₙ, N::Int,
          m::Int, T = Tₘ, M::Int) -> Real

    Модельная функция ``(1-2t)\sin(\pi x)``.
    `n` — номер узла в сетке по X. `m` — номер шага в сетке по T.
    `X` — сетка по X, `N` — кол-во **интервалов** сетки.
    `T` — сетка по T, `M` — кол-во **интервалов** сетки.
    """
    function g(n::Int, X, N::Int,
               m::Int, T, M::Int)
        t = T[m]
        x = X[n]
        return (1 - 2t)*sin(π*x)
    end

    raw"""
        g_d(n::Int, Xₙ::Vector, N::Int,
            m::Int, Tₘ::Vector, M::Int,
            ε::Real, qₙ::Vector,
            u::Matrix, f1::Vector, f2::Vector) -> Real

    Вычисляет невязку, т.е. рез ультат подстановки ``g`` в постановку сопряженной задачи.
    `n` — номер узла в сетке по X. `m` — номер шага в сетке по T.
    `X` — сетка по X, `N` — кол-во **интервалов** сетки.
    `T` — сетка по T, `M` — кол-во **интервалов** сетки.
    """
    function g_d(n::Int, Xₙ::Vector, N::Int,
                 m::Int, Tₘ::Vector, M::Int,
                 ε::Real, qₙ::Vector,
                 u::Matrix, f1::Vector, f2::Vector)
        x = Xₙ[n];
        t = Tₘ[m];
        out  = 2 * sin(π * x) - ( - ε * π^2 * (1 - 2t) * sin(π * x)) + π * (1 - 2t) * cos(π * x) * u[n,m] + qₙ[n] * (1 - 2t) * sin(π * x) - heterogenety(n, m, Xₙ, N, Tₘ, M, u, f1, f2)
        return out;
    end
    #######################################################################################

    Uₙₘ = u;                                        # Сохраним матрицу решения прямой задачи
    ψl = [0.0 for i in 1:M+1];                      # Левые  ГУ для модельной функции
    ψr = [0.0 for i in 1:M+1];                      # Правые ГУ для модельной функции
    y₀ = [g(n, Xₙ, N, M+1, Tₘ, M) for n in 1:N+1];  # Начальные условия для модельной функции
                                                    # Внимание! сопряженная задача — ретроспективная
                                                    # ``y₀ = g(x, T)``
    # Модельное решение найденное с помощью известного аналитического решения
    # Для его генерации используем инвертированную сетку по времени, чтобы
    # массив `ψ_model` и `ψ` имели одно и тоже направление хода времени.
    ψ_model = [ g(n, Xₙ, N, m, Tₘ, M) for n in 1:N+1, m in M+1:-1:1];

    # Создадим функцию, которая будет вычислять вектор правой части с добавлением невязки
    raw"""
        RP(y, m, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2) -> Vector

    # Return
        Вектор размера `N-1`, сеточные значения правой части уравнения,
        для которого `g(x,t)` будет являться решением.

    !!! warning
        Функция не предназначена для самостоятельного использования,
        она передается в качестве аргумента в `solve`, внутри которой
        для нее сформуются аргументы нужной длины, для которых не
        потребуется смещений индексов.
    """
    function RP(y, m, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, u, f1, f2)
        arp = adjointRP(y, m, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, u, f1, f2)
        d = [ g_d(n, Xₙ, N, m, Tₘ, M, ε, qₙ, u, f1, f2) for n in 1:N-1]
        return arp .- d
    end
    # Якобиан сконструируем тем же способом, что и внутри пакета в `src/adjoint.jl`.
    j(y, m, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, U, f1, f2) = ForwardDiff.jacobian( z -> adjointRP(z, m, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, U, f1, f2), y)

    # С использованием автоматического дифференцирования
    ψ = solve_adjoint(y₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, Uₙₘ, f1, f2, RP, j)
    @test all(isapprox.(ψ_model, ψ, atol = 0.01))

    # С использованием трехдиагонального якобиана
    ψ = solve_adjoint(y₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, Uₙₘ, f1, f2, RP, ∂ARP_∂y)
    @test all(isapprox.(ψ_model, ψ, atol = 0.01))

    return (ψ, ψ_model, Xₙ, Tₘ)
end

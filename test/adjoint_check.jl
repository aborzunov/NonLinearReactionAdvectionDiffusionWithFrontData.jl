
@testset "Testing adjoint problem solving with model function" begin

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
    ulₙ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
    urₙ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
    y₀ = u_init.(Xₙ);               # Начальные условия
    nothing #hide

    u = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₙ, urₙ, qₙ);

    # Вырожденные корни
    # TODO: FIX PHIDETERMINATION
    y = u[:,1]
    ϕl = phidetermination(qₙ, y, ulₙ, Xₙ, N::Int);
    # Разворачиваем сетку по Х, а после — решение
    ϕr = phidetermination(qₙ, y, urₙ, Xₙ[end:-1:1], N::Int);
    ϕr = ϕr[end:-1:1];
    # Полуразность вырожденных корней
    ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N);
    # Положение переходного слоя
    f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);
    # Значение функции на переходном слое
    f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);
    nothing #hide

    #' ## Решение сопряженной задачи


    #######################################################################################
    # Теперь переходим непосредственно к проверки корректности решения сопряженной задачи #
    #######################################################################################
    #
    ################################ Модельная функция ####################################
    """
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
    """
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
        out  = 2 * sin(π * x) - ( - ε * π^2 * (1 - 2t) * sin(π * x) ) + π * (1 - 2t) * cos(π * x) * u[n,m] + qₙ[n] * (1 - 2t) * sin(π * x) - 2 * delta(x, Xₙ, f1[m]) * (u[n,m] - f2[m])
        return out;
    end
    #######################################################################################

    Uₙₘ = u;                        # Сохраним матрицу решения прямой задачи
    ψl = [0.0 for i in 1:M+1];      # Левые  ГУ для модельной функции
    ψr = [0.0 for i in 1:M+1];      # Правые ГУ для модельной функции
    y₀ = [g(x,M+1) for x in Xₙ];    # Начальные условия для модельной функции
                                    # Внимание! сопряженная задача — ретроспективная
                                    # ``y₀ = g(x, T)``
    # Модельное решение найденное с помощью известного аналитического решения
    # Для его генерации используем инвертированную сетку по времени, чтобы
    # массив `model` и `ψ` имели одно и тоже направление хода времени.
    model = [ g(n, Xₙ, N, m, Tₘ, M) for n in 1:N+1, m in M+1:-1:1];

    # Создадим функцию, которая будет вычислять вектор правой части с добавлением невязки
    """
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
    function RP(y, m, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2)
        arp = adjointRP(y, M, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2)
        d = [ g_d(n, Xₙ, N, m, Tₘ, M, ε, qₙ, u, f1, f2) for n in 1:N-1]
        return arp .- d
    end
    # Якобиан сконструируем тем же способом, что и внутри пакета в `src/adjoint.jl`.
    j(y, m, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, U, f1, f2) = ForwardDiff.jacobian( z -> adjointRP(z, m, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, U, f1, f2), y)

    ψ = solve_adjoint(y₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, Uₙₘ, f1, f2, RP, j)

    @test_broken isapprox(model, u, rtol = 1E-3)

    return (ψ, model, Xₙ, Tₘ)
end

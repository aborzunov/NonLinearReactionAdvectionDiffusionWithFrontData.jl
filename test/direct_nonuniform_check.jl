@testset "Якобиан прямой задачи в техдиагональном виде на кусчно-раномерной сетке." begin

    using NonLinearReactionAdvectionDiffusionWithFrontData: shishkin_mesh;
    u_l(t) = -8;                    # ГУ
    u_r(t) =  4;                    #
    qf(x) = 4*sin(3 * π * x);       # Коэффициент линейного усиления, который в обратной
    # задаче необходимо определить, но при генерации априорной
    # информации мы задаем некоторый коэффициент, который,
    # собственно, после имея априорную информацию и будем определять.
    ε = 0.01;                       # Малый параметр при старшей производной
    a, b = 0, 1;                    # Область по X
    t₀, T = 0, 0.50;                # Область по T
    x_tp = 0.12;                    # Положение переходного слоя
    M = 500;                        # Кол-во разбиений по X, T
    τ = (T-t₀)/M;                   # шаг по T
    Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
    ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
    urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
    nothing #hide

    #' Новое поведение будет реализоваться с помощью передачи функции создания сетки внутрь `solve` в качестве аргумента.
    #' Эта функция должна принимать только один аргумент — положение переходного слоя и формировать соответствующую сетку.
    #' В пакете есть формирование кусочно-равномерной сетки со сгущениями на границе и на переходном слое.
    #' Создадим замыкание этой функции, которое будет иметь нужную сигнатуру.
    mshfrm(x_tp) = shishkin_mesh(a, b, x_tp, ε, 40, 0.5, 1.0, 1.0, 0.25);
    Xₙ = mshfrm(x_tp); ;                            # Сетка по Х
    N = length(Xₙ) - 1                              # Примем за N длину сетки, что получилась в итоге.
    qₙ =    qf.(Xₙ);                                # Сеточные значения коэффициента линейного усиления
    y₀ = u_init.(Xₙ, ε=ε, x_tp = x_tp);             # Начальные условия

    # Подготовим массивы, выбросив граничные точки
    # ведь тестируемая функция — для внутреннего использования
    qq = qₙ[2:N];
    y = y₀[2:N];

    dl, d, du = NonLinearReactionAdvectionDiffusionWithFrontData.DRP_y(y, 1, Xₙ, N, ε, ulₘ, urₘ, qq)

    @test length(dl) == N - 2
    @test length(d)  == N - 1
    @test length(du) == N - 2

    @test count(x -> x == 0, dl) == 0
    @test count(x -> x == 0, d)  == 0
    @test count(x -> x == 0, du) == 0

    tdjac = NonLinearReactionAdvectionDiffusionWithFrontData.∂DRP_∂y(y, 1, Xₙ, N, ε, ulₘ, urₘ, qq)
    @test typeof(tdjac) <: Tridiagonal
    @test tdjac == Tridiagonal( dl, d, du )

    # Сравним с якобианом автодифференцирования
    jac = NonLinearReactionAdvectionDiffusionWithFrontData.∂directRP_∂y(y, 1, Xₙ, N, ε, ulₘ, urₘ, qq)
    @test isapprox(tdjac, jac)

end

# Решение с динамической сеткой опирается на поиск положения переходного слоя.
# Такого в модельной задаче нет!
#=
# Тест проверяет корректность решения прямой задачи. Алгоритм описан в /docs/src/direct/direct_check.md
# Возвращает решение, аналитическое решение, сетку по X, сетку по T.
# ДОЛЖЕН быть последним в этом файле!
@testset "Модельная невзяка для прямой задачи" begin

    # Зададим параметры для прямой задачи
    u_l(t) = 0;                     # ГУ удовлетворяющие модельной функции
    u_r(t) = 0;                     # ГУ удовлетворяющие модельной функции
    q(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления
    ε = 0.02;                        # Малый параметр при старшей производной
    a, b = 0, 1;                    # Область по X
    t₀, T = 0, 1;                   # Область по T
    M = 800;                  # Кол-во разбиений по X, T
    h = (b-a)/N;                    # шаг по X
    τ = (T-t₀)/M;                   # шаг по T
    Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
    ulₘ= u_l.(Tₘ);                  # Сеточные значения левого  ГУ
    urₘ= u_r.(Tₘ);                  # Сеточные значения правого ГУ

    #' Новое поведение будет реализоваться с помощью передачи функции создания сетки внутрь `solve` в качестве аргумента.
    #' Эта функция должна принимать только один аргумент — положение переходного слоя и формировать соответствующую сетку.
    #' В пакете есть формирование кусочно-равномерной сетки со сгущениями на границе и на переходном слое.
    #' Создадим замыкание этой функции, которое будет иметь нужную сигнатуру.
    mshfrm(x_tp) = shishkin_mesh(a, b, x_tp, ε, 40, 0.5, 1.0, 1.0, 0.25);
    Xₙ = mshfrm(x_tp); ;                            # Сетка по Х
    N = length(Xₙ) - 1                              # Примем за N длину сетки, что получилась в итоге.
    qₙ =    qf.(Xₙ);                                # Сеточные значения коэффициента линейного усиления

    # Зададим модельную функцию и невязку, получаемую после подстановки `g` в исходное уравнение
    function g(x, m)
        t = Tₘ[m]
        (1 - 2t)*sin(π*x)
    end
    function g_d(x::Real, m::Int)
        t = Tₘ[m];
        - ε * π^2 * (1 - 2t) * sin(π * x) + π * (1 - 2t)^2 * sin(π * x) * cos(π * x) - q(x) * (1 -2t) * sin(π * x) + 2 * sin(π * x)
    end

    y₀ = g.(Xₙ, 1);               # Начальные условия
    # Модельное решение найденное с помощью известного аналитического решения
    u_model = [ g(x, m) for x in Xₙ, m in 1:M+1];

    # Создадим функцию, которая будет вычислять вектор правой части с добавлением невязки
    function RP(y, m, Xₙ, N, ε, ulₘ, urₘ, qₙ)
        d = [ g_d(x, m) for x in Xₙ[2:N] ]
        NonLinearReactionAdvectionDiffusionWithFrontData.directRP(y, m, Xₙ, N, ε, ulₘ, urₘ, qₙ) - d
    end
    # Хоть мы и конструируем якобиан с помощью автоматического дифференцирования, примите во внимание, что
    # Якобиан ``f_y`` при добавлении `g_d` останется без изменений, т.к. `g_d` зависит только от ``x,t``.
    # То, что он не зависит от добавления `g_d` можно убедиться изменением порядка этих двух строк, ну а так же на бумаге.
    j(y, m, Xₙ, N, ε, ulₘ, urₘ, qₙ) = ForwardDiff.jacobian( z -> RP(z, m, Xₙ, N, ε, ulₘ, urₘ, qₙ), y)

    # С использованием автоматического дифференцирования
    u, XX, TP = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, RP, j, create_mesh = mshfrm);
    @test all(isapprox.(u_model, u, atol = 0.01))

    # С использованием трехдиагонального якобиана
    u, XX, TP = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, RP, NonLinearReactionAdvectionDiffusionWithFrontData.∂DRP_∂y);
    @test all(isapprox.(u_model, u, atol = 0.01))

    return (u, u_model, Xₙ, Tₘ)
end
=#

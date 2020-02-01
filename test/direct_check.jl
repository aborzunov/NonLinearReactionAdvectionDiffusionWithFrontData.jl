@testset "Direct check" begin

    # Зададим параметры для прямой задачи
    u_l(t) = 0;                     # ГУ удовлетворяющие модельной функции
    u_r(t) = 0;                     # ГУ удовлетворяющие модельной функции
    q(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления
    ε = 0.2;                        # Малый параметр при старшей производной
    a, b = 0, 1;                    # Область по X
    t₀, T = 0, 1;                # Область по T
    N, M = 50, 80;                  # Кол-во разбиений по X, T
    h = (b-a)/N;                    # шаг по X
    τ = (T-t₀)/M;                   # шаг по T
    Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
    Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
    qₙ =      q.(Xₙ);               # Сеточные значения коэффициента линейного усиления
    ulₙ= u_l.(Tₘ);                  # Сеточные значения левого  ГУ
    urₙ= u_r.(Tₘ);                  # Сеточные значения правого ГУ


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
    function RP(y, m, Xₙ, N, ε, ulₙ, urₙ, qₙ)
        NonLinearReactionAdvectionDiffusionWithFrontData.directRP(y, m, Xₙ, N, ε, ulₙ, urₙ, qₙ) - g_d.(Xₙ, m)
    end
    # Хоть мы и конструируем якобиан с помощью автоматического дифференцирования, примите во внимание, что
    # Якобиан ``f_y`` при добавлении `g_d` останется без изменений, т.к. `g_d` зависит только от ``x,t``.
    # То, что он не зависит от добавления `g_d` можно убедиться изменением порядка этих двух строк, ну а так же на бумаге.
    j(y, m, Xₙ, N, ε, ulₙ, urₙ, qₙ) = ForwardDiff.jacobian( z -> RP(z, m, Xₙ, N, ε, ulₙ, urₙ, qₙ), y)

    u = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₙ, urₙ, qₙ, RP, j);

    # md d = [missing, missing];
    # md make_gif(u, Xₙ, Tₘ, d, d, d, d, u_model; frame_skip = div(M,50), frames_to_write=M+1, name="direct_check.gif")

    # md using LaTeXStrings
    # md err = u - u_model
    # md heatmap(Xₙ, Tₘ, err', xlabel=L"X_n", ylabel=L"T_m", title="Absolute Error", size=(1200, 800))

    @test all(isapprox.(u_model, u, atol = 0.01))

end

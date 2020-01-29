
@testset "Testing adjoint problem solving with model function" begin

    # Сначала сгенерируем априорные данные f1, f2 и сеточные значения некоторой u
    # ведь подойдут ЛЮБЫЕ
    # Зададим параметры для прямой задачи
    u_l(t) = -8#*cos(2*π * t);
    u_r(t) = 4# * (1 + sin(2*π * t));
    ε = 0.2;
    a, b = 0, 1; # Область по Х
    t₀, T = 0, 0.20;
    N, M = 500, 800; # Делаем сетку по `x` достаточно густой.
    h = (b-a)/N;
    τ = (T-t₀)/M;
    Xₙ = [a  + n*h for n in 0:N];
    Tₘ = [t₀ + m*τ for m in 0:M];
    q(x) = 4*sin(3 * π * x);
    u = zeros(M+1, N+1);

    # Начальные условия
    y = u_init.( Xₙ[n] for n in 2:N );
    # Некоторая функция q, т.к. мы моделируем априорные данные f1, f2,
    # то она нам известна
    qₙ = [ q(x) for x in Xₙ[2:N] ];
    u= solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, qₙ)

    # Вырожденные корни
    qn = [q(x) for x in Xₙ]
    q_low = NonLinearReactionAdvectionDiffusionWithFrontData.phidetermination(qn, y, u_l(0), Xₙ, N::Int);
    q_top = NonLinearReactionAdvectionDiffusionWithFrontData.phidetermination(qn, y, u_r(0), Xₙ[end:-1:1], N::Int);
    q_top = q_top[end:-1:1];
    # Полуразность вырожденных корней
    ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(q_low, q_top, N);
    # Положение переходного слоя
    f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);
    # Значение функции на переходном слое
    f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);

    #
    #######################################################################################
    # Теперь переходим непосредственно к проверки корректности решения сопряженной задачи #
    #######################################################################################
    #


    # Зададим модельную функцию и невязку, получаемую после подстановки `g` в исходное уравнение
    g(x, t) = (1 - 2t)*sin(π*x)
    function g_d(n::Int, Xₙ::Vector, N::Int,
                 m::Int, Tₘ::Vector, M::Int,
                 ε::Real, qₙ::Vector,
                 u::Matrix, f1::Vector, f2::Vector)
        x = Xₙ[n];
        t = Tₘ[m];
        out = 0.0
        out  = 2 * sin(π * x) - ( - ε * π^2 * (1 - 2t) * sin(π * x) ) + π * (1 - 2t) * cos(π * x) * u[n,m] - qₙ[n] * (1 - 2t) * sin(π * x) - 2 * delta(x, Xₙ, f1[m]) * (u[n,m] - f2[m])
        return out;
    end
    y₀ = [ g(x, 0) for x in Xₙ[2:N]];
    # Модельное решение найденное с помощью известного аналитического решения
    model = [ g(x,t) for x in Xₙ, t in Tₘ];

    # Создадим функцию, которая будет вычислять вектор правой части с добавлением невязки
    function RP(y, m, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2)
        arp = adjointRP(y, M, Xₙ, N, Tₘ, M, ε, qₙ, u, f1, f2)
        d = [ g_d(n, Xₙ, N, m, Tₘ, M, ε, qₙ, u, f1, f2) for n in 1:N-1]
        return arp .- d
    end

    make_gif(ψ[:,end:-1:1], Xₙ, Tₘ[end:-1:1]; frame_skip = div(M,500), frames_to_write=81, name="adjoint.gif", convert2mp4 = true)

    @test_broken isapprox(model, u, rtol = 1E-3)
end

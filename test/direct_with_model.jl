@testset "Testing direct problem solving with model function" begin

    u_l(t) = 0
    u_r(t) = 0
    ε = 0.1;
    a, b = 0, 1;
    t₀, T = 0, 1;
    N, M = 40, 80;
    h = (b-a)/N;
    τ = (T-t₀)/M;
    Xₙ = [a  + n*h for n in 0:N];
    Tₘ = [t₀ + m*τ for m in 0:M];
    q(x) = sin(3 * π * x);
    u = zeros(M+1, N+1);
    model = similar(u);
    model = [ g(x,t) for x in Xₙ, t in Tₘ]

    # Зададим модельную функцию и невязку, получаемую после подстановки `g` в исходное уравнение
    g(x, t) = (1 - 2t)*sin(π*x)
    g_d(x,t) =  - ε * π^2 * (1 - 2t) * sin(π * x) + π * (1 - 2t)^2 * sin(π * x) * cos(π * x) - q(x) * (1 -2t) * sin(π * x) + 2 * sin(π * x)

    y = g.( (Xₙ[n] for n in 2:N), 0 );
    qₙ = [ q(x) for x in Xₙ[2:N] ];

    # Создадим функцию, которая будет вычислять вектор правой части с добавлением невязки
    RP(y, t, Xₙ, N, ε, u_l, u_r, qₙ) = f(y, t, Xₙ, N, ε, u_l, u_r, qₙ) - g_d.(Xₙ[2:N], t)
    # Хоть мы и конструируем якобиан с помощью автоматического дифференцирования, примите во внимание, что
    # Якобиан ``f_y`` при добавлении `g_d` останется без изменений, т.к. `g_d` зависит только от ``x,t``.
    # То, что он не зависит от добавления `g_d` можно убедиться изменением порядка этих двух строк, ну а так же на бумаге.
    j(y, t, Xₙ, N, ε, u_l, u_r, q) = ForwardDiff.jacobian( z -> RP(z, t, Xₙ, N, ε, u_l, u_r, q), y)

    u= solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, qₙ, RP, j)
    @test isapprox(model, u, rtol = 1E-3)
end

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
    α = complex(0.5, 0.5);

    g(x, t) = (1 - 2t)*sin(π*x)
    g_d(x,t) =  - ε * π^2 * (1 - 2t) * sin(π * x) + π * (1 - 2t)^2 * sin(π * x) * cos(π * x) - q(x) * (1 -2t) * sin(π * x) + 2 * sin(π * x)

    RP = NonLinearReactionAdvectionDiffusionWithFrontData.f
    j(y, t, Xₙ, N, ε, u_l, u_r, q) = ForwardDiff.jacobian( z -> RP(z, t, Xₙ, N, ε, u_l, u_r, q), y)

    y = g.( (Xₙ[n] for n in 2:N), 0 );
    # `u` – матрица содержащая искомую функцию на каждом временном шаге
    u = zeros(N+1, M+1);
    # Запишем граничные условия в матрицу `u` для нулевого шага по времени.
    u[1, 1]   = u_l(Tₘ[1]);
    u[N+1, 1] = u_r(Tₘ[1]);
    # Запишем искомый вектор, здесь он соответствует начальным условиям переданным внутрь функции
    u[2:N, 1] = y;
    for m in 1:M
        W = (I - α * (Tₘ[m+1] - Tₘ[m]) * j(y, Tₘ[m], Xₙ, N, ε, u_l, u_r, q)) \ ( f(y, (Tₘ[m+1] + Tₘ[m])/2, Xₙ, N, ε, u_l, u_r, q) - g_d.(Xₙ[2:N], (Tₘ[m+1] + Tₘ[m])/2))
        y = y .+ (Tₘ[m+1] - Tₘ[m])  * real(W);
        # Запишем найденный вектор.
        # Запишем граничные условия в матрицу `u` для нулевого шага по времени.
        # Т.к. `u` имеет размеры (N+1, M+1), то как и для Tₘ не забудем сместить нумерацию на 1.
        u[1, m+1]   = u_l(Tₘ[m+1]);
        u[N+1, m+1] = u_r(Tₘ[m+1]);
        u[2:N, m+1] = y
    end

    model = similar(u);
    model = [ g(x,t) for x in Xₙ, t in Tₘ]

    @test isapprox(model, u, rtol = 1E-3)
end

@testset "Различные вспомогательные функции                " begin

    N = 20;
    X = [i/N for i in 0:N];
    Y = u_init.(X);

    @testset "Initial condition" begin

        @test length(u_init.(X)) == N+1
        @test_broken isapprox(u_init(X[1]), -8, rtol=0.001)
        @test isapprox(u_init(X[end]), 4, rtol=0.001)
    end

    @testset "`strip_borderPoints`" begin
        B = copy(Y)
        B = NonLinearReactionAdvectionDiffusionWithFrontData.strip_borderPoints(B, N)
        @test length(B) == N-1
        @test B == Y[2:N]
        @test_throws AssertionError NonLinearReactionAdvectionDiffusionWithFrontData.strip_borderPoints(B, N)
    end


end

@testset "Формирование сетки                               " begin

    using NonLinearReactionAdvectionDiffusionWithFrontData: shishkin_mesh;

    e = 0.02;
    a, b = 0, 1;
    xtp = 0.3;
    N = 50;
    X = shishkin_mesh(a, b, xtp, e, N);

    @test sort(X) == X
    @test unique(X) == X
    @test length(X) == N + N + 2N + 1

    # Проверка упорядоченности a, x_tp, b
    @test_throws ArgumentError shishkin_mesh(b, a, 0.5, 0.1, N)
    @test_throws ArgumentError shishkin_mesh(a, b, 1.5, 0.1, N)
    @test_throws ArgumentError shishkin_mesh(a, b, -10, 0.1, N)

end

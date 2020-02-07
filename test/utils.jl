@testset "Utils" begin

    N = 20;
    X = [i/N for i in 0:N];
    Y = u_init.(X);

    @testset "Initial condition" begin

        @test length(u_init.(X)) == N+1
        @test isapprox(u_init(X[1]), -8, rtol=0.001)
        @test isapprox(u_init(X[end]), 4, rtol=0.001)
    end

    @testset "`strip_borderPoints`" begin
        B = copy(Y)
        B = NonLinearReactionAdvectionDiffusionWithFrontData.strip_borderPoints(B, N)
        @test length(B) == N-1
        @test B == Y[2:N]
        @test_throws AssertionError NonLinearReactionAdvectionDiffusionWithFrontData.strip_borderPoints(B, N)
    end

    @testset "`delta` function" begin
        d = [ delta(x, X, 0.5) for x in X ]

        @testset "формула прямоугольников" begin
        I = 0.0;
            for n in 1:N
                δx = X[n+1] - X[n]
                I += d[n] * δx
            end
            @test isapprox(I, 1.0)
        end

        @testset "формула трапеций" begin
        I = 0.0;
            for n in 1:N
                δx = X[n+1] - X[n]
                I += ((d[n] + d[n+1]) / 2.0) * δx
            end
            @test isapprox(I, 1.0)
        end

        # Возвращает только один ненулевой элемент
        @test length( filter( x -> x != 0, d) ) == 1

        @testset "Проверка области определения" begin
            @test_throws DomainError delta(-2, X, 0.5)
            @test_throws DomainError delta(2, X, 0.5)
            @test_throws DomainError delta(0.5, X, -2)
            @test_throws DomainError delta(0.5, X, 2)
            @test_throws DomainError delta(0.5, X, 1.0)
        end

    end

end

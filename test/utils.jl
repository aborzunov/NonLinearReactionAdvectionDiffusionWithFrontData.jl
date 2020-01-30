@testset "Utils" begin

    N = 20;
    X = [i/N for i in 0:N];
    Y = u_init.(X);

    @testset "Initial condition" begin

        @test length(Y) == N+1
        @test isapprox(Y[1], -8, rtol=0.001)
        @test isapprox(Y[end], 4, rtol=0.001)
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

        I = 0.0;
        @testset "формула прямоугольников" begin
            for n in 1:N
                δx = X[n+1] - X[n]
                I += d[n] * δx
            end
            @test isapprox(I, 1.0)
        end
        @testset "формула трапеций" begin
            for n in 1:N
                δx = X[n+1] - X[n]
                I += ((d[n] + d[n+1]) / 2.0) * δx
            end
            @test_broken isapprox(I, 1.0)
        end

        @test length( filter( x -> x != 0, d) ) == 1
        @test_throws DomainError delta(-2, X, 0.5)
        @test_throws DomainError delta(2, X, 0.5)
        @test_throws DomainError delta(0.5, X, -2)
        @test_throws DomainError delta(0.5, X, 2)

    end

end

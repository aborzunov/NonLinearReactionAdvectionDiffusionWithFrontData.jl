using NonLinearReactionAdvectionDiffusionWithFrontData;

@testset "Дельта-функция                                   " begin
    N = 20;
    X = [i/N for i in 0:N];

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

    @testset "Проверка области определения дельта-функции" begin
        @test_throws DomainError delta(-2, X, 0.5)
        @test_throws DomainError delta(2, X, 0.5)
        @test_throws DomainError delta(0.5, X, -2)
        @test_throws DomainError delta(0.5, X, 2)
        @test_throws DomainError delta(0.5, X, 1.0)
    end

end

@testset "Гладкая аппроксимация дельта-функции (deltaw/δw)    " begin
    using NonLinearReactionAdvectionDiffusionWithFrontData: deltaw

    N = 100;
    h = 1.0 / N;
    # deltaw принимает Xₙ размера N-1 (внутренние узлы, без граничных точек)
    X_int = [n * h for n in 1:N-1];
    x₀ = 0.5;
    w  = 0.1;

    d = [ deltaw(n, x₀, X_int, N, w) for n in 1:N-1 ]

    @testset "нормировка ≈ 1 (формула трапеций)" begin
        I = 0.0
        for n in 1:N-2
            I += (d[n] + d[n+1]) / 2.0 * h
        end
        @test isapprox(I, 1.0, atol=0.02)
    end

    @testset "носитель — только вблизи x₀" begin
        # Узлы дальше чем w от x₀ должны давать 0
        far = filter(n -> abs(X_int[n] - x₀) > w, 1:N-1)
        @test all(d[far] .== 0.0)
    end

    @testset "максимум — в узле ближайшем к x₀" begin
        imax = argmax(d)
        @test abs(X_int[imax] - x₀) <= w
    end

end

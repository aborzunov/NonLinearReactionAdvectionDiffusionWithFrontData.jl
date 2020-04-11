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

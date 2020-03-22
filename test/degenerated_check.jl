@testset "Поиск нуля функции" begin

    using NonLinearReactionAdvectionDiffusionWithFrontData: find_f_zeros;

    g(x) = (x - 1.1)*(x + 1.1) # функция

    X1 = collect(range(-2, 0, length=20));
    X2 = collect(range(0, 2, length=20));

    @test isapprox( find_f_zeros( g.(X1), X1), -1.1, atol=1e-4)
    @test isapprox( find_f_zeros( g.(X2), X2), 1.1, atol=1e-4)

    h(x) = (x - 1) * (x - 2) * (x - 3)

    # Введем недостаточно густую сетку.
    # Выбранное окно функции для аппроксимации будет немонотонным.
    X = collect( range(-4, 4, length = 50))
    @test_throws ArgumentError find_f_zeros( h.(X), X)

    X = collect( range(-4, 4, length = 100))
    @test isapprox( find_f_zeros( h.(X), X), 1.0, atol=1e-4)

    # Не работает на неотрицательных функция
    j(x) = sin(x)*sin(2x)
    X = collect( range(-0.5, 0.5, length = 100))
    @test_throws ArgumentError find_f_zeros( j.(X), X)

    # Подадим немонотонные значения функции: см первый массив
    @test_throws ArgumentError find_f_zeros( [-1, 1, 3, 2, 4, 5, 6, 8, 7], [-1, 1, 2, 3, 4, 5, 6, 7, 8] )

end

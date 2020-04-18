using NonLinearReactionAdvectionDiffusionWithFrontData;
using NonLinearReactionAdvectionDiffusionWithFrontData: phidetermination, Φ;
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using NonLinearReactionAdvectionDiffusionWithFrontData: shishkin_mesh;
using NonLinearReactionAdvectionDiffusionWithFrontData: apply_on_dynamic_mesh;

using LinearAlgebra;

# Регрессивный тест
# На точных данных, градиент должен быть мал
@testset "Градиент на статической сетке, точные данные    " begin

    u_l(t) = -8
    u_r(t) =  4
    qf(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
                                    # задаче необходимо определить, но при генерации априорной
                                    # информации мы задаем некоторый коэффициент, который,
                                    # собственно, после имея априорную информацию и будем определять.
    ε = 0.2;                        # Малый параметр при старшей производной
    a, b = 0, 1;                    # Область по X
    t₀, T = 0, 0.28;                # Область по T
    N, M = 50, 80;                  # Увеличенное Кол-во разбиений по X, T
    h = (b-a)/N;                    # шаг по X
    τ = (T-t₀)/M;                   # шаг по T
    Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
    Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
    qₙ =      qf.(Xₙ);               # Сеточные значения коэффициента линейного усиления
    ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
    urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
    y₀ = u_init.(Xₙ);               # Начальные условия
    #
    u, XX, TP = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
    nothing #hide

    #' ## Генерация априорной информации
    ϕl      = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                  # Левый вырожденный корень
    ϕr      = phidetermination(qₙ, urₘ, Xₙ, N, Tₘ, M, reverseX = true); # Правый вырожденный корень
    ϕ       = Φ(ϕl, ϕr, N, M);                                          # Серединный корень
    f1_data = f1(ϕ, u, Xₙ, N, M);                                       # Положение переходного слоя
    f2_data = f2(f1_data, u, Xₙ, N, M);                                 # Значение функции на переходном слое
    nothing #hide

    #' ## Решение сопряженной задачи

    Uₙₘ = u;                        # Сохраним старую матрицу
    y₀ = [0.0 for i in 1:N+1];      # Нулевые начальные условия
    ψl = [0.0 for i in 1:M+1];      # Нулевые ГУ
    ψr = [0.0 for i in 1:M+1];      # Нулевые ГУ

    ψ = solve_adjoint(y₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, Uₙₘ, f1_data, f2_data, w = 0.0005)

    grad = J_q(Uₙₘ, ψ, XX, N, Tₘ, M)

    # Значение для этих данных чисто эмпирическое
    # Мы просто проверяем, что градиент достаточно мал
    # "Достаточно" — понятие субъективное
    # Поэтому при тестировании, мы можем лишь удостовериться,
    # Что не произошло регрессии
    @test norm(grad) < 0.05             # L2 норма
    @test norm(grad, Inf) < 0.01        # Максимальный элемент в векторе max(abs.(v))

end

# Регрессивный тест
# На точных данных, градиент должен быть мал
@testset "Градиент на динамической сетке, точные данные   " begin

    u_l(t) = -8;                    # ГУ
    u_r(t) =  4;                    #
    qf(x) = 4*sin(3 * π * x);       # Коэффициент линейного усиления, который в обратной
                                    # задаче необходимо определить, но при генерации априорной
                                    # информации мы задаем некоторый коэффициент, который,
                                    # собственно, после имея априорную информацию и будем определять.
    ε = 0.03;                       # Малый параметр при старшей производной
    a, b = 0, 1;                    # Область по X
    t₀, T = 0, 0.20;                # Область по T
    x_tp = 0.22;                    # Положение переходного слоя
    M = 500;                        # Кол-во разбиений по X, T
    τ = (T-t₀)/M;                   # шаг по T
    Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
    ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
    urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
    nothing #hide

    #' Новое поведение будет реализоваться с помощью передачи функции создания сетки внутрь `solve` в качестве аргумента.
    #' Эта функция должна принимать только один аргумент — положение переходного слоя и формировать соответствующую сетку.
    #' В пакете есть формирование кусочно-равномерной сетки со сгущениями на границе и на переходном слое.
    #' Создадим замыкание этой функции, которое будет иметь нужную сигнатуру.
    mshfrm(x_tp) = shishkin_mesh(a, b, x_tp, ε, 40, 1.0, 1.0, 0.2, 1.0);
    Xₙ = mshfrm(x_tp); ;                            # Сетка по Х
    N = length(Xₙ) - 1                              # Примем за N длину сетки, что получилась в итоге.
    qₙ =  1.0 .* qf.(Xₙ);                            # Сеточные значения коэффициента линейного усиления
                                                    # Измените линейный коэффициент, чтобы посмотреть
                                                    # на поведение сопряженной задачи (а следовательно градиента),
                                                    # когда прямая задача не соответствует наблюдаемым данным.
    u₀ = u_init.(Xₙ, ε=ε, x_tp = x_tp);             # Начальные условия

    u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, create_mesh = mshfrm);
    #' ## Генерация априорной информации
    ϕl      = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                          # Левый вырожденный корень
    ϕr      = phidetermination(qₙ, urₘ, Xₙ, N, Tₘ, M, reverseX = true);         # Правый вырожденный корень
    ϕ       = Φ(ϕl, ϕr, N, M);                                                  # Полуразность вырожденных корней
    ϕ       = apply_on_dynamic_mesh(ϕ, XX, N, M);                               # Аппроксимация на переменную сетку
    ϕl      = apply_on_dynamic_mesh(ϕl, XX, N, M);                              # Аппроксимация на переменную сетку
    ϕr      = apply_on_dynamic_mesh(ϕr, XX, N, M);                              # Аппроксимация на переменную сетку
    f1_data = f1(ϕ, u, XX, N, M);                                               # Положение переходного слоя
    f2_data = f2(f1_data, u, XX, N, M);                                         # Значение функции на переходном слое
    #
    w = 0.000011;                                                               # Априорный параметр апроксимации дельта-функции

    #' ## Решение сопряженной задачи

    Uₙₘ = u;                        # Сохраним старую матрицу
    ψ₀ = [0.0 for i in 1:N+1];      # Нулевые начальные условия
    ψl = [0.0 for i in 1:M+1];      # Нулевые ГУ
    ψr = [0.0 for i in 1:M+1];      # Нулевые ГУ
    ψ = solve_adjoint(ψ₀, XX, N, Tₘ, M, ε, ψl, ψr, qₙ, Uₙₘ, f1_data, f2_data, w = w)

    # Ошибка, если не передадим вспомогательный массив
    @test_throws ArgumentError J_q(Uₙₘ, ψ, XX, N, Tₘ, M)

    X_q = [i/50.0 for i in 0:50]
    grad = J_q(Uₙₘ, ψ, XX, N, Tₘ, M, X_q)

    #=

    plot(Xₙ, grad)


    =#

    # Значение для этих данных чисто эмпирическое
    # Мы просто проверяем, что градиент достаточно мал
    # "Достаточно" — понятие субъективное
    # Поэтому при тестировании, мы можем лишь удостовериться,
    # Что не произошло регрессии
    #
    # Тесты помечены broken, потому что алгоритм неверен
    @test norm(grad) < 0.05             # L2 норма
    @test norm(grad, Inf) < 0.01        # Максимальный элемент в векторе max(abs.(v))

end
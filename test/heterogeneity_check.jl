using  NonLinearReactionAdvectionDiffusionWithFrontData;
using  NonLinearReactionAdvectionDiffusionWithFrontData: heterogeneity;
using  NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using  NonLinearReactionAdvectionDiffusionWithFrontData: apply_on_dynamic_mesh;

@testset "Градиент на точных данных и статической сетке    " begin

    # Создадим априорную информацию и точное решение для точного q(x)
    a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = NonLinearReactionAdvectionDiffusionWithFrontData.dparams();
    #
    u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
    #' ## Генерация априорной информации
    #########################################################################################
    ϕl      = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                      # Левый вырожденный корень
    ϕr      = phidetermination(qₙ, urₘ, Xₙ, N, Tₘ, M, reverseX = true);     # Правый вырожденный корень
    ϕ       = Φ(ϕl, ϕr, N, M);                                              # Серединный корень
    f1_data = f1(ϕ, u, Xₙ, N, M);                                           # Положение переходного слоя
    f2_data = f2(f1_data, u, Xₙ, N, M);                                     # Значение функции на переходном слое
    #########################################################################################
    nothing #hide

    # Функция `heterogeneity` — для внутреннего использования,
    # она принимает аргументы без граничных точек по x.
    Uₙₘ = u[2:N, :];
    X = XX[2:N, :];
    # Создадим замыкания, для удобства
    hmap(z) = [  - heterogeneity(n, m, X[:, m], N, Tₘ, M, Uₙₘ, f1_data, f2_data, z) for n in 1:N-1, m in 1:M+1]
    nnz(arr) = length( findall( x -> ! isapprox(x, 0), arr ))

    # Количество ненулевых элементов — 0
    @test       nnz(hmap(0.0001)) < M/4        # 0, Не рекомендовано
    @test 0 < nnz(hmap(0.0004)) < M/10    # 5, норм!

    # Для визуального контроля, можно использовать следующий код дающий анимацию
    #=
    # {{{
    using Plots, LaTeXStrings;
    pyplot();

    heatmap(h(0.0004)', label=L"- 2 \delta(x - f_1(t)) ( u^s(x,t) - f_2(t))")

    d = [ deltaw( n, f1[m], Xₙ[2:N], N) for n in 1:N-1, m in 1:M+1];    # Значения первого множителя
    uf = [ ( Uₙₘ[n+1, m] - f2[m] ) for n in 1:N-1, m in 1:M+1];         # Значения второго множителя

    yl = extrema( [H; d; uf]  );
    a = Animation()
    for m in 1:M
        nnd = findfirst( x-> x!= 0, d[:,m]);                # Первый ненулевой элемент на m-ом шаге по T
        nnh = findfirst( x-> x!= 0, H[:,m]);                # Первый ненулевой элемент на m-ом шаге по T
        ym = d[nnd, m] > H[nnh, m] ? d[nnd,m] : H[nnh, m]   # Выбор наибольшего ненулевого элемента
                                                            # чтобы нарисовать красивую вертикальную лиинию
        scatter(X,  H[:, m],    label=L"- 2 \delta(x - f_1(t)) ( u^s(x,t) - f_2(t))", ylims = yl)
        scatter!(X, d[:, m],    label=L"\delta(x - f_1(t))")
        scatter!(X, uf[:, m],   label=L"( u^s(x,t) - f_2(t))")
        plot!([Xₙ[nnd+1], Xₙ[nnd+1]], [ym, 0], line=:dash, label="") # Вертикальная линия
        frame(a)
    end
    g = gif(a, "heterogenyry.gif")

    # Карта поверхности неоднородности
    heatmap(X, Tₘ, H',  xlabel = "x", ylabel="t", title=L"- 2 \delta(x - f_1(t)) ( u^s(x,t) - f_2(t))")
    savefig("heterogeneity.png")
    # }}}
    =#

end

@testset "Градиент на точных данных и динамической сетке   " begin

    using NonLinearReactionAdvectionDiffusionWithFrontData;
    using NonLinearReactionAdvectionDiffusionWithFrontData: Φ, apply_on_dynamic_mesh;
    using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;

    a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀, mshfrm = NonLinearReactionAdvectionDiffusionWithFrontData.dparams_nonuniform();
    #
    u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, create_mesh = mshfrm);
    #' ## Генерация априорной информации
    #########################################################################################
    ϕl      = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                          # Левый вырожденный корень
    ϕr      = phidetermination(qₙ, urₘ, Xₙ, N, Tₘ, M, reverseX = true);         # Правый вырожденный корень
    ϕ       = Φ(ϕl, ϕr, N, M);                                                  # Полуразность вырожденных корней
    ϕ       = apply_on_dynamic_mesh(ϕ, XX, N, M);                               # Аппроксимация на переменную сетку
    ϕl      = apply_on_dynamic_mesh(ϕl, XX, N, M);                              # Аппроксимация на переменную сетку
    ϕr      = apply_on_dynamic_mesh(ϕr, XX, N, M);                              # Аппроксимация на переменную сетку
    f1_data = f1(ϕ, u, XX, N, M);                                               # Положение переходного слоя
    f2_data = f2(f1_data, u, XX, N, M);                                         # Значение функции на переходном слое
    #########################################################################################


    using  NonLinearReactionAdvectionDiffusionWithFrontData: heterogeneity
    # Функция `heterogeneity` — для внутреннего использования,
    # она принимает аргументы без граничных точек по x.
    Uₙₘ = u[2:N, :];
    X = XX[2:N, :];
    hmap(z) = [  - heterogeneity(n, m, X[:, m], N, Tₘ, M, Uₙₘ, f1_data, f2_data, z) for n in 1:N-1, m in 1:M+1]
    nnz(arr) = length( findall( x -> ! isapprox(x, 0), arr ));

    @test nnz(hmap(0.00002)) > M * 0.05         # Много ненулевых элементов,
                                                # значит итерационный процесс продолжится и уйдет от решения
    @test nnz(hmap(0.0000001)) < M * 0.05        # 0, Не рекомендовано
    @test 0 < nnz(hmap(0.000003)) < M * 0.01    # 4, Должно быть норм!

end

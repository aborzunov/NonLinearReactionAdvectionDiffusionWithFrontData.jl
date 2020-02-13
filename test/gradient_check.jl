@testset "Градиент на точных данных" begin

    # Создадим априорную информацию и точное решение для точного q(x)
    #########################################################################################
    u_l(t) = -8.0
    u_r(t) =  4.0
    q(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления, который в обратной
    # задаче необходимо определить, но при генерации априорной
    # информации мы задаем некоторый коэффициент, который,
    # собственно, после имея априорную информацию и будем определять.
    ε = 0.2;                        # Малый параметр при старшей производной
    a, b = 0, 1;                    # Область по X
    t₀, T = 0, 0.36;                # Область по T
    N, M = 100, 80;                 # Кол-во разбиений по X, T
    h = (b-a)/N;                    # шаг по X
    τ = (T-t₀)/M;                   # шаг по T
    Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
    Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
    qₙ =      q.(Xₙ);               # Сеточные значения коэффициента линейного усиления
    ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
    urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
    u₀ = u_init.(Xₙ);               # Начальные условия
    nothing #hide

    u = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
    ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
    ϕr = phidetermination(qₙ, urₘ, reverse(Xₙ), N, Tₘ, M);                      # Нужно подать инвертированную сетку
    ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
    ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(ϕl, ϕr, N, M);       # Серединный корень
    f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);   # Положение переходного слоя
    f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);  # Значение функции на переходном слое
    nothing #hide
    #########################################################################################

    using  NonLinearReactionAdvectionDiffusionWithFrontData: heterogenety
    # Функция `heterogenety` — для внутреннего использования,
    # она принимает аргументы без граничных точек по x.
    Uₙₘ = u[2:N, :];
    X = Xₙ[2:N];
    H = [  - heterogenety(n, m, X, N, Tₘ, M, Uₙₘ, f1, f2) for n in 1:N-1, m in 1:M+1];

    # Количество ненулевых элементов меньше 25%
    @test length( findall( x -> x!= 0, H) ) < M/4

    # Для визуального контроля, можно использовать следующий код
    #=
    using Plots, LaTeXStrings

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
    savefig("heterogenety.png")
    =#
end

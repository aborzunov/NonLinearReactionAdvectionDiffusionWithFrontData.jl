# # Подготовка данных для тестирования алгоритма выбора начального приближения

# Здесь мы готовим данные, для тестирования алгоритмов.
# Подразумевается, что скрипт будет исполнятся очень редко. Все свои результаты
# он запишет в соответствующие файлы в текущей папке.

using NonLinearReactionAdvectionDiffusionWithFrontData;
using DelimitedFiles;
using Plots; pyplot();


# Выберем параметры таким образом, чтобы решение точно было хорошим.
# Здесь мы генерируем экспериментальные данные, поэтому решим прямую задачу
# всего раз для каждого случая.

x_tp    = 0.05;         # Стартовое местоположение фронта
T_end   = 0.48;         # Регулируем конечное местоположение фронта
Nx      = 18000;         # Число интервалов по ``X``
Mt      = 27000;         # Число интервалов по ``T``


# Сгенерируем данные для набора малых параметров
for eps in [0.03, 0.01, 0.001]
    @info "Генерируем данные для ε=$(eps)"
    ε       = eps;          # Крутизна фронта

    a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
                                                            ε = ε,
                                                            Nx = Nx,
                                                            Mt = Mt,
                                                            T_end = T_end,
                                                            qfunc = x -> sin(3pi*x)
                                                            );
    u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
    ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);

    # Убеждаемся в корректности прямого решения на глаз и пишем данные в txt
    directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения")
    savefig("direct_eps_$(ε)_sin.png");
    writedlm("tm_f1_m_eps_$(ε)_sin.txt",[Tₘ f1_data])
    writedlm("xn_qn_eps_$(ε).txt",[Xₙ qₙ])

    # Два гаусиана разной высоты, которые частично перекрываются
    gauss_init(x) = 2 * exp(-((x-5.0/13.0)^2)/(1.5/13.0^2)) +
        1.2 * exp((-(x-8.0/13.0)^2)/(3.0/13.0^2));
    qₙ = [gauss_init(x) for x in Xₙ]
    plot(Xₙ, qₙ)

    u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
    ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);

    # Убеждаемся в корректности прямого решения на глаз и пишем данные в txt
    directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения")
    savefig("direct_eps_$(ε)_gauss.png");
    writedlm("tm_f1_m_eps_$(ε)_gauss.txt",[Tₘ f1_data])
    writedlm("xn_qn_eps_$(ε)_gauss.txt",[Xₙ qₙ])
end

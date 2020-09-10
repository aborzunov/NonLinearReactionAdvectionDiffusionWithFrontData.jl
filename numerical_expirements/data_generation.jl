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
T_end   = 0.58         # Регулируем конечное местоположение фронта
Nx      = 15000;         # Число интервалов по ``X``
Mt      = 30000;        # Число интервалов по ``T``
α       = 0.005;        # Параметр регуляризации в нахождении начального приближения

# Генерируем данные для набора малых параметров

plts = Any[];           # Вектор для хранения графиков Plots

for (eps, x_tp) in [(0.03, 0.04), (0.01, 0.03), (0.001,0.02)]
    @info "Генерируем данные для ε = $(eps), x_tp = $(x_tp)"
    ε       = eps;          # Крутизна фронта
#+
    a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams(x_tp = x_tp,
                                                            ε = ε,
                                                            Nx = Nx,
                                                            Mt = Mt,
                                                            T_end = T_end,
                                                            qfunc = x -> sin(3pi*x)
                                                            );
    u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
    ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
#+
    # Убеждаемся в корректности прямого решения на глаз и пишем данные в txt
    directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения, \\varepsilon = $(ε)")
    savefig("direct_eps_$(ε)_sin.png");
    writedlm("tm_f1_m_eps_$(ε)_sin.txt",[Tₘ f1_data])
    writedlm("xn_qn_eps_$(ε).txt",[Xₙ qₙ])
#+
    # Находим начальные приближения
    #  q_guess = initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ, α);
    #  p = plot(Xₙ, qₙ, label="Истинное",  title = "Начальное приближение, \\varepsilon = $(ε)");
    #  p = plot!(Xₙ, q_guess, label="Найденное")
    #  savefig("guess_eps_$(ε)_sin.png");
    #  push!(plts, p);
#+
    # Два гаусиана разной высоты, которые частично перекрываются
    gauss_init(x) = 2 * exp(-((x-5.0/13.0)^2)/(1.5/13.0^2)) +
        1.2 * exp((-(x-8.0/13.0)^2)/(3.0/13.0^2));
    qₙ = [gauss_init(x) for x in Xₙ]
    plot(Xₙ, qₙ)
#+
    u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ);
    ϕl, ϕr, ϕ, f1_data, f2_data = generate_obs_data(u, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ);
#+
    # Убеждаемся в корректности прямого решения на глаз и пишем данные в text
    directP = draft(u, Xₙ, N, Tₘ, M, title = "Эскиз прямого решения, \\varepsilon = $(ε)")
    savefig("direct_eps_$(ε)_gauss.png");
    writedlm("tm_f1_m_eps_$(ε)_gauss.txt",[Tₘ f1_data])
    writedlm("xn_qn_eps_$(ε)_gauss.txt",[Xₙ qₙ])
#+
    # Находим начальные приближения
    #  q_guess = initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ, α);
    #  p = plot(Xₙ, qₙ, label="Истинное",  title = "Начальное приближение, \\varepsilon = $(ε)");
    #  p = plot!(Xₙ, q_guess, label="Найденное")
    #  savefig("guess_eps_$(ε)_gauss.png");
    #  push!(plts, p);

end

# График
plts[1]

# График
plts[2]

# График
plts[3]

# График
plts[4]

# График
plts[5]

# График
plts[6]

# Файл с некоторыми скриптами,
# Которые были полезны во время работы, но не могут быть внесены на
# постоянной основе в скрипты, которые представляют реальный интерес.


# В одном из экспериментов потребовалось запустить вычисления с сохраненных
# ранее сериализованных результатов, здесь их чтение.
Js = deserialize("Js_zero.jld")
Qs = deserialize("QQs_zero.jld")
qs = Qs_old[:, end];

Js = deserialize("Js.jld")
Qs = deserialize("QQs.jld")
qs = Qs_old[:, end];

Js_old_1 = Js;
Qs_old_1 = Qs;

q₀ = Qs[:,end];

SS = 100000;

SS = 3000;

# Отрисовка больших gif
frames = collect(1:div(S, 500):S);
make_minimzation_gif(
    Js, Qs, qₙ, Xₙ,
    name = "Minimization_noguess.gif",
    frames_to_write = frames,
    β = β
)

# Отрисовка больших gif
frames = collect(1:div(S, 500):S);
make_minimzation_gif(
    Js, Qs, qₙ, Xₙ,
    name = "Minimization_withguess.gif",
    frames_to_write = frames,
    β = β
)

initial_guessP = plot(Xₙ, q₀, xlabel = "X", ylabel = "q(x)", label="Начальное приближение");
initial_guessP = plot!(Xₙ, Qs[:, 3000], label="Найденное решение")
savefig("results.svg");

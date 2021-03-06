# # Пример поиска начального приближения
#
# На этой странице мы делаем `include` примера, который находится
# в корне пакета в каталоге `"examples/"`.
#
# [Literate.jl](https://github.com/fredrikekre/Literate.jl) включает
# текст этих примеров и в md файл.

include("examples/example_initial_guess.jl")

# ### Визуализация

# Набросаем эскиз решения прямой задачи
psol = plot(title="Эскиз решения");
plot!(Xₙ, u[:, 1], label="u(x, 0)");
plot!(Xₙ, u[:, div(end,2)], label="u(x, T/2)");
plot!(Xₙ, u[:, end], label="u(x, T)");
psol

# Найденное начальное приближение.
#
# Начальное приближение не может быть определено в точках сетки ``X_N``,
# которые не пересек переходный слой, поэтому там, начальное приближение
# продолжено константой.
plot(Xₙ, qₙ, label="Истинное")
plot!(Xₙ, q_guess, label="Начальное приближение")

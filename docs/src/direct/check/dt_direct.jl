
# ## Непосредственная реализация проверки
#
using Test;
using LaTeXStrings;
using Plots;

include("test/direct_check.jl")

# Посмотрим на результат решения в сравнении с аналитическим решением

d = [missing, missing];
dd = [missing missing; missing missing];
make_gif(u, Xₙ, Tₘ, dd, dd, dd, d, d, u_model; name = "dicrect_check.gif")

# Найдем абсолютную погрешность численного решения

err = u .- u_model
heatmap(Xₙ, Tₘ, err', xlabel=L"X_n", ylabel=L"T_m", title="Absolute Error", size=(1200, 800))

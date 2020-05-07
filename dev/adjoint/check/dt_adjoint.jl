
## Непосредственная реализация  проверки сопряженной задачи
#
using Test;
using LaTeXStrings;
using Plots;

include("test/adjoint_check.jl")

# Посмотрим на результат решения в сравнении с аналитическим решением

d = [missing, missing];
dd = [missing missing; missing missing];
make_gif(ψ, Xₙ, Tₘ, dd, dd, dd, d, d, ψ_model; name = "direct_check.gif")

# Найдем абсолютную погрешность численного решения

err = ψ .- ψ_model
heatmap(Xₙ, Tₘ, err', xlabel=L"X_n", ylabel=L"T_m", title="Absolute Error", size=(1200, 800))


# Проверка решения сопряженной задачи

Зададим следующую пробную функцию: ``g(x,t) = (1-2t) \sin(\pi x)`` и найдем её производные:

* ``\frac{\partial g}{\partial x} = \pi (1-2t) \cos(\pi x)``
* ``\frac{\partial^2 g}{\partial x^2} = - \pi^2 (1 - 2t) \sin(\pi x)``
* ``\frac{\partial g}{\partial t} = -2 \sin(\pi x)``.

Подставим найденные производные в исходное уравнение (см.[Сопряженная задача](@ref)) и перенеся все в правую часть, определим ``g_d``.
```math
g_d(x,t) = 2 \sin(\pi x) - ( - \varepsilon \pi^2 (1 - 2t) \sin(\pi x) )+
u(x,t) \pi (1 - 2t) \cos(\pi x) + q(x) (1 - 2t) \sin(\pi x) - 2 \delta( x - f_1(t) ) (u(x,t) - f_2(t))
```

Таким образом, пробная функция ``g`` будет являться решением уравнения

```math
\left\{
\begin{aligned}
    & \frac{\partial \psi}{\partial t} = - \varepsilon\frac{\partial^2 \psi}{\partial x^2} + u \frac{\partial \psi}{\partial x} + q(x)\,\psi  -\\
    & \qquad  - 2\delta(x - f_1(t))(u(x,t) - f_2(t)) - g_d, \quad x \in (0,1), \quad t \in (0,T), \\
    &\psi(0,t) = 0, \quad \psi(1,t) = 0, \quad t \in (0,T), \\
    &\psi(x,T) = 0, \qquad x \in [0,1].
\end{aligned}
\right.
```


Такая проверка корректности решения применяется в юнит тесте `"tests/adjoint_check.jl"`.
Файл содержит один `@testset`, внутри него реализовано решение вышеописанной системы, проверка его корректности
через `@test`. А так же, `@testset` возвращает `ψ, ψ_model, Xₙ, Tₘ`, что соответствует решению, аналитическому
решению, сетке по X, T.

```@example test_direct_check
using NonLinearReactionAdvectionDiffusionWithFrontData, Test, ForwardDiff
using LaTeXStrings, Plots
ψ, ψ_model, Xₙ, Tₘ = include("../../../test/adjoint_check.jl")
nothing #hide
```

Анимация решения нарисуем с неравномерным течением времени: первые 80 кадров, а потом каждый десятый.
```@example test_direct_check
d = [missing, missing];
@info "$( splitdir(@__FILE__)[2] ) Рисует решение сопряженной задачи на модельной функции."
make_gif(ψ, Xₙ, Tₘ[end:-1:1], d, d, d, d, ψ_model;
            name="adjoint_check.gif", label="\\psi", frames_to_write=[1:80; 81:10:length(Tₘ)], convert2mp4=true)
```

```@example test_direct_check
using LaTeXStrings, Plots
err = ψ .- ψ_model
heatmap(Xₙ, Tₘ, err', xlabel=L"X_n", ylabel=L"T_m", title="Absolute Error", size=(1200, 800))
```

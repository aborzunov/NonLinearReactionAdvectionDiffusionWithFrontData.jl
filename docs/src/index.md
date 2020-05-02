```@meta
DocTestSetup = quote
    using NonLinearReactionAdvectionDiffusionWithFrontData
    using Plots
end

```

# Главная

## Постановка обратной задачи

Рассмотрим прямую задачу для сингулярно возмущенного уравнения типа Бюргерса:
```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 u}{\partial x^2} - \frac{\partial u}{\partial t} = -u \frac{\partial u}{\partial x} +  q(x)\,u, \quad x \in (0,1), \quad t \in (0,T), \\
    &u(0,t) = u_{left}(t), \quad u(1,t) = u_{right}(t), \quad t \in (0,T), \\
    &u(x,t) = u_{init}(x), \qquad x \in [0,1].
\end{aligned}
\right.
```

Решение этой задачи имеет движущийся слой, положение которого во времени
описывает ``x = x_{t.p}(t)``.

Обратная задача состоит в определении коэффициента линейного усиления ``q(x)``,
``x \in [0,1]``, по известной дополнительной информации о положении переходного
слоя и значения функции на переходном слое:
```math
x_{t.p} = f_1(t), \qquad u(x_{t.p}(t),t) = f_2(t), \qquad t \in [0, T].
```

## Содержание

### Прямая задача

```@contents
Pages = [
    "direct/direct.md",
    "direct/experimental_data.md",
    "generated/docexample_direct.md",
    "generated/direct_check.md",
]
Depth = 1
```

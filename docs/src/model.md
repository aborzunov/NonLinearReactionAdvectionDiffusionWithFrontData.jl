# Тестирование

## Проверка корректности решения прямой задачи

Зададим следующую пробную функцию: ``g(x,t) = (1-2t) \sin(\pi x)`` и найдем её производные:

* ``\frac{\partial g}{\partial x} = \pi (1-2t) \cos(\pi x)``
* ``\frac{\partial^2 g}{\partial x^2} = - \pi^2 (1 - 2t) \sin(\pi x)``
* ``\frac{\partial g}{\partial t} = -2 \sin(\pi x)``.

Подставим найденные производные в исходное уравнение (см.[Постановка задачи](@ref)), таким образом определив ``g_d`` (от англ. discrepancy).
```math
g_d(x,t) = 2 \sin(\pi x) - \varepsilon \pi^2 (1 - 2t) \sin(\pi x) +
\pi (1 - 2t)^2 \sin(\pi x) \cos(\pi x) - q(x) (1 -2t) \sin(\pi x)
```

Таким образом, пробная функция ``g`` будет являться решением уравнения

```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 u}{\partial x^2} - \frac{\partial u}{\partial t} = -u \frac{\partial u}{\partial x} +  q(x)\,u + g_d(x,t), \quad x \in (0,1), \quad t \in (0,T), \\
    &u(0,t) = u_{left}(t), \quad u(1,t) = u_{right}(t), \quad t \in (0,T), \\
    &u(x,t) = u_{init}(x), \qquad x \in [0,1].
\end{aligned}
\right.
```

Найдем решение следующим образом:
```math
    \begin{aligned}
        &\mathbf{\textbf{y}}(t_{m + 1}) = \mathbf{\textbf{y}}(t_m) + (t_{m + 1} - t_m) \, \mathrm{Re} \, \mathbf{\textbf{w}}_1 \,\\
        &\left[\mathbf{\textbf{E}} - \dfrac{1 + i}{2} \, (t_{m + 1} - t_m) \, \mathbf{\textbf{f}}_\mathbf{\textbf{y}}\Big(\mathbf{\textbf{y}}(t_m),t_m\Big)\right] \, \mathbf{\textbf{w}}_1 = \\
        &\qquad\qquad\qquad\qquad\quad = \mathbf{\textbf{f}} \, \Big(\mathbf{\textbf{y}}(t_m),\frac{t_{m + 1} + t_m}{2}\Big) + g_d(\mathbf{X_n},\frac{t_{m + 1} + t_m}{2}).
    \end{aligned}
```
где ``g_d(\mathbf{X_n})`` — значения ``g_d`` на сетке ``x_1, x_x, \ldots, x_N``, т.е. на сетке ``\mathbf{X\_n}`` без граничных точек.

Такая проверка корректности решения применяется в юнит тесте `"tests/direct_with_model.jl"`.

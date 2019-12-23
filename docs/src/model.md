# Проверка корректности решения с помощью модельной функции.

В качестве пробной функции возьмем ``g = (1-2t) \sin(\pi x)``.

* ``\frac{\partial g}{\partial x} = \pi (1-2t) \cos(\pi x)``
* ``\frac{\partial^2 g}{\partial x^2} = - \pi^2 (1 - 2t) \sin(\pi x)``
* ``\frac{\partial g}{\partial t} = -2 \sin(\pi x) ``

Подставим найденные производные в исходное уравнение (см.[Постановка задачи](@ref)), таким образом определив ``g_d`` (от англ. discrepancy).
```math
g_d = 2 \sin(\pi x) - \varepsilon \pi^2 (1 - 2t) \sin(\pi x) +
\pi (1 - 2t)^2 \sin(\pi x) \cos(\pi x) - q(x) (1 -2t) \sin(\pi x)
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

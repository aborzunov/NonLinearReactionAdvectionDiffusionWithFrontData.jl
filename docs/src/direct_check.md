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

```@example test
using Test, ForwardDiff
using NonLinearReactionAdvectionDiffusionWithFrontData

u_l(t) = 0
u_r(t) = 0
ε = 0.1;
a, b = 0, 1;
t₀, T = 0, 1;
N, M = 40, 80;
h = (b-a)/N;
τ = (T-t₀)/M;
Xₙ = [a  + n*h for n in 0:N];
Tₘ = [t₀ + m*τ for m in 0:M];
q(x) = sin(3 * π * x);
u = zeros(M+1, N+1);

# Зададим модельную функцию и невязку, получаемую после подстановки `g` в исходное уравнение
g(x, t) = (1 - 2t)*sin(π*x);
g_d(x,t) =  - ε * π^2 * (1 - 2t) * sin(π * x) + π * (1 - 2t)^2 * sin(π * x) * cos(π * x) - q(x) * (1 -2t) * sin(π * x) + 2 * sin(π * x);

y = g.( (Xₙ[n] for n in 2:N), 0 );
qₙ = [ q(x) for x in Xₙ[2:N] ];

# Модельное решение найденное с помощью известного аналитического решения
model = [ g(x,t) for x in Xₙ, t in Tₘ];

# Создадим функцию, которая будет вычислять вектор правой части с добавлением невязки
RP(y, t, Xₙ, N, ε, u_l, u_r, qₙ) = f(y, t, Xₙ, N, ε, u_l, u_r, qₙ) - g_d.(Xₙ[2:N], t)
# Хоть мы и конструируем якобиан с помощью автоматического дифференцирования, примите во внимание, что
# Якобиан ``f_y`` при добавлении `g_d` останется без изменений, т.к. `g_d` зависит только от ``x,t``.
# То, что он не зависит от добавления `g_d` можно убедиться изменением порядка этих двух строк, ну а так же на бумаге.
j(y, t, Xₙ, N, ε, u_l, u_r, qₙ) = ForwardDiff.jacobian( z -> RP(z, t, Xₙ, N, ε, u_l, u_r, qₙ), y)

u= solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, qₙ, RP, j);
@test isapprox(model, u, rtol = 1E-3)
```

Сохраним анимацию решения на модельной функции.
```@example test
d = [missing, missing]
make_gif(u, Xₙ, Tₘ, d, d, d, d, g; frame_skip = div(M,30), frames_to_write=80, name="solution_model.gif")
nothing #hide
```

![](assets/solution_model.mp4)

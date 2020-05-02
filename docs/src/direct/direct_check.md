# Проверка корректности решения прямой задачи

Применим метод пробной функции или контрольной невязки.
Это позволит сравнивать решение, получаемое нашей схемой, с аналитическим
решением, даже если мы не можем найти аналитическое решение исходного
уравнения.

Идея проверки заключается в выборе некоторой известной пробной функции,
например, ``g(x,t)``, которая не является решением исходного уравнения.
Подставив пробную функцию уравнение, оно не обращается в тождество.
Если преобразовать уравнение путем вычитания невязки, полученной в прошлом
действии, то мы получим новое уравнение, для которого уже ``g(x,t)`` ---
решение.

После этого, мы решаем решаем модифицированное уравнение с помощью нашей схемы,
в которой нужно лишь подправить функцию возвращающую вектор правой части ---
нужно дополнительно вычесть невязку, формулу для которой мы нашли аналитически.

## Алгоритм

Зададим следующую пробную функцию ``g(x,t) = (1-2t) \sin(\pi x)``.
Найдем её производные:

* ``\frac{\partial g}{\partial x} = \pi (1-2t) \cos(\pi x)``
* ``\frac{\partial^2 g}{\partial x^2} = - \pi^2 (1 - 2t) \sin(\pi x)``
* ``\frac{\partial g}{\partial t} = -2 \sin(\pi x)``.

Подстановка пробной функции в исходной уравнение
(см. [Постановка обратной задачи](@ref)), позволит нам определить
``g_d``
```math
g_d(x,t) = 2 \sin(\pi x) - \varepsilon \pi^2 (1 - 2t) \sin(\pi x) +
\pi (1 - 2t)^2 \sin(\pi x) \cos(\pi x) - q(x) (1 -2t) \sin(\pi x)
```

Таким образом, пробная функция ``g`` будет являться решением уравнения

```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 u}{\partial x^2} -
    \frac{\partial u}{\partial t} = -u \frac{\partial u}{\partial x} +
    q(x)\,u - g_d(x,t), \quad x \in (0,1), \quad t \in (0,T], \\
    &u(0,t) = u_{left}(t), \quad u(1,t) = u_{right}(t), \quad t \in (0,T], \\
    &u(x,t) = u_{init}(x), \qquad x \in [0,1], t = 0.
\end{aligned}
\right.
```

Найдем решение следующим образом:
```math
    \begin{aligned}
        &\mathbf{\textbf{y}}(t_{m + 1}) = \mathbf{\textbf{y}}(t_m) +
        (t_{m + 1} - t_m) \, \mathrm{Re} \, \mathbf{\textbf{w}}_1 \,\\
        &\left[\mathbf{\textbf{E}} - \dfrac{1 + i}{2} \, (t_{m + 1} - t_m) \,
        \mathbf{\textbf{f}}_\mathbf{\textbf{y}}\Big(\mathbf{\textbf{y}}(t_m),
        t_m\Big)\right] \, \mathbf{\textbf{w}}_1 = \\
        &\qquad\qquad\qquad\qquad\quad = \mathbf{\textbf{f}} \,
        \Big(\mathbf{\textbf{y}}(t_m),\frac{t_{m + 1} + t_m}{2}\Big) -
        g_d(\mathbf{X_n},\frac{t_{m + 1} + t_m}{2}).
    \end{aligned}
```
где ``g_d(\mathbf{X_n})`` — значения ``g_d`` на сетке
``x_1, x_x, \ldots, x_N``, т.е. на сетке ``\mathbf{X_n}`` без граничных точек.
Решение будем находить с помощью функции [`solve`](@ref), но перед этим
сконструируем функцию правой части и её якобиан и передадим их в качестве
аргументов.

Такая проверка корректности решения применяется в unit тесте
`"tests/direct_with_model.jl"`. Здесь мы сделаем `include`,
[Literate.jl](https://github.com/fredrikekre/Literate.jl) включит
текст теста в md файл.

Перед запуском, необходимо подключить *все* пакеты, что используются при
тестировании.

```@meta
EditURL = "<unknown>/src/examples/dt_direct.jl"
```

## Непосредственная реализация проверки

```@example doctest_direct
using Test;
using LaTeXStrings;
using Plots;
nothing #hide
```

Тест проверяет корректность решения прямой задачи. Алгоритм описан в /docs/src/direct/direct_check.md
Возвращает решение, аналитическое решение, сетку по X, сетку по T.

```@example doctest_direct
using NonLinearReactionAdvectionDiffusionWithFrontData;
using NonLinearReactionAdvectionDiffusionWithFrontData: phidetermination, Φ;
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using ForwardDiff;
nothing #hide
```

Зададим параметры для прямой задачи

```@example doctest_direct
u_l(t) = 0;                     # ГУ удовлетворяющие модельной функции
u_r(t) = 0;                     # ГУ удовлетворяющие модельной функции
qf(x) = 4*sin(3 * π * x);        # Коэффициент линейного усиления
ε = 0.2;                        # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 1;                   # Область по T
N, M = 50, 80;                  # Кол-во разбиений по X, T
h = (b-a)/N;                    # шаг по X
τ = (T-t₀)/M;                   # шаг по T
Xₙ = [a  + n*h for n in 0:N];   # Сетка по Х
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
qₙ =      qf.(Xₙ);               # Сеточные значения коэффициента линейного усиления
ulₘ= u_l.(Tₘ);                  # Сеточные значения левого  ГУ
urₘ= u_r.(Tₘ);                  # Сеточные значения правого ГУ
nothing #hide
```

Зададим модельную функцию и невязку, получаемую после подстановки `g` в исходное уравнение

```@example doctest_direct
function g(x, m)
    t = Tₘ[m]
    (1 - 2t)*sin(π*x)
end
function g_d(x::Real, m::Int)
    t = Tₘ[m];
    - ε * π^2 * (1 - 2t) * sin(π * x) + π * (1 - 2t)^2 * sin(π * x) * cos(π * x) - qf(x) * (1 -2t) * sin(π * x) + 2 * sin(π * x)
end

y₀ = g.(Xₙ, 1);               # Начальные условия
nothing #hide
```

Модельное решение найденное с помощью известного аналитического решения

```@example doctest_direct
u_model = [ g(x, m) for x in Xₙ, m in 1:M+1];
nothing #hide
```

Создадим функцию, которая будет вычислять вектор правой части с добавлением невязки

```@example doctest_direct
function RP(y, m, Xₙ, N, ε, ulₘ, urₘ, qₙ)
    d = [ g_d(x, m) for x in Xₙ[2:N] ]
    NonLinearReactionAdvectionDiffusionWithFrontData.directRP(y, m, Xₙ, N, ε, ulₘ, urₘ, qₙ) - d
end
```

Хоть мы и конструируем якобиан с помощью автоматического дифференцирования, примите во внимание, что
Якобиан ``f_y`` при добавлении `g_d` останется без изменений, т.к. `g_d` зависит только от ``x,t``.
То, что он не зависит от добавления `g_d` можно убедиться изменением порядка этих двух строк, ну а так же на бумаге.

```@example doctest_direct
j(y, m, Xₙ, N, ε, ulₘ, urₘ, qₙ) = ForwardDiff.jacobian( z -> RP(z, m, Xₙ, N, ε, ulₘ, urₘ, qₙ), y)
```

С использованием автоматического дифференцирования

```@example doctest_direct
u, XX, TP = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, RP, j);
@test all(isapprox.(u_model, u, atol = 0.01))
```

С использованием трехдиагонального якобиана

```@example doctest_direct
u, XX, TP = solve(y₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, RP, NonLinearReactionAdvectionDiffusionWithFrontData.∂DRP_∂y);
@test all(isapprox.(u_model, u, atol = 0.01))

(u, u_model, Xₙ, Tₘ)
```

Посмотрим на результат решения в сравнении с аналитическим решением

```@example doctest_direct
d = [missing, missing];
dd = [missing missing; missing missing];
make_gif(u, Xₙ, Tₘ, dd, dd, dd, d, d, u_model; name = "dicrect_check.gif")
```

Найдем абсолютную погрешность численного решения

```@example doctest_direct
err = u .- u_model
heatmap(Xₙ, Tₘ, err', xlabel=L"X_n", ylabel=L"T_m", title="Absolute Error", size=(1200, 800))
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


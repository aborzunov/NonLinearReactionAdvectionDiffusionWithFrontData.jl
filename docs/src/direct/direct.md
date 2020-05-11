# Прямая задача

## [Постановка прямой задачи](@id direct-formulation)

Прямая задача представляет из себя параболическое уравнение типа
реакция-адвекция-диффузия для отрезка ``[0,1]`` с граничными условиями
первого рода, где искомая функция ``u(x,y)`` подлежит определению на
полуинтервале ``(0, T]``
```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 u}{\partial x^2} -
    \frac{\partial u}{\partial t} = -u \frac{\partial u}{\partial x} +
    q(x)\,u, \quad x \in (0,1), \quad t \in (0,T], \\
    &u(0,t) = u_l(t), \quad u(1,t) = u_r(t), \quad t \in (0,T], \\
    &u(x,t) = u_i(x), \qquad x \in [0,1], t=0.
\end{aligned}
\right.
```

## Схема решения

Схема решения будет подробно объяснена на равномерной сетке.
Формулы для неравномерной сетки будет даны без объяснений.
В обоих случаях будет использоваться метод жестких прямых для сведения
уравнения в частных производных к системе обыкновенных дифференциальных
уравнений. После чего, последняя будет решаться с помощью одностадийной схемы Розенброка с
комплексными коэффициентами (CROS1).

Введём равномерную сетку ``X_N`` с ``N`` интервалами по пространственной
переменной ``x`` с шагом `` h = 1/N``: ``X_N = \lbrace x_n, 0 \le n \le N:
x_n = n h \rbrace``.
Произведём конечно-разностную аппроксимацию пространственных производных,
перенесём всё в правую часть, так, чтобы слева осталась только производная по
времени и получим дифференциально-алгебраическую систему, где количество
неизвестных ``u_n \equiv u_n(t) \equiv u(x_n, t), n = \overline{0, N}``
равняется ``N+1``.

```math
\left\{
    \begin{aligned}
    & \frac{d u}{d t} = \varepsilon \frac{u_{n+1} - 2 u_n + u_{n-1}}{h^2} +
    u_n \frac{u_{n+1} - u_{n-1}}{2h} - q_n u_n, \quad n = \overline{1, N-1}, t \in
    (0, T] \\
    & u_0 = u_l(t), u_N = u_r(t), \quad t \in (0, T], \\
    & u(x_n, 0) = u_i(x_n), \quad x \in [0, 1], t =0.
    \end{aligned}
\right.
```

Теперь сведём полученную дифференциально-алгебраическую систему к системе
обыкновенных дифференциальных уравнений их ``N-1`` уравнений с ``N-1``
неизвестной путем подстановки алгебраических
выражений для ``u_0, u_N`` в правую часть дифференциальных уравнений
```math
\left\{
\begin{aligned}
    & \frac{du}{dt} = \frac{u_{2} - 2 u_1 + u_l(t)}{h^2} + u_1
    \frac{u_2 - u_l(t)}{2h} - q_1 u_1, \quad t \in (0, T] \\
    & \frac{d u}{d t} = \varepsilon \frac{u_{n+1} - 2 u_n + u_{n-1}}{h^2} +
    u_n \frac{u_{n+1} - u_{n-1}}{2h} - q_n u_n, \quad n = \overline{2, N-2},
    t \in (0, T] \\
    & \frac{du}{dt} = \varepsilon \frac{u_r(t) - 2 u_{N-1} + u_{N-2}}{h^2} +
    u_{N-1} \frac{u_r(t) - u_{N-2}}{2h} - q_{N-1} u_{N-1}, \quad t \in (0, T]\\
    & u_n(0) = u_i(x_n), \quad n = \overline{0, N} .
\end{aligned}
\right.
```

Введём новый вектор искомой функции ``\mathbf{y} =
(u_1, u_2, \dots, u_{N-1})^T``, тогда можно записать систему в следующем виде:
```math
    \left\{
    \begin{aligned}
        &\frac{d \mathbf{y}}{d t} = \mathbf{f} \, (\mathbf{y},t), \quad t \in (0,T],\\
        &\mathbf{y}(0) = \mathbf{y}_i,
    \end{aligned}
    \right.
```
где
```math
    \begin{aligned}
        &\mathbf{y} = \big(u_1 \; u_2 \;  \ldots \; u_{N - 1} \big)^T, \\
        &\mathbf{f} = \big(f_1 \; f_2 \; \ldots \; f_{N - 1}\big)^T, \\
        &\mathbf{y}_i = \big(u_i (x_1) \; u_i (x_2) \; \ldots \; u_i (x_{N - 1}) \big)^T,
    \end{aligned}
```
вектор-функция правой части ``\mathbf{f}`` определяется следующим
образом:
```math
\mathbf{f} = \left\{
    \begin{aligned}
        &f_1 =       \varepsilon \frac{y_{2}     - 2y_1  + u_l(t)}{h^2}
        + y_1       \frac{y_{2}        - u_l(t)}{2h} - q(x_1) y_1, \\

        &f_n =       \varepsilon \frac{y_{n + 1} - 2y_n  + y_{n - 1}}{h^2}
        + y_n       \frac{y_{n + 1}    - y_{n - 1}}{2h}   - q(x_n) u_n,
        \quad n=\overline{2, N-2}, \\

        &f_{N - 1} = \varepsilon \frac{u_r(t) - 2y_{N - 1} + y_{N - 2}}{h^2}
        + y_{N - 1} \frac{u_r(t) - y_{N - 2}}{2h}   - q(x_{N-1}) y_{N-1}.
    \end{aligned}
    \right.
```

Введём равномерную сетку по времени ``T_M`` с шагом ``\tau = T/M``, состоящую
из ``M`` интервалов и ``M+1`` точек:
``T_M = \lbrace t_m, 0 \le m \le M: t_m = \tau m \rbrace ``.

Записав однастадийную схему Розенброка, ешение на следующем временном определяется как
```math
    \begin{aligned}
        &\mathbf{y}(t_{m + 1}) = \mathbf{y}(t_m) + (t_{m + 1} - t_m) \, \mathrm{Re} \, \mathbf{w},\\
    \end{aligned}
```
где ``w`` находится
```math
\begin{aligned}
    &\left[\mathbf{E} - \alpha \, (t_{m + 1} - t_m) \, \mathbf{f}_\mathbf{y}\Big(\mathbf{y}(t_m),t_m\Big)\right] \, \mathbf{w}_1 = \\
    &\qquad\qquad\qquad\qquad\quad = \mathbf{f} \, \Big(\mathbf{y}(t_m),\frac{t_{m + 1} + t_m}{2}\Big).
\end{aligned}
```

``\mathbf{f}_\mathbf{y}(\mathbf{y}(t_m), t_m)`` — якобиан функции ``f``
по вектору ``y(t_m)``, где члены явно зависящие от времени взяты при ``t\_m``.
Эта матрица Якоби имеет следующие ненулевые элементы.
```math
\begin{aligned}
    & \left(f_y\right)_{1,1}  & \equiv & \frac{\partial f_1}{\partial y_1} & = & \varepsilon\frac{-2}{h^2} - \frac{y_{2} - u_l(t)}{2h} + q(x_1), \\

     & \left(f_y\right)_{n,n - 1}  & \equiv & \frac{\partial f_n}{\partial y_{n - 1}} & = & \varepsilon \frac{1}{h^2} + \frac{y_{n}}{2h}, \quad n=\overline{2, N-1},\\

     & \left(f_y\right)_{n,n}  & \equiv & \frac{\partial f_n}{\partial y_{n}} & = &  -\varepsilon \frac{2}{h^2} - \frac{y_{n+1} - y_{n-1}}{2h} + q(x_n), \quad n=\overline{2, N-2},\\

     & \left(f_y\right)_{n,n + 1}  & \equiv & \frac{\partial f_n}{\partial y_{n + 1}} & = & \varepsilon \frac{1}{h^2} - \frac{y_{n}}{2h}, \quad n=\overline{1, N-2},\\

     & \left(f_y\right)_{N - 1,N - 1}  & \equiv & \frac{\partial f_{N - 1}}{\partial y_{N - 1}} & = &  \varepsilon \frac{-2}{h^2} - \frac{u_r(t) - y_{N - 2}}{2h} + q(x_{N-1}).
\end{aligned}
```

Если от нас потребуется использовать сеточные значения ``u_l(t), u_r(t), q(x)``,
то используются соответствующие сетки заменим их на ``u_l^m, u_r^m, q_n``
соответственно.

## Случай неравномерной сетки

Введение неравномерной сетки необходимо для эффективной обработки случаев,
когда малый параметр ``\varepsilon`` начинает приобретать достаточно большой
порядок малости, и вычисления на равномерной сетке становятся слишком долгими.
Тогда мы сгустим сетку в окрестностях особенностей решения.
В дальнейшем мы будем использовать кусочно-равномерную сетку, но формулы
запишем в общем виде для неравномерных сеток.

Подход к решению остаётся неизменным, нужно лишь модифицировать формулы
``\mathbf{f}`` и ``\mathbf{F}_{\mathbf{y}}`` для неравномерной сетки.
Запишем их с использованием сеточных значений ``u_l^m, u_r^m, q_n``.

```math
\mathbf{f} = \left\{
\begin{aligned}
        &f_1 =  \frac{2 \varepsilon}{x_2 - x_0} \left(
              \frac{y_{2} - y_1}{x_2 - x_1}
            - \frac{y_{1} + u_l^m}{x_1 - x_0}
        \right)
        + y_1 \frac{y_{2} - u_l^m}{x_2 - x_0} - q_{1} y_1, \\

        &f_n =  \frac{2 \varepsilon}{x_{n+1} - x_{n-1}}
        \left(
              \frac{y_{n+1} - y_{n}}{x_{n+1} - x_{n}}
            - \frac{y_{n} + y_{n-1}}{x_{n} - x_{n-1}}
        \right)
        + y_n  \frac{y_{n + 1} - y_{n - 1}}{x_{n+1} - x_{n-1}}
        - q_{n} u_n, \quad n=\overline{2, N-2}, \\

        &f_{N - 1} =  \frac{2 \varepsilon}{x_{N} - x_{N-1}}
        \left(
              \frac{u_r^m - y_{N - 1}}{x_{N} - x_{N-1}}
            - \frac{y_{N-1} + y_{N - 2}}{x_{N-1} - x_{N-2}}
        \right)
        + y_{N - 1} \frac{u_r^m - y_{N - 2}}{x_{N} - x_{N-2}}
        - q_{N-1} y_{N-1}.
\end{aligned}
\right.
```

Ненулевые элементы матрицы Якоби
```math
\begin{aligned}
    (f_y)_{1, 1}          & =
    \frac{2 \varepsilon}{x_2 - x_0}
    \left(
        \frac{-1}{x_{2} - x_{1}} - \frac{1}{x_{1} - x_{0}}
    \right)
    + \frac{y_{2} - u_l^m}{x_2 - x_0} - q_1, \\

    (f_y)_{n, n - 1}      & =
    \frac{2 \varepsilon}{x_{n+1} - x_{n-1}}
    \left(
        \frac{1}{x_{n} - x_{n-1}}
    \right)
    - \frac{y_{n}}{x_{n+1} - x_{n-1}}, \quad n=\overline{2, N-1},\\

    (f_y)_{n, n}          & =
    \frac{2 \varepsilon}{x_{n+1} - x_{n-1}}
    \left(
        \frac{-1}{x_{n+1} - x_{n}} - \frac{1}{x_{n} - x_{n-1}}
    \right)
    + \frac{y_{n+1} - y_{n-1}}{x_{n+1} - x_{n-1}} - q_n, \quad n=\overline{2, N-2},\\

    (f_y)_{n, n + 1}      & =
    \frac{2 \varepsilon}{x_{n+1} - x_{n-1}}
    \left(
        \frac{1}{x_{n+1} - x_{n}}
    \right)
    + \frac{y_{n}}{x_{n+1} - x_{n-1}}, \quad n=\overline{1, N-2},\\

    (f_y)_{N - 1,N - 1} & =
    \frac{2 \varepsilon}{x_{N} - x_{N-2}}
    \left(
        \frac{-1}{x_{N} - x_{N-1}} - \frac{1}{x_{N-1} - x_{N-2}}
    \right)
    + \frac{u_r^m - y_{N - 2}}{x_{N} - x_{N-2}} - q_{N-1}.
\end{aligned}
```

## Случай динамической сетки

В дальнейшем, мы построим алгоритм определения положения переходного слоя.
!!! tip
    Изменение алгоритма при динамической сетке
    Ну собственно здесь нужно нормально всё описать

## Программная реализация

*   Функция правой части
    [`NonLinearReactionAdvectionDiffusionWithFrontData.directRP`](@ref).

*   Функция якобиана
    [`NonLinearReactionAdvectionDiffusionWithFrontData.∂DRP_∂y`](@ref),
    возвращает матрицу типа `Tridiagonal` (см. [официальную
    документацию](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.Tridiagonal))

*   Функция якобиана, вычисляемого автоматическим дифференцированием
    [`NonLinearReactionAdvectionDiffusionWithFrontData.∂directRP_∂y`](@ref)
    [ForwardDiff.jl](http://www.juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.jacobian)

*   Функция поиска решения по схеме CROS1 [`solve`](@ref).

# [Сопряженная задача](@id adjoint)

## Постановка сопряженной задачи


Сопряженная задача формулируется следующим образом:
```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 \psi}{\partial x^2} +
    \frac{\partial \psi}{\partial t} =
    u \frac{\partial \psi}{\partial x} + q(x)\,\psi  -\\
    & \qquad  - 2\delta(x - f_1(t))(u(x,t) - f_2(t)),
    \quad x \in (0,1), \quad t \in (0,T], \\

    &\psi(0,t) = 0, \quad \psi(1,t) = 0, \quad t \in (0,T], \\

    &\psi(x,T) = 0, \qquad x \in [0,1].
\end{aligned}
\right.
```
Причем, важной особенностью численной реализации решения сопряженной задачи,
является динамическая, неравномерная по ``x`` сетка, уникальная для
каждого шага по времени.
Эту сетку возвращается нам функция [`solve`](@ref) решения прямой задачи.

Сопряженная задача является ретроспективной. Она должна использоваться
переданную ей сетку ``X_N^M`` в обратном по времени направлении. Найдя решение
на следующем временном слое ``\psi^{m-1}`` мы должны аппроксимировать его на
соответствующую новому временному слою сетку ``X_N^{m-1}``.

Введем следующие обозначения

|                        Описание                          |                          Обозначение                          |
|----------------------------------------------------------|---------------------------------------------------------------|
| вектор столбец искомой функции размерностью ``N-1``:     |   ``\mathbf{y} = (\psi_1, \psi_2, \dots, \psi_{N-1})^T``.     |
| вектор столбец начальных значений размерностью ``N-1``:  |   ``\mathbf{y_0} = (0, 0, \dots, 0)^T``.                      |
| вектор столбец правой части размерность ``N-1``:         |   ``\mathbf{f}(\mathbf{y}, t)``                               |
Приведём только формулы для неравномерной сетки.

Сопряженная задача легко приводится к виду
```math
    \left\{
    \begin{aligned}
        &\frac{d \mathbf{y}}{d t} = \mathbf{f} \, (\mathbf{y},t), \quad t \in [t_0,T),\\
        &\mathbf{y}(T) = \mathbf{y}_{init},
    \end{aligned}
    \right.
```

Хоть граничные условия у нас и нулевые ``\psi_r(t) = \psi_l(t) = 0``, в
формулах выпишем их явно. Причем ``f^m`` вычисляется с использованием ``X_N^m``.
```math
\mathbf{f} = \left\{
\begin{aligned}
        &f_1 =  \frac{-2 \varepsilon}{x_2 - x_0} \left(
              \frac{y_{2} - y_1}{x_2 - x_1}
            - \frac{y_{1} + \psi_l^m}{x_1 - x_0}
        \right)
        + y_1 \frac{y_{2} - \psi_l^m}{x_2 - x_0}
        + q_{1} y_1 - 2 \delta( x_1 - f_1^m)) ( u_n^m - f_2^m), \\

        &f_n =  \frac{-2 \varepsilon}{x_{n+1} - x_{n-1}}
        \left(
              \frac{y_{n+1} - y_{n}}{x_{n+1} - x_{n}}
            - \frac{y_{n} + y_{n-1}}{x_{n} - x_{n-1}}
        \right)
        + y_n  \frac{y_{n + 1} - y_{n - 1}}{x_{n+1} - x_{n-1}}
        + q_{n} u_n - 2 \delta( x_n - f_1^m)) ( u_n^m - f_2^m),
        \quad n=\overline{2, N-2}, \\

        &f_{N - 1} =  \frac{-2 \varepsilon}{x_{N} - x_{N-1}}
        \left(
              \frac{\psi_r^m - y_{N - 1}}{x_{N} - x_{N-1}}
            - \frac{y_{N-1} + y_{N - 2}}{x_{N-1} - x_{N-2}}
        \right)
        + y_{N - 1} \frac{\psi_r^m - y_{N - 2}}{x_{N} - x_{N-2}}
        + q_{N-1} y_{N-1} - 2 \delta( x_{N-1} - f_1^m)) ( u_n^m - f_2^m).
\end{aligned}
\right.
```
Ненулевые элементы Якобиана ``\mathbf{F}_{\mathbf{y}}``
```math
\begin{aligned}
    (f_y)_{1, 1}          & =
    \frac{-2 \varepsilon}{x_2 - x_0}
    \left(
        \frac{-1}{x_{2} - x_{1}} - \frac{1}{x_{1} - x_{0}}
    \right)
    + q_1, \\

    (f_y)_{n, n - 1}      & =
    \frac{-2 \varepsilon}{x_{n+1} - x_{n-1}}
    \left(
        \frac{1}{x_{n} - x_{n-1}}
    \right)
    - \frac{u_{n}^m}{x_{n+1} - x_{n-1}}, \quad n=\overline{2, N-1},\\

    (f_y)_{n, n}          & =
    \frac{-2 \varepsilon}{x_{n+1} - x_{n-1}}
    \left(
        \frac{-1}{x_{n+1} - x_{n}} - \frac{1}{x_{n} - x_{n-1}}
    \right)
    + q_n, \quad n=\overline{2, N-2},\\

    (f_y)_{n, n + 1}      & =
    \frac{-2 \varepsilon}{x_{n+1} - x_{n-1}}
    \left(
        \frac{1}{x_{n+1} - x_{n}}
    \right)
    + \frac{u_{n}^m}{x_{n+1} - x_{n-1}}, \quad n=\overline{1, N-2},\\

    (f_y)_{N - 1,N - 1} & =
    \frac{-2 \varepsilon}{x_{N} - x_{N-2}}
    \left(
        \frac{-1}{x_{N} - x_{N-1}} - \frac{1}{x_{N-1} - x_{N-2}}
    \right)
    + q_{N-1}.
\end{aligned}
```

-------------------------------------------------------------------------------

!!! note
    Правая часть в уравнении для ``\mathbf{w}`` в явном виде не зависит от
    ``t``, а все сеточные функции `u, f_1, f_2` берутся в момент времени
    ``t_m``, вместо ``\frac{ t_{m+1} + t_m}{2}``.

Найдем решение сопряженной задачи следующим образом:
```math
    \begin{aligned}
        & \mathbf{y}(t_{m - 1}) = \mathbf{y}(t_m) + (t_{m} - t_{m-1}) \,
        \mathrm{Re} \, \mathbf{w} \,\\

        & \left[\mathbf{E} - \frac{1 + i}{2} \, (t_{m} - t_{m-1}) \,
        \mathbf{F}_\mathbf{y}\Big(\mathbf{y}(t_m),t_m\Big)\right] \,
        \mathbf{w} = \\

        & \qquad\qquad\qquad\qquad\quad = \mathbf{f} \,
        \Big(\mathbf{y}(t_m), t_{m}\Big).
    \end{aligned}
```

## Программная реализация

*   Функция неоднородности
    [`NonLinearReactionAdvectionDiffusionWithFrontData.heterogeneity`](@ref).

*   Конечно-разностная аппроксимация дельта функции
    [`NonLinearReactionAdvectionDiffusionWithFrontData.deltaw`](@ref).

*   Функция правой части
    [`NonLinearReactionAdvectionDiffusionWithFrontData.adjointRP`](@ref).

*   Функция якобиана
    [`NonLinearReactionAdvectionDiffusionWithFrontData.∂ARP_∂y`](@ref), возвращает
    матрицу типа `Tridiagonal` (см. [официальную
    документацию](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.Tridiagonal))

*   Функция якобиана, вычисляемого автоматическим дифференцированием
    [`NonLinearReactionAdvectionDiffusionWithFrontData.∂adjointRP_∂y`](@ref)
    [ForwardDiff.jl](http://www.juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.jacobian)

*   Функция поиска решение по схеме CROS1 [`solve_adjoint`](@ref).

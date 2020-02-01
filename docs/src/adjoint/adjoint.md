# Сопряженная задача

Сопряженная задача формулируется следующим образом:
```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial \psi}{\partial t} = u \frac{\partial \psi}{\partial x} + q(x)\,\psi  -\\
    & \qquad  - 2\delta(x - f_1(t))(u(x,t) - f_2(t)), \quad x \in (0,1), \quad t \in (0,T), \\
    &\psi(0,t) = 0, \quad \psi(1,t) = 0, \quad t \in (0,T), \\
    &\psi(x,T) = 0, \qquad x \in [0,1].
\end{aligned}
\right.
```

Перепишем её в удобном виде для применений метода жесткий прямых
```math
\left\{
\begin{aligned}
    &\frac{\partial \psi}{\partial t} = - \varepsilon\frac{\partial^2 \psi}{\partial x^2} +  u \frac{\partial \psi}{\partial x} + q(x)\,\psi  -\\
    & \qquad  - 2\delta(x - f_1(t))(u(x,t) - f_2(t)), \quad x \in (0,1), \quad t \in (0,T), \\
    &\psi(0,t) = 0, \quad \psi(1,t) = 0, \quad t \in (0,T), \\
    &\psi(x,T) = 0, \qquad x \in [0,1].
\end{aligned}
\right.
```

Дальше будем использовать уже введенную в решении прямой задачи сетку, ``x_n \in X_n``, ``\psi_n \equiv \psi_n(t) \equiv \psi(x_n, t)``, ``u_n \equiv u_n(t) \equiv u(x_n, t)``, ``q_n \equiv q(x_n)``.
Произведем аппроксимацию частных производных ``\frac{\partial }{\partial x}``, ``\frac{\partial^2 }{\partial x^2}`` по ``x`` с помощью конечных разностей.
```math
\left\{
\begin{aligned}
    &\frac{\partial \psi_n}{\partial t} = - \varepsilon \frac{ \psi_{n+1} - 2 \psi_n + \psi_{n-1} }{h^2} +  u_n \frac{\psi_{n+1} - \psi_{n-1}}{2h} + q_n\,\psi_n  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t)), \quad n = \overline{1, N-1}, \quad t \in (0,T), \\
    &\psi(0,t) = 0, \quad \psi(1,t) = 0, \quad t \in (0,T), \\
    &\psi_n(T) = 0, \quad x \in [0, 1].
\end{aligned}
\right.
```

Сведем дифференциально-алгебраическую систему к дифференциальной, путем подстановки нулевых граничных условий.

```math
\left\{
\begin{aligned}
    &\frac{\partial \psi_1     }{\partial t} = - \varepsilon \frac{ \psi_{2}     - 2 \psi_1              }{h^2} +  u_1 \frac{\psi_{2} }{2h} + q_n\,\psi_1  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t)),  \quad t \in (0,T), \\
    &\frac{\partial \psi_n     }{\partial t} = - \varepsilon \frac{ \psi_{n+1}   - 2 \psi_n + \psi_{n-1} }{h^2} +  u_n \frac{\psi_{n+1} - \psi_{n-1}}{2h} + q_n\,\psi_n  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t)), \quad n = \overline{2, N-2}, \quad t \in (0,T), \\
    &\frac{\partial \psi_{N-1} }{\partial t} = - \varepsilon \frac{          -2 \psi_{N-1} + \psi_{N-2}  }{h^2} +  u_{N-1} \frac{ - \psi_{N-2}}{2h} + q_{N-1}\,\psi_{N-1}  -\\
    & \qquad  - 2\delta(x_{N-1} - f_1(t))(u_{N-1}(t) - f_2(t)),  \quad t \in (0,T), \\
    &\psi_n(T) = 0, \quad n = \overline{0, N}
\end{aligned}
\right.
```

Введем
 - вектор столбец искомой функции размерностью ``N-1``: ``\mathbf{y} = (\psi_1, \psi_2, \dots, \psi_{N-1})^T``.
 - вектор столбец начальных значений размерностью ``N-1``: ``\mathbf{y_i} = (0, 0, \dots, 0)^T``.
 - вектор столбец ``\mathbf{ARP}(\mathbf{y}, t)`` который будет представлять правую часть вышеприведенной системы.

!!! note
    **Замечание, о технической реализации в коде `adjointRP`**.

    `u`, `Xₙ` передаются как есть, вместе с граничными точками!
    Внутри функции они локально модифицируются, для сохранения индексации.

```math
\begin{aligned}
    &ARP_1 = - \varepsilon \frac{ y_{2} - 2 y_1 }{h^2} +  u_1 \frac{y_{2} }{2h} + q_1\,y_1  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t))\\
    &ARP_n = - \varepsilon \frac{ y_{n+1} - 2 y_n + y_{n-1} }{h^2} +  u_n \frac{y_{n+1} - y_{n-1}}{2h} + q_n\,y_n  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t)), \quad n = \overline{2, N-2}\\
    &ARP_{N-1} = - \varepsilon \frac{  2 y_{N-1} + y_{N-2} }{h^2} +  u_{N-1} \frac{ - y_{N-2}}{2h} + q_{N-1}\,y_{N-1}  -\\
    & \qquad  - 2\delta(x_{N-1} - f_1(t))(u_{N-1}(t) - f_2(t)) \\
\end{aligned}
```

!!! warning
    Сетку нужно развернуть, в тексте об этом дописать.

Используя уже введенную временную сетку по времени ``t_m \in T_m``, ``ARP_n^m \equiv ARP_n(t_m) \equiv ARP(x_n, t_m)``.
```math
\begin{aligned}
    &ARP_1^m = - \varepsilon \frac{ y_{2} - 2 y_1 }{h^2} +  u_1^m \frac{y_{2} }{2h} + q_n\,y_1  -\\
    & \qquad  - 2\delta(x_1 - f_1^m)(u_1^m - f_2^m)\\
    &ARP_n^m = - \varepsilon \frac{ y_{n+1} - 2 y_n + y_{n-1} }{h^2} +  u_n^m \frac{y_{n+1} - y_{n-1}}{2h} + q_n\,y_n  -\\
    & \qquad  - 2\delta(x_n - f_1^m)(u_n^m - f_2^m), \quad n = \overline{2, N-2}\\
    &ARP_{N-1}^m = - \varepsilon \frac{  2 y_{N-1} + y_{N-2} }{h^2} +  u_{N-1}^m \frac{ - y_{N-2}}{2h} + q_{N-1}\,y_{N-1}  -\\
    & \qquad  - 2\delta(x_{N-1} - f_1^m)(u_{N-1}^m - f_2^m) \\
\end{aligned}
```

Теперь, систему ОДУ можно записать следующим образом:
```math
    \left\{
    \begin{aligned}
        &\dfrac{d \mathbf{\textbf{y}}}{d t} = \mathbf{\textbf{ARP}} \, (\mathbf{\textbf{y}},t), \quad t \in [t_0,T),\\
        &\mathbf{\textbf{y}}(T) = \mathbf{\textbf{y}}_{init},
    \end{aligned}
    \right.
```


!!! note
    Правая часть в уравнении для ``\mathbf{W}`` в явном виде не зависит от ``t``, а все сеточные функции `u, f_1, f_2` берутся в момент времени ``t_m``, вместо ``\frac{ t_{m+1} + t_m}{2}``.

Найдем решение сопряженной задачи следующим образом:
```math
    \begin{aligned}
        &\mathbf{\textbf{y}}(t_{m + 1}) = \mathbf{\textbf{y}}(t_m) + (t_{m + 1} - t_m) \, \mathrm{Re} \, \mathbf{\textbf{W}} \,\\
        &\left[\mathbf{\textbf{E}} - \dfrac{1 + i}{2} \, (t_{m + 1} - t_m) \, \mathbf{\textbf{ARP}}_\mathbf{\textbf{y}}\Big(\mathbf{\textbf{y}}(t_m),t_m\Big)\right] \, \mathbf{\textbf{W}} = \\
        &\qquad\qquad\qquad\qquad\quad = \mathbf{\textbf{ARP}} \, \Big(\mathbf{\textbf{y}}(t_m), t_{m}\Big).
    \end{aligned}
```

![](../assets/adjoint.mp4)

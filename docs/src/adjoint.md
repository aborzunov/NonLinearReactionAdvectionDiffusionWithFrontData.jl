# Сопряженная задача

Сопряженная задача формулируется следующим образом:
```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial \psi}{\partial t} = u \frac{\partial \psi}{\partial x} + q(x)\,\psi  -\\
    & \qquad  - 2\delta(x - f_1(t))(u(x,t) - f_2(t)), \quad x \in (0,1), \quad t \in (0,T), \\
    &\psi(0,t) = 0, \quad \psi(1,t) = 0, \quad t \in (0,T), \\
    &\psi(x,t) = 0, \qquad x \in [0,1].
\end{aligned}
\right.
```

Перепишем её в удобном виде для применений метода жесткий прямых
```math
\left\{
\begin{aligned}
    &\frac{\partial \psi}{\partial t} = \varepsilon\frac{\partial^2 \psi}{\partial x^2} +  u \frac{\partial \psi}{\partial x} + q(x)\,\psi  -\\
    & \qquad  - 2\delta(x - f_1(t))(u(x,t) - f_2(t)), \quad x \in (0,1), \quad t \in (0,T), \\
    &\psi(0,t) = 0, \quad \psi(1,t) = 0, \quad t \in (0,T), \\
    &\psi(x,t) = 0, \qquad x \in [0,1].
\end{aligned}
\right.
```

Дальше будем использовать уже введенную в решении прямой задачи сетку, ``x_n \in X_n``, ``\psi_n \equiv \psi_n(t) \equiv \psi(x_n, t)``, ``u_n \equiv u_n(t) \equiv u(x_n, t)``, ``q_n \quiv q(x_n)``.
Произведем аппроксимацию частных производных ``\frac{\partial }{\partial x}``, ``\frac{\partial^2 }{\partial x^2}`` по ``x`` с помощью конечных разностей.
```math
\left\{
\begin{aligned}
    &\frac{\partial \psi_n}{\partial t} = \varepsilon \frac{ \psi_{n+1} - 2 \psi_n + \psi_{n-1} }{h^2} +  u_n \frac{\psi_{n+1} - \psi_{n-1}}{h} + q_n\,\psi_n  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t)), \quad n = \overline{1, N-1}, \quad t \in (0,T), \\
    &\psi(0,t) = 0, \quad \psi(1,t) = 0, \quad t \in (0,T), \\
    &\psi_n(t) = 0
\end{aligned}
\right.
```

Сведем дифференциально-алгебраическую систему к дифференциальной, путем подстановки нулевых граничных условий.

```math
\left\{
\begin{aligned}
    &\frac{\partial \psi_1}{\partial t} = \varepsilon \frac{ \psi_{2} - 2 \psi_1 }{h^2} +  u_n \frac{\psi_{2} }{h} + q_n\,\psi_1  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t)),  \quad t \in (0,T), \\
    &\frac{\partial \psi_n}{\partial t} = \varepsilon \frac{ \psi_{n+1} - 2 \psi_n + \psi_{n-1} }{h^2} +  u_n \frac{\psi_{n+1} - \psi_{n-1}}{h} + q_n\,\psi_n  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t)), \quad n = \overline{2, N-2}, \quad t \in (0,T), \\
    &\frac{\partial \psi_{N-1}}{\partial t} = \varepsilon \frac{  2 \psi_{N-1} + \psi_{N-2} }{h^2} +  u_n \frac{ - \psi_{N-2}}{h} + q_{N-1}\,\psi_{N-1}  -\\
    & \qquad  - 2\delta(x_{N-1} - f_1(t))(u_{N-1}(t) - f_2(t)),  \quad t \in (0,T), \\
    &\psi_n(t) = 0, \quad n = \overline{0, N}
\end{aligned}
\right.
```

Введем
 - вектор столбец искомой функции размерностью ``N-1``: ``\mathbf{y} = (u_1, u_2, \dots, u_{N-1})^T``.
 - вектор столбец начальных значений размерностью ``N-1``: ``\mathbf{y_i} = (0, 0, \dots, 0)^T``.
 - вектор столбец ``\mathbf{ARP}(y, t)`` который будет представлять правую часть вышеприведенной системы.
```math
\begin{aligned}
    &ARP_1 = \varepsilon \frac{ \psi_{2} - 2 \psi_1 }{h^2} +  u_n \frac{\psi_{2} }{h} + q_n\,\psi_1  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t))\\
    &ARP_1 = \varepsilon \frac{ \psi_{n+1} - 2 \psi_n + \psi_{n-1} }{h^2} +  u_n \frac{\psi_{n+1} - \psi_{n-1}}{h} + q_n\,\psi_n  -\\
    & \qquad  - 2\delta(x_n - f_1(t))(u_n(t) - f_2(t)), \quad n = \overline{2, N-2}\\
    &ARP_1 = \varepsilon \frac{  2 \psi_{N-1} + \psi_{N-2} }{h^2} +  u_n \frac{ - \psi_{N-2}}{h} + q_{N-1}\,\psi_{N-1}  -\\
    & \qquad  - 2\delta(x_{N-1} - f_1(t))(u_{N-1}(t) - f_2(t)) \\
\end{aligned}
```

Используя уже введенную временную сетку по времени ``t_m \in T_m``.


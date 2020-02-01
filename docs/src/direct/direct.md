# Решение прямой задачи

Если записать систему в следующем виде:
```math
    \left\{
    \begin{aligned}
        &\dfrac{d \mathbf{\textbf{y}}}{d t} = \mathbf{\textbf{f}} \, (\mathbf{\textbf{y}},t), \quad t \in (t_0,T],\\
        &\mathbf{\textbf{y}}(t_0) = \mathbf{\textbf{y}}_{init},
    \end{aligned}
    \right.
```
где
```math
    \begin{aligned}
        &\mathbf{\textbf{y}} = \big(u_1 \; u_2 \;  \ldots \; u_{N - 1} \big)^T, \\
        &\mathbf{\textbf{f}} = \big(f_1 \; f_2 \; \ldots \; f_{N - 1}\big)^T, \\
        &\mathbf{\textbf{y}}_{init} = \big(u_{init} (x_1) \; u_{init} (x_2) \; \ldots \; u_{init} (x_{N - 1}) \big)^T.
    \end{aligned}
```
``u_init(x_n)`` вычисляется с помощью [`u_init(x)`](@ref).


То текущая функция определяет вектор-функцию $\mathbf{\textbf{f}}$ следующим образом:
```math
    \begin{aligned}
        &f_1 =       \varepsilon \dfrac{y_{2}        - 2y_1       + u_{left}(t)}{h^2} + y_1       \dfrac{y_{2}        - u_{left}(t)}{2h} - q(x_1) y_1, \\
        &f_n =       \varepsilon \dfrac{y_{n + 1}    - 2y_n       + y_{n - 1}}{h^2}   + y_n       \dfrac{y_{n + 1}    - y_{n - 1}}{2h}   - q(x_n) u_n, \quad n=\overline{2, N-2}, \\
        &f_{N - 1} = \varepsilon \dfrac{u_{right}(t) - 2y_{N - 1} + y_{N - 2}}{h^2}   + y_{N - 1} \dfrac{u_{right}(t) - y_{N - 2}}{2h}   - q(x_{N-1}) y_{N-1}.
    \end{aligned}
```
На каждом временном шаге, решение находится как:
```math
    \begin{aligned}
        &\mathbf{\textbf{y}}(t_{m + 1}) = \mathbf{\textbf{y}}(t_m) + (t_{m + 1} - t_m) \, \mathrm{Re} \, \mathbf{\textbf{w}}_1,\\
    \end{aligned}
```
где ``W_1`` находится из
```math
\begin{aligned}
    &\left[\mathbf{\textbf{E}} - \alpha \, (t_{m + 1} - t_m) \, \mathbf{\textbf{f}}_\mathbf{\textbf{y}}\Big(\mathbf{\textbf{y}}(t_m),t_m\Big)\right] \, \mathbf{\textbf{w}}_1 = \\
    &\qquad\qquad\qquad\qquad\quad = \mathbf{\textbf{f}} \, \Big(\mathbf{\textbf{y}}(t_m),\frac{t_{m + 1} + t_m}{2}\Big).
\end{aligned}
```

``\mathbf{f}_\mathbf{y}(\mathbf{y}(t_m), t_m)`` — якобиан функции [`directRP(...)`](@ref) по вектору ``y`` (в момент времени ``t_m``) в момент времени ``t\_m``.
Эта матрица Якоби имеет следущие ненулевые элементы.
```math
\begin{aligned}
    & \left(f_y\right)_{1,1}  & \equiv & \frac{\partial f_1}{\partial y_1} & = & \varepsilon\dfrac{-2}{h^2} - \dfrac{y_{2} - u_{left}(t)}{2h} + q(x_1), \\

     & \left(f_y\right)_{n,n - 1}  & \equiv & \frac{\partial f_n}{\partial y_{n - 1}} & = & \varepsilon \dfrac{1}{h^2} + \dfrac{y_{n}}{2h}, \quad n=\overline{2, N-1},\\

     & \left(f_y\right)_{n,n}  & \equiv & \frac{\partial f_n}{\partial y_{n}} & = &  -\varepsilon \dfrac{2}{h^2} - \dfrac{y_{n+1} - y_{n-1}}{2h} + q(x_n), \quad n=\overline{2, N-2},\\

     & \left(f_y\right)_{n,n + 1}  & \equiv & \frac{\partial f_n}{\partial y_{n + 1}} & = & \varepsilon \dfrac{1}{h^2} - \dfrac{y_{n}}{2h}, \quad n=\overline{1, N-2},\\

     & \left(f_y\right)_{N - 1,N - 1}  & \equiv & \frac{\partial f_{N - 1}}{\partial y_{N - 1}} & = &  \varepsilon \dfrac{-2}{h^2} - \dfrac{u_{right}(t) - y_{N - 2}}{2h} + q(x_N).
\end{aligned}
```

![](../assets/solution.mp4)

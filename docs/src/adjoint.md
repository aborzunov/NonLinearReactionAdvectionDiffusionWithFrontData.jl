# Сопряженная задача

Сопряженная задача формулируется следующим образом:
```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 \psi}}{\partial x^2} + \frac{\partial \psi}}{\partial t} = u}\frac{\partial \psi}}{\partial x} +  q}(x)\,\psi} - \\
    &\hspace{1.2cm} -2\delta(x-f_1(t))(u}(x,t)-f_2(t))=0, \quad x \in (0,1), \quad t \in [0,T), \\
    &\psi}(0,t) = 0, \quad \psi}(1,t) = 0, \quad t \in [0,T), \\
    &\psi}(x,T) = 0, \qquad x \in (0,1).
\end{aligned}
\right.
```

Перепишем её в удобном виде для применений метода жесткий прямых
```math
\left\{
\begin{aligned}
    \frac{\partial \psi}}{\partial t} = &\varepsilon\frac{\partial^2 \psi}}{\partial x^2} + u}\frac{\partial \psi}}{\partial x} +  q}(x)\,\psi} - \\
    &\hspace{1.2cm} -2\delta(x-f_1(t))(u}(x,t)-f_2(t)), \quad x \in (0,1), \quad t \in [0,T), \\
    &\psi}(0,t) = 0, \quad \psi}(1,t) = 0, \quad t \in [0,T), \\
    &\psi}(x,T) = 0, \qquad x \in (0,1).
\end{aligned}
\right.
`

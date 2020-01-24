```@meta
DocTestSetup = quote
    using NonLinearReactionAdvectionDiffusionWithFrontData
end

```# Главная

```@contents
```

# Постановка задачи

```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 u}{\partial x^2} - \frac{\partial u}{\partial t} = -u \frac{\partial u}{\partial x} +  q(x)\,u, \quad x \in (0,1), \quad t \in (0,T), \\
    &u(0,t) = u_{left}(t), \quad u(1,t) = u_{right}(t), \quad t \in (0,T), \\
    &u(x,t) = u_{init}(x), \qquad x \in [0,1].
\end{aligned}
\right.
```

# Пример решения прямой задачи

Зададим параметры.
```@example demo
using NonLinearReactionAdvectionDiffusionWithFrontData
using Plots

u_l(t) = -8
u_r(t) = 4
ε = 0.2;
a, b = 0, 1;
t₀, T = 0, 1;
N, M = 40, 80;
h = (b-a)/N;
τ = (T-t₀)/M;
Xₙ = [a  + n*h for n in 0:N];
Tₘ = [t₀ + m*τ for m in 0:M];
q(x) = sin(3 * π * x);
u = zeros(M+1, N+1);
nothing #hide
```

Зададим начальное приближение [`u_init`](@ref)
```@example demo
y = u_init.( Xₙ[n] for n in 2:N );
qₙ = [ q(x) for x in Xₙ[2:N] ];
```

Решим задачу [`solve!`](@ref) и построим gif решения [`make_gif`](@ref)
```@example demo
u= solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, qₙ)
make_gif(u, Xₙ, Tₘ; frame_skip = div(M,50), frames_to_write=80, name="solution.gif");
nothing # hide
```

![](solution.gif)

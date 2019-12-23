# Главная

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

# Пример решения

Зададим параметры.
```jldoctest demo
julia> using NonLinearReactionAdvectionDiffusionWithFrontData

julia> u_l(t) = -8#*cos(2*π * t);


julia> u_r(t) = 4# * (1 + sin(2*π * t));


julia> ε = 0.2;


julia> a, b = 0, 1;


julia> t₀, T = 0, 1;


julia> N, M = 40, 80;


julia> h = (b-a)/N;


julia> τ = (T-t₀)/M;


julia> Xₙ = [a  + n*h for n in 0:N];


julia> Tₘ = [t₀ + m*τ for m in 0:M];


julia> q(x) = sin(3 * π * x)

q (generic function with 1 method)

julia> u = zeros(M+1, N+1);
```

Зададим начальное приближение
```jldoctest demo
julia> y = u_init.( Xₙ[n] for n in 2:N )
39-element Array{Float64,1}:
 -8.010340876819335
 -8.017828522120382
 -8.006733491677299
 -7.958156688432881
 -7.833646561079693
 -7.5583895218691985
 -7.000181421210685
 -5.970893714323723
 -4.324519390104715
 -2.1875
 -0.049230609895282385
  1.6008937143237265
  2.636431421210687
  3.203389521869199
  3.489896561079693
  3.628156688432881
  3.6929834916772997
  3.7228285221203827
  3.7365908768193337
  3.7433666563569172
  3.747490716170851
  3.7510192650881646
  3.754925504812066
  3.759669571706625
  3.7654689144584026
  3.7724262699047735
  3.780590172256178
  3.7899835485135176
  3.8006172288624214
  3.8124963291732765
  3.8256232660239555
  3.8399991809276504
  3.8556246130976035
  3.872499817240246
  3.890624913670404
  3.9099999592207864
  3.930624980737263
  3.9524999909009275
  3.975624995701903
```
Решим задачу
```jldoctest demo
julia> u = solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, q; α = complex(0.5, 0.5))

41×81 Array{Float64,2}:
 [...]
```

Построим gif решения
```jldoctest demo
julia> make_gif(u, Xₙ, Tₘ; frame_skip = div(M,30), frames_to_write=80, name="build/solution.gif"); nothing # hide
```
![](solution.gif)

# Решение прямой задачи

```@example demo
using NonLinearReactionAdvectionDiffusionWithFrontData
using Plots

# Зададим параметры для прямой задачи
u_l(t) = -8#*cos(2*π * t);
u_r(t) = 4# * (1 + sin(2*π * t));
ε = 0.2;
a, b = 0, 1; # Область по Х
t₀, T = 0, 0.28;
N, M = 40, 80;
h = (b-a)/N;
τ = (T-t₀)/M;
Xₙ = [a  + n*h for n in 0:N];
Tₘ = [t₀ + m*τ for m in 0:M];
q(x) = 4*sin(3 * π * x);
u = zeros(M+1, N+1);

# Начальные условия
y = u_init.( Xₙ[n] for n in 2:N );
# Некоторая функция q, т.к. мы моделируем априорные данные f1, f2,
# то она нам известна
qₙ = [ q(x) for x in Xₙ[2:N] ];
u = solve!(y, Xₙ, Tₘ, N, M, ε, u_l, u_r, qₙ);
nothing #hide
```

# Генерирование априорной информации

```@example demo
# Вырожденные корни
qn = [q(x) for x in Xₙ]
q_low = NonLinearReactionAdvectionDiffusionWithFrontData.phidetermination(qn, y, u_l(0), Xₙ, N::Int);
q_top = NonLinearReactionAdvectionDiffusionWithFrontData.phidetermination(qn, y, u_r(0), Xₙ[end:-1:1], N::Int);
q_top = q_top[end:-1:1];
# Полуразность вырожденных корней
ϕ = NonLinearReactionAdvectionDiffusionWithFrontData.Φ(q_low, q_top, N);
# Положение переходного слоя
f1 = NonLinearReactionAdvectionDiffusionWithFrontData.f1(ϕ, u, Xₙ, N, M);
# Значение функции на переходном слое
f2 = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, u, Xₙ, N, M);
nothing #hide
```

# График и анимация решения

```@xample demo
make_plot(u, Xₙ, Tₘ, 50, q_low, q_top, f1, f2)
```

```@example demo
make_gif(u, Xₙ, Tₘ, q_low, q_top, f1, f2; frame_skip = div(M,50), frames_to_write=80, convert2mp4 = true)
nothing #hide
```

![](./assets/solution.mp4)
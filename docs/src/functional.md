# Вычисление функционала

```math
J[q] = \int\limits^T_0 \left( u(f_1(t), t; q) - f_2(t) \right)^2 dt + \alpha \int\limits_0^1 \left( q(x) \right)^2 dx,
```
где ``\alpha`` — параметр регуляризации, пока равный нулю.

# Вычисление градиента функционала

```math
J'[q^s](x) = \int\limits_0^T u^{(s)} (x,t) \psi^{(s)}(x,t) dt
```

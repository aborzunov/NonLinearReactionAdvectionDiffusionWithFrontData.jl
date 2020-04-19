# ## [Пример №4, краткий, динамическая сетка](@id de_4)
#
# Делаем тоже самое, что в [развернотом примере](@ref de_3),
# только с использованием функций сокращений и стандартных параметров.
using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: shishkin_mesh;

## Задача параметров и решение в коротком виде
a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀, mshfrm_std = NonLinearReactionAdvectionDiffusionWithFrontData.dparams_nonuniform();
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, create_mesh = mshfrm_std);
nothing #hide


## Генерация априорной информации
ϕl      = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                  # Левый вырожденный корень
ϕr      = phidetermination(qₙ, urₘ, Xₙ, N, Tₘ, M, reverseX = true); # Правый вырожденный корень
ϕ       = Φ(ϕl, ϕr, N, M);                                          # Полуразность вырожденных корней
ϕ       = apply_on_dynamic_mesh(ϕ, XX, N, M);                       # Аппроксимация на переменную сетку
ϕl      = apply_on_dynamic_mesh(ϕl, XX, N, M);                      # Аппроксимация на переменную сетку
ϕr      = apply_on_dynamic_mesh(ϕr, XX, N, M);                      # Аппроксимация на переменную сетку
f1_data = f1(ϕ, u, XX, N, M);                                       # Положение переходного слоя
f2_data = f2(f1_data, u, XX, N, M);                                 # Значение функции на переходном слое
nothing # hide

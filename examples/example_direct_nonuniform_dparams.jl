using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: shishkin_mesh;
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using NonLinearReactionAdvectionDiffusionWithFrontData: phidetermination, apply_on_dynamic_mesh, Φ;

# ## Задача параметров и решение в коротком виде
#
a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀, mshfrm = NonLinearReactionAdvectionDiffusionWithFrontData.dparams_nonuniform();
u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, create_mesh = mshfrm);
nothing #hide


# ## Генерация априорной информации
ϕl      = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                  # Левый вырожденный корень
ϕr      = phidetermination(qₙ, urₘ, Xₙ, N, Tₘ, M, reverseX = true); # Правый вырожденный корень
ϕ       = Φ(ϕl, ϕr, N, M);                                          # Полуразность вырожденных корней
ϕ       = apply_on_dynamic_mesh(ϕ, XX, N, M);                       # Аппроксимация на переменную сетку
ϕl      = apply_on_dynamic_mesh(ϕl, XX, N, M);                      # Аппроксимация на переменную сетку
ϕr      = apply_on_dynamic_mesh(ϕr, XX, N, M);                      # Аппроксимация на переменную сетку
f1_data = f1(ϕ, u, XX, N, M);                                       # Положение переходного слоя
f2_data = f2(f1_data, u, XX, N, M);                                 # Значение функции на переходном слое
nothing # hide

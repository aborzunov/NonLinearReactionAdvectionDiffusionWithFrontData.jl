using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: shishkin_mesh;
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using NonLinearReactionAdvectionDiffusionWithFrontData: phidetermination, apply_on_dynamic_mesh, Φ;

u_l(t) = -8;                    # ГУ
u_r(t) =  4;                    #
qf(x) = 4*sin(3 * π * x);       # Коэффициент линейного усиления, который в обратной
                                # задаче необходимо определить, но при генерации априорной
                                # информации мы задаем некоторый коэффициент, который,
                                # собственно, после имея априорную информацию и будем определять.
ε = 0.01;                       # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 0.50;                # Область по T
x_tp = 0.12;                    # Положение переходного слоя
M = 500;                        # Кол-во разбиений по X, T
τ = (T-t₀)/M;                   # шаг по T
Tₘ = [t₀ + m*τ for m in 0:M];   # Сетка по Т
ulₘ=    u_l.(Tₘ);               # Сеточные значения левого  ГУ
urₘ=    u_r.(Tₘ);               # Сеточные значения правого ГУ
nothing #hide

#' Новое поведение будет реализоваться с помощью передачи функции создания сетки внутрь `solve` в качестве аргумента.
#' Эта функция должна принимать только один аргумент — положение переходного слоя и формировать соответствующую сетку.
#' В пакете есть формирование кусочно-равномерной сетки со сгущениями на границе и на переходном слое.
#' Создадим замыкание этой функции, которое будет иметь нужную сигнатуру.
mshfrm(x_tp) = shishkin_mesh(a, b, x_tp, ε, 40, 0.5, 1.0, 1.0, 0.25);
Xₙ = mshfrm(x_tp); ;                            # Сетка по Х
N = length(Xₙ) - 1                              # Примем за N длину сетки, что получилась в итоге.
qₙ =    qf.(Xₙ);                                # Сеточные значения коэффициента линейного усиления
u₀ = u_init.(Xₙ, ε=ε, x_tp = x_tp);             # Начальные условия

u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, create_mesh = mshfrm);
nothing #hide


#' ## Генерация априорной информации
ϕl = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                               # Левый вырожденный корень
ϕr = phidetermination(reverse(qₙ), urₘ, reverse(Xₙ), N, Tₘ, M);             # Нужно подать инвертированную сетку
ϕr = reverse(ϕr, dims=1);                                                   # А после — инвертировать решение по X
ϕ = Φ(ϕl, ϕr, N, M);                                                        # Полуразность вырожденных корней
ϕ = apply_on_dynamic_mesh(ϕ, XX, N, M);                                     # Аппроксимация на переменную сетку
ϕl = apply_on_dynamic_mesh(ϕl, XX, N, M);                                   # Аппроксимация на переменную сетку
ϕr = apply_on_dynamic_mesh(ϕr, XX, N, M);                                   # Аппроксимация на переменную сетку
f1_data = f1(ϕ, u, XX, N, M);                                               # Положение переходного слоя
f2_data = f2(f1_data, u, XX, N, M);                                         # Значение функции на переходном слое
nothing # hide

#' ## Визуализация
#' Запись **только** mp4 вместе с вырожденными корнями
make_gif(u, XX, Tₘ, ϕl, ϕr, ϕ, f1_data, f2_data; convert2mp4 = true, name="example_direct_nonuniform_with_f1_f2.gif")

#' ## Решение прямой задачи для генерации априорной информации
#' Сопряженная задача определяется для сеточного решения некоторой прямой задачи,
#' и некоторых априорных данных. Поэтому, сначала необходимо сгенерировать
#' априорные данные и некоторую матрицу `u`.
#' Сделаем это на увеличенном числе интервалов, т.к. решение сопряженной задачи
#' менее гладкое.
using NonLinearReactionAdvectionDiffusionWithFrontData
using NonLinearReactionAdvectionDiffusionWithFrontData: f1, f2;
using NonLinearReactionAdvectionDiffusionWithFrontData: shishkin_mesh;
using NonLinearReactionAdvectionDiffusionWithFrontData: phidetermination, Φ;
using NonLinearReactionAdvectionDiffusionWithFrontData: apply_on_dynamic_mesh;

u_l(t) = -8;                    # ГУ
u_r(t) =  4;                    #
qf(x) = 4*sin(3 * π * x);       # Коэффициент линейного усиления, который в обратной
                                # задаче необходимо определить, но при генерации априорной
                                # информации мы задаем некоторый коэффициент, который,
                                # собственно, после имея априорную информацию и будем определять.
ε = 0.03;                       # Малый параметр при старшей производной
a, b = 0, 1;                    # Область по X
t₀, T = 0, 0.20;                # Область по T
x_tp = 0.22;                    # Положение переходного слоя
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
mshfrm(x_tp) = shishkin_mesh(a, b, x_tp, ε, 40, 1.0, 1.0, 0.2, 1.0);
Xₙ = mshfrm(x_tp); ;                            # Сетка по Х
N = length(Xₙ) - 1                              # Примем за N длину сетки, что получилась в итоге.
qₙ =  1.0 * qf.(Xₙ);                            # Сеточные значения коэффициента линейного усиления
                                                # Измените линейный коэффициент, чтобы посмотреть
                                                # на поведение сопряженной задачи (а следовательно градиента),
                                                # когда прямая задача не соответствует наблюдаемым данным.
u₀ = u_init.(Xₙ, ε=ε, x_tp = x_tp);             # Начальные условия

u, XX, TP = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ, create_mesh = mshfrm);
#' ## Генерация априорной информации
ϕl      = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                          # Левый вырожденный корень
ϕr      = phidetermination(reverse(qₙ), urₘ, reverse(Xₙ), N, Tₘ, M);        # Нужно подать инвертированную сетку и q
ϕr      = reverse(ϕr, dims=1);                                              # А после — инвертировать решение по X
ϕ       = Φ(ϕl, ϕr, N, M);                                                  # Полуразность вырожденных корней
ϕ       = apply_on_dynamic_mesh(ϕ, XX, N, M);                               # Аппроксимация на переменную сетку
ϕl      = apply_on_dynamic_mesh(ϕl, XX, N, M);                              # Аппроксимация на переменную сетку
ϕr      = apply_on_dynamic_mesh(ϕr, XX, N, M);                              # Аппроксимация на переменную сетку
f1_data = f1(ϕ, u, XX, N, M);                                               # Положение переходного слоя
f2_data = f2(f1_data, u, XX, N, M);                                         # Значение функции на переходном слое
#
w = 0.000051;                                                               # Априорный параметр апроксимации дельта-функции

#' Проведем контроль выбора априорного параметра в аппроксимации дельта-функции
#' На точных данных, эта неоднородность должна быть нулевой почти везде
#' Ведь сопряженная задача — ретроспективная. Её начальные условие — нулевые.
#' Нулевое решение сопряженной задачи даст нам нулевой градиент,
#' а значит мы нашли решение.
using  NonLinearReactionAdvectionDiffusionWithFrontData: heterogeneity;
using LaTeXStrings, Plots;
Uₙₘ = u[2:N, :];
X = XX[2:N, :];
H = [  - heterogeneity(n, m, X[:, m], N, Tₘ, M, Uₙₘ, f1_data, f2_data, w) for n in 1:N-1, m in 1:M+1]
heatmap(H', title=L"-2 δ( x - f_1(t))(u(f1(t), t) - f2(t) ")

#' ## Решение сопряженной задачи

Uₙₘ = u;                        # Сохраним старую матрицу
ψ₀ = [0.0 for i in 1:N+1];      # Нулевые начальные условия
ψl = [0.0 for i in 1:M+1];      # Нулевые ГУ
ψr = [0.0 for i in 1:M+1];      # Нулевые ГУ
ψ = solve_adjoint(ψ₀, XX, N, Tₘ, M, ε, ψl, ψr, qₙ, Uₙₘ, f1_data, f2_data, w = w)
nothing #hide

#' ## Визуализация

# На отрисовку, решение сопряженной задачи передадим в инвертированном времени.
# Передадим свою подпись к графиками с помощью keyword `label="\\psi"`,
# не забыть про экранировку спец символа
# Нарисует гиф с неравномерной скоростью по оси T, первые 50 шагов
# отрисуем полностью, а из последующих выберем каждый десятый и последние — опять полностью.
# С помощью keyword `frames_to_write`
ftw = [1:50; 51:10:M-50; M-49:M+1];
make_gif(ψ[:,end:-1:1], XX, Tₘ; frames_to_write=ftw, label="\\psi",
         name="adjoint_example_nonuniform.gif", convert2mp4 = true)

#' Результат должен быть около нулевой, ведь в качестве текущего приближения `q` мы взяли искомое,
#' а при нем — градиент должен обнуляться.

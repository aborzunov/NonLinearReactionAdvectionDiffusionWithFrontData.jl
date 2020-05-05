# Вычисление вырожденных корней

@doc raw"""
    generate_obs_data(u::Matrix, Xₙ::Vector, N::Int,
                      Tₘ::Vector, M::Int,
                      qₙ::Vector,
                      ulₘ::Vector, urₘ::Vector)

Фнукция-сокращение.

#Return
Левый вырожденный корень, правый, их полуразность, положение переходного слоя,
значение `u` на переходном слое.
`ϕl, ϕr, ϕ, f1_data, f2_data`

See also: [`phidetermination`](@ref), [`Φ`](@ref), [`f1`](@ref), [`f2`](@ref).
"""
function generate_obs_data(u::Matrix, Xₙ::Vector, N::Int,
                           Tₘ::Vector, M::Int,
                           qₙ::Vector,
                           ulₘ::Vector, urₘ::Vector)

    ϕl      = phidetermination(qₙ, ulₘ, Xₙ, N, Tₘ, M);                      # Левый вырожденный корень
    ϕr      = phidetermination(qₙ, urₘ, Xₙ, N, Tₘ, M, reverseX = true);     # Правый вырожденный корень
    ϕ       = Φ(ϕl, ϕr, N, M);                                              # Серединный корень
    f1_data = f1(ϕ, u, Xₙ, N, M);                                           # Положение переходного слоя
    f2_data = f2(f1_data, u, Xₙ, N, M);                                     # Значение функции на переходном слое

    return ϕl, ϕr, ϕ, f1_data, f2_data
end

@doc raw"""
    phidetermination(q::Vector, ub::Vector,
                     Xₙ::Vector, N::Int,
                     Tₘ::Vector, M::Int;
                     reverseX = false)

Решает ОДУ для нахождения вырожденного корня.

`reverseX` флаг обозначающий обратное направление интегрирования по оси ``X``.
"""
function phidetermination(qₙ::Vector, ub::Vector,
                          Xₙ::Vector, N::Int,
                          Tₘ::Vector, M::Int;
                          reverseX::Bool = false)

    @assert length(qₙ) == N + 1
    @assert length(Xₙ) == N + 1

    @assert length(ub) == M + 1
    @assert length(Tₘ) == M + 1

    if (Xₙ[end] - Xₙ[1]) < 0
        throw(ArgumentError("Сетку следует передавать в нормальном виде," *
                            " а не инвертированном. Исльзуйте kw `reverseX`."))
    end

    X = reverseX ? reverse(Xₙ) : Xₙ;
    q = reverseX ? reverse(qₙ) : qₙ;

    phi = zeros(N + 1, M + 1);

    for m in 1:M+1

        # Граничное условие на ϕ
        phi[1,m] = ub[m];

        for n = 1:N
            phi[n+1, m ] = phi[n, m ] + (X[n+1] - X[n]) * ( q[n+1] + q[n] ) / 2;
        end
    end

    return reverseX ? reverse(phi, dims=1) : phi
end

@doc raw"""
    Φ(ϕ_l::Matrix, ϕ_r::Matrix, N::Int, M::Int) -> ::Vector

Вычисляет полуразность вырожденных решений
``|\phi_l^m - \phi_r^m|/2 + \phi_l^m`` на каждом шаге по времени `m`
с помощью матриц вырожденных решений ``\phi_l`` и ``\phi_r`` вырожденного корня.

"""
function Φ(ϕ_l::Matrix, ϕ_r::Matrix, N::Int, M::Int)
    @assert size(ϕ_l) == (N+1, M+1)
    @assert size(ϕ_r) == (N+1, M+1)

    Φ = zeros(size(ϕ_l))
    Φ = [ abs(ϕ_l[n, m] - ϕ_r[n, m])/2 + ϕ_l[n, m] for n in 1:N+1, m in 1:M+1]
    return Φ;
end

@doc raw"""
    apply_on_dynamic_mesh(ϕ::Matrix, XX::Matrix,
                      N::Int, M::Int) -> Matrix

Аппроксимирует `ϕ` на динамичски изменяющуюся на каждом временном шаге сетку
`XX`.

Функция [`phidetermination`](@ref) определяют вырожденные корни на стартовой
сетке.  [`Φ`](@ref) проиозводит линейные операции над векторами на каждом
временном шаге, опять же на стартовой сетке. Поэтому если в расчетах была
использована динамическая сетка, то нужно переопределить `ϕ` на этой
динамической сетке.
"""
function apply_on_dynamic_mesh(ϕ::Matrix, XX::Matrix,
                           N::Int, M::Int)
    out = zero(ϕ)
    size(ϕ) == size(XX) || throw(ArgumentError("Входные матрицы различной длины"))

    size(ϕ) == size(XX) == (N+1, M+1) || throw(ArgumentError("Входные матрицы неверной длины" *
                        "size(ϕ) = $(size(ϕ)), size(XX) = $(size(XX)), (N+1, M+1) = $((N+1, M+1))"))

    # Функция `phidetermination` работает на статичной сетке.
    # Функция `Φ` просто производит линейные операции.
    # Поэтому `ϕ`, наш входной аргумент определен на статичной сетке на первом временном шаге.
    for m in 1:M+1
        # Так, мы будем брать `ϕ` на шаге m, а сетку — из самого первого.
        # После, интерполировать на текущую сетку.
        spl = Spline1D( XX[:, 1], ϕ[:, m ] )
        out[:, m] = spl(XX[:, m])
    end

    return out
end

@doc raw"""
    find_f_zeros(f::Vector, Xₙ::Vector)

Находит такой ``x``, что ``f(x) = 0``.
`f` — сеточные значения функции на сетке `Xₙ`.

!!! warning
     - Функция обязана пересекать ноль. Не сработает на неотрицательных функциях.
     - Возвращает только аргуент реализующий **первый** ноль.
     - Решение ищется аппроксимацией.
"""
function find_f_zeros(f::Vector, Xₙ::Vector)

    N = length(f);
    @assert length(f) == length(Xₙ) "Массивы f, Xₙ разной длины: $((length(f), length(Xₙ)))"

    k = 1;
    # Найдем номер элемента, после которого стоит положительный элемент
    # Если вначале функция отрицательная, то произведение больше нуля
    # Если вначале функция положительная, то произведение больше нуля
    for n in 1:N-1
        k = n;
        if f[n]*f[n+1] < 0
            break;
        end
    end

    if k == length(f) - 1
        throw(ArgumentError("Функция не пересекает нуля"))
    end

    # Нормальный случай, когда мы можем выбрать по 4 элемента с каждой стороны от k-ой точки.
    if 4 <= k <= length(f) - 4
        # Выберем 8 точек, по 4 с каждой стороны от нуля функции.
        x = Xₙ[k-3:k+4]
        y = f[k-3:k+4]
    elseif k < 4
        # Если ноль лежит вблизи левой границы, возьмем первых 8 точек
        x = Xₙ[1:8]
        y = f[1:8]
    else
        # Если ноль лежит вблизи правой границы, возьмем последний 8 точек
        x = Xₙ[end - 7:end]
        y = f[end - 7:end]

    end

    # Проверим функцию на монотонность.
    # Spline1D требует упорядоченной возрастающей по ``x`` сетки.
    # В лучшем случае, `y` — упорядоченный массив по возрастанию.
    if !issorted(y)

        # Если он упорядочен по убыванию, то развернем его.
        y = reverse(y);
        x = reverse(x);

        # Проверим, что он упорядочен по возрастанию.
        if !issorted(y) # Вектор не упорядочен ни в одном из двух направлений, f(x) немонотонна
            errstr =    "Окно из 8 значений функции ``f`` немонотонно! Нужно увеличить плотность сетки." *
                        "\nx = $(x)" *
                        "\ny = $(y)";
            throw(ArgumentError(errstr))
        end
    end

    # Теперь построим функцию ``x = f^-1(y)``
    # С помощью интерполяции найдем ``f^-1(0)``
    # Внимание!
    # Интерполяция производится на неравномерной сетке: см массив `y`
    spl = Spline1D(y, x)    # это аппроксимация ``f^-1(y) = x``.

    return spl(0);
end

@doc raw"""
    f1(ϕ::Matrix, u::Matrix, Xₙ::Array, N::Int, M::Int) -> Vector

Находит положение переходного слоя ``f_1(t) = x_{t.p.}(t)``
путем поиска точки пересечения аргумента, при котором пересекаются ``u(x,t)`` и ``ϕ(x)``.

Точка пересечения находится путем интерполяции функции обратной к ``u(x,t) - ϕ(x) = 0``,
передавай аргументы в инвертированном виде в [`find_f_zeros`](@ref).

See also: [`find_f_zeros`](@ref)
"""
function f1(ϕ::Matrix, u::Matrix, Xₙ::Array, N::Int, M::Int)
    @assert size(u) == (N+1, M+1)
    @assert size(ϕ) == (N+1, M+1)

    f1 = zeros(M+1);

    # Необходимо найти абсциссу пересечения Φ и u(x, t), т.е. аргумент х на каждом временном шаге
    if typeof(Xₙ) <: Vector
        for m in 1:M+1
            f1[m] = find_f_zeros(u[:,m] - ϕ[:,m], Xₙ);
        end
    elseif typeof(Xₙ) <: Matrix
        for m in 1:M+1
            # Только в отличии от прошлого случая, здесь `ϕ, u` — определены на динамической сетке
            # Мы должны подавать соответствующую сетку на каждом шаге
            f1[m] = find_f_zeros(u[:,m] - ϕ[:,m], Xₙ[:,m]);
        end
    else
        throw(ArgumentError("Функция принимает `Xₙ` только в виде вектора или матрицы"))
    end

    return f1
end

@doc raw"""
    f2(ϕ::Matrix, u::Matrix, Xₙ::Array, N::Int, M::Int) -> Vector

Находит значение искомой функции на переходном слое ``f_2(t) = u(x_{t.p.}, t)``.
Находится путем интерполяции функции ``u(x - f1(t), t) = 0``.

See also: [`find_f_zeros`](@ref)
"""
function f2(f1::Vector, u::Matrix, Xₙ::Array, N::Int, M::Int)
    @assert size(u) == (N+1, M+1)
    @assert length(f1) == M+1;

    f2 = zeros(M+1);

    # находим значение функции на каждом временном шаге
    if typeof(Xₙ) <: Vector
        for m in 1:M+1
            f2[m] = find_f_zeros(Xₙ .- f1[m], u[:, m]); # Передадим аргументы в обратном порядке
        end
    elseif typeof(Xₙ) <: Matrix
        for m in 1:M+1
            # Только в отличии от прошлого случая, здесь `ϕ, u` — определены на динамической сетке
            # Мы должны подавать соответствующую сетку на каждом шаге
            f2[m] = find_f_zeros(Xₙ[:, m] .- f1[m], u[:, m]); # Передадим аргументы в обратном порядке
        end
    else
        throw(ArgumentError("Функция принимает `Xₙ` только в виде вектора или матрицы"))
    end

    return f2
end

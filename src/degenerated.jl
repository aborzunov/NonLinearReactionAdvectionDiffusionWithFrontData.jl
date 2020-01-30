# Вычисление вырожденных корней

@doc raw"""
    phidetermination(q, u, u_border, Xₙ, N::Int)

Решает ОДУ для нахождения вырожденных корней.

!!! danger
    Функция некорректно обрабатыет непостоянные ГУ.
"""
function phidetermination(q::Vector, u::Vector,
                          u_border::Vector, Xₙ::Vector, N::Int)
    @assert length(Xₙ) == N + 1

    phi = similar(Xₙ);
    phi[1] = u_border[1];

    for n = 1:N
        phi[n+1] = phi[n] + (Xₙ[n+1] - Xₙ[n]) * ( q[n+1] + q[n] ) / 2;
    end

    return phi
end

@doc raw"""
    Φ(ϕ_l::Vector, ϕ_r::Vector, N::Int) -> ::Vector

Вычисляет значение функции на переходном слое ``|\phi_l - \phi_r|/2`` с помощью левого ``\phi_l`` и правого ``\phi_r`` вырожденного корня.
"""
function Φ(ϕ_l::Vector, ϕ_r::Vector, N::Int)
    Φ = similar(ϕ_l)
    Φ = abs.(ϕ_l - ϕ_r)/2 + ϕ_l
    return Φ;
end

@doc raw"""
    f1(ϕ::Vector, u::Matrix, Xₙ::Vector, N::Int, M::Int)

Находит значение искомой функции на переходном слое ``f_1(t) = u(x_{t.p.}, t)`` путем поиска точки пересечения ``u(x,t)`` и ``ϕ(x)``.

Точка пересечения находится путем интерполяции функции ``u(x,t) - ϕ(x) = 0``.
"""
function f1(ϕ::Vector, u::Matrix, Xₙ::Vector, N::Int, M::Int)
    @assert length(Xₙ) == N+1
    @assert size(u) == (N+1, M+1)


    f1 = zeros(M+1);

    # Необходимо найти абсциссу пересечения Φ и u(x, t), т.е. аргумент х на каждом временном шаге
    for m in 1:M+1
        # Начальные условия при разработки программы заданы таким образом,
        # что в начале массива `V` образуются отрицательные элементы,
        # в конце — положительные
        V = u[:,m] - ϕ
        @assert V[1] < 0
        @assert V[end] > 0

        k = 1;
        # Найдем номер элемента, после которого стоит положительный элемент
        for n in 1:N
            if V[n]*V[n+1] < 0
                k = n;
                break;
            end
        end

        # try-catch Заглушка, на случай если разбиений по Х мало,
        # а переходный слой подошел слишком близко к границе
        # тогда k + 3 > N
        try
            # Выберем 8 ближайших элементов, по К.Р.А. от нуля функции
            x = Xₙ[k-3:k+4]
            y = V[k-3:k+4]

            # Теперь построим функцию ``x = v(y)``, она задана сеточными значениями выше
            # С помощью интерполяции найдем ``v(0)``
            # WARN интерполяция производится на неравномерной сетке: см массив `y`

            spl = Spline1D(y, x)

            f1[m] =spl(0);

        catch e
            isa(e, BoundsError) || rethrow(e)
            @warn "Bound Error" m k
            break
        end

    end

    return f1
end

function f2(f1::Vector, u::Matrix, Xₙ::Vector, N::Int, M::Int)
    @assert size(u) == (N+1, M+1)
    @assert length(f1) == M+1;

    f2 = similar(f1);

    # находим ординату на каждом временном шаге
    for m in 1:M+1
        k = 1;

        # ищем место в сетке по X
        for n in 1:N+1
            if Xₙ[n] > f1[m]
                k = n
                break
            end
        end

        # try-catch Заглушка, на случай если разбиений по Х мало,
        # а переходный слой подошел слишком близко к границе
        # тогда k + 3 > N
        try
        x = Xₙ[k-3:k+4]
        y = u[k-3:k+4,m]

        spl = Spline1D(x, y)

        f2[m] = spl(f1[m])

        catch e
            isa(e, BoundsError) || rethrow(e)
            @warn "Bound Error" m k
            break
        end
    end

    return f2
end



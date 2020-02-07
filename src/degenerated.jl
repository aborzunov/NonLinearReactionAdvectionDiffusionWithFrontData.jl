# Вычисление вырожденных корней

@doc raw"""
    phidetermination(q::Vector, ub::Vector,
                     Xₙ::Vector, N::Int,
                     Tₘ::Vector, M::Int)

Решает ОДУ для нахождения вырожденного корня.
"""
function phidetermination(q::Vector, ub::Vector,
                          Xₙ::Vector, N::Int,
                          Tₘ::Vector, M::Int)
    @assert length(q ) == N + 1
    @assert length(Xₙ) == N + 1

    @assert length(ub) == M + 1
    @assert length(Tₘ) == M + 1

    phi = zeros(N + 1, M + 1);

    for m in 1:M+1

        # Граничное условие на ϕ
        phi[1,m] = ub[m];

        for n = 1:N
            phi[n+1, m ] = phi[n, m ] + (Xₙ[n+1] - Xₙ[n]) * ( q[n+1] + q[n] ) / 2;
        end
    end

    return phi
end

@doc raw"""
    Φ(ϕ_l::Matrix, ϕ_r::Matrix, N::Int, M::Int) -> ::Vector

Вычисляет значение функции на переходном слое на каждом шаге по времени `m`
``|\phi_l^m - \phi_r^m|/2`` с помощью матриц вырожденных решений ``\phi_l`` и
``\phi_r`` вырожденного корня.

"""
function Φ(ϕ_l::Matrix, ϕ_r::Matrix, N::Int, M::Int)
    @assert size(ϕ_l) == (N+1, M+1)
    @assert size(ϕ_r) == (N+1, M+1)

    Φ = similar(ϕ_l)
    Φ = [ abs(ϕ_l[n, m] - ϕ_r[n, m])/2 + ϕ_l[n, m] for n in 1:N+1, m in 1:M+1]
    return Φ;
end

@doc raw"""
    f1(ϕ::Matrix, u::Matrix, Xₙ::Vector, N::Int, M::Int)

Находит значение искомой функции на переходном слое ``f_1(t) = u(x_{t.p.}, t)`` путем поиска точки пересечения ``u(x,t)`` и ``ϕ(x)``.

Точка пересечения находится путем интерполяции функции ``u(x,t) - ϕ(x) = 0``.
"""
function f1(ϕ::Matrix, u::Matrix, Xₙ::Vector, N::Int, M::Int)
    @assert length(Xₙ) == N+1
    @assert size(u) == (N+1, M+1)
    @assert size(ϕ) == (N+1, M+1)


    f1 = zeros(M+1);

    # Необходимо найти абсциссу пересечения Φ и u(x, t), т.е. аргумент х на каждом временном шаге
    for m in 1:M+1
        # Начальные условия при разработки программы заданы таким образом,
        # что в начале массива `V` образуются отрицательные элементы,
        # в конце — положительные
        V = u[:,m] - ϕ[:,m]
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



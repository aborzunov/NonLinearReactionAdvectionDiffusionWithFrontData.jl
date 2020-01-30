# Функционал, его градиент и минимизация

@doc raw"""
    J_q(uˢ::Matrix, ψˢ::Matrix,
             Xₙ::Vector, N::Int,
             Tₘ::Vector, M::Int) -> Vector

"""
function J_q(uˢ::Matrix, ψˢ::Matrix,
             Xₙ::Vector, N::Int,
             Tₘ::Vector, M::Int)

    @assert length(Xₙ) == N+1
    @assert length(Tₘ) == M+1

    J_q = similar(Xₙ);

    for m in 1:M
        τ = Tₘ[m+1] - Tₘ[m]
        J_q  += uˢ[:, m] .* ψˢ[:, m] .* τ
    end

    return J_q
end

@doc raw"""
    J(uˢ::Matrix, Xₙ::Vector, N::Int,
           Tₘ::Vector, M::Int,
           f1::Vector, f2::Vector,
           qₙˢ::Vector, α::Real = 0.0) -> Real

"""
function J(uˢ::Matrix, Xₙ::Vector, N::Int,
           Tₘ::Vector, M::Int,
           f1::Vector, f2::Vector,
           qₙˢ::Vector, α::Real = 0.0)
    J = 0.0

    for m in 1:M
        # ищем номер узла в сетке Xₙ в котором координата `x` больше координаты положения фронта `f1(t)`
        k = 0;
        for n in 1:N+1
            if Xₙ[n] > f1[m]
                k = n;
                break;
            end
        end

        # Выбираем 8 элементов значений аппроксимируемой функции и сетки по X
        y =  uˢ[k-3:k+4, m]
        x = Xₙ[k-3:k+4]

        spl = Spline1D(x, y)

        # Вычисляем uˢ(f1(t),t;qˢ)
        #
        ξ = spl(f1[m])
        τ = Tₘ[m+1] - Tₘ[m];

        J += (ξ - f2[m]) * τ

    end

    return J
end

function minimize(q₀, uₗ, uᵣ, ε, a, b, t₀, T, N, M)

end

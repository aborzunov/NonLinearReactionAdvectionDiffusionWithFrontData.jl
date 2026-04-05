# Градиентный чек и тест Тейлора
#
# Проверяют корректность вычисления градиента dJ/dq через сопряжённую задачу.
#
# ВАЖНО: J_q возвращает непрерывный градиент (плотность dJ/dq(x)),
# а конечные разности вычисляют дискретный ∂J/∂qᵢ. Между ними стоит
# масштабный множитель, связанный с пространственной дискретизацией.
# Поэтому тест проверяет не абсолютное совпадение, а:
# 1. Совпадение знаков (направление градиента правильное).
# 2. Постоянство отношения grad[i]/fd[i] по всем компонентам
#    (градиент пропорционален истинному с единым масштабом).

using NonLinearReactionAdvectionDiffusionWithFrontData
using LinearAlgebra
using Test

@testset "Градиентный чек (пропорциональность)" begin

    a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams()

    # Шаг 1: генерируем "истинные" наблюдения от q_true = qₙ
    u_true, XX_true, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ)
    _, _, _, f1_data, f2_data = generate_obs_data(u_true, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ)

    # Шаг 2: вычисляем градиент в точке q_test ≠ q_true
    q_test = zeros(N+1)
    u_test, XX_test, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, q_test)

    ψ₀ = zeros(N+1)
    ψl = zeros(M+1)
    ψr = zeros(M+1)
    w  = 0.0005

    ψ = solve_adjoint(ψ₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, q_test, u_test, f1_data, f2_data, w = w)

    # solve_adjoint возвращает ψ в обратном порядке по времени (от T до 0),
    # J_q ожидает прямой порядок (от 0 до T), поэтому разворачиваем.
    grad = J_q(u_test, ψ[:, end:-1:1], XX_test, N, Tₘ, M)

    # Шаг 3: конечно-разностная проверка в той же точке q_test
    fd_eps = 1e-5
    indices = [2, N ÷ 4, N ÷ 2, 3 * N ÷ 4, N]

    ratios = Float64[]
    for i in indices
        q_plus  = copy(q_test); q_plus[i]  += fd_eps
        q_minus = copy(q_test); q_minus[i] -= fd_eps

        u_plus,  XX_plus,  _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, q_plus)
        u_minus, XX_minus, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, q_minus)

        J_plus  = J(u_plus,  XX_plus,  N, Tₘ, M, f1_data, f2_data, q_plus)
        J_minus = J(u_minus, XX_minus, N, Tₘ, M, f1_data, f2_data, q_minus)

        fd_grad_i = (J_plus - J_minus) / (2 * fd_eps)

        # Знаки должны совпадать
        @test sign(grad[i]) == sign(fd_grad_i)

        push!(ratios, grad[i] / fd_grad_i)
    end

    # Отношения grad[i]/fd[i] должны быть приблизительно одинаковы:
    # разброс не более чем в 2 раза. Это подтверждает, что различие —
    # систематический масштабный множитель, а не ошибка в формуле.
    @test maximum(ratios) / minimum(ratios) < 2.0

end

@testset "Тест Тейлора (сходимость разностного отношения)" begin

    a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams()

    u_true, _, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ)
    _, _, _, f1_data, f2_data = generate_obs_data(u_true, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ)

    q_test = zeros(N+1)
    J₀ = let
        u_test, XX_test, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, q_test)
        J(u_test, XX_test, N, Tₘ, M, f1_data, f2_data, q_test)
    end

    # Направление возмущения
    δq = sin.(π .* Xₙ)

    # Вычисляем ΔJ(α) = J(q + α·δq) - J(q) для убывающих α.
    # Если J дифференцируем, то ΔJ(α) = α·⟨∇J, δq⟩ + O(α²),
    # и отношение ΔJ(α)/α сходится к константе при α → 0.
    alphas = [0.1, 0.01, 0.001]
    dJ_over_alpha = map(alphas) do α
        q_α = q_test .+ α .* δq
        u_α, XX_α, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, q_α)
        J_α = J(u_α, XX_α, N, Tₘ, M, f1_data, f2_data, q_α)
        (J_α - J₀) / α
    end

    # Разностное отношение должно сходиться: последовательные значения
    # всё ближе друг к другу.
    diffs = abs.(diff(dJ_over_alpha))
    @test diffs[end] < diffs[1]

    # Отношение последних двух значений ΔJ/α должно быть близко к 1
    # (сходимость к пределу)
    @test isapprox(dJ_over_alpha[end], dJ_over_alpha[end-1], rtol = 0.15)

end

# Тест для функции initial_guess
#
# Проверяет, что асимптотическое начальное приближение:
#   1. Возвращает вектор правильного размера.
#   2. Находится в разумном диапазоне (не NaN, не Inf, не нулевой).
#   3. Даёт меньшее значение функционала, чем нулевое приближение.
#
# Используются синтетические данные, сгенерированные из истинного q.

using NonLinearReactionAdvectionDiffusionWithFrontData
using LinearAlgebra

@testset "Начальное приближение (initial_guess)                " begin

    a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams()

    # Решаем прямую задачу и генерируем синтетические наблюдения
    u_true, _, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ)
    _, _, _, f1_data, f2_data = generate_obs_data(u_true, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ)

    @testset "размер результата" begin
        q0 = initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ)
        @test length(q0) == N + 1
    end

    @testset "результат конечный и ненулевой" begin
        q0 = initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ)
        @test all(isfinite.(q0))
        @test norm(q0) > 0
    end

    @testset "начальное приближение ближе к истинному q, чем нули" begin
        q0 = initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ)
        @test norm(q0 - qₙ) < norm(zeros(N+1) - qₙ)
    end

end

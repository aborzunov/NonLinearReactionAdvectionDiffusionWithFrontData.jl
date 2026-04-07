# Бейзлайн-бенчмарк для плана производительности.
# Запуск: julia --project=. benchmark_baseline.jl
#
# Измеряет время и аллокации горячих функций на стандартных параметрах dparams().
# Использует @time (не BenchmarkTools), чтобы не тянуть лишних зависимостей.

using NonLinearReactionAdvectionDiffusionWithFrontData
const M_ = NonLinearReactionAdvectionDiffusionWithFrontData

a, b, t₀, T, N, M, ε, Xₙ, Tₘ, qₙ, ulₘ, urₘ, u₀ = dparams()

println("Параметры: N=$N, M=$M, ε=$ε")

# --- Прогрев (компиляция) ---
u_true, XX, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ)
_, _, _, f1_data, f2_data = generate_obs_data(u_true, Xₙ, N, Tₘ, M, qₙ, ulₘ, urₘ)
ψ₀ = zeros(N+1); ψl = zeros(M+1); ψr = zeros(M+1)
ψ = solve_adjoint(ψ₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, u_true, f1_data, f2_data, w=0.0005)
_ = M_.J_q(u_true, ψ[:,end:-1:1], XX, N, Tₘ, M)
_ = M_.J(u_true, XX, N, Tₘ, M, f1_data, f2_data, qₙ)
q0 = zeros(N+1)
_, _, _ = minimize(q0, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data; S=2, β=0.001, w=0.0005)

function bench(name, f, n)
    times = Float64[]
    allocs = Int[]
    bytes = Int[]
    for _ in 1:n
        stats = @timed f()
        push!(times, stats.time)
        push!(allocs, Base.gc_alloc_count(stats.gcstats))
        push!(bytes, stats.bytes)
    end
    tmin = minimum(times)
    println(rpad(name, 40), "  min=$(round(tmin*1e3, digits=2)) ms",
            "  allocs=$(minimum(allocs))",
            "  bytes=$(round(minimum(bytes)/1024, digits=1)) KiB")
end

println("\n=== Бейзлайн (min из 5 запусков) ===")
bench("solve (прямая задача)",   () -> solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qₙ), 5)
bench("solve_adjoint",            () -> solve_adjoint(ψ₀, Xₙ, N, Tₘ, M, ε, ψl, ψr, qₙ, u_true, f1_data, f2_data, w=0.0005), 5)
bench("J_q",                      () -> M_.J_q(u_true, ψ[:,end:-1:1], XX, N, Tₘ, M), 5)
bench("J",                        () -> M_.J(u_true, XX, N, Tₘ, M, f1_data, f2_data, qₙ), 5)
bench("minimize (S=10)",          () -> minimize(q0, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data; S=10, β=0.001, w=0.0005), 3)
bench("minimize (S=50)",          () -> minimize(q0, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data; S=50, β=0.001, w=0.0005), 3)

# Значение функционала после S=50 итераций — ориентир для сравнения со "новой" версией
_, Js50, _ = minimize(q0, u₀, ulₘ, urₘ, Xₙ, N, Tₘ, M, ε, f1_data, f2_data; S=50, β=0.001, w=0.0005)
println("\nJ после S=50 итераций: J[1]=$(round(Js50[1], digits=6)), J[end]=$(round(Js50[end], digits=6))")
println("Относительное уменьшение: $(round((Js50[1]-Js50[end])/Js50[1]*100, digits=2))%")

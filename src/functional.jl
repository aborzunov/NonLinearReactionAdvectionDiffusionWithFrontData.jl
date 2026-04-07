# Функционал, его градиент и минимизация

@doc raw"""
    J_q(uˢ::Matrix, ψˢ::Matrix,
        Xₙ::Array, N::Int,
        Tₘ::Vector, M::Int,
        X_q::Vector = [NaN, NaN]) -> Vector

Вычисляет **непрерывную плотность градиента** функционала по параметру ``q``:
```math
\frac{dJ}{dq}(x) = -\int_0^T u^s(x, t)\, \psi^s(x, t)\, dt
```

!!! warning "Масштаб непрерывный vs дискретный"
    Возвращается **непрерывная плотность** ``dJ/dq(x)``, а не дискретный вектор ``\partial J/\partial q_i``.
    Конечно-разностная проверка даёт систематический множитель ~``h`` (шаг сетки) между ними.
    Тесты в `gradient_check.jl` проверяют **пропорциональность** и совпадение знаков, а не абсолютное равенство.

!!! warning "Направление времени ψ"
    `solve_adjoint` возвращает ``\psi`` в **обратном** порядке по времени.
    Перед передачей в `J_q` необходимо инвертировать: `ψ[:, end:-1:1]`.

# Arguments
- `uˢ::Matrix`:  Решение прямой задачи, размер `(N+1, M+1)`.
- `ψˢ::Matrix`:  Решение сопряженной задачи (в **прямом** порядке по времени), размер `(N+1, M+1)`.
- `Xₙ::Array`:   Сетка по пространству — вектор `(N+1,)` или матрица `(N+1, M+1)` для динамической сетки.
- `N::Int`:       Число интервалов по пространству.
- `Tₘ::Vector`:  Сетка по времени, размер `M+1`.
- `M::Int`:       Число интервалов по времени.
- `X_q::Vector`: Вспомогательная сетка для ``q`` при динамической сетке `Xₙ`. При статической — не нужна.

# Return
Вектор — непрерывная плотность градиента ``dJ/dq`` на сетке `Xₙ[:,1]` (или `X_q`).
"""
function J_q(uˢ::AbstractMatrix, ψˢ::AbstractMatrix,
             Xₙ::Matrix, N::Int,
             Tₘ::Vector, M::Int,
             X_q::Vector = [NaN, NaN])

    size(Xₙ) == (N+1, M+1) ||
    throw(ArgumentError("size(Xₙ) == $(size(Xₙ)), N == $(N), M == $(M) " *
                        "Массив Xₙ должен иметь размерность (N+1, M+1)."))

    isDynamicMesh = ! isapprox(Xₙ[:, 1], Xₙ[:, end])
    # Когда мы понимаем, что у нас динамическая сетка,
    # мы должны убедиться, что нам передали необходимую вспомогательную сетку.
    if isDynamicMesh
        if all(isnan.(X_q))
            throw(ArgumentError("Xₙ передан в виде матрицы, значит подразумевается," *
                                "\nчто используется динамическая сетка. Это требует" *
                                "\nпередачи вспомогательной сетки `X_q`, для ``q``"))
        end
    end
    @assert length(Tₘ) == M+1


    # заполним рабочие массивы специальным образом
    # Если сетка динамическая,
    # Нам нужна интерполяция решений на вспомогательную сетку `X_q`
    if isDynamicMesh
        J_q = zero(X_q);                # Градиент будет определен на вспомогательной сетке
        u = zeros(length(X_q), M+1);
        ψ = zeros(length(X_q), M+1);
        for m in 1:M+1
            X = Xₙ[:, m];
            uspl = Spline1D(X, uˢ[:, m]);
            ψspl = Spline1D(X, ψˢ[:, m]);
            u[:, m] = uspl(X_q);
            ψ[:, m] = ψspl(X_q);
        end
    else                        # Если все вычислялось на статических сетках
        u = uˢ;                 # То будем работать с исходными массивами
        ψ = ψˢ;
        J_q = zeros(N+1);       # Градиент определяется на сетке Xₙ
    end

    #  Здесь уже не важно, на какой сетке определены `u`, `ψ`
    #  Вычислим интеграл по формуле трапеций.
    #  ``\int_0^T u^s(x, t) \psi^s(x, t) dt ``
    #  Параллелим по времени: каждый @spawn-таск получает свой chunk и локальный буфер.
    #  Избегаем threadid() — в Julia 1.9+ он нестабилен при :dynamic scheduling.
    #  grad_len вычисляется до замыканий, чтобы не захватывать J_q (одноимённая функция).
    grad_len = length(J_q)
    nt = Threads.nthreads()
    chunk = max(1, div(M, nt))
    tasks = map(1:nt) do t
        m_lo = (t - 1) * chunk + 1
        m_hi = t == nt ? M : t * chunk
        Threads.@spawn begin
            local_buf = zeros(grad_len)
            @inbounds for m in m_lo:m_hi
                τ = Tₘ[m+1] - Tₘ[m]
                @views @. local_buf += (u[:, m] * ψ[:, m] + u[:, m+1] * ψ[:, m+1]) * τ / 2
            end
            local_buf
        end
    end
    for task in tasks
        J_q .+= fetch(task)
    end

    @. J_q = -J_q
    return J_q
end

@doc raw"""
    J(uˢ::Matrix, Xₙ::Array, N::Int,
      Tₘ::Vector, M::Int,
      f1::Vector, f2::Vector,
      qₙˢ::Vector, α::Real = 0.0) -> Real

Вычисляет функционал `` J(\mathbf{x}) = \int_0^T \left( u(f_1(t), t; q^s) - f_2(t) \right)^2 +
\alpha \int_0^1 q^2(x) dx ``.

Использует [`f2`](@ref) для вычисления ``u(f_1(t), t)`` и после по формуле трапеций.
"""
function J(uˢ::Matrix, Xₙ::Array, N::Int,
           Tₘ::Vector, M::Int,
           f1::Vector, f2::Vector,
           qₙˢ::Vector, α::Real = 0.0)
    J = 0.0
    ξ = NonLinearReactionAdvectionDiffusionWithFrontData.f2(f1, uˢ, Xₙ, N, M);

    @inbounds for m in 1:M
        τ = Tₘ[m+1] - Tₘ[m];
        J += ( (ξ[m] - f2[m])^2 + (ξ[m+1] - f2[m+1])^2) * τ / 2
    end

    # Тихоновская регуляризация: α ∫ q²(x) dx по формуле трапеций на сетке Xₙ[:,1]
    if α > 0
        X = Xₙ isa Matrix ? view(Xₙ, :, 1) : Xₙ
        reg = 0.0
        @inbounds for n in 1:N
            h = X[n+1] - X[n]
            reg += (qₙˢ[n]^2 + qₙˢ[n+1]^2) * h / 2
        end
        J += α * reg
    end

    return J
end

@doc raw"""
    minimize(q₀::Vector, u₀::Vector,
             ulₘ::Vector, urₘ::Vector,
             Xₙ, N,
             Tₘ::Vector, M,
             ε,
             f1_data::Vector, f2_data::Vector;
             S::Int = 10,
             β::Real = 0.01,
             w::Real = 0.0001,
             α::Real = 0.0,
             linesearch::Bool = true,
             tol_grad::Real = 0.0,
             tol_J::Real = 0.0,
             ls_c::Real = 1e-4,
             ls_shrink::Real = 0.5,
             ls_max::Int = 10,
             β_min::Real = 1e-12,
             create_mesh::Function = x -> [NaN, NaN]) -> Vector, Vector, Matrix


Градиентный спуск с опциональным backtracking line search (Armijo) для минимизации функционала.

# Keyword Arguments
- `S`:           Максимальное число итераций.
- `β`:           Начальный размер шага градиентного спуска.
- `w`:           Ширина аппроксимации δ-функции для сопряжённой задачи.
- `α`:           Параметр тихоновской регуляризации (добавляет `α ∫ q² dx` к функционалу и `2α·q` к градиенту).
- `linesearch`:  Если `true` (по умолчанию), используется backtracking Armijo для подбора `β` на каждом шаге.
- `tol_grad`:    Остановка по норме градиента ``\|∇J\| < tol_grad``. По умолчанию `0` — не использовать.
- `tol_J`:       Остановка по относительному изменению функционала ``|ΔJ| < tol_J``. По умолчанию `0`.
- `ls_c`:        Константа в условии Armijo (по умолчанию 1e-4).
- `ls_shrink`:   Множитель уменьшения шага при backtracking (по умолчанию 0.5).
- `ls_max`:      Максимальное число backtracking-попыток на одной итерации.
- `β_min`:       Минимальный допустимый шаг; если `β` упал ниже — считается, что спуск невозможен.

# Return
Кортеж `(qˢ, J_values, Q_values)`. При ранней остановке `J_values`/`Q_values` обрезаны до фактического числа итераций.
"""
function minimize(q₀::Vector, u₀::Vector,
                  ulₘ::Vector, urₘ::Vector,
                  Xₙ, N,
                  Tₘ::Vector, M,
                  ε,
                  f1_data::Vector, f2_data::Vector;
                  S::Int = 10,
                  β::Real = 0.01,
                  w::Real = 0.0001,
                  α::Real = 0.0,
                  linesearch::Bool = true,
                  tol_grad::Real = 0.0,
                  tol_J::Real = 0.0,
                  ls_c::Real = 1e-4,
                  ls_shrink::Real = 0.5,
                  ls_max::Int = 10,
                  β_min::Real = 1e-12,
                  create_mesh::Function = x -> [NaN, NaN],
                  showProgress = false)

    ψ₀ = zero(Xₙ);
    ψl = zero(Tₘ);
    ψr = zero(Tₘ);

    isDynamicMesh = all(isnan.(create_mesh(Xₙ[end - div(end,2)])))
    @info isDynamicMesh ? "Используем оригинальную сетку Xₙ для q₀" : "Используем отдельную сетку для q₀"

    # Создадим интерполяционный объект для q
    qspl = Spline1D(Xₙ, q₀)
    k = 100;                                 # Кол-во интервалов в вспомогательной сетке
    if isDynamicMesh
        X_q = Xₙ;
    else
        X_q = [ first(Xₙ) + n * (last(Xₙ) - first(Xₙ))/Float64(k) for n in 0:k]
    end
    q_aux = qspl(X_q);

    if ! isapprox(qspl.(Xₙ), q₀)
        throw(ArgumentError("Сеточные значения q₀ определены на какой-то другой сетке, а не Xₙ"))
    end

    #' ## Подготовка к итерационному процессу
    qˢ = copy(q₀);
    J_values = zeros(S);
    Q_values = zeros(N+1, S);

    p = showProgress ? Progress(S, 5, "Iterating minimization loop S=$(S)... ") : nothing;

    # Начальное решение прямой задачи и значение функционала в точке q₀
    uˢ, XXˢ, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qˢ, create_mesh = create_mesh);
    J_current = J(uˢ, XXˢ, N, Tₘ, M, f1_data, f2_data, qˢ, α);

    last_s = 0
    q_aux_old = similar(q_aux)
    last_β = float(β)          # фактический шаг, принятый line search на прошлой итерации

    for s in 1:S
        showProgress && next!(p; showvalues = [
            (:iter, "$(s)/$(S)"),
            (:J, J_current),
            (:β, last_β),
        ]);

        Q_values[:, s] .= qˢ
        J_values[s] = J_current
        last_s = s

        # Решаем сопряженную задачу в текущей точке
        ψˢ = solve_adjoint(ψ₀, XXˢ, N, Tₘ, M, ε, ψl, ψr, qˢ, uˢ, f1_data, f2_data, w=w)

        # Градиент: J_q принимает ψ в прямом порядке по времени.
        # Передаём view, чтобы избежать копии полной матрицы.
        ∇J = J_q(uˢ, @view(ψˢ[:, end:-1:1]), XXˢ, N, Tₘ, M, X_q)

        # Добавляем градиент тихоновской регуляризации: d/dq (α ∫ q² dx) = 2α·q
        if α > 0
            @. ∇J += 2 * α * q_aux
        end

        grad_norm_sq = sum(abs2, ∇J)
        if tol_grad > 0 && sqrt(grad_norm_sq) < tol_grad
            break
        end

        # --- Backtracking Armijo line search ---
        copyto!(q_aux_old, q_aux)
        β_try = float(β)
        local uˢ_new, XXˢ_new, qˢ_new, J_new
        uˢ_new = uˢ; XXˢ_new = XXˢ; qˢ_new = qˢ; J_new = J_current

        if linesearch && Threads.nthreads() > 1
            # Parallel Armijo probe: запускаем до ls_max кандидатов одновременно
            # β, β·shrink, β·shrink², … — каждый в отдельном таске.
            # Выбираем наибольший β, при котором условие Armijo выполнено.
            β_candidates = Float64[]
            β_k = β_try
            while β_k >= β_min && length(β_candidates) < ls_max
                push!(β_candidates, β_k)
                β_k *= ls_shrink
            end
            if isempty(β_candidates)
                copyto!(q_aux, q_aux_old)
                @goto finish
            end

            # Каждый таск получает свою копию q_try и не трогает общее состояние
            probe_tasks = map(β_candidates) do β_k
                q_try = q_aux_old .- β_k .* ∇J
                Threads.@spawn begin
                    spl   = Spline1D(X_q, q_try)
                    q_s   = spl(Xₙ)
                    u_s, XX_s, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, q_s,
                                         create_mesh = create_mesh)
                    J_s   = J(u_s, XX_s, N, Tₘ, M, f1_data, f2_data, q_s, α)
                    (J_s, q_s, u_s, XX_s, q_try)
                end
            end

            probe_results = map(fetch, probe_tasks)

            # Наибольший β, удовлетворяющий условию Armijo (кандидаты в порядке убывания β)
            accepted = findfirst(eachindex(β_candidates)) do k
                probe_results[k][1] ≤ J_current - ls_c * β_candidates[k] * grad_norm_sq
            end

            if accepted === nothing
                copyto!(q_aux, q_aux_old)
                @goto finish
            end
            J_new, qˢ_new, uˢ_new, XXˢ_new, q_aux_accepted = probe_results[accepted]
            copyto!(q_aux, q_aux_accepted)
            β_try = β_candidates[accepted]

        else
            # Последовательный backtracking Armijo (фоллбэк при 1 потоке)
            while true
                @. q_aux = q_aux_old - β_try * ∇J
                qs_spl = Spline1D(X_q, q_aux)
                qˢ_new = qs_spl(Xₙ)

                uˢ_new, XXˢ_new, _ = solve(u₀, Xₙ, N, Tₘ, M, ε, ulₘ, urₘ, qˢ_new,
                                            create_mesh = create_mesh)
                J_new = J(uˢ_new, XXˢ_new, N, Tₘ, M, f1_data, f2_data, qˢ_new, α)

                if !linesearch
                    break  # простой шаг без проверки
                end
                # Условие Armijo: J(q_new) ≤ J(q) − c·β·‖∇J‖²
                if J_new ≤ J_current - ls_c * β_try * grad_norm_sq
                    break
                end
                β_try *= ls_shrink
                if β_try < β_min
                    copyto!(q_aux, q_aux_old)
                    @goto finish
                end
            end
        end

        # Принимаем шаг
        qˢ = qˢ_new
        uˢ = uˢ_new
        XXˢ = XXˢ_new
        ΔJ = J_current - J_new
        J_current = J_new
        last_β = β_try

        if tol_J > 0 && abs(ΔJ) < tol_J
            break
        end
    end
    @label finish

    # Обрезаем историю до фактического числа итераций
    if last_s < S
        J_values = J_values[1:last_s]
        Q_values = Q_values[:, 1:last_s]
    end

    return (qˢ, J_values, Q_values)
end

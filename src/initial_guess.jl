# Вычисление начального приближения с использованием асимптотических методов

function front_velocity(f1_data, Tₘ, M)

    df1dt = zero(f1_data);
    τ = Tₘ[2] - Tₘ[1];

    df1dt[1] = (-3 * f1_data[1] + 2 * f1_data[2] - f1_data[3]) / (2 * τ);
    for m in 2:M-2
        df1dt[m] = (f1_data[m+1] - f1_data[m-1]) / (2 * τ);
    end
    df1dt[M+1] = (f1_data[M+1 - 2] - 2 * f1_data[M+1 - 1] + 3 * f1_data[M+1]) / (2 * τ);

    return df1dt;
end

# Найдем минимальный и максимальный номер узлов сетки Xₙ, которые пересекает
# переходный слой
# .     Узлы сетки
# x     f1(t) положение переходного слоя
#                     k                                   K
# .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
#                   x x x x x x x x x x x x x x x x x x x x x
function get_nodes_range_where_fron_passed(x, Xₙ)
    k = findfirst(z -> z > minimum(x), Xₙ)
    K = findfirst(z -> z > maximum(x), Xₙ) - 1

    # Первые и последние точки, часто оказываются плохими
    # Отрежем на всякий случай
    return k, K
end

function get_A_matrix(k, K, Xₙ, N)
    h = Xₙ[2] - Xₙ[1];

    Al = LowerTriangular(h * ones(K-k+1, K-k+1));
    Ar = UpperTriangular(-h * ones(K-k+1, K-k+1));
    Al[:, 1]    .= (k-1)*h;
    Ar[:, end]  .= (N+1-K)*h;
    # Повторим чертеж для объяснения двух прошлых строк
    # N = 17, тогда нумерация массива будет следующей
    # k = 6, K = 16
    # Интервалов лежащих слева от k -- 5, (k-1)
    # Интервалов лежащих справа от K -- 2, (N+1 - K)
    # 1   2   3   4   5   6   7   8   9   11  12  13  14  15  16  17  18
    #                     k                                   K
    # .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
    #                   x x x x x x x x x x x x x x x x x x x x x

    return A = -(Al + Ar) / 2.0;
end

function get_B_vector(k, K, Vn, ulₘ, urₘ)
    # `B` составим из ненулевых v_n
    B = Vn[k:K];
    B = B .- ( -(ulₘ[1] + ulₘ[1])  / 2.0);

    return B;
end

function solve_system_for_initial_guess(k, K, N, A, B, α = 0.001)

    qx = ( transpose(A) * A + α*I) \ transpose(A)*B

    l = [qx[2] for i in 1:k-1];     # Дополняем начальное приближение константами
    r = [qx[end-1] for i in K+1:N+1]; # Дополняем Справа
    qx[end] = qx[end-1];
    qx[1] = qx[2];

    q_guess = [l; qx; r];
end

function initial_guess(f1_data, Xₙ, N, Tₘ, M, ulₘ, urₘ, α = 0.001)

    df1dt = front_velocity(f1_data, Tₘ, M);

    # Последние точки некрасивые, обрежем их
    x = f1_data[2:end-3];
    y = df1dt[2:end-3];
    tm = Tₘ[2:end-3];

    k, K = get_nodes_range_where_fron_passed(x, Xₙ);
    @info "Фронт прошел точки $((k, K)), $(div((100*(K-k)),(N+1)))% кадра"

    # Найдем скорость фронта как функцию координат
    # Иначе говоря, найдем скорость, которой обладал фронт проходя узел сетки
    vspl = Spline1D(x, y);
    Vn = zero(Xₙ);
    for n in k:K
            Vn[n] = vspl(Xₙ[n])
    end

    A = get_A_matrix(k, K, Xₙ, N);
    B = get_B_vector(k, K, Vn, ulₘ, urₘ);
    q_guess = solve_system_for_initial_guess(k, K, N, A, B, α);

    return q_guess;
end

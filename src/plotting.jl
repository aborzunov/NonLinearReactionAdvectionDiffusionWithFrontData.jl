# Файл содержит функции, помогающие быстро нарисовать полученное решение

@doc raw"""
    make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector,
                  ϕ_l::Matrix = missings(2), ϕ_r::Matrix = missings(2),
                  f1::Vector = missings(2), f2::Vector = missings(2),
                  analitical = nothing;
                  frames_to_write::Int = -1,
                  name = "solution.gif",
                  convert2mp4 = false)

Рисует gif анимацию решения.
вплоть по `frames_to_write`-ый кадр, сохраняет как "results/`name`".

Так же рисует аналитическое решение `analitic(x,t)`, если таково передано.
`label`::String — LaTeX строка подписи искомой функции с экранирование спец. символов.

!!! tip
    Pass an empty string to avoid saving at disk.

TODO: Fix doc
"""
function make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector,
                  ϕ_l::Matrix = missings(2), ϕ_r::Matrix = missings(2),
                  f1::Vector = missings(2), f2::Vector = missings(2),
                  analitical = nothing;
                  frames_to_write::Vector = Vector(),
                  name = "solution.gif", convert2mp4 = false, label="u")
    N,M = size(u) .-1
    if frames_to_write == Vector()
        frames_to_write = collect(1:M+1);
    end

    a = Animation()
    for m in frames_to_write
        make_plot(u, Xₙ, Tₘ, m, ϕ_l, ϕ_r, f1, f2, analitical; label=label)
        frame(a)
    end

    if convert2mp4
        g = mp4(a, replace(name, "gif" => "mp4"), show_msg=false)
    else
        g = gif(a, name, show_msg=false)
    end

    return g
end

@doc raw"""
    make_plot(u::Matrix, Xₙ::Vector, Tₘ::Vector, m::Int,
                   ϕ_l::Matrix = missings(2), ϕ_r::Matrix = missings(2),
                   f1::Vector = missings(2), f2::Vector = missings(2),
                   analitical = nothing;
                   label::String = "u")

Рисует `m`-ый кадр решения `u`. `Xₙ, Tₘ` — сетки.
`ϕ_l, ϕ_r` — вырожденные решения.
`f1, f2` — сгенерированные априорные данные.
`analitical` — или функция или сеточные значения аналитического решения.
`label`::String — LaTeX строка подписи искомой функции с экранирование спец. символов.
"""
function make_plot(u::Matrix, Xₙ::Vector, Tₘ::Vector, m::Int,
                   ϕ_l::Matrix = missings(2), ϕ_r::Matrix = missings(2),
                   f1::Vector = missings(2), f2::Vector = missings(2),
                   analitical = nothing;
                   label = "u")

    N, M = size(u) .- 1;
    yl = extrema(u[:,:]).*1.05;

    # График, оси, подписи осей и пр.
    pl = plot(size=(800, 600), xlabel="x", ylabel=latexstring(label, "(x)"), ylims=yl)

    # Начальное условие
    plot!(Xₙ, u[:,1], line=:dash, label=latexstring(label,"(x, 0)"))

    # Найденное решение
    plot!(Xₙ, u[:,m], label=latexstring(label,"(x, t)"), color=:blue)

    # Точки сетки на найденной функции и их проекция на ось Х
    scatter!(Xₙ, u[:,m], color=:blue, label="", markersize=3)
    scatter!(Xₙ, [0 for x in Xₙ], color=:black, label="", markersize=2)

    # Надпись слева внизу с текущим временем
    annotate!(0.0, 0.85*first(yl), Plots.text(@sprintf("t = %.2f",Tₘ[m]), 14, :left ))

    check(x::Array) = !any(ismissing.(x))
    if ( check(ϕ_l) ) && ( check(ϕ_r)) && ( check(f1) ) && ( check(f2) )
        ϕ = Φ(ϕ_l, ϕ_r, N, M);
        plot!(Xₙ, ϕ_l[:, m], label=L"\phi_l", color=:darkgoldenrod)
        plot!(Xₙ, ϕ_r[:, m], label=L"\phi_r", color=:darkgoldenrod)
        plot!(Xₙ, ϕ[:, m], label=L"\widetilde{\Phi}", color=:gold)

        # Плавающая по оси Y вслед за пунктиром надпись
        buff = latexstring("f_2", @sprintf("(t) = %.2f",f2[m]));
        annotate!(0.0, f2[m] * 1.25, Plots.text(buff, :left))
        # Горизонтальный пунктир
        plot!([0, f1[m]], [f2[m], f2[m]], line=:dash, color=:black, label="")
        # Красная точка слева, около подписи f2
        scatter!( [0], [f2[m]], color=:red, label="")

        # Плавающая по оси X вслед за пунктиром надпись
        buff = latexstring("f_1" , @sprintf("(t) = %.2f",f1[m]));
        annotate!(f1[m] * 1.05, 0.95 * first(yl), Plots.text(buff, :left))
        # Вертикальная линия
        plot!([f1[m], f1[m]], [yl[1], f2[m]], line=:dash, color=:black, label="")
        # Красная точка снизу, около надписи f_1
        # Значение по Y умножаем на нормировочный коэффициент, чтобы точка полностью влезла в кадр
        scatter!( [f1[m]], [yl[1]*0.99], color=:red, label="")

        # красная точка на пересечении пунктиров априорной информации
        scatter!( [f1[m]], [f2[m]], color=:red, label="")
    end

    # Если передали аналитическое решение
    if analitical != nothing
        # В виде функции, применяем её к нашим сеткам
        if typeof( analitical ) <: Function
            plot!(Xₙ, analitical.(Xₙ, Tₘ[m]), label="analitical(x,t)", color=:green, linewidth = 5, alpha=0.3)
        elseif typeof( analitical) <: Matrix
            # В виде матрицы сеточных значений, рисуем m-ый шаг по времени их неё
            if size(analitical) != size(u)
                throw(ArgumentError("Размер `analitical` != $(size(u))"))
            end
            plot!(Xₙ, analitical[:,m], label="analitical(x,t)", color=:green, linewidth = 5, alpha=0.3)
        else
            throw(ArgumentError("Unknown typeof(analitical): $(typeof(analitical))"))
        end
    end

    return pl
end

# Файл содержит функции, помогающие быстро нарисовать полученное решение

@doc raw"""
    make_draft(

Рисует экскиз решения: в первый, в последний момент времени и момент вермени
посередине.
"""
function draft(u, Xₙ, N, Tₘ, M; title = "")

    # Набросаем эскиз решения прямой задачи
    psol = plot(title=title, size = (800, 800));
    plot!(Xₙ, u[:, 1], label="u(x, 0)");
    plot!(Xₙ, u[:, div(end,2)], label="u(x, T/2)");
    plot!(Xₙ, u[:, end], label="u(x, T)");

    return psol
end

@doc raw"""
    minimization_draft(qₙ::Vector, Qs::Matrix,
                       Xₙ::Vector, N::Integer,
                       Js::Vector;
                       zoom::Bool = false,
                       zoom_chunk::Rational = 7//8,
                       annotate_string::String = "")

Рисует один график с результатами минимизации.

Строит истинное, начальное, среднее и последнее приближение ``q``.
Значение функционала на всём пути. По требованию `zoom = true` строит функционал
ближе к концу процесса минимизации. Доля увеличинного графика устанавливается
`zoom_chunk`-ом, график будет построен от `zoom_chunk * S` до конца.
"""
function minimization_draft(qₙ::Vector, Qs::Matrix,
                                 Xₙ::Vector, N::Integer,
                                 Js::Vector;
                                 zoom::Bool = false,
                                 zoom_chunk::Rational = 7//8,
                                 annotate_string::String = "")

    S = size(Js, 1)
    if (zoom_chunk <= 0.) || (zoom_chunk >= 1.0)
        throw(ArgumentError("`zoom_chunk` должен лежать в интервале [0,1]"));
    end

    # Восстанавливаемая функция, будет слева
    pQs = plot(Xₙ, qₙ, linewidth=2, label="Истинное", title = "q(x)", size = (800, 800));
    pQs = plot!(Xₙ, Qs[:, 1], line=:dash, color=:green, label="Начальное");
    pQs = plot!(Xₙ, Qs[:, div(S, 2)], label="на полпути");
    pQs = plot!(Xₙ, Qs[:, end], label="Последнее")

    # Функционал на всём пути минимизации
    pJs = plot(Js, ylims = (0, maximum(Js)), title = "Значение функционала", size = (800, 800))
    pJs = annotate!( div(S, 10), 7*maximum(Js)/8, Plots.text(annotate_string, :left))

    # Если нужно, добавляем третий график с функционалам в конце пути
    if zoom
        range_to_zoom = div(numerator(zoom_chunk) * S,
                            denominator(zoom_chunk)):S;     # Выбираем участок
        zoommed_data = Js[range_to_zoom];                   # И данные на нём

        pJs_zoom = plot(range_to_zoom, zoommed_data, ylims = (0, maximum(zoommed_data)), size = (800, 800));
        return pQs, pJs, pJs_zoom;
    else
        return pQs, pJs;
    end
end


@doc raw"""
    make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector,
                  ϕ_l::Matrix = missings(2), ϕ_r::Matrix = missings(2),
                  ϕ::Matrix = missings(2,2),
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
function make_gif(u::Matrix, Xₙ::Array, Tₘ::Vector,
                  ϕ_l::Matrix = missings(2,2), ϕ_r::Matrix = missings(2,2),
                  ϕ::Matrix = missings(2,2),
                  f1::Vector = missings(2), f2::Vector = missings(2),
                  analitical = nothing;
                  frames_to_write::Vector = Vector(),
                  name = "solution.gif", convert2mp4 = false, label="u",
                  title::String = "")

    N,M = size(u) .-1
    if frames_to_write == Vector()
        frames_to_write = collect(1:div(M,80):M+1);
    end

    a = Animation()
    @showprogress "Plotting animation..." for m in frames_to_write
        X = typeof(Xₙ) <: Vector ? Xₙ : Xₙ[:, m]

        make_plot(u, X, Tₘ, m, ϕ_l, ϕ_r, ϕ, f1, f2, analitical; label=label,
                  y_lim = extrema( u[:,1:last(frames_to_write)] * 1.05),
                  title = title)

        frame(a)
    end

    convert2mp4 = endswith(name, ".mp4");

    g = convert2mp4 ? mp4(a, name, show_msg=false) : gif(a, name, show_msg=false)

    return g
end

@doc raw"""
    make_plot(u::Matrix, X::Vector, Tₘ::Vector, m::Int,
                   ϕ_l::Matrix = missings(2), ϕ_r::Matrix = missings(2),
                   f1::Vector = missings(2), f2::Vector = missings(2),
                   analitical = nothing;
                   label::String = "u")

Рисует `m`-ый кадр решения `u`. `X, Tₘ` — сетки.
`ϕ_l, ϕ_r` — вырожденные решения.
`f1, f2` — сгенерированные априорные данные.
`analitical` — или функция или сеточные значения аналитического решения.
`label`::String — LaTeX строка подписи искомой функции с экранирование спец. символов.
"""
function make_plot(u::Matrix, X::Vector, Tₘ::Vector, m::Int,
                   ϕ_l::Matrix = missings(2, 2), ϕ_r::Matrix = missings(2, 2),
                   ϕ::Matrix = missings(2,2),
                   f1::Vector = missings(2), f2::Vector = missings(2),
                   analitical = nothing;
                   label = "u",
                   y_lim = missings(2),
                   title::String)

    check(x::Array) = !any(ismissing.(x))
    check(x::Tuple) = !any(ismissing.(x))

    N, M = size(u) .- 1;
    if check(y_lim)
        yl = y_lim;
    else
        yl = extrema(u[:,:]).*1.05;
    end

    # График, оси, подписи осей и пр.
    pl = plot(size=(800, 600), xlabel="x",
              ylabel=latexstring(label, "(x)"),
              ylims=yl,
             title = title)

    # Начальное условие
    #plot!(X, u[:,1], line=:dash, label=latexstring(label,"(x, 0)"))

    # Найденное решение
    plot!(X, u[:,m], label=latexstring(label,"(x, t)"), color=:blue)


    scatter!(X, u[:,m], color=:blue, label="", markersize=3)               # Точки сетки на найденной функции
    scatter!(X, [0 for x in X], color=:black, label="", markersize=2)     # Сетка по Х

    # Надпись слева внизу с текущим временем
    annotate!(0.0, 0.85*first(yl), Plots.text(@sprintf("t = %.2f",Tₘ[m]), 14, :left ))

    if ( check(ϕ_l) ) && ( check(ϕ_r)) && check(ϕ) && ( check(f1) ) && ( check(f2) )

        # Рисуем вырожденные корни
        plot!(X, ϕ_l[:, m], label=L"\phi_l", color=:darkgoldenrod)
        plot!(X, ϕ_r[:, m], label=L"\phi_r", color=:darkgoldenrod)
        plot!(X, ϕ[:, m], label=L"\widetilde{\Phi}", color=:gold)

        # Плавающая по оси Y вслед за пунктиром надпись
        buff = latexstring("f_2", @sprintf("(t) = %.2f",f2[m]));
        annotate!(0.0, f2[m] * 1.25, Plots.text(buff, :left))
        # Горизонтальный пунктир
        plot!([0, f1[m]], [f2[m], f2[m]], line=:dash, color=:black, label="")
        # Красная точка слева, около подписи f2
        scatter!( [0], [f2[m]], color=:red, label="")

        # история f1, f2
        scatter!(f1[1:div(M+1, 80)+1:m], f2[1:div(M+1, 80)+1:m], color=:red, alpha = 0.3, label="")

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
            plot!(X, analitical.(X, Tₘ[m]), label="analitical(x,t)", color=:green, linewidth = 5, alpha=0.3)
        elseif typeof( analitical) <: Matrix
            # В виде матрицы сеточных значений, рисуем m-ый шаг по времени их неё
            if size(analitical) != size(u)
                throw(ArgumentError("Размер `analitical` != $(size(u))"))
            end
            plot!(X, analitical[:,m], label="analitical(x,t)", color=:green, linewidth = 5, alpha=0.3)
        else
            throw(ArgumentError("Unknown typeof(analitical): $(typeof(analitical))"))
        end
    end

    return pl
end

@doc raw"""
    make_minimzation_gif(Js::Vector, Qs::Matrix,
                         qₙ::Vector, Xₙ::Vector;
                         frames_to_write::Vector = Vector(),
                         name = "solution.gif", convert2mp4 = false,
                         β::Real = 0.0)

Сохраняет анимацию процесса минимизации.
"""
function make_minimzation_gif(Js::Vector, Qs::Matrix,
                              qₙ::Vector, Xₙ::Vector;
                              frames_to_write::Vector = Vector(),
                              name = "solution.gif", convert2mp4 = false,
                              β::Real = 0.0)

    S = length(Js)
    if frames_to_write == Vector()
        frames_to_write = collect(1:S);
    end

    yl = extrema([Qs qₙ])

    a = Animation()
    @showprogress "Composing minimization.." for s in frames_to_write
        pQs = plot(xlabel = "x", ylabel="q(x)", ylims=yl, size = (800, 800) )
        pQs = plot!(Xₙ, qₙ, label="q(x)")
        pQs = scatter!(Xₙ, Qs[:,s], ylims = (-1.1, 1.1), xlims = (0, 1), title="Искомая qˢ(x) при s = $(s)", label=L"q^s(x)")
        pQs = plot!(Xₙ, Qs[:,1], line = :dash, label=L"q^0(x)")

        pJs = plot(title="Значение функционала на шаге s = $(s)", size = (800, 800),
                   xlabel = "s", ylabel=L"J", xlims = (1,S), ylimits = extrema(Js))
        pJs = plot!(yaxis=:log)
        pJs = plot!(1:s, Js[1:s], label=L"J(q^s)")

        β == 0 || (pJs = annotate!(S/10, minimum(Js)*1.2, "\\beta = $(β)"))
        #pJs = annotate!(S/10, minimum(Js)*1.5, "S = $(S)")

        p = plot(pQs, pJs, size = (1600, 800) )
        frame(a);
    end

    if convert2mp4
        g = mp4(a, replace(name, "gif" => "mp4"), show_msg=false)
    else
        g = gif(a, name, show_msg=false)
    end
    return g

end

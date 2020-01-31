# Файл содержит функции, помогающие быстро нарисовать полученное решение

@doc raw"""
    make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector,
                  ϕ_l::Vector = missings(2), ϕ_r::Vector = missings(2),
                  f1::Vector = missings(2), f2::Vector = missings(2),
                  analitical = nothing;
                  frames_to_write::Int = -1, frame_skip::Int=-1,
                  name = "solution.gif", convert2mp4 = false)

Рисует gif анимацию решения каждый `frame_skip` кадр, вплоть по `frames_to_write`-ый кадр, сохраняет как "results/`name`".

Так же рисует аналитическое решение `analitic(x,t)`, если таково передано.

TODO: Fix doc
"""
function make_gif(u::Matrix, Xₙ::Vector, Tₘ::Vector,
                  ϕ_l::Vector = missings(2), ϕ_r::Vector = missings(2),
                  f1::Vector = missings(2), f2::Vector = missings(2),
                  analitical = nothing;
                  frames_to_write::Int = -1, frame_skip::Int=-1,
                  name = "solution.gif", convert2mp4 = false)
    N,M = size(u) .-1
    if frames_to_write < 0
        frames_to_write = M;
    end
    if frame_skip < 0
        frame_skip = div(frames_to_write, 40) + 1
    end

    a = Animation()
    for m in 1:frame_skip:frames_to_write
        make_plot(u, Xₙ, Tₘ, m, ϕ_l, ϕ_r, f1, f2, analitical)
        frame(a)
    end

    path = joinpath(mkpath("results"), name)
    if convert2mp4
        run(`gif2mp4 $(path) $(replace(path, "gif" => "mp4")) \&`)
        path = replace(path, "gif" => "mp4")
    else
        g = gif(a, path, show_msg=false)
    end

    return path
end

@doc raw"""
    make_plot(u::Matrix, Xₙ::Vector, Tₘ::Vector, m::Int,
                   ϕ_l::Vector = missings(2), ϕ_r::Vector = missings(2),
                   f1::Vector = missings(2), f2::Vector = missings(2),
                   analitical = nothing)

Рисует `m`-ый кадр решения `u`. `Xₙ, Tₘ` — сетки.
`ϕ_l, ϕ_r` — вырожденные решения.
`f1, f2` — сгенерированные априорные данные.
`analitical` — или функция или сеточные значения аналитического решения.
"""
function make_plot(u::Matrix, Xₙ::Vector, Tₘ::Vector, m::Int,
                   ϕ_l::Vector = missings(2), ϕ_r::Vector = missings(2),
                   f1::Vector = missings(2), f2::Vector = missings(2),
                   analitical = nothing)

    yl = extrema(u[:,:]).*1.05;

    # График, оси, подписи осей и пр.
    pl = plot(size=(800, 600), xlabel="x", ylabel="u(x)", ylims=yl)

    # Начальное условие
    plot!(Xₙ, u[:,1], line=:dash, label="u_inital")

    # Найденное решение
    plot!(Xₙ, u[:,m], label="u(x,t)", color=:blue)

    # Точки сетки на найденной функции и их проекция на ось Х
    scatter!(Xₙ, u[:,m], color=:blue, label="", markersize=3)
    scatter!(Xₙ, [0 for x in Xₙ], color=:black, label="", markersize=2)

    # Надпись слева внизу с текущим временем
    annotate!(0.0, 0.9*first(yl), Plots.text(@sprintf("t = %.2f",Tₘ[m]), 16, :left ))

    check(x::Vector) = !any(ismissing.(x))
    if ( check(ϕ_l) ) && ( check(ϕ_r)) && ( check(f1) ) && ( check(f2) )
        ϕ = Φ(ϕ_l, ϕ_r, length(Xₙ)-1);
        plot!(Xₙ, ϕ_l, label=L"\phi_l", color=:darkgoldenrod)
        plot!(Xₙ, ϕ_r, label=L"\phi_r", color=:darkgoldenrod)
        plot!(Xₙ, ϕ, label=L"\widetilde{\Phi}", color=:gold)

        # Плавающая вслед за пунктиром надпись
        annotate!(0.0, f2[m] + 0.5, Plots.text(@sprintf("f2(t) = %.2f",f2[m]), 16, :left ))
        plot!([0, f1[m]], [f2[m], f2[m]], line=:dash, color=:black, label="")
        # красная точка слева, около подписи f2
        scatter!( [0], [f2[m]], color=:red, label="")

        # Плавающая вслед за пунктиром надпись
        annotate!(f1[m] + 0.01, 0.95 * first(yl), Plots.text(@sprintf("f1(t) = %.2f",f1[m]), 16, :left ))
        plot!([f1[m], f1[m]], [yl[1], f2[m]], line=:dash, color=:black, label="")
        # красная точка, нормировочный коэффициент, чтобы влезло в кадр
        scatter!( [f1[m]], [yl[1]*0.99], color=:red, label="")

        # красная точка на пересечении пунктиров априорной информации
        scatter!( [f1[m]], [f2[m]], color=:red, label="")
    end

    # Аналитическое решение
    if analitical != nothing
        if typeof( analitical ) <: Function
            plot!(Xₙ, analitical.(Xₙ, Tₘ[m]), label="analitical(x,t)", color=:green, linewidth = 5, alpha=0.3)
        elseif typeof( analitical) <: Matrix
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

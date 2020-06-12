"""
    module NonLinearReactionAdvectionDiffusionWithFrontData

Пакет реализует исследование в рамках гранта РФФИ #20-31-70016 "Численные методы
решения обратных задач для нелинейных сингулярно возмущённых уравнений типа
реакция-диффузия-адвекция с данными о положении фронта реакции".

Автор исходного кода пакета:
- [Андрей Борзунов](mailto:aborzunov@physics.msu.ru), Кафедра математики физического факультета МГУ им. Ломоносова.

Исполнители гранта:
 - Мельникова Алина Александровна
 - Левашова Наталия Тимуровна
 - Лукьяненко Дмитрий Витальевич (Руководитель)
 - Быцюра Светлана Владимировна
 - Борзунов Андрей Анатольевич
 - Исаев Темур Фуркатович
 - Аргун Рауль Ларикович
 - Горбачев Александр Викторович

Документация: https://github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl
"""
module NonLinearReactionAdvectionDiffusionWithFrontData

using ForwardDiff
using LinearAlgebra
using Dierckx
using LaTeXStrings
using Plots
using Printf
using Missings
using ProgressMeter

export u_init, solve;
export make_plot, make_gif, make_minimzation_gif;
export phidetermination, Φ;
export apply_on_dynamic_mesh;
export f1, f2;
export generate_obs_data;
export delta, solve_adjoint;
export J, J_q, minimize;

ENV["GKS_ENCODING"] = "utf-8"

include("utils.jl")
include("direct.jl")
include("degenerated.jl")
include("initial_guess.jl")
include("adjoint.jl")
include("functional.jl")
include("plotting.jl")

end # module

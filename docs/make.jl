using Documenter, NonLinearReactionAdvectionDiffusionWithFrontData
using ForwardDiff, LaTeXStrings, Plots

# Autogenerate documentation markdown for for all scripts
# in "examples/" package subfolder
# based on https://github.com/Evizero/Augmentor.jl
include("generatemd.jl")

# workaround GR warnings
# activates null device as output for GR
ENV["GKSwstype"] = "100"
GenerateMD.weave(overwrite=true)

DocMeta.setdocmeta!( NonLinearReactionAdvectionDiffusionWithFrontData, :DocTestSetup, :(using NonLinearReactionAdvectionDiffusionWithFrontData); recursive=true)
makedocs(
    modules=[NonLinearReactionAdvectionDiffusionWithFrontData],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages= Any[
        "Главная" => "index.md",
        "Прямая задача" => Any["Прямая задача" => "direct/direct.md",
                               "Генерирование априорной информации" => "direct/apriordata.md",
                               "Пример" => "generated/example_direct.md",
                               "Проверка на модельном решении" => "direct/direct_check.md",
                              ],
        "Сопряженная задача" => Any["Сопряженная задача" => "adjoint/adjoint.md",
                                    "Пример" => "generated/example_adjoint.md",
                                    "Проверка корректности решения сопряженной задачи" => "adjoint/adjoint_check.md",
                                   ],
        "Функционал, Градиент, Минимизация" => Any["Функционал" => "functional/functional.md",
                                                   "Пример" => "generated/example_functional.md",
                                                  ],
        "Справочник" => "reference.md",
    ],
    repo="https://github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl/blob/{commit}{path}#L{line}",
    sitename="NonLinearReactionAdvectionDiffusionWithFrontData.jl",
    authors="Andrey Borzunov",
)

deploydocs(;
    repo="github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl",
)

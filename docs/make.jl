using Documenter;
using NonLinearReactionAdvectionDiffusionWithFrontData;
using Literate;

using ForwardDiff, LaTeXStrings, Plots

# Перейдем в каталог текущего файла
cd(@__DIR__)
# activates null device as output for GR
ENV["GKSwstype"] = "100"
istravis = in("TRAVIS", keys(ENV))

# Autogenerate documentation markdown for for all scripts
# in "examples/" package subfolder
# based on https://github.com/Evizero/Augmentor.jl
#include("generatemd.jl")

# Функция реализующая предобработку
function replace_includes(str)

    included = [
                "example_direct.jl",
                "example_direct_dparams.jl",
                "example_direct_nonuniform.jl",
                "example_direct_nonuniform_dparams.jl",
                ]

    # Пусть к каталогу с примерами
    path = dirname(dirname(pathof(NonLinearReactionAdvectionDiffusionWithFrontData)))*"/examples/"

    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end

    return str
end

Literate.markdown("src/examples/de_direct.jl", "src/generated/"; name = "docexample_direct", preprocess = replace_includes, documenter = true)
#Literate.markdown("src/examples/de_direct_dparams.jl", "src/generated/"; name = "docexample_direct_dparams", preprocess = replace_includes, documenter = true)
#Literate.markdown("src/examples/de_direct_nonuniform.jl", "src/generated/"; name = "docexample_direct_nonuniform", preprocess = replace_includes, documenter = true)
#Literate.markdown("src/examples/de_direct_nonuniform_dparams.jl", "src/generated/"; name = "docexample_direct_nonuniform_dparams", preprocess = replace_includes, documenter = true)

DocMeta.setdocmeta!( NonLinearReactionAdvectionDiffusionWithFrontData, :DocTestSetup, :(using NonLinearReactionAdvectionDiffusionWithFrontData); recursive=true)
makedocs(
    modules=[NonLinearReactionAdvectionDiffusionWithFrontData],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages= Any[
        "index.md",
        "Прямая задача" => Any["direct/direct.md",
                               "direct/experimental_data.md",
                               "generated/docexample_direct.md",
                               #"generated/docexample_direct_dparams.md",
                               #"generated/docexample_direct_nonuniform.md",
                               #"generated/docexample_direct_nonuniform_dparams.md",
                               "direct/direct_check.md",
                              ],
        "Сопряженная задача" => Any["Сопряженная задача" => "adjoint/adjoint.md",
                                    "Проверка корректности решения сопряженной задачи" => "adjoint/adjoint_check.md",
                                   ],
        "Функционал, Градиент, Минимизация" => Any["Функционал" => "functional/functional.md",
                                                   #"Пример" => "generated/example_functional.md",
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

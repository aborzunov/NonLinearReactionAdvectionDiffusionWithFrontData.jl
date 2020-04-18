using Documenter;
using NonLinearReactionAdvectionDiffusionWithFrontData;
using Literate;

using ForwardDiff, LaTeXStrings, Plots

# Перейдем в каталог текущего файла
cd(@__DIR__)

# Autogenerate documentation markdown for for all scripts
# in "examples/" package subfolder
# based on https://github.com/Evizero/Augmentor.jl
#include("generatemd.jl")

# Функция реализующая предобработку
function replace_includes(str)

    included = [
                "example_direct.jl",
                "example_direct_nonuniform.jl",
                ]

    # Пусть к каталогу с примерами
    path = dirname(dirname(pathof(NonLinearReactionAdvectionDiffusionWithFrontData)))*"/examples/"

    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end

    return str
end

Literate.markdown("src/examples/doc_example_direct.jl", "src/examples/"; name = "docexample_direct", preprocess = replace_includes, documenter = true)

# workaround GR warnings
# activates null device as output for GR
ENV["GKSwstype"] = "100"
#istravis = in("TRAVIS", keys(ENV))
#istravis || GenerateMD.weave(overwrite=true) # Do not generate mds from examples

DocMeta.setdocmeta!( NonLinearReactionAdvectionDiffusionWithFrontData, :DocTestSetup, :(using NonLinearReactionAdvectionDiffusionWithFrontData); recursive=true)
makedocs(
    modules=[NonLinearReactionAdvectionDiffusionWithFrontData],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages= Any[
        "Главная" => "index.md",
        "Прямая задача" => Any["Прямая задача" => "direct/direct.md",
                               "Генерирование априорной информации" => "direct/apriordata.md",
                               "Пример на статической сетке" => "examples/docexample_direct.md",
                               "Проверка на модельном решении" => "direct/direct_check.md",
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

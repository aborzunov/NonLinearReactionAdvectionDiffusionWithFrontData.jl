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

# Функция реализующая предобработку примеров
function replace_includes_examples(str)

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

# Функция реализующая предобработку тестов
function replace_includes_test(str)

    included = [
                "direct_check.jl",
                ]

    # Пусть к каталогу с примерами
    path = dirname(dirname(pathof(NonLinearReactionAdvectionDiffusionWithFrontData)))*"/test/"

    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end

    return str
end

Literate.markdown("src/examples/de_direct.jl", "src/generated/"; name = "docexample_direct", preprocess = replace_includes_examples, documenter = true)
Literate.markdown("src/examples/dt_direct.jl", "src/generated/"; name = "doctest_direct", preprocess = replace_includes_test, documenter = true)

# Нам нужно, чтобы Literate сделал предобработку, но это возможно только
# для jl скриптов, а мы хотим оставить весь текст внутри md страницы
# документации. А подключать `doctest_direct.md` в качестве отдельной
# страницы не хотим. Поэтому сделаем append сгененированной md страницы
# к написанной.
content1 = read(dirname(dirname(pathof(NonLinearReactionAdvectionDiffusionWithFrontData)))*"/docs/src/direct/direct_check.md", String)
content2 = read(dirname(dirname(pathof(NonLinearReactionAdvectionDiffusionWithFrontData)))*"/docs/src/generated/doctest_direct.md", String)

io = open(dirname(dirname(pathof(NonLinearReactionAdvectionDiffusionWithFrontData)))*"/docs/src/direct/direct_check.md", "w")
print(io, content1 * content2)
close(io)


DocMeta.setdocmeta!( NonLinearReactionAdvectionDiffusionWithFrontData, :DocTestSetup, :(using NonLinearReactionAdvectionDiffusionWithFrontData); recursive=true)
makedocs(
    modules=[NonLinearReactionAdvectionDiffusionWithFrontData],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages= Any[
        "index.md",
        "Прямая задача" => Any["direct/direct.md",
                               "direct/experimental_data.md",
                               "generated/docexample_direct.md",
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

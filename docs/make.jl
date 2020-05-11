using Documenter;
using NonLinearReactionAdvectionDiffusionWithFrontData;
using Literate;

using ForwardDiff, LaTeXStrings, Plots

# Перейдем в каталог текущего файла
cd(@__DIR__)
# activates null device as output for GR
ENV["GKSwstype"] = "100"
istravis = in("TRAVIS", keys(ENV))


# package prefix
prefix =
dirname(dirname(pathof(NonLinearReactionAdvectionDiffusionWithFrontData)));

@info "\tGlobal package prefix: `$(prefix)`"

# path to folder, where generated files should be placed
gh_path = joinpath(prefix, "docs", "src", "generated", "helpers");

@info "\tSetting up \"docs/src/generated/helpers\" folder" gh_path
mkpath(gh_path);

# Функция реализующая предобработку тестов
# Ключевое отличие этой функции, от приведенной в примере Literate.jl
# в перечислении файлов вместе с их каталогом, т.к. у нас они берутся из двух
# разных каталогов. Городить что-то более сложное лень, поэтому здесь и в
# `incldue` примеров в документации, просто укажем папку откуда их берём.
function replace_includes(str)

    # Перечислим все файлы, содержимое которых эта функция может подставить
    included = [
                "examples/example_direct.jl",
                "examples/example_direct_dparams.jl",
                "examples/example_direct_nonuniform.jl",
                "examples/example_direct_nonuniform_dparams.jl",
                "test/direct_check.jl",
                "test/adjoint_check.jl",
                "examples/example_adjoint.jl",
                "examples/example_adjoint_nonuniform.jl",
                ]

    # `prefix` -- путь к нашему пакету
    # Все файлы в прошлом массиве должны быть определены через пути,
    # относительно prefix
    for ex in included
        # Начинаем чтение, только когда `ex` действительно инклюдится внутри str
        if occursin("include(\"$(ex)\")", str)
            content = read(joinpath(prefix, ex), String)
            str = replace(str, "include(\"$(ex)\")" => content)
        end
    end
return str
end

@info "\tGenerating sripts from `examples/` folder"
Literate.markdown("src/direct/direct_examples.jl",
                  "src/generated/";
                  name = "docexample_direct",
                  preprocess = replace_includes, documenter = true)
Literate.markdown("src/adjoint/adjoint_examples.jl",
                  "src/generated/";
                  name = "docexample_adjoint",
                  preprocess = replace_includes, documenter = true)

@info "\tGenerating scripts from `tests/` folder"
Literate.markdown("src/direct/check/dt_direct.jl",
                  "src/generated/helpers";
                  name = "doctest_direct",
                  preprocess = replace_includes, documenter = true)
Literate.markdown("src/adjoint/check/dt_adjoint.jl",
                  "src/generated/helpers";
                  name = "doctest_adjoint",
                  preprocess = replace_includes, documenter = true)

@info "\tComposing final check .md files"
# Нам нужно, чтобы Literate сделал предобработку, но это возможно только
# для jl скриптов, а мы хотим оставить весь текст с фомулами внутри md
# страницы документации. А создавать отдельную страницу для формул и отдельную
# страницу для тестов не хотим. Поэтому сделаем append сгененированных md
# страниц содержащих md unit-тестов к написанным страницами md с формулами.
#
# Пути к md страницам формул
# docs source file prefix
# -----------------------------------------------------------------------------
doc_prefix = "docs/src/";
generated_prefix = "docs/src/generated/helpers";
output_prefix = "docs/src/generated"
#
formulas_pages = ["direct/check/direct_check.md"; "adjoint/check/adjoint_check.md"]
generated_pages = ["doctest_direct.md"; "doctest_adjoint.md"]
output_pages = ["final_direct_check.md"; "final_adjoint_check.md"]
#
fp = map( x-> joinpath(prefix, doc_prefix, x), formulas_pages);
gp = map( x-> joinpath(prefix, generated_prefix, x), generated_pages);
op = map( x-> joinpath(prefix, output_prefix, x), output_pages);


for (f, g, o) in zip(fp, gp, op)

    content1 = read(f, String)
    content2 = read(g, String)

    io = open(o, create=true, truncate=true, write=true)
    print(io, content1 * "\n" * content2)
    close(io)

    @info "writing composed file to $(o)"
end
@info "\tLiterate stage has finished"
# -----------------------------------------------------------------------------


DocMeta.setdocmeta!( NonLinearReactionAdvectionDiffusionWithFrontData,
                    :DocTestSetup,
                    :(using NonLinearReactionAdvectionDiffusionWithFrontData);
                    recursive=true)

makedocs(
    modules=[NonLinearReactionAdvectionDiffusionWithFrontData],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages= Any[
        "index.md",
        "Прямая задача" => Any["direct/direct.md",
                               "direct/experimental_data.md",
                               "generated/docexample_direct.md",
                               "generated/final_direct_check.md",
                              ],
        "Сопряженная задача" => Any["adjoint/adjoint.md",
                                    "generated/final_adjoint_check.md",
                                    "generated/docexample_adjoint.md",
                                   ],
        "Обратная задача" => Any["functional/functional.md",
                                                   #"Пример" => "generated/example_functional.md",
                                                  ],
        "Справочник" => "reference.md",
    ],
    sitename    = "NonLinearReactionAdvectionDiffusionWithFrontData.jl",
    authors     = "Andrey Borzunov",
    clean       = true
)

deploydocs(;
    repo="github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl",
)

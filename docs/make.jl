using Documenter;
using NonLinearReactionAdvectionDiffusionWithFrontData;
using Literate;

using ForwardDiff, LaTeXStrings, Plots

# Перейдем в каталог текущего файла
cd(@__DIR__)
# activates null device as output for GR
ENV["GKSwstype"] = "100"
ENV["GKS_ENCODING"] = "utf-8"

# Режим сборки документации:
#   "ci"   (по умолчанию) — все примеры как статический код, без исполнения
#   "full" — примеры исполняются, генерируются графики и GIF
const FULL_BUILD = get(ENV, "DOCS_MODE", "ci") == "full"
@info "\tDocs build mode: $(FULL_BUILD ? "full" : "ci")"


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
                "examples/example_initial_guess.jl",
                "examples/example_functional.jl",
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

@info "\tGenerating md for numerical expirements"
# documenter=false всегда: это архивные расчёты (Nx=15000, S=35000),
# предназначены для ручного запуска через reports/run.jl
Literate.markdown(joinpath(prefix, "reports/experiments/data_generation.jl"),
                  "src/generated/";
                  name = "data_generation",
                  documenter = false)
Literate.markdown(joinpath(prefix, "reports/experiments/same_params.jl"),
                  "src/generated/";
                  name = "same_params",
                  documenter = false)

@info "\tGenerating md for functional section"
Literate.markdown("src/functional/numerical_expirements.jl",
                  "src/generated/";
                  name = "example_functional",
                  preprocess = replace_includes, documenter = FULL_BUILD)

@info "\tGenerating sripts from `examples/` folder"
Literate.markdown("src/direct/direct_examples.jl",
                  "src/generated/";
                  name = "docexample_direct",
                  preprocess = replace_includes, documenter = FULL_BUILD)
Literate.markdown("src/adjoint/adjoint_examples.jl",
                  "src/generated/";
                  name = "docexample_adjoint",
                  preprocess = replace_includes, documenter = FULL_BUILD)
Literate.markdown("src/asymptotics/initial_guess_example.jl",
                  "src/generated/";
                  name = "docexample_initial_guess",
                  preprocess = replace_includes, documenter = FULL_BUILD)

@info "\tGenerating scripts from `tests/` folder"
Literate.markdown("src/direct/check/dt_direct.jl",
                  "src/generated/helpers";
                  name = "doctest_direct",
                  preprocess = replace_includes, documenter = FULL_BUILD)
Literate.markdown("src/adjoint/check/dt_adjoint.jl",
                  "src/generated/helpers";
                  name = "doctest_adjoint",
                  preprocess = replace_includes, documenter = FULL_BUILD)

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
formulas_pages = [
                  "direct/check/direct_check.md";
                  "adjoint/check/adjoint_check.md"
                 ];
generated_pages = [
                   "doctest_direct.md";
                   "doctest_adjoint.md"
                  ];
output_pages = [
                "final_direct_check.md";
                "final_adjoint_check.md"
               ];
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
    warnonly = [:example_block, :autodocs_block],
    pages= Any[
        "index.md",
        "Прямая задача" =>
            Any[
                    "direct/direct.md",
                    "direct/experimental_data.md",
                    "generated/docexample_direct.md",
                    "generated/final_direct_check.md",
                ],

        "Асимптотические методы" =>
            Any[
                    "asymptotics/initial_guess.md",
                    "generated/docexample_initial_guess.md",
               ],

        "Сопряженная задача" =>
            Any[
                    "adjoint/adjoint.md",
                    "generated/final_adjoint_check.md",
                    "generated/docexample_adjoint.md",
               ],

        "Обратная задача" =>
            Any[
                    "functional/functional.md",
                    "generated/example_functional.md",
               ],

        "Эксперименты" =>
            Any[
                    "generated/data_generation.md",
                    "generated/same_params.md",
               ],

        "Методология тестирования" =>
            Any[
                    "testing/testing_mms.md",
                    "testing/testing_gradient.md",
                    "testing/testing_adjoint.md",
                    "testing/testing_convergence.md",
               ],

        "reference.md",
    ],
    sitename    = "NonLinearReactionAdvectionDiffusionWithFrontData.jl",
    authors     = "Andrey Borzunov",
    clean       = true
)

deploydocs(;
    repo="github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl",
    devbranch = "master",
    push_preview = true,
)

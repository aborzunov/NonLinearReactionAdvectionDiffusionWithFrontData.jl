using Documenter, NonLinearReactionAdvectionDiffusionWithFrontData

DocMeta.setdocmeta!( NonLinearReactionAdvectionDiffusionWithFrontData, :DocTestSetup, :(using NonLinearReactionAdvectionDiffusionWithFrontData); recursive=true)
makedocs(
    modules=[NonLinearReactionAdvectionDiffusionWithFrontData],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages= Any[
        "Главная" => "index.md",
        "Прямая заадча" => "direct.md",
        "Проверка на модельном решении" => "direct_check.md",
        "Генерирование априорной информации" => "apriordata.md",
        "Сопряженная задача" => "adjoint.md",
        "Проверка корректности решения сопряженной задачи" => "adjoint_check.md",
        "Справочник" => "reference.md",
    ],
    repo="https://github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl/blob/{commit}{path}#L{line}",
    sitename="NonLinearReactionAdvectionDiffusionWithFrontData.jl",
    authors="Andrey Borzunov",
)

deploydocs(;
    repo="github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl",
)

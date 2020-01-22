using Documenter, NonLinearReactionAdvectionDiffusionWithFrontData

DocMeta.setdocmeta!( NonLinearReactionAdvectionDiffusionWithFrontData, :DocTestSetup, :(using NonLinearReactionAdvectionDiffusionWithFrontData); recursive=true)
makedocs(;
    modules=[NonLinearReactionAdvectionDiffusionWithFrontData],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=[
        "Главная" => "index.md",
        "Проверка на модельном решении" => "model.md",
        "Вырожденные решения" => "degenerated.md",
        "Сопряженная задача" => "adjoint.md",
        "Справочник" => "reference.md",
    ],
    repo="https://github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl/blob/{commit}{path}#L{line}",
    sitename="NonLinearReactionAdvectionDiffusionWithFrontData.jl",
    authors="Andrey Borzunov",
)

deploydocs(;
    repo="github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl",
)

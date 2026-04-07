#!/usr/bin/env julia
#
# Запуск численных экспериментов (тяжёлые расчёты).
#
# Использование:
#   julia --project=reports reports/run.jl                    # все эксперименты
#   julia --project=reports reports/run.jl data_generation    # один конкретный
#   julia --project=reports reports/run.jl same_params        # другой конкретный
#
# Перед первым запуском:
#   julia --project=reports -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

using Dates

const EXPERIMENTS_DIR = joinpath(@__DIR__, "experiments")
const OUTPUT_DIR = joinpath(@__DIR__, "output")

# Список доступных экспериментов (имя => файл)
const AVAILABLE = Dict(
    "data_generation" => "data_generation.jl",
    "same_params"     => "same_params.jl",
    "noised"          => "noised.jl",
    "non_noised"      => "non_noised.jl",
    "perfect_example" => "perfect_example.jl",
    "qualitative"     => "qualitative.jl",
)

function run_experiment(name::String)
    if !haskey(AVAILABLE, name)
        error("Неизвестный эксперимент: \"$name\". Доступные: $(join(sort(collect(keys(AVAILABLE))), ", "))")
    end

    script = joinpath(EXPERIMENTS_DIR, AVAILABLE[name])
    stamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    outdir = joinpath(OUTPUT_DIR, "$(name)_$(stamp)")
    mkpath(outdir)

    @info "Запуск эксперимента \"$name\"" script outdir
    cd(outdir) do
        include(script)
    end
    @info "Эксперимент \"$name\" завершён. Результаты в: $outdir"
end

# --- main ---
if isempty(ARGS)
    @info "Запуск всех экспериментов"
    for name in sort(collect(keys(AVAILABLE)))
        run_experiment(name)
    end
else
    for name in ARGS
        run_experiment(name)
    end
end

```@meta
DocTestSetup = quote
    using NonLinearReactionAdvectionDiffusionWithFrontData
    using Plots
end
```

# Главная

## Постановка задачи

Рассмотрим прямую задачу для сингулярно возмущенного уравнения типа Бюргерса:
```math
\left\{
\begin{aligned}
    &\varepsilon\frac{\partial^2 u}{\partial x^2} -
    \frac{\partial u}{\partial t} = -u \frac{\partial u}{\partial x} +
    q(x)\,u, \quad x \in (0,1), \quad t \in (0,T], \\
    &u(0,t) = u_l(t), \quad u(1,t) = u_r(t), \quad t \in (0,T], \\
    &u(x,t) = u_i(x), \qquad x \in [0,1], t=0.
\end{aligned}
\right.
```

Решение этой задачи имеет движущийся слой, положение которого во времени
описывает ``x = x_{t.p}(t)``.

Обратная задача состоит в определении коэффициента линейного усиления ``q(x)``,
``x \in [0,1]``, по известной дополнительной информации о положении переходного
слоя и значения функции на переходном слое:
```math
x_{t.p} = f_1(t), \qquad u(x_{t.p}(t),t) = f_2(t), \qquad t \in [0, T].
```

## Содержание

### Прямая задача

```@contents
Pages = [
    "direct/direct.md",
    "direct/experimental_data.md",
    "generated/docexample_direct.md",
    "generated/final_direct_check.md",
]
Depth = 1
```

### Сопряженная задача

```@contents
Pages = [
    "adjoint/adjoint.md",
    "generated/final_adjoint_check.md",
    "generated/docexample_adjoint.md",
]
Depth = 1
```

### Обратная задача


```@contents
Pages = [
    "functional/functional.md",
    "functional/numerical_expirements.md",
]
Depth = 1
```

## Локальная документация

Документация на `gh-pages` собирается в режиме `ci` — примеры показываются
как статический код без графиков и анимаций.

Чтобы сгенерировать полную документацию с графиками и GIF локально:
```bash
git clone https://github.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl
cd NonLinearReactionAdvectionDiffusionWithFrontData.jl
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
DOCS_MODE=full julia --project=docs docs/make.jl
python3 -m http.server 8000 --directory docs/build
```
После, откройте [http://localhost:8000](http://localhost:8000) в браузере.

Тяжёлые численные эксперименты (для статьи) вынесены в отдельный каталог `reports/`
и запускаются отдельно:
```bash
julia --project=reports -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=reports reports/run.jl
```

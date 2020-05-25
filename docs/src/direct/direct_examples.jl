# # Примеры
#
# На этой странице мы делаем `include` примеров,
# которые лежат в корневом каталоге пакета в папке `"examples/"`.
#
# [Literate.jl](https://github.com/fredrikekre/Literate.jl) включает
# текст этих примеров и в md файл.
#
#src `#src` это комментарии даже для Literate.jl
#src `#md` Строки Literate.jl включит только в md файлы.
#
#src Перед всматриванием в следующий блок, посмотри `docs/make.jl`.
#src Там, скормив этот скрипт Literate.jl мы получим страницу
#src `docexample_direct.md`. Именно её будет обрабатывать Documenter.jl
#
#md # ```@contents
#md # Pages = [ "docexample_direct.md" ]
#md # Depth = 5
#md # ```

include("examples/example_direct.jl")

# ### Визуализация

# Мы не можем строить большие анимации на стороне `Travis`-a.
# Если мы на CI, то будем рисовать только 10 кадров.
# Если мы генерируем документацию локально, то рисуем 80 кадров.
isTravis = in("Travis", keys(ENV))
ftw = isTravis ? range(1, stop = M+1, length=9) : [1; 2:div(M+1, 80):M+1];

## Запись gif одного только решения
make_gif(u, XX, Tₘ; name="solution_direct_ex1.gif", frames_to_write = ftw)

## Запись **только** mp4 вместе с вырожденными корнями и переходным слоем.
make_gif(u, XX, Tₘ, ϕl, ϕr, ϕ, f1_data, f2_data;
         name="solution_direct_ex11.mp4", frames_to_write = ftw)






include("examples/example_direct_dparams.jl")

# ### Визуализация
# Так выглядит решение для стандартных данных.
make_gif(u, XX, Tₘ, ϕl, ϕr, ϕ, f1_data, f2_data; name="solution_direct_ex2.gif", frames_to_write = ftw)






include("examples/example_direct_nonuniform.jl")

# ### Визуализация
# Так выглядит решение для динамической сетке
ftw = isTravis ? range(1, stop = M+1, length=15) : [1; 2:div(M+1, 80):M+1];
make_gif(u, XX, Tₘ, ϕl, ϕr, ϕ, f1_data, f2_data; name="solution_direct_ex3.gif", frames_to_write = ftw)





include("examples/example_direct_nonuniform_dparams.jl")

# ### Визуализация
# Так выглядит решение для стандартных данных динамической сетке.
ftw = [1; 2:div(M+1, 80):M+1];
make_gif(u, XX, Tₘ, ϕl, ϕr, ϕ, f1_data, f2_data; name="solution_direct_ex4.gif", frames_to_write = ftw)

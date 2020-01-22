# NonLinearReactionAdvectionDiffusionWithFrontData.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aborzunov.github.io/NonLinearReactionAdvectionDiffusionWithFrontData.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aborzunov.github.io/NonLinearReactionAdvectionDiffusionWithFrontData.jl/dev)
[![Build Status](https://travis-ci.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl.svg?branch=master)](https://travis-ci.com/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl)
[![Codecov](https://codecov.io/gh/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl)
[![Coveralls](https://coveralls.io/repos/github/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl/badge.svg?branch=master)](https://coveralls.io/github/aborzunov/NonLinearReactionAdvectionDiffusionWithFrontData.jl?branch=master)

# Численные методы решения обратных задач для нелинейных сингулярно возмущённых уравнений типа реакция-диффузия-адвекция с данными о положении фронта реакции


## Постановка задачи

![statement](http://www.sciweavers.org/tex2img.php?eq=%5Cleft%5C%7B%0A%5Cbegin%7Baligned%7D%0A%20%20%20%20%26%5Cvarepsilon%5Cfrac%7B%5Cpartial%5E2%20u%7D%7B%5Cpartial%20x%5E2%7D%20-%20%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20t%7D%20%3D%20-u%20%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20x%7D%20%2B%20%20q%28x%29%5C%2Cu%2C%20%5Cquad%20x%20%5Cin%20%280%2C1%29%2C%20%5Cquad%20t%20%5Cin%20%280%2CT%29%2C%20%5C%5C%0A%20%20%20%20%26u%280%2Ct%29%20%3D%20u_%7Bleft%7D%28t%29%2C%20%5Cquad%20u%281%2Ct%29%20%3D%20u_%7Bright%7D%28t%29%2C%20%5Cquad%20t%20%5Cin%20%280%2CT%29%2C%20%5C%5C%0A%20%20%20%20%26u%28x%2Ct%29%20%3D%20u_%7Binit%7D%28x%29%2C%20%5Cqquad%20x%20%5Cin%20%5B0%2C1%5D.%0A%5Cend%7Baligned%7D%0A%5Cright.%0A&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

Решение этой задачи имеет движущийся слой, положение которого во времени описывается функцией `x = x_{t.p.}(t)`.

Обратная задача состоит в определении коэффициента линейного усиления `q(x)$, $x \in [0,1]`,
по известной дополнительной информации о положении переходного слоя и значении функции
(если необходимо, но желательно без дополнительной функции) на переходном слое:

![img](http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Bequation%2A%7D%0A%20%20%20%20x_%7Bt.p.%7D%28t%29%20%3D%20f_%7B1%7D%28t%29%2C%20%5Cquad%20u%5Cbig%28x_%7Bt.p.%7D%28t%29%2Ct%5Cbig%29%20%3D%20f_%7B2%7D%28t%29%2C%20%5Cquad%20t%20%5Cin%20%5B0%2CT%5D.%0A%5Cend%7Bequation%7D%0A&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

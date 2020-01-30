module NonLinearReactionAdvectionDiffusionWithFrontData

using ForwardDiff
using LinearAlgebra
using Dierckx
using LaTeXStrings
using Plots
using Printf
using Missings

export u_init, directRP, solve;
export make_plot, make_gif;
export phidetermination, Φ;
export delta, adjointRP, solve_adjoint;
export J, J_q;


include("utils.jl")
include("direct.jl")
include("degenerated.jl")
include("adjoint.jl")
include("functional.jl")
include("plotting.jl")

end # module

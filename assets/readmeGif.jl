# automation script for generating READEME gif
# Runs example script for direct problem
# Generates gif
# Saves it in "assets/" folder
# Later, you should manually upload it via Git-LFS


#package prefix
prefix =
dirname(dirname(pathof(NonLinearReactionAdvectionDiffusionWithFrontData)));
assets_path =  joinpath(prefix, "assets/");
gif_path = joinpath(assets_path, "readme_direct.gif");

include( joinpath(prefix, "examples/example_direct_nonuniform_dparams.jl"))

ENV["GKS_ENCODING"] = "utf-8"
@info "Writing README gif" gif_path
make_gif(u, XX, Tₘ, ϕl, ϕr, ϕ, f1_data, f2_data;
         name=gif_path, title = "Прямая задача и экспериментальные данные")

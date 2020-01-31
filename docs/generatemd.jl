"""
    module GenerateMD

Generates `Documenter` compatible `.md`  files using `Weave` from
`.jl` scripts in "examples/" folder of `NonLinearReactionAdvectionDiffusionWithFrontData`
package. Based on [ExampleWeaver](https://github.com/Evizero/Augmentor.jl).
"""
module GenerateMD

##########################################################################################
const EXAMPLES_DIR  = abspath(joinpath(@__DIR__, "..", "examples"))
const GENERATED_DIR = abspath(joinpath(@__DIR__, "src", "generated"))
function _listfiles(dir, ext, fullpath=false)
    fnames = filter(fname->splitext(fname)[2]==ext, readdir(dir))
    fullpath ? map(fname->joinpath(dir, fname), fnames) : fnames
end
listexamples(fullpath=false) = _listfiles(EXAMPLES_DIR, ".jl", fullpath)
listmarkdown(fullpath=false) = _listfiles(GENERATED_DIR, ".md", fullpath)
##########################################################################################

function weave_markdown(scriptname; overwrite=false)
    splitext(scriptname)[2] == ".jl" || return
    name = splitext(scriptname)[1]
    # define all required paths
    scriptpath = joinpath(EXAMPLES_DIR, scriptname)
    processed_scriptpath = joinpath(GENERATED_DIR, name * ".jl")
    jmdpath = joinpath(GENERATED_DIR, name * ".jmd")
    mdpath = joinpath(GENERATED_DIR, name * ".md")
    # if markdown file already exists, only overwrite if requested
    if isfile(mdpath) && !overwrite
        @info "skipping markdown generation for \"$scriptname\" (file already exists)"
        return mdpath
    else
        @info "generating markdown \"$(name*".md")\" for \"$scriptname\""
        mkpath(GENERATED_DIR)
    end
    convert_doc(scriptpath, jmdpath)
    # posprocess the .jmd and save it as .md for documenter
    str_md = read(jmdpath, String)
    str_md = replace(str_md, "```julia" => "```@example $name")
    write(mdpath, str_md)
    # cleanup temporary files
    rm(jmdpath)
    # return path to final .md file
    mdpath
end

function weave(; overwrite=false)
    mds = String[];
    for scriptname in listexamples()
        md = weave_markdown(scriptname; overwrite=overwrite)
        push!(mds, md)
    end
    mds
end

end # module

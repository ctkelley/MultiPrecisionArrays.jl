push!(LOAD_PATH, "../src/")
using Documenter, MultiPrecisionArrays, DocumenterTools

struct LaTeXEquation
    content::String
end

function Base.show(io::IO, ::MIME"text/latex", x::LaTeXEquation)
    # Wrap in $$ for display math printing
    return print(io, "\$\$ " * x.content * " \$\$")
end

makedocs(
    sitename = "MultiPrecisionArrays.jl",
    authors = "C. T. Kelley",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = Any[
        "Home"=>"index.md",
        "Factorizations"=>Any["functions/mplu!.md", "functions/hlu!.md"],
    ],
)
deploydocs(repo = "github.com/ctkelley/MultiPrecisionArrays.jl.git")

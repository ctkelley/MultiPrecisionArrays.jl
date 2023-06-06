using Documenter, MultiPrecisionArrays, DocumenterTools
push!(LOAD_PATH,"../src/")
makedocs(sitename="MPArrays.jl",
authors="C. T. Kelley",
format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true"
           ),
pages = Any[
     "Home" => "index.md",
]
)
deploydocs(
     repo="github.com/ctkelley/MultiPrecisionArrays.jl.git"
)

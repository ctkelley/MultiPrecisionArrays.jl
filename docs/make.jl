push!(LOAD_PATH,"../src/")
using Documenter, MultiPrecisionArrays, DocumenterTools
makedocs(sitename="MultiPrecisionArrays.jl",
authors="C. T. Kelley",
format = Documenter.HTML(
               prettyurls = get(ENV, "CI", nothing) == "true"
           ),
pages = Any[
     "Home" => "index.md",
     "Factorizations" => Any[
     "functions/hlu!.md",
     ]
]
)
deploydocs(
     repo="github.com/ctkelley/MultiPrecisionArrays.jl.git"
)

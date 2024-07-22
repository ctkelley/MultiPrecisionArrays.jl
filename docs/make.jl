push!(LOAD_PATH, "../src/")
using Documenter, MultiPrecisionArrays, DocumenterTools, DocumenterCitations
import MultiPrecisionArrays.MPArray
import MultiPrecisionArrays.MPHArray
import MultiPrecisionArrays.MPGArray
import MultiPrecisionArrays.MPBArray
import MultiPrecisionArrays.MPBFact
import MultiPrecisionArrays.MPKFact
import MultiPrecisionArrays.MPLFact
import MultiPrecisionArrays.MPGEFact
import MultiPrecisionArrays.MPFact
import MultiPrecisionArrays.mpgeslir
import MultiPrecisionArrays.mpkrir


mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "MPArray.bib"),
#    style=:numeric  # default
    style=:authoryear # My way!
)


makedocs(
    sitename = "MultiPrecisionArrays.jl",
    authors = "C. T. Kelley",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = Any[
        "Home"=>"index.md",
        "Half Precision and Krylov-IR" => Any["Half_1.md",], 
        "More than you want to know" => 
           Any["Details/Termination.md", "Details/N2Work.md",
         "Details/Interprecision_1.md", "Details/Extended.md",],
        "MPArray Constructors" => 
           Any["functions/MPArray.md","functions/MPGArray.md",
               "functions/MPBArray.md",],
        "Factorizations"=>
          Any["functions/hlu!.md", "functions/mplu!.md", "functions/mplu.md", 
                 "functions/mpglu!.md", "functions/mpglu.md",
                 "functions/mpblu!.md", "functions/mpblu.md",],
        "Iteration Statistics"=>
          Any["Details/Stats.md",],
        "Solvers"=>
             Any["functions/mpgeslir.md", "functions/mpkrir.md"],
        "References" => ["References.md",],
],
; plugins=[bib]
)
deploydocs(repo = "github.com/ctkelley/MultiPrecisionArrays.jl.git")

push!(LOAD_PATH, "../src/")
using Documenter, MultiPrecisionArrays, DocumenterTools

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

makedocs(
    sitename = "MultiPrecisionArrays.jl",
    authors = "C. T. Kelley",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = Any[
        "Home"=>"index.md",
        "More than you want to know" => "Details.md",
        "Factorizations"=>Any["functions/MPArray.md", "functions/mplu!.md", "functions/hlu!.md"],
    ],
)
deploydocs(repo = "github.com/ctkelley/MultiPrecisionArrays.jl.git")

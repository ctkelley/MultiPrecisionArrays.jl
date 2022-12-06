module TestMP

using MPArrays
using LinearAlgebra
using SparseArrays
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using LaTeXStrings
using Printf
using PyPlot

export jacmp!
export basic2d!
export jbasic2d!
export basicnewt
export htest
export plot_its_funs
export nl_stats!
export ktst
export dstest
export Gmat


include("ntest.jl")
include("plot_its_funs.jl")
include("greens.jl")
include("../Examples/Gmat.jl")

end

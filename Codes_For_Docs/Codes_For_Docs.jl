module Codes_For_Docs

using MultiPrecisionArrays
using LinearAlgebra
using LinearAlgebra.BLAS
using BenchmarkTools

include("IR.jl")
include("HalfTime.jl")

export IR
export HalfTime

end

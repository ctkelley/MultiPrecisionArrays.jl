module Codes_For_Docs

using MultiPrecisionArrays
using MultiPrecisionArrays.Examples
using LinearAlgebra
using LinearAlgebra.BLAS
using BenchmarkTools
using Printf

include("IR.jl")
include("HalfTime.jl")
include("fprintTeX.jl")

export IR
export HalfTime

end

module MPArrays

using LaTeXStrings
using LinearAlgebra
using SparseArrays
using Printf

struct MPArray
   AD::Array{Float64,2}
   AS::Array{Float32,2}
end

struct MPFArray
   AD::Array{Float64,2}
   AFS::LU{Float32, Matrix{Float32}}
end

function MPArray(AD::Array{Float64,2})
   AS=Float32.(AD)
   MPB=MPArray(AD,AS)
   return MPB
end

function MPFArray(AD::Array{Float64,2})
   AS=Float32.(AD)
   ASF=lu!(AS)
   MPB=MPFArray(AD,ASF)
   return MPB
end

import Base.\
function \(AF::MPFArray, b)
xi = mpgesl(AF,b)
return xi
end

export MPArray
export MPFArray
export mpgesl
export testback

include("Solvers/mpgesl.jl")
include("testback.jl")

end

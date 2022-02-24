module MPArrays

using LaTeXStrings
using LinearAlgebra
using SparseArrays
using Printf

struct MPTest
  AH::Array
  AL::Array
  AF::Factorization
end

function MPTest(AD::Array{Float64,2})
   AL=Float32.(AD);
   AF=lu!(AL);
   MPTest(AD, AL, AF)
end

function MPTest(AD::Array{Float32,2})
   AL=Float16.(AD);
   AF=lu!(AL);
   MPTest(AD, AL, AF)
end

function MPTest!(MTF::MPTest, AD::Array{Float64,2})
   Tin=eltype(AD)
   Tout=eltype(MTF.AH)
   Tlow=eltype(MTF.AL)
#   println(Tin, "  ", Tout,"  ",Tlow)
   (Tin == Tout) || error("Types of AD and MPArray must agree")
   MTF.AH .= AD
   MTF.AL .= Tlow.(AD);
   ZF = lu!(MTF.AL);
   MPTest(MTF.AH, MTF.AL, ZF)
end

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

import Base.\
function \(AF::MPTest, b)
xi = mpgesl2(AF,b)
return xi
end

export MPArray
export MPFArray
export MPTest
export MPTest!
export mpgesl
export mpgesl2
export testback

include("Solvers/mpgesl.jl")
include("Solvers/mpgesl2.jl")
include("testback.jl")

end

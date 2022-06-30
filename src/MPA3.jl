module MPA3

using LinearAlgebra
using SparseArrays
using SIAMFANLEquations

struct MPTest3{TH<:AbstractFloat, TL<:AbstractFloat, TF<:Factorization}
  AH::Array{TH,2}
  AL::Array{TL,2}
  AF::TF
end

struct MPArray3{TH<:AbstractFloat, TL<:AbstractFloat}
   AH::Array{TH,2}
   AL::Array{TL,2}
end

function MPArray3(AH::Array{Float64,2})
AL=Float32.(AH)
MPA=MPArray3(AH,AL)
end

function mplu3!(MPA::MPArray3)
AH=MPA.AH
AL=MPA.AL
AF=lu!(AL)
MPF=MPTest3(AH, AL, AF)
return MPF
end

export MPArray3
export mplu3!

end


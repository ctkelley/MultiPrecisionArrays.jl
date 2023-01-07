module MPArrays

using LinearAlgebra 
using SparseArrays
using SIAMFANLEquations

#
# This is MPArrays.jl
# The package has data structures and algorithms to manage several
# variations on iterative refinement.
#

include("MPStructs/MPLight.jl")
include("MPStructs/MPHeavy.jl")

struct MPGArray
    AH::Array
    AH2::Array
    VH::Array
    AS::Array
end

struct MPGFact
  AH::Array
  AL::Array
  VH::Array
  AS::Array
  AF::Factorization
end


MPFact=Union{MPLFact, MPLEFact, MPHFact, MPGFact}

MPGMFact=Union{MPHFact, MPGFact}

function MPhatv(x, MPHF::MPGMFact)
atv=MPHF.AH*x
return atv
end

function MPhptv(x, MPHF::MPGMFact)
ptv = MPHF.AF\x
return ptv
end


MPIRArray=Union{MPArray,MPHArray}

function MPGArray(AH::Array{Float64,2}, TL=Float32, itg=5)
(ma, na)=size(AH)
T=eltype(AH)
AH2=copy(AH)
AS=TL.(AH)
VH=zeros(T, na, itg)
MPG=MPGArray(AH, AH2, VH, AS)
end

function MPGArray(AH::Array{Float32,2}, TL=Float16, itg=5)
(ma, na)=size(AH)
T=eltype(AH)
AH2=copy(AH)
AS=TL.(AH)
VH=zeros(T, na, itg)
MPG=MPGArray(AH, AH2, VH, AS)
end

function mpglu!(MPG::MPGArray)
AH=MPG.AH
VH=MPG.VH
TD=eltype(AH)
AH2=MPG.AH2
AS=MPG.AS
ASF=lu!(AS)
AH2 .= TD.(AS)
AF = LU(AH2, ASF.ipiv, ASF.info)
MPF=MPGFact(AH, AH2, VH, AS, AF)
end

import Base.eltype
function eltype(MP::MPArray)
TP=eltype(MP.AH)
return TP
end

import Base.\
function \(AF::MPLFact, b; verbose=false, reporting=false)
xi = mpgesl2(AF,b; verbose=verbose, reporting=reporting)
return xi
end

function \(AF::MPLEFact, b; verbose=false, reporting=false)
xi = mpgesl2(AF,b; verbose=verbose, reporting=reporting)
return xi
end

function \(AF::MPHFact, b; verbose=false, reporting=false)
xi = mpgesl2(AF,b; verbose=verbose, reporting=reporting)
return xi
end

function \(AF::MPGFact, b; verbose=false, reporting=false)
xi = mpgesl2(AF,b; verbose=verbose, reporting=reporting)
return xi
end

function promotelu!(AL, AH)
TH=eltype(AH)
ALF=lu!(AL)
AH .= TH.(AL)
AHF = LU(AH, ALF.ipiv, ALF.info)
return AHF
end

function promotelu(A, T=Float32)
AF=lu!(A)
AFHigh=LU(T.(A), AF.ipiv, AF.info)
return AFHigh
end

export mplu!
export mphlu!
export mpglu!
export mpqr!
export mpcholesky!
export MPArray
export MPEArray
export MPFArray
export MPHArray
export MPGArray
export MPLFact
export MPLEFact
export MPHFact
export MPGFact
export mpgesl2
export promotelu
export promotelu!
export MPhatv
export MPhptv

export MPGStats
export MPIRStats

include("Solvers/mpgesl2.jl")
include("Solvers/IRTriangle.jl")
include("MPAStats.jl")

module Examples
export Gmat

include("Examples/Gmat.jl")

end

end # module MPArrays

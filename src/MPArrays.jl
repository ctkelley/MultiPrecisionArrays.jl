module MPArrays

using LinearAlgebra 
using SparseArrays
using SIAMFANLEquations

struct MPTest
  AH::Array
  AL::Array
  AF::Factorization
end

struct MPHTest
  AH::Array
  AL::Array
  AS::Array
  AF::Factorization
end

struct MPGArray
    AH::Array
    AH2::Array
    VH::Array
    AS::Array
end

struct MPGTest
  AH::Array
  AL::Array
  VH::Array
  AS::Array
  AF::Factorization
end

MPFact=Union{MPTest, MPHTest, MPGTest}

MPHFact=Union{MPHTest, MPGTest}

function MPhatv(x, MPHF::MPHFact)
atv=MPHF.AH*x
return atv
end

function MPhptv(x, MPHF::MPHFact)
ptv = MPHF.AF\x
return ptv
end


struct MPArray
   AH::Array
   AL::Array
end

struct MPHArray
    AH::Array
    AH2::Array
    AS::Array
end

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


function MPHArray(AH::Array{Float64,2}, TL=Float32)
AH2=copy(AH)
AS=TL.(AH)
MPH=MPHArray(AH, AH2, AS)
end

function MPHArray(AH::Array{Float32,2}, TL=Float16)
AH2=copy(AH)
AS=TL.(AH)
MPH=MPHArray(AH, AH2, AS)
end

function mphlu!(MPH::MPHArray)
AH=MPH.AH
TD=eltype(AH)
AH2=MPH.AH2
AS=MPH.AS
ASF=lu!(AS)
AH2 .= TD.(AS)
AF = LU(AH2, ASF.ipiv, ASF.info)
MPF=MPHTest(AH, AH2, AS, AF)
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
MPF=MPGTest(AH, AH2, VH, AS, AF)
end

function MPArray(AH::Array{Float32,2})
AL=Float16.(AH)
MPA=MPArray(AH,AL)
end

function MPArray(AH::Array{Float64,2})
AL=Float32.(AH)
MPA=MPArray(AH,AL)
end

function mplu!(MPA::MPArray)
AH=MPA.AH
AL=MPA.AL
AF=lu!(AL)
MPF=MPTest(AH, AL, AF)
return MPF
end

function mpqr!(MPA::MPArray)
AH=MPA.AH
AL=MPA.AL
AF=qr!(AL)
MPF=MPTest(AH, AL, AF)
return MPF
end

function mpcholesky!(MPA::MPArray)
AH=MPA.AH
AL=MPA.AL
AF=cholesky!(AL)
MPF=MPTest(AH, AL, AF)
return MPF
end

import Base.eltype
function eltype(MP::MPArray)
TP=eltype(MP.AH)
return TP
end

import Base.\
function \(AF::MPTest, b)
xi = mpgesl2(AF,b)
return xi
end

function \(AF::MPHTest, b)
xi = mpgesl2(AF,b)
return xi
end

function \(AF::MPGTest, b)
xi = mpgesl2(AF,b)
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
export MPFArray
export MPHArray
export MPGArray
export MPTest
export MPHTest
export MPGTest
export mpgesl2
export promotelu
export MPhatv
export MPhptv

include("Solvers/mpgesl2.jl")

end

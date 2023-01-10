struct MPHArray
    AH::Array
    AH2::Array
    AS::Array
end

struct MPHFact
  AH::Array
  AL::Array
  AS::Array
  AF::Factorization
end

function MPHArray(AH::Array{Float64,2}; TL=Float32)
AH2=copy(AH)
AS=TL.(AH)
MPH=MPHArray(AH, AH2, AS)
end

function MPHArray(AH::Array{Float32,2}; TL=Float16)
AH2=copy(AH)
AS=TL.(AH)
MPH=MPHArray(AH, AH2, AS)
end

function mphlu!(MPH::MPHArray)
AH=MPH.AH
TD=eltype(AH)
AH2=MPH.AH2
AS=MPH.AS
#
# Factor in low precision
#
ASF=lu!(AS)
#
# Promote he low-precision lu
#
AH2 .= TD.(AS)
AF = LU(AH2, ASF.ipiv, ASF.info)
MPF=MPHFact(AH, AH2, AS, AF)
end


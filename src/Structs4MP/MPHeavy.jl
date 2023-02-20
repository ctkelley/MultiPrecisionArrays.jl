struct MPHArray
    AH::Array
    AStore::Array
    ALow::Array
end

struct MPHFact
  AH::Array
  AStore::Array
  ALow::Array
  AF::Factorization
end

#
# Heavy factorization for GMRES-IR
#
struct MPGHFact
  AH::Array
  AStore::Array
  ALow::Array
  AF::Factorization
end


function MPHArray(AH::Array{Float64,2}; TL=Float32)
AStore=copy(AH)
ALow=TL.(AH)
MPH=MPHArray(AH, AStore, ALow)
end

function MPHArray(AH::Array{Float32,2}; TL=Float16)
AStore=copy(AH)
ALow=TL.(AH)
MPH=MPHArray(AH, AStore, ALow)
end

function mphlu!(MPH::MPHArray)
AH=MPH.AH
TD=eltype(AH)
AStore=MPH.AStore
ALow=MPH.ALow
#
# Factor in low precision
#
ALowF=lu!(ALow)
#
# Promote the low-precision lu
#
AStore .= TD.(ALow)
AF = LU(AStore, ALowF.ipiv, ALowF.info)
MPF=MPHFact(AH, AStore, ALow, AF)
end

function mpglu!(MPH::MPHArray)
AH=MPH.AH
TD=eltype(AH)
AStore=MPH.AStore
ALow=MPH.ALow
#
# Factor in low precision
#
ALowF=lu!(ALow)
#
# Promote the low-precision lu
#
AStore .= TD.(ALow)
AF = LU(AStore, ALowF.ipiv, ALowF.info)
MPF=MPGHFact(AH, AStore, ALow, AF)
end




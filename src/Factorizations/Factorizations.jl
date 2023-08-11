# 
# A few factorizations
# 

function mphlu!(MPH::MPHArray)
AH=MPH.AH
TD=eltype(AH)
AStore=MPH.AStore
AL=MPH.AL
TL=eltype(AL)
(TL==Float16) ? ALF=hlu!(AL) : ALF=lu!(AL)
#
# Factor in low precision
#
#ALF=lu!(AL)
#
# Promote the low-precision lu
#
AStore .= TD.(AL)
AF = LU(AStore, ALF.ipiv, ALF.info)
MPF=MPHFact(AH, AStore, AL, AF)
end

function mpglu!(MPH::MPHArray)
AH=MPH.AH
TD=eltype(AH)
AStore=MPH.AStore
AL=MPH.AL
TL=eltype(AL)
(TL==Float16) ? ALF=hlu!(AL) : ALF=lu!(AL)
#
# Factor in low precision
#
#ALF=lu!(AL)
#
# Promote the low-precision lu
#
AStore .= TD.(AL)
AF = LU(AStore, ALF.ipiv, ALF.info)
MPF=MPGHFact(AH, AStore, AL, AF)
end

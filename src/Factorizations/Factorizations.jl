# 
# A few factorizations
# 

function mplu!(MPA::MPEArray)
AH=MPA.AH
AL=MPA.AL
TL=eltype(AL)
(TL == Float16) ? AF=hlu!(AL) : AF=lu!(AL)
#AF=lu!(AL)
MPF=MPLEFact(AH, AL, AF)
return MPF
end

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

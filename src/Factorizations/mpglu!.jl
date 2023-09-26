#
# Factor a MPErray and set it up for GMRES with \
#
function mpglu!(MPH::MPArray; basissize=10)
AH = MPH.AH
TD = eltype(AH)
res = MPH.residual
n=length(res)
AL = MPH.AL
    TL = eltype(AL)
    #
    # Factor in low precision
    #
    (TL == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
    #
VStore=zeros(TD, n, basissize)
KStore=kstore(n,"gmres")
MPF=MPGEFact(AH, AL, ALF, VStore, KStore, res, true)
return MPF
end

#
# Factor a heavy MPArray and set it up for GMRES with \
# If you want to use it with IR (why?) then set gmresok=false
#
function mpglu!(MPH::MPHArray; gmresok = true, basissize=10)
    AH = MPH.AH
    TD = eltype(AH)
    res = MPH.residual
    n=length(res)
    AStore = MPH.AStore
    AL = MPH.AL
    TL = eltype(AL)
    #
    # Factor in low precision
    #
    (TL == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
    #
    # Promote the low-precision lu
    #
    AStore .= TD.(AL)
    AF = LU(AStore, ALF.ipiv, ALF.info)
    if gmresok
        VStore=zeros(TD, n, basissize)
        KStore=kstore(n,"gmres")
        MPF = MPGHFact(AH, AL, AF, VStore, KStore, res, true)
    else
        MPF = MPHFact(AH, AL, AF, res, true)
    end
end

#
# Using a heavy MPArray with normal IR is insane. I use it only
# for CI.
#
function mphlu!(MPH::MPHArray)
    MPF = mpglu!(MPH; gmresok = false)
end

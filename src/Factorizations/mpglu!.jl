"""
mpglu!(MPH::MPArray; basissize=10)
Factor a MPArray and set it up for GMRES by allocating room
for Krylov vectors etc
"""
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

"""
mpglu(A::Array{TH,2}; TL=Float32, basissize=10) where TH <: Real

Combines the constructor of the multiprecision GMRES-ready array with the
factorization.
"""
function mpglu(A::Array{TH,2}; TL=Float32, basissize=10) where TH <: Real
#
# If the high precision matrix is single, the low precision must be half.
#
(TH == Float32) && (TL = Float16)
#
# Unless you tell me otherwise, onthefly is true if low precision is half
# and false if low precision is single.
#
MPA=MPArray(A; TL=TL, onthefly=true)
MPGF=mpglu!(MPA; basissize=basissize)
return MPGF
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

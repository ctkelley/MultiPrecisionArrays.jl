"""
mpglu!(MPGA::MPGArray)
Factor a MPGArray and set it up for GMRES by allocating room
for Krylov vectors etc

This function factors the low precision copy
and leaves the high precision matrix alone. ```mpglu!``` allocates
storage for ```basissize``` Kylov vectors and some other things
GMRES needs.
You get a factorization
object as output and can use ```\\``` to solve linear systems.
"""
function mpglu!(MPGA::MPGArray)
AL=MPGA.AL
AH=MPGA.AH
VStore=MPGA.VStore
KStore=MPGA.KStore
res=MPGA.residual
TL=eltype(AL)
(TL == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
MPF=MPGEFact(AH, AL, ALF, VStore, KStore, res, true)
return MPF
end

"""
mpglu(A::Array{TH,2}; TL=Float32, basissize=10) where TH <: Real

Combines the constructor of the multiprecision GMRES-ready array with the
factorization.

Step 1: build the MPGArray
Step 2: Call mpglu! to build the factorization object
"""
function mpglu(A::Array{TH,2}; TL=Float32, basissize=10) where TH <: Real
(TH==Float32) ? TL=Float16 : TL=TL
MPGA=MPGArray(A; basissize=basissize, TL=TL)
MPGF=mpglu!(MPGA)
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

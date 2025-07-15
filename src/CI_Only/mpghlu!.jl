# Factor a heavy MPArray and set it up for GMRES with \
# If you want to use it with IR (why?) then set gmresok=false
#
function mpglu!(MPH::MPHArray; gmresok = true, basissize = 10, residterm = residtermdefault)
    AH = MPH.AH
    TD = eltype(AH)
    res = MPH.residual
    sol = MPH.sol
    n = length(res)
    AStore = MPH.AStore
    AL = MPH.AL
    TF = eltype(AL)
    #
    # Factor in low precision
    #
    (TF == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
    #
    # Promote the low-precision lu
    #
    AStore .= TD.(AL)
    AF = LU(AStore, ALF.ipiv, ALF.info)
    anrm = TF(0.0)
    if gmresok
        VStore = zeros(TD, n, basissize)
        KStore = KIRstore(n, "gmres", TD)
        MPF = MPGHFact(AH, AL, AF, VStore, KStore, res, sol, true, residterm, anrm)
    else
        MPF = MPHFact(AH, AL, AF, res, sol, true, residterm, anrm)
    end
end

#
# Using a heavy MPArray with normal IR is insane. I use it only
# for CI.
#
function mphlu!(MPH::MPHArray)
    MPF = mpglu!(MPH; gmresok = false, residterm = residtermdefault)
end

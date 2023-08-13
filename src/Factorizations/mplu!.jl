"""
mplu!(MPA::MPArray)

Plain vanilla MPArray factorization.
Economy trianglar solve without interprecision transfers on the fly.
"""
function mplu!(MPA::MPArray)
    AH = MPA.AH
    AL = MPA.AL
    TL = eltype(AL)
    (TL == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
    MPF = MPLFact(AH, AL, AF)
    return MPF
end

function mplu!(MPA::MPEArray)
    AH = MPA.AH
    AL = MPA.AL
    TL = eltype(AL)
    (TL == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
    #AF=lu!(AL)
    MPF = MPLEFact(AH, AL, AF)
    return MPF
end

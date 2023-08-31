"""
mplu!(MPA::Union{MPArray,MPEArray})

Plain vanilla MPArray factorization.

The story on interprecision transfers is that 

- MPLFact downcasts the residual before the solve and avoids N^2 
  interprecision transfers. MPLFact factors MPArrays.

- MPLEFact factors MPEArrays and therefore does interprecision transfers 
  on the fly and incurs the N^2 interprecision transfer cost for that. 

  MPLEFact is what you must use if you plan to use the low precision 
  factorization as a preconditioner in IR-GMRES or you're working in 
  Float16 and the matrix is very ill-conditioned. MPLEFact factors 
  MPEArrays, which know to do interprecision transfers on-the-fly.

The 

Union{MPArray,MPEArray}

lets me use the ```on_the_fly``` trait to figure out what do to.

"""
function mplu!(MPA::Union{MPArray,MPEArray})
    AH = MPA.AH
    AL = MPA.AL
    TL = eltype(AL)
    residual=MPA.residual
    (TL == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
    # For the MPEArray
    if on_the_fly(MPA)
        MPF = MPLEFact(AH, AL, AF, residual)
    else
    # For the plain vanilla MPArray
        MPF = MPLFact(AH, AL, AF, residual)
    end
    return MPF
end

"""
mplu!(MPA::MPArray)

Plain vanilla MPArray factorization.

The story on interprecision transfers is that 

- MPLFact downcasts the residual before the solve and avoids N^2 
  interprecision transfers

- MPLEFact does interprecision transfers on the fly and incurs that N^2
  interprecision transfer cost. MPLEFact is what you must use if you plat
  to use the low precision factorization as a preconditioner in IR-GMRES.

MPLArray = Union{MPArray,MPEArray}

lets me use on_the_fly to figure out what do to.

"""
function mplu!(MPA::MPLArray)
    AH = MPA.AH
    AL = MPA.AL
    TL = eltype(AL)
    (TL == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
    if on_the_fly(MPA)
         MPF = MPLEFact(AH, AL, AF)
    else
         MPF = MPLFact(AH, AL, AF)
    end
    return MPF
end

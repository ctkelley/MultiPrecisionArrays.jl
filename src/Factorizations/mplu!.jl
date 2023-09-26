"""
mplu!(MPA::MPArray)

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
function mplu!(MPA::MPArray)
    AH = MPA.AH
    AL = MPA.AL
    TL = eltype(AL)
    residual=MPA.residual
    (TL == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
    # For the MPEArray
    on_the_fly=MPA.onthefly
    MPF = MPLFact(AH, AL, AF, residual, on_the_fly)
    return MPF
end


"""
mplu(A::Array{Float64,2}; TL=Float32, onthefly=false)

Combines the constructor of the multiprecision array with the
factorization.
"""
function mplu(A::Array{TH,2}; TL=Float32, onthefly=nothing) where TH <: Real
#
# If the high precision matrix is single, the low precision must be half.
#
(TH == Float32) && (TL = Float16)
#
# Unless you tell me otherwise, onthefly is true if low precision is half
# and false if low precision is single.
#
(onthefly == nothing ) && (onthefly = (TL==Float16))
MPA=MPArray(A; TL=TL, onthefly=onthefly)
MPF=mplu!(MPA)
return MPF
end

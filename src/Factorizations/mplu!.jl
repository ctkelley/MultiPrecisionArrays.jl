"""
mplu!(MPA::MPArray)

Plain vanilla MPArray factorization: Factor the low precision copy
and leave the high precision matrix alone. You get a factorization
object as output and can use ```\\``` to solve linear systems.

The story on interprecision transfers is that you can set the Boolean
```onthefly``` when you construct the MPArray. If you use ```mplu```
then you get the defaults

- If ```onthefly == false ``` then the solver downcasts the residual 
before the solve and avoids N^2 interprecision transfers.

- If ```onthefly == true``` then the solver does interprecision transfers 
  on the fly and incurs the N^2 interprecision transfer cost for that. 

  ```onthefly == true``` is what you must use if you plan to use 
  the low precision 
  factorization as a preconditioner in IR-GMRES or you're working in 
  Float16 and the matrix is very ill-conditioned. 

  ```onthefly == nothing``` means you take the defaults.

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
mplu!(MPF::MPLFact,A::AbstractArray{TH,2}) where TH <: Real

Overwrite a multiprecision factorization MPF to reuse the
storage to make a multiprecision of a new matrix A.

This will, of course, trash the original factorization.

To use this do
```
MPF=mplu!(MPF,A)
```
Simply using 
```
mplu!(MPF,A)
```
(ie without explicitly returning MPF)

may not do what you want because the multiprecision factorization
structure is immutable and MPF.AF.info cannot be changed.

Reassigning MPF works and resuses almost all of the storage in the 
original array.
"""
function mplu!(MPF::MPLFact,A::AbstractArray{TH,2}) where TH
TF=eltype(MPF.AH)
(TF == TH) || error("Precision error in mplu!")
AH=MPF.AH
AH = A
TL = eltype(MPF.AL)
AL=MPF.AL
AL .= TL.(A)
(TL == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
MPF.AF.ipiv .= AF.ipiv
MPF = MPLFact(A, AL, AF, MPF.residual, MPF.onthefly)
return MPF
end



"""
mplu(A::AbstractArray{TH,2}; TL=Float32, onthefly=nothing) where TH <: Real

Combines the constructor of the multiprecision array with the
factorization. 

Step 1: build the MPArray 

Step 2: factor the low precision copy and return the factorization object
"""
function mplu(A::AbstractArray{TH,2}; TL=Float32, onthefly=nothing) where TH <: Real
#
# If the high precision matrix is single, the low precision must be half.
#
(TH == Float32) && (TL = Float16)
#
# Unless you tell me otherwise, onthefly is true if low precision is half
# and false if low precision is single.
#
(onthefly == nothing ) && (onthefly = (TL==Float16))
#
# Build the multiprecision array MPA
#
MPA=MPArray(A; TL=TL, onthefly=onthefly)
#
# Factor the low precision copy to get the factorization object MPF
#
MPF=mplu!(MPA)
return MPF
end

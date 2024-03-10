"""
mpblu!(MPBA::MPBArray)
Factor a MPBArray and set it up for BiCGSTAB-IR

This function factors the low precision copy
and leaves the high precision matrix alone. The constructor for
MPBArray allocates
storage for the things BiCGSTAB needs.

You get a factorization
object as output and can use ```\\``` to solve linear systems.
"""
function mpblu!(MPBA::MPBArray)
AL=MPBA.AL
AH=MPBA.AH
VStore=MPBA.VStore
KStore=MPBA.KStore
res=MPBA.residual
TF=eltype(AL)
(TF == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
MPF=MPBFact(AH, AL, ALF, VStore, KStore, res, true)
return MPF
end

"""
mpblu!(MPG::MPBFact, A::AbstractArray{TW,2}) where TW <: Real
Overwrite a multiprecision factorization MPF to reuse the
storage to make a multiprecision factorization of a new matrix A.

This will, of course, trash the original factorization.

To use this do
```
MPG=mpblu!(MPF,A)
```
Simply using
```
mpblu!(MPF,A) # Don't do this.
```
(ie without explicitly returning MPG)

may not do what you want because the multiprecision factorization
structure is immutable and MPF.AF.info cannot be changed.

Reassigning MPG works and resuses almost all of the storage in the
original array.
"""
function mpblu!(MPG::MPBFact, A::AbstractArray{TW,2}) where TW <: Real
TF=eltype(MPG.AH)
(TF == TW) || error("Precision error in mplu!")
AH=MPG.AH
AH = A
TF = eltype(MPG.AL)
AL=MPG.AL
AL .= TF.(A)
(TF == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
MPG.AF.ipiv .= AF.ipiv
VStore=MPG.VStore 
KStore=MPG.KStore
MPG=MPBFact(AH, AL, AF, VStore, KStore, MPG.residual, true)
return MPG
end



"""
mpblu(A::AbstractArray{TW,2}; TF=Float32) where TW <: Real

Combines the constructor of the multiprecision BiCGSTAB-ready array with the
factorization.

Step 1: build the MPBArray

Step 2: Call mpblu! to build the factorization object
"""
function mpblu(A::AbstractArray{TW,2}; 
          TF=Float32) where TW <: Real
(TW==Float32) ? TF=Float16 : TF=TF
MPBA=MPBArray(A; TF=TF)
MPBF=mpblu!(MPBA)
return MPBF
end

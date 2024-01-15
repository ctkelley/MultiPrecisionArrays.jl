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
TL=eltype(AL)
(TL == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
MPF=MPBFact(AH, AL, ALF, VStore, KStore, res, true)
return MPF
end

"""
mpblu!(MPG::MPBFact, A::AbstractArray{TH,2}) where TH <: Real
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
function mpblu!(MPG::MPBFact, A::AbstractArray{TH,2}) where TH <: Real
TF=eltype(MPG.AH)
(TF == TH) || error("Precision error in mplu!")
AH=MPG.AH
AH = A
TL = eltype(MPG.AL)
AL=MPG.AL
AL .= TL.(A)
(TL == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
MPG.AF.ipiv .= AF.ipiv
VStore=MPG.VStore 
KStore=MPG.KStore
MPG=MPBFact(AH, AL, AF, VStore, KStore, MPG.residual, true)
return MPG
end



"""
mpblu(A::AbstractArray{TH,2}; TL=Float32) where TH <: Real

Combines the constructor of the multiprecision BiCGSTAB-ready array with the
factorization.

Step 1: build the MPBArray

Step 2: Call mpblu! to build the factorization object
"""
function mpblu(A::AbstractArray{TH,2}; 
          TL=Float32) where TH <: Real
(TH==Float32) ? TL=Float16 : TL=TL
MPBA=MPBArray(A; TL=TL)
MPBF=mpblu!(MPBA)
return MPBF
end

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
function mpblu!(MPBA::MPBArray; residterm=residtermdefault)
AL=MPBA.AL
AH=MPBA.AH
VStore=MPBA.VStore
KStore=MPBA.KStore
res=MPBA.residual
sol=MPBA.sol
TF=eltype(AL)
(TF == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
MPF=MPBFact(AH, AL, ALF, VStore, KStore, res, sol, true, residterm)
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
residterm=MPG.residterm 
MPG=MPBFact(AH, AL, AF, VStore, KStore, MPG.residual, MPG.sol, true, residterm)
return MPG
end



"""
mpblu(A::AbstractArray{TW,2}; TF=Float32, TR=nothing, i
          residterm=residtermdefault) where TW <: Real

Combines the constructor of the multiprecision BiCGSTAB-ready array with the
factorization.

Step 1: build the MPBArray

Step 2: Call mpblu! to build the factorization object

Step 3: Allocate the internal vectors that BiCGSTAB needs.

IR-Krylov methods were designed as saviors of half precision. So, as we
did with mpglu, we will run through some examples. The difference here
is that there is no Kylov basis to store. 

Here's the example from ```mpglu```

You may not get exactly the same results for this example on
different hardware, BLAS, versions of Julia or this package. 
I am still playing with the termination criteria and the iteration
count could grow or shrink as I do that.

#Example
```jldoctest
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> N=4096; alpha=799.0; AD = I - alpha*Gmat(N); A = Float32.(AD);

# Try vanilla IR with TF=Float16

julia> xe=ones(Float32,N); b=A*xe; AF=mplu(A);

julia> mout=\\(AF, b; reporting=true);

julia> mout.rhist
5-element Vector{Float64}:
 9.88750e+01
 3.92435e+00
 3.34373e-01
 2.02045e-01
 2.24720e-01

# That does not look like convergence. What about TR=Float64?

julia> BF=mplu(A; TR=Float64);

julia> mout2=\\(BF, b; reporting=true);

julia> mout2.rhist
5-element Vector{Float64}:
 9.88750e+01
 3.92614e+00
 3.34301e-01
 2.01975e-01
 2.24576e-01

# This is where we were in the docs for ```mpglu```. I'll try IR-BiCGSTAB

julia> GF=mpblu(A; TR=Float64);

julia> mout3=\\(GF, b; reporting=true);

julia> mout3.rhist
4-element Vector{Float64}:
 9.88750e+01
 2.16858e-11
 8.38440e-13
 8.81073e-13

# That is very different from IR-GMRES. The reason is that there is no
# limit on Krylov iterations because there is no Krylov subspace.
# So, the linear work produces a large reduction in the residual in
# the first iterate.

# Of course, we got it right and solved the promoted problem.

julia> xp=Float64.(A)\\b; norm(xp-mout3.sol,Inf)
1.63376e-11
```

"""
function mpblu(A::AbstractArray{TW,2}; residterm=residtermdefault,
          TF=Float32, TR=nothing) where TW <: Real
(TW==Float32) ? TF=Float16 : TF=TF
MPBA=MPBArray(A; TF=TF, TR=TR)
MPBF=mpblu!(MPBA; residterm=residterm)
return MPBF
end

"""
mpglu!(MPGA::MPGArray; residterm=residtermdefault)
Factor a MPGArray using the allocations from the structure.

This function factors the low precision copy
and leaves the high precision matrix alone. The constructor 
for MPGArray allocates
storage for ```basissize``` Kylov vectors and some other things
GMRES needs.
You get a factorization
object as output and can use ```\\``` to solve linear systems.

The kwarg ```residterm``` sets the termination criterion. 
```residterm == true``` (default) terminates the iteration on 
small residuals.  ```residterm == false``` terminates the iteration on
small normwise backward errors. Look at the docs for details.
"""
function mpglu!(MPGA::MPGArray; residterm = residtermdefault)
    AL = MPGA.AL
    TF = eltype(AL)
    AH = MPGA.AH
    (TF == Float16) ? AX=AH : AX=AL
    residterm ? anrm=TF.(0.0) : anrm=opnorm(AX,1)
    VStore = MPGA.VStore
    KStore = MPGA.KStore
    res = MPGA.residual
    sol = MPGA.sol
    (TF == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
    MPF = MPGEFact(AH, AL, ALF, VStore, KStore, res, sol, true, 
                 residterm, anrm)
    return MPF
end

"""
mpglu!(MPG::MPGEFact, A::AbstractArray{TW,2}) where TW <: Real
Overwrite a multiprecision factorization MPF to reuse the
storage to make a multiprecision array with a new matrix A.

This will, of course, trash the original factorization.

To use this do
```
MPG=mpglu!(MPG,A)
```
Simply using
```
mpglu!(MPG,A) # Don't do this.
```
(ie without explicitly returning MPG)

may not do what you want because the multiprecision factorization
structure is immutable and MPF.AF.info cannot be changed.

Reassigning MPG works and resuses almost all of the storage in the
original array
"""
function mpglu!(MPG::MPGEFact, A::AbstractArray{TW,2}) where {TW<:Real}
    TF = eltype(MPG.AH)
    (TF == TW) || error("Precision error in mplu!")
    AH = MPG.AH
    AH = A
    TF = eltype(MPG.AL)
    AL = MPG.AL
    AL .= TF.(A)
    (TF == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
    MPG.AF.ipiv .= AF.ipiv
    VStore = MPG.VStore
    KStore = MPG.KStore
    residterm = MPG.residterm
    anrm = MPG.anrm
    MPG = MPGEFact(AH, AL, AF, VStore, KStore, MPG.residual, MPG.sol, true, residterm, anrm)
    return MPG
end



"""
mpglu(A::AbstractArray{TW,2}; TF=Float32, TR=nothing, 
    residterm=residtermdefault1, basissize=10) where TW <: Real

Combines the constructor of the multiprecision GMRES-ready array with the
factorization.

Step 1: build the MPGArray

      (a) Store A and b in precision TW.

      (b) Store the factorization (copy of A) in precision TF

      (c) Preallocate storage for the residual, a local copy of the
      solution, and the Krylov vectors in precision TR


Step 2: Call mpglu! to build the factorization object

The ```TR``` kwarg is the residual precision. Leave this alone unless you know
what you are doing. The default is ```nothing``` which tells the solver to
set ```TR = TW```. If you set ```TR``` it must be
a higher precision than TW and
you are essentially solving ```TR.(A) x = TR.(b)```
with IR with the factorization in ```TF```.

The case of interest here is ```TW = Float32; TF = Float16; TR = Float64```.
IR-GMRES with those precisions is sometimes the savior of half precision.

The default basis size is 10. You can (and should) play with that if you
are not happy with the results. If you are having trouble storing enough
vectors, IR-BiCGSTAB ```mpblu``` is worth a shot.

Take this example, please.

You may not get exactly the same results for this example on
different hardware, BLAS, versions of Julia or this package. 
I am still playing with the termination criteria and the iteration
count could grow or shrink as I do that.

## Example
```jldoctest
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> N=4096; alpha=799.0; AD = I - alpha*Gmat(N); A = Float32.(AD);

# Try vanilla IR with TF=Float16

julia> xe=ones(Float32,N); b=A*xe; AF=mplu(A); 

julia> mout=\\(AF, b; reporting=true);

julia> mout.rhist
4-element Vector{Float64}:
 9.88750e+01
 3.92435e+00
 3.34373e-01
 2.02045e-01

# That does not look like convergence. What about TR=Float64?

julia> BF=mplu(A; TR=Float64);

julia> mout2=\\(BF, b; reporting=true);

julia> mout2.rhist
5-element Vector{Float64}:
4-element Vector{Float64}:
 9.88750e+01
 3.92451e+00
 3.34292e-01
 2.02204e-01

# Can Float16 be saved? How 'bout IR-GMRES with a generous basissize.

julia> GF=mpglu(A; TR=Float64, basissize=20);

julia> mout3=\\(GF, b; reporting=true)

julia> mout3.rhist
4-element Vector{Float64}:
 9.88750e+01
 7.59075e-03
 1.48842e-05
 2.17281e-07

# Shazam! Did we get the solution of the promoted problem?

julia> xp=Float64.(A)\\b; norm(xp-mout3.sol,Inf)
9.02494e-06

# Maybe. 

```
"""
function mpglu(
    A::AbstractArray{TW,2};
    TF = Float32,
    TR = nothing,
    residterm = residtermdefault,
    basissize = 10,
) where {TW<:Real}
    (TW == Float32) ? TF = Float16 : TF = TF
    MPGA = MPGArray(A; basissize = basissize, TF = TF, TR = TR)
    MPGF = mpglu!(MPGA; residterm = residterm)
    return MPGF
end



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
    residterm ? anrm=TF.(0.0) : anrm=norm(AL,1)
    VStore = MPGA.VStore
    KStore = MPGA.KStore
    res = MPGA.residual
    sol = MPGA.sol
    (TF == Float16) ? ALF = hlu!(AL) : ALF = lu!(AL)
#    anrm=TF.(0.0)
    MPF = MPGEFact(AH, AL, ALF, VStore, KStore, res, sol, true, 
                 residterm, anrm)
    return MPF
end

"""
mpglu!(MPG::MPGEFact, A::AbstractArray{TW,2}) where TW <: Real
Overwrite a multiprecision factorization MPF to reuse the
storage to make a multiprecision of a new matrix A.

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

# Can Float16 be saved? How 'bout IR-GMRES?

julia> GF=mpglu(A; TR=Float64);

julia> mout3=\\(GF, b; reporting=true);

julia> mout3.rhist
6-element Vector{Float64}:
 9.88750e+01
 2.23211e-04
 9.61252e-09
 1.26477e-12
 8.10019e-13
 8.66862e-13

# Shazam! Did we get the solution of the promoted problem?

julia> xp=Float64.(A)\\b; norm(xp-mout3.sol,Inf)
9.43157e-12

# Yup.

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


#
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
        KStore = kstore(n, "gmres")
        MPF = MPGHFact(AH, AL, AF, VStore, KStore, res, sol, true, 
               residterm,anrm)
    else
        MPF = MPHFact(AH, AL, AF, res, sol, true, residterm,anrm)
    end
end

#
# Using a heavy MPArray with normal IR is insane. I use it only
# for CI.
#
function mphlu!(MPH::MPHArray)
    MPF = mpglu!(MPH; gmresok = false, residterm = residtermdefault)
end

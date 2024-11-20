"""
mplu!(MPA::MPArray; residterm=residtermdefault)

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

The kwarg ```residterm``` sets the termination criterion. 
```residterm == true``` (default) terminates the iteration on 
small residuals.  ```residterm == false``` terminates the iteration on
small normwise backward errors. Look at the docs for details.

If you want to use static arrays with this stuff, use the 
mutable @MArray constructor

"""
function mplu!(MPA::MPArray; residterm = residtermdefault)
    AH = MPA.AH
    AL = MPA.AL
    TF = eltype(AL)
    residual = MPA.residual
    sol = MPA.sol
    (TF == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
    # For the MPEArray
    on_the_fly = MPA.onthefly
    anrm=TF.(0.0)
    MPF = MPLFact(AH, AL, AF, residual, sol, on_the_fly, residterm, anrm)
    return MPF
end

"""
mplu!(MPF::MPLFact,A::AbstractArray{TW,2}) where TW <: Real

Overwrite a multiprecision factorization MPF to reuse the
storage to make a multiprecision of a new matrix A.

This will, of course, trash the original factorization.

To use this do
```
MPF=mplu!(MPF,A)
```
Simply using 
```
mplu!(MPF,A) # Don't do this!
```
(ie without explicitly returning MPF)

may not do what you want because the multiprecision factorization
structure is immutable and MPF.AF.info cannot be changed.

Reassigning MPF works and resuses almost all of the storage in the 
original array.

If you want to use static arrays with this stuff, use the 
mutable @MArray constructor
"""
function mplu!(MPF::MPLFact, A::AbstractArray{TW,2}) where {TW}
    TH = eltype(MPF.AH)
    (TH == TW) || error("Precision error in mplu!")
    AH = MPF.AH
    AH = A
    TF = eltype(MPF.AL)
    AL = MPF.AL
    AL .= TF.(A)
    (TF == Float16) ? AF = hlu!(AL) : AF = lu!(AL)
    MPF.AF.ipiv .= AF.ipiv
    residterm = MPF.residterm
    anrm=TF.(0.0)
    MPF = MPLFact(A, AL, AF, MPF.residual, MPF.sol, MPF.onthefly, 
          residterm, anrm)
    return MPF
end



"""
mplu(A::AbstractArray{TW,2}; TF=nothing, TR=nothing, residterm=residtermdefault,
                    onthefly=nothing) where TW <: Real

Combines the constructor of the multiprecision array with the
factorization. 

Step 1: build the MPArray. 

      (a) Store A and b in precision TW. 

      (b) Store the factorization (copy of A) in precision TF

      (c) Preallocate storage for the residual and a local copy of the
      solution in precision TR

Step 2: factor the low precision copy and return the factorization object


The ```TR``` kwarg is the residual precision. Leave this alone unless you know
what you are doing. The default is ```nothing``` which tells the solver to
set ```TR = TW```. If you set ```TR``` it must be
a higher precision than TW and
you are essentially solving ```TR.(A) x = TR.(b)```
with IR with the factorization in ```TF```. The MPArray structure stores
the solution and the residual in precision ```TR``` and so
the residual computation is done via ```TR.(A) x```. The
interprecision transfers are on the fly. So, the storage cost is the matrix,
and the copy in the factorization precision.

 The classic case is ```TW = TF = Float32``` and ```TR = Float64```. The nasty
part of this is that you must store TWO copies of the matrix. One for
the residual computation and the other to overwrite with the factors.
I do not think this is a good deal unless A is seriously ill-conditioned.
My support for this is through ```mplu```. To do this you must put the
```TR``` kwarg explicitly in your call to ```mplu```.   

The kwarg ```residterm``` sets the termination criterion.
```residterm == true``` (default) terminates the iteration on
small residuals.  ```residterm == false``` terminates the iteration on
small normwise backward errors. Look at the docs for details.

You may not get exactly the same results for this example on
different hardware, BLAS, versions of Julia or this package. 
I am still playing with the termination criteria and the iteration
count could grow or shrink as I do that.

## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> n=31; alpha=Float32(1.0);

julia> G=Gmat(n, Float32);

julia> A = I + alpha*G;

julia> b = A*ones(Float32,n);

# use mpa with TF=TW=Float32 and TR=Float64

julia> AF = mplu(A; TF=Float32, TR=Float64, onthefly=true);

# Solve and save the iteration history

julia> mout = \\(AF, b; reporting=true);

julia> mout.rhist
5-element Vector{Float64}:
 1.12500e+00
 1.74765e-07
 3.10862e-14
 6.66134e-16
 2.22045e-16

# What does this mean. I'll solve the promoted problem. TR.(A) x = b

julia> AD=Float64.(A);

julia> xd = AD\\b;

julia> norm(xd - mout.sol,Inf)
1.11022e-15

# So IR with TR > TW solves a promoted problem.
```


"""
function mplu(
    A::AbstractArray{TW,2};
    TF = nothing,
    TR = nothing,
    residterm = residtermdefault,
    onthefly = nothing,
) where {TW<:Real}
    #
    # If the high precision matrix is single, the low precision must be half
    # unless you're planning on using a high-precision residual where TR > TW
    # and also factoring in the working precision, so TW == TF.
    #
    #
    (TR == nothing) && (TR = TW)
    TFdef = Float32
    (TW == Float32) && (TFdef = Float16)
    (TF == nothing) && (TF = TFdef)
    #
    # Unless you tell me otherwise, onthefly is true if low precision is half
    # and false if low precision is single.
    #
    (onthefly == nothing) && (onthefly = (TF == Float16))
    #
    # IF TF = TW then something funny is happening with the residual precision.
    #
    (TF == TW) && (onthefly = true)
    #
    # Build the multiprecision array MPA
    #
    MPA = MPArray(A; TF = TF, TR = TR, onthefly = onthefly)
    #
    # Factor the low precision copy to get the factorization object MPF
    #
    MPF = mplu!(MPA; residterm = residterm)
    return MPF
end

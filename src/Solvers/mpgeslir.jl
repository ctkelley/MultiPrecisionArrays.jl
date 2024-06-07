"""
mpgeslir(MPA::MPArray, b; TR=Float16, reporting = false, verbose = true)

I do not export this function. The idea is that you use ```mpglu```
and do not touch either the constructor or the solver directly.

Use a multi-precision factorization to solve a linear system with
plain vanilla iterative refinement.

This version is analogous to ```A\\b``` and combines the factorization
and the solve. You start with MPA=MPArray(A) and then pass MPA
to mpgeslir and combine the factorization and the solve. 

You can also get the multiprecision factorization directly with
```
MPF=mplu(A)
```
and then pass ```MPF``` to mpgeslir.

I use this to get some timing results and it's also convenient
if you want to do factor and solve in one statement. 

You can also get this with ```x = MPA\\b```.

If you set the kwarg ```reporting``` to true you can get the IR
residual history. The output of 
```
x = MPA\\b
```
or
```
x=MPF\\b
```
is the solition. The output of 
```
mout = \\(MPA,b; reporting=true)
```
or
```
mout = \\(MPF,b; reporting=true)
```
is a structure. ```mpout.sol``` is the solution. ```mpout.rhist```
is the residual history. mpout also contains the datatypes TW for
high precision and TF for low precision.

## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> N=4096; A = I - 800.0 * Gmat(N); b=ones(N);

julia> MPF=mplu(A);

julia> mout=\\(MPF, b; reporting=true);

julia> mout.rhist
6-element Vector{Float64}:
 1.00000e+00
 5.36483e-02
 1.57977e-05
 5.10232e-09
 7.76756e-12
 9.90008e-12

# Stagnation after four IR iterations

julia> [mout.TW mout.TF]
1×2 Matrix{DataType}:
 Float64  Float32

```

The ```TR``` kwarg is the residual precision. Leave this alone unless you know
what you are doing. The default is ```Float16``` which tells the solver to
set ```TR = TW```. If you use this option, then
```TR``` is a higher precision than TW and when you set TR
you are essentially solving ```TR.(A) x = TR.(b)``` 
with IR with the factorization
in ```TF``` and the residual computation done via ```A TR.(x)``` and 
interprecision transfers on the fly. So, the storage cost is the matrix,
and the copy in the factorization precision.

 The classic case is ```TW = TF = Float32``` and ```TR = Float64```. The nasty
part of this is that you must store TWO copies of the matrix. One for
the residual computation and the other to overwrite with the factors.
I do not think this is a good deal unless A is seriously ill-conditioned.
My support for this through ```mplu```. To do this you must put the 
```TR``` kwarg explicitly in your call to the solver.

## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> n=31; alpha=Float32(1.0);

julia> G=Gmat(n, Float32);

julia> A = I + alpha*G;

julia> b = A*ones(Float32,n);

# use mpa with TF=TW=Float32

AF = mplu(A; TF=Float32, onthefly=true);

# now set TR=Float64 with the kwarg and solve

mout = \\(AF, b; TR=Float64, reporting=true);

# The solution and the residuals are in double. The iteration drives
# the residual (evaluated in double) to close to double precision roundoff

julia> mout.rhist
4-element Vector{Float64}:
 1.12500e+00
 2.65153e-07
 6.90559e-14
 6.66134e-16

# What does this mean. I'll solve the promoted problem. TR.(A) x = b

julia> AD=Float64.(A);

julia> xd = AD\\b;

julia> norm(xd - mout.sol,Inf)
1.11022e-15

# Is that better?
```

"""
function mpgeslir(MPA::MPArray, b; TR=Float16, reporting = false, verbose = true)
# Factor MPA and return Factorization object
MPF=mplu!(MPA)
# Call mpgeslir for the solve
xi=\(MPF, b; TR=TR, reporting=reporting, verbose=verbose)
return xi
end

"""
mpgeslir(AF::MPFact, b; TR=Float16, reporting=false, verbose=true)

I do not export this function. The idea is that you use ```mpglu```
and do not touch either the constructor or the solver directly.

Use a multi-precision factorization to solve a linear system with
plain vanilla iterative refinement.

This version is analogous to ```A\\b``` and combines the factorization
and the solve. You start with MPA=MPArray(A) and then pass MPA
to mpgeslir and combine the factorization and the solve. 

You can also get the multiprecision factorization directly with
```
MPF=mplu(A)
```
and then pass ```MPF``` to mpgeslir.

I use this to get some timing results and it's also convenient
if you want to do factor and solve in one statement. 

You can also get this with ```x = MPA\\b```.

If you set the kwarg ```reporting``` to true you can get the IR
residual history. The output of 
```
x = MPA\\b
```

Use a multi-precision factorization to solve a linear system with
plain vanilla iterative refinement.

MPFact is a union of all the MultiPrecision factorizations in the package. 
The triangular solver will dispatch on the various types depending on
how the interprecision transfers get done.
"""
function mpgeslir(AF::MPFact, b; TR=Float16, reporting = false, verbose = true)
    #
    # What kind of problem are we dealing with?
    #
    mpdebug = false
    normtype = Inf
    TB = eltype(b)
    MPStats = getStats(AF)
    TF = MPStats.TF
    TW = MPStats.TW
    r = AF.residual
    onthefly=AF.onthefly
    (TR == Float16) && (TR=TW)
    # If I'm computing a high precision residual, TS=TR
    # and I must do interprecision transfers on the fly.
    HiRes = (eps(TR) < eps(TW))
    HiRes && (onthefly=true)
    #
    # TFact is the precision of the factors; should be TF
    # unless we're dealing with a heavy MPArray for CI
    #
    TFact = MPStats.TFact
    TFok = ((TF == TFact) || (typeof(AF) == MPHFact))
    TFok  || error("TF is supposed to be TFact")
    #
    # Are the precisions consistent? If not, I have a bug somewhere.
    # Otherwise, set the tolerance on the iteration to 10*eps.
    # If the iteration can't meet the tolerance, terminate when
    # the residual norms stagnate (res_old > .9 res_new)
    #
    (TW == TB) || error("inconsistent precisions; A and b must have same type")
    tolf = eps(TR)*TR.(10.0)
    #
    # Keep the records and accumulate the statistics. 
    #
    Meth = MPStats.Meth
    verbose && println(
        Meth,
        ": High precision = $TW, Low precision = $TF, Factorization storage precision = $TFact, Residual precision = $TR"
    )
    #
    # Showtime!
    #
    AD = AF.AH
    bnrm = norm(b, normtype) 
    anrm = opnorm(AD, normtype)
    bsc = b
    AFS = AF.AF
    bS = TFact.(bsc)
    #
    # Initialize the iteration. I initialize to zero. That makes the
    # iteration count the same as the high precision matvec and the 
    # triangular sovles
    #
    x = zeros(TR, size(b))
    xnrm = norm(x, normtype)
    #
    # Initial residual
    #
    r .= b 
    tol = tolf * bnrm
    rs = bS
#
#   Keep the books. Test for excessive residual precision.
#  
    ERes = (eps(TR) < eps(Float64))
    HiRes ? (THist = TR) : (THist=Float64)
    rhist = Vector{THist}()
    rnrm = TR(norm(r, normtype))
    rnrmx = rnrm * TR(2.0)
    oneb = TR(1.0)
    itc = 0
    #
    # Put initial residual norm into the history and iterate.
    #
    push!(rhist, rnrm)
    # Store r and x in the residual precision if TR is not TW
    HiRes ? rloop=TR.(r) : rloop=r
    HiRes ? xloop=TR.(x) : xloop=x
    # Solve loop
    # while (rnrm > tol) && (rnrm <= .9*rnrmx)
    while (rnrm > (anrm * xnrm + bnrm) *tolf) && (rnrm <= .9*rnrmx)
        #
        # Scale the residual
        #
        rloop ./= rnrm
        #
        # Use the low-precision factorization
        #
        rloop .= IRTriangle!(AF, rloop, rs, verbose)
        #
        # Undo the scaling
        #
        rloop .*= rnrm
        #
        # Update the solution and residual
        #
        xloop .+= rloop
        mul!(rloop, AD, xloop)
        rloop .*= -oneb
        axpy!(oneb, bsc, rloop)
        rnrmx = rnrm
        rnrm = norm(rloop, normtype)
        itc += 1
        push!(rhist, rnrm)
        mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
        #
        # If the residual norm increased, complain.
        #
#        complain_resid = (rnrm >= rnrmx) && (rnrm > 1.e3 * tol)
#        complain_resid && println("IR Norm increased: $rnrm, $rnrmx, $tol")
         xnrm=norm(xloop,normtype)
    end
    x = xloop
    verbose && println("Residual history = $rhist")
    if reporting
        return (rhist = rhist, sol = x, TW = TW, TF = TF, TFact = TFact)
    else
        return x
    end
end

function getTF(AF::MPFact)
    TF = eltype(AF.AL)
    if is_heavy(AF)
    TFact = eltype(AF.AH)
    else
    TFact = eltype(AF.AL)
    end
    return (TF, TFact)
end

function getStats(AF)
    TW = eltype(AF.AH)
    (TF, TFact) = getTF(AF)
    MPStats = MPIRStats(TW, TF, TFact)
    return MPStats
end

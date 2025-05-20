"""
mpgeslir(MPA::MPArray, b; reporting = false, verbose = true)

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

You may not get exactly the same results for this example on
different hardware, BLAS, versions of Julia or this package. 
I am still playing with the termination criteria and the iteration
count could grow or shrink as I do that.

## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> N=4096; A = I - 800.0 * Gmat(N); b=ones(N);

julia> MPF=mplu(A);

julia> mout=\\(MPF, b; reporting=true);

julia> mout.rhist
6-element Vector{Float64}:
 9.90000e+01
 3.65823e-03
 6.17917e-07
 8.74678e-11
 2.04636e-12
 2.03215e-12

# Stagnation after four IR iterations

julia> [mout.TW mout.TF]
1×2 Matrix{DataType}:
 Float64  Float32

```

"""
function mpgeslir(MPA::MPArray, b; reporting = false, verbose = true)
    # Factor MPA and return Factorization object
    MPF = mplu!(MPA)
    # Call mpgeslir for the solve
    xi = \(MPF, b; reporting = reporting, verbose = verbose)
    return xi
end

"""
mpgeslir(AF::MPFact, b; reporting=false, verbose=true)

I do not export this function. The idea is that you use ```mplu```
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
function mpgeslir(AF::MPFact, b; reporting = false, verbose = true)
    #
    # What kind of problem are we dealing with?
    #
    mpdebug = false
    N = length(b)
    normtype = Inf
    #    normtype = 1
    TB = eltype(b)
    MPStats = getStats(AF)
    TF = MPStats.TF
    TW = MPStats.TW
    r = AF.residual
    TR = eltype(r)
    x = AF.sol
    x .*= TW(0.0)
    #    xa = AF.sol .* TRA(0.0)
    #    xr = AF.sol .* TR(0.0)
    #    x=xr
    onthefly = AF.onthefly
    # If I'm computing a high precision residual, TS=TR
    # and I must do interprecision transfers on the fly.
    HiRes = (eps(TR) < eps(TW))
    HiRes && (onthefly = true)
    #
    # TFact is the precision of the factors; should be TF
    # unless we're dealing with a heavy MPArray for CI
    #
    TFact = MPStats.TFact
    TFok = ((TF == TFact) || (typeof(AF) == MPHFact))
    TFok || error("TF is supposed to be TFact")
    #
    # Are the precisions consistent? If not, I have a bug somewhere.
    # Otherwise, set the tolerance on the iteration to 10*eps.
    # If the iteration can't meet the tolerance, terminate when
    # the residual norms stagnate (res_old > rrf * res_new)
    #
    (TW == TB) || error("inconsistent precisions; A and b must have same type")
    residterm = AF.residterm
    term_data = termination_settings(TW, residterm)
    tolf = term_data.tolf
    rrf = term_data.redmax
    litmax = term_data.litmax
    AD = AF.AH
    #
    # I compute the norm of AF if needed in single
    # Half is still to slow.
    #
    anrm = AF.anrm
    #
    # Keep the records and accumulate the statistics. 
    #
    Meth = MPStats.Meth
    verbose && println(
        Meth,
        ": High precision = $TW, Low precision = $TF, Factorization storage precision = $TFact, Residual precision = $TR",
    )
    #
    # Showtime!
    #
    AD = AF.AH
    bnrm = norm(b, normtype)
    bsc = b
    AFS = AF.AF
    #
    # Initialize the iteration. I initialize to zero. That makes the
    # iteration count the same as the high precision matvec and the 
    # triangular sovles
    #
    xnrm = norm(x, normtype)
    #
    # Initial residual
    #
    r .= b
    tol = tolf
    onthefly ? (rs = ones(TF, 1)) : (rs = zeros(TF, size(b)))
#
# If TR > TW then do the solve in TW after computing r in TR
#
    (TR == TW) || (rs = zeros(TW,size(b)))
    #
    #   Keep the books. Test for excessive residual precision.
    #  
    ERes = (eps(TR) < eps(Float64))
    HiRes ? (THist = TR) : (THist = Float64)
    rhist = Vector{THist}()
    rnrm = TR(norm(r, normtype))
    rnrmx = rnrm * 1.e6
    oneb = TR(1.0)
    itc = 0
    #
    # Put initial residual norm into the history and iterate.
    #
    push!(rhist, rnrm)
    # Store r and x in the residual precision if TR is not TW
    HiRes ? rloop = TR.(r) : rloop = r
    HiRes ? xloop = TR.(x) : xloop = x
    #    rrf = 0.5
    #    rrf = term_data.redmax
    # Solve loop
    tol=(anrm * xnrm + bnrm) * tolf
    while (rnrm > tol) && (rnrm <= rrf*rnrmx) && (itc < litmax)
        #
        # Scale the residual
        #
        rloop ./= rnrm
        #
        # Use the low-precision factorization
        # The residual is overwritten with the correction here.
        #
        rloop .= IRTriangle!(AF, rloop, rs, verbose)
        #
        # Undo the scaling
        #
        rloop .*= rnrm
        #
        # Update the solution and residual
        #
        x .+= rloop
        xloop .= TR.(x)
#        xloop .+= rloop
        mul!(rloop, AD, xloop)
        #
        # After mul! the residual is overwritten with Ax
        #
        rloop .*= -oneb
        axpy!(oneb, bsc, rloop)
        #
        # and now the residual is b-Ax like it needs to be
        #
        rnrmx = rnrm
        rnrm = norm(rloop, normtype)
        itc += 1
        push!(rhist, rnrm)
        xnrm = norm(xloop, normtype)
        tol = tolf * (anrm * xnrm + bnrm)
        mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
        #
        # If the residual norm increased, complain.
        #
        complain_resid = mpdebug && (rnrm >= rnrmx) && (rnrm > 1.e3 * tol)
        complain_resid && println("IR Norm increased: $rnrm, $rnrmx, $tol")
    end
#    x = xloop
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

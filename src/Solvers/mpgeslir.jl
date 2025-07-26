"""
mpgeslir(MPA::MPArray, b; reporting = false, verbose = false,
         term_parms=term_parms_default)

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
is the residual history. ```mpout.dhist``` is the history of the
norms of the corrections.  mpout also contains the datatypes TW for
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
7-element Vector{Float64}:
 1.00000e+00
 4.24311e-02
 9.03500e-05
 9.59277e-08
 6.12346e-11
 9.75675e-12
 8.49343e-12

# Stagnation after six IR iterations

julia> mout.dhist
6-element Vector{Float64}:
 2.00420e+02
 3.15198e-01
 4.48500e-04
 4.57912e-07
 5.88713e-10
 3.62408e-11

# No correction for the final itertation.

julia> [mout.TW mout.TF]
1×2 Matrix{DataType}:
 Float64  Float32

```

"""
function mpgeslir(
    MPA::MPArray,
    b;
    reporting = false,
    verbose = false,
    term_parms = term_parms_default,
)
    # Factor MPA and return Factorization object
    MPF = mplu!(MPA)
    # Call mpgeslir for the solve
    xi = \(MPF, b; reporting = reporting, verbose = verbose)
    return xi
end

"""
mpgeslir(AF::MPFact, b; reporting=false, verbose=false, 
         term_parms=term_parms_default)

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
function mpgeslir(
    AF::MPFact,
    b;
    reporting = false,
    verbose = false,
    term_parms = term_parms_default,
)
    #
    # What kind of problem are we dealing with?
    #
    mpdebug = false
    N = length(b)
    normtype = Inf
    (TW, TF, TR, TFact) = Types_IR_Init(AF, b, normtype)
    #
    # Initialize the iteration. I initialize to zero. That makes the
    # iteration count the same as the high precision matvec and the 
    # triangular solves
    #
    (x, r, rs, xnrm, bnrm, anrm) = Solver_IR_Init(AF, b, normtype)
    rnrm = bnrm
    #
    #  get the termination data 
    #
    residterm = AF.residterm
    tolf = termination_settings(TW, term_parms, residterm)
    Rmax = term_parms.Rmax
    litmax = term_parms.litmax
    AD = AF.AH
    #
    # Keep the records and accumulate the statistics. 
    #
    verbose && println(
        "High precision = $TW, Low precision = $TF, Factorization storage precision = $TFact, Residual precision = $TR",
    )
    #
    # Showtime!
    #
    AD = AF.AH
    bsc = b
    tol = tolf
    #
    #   Keep the books. 
    #  
    rhist = Vector{TR}()
    dhist = Vector{TW}()
    rnrmx = rnrm * 1.e6
    oneb = TR(1.0)
    itc = 0
    #
    # Put initial residual norm into the history and iterate.
    #
    push!(rhist, rnrm)
    # Store r and x in the residual precision if TR is not TW
    (TR == TW) ? rloop = r : rloop = TR.(r)
    (TR == TW) ? xloop = x : xloop = TR.(x)
    # Solve loop
    tol=(anrm * xnrm + bnrm) * tolf
    dnormold=1.0
    etest=true
    while (rnrm > tol) && (rnrm <= Rmax*rnrmx) && (itc < litmax) || etest
        #
        # Scale the residual
        #
        rloop ./= rnrm
        #
        # Use the low-precision factorization
        # The residual is overwritten with the correction here.
        # If TR > TW, then rs = TW.(r) and I use that as the rhs
        # for the working precision solve.
        #
        rloop .= IRTriangle!(AF, rloop, rs, verbose)
        #
        # Undo the scaling
        #
        rloop .*= rnrm
        #
        # Update the solution and residual
        #
        dnorm = norm(rloop, normtype)
        push!(dhist, dnorm)
        drat=dnorm/dnormold
        dnormold=dnorm
        # High precision residual? Use ||d|| in termination.
        etest = (eps(TR) < eps(TW)) && (drat < Rmax) || (itc==0)
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
        return (rhist = rhist, dhist = dhist, sol = x, TW = TW, TF = TF, TFact = TFact)
    else
        return x
    end
end

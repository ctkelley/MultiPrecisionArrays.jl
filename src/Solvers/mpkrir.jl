"""
mpkrir(AF::MPKFact, b; reporting=false, verbose=false, 
           mpdebug=false, term_parms=term_parms_default)

Krylov-IR solver 

I do not export this function. The idea is that you use ```mpglu```
and do not touch either the constructor or the solver directly.

Use a multi-precision factorization to solve a linear system with
IR-Krylov metods

This is the generic solver used by GMRES-IR and BiCGSTAB-IR. You use
the correct MPKFact = Union{MPGFact,MPBFact} structure and mpkrir
will do the right thing. 

You should not be calling this directly. Use ```\\``` to solve
linear systems with the multiprecision factorization and to 
use the optional kwargs.

We overload the backslash operator to call mpkrir for a multiprecision
MPGFact factorization. So if ```MPA``` is an MPGArray and  
```
AF = mpglu!(MPA)
```
Then ```AF\\b``` maps to
```
mpkrir(AF, b)
```
which does the GMRES-IR solve. You can also get the multiprecision
factoriztion directly with
```
AF = mpglu(A)
```
which builds the mutliprcision MPGArray and then factors the low
preicsion copy.

Similarly if  ```MPA``` is an MPBArray. Then
```
AF = mpblu!(MPA)
```
Then ```AF\\b``` maps to
```
mpkrir(AF, b)
```
which does the BiCGSTAB-IR solve.

You can also use the ```\\``` operator to harvest iteration statistics.

You may not get exactly the same results for this example on
different hardware, BLAS, versions of Julia or this package. 
I am still playing with the termination criteria and the iteration
count could grow or shrink as I do that.

## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> N=4096; A = I - 800.0 * Gmat(N); b=ones(N);

julia> AF=mpglu(A);

julia> solout=\\(AF, b; reporting=true);

# Correct result?

julia> x=solout.sol; norm(b-A*x,Inf)
8.20877e-12

# Look at the residual history

julia> solout.rhist

4-element Vector{Float64}:
 1.00000e+00
 1.25881e-10
 8.58902e-12
 8.20877e-12

# and the correction norm history. No correction for the final iteration.

julia> solout.dhist
3-element Vector{Float64}:
 2.00129e+02
 1.05241e-10
 2.28629e-11


# Stagnation after the 2nd iteration. Now the Krylovs/iteration

julia> solout.khist
3-element Vector{Int64}:
 4
 5
 5
# 4-5 Krylovs per iteration.

BiCGSTAB works the same way.


"""
function mpkrir(
    AF::MPKFact,
    b;
    reporting = false,
    verbose = false,
    mpdebug = false,
    term_parms = term_parms_default,
)
    #
    # Which Krylov method are we talking about?
    ktype = kmeth(AF)
    krylov_ok = (ktype == "GMRES") || (ktype == "BiCGSTAB")
    krylov_ok || error("$ktype is not supported")
    #
    normtype = Inf
    # Harvest the types of everyhting
    (TW, TF, TR, TFact) = Types_IR_Init(AF, b)
    # 
    # Initialize the iteration. I initialize to zero. That makes the
    # iteration count the same as the high precision matvec and the
    # solves for the defect
    #
    (x, xres, r, rs, bnrm, rhist, dhist) = Solver_IR_Init(AF, b, normtype)
    rnrm = bnrm
    rnrmx = rnrm * 1.e6
    #
    #  get the termination data
    #
    (tolf, Rmax, litmax, tol) = Term_Init(AF, term_parms, bnrm)
    #   
    # Tell 'em more than they need to know. 
    #  
    ir_vmsg(TW, TF, TFact, TR, verbose)
    # Copy x to the residual precision if TR is not TW
    #    (TR == TW) ? xres=x : xres = copy(r)
    #
    # Initialize Krylov-IR
    #
    khist = Vector{Int64}()
    atvd = copy(x)
    x0 = zeros(TW, size(b))
    MP_Data = (MPF = AF, atv = atvd, x0 = x0, ktype = ktype, eta = tolf)
    # 
    dnormold=1.0
    etest=true
    itc = 0
    #
    # Krylov-IR loop
    #
    while (rnrm > tol) && (rnrm <= Rmax * rnrmx) && (itc < litmax) || etest
        #
        # If TR > TW, then rs = TW.(r) and I use that as the rhs
        # for the working precision solve.
        #
        kout = IRKsolve!(AF, r, rs, rnrm; itc = itc, verbose = verbose, MP_Data = MP_Data)
        r .= TR.(kout.sol)*rnrm
        #
        # Manage the results and keep the books
        push!(khist, length(kout.reshist))
        dnorm=norm(r, normtype)
        drat=dnorm/dnormold
        dnormold=dnorm
        push!(dhist, dnorm)
        #
        # High precision residual? Use ||d|| in termination.
        #
        etest = (eps(TR) < eps(TW)) && (drat < Rmax) || (itc==0)
        #
        # Update the solution and residual
        # xres is in the residual precision to get a
        # high precision residual.
        #
        d=kout.sol
        d .*= rnrm
        x .+= d
        xnrm = norm(x, normtype)
        r = Resid_IR(r, x, xres, b, AF)
        rnrmx = rnrm
        rnrm = norm(r, normtype)
        itc += 1
        push!(rhist, rnrm)
        tol = tolf * (bnrm + AF.anrm * xnrm)
        #
        # Debugging? Report iteration data
        #
        ir_debug_msg(mpdebug, itc, tol, rnrm, rnrmx, dnorm, drat)
        #
    end
    outdat=closeout(AF, rhist, dhist, khist, x, TW, TF, verbose)
    if reporting
        return outdat
    else
        return x
    end
end

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
    (TW, TF, TR, TFact) = Types_IR_Init(AF, b, normtype)
    # 
    # Initialize the iteration. I initialize to zero. That makes the
    # iteration count the same as the high precision matvec and the
    # solves for the defect
    #
    (x, r, rs, xnrm, bnrm, anrm) = Solver_IR_Init(AF,b,normtype)
    rnrm = bnrm
    rnrmx = rnrm * 1.e6
    itc = 0
    #
    #  get the termination data
    #
    tolf = termination_settings(AF, term_parms)
    Rmax = term_parms.Rmax
    litmax = term_parms.litmax
    AD = AF.AH
    bsc = TR.(b)
    #
    # Initialize Krylov-IR
    #
    rhist = Vector{TR}()
    dhist = Vector{TW}()
    onetb = TR(1.0)
    khist = Vector{Int64}()
    push!(rhist, rnrm)
    #
    # Krylov-IR loop
    #
    itc = 0
    atvd = copy(r)
    xloop=copy(r)
    x0 = zeros(TW, size(b))
    eta = tolf
    MP_Data = (MPF = AF, atv = atvd, x0 = x0,
                ktype = ktype, eta = eta)
    tol = tolf * (bnrm + anrm * xnrm)
    dnormold=1.0
    etest=true
    while (rnrm > tol) && (rnrm <= Rmax * rnrmx) && (itc < litmax) || etest
        #
        # Scale the residual 
        #
        r ./= rnrm
        # Solve the correction equation
        # If TR > TW, then rs = TW.(r) and I use that as the rhs
        # for the working precision solve.
        #
        kout = IRKsolve!(AF, r, rs; MP_Data = MP_Data)
        # Manage the results and keep the books
        push!(khist, length(kout.reshist))
        itcp1 = itc + 1
#        winner = kout.idid ? " $ktype converged" : " $ktype failed"
        #
        # Make some noise
        #
        irk_msg(itcp1, kout, ktype, verbose)
#        verbose && (println(
#            "Krylov stats: Iteration $itcp1 :",
#            length(kout.reshist),
#            " iterations",
#            "  ",
#            winner,
#        ))
        #
        # Overwrite the residual with the correction
        #
        r .= TR.(kout.sol)
        #
        # Undo the scaling
        #
        r .*= rnrm
        #
        # Update the solution and residual
        # xloop is in the residual precision to get a
        # high precision residual.
        #
        x .+= r
        dnorm=norm(r, normtype)
        drat=dnorm/dnormold
        dnormold=dnorm
        push!(dhist, dnorm)
        # High precision residual? Use ||d|| in termination.
        etest = (eps(TR) < eps(TW)) && (drat < Rmax) || (itc==0)
        xloop .= TR.(x)
        mul!(r, AD, xloop)
        r .*= -onetb
        axpy!(1.0, bsc, r)
        rnrmx = rnrm
        rnrm = norm(r, normtype)
        itc += 1
        push!(rhist, rnrm)
        #        tol = tolf * bnrm
        mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
        #
        # If the residual norm increased, complain.
        #
        complain_resid = mpdebug && (rnrm >= rnrmx) && (rnrm > 1.e3 * tol)
        complain_resid && println("IR Norm increased: $rnrm, $rnrmx, $tol")
        xnrm = norm(x, normtype)
        tol = tolf * (bnrm + anrm * xnrm)
    end
    verbose && println("Residual history = $rhist")
    if reporting
        #        TF = eltype(AF.AL); TW = eltype(AF.AH);
        return (rhist = rhist, dhist = dhist, khist = khist, sol = x, TW = TW, TF = TF)
    else
        return x
    end
end

function irk_msg(itcp1, kout, ktype, verbose)
winner = kout.idid ? " $ktype converged" : " $ktype failed"
verbose && (println(
            "Krylov stats: Iteration $itcp1 :",
            length(kout.reshist),
            " iterations",
            "  ",
            winner,
        ))
end



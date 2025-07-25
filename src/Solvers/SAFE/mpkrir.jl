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
    (x, r, rs, anrm, onthefly, HiRes) = Solver_IR_Init(AF, b)
    #    x = AF.sol
    AD = AF.AH
    # remember that eps(TR) = 2 * unit roundoff
    residterm = AF.residterm
    tolf = termination_settings(TW, term_parms, residterm)
    Rmax = term_parms.Rmax
    litmax = term_parms.litmax
    #
    # I compute the norm of AF if needed in single
    # Half is still too slow.
    #
    anrm = AF.anrm
    #
    n = length(b)
    onetb = TR(1.0)
    bsc = TR.(b)
    bnorm = norm(b, normtype)
    xnorm = TR(0.0)
    AD = AF.AH
    #
    # Initialize Krylov-IR
    #
    #    r = AF.residual
    r .= TR.(b)
    rnrm = norm(r, normtype)
    bnrm = norm(b, normtype)
    rnrmx = rnrm * TR(2.0)
    rhist = Vector{TR}()
    dhist = Vector{TW}()
    khist = Vector{Int64}()
    push!(rhist, rnrm)
    eta = tolf
    #
    # Krylov-IR loop
    #
    itc = 0
    VF = AF.VStore
    kl_store = AF.KStore
    atvd = copy(r)
    xloop=copy(r)
    MP_Data = (MPF = AF, atv = atvd, TF = TF, TW = TW)
    #    rrf = 0.5
    #    rrf = term_data.Rmax
    tol = tolf * (bnorm + anrm * xnorm)
    dnormold=1.0
    etest=true
    while (rnrm > tol) && (rnrm <= Rmax * rnrmx) && (itc < litmax) || etest
        x0 = zeros(TR, n)
        #
        # Scale the residual 
        #
        r ./= rnrm
        #
        # Solve the correction equation with a Krylov method
        #
        if ktype == "GMRES"
            kout = kl_gmres(
                x0,
                r,
                MPhatv,
                VF,
                eta,
                MPhptv;
                pdata = MP_Data,
                side = "left",
                kl_store = kl_store,
            )
        elseif ktype == "BiCGSTAB"
            kout = kl_bicgstab(
                x0,
                r,
                MPhatv,
                VF,
                eta,
                MPhptv;
                pdata = MP_Data,
                side = "left",
                kl_store = kl_store,
            )
        end
        push!(khist, length(kout.reshist))
        itcp1 = itc + 1
        winner = kout.idid ? " $ktype converged" : " $ktype failed"
        #
        # Make some noise
        #
        verbose && (println(
            "Krylov stats: Iteration $itcp1 :",
            length(kout.reshist),
            " iterations",
            "  ",
            winner,
        ))
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
        xloop=TR.(x)
        mul!(r, AD, xloop)
        r .*= -onetb
        axpy!(1.0, bsc, r)
        rnrmx = rnrm
        rnrm = norm(r, normtype)
        itc += 1
        push!(rhist, rnrm)
        #        tol = tolf * bnorm
        mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
        #
        # If the residual norm increased, complain.
        #
        complain_resid = mpdebug && (rnrm >= rnrmx) && (rnrm > 1.e3 * tol)
        complain_resid && println("IR Norm increased: $rnrm, $rnrmx, $tol")
        xnorm = norm(x, normtype)
        tol = tolf * (bnorm + anrm * xnorm)
    end
    verbose && println("Residual history = $rhist")
    if reporting
        #        TF = eltype(AF.AL); TW = eltype(AF.AH);
        return (rhist = rhist, dhist = dhist, khist = khist, sol = x, TW = TW, TF = TF)
    else
        return x
    end
end

# Matrix-vector product for Krylov-IR
function MPhatv(x, pdata)
    atv = pdata.atv
    mul!(atv, pdata.MPF.AH, x)
    return atv
end

# Preconditioner-vector product for Krylov-IR
function MPhptv(x, pdata)
    TF=pdata.TF
    TW=pdata.TW
    #
    # ldiv! is not doing well for one case, so I hide from it.
    #
    if (TF == Float16) && (TW == Float32)
        x .= pdata.MPF.AF\x;
    else
        ldiv!(pdata.MPF.AF, x)
    end
    return x
end

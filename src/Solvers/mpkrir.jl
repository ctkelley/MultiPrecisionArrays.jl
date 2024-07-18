"""
mpkrir(AF::MPKFact, b; reporting=false, 
                               verbose=false, mpdebug=false)

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

Krylov-IR solver 

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

## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> N=4096; A = I - 800.0 * Gmat(N); b=ones(N);

julia> AF=mpglu(A);

julia> solout=\\(AF, b; reporting=true);

# Correct result?

julia> x=solout.sol; norm(b-A*x,Inf)
9.12193e-12

# Look at the residual history

julia> solout.rhist
5-element Vector{Float64}:
 1.00000e+00
 1.16784e-10
 8.47566e-12
 7.95053e-12
 9.12193e-12
# Stagnation after the 3rd iteration. Now the Krylovs/iteration

julia> solout.khist
4-element Vector{Int64}:
 4
 5
 4
 5
# 4-5 Krylovs per iteration.

BiCGSTAB works the same way.


"""
function mpkrir(AF::MPKFact, b; reporting = false, 
                verbose = false, mpdebug = false)
    #
    # Which Krylov method are we talking about?
    ktype=kmeth(AF)
    krylov_ok = (ktype=="GMRES") || (ktype=="BiCGSTAB")
    krylov_ok || error("$ktype is not supported")
    #
    normtype = Inf
    x = AF.sol
    TR = eltype(x)
    x .*= TR(0.0)
    # remember that eps(TR) = 2 * unit roundoff
    residterm=AF.residterm
    residterm ? tf=1.0 : tf=.9
    tolf = tf*eps(TR)
    n = length(b)
    onetb = TR(1.0)
    bsc = copy(b)
#    x = zeros(TR, size(b))
    bnorm = norm(b, normtype)
    xnorm = TR(0.0)
    #
    AFS = AF.AF
    AD = AF.AH
    residterm ?  anrm = 0.0 : anrm = opnorm(AD, 1)
#    anorm = opnorm(AD,normtype)
    #
    # Initialize Krylov-IR
    #
    r = AF.residual
    r .= b
    rnrm = norm(r, normtype)
    rnrmx = rnrm * TR(2.0)
    rhist = Vector{TR}()
    khist = Vector{Int64}()
    push!(rhist, rnrm)
    eta = tolf
    #
    # Krylov-IR loop
    #
    itc = 0
    VF=AF.VStore
    normdec = true
    kl_store = AF.KStore
    atvd=copy(r)
    MP_Data = (MPF = AF, atv = atvd)
    tol = tolf *(bnorm + anrm *xnorm)
    rrf = .9
    while (rnrm > tol) && ( rnrm <= rrf * rnrmx )
        x0 = zeros(TR, n)
        #
        # Scale the residual 
        #
        r ./= rnrm
        #
        # Solve the correction equation with a Krylov method
        #
        if ktype=="GMRES"
        kout = kl_gmres(x0, r, MPhatv, VF, eta, MPhptv; 
               pdata = MP_Data, side = "left", kl_store=kl_store)
        elseif ktype=="BiCGSTAB"
        kout = kl_bicgstab(x0, r, MPhatv, VF, eta, MPhptv;
               pdata = MP_Data, side = "left", kl_store=kl_store)
        end
        #
        # Make some noise
        #
        push!(khist,length(kout.reshist))
        itcp1 = itc + 1
        winner = kout.idid ? " $ktype converged" : " $ktype failed"
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
        #
        x .+= r
        mul!(r, AD, x)
        r .*= -onetb
        axpy!(1.0, bsc, r)
        rnrmx = rnrm
        rnrm = norm(r, normtype)
        itc += 1
        push!(rhist, rnrm)
        tol = tolf * bnorm
        mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
        #
        # If the residual norm increased, complain.
        #
        (rnrm >= rnrmx) && (normdec = false)
        ~normdec && mpdebug && (rnrm >= rnrmx) && println("Residual norm increased")
    xnorm=norm(x,normtype)
    tol = tolf *(bnorm + anrm *xnorm)
    end
    verbose && println("Residual history = $rhist")
    if reporting
        TF = eltype(AF.AL)
        TW = eltype(AF.AH)
        return (rhist = rhist, khist = khist,
               sol = x, TW = TW, TF = TF)
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
    ldiv!(pdata.MPF.AF, x)
    return x
end


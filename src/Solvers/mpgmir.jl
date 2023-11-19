"""
mpgmir(AF::MPGFact, b; reporting=false, verbose=false, mpdebug=false)

GMRES-IR solver 

We overload the backslash operator to call mprmir for a multiprecision
MPGFact factorization. So if ```MPA``` is an MPArray and
```
AF = mpglu!(MPA)
```
Then ```AF\\b``` maps to
```
mpgmir(AF, b)
```
You should call ```mpgir``` explicitly if you want the iteration
statistics.

The three kwargs are ```reporting```, ```verbose```, and ```mpdebug```.
Of these ```reporting``` is most useful. The other two are for my
debugging pleasure.

When you do
```
mpout = mpgmir(NAF, b; reporting=true)
```
You get a structure where

```mpout.sol``` is the solution

```mpout.rhist``` is the residual history for IR

```mpout.khist``` is the krylovs/IR iteration 

Other parts of ```mpout``` are the high and low precisions 
```TH``` and  ```TL```.


## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> N=4096; A = I - 800.0 * Gmat(N); b=ones(N);

julia> MPA=MPArray(A); AF=mpglu!(MPA);

julia> mpout=mpgmir(AF, b; reporting=true);

julia> x=mpout.sol; norm(b-A*x,Inf)
8.92664e-12

julia> mpout.rhist
4-element Vector{Float64}:
 6.40000e+01
 3.32046e-09
 1.17624e-10
 1.17893e-10
# Stagnation after the second iteration

julia> mpout.khist
3-element Vector{Int64}:
 4
 5
 4
# 4-5 Krylovs per iteration.

julia> mpout.TH
Float64

julia> mpout.TL
Float32

# Repeat the experiment with low precision TL=Float16 (half)

julia> MPA=MPArray(A; TL=Float16); AF=mpglu!(MPA); 

# You can use backslash too

julia> mpout=\\(AF, b; reporting=true);

julia> x=mpout.sol; norm(b-A*x,Inf)
8.65075e-12
# Residual is as good as the TL=Float32 case.

julia> mpout.rhist
5-element Vector{Float64}:
 6.40000e+01
 2.00140e-03
 2.05307e-07
 1.16612e-10
 1.18166e-10

# Stagnaton after 3 iterations

julia> mpout.khist
5-element Vector{Int64}:
 10
 10
 10
 10
# The default basissize=10 so we are taking all the GMRES iterations we
# can at each iteration.

julia> 

julia> mpout.TH
Float64

julia> mpout.TL
Float16


```

"""
function mpgmir(AF::MPGFact, b; reporting = false, 
                verbose = false, mpdebug = false)
    #
    normtype = 2
    TB = eltype(b)
    tolf = TB(10.0)*eps(TB)
    n = length(b)
    onetb = TB(1.0)
    bsc = copy(b)
    x = zeros(TB, size(b))
    bnorm = norm(b, normtype)
    #
    AFS = AF.AF
    AD = AF.AH
    #
    # Initialize GMRES-IR
    #
    r = AF.residual
    r .= b
#    r = copy(b)
    rnrm = norm(r, normtype)
    rnrmx = rnrm * TB(2.0)
    rhist = Vector{TB}()
    khist = Vector{Int64}()
    push!(rhist, rnrm)
#    eta = TB(1.e-8)
    eta = tolf
    #
    # GMRES-IR loop
    #
    itc = 0
    VF=AF.VStore
    normdec = true
#    kl_store=kstore(n,"gmres")
    kl_store = AF.KStore
    atvd=copy(r)
    MP_Data = (MPF = AF, atv = atvd)
    while (rnrm > tolf * bnorm) && ( rnrm <= .99 * rnrmx )
        x0 = zeros(TB, n)
        #
        # Scale the residual 
        #
        r ./= rnrm
        #
        # Solve the correction equation with GMRES
        #
        kout = kl_gmres(x0, r, MPhatv, VF, eta, MPhptv; 
               pdata = MP_Data, side = "left", kl_store=kl_store)
        #
        # Make some noise
        #
        push!(khist,length(kout.reshist))
        itcp1 = itc + 1
        winner = kout.idid ? " GMRES converged" : " GMRES failed"
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
        r .= TB.(kout.sol)
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
    end
    verbose && println("Residual history = $rhist")
    if reporting
        TL = eltype(AF.AL)
        return (rhist = rhist, khist = khist,
               sol = x, TH = TB, TL = TL)
    else
        return x
    end
end

#function MPhatv(x, MPF::MPGFact)
function MPhatv(x, pdata)
#    atv = MPF.AH * x
    atv = pdata.atv
    mul!(atv, pdata.MPF.AH, x)
#    atv = pdata.MPF.AH * x
    return atv
end

function MPhptv(x, pdata)
#function MPhptv(x, MPF::MPGFact)
#    ptv = MPF.AF \ x
#    ptv = pdata.ptv
#    ptv .= x
     ldiv!(pdata.MPF.AF, x)
#     ldiv!(pdata.MPF.AF, ptv)
#    ptv .= pdata.MPF.AF \ x
#    return ptv
    return x
end

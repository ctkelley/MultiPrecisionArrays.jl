"""
mpkrir(AF::MPKFact, b; reporting=false, verbose=false, mpdebug=false)

Krylov-IR solver 

This is the generic solver used by GMRES-IR and BiCGSTAB-IR. You use
the correct MPKFact = Union{MPGFact,MPBFact} structure and mpkrir
will do the right thing. 

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
    # Initialize Krylov-IR
    #
    r = AF.residual
    r .= b
    rnrm = norm(r, normtype)
    rnrmx = rnrm * TB(2.0)
    rhist = Vector{TB}()
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
    while (rnrm > tolf * bnorm) && ( rnrm <= .99 * rnrmx )
        x0 = zeros(TB, n)
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


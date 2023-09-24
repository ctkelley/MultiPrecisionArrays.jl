"""
mpgmir(AF::MPGFact, b; reporting=false, verbose=false, mpdebug=false)

GMRES-IR solver 
"""
function mpgmir(AF::MPGFact, b; reporting = false, 
                verbose = false, mpdebug = false)
    #
    normtype = Inf
    TB = eltype(b)
    irtol = (TB == Float64) ? 1.e-14 : 1.e-7
    tolf = TB.(10.0)*eps(TB)
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
    while (rnrm > tolf * bnorm) && ( rnrm <= .9 * rnrmx )
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
        tol = irtol * bnorm
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
        TFact = eltype(AF.AL)
        return (rhist = rhist, sol = x, TH = TB, TL = TL, TFact = TFact)
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

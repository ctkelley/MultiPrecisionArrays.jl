"""
mpgmir(AF::MPGHFact, b; reporting=false, verbose=true, mpdebug=false)

Prototype GMRES-IR solver
"""
function mpgmir(AF::MPGHFact, b; reporting = false, 
                verbose = false, mpdebug = false)
    #
    normtype = Inf
    TB = eltype(b)
    irtol = (TB == Float64) ? 1.e-14 : 1.e-7;
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
    r = copy(x)
    mul!(r, AD, x)
    r .*= -onetb
    axpy!(onetb, bsc, r)
    rnrm = norm(r, normtype)
    rnrmx = rnrm * TB(1.1)
    rhist = Vector{TB}()
    push!(rhist, rnrm)
    eta = TB(1.e-6)
    #
    # GMRES-IR loop
    #
    itc = 0
    VF = zeros(TB, n, 20)
    normdec=true
    while (rnrm > irtol*bnorm) && (itc < 10) && normdec
        x0 = zeros(TB, n)
        #
        # Scale the residual 
        #
        r ./= rnrm
        #
        # Solve the correction equation with GMRES
        #
        kout = kl_gmres(x0, r, MPhatv, VF, eta, MPhptv; 
               pdata = AF, side = "left")
        #
        # Make some noise
        #
        verbose && (itc +=1; println("Krylov stats: Iteration $itcc ", 
                     kout.reshist, "  ", kout.idid))
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
        mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
        #
        # If the residual norm increased, complain.
        #
        (rnrm >= rnrmx) && (println("Residual norm increased"); normdec=false)
    end
    verbose && println("Residual history = $rhist")
    if reporting
        TL = eltype(AF.ALow)
        TFact = eltype(AF.ALow)
        return (rhist = rhist, sol = x, TH = TB, TL = TL, TFact = TFact)
    else
        return x
    end
end

function MPhatv(x, MPHF::MPGHFact)
    atv = MPHF.AH * x
    return atv
end

function MPhptv(x, MPHF::MPGHFact)
    ptv = MPHF.AF \ x
    return ptv
end

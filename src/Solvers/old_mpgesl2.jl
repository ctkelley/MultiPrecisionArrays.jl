"""
mpgesl2(AF::MPFact, b; reporting=false, verbose=true, expensive=false)

Use a multi-precision factorization to solve a linear system.
"""
function mpgesl2(AF::MPFact, b; reporting = false, verbose = true)
    # # What kind of problem are we dealing with?
    #
    mpdebug = false
    normtype = Inf
    TB = eltype(b)
    MPStats = getStats(AF)
    TL = MPStats.TL
    TH = MPStats.TH
    #
    # TFact is the precision of the factors
    #
    TFact = MPStats.TFact
    #
    # Are the precisions consistent? If not, I have a bug somewhere.
    # If so, set the completely arbitrary tolerances
    #
    (TH == TB) || error("inconsistent precisions")
    (TH == Float64) ? tolf = 1.e-13 : tolf = 1.e-6
    #
    # Keep the records and accumulate the statistics. 
    #
    Meth = MPStats.Meth
    if verbose
        println( Meth,": High precision = $TH, Low precision = $TL,
             Factorization storage precision = $TFact",
        )
    end
    #
    # Showtime!
    #
    AD = AF.AH
    bnrm = norm(b, normtype)
    bsc = b
    AFS = AF.AF
    bS = TFact.(bsc)
    #
    # Initialize the iteration. I am still thinking about how I want
    # to do this. For now I initialize to zero.
    #
    x = zeros(size(b))
#    if (typeof(AF) == MPGFact)
#        x = zeros(size(b))
#    else
#        x = zeros(size(b))
#    end
    #
    # Initial residual
    #
    r = copy(x)
    mul!(r, AD, x)
    r .*= -1.0
    axpy!(1.0, bsc, r)
    tol = tolf * bnrm
    rs = bS
    rhist = Vector{Float64}()
    rnrm = norm(r, normtype)
    rnrmx = rnrm * 1.1
    itc = 0
    #
    # Put initial residual norm into the history and iterate.
    #
    push!(rhist, rnrm)
    while (rnrm > tol) && (rnrm < rnrmx)
        r ./= rnrm
        #
        # Hungry IR?
        #
        if (typeof(AF) <: MPHFact)
#            ldiv!(AFS, r)
            r .= IRTriangle!(AF, r, rs)
        elseif (typeof(AF) <: MPLFact) || (typeof(AF) <: MPLEFact)
            #
            # Solver for the correction; AFS in factorization precision
            #
            r .= IRTriangle!(AF, r, rs)
#            rs .= TFact.(r)
#            ldiv!(AFS, rs)
#            r .= TH.(rs)
            #
            # IR-GMRES? call my kl_gmres code. I need to make a special-purpose
            # code for this.
            #
        elseif (typeof(AF) <: MPGFact)
            x0 = zeros(size(r))
            eta = 1.e-4
            V = AF.VH
            kout = kl_gmres(x0, r, MPhatv, V, eta, MPhptv; pdata = AF, side = "left")
            #
            # Make noise?
            #
            if verbose
                itcc = itc + 1
                println("Krylov stats: Iteration $itcc ", kout.reshist, "  ", kout.idid)
            end
            #
            # Overwite r with d
            #
            r .= kout.sol
        else
            error("missing MP Fact type")
        end
        #
        # Undo the scaling
        #
        r .*= rnrm
        #
        # Update the solution and residual
        #
        x .+= r
        mul!(r, AD, x)
        r .*= -1.0
        axpy!(1.0, bsc, r)
        rnrmx = rnrm
        rnrm = norm(r, normtype)
        itc += 1
        push!(rhist, rnrm)
        mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
        #
        # If the residual norm increased, complain.
        #
        (rnrm >= rnrmx) && println("Norm increased")
    end
    if verbose
        println("Residual history = $rhist")
    end
    if reporting
        return (rhist = rhist, sol = x, TH = TH, TL = TL, TFact = TFact)
    else
        return x
    end
end

function getTL(AF)
    if (typeof(AF) <: MPLFact)
        TL = eltype(AF.AL)
    elseif (typeof(AF) <: MPLEFact)
        TL = eltype(AF.AL)
    elseif (typeof(AF) <: MPHFact)
        TL = eltype(AF.AS)
    elseif (typeof(AF) <: MPGFact)
        TL = eltype(AF.AS)
    else
        TX = typeof(AF)
        error("illegal MPFact type $TX")
    end
    return TL
end

function getMeth(AF)
    if (typeof(AF) <: MPLFact) || (typeof(AF) <: MPLEFact)
        Meth = "IR"
    elseif (typeof(AF) <: MPHFact)
        Meth = "Hungry IR"
    elseif (typeof(AF) <: MPGFact)
        Meth = "IRGM"
    else
        TX = typeof(AF)
        error("illegal MPFact type $TX")
    end
    return Meth
end

function getStats(AF)
    TH = eltype(AF.AH)
    TL = getTL(AF)
    TFact = eltype(AF.AL)
    if (typeof(AF) == MPGFact)
        MPStats = MPGStats()
    else
        MPStats = MPIRStats(TH, TL, TFact)
    end
    return MPStats
end

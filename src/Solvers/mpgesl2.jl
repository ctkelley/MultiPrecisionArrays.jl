"""
mpgesl2(AF::MPFact, b; reporting=false, verbose=true)

Use a multi-precision factorization to solve a linear system with
plain vanilla iterative refinement.
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
    verbose &&
        println( Meth,": High precision = $TH, Low precision = $TL, Factorization storage precision = $TFact")
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
    x = zeros(TB,size(b))
#    if (typeof(AF) == MPGFact)
#        x = zeros(size(b))
#    else
#        x = zeros(size(b))
#    end
    #
    # Initial residual
    #
    oneb=TB(1.0)
    r = copy(x)
    mul!(r, AD, x)
    r .*= -oneb
    axpy!(oneb, bsc, r)
    tol = tolf * bnrm
    rs = bS
    rhist = Vector{Float64}()
    rnrm = norm(r, normtype)
    rnrmx = rnrm * TB(1.1)
    itc = 0
    #
    # Put initial residual norm into the history and iterate.
    #
    push!(rhist, rnrm)
    while (rnrm > tol) && (rnrm < rnrmx)
        #
        # Scale the residual
        #
        r ./= rnrm
        #
        # Use the low-precision factorization
        #
        r .= IRTriangle!(AF, r, rs, verbose)
        #
        # Undo the scaling
        #
        r .*= rnrm
        #
        # Update the solution and residual
        #
        x .+= r
        mul!(r, AD, x)
        r .*= -oneb
        axpy!(oneb, bsc, r)
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
    verbose && println("Residual history = $rhist")
    if reporting
        return (rhist = rhist, sol = x, TH = TH, TL = TL, TFact = TFact)
    else
        return x
    end
end

function getTL(AF::MPLFacts)
TL = eltype(AF.AL)
return TL 
end

function getTL(AF::MPHFact)
TL = eltype(AF.AS)
return TL 
end


#function getTL(AF)
#    if (typeof(AF) <: MPLFact)
#        TL = eltype(AF.AL)
#    elseif (typeof(AF) <: MPLEFact)
#        TL = eltype(AF.AL)
#    elseif (typeof(AF) <: MPHFact)
#        TL = eltype(AF.AS)
#    elseif (typeof(AF) <: MPGFact)
#        TL = eltype(AF.AS)
#    else
#        TX = typeof(AF)
#        error("illegal MPFact type $TX")
#    end
#    return TL
#end

function getStats(AF)
    TH = eltype(AF.AH)
    TL = getTL(AF)
    TFact = eltype(AF.AL)
#    if (typeof(AF) == MPGFact)
#        MPStats = MPGStats()
#    else
        MPStats = MPIRStats(TH, TL, TFact)
#    end
    return MPStats
end

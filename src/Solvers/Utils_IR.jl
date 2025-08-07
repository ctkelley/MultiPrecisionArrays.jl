#
# closeout gathers the data for output if reporting=true
#
function closeout(AF::MPFact, rhist, dhist, x, TW, TF, verbose)
    verbose && println("Residual history = $rhist")
    return (rhist = rhist, dhist = dhist, sol = x, TW = TW, TF = TF)
end
#
function closeout(AF::MPKFact, rhist, dhist, khist, x, TW, TF, verbose)
    verbose && println("Residual history = $rhist")
    return (rhist = rhist, dhist = dhist, khist = khist, sol = x, TW = TW, TF = TF)
end
#
# These functions print various messages if
# verbose or mpdebug are set to true.
#
function ir_debug_msg(mpdebug, itc, tol, rnrm, rnrmx, dnorm, drat)
    mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
    complain_resid = mpdebug && (rnrm >= rnrmx) && (rnrm > 1.e3 * tol)
    complain_resid && println("IR resid Norm increased: $rnrm, $rnrmx, $tol")
    complain_resid && println("Correction norm $dnorm. Ratio $drat.")
end

function ir_vmsg(TW, TF, TFact, TR, verbose)
    verbose && println(
        "High precision = $TW, Low precision = $TF, Factorization storage precision = $TFact, Residual precision = $TR",
    )
end

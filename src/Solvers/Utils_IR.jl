function ir_debug_msg(mpdebug, itc, tol, rnrm, rnrmx)
    mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
    complain_resid = mpdebug && (rnrm >= rnrmx) && (rnrm > 1.e3 * tol)
    complain_resid && println("IR Norm increased: $rnrm, $rnrmx, $tol")
end

function ir_vmsg(TW, TF, TFact, TR, verbose)
    verbose && println(
        "High precision = $TW, Low precision = $TF, Factorization storage precision = $TFact, Residual precision = $TR",
    )
end

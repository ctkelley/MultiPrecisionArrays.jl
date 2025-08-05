function ir_debug_msg(mpdebug, itc, tol, rnrm, rnrmx)
    mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
    complain_resid = mpdebug && (rnrm >= rnrmx) && (rnrm > 1.e3 * tol)
    complain_resid && println("IR Norm increased: $rnrm, $rnrmx, $tol")
end

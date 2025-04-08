function termination_settings(TR, residterm)
    #
    # I am continually tweaking this stuff
    #
    # I will terminate IR using either small relative residual
    # (residterm = true, the default)
    # || r || < Cr eps(TR) || b || 
    #
    # or small normwise backward error
    # (residterm = false)
    # || r || < Ce eps(TR) (|| A || || x || + || b ||) 
    #
    # I use the L1 norm for A because that is less expensive that 
    # Linfty.  || A || is as expensive as a few IR iterations, 
    # so it's not the default.
    #
    # I also stop when I see stagnation, ie
    # || r_new || > redmax || r_old ||
    #
    (Cr, Ce, redmax) = termdata()
#    (Cr, Ce, Rr, Re) = termdata()
#    residterm ? redmax = Rr : redmax = Re
    residterm ? tf = Cr : tf = Ce
    tolf = eps(TR) * tf
    #
#    anrm = 0.0
#    term_out = (tolf = tolf, anrm = anrm, redmax = redmax)
    term_out = (tolf = tolf, redmax = redmax)
    return term_out
end

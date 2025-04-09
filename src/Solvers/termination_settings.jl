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
    # || r_new || > Rmax || r_old ||
    #
    # Cr, Ce, and Rmax are set in the main module MultiPrecisionArrays.jl
    # as fields in a mutable struct. You can change them, but if you do
    # that (not recommended) understand that the parameters are global
    # within the module. Take care with changing the parameters while
    # using the solvers in a multi-threaded code. 
    #
    (Cr, Ce, Rmax) = termdata()
    residterm ? tf = Cr : tf = Ce
    tolf = eps(TR) * tf
    term_out = (tolf = tolf, redmax = Rmax)
    return term_out
end

function restore_default_parms(t::TERM = term_parms)
    #
    # The default parameters are stored as a constant mutable struct
    # term_parms_default in MultiPrecisionArrays.jl
    #
    t.Cr = term_parms_default.Cr
    t.Ce = term_parms_default.Ce
    t.Rmax = term_parms_default.Rmax
    return t
end

function update_parms(t::TERM = term_parms; Cr = 20.0, Ce = 1.0, Rmax = 0.5)
    #
    # The parameters live in a mutable struct term_parms in 
    # MultiPrecisionArrays.jl and are initialized to the default values. 
    #
    t.Cr = Cr
    t.Ce = Ce
    t.Rmax = Rmax
    return t
end

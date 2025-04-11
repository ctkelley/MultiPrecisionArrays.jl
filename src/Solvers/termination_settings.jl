function termination_settings(TR, residterm)
    #
    # I am continually tweaking this stuff. 
    # I do not export this and do not encourage you to play with it.
    # If you want to mess with it anyway, the functions in this file
    # and the ***Termination criteria defaults*** in MultiPrecisionArrays.jl
    # will let you get started.
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
    (Cr, Ce, Rmax, litmax) = termdata()
    residterm ? tf = Cr : tf = Ce
    tolf = eps(TR) * tf
    term_out = (tolf = tolf, redmax = Rmax, litmax=litmax)
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
    t.litmax = term_parms_default.litmax
    return t
end

"""
    update_parms(t::TERM = term_parms; Cr = 20.0, Ce = 1.0,
           Rmax = 0.5, litmax=1000
)

C. T. Kelley, 2025

Update the termination parameters in MultiPrecisionArrays.

This changes the values of the parameters. I do not recommend that
you mess with this unless you have a good reason. One such reason
might be that the default limit on the number of iterations (1000,
essentially infinite) is giving you many iterations that make very
little progress. Changing Rmax is probably a better way to control this.

Note that the default values for the kwargs are the defaults, so
unless you specify otherwise, any parameter will take its default value.

We store the parameters in a mutable structure TERM
```
mutable struct TERM
       Cr::Real
       Ce::Real
       Rmax::Real
       litmax::Int
end
```
where the fields are the parameters. The parameters we use in 
the solvers are in a global TERM structure defined in the main module.
```
term_parms=TERM(20.0, 1.0, .5, 1000)
```
As you can see, I start with the defaults. When you make changes, you
write into this structure and it is your job to keep track of what you
did.

To query the parameters type the name of the structure
```
julia> term_parms
MultiPrecisionArrays.TERM(2.00000e+01, 1.00000e+00, 5.00000e-01, 1000)
```
To change Cr from 20 to 40 type
```
julia> update_parms(; Cr=40.0)
MultiPrecisionArrays.TERM(4.00000e+01, 1.00000e+00, 5.00000e-01, 1000)
```
To change Ce to 5 and revert Cr back to the default
```
julia> update_parms(; Ce=5.0)
MultiPrecisionArrays.TERM(2.00000e+01, 5.00000e+00, 5.00000e-01, 1000)
```
Finally, to get the defaults back, enter no changes.
```
julia> update_parms(;)
MultiPrecisionArrays.TERM(2.00000e+01, 1.00000e+00, 5.00000e-01, 1000)
```


"""
function update_parms(t::TERM = term_parms; Cr = 20.0, Ce = 1.0,
         Rmax = 0.5, litmax=1000)
    #
    # The parameters live in a mutable struct term_parms in
    # MultiPrecisionArrays.jl and are initialized to the default values.
    #
    t.Cr = Cr
    t.Ce = Ce
    t.Rmax = Rmax
    t.litmax = litmax
    return t
end


"""
    update_parms(t::TERM = term_parms; Cr = 1.0, Ce = 1.0,
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
term_parms=TERM(1.0, 1.0, .5, 1000)
```
As you can see, I start with the defaults. When you make changes, you
write into this structure and it is your job to keep track of what you
did.

To query the parameters type the name of the structure
```
julia> term_parms
TERM(1.00000e+00, 1.00000e+00, 5.00000e-01, 1000)
```
To change Cr from 1 to 40 type
```
julia> update_parms(; Cr=40.0)
TERM(4.00000e+01, 1.00000e+00, 5.00000e-01, 1000)
```
To change Ce to 5 and revert Cr back to the default
```
julia> update_parms(; Ce=5.0)
TERM(1.00000e+00, 5.00000e+00, 5.00000e-01, 1000)

```
Finally, to get the defaults back, enter no changes.
```
julia> update_parms(;)
TERM(1.00000e+00, 1.00000e+00, 5.00000e-01, 1000)

```


"""
function update_parms(t::TERM = term_parms; Cr = 1.0, Ce = 1.0,
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


"""
update_parms(; Cr=Cr_default, Ce=Ce_default, 
      Rmax=Rmax_default, litmax=litmax_default)

C. T. Kelley 2025

Update the termination parameters in MultiPrecisionArrays.

This creates a new TERM data structure that you send to the solver
as a kwarg. 

This changes the values of the parameters. I do not recommend that
you mess with this unless you have a good reason. One such reason
might be that the default limit on the number of iterations (10)
is not working for you.

I start with the default values, so
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

and that is passed to the solvers. So, if AF is a multiprecision
you use the optional argument term_parms.


"""
function update_parms(; Cr=Cr_default, Ce=Ce_default, 
      Rmax=Rmax_default, litmax=litmax_default)
term_parms=TERM(Cr, Ce, Rmax, litmax)
return term_parms
end

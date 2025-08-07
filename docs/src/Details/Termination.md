# Terminating the while loop

Today's values are
```math
C_e = 1.0, C_r = 1.0, R_{max} = .5, litmax = 10
```

Floating point roundoff is 
```math
u_w = 0.5 * eps(TW)
```

We can terminate on a small residual
```math
\| r \| < C_r u_w \| b \|
```
or a small relative normwise backward error
```math
\| r \| < C_e u_w (\| b \| + \| A \| \| x \|)
```

Termination on small residual is the default because computing $\| A \|$
is $N^2$ work and is more expensive that a few IR iterations. I am using
$\| A \|_1$ because that takes less time (data locality) than 
$\| A \|_\infty$ but it is still a pain. I also compute $\| A \|_1$
in ```TF``` if ```TF = Float32``` or ```TF = Float64```. 
Otherwise I use ```TW```. 
LAPACK does not support half precision so that's out.

The problem with these criteria is
that IR can stagnate, especially for ill-conditioned problems, before
the termination criterion is attained. We detect stagnation by looking
for a unacceptable decrease (or increase) in the residual norm. So we will
also terminate the iteration if
```math
\| r_{new} \| \ge R_{max} \| r_{old} \|
```
even if the small residual condition is not satisfied. You can also 
limit the number of IR iterations to manage stagnation. 
Higham [higham97]@cite recomments a limit of ```litmax = 5```. Our default
is ```litmax = 10```, which I may change at any time.

If ```TR > TW``` then I assume you are trying to address extreme
ill-conditioning. In that case I terminate when the norm of the
correction seems to stagnate. 
```math
\| d_{new} \| \ge R_{max} \| d_{old} \|
```

This approach may take one more iteration than the one
recommended in [demmelir](@cite) but is simpler and will 
get you to the theoretical error bould of working precision if
$u_r = u_w^2$.

I am still playing with the termination criteria and the iteration
counts and timings could grow or shrink as I do that. 

We store the parameters in a TERM structure, which define in the main file
for __MultiPrecisionArrays.jl__.
``` 
struct TERM  
       Cr::Real
       Ce::Real
       Rmax::Real
       litmax::Int
end
```
We create one ```TERM``` structure for the defaults
```term_parms_default```. The solvers take a ```TERM``` structure
as a kwarg, with ```term_parms_default```
as the default.


To change the parameters use the __update_parms__ function.
This function makes a TERM structure for you to pass to the solvers.
Herewith the docstrings
```
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
```

Here is an example with the ill-conditioned problem we've been using
in other examples. In this example we decrease Rmax and see that the
solve takes fewer iterations with no change in the residual quality.

```
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> N=512; A=I - 799.0*Gmat(N); AF=mplu(A);

julia> b=ones(N);

julia> mout=\(AF,b; reporting=true);

# Termination on a residual norm increase.

julia> mout.rhist
6-element Vector{Float64}:
 1.00000e+00
 4.39096e-03
 2.85170e-07
 4.30167e-11
 6.05982e-12
 6.34648e-12

julia> term_parms=update_parms(;Rmax=.1);

julia> mout2=\(AF,b; reporting=true, term_parms=term_parms);

# Termination on minimal decrease in the residual norm.

julia> mout2.rhist
5-element Vector{Float64}:
 1.00000e+00
 4.39096e-03
 2.85170e-07
 4.30167e-11
 6.05982e-12
```

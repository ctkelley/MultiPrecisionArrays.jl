# Terminating the while loop

Today's values are
```math
C_e = 1.0, C_r = 20.0, r_{max} = .5, litmax = 1000
```

Floating point roundoff is 
```math
u_r = 0.5 * eps(TR)
```

We can terminate on a small residual
```math
\| r \| < C_r u_r \| b \|
```
or a small relative normwise backward error
```math
\| r \| < C_e u_r (\| b \| + \| A \| \| x \|)
```

Termination on small residual is the default because computing $\| A \|$
is $N^2$ work and is more expensive that a few IR iterations. I am using
$\| A \|_1$ because that takes less time (data locality) than 
$\| A \|_\infty$ but it is still a pain. I also compute $\| A \|_1$
in ```TF``` if ```TF = Float32``` or ```TF = Float64```. 
Otherwise I use ```TW```. 
LAPACK does not support half precision so that's out.

The reason $C_r > C_e$ is
that putting $\| A \|$ in the termination criterion is very useful
and making $C_r$ large helps compensate for that.

The problem with these criteria is
that IR can stagnate, especially for ill-conditioned problems, before
the termination criterion is attained. We detect stagnation by looking
for a unacceptable decrease (or increase) in the residual norm. So we will
also terminate the iteration if
```math
\| r_{new} \| \ge r_{max} \| r_{old} \|
```
even if the small residual condition is not satisfied. You can also 
limit the number of IR iterations to manage stagnation. 
Higham [higham97]@cite recomments a limit of ```litmax = 5```. Our default
is ```litmax = 1000``, which is essentially infinite in this context.

I am still playing with the termination criteria and the iteration
counts and timings could grow or shrink as I do that. 

You can use the __update_parms__ command to
change $C_r$, $C_e$ and $r_{max}$ if you must. I do not advise that.
Anyhow, here are the docstrings.
```
  update_parms(t::TERM = term_parms; Cr = 20.0, Ce = 1.0,
         Rmax = 0.5, litmax=1000

  )

  C. T. Kelley, 2025

  Update the termination parameters in MultiPrecisionArrays.

  This changes the values of the parameters. I do not recommend that you mess
  with this unless you have a good reason. One such reason might be that the
  default limit on the number of iterations (1000, essentially infinite) is
  giving you many iterations that make very little progress. Changing Rmax is
  probably a better way to control this.

  Note that the default values for the kwargs are the defaults, so unless you
  specify otherwise, any parameter will take its default value.

  We store the parameters in a mutable structure TERM

  mutable struct TERM
         Cr::Real
         Ce::Real
         Rmax::Real
         litmax::Int
  end

  where the fields are the parameters. The parameters we use in the solvers
  are in a global TERM structure defined in the main module.

  term_parms=TERM(20.0, 1.0, .5, 1000)

  As you can see, I start with the defaults. When you make changes, you write
  into this structure and it is your job to keep track of what you did.

  To query the parameters type the name of the structure

  julia> term_parms
  MultiPrecisionArrays.TERM(2.00000e+01, 1.00000e+00, 5.00000e-01, 1000)

  To change Cr from 20 to 40 type

  julia> update_parms(; Cr=40.0)
  MultiPrecisionArrays.TERM(4.00000e+01, 1.00000e+00, 5.00000e-01, 1000)

  To change Ce to 5 and revert Cr back to the default

  julia> update_parms(; Ce=5.0)
  MultiPrecisionArrays.TERM(2.00000e+01, 5.00000e+00, 5.00000e-01, 1000)

  Finally, to get the defaults back, enter no changes.

  julia> update_parms(;)
  MultiPrecisionArrays.TERM(2.00000e+01, 1.00000e+00, 5.00000e-01, 1000)
```

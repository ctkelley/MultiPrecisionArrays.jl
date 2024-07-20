# Terminating the while loop

Todays values are
```math
C_e = 1.0, C_r = 20.0, r_{max} = .1
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

The problem with these criteria is
that IR can stagnate, especially for ill-conditioned problems, before
the termination criterion is attained. We detect stagnation by looking
for a unacceptable decrease (or increase) in the residual norm. So we will
also terminate the iteration if
```math
\| r_{new} \| \ge r_{max} \| r_{old} \|
```
even if the small residual condition is not satisfied.

I am still playing with the termination criteria and the iteration
counts and timings could grow or shrink as I do that. 


# Terminating the while loop

Today's values are
```math
C_e = 1.0, C_r = 20.0, r_{max} = .5
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
is $N^2$ work is more expensive that a few IR iterations. I am using
$\| A \|_1$ because that takes less time (data locality) than 
$\| A \|_\infty$ but it is still a pain. The reason $C_r > C_e$ is
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
even if the small residual condition is not satisfied.

I am still playing with the termination criteria and the iteration
counts and timings could grow or shrink as I do that. 


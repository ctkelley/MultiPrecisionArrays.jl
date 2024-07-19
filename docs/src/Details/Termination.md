# Terminating the while loop
Todays values are

```math
C_e = .5, C_r = 1.0
```

We terminate the loop when 
```math
\| r \| < C_e \tau (\| b \| + \| A \| \| x \|)
```
where we use $\tau = 0.5 * eps(TW)$, where $eps(TW)$ is the working
precision floating
point machine epsilon.  The problem with this criterion is
that IR can stagnate, especially for ill-conditioned problems, before
the termination criterion is attained. We detect stagnation by looking
for a unacceptable decrease (or increase) in the residual norm. So we will
terminate the iteration if
```math
\| r_{new} \| \ge .9 \| r_{old} \|
```
even if the small residual condition is not satisfied.

I am still playing with the termination criteria and the iteration
counts and timings could grow or shrink as I do that.


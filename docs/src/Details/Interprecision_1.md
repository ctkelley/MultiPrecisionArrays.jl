# Interprecision Transfers: Part I

The discussion in this section is for the most useful case where
high precision is ```Float64``` and low precision is ```Float32```.
Things are different if low precision is ```Float16```.

Recall that the default way to use the low precision factorization 
is to copy $r$ into low precision, scale it, perform the solve in 
low precision, and then reverse the scaling and promote the 
correction $d$. So if $AF = lu!(A32)$ is 
the factorization object for the low precision factorization, then we
compute $d$ via

```
d = norm(r)* Float64.( AF\ (Float32.(r / norm(r))))
```
We will refer to this approach as the low precision solve (LPS). 
As we said earlier, if one simply does
```
d = AF\r
```
the elements of the triangular matrices are promoted to double as
the solves take place. We will refer to this as a mixed precision
solve (MPS). In the table below we report
timings from Julia's  __BenchmarkTools__ package for double precision
matrix vector multiply (MV64),
single precision LU factorization (LU32) and three approaches
for using the factors to solve a linear system. HPS is the time for
a fully double precision triangular solved and MPS and LPS are the
mixed precision solve and the fully low precision solve.
IR will use a high precision
matrix vector multiply to compute the residual and a solve to
compute the correction for each iteration. The low precision
factorization is done only once.

In this example $A = I + 800 G(N)$ and we look at several values of $N$.

```
    N      MV64       LU32       HPS        MPS        LPS   LU32/MPS
  512    2.8e-05    7.7e-04    5.0e-05    1.0e-04    2.8e-05 7.8e+00
 1024    1.1e-04    2.6e-03    1.9e-04    7.7e-04    1.0e-04 3.4e+00
 2048    6.1e-04    1.4e-02    8.8e-04    3.5e-03    4.0e-04 4.0e+00
 4096    1.9e-03    8.4e-02    4.7e-03    1.4e-02    2.2e-03 5.8e+00
 8192    6.9e-03    5.9e-01    1.9e-02    5.9e-02    9.7e-03 9.9e+00
```

The last column of the table is the ratio of timings for the low precision
factorization and the mixed precision solve. Keeping in mind that at least
two solves will be needed in IR, the table shows that MPS can be
a significant fraction of the cost of the solve for smaller problems and
that LPS is at least 4 times less costly. This is a compelling case
for using LPS in case considered in this section, where high precision
is double and low precision is single, provided the performance of IR
is equally good.

If one is solving $\ma \vx = \vb$ for multiple right hand sides, as one
would do for nonlinear equations in many cases, then
LPS is significantly faster for small and moderately large problems. For
example, for $N=4096$ the cost of MPS is roughly $15\%$ of the low precision
LU factorization, so if one does more than 6 solves with the same
factorization, the solve cost would be more than the factorization cost.
LPS is five times faster and we saw this effect while preparing our
our nonlinear solver package __SIAMFANL.jl__.
The situation for IR is similar, but one must consider
the cost of the high precision matrix-vector multiply, which is about
the same as LPS.

We make LPS the default for IR if high precision is double and low precision
is single. This decision is good for desktop computing. If low precision
is half, then the LPS vs MPS decision needs more scrutiny.

Since MPS does the triangular solves in high precision, one should expect
that the results will be more accurate and that the improved accuracy
might enable the IR loop to terminate earlier \cite{CarsonHigham}.
We should be able to see
that by timing the IR loop after computing the factorization. One should
also verify that the residual norms are equally good.

We will conclude this section with two final tables for the results of IR
with $A = I + \alpha G(N)$. We compare the well
conditioned case ($\alpha=1$) and the ill-conditioned case ($\alpha=800$)
for a few values of $N$. We will look at residual and error norms
for both approaches to interprecision transfer. The conclusion is that
if high precision is double and low is single, the two approaches give
equally good results. 

The columns of the tables are the dimensions, the
$\ell^\infty$ relative error norms for both
LP and MP interprecision transfers (ELP and EMP) and the corresponding
relative residual norms (RLP and RMP).

The results for $\alpha=1$ took 5 IR iterations for all cases. As expected
the LPS iteration was faster than MPS.
However,
for the ill-conditioned $\alpha=800$ case, MPS took one fewer iteration
(5 vs 6)
than EPS
for all but the smallest problem.
Even so, the overall solve
times were essentially the same.

$\alpha=1$
```
    N      ELP        EMP        RLP         RMP        TLP       TMP 
  512    4.4e-16    5.6e-16    3.9e-16    3.9e-16    2.8e-04   3.6e-04 
 1024    6.7e-16    4.4e-16    3.9e-16    3.9e-16    1.2e-03   1.5e-03 
 2048    5.6e-16    4.4e-16    3.9e-16    3.9e-16    5.8e-03   6.2e-03 
 4096    1.1e-15    1.1e-15    7.9e-16    7.9e-16    1.9e-02   2.4e-02 
 8192    8.9e-16    6.7e-16    7.9e-16    5.9e-16    7.0e-02   8.9e-02 
```


$\alpha=800$
```
    N      ELP        EMP        RLP         RMP        TLP       TMP 
  512    6.3e-13    6.2e-13    2.1e-15    1.8e-15    3.3e-04   3.8e-04 
 1024    9.6e-13    1.1e-12    3.4e-15    4.8e-15    1.3e-03   1.4e-03 
 2048    1.0e-12    1.2e-12    5.1e-15    4.5e-15    7.2e-03   6.8e-03 
 4096    2.1e-12    2.1e-12    6.6e-15    7.5e-15    2.4e-02   2.5e-02 
 8192    3.3e-12    3.2e-12    9.0e-15    1.0e-14    8.4e-02   8.9e-02 
```



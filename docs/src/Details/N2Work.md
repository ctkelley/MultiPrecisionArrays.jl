# Is O(N^2) work negligible?

In this sectiion ```TR = TW = Float64``` and
```TF = TS = Float32```, which means that the iterprecision transfers
in the triangular solvers are done in-place. We terminate on small residuals.

The premise behind IR is that reducing the $O(N^3)$ cost of the
factorization will make the solve faster because everything else
is $O(N^2)$ work. It's worth looking to this.

We will use the old-fashioned defintion of a FLOP as an add, a multiply
and a bit of address computation. So we have $N^2$ flops for any of

 - matrix-vector multiply $A x$,
 - the two triangular solves with the LU factors $(LU)^{-1} b$, and
 - computation of the $\ell^1$ or $\ell^\infty$ matrix operator norms $|| A ||_{1,\infty}$.

A linear solve with an LU factorization and the standard triangular
solve has a cost of $(N^3/3) + N^2$ TR-FLOPS. The factorization for IR
has a cost of $N^3/3$ TF-FLOPS or $N^3/6$ TR-FLOPS.

A single IR iteration costs a matrix-vector product in precision TR
and a triangular solve in precision TF for a total of
$3 N^2/2$ TR-FLOPS. Hence a linear solve with IR that needs $n_I$ iterations
costs
```math
\frac{N^3}{6} + 3 n_I N^2/2
```
TR-FLOPS if one terminates on small residuals and an extra $N^2$ TR-FLOPS
if one computes the norm of $\ma$ in precision TR.

IR will clearly be better for large values of $N$. How large is that?
In this example we compare the cost of factorization and solve
using ```lu!``` and ```ldiv``` (cols 2-4) with the equivalent multiprecision
commands $mplu$ and $\backslash$

The operator is ```A = I - Gmat(N)```. We tabulate

 - LU: time for ```AF=lu!(A)```

 - TS: time for ```ldiv!(AF,b)```

 - TOTL = LU+TS

 - MPLU: time for ```MPF=mplu(A)```

 - MPS: time for ```MPF\b```

 - TOT: MPLU+MPS

 - OPNORM: Cost for $\| A \|_1$, which one needs to terminate on small normwise backward error.

The message from the tables is that the triangular solves are more costly
than operation counts might indicate. One reason for this is that the
LU factorization exploits multi-core computing better than a triangular
solve. It is also interesting to see how the choice of BLAS affects the
results and how the cost of the operator norm of the matrix is more than
a triangular solve.

For both cases, multiprecision arrays perform better when $N \ge 2048$
with the difference becoming larger as $N$ increases.

openBLAS
```
   N        LU        TS       TOTL      MPLU       MPS       TOT    OPNORM  
   512  9.87e-04  5.35e-05  1.04e-03  7.74e-04  2.90e-04  1.06e-03  1.65e-04 
  1024  3.67e-03  2.14e-04  3.88e-03  3.21e-03  8.83e-04  4.10e-03  7.80e-04 
  2048  2.10e-02  1.20e-03  2.22e-02  1.54e-02  4.72e-03  2.02e-02  3.36e-03 
  4096  1.46e-01  5.38e-03  1.51e-01  8.93e-02  1.91e-02  1.08e-01  1.43e-02 
  8192  1.10e+00  2.01e-02  1.12e+00  5.98e-01  6.89e-02  6.67e-01  5.83e-02 
```

AppleAcclerateBLAS
```
   N        LU        TS       TOTL      MPLU       MPS       TOT    OPNORM  
   512  7.86e-04  1.10e-04  8.96e-04  4.76e-04  3.58e-04  8.35e-04  1.65e-04 
  1024  3.28e-03  4.54e-04  3.73e-03  2.44e-03  1.30e-03  3.74e-03  7.80e-04 
  2048  2.27e-02  2.46e-03  2.52e-02  1.32e-02  1.13e-02  2.45e-02  3.36e-03 
  4096  1.52e-01  1.13e-02  1.63e-01  6.51e-02  5.62e-02  1.21e-01  1.40e-02 
  8192  1.17e+00  5.35e-02  1.23e+00  6.27e-01  2.48e-01  8.75e-01  5.67e-02 
```



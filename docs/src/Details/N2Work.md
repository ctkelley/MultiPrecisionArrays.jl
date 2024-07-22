# Is O(N^2) work negligible?

In this example we compare the cost of factorization and solve
using $lu!$ and $ldiv$ (cols 2-4) with the equivalent multiprecision
commands $mplu$ and $\\$
$A = I - G(N)$ 

openBLAS
```
   N        LU        TS       TOT      MPLU       MPS       TOT    OPNORM  
   512  9.87e-04  5.35e-05  1.04e-03  7.74e-04  2.90e-04  1.06e-03  1.65e-04 
  1024  3.67e-03  2.14e-04  3.88e-03  3.21e-03  8.83e-04  4.10e-03  7.80e-04 
  2048  2.10e-02  1.20e-03  2.22e-02  1.54e-02  4.72e-03  2.02e-02  3.36e-03 
  4096  1.46e-01  5.38e-03  1.51e-01  8.93e-02  1.91e-02  1.08e-01  1.43e-02 
  8192  1.10e+00  2.01e-02  1.12e+00  5.98e-01  6.89e-02  6.67e-01  5.83e-02 
```

AppleAcclerateBLAS



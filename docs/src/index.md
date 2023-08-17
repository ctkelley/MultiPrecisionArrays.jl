# MultiPrecisionArrays.jl v0.0.4

[C. T. Kelley](https://ctk.math.ncsu.edu)

[MultiPrecisionArrays.jl](https://github.com/ctkelley/MultiPrecisionArrays.jl) is a package for iterative refinement. 

This is mostly for research at this point. It will get more interesting as hardware support for Float16 comes online.

This package will make solving dense systems of linear equations faster by using the LU factorization and iterative refinement. It is limited to LU for now. A very generic description of this for solving a linear systrem $A x = b$ is

__IR(A, b, x)__
- $r = b - Ax$
- Factor $A = LU$ in a lower precision
- While $\|| r \||$ is too large
  - $d = (LU)^{-1} r$
  - $x = x + d$
  - $r = b - Ax$
- end

In Julia, a code to do this would solve the linear system $A x = b$ in double precision by using a
factorization in a lower precision, say single, within a residual correction iteration. This means that one would need
to allocate storage for a copy of $A$ is the lower precision and factor that copy. The one has to determine what the line
$d = (LU)^{-1} r$ means. Do you case $r$ into the lower precison before the solve or not? __MultiPrecisionArrays.jl__ provides
data structures and solvers to manage this. The __MPArray__ structure lets you preallocate $A$ and the low precision copy.



# MultiPrecisionArrays.jl v0.0.4

[C. T. Kelley](https://ctk.math.ncsu.edu)

[MultiPrecisionArrays.jl](https://github.com/ctkelley/MultiPrecisionArrays.jl) is a package for iterative refinement. 

This package provides data atructures and solvers for several variants of iterative refinement. It will become much more useful when half precision (aka ```Float16```) is fully supported in LAPACK/BLAS. For now, it's only general-purpose
application is classical iterative refinement with double precision equations and single precision factorizations.

The half precision stuff is good for those of us doing research in this field. Half precision performace has progressed to the point where you can acutally get things done. On an Apple M2-Pro, a half precision LU only costs 3--5 times
what a double precision LU costs. This may be as good as it gets unless someone wants to duplicate the LAPACK implementation and get the benefits from blocking, recursion, and clever cache management.

We use a hack-job LU factorization for half precision. Look at the source
for __hlu!.jl__.

## What is iterative refinement.

The idea is to solve $Ax=b$ in high precision using a factorization
in lower precision. 

__IR(A, b)__

- $x = 0$
- $r = b$
- Factor $A = LU$ in a lower precision
- While $\| r \|$ is too large
  - $d = (LU)^{-1} r$
  - $x = x + d$
  - $r = b - Ax$
- end

In Julia, a code to do this would solve the linear system $A x = b$ in double precision by using a
factorization in a lower precision, say single, within a residual correction iteration. This means that one would need
to allocate storage for a copy of $A$ is the lower precision and factor that copy. The one has to determine what the line
$d = (LU)^{-1} r$ means. Do you case $r$ into the lower precison before the solve or not? __MultiPrecisionArrays.jl__ provides
data structures and solvers to manage this. The __MPArray__ structure lets you preallocate $A$ and the low precision copy.

The submodule __MultiPrecisionArrays.Examples__ has an example which we will use for most of the documentation. The function ```Gmat(N)``` returns
the trapeziod rule discretization of the Greens operator for $-d^2/dx^2$ on $[0,1]$


```math 
G u(x) = \int_0^1 g(x,y) u(y) \, dy 
```


where


```math 
g(x,y) = \left( \begin{array}{l} 
y (1 - x) \ \mbox{if x > y} \\
x (1 -y ) \ \mbox{otherwise}
\end{array}
\right.
```

LaTeXEquation(raw"""
    \left[\begin{array}{c}
        x \\
        y
    \end{array}\right]
""")

- $g(x,y) = y (1 - x)$  if x > y
- $g(x,y) = x (1 -y )$   otherwise

The code for this is in the __/src/Examples__ directory. The file is __Gmat.jl__.



# MultiPrecisionArrays.jl v0.0.5

[C. T. Kelley](https://ctk.math.ncsu.edu)

[MultiPrecisionArrays.jl](https://github.com/ctkelley/MultiPrecisionArrays.jl) is a package for iterative refinement. 

This package provides data structures and solvers for several variants of iterative refinement. It will become much more useful when half precision (aka ```Float16```) is fully supported in LAPACK/BLAS. For now, it's only general-purpose
application is classical iterative refinement with double precision equations and single precision factorizations.

The half precision stuff is good for those of us doing research in this field. Half precision performance has progressed to the point where you can actually get things done. On an Apple M2-Pro, a half precision LU only costs 3--5 times
what a double precision LU costs. This may be as good as it gets unless someone wants to duplicate the LAPACK implementation and get the benefits from blocking, recursion, and clever cache management.

We use a hack-job LU factorization for half precision. Look at the source
for __hlu!.jl__.

## What is iterative refinement.

The idea is to solve $Ax=b$ in high precision using a factorization
in lower precision. 

__IR(A, b)__

- Initialize: $x = 0$,  $r = b$
- Factor $A = LU$ in a lower precision
- While $\| r \|$ is too large
  - Compute the defect $d = (LU)^{-1} r$
  - Correct the solution $x = x + d$
  - Update the residual $r = b - Ax$
- end

In Julia, a code to do this would solve the linear system $A x = b$ in double precision by using a
factorization in a lower precision, say single, within a residual correction iteration. This means that one would need
to allocate storage for a copy of $A$ is the lower precision and factor that copy. 

Then one has to determine what the line
$d = (LU)^{-1} r$ means. Do you cast $r$ into the lower precision before the solve or not? __MultiPrecisionArrays.jl__ provides
data structures and solvers to manage this. 

Here's a simple Julia function for IR that
```
"""
IR(A,b)
Simple minded iterative refinement
Solve Ax=b
"""
function IR(A, b)
    x = zeros(length(b))
    r = copy(b)
    tol = 100.0 * eps(Float64)
    #
    # Allocate a single precision copy of A and factor in place
    #
    A32 = Float32.(A)
    AF = lu!(A32)
    #
    # Give IR at most ten iterations, which it should not need
    # in this case
    #
    itcount = 0
    # The while loop will get more attention later.
    while (norm(r) > tol * norm(b)) && (itcount < 10)
        #
        # Store r and d = AF\r in the same place.
        #
        ldiv!(AF, r)
        x .+= r
        r .= b - A * x
        itcount += 1
    end
    return x
end
```

The __MPArray__ structure contains both $A$, the low precision copy,
and a vector for the residual. 
This lets you allocate the data in advance and reuse the structure
for other right hand sides without rebuilding (or refactoring!) the
low precision copy. 

As written in the function, the defect uses ```ldiv!``` to compute
```AF\r```. This means that the two triangular factors are stored in
single precision and interprecision transfers are done with each
step in the factorization. While that ```on the fly``` interprecision 
transfer is an option, and is needed in many situations, the
default is to downcast $r$ to low precision, do the solve entirely in
low precision, and the upcast the result. The code for that looks like
```
normr=norm(r)
ds=Float32.(r)/normr
ldiv!(AF, ds)
r .= Float64.(ds)*normr
```
The scaling by ```1.0/normr``` avoids underflow, which is most important
when the low precision is ```Float16```. We will discuss interprecision 
transfer costs later.

## Integral Equations Example

The submodule __MultiPrecisionArrays.Examples__ has an example which we will 
use for most of the documentation. The function ```Gmat(N)``` returns
the trapeziod rule discretization of the Greens operator 
for $-d^2/dx^2$ on $[0,1]$ with homogeneous Dirichlet boundary conditions.

```math 
G u(x) = \int_0^1 g(x,y) u(y) \, dy 
```

where


```math
g(x,y) = 
    \left\{\begin{array}{c}
        y (1-x) ; \ x > y\\
        x (1-y) ; \ x \le y
    \end{array}\right.
```

The eigenvalues of $G$ are $1/(n^2 \pi^2)$ for $n = 1, 2, \dots$.

The code for this is in the __/src/Examples__ directory. 
The file is __Gmat.jl__.

In the example we will build a matrix $A = I - \alpha G$. In the examples
we will use $\alpha=1.0$, a very well conditioned case, and $\alpha=800.0$
This latter case is very near singularity.

We will solve a linear system with both double precision $LU$ and an MPArray
and compare execution time and the quality of the results.
The problem setup is pretty simple
```
julia> using MultiPrecisionArrays

julia> using BenchmarkTools

julia> using MultiPrecisionArrays.Examples

julia> N=4096; G=Gmat(N); A=I - G; x=ones(N); b=A*x;

julia> @belapsed lu!(AC) setup=(AC=copy($A))
1.43148e-01
```
At this point we have timed ```lu!```. The next step is to construct
an MPArray and factor the low precision matrix. We use the
constructor ```MPArray``` to store $A$, the low precision copy
and the residual. The we apply
the function ```mplu!``` to factor the low precision copy in place.
```
julia> MPA=MPArray(A);

julia> @belapsed mplu!(MPAC) setup=(MPAC=deepcopy($MPA))
8.02158e-02
```
So the single precision factorization is roughly half the cost of the
double precision one. Now for the solves. Both ```lu!``` and ```mplu!```
produce a factorization object and ```\``` works with both.
You have to be a bit careful because MPA and A share storage. So
I will use ```lu``` instead of ```lu!``` when factoring $A$.
```
julia> AF=lu(A); xf = AF\b;

julia> MPAF=mplu!(MPA); xmp=MPAF\b;

julia> println(norm(xf-x,Inf),"  ",norm(xmp-x,Inf))
7.41629e-14  8.88178e-16
```
You can see that the solutions are equally good.



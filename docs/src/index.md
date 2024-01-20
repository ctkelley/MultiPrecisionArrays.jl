# MultiPrecisionArrays.jl v0.0.10

[C. T. Kelley](https://ctk.math.ncsu.edu)

[MultiPrecisionArrays.jl](https://github.com/ctkelley/MultiPrecisionArrays.jl) 
[ctk:mparrays](@cite)
is a package for iterative refinement. 

These docs are enough to get you started. The complete version with
a better account of the theory is [ctk:mparraysdocs](@cite). 

This package uses __SIAMFANLEquations.jl__ [ctk:siamfanl](@cite), the
solver package associated with a book [ctk:fajulia](@cite) and
suite of IJulia notebooks [ctk:notebooknl](@cite).

This package provides data structures and solvers for several variants of iterative refinement (IR). It will become much more useful when half precision (aka ```Float16```) is fully supported in LAPACK/BLAS. For now, it's only general-purpose
application is classical iterative refinement with double precision equations and single precision factorizations.

The half precision stuff is good for those of us doing research in this field. Half precision performance has progressed to the point where you can actually get things done. On an Apple M2-Pro, a half precision LU only costs 3--5 times
what a double precision LU costs. This may be as good as it gets unless someone wants to duplicate the LAPACK implementation and get the benefits from blocking, recursion, and clever cache management.

We use a hack-job LU factorization for half precision. Look at the source
for __hlu!.jl__.

## What is iterative refinement.

The idea is to solve $Ax=b$ in high precision using a factorization
in lower precision. 

IR is a perfect example of a storage/time tradeoff. To solve a linear system  
$Ax=b$ with IR, one incurs the storage penalty of making a low 
precision copy of  and reaps the benefit of only having to factor the 
low precision copy.

Here is the textbook 
version [higham](@cite) using the LU factorization.

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

The ```MPArray``` structure contains both $A$, the low precision copy,
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

The example below compares the cost of a double precision factorization to a MPArray factorization. The ```MPArray``` structure has a high precision and a low precision matrix. The structure we will start with 
is
```
struct MPArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::AbstractArray{TH,2}
    AL::AbstractArray{TL,2}
    residual::Vector{TH}
    onthefly::Bool
end
```
The structure also stores the residual. The ```onthefly``` Boolean tells the solver how to do the interprecision transfers. The easy way to get started is to use the ```mplu``` 
command directly on the matrix. That will build the MPArray, follow that with the factorization, and put in all in a structure
that you can use with ```\```.

Now we will see how the results look. In this example we compare the result with iterative refinement with ```A\b```, which is LAPACK's LU. 
As you can see the results are equally good. Note that the factorization object ```MPF``` is the
output of ```mplu```. This is analogous to ```AF=lu(A)``` in LAPACK.

Now we will see how the results look. In this example we compare the result with iterative refinement with ```A\b```, which is LAPACK's LU. 
As you can see the results are equally good. Note that the factorization object ```MPF``` is the
output of ```mplu```. This is analogous to ```AF=lu(A)``` in LAPACK.

```
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> using BenchmarkTools

julia> N=4096;

julia> G=Gmat(N);

julia> A = I - G;

julia> MPF=mplu(A); AF=lu(A);

julia> z=MPF\b; w=AF\b;

julia> ze=norm(z-x,Inf); zr=norm(b-A*z,Inf)/norm(b,Inf);

julia> we=norm(w-x,Inf); wr=norm(b-A*w,Inf)/norm(b,Inf);

julia> println("Errors: $ze, $we. Residuals: $zr, $wr")
Errors: 8.88178e-16, 7.41629e-14. Residuals: 1.33243e-15, 7.40609e-14

```

So the resuts are equally good.

The compute time for ```mplu``` should be a bit more than half that of ```lu!```. The reason is
that ```mplu``` factors a low precision array, so the factorization cost is cut in half. Memory
is a different story because. The reason
is that both ```mplu``` and ```lu!``` do not allocate storage for a new high precision array,
but ```mplu``` allocats for a low precision copy, so the memory and allocation cost for ```mplu```
is 50% more than ```lu```. 

```
julia> @belapsed mplu($A)
8.60945e-02

julia> @belapsed lu!(AC) setup=(AC=copy($A))
1.42840e-01

```
It is no surprise that the factorization in single precision took roughly half as long as the one in double. In the double-single precision case, iterative refinement is a great
expample of a time/storage tradeoff. You have to store a low precision copy of $A$, so the storage burden increases by 50\% and the factoriztion time is cut in half.
The advantages of IR increase as the dimension increases. IR is less impressive for smaller problems and can even be slower
```
julia> N=30; A=I + Gmat(N); 

julia> @belapsed mplu($A)
4.19643e-06

julia> @belapsed lu!(AC) setup=(AC=copy($A))
3.70825e-06
```

## Options and data structures for mplu

Here is the source for ```mplu```
```
"""
mplu(A::AbstractArray{Float64,2}; TL=Float32, onthefly=false)

Combines the constructor of the multiprecision array with the
factorization.
"""
function mplu(A::AbstractArray{TH,2}; TL=Float32, onthefly=nothing) where TH <: Real
#
# If the high precision matrix is single, the low precision must be half.
#
(TH == Float32) && (TL = Float16)
#
# Unless you tell me otherwise, onthefly is true if low precision is half
# and false if low precision is single.
#
(onthefly == nothing ) && (onthefly = (TL==Float16))
MPA=MPArray(A; TL=TL, onthefly=onthefly)
MPF=mplu!(MPA)
return MPF
end
```

The function ```mplu``` has two keyword arguments. The easy one to understand is ```TL``` which is the precision of the factoriztion. Julia has support for single (```Float32```) and half (```Float16```)
precisions. If you set ```TL=Float16``` then low precision will be half. Don't do that unless you know what you're doing. Using half precision is a fast way to get incorrect results. Look at the section on [half precision](#half-Precision) in this Readme for a bit more bad news.

The other keyword arguemnt is __onthefly__. That keyword controls how the triangular solvers from the factorization work. When you solve

```math
LU d = r
```

The LU factors are in low precision and the residual $r$ is in high precision. If you let Julia and LAPACK figure out what to do, then the solves will be done in high precision and
the entries in the LU factors will be comverted to high precision with each binary operation. The output $d$ will be in high precision. This is called interprecision transfer on-the-fly
and ```onthefly = true``` will tell the solvers to do it that way. You have $N^2$ interprecsion transfers with each solve and, as we will see, that can have a non-trivial cost.

When low precision is Float32, then the default is (```onthefly = false```). This converts $r$ to low precision, does the solve entirely in low precision, and then promotes $d$ to high precision. You need to be careful to avoid
overflow and, more importantly, underflow when you do that and we scale $r$ to be a unit vector before conversion to low precisiion and reverse the scaling when we promote $d$. We take care of this for you.

```mplu``` calls the constructor for the multiprecision array and then factors the low precision matrix. In some cases, such as nonlinear solvers, you will want to separate the constructor and the factorization. When you do that
remember that ```mplu!``` overwrites the low precision copy of A 
with the factors. The factorizaton object is different from the multiprecision
array, even though they share storage. This is just like ```lu!```.

## Memory Allocations for mplu

The memory footprint of a multiprecision array is dominated by
the high precision array and the low precision copy. The allocations of
```
AF1=lu(A)
```
and
```
AF2=mplu(A)
```
are very different. Typically ```lu``` makes a high precision copy of $\ma$ and
factors that with ```lu!```. ```mplu``` on the other hand, uses $\ma$
as the high precision matrix in the multiprecision array structure and
the makes a low precision copy to send to ```lu!```. Hence ```mplu```
has half the allocation burden of ```lu```.

That is, of course misleading. The best way to apply ```lu``` is to
overwrite $A$ with the factorization using
```
AF1=lu!(A).
```
The analog of this approach with a multiprecision array would be to
first build an ```MPArray``` structure with
```
MPA = MPArray(A)
```
which makes $A$ the high precision matrix and also makes a low
precision copy. This is the stage where the extra memory is allocated
for the the low precision copy. One follows that with the factorization
of the low precision matrix to construct the factorization object.
```
MPF = mplu!(MPA).
```
The function ```mplu``` simply applies ```MPArray``` and follows that with
```mplu!```.

## Other IR software in Julia

The package
[IterativeRefinement.jl](https://github.com/RalphAS/IterativeRefinement.jl)
is an implementation of the IR method from
[dongarra:1983](@cite).

The unregistered
package [Itref.jl](https://github.com/bvieuble/Itref.jl) implements
IR and the GMRES-IR method from 
[amestoy:2023](@cite)
and was used to obtain
the numerical results in that paper. It does not provide the
data structures for preallocation that we do and does not seem to have
been updated lately.




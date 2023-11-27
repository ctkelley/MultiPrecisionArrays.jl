[![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url]
[![][build-status-img]][build-status-url]
[![][codecov-img]][codecov-url]
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7521427.svg)](https://doi.org/10.5281/zenodo.7521427)

# MultiPrecisionArrays.jl v0.0.8

## [C. T. Kelley](https://ctk.math.ncsu.edu)

This is the start of a package to support multiprecision arrays. This is mostly for my own research right now.

This package provides data structures and solvers for several variants of iterative refinement (IR). It will become much more useful when half precision (aka ```Float16```) is fully supported in LAPACK/BLAS. For now, its only general-purpose
application is classical iterative refinement with double precision equations and single precision factorizations. 

The half precision stuff is good for those of us doing research in this field. Half precision performace has progressed to the point where you can actually get things done. On an Apple M2-Pro, a half precision LU only costs 3--5 times
what a double precision LU costs. This may be as good as it gets unless someone wants to duplicate the LAPACK implementation and get the benefits from blocking, recursion, and clever cache management.

__The half precision LU is much faster (more than 10x) as of v0.0.3. Look at [src/Factorizations/hlu!.jl](src/Factorizations/hlu!.jl) to see the hack job I did to [generic_lu](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/lu.jl).__

__The ArXiv paper I just submitted was for v0.0.7.__

## What parts of MultiPrecisionArrays.jl are ready for prime time?

- Using ```mplu``` to cut your solve time in half for double precision matrices is simple to use and works well.
- The API for harvesting iteration statistics is stable.
- If you're a half precision person, the new factorizatino ```hlu!``` is __only__ 3-5x slower than a double precision ```lu```.

## I found this repo by accident on the internet. Can I complain about it now?

__I have registered this pacakge so a few people who work in this field can play with it. I'll announce it on disocourse when I have the API more stable than it is now.__

__Please do not make PRs. If you stumble on this mess and have questions/ideas ..., raise an issue or email me at tim_kelley@ncsu.edu__

Nothing is in final form and I am changing the API, internal structures, exported functions/structs and all kinds of other stuff frequently. When/if I announce this on Discourse, then it will be time for complaints.

## Coming Attractions

- Nonlinear solver applications
- Three precision IR
- BiCGSTAB-IR

## Documentation

- This README file
- The [docs for the package](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev)
- [The ArXiv paper for v0.0.7](https://arxiv.org/abs/2311.14616)
- [The updated ArXiv paper for v0.0.8](https://github.com/ctkelley/MultiPrecisionArrays.jl/blob/main/Publications_and_Presentations/MPArray.pdf)


## Readme Contents:

- [Algorithms](#algorithms)
- [Example](#example)
  - [Subtleties in the example](#a-few-subtleties-in-the-example)
- [Be Careful with Half Precision](#half-precision)
- [Harvesting Iteration Statistics](#harvesting-iteration-statistics)
- [Dependencies](#dependencies)
- [Endorsement](#endorsement)
- [Funding](#funding)


## Algorithms

### What is iterative refinement?

This package will make solving dense systems of linear equations faster by using the LU factorization and IR. It is limited to LU for now. A very generic description of this for solving a linear systrem $A x = b$ is

__IR(A, b)__
- $x = 0$
- $r = b$
- Factor $A = LU$ in a lower precision
- While $\|| r \||$ is too large
  - $d = (LU)^{-1} r$
  - $x = x + d$
  - $r = b - Ax$
- end

In Julia, a code to do this would solve the linear system $A x = b$ in double precision by using a
factorization in a lower precision, say single, within a residual correction iteration. This means that one would need
to allocate storage for a copy of $A$ is the lower precision and factor that copy. The one has to determine what the line
$d = (LU)^{-1} r$ means. Do you cast $r$ into the lower precison before the solve or not? __MultiPrecisionArrays.jl__ provides
data structures and solvers to manage this. The __MPArray__ structure lets you preallocate $A$, the low precision copy, and the residual $r$.

IR is a perfect example of a storage/time tradeoff.
To solve a linear system $A x = b$ in $R^N$ with IR,
one incurs the storage penalty of making a low
precision copy of $A$ and reaps the benefit of only having to
factor the low precision copy.


### Termination of the IR loop.

We terminate the while loop when 

$$
\\| r \\| < \tau \\| b \\|
$$

where we use $\tau = 10 * {\mbox eps}(TH)$. Here $\mbox{eps}(TH)$
is high precision
machine epsilon.  The problem with this criterion is
that IR can stagnate, especially for ill-conditioned problems, before
the termination criterion is attained. We detect stagnation by looking
for a unacceptable decrease (or increase) in the residual norm. So we will
terminate the iteration if

$$
\\| r_{new} \\| \ge .9 \\| vr_{old} \\|
$$

even if the residual decrease criterion is not satisfied.




## Example

### An example to get started

Herewith, the world's most simple example to show how iterative refienment works. We will follow that with some benchmarking on the cost of factorizations.
The functions we use are __MPArray__ to create the structure and __mplu!__ to factor the low precision copy. In this example high precision is ```Float64``` and low
precision is ```Float32```. The matrix is the sum of the identity and a constant multiple of the trapezoid rule discretization of the Greens operator for $-d^2/dx^2$ on $[0,1]$

$$
G u(x) = \int_0^1 g(x,y) u(y) \, dy 
$$

where

$$
g(x,y) = \\left\\{ \begin{array}{l} 
y (1 - x) \ \mbox{if x > y} \\
x (1 -y ) \ \mbox{otherwise}
\end{array}
\right.
$$

The code for this is in the __/src/Examples__ directory. The file is __Gmat.jl__. You need to do 
```
using MultiPrecisionArrays.Examples
```
to get to it.

The example below compares the cost of a double precision factorization to a MPArray factorization. The ```MPArray``` structure has a high precision and a low precision matrix. The structure we will start with 
is
```
struct MPArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::Array{TH,2}
    AL::Array{TL,2}
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

The compute time for ```mplu``` should be half that of ```lu```.



```
julia> @belapsed mplu($A)
8.55328e-02

julia> @belapsed lu($A)
1.49645e-01

```

It is no surprise that the factorization in single precision took roughly half as long as the one in double. In the double-single precision case, iterative refinement is a great
expample of a time/storage tradeoff. You have to store a low precision copy of $A$, so the storage burden increases by 50\% and the factoriztion time is cut in half.

### A few subtleties in the example

Here is the source for ```mplu```
```
"""
mplu(A::Array{Float64,2}; TL=Float32, onthefly=false)

Combines the constructor of the multiprecision array with the
factorization.
"""
function mplu(A::Array{TH,2}; TL=Float32, onthefly=nothing) where TH <: Real
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

The other keyword arguement is __onthefly__. That keyword controls how the triangular solvers from the factorization work. When you solve

$$ 
LU d = r
$$

The LU factors are in low precision and the residual $r$ is in high precision. If you let Julia and LAPACK figure out what to do, then the solves will be done in high precision and
the entries in the LU factors will be comverted to high precision with each binary operation. The output $d$ will be in high precision. This is called interprecision transfer on-the-fly
and ```onthefly = true``` will tell the solvers to do it that way. You have $N^2$ interprecsion transfers with each solve and, as we will see, that can have a non-trivial cost.

When low precision is Float32, then the default is (```onthefly = false```). This converts $r$ to low precision, does the solve entirely in low precision, and then promotes $d$ to high precision. You need to be careful to avoid
overflow and, more importantly, underflow when you do that and we scale $r$ to be a unit vector before conversion to low precisiion and reverse the scaling when we promote $d$. We take care of this for you.

```mplu``` calls the constructor for the multiprecision array and then factors the low precision matrix. In some cases, such as nonlinear solvers, you will want to separate the constructor and the factorization. When you do that
remember that ```mplu!``` overwrites the low precision copy of A with the factors, so you can't resuse the multiprecision array for other problems unless you restore the low precision copy.


__MultiPrecisionArrays.jl__ supports many variations of iterative refinement and we explain all that in the docs and a [users guide](https://github.com/ctkelley/MultiPrecisionArrays.jl/blob/main/Publications_and_Presentations/MPArray.pdf).

## Half Precision

Bottom line: don't use it and expect good performance. See [this page](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Half_1/) in the docs for details.

It's a good idea to try GMRES-IR if you are playing with half precision. GMRES-IR uses a different factorization ```mpglu``` which factors the 
low precision matrix, allocates room for the Krylov basis, and builds a factorization object. GMRES-IR uses the low precision factorization as
a preconditioner for GMRES. The call looks like ```MPGF=mpglu(A)```. You solve $A x = b$ with
```x = MPGF\b```. More details [here](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Half_1/#GMRES-IR).

## Harvesting Iteration Statistics

You can get some iteration statistics by using the
__reporting__ keyword argument to the solvers. The easiest way
to do this is with the backslash command. When you use this option you
get a data structure with the solution and the residual history.

```
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> N=4096; A = I - Gmat(N); x=ones(N); b=A*x;

julia> MPF=mplu(A);

julia> # Use \ with reporting=true

julia> mpout=\(MPF, b; reporting=true);

julia> norm(b-A*mpout.sol, Inf)
1.33227e-15

julia> # Now look at the residual history

julia> mpout.rhist
5-element Vector{Float64}:
 9.99878e-01
 1.21892e-04
 5.25805e-11
 2.56462e-14
 1.33227e-15
```
As you can see, IR does well for this problem. The package uses an initial
iterate of $x = 0$ and so the initial residual is simply $r = b$
and the first entry in the residual history is $|| b ||_\infty$. The
iteration terminates successfully after four matrix-vector products.

There are more examples for this in the [paper](https://github.com/ctkelley/MultiPrecisionArrays.jl/blob/main/Publications_and_Presentations/MPArray.pdf).



## Dependencies

As of now you need to install these packages

- Polyester
- SIAMFANLEquations

and use LinearAlgebra and SparseArrays from Base. I use the Krylov solvers and examples from SIAMFANLEquations. Polyester is here because 
threading with ```Polyester.@batch``` makes the LU for Float16 run much faster. Once LU for Float16 is in LAPACK/BLAS, I will eliminate that dependency.

## Endorsement

I have used this in my own work and will be adding links to that stuff as I finish it. 

I started on this package after finishing

(KEL22a) C. T. Kelley, [__Newton's Method in Mixed Precision__](https://epubs.siam.org/doi/10.1137/20M1342902), SIAM Review 35 (1998), pp 191-211. 

(KEL22b) C. T. Kelley, [__Solving Nonlinear Equations with Iterative Methods: Solvers and Examples in Julia__](https://my.siam.org/Store/Product/viewproduct/?ProductId=44313635), SIAM, Philadelphia, 2022. 

A new paper exploits most of the package. I use the examples in that paper for CI. If you do iterative refinement well, you can make half precision work far better than it did in (KEL22a). MultiPrecisionArrays drop right into the solvers in [SIAMFANLEquations.jl](https://github.com/ctkelley/SIAMFANLEquations.jl). 

(KEL23a) C. T. Kelley, [__Newton's method in three precisions__](https://arxiv.org/abs/2307.16051) To appear in Pacific Journal of Optimization.

This paper has a [repo](https://github.com/ctkelley/N3Presults) for reproducing the results with an early version of this package.


## __Interprecision transfers__:

This is not a trivial matter. I gave a talk at the XSDK-MULTIPRECISION meeting on June 15, 2023, about this issue. See [the docs](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Details/Interprecision_1/) for the story on this. 
Most users of this package can ignore this issue.

  - [Interprecision Transfers in Iterative Refinement: Making Half Precision on Desktops Less Painful](Publications_and_Presentations/MPArrays_XSDK-MULTIPRECISION_June_15.pdf).



## Funding

This project was partially supported by

1.  [National Science Foundation Grant DMS-1906446](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1906446&HistoricalAwards=false) and
2. [Department of Energy grant DE-NA003967](https://cement-psaap.github.io)

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://ctkelley.github.io/MultiPrecisionArrays.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ctkelley.github.io/MultiPrecisionArrays.jl/dev

[build-status-img]: https://github.com/ctkelley/MultiPrecisionArrays.jl/workflows/CI/badge.svg
[build-status-url]: https://github.com/ctkelley/MultiPrecisionArrays.jl/actions

[codecov-img]: https://codecov.io/gh/ctkelley/MultiPrecisionArrays.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/ctkelley/MultiPrecisionArrays.jl


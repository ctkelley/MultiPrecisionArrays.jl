[![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url]
[![][build-status-img]][build-status-url]
[![][codecov-img]][codecov-url]
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7521427.svg)](https://doi.org/10.5281/zenodo.7521427)

# MultiPrecisionArrays.jl v0.0.4

## [C. T. Kelley](https://ctk.math.ncsu.edu)

This is the start of a package to support multiprecision arrays. This is for my own research right now.

This package provides data structures and solvers for several variants of iterative refinement. It will become much more useful when half precision (aka ```Float16```) is fully supported in LAPACK/BLAS. For now, its only general-purpose
application is classical iterative refinement with double precision equations and single precision factorizations. 

The half precision stuff is good for those of us doing research in this field. Half precision performace has progressed to the point where you can acutally get things done. On an Apple M2-Pro, a half precision LU only costs 3--5 times
what a double precision LU costs. This may be as good as it gets unless someone wants to duplicate the LAPACK implementation and get the benefits from blocking, recursion, and clever cache management.

__The half precision LU is much faster (more than 10x) as of v0.0.3. Look at [src/Factorizations/hlu.jl](src/Factorizations/hlu.jl) to see the hack job I did to [generic_lu](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/lu.jl).__


## Is MultiPrecisionArrays.jl ready for prime time?

No

__This package is nowhere close to ready for registration or release. It's public only to help me do CI and clean up the docs.__

__Please do not make PRs. If you stumble on this mess and have questions/ideas ..., raise an issue or email me at tim_kelley@ncsu.edu__

Nothing is in final form and I am changing the API, internal structures, exported functions/structs and all kinds of other stuff frequently. When/if I register this and announce it, then it will be time for complaints and offers to collaborate. 



## Readme Contents:
- [Algorithms](#algorithms)
- [Example](#example)
  - [Subtleties in the example](#a-few-subtleties-in-the-example)
- [Be Careful with Half Precision](#half-precision)    
- [Dependencies](#dependencies)
- [Endorsement](#endorsement)
- [Funding](#funding)

## Algorithms

### What is iterative refinement?

This package will make solving dense systems of linear equations faster by using the LU factorization and iterative refinement. It is limited to LU for now. A very generic description of this for solving a linear systrem $A x = b$ is

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
data structures and solvers to manage this. The __MPArray__ structure lets you preallocate $A$ and the low precision copy. 

## Example

### An example to get started

Herewith, the world's most simple example to show how iterative refienment works. We will follow that with some benchmarking on the cost of factorizations.
The functions we use are __MPArray__ to create the structure and __mplu!__ to factor the low precision copy. In this example high precision is ```Float64``` and low
precision is ```Float16```. The matrix is the sum of the identity and a constant multiple of the trapezoid rule discretization of the Greens operator for $-d^2/dx^2$ on $[0,1]$


$$
G u(x) = \int_0^1 g(x,y) u(y) \, dy 
$$

where

$$
g(x,y) = \left( \begin{array}{l} 
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

The example below compares the cost of a double precision factorization to a MPArray factorization. The ```MPArray``` structure has high precision and low precision matrix. The structure we will start with 
is
```
struct MPArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::Array{TH,2}
    AL::Array{TL,2}
end
```

We will be using the constructor
```
function MPArray(AH::Array{Float64,2}; TL = Float32, onthefly=false)
    AL = TL.(AH)
    onthefly ?  MPA = MPEArray(AH, AL) : MPA = MPArray(AH, AL)
end
```
The MPArray uses the storage for A, so be carefull about using A for anything else.

The factorization ```mplu!``` simiply applies ```lu!``` to the low-precision matrix and the solver, which we have overloaded with ```\``` runs iterative refinement. Here are some results
```
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> using BenchmarkTools

julia> N=1024;

julia> G=Gmat(N);

julia> A = I - G;

julia> MPA=MPArray(A);

julia> @belapsed lu!(AC) setup = (AC = copy($A))
3.88038e-03

julia> @belapsed mplu!(MPC) setup = (MPC=deepcopy($MPA))
2.13324e-03
```
It is no surprise that the factorization in single precision took half as long as the one in double. In the double-single precision case, iterative refinement is a great
expample of a time/storage tradeoff. You have to store a low precision copy of $A$, so the storage burden increases by 50\% and the factoriztion time is cut in half.

Now we will see how the results look. In this example we compare the result with iterative refinement with ```A\b```, which is LAPACK's LU. 
As you can see the results are equally good. Note that the factorization object ```MPAF``` is the
output of ```mplu!```. This is analogous to ```AF=lu!(A)``` in LAPACK.

```
julia> N=1024; G=Gmat(N); A=I-G;

julia> MPA=MPArray(A); 

julia> MPF=mplu!(MPA);

julia> b=ones(N); x=MPF\b; y=A\b;

julia> norm(x-y,Inf)
4.88498e-15

julia> norm(b-A*x, Inf)
2.10942e-15

julia> norm(b-A*y,Inf)
4.44089e-15

```

As of today, you'll need to manage the factorization and the solve separately. One reason for this is that we provide several variations of iterative refinement and the solvers dispatch on the way we configure
the multiprecision array. I do not expect this to change.

### A few subtleties in the example

The constructor ```MPArray``` has two keyword arguments. The easy one to understand is ```TL``` which is the precision of the factoriztion. Julia has support for single (```Float32```) and half (```Float16```)
precisions. If you set ```TL=Float16``` then low precision will be half. Don't do that unless you know what you're doing. Using half precision is a fast way to get incorrect results. Look at the section on [half precision](#half-Precision) in this Readme for a bit more bad news.

The other keyword arguemnt is __onthefly__. That keyword controls how the triangular solvers from the factorization work. When you solve

$$ 
LU d = r
$$

The LU factors are in low precision and the residual $r$ is in high precision. If you let Julia and LAPACK figure out what to do, then the solves will be done in high precision and
the entries in the LU factors will be comverted to high precision with each binary operation. The output $d$ will be in high precision. This is called interprecision transfer on-the-fly
and ```onthefly = true``` will tell the solvers to do it that way. You have $N^2$ interprecsion transfers with each solve and, as we will see, that can have a non-trivial cost.

The default (```onthefly = false```) converts $r$ to low precision, does the solve entirely in low precision, and then promotes $d$ to high precision. You need to be careful to avoid
overflow and, more importantly, underflow when you do that and we scale $r$ to be a unit vector before conversion to low precisiion and reverse the scaling when we promote $d$. 


__MultiPrecisionArrays.jl__ supports many variations of iterative refinement and we will explain all that in the docs and in a paper in the works.

## Half Precision

Using half precision will not speed anything up, in fact will make the solver slower. The reason for this is that LAPACK and the BLAS do not (___YET__) support half precision, so all the clever stuff in 
there is missing. We proved a half precision LU factorization __/src/Factorizations/hlu!.jl__ that is better than nothing. It's a hack of Julia's  ```generic_lu!``` with threading and a couple 
complier directives. Even so, it's 2.5 -- 5 x __slower__ that a double precision LU. Half precision suppor is coming (Julia and Apple support it in hardware!) but for now, half precision is for
research in iterative refinement, not applications. Here's a table (created with  __/Code_For_Docs/HalfTime.jl__ ) that illustrates the point.

Half precision is also difficult to use properly. __Kids, don't try this at home!__. The low precsion can make iterative refinement fail because the half precision factorization can have
a large error. Here is an example to illustrate this point. The matrix here is modeslty ill-conditioned and you can see that in the error from a direct solve in double precision.

```
julia> A=I - 800.0*G;

julia> x=ones(N);

julia> b=A*x;

julia> xd=A\b;

julia> norm(b-A*xd,Inf)
6.96332e-13

julia> norm(xd-x,Inf)
2.30371e-12
```
Now, if we downcast things to half precision, nothing good happens.
```
julia> AH=Float16.(A);

julia> AHF=hlu!(AH);

julia> z=AHF\b;

julia> norm(b-A*z,Inf)
6.25650e-01

julia> norm(z-xd,Inf)
2.34975e-01
```
So you get very poor, but unsurprising, results. While __MultiPrecisionArrays.jl__ supports half precision and I use it all the time, it is not something you would use in your own
work without looking at the literature and makeing certain you are prepared for strange results. Gettting good results consistently from half precision is an active research area.


## Dependencies

As of now you need to install these packages

- Polyester
- SIAMFANLEquations

and use LinearAlgebra and SparseArrays from Base. I will remove SIAMFANLEquations from the list once I get a custom (for GMRES-IR) GMRES solver in here. Polyester is here to stay because 
threading with ```Polyester.@batch``` makes the LU for Float16 run much faster.

## Endorsement

I have used this in my own work and will be adding links to that stuff as I finish it. 

I started on this package after finishing

(KEL22a) C. T. Kelley, [__Newton's Method in Mixed Precision__](https://epubs.siam.org/doi/10.1137/20M1342902), SIAM Review 35 (1998), pp 191-211. 

(KEL22b) C. T. Kelley, [__Solving Nonlinear Equations with Iterative Methods: Solvers and Examples in Julia__](https://my.siam.org/Store/Product/viewproduct/?ProductId=44313635), SIAM, Philadelphia, 2022. 

A new paper exploits most of the package. I use the examples in that paper for CI. If you do iterative refinement well, you can make half precision work far better than it did in (KEL22a). MultiPrecisionArrays drop right into the solvers in [SIAMFANLEquations.jl](https://github.com/ctkelley/SIAMFANLEquations.jl). 

(KEL23a) C. T. Kelley, [__Newton's method in three precisions__](https://arxiv.org/abs/2307.16051) 

This paper has a [repo](https://github.com/ctkelley/N3Presults) for reproducing the results with an early version of this package.


## __Interprecision transfers__:

This is not a trivial matter. I gave a talk at the XSDK-MULTIPRECISION meeting on June 15, 2023, about this issue.

  - [Interprecision Transfers in Iterative Refinement: Making Half Precision on Desktops Less Painful](Publications_and_Presentations/MPArrays_XSDK-MULTIPRECISION_June_15.pdf).



## Funding

This project was partially supported by

1. National Science Foundation Grant DMS-1906446 and
2. Department of Energy grant DE-NA003967

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://ctkelley.github.io/MultiPrecisionArrays.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ctkelley.github.io/MultiPrecisionArrays.jl/dev

[build-status-img]: https://github.com/ctkelley/MultiPrecisionArrays.jl/workflows/CI/badge.svg
[build-status-url]: https://github.com/ctkelley/MultiPrecisionArrays.jl/actions

[codecov-img]: https://codecov.io/gh/ctkelley/MultiPrecisionArrays.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/ctkelley/MultiPrecisionArrays.jl


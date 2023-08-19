[![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url]
[![][build-status-img]][build-status-url]
[![][codecov-img]][codecov-url]
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7521427.svg)](https://doi.org/10.5281/zenodo.7521427)

# MultiPrecisionArrays.jl v0.0.4

## [C. T. Kelley](https://ctk.math.ncsu.edu)

This package provides data structures and solvers for several variants of iterative refinement. It will become much more useful when half precision (aka ```Float16```) is fully supported in LAPACK/BLAS. For now, its only general-purpose
application is classical iterative refinement with double precision equations and single precision factorizations. 

The half precision stuff is good for those of us doing research in this field. Half precision performace has progressed to the point where you can acutally get things done. On an Apple M2-Pro, a half precision LU only costs 3--5 times
what a double precision LU costs. This may be as good as it gets unless someone wants to duplicate the LAPACK implementation and get the benefits from blocking, recursion, and clever cache management.

You might take a look at [hlu!.jl](src/Factorizations/hlu!.jl), my shameless hack-job LU for half precision.

__This package is nowhere close to ready for registration or release. It's public only to help me do CI and clean up the docs.__

__Please do not make PRs. If you stumble on this mess and have questions/ideas ..., raise an issue or email me at tim_kelley@ncsu.edu__

__The half precision lu is much faster (more than 10x) as of v0.0.3. Look at [src/Factorizations/hlu.jl](src/Factorizations/hlu.jl) to see the hack job I did to [generic_lu](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/lu.jl).__

Nothing is in final form and I am changing the API, internal structures, exported functions/structs and all kinds of other stuff frequently. When/if I register this and announce it, then it will be time for complaints and offers to collaborate. 

This is the start of a package to support multiprecision arrays. This is for my own research right now.

## Readme Contents:
- [Algorithms](#What-is-iterative-refinement?)
- [Endorsement](#I'm-using-this-myself)
- [Funding](#Funding)

## What is iterative refinement?

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

Herewith, the world's most simple example to show how iterative refienment works. We will follow that with some benchmarking on the cost of factorizations.
The functions we use are __MPArray__ to create the structure and __mplu!__ to factor the low precision copy. In this example high precision is ```Float64``` and low
precision is ```Float16```.
```
julia> using MultiPrecisionArrays

julia> A=rand(5000,5000);

julia> x=ones(5000);

julia> b=A*x;

julia> MPA=MPArray(A); MPAF=mplu!(MPA); y=MPAF\b;

julia> norm(y-x,Inf)
6.65401e-11

julia> z=A\b; norm(z-x,Inf)
7.65799e-11
```
In this example we compare the result with iterative refinement with ```A\b```, which is LAPACK's LU. As you can see the results are equally good. Note that the factorization object ```MPAF``` is the
output of ```mplu!```. This is analogous to ```AF=lu!(A)``` in LAPACK.

As of today, you'll need to manage the factorization and the solve separately. One reason for this is that we provide several variations of iterative refinement and the solvers dispatch on the way we configure
the multiprecision array. I do not expect this to change.

As for the cost.
```
julia> using MultiPrecisionArrays

julia> using BenchmarkTools

julia> A=rand(5000,5000);

julia> @belapsed lu!(AA) setup = (AA = copy($A))
2.85270e-01

julia> @belapsed mplu!(MPAA) setup = (MPAA = deepcopy($MPA))
1.32841e-01
```
It is no surprise that the factorization in single precision took half as long as the one in double. In the double-single precision case, iterative refinement is a great
expample of a time/storage tradeoff. You have to store a low precision copy of $A$, so the storage burden increases by 50\% and the factoriztion time is cut in half.

__MultiPrecisionArrays.jl__ supports many variations of iterative refinement and we will explain all that in the docs and in a paper in the works.


## Dependencies

As of now you need to install these packages

- Polyester
- SIAMFANLEquations

and use LinearAlgebra and SparseArrays from Base. I will remove SIAMFANLEquations from the list once I get a custom (for GMRES-IR) GMRES solver in here. Polyester is here to stay because 
threading with ```Polyester.@batch``` makes the LU for Float16 run much faster.


## I'm using this myself.

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


[![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url]
[![][build-status-img]][build-status-url]
[![][codecov-img]][codecov-url]
[![DOI](https://zenodo.org/badge/462014459.svg)](https://doi.org/10.5281/zenodo.13851500)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06698/status.svg)](https://doi.org/10.21105/joss.06698)


<!--- [![MultiPrecisionArrays Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/MultiPrecisionArrays)](https://pkgs.genieframework.com?packages=MultiPrecisionArrays) --->


# MultiPrecisionArrays.jl v0.1.4
  
## [C. T. Kelley](https://ctk.math.ncsu.edu)

This package provides data structures and solvers for several variants of iterative refinement (IR). It will become much more useful when half precision (aka ```Float16```) is fully supported in LAPACK/BLAS. 

For now, its only general-purpose
applications are classical iterative refinement with double precision equations and single precision factorizations or, for you Wilkinsonian ill-conditioned problem connisseurs, single precision equations and factorizations with double precision residuals.

The half precision stuff is good for those of us doing research in this field. Half precision performance has progressed to the point where you can actually get things done. On an Apple M2-Pro, a half precision LU only costs 3--8 times
what a double precision LU costs. This may be as good as it gets unless someone wants to duplicate the LAPACK implementation and get the benefits from blocking, recursion, and clever cache management.

__The half precision LU for Float16 in this package is much faster (more than 10x on my M2 Pro) than the default in Julia. Look at [src/Factorizations/hlu!.jl](src/Factorizations/hlu!.jl) to see the hack job I did to [generic_lu](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/lu.jl).__

## What parts of MultiPrecisionArrays.jl are ready for prime time?

- Using [```mplu```](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/functions/mplu/) and [```mplu!```](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/functions/mplu!/) to cut your factorization time in half for double precision matrices is simple and works well.
- The API for [harvesting iteration statistics](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Details/Stats/) is stable.
- [Using extended precision for IR](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Details/Extended) works the way Wilksinson said it did. It's fun (but slow) to play with ```TR=Float64, TW=Float32, TF=Float16```. The old-time religion works with "TW=TF=Float32" and "TR= Float64". You might try "TW=TF=Float64" and "TR=Double64" to see what happens. 
- If you're a half precision person,
   - GMRES-IR works with [```mpglu```](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/functions/mpglu/) and [```mpglu!```](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/functions/mpglu!/)
   - BiCGSTAB-IR works with  [```mpblu```](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/functions/mpblu/) and [```mpblu!```](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/functions/mpblu!/)
   - the new half precision LU factorization [```hlu!```](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/functions/hlu!/) is __only__ 3-5x slower than a double precision ```lu```.

## What's new?

- v0.1.3: Still better docs and ..
   - Fixing a performance bug.
   - Add options to termination criterion. __Change default back to small residuals.__
   - v0.1.3 is the version for the [__JOSS paper__](https://joss.theoj.org/papers/10.21105/joss.06698). I gave that version the v0.1.3-joss tag.
 
 - v0.1.4: Continuous improvement for the docs and ...
    - Enable fine control of termination criteria parameters 

 - v0.1.5: QR and linear least squares
      
##  Can I complain about this package?

Yes, but ...

__Please do not make PRs. If you have questions/ideas ..., raise an issue or email me at tim_kelley@ncsu.edu__

Since log(version number) < 0, you can expect changes in API, internal structures, and exported functions.

## Coming Attractions

- Nonlinear solver applications
- More factorizations: cholesky, qr
- BFloat16 when/if the support is there. You can use it now, but it is pathologically slow.
  - You'll need [BFloat16s.jl](https://github.com/JuliaMath/BFloat16s.jl) because Julia does not support BFloat16 natively.
  - There's work in progress. See issue [41075](https://github.com/JuliaLang/julia/issues/41075) and PRs [51470](https://github.com/JuliaLang/julia/pull/51470) and [53059](https://github.com/JuliaLang/julia/pull/53059).

## Documentation

- This README file
- The [docs for the package](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev)
- [The ArXiv paper for v0.1.3](https://arxiv.org/pdf/2311.14616.pdf)
- [The working paper for v0.1.4](https://github.com/ctkelley/MultiPrecisionArrays.jl/blob/main/Publications_and_Presentations/MPArray.pdf) 
- [The JOSS paper](https://joss.theoj.org/papers/10.21105/joss.06698)
## Readme Contents:

- [Algorithms](#algorithms)
- [Example](#example)
  - [A few details](#a-few-details)
- [Be Careful with Half Precision](#half-precision)
- [Harvesting Iteration Statistics](#harvesting-iteration-statistics)
- [Residuals in Higher Precision](#evaluating-residuals-in-higher-precision-than-the-working-precision)
- [Dependencies](#dependencies)
- [Endorsement](#endorsement)
- [Other Stuff in Julia](#other-ir-in-julia)
- [Funding](#funding)


## Algorithms

### What is iterative refinement?

This package will make solving dense systems of linear equations faster or more accurate by using the LU factorization and IR. It is limited to LU for now. A very generic description of this for solving a linear system $A x = b$ faster is

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
to allocate storage for a copy of $A$ in the lower precision and factor that copy. Then one has to determine what the line
$d = (LU)^{-1} r$ means. Do you cast $r$ into the lower precision before the solve or not? __MultiPrecisionArrays.jl__ provides
data structures and solvers to manage this. The __MPArray__ structure lets you preallocate $A$, the low precision copy, and the residual $r$.
The factorizations factor the low-precision copy and the solvers use that factorization and the original high-precision matrix to run
the while loop in the algorithm. We encode all this is functions like ```mplu```, which builds a factorization object that does IR 
when you use the backslash operator ```\```.

IR is a perfect example of a storage/time tradeoff.
To solve a linear system $A x = b$ in $R^N$ with IR,
one incurs the storage penalty of making a low
precision copy of $A$ and reaps the benefit of only having to
factor the low precision copy.

While you might think this is a good idea for all problems, it is not so good for smaller problems. The reason is that IR swaps the factoriation cost for
a matrix-vector multiply and the two triangular solves for LU __in each IR iteration__. Triangular solves do not thread as well as factorizations or matrix-vector
multiplies and that can affect the performance in a significant way, even though it is only $N^2$ work. The details are [in the docs](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Details/N2Work).



## Example

### An example to get started

Herewith, the world's most simple example to show how iterative refinement works. We will follow that with some benchmarking on the cost of factorizations.

The results in this README use Julia 1.10.4 and openBLAS. You should consider trying MKL (Intel machines) or [AppleAccelerate](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl). I have switched to AppleAccelerate myself. It gets faster with every new OS release, but the triangular solvers are still slower than openBLAS. The performance (iteration counts and timings) of IR depend on versions of OS/Julia/MultiPrecsionArrays and hour choice of BLAS. I have used several different versions in the docs and the ArXiV paper and you should try the
examples in your own environment.

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

The example below compares the cost of a double precision factorization to a MPArray factorization. The ```MPArray``` structure has a working
precision (```TW```) and a factorization (low)  precision (```TF```) matrix. The structure we will start with 
is
```
struct MPArray{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat}
    AH::Array{TW,2}
    AL::Array{TF,2}
    residual::Vector{TR}
    sol::Vector{TR}
    onthefly::Bool
end
```
The structure also stores the residual and the solution in the residual precision ```TR```. Most of the time ```TR = TW```.
The ```onthefly``` Boolean tells the solver how to do the interprecision transfers. The easy way to get started is to use the ```mplu``` 
command directly on the matrix. That will build the MPArray, follow that with the factorization of ```AL```, and put in all in a structure
that you can use as a factorization object with ```\```. I do not export the constructor for this or any other multiprecision array structure. You should use functions like ```mplu``` to build 
multiprecision factorization objects.

Now we will see how the results look. In this example we compare the result with iterative refinement with ```A\b```, which is LAPACK's LU. 
As you can see the results are equally good. Note that the factorization object ```MPF``` is the
output of ```mplu```. This is analogous to ```AF=lu(A)``` in LAPACK.

You may not get exactly the same results for this example on
different hardware, BLAS, number of cores, versions of Julia/OS/MulitPrecisionArrays.
I am still playing with the termination criteria and the iteration
count could grow or shrink as I do that as could the residual for the converged result.

```
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> using BenchmarkTools
julia> N=4096;

julia> G=Gmat(N);

julia> A = I - G;

julia> x=ones(N); b=A*x;

julia> MPF=mplu(A); AF=lu(A);

julia> z=MPF\b; w=AF\b;

julia> ze=norm(z-x,Inf); zr=norm(b-A*z,Inf)/norm(b,Inf);

julia> we=norm(w-x,Inf); wr=norm(b-A*w,Inf)/norm(b,Inf);

julia> println("Errors: $ze, $we. Residuals: $zr, $wr")
Errors: 5.55112e-16, 6.68354e-14. Residuals: 6.66134e-16, 6.68354e-14
```

So the results are equally good.

The compute time for ```mplu``` should be a bit more than half that of ```lu!```. The reason is
that ```mplu``` factors a low precision array, so the factorization cost is cut in half. Memory
is a different story because. The reason
is that both ```mplu``` and ```lu!``` do not allocate storage for a new high precision array,

but ```mplu``` allocates for a low precision copy, so the memory and allocation cost for ```mplu```
is 50% more than ```lu```. 

One issue with smaller problems is that the triangular solve
does not parallelize as well as the factorization, so does not exploit
multi-core processor as well. We can see this in the IR solver times because
each iteration of IR needs a matrix-vector multiply and a triangular solve.

```
julia> @belapsed mplu($A)
8.59528e-02

julia> @belapsed lu!(AC) setup=(AC=copy($A))
1.42112e-01

# And now for the solve times.

julia> @belapsed ldiv!($AF,bb) setup=(bb = copy($b))
4.79117e-03

julia> @belapsed $MPF\$b
2.01195e-02
```
So the total solve time is less, but the $O(N^2)$ work is not zero.

It is no surprise that the factorization in single precision took roughly half as long as the one in double. In the double-single precision case, iterative refinement is a great
example of a time/storage tradeoff. You have to store a low precision copy of $A$, so the storage burden increases by 50\% and the factorization time is cut in half.
The advantages of IR increase as the dimension increases. IR is less impressive for smaller problems and can even be slower for the factorization.
```
julia> N=30; A=I + Gmat(N); 

julia> @belapsed mplu($A)
5.22217e-06

julia> @belapsed lu!(AC) setup=(AC=copy($A))
3.64062e-06
```

### A few details

Look at the docs for things like
- [Memory allocation costs](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/#Memory-Allocations-for-mplu)
- [Terminating the while loop](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Details/Termination/)
- [Is O(N^2) work really negligible?](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Details/N2Work)
- [Options and data structures for mplu](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/#Options-and-data-structures-for-mplu)



__MultiPrecisionArrays.jl__ supports many variations of iterative refinement and we explain all that in the docs and a [users guide](https://github.com/ctkelley/MultiPrecisionArrays.jl/blob/main/Publications_and_Presentations/MPArray.pdf).

## Half Precision

Bottom line: don't use it and expect good performance. See [this page](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Half_1/) in the docs for details.

It's a good idea to try GMRES-IR if you are playing with half precision. GMRES-IR uses a different factorization ```mpglu``` which factors the 
low precision matrix, allocates room for the Krylov basis, and builds a factorization object. GMRES-IR uses the low precision factorization as
a preconditioner for GMRES. The call looks like ```MPGF=mpglu(A)```. You solve $A x = b$ with
```x = MPGF\b```. More details [here](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Half_1/#GMRES-IR).

## Evaluating residuals in higher precision than the working precision.

This is the classical application. __MultiPrecisionArrays.jl__ supports that. It's all explained [here](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Details/Extended/) in the docs. Doing this can solve some ill-conditioning problems, 
but it is not something you should do routinely.

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
2.22045e-16

julia> # Now look at the residual history

julia> mpout.rhist
6-element Vector{Float64}:
 1.00000e+00
 1.89553e-05
 4.56056e-11
 3.08642e-14
 4.44089e-16
 2.22045e-16
```
As you can see, IR does well for this problem. The package uses an initial
iterate of $x = 0$ and so the initial residual is simply $r = b$
and the first entry in the residual history is $|| b ||_\infty$. The
iteration terminates successfully after five matrix-vector products.

You may wonder why the residual after the first iteration was so much
larger than single precision roundoff. The reason is that the default 
when the low precision is single is to scale and downcast the residual
for the triangular solves. This is faster for medium sized problems.

One can enable interprecision transfers on the fly and see the difference.

 ```
julia> MPF2=mplu(A; onthefly=true);

julia> @belapsed $MPF2\$b
2.47693e-02

julia> mpout2=\(MPF2, b; reporting=true);

julia> mpout2.rhist
5-element Vector{Float64}:
 5-element Vector{Float64}:
 1.00000e+00
 6.29385e-07
 3.98570e-13
 5.21805e-15
 2.22045e-16
```
So the second iteration took fewer iterations but a bit more time.



There are more examples for this in the [paper](https://github.com/ctkelley/MultiPrecisionArrays.jl/blob/main/Publications_and_Presentations/MPArray.pdf).



## Dependencies

As of now these packages are the only direct dependencies

- OhMyThreads
- Reexport
- SIAMFANLEquations

as we use LinearAlgebra and SparseArrays from Base. I use the Krylov solvers and examples from SIAMFANLEquations. OhMyThreads is here because 
threading with ```OhMyThreads.@tasks``` makes the LU for Float16 run much faster. Once LU for Float16 is in LAPACK/BLAS, I will eliminate that dependency.

I use Reexport to make the identity ```I``` part of the Examples submodule. This makes it easier to run the examples in the docs.

## Endorsement

I have used this in my own work and will be adding links to that stuff as I finish it. 

I started on this package after finishing

- (KEL22a) C. T. Kelley, [__Newton's Method in Mixed Precision__](https://epubs.siam.org/doi/10.1137/20M1342902), SIAM Review 35 (1998), pp 191-211. 

- (KEL22b) C. T. Kelley, [__Solving Nonlinear Equations with Iterative Methods: Solvers and Examples in Julia__](https://my.siam.org/Store/Product/viewproduct/?ProductId=44313635), SIAM, Philadelphia, 2022. 

A new paper exploits most of the package. I use the examples in that paper for CI. If you do iterative refinement well, you can make half precision work far better than it did in (KEL22a). MultiPrecisionArrays drop right into the solvers in [SIAMFANLEquations.jl](https://github.com/ctkelley/SIAMFANLEquations.jl). 

- (KEL23a) C. T. Kelley, [__Newton's method in three precisions__](https://arxiv.org/abs/2307.16051) To appear in Pacific Journal of Optimization.
  - This paper has a [repo](https://github.com/ctkelley/N3Presults) for reproducing the results with an early version of this package.

I just submitted a paper about interprecison transfers and used MulitPrecisionArrays.jl for the computations.

- (Kel24a) C. T. Kelley, [__Interprecision transfers in iterative refinement__](https://arxiv.org/pdf/2407.00827)

Here's a new JOSS paper about the package.

- (Kel24b) C. T. Kelley, [__MultiPrecisionArrays.jl: A Julia package for iterative refinement__](https://joss.theoj.org/papers/10.21105/joss.06698)

## __Interprecision transfers__:

This is not a trivial matter. I gave a talk at the XSDK-MULTIPRECISION meeting on June 15, 2023, about this issue. See [the docs](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/Details/Interprecision_1/) for the story on this. 
Most users of this package can ignore this issue.

  - [Interprecision Transfers in Iterative Refinement: Making Half Precision on Desktops Less Painful](Publications_and_Presentations/MPArrays_XSDK-MULTIPRECISION_June_15.pdf).

Also see ref (Kel 24a) [__Interprecision transfers in iterative refinement__](https://arxiv.org/pdf/2407.00827)

## Other IR in Julia

See [the docs](https://ctkelley.github.io/MultiPrecisionArrays.jl/dev/#Other-IR-software-in-Julia) for a couple citations. 

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


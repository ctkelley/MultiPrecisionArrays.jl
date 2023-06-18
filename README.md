[![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url]
[![][build-status-img]][build-status-url]
[![][codecov-img]][codecov-url]
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7521427.svg)](https://doi.org/10.5281/zenodo.7521427)

# MultiPrecisionArrays.jl v0.0.3

## [C. T. Kelley](https://ctk.math.ncsu.edu)

__This package is nowhere close to ready for registration or release. It's public only to help me do CI and clean up the docs.__

__Please do not make PRs. If you stumble on this mess and have questions/ideas ..., raise an issue or email me at tim_kelley@ncsu.edu__

Nothing is in final form and I am changing the API, internal structures, exported functions/structs and all kinds of other stuff frequently. When/if I register this and announce it, then it will be time for complaints and offers to collaborate. 

This is the start of a package to support multiprecision arrays. This is for my own research right now.

## Readme Contents:
- [Endorsement](##I'm-using-this-myself)
- [Algorithms](##What-is-iterative-refinement?)
- [Funding](##Funding)

## I'm using this myself.

I have used this in my own work and will be adding links to that stuff as I finish it. 

I started on this package after finishing

(KEL22a) C. T. Kelley, [__Newton's Method in Mixed Precision__](https://epubs.siam.org/doi/10.1137/20M1342902), SIAM Review 35 (1998), pp 191-211. 

(KEL22b) C. T. Kelley, [__Solving Nonlinear Equations with Iterative Methods: Solvers and Examples in Julia__](https://my.siam.org/Store/Product/viewproduct/?ProductId=44313635), SIAM, Philadelphia, 2022. 

- __Interprecision transfers__:This is not a trivial matter. I gave [this talk](Publications_and_Presentations/MPArrays_XSDK-MULTIPRECISION_June_15.pdf) at the XSDK-MULTIPRECSSION meeting on June 15, 2023, about this issue.
- __Newton's method in three precisions:__ If you do iterative refinement well, you can make half precision work far better than it did in (KEL22a). 

## What is iterative refinement?

One example is a struct with a Float64 array A and Float32.(A). To do factor/solve you factor the Float32 array and solve with iterative refinement.

I have this working to the point where ```\``` does the right thing.

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


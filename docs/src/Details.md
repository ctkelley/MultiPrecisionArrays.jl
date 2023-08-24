# More than you want to know

## Interprecision Transfers: Part I


## Half Precision

Using half precision will not speed anything up, in fact will make the solver slower. The reason for this is that LAPACK and the BLAS do not (__YET__) support half precision, so all the clever stuff in
there is missing. We provide a half precision LU factorization __/src/Factorizations/hlu!.jl__ that is better than nothing. It's a hack of Julia's  ```generic_lu!``` with threading and a couple
complier directives. Even so, it's 2.5 -- 5 x __slower__ that a double precision LU. Half precision suppor is coming (Julia and Apple support it in hardware!) but for now, half precision is for
research in iterative refinement, not applications. Here's a table (created with  __/Code_For_Docs/HalfTime.jl__ ) that illustrates the point. In the table we compare LAPACK's LU to the LU we compute with ```hlu!.jl```.

Half precision is also difficult to use properly. __Kids, don't try this at home!__. The low precsion can make iterative refinement fail because the half precision factorization can have
a large error. Here is an example to illustrate this point. The matrix here is modestly ill-conditioned and you can see that in the error from a direct solve in double precision.

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


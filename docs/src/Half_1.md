# Half Precision and Krylov-IR

There are two half precision (16 bit) formats. Julia has native
support for IEEE 16 bit floats (Float16). A second format
(BFloat16) has a larger exponent field and a smaller 
significand (mantissa), thereby
trading precision for range. In fact, the exponent field in BFloat is
the same size (8 bits) as that for single precision (Float32). The 
significand, however, is only 8 bits. Which is less than that for 
Float16 (11 bits) and single (24 bits). The size of the significand
means that you can get in real trouble with half precision in either
format.  

At this point we offer no support for BFloat16. Progress is being
made.

Using half precision will not speed anything up, in fact it will make 
the solver slower. The reason for this is that LAPACK and the BLAS 
do not (__YET__) support half precision, so all the clever stuff in
there is missing. We provide a half precision LU
factorization for Float16
 __/src/Factorizations/hlu!.jl__ that is better than nothing. 
It's a hack of Julia's  ```generic_lu!``` with threading and a couple
compiler directives. Even so, it's 2.5 -- 5 x __slower__ than a 
double precision LU. Half precision support for both Float16
and BFloat16 is coming 
(Julia and Apple support it in hardware!) but for now, at least for desktop
computing, half precision is for
research in iterative refinement, not applications. 


Here's a table (created with  __/Code_For_Docs/HalfTime.jl__ ) that illustrates the point. In the table we compare timings for
LAPACK's LU to the LU we compute with ```hlu!.jl```. The matrix is 
$I-800.0*G$.

```
      N       F64       F32       F16     F16/F64 
     1024  3.65e-03  2.65e-03  5.26e-03  1.44e+00 
     2048  2.26e-02  1.41e-02  3.70e-02  1.64e+00 
     4096  1.55e-01  8.53e-02  2.55e-01  1.65e+00 
     8192  1.15e+00  6.05e-01  4.23e+00  3.69e+00 
```
The columns of the table are the dimension of the problem, timings
for double, single, and half precision, and the ratio of the half
precision timings to double. The timings came from Julia 1.10-beta2
running on an Apple M2 Pro with 8 performance cores.

## Half Precision is Subtle

Half precision is also difficult to use properly. The low precision can 
make iterative refinement fail because the half precision factorization 
can have a large error. Here is an example to illustrate this point. 
The matrix here is modestly ill-conditioned and you can see that in the 
error from a direct solve in double precision.

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
work without looking at the literature and making certain you are prepared for strange results. Getting good results consistently from half precision is an active research area.

So, it should not be a surprise that IR also struggles with half precision.
We will illustrate this with one simple example. In this example high
precision will be single and low will be half. Using {\bf MPArray} with
a single precision matrix will automatically make the low precision matrix
half precision.
```
julia> N=4096; G=800.0*Gmat(N); A=I - Float32.(G);

julia> x=ones(Float32,N); b=A*x;

julia> MPF=mplu(A; onthefly=false);

julia> y=MPF\b;

julia> norm(b - A*y,Inf)
1.05272e+02
```
So, IR completely failed for this example. We will show how to extract
the details of the iteration in a later section.

It is also worthwhile to see if doing the triangular solves on-the-fly
(MPS) helps. 

```
julia> MPBF=mplu(A);

julia> z=MPBF\b;

julia> norm(b-A*z,Inf)
1.28174e-03
```
So, MPS is better in the half precision case. Moreover, it is also less
costly thanks to the limited support for half precision computing.
For that reason, MPS is the default when high precision is single.

However, on-the-fly solves are not enough to get good results and IR
still terminates too soon.

## GMRES-IR

GMRES-IR solves the correction equation
with a preconditioned GMRES [gmres](@cite) iteration. One way to think of this
is that the solve in the IR loop is an approximate solver for the
correction equation
```math
A d = r
```
where one replaces $A$ with the low precision factors
$LU$. In GMRES-IR one solves the correction
equation with a left-preconditioned GMRES iteration using
$U^{-1} L^{-1}$ as
the preconditioner. The preconditioned equation is
```math
U^{-1} L^{-1}  A d = U^{-1} L^{-1} r.
```

GMRES-IR will not be as efficient as IR because each iteration is itself
an GMRES iteration and application of the preconditioned matrix-vector
product has the same cost (solve + high precision matrix vector product)
as a single IR iteration. However, if low precision is half, this approach
can recover the residual norm one would get from a successful IR iteration.

There is also a storage problem. One should allocate storage for the Krylov
basis vectors and other vectors that GMRES needs internally. We do that
in the factorization phase. So the structure __MPGEFact__ has the 
factorization of the low precision matrix, the residual, the Krylov
basis and some other vectors needed in the solve. 

The Julia function
```mpglu``` constructs the data structure and factors the low precision
copy of the matrix. The output, like that of ```mplu``` is a factorization
object that you can use with backslash.

Here is a well conditioned example. Both IR and GMRES-IR perform well, with
GMRES-IR taking significantly more time.  

```
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> using BenchmarkTools

julia> N=4069; AD= I - Gmat(N); A=Float32.(AD); x=ones(Float32,N); b=A*x;

julia> MPF=mplu(A); MPF2=mpglu(A);

julia> z=MPF\b; y=MPF2\b; println(norm(z-x,Inf),"  ",norm(y-x,Inf))
5.9604645e-7  4.7683716e-7

julia> @btime $MPF\$b;
  13.582 ms (4 allocations: 24.33 KiB)

julia> @btime $MPF2\$b;
  40.344 ms (183 allocations: 90.55 KiB)
```

If you dig into the iterations statistics (more on that later) you
will see that the GMRES-IR iteration took almost exactly four times
as many solves and residual computations as the simple IR solve.

We will repeat this experiment on the ill-conditioned example. In this
example, as we saw earlier, IR fails to converge.

```
julia> N=4069; AD= I - 800.0*Gmat(N); A=Float32.(AD); x=ones(Float32,N); b=A*x;
        
julia> MPF=mplu(A); MPF2=mpglu(A);

julia> z=MPF\b; y=MPF2\b; println(norm(z-x,Inf),"  ",norm(y-x,Inf))
0.2875508  0.004160166

julia> println(norm(b-A*z,Inf)/norm(b,Inf),"  ",norm(b-A*y,Inf)/norm(b,Inf))
0.0012593127  1.4025759e-5
```
So, the relative error and relative residual norm for GMRES-IR
is much smaller than that for IR.


## BiCGSTAB-IR

__MultiPrecisionArrays.jl__ also supports BiCGSTAB [bicgstab](@cite). The
API is the same as for GMRES-IR with the swap of ```b``` for ```g```. So the
data structure is __MPBArray__ and the factorizations are
```mpblu``` and ```mpblu!```. Keep in mind that a BiCGSTAB iteration costs
two matrix-vector and two preconditioner-vector products.

We will do the same examples we did for GMRES-IR.
```
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> using BenchmarkTools

julia> N=4069; AD= I - Gmat(N); A=Float32.(AD); x=ones(Float32,N); b=A*x;

julia> MPFB=mpblu(A); MPF=mplu(A);

julia> z=MPF\b; y=MPFB\b; println(norm(z-x,Inf),"  ",norm(y-x,Inf))
5.9604645e-7  4.172325e-7

julia> @btime $MPF\$b;
  13.371 ms (5 allocations: 24.38 KiB)

julia> @btime $MPFB\$b;
  77.605 ms (178 allocations: 821.45 KiB)
```
The only new thing is that the solve for BiCGSAB-IR took roughly twice
as long. That is no surprise since the cost of a BiCGSTAB iteration is
about double that of a GMRES iteration.

The ill-conditioned example tells the same story.
```
julia> N=4069; AD= I - 800.0*Gmat(N); A=Float32.(AD); x=ones(Float32,N); b=A*x;

julia> MPF=mplu(A); MPFB2=mpblu(A);

julia> z=MPF\b; y=MPFB2\b; println(norm(z-x,Inf),"  ",norm(y-x,Inf))
0.2875508  0.0026364326

julia> println(norm(b-A*z,Inf)/norm(b,Inf),"  ",norm(b-A*y,Inf)/norm(b,Inf))
0.0012593127  7.86059e-6
```
So, as was the case with GMRES, BiCGSTAB-IR is a much better solver
that IR alone.


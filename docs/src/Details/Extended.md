# Evaluating residuals in higher precision

This idea comes from [Wilkinson48](@cite) and I am using the notation
from [demmelir](@cite) and [amestoy:2024](@cite).
The idea is to evaluate the
residual in a precision ```TR``` 
higher than the working precision. If you do this,
you should store both the solution and the residual in precision
TR and to interprecision transfers on the fly. In that case you are really
solving a promoted problem [ctk:irnote](@cite)
```math
(I_W^R A) x = I_W^R b
```
and, by driving the residual to a small value can mitigate ill-conditioning
to some degree. Here $I_P^Q$ is the interprecision transfer from precision
$P$ to precision $Q$.
__MultiPrecisionArrays__  
allows you to do this with the multiprecision
factorization you get from ```mplu```. 

The classic example is to let ```TR``` 
and ```TF``` be single precision and ```TR``` be double.
The storage penalty is that you must store two copies of $A$, one for the
residual computation and the other for the factorization.

Here is an example with a badly conditioned matrix. You must tell
```mplu``` to factor in the working precision and use the
```kwargs``` in ```mplu``` to set ```TR```.

You may not get exactly the same results for this example on
different hardware, BLAS, versions of Julia or this package.
I am still playing with the termination criteria and the iteration
count could grow or shrink as I do that.


```
julia> using MultiPrecisionArrays

julia> using MultiPrecisionArrays.Examples

julia> N=4096; alpha=799.0; AD=I - alpha*Gmat(N);

# conditioning is bad

julia> cond(AD,Inf)
2.35899e+05

# Set up the single precision computation

julia> A = Float32.(AD); xe=ones(Float32,N); b=A*xe;

# Make sure TF is what it needs to be for this example

julia> AF = mplu(A; TF=Float32, TR=Float64);

# Use the multiprecision array to solve the problem, set TR.
julia> mrout = \(AF, b; reporting=true);

# look at the residual history

julia> mrout.rhist
6-element Vector{Float64}:
 9.88750e+01
 3.94930e-03
 7.66242e-07
 1.13667e-07
 7.53089e-08
 7.53089e-08

# Compare results with LU and the exact(?) solution

julia> xr=Float32.(mrout.sol); xs = A\b;

julia> [norm(b - A*xr, Inf) norm(b - A*xs, Inf)]
1×2 Matrix{Float32}:
 1.69373e-03  2.53296e-03

# So the residual is bit better. What about the difference from xe?

julia> [norm(xr - xe, Inf) norm(xs - xe, Inf)]
1×2 Matrix{Float32}:
 2.27022e-03  7.20143e-03

# Nothing very exciting here. You have to wonder what this all means.
# Finally, how did we do with the promoted problem?

julia> AP=Float64.(A);

julia> norm(b - AP*mrout.sol, Inf)
7.53089e-08

# Which is consistent with theory.
```

So, is the solution to the promoted problem better than the exact solution
I used to build the problem?

The reader might try this with ```TF=Float16```, the default when
```TW = Float32```. All that you'll need to do is replace
```
AF = mplu(A; TF=Float32, TR=Float64);
```
with
```
AF = mplu(A; TR=Float64);
```

What goes wrong and why? Fixup to follow.

The advantages of evaluating the residual in extended precision grow
when $A$ is extremely ill-conditioned. Of course, in this case the
factorization in the factorization precision could be so inaccurate that
IR will fail to converge. One approach to respond to this, as you
might expect, is to use the factorization as a preconditioner and not
a solver [amestoy:2024](@cite). We will support this in a later version of
__MultiPrecisionArrays.jl__.

## IR-Krylov with high precision residuals

Half precision factorizations can lead to failure in IR. IR-Krylov
methods try to fix this by using the low precision factorization as
a preconditioner, not a solver. __MultiPrecisionArrays.jl__ has two
IR-Krylov methods and one can adjust the residual precision just as one
does with ```mplu```. Here is an example of this using the two 
multiprecision IR-Krylov factorizations ```mpglu``` (GMRES) and 
```mpblu``` (BiCGSTAB).

Using the same example, we examine the results.

```
julia> AF = mplu(A; TR=Float64);

julia> mout16=\(AF, b; reporting=true);

julia> mout16.rhist
4-element Vector{Float64}:
 9.88750e+01
 3.92451e+00
 3.34292e-01
 2.02204e-01
```
so plain vanilla IR with ```TF=Float16```, ```TW=Float32```, and
```TR=Float64``` fails to converge. 

We support both IR-GMRES and IR-BiCGSTAB for ```TR > TW```. You get this to
work just like in {\tt mplu} by using the keyword argument {\tt TR}.
We will continue with the example in this section and do that. For this
example the default basis size of $10$ is not enough, so we use 20.

```
julia> GF = mpglu(A; TR=Float64, basissize=20);

julia> moutG = \(GF, b; reporting=true);

julia> moutG.rhist
8-element Vector{Float64}:
 9.88750e+01
 7.59075e-03
 1.48842e-05
 2.17281e-07
 8.60429e-08
 7.45077e-08
 7.91866e-08
 7.53089e-08

julia> moutG.dhist
7-element Vector{Float32}:
 1.03306e+00
 3.23734e-02
 2.83073e-04
 8.92629e-06
 1.55432e-07
 6.05685e-08
 5.96633e-08

julia> xp = Float64.(A)\b;

julia> norm(xp-moutG.sol, Inf)
5.95825e-08
```

IR-BiCGSTAB whould take fewer iterations than IR-GMRES had we used
the default basis size because there's no storage
issue. But remember that BiCGSTAB has a higher cost per linear iteration.

```
julia> BF = mpblu(A; TR=Float64);

julia> moutB = \(BF, b; reporting=true);

julia> moutB.rhist
5-element Vector{Float64}:
 9.88750e+01
 1.86437e-07
 7.53089e-08
 7.53089e-08
 7.53089e-08

julia> moutB.dhist
4-element Vector{Float32}:
 1.00227e+00
 2.05457e-06
 5.95821e-08
 5.95821e-08

julia> norm(xp - moutB.sol, Inf)
5.95825e-08
```

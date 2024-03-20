# Evaluating residuals in higher precision

This idea comes from [Wilkinson48](@cite) and I am using the notation
from [demmelir](@cite) and [amestoy:2023](@cite).
The idea is to evaluate the
residual in a precision ```TR``` 
higher than the working precision. If you do this,
you should store both the solution and the residual in precision
TR and to interprecision transfers on the fly. In that case you are really
solving a promoted problem
```math
(I_W^R A) x = I_W^R b
```
and, by driving the residual to a small value can mitigate ill-conditioning
to some degree. Here $I_P^Q$ is the interprecision transfer from precision
$P$ to precision $Q$.
__MultiPrecisionArrays__  
allows you to do this with the multiprecision
factorization you get from ```mplu```. I may add this feature to the
Krylov-IR solvers later.

The classic example is to let ```TR``` 
and ```TF``` be single precision and ```TR``` be double.
The storage penalty is that you must store two copies of $A$, one for the
residual computation and the other for the factorization.

Here is an example with a badly conditioned matrix. You must tell
```mplu``` to factor in the working precision and then use the
```kwargs``` in the solver to set ```TR```.

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

julia> AF = mplu(A; TF=Float32);

# Use the multiprecision array to solve the problem, set TR.

julia> mrout = \(AF, b; reporting=true, TR=Float64);

# look at the residual history

julia> mrout.rhist
5-element Vector{Float64}:
 9.88750e+01
 1.67735e-05
 9.23976e-10
 9.37916e-13
 9.09495e-13

# Compare results with LU and the exact(?) solution

julia> xr=Float32.(mrout.sol); xs = A\b;

julia> [norm(b - A*xr, Inf) norm(b - A*xs, Inf)]
1×2 Matrix{Float32}:
 1.29700e-04  1.41144e-03

# So the residual is better. What about the difference from xe?

julia> [norm(xr - xe, Inf) norm(xs - xe, Inf)]
1×2 Matrix{Float32}:
 8.47816e-04  8.82089e-04

# Nothing exciting here. You have to wonder what this all means.
# Finally, how did we do with the promoted problem?

julia> AP=Float64.(A);

julia> norm(b - AP*mrout.sol, Inf)
7.10543e-13

# Which is what I said it was above.
```

So, is the solution to the promoted problem better than the exact solution
I used to build the problem?

The reader might try this with ```TF=Float16```, the default when
```TW = Float32```. All that you'll need to do is replace
```
AF = mplu(A; TF=Float32);
```
with
```
AF = mplu(A);
```

What goes wrong and why?

The advantages of evaluating the residual in extended precision grow
when $A$ is extremely ill-conditioned. Of course, in this case the
factorization in the working precision could be so inaccurate that
IR will fail to converge. One approach to respond to this, as you
might expect, is to use the factorization as a preconditioner and not
a solver [amestoy:2023](@cite). We will support this in a later version of
__MultiPrecisionArrays.jl__.


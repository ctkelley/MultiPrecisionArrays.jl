# Harvesting Iteration Statistics

You can get some iteration statistics by using the
__reporting__ keyword argument to the solvers. The easiest way
to do this is with the backslash command. When you use this option you
get a data structure with the solution and the residual history.

You may not get exactly the same results for this example on
different hardware, BLAS, versions of Julia or this package.
I am still playing with the termination criteria and the iteration
count could grow or shrink as I do that.

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
 1.64859e-05
 1.02649e-10
 6.35048e-14
 4.44089e-16
 2.22045e-16
```
As you can see, IR does well for this problem. The package uses an initial
iterate of $x = 0$ and so the initial residual is simply $r = b$
and the first entry in the residual history is $|| b ||_\infty$. The
iteration terminates successfully after five matrix-vector products.

You can also look at the norm of the defect.
```
julia> mpout.dhist
5-element Vector{Float64}:
 1.00002e+00
 1.62124e-05
 1.02649e-10
 6.35100e-14
 4.45650e-16
```

You may wonder why the residual after the first iteration was so much
larger than single precision roundoff. The reason is that the default 
when the low precision is single is to downcast the residual to single before 
the solve (onthefly=false).
One can enable interprecision transfers on the fly and see the difference.

```
julia> MPF2=mplu(A; onthefly=true);

julia> mpout2=\(MPF2, b; reporting=true);

julia> mpout2.rhist
6-element Vector{Float64}:
 1.00000e+00
 3.10192e-06
 1.04403e-11
 6.13953e-14
 4.44089e-16
 2.22045e-16
```
So the second iteration is somewhat better, but the iteration terminated after
five iterations in both cases.


There are more examples for this in [ctk:mparraysdocs](@cite).

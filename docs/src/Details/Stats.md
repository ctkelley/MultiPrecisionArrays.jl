# Harvesting Iteration Statistics

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

You may wonder why the residual after the first iteration was so much
larger than single precision roundoff. The reason is that the default 
when the low precision is single is to downcast the residual to single before 
the solve (onthefly=false).
One can enable interprecision transfers on the fly and see the difference.

```
julia> MPF2=mplu(A; onthefly=true);

julia> mpout2=\(MPF2, b; reporting=true);

julia> mpout2.rhist
5-element Vector{Float64}:
 9.99878e-01
 6.17721e-07
 3.84581e-13
 7.99361e-15
 8.88178e-16
```
So the second iteration is much better, but the iteration terminated after
four iterations in both cases.


There are more examples for this in [ctk:mparraysdocs](@cite).

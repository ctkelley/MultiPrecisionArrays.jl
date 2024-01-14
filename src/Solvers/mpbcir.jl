"""
mpbcir(AF::MPBFact, b; reporting=false, verbose=false, mpdebug=false)

BiCGSTAB-IR solver 

We overload the backslash operator to call mprmir for a multiprecision
MPBFact factorization. So if ```MPA``` is an MPB:Array and
```
AF = mpblu!(MPA)
```
Then ```AF\\b``` maps to
```
mpbcir(AF, b)
```
You should call ```mpgir``` explicitly if you want the iteration
statistics or use the kwargs for the backslash operator.

You can do the construction and factorization all at once with

```
MPF=mpblu(A)
```
Which is what I will do in the examples.

The three kwargs are ```reporting```, ```verbose```, and ```mpdebug```.
Of these ```reporting``` is most useful. The other two are for my
debugging pleasure.

When you do
```
mpout = mpbcir(NAF, b; reporting=true)
```
You get a structure where

```mpout.sol``` is the solution

```mpout.rhist``` is the residual history for IR

```mpout.khist``` is the krylovs/IR iteration 

Other parts of ```mpout``` are the high and low precisions 
```TH``` and  ```TL```.


## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> N=4096; A = I - 800.0 * Gmat(N); b=ones(N);

julia> MBF=mpblu(A);

julia> mbout=\\(MBF, b; reporting=true);

julia> x=mbout.sol; norm(b-A*x,Inf)
8.88534e-12

julia> mpout.rhist
5-element Vector{Float64}:
 1.00000e+00
 1.54060e-10
 8.43281e-12
 8.24296e-12
 8.88534e-12
# Stagnation after the fourth iteration and no progress after the second

julia> mpout.khist
4-element Vector{Int64}:
 3
 3
 3
 3
# 3 Krylovs per iteration. This is worse that GMRES because BiCGSTAB
# needs two matrix-vector products per iteration

julia> mpout.TH
Float64

julia> mpout.TL
Float32

# Repeat the experiment with low precision TL=Float16 (half)

julia> MPA=MPBArray(A; TL=Float16); AF=mpblu!(MPA); 

# You can use backslash too

julia> mpout=\\(AF, b; reporting=true);

julia> x=mpout.sol; norm(b-A*x,Inf)
9.34008e-12
# Residual is as good as the TL=Float32 case.

julia> mpout.rhist
5-element Vector{Float64}:
4-element Vector{Float64}:
 1.00000e+00
 1.11821e-11
 8.88811e-12
 9.34008e-12

# Stagnaton after 3 iterations

julia> mpout.khist
3-element Vector{Int64}:
 11
 11
 11
# No storage constraints.

julia> 

julia> mpout.TH
Float64

julia> mpout.TL
Float16


```

"""
function mpbcir(AF::MPBFact, b; reporting = false, 
                verbose = false, mpdebug = false)
     sol_out=mpkrir(AF, b; reporting=reporting,
                verbose = verbose, mpdebug = mpdebug)
    return sol_out
end

"""
mpgmir(AF::MPGFact, b; reporting=false, verbose=false, mpdebug=false)

GMRES-IR solver 

We overload the backslash operator to call mprmir for a multiprecision
MPGFact factorization. So if ```MPA``` is an MPArray and
```
AF = mpglu!(MPA)
```
Then ```AF\\b``` maps to
```
mpgmir(AF, b)
```
You should call ```mpgir``` explicitly if you want the iteration
statistics.

The three kwargs are ```reporting```, ```verbose```, and ```mpdebug```.
Of these ```reporting``` is most useful. The other two are for my
debugging pleasure.

When you do
```
mpout = mpgmir(NAF, b; reporting=true)
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

julia> MPA=MPArray(A); AF=mpglu!(MPA);

julia> mpout=mpgmir(AF, b; reporting=true);

julia> x=mpout.sol; norm(b-A*x,Inf)
8.92664e-12

julia> mpout.rhist
4-element Vector{Float64}:
 6.40000e+01
 3.32046e-09
 1.17624e-10
 1.17893e-10
# Stagnation after the second iteration

julia> mpout.khist
3-element Vector{Int64}:
 4
 5
 4
# 4-5 Krylovs per iteration.

julia> mpout.TH
Float64

julia> mpout.TL
Float32

# Repeat the experiment with low precision TL=Float16 (half)

julia> MPA=MPArray(A; TL=Float16); AF=mpglu!(MPA); 

# You can use backslash too

julia> mpout=\\(AF, b; reporting=true);

julia> x=mpout.sol; norm(b-A*x,Inf)
8.65075e-12
# Residual is as good as the TL=Float32 case.

julia> mpout.rhist
5-element Vector{Float64}:
 6.40000e+01
 2.00140e-03
 2.05307e-07
 1.16612e-10
 1.18166e-10

# Stagnaton after 3 iterations

julia> mpout.khist
5-element Vector{Int64}:
 10
 10
 10
 10
# The default basissize=10 so we are taking all the GMRES iterations we
# can at each iteration.

julia> 

julia> mpout.TH
Float64

julia> mpout.TL
Float16


```

"""
function mpgmir(AF::MPGFact, b; reporting = false, 
                verbose = false, mpdebug = false)
    sol_out=mpkrir(AF, b; reporting=reporting, 
                verbose = verbose, mpdebug = mpdebug)
    return sol_out
end

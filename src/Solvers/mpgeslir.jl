"""
mpgeslir(MPA::MPArray, b; reporting = false, verbose = true)

Use a multi-precision factorization to solve a linear system with
plain vanilla iterative refinement.

This version is analogous to ```A\\b``` and combines the factorization
and the solve. You start with MPA=MPArray(A) and then pass MPA
to mpgeslir and combine the factorization and the solve. 

You can also get the multiprecision factorization directly with
```
MPF=mplu!(A)
```
and then pass ```MPF``` to mpgeslir.

Unlike lu, this does overwrite the low precision part of MPA.
I use this to get some timing results and it's also convenient
if you want to do factor and solve in one statement. 

You can also get this with ```x = MPA\\b```.

If you set the kwarg ```reporting``` to true you can get the IR
residual history. The output of 
```
x = MPA\\b
```
or
```
x=MPF\\b
```
is the solition. The output of 
```
mout = \\(MPA,b; reporting=true)
```
or
```
mout = \\(MPF,b; reporting=true)
```
is a structure. ```mpout.sol``` is the solution. ```mpout.rhist```
is the residual history. mpout also contains the datatypes TH for
high precision and TL for low precision.

## Example
```jldoctest
julia> using MultiPrecisionArrays.Examples

julia> N=4096; A = I - 800.0 * Gmat(N); b=ones(N);

julia> MPF=mplu(A);

julia> mout=\\(MPF, b; reporting=true);

julia> mout.rhist
6-element Vector{Float64}:
 1.00000e+00
 5.36483e-02
 1.57977e-05
 5.10232e-09
 7.76756e-12
 9.90008e-12

# Stagnation after four IR iterations

julia> [mout.TH mout.TL]
1×2 Matrix{DataType}:
 Float64  Float32

```
"""
function mpgeslir(MPA::MPArray, b; reporting = false, verbose = true)
# Factor MPA and return Factorization object
MPF=mplu!(MPA);
# Call mpgeslir for the solve
xi=\(MPF, b; reporting=reporting, verbose=verbose)
return xi
end

"""
mpgeslir(AF::MPFact, b; reporting=false, verbose=true)

Use a multi-precision factorization to solve a linear system with
plain vanilla iterative refinement.

MPFact is a union of all the MultiPrecision factorizations in the package. 
The triangular solver will dispatch on the various types depending on
how the interprecision transfers get done.
"""
function mpgeslir(AF::MPFact, b; reporting = false, verbose = true)
    #
    # What kind of problem are we dealing with?
    #
    mpdebug = false
    normtype = Inf
    TB = eltype(b)
    MPStats = getStats(AF)
    TL = MPStats.TL
    TH = MPStats.TH
    r = AF.residual
    onthefly=AF.onthefly
    #
    # TFact is the precision of the factors
    #
    TFact = MPStats.TFact
    #
    # Are the precisions consistent? If not, I have a bug somewhere.
    # Otherwise, set the tolerance on the iteration to 10*eps.
    # If the iteration can't meet the tolerance, terminate when
    # the residual norms stagnate (res_old > .9 res_new)
    #
    (TH == TB) || error("inconsistent precisions")
    tolf = eps(TH)*TH.(10.0)
    #
    # Keep the records and accumulate the statistics. 
    #
    Meth = MPStats.Meth
    verbose && println(
        Meth,
        ": High precision = $TH, Low precision = $TL, Factorization storage precision = $TFact",
    )
    #
    # Showtime!
    #
    AD = AF.AH
    bnrm = norm(b, normtype)
    bsc = b
    AFS = AF.AF
    bS = TFact.(bsc)
    #
    # Initialize the iteration. I initialize to zero. That makes the
    # iteration count the same as the high precision matvec and the 
    # triangular sovles
    #
    x = zeros(TB, size(b))
    #
    # Initial residual
    #
    r .= b
    tol = tolf * bnrm
    rs = bS
#
#
    rhist = Vector{Float64}()
    rnrm = norm(r, normtype)
    rnrmx = rnrm * TB(2.0)
    oneb = TB(1.0)
    itc = 0
    #
    # Put initial residual norm into the history and iterate.
    #
    push!(rhist, rnrm)
    while (rnrm > tol) && (rnrm <= .9*rnrmx)
        #
        # Scale the residual
        #
        r ./= rnrm
        #
        # Use the low-precision factorization
        #
        r .= IRTriangle!(AF, r, rs, verbose)
        #
        # Undo the scaling
        #
        r .*= rnrm
        #
        # Update the solution and residual
        #
        x .+= r
        mul!(r, AD, x)
        r .*= -oneb
        axpy!(oneb, bsc, r)
        rnrmx = rnrm
        rnrm = norm(r, normtype)
        itc += 1
        push!(rhist, rnrm)
        mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol")
        #
        # If the residual norm increased, complain.
        #
#        complain_resid = (rnrm >= rnrmx) && (rnrm > 1.e3 * tol)
#        complain_resid && println("IR Norm increased: $rnrm, $rnrmx, $tol")
    end
    verbose && println("Residual history = $rhist")
    if reporting
        return (rhist = rhist, sol = x, TH = TH, TL = TL, TFact = TFact)
    else
        return x
    end
end

function getTL(AF::MPFact)
    TL = eltype(AF.AL)
    if is_heavy(AF)
    TFact = eltype(AF.AH)
    else
    TFact = eltype(AF.AL)
    end
    return (TL, TFact)
end

function getStats(AF)
    TH = eltype(AF.AH)
    (TL, TFact) = getTL(AF)
    MPStats = MPIRStats(TH, TL, TFact)
    return MPStats
end

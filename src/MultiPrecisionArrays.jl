module MultiPrecisionArrays

using LinearAlgebra: LinearAlgebra, LAPACK, LU, RowMaximum
using LinearAlgebra: Factorization, SingularException
using LinearAlgebra: axpy!, ldiv!, lu!, mul!, norm, opnorm

using SparseArrays: SparseArrays
using SIAMFANLEquations: SIAMFANLEquations, kl_bicgstab, kl_gmres
using SIAMFANLEquations: kstore
#using FLoops: @floop
using OhMyThreads: tforeach, @tasks
using Base.Threads: nthreads

#
# This is MultiPrecisionArrays.jl
# The package has data structures and algorithms to manage several
# variations on iterative refinement.
#

include("Structs4MP/MPBase.jl")
include("Structs4MP/MPArray.jl")
include("Structs4MP/MPGArray.jl")
include("Structs4MP/MPBArray.jl")
#
# MPHArrays may go away soon. They are only used for CI
include("CI_Only/MPHeavy.jl")
is_heavy(x::MPHFact) = true
is_heavy(x::MPGHFact) = true
is_heavy(x::MPLFact) = false
is_heavy(x::MPBFact) = false
is_heavy(x::MPGEFact) = false

MPSFact = Union{MPLFact,MPHFact,MPGEFact,MPGHFact}
MPFact = Union{MPLFact,MPHFact}
MPGFact = Union{MPGEFact,MPGHFact}
MPLFacts = Union{MPLFact}
MPKFact = Union{MPGFact,MPBFact}
MPHHFact = Union{MPGHFact,MPHFact}
#is_heavy(x::MPHHFact) = true
#
# Termination criteria defaults
#
const residtermdefault = true
#
# The termparms live in the module and termparms 
# is therefore global. Don't change this outside
# of the main thread.
#
struct TERM
    Cr::Real
    Ce::Real
    Rmax::Real
    litmax::Int
end
const Rmax_default = 0.5
const Cr_default = 1.0
const Ce_default = 1.0
const litmax_default = 10
term_parms_default=TERM(Cr_default, Ce_default, Rmax_default, litmax_default)


export termdata, term_parms_default, term_parms, update_parms
export update_parms, TERM


kmeth(x::MPGFact) = "GMRES"
kmeth(x::MPBFact) = "BiCGSTAB"

#export kmeth

import Base.eltype

function eltype(MP::Union{MPArray,MPHArray,MPGArray,MPBArray})
    TP = eltype(MP.AH)
    return TP
end

import Base.\
function \(
    AF::MPFact,
    b;
    verbose = false,
    reporting = false,
    term_parms = term_parms_default,
)
    xi = mpgeslir(AF, b; verbose = verbose, reporting = reporting, term_parms = term_parms)
    return xi
end

function \(
    AF::MPGFact,
    b;
    verbose = false,
    reporting = false,
    mpdebug = false,
    term_parms = term_parms_default,
)
    xi = mpkrir(
        AF,
        b;
        verbose = verbose,
        reporting = reporting,
        mpdebug = false,
        term_parms = term_parms,
    )
    return xi
end

function \(
    MPA::Union{MPArray},
    b;
    verbose = false,
    reporting = false,
    term_parms = term_parms_default,
)
    xi = mpgeslir(MPA, b; verbose = verbose, reporting = reporting, term_parms = term_parms)
    return xi
end

function \(
    AF::MPBFact,
    b;
    verbose = false,
    reporting = false,
    term_parms = term_parms_default,
)
    xi = mpkrir(AF, b; verbose = verbose, reporting = reporting, term_parms = term_parms)
    return xi
end


#
# hlu and hlu! are hacks of the Julia generic_lu source. I've used
# Polyester.jl to thread them and done some tweaking. On
# my Mac M2 Pro with 8 performance cores I'm seeing a 10x speedup
# over where I was before and seem to be using the hardware support
# for half precision as well as I can. I am still 5x slower than
# the LAPACK double precision LU, so there's a factor of 20 out there
# somewhere.
#
export hlu!
export hlu
#
# I'm exporting the multiprecision factorizations so you can pass
# them to nonlinear solvers and eigen solvers.
#
export mplu!
export mplu
export mphlu!
export mpblu
export mpblu!
export mpglu!
export mpghlu!
export mpglu
export Types_IR_Init
#export MPStats
#export mpqr!
#export mpcholesky!
#
# The solvers are mpgeslir (iterative refinement) and mpkrir (IR-Krylov).
# I'm not working on more than that right now. I have overloaded
# \ with these so you should not have to call them directly.
#
#export mpgeslir
#export mpgmir
#export mpkrir
#export mpbcir
#
# Each MPArray data structure comes with a structure to store a factorization.
# The differences are whether one does on-the-fly interprecision transfers
# of not. For plain IR with high=double and low=single, I think the answer 
# is clear (NO) and you should use MPArray with onthefly = false 
# (the default).
# and MPLEFact. If low precision is half, it's not so clear and the 
# documentation has an example to illustrate that.
#
# The factorization structures should be invisible to most people 
# and I may stop exporting them. 
#
# For IR-GMRES, it's more subtle.  The cost of Heavy IR with MPHArray
# and MPGHFact is an extra high precision matrix. I only use this for CI
# and am thinking about removing this stuff. 
#
# The data structures for statistics are almost surely not 
# interesting to anyone but me. Look at "Harvesting Iteration Statistics"
# in the docs for details. 
#

include("Solvers/mpkrir.jl")
include("Solvers/mpgeslir.jl")
include("Solvers/IRTriangle.jl")
include("Solvers/termination_settings.jl")
include("Solvers/update_parms.jl")
include("Solvers/Solvers_IR_Init.jl")
include("Solvers/Resid_IR.jl")
include("Solvers/IRKsolve.jl")
#include("Structs4MP/MPStats.jl")
include("Factorizations/hlu!.jl")
include("Factorizations/mplu!.jl")
include("Factorizations/mpglu!.jl")
include("Factorizations/mpblu!.jl")
include("CI_Only/mpghlu!.jl")

module Examples
using MultiPrecisionArrays
using LinearAlgebra: LinearAlgebra
using Reexport;
@reexport import LinearAlgebra.I

export Gmat

include("Examples/Gmat.jl")

end

end # module MPArrays

module MultiPrecisionArrays

using LinearAlgebra
using LinearAlgebra.BLAS
using SparseArrays
using SIAMFANLEquations
using Polyester

#
# This is MultiPrecisionArrays.jl
# The package has data structures and algorithms to manage several
# variations on iterative refinement.
#

include("Structs4MP/MPBase.jl")
include("Structs4MP/MPArray.jl")
include("Structs4MP/MPGArray.jl")
include("Structs4MP/MPHeavy.jl")

MPFact = Union{MPLFact,MPHFact}
MPLFacts = Union{MPLFact}
MPGFact = Union{MPGEFact, MPGHFact}

export MPGFact

#MPIRArray = Union{MPArray,MPHArray}


is_heavy(x::MPHFact) = true
is_heavy(x::MPLFact) = false




import Base.eltype

function eltype(MP::Union{MPArray,MPHArray})
    TP = eltype(MP.AH)
    return TP
end

#function eltype(MPH::MPHArray)
#    TP = eltype(MPH.AH)
#    return TP
#end

import Base.\
function \(AF::MPFact, b; verbose = false, reporting = false)
    xi = mpgeslir(AF, b; verbose = verbose, reporting = reporting)
    return xi
end

function \(AF::MPGFact, b; verbose = false, reporting = false, mpdebug=false) 
    xi = mpgmir(AF, b; verbose = verbose, reporting = reporting, mpdebug=mpdebug)
    return xi
end

function \(MPA::Union{MPArray}, b; verbose=false, reporting=false)
          xi = mpgeslir(MPA, b; verbose = verbose, reporting = reporting)
          return xi
end

#
# hlu and hlu! are hacks of the Julia generic_lu source. I've used
# Polyester.jl to thread them and done some tweaking. On
# my Mac M2 Pro with 8 performace cores I'm seeing a 10x speedup
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
export mpglu
export mpglu!
export mpqr!
export mpcholesky!
#
# The solvers are mpgeslir (iterative refinement) and mpgmir (IR-GMRES).
# I'm not working on more than that right now. I have overloaded
# \ with these so you should not have to call them directly unless
# you are looking at iteration statistics
#
export mpgeslir
export mpgmir
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
# and MPGHFact is an extra high precision matrix. If you can afford the 
# storage and communication burden, it's a reasonable thing to do. 
# If you can't, on-the-fly is your only option. The differences in time
# become less significant as the problem size gets large and the O(N^2)
# interprecision transfer cost is completely dominated by the O(N^3) 
# factorization cost. 
#
export MPArray
export MPHArray
#
#
#
export MPLFact
export MPGHFact
export MPGEFact
export MPFact
#
#
#
#export MPEArray
export MPFArray
export MPGArray
export MPHFact
export MPhatv
export MPhptv
#
# The data structres for statistics are almost surely not 
# interesting to anyone but me. I will document how one can
# use the solvers to get the stats sooner or later.
#
export MPGStats
export MPIRStats

include("Solvers/mpgmir.jl")
include("Solvers/mpgeslir.jl")
include("Solvers/IRTriangle.jl")
include("Structs4MP/MPStats.jl")
include("Factorizations/hlu!.jl")
include("Factorizations/mplu!.jl")
include("Factorizations/mpglu!.jl")

module Examples
using MultiPrecisionArrays
using LinearAlgebra
using LinearAlgebra.BLAS

export Gmat

include("Examples/Gmat.jl")

end

end # module MPArrays

module MultiPrecisionArrays

using LinearAlgebra 
using SparseArrays
using SIAMFANLEquations

#
# This is MPArrays.jl
# The package has data structures and algorithms to manage several
# variations on iterative refinement.
#

include("Structs4MP/MPLight.jl")
include("Structs4MP/MPHeavy.jl")

MPIRArray=Union{MPArray,MPHArray}


MPFact=Union{MPLFact, MPLEFact, MPHFact}

MPLFacts=Union{MPLFact, MPLEFact}

import Base.eltype
function eltype(MP::MPArray)
TP=eltype(MP.AH)
return TP
end

function eltype(MPH::MPHArray)
TP = eltype(MPH.AH)
return TP
end

import Base.\
function \(AF::MPLFact, b; verbose=false, reporting=false)
xi = mpgesl2(AF,b; verbose=verbose, reporting=reporting)
return xi
end

function \(AF::MPLEFact, b; verbose=false, reporting=false)
xi = mpgesl2(AF,b; verbose=verbose, reporting=reporting)
return xi
end

function \(AF::MPHFact, b; verbose=false, reporting=false)
xi = mpgesl2(AF,b; verbose=verbose, reporting=reporting)
return xi
end

function \(AF::MPGHFact, b; verbose=false, reporting=false)
xi = mpgmir(AF,b; verbose=verbose, reporting=reporting)
return xi
end

#
# I'm exporting the multiprecision factorizations so you can pass
# them to nonlinear solvers and eigen solvers.
#
export mplu!
export mphlu!
export mpglu!
export mpqr!
export mpcholesky!
#
# The solvers are mpgesl2 (iterative refinement) and mpgmir (IR-GMRES).
# I'm not working on more than that right now. I have overloaded
# \ with these so you should not have to call them directly unless
# you are looking at iteration statistics
#
export mpgesl2
export mpgmir
#
# Each MPArray data structure comes with a structure to store a factorization.
# The ones you care about, unless you are a real geek, are
# MPArray and MPLFact for IR and MPHArray and MPGHFact for IR-GMRES
#
export MPArray
export MPHArray
#
#
#
export MPLFact
export MPGHFact
#
#
#
export MPEArray
export MPFArray
export MPHFact
export MPLEFact
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
include("Solvers/mpgesl2.jl")
include("Solvers/IRTriangle.jl")
include("Structs4MP/MPStats.jl")

module Examples
export Gmat

include("Examples/Gmat.jl")

end

end # module MPArrays

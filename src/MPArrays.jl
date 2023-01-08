module MPArrays

using LinearAlgebra 
using SparseArrays
using SIAMFANLEquations

#
# This is MPArrays.jl
# The package has data structures and algorithms to manage several
# variations on iterative refinement.
#

include("MPStructs/MPLight.jl")
include("MPStructs/MPHeavy.jl")

MPIRArray=Union{MPArray,MPHArray}

#MPFact=Union{MPLFact, MPLEFact, MPHFact, MPGFact}

MPFact=Union{MPLFact, MPLEFact, MPHFact}

import Base.eltype
function eltype(MP::MPArray)
TP=eltype(MP.AH)
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

export mplu!
export mphlu!
export mpglu!
export mpqr!
export mpcholesky!
export MPArray
export MPEArray
export MPFArray
export MPHArray
export MPLFact
export MPLEFact
export MPHFact
export mpgesl2
export MPhatv
export MPhptv

export MPGStats
export MPIRStats

include("Solvers/mpgesl2.jl")
include("Solvers/IRTriangle.jl")
include("MPAStats.jl")

module Examples
export Gmat

include("Examples/Gmat.jl")

end

end # module MPArrays

module MPArrays

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


#MPFact=Union{MPLFact, MPLEFact, MPHFact, MPGFact}

MPFact=Union{MPLFact, MPLEFact, MPHFact}

MPLFacts=Union{MPLFact, MPLEFact}

import Base.eltype
function eltype(MP::MPArray)
TP=eltype(MP.AH)
return TP
end

#function eltype(MPH::MPHArray)
#TP = eltype(MPH.AH)
#return TP
#end

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

#function \(AF::MPGHFact, b; verbose=false, reporting=false)
#xi = mpgmir(AF,b; verbose=verbose, reporting=reporting)
#return xi
#end

export mplu!
export mphlu!
export mpglu!
#export mpglu!
export mpqr!
export mpcholesky!
export mpgmir
export MPArray
export MPEArray
export MPFArray
export MPHArray
export MPLFact
export MPLEFact
export MPHFact
export MPGHFact
export mpgesl2
export MPhatv
export MPhptv

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

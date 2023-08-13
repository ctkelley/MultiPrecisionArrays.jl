struct MPHArray
    AH::Array
    AStore::Array
    AL::Array
end

struct MPHFact
    AH::Array
    AStore::Array
    AL::Array
    AF::Factorization
end

#
# Heavy factorization for GMRES-IR
#
struct MPGHFact
    AH::Array
    AStore::Array
    AL::Array
    AF::Factorization
end


function MPHArray(AH::Array{Float64,2}; TL = Float32)
    AStore = copy(AH)
    AL = TL.(AH)
    MPH = MPHArray(AH, AStore, AL)
end

function MPHArray(AH::Array{Float32,2}; TL = Float16)
    AStore = copy(AH)
    AL = TL.(AH)
    MPH = MPHArray(AH, AStore, AL)
end

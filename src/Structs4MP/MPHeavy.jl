#
# MPEArray factorization for GMRES-IR
#
struct MPGEFact
    AH::AbstractArray
    AL::AbstractArray
    AF::Factorization
    VStore:: Array
    KStore:: NTuple
    residual::Vector
    onthefly::Bool
end

struct MPHArray
    AH::AbstractArray
    AStore::AbstractArray
    AL::AbstractArray
    residual::Vector
    onthefly::Bool
end

struct MPHFact
    AH::AbstractArray
    AL::AbstractArray
    AF::Factorization
    residual::Vector
    onthefly::Bool
end

#
# Heavy factorization for GMRES-IR
#
struct MPGHFact
    AH::AbstractArray
    AL::AbstractArray
    AF::Factorization
    VStore:: Array
    KStore:: NTuple
    residual::Vector
    onthefly::Bool
end

function MPHArray(AH::AbstractArray{Float64,2}; TL = Float32)
    AStore = copy(AH)
    AL = TL.(AH)
    (m,n)=size(AH); res=ones(eltype(AH),n)
    MPH = MPHArray(AH, AStore, AL, res, true)
end

function MPHArray(AH::AbstractArray{Float32,2}; TL = Float16)
    AStore = copy(AH)
    AL = TL.(AH)
    (m,n)=size(AH); res=ones(eltype(AH),n)
    MPH = MPHArray(AH, AStore, AL, res, true)
end

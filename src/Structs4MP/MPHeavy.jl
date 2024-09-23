#
# MPEArray factorization for GMRES-IR
# This stuff is only used for CI
#
struct MPGEFact
    AH::AbstractArray
    AL::AbstractArray
    AF::Factorization
    VStore::Array
    KStore::NTuple
    residual::Vector
    sol::Vector
    onthefly::Bool
    residterm::Bool
end

struct MPHArray
    AH::AbstractArray
    AStore::AbstractArray
    AL::AbstractArray
    residual::Vector
    sol::Vector
    onthefly::Bool
end

struct MPHFact
    AH::AbstractArray
    AL::AbstractArray
    AF::Factorization
    residual::Vector
    sol::Vector
    onthefly::Bool
    residterm::Bool
end

#
# Heavy factorization for GMRES-IR
#
struct MPGHFact
    AH::AbstractArray
    AL::AbstractArray
    AF::Factorization
    VStore::Array
    KStore::NTuple
    residual::Vector
    sol::Vector
    onthefly::Bool
    residterm::Bool
end

function MPHArray(AH::AbstractArray{Float64,2}; TF = Float32)
    AStore = copy(AH)
    AL = TF.(AH)
    (m, n) = size(AH)
    res = ones(eltype(AH), n)
    sol = ones(eltype(AH), n)
    MPH = MPHArray(AH, AStore, AL, res, sol, true)
end

function MPHArray(AH::AbstractArray{Float32,2}; TF = Float16)
    AStore = copy(AH)
    AL = TF.(AH)
    (m, n) = size(AH)
    res = ones(eltype(AH), n)
    sol = ones(eltype(AH), n)
    MPH = MPHArray(AH, AStore, AL, res, sol, true)
end

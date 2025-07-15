# This file is MPBase.jl 
# Plain vanilla iterative refinement. The structures store
# High and low precision arrays and factorizations of these arrays
#
struct MPArray{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    residual::Vector{TR}
    sol::Vector{TW}
    onthefly::Bool
end

# Factorizations. All in one place for now as I tweak the structure

struct MPLFact{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat,ATF<:Factorization}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    AF::ATF
    residual::Vector{TR}
    sol::Vector{TW}
    onthefly::Bool
    residterm::Bool
    anrm::AbstractFloat
end

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
    anrm::AbstractFloat
end

struct MPBFact
    AH::AbstractArray
    AL::AbstractArray
    AF::Factorization
    VStore::Vector
    KStore::NTuple
    residual::Vector
    sol::Vector
    onthefly::Bool
    residterm::Bool
    anrm::AbstractFloat
end

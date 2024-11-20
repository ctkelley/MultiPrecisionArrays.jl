# This file is MPBase.jl 
# Plain vanilla iterative refinement. The structures store
# High and low precision arrays and factorizations of these arrays
#
struct MPArray{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    residual::Vector{TR}
    sol::Vector{TR}
    onthefly::Bool
end

struct MPLFact{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat,ATF<:Factorization}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    AF::ATF
    residual::Vector{TR}
    sol::Vector{TR}
    onthefly::Bool
    residterm::Bool
    anrm::TF
end

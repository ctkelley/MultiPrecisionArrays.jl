# This file is MPBase.jl 
# Plain vanilla iterative refinement. The structures store
# High and low precision arrays and factorizations of these arrays
#
struct MPArray{TH<:AbstractFloat,TF<:AbstractFloat}
    AH::AbstractArray{TH,2}
    AL::AbstractArray{TF,2}
    residual::Vector{TH}
    onthefly::Bool
end

struct MPLFact{TH<:AbstractFloat,TF<:AbstractFloat,ATF<:Factorization}
    AH::AbstractArray{TH,2}
    AL::AbstractArray{TF,2}
    AF::ATF
    residual::Vector{TH}
    onthefly::Bool
end

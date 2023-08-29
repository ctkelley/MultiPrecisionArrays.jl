#
# This file is MPArray.jl 
# Plain vanilla iterative refinement. The structures store
# High and low precision arrays
# Factorizations of these arrays
#
struct MPArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::Array{TH,2}
    AL::Array{TL,2}
    residual::Vector{TH}
end
#
# The constructors for the multi-precision arrays
#
"""
MPArray(AH::Array{Float64,2}; TL = Float32, onthefly=false)
Default constructor for MPArray. The difference between
and MPArray and an MPEArray is that the MPEArray does interprecision
transfres on the fly.

The MPArray data structure is

struct MPEArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::Array{TH,2}
    AL::Array{TL,2}
    residual::Vector{TH}
end

MPEArray is exactly the same but the triangular solver dispatches
differently.
"""
function MPArray(AH::Array{Float64,2}; TL = Float32, onthefly=false)
    AL = TL.(AH)
    (m,n)=size(AH); res=ones(eltype(AH),n)
    onthefly ?  MPA = MPEArray(AH, AL, res) : MPA = MPArray(AH, AL, res)
end


struct MPLFact{TH<:AbstractFloat,TL<:AbstractFloat,TF<:Factorization}
    AH::Array{TH,2}
    AL::Array{TL,2}
    AF::TF
    residual::Vector{TH}
end

#
# The MPAE and MPLEFact structures tell the solver to 
# do the triangular solve in high precision. The E stands
# for expensive
#

struct MPEArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::Array{TH,2}
    AL::Array{TL,2}
    residual::Vector{TH}
end


struct MPLEFact{TH<:AbstractFloat,TL<:AbstractFloat,TF<:Factorization}
    AH::Array{TH,2}
    AL::Array{TL,2}
    AF::TF
    residual::Vector{TH}
end


function MPArray(AH::Array{Float32,2}; TL = Float16, onthefly=false)
    AL = TL.(AH)
    (m,n)=size(AH); res=ones(eltype(AH),n)
    onthefly ?  MPA = MPEArray(AH, AL, res) : MPA = MPArray(AH, AL, res)
end

function MPEArray(AH::Array{Float32,2}; TL = Float16)
    AL = TL.(AH)
    (m,n)=size(AH); res=ones(eltype(AH),n)
    MPA = MPEArray(AH, AL, res)
end

function MPEArray(AH::Array{Float64,2}; TL = Float32)
    (m,n)=size(AH); res=ones(eltype(AH),n)
    AL = TL.(AH)
    MPA = MPEArray(AH, AL, res)
end

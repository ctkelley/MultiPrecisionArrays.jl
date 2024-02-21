module IR3P

using LinearAlgebra
using LinearAlgebra.BLAS
using SparseArrays
using SIAMFANLEquations
using Printf
using Polyester

#
# Test of three precision IR
#

struct M2EPArray{TH<:AbstractFloat,TE<:AbstractFloat}
    AH::AbstractArray{TH,2}
    AH2::AbstractArray{TH,2}
    residual::Vector{TE}
    defect::Vector{TE}
    sole::Vector{TE}
end

function M2EPArray(AH::AbstractArray{TH,2}; TE=Float64) where TH
AH2 = copy(AH)
(m,n) = size(AH)
res=zeros(TE,n)
def=zeros(TE,n)
sole=zeros(TE,n)
M2EPA = M2EPArray(AH, AH2, res, def, sole)
end

include("tir3p.jl")
include("wilk.jl")


export tir3p
export wilk

end

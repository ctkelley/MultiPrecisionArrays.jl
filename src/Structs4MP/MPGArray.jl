"""
MPGArray(AH::AbstractArray{Float64,2}; TF = Float32, basissize=10)
Default constructor for MPGArray. Allocate the storage for 
GMRES-IR

C. T. Kelley 2023


The MPGArray data structure is

```
struct MPGArray{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    VStore::AbstractArray{TR,2}
    KStore::NTuple
    residual::Vector{TR}
    sol::Vector{TR}
    onthefly::Bool
end
```
The constructor just builds an MPGArray with TW=Float64. Set TF=Float16
to get double/half IR.
"""

struct MPGArray{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    VStore::AbstractArray{TR,2}
    KStore::NTuple
    residual::Vector{TR}
    sol::Vector{TR}
    onthefly::Bool
end

"""
MPGArray(AH::AbstractArray{Float64,2}; basissize=10, TF=Float32, TR=nothing)

An MPGArray stores the high precision matrix, the low precision factorization,
the Krylov basis, and a few other things GMRES needs. If the high precision
matrix is double, the low precision is single by default. Half is an option
which you get with TF=Float16.
"""
function MPGArray(AH::AbstractArray{Float64,2}; TR=nothing,
                 basissize=10, TF=Float32)
AL=TF.(AH)
(m,n)=size(AH)
(TR == nothing) ? TRR=eltype(AH) : TRR = TR
res=ones(TRR,n)
sol=ones(TRR,n)
VStore=zeros(TRR,n,basissize)
KStore=KIRstore(n,"gmres")
MPGA=MPGArray(AH, AL, VStore, KStore, res, sol, true)
end


"""
MPGArray(AH::AbstractArray{Float32,2}; basissize=10, TF=Float16, TR=nothing)

An MPGArray stores the high precision matrix, the low precision factorization,
the Krylov basis, and a few other things GMRES needs. Since High precision is 
single, low is half. I'm leaving the kwarg for TF in there because it makes
is easier to cut/paste calls to MPGArray different precisions into a CI loop.
"""
function MPGArray(AH::AbstractArray{Float32,2}; basissize=10, TF=Float16,
                  TR=nothing)
AL=TF.(AH)
(m,n)=size(AH)
(TR == nothing) ? TRR=eltype(AH) : TRR = TR
res=ones(TRR,n)
sol=ones(TRR,n)
VStore=zeros(TRR,n,basissize)
KStore=KIRstore(n,"gmres")
MPGA=MPGArray(AH, AL, VStore, KStore, res, sol, true)
return MPGA
end

"""     
KIRstore(n, lsolver, TR=Float64)

Preallocates the vectors a Krylov method uses internally.
"""
function KIRstore(n, lsolver, TR=Float64)
    tmp1 = zeros(TR,n)
    tmp2 = zeros(TR,n)
    tmp3 = zeros(TR,n)
    tmp4 = zeros(TR,n)
    if lsolver == "gmres"
        return (tmp1, tmp2, tmp3, tmp4)
    else
        tmp5 = zeros(TR,n)
        tmp6 = zeros(TR,n)
        tmp7 = zeros(TR,n)
        return (tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)
    end
end


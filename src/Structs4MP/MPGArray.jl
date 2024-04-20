"""
MPGArray(AH::AbstractArray{Float64,2}; TF = Float32, basissize=10)
Default constructor for MPGArray. Allocate the storage for 
GMRES-IR

C. T. Kelley 2023


The MPGArray data structure is

```
struct MPGArray{TW<:AbstractFloat,TF<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    VStore::AbstractArray{TW,2}
    KStore::NTuple
    residual::Vector{TW}
    onthefly::Bool
end
```
The constructor just builds an MPGArray with TW=Float64. Set TF=Float16
to get double/half IR.
"""

struct MPGArray{TW<:AbstractFloat,TF<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    VStore::AbstractArray{TW,2}
    KStore::NTuple
    residual::Vector{TW}
    onthefly::Bool
end

"""
MPGArray(AH::AbstractArray{Float64,2}; basissize=10, TF=Float32)

An MPGArray stores the high precision matrix, the low precision factorization,
the Krylov basis, and a few other things GMRES needs. If the high precision
matrix is double, the low precision is single by default. Half is an option
which you get with TF=Float16.
"""
function MPGArray(AH::AbstractArray{Float64,2}; basissize=10, TF=Float32)
AL=TF.(AH)
(m,n)=size(AH)
res=ones(eltype(AH),n)
VStore=zeros(eltype(AH),n,basissize)
KStore=kstore(n,"gmres")
MPGA=MPGArray(AH, AL, VStore, KStore, res, true)
end


"""
MPGArray(AH::AbstractArray{Float32,2}; basissize=10, TF=Float16)

An MPGArray stores the high precision matrix, the low precision factorization,
the Krylov basis, and a few other things GMRES needs. Since High precision is 
single, low is half. I'm leaving the kwarg for TF in there because it makes
is easier to cut/paste calls to MPGArray different precisions into a CI loop.
"""
function MPGArray(AH::AbstractArray{Float32,2}; basissize=10, TF=Float16)
AL=TF.(AH)
(m,n)=size(AH)
res=ones(eltype(AH),n)
VStore=zeros(eltype(AH),n,basissize)
KStore=kstore(n,"gmres")
MPGA=MPGArray(AH, AL, VStore, KStore, res, true)
return MPGA
end

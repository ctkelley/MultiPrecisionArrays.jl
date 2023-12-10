"""
MPGArray(AH::AbstractArray{Float64,2}; TL = Float32, basissize=10)
Default constructor for MPGArray. Allocate the storage for 
GMRES-IR

C. T. Kelley 2023


The MPGArray data structure is

```
struct MPGArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::AbstractArray{TH,2}
    AL::AbstractArray{TL,2}
    VStore::AbstractArray{TH,2}
    KStore::NTuple
    residual::Vector{TH}
    onthefly::Bool
end
```
The constructor just builds an MPGArray with TH=Float64. Set TL=Float16
to get double/half IR.
"""

struct MPGArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::AbstractArray{TH,2}
    AL::AbstractArray{TL,2}
    VStore::AbstractArray{TH,2}
    KStore::NTuple
    residual::Vector{TH}
    onthefly::Bool
end

"""
MPGArray(AH::AbstractArray{Float64,2}; basissize=10, TL=Float32)

An MPGArray stores the high precision matrix, the low precision factorization,
the Krylov basis, and a few other things GMRES needs. If the high precision
matrix is double, the low precision is single by default. Half is an optioin
which you get with TL=Float16.
"""
function MPGArray(AH::AbstractArray{Float64,2}; basissize=10, TL=Float32)
AL=TL.(AH)
(m,n)=size(AH)
res=ones(eltype(AH),n)
VStore=zeros(eltype(AH),n,basissize)
KStore=kstore(n,"gmres")
MPGA=MPGArray(AH, AL, VStore, KStore, res, true)
end


"""
MPGArray(AH::AbstractArray{Float32,2}; basissize=10, TL=Float16)

An MPGArray stores the high precision matrix, the low precision factorization,
the Krylov basis, and a few other things GMRES needs. Since High precision is 
single, low is half. I'm leaving the kwarg for TL in there because it makes
is easier to cut/paste calls to MPGArray different precisions into a CI loop.
"""
function MPGArray(AH::AbstractArray{Float32,2}; basissize=10, TL=Float16)
AL=TL.(AH)
(m,n)=size(AH)
res=ones(eltype(AH),n)
VStore=zeros(eltype(AH),n,basissize)
KStore=kstore(n,"gmres")
MPGA=MPGArray(AH, AL, VStore, KStore, res, true)
return MPGA
end

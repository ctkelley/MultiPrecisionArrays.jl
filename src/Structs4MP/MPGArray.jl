"""
MPGArray(AH::Array{Float64,2}; TL = Float32, basissize=10)
Default constructor for MPGArray. Allocate the storage for 
GMRES-IR

C. T. Kelley 2023


The MPGArray data structure is

```
struct MPGArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::Array{TH,2}
    AL::Array{TL,2}
    VStore::Array{TH,2}
    KStore::NTuple
    residual::Vector{TH}
    onthefly::Bool
end
```
The constructor just builds an MPGArray with TH=Float64. Set TL=Float16
to get double/half IR.
"""

struct MPGArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::Array{TH,2}
    AL::Array{TL,2}
    VStore::Array{TH,2}
    KStore::NTuple
    residual::Vector{TH}
    onthefly::Bool
end

"""
MPGArray(AH::Array{Float64,2}; basissize=10, TL=Float32)
"""
function MPGArray(AH::Array{Float64,2}; basissize=10, TL=Float32)
AL=TL.(AH)
(m,n)=size(AH)
res=ones(eltype(AH),n)
VStore=zeros(eltype(AH),n,basissize)
KStore=kstore(n,"gmres")
MPGA=MPGArray(AH, AL, VStore, KStore, res, true)
end


"""
MPGArray(AH::Array{Float32,2}; basissize=10, TL=Float16)
"""
function MPGArray(AH::Array{Float32,2}; basissize=10, TL=Float16)
AL=TL.(AH)
(m,n)=size(AH)
res=ones(eltype(AH),n)
VStore=zeros(eltype(AH),n,basissize)
KStore=kstore(n,"gmres")
MPGA=MPGArray(AH, AL, VStore, KStore, res, true)
return MPGA
end

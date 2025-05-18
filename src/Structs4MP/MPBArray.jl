"""
MPBArray(AH::AbstractArray{Float64,2}; TF = Float32)
Default constructor for MPBArray. Allocate the storage for 
BiCGSTAB-IR

C. T. Kelley 2023


The MPBArray data structure is

```
struct MPBArray{TW<:AbstractFloat,TF<:AbstractFloat, TR<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    VStore::Vector
    KStore::NTuple
    residual::Vector{TR}
    sol::Vector{TW}
    onthefly::Bool
end
```
The constructor just builds an MPBArray with TW=Float64. Set TF=Float16
to get double/half IR.
"""

struct MPBArray{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    VStore::Vector
    KStore::NTuple
    residual::Vector{TR}
    sol::Vector{TW}
    onthefly::Bool
end

"""
MPBArray(AH::AbstractArray{Float64,2}; TF=Float32, TR=nothing)

An MPBArray stores the high precision matrix, the low precision factorization
and a few other things BiCGSTAB needs. If the high precision
matrix is double, the low precision is single by default. Half is an option
which you get with TF=Float16.
"""
function MPBArray(AH::AbstractArray{Float64,2}; TF = Float32, TR = nothing)
    AL = TF.(AH)
    TW = eltype(AH)
    (m, n) = size(AH)
    (TR == nothing) && (TR = TW)
    res = ones(TR, n)
    sol = ones(TW, n)
    VStore = copy(res)
    KStore = KIRstore(n, "bicgstab", TW)
    MPBA = MPBArray(AH, AL, VStore, KStore, res, sol, true)
end


"""
MPBArray(AH::AbstractArray{Float32,2}; TF=Float16, TR=nothing)

An MPBArray stores the high precision matrix, the low precision factorization
and a few other things BiCGSTAB needs. Since High precision is 
single, low is half. I'm leaving the kwarg for TF in there because it makes
is easier to cut/paste calls to MPBArray different precisions into a CI loop.
"""
function MPBArray(AH::AbstractArray{Float32,2}; TF = Float16, TR = nothing)
    AL = TF.(AH)
    (m, n) = size(AH)
    TW = eltype(AH)
    (TR == nothing) && (TR = TW)
    res = ones(TR, n)
    sol = ones(TW, n)
    VStore = copy(res)
    KStore = KIRstore(n, "bicgstab", TW)
    MPBA = MPBArray(AH, AL, VStore, KStore, res, sol, true)
    return MPBA
end

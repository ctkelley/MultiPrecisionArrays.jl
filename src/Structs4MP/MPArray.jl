"""
MPArray(AH::AbstractArray{Float64,2}; TF = Float32, onthefly=false)
Default constructor for MPArray. 

C. T. Kelley 2024

The MPArray data structure is

```
struct MPArray{TW<:AbstractFloat,TF<:AbstractFloat,TR<:AbstractFloat}
    AH::AbstractArray{TW,2}
    AL::AbstractArray{TF,2}
    residual::Vector{TR}
    sol::Vector{TR}
    onthefly::Bool
end
```
The constructor just builds an MPArray with TW=Float64. Set TF=Float16
to get double/half IR.
"""
function MPArray(AH::AbstractArray{Float64,2}; 
      TR=nothing, TF = Float32, onthefly=nothing)
    AL = TF.(AH)
    TH = eltype(AH)
# Default is interprecision on the fly if TF = Float32
    (onthefly==nothing) && (onthefly = (TF==Float16))
# If TR = nothing, that's a signal to set TR=TH
    (TR==nothing) ? TRR=TH : TRR=TR
    (m,n)=size(AH); res=ones(TRR,n); sol=zeros(TRR,n)
    MPA = MPArray(AH, AL, res, sol, onthefly) 
end
"""
MPArray(AH::AbstractArray{Float32,2}; 
             TR = nothing, TF = Float16, onthefly=true)
Default single precision constructor for MPArray with TF=Float16

If your high precision array is single, then your low precision
array is half (Duh!). 

We do the triangular
solves with on-the-fly interprecision transfer in this case because
the bit of extra accuracy makes a difference and, at least for now,
on-the-fly interprecision transfers are cheaper.

Data structures etc are the same as in the 
double-single/half case, but you don't have the option to go lower than
half.
"""
function MPArray(AH::AbstractArray{Float32,2}; 
               TR=nothing, TF = Float16, onthefly=true)
    AL = TF.(AH)
    TH = eltype(AH)
# If TR = nothing, that's a signal to set TR=TH
    (TR==nothing) ? TRR=TH : TRR=TR
    (m,n)=size(AH); res=ones(TRR,n); sol=zeros(TRR,n)
    MPA = MPArray(AH, AL, res, sol, onthefly)
end



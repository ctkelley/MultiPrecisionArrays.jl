"""
MPArray(AH::Array{Float64,2}; TL = Float32, onthefly=false)
Default constructor for MPArray. 

C. T. Kelley 2023

The difference between and MPArray and an MPEArray is that 
MPEArray does interprecision transfers on the fly. Set onthefly = true
and get an MPEArray.

The MPArray data structure is

```
struct MPArray{TH<:AbstractFloat,TL<:AbstractFloat}
    AH::Array{TH,2}
    AL::Array{TL,2}
    residual::Vector{TH}
end
```
The constructor just builds an MPArray with TH=Float64. Set TL=Float16
to get double/half IR.

MPEArray is exactly the same but the triangular solver dispatches
differently.
"""
function MPArray(AH::Array{Float64,2}; TL = Float32, onthefly=false)
    AL = TL.(AH)
    (m,n)=size(AH); res=ones(eltype(AH),n)
    onthefly ?  MPA = MPEArray(AH, AL, res) : MPA = MPArray(AH, AL, res)
end
"""
MPArray(AH::Array{Float32,2}; TL = Float16, onthefly=true)
Default single precision constructor for MPArray. 

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
function MPArray(AH::Array{Float32,2}; TL = Float16, onthefly=true)
    AL = TL.(AH)
    (m,n)=size(AH); res=ones(eltype(AH),n)
    onthefly ?  MPA = MPEArray(AH, AL, res) : MPA = MPArray(AH, AL, res)
end


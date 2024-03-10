# mpglu!: Factor a MPGArray and set it up for GMRES by allocating room for Krylov vectors etc
```@docs
mpglu!(MPGA::MPGArray)
mpglu!(MPG::MPGEFact, A::AbstractArray{TW,2}) where TW <: Real
```

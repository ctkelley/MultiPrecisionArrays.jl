# mpglu!: Factor a MPGArray and set it up for GMRES by allocating room for Krylov vectors etc
```@docs
mpglu!(MPGA::MPGArray)
mpglu!(MPG::MPGEFact, A::AbstractArray{TH,2}) where TH <: Real
```

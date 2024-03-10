# mpblu!: Factor a MPBArray and set it up for BiCGSTAB by allocating room for a few vectors
```@docs
mpblu!(MPBA::MPBArray)
mpblu!(MPG::MPBFact, A::AbstractArray{TW,2}) where TW <: Real
```

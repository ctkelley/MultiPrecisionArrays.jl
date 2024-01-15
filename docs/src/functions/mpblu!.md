# mpblu!: Factor a MPBArray and set it up for BiCGSTAB by allocating room for a few vectors
```@docs
mpblu!(MPBA::MPBArray)
mpblu!(MPG::MPBFact, A::AbstractArray{TH,2}) where TH <: Real
```

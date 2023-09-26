"""
IRTriangle!(AF::Union{MPLFact,MPLEFact,MPHFact}, r, rs, verbose)
This is the solve phase using the factorization object (AFS) you get
from a multiprecision LU factorization.

Solve for the defect by quering on_the_fly
to figure out if we can do the triangular solves entirely in low precision.

If on_the_fly(AF) == false, then demote the residual and solve in low
precision.

The solve overwrites the residual with the defect.
"""
function IRTriangle!(AF::Union{MPLFact,MPHFact}, r, rs, verbose)
    AFS = AF.AF
    on_the_fly=AF.onthefly
    if on_the_fly
        ldiv!(AFS, r)
    else
        TFact = eltype(AFS)
        TH = eltype(r)
        rs .= TFact.(r)
        ldiv!(AFS, rs)
        r .= TH.(rs)
    end
    return r
end

"""
IRTriangle!(AF::Union{MPLFact,MPHFact}, r, rs)
This is the solve phase using the factorization object (AFS) you get
from a multiprecision LU factorization.

Solve for the defect by querying on_the_fly
to figure out if we can do the triangular solves entirely in low precision.

If on_the_fly(AF) == false, then demote the residual and solve in low
precision.

The solve overwrites the residual with the defect if TR=TW.
"""
function IRTriangle!(AF::Union{MPLFact,MPHFact}, r, rs)
    AFS = AF.AF
    on_the_fly = AF.onthefly
    TR=eltype(r)
    TW=eltype(AF.AH)
    if on_the_fly && (TR == TW)
        ldiv!(AFS, r)
    else
        rs .= r
        ldiv!(AFS, rs)
        r .= rs
    end
    return r
end

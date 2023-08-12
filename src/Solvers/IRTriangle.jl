"""
IRTriangle!(AF::MPLFact, r, rs, verbose)

Solve for the defect by demoting the residual and using a low-precision
factorization. Triangular solves entirely in low precision.
"""
function IRTriangle!(AF::MPLFact, r, rs, verbose)
AFS = AF.AF
#MPStats = getStats(AF)
TFact = eltype(AFS)
TH = eltype(r)
rs .= TFact.(r)
ldiv!(AFS, rs)
r .= TH.(rs)
return r
end

"""
IRTriangle!(AF::MPLEFact, r, rs, verbose)

Solve for the defect using a low-precision factorization and the HIGH-PRECISION residual. Triangular solves entirely in high precision with expensive promotion from high to low for the triangular matrices with each axpy. This is only for CI and making points in papers.
"""
function IRTriangle!(AF::MPLEFact, r, rs, verbose)
AFS=AF.AF
ldiv!(AFS, r)
return r
end

"""
IRTriangle!(AF::MPHFact, r, rs, verbose)

If you want to use the low-precision factorization as a preconditioner, you are best advised to use heavy MPArrays and promote the low-precision factorization to high before doing the solve. The solve in this version is entirely in high.
"""
function IRTriangle!(AF::MPHFact, r, rs, verbose)
AFS=AF.AF
ldiv!(AFS, r)
return r
end

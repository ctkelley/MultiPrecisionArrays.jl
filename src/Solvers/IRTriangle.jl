function IRTriangle!(AF::MPLFact, r, rs, verbose)
expensive=AF.expensive
AFS = AF.AF
if expensive
ldiv!(AFS,r)
else
MPStats = getStats(AF)
TFact = eltype(AFS)
TH = eltype(r)
rs .= TFact.(r)
ldiv!(AFS, rs)
r .= TH.(rs)
end
return r
end

function IRTriangle!(AF::MPLEFact, r, rs, verbose)
AFS=AF.AF
ldiv!(AFS, r)
return r
end

function IRTriangle!(AF::MPHFact, r, rs, verbose)
AFS=AF.AF
ldiv!(AFS, r)
return r
end

function IRTriangle(AF::MPGFact, r, rs, verbose)
x0 = zeros(size(r))
eta = 1.e-4
V = AF.VH
kout = kl_gmres(x0, r, MPhatv, V, eta, MPhptv; pdata = AF, side = "left")
#
# Tell them more than they need to know?
#
if verbose
   itcc = itc + 1
   println("Krylov stats: Iteration $itcc ", kout.reshist, "  ", kout.idid)
end
r .= kout.sol
return r
end


function mpgesl2(AF::MPTest, b)
TH = eltype(b);
(TH == Float64) ? TL = Float32 : TL = Float16
(TH == Float64) ? tolf=1.e-13 : tolf=1.e-6
TL = eltype(AF.AL)
AD=AF.AH
bnrm=norm(b,Inf)
bsc=b
AFS=AF.AF
bS=TL.(bsc)
ldiv!(AFS,bS)
println(typeof(bS))
x = TH.(bS)
r = copy(x)
mul!(r,AD,x)
r .*= -1.0
axpy!(1.0, bsc, r)
tol = tolf * bnrm
rs=bS
rnrm=norm(r, Inf)
itc=0
while rnrm > tol
r ./=rnrm
rs .= TL.(r)
ldiv!(AFS,rs)
r .= TH.(rs)
r .*= rnrm
x .+= r
mul!(r,AD,x)
r .*= -1.0
axpy!(1.0, bsc, r)
rnrm=norm(r,Inf)
itc += 1
println(itc, "  ", rnrm, "  ",tol)
end
return x
end


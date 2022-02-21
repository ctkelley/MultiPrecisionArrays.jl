function mpgesl(AF::MPFArray, b)
AD=AF.AD
bnrm=norm(b,Inf)
bsc=b
AFS=AF.AFS
bS=Float32.(bsc)
ldiv!(AFS,bS)
x = Float64.(bS)
r = copy(x)
mul!(r,AD,x)
r .*= -1.0
axpy!(1.0, bsc, r)
tol = 1.e-13 * bnrm
rs=bS
rnrm=norm(r, Inf)
while rnrm > tol
r ./=rnrm
rs .= Float32.(r)
ldiv!(AFS,rs)
r .= Float64.(rs)
r .*= rnrm
x .+= r
mul!(r,AD,x)
r .*= -1.0
axpy!(1.0, bsc, r)
rnrm=norm(r,Inf)
end
return x
end


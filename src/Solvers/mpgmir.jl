"""
mpgmir(AF::MPGHFact, b; reporting=false, verbose=true)

Prototype GMRES-IR solver
"""
function mpgmir(AF::MPGHFact, b; reporting=false, verbose=true)
#
normtype=Inf
TB=eltype(b)
onetb=TB(1.0)
x=zeros(TB,size(b))
brnm=norm(b,normtype)
#
AFS=AF.AF
AH=AF.AH
Gdata=(AFS=AFS,AH=AH)
#
# Initialize GMRES-IR
#
r=copy(x)
mul!(r,AD,x)
r .*= -onetb
axpy!(onetb, bsc, r)
rnrm=norm(r,normtype)
rnrmx=rnrm*TB(1.1)
rhist=Vector(TB){}
push!(rhist,rnrm)
#
# GMRES-IR loop
#
end

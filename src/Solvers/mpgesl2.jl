function mpgesl2(AF::MPFact, b)
mpdebug=false
normtype=Inf
TB = eltype(b)
TH = eltype(AF.AH)
(TH == TB) || error("inconsistent precisions")
(TH == Float64) ? tolf=1.e-13 : tolf=1.e-6
TFact = eltype(AF.AL)
TB = eltype(b)
AD=AF.AH
bnrm=norm(b,normtype)
bsc=b
AFS=AF.AF
bS=TFact.(bsc)
if (typeof(AF)==MPGTest)
x=zeros(size(b))
else
ldiv!(AFS,bS); x = TH.(bS)
end
#n=length(b); x=zeros(TH,n)
#
r = copy(x)
mul!(r,AD,x)
r .*= -1.0
axpy!(1.0, bsc, r)
tol = tolf * bnrm
rs=bS
rhist=Vector{Float64}()
rnrm=norm(r, normtype)
rnrmx = rnrm*1.1
itc=0
println("High precision = $TH, Factorization precision = $TFact")
#println("current rnrm($itc) = $rnrm")
push!(rhist,rnrm)
while (rnrm > tol) && (rnrm < rnrmx)
r ./=rnrm
#
# Hungry IR?
#
if (typeof(AF) == MPHTest)
ldiv!(AFS,r)
elseif (typeof(AF) == MPTest)
#
# Solver for the correction; AFS in factorization precision
#
rs .= TFact.(r)
ldiv!(AFS,rs)
r .= TH.(rs)
elseif (typeof(AF) == MPGTest)
x0=zeros(size(r))
eta=1.e-4
V=AF.VH
kout=kl_gmres(x0, r, MPhatv, V, eta, MPhptv; pdata=AF, side="left");
#println(kout.reshist,"  ",kout.idid)
r .= kout.sol
else
error("missing MP Fact type")
end
r .*= rnrm
#
#
#
x .+= r
mul!(r,AD,x)
r .*= -1.0
axpy!(1.0, bsc, r)
rnrmx = rnrm
rnrm=norm(r,normtype)
itc += 1
#println("current rnrm($itc) = $rnrm")
push!(rhist,rnrm)
mpdebug && println("Iteration $itc: rnorm = $rnrm, tol = $tol");
(rnrm >= rnrmx) && println("Norm increased")
end
println("Residual history = $rhist")
return x
end

#function mpgesl2(MPGF::xMPGTest, b)
#TB = eltype(b)
#(TB == Float64) ? eta=1.e-13 : eta=1.e-6
#n=length(b)
#x0=zeros(n)
#V=MPGF.VH
#kout=kl_gmres(x0, b, MPhatv, V, eta, MPhptv; pdata=MPGF, side="left");
#println(kout.reshist)
#return kout.sol
#end

function mpgesl2(AF::MPFact, b; reporting=false, verbose=true)
mpdebug=false
normtype=Inf
TB = eltype(b)
MPStats=getStats(AF)
TL=MPStats.TL
TH=MPStats.TH
TFact=MPStats.TFact
#TH = eltype(AF.AH); TL=getTL(AF); TFact = eltype(AF.AL); 
(TH == TB) || error("inconsistent precisions")
(TH == Float64) ? tolf=1.e-13 : tolf=1.e-6
MPStats=getStats(AF)
Meth=MPStats.Meth
if verbose
println(Meth, ": High precision = $TH, Low precision = $TL,
Factorization storage precision = $TFact")
end
AD=AF.AH
bnrm=norm(b,normtype)
bsc=b
AFS=AF.AF
bS=TFact.(bsc)
if (typeof(AF)==MPGTest)
x=zeros(size(b))
else
#ldiv!(AFS,bS); x = TH.(bS)
n=length(b); x=zeros(TH,n)
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
if verbose
itcc=itc+1
println("Krylov stats: Iteration $itcc ",  kout.reshist,"  ",kout.idid)
end
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
if verbose
println("Residual history = $rhist")
end
if reporting
return (rhist = rhist, sol = x, TH=TH, TL=TL, TFact=TFact)
else
return x
end
end

function getTL(AF)
if (typeof(AF) == MPTest)
TL = eltype(AF.AL)
elseif (typeof(AF) == MPHTest)
TL = eltype(AF.AS)
elseif (typeof(AF) == MPGTest)
TL = eltype(AF.AS)
else
TX=typeof(AF)
error("illegal MPTest type $TX")
end
return TL
end

function getMeth(AF)
if (typeof(AF) == MPTest)
Meth = "IR";
elseif (typeof(AF) == MPHTest)
Meth = "Hungry IR";
elseif (typeof(AF) == MPGTest)
Meth = "IRGM";
else
TX=typeof(AF)
error("illegal MPTest type $TX")
end
return Meth
end

function getStats(AF)
TH = eltype(AF.AH); TL=getTL(AF); TFact = eltype(AF.AL); 
if (typeof(AF) == MPGTest)
   MPStats=MPGStats()
else
   MPStats=MPIRStats(TH, TL, TFact)
end
return MPStats
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

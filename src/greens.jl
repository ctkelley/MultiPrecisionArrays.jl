function ktst(n=1000,kappa=1.0; hungry=false, verbose=false, 
              reporting=false)
G=Gmat(n);
AD=I - kappa*G;
cd = cond(AD);
h=1.0/(n+1.0);
X=collect(h:h:1.0-h);
xe=sin.(pi*X) .* exp.(X .- .5);
b = AD*xe;
jacobipc=false
if jacobipc
dd=1.0./diag(AD); D=Diagonal(dd); lmul!(D,b); lmul!(D,AD);
end
xd=AD\b;
println("Solution error = ",norm(xd .- xe,Inf)/norm(xe,Inf),
        " Residual norm = ", norm(b - AD*xd, Inf)/norm(b,Inf))
println("cond(A) = $cd")
AS=Float32.(AD);
bs=Float32.(b);
halfok=true
xs=AS\bs;
nsd=norm(xs-xd,Inf)/norm(xd,Inf);
println("Single precision Solution error = $nsd")
spacepad()
if halfok
AMP=MPGArray(AD, Float16);
AMF = mpglu!(AMP);
#xmp=AMF\b;
solout=\(AMF,b; reporting=reporting, verbose=verbose)
if reporting
xmp=solout.sol
else
xmp=solout
end
nshr= norm(b-AD*xmp,Inf)/norm(b,Inf)
nsm=norm(xmp-xd,Inf)/norm(xd,Inf); 
end
spacepad()
if hungry
AMD=MPHArray(AD)
AMFD=mphlu!(AMD)
else
AMD=MPArray(AD)
AMFD=mplu!(AMD)
end
#xmpd=AMFD\b;
solout=\(AMFD,b; reporting=reporting, verbose=verbose)
if reporting
xmpd=solout.sol
else
xmpd=solout
end
nsdr= norm(b-AD*xmpd,Inf)/norm(b,Inf)
nsdd = norm(xmpd-xd,Inf)/norm(xd,Inf)
spacepad()
if halfok
println("IR error = $nsdd, IRGM error = $nsm")
println("Residuals: IRGM: $nsdr, IR: $nsdr")
return [nsdd nsm]
else
println("outgoing relative residual norm = ",
          norm(b-AD*xmpd,Inf)/norm(b,Inf))
println("IR error = $nsdd")
return [nsdd nsd]
end
end

function dstest(n=1000,kappa=1.0;TL=Float32, TH=Float64)
#
# spectral radius of iteration matrix
#
G=Gmat(n);
AD=I - kappa*G;
cd = cond(AD);
AS=TL.(AD);
ASF=lu!(AS);
ADF=LU(TH.(AS), ASF.ipiv, ASF.info)
A1=ADF\AD;
c1=cond(A1)
M=I - A1;
sprad=opnorm(M,2)
println("condition = $cd, NormM = $sprad, modified condition = $c1")
end



function gtst(n)
G = Gmat(n);
z=rand(n); 
D2=Lap1d(n);
P=G*z;
Q=D2\z;
w=D2*P;
#println(norm(w-z,Inf),"  ",norm(P-Q,Inf))
end

function Gmat(n)
h=1.0/(n+1.0);
X=collect(h:h:1.0-h);
G=[greens(x,y) for x in X, y in X]
G.*=h;
return G
end


function ker(x,y)
ker = x*(x-y)^2
end

function greens(x,y)
if x > y
   gf=y*(1.0-x)
else
   gf=x*(1.0-y)
end
return gf
end

function spacepad()
println("  ")
println(" ------------------- ")
println("  ")
end

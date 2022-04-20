function ktst(n=1000,kappa=1.0; hungry=false)
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
#b=rand(n);
#b=sin.(X);
xd=AD\b;
println(norm(xd .- xe,Inf)/norm(xe,Inf),
        "   ", norm(b - AD*xd, Inf)/norm(b,Inf))
AS=Float32.(AD);
bs=Float32.(b);
halfok=true
xs=AS\bs;
nsd=norm(xs-xd,Inf)/norm(xd,Inf);
if halfok
println("Double-Half IRGM")
AMP=MPGArray(AD, Float16);
AMF = mpglu!(AMP);
xmp=AMF\b;
nsm=norm(xmp-xd,Inf)/norm(xd,Inf); 
end
if hungry
println("Double-Single Hungry")
AMD=MPHArray(AD)
AMFD=mphlu!(AMD)
else
println("Double-Single IR")
AMD=MPArray(AD)
AMFD=mplu!(AMD)
end
xmpd=AMFD\b;
nsdd=norm(xmpd-xd,Inf)/norm(xd,Inf)
if halfok
println("cond=$cd; s2derr=$nsd, mpderr=$nsdd, mpherr=$nsm")
return [cd nsdd nsd nsm]
else
println("outgoing relative residual norm = ",
          norm(b-AD*xmpd,Inf)/norm(b,Inf))
println("s2derr=$nsd, mpderr=$nsdd")
return [cd nsdd nsd]
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


"""
eir(n=32)

Make sure that Heavy IR and Expensive IR give the same results which
are not the same as normal IR.
"""
function eir(n=512, kappa=.1; TL=Float32, TH=Float64)
G = Gmat(n);
A = I + .1*TH.(G);
AC = copy(A);
AHC = copy(A);
b = ones(TH,n);
c = ones(TH,n);
d = ones(TH,n)
xe = A\b;
MPA=MPArray(A);
MPB=MPArray(AC);
MPHA = MPHArray(AHC)
MPF=mplu!(MPA);
MPFE=mplu!(MPB; expensive=true);
MPFH=mphlu!(MPHA)
normout=\(MPF,b; reporting=true);
Eout=\(MPFE,c; reporting=true);
hout=\(MPFH,d; reporting=true);
Ldiff = length(normout.rhist) - length(Eout.rhist)
Lengthok = (Ldiff != 0)
Hdiff = length(hout.rhist) - length(Eout.rhist)
if Hdiff == 0
Rdiff = norm(hout.rhist - Eout.rhist,Inf)
HEok = (Rdiff < 1.e-15)
else
println("HIR and EIR not the same. Residual histories have different size.")
return false
end
eirok=HEok && Lengthok
println(eirok)
return(Eout=Eout, normout=normout, hout=hout)
end



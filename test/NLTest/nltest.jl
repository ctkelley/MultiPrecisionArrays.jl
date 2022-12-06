"""
nltest()

Use MPArrays to solve H-equation
"""
function nltest()
hpass=nlhtest();
end

#
# H-equation test: Float32 Jacobian vs MP(32,16) MPArray
#
function nlhtest(n=32, c=.999)
FV=zeros(n);
JV=zeros(Float32,n,n);
JVD=zeros(Float64,n,n);
JV16=zeros(Float16,n,n);
JVMP=MPArray(JV);
JVMPD=MPArray(JVD);
x0=ones(n);
tol=1.e-14
dohalf=false
hdata = heqinit(x0, c);
nout=nsol(heqf!, x0, FV, JV, heqJ!;
          rtol=tol, atol=tol, pdata = hdata, sham = 1, jfact=lu!)
mpnout=nsol(heqf!, x0, FV, JVMP, jheqmp!;
          rtol=tol, atol=tol, pdata = hdata, sham = 1, jfact=mplu!)
llu=length(nout.history)
lmp=length(mpnout.history)
hpass=(llu==lmp)
end

function jheqmp!(JVMP, FV, x,hdata)
JV=JVMP.AH
JL=JVMP.AL
JV = heqJ!(JV, FV, x, hdata)
TH=eltype(JVMP)
TH == Float32 ? TL = Float16 : TL = Float32
JL .= TL.(JV)
return JVMP
end


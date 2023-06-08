"""
nltest()

Use MPArrays to solve H-equation
"""
function nltest()
hpass=nlhtest();
hpass16=nlhtest
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
JVMPH=MPHArray(JV);
JVMPD=MPArray(JVD);
x0=ones(n);
tol=1.e-14
hdata = heqinit(x0, c);
nout=nsol(heqf!, x0, FV, JV, heqJ!;
          rtol=tol, atol=tol, pdata = hdata, sham = 1, jfact=lu!)
mpnout=nsol(heqf!, x0, FV, JVMP, jheqmp!;
          rtol=tol, atol=tol, pdata = hdata, sham = 1, jfact=mplu!)
mphnout=nsol(heqf!, x0, FV, JVMPH, jheqmp!;
          rtol=tol, atol=tol, pdata = hdata, sham = 1, jfact=mpglu!)
llu=length(nout.history)
lmp=length(mpnout.history)
lmph=length(mphnout.history)
hpass=(llu==lmp) && (lmph == lmp)
end

function jheqmp!(JVMP, FV, x,hdata)
JV=JVMP.AH
JL=JVMP.AL
JV = heqJ!(JV, FV, x, hdata)
TH=eltype(JVMP)
TL=eltype(JL)
JL .= TL.(JV)
return JVMP
end

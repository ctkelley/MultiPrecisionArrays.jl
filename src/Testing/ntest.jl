struct Data_4_Plots
relreshist::Vector{Float64}
itdata::Vector{Int}
fdata::Vector{Int}
legend::String
end

function htest(n=1024,c=.99)
FV=zeros(n);
JV=zeros(Float32,n,n);
JVD=zeros(Float64,n,n);
JV16=zeros(Float16,n,n);
JVMP=MPArray(JV);
JVMPD=MPArray(JVD);
x0=ones(n);
tol=1.e-14
dohalf=false
#
# Set up the data for the plots. Look at the top of this file for the
# definition of the Data_4_Plots structure.
#
plot_hist = Vector{Data_4_Plots}()
#
hdata = heqinit(x0, c);
nout=nsol(heqf!, x0, FV, JV, heqJ!;
          rtol=tol, atol=tol, pdata = hdata, sham = 1, jfact=lu!)
nl_stats!(plot_hist, nout, "Single")
println("IR S-H")
mpnout=nsol(heqf!, x0, FV, JVMP, jheqmp!; 
          rtol=tol, atol=tol, pdata = hdata, sham = 1, jfact=mplu!)
nl_stats!(plot_hist, mpnout, "IR S-H")
println("IR D-S")
mpnoutd=nsol(heqf!, x0, FV, JVMPD, jheqmp!; 
          rtol=tol, atol=tol, pdata = hdata, sham = 1, jfact=mplu!)
nl_stats!(plot_hist, mpnoutd, "IR D-S")
if dohalf
nout16=nsol(heqf!, x0, FV, JV16, heqJ!; 
          rtol=1.e-10, atol=1.e-10, pdata = hdata, sham = 1, jfact=lu!)
nhist16=nout16.history/nout16.history[1]
nl_stats!(plot_hist, nout16, "Half")
end
caption=@sprintf("c = %5.4f; n = %5d", c, n)
plot_its_newton(plot_hist, caption)
if dohalf
return hout=(nout=nout, nout16=nout16, mpnoutd=mpnoutd, mpnout=mpnout)
else
return hout=(nout=nout, mpnoutd=mpnoutd, mpnout=mpnout)
end
end

function basicnewt(x)
FV=zeros(2)
JV=zeros(Float32,2,2);
JVMP=MPArray(JV)
nout=nsol(basic2d!, x, FV, JV, jbasic2d!; sham=1, jfact=lu!)
mpout=nsol(basic2d!, x, FV, JVMP, jacmp!; sham=1, jfact=mplu!)
return (nout=nout, mpout=mpout)
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

function jacmp!(JVMP, FV, x)
JV=JVMP.AH
JL=JVMP.AL
JV = jbasic2d!(JV, FV, x)
JL .= Float16.(JV)
return JVMP
end

function basic2d!(FV, x)
    FV[1] = x[1] * x[1] - 2.0 + x[2]
    FV[2] = exp(x[1] - 1) + x[2] * x[2] - 2.0
    return FV
end

function jbasic2d!(JV, FV, x)
    JV[1, 1] = 2 * x[1]
    JV[1, 2] = 1.0
    JV[2, 1] = exp(x[1] - 1)
    JV[2, 2] = 2 * x[2]
    return JV
end


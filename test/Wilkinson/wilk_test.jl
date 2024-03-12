function hi_precision_res()
a1 = wilk_test()
a8 = wilk_test(31, 800.0)
return a1 && a8
end


function wilk_test(n=31, a=1.0)
G=Gmat(n, Float32);
#
# Well conditioned case
#
alpha=Float32(a)
A=I + alpha*G;
b=ones(Float32,n);
AD=Float64.(A);
bd=Float64.(b);
xe=AD\bd;
# set up for high precision residual
AF = mplu(A; TF=Float32, onthefly=true);
mout = \(AF, b; TR=Float64, reporting=true);
lres=length(mout.rhist)
ndel=norm(mout.sol - xe,Inf)
l1ok=(lres == 4)
s1ok=(ndel < 1.e-10);
A1ok=l1ok && s1ok
A1ok || println("32-32 Failure for alpha=$alpha:   ", lres,"  ",ndel)
BF=mplu(A; onthefly=true);
mout = \(BF, b; TR=Float64, reporting=true);
lres=length(mout.rhist)
ndel=norm(mout.sol - xe,Inf)
(alpha < 500.0) ? lt=6 : lt=8
l2ok=(lres == lt) 
s2ok=(ndel < 1.e-10);
A2ok=l2ok && s1ok
A2ok || println("32-16 Failure for alpha=$alpha:   ", lres,"  ",ndel)
BF=mplu(A; onthefly=true);
return A2ok && s1ok
end

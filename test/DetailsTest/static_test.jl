function static_test(N=20; alpha=1.0)
AS = @MArray zeros(N,N);
bs = @MArray zeros(N);
AD = I - alpha*Gmat(N);
AS .= AD;
x=ones(N);
b=AD*x;
bs .= b;
xd=AD\b;
xs=AS\bs;
solok = (norm(xd-xs,Inf) < 1.e-14) && norm(xd-x,Inf) < 1.e-14
solok || println("Static sol problem")
MAF = mplu(AS)
MDF = mplu(AD)
xmd=MDF\b;
xms=MDF\bs;
mpsolok = (norm(xmd-xms,Inf) < 1.e-14) && norm(xmd-x,Inf) < 1.e-14
mpsolok || println("Static mpsol problem")
MAF2 = mplu(AS; TF=Float16)
MDF2 = mplu(AD; TF=Float16)
xmd2=MDF2\b;
xms2=MDF2\bs;
mpsolok2 = (norm(xmd2-xms2,Inf) < 1.e-14) && norm(xmd2-x,Inf) < 1.e-14
mpsolok2 || println("Static mpsol16 problem")
stok= mpsolok && solok && mpsolok2
end

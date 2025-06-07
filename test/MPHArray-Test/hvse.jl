#
# Things are getting strange with this test. It seems that
# ldiv!, which I use, does not give the same results as \\
# This is trouble if TF=Float16 and TW=Float32.
# This happens for both OpenBlas and AppleAccelerate.
#
# Everything is fine for Float32-Float64 and for regular mplu.
#
function hvse(N=128)
(MPF, MPHF, x, b) = unitt(N);
eout=mpkrir(MPF, b; reporting=true);
hout=mpkrir(MPHF, b; reporting=true);
leno=length(eout.rhist); heno=length(hout.rhist);
mlen=min(leno,heno)
nres=norm(eout.rhist[1:mlen] - hout.rhist[1:mlen], Inf);
sres=norm(eout.sol - hout.sol, Inf);
tolok=1.e-14
hvseok= ( (nres < tolok) && (sres < tolok) )
hvseok || println("MGF problem resdiff=$nres or soldiff=$sres too large")
return hvseok
end


function unitt(N; TW=Float32, TF=Float16)
A64=I - 800.0*Gmat(N);
A=TW.(A64);
MPA=MPGArray(A);
MPF=mpglu!(MPA);
MPH=MPHArray(A);
MPHF=mpglu!(MPH);
x=ones(TW,N)
b=A*x
return (MPF=MPF, MPHF=MPHF, x=x, b=b)
end


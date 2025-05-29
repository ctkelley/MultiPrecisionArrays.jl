#
# Things are getting strange with this test. I'm now seeing
# for lower triangular L in Float 16 and b in Float32 that
# L\b is not the same as Float32.(L)\b which seemed to be fine 
# before. This happens for both OpenBlas and AppleAccelerate.
#
# My suspicion is that the Float16 bit is handled with the Julia
# internal solver and that the Float32 part uses BLAS, so they are not
# the same operations. 
#
# Everything is fine for Float32-Float64.
#
function hvse(N=128)
(MPF, MPHF, x, b) = unitt(N);
eout=mpkrir(MPF, b; reporting=true);
hout=mpkrir(MPHF, b; reporting=true);
leno=length(eout.rhist); heno=length(hout.rhist);
mlen=min(leno,heno)
nres=norm(eout.rhist[1:mlen] - hout.rhist[1:mlen], Inf);
sres=norm(eout.sol - hout.sol, Inf);
tolok=eps(Float16)
hvseok= ( (nres < tolok) && (sres < tolok) )
hvseok || println("MGF problem resdiff=$nres or soldiff=$sres too large")
return hvseok
end


function unitt(N)
A64=I - 800.0*Gmat(N);
A=Float32.(A64);
MPA=MPGArray(A);
MPF=mpglu!(MPA);
MPH=MPHArray(A);
MPHF=mpglu!(MPH);
x=ones(Float32,N)
b=A*x
return (MPF=MPF, MPHF=MPHF, x=x, b=b)
end


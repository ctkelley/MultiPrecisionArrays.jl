function hvse(N=128)
(MPF, MPHF, x, b) = unitt(N)
eout=mpgmir(MPF, b; reporting=true)
hout=mpgmir(MPHF, b; reporting=true)
nres=norm(eout.rhist - hout.rhist, Inf)
sres=norm(eout.sol - hout.sol, Inf)
hvseok= ( (nres < 1.e-15) && (sres < 1.e-15) )
hvseok || println("MGF problem $nres or $sres too large")
return hvseok
end


function unitt(N)
A64=I - 800.0*Gmat(N)
A=Float32.(A64)
MPA=MPGArray(A)
MPF=mpglu!(MPA)
MPH=MPHArray(A)
MPHF=mpglu!(MPH)
x=ones(Float32,N)
b=A*x
return (MPF=MPF, MPHF=MPHF, x=x, b=b)
end


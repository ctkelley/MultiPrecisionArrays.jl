function mpgmtest(n=100)
normok=true
for TA in [Float64, Float32]
for alpha in [1.0, 10.0, 800.0]
    AD = TA.(I - alpha*Gmat(n))
    xs = ones(TA,n)
    b = AD*xs
    ADM = MPHArray(AD)
    ADF = mpglu!(ADM)
    z = mpgmir(ADF, b)
    delnorm=norm(z - xs, Inf)
    normok=normok && testnorm(delnorm, TA, alpha)
end
end
return normok
end

function testnorm(delnorm, TA, alpha)
if TA == Float64
   normok = (alpha < 750.0) ? (delnorm < 1.e-13) : (delnorm < 1.e-11)
else
   normok = (alpha < 750.0) ? (delnorm < 1.e-6) : (delnorm < 1.e-3)
end
return normok
end


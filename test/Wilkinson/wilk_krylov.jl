function wilk_krylov(n=4096,a=800.0; basissize=10)
G=Gmat(n, Float32);
#
# alpha=1.0 = well conditioned
# alpha=800.0 = moderately ill conditioned
#
alpha=Float32(a)
A=I + alpha*G;
b=ones(Float32,n);
AD=Float64.(A);
bd=Float64.(b);
# solve the promoted problem
xp = AD\bd;
# Set it up for IR-GMRES
AFG=mpglu(A; TR=Float64, basissize=basissize)
mgout=\(AFG, b; reporting=true);
lmg=length(mgout.rhist)
promdiffg = norm(xp-mgout.sol,Inf)
#
# I want to see somewhere between 7 and 9 iterations
# and an error < 10^{-13}
#
gmresok= (7 <= lmg <= 9) && (promdiffg < 1.e-13)
#println(lmg,"  ",promdiffg,"  ",gmresok)
gmresok || println("IRGMRES error. lmg=$lmg, promdiffg=$promdiffg")
ABG=mpblu(A; TR=Float64)
mbout=\(ABG, b; reporting=true);
#
# BiCGSTAB needs fewer iterations because I limit the dimension of the
# Krylov basis to 10 for GMRES
#
lmb=length(mbout.rhist)
promdiffb = norm(xp-mbout.sol,Inf)
bicgstabok= (4 <= lmb <= 6) && (promdiffb < 1.e-13)
bicgstabok || println("IRBiCGTAB error. lmb=$lmb, promdiffb=$promdiffb")
#println(lmb,"  ",promdiffb,"  ",bicgstabok)
#results=(mgout=mgout, mbout=mbout)
IRKrylovok = (gmresok && bicgstabok)
return IRKrylovok
end





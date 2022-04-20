function testgm(n=100, kappa=1.0, eta=1.e-8, dpre=false)
FT=Float16
G=rand(n,n)
G .-= .5
G .*= 1.0/n
A0 = I - kappa*G
x=ones(n); b0=A0*x; 
if dpre
D=Diagonal(A0)
A=D\A0
b=D\b0
else
A=A0
b=b0
end
MPA=MPHArray(A, FT);
MPHF=mphlu!(MPA);
x0=zeros(n)
V=zeros(n,20)
println(cond(A),"   ",cond(MPHF.AF\A))
kout=kl_gmres(x0, b, MPhatv, V, eta, MPhptv; pdata=MPHF, side="left");
return kout
end

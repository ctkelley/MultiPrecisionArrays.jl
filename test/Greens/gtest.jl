function greensok(n=31)
G=Gmat(n)
(e32,l32) = gtest(G);
ok32 = (e32 < 1.e-13) && (l32 == 3)
ok32 || println("Greens fail at TL=Float64-32")
(e16,l16) = gtest(G; TL=Float16);
ok16 = (e16 < 1.e-13) && (l16 == 6)
ok16 || println("Greens fail at TL=Float64-16")
(e3216,l3216) = gtest(G; TL=Float16,TH=Float32);
ok3216 = (e3216 < 1.e-6) && (l3216 == 4)
ok3216 || println("Greens fail at TL=Float32-16")
lightok=(ok16 && ok32 && ok3216)
return lightok
end

function greensHok(n=31)
G=Gmat(n)
(e32,l32) = gtestH(G);
ok32 = (e32 < 1.e-13) && (l32 == 3)
ok32 || println("GreensH fail at TL=Float64-32")
(e16,l16) = gtestH(G; TL=Float16);
ok16 = (e16 < 1.e-14) && (l16 == 6)
ok16 || println("GreensH fail at TL=Float64-16")
(e3216,l3216) = gtestH(G; TL=Float16,TH=Float32);
ok3216 = (e3216 < 1.e-6) && (l3216 == 3)
ok3216 || println("GreensH fail at TL=Float32-16")
heavyok = (ok16 && ok32 && ok3216)
return heavyok
end


function gtest(G; TL=Float32,TH=Float64)
#G=Gmat(n)
(n,n)=size(G)
A = I + TH.(G)
b=ones(TH,n)
xe = A\b
MPA=MPArray(A;TL=TL)
MPF=mplu!(MPA)
soldata=\(MPF,b;reporting=true)
xm=soldata.sol
rhist=soldata.rhist
nerr=norm(xm-xe,Inf)
lenit=length(rhist)
return (nerr, lenit)
end

function gtestH(G; TL=Float32,TH=Float64)
#G=Gmat(n)
(n,n)=size(G)
A = I + TH.(G)
b=ones(TH,n)
xe = A\b
MPA=MPHArray(A;TL=TL)
MPF=mphlu!(MPA)
soldata=\(MPF,b;reporting=true)
xm=soldata.sol
rhist=soldata.rhist
nerr=norm(xm-xe,Inf)
lenit=length(rhist)
return (nerr, lenit)
end

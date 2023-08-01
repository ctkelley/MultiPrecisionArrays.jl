"""
greensok(n=31)

Test vanilla IR with the inverse Laplacian. 
"""
function greensok(n=31)
G=Gmat(n)
(e32,l32) = gtest(G);
ok32 = (e32 < 1.e-13) && (l32 <= 4)
ok32 || println("Greens fail at TL=Float64-32, $e32, $l32")
(e16,l16) = gtest(G; TL=Float16);
ok16 = (e16 < 1.e-13) && (l16 == 6)
ok16 || println("Greens fail at TL=Float64-16")
(e3216,l3216) = gtest(G; TL=Float16,TH=Float32);
ok3216 = (e3216 < 1.e-6) && (l3216 == 4)
ok3216 || println("Greens fail at TL=Float32-16: $e3216, $l3216")
lightok=(ok16 && ok32 && ok3216)
return lightok
end

"""
greenHsok(n=31)

Test heavy IR with the inverse Laplacian. 
"""
function greensHok(n=31)
G=Gmat(n)
(e32,l32) = gtestH(G);
ok32 = (e32 < 1.e-13) && (l32 == 3)
ok32 || println("GreensH fail at TL=Float64-32, $ok32, $e32, $le32")
(e16,l16) = gtestH(G; TL=Float16);
ok16 = (e16 < 1.e-14) && (l16 == 6)
ok16 || println("GreensH fail at TL=Float64-16")
(e3216,l3216) = gtestH(G; TL=Float16,TH=Float32);
ok3216 = (e3216 < 1.e-6) && (l3216 <= 4)
ok3216 || println("GreensH fail at TL=Float32-16")
heavyok = (ok16 && ok32 && ok3216)
return heavyok
end

"""
greensEvsH(n=31)

Make sure heavy IR and expensive IR are give identical results
"""
function greensEvsH(n=31)
G=Gmat(n)
(ee32, le32, he32) =gtestE(G; reshistout=true)
(e32, l32, h32) = gtestH(G; reshistout=true)
edel=abs(e32-ee32)/e32
hdel=norm(h32-he32,Inf)/norm(he32,Inf)
#errok=(e32==ee32)
#histok=(h32==he32)
errok=(edel < 1.e-2)
histok=(hdel < 1.e-15)
EvsHok32=(errok && histok)
EvsHok32 || println("TL=F32: heavy IR and expensive IR differ")
(ee16, le16, he16) =gtestE(G; TL=Float16, reshistout=true)
(e16, l16, h16) = gtestH(G; TL=Float16, reshistout=true)
errok=(e16==ee16)
histok=(h16==he16)
histok || println("he16 - h17 error  ", norm(he16-h16,Inf))
EvsHok16=(errok && histok)
(ee3216, le3216, he3216) =gtestE(G; TL=Float16, TH=Float32, reshistout=true)
(e3216, l3216, h3216) = gtestH(G; TL=Float16, TH=Float32, reshistout=true)
errok=(e3216==ee3216)
histok=(h3216==he3216)
histok=(norm(he3216-he3216,Inf) < 1.e-14)
EvsHok3216=(errok && histok)
EvsHPass=EvsHok32 && EvsHok16 && EvsHok3216
EvsHPass || println("heavy vs expensive too far apart")
return EvsHPass
end

"""
gtest(G; TL=Float32,TH=Float64)

Solve the inverse Laplacian problem with IR
"""
function gtest(G; TL=Float32,TH=Float64)
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

"""
gtestH(G; TL=Float32,TH=Float64, reshistout=false)

Solve the inverse Laplacian problem with heavy IR
"""
function gtestH(G; TL=Float32,TH=Float64, reshistout=false)
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
if reshistout
return (nerr, lenit, rhist)
else
return (nerr, lenit)
end
end

"""
gtestE(G; TL=Float32,TH=Float64)

Solve the inverse Laplacian problem with expensive IR
This is only for CI
"""
function gtestE(G; TL=Float32,TH=Float64, reshistout=false)
#G=Gmat(n)
(n,n)=size(G)
A = I + TH.(G)
b=ones(TH,n)
xe = A\b
MPA=MPEArray(A;TL=TL)
MPF=mplu!(MPA)
soldata=\(MPF,b;reporting=true)
xm=soldata.sol
rhist=soldata.rhist
nerr=norm(xm-xe,Inf)
lenit=length(rhist)
if reshistout
return (nerr, lenit, rhist)
else
return (nerr, lenit)
end
end



function vtimes(p=4)
pv=collect(1:1:p);
ppv=2 .^ pv;
nv=512 .* ppv;
Time_LU=zeros(p,3)
for p in pv
n=nv[p]
AD=rand(n,n)
AS=Float32.(AD)
AH=Float16.(AD)
td=@belapsed(lu($AD))
ts=@belapsed(lu($AS))
th=@belapsed(lu($AH))
Time_LU[p,:]=[td,ts,th]
end
return Time_LU
end


function HalfTime(p=3)
pv = collect(1:1:p)
ppv = 2 .^ pv
nv = 512 .* ppv
Time_LU = zeros(p, 5)
headers=["N", "LU64", "LU32", "LU16", "LU16/LU64"]
for p in pv
N=nv[p]
G=Gmat(N);
AD=I - G;
AS=Float32.(AD)
AH=Float16.(AD)
td = @belapsed lu($AD)
ts = @belapsed lu($AS)
th = @belapsed hlu($AH)
rt = th / td
Time_LU[p, :] = [N, td, ts, th, rt]
end
dformat = "%9d %9.2e %9.2e %9.2e %9.2e \n"
printf(fmt::String, args...) = @eval @printf($fmt, $(args...))
for i = 1:5
printf("%11s", headers[i])
end
printf("\n")
for p in pv
printf(dformat,Time_LU[p,:]...)
end
return Time_LU
end






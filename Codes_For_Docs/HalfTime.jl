function HalfTime(p=3)
pv = collect(1:1:p)
ppv = 2 .^ pv
nv = 512 .* ppv
Time_LU = zeros(p, 5)
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
return Time_LU
end






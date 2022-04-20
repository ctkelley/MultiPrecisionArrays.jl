function pitch(n)
AH=rand(n,n); AS=Float32.(AH); AS2=copy(AS);
bh=rand(n); bs=Float32.(bh); 
AHF=lu(AH); ASF=lu(AS);
tluh=@belapsed lu!($AH);
println("High precsion LU time = $tluh")
tlus=@belapsed lu!($AS);
println("Low precision LU time = $tlus")
thsol=@belapsed $AHF\$bh;
println("High precision solve time = $thsol")
tssol=@belapsed $ASF\$bh;
println("Low precision solve time = $tssol")
tgssol=@belapsed $ASF\$bs;
println("Better Low precision solve time = $tgssol")
rvec = [tluh tlus thsol tssol tgssol]
return rvec
end

function makemptab(m=4, texok=false)
AT=zeros(m,6)
p=collect(0:1:m-1)
tp = 2 .^ p
np = 32 .* tp
for idim=1:m
n=np[idim]
AT[idim,2:6] = pitch(n);
AT[idim,1]=n;
end
@printf("%5s     %5s      %5s       %5s        %5s        %5s \n",
         "n", "luh", "lus", "A\\b", "AS\\b", "AS\\bs")
for idim=1:m
@printf("%5d    %5.2e    %5.2e    %5.2e    %5.2e    %5.2e \n", 
       AT[idim,1], AT[idim,2], AT[idim,3], AT[idim,4], AT[idim,5],
       AT[idim,6])
end
headers=["N", "\$lu(\\ma_h)\$", "\$lu(\\ma_s)\$", 
"\$\\ma_h \\backslash \\vb",
"\$\\ma_s \\backslash \\vb",
"\$\\ma_s \\backslash \\vb_s"];
formats=("%5d  &  %5.2e &   %5.2e  &  %5.2e &  %5.2e  & %5.2e");
if texok
fprintTeX(headers, formats, AT)
end
return AT
end

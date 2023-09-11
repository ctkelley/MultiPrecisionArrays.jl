function slashtest()
G=Gmat(10)
A=I+G
x=ones(10)
b=A*x
MPA=MPArray(A)
MPB=deepcopy(MPA)
MPF=mplu!(MPA)
out1=\(MPF, b; reporting=true)
out2=\(MPB, b; reporting=true)
doubleok=(out1==out2)
doubleok || println("64-32 fails in slashtest")
As=Float32.(A)
xs=Float32.(x)
bs=As*xs;
MPAS=MPArray(As)
MPBS=deepcopy(MPAS)
MPFS=mplu!(MPAS)
out3=\(MPFS, bs; reporting=true)
out4=\(MPBS, bs; reporting=true)
singleok=(out3==out4)
singleok || println("32-16 fails in slashtest")
slashok=(doubleok && singleok)
return slashok
end


function conv_tab(plim=4)
c_data=zeros(plim,6)
for p=1:plim
    N=512*2^p
    c_data[p,:].=Convergence(N)
end
conv_hist(c_data)
end

 
function Convergence(N)
G=Gmat(N)
A = I - 800.0*G
A32 = Float32.(A)
A16 = Float16.(A)
AF=lu(A)
A32F=lu(A32)
A16F=lu(A16)
P = AF\A
P32 = A32F\A
P16 = A16F\A
Psh = A16F\A32
Nd = opnorm(I - P,1)
Ndsh = opnorm(I - Psh,1)
Nd32 = opnorm(I - P32,1)
Nd16 = opnorm(I - P16,1)
condA=cond(A,1)
dreport=[N, Nd, Nd32, Nd16, Ndsh, condA]
return dreport
#println("||M|| for A = I - 800*G; N = $N; cond(A,1)= $condA")
#println(" 64-64: ", Nd,
#        " |  64-32: ", Nd32,"  |  64-16: ", Nd16, " |  32-16: ", Ndsh)
#hformat="%7s %30s \n"
#printf(fmt::String, args...) = @eval @printf($fmt, $(args...))
#headers=["   ", "opnorm(M,1)"]
#dlab=["N/TH-TL","64-64","64-32","64-16","32-16"]
#labformat=" %6s %9s %9s %9s %9s \n"
#dformat="%5d  %13.2e %9.2e %9.2e %9.2e \n"
#printf(hformat, headers...)
#println("______________________________________________________________")
#printf(labformat,dlab...)
#printf("%7s","  ")
#printf(dformat,dreport...)
end

function conv_hist(c_data)
(m,n)=size(c_data)
for ir=1:m
println(c_data[ir,:])
end
end

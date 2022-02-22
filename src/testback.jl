function testback(n=100, Tlow=Float16)
#B=rand(n,n); B ./= opnorm(B,1); A = I + B'*B;
B=rand(n,n); B .-= .5*rand(n,n); B ./= n; A = I + 60.0*B;
AH=Tlow.(A);
AHF=lu(AH);
blow=norm(AHF.P*AHF.L*AHF.U - AH,Inf);
err16=(Float64.(AH) - A)
bhigh=norm(err16,Inf)
return (blow, bhigh, cond(A,1), opnorm(err16,1))
end


function Quality(m = 4, alpha=1.0, texok = false; T = Float64)
    p = collect(0:1:m-1)
    tp = 2 .^ p
    np = 512 .* tp
    itcount=zeros(Int64,m,2)
    IRData=zeros(m,7)
for idim=1:m
    N=np[idim]
    G=Gmat(N); 
    A = T.(I + alpha*G);
    xe=ones(T,N); b=A*xe;
    AC=copy(A)
    MPA=MPArray(A);
    MPAF=mplu!(MPA);
    MPE=MPArray(A; onthefly=true);
    MPEF=mplu!(MPE);
#
# Timings
#
    MPB=MPArray(A);
    MPBF=mplu!(MPB)
    MPEB=MPArray(A; onthefly=true);
    MPEBF=mplu!(MPEB)
#    TLP=@belapsed MPC\$b setup=(MPC=deepcopy($MPB)) evals =1
#    TMP=@belapsed MPEC\$b setup=(MPEC=deepcopy($MPEB)) evals=1
    TLP=@belapsed MPC\$b setup=(MPC=deepcopy($MPBF)) evals =1
    TMP=@belapsed MPEC\$b setup=(MPEC=deepcopy($MPEBF)) evals=1
#
    olps=\(MPAF,b; reporting=true)
    xlps=olps.sol    
    omps=\(MPEF,b; reporting=true)
    xmps=omps.sol
    llps=length(olps.rhist)
    lmps=length(omps.rhist)
    itcount[idim,:]=[llps, lmps]
#    println(llps,"  ",lmps)
    nerra=norm(xe-xlps, Inf)/norm(xe,Inf)
    nerre=norm(xe-xmps, Inf)/norm(xe,Inf)
    nresa=norm(b-A*xlps,Inf)/norm(b,Inf)
    nrese=norm(b-A*xmps,Inf)/norm(b,Inf)
    IRData[idim,:]=[N nerra nerre nresa nrese TLP TMP]
#    println("errLP = $nerra,  errMP=$nerre, ResLP=$nresa, ResMP=$nrese")
end
if texok
headers = ["N", "ELP", "EMP", "RLP", "RMP", "TLP", "TMP"]
formats = ("%5d  &  %5.1e &   %5.1e  &  %5.1e &  %5.1e  & %5.1e & %5.1e")
fprintTeX(headers,formats,IRData)
else
println(itcount)
@printf("%s  %s   %s   %s  %s   %s   %s \n", "    N", "    ELP", 
             "     EMP", "     RLP", "       RMP","     TLP", "    TMP")
for idim=1:m
    @printf("%5d    %5.1e    %5.1e    %5.1e    %5.1e    %5.1e   %5.1e \n", IRData[idim,:]...)
end
end
return
end
      



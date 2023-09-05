function Quality(m = 4, alpha=1.0, texok = false; T = Float64)
    p = collect(0:1:m-1)
    tp = 2 .^ p
    np = 512 .* tp
    @printf("%s  %s   %s   %s  %s \n", "    N", "    ELP", 
             "     EMP", "     RLP", "       RMP")
for idim=1:m
    N=np[idim]
    G=Gmat(N); 
    A = T.(I + alpha*G);
    xe=ones(T,N); b=A*xe;
    MPA=MPArray(A);
    MPAF=mplu!(MPA)
    MPE=MPArray(A; onthefly=true);
    MPEF=mplu!(MPE)
    xlps=MPAF\b
    xmps=MPEF\b
    nerra=norm(xe-xlps, Inf)/norm(xe,Inf)
    nerre=norm(xe-xmps, Inf)/norm(xe,Inf)
    nresa=norm(b-A*xlps,Inf)/norm(b,Inf)
    nrese=norm(b-A*xmps,Inf)/norm(b,Inf)
    @printf("%5d    %5.1e    %5.1e    %5.1e    %5.1e \n", 
             N, nerra, nerre, nresa, nrese)
#    println("errLP = $nerra,  errMP=$nerre, ResLP=$nresa, ResMP=$nrese")
end
return
end
      



function low2basis(AL,AH)
TL=eltype(AL)
TH=eltype(AH)
B=reinterpret(TH,AL)
(ml,nl)=size(AL)
(ml == nl) || error("AL not square")
(mk,nk)=size(B)
(nk == nl) || error("Dimension problem in low2basis")
BK=reshape(B,nl,mk)
BK .*= 0.0
return BK
end


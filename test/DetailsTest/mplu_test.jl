"""
mplu_test()

Make sure that mplu and mplu! do what they are supposed to do
"""
function mplu_test()
AD=rand(10,10); MPD=MPArray(AD); MPF1=mplu!(MPD); MPF2=mplu(AD);
eq64=test_eq(MPF1,MPF2)
ADx=rand(10,10); MPDx=MPArray(ADx; TL=Float16); 
MPF1x=mplu!(MPDx); MPF2x=mplu(ADx; TL=Float16);
eq64x=test_eq(MPF1x,MPF2x)
AS=Float32.(AD); MPS=MPArray(AS); MSF1=mplu!(MPS); MSF2=mplu(AS)
eq32=test_eq(MSF1,MSF2)
mpluok = (eq64 && eq64x && eq32)
mpluok || println("mplu failure")
return mpluok
end


function test_eq(MF1,MF2)
eqok=true
for nf in fieldnames(MPLFact)
gx=getfield(MF1,nf); hx =getfield(MF2,nf)
eqok= ((gx==hx) && eqok)
end
return eqok
end

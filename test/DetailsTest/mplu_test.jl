"""
mplu_test()

Make sure that mplu and mplu! do what they are supposed to do
"""
function mplu_test()
AD=rand(10,10); MPD=MPArray(AD); MPF1=mplu!(MPD); MPF2=mplu(AD);
eq64=test_eq(MPF1,MPF2)
eq64 || println("mplu t1 fails")
ADx=rand(10,10); MPDx=MPArray(ADx; TL=Float16); 
MPF1x=mplu!(MPDx); MPF2x=mplu(ADx; TL=Float16);
eq64x=test_eq(MPF1x,MPF2x)
eq64x || println("mplu t2 fails")
AS=Float32.(AD); MPS=MPArray(AS); MSF1=mplu!(MPS); MSF2=mplu(AS)
eq32=test_eq(MSF1,MSF2)
eq32 || println("mplu t3 fails")
mpluok = (eq64 && eq64x && eq32)
mpluok || println("mplu failure")
return mpluok
end


"""
mpglu_test()

Make sure that mpglu and mpglu! do what they are supposed to do
"""
function mpglu_test()
AD=rand(10,10); MPD=MPGArray(AD); MPF1=mpglu!(MPD); MPF2=mpglu(AD);
eq64=test_eq(MPF1,MPF2)
eq64 || println("mpglu t1 fails")
ADx=rand(10,10); MPDx=MPGArray(ADx; TL=Float16);
MPF1x=mpglu!(MPDx); MPF2x=mpglu(ADx; TL=Float16);
eq64x=test_eq(MPF1x,MPF2x)
eq64x || println("mpglu t2 fails")
AS=Float32.(AD); MPS=MPGArray(AS); MSF1=mpglu!(MPS); MSF2=mpglu(AS)
eq32=test_eq(MSF1,MSF2)
eq32 || println("mpglu t3 fails")
mpgluok = (eq64 && eq64x && eq32)
mpgluok || println("mpglu failure")
return mpgluok
end


function test_eq(MF1,MF2)
eqok=true
for nf in fieldnames(MPLFact)
gx=getfield(MF1,nf); hx =getfield(MF2,nf)
eqok= ((gx==hx) && eqok)
eqok || println(nf)
end
return eqok
end

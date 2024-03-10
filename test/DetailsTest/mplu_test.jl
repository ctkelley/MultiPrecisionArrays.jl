"""
mplu_test()

Make sure that mplu and mplu! do what they are supposed to do
"""
function mplu_test()
AD=rand(10,10); MPD=MPArray(AD); MPF1=mplu!(MPD); MPF2=mplu(AD);
eq64=test_eq(MPF1,MPF2)
eq64 || println("mplu t1 fails")
ADx=rand(10,10); MPDx=MPArray(ADx; TF=Float16); 
MPF1x=mplu!(MPDx); MPF2x=mplu(ADx; TF=Float16);
eq64x=test_eq(MPF1x,MPF2x)
eq64x || println("mplu t2 fails")
AS=Float32.(AD); MPS=MPArray(AS); MSF1=mplu!(MPS); MSF2=mplu(AS)
eq32=test_eq(MSF1,MSF2)
eq32 || println("mplu t3 fails")
BD=I + AD; MPB=mplu(BD); MPA=mplu(AD); 
MPX=mplu!(MPA,BD);
equp1=test_eq(MPX,MPB)
equp1 || println("mplu! failure, updating factorization")
mpluok = (eq64 && eq64x && eq32 && equp1)
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
ADx=rand(10,10); MPDx=MPGArray(ADx; TF=Float16);
MPF1x=mpglu!(MPDx); MPF2x=mpglu(ADx; TF=Float16);
eq64x=test_eq(MPF1x,MPF2x)
eq64x || println("mpglu t2 fails")
AS=Float32.(AD); MPS=MPGArray(AS); MSF1=mpglu!(MPS); MSF2=mpglu(AS)
eq32=test_eq(MSF1,MSF2)
eq32 || println("mpglu t3 fails")
BDx=I + ADx; MPBx=mpglu(BDx; TF=Float16);
MPTx=mpglu!(MPBx,ADx); 
equp2=test_eq(MPTx,MPF1x)
equp2 || println("mpglu! failure, updating factorization")
mpgluok = (eq64 && eq64x && eq32 && equp2)
mpgluok || println("mpglu failure")
return mpgluok
end

"""
mpblu_test()

Make sure that mpblu and mpblu! do what they are supposed to do
"""
function mpblu_test()
AD=rand(10,10); MPD=MPBArray(AD); MPF1=mpblu!(MPD); MPF2=mpblu(AD);
eq64=test_eq(MPF1,MPF2)
eq64 || println("mpblu t1 fails")
ADx=rand(10,10); MPDx=MPBArray(ADx; TF=Float16);
MPF1x=mpblu!(MPDx); MPF2x=mpblu(ADx; TF=Float16);
eq64x=test_eq(MPF1x,MPF2x)
eq64x || println("mpblu t2 fails")
AS=Float32.(AD); MPS=MPBArray(AS); MSF1=mpblu!(MPS); MSF2=mpblu(AS)
eq32=test_eq(MSF1,MSF2)
eq32 || println("mpblu t3 fails")
BDx=I + ADx; MPBx=mpblu(BDx; TF=Float16);
MPTx=mpblu!(MPBx,ADx);
equp2=test_eq(MPTx,MPF1x)
equp2 || println("mpblu! failure, updating factorization")
mpbluok = (eq64 && eq64x && eq32 && equp2)
mpbluok || println("mpblu failure")
return mpbluok
end



function test_eq(MF1,MF2)
eqok=true
for nf in fieldnames(MPLFact)
gx=getfield(MF1,nf); hx =getfield(MF2,nf)
eqok= ((gx==hx) && eqok)
(gx==hx) || println(nf)
end
return eqok
end



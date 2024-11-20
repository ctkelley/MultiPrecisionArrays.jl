function term_test(N=100,TF=Float32)
G=Gmat(N)
A=I - 800.0*G
b=ones(N)
# mplu test
AF=mplu(A; TF=TF, onthefly=true, residterm=true)
BF=mplu(A; TF=TF, onthefly=true, residterm=false)
mpluok = lentest(AF, BF, b, "mplu")
# mpglu test
AF = mpglu(A; TF=TF, residterm=true)
BF = mpglu(A; TF=TF, residterm=false)
mpgluok = lentest(AF, BF, b, "mpglu")
# mpblu test
AF = mpblu(A; TF=TF, residterm=true)
BF = mpblu(A; TF=TF, residterm=false)
mpbluok = lentest(AF, BF, b, "mpblu")
termtestok = (mpluok & mpgluok & mpbluok)
return termtestok
end

function lentest(AF, BF, b, fname)
resout=\(AF, b; reporting=true)
bresout=\(BF, b; reporting=true)
rlen=length(resout.rhist)
brlen=length(bresout.rhist)
termtestok=(rlen > brlen)
termtestok || println("term_test failure in $fname")
#println(termtestok)
#return (rhist=resout.rhist, brhist=bresout.rhist)
return termtestok
end



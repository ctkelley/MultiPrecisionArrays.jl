function term_test(N=100,TF=Float32)
G=Gmat(N)
A=I - 800.0*G
b=ones(N)
AF=mplu(A; TF=TF, onthefly=true, residterm=true)
BF=mplu(A; TF=TF, onthefly=true, residterm=false)
resout=\(AF, b; reporting=true)
bresout=\(BF, b; reporting=true)
rlen=length(resout.rhist)
brlen=length(bresout.rhist)
termtestok=(rlen > brlen)
termtestok || println("term_test failure")
#println(termtestok)
#return (rhist=resout.rhist, brhist=bresout.rhist)
return termtestok
end

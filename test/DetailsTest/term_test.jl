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
resout=\(AF, b; reporting=true);
bresout=\(BF, b; reporting=true);
rlen=length(resout.rhist)
brlen=length(bresout.rhist)
mlen=min(rlen,brlen)
termtestok=(norm(resout.rhist[1:mlen]-bresout.rhist[1:mlen]) < 1.e-10)
#termtestok=(rlen > brlen)
termtestok || println("term_test failure in $fname")
#println(termtestok)
#return (rhist=resout.rhist, brhist=bresout.rhist)
return termtestok
end

function test_term_parms(N=100, TF=Float32)
N=100; G=Gmat(N); A=I - 9.0*G;
AF=mplu(A); b=ones(N);
restore_default_parms()
plain_out = \(AF, b; reporting=true);
rhist1=plain_out.rhist;
lhist1=length(rhist1)
update_parms(; Cr=400.0);
swap_out = \(AF, b; reporting=true);
rhist2=swap_out.rhist;
lhist2=length(rhist2)
swapok = (lhist1 > lhist2)
swapok || println("Error in test_term_parms", lhist1,"  ",lhist2)
lmsx_out = update_parms(; litmax=2)
new_out=\(AF, b; reporting=true);
lhist3=length(new_out.rhist)
lhistok=(lhist3 < lhist1)
restore_default_parms()
return swapok && lhistok
end

function swap_test()
update_parms(; Cr=100.0, Ce=10.0, Rmax=.01)
swap1ok = terms_the_same(term_parms_default, term_parms)
restore_default_parms()
swap2ok = terms_the_same(term_parms_default, term_parms)
swap_working = swap2ok && ~swap1ok
return swap_working
end

function terms_the_same(term1::TERM, term2::TERM)
CrOK = (term1.Cr == term2.Cr)
CeOK = (term1.Ce == term2.Ce)
RmaxOK = (term1.Rmax == term2.Rmax)
sameOK = (CrOK && CeOK && RmaxOK)
return sameOK
end


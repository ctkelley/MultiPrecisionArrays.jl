function AbsArray()
A=rand(10,10); x=ones(10); b=A'*x;
B=A';
MPB=mplu(B)
z=MPB\b
absok = (norm(x - z,Inf) < 1.e-11)
absok || println(norm(x-z,Inf))
return absok
end

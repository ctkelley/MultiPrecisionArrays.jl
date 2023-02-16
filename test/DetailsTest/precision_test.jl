"""
precision_test()

Check for promotions to Float64 that should not be there.
"""
function precision_test()
#
# Get organized
#
A=rand(10,10); A32=Float32.(A);
b=rand(10); b32=Float32.(b);
xe=ones(10); xe32=ones(Float32,10); 
b .= A*xe; b32 .= A32*xe32;
#
# Cycle through all combinations and make sure the precisions line up.
#
MPA=MPArray(A); MPF=mplu!(MPA);
x=MPF\b;
ok64 = (eltype(x) == eltype(xe))
#
MPHA=MPHArray(A); MPHF=mphlu!(MPHA);
xh=MPHF\b;
ok64h = (eltype(xh) == eltype(xe))
doubleok=ok64 && ok64h
#
MPAS=MPArray(A32); MPFS=mplu!(MPAS);
xs=MPFS\b32; 
ok32 = (eltype(xs) == eltype(xe32))
#
MPHAS=MPHArray(A32); MPFHS=mphlu!(MPHAS);
xsh=MPFHS\b32; 
ok32h = (eltype(xsh) == eltype(xe32))
singleok = (ok32 && ok32h)
#
MPAx=MPArray(A;TL=Float16); MPFx=mplu!(MPAx);
x16=MPFx\b;
ok6416 = (eltype(x16) == eltype(xe))
#
MPHAx=MPHArray(A;TL=Float16); MPHFx=mphlu!(MPHAx);
xh16=MPHFx\b;
okh6416 = (eltype(xh16) == eltype(xe))
doubleok = (doubleok && ok6416 && okh6416)
#
MPE=MPEArray(A); MPEF=mplu!(MPE);
xee=MPEF\b;
oke= (eltype(xee) == eltype(xe))
#
MPEx=MPEArray(A; TL=Float16); MPEFx=mplu!(MPEx);
xeex=MPEFx\b;
okex= (eltype(xeex) == eltype(xe))
doubleok = (doubleok && oke && okex)
#
MPE32=MPEArray(A32); MPEF32=mplu!(MPE32);
xee32=MPEF32\b32;
oke32 = (eltype(xee32) == eltype(xe32))
singleok = (singleok && oke32)

precisionok= (singleok && doubleok)
end


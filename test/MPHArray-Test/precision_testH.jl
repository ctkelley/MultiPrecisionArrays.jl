"""
precision_test()

Check for promotions to Float64 that should not be there.
"""
function precision_testH()
    #
    # Get organized
    #
    A = rand(10, 10)
    A32 = Float32.(A)
    b = rand(10)
    b32 = Float32.(b)
    xe = ones(10)
    xe32 = ones(Float32, 10)
    b .= A * xe
    b32 .= A32 * xe32
    #
    # Cycle through all combinations and make sure the precisions line up.
    #
    MPHA = MPHArray(A)
    MPHF = mphlu!(MPHA)
    xh = MPHF \ b
    ok64h = (eltype(xh) == eltype(xe))
    doubleok = ok64h
    #
    #
    MPHAS = MPHArray(A32)
    MPFHS = mphlu!(MPHAS)
    xsh = MPFHS \ b32
    ok32h = (eltype(xsh) == eltype(xe32))
    singleok = ok32h
    #
    #
    MPHAx = MPHArray(A; TF = Float16)
    MPHFx = mphlu!(MPHAx)
    xh16 = MPHFx \ b
    okh6416 = (eltype(xh16) == eltype(xe))
    doubleok = (doubleok && okh6416)
    #
    #
    precisionokH = (singleok && doubleok)
end

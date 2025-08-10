"""
nltest()

Use MPArrays to solve H-equation
"""
function nltest()
    hpass = nlhtest()
    hpass16 = nlhtest16()
    return hpass && hpass16
end

#
# H-equation test: Float32 Jacobian vs MP(32,16) MPArray
#
function nlhtest(n = 32, c = 0.999)
    FV = zeros(n)
    JV = zeros(Float32, n, n)
    JVD = zeros(Float64, n, n)
    JV16 = zeros(Float16, n, n)
    JVMP = MPArray(JV)
    JVMPH = MPHArray(JV)
    JVMPG = MPGArray(JVD)
    x0 = ones(n)
    tol = 1.e-14
    hdata = heqinit(x0, c)
    test_data=(FV = FV, x0 = x0, tol = tol, hdata = hdata)
    # Normal Newton
    nout = nsol_test(JV, heqJ!, lu!, test_data)
    # MPLU Newton
    mpnout = nsol_test(JVMP, jheqmp!, mplu!, test_data)
    # IR-G Heavy Newton
    mphnout = nsol_test(JVMPH, jheqmp!, mpglu!, test_data)
    # IR-G Newton
    mpgnout = nsol_test(JVMPG, jheqmp!, mpglu!, test_data)
    #
    llu = length(nout.history)
    lmp = length(mpnout.history)
    lmph = length(mphnout.history)
    lmpg = length(mpgnout.history)
    hpass = (llu == lmp) && (lmph == lmp) && (lmp == lmpg)
end

#
# H-equation test: Float16 Jacobian vs IR-GMRES for singular Jacobian
#
function nlhtest16(n = 32, c = 1.0)
    FV = zeros(n)
    JV = zeros(Float32, n, n)
    JV16 = zeros(Float16, n, n)
    JVMP = MPArray(JV; onthefly = false)
    JVME = MPArray(JV; onthefly = true)
    JVMPG = MPGArray(JV)
    JVMPH = MPHArray(JV)
    x0 = ones(n)
    tol = 1.e-8
    hdata = heqinit(x0, c)
    test_data=(FV = FV, x0 = x0, tol = tol, hdata = hdata)
    # Normal Newton
    nout = nsol_test(JV, heqJ!, lu!, test_data)
    # Normal Newton F16 Jacobian
    nout16 = nsol_test(JV16, heqJ!, lu!, test_data)
    # MP Newton
    mpnout = nsol_test(JVMP, jheqmp!, mplu!, test_data)
    # MPH Newton
    mphnout = nsol_test(JVMPH, jheqmp!, mpglu!, test_data)
    # MPG Newton
    mpgmeout = nsol_test(JVMPG, jheqmp!, mpglu!, test_data)
    #
    llu = length(nout.history)
    ll16 = length(nout16.history)
    lmp = length(mpnout.history)
    lmph = length(mphnout.history)
    lmpge = length(mpgmeout.history)
    halfpass = (ll16 == 21) && (llu == lmph) && (lmp > llu) && (lmpge == llu)
    return halfpass
end

function nsol_test(JStore, jeval!, lsol!, test_data)
    FV=test_data.FV
    x0=test_data.x0
    tol=test_data.tol
    hdata=test_data.hdata
    test_out = nsol(
        heqf!,
        x0,
        FV,
        JStore,
        jeval!;
        rtol = tol,
        atol = tol,
        pdata = hdata,
        sham = 1,
        jfact = lsol!,
    )
    return test_out
end




function jheqmp!(JVMP, FV, x, hdata)
    JV = JVMP.AH
    JL = JVMP.AL
    JV = heqJ!(JV, FV, x, hdata)
    TW = eltype(JVMP)
    TF = eltype(JL)
    JL .= TF.(JV)
    return JVMP
end

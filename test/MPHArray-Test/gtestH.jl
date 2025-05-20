"""
greenHsok(n=31)

Test heavy IR with the inverse Laplacian. 
"""
function greensHok(n = 31)
    G = Gmat(n)
    (e32, l32) = gtestH(G)
    ok32 = (e32 < 1.e-13) && (l32 <= 5)
    ok32 || println("GreensH fail at TF=Float64-32, $ok32, $e32, $l32")
    (e16, l16) = gtestH(G; TF = Float16)
    ok16 = (e16 < 1.e-14) && (l16 <= 7)
    ok16 || println("GreensH fail at TF=Float64-16")
    (e3216, l3216) = gtestH(G; TF = Float16, TW = Float32)
    ok3216 = (e3216 < 1.e-6) && (l3216 <= 4)
    ok3216 || println("GreensH fail at TF=Float32-16")
    heavyok = (ok16 && ok32 && ok3216)
    return heavyok
end

"""
greensEvsH(n=31)

Make sure heavy IR and expensive IR are give identical results
"""
function greensEvsH(n = 31)
    G = Gmat(n)
    (ee32, le32, he32) = gtestE(G; reshistout = true)
    (e32, l32, h32) = gtestH(G; reshistout = true)
    edel = abs(e32 - ee32) / e32
    hdel = norm(h32 - he32, Inf) / norm(he32, Inf)
    #errok=(e32==ee32)
    #histok=(h32==he32)
    errok = (edel < 1.e-2)
    histok = (hdel < 1.e-15)
    EvsHok32 = (errok && histok)
    EvsHok32 || println("TF=F32: heavy IR and expensive IR differ")
    (ee16, le16, he16) = gtestE(G; TF = Float16, reshistout = true)
    (e16, l16, h16) = gtestH(G; TF = Float16, reshistout = true)
    errok = (e16 == ee16)
    histok = (norm(h16 - he16,Inf) < 1.e-14)
    histok || println("he16 - h17 error  ", norm(he16 - h16, Inf))
    EvsHok16 = (errok && histok)
    (ee3216, le3216, he3216) = gtestE(G; TF = Float16, TW = Float32, reshistout = true)
    (e3216, l3216, h3216) = gtestH(G; TF = Float16, TW = Float32, reshistout = true)
    errok = (e3216 == ee3216)
    histok = (h3216 == he3216)
    histok = (norm(he3216 - he3216, Inf) < 1.e-14)
    EvsHok3216 = (errok && histok)
    EvsHPass = EvsHok32 && EvsHok16 && EvsHok3216
    EvsHPass || println("heavy vs expensive too far apart")
    return EvsHPass
end

"""
gtestH(G; TF=Float32,TW=Float64, reshistout=false)

Solve the inverse Laplacian problem with heavy IR
"""
function gtestH(G; TF = Float32, TW = Float64, reshistout = false)
    #G=Gmat(n)
    (n, n) = size(G)
    A = I + TW.(G)
    b = ones(TW, n)
    xe = A \ b
    MPA = MPHArray(A; TF = TF)
    MPF = mphlu!(MPA)
    soldata = \(MPF, b; reporting = true)
    xm = soldata.sol
    rhist = soldata.rhist
    nerr = norm(xm - xe, Inf)
    lenit = length(rhist)
    if reshistout
        return (nerr, lenit, rhist)
    else
        return (nerr, lenit)
    end
end

"""
greensok(n=31)

Test vanilla IR with the inverse Laplacian. 
"""
function greensok(n = 31)
    G = Gmat(n)
    (e32, l32) = gtest(G)
    ok32 = (e32 < 1.e-13) && (l32 <= 5)
    ok32 || println("Greens fail at TF=Float64-32, $e32, $l32")
    (e16, l16) = gtest(G; TF = Float16)
    ok16 = (e16 < 1.e-13) && (l16 <= 7)
    ok16 || println("Greens fail at TF=Float64-16: $e16, $l16")
    (e3216, l3216) = gtest(G; TF = Float16, TW = Float32)
    ok3216 = (e3216 < 5.e-6) && (l3216 <= 4)
    ok3216 || println("Greens fail at TF=Float32-16: $e3216, $l3216")
    lightok = (ok16 && ok32 && ok3216)
    return lightok
end


"""
gtest(G; TF=Float32,TW=Float64)

Solve the inverse Laplacian problem with IR
"""
function gtest(G; TF = Float32, TW = Float64)
    (n, n) = size(G)
    A = I + TW.(G)
    b = ones(TW, n)
    xe = A \ b
    MPA = MPArray(A; TF = TF)
    MPF = mplu!(MPA)
    soldata = \(MPF, b; reporting = true)
    xm = soldata.sol
    rhist = soldata.rhist
    nerr = norm(xm - xe, Inf)
    lenit = length(rhist)
    return (nerr, lenit)
end

"""
gtestE(G; TF=Float32,TW=Float64)

Solve the inverse Laplacian problem with expensive IR
This is only for CI
"""
function gtestE(G; TF = Float32, TW = Float64, reshistout = false)
    #G=Gmat(n)
    (n, n) = size(G)
    A = I + TW.(G)
    b = ones(TW, n)
    xe = A \ b
    MPA = MPArray(A; TF = TF, onthefly=true)
    MPF = mplu!(MPA)
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

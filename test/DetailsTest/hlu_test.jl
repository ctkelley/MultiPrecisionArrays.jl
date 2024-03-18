function hlu_test()
    N = 128
    h = 1.0 / (N + 1)
    X = collect(h:h:1-h)
    K = [ker(x, y) for x in X, y in X]
    A = I + 10.0* h * K
    Ah = Float16.(A)
    AFt = hlu(A)
    delnormA=norm(I - AFt\A)
    dok = (delnormA < 1.e-12)
    dok || println("delnormA = $delnormA")
    Aht = hlu(Ah)
    delah = norm(I - Aht \ Ah,Inf)*h
    hok = (delah < 1.e-5)
    hok || println("delah = $delah")
    #
    failout = testhlufail()
    failok = (failout == "Fact fails")
    #
    ndifflh = samefordouble(A)
    doubok = (ndifflh < 1.e-14)
    doubok || println("samefordouble fails $ndifflh")
    hluok = dok && hok && failok && doubok
    return hluok
end

function ker(x, y)
    ker = sin.(x - y)
    return ker
end

function testhlufail()
    A = ones(Float16,128, 128)
    try
        hlu!(A)
    catch
        return "Fact fails"
    end
end

function samefordouble(A)
    B = copy(A)
    BF = hlu!(B)
    AF = lu!(A)
    return norm(A - B, Inf)
end

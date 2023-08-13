function hlu_test()
    N = 100
    h = 1.0 / (N + 1)
    X = collect(h:h:1-h)
    K = [ker(x, y) for x in X, y in X]
    A = I + 0.1 * K
    Ah = Float16.(A)
    AFt = hlu(A)
    dok = (norm(I - AFt \ A) < 1.e-12)
    Aht = hlu(Ah)
    hok = (norm(I - Aht \ Ah) < 1.e-3)
    #
    failout = testhlufail()
    failok = (failout == "Fact fails")
    #
    ndifflh = samefordouble(A)
    doubok = (ndifflh < 1.e-14)
    hluok = dok && hok && failok && doubok
    return hluok
end

function ker(x, y)
    ker = sin.(x - y)
    return ker
end

function testhlufail()
    A = ones(2, 2)
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

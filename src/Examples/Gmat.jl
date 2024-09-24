"""
Gmat: discrete Green's operator for -d^2/dx^2 in 1D
"""
function Gmat(n, T = Float64)
    onet = T(1.0)
    h = onet / (n - onet)
    Ns = collect(0:1:n-1)
    X = h * Ns
    X[n] = onet
    G = [greens(x, y, onet) for x in X, y in X]
    G .*= h
    return G
end

function greens(x, y, onet)
    if x > y
        gf = y * (onet - x)
    else
        gf = x * (onet - y)
    end
    return gf
end

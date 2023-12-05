"""
Gmat: discrete Green's operator for -d^2/dx^2 in 1D
"""
function Gmat(n,T=Float64)
    onet=T(1.0)
# fudge fixes a corner case if T=Float32
    fudge=eps(T)*T(5.0)
    h = onet/ (n + onet)
    X = collect(T, h:h:1.0-h+fudge)
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

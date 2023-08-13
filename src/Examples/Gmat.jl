"""
Gmat: discrete Green's operator for -d^2/dx^2 in 1D
"""
function Gmat(n)
    h = 1.0 / (n + 1.0)
    X = collect(h:h:1.0-h)
    G = [greens(x, y) for x in X, y in X]
    G .*= h
    return G
end

function greens(x, y)
    if x > y
        gf = y * (1.0 - x)
    else
        gf = x * (1.0 - y)
    end
    return gf
end

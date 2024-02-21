function tir3p(A::AbstractArray{Float32,2}, xe::Vector{Float32}; otf = true)
    TL=Float16
#    TL = Float32
    TW = Float32
    TH = Float64
    nrmtype = Inf
    AL = TL.(A)
    ALF = lu!(AL)
    tol = 1.e-14
    itc = 0
    itmax = 20
    itmin = 4
    b = A * xe
    xl = ALF \ b
    AH = TH.(A)
    xh = AH\b;
    if otf
        x = TH.(xl)
        #    r = b - A * x
        r = TH.(b) - A * TH.(x)
    else
        x = xl
        r = b - A * x
    end
    d = copy(r)
    resnorm = norm(r, nrmtype)
    # IR Loop
    resnorm = norm(r, nrmtype) / norm(b, nrmtype)
    errnorm = norm(x - xh, nrmtype) / norm(xe, nrmtype)
    @printf("%1.5e   %1.5e \n", resnorm, errnorm)
    while ((resnorm > tol) || (itc < itmin)) && itc < itmax
        if otf
            nrmr = norm(r, nrmtype)
            rs = r / nrmr
            #        d = ALF \ rs
            ldiv!(d, ALF, rs)
            d .*= nrmr
            x .+= d
            #       r .= b - A * x
            mul!(r, A, -x)
            r .+= b
        else
            nrmr = norm(r, nrmtype)
            rs = r / nrmr
            d .= ALF \ rs
            d .*= nrmr
            x .+= d
            r .= b - A * x
        end
        resnorm = norm(r, nrmtype) / norm(b, nrmtype)
        dnorm = norm(d, nrmtype) / norm(x, nrmtype)
        errnorm = norm(x - xh, nrmtype) / norm(xe, nrmtype)
        @printf("%1.5e   %1.5e   %1.5e \n", resnorm, errnorm, dnorm)
        itc += 1
    end
    return x
end

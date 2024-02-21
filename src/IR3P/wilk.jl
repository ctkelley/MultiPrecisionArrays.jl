function wilk(A::AbstractArray{Float32,2}, xe::Vector{Float32})
    TE = Float64
    TW = Float32
    nrmtype = Inf
    tol = 1.e-14
    ALF = lu(A)
    tol = 1.e-14
    itc = 0
    itmax = 10
    itmin = 2
    b = A * xe
    AH=TE.(A)
    bh = AH*xe
    xeh = AH\b
    xl = ALF \ b
    # Initialize
    x = TE.(xl)
#    x = xl 
#    r = TE.(b) - A * TE.(x)
#    r = muladd(A, -x, TE.(b))
    r = b - A*x
    d = copy(r)
    resnorm = norm(r, nrmtype) / norm(b, nrmtype)
    errnorm = norm(x - xe, nrmtype) / norm(xe, nrmtype)
    resnorm = norm(r, nrmtype)/norm(b,nrmtype)
    @printf("%1.5e   %1.5e \n", resnorm, errnorm)
    while ((resnorm > tol) || (itc < itmin)) && itc < itmax
        ldiv!(d, ALF, r)
        x .+= TW.(d)
        mul!(r, A, -x)
        r .+= b
        resnorm = norm(r, nrmtype) / norm(b, nrmtype)
        dnorm = norm(d, nrmtype) / norm(x, nrmtype)
        errnorm = norm(x - xeh, nrmtype) / norm(xe, nrmtype)
        @printf("%1.5e   %1.5e   %1.5e \n", resnorm, errnorm, dnorm)
        itc += 1
    end
    println(norm(x-xe,Inf))
    return x
end

function resid(A, x, b)
n=length(x)
r = Float64.(b)
for j=1:n
    for i=1:n
       r[i] -= Float64(A[i,j])*x[j]
    end
end
return r
end

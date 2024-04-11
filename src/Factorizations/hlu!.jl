"""
hlu!(A::AbstractMatrix{T}) where {T}
Return LU factorization of A

C. T. Kelley, 2023

This function is a hack of generic_lufact! which is part of

https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/lu.jl

I "fixed" the code to be Float16 only and fixed pivoting
to only MaxRow.

All I did in the factorization
was thread the critical loop with FLoops.@floop and
put @simd in the inner loop. For larger problems (n > 128)
these changes got me a 2-10x speedup
on my Mac M2 Pro with 8 performance cores. I'm happy.

"""
function hlu!(A::AbstractMatrix{T}) where {T}
    pivot = RowMaximum()
    T == Float16 || @warn("Use hlu for half precision only!")
    LAPACK.chkfinite(A)
    # Extract values and make sure the problem is square
    m, n = size(A)
    # Small n? Revert to normal lu
    (n < 128) && (AF=lu!(A); return AF )
    minmn = min(m, n)
    # Initialize variables
    info = 0
    BlasInt = LinearAlgebra.BLAS.BlasInt
    ipiv = Vector{BlasInt}(undef, minmn)
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if k < m
                amax = abs(A[k, k])
                for i = k+1:m
                    absi = abs(A[i, k])
                    if absi > amax
                        kp = i
                        amax = absi
                    end
                end
                #            elseif pivot === RowNonZero()
                #                for i = k:m
                #                    if !iszero(A[i, k])
                #                        kp = i
                #                        break
                #                    end
                #                end
            end
            ipiv[k] = kp
            if !iszero(A[kp, k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k, i]
                        A[k, i] = A[kp, i]
                        A[kp, i] = tmp
                    end
                end
                # Scale first column
                Akkinv = inv(A[k, k])
                for i = k+1:m
                    A[i, k] *= Akkinv
                end
            elseif info == 0
                info = k
            end
            # Update the rest
            ntasks=min(nthreads(), 1 + floor(Int,(n-k)/8))
                  tforeach(k+1:n; ntasks=ntasks) do j
                    Akj = -A[k, j]
                    @inbounds @simd ivdep for i = k+1:m
                        A[i, j] += A[i, k] * Akj
                    end # i loop
                end #j loop
        end
    end
    checknonsingular(info, pivot)
    return LU{T,typeof(A),typeof(ipiv)}(A, ipiv, convert(BlasInt, info))
end

function hlu(A)
    C = copy(A)
    AF = hlu!(C)
    return AF
end

# More stuff I got from Base
checknonsingular(info, ::RowMaximum) = info == 0 || throw(SingularException(info))

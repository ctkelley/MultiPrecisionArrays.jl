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
put @simd in the inner loop. These changes got me a 10x speedup
on my Mac M2 Pro with 8 performance cores. I'm happy.

"""
function hlu!(A::AbstractMatrix{T}) where {T}
    pivot = RowMaximum()
    T == Float16 || @warn("Use hlu for half precision only!")
    LAPACK.chkfinite(A)
    # Extract values and make sure the problem is square
    m, n = size(A)
    minmn = min(m, n)
    (m == n) || @warn("hlu is only tested for square problems")
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
            #            @tasks for j = k+1:n
            # This is where I need the problem to be square
            ntasks=nthreads()
#            if k < minmn
            ntasks=min(nthreads(), 1 + floor(Int,(n-k)/8))
                @tasks for j = k+1:n
                   @set begin
                        scheduler=DynamicScheduler(; nchunks=ntasks)
                    end
                    Akj = -A[k, j]
                    @simd ivdep for i = k+1:m
                        @inbounds A[i, j] += A[i, k] * Akj
                        # @inbounds A[i, j] -= A[i, k] * A[k, j]
                        # @inbounds A[i,j] = muladd(A[i,k],Akj,A[i,j])
                    end # i loop
                end #j loop
#            end # k < minmn if statement
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

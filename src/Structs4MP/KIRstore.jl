"""     
KIRstore(n, lsolver, TR=Float64)

Preallocates the vectors a Krylov method uses internally.
"""
function KIRstore(n, lsolver, TR=Float64)
    tmp1 = zeros(TR,n)
    tmp2 = zeros(TR,n)
    tmp3 = zeros(TR,n)
    tmp4 = zeros(TR,n)
    if lsolver == "gmres"
        return (tmp1, tmp2, tmp3, tmp4)
    else
        tmp5 = zeros(TR,n)
        tmp6 = zeros(TR,n)
        tmp7 = zeros(TR,n)
        return (tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)
    end
end


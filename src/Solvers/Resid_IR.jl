"""

Residual computation for MultiPrecisionArrays. AF is the multiprecision
factorization. 

xres is TR.(current iteration). I only allocate storage for this
if TR > TW

oneb = TR.(1.0)
"""
function Resid_IR(r, x, xres, b, AF)
    TR = eltype(r)
    TW = eltype(b)
    #
    if TR==TW
        mul!(r, AF.AH, x)
    else
        xres .= TR.(x)
        mul!(r, AF.AH, xres)
    end
    #
    # After mul! the residual is overwritten with Ax 
    #
    oneb = TR(1.0)
    r .*= -oneb
    #
    # Now r = - Ax 
    axpy!(oneb, b, r)
    # and now the residual is b-Ax like it needs to be
    #
    return r
end

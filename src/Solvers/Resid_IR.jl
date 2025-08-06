"""

Residual computation for MultiPrecisionArrays. AF is the multiprecision
factorization. 

rloop and xloop are the residual vector and TR.(current iteration). 

bsc = TR.(b), oneb = TR.(1.0)
"""
function Resid_IR(r, xloop, bsc, AF)
    TR = eltype(r)
    #
    mul!(r, AF.AH, xloop)
    #
    # After mul! the residual is overwritten with Ax 
    #
    oneb = TR(1.0)
    r .*= -oneb
    #
    # Now r = - Ax 
    axpy!(oneb, bsc, r)
    # and now the residual is b-Ax like it needs to be
    #
    return r
end

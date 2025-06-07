function testnorm(delnorm, TA, alpha)
#
# This function is used for a couple problems where the convergence
# of GMRES-IR is poor if TW=TR=Float32 and TF=Float64. That is why
# the tolerances are large in that case. 
#
    if TA == Float64
        normok = (alpha < 750.0) ? (delnorm < 1.e-13) : (delnorm < 1.e-10)
    else
        normok = (alpha < 750.0) ? (delnorm < 1.e-5) : (delnorm < 5.e-3)
    end
    normok || println("delnorm= $delnorm ; precision= $TA, alpha = $alpha")
    return normok
end


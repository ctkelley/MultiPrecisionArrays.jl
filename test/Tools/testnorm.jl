function testnorm(delnorm, TA, alpha)
    if TA == Float64
        normok = (alpha < 750.0) ? (delnorm < 1.e-13) : (delnorm < 1.e-10)
    else
        normok = (alpha < 750.0) ? (delnorm < 1.e-6) : (delnorm < 5.e-3)
    end
    normok || println("delnorm= $delnorm ; precision= $TA, alpha = $alpha")
    return normok
end


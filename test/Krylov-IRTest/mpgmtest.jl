function mpgmtest(n = 100)
    normok = true
    for TA in [Float64, Float32]
        normok = normok && get_gmir_results(n, TA)
    end
    normok = normok && get_gmir_results(n, Float64; test6416 = true)
    return normok
end

function get_gmir_results(n, TA; test6416 = false, debug = false)
    TAok = true
    for alpha in [1.0, 10.0, 800.0]
        AD = TA.(I - alpha * Gmat(n))
        xs = ones(TA, n)
        b = AD * xs
        if test6416
            ADM = MPHArray(AD; TL = Float16)
        else
            ADM = MPHArray(AD)
        end
        ADF = mpglu!(ADM)
        zt = mpgmir(ADF, b; reporting = true)
        z = zt.sol
        delnorm = norm(z - xs, Inf)
        if debug
            its = length(zt.rhist)
            rnorm = zt.rhist[its]
            println("alpha=$alpha, delnorm=$delnorm, rnorm=$rnorm, its = $its")
        end
        TAok = TAok && testnorm(delnorm, TA, alpha)
    end
    return TAok
end

function testnorm(delnorm, TA, alpha)
    if TA == Float64
        normok = (alpha < 750.0) ? (delnorm < 1.e-13) : (delnorm < 1.e-10)
    else
        normok = (alpha < 750.0) ? (delnorm < 1.e-5) : (delnorm < 1.e-2)
    end
    normok || println("delnorm= $delnorm ; precision= $TA, alpha = $alpha")
    return normok
end

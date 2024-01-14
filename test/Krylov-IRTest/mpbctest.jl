function mpbctest(n = 100)
    normok = true
    for TA in [Float64, Float32]
        normok = normok && get_bcir_results(n, TA)
    end
    normok = normok && get_bcir_results(n, Float64; test6416 = true)
    return normok
end

function get_bcir_results(n, TA; test6416 = false, debug = false)
    TAok = true
    for alpha in [1.0, 10.0, 800.0]
        AD = TA.(I - alpha * Gmat(n));
        xs = ones(TA, n);
        b = AD * xs;
        if test6416
            ADM = MPBArray(AD; TL = Float16)
        else
            ADM = MPBArray(AD)
        end
        ADF = mpblu!(ADM)
        zt = mpbcir(ADF, b; reporting = true)
        wt = ADF\b;
        slasherr=norm(zt.sol-wt,Inf)
        slashok = (slasherr < 1.e-17)
        z = zt.sol
        delnorm = norm(z - xs, Inf)
        if debug
            its = length(zt.rhist)
            rnorm = zt.rhist[its]
            println("alpha=$alpha, delnorm=$delnorm, rnorm=$rnorm, its = $its")
        end
        TAok = TAok && testnorm(delnorm, TA, alpha) && slashok
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

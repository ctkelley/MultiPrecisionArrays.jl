function mpgmtest(n = 100)
    normok = true
    for TA in [Float64, Float32]
        normok = normok && get_gmir_results(n, TA)
        normok || println("mpgm failure TA = $TA")
    end
    ok6416 = get_gmir_results(n, Float64; test6416 = true)
    ok6416 || println("6416 failure")
    normok = normok && ok6416
    return normok
end

function get_gmir_results(n, TA; test6416 = false, debug = false)
    TAok = true
    for alpha in [1.0, 10.0, 800.0]
        AD = TA.(I - alpha * Gmat(n));
        xs = ones(TA, n);
        b = AD * xs;
        if test6416
            ADM = MPHArray(AD; TF = Float16);
        else
            ADM = MPHArray(AD);
        end
        ADF = mpglu!(ADM);
        zt = mpkrir(ADF, b; reporting = true);
        z = zt.sol;
        delnorm = norm(z - xs, Inf)
        if debug
            its = length(zt.rhist)
            rnorm = zt.rhist[its]
            println("alpha=$alpha, delnorm=$delnorm, rnorm=$rnorm, its = $its")
        end
        TAok = TAok && testnorm(delnorm, TA, alpha)
        TAok || println("delnorm = $delnorm")
    end
    return TAok
end

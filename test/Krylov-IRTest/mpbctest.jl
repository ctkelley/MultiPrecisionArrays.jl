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
            ADM = MPBArray(AD; TF = Float16)
        else
            ADM = MPBArray(AD)
        end
        ADF = mpblu!(ADM);
        zt = mpkrir(ADF, b; reporting = true);
        wt = ADF\b;
        slasherr=norm(zt.sol-wt,Inf)
        slashok = (slasherr < 1.e-17)
        z = zt.sol;
        delnorm = norm(z - xs, Inf)
        if debug
            its = length(zt.rhist)
            rnorm = zt.rhist[its]
            println("alpha=$alpha, delnorm=$delnorm, rnorm=$rnorm, its = $its")
        end
        TAok = TAok && testnorm(delnorm, TA, alpha) && slashok
        TAok || println("mpkrir fails; $TA, $delnorm, $alpha")
    end
    return TAok
end

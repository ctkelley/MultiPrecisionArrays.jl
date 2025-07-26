function IRKsolve(x0, r, MPhatv, AF, eta, MP_Data, ktype)
    #
    # Solve the correction equation with a Krylov method
    #
    VF = AF.VStore
    kl_store = AF.KStore
    if ktype == "GMRES"
        kout = kl_gmres(
            x0,
            r,
            MPhatv,
            VF,
            eta,
            MPhptv;
            pdata = MP_Data,
            side = "left",
            kl_store = kl_store,
        )
    elseif ktype == "BiCGSTAB"
        kout = kl_bicgstab(
            x0,
            r,
            MPhatv,
            VF,
            eta,
            MPhptv;
            pdata = MP_Data,
            side = "left",
            kl_store = kl_store,
        )
    end
    return kout
end

# Matrix-vector product for Krylov-IR
function MPhatv(x, pdata)
    atv = pdata.atv
    mul!(atv, pdata.MPF.AH, x)
    return atv
end

# Preconditioner-vector product for Krylov-IR
function MPhptv(x, pdata)
    TF=pdata.TF
    TW=pdata.TW
    #
    # ldiv! is not doing well for one case, so I hide from it.
    #
    if (TF == Float16) && (TW == Float32)
        x .= pdata.MPF.AF\x;
    else
        ldiv!(pdata.MPF.AF, x)
    end
    return x
end




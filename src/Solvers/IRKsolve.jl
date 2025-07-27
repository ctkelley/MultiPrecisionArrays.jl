function IRKsolve!(AF::Union{MPKFact,MPGHFact}, r, rs ; MP_Data =[])
#function IRKsolve!(AF::Union{MPKFact,MPGHFact}, r, rs, x0, eta, MP_Data, ktype)
    #
    # Solve the correction equation in working precision with a Krylov method
    #
    TR = eltype(r)
    TW = eltype(AF.sol)
    (TR == TW) ? zs=r : (rs .= r; zs = rs)
    #
    ktype = MP_Data.ktype
    eta = MP_Data.eta
    x0 = MP_Data.x0
    VF = AF.VStore
    kl_store = AF.KStore
    if ktype == "GMRES"
        kout = kl_gmres(
            x0,
            zs,
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
            zs,
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
    TW = eltype(pdata.x0)
    TF = eltype(pdata.MPF.AF)
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




function Types_IR_Init(AF, b)
    TB=eltype(b)
    TF=eltype(AF.AF)
    TW=eltype(AF.AH)
    r = AF.residual
    TR=eltype(r)
    #
    # Are the precisions consistent? If not, I have a bug somewhere.
    # Otherwise, set the tolerance on the iteration to 10*eps.
    # If the iteration can't meet the tolerance, terminate when
    # the residual norms stagnate (res_old > Rmax * res_new)
    #
    TFact = consistency(AF, TF, TW, TB)
    return (TW, TF, TR, TFact)
end

function Term_Init(AF, term_parms, bnrm)
    #   
    #  get the termination data
    #   
    TR = eltype(AF.residual)
    TW = eltype(bnrm)
    tolf = termination_settings(AF, term_parms)
    Rmax = term_parms.Rmax
    litmax = term_parms.litmax
    # norm(x) = 0 on initialization
    tol = tolf * bnrm
    return (tolf, Rmax, litmax, tol)
end


function Solver_IR_Init(AF, b, normtype)
    r = AF.residual
    TR=eltype(r)
    TW=eltype(b)
    r .= TR.(b)
    x = AF.sol
    x .*= TW(0.0)
    xnrm = norm(x, normtype)
    bnrm = norm(b, normtype)
    TF=eltype(AF.AF)
    onthefly = AF.onthefly
    HiRes = (eps(TR) < eps(TW))
    HiRes && (onthefly = true)
    anrm = AF.anrm
    #
    # rs goes into the triangular solve to compute the correction
    # So if TW==TR and onthefly=true, you can use rs = r
    # If TW==TR and onthefly=false, then rs = TF.(r)
    # If TR > TW then rs = TW.(r) 
    #
    rs = r
    (TR == TW) || (rs = TW.(r))
    onthefly || (rs = TF.(r))
    rhist = Vector{TR}()
    dhist = Vector{TW}()
    push!(rhist, bnrm)
    return (x, r, rs, bnrm, anrm, rhist, dhist)
end



function consistency(AF, TF, TW, TB)
    #
    # Test for errors in the structures. If this fails I have done something
    # fundamentally wrong.
    #
    if is_heavy(AF)
        TFact = eltype(AF.AH)
    else
        TFact = eltype(AF.AL)
    end
    # TFact is the precision of the factors; should be TF
    # unless we're dealing with a heavy MPArray for CI
    #
    TFok = ((TF == TFact) || (typeof(AF) == MPHHFact))
    TFok || error("TF is supposed to be TFact")
    (TW == TB) || error("inconsistent precisions; A and b must have same type")
    return TFact
end

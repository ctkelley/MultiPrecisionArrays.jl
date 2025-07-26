function Types_IR_Init(AF, b, normtype = Inf)
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

function Solver_IR_Init(AF, b)
    r = AF.residual
    TR=eltype(r)
    TW=eltype(b)
    r .= TR.(b)
    x = AF.sol
    x .*= TW(0.0)
    TF=eltype(AF.AF)
    onthefly = AF.onthefly
    HiRes = (eps(TR) < eps(TW))
    HiRes && (onthefly = true)
    anrm = AF.anrm
    onthefly ? (rs = ones(TF, 1)) : (rs = zeros(TF, size(b)))
    #
    # If TR > TW then do the solve in TW after computing r in TR
    #
    (TR == TW) || (rs = zeros(TW, size(b)))
    return (x, r, rs, anrm, onthefly, HiRes)
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

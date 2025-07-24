function IR_Init(AF,b,normtype=Inf)
TB=eltype(b)
TF=eltype(AF.AF)
TW=eltype(AF.AH)
r = AF.residual
TR=eltype(r)
x = AF.sol
x .*= TW(0.0)
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
  #
    # Are the precisions consistent? If not, I have a bug somewhere.
    # Otherwise, set the tolerance on the iteration to 10*eps.
    # If the iteration can't meet the tolerance, terminate when
    # the residual norms stagnate (res_old > Rmax * res_new)
    #
    (TW == TB) || error("inconsistent precisions; A and b must have same type")
return (x, r, TB, TW, TF, TR, TFact)
end

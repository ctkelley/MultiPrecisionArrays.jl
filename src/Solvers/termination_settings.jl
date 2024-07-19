function termination_settings(TR, residterm)
#
# I am continually tweaking this stuff
#
Cr=1.0; Ce=.5;
residterm ? tf=Cr : tf=Ce
tolf = eps(TR)*tf
#
#   I'm using the L1 norm because it's much faster.
#
residterm ?  anrm = 0.0 : anrm = opnorm(AD, 1)
term_out = (tolf=tolf, anrm=anrm)
return term_out
end

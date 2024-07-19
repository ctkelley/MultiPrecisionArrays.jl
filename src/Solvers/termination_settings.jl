function termination_settings(TR, residterm)
#
# I am continually tweaking this stuff
#
residterm ? tf=1.0 : tf=.5
tolf = eps(TR)*tf
#
#   I'm using the L1 norm because it's much faster.
#
residterm ?  anrm = 0.0 : anrm = opnorm(AD, 1)
term_out = (tolf=tolf, anrm=anrm)
return term_out
end

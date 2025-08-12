# Interprecision Transfers: the triangular solve

This default way to use the low precision factorization is to 
simply do
```
d = norm(r) *TW.( AF\(r /norm(r)))
```
Since ```AF``` is stored in ```TF``` and the residual is stored in
precision ```TR```, an interprecision transfer happens during the triangular
solve and the result is in precision ```TR```. We refer to this approach
as interprecision transfer on the fly or a mixed precision
solve (MPS). 

Another way to do this is to downcase the residual before the solve. In
this low precision solve (LPS) case have
```
d = norm(r) *TW.( AF\ TF.(r /norm(r))))
```

LPS saves $N^2$ interprecision transfers and one might think that
this savings could be useful. We look into this in the appendix of
[ctk:mparraysdocs](@cite) and find that LPS has no meaningful
advantage. This is a change in my long-held incorrect opinion and
MPS is the default. You might read the discussion on this issue and
see if using LPS would be useful for you. The ```onthefly``` keyword
argument of ```mplu``` will do this.


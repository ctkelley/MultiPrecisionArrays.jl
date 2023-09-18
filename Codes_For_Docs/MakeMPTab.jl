#
# This file makes the table comparing the triangular solve options
# for computer time.
#
"""
MakeMPTab(m = 4, texok = false; T = Float64)
Make a nice table of the timings for the solves.
"""
function MakeMPTab(m = 4, texok = false; T = Float64)
    AT = zeros(m, 7)
    p = collect(0:1:m-1)
    tp = 2 .^ p
    np = 512 .* tp
    for idim = 1:m
        n = np[idim]
        AT[idim, 2:7] = pitch(n; T = T)
        AT[idim, 1] = n
    end
    T==Float64 ?  LULow = "LU32" : LULow="LU16"
    T==Float64 ?  MVHigh = "MV64" : MVHigh="MV32"
    @printf(
        "%5s     %5s      %5s       %3s        %3s        %3s   %8s\n",
        "n",
        "$MVHigh",
        "$LULow",
#        "LU32",
        "HPS",
        "MPS",
        "LPS",
        "$LULow/MPS"
    )
    for idim = 1:m
        @printf(
            "%5d    %5.1e    %5.1e    %5.1e    %5.1e    %5.1e %5.1e\n",
            AT[idim, 1],
            AT[idim, 2],
            AT[idim, 3],
            AT[idim, 4],
            AT[idim, 5],
            AT[idim, 6],
            AT[idim, 7]
        )
    end
    headers = ["N", "$MVHigh", "$LULow", "HPS", "MPS", "LPS", "$LULow/MPS"]
    formats = ("%5d  &  %5.1e &   %5.1e  &  %5.1e &  %5.1e  & %5.1e & %5.1e")
    if texok
        fprintTeX(headers, formats, AT)
    end
    return AT
end

"""
pitch(n; T = Float64)
Collect the timings for the triangular solve options.
"""
function pitch(n; T = Float64)
    G=Gmat(n); 
    x=ones(n); 
    A = I + 800.0*G
    if T == Float64
        x=ones(T,n)
        bh=A*x
        AH = T.(A)
        AS = Float32.(AH)
        AS2 = copy(AS)
        bs = Float32.(bh)
    else
        x=ones(T,n)
        bh=A*x
        AH = T.(A)
        AS = Float16.(AH)
        AS2 = copy(AS)
        bs = Float16.(bh)
    end
    AHF = lu(AH)
    ASF = lu(AS)
#
# Prints the 
#
#    tluh = @belapsed lu(AV) setup = (AV = copy($AH)) evals=1
    tluh = @belapsed $AH*$bh 
    tlus = @belapsed lu(AVS) setup = (AVS = copy($AS)) evals=1
    thsol = @belapsed $AHF \ $bh
    tssol = @belapsed $ASF \ $bh
    tgssol = @belapsed $ASF \ $bs
    ratio=tlus/tssol
#
# Prints the play-by-play to the REPL for my debugging.
#
playbyplay=false
if playbyplay
    println("High precsion LU time = $tluh")
    println("Low precision LU time = $tlus")
    println("High precision solve time = $thsol")
    println("MPS time = $tssol")
    println("LPS time = $tgssol")
end
    rvec = [tluh tlus thsol tssol tgssol ratio]
    return rvec
end

#
# This file makes the table comparing the triangular solve options
# It's the first table in the talk and Table 1.xx in the paper
#

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
    @printf(
        "%5s     %5s      %5s       %3s        %3s        %3s   %8s\n",
        "n",
        "LU64",
        "LU32",
        "HPS",
        "MPS",
        "LPS",
        "LU32/MPS"
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
    headers = ["N", "LU64", "LU32", "HPS", "MPS", "LPS", "LU32/MPS"]
    formats = ("%5d  &  %5.1e &   %5.1e  &  %5.1e &  %5.1e  & %5.1e & %5.1e")
    if texok
        fprintTeX(headers, formats, AT)
    end
    return AT
end

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
    tluh = @belapsed lu(AV) setup = (AV = copy($AH))
    println("High precsion LU time = $tluh")
    tlus = @belapsed lu(AVS) setup = (AVS = copy($AS))
    println("Low precision LU time = $tlus")
    thsol = @belapsed $AHF \ $bh
    println("High precision solve time = $thsol")
    tssol = @belapsed $ASF \ $bh
    println("MPS time = $tssol")
    tgssol = @belapsed $ASF \ $bs
    println("LPS time = $tgssol")
    ratio=tlus/tssol
    rvec = [tluh tlus thsol tssol tgssol ratio]
    return rvec
end

function HalfTime(p = 3)
    pv = collect(1:1:p)
    ppv = 2 .^ pv
    nv = 512 .* ppv
    Time_F = zeros(p, 5)
    headers = ["N", "F64", "F32", "F16", "F16/F64"]
    for p in pv
        N = nv[p]
        G = Gmat(N)
        AD = I - 800.0*G
        AS = Float32.(AD)
        AH = Float16.(AD)
        td = @belapsed lu($AD)
        ts = @belapsed lu($AS)
        th = @belapsed hlu($AH)
        rt = th / td
        Time_F[p, :] = [N, td, ts, th, rt]
    end
    dformat = "%9d %9.2e %9.2e %9.2e %9.2e \n"
    hformat = "%7s %9s %9s %9s %11s \n"
    printf(fmt::String, args...) = @eval @printf($fmt, $(args...))
    printf(hformat,headers...)
#    for i = 1:5
#        printf("%9s", headers[i])
#    end
#    printf("\n")
    for p in pv
        printf(dformat, Time_F[p, :]...)
    end
    return Time_F
end

function timetab(TabMac)
    headers = ["N", "Double", "Single", "Half", "Ratio"]
    formats = "%d & %7.2e & %7.2e & %7.2e & %7.2e"
    fprintTeX(headers, formats, TabMac)
end

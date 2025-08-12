function wilk_krylov(n = 31, a = 800.0; basissize = 3)
    #
    # Give GMRES a very small Krylov space to make it take more iterations
    #
    G=Gmat(n, Float32);
    #
    # alpha=1.0 = well conditioned
    # alpha=800.0 = moderately ill conditioned
    #
    alpha=Float32(a)
    A=I + alpha*G;
    b=ones(Float32, n);
    AD=Float64.(A);
    bd=Float64.(b);
    # solve the promoted problem
    xp = AD\bd;
    # Set it up for IR-GMRES
    AFG=mpglu(A; TR = Float64, basissize = basissize)
    mgout=\(AFG, b; reporting = true);
    lmg=length(mgout.rhist)
    promdiffg = norm(xp-mgout.sol, Inf)
    #
    # I want to see somewhere between 3 and 7 iterations
    # and an error < 10^{-7}
    #
    gmresok = (3 <= lmg <= 7) && (promdiffg < 1.e-7)
    #println(lmg,"  ",promdiffg,"  ",gmresok)
    gmresok || println("IRGMRES error. lmg=$lmg, promdiffg=$promdiffg")
    ABG=mpblu(A; TR = Float64);
    mbout=\(ABG, b; reporting = true);
    #
    # BiCGSTAB will convege to the defect so you'll need only one iteration.
    # Making the problem larger will result in a more interesting computation,
    # but I'm doing CI here and it must be fast.
    #
    lmb=length(mbout.rhist)
    promdiffb = norm(xp-mbout.sol, Inf)
    bicgstabok = (1 <= lmb <= 5) && (promdiffb < 1.e-7)
    bicgstabok || println("IRBiCGSTAB error. lmb=$lmb, promdiffb=$promdiffb")
    #println(lmb,"  ",promdiffb,"  ",bicgstabok)
    #results=(mgout=mgout, mbout=mbout)
    IRKrylovok = (gmresok && bicgstabok)
    return IRKrylovok
end

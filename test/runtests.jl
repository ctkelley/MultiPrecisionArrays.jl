using Test, MultiPrecisionArrays
using MultiPrecisionArrays.Examples: Gmat
import MultiPrecisionArrays: TERM
using LinearAlgebra: I, lu!, norm
using SIAMFANLEquations: nsol
using SIAMFANLEquations.TestProblems: heqinit, heqf!, heqJ!
using StaticArrays

import MultiPrecisionArrays.MPArray
import MultiPrecisionArrays.MPHArray
import MultiPrecisionArrays.MPGArray
import MultiPrecisionArrays.MPBArray
import MultiPrecisionArrays.MPLFact
import MultiPrecisionArrays.mpkrir

include("Tools/testnorm.jl")
include("Greens/gtest.jl")
include("NLTest/nltest.jl")
include("DetailsTest/precision_test.jl")
include("DetailsTest/hlu_test.jl")
include("DetailsTest/slashtest.jl")
include("DetailsTest/mplu_test.jl")
include("DetailsTest/AbsArray.jl")
include("DetailsTest/static_test.jl")
include("DetailsTest/term_test.jl")
include("Krylov-IRTest/mpbctest.jl")
include("Wilkinson/wilk_test.jl")
include("Wilkinson/wilk_krylov.jl")
#
# Using MPHArray for CI only. This may go soon.
#
include("MPHArray-Test/gtesth.jl")
include("MPHArray-Test/mpgmtest.jl")
include("MPHArray-Test/hvse.jl")
include("MPHArray-Test/precision_testH.jl")

println("starting")

@testset "Greens Functions" begin
    @test greensok()
#
# Using MPHArray for CI only. This may go soon.
#
    @test greensHok()
    @test greensEvsH()
end

@testset "Nonlinear Equations" begin
    @test nltest()
end

@testset "Details" begin
    @test precision_test()
    @test hlu_test()
    @test slashtest()
    @test mplu_test()
    @test mpglu_test()
    @test mpblu_test()
    @test AbsArray()
    @test static_test()
    @test term_test()
    @test test_term_parms()
#
# Using MPHArray for CI only. This may go soon.
#
    @test precision_testH()
end

@testset "Krylov-IR" begin
    @test mpbctest(1000)
end
#
# Using MPHArray for CI only. This may go soon.
#
@testset "Krylov-IR2" begin
    @test mpgmtest(1000)
end
@testset "Krylov-IR3" begin
    @test hvse(128)
end

@testset "High Residual Precision" begin
   @test hi_precision_res()
   @test wilk_krylov()
end

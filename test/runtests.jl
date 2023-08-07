using MultiPrecisionArrays
using MultiPrecisionArrays.Examples
using LinearAlgebra
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using Test

include("Greens/gtest.jl")
include("NLTest/nltest.jl")
include("DetailsTest/precision_test.jl")
include("DetailsTest/hlu_test.jl")
include("GM-IRTest/mpgmtest.jl")

println("starting")

@testset "Greens Functions" begin
   @test greensok();
   @test greensHok();
   @test greensEvsH();
end

@testset "Nonlinear Equations" begin
   @test nltest();
end

@testset "Details" begin
   @test precision_test();
   @test hlu_test();
end

@testset "GM-IR" begin
   @test mpgmtest(1000)
end

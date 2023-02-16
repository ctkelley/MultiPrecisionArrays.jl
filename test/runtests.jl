using MPArrays
using MPArrays.Examples
using LinearAlgebra
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using Test

include("Greens/gtest.jl")
include("NLTest/nltest.jl")
include("DetailsTest/precision_test.jl")

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
end

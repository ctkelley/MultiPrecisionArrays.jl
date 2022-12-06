using MPArrays
using MPArrays.Examples
using LinearAlgebra
using SIAMFANLEquations
using SIAMFANLEquations.TestProblems
using Test

include("Greens/gtest.jl")
include("NLTest/nltest.jl")

@testset "Greens Functions" begin
   @test greensok();
   @test greensHok();
end

@testset "Nonlinear Equations" begin
   @test nltest();
end

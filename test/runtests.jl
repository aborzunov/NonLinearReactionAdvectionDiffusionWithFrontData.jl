using NonLinearReactionAdvectionDiffusionWithFrontData
using Test
using Documenter

@testset "Initial condition" begin
    N = 20;
    X = [i/N for i in 0:N];
    Y = NonLinearReactionAdvectionDiffusionWithFrontData.u_init.(X);

    @test length(Y) == N+1
    @test isapprox(Y[1], -8, rtol=0.001)
    @test isapprox(Y[end], 4, rtol=0.001)
end

@testset "DocTests" begin
#    doctest(NonLinearReactionAdvectionDiffusionWithFrontData)
end

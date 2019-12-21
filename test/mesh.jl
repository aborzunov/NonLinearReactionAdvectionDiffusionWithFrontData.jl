@testset "Mesh Creation" begin

    @testset "Uniform" begin
        N = 10;
        l = -1;
        r = 1;

        X, h = NonLinearReactionAdvectionDiffusionWithFrontData.create_mesh(l, r, N);

        @test length(X) == N+1
        @test isapprox(X, collect(range(l,r, length=N+1)))
        @test isapprox(h, (r - l)/N)
    end
end



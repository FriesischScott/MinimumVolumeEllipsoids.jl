@testset "rand" begin
    ϵ = Ellipsoid(PDMat([2 0; 0 2]), [0, 0])

    Random.seed!(8128)
    U = rand(ϵ, 1000000)
    Random.seed!()

    @test size(U) == (2, 1000000)
    @test maximum(U; dims=2) ≈ [1.0, 1.0] atol = 1e-3
    @test minimum(U; dims=2) ≈ [-1.0, -1.0] atol = 1e-3
    @test mean(U; dims=2) ≈ [0.0; 0.0] atol = 1e-3

    ϵ = Ellipsoid(PDMat([2 0; 0 2]), [0.5, 0.5])

    Random.seed!(8128)
    U = rand(ϵ, 1000000)
    Random.seed!()

    @test mean(U; dims=2) ≈ [0.5, 0.5] atol = 1e-2
end

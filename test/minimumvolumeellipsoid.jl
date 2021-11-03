@testset "Minimum Volume Ellipsoid" begin
    X = [
        -1 -1 1 2
        1 -1 -1 2
    ]

    @testset "Centered" begin
        ϵ = minimum_volume_ellipsoid(X, 1e-10; centered=true)

        @test ϵ.x == [0, 0]
        @test ϵ.H ≈ [5/8 -3/8; -3/8 5/8] atol = 1e-10
    end

    @testset "Not Centered" begin
        ϵ = minimum_volume_ellipsoid(X, 1e-10)

        @test ϵ.x ≈ [0.5, 0.5] atol = 1e-10
        @test ϵ.H ≈ [2/3 -2/9; -2/9 2/3] atol = 1e-10
    end
end

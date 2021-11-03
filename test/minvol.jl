@testset "_minvol" begin
    Y = [
        -1 -1 1 2
        1 -1 -1 2
    ]

    @testset "centered" begin
        u, R = MinimumVolumeEllipsoids._minvol(Y, 1e-10)

        @test u == [0.5, 0, 0, 0.5]
        @test R.L * R.L' ≈ [5/2 3/2; 3/2 5/2] atol = 1e-10
    end

    @testset "not centered" begin
        X = [Y; ones(1, 4)]

        u, R = MinimumVolumeEllipsoids._minvol(X, 1e-10)

        @test u ≈ [9 / 32, 4 / 32, 9 / 32, 10 / 32] atol = 1e-10
        @test X * Diagonal(u) * X' ≈ [31/16 13/16 1/2; 13/16 31/16 1/2; 1/2 1/2 1] atol =
            1e-10
    end
end

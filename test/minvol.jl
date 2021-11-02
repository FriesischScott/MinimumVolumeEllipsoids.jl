@testset "Minimum Volume" begin
    Y = [
        -1 -1 1 2
        1 -1 -1 2
    ]

    @testset "Centered" begin
        u, R = minvol(Y, 1e-10)

        @test u == [0.5, 0, 0, 0.5]
        @test R.L * R.L' ≈ [2.5 1.5; 1.5 2.5] atol = 1e-10
    end

    @testset "Not Centered" begin
        X = [Y; ones(1, 4)]

        u, R = minvol(X, 1e-10)

        H = Y * Diagonal(u) * Y' - Y * u * u' * Y'
        H = (H + H') / 2 # symmetry!
        H = cholesky(H)

        @test Y * u ≈ [0.5, 0.5] atol = 1e-10
        @test H.L * H.L' ≈ [1.6875 0.5625; 0.5625 1.6875]
    end
end

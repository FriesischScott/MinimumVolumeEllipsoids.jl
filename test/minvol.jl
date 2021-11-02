@testset "Minimum Volume" begin
    Y = [
        -1 -1 1 2
        1 -1 -1 2
    ]

    @testset "Centered" begin
        u, R = minvol(Y, 1e-10)
        H = inv(R)

        @test u == [0.5, 0, 0, 0.5]
        @test H ≈ [5/8 -3/8; -3/8 5/8] atol = 1e-10
    end

    @testset "Not Centered" begin
        X = [Y; ones(1, 4)]

        u, R = minvol(X, 1e-10)

        H = Y * Diagonal(u) * Y' - Y * u * u' * Y'
        H = (H + H') / 2 # symmetry!
        H = inv(H)

        @test Y * u ≈ [0.5, 0.5] atol = 1e-10
        @test H ≈ [2/3 -2/9; -2/9 2/3] atol = 1e-10
    end
end

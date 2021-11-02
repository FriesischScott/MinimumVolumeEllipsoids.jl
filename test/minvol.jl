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
end

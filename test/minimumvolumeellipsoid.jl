@testset "Minimum Volume Ellipsoid" begin
    @testset "Centered" begin
        X = [
            -1 -1 1 2
            1 -1 -1 2
        ]
        ϵ = minimum_volume_ellipsoid(X, 1e-10; centered=true)

        @test ϵ.c == [0, 0]
        @test ϵ.H ≈ [5/8 -3/8; -3/8 5/8] atol = 1e-10

        for x in eachcol(X)
            @test x in ϵ
        end
    end

    @testset "Not Centered" begin
        X = [
            -1 -1 1 2
            1 -1 -1 2
        ]
        ϵ = minimum_volume_ellipsoid(X, 1e-10)

        @test ϵ.c ≈ [0.5, 0.5] atol = 1e-10
        @test ϵ.H ≈ [2/3 -2/9; -2/9 2/3] atol = 1e-10

        for x in eachcol(X)
            @test x in ϵ
        end
    end

    @testset "KKY" begin
        X = [
            -1 -1 1 2
            1 -1 -1 2
        ]
        ϵ = minimum_volume_ellipsoid(X, 1e-10, 1)

        @test ϵ.c ≈ [0.5, 0.5] atol = 1e-10
        @test ϵ.H ≈ [2/3 -2/9; -2/9 2/3] atol = 1e-10

        for x in eachcol(X)
            @test x in ϵ
        end
    end

    @testset "Datasets" begin
        @testset "Iris" begin
            using RDatasets

            iris = dataset("datasets", "iris")

            X = permutedims([iris.SepalLength iris.SepalWidth])

            ϵ = minimum_volume_ellipsoid(X, 1e-10)

            for x in eachcol(X)
                @test x in ϵ
            end
        end

    end

    @testset "rand" begin
        X = [
            -1 -1 1 2
            1 -1 -1 2
        ]
        ϵ = minimum_volume_ellipsoid(X, 1e-10)

        U = rand(ϵ, 10^4)

        for x in eachcol(U)
            @test x in ϵ
        end
    end
end

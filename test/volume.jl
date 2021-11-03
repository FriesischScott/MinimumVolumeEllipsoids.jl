@testset "_unit_ball_volume" begin
    @test MinimumVolumeEllipsoids._unit_ball_volume(0) == 1.0
    @test MinimumVolumeEllipsoids._unit_ball_volume(1) ≈ 2.0 atol = 1e-10
    @test MinimumVolumeEllipsoids._unit_ball_volume(2) ≈ π atol = 1e-10
    @test MinimumVolumeEllipsoids._unit_ball_volume(3) ≈ 4π / 3 atol = 1e-10
    @test MinimumVolumeEllipsoids._unit_ball_volume(4) ≈ π^2 / 2 atol = 1e-10
    @test MinimumVolumeEllipsoids._unit_ball_volume(5) ≈ 8π^2 / 15 atol = 1e-10
end

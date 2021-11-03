using MinimumVolumeEllipsoids
using Plots

X = [
    -1 -1 1 2
    1 -1 -1 2
]

ϵ = minimum_volume_ellipsoid(X, centered=true)

U = rand(ϵ, 5000)

p = scatter(U[1, :], U[2, :]; legend=:none)
scatter!(p, X[1, :], X[2, :]; legend=:none)

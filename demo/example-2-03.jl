using MinimumVolumeEllipsoids
using LinearAlgebra
using Plots

X = transpose([-1 1; -1 -1; 1 -1; 2 2])

u, R = minvol(X, 1e-10)

ϵ = Ellipsoid(R, zeros(2))

U = rand(ϵ, 5000)

p = scatter(U[1, :], U[2, :]; legend=:none)
scatter!(p, X[1, :], X[2, :]; legend=:none)

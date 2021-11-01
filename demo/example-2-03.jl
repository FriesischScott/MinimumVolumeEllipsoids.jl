using MinimumVolumeEllipsoids
using LinearAlgebra
using Plots

X = transpose([-1 1; -1 -1; 1 -1; 2 2])

u, R = minvol(X, 1e-10)
U = rand(R.L, [0, 0], 1000)

p = scatter(U[1, :], U[2, :]; legend=:none, aspect_ratio=:equal)
scatter!(p, X[1, :], X[2, :]; legend=:none)

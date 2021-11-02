using MinimumVolumeEllipsoids
using LinearAlgebra
using Plots

Y = transpose([-1 1; -1 -1; 1 -1; 2 2])
X = [Y; ones(1, 4)]

u, R = minvol(X, 1e-10)

H = Y * Diagonal(u) * Y' - Y * u * u' * Y'
H = (H + H') / 2 # symmetry!
H = PDMat(inv(H))

ϵ = Ellipsoid(H, Y * u)

U = rand(ϵ, 5000)

p = scatter(U[1, :], U[2, :]; legend=:none)
scatter!(p, X[1, :], X[2, :]; legend=:none)

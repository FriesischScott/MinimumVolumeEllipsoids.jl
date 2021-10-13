using MinimumVolumeEllipsoids
using LinearAlgebra
using Plots
#= X = [-1 1; -1 -1; 1 -1; 2 2] |> transpose
u, R = minvol(X, 1e-10, 2)
H = X * Diagonal(u) * X' |> inv
U = rand(inv(H), [0.0, 0.0], 500)
p = scatter(U[1, :], U[2, :], aspect_ratio=:equal, legend=:none) =#

# X = [-1 1 1; -1 -1 1; 1 -1 1; 2 2 1] |> transpose
X = Matrix{Float64}([-2 0 1; 0 1 1; 0 -1 1; 2 0 1]) |> transpose

X[1:2, :] = [sqrt(2.0) / 2.0 sqrt(2.0) / 2.0; sqrt(2.0) / 2.0 -sqrt(2.0) / 2.0] * X[1:2, :]
u, R = minvol(X, 1e-10)
Y = X[1:2, :]
y = Y * u
# H = inv(R.U' * R.U)[1:2, 1:2]
H = Y * Diagonal(u) * Y' - Y * u * u' * Y' |> inv
H = (H' + H) / 2 # * symmetry!
# H[1, 2] *= -1
# H[2, 1] *= -1

p = scatter(Y[1, :], Y[2, :], aspect_ratio=:equal, legend=:none)

U = rand(inv(H), vec(y), 10000)
scatter!(p, U[1, :], U[2, :], aspect_ratio=:equal, legend=:none)

#= x
H

(xᵢ - x)'H(xᵢ - x) < γ
x'inv(H)x < γ =#

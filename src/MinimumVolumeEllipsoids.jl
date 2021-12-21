module MinimumVolumeEllipsoids

using Distributions
using LinearAlgebra
using PDMats
using SparseArrays
using SpecialFunctions

export Ellipsoid

export minimum_volume_ellipsoid
export rand
export volume

struct Ellipsoid
    H::AbstractPDMat
    x::AbstractVector
end

function minimum_volume_ellipsoid(
    X::AbstractMatrix,
    tol::Real=1e-7,
    KKY::Integer=0,
    maxit::Integer=100000;
    centered::Bool=false,
)
    n, m = size(X)

    if centered
        _, R, ϕ = _minvol(X, tol, KKY, maxit)

        H = inv(R.L * R.L' / ϕ)
        x = zeros(n)
    else
        Y = [X; ones(1, m)]

        u, R, ϕ = _minvol(Y, tol, KKY, maxit)

        H = inv(R.L * R.L' / ϕ)[1:n, 1:n]
        x = X * u
    end

    H = (H + H') / 2
    H = PDMat(H)

    return Ellipsoid(H, x)
end

function _minvol(X::AbstractMatrix, tol::Real=1e-7, KKY::Integer=0, maxit::Integer=100000)
    n, m = size(X)

    n100 = max(n, 100)
    n50000 = max(n, 50000)

    iter = 0

    if KKY >= 1
        u = (1 / m) * ones(m)
    else
        u = initwt(X)
    end

    # Initialize Cholesky Factor
    upos = findall(u .> 0)
    R, var = _compute_R_and_var(u, X, upos)
    ϕ = 1

    ω₊ = maximum(var)

    # Use the Harman-Pronzato test to see if columns of X can be eliminated.
    essential_indices = _harman_pronzato_elimination(ω₊, n, var, u)
    XX = X[:, essential_indices]

    # If only n columns remain, recompute u and R
    if length(essential_indices) == n
        u = (1 / n) * ones(n)
        upos = findall(u .> 1e-8)
        R, var = _compute_R_and_var(u, XX, upos)
    else
        var = var[essential_indices]
        u = u[essential_indices] / sum(u[essential_indices])
        upos = findall(u .> 1e-8)
    end

    # Find "furthest" and "closest" points
    ω₊, j = findmax(var)
    ω₋, i = findmin(var[upos])
    i = upos[i]

    while ((ω₊ > (1 + tol) * n) || (ω₋ < (1 - tol) * n)) && iter < maxit
        iter += 1

        (j, ω) = ω₊ + ω₋ > 2 * n ? (j, ω₊) : (i, ω₋)

        # compute new λ
        λ = (ω - n) / ((n - 1) * ω)
        λ = max(λ, -u[j])

        # Update u and make sure it stays nonnegative
        uold = u
        u[j] = max(u[j] + λ, 0)
        u = (1 / (1 + λ)) * u
        upos = findall(u .> 1e-8)

        Rxj = R.U' \ XX[:, j]
        x = ϕ * (R.U \ Rxj)
        ωⱼ = ϕ * (Rxj' * Rxj)

        recompute_R = abs(ωⱼ - ω) / max(1, ω) > 1e-8

        ω = ωⱼ

        # Every 50.000 iterations: Check if R should be recomputed
        if mod(iter, n50000) == 0
            upos = findall(uold .> 0)
            M = XX[:, upos] * Diagonal(uold[upos]) * XX[:, upos]'
            if norm(ϕ * M - R.L * R.L') / (ϕ * norm(M)) > 1e-8
                recompute_R = true
            end
        end

        if recompute_R
            upos = findall(u .> 0)
            R, var = _compute_R_and_var(u, XX, upos)
            ϕ = 1
        else
            xx = sqrt(abs(λ) * ϕ) * XX[:, j]
            if λ > 0
                lowrankupdate!(R, xx)
            else
                lowrankdowndate!(R, xx)
            end
            ϕ *= (1 + λ)
            mult = λ / (1 + λ * ω)
            var = (1 + λ) * (var - mult * transpose((transpose(x) * XX)) .^ 2)
        end

        ω₊, j = findmax(var)

        # Every 100 iterations: Check if more points can be eliminated
        if mod(iter, n100) == 0
            essential = _harman_pronzato_elimination(ω₊, n, var, u)
            if length(essential) < length(essential_indices)
                essential_indices = essential_indices[essential]
                XX = XX[:, essential]
                if length(essential_indices) == n
                    u = (1 / n) * ones(n, 1)
                    R, var = _compute_R_and_var(u, XX, upos)
                    ϕ = 1
                else
                    var = var[essential]
                    u = u[essential] / sum(u[essential])
                end
                ω₊, j = findmax(var)
            end
        end

        upos = findall(u .> 0)
        ω₋, i = findmin(var[upos])
        i = upos[i]
    end

    uu = zeros(m)
    uu[essential_indices] = u
    u = uu

    return u, R, ϕ
end

function _harman_pronzato_elimination(
    ω₊::Float64, n::Integer, var::AbstractVector, u::AbstractVector
)
    δn = ω₊ - n
    # TODO: Figure out when and why this happens
    if δn < 0
        δn = 0
    end
    threshold = n * (1 + δn / 2 - √(δn - δn / n + ((δn / n)^2 * n^2) / 4))
    return findall((var .> threshold) .| (u .> 1e-8))
end

function _compute_R_and_var(u::AbstractVector, X::AbstractMatrix, upos::AbstractVector)
    A = spdiagm(sqrt.(u[upos])) * X[:, upos]'
    F = qr(A)
    R = Cholesky(F.R, :U, 0)
    RX = R.U' \ X
    var = vec(sum(RX .* RX; dims=1))
    return R, var
end

"""
    initwt(X)

Obtain the initial weights `u` using the Kumar-Yildirim algorithm, taking into account that
`X` represents [X, -X].
"""
function initwt(X::AbstractMatrix)
    n, m = size(X)
    u = zeros(m)
    Q = 1.0 * I(n)
    d = Q[:, 1]

    for j in 1:n
        # compute the maximizer of | d'*x | over the columns of X.
        dX = vec(abs.(d' * X))
        _, ind = findmax(dX)
        u[ind] = 1

        j == n && break

        # update Q
        y = X[:, ind]
        z = Q' * y

        if j > 1
            z[1:(j - 1)] = zeros(j - 1, 1)
        end

        ζ = norm(z)
        zj = z[j]
        if zj < 0
            ζ = -ζ
        end
        zj += ζ
        z[j] = zj
        Q = Q - (Q * z) * ((1 / (ζ * zj)) * z')
        d = Q[:, j + 1]
    end

    u /= n
    return u
end

function Base.rand(ϵ::Ellipsoid, m::Integer)
    n = size(ϵ.H, 1)
    X = randn(n, m)
    X = X ./ kron(ones(n, 1), sqrt.(sum(X .^ 2; dims=1)))
    R = ones(n, 1) * rand(1, m) .^ (1 / n)
    sphere = R .* X
    ellipsoid = sqrt(n) * cholesky(inv(ϵ.H)).L * sphere + ϵ.x .* ones(1, m)
    return ellipsoid
end

include("volume.jl")

end # module

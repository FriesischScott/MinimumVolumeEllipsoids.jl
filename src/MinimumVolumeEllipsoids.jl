module MinimumVolumeEllipsoids

using Distributions
using LinearAlgebra
using PDMats

export Ellipsoid

export minimum_volume_ellipsoid
export rand
export volume

struct Ellipsoid
    H::PDMat
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
        u, R = _minvol(X, tol, KKY, maxit)

        H = PDMat(inv(R))
        return Ellipsoid(H, zeros(n))
    else
        Y = [X; ones(1, m)]

        u, R = _minvol(Y, tol, KKY, maxit)

        H = X * Diagonal(u) * X' - X * u * u' * X'
        H = (H + H') / 2 # symmetry!
        H = PDMat(inv(H))

        return Ellipsoid(H, X * u)
    end
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

    act = fill(true, m)
    XX = copy(X)

    # Use the Harman-Pronzato test to see if columns of X can be eliminated.
    δn = ω₊ - n
    thresh = n * (1 + δn / 2 - √(δn - δn / n + ((δn / n)^2 * n^2) / 4))
    essential = (var .> thresh) .| (u .> 1e-8)
    act = act .& essential
    XX = X[:, essential]

    # If only n columns remain, recompute u and R
    if sum(act) == n
        u = (1 / n) * ones(n)
        upos = findall(u .> 1e-8)
        R, var = _compute_R_and_var(u, XX, upos)
    else
        var = var[essential]
        u = u[essential] / sum(u[essential])
        upos = findall(u .> 1e-8)
    end

    # Find "furthest" and "closest" points
    ω₊, j = findmax(var)
    ω₋, i = findmin(var[upos])
    i = upos[i]

    ϵ₊ = (ω₊ - n) / n
    ϵ₋ = (ω₋ - n) / n

    while max(ϵ₊, ϵ₋) > tol && iter < maxit
        iter += 1

        (j, ω) = ϵ₊ > ϵ₋ ? (j, ω₊) : (i, ω₋)

        # compute new λ
        λ = (ω - n) / ((n - 1) * ω)
        λ = max(λ, -u[j])

        # Update u and make sure it stays nonnegative
        uold = u
        u[j] = max(u[j] + λ, 0)
        u = (1 / (1 + λ)) * u
        upos = findall(u .> 1e-8)

        z = inv(R.L) * XX[:, j]
        ωⱼ = ϕ * transpose(z) * z
        x = ϕ * transpose(inv(R.L)) * z

        recompute_R = abs(ωⱼ - ω) / max(1, ω) > 1e-8

        ω = ωⱼ

        # Every 50.000 iterations: Check if R should be recomputed
        if mod(iter, n50000) == 0
            upos = findall(uold .> 0)
            M = XX[:, upos] * Diagonal(uold[upos]) * XX[:, upos]'
            if norm(ϕ * M - R) / (ϕ * norm(M)) > 1e-8
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
            δn = ω₊ - n
            thresh = n * (1 + δn / 2 - √(δn - δn / n + ((δn / n)^2 * n^2) / 4))
            essential = (var .> thresh) .| (u .> 1e-8)
            if sum(essential) < sum(act)
                act = act .& essential
                XX = XX[:, essential]
                if sum(act) == n
                    u = (1 / n) * ones(n, 1)
                    uold = u
                    upos = findall(u .> 1e-8)
                    R, var = _compute_R_and_var(u, XX, upos)
                    ϕ = 1
                else
                    var = var[essential]
                    u = u[essential] / sum(u[essential])
                    uold = uold[essential] / sum(uold[essential])
                    upos = findall(u .> 1e-8)
                end
                ω₊, j = findmax(var)
            end
        end

        upos = findall(u .> 0)
        ω₋, i = findmin(var[upos])
        i = upos[i]

        ϵ₊ = (ω₊ - n) / n
        ϵ₋ = (ω₋ - n) / n
    end

    uu = zeros(m)
    uu[act] = u
    u = uu

    return u, R
end

function _compute_R_and_var(u::AbstractVector, X::AbstractMatrix, upos::AbstractVector)
    A = Diagonal(sqrt.(u[upos])) * transpose(X[:, upos])
    _, R = qr(A)
    R = Cholesky(R, :U, 0)
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

function volume(ϵ::Ellipsoid)
    n = size(ϵ.H, 1)
    return (n^(n / 2) * _unit_ball_volume(n)) / sqrt(det(ϵ.H))
end

function _unit_ball_volume(n::Integer)
    return π^(n / 2) / gamma(n / 2 + 1)
end

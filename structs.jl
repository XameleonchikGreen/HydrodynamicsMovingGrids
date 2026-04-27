# primitive variables struct
struct PrimVars
    ρ::AbstractFloat
    u::AbstractFloat
    p::AbstractFloat
end

function soundSpeed(u::PrimVars, γ::BigFloat)::BigFloat
    return sqrt(γ * u.p / u.ρ)
end

function calculateEnergy(u::PrimVars, γ::BigFloat)::BigFloat
    return u.p / (γ - 1) / u.ρ
end

PrimVars((ρ, u, p)) = PrimVars(ρ, u, p)

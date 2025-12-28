# primitive variables struct
struct PrimVars
    ρ::AbstractFloat
    u::AbstractFloat
    p::AbstractFloat
end

function sound_speed(u::PrimVars, γ::BigFloat)::BigFloat
    return sqrt(γ * u.p / u.ρ)
end

PrimVars((ρ, u, p)) = PrimVars(ρ, u, p)

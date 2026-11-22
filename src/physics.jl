function sound_speed(u::PrimVars, γ::BigFloat)::BigFloat
    return sqrt(γ * u.p / u.ρ)
end

function calculate_energy(u::PrimVars, γ::BigFloat)::BigFloat
    return u.p / (γ - 1) / u.ρ
end
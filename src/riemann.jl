function gₖ(p::BigFloat, pₖ::BigFloat, ρₖ::BigFloat, γ::BigFloat)
    Aₖ = 2 / (γ + 1)ρₖ
    Bₖ = (γ - 1) / (γ + 1) * pₖ
    return sqrt(Aₖ / (p + Bₖ))
end

function fₖ(p::BigFloat, Wₖ::PrimVars, γ::BigFloat)
    ρₖ = Wₖ.ρ
    pₖ = Wₖ.p
    if p > pₖ
        return (p - pₖ) * gₖ(p, pₖ, ρₖ, γ)
    else
        aₖ = sound_speed(Wₖ, γ)
        return 2aₖ / (γ - 1) * ((p / pₖ) ^ ((γ - 1) / 2γ) - 1)
    end
end

function fꜝₖ(p::BigFloat, Wₖ::PrimVars, γ::BigFloat)
    ρₖ = Wₖ.ρ
    pₖ = Wₖ.p
    if p > pₖ
        Bₖ = (γ - 1) / (γ + 1) * pₖ
        return (1 - (p - pₖ) / (2 * (Bₖ + p))) * gₖ(p, pₖ, ρₖ, γ)
    else
        aₖ = sound_speed(Wₖ, γ)
        return 1 / (ρₖ * aₖ) * (p / pₖ) ^ (-(γ + 1) / 2γ)
    end
end    

f(p, Wₗ, Wᵣ, γ) = fₖ(p, Wₗ, γ) + fₖ(p, Wᵣ, γ) + Wᵣ.u - Wₗ.u
fꜝ(p, Wₗ, Wᵣ, γ) = fꜝₖ(p, Wₗ, γ) + fꜝₖ(p, Wᵣ, γ)

function newton_step(p, Wₗ, Wᵣ, γ)
    return p - f(p, Wₗ, Wᵣ, γ) / fꜝ(p, Wₗ, Wᵣ, γ)
end

function newton_method(p₀, Wₗ, Wᵣ, γ, ε)
    p₁ = newton_step(p₀, Wₗ, Wᵣ, γ)
    while 2abs(p₁ - p₀) / (p₁ + p₀) > ε
        p₀ = p₁
        p₁ = newton_step(p₀, Wₗ, Wᵣ, γ)
    end
    return p₁
end

function Wₗᶠᵃⁿ(Wₗ::PrimVars, γ::BigFloat, S::BigFloat)::PrimVars
    ρₗ = Wₗ.ρ
    uₗ = Wₗ.u
    pₗ = Wₗ.p
    aₗ = sound_speed(Wₗ, γ)
    ρ = ρₗ * (2 / (γ + 1) + (γ - 1) / (γ + 1) / aₗ * (uₗ - S)) ^ (2 / (γ - 1))
    u = 2 / (γ + 1) * (aₗ + (γ - 1) / 2 * uₗ + S)
    p = pₗ * (2 / (γ + 1) + (γ - 1) / (γ + 1) / aₗ * (uₗ - S)) ^ (2γ / (γ - 1))
    return PrimVars(ρ, u, p)
end

function Wᵣᶠᵃⁿ(Wᵣ::PrimVars, γ::BigFloat, S::BigFloat)::PrimVars
    ρᵣ = Wᵣ.ρ
    uᵣ = Wᵣ.u
    pᵣ = Wᵣ.p
    aᵣ = sound_speed(Wᵣ, γ)
    ρ = ρᵣ * (2 / (γ + 1) - (γ - 1) / (γ + 1) / aᵣ * (uᵣ - S)) ^ (2 / (γ - 1))
    u = 2 / (γ + 1) * (-aᵣ + (γ - 1) / 2 * uᵣ + S)
    p = pᵣ * (2 / (γ + 1) - (γ - 1) / (γ + 1) / aᵣ * (uᵣ - S)) ^ (2γ / (γ - 1))
    return PrimVars(ρ, u, p)
end

function W_non_vacuum(Wₗ::PrimVars,
                      Wᵣ::PrimVars,
                      W⁺::Tuple{BigFloat, BigFloat}, 
                      γ::BigFloat, 
                      S::BigFloat,
                     )::PrimVars
    ρₗ = Wₗ.ρ
    uₗ = Wₗ.u
    pₗ = Wₗ.p
    aₗ = sound_speed(Wₗ, γ)

    ρᵣ = Wᵣ.ρ
    uᵣ = Wᵣ.u
    pᵣ = Wᵣ.p
    aᵣ = sound_speed(Wᵣ, γ)

    u⁺ = W⁺[1]
    p⁺ = W⁺[2]

    if S < u⁺
        if p⁺ > pₗ # left shock wave
            Sₗ = uₗ - aₗ * sqrt(p⁺ * (γ + 1) / 2γ / pₗ + (γ - 1) / 2γ)
            if Sₗ ≤ S
                ρ⁺ₗ = ρₗ * (p⁺ / pₗ + (γ - 1) / (γ + 1)) / ((γ - 1) * p⁺ / (γ + 1) / pₗ + 1)
                return PrimVars(ρ⁺ₗ, u⁺, p⁺)
            else
                return Wₗ
            end
        else # left rarefaction
            a⁺ₗ = aₗ * (p⁺ / pₗ) ^ ((γ - 1) / 2γ)
            Sₕₗ = uₗ - aₗ
            Sₜₗ = u⁺ - a⁺ₗ
            if S < Sₕₗ
                return Wₗ
            elseif S ≤ Sₜₗ
                return Wₗᶠᵃⁿ(Wₗ, γ, S)
            else
                ρ⁺ₗ = ρₗ * (p⁺ / pₗ) ^ (1 / γ)
                return PrimVars(ρ⁺ₗ, u⁺, p⁺)
            end
        end
    else
        if p⁺ > pᵣ # right shock wave
            Sᵣ = uᵣ + aᵣ * sqrt(p⁺ * (γ + 1) / 2γ / pᵣ + (γ - 1) / 2γ)
            if S ≤ Sᵣ
                ρ⁺ᵣ = ρᵣ * (p⁺ / pᵣ + (γ - 1) / (γ + 1)) / ((γ - 1) * p⁺ / (γ + 1) / pᵣ + 1)
                return PrimVars(ρ⁺ᵣ, u⁺, p⁺)
            else
                return Wᵣ
            end
        else # right rarefaction
            a⁺ᵣ = aᵣ * (p⁺ / pᵣ) ^ ((γ - 1) / 2γ)
            Sₕᵣ = uᵣ + aᵣ
            Sₜᵣ = u⁺ + a⁺ᵣ
            if S < Sₜᵣ
                ρ⁺ᵣ = ρᵣ * (p⁺ / pᵣ) ^ (1 / γ)
                return PrimVars(ρ⁺ᵣ, u⁺, p⁺)
            elseif S ≤ Sₕᵣ
                return Wᵣᶠᵃⁿ(Wᵣ, γ, S)
            else
                return Wᵣ
            end
        end
    end
end

function W_vacuum(Wₗ::PrimVars,
                  Wᵣ::PrimVars, 
                  γ::BigFloat,
                  ε::BigFloat, 
                  S::BigFloat,
                 )::PrimVars

    uₗ = Wₗ.u
    aₗ = sound_speed(Wₗ, γ)

    uᵣ = Wᵣ.u
    aᵣ = sound_speed(Wᵣ, γ)

    if abs(ρᵣ) < ε # right vacuum
        u₀ = uₗ + 2aₗ / (γ - 1)
        if S ≤ uₗ - aₗ
            return Wₗ
        elseif S < u₀
            return Wₗᶠᵃⁿ(Wₗ, γ, S)
        else
            return PrimVars(BigFloat(0), u₀, BigFloat(0))
        end
    elseif abs(ρₗ) < ε # left vacuum
        u₀ = uᵣ - 2aᵣ / (γ - 1)
        if S ≤ u₀
            return PrimVars(BigFloat(0), u₀, BigFloat(0))
        elseif S < uᵣ + aᵣ
            return Wᵣᶠᵃⁿ(Wᵣ, γ, S)
        else
            return Wᵣ
        end
    else # vacuum generated
        S⁺ₗ = uₗ + 2aₗ / (γ - 1)
        S⁺ᵣ = uᵣ - 2aᵣ / (γ - 1)
        if S < S⁺ₗ
            if S ≤ uₗ - aₗ
                return Wₗ
            else
                return Wₗᶠᵃⁿ(Wₗ, γ, S)
            end
        elseif S⁺ₗ ≤ S ≤ S⁺ᵣ
            return PrimVars(BigFloat(0), S⁺ₗ, BigFloat(0))
        else
            if S < uᵣ + aᵣ
                return Wᵣᶠᵃⁿ(Wᵣ, γ, S)
            else
                return Wᵣ
            end
        end
    end
end

function W⁺(Wₗ::PrimVars,
            Wᵣ::PrimVars, 
            γ::BigFloat, 
            ε::BigFloat;
            p⁺::Union{BigFloat, Nothing} = nothing,
            )::Function
    
    final_p⁺ = isnothing(p⁺) ? find_p⁺(Wₗ, Wᵣ, γ, ε) : p⁺
    uₗ = Wₗ.u
    uᵣ = Wᵣ.u
    if isnothing(final_p⁺)
        return S -> W_vacuum(Wₗ, Wᵣ, γ, ε, S)
    else
        u⁺ = 1 / 2 * (uₗ + uᵣ) + 1 / 2 * (fₖ(final_p⁺, Wᵣ, γ) - fₖ(final_p⁺, Wₗ, γ))
        W⁺ = (u⁺, final_p⁺)
        return S -> W_non_vacuum(Wₗ, Wᵣ, W⁺, γ, S)
    end
end

function find_p⁺(
                Wₗ::PrimVars,
                Wᵣ::PrimVars, 
                γ::BigFloat, 
                ε::BigFloat,
              )::BigFloat
    ρₗ = Wₗ.ρ
    uₗ = Wₗ.u
    pₗ = Wₗ.p
    aₗ = sound_speed(Wₗ, γ)

    ρᵣ = Wᵣ.ρ
    uᵣ = Wᵣ.u
    pᵣ = Wᵣ.p
    aᵣ = sound_speed(Wᵣ, γ)

    if 2aₗ / (γ - 1) + 2aᵣ / (γ - 1) > uᵣ - uₗ # non-vacuum
        pₚᵥ = 1 / 2 * (pₗ + pᵣ) - 1 / 8 * (uᵣ - uₗ) * (ρₗ + ρᵣ) * (aₗ + aᵣ)
        p̂ = max(ε, pₚᵥ)
        gₗ = gₖ(p̂, pₗ, ρₗ, γ)
        gᵣ = gₖ(p̂, pᵣ, ρᵣ, γ)
        pₜₛ = (gₗ * pₗ + gᵣ * pᵣ - uᵣ + uₗ) / (gₗ + gᵣ)
        p₀ = max(ε, pₜₛ)

        p⁺ = newton_method(p₀, Wₗ, Wᵣ, γ, ε)
        return p⁺
    else # vacuum
        return nothing
    end
end

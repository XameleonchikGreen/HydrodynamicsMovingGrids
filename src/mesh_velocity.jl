function find_Q_by_pressure(
                Wₗ::PrimVars,
                Wᵣ::PrimVars, 
                γ::BigFloat, 
                p⁺::BigFloat,
                ε₁::BigFloat,
                ε₂::BigFloat
              )::BigFloat
    ρₗ = Wₗ.ρ
    uₗ = Wₗ.u
    pₗ = Wₗ.p
    aₗ = sound_speed(Wₗ, γ)

    ρᵣ = Wᵣ.ρ
    uᵣ = Wᵣ.u
    pᵣ = Wᵣ.p
    aᵣ = sound_speed(Wᵣ, γ)

    u⁺ = 1 / 2 * (uₗ + uᵣ) + 1 / 2 * (fₖ(p⁺, Wᵣ, γ) - fₖ(p⁺, Wₗ, γ))

    if abs(pₗ - pᵣ) < ε₁ && abs(ρₗ - ρᵣ) > ε₂
        return - u⁺
    elseif abs(pₗ - pᵣ) < ε₁
        return BigFloat("0")
    else
        if abs(pₗ - p⁺) > abs(pᵣ - p⁺) # left wave dominates right
            if p⁺ > pₗ                      # left wave is shock
                return - (uₗ - aₗ * sqrt(p⁺ * (γ + 1) / 2γ / pₗ + (γ - 1) / 2γ))
            else                            # left wave is rarefaction
                return - ((uₗ + uᵣ) / 2 - (aₗ + aᵣ) / 2)
            end
        elseif abs(pₗ - p⁺) < abs(pᵣ - p⁺) # right wave dominates left
            if p⁺ > pᵣ                      # right wave is shock
                return - (uᵣ + aᵣ * sqrt(p⁺ * (γ + 1) / 2γ / pᵣ + (γ - 1) / 2γ))
            else                            # right wave is rarefaction
                return - ((uₗ + uᵣ) / 2 + (aₗ + aᵣ) / 2)
            end
        else
            return BigFloat("0")
        end
    end
end

function find_Q_by_density(
                Wₗ::PrimVars,
                Wᵣ::PrimVars, 
                γ::BigFloat, 
                p⁺::BigFloat,
                ε₁::BigFloat,
                ε₂::BigFloat
              )::BigFloat
    ρₗ = Wₗ.ρ
    uₗ = Wₗ.u
    pₗ = Wₗ.p
    aₗ = sound_speed(Wₗ, γ)

    ρᵣ = Wᵣ.ρ
    uᵣ = Wᵣ.u
    pᵣ = Wᵣ.p
    aᵣ = sound_speed(Wᵣ, γ)

    u⁺ = 1 / 2 * (uₗ + uᵣ) + 1 / 2 * (fₖ(p⁺, Wᵣ, γ) - fₖ(p⁺, Wₗ, γ))

    if p⁺ > pₗ # left shock wave
        ρ⁺ₗ = ρₗ * (p⁺ / pₗ + (γ - 1) / (γ + 1)) / ((γ - 1) * p⁺ / (γ + 1) / pₗ + 1)
    else # left rarefaction
        ρ⁺ₗ = ρₗ * (p⁺ / pₗ) ^ (1 / γ)
    end
    if p⁺ > pᵣ # right shock wave
        ρ⁺ᵣ = ρᵣ * (p⁺ / pᵣ + (γ - 1) / (γ + 1)) / ((γ - 1) * p⁺ / (γ + 1) / pᵣ + 1)
    else # right rarefaction
        ρ⁺ᵣ = ρᵣ * (p⁺ / pᵣ) ^ (1 / γ)
    end

    Δρₗ = abs(ρₗ - ρ⁺ₗ)
    Δρ = abs(ρ⁺ₗ - ρ⁺ᵣ)
    Δρᵣ = abs(ρᵣ - ρ⁺ᵣ)

    # if abs(pₗ - pᵣ) < ε₁ && abs(ρₗ - ρᵣ) > ε₂
    #     return - u⁺
    # end

    if max(Δρₗ, Δρ, Δρᵣ) < 1e-5
        return BigFloat(0)
    end

    if Δρₗ > Δρ && Δρₗ > Δρᵣ            # left wave dominates                        
        if p⁺ > pₗ                      # left wave is shock
            return - (uₗ - aₗ * sqrt(p⁺ * (γ + 1) / 2γ / pₗ + (γ - 1) / 2γ))
        else                            # left wave is rarefaction
            return - ((uₗ + uᵣ) / 2 - (aₗ + aᵣ) / 2)
        end
    elseif Δρᵣ > Δρ && Δρᵣ > Δρₗ        # right wave dominates
        if p⁺ > pᵣ                      # right wave is shock
                return - (uᵣ + aᵣ * sqrt(p⁺ * (γ + 1) / 2γ / pᵣ + (γ - 1) / 2γ))
        else                            # right wave is rarefaction
            return - ((uₗ + uᵣ) / 2 + (aₗ + aᵣ) / 2)
        end        
    elseif Δρ > Δρₗ && Δρ > Δρᵣ
        return - u⁺
    else
        return BigFloat(0)
    end
end

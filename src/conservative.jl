function ρ_u_p_E(W::PrimVars, γ::BigFloat)
    ρ = W.ρ
    u = W.u
    p = W.p
    E = p / (γ - 1) + 1 / 2 * ρ * u ^ 2
    return ρ, u, p, E
end

function ρ_u_p_E(states::AbstractVector{<:PrimVars}, γ::BigFloat)
    ρ = getfield.(states, :ρ)
    u = getfield.(states, :u)
    p = getfield.(states, :p)
    E = @. p / (γ - 1) + 1 / 2 * ρ * u ^ 2
    return ρ, u, p, E
end

function compute_static_update(
    ρᵏ::AbstractVector, uᵏ::AbstractVector, Eᵏ::AbstractVector,
    ρ::AbstractVector, u::AbstractVector, p::AbstractVector, E::AbstractVector,
    Δt::BigFloat, Δy::AbstractVector, γ::BigFloat
)
    ρᵏ⁺¹ = @. ρᵏ - Δt / Δy * (ρ[2:end] * u[2:end] -
                            ρ[1:end-1] * u[1:end-1])
    ρuᵏ⁺¹ = @. ρᵏ * uᵏ - Δt / Δy * (ρ[2:end] * u[2:end] ^ 2 + p[2:end] -
                                    ρ[1:end-1] * u[1:end-1] ^ 2 - p[1:end-1])
    Eᵏ⁺¹ = @. Eᵏ - Δt / Δy * ((E[2:end] + p[2:end]) * u[2:end] -
                            (E[1:end-1] + p[1:end-1]) * u[1:end-1])
    uᵏ⁺¹ = @. ρuᵏ⁺¹ / ρᵏ⁺¹
    pᵏ⁺¹ = @. (γ - 1) * (Eᵏ⁺¹ - 1 / 2 * ρᵏ⁺¹ * uᵏ⁺¹ ^ 2)

    return ρᵏ⁺¹, uᵏ⁺¹, pᵏ⁺¹
end

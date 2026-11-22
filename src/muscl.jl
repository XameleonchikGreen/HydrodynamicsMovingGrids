function muscl_predictor_interface_states(
    Ws::AbstractVector{<:PrimVars},
    γ::BigFloat,
    Δt::BigFloat,
    Δy₁::BigFloat
)
    ρᵏ, uᵏ, pᵏ, _ = ρ_u_p_E(Ws, γ)

    ρˡ, ρʳ = muscl_reconstruct(ρᵏ)
    uˡ, uʳ = muscl_reconstruct(uᵏ)
    pˡ, pʳ = muscl_reconstruct(pᵏ)

    Eˡ = @. pˡ / (γ - 1) + 1 / 2 * ρˡ * uˡ^2
    Eʳ = @. pʳ / (γ - 1) + 1 / 2 * ρʳ * uʳ^2

    half_dt_over_dy = Δt / (2 * Δy₁)

    mass_flux_diff = @. ρˡ * uˡ - ρʳ * uʳ

    momentum_flux_diff = @. ρˡ * uˡ^2 + pˡ - ρʳ * uʳ^2 - pʳ

    energy_flux_diff = @. (Eˡ + pˡ) * uˡ - (Eʳ + pʳ) * uʳ

    ρuˡ_new = @. ρˡ * uˡ + half_dt_over_dy * momentum_flux_diff
    ρuʳ_new = @. ρʳ * uʳ + half_dt_over_dy * momentum_flux_diff

    ρˡ_new = @. ρˡ + half_dt_over_dy * mass_flux_diff
    ρʳ_new = @. ρʳ + half_dt_over_dy * mass_flux_diff

    Eˡ_new = @. Eˡ + half_dt_over_dy * energy_flux_diff
    Eʳ_new = @. Eʳ + half_dt_over_dy * energy_flux_diff

    uˡ_new = @. ρuˡ_new / ρˡ_new
    uʳ_new = @. ρuʳ_new / ρʳ_new

    pˡ_new = @. (γ - 1) * (Eˡ_new - 1 / 2 * ρˡ_new * uˡ_new^2)
    pʳ_new = @. (γ - 1) * (Eʳ_new - 1 / 2 * ρʳ_new * uʳ_new^2)

    Fˡ = PrimVars.(ρˡ_new[2:end], uˡ_new[2:end], pˡ_new[2:end])
    Fʳ = PrimVars.(ρʳ_new[1:end-1], uʳ_new[1:end-1], pʳ_new[1:end-1])

    return Fʳ, Fˡ
end

function MUSCL_step(
        solution::Solution{Static},
        config::ProblemConfig,
        Δt::BigFloat
    )
    Ws = solution.Ws
    phys = config.physics
    γ = phys.γ
    grid = solution.grid

    Fʳ, Fˡ = muscl_predictor_interface_states(
        Ws,
        γ,
        Δt,
        grid.Δy[1]
    )

    p⁺ = find_p⁺.(Fʳ, Fˡ, γ, phys.ε_newton)
    Wᶠˡᵃˣ = [
        W⁺(Ws1, Ws2, γ, phys.ε_newton; p⁺=p1)(BigFloat(0))
        for (Ws1, Ws2, p1) in zip(Fʳ, Fˡ, p⁺)
    ]

    ρᵏ, uᵏ, _, Eᵏ = ρ_u_p_E(Ws[2:end-1], γ)
    ρ, u, p, E = ρ_u_p_E(Wᶠˡᵃˣ, γ)

    ρᵏ⁺¹, uᵏ⁺¹, pᵏ⁺¹ = compute_static_update(
        ρᵏ, uᵏ, Eᵏ,
        ρ, u, p, E,
        Δt, grid.Δy, γ
    )

    solution.t += Δt
    solution.Ws[2:end-1] .= PrimVars.(ρᵏ⁺¹, uᵏ⁺¹, pᵏ⁺¹)
    apply_bc!(solution, config)
    
    return solution
end

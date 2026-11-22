function godunov_step(
        solution::Solution{Static},
        config::ProblemConfig,
        Δt::BigFloat
    )
    Ws = solution.Ws
    phys = config.physics
    γ = phys.γ
    grid = solution.grid

    p⁺ = find_p⁺.(Ws[1:end-1], Ws[2:end], γ, phys.ε_newton)
    Wᶠˡᵃˣ = [
        W⁺(Ws1, Ws2, γ, phys.ε_newton; p⁺=p1)(BigFloat(0))
        for (Ws1, Ws2, p1) in zip(Ws[1:end-1], Ws[2:end], p⁺)
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

function godunov_step(
        solution::Solution{Moving},
        config::ProblemConfig,
        Δt::BigFloat
    )
    Ws = solution.Ws
    phys = config.physics
    γ = phys.γ
    grid = solution.grid

    p⁺ = find_p⁺.(Ws[1:end-1], Ws[2:end], γ, phys.ε_newton)
    find_Q = config.mesh_velocity_alg == :press ? find_Q_by_pressure : find_Q_by_density
    grid.Qs[2:end-1] = find_Q.(Ws[2:end-2], Ws[3:end-1], γ, p⁺[2:end-1], 
                        phys.ε_pressure, phys.ε_density) .+ 
                        phys.α ./ max.(grid.Δy[2:end], grid.Δy[1:end-1]) .* 
                        (1 ./ grid.ψs[2:end] .- 1 ./ grid.ψs[1:end-1])
    grid.nodes[2:end-1] = grid.nodes[2:end-1] .- grid.Qs[2:end-1] .* Δt
    if any(grid.nodes[2:end] .< grid.nodes[1:end-1])
        @warn "Grid degenerated"
    end
    
    ψs_n = @. grid.ψs - Δt / grid.Δy * (grid.Qs[2:end] - grid.Qs[1:end-1])

    Wᶠˡᵃˣ = [W⁺(Ws1, Ws2, γ, phys.ε_newton; p⁺=p1)(-Qs) for (Ws1, Ws2, Qs, p1) in 
        zip(Ws[1:end-1], Ws[2:end], grid.Qs, p⁺)]
    ρᵏ, uᵏ, _, Eᵏ = ρ_u_p_E(Ws[2:end-1], γ)
    ρ, u, p, E = ρ_u_p_E(Wᶠˡᵃˣ, γ)

    ρᵏ⁺¹ = @. (
                grid.ψs * ρᵏ - Δt / grid.Δy * (
                    grid.Qs[2:end] * ρ[2:end] + ρ[2:end] * u[2:end] -
                    grid.Qs[1:end-1] * ρ[1:end-1] - ρ[1:end-1] * u[1:end-1]
                )
              ) / ψs_n
    ρuᵏ⁺¹ = @. (
                grid.ψs * ρᵏ * uᵏ - Δt / grid.Δy * (
                    grid.Qs[2:end] * ρ[2:end] * u[2:end] + ρ[2:end] * u[2:end] ^ 2 + p[2:end] -
                    grid.Qs[1:end-1] * ρ[1:end-1] * u[1:end-1] - ρ[1:end-1] * u[1:end-1] ^ 2 - p[1:end-1]
                )
               ) / ψs_n
    Eᵏ⁺¹ = @. (
                grid.ψs * Eᵏ - Δt / grid.Δy * (
                    grid.Qs[2:end] * E[2:end] + (E[2:end] + p[2:end]) * u[2:end] -
                    grid.Qs[1:end-1] * E[1:end-1] - (E[1:end-1] + p[1:end-1]) * u[1:end-1]
                )
              ) / ψs_n
    uᵏ⁺¹ = @. ρuᵏ⁺¹ / ρᵏ⁺¹
    pᵏ⁺¹ = @. (γ - 1) * (Eᵏ⁺¹ - 1 / 2 * ρᵏ⁺¹ * uᵏ⁺¹ ^ 2)

    grid.ψs = ψs_n
    solution.t += Δt
    solution.Ws[2:end-1] .= PrimVars.(ρᵏ⁺¹, uᵏ⁺¹, pᵏ⁺¹)
    apply_bc!(solution, config)
    
    return solution
end

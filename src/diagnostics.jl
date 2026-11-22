function spatial_error(
        solution::Solution{T},
        exact_solution::Function,
        γ::BigFloat,
        norm::ErrorNorm
    ) where T
    if iszero(solution.t)
        z = zero(BigFloat)
        return (ρ = z, u = z, p = z, energy = z)
    end

    grid = solution.grid
    Ws = solution.Ws
    centers = cell_centers(grid)
    
    err_ρ = BigFloat(0)
    err_u = BigFloat(0)
    err_p = BigFloat(0)
    err_energy = BigFloat(0)
    
    for i in eachindex(centers)
        exact = exact_solution(centers[i] / solution.t)
        
        e_ρ = abs(Ws[i + 1].ρ - exact.ρ)
        e_u = abs(Ws[i + 1].u - exact.u)
        e_p = abs(Ws[i + 1].p - exact.p)
        e_energy = abs(calculate_energy(Ws[i + 1], γ) - calculate_energy(exact, γ))
        
        Δx = grid.nodes[i + 1] - grid.nodes[i]
        if norm == L1
            err_ρ += e_ρ * Δx
            err_u += e_u * Δx
            err_p += e_p * Δx
            err_energy += e_energy * Δx
        elseif norm == L2
            err_ρ += e_ρ^2 * Δx
            err_u += e_u^2 * Δx
            err_p += e_p^2 * Δx
            err_energy += e_energy^2 * Δx
        elseif norm == Linf
            err_ρ = max(err_ρ, e_ρ)
            err_u = max(err_u, e_u)
            err_p = max(err_p, e_p)
            err_energy = max(err_energy, e_energy)
        end
    end
    
    if norm == L2
        err_ρ = sqrt(err_ρ)
        err_u = sqrt(err_u)
        err_p = sqrt(err_p)
        err_energy = sqrt(err_energy)
    end
    
    return (ρ = err_ρ, u = err_u, p = err_p, energy = err_energy)
end

function accumulate_error!(
    acc::ErrorAccumulator,
    err,
    Δt_prev::BigFloat,
    Δt::BigFloat
)
    acc.ρ += err.ρ * (Δt_prev + Δt) / 2
    acc.u += err.u * (Δt_prev + Δt) / 2
    acc.p += err.p * (Δt_prev + Δt) / 2
    acc.energy += err.energy * (Δt_prev + Δt) / 2
end

function make_error_hooks!(
    hooks::HookSet,
    exact_solution::Function,
    norm::ErrorNorm
)
    acc = ErrorAccumulator()

    # push!(hooks.after_birth, function (sol, cfg, Δt_prev, Δt)
    #     err = spatial_error(sol, exact_solution, cfg.physics.γ, norm)
    #     accumulate_error!(acc, err, Δt_prev, Δt)
    #     return nothing
    # end)

    push!(hooks.after_step, function (sol, cfg, Δt_prev)
        err = spatial_error(sol, exact_solution, cfg.physics.γ, norm)
        accumulate_error!(acc, err, Δt_prev, Δt_prev)
        return nothing
    end)

    return acc
end

function make_trajectory_hooks!(
    hooks::HookSet,
    solution::Solution,
    snapshot_step::BigFloat
)
    snapshots = GridSnapshot[]

    grid = solution.grid

    push!(
        snapshots,
        GridSnapshot(grid.nodes, grid.ids, grid.max_id, solution.t)
    )

    next_snapshot_time = Ref(snapshot_step)

    function push_snapshot_unique!(sol::Solution)
        if isempty(snapshots) || snapshots[end].t < sol.t
            grid = sol.grid
            push!(
                snapshots,
                GridSnapshot(grid.nodes, grid.ids, grid.max_id, sol.t)
            )
        end
        return nothing
    end

    push!(hooks.after_birth, function (sol, cfg, Δt_prev, Δt)
        if sol.t ≥ next_snapshot_time[]
            push_snapshot_unique!(sol)

            while next_snapshot_time[] ≤ sol.t
                next_snapshot_time[] += snapshot_step
            end
        end

        return nothing
    end)

    push!(hooks.after_loop, function (sol, cfg)
        push_snapshot_unique!(sol)

        return nothing
    end)

    return snapshots
end

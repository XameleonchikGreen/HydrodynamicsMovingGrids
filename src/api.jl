function simulate(
    grid::Grid{T},
    config::ProblemConfig,
    method::NumericalMethod;
    birth_nodes::Bool = false,
    exact_solution::Union{Nothing,Function} = nothing,
    error_norm::Union{Nothing,ErrorNorm} = nothing,
    snapshot_step::Union{Nothing,BigFloat} = nothing
) where T
    if (exact_solution === nothing) != (error_norm === nothing)
        throw(ArgumentError(
            "exact_solution and error_norm must be provided together."
        ))
    end

    if snapshot_step !== nothing && snapshot_step ≤ 0
        throw(ArgumentError(
            "snapshot_step must be positive if snapshots are requested."
        ))
    end

    solution = Solution{T}(grid, config.ic, BigFloat(0))

    hooks = HookSet()

    err_acc = nothing

    if exact_solution !== nothing && error_norm !== nothing
        err_acc = make_error_hooks!(
            hooks, exact_solution, error_norm
        )
    end

    snaps = nothing

    if snapshot_step !== nothing
        snaps = make_trajectory_hooks!(
            hooks, solution, snapshot_step
        )
    end

    solve!(
        solution,
        config,
        method;
        birth_nodes = birth_nodes,
        hooks = hooks
    )

    return (
        solution=solution,
        error=err_acc,
        snapshots=snaps
    )
end


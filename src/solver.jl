function create_simulation_progress(
    grid::Grid{T},
    method::NumericalMethod
) where T <: GridType
    return ProgressUnknown(
        dt=0.2,
        desc=method_name(method) * ", " * 
            grid_name(T) * ", N = $(length(grid.nodes) - 1), iteration:",
        output=stdout,
        showspeed=true
    )
end

function step!(
    ::GodunovScheme,
    solution::Solution,
    config::ProblemConfig,
    Δt::BigFloat
)
    return godunov_step(solution, config, Δt)
end

function step!(
    ::MUSCLScheme,
    solution::Solution,
    config::ProblemConfig,
    Δt::BigFloat
)
    return MUSCL_step(solution, config, Δt)
end

function solve!(
    solution::Solution,
    config::ProblemConfig,
    method::NumericalMethod;
    birth_nodes::Bool=false,
    hooks::HookSet = HookSet()
)
    add_ghost_cells!(solution, config)

    Δx = min_step(solution.grid)
    Δt = calculate_Δt(solution, config)

    progress = create_simulation_progress(solution.grid, method)

    while solution.t < config.t_max
        solution = step!(method, solution, config, Δt)
        Δt_prev = Δt

        run_after_step(hooks, solution, config, Δt_prev)

        ProgressMeter.next!(
            progress;
            showvalues = [
                (:Прогресс, @sprintf("%3.1f%%", solution.t / config.t_max * 100)),
                (:Ячеек, @sprintf("%3d", length(solution.grid.nodes) - 1)),
                (:Δt, @sprintf("%.2e", Δt))
            ]
        )

        birth_nodes && birth_nodes!(solution, config, Δx)

        Δt = calculate_Δt(solution, config)
        Δt = limit_time_step(Δt, solution.t, config)

        run_after_birth(hooks, solution, config, Δt_prev, Δt)
    end

    run_after_loop(hooks, solution, config)

    finish!(progress)

    remove_ghost_cells!(solution)

    return solution
end

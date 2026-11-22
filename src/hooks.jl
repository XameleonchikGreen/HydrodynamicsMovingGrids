function run_after_step(
    hooks::HookSet,
    solution,
    config,
    Δt_prev
)
    for hook in hooks.after_step
        hook(solution, config, Δt_prev)
    end

    return nothing
end

function run_after_birth(
    hooks::HookSet,
    solution,
    config,
    Δt_prev,
    Δt
)
    for hook in hooks.after_birth
        hook(solution, config, Δt_prev, Δt)
    end

    return nothing
end

function run_after_loop(
    hooks::HookSet,
    solution,
    config,
)
    for hook in hooks.after_loop
        hook(solution, config)
    end

    return nothing
end

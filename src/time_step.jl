function calculate_Δt(solution::Solution{T}, config::ProblemConfig) where T
    γ = config.physics.γ
    CFL = config.physics.CFL
    grid = solution.grid
    
    Δx = diff(grid.nodes)
    c = sound_speed.(solution.Ws[2:end-1], γ)
    u = [W.u for W in solution.Ws[2:end-1]]
    Q̇ = abs.(diff(grid.Qs))
    
    λ = c .+ abs.(u) .+ Q̇
    
    return CFL * minimum(Δx ./ λ)
end

function limit_time_step(Δt::BigFloat, t::BigFloat, config::ProblemConfig)
    return max(min(config.t_max - t, Δt), config.min_time_step)
end

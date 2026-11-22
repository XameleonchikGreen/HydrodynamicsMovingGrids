function add_ghost_cells!(solution::Solution, config::ProblemConfig)
    @assert length(solution.Ws) == length(solution.grid.nodes) - 1
        "add_ghost_cells! expects a solution without ghost cells"
    
    Ws = solution.Ws
    t = solution.t

    pushfirst!(Ws, first(Ws))
    push!(Ws, last(Ws))

    left = config.bc.left(t, Ws)
    right = config.bc.right(t, Ws)

    Ws[1] = left
    Ws[end] = right

    return solution
end

function remove_ghost_cells!(solution::Solution)
    @assert length(solution.Ws) == length(solution.grid.nodes) + 1
        "remove_ghost_cells! expects a solution with ghost cells"
    
    pop!(solution.Ws)
    popfirst!(solution.Ws)

    return solution
end

function apply_bc!(solution::Solution, config::ProblemConfig)
    @assert length(solution.Ws) == length(solution.grid.nodes) + 1
        "apply_bc! expects a solution with ghost cells"
    
    Ws = solution.Ws
    t = solution.t

    left = config.bc.left(t, Ws)
    right = config.bc.right(t, Ws)

    Ws[1] = left
    Ws[end] = right

    return solution
end

mirror(W::PrimVars) = PrimVars(W.ρ, -W.u, W.p)

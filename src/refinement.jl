function birth_nodes!(solution::Solution{T}, config::ProblemConfig, Δx::BigFloat) where T
    grid = solution.grid
    Ws = solution.Ws
    γ = config.physics.γ
    ε = config.physics.ε_newton
    i = 2
    while i < length(grid.nodes)
        if grid.nodes[i + 1] - grid.nodes[i] > config.refinement_factor * Δx
            insert_node!(grid, i)
            insert!(Ws, i + 1, Ws[i])
            Ws[i + 1] = W⁺(Ws[i + 1], W⁺(Ws[i + 1], Ws[i + 2], γ, ε)(- grid.Qs[i + 2]), 
                γ, ε)(- (grid.Qs[i + 1] + grid.Qs[i + 2]) / 2)
            Ws[i] = W⁺(W⁺(Ws[i - 1], Ws[i], γ, ε)(- grid.Qs[i]), Ws[i],
                γ, ε)(- (grid.Qs[i] + grid.Qs[i + 1]) / 2)
        else
            i += 1
        end
    end
end

function insert_node!(grid::Grid{Moving}, idx::Integer)
    insert!(grid.nodes, idx + 1, (grid.nodes[idx + 1] + grid.nodes[idx]) / 2)
    grid.max_id += 1
    insert!(grid.ids, idx + 1, grid.max_id)
    insert!(grid.Qs, idx + 1, (grid.Qs[idx + 1] + grid.Qs[idx]) / 2)
    insert!(grid.ψs, idx + 1, grid.ψs[idx])
    grid.Δy[idx] /= 2
    insert!(grid.Δy, idx + 1, grid.Δy[idx])
end

cell_centers(grid::Grid{T}) where T = (grid.nodes[2:end] .+ grid.nodes[1:end-1]) ./ 2
min_step(grid::Grid{T}) where T = minimum(grid.Δy)

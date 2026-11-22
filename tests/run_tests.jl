using Plots
using LaTeXStrings
using CSV
using DataFrames
using Polynomials

include(joinpath(@__DIR__, "..", "GasDynamics.jl"))
using .GasDynamics

include("setup_tests.jl")
include("plotting.jl")
include("plotting_err.jl")

function run_plot(
    method::NumericalMethod,
    ::Type{T},
    n_cells::Int,
    test_case::Union{Int,Symbol},
    mesh_velocity_alg::Symbol,
    birth_nodes::Bool,
    take_snapshots::Bool;
    save_tiff::Bool=false
) where T <: GridType

    outdir = joinpath(@__DIR__, "..", "images", "$(test_case)_$(mesh_velocity_alg)")

    simulate_and_plot(
        T,
        test_case,
        n_cells,
        method,
        mesh_velocity_alg;
        birth_nodes=birth_nodes,
        snapshots=take_snapshots,
        outdir = outdir,
        save_tiff = save_tiff
    )

    return nothing
end

function run_err(
    ns_cells::Vector{Int},
    test_case::Int,
    mesh_velocity_alg::Symbol;
    save_tiff::Bool=false
)

    outdir = joinpath(@__DIR__, "..", "images", "$(test_case)_$(mesh_velocity_alg)")

    simulate_and_plot_err(
        test_case,
        ns_cells,
        mesh_velocity_alg;
        outdir = outdir,
        save_tiff = save_tiff
    )

    return nothing
end

function run_and_print_err(
    method::NumericalMethod,
    ::Type{T},
    n_cells::Int,
    test_case::Int,
    mesh_velocity_alg::Symbol,
    birth_nodes::Bool
) where T <: GridType


    simulate_and_print_err(
        T,
        test_case,
        n_cells,
        method,
        mesh_velocity_alg;
        birth_nodes=birth_nodes
    )

    return nothing
end

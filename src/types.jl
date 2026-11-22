struct PrimVars
    ρ::BigFloat
    u::BigFloat
    p::BigFloat
end

PrimVars((ρ, u, p)) = PrimVars(ρ, u, p)


struct PhysicsParams
    γ::BigFloat
    ε_newton::BigFloat
    ε_pressure::BigFloat
    ε_density::BigFloat
    α::BigFloat
    CFL::BigFloat
end


abstract type GridType end
struct Static <: GridType end
struct Moving <: GridType end

grid_name(::Type{Static}) = "Static"
grid_name(::Type{Moving}) = "Moving"


abstract type NumericalMethod end
struct GodunovScheme <: NumericalMethod end
struct MUSCLScheme <: NumericalMethod end

method_name(::GodunovScheme) = "Godunov"
method_name(::MUSCLScheme) = "MUSCL"


mutable struct Grid{T<:GridType}
    nodes::Vector{BigFloat}
    ids::Vector{Int}
    max_id::Int
    Δy::Vector{BigFloat}
    ψs::Vector{BigFloat}
    Qs::Vector{BigFloat}

    function Grid(::Type{T}, N::Integer, x_min::BigFloat, x_max::BigFloat) where T
        @assert N > 1
        @assert x_max > x_min
        nodes = collect(range(x_min, x_max, length=N+1))
        max_id = N+1
        ids = collect(1:(N+1))
        Δy = diff(nodes)
        ψs = ones(BigFloat, N)
        Qs = zeros(BigFloat, N + 1)
        return new{T}(nodes, ids, max_id, Δy, ψs, Qs)
    end
end

struct BoundaryConditions
    left::Function   
    right::Function
end

struct InitialCondition
    func::Function
end


mutable struct Solution{T<:GridType}
    grid::Grid{T}
    Ws::Vector{PrimVars}
    t::BigFloat

    function Solution{T}(grid::Grid{T}, Ws::Vector{PrimVars}, t::BigFloat) where T
        @assert length(Ws) == length(grid.nodes) - 1
        return new{T}(grid, Ws, t)
    end
end

function Solution{T}(grid::Grid{T}, ic::InitialCondition, t::BigFloat) where T<:GridType
    return Solution{T}(grid, ic.func.(cell_centers(grid)), t)
end


struct ProblemConfig
    physics::PhysicsParams
    bc::BoundaryConditions
    ic::InitialCondition
    t_max::BigFloat
    refinement_factor::BigFloat
    min_time_step::BigFloat
    mesh_velocity_alg::Symbol
end

function ProblemConfig(;
    γ::BigFloat,
    ε_newton::BigFloat,
    ε_pressure::BigFloat,
    ε_density::BigFloat,
    α::BigFloat,
    CFL::BigFloat,
    t_max::BigFloat,
    BCₗ::Function,
    BCᵣ::Function,
    IC::Function,
    refinement_factor::BigFloat,
    min_time_step::BigFloat,
    mesh_velocity_alg::Symbol
)
    physics = PhysicsParams(γ, ε_newton, ε_pressure, ε_density, α, CFL)
    bc = BoundaryConditions(BCₗ, BCᵣ)
    ic = InitialCondition(IC)
    
    return ProblemConfig(
        physics,
        bc,
        ic,
        t_max,
        refinement_factor,
        min_time_step,
        mesh_velocity_alg
    )
end


@enum ErrorNorm L1 L2 Linf


mutable struct ErrorAccumulator
    ρ::BigFloat
    u::BigFloat
    p::BigFloat
    energy::BigFloat
end

function ErrorAccumulator()
    return ErrorAccumulator(0.0, 0.0, 0.0, 0.0)
end

struct GridSnapshot
    nodes::Vector{BigFloat}
    ids::Vector{Int}
    max_id::Int
    t::BigFloat

    function GridSnapshot(
        nodes::Vector{BigFloat},
        ids::Vector{Int},
        max_id::Int,
        t::BigFloat
    )
        return new(copy(nodes), copy(ids), max_id, t)
    end
end


struct HookSet
    after_step::Vector{Function}
    after_birth::Vector{Function}
    after_loop::Vector{Function}
end

function HookSet()
    return HookSet(
        Function[],
        Function[],
        Function[]
    )
end

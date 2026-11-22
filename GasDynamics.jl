module GasDynamics

using ProgressMeter
using Printf

include("src/types.jl")
include("src/hooks.jl")
include("src/grid.jl")
include("src/physics.jl")
include("src/riemann.jl")
include("src/mesh_velocity.jl")
include("src/conservative.jl")
include("src/reconstruction.jl")
include("src/boundary.jl")
include("src/time_step.jl")
include("src/godunov.jl")
include("src/muscl.jl")
include("src/refinement.jl")
include("src/diagnostics.jl")
include("src/solver.jl")
include("src/api.jl")

export
    # Типы
    PrimVars,
    PhysicsParams,
    GridType,
    Static,
    Moving,
    Grid,
    BoundaryConditions,
    InitialCondition,
    Solution,
    ProblemConfig,
    ErrorNorm,
    L1,
    L2,
    Linf,
    ErrorAccumulator,
    GridSnapshot,
    GodunovScheme,
    MUSCLScheme,
    NumericalMethod,

    # Сетка
    cell_centers,
    min_step,
    insert_node!,

    # Граничные условия
    mirror,

    # Физика
    sound_speed,
    calculate_energy,

    # Запуск методов
    simulate,

    # Точное решение
    W⁺

end
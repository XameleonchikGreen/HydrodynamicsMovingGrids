function load_tests(filepath::String)::DataFrame
    df = CSV.read(filepath, DataFrame)

    float_cols = [:ρₗ, :uₗ, :pₗ, :ρᵣ, :uᵣ, :pᵣ, :t]

    for col in float_cols
        if hasproperty(df, col)
            df[!, col] = BigFloat.(string.(df[!, col]))
        end
    end

    if hasproperty(df, :Test)
        df[!, :Test] = parse.(Int, string.(df[!, :Test]))
    end

    return df
end

function test_row(tests::DataFrame, case_id::Int)
    idx = findfirst(==(case_id), tests.Test)

    if idx === nothing
        error("Test $case_id not found in CSV")
    end

    return tests[idx, :]
end

function make_problem(row)
    Wₗ = PrimVars(row[:ρₗ], row[:uₗ], row[:pₗ])
    Wᵣ = PrimVars(row[:ρᵣ], row[:uᵣ], row[:pᵣ])

    BCₗ(t, Ws) = Wₗ
    BCᵣ(t, Ws) = Wᵣ

    IC(x) = x >= 0 ? Wᵣ : Wₗ

    t_max = row[:t]

    return BCₗ, BCᵣ, IC, Wₗ, Wᵣ, t_max
end

function make_config(BCₗ, BCᵣ, IC, t_max, grid, mesh_velocity_alg)
    bc = BoundaryConditions(BCₗ, BCᵣ)
    ic = InitialCondition(IC)

    γ = BigFloat("1.4")

    W_ic = IC.(cell_centers(grid))

    char_speed = 3 * (maximum(sound_speed.(W_ic, γ)) + 
                      maximum(abs, getfield.(W_ic, :u)))

    physics = PhysicsParams(
        γ,                     # γ
        BigFloat("1e-15"),     # ε_newton
        BigFloat("0.002"),     # ε_pressure
        BigFloat("0.0002"),    # ε_density
        char_speed * min_step(grid) / 10,
        BigFloat("0.4")        # CFL
    )

    return ProblemConfig(
        physics,
        bc,
        ic,
        t_max,                 # t_max
        BigFloat("2"),         # refinement_factor
        BigFloat("1e-6"),      # min_time_step
        mesh_velocity_alg      # choose mesh_velocity by density or pressure
    )
end

function exact_solution_arrays(Wₗ, Wᵣ, config, cells)
    γ = config.physics.γ
    ε = config.physics.ε_newton
    t = config.t_max

    W_exact = W⁺(Wₗ, Wᵣ, γ, ε)

    N = length(cells)

    ρ_ex = Vector{BigFloat}(undef, N)
    u_ex = Vector{BigFloat}(undef, N)
    p_ex = Vector{BigFloat}(undef, N)
    E_ex = Vector{BigFloat}(undef, N)

    for i in eachindex(cells)
        S = cells[i] / t
        W = W_exact(S)

        ρ_ex[i] = W.ρ
        u_ex[i] = W.u
        p_ex[i] = W.p
        E_ex[i] = W.p / (γ - 1) / W.ρ
    end

    return ρ_ex, u_ex, p_ex, E_ex
end

function solution_arrays(solution, γ)
    Ws = solution.Ws

    ρ = [W.ρ for W in Ws]
    u = [W.u for W in Ws]
    p = [W.p for W in Ws]

    E = @. p / (γ - 1) / ρ

    return ρ, u, p, E
end

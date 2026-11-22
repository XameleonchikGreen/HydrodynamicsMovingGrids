function setup_plots!()
    default(:linewidth, 2)
    default(:c, "#000000")
    default(:dpi, 600)
    default(:xlabel, L"x")
    default(:legend, false)
    default(:markersize, 2.5)
    default(:markerstrokewidth, 0)
    default(:xguidefontsize, 14)
    default(:yguidefontsize, 14)
    default(:xtickfontsize, 12)
    default(:ytickfontsize, 12)
end

function save_plot_image!(
    p,
    method::AbstractString,
    N::Integer,
    ::Type{T},
    birth_nodes::Bool,
    outdir::AbstractString,
    save_tiff::Bool
) where T <: GridType
    mkpath(outdir)

    safe(s) = replace(string(s), r"[^A-Za-z0-9_-]+" => "_")

    method_name = safe(method)
    grid_name   = safe(nameof(T))
    birth_name  = birth_nodes ? "birth" : "nobirth"

    filename = "$(method_name)_N$(N)_$(grid_name)_$(birth_name)"
    path = joinpath(outdir, filename)

    png_path = path * ".png"

    savefig(p, png_path)
    println("  saved: ", path)

    if save_tiff
        tif_path = path * ".tif"
        run(`convert $png_path -depth 8 -compress LZW $tif_path`)
    end

    return path
end

function plot_solution_fields(cells_ex, cells, exact, numerical, ψ)
    ρ_ex, u_ex, p_ex, E_ex = exact
    ρ_n, u_n, p_n, E_n = numerical

    ψ_min = minimum(ψ)
    ψ_max = maximum(ψ)

    x_ex = Float64.(cells_ex)
    x = Float64.(cells)

    p1 = plot(x_ex, Float64.(ρ_ex); ylabel=L"\rho")
    if ψ_max - ψ_min > BigFloat("1e-5")
        ρ_min = min(minimum(ρ_ex), minimum(ρ_n))
        ρ_max = max(maximum(ρ_ex), maximum(ρ_n))
        ψ_ρ = (ψ .- ψ_min) ./ (ψ_max - ψ_min) .* (ρ_max - ρ_min) .+ ρ_min
        plot!(p1, x, Float64.(ψ_ρ); c="#0072B2", ls=:dash, lw=1)
    end
    scatter!(p1, x, Float64.(ρ_n); c="#CC79A7")

    p2 = plot(x_ex, Float64.(u_ex); ylabel=L"u")
    if ψ_max - ψ_min > BigFloat("1e-5")
        u_min = min(minimum(u_ex), minimum(u_n))
        u_max = max(maximum(u_ex), maximum(u_n))
        ψ_u = (ψ .- ψ_min) ./ (ψ_max - ψ_min) .* (u_max - u_min) .+ u_min
        plot!(p2, x, Float64.(ψ_u); c="#0072B2", ls=:dash, lw=1)
    end
    scatter!(p2, x, Float64.(u_n); c="#CC79A7")

    p3 = plot(x_ex, Float64.(p_ex); ylabel=L"p")
    if ψ_max - ψ_min > BigFloat("1e-5")
        p_min = min(minimum(p_ex), minimum(p_n))
        p_max = max(maximum(p_ex), maximum(p_n))
        ψ_p = (ψ .- ψ_min) ./ (ψ_max - ψ_min) .* (p_max - p_min) .+ p_min
        plot!(p3, x, Float64.(ψ_p); c="#0072B2", ls=:dash, lw=1)
    end
    scatter!(p3, x, Float64.(p_n); c="#CC79A7")

    p4 = plot(x_ex, Float64.(E_ex); ylabel=L"E")
    if ψ_max - ψ_min > BigFloat("1e-5")
        E_min = min(minimum(E_ex), minimum(E_n))
        E_max = max(maximum(E_ex), maximum(E_n))
        ψ_E = (ψ .- ψ_min) ./ (ψ_max - ψ_min) .* (E_max - E_min) .+ E_min
        plot!(p4, x, Float64.(ψ_E); c="#0072B2", ls=:dash, lw=1)
    end
    scatter!(p4, x, Float64.(E_n); c="#CC79A7")

    return plot(p1, p2, p3, p4; layout=(2, 2), size=(800, 450))
end

function save_quants(cells_ex, cells, exact, numerical, ψ, path)
    ρ_ex, u_ex, p_ex, E_ex = exact
    ρ_n, u_n, p_n, E_n = numerical

    x_ex = Float64.(cells_ex)
    x = Float64.(cells)
    
    df_n = DataFrame(
        x = x,
        rho = ρ_n, u = u_n, p = p_n, E = E_n,
        psi = ψ
    )
    CSV.write(path * ".csv", df_n)

    df_ex = DataFrame(
        x = x_ex,
        rho = ρ_ex, u = u_ex, p = p_ex, E = E_ex
    )
    CSV.write(path * "_exact.csv", df_ex)
end

function plot_snapshots(snapshots, max_id, N)
    p = plot(;
        xlabel=L"x",
        ylabel=L"t",
        margin=5Plots.mm,
        size=(900, 500)
    )

    for id in 1:max_id
        t = Float64[]
        x = Float64[]
        c = id < N+2 ? "#0072B2" : "#E69F00"

        for snapshot in snapshots
            snapshot.max_id < id && continue

            push!(t, Float64.(snapshot.t))

            i = findfirst(==(id), snapshot.ids)
            push!(x, Float64.(snapshot.nodes[i]))
        end

        if !isempty(x)
            plot!(
                p,
                x,
                t;
                c=c
            )
        end
    end

    return p
end

function save_trajectories(snapshots, max_id, N, path)
    t = Float64[]
    for snapshot in snapshots
        push!(t, Float64.(snapshot.t))
    end
    df = DataFrame(t = t)
    for id in 1:(N+1)
        x = Float64[]
        for snapshot in snapshots
            i = findfirst(==(id), snapshot.ids)
            push!(x, Float64.(snapshot.nodes[i]))
        end
        df[!, Symbol(id)] = x
    end
    CSV.write(path * "_traj_old.csv", df)

    df = DataFrame(t = t)
    for id in N+2:max_id
        x = []
        for snapshot in snapshots
            if snapshot.max_id < id
                push!(x, missing)
            else
                i = findfirst(==(id), snapshot.ids)
                push!(x, Float64.(snapshot.nodes[i]))
            end
        end
        df[!, Symbol(id)] = x
    end
    CSV.write(path * "_traj_new.csv", df)
end

function simulate_and_plot(
    ::Type{T},
    case_id::Int,
    N::Integer,
    method::GasDynamics.NumericalMethod,
    mesh_velocity_alg::Symbol;
    birth_nodes::Bool=false,
    snapshots::Bool=false,
    snapshot_step::BigFloat=BigFloat("2e-3"),
    outdir::AbstractString,
    save_tiff::Bool=false
) where T <: GridType

    setup_plots!()

    grid = Grid(T, N, BigFloat("-1"), BigFloat("1"))

    row = test_row(load_tests(joinpath(@__DIR__, "./tests.csv")), case_id)

    BCₗ, BCᵣ, IC, Wₗ, Wᵣ, t_max = make_problem(row)
    config = make_config(BCₗ, BCᵣ, IC, t_max, grid, mesh_velocity_alg)

    @time result =
        if snapshots
            GasDynamics.simulate(
                grid,
                config,
                method;
                birth_nodes=birth_nodes,
                snapshot_step=snapshot_step
            )
        else
            GasDynamics.simulate(
                grid,
                config,
                method;
                birth_nodes=birth_nodes
            )
        end

    γ = config.physics.γ
    cells = cell_centers(grid)
    cells_ex = collect(range(BigFloat("-1"), BigFloat("1"), length=2N+1))

    exact = exact_solution_arrays(Wₗ, Wᵣ, config, cells_ex)
    numerical = solution_arrays(result.solution, γ)
    ψ = 1 ./ result.solution.grid.ψs

    p = plot_solution_fields(cells_ex, cells, exact, numerical, ψ)

    mname = GasDynamics.method_name(method)

    path = save_plot_image!(
        p,
        mname,
        N,
        T,
        birth_nodes,
        outdir,
        save_tiff
    )

    save_quants(cells_ex, cells, exact, numerical, ψ, path)

    display(p)

    # Если снапшоты собирались, рисуем их отдельно.
    if snapshots
        snaps = get(result, :snapshots, nothing)

        if snaps !== nothing && !isempty(snaps)
            p_traj = plot_snapshots(snaps, grid.max_id, N)

            save_plot_image!(
                p_traj,
                mname * "_Trajectories",
                N,
                T,
                birth_nodes,
                outdir,
                save_tiff
            )

            save_trajectories(snaps, grid.max_id, N, path)

            display(p_traj)
        end
    end

    return nothing
end

function simulate_and_plot(
    ::Type{T},
    case_id::Symbol,
    N::Integer,
    method::GasDynamics.NumericalMethod,
    mesh_velocity_alg::Symbol;
    birth_nodes::Bool=false,
    snapshots::Bool=false,
    snapshot_step::BigFloat=BigFloat("2e-3"),
    outdir::AbstractString,
    save_tiff::Bool=false
) where T <: GridType

    setup_plots!()

    if case_id == :shu_osher

        grid = Grid(T, N, BigFloat("-5"), BigFloat("5"))

        function BCₗ(t::BigFloat, Ws::Vector{PrimVars})
            return PrimVars(BigFloat("3.857143"), BigFloat("2.629369"), BigFloat("10.33333"))
        end

        function BCᵣ(t::BigFloat, Ws::Vector{PrimVars})
            PrimVars(BigFloat("1") + BigFloat("0.2") * sin(25), BigFloat("0"), BigFloat("1"))
        end
        
        function IC(x::BigFloat)
            if x ≤ BigFloat("-4")
                return PrimVars(BigFloat("3.857143"), BigFloat("2.629369"), BigFloat("10.33333"))
            end
            return PrimVars(BigFloat("1") + BigFloat("0.2") * sin(5x), BigFloat("0"), BigFloat("1"))
        end

        config = make_config(BCₗ, BCᵣ, IC, BigFloat("1.8"), grid, mesh_velocity_alg)

        result =
            if snapshots
                GasDynamics.simulate(
                    grid,
                    config,
                    method;
                    birth_nodes=birth_nodes,
                    snapshot_step=snapshot_step
                )
            else
                GasDynamics.simulate(
                    grid,
                    config,
                    method;
                    birth_nodes=birth_nodes
                )
            end

        γ = config.physics.γ
        numerical = solution_arrays(result.solution, γ)
        ψ = 1 ./ result.solution.grid.ψs

        path_to_exact = joinpath(@__DIR__, "..", "..", "CFD_Julia", "09_Euler_1D_Roe")

        df_ρ = CSV.read(
            joinpath(path_to_exact, "solution_d.txt"),
            DataFrame;
            delim=' ',
            header=false,
            ignorerepeated=true
        )
        rename!(df_ρ, [:x, :rho_0, :rho])

        df_ρu = CSV.read(
            joinpath(path_to_exact, "solution_v.txt"),
            DataFrame;
            delim=' ',
            header=false,
            ignorerepeated=true
        )
        rename!(df_ρu, [:x, :rho_u_0, :rho_u])

        df_e = CSV.read(
            joinpath(path_to_exact, "solution_e.txt"),
            DataFrame;
            delim=' ',
            header=false,
            ignorerepeated=true
        )
        rename!(df_e, [:x, :e_0, :e])

        df_ρu[:, :u] = df_ρu[:, :rho_u] ./ df_ρ[:, :rho]
        df_e[:, :p] = (df_e[:, :e] - df_ρu[:, :rho_u] .* df_ρu[:, :u] ./ 2) .* (γ - 1)

        cells_ex = df_ρ[:, :x]
        exact = (df_ρ[:, :rho], df_ρu[:, :u], df_e[:, :p], df_e[:, :e])

        cells = cell_centers(grid)

        p = plot_solution_fields(cells_ex, cells, exact, numerical, ψ)

        mname = GasDynamics.method_name(method)

        save_plot_image!(
            p,
            mname,
            N,
            T,
            birth_nodes,
            outdir,
            save_tiff
        )

        display(p)

        # Если снапшоты собирались, рисуем их отдельно.
        if snapshots
            snaps = get(result, :snapshots, nothing)

            if snaps !== nothing && !isempty(snaps)
                p_traj = plot_snapshots(snaps, grid.max_id, N)

                save_plot_image!(
                    p_traj,
                    mname * "_Trajectories",
                    N,
                    T,
                    birth_nodes,
                    outdir,
                    save_tiff
                )

                display(p_traj)
            end
        end
    end

    return nothing
end


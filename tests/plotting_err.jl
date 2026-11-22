function setup_err_plots!()
    default(:linewidth, 2)
    default(:c, "#000000")
    default(:dpi, 600)
    default(:xlabel, L"h")
    default(:legend, false)
    default(:xscale, :log)
    default(:yscale, :log)
    default(:markersize, 6)
    default(:markerstrokewidth, 0)
    default(:xguidefontsize, 14)
    default(:yguidefontsize, 14)
    default(:xtickfontsize, 14)
    default(:ytickfontsize, 14)
    default(:legendfontsize, 14)
    default(:framestyle, :box)
    default(:fontfamily, :times)
end

function calculate_degree(acc, h)
    p1 = fit(log10.(h), log10.(getfield.(acc, :ρ)), 1)[1]
    p1_text = string(round(Float64(p1), digits=4))

    p2 = fit(log10.(h), log10.(getfield.(acc, :u)), 1)[1]
    p2_text = string(round(Float64(p2), digits=4))

    p3 = fit(log10.(h), log10.(getfield.(acc, :p)), 1)[1]
    p3_text = string(round(Float64(p3), digits=4))

    p4 = fit(log10.(h), log10.(getfield.(acc, :energy)), 1)[1]
    p4_text = string(round(Float64(p4), digits=4))

    return p1_text, p2_text, p3_text, p4_text
end

function unpack_error_accumulator(vec_acc)
    ρ = [acc.ρ for acc in vec_acc]
    u = [acc.u for acc in vec_acc]
    p = [acc.p for acc in vec_acc]
    E = [acc.energy for acc in vec_acc]

    return ρ, u, p, E
end

function plot_errs(
    acc_st,
    acc_mv,
    acc_bn,
    acc_mc,
    Ns,
    final_N
)
    ρ_st, u_st, p_st, E_st = unpack_error_accumulator(acc_st)
    ρ_mv, u_mv, p_mv, E_mv = unpack_error_accumulator(acc_mv)
    ρ_bn, u_bn, p_bn, E_bn = unpack_error_accumulator(acc_bn)
    ρ_mc, u_mc, p_mc, E_mc = unpack_error_accumulator(acc_mc)

    h = 2 ./ Ns
    h_float = Float64.(h)

    p1 = plot(
        h_float,
        Float64.(ρ_st);
        ylabel=L"\rho_{err}",
        c="#0072B2",
        marker=:circle
    )
    annotate!(p1, h .* 0.9, ρ_st .* 1.15, text.(string.(Ns), 10, "#0072B2", :left))
    plot!(p1, h, ρ_mv, c="#E69F00", marker=:diamond)
    plot!(p1, h, ρ_bn, c="#009E73", marker=:utriangle)
    annotate!(p1, h .* 0.9, ρ_bn * 1.15, text.(string.(final_N), 10, "#009E73", :left))
    plot!(p1, h, ρ_mc, c="#CC79A7", marker=:rect)

    p2 = plot(
        h_float,
        Float64.(u_st);
        ylabel=L"u_{err}",
        c="#0072B2",
        marker=:circle
    )
    annotate!(p2, h .* 0.9, u_st .* 1.15, text.(string.(Ns), 10, "#0072B2", :left))
    plot!(p2, h, u_mv, c="#E69F00", marker=:diamond)
    plot!(p2, h, u_bn, c="#009E73", marker=:utriangle)
    annotate!(p2, h .* 0.9, u_bn * 1.15, text.(string.(final_N), 10, "#009E73", :left))
    plot!(p2, h, u_mc, c="#CC79A7", marker=:rect)

    p3 = plot(
        h_float,
        Float64.(p_st);
        ylabel=L"p_{err}",
        c="#0072B2",
        marker=:circle
    )
    annotate!(p3, h .* 0.9, p_st .* 1.15, text.(string.(Ns), 10, "#0072B2", :left))
    plot!(p3, h, p_mv, c="#E69F00", marker=:diamond)
    plot!(p3, h, p_bn, c="#009E73", marker=:utriangle)
    annotate!(p3, h .* 0.9, p_bn * 1.15, text.(string.(final_N), 10, "#009E73", :left))
    plot!(p3, h, p_mc, c="#CC79A7", marker=:rect)

    p4 = plot(
        h_float,
        Float64.(E_st);
        ylabel=L"E_{err}",
        c="#0072B2",
        marker=:circle
    )
    annotate!(p4, h .* 0.9, E_st .* 1.15, text.(string.(Ns), 10, "#0072B2", :left))
    plot!(p4, h, E_mv, c="#E69F00", marker=:diamond)
    plot!(p4, h, E_bn, c="#009E73", marker=:utriangle)
    annotate!(p4, h .* 0.9, E_bn * 1.15, text.(string.(final_N), 10, "#009E73", :left))
    plot!(p4, h, E_mc, c="#CC79A7", marker=:rect)

    for p in [p1, p2, p3, p4]
        y_min, y_max = ylims(p)
        x_min, x_max = xlims(p)
        plot!(
            p,
            grid=true, minorgrid=true,
            xlims=(x_min * 0.9, x_max * 1.2),
            ylims=(y_min * 0.95, y_max * 1.2)
        )
    end

    legend_plot = plot(
        grid=false,
        showaxis=false,
        framestyle=:none,
        legend=:bottom,
        size=(800, 80),
        ticks=nothing,
        xticks=nothing,
        yticks=nothing,
        xlims = (1, 2),
        ylims = (1, 2),
        foreground_color_legend = RGBA(0, 0, 0, 0)
    )

    plot!(legend_plot, [0], [0], label="Равномерная сетка", c="#0072B2", marker=:circle)
    plot!(legend_plot, [0], [0], label="Адаптивная сетка", c="#E69F00", marker=:diamond)
    plot!(legend_plot, [0], [0], label="Измельчающаяся сетка", c="#009E73", marker=:utriangle)
    plot!(legend_plot, [0], [0], label="MUSCL", c="#CC79A7", marker=:rect)

    layout = @layout [a b; c d; e{0.15h}]

    return plot(p1, p2, p3, p4, legend_plot, layout=layout, size=(800, 700))
end

function save_errs(
    acc_st,
    acc_mv,
    acc_bn,
    acc_mc,
    Ns,
    final_N,
    path
)
    ρ_st, u_st, p_st, E_st = unpack_error_accumulator(acc_st)
    ρ_mv, u_mv, p_mv, E_mv = unpack_error_accumulator(acc_mv)
    ρ_bn, u_bn, p_bn, E_bn = unpack_error_accumulator(acc_bn)
    ρ_mc, u_mc, p_mc, E_mc = unpack_error_accumulator(acc_mc)

    df = DataFrame(
        Ns = Ns, final_N = final_N,
        rho_st = ρ_st, u_st = u_st, p_st = p_st, E_st = E_st,
        rho_mv = ρ_mv, u_mv = u_mv, p_mv = p_mv, E_mv = E_mv,
        rho_bn = ρ_bn, u_bn = u_bn, p_bn = p_bn, E_bn = E_bn,
        rho_mc = ρ_mc, u_mc = u_mc, p_mc = p_mc, E_mc = E_mc,
    )
    CSV.write(joinpath(path, "plot_errs.csv"), df)
end

function simulate_and_plot_err(
    case_id::Int,
    Ns::Vector{Int},
    mesh_velocity_alg::Symbol;
    norm::GasDynamics.ErrorNorm=GasDynamics.L1,
    outdir::AbstractString,
    save_tiff::Bool=false
)
    setup_err_plots!()

    row = test_row(load_tests(joinpath(@__DIR__, "./tests.csv")), case_id)

    BCₗ, BCᵣ, IC, Wₗ, Wᵣ = make_problem(row)
    
    acc_st = []
    acc_mv = []
    acc_bn = []
    acc_mc = []
    final_N = []

    for N in Ns
        grid_st = Grid(Static, N, BigFloat("-1"), BigFloat("1"))
        config_st = make_config(BCₗ, BCᵣ, IC, t_max, grid_st, mesh_velocity_alg)
        W_st = GasDynamics.W⁺(Wₗ, Wᵣ, config_st.physics.γ, config_st.physics.ε_newton)
        res_st = GasDynamics.simulate(
            grid_st,
            config_st,
            GodunovScheme();
            exact_solution=W_st,
            error_norm=norm
        )
        push!(acc_st, res_st.error)

        grid_mv = Grid(Moving, N, BigFloat("-1"), BigFloat("1"))
        config_mv = make_config(BCₗ, BCᵣ, IC, t_max, grid_mv, mesh_velocity_alg)
        W_mv = GasDynamics.W⁺(Wₗ, Wᵣ, config_mv.physics.γ, config_mv.physics.ε_newton)
        res_mv = GasDynamics.simulate(
            grid_mv,
            config_mv,
            GodunovScheme();
            exact_solution=W_mv,
            error_norm=norm
        )
        push!(acc_mv, res_mv.error)

        grid_bn = Grid(Moving, N, BigFloat("-1"), BigFloat("1"))
        config_bn = make_config(BCₗ, BCᵣ, IC, t_max, grid_bn, mesh_velocity_alg)
        W_bn = GasDynamics.W⁺(Wₗ, Wᵣ, config_bn.physics.γ, config_bn.physics.ε_newton)
        res_bn = GasDynamics.simulate(
            grid_bn,
            config_bn,
            GodunovScheme();
            birth_nodes=true,
            exact_solution=W_bn,
            error_norm=norm
        )
        push!(acc_bn, res_bn.error)
        push!(final_N, length(grid_bn.nodes) - 1)

        grid_mc = Grid(Static, N, BigFloat("-1"), BigFloat("1"))
        config_mc = make_config(BCₗ, BCᵣ, IC, t_max, grid_mc, mesh_velocity_alg)
        W_mc = GasDynamics.W⁺(Wₗ, Wᵣ, config_mc.physics.γ, config_mc.physics.ε_newton)
        res_mc = GasDynamics.simulate(
            grid_mc,
            config_mc,
            MUSCLScheme();
            exact_solution=W_mc,
            error_norm=norm
        )
        push!(acc_mc, res_mc.error)

    end

    p = plot_errs(acc_st, acc_mv, acc_bn, acc_mc, Ns, final_N)

    filename = "Approx_degrees_N" * join(Ns, "_")
    path = joinpath(outdir, filename)

    png_path = path * ".png"

    savefig(p, png_path)
    save_errs(acc_st, acc_mv, acc_bn, acc_mc, Ns, final_N, outdir)
    println("  saved: ", path)

    if save_tiff
        tif_path = path * ".tif"
        run(`convert $png_path -depth 8 -compress LZW $tif_path`)
    end
    
    display(p)

    h = 2 ./ Ns

    results = [
        calculate_degree(acc_st, h),
        calculate_degree(acc_mv, h),
        calculate_degree(acc_bn, h),
        calculate_degree(acc_mc, h)
    ]

    open(path * ".txt", "w") do file
        for res in results
            println(file, join(res, ",")) 
        end
    end

    return nothing
end

function simulate_and_print_err(
    ::Type{T},
    case_id::Int,
    N::Integer,
    method::GasDynamics.NumericalMethod,
    mesh_velocity_alg::Symbol;
    norm::GasDynamics.ErrorNorm=GasDynamics.L1,
    birth_nodes::Bool=false
) where T <: GridType

    setup_plots!()

    grid = Grid(T, N, BigFloat("-1"), BigFloat("1"))

    row = test_row(load_tests(joinpath(@__DIR__, "./tests.csv")), case_id)

    BCₗ, BCᵣ, IC, Wₗ, Wᵣ, t_max = make_problem(row)
    config = make_config(BCₗ, BCᵣ, IC, t_max, grid, mesh_velocity_alg)
    W = GasDynamics.W⁺(Wₗ, Wᵣ, config.physics.γ, config.physics.ε_newton)

    result = GasDynamics.simulate(
        grid,
        config,
        method;

        birth_nodes=birth_nodes,
        exact_solution=W,
        error_norm=norm
    )
    err = result.error
    println("Final error: $(Float64(round(err.ρ+err.p+err.u, digits=8)))")

    return nothing
end

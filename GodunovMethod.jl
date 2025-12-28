include("RiemannSolver.jl")

function find_Δt(CourantNumber::BigFloat, Δx::BigFloat, a::BigFloat, t::BigFloat)
    Δt = CourantNumber * Δx / a
    numOfSteps = div(t, Δt, RoundUp)
    return t / numOfSteps
end

function ρupE(W::PrimVars, γ::BigFloat)
    ρ = W.ρ
    u = W.u
    p = W.p
    E = p / (γ - 1) + 1 / 2 * ρ * u ^ 2
    return ρ, u, p, E
end

function GodunovMethodOneStep(
                              Ws::Vector{PrimVars},
                              xs::Vector{BigFloat},
                              Δx::BigFloat,
                              Δt::BigFloat,
                              γ::BigFloat,
                              ε::BigFloat
                             )
    Ws_n = Vector{PrimVars}()
    W₁ = W⁺(Ws[1], Ws[2], γ, ε)(BigFloat(0))
    for x in 2:(length(xs)+1)
        Wₗ = Ws[x]
        Wᵣ = Ws[x + 1]
        W₂ = W⁺(Wₗ, Wᵣ, γ, ε)(BigFloat(0))

        ρᴷ, uᴷ, _, Eᴷ = ρupE(Wₗ, γ)
        ρ₁, u₁, p₁, E₁ = ρupE(W₁, γ)
        ρ₂, u₂, p₂, E₂ = ρupE(W₂, γ)

        ρᴷ⁺¹ = ρᴷ - Δt / Δx * (ρ₂ * u₂ - ρ₁ * u₁)
        ρuᴷ⁺¹ = ρᴷ * uᴷ - Δt / Δx * (ρ₂ * u₂ ^ 2 + p₂ - ρ₁ * u₁ ^ 2 - p₁)
        Eᴷ⁺¹ = Eᴷ - Δt / Δx * ((E₂ + p₂) * u₂ - (E₁ + p₁) * u₁)

        uᴷ⁺¹ = ρuᴷ⁺¹ / ρᴷ⁺¹
        pᴷ⁺¹ = (γ - 1) * (Eᴷ⁺¹ - 1 / 2 * ρᴷ⁺¹ * uᴷ⁺¹ ^ 2)

        Wᴷ⁺¹ = PrimVars(ρᴷ⁺¹, uᴷ⁺¹, pᴷ⁺¹)
        push!(Ws_n, Wᴷ⁺¹)

        W₁ = W₂
    end
    return Ws_n
end

function GodunovMethod(xs::Vector{BigFloat}, γ::BigFloat, t_max::BigFloat, ε::BigFloat,
                       IC::Function, BCₗ::Function, BCᵣ::Function)
    a_max = 0
    for x in xs
        W = IC(x)
        a = sound_speed(W, γ)
        a_max = max(a, a_max)
    end
    Δx = xs[2] - xs[1]
    Δt = find_Δt(BigFloat("0.125"), Δx, a_max, t_max)
    Ws = IC.(xs)
    pushfirst!(Ws, BCₗ(BigFloat(0)))
    push!(Ws, BCᵣ(BigFloat(0)))
    for t in BigFloat("0"):Δt:t_max
        Ws = GodunovMethodOneStep(Ws, xs, Δx, Δt, γ, ε)
        pushfirst!(Ws, BCₗ(t))
        push!(Ws, BCᵣ(t))
    end
    pop!(Ws)
    popfirst!(Ws)
    return Ws
end

function MovingGridsOneStep(
                            Ws::Vector{PrimVars}, # length n + 1
                            ψs::Vector{BigFloat}, # length n - 1
                            xs::Vector{BigFloat}, # length n
                            Δy::BigFloat,
                            Δt::BigFloat,
                            γ::BigFloat,
                            ε::BigFloat,
                            ε₁::BigFloat,
                            ε₂::BigFloat,
                            α::BigFloat
                           )
    Qs = zeros(BigFloat, length(xs))
    for i in 2:(length(Ws)-2)
        Wₗ = Ws[i]
        Wᵣ = Ws[i + 1]
        Qs[i] = findQ(Wₗ, Wᵣ, γ, ε, ε₁, ε₂) + α / Δy * (1 / ψs[i] - 1 / ψs[i - 1])
        xs[i] = xs[i] - Qs[i] * Δt
        if xs[i] < xs[i - 1]
            println("   $i: Grid degenerated")
        end
    end
    ψs_n = zeros(BigFloat, length(ψs))
    for i in 1:length(ψs)
        ψs_n[i] = ψs[i] - Δt / Δy * (Qs[i + 1] - Qs[i])
    end
    Ws_n = Vector{PrimVars}()
    
    W₁ = W⁺(Ws[1], Ws[2], γ, ε)(- Qs[1])
    for i in 2:(length(Ws)-1)
        Wₗ = Ws[i]
        Wᵣ = Ws[i + 1]
        W₂ = W⁺(Wₗ, Wᵣ, γ, ε)(- Qs[i])

        ρᴷ, uᴷ, _, Eᴷ = ρupE(Wₗ, γ)
        ρ₁, u₁, p₁, E₁ = ρupE(W₁, γ)
        ρ₂, u₂, p₂, E₂ = ρupE(W₂, γ)

        ρᴷ⁺¹ = (
            ψs[i - 1] * ρᴷ - Δt / Δy * (
                Qs[i] * ρ₂ + ρ₂ * u₂ - 
                Qs[i - 1] * ρ₁ - ρ₁ * u₁
            )
        ) / ψs_n[i - 1]
        ρuᴷ⁺¹ = (
            ψs[i - 1] * ρᴷ * uᴷ - Δt / Δy * (
                Qs[i] * ρ₂ * u₂ + ρ₂ * u₂ ^ 2 + p₂ - 
                Qs[i - 1] * ρ₁ * u₁ -  ρ₁ * u₁ ^ 2 - p₁
            )
        ) / ψs_n[i - 1]
        Eᴷ⁺¹ = (
            ψs[i - 1] * Eᴷ - Δt / Δy * (
                Qs[i] * E₂ + (E₂ + p₂) * u₂ - 
                Qs[i - 1] * E₁ - (E₁ + p₁) * u₁
            )
        ) / ψs_n[i - 1]

        uᴷ⁺¹ = ρuᴷ⁺¹ / ρᴷ⁺¹
        pᴷ⁺¹ = (γ - 1) * (Eᴷ⁺¹ - 1 / 2 * ρᴷ⁺¹ * uᴷ⁺¹ ^ 2)

        Wᴷ⁺¹ = PrimVars(ρᴷ⁺¹, uᴷ⁺¹, pᴷ⁺¹)
        push!(Ws_n, Wᴷ⁺¹)

        W₁ = W₂
    end
    return Ws_n, ψs_n, xs
end

function MovingGrids(xs::Vector{BigFloat}, γ::BigFloat, t_max::BigFloat, ε::BigFloat,
                     ε₁::BigFloat, ε₂::BigFloat, α::BigFloat, CourantNumber::BigFloat,
                     IC::Function, BCₗ::Function, BCᵣ::Function)
    a_max = 0
    for x in xs
        W = IC(x)
        a = sound_speed(W, γ)
        a_max = max(a, a_max)
    end
    Δx = xs[2] - xs[1]
    Δt = find_Δt(CourantNumber, Δx, a_max, t_max)
    Ws = Vector{PrimVars}()
    for i in 1:(length(xs)-1)
        push!(Ws, IC((xs[i] + xs[i+1]) / 2))
    end
    ψs = ones(BigFloat, length(Ws))

    pushfirst!(Ws, BCₗ(BigFloat(0)))
    push!(Ws, BCᵣ(BigFloat(0)))
    for t in BigFloat("0"):Δt:t_max
        println("t = $t")
        Ws, ψs, xs = MovingGridsOneStep(Ws, ψs, xs, Δx, Δt, γ, ε, ε₁, ε₂, α)
        pushfirst!(Ws, BCₗ(t))
        push!(Ws, BCᵣ(t))
    end
    pop!(Ws)
    popfirst!(Ws)
    return Ws, xs
end

function MovingGridsDynamicStep(xs::Vector{BigFloat}, γ::BigFloat, t_max::BigFloat, ε::BigFloat,
                     ε₁::BigFloat, ε₂::BigFloat, α::BigFloat, CourantNumber::BigFloat,
                     IC::Function, BCₗ::Function, BCᵣ::Function)
    Δt = t_max
    t = BigFloat(0)
    Ws = Vector{PrimVars}()
    for i in 1:(length(xs)-1)
        push!(Ws, IC((xs[i] + xs[i+1]) / 2))
    end
    ψs = ones(BigFloat, length(Ws))
    Δx = xs[2] - xs[1]

    for i in 1:(length(xs)-1)
        h = xs[i+1] - xs[i]
        a = sound_speed(Ws[i], γ)
        Δt = min(Δt, CourantNumber * h / (a + abs(Ws[i].u)))
    end

    pushfirst!(Ws, BCₗ(BigFloat(0)))
    push!(Ws, BCᵣ(BigFloat(0)))
    while t < t_max
        Ws, ψs, xs = MovingGridsOneStep(Ws, ψs, xs, Δx, Δt, γ, ε, ε₁, ε₂, α)
        t += Δt

        Δt = t_max - t
        for i in 1:(length(xs)-1)
            h = xs[i+1] - xs[i]
            a = sound_speed(Ws[i], γ)
            Δt = min(Δt, CourantNumber * h / (a + abs(Ws[i].u)))
        end

        pushfirst!(Ws, BCₗ(t))
        push!(Ws, BCᵣ(t))
    end
    pop!(Ws)
    popfirst!(Ws)
    return Ws, ψs, xs
end

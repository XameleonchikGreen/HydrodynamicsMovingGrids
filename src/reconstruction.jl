function limiter(σ₁, σ₂)
    r = iszero(σ₂) ? BigFloat("0") : σ₁ / σ₂
    return max(0, min(1, 2r), min(2, r)) * σ₂
end

function limited_slopes(values::AbstractVector)
    σ = diff(values)

    pushfirst!(σ, zero(eltype(σ)))
    push!(σ, zero(eltype(σ)))

    return limiter.(σ[1:end-1], σ[2:end])
end

function muscl_reconstruct(values::AbstractVector)
    σ = limited_slopes(values)

    left = @. values - σ / 2
    right = @. values + σ / 2
    
    return left, right
end

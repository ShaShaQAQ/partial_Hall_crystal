using LinearAlgebra

struct ResponseLanczosKernel
    alpha::Vector{Float64}
    beta::Vector{Float64}
    source_norm2::Float64
    breakdown::Bool
end


function projected_response_source(
        J::AbstractMatrix,
        ground_state::AbstractVector;
        ground_subspace::Union{Nothing,AbstractMatrix}=nothing,
        atol::Float64=1e-12)
    g = ComplexF64.(ground_state)
    f = ComplexF64.(J * g)
    states = ground_subspace === nothing ? reshape(g, :, 1) : ground_subspace
    for state in eachcol(states)
        abs(norm(state) - 1) <= 1e-10 ||
            throw(ArgumentError("ground projector states must be normalized"))
        f .-= state .* dot(state, f)
    end
    norm(f) > atol ||
        throw(ArgumentError("projected current source has zero norm"))
    return f
end


function response_lanczos(
        H,
        E0::Real,
        source::AbstractVector;
        mmax::Int=length(source),
        breakdown_tol::Float64=1e-13)
    n = length(source)
    1 <= mmax <= n ||
        throw(ArgumentError("mmax must satisfy 1 <= mmax <= length(source)"))
    source_norm = norm(source)
    source_norm > breakdown_tol ||
        throw(ArgumentError("response source has zero norm"))

    q_prev = zeros(ComplexF64, n)
    q = ComplexF64.(source) ./ source_norm
    w = similar(q)
    alpha = Float64[]
    beta = Float64[]
    beta_prev = 0.0
    did_breakdown = false

    for step in 1:mmax
        mul!(w, H, q)
        @. w = w - E0 * q - beta_prev * q_prev
        a = real(dot(q, w))
        push!(alpha, a)
        @. w = w - a * q

        # Remove roundoff components along the local recurrence vectors.
        w .-= q .* dot(q, w)
        if step > 1
            w .-= q_prev .* dot(q_prev, w)
        end
        b = norm(w)

        if b <= breakdown_tol
            did_breakdown = true
            break
        end
        step == mmax && break

        push!(beta, b)
        q_prev, q, w = q, w, q_prev
        q ./= b
        beta_prev = b
    end

    return ResponseLanczosKernel(alpha, beta, source_norm^2, did_breakdown)
end


function lanczos_resolvent(kernel::ResponseLanczosKernel, z::Number)
    isempty(kernel.alpha) && throw(ArgumentError("empty Lanczos kernel"))
    denominator = complex(z) - kernel.alpha[end]
    for j in length(kernel.alpha)-1:-1:1
        denominator = complex(z) - kernel.alpha[j] -
                      kernel.beta[j]^2 / denominator
    end
    return kernel.source_norm2 / denominator
end


function truncate_kernel(kernel::ResponseLanczosKernel, m::Int)
    1 <= m <= length(kernel.alpha) ||
        throw(ArgumentError("invalid truncated dimension"))
    beta = m == 1 ? Float64[] : copy(kernel.beta[1:m-1])
    return ResponseLanczosKernel(
        copy(kernel.alpha[1:m]), beta, kernel.source_norm2,
        kernel.breakdown && m == length(kernel.alpha))
end

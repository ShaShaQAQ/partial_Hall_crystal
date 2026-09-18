using LinearAlgebra
using SparseArrays

"""Row-oriented sparse matrix with disjoint rows for threaded multiplication."""
struct ThreadedCSR
    nrows::Int
    ncols::Int
    rowptr::Vector{Int32}
    colind::Vector{Int32}
    nzval::Vector{ComplexF64}
end

Base.size(matrix::ThreadedCSR) = (matrix.nrows, matrix.ncols)
SparseArrays.nnz(matrix::ThreadedCSR) = length(matrix.nzval)

function threaded_csr(
        matrix::SparseMatrixCSC{ComplexF64,I}) where {I<:Integer}
    nrows, ncols = size(matrix)
    row_counts = [zeros(Int32, nrows) for _ in 1:Threads.nthreads()]
    Threads.@threads :static for column in 1:ncols
        counts = row_counts[Threads.threadid()]
        first_index = Int(matrix.colptr[column])
        last_index = Int(matrix.colptr[column + 1]) - 1
        @inbounds for p in first_index:last_index
            counts[Int(matrix.rowval[p])] += 1
        end
    end

    rowptr = Vector{Int32}(undef, nrows + 1)
    rowptr[1] = Int32(1)
    @inbounds for row in 1:nrows
        count = sum(Int(counts[row]) for counts in row_counts)
        rowptr[row + 1] = rowptr[row] + Int32(count)
    end

    colind = Vector{Int32}(undef, nnz(matrix))
    nzval = Vector{ComplexF64}(undef, nnz(matrix))
    cursor = copy(rowptr[1:nrows])
    @inbounds for column in 1:ncols
        first_index = Int(matrix.colptr[column])
        last_index = Int(matrix.colptr[column + 1]) - 1
        for p in first_index:last_index
            row = Int(matrix.rowval[p])
            destination = Int(cursor[row])
            colind[destination] = Int32(column)
            nzval[destination] = matrix.nzval[p]
            cursor[row] += 1
        end
    end
    return ThreadedCSR(nrows, ncols, rowptr, colind, nzval)
end

function LinearAlgebra.mul!(
        output::Vector{ComplexF64},
        matrix::ThreadedCSR,
        input::Vector{ComplexF64})
    length(input) == matrix.ncols ||
        throw(DimensionMismatch("input has wrong length"))
    length(output) == matrix.nrows ||
        throw(DimensionMismatch("output has wrong length"))
    Threads.@threads :static for row in 1:matrix.nrows
        value = 0.0 + 0.0im
        first_index = Int(matrix.rowptr[row])
        last_index = Int(matrix.rowptr[row + 1]) - 1
        @inbounds for p in first_index:last_index
            value += matrix.nzval[p] * input[Int(matrix.colind[p])]
        end
        output[row] = value
    end
    return output
end

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
        breakdown_tol::Float64=1e-13,
        ground_state::Union{Nothing,AbstractVector}=nothing,
        progress::Function=(args...)->nothing)
    n = length(source)
    1 <= mmax <= n ||
        throw(ArgumentError("mmax must satisfy 1 <= mmax <= length(source)"))
    source_norm = norm(source)
    source_norm > breakdown_tol ||
        throw(ArgumentError("response source has zero norm"))
    ground_state === nothing ||
        (length(ground_state) == n &&
         abs(norm(ground_state) - 1) <= 1e-8) ||
        throw(ArgumentError(
            "ground_state must be normalized and match source"))

    q_prev = zeros(ComplexF64, n)
    q = ComplexF64.(source) ./ source_norm
    w = similar(q)
    alpha = Float64[]
    beta = Float64[]
    beta_prev = 0.0
    did_breakdown = false

    for step in 1:mmax
        mul!(w, H, q)
        if ground_state !== nothing
            w .-= ground_state .* dot(ground_state, w)
        end
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
        progress(step, alpha, beta, q_prev, q, beta_prev)
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


struct ExactSpectralData
    deltas::Vector{Float64}
    weights::Vector{Float64}
end


function exact_spectral_data(
        eigenvalues::AbstractVector,
        eigenvectors::AbstractMatrix,
        E0::Real,
        source::AbstractVector;
        gap_tol::Float64=1e-10)
    deltas_all = real.(eigenvalues .- E0)
    amplitudes = eigenvectors' * source
    keep = findall(delta -> delta > gap_tol, deltas_all)
    return ExactSpectralData(deltas_all[keep], abs2.(amplitudes[keep]))
end


function exact_resolvent(data::ExactSpectralData, z::Number)
    return sum(data.weights ./ (complex(z) .- data.deltas))
end


function system_area(lat::GenLat)
    a1 = collect(lat.a1)
    a2 = collect(lat.a2)
    T1 = lat.T1[1] .* a1 .+ lat.T1[2] .* a2
    T2 = lat.T2[1] .* a1 .+ lat.T2[2] .* a2
    return abs(T1[1] * T2[2] - T1[2] * T2[1])
end


function lanczos_drude_weight(kernel::ResponseLanczosKernel, Kexp::Real)
    return Float64(Kexp + 2real(lanczos_resolvent(kernel, 0.0 + 0.0im)))
end


function conductivity_parts(
        G,
        Kexp::Real,
        area::Real,
        omega::Real,
        eta::Real)
    area > 0 || throw(ArgumentError("system area must be positive"))
    eta > 0 || throw(ArgumentError("eta must be positive"))
    z = complex(omega, eta)
    chi = G(z) + G(-z)
    chi0 = 2real(G(0.0 + 0.0im))
    D = Float64(Kexp + chi0)
    prefactor = 2pi * im / area
    drude = prefactor * D / z
    regular = prefactor * (chi - chi0) / z
    return (total=drude + regular,
            drude=drude,
            regular=regular,
            drude_weight=D)
end


function scaled_max_error(
        actual::AbstractVector,
        reference::AbstractVector)
    length(actual) == length(reference) ||
        throw(DimensionMismatch("curve lengths differ"))
    scale = max(1.0, maximum(abs, reference))
    return maximum(abs, actual .- reference) / scale
end

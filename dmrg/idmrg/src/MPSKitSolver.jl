struct MPSKitSolverStageRecord
    stage::Int
    requested_maxdim::Int
    actual_maxdim::Int
    cutoff::Float64
    idmrg_maxiter::Int
    vumps_maxiter::Int
    energy_per_site::Float64
    energy_imaginary::Float64
    idmrg_diagnostic::Float64
    vumps_residual::Float64
    recomputed_galerkin_residual::Float64
    elapsed_seconds::Float64
end

struct MPSKitSolverResult{S,E}
    state::S
    environments::E
    records::Vector{MPSKitSolverStageRecord}
    energy_per_site::Float64
    energy_imaginary::Float64
    galerkin_residual::Float64
    recomputed_galerkin_residual::Float64
    link_dimensions::Vector{Int}
    converged::Bool
end

function mpskit_solver_converged(;
    idmrg_diagnostic::Real,
    recomputed_galerkin_residual::Real,
    energy_per_site::Real,
    energy_imaginary::Real,
    final_stage_reached::Bool,
    galerkin_tol::Real,
    energy_imag_tol::Real,
)
    _ = idmrg_diagnostic
    return final_stage_reached &&
        isfinite(recomputed_galerkin_residual) &&
        recomputed_galerkin_residual <= galerkin_tol &&
        isfinite(energy_per_site) &&
        isfinite(energy_imaginary) &&
        abs(energy_imaginary) <= energy_imag_tol
end

function _mpskit_solver_link_dimensions(state)
    return [
        TensorKit.dim(MPSKit.right_virtualspace(state, site)) for
        site in 1:length(state)
    ]
end

function _mpskit_validate_solver_settings(
    maxdim_schedule,
    cutoff::Real,
    idmrg_maxiter::Int,
    vumps_maxiter::Int,
    galerkin_tol::Real,
    energy_imag_tol::Real,
)
    schedule = Int[dimension for dimension in maxdim_schedule]
    isempty(schedule) && throw(ArgumentError("maxdim_schedule cannot be empty"))
    all(>(0), schedule) || throw(
        ArgumentError("maxdim_schedule entries must be positive")
    )
    all(schedule[index] < schedule[index + 1] for index in 1:(length(schedule) - 1)) || throw(
        ArgumentError("maxdim_schedule must be strictly increasing")
    )
    isfinite(cutoff) && cutoff > 0 || throw(
        ArgumentError("cutoff must be finite and positive")
    )
    idmrg_maxiter > 0 || throw(ArgumentError("idmrg_maxiter must be positive"))
    vumps_maxiter > 0 || throw(ArgumentError("vumps_maxiter must be positive"))
    isfinite(galerkin_tol) && galerkin_tol > 0 || throw(
        ArgumentError("galerkin_tol must be finite and positive")
    )
    isfinite(energy_imag_tol) && energy_imag_tol >= 0 || throw(
        ArgumentError("energy_imag_tol must be finite and nonnegative")
    )
    return schedule
end

function run_mpskit_idmrg(
    hamiltonian,
    initial_state;
    maxdim_schedule,
    cutoff::Real,
    idmrg_maxiter::Int,
    vumps_maxiter::Int,
    galerkin_tol::Real,
    energy_imag_tol::Real=1e-10,
    verbosity::Int=0,
    progress_callback=(state, environments, records) -> nothing,
)
    length(initial_state) == length(hamiltonian) || throw(
        ArgumentError("state and Hamiltonian unit-cell lengths do not match")
    )
    schedule = _mpskit_validate_solver_settings(
        maxdim_schedule,
        cutoff,
        idmrg_maxiter,
        vumps_maxiter,
        galerkin_tol,
        energy_imag_tol,
    )

    state = initial_state
    final_environments = MPSKit.environments(state, hamiltonian, state)
    records = MPSKitSolverStageRecord[]
    local_eigsolve = KrylovKit.Lanczos(;
        tol=min(Float64(cutoff), Float64(galerkin_tol)),
        maxiter=200,
        eager=true,
        krylovdim=30,
        verbosity=0,
    )
    for (stage, requested_maxdim) in enumerate(schedule)
        started = time_ns()
        idmrg_environments = MPSKit.environments(state, hamiltonian, state)
        idmrg_algorithm = MPSKit.IDMRG2(;
            tol=Float64(cutoff),
            maxiter=idmrg_maxiter,
            verbosity,
            trunc=MatrixAlgebraKit.truncrank(requested_maxdim) &
                MatrixAlgebraKit.trunctol(; atol=Float64(cutoff)),
            alg_eigsolve=local_eigsolve,
        )
        state, _, idmrg_diagnostic = MPSKit.find_groundstate(
            state,
            hamiltonian,
            idmrg_algorithm,
            idmrg_environments,
        )

        vumps_environments = MPSKit.environments(state, hamiltonian, state)
        vumps_algorithm = MPSKit.VUMPS(;
            tol=Float64(galerkin_tol),
            maxiter=vumps_maxiter,
            verbosity,
            alg_eigsolve=local_eigsolve,
        )
        state, _, vumps_residual = MPSKit.find_groundstate(
            state,
            hamiltonian,
            vumps_algorithm,
            vumps_environments,
        )

        final_environments = MPSKit.environments(state, hamiltonian, state)
        recomputed_residual = MPSKit.calc_galerkin(
            state,
            hamiltonian,
            state,
            final_environments,
        )
        energy = ComplexF64(
            MPSKit.expectation_value(state, hamiltonian, final_environments) /
            length(state),
        )
        link_dimensions = _mpskit_solver_link_dimensions(state)
        elapsed_seconds = (time_ns() - started) / 1e9
        push!(
            records,
            MPSKitSolverStageRecord(
                stage,
                requested_maxdim,
                maximum(link_dimensions),
                Float64(cutoff),
                idmrg_maxiter,
                vumps_maxiter,
                real(energy),
                imag(energy),
                Float64(idmrg_diagnostic),
                Float64(vumps_residual),
                Float64(recomputed_residual),
                elapsed_seconds,
            ),
        )
        progress_callback(state, final_environments, records)
    end

    final_record = last(records)
    link_dimensions = _mpskit_solver_link_dimensions(state)
    converged = mpskit_solver_converged(;
        idmrg_diagnostic=final_record.idmrg_diagnostic,
        recomputed_galerkin_residual=final_record.recomputed_galerkin_residual,
        energy_per_site=final_record.energy_per_site,
        energy_imaginary=final_record.energy_imaginary,
        final_stage_reached=length(records) == length(schedule),
        galerkin_tol,
        energy_imag_tol,
    )
    return MPSKitSolverResult(
        state,
        final_environments,
        records,
        final_record.energy_per_site,
        final_record.energy_imaginary,
        final_record.vumps_residual,
        final_record.recomputed_galerkin_residual,
        link_dimensions,
        converged,
    )
end

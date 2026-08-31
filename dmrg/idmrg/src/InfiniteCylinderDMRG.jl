module InfiniteCylinderDMRG

include("../../CylinderModelCore.jl")
using .CylinderModelCore:
    CylinderModelParams,
    NB_DISPS,
    reciprocal_vectors,
    hopping_matrices_from_get_Hk

using ITensors
using ITensorMPS
using ITensorInfiniteMPS
using KrylovKit
using HDF5
using TOML
using Random
using SHA
using SparseArrays
using LinearAlgebra

import JLD2
import BlockTensorKit
import MatrixAlgebraKit
import MPSKit
import TensorKit
import TensorKitTensors

include("Configuration.jl")
include("PaperGeometry.jl")
include("InitialCandidates.jl")
include("TransverseTranslation.jl")
include("Hamiltonian.jl")
include("VUMPSRunner.jl")
include("Observables.jl")
include("MomentumEntanglement.jl")
include("BranchTracking.jl")
include("Output.jl")
include("CLI.jl")
include("Fig2Benchmark.jl")
include("MPSKitBackend.jl")
include("MPSKitHamiltonian.jl")

export CylinderModelParams,
    NB_DISPS,
    reciprocal_vectors,
    hopping_matrices_from_get_Hk,
    InfiniteCylinderConfig,
    SUPPORTED_GEOMETRIES,
    circumference_unit_cells,
    physical_circumference_sites,
    sites_per_cell,
    unit_cells_per_cell,
    particles_per_cell,
    physical_site_density,
    charge_scale,
    linear_site,
    site_coordinates,
    canonical_site,
    paper_site,
    paper_coordinates,
    with_flux,
    configuration_signature,
    default_occupied_sites,
    PAPER_A1,
    PAPER_A2,
    PAPER_SUBLATTICE_POSITIONS,
    PaperHoppingOrbit,
    paper_triangular_coordinates,
    paper_physical_position,
    paper_periodic_distance_squared,
    paper_hopping_orbits,
    paper_orbits_are_hermitian,
    paper_bloch_from_orbits,
    paper_reference_Hk,
    InitialCandidate,
    paper_cdw_color,
    fig2_initial_candidates,
    translated_fock_basis,
    fermionic_ring_translation,
    apply_fermionic_ring_translation,
    InfiniteHoppingTerm,
    InfiniteInteractionTerm,
    normalize_directed,
    normalize_undirected,
    build_infinite_hoppings,
    build_infinite_interactions,
    build_infinite_model_terms,
    hamiltonian_is_hermitian,
    canonical_term_dict,
    tile_to_finite_window,
    finite_term_parity,
    neighbor_shell_counts,
    straight_seam_is_connected,
    paper_window_parity,
    infinite_opsum,
    build_infinite_mpo,
    initial_infinite_mps,
    centered_site_charges,
    validate_centered_site_charges,
    link_dimensions,
    vumps_converged,
    VUMPSRecord,
    VUMPSProgressEvent,
    SubspaceExpansionRecord,
    VUMPSResult,
    unit_cell_energy,
    vumps_iteration,
    expand_subspace,
    run_vumps,
    EnergyData,
    DensityRow,
    EntanglementLevel,
    EntanglementSector,
    EntanglementData,
    NeutralTransferData,
    normalize_energy,
    energy_data,
    density_data,
    entanglement_data,
    neutral_transfer_data,
    MomentumValidation,
    MomentumEntanglementLevel,
    MomentumEntanglementData,
    validate_momentum_blocks,
    momentum_entanglement_data,
    bhattacharyya_distance,
    best_sector_shift,
    charge_matched_spectrum_distance,
    branch_flags,
    MixedTransferFidelity,
    mixed_transfer_fidelity,
    BranchCandidate,
    CandidateSelection,
    select_candidate,
    InfiniteCanonicalMPS,
    CheckpointFormatError,
    CheckpointCompatibilityError,
    save_checkpoint,
    load_checkpoint,
    write_output_files,
    SinglePointSettings,
    parse_single_point_args,
    FluxScanSettings,
    parse_flux_scan_args,
    flux_grid,
    FluxCandidateData,
    select_flux_candidate_data,
    flux_scan_row,
    flux_branch_event,
    WorkflowValidationError,
    SinglePointOperations,
    SinglePointRunResult,
    configure_cli_threads,
    run_single_point,
    SCAN_SUMMARY_HEADER,
    BRANCH_EVENTS_HEADER,
    SECTOR_GAUGE_HEADER,
    FluxScanRunResult,
    run_flux_scan,
    Fig2BenchmarkSpec,
    Fig2CandidateEvidence,
    Fig2BenchmarkOperations,
    Fig2Selection,
    Fig2BenchmarkRun,
    Fig2BenchmarkSettings,
    parse_fig2_benchmark_args,
    load_fig2_benchmark,
    run_fig2_benchmark,
    write_fig2_acceptance_report!

export BLOCKTENSORKIT_COMMIT,
    MPSKIT_BACKEND_COMMIT,
    TENSORKITTENSORS_COMMIT,
    mpskit_backend_provenance,
    mpskit_centered_charges,
    mpskit_total_raw_charge,
    mpskit_physical_spaces,
    mpskit_physical_raw_charges,
    mpskit_product_state,
    mpskit_product_occupations,
    mpskit_finite_hamiltonian,
    mpskit_dense_matrix

end

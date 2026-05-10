"""Measurement-based-QC routing demo on a 2D cluster state.

The 2D cluster state is the existing Genuine-Multipartite-Entanglement resource.
Beyond certifying it as a global entangled state via the stabilizer witness,
this module *uses* it: by measuring off-path qubits in Z and internal-path
qubits in X, we localize entanglement onto two chosen distant endpoints A and
B. The endpoint pair is then verified by Bell fidelity (F > 1/2) and
optionally by a CHSH game (S > 2).
"""

from .cluster_routing import (
    build_cluster_routing_circuit,
    build_routing_bell_circuits,
    build_routing_chsh_circuits,
    byproduct_signs_from_outcome,
    compute_bell_fidelity,
    compute_chsh_game,
    corrected_endpoint_correlation,
    derive_byproduct_table_simulator,
    grid_edges,
    grid_neighbors,
    run_routing_demo,
    save_results_json,
    sweep_routes,
    uncorrected_endpoint_correlation,
    validate_path_is_connected,
)

__all__ = [
    "build_cluster_routing_circuit",
    "build_routing_bell_circuits",
    "build_routing_chsh_circuits",
    "byproduct_signs_from_outcome",
    "compute_bell_fidelity",
    "compute_chsh_game",
    "corrected_endpoint_correlation",
    "derive_byproduct_table_simulator",
    "grid_edges",
    "grid_neighbors",
    "run_routing_demo",
    "save_results_json",
    "sweep_routes",
    "uncorrected_endpoint_correlation",
    "validate_path_is_connected",
]

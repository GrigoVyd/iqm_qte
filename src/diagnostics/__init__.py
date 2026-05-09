"""Hardware-diagnostic layer: experimental entanglement-aware audit of qubit
selectors and physical couplers."""

from .edge_bell_map import (
    EdgeResult,
    audit_layout,
    build_matching_graph_state_circuit,
    compute_edge_fidelity,
    edge_color_matchings,
    edge_residuals,
    expectation_from_counts,
    get_coupling_edges,
    isolated_vs_parallel,
    make_undirected_edges,
    node_diagnostics,
    plot_edge_map,
    recommend_best_chain,
    recommend_best_patch,
    run_edge_map,
    run_isolated_edge_map,
    save_results_json_csv,
)

__all__ = [
    "EdgeResult",
    "audit_layout",
    "build_matching_graph_state_circuit",
    "compute_edge_fidelity",
    "edge_color_matchings",
    "edge_residuals",
    "expectation_from_counts",
    "get_coupling_edges",
    "isolated_vs_parallel",
    "make_undirected_edges",
    "node_diagnostics",
    "plot_edge_map",
    "recommend_best_chain",
    "recommend_best_patch",
    "run_edge_map",
    "run_isolated_edge_map",
    "save_results_json_csv",
]

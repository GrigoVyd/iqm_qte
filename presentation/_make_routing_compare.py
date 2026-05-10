"""Make a side-by-side chip-layout comparison of the IQM Selector chain
vs our beam-search chain at n=15 on Garnet. Saves to
presentation/assets/results/routing_chains_garnet.png."""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib.pyplot as plt
import os as _os
_os.environ.pop("IQM_TOKEN", None)  # use Aer fallback so get_backend works offline
from src.backend import get_backend
from src.visualization import plot_device_topology, GARNET_POSITIONS

# We can't get_backend(device='garnet') without a token, so build a minimal
# stand-in: we only need num_qubits + coupling_map for the visualizer. Use a
# tiny shim.
class _ShimBackend:
    num_qubits = 20
    @property
    def coupling_map(self):
        # Edges as a list of pairs (both directions OK; visualizer dedups)
        edges = [
            (0,1),(0,3),(1,4),(2,3),(2,7),(3,4),(3,8),(4,5),(4,9),
            (5,6),(5,10),(6,11),(7,8),(7,12),(8,9),(8,13),(9,10),(9,14),
            (10,11),(10,15),(11,16),(12,13),(13,14),(13,17),(14,15),(14,18),
            (15,16),(15,19),(16,19),(17,18),(18,19),
        ]
        # CouplingMap-like: iterable of pairs
        return [tuple(e) for e in edges]
    def __repr__(self):
        return "Garnet(stub)"

backend = _ShimBackend()
data = json.load(open(os.path.join(os.path.dirname(__file__), "..",
                                     "experiments", "consolidated_results",
                                     "routing_n15_garnet.json")))
sel_layout = data["selector"]["layout"]
beam_layout = data["beam"]["layout"]

pos = {q: tuple(map(float, GARNET_POSITIONS[q])) for q in GARNET_POSITIONS}

fig, axes = plt.subplots(1, 2, figsize=(15, 7))

for ax, layout, title, F_z in zip(
        axes,
        [sel_layout, beam_layout],
        ["IQM Qubit Selector layout", "Our beam-search Hamiltonian-path layout"],
        [0.31, 0.68],
    ):
    edges = [(layout[i], layout[i + 1]) for i in range(len(layout) - 1)]
    plot_device_topology(
        backend, metrics=None, cz_fidelities=None,
        color_by="bipartite",
        highlight_qubits=layout, highlight_edges=edges,
        title=f"{title}\n(n=15 W-state on Garnet · F_z = {F_z:.2f})",
        ax=ax, pos=pos, show_labels=True, spotlight=True,
    )

fig.suptitle("Same circuit, two routings — same hardware, very different fidelity",
              fontsize=15, fontweight="bold", y=1.005)
plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "assets", "results",
                    "routing_chains_garnet.png")
plt.savefig(out, dpi=130, bbox_inches="tight")
print("saved", out)

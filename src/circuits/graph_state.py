"""Graph state circuit on an arbitrary connected bipartite graph.

A graph state on graph G = (V, E) is defined by:
    |G⟩ = ∏_{(i,j) ∈ E} CZ_{ij} · H^⊗n |0⟩^⊗n

For any connected bipartite (2-colorable) graph, this state is genuinely
multipartitely entangled and the witness W = Σ_i ⟨g_i⟩ where g_i = X_i ⊗ ⊗_{j∈N(i)} Z_j
satisfies W ≤ n−1 for any biseparable state (Tóth & Gühne 2005).

Compared to a rectangular cluster state on the same n qubits:
- A spanning tree has n−1 edges (vs ~2n for a 2D grid)
- Fewer CZ gates → higher prep fidelity → certifies more qubits
"""

from __future__ import annotations

from collections import defaultdict, deque

from qiskit import QuantumCircuit


def build_graph_state(n: int, edges: list[tuple[int, int]]) -> QuantumCircuit:
    """Build the graph-state preparation circuit on n logical qubits.

    edges: list of (i, j) pairs with 0 <= i, j < n. Self-loops and duplicates
    are tolerated (deduped).
    """
    qc = QuantumCircuit(n, name=f"GraphState-{n}-{len(edges)}e")
    qc.h(range(n))
    qc.barrier()
    seen: set[tuple[int, int]] = set()
    for a, b in edges:
        if a == b:
            continue
        e = (min(a, b), max(a, b))
        if e in seen:
            continue
        seen.add(e)
        qc.cz(a, b)
    return qc


def two_coloring(n: int, edges: list[tuple[int, int]]) -> list[int] | None:
    """BFS 2-coloring of the graph. Returns list of {0, 1} per node, or None
    if the graph is not bipartite. Disconnected components are handled.
    """
    adj: dict[int, list[int]] = defaultdict(list)
    for a, b in edges:
        if a == b:
            continue
        adj[a].append(b)
        adj[b].append(a)

    color = [-1] * n
    for start in range(n):
        if color[start] != -1:
            continue
        color[start] = 0
        q = deque([start])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if color[v] == -1:
                    color[v] = 1 - color[u]
                    q.append(v)
                elif color[v] == color[u]:
                    return None  # odd cycle → not bipartite
    return color


def neighbours_from_edges(n: int, edges: list[tuple[int, int]]) -> dict[int, list[int]]:
    """Adjacency list (sorted, deduped) from an edge list."""
    adj: dict[int, set[int]] = defaultdict(set)
    for a, b in edges:
        if a == b:
            continue
        adj[a].add(b)
        adj[b].add(a)
    return {q: sorted(adj.get(q, [])) for q in range(n)}

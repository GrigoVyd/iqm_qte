// chip.js — render Garnet's coupling map as an SVG chip layout.
// Coordinates transcribed from src/visualization.py:GARNET_POSITIONS so the
// website's chip drawing matches the matplotlib figures exactly.

const GARNET_POSITIONS = {
  19: [2, 6], 16: [4, 6],
  18: [1, 5], 15: [3, 5], 11: [5, 5],
  17: [0, 4], 14: [2, 4], 10: [4, 4],  6: [6, 4],
  13: [1, 3],  9: [3, 3],  5: [5, 3],
  12: [0, 2],  8: [2, 2],  4: [4, 2],
   7: [1, 1],  3: [3, 1],  1: [5, 1],
   2: [2, 0],  0: [4, 0],
};

// All native CZ edges on Garnet (derived from get_coupling_edges output, both directions deduplicated).
const GARNET_EDGES = [
  [0,1],[0,3],
  [1,4],
  [2,3],[2,7],
  [3,4],[3,8],
  [4,5],[4,9],
  [5,6],[5,10],
  [6,11],
  [7,8],[7,12],
  [8,9],[8,13],
  [9,10],[9,14],
  [10,11],[10,15],
  [11,16],
  [12,13],
  [13,14],[13,17],
  [14,15],[14,18],
  [15,16],[15,19],
  [16,19],   // Garnet has the (16,19) closure
  [17,18],
  [18,19],
];

// Render the chip into a host element. Returns handles for the animation step.
function renderGarnetChip(host, opts = {}) {
  const margin = 30;
  const cellSize = 64;
  const cols = 7, rows = 7;
  const W = cellSize * (cols - 1) + 2 * margin;
  const H = cellSize * (rows - 1) + 2 * margin;
  const xy = (q) => {
    const [cx, cy] = GARNET_POSITIONS[q];
    return [margin + cx * cellSize, margin + (rows - 1 - cy) * cellSize];
  };

  const svgNS = "http://www.w3.org/2000/svg";
  const svg = document.createElementNS(svgNS, "svg");
  svg.classList.add("chip-svg");
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  svg.setAttribute("xmlns", svgNS);

  // Background edges (all coupling-map edges)
  const edgeRefs = {};
  for (const [a, b] of GARNET_EDGES) {
    const [x1, y1] = xy(a);
    const [x2, y2] = xy(b);
    const e = document.createElementNS(svgNS, "line");
    e.setAttribute("x1", x1); e.setAttribute("y1", y1);
    e.setAttribute("x2", x2); e.setAttribute("y2", y2);
    e.classList.add("edge");
    svg.appendChild(e);
    edgeRefs[edgeKey(a, b)] = e;
  }

  // Tree-edge overlays (drawn over the dim baseline edges, animated separately)
  const treeEdgeRefs = {};
  if (opts.treeEdges) {
    for (const [a, b] of opts.treeEdges) {
      const [x1, y1] = xy(a);
      const [x2, y2] = xy(b);
      const e = document.createElementNS(svgNS, "line");
      e.setAttribute("x1", x1); e.setAttribute("y1", y1);
      e.setAttribute("x2", x2); e.setAttribute("y2", y2);
      e.classList.add("edge", "tree");
      // Stroke-dash for draw-on animation
      const len = Math.hypot(x2 - x1, y2 - y1);
      e.setAttribute("stroke-dasharray", len);
      e.setAttribute("stroke-dashoffset", len);
      e.dataset.length = len;
      svg.appendChild(e);
      treeEdgeRefs[edgeKey(a, b)] = e;
    }
  }

  // Nodes
  const nodeRefs = {};
  const labelRefs = {};
  for (const q of Object.keys(GARNET_POSITIONS).map(Number)) {
    const [x, y] = xy(q);
    const n = document.createElementNS(svgNS, "circle");
    n.setAttribute("cx", x); n.setAttribute("cy", y);
    n.setAttribute("r", 18);
    n.classList.add("node");
    svg.appendChild(n);
    const t = document.createElementNS(svgNS, "text");
    t.setAttribute("x", x); t.setAttribute("y", y);
    t.classList.add("node-label");
    t.textContent = `Q${q + 1}`;
    svg.appendChild(t);
    nodeRefs[q] = n;
    labelRefs[q] = t;
  }

  host.appendChild(svg);
  return { svg, nodeRefs, labelRefs, edgeRefs, treeEdgeRefs };
}

function edgeKey(a, b) {
  return a < b ? `${a}_${b}` : `${b}_${a}`;
}

// Reset all tree-related state to the initial (un-grown) state.
function resetTree(handles) {
  for (const k in handles.treeEdgeRefs) {
    const e = handles.treeEdgeRefs[k];
    e.setAttribute("stroke-dashoffset", e.dataset.length);
    e.style.opacity = 0;
  }
  for (const q in handles.nodeRefs) {
    handles.nodeRefs[q].classList.remove("tree");
    handles.labelRefs[q].setAttribute("fill", "");
  }
}

window.QteChip = { renderGarnetChip, edgeKey, resetTree, GARNET_POSITIONS, GARNET_EDGES };

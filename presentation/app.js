// app.js — animations, dot nav, keyboard, fullscreen.

(function () {
  const deck = document.getElementById("deck");
  const sections = Array.from(deck.querySelectorAll("section.panel"));
  const dotsHost = document.getElementById("dots");
  const fsBtn = document.getElementById("fs");

  // ---------- dots ----------
  sections.forEach((sec, i) => {
    const d = document.createElement("button");
    d.className = "dot";
    d.title = `Slide ${i + 1}`;
    d.addEventListener("click", () => stepTo(i));
    dotsHost.appendChild(d);
  });
  const dots = Array.from(dotsHost.children);

  // ---------- IntersectionObserver: which panel is active ----------
  let currentIdx = 0;
  const triggered = new Set();
  const io = new IntersectionObserver(
    (entries) => {
      entries.forEach((e) => {
        if (e.isIntersecting && e.intersectionRatio > 0.55) {
          const idx = sections.indexOf(e.target);
          currentIdx = idx;
          dots.forEach((d, i) => d.classList.toggle("active", i === idx));
          if (!triggered.has(idx)) {
            triggered.add(idx);
            runAnimationFor(idx);
          }
        }
      });
    },
    { threshold: [0, 0.55, 1] }
  );
  sections.forEach((s) => io.observe(s));

  // ---------- step helpers ----------
  function stepTo(idx) {
    idx = Math.max(0, Math.min(sections.length - 1, idx));
    sections[idx].scrollIntoView({ behavior: "smooth", block: "start" });
  }

  // ---------- keyboard ----------
  document.addEventListener("keydown", (e) => {
    if (e.key === "ArrowRight" || e.key === "ArrowDown" || e.key === " " || e.key === "PageDown") {
      e.preventDefault(); stepTo(currentIdx + 1);
    } else if (e.key === "ArrowLeft" || e.key === "ArrowUp" || e.key === "PageUp") {
      e.preventDefault(); stepTo(currentIdx - 1);
    } else if (e.key === "Home") { e.preventDefault(); stepTo(0); }
    else if (e.key === "End") { e.preventDefault(); stepTo(sections.length - 1); }
    else if (e.key.toLowerCase() === "f") { toggleFullscreen(); }
  });

  // ---------- fullscreen ----------
  fsBtn.addEventListener("click", toggleFullscreen);
  function toggleFullscreen() {
    if (!document.fullscreenElement) {
      document.documentElement.requestFullscreen?.();
    } else {
      document.exitFullscreen?.();
    }
  }

  // ---------- per-panel animations ----------
  function runAnimationFor(idx) {
    const id = sections[idx].id;
    if (id === "s1-cover")        animateCover();
    else if (id === "s2-challenge")    animateChallenge();
    else if (id === "s3-story")        animateStory();
    else if (id === "s4-routing")      animateNumbers(sections[idx]);
    else if (id === "s5-w-sweep")      animateNumbers(sections[idx]);
    else if (id === "s6-tree")         animateTree();
    else if (id === "s7-mitigation")   animateMitigation();
    else if (id === "s8-bottleneck")   animateNumbers(sections[idx]);
    else if (id === "s9-mbqc")         animateMBQC();
    else if (id === "s10-summary")     animateSummary();
  }

  // Panel 1 — Cover: typewriter the subtitle
  function animateCover() {
    const tag = document.querySelector("#s1-cover .cover-tag");
    const sub = document.querySelector("#s1-cover .cover-sub");
    const typer = document.querySelector("#s1-cover .cover-typer");
    if (!typer) return;
    anime.timeline({ easing: "easeOutCubic" })
      .add({ targets: tag, opacity: [0, 1], translateY: [10, 0], duration: 600 })
      .add({ targets: sub, opacity: [0, 1], duration: 200 }, "+=300")
      .add({ targets: typer, width: ["0px", typer.scrollWidth + "px"], duration: 1400, easing: "steps(28)" });
  }

  // Panel 2 — Challenge: stagger cards, count percentages
  function animateChallenge() {
    const cards = document.querySelectorAll("#s2-challenge .card");
    anime({ targets: cards, opacity: [0, 1], translateY: [22, 0], duration: 700, delay: anime.stagger(120), easing: "easeOutQuart" });
    cards.forEach((c) => {
      const numEl = c.querySelector(".num");
      if (!numEl) return;
      const target = parseFloat(numEl.textContent);
      const obj = { v: 0 };
      anime({
        targets: obj, v: target, duration: 1400, easing: "easeOutQuart",
        update: () => { numEl.firstChild.nodeValue = Math.round(obj.v); }
      });
    });
  }

  // Panel 3 — Story cards stagger in
  function animateStory() {
    const cards = document.querySelectorAll("#s3-story .story-card");
    anime({
      targets: cards, opacity: [0, 1], translateX: [-30, 0],
      duration: 700, delay: anime.stagger(150), easing: "easeOutQuart"
    });
    anime({ targets: "#s3-story .fade-in", opacity: [0, 1], duration: 800, delay: 600 });
  }

  // Panels 4, 5, 8 — number rolls on .counter elements
  function animateNumbers(panel) {
    const counters = panel.querySelectorAll(".counter");
    counters.forEach((el) => {
      const end = parseFloat(el.dataset.end);
      const dec = parseInt(el.dataset.decimals || "0", 10);
      const obj = { v: 0 };
      anime({
        targets: obj, v: end, duration: 1200, easing: "easeOutCubic",
        update: () => { el.textContent = obj.v.toFixed(dec); }
      });
    });
    // Also reveal figures + table rows
    anime({ targets: panel.querySelectorAll(".figure"), opacity: [0, 1], duration: 700, easing: "easeOutQuad" });
    anime({ targets: panel.querySelectorAll("tr"), opacity: [0, 1], translateX: [-12, 0], delay: anime.stagger(70, { start: 200 }), duration: 500, easing: "easeOutQuad" });
  }

  // Panel 6 — Garnet spanning-tree growth (Prim's-style draw-on)
  function animateTree() {
    const host = document.getElementById("garnet-chip-host");
    if (host.dataset.rendered) return;
    host.dataset.rendered = "1";
    const data = JSON.parse(document.getElementById("data-tree-garnet").textContent);
    const handles = window.QteChip.renderGarnetChip(host, { treeEdges: data.edges_in_order });
    const tree = handles.treeEdgeRefs;
    const nodes = handles.nodeRefs;

    const wEl = document.getElementById("prim-W");
    const eEl = document.getElementById("prim-edges");
    const qEl = document.getElementById("prim-qubits");

    // Step through edges in order, growing the tree.
    const edges = data.edges_in_order;
    const qubitsAdded = new Set();
    qubitsAdded.add(edges[0][0]); // seed

    edges.forEach(([a, b], i) => {
      const e = tree[window.QteChip.edgeKey(a, b)];
      if (!e) return;
      const len = parseFloat(e.dataset.length);
      const startMs = 600 + i * 320;
      // edge draw-on
      anime({
        targets: e,
        strokeDashoffset: [len, 0],
        opacity: [0, 1],
        duration: 380,
        easing: "easeOutQuad",
        delay: startMs,
      });
      // grow nodes on either side
      [a, b].forEach((q) => {
        if (qubitsAdded.has(q)) return;
        qubitsAdded.add(q);
        anime({
          targets: nodes[q],
          duration: 280,
          delay: startMs + 100,
          easing: "easeOutBack",
          begin: () => nodes[q].classList.add("tree"),
        });
      });
      // counters
      setTimeout(() => {
        eEl.textContent = i + 1;
        qEl.textContent = qubitsAdded.size;
        const obj = { v: parseFloat(wEl.textContent) };
        const target = ((i + 1) / edges.length) * data.predicted_W;
        anime({ targets: obj, v: target, duration: 280, easing: "easeOutCubic",
                update: () => wEl.textContent = obj.v.toFixed(2) });
      }, startMs);
    });
  }

  // Panel 7 — pulse the n=20 row
  function animateMitigation() {
    animateNumbers(sections.find((s) => s.id === "s7-mitigation"));
    const row = document.getElementById("qrem20-sigma");
    if (!row) return;
    setTimeout(() => {
      anime({ targets: row, scale: [1, 1.18, 1], duration: 800, easing: "easeOutQuad", loop: 2 });
    }, 1400);
  }

  // Panel 9 — MBQC sequence frames
  function animateMBQC() {
    const data = JSON.parse(document.getElementById("data-mbqc").textContent);
    const labelEl = document.getElementById("mbqc-label");

    const links  = ["link-0","link-1","link-2","link-3","link-4","link-5"].map(id => document.getElementById(id));
    const halos  = [0,1,2,3,4,5,6].map(i => document.getElementById(`halo-${i}`));
    const qubits = [0,1,2,3,4,5,6].map(i => document.getElementById(`q-${i}`));
    const outs   = [1,2,3,4,5].map(i => document.getElementById(`out-${i}`));
    const pauli  = document.getElementById("pauli-string");
    const bell   = document.getElementById("bell-state");
    const fcap   = document.getElementById("fcap");

    // Reset
    [...halos, ...outs].forEach(el => el.style.opacity = 0);
    links.forEach(l => l.classList.remove("cz"));
    qubits.forEach(q => q.classList.remove("measured"));
    [pauli, bell, fcap].forEach(el => el.style.opacity = 0);

    const tl = anime.timeline({ easing: "easeOutQuad" });
    // Frame 1: cluster prep — H halos pulse, CZ links draw left-to-right
    tl.add({ targets: halos, opacity: [0, 0.7, 0], duration: 800, delay: anime.stagger(70) });
    tl.add({ targets: links, opacity: [0.25, 1], duration: 60, delay: anime.stagger(120),
             begin: () => labelEl.textContent = "1 · cluster prep — CZ links" }, "-=400");
    links.forEach((l, i) => tl.add({ duration: 1, begin: () => l.classList.add("cz") }, 800 + i * 120));

    // Frame 2: internal X-measurements
    tl.add({ duration: 1, begin: () => labelEl.textContent = "2 · internal X-measurements" }, "+=400");
    data.outcomes.forEach((m, i) => {
      const sign = m === 0 ? "+1" : "−1";
      tl.add({ duration: 1, begin: () => { outs[i].textContent = sign; } }, "+=0");
      tl.add({ targets: outs[i], opacity: [0, 1], translateY: [-8, 0], duration: 350 }, "-=1");
      tl.add({ targets: qubits[i + 1], duration: 350, begin: () => qubits[i + 1].classList.add("measured") }, "-=300");
    });

    // Frame 3: byproduct overlay
    tl.add({ duration: 1, begin: () => labelEl.textContent = "3 · byproduct correction" }, "+=300");
    const a = data.outcomes.filter((_, i) => i % 2 === 0).reduce((s, x) => s + x, 0);
    const b = data.outcomes.filter((_, i) => i % 2 === 1).reduce((s, x) => s + x, 0);
    pauli.textContent = `byproduct = Z^${a} X^${b}`;
    tl.add({ targets: pauli, opacity: [0, 1], translateY: [12, 0], duration: 600 }, "-=1");

    // Frame 4: |Φ+⟩ on endpoints
    tl.add({ duration: 1, begin: () => labelEl.textContent = "4 · |Φ⁺⟩ recovered on (A, B)" }, "+=400");
    tl.add({ targets: bell, opacity: [0, 1], translateY: [10, 0], duration: 500 }, "-=1");
    tl.add({ targets: fcap, opacity: [0, 1], duration: 500 }, "-=200");
  }

  // Panel 10 — variety table stagger + thanks
  function animateSummary() {
    anime({
      targets: "#s10-summary table tr",
      opacity: [0, 1], translateX: [-12, 0],
      duration: 500, delay: anime.stagger(60), easing: "easeOutQuart",
    });
    anime({
      targets: "#s10-summary li",
      opacity: [0, 1], translateY: [10, 0],
      duration: 400, delay: anime.stagger(50, { start: 200 }), easing: "easeOutQuart",
    });
  }

  // Trigger panel 1 animation immediately on load
  setTimeout(() => runAnimationFor(0), 200);
  triggered.add(0);
  dots[0]?.classList.add("active");
})();

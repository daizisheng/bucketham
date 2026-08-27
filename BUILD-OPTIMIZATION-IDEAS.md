# Build-phase optimization backlog

The build cost ≈ **peak state count × columns swept × per-state cost**.
Measured split (m=6, AGGR=2): **56% expand + 43% sort/reduce + ~0% compress**.
Scaling: peak state ~×32 per row → m6 27s, m7 ~15min, m8 ~8h / ~350GB–1TB (memory wall).

Goal to keep in mind: ultimately we want **a₈,b₈** in `(a+bn)·ρ^n` (F8A ex210), and
exact open m×n counts along the way.

Status legend: ⬜ todo · 🔬 experimenting · ✅ done/adopted · ❌ tried, rejected · 💤 parked

---

## A. Kill redundant re-expansion (algorithmic — biggest lever)

- **A1. Memoization / closure BFS** 🔬 ← **chosen; committing to the PARALLEL version**

  **Cost analysis (2026-08-27) — must be PARALLEL to win:** the current sweep is already
  16-core parallel. Effective wall (in "stable-set×" state-expansion units):
  - current parallel sweep: ~8× expansions / 16 cores ≈ **0.5×**
  - serial BFS: ~1.1× / 1 core ≈ **1.1×** → *2× SLOWER than the sweep* (loses the 16×)
  - parallel BFS: ~1.1× / 16 cores ≈ **0.07×** → ~8× faster (if the cache parallelizes)

  So A1 needs a **thread-safe sharded cache** (same hard problem as the parallel hash-agg).
  The 43% sort is NOT helped (successors still deduped). Boundary caveat: cols 0,1 have a
  different (boundary) transfer, don't cache them.

  **Design (sharded parallel BFS / closure):** partition the state space by hash(key)%P;
  each shard owned by one thread with its own hash(key→id)+worklist. Pop→expand(add_derived)
  →successors; same-shard insert+enqueue, cross-shard route via an outbox; exchange outboxes;
  repeat to global closure. Each (substep,state) expanded ONCE; edges collected per shard.
  BFS yields the bulk stable set + periodic block; the transient (cols 0..R) still uses the
  short sweep, then junction+replay as now.

  **Plan:** (1) serial BFS closure prototype → validate it reproduces the periodic block /
  OEIS on m=5/6 (correctness of "each state once = correct transfer"); (2) shard it for
  parallelism; (3) integrate with transient + junction.

  Memoization cache ≈ the periodic block (already stored under aggressive) → not much extra
  memory for m=7 (~fits 128GB); for m=8 the cache is ~the full transfer (big) — A1 mainly a
  build-SPEED lever, less an m=8-memory one.
  Build the transfer matrix by BFS from the seed: expand each *reachable* state ONCE,
  store its successor edges in a growing edge table (hash state→id + edges), follow to
  new states until closure. The sweep currently re-expands the accumulated set every
  column (~5–8×, since Sₖ₋₁ ⊆ Sₖ). BFS expands each reachable state exactly once.
  - Payoff: expand ~5–8× → **build ~2–3×**; also removes the sort/dedup (states are
    discovered distinct). Cost: holds the edge table in RAM → **+memory** (pair with B1).
  - The state graph is layered by substep mod m (periodic); BFS fills each (phase,state)
    once. The edge table IS the saved transfer matrix — replay runs on it directly.

- **A2. Direct enumeration of the state universe** 💤 (calibrated; A1 is strictly better)
  Idea was to enumerate the frontier-state universe and compute edges directly, no sweep.
  - **Calibration (2026-08-27):** states are NOT non-crossing (knight paths cross in the
    plane), so the universe = ALL partial matchings + bare/inner labels. Exact count
    `U(q)=Σ_k C(q,2k)(2k−1)!!·2^{q−2k}`. Verified on m=5 (q=11): U=538,078, reachable=143,448,
    **reachable ⊆ universe ✓**, ratio **U/R = 3.75** (only 26.7% of the universe reachable).
    Ratio falls with m: 3.75(m5) → 2.6(m6) → 1.7(m7) → ~1.4(m8).
  - **Verdict:** the universe is only O(1)× the reachable set (good), BUT A1's BFS expands
    ONLY reachable states (1×) with no over-generation, and stores only the reachable edge
    table — strictly better than enumerating+pruning the universe (1.7–3.75×). So do A1.
    (A2 stays parked as the fallback if BFS turns out awkward.)

## B. Reduce the state count / its footprint

- **B1. In-RAM state-set compression (delta keys)** ⬜
  The state set is sorted keys with high prefix redundancy → delta+varint to ~5 B/key
  (vs 15–20 B). Cuts state-set memory; offsets A1's caching cost.
  - Payoff: ~3× on the state-set memory. Low-ish risk.

- **B2. State minimization (merge count-equivalent states)** ⬜
  If two mate-patterns have identical future-count behaviour (Nerode/bisimulation),
  merge them → shrink the transfer-matrix dimension itself (the fundamental driver).
  - Payoff: uncertain (frontier states may already be near-minimal); if it bites, it's
    fundamental. Risk/cost: computing equivalence is expensive.

- **B3. More symmetry** ❌ (essentially exhausted)
  Row reflection already used (2×). Column reflection / rotation are global, not
  transfer-compatible. Nothing more local to exploit.

- **B4. Better cell-processing order** 💤
  Outside-in (OUTIN) already compacted the profile 1.34×. A different order might lower
  the mid-column peak further, but marginal.

## C. Hardware

- **C1. GPU expand + dedup** ⬜
  Expand (one thread per state) and dedup (GPU hash/radix) are massively regular-parallel.
  - Payoff: **10–100×** on the parallel parts — the single biggest raw lever, esp. for m=8.
  - Cost: large (state encoding, GPU hashing, memory management).

## D. m=8 memory wall (memory ~350GB–1TB is the real wall; time ~8h is secondary)

- **D1. External-memory sweep** ⬜
  Stream successors/states to disk, external merge-sort; RAM = buffers only.
  - Payoff: **the only way to fit m=8 in 128GB.** Disk-bound (slower) but fits.

## E. Constant factors (safe, low risk)

- **E1. Bit-pack the state key** ⬜  (~5 bits/position not 1 byte → ~35% key memory + faster memcmp/hash)
- **E2. Radix sort for the reduce** ⬜  (O(n) vs O(n log n) → ~1.3–1.5× on the 43%; fixed-length keys suit it; needs a temp buffer = +memory)
- **E3. SIMD / cache-friendly layout** 💤  (vectorize mate-splice + key compare)

## F. Reframe the goal (helps overall, not build)

- **F1. Spectral a,b (fit in the asymptotic regime)** ⬜
  a,b,ρ are the transfer matrix's dominant eigenvalue + sub-leading structure. Once the
  matrix is built, a,b need only ~tens of columns in the asymptotic regime to fit
  `(a+bn)·ρ^n` — no need to iterate exact to n=200. Shortens the REPLAY, not the build.

---

## Priority (my recommendation)

1. **A1 memoization + B1 state compression** — the practical combo to make m=7 / near-m=8
   build fast within 128GB (A1 speeds, B1 offsets its memory cost).
2. **A2 direct enumeration** — potentially removes the sweep entirely; do the cheap
   universe/reachable ratio experiment first.
3. **C1 GPU** or **D1 external memory** — the big mid-term levers (m=8).

## Already done this session (context)
weight-free build · aggressive recording · edge compression (parallel) · exact CRT
mode · parallel split-point replay decode · exact period detection (count→full compare).
Branch: `weightfree-wip`. Validated m=5(n≤50)/m=6(n≤32)/m=7(n≤12) vs OEIS b-files.

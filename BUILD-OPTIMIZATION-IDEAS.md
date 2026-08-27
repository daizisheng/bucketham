# Build-phase optimization backlog

The build cost ≈ **peak state count × columns swept × per-state cost**.
Measured split (m=6, AGGR=2): **56% expand + 43% sort/reduce + ~0% compress**.
Scaling: peak state ~×32 per row → m6 27s, m7 ~15min, m8 ~8h / ~350GB–1TB (memory wall).

Goal to keep in mind: ultimately we want **a₈,b₈** in `(a+bn)·ρ^n` (F8A ex210), and
exact open m×n counts along the way.

Status legend: ⬜ todo · 🔬 experimenting · ✅ done/adopted · ❌ tried, rejected · 💤 parked

---

## A. Kill redundant re-expansion (algorithmic — biggest lever)

- **A1. Memoization / closure BFS** ⬜
  Expand each *distinct* state ONCE (cache its successor edges); later columns hit
  the cache instead of re-running `add_derived`. The sweep currently re-expands the
  accumulated set every column (~8× total, since Sₖ₋₁ ⊆ Sₖ).
  - Payoff: expand ~8× → **build ~2–3×**. Cost: caches the edge table → **+memory**.
  - Depends/pairs with: B1 (compress the cached state/edges).

- **A2. Direct enumeration of the state universe** ⬜
  Frontier states = non-crossing partial matchings on 2m+1 positions (Catalan family).
  Enumerate them combinatorially and compute each state's transfer edges directly —
  **no growth sweep, no re-expansion, transfer matrix in one shot**. Reachability is
  then left to the SpMV (zeros for unreachable).
  - Payoff: eliminates the discovery sweep + all re-expansion. Bigger than A1.
  - Risk: the "valid non-crossing" universe may be >> the reachable set (over-generation).
  - **Cheap pre-experiment: measure |universe| / |reachable| before committing.**

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

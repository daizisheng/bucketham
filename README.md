# bucketham

A from-scratch, **parallel-friendly** transfer engine for counting **open
Hamiltonian paths (open "tours") of the m×n knight graph**, written as a
[CWEB](https://www-cs-faculty.stanford.edu/~knuth/cweb.html) literate program.

It is a bucket/sort reformulation of the broken-profile (frontier) method used
by Knuth's `DYNAHAM`. Where `DYNAHAM` stores frontier classes in a **trie**
(compact, but a serial, latency-bound, pointer-chasing update), bucketham
**emits every successor as a record, sorts by a canonical key, then reduces
equal keys**. Sorting is bandwidth-bound and embarrassingly parallel, so the
whole transfer parallelizes — the point being to scale to `m=8` and produce
many terms of the open-tour counts.

## Build

```sh
make          # ctangle + compile  ->  ./bucketham
make pdf      # cweave -> bucketham.pdf   (the readable literate document)
make check    # run the built-in self-tests
```

Needs a C compiler and (for `pdf`) a CWEB + TeX installation
(`ctangle`, `cweave`, `pdftex`).

## Use

```sh
./bucketham            # self-tests: open m×n for known m,n
./bucketham 5 8        # open 5×8  = 18061054
./bucketham w 5 12     # open 5×c for every c=1..12, in one strip sweep
./bucketham s 5 11     # stabilization diagnostic (translation-invariant states)
```

## Programs

- **`dualham.w`** — the current engine. A **dual-frontier** transfer (track the
  *unprocessed* boundary + an apex vertex; open tours = cycles through the apex),
  whose state counts match Knuth's `DYNAHAM` exactly. It sweeps once to count
  open `m×c` for every `c`, detects the periodic transfer once the frontier
  stabilizes, extracts the integer edge tables, and then reaches far columns by
  a **sparse matrix–vector product** — `build the stable transfer once, iterate
  cheaply`. Validated: open `5×c` matches known values, and the SpMV agrees with
  the direct sweep and extends correctly beyond the built range.
- `bucketham.w` — the first version (tracks *processed* cells). Correct but its
  mid-column state set is ~10× larger; kept for reference. `dual.c`,
  `periodic.c` are the C prototypes the CWEB was derived and validated from.

```sh
make dualham                     # tangle + compile the dual engine
./dualham                        # self-checks
./dualham 5 20 30                # build, then SpMV-extend to col 30 (one process)

# crash-resilient, two-phase (recommended for big m):
./dualham build 8 40 t.bin ck.bin 8 # build+extract, checkpoint every 8 columns to ck.bin,
                                     #   dump the periodic tables to t.bin
./dualham resume ck.bin 40 t.bin     # a crashed build resumes from ck.bin
./dualham run t.bin 300              # reload tables, SpMV open 8xc out to col 300
```

Weights are exact `u64` by default; pass a prime to `run` for one residue, and
`crt.py` combines several residues into exact big integers:

```sh
./dualham build 6 20 t6.bin ck6.bin   # build once (Wb=20 >= stab col 9 + margin)
python3 crt.py ./dualham t6.bin 6 30 8 # 8 primes -> exact open 6xc to col 30
```

The build (`expand` + `sort/reduce` + recording) and the SpMV are OpenMP-parallel.
A background watcher samples the **resident set size** and exits gracefully
(`exit 3`) if it exceeds 85% of RAM (override with `MEMCAP_GB=N`); every allocation
is checked too. Bounding physical memory (not the virtual address space) is what
keeps a run from ever triggering the kernel OOM killer.

**Choosing `Wb`.** `Wb` is the strip width the build sweeps to discover the
periodic transfer. It must be a few columns wider than the stabilization column
so the *recorded* columns sit in the bulk, away from the strip's right edge —
concretely `Wb >= (stabilization col) + 4`. The build reports the stabilization
column (`stable col C, period P (confirmed)`); if `Wb` is too small it prints
`no confirmed period ... : direct counts only` and simply reports the directly
counted columns (no extrapolation). For `m=6` (stabilizes ~col 9) use `Wb>=20`;
larger `m` stabilize later and need a larger `Wb`.

## Status

Validated: dual state machine (state counts match DYNAHAM); whole-table sweep;
translation-invariant key + period detection; periodic edge-table extraction;
SpMV extension. Next: OpenMP over the bucket sort and SpMV; modular/CRT (or
double) arithmetic for large `n`; the `m=8` run.

Background: this targets reliable `a,b` in the `8×n` open-tour asymptotics
`(a+b·n)·ρ^n`, ρ ≈ 526.458 — Knuth *TAOCP* Pre-Fascicle 8A, Exercise 210.

## Background

This targets reliable estimates of the leading coefficients in the `8×n`
open-tour asymptotics `(a + b·n)·ρ^n`, `ρ ≈ 526.458` — a question raised in
Knuth's *TAOCP* Pre-Fascicle 8A, Exercise 210.

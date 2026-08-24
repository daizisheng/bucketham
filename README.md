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

## Status

Validated:

- **State machine** (frontier canon/decode, edge splice, freeze accounting):
  open `m×n` correct for `m = 4,5,6,7` — built-in checks, exact agreement with
  known open 5×n values, and transpose symmetry `open(m,n)=open(n,m)` computed
  via a different vertex order.
- **Whole-table sweep**: open `m×c` for every `c` from one strip run.
- **Translation-invariant key** + stabilization detection (the prerequisite for
  extracting the reusable periodic transfer).

In progress: tightening the canonical form to the minimal frontier code;
extracting the per-substep edge tables and the parallel sparse matrix-vector
product that scales to large `n`; OpenMP; modular/CRT exact arithmetic.

## Background

This targets reliable estimates of the leading coefficients in the `8×n`
open-tour asymptotics `(a + b·n)·ρ^n`, `ρ ≈ 526.458` — a question raised in
Knuth's *TAOCP* Pre-Fascicle 8A, Exercise 210.

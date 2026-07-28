# Verification scripts

Repository: https://github.com/thawatchai2799/DR_Code_Clustering_20260728_0828

Reproduces every numerical claim in

> W. Sriphum and T. Chomsiri,
> *Symmetry-Reduced Enumeration and Canonical Forms for Grid-Based Density
> Clustering under the Hyperoctahedral Group*, Symmetry (2026).

## Requirements

Python 3.9 or later. Both scripts use nothing outside the standard library, so
there is nothing to install.

## Run

Everything at once:

```
./run_all.sh          # macOS / Linux
run_all.bat           # Windows
```

Or individually:

```
python3 verify_all.py
python3 verify_dedup_independent.py
```

Combined runtime is under a minute (about 50 s for `verify_all.py` and about 6 s
for `verify_dedup_independent.py` on a 2024 laptop).

## What is checked

### `verify_all.py`

Prints an `OK` for each checked claim:

- **T1** the hyperoctahedral group `B_d` acts on `[m]^d` by adjacency-preserving
  bijections, the wreath-product composition rule agrees with composition of the
  induced bin-maps, and the inverse formula is correct;
- **T2** grid density clustering is `B_d`-equivariant;
- **T3** the Burnside count `N(d,m,q)` equals exhaustive orbit enumeration;
- **T4** the general cycle index (product-action / gcd–lcm formula) reproduces the
  cycle count `c(g)`;
- **T5** the canonical form is orbit-invariant and separates orbits;
- **C4** the orbit-membership test: `canon(w) == canon(w')` iff `w` and `w'` share
  an orbit;
- the closed form for `N(2,m,q)` and the parity-split closed form for `N(3,m,q)`
  (Corollary 2) match Burnside for all tested `m, q`;
- the worked micro-example (`d=2, m=2, q=2`) reproduces the six orbits with the
  stated sizes and stabilisers;
- **Figure 6a** the fraction of random occupancy functions with a non-trivial
  stabiliser, which reaches 0% at `d=3, m=3, q=4` over the sampled instances;
- **Section 7** the constructed `B_d`-closed library (400 instances collapsing to
  50 canonical forms, 8.00x) and the measured canonicalisation work as a fraction
  of the naive bound (0.51–0.64).

### `verify_dedup_independent.py`

Reproduces **Figure 6b**. The 400 → 50 figure above is a construction: the
library was built from complete orbits, so it demonstrates only that the bound of
Definition 4 is attainable. This script measures the complementary quantity — on
instances drawn *independently*, how much deduplication does the canonical form
add beyond what hashing the occupancy function already achieves?

```
extra factor = (number of distinct w) / (number of distinct canon(w))
```

A value of 1.000 means canonicalisation removed nothing that a plain hash of `w`
had not already removed. Over 20,000 draws per setting:

| d | m | q | distinct w | distinct canon | extra factor |
|---|---|---|-----------:|---------------:|-------------:|
| 2 | 3 | 2 |        512 |            102 |       5.020x |
| 2 | 4 | 2 |     17,238 |          7,671 |       2.247x |
| 2 | 5 | 2 |     19,994 |         19,945 |       1.002x |
| 2 | 6 | 2 |     20,000 |         20,000 |       1.000x |
| 3 | 2 | 2 |        256 |             22 |      11.636x |
| 3 | 3 | 2 |     19,998 |         19,927 |       1.004x |
| 2 | 4 | 4 |     20,000 |         19,999 |       1.000x |
| 3 | 3 | 4 |     20,000 |         20,000 |       1.000x |

The first two rows are an artefact of a workload larger than the instance space:
at `d=2, m=3, q=2` the 20,000 draws saturate all 512 patterns, which collapse to
the 102 orbits of Table 5. Once the space is large relative to the workload —
the case at every resolution of practical interest — the extra factor is 1.00.
The script self-checks each value against the figure quoted in the paper.

## Files

| File | Purpose |
|---|---|
| `verify_all.py` | all theorem, enumeration and canonical-form checks |
| `verify_dedup_independent.py` | the independent-instance deduplication measurement (Figure 6b) |
| `run_all.sh` / `run_all.bat` | run both scripts in order |
| `requirements.txt` | records that there are no dependencies |
| `LICENSE` | MIT |
| `README.md` / `README.html` | this document |

## Licence

MIT. See `LICENSE`.

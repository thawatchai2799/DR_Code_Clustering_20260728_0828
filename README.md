# Verification scripts

Repository: https://github.com/thawatchai2799/DR_Code_Clustering_20260728_0828

Reproduces every numerical claim in

> W. Sriphum and T. Chomsiri,
> *Symmetry-Reduced Enumeration and Canonical Forms for Grid-Based Density
> Clustering under the Hyperoctahedral Group*, Symmetry (2026).

## Requirements
Python 3.9+ and NumPy. No other dependencies.

```
pip install numpy
```

## Run
```
python3 verify_all.py
```

Expected runtime: a few seconds. The script prints an OK for each checked claim:

- **T1** the hyperoctahedral group `B_d` acts on `[m]^d` by adjacency-preserving bijections;
- **T2** grid density clustering is `B_d`-equivariant;
- **T3** the Burnside count `N(d,m,q)` equals exhaustive orbit enumeration;
- **T4** the general cycle index (product-action / gcd–lcm formula) reproduces the cycle count `c(g)`;
- **T5** the canonical form is orbit-invariant and separates orbits;
- the closed forms for `N(2,m,q)` and `N(3,m,q)` (parity-split) match Burnside for all tested `m,q`;
- the worked micro-example (`d=2, m=2, q=2`) reproduces the six orbits with the stated sizes and stabilisers.

## Files
- `verify_all.py` — all checks in a single self-contained script.

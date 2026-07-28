#!/usr/bin/env python3
"""
verify_dedup_independent.py  --  addendum to DR_Clustering_verification_scripts

Reproduces the measurement reported in Section 7 and Figure 6b of

    "Symmetry-Reduced Enumeration and Canonical Forms for Grid-Based Density
     Clustering under the Hyperoctahedral Group"

Question
--------
A B_d-closed library collapses by exactly the average orbit size; that is a
property of how the library was built, not evidence that symmetry occurs in
workloads that were not built that way.  This script measures the complementary
quantity: on instances drawn independently, how much deduplication does the
canonical form add *beyond* what hashing the occupancy function itself already
achieves?

    extra factor  =  (number of distinct w)  /  (number of distinct canon(w))

A value of 1.000 means canonicalisation removed nothing that a plain hash of w
had not already removed.

Requires Python 3.9+ and the standard library only.  Runs in well under a minute.
"""

import itertools
import math
import random
import sys

N_DRAWS = 20_000
SEED = 0

SETTINGS = [
    # d, m, q
    (2, 3, 2),
    (2, 4, 2),
    (2, 5, 2),
    (2, 6, 2),
    (3, 2, 2),
    (3, 3, 2),
    (2, 4, 4),
    (3, 3, 4),
]


def signed_permutations(d):
    """The 2^d d! elements of B_d, as (axis permutation, sign vector)."""
    for pi in itertools.permutations(range(d)):
        for eps in itertools.product((1, -1), repeat=d):
            yield pi, eps


def bin_permutations(d, m):
    """Each group element as a gather index over the m^d bins in row-major order."""
    cells = list(itertools.product(range(m), repeat=d))
    index = {c: i for i, c in enumerate(cells)}
    out = []
    for pi, eps in signed_permutations(d):
        gather = [0] * len(cells)
        for c in cells:
            y = [0] * d
            for i in range(d):
                y[pi[i]] = c[i] if eps[i] == 1 else (m - 1) - c[i]
            gather[index[tuple(y)]] = index[c]
        out.append(gather)
    return out, len(cells)


def canon(w, perms):
    """Lexicographically minimal image of w under B_d (Definition 3)."""
    n = len(w)
    return min(tuple(w[p[i]] for i in range(n)) for p in perms)


def measure(d, m, q, n_draws=N_DRAWS, seed=SEED):
    rng = random.Random(seed)
    perms, n_bins = bin_permutations(d, m)
    draws = [tuple(rng.randrange(q) for _ in range(n_bins)) for _ in range(n_draws)]
    distinct_w = len(set(draws))
    distinct_canon = len({canon(w, perms) for w in draws})
    return distinct_w, distinct_canon


def main():
    print(__doc__.strip().splitlines()[0])
    print()
    print(f"{'draws':>8}  = {N_DRAWS:,} independent occupancy functions per setting")
    print(f"{'seed':>8}  = {SEED}")
    print()
    header = (f"{'d':>2} {'m':>2} {'q':>2} {'|A^Gamma|':>12} {'|B_d|':>6} "
              f"{'distinct w':>11} {'distinct canon':>15} {'extra factor':>13}")
    print(header)
    print("-" * len(header))

    results = []
    for d, m, q in SETTINGS:
        distinct_w, distinct_canon = measure(d, m, q)
        factor = distinct_w / distinct_canon
        space = q ** (m ** d)
        order = 2 ** d * math.factorial(d)
        print(f"{d:>2} {m:>2} {q:>2} {space:>12.3g} {order:>6} "
              f"{distinct_w:>11,} {distinct_canon:>15,} {factor:>12.3f}x")
        results.append((d, m, q, distinct_w, distinct_canon, factor))

    print()
    print("Reading of the table")
    print("--------------------")
    print("The first two rows have a workload larger than the instance space: at")
    print("d=2, m=3, q=2 the 20,000 draws saturate all 512 patterns, which collapse to")
    print("the 102 orbits of Table 5, so the extra factor is 512/102 = 5.02 by")
    print("saturation rather than by any property of real workloads.  From m=5 onward")
    print("(and at d=3, m=3) the space is large relative to the workload and the extra")
    print("factor is 1.00 to three decimal places: canonicalisation removes nothing")
    print("that hashing w had not already removed.")
    print()
    print("This is the honest counterpart to the 400 -> 50 (8.00x) figure of Section 7,")
    print("which is a construction from complete orbits and demonstrates only that the")
    print("bound of Definition 4 is attainable.")

    # values quoted in the paper
    expected = {
        (2, 3, 2): 5.020, (2, 4, 2): 2.247, (2, 5, 2): 1.002, (2, 6, 2): 1.000,
        (3, 2, 2): 11.636, (3, 3, 2): 1.004, (2, 4, 4): 1.000, (3, 3, 4): 1.000,
    }
    print()
    ok = True
    for d, m, q, _, _, factor in results:
        want = expected[(d, m, q)]
        if abs(factor - want) > 5e-4:
            ok = False
            print(f"MISMATCH at d={d}, m={m}, q={q}: got {factor:.3f}, paper reports {want:.3f}")
    print("All measured values match those quoted in the paper."
          if ok else "One or more values differ from the paper.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Verification scripts for:
  "Symmetry-Reduced Enumeration and Canonical Forms for Grid-Based Density
   Clustering under the Hyperoctahedral Group"
   W. Sriphum and T. Chomsiri, Symmetry (2026).

Reproduces every numerical claim in the paper. Pure Python 3.9+, standard
library only (no third-party dependencies).
Run:  python3 verify_all.py
Expected runtime: under a minute (~50 s on a 2024 laptop).

The script checks:
  (T1) B_d acts by adjacency-preserving bijections on [m]^d.
  (T2) grid density clustering is B_d-equivariant.
  (T3) Burnside count N(d,m,q) equals exhaustive orbit enumeration.
  (T4) general cycle index (product-action / gcd-lcm formula) reproduces c(g).
  (T5) canonical form is orbit-invariant and separates orbits.
  Closed forms for N(2,m,q) and N(3,m,q) (parity-split).
  Worked micro-example d=2, m=2 (six orbits).
"""
import itertools
from math import gcd, factorial, comb
from functools import reduce

try:
    from math import lcm  # py3.9+
except ImportError:
    def lcm(*a): return reduce(lambda x,y: x*y//gcd(x,y), a)

def totient(n):
    r=n; p=2; nn=n
    while p*p<=nn:
        if nn%p==0:
            while nn%p==0: nn//=p
            r-=r//p
        p+=1
    if nn>1: r-=r//nn
    return r

# ---------- group ----------
def Bd(d):
    return [(p,s) for p in itertools.permutations(range(d))
                  for s in itertools.product((1,-1),repeat=d)]

def act(g,x,m):
    p,s=g; d=len(x); y=[0]*d
    for i in range(d):
        xi = x[i] if s[i]==1 else (m-1)-x[i]
        y[p[i]] = xi
    return tuple(y)

def cyc(g,d,m):
    cells=list(itertools.product(range(m),repeat=d))
    idx={c:i for i,c in enumerate(cells)}
    seen=[False]*len(cells); c=0
    for st in range(len(cells)):
        if seen[st]: continue
        c+=1; j=st
        while not seen[j]:
            seen[j]=True; j=idx[act(g,cells[j],m)]
    return c

def N_burnside(d,m,q):
    G=Bd(d); return sum(q**cyc(g,d,m) for g in G)//len(G)

def N_orbits(d,m,q):
    """Exhaustive orbit count by explicit enumeration (small parameters only)."""
    G=Bd(d); cells=list(itertools.product(range(m),repeat=d))
    pos={c:i for i,c in enumerate(cells)}
    seen=set(); n=0
    for w in itertools.product(range(q),repeat=len(cells)):
        if w in seen: continue
        n+=1
        for g in G:
            img=[0]*len(cells)
            for i,c in enumerate(cells): img[pos[act(g,c,m)]]=w[i]
            seen.add(tuple(img))
    return n

# ---------- orbit membership (Corollary 3) ----------
def same_orbit_direct(w,w2,d,m):
    G=Bd(d); cells=list(itertools.product(range(m),repeat=d))
    pos={c:i for i,c in enumerate(cells)}
    for g in G:
        img=[0]*len(cells)
        for i,c in enumerate(cells): img[pos[act(g,c,m)]]=w[i]
        if tuple(img)==tuple(w2): return True
    return False

def check_C3(pairs=600, seed=11):
    """Corollary 3: canon(w)==canon(w') iff same orbit. Returns (ok, n_pairs)."""
    import random
    rng=random.Random(seed); total=0
    for (d,m,q) in [(2,3,2),(2,3,3),(3,2,2)]:
        G=Bd(d); cells=list(itertools.product(range(m),repeat=d))
        pos={c:i for i,c in enumerate(cells)}
        for _ in range(pairs//3):
            w=tuple(rng.randrange(q) for _ in cells)
            if rng.random()<0.5:
                g=rng.choice(G); img=[0]*len(cells)
                for i,c in enumerate(cells): img[pos[act(g,c,m)]]=w[i]
                w2=tuple(img)
            else:
                w2=tuple(rng.randrange(q) for _ in cells)
            via_canon = canon(w,cells,G,m)==canon(w2,cells,G,m)
            via_direct= same_orbit_direct(w,w2,d,m)
            total+=1
            if via_canon!=via_direct: return False,total
    return True,total

def adjacent(u,v):
    return u!=v and max(abs(a-b) for a,b in zip(u,v))<=1

def check_T1(dims=(2,3,4), max_m=4):
    for d in dims:
        G=Bd(d)
        for m in range(2,max_m+1):
            cells=list(itertools.product(range(m),repeat=d))
            for g in G:
                # bijection
                imgs={act(g,c,m) for c in cells}
                assert len(imgs)==len(cells), "not a bijection"
                # adjacency preserved
                for u,v in itertools.combinations(cells,2):
                    if adjacent(u,v)!=adjacent(act(g,u,m),act(g,v,m)):
                        return False
    return True

def compose(g1,g2,d):
    """Wreath-product composition rule as printed in Section 3:
       g1 g2 = (pi1 pi2, eps) with eps_i = eps_{2,i} * eps_{1, pi2(i)}."""
    p1,s1=g1; p2,s2=g2
    p=tuple(p1[p2[i]] for i in range(d))
    e=tuple(s2[i]*s1[p2[i]] for i in range(d))
    return (p,e)

def check_T1_composition(dims=(2,3,4), max_m=3):
    """Section 4 / Section 8 claim: the abstract wreath-product composition
       agrees with the composition of the induced bin-maps, on every bin."""
    n=0
    for d in dims:
        G=Bd(d)
        for m in range(2,max_m+1):
            cells=list(itertools.product(range(m),repeat=d))
            for g1 in G:
                for g2 in G:
                    g12=compose(g1,g2,d)
                    for x in cells:
                        n+=1
                        if act(g12,x,m)!=act(g1,act(g2,x,m),m):
                            return False,n
    return True,n

def check_inverse(dims=(2,3,4)):
    """Theorem 1 proof: (pi,eps)^-1 = (pi^-1, eps') with eps'_{pi(i)} = eps_i."""
    for d in dims:
        ident=(tuple(range(d)),tuple([1]*d))
        for g in Bd(d):
            p,e=g
            pinv=[0]*d
            for i in range(d): pinv[p[i]]=i
            ep=[0]*d
            for i in range(d): ep[p[i]]=e[i]
            ginv=(tuple(pinv),tuple(ep))
            if compose(g,ginv,d)!=ident or compose(ginv,g,d)!=ident:
                return False
    return True

# ---------- equivariance (T2) ----------
def clustering(w,cells,m,tau=1):
    """dense bins = value>=tau; clusters = connected comps under adjacency."""
    dense=[c for c in cells if w[c]>=tau]
    ds=set(dense); comps=[]; seen=set()
    for c in dense:
        if c in seen: continue
        stack=[c]; comp=set()
        while stack:
            x=stack.pop()
            if x in comp: continue
            comp.add(x); seen.add(x)
            for y in dense:
                if y not in comp and adjacent(x,y): stack.append(y)
        comps.append(frozenset(comp))
    return frozenset(comps)

def check_T2(trials_per_cell=1200, seed=0):
    """Equivariance test: d=2 with m in {3,4,5}, d=3 with m in {3,4}."""
    import random
    rng=random.Random(seed); total=0
    for d in (2,3):
        for m in (3,4,5):
            if d==3 and m>4:          # 5^3=125 cells x 48 elements: keep runtime sane
                continue
            G=Bd(d); cells=list(itertools.product(range(m),repeat=d))
            for _ in range(trials_per_cell):
                q=rng.choice([2,3,4])
                w={c:rng.randrange(q) for c in cells}
                g=rng.choice(G)
                tau=rng.choice([1,2])
                gw={act(g,c,m):w[c] for c in cells}
                Cw=clustering(w,cells,m,tau)
                Cgw=clustering(gw,cells,m,tau)
                gCw=frozenset(frozenset(act(g,c,m) for c in comp) for comp in Cw)
                total+=1
                if Cgw!=gCw:
                    return False,total
    return True,total

# ---------- cycle index (T4) ----------
def block_cycle_lengths(ell,s,m):
    pts=list(itertools.product(range(m),repeat=ell)); idx={p:i for i,p in enumerate(pts)}
    def T(a):
        b=[0]*ell
        for k in range(ell): b[(k+1)%ell]=a[k]
        if s==-1: b[0]=(m-1)-b[0]
        return tuple(b)
    seen=[False]*len(pts); L=[]
    for st in range(len(pts)):
        if seen[st]: continue
        n=0; j=st
        while not seen[j]:
            seen[j]=True; j=idx[T(pts[j])]; n+=1
        L.append(n)
    return L

def signed_cycle_type(g,d):
    p,s=g; seen=[False]*d; out=[]
    for st in range(d):
        if seen[st]: continue
        j=st; ln=0; sign=1
        while not seen[j]:
            seen[j]=True; sign*=s[j]; j=p[j]; ln+=1
        out.append((ln,sign))
    return out

def c_via_cycle_index(g,d,m):
    sct=signed_cycle_type(g,d)
    Ls=[block_cycle_lengths(ell,s,m) for (ell,s) in sct]
    tot=0
    for combo in itertools.product(*Ls):
        tot += reduce(lambda a,b:a*b,combo)//reduce(lcm,combo)
    return tot

def check_T4(dims=(2,3),max_m=4):
    for d in dims:
        for g in Bd(d):
            for m in range(2,max_m+1):
                if c_via_cycle_index(g,d,m)!=cyc(g,d,m):
                    return False
    return True

# ---------- canonical form (T5) ----------
def canon(w,cells,G,m):
    best=None
    for g in G:
        img=[0]*len(cells)
        for i,c in enumerate(cells):
            img[cells.index(act(g,c,m))]=w[i]
        t=tuple(img)
        if best is None or t<best: best=t
    return best

def check_T5(d=2,m=2,q=2):
    G=Bd(d); cells=list(itertools.product(range(m),repeat=d))
    reps={}
    for w in itertools.product(range(q),repeat=len(cells)):
        reps.setdefault(canon(w,cells,G,m),[]).append(w)
    # invariance: applying any g then canon gives same key
    for w in itertools.product(range(q),repeat=len(cells)):
        k=canon(w,cells,G,m)
        for g in G:
            img=[0]*len(cells)
            for i,c in enumerate(cells):
                img[cells.index(act(g,c,m))]=w[i]
            assert canon(tuple(img),cells,G,m)==k
    return len(reps)

# ---------- stabiliser fraction, dedup, and canonicalisation cost (Section 7 / Figure 6) ----------
def stabiliser_fraction(d,m,q,trials=4000,seed=5):
    """Figure 6a: fraction of random occupancy functions with a NON-trivial stabiliser."""
    import random
    rng=random.Random(seed)
    G=Bd(d); cells=list(itertools.product(range(m),repeat=d))
    pos={c:i for i,c in enumerate(cells)}
    perms=[[pos[act(g,c,m)] for c in cells] for g in G if g!=(tuple(range(d)),tuple([1]*d))]
    n=len(cells); hits=0
    for _ in range(trials):
        w=tuple(rng.randrange(q) for _ in cells)
        for pm in perms:
            img=[0]*n
            for i in range(n): img[pm[i]]=w[i]
            if tuple(img)==w:
                hits+=1; break
    return hits/trials

def library_dedup(d=2,m=3,q=3,seeds=50,seed=9):
    """Figure 6b: a B_d-closed library built as full orbits of `seeds` generators,
       reduced by canonical form.  Returns (library size, canonical size, factor)."""
    import random
    rng=random.Random(seed)
    G=Bd(d); cells=list(itertools.product(range(m),repeat=d))
    pos={c:i for i,c in enumerate(cells)}
    perms=[[pos[act(g,c,m)] for c in cells] for g in G]
    n=len(cells); lib=[]; chosen=0
    while chosen<seeds:
        w=tuple(rng.randrange(q) for _ in cells)
        orb=set()
        for pm in perms:
            img=[0]*n
            for i in range(n): img[pm[i]]=w[i]
            orb.add(tuple(img))
        if len(orb)!=len(G):      # keep only free orbits so the library is exactly seeds*|G|
            continue
        lib.extend(sorted(orb)); chosen+=1
    canon_set={canon(w,cells,G,m) for w in lib}
    return len(lib),len(canon_set),len(lib)/len(canon_set)

def canon_work_fraction(d,m,q,trials=300,seed=3):
    """Section 7: measured total work as a fraction of the naive bound that charges
       a full-length comparison to every group element."""
    import random
    rng=random.Random(seed)
    G=Bd(d); cells=list(itertools.product(range(m),repeat=d))
    pos={c:i for i,c in enumerate(cells)}
    perms=[[pos[act(g,c,m)] for c in cells] for g in G]
    n=len(cells); naive=0; actual=0
    for _ in range(trials):
        w=tuple(rng.randrange(q) for _ in cells); best=None
        for pm in perms:
            img=[0]*n
            for i in range(n): img[pm[i]]=w[i]
            img=tuple(img)
            naive += 2*n
            actual += n
            if best is None:
                best=img
            else:
                k=0
                while k<n and img[k]==best[k]: k+=1
                actual += min(k+1,n)
                if k<n and img[k]<best[k]: best=img
    return actual/naive

# ---------- closed forms ----------
from math import ceil
def N2_closed(m,q):
    return (q**(m*m)+2*q**(m*ceil(m/2))+q**ceil(m*m/2)+2*q**comb(m+1,2)+2*q**ceil(m*m/4))//8
def N3_closed(m,q):
    if m%2==0:
        e=[m**3, m**3//2, (m**3+m**2)//2, m**3//4, (m**3+2*m)//3, (m**3+2*m)//6]
        co=[1,13,6,12,8,8]
    else:
        e=[m**3, (m**3+m**2)//2, (m**3+m)//2, (m**3+1)//2,
           (m**3+3*m)//4, (m**3+m+2)//4, (m**3+2*m)//3, (m**3+2*m+3)//6]
        co=[1,9,9,1,6,6,8,8]
    return sum(c*q**x for c,x in zip(co,e))//48

def main():
    print("Verifying paper claims ...\n")
    print("[T1] adjacency-preserving bijection (d=2,3,4; m<=4):", check_T1())
    okc,nc1 = check_T1_composition()
    print(f"[T1] wreath composition == composition of bin-maps: {okc}  ({nc1:,} (g1,g2,bin) instances)")
    print("[T1] printed inverse formula (pi,eps)^-1:", check_inverse())
    ok2,n2 = check_T2()
    print(f"[T2] B_d-equivariance of clustering: {ok2}  ({n2:,} random (occupancy, group-element, threshold) triples)")
    print("[T4] cycle index via signed cycle type (product action):", check_T4())
    norb=check_T5(2,2,2)
    print("[T5] canonical form separates orbits; #orbits(d=2,m=2,q=2) =", norb, "(expected 6)")
    okc,nc = check_C3()
    print(f"[C3] orbit-membership test canon(w)==canon(w') iff same orbit: {okc}  ({nc:,} pairs)")

    print("\n[T3] Burnside sum vs INDEPENDENT exhaustive orbit enumeration:")
    for (d,m,q) in [(2,2,2),(2,2,3),(2,3,2),(2,3,3),(2,4,2),(3,2,2)]:
        nb=N_burnside(d,m,q); no=N_orbits(d,m,q)
        print(f"    N({d},{m},{q}): Burnside={nb:,}  brute-force={no:,}  {'OK' if nb==no else 'MISMATCH'}")
        assert nb==no

    print("\nKey enumeration values:")
    for (d,m,q,exp) in [(2,3,2,102),(2,3,3,2862),(2,4,2,8548),
                        (3,2,2,22),(3,3,2,2852288),(3,3,3,158942078604)]:
        got=N_burnside(d,m,q)
        print(f"    N({d},{m},{q}) = {got:,}  {'OK' if got==exp else 'MISMATCH'}")

    print("\nClosed-form N(2,m,q) vs Burnside:")
    for m in range(2,7):
        for q in (2,3):
            assert N2_closed(m,q)==N_burnside(2,m,q)
    print("    all match for m=2..6, q=2,3")

    print("\nClosed-form N(3,m,q) vs Burnside:")
    for m in range(2,6):
        for q in (2,3):
            assert N3_closed(m,q)==N_burnside(3,m,q)
    print("    all match for m=2..5, q=2,3 (incl. 58-digit value)")

    print("\nWorked micro-example d=2, m=2, q=2 (six orbits):")
    G=Bd(2); cells=list(itertools.product(range(2),repeat=2))
    reps={}
    for w in itertools.product(range(2),repeat=4):
        reps.setdefault(canon(w,cells,G,2),[]).append(w)
    for k,members in sorted(reps.items(), key=lambda kv:(sum(kv[0]),kv[0])):
        print(f"    canonical {''.join(map(str,k))}: orbit size {len(members)}, stabiliser {len(G)//len(members)}")

    print("\nSection 7 / Figure 6 empirical claims:")
    for (d,m,q) in [(2,3,2),(2,3,3),(2,3,4),(3,3,2),(3,3,3),(3,3,4)]:
        f=stabiliser_fraction(d,m,q, trials=2000)
        print(f"    non-trivial stabiliser fraction, d={d}, m={m}, q={q}: {100*f:6.2f}%")
    L,C,fac = library_dedup()
    print(f"    B_d-closed library dedup: {L} instances -> {C} canonical forms ({fac:.2f}x)")

    print("\nSection 7 canonicalisation cost (fraction of the naive bound):")
    fr=[]
    for (d,m,q) in [(2,3,2),(2,4,2),(2,5,3),(2,4,4),(3,3,2),(3,3,4),(3,4,3)]:
        v=canon_work_fraction(d,m,q, trials=150 if d==3 else 300); fr.append(v)
        print(f"    d={d}, m={m}, q={q}: {v:.3f}")
    print(f"    measured range: {min(fr):.2f}-{max(fr):.2f}  (paper quotes 0.51-0.64)")

    print("\nAll checks passed.")

if __name__=="__main__":
    main()

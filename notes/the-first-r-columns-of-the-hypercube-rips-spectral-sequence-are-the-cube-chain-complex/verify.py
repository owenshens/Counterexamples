#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- checks the computational claims of the accompanying paper

    "The First r+1 Columns of the Cubic-Dimension Spectral Sequence of a
     Hypercube Vietoris-Rips Complex"

Python 3.9+, STANDARD LIBRARY ONLY (itertools, sys, time, math).  No numpy, no
sympy, no networkx, no external data file, no network.  All arithmetic is exact
integer arithmetic; no decision anywhere is taken on a floating-point value.

WHAT IT READS.  Every input is transcribed below from the paper, character for
character, in the block marked FROM THE PAPER: the dimension vector (3.2) of the
63,775-generator complex, the inclusion-exclusion of Remark 3.2, the sixteen
signed terms of the cycle z printed after (3.3), and the dimension vectors that
Galetto-Montano-Wellner publish for their own scale-2 and scale-3 computations
and that the paper reproduces in Section 7.  Nothing else is fed in: the
complexes themselves are rebuilt from the definitions of Section 1.

HOW RANKS ARE TAKEN, AND WHY THAT IS SOUND OVER A FIELD OF CHARACTERISTIC ZERO.
Boundary-matrix ranks are computed by exact Markowitz-pivoted sparse Gaussian
elimination over a prime field GF(P).  A rank mod P never exceeds the rational
rank, so the mod-P Betti numbers are UPPER BOUNDS for the characteristic-zero
ones; the alternating sum sum_k (-1)^k dim C_k is field-independent and equals
sum_k (-1)^k betti_k over every field.  Hence whenever the mod-P Betti vector is
concentrated in one degree with value 1, the characteristic-zero Betti vector is
forced to be exactly the same vector.  That implication is CHECKED as its own
check (`char0-pin-*`), and it is the only bridge used between the machine's
field and the paper's.

Exit status is 0 if and only if every check passed.
"""

import itertools
import math
import sys
import time

# ===========================================================================
# FROM THE PAPER -- transcribed inputs.  Nothing here is computed.
# ===========================================================================

# Equation (3.2): dim C^{4,4}_{4+.,4} by simplicial dimension 1..15.
PAPER_DIMS_44 = [8, 208, 1284, 3920, 7784, 11376, 12862, 11440,
                 8008, 4368, 1820, 560, 120, 16, 1]
PAPER_TOTAL_44 = 63775           # Remark 3.2: 65535 - 1760
PAPER_EULER_44 = 1               # Remark 3.2: sum_k (-1)^k dim = +1

# Remark 3.2, the inclusion-exclusion for |U| at p = 4:
#   1760 = 8*255 - 24*15 + 32*3 - 16*1 = 2040 - 360 + 96 - 16
PAPER_IE_TERMS = [(8, 255), (24, 15), (32, 3), (16, 1)]
PAPER_IE_PARTIALS = [2040, 360, 96, 16]
PAPER_U4_FACES = 1760
PAPER_DELTA4_FACES = 65535       # 2^(2^4) - 1

# The sixteen signed terms of the cycle z at p = 4, exactly as the table in
# Section 3 prints them: (coefficient, eps_1 eps_2 eps_3 eps_4).  The facet is
# {(1,eps_1),(2,eps_2),(3,eps_3),(4,eps_4)}.
PAPER_Z_P4 = [
    (+1, "0000"), (-1, "0001"), (-1, "0010"), (+1, "0011"),
    (-1, "0100"), (+1, "0101"), (+1, "0110"), (-1, "0111"),
    (-1, "1000"), (+1, "1001"), (+1, "1010"), (-1, "1011"),
    (+1, "1100"), (-1, "1101"), (-1, "1110"), (+1, "1111"),
]

# Section 7, the published dimension vectors used as controls, quoted there
# from Galetto-Montano-Wellner.  Keys are (p, r); values are
# (first simplicial dimension, [dims]).
PAPER_PUBLISHED_DIMS = {
    (2, 2): (1, [2, 4, 1]),
    (3, 3): (1, [4, 32, 64, 56, 28, 8, 1]),
    (4, 3): (2, [96, 584, 1344, 1568, 960, 248]),
}
# Section 6(3): the two published off-range facts.
PAPER_OFFRANGE_43_DEGREES = [4, 7]     # H_4 = {0;1^4} and H_7 = {4;0}
PAPER_OFFRANGE_32_DIM_H3 = 2           # H_3 = {1^3;0} + {0;1^3}, so 2-dimensional

# Section 1 / Proposition 2.3: dim F_p = C(n,p) 2^{n-p}.  Section 5 / D1: at
# n = 6, dim F_4 = 60.
PAPER_N = 6
PAPER_DIM_F = {0: 64, 1: 192, 2: 240, 3: 160, 4: 60}
PAPER_INDEX_60 = 60

# The two primes.  46337 is prime and coprime to |H_4| = 384 = 2^7 * 3.
P1 = 46337
P2 = 1000003

# ===========================================================================
# harness
# ===========================================================================

T0 = time.time()
_STATE = {"pass": 0, "fail": 0}
_SCOPE = []


def ok(name, detail=""):
    _STATE["pass"] += 1
    print("PASS %s%s" % (name, (" " + detail) if detail else ""))
    sys.stdout.flush()


def bad(name, detail=""):
    _STATE["fail"] += 1
    print("FAIL %s%s" % (name, (" " + detail) if detail else ""))
    sys.stdout.flush()


def check(name, got, want, extra=""):
    if got == want:
        ok(name, ("[%s]" % (extra if extra else _short(got))))
    else:
        bad(name, "[got %s want %s %s]" % (_short(got), _short(want), extra))


def _short(v):
    s = repr(v)
    return s if len(s) <= 110 else s[:107] + "..."


def scope(line):
    _SCOPE.append(line)


def note(line):
    print("-- %s" % line)
    sys.stdout.flush()


# ===========================================================================
# exact linear algebra over GF(P)
# ===========================================================================

def rank_mod(rows, P):
    """Exact rank over GF(P) of a sparse matrix given as a list of
    {col: value} dicts.  Markowitz pivoting: the pivot minimises
    (|row|-1)*(|col|-1), searched over the 40 sparsest rows.  Pure integers."""
    Rw = {i: dict(d) for i, d in enumerate(rows) if d}
    Cl = {}
    for i, d in Rw.items():
        for c in d:
            Cl.setdefault(c, set()).add(i)
    Cl = {c: s for c, s in Cl.items() if s}
    rank = 0
    while Rw:
        best = None
        bestcost = None
        order = sorted(Rw, key=lambda i: len(Rw[i]))
        for i in order[:40]:
            for c in Rw[i]:
                cost = (len(Rw[i]) - 1) * (len(Cl[c]) - 1)
                if bestcost is None or cost < bestcost:
                    bestcost, best = cost, (i, c)
                    if cost == 0:
                        break
            if bestcost == 0:
                break
        pr, pc = best
        inv = pow(Rw[pr][pc], P - 2, P)
        prow = {c: (v * inv) % P for c, v in Rw[pr].items() if c != pc}
        for i in list(Cl[pc]):
            if i == pr:
                continue
            f = Rw[i].pop(pc, 0)
            Cl[pc].discard(i)
            if f:
                ri = Rw[i]
                for c, v in prow.items():
                    nv = (ri.get(c, 0) - f * v) % P
                    if nv:
                        if c not in ri:
                            Cl.setdefault(c, set()).add(i)
                        ri[c] = nv
                    elif c in ri:
                        del ri[c]
                        Cl[c].discard(i)
            if not Rw[i]:
                del Rw[i]
        for c in Rw[pr]:
            Cl[c].discard(pr)
        del Rw[pr]
        Cl.pop(pc, None)
        Cl = {c: s for c, s in Cl.items() if s}
        rank += 1
    return rank


def boundary_rows(bydeg, k, P):
    """The simplicial boundary d_k as {row: {col: coeff}}, rows indexed by the
    (k-1)-simplices, columns by the k-simplices.  Faces are sorted tuples."""
    cols = bydeg.get(k, [])
    tgt = bydeg.get(k - 1, [])
    idx = {s: i for i, s in enumerate(tgt)}
    rows = [dict() for _ in range(len(tgt))]
    for j, s in enumerate(cols):
        for a in range(len(s)):
            i = idx.get(s[:a] + s[a + 1:])
            if i is not None:
                rows[i][j] = ((-1) ** a) % P
    return rows


def ranks_and_betti(bydeg, P):
    """-> (dims, ranks, betti) of the chain complex `bydeg` over GF(P)."""
    dims = {k: len(v) for k, v in bydeg.items() if v}
    ranks = {}
    for k in sorted(dims):
        ranks[k] = rank_mod(boundary_rows(bydeg, k, P), P)
    betti = {k: dims[k] - ranks.get(k, 0) - ranks.get(k + 1, 0)
             for k in sorted(dims)}
    return dims, ranks, betti


def sgnk(k):
    """(-1)**k as an exact int.  Python's (-1)**-1 is a FLOAT, and this program
    takes no decision on a float, so the sign is written out."""
    return 1 if k % 2 == 0 else -1


def euler(dims):
    return sum(sgnk(k) * d for k, d in dims.items())


def char0_pin(name, dims, betti_modp, expect_degree):
    """The one bridge from GF(P) to characteristic zero, checked explicitly.

    rank_P <= rank_Q, hence betti_Q[k] <= betti_P[k] for every k.  The
    alternating sum of the dims is field-independent.  So if betti_P is 1 in a
    single degree m and 0 elsewhere, then betti_Q[k] = 0 for k != m and
    betti_Q[m] <= 1, and euler = (-1)^m pins betti_Q[m] = 1.
    """
    nz = sorted(k for k, v in betti_modp.items() if v)
    e = euler(dims)
    good = (nz == [expect_degree]
            and betti_modp[expect_degree] == 1
            and e == sgnk(expect_degree))
    check("char0-pin-" + name, good, True,
          "nonzero_mod_P=%s euler=%d (-1)^%d=%d => betti_Q = 1 in degree %d only"
          % (nz, e, expect_degree, sgnk(expect_degree), expect_degree))


# ===========================================================================
# the objects of Section 1 and Section 3, rebuilt from their definitions
# ===========================================================================

def popcount(x):
    return bin(x).count("1")


def cdim(p, verts):
    """Section 1: the number of coordinates in which the vertex set is not
    constant."""
    orv = 0
    andv = (1 << p) - 1
    for v in verts:
        orv |= v
        andv &= v
    return popcount(orv & ~andv)


def diam(verts):
    d = 0
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            d = max(d, popcount(verts[i] ^ verts[j]))
    return d


_FACE_CACHE = {}


def all_faces(p, predicate):
    """All non-empty S subset {0,1}^p (as sorted tuples) with predicate(S),
    graded by simplicial dimension |S| - 1."""
    N = 1 << p
    bydeg = {}
    for mask in range(1, 1 << N):
        verts = tuple(v for v in range(N) if (mask >> v) & 1)
        if predicate(verts):
            bydeg.setdefault(len(verts) - 1, []).append(verts)
    return bydeg


def _cached(key, build):
    if key not in _FACE_CACHE:
        _FACE_CACHE[key] = build()
    return _FACE_CACHE[key]


def rips_faces(p, r):
    """X^{p,r}: the non-empty subsets of diameter <= r."""
    return _cached(("rips", p, r),
                   lambda: all_faces(p, lambda S: diam(S) <= r))


def local_complex(p, r):
    """C^{p,r}_{p+.,p}: faces of X^{p,r} of cubic dimension exactly p, graded by
    simplicial dimension.  This is E^0_{p,.} at n = p."""
    return _cached(("loc", p, r), lambda: all_faces(
        p, lambda S: cdim(p, S) == p and diam(S) <= r))


def U_complex(p, r):
    """X^{p,r}_{p-1}: faces of cubic dimension < p."""
    return _cached(("U", p, r), lambda: all_faces(
        p, lambda S: cdim(p, S) < p and diam(S) <= r))


def halfcube_cover(p):
    """The 2p sets F_{i,eps} = {x : x_i = eps} of Lemma 3.1, keyed (i, eps)."""
    N = 1 << p
    cov = {}
    for i in range(p):
        for eps in (0, 1):
            cov[(i, eps)] = frozenset(
                v for v in range(N) if ((v >> i) & 1) == eps)
    return cov


def augment(bydeg):
    """Add the empty simplex in degree -1, so that the homology of the result is
    the REDUCED homology of the original complex."""
    out = {k: list(v) for k, v in bydeg.items()}
    out[-1] = [()]
    return out


# --- the cross-polytope boundary Diamond_p ---------------------------------

def diamond(p):
    """Section 3: the complex on the 2p symbols (i,eps), encoded 2*i+eps, whose
    faces are the subsets carrying at most one symbol per i."""
    bydeg = {}
    for mask in range(1, 1 << (2 * p)):
        S = tuple(v for v in range(2 * p) if (mask >> v) & 1)
        seen = set()
        good = True
        for v in S:
            if v // 2 in seen:
                good = False
                break
            seen.add(v // 2)
        if good:
            bydeg.setdefault(len(S) - 1, []).append(S)
    return bydeg


# --- the hyperoctahedral group ---------------------------------------------

def hyperoctahedral(p):
    """H_p = (Z/2)^p semidirect S_p, as a list of (s, sigma) with s a tuple of
    p bits and sigma a tuple with sigma[i] = image of i."""
    out = []
    for sigma in itertools.permutations(range(p)):
        for s in itertools.product((0, 1), repeat=p):
            out.append((s, sigma))
    return out


def perm_sign(seq):
    """Sign of the permutation sorting `seq` (a sequence of distinct items)."""
    a = list(seq)
    n = len(a)
    sign = 1
    for i in range(n):
        m = i
        for j in range(i + 1, n):
            if a[j] < a[m]:
                m = j
        if m != i:
            a[i], a[m] = a[m], a[i]
            sign = -sign
    return sign


def det_char(g):
    """det of the signed permutation matrix: sgn(sigma) * (-1)^{sum s}."""
    s, sigma = g
    return perm_sign(list(sigma)) * (-1) ** (sum(s) % 2)


def vertex_action(p, g):
    """Section 1 and the proof of Lemma 4.1: (g x)_{sigma(i)} = x_i + s_{sigma(i)}.
    -> a list vmap with vmap[v] = g.v on {0,1}^p encoded as an int."""
    s, sigma = g
    N = 1 << p
    vmap = [0] * N
    for v in range(N):
        w = 0
        for i in range(p):
            bit = (v >> i) & 1
            j = sigma[i]
            w |= (bit ^ s[j]) << j
        vmap[v] = w
    return vmap


def symbol_action(p, g):
    """Section 3: (i,eps) -> (sigma(i), eps + s_{sigma(i)}), encoded 2i+eps."""
    s, sigma = g
    out = [0] * (2 * p)
    for i in range(p):
        for eps in (0, 1):
            j = sigma[i]
            out[2 * i + eps] = 2 * j + (eps ^ s[j])
    return out


# ===========================================================================
# G1  the collapse:  Lemma 3.1
# ===========================================================================

def g1_collapse():
    note("G1  Lemma 3.1, the collapse:  C^{p,r}_{p+.,p} = C_{p+.}(Delta, U)")
    for p in (1, 2, 3, 4):
        r = p                            # the extreme in-range case p = r
        N = 1 << p
        # (a) the Rips condition is vacuous on V(Q_p) because diam Q_p = p <= r
        worst = max(popcount(a ^ b) for a in range(N) for b in range(N))
        check("collapse-rips-vacuous-p%d" % p, worst, r,
              "diam Q_%d = %d <= r = %d, so every non-empty subset is a face"
              % (p, worst, r))
        # (b) cdim = p  <=>  full support  <=>  in Delta but not in U
        delta = rips_faces(p, r)
        U = U_complex(p, r)
        loc = local_complex(p, r)
        dset = set(s for v in delta.values() for s in v)
        uset = set(s for v in U.values() for s in v)
        lset = set(s for v in loc.values() for s in v)
        allsub = set(tuple(v for v in range(N) if (m >> v) & 1)
                     for m in range(1, 1 << N))
        full = set(S for S in allsub
                   if all(any(((v >> i) & 1) != ((S[0] >> i) & 1) for v in S)
                          for i in range(p)))
        check("collapse-cdim-is-full-support-p%d" % p,
              (dset == allsub, lset == full, lset == dset - uset), (True,) * 3,
              "|Delta|=%d |U|=%d |Delta\\U|=%d" % (len(dset), len(uset), len(lset)))
        # (c) U is exactly the union of the 2p half-cube simplices
        cov = halfcube_cover(p)
        union = set()
        for F in cov.values():
            for S in allsub:
                if set(S) <= F:
                    union.add(S)
        check("collapse-U-is-union-of-2p-halfcubes-p%d" % p, uset, union,
              "2p = %d half-cube simplices, |U| = %d faces" % (2 * p, len(uset)))
        # (d) the differentials agree: the relative complex has the simplicial
        #     boundary, i.e. dropping a vertex of a full-support face gives
        #     either a full-support face or a face of U.
        okdiff = True
        for S in full:
            for a in range(len(S)):
                T = S[:a] + S[a + 1:]
                if T and not (T in full or T in uset):
                    okdiff = False
        check("collapse-differential-p%d" % p, okdiff, True,
              "every codim-1 face of a full-support face is full-support or in U")


# ===========================================================================
# G2  the counts of Remark 3.2
# ===========================================================================

def g2_counts():
    note("G2  Remark 3.2, the counts at p = 4 and the dimension vector (3.2)")
    for p in (1, 2, 3, 4):
        N = 1 << p
        U = U_complex(p, p)
        loc = local_complex(p, p)
        nu = sum(len(v) for v in U.values())
        nl = sum(len(v) for v in loc.values())
        check("count-complement-p%d" % p, nu + nl, (1 << N) - 1,
              "|U| = %d, |full support| = %d, 2^(2^%d) - 1 = %d"
              % (nu, nl, p, (1 << N) - 1))
    # the inclusion-exclusion printed in Remark 3.2
    parts = [a * b for a, b in PAPER_IE_TERMS]
    check("count-ie-partials", parts, PAPER_IE_PARTIALS,
          "8*255, 24*15, 32*3, 16*1 as printed")
    check("count-ie-alternating", parts[0] - parts[1] + parts[2] - parts[3],
          PAPER_U4_FACES, "2040 - 360 + 96 - 16 = 1760")
    U4 = U_complex(4, 4)
    check("count-U4-faces", sum(len(v) for v in U4.values()), PAPER_U4_FACES)
    # and the inclusion-exclusion is the right one: multiplicities and simplex
    # sizes of the 8-fold half-cube cover of Q_4
    cov = halfcube_cover(4)
    sizes = {}
    for k in range(1, 5):
        for combo in itertools.combinations(sorted(cov), k):
            inter = frozenset.intersection(*[cov[c] for c in combo])
            if inter:
                sizes.setdefault(k, []).append(len(inter))
    got = [(len(sizes[k]), (1 << sizes[k][0]) - 1) for k in range(1, 5)]
    check("count-ie-cover-structure", got, PAPER_IE_TERMS,
          "k-fold non-empty intersections of the 8 half-cubes, and 2^|.|-1")
    # equation (3.2)
    loc44 = local_complex(4, 4)
    dims = [len(loc44.get(k, [])) for k in range(1, 16)]
    check("dims-3.2", dims, PAPER_DIMS_44, "degrees 1..15 of C^{4,4}_{4+.,4}")
    check("dims-3.2-total", sum(dims), PAPER_TOTAL_44)
    check("dims-3.2-complement", PAPER_DELTA4_FACES - PAPER_U4_FACES,
          PAPER_TOTAL_44, "65535 - 1760")
    check("dims-3.2-euler", sum((-1) ** k * d for k, d in
                                zip(range(1, 16), dims)), PAPER_EULER_44,
          "sum_k (-1)^k dim = +1")


# ===========================================================================
# G3  the cross-polytope nerve and the cycle z
# ===========================================================================

def g3_nerve():
    note("G3  Section 3, Diamond_p = (S^0)^{*p}, the cycle z, and its character")
    for p in range(1, 8):
        D = diamond(p)
        f = [len(D.get(s, [])) for s in range(p)]
        want = [math.comb(p, s + 1) * 2 ** (s + 1) for s in range(p)]
        check("nerve-fvector-p%d" % p, f, want,
              "faces of dimension 0..%d of Diamond_%d" % (p - 1, p))
        dims, ranks, betti = ranks_and_betti(augment(D), P1)
        red = {k: v for k, v in betti.items() if k >= 0}
        char0_pin("nerve-p%d" % p, dims, betti, p - 1)
        check("nerve-is-sphere-p%d" % p,
              sorted(k for k, v in red.items() if v), [p - 1],
              "reduced homology of Diamond_%d = that of S^%d" % (p, p - 1))
    # Diamond_p really is the nerve of the 2p half-cube cover of U_p
    for p in (2, 3, 4):
        cov = halfcube_cover(p)
        keys = sorted(cov)
        nerve = set()
        for k in range(1, 2 * p + 1):
            for combo in itertools.combinations(keys, k):
                if frozenset.intersection(*[cov[c] for c in combo]):
                    nerve.add(tuple(sorted(2 * i + e for (i, e) in combo)))
        D = diamond(p)
        check("nerve-is-the-cover-nerve-p%d" % p,
              nerve, set(s for v in D.values() for s in v),
              "%d non-empty intersection patterns" % len(nerve))
        # and every non-empty intersection is a full simplex, hence contractible
        allfull = True
        for k in range(1, 2 * p + 1):
            for combo in itertools.combinations(keys, k):
                inter = frozenset.intersection(*[cov[c] for c in combo])
                if inter and len(inter) != (1 << max(0, p - len(
                        set(i for (i, _) in combo)))):
                    allfull = False
        check("nerve-intersections-are-simplices-p%d" % p, allfull, True,
              "each non-empty intersection is a full simplex on 2^(p-|I|) "
              "vertices, so the nerve theorem applies")
    # the printed z, term for term
    gen = [((-1) ** sum(int(c) for c in w), w)
           for w in ("".join(t) for t in
                     itertools.product("01", repeat=4))]
    check("z-table-matches-formula-3.3", PAPER_Z_P4, gen,
          "the 16 printed terms are exactly (-1)^{eps_1+...+eps_4}")
    for p in range(1, 7):
        z = z_chain(p)
        check("z-is-a-cycle-p%d" % p, boundary_of(z), {},
              "dz = 0 in C_{%d}(Diamond_%d), %d terms" % (p - 1, p, len(z)))
        D = augment(diamond(p))
        top = D[p - 1]
        # p-1 is the top dimension of Diamond_p, so H_{p-1} = ker d_{p-1}
        # (computed on the AUGMENTED complex, so that p = 1 gives reduced H_0)
        rows = boundary_rows(D, p - 1, P1)
        r = rank_mod(rows, P1)
        check("z-spans-top-homology-p%d" % p, len(top) - r, 1,
              "dim ker d_%d = %d - %d = 1, and z is a non-zero element of it"
              % (p - 1, len(top), r))
    for p in range(1, 6):
        G = hyperoctahedral(p)
        z = z_chain(p)
        allok = True
        vals = set()
        for g in G:
            sm = symbol_action(p, g)
            img = {}
            for facet, c in z.items():
                tgt = [sm[v] for v in facet]
                key = tuple(sorted(tgt))
                img[key] = img.get(key, 0) + c * perm_sign(tgt)
            img = {k: v for k, v in img.items() if v}
            d = det_char(g)
            vals.add(d)
            if img != {k: d * v for k, v in z.items()}:
                allok = False
        check("z-character-is-det-p%d" % p, (allok, sorted(vals)),
              (True, [-1, 1] if p >= 1 else [1]),
              "g.z = det(g) z on all %d elements of H_%d" % (len(G), p))


def z_chain(p):
    """Equation (3.3) as a dict {sorted facet: coefficient}."""
    z = {}
    for eps in itertools.product((0, 1), repeat=p):
        facet = tuple(2 * i + eps[i] for i in range(p))
        z[facet] = (-1) ** (sum(eps) % 2)
    return z


def boundary_of(chain):
    """Simplicial boundary of a dict {sorted face: coeff}, vertices ordered by
    their integer code."""
    out = {}
    for S, c in chain.items():
        for a in range(len(S)):
            T = S[:a] + S[a + 1:]
            out[T] = out.get(T, 0) + ((-1) ** a) * c
    return {k: v for k, v in out.items() if v}


# ===========================================================================
# G4  Lemma 4.1: the homology of the pair, both routes, and their agreement
# ===========================================================================

def g4_pair():
    note("G4  Lemma 4.1 and Remark 4.2(a): H_k(Delta,U) = k in degree p only")
    for p in (1, 2, 3, 4):
        # route (a): reduced homology of U, then the long exact sequence
        U = U_complex(p, p)
        dimsU, ranksU, bettiU = ranks_and_betti(augment(U), P1)
        char0_pin("U-p%d" % p, dimsU, bettiU, p - 1)
        red = sorted(k for k, v in bettiU.items() if v and k >= 0)
        check("U-is-a-sphere-p%d" % p, red, [p - 1],
              "reduced homology of U_%d concentrated in degree %d "
              "(U_%d has %d faces)"
              % (p, p - 1, p, sum(len(v) for v in U.values())))
        # route (b): the relative complex itself
        loc = local_complex(p, p)
        dims, ranks, betti = ranks_and_betti(loc, P1)
        char0_pin("local-p%d" % p, dims, betti, p)
        check("local-homology-p%d" % p,
              sorted(k for k, v in betti.items() if v), [p],
              "H_k(C^{%d,%d}_{%d+.,%d}) = k for k = %d only, %d generators"
              % (p, p, p, p, p, sum(len(v) for v in loc.values())))
        # the two routes agree: H_k(Delta,U) = reduced H_{k-1}(U)
        shifted = {k: bettiU.get(k - 1, 0) for k in betti}
        check("two-routes-agree-p%d" % p, betti, shifted,
              "H_k(Delta,U) = reduced H_{k-1}(U) for every k")
        # a second prime, on the cells small enough for it
        if sum(len(v) for v in loc.values()) <= 6000:
            _, _, b2 = ranks_and_betti(loc, P2)
            check("local-second-prime-p%d" % p, b2, betti,
                  "GF(%d) agrees with GF(%d)" % (P2, P1))
    # the character of the local homology, computed EQUIVARIANTLY by the
    # Hopf-Lefschetz trace on the chain level, at the target cell p = r = 4
    lefschetz_character("target-4-4", 4, 4, expect_degrees={4: "det"})


def lefschetz_character(tag, p, r, expect_degrees):
    """L(g) = sum_k (-1)^k tr(g | C_k) = sum_k (-1)^k tr(g | H_k).

    Only g-invariant faces contribute to tr(g | C_k), and a g-invariant face is
    a union of orbits of g on V(Q_p); so the trace is computed by enumerating
    unions of orbits, never all faces.  For g = 1 the Lefschetz number is the
    Euler characteristic, which is taken from the dimension vector instead.
    """
    G = hyperoctahedral(p)
    loc = local_complex(p, r)
    dims = {k: len(v) for k, v in loc.items()}
    e = euler(dims)
    predicted = {}
    for g in G:
        vmap = vertex_action(p, g)
        if all(vmap[v] == v for v in range(1 << p)):
            predicted[g] = e
            continue
        # orbits of g on V
        seen = [False] * (1 << p)
        orbits = []
        for v in range(1 << p):
            if seen[v]:
                continue
            o = []
            w = v
            while not seen[w]:
                seen[w] = True
                o.append(w)
                w = vmap[w]
            orbits.append(tuple(sorted(o)))
        tot = 0
        m = len(orbits)
        for mask in range(1, 1 << m):
            S = []
            for i in range(m):
                if (mask >> i) & 1:
                    S.extend(orbits[i])
            S = tuple(sorted(S))
            if cdim(p, S) != p:
                continue
            if r < p and diam(S) > r:
                continue
            tot += ((-1) ** (len(S) - 1)) * perm_sign([vmap[v] for v in S])
        predicted[g] = tot
    # what the paper's answer predicts
    def want(g):
        s = 0
        for k, kind in expect_degrees.items():
            s += ((-1) ** k) * (det_char(g) if kind == "det" else 1)
        return s
    allok = all(predicted[g] == want(g) for g in G)
    check("lefschetz-character-" + tag, allok, True,
          "sum_k (-1)^k tr(g|C_k) = %s on all %d elements of H_%d"
          % (" + ".join("(-1)^%d*%s" % (k, v)
                        for k, v in sorted(expect_degrees.items())),
             len(G), p))


# ===========================================================================
# G5  the cube's own chain complex, and the flag transitivity of D4
# ===========================================================================

def cells(n, p):
    """The p-cells of the n-cube as words in {0,1,2}^n with exactly p twos
    (2 = a free coordinate)."""
    out = []
    for A in itertools.combinations(range(n), p):
        for vals in itertools.product((0, 1), repeat=n - p):
            w = [None] * n
            it = iter(vals)
            for i in range(n):
                w[i] = 2 if i in A else next(it)
            out.append(tuple(w))
    return out


def cube_boundary(w):
    """The cellular boundary of a cell, as {facet: coefficient}."""
    stars = [i for i, x in enumerate(w) if x == 2]
    out = {}
    for j, i in enumerate(stars):
        for eps in (0, 1):
            f = list(w)
            f[i] = eps
            out[tuple(f)] = ((-1) ** j) * (1 if eps == 1 else -1)
    return out


def cube_group_action(n, g, w):
    s, sigma = g
    out = [None] * n
    for i in range(n):
        j = sigma[i]
        out[j] = 2 if w[i] == 2 else (w[i] ^ s[j])
    return tuple(out)


def g5_cube():
    note("G5  Proposition 2.3 and step D4: the cube's chain complex at n = 6")
    n = PAPER_N
    for p in range(0, 5):
        C = cells(n, p)
        check("cube-dim-F%d" % p, len(C), PAPER_DIM_F[p],
              "C(%d,%d) 2^%d = %d p-cells of Q_%d" % (n, p, n - p, len(C), n))
    # dim F_p is the index [H_n : S_{n-p} x H_p] of Proposition 2.3
    Hn = 2 ** n * math.factorial(n)
    for p in range(0, 5):
        sub = math.factorial(n - p) * (2 ** p * math.factorial(p))
        check("cube-index-F%d" % p, Hn // sub, PAPER_DIM_F[p],
              "[H_%d : S_%d x H_%d] = %d / %d" % (n, n - p, p, Hn, sub))
    check("cube-dim-F4-is-60", PAPER_DIM_F[4], PAPER_INDEX_60,
          "dim E^1_{4,0}(X^{6,4}) = dim F_4 = 60")
    # every p-cell has exactly 2p facets, all coefficients units, and dd = 0
    for p in range(1, 5):
        C = cells(n, p)
        allok = all(len(cube_boundary(w)) == 2 * p
                    and set(abs(v) for v in cube_boundary(w).values()) == {1}
                    for w in C)
        check("cube-boundary-units-p%d" % p, allok, True,
              "each of the %d %d-cells has %d facets, coefficients +-1"
              % (len(C), p, 2 * p))
        dd = {}
        for w in C:
            for f, c in cube_boundary(w).items():
                for h, c2 in cube_boundary(f).items():
                    dd[(w, h)] = dd.get((w, h), 0) + c * c2
        check("cube-boundary-squares-to-zero-p%d" % p,
              all(v == 0 for v in dd.values()), True,
              "d_%d d_%d = 0" % (p - 1, p))
    # D4: H_n is transitive on flags (K' facet of K)
    gens = []
    for i in range(n - 1):
        sig = list(range(n))
        sig[i], sig[i + 1] = sig[i + 1], sig[i]
        gens.append((tuple([0] * n), tuple(sig)))
    for i in range(n):
        s = [0] * n
        s[i] = 1
        gens.append((tuple(s), tuple(range(n))))
    for p in range(1, 5):
        C = cells(n, p)
        flags = set((w, f) for w in C for f in cube_boundary(w))
        start = next(iter(sorted(flags)))
        seen = {start}
        frontier = [start]
        while frontier:
            (w, f) = frontier.pop()
            for g in gens:
                nxt = (cube_group_action(n, g, w), cube_group_action(n, g, f))
                if nxt not in seen:
                    seen.add(nxt)
                    frontier.append(nxt)
        check("cube-flag-transitive-p%d" % p, (len(seen), len(flags)),
              (len(flags), len(flags)),
              "H_%d is transitive on the %d flags (K' facet of K)"
              % (n, len(flags)))


# ===========================================================================
# G6  controls, both polarities, on the source's own published integers
# ===========================================================================

def g6_controls():
    note("G6  controls on the published integers of Galetto-Montano-Wellner")
    for (p, r), (first, pub) in sorted(PAPER_PUBLISHED_DIMS.items()):
        loc = local_complex(p, r)
        got = [len(loc.get(k, [])) for k in range(first, first + len(pub))]
        check("published-dims-p%d-r%d" % (p, r), got, pub,
              "degrees %d..%d of C^{%d,%d}_{.,%d}"
              % (first, first + len(pub) - 1, p, r, p))
        check("published-dims-p%d-r%d-exhaustive" % (p, r),
              sorted(loc), list(range(first, first + len(pub))),
              "no degree outside the published range is occupied")
        dims, ranks, betti = ranks_and_betti(loc, P1)
        _, _, b2 = ranks_and_betti(loc, P2)
        check("published-two-primes-p%d-r%d" % (p, r), b2, betti,
              "GF(%d) agrees with GF(%d)" % (P2, P1))
        nz = sorted(k for k, v in betti.items() if v)
        if p <= r:
            char0_pin("control-p%d-r%d" % (p, r), dims, betti, p)
            check("control-in-range-silent-p%d-r%d" % (p, r), nz, [p],
                  "in range p <= r: homology in degree %d only, as published"
                  % p)
        else:
            check("control-off-range-FIRES-p%d-r%d" % (p, r), nz,
                  PAPER_OFFRANGE_43_DEGREES,
                  "p = %d > r = %d: the mod-P Betti vector is non-zero in TWO "
                  "degrees, reproducing the anomaly published at source line "
                  "2309 -- the negative control fires. NB a mod-P Betti number "
                  "only BOUNDS the characteristic-zero one from above, so this "
                  "reproduces the published answer, it does not derive it"
                  % (p, r))
            check("control-off-range-euler-p%d-r%d" % (p, r), euler(dims), 0,
                  "and the alternating sum is 0, consistent with one class in "
                  "degree 4 and one in degree 7")
        if p <= 3 or sum(len(v) for v in loc.values()) <= 6000:
            kinds = ({4: "det", 7: "triv"} if (p, r) == (4, 3)
                     else {p: "det"})
            lefschetz_character("control-p%d-r%d" % (p, r), p, r, kinds)
    # the sharper published negative: at (p,r) = (3,2) it is bullet (i) that
    # fails, and H_3 is 2-dimensional
    loc32 = local_complex(3, 2)
    dims32, ranks32, betti32 = ranks_and_betti(loc32, P1)
    _, _, b32b = ranks_and_betti(loc32, P2)
    check("offrange-3-2-two-primes", b32b, betti32,
          "GF(%d) agrees with GF(%d)" % (P2, P1))
    check("offrange-3-2-bullet-i-FAILS", betti32.get(3), PAPER_OFFRANGE_32_DIM_H3,
          "dim_GF(P) H_3(C^{3,2}_{.,3}) = 2 > 1, over both primes: bullet (i) "
          "fails off range while bullet (ii) survives. The characteristic-zero "
          "value 2 and the splitting {1^3;0}+{0;1^3} are published at source "
          "lines 1760-1780 and are reproduced, not derived, here")
    check("offrange-3-2-bullet-ii-survives",
          sorted(k for k, v in betti32.items() if v), [3],
          "the vanishing for q != 0 still holds at (p,r) = (3,2)")
    # and the boundary case of the truncation: at p = r+1 the minimal
    # full-support subsets are the antipodal pairs, at distance p > r
    for r in (2, 3, 4):
        p = r + 1
        N = 1 << p
        anti = [(v, N - 1 - v) for v in range(N // 2)]
        check("truncation-boundary-r%d" % r,
              (len(anti), sorted(set(popcount(a ^ b) for a, b in anti)),
               p > r), (N // 2, [p], True),
              "at p = r+1 = %d the %d antipodal pairs are the minimal "
              "full-support subsets and sit at distance %d > r = %d"
              % (p, N // 2, p, r))


# ===========================================================================
# main
# ===========================================================================

def main():
    print("verify.py -- the first r+1 columns of the cubic-dimension spectral "
          "sequence")
    print("exact integer arithmetic only; ranks over GF(%d) and GF(%d)"
          % (P1, P2))
    print("")
    g1_collapse()
    g2_counts()
    g3_nerve()
    g4_pair()
    g5_cube()
    g6_controls()

    print("")
    scope("NOT RE-RUN: the two results quoted from Galetto-Montano-Wellner and "
          "not reproved in the paper -- the n-freeness of the p-th column "
          "(source lines 930-939, label pro:1) and the induced description of "
          "F_i with its determinant character (lines 1029-1047, lem:1). Every "
          "check above is therefore made at n = p, which is exactly the case "
          "those results reduce to; nothing here verifies the reduction "
          "itself.")
    scope("NOT RE-RUN: the standard topology the proof cites -- the long exact "
          "sequence of a pair, the nerve theorem, the Mayer-Vietoris spectral "
          "sequence of a cover, and the Kunneth formula for a join. What is "
          "checked is that their CONCLUSIONS agree with direct computation at "
          "p <= 4 (check two-routes-agree-p*), not the theorems.")
    scope("NOT RE-RUN: Lemma 4.1 itself. The marginal map argument is a "
          "homotopy of topological pairs and is not modelled here; what is "
          "checked is its numerical consequence, degree by degree, at p <= 4, "
          "plus the H_p-character of the answer.")
    scope("NOT RE-RUN: any cell with p = r >= 5. That would need 2^32 subsets "
          "of V(Q_5). The proof of the paper is uniform in n, r and p; the "
          "machine evidence here stops at p = 4, and the equivariant character "
          "check stops at p = 5 (nerve) and p = 4 (chain level).")
    scope("NOT RE-RUN: the flag-transitivity and unit-coefficient steps of D4 "
          "at general n. They are checked only at n = 6 and p <= 4.")
    scope("NOT RE-RUN: the characteristic-zero labelling {lambda;mu} of "
          "H_m-irreducibles. This program computes characters as integer "
          "class functions and compares them with det and with the trivial "
          "character; it does not decompose any representation into the "
          "source's labelled irreducibles.")
    scope("NOT RE-RUN: the two off-range facts of Section 6(3) as statements "
          "about H_m-isomorphism TYPES. Their dimensions and Lefschetz "
          "characters are reproduced; the identifications {1^3;0}+{0;1^3} at "
          "(3,2) and {4;0} at (4,3) are quoted from the source.")
    scope("NOT RE-RUN: the literature. This program performs no search, "
          "fetches nothing, and checks no citation.")
    for s in _SCOPE:
        print(s)
        print("")

    n = _STATE["pass"]
    bad_n = _STATE["fail"]
    print("[%.1fs] %d checks ran, %d failed" % (time.time() - T0, n + bad_n, bad_n))
    if bad_n:
        print("VERDICT: %d OF %d CHECKS FAILED" % (bad_n, n + bad_n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())

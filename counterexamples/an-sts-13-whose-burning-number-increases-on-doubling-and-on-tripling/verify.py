#!/usr/bin/env python3
"""Verification program for the accompanying paper.

It reads the object printed in the paper -- the 26 blocks of the Steiner triple system H of
order 13 -- and re-derives every quantity the paper asserts.  Python 3.9 or later, standard
library only, no external data file.  All arithmetic is exact integer arithmetic; in
particular the ceiling of a base-two logarithm is computed with `int.bit_length()` and never
with a floating-point log, so no decision in this program depends on a float.

Two independent routes to the conclusion are both re-derived here:

  ROUTE A (self-contained).  The burning number is computed DIRECTLY from the definition of
  the burning process, by a complete breadth-first enumeration of every reachable burnt set,
  round by round, deduplicated by bitmask.  The deduplication is sound because the burnt set
  F(r) is a sufficient state: both the next propagation step and the legality of the next
  source depend on F(r) alone.  The enumeration also permits the arsonist to pass, so a
  lower bound produced here holds for generalized burning sequences as well as ordinary ones.
  This route uses no theorem of the target paper at all.

  ROUTE B (the paper's own two unconditional bounds).  The recursion h(1) = 1,
  h(r) = C(h(r-1),2) - 2h(r-1) + 3r - 2 is recomputed and compared with the nine integers
  the target prints; rho_v (least r with v <= h(r)) and U(v) = ceil(log2 v) + 1 are then
  evaluated exactly.

Every published integer this program tests itself against belongs to somebody else: the three
burning numbers the target computes at orders 3, 7 and 9, its two closed forms
b(PG(n,2)) = n+2 and b(AG(m,3)) = m+3 evaluated at n = 2, 3 and m = 2, 3, its printed table
of h, and the classical fact that tripling AG(2,3) yields AG(3,3).  (The target also computes
b at order 31; that order is out of reach of the enumeration here and is not used.)

Run:  python3 verify.py
"""

import sys
from itertools import combinations, permutations

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append(bool(ok))
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, (" " + detail) if detail else ""))


# ---------------------------------------------------------------------------
# 0.  THE OBJECT, EXACTLY AS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
# The 26 blocks of H, transcribed from the display in Section 2 of the paper.
H_PRINTED = [
    (0, 1, 4), (0, 2, 7), (0, 3, 12), (0, 5, 11), (0, 6, 8), (0, 9, 10),
    (1, 2, 5), (1, 3, 8), (1, 6, 12), (1, 7, 9), (1, 10, 11),
    (2, 3, 6), (2, 4, 9), (2, 8, 10), (2, 11, 12),
    (3, 4, 7), (3, 5, 10), (3, 9, 11),
    (4, 5, 8), (4, 6, 11), (4, 10, 12),
    (5, 6, 9), (5, 7, 12),
    (6, 7, 10),
    (7, 8, 11),
    (8, 9, 12),
]
# The burning sequence printed in the paper.
SEQ_PRINTED = (0, 1, 2, 3, 6)
# The nine values of h printed in the target paper's own table.
H_TABLE_PRINTED = [1, 2, 4, 8, 25, 266, 34732, 603069104, 181846170592008673]


# ---------------------------------------------------------------------------
# 1.  STEINER TRIPLE SYSTEMS
# ---------------------------------------------------------------------------
def is_sts(v, blocks):
    """True iff `blocks` is the block set of a Steiner triple system of order v: every
    unordered pair of the v points lies in exactly one block."""
    seen = set()
    for B in blocks:
        if len(set(B)) != 3:
            return False
        for p in B:
            if not (0 <= p < v):
                return False
        for pr in combinations(sorted(B), 2):
            if pr in seen:
                return False
            seen.add(pr)
    return len(seen) == v * (v - 1) // 2 and len(blocks) == v * (v - 1) // 6


def norm(blocks):
    return sorted(tuple(sorted(B)) for B in blocks)


def cyclic_sts13():
    """H as the development of the base blocks {0,1,4} and {0,2,7} modulo 13."""
    return norm(((a + i) % 13, (b + i) % 13, (c + i) % 13)
                for (a, b, c) in [(0, 1, 4), (0, 2, 7)] for i in range(13))


def one_factors_k14():
    """The 1-factorization of the complete graph on Z_13 u {inf} used by the paper:
    F_i = {{inf,i}} u {{i+j, i-j} : j = 1..6}.  `inf` is the token 13."""
    F = []
    for i in range(13):
        f = [(min(13, i), max(13, i))]
        for j in range(1, 7):
            a, b = (i + j) % 13, (i - j) % 13
            f.append((min(a, b), max(a, b)))
        F.append(f)
    return F


def double_sts(v, blocks, factors, w_of, w_inf):
    """Double-and-add-one: an STS(v) plus a 1-factorization of K_{v+1} on the new points,
    giving an STS(2v+1).  `factors[i]` is the 1-factor indexed by the old point i."""
    D = norm(blocks)
    for i in range(v):
        for (a, b) in factors[i]:
            D.append(tuple(sorted((i, w_inf if a == v else w_of(a),
                                   w_inf if b == v else w_of(b)))))
    return sorted(D)


def triple_sts(v, blocks):
    """The tripling construction: points V x Z_3, written x + v*i.  Blocks are the v vertical
    triples {(x,0),(x,1),(x,2)} together with, for each block {x,y,z} and each
    (i,j,k) in Z_3^3 with i+j+k = 0, the triple {(x,i),(y,j),(z,k)}."""
    P = lambda x, i: x + v * (i % 3)
    T = [tuple(sorted((P(x, 0), P(x, 1), P(x, 2)))) for x in range(v)]
    for (x, y, z) in blocks:
        for i in range(3):
            for j in range(3):
                T.append(tuple(sorted((P(x, i), P(y, j), P(z, (-i - j) % 3)))))
    return sorted(T)


def affine_sts(m):
    """AG(m,3): the points of F_3^m, blocks the triples of distinct points summing to zero.
    Points are labelled by the base-3 digit string read with digit t as the coefficient of
    3^t, which for m = 2,3 makes the labelling of AG(3,3) agree with V x Z_3 for
    V = AG(2,3) -- exactly the identification the tripling construction uses."""
    pts = []
    for code in range(3 ** m):
        d, c = [], code
        for _ in range(m):
            d.append(c % 3)
            c //= 3
        pts.append(tuple(d))
    lab = {p: i for i, p in enumerate(pts)}
    B = set()
    for a, b in combinations(pts, 2):
        c = tuple((-a[t] - b[t]) % 3 for t in range(m))
        if c != a and c != b:
            B.add(tuple(sorted((lab[a], lab[b], lab[c]))))
    return sorted(B)


def projective_sts(n):
    """PG(n,2): the nonzero vectors of F_2^{n+1}, blocks {x, y, x+y}.  Point x is labelled
    x-1 so that the point set is {0,...,2^{n+1}-2}."""
    N = 1 << (n + 1)
    B = set()
    for x in range(1, N):
        for y in range(x + 1, N):
            B.add(tuple(sorted((x - 1, y - 1, (x ^ y) - 1))))
    return sorted(B)


# ---------------------------------------------------------------------------
# 2.  THE BURNING PROCESS, STRAIGHT FROM THE DEFINITION
# ---------------------------------------------------------------------------
def block_masks(blocks):
    return [(1 << a) | (1 << b) | (1 << c) for (a, b, c) in blocks]


def spread(S, ms):
    """ONE propagation step, not a closure: the set of points not in S that lie in a block
    whose other two points are both in S."""
    new = 0
    for m in ms:
        r = m & ~S
        if r and not (r & (r - 1)):        # exactly one point of this block is unburnt
            new |= r
    return new


def burns_in(v, blocks, seq):
    """Simulate an ordinary burning sequence.  Returns (burnt_everything, trace)."""
    ms = block_masks(blocks)
    S = 0
    trace = []
    for u in seq:
        if (S >> u) & 1:                  # a source must be unburnt at the end of the
            return False, trace           # previous round
        S = S | spread(S, ms) | (1 << u)
        trace.append(sorted(i for i in range(v) if (S >> i) & 1))
    return S == (1 << v) - 1, trace


def burning_number(v, blocks, kmax):
    """Exact burning number by complete level-by-level enumeration of reachable burnt sets,
    or None if b > kmax.  The arsonist may also pass, so `None` is a lower bound valid for
    generalized burning sequences too.  Levels beyond kmax-1 are never expanded, which is
    what keeps the kmax = 5 lower-bound run at order 39 small."""
    ms = block_masks(blocks)
    FULL = (1 << v) - 1
    level = {0}
    sizes = []
    for k in range(1, kmax + 1):
        for S in level:
            rest = FULL & ~(S | spread(S, ms))
            if rest == 0 or not (rest & (rest - 1)):
                return k, sizes           # one more source (or a pass) finishes
        if k == kmax:
            return None, sizes
        nxt = set()
        for S in level:
            T = S | spread(S, ms)
            nxt.add(T)                    # the arsonist passes
            u = FULL & ~S
            while u:
                b = u & (-u)
                nxt.add(T | b)
                u ^= b
        level = nxt
        sizes.append(len(level))
    return None, sizes


# ---------------------------------------------------------------------------
# 3.  SUBSYSTEMS AND LAZY BURNING
# ---------------------------------------------------------------------------
def third_map(blocks):
    t = {}
    for B in blocks:
        for x, y in combinations(sorted(B), 2):
            t[(x, y)] = [p for p in B if p not in (x, y)][0]
    return t


def all_subsystems(v, blocks):
    """Every subset of the point set closed under `the third point of the block through a
    pair`.  Exhaustive over all 2^v subsets; v = 13 here, so 8192 subsets."""
    t = third_map(blocks)
    out = []
    for mask in range(1 << v):
        pts = [i for i in range(v) if (mask >> i) & 1]
        ok = True
        for x, y in combinations(pts, 2):
            if not ((mask >> t[(x, y)]) & 1):
                ok = False
                break
        if ok:
            out.append(frozenset(pts))
    return out


def closure(seed, t):
    S = set(seed)
    grew = True
    while grew:
        grew = False
        for x, y in combinations(sorted(S), 2):
            p = t[(x, y)]
            if p not in S:
                S.add(p)
                grew = True
    return frozenset(S)


def lazy_burning_number(points, blocks):
    """b_L = the least size of a subset whose closure is the whole point set.  (Lazy burning
    has no arsonist after round 1: the fire spreads to a fixpoint, so b_L is the least size
    of a generating set.)  Exhaustive over subsets in increasing size."""
    pts = sorted(points)
    sub = [B for B in blocks if set(B) <= set(pts)]
    t = third_map(sub)
    full = frozenset(pts)
    for k in range(1, len(pts) + 1):
        for C in combinations(pts, k):
            if closure(C, t) == full:
                return k
    return None


# ---------------------------------------------------------------------------
# 4.  THE PAPER'S TWO UNCONDITIONAL BOUNDS
# ---------------------------------------------------------------------------
def h_table(n):
    """h(1) = 1, h(r) = C(h(r-1),2) - 2h(r-1) + 3r - 2, exact integers."""
    h = [None, 1]
    for r in range(2, n + 1):
        p = h[r - 1]
        h.append(p * (p - 1) // 2 - 2 * p + 3 * r - 2)
    return h


HT = h_table(9)


def rho(v):
    """The least r with v <= h(r).  Lower bound of the target's Theorem `small_lem`."""
    r = 1
    while HT[r] < v:
        r += 1
    return r


def upper(v):
    """ceil(log2 v) + 1, exactly, with no float: ceil(log2 v) = (v-1).bit_length()."""
    return (v - 1).bit_length() + 1


def admissible(v):
    return v % 6 in (1, 3)


# ===========================================================================
#                                THE CHECKS
# ===========================================================================
print("verify.py -- an STS(13) whose burning number increases on doubling and on tripling")
print("python %s" % sys.version.split()[0])
print()

H = norm(H_PRINTED)

print("-- the objects --")
check("sts13-printed-blocks", is_sts(13, H),
      "the 26 blocks printed in the paper cover each of the 78 pairs of a 13-set exactly once")
check("sts13-cyclic-development", H == cyclic_sts13(),
      "and they are label-for-label the development of {0,1,4} and {0,2,7} modulo 13")

broken = [b for b in H if b != (0, 1, 4)] + [(0, 1, 5)]
check("sts13-checker-is-not-vacuous", not is_sts(13, broken),
      "the same test REJECTS H with the block {0,1,4} replaced by {0,1,5}")

F = one_factors_k14()
edges = [e for f in F for e in f]
ok_f = (len(F) == 13 and all(len(f) == 7 for f in F)
        and all(sorted(p for e in f for p in e) == list(range(14)) for f in F)
        and len(set(edges)) == 91 and len(edges) == 91)
check("one-factorization-k14", ok_f,
      "F_0..F_12 are 13 perfect matchings of K_14 on Z_13 u {inf} partitioning all 91 edges")

D = double_sts(13, H, F, lambda k: 13 + (k % 13), 26)
check("sts27-double", is_sts(27, D),
      "the 117 blocks of the double-and-add-one image cover each of the 351 pairs exactly once")
old = [B for B in D if max(B) <= 12]
new = [B for B in D if max(B) > 12]
check("sts27-restricts-to-H",
      norm(old) == H and len(new) == 91 and all(len([p for p in B if p <= 12]) == 1 for B in new),
      "D induces exactly H on {0,...,12}, and each of its 91 other blocks meets {0,...,12} once")

T = triple_sts(13, H)
check("sts39-triple", is_sts(39, T),
      "the 247 blocks of the tripling image cover each of the 741 pairs of a 39-set exactly once")
copies_ok = True
for i in range(3):
    inside = norm(tuple(p - 13 * i for p in B) for B in T if all(13 * i <= p < 13 * (i + 1) for p in B))
    if inside != H:
        copies_ok = False
check("sts39-three-copies-of-H", copies_ok,
      "T induces a copy of H on each of the three classes V x {0}, V x {1}, V x {2}")

AG2, AG3 = affine_sts(2), affine_sts(3)
check("tripling-control-affine", triple_sts(9, AG2) == AG3,
      "the same tripling rule applied to AG(2,3) returns AG(3,3) block for block, the "
      "classical fact the target relies on")

print()
print("-- H is non-degenerate --")
subs = all_subsystems(13, H)
orders = sorted(set(len(S) for S in subs))
check("subsystems-of-H-exhaustive",
      len(subs) == 41 and orders == [0, 1, 3, 13]
      and sum(1 for S in subs if len(S) == 3) == 26,
      "all 2^13 subsets tested: exactly 41 subsystems, of orders 0 (1), 1 (13), 3 (26), 13 (1)")
proper = [S for S in subs if len(S) < 13]
maximal = [S for S in proper if not any(S < X for X in proper)]
check("maximal-proper-subsystems",
      norm(tuple(sorted(S)) for S in maximal) == H,
      "the maximal proper subsystems of H are exactly its 26 blocks; H has no subsystem of "
      "order 7 or 9, so it is subsystem-free")
bL_H = lazy_burning_number(range(13), H)
check("lazy-burning-of-H", bL_H == 3,
      "b_L(H) = 3: no 2-subset of the 78 generates H, and some 3-subset does")
bL_max = sorted(set(lazy_burning_number(S, H) for S in maximal))
check("lazy-burning-of-maximal-subsystems", bL_max == [2],
      "every one of the 26 maximal proper subsystems has b_L = 2")
check("H-is-non-degenerate", bL_max == [bL_H - 1],
      "every maximal proper subsystem has b_L exactly b_L(H) - 1 = 2, which is the target's "
      "definition of non-degenerate")

print()
print("-- the LAZY burning number does not increase, so Zeitler's theorem is not contradicted --")
bL_D = lazy_burning_number(range(27), D)
check("lazy-consistency-doubling", bL_D == bL_H == 3,
      "b_L(D) = b_L(H) = 3: the LAZY burning number does not increase on doubling. Zeitler "
      "proved the projective systems are the only non-degenerate systems for which it does, "
      "and H is non-degenerate and not projective, so a lazy increase here would have been a "
      "contradiction rather than a result")
bL_T = lazy_burning_number(range(39), T)
check("lazy-consistency-tripling", bL_T == bL_H == 3,
      "b_L(T) = b_L(H) = 3 as well, so the separate lazy-burning conjecture of the target, "
      "which allows an increase on tripling only for affine systems, is also untouched")

print()
print("-- H is neither projective nor affine --")
check("order-13-not-projective",
      all(13 != (1 << (n + 1)) - 1 for n in range(1, 40)),
      "13 is not of the form 2^(n+1) - 1; the neighbouring values are 7 and 15")
check("order-13-not-affine",
      all(13 != 3 ** m for m in range(1, 40)),
      "13 is not of the form 3^m; the neighbouring values are 9 and 27")

print()
print("-- ROUTE B: the target's own two unconditional bounds --")
check("h-recursion-vs-printed-table", HT[1:10] == H_TABLE_PRINTED,
      "h(1..9) = %s, reproducing the nine integers the target prints" % (HT[1:10],))
check("rho-values", (rho(13), rho(27), rho(39)) == (5, 6, 6),
      "rho_13 = 5 (h(4) = 8 < 13 <= 25), rho_27 = 6 (h(5) = 25 < 27), rho_39 = 6")
check("upper-bound-values", (upper(13), upper(27), upper(39)) == (5, 6, 7),
      "ceil(log2 v) + 1 = 5, 6, 7 at v = 13, 27, 39, computed by bit_length with no float")
check("order-13-pinned", rho(13) == upper(13) == 5,
      "b(K) = 5 for EVERY Steiner triple system K of order 13")
check("order-27-pinned", rho(27) == upper(27) == 6,
      "b(K) = 6 for EVERY Steiner triple system K of order 27, hence for every 1-factorization "
      "used in the doubling")
check("routeB-doubling-increase", upper(13) < rho(27),
      "5 = b(H) < 6 = b(D): the burning number increases on doubling")
check("routeB-tripling-increase", upper(13) < rho(39),
      "5 = b(H) < 6 <= b(T): the burning number increases on tripling, though the bounds "
      "leave b(T) in {6,7} and so do not pin the increase to exactly one")

print()
print("-- ROUTE A: burning numbers computed from the definition, using no theorem of the target --")
b3, _ = burning_number(3, [(0, 1, 2)], 5)
check("burn-control-order-3", b3 == 3, "b(STS(3)) = 3, the value the target computes")
bf, _ = burning_number(7, projective_sts(2), 6)
check("burn-control-fano", bf == 4,
      "b(PG(2,2)) = 4, the value the target computes and its closed form n+2 at n = 2")
ba, _ = burning_number(9, AG2, 6)
check("burn-control-ag23", ba == 5,
      "b(AG(2,3)) = 5, the value the target computes and its closed form m+3 at m = 2")
bp, _ = burning_number(15, projective_sts(3), 6)
check("burn-control-pg32", bp == 5, "b(PG(3,2)) = 5, the target's closed form n+2 at n = 3")
bag3, _ = burning_number(27, AG3, 7)
check("burn-control-ag33", bag3 == 6, "b(AG(3,3)) = 6, the target's closed form m+3 at m = 3")

bH, szH = burning_number(13, H, 6)
check("burn-H-exact", bH == 5,
      "b(H) = 5 from the definition alone; no reachable burnt set after 3 rounds leaves at "
      "most one point unburnt (levels: %s)" % (szH,))
ok_seq, trace = burns_in(13, H, SEQ_PRINTED)
check("burn-H-printed-sequence", ok_seq and len(trace) == 5 and len(trace[-1]) == 13,
      "the sequence %s printed in the paper burns H; burnt sets have sizes %s"
      % (list(SEQ_PRINTED), [len(x) for x in trace]))
bD, szD = burning_number(27, D, 7)
check("burn-double-exact", bD == 6,
      "b(D) = 6 from the definition alone: every reachable burnt set after 4 rounds leaves "
      "at least two points unburnt (levels: %s)" % (szD,))
bT, szT = burning_number(39, T, 5)
check("burn-triple-lower", bT is None,
      "b(T) > 5 from the definition alone: none of the %d reachable burnt sets after 4 "
      "rounds leaves at most one point unburnt" % (szT[-1],))
check("routeA-doubling-increase", bH == 5 and bD == 6,
      "5 = b(H) < 6 = b(D), with no appeal to any theorem of the target paper")
check("routeA-tripling-increase", bH == 5 and bT is None,
      "5 = b(H) < b(T), with no appeal to any theorem of the target paper")

print()
print("-- the census that decides how much of the question survives --")
adm = [v for v in range(3, 10000) if admissible(v)]
pinned = [v for v in adm if rho(v) == upper(v)]
check("census-pinned-orders", pinned == [3, 7, 9, 13, 15, 27, 31],
      "over all %d admissible orders 3 <= v < 10000, the two bounds coincide exactly at "
      "v = 3, 7, 9, 13, 15, 27, 31" % len(adm))
fd = [v for v in adm if rho(2 * v + 1) > upper(v)]
check("census-forced-doubling", fd == [3, 7, 13, 15],
      "the bounds force an increase on doubling exactly at source orders 3, 7, 13, 15; "
      "3 and 7 are projective, so 13 and 15 are the only ones they leave")
ft = [v for v in adm if rho(3 * v) > upper(v)]
check("census-forced-tripling", ft == [3, 7, 9, 13, 15],
      "and on tripling exactly at 3, 7, 9, 13, 15; 3 and 7 are projective and 9 affine, so "
      "again only 13 and 15 remain")
anti = [(v, rho(v), upper(v)) for v in (19, 21, 25, 63, 81)]
check("census-anti-control", all(r != u for (_, r, u) in anti),
      "the same sandwich pins nothing at v = 19, 21, 25 (bounds [5,6]) nor at PG(5,2), "
      "v = 63 ([6,7]), nor at AG(4,3), v = 81 ([6,8]): %s" % (anti,))
HT20 = h_table(20)
check("census-induction-range", all(HT20[r] > 3 * 2 ** r for r in range(6, 21)),
      "h(r) > 3*2^r for r = 6..20, which is exactly the inequality the paper's Proposition "
      "proves for every r >= 6 by induction; its base case r = 6 is h(6) = %d against "
      "3*2^6 = %d, and it is what carries the previous three censuses from v < 10000 to every "
      "larger order" % (HT20[6], 3 * 2 ** 6))

print()
print("NOT RE-RUN: the target e-print is not fetched or parsed here. The verbatim question, its")
print("NOT RE-RUN: byte offset, the theorem labels and the nine printed values of h are quoted in")
print("NOT RE-RUN: the paper from the source; only the h values are re-derived above.")
print("NOT RE-RUN: Constructions 2.15 and 2.17 of the Colbourn-Dinitz Handbook are not consulted.")
print("NOT RE-RUN: The doubling and tripling rules are implemented from the standard statements")
print("NOT RE-RUN: reproduced in the paper, and are validated only by (i) the resulting designs")
print("NOT RE-RUN: being Steiner triple systems and (ii) the control triple(AG(2,3)) = AG(3,3).")
print("NOT RE-RUN: No isomorphism census is performed. The remark that there are exactly two")
print("NOT RE-RUN: STS(13) is neither used nor checked, and the non-projective non-degenerate")
print("NOT RE-RUN: STS(15) that would give a second witness order are not enumerated.")
print("NOT RE-RUN: b(T) is bounded below only; its exact value, 6 or 7, is not determined here.")
print("NOT RE-RUN: The separate observation about the target's Corollary with a missing rho >= 9")
print("NOT RE-RUN: hypothesis is not tested; nothing above depends on that Corollary.")
print()

n = len(CHECKS)
if all(CHECKS):
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)
print("VERDICT: %d of %d checks FAILED" % (sum(1 for c in CHECKS if not c), n))
sys.exit(1)

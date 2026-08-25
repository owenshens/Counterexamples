#!/usr/bin/env python3
"""Verification program for a refutation of the injectivity conjecture for
numerical labels on generalized Markov trees.

The generalized Markov equation, for k = (k1, k2, k3) in Z_{>=0}^3, is

    x^2 + y^2 + z^2 + k1*y*z + k2*z*x + k3*x*y = (3 + k1 + k2 + k3)*x*y*z.

For sigma = id the generalized Markov tree MT(k1,k2,k3,id) has root
((1,1), (k2+2,2), (1,3)); a node v = ((a,h),(b,i),(c,j)) has children

    L(v) = ((a,h), ((a^2 + k_j*a*b + b^2)/c, j), (b,i)),
    R(v) = ((b,i), ((b^2 + k_h*b*c + c^2)/a, h), (c,j)),

and it is synchronized with the Farey tree, whose root triple is
(0/1, 1/1, 1/0) and whose children of (r,t,s) are (r, r(+)t, t) and
(t, t(+)s, s) for the mediant (+). The middle pair of the node sitting over
the Farey triple (r,t,s) is written (n_t, i_t).

WHAT IS INPUT AND WHAT IS DERIVED.  The only things taken from the paper are
the exhibited data: the parameter triple (22,0,5); the Farey labels 1/4 and
2/3; the printed node values 15, 113, 889 and the printed substitution totals
3013710 and 800100; the family d_m = 30m+5, q_m = (F_{60m+11}-29)/10; and the
claimed labels/tags. Every tree, every node, every mutation, every Farey
label, every path, every Fibonacci number and every equation evaluation below
is recomputed here from the recursions above and only then compared with the
paper's printed data.

Exact arithmetic only: Python integers and fractions.Fraction. No floats, no
randomness, no third-party modules, no external data.
"""

from fractions import Fraction
import sys

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         ("  --  " + detail) if detail else ""))
    return bool(ok)


# --------------------------------------------------------------------------
# exhibited data taken from the paper (INPUT)
# --------------------------------------------------------------------------

K_MAIN = (22, 0, 5)                     # (k1, k2, k3) of the main example
LABEL_A = Fraction(1, 4)                # first colliding Farey label
LABEL_B = Fraction(2, 3)                # second colliding Farey label
CLAIM_A = (889, 3)                      # claimed (n_{1/4}, i_{1/4})
CLAIM_B = (889, 1)                      # claimed (n_{2/3}, i_{2/3})
SPINE_PRINTED = [15, 113, 889]          # printed left-spine mutation values
NODE_A_PRINTED = ((1, 1), (889, 3), (113, 2))
NODE_B_PRINTED = ((15, 3), (889, 1), (2, 2))
SOL_A_PRINTED = (1, 113, 889)           # ordered solution at the 1/4 node
SOL_B_PRINTED = (889, 2, 15)            # ordered solution at the 2/3 node
TOTAL_A_PRINTED = 3013710
TOTAL_B_PRINTED = 800100
WEIGHT_PRINTED = 5                      # 1+4 = 2+3
DEPTH_A_PRINTED = 3
DEPTH_B_PRINTED = 2
N23_PRINTED = "29 + 10*q"               # Lemma: n_{2/3} in MT(q,0,0,id)
M_RANGE_PRINTED = 11                    # paper's independent recomputation: 0 <= m <= 10

# ranges this program actually covers
BFS_DEPTH = 7                           # full binary tree depth for MT(22,0,5,id)
LEMMA_QS = (0, 1, 2, 3, 6, 7, 12, 35, 100)
LEMMA_DMAX = 30
Q_MAX = 40                              # q range for the n_{2/3} half of the lemma
CATALAN_DMAX = 60
PISANO_RMAX = 700
M_MAX_CONGRUENCE = 40
M_MAX_MUTATION = M_RANGE_PRINTED - 1    # 0..10, the paper's whole stated range
FIB_MAX = 60 * M_MAX_CONGRUENCE + 200


# --------------------------------------------------------------------------
# the generalized Markov equation
# --------------------------------------------------------------------------

def kof(k, idx):
    """k_idx for idx in {1,2,3}."""
    return k[idx - 1]


def gm_lhs(k, x, y, z):
    k1, k2, k3 = k
    return x * x + y * y + z * z + k1 * y * z + k2 * z * x + k3 * x * y


def gm_rhs(k, x, y, z):
    k1, k2, k3 = k
    return (3 + k1 + k2 + k3) * x * y * z


def gm_holds(k, x, y, z):
    return gm_lhs(k, x, y, z) == gm_rhs(k, x, y, z)


# --------------------------------------------------------------------------
# the generalized Markov tree
# --------------------------------------------------------------------------

class NotExact(Exception):
    pass


def gm_root(k):
    return ((1, 1), (kof(k, 2) + 2, 2), (1, 3))


def gm_left(k, v):
    (a, h), (b, i), (c, j) = v
    num = a * a + kof(k, j) * a * b + b * b
    q, r = divmod(num, c)
    if r != 0:
        raise NotExact("left mutation not exact")
    return ((a, h), (q, j), (b, i))


def gm_right(k, v):
    (a, h), (b, i), (c, j) = v
    num = b * b + kof(k, h) * b * c + c * c
    q, r = divmod(num, a)
    if r != 0:
        raise NotExact("right mutation not exact")
    return ((b, i), (q, h), (c, j))


def gm_step(k, v, letter):
    return gm_left(k, v) if letter == "L" else gm_right(k, v)


def gm_walk(k, path):
    v = gm_root(k)
    for letter in path:
        v = gm_step(k, v, letter)
    return v


def middle(v):
    return v[1]


def tags(v):
    return tuple(t for (_, t) in v)


def ordered_solution(v):
    """Place each node value at the coordinate position given by its tag."""
    slot = [None, None, None]
    for (val, tag) in v:
        if tag not in (1, 2, 3) or slot[tag - 1] is not None:
            return None
        slot[tag - 1] = val
    return tuple(slot)


# --------------------------------------------------------------------------
# the Farey (Stern-Brocot) tree
# --------------------------------------------------------------------------

FAREY_ROOT = ((0, 1), (1, 1), (1, 0))


def mediant(f, g):
    return (f[0] + g[0], f[1] + g[1])


def farey_left(tri):
    r, t, s = tri
    return (r, mediant(r, t), t)


def farey_right(tri):
    r, t, s = tri
    return (t, mediant(t, s), s)


def farey_step(tri, letter):
    return farey_left(tri) if letter == "L" else farey_right(tri)


def as_fraction(f):
    """The mediant of a Farey triple always has positive denominator."""
    return Fraction(f[0], f[1])


def farey_path(target, limit=2000):
    """Descend the Farey tree to the node whose middle label is target."""
    tri = FAREY_ROOT
    path = []
    for _ in range(limit):
        m = as_fraction(middle_of_triple(tri))
        if m == target:
            return path, tri
        if target < m:
            path.append("L")
            tri = farey_left(tri)
        else:
            path.append("R")
            tri = farey_right(tri)
    raise ValueError("Farey label not reached within the step limit")


def middle_of_triple(tri):
    return tri[1]


def label_at(path):
    tri = FAREY_ROOT
    for letter in path:
        tri = farey_step(tri, letter)
    return as_fraction(middle_of_triple(tri))


def node_at_label(k, target):
    """Return (path, gm node) for the node whose middle Farey label is target."""
    path, _tri = farey_path(target)
    return path, gm_walk(k, path)


# --------------------------------------------------------------------------
# Fibonacci numbers
# --------------------------------------------------------------------------

def fib_list(n):
    f = [0, 1]
    while len(f) <= n:
        f.append(f[-1] + f[-2])
    return f[:n + 1]


FIB = fib_list(FIB_MAX)


# --------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------

def check_equation_specialization():
    """The displayed specialization x^2+y^2+z^2+22yz+5xy = 30xyz is the
    equation of the paper's general form at (k1,k2,k3) = (22,0,5)."""
    pts = [(x, y, z)
           for x in (1, 2, 3, 7)
           for y in (1, 4, 113)
           for z in (2, 15, 889)]
    ok_l = all(gm_lhs(K_MAIN, x, y, z)
               == x * x + y * y + z * z + 22 * y * z + 5 * x * y
               for (x, y, z) in pts)
    ok_r = all(gm_rhs(K_MAIN, x, y, z) == 30 * x * y * z for (x, y, z) in pts)
    ck("displayed equation is the general form at (22,0,5), on %d integer points"
       % len(pts), ok_l and ok_r,
       "3+k1+k2+k3 = %d" % (3 + sum(K_MAIN)))


def check_parameters_distinct():
    k1, k2, k3 = K_MAIN
    ok = (len({k1, k2, k3}) == 3) and all(v >= 0 for v in K_MAIN)
    ck("parameter triple (22,0,5) is pairwise distinct and nonnegative", ok,
       "sorted = %s" % (tuple(sorted(K_MAIN)),))


def check_root():
    v = gm_root(K_MAIN)
    sol = ordered_solution(v)
    ok = (v == ((1, 1), (2, 2), (1, 3))
          and sol is not None
          and gm_holds(K_MAIN, *sol))
    ck("derived root of MT(22,0,5,id) is ((1,1),(2,2),(1,3)) and solves the equation",
       ok, "root = %s, ordered = %s" % (v, sol))


def check_farey_paths():
    pa, _ = farey_path(LABEL_A)
    pb, _ = farey_path(LABEL_B)
    mids = [label_at(["L"] * n) for n in (1, 2, 3)]
    ok = (pa == ["L", "L", "L"]
          and pb == ["L", "R"]
          and mids == [Fraction(1, 2), Fraction(1, 3), Fraction(1, 4)]
          and len(pa) == DEPTH_A_PRINTED
          and len(pb) == DEPTH_B_PRINTED)
    ck("Farey routing: 1/4 is three left steps, 2/3 is the right child of the 1/2 node",
       ok, "left-spine middles %s ; depths %d and %d"
       % ([str(m) for m in mids], len(pa), len(pb)))


def check_left_spine_values():
    v = gm_root(K_MAIN)
    got = []
    nodes = []
    for _ in range(3):
        v = gm_left(K_MAIN, v)
        got.append(middle(v)[0])
        nodes.append(v)
    ok = (got == SPINE_PRINTED
          and nodes[0] == ((1, 1), (15, 3), (2, 2))
          and nodes[1] == ((1, 1), (113, 2), (15, 3))
          and nodes[2] == NODE_A_PRINTED)
    ck("three derived left mutations reproduce the printed 15, 113, 889 and their nodes",
       ok, "derived = %s" % (got,))


def check_right_child_of_half():
    v = gm_walk(K_MAIN, ["L"])
    w = gm_right(K_MAIN, v)
    ok = (w == NODE_B_PRINTED and middle(w)[0] == 889)
    ck("derived right child of the 1/2 node is ((15,3),(889,1),(2,2))", ok,
       "derived = %s" % (w,))


def check_label_a():
    path, v = node_at_label(K_MAIN, LABEL_A)
    ok = (middle(v) == CLAIM_A and label_at(path) == LABEL_A)
    ck("derived (n_{1/4}, i_{1/4}) equals the claimed (889,3)", ok,
       "path %s -> %s" % ("".join(path), middle(v)))


def check_label_b():
    path, v = node_at_label(K_MAIN, LABEL_B)
    ok = (middle(v) == CLAIM_B and label_at(path) == LABEL_B)
    ck("derived (n_{2/3}, i_{2/3}) equals the claimed (889,1)", ok,
       "path %s -> %s" % ("".join(path), middle(v)))


def check_collision():
    _, va = node_at_label(K_MAIN, LABEL_A)
    _, vb = node_at_label(K_MAIN, LABEL_B)
    na, ia = middle(va)
    nb, ib = middle(vb)
    ok = (na == nb and LABEL_A != LABEL_B and ia != ib and va != vb)
    ck("t -> n_t is NOT injective in MT(22,0,5,id): equal values, distinct labels, distinct tags",
       ok, "n = %d at both 1/4 and 2/3, tags %d vs %d, nodes distinct: %s"
       % (na, ia, ib, va != vb))


def check_nodes_distinct_and_middle_slot():
    """The Status section's claim: the two nodes are distinct triples, yet 889
    stands in the second component -- the one recorded by the Farey label."""
    _, va = node_at_label(K_MAIN, LABEL_A)
    _, vb = node_at_label(K_MAIN, LABEL_B)
    ok = (va != vb
          and va[1][0] == 889 and vb[1][0] == 889
          and va[0][0] != 889 and va[2][0] != 889
          and vb[0][0] != 889 and vb[2][0] != 889)
    ck("the 1/4 and 2/3 nodes are distinct triples and 889 is the middle entry of each",
       ok, "%s / %s" % (va, vb))


def check_substitution_a():
    _, va = node_at_label(K_MAIN, LABEL_A)
    sol = ordered_solution(va)
    lhs = gm_lhs(K_MAIN, *sol)
    rhs = gm_rhs(K_MAIN, *sol)
    ok = (sol == SOL_A_PRINTED and lhs == rhs and lhs == TOTAL_A_PRINTED)
    ck("substitution at the 1/4 node: (1,113,889) solves the equation with total 3013710",
       ok, "ordered %s, lhs %d, rhs %d" % (sol, lhs, rhs))


def check_substitution_b():
    _, vb = node_at_label(K_MAIN, LABEL_B)
    sol = ordered_solution(vb)
    lhs = gm_lhs(K_MAIN, *sol)
    rhs = gm_rhs(K_MAIN, *sol)
    ok = (sol == SOL_B_PRINTED and lhs == rhs and lhs == TOTAL_B_PRINTED)
    ck("substitution at the 2/3 node: (889,2,15) solves the equation with total 800100",
       ok, "ordered %s, lhs %d, rhs %d" % (sol, lhs, rhs))


def check_ordering_not_free():
    """Of the six orderings of {2,15,889}, only (889,2,15) satisfies the
    equation; and 889 is the largest entry of both exhibited solutions."""
    vals = (2, 15, 889)
    good = []
    for x in vals:
        for y in vals:
            for z in vals:
                if sorted((x, y, z)) != sorted(vals):
                    continue
                if gm_holds(K_MAIN, x, y, z):
                    good.append((x, y, z))
    largest = (max(SOL_A_PRINTED) == 889 and max(SOL_B_PRINTED) == 889
               and SOL_A_PRINTED.index(889) == 2 and SOL_B_PRINTED.index(889) == 0)
    ok = (good == [SOL_B_PRINTED]) and largest
    ck("exactly one of the six orderings of {2,15,889} solves the equation, namely (889,2,15)",
       ok, "solutions found: %s ; 889 largest in both, at positions 3 and 1: %s"
       % (good, largest))


def check_tree_consistency():
    """Every mutation in the full binary tree to depth BFS_DEPTH is an exact
    division, every node's tags are a permutation of {1,2,3}, every node is an
    ordered solution of the equation, and every node's values are positive."""
    frontier = [(gm_root(K_MAIN), FAREY_ROOT, [])]
    nodes = 0
    bad = []
    labels = {}
    for depth in range(BFS_DEPTH + 1):
        nxt = []
        for (v, tri, path) in frontier:
            nodes += 1
            sol = ordered_solution(v)
            if sol is None or not gm_holds(K_MAIN, *sol):
                bad.append(("equation", path))
            if sorted(tags(v)) != [1, 2, 3]:
                bad.append(("tags", path))
            if any(val <= 0 for (val, _t) in v):
                bad.append(("positivity", path))
            lab = as_fraction(middle_of_triple(tri))
            prev = labels.get(lab)
            if prev is not None and prev != middle(v):
                bad.append(("label-collision", path))
            labels[lab] = middle(v)
            if depth < BFS_DEPTH:
                for letter in ("L", "R"):
                    try:
                        nxt.append((gm_step(K_MAIN, v, letter),
                                    farey_step(tri, letter),
                                    path + [letter]))
                    except NotExact:
                        bad.append(("inexact", path + [letter]))
        frontier = nxt
    ok = (not bad) and nodes == 2 ** (BFS_DEPTH + 1) - 1
    ck("all %d nodes of MT(22,0,5,id) to depth %d: exact mutations, tag permutations, exact solutions"
       % (nodes, BFS_DEPTH), ok,
       "distinct Farey labels %d, anomalies %d" % (len(labels), len(bad)))


def check_equal_weight_remark():
    """The least weight (numerator + denominator) admitting two distinct
    irreducible fractions in [0,1] is 5, realized exactly by {1/4, 2/3}."""
    from math import gcd
    first_w = None
    first_set = None
    for w in range(1, 40):
        fr = set()
        for den in range(1, w + 1):
            num = w - den
            if 0 <= num <= den and gcd(num, den) == 1:
                fr.add(Fraction(num, den))
        if len(fr) >= 2 and first_w is None:
            first_w = w
            first_set = fr
    pa, _ = farey_path(LABEL_A)
    pb, _ = farey_path(LABEL_B)
    wa = LABEL_A.numerator + LABEL_A.denominator
    wb = LABEL_B.numerator + LABEL_B.denominator
    ok = (first_w == WEIGHT_PRINTED
          and first_set == {LABEL_A, LABEL_B}
          and wa == wb == WEIGHT_PRINTED
          and len(pa) != len(pb))
    ck("least equal-weight pair of distinct irreducible fractions in [0,1] is {1/4,2/3} at weight 5",
       ok, "first duplicated weight %s realized by %s ; Farey depths %d vs %d"
       % (first_w, sorted(str(f) for f in (first_set or [])), len(pa), len(pb)))


def check_lemma_spine():
    """n_{1/d} = F_{2d+1} and i_{1/d} = 2 (d odd) / 3 (d even) in MT(q,0,0,id),
    obtained by mutating from the root rather than by quoting the lemma."""
    bad = []
    for q in LEMMA_QS:
        k = (q, 0, 0)
        for d in range(1, LEMMA_DMAX + 1):
            path, _ = farey_path(Fraction(1, d))
            if path != ["L"] * (d - 1):
                bad.append(("path", q, d))
            v = gm_walk(k, path)
            n, i = middle(v)
            if n != FIB[2 * d + 1]:
                bad.append(("value", q, d))
            if i != (2 if d % 2 == 1 else 3):
                bad.append(("tag", q, d))
    ok = not bad
    ck("spine lemma n_{1/d}=F_{2d+1}, i_{1/d} alternating, for %d values of q and d<=%d"
       % (len(LEMMA_QS), LEMMA_DMAX), ok,
       "q in %s, anomalies %d" % (LEMMA_QS, len(bad)))


def check_lemma_n23():
    """n_{2/3} = 29 + 10q and i_{2/3} = 1 in MT(q,0,0,id)."""
    bad = []
    path, _ = farey_path(LABEL_B)
    for q in range(0, Q_MAX + 1):
        v = gm_walk((q, 0, 0), path)
        n, i = middle(v)
        if n != 29 + 10 * q or i != 1:
            bad.append(q)
        if v != ((5, 3), (29 + 10 * q, 1), (2, 2)):
            bad.append(("node", q))
    ok = not bad
    ck("spine lemma %s with i_{2/3}=1 in MT(q,0,0,id) for 0<=q<=%d" % (N23_PRINTED, Q_MAX),
       ok, "anomalies %d" % len(bad))


def check_catalan_instance():
    """F_{2d-1} F_{2d+3} = F_{2d+1}^2 + 1, the identity the induction uses."""
    bad = [d for d in range(1, CATALAN_DMAX + 1)
           if FIB[2 * d - 1] * FIB[2 * d + 3] != FIB[2 * d + 1] ** 2 + 1]
    ok = not bad
    ck("Fibonacci identity F_{2d-1}F_{2d+3} = F_{2d+1}^2 + 1 for 1<=d<=%d" % CATALAN_DMAX,
       ok, "failures %d" % len(bad))


def check_pisano():
    """(F_60, F_61) = (0,1) mod 10, hence F_{r+60} = F_r mod 10."""
    pair_ok = (FIB[60] % 10 == 0 and FIB[61] % 10 == 1)
    bad = [r for r in range(0, PISANO_RMAX + 1) if FIB[r + 60] % 10 != FIB[r] % 10]
    ok = pair_ok and not bad
    ck("period-60 congruence: (F_60,F_61)=(0,1) mod 10 and F_{r+60}=F_r mod 10 for r<=%d"
       % PISANO_RMAX, ok,
       "F_60 mod 10 = %d, F_61 mod 10 = %d, failures %d"
       % (FIB[60] % 10, FIB[61] % 10, len(bad)))


def check_family_integrality():
    """q_m = (F_{60m+11} - 29)/10 is a nonnegative integer, q_m >= 6, strictly
    increasing; d_m = 30m+5 is odd and 2 d_m + 1 = 60m + 11."""
    bad = []
    prev = None
    for m in range(0, M_MAX_CONGRUENCE + 1):
        f = FIB[60 * m + 11]
        d = 30 * m + 5
        if f % 10 != 9:
            bad.append(("mod10", m))
        num, rem = divmod(f - 29, 10)
        if rem != 0 or num < 0 or num < 6:
            bad.append(("q", m))
        if prev is not None and not num > prev:
            bad.append(("monotone", m))
        prev = num
        if d % 2 != 1 or 2 * d + 1 != 60 * m + 11 or d < 5:
            bad.append(("d", m))
    ok = not bad
    ck("q_m integral, >=6 and strictly increasing, d_m odd with 2d_m+1=60m+11, for 0<=m<=%d"
       % M_MAX_CONGRUENCE, ok,
       "F_11 = %d, q_0 = %d, q_%d has %d digits, anomalies %d"
       % (FIB[11], (FIB[11] - 29) // 10, M_MAX_CONGRUENCE,
          len(str((FIB[60 * M_MAX_CONGRUENCE + 11] - 29) // 10)), len(bad)))


def check_family_by_mutation():
    """For each m in the paper's range, mutate MT(q_m,0,0,id) from its root and
    confirm the collision n_{1/d_m} = n_{2/3} = F_{60m+11} with tags 2 and 1."""
    bad = []
    path_b, _ = farey_path(LABEL_B)
    for m in range(0, M_MAX_MUTATION + 1):
        f = FIB[60 * m + 11]
        q = (f - 29) // 10
        d = 30 * m + 5
        k = (q, 0, 0)
        label_dm = Fraction(1, d)
        path_a, _ = farey_path(label_dm)
        va = gm_walk(k, path_a)
        vb = gm_walk(k, path_b)
        na, ia = middle(va)
        nb, ib = middle(vb)
        if (na, ia) != (f, 2):
            bad.append(("1/d_m", m, na == f, ia))
        if (nb, ib) != (f, 1):
            bad.append(("2/3", m, nb == f, ib))
        if na != nb or label_dm == LABEL_B or ia == ib or va == vb:
            bad.append(("collision", m))
        for v in (va, vb):
            sol = ordered_solution(v)
            if sol is None or not gm_holds(k, *sol):
                bad.append(("equation", m))
    ok = not bad
    ck("family verified by mutation from the root for every m with 0<=m<=%d"
       % M_MAX_MUTATION, ok,
       "m=0 gives q=%d, d=5, n=%d ; m=%d gives a %d-digit label ; anomalies %d"
       % ((FIB[11] - 29) // 10, FIB[11], M_MAX_MUTATION,
          len(str(FIB[60 * M_MAX_MUTATION + 11])), len(bad)))


def check_family_first_member():
    """The stated first member: n_{1/5} = n_{2/3} = 89 in MT(6,0,0,id)."""
    k = (6, 0, 0)
    pa, va = node_at_label(k, Fraction(1, 5))
    pb, vb = node_at_label(k, LABEL_B)
    q0 = (FIB[11] - 29) // 10
    ok = (q0 == 6
          and middle(va) == (89, 2)
          and middle(vb) == (89, 1)
          and FIB[11] == 89
          and va != vb)
    ck("first family member: derived n_{1/5}=(89,2) and n_{2/3}=(89,1) in MT(6,0,0,id)",
       ok, "q_0 = %d, nodes %s and %s" % (q0, va, vb))


def check_repeated_values_have_distinct_tags():
    """Scope statement: the collisions exhibited are collisions of t -> n_t and
    leave t -> (n_t,i_t) untouched. Derived form: enumerate the searched region,
    collect every middle value that occurs at more than one Farey label, and
    confirm that each such value carries pairwise distinct position tags -- and
    that 889 is one of them."""
    occ = {}
    frontier = [(gm_root(K_MAIN), FAREY_ROOT)]
    for depth in range(BFS_DEPTH + 1):
        nxt = []
        for (v, tri) in frontier:
            n, i = middle(v)
            occ.setdefault(n, []).append((as_fraction(middle_of_triple(tri)), i))
            if depth < BFS_DEPTH:
                for letter in ("L", "R"):
                    nxt.append((gm_step(K_MAIN, v, letter),
                                farey_step(tri, letter)))
        frontier = nxt
    repeated = {n: lst for (n, lst) in occ.items() if len(lst) > 1}
    tag_clash = [n for (n, lst) in repeated.items()
                 if len(set(i for (_f, i) in lst)) != len(lst)]
    ok = (889 in repeated) and (not tag_clash)
    ck("every repeated middle value to depth %d carries pairwise distinct tags, 889 included"
       % BFS_DEPTH, ok,
       "distinct values %d, repeated values %s, tag clashes %d"
       % (len(occ), sorted(repeated), len(tag_clash)))


def check_value_collisions_in_tree():
    """Independently of the Farey routing, the value 889 occurs as a middle
    entry at exactly two distinct labels in the searched region of the tree,
    and those labels are 1/4 and 2/3."""
    occurrences = []
    frontier = [(gm_root(K_MAIN), FAREY_ROOT)]
    for depth in range(BFS_DEPTH + 1):
        nxt = []
        for (v, tri) in frontier:
            n, i = middle(v)
            if n == 889:
                occurrences.append((as_fraction(middle_of_triple(tri)), i))
            if depth < BFS_DEPTH:
                for letter in ("L", "R"):
                    nxt.append((gm_step(K_MAIN, v, letter),
                                farey_step(tri, letter)))
        frontier = nxt
    labels = sorted((f for (f, _i) in occurrences))
    ok = (len(occurrences) == 2
          and labels == sorted([LABEL_A, LABEL_B])
          and sorted(i for (_f, i) in occurrences) == [1, 3])
    ck("independent scan to depth %d finds the value 889 at exactly the labels 1/4 and 2/3"
       % BFS_DEPTH, ok,
       "occurrences %s" % ([(str(f), i) for (f, i) in occurrences],))


def summarize():
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
    else:
        print("VERDICT: ALL %d CHECKS PASS" % n)
    return 1 if bad else 0


NOT_RE_RUN = (
    "NOT RE-RUN HERE: "
    "(1) the infinitude of the family -- q_m integrality and monotonicity were "
    "re-derived only for 0 <= m <= %d, and the full mutation of MT(q_m,0,0,id) "
    "from its root only for 0 <= m <= %d (the paper's own independently "
    "recomputed range), so the statement for all m >= 0, which rests on the "
    "induction and the mod-10 period argument, is not machine-checked; "
    "(2) the spine lemma as a theorem -- n_{1/d}=F_{2d+1} and the alternating "
    "tags were mutated out only for q in %s and d <= %d, and n_{2/3}=29+10q only "
    "for 0 <= q <= %d, with the Fibonacci identity instance checked for d <= %d "
    "and the mod-10 period for r <= %d; the general proofs by induction and "
    "Catalan's identity in general are not verified; "
    "(3) tree-wide facts for MT(22,0,5,id) -- exactness of mutations, tag "
    "permutations, the equation at every node, and the scans for repeated values "
    "or repeated (n_t,i_t) pairs cover only the full binary tree to depth %d, not "
    "the whole infinite tree, and therefore establish nothing about collisions "
    "deeper than that; and the Status section's algebraic remark that the "
    "equation is symmetric in x,y,z only when k1=k2=k3 is exercised here only at "
    "the single triple (22,0,5), through the scan of the six orderings of "
    "{2,15,889}, and is not verified for general k; "
    "(4) everything external to the paper's own arithmetic -- that the quoted "
    "conjecture is Conjecture 7.6 of the cited reference and is stated there as "
    "quoted, that the cited reference records that no counterexample was known, "
    "the definitions of the generalized Markov tree, of the root, of the mutation "
    "rules and of the synchronization with the Farey tree (all taken here from the "
    "paper's own recollection of them and not from the sources), the relation to "
    "Question 7.3 and to the separate Chen-Jia Conjecture 8.2, the posting dates "
    "of the cited preprints, the literature search, and the novelty/priority claim "
    "that no earlier refutation or earlier appearance of the collision at 889 "
    "exists, and the paper's account of its own provenance -- that the printed "
    "data were recomputed once independently, on an earlier occasion, with every "
    "check agreeing -- which this program can neither confirm nor deny, though "
    "the recomputation the paper describes there (the six checks of Section 2, "
    "and for each m with 0 <= m <= 10 the integrality of q_m together with the "
    "two labels obtained by mutating MT(q_m,0,0,id) from its root) is carried out "
    "among the checks above, so what is unverifiable is the history of that "
    "earlier run and not the facts it is said to have confirmed; "
    "(5) claims the paper explicitly does not make and this program does not "
    "supply -- global minimality of the counterexample, any statement about the "
    "classical case (0,0,0) or the symmetric locus k1=k2=k3, trees with sigma "
    "other than the identity, injectivity of t -> (n_t,i_t) as a theorem, and "
    "infinitely many collisions inside one fixed tree."
) % (M_MAX_CONGRUENCE, M_MAX_MUTATION, LEMMA_QS, LEMMA_DMAX, Q_MAX,
     CATALAN_DMAX, PISANO_RMAX, BFS_DEPTH)


def run(fn):
    """A check that raises (e.g. an inexact mutation, an unreachable Farey
    label) counts as a failed check, not as a crash."""
    before = len(CHECKS)
    try:
        fn()
    except Exception as exc:              # noqa: BLE001 - deliberate
        if len(CHECKS) == before:
            ck(fn.__name__.replace("check_", "").replace("_", " ")
               + " (raised an exception)", False, repr(exc))


def main():
    for fn in (
        check_equation_specialization,
        check_parameters_distinct,
        check_root,
        check_farey_paths,
        check_left_spine_values,
        check_right_child_of_half,
        check_label_a,
        check_label_b,
        check_collision,
        check_nodes_distinct_and_middle_slot,
        check_substitution_a,
        check_substitution_b,
        check_ordering_not_free,
        check_tree_consistency,
        check_repeated_values_have_distinct_tags,
        check_value_collisions_in_tree,
        check_equal_weight_remark,
        check_lemma_spine,
        check_lemma_n23,
        check_catalan_instance,
        check_pisano,
        check_family_integrality,
        check_family_by_mutation,
        check_family_first_member,
    ):
        run(fn)
    rc = summarize()
    print(NOT_RE_RUN)
    return rc


if __name__ == "__main__":
    sys.exit(main())

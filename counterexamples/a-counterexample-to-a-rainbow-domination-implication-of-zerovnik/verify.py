#!/usr/bin/env python3
"""Verification program for a note on k-rainbow domination of the
Tutte-Coxeter (Levi) graph.

Decisive claim under test:
    (gamma_r3, gamma_r4, gamma_r5) of the Tutte-Coxeter graph = (16, 21, 25),
    with |V| = 30, so gamma_r5 <= 5|V|/6 while gamma_r4 > 2|V|/3.

Standard library only; exact integer / Fraction arithmetic throughout.
The graph is DERIVED from its incidence definition, not typed in: the only
inputs taken from the paper are the duad-syntheme description of the graph
and the three exhibited labellings f_3, f_4, f_5.

Upper bounds: each exhibited labelling is checked to be a k-RDF and its
weight is summed, giving gamma_r3 <= 16, gamma_r4 <= 21, gamma_r5 <= 25.

Lower bounds, established here without any solver:
  (i)   an exhaustive integer scan of the paper's counting inequalities shows
        that a k-RDF of weight at most 5k forces weight exactly 5k, exactly
        15 nonempty vertices and no edge among them (so gamma_rk >= 5k);
  (ii)  an exhaustive branch-and-bound shows the only independent 15-sets are
        the two bipartition classes, so the nonempty vertices form one class;
  (iii) an exhaustive scan of all 2^15 subsets of a class shows every subset
        of one class dominating the other has size at least 5, and that those
        of size exactly 5 are precisely six "perfect classes";
  (iv)  hence at weight 5k the k colour supports are k perfect classes whose
        union is a whole class; an exhaustive scan of all multisets of perfect
        classes shows no 3 and no 4 of them cover a class, which rules out
        weights 15 and 20 and gives gamma_r3 >= 16, gamma_r4 >= 21.
"""

import sys
import itertools
from fractions import Fraction

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         ("  --  " + detail) if detail else ""))
    return bool(ok)


def finish(disclosure):
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
    else:
        print("VERDICT: ALL %d CHECKS PASS" % n)
    print(disclosure)
    sys.exit(1 if bad else 0)


# ---------------------------------------------------------------------------
# 1.  The duad-syntheme model, derived from Omega = {1,...,6}
# ---------------------------------------------------------------------------

OMEGA = tuple(range(1, 7))

DUADS = sorted((frozenset(p) for p in itertools.combinations(OMEGA, 2)),
               key=lambda d: sorted(d))


def _matchings(points):
    """All partitions of the tuple `points` into 2-element blocks."""
    if not points:
        yield frozenset()
        return
    a = points[0]
    for j in range(1, len(points)):
        b = points[j]
        rest = points[1:j] + points[j + 1:]
        for tail in _matchings(rest):
            yield tail | {frozenset((a, b))}


SYNTHEMES = sorted(_matchings(OMEGA), key=lambda s: sorted(sorted(d) for d in s))

# Vertices: ('D', duad) and ('S', syntheme).  Adjacency: d in s.
VD = [('D', d) for d in DUADS]
VS = [('S', s) for s in SYNTHEMES]
VERTS = VD + VS
IDX = {v: i for i, v in enumerate(VERTS)}

ADJ = {v: set() for v in VERTS}
for s in SYNTHEMES:
    for d in s:
        ADJ[('D', d)].add(('S', s))
        ADJ[('S', s)].add(('D', d))

EDGES = sorted(((IDX[u], IDX[v]) for u in VERTS for v in ADJ[u] if IDX[u] < IDX[v]))


def ck_model():
    ok_d = (len(DUADS) == 15 and all(len(d) == 2 for d in DUADS)
            and len(set(DUADS)) == 15)
    ok_s = (len(SYNTHEMES) == 15 and len(set(SYNTHEMES)) == 15
            and all(len(s) == 3 and set().union(*s) == set(OMEGA)
                    and sum(len(d) for d in s) == 6 for s in SYNTHEMES))
    ck("model: 15 duads and 15 synthemes, each syntheme a partition of Omega",
       ok_d and ok_s, "|D|=%d |S|=%d" % (len(DUADS), len(SYNTHEMES)))

    per_duad = sorted(sum(1 for s in SYNTHEMES if d in s) for d in DUADS)
    ck("model: every duad lies in exactly three synthemes",
       per_duad == [3] * 15, "multiset of duad-degrees = %s" % (set(per_duad),))


def ck_graph():
    degs = sorted(len(ADJ[v]) for v in VERTS)
    ck("graph: order 30, size 45, 3-regular",
       len(VERTS) == 30 and len(EDGES) == 45 and degs == [3] * 30,
       "n=%d m=%d degrees=%s" % (len(VERTS), len(EDGES), set(degs)))

    # connectivity and 2-colouring by BFS from one vertex
    start = VERTS[0]
    colour = {start: 0}
    order = [start]
    qi = 0
    while qi < len(order):
        u = order[qi]
        qi += 1
        for w in ADJ[u]:
            if w not in colour:
                colour[w] = 1 - colour[u]
                order.append(w)
    connected = len(colour) == 30
    proper = all(colour[a] != colour[b] for a, b in
                 ((VERTS[i], VERTS[j]) for i, j in EDGES))
    classes = (frozenset(v for v in VERTS if colour.get(v) == 0),
               frozenset(v for v in VERTS if colour.get(v) == 1))
    matches = set(classes) == {frozenset(VD), frozenset(VS)}
    ck("graph: connected and bipartite, and the unique bipartition is "
       "{duads, synthemes}",
       connected and proper and matches,
       "connected=%s proper2colouring=%s classes=D/S:%s"
       % (connected, proper, matches))


def girth():
    best = None
    for (i, j) in EDGES:
        u, v = VERTS[i], VERTS[j]
        dist = {u: 0}
        q = [u]
        qi = 0
        while qi < len(q):
            x = q[qi]
            qi += 1
            for y in ADJ[x]:
                if x == u and y == v:
                    continue
                if y == v and x == u:
                    continue
                if y not in dist:
                    dist[y] = dist[x] + 1
                    q.append(y)
        # shortest u-v path avoiding the edge uv
        if v in dist and dist[v] > 1:
            cand = dist[v] + 1
            if best is None or cand < best:
                best = cand
    return best


def ck_girth():
    g = girth()
    ck("graph: girth is 8 (the Tutte 8-cage)", g == 8, "girth=%s" % (g,))


# ---------------------------------------------------------------------------
# 2.  The three exhibited functions, transcribed from the paper's definitions
#     and applied to the DERIVED duads.  Synthemes always get the empty label.
# ---------------------------------------------------------------------------

D56 = frozenset((5, 6))
D16 = frozenset((1, 6))


def f5_label(d):
    return frozenset(x for x in d if x <= 5)


def f4_label(d):
    if d == D56:
        return frozenset((1,))
    return frozenset(x for x in d if x <= 4)


def f3_label(d):
    if d == D16:
        return frozenset((1, 2))
    if 6 in d and 1 not in d:
        return frozenset((1,))
    if 1 in d and 6 not in d:
        return frozenset((2,))
    if not (d & {1, 6}):
        return frozenset((3,))
    raise AssertionError("f3 is undefined on a duad: the transcription is "
                         "not a total function")


def build(labeller):
    f = {}
    for d in DUADS:
        f[('D', d)] = labeller(d)
    for s in SYNTHEMES:
        f[('S', s)] = frozenset()
    return f


def is_krdf(f, k):
    """Return (ok, list of violations).  f maps every vertex to a subset of [k]."""
    bad = []
    colours = set(range(1, k + 1))
    for v in VERTS:
        if not f[v] <= colours:
            bad.append(('label out of range', v, sorted(f[v])))
    for v in VERTS:
        if f[v]:
            continue
        seen = set()
        for u in ADJ[v]:
            seen |= f[u]
        if seen != colours:
            bad.append(('undominated', v, sorted(colours - seen)))
    return (not bad), bad


def weight(f):
    return sum(len(f[v]) for v in VERTS)


def ck_upper(k, labeller, claimed_weight, name):
    f = build(labeller)
    parsed = (len(f) == 30
              and all(f[('D', d)] for d in DUADS)
              and all(not f[('S', s)] for s in SYNTHEMES))
    ck("f_%d parses: 30 labels, every duad nonempty, every syntheme empty"
       % k, parsed,
       "labels=%d distinct duad labels=%d"
       % (len(f), len(set(f[('D', d)] for d in DUADS))))

    ok, bad = is_krdf(f, k)
    ck("f_%d (%s) really is a %d-rainbow dominating function: every empty "
       "vertex sees all %d colours" % (k, name, k, k), ok,
       "violations=%d%s" % (len(bad), ("  first=%s" % (bad[0],)) if bad else ""))

    w = weight(f)
    ck("w(f_%d) = %d as claimed (derived by summing the label sizes)"
       % (k, claimed_weight), w == claimed_weight, "derived weight=%d" % w)
    return f


def supports(f, k):
    return [frozenset(v for v in VERTS if i in f[v]) for i in range(1, k + 1)]


def ck_support_sizes():
    f3 = build(f3_label)
    sizes = sorted(len(q) for q in supports(f3, 3))
    ck("f_3 colour supports have sizes 5, 5, 6 (the paper's 5, 5, C(4,2))",
       sizes == [5, 5, 6], "derived sizes=%s, sum=%d" % (sizes, sum(sizes)))

    f5 = build(f5_label)
    s5 = sorted(len(q) for q in supports(f5, 5))
    ck("f_5 colour supports all have size 5 (each colour on five duads)",
       s5 == [5] * 5, "derived sizes=%s" % (s5,))


def ck_exception_needed():
    lone = [d for d in DUADS if not (d & {1, 2, 3, 4})]
    unique = (lone == [D56])
    naive = build(lambda d: frozenset(x for x in d if x <= 4))
    ok_naive, bad_naive = is_krdf(naive, 4)
    ck("the exceptional value of f_4 is necessary: {5,6} is the only duad "
       "disjoint from [4], and d -> d cap [4] is NOT a 4-RDF",
       unique and (not ok_naive),
       "duads disjoint from [4]=%s ; violations of the naive map=%d"
       % ([sorted(d) for d in lone], len(bad_naive)))


# ---------------------------------------------------------------------------
# 3.  Perfect classes, by exhaustive enumeration over one bipartition class
# ---------------------------------------------------------------------------

def enumerate_classes(side, other):
    """Exhaustive scan of all 2^15 subsets of `side`.

    Returns (perfect, dominating) where `perfect` are the subsets Q such that
    every vertex of `other` has exactly one neighbour in Q, and `dominating`
    are those where every vertex of `other` has at least one neighbour in Q.
    """
    n = len(side)
    masks = []
    for v in other:
        m = 0
        for j, u in enumerate(side):
            if u in ADJ[v]:
                m |= 1 << j
        masks.append(m)
    perfect, dominating = [], []
    for bits in range(1 << n):
        good_all = True
        exact = True
        for m in masks:
            c = bin(bits & m).count('1')
            if c == 0:
                good_all = False
                break
            if c != 1:
                exact = False
        if good_all:
            q = frozenset(side[j] for j in range(n) if bits >> j & 1)
            dominating.append(q)
            if exact:
                perfect.append(q)
    return perfect, dominating


PERF_D, DOM_D = enumerate_classes(VD, VS)
PERF_S, DOM_S = enumerate_classes(VS, VD)


def ck_perfect_classes():
    stars = [frozenset(('D', d) for d in DUADS if i in d) for i in OMEGA]
    ck("perfect classes on the duad side: exactly the six stars D_i, each of "
       "size 5 (exhaustive over all 2^15 subsets)",
       len(PERF_D) == 6 and set(PERF_D) == set(stars)
       and all(len(q) == 5 for q in PERF_D),
       "found %d perfect classes, sizes %s"
       % (len(PERF_D), sorted(len(q) for q in PERF_D)))

    ck("perfect classes on the syntheme side: exactly six, each of size 5 "
       "(verified directly, so the paper's self-duality is not needed here)",
       len(PERF_S) == 6 and all(len(q) == 5 for q in PERF_S),
       "found %d perfect classes, sizes %s"
       % (len(PERF_S), sorted(len(q) for q in PERF_S)))

    md = min(len(q) for q in DOM_D)
    ms = min(len(q) for q in DOM_S)
    minimal_d = [q for q in DOM_D if len(q) == md]
    minimal_s = [q for q in DOM_S if len(q) == ms]
    ck("a subset of one class dominating the other has size >= 5, and every "
       "such subset of size 5 is a perfect class",
       md == 5 and ms == 5
       and set(minimal_d) == set(PERF_D) and set(minimal_s) == set(PERF_S),
       "min dominating sizes = %d (duads) and %d (synthemes); minimal "
       "families of sizes %d and %d" % (md, ms, len(minimal_d), len(minimal_s)))


def ck_cover_lemma():
    fails = []
    for side, perf, name in ((VD, PERF_D, 'duads'), (VS, PERF_S, 'synthemes')):
        full = frozenset(side)
        for r in range(0, 5):
            for combo in itertools.combinations(perf, r):
                u = frozenset().union(*combo) if combo else frozenset()
                if u == full:
                    fails.append((name, r))
        got5 = any(frozenset().union(*combo) == full
                   for combo in itertools.combinations(perf, 5))
        if not got5:
            fails.append((name, 'no five cover'))
    ck("cover lemma: no four or fewer perfect classes cover a bipartition "
       "class, while some five do (both sides, exhaustive)",
       not fails, "counterexamples to the lemma: %s" % (fails or 'none',))


# ---------------------------------------------------------------------------
# 4.  The counting bound, verified as an exhaustive integer feasibility scan
# ---------------------------------------------------------------------------
#
# Let f be a k-RDF of weight w on our graph (order 30, 3-regular, 45 edges),
# let C be the set of nonempty vertices, c = |C| and e the number of edges
# inside C.  The paper's double counting of the triples (x, y, i) with
# x outside C, xy an edge and i in f(y) yields
#
#     (A)  k(30 - c) <= tau            (rainbow domination at each empty x)
#     (B)  tau       <= 3w - 2e        (each y contributes |f(y)|(3-deg_C y))
#     (C)  3c - 2e   <= 3(30 - c)      (the C-to-complement cut is at most
#                                       3 times the number of empty vertices)
#     (D)  c <= w                      (every vertex of C has |f(v)| >= 1)
#
# Everything below scans ALL integer triples (w, c, e) in range that satisfy
# (A)+(B), (C) and (D), so the conclusion is a finite verified fact about the
# inequality system rather than a re-typed algebraic manipulation.

N_ORDER = len(VERTS)
N_EDGES = len(EDGES)
DEGREE = 3


def feasible_triples(k, wmax):
    out = []
    for w in range(0, wmax + 1):
        for c in range(0, N_ORDER + 1):
            if c > w:
                continue
            for e in range(0, N_EDGES + 1):
                if k * (N_ORDER - c) > DEGREE * w - 2 * e:
                    continue
                if DEGREE * c - 2 * e > DEGREE * (N_ORDER - c):
                    continue
                out.append((w, c, e))
    return out


def ck_counting_bound():
    for k in (3, 4, 5):
        wmax = 5 * k
        tri = feasible_triples(k, wmax)
        ck("counting bound for k=%d: among all integer (w,c,e) with w <= %d "
           "the inequality system admits only (w,c,e)=(%d,15,0); hence "
           "gamma_r%d >= %d and any function of weight %d has C independent "
           "of size 15" % (k, wmax, wmax, k, wmax, wmax),
           tri == [(wmax, 15, 0)],
           "feasible triples=%s" % (tri if len(tri) <= 4 else
                                    "%d triples, e.g. %s" % (len(tri), tri[:4])))


def independent_sets_at_least(target):
    nb = [0] * N_ORDER
    for i, j in EDGES:
        nb[i] |= 1 << j
        nb[j] |= 1 << i
    found = []
    best = [0]
    full = (1 << N_ORDER) - 1

    def rec(i, mask, count, blocked):
        if i == N_ORDER:
            if count >= target:
                found.append((count, mask))
            if count > best[0]:
                best[0] = count
            return
        # upper bound: only unblocked vertices from i onwards can still be used
        avail = bin(((~blocked) & full) >> i).count('1')
        if count + avail < target:
            return
        if not (blocked >> i & 1):
            rec(i + 1, mask | (1 << i), count + 1, blocked | nb[i])
        rec(i + 1, mask, count, blocked)

    rec(0, 0, 0, 0)
    return found, best[0]


_IND_CACHE = {}


def independent_15():
    if 'r' not in _IND_CACHE:
        _IND_CACHE['r'] = independent_sets_at_least(15)
    return _IND_CACHE['r']


def mask_of(vs):
    m = 0
    for v in vs:
        m |= 1 << IDX[v]
    return m


def bipartition_classes_are_the_only_independent_15_sets():
    found, _ = independent_15()
    at15 = set(m for (sz, m) in found if sz == 15)
    bigger = [sz for (sz, m) in found if sz > 15]
    return (not bigger) and at15 == {mask_of(VD), mask_of(VS)}


def minimum_dominating_sets_are_perfect():
    ok = True
    for dom, perf in ((DOM_D, PERF_D), (DOM_S, PERF_S)):
        m = min(len(q) for q in dom)
        ok = ok and m == 5 and set(q for q in dom if len(q) == m) == set(perf)
    return ok


def ck_independent_sets():
    found, best = independent_15()
    at15 = [m for (sz, m) in found if sz == 15]
    bigger = [sz for (sz, m) in found if sz > 15]
    want = {mask_of(VD), mask_of(VS)}
    ck("independence: the graph has no independent set of 16 vertices, and "
       "its independent 15-sets are exactly the duad class and the syntheme "
       "class (exhaustive branch-and-bound over all 30 vertices)",
       not bigger and len(at15) == 2 and set(at15) == want,
       "max independent size=%d, number of independent 15-sets=%d"
       % (best, len(at15)))


def ck_extremal_impossible():
    """No k-RDF of weight 5k for k = 3, 4: the k colour supports would have to
    be perfect classes covering a whole bipartition class."""
    bad = []
    for k in (3, 4):
        for side, perf, name in ((VD, PERF_D, 'duads'), (VS, PERF_S, 'synthemes')):
            full = frozenset(side)
            for combo in itertools.combinations_with_replacement(perf, k):
                if frozenset().union(*combo) == full:
                    bad.append((k, name))
    ck("weights 15 and 20 are impossible: no multiset of 3 (resp. 4) perfect "
       "classes has union a whole bipartition class, on either side",
       not bad, "coverings found: %s" % (bad or 'none',))

    # the k = 5 case is consistent, and f_5 realises it
    f5 = build(f5_label)
    sup5 = supports(f5, 5)
    ok5 = (all(q in set(PERF_D) for q in sup5)
           and frozenset().union(*sup5) == frozenset(VD)
           and len(set(sup5)) == 5)
    ck("for k=5 the extremal configuration exists and f_5 realises it: its "
       "five supports are five distinct perfect classes covering the duads",
       ok5, "distinct supports=%d, all perfect=%s"
       % (len(set(sup5)), all(q in set(PERF_D) for q in sup5)))


def extremal_covering_exists(k):
    for side, perf in ((VD, PERF_D), (VS, PERF_S)):
        full = frozenset(side)
        for combo in itertools.combinations_with_replacement(perf, k):
            if frozenset().union(*combo) == full:
                return True
    return False


# ---------------------------------------------------------------------------
# 5.  The 0/1 program quoted in the paper, and the elementary condition count
# ---------------------------------------------------------------------------

def ck_program_model():
    """The paper states that the feasible points of

           sum_{u in N(v)} x_{u,i} + sum_{j in [k]} x_{v,j} >= 1   (v, i)

    are exactly the k-RDFs under f(v) = {i : x_{v,i} = 1}.  Both sides of
    that equivalence depend only on the labels on a closed neighbourhood, so
    for a 3-regular graph the claim is settled by scanning every tuple
    (f(v), f(u1), f(u2), f(u3)) of subsets of [k]."""
    mismatches = 0
    tuples_scanned = 0
    for k in (3, 4, 5):
        colours = list(range(1, k + 1))
        subsets = [frozenset(c) for r in range(k + 1)
                   for c in itertools.combinations(colours, r)]
        for nbrs in itertools.product(subsets, repeat=3):
            union = frozenset().union(*nbrs)
            # the left sums of the quoted program, per colour
            counts = [sum(1 for t in nbrs if i in t) for i in colours]
            for lab in subsets:
                tuples_scanned += 1
                rdf_ok = bool(lab) or (union == frozenset(colours))
                lp_ok = all(counts[i - 1] + len(lab) >= 1 for i in colours)
                if rdf_ok != lp_ok:
                    mismatches += 1
    ck("the 0/1 program quoted in the paper has exactly the k-RDFs as its "
       "feasible points (exhaustive over every closed-neighbourhood label "
       "tuple, k=3,4,5, degree 3)",
       mismatches == 0,
       "label tuples scanned=%d, disagreeing with the rainbow condition=%d"
       % (tuples_scanned, mismatches))

    obj = []
    for k, lab in ((3, f3_label), (4, f4_label), (5, f5_label)):
        f = build(lab)
        x = {(v, i): (1 if i in f[v] else 0)
             for v in VERTS for i in range(1, k + 1)}
        feas = all(sum(x[(u, i)] for u in ADJ[v]) + sum(x[(v, j)] for j in
                                                        range(1, k + 1)) >= 1
                   for v in VERTS for i in range(1, k + 1))
        obj.append((k, feas, sum(x.values())))
    ck("each exhibited labelling gives a feasible 0/1 point whose objective "
       "equals its weight (objectives 16, 21, 25 for k=3,4,5)",
       all(f for (_, f, _) in obj)
       and [o for (_, _, o) in obj] == [16, 21, 25],
       "(k, feasible, objective) = %s" % (obj,))


def ck_ninetythree():
    """The paper reports 93 elementary conditions: 45 syntheme unions,
    45 nonemptiness conditions and 3 weight sums."""
    conditions = []
    for k, lab in ((5, f5_label), (4, f4_label), (3, f3_label)):
        f = build(lab)
        for s in SYNTHEMES:
            u = frozenset().union(*(f[('D', d)] for d in s))
            conditions.append(('union', k, u == frozenset(range(1, k + 1))))
        for d in DUADS:
            conditions.append(('nonempty', k, bool(f[('D', d)])))
        conditions.append(('weight', k, weight(f) == {3: 16, 4: 21, 5: 25}[k]))
    unions = sum(1 for c in conditions if c[0] == 'union')
    nonempty = sum(1 for c in conditions if c[0] == 'nonempty')
    weights = sum(1 for c in conditions if c[0] == 'weight')
    failed = [c for c in conditions if not c[2]]
    ck("the paper's 93 elementary conditions (45 syntheme unions, 45 "
       "nonemptiness conditions, 3 weight sums) all hold",
       len(conditions) == 93 and unions == 45 and nonempty == 45
       and weights == 3 and not failed,
       "conditions=%d (%d/%d/%d), failing=%d"
       % (len(conditions), unions, nonempty, weights, len(failed)))


# ---------------------------------------------------------------------------
# 6.  Assembling the three exact values and the refutation
# ---------------------------------------------------------------------------

DERIVED = {}


def ck_values():
    # the three structural facts the weight-5k exclusion rests on, recomputed
    struct = (bipartition_classes_are_the_only_independent_15_sets()
              and minimum_dominating_sets_are_perfect())
    upper, lower = {}, {}
    for k, lab in ((3, f3_label), (4, f4_label), (5, f5_label)):
        f = build(lab)
        upper[k] = weight(f) if is_krdf(f, k)[0] else None
        tri = feasible_triples(k, 5 * k)
        if tri != [(5 * k, 15, 0)]:
            lower[k] = None
            continue
        if k < 5 and struct and not extremal_covering_exists(k):
            lower[k] = 5 * k + 1
        else:
            lower[k] = 5 * k
    for k in (3, 4, 5):
        DERIVED[k] = upper[k] if (upper[k] is not None
                                  and upper[k] == lower[k]) else None
    ck("lower and upper bounds meet, giving (gamma_r3, gamma_r4, gamma_r5) "
       "= (16, 21, 25) as the paper claims",
       (DERIVED[3], DERIVED[4], DERIVED[5]) == (16, 21, 25),
       "derived lower=%s upper=%s" % (lower, upper))

    ceilings = {k: -((-k * 30) // 6) for k in (3, 4, 5)}
    ck("the derived weight floor 5k agrees with the cited cubic bound "
       "ceil(kN/6) for N=30, k=3,4,5",
       ceilings == {3: 15, 4: 20, 5: 25},
       "ceil(kN/6) = %s" % (ceilings,))


def ck_refutation():
    n = Fraction(len(VERTS))
    g3, g4, g5 = DERIVED[3], DERIVED[4], DERIVED[5]
    hyp = (g5 is not None and g5 == Fraction(5, 6) * n)
    concl = (g4 is not None and g4 == Fraction(2, 3) * n)
    strict = (g4 is not None and g4 > Fraction(2, 3) * n)
    ck("the first implication of the question fails: the hypothesis "
       "gamma_r5 = 5N/6 holds while gamma_r4 = 2N/3 fails, strictly above",
       hyp and (not concl) and strict,
       "N=%s, 5N/6=%s, gamma_r5=%s ; 2N/3=%s, gamma_r4=%s"
       % (n, Fraction(5, 6) * n, g5, Fraction(2, 3) * n, g4))

    ck("the second implication is untouched by this graph: its hypothesis "
       "gamma_r4 = 2N/3 fails here, and gamma_r3 = 16 differs from N/2 = 15",
       (not concl) and g3 is not None and g3 != Fraction(1, 2) * n,
       "gamma_r4=%s vs 2N/3=%s ; gamma_r3=%s vs N/2=%s"
       % (g4, Fraction(2, 3) * n, g3, Fraction(1, 2) * n))


DISCLOSURE = (
    "NOT RE-RUN HERE: "
    "(1) the three double-counting premises behind the lower bounds -- "
    "tau >= k(30-c) at the empty vertices, tau <= 3w-2e, and the cut bound "
    "3c-2e <= 3(30-c) -- are transcribed from the paper as premises, and only "
    "their integer consequences are scanned exhaustively; they are not "
    "re-derived symbolically here. "
    "(2) No solver and no exhaustive search over all labellings was run, so "
    "the paper's report that its 0/1 program was computed to optimality with "
    "optima 16, 21 and 25 is not reproduced by an independent optimisation: "
    "the upper bounds here come from the three exhibited labellings and the "
    "lower bounds from those transcribed premises plus finite enumeration. "
    "(3) The cited cubic bound gamma_rk(G) >= ceil(k|V(G)|/6) for all cubic G "
    "and all k <= 6, and the cited equality lemma in its published "
    "generality, are not verified beyond the order-30, k=3,4,5 instances used "
    "here. "
    "(4) Of the identification as the Tutte 8-cage, only girth 8 is verified; "
    "the quoted uniqueness of the cubic graph of girth 8 on 30 vertices is "
    "not. "
    "(5) The self-duality of the duad-syntheme configuration is not verified; "
    "its only use, the perfect-class statement on the syntheme side, is "
    "verified directly instead. "
    "(6) No bibliographic or attribution claim is checked: the wording of the "
    "question being refuted and the |V|=2n reading of it, the cited remark "
    "that this graph is not 3-rainbow domination regular, the claim that no "
    "earlier work determines gamma_r4 of this graph, and the paper's own "
    "acknowledgement that its literature search was not exhaustive. "
    "(7) The paper's account of its own process is not checked: that the 93 "
    "elementary conditions were done by hand and independently re-done, and "
    "that the machine cross-checks the paper reports were carried out there as "
    "described. This program is an independent re-derivation of those "
    "cross-checks, not a transcript of the paper's own run, and no such "
    "transcript accompanies it. "
    "(8) Nothing is checked about graphs other than this one, including the "
    "second implication of the question for other cubic graphs and any "
    "characterisation of the cubic graphs attaining the bound."
)


def main():
    ck_model()
    ck_graph()
    ck_girth()
    ck_upper(5, f5_label, 25, "d -> d cap [5]")
    ck_upper(4, f4_label, 21, "d -> d cap [4], with one exceptional duad")
    ck_upper(3, f3_label, 16, "the four-case labelling")
    ck_support_sizes()
    ck_exception_needed()
    ck_perfect_classes()
    ck_cover_lemma()
    ck_counting_bound()
    ck_independent_sets()
    ck_extremal_impossible()
    ck_program_model()
    ck_ninetythree()
    ck_values()
    ck_refutation()
    finish(DISCLOSURE)


if __name__ == "__main__":
    main()

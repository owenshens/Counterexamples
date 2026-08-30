#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- re-derives every numerical claim of

    "A 24-cell weak 4-resolving set of K_11 x K_11 refutes a surmised exact value"

from the objects PRINTED IN THE PAPER and from nothing else.

Python 3.9 or later, standard library only: no third-party package, no external
data file, no solver, no network.  All arithmetic is on Python integers; there
is no floating point anywhere and no decision is taken on a float.

WHAT IT READS.  The cell lists below are transcribed literally from the paper:
S (24 cells, n = 11) from the statement of Theorem 1.1, and S_12, S_13, S_14
from Section 5.  The block-diagonal sets S(m) are rebuilt here from S by the
rule stated in Section 5.  The one number quoted from the source paper that this program takes
on trust is wdim_4(K_3 x K_3) = 6, used only in the sharpness remark of
Section 4; every other quantity is recomputed.

HOW IT CHECKS.  The primary route is the RAW DISTANCE RULE of Section 1 --
d(u,v) = 0 if u = v, 1 if both coordinates differ, 2 otherwise -- summed term
by term over every unordered pair of vertices.  Lemma 2.1 (the closed formula
for Delta_S) and Lemma 2.2 (the two-line criterion) are then checked AGAINST
that raw route rather than assumed, so a reader who distrusts either lemma can
read only the raw checks and still see the result.

Output: one `PASS <name> [detail]` line per check, then

    VERDICT: ALL <n> CHECKS PASS

with exit status 0 if and only if every check passed.  The program closes with
its own statement of what it does NOT cover (`NOT RE-RUN: ...`), repeated in
REVIEW_NOTE.md under `## Scope`.
"""

import sys
from collections import deque
from itertools import combinations

# ---------------------------------------------------------------------------
# the check harness
# ---------------------------------------------------------------------------
_PASSES = 0
_FAILS = []


def chk(name, got, want, detail=''):
    """One check.  Prints PASS/FAIL and records it."""
    global _PASSES
    if got == want:
        _PASSES += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILS.append(name)
        print('FAIL %s got=%r want=%r%s'
              % (name, got, want, (' [%s]' % detail) if detail else ''))


# ---------------------------------------------------------------------------
# THE OBJECTS, exactly as printed in the paper
# ---------------------------------------------------------------------------
# Theorem 1.1: the 24 cells of S in K_11 x K_11
S11 = [(1, 1), (1, 10), (2, 1), (2, 11), (3, 2), (3, 10), (4, 2), (4, 11),
       (5, 3), (5, 10), (6, 3), (6, 11), (7, 4), (7, 5), (8, 6), (8, 7),
       (9, 8), (9, 9), (10, 4), (10, 6), (10, 8), (11, 5), (11, 7), (11, 9)]

# Section 5: S_12
S12 = [(1, 11), (1, 12), (2, 5), (2, 9), (3, 7), (3, 10), (3, 12), (4, 2),
       (4, 4), (5, 1), (5, 3), (5, 6), (6, 6), (6, 8), (7, 2), (7, 8),
       (7, 11), (8, 2), (8, 9), (9, 4), (9, 9), (10, 1), (10, 10), (11, 3),
       (11, 7), (12, 2), (12, 5)]

# Section 5: S_13
S13 = [(1, 6), (1, 12), (2, 7), (2, 13), (3, 1), (3, 2), (3, 5), (4, 6),
       (4, 10), (5, 3), (5, 5), (6, 10), (6, 12), (7, 1), (7, 7), (8, 8),
       (8, 11), (9, 9), (9, 10), (10, 4), (10, 8), (10, 9), (11, 3), (11, 7),
       (11, 11), (12, 2), (12, 4), (13, 12), (13, 13)]

# Section 5: S_14
S14 = [(1, 12), (1, 13), (2, 3), (2, 8), (3, 1), (3, 8), (4, 6), (4, 9),
       (4, 14), (5, 8), (5, 10), (6, 2), (6, 4), (7, 4), (7, 10), (8, 1),
       (8, 2), (8, 5), (9, 7), (9, 11), (9, 12), (10, 5), (10, 11), (11, 4),
       (11, 9), (12, 13), (12, 14), (13, 3), (13, 13), (14, 6), (14, 7)]

PRINTED = [('S11', 11, S11), ('S12', 12, S12), ('S13', 13, S13),
           ('S14', 14, S14)]

WDIM4_K3xK3_FROM_SOURCE = 6      # the one value quoted from arXiv:2605.22307v1


# ---------------------------------------------------------------------------
# the graph, the distance rule, and the weak-k predicate
# ---------------------------------------------------------------------------
def cells(n):
    return [(i, j) for i in range(1, n + 1) for j in range(1, n + 1)]


def dist_rule(u, v):
    """The distance function of K_n x K_n, Section 1 of the paper."""
    if u == v:
        return 0
    return 1 if (u[0] != v[0] and u[1] != v[1]) else 2


def bfs_distances(n):
    """Distances in K_n x K_n built edge by edge from the DEFINITION of the
    direct product, with no closed formula: (i,j) ~ (i',j') iff i != i' and
    j != j'.  Returns {vertex: {vertex: distance}}."""
    V = cells(n)
    adj = dict((u, []) for u in V)
    for u in V:
        for v in V:
            if u[0] != v[0] and u[1] != v[1]:
                adj[u].append(v)
    out = {}
    for s in V:
        d = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for w in adj[u]:
                if w not in d:
                    d[w] = d[u] + 1
                    q.append(w)
        out[s] = d
    return out


def degrees(n, S):
    a = dict((i, 0) for i in range(1, n + 1))
    b = dict((j, 0) for j in range(1, n + 1))
    for (i, j) in S:
        a[i] += 1
        b[j] += 1
    return a, b


def min_delta_raw(n, S):
    """(min Delta_S over all unordered pairs, number of pairs attaining it,
    number of pairs, one argmin) -- computed from dist_rule alone."""
    V = cells(n)
    prof = dict((v, tuple(dist_rule(v, s) for s in S)) for v in V)
    best = None
    tight = 0
    arg = None
    for x, y in combinations(V, 2):
        px = prof[x]
        py = prof[y]
        t = 0
        for p, q in zip(px, py):
            t += p - q if p > q else q - p
        if best is None or t < best:
            best = t
            tight = 1
            arg = (x, y)
        elif t == best:
            tight += 1
    return best, tight, len(V) * (len(V) - 1) // 2, arg


def delta_formula(n, S, x, y):
    """Lemma 2.1 of the paper."""
    Ss = set(S)
    a, b = degrees(n, S)
    (i, j), (i2, j2) = x, y
    al = (x in Ss) + (y in Ss)
    if i == i2:
        return b[j] + b[j2] + al
    if j == j2:
        return a[i] + a[i2] + al
    be = ((i, j2) in Ss) + ((i2, j) in Ss)
    return a[i] + a[i2] + b[j] + b[j2] - al - 2 * be


def formula_mismatches(n, S):
    """Number of unordered pairs where Lemma 2.1 disagrees with the raw sum."""
    Ss = set(S)
    a, b = degrees(n, S)
    V = cells(n)
    bad = 0
    for x, y in combinations(V, 2):
        (i, j), (i2, j2) = x, y
        al = (x in Ss) + (y in Ss)
        if i == i2:
            f = b[j] + b[j2] + al
        elif j == j2:
            f = a[i] + a[i2] + al
        else:
            be = ((i, j2) in Ss) + ((i2, j) in Ss)
            f = a[i] + a[i2] + b[j] + b[j2] - al - 2 * be
        raw = 0
        for s in S:
            p = dist_rule(x, s)
            q = dist_rule(y, s)
            raw += p - q if p > q else q - p
        if raw != f:
            bad += 1
    return bad


def criterion(n, S):
    """Lemma 2.2: with all degrees >= 2, S is weak 4-resolving iff every
    three-edge path of G_S has degree sum >= 9 and every 4-cycle has degree sum
    >= 10.  Returns (verdict, witnesses-of-failure)."""
    a, b = degrees(n, S)
    bad = []
    if min(a.values()) < 2 or min(b.values()) < 2:
        return False, [('degree<2',)]
    rows = dict((i, set()) for i in range(1, n + 1))
    cols = dict((j, set()) for j in range(1, n + 1))
    for (i, j) in S:
        rows[i].add(j)
        cols[j].add(i)
    for i, i2 in combinations(range(1, n + 1), 2):          # 4-cycles
        common = sorted(rows[i] & rows[i2])
        for j, j2 in combinations(common, 2):
            T = a[i] + a[i2] + b[j] + b[j2]
            if T < 10:
                bad.append(('C4', i, i2, j, j2, T))
    for i in range(1, n + 1):                               # three-edge paths
        for j, j2 in combinations(sorted(rows[i]), 2):
            for jm, jo in ((j, j2), (j2, j)):
                for i2 in cols[jo]:
                    if i2 == i:
                        continue
                    T = a[i] + a[i2] + b[jm] + b[jo]
                    if T < 9:
                        bad.append(('P4', i, jm, jo, i2, T))
    return (not bad), bad[:4]


def light_components(n, S):
    """Components of the subgraph of G_S induced by the degree-2 vertices."""
    a, b = degrees(n, S)
    edges = [(i, j) for (i, j) in S if a[i] == 2 and b[j] == 2]
    par = {}

    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    # every light vertex is a component of its own until an edge joins it, so
    # that a light vertex both of whose edges go to heavy vertices is counted
    for i in range(1, n + 1):
        if a[i] == 2:
            find(('R', i))
    for j in range(1, n + 1):
        if b[j] == 2:
            find(('C', j))
    for (i, j) in edges:
        ri, cj = find(('R', i)), find(('C', j))
        if ri != cj:
            par[ri] = cj
    comp = {}
    for v in list(par):
        comp.setdefault(find(v), []).append(v)
    return edges, sorted(len(c) for c in comp.values())


def four_cycles(n, S):
    rows = dict((i, set()) for i in range(1, n + 1))
    for (i, j) in S:
        rows[i].add(j)
    out = []
    for i, i2 in combinations(range(1, n + 1), 2):
        common = sorted(rows[i] & rows[i2])
        for j, j2 in combinations(common, 2):
            out.append((i, i2, j, j2))
    return out


def ceil_div(p, q):
    return -((-p) // q)


def lower_bound(n):
    """Theorem 1.2 of the paper, for n >= 4."""
    return 2 * n + ceil_div(2 * n, 11)


def surmise(n):
    """The value surmised in Section 1: 2n + 1 + floor(n/4)."""
    return 2 * n + 1 + n // 4


def block_diagonal(m):
    """S(m), the m-fold block-diagonal copy of S, in K_(11m) x K_(11m)."""
    return sorted((i + 11 * t, j + 11 * t) for t in range(m) for (i, j) in S11)


def thm33_set(n):
    """The set constructed in the proof of Theorem 3.3 of arXiv:2605.22307v1,
    rebuilt from the description quoted in Section 6 of the paper: the n
    diagonal cells, (i, i+2) for 1 <= i <= n-2, the four cells (2,1), (3,2),
    (n-1,n-2), (n,n-1), and (i, i-1) for 6 <= i <= n-2 with i = 2 mod 4."""
    T = set((i, i) for i in range(1, n + 1))
    T |= set((i, i + 2) for i in range(1, n - 1))
    T |= {(2, 1), (3, 2), (n - 1, n - 2), (n, n - 1)}
    T |= set((i, i - 1) for i in range(6, n - 1) if i % 4 == 2)
    return sorted(T)


def k3_basis(n):
    """The weak 3-resolving set of size 2n from Theorem 3.2 of
    arXiv:2605.22307v1: the diagonal, the subdiagonal, and (1,n)."""
    return sorted(set([(i, i) for i in range(1, n + 1)]
                      + [(i + 1, i) for i in range(1, n)] + [(1, n)]))


# ===========================================================================
print('=== A: the distance rule of Section 1 is the true distance function ===')
# The rule is checked against breadth-first distances in the product built edge
# by edge from the definition -- the one place a Cartesian/direct mix-up would
# show up.
for n in (3, 4, 5, 6, 7):
    D = bfs_distances(n)
    V = cells(n)
    bad = sum(1 for u in V for v in V if D[u].get(v) != dist_rule(u, v))
    reach = min(len(D[u]) for u in V)
    chk('dist-rule-vs-bfs-n%d' % n, (bad, reach), (0, n * n),
        '%d ordered pairs, connected, diameter 2' % (n * n * n * n))

print()
print('=== B: Lemma 2.1 (the closed form for Delta_S) against raw sums ===')
for label, n, S in (('S11', 11, S11), ('S12', 12, S12), ('S13', 13, S13),
                    ('k3basis-n5', 5, k3_basis(5)),
                    ('full-n6', 6, cells(6))):
    npairs = (n * n) * (n * n - 1) // 2
    chk('delta-formula-%s' % label, formula_mismatches(n, S), 0,
        'n=%d |S|=%d, %d pairs, 0 mismatches' % (n, len(S), npairs))

print()
print('=== C: the printed sets are weak 4-resolving and beat the surmise ===')
RAW = {}
for label, n, S in PRINTED:
    chk('cells-distinct-%s' % label, len(set(S)), len(S), '%d cells' % len(S))
    chk('cells-in-range-%s' % label,
        all(1 <= i <= n and 1 <= j <= n for (i, j) in S), True, 'n=%d' % n)
    a, b = degrees(n, S)
    chk('min-degree-%s' % label, min(min(a.values()), min(b.values())), 2,
        'row degrees %s, column degrees %s'
        % (sorted(a.values()), sorted(b.values())))
    best, tight, npairs, arg = min_delta_raw(n, S)
    RAW[label] = (best, tight, npairs, arg)
    chk('mindelta-raw-%s' % label, best, 4,
        'min Delta_S = %d over all %d pairs, %d tight, argmin %s'
        % (best, npairs, tight, arg))
    ok, bad = criterion(n, S)
    chk('criterion-agrees-%s' % label, (ok, bad), (best >= 4, []),
        'Lemma 2.2 verdict = raw verdict')
    chk('beats-surmise-%s' % label,
        (len(S), 2 * n + n // 4, surmise(n)),
        (2 * n + n // 4, 2 * n + n // 4, len(S) + 1),
        '|S| = %d = 2n+floor(n/4), surmised %d, short by 1'
        % (len(S), surmise(n)))

print()
print('=== D: the structure Section 3 argues from, recounted ===')
a11, b11 = degrees(11, S11)
chk('degseq-rows-S11', [a11[i] for i in range(1, 12)],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3], 'as printed in Section 3')
chk('degseq-cols-S11', [b11[j] for j in range(1, 12)],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3], 'as printed in Section 3')
chk('degsum-S11', (sum(a11.values()), sum(b11.values())), (24, 24),
    '9*2 + 2*3 = 24 both ways')
chk('no-4-cycle-S11', four_cycles(11, S11), [],
    'no two rows share two columns, so condition (C) is vacuous')
_edges, _sizes = light_components(11, S11)
chk('light-edge-count-S11', len(_edges), 12,
    'cells in a light row and a light column')
chk('light-components-S11', _sizes, [3, 3, 3, 3, 3, 3],
    'six two-edge paths, so no three-edge path is all-light')
chk('heavy-vertices-S11',
    sorted([('R', i) for i in range(1, 12) if a11[i] >= 3]
           + [('C', j) for j in range(1, 12) if b11[j] >= 3]),
    [('C', 10), ('C', 11), ('R', 10), ('R', 11)], 'r10 r11 c10 c11')

print()
print('=== E: the tight pairs quoted in the Remark of Section 3 ===')
chk('tight-pair-count-S11', RAW['S11'][1], 618,
    '618 of 7260 pairs attain Delta_S = 4, none below')
chk('tight-pair-1-1-vs-2-10',
    (delta_formula(11, S11, (1, 1), (2, 10)),
     sum(abs(dist_rule((1, 1), s) - dist_rule((2, 10), s)) for s in S11)),
    (4, 4), 'alpha=1, beta=2, T=9, Delta=4')
chk('tight-pair-invariants-1-1-2-10',
    ((1, 1) in set(S11), (2, 10) in set(S11),
     ((1, 10) in set(S11)) + ((2, 1) in set(S11)),
     a11[1] + a11[2] + b11[1] + b11[10]),
    (True, False, 2, 9), 'the (alpha,beta,T) = (1,2,9) case of Lemma 2.2')
chk('same-column-pair-7-4-vs-8-4',
    (delta_formula(11, S11, (7, 4), (8, 4)),
     sum(abs(dist_rule((7, 4), s) - dist_rule((8, 4), s)) for s in S11)),
    (5, 5), 'a_7 + a_8 + alpha = 2 + 2 + 1')

print()
print('=== F: Theorem 1.2, the lower bound: its arithmetic ===')
for n in (9, 10, 11, 12, 13, 14):
    chk('lower-bound-n%d' % n, lower_bound(n),
        {9: 20, 10: 22, 11: 24, 12: 27, 13: 29, 14: 31}[n],
        '2n + ceil(2n/11) = %d, surmised %d' % (lower_bound(n), surmise(n)))
chk('step2-case-dominates-bound',
    all(3 * n - 2 >= lower_bound(n) for n in range(4, 301)), True,
    'a degree <= 1 forces |S| >= 3n-2 >= 2n+ceil(2n/11), n = 4..300')
_bad = []
for n in range(4, 301):
    for E in range(0, ceil_div(2 * n, 11)):
        # Step 6: h <= 2E and 2n <= 4h + 3E would be needed; the maximum of
        # 4h + 3E over h <= 2E is 11E, and 11E < 2n whenever E < ceil(2n/11).
        if 11 * E >= 2 * n:
            _bad.append((n, E))
chk('step6-excludes-every-smaller-excess', _bad, [],
    'for n = 4..300 and every E < ceil(2n/11): max(4h+3E : h <= 2E) = 11E < 2n')
chk('hypothesis-n-at-least-4-is-needed',
    (lower_bound(3), WDIM4_K3xK3_FROM_SOURCE, min(3, 4), min(4, 4)),
    (7, 6, 3, 4),
    'at n=3 the bound would read 7 > 6 = wdim_4(K_3 x K_3) as reported by the '
    'source, and Step 1 there yields only min(n,4) = 3, not 4')

print()
print('=== G: Theorem 1.3, the exact values ===')
for label, n, S in PRINTED:
    chk('exact-value-n%d' % n, (lower_bound(n), len(S)),
        (len(S), 2 * n + n // 4),
        'lower bound = |S| = %d = 2n+floor(n/4); surmised %d'
        % (len(S), surmise(n)))
chk('exact-values-vs-surmised',
    ([lower_bound(n) for n in (11, 12, 13, 14)],
     [surmise(n) for n in (11, 12, 13, 14)]),
    ([24, 27, 29, 31], [25, 28, 30, 32]),
    'each surmised value is too large by exactly 1')

print()
print('=== H: the block-diagonal family S(m) ===')
for m in (1, 2, 3, 4, 5, 6):
    BD = block_diagonal(m)
    n = 11 * m
    chk('bd-size-m%d' % m, (len(BD), len(set(BD))), (24 * m, 24 * m),
        'n=%d, |S(m)| = 24m = %d distinct cells' % (n, 24 * m))
    ok, bad = criterion(n, BD)
    chk('bd-criterion-m%d' % m, (ok, bad), (True, []),
        'Lemma 2.2 holds: no 4-cycle, every three-edge path meets a degree-3 '
        'vertex')
    chk('bd-meets-lower-bound-m%d' % m, lower_bound(n), 24 * m,
        '2(11m) + ceil(22m/11) = 24m, so the construction is optimal')
    chk('bd-beats-surmise-m%d' % m, (24 * m < surmise(n), surmise(n) - 24 * m),
        (True, 1 + (11 * m) // 4 - 2 * m),
        'surmised %d, true %d, gap %d' % (surmise(n), 24 * m,
                                          surmise(n) - 24 * m))
for m in (1, 2, 3):
    BD = block_diagonal(m)
    n = 11 * m
    best, tight, npairs, arg = min_delta_raw(n, BD)
    chk('bd-mindelta-raw-m%d' % m, best, 4,
        'n=%d |S|=%d: min Delta = %d over all %d pairs, %d tight'
        % (n, len(BD), best, npairs, tight))
chk('bd-m1-is-S11', block_diagonal(1), sorted(S11),
    'S(1) is literally the set of Theorem 1.1')
chk('bd-m2-is-S11-plus-shift', block_diagonal(2),
    sorted(list(S11) + [(i + 11, j + 11) for (i, j) in S11]),
    'S(2) = S union (S shifted by (+11,+11)), 48 cells, no cross-block cell')
chk('bd-no-cross-block-cell-m3',
    sorted(set((i - 1) // 11 == (j - 1) // 11 for (i, j) in block_diagonal(3))),
    [True], 'every cell of S(3) lies in a diagonal 11x11 block')

print()
print('=== I: the published upper bound of arXiv:2605.22307v1 is TRUE ===')
# Rebuilt from the construction in the proof of its Theorem 3.3 and checked
# against the raw distance rule.  This is the check that says the present note
# is not an erratum.
for n in range(9, 15):
    T = thm33_set(n)
    chk('thm33-size-n%d' % n, len(T), surmise(n),
        'the paper\'s own set has 2n+1+floor(n/4) = %d cells' % surmise(n))
for n in range(9, 15):
    T = thm33_set(n)
    best = min_delta_raw(n, T)[0]
    chk('thm33-resolving-n%d' % n, best, 4,
        'min Delta = %d, so wdim_4 <= %d holds as published' % (best, len(T)))

print()
print('=== J: controls -- the same code must accept a proved positive and '
      'reject a proved negative ===')
for n in (4, 5, 6, 8):
    best = min_delta_raw(n, cells(n))[0]
    chk('control-full-set-n%d' % n, best, min(2 * n + 2, 4 * n - 6),
        'S = V gives min Delta = %d = kappa(K_n x K_n) = min(2n+2, 4n-6)'
        % best)
for n in (5, 7, 9, 11):
    B = k3_basis(n)
    best = min_delta_raw(n, B)[0]
    chk('control-k3-basis-n%d' % n, (len(B), best), (2 * n, 3),
        'the source\'s proved weak 3-basis: |S| = 2n, min Delta = 3 exactly, '
        'so it PASSES k=3 and FAILS k=4')

# ---------------------------------------------------------------------------
print()
print('NOT RE-RUN: n = 9 and n = 10 are NOT decided here.  This program '
      'computes the lower bounds 20 and 22 and confirms the published upper '
      'bounds 21 and 23, and performs no search in between; the surmised '
      'equality may well be correct at those two values of n.')
print('NOT RE-RUN: no exhaustive or solver search over candidate sets is '
      'performed at any n.  Minimality comes only from Theorem 1.2, whose '
      'combinatorial steps are hand proofs; this program checks that '
      'theorem\'s ARITHMETIC (the chain 2n <= 4h+3E <= 11E for n = 4..300) '
      'and not its Steps 1, 4 and 5.')
print('NOT RE-RUN: the block-diagonal family is verified against the raw '
      'distance rule only at m = 1, 2, 3, and against the Lemma 2.2 criterion '
      'at m = 1..6.  The statement for all m >= 1 is the hand proof of '
      'Section 5.')
print('NOT RE-RUN: no text of arXiv:2605.22307v1 is fetched or hashed here, '
      'so the quotations, line and byte locators, statement numbers and file '
      'digests of Section 1 are NOT machine-checked by this program.  The one '
      'numerical value taken from that source on trust is '
      'wdim_4(K_3 x K_3) = 6.')
print('NOT RE-RUN: the constraint-programming runs that originally found '
      'S_12, S_13, S_14 and a second 24-cell set at n = 11 are not repeated; '
      'they are needed only as a source of upper-bound witnesses, and the '
      'witnesses shipped are re-verified above from scratch.  No solver is '
      'imported.')
print('NOT RE-RUN: exact values for n outside {11, 12, 13, 14} and outside '
      'the multiples of 11 are neither computed nor bounded better than by '
      '2n+ceil(2n/11) <= wdim_4 <= 2n+1+floor(n/4).')

print()
if _FAILS:
    print('VERDICT: %d CHECKS FAILED: %s' % (len(_FAILS), ', '.join(_FAILS)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _PASSES)
sys.exit(0)

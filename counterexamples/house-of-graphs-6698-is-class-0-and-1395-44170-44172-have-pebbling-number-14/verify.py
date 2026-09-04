#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Verification program for "Pebbling numbers of four cubic 12-vertex graphs of
girth 3 and diameter 3: one is Class 0, the other three are 14".

Python 3.9+, STANDARD LIBRARY ONLY (no numpy / sympy / networkx and no external
data file), exact integer arithmetic throughout.  No floating point is used for
any decision; the only floats printed are wall-clock timings.

--------------------------------------------------------------------------
INPUTS EMBEDDED IN THIS FILE (these are NOT checks by themselves)
--------------------------------------------------------------------------
  ADJ            the four adjacency lists printed in Section 2 of the paper
                 (House of Graphs #1395, #6698, #44170, #44172; 0-indexed with
                 API_index = HoG_label - 1).
  J3_FIG_EDGES   the edge list of the Flower snark J_3 read off the source
                 paper's own Figure (labels A..L -> 0..11); J_3 is a CONTROL,
                 not a claim of this paper.
  CERTS          the eight exhibited configurations: the three 13-pebble
                 witnesses W1/W2/W3 (Theorem 1), two further 12-pebble
                 certificates, and the source paper's own three 12-pebble
                 configurations -- each with its target r and the claimed size
                 of its move-closed reachable set.
  CLAIM_*        the invariant table (triangles, 4-cycles, |Aut|, vertex
                 orbits), the per-target unsolvable-configuration counts at
                 totals 12 and 13, and the pebbling numbers, all taken from the
                 paper.
  CENSUS         16 graph6 strings, transcribed into this file only (the paper
                 prints no graph6 strings and claims no sixteen-member family),
                 each with the pebbling number this program re-derives for it.
Every one of these is compared against an independently recomputed value.  The
graphs themselves are inputs (they are the objects under study, printed in the
paper); everything asserted ABOUT them is recomputed here from the adjacency
lists, from scratch.

--------------------------------------------------------------------------
WHAT IS DERIVED HERE (the checks)
--------------------------------------------------------------------------
  * each printed adjacency list is validated as a simple symmetric cubic
    12-vertex graph and its connectivity, girth, radius and diameter are
    recomputed (every ambient clause of the source's sentence);
  * triangles and 4-cycles are counted TWICE by structurally different routes
    (direct enumeration over all C(12,4)=495 four-subsets x 3 pairings, and the
    trace formula (tr A^4 - 2m - 4*sum_i C(d_i,2))/8), and the "6 and 8" that a
    naive pair-sum reports is exhibited as exactly twice the true count;
  * |Aut| and the vertex orbits are computed by explicit backtracking
    enumeration of all automorphisms; the five graphs are shown pairwise
    non-isomorphic by exact backtracking isomorphism tests;
  * #1395 is shown isomorphic to a from-scratch truncation of K_4 and the
    figure's J_3 to a from-scratch Flower-snark construction AND to Petersen
    with one vertex expanded into a triangle (Tietze's graph);
  * every exhibited configuration is certified r-unsolvable with NO recourse to
    the exhaustive sweep: the full forward reachable set is generated, shown
    closed under pebbling moves, shown to contain no configuration with a
    pebble on r, and its size compared with the paper's;
  * the exact pebbling numbers are re-derived by an exhaustive level-indexed
    search that uses only the definition plus the two lemmas proved in the
    paper (Lemma 1 downward closure, Lemma 2 the one-move recurrence), and the
    completeness of the level-by-level generation is cross-checked at low
    levels against literal brute-force enumeration of EVERY configuration,
    whose visited count is compared with C(s+10,10) exactly;
  * the per-target unsolvable-count vectors, the totals, the pebbling numbers,
    the emptiness of level 14 for the three graphs, the emptiness of level 12
    for #6698, and the pi-distribution {12:2, 13:5, 14:9} over the sixteen
    graph6 strings transcribed into this file are all recomputed.

WHAT IS **NOT** DERIVED HERE (stated again at the end of the run):
  * that the cubic 12-vertex girth-3 diameter-3 graphs number EXACTLY sixteen.
    This program re-derives only that the sixteen graph6 strings transcribed
    into it are sixteen pairwise non-isomorphic cubic 12-vertex girth-3
    diameter-3 graphs -- i.e. "at least sixteen".  No completeness claim is
    made or supported here, and the note claims none.
  * what the House of Graphs database holds; only the four adjacency lists it
    served are used, and they are printed in the paper.

USAGE
  python3 verify.py            full run (the exhaustive sweep included)
  python3 verify.py --quick    skip the exhaustive sweep (Steps 5 and 6b only);
                               NOT a complete verification -- it drops the
                               upper bounds pi <= 14 and pi(#6698) <= 12.
The full run is parallel over (graph, target) tasks via multiprocessing from
the standard library; on 32 cores it takes a few minutes, single-core rather
longer.  Use -jN to fix the worker count.
"""

import itertools
import math
import os
import sys
import time

RESULTS = []


def check(name, cond, detail=''):
    cond = bool(cond)
    RESULTS.append((name, cond))
    print('%s %s%s' % ('PASS' if cond else 'FAIL', name, (' [%s]' % detail) if detail else ''))
    sys.stdout.flush()
    return cond


def note(s):
    print('NOTE ' + s)
    sys.stdout.flush()


def section(t):
    print('')
    print('=== ' + t)
    sys.stdout.flush()


# ======================================================================
# INPUTS, copied from the paper
# ======================================================================

ADJ = {
    '#1395': [[1, 2, 3], [0, 8, 9], [0, 3, 11], [0, 2, 10],
              [5, 6, 10], [4, 7, 11], [4, 8, 10], [5, 9, 11],
              [1, 6, 9], [1, 7, 8], [3, 4, 6], [2, 5, 7]],
    '#6698': [[1, 2, 3], [0, 8, 9], [0, 3, 11], [0, 2, 10],
              [5, 8, 10], [4, 8, 11], [7, 9, 11], [6, 9, 10],
              [1, 4, 5], [1, 6, 7], [3, 4, 7], [2, 5, 6]],
    '#44170': [[1, 2, 3], [0, 10, 11], [0, 10, 11], [0, 8, 9],
               [5, 6, 8], [4, 7, 8], [4, 9, 10], [5, 9, 11],
               [3, 4, 5], [3, 6, 7], [1, 2, 6], [1, 2, 7]],
    '#44172': [[1, 2, 3], [0, 9, 11], [0, 10, 11], [0, 7, 8],
               [5, 6, 8], [4, 9, 10], [4, 9, 10], [3, 8, 11],
               [3, 4, 7], [1, 5, 6], [2, 5, 6], [1, 2, 7]],
}

# The source paper's own Figure for the Flower snark J_3, labels A..L -> 0..11.
J3_FIG_EDGES = ['AB', 'AC', 'AD', 'BE', 'BF', 'CG', 'CH', 'DI', 'DJ',
                'EK', 'FL', 'EF', 'GK', 'IK', 'HL', 'JL', 'HI', 'GJ']

NAMES5 = ['#1395', '#6698', '#44170', '#44172', 'J3']

# claimed invariants: (triangles, 4-cycles, |Aut|, number of vertex orbits)
CLAIM_INV = {
    '#1395': (4, 0, 24, 1),
    '#6698': (3, 0, 36, 2),
    '#44170': (1, 3, 12, 5),
    '#44172': (1, 4, 8, 7),
    'J3': (1, 0, 12, 3),
}
CLAIM_ORBITS = {
    '#44170': [[0, 10, 11], [1, 2], [3, 6, 7], [4, 5, 8], [9]],
    '#44172': [[0, 11], [1, 2], [3, 7], [4], [5, 6], [8], [9, 10]],
}

# (id, graph, target, configuration, claimed |Reach(C)|, what it is)
CERTS = [
    ('W1', '#1395', 0, [0, 0, 0, 0, 3, 4, 1, 3, 0, 1, 1, 0], 6358, 'refutation witness, total 13'),
    ('W2', '#44170', 1, [0, 0, 1, 1, 5, 5, 0, 0, 1, 0, 0, 0], 1912, 'refutation witness, total 13'),
    ('W3', '#44172', 5, [3, 1, 1, 0, 0, 0, 1, 7, 0, 0, 0, 0], 1711, 'refutation witness, total 13'),
    ('X1', '#44170', 9, [1, 5, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1], 198, 'further certificate, total 12'),
    ('X2', '#44172', 3, [0, 0, 1, 0, 1, 1, 1, 1, 0, 7, 0, 0], 313, 'further certificate, total 12'),
    ('S1', '#1395', 0, [0, 0, 0, 0, 0, 7, 3, 0, 0, 1, 1, 0], 1262, "source's own configuration, total 12"),
    ('S2', '#44170', 5, [5, 1, 1, 0, 0, 0, 0, 0, 0, 0, 5, 0], 1251, "source's own configuration, total 12"),
    ('S3', '#44172', 0, [0, 0, 0, 0, 7, 1, 1, 0, 0, 1, 1, 1], 483, "source's own configuration, total 12"),
]

CLAIM_PI = {'#1395': 14, '#6698': 12, '#44170': 14, '#44172': 14, 'J3': 12}

CLAIM_T12 = {
    '#1395': [165] * 12,
    '#44170': [0, 229, 229, 0, 91, 91, 0, 0, 91, 3, 0, 0],
    '#44172': [8, 0, 0, 94, 4, 153, 153, 94, 18, 5, 5, 8],
    '#6698': [0] * 12,
    'J3': [0] * 12,
}
CLAIM_T13 = {
    '#1395': [2] * 12,
    '#44170': [0, 20, 20, 0, 4, 4, 0, 0, 4, 0, 0, 0],
    '#44172': [0, 0, 0, 0, 0, 6, 6, 0, 0, 0, 0, 0],
    '#6698': [0] * 12,
    'J3': [0] * 12,
}
CLAIM_TOTALS = {'#1395': (1980, 24), '#44170': (734, 52), '#44172': (542, 12),
                '#6698': (0, 0), 'J3': (0, 0)}

# sixteen graph6 strings transcribed into this file: (graph6, pi, identification or None)
CENSUS = [
    ('K?`@E`gFCKEO', 12, 'J3'),
    ('K?`@Cr_T?[EO', 12, '#6698'),
    ('K?ABE`Ke@WEO', 13, None),
    ('K?ABEESU@WF?', 13, None),
    ('K?`@Eb_F?[EO', 13, None),
    ('K?`@E`gFCcCo', 13, None),
    ('K?`D@`gE_iW_', 13, None),
    ('K?`D@`ae?iH_', 14, '#1395'),
    ('K?AEDDWX@oB_', 14, '#44170'),
    ('K?ABCN_U?wF?', 14, '#44172'),
    ('K?ABE`Ke@gDO', 14, None),
    ('K?ABArCe@gDO', 14, None),
    ('K?ABCNOY?wF?', 14, None),
    ('K?BD?pck?sHG', 14, None),
    ('K?`@E`Wh?kDO', 14, None),
    ('K?`D@bAJAgBG', 14, None),
]
CLAIM_PI_DIST = {12: 2, 13: 5, 14: 9}

N = 12

# ======================================================================
# graph utilities (all from scratch)
# ======================================================================


def edges_of(adj):
    return sorted((min(u, v), max(u, v)) for u in range(len(adj)) for v in adj[u] if u < v)


def is_simple_symmetric(adj):
    n = len(adj)
    for u in range(n):
        if u in adj[u]:
            return False, 'loop at %d' % u
        if len(set(adj[u])) != len(adj[u]):
            return False, 'repeated neighbour at %d' % u
        for v in adj[u]:
            if not (0 <= v < n):
                return False, 'out-of-range neighbour %d' % v
            if u not in adj[v]:
                return False, 'asymmetric pair (%d,%d)' % (u, v)
    return True, ''


def bfs(adj, s):
    n = len(adj)
    d = [-1] * n
    d[s] = 0
    q = [s]
    while q:
        nxt = []
        for u in q:
            for v in adj[u]:
                if d[v] < 0:
                    d[v] = d[u] + 1
                    nxt.append(v)
        q = nxt
    return d


def connected(adj):
    return all(x >= 0 for x in bfs(adj, 0))


def ecc_profile(adj):
    """(radius, diameter)."""
    e = []
    for s in range(len(adj)):
        d = bfs(adj, s)
        if min(d) < 0:
            return None, None
        e.append(max(d))
    return min(e), max(e)


def girth(adj):
    """Shortest cycle length by BFS from every vertex over every removed edge."""
    n = len(adj)
    best = 10 ** 9
    for s in range(n):
        d = [-1] * n
        par = [-1] * n
        d[s] = 0
        q = [s]
        while q:
            nxt = []
            for u in q:
                for v in adj[u]:
                    if d[v] < 0:
                        d[v] = d[u] + 1
                        par[v] = u
                        nxt.append(v)
                    elif v != par[u]:
                        best = min(best, d[u] + d[v] + 1)
            q = nxt
    return best


def triangles(adj):
    nb = [set(a) for a in adj]
    c = 0
    for a, b, d in itertools.combinations(range(len(adj)), 3):
        if b in nb[a] and d in nb[a] and d in nb[b]:
            c += 1
    return c


def c4_direct(adj):
    """4-cycles counted by literal enumeration: every 4-subset, every 3 pairings."""
    nb = [set(a) for a in adj]
    c = 0
    for q in itertools.combinations(range(len(adj)), 4):
        a, b, x, y = q
        for cyc in ((a, b, x, y), (a, b, y, x), (a, x, b, y)):
            p, r, s, t = cyc
            if r in nb[p] and s in nb[r] and t in nb[s] and p in nb[t]:
                c += 1
    return c


def c4_trace(adj):
    """(tr A^4 - 2m - 4*sum_i C(d_i,2)) / 8, exact integer arithmetic."""
    n = len(adj)
    A = [[1 if v in adj[u] else 0 for v in range(n)] for u in range(n)]

    def mul(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

    A2 = mul(A, A)
    A4 = mul(A2, A2)
    tr = sum(A4[i][i] for i in range(n))
    m = sum(len(a) for a in adj) // 2
    s = sum(math.comb(len(adj[i]), 2) for i in range(n))
    num = tr - 2 * m - 4 * s
    assert num % 8 == 0, num
    return num // 8


def c4_pairsum(adj):
    """sum over unordered pairs {a,b} of C(|N(a) cap N(b)|, 2) -- the UNDIVIDED count."""
    nb = [set(a) for a in adj]
    return sum(math.comb(len(nb[a] & nb[b]), 2)
               for a, b in itertools.combinations(range(len(adj)), 2))


def wl_colors(adj):
    """1-WL refinement seeded with (degree, triangles through the vertex)."""
    nb = [set(a) for a in adj]
    col = []
    for u in range(len(adj)):
        t = sum(1 for x, y in itertools.combinations(sorted(nb[u]), 2) if y in nb[x])
        col.append((len(adj[u]), t))
    col = _normalise(col)
    for _ in range(len(adj)):
        new = [(col[u], tuple(sorted(col[v] for v in adj[u]))) for u in range(len(adj))]
        new = _normalise(new)
        if new == col:
            break
        col = new
    return col


def _normalise(col):
    order = {c: i for i, c in enumerate(sorted(set(col)))}
    return [order[c] for c in col]


def _search(adjG, adjH, want_all):
    """Backtracking search for isomorphisms G -> H.  Yields complete mappings."""
    n = len(adjG)
    if len(adjH) != n:
        return
    cg, ch = wl_colors(adjG), wl_colors(adjH)
    if sorted(cg) != sorted(ch):
        return
    nbG = [set(a) for a in adjG]
    nbH = [set(a) for a in adjH]
    # order the vertices of G so that each one after the first touches an earlier one
    order = [0]
    seen = {0}
    while len(order) < n:
        for u in order:
            for v in adjG[u]:
                if v not in seen:
                    seen.add(v)
                    order.append(v)
        for v in range(n):
            if v not in seen:
                seen.add(v)
                order.append(v)
                break
    cand = [[h for h in range(n) if ch[h] == cg[g]] for g in range(n)]
    img = [-1] * n
    used = [False] * n
    out = []

    def rec(i):
        if i == n:
            out.append(list(img))
            return not want_all       # stop after the first if we only want one
        g = order[i]
        for h in cand[g]:
            if used[h]:
                continue
            ok = True
            for j in range(i):
                gj = order[j]
                if (gj in nbG[g]) != (img[gj] in nbH[h]):
                    ok = False
                    break
            if ok:
                img[g] = h
                used[h] = True
                if rec(i + 1):
                    used[h] = False
                    img[g] = -1
                    return True
                used[h] = False
                img[g] = -1
        return False

    rec(0)
    for m in out:
        yield m


def automorphisms(adj):
    return list(_search(adj, adj, True))


def orbits_of(adj):
    auts = automorphisms(adj)
    parent = list(range(len(adj)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in auts:
        for v in range(len(adj)):
            a, b = find(v), find(p[v])
            if a != b:
                parent[a] = b
    groups = {}
    for v in range(len(adj)):
        groups.setdefault(find(v), []).append(v)
    return len(auts), sorted(sorted(g) for g in groups.values())


def isomorphic(adjG, adjH):
    for _ in _search(adjG, adjH, False):
        return True
    return False


def g6_decode(s):
    """graph6 -> adjacency lists (n < 63 only, which is all we need)."""
    b = [ord(c) - 63 for c in s]
    n = b[0]
    if not (0 <= n < 63):
        raise ValueError('graph6 order byte out of the supported range')
    bits = []
    for x in b[1:]:
        for k in range(5, -1, -1):
            bits.append((x >> k) & 1)
    adj = [[] for _ in range(n)]
    i = 0
    for col in range(1, n):
        for row in range(col):
            if i < len(bits) and bits[i]:
                adj[row].append(col)
                adj[col].append(row)
            i += 1
    return [sorted(a) for a in adj]


def build_j3_from_figure():
    adj = [[] for _ in range(N)]
    for e in J3_FIG_EDGES:
        u, v = ord(e[0]) - 65, ord(e[1]) - 65
        adj[u].append(v)
        adj[v].append(u)
    return [sorted(a) for a in adj]


def build_flower_snark_j3():
    """a_i ~ b_i,c_i,d_i;  triangle b0b1b2;  6-cycle c0 c1 c2 d0 d1 d2."""
    idx = {}
    for i in range(3):
        for k, nm in enumerate('abcd'):
            idx[(nm, i)] = 4 * i + k
    E = []
    for i in range(3):
        for nm in 'bcd':
            E.append((idx[('a', i)], idx[(nm, i)]))
    for i in range(3):
        E.append((idx[('b', i)], idx[('b', (i + 1) % 3)]))
    cyc = [('c', 0), ('c', 1), ('c', 2), ('d', 0), ('d', 1), ('d', 2)]
    for i in range(6):
        E.append((idx[cyc[i]], idx[cyc[(i + 1) % 6]]))
    adj = [[] for _ in range(N)]
    for u, v in E:
        adj[u].append(v)
        adj[v].append(u)
    return [sorted(a) for a in adj]


def build_petersen():
    adj = [[] for _ in range(10)]

    def add(u, v):
        adj[u].append(v)
        adj[v].append(u)

    for i in range(5):
        add(i, (i + 1) % 5)              # outer 5-cycle
        add(5 + i, 5 + (i + 2) % 5)      # inner pentagram
        add(i, 5 + i)                    # spokes
    return [sorted(set(a)) for a in adj]


def build_tietze():
    """Petersen with one vertex expanded into a triangle -- Tietze's graph."""
    P = build_petersen()
    v = 0
    nbrs = list(P[v])
    keep = [u for u in range(10) if u != v]
    ren = {u: i for i, u in enumerate(keep)}
    adj = [[] for _ in range(12)]

    def add(a, b):
        adj[a].append(b)
        adj[b].append(a)

    for u in keep:
        for w in P[u]:
            if w != v and u < w:
                add(ren[u], ren[w])
    t = [9, 10, 11]
    add(t[0], t[1])
    add(t[1], t[2])
    add(t[2], t[0])
    for k, u in enumerate(nbrs):
        add(t[k], ren[u])
    return [sorted(a) for a in adj]


def build_truncated_tetrahedron():
    """Truncation of K_4: vertices (u,v), u != v; a triangle on each (u,*);
    the matching (u,v)-(v,u) across each K_4 edge."""
    verts = [(u, v) for u in range(4) for v in range(4) if u != v]
    idx = {p: i for i, p in enumerate(verts)}
    adj = [[] for _ in range(12)]

    def add(a, b):
        adj[a].append(b)
        adj[b].append(a)

    for u in range(4):
        trio = [idx[(u, v)] for v in range(4) if v != u]
        add(trio[0], trio[1])
        add(trio[1], trio[2])
        add(trio[2], trio[0])
    for u in range(4):
        for v in range(u + 1, 4):
            add(idx[(u, v)], idx[(v, u)])
    return [sorted(a) for a in adj]


# ======================================================================
# pebbling
# ======================================================================

SH = 5                        # 5 bits per vertex; totals here never exceed 15
MASK = (1 << SH) - 1


def key_of(cfg):
    k = 0
    for i, c in enumerate(cfg):
        k |= c << (SH * i)
    return k


def cfg_of(key):
    return [(key >> (SH * i)) & MASK for i in range(N)]


def reach_closure(adj, cfg, r, cap=5 * 10 ** 6):
    """The full set of configurations reachable from cfg by pebbling moves.

    Returns (size, hits_target, closed).  `closed` re-checks, after the fact,
    that every move out of every member of the set lands back inside it --
    which together with hits_target=False is a self-contained proof that cfg is
    r-unsolvable, using no recurrence and no lemma.
    """
    start = key_of(cfg)
    seen = {start}
    stack = [start]
    hits = cfg[r] >= 1
    while stack:
        k = stack.pop()
        c = cfg_of(k)
        for v in range(N):
            if c[v] >= 2:
                for u in adj[v]:
                    k2 = k - (2 << (SH * v)) + (1 << (SH * u))
                    if k2 not in seen:
                        if u == r:
                            hits = True
                        seen.add(k2)
                        if len(seen) > cap:
                            raise RuntimeError('reachable set exceeded the cap')
                        stack.append(k2)
    closed = True
    for k in seen:
        c = cfg_of(k)
        for v in range(N):
            if c[v] >= 2:
                for u in adj[v]:
                    if (k - (2 << (SH * v)) + (1 << (SH * u))) not in seen:
                        closed = False
    return len(seen), hits, closed


def _unsolvable(k, c, r, adj, prev):
    """Lemma 2, applied to one configuration with c[r] == 0 at level s, given
    the complete unsolvable set `prev` at level s-1."""
    for v in range(N):
        if c[v] >= 2:
            for u in adj[v]:
                if u == r:
                    return False
                if (k - (2 << (SH * v)) + (1 << (SH * u))) not in prev:
                    return False
    return True


def sweep_target(adj, r, hard_stop=15):
    """Exhaustive level-indexed sweep for one target.

    Returns (first_empty_level, {level: count}).  Level s is generated from
    level s-1 by adding one pebble anywhere but r, which is COMPLETE by Lemma 1
    (downward closure), and each candidate is decided by Lemma 2.
    """
    prev = {0}                                  # level 0: the empty configuration
    counts = {0: 1}
    s = 0
    while prev:
        s += 1
        if s > hard_stop:
            raise RuntimeError('sweep did not terminate by level %d' % hard_stop)
        cur = set()
        tried = set()
        for k in prev:
            for v in range(N):
                if v == r:
                    continue
                k2 = k + (1 << (SH * v))
                if k2 in tried:
                    continue
                tried.add(k2)
                if _unsolvable(k2, cfg_of(k2), r, adj, prev):
                    cur.add(k2)
        counts[s] = len(cur)
        prev = cur
    return s, counts


def brute_level(adj, r, s):
    """Every configuration of total s with C(r)=0, decided by Lemma 2 from the
    brute-forced level below.  Returns (set, number of configurations visited),
    the visited count being comparable with C(s+10,10)."""
    free = [v for v in range(N) if v != r]
    levels = {0: {0}}
    visited = {}
    for t in range(1, s + 1):
        cur = set()
        seen = 0
        for comp in _compositions(t, len(free)):
            seen += 1
            k = 0
            for i, x in enumerate(comp):
                k |= x << (SH * free[i])
            if _unsolvable(k, cfg_of(k), r, adj, levels[t - 1]):
                cur.add(k)
        levels[t] = cur
        visited[t] = seen
    return levels, visited


def _compositions(total, parts):
    """All ordered tuples of `parts` non-negative integers summing to `total`."""
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, parts - 1):
            yield (first,) + rest


def _task(arg):
    name, adj, r, hard_stop = arg
    t0 = time.time()
    first_empty, counts = sweep_target(adj, r, hard_stop)
    return (name, r, first_empty, counts, time.time() - t0)


# ======================================================================
# the steps
# ======================================================================

GRAPHS = {}


def step1_objects(ctx):
    section('Step 1: the exhibited graphs, validated as objects')
    GRAPHS.update({k: [sorted(a) for a in v] for k, v in ADJ.items()})
    GRAPHS['J3'] = build_j3_from_figure()
    for nm in NAMES5:
        adj = GRAPHS[nm]
        ok, why = is_simple_symmetric(adj)
        m = sum(len(a) for a in adj) // 2
        degs = sorted(len(a) for a in adj)
        check('%s_is_a_simple_symmetric_cubic_graph_on_12_vertices_with_18_edges' % _sn(nm),
              ok and len(adj) == 12 and m == 18 and degs == [3] * 12,
              'n=%d m=%d degrees=%s%s' % (len(adj), m, set(degs), '' if ok else ' ' + why))
        rad, dia = ecc_profile(adj)
        g = girth(adj)
        check('%s_is_connected_with_girth_3_radius_3_diameter_3' % _sn(nm),
              connected(adj) and g == 3 and rad == 3 and dia == 3,
              'connected=%s girth=%d radius=%d diameter=%d' % (connected(adj), g, rad, dia))
    note('every ambient clause of the source sentence (cubic, 12 vertices, girth 3, diameter 3) '
         'holds for all five graphs, recomputed above')


def _sn(nm):
    return nm.replace('#', 'hog')


def step2_invariants(ctx):
    section('Step 2: invariants, each by two structurally different routes')
    for nm in NAMES5:
        adj = GRAPHS[nm]
        tri = triangles(adj)
        d1 = c4_direct(adj)
        d2 = c4_trace(adj)
        ps = c4_pairsum(adj)
        ctri, cc4, caut, corb = CLAIM_INV[nm]
        check('%s_triangle_count_is_%d' % (_sn(nm), ctri), tri == ctri, 'triangles=%d' % tri)
        check('%s_4_cycle_count_is_%d_by_direct_enumeration_and_by_the_trace_formula'
              % (_sn(nm), cc4), d1 == d2 == cc4,
              'direct=%d trace=%d pair_sum=%d' % (d1, d2, ps))
        check('%s_the_naive_pair_sum_is_exactly_twice_the_true_4_cycle_count' % _sn(nm),
              ps == 2 * d1, 'pair_sum=%d = 2 x %d' % (ps, d1))
        naut, orbs = orbits_of(adj)
        check('%s_has_automorphism_group_of_order_%d_and_%d_vertex_orbits'
              % (_sn(nm), caut, corb), naut == caut and len(orbs) == corb,
              '|Aut|=%d orbits=%d %s' % (naut, len(orbs), orbs))
        ctx.setdefault('orbits', {})[nm] = orbs
    for nm, claim in CLAIM_ORBITS.items():
        check('%s_vertex_orbit_partition_is_as_printed' % _sn(nm),
              ctx['orbits'][nm] == [sorted(g) for g in claim], str(ctx['orbits'][nm]))
    note('the 4-cycle counts of #44170 and #44172 are 3 and 4; the value 6 and 8 that a pair-sum '
         'reports is the undivided count, each 4-cycle being counted once per diagonal pair')


def step3_identify(ctx):
    section('Step 3: identification of the objects, and pairwise non-isomorphism')
    tt = build_truncated_tetrahedron()
    ok, _ = is_simple_symmetric(tt)
    check('the_from_scratch_truncation_of_k4_is_a_cubic_12_vertex_graph',
          ok and sorted(len(a) for a in tt) == [3] * 12, '%d edges' % (sum(len(a) for a in tt) // 2))
    check('hog1395_is_isomorphic_to_the_truncated_tetrahedron_built_from_scratch',
          isomorphic(GRAPHS['#1395'], tt))
    fs = build_flower_snark_j3()
    check('the_figure_j3_is_isomorphic_to_the_flower_snark_j3_built_from_the_standard_recipe',
          isomorphic(GRAPHS['J3'], fs), 'a_i ~ b_i,c_i,d_i; triangle on the b_i; 6-cycle c0c1c2d0d1d2')
    tz = build_tietze()
    check('the_figure_j3_is_isomorphic_to_petersen_with_one_vertex_expanded_to_a_triangle',
          isomorphic(GRAPHS['J3'], tz), "Tietze's graph")
    bad = [(a, b) for a, b in itertools.combinations(NAMES5, 2)
           if isomorphic(GRAPHS[a], GRAPHS[b])]
    check('the_five_graphs_are_pairwise_non_isomorphic', not bad,
          '10 pairs tested, %d isomorphic' % len(bad))
    note('so #6698 is not $J_3$ relabelled: its Class-0 status is a fact about a different graph')


def step4_certificates(ctx):
    section('Step 4: the exhibited configurations, certified with no exhaustive search')
    lower = {}
    for cid, nm, r, cfg, claim_reach, what in CERTS:
        adj = GRAPHS[nm]
        tot = sum(cfg)
        size, hits, closed = reach_closure(adj, cfg, r)
        check('%s_on_%s_has_total_%d_and_no_pebble_on_the_target' % (cid, _sn(nm), tot),
              len(cfg) == 12 and cfg[r] == 0, 'r=%d C=%s (%s)' % (r, cfg, what))
        check('%s_reachable_set_is_move_closed_misses_the_target_and_has_size_%d'
              % (cid, claim_reach), closed and not hits and size == claim_reach,
              '|Reach|=%d closed=%s pebble_on_target_ever=%s' % (size, closed, hits))
        lower[nm] = max(lower.get(nm, 0), tot + 1)
    for nm in ('#1395', '#44170', '#44172'):
        check('pi_%s_is_at_least_14_from_an_exhibited_13_pebble_certificate_alone' % _sn(nm),
              lower.get(nm) == 14, 'pi >= %d' % lower.get(nm))
    note('Theorem 1 (the refutation of the suggested value 13) is complete at this point and needs '
         'no exhaustive search: three closed reachable sets, checkable by hand')


def step5_lemmas(ctx):
    section('Step 5: the two lemmas the exhaustive sweep rests on, tested on data')
    # Lemma 1, downward closure, tested by sampling every one-pebble removal from a
    # certified unsolvable configuration and re-certifying it by forward reachability.
    bad = 0
    tested = 0
    for cid, nm, r, cfg, _cr, _w in CERTS[:5]:
        adj = GRAPHS[nm]
        for v in range(N):
            if cfg[v] >= 1:
                c2 = list(cfg)
                c2[v] -= 1
                _s2, hits, _c2 = reach_closure(adj, c2, r)
                tested += 1
                if hits:
                    bad += 1
    check('lemma_1_downward_closure_holds_on_every_one_pebble_removal_from_the_certificates',
          bad == 0, '%d removals re-certified unsolvable, %d failures' % (tested, bad))
    # Lemma 2, the one-move recurrence, tested against forward reachability at low levels.
    adj = GRAPHS['#44172']
    r = 5
    levels, visited = brute_level(adj, r, 6)
    ok_odo = all(visited[s] == math.comb(s + 10, 10) for s in visited)
    check('the_brute_force_odometer_visits_exactly_binom_s_plus_10_choose_10_configurations',
          ok_odo, ', '.join('s=%d:%d' % (s, visited[s]) for s in sorted(visited)))
    mism = 0
    for s in range(1, 6):
        for k in levels[s]:
            _sz, hits, _c = reach_closure(adj, cfg_of(k), r)
            if hits:
                mism += 1
    check('lemma_2_recurrence_and_forward_reachability_agree_on_every_configuration_to_level_5',
          mism == 0, '%d configurations, %d disagreements (graph #44172, target 5)'
          % (sum(len(levels[s]) for s in range(1, 6)), mism))
    ctx['brute'] = (levels, visited)


def step6a_sweep(ctx):
    section('Step 6a: the exhaustive level-indexed sweep (the exact pebbling numbers)')
    tasks = [(nm, GRAPHS[nm], r, 15) for nm in NAMES5 for r in range(12)]
    res = _run(tasks, ctx)
    per = {}
    for nm, r, first_empty, counts, _dt in res:
        per.setdefault(nm, {})[r] = (first_empty, counts)
    # completeness cross-check: the closure-generated sets equal the brute-forced ones
    levels, _visited = ctx['brute']
    _fe, counts = sweep_target(GRAPHS['#44172'], 5, 15)
    check('the_level_generation_reproduces_the_brute_force_level_sizes_up_to_level_6',
          all(counts[s] == len(levels[s]) for s in range(0, 7)),
          'sweep=%s brute=%s' % ([counts[s] for s in range(7)], [len(levels[s]) for s in range(7)]))
    for nm in NAMES5:
        pi = max(per[nm][r][0] for r in range(12))
        v12 = [per[nm][r][1].get(12, 0) for r in range(12)]
        v13 = [per[nm][r][1].get(13, 0) for r in range(12)]
        check('per_target_unsolvable_counts_of_%s_at_total_12_match_the_transcribed_table' % _sn(nm),
              v12 == CLAIM_T12[nm], str(v12))
        check('per_target_unsolvable_counts_of_%s_at_total_13_match_the_transcribed_table' % _sn(nm),
              v13 == CLAIM_T13[nm], str(v13))
        check('summed_unsolvable_counts_of_%s_are_%s_at_totals_12_and_13'
              % (_sn(nm), str(CLAIM_TOTALS[nm])),
              (sum(v12), sum(v13)) == CLAIM_TOTALS[nm], '(%d, %d)' % (sum(v12), sum(v13)))
        check('pi_%s_equals_%d_exactly' % (_sn(nm), CLAIM_PI[nm]), pi == CLAIM_PI[nm], 'pi=%d' % pi)
        empt = {r: per[nm][r][0] for r in range(12)}
        check('%s_every_target_reaches_an_empty_unsolvable_level_at_or_below_%d'
              % (_sn(nm), CLAIM_PI[nm]), max(empt.values()) == CLAIM_PI[nm],
              'first empty level per target = %s' % [empt[r] for r in range(12)])
        orbs = ctx['orbits'][nm]
        const = all(len(set(v12[x] for x in g)) == 1 and len(set(v13[x] for x in g)) == 1
                    for g in orbs)
        check('%s_per_target_counts_are_constant_on_automorphism_orbits' % _sn(nm), const,
              '%d orbits' % len(orbs))
    for nm in ('#1395', '#44170', '#44172'):
        allempty = all(per[nm][r][0] <= 14 for r in range(12))
        check('%s_has_no_unsolvable_configuration_of_total_14_for_any_target' % _sn(nm), allempty,
              'first empty level <= 14 for all 12 targets, so U_14 = U_15 = 0 by Lemma 1')
    for nm in ('#6698', 'J3'):
        z = all(per[nm][r][1].get(12, 0) == 0 for r in range(12))
        nz11 = all(per[nm][r][1].get(11, 0) > 0 for r in range(12))
        check('%s_has_no_unsolvable_configuration_of_total_12_for_any_target_so_it_is_class_0'
              % _sn(nm), z, '0 of 12 x C(22,10) = %d configurations' % (12 * math.comb(22, 10)))
        check('%s_does_have_unsolvable_configurations_of_total_11_so_pi_is_at_least_12' % _sn(nm),
              nz11, 'per-target counts at total 11 = %s'
              % [per[nm][r][1].get(11, 0) for r in range(12)])
    ctx['per'] = per


def step6b_census_structure(ctx):
    section('Step 6b: the sixteen transcribed graph6 strings, as objects')
    dec = []
    for g6, pi, ident in CENSUS:
        adj = g6_decode(g6)
        ok, why = is_simple_symmetric(adj)
        rad, dia = ecc_profile(adj)
        good = (ok and len(adj) == 12 and sorted(len(a) for a in adj) == [3] * 12
                and connected(adj) and girth(adj) == 3 and dia == 3)
        dec.append((g6, pi, ident, adj, good))
    check('all_sixteen_graph6_strings_decode_to_cubic_connected_12_vertex_girth_3_diameter_3_graphs',
          all(d[4] for d in dec), '%d of %d' % (sum(1 for d in dec if d[4]), len(dec)))
    bad = [(a[0], b[0]) for a, b in itertools.combinations(dec, 2) if isomorphic(a[3], b[3])]
    check('the_sixteen_decoded_graph6_graphs_are_pairwise_non_isomorphic', not bad,
          '%d pairs tested, %d isomorphic' % (len(dec) * (len(dec) - 1) // 2, len(bad)))
    hit = []
    for g6, pi, ident, adj, _g in dec:
        if ident:
            hit.append((ident, isomorphic(GRAPHS[ident], adj)))
    check('the_five_named_graphs_each_appear_among_the_sixteen_decoded_graph6_graphs',
          len(hit) == 5 and all(h[1] for h in hit), str([h[0] for h in hit]))
    check('the_transcribed_pi_values_of_the_sixteen_graph6_strings_form_the_multiset_12x2_13x5_14x9',
          {k: sum(1 for d in dec if d[1] == k) for k in CLAIM_PI_DIST} == CLAIM_PI_DIST
          and len(dec) == sum(CLAIM_PI_DIST.values()), str(CLAIM_PI_DIST))
    ctx['census'] = dec


def step6c_census_pi(ctx):
    section('Step 6c: the pebbling number of every one of the sixteen transcribed graph6 graphs')
    dec = ctx['census']
    tasks = []
    reps = {}
    for g6, pi, ident, adj, _g in dec:
        if ident:
            continue                      # pi is an isomorphism invariant; done in Step 6a
        _n, orbs = orbits_of(adj)
        reps[g6] = [g[0] for g in orbs]
        for r in reps[g6]:
            tasks.append((g6, adj, r, 15))
    note('pi is an isomorphism invariant and |U_s(phi(r))| = |U_s(r)| for phi in Aut(G), so the '
         'eleven unnamed transcribed graph6 graphs are swept on one target per vertex orbit '
         '(%d sweeps)' % len(tasks))
    res = _run(tasks, ctx)
    got = {}
    for g6, r, first_empty, _counts, _dt in res:
        got[g6] = max(got.get(g6, 0), first_empty)
    okall = True
    for g6, pi, ident, _adj, _g in dec:
        if ident:
            continue
        if got.get(g6) != pi:
            okall = False
    check('the_pi_value_of_each_of_the_eleven_unnamed_graph6_graphs_matches_the_transcribed_value', okall,
          ', '.join('%s:%d' % (g6[-4:], got[g6]) for g6, pi, i, a, g in dec if not i))
    check('the_five_named_graph6_strings_carry_the_pi_values_re_derived_in_step_6a',
          all(pi == CLAIM_PI[ident] for g6, pi, ident, a, g in dec if ident),
          str({ident: pi for g6, pi, ident, a, g in dec if ident}))
    dist = {}
    for g6, pi, ident, _a, _g in dec:
        dist[pi] = dist.get(pi, 0) + 1
    check('the_re_derived_pi_distribution_over_the_sixteen_graph6_graphs_is_12_twice_13_five_times_14_nine',
          dist == CLAIM_PI_DIST, str(dist))


def step7_arithmetic(ctx):
    section('Step 7: the arithmetic the paper prints')
    check('the_search_space_per_target_at_total_12_is_binom_22_10_equals_646646',
          math.comb(22, 10) == 646646, '%d' % math.comb(22, 10))
    check('the_search_space_per_target_at_total_13_is_binom_23_10_equals_1144066',
          math.comb(23, 10) == 1144066, '%d' % math.comb(23, 10))
    check('twelve_targets_times_binom_22_10_is_7759752',
          12 * math.comb(22, 10) == 7759752, '%d' % (12 * math.comb(22, 10)))
    check('three_graphs_times_twelve_targets_times_binom_23_10_is_41186376',
          3 * 12 * math.comb(23, 10) == 41186376, '%d' % (3 * 12 * math.comb(23, 10)))
    check('the_printed_row_sums_of_the_per_target_tables_are_arithmetically_right',
          sum(CLAIM_T12['#1395']) == 1980 and sum(CLAIM_T12['#44170']) == 734
          and sum(CLAIM_T12['#44172']) == 542 and sum(CLAIM_T13['#44170']) == 52
          and sum(CLAIM_T13['#44172']) == 12 and sum(CLAIM_T13['#1395']) == 24,
          '1980 / 734 / 542 and 24 / 52 / 12')
    check('each_of_the_three_refutation_witnesses_has_exactly_thirteen_pebbles',
          all(sum(c[3]) == 13 for c in CERTS[:3]),
          str([sum(c[3]) for c in CERTS[:3]]))


def _run(tasks, ctx):
    """Run the sweep tasks, in parallel when a pool is available."""
    t0 = time.time()
    jobs = ctx.get('jobs')
    if jobs and jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs) as pool:
            res = pool.map(_task, tasks, chunksize=1)
    else:
        res = [_task(t) for t in tasks]
    note('%d sweeps in %.1f s wall (%d worker%s), slowest single sweep %.1f s'
         % (len(tasks), time.time() - t0, jobs or 1, '' if (jobs or 1) == 1 else 's',
            max(r[4] for r in res)))
    return res


STEPS_FAST = [step1_objects, step2_invariants, step3_identify, step4_certificates,
              step5_lemmas, step6b_census_structure, step7_arithmetic]
STEPS_FULL = [step1_objects, step2_invariants, step3_identify, step4_certificates,
              step5_lemmas, step6a_sweep, step6b_census_structure, step6c_census_pi,
              step7_arithmetic]


def verdict():
    n = len(RESULTS)
    bad = [nm for nm, ok in RESULTS if not ok]
    print('')
    if bad:
        print('VERDICT: %d OF %d CHECKS FAILED' % (len(bad), n))
        for nm in bad:
            print('  failed: %s' % nm)
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % n)
    return 0


def main(argv):
    quick = '--quick' in argv
    jobs = 0
    for a in argv:
        if a.startswith('-j'):
            jobs = int(a[2:])
    if not jobs:
        try:
            jobs = min(64, os.cpu_count() or 1)
        except Exception:                                   # noqa: BLE001
            jobs = 1
    ctx = {'jobs': 1 if quick else jobs}
    steps = STEPS_FAST if quick else STEPS_FULL
    print('verification of the note: pebbling numbers of four cubic 12-vertex graphs')
    print('of girth 3 and diameter 3 -- #6698 is Class 0, and')
    print('pi(#1395) = pi(#44170) = pi(#44172) = 14, not the suggested 13')
    print('python %s, exact integer arithmetic only' % sys.version.split()[0])
    print('mode: %s, workers: %d' % ('QUICK (exhaustive sweep skipped)' if quick else 'FULL',
                                     ctx['jobs']))
    for step in steps:
        try:
            step(ctx)
        except SystemExit:
            raise
        except Exception:                                    # noqa: BLE001
            import traceback
            traceback.print_exc()
            check('step_%s_completed_without_an_exception' % step.__name__, False,
                  'see traceback above')
    section('scope of this run')
    if quick:
        note('SCOPE: --quick was used, so the exhaustive sweep did NOT run: the upper bounds '
             'pi <= 14 and pi(#6698) <= 12 are NOT verified by this transcript.')
    note('SCOPE: this program recomputes, from the adjacency lists and graph6 strings transcribed '
         'into it, each quantity listed in the CLAIM_* and CENSUS tables above (invariants, '
         'orbits, certificate reachable-set sizes, per-target unsolvable counts at totals 12 and '
         '13, the pebbling numbers of all five named graphs and of all sixteen transcribed '
         'graph6 graphs). '
         'The paper itself is not parsed, so the transcription of its printed objects and figures '
         'into the embedded inputs, and the coverage of every quantity the paper asserts, are not '
         'checked here.')
    note('NOT RE-RUN: any count of how many cubic 12-vertex girth-3 diameter-3 graphs there are. '
         'This run re-derives only that the sixteen graph6 strings transcribed into this file are '
         'sixteen pairwise non-isomorphic such graphs, i.e. at least sixteen; no completeness '
         'statement is verified here, and the note claims none.')
    note('NOT RE-RUN: what the House of Graphs database holds, and the identification of the four '
         'adjacency lists with the database ids #1395, #6698, #44170 and #44172 -- the lists are '
         'inputs, printed in the paper, and the paper pins them to the source figure by an edge-set '
         'identity rather than to the database.')
    return verdict()


if __name__ == '__main__':
    try:
        code = main(sys.argv[1:])
    except SystemExit:
        raise
    except Exception:                                        # noqa: BLE001
        import traceback
        traceback.print_exc()
        check('verifier_ran_to_completion_without_an_exception', False)
        code = verdict()
    sys.exit(code)

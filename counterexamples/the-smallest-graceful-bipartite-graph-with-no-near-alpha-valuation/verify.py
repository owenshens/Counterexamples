#!/usr/bin/env python3
"""verify.py -- re-derives every quantity claimed in paper.tex, from the objects the paper prints.

    python3 verify.py            # Python 3.9+, standard library only, no external data file

WHAT IT READS. Everything the program checks against is transcribed below from the paper itself:
the edge list of G = C_4 u P_3 (paper, Section 2), the beta-valuation exhibited there, the 15-row
case table of the proof of Theorem 1, the two further witnesses of Section 3 with their printed
labellings, and the census table of Section 3. No value is imported from the run that produced the
paper; each is recomputed here and compared.

ARITHMETIC. Every quantity is an integer. There is no floating point anywhere in this file and no
decision is taken on an inexact value.

DEFINITIONS USED, verbatim from the source paper (arXiv:2603.05662, lines 104 and 317-326).
  * A beta-valuation (graceful labelling) of a graph with m edges is an INJECTIVE map
    b : V -> {0, 1, ..., m} whose induced edge labels |b(u) - b(v)| are exactly {1, ..., m}.
  * A near alpha-valuation is a beta-valuation for which V admits a partition into V_small and
    V_large such that every vertex of V_small carries a label below all of its neighbours' and
    every vertex of V_large a label above all of its neighbours'.
  * An alpha-valuation is a beta-valuation for which some integer x satisfies
    min(b(u), b(v)) <= x < max(b(u), b(v)) for every edge uv.

WHAT IT DOES NOT COVER is printed by the program itself, in the NOT RE-RUN block before the verdict,
and is repeated in REVIEW_NOTE.md's `## Scope`.
"""

import itertools
import math
import sys

# =====================================================================================
# 0. THE OBJECTS AS PRINTED IN THE PAPER
# =====================================================================================

# Section 2. G = C_4 u P_3 on v1,...,v7, stored 0-based: vertex i is v_{i+1}.
G1_N = 7
G1_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6)]
G1_PARTS = ([0, 2, 4, 6], [1, 3, 5])                       # {v1,v3,v5,v7} and {v2,v4,v6}
# The beta-valuation exhibited in Section 2: b(v1..v4) = 0,6,2,5 on the cycle, b(v5..v7) = 1,3,4.
G1_LABELLING = (0, 6, 2, 5, 1, 3, 4)
G1_EDGE_LABELS = [6, 4, 3, 5, 2, 1]                        # in the order of G1_EDGES

# Section 3, the two further witnesses.
G2_N = 11
G2_EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),          # C_6
            (6, 7), (6, 8), (6, 9), (6, 10)]                          # K_{1,4}
# The labelling printed for C_6 u K_{1,4}: cycle 0-9-5-3-2-10-0, star centre 1 with leaves 4,6,7,8.
G2_LABELLING = (0, 9, 5, 3, 2, 10, 1, 4, 6, 7, 8)
G2_MIDDLE = (3, 2, 5)                                      # label 3 lies strictly between 2 and 5

G3_N = 12
G3_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0),                           # C_4
            (4, 5), (5, 6), (6, 7), (7, 4),                           # C_4
            (8, 9), (8, 10), (8, 11)]                                 # K_{1,3}
# The labelling printed for 2C_4 u K_{1,3}: cycles 0-10-2-11-0 and 1-7-3-8-1, star centre 6 (4,5,9).
G3_LABELLING = (0, 10, 2, 11, 1, 7, 3, 8, 6, 4, 5, 9)
G3_MIDDLE = (6, 5, 9)                                      # label 6 lies strictly between 5 and 9

# Section 2, the 15-row case table of the proof of Theorem 1. Each row is
#   (case, cycle labels, cycle edge labels D, path labels, required path edge labels, reason)
# with reason codes  R = two cycle edge labels coincide,
#                    N = no choice of path centre realises the required pair,
#                    X = the only centre realising it is not a local extremum.
PAPER_TABLE = [
    ('A', (0, 1, 2, 6), (1, 2, 5, 6), (3, 4, 5), (3, 4), 'N'),
    ('A', (0, 1, 3, 6), (2, 3, 5, 6), (2, 4, 5), (1, 4), 'N'),
    ('A', (0, 1, 4, 6), (3, 4, 5, 6), (2, 3, 5), (1, 2), 'X'),
    ('A', (0, 1, 5, 6), (4, 5, 5, 6), (2, 3, 4), (),     'R'),
    ('A', (0, 2, 3, 6), (1, 3, 4, 6), (1, 4, 5), (2, 5), 'N'),
    ('A', (0, 2, 4, 6), (2, 4, 4, 6), (1, 3, 5), (),     'R'),
    ('A', (0, 2, 5, 6), (3, 4, 5, 6), (1, 3, 4), (1, 2), 'X'),
    ('A', (0, 3, 4, 6), (1, 3, 4, 6), (1, 2, 5), (2, 5), 'N'),
    ('A', (0, 3, 5, 6), (2, 3, 5, 6), (1, 2, 4), (1, 4), 'N'),
    ('A', (0, 4, 5, 6), (1, 2, 5, 6), (1, 2, 3), (3, 4), 'N'),
    ('B', (1, 2, 3, 4), (1, 2, 2, 3), (0, 5, 6), (),     'R'),
    ('B', (1, 2, 3, 5), (1, 2, 3, 4), (0, 4, 6), (5, 6), 'N'),
    ('B', (1, 2, 4, 5), (2, 3, 3, 4), (0, 3, 6), (),     'R'),
    ('B', (1, 3, 4, 5), (1, 2, 3, 4), (0, 2, 6), (5, 6), 'N'),
    ('B', (2, 3, 4, 5), (1, 2, 2, 3), (0, 1, 6), (),     'R'),
]

# Section 3, the census table. m -> (isomorphism classes of graceful graphs,
#                                   of which bipartite,
#                                   of which with no near alpha-valuation).
# The m = 11 row is printed in the paper and is NOT re-run by this program (see NOT RE-RUN below).
PAPER_CENSUS = {
    1: (1, 1, 0), 2: (1, 1, 0), 3: (3, 2, 0), 4: (5, 4, 0), 5: (12, 7, 0),
    6: (37, 18, 1), 7: (112, 44, 0), 8: (340, 119, 0), 9: (1078, 323, 0), 10: (3620, 919, 1),
}
CENSUS_M = 10                                  # the range this program re-derives
PAPER_A006967 = [1, 2, 2, 4, 12, 16, 20, 60, 148, 324]     # A006967(m+1)/2, m = 1..10
PAPER_A000055 = [1, 1, 2, 3, 6, 11, 23, 47, 106, 235]      # A000055(m+1),   m = 1..10

# Section 3 / Section 5, the named controls. name -> (n, edges, graceful, near-alpha, alpha)
PAPER_CONTROLS = [
    ('K_{2,3}', 5, [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)], 48, 48, 48),
    ('C_8', 8, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0)], 384, 96, 96),
    ('S_{3,2}', 7, [(0, 1), (1, 2), (0, 3), (3, 4), (0, 5), (5, 6)], 60, 24, 0),
    ('C_3', 3, [(0, 1), (1, 2), (2, 0)], 12, 0, 0),
    ('C_5', 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)], 0, 0, 0),
    ('C_6', 6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)], 0, 0, 0),
    ('C_7', 7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0)], 168, 0, 0),
    ('C_4 u K_2', 6, [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5)], 0, 0, 0),
    ('P_7', 7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)], 32, 16, 16),
    ('2C_4 u K_2', 10, [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (8, 9)],
     512, 512, 512),
]

# =====================================================================================
# 1. THE CHECK HARNESS
# =====================================================================================

_passed = 0
_failed = 0


def ok(name, detail=''):
    global _passed
    _passed += 1
    print('PASS %s%s' % (name, ('  ' + detail) if detail else ''))
    sys.stdout.flush()


def check(name, condition, detail=''):
    global _failed
    if condition:
        ok(name, detail)
    else:
        _failed += 1
        print('FAIL %s  %s' % (name, detail))
        sys.stdout.flush()


def eq(name, got, want, extra=''):
    check(name, got == want, 'got %r, paper says %r%s' % (got, want, ('; ' + extra) if extra else ''))


# =====================================================================================
# 2. GRAPH PRIMITIVES
# =====================================================================================

def adjacency(n, edges):
    nbrs = [[] for _ in range(n)]
    for u, v in edges:
        nbrs[u].append(v)
        nbrs[v].append(u)
    return nbrs


def components(n, nbrs):
    seen = [False] * n
    out = []
    for s in range(n):
        if seen[s]:
            continue
        seen[s] = True
        stack = [s]
        comp = []
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in nbrs[x]:
                if not seen[y]:
                    seen[y] = True
                    stack.append(y)
        out.append(sorted(comp))
    return out


def two_colouring(n, nbrs):
    """(is_bipartite, colour list). Colour -1 on isolated vertices cannot occur here."""
    col = [-1] * n
    bip = True
    for s in range(n):
        if col[s] != -1:
            continue
        col[s] = 0
        stack = [s]
        while stack:
            x = stack.pop()
            for y in nbrs[x]:
                if col[y] == -1:
                    col[y] = 1 - col[x]
                    stack.append(y)
                elif col[y] == col[x]:
                    bip = False
    return bip, col


def is_graceful(n, edges, lab):
    m = len(edges)
    if len(set(lab)) != n or min(lab) < 0 or max(lab) > m:
        return False
    return sorted(abs(lab[u] - lab[v]) for u, v in edges) == list(range(1, m + 1))


def is_near_alpha(nbrs, lab):
    """The paper's definition. With no isolated vertex, the V_small/V_large partition exists exactly
    when every vertex is strictly below all its neighbours or strictly above all of them: take
    V_small to be the local minima. A vertex could be both only if it had no neighbour."""
    for v, nb in enumerate(nbrs):
        if not nb:
            return False
        if not (all(lab[w] > lab[v] for w in nb) or all(lab[w] < lab[v] for w in nb)):
            return False
    return True


def is_alpha(edges, lab):
    lo = max(min(lab[u], lab[v]) for u, v in edges)
    hi = min(max(lab[u], lab[v]) for u, v in edges)
    return lo < hi


def order_for_search(n, nbrs):
    """Vertex order in which each vertex after the first of its component has at least one
    already-ordered neighbour: every edge then fixes its difference exactly once."""
    order = []
    placed = [False] * n
    while len(order) < n:
        best, key = None, None
        for v in range(n):
            if placed[v]:
                continue
            k = (sum(1 for u in nbrs[v] if placed[u]), len(nbrs[v]))
            if key is None or k > key:
                key, best = k, v
        order.append(best)
        placed[best] = True
    return order


def all_graceful_labellings(n, edges):
    """Every beta-valuation of the given graph, as a tuple indexed by vertex. Exhaustive: the search
    tries every label for every vertex and prunes only on the two conditions in the definition
    (labels distinct, edge differences distinct)."""
    m = len(edges)
    nbrs = adjacency(n, edges)
    order = order_for_search(n, nbrs)
    seen = set()
    back = []
    for v in order:
        back.append([u for u in nbrs[v] if u in seen])
        seen.add(v)
    lab = [-1] * n
    usedlab = [False] * (m + 1)
    useddiff = [False] * (m + 1)
    out = []

    def bt(i):
        if i == n:
            out.append(tuple(lab))
            return
        v = order[i]
        for L in range(m + 1):
            if usedlab[L]:
                continue
            ds = []
            good = True
            for u in back[i]:
                d = abs(L - lab[u])
                if useddiff[d] or d in ds:
                    good = False
                    break
                ds.append(d)
            if not good:
                continue
            lab[v] = L
            usedlab[L] = True
            for d in ds:
                useddiff[d] = True
            bt(i + 1)
            for d in ds:
                useddiff[d] = False
            usedlab[L] = False
            lab[v] = -1

    bt(0)
    return out


def automorphism_count(n, edges):
    """|Aut(G)|, by exhaustive backtracking over adjacency-preserving bijections."""
    nbrs = adjacency(n, edges)
    adm = [0] * n
    for v in range(n):
        for w in nbrs[v]:
            adm[v] |= 1 << w
    deg = [len(nbrs[v]) for v in range(n)]
    order = sorted(range(n), key=lambda v: -deg[v])
    mp = [-1] * n
    total = 0

    def bt(i, img, used):
        nonlocal total
        if i == n:
            total += 1
            return
        v = order[i]
        want = 0
        for u in nbrs[v]:
            if mp[u] >= 0:
                want |= 1 << mp[u]
        for w in range(n):
            b = 1 << w
            if used & b or deg[w] != deg[v]:
                continue
            if adm[w] & img != want:
                continue
            mp[v] = w
            bt(i + 1, img | b, used | b)
            mp[v] = -1

    bt(0, 0, 0)
    return total


def refine(n, nbrs):
    """An isomorphism-invariant vertex colouring, used only to bucket and to restrict candidates in
    the exact isomorphism test below. It is never used as a substitute for that test."""
    inv = [len(nbrs[v]) for v in range(n)]
    for _ in range(n):
        new = [(inv[v],) + tuple(sorted(inv[u] for u in nbrs[v])) for v in range(n)]
        idx = {x: i for i, x in enumerate(sorted(set(new)))}
        nxt = [idx[x] for x in new]
        if nxt == inv:
            break
        inv = nxt
    return inv


def isomorphic(n, nbrsA, invA, admB, invB):
    """Exact isomorphism test between two graphs on the same vertex count."""
    order = sorted(range(n), key=lambda v: (-len(nbrsA[v]), invA[v]))
    cand = [[w for w in range(n) if invB[w] == invA[v]] for v in order]
    mp = [-1] * n

    def bt(i, img, used):
        if i == n:
            return True
        v = order[i]
        want = 0
        for u in nbrsA[v]:
            if mp[u] >= 0:
                want |= 1 << mp[u]
        for w in cand[i]:
            b = 1 << w
            if used & b or admB[w] & img != want:
                continue
            mp[v] = w
            if bt(i + 1, img | b, used | b):
                return True
            mp[v] = -1
        return False

    return bt(0, 0, 0)


def graph_record(n, edges):
    nbrs = adjacency(n, edges)
    adm = [0] * n
    for v in range(n):
        for w in nbrs[v]:
            adm[v] |= 1 << w
    return nbrs, adm, refine(n, nbrs)


# =====================================================================================
# 3. THE OBJECT OF SECTION 2
# =====================================================================================

print('=== 1. THE OBJECT PRINTED IN SECTION 2: G = C_4 u P_3 ===')
g1_nbrs = adjacency(G1_N, G1_EDGES)
m1 = len(G1_EDGES)
eq('object-order-and-size', (G1_N, m1), (7, 6))
check('object-vertices-is-edges-plus-one', G1_N == m1 + 1,
      'n = %d, m = %d, so every beta-valuation of G is a bijection onto {0,...,%d}' % (G1_N, m1, m1))
check('object-simple', all(u != v for u, v in G1_EDGES)
      and len({frozenset(e) for e in G1_EDGES}) == m1, 'no loop, no repeated edge')
check('object-no-isolated-vertex', all(g1_nbrs[v] for v in range(G1_N)),
      'minimum degree %d' % min(len(x) for x in g1_nbrs))
_bip, _col = two_colouring(G1_N, g1_nbrs)
check('object-bipartite', _bip and (sorted(v for v in range(G1_N) if _col[v] == _col[0])
                                    in (sorted(G1_PARTS[0]), sorted(G1_PARTS[1]))),
      'parts {v1,v3,v5,v7} and {v2,v4,v6} as printed')
_comps = components(G1_N, g1_nbrs)
eq('object-two-components', [len(c) for c in _comps], [4, 3], 'so G is DISCONNECTED')
eq('object-degree-sequence', sorted(len(x) for x in g1_nbrs), [1, 1, 2, 2, 2, 2, 2],
   'a 4-cycle and a path on three vertices')
check('object-cycle-is-C4', all(len(g1_nbrs[v]) == 2 for v in _comps[0]) and len(_comps[0]) == 4,
      'component {v1,v2,v3,v4} is 2-regular on 4 vertices')

print()
print('=== 2. THE BETA-VALUATION EXHIBITED IN SECTION 2 ===')
eq('exhibited-labels-injective', len(set(G1_LABELLING)), 7)
eq('exhibited-label-set', sorted(G1_LABELLING), list(range(7)))
eq('exhibited-edge-labels', [abs(G1_LABELLING[u] - G1_LABELLING[v]) for u, v in G1_EDGES],
   G1_EDGE_LABELS, 'in the order the paper lists the edges')
check('exhibited-is-graceful', is_graceful(G1_N, G1_EDGES, G1_LABELLING),
      'edge labels are exactly {1,...,6}, so G is graceful')
check('exhibited-is-not-near-alpha', not is_near_alpha(g1_nbrs, G1_LABELLING),
      'v6 carries 3 with neighbours 1 and 4, so it is neither a local minimum nor a local maximum')

# =====================================================================================
# 4. THE PROOF OF THEOREM 1, RE-DERIVED
# =====================================================================================

print()
print('=== 3. THE PROOF OF THEOREM 1: THE 15-CASE TABLE ===')

# The equivalence the proof uses in place of the partition.
_bad = []
for lab in itertools.permutations(range(7), 7):
    partition_exists = False
    for mask in range(1 << 7):
        small = [v for v in range(7) if mask >> v & 1]
        large = [v for v in range(7) if not mask >> v & 1]
        if all(all(lab[w] > lab[v] for w in g1_nbrs[v]) for v in small) and \
           all(all(lab[w] < lab[v] for w in g1_nbrs[v]) for v in large):
            partition_exists = True
            break
    if partition_exists != is_near_alpha(g1_nbrs, lab):
        _bad.append(lab)
eq('near-alpha-equals-every-vertex-extremal', len(_bad), 0,
   'over all 5040 labellings of G, the V_small/V_large partition exists exactly when every vertex '
   'is a local minimum or a local maximum')

# The C_4 identity.
_ident = [(p, q, r, s) for p, q, r, s in itertools.combinations(range(7), 4)
          if (s - p) + (r - q) != (r - p) + (s - q)]
eq('c4-identity', len(_ident), 0,
   '(s-p)+(r-q) = (r-p)+(s-q) for every p<q<r<s in {0,...,6}: %d quadruples checked'
   % len(list(itertools.combinations(range(7), 4))))
_extreme = [(p, q, r, s) for p, q, r, s in itertools.combinations(range(7), 4)
            if not (s - p == max(s - p, r - p, s - q, r - q)
                    and r - q == min(s - p, r - p, s - q, r - q))]
eq('c4-identity-extremes', len(_extreme), 0,
   's-p is the largest and r-q the smallest of the four cycle edge labels')

# On a 4-cycle, the minima must be a bipartition class.
_splits = []
cyc = [(0, 1), (1, 2), (2, 3), (3, 0)]
for mask in range(1 << 4):
    S = {v for v in range(4) if mask >> v & 1}
    if all((u in S) != (v in S) for u, v in cyc):
        _splits.append(tuple(sorted(S)))
eq('c4-forced-split', sorted(_splits), [(0, 2), (1, 3)],
   'the only ways to split the four cycle vertices into minima and maxima are the two bipartition '
   'classes, so the cycle labels are a low pair {p<q} and a high pair {r<s} with q < r')


def path_rows(cycle_labels):
    """For a 4-set of cycle labels split low/high, the paper's row: cycle edge labels D, the three
    remaining labels, the required pair of path edge labels, and the reason the row dies."""
    p, q, r, s = cycle_labels
    D = sorted([s - p, r - p, s - q, r - q])
    path = tuple(sorted(set(range(7)) - set(cycle_labels)))
    if len(set(D)) < 4:
        return tuple(D), path, (), 'R', []
    need = tuple(sorted(set(range(1, 7)) - set(D)))
    hits = []
    for c in path:
        other = [x for x in path if x != c]
        labels = tuple(sorted([abs(c - other[0]), abs(c - other[1])]))
        if labels == need:
            extremal = (c < other[0] and c < other[1]) or (c > other[0] and c > other[1])
            hits.append((c, extremal))
    if not hits:
        return tuple(D), path, need, 'N', hits
    if any(e for _, e in hits):
        return tuple(D), path, need, 'SURVIVES', hits
    return tuple(D), path, need, 'X', hits


rows = []
for q, r in itertools.combinations(range(1, 6), 2):          # CASE A: 0 and 6 in the cycle
    cyc4 = (0, q, r, 6)
    D, path, need, why, _hits = path_rows(cyc4)
    rows.append(('A', cyc4, D, path, need, why))
for L in itertools.combinations(range(1, 6), 4):             # CASE B: 0 and 6 in the path
    D, path, need, why, _hits = path_rows(L)
    rows.append(('B', tuple(L), D, path, need, why))

eq('case-table-row-count', len(rows), 15, '10 rows in case A and 5 in case B')
eq('case-table-case-A', sum(1 for x in rows if x[0] == 'A'), 10,
   'the pairs 0<q<r<6 giving cycle labels {0,q,r,6}')
eq('case-table-case-B', sum(1 for x in rows if x[0] == 'B'), 5,
   'the 4-subsets of {1,...,5} giving the cycle labels')
check('case-table-matches-paper', rows == PAPER_TABLE,
      '%d of the %d printed rows disagree; each row carries the cycle labels, the cycle edge labels, '
      'the path labels, the required path edge labels and the reason, and all of them are '
      're-derived here' % (sum(1 for a, b in zip(rows, PAPER_TABLE) if a != b), len(PAPER_TABLE)))
eq('case-table-no-survivor', [x for x in rows if x[5] == 'SURVIVES'], [],
   'no row admits a near alpha-valuation, which proves Theorem 1')
eq('case-table-reason-tally',
   tuple(sum(1 for x in rows if x[5] == c) for c in 'RNX'), (5, 8, 2),
   '5 rows die because two cycle edge labels coincide, 8 because no path centre realises the '
   'required pair, 2 because the only centre that realises it is not a local extremum')

# =====================================================================================
# 5. BRUTE FORCE OVER THE WHOLE LABELLING SPACE OF G, INDEPENDENT OF SECTION 4
# =====================================================================================

print()
print('=== 4. BRUTE FORCE OVER ALL 5040 LABELLINGS OF G (independent of the case table) ===')
_tot = _g = _na = _al = 0
_gs = []
for lab in itertools.permutations(range(7), 7):
    _tot += 1
    if not is_graceful(G1_N, G1_EDGES, lab):
        continue
    _g += 1
    _gs.append(lab)
    if is_near_alpha(g1_nbrs, lab):
        _na += 1
    if is_alpha(G1_EDGES, lab):
        _al += 1
eq('brute-labellings-examined', _tot, 5040, '7! bijections of V onto {0,...,6}')
eq('brute-graceful', _g, 32)
eq('brute-near-alpha', _na, 0, 'THE REFUTATION: no beta-valuation of G is a near alpha-valuation')
eq('brute-alpha', _al, 0, 'and none is an alpha-valuation either')
check('brute-contains-exhibited', tuple(G1_LABELLING) in _gs,
      'the labelling printed in Section 2 is among the 32 found')
_aut1 = automorphism_count(G1_N, G1_EDGES)
eq('brute-automorphisms', _aut1, 16, '|Aut(C_4)| * |Aut(P_3)| = 8 * 2')
eq('brute-orbits', _g // _aut1, 2, '32 labellings in 2 orbits, so one graph and not several')

# =====================================================================================
# 6. THE TWO FURTHER WITNESSES OF SECTION 3
# =====================================================================================

print()
print('=== 5. THE TWO FURTHER WITNESSES OF SECTION 3 ===')
for tag, n, edges, lab, middle, want_nm, want_g, want_aut, want_orb, shape in (
        ('C_6uK_{1,4}', G2_N, G2_EDGES, G2_LABELLING, G2_MIDDLE, (11, 10), 576, 288, 2, [5, 6]),
        ('2C_4uK_{1,3}', G3_N, G3_EDGES, G3_LABELLING, G3_MIDDLE, (12, 11), 4608, 768, 6,
         [4, 4, 4])):
    nb = adjacency(n, edges)
    m = len(edges)
    eq('%s-order-and-size' % tag, (n, m), want_nm, 'vertices and edges, as printed')
    check('%s-vertices-is-edges-plus-one' % tag, n == m + 1, 'n = %d = m + 1' % n)
    _b, _c = two_colouring(n, nb)
    check('%s-bipartite' % tag, _b, 'two-colourable')
    _cm = components(n, nb)
    eq('%s-components' % tag, sorted(len(c) for c in _cm), sorted(shape),
       'DISCONNECTED: %d components' % len(_cm))
    check('%s-printed-labelling-graceful' % tag, is_graceful(n, edges, lab),
          'edge labels are exactly {1,...,%d}' % m)
    _mid, _lo, _hi = middle
    _vs = [v for v in range(n) if lab[v] == _mid]
    _nbl = sorted(lab[w] for w in nb[_vs[0]])
    check('%s-printed-middle-vertex' % tag, _lo in _nbl and _hi in _nbl and _lo < _mid < _hi,
          'the vertex labelled %d has neighbours %s, and %d < %d < %d, so it is neither a local '
          'minimum nor a local maximum' % (_mid, _nbl, _lo, _mid, _hi))
    labs = all_graceful_labellings(n, edges)
    eq('%s-graceful-count' % tag, len(labs), want_g, 'every beta-valuation enumerated')
    check('%s-printed-labelling-found' % tag, tuple(lab) in set(labs),
          'the printed labelling is one of the %d' % want_g)
    eq('%s-near-alpha-count' % tag, sum(1 for L in labs if is_near_alpha(nb, L)), 0,
       'A COUNTEREXAMPLE: no beta-valuation of this graph is a near alpha-valuation')
    _aut = automorphism_count(n, edges)
    eq('%s-automorphisms' % tag, _aut, want_aut)
    eq('%s-orbits' % tag, len(labs) // _aut, want_orb,
       '%d labellings in %d orbits' % (want_g, want_orb))

# =====================================================================================
# 7. THE CENSUS
# =====================================================================================

print()
print('=== 6. THE CENSUS: every graph with at most %d edges, up to isomorphism ===' % CENSUS_M)
print('    (enumerating LABELLINGS, which is legitimate because a graceful graph with m edges has')
print('     at most m+1 vertices: for each difference d in {1,...,m} the edge realising it is one of')
print('     the m-d+1 pairs {a, a+d}, so there are exactly m! graceful labelled edge sets.)')
sys.stdout.flush()

buckets = {}
census_rows = {}
census_paths = []
census_trees = []
total_labellings = 0
total_processed = 0
for m in range(1, CENSUS_M + 1):
    ds = list(range(1, m + 1))
    ranges = [range(0, m - d + 1) for d in ds]
    npaths = 0
    for choice in itertools.product(*ranges):
        total_labellings += 1
        nb = {}
        for d, a in zip(ds, choice):
            nb.setdefault(a, []).append(a + d)
            nb.setdefault(a + d, []).append(a)
        n = len(nb)
        if n == m + 1 and sorted(len(x) for x in nb.values()) == [1, 1] + [2] * (m - 1):
            start = next(iter(nb))
            seen = {start}
            stack = [start]
            while stack:
                x = stack.pop()
                for y in nb[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            if len(seen) == n:
                npaths += 1
        # b -> m - b maps graceful labellings to graceful labellings of an isomorphic graph and
        # preserves "every vertex is a local extremum", so one of each pair suffices here.
        if tuple(m - d - a for d, a in zip(ds, choice)) < choice:
            continue
        total_processed += 1
        near = True
        for v, lst in nb.items():
            lo = hi = True
            for w in lst:
                if w > v:
                    lo = False
                else:
                    hi = False
            if not (lo or hi):
                near = False
                break
        verts = sorted(nb)
        pos = {x: i for i, x in enumerate(verts)}
        nbrs = [[] for _ in range(n)]
        adm = [0] * n
        for v in verts:
            i = pos[v]
            for w in nb[v]:
                j = pos[w]
                nbrs[i].append(j)
                adm[i] |= 1 << j
        inv = refine(n, nbrs)
        key = (n, tuple(sorted(len(x) for x in nbrs)), tuple(sorted(inv)))
        lst = buckets.setdefault(key, [])
        rec = None
        for cand in lst:
            if isomorphic(n, nbrs, inv, cand['adm'], cand['inv']):
                rec = cand
                break
        if rec is None:
            bip, _ = two_colouring(n, nbrs)
            ncomp = len(components(n, nbrs))
            rec = {'adm': adm, 'inv': inv, 'n': n, 'm': m, 'near': False, 'bip': bip,
                   'conn': ncomp == 1, 'tree': ncomp == 1 and n == m + 1,
                   'edges': tuple(sorted((verts[i], verts[j])
                                         for i in range(n) for j in nbrs[i] if i < j))}
            lst.append(rec)
        if near:
            rec['near'] = True
    recs = [r for b in buckets.values() for r in b if r['m'] == m]
    bip = [r for r in recs if r['bip']]
    ref = [r for r in bip if not r['near']]
    trees = [r for r in recs if r['tree']]
    census_rows[m] = (len(recs), len(bip), len(ref))
    census_paths.append(npaths)
    census_trees.append(len(trees))
    print('    m=%-3d classes=%-6d bipartite=%-5d connected-bipartite=%-5d no-near-alpha=%d'
          % (m, len(recs), len(bip), sum(1 for r in bip if r['conn']), len(ref)))
    sys.stdout.flush()
    for r in ref:
        print('       COUNTEREXAMPLE at m=%d: %d vertices, connected=%s, labelled edges %s'
              % (m, r['n'], r['conn'], list(r['edges'])))
        sys.stdout.flush()

allrecs = [r for b in buckets.values() for r in b]
bipall = [r for r in allrecs if r['bip']]
refall = [r for r in bipall if not r['near']]
eq('census-labellings-enumerated', total_labellings,
   sum(math.factorial(m) for m in range(1, CENSUS_M + 1)),
   'the sum of m! for m = 1..%d; %d of them classified after the b -> m-b reduction'
   % (CENSUS_M, total_processed))
for m in range(1, CENSUS_M + 1):
    eq('census-row-m%d' % m, census_rows[m], PAPER_CENSUS[m],
       '(classes, bipartite, with no near alpha-valuation)')
eq('census-counterexample-edge-counts', sorted(r['m'] for r in refall), [6, 10],
   'these are the only edge counts at most %d at which a bipartite graceful graph has no near '
   'alpha-valuation' % CENSUS_M)
eq('census-counterexample-count', len(refall), 2, 'one at m=6 and one at m=10, up to isomorphism')
_r6 = [r for r in refall if r['m'] == 6][0]
_r10 = [r for r in refall if r['m'] == 10][0]
_n1, _a1, _i1 = graph_record(G1_N, G1_EDGES)
check('census-m6-counterexample-is-G', isomorphic(7, _n1, _i1, _r6['adm'], _r6['inv']),
      'the unique 6-edge counterexample is C_4 u P_3, which is the minimality and uniqueness '
      'statement of Section 3 of the paper')
_n2, _a2, _i2 = graph_record(G2_N, G2_EDGES)
check('census-m10-counterexample-is-C6uK14', isomorphic(11, _n2, _i2, _r10['adm'], _r10['inv']),
      'the unique 10-edge counterexample is C_6 u K_{1,4}')
eq('census-none-below-6-edges', [r['m'] for r in refall if r['m'] < 6], [],
   'MINIMALITY: every bipartite graceful graph with at most 5 edges has a near alpha-valuation '
   '(%d isomorphism classes)' % len([r for r in bipall if r['m'] <= 5]))
eq('census-all-counterexamples-disconnected', [r['m'] for r in refall if r['conn']], [],
   'THE SCOPE LIMIT: both counterexamples in range are disconnected')
eq('census-connected-counterexamples', sum(1 for r in bipall if r['conn'] and not r['near']), 0,
   'all %d connected bipartite graceful classes with at most %d edges HAVE a near alpha-valuation'
   % (sum(1 for r in bipall if r['conn']), CENSUS_M))
_conn = sum(1 for r in bipall if r['conn'])
_disc = sum(1 for r in bipall if not r['conn'])
_nonbip = sum(1 for r in allrecs if not r['bip'])
eq('census-crossfoot', (_conn + _disc, len(bipall) + _nonbip, sum(census_rows[m][0]
                                                                 for m in census_rows)),
   (len(bipall), len(allrecs), len(allrecs)),
   'connected %d + disconnected %d = %d bipartite classes; bipartite + non-bipartite = %d classes, '
   'which is also the sum of the per-m rows' % (_conn, _disc, len(bipall), len(allrecs)))
eq('control-A006967-paths', census_paths, PAPER_A006967,
   'graceful labelled edge sets of the path on m+1 vertices, against A006967(m+1)/2')
eq('control-A000055-trees', census_trees, PAPER_A000055,
   'graceful tree classes with m edges -- every tree on at most %d vertices is graceful -- against '
   'A000055(m+1)' % (CENSUS_M + 1))
eq('control-every-tree-has-near-alpha',
   sum(1 for r in allrecs if r['tree'] and not r['near']), 0,
   'all %d tree classes in range have a near alpha-valuation, as El-Zanati-Kenig-Vanden Eynden and '
   'Grannell-Griggs-Holroyd conjecture in general' % sum(1 for r in allrecs if r['tree']))
eq('control-no-non-bipartite-has-near-alpha',
   sum(1 for r in allrecs if not r['bip'] and r['near']), 0,
   'none of the %d non-bipartite classes has one, as the source paper proves'
   % sum(1 for r in allrecs if not r['bip']))
_cyc = sorted(r['m'] for r in allrecs
              if r['conn'] and r['n'] == r['m']
              and all(bin(a).count('1') == 2 for a in r['adm']))
eq('control-cycles-graceful-iff-0-or-3-mod-4', _cyc,
   [m for m in range(3, CENSUS_M + 1) if m % 4 in (0, 3)],
   "Rosa's theorem: C_m is graceful exactly for m = 0 or 3 mod 4, and the census finds a graceful "
   "2-regular connected class at exactly those m")

# =====================================================================================
# 8. NAMED CONTROLS
# =====================================================================================

print()
print('=== 7. NAMED CONTROLS, both polarities ===')
for name, n, edges, wg, wna, wal in PAPER_CONTROLS:
    nb = adjacency(n, edges)
    labs = all_graceful_labellings(n, edges)
    g = len(labs)
    na = sum(1 for L in labs if is_near_alpha(nb, L))
    al = sum(1 for L in labs if is_alpha(edges, L))
    eq('control-%s' % name.replace(' ', ''), (g, na, al), (wg, wna, wal),
       '(beta-valuations, of which near alpha, of which alpha)')

# =====================================================================================
# 9. SCOPE AND VERDICT
# =====================================================================================

print()
print('NOT RE-RUN: the census at m = 11. This program re-derives the classification for m <= %d '
      'only, so the paper\'s statement that C_4 u P_3, C_6 u K_{1,4} and 2C_4 u K_{1,3} are the '
      'ONLY counterexamples with at most 11 edges is verified here only up to 10 edges. The 11-edge '
      'witness 2C_4 u K_{1,3} is itself fully verified above, over all 4608 of its '
      'beta-valuations.' % CENSUS_M)
print('NOT RE-RUN: any m >= 12. Nothing in this program, and nothing in the paper, bounds the number '
      'of counterexamples with 12 or more edges.')
print('NOT RE-RUN: the connected form of the question. It is OPEN. What is checked here is only that '
      'every connected bipartite graceful graph with at most %d edges has a near alpha-valuation; '
      'its tree sub-case is an open conjecture of El-Zanati, Kenig and Vanden Eynden.' % CENSUS_M)
print('NOT RE-RUN: vertex-minimality. Graphs on at most 6 vertices with more than 6 edges were never '
      'examined, so no claim is made or checked that 7 is the least possible number of vertices.')
print('NOT RE-RUN: the literature. That the non-existence half is new, and that the gracefulness of '
      'C_4 u P_3 is already published, are bibliographic statements no program can settle.')
print()
if _failed:
    print('VERDICT: %d CHECKS RUN, %d FAILED' % (_passed + _failed, _failed))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _passed)
sys.exit(0)

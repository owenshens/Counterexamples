#!/usr/bin/env python3
"""verify.py -- checks computational claims of

    "A Matching Labelling for the Cyclic Antibandwidth of the 2 x 2 x n Mesh"
    CAB(P_2 x P_2 x P_n) = 2(n-1)   for every integer n >= 2

Python 3.9+, STANDARD LIBRARY ONLY (re, sys, itertools, collections).  No third-party package, no
external data file, and no floating point anywhere: every quantity is an exact integer.

The program reads the two label tables printed in the note -- the n=3 and the n=4 table, pasted in
below character for character and parsed, not regenerated -- and then re-derives the numbers stated
beside them.

WHAT IS CHECKED, section by section:
  A  the printed finite objects: bijectivity, the edge count, the full edge-distance multiset, the
     minimum, and that the printed n=3 / n=4 tables agree with the closed form (2) of
     Proposition 1.  Also that P_2 x P_2 x P_n really is the cylinder C_4 box P_n, and that the
     O(n) neighbour rule used by the sweep produces exactly the edge set of the LITERAL product
     predicate ("exactly one coordinate differs, and differs there by exactly 1") tested over all
     pairs of vertices.
  B  the sweep n = 2..500 -- all 498 tabulated cells n_3 = 3..500 plus the endpoint n = 2:
     |V| = 4n, |E| = 8n-4, the degree multiset, bijectivity of the closed form, and the exact
     edge-distance multiset {2n-2 : 4(n-1), 2n-1 : 2n, 2n : 2n}.
  C  the counting upper bound of Lemma 2, evaluated by brute force in Z_{4n} for every swept n,
     and the resulting exact values, including the four benchmark cells 3, 168, 335, 500.
  D  exhaustive corroboration on small cells: a complete rotation-pinned depth-first search whose
     ONLY pruning is the cyclic-distance constraint itself -- so the upper bound is TESTED rather
     than assumed -- plus the endpoints n = 1 and n = 2 decided by unpruned enumeration.
  E  negative controls in both polarities: the checker and the search must be able to say NO.
  F  the arithmetic of the two-dimensional route recorded in section 5 of the note.
  G  the two tori C_4 box C_3 and C_4 box C_4, by complete search.

Runtime is a few seconds; the sweep touches about a million edges.
Exit status is 0 if and only if every check passed.
"""

import re
import sys
from collections import Counter
from itertools import permutations

sys.setrecursionlimit(10000)

# ---------------------------------------------------------------------------
# 0.  The label tables, pasted from the note
# ---------------------------------------------------------------------------


def formula(u, v, i, n):
    """Proposition 1, equation (2).  Vertices (u,v,i); labels 1..4n."""
    return (1 + 2 * i + u) + 2 * n * ((u + v + i) % 2)


PAPER_N3 = r"""
 1=(0,0,0)   2=(1,1,0)   3=(0,1,1)   4=(1,0,1)
 5=(0,0,2)   6=(1,1,2)   7=(0,1,0)   8=(1,0,0)
 9=(0,0,1)  10=(1,1,1)  11=(0,1,2)  12=(1,0,2)
"""

PAPER_N4 = r"""
 1=(0,0,0)   2=(1,1,0)   3=(0,1,1)   4=(1,0,1)
 5=(0,0,2)   6=(1,1,2)   7=(0,1,3)   8=(1,0,3)
 9=(0,1,0)  10=(1,0,0)  11=(0,0,1)  12=(1,1,1)
13=(0,1,2)  14=(1,0,2)  15=(0,0,3)  16=(1,1,3)
"""

# ---------------------------------------------------------------------------
# 1.  Check bookkeeping
# ---------------------------------------------------------------------------

_passes = 0
_fails = 0


def ck(name, cond, detail=''):
    """One check.  Exactly one PASS/FAIL line, so the verdict count and the PASS lines agree."""
    global _passes, _fails
    if cond:
        _passes += 1
        sys.stdout.write('PASS %s%s\n' % (name, ('   ' + detail) if detail else ''))
    else:
        _fails += 1
        sys.stdout.write('FAIL %s%s\n' % (name, ('   ' + detail) if detail else ''))


def head(text):
    sys.stdout.write('\n--- %s\n' % text)


# ---------------------------------------------------------------------------
# 2.  Parsers for the printed tables
# ---------------------------------------------------------------------------

def parse_triples(block):
    """`label=(u,v,i)` -> {label: (u,v,i)}.  Refuses a repeated label."""
    out = {}
    for a, b, c, d in re.findall(r'(\d+)=\((\d+),(\d+),(\d+)\)', block):
        lab = int(a)
        if lab in out:
            raise ValueError('label %d printed twice' % lab)
        out[lab] = (int(b), int(c), int(d))
    return out


# ---------------------------------------------------------------------------
# 3.  The graph P_2 x P_2 x P_n and its cyclic distances
# ---------------------------------------------------------------------------

def mesh_vertices(n):
    return [(u, v, i) for u in (0, 1) for v in (0, 1) for i in range(n)]


def mesh_edges(n):
    """The three edge families named in the paper, each edge listed once.  |E| = 8n-4."""
    E = []
    for i in range(n):
        for u in (0, 1):                                     # v-edges
            E.append(((u, 0, i), (u, 1, i)))
        for v in (0, 1):                                     # u-edges
            E.append(((0, v, i), (1, v, i)))
        if i + 1 < n:                                        # i-edges
            for u in (0, 1):
                for v in (0, 1):
                    E.append(((u, v, i), (u, v, i + 1)))
    return E


def mesh_edges_literal(n):
    """The LITERAL Cartesian-product predicate over all pairs: the coordinatewise absolute
    differences sum to 1, i.e. exactly one coordinate differs and differs by exactly 1.  O(N^2)."""
    V = mesh_vertices(n)
    E = set()
    for a in range(len(V)):
        x = V[a]
        for b in range(a + 1, len(V)):
            y = V[b]
            if abs(x[0] - y[0]) + abs(x[1] - y[1]) + abs(x[2] - y[2]) == 1:
                E.add((x, y) if x < y else (y, x))
    return E


def cylinder_edges(n):
    """C_4 box P_n on vertices (a,i), a in Z_4, i in 0..n-1."""
    E = set()
    for i in range(n):
        for a in range(4):
            E.add(tuple(sorted(((a, i), ((a + 1) % 4, i)))))
        if i + 1 < n:
            for a in range(4):
                E.add(tuple(sorted(((a, i), (a, i + 1)))))
    return E


def cyc(a, b, N):
    """min(|a-b|, N-|a-b|).  Exact integers only."""
    d = a - b
    if d < 0:
        d = -d
    return d if d + d <= N else N - d


def multiset_from_labelling(n, lab):
    """lab: {(u,v,i): label}.  -> Counter of cyclic distances over E(P_2 x P_2 x P_n)."""
    N = 4 * n
    return Counter(cyc(lab[x], lab[y], N) for x, y in mesh_edges(n))


# ---------------------------------------------------------------------------
# 4.  A complete rotation-pinned search whose only prune is the constraint
# ---------------------------------------------------------------------------

def graph_mesh(n):
    V = mesh_vertices(n)
    idx = {x: j for j, x in enumerate(V)}
    adj = [[] for _ in V]
    for x, y in mesh_edges(n):
        adj[idx[x]].append(idx[y])
        adj[idx[y]].append(idx[x])
    return len(V), adj, idx


def graph_torus(m, n):
    """C_m box C_n, m, n >= 3.  -> (N, adj, idx, sorted edge list on (a,i) labels)."""
    V = [(a, i) for a in range(m) for i in range(n)]
    idx = {x: j for j, x in enumerate(V)}
    E = set()
    for a in range(m):
        for i in range(n):
            for y in (((a + 1) % m, i), (a, (i + 1) % n)):
                E.add(tuple(sorted(((a, i), y))))
    adj = [[] for _ in V]
    for x, y in E:
        adj[idx[x]].append(idx[y])
        adj[idx[y]].append(idx[x])
    return len(V), adj, idx, sorted(E)


def graph_cycle(m):
    return m, [[(a - 1) % m, (a + 1) % m] for a in range(m)]


def bfs_order(N, adj, root):
    seen = [False] * N
    order = [root]
    seen[root] = True
    qi = 0
    while qi < len(order):
        x = order[qi]
        qi += 1
        for y in adj[x]:
            if not seen[y]:
                seen[y] = True
                order.append(y)
    return order


def search(N, adj, k, root=None, count_all=False):
    """Complete depth-first search for a labelling V -> Z_N with every edge at cyclic distance >= k.

    Labels are read in Z_N as 0..N-1.  The paper writes 1..N; cyclic distance is invariant under
    adding a constant, so the conventions agree.  The label of order[0] is PINNED to 0, for the same
    reason: rotation preserves every cyclic distance and every rotation orbit has exactly N members,
    so the number of labellings is N times the number of pinned ones.

    The ONLY pruning is the constraint itself -- no degree argument, no window lemma, no symmetry
    quotient beyond that single rotation.  A count of 0 therefore TESTS the upper bound instead of
    assuming it.

    -> (number of pinned solutions, one witness as a label-per-vertex list, or None)
    """
    if root is None:
        root = max(range(N), key=lambda x: len(adj[x]))
    order = bfs_order(N, adj, root)
    if len(order) != N:
        raise ValueError('the graph is not connected; the search assumes it is')
    allowed = [frozenset(b for b in range(N) if cyc(a, b, N) >= k) for a in range(N)]
    placed = set()
    prev_nb = []
    for x in order:
        prev_nb.append([y for y in adj[x] if y in placed])
        placed.add(x)
    lab = [None] * N
    used = [False] * N
    total = [0]
    witness = [None]

    def rec(t):
        if t == N:
            total[0] += 1
            if witness[0] is None:
                witness[0] = list(lab)
            return not count_all
        x = order[t]
        if t == 0:
            cands = (0,)
        else:
            nb = prev_nb[t]
            s = allowed[lab[nb[0]]]
            for y in nb[1:]:
                s = s & allowed[lab[y]]
            cands = [c for c in s if not used[c]]
        for c in cands:
            lab[x] = c
            used[c] = True
            if rec(t + 1):
                lab[x] = None
                used[c] = False
                return True
            lab[x] = None
            used[c] = False
        return False

    rec(0)
    return total[0], witness[0]


def min_dist(N, adj, lab):
    best = None
    for x in range(N):
        for y in adj[x]:
            if y > x:
                d = cyc(lab[x], lab[y], N)
                if best is None or d < best:
                    best = d
    return best


def brute_force_max_min(N, adj):
    """Exhaustive over all (N-1)! rotation-pinned labellings, with no pruning whatsoever."""
    best = 0
    for perm in permutations(range(1, N)):
        m = min_dist(N, adj, (0,) + perm)
        if m > best:
            best = m
    return best


def brute_force_count(N, adj, k):
    """Exhaustive count of rotation-pinned labellings with minimum >= k, no pruning."""
    c = 0
    for perm in permutations(range(1, N)):
        if min_dist(N, adj, (0,) + perm) >= k:
            c += 1
    return c


# ===========================================================================
#  A.  THE PRINTED FINITE OBJECTS
# ===========================================================================
sys.stdout.write('verify.py -- a matching labelling for CAB(P_2 x P_2 x P_n) = 2(n-1)\n')
sys.stdout.write('interpreter: %s\n' % ' '.join(sys.version.split()))

head('A. the label tables printed in the note, read literally')

for n, block, want_mult, tag in ((3, PAPER_N3, {4: 8, 5: 6, 6: 6}, 'n3'),
                                 (4, PAPER_N4, {6: 12, 7: 8, 8: 8}, 'n4')):
    lab = parse_triples(block)
    N = 4 * n
    ck('%s-table-is-a-labelling' % tag,
       sorted(lab) == list(range(1, N + 1))
       and sorted(lab.values()) == sorted(mesh_vertices(n)),
       'labels 1..%d used once each, on all %d vertices' % (N, N))
    ck('%s-table-matches-the-closed-form' % tag,
       all(formula(u, v, i, n) == l for l, (u, v, i) in lab.items()),
       'every printed label equals the closed form at n=%d' % n)
    mult = multiset_from_labelling(n, {t: l for l, t in lab.items()})
    ck('%s-edge-count' % tag, sum(mult.values()) == 8 * n - 4,
       '|E| = %d = 8*%d-4' % (sum(mult.values()), n))
    ck('%s-distance-multiset' % tag, dict(mult) == want_mult,
       'multiset %s' % sorted(mult.items()))
    ck('%s-minimum-is-2n-2' % tag, min(mult) == 2 * n - 2,
       'min = %d = 2*%d-2' % (min(mult), n))

_fmap = {t: formula(t[0], t[1], t[2], 4) for t in mesh_vertices(4)}

_PHI = {(0, 0): 0, (1, 0): 1, (1, 1): 2, (0, 1): 3}
_cyl_bad = 0
for n in range(2, 31):
    img = set(tuple(sorted(((_PHI[(x[0], x[1])], x[2]), (_PHI[(y[0], y[1])], y[2]))))
              for x, y in mesh_edges(n))
    if img != cylinder_edges(n):
        _cyl_bad += 1
ck('mesh-is-the-cylinder-C4-box-Pn-for-n-2-to-30', _cyl_bad == 0,
   'edge sets coincide under (u,v) -> 0,1,2,3 around the 4-cycle, n = 2..30')

_lit_bad = 0
for n in range(2, 15):
    if set(tuple(sorted(e)) for e in mesh_edges(n)) != mesh_edges_literal(n):
        _lit_bad += 1
ck('neighbour-rule-equals-the-literal-product-predicate-n-2-to-14', _lit_bad == 0,
   'the three O(n) families == the O(N^2) all-pairs product test, n = 2..14')

# ===========================================================================
#  B.  THE SWEEP OVER EVERY TABULATED CELL
# ===========================================================================
head('B. the sweep n = 2..500: the 498 tabulated cells n_3 = 3..500, plus the endpoint n = 2')

NLO, NHI = 2, 500
bad_v = bad_e = bad_deg = bad_bij = bad_mult = bad_min = 0
cells = 0
big_mult = {}
for n in range(NLO, NHI + 1):
    cells += 1
    N = 4 * n
    two_n = 2 * n
    lab = [0] * N
    for i in range(n):
        base = 4 * i
        for u in (0, 1):
            for v in (0, 1):
                lab[base + 2 * u + v] = formula(u, v, i, n)
    if sorted(lab) != list(range(1, N + 1)):
        bad_bij += 1
    deg = [0] * N
    cnt = {}
    for i in range(n):
        base = 4 * i
        for u in (0, 1):                                     # v-edges
            p, q = base + 2 * u, base + 2 * u + 1
            d = cyc(lab[p], lab[q], N)
            cnt[d] = cnt.get(d, 0) + 1
            deg[p] += 1
            deg[q] += 1
        for v in (0, 1):                                     # u-edges
            p, q = base + v, base + 2 + v
            d = cyc(lab[p], lab[q], N)
            cnt[d] = cnt.get(d, 0) + 1
            deg[p] += 1
            deg[q] += 1
        if i + 1 < n:                                        # i-edges
            for j in range(4):
                p, q = base + j, base + 4 + j
                d = cyc(lab[p], lab[q], N)
                cnt[d] = cnt.get(d, 0) + 1
                deg[p] += 1
                deg[q] += 1
    if len(lab) != 4 * n or len(deg) != 4 * n:
        bad_v += 1
    if sum(cnt.values()) != 8 * n - 4:
        bad_e += 1
    want_deg = {3: 8} if n == 2 else {3: 8, 4: 4 * (n - 2)}
    if dict(Counter(deg)) != want_deg:
        bad_deg += 1
    if cnt != {two_n - 2: 4 * (n - 1), two_n - 1: two_n, two_n: two_n}:
        bad_mult += 1
    if min(cnt) != two_n - 2:
        bad_min += 1
    if n in (335, 500):
        big_mult[n] = dict(cnt)

ck('sweep-covers-the-whole-tabulated-row',
   cells == 499 and (500 - 3 + 1) == 498 and NLO == 2 and NHI == 500,
   '%d values of n swept, containing all 498 cells n_3 = 3..500' % cells)
ck('sweep-vertex-and-edge-counts', bad_v == 0 and bad_e == 0,
   '|V| = 4n and |E| = 8n-4 on every cell')
ck('sweep-degree-multiset', bad_deg == 0,
   'eight 3s and 4(n-2) 4s for n >= 3; 3-regular at n = 2; every cell')
ck('sweep-closed-form-is-a-bijection', bad_bij == 0,
   'eq. (2) hits each of 1..4n exactly once, every cell')
ck('sweep-distance-multiset', bad_mult == 0,
   'exactly {2n-2: 4(n-1), 2n-1: 2n, 2n: 2n} on every cell')
ck('sweep-minimum-is-2n-2', bad_min == 0,
   'so CAB(P_2xP_2xP_n) >= 2(n-1) on every cell, n = 2..500')
ck('cell-335-multiset', big_mult.get(335) == {668: 1336, 669: 670, 670: 670},
   'P_2x2x335, 1340 vertices, 2676 edges: %s' % sorted(big_mult.get(335, {}).items()))
ck('cell-500-multiset', big_mult.get(500) == {998: 1996, 999: 1000, 1000: 1000},
   'P_2x2x500, 2000 vertices, 3996 edges: %s' % sorted(big_mult.get(500, {}).items()))

# ===========================================================================
#  C.  THE COUNTING UPPER BOUND
# ===========================================================================
head('C. the counting upper bound of Lemma 2, by brute force in Z_{4n}')

bad_win = bad_mono = 0
for n in range(NLO, NHI + 1):
    N = 4 * n
    k = 2 * n - 1
    win = sum(1 for b in range(N) if cyc(0, b, N) >= k)
    if win != 3 or win != N - 2 * k + 1:
        bad_win += 1
    if sum(1 for b in range(N) if cyc(0, b, N) >= k + 1) > win:
        bad_mono += 1
ck('window-count-is-three', bad_win == 0,
   '#{b in Z_4n : D_c(0,b) >= 2n-1} = 3 = N-2k+1, every n = 2..500')
ck('window-count-is-monotone-in-k', bad_mono == 0,
   'the count cannot grow with k, so ruling out k = 2n-1 rules out every larger k')
ck('lemma-2-side-condition-holds',
   all(2 * n - 1 <= (4 * n) // 2 and 2 * ((2 * n - 1) - 1) < 4 * n
       for n in range(NLO, NHI + 1)),
   'k = 2n-1 <= N/2 and 2(k-1) < N at every swept n')
ck('maximum-degree-is-four-for-n-ge-3',
   all(max(Counter(x for e in mesh_edges(n) for x in e).values()) == 4 for n in range(3, 31)),
   'Delta(G_n) = 4 for n = 3..30, and by the sweep degree multiset for all n >= 3')
ck('upper-bound-follows-from-lemma-2',
   all(4 > 4 * n - 2 * (2 * n - 1) + 1 for n in range(3, NHI + 1)),
   'a degree-4 vertex contradicts Lemma 2 at k = 2n-1: 4 > 3, every n = 3..500')
ck('both-bounds-agree-at-2n-2-on-every-swept-cell',
   bad_min == 0 and bad_win == 0 and all(2 * (n - 1) == 2 * n - 2 for n in range(3, NHI + 1)),
   'the attained minimum and the Lemma 2 ceiling both read 2(n-1), n = 3..500')
ck('formula-values-at-four-cells', [2 * (n - 1) for n in (3, 168, 335, 500)] == [4, 334, 668, 998],
   '2(n-1) at n_3 = 3, 168, 335, 500 is 4, 334, 668, 998')
ck('668-exceeds-667-and-998-exceeds-997', 668 > 667 and 998 > 997,
   'the two quoted integers 667 and 997 are strictly below 668 and 998')

_v333 = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
_e333 = [(x, y) for j, x in enumerate(_v333) for y in _v333[j + 1:]
         if abs(x[0] - y[0]) + abs(x[1] - y[1]) + abs(x[2] - y[2]) == 1]
_deg333 = max(Counter(x for e in _e333 for x in e).values())
ck('lemma-2-is-not-tight-in-general',
   (len(_v333), _deg333, (len(_v333) - _deg333 + 1) // 2) == (27, 6, 11),
   'P_3xP_3xP_3: N=27, Delta=6, the Lemma 2 ceiling is 11')

# ===========================================================================
#  D.  EXHAUSTIVE CORROBORATION
# ===========================================================================
head('D. complete rotation-pinned search, with the constraint as the only prune')

_no_sol = []
for n in range(3, 9):
    N, adj, _ = graph_mesh(n)
    _no_sol.append((n, search(N, adj, 2 * n - 1)[0]))
ck('exhaustive-nothing-attains-2n-1', all(c == 0 for _, c in _no_sol),
   'pinned solution counts at k = 2n-1 for n = 3..8: %s' % [c for _, c in _no_sol])

_counts = []
for n in (3, 4, 5):
    N, adj, _ = graph_mesh(n)
    c = search(N, adj, 2 * n - 2, count_all=True)[0]
    _counts.append((n, c, N * c))
ck('exhaustive-optimum-counts', _counts == [(3, 8, 96), (4, 8, 128), (5, 8, 160)],
   '(n, pinned, total) optimal labellings: %s' % _counts)

_N1, _adj1 = graph_cycle(4)
ck('endpoint-n-1-CAB-C4-is-1', brute_force_max_min(_N1, _adj1) == 1,
   'unpruned enumeration of all 3! = 6 pinned labellings of C_4; 2(n-1) = 0, so n=1 fails')

_N2, _adj2, _ = graph_mesh(2)
ck('endpoint-n-2-CAB-Q3-is-2',
   brute_force_max_min(_N2, _adj2) == 2 and search(_N2, _adj2, 3)[0] == 0,
   'unpruned enumeration of all 7! = 5040 pinned labellings of Q_3 gives 2; '
   'the search finds nothing at k = 3')
ck('search-engine-agrees-with-unpruned-enumeration',
   search(_N2, _adj2, 2, count_all=True)[0] == brute_force_count(_N2, _adj2, 2),
   'pinned count at k = 2 on Q_3, pruned search vs full enumeration')

# ===========================================================================
#  E.  NEGATIVE CONTROLS
# ===========================================================================
head('E. negative controls -- the checker and the search must be able to say NO')

_ident = [(n, min(multiset_from_labelling(
    n, {(u, v, i): 1 + u + 2 * v + 4 * i for (u, v, i) in mesh_vertices(n)})))
    for n in (3, 4, 5, 50)]
ck('control-identity-labelling-scores-1', all(m == 1 for _, m in _ident),
   'naive layer-by-layer labelling, (n, min): %s -- far below 2n-2' % _ident)

_swap = dict(_fmap)
_a = [t for t, l in _fmap.items() if l == 1][0]
_b = [t for t, l in _fmap.items() if l == 2][0]
_swap[_a], _swap[_b] = _swap[_b], _swap[_a]
_swapmin = min(multiset_from_labelling(4, _swap))
ck('control-one-transposition-breaks-the-witness', _swapmin < 6,
   'swapping labels 1 and 2 in the n=4 witness drops the minimum to %d' % _swapmin)

ck('control-the-witness-does-not-over-attain',
   all(min(multiset_from_labelling(
       n, {t: formula(t[0], t[1], t[2], n) for t in mesh_vertices(n)})) < 2 * n - 1
       for n in (3, 4, 5, 10, 100)),
   'the closed form attains 2n-2 and not 2n-1 at n = 3, 4, 5, 10, 100')

# ===========================================================================
#  F.  THE PUBLISHED TWO-DIMENSIONAL ROUTE
# ===========================================================================
head('F. the arithmetic of the two-dimensional route recorded in the note')

ck('floor-is-vacuous-at-n2-equals-4',
   all((4 * (m - 1)) // 2 == 2 * (m - 1) for m in range(4, 501)),
   'floor(n_2(n_1-1)/2) = 2(n_1-1) at n_2 = 4, n_1 = 4..500')

_span_bad = 0
for n in range(2, 41):
    cyl = cylinder_edges(n)
    grid = set(tuple(sorted(((a, i), (a + 1, i)))) for a in range(3) for i in range(n))
    grid |= set(tuple(sorted(((a, i), (a, i + 1)))) for a in range(4) for i in range(n - 1))
    if not (grid < cyl and len(cyl) - len(grid) == n and len(grid) == 7 * n - 4):
        _span_bad += 1
ck('grid-is-a-spanning-subgraph-of-the-cylinder-n-2-to-40', _span_bad == 0,
   'P_4 box P_n is a proper subgraph of C_4 box P_n with exactly n fewer edges, n = 2..40')
ck('edge-count-difference-is-n-at-every-cell',
   all((8 * n - 4) - (7 * n - 4) == n for n in range(2, 501)),
   '(8n-4) - (7n-4) = n at every n = 2..500')

# ===========================================================================
#  G.  THE TORUS REMARK
# ===========================================================================
head('G. the two tori C_4 box C_3 and C_4 box C_4, by complete search')

_N3t, _adj3t, _idx3t, _E3t = graph_torus(4, 3)
ck('torus3-nothing-attains-5',
   max(len(a) for a in _adj3t) == 4 and (12 - 4 + 1) // 2 == 4
   and search(_N3t, _adj3t, 5)[0] == 0,
   '4-regular on 12 vertices: Lemma 2 caps at 4, and the complete search finds nothing at k = 5')

_N4t, _adj4t, _idx4t, _E4t = graph_torus(4, 4)
ck('torus4-does-not-attain-6', search(_N4t, _adj4t, 6)[0] == 0,
   'C_4 box C_4 has no labelling with every edge at cyclic distance >= 6 = 2*4-2')

# ===========================================================================
#  SCOPE, then the verdict
# ===========================================================================
sys.stdout.write('\n')
for line in (
    'NOT RE-RUN: the quoted source statement itself.  The conjectured row, the certified optima',
    'NOT RE-RUN:   4 and 334, and the quoted best-known values 667 and 997 are facts about',
    'NOT RE-RUN:   external documents.  This program checks the MATHEMATICS against those',
    'NOT RE-RUN:   integers; it does not fetch them.',
    'NOT RE-RUN: the priority question.  Nothing here can tell whether an extremal labelling of',
    'NOT RE-RUN:   the 4 x n grid, or a toroidal formula, already yields Proposition 1.',
    'NOT RE-RUN: the exhaustive census beyond n = 8.  Section D decides n = 3..8 and counts the',
    'NOT RE-RUN:   optima at n = 3,4,5 only.  The claim on all 498 cells rests on the PROOF',
    'NOT RE-RUN:   (Proposition 1 and Lemma 2), checked cell by cell in sections B and C; the',
    'NOT RE-RUN:   searches are corroboration, not the load-bearing step.',
    'NOT RE-RUN: the count 32n of optimal labellings.  It is computed at n = 3,4,5 only.',
    'NOT RE-RUN: CAB(C_4 box C_n) beyond n = 4.  Section G searches n = 3 at k = 5 and n = 4 at',
    'NOT RE-RUN:   k = 6 only.',
    'NOT RE-RUN: the other four rows of the conjectured table (3x3, 4x4, 5x5, 6x6).  Only the',
    'NOT RE-RUN:   Lemma 2 ceiling at P_3xP_3xP_3 is computed, as a witness that the bound is not',
    'NOT RE-RUN:   tight in general.',
):
    sys.stdout.write(line + '\n')

sys.stdout.write('\n')
if _fails:
    sys.stdout.write('VERDICT: %d OF %d CHECKS FAILED\n' % (_fails, _passes + _fails))
    sys.exit(1)
sys.stdout.write('VERDICT: ALL %d CHECKS PASS\n' % _passes)
sys.exit(0)

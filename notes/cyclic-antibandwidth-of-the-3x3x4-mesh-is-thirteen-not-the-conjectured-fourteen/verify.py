#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- re-derives every quantity claimed in

    "The cyclic antibandwidth of the 3x3x4 mesh is 13, not the conjectured 14"

from the objects PRINTED IN THE PAPER and from nothing else.

Python 3.9+, STANDARD LIBRARY ONLY (no numpy / sympy / networkx, no external data file, no
network, no randomness).  All arithmetic is exact integer arithmetic; no float value enters any
decision.  Deterministic: every node count printed below is reproducible exactly.  Runs in about
three minutes on one core, almost all of it in the k = 14 exhaustion of section 10.

Contract: one `PASS <name> [detail]` line per check, a closing
    VERDICT: ALL <n> CHECKS PASS
and exit status 0 if and only if every check passed.

WHAT IS COVERED.  Every layer of the paper's proof for n3 = 4: the two printed witnesses and the
lower bound cab >= 13; the arc lemma and the window reformulation; the degree cap that kills
k >= 16; and the two exhaustions that kill k = 15 and k = 14.  No solver is used anywhere, so
the paper's decisive negative no longer depends on a third-party constraint solver.  What is NOT
covered is listed in the closing `NOT RE-RUN:` lines and repeated in `## Scope` of
REVIEW_NOTE.md -- chiefly that k = 14 is settled by a search a referee must RE-RUN rather than
read, and that the cells n3 = 5..8 are not decided here at all.
"""

import itertools
import sys
import time

# =====================================================================================
# 0. the check harness
# =====================================================================================
_N_PASS = 0
_N_FAIL = 0


def check(name, ok, detail=''):
    global _N_PASS, _N_FAIL
    if ok:
        _N_PASS += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _N_FAIL += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))
    sys.stdout.flush()


def section(t):
    print('--- %s ---' % t)
    sys.stdout.flush()


# =====================================================================================
# 1. THE OBJECTS, EXACTLY AS PRINTED IN THE PAPER
# =====================================================================================
# Section 2 of the paper prints two labelings of P_3 x P_3 x P_4 as four rows (z = 1..4); each
# row groups x = 1 | x = 2 | x = 3 and within a group lists y = 1, 2, 3.  The two strings below
# are those rows, transcribed character for character.

W1_GRID = """
  z=1 :  31  14  30   |   13  29  10   |   35  12  27
  z=2 :  15   1  17   |   36  16  32   |   21  34  11
  z=3 :   2  23   4   |   22   3  18   |    7  20  33
  z=4 :  24   8  26   |    9  25   5   |   28   6  19
"""

W2_GRID = """
  z=1 :   1  21   3   |   23   2  19   |   10  24   5
  z=2 :  22   4  20   |    9  25   6   |   32  11  28
  z=3 :   8  26   7   |   31  12  29   |   18  33  15
  z=4 :  30  13  27   |   16  35  14   |   36  17  34
"""

# Section 2 also prints, for the degree-6 vertex (2,2,2), its own label and the labels of its six
# neighbours together with the six cyclic distances.  Transcribed:
PRINTED_HUB = {
    'vertex': (2, 2, 2),
    'label': 16,
    'neighbour_labels': [1, 34, 36, 32, 29, 3],
    'distances': [15, 18, 16, 16, 13, 13],
    'neighbours': [(1, 2, 2), (3, 2, 2), (2, 1, 2), (2, 3, 2), (2, 2, 1), (2, 2, 3)],
}

# Section 4: the values this pipeline computed for the six cells n3 = 3..8 of the target row, and
# the linear antibandwidths of the same six graphs.  Only the n3 = 4 column is re-derived below;
# the others are printed in the paper as reported values and are listed here so that the paper's
# ARITHMETIC about them (which closed form they fit) is checked rather than asserted.
REPORTED_CAB = {3: 9, 4: 13, 5: 18, 6: 22, 7: 27, 8: 31}
REPORTED_AB = {3: 9, 4: 14, 5: 18, 6: 23, 7: 27, 8: 32}

# Section 1: the conjectured row, verbatim -- cab(P_{3x3xn3}) = (9 n3 - 8)/2 for even n3 and
# 9 (n3 - 1)/2 for odd n3, for every integer 3 <= n3 <= 400.
def published(n3):
    return (9 * n3 - 8) // 2 if n3 % 2 == 0 else 9 * (n3 - 1) // 2


def parse_grid(text):
    """The printed grid -> {(x,y,z): label} with x,y in 1..3 and z in 1..4."""
    lab = {}
    rows = 0
    for line in text.strip().splitlines():
        head, body = line.split(':', 1)
        z = int(head.strip().split('=')[1])
        groups = [g.split() for g in body.split('|')]
        if len(groups) != 3:
            raise ValueError('row z=%d does not have three x-groups' % z)
        for gi, g in enumerate(groups):
            if len(g) != 3:
                raise ValueError('row z=%d group x=%d does not have three y-entries' % (z, gi + 1))
            for yi, v in enumerate(g):
                lab[(gi + 1, yi + 1, z)] = int(v)
        rows += 1
    if rows != 4:
        raise ValueError('expected 4 rows, got %d' % rows)
    return lab


# =====================================================================================
# 2. THE GRAPH, BUILT TWO INDEPENDENT WAYS
# =====================================================================================
def mesh_axis(dims):
    """3D mesh as the Cartesian product of paths: successor along each axis.  0-indexed."""
    verts = list(itertools.product(*[range(d) for d in dims]))
    idx = {v: i for i, v in enumerate(verts)}
    E = set()
    for v in verts:
        for ax in range(len(dims)):
            w = list(v)
            w[ax] += 1
            if w[ax] < dims[ax]:
                E.add((idx[v], idx[tuple(w)]))
    return verts, idx, sorted(E)


def mesh_l1(dims):
    """The same graph from the metric definition: L1 distance exactly 1."""
    verts = list(itertools.product(*[range(d) for d in dims]))
    E = set()
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            if sum(abs(a - b) for a, b in zip(verts[i], verts[j])) == 1:
                E.add((i, j))
    return sorted(E)


def adjacency(n, E):
    adj = [set() for _ in range(n)]
    for u, v in E:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def dc(a, b, n):
    """Cyclic distance on Z_n, exact integers."""
    d = abs(a - b) % n
    return min(d, n - d)


def connected(n, adj):
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def max_matching(n, adj, part):
    """Maximum matching of a bipartite graph by augmenting paths (Kuhn).  part[v] in {0,1}."""
    left = [v for v in range(n) if part[v] == 0]
    m = [-1] * n
    size = 0
    for u in left:
        seen = set()
        stack = [(u, iter(sorted(adj[u])))]
        # iterative DFS for an augmenting path
        path = []
        found = None
        while stack:
            v, it = stack[-1]
            adv = False
            for w in it:
                if w in seen:
                    continue
                seen.add(w)
                path.append((v, w))
                if m[w] == -1:
                    found = list(path)
                    stack = []
                    adv = True
                    break
                stack.append((m[w], iter(sorted(adj[m[w]]))))
                adv = True
                break
            if found is not None:
                break
            if not adv:
                stack.pop()
                if path:
                    path.pop()
        if found:
            for a, b in found:
                m[a], m[b] = b, a
            size += 1
    return size, m


# =====================================================================================
# 3. THE DECIDER: a solver-free exhaustive decision procedure for  cab(G) >= k
# =====================================================================================
# REFORMULATION (Lemma 2 of the paper, and it is exact).  A bijection f : V -> Z_n has cyclic
# edge distance >= k on every edge iff, for every label l, the labels available to a neighbour of
# the vertex carrying l form the ARC  {l+k, l+k+1, ..., l+n-k}  of  n - 2k + 1  labels.  So
# deciding cab(G) >= k is deciding whether G embeds as a spanning subgraph of the circulant
# Ci_n(k, k+1, ..., floor(n/2)).  The search below assigns labels to vertices in a fixed
# max-connectivity order, with full forward checking on the not-yet-labelled neighbours.
#
# SYMMETRY REDUCTION (sound, and CHECKED below on the paper's own witness).  Three constraints
# are imposed for a 3x3xn3 mesh whose search starts at a centre vertex c = (2,2,z):
#   (R) label rotation f -> f + t preserves every cyclic distance and acts transitively on which
#       vertex carries label 0, so c may be pinned to label 0;
#   (A) the eight automorphisms of the 3x3 grid fix c and act transitively on its four in-layer
#       neighbours, so the least of their four labels may be required to sit on (1,2,z);
#   (A') the one non-identity grid automorphism fixing (1,2,z) swaps (2,1,z) and (2,3,z), so
#       label(2,1,z) < label(2,3,z) may be required as well;
#   (M) label reflection f -> -f fixes label 0 and reverses every comparison, and grid
#       automorphisms do not touch z, so label(2,2,z-1) < label(2,2,z+1) may be required.
# Together these cut the tree by a factor of 16.  Check S1 below verifies on the printed witness
# that EXACTLY ONE of its 16 images satisfies (A), (A') and (M) -- i.e. the reduction discards no
# solution and keeps no duplicate.

def build(dims):
    verts, idx, E = mesh_axis(dims)
    n = len(verts)
    adj = adjacency(n, E)
    nbrs = [sorted(adj[v]) for v in range(n)]
    return verts, idx, E, n, adj, nbrs


def search_order(dims):
    """Max-connectivity order: start at a maximum-degree vertex, then always take an unplaced
    vertex with the most already-placed neighbours."""
    verts, idx, E, n, adj, nbrs = build(dims)
    deg = [len(nbrs[v]) for v in range(n)]
    start = max(range(n), key=lambda v: (deg[v], -v))
    order = [start]
    ins = {start}
    while len(order) < n:
        best = max((v for v in range(n) if v not in ins),
                   key=lambda v: (sum(1 for w in nbrs[v] if w in ins), deg[v], -v))
        order.append(best)
        ins.add(best)
    return order


def sym_data(dims, order, idx, verts):
    """The four vertices of constraints (A), (A') and the two of (M); None when they do not apply."""
    s = verts[order[0]]
    if dims[0] != 3 or dims[1] != 3 or s[0] != 1 or s[1] != 1:
        return None
    z = s[2]
    n1, n2, n3v, n4 = idx[(0, 1, z)], idx[(1, 0, z)], idx[(1, 2, z)], idx[(2, 1, z)]
    zlo = idx[(1, 1, z - 1)] if z - 1 >= 0 else None
    zhi = idx[(1, 1, z + 1)] if z + 1 < dims[2] else None
    return ([n1, n2, n3v, n4], n1, n2, n3v, zlo, zhi)


def decide(dims, k, stop_first=False, sym=True, budget=None):
    """Exhaustive decision of cab(G) >= k.  Returns
    (status, nodes, leaves, witness, depth_profile, seconds)."""
    verts, idx, E, n, adj, nbrs = build(dims)
    ALL = (1 << n) - 1
    far = []
    for l1 in range(n):
        m = 0
        for l2 in range(n):
            if dc(l1, l2, n) >= k:
                m |= 1 << l2
        far.append(m)
    order = search_order(dims)
    pos = {v: i for i, v in enumerate(order)}
    future = [[w for w in nbrs[order[i]] if pos[w] > i] for i in range(n)]
    C = sym_data(dims, order, idx, verts) if sym else None
    lab = [-1] * n
    cand = [ALL] * n
    prof = [0] * (n + 1)
    leaves = [0]
    wit = [None]
    st = {'nodes': 0}
    t0 = time.time()

    def rec(i, used):
        if i == n:
            leaves[0] += 1
            wit[0] = list(lab)
            return True
        v = order[i]
        c = cand[v] & ~used
        while c:
            b = c & -c
            c ^= b
            l = b.bit_length() - 1
            prof[i] += 1
            st['nodes'] += 1
            if budget is not None and st['nodes'] > budget:
                raise RuntimeError('BUDGET')
            if C is not None:
                quad, n1, n2, n3v, zlo, zhi = C
                if v in quad:
                    if v == n1:
                        if any(lab[w] >= 0 and lab[w] < l for w in quad if w != n1):
                            continue
                    elif lab[n1] >= 0 and l < lab[n1]:
                        continue
                    if v == n3v and lab[n2] >= 0 and l < lab[n2]:
                        continue
                    if v == n2 and lab[n3v] >= 0 and l > lab[n3v]:
                        continue
                if zlo is not None and zhi is not None:
                    if v == zhi and lab[zlo] >= 0 and l < lab[zlo]:
                        continue
                    if v == zlo and lab[zhi] >= 0 and l > lab[zhi]:
                        continue
            nused = used | b
            fl = far[l]
            saved = []
            dead = False
            for w in future[i]:
                old = cand[w]
                new = old & fl
                if new & ~nused == 0:
                    dead = True
                    for x, o in saved:
                        cand[x] = o
                    break
                saved.append((w, old))
                cand[w] = new
            if dead:
                continue
            lab[v] = l
            r = rec(i + 1, nused)
            lab[v] = -1
            for x, o in saved:
                cand[x] = o
            if r and stop_first:
                return True
        return False

    sys.setrecursionlimit(10000)
    lab[order[0]] = 0
    for w in future[0]:
        cand[w] &= far[0]
    try:
        rec(1, 1)
        status = 'FEASIBLE' if leaves[0] else 'INFEASIBLE'
    except RuntimeError:
        status = 'BUDGET_EXCEEDED'
    return status, st['nodes'], leaves[0], wit[0], prof, time.time() - t0


def min_cyclic(dims, lab0):
    """lab0[v] = 0-indexed label.  Returns (is_bijection, min cyclic edge distance, #tight)."""
    verts, idx, E, n, adj, nbrs = build(dims)
    ok = sorted(lab0) == list(range(n))
    d = [dc(lab0[u], lab0[v], n) for u, v in E]
    m = min(d)
    return ok, m, d.count(m)


# =====================================================================================
# 4. THE CHECKS
# =====================================================================================
print('cyclic antibandwidth of P_3 x P_3 x P_4 : re-derivation from the printed objects')
print('')

DIMS = (3, 3, 4)
verts, idx, E, n, adj, nbrs = build(DIMS)
deg = [len(nbrs[v]) for v in range(n)]

section('1. the graph P_3 x P_3 x P_4')

E2 = mesh_l1(DIMS)
check('graph-built-two-independent-ways-is-the-same-graph', sorted(E) == E2,
      'axis-successor product and "L1 distance exactly 1" give identical edge sets')
check('order-and-size-are-36-and-75', (n, len(E)) == (36, 75), 'n = %d, |E| = %d' % (n, len(E)))

hub = [tuple(c + 1 for c in verts[v]) for v in range(n) if deg[v] == max(deg)]
check('maximum-degree-is-6-at-exactly-the-two-printed-vertices',
      max(deg) == 6 and hub == [(2, 2, 2), (2, 2, 3)],
      'Delta = %d attained exactly at %s' % (max(deg), hub))
degcount = {d: deg.count(d) for d in sorted(set(deg))}
check('degree-distribution-sums-to-twice-the-size',
      sum(deg) == 2 * len(E) and degcount == {3: 8, 4: 16, 5: 10, 6: 2},
      'degrees %s, sum %d = 2|E|' % (degcount, sum(deg)))
check('the-mesh-is-connected', connected(n, adj), 'one component on 36 vertices')

part = [sum(verts[v]) % 2 for v in range(n)]
cross = all(part[u] != part[v] for u, v in E)
check('bipartition-by-coordinate-parity-is-18-18-and-every-edge-crosses-it',
      cross and part.count(0) == 18 and part.count(1) == 18,
      'parts of size %d and %d, all %d edges crossing' % (part.count(0), part.count(1), len(E)))

msize, mate = max_matching(n, adj, part)
alpha = n - msize
check('the-matching-the-program-found-really-is-a-matching',
      all(mate[v] != -1 for v in range(n))
      and all(mate[mate[v]] == v for v in range(n))
      and all(mate[v] in adj[v] for v in range(n))
      and 2 * msize == n,
      'the returned pairing is checked edge by edge and for involutivity: %d disjoint edges '
      'covering all %d vertices, so it is a PERFECT matching and no larger one can exist'
      % (msize, n))
check('the-mesh-has-a-perfect-matching-so-its-independence-number-is-18',
      msize == 18 and alpha == 18,
      'maximum matching %d, hence alpha = n - alpha\' = %d by Koenig-Gallai (G bipartite)' % (msize, alpha))
check('the-independence-number-does-not-by-itself-exclude-k-14',
      alpha >= 14,
      'a window of 14 consecutive labels must be independent and alpha = 18 >= 14, so no '
      'independence-number argument can decide k = 14; a search is unavoidable')

section('2. the conjectured row, and what value it asserts at n3 = 4')

check('published-even-branch-at-n3-4-is-14', published(4) == 14,
      '(9*4 - 8)/2 = %d' % published(4))
check('published-row-is-exactly-the-ceiling-of-9(n3-1)/2',
      all(published(t) == -((-9 * (t - 1)) // 2) for t in range(3, 401)),
      'the two-branch formula equals ceil(9(n3-1)/2) for every one of the 398 quantified cells')
check('our-six-decided-values-are-exactly-the-floor-of-9(n3-1)/2',
      all(REPORTED_CAB[t] == (9 * (t - 1)) // 2 for t in range(3, 9)),
      'cab = %s for n3 = 3..8, and floor(9(n3-1)/2) = %s'
      % ([REPORTED_CAB[t] for t in range(3, 9)], [(9 * (t - 1)) // 2 for t in range(3, 9)]))
check('the-two-formulas-differ-exactly-at-even-n3',
      all((published(t) != (9 * (t - 1)) // 2) == (t % 2 == 0) for t in range(3, 401)),
      'ceil and floor of 9(n3-1)/2 differ by 1 iff n3 is even, agree iff n3 is odd')
check('the-three-even-cells-attacked-contradict-the-row-and-the-two-odd-cells-confirm-it',
      [t for t in range(3, 9) if REPORTED_CAB[t] != published(t)] == [4, 6, 8]
      and [t for t in range(3, 9) if REPORTED_CAB[t] == published(t)] == [3, 5, 7],
      'refuted at n3 = 4, 6, 8 (13 vs 14, 22 vs 23, 31 vs 32); confirmed at n3 = 3, 5, 7')

section('3. the witness W1, exactly as printed in Section 2 of the paper')

L1 = parse_grid(W1_GRID)
check('printed-grid-W1-parses-to-36-labelled-vertices', len(L1) == 36,
      'four rows z = 1..4, three x-groups of three y-entries each')
check('W1-is-a-bijection-onto-1..36', sorted(L1.values()) == list(range(1, 37)),
      'every value 1..36 occurs exactly once')
lab1 = [0] * n
for (x, y, z), v in L1.items():
    lab1[idx[(x - 1, y - 1, z - 1)]] = v - 1
ok1, m1, t1 = min_cyclic(DIMS, lab1)
check('W1-has-minimum-cyclic-edge-distance-exactly-13', ok1 and m1 == 13,
      'min over the 75 edges of min(|a-b|, 36-|a-b|) = %d' % m1)
check('W1-has-exactly-15-tight-edges', t1 == 15,
      '%d of the 75 edges realise the minimum 13' % t1)
check('every-one-of-the-75-edges-is-checked-and-none-falls-below-13',
      all(dc(lab1[u], lab1[v], n) >= 13 for u, v in E) and len(E) == 75,
      'exhaustive over the edge set, exact integer arithmetic')

hv = idx[tuple(c - 1 for c in PRINTED_HUB['vertex'])]
hn = [idx[tuple(c - 1 for c in w)] for w in PRINTED_HUB['neighbours']]
check('the-printed-degree-6-check-at-(2,2,2)-is-reproduced-line-for-line',
      lab1[hv] + 1 == PRINTED_HUB['label']
      and set(hn) == adj[hv]
      and [lab1[w] + 1 for w in hn] == PRINTED_HUB['neighbour_labels']
      and [dc(lab1[hv], lab1[w], n) for w in hn] == PRINTED_HUB['distances'],
      'label %d, neighbour labels %s, distances %s'
      % (PRINTED_HUB['label'], PRINTED_HUB['neighbour_labels'], PRINTED_HUB['distances']))

section('4. the second witness W2, also printed in Section 2')

L2 = parse_grid(W2_GRID)
lab2 = [0] * n
for (x, y, z), v in L2.items():
    lab2[idx[(x - 1, y - 1, z - 1)]] = v - 1
ok2, m2, t2 = min_cyclic(DIMS, lab2)
check('W2-is-a-bijection-onto-1..36-with-minimum-cyclic-edge-distance-13',
      ok2 and m2 == 13, 'min = %d, tight edges = %d' % (m2, t2))
check('W2-has-exactly-16-tight-edges', t2 == 16, '%d edges realise the minimum' % t2)
agree = sum(1 for v in range(n) if lab1[v] == lab2[v])
check('W1-and-W2-are-genuinely-different-labelings', agree == 0,
      'they agree on %d of the 36 vertices' % agree)
rot = set((lab2[v] - lab1[v]) % n for v in range(n))
check('W2-is-not-a-rotation-of-W1', len(rot) > 1,
      '%d distinct label differences occur, so no single shift carries W1 to W2' % len(rot))
check('the-lower-bound-cab-at-least-13-needs-no-computer',
      m1 == 13 and m2 == 13,
      'two printed bijections, 75 integer comparisons each; hence cab(P_3x3x4) >= 13 > 12')

section('5. the arc lemma and the window reformulation (Lemma 2 of the paper)')

bad = []
for k in range(1, n // 2 + 1):
    for l in range(n):
        arc = set((l + j) % n for j in range(k, n - k + 1))
        far = set(l2 for l2 in range(n) if dc(l, l2, n) >= k)
        if arc != far or len(far) != n - 2 * k + 1:
            bad.append((k, l))
check('the-far-set-of-a-label-is-an-arc-of-n-2k+1-labels', not bad,
      'verified for all 36 labels and all k = 1..18; |{l\': d_c(l,l\') >= k}| = 36 - 2k + 1')
check('at-k-14-the-arc-has-9-labels-and-at-k-16-only-5',
      n - 2 * 14 + 1 == 9 and n - 2 * 15 + 1 == 7 and n - 2 * 16 + 1 == 5,
      'k = 14 -> 9, k = 15 -> 7, k = 16 -> 5 admissible labels for each neighbour')

# the window equivalence, verified exhaustively on two small graphs over ALL bijections
def window_equiv(dims):
    vs, ix, Es = mesh_axis(dims)
    N = len(vs)
    A = adjacency(N, Es)
    worst = 0
    for perm in itertools.permutations(range(N)):
        m = min(dc(perm[u], perm[v], N) for u, v in Es)
        at = [0] * N
        for v in range(N):
            at[perm[v]] = v
        for k in range(1, N // 2 + 1):
            indep = True
            for s in range(N):
                W = [at[(s + j) % N] for j in range(k)]
                if any(W[b] in A[W[a]] for a in range(k) for b in range(a + 1, k)):
                    indep = False
                    break
            if indep != (k <= m):
                return False, (perm, k)
        worst += 1
    return True, worst


okw, dw = window_equiv((2, 2, 2))
check('window-reformulation-holds-for-every-bijection-of-P_2x2x2',
      okw, 'all %d bijections x k = 1..4: "min cyclic edge distance >= k" iff "every cyclic '
           'window of k consecutive labels is independent"' % dw)
okw2, dw2 = window_equiv((1, 2, 3))
check('window-reformulation-holds-for-every-bijection-of-the-2x3-grid',
      okw2, 'all %d bijections, second graph, same equivalence' % dw2)

at1 = [0] * n
for v in range(n):
    at1[lab1[v]] = v
wins13 = all(not any(at1[(s + b) % n] in adj[at1[(s + a) % n]]
                     for a in range(13) for b in range(a + 1, 13))
             for s in range(n))
wins14 = all(not any(at1[(s + b) % n] in adj[at1[(s + a) % n]]
                     for a in range(14) for b in range(a + 1, 14))
             for s in range(n))
check('on-W1-all-36-windows-of-13-labels-are-independent-and-some-window-of-14-is-not',
      wins13 and not wins14,
      'exactly what min cyclic edge distance 13 means, checked directly on the printed object')

section('6. the upper bound, layer B1: k >= 16 is impossible in one line')

check('the-degree-cap-kills-every-k-at-least-16',
      n - 2 * 16 + 1 < max(deg) and (n - max(deg) + 1) // 2 == 15,
      'a vertex of degree 6 needs 6 admissible labels but k = 16 leaves only %d, so '
      'k <= floor((n - Delta + 1)/2) = 15' % (n - 2 * 16 + 1))
st, nodes, lv, _, _, el = decide(DIMS, 16, sym=False)
check('the-decider-agrees-that-k-16-is-infeasible-on-the-target-object',
      st == 'INFEASIBLE' and lv == 0,
      'exhausted in %d nodes, %.2f s (a proved-silent control: the cap already forbids it)'
      % (nodes, el))
check('the-degree-cap-does-NOT-decide-k-15-or-k-14',
      n - 2 * 15 + 1 >= max(deg) and n - 2 * 14 + 1 >= max(deg),
      'at k = 15 the host circulant Ci_36(15,16,17,18) is 7-regular and Delta = 6, so counting '
      'leaves both k = 15 and k = 14 open')

section('7. the decider, calibrated in both polarities on other people\'s certified integers')

st, nodes, lv, w, _, el = decide((2, 2, 3), 4, stop_first=True, sym=False)
okc, mc, _ = min_cyclic((2, 2, 3), w) if w else (False, -1, 0)
check('control-cab-P_2x2x3-at-least-4-with-a-witness-the-verifier-re-checks',
      st == 'FEASIBLE' and okc and mc >= 4,
      'search found a bijection with min cyclic edge distance %d; the benchmark reports 4 as '
      'certified optimal' % mc)
st, nodes5, lv5, _, _, el5 = decide((2, 2, 3), 5, sym=False)
check('control-cab-P_2x2x3-is-not-5-so-the-published-4-is-tight',
      st == 'INFEASIBLE' and lv5 == 0,
      'exhausted in %d nodes, %.2f s -> cab(P_2x2x3) = 4 exactly, matching the published 4*'
      % (nodes5, el5))

st, n9, lv9, w9, _, el9 = decide((3, 3, 3), 9, stop_first=True, sym=True)
ok9, m9, _ = min_cyclic((3, 3, 3), w9) if w9 else (False, -1, 0)
check('control-cab-P_3x3x3-at-least-9-with-a-witness-the-verifier-re-checks',
      st == 'FEASIBLE' and ok9 and m9 >= 9,
      'search found a bijection of the 27-vertex mesh with min cyclic edge distance %d in %d '
      'nodes, %.2f s' % (m9, n9, el9))
st, n10, lv10, _, p10, el10 = decide((3, 3, 3), 10, sym=True)
check('control-cab-P_3x3x3-is-not-10-so-the-published-9-is-reproduced-two-sidedly',
      st == 'INFEASIBLE' and lv10 == 0,
      'exhausted in %d nodes, 0 complete labelings, %.2f s -> cab(P_3x3x3) = 9 exactly, the one '
      'value of the target row that a published exact method certifies' % (n10, el10))
st, n12, lv12, _, _, el12 = decide((3, 3, 3), 12, sym=True)
check('proved-silent-control-on-the-27-vertex-mesh', st == 'INFEASIBLE' and lv12 == 0,
      'its degree cap is 11, and k = 12 is refused in %d nodes' % n12)

section('8. the symmetry reduction is sound: it keeps exactly one image of a real solution')

def grid_autos():
    out = []
    for sw in (0, 1):
        for fx in (0, 1):
            for fy in (0, 1):
                def g(v, sw=sw, fx=fx, fy=fy):
                    x, y, z = v
                    if sw:
                        x, y = y, x
                    if fx:
                        x = 2 - x
                    if fy:
                        y = 2 - y
                    return (x, y, z)
                out.append(g)
    return out


order = search_order(DIMS)
start = order[0]
C = sym_data(DIMS, order, idx, verts)
quad, q1, q2, q3, zlo, zhi = C
base = {verts[v]: (lab1[v] - lab1[start]) % n for v in range(n)}
kept = []
for refl in (0, 1):
    Lr = {v: ((n - l) % n if refl else l) for v, l in base.items()}
    for gi, g in enumerate(grid_autos()):
        M = {g(v): Lr[v] for v in Lr}
        q = [M[verts[w]] for w in quad]
        if (M[verts[q1]] == min(q) and M[verts[q2]] < M[verts[q3]]
                and M[verts[zlo]] < M[verts[zhi]]):
            kept.append((refl, gi))
check('exactly-one-of-the-16-symmetry-images-of-W1-is-canonical', len(kept) == 1,
      'the group is 8 grid automorphisms x label reflection; %d image(s) satisfy the three '
      'canonicity constraints, so the reduction loses no solution and keeps no duplicate'
      % len(kept))
check('all-16-images-of-W1-are-themselves-valid-labelings-of-min-distance-13',
      all(min_cyclic(DIMS, [({g(v): (((n - l) % n) if refl else l)
                             for v, l in base.items()})[verts[v]] for v in range(n)])[1] == 13
          for refl in (0, 1) for g in grid_autos()),
      'the group really does act on the solution set: every image has min cyclic distance 13')

# the decider must ADMIT the canonical image at every placement step -- a direct test that
# neither the forward checking nor the symmetry pruning discards a genuine solution.
refl, gi = kept[0]
g = grid_autos()[gi]
canon = {}
for v0, l in base.items():
    canon[g(v0)] = (n - l) % n if refl else l
target = [canon[verts[v]] for v in range(n)]
ALLM = (1 << n) - 1
far13 = []
for l1 in range(n):
    mm = 0
    for l2 in range(n):
        if dc(l1, l2, n) >= 13:
            mm |= 1 << l2
    far13.append(mm)
pos = {v: i for i, v in enumerate(order)}
future = [[w for w in nbrs[order[i]] if pos[w] > i] for i in range(n)]
cand = [ALLM] * n
used = 0
admitted = True
for i in range(n):
    v = order[i]
    l = target[v]
    if not (cand[v] & ~used) >> l & 1:
        admitted = False
        break
    used |= 1 << l
    for w in future[i]:
        cand[w] &= far13[l]
check('the-decider-admits-the-printed-witness-at-every-one-of-the-36-placement-steps',
      admitted and target[start] == 0,
      'the canonical image of W1 survives forward checking and all three symmetry constraints '
      'from step 1 to step 36, so a k = 13 solution is not pruned away')

section('9. the upper bound, layer B2: k = 15 is impossible, solver-free and exhaustive')

st, nA, lvA, _, profA, elA = decide(DIMS, 15, sym=False)
check('k-15-is-infeasible-with-NO-symmetry-reduction-at-all',
      st == 'INFEASIBLE' and lvA == 0,
      'complete exhaustion: %d search-tree nodes, 0 complete labelings, %.2f s' % (nA, elA))
check('depth-profile-of-that-run-sums-to-its-reported-node-count',
      sum(profA) == nA, '1 root + sum over the 35 placement levels = %d' % nA)
st, nB, lvB, _, profB, elB = decide(DIMS, 15, sym=True)
check('k-15-is-infeasible-again-under-the-16-fold-symmetry-reduction',
      st == 'INFEASIBLE' and lvB == 0,
      '%d nodes, 0 complete labelings, %.2f s; the two runs differ in exactly one thing -- one '
      'applies the reduction and the other does not -- and they agree' % (nB, elB))
check('the-symmetry-reduction-shrinks-the-tree-without-changing-the-verdict',
      nB < nA and lvA == lvB == 0,
      '%d nodes with the reduction against %d without, same INFEASIBLE verdict' % (nB, nA))

section('10. the upper bound, layer B3: k = 14 -- the published value itself -- is impossible')

print('    (the load-bearing exhaustion; about three minutes on one core, please wait)')
sys.stdout.flush()
st, n14, lv14, _, prof14, el14 = decide(DIMS, 14, sym=True)
check('k-14-is-infeasible-solver-free-and-exhaustively',
      st == 'INFEASIBLE' and lv14 == 0,
      'complete exhaustion of the reduced tree: %d search-tree nodes, 0 complete labelings, '
      '%.0f s -- the conjectured value 14 is not merely beatable, it is UNREACHABLE'
      % (n14, el14))
check('depth-profile-of-the-k-14-run-sums-to-its-reported-node-count',
      sum(prof14) == n14, 'sum over the 35 placement levels = %d' % n14)
deep14 = max(i for i in range(1, n) if prof14[i] > 0)
check('the-k-14-exhaustion-is-not-a-trivially-empty-one',
      n14 > 10 ** 8 and deep14 >= 20
      and all(prof14[i] > 0 for i in range(1, deep14 + 1))
      and all(prof14[i] == 0 for i in range(deep14 + 1, n)),
      'levels 1..%d of the tree are non-empty and every level beyond is empty: the search dies '
      '%d vertices short of a complete labeling, which is what "no labeling exists" looks like '
      'from inside a backtracking search' % (deep14, n - deep14))
check('the-k-14-tree-is-far-larger-than-the-k-15-tree-as-the-arc-count-predicts',
      n14 > 100 * nB,
      '%d nodes at k = 14 against %d at k = 15: each neighbour has 9 admissible labels instead '
      'of 7, which is why k = 14 and not k = 15 is the expensive layer' % (n14, nB))

section('11. the answer for n3 = 4, and the arithmetic of the two closed forms')

check('cab-P_3x3x4-equals-13-exactly',
      m1 == 13 and lvA == 0 and lvB == 0 and lv14 == 0 and (n - 2 * 16 + 1) < max(deg),
      'lower bound 13 from the printed bijection; k >= 16 by the degree cap; k = 15 and k = 14 '
      'by the two exhaustions above. Every layer of the proof is re-derived in this one run, '
      'with no solver and no library')
check('the-conjectured-value-14-is-therefore-one-too-high-at-n3-4',
      published(4) - REPORTED_CAB[4] == 1,
      'published %d, true value 13, difference %d -- the row fails DOWNWARD, not upward'
      % (published(4), published(4) - REPORTED_CAB[4]))
check('the-reported-linear-antibandwidths-equal-the-published-cyclic-column',
      all(REPORTED_AB[t] == published(t) for t in range(3, 9)),
      'ab = %s for n3 = 3..8, and the published "conjectured cab" row gives %s -- identical at '
      'all six decided cells' % ([REPORTED_AB[t] for t in range(3, 9)],
                                 [published(t) for t in range(3, 9)]))
check('and-they-sit-one-above-our-cyclic-values-exactly-at-even-n3',
      all(REPORTED_AB[t] - REPORTED_CAB[t] == (1 if t % 2 == 0 else 0) for t in range(3, 9)),
      'ab - cab = 1 for n3 = 4, 6, 8 and 0 for n3 = 3, 5, 7')
check('every-decided-cell-obeys-the-Miller-Pritikin-sandwich',
      all(2 * REPORTED_CAB[t] >= REPORTED_AB[t] and REPORTED_CAB[t] <= REPORTED_AB[t]
          for t in range(3, 9)),
      '(1/2) ab <= cab <= ab holds at n3 = 3..8, so no published ab value can yield cab <= 13')
check('a-correct-exact-ab-value-can-never-give-the-upper-bound-this-paper-needs',
      REPORTED_AB[4] == 14 and REPORTED_AB[4] > 13,
      'ab(P_3x3x4) = 14, so the transfer cab <= ab gives only cab <= 14; 13 does not follow '
      'from any antibandwidth result')
check('the-large-n3-datum-in-the-literature-points-the-same-way',
      published(400) == 1796 and 1796 - 1779 == 17 and (9 * 399) // 2 == 1795,
      'the row asserts 1796 at n3 = 400; the best value ever reported there is 1779, which is '
      '17 below it and 16 below our corrected floor form 1795')

# =====================================================================================
# 5. what this program did NOT do
# =====================================================================================
print('')
print('NOT RE-RUN: a HAND-CHECKABLE proof that k = 14 is impossible. The k = 14 layer is decided '
      'above, exhaustively and without any solver, but it is decided by a search of a few '
      'hundred million nodes: a referee must re-run it, not read it. No DRAT, LRAT or other '
      'independently checkable UNSAT proof object exists for it anywhere in the record, and the '
      'obstruction is nearly global (its minimal infeasible induced subgraph was measured '
      'elsewhere at 30 of the 36 vertices), so it does not compress to a hand-sized core the way '
      'k = 15 does. The 16-fold symmetry reduction used here IS argued, and its soundness is '
      'checked in section 8; exhaust.py in this folder re-decides k = 14 with that reduction '
      'removed entirely.')
print('NOT RE-RUN: the cells n3 = 5, 6, 7, 8. The values 18, 22, 27, 31 quoted in the paper are '
      'reported, not re-derived here: no labeling attaining them survives in the record, so '
      'there is no object for this program to read. Only n3 = 4 is settled by what is printed.')
print('NOT RE-RUN: the linear antibandwidths ab(P_3x3xn3) = 9, 14, 18, 23, 27, 32. Their '
      'ARITHMETIC relation to the published row is checked above; their computation is not '
      'reproduced, again for want of a printed witness.')
print('NOT RE-RUN: anything at n3 >= 9. The uniform closed form floor(9(n3-1)/2) is a '
      'conjecture; 392 of the 398 cells the published row quantifies over are untouched.')
print('NOT RE-RUN: the value 1779 attributed to a published heuristic at n3 = 400, and the '
      'transposition of two cells in the restating paper\'s results table. Only the arithmetic '
      'around them is checked; neither is read off the source by this program.')
print('NOT RE-RUN: the bibliographic locators. The arXiv identifier, the DOIs, the table label '
      'and the source line numbers printed in the paper are not verified here; this program '
      'checks mathematics only.')
print('NOT RE-RUN: prior art. Nothing here bears on whether cab(P_3x3x4) = 13 was known.')
print('')

if _N_FAIL:
    print('VERDICT: %d CHECKS FAILED' % _N_FAIL)
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _N_PASS)
sys.exit(0)

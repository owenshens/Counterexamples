#!/usr/bin/env python3
"""verify.py -- independent check of

    "A four-vertex counterexample to the Clow--van Bommel eternal distance-k
     domination conjecture"

Python 3.9+, STANDARD LIBRARY ONLY, exact integer arithmetic throughout (there is no float
anywhere in this file and no decision is taken on one).  No external data file: every object
it checks is the object PRINTED IN THE NOTE, entered below as an explicit adjacency list or
built from the printed path labelling.

What it re-derives, from the objects alone:
  * the exhibited trees are trees, with the printed order, size, diameter and radius;
  * gamma_j for each j the note names, by EXACT minimum set cover (branch on the lowest
    uncovered vertex; no heuristic, no bound taken on faith);
  * gamma_{all,k}^infty for each witness, by a greatest-fixed-point solution of the game as
    DEFINED IN SECTION 1 of the note, over the FULL configuration space (all multisets of
    the given size);
  * the premise and the conclusion of Conjecture 8.2 on each witness;
  * instances of Proposition 2 (odd k) at k=3,5,7 and on one non-path tree.

The eternal game is solved from the definition and not from any closed form.

Every check prints one `PASS <name>` line.  The program exits 0 iff every check passed and
closes with `VERDICT: ALL <n> CHECKS PASS`.
"""
import sys
from itertools import combinations_with_replacement, permutations

CHECKS = []
FAILED = []


def PASS(name, detail=''):
    CHECKS.append(name)
    print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))


def CHECK(name, cond, detail=''):
    if cond:
        PASS(name, detail)
    else:
        FAILED.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def NOTE(text):
    print('NOTE ' + text)


# ---------------------------------------------------------------------------
# 0.  THE OBJECTS EXHIBITED IN THE PAPER
# ---------------------------------------------------------------------------
# Section 2, Theorem 1: T = P_4, the path 1-2-3-4.  It is entered here as a 0-indexed
# adjacency list, i.e. note vertex i is index i-1.
W1_ADJ = [[1], [0, 2], [1, 3], [2]]                      # P_4


def path_adj(n):
    """The path 1-2-...-n as printed in the paper, 0-indexed."""
    adj = [[] for _ in range(n)]
    for i in range(n - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    return adj


def is_tree(adj):
    n = len(adj)
    if n == 0:
        return False
    edges = sum(len(a) for a in adj)
    if edges % 2:
        return False
    if edges // 2 != n - 1:
        return False
    for v in range(n):
        if v in adj[v] or len(set(adj[v])) != len(adj[v]):
            return False
        for u in adj[v]:
            if v not in adj[u]:
                return False
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for u in adj[v]:
            if u not in seen:
                seen.add(u)
                stack.append(u)
    return len(seen) == n


def edge_set(adj):
    return sorted({(min(u, v), max(u, v)) for u in range(len(adj)) for v in adj[u]})


def dists(adj):
    """All-pairs distances by BFS.  Integers only."""
    n = len(adj)
    D = []
    for s in range(n):
        d = [-1] * n
        d[s] = 0
        frontier = [s]
        while frontier:
            nxt = []
            for v in frontier:
                for u in adj[v]:
                    if d[u] < 0:
                        d[u] = d[v] + 1
                        nxt.append(u)
            frontier = nxt
        D.append(d)
    return D


def diam_rad(D):
    ecc = [max(row) for row in D]
    return max(ecc), min(ecc)


def balls(D, r):
    """N_r[v] as a bitmask, for every v."""
    n = len(D)
    return [sum(1 << u for u in range(n) if D[v][u] <= r) for v in range(n)]


# ---------------------------------------------------------------------------
# 1.  gamma_j  --  EXACT minimum set cover
# ---------------------------------------------------------------------------
def gamma_dist(D, j):
    """gamma_j(G) = fewest balls N_j[v] whose union is V(G).  Exact: iterative deepening on
    the budget, branching only on balls that cover the LOWEST still-uncovered vertex (which
    is complete, since some chosen ball must cover it)."""
    n = len(D)
    B = balls(D, j)
    full = (1 << n) - 1
    cov_by = [[v for v in range(n) if (B[v] >> u) & 1] for u in range(n)]

    def dfs(cov, budget):
        if cov == full:
            return True
        if budget == 0:
            return False
        u = 0
        while (cov >> u) & 1:
            u += 1
        for v in cov_by[u]:
            if dfs(cov | B[v], budget - 1):
                return True
        return False

    for size in range(1, n + 1):
        if dfs(0, size):
            return size
    return None            # unreachable for a connected graph


def dominating(D, S, j):
    """Is S a distance-j dominating set?  (Used to certify the explicit sets the paper names.)"""
    n = len(D)
    return all(any(D[s][u] <= j for s in S) for u in range(n))


# ---------------------------------------------------------------------------
# 2.  gamma_{all,k}^infty  --  the game of Section 1, solved as a greatest fixed point
# ---------------------------------------------------------------------------
def game_wins(D, k, g):
    """True iff g guards admit an eternal distance-k dominating family, in the model of
    Section 1: configurations are MULTISETS of size g; EVERY member must be distance-k
    dominating (condition (i)); a response to an attack at v is a member containing v
    reachable by a bijection of guards each moving at most k (condition (ii), all guards
    move).

    Correctness: the set of configurations that survive is the GREATEST fixed point of
    "for every attack v there is a live successor containing v", and any eternal family is
    contained in it, so the answer is exact -- not a strategy search with a horizon.
    """
    n = len(D)
    B = balls(D, k)
    full = (1 << n) - 1
    cfg = []
    for c in combinations_with_replacement(range(n), g):
        m = 0
        for v in c:
            m |= B[v]
        if m == full:
            cfg.append(c)
    if not cfg:
        return False
    N = len(cfg)
    perms = list(permutations(range(g)))
    succ = [0] * N
    for i in range(N):
        A = cfg[i]
        bits = 0
        for j2 in range(N):
            C = cfg[j2]
            for p in perms:
                if all(D[A[t]][C[p[t]]] <= k for t in range(g)):
                    bits |= (1 << j2)
                    break
        succ[i] = bits
    holds = [0] * n                     # holds[v] = configurations containing v
    for j2 in range(N):
        for v in set(cfg[j2]):
            holds[v] |= (1 << j2)
    alive = (1 << N) - 1
    while True:
        drop = 0
        for i in range(N):
            if not (alive >> i) & 1:
                continue
            s = succ[i] & alive
            for v in range(n):
                if not (s & holds[v]):
                    drop |= (1 << i)
                    break
        if not drop:
            return True
        alive &= ~drop
        if alive == 0:
            return False


def gamma_all(D, k, cap=4):
    """Least g <= cap with game_wins, or None."""
    for g in range(1, cap + 1):
        if game_wins(D, k, g):
            return g
    return None


# ===========================================================================
print('verification of "A four-vertex counterexample to the Clow--van Bommel eternal')
print('distance-k domination conjecture" -- objects: P_4 at k=3, P_{k+1} for odd k')
print('python %d.%d, exact integer arithmetic only' % sys.version_info[:2])
print('')
print('=== Step 1: the exhibited object of Theorem 1 is P_4, and its printed invariants')

D1 = dists(W1_ADJ)
CHECK('W1_is_a_tree', is_tree(W1_ADJ), 'n=%d, %d edges' % (len(W1_ADJ), len(edge_set(W1_ADJ))))
CHECK('W1_order_is_4', len(W1_ADJ) == 4)
CHECK('W1_is_the_printed_path_1_2_3_4',
      edge_set(W1_ADJ) == edge_set(path_adj(4)),
      'edges (0-indexed) = %s' % (edge_set(W1_ADJ),))
d1, r1 = diam_rad(D1)
CHECK('W1_diameter_is_3_equals_k', d1 == 3, 'diam=%d, k=3' % d1)
CHECK('W1_radius_is_2_exceeds_floor_k_over_2', r1 == 2 and r1 > 3 // 2,
      'rad=%d > floor(3/2)=%d' % (r1, 3 // 2))

print('')
print('=== Step 2: Theorem 1 -- gamma_3(P_4) = gamma_all,3^inf(P_4) = 1 and gamma(P_4) = 2')

g3_w1 = gamma_dist(D1, 3)
CHECK('W1_gamma_3_is_1', g3_w1 == 1, 'gamma_3(P_4)=%d' % g3_w1)
CHECK('W1_every_singleton_distance_3_dominates',
      all(dominating(D1, [v], 3) for v in range(4)),
      'all 4 singletons are distance-3 dominating, as the proof says')
g1_w1 = gamma_dist(D1, 1)
CHECK('W1_gamma_1_is_2', g1_w1 == 2, 'gamma(P_4)=%d' % g1_w1)
CHECK('W1_no_single_vertex_dominates',
      not any(dominating(D1, [v], 1) for v in range(4)),
      'N_1[v] for v=1..4 = %s' % ([sorted(u + 1 for u in range(4) if D1[v][u] <= 1)
                                   for v in range(4)],))
CHECK('W1_the_printed_dominating_set_2_3_works', dominating(D1, [1, 2], 1),
      '{2,3} in paper labels')
CHECK('W1_gamma_all_3_is_1_by_the_game_solver', game_wins(D1, 3, 1),
      'one guard survives the greatest-fixed-point iteration')
ga_w1 = gamma_all(D1, 3)
CHECK('W1_gamma_all_3_least_value_is_1', ga_w1 == 1, 'gamma_all,3^inf(P_4)=%d' % ga_w1)
CHECK('W1_premise_of_conjecture_8_2_holds', g3_w1 == ga_w1, '%d = %d' % (g3_w1, ga_w1))
CHECK('W1_conclusion_of_conjecture_8_2_fails', ga_w1 != g1_w1,
      'conclusion would need %d = %d' % (ga_w1, g1_w1))
CHECK('W1_sandwich_of_CMM_respected', g3_w1 <= ga_w1 <= g1_w1,
      '%d <= %d <= %d' % (g3_w1, ga_w1, g1_w1))
CHECK('W1_refutes_conjecture_8_2', (g3_w1 == ga_w1) and (ga_w1 != g1_w1) and 3 > 2,
      'k=3>2, tree, premise true, conclusion false')

print('')
print('=== Step 3: Proposition 2 -- instances at odd k = 3, 5, 7 and one non-path tree')

for k in (3, 5, 7):
    T = path_adj(k + 1)
    D = dists(T)
    dm, rd = diam_rad(D)
    gk = gamma_dist(D, k)
    gf = gamma_dist(D, k // 2)
    ga = gamma_all(D, k)
    CHECK('oddk_instance_P%d_at_k%d_is_a_counterexample' % (k + 1, k),
          dm == k and rd == (k + 1) // 2 and gk == 1 and ga == 1 and gf >= 2,
          'diam=%d rad=%d gamma_%d=%d gamma_all=%d gamma_%d=%d' % (dm, rd, k, gk, ga, k // 2, gf))

for k in (4, 6, 8):
    T = path_adj(k + 1)
    D = dists(T)
    dm, rd = diam_rad(D)
    gf = gamma_dist(D, k // 2)
    ga = gamma_all(D, k)
    CHECK('evenk_instance_P%d_at_k%d_is_NOT_a_counterexample' % (k + 1, k),
          dm == k and rd == k // 2 and gf == 1 and ga == 1,
          'diam=%d rad=%d=floor(k/2) gamma_%d=%d gamma_all=%d -- conclusion HOLDS'
          % (dm, rd, k // 2, gf, ga))

# The parity identity the paper names as the mechanism.
CHECK('parity_identity_2_floor_k_over_2_plus_1',
      all((2 * (k // 2) + 1) == (k + 1 if k % 2 == 0 else k) for k in range(2, 200)),
      '2*floor(k/2)+1 = k+1 for even k and k for odd k, all k in 2..199')

# Non-path trees of diameter exactly k, to show Proposition 2 is not about paths.
DS = [[1, 2, 3], [0, 4, 5], [0], [0], [1], [1]]        # double star S(2,2), diameter 3
Dds = dists(DS)
dm, rd = diam_rad(Dds)
CHECK('oddk_double_star_S22_at_k3_is_a_counterexample',
      is_tree(DS) and dm == 3 and gamma_dist(Dds, 3) == 1 and gamma_all(Dds, 3) == 1
      and gamma_dist(Dds, 1) == 2,
      'n=6 double star, diam=%d rad=%d, (1,1,2)' % (dm, rd))

print('')
NOTE('SCOPE: Proposition 2 is stated for ALL odd k >= 3 and ALL trees of diameter exactly '
     'k; Step 3 instantiates it at k=3,5,7 on P_{k+1} and on one non-path tree. The general '
     'statement rests on the proof in the note, not on these instances. Nothing here bears '
     'on the even-k restriction of Conjecture 8.2, and no claim about the contents of the '
     'cited papers is checked by this program.')

print('')
if FAILED:
    print('VERDICT: %d CHECK(S) FAILED: %s' % (len(FAILED), FAILED))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(CHECKS))
sys.exit(0)

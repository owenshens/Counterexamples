#!/usr/bin/env python3
"""verify.py -- re-derives every quantity claimed in paper.tex, from the object printed there.

Python 3.9+, STANDARD LIBRARY ONLY (itertools, fractions, sys).  No external data file: every input is a
literal below, copied from the paper, so the program and the paper can be diffed by eye.

Exact arithmetic throughout: pebble counts are integers, potentials and losses are fractions.Fraction.
No floating point is used in any decision.

WHAT IT DOES
  1. rebuilds the weighted triangle G = K_3 with w(1,2)=2, w(1,3)=3, w(2,3)=2 and its square G x G;
  2. recomputes pi(G,v) at all three roots and pi(G) by exhaustive search;
  3. re-decodes the graph6 string printed in the paper and checks it is the product's skeleton;
  4. checks the printed distribution D is (1,1)-unsolvable TWICE -- once by a plain forward search over
     every state reachable from D (no potential, no box, no closure), once by the paper's own potential
     argument, mechanised: every move of a successful play has loss <= 1/9, only 14 of the 36 directed
     product moves do, and the root is unreachable from D using only those;
  5. recomputes the exact value pi(G x G,(1,1)) = 10, the level counts of unsolvable distributions, and
     pi(G x G, r) at all nine roots;
  6. re-runs the controls: five published unweighted pebbling numbers at w == 2, Chung's Theorem 3 over a
     12 x 12 grid of weighted K_2 products, the source paper's own weighted-K_4 cover-pebbling example,
     and three anti-controls that must NOT produce a violation;
  7. re-checks the two further witnesses and the exact sense in which they are weaker than the first.

Exits 0 if and only if every check passes.
"""

import itertools
import sys
from fractions import Fraction as Fr

# ---------------------------------------------------------------------------
# THE OBJECT, EXACTLY AS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
W_G = {(1, 2): 2, (1, 3): 3, (2, 3): 2}        # paper, Section 2, first line
VERTEX_ORDER = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
D_PRINTED = (0, 0, 0, 0, 0, 1, 0, 1, 7)        # 7 on (3,3), one each on (2,3) and (3,2)
D_VEC_PRINTED = (1, 2, 3, 2, 4, 6, 3, 6, 9)    # min-product distances to the root (1,1)
G6_SKELETON = 'H{S{aSf'                        # graph6 of the product's unweighted skeleton
PI_FACTOR_PRINTED = (3, 3, 3)                  # pi(G,1), pi(G,2), pi(G,3)
LARGEST_UNSOLVABLE_PRINTED = {1: [(0, 0, 2), (0, 1, 1)],      # every unsolvable distribution of maximum
                              2: [(1, 0, 1)],                 # size, at each of the three roots, as the
                              3: [(1, 1, 0), (2, 0, 0)]}      # paper's Lemma 1 lists them
PI_PRODUCT_ROOT_PRINTED = 10
LEVEL_COUNTS_PRINTED = (1, 8, 34, 98, 205, 302, 275, 121, 24, 4, 0)
SIZE9_FAILURES_PRINTED = {(0, 0, 0, 0, 0, 1, 0, 1, 7),
                          (0, 0, 0, 0, 1, 0, 0, 0, 8),
                          (0, 0, 0, 0, 1, 0, 0, 1, 7),
                          (0, 0, 0, 0, 1, 1, 0, 0, 7)}
PI_ALL_ROOTS_PRINTED = {(1, 1): 10, (1, 2): 9, (1, 3): 10,
                        (2, 1): 9, (2, 2): 9, (2, 3): 9,
                        (3, 1): 10, (3, 2): 9, (3, 3): 10}
VIOLATING_ROOTS_PRINTED = [(1, 1), (1, 3), (3, 1), (3, 3)]
BUDGET_PRINTED = Fr(1, 9)                      # RESIDUAL + TOTAL LOSS, from VAL(D) = 10/9
FACTOR_LOSS_PRINTED = {(1, 2): Fr(3, 2), (1, 3): Fr(8, 3), (2, 1): Fr(0),
                       (2, 3): Fr(2, 3), (3, 1): Fr(0), (3, 2): Fr(1, 6)}
W_G243 = {(1, 2): 2, (1, 3): 4, (2, 3): 3}     # the second weighted triangle of Section 3
W_G223 = {(1, 2): 2, (1, 3): 2, (2, 3): 3}     # the same total weight 7 as W_G, and it holds
W_K4 = {(1, 2): 2, (3, 4): 2, (1, 3): 5, (1, 4): 5, (2, 3): 5, (2, 4): 5}   # source paper, lines 658-665

_ran = 0
_bad = 0


def check(name, ok, detail=''):
    global _ran, _bad
    _ran += 1
    if ok:
        print('PASS %s %s' % (name, detail))
    else:
        print('FAIL %s %s' % (name, detail))
        _bad += 1


# ---------------------------------------------------------------------------
# WEIGHTED GRAPHS, THE MOVE, AND pi
# ---------------------------------------------------------------------------
def mkgraph(n, wdict):
    """Symmetric weight matrix on 0..n-1 from a 1-indexed dict of undirected edge weights."""
    W = [[0] * n for _ in range(n)]
    for (u, v), w in wdict.items():
        W[u - 1][v - 1] = w
        W[v - 1][u - 1] = w
    return W, n


def product(WG, nG, WH, nH):
    """The paper's weighted Cartesian product: vertex (i,j) -> index i*nH + j; a G-direction edge keeps its
    G-weight, an H-direction edge its H-weight."""
    n = nG * nH
    W = [[0] * n for _ in range(n)]
    for i in range(nG):
        for j in range(nH):
            u = i * nH + j
            for i2 in range(nG):
                if i2 != i and WG[i][i2]:
                    W[u][i2 * nH + j] = WG[i][i2]
            for j2 in range(nH):
                if j2 != j and WH[j][j2]:
                    W[u][i * nH + j2] = WH[j][j2]
    return W, n


def distances(W, n, root):
    """d(v) = min over paths v -> root of the product of the weights; d(root) = 1."""
    d = [None] * n
    d[root] = 1
    changed = True
    while changed:
        changed = False
        for u in range(n):
            if d[u] is None:
                continue
            for v in range(n):
                if W[u][v] and (d[v] is None or d[u] * W[u][v] < d[v]):
                    d[v] = d[u] * W[u][v]
                    changed = True
    return d


def is_tight(W, n):
    """w(u,v) = min over paths from u to v of the product of the weights, for every edge -- the source
    paper's normal form (its lines 640-646)."""
    for u in range(n):
        du = distances(W, n, u)
        for v in range(n):
            if u != v and W[u][v] and W[u][v] != du[v]:
                return False
    return True


def census(W, n, root):
    """Exhaustive backward closure over the box {D : 0 <= D(v) < d(v)}.

    Two facts make the box legitimate and are checked separately below: (a) if D(v) >= d(v) for some v then
    a single stack at v reaches the root along a minimising path, so D is solvable; (b) solvability is
    monotone in D.  Hence every unsolvable D lies in the box, and pi = 1 + (largest unsolvable size).

    -> (pi, sorted list of the largest unsolvable distributions, d, {size: count})
    """
    d = distances(W, n, root)
    moves = [(u, v, W[u][v], d[v] - 1) for u in range(n) for v in range(n) if W[u][v]]
    memo = {}

    def solvable(D):
        """D is inside the box, so D(root) = 0 (d(root) = 1) and no coordinate is saturated."""
        r = memo.get(D)
        if r is not None:
            return r
        memo[D] = False                      # a state cannot help solve itself
        res = False
        for (u, v, w, cap) in moves:
            if D[u] >= w:
                if D[v] >= cap:              # the move leaves the box: target reaches d(v), so solvable
                    res = True
                    break
                E = list(D)
                E[u] -= w
                E[v] += 1
                if solvable(tuple(E)):
                    res = True
                    break
        memo[D] = res
        return res

    bysize = {}
    for D in itertools.product(*[range(x) for x in d]):
        if not solvable(D):
            bysize.setdefault(sum(D), []).append(D)
    worst = max(bysize)
    counts = dict((s, len(v)) for s, v in bysize.items())
    return worst + 1, sorted(bysize[worst]), d, counts


def pi_rooted(W, n, root):
    return census(W, n, root)[0]


def pi_unrooted(W, n):
    """pi(G) = max over v of pi(G,v) -- the source paper's Proposition, its lines 215-216."""
    return max(pi_rooted(W, n, r) for r in range(n))


def reachable_forward(W, n, root, D):
    """Plain forward search: enumerate EVERY state reachable from D and report whether any of them carries a
    pebble on the root.  No distances, no box, no downward closure, no potential -- an independent path to
    the same conclusion.  -> (root_reached, number of states explored)"""
    seen = {tuple(D)}
    stack = [tuple(D)]
    while stack:
        S = stack.pop()
        if S[root] >= 1:
            return True, len(seen)
        for u in range(n):
            if S[u] < 2:
                continue
            for v in range(n):
                w = W[u][v]
                if w and S[u] >= w:
                    E = list(S)
                    E[u] -= w
                    E[v] += 1
                    T = tuple(E)
                    if T not in seen:
                        seen.add(T)
                        stack.append(T)
    return False, len(seen)


def cover_reachable(W, n, D):
    """Cover pebbling: can a state with at least one pebble on EVERY vertex be reached from D?"""
    seen = {tuple(D)}
    stack = [tuple(D)]
    while stack:
        S = stack.pop()
        if all(x >= 1 for x in S):
            return True
        for u in range(n):
            for v in range(n):
                w = W[u][v]
                if w and S[u] >= w:
                    E = list(S)
                    E[u] -= w
                    E[v] += 1
                    T = tuple(E)
                    if T not in seen:
                        seen.add(T)
                        stack.append(T)
    return False


def g6_decode(s):
    """graph6 -> (n, set of edges), for n < 63.  Written out so the program needs no graph library."""
    b = [ord(c) - 63 for c in s]
    n = b[0]
    bits = ''.join(format(x, '06b') for x in b[1:])
    E = set()
    k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k] == '1':
                E.add((i, j))
            k += 1
    return n, E


# ---------------------------------------------------------------------------
# 1. THE FACTOR
# ---------------------------------------------------------------------------
print('--- 1. the weighted triangle G: w(1,2)=2, w(1,3)=3, w(2,3)=2 ---')
WG, nG = mkgraph(3, W_G)
check('factor-is-tight', is_tight(WG, nG),
      '[every edge weight equals the least path product: 2 <= 3*2, 3 <= 2*2, 2 <= 2*3]')

d1 = distances(WG, nG, 0)
check('factor-distances-root-1', d1 == [1, 2, 3], '[d = %s]' % d1)

pf, argf = [], {}
for r in (0, 1, 2):
    p, worst, dd, _ = census(WG, nG, r)
    pf.append(p)
    argf[r + 1] = worst
    check('pi-factor-root-%d' % (r + 1),
          p == PI_FACTOR_PRINTED[r] and worst == sorted(LARGEST_UNSOLVABLE_PRINTED[r + 1]),
          '[pi = %d, d = %s, every largest unsolvable distribution: %s]'
          % (p, dd, ' '.join(str(x) for x in worst)))

check('pi-factor-unrooted', pi_unrooted(WG, nG) == 3,
      '[pi(G) = max_v pi(G,v) = 3, so the unrooted right-hand side is 3*3 = 9]')
check('rhs-rooted', pf[0] * pf[0] == 9, '[pi(G,1) pi(G,1) = 3*3 = 9]')

# ---------------------------------------------------------------------------
# 2. THE PRODUCT
# ---------------------------------------------------------------------------
print('--- 2. the product G x G ---')
WP, nP = product(WG, nG, WG, nG)
edges = [(u, v, WP[u][v]) for u in range(nP) for v in range(u + 1, nP) if WP[u][v]]
check('product-order-and-size', (nP, len(edges)) == (9, 18), '[9 vertices, 18 edges]')

want_w = {}
for i in range(3):
    for j in range(3):
        for i2 in range(3):
            if i2 != i:
                want_w[tuple(sorted((i * 3 + j, i2 * 3 + j)))] = W_G[tuple(sorted((i + 1, i2 + 1)))]
        for j2 in range(3):
            if j2 != j:
                want_w[tuple(sorted((i * 3 + j, i * 3 + j2)))] = W_G[tuple(sorted((j + 1, j2 + 1)))]
check('product-weights-match-the-definition',
      dict(((u, v), w) for (u, v, w) in edges) == want_w,
      '[each G-direction edge carries its G-weight, each H-direction edge its H-weight]')

deg = [sum(1 for v in range(nP) if WP[u][v]) for u in range(nP)]
check('product-skeleton-is-4-regular', deg == [4] * 9, '[the 3x3 rook graph K_3 box K_3]')

n6, E6 = g6_decode(G6_SKELETON)
check('graph6-decodes-to-the-skeleton',
      n6 == 9 and E6 == set((u, v) for (u, v, _) in edges),
      '[%r -> 9 vertices, %d edges, identical to the product skeleton]' % (G6_SKELETON, len(E6)))

dP = distances(WP, nP, 0)
check('distance-vector-as-printed', tuple(dP) == D_VEC_PRINTED, '[d = %s]' % (tuple(dP),))
check('distance-product-formula',
      all(dP[i * 3 + j] == d1[i] * d1[j] for i in range(3) for j in range(3)),
      '[d((i,j)) = d(i) d(j) at all nine vertices]')

# ---------------------------------------------------------------------------
# 3. THE WITNESS
# ---------------------------------------------------------------------------
print('--- 3. the distribution D: 7 on (3,3), one on (2,3), one on (3,2) ---')
check('witness-size', sum(D_PRINTED) == 9, '[|D| = 9 = pi(G,1) pi(G,1)]')
check('witness-inside-the-box', all(D_PRINTED[v] < dP[v] for v in range(nP)),
      '[D(v) < d(v) at every vertex, so no single stack reaches the root by itself]')

hit, states = reachable_forward(WP, nP, 0, D_PRINTED)
check('witness-unsolvable-by-plain-forward-search', hit is False,
      '[%d states reachable from D, none of them carries a pebble on (1,1)]' % states)

piP, worstP, _, countsP = census(WP, nP, 0)
check('pi-product-at-root-(1,1)', piP == PI_PRODUCT_ROOT_PRINTED,
      '[pi(G x G,(1,1)) = %d, over a box of %d states]' % (piP, 46656))
check('violation-rooted', piP > pf[0] * pf[0],
      '[%d > %d = pi(G,1) pi(G,1): the rooted conjecture fails]' % (piP, pf[0] * pf[0]))

lvl = tuple(countsP.get(s, 0) for s in range(11))
check('unsolvable-level-counts', lvl == LEVEL_COUNTS_PRINTED, '[by size: %s]' % (lvl,))
check('exactly-four-failures-of-size-9',
      set(worstP) == SIZE9_FAILURES_PRINTED and len(worstP) == 4,
      '[4 of the C(17,8) = 24310 distributions of 9 pebbles fail, and D is one of them]')

swap = lambda D: tuple(D[(k % 3) * 3 + k // 3] for k in range(9))
fixed = [D for D in sorted(worstP) if swap(D) == D]
moved = [D for D in sorted(worstP) if swap(D) != D]
check('size-9-failures-under-the-factor-swap',
      set(map(swap, worstP)) == set(worstP) and len(fixed) == 2 and len(moved) == 2 and
      swap(moved[0]) == moved[1] and sorted(worstP)[:2] == fixed,
      '[exchanging the two factors permutes the four failures: it fixes the first two, %s and %s, and '
      'exchanges the last two, so there are three up to that symmetry]' % (fixed[0], fixed[1]))

allroots = {}
for r in range(nP):
    allroots[(r // 3 + 1, r % 3 + 1)] = pi_rooted(WP, nP, r)
check('pi-product-at-all-nine-roots', allroots == PI_ALL_ROOTS_PRINTED, '[%s]' % (
    ', '.join('%s:%d' % (k, allroots[k]) for k in sorted(allroots))))
viol = sorted(k for k in allroots
              if allroots[k] > pf[k[0] - 1] * pf[k[1] - 1])
check('violating-roots', viol == VIOLATING_ROOTS_PRINTED, '[%s, all at 10 against 9]' % (viol,))
check('pi-product-unrooted', max(allroots.values()) == 10,
      '[pi(G x G) = 10 > 9 = pi(G) pi(G): the unrooted conjecture fails too]')

# ---------------------------------------------------------------------------
# 4. THE HAND PROOF, MECHANISED
# ---------------------------------------------------------------------------
print('--- 4. the potential argument of Section 2, step by step ---')
phi = [Fr(1, dP[v]) for v in range(nP)]
dirmoves = [(u, v, WP[u][v]) for u in range(nP) for v in range(nP) if WP[u][v]]
check('potential-is-valid',
      all(phi[v] <= w * phi[u] for (u, v, w) in dirmoves),
      '[phi = 1/d satisfies phi(u\') <= w phi(u) on all %d directed moves, so VAL never rises]'
      % len(dirmoves))

tight_chain = True
for v in range(nP):
    cur, ok = v, (v == 0)
    for _ in range(nP):
        if cur == 0:
            ok = True
            break
        nxt = [u for u in range(nP) if WP[cur][u] and dP[cur] == WP[cur][u] * dP[u]]
        if not nxt:
            break
        cur = nxt[0]
    tight_chain = tight_chain and ok
check('potential-is-pointwise-minimal', tight_chain,
      '[every vertex has a root path on which phi(u) = w phi(u\') is tight, so any valid potential '
      'normalised at the root has phi >= 1/d: 1/d is the BEST such certificate]')

VAL = sum(D_PRINTED[v] * phi[v] for v in range(nP))
check('val-of-the-witness', VAL == Fr(10, 9), '[VAL(D) = 7/9 + 1/6 + 1/6 = %s > 1]' % VAL)
check('potential-alone-does-not-suffice', VAL > 1,
      '[the obstruction is integrality, not the potential: VAL(D) = 10/9 exceeds 1]')
check('budget-identity', VAL - 1 == BUDGET_PRINTED,
      '[stopping when the root is first reached: RESIDUAL + TOTAL LOSS = 10/9 - 1 = 1/9]')

floss = {}
for a in (1, 2, 3):
    for b in (1, 2, 3):
        if a != b:
            floss[(a, b)] = Fr(W_G[tuple(sorted((a, b)))], d1[a - 1]) - Fr(1, d1[b - 1])
check('factor-loss-table', floss == FACTOR_LOSS_PRINTED,
      '[1->2: 3/2, 1->3: 8/3, 2->1: 0, 2->3: 2/3, 3->1: 0, 3->2: 1/6]')

loss = dict(((u, v), w * phi[u] - phi[v]) for (u, v, w) in dirmoves)
check('product-loss-is-factor-loss-over-the-other-coordinate',
      all(loss[(i * 3 + j, i2 * 3 + j)] == floss[(i + 1, i2 + 1)] / d1[j] and
          loss[(i * 3 + j, i * 3 + j2)] == floss[(j + 1, j2 + 1)] / d1[i]
          for i in range(3) for j in range(3) for i2 in range(3) for j2 in range(3)
          if i2 != i and j2 != j),
      '[LOSS((i,j)->(i\',j)) = L(i->i\')/d(j), and symmetrically]')

check('minimum-off-root-potential', min(phi[1:]) == Fr(1, 9),
      '[phi >= 1/9 off the root, so RESIDUAL is 0 or at least 1/9: either TOTAL LOSS = 0 (case B) or '
      'RESIDUAL = 0 and TOTAL LOSS = 1/9 (case A)]')

admissible = sorted(k for k in loss if loss[k] <= BUDGET_PRINTED)
lossless = sorted(k for k in loss if loss[k] == 0)
lossy_adm = sorted(k for k in admissible if loss[k] > 0)
want_lossless = sorted([(i * 3 + j, j) for i in (1, 2) for j in range(3)] +
                       [(i * 3 + j, i * 3) for j in (1, 2) for i in range(3)])
check('admissible-move-set',
      len(admissible) == 16 and lossless == want_lossless and
      [(k, loss[k]) for k in lossy_adm] == [((5, 4), Fr(1, 12)), ((7, 4), Fr(1, 12)),
                                            ((8, 5), Fr(1, 18)), ((8, 7), Fr(1, 18))],
      '[16 of the %d directed moves cost at most the whole budget 1/9: the 12 lossless ones (a coordinate '
      'dropping to 1), (2,3)->(2,2) and (3,2)->(2,2) at 1/12, and (3,3)->(2,3) and (3,3)->(3,2) at 1/18]'
      % len(dirmoves))

vals = sorted(set(loss.values()) - {Fr(0)})
sums = [m for k in range(1, 6)
        for m in itertools.combinations_with_replacement(vals, k)
        if sum(m) == BUDGET_PRINTED]
check('only-two-eighteenths-fill-the-budget', sums == [(Fr(1, 18), Fr(1, 18))],
      '[over every multiset of positive losses %s, the budget 1/9 is filled in exactly one way: two 1/18 '
      'moves, both out of (3,3). 1/12 = 1.5/18 leaves an unachievable half, so case A admits no 1/12 '
      'move and case B no lossy move at all]' % ('{' + ', '.join(str(x) for x in vals) + '}'))

check('caseB-strands-(3,2)',
      not any(k[1] == 7 for k in lossless) and min(WP[7][v] for v in range(nP) if WP[7][v]) == 2
      and phi[7] == Fr(1, 6) > BUDGET_PRINTED,
      '[with TOTAL LOSS = 0 only the 12 lossless moves occur; none lands on (3,2), and its single pebble '
      'cannot pay for any move (all weights are at least 2), so RESIDUAL >= phi(3,2) = 1/6 > 1/9]')

CASEA = [k for k in loss if loss[k] in (Fr(0), Fr(1, 18))]   # the 14 moves case A can use
check('caseA-move-set', len(CASEA) == 14 and sorted(k for k in admissible if k not in CASEA)
      == [(5, 4), (7, 4)],
      '[case A drops the two 1/12 moves from the 16, leaving 14]')

into33 = [k for k in loss if k[1] == 8]
check('caseA-vertex-(3,3)-never-gains',
      all(loss[k] > BUDGET_PRINTED for k in into33),
      '[every move landing on (3,3) costs at least %s > 1/9]' % min(loss[k] for k in into33))
out33 = dict((k, WP[8][k[1]]) for k in CASEA if k[0] == 8)
check('caseA-emptying-(3,3)',
      sorted(set(out33.values())) == [2, 3] and
      [(a, b) for a in range(4) for b in range(4) if 3 * a + 2 * b == 7 and b == 2] == [(1, 2)],
      '[its 7 pebbles leave in weight-3 and weight-2 moves only; with exactly two weight-2 moves, '
      '3a + 4 = 7 forces a = 1]')

into23 = [k for k in CASEA if k[1] == 5]
out23 = sorted(set(WP[5][k[1]] for k in CASEA if k[0] == 5))
check('caseA-forces-one-lossy-move-each-way',
      into23 == [(8, 5)] and out23 == [2, 3] and
      not any(2 * a + 3 * b == 1 for a in range(3) for b in range(3)),
      '[(2,3) gains only from (3,3), and must end empty; 1 pebble cannot be spent by any move, so it '
      'must receive one, and symmetrically for (3,2): the two lossy moves split 1 and 1]')

out13 = dict((k[1], WP[2][k[1]]) for k in CASEA if k[0] == 2)
check('caseA-endgame-is-stuck',
      out13 == {0: 3} and 3 * Fr(1, 3) == 1,
      '[from (1,3) the only admissible move is the weight-3 move to the root; the three surviving '
      'pebbles split 2 + 1 across (1,3) and (3,1), so none can move and RESIDUAL = 1 > 0]')

adm_by_u = {}
for (u, v) in admissible:
    adm_by_u.setdefault(u, []).append((v, WP[u][v]))
seen = {D_PRINTED}
stack = [D_PRINTED]
root_hit = False
while stack and not root_hit:
    S = stack.pop()
    if S[0] >= 1:
        root_hit = True
        break
    for u, outs in adm_by_u.items():
        for (v, w) in outs:
            if S[u] >= w:
                E = list(S)
                E[u] -= w
                E[v] += 1
                T = tuple(E)
                if T not in seen:
                    seen.add(T)
                    stack.append(T)
check('hand-proof-closure', root_hit is False,
      '[every move of a successful play has loss <= the total loss, hence <= 1/9, hence is one of the 16; '
      'exhausting the %d states reachable from D by those 16 moves never reaches (1,1). This is the whole '
      'hand proof in one search, independent of the case split above]' % len(seen))
check('hand-proof-gives-the-bound', not hit and sum(D_PRINTED) == 9 and 10 > pf[0] * pf[0],
      '[an unsolvable distribution of 9 pebbles gives pi(G x G,(1,1)) >= 10 > 9 without any census]')

# ---------------------------------------------------------------------------
# 5. CONTROLS
# ---------------------------------------------------------------------------
print('--- 5. controls: the semantics reproduce values published by other authors ---')


def unweighted(n, E):
    W = [[0] * n for _ in range(n)]
    for (u, v) in E:
        W[u][v] = 2
        W[v][u] = 2
    return W, n


cyc = lambda m: unweighted(m, [(i, (i + 1) % m) for i in range(m)])
kn = lambda m: unweighted(m, [(i, j) for i in range(m) for j in range(i + 1, m)])
cube = unweighted(8, [(x, x ^ (1 << b)) for x in range(8) for b in range(3) if x < x ^ (1 << b)])
pete = unweighted(10, [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] +
                  [(5 + i, 5 + (i + 2) % 5) for i in range(5)])
for name, (W, n), want in (('C_5', cyc(5), 5), ('C_6', cyc(6), 8), ('Q_3', cube, 8),
                           ('K_6', kn(6), 6), ('Petersen', pete, 10)):
    got = pi_unrooted(W, n)
    check('control-published-pi-%s' % name, got == want,
          '[all weights 2: pi = %d, published value %d]' % (got, want))

bad = []
for a in range(1, 13):
    for b in range(1, 13):
        WA, _ = mkgraph(2, {(1, 2): a})
        WB, _ = mkgraph(2, {(1, 2): b})
        W, n = product(WA, 2, WB, 2)
        if pi_rooted(W, n, 0) != a * b:
            bad.append((a, b))
check('control-chung-theorem-3-grid', bad == [],
      '[pi((K_2,a) x (K_2,b), corner) = ab on all 144 pairs with a,b <= 12, 0 failures]')

WA, _ = mkgraph(2, {(1, 2): 3})
WB, _ = mkgraph(2, {(1, 2): 2})
W, n = product(WA, 2, WB, 2)
check('control-weight-blindness-discriminator', pi_rooted(W, n, 0) == 6,
      '[pi((K_2,3) x (K_2,2)) = 6; a weight-blind engine would print 4]')

WK4, n4 = mkgraph(4, W_K4)
check('control-source-paper-K4-cover-13',
      all(cover_reachable(WK4, n4, [13 if k == s else 0 for k in range(4)]) for s in range(4)),
      '[the source paper\'s own example, its lines 658-665: thirteen pebbles on any single vertex cover '
      'the weighted K_4 with w(x1,x2) = w(x3,x4) = 2 and every other weight 5]')
check('control-source-paper-K4-cover-9-plus-4',
      not cover_reachable(WK4, n4, [9, 4, 0, 0]),
      '[nine on x_1 and four on x_2 do NOT cover it -- the paper\'s printed negative]')
check('control-source-paper-K4-cover-12-fails',
      not cover_reachable(WK4, n4, [12, 0, 0, 0]),
      '[twelve on x_1 do not cover, so the threshold is exactly the printed thirteen]')

check('outside-chungs-proved-region',
      d1[2] < W_G[(2, 3)] * d1[1] and d1[1] < W_G[(2, 3)] * d1[2],
      '[in a weighted power of K_2 with parallel weights equal every edge is geodesic to the root; here '
      'd(3) = %d < %d = w(2,3) d(2) and d(2) = %d < %d = w(2,3) d(3), so the edge {2,3} is geodesic in '
      'neither direction and G is not such a cube]'
      % (d1[2], W_G[(2, 3)] * d1[1], d1[1], W_G[(2, 3)] * d1[2]))

print('--- 6. anti-controls: the same code must NOT refute these ---')
for label, wd in (('unweighted-K3xK3', {(1, 2): 2, (1, 3): 2, (2, 3): 2}),
                  ('K3(2,2,3)-same-total-weight-7', W_G223)):
    WA, na = mkgraph(3, wd)
    W, n = product(WA, na, WA, na)
    p0 = pi_rooted(WA, na, 0)
    got = pi_rooted(W, n, 0)
    check('anti-control-%s' % label, got == p0 * p0 == 9,
          '[pi = %d = %d * %d: no violation, so the phenomenon is weight ASYMMETRY at the root and not '
          'total weight]' % (got, p0, p0))

# ---------------------------------------------------------------------------
# 6. THE TWO FURTHER WITNESSES
# ---------------------------------------------------------------------------
print('--- 7. the two further witnesses of Section 3 ---')
W243, n243 = mkgraph(3, W_G243)
p243 = [pi_rooted(W243, n243, r) for r in range(3)]
check('second-triangle-pi-at-its-three-roots', p243 == [4, 4, 5],
      '[K_3(2,4,3): pi = 4, 4, 5, so pi(K_3(2,4,3)) = 5 and vertex 1 is NOT a worst root]')

W2, m2 = product(WG, nG, W243, n243)
r2 = pi_rooted(W2, m2, 0)
u2 = pi_unrooted(W2, m2)
check('witness-2-rooted', r2 == 13 and r2 > pf[0] * p243[0],
      '[pi(K_3(2,3,2) x K_3(2,4,3),(1,1)) = 13 > 12 = 3*4]')
check('witness-2-unrooted-not-violated', u2 == 14 and u2 <= 3 * 5,
      '[its unrooted value is 14, at or below 3*5 = 15: no unrooted violation]')

W3, m3 = product(W243, n243, W243, n243)
r3 = pi_rooted(W3, m3, 0)
u3 = pi_unrooted(W3, m3)
check('witness-3-rooted', r3 == 18 and r3 > p243[0] * p243[0],
      '[pi(K_3(2,4,3)^2,(1,1)) = 18 > 16 = 4*4]')
check('witness-3-unrooted-not-violated', u3 == 19 and u3 <= 5 * 5,
      '[its unrooted value is 19, well below 5*5 = 25: no unrooted violation]')

# ---------------------------------------------------------------------------
# SCOPE, then the verdict
# ---------------------------------------------------------------------------
print()
print('NOT RE-RUN: the search that FOUND this witness -- a census of 47,336 tight weighted instances over '
      'six cells -- is not repeated here, and no minimality claim of any kind is made or checked. The '
      'smallest counterexample, over weighted graphs or even over weighted triangles, is not determined '
      'by this program.')
print('NOT RE-RUN: the printed text of Discrete Math. 312 (2012) 2286-2293 is paywalled and was never '
      'read, so the printed CONJECTURE NUMBERS are unverified; the paper quotes the statements from the '
      'e-print source by LaTeX label and line range instead, and this program checks mathematics only, '
      'never bibliography.')
print('NOT RE-RUN: nothing here bears on Graham\'s (unweighted) conjecture, on the unweighted rooted '
      'conjecture, or on Chung\'s proved region beyond the 144-instance grid check above.')
print()
if _bad:
    print('%d CHECK(S) FAILED of %d' % (_bad, _ran))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _ran)
sys.exit(0)

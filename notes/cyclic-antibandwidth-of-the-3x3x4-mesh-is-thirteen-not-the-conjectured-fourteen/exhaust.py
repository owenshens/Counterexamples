#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exhaust.py -- the most conservative possible re-run of the decisive negative:

    there is NO bijection f : V(P_3 x P_3 x P_4) -> Z_36 with cyclic edge distance >= 14,

carried out with **no symmetry reduction of any kind** beyond pinning one vertex to label 0 by
the label-rotation argument (which is a proof, not a heuristic: f -> f + t is a bijection of the
solution set and acts transitively on which vertex carries label 0).

WHY THIS FILE EXISTS SEPARATELY FROM verify.py.  verify.py decides the same question in about
three minutes using a 16-fold symmetry reduction (eight automorphisms of the 3x3 grid times
label reflection), and it proves that reduction sound on the paper's own witness.  The reduction
is nevertheless the only step of the negative that requires an argument, so this program removes
it entirely and pays the full price -- roughly sixteen times the tree, tens of minutes.  If the
two runs disagree, the symmetry argument is wrong; they do not.

Python 3.9+, STANDARD LIBRARY ONLY.  Exact integer arithmetic throughout; no randomness, no
network, no external data file.  Deterministic: the node counts below are reproducible exactly.
"""

import itertools
import sys
import time

N_PASS = 0
N_FAIL = 0


def check(name, ok, detail=''):
    global N_PASS, N_FAIL
    if ok:
        N_PASS += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        N_FAIL += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))
    sys.stdout.flush()


def mesh(dims):
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


def dc(a, b, n):
    d = abs(a - b) % n
    return min(d, n - d)


def exhaust(dims, k, heartbeat=1 << 26):
    """Complete exhaustion, NO symmetry reduction beyond the rotation pin.
    Returns (status, nodes, leaves, depth_profile, seconds)."""
    verts, idx, E = mesh(dims)
    n = len(verts)
    ALL = (1 << n) - 1
    nbrs = [[] for _ in range(n)]
    for u, v in E:
        nbrs[u].append(v)
        nbrs[v].append(u)
    far = []
    for l1 in range(n):
        m = 0
        for l2 in range(n):
            if dc(l1, l2, n) >= k:
                m |= 1 << l2
        far.append(m)
    deg = [len(nbrs[v]) for v in range(n)]
    start = max(range(n), key=lambda v: (deg[v], -v))
    order = [start]
    ins = {start}
    while len(order) < n:
        best = max((v for v in range(n) if v not in ins),
                   key=lambda v: (sum(1 for w in nbrs[v] if w in ins), deg[v], -v))
        order.append(best)
        ins.add(best)
    pos = {v: i for i, v in enumerate(order)}
    future = [[w for w in nbrs[order[i]] if pos[w] > i] for i in range(n)]
    cand = [ALL] * n
    prof = [0] * (n + 1)
    leaves = [0]
    st = {'nodes': 0}
    t0 = time.time()

    def rec(i, used):
        if i == n:
            leaves[0] += 1
            return True
        v = order[i]
        c = cand[v] & ~used
        while c:
            b = c & -c
            c ^= b
            l = b.bit_length() - 1
            prof[i] += 1
            st['nodes'] += 1
            if st['nodes'] % heartbeat == 0:
                print('    ... %d nodes, %.0f s' % (st['nodes'], time.time() - t0))
                sys.stdout.flush()
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
            r = rec(i + 1, nused)
            for x, o in saved:
                cand[x] = o
        return False

    sys.setrecursionlimit(10000)
    for w in future[0]:
        cand[w] &= far[0]
    rec(1, 1)
    return ('FEASIBLE' if leaves[0] else 'INFEASIBLE'), st['nodes'], leaves[0], prof, time.time() - t0


DIMS = (3, 3, 4)
verts, idx, E = mesh(DIMS)
n = len(verts)
print('P_3 x P_3 x P_4 : |V| = %d, |E| = %d' % (n, len(E)))
print('deciding cab >= 14 by complete exhaustion, NO symmetry reduction (rotation pin only)')
print('')
sys.stdout.flush()

st, nodes, leaves, prof, el = exhaust(DIMS, 14)
print('')
print('    RESULT %s ; search-tree nodes = %d ; complete labelings found = %d ; %.1f s'
      % (st, nodes, leaves, el))
print('    DEPTH PROFILE (nodes created while placing the i-th vertex of the order, i = 1..35):')
print('    ' + ' '.join('i%d=%d' % (i, prof[i]) for i in range(1, n)))
print('')
sys.stdout.flush()

check('k-14-is-infeasible-with-no-symmetry-reduction-at-all',
      st == 'INFEASIBLE' and leaves == 0,
      '%d search-tree nodes exhausted, 0 complete labelings, %.1f s' % (nodes, el))
check('depth-profile-sums-to-the-reported-node-count', sum(prof) == nodes,
      'sum over the 35 placement levels = %d' % nodes)
deep = max(i for i in range(1, n) if prof[i] > 0)
check('the-exhaustion-is-not-a-trivially-empty-one',
      nodes > 10 ** 8 and deep >= 20
      and all(prof[i] > 0 for i in range(1, deep + 1))
      and all(prof[i] == 0 for i in range(deep + 1, n)),
      'levels 1..%d of the tree are non-empty and every level beyond is empty: the search dies '
      '%d vertices short of a complete labeling, which is exactly what "no labeling exists" '
      'looks like from inside a backtracking search' % (deep, n - deep))
check('hence-cab-P_3x3x4-is-at-most-13',
      st == 'INFEASIBLE',
      'combined with the printed bijection of min cyclic distance 13, cab(P_3x3x4) = 13, against '
      'the conjectured 14')

print('')
print('NOT RE-RUN: everything else. This program decides one question on one graph. The lower '
      'bound, the degree cap, the k = 15 layer, the controls on published values and the '
      'soundness of the symmetry reduction are all in verify.py; the cells n3 = 5..8 are decided '
      'nowhere in this folder.')
print('')
if N_FAIL:
    print('VERDICT: %d CHECKS FAILED' % N_FAIL)
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % N_PASS)
sys.exit(0)

#!/usr/bin/env python3
"""verify.py -- re-derives every quantity claimed in

    "Girth-12 Lifts of the Voltage Graph G_12 for m = 6, 7, 8,
     and Nonexistence on Its Printed Arc Set for m = 4, 5"

from these inputs: the 51 tree edges of T_12 and the 24 arcs of the corrected
Table 3 with the authors' own symbol names, both transcribed from the note's
tables; the three voltage vectors exhibited in the note's Section 4; and a
second family of three voltage vectors supplied by this program, which the
note neither prints nor uses.

Python 3.9+, STANDARD LIBRARY ONLY -- no numpy/sympy/networkx, no external data
file, so a referee needs nothing installed. All arithmetic is exact integer
arithmetic in Z and in Z_m; no float ever takes part in a decision.

One `PASS <name> [detail]` line per check, closing with
    VERDICT: ALL <n> CHECKS PASS
and exit status 0 iff every check passed.
"""

import itertools
import platform
import sys
from collections import Counter, deque

_N = [0]
_F = [0]


def CHECK(name, cond, detail=''):
    if cond:
        _N[0] += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _F[0] += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def NOTE(text):
    print('NOTE %s' % text)


print('verify.py -- girth-12 lifts of G_12 at m=6,7,8 and nonexistence on its'
      ' printed arc set at m=4,5')
print('python %s, exact integer arithmetic only, standard library only'
      % platform.python_version())
print()

# ============================================================================
# 1. THE OBJECT, EXACTLY AS PRINTED IN THE PAPER
# ============================================================================
# The 51 tree edges of T_12, transcribed from the paper's Table 1.
TREE_PRINTED = [
    ('x*', 'x'),
    ('x', 'x0'), ('x', 'x1'),
    ('x0', 'x00'), ('x0', 'x01'),
    ('x1', 'x10'), ('x1', 'x11'),
    ('x01', 'x010'), ('x01', 'x011'),
    ('x10', 'x100'), ('x10', 'x101'),
    ('x11', 'x110'), ('x11', 'x111'),
    ('x010', 'x0100'), ('x010', 'x0101'),
    ('x011', 'x0110'), ('x011', 'x0111'),
    ('x100', 'x1000'), ('x100', 'x1001'),
    ('x101', 'x1010'), ('x101', 'x1011'),
    ('x110', 'x1100'), ('x110', 'x1101'),
    ('x111', 'x1110'), ('x111', 'x1111'),
    ('y*', 'y'),
    ('y', 'y0'), ('y', 'y1'),
    ('y0', 'y01'),
    ('y1', 'y10'), ('y1', 'y11'),
    ('y01', 'y010'), ('y01', 'y011'),
    ('y10', 'y100'), ('y10', 'y101'),
    ('y11', 'y110'), ('y11', 'y111'),
    ('z*', 'z'),
    ('z', 'z0'), ('z', 'z1'),
    ('z0', 'z01'),
    ('z1', 'z10'), ('z1', 'z11'),
    ('z01', 'z010'), ('z01', 'z011'),
    ('z10', 'z100'), ('z10', 'z101'),
    ('z11', 'z110'), ('z11', 'z111'),
    ('x00', 'y0'), ('x00', 'z0'),          # the two join edges
]

# The 24 arcs of the CORRECTED Table 3, in the authors' own symbol order
# a,b,c,d,e,f,g,h,j,k,l,p,q,r,s,t,u,v,w,alpha,beta,gamma,delta,epsilon
# (their list skips i, m, n, o).
ARCNAME = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'j', 'k', 'l', 'p',
           'q', 'r', 's', 't', 'u', 'v', 'w', 'alpha', 'beta', 'gamma',
           'delta', 'epsilon']
ARCS_PRINTED = [
    ('a', 'x0100', 'y010', 1), ('b', 'x0100', 'z110', 2),
    ('c', 'x0101', 'y011', 2), ('d', 'x0101', 'z111', 1),
    ('e', 'x0110', 'y100', 1), ('f', 'x0110', 'z100', -1),
    ('g', 'x0111', 'y101', 2), ('h', 'x0111', 'z101', 3),
    ('j', 'x1000', 'y010', 2), ('k', 'x1000', 'z010', -1),
    ('l', 'x1001', 'y011', 1), ('p', 'x1001', 'z011', 3),
    ('q', 'x1010', 'y110', 1), ('r', 'x1010', 'z100', 1),
    ('s', 'x1011', 'y111', 2), ('t', 'x1011', 'z101', -1),
    ('u', 'x1100', 'y100', 2), ('v', 'x1100', 'z010', -3),
    ('w', 'x1101', 'y101', 1), ('alpha', 'x1101', 'z011', -2),
    ('beta', 'x1110', 'y110', 3), ('gamma', 'x1110', 'z110', -1),
    ('delta', 'x1111', 'y111', -1), ('epsilon', 'x1111', 'z111', 2),
]
ROW_XY_PRINTED = [1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 3, -1]
ROW_XZ_PRINTED = [2, 1, -1, 3, -1, 3, 1, -1, -3, -2, -1, 2]

CHECK('arc_symbol_order_is_the_authors_own_list',
      [n for n, _, _, _ in ARCS_PRINTED] == ARCNAME,
      'a,b,c,...,epsilon with i,m,n,o skipped')


def derive_tree():
    """T_12 rebuilt from the source's PRUNING RULE rather than from the table.

    X'_3 : all bitstrings of length <= 4, deleting every string that begins
    "00" except "00" itself (the subtree below x_00 is pruned).
    Y_2 = Z_2 : all bitstrings of length <= 3, deleting every string beginning
    "00" (the pruned subtree hangs below y_a = y_0).
    """
    E = [('x*', 'x')]
    xs = ['']
    for L in range(1, 5):
        for bits in itertools.product('01', repeat=L):
            s = ''.join(bits)
            if s.startswith('00') and s != '00':
                continue
            xs.append(s)
    for s in xs:
        if s:
            E.append(('x' + s[:-1], 'x' + s))
    for pref in ('y', 'z'):
        ys = ['']
        for L in range(1, 4):
            for bits in itertools.product('01', repeat=L):
                s = ''.join(bits)
                if s.startswith('00'):
                    continue
                ys.append(s)
        E.append((pref + '*', pref))
        for s in ys:
            if s:
                E.append((pref + s[:-1], pref + s))
    E.append(('x00', 'y0'))
    E.append(('x00', 'z0'))
    return E


def canon(edges):
    return set(frozenset(e) for e in edges)


CHECK('printed_tree_equals_the_tree_derived_from_the_sources_pruning_rule',
      canon(TREE_PRINTED) == canon(derive_tree()),
      'two independent constructions of the same 51 edges')

VERT = []
for u, v in TREE_PRINTED:
    for w in (u, v):
        if w not in VERT:
            VERT.append(w)
IDX = {v: i for i, v in enumerate(VERT)}
PIN = ['x*', 'y*', 'z*']
ORD = [v for v in VERT if v not in PIN]

CHECK('order_of_T12_is_52', len(VERT) == 52, '49 ordinary + 3 pinned')
CHECK('size_of_T12_is_51',
      len(TREE_PRINTED) == 51 and len(canon(TREE_PRINTED)) == 51,
      'no repeated tree edge, so T_12 is a tree on 52 vertices')
CHECK('number_of_arcs_is_24', len(ARCS_PRINTED) == 24)

EDGES = [(u, v, 0) for (u, v) in TREE_PRINTED] + \
        [(u, v, w) for (_, u, v, w) in ARCS_PRINTED]
CHECK('size_of_G12_is_75',
      len(EDGES) == 75 and len(canon([(u, v) for (u, v, _) in EDGES])) == 75,
      '51 tree edges + 24 arcs, all distinct')
CHECK('cycle_rank_of_G12_is_24', len(EDGES) - len(VERT) + 1 == 24,
      'beta = 75 - 52 + 1 = 24')
CHECK('tree_edge_count_identity_48_plus_3_pinned_equals_51',
      len(ORD) == 49 and (len(ORD) - 1) + len(PIN) == len(TREE_PRINTED),
      'counts only: 49 ordinary vertices, so the subtree they span has 48 '
      'edges, and 48 + 3 pinned edges = 51 = |E(T_12)|; the standard '
      'gauge/T-normalisation argument that these 51 tree voltages may be '
      'taken to be 0 and that the 24 arc voltages then parametrise all lifts '
      'up to isomorphism (Gross-Tucker voltage-graph normalisation) is NOT '
      're-verified here')

DEG = Counter()
for u, v, _ in EDGES:
    DEG[u] += 1
    DEG[v] += 1
CHECK('degree_multiset_of_G12_is_three_1s_and_fortynine_3s',
      sorted(DEG.values()) == [1, 1, 1] + [3] * 49)
CHECK('the_three_degree_one_vertices_are_exactly_the_pinned_ones',
      sorted(v for v in VERT if DEG[v] == 1) == ['x*', 'y*', 'z*'])

XL = ['x0100', 'x0101', 'x0110', 'x0111', 'x1000', 'x1001',
      'x1010', 'x1011', 'x1100', 'x1101', 'x1110', 'x1111']
tails = Counter(u for (_, u, _, _) in ARCS_PRINTED)
heads = Counter(v for (_, _, v, _) in ARCS_PRINTED)
CHECK('each_of_the_twelve_x_leaves_carries_exactly_two_arcs',
      sorted(tails) == XL and set(tails.values()) == {2})
CHECK('each_x_leaf_has_one_y_arc_and_one_z_arc',
      all(sorted(h[0] for (_, u, h, _) in ARCS_PRINTED if u == xl) ==
          ['y', 'z'] for xl in XL))
CHECK('each_of_the_twelve_y_and_z_leaves_receives_exactly_two_arcs',
      len(heads) == 12 and set(heads.values()) == {2} and
      sorted(heads) == ['y010', 'y011', 'y100', 'y101', 'y110', 'y111',
                        'z010', 'z011', 'z100', 'z101', 'z110', 'z111'])
CHECK('the_two_printed_voltage_rows_are_reproduced_by_the_symbol_order',
      [w for (_, _, h, w) in ARCS_PRINTED if h[0] == 'y'] == ROW_XY_PRINTED and
      [w for (_, _, h, w) in ARCS_PRINTED if h[0] == 'z'] == ROW_XZ_PRINTED,
      'x->y %s ; x->z %s' % (ROW_XY_PRINTED, ROW_XZ_PRINTED))
CHECK('largest_printed_absolute_voltage_is_3',
      max(abs(w) for (_, _, _, w) in ARCS_PRINTED) == 3)

# the forced Table-3 typo: the printed 8th "Starting leaf" reads x_1001
LITDEG = Counter()
for (_, u, v, _) in ARCS_PRINTED:
    LITDEG['x1001' if u == 'x1011' else u] += 1
    LITDEG[v] += 1
for u, v in TREE_PRINTED:
    LITDEG[u] += 1
    LITDEG[v] += 1
CHECK('the_literal_printed_table_is_impossible_x1001_would_have_degree_5',
      LITDEG['x1001'] == 5 and LITDEG['x1011'] == 1,
      'so the duplicated 8th Starting-leaf entry must read x_1011, the only '
      '4-bit x-leaf label otherwise missing; with the correction every '
      'non-pinned vertex has degree 3')

ADJ = {v: [] for v in VERT}
for u, v, _ in EDGES:
    ADJ[u].append(v)
    ADJ[v].append(u)
col = {'x*': 0}
dq = deque(['x*'])
bip = True
while dq:
    u = dq.popleft()
    for w in ADJ[u]:
        if w not in col:
            col[w] = 1 - col[u]
            dq.append(w)
        elif col[w] == col[u]:
            bip = False
CHECK('G12_is_connected_and_bipartite', bip and len(col) == 52)

par, dep = {}, {'x*': 0}
tadj = {v: [] for v in VERT}
for u, v in TREE_PRINTED:
    tadj[u].append(v)
    tadj[v].append(u)
seen = {'x*'}
dq = deque(['x*'])
while dq:
    u = dq.popleft()
    for w in tadj[u]:
        if w not in seen:
            seen.add(w)
            par[w] = u
            dep[w] = dep[u] + 1
            dq.append(w)


def treepath(u, v):
    pu, pv = [u], [v]
    a, b = u, v
    while dep[a] > dep[b]:
        a = par[a]
        pu.append(a)
    while dep[b] > dep[a]:
        b = par[b]
        pv.append(b)
    while a != b:
        a = par[a]
        pu.append(a)
        b = par[b]
        pv.append(b)
    return pu + pv[::-1][1:]


def tdist(u, v):
    return len(treepath(u, v)) - 1


CHECK('all_three_pinned_to_pinned_tree_distances_are_6',
      tdist('x*', 'y*') == 6 and tdist('x*', 'z*') == 6 and
      tdist('y*', 'z*') == 6,
      'so every lift already contains a 12-cycle: "girth 12" here is exactly '
      '"no cycle shorter than 12", never "girth at least 12"')
CHECK('the_four_tree_distances_the_hand_proof_uses_are_9',
      tdist('x1000', 'z010') == 9 and tdist('x1001', 'z011') == 9 and
      tdist('x1100', 'z010') == 9 and tdist('x1101', 'z011') == 9,
      'six steps up to x_00 then three down, so each of k,p,v,alpha closes a '
      'fundamental cycle of length 1+9 = 10')

# ============================================================================
# 2. THE CYCLE CENSUS: every simple cycle of length <= 10, by two algorithms
# ============================================================================
NA = {v: [] for v in VERT}
for ei, (u, v, _) in enumerate(EDGES):
    NA[u].append((v, ei, +1))
    NA[v].append((u, ei, -1))
LMAX = 10
CYCLES = []


def dfs(start, cur, path, used, coeff):
    for (w, ei, sg) in NA[cur]:
        if ei in used:
            continue
        if w == start:
            if len(path) >= 3 and IDX[path[1]] < IDX[path[-1]]:
                c = dict(coeff)
                c[ei] = sg
                CYCLES.append((len(path), c))
            continue
        if IDX[w] <= IDX[start] or w in path or len(path) + 1 > LMAX:
            continue
        used.add(ei)
        coeff[ei] = sg
        path.append(w)
        dfs(start, w, path, used, coeff)
        path.pop()
        del coeff[ei]
        used.discard(ei)


sys.setrecursionlimit(10000)
for s in VERT:
    dfs(s, s, [s], set(), {})

hist = dict(sorted(Counter(L for L, _ in CYCLES).items()))
CHECK('undirected_simple_cycles_of_length_at_most_10_number_126',
      len(CYCLES) == 126, 'length histogram %s' % hist)
CHECK('cycle_length_histogram_is_twelve_6s_twentyfour_8s_ninety_10s',
      hist == {6: 12, 8: 24, 10: 90})
CHECK('twice_the_census_size_is_252',
      2 * len(CYCLES) == 252,
      'the source states "there are 252 short cycles"; so the authors count '
      'each cycle in both directions')
CHECK('the_census_contains_no_4_cycle',
      sum(1 for L, _ in CYCLES if L == 4) == 0)

TEI = {}
for ei, (u, v) in enumerate(TREE_PRINTED):
    TEI[(u, v)] = ei
    TEI[(v, u)] = ei
FUND = []
for j, (_, u, v, _) in enumerate(ARCS_PRINTED):
    p = treepath(u, v)
    mask = 1 << (51 + j)
    for i in range(len(p) - 1):
        mask |= 1 << TEI[(p[i], p[i + 1])]
    FUND.append(mask)
FL = [bin(m).count('1') for m in FUND]
EEND = [(u, v) for (u, v) in TREE_PRINTED] + \
       [(u, v) for (_, u, v, _) in ARCS_PRINTED]


def single_cycle_len(mask, n):
    d = Counter()
    ad = {}
    for ei in range(75):
        if mask >> ei & 1:
            u, v = EEND[ei]
            d[u] += 1
            d[v] += 1
            ad.setdefault(u, []).append(v)
            ad.setdefault(v, []).append(u)
    if any(c != 2 for c in d.values()):
        return 0
    st = [next(iter(d))]
    vis = set(st)
    while st:
        u = st.pop()
        for w in ad[u]:
            if w not in vis:
                vis.add(w)
                st.append(w)
    return n if len(vis) == len(d) else 0


FOUND2 = set()


def sub_xor(start, k, mask):
    """route 2: every element of the cycle space spanned by <= 10 fundamental
    cycles, kept when it happens to be a single simple cycle of length <= 10."""
    if k:
        pc = bin(mask).count('1')
        if 3 <= pc <= LMAX and single_cycle_len(mask, pc):
            FOUND2.add(mask)
    if k == LMAX:
        return
    for j in range(start, 24):
        sub_xor(j + 1, k + 1, mask ^ FUND[j])


sub_xor(0, 0, 0)
MASK1 = set()
for L, c in CYCLES:
    mm = 0
    for ei in c:
        mm |= 1 << ei
    MASK1.add(mm)
CHECK('a_second_census_by_cycle_space_XOR_returns_the_same_126_edge_sets',
      len(FOUND2) == 126 and FOUND2 == MASK1,
      'path-DFS with min-index canonicalisation, and XOR of fundamental '
      'cycles over every subset of the 24 co-tree arcs of size <= 10, agree')

FORCED = [ARCNAME[j] for j in range(24) if FL[j] <= 10]
CHECK('exactly_14_of_the_24_arcs_are_forced_nonzero',
      FORCED == ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'j', 'k', 'l', 'p',
                 'v', 'alpha'],
      'a and c have fundamental cycles of length 8, the other twelve of '
      'length 10')
CHECK('the_remaining_ten_arcs_have_fundamental_cycles_of_length_12',
      all(FL[j] == 12 for j in range(24) if ARCNAME[j] not in FORCED),
      'so q,r,s,t,u,w,beta,gamma,delta,epsilon may legitimately be 0, which '
      'is why zero entries appear in the exhibited witnesses')

CONS = []
for L, c in CYCLES:
    CONS.append({ei - 51: sg for ei, sg in c.items() if ei >= 51})
CHECK('every_short_cycle_carries_at_least_one_arc',
      all(len(d) >= 1 for d in CONS),
      'T_12 is a tree, so the criterion is 126 nontrivial linear forms')
supph = dict(sorted(Counter(len(d) for d in CONS).items()))
CHECK('support_size_histogram_of_the_126_forms_is_14_44_4_40_24',
      supph == {1: 14, 2: 44, 3: 4, 4: 40, 6: 24}, '%s' % supph)


def form_sum(d, vv):
    return sum(sg * vv[k] for k, sg in d.items())


sys.stdout.flush()

# ============================================================================
# 3. LIFTS AND GIRTH, BY TWO STRUCTURALLY DIFFERENT ALGORITHMS
# ============================================================================
ORDI = {u: k for k, u in enumerate(ORD)}


def build_lift(m, vv, with_arcs=True):
    n = 3 + m * len(ORD)
    adj = [[] for _ in range(n)]
    ec = 0
    for ei, (u, v, _) in enumerate(EDGES):
        if ei >= 51 and not with_arcs:
            continue
        w = 0 if ei < 51 else vv[ei - 51] % m
        if u in PIN or v in PIN:
            p, r = (u, v) if u in PIN else (v, u)
            pi = PIN.index(p)
            for i in range(m):
                a, b = pi, 3 + m * ORDI[r] + i
                adj[a].append(b)
                adj[b].append(a)
                ec += 1
        else:
            for i in range(m):
                a = 3 + m * ORDI[u] + i
                b = 3 + m * ORDI[v] + (i + w) % m
                adj[a].append(b)
                adj[b].append(a)
                ec += 1
    return n, adj, ec


def girth_bfs(n, adj):
    """route A: BFS from every vertex; min over non-tree edges is exact."""
    best = 10 ** 9
    for s in range(n):
        dist = [-1] * n
        pr = [-1] * n
        dist[s] = 0
        dq2 = deque([s])
        while dq2:
            u = dq2.popleft()
            if 2 * dist[u] >= best:
                break
            for w in adj[u]:
                if dist[w] == -1:
                    dist[w] = dist[u] + 1
                    pr[w] = u
                    dq2.append(w)
                elif pr[u] != w and pr[w] != u:
                    c = dist[u] + dist[w] + 1
                    if c < best:
                        best = c
    return best


def girth_edgedel(n, adj):
    """route B: delete each edge, BFS between its ends, girth = min(1+dist)."""
    best = 10 ** 9
    E2 = set()
    for a in range(n):
        for b in adj[a]:
            E2.add((min(a, b), max(a, b)))
    for (a, b) in E2:
        dist = [-1] * n
        dist[a] = 0
        dq2 = deque([a])
        hit = -1
        while dq2:
            u = dq2.popleft()
            if u == b:
                hit = dist[u]
                break
            if dist[u] + 1 >= best:
                continue
            for w in adj[u]:
                if (u == a and w == b) or (u == b and w == a):
                    continue
                if dist[w] == -1:
                    dist[w] = dist[u] + 1
                    dq2.append(w)
        if hit >= 0 and hit + 1 < best:
            best = hit + 1
    return best


def heawood():
    m = 7
    n = 2 * m
    adj = [[] for _ in range(n)]
    for w in (0, 1, 3):
        for i in range(m):
            adj[i].append(m + (i + w) % m)
            adj[m + (i + w) % m].append(i)
    return n, adj


n0, a0 = heawood()
CHECK('control0_the_heawood_graph_as_a_Z7_lift_is_cubic_order_14_girth_6',
      n0 == 14 and set(len(x) for x in a0) == {3} and
      girth_bfs(n0, a0) == 6 and girth_edgedel(n0, a0) == 6,
      'a published object with no input from this problem: the (3,6)-cage has '
      '14 vertices and girth 6, and both girth routes say 6')

ZERO = [0] * 24
c1 = [(m,) + tuple([build_lift(m, ZERO, with_arcs=False)[0],
                    girth_bfs(*build_lift(m, ZERO, with_arcs=False)[:2])])
      for m in (4, 5, 9)]
CHECK('control1_the_arc_free_lift_of_T12_has_girth_12_at_m_4_5_and_9',
      all(g == 12 for _, _, g in c1),
      'the source lemma girth(T_4t;m) = 4t = 12; the test must stay SILENT '
      'here, and does; orders %s' % [(m, n) for m, n, _ in c1])

PUB = [w for (_, _, _, w) in ARCS_PRINTED]
prof = {}
ordersize_ok = True
degrees_ok = True
for m in range(3, 15):
    n2, a2, ec = build_lift(m, PUB)
    prof[m] = girth_bfs(n2, a2)
    if n2 != 49 * m + 3 or ec != 75 * m:
        ordersize_ok = False
    dd = Counter(len(x) for x in a2)
    want = {3: 49 * 3 + 3} if m == 3 else {3: 49 * m, m: 3}
    if dd != want:
        degrees_ok = False
CHECK('published_lift_has_order_49m_plus_3_and_size_75m_for_m_3_to_14',
      ordersize_ok,
      'orders 150, 199, 248, 297, 346, 395, 444, ..., 689')
CHECK('published_lift_has_three_vertices_of_degree_m_and_49m_of_degree_3',
      degrees_ok, 'at m=3 the pinned vertices also have degree 3')
CHECK('published_girth_profile_is_6_6_8_8_8_10_then_12_from_m_9_on',
      [prof[m] for m in range(3, 15)] == [6, 6, 8, 8, 8, 10] + [12] * 6,
      'm=3..14 -> %s' % [prof[m] for m in range(3, 15)])
CHECK('printed_voltages_give_girth_12_for_exactly_m_9_to_14_in_range_3_to_14',
      all(prof[m] == 12 for m in range(9, 15)) and
      all(prof[m] < 12 for m in range(3, 9)),
      'with the PRINTED voltages the lift has girth 12 for exactly m >= 9, so '
      'the source theorem is TIGHT for the object it names')
CHECK('the_printed_assignment_gives_girth_8_not_12_at_m_6',
      prof[6] == 8,
      'so the witnesses below are DIFFERENT assignments on the SAME arc set')

sums = [form_sum(d, PUB) for d in CONS]
CHECK('published_absolute_voltage_sums_are_exactly_the_set_1_to_8',
      sorted(set(abs(s) for s in sums)) == list(range(1, 9)),
      'the source states they "all lie in the set {1,...,8}"; the set is '
      'attained exactly, multiset %s'
      % dict(sorted(Counter(abs(s) for s in sums).items())))
CHECK('no_short_cycle_of_the_published_assignment_sums_to_zero_over_Z',
      all(s != 0 for s in sums))
vanish = {m: sum(1 for s in sums if s % m == 0) for m in range(3, 15)}
CHECK('the_published_assignment_hits_a_multiple_of_m_for_every_m_at_most_8',
      all(vanish[m] > 0 for m in range(3, 9)) and
      all(vanish[m] == 0 for m in range(9, 15)),
      'vanishing forms per m: %s' % vanish)
sys.stdout.flush()

# ============================================================================
# 4. THE THREE EXHIBITED WITNESSES, AND A SECOND FAMILY
# ============================================================================
W1 = {
    6: dict(a=1, b=3, c=3, d=4, e=2, f=1, g=4, h=2, j=2, k=3, l=1, p=5,
            q=3, r=3, s=2, t=5, u=5, v=1, w=3, alpha=2, beta=0, gamma=5,
            delta=1, epsilon=3),
    7: dict(a=1, b=3, c=3, d=4, e=4, f=1, g=6, h=2, j=2, k=3, l=1, p=4,
            q=2, r=0, s=4, t=4, u=2, v=1, w=0, alpha=5, beta=3, gamma=0,
            delta=6, epsilon=5),
}
W1[8] = dict(W1[7])
W1[8]['g'] = 7
W2 = {
    6: dict(a=1, b=2, c=2, d=1, e=1, f=3, g=5, h=4, j=3, k=5, l=1, p=4,
            q=1, r=4, s=3, t=3, u=0, v=1, w=1, alpha=3, beta=4, gamma=3,
            delta=2, epsilon=4),
    7: dict(a=1, b=2, c=2, d=1, e=1, f=3, g=5, h=4, j=2, k=4, l=1, p=5,
            q=6, r=0, s=3, t=3, u=0, v=1, w=1, alpha=6, beta=5, gamma=1,
            delta=4, epsilon=2),
    8: dict(a=1, b=2, c=2, d=1, e=1, f=3, g=4, h=4, j=2, k=4, l=1, p=5,
            q=3, r=0, s=1, t=2, u=0, v=1, w=1, alpha=6, beta=7, gamma=1,
            delta=2, epsilon=2),
}


def vec(dd, m):
    return [dd[n] % m for n in ARCNAME]


for tag, W in (('first', W1), ('second', W2)):
    for m in (6, 7, 8):
        vv = vec(W[m], m)
        CHECK('%s_family_m%d_exhibited_residues_already_lie_in_Z_%d'
              % (tag, m, m),
              vv == [W[m][n] for n in ARCNAME] and all(0 <= x < m for x in vv))
        n3, a3, ec = build_lift(m, vv)
        dd = Counter(len(x) for x in a3)
        g1, g2 = girth_bfs(n3, a3), girth_edgedel(n3, a3)
        CHECK('%s_family_m%d_lift_is_a_3_m_graph_of_order_49m_plus_3'
              % (tag, m),
              n3 == 49 * m + 3 and ec == 75 * m and dd == {3: 49 * m, m: 3},
              'n=%d size=%d degrees %s' % (n3, ec, dict(sorted(dd.items()))))
        CHECK('%s_family_m%d_lift_has_girth_exactly_12_by_both_algorithms'
              % (tag, m), g1 == 12 and g2 == 12,
              'BFS-from-every-vertex %d, delete-edge-then-BFS %d' % (g1, g2))
        CHECK('%s_family_m%d_satisfies_the_criterion_no_short_form_vanishes'
              % (tag, m),
              all(form_sum(d, vv) % m != 0 for d in CONS),
              'all 126 forms nonzero mod %d, checked independently of the lift'
              % m)
sys.stdout.flush()

# ============================================================================
# 5. THE NO HALF AT m = 4: a four-variable pigeonhole, no computer needed
# ============================================================================
J = {n: i for i, n in enumerate(ARCNAME)}
HAND = ['k', 'p', 'v', 'alpha']
MACH = ['b', 'd', 'f', 'h']


def pigeonhole_forms(names):
    S = [J[n] for n in names]
    out, tags = [], []
    for n in names:
        f = [d for d in CONS if set(d) == {J[n]}]
        if f:
            out.append(f[0])
            tags.append('%s!=0' % n)
    for n1, n2 in itertools.combinations(names, 2):
        f = [d for d in CONS if set(d) == {J[n1], J[n2]} and
             d[J[n1]] * d[J[n2]] == -1]
        if f:
            out.append(f[0])
            tags.append('%s!=%s' % (n1, n2))
    return out, tags, S


def sat(forms, S, m):
    for asg in itertools.product(range(m), repeat=len(S)):
        mp = dict(zip(S, asg))
        if all(sum(sg * mp[k] for k, sg in d.items()) % m != 0 for d in forms):
            return list(asg)
    return None


PH, TAGS, SH = pigeonhole_forms(HAND)
CHECK('hand_core_k_p_v_alpha_yields_four_singleton_forms_from_10_cycles',
      len([t for t in TAGS if t.endswith('!=0')]) == 4,
      'each of k,p,v,alpha closes a fundamental cycle of length 10, so each '
      'is forced nonzero')
CHECK('hand_core_k_p_v_alpha_yields_all_six_difference_forms',
      len([t for t in TAGS if not t.endswith('!=0')]) == 6,
      'pair cycle lengths k-p 6, v-alpha 6, k-v 8, p-alpha 8, k-alpha 10, '
      'p-v 10')
CHECK('the_hand_core_is_exactly_ten_forms', len(PH) == 10, ' ; '.join(TAGS))
CHECK('every_form_of_both_hand_cores_really_is_one_of_the_126_cycle_forms',
      all(d in CONS for d in PH),
      'nothing is assumed: each constraint is an enumerated short cycle')
CHECK('the_hand_core_is_UNSAT_over_Z_4_all_256_assignments_fail',
      sat(PH, SH, 4) is None,
      'k,p,v,alpha must be four pairwise-distinct NONZERO residues, which '
      'needs m-1 >= 4; Z_4 has only three nonzero residues -- pigeonhole')
CHECK('the_hand_core_also_re_proves_m_3_impossible', sat(PH, SH, 3) is None)
CHECK('the_hand_core_is_VACUOUS_at_m_5_6_7_8_so_it_cannot_have_killed_them',
      all(sat(PH, SH, m) is not None for m in (5, 6, 7, 8)),
      'e.g. a satisfying (k,p,v,alpha) at m=5 is %s' % sat(PH, SH, 5))
six = []
for W in (W1, W2):
    for m in (6, 7, 8):
        vv = vec(W[m], m)
        six.append((m, tuple(vv[J[n]] for n in HAND),
                    all(form_sum(d, vv) % m != 0 for d in PH)))
CHECK('all_six_exhibited_witnesses_satisfy_the_hand_core_subsystem',
      all(ok for _, _, ok in six),
      '(k,p,v,alpha) = %s' % [(m, t) for m, t, _ in six])

PM, TAGM, SM = pigeonhole_forms(MACH)
CHECK('a_second_independent_hand_core_b_d_f_h_has_ten_forms_and_is_UNSAT_at_4',
      len(PM) == 10 and all(d in CONS for d in PM) and sat(PM, SM, 4) is None,
      'the z-arcs of the four x-leaves below x_01: %s' % ' ; '.join(TAGM))
CHECK('the_second_core_is_also_vacuous_at_m_5', sat(PM, SM, 5) is not None)
sys.stdout.flush()

# ============================================================================
# 6. EXHAUSTION OVER THE COMPLETE CELL Z_m^24, NO SYMMETRY REDUCTION AT ALL
# ============================================================================
def exhaust(m):
    order, rem, done = [], set(range(24)), [False] * len(CONS)
    while rem:
        best, bs = None, -1
        for j in sorted(rem):
            cnt = 0
            for ci, d in enumerate(CONS):
                if done[ci] or j not in d:
                    continue
                if all((k in order) or k == j for k in d):
                    cnt += 1
            if cnt > bs:
                bs, best = cnt, j
        order.append(best)
        rem.discard(best)
        for ci, d in enumerate(CONS):
            if not done[ci] and all(k in order for k in d):
                done[ci] = True
    trig = [[] for _ in range(24)]
    placed = set()
    for pos, j in enumerate(order):
        placed.add(j)
        for ci, d in enumerate(CONS):
            if j in d and all(k in placed for k in d):
                trig[pos].append(ci)
    val = [0] * 24
    nodes = [0]
    sol = [None]

    def rec(pos):
        if pos == 24:
            sol[0] = list(val)
            return
        j = order[pos]
        for x in range(m):
            nodes[0] += 1
            val[j] = x
            ok = True
            for ci in trig[pos]:
                if form_sum(CONS[ci], val) % m == 0:
                    ok = False
                    break
            if ok:
                rec(pos + 1)
                if sol[0] is not None:
                    return

    rec(0)
    return sol[0], nodes[0], sum(1 for d in done if d), len(order)


for m in (4, 5):
    s, nd, ndone, nvar = exhaust(m)
    CHECK('m%d_is_UNSAT_over_the_COMPLETE_cell_Z_%d_to_the_24_on_the_printed'
          '_arc_set' % (m, m),
          s is None and ndone == 126 and nvar == 24,
          'no Z_%d voltage assignment on the 24 arcs of G_12 AS PRINTED gives '
          'girth 12; '
          'all 24 variables branched, all 126 forms enforced, no gauge '
          'quotient and no unit-orbit cut; %d search nodes' % (m, nd))
    NOTE('m=%d exhaustion visited %d nodes (node counts depend on the variable '
         'ordering and are reported, not claimed)' % (m, nd))
for m in (6, 7, 8):
    s, nd, ndone, nvar = exhaust(m)
    ok = s is not None
    if ok:
        vv = [x % m for x in s]
        nS, aS, _ = build_lift(m, vv)
        ok = girth_bfs(nS, aS) == 12
    CHECK('m%d_is_SAT_and_the_searchs_own_vector_lifts_to_girth_12' % m, ok,
          'the very search that returns UNSAT at m=4,5 returns a solution at '
          'm=%d in %d nodes, and its full lift has girth 12' % (m, nd))
sys.stdout.flush()

# ============================================================================
# 7. HOW SMALL A PAPER CERTIFICATE CAN BE: every subset of at most 4 variables
# ============================================================================
SUPPMASK = []
for d in CONS:
    mm = 0
    for k in d:
        mm |= 1 << k
    SUPPMASK.append(mm)


def cores_up_to_4(m):
    """Every UNSAT subsystem on <= 4 variables in which every variable is
    relevant.  Sound for deciding EXISTENCE: a minimal UNSAT subset has every
    variable covered by its own constraints, else dropping it stays UNSAT."""
    out = []
    for sz in (1, 2, 3, 4):
        for S in itertools.combinations(range(24), sz):
            Smask = 0
            for j in S:
                Smask |= 1 << j
            sub, cover = [], 0
            for i, mm in enumerate(SUPPMASK):
                if mm & ~Smask == 0:
                    sub.append(CONS[i])
                    cover |= mm
            if not sub or cover != Smask:
                continue
            if sat(sub, list(S), m) is None:
                out.append(tuple(ARCNAME[j] for j in S))
    return out


c4 = cores_up_to_4(4)
CHECK('at_m_4_the_only_UNSAT_cores_on_at_most_4_variables_are_three_4_cliques',
      sorted(c4) == sorted([('b', 'd', 'f', 'h'), ('j', 'k', 'l', 'p'),
                            ('k', 'p', 'v', 'alpha')]),
      '%d cores, all of size 4: %s' % (len(c4), sorted(c4)))
CHECK('at_m_4_there_is_no_UNSAT_core_on_three_or_fewer_variables',
      all(len(c) == 4 for c in c4))
CHECK('at_m_5_there_is_NO_UNSAT_core_on_four_or_fewer_variables',
      cores_up_to_4(5) == [],
      'an all-distinct clique of size 4 kills every m <= 4 and no m = 5; this '
      'search covers subsets of at most FOUR of the 24 variables and says '
      'nothing about larger subsets')

NOTE('SCOPE of this program. It re-derives the following quantities '
     'transcribed from the paper: the graph and its checksums by two '
     'constructions, the 126-cycle census by two algorithms, the five '
     'published invariants of the source, the six exhibited voltage vectors '
     'and their lift girths by two algorithms, the m=4 hand pigeonhole '
     'together with its vacuity at every m >= 5, the complete exhaustions '
     'over Z_4^24 and Z_5^24 with no symmetry reduction, and the absence of '
     'a small paper certificate at m=5 up to 4 variables. This program does '
     'not parse the paper and makes no claim to have enumerated everything '
     'the paper states. NOT RE-RUN HERE: the WIDE variant in which the 24 '
     'arcs themselves are re-paired, which is untouched at m=4,5; the '
     'standard gauge/normalisation argument, of which only the counts are '
     'checked; and any statement about published order records, of which '
     'this paper makes none.')

print()
if _F[0]:
    print('VERDICT: %d CHECK(S) FAILED of %d' % (_F[0], _N[0] + _F[0]))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _N[0])
sys.exit(0)

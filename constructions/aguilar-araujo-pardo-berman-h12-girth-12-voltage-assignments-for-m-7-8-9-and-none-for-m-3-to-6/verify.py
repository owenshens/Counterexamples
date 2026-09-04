#!/usr/bin/env python3
"""verify.py -- verification program for the note

    "Voltage assignments on H_12 giving girth 12 for m = 7, 8, 9, and none for 3 <= m <= 6"

Python 3.9 or later, standard library only: no third-party package and no external data file.

Everything this program consumes is PRINTED IN THE PAPER: the 44 vertices of the base graph
H_12, its 43 tree edges, its 20 labelled arcs in the fixed coordinate order, the deletion list
that turns T_12 into H_12, the six voltage 20-tuples, the three shortest cycles exhibited for
the three headline witnesses, and the core variable sets of the four emptiness claims.  Nothing
is read from a file and nothing is imported beyond the standard library.

All arithmetic is exact: Python integers, plus fractions.Fraction where the source's own
Corollary 2.4 has thirds in it.  No floating point enters any decision.

It prints one `PASS <name> [detail]` line per check and closes with

    VERDICT: ALL <n> CHECKS PASS

exiting 0 if and only if every check passed.
"""

import sys
from collections import deque
from fractions import Fraction
from itertools import product

_PASS = 0
_FAIL = 0


def ck(name, cond, detail=''):
    """One check.  Prints PASS or FAIL and nothing else."""
    global _PASS, _FAIL
    if cond:
        _PASS += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAIL += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    """An informational line.  Deliberately NOT a check: it asserts nothing."""
    print('    -- %s' % text)


print('verification of the note: voltage assignments on H_12 giving girth 12 for m = 7, 8, 9,')
print('and none for 3 <= m <= 6 -- everything below is derived from the objects printed in the')
print('paper; exact integer and Fraction arithmetic only, no floating point in any decision.')
print('python %s' % sys.version.split()[0])


# =====================================================================================
# 1.  THE DATA, EXACTLY AS PRINTED IN THE PAPER
# =====================================================================================

# Section 2 of the paper: the 44 vertices of H_12.  x*, y*, z* are the three pinned vertices.
VERTS_PRINTED = """
x* x x0 x1 x00 x01 x10 x11 x010 x011 x100 x101 x110 x111
x0100 x0101 x0110 x0111 x1010 x1011 x1110 x1111
y* y y0 y1 y01 y10 y11 y010 y011 y101 y110
z* z z0 z1 z01 z10 z11 z010 z011 z101 z110
""".split()

PINNED = ('x*', 'y*', 'z*')

# The 43 tree edges, all at voltage 0.
TREE = [
    ('x*', 'x'), ('x', 'x0'), ('x', 'x1'), ('x0', 'x00'), ('x0', 'x01'),
    ('x1', 'x10'), ('x1', 'x11'), ('x01', 'x010'), ('x01', 'x011'),
    ('x10', 'x100'), ('x10', 'x101'), ('x11', 'x110'), ('x11', 'x111'),
    ('x010', 'x0100'), ('x010', 'x0101'), ('x011', 'x0110'), ('x011', 'x0111'),
    ('x101', 'x1010'), ('x101', 'x1011'), ('x111', 'x1110'), ('x111', 'x1111'),
    ('y*', 'y'), ('y', 'y0'), ('y', 'y1'), ('y0', 'y01'), ('y1', 'y10'), ('y1', 'y11'),
    ('y01', 'y010'), ('y01', 'y011'), ('y10', 'y101'), ('y11', 'y110'),
    ('z*', 'z'), ('z', 'z0'), ('z', 'z1'), ('z0', 'z01'), ('z1', 'z10'), ('z1', 'z11'),
    ('z01', 'z010'), ('z01', 'z011'), ('z10', 'z101'), ('z11', 'z110'),
    ('x00', 'y0'), ('x00', 'z0'),
]

# The 20 labelled arcs, in THE fixed coordinate order of every voltage tuple below:
#   [a, b, c, d, e, j, q, u, w, s, l, g, f, k, r, v, al, t, p, h]
ARCS = [
    ('a', 'x100', 'y10'), ('b', 'x100', 'z10'), ('c', 'x110', 'y11'), ('d', 'x110', 'z11'),
    ('e', 'x0100', 'y010'), ('j', 'x0101', 'y011'), ('q', 'x0110', 'y101'),
    ('u', 'x0111', 'y110'), ('w', 'x1010', 'y110'), ('s', 'x1011', 'y101'),
    ('l', 'x1110', 'y011'), ('g', 'x1111', 'y010'),
    ('f', 'x0100', 'z101'), ('k', 'x0101', 'z110'), ('r', 'x0110', 'z010'),
    ('v', 'x0111', 'z011'), ('al', 'x1010', 'z010'), ('t', 'x1011', 'z011'),
    ('p', 'x1110', 'z101'), ('h', 'x1111', 'z110'),
]
ARC_NAMES = [a[0] for a in ARCS]
IDX = {n: i for i, n in enumerate(ARC_NAMES)}

# The eight leaves of T_12 deleted on p. 27 of the source, with their parents in T_12.
DELETED = [('x100', 'x1000'), ('x100', 'x1001'), ('x110', 'x1100'), ('x110', 'x1101'),
           ('y10', 'y100'), ('y11', 'y111'), ('z10', 'z100'), ('z11', 'z111')]

# The six voltage tuples of the paper.
W = {7: [1, 2, 2, 1, 2, 4, 6, 3, 0, 5, 5, 6, 1, 5, 3, 1, 1, 4, 0, 3],
     8: [1, 2, 2, 1, 3, 7, 3, 1, 0, 4, 1, 7, 1, 3, 7, 3, 1, 7, 0, 4],
     9: [1, 2, 2, 1, 7, 6, 3, 1, 0, 5, 1, 3, 1, 3, 6, 7, 1, 3, 0, 5]}
# The restricted witnesses: a, b, c, d pinned to the source's published Table-4 values 1, -1, -1, 1.
R = {7: [1, 6, 6, 1, 6, 5, 2, 5, 0, 5, 3, 2, 1, 6, 4, 6, 5, 4, 2, 5],
     8: [1, 7, 7, 1, 6, 5, 5, 1, 0, 3, 6, 1, 1, 4, 7, 5, 1, 2, 0, 2],
     9: [1, 8, 8, 1, 1, 6, 3, 1, 0, 7, 1, 6, 2, 4, 5, 8, 1, 2, 0, 3]}

# The three 12-cycles exhibited in the paper, as lift vertices (base vertex, fibre index).
SHORT_CYCLE = {
    7: [('y11', 5), ('y1', 5), ('y', 5), ('y0', 5), ('y01', 5), ('y010', 5),
        ('x0100', 3), ('z101', 4), ('z10', 4), ('z1', 4), ('z11', 4), ('x110', 3)],
    8: [('x100', 7), ('y10', 0), ('y1', 0), ('y', 0), ('y0', 0), ('y01', 0),
        ('y011', 0), ('x1110', 7), ('x111', 7), ('x11', 7), ('x1', 7), ('x10', 7)],
    9: [('x0', 7), ('x00', 7), ('z0', 7), ('z', 7), ('z1', 7), ('z10', 7),
        ('z101', 7), ('x1110', 7), ('x111', 7), ('x11', 7), ('x1', 7), ('x', 7)],
}

# The core variable sets of the four emptiness claims (paper, Table of cores).
CORES = {3: ['g', 'f', 'k', 'p', 'h'],
         4: ['u', 'w', 's', 'r', 'v', 'al', 't'],
         5: ['q', 'u', 'w', 's', 'r', 'v', 'al', 't'],
         6: ['q', 'u', 'w', 's', 'r', 'v', 'al', 't']}

# The ten relations the paper prints for m = 3, as coefficient vectors over (g, f, k, p, h).
M3_RELATIONS = [
    ('g_ne_0', (1, 0, 0, 0, 0)),
    ('f_ne_0', (0, 1, 0, 0, 0)),
    ('k_ne_0', (0, 0, 1, 0, 0)),
    ('h_ne_g', (-1, 0, 0, 0, 1)),
    ('h_ne_k', (0, 0, -1, 0, 1)),
    ('h_ne_g_plus_k', (-1, 0, -1, 0, 1)),
    ('f_ne_k', (0, 1, -1, 0, 0)),
    ('p_ne_f', (0, -1, 0, 1, 0)),
    ('p_ne_h', (0, 0, 0, 1, -1)),
    ('f_minus_k_minus_p_plus_h_ne_0', (0, 1, -1, -1, 1)),
]
M3_VARS = ('g', 'f', 'k', 'p', 'h')


# =====================================================================================
# 2.  THE BASE GRAPH
# =====================================================================================

def build_base():
    """-> (vertices, adjacency, edge_label) for H_12.

    `edge_label[(u, v)]` is None for a tree edge and (arc_index, +1) when the edge is the arc
    traversed from tail to head (so (arc_index, -1) for the reverse orientation).
    """
    verts, adj, lab = [], {}, {}

    def add(u, v, label):
        for x in (u, v):
            if x not in adj:
                adj[x] = []
                verts.append(x)
        adj[u].append(v)
        adj[v].append(u)
        lab[(u, v)] = label
        lab[(v, u)] = None if label is None else (label[0], -label[1])

    for u, v in TREE:
        add(u, v, None)
    for i, (_n, tail, head) in enumerate(ARCS):
        add(tail, head, (i, 1))
    return verts, adj, lab


VERTS, ADJ, LAB = build_base()

ck('base_vertex_set_is_the_44_printed_vertices',
   sorted(VERTS) == sorted(VERTS_PRINTED) and len(VERTS) == 44,
   '%d vertices, and the tree edges span exactly the printed list' % len(VERTS))
ck('base_tree_edge_count_43', len(TREE) == 43, '43 tree edges, all at voltage 0')
ck('base_arc_count_20', len(ARCS) == 20 and len(set(ARC_NAMES)) == 20,
   'arc order [%s]' % ','.join(ARC_NAMES))
ck('base_edge_count_63', len(TREE) + len(ARCS) == 63, '43 + 20 = 63')

_undirected = set(frozenset(e) for e in TREE) | set(frozenset((t, h)) for _n, t, h in ARCS)
ck('base_has_no_loop_and_no_repeated_edge',
   len(_undirected) == 63 and all(len(e) == 2 for e in _undirected),
   '63 distinct unordered pairs')

_deg = {v: len(ADJ[v]) for v in VERTS}
_hist = {}
for v in VERTS:
    _hist[_deg[v]] = _hist.get(_deg[v], 0) + 1
ck('base_degree_histogram_is_1x3_and_3x41', _hist == {1: 3, 3: 41}, 'got %s' % _hist)
ck('base_degree_one_vertices_are_exactly_the_three_pinned',
   sorted(v for v in VERTS if _deg[v] == 1) == sorted(PINNED), 'x*, y*, z*')
ck('base_non_pinned_count_41_reproduces_theorem_5_6',
   len(VERTS) - 3 == 41, '41 vertices of degree 3, i.e. order 41m + 3')


def bfs_component(adj, start):
    seen, q = {start}, deque([start])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                q.append(w)
    return seen


_tree_adj = {v: [] for v in VERTS}
for u, v in TREE:
    _tree_adj[u].append(v)
    _tree_adj[v].append(u)
ck('the_43_tree_edges_form_a_spanning_tree',
   len(bfs_component(_tree_adj, 'x*')) == 44 and len(TREE) == 43,
   'BFS from x* reaches all 44 vertices on 44-1 = 43 edges')
ck('cycle_rank_of_the_base_is_20', 63 - 44 + 1 == 20,
   '|E| - |V| + 1 = 63 - 44 + 1 = 20, so the 20 arcs are a cotree')


def two_colour(adj, start):
    col = {start: 0}
    q = deque([start])
    bad = 0
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in col:
                col[w] = 1 - col[u]
                q.append(w)
            elif col[w] == col[u]:
                bad += 1
    return col, bad


_col, _bad = two_colour(ADJ, 'x*')
_sizes = sorted([sum(1 for v in VERTS if _col[v] == 0), sum(1 for v in VERTS if _col[v] == 1)])
ck('base_is_bipartite_with_classes_of_size_21_and_23',
   len(_col) == 44 and _bad == 0 and _sizes == [21, 23],
   'zero monochromatic edges, class sizes %s' % _sizes)
ck('every_one_of_the_20_arcs_joins_opposite_classes',
   all(_col[t] != _col[h] for _n, t, h in ARCS),
   'so every lift (H_12, m) is bipartite as well')

# T_12, recovered by putting the eight deleted leaves back.
_t12_v = set(VERTS) | set(c for _p, c in DELETED)
ck('the_eight_deleted_leaves_are_absent_from_H12',
   all(c not in VERTS for _p, c in DELETED) and all(_p in VERTS for _p, _c in DELETED),
   ', '.join(c for _p, c in DELETED))
ck('T12_has_52_vertices_and_51_tree_edges',
   len(_t12_v) == 52 and len(TREE) + len(DELETED) == 51,
   '44 + 8 = 52 and 43 + 8 = 51')
ck('T12_non_pinned_count_49_reproduces_the_sources_49m_plus_3',
   52 - 3 == 49, "the source's own 49m + 3 and its 441 = 49 * 9")


# =====================================================================================
# 3.  THE SHORT-CYCLE CENSUS AND THE CONSTRAINT SYSTEM
# =====================================================================================

def cycles_upto(limit):
    """Every simple cycle of length <= limit, each once, as a vertex list."""
    order = {v: i for i, v in enumerate(sorted(VERTS))}
    out = []

    def dfs(start, path, on):
        u = path[-1]
        for w in ADJ[u]:
            if w == start:
                if len(path) >= 3 and order[path[1]] < order[path[-1]]:
                    out.append(list(path))
            elif order[w] > order[start] and w not in on and len(path) < limit:
                on.add(w)
                path.append(w)
                dfs(start, path, on)
                path.pop()
                on.discard(w)

    for s in sorted(VERTS, key=lambda v: order[v]):
        dfs(s, [s], {s})
    return out


def constraint_of(cycle):
    """The 20-vector of arc coefficients picked up by traversing `cycle` once."""
    co = [0] * 20
    n = len(cycle)
    for i in range(n):
        e = LAB[(cycle[i], cycle[(i + 1) % n])]
        if e is not None:
            co[e[0]] += e[1]
    return tuple(co)


CYC = cycles_upto(11)
_len_hist = {}
for c in CYC:
    _len_hist[len(c)] = _len_hist.get(len(c), 0) + 1

ck('shortest_cycle_of_the_base_graph_has_length_6', min(_len_hist) == 6,
   'so the list of lengths below 12 is complete at 6, 8, 10')
ck('base_has_no_cycle_of_odd_length_below_12',
   all(L % 2 == 0 for L in _len_hist), 'lengths present: %s' % sorted(_len_hist))
ck('base_6_cycles_number_18', _len_hist.get(6) == 18, 'got %s' % _len_hist.get(6))
ck('base_8_cycles_number_27', _len_hist.get(8) == 27, 'got %s' % _len_hist.get(8))
ck('base_10_cycles_number_82', _len_hist.get(10) == 82, 'got %s' % _len_hist.get(10))
ck('base_cycles_shorter_than_12_number_127', len(CYC) == 127,
   '18 + 27 + 82 = %d undirected' % len(CYC))
ck('directed_short_cycles_number_254_the_sources_published_integer',
   2 * len(CYC) == 254, "reproduces 'There are 254 such cycles', source p. 27")

CONS = [constraint_of(c) for c in CYC]
ck('every_constraint_coefficient_is_plus_or_minus_one',
   all(all(x in (-1, 0, 1) for x in co) for co in CONS),
   'no arc is traversed twice by a simple cycle of length below 12')

_ar = {}
for co in CONS:
    a = sum(1 for x in co if x)
    _ar[a] = _ar.get(a, 0) + 1
ck('constraint_arity_histogram_is_16_44_4_55_6_2',
   _ar == {1: 16, 2: 44, 3: 4, 4: 55, 6: 6, 8: 2}, 'got %s' % _ar)

_forced = sorted(ARC_NAMES[i] for co in CONS if sum(1 for x in co if x) == 1
                 for i in range(20) if co[i])
ck('exactly_16_arcs_are_forced_nonzero_by_a_single_arc_short_cycle',
   len(set(_forced)) == 16, '%d distinct arcs' % len(set(_forced)))
ck('the_four_unforced_arcs_are_w_s_p_h',
   sorted(set(ARC_NAMES) - set(_forced)) == sorted(['w', 's', 'p', 'h']),
   'got %s' % sorted(set(ARC_NAMES) - set(_forced)))


def satisfies_all(phi, m):
    """True iff every base cycle shorter than 12 has voltage sum != 0 mod m."""
    for co in CONS:
        if sum(c * phi[i] for i, c in enumerate(co)) % m == 0:
            return False
    return True


# =====================================================================================
# 4.  THE LIFTS
# =====================================================================================

def build_lift(phi, m):
    """The derived graph (H_12, phi) over Z_m: adjacency dict, plus loop and multi-edge counts."""
    adj = {}

    def V(v, i):
        return v if v in PINNED else (v, i)

    for v in VERTS:
        if v in PINNED:
            adj[v] = []
        else:
            for i in range(m):
                adj[(v, i)] = []

    pairs = []
    for u, v in TREE:
        for i in range(m):
            pairs.append((V(u, i), V(v, i)))
    for idx, (_n, t, h) in enumerate(ARCS):
        c = phi[idx] % m
        for i in range(m):
            pairs.append(((t, i), (h, (i + c) % m)))

    loops = sum(1 for a, b in pairs if a == b)
    seen = set()
    multi = 0
    for a, b in pairs:
        key = frozenset((a, b)) if a != b else (a,)
        if key in seen:
            multi += 1
        seen.add(key)
        adj[a].append(b)
        adj[b].append(a)
    return adj, len(pairs), loops, multi


def girth_at_most_12(adj):
    """The exact girth when it is <= 12, else None.  BFS from every vertex, depth capped at 6.

    For a graph of girth g <= 12 there is a root on a shortest cycle from which the two halves
    meet at depth floor(g/2) <= 6, so the cap loses nothing.
    """
    best = None
    for root in adj:
        dist = {root: 0}
        par = {root: None}
        q = deque([root])
        while q:
            u = q.popleft()
            du = dist[u]
            for w in adj[u]:
                if w not in dist:
                    if du < 6:
                        dist[w] = du + 1
                        par[w] = u
                        q.append(w)
                elif w != par[u]:
                    cand = du + dist[w] + 1
                    if best is None or cand < best:
                        best = cand
    return best


for tag, table in (('W', W), ('R', R)):
    for m in (7, 8, 9):
        phi = table[m]
        adj, ne, loops, multi = build_lift(phi, m)
        pre = '%s%d' % (tag, m)
        ck('%s_lift_order_is_41m_plus_3' % pre, len(adj) == 41 * m + 3,
           'order %d = 41*%d + 3' % (len(adj), m))
        ck('%s_lift_edge_count_is_63m' % pre, ne == 63 * m, '%d = 63*%d' % (ne, m))
        ck('%s_lift_has_no_loop_and_no_multi_edge' % pre, loops == 0 and multi == 0,
           'loops %d, repeated pairs %d' % (loops, multi))
        h = {}
        for v in adj:
            h[len(adj[v])] = h.get(len(adj[v]), 0) + 1
        ck('%s_lift_degree_histogram_is_3x41m_and_mx3' % pre, h == {3: 41 * m, m: 3},
           'got %s' % h)
        ck('%s_lift_is_connected' % pre, len(bfs_component(adj, next(iter(adj)))) == len(adj),
           'one component of %d' % len(adj))
        g = girth_at_most_12(adj)
        ck('%s_lift_girth_is_exactly_12' % pre, g == 12, 'girth = %s' % g)
        ck('%s_satisfies_all_127_necessary_short_cycle_constraints' % pre,
           satisfies_all(phi, m), 'every base cycle of length 6, 8, 10 has nonzero voltage sum')

# The three exhibited 12-cycles, checked edge by edge -- the reader's hand-sized certificate.
for m in (7, 8, 9):
    adj, _ne, _lo, _mu = build_lift(W[m], m)
    cyc = [(v if v in PINNED else (v, i)) for v, i in SHORT_CYCLE[m]]
    ok = (len(cyc) == 12 and len(set(cyc)) == 12
          and all(c in adj for c in cyc)
          and all(cyc[(i + 1) % 12] in adj[cyc[i]] for i in range(12)))
    ck('W%d_printed_shortest_cycle_is_a_genuine_12_cycle_of_the_lift' % m, ok,
       '12 distinct lift vertices, all 12 consecutive pairs are edges')

# The restricted witnesses keep the source's own a, b, c, d = 1, -1, -1, 1.
for m in (7, 8, 9):
    ck('R%d_keeps_the_published_table_4_values_of_a_b_c_d' % m,
       [R[m][IDX[n]] % m for n in ('a', 'b', 'c', 'd')] == [1 % m, -1 % m, -1 % m, 1 % m],
       'a, b, c, d = 1, %d, %d, 1 = 1, -1, -1, 1 mod %d' % (-1 % m, -1 % m, m))
    ck('R%d_differs_from_W%d_so_the_restricted_cell_is_separately_non_empty' % (m, m),
       R[m] != W[m], 'the sub-bullet prescribes keeping a, b, c, d')


# =====================================================================================
# 5.  THE EMPTINESS OF m = 3, 4, 5, 6
# =====================================================================================
# For each m the paper names a set of core arcs.  Every constraint whose support lies inside
# that set is collected, and EVERY tuple of Z_m over the core is enumerated flat -- no search,
# no propagation, no heuristic.  Zero survivors decides the whole cell m^20, because the
# constraints are NECESSARY for girth >= 12 and a partial assignment that already violates one
# cannot be completed.

_total_cells = 0
for m in (3, 4, 5, 6):
    core = CORES[m]
    ids = [IDX[n] for n in core]
    inside = []
    for co in CONS:
        if all(co[i] == 0 for i in range(20) if i not in ids):
            inside.append(tuple(co[i] for i in ids))
    inside.sort(key=lambda v: sum(1 for x in v if x))
    survivors = 0
    seen = 0
    for tup in product(range(m), repeat=len(core)):
        seen += 1
        if all(sum(c * tup[j] for j, c in enumerate(co)) % m for co in inside):
            survivors += 1
    ck('m%d_core_enumeration_is_the_full_flat_cell' % m, seen == m ** len(core),
       'core (%s): %d = %d^%d tuples enumerated, no search' % (','.join(core), seen, m, len(core)))
    ck('m%d_core_constraints_are_short_cycle_constraints_supported_in_the_core' % m,
       len(inside) > 0 and all(sum(1 for x in c if x) == sum(1 for x in co if x)
                               for c, co in zip(inside, sorted(
                                   (co for co in CONS
                                    if all(co[i] == 0 for i in range(20) if i not in ids)),
                                   key=lambda v: sum(1 for x in v if x)))),
       '%d of the 127 constraints have their whole support inside the core' % len(inside))
    ck('m%d_core_has_zero_survivors' % m, survivors == 0,
       'so no assignment in Z_%d^20 gives girth >= 12' % m)
    ck('m%d_lift_order_would_have_been_41m_plus_3' % m, 41 * m + 3 == (126, 167, 208, 249)[m - 3],
       'order %d' % (41 * m + 3))
    _total_cells += m ** 20
    note('m=%d: cell size %d^20 = %d, decided in full' % (m, m, m ** 20))

ck('the_four_empty_cells_total_3752628870115778_assignments',
   _total_cells == 3752628870115778, '%d' % _total_cells)

# ---- m = 3: the ten printed relations, and the four-case pencil proof -----------------
_ids3 = [IDX[n] for n in M3_VARS]
_derived3 = set()
for co in CONS:
    if all(co[i] == 0 for i in range(20) if i not in _ids3):
        v = tuple(co[i] for i in _ids3)
        first = next(x for x in v if x)
        _derived3.add(v if first > 0 else tuple(-x for x in v))


def canon(v):
    first = next((x for x in v if x), 0)
    return tuple(-x for x in v) if first < 0 else tuple(v)


_printed3 = set(canon(v) for _n, v in M3_RELATIONS)
ck('the_ten_relations_printed_for_m3_are_all_genuine_constraints_of_H12',
   _printed3 <= _derived3,
   '%d printed relations, all among the %d constraints supported on (g,f,k,p,h)'
   % (len(_printed3), len(_derived3)))
note('constraints of H_12 supported on (g,f,k,p,h): %d' % len(_derived3))

_surv3 = 0
for tup in product(range(3), repeat=5):
    if all(sum(c * tup[j] for j, c in enumerate(v)) % 3 for _n, v in M3_RELATIONS):
        _surv3 += 1
ck('the_ten_printed_relations_alone_have_no_solution_over_Z3', _surv3 == 0,
   'flat enumeration of all 3^5 = 243 tuples, 0 survivors')

# The four-case argument of the paper, re-derived case by case.
_cases = 0
for g in (1, 2):
    for k in (1, 2):
        f = 3 - k                                   # f != 0 and f != k force f = the other of {1,2}
        hs = [x for x in range(3) if x != g % 3 and x != k % 3 and x != (g + k) % 3]
        if not hs:
            _cases += 1                             # no value left for h
            continue
        h = hs[0]
        ps = [x for x in range(3) if x != f % 3 and x != h % 3]
        if len(ps) == 1 and (f - k - ps[0] + h) % 3 == 0:
            _cases += 1                             # the arity-4 relation is violated
ck('the_four_case_pencil_proof_for_m3_closes_all_four_branches', _cases == 4,
   'each of (g,k) in {1,2}^2 ends in a contradiction, so no Z_3 assignment exists')
ck('the_m3_lift_would_have_been_cubic_on_126_vertices',
   41 * 3 + 3 == 126, 'at m = 3 the three pinned vertices also have degree 3')


# =====================================================================================
# 6.  THE ARITHMETIC OF THE SOURCE AND OF THE RECORD TABLE
# =====================================================================================

def lb(m):
    """The published lower bound of Lemma 3.4 of [ABGMV]: ceil(109m/3) + 17."""
    return -((-109 * m) // 3) + 17


ck('lower_bound_formula_reproduces_the_sources_own_four_values',
   (lb(3), lb(7), lb(9), lb(10)) == (126, 272, 344, 381),
   'ceil(109m/3)+17 = 126, 272, 344, 381 at m = 3, 7, 9, 10')
ck('lower_bound_at_m8_is_308', lb(8) == 308, 'ceil(872/3) + 17 = 291 + 17')
ck('order_of_the_m7_witness_is_290_equals_41x7_plus_3', 41 * 7 + 3 == 290)
ck('order_of_the_m8_witness_is_331_equals_41x8_plus_3', 41 * 8 + 3 == 331)
ck('order_of_the_m9_witness_is_372_equals_41x9_plus_3', 41 * 9 + 3 == 372,
   'the order the source names on its p. 28')
ck('the_sources_own_413_and_444_reproduce',
   41 * 10 + 3 == 413 and 49 * 9 + 3 == 444, '(H_12, 10) and (G_12, 9)')
ck('290_lies_strictly_inside_the_published_window', lb(7) < 290 < 334,
   '272 < 290 < 334, remaining gap %d' % (290 - lb(7)))
ck('331_lies_strictly_inside_the_published_window', lb(8) < 331 < 374,
   '308 < 331 < 374, remaining gap %d' % (331 - lb(8)))
ck('372_does_not_improve_the_published_360_at_m9', 372 > 360,
   'so the m = 9 witness carries NO record claim')


def cor24(m):
    """Corollary 2.4 of the source, verbatim, in exact rational arithmetic."""
    if m % 3 == 0:
        return 41 * m + Fraction(m, 3) + 2
    if m % 3 == 1:
        return 41 * m + Fraction(m, 3) + 86 + Fraction(2, 3)
    return 41 * m + Fraction(m, 3) + 43 + Fraction(1, 3)


ck('corollary_2_4_reproduces_five_of_the_six_table_1_entries',
   (cor24(4), cor24(5), cor24(6), cor24(7), cor24(9)) == (252, 250, 250, 376, 374),
   'm = 4, 5, 6, 7, 9 give 252, 250, 250, 376, 374 exactly')
ck('corollary_2_4_at_m8_is_374_not_the_304_printed_in_table_1', cor24(8) == 374,
   '328 + 8/3 + 43 + 1/3 = 374, and Table 1 prints 304')
ck('the_printed_304_is_below_the_lower_bound_of_its_own_row', 304 < lb(8),
   '304 < 308, so the printed 304 is inconsistent with the lower bound of its own row; '
   "the corollary's own value at m=8 is 374 (checked above), but no erratum is verified here")
ck('table_1_lower_bounds_at_m5_and_m6_also_disagree_with_the_papers_own_formula',
   (lb(5), lb(6)) == (199, 235), 'the formula gives 199 and 235; Table 1 prints 202 and 237')

# GJV Theorem 4.2 with r = 3, g = 12, n_12 = 126, and s = 9 pinned by their own two published
# improvements.  Cases 2 and 3 additionally require s even, so there the usable s is at most 8.
NG, S = 126, 9


def th42_case1(k):
    return k * (NG - S) + S


def th42_case2(k, s):
    return k * (NG - s) + NG + s


def th42_case3(k, l, s):
    return k * (NG - s) + NG + (l - 1) * (NG - s)


ck('gjv_theorem_4_2_case_1_reproduces_their_published_360_at_m9', th42_case1(3) == 360,
   'k = 3, s = 9: 3*117 + 9 = 360, which pins s = 9')
ck('gjv_theorem_4_2_case_1_reproduces_their_published_243_at_m6', th42_case1(2) == 243,
   'k = 2, s = 9: 2*117 + 9 = 243, the same s again')
ck('gjv_theorem_4_2_selected_case_values_have_minimum_362_at_m7',
   min(th42_case2(2, 8), th42_case3(1, 2, 8)) == 362,
   'for case 2 at (k,s)=(2,8), value %d; for case 3 at (k,l,s)=(1,2,8), value %d; NOT '
   'CHECKED here: any sweep over all admissible cases and parameters'
   % (th42_case2(2, 8), th42_case3(1, 2, 8)))
ck('gjv_theorem_4_2_case_3_with_k2_l1_s8_gives_362_at_m8', th42_case3(2, 1, 8) == 362,
   'th42_case3(2, 1, 8) = 378 - 2*8 = 362; NOT CHECKED here: the applicability conditions '
   'of cases 1 and 2 at m = 8, nor a sweep over all admissible (k, l, s), so this is one '
   'value of case 3, not a proven floor')
ck('290_is_outside_the_reach_of_theorem_4_2_by_72', 362 - 290 == 72)
ck('331_is_outside_the_reach_of_theorem_4_2_by_31', 362 - 331 == 31)

print('NOTE SCOPE: this program checks the objects and the arithmetic printed in the paper. '
      'NOT RE-RUN here: the uniqueness of the (3,12)-cage, which the paper quotes from the '
      'literature and flags as an inherited external fact; the existence of the s/2 = 4 edges '
      'at pairwise distance >= 6 that case 3 of GJV Theorem 4.2 needs, nor the applicability '
      'conditions of its cases 1 and 2, nor any minimisation over its parameters -- the two '
      'case-value checks above are single evaluations, not floors; the numerical values 334, '
      '374 and 360 '
      'of the published record table, which are transcribed from the literature and not '
      'recomputed; and the total number of girth-12 assignments at any m, which the paper does '
      'not claim. The four emptiness claims are re-proved here in full, by flat enumeration '
      'over the printed cores with no search.')

if _FAIL:
    print('VERIFICATION FAILED: %d of %d checks did not pass' % (_FAIL, _PASS + _FAIL))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _PASS)
sys.exit(0)

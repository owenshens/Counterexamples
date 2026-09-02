#!/usr/bin/env python3
"""Verification program for

    "A 15-Vertex Counterexample to the Minimum-Degree-Two Question of
     Mynhardt and Neilson"

Python 3.9+, STANDARD LIBRARY ONLY (no numpy/sympy/networkx), no external data
file: every input is the object printed in the paper, transcribed below.  All
arithmetic is exact integer / bitmask arithmetic; no floating point is used and
no decision depends on one.

What it does, in order:

  A. rebuilds G' from the 26 named edges printed in the paper, and checks that
     the paper's graph6 string decodes to exactly that LABELLED graph;
  B. checks the doubling structure: deleting the four primed vertices returns
     the paper's G_1, the four doubled vertices are end-vertices of G_1, each
     primed vertex is a true twin, and the doubling changes no distance and no
     eccentricity;
  C. re-derives, for the witness broadcast f of the paper, every set the paper
     prints (N_f, B_f, PB_f), its weight, its legality, its irredundance, its
     bn-independence, and the fact that it is NOT dominating;
  D. runs its own exhaustive census of the broadcast space and re-derives
     Gamma_b(G') = 5, alpha_bnr(G') = 6, alpha_bn(G') = 7 together with the
     three counts 283 / 1140 / 1416, and the same four numbers for G_1;
  E. validates the census enumerator itself against an UNPRUNED brute force
     over all 3,686,400 broadcasts of G_1, and against published values for
     P_6, P_7, C_6 and two grids;
  F. re-runs the end-vertex doubling sweep: every connected labelled graph on
     3..6 vertices, every end-vertex, 22,654 pairs, 0 mismatches.

One `PASS <name> [detail]` line per check; `NOTE` lines carry data that is not
itself a check.  Closes with `VERDICT: ALL <n> CHECKS PASS` and exits 0 iff
every check passed.
"""

import functools
import itertools
import sys
import time

# ----------------------------------------------------------------------------
# check harness
# ----------------------------------------------------------------------------
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


def note(s):
    print('NOTE %s' % s)


# ----------------------------------------------------------------------------
# THE OBJECT, EXACTLY AS PRINTED IN THE PAPER
# ----------------------------------------------------------------------------
# Section 2 of the paper ("The graph"): the vertex order (index 0..14) ...
VERTS = ['u1', 'u2', 'w1', 'w2', 'x1', 'x2', 'y1', 'y2', 'z1', 'z2',
         'v', 'w1p', 'w2p', 'x1p', 'x2p']
# ... and the 26 edges, in the paper's own grouping.
EDGES = [('u1', 'u2'), ('u1', 'v'), ('u2', 'v'),
         ('z1', 'z2'), ('u1', 'z1'), ('u1', 'z2'), ('u2', 'z1'), ('u2', 'z2'),
         ('y1', 'z1'), ('y1', 'z2'), ('y2', 'z1'), ('y2', 'z2'),
         ('u1', 'y1'), ('u2', 'y2'),
         ('x1', 'y1'), ('x2', 'y2'),
         ('x1p', 'x1'), ('x1p', 'y1'), ('x2p', 'x2'), ('x2p', 'y2'),
         ('w1', 'z1'), ('w2', 'z2'),
         ('w1p', 'w1'), ('w1p', 'z1'), ('w2p', 'w2'), ('w2p', 'z2')]
# ... the graph6 string of the same labelled graph, as printed in the paper.
G6 = 'N_?CQJbs}?CCCCA_?g?'
# ... the primed (added) vertices, each a true twin of the vertex it doubles.
TWINS = [('x1p', 'x1'), ('x2p', 'x2'), ('w1p', 'w1'), ('w2p', 'w2')]
# ... the witness broadcast of Section 3 ("The counterexample").
WITNESS = {'x1': 2, 'x2': 2, 'w1': 1, 'w2': 1}
# ... and the sets the paper prints for it (names, sorted).
WITNESS_SETS = {
    'x1': (['u1', 'x1', 'x1p', 'y1', 'z1', 'z2'], ['u1', 'z1', 'z2'], ['u1']),
    'x2': (['u2', 'x2', 'x2p', 'y2', 'z1', 'z2'], ['u2', 'z1', 'z2'], ['u2']),
    'w1': (['w1', 'w1p', 'z1'], ['w1p', 'z1'], ['w1', 'w1p']),
    'w2': (['w2', 'w2p', 'z2'], ['w2p', 'z2'], ['w2', 'w2p']),
}
# ... the degree table of Section 2.
DEGREES = {'v': 2, 'u1': 5, 'u2': 5, 'z1': 7, 'z2': 7, 'y1': 5, 'y2': 5,
           'w1': 2, 'w2': 2, 'x1': 2, 'x2': 2,
           'w1p': 2, 'w2p': 2, 'x1p': 2, 'x2p': 2}
# ... the eccentricity table of Section 2.
ECC_TABLE = {'u1': 3, 'u2': 3, 'w1': 3, 'w2': 3, 'x1': 4, 'x2': 4,
             'y1': 3, 'y2': 3, 'z1': 2, 'z2': 2, 'v': 3,
             'w1p': 3, 'w2p': 3, 'x1p': 4, 'x2p': 4}
CUT_VERTICES = ['y1', 'y2', 'z1', 'z2']
# ... the minimum dominating broadcast of weight 5 printed in Section 3.
GAMMA_WITNESS = {'v': 1, 'x1': 1, 'x2': 1, 'w1': 1, 'w2': 1}


# ----------------------------------------------------------------------------
# graph machinery -- integers and bitmasks only
# ----------------------------------------------------------------------------
def build(names, edges):
    """-> (n, index map, distance matrix).  BFS distances; -1 if unreachable."""
    n = len(names)
    ix = {x: i for i, x in enumerate(names)}
    adj = [[0] * n for _ in range(n)]
    for a, b in edges:
        adj[ix[a]][ix[b]] = adj[ix[b]][ix[a]] = 1
    return n, ix, dists(n, adj), adj


def dists(n, adj):
    D = []
    for s in range(n):
        d = [-1] * n
        d[s] = 0
        q = [s]
        while q:
            nq = []
            for u in q:
                for w in range(n):
                    if adj[u][w] and d[w] < 0:
                        d[w] = d[u] + 1
                        nq.append(w)
            q = nq
        D.append(d)
    return D


def adj_from(n, edges):
    a = [[0] * n for _ in range(n)]
    for i, j in edges:
        a[i][j] = a[j][i] = 1
    return a


def decode_g6(s):
    """graph6 -> (order, set of index pairs i<j).  Upper triangle, column major."""
    n = ord(s[0]) - 63
    bits = ''.join(bin(ord(c) - 63)[2:].zfill(6) for c in s[1:])
    need = n * (n - 1) // 2
    if len(bits) < need:
        return n, None
    ed = set()
    p = 0
    for j in range(1, n):
        for i in range(j):
            if bits[p] == '1':
                ed.add((i, j))
            p += 1
    return n, ed


# ----------------------------------------------------------------------------
# broadcasts.  Definitions transcribed from the source paper:
#   a broadcast is f: V -> {0,...,diam G} with f(v) <= e(v);
#   N_f(v) = ball of radius f(v) about v, B_f(v) = sphere of radius f(v);
#   PB_f(v) = own(v) minus the union of N_f(z) over broadcasters z != v, where
#             own(v) = B_f(v) if f(v) >= 2 and own(v) = N_f(v) if f(v) = 1
#             (the source puts a broadcaster in its own private boundary
#             exactly when f(v) = 1);
#   f irredundant  <=> PB_f(v) non-empty for every broadcaster v;
#   f bn-independent <=> N_f(u) & N_f(v) subset of B_f(u) & B_f(v) for all
#             distinct broadcasters u, v;
#   Gamma_b   = max weight of a MINIMAL DOMINATING broadcast
#             = max weight of a dominating irredundant broadcast;
#   alpha_bn  = max weight of a bn-independent broadcast;
#   alpha_bnr = max weight of a bn-independent irredundant broadcast.
# ----------------------------------------------------------------------------
def tables(n, D):
    ecc = [max(D[s]) for s in range(n)]
    ball = [[sum(1 << w for w in range(n) if D[v][w] <= r) for r in range(ecc[v] + 1)]
            for v in range(n)]
    bnd = [[sum(1 << w for w in range(n) if D[v][w] == r) for r in range(ecc[v] + 1)]
           for v in range(n)]
    own = [[(bnd[v][r] if r >= 2 else ball[v][1]) for r in range(ecc[v] + 1)]
           for v in range(n)]
    return ecc, ball, bnd, own


def census(n, D, collect_irr=False):
    """Exhaustive census of the whole broadcast space of a connected graph.

    Depth-first over f(0), f(1), ..., f(n-1) with f(i) in 0..e(i).  TWO prunes,
    both MONOTONE and therefore complete:

      * irredundance: for a fixed broadcaster v, own(v) does not change as
        later vertices are assigned, while the union of the other N_f(z) only
        GROWS; so once PB_f(v) is empty on a partial assignment it is empty on
        every extension.  Any irredundant full broadcast therefore has every
        prefix irredundant, and no irredundant broadcast is missed.
      * bn-independence: the condition is a conjunction over PAIRS of
        broadcasters, and a violated pair stays violated under extension.

    The branch is abandoned only when BOTH flags have died, so the maxima over
    the irredundant, the bn-independent and the bn-independent-and-irredundant
    families are all exact.  Validated against an unpruned brute force below.
    """
    ecc, ball, bnd, own = tables(n, D)
    full = (1 << n) - 1
    res = {'gam': -1, 'gam_w': None, 'bnr': -1, 'bnr_w': None, 'bn': -1,
           'bn_w': None, 'n_irr': 0, 'n_mindom': 0, 'n_bn': 0, 'n_bnr': 0,
           'nodes': 0, 'irr': []}
    cur = [0] * n

    def rec(i, sup, w, cov, irr, bni):
        res['nodes'] += 1
        if i == n:
            if irr:
                res['n_irr'] += 1
                if collect_irr:
                    res['irr'].append(tuple(cur))
                if cov == full:
                    res['n_mindom'] += 1
                    if w > res['gam']:
                        res['gam'], res['gam_w'] = w, tuple(cur)
            if bni:
                res['n_bn'] += 1
                if w > res['bn']:
                    res['bn'], res['bn_w'] = w, tuple(cur)
                if irr:
                    res['n_bnr'] += 1
                    if w > res['bnr']:
                        res['bnr'], res['bnr_w'] = w, tuple(cur)
            return
        for r in range(ecc[i] + 1):
            cur[i] = r
            if r == 0:
                rec(i + 1, sup, w, cov, irr, bni)
                continue
            nsup = sup + [i]
            nirr = irr
            if irr:
                for j in nsup:
                    others = 0
                    for u in nsup:
                        if u != j:
                            others |= ball[u][cur[u]]
                    if own[j][cur[j]] & ~others == 0:
                        nirr = False
                        break
            nbni = bni
            if bni:
                for u in sup:
                    nu, nv = ball[u][cur[u]], ball[i][r]
                    if (nu & nv) & ~(bnd[u][cur[u]] & bnd[i][r]) != 0:
                        nbni = False
                        break
            if nirr or nbni:
                rec(i + 1, nsup, w + r, cov | ball[i][r], nirr, nbni)
        cur[i] = 0

    rec(0, [], 0, 0, True, True)
    return res


def census_brute(n, D):
    """The SAME quantities with NO pruning at all: every f in the product space."""
    ecc, ball, bnd, own = tables(n, D)
    full = (1 << n) - 1
    out = {'gam': -1, 'bnr': -1, 'bn': -1, 'n_mindom': 0, 'n_bn': 0,
           'n_bnr': 0, 'n_irr': 0, 'total': 0}
    for f in itertools.product(*[range(e + 1) for e in ecc]):
        out['total'] += 1
        sup = [v for v in range(n) if f[v] > 0]
        cov = 0
        for v in sup:
            cov |= ball[v][f[v]]
        irr = True
        for j in sup:
            others = 0
            for u in sup:
                if u != j:
                    others |= ball[u][f[u]]
            if own[j][f[j]] & ~others == 0:
                irr = False
                break
        bni = True
        for a in range(len(sup)):
            for b in range(a + 1, len(sup)):
                u, v = sup[a], sup[b]
                if (ball[u][f[u]] & ball[v][f[v]]) & ~(bnd[u][f[u]] & bnd[v][f[v]]):
                    bni = False
                    break
            if not bni:
                break
        w = sum(f)
        if irr:
            out['n_irr'] += 1
            if cov == full:
                out['n_mindom'] += 1
                out['gam'] = max(out['gam'], w)
        if bni:
            out['n_bn'] += 1
            out['bn'] = max(out['bn'], w)
            if irr:
                out['n_bnr'] += 1
                out['bnr'] = max(out['bnr'], w)
    return out


def names_of(mask, names):
    return sorted(names[i] for i in range(len(names)) if mask >> i & 1)


# ============================================================================
print('verification of the note: a 15-vertex counterexample to the '
      'minimum-degree-two question')
print('of Mynhardt and Neilson -- G\' with alpha_bnr = 6 > 5 = Gamma_b and delta = 2')
print('python %s, exact integer arithmetic only' % sys.version.split()[0])
T0 = time.time()

# ---------------------------------------------------------------- A. the object
print('')
print('=== A. the object exhibited in the paper')
n, ix, D, adj = build(VERTS, EDGES)
pairs = set(frozenset(e) for e in EDGES)
check('order_is_15', n == 15, 'n = %d' % n)
check('edge_list_is_26_distinct_pairs_of_distinct_vertices',
      len(EDGES) == 26 and len(pairs) == 26 and all(len(p) == 2 for p in pairs),
      'm = %d, distinct = %d' % (len(EDGES), len(pairs)))
check('every_endpoint_is_a_declared_vertex',
      all(a in ix and b in ix for a, b in EDGES), '15 declared names')
deg = {VERTS[i]: sum(1 for j in range(n) if D[i][j] == 1) for i in range(n)}
check('degree_table_matches_the_paper', deg == DEGREES,
      ' '.join('%s:%d' % (k, deg[k]) for k in VERTS))
check('degree_sum_is_52_i_e_26_edges', sum(deg.values()) == 52,
      'sum deg = %d = 2*%d' % (sum(deg.values()), len(EDGES)))
check('minimum_degree_is_2', min(deg.values()) == 2,
      'delta = %d, attained at %s' % (min(deg.values()),
                                      ','.join(sorted(k for k in deg if deg[k] == 2))))
check('graph_is_connected', all(D[0][j] >= 0 for j in range(n)))
ecc = {VERTS[i]: max(D[i]) for i in range(n)}
check('eccentricity_table_matches_the_paper', ecc == ECC_TABLE,
      ' '.join('%s:%d' % (k, ecc[k]) for k in VERTS))
check('diameter_is_4', max(ecc.values()) == 4, 'diam = %d' % max(ecc.values()))

g6n, g6e = decode_g6(G6)
mine = set(tuple(sorted((ix[a], ix[b]))) for a, b in EDGES)
check('graph6_string_decodes_to_the_same_LABELLED_graph',
      g6n == n and g6e == mine,
      'graph6 %s -> order %d, %d edges, label-equal to the printed edge list'
      % (G6, g6n, len(g6e or ())))


def cut_vertices(nn, dd, aa):
    out = []
    for s in range(nn):
        rest = [i for i in range(nn) if i != s]
        sub = dists(nn, [[aa[i][j] if i != s and j != s else 0 for j in range(nn)]
                         for i in range(nn)])
        if any(sub[rest[0]][j] < 0 for j in rest):
            out.append(s)
    return out


cuts = sorted(VERTS[i] for i in cut_vertices(n, D, adj))
check('cut_vertices_are_exactly_y1_y2_z1_z2', cuts == sorted(CUT_VERTICES),
      'cut vertices %s -- so G\' is NOT 2-connected' % ' '.join(cuts))
space = functools.reduce(lambda a, b: a * b, [ecc[t] + 1 for t in VERTS])
check('broadcast_space_has_1474560000_elements', space == 1474560000,
      'prod_t (e(t)+1) = %d' % space)

# ------------------------------------------------- B. the doubling structure
print('')
print('=== B. G\' is G_1 with its four end-vertices doubled into true twins')
V1 = [t for t in VERTS if t not in [p for p, _ in TWINS]]
E1 = [(a, b) for a, b in EDGES if a in V1 and b in V1]
n1, ix1, D1, adj1 = build(V1, E1)
check('deleting_the_four_primed_vertices_leaves_11_vertices_and_18_edges',
      n1 == 11 and len(E1) == 18, 'G_1: n = %d, m = %d' % (n1, len(E1)))
check('G_1_is_connected', all(D1[0][j] >= 0 for j in range(n1)))
deg1 = {V1[i]: sum(1 for j in range(n1) if D1[i][j] == 1) for i in range(n1)}
ends = sorted(t for t in V1 if deg1[t] == 1)
check('the_end_vertices_of_G_1_are_exactly_the_four_doubled_vertices',
      ends == sorted(q for _, q in TWINS),
      'end-vertices of G_1: %s' % ' '.join(ends))
ok = True
for p, q in TWINS:
    cp = set(names_of(sum(1 << ix[t] for t in VERTS if D[ix[p]][ix[t]] <= 1), VERTS))
    cq = set(names_of(sum(1 << ix[t] for t in VERTS if D[ix[q]][ix[t]] <= 1), VERTS))
    ok = ok and cp == cq and len(cp) == 3
check('each_primed_vertex_is_a_true_twin_closed_neighbourhoods_equal', ok,
      '; '.join('N[%s] = N[%s] = {%s}'
                % (p, q, ','.join(names_of(sum(1 << ix[t] for t in VERTS
                                               if D[ix[q]][ix[t]] <= 1), VERTS)))
                for p, q in TWINS))
check('the_doubling_changes_no_distance_inside_V_G_1',
      all(D[ix[a]][ix[b]] == D1[ix1[a]][ix1[b]] for a in V1 for b in V1),
      '11x11 distance matrices agree entrywise')
check('the_doubling_changes_no_eccentricity_inside_V_G_1',
      all(max(D[ix[a]]) == max(D1[ix1[a]]) for a in V1),
      'so every broadcast of G_1 is a broadcast of G\' and conversely')
check('each_primed_vertex_inherits_the_eccentricity_of_its_partner',
      all(max(D[ix[p]]) == max(D[ix[q]]) for p, q in TWINS),
      '; '.join('e(%s) = e(%s) = %d' % (p, q, max(D[ix[p]])) for p, q in TWINS))

# ------------------------------------------------------------- C. the witness
print('')
print('=== C. the witness broadcast f of Section 3 (the hand lower bound)')
ecc_l, ball, bnd, own = tables(n, D)
sup = sorted(WITNESS, key=lambda t: ix[t])
check('f_is_a_legal_broadcast_f_t_at_most_e_t',
      all(WITNESS[t] <= ecc[t] for t in WITNESS),
      ', '.join('f(%s)=%d <= e(%s)=%d' % (t, WITNESS[t], t, ecc[t]) for t in sup))
check('the_weight_of_f_is_6', sum(WITNESS.values()) == 6,
      'sigma(f) = %d' % sum(WITNESS.values()))
allmask = 0
for t in sup:
    allmask |= ball[ix[t]][WITNESS[t]]
for t in sup:
    i = ix[t]
    r = WITNESS[t]
    others = 0
    for u in sup:
        if u != t:
            others |= ball[ix[u]][WITNESS[u]]
    Nf = names_of(ball[i][r], VERTS)
    Bf = names_of(bnd[i][r], VERTS)
    Pf = names_of(own[i][r] & ~others, VERTS)
    want = WITNESS_SETS[t]
    check('N_f_B_f_and_PB_f_of_%s_match_the_paper' % t,
          (Nf, Bf, Pf) == (sorted(want[0]), sorted(want[1]), sorted(want[2])),
          'f(%s)=%d  N_f={%s}  B_f={%s}  PB_f={%s}'
          % (t, r, ','.join(Nf), ','.join(Bf), ','.join(Pf)))
check('f_is_irredundant_every_private_boundary_is_non_empty',
      all((own[ix[t]][WITNESS[t]]
           & ~functools.reduce(lambda a, b: a | b,
                               [ball[ix[u]][WITNESS[u]] for u in sup if u != t], 0)) != 0
          for t in sup), 'all 4 private boundaries non-empty')
bad = [(u, w) for a, u in enumerate(sup) for w in sup[a + 1:]
       if (ball[ix[u]][WITNESS[u]] & ball[ix[w]][WITNESS[w]])
       & ~(bnd[ix[u]][WITNESS[u]] & bnd[ix[w]][WITNESS[w]])]
check('f_is_bn_independent', not bad,
      'all %d pairs satisfy N_f(u) & N_f(w) subset of B_f(u) & B_f(w)'
      % (len(sup) * (len(sup) - 1) // 2))
dd = {('x1', 'x2'): 4, ('x1', 'w1'): 3, ('x1', 'w2'): 3,
      ('x2', 'w1'): 3, ('x2', 'w2'): 3, ('w1', 'w2'): 3}
check('the_pairwise_distances_of_the_support_match_the_paper',
      all(D[ix[a]][ix[b]] == v for (a, b), v in dd.items())
      and all(D[ix[a]][ix[b]] >= WITNESS[a] + WITNESS[b] for a, b in dd),
      ', '.join('d(%s,%s)=%d>=%d' % (a, b, D[ix[a]][ix[b]], WITNESS[a] + WITNESS[b])
                for a, b in sorted(dd)))
unc = names_of(((1 << n) - 1) & ~allmask, VERTS)
check('f_is_NOT_dominating_and_the_only_unheard_vertex_is_v', unc == ['v'],
      'uncovered = {%s}; so f does not compete in Gamma_b' % ','.join(unc))
check('alpha_bnr_of_G_prime_is_at_least_6_by_this_witness_alone',
      sum(WITNESS.values()) == 6,
      'f bn-independent + irredundant of weight 6 -- no computation needed')

gw = sorted(GAMMA_WITNESS, key=lambda t: ix[t])
gcov = 0
for t in gw:
    gcov |= ball[ix[t]][GAMMA_WITNESS[t]]
gpb = {}
for t in gw:
    o = 0
    for u in gw:
        if u != t:
            o |= ball[ix[u]][GAMMA_WITNESS[u]]
    gpb[t] = names_of(own[ix[t]][GAMMA_WITNESS[t]] & ~o, VERTS)
check('the_weight_5_broadcast_of_Section_4_is_minimal_dominating',
      sum(GAMMA_WITNESS.values()) == 5 and gcov == (1 << n) - 1
      and all(gpb[t] for t in gw),
      'sigma = 5, dominating, ' + '; '.join('PB(%s)={%s}' % (t, ','.join(gpb[t]))
                                            for t in gw))

# ------------------------------------------------------------- D/E. the census
print('')
print('=== E. the census enumerator, validated against an UNPRUNED brute force')
t = time.time()
r1 = census(n1, D1)
b1 = census_brute(n1, D1)
check('pruned_and_unpruned_censuses_of_G_1_agree_on_all_seven_quantities',
      all(r1[k] == b1[k] for k in ('gam', 'bnr', 'bn', 'n_mindom', 'n_bn',
                                   'n_bnr', 'n_irr')),
      'brute force visited all %d broadcasts of G_1 (%d DFS nodes pruned to); '
      'Gamma_b=%d alpha_bnr=%d alpha_bn=%d, counts %d/%d/%d/%d'
      % (b1['total'], r1['nodes'], r1['gam'], r1['bnr'], r1['bn'],
         r1['n_mindom'], r1['n_bnr'], r1['n_bn'], r1['n_irr']))
note('the brute force took %.1f s; it is the only check that touches the whole '
     'product space of a graph in this note' % (time.time() - t))

CONTROLS = [('P6', 6, [(i, i + 1) for i in range(5)], 5, 5, 5),
            ('P7', 7, [(i, i + 1) for i in range(6)], 6, 6, 6),
            ('C6', 6, [(i, (i + 1) % 6) for i in range(6)], 4, 3, 3),
            ('grid3x3', 9, [(i * 3 + j, i * 3 + j + 1) for i in range(3) for j in range(2)]
             + [(i * 3 + j, i * 3 + j + 3) for i in range(2) for j in range(3)], 6, 5, 5),
            ('grid3x4', 12, [(i * 4 + j, i * 4 + j + 1) for i in range(3) for j in range(3)]
             + [(i * 4 + j, i * 4 + j + 4) for i in range(2) for j in range(4)], 9, 6, 6)]
for nm, nn, ee, g, bnr, bn in CONTROLS:
    rr = census(nn, dists(nn, adj_from(nn, ee)))
    check('control_%s_reproduces_the_published_values' % nm,
          (rr['gam'], rr['bnr'], rr['bn']) == (g, bnr, bn),
          'Gamma_b=%d alpha_bnr=%d alpha_bn=%d' % (rr['gam'], rr['bnr'], rr['bn']))

print('')
print('=== D. the exhaustive census of G_1 and of G\'')
check('Gamma_b_of_G_1_is_5_the_published_value_2k+3_at_k=1', r1['gam'] == 5,
      'Gamma_b(G_1) = %d, %d minimal dominating broadcasts' % (r1['gam'], r1['n_mindom']))
check('alpha_bnr_of_G_1_is_6_the_published_value', r1['bnr'] == 6,
      'alpha_bnr(G_1) = %d > 5 = Gamma_b(G_1), but delta(G_1) = %d'
      % (r1['bnr'], min(deg1.values())))
check('G_1_has_minimum_degree_1_so_it_is_not_itself_an_answer',
      min(deg1.values()) == 1, 'delta(G_1) = 1')
t = time.time()
R = census(n, D, collect_irr=True)
note('the census of G\' visited %d DFS nodes in %.2f s and found %d irredundant '
     'broadcasts' % (R['nodes'], time.time() - t, R['n_irr']))
check('Gamma_b_of_G_prime_is_5', R['gam'] == 5, 'Gamma_b(G\') = %d' % R['gam'])
check('G_prime_has_exactly_283_minimal_dominating_broadcasts',
      R['n_mindom'] == 283, 'n = %d' % R['n_mindom'])
check('alpha_bnr_of_G_prime_is_6', R['bnr'] == 6, 'alpha_bnr(G\') = %d' % R['bnr'])
check('G_prime_has_exactly_1140_bn_independent_irredundant_broadcasts',
      R['n_bnr'] == 1140, 'n = %d' % R['n_bnr'])
check('alpha_bn_of_G_prime_is_7', R['bn'] == 7, 'alpha_bn(G\') = %d' % R['bn'])
check('G_prime_has_exactly_1416_bn_independent_broadcasts',
      R['n_bn'] == 1416, 'n = %d' % R['n_bn'])
check('Gamma_b_is_unchanged_by_the_four_doublings',
      R['gam'] == r1['gam'] == 5,
      'Gamma_b(G_1) = Gamma_b(G\') = 5, the conclusion of the doubling lemma at '
      'this graph')
mw = {VERTS[i]: R['gam_w'][i] for i in range(n) if R['gam_w'][i]}
check('the_optimal_Gamma_b_broadcast_found_is_the_twin_image_of_the_printed_one',
      sorted(mw.values()) == [1] * 5
      and sorted(mw) == sorted(['v', 'w1p', 'w2p', 'x1p', 'x2p']),
      ' '.join('%s:%d' % (k, mw[k]) for k in sorted(mw)))
bw = {VERTS[i]: R['bnr_w'][i] for i in range(n) if R['bnr_w'][i]}
check('the_optimal_alpha_bnr_broadcast_found_is_the_twin_image_of_the_witness',
      bw == {'x1p': 2, 'x2p': 2, 'w1p': 1, 'w2p': 1},
      ' '.join('%s:%d' % (k, bw[k]) for k in sorted(bw)))
tw = [(ix[p], ix[q]) for p, q in TWINS]
check('no_irredundant_broadcast_of_G_prime_weights_both_of_a_twin_pair',
      not any(f[a] > 0 and f[b] > 0 for f in R['irr'] for a, b in tw),
      'checked over all %d irredundant broadcasts and all 4 twin pairs '
      '-- the key step of the doubling lemma' % R['n_irr'])
check('THE_RESULT_alpha_bnr_exceeds_Gamma_b_at_minimum_degree_2',
      R['bnr'] > R['gam'] and min(deg.values()) == 2,
      'alpha_bnr(G\') = %d > %d = Gamma_b(G\'), delta(G\') = 2, and G\' is not '
      '2-connected' % (R['bnr'], R['gam']))

# --------------------------------------------------- F. the doubling sweep
print('')
print('=== F. the end-vertex doubling lemma, swept over all small graphs')
t = time.time()
counts, pairs_n, mism = {}, 0, []
for nn in range(3, 7):
    pos = [(i, j) for i in range(nn) for j in range(i + 1, nn)]
    cg = 0
    for mask in range(1 << len(pos)):
        es = [pos[k] for k in range(len(pos)) if mask >> k & 1]
        aa = adj_from(nn, es)
        DD = dists(nn, aa)
        if any(DD[0][j] < 0 for j in range(nn)):
            continue
        cg += 1
        endv = [a for a in range(nn) if sum(aa[a]) == 1]
        if not endv:
            continue
        gH = census(nn, DD)['gam']
        for p in endv:
            q = aa[p].index(1)
            D2 = dists(nn + 1, adj_from(nn + 1, es + [(p, nn), (q, nn)]))
            g2 = census(nn + 1, D2)['gam']
            pairs_n += 1
            if g2 != gH:
                mism.append((nn, mask, p, gH, g2))
    counts[nn] = cg
check('the_connected_labelled_graph_counts_match_OEIS_A001187',
      [counts[k] for k in (3, 4, 5, 6)] == [4, 38, 728, 26704],
      'n=3..6: %s, total %d' % (', '.join(str(counts[k]) for k in (3, 4, 5, 6)),
                                sum(counts.values())))
check('Gamma_b_is_unchanged_by_doubling_an_end_vertex_on_every_small_instance',
      pairs_n == 22654 and not mism,
      '%d (H,p) pairs with H connected on 3..6 vertices, %d mismatches, %.1f s'
      % (pairs_n, len(mism), time.time() - t))

# ------------------------------------------------------------------ scope
print('')
note('SCOPE -- what this program does NOT establish.  (1) It settles only the '
     'SECOND sentence of Question 4 of Mynhardt and Neilson, the delta(G) >= 2 '
     'one, in its repaired reading; the FIRST sentence, for 2-connected G, is '
     'untouched and remains open, and G\' is deliberately not 2-connected (four '
     'cut vertices, checked above). (2) NOT RE-RUN: any claim of minimality -- no '
     'census of graphs of order below 15 with delta >= 2 is performed here, so '
     'nothing above says G\' is a smallest witness. (3) NOT RE-RUN: the family '
     'G_k for k >= 2 -- this program reads only the 15-vertex object printed in '
     'the paper and the 11-vertex graph G_1 obtained from it, so the unbounded-gap '
     'statement is outside its scope. (4) GAPS NOT COVERED: the definitions of '
     'bn-independence, irredundance, alpha_bnr and Gamma_b, and the '
     'characterisation of a minimal dominating broadcast as a dominating '
     'irredundant one, are TRANSCRIBED from the source paper and not re-derived; '
     'the evidence that they were transcribed correctly is that the census '
     'reproduces the source\'s own published Gamma_b(G_1) = 5 and '
     'alpha_bnr(G_1) = 6 and the five independent published control values above. '
     '(5) The census of G\' uses two monotone prunes rather than a full walk of '
     'its 1,474,560,000 broadcasts; the prunes are validated by exact agreement '
     'with an unpruned walk of all 3,686,400 broadcasts of G_1 and with all five '
     'controls, but no unpruned walk of G\' itself is performed here. (6) The '
     'end-vertex doubling sweep is exhaustive only for H on at most 6 vertices; '
     'the general lemma rests on the hand proof in the paper, not on the sweep.')
print('')
note('elapsed %.1f s' % (time.time() - T0))
if _N_FAIL:
    print('VERDICT: %d CHECK(S) FAILED of %d' % (_N_FAIL, _N_FAIL + _N_PASS))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _N_PASS)
sys.exit(0)

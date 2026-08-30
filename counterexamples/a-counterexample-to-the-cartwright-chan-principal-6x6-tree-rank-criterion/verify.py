#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- checks the claims of the note

    "A Counterexample to Cartwright--Chan's Principal 6x6 Criterion for Tree Rank"

Python 3.9+, STANDARD LIBRARY ONLY (itertools, functools, sys, collections). No numpy, no
solver, no external data file, no randomness, and no floating point in any decision: every
comparison below is between exact integers or between sets of integers.

It reads the object AS PRINTED IN THE PAPER -- the 7x7 matrix is transcribed from the paper's
display, and so are the three-element cover, the seven two-element covers, the edge and
non-edge lists, and every number the paper states -- and re-derives all of them.

Prints one `PASS <name> [detail]` line per check and closes with
    VERDICT: ALL <n> CHECKS PASS
exiting 0 if and only if every check passed. What is NOT covered is stated in the output
itself, on the `NOT RE-RUN:` lines at the end.
"""

import itertools
import sys
from collections import Counter
from functools import lru_cache

# =====================================================================================
# PART 0.  EVERYTHING THE PAPER PRINTS, TRANSCRIBED.  Nothing below this block is input.
# =====================================================================================

# The matrix M of the paper's display (rows and columns 1..7, `*` on the absent diagonal).
PAPER_MATRIX = """
* 0 1 0 0 0 0
0 * 1 0 0 0 0
1 1 * 0 0 0 0
0 0 0 * 1 1 0
0 0 0 1 * 0 1
0 0 0 1 0 * 1
0 0 0 0 1 1 *
"""

PAPER_EDGES = "12 14 15 16 17 24 25 26 27 34 35 36 37 47 56"      # E(G), fifteen of them
PAPER_NONEDGES = "13 23 45 46 57 67"                              # the six one-entries of M
PAPER_A, PAPER_B = [1, 2, 3], [4, 5, 6, 7]                        # the paper's bipartition
PAPER_EXTRA_A, PAPER_EXTRA_B = ["12"], ["47", "56"]               # edges inside A, inside B

# The three-element cover of the upper bound, written as its parts.
PAPER_COVER3 = [
    [[1, 2, 3], [4, 5, 6, 7]],
    [[1], [2], [4], [7]],
    [[5], [6]],
]
PAPER_COVER3_SIZES = [12, 6, 1]          # edges contributed by each member, as printed
PAPER_COVER3_SLOTS = 19                  # 12 + 6 + 1
PAPER_COVER3_DUPLICATES = "14 17 24 27"  # the four repeats, as printed
PAPER_COVER3_UNION = 15                  # 12 + 2 + 1

# The seven two-element covers for the principal 6x6 submatrices, written as their parts.
PAPER_COVER2 = {
    1: ([[4], [7], [2, 3]], [[5], [6], [2, 3]]),
    2: ([[4], [7], [1, 3]], [[5], [6], [1, 3]]),
    3: ([[4], [7], [1, 2]], [[1], [2], [5], [6]]),
    4: ([[1, 2, 3], [5, 6, 7]], [[1], [2], [5], [6]]),
    5: ([[1, 2, 3], [4, 6, 7]], [[1], [2], [4], [7]]),
    6: ([[1, 2, 3], [4, 5, 7]], [[1], [2], [4], [7]]),
    7: ([[1, 2, 3], [4, 5, 6]], [[1], [2], [5], [6]]),
}
PAPER_DELETION_EDGE_COUNTS = {1: 10, 2: 10, 3: 11, 4: 11, 5: 11, 6: 11, 7: 11}

PAPER_TREE_RANK_M = 3                    # tree rank of M
PAPER_TREE_RANK_PRINCIPAL = 2            # tree rank of every principal 6x6 submatrix of M
PAPER_MIN_DEGREE_G = 4                   # the paper's proof of the uniform theorem uses this
PAPER_UNIFORM_FROM, PAPER_UNIFORM_TO = 8, 14   # proved for all n >= 7; re-run here to n = 14

# Values quoted from the source paper (Cartwright--Chan), used here as controls.
SOURCE_MAX_TREE_RANK = {3: 1, 4: 2, 5: 3, 6: 3}     # their table of maxima
SOURCE_C5_TREE_RANK = 3                             # their value for the 5-cycle
SOURCE_CM_COUNTS = {5: 171, 6: 813, 7: 4012}        # nonempty complete multipartite graphs

# The census of the whole 0/1 cell at n = 7, as the paper reports it.
PAPER_CELL_SIZE = 2097152           # 2^21
PAPER_R_LE_1 = 4013                 # complete multipartite graphs on [7], the empty one included
PAPER_R_LE_2 = 502511
PAPER_TR7_LE_2 = 486761
PAPER_TR7_GE_3 = 1610391
PAPER_TR6_LE_2 = 21044              # of 2^15
PAPER_ALL_SEVEN_OK = 488336
PAPER_COUNTEREXAMPLES = 1575
PAPER_MONOTONICITY_VIOLATIONS = 0
PAPER_CLASSES = 2
PAPER_ORBIT_SIZES = [315, 1260]
PAPER_OTHER_CLASS_ONES = 5          # the second class has five one-entries; M has six
PAPER_OTHER_CLASS_IS_G_PLUS = (4, 5)   # ... and its zero-entry graph is G with 45 added
PAPER_CM_SUBGRAPHS_OF_G = 209       # nonempty complete multipartite subgraphs of G
PAPER_CM_SUBGRAPHS_AT_8 = 1004      # ... and of the order-8 member of the uniform family
PAPER_N5_TREE_RANK_3 = 12           # labelled objects of tree rank 3 at n = 5
PAPER_N5_CLASSES = 1                # ... forming one class, namely the 5-cycle
PAPER_A000088_7 = 1044              # graphs on 7 unlabelled nodes, OEIS A000088

# =====================================================================================
# PART 1.  the harness
# =====================================================================================

_PASSED = 0
_FAILED = 0


def ck(name, cond, detail=''):
    global _PASSED, _FAILED
    if cond:
        _PASSED += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILED += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(s):
    print('NOTE %s' % s)


def head(s):
    print()
    print('=== %s' % s)


# =====================================================================================
# PART 2.  graphs as edge bitmasks over the lexicographic pair order
# =====================================================================================

@lru_cache(maxsize=None)
def pair_list(n):
    """the C(n,2) pairs of {1,...,n} in lexicographic order; bit k of a mask is pair k."""
    return tuple((i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1))


@lru_cache(maxsize=None)
def pair_index(n):
    return {p: k for k, p in enumerate(pair_list(n))}


@lru_cache(maxsize=None)
def incidence(n):
    inc = [0] * (n + 1)
    for k, (i, j) in enumerate(pair_list(n)):
        inc[i] |= 1 << k
        inc[j] |= 1 << k
    return tuple(inc)


def bit(n, i, j):
    return 1 << pair_index(n)[(min(i, j), max(i, j))]


def mask_of(n, edges):
    ix = pair_index(n)
    m = 0
    for (i, j) in edges:
        m |= 1 << ix[(min(i, j), max(i, j))]
    return m


def edges_of(n, mask):
    return [p for k, p in enumerate(pair_list(n)) if (mask >> k) & 1]


def label(edges):
    return ','.join('%d%d' % e for e in sorted(edges))


def parse_edge_words(s):
    """'12 14 15' -> [(1,2),(1,4),(1,5)]; only single-digit vertex labels occur in the paper."""
    return [(int(w[0]), int(w[1])) for w in s.replace(',', ' ').split()]


def isolated_vertices(n, mask, verts=None):
    inc = incidence(n)
    return [v for v in (verts if verts is not None else range(1, n + 1)) if not (mask & inc[v])]


def popcount(x):
    return bin(x).count('1')


def neighbourhoods(n, mask, verts=None):
    """N(v) as a bitmask over vertex labels, restricted to `verts`."""
    V = tuple(verts if verts is not None else range(1, n + 1))
    ix = pair_index(n)
    nbr = {v: 0 for v in V}
    for a in range(len(V)):
        for b in range(a + 1, len(V)):
            u, v = V[a], V[b]
            if (mask >> ix[(min(u, v), max(u, v))]) & 1:
                nbr[u] |= 1 << v
                nbr[v] |= 1 << u
    return V, nbr


def is_complete_multipartite(n, mask, verts=None):
    """Cartwright--Chan's characterisation, in the equivalent neighbourhood form.

    Their form: AMONG THE VERTICES INCIDENT TO SOME EDGE, non-adjacency is transitive (vertices
    in no part are isolated and are ignored).  Equivalent form used here: for all distinct
    non-isolated u, v, if u and v are non-adjacent then they have the same neighbourhood.
    (=>) In a complete multipartite graph non-adjacent vertices lie in one part and so have the
    same neighbourhood.  (<=) If the implication holds then non-adjacency is transitive on the
    non-isolated vertices -- u!~v, v!~w and u~w would put w in N(u)=N(v), i.e. v~w -- and the
    classes of that equivalence are parts between which every pair is an edge.
    The two forms are checked against each other exhaustively at n = 5 below.
    """
    V, nbr = neighbourhoods(n, mask, verts)
    W = [v for v in V if nbr[v]]
    for a in range(len(W)):
        for b in range(a + 1, len(W)):
            u, v = W[a], W[b]
            if not (nbr[u] >> v) & 1 and nbr[u] != nbr[v]:
                return False
    return True


def is_complete_multipartite_literal(n, mask, verts=None):
    """the source's wording, coded literally: transitivity of non-adjacency on the non-isolated
    vertices.  Quadratically slower, so used only as a control on the routine above."""
    V, nbr = neighbourhoods(n, mask, verts)
    W = [v for v in V if nbr[v]]

    def nonadj(x, y):
        return x != y and not (nbr[x] >> y) & 1

    for u in W:
        for v in W:
            if not nonadj(u, v):
                continue
            for w in W:
                if not nonadj(v, w) or w == u:
                    continue
                if not nonadj(u, w):
                    return False
    return True


def parts_to_mask(n, parts):
    """the complete multipartite graph with these DISJOINT parts, as an edge mask."""
    seen = set()
    for p in parts:
        for v in p:
            if v in seen:
                raise ValueError('parts are not disjoint')
            seen.add(v)
    m = 0
    ix = pair_index(n)
    for a in range(len(parts)):
        for b in range(a + 1, len(parts)):
            for u in parts[a]:
                for v in parts[b]:
                    m |= 1 << ix[(min(u, v), max(u, v))]
    return m


def set_partitions(lst):
    if not lst:
        yield []
        return
    first, rest = lst[0], lst[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


@lru_cache(maxsize=None)
def cm_masks_by_definition(n):
    """Every complete multipartite graph on {1,...,n} as an edge mask, from the DEFINITION:
    disjoint sets I_1,...,I_k with k >= 2 and an edge between elements of distinct sets.  The
    empty graph is included as mask 0."""
    out = {0}
    for size in range(2, n + 1):
        for A in itertools.combinations(range(1, n + 1), size):
            for part in set_partitions(list(A)):
                if len(part) >= 2:
                    out.add(parts_to_mask(n, part))
    return tuple(sorted(out))


@lru_cache(maxsize=None)
def cover_layers(n, kmax):
    """L[k] = the masks that are unions of at most k complete multipartite subgraphs."""
    cm = cm_masks_by_definition(n)
    cur = {0}
    L = [frozenset(cur)]
    for _ in range(kmax):
        nxt = set(cur)
        for a in cur:
            nxt.update(a | c for c in cm)
        cur = nxt
        L.append(frozenset(cur))
    return cm, tuple(L)


def cover_number(n, mask, kmax=3):
    """r(G): the least number of complete multipartite SUBGRAPHS of G whose union is E(G).
    A breadth-first search with an early exit; used only where the answer is small."""
    if mask == 0:
        return 0
    sub = [c for c in cm_masks_by_definition(n) if c and not (c & ~mask)]
    reach = {0}
    for k in range(1, kmax + 1):
        nxt = set()
        for a in reach:
            for c in sub:
                u = a | c
                if u == mask:
                    return k
                nxt.add(u)
        reach = nxt
    return None


def tree_rank_small(n, mask, L):
    """Cartwright--Chan's proposition for 0/1 dissimilarity matrices: the tree rank is r, except
    that it is r + 1 when the zero-entry graph has at least two isolated vertices."""
    r = None
    for k in range(len(L)):
        if mask in L[k]:
            r = k
            break
    if r is None:
        return None
    return r + (1 if len(isolated_vertices(n, mask)) >= 2 else 0)


def bitset(size, masks):
    """one BYTE per mask position, packed into a single integer, so that whole-cell set algebra
    runs at C speed; every byte is 0 or 1, so &, | and ^ act pointwise."""
    ba = bytearray(size)
    for m in masks:
        ba[m] = 1
    return int.from_bytes(bytes(ba), 'little')


def bell(k):
    b = [[0] * (k + 1) for _ in range(k + 1)]
    b[0][0] = 1
    for i in range(1, k + 1):
        b[i][0] = b[i - 1][i - 1]
        for j in range(1, i + 1):
            b[i][j] = b[i - 1][j - 1] + b[i][j - 1]
    return b[k][0]


def binom(a, b):
    if b < 0 or b > a:
        return 0
    r = 1
    for t in range(b):
        r = r * (a - t) // (t + 1)
    return r


def orbits_under_relabelling(n, masks):
    """the orbits of `masks` under the symmetric group on the n vertices, and whether the set is
    closed; generated by the n-1 adjacent transpositions."""
    ix = pair_index(n)
    gens = []
    for t in range(1, n):
        perm = {v: v for v in range(1, n + 1)}
        perm[t], perm[t + 1] = t + 1, t
        tab = [0] * len(ix)
        for (a, b), k in ix.items():
            u, v = perm[a], perm[b]
            tab[k] = 1 << ix[(min(u, v), max(u, v))]
        gens.append(tab)

    def act(tab, g):
        m, z = 0, g
        while z:
            lb = z & -z
            m |= tab[lb.bit_length() - 1]
            z ^= lb
        return m

    pool, seen, out = set(masks), set(), []
    for g in sorted(pool):
        if g in seen:
            continue
        orb, stack = {g}, [g]
        while stack:
            x = stack.pop()
            for tab in gens:
                y = act(tab, x)
                if y not in orb:
                    orb.add(y)
                    stack.append(y)
        seen |= orb
        out.append(sorted(orb))
    return out, seen == pool


# =====================================================================================
# STEP 1.  the object, read off the paper's printed matrix
# =====================================================================================

head('Step 1: the matrix M printed in the paper, and its zero-entry graph G')

rows = [r.split() for r in PAPER_MATRIX.strip().splitlines()]
ck('printed_matrix_is_7_by_7', len(rows) == 7 and all(len(r) == 7 for r in rows),
   'n=%d' % len(rows))
ck('printed_matrix_has_a_star_exactly_on_the_diagonal',
   all((rows[i][j] == '*') == (i == j) for i in range(7) for j in range(7)),
   'a dissimilarity matrix has no diagonal entry')
ck('printed_off_diagonal_entries_are_all_0_or_1',
   all(rows[i][j] in ('0', '1') for i in range(7) for j in range(7) if i != j))
ck('printed_matrix_is_symmetric',
   all(rows[i][j] == rows[j][i] for i in range(7) for j in range(7)))

G_edges = [(i + 1, j + 1) for i in range(7) for j in range(i + 1, 7) if rows[i][j] == '0']
G_ones = [(i + 1, j + 1) for i in range(7) for j in range(i + 1, 7) if rows[i][j] == '1']
G = mask_of(7, G_edges)
note('zero-entry graph G read off the printed matrix: {%s}' % label(G_edges))
note('one-entries of M: {%s}' % label(G_ones))

ck('zero_entry_graph_matches_the_edge_list_printed_in_the_paper',
   sorted(G_edges) == sorted(parse_edge_words(PAPER_EDGES)), '%d edges' % len(G_edges))
ck('one_entries_match_the_non_edge_list_printed_in_the_paper',
   sorted(G_ones) == sorted(parse_edge_words(PAPER_NONEDGES)), '{%s}' % label(G_ones))
ck('edges_and_non_edges_partition_all_21_pairs',
   len(G_edges) + len(G_ones) == 21 and not (set(G_edges) & set(G_ones))
   and set(G_edges) | set(G_ones) == set(pair_list(7)),
   '15 + 6 = 21, intersection empty')

cross = [(a, b) for a in PAPER_A for b in PAPER_B]
ck('G_is_complete_bipartite_between_A_and_B_plus_the_printed_inside_edges',
   set(G_edges) == set(cross) | set(parse_edge_words(' '.join(PAPER_EXTRA_A + PAPER_EXTRA_B))),
   'A={%s} B={%s}: %d cross edges, plus %s inside A, plus %s inside B'
   % (','.join(map(str, PAPER_A)), ','.join(map(str, PAPER_B)), len(cross),
      ','.join(PAPER_EXTRA_A), ','.join(PAPER_EXTRA_B)))

deg = {v: sum(1 for e in G_edges if v in e) for v in range(1, 8)}
note('degrees in G: %s' % dict(sorted(deg.items())))
ck('G_has_no_isolated_vertex', not isolated_vertices(7, G), 'minimum degree %d' % min(deg.values()))
ck('minimum_degree_of_G_is_the_value_the_paper_uses',
   min(deg.values()) == PAPER_MIN_DEGREE_G, 'delta(G) = %d' % min(deg.values()))

Gdel = {}
for i in range(1, 8):
    Gdel[i] = (mask_of(7, [e for e in G_edges if i not in e]),
               tuple(v for v in range(1, 8) if v != i))
ck('no_principal_6x6_zero_entry_graph_has_an_isolated_vertex_either',
   all(not isolated_vertices(7, m, V) for m, V in Gdel.values()),
   'so the +1 branch of the source proposition never fires, for G or for any G-i')
ck('edge_counts_of_the_seven_deletions_are_as_printed',
   all(popcount(Gdel[i][0]) == PAPER_DELETION_EDGE_COUNTS[i] for i in range(1, 8)),
   ' '.join('|E(G-%d)|=%d' % (i, popcount(Gdel[i][0])) for i in range(1, 8)))

# =====================================================================================
# STEP 2.  the complete-multipartite machinery, validated three ways
# =====================================================================================

head('Step 2: complete multipartite graphs -- definition, transitivity test, closed form')

CM = {n: cm_masks_by_definition(n) for n in (3, 4, 5, 6, 7, 8)}
for n in (5, 6, 7):
    ck('definition_enumeration_counts_the_nonempty_complete_multipartite_graphs_at_n_%d' % n,
       len(CM[n]) - 1 == SOURCE_CM_COUNTS[n],
       '%d nonempty; the closed form below and the source agree at %d'
       % (len(CM[n]) - 1, SOURCE_CM_COUNTS[n]))
    closed = sum(binom(n, a) * (bell(a) - 1) for a in range(2, n + 1))
    ck('closed_form_sum_of_C_n_a_times_Bell_a_minus_1_agrees_at_n_%d' % n,
       closed == len(CM[n]) - 1, 'sum_{a>=2} C(%d,a)(Bell(a)-1) = %d' % (n, closed))

for n in (5, 6):
    brute = tuple(sorted(g for g in range(1 << len(pair_list(n)))
                         if is_complete_multipartite(n, g)))
    ck('the_recogniser_agrees_with_the_definition_over_every_graph_at_n_%d' % n,
       brute == CM[n], '%d masks of 2^%d, identical lists' % (len(brute), len(pair_list(n))))

ck('the_recogniser_agrees_with_the_sources_literal_transitivity_wording_at_n_5',
   all(is_complete_multipartite(5, g) == is_complete_multipartite_literal(5, g)
       for g in range(1 << 10)), 'all 1024 graphs on 5 vertices')

C5 = mask_of(5, [(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)])
ck('control_the_5_cycle_is_not_complete_multipartite', not is_complete_multipartite(5, C5))
ck('control_every_complete_graph_is_complete_multipartite',
   all(is_complete_multipartite(n, (1 << len(pair_list(n))) - 1) for n in (3, 4, 5, 6, 7)))
ck('control_two_disjoint_edges_are_not_but_one_edge_with_isolated_vertices_is',
   not is_complete_multipartite(4, mask_of(4, [(1, 2), (3, 4)]))
   and is_complete_multipartite(4, mask_of(4, [(1, 2)])),
   'K_2 + K_2 no; K_2 plus two vertices in no part yes')

# =====================================================================================
# STEP 3.  restriction, hence monotonicity of r under principal submatrices
# =====================================================================================

head('Step 3: restriction, hence monotonicity of r under principal submatrices')

viol = 0
for c in CM[6]:
    for size in range(1, 7):
        for S in itertools.combinations(range(1, 7), size):
            sub = mask_of(6, [e for e in edges_of(6, c) if e[0] in S and e[1] in S])
            if not is_complete_multipartite(6, sub, S):
                viol += 1
ck('an_induced_subgraph_of_a_complete_multipartite_graph_is_complete_multipartite',
   viol == 0, 'exhaustive at n=6: all %d such graphs x all 62 nonempty proper subsets, %d '
              'violations -- this is what makes r monotone' % (len(CM[6]), viol))

cm6, L6 = cover_layers(6, 3)
r6 = {}
for g in range(1 << 15):
    for k in range(4):
        if g in L6[k]:
            r6[g] = k
            break
mono = 0
for g in range(1 << 15):
    ge = edges_of(6, g)
    for i in range(1, 7):
        if r6[mask_of(6, [e for e in ge if i not in e])] > r6[g]:
            mono += 1
ck('r_never_increases_when_a_vertex_is_deleted', mono == 0,
   'exhaustive at n=6: 32768 graphs x 6 deletions, %d violations' % mono)

# =====================================================================================
# STEP 4.  the upper bound r(G) <= 3
# =====================================================================================

head('Step 4: the printed three-element cover, so r(G) <= 3')

cov3 = [parts_to_mask(7, p) for p in PAPER_COVER3]
for t, (parts, m) in enumerate(zip(PAPER_COVER3, cov3), start=1):
    ck('cover_member_%d_is_complete_multipartite_and_a_subgraph_of_G' % t,
       is_complete_multipartite(7, m) and not (m & ~G),
       'parts %s -> {%s}' % (parts, label(edges_of(7, m))))
ck('the_three_members_contribute_the_printed_edge_counts',
   [popcount(m) for m in cov3] == PAPER_COVER3_SIZES,
   ' + '.join(str(popcount(m)) for m in cov3))
ck('the_printed_edge_slot_total_is_the_sum_of_those_counts',
   sum(PAPER_COVER3_SIZES) == PAPER_COVER3_SLOTS,
   '12 + 6 + 1 = %d edge slots, which is not 15' % PAPER_COVER3_SLOTS)
dups = [e for e, c in Counter(e for m in cov3 for e in edges_of(7, m)).items() if c > 1]
ck('the_duplicated_edges_are_exactly_the_four_the_paper_names',
   sorted(dups) == sorted(parse_edge_words(PAPER_COVER3_DUPLICATES)), '{%s}' % label(dups))
ck('slots_minus_duplicates_is_the_printed_union_size',
   PAPER_COVER3_SLOTS - len(dups) == PAPER_COVER3_UNION == popcount(G),
   '%d - %d = %d = |E(G)|' % (PAPER_COVER3_SLOTS, len(dups), PAPER_COVER3_UNION))
ck('the_three_members_union_to_exactly_E_of_G', cov3[0] | cov3[1] | cov3[2] == G,
   '15 edges, and no edge outside G')

# =====================================================================================
# STEP 5.  the lower bound r(G) >= 3, exhaustively
# =====================================================================================

head('Step 5: r(G) >= 3 -- exhaustively, over the complete multipartite subgraphs of G')

sub7 = [c for c in CM[7] if c and not (c & ~G)]
ck('the_number_of_complete_multipartite_subgraphs_of_G_is_as_the_paper_states',
   len(sub7) == PAPER_CM_SUBGRAPHS_OF_G, '%d nonempty ones' % len(sub7))
ck('no_single_complete_multipartite_subgraph_of_G_covers_E_of_G', all(c != G for c in sub7),
   'so r(G) >= 2; %d candidates examined' % len(sub7))
b47, b56 = bit(7, 4, 7), bit(7, 5, 6)
ck('no_complete_multipartite_subgraph_of_G_contains_both_47_and_56',
   not any((c & b47) and (c & b56) for c in sub7),
   'step (i) of the paper, which is also what excludes r(G) <= 1')
covering_pairs = [(a, b) for a in sub7 for b in sub7 if (a | b) == G]
ck('no_pair_of_complete_multipartite_subgraphs_of_G_covers_E_of_G', not covering_pairs,
   '%d ordered pairs examined, 0 of them a cover; hence r(G) >= 3' % (len(sub7) ** 2))
ck('r_of_G_is_therefore_exactly_3',
   (not covering_pairs) and (cov3[0] | cov3[1] | cov3[2] == G),
   'lower bound from the pair search, upper bound from the printed cover')

# =====================================================================================
# STEP 6.  the seven principal 6x6 submatrices
# =====================================================================================

head('Step 6: each principal 6x6 submatrix, from the printed pair of covers')

for i in range(1, 8):
    m, V = Gdel[i]
    p1, p2 = PAPER_COVER2[i]
    c1, c2 = parts_to_mask(7, p1), parts_to_mask(7, p2)
    ck('printed_two_cover_is_valid_for_the_deletion_of_vertex_%d' % i,
       is_complete_multipartite(7, c1, V) and is_complete_multipartite(7, c2, V)
       and not (c1 & ~m) and not (c2 & ~m) and (c1 | c2) == m,
       '%s and %s, union E(G-%d) with %d edges' % (p1, p2, i, popcount(m)))

ck('no_principal_6x6_zero_entry_graph_is_complete_multipartite',
   all(not is_complete_multipartite(7, m, V) for m, V in Gdel.values()),
   'so every r(G-i) is exactly 2, never 1')
rvals = {i: cover_number(7, Gdel[i][0]) for i in range(1, 8)}
ck('an_independent_minimum_cover_search_returns_2_for_every_deletion',
   all(v == 2 for v in rvals.values()),
   ' '.join('r(G-%d)=%d' % (i, rvals[i]) for i in range(1, 8)))

# =====================================================================================
# STEP 7.  the source proposition, and the refutation
# =====================================================================================

head('Step 7: the source proposition applied, and the refutation')

cm7, L7 = cover_layers(7, 2)
ck('tree_rank_of_M_is_the_value_the_paper_claims',
   (G not in L7[2]) and (cov3[0] | cov3[1] | cov3[2] == G)
   and not isolated_vertices(7, G) and PAPER_TREE_RANK_M == 3,
   'no isolated vertex, so tree rank = r = 3')
ck('tree_rank_of_every_principal_6x6_submatrix_is_the_value_the_paper_claims',
   all((Gdel[i][0] in L7[2]) and (Gdel[i][0] not in L7[1])
       and not isolated_vertices(7, Gdel[i][0], Gdel[i][1]) for i in range(1, 8))
   and PAPER_TREE_RANK_PRINCIPAL == 2,
   'each is a union of two complete multipartite subgraphs and of no single one')
ck('REFUTATION_M_has_tree_rank_3_while_every_principal_6x6_has_tree_rank_at_most_2',
   (G not in L7[2]) and all(Gdel[i][0] in L7[2] for i in range(1, 8)),
   'the right-hand side of the conjecture holds for M and the left-hand side fails')

# =====================================================================================
# STEP 8.  controls quoted from the source paper
# =====================================================================================

head('Step 8: controls -- values published in the source paper, recomputed')

for n in (3, 4, 5, 6):
    _, Ln = cover_layers(n, 4 if n <= 5 else 3)
    trn = {g: tree_rank_small(n, g, Ln) for g in range(1 << len(pair_list(n)))}
    ck('maximum_tree_rank_over_the_whole_0_1_cell_at_n_%d_matches_the_source_table' % n,
       max(trn.values()) == SOURCE_MAX_TREE_RANK[n],
       'max %d, source says %d; distribution %s'
       % (max(trn.values()), SOURCE_MAX_TREE_RANK[n],
          dict(sorted(Counter(trn.values()).items()))))
    if n == 5:
        ck('control_the_5_cycle_has_the_tree_rank_the_source_publishes',
           trn[C5] == SOURCE_C5_TREE_RANK, 'tree rank(C_5) = %d' % trn[C5])
        three = [g for g in trn if trn[g] == 3]
        ck('the_number_of_labelled_tree_rank_3_objects_at_n_5_is_as_the_paper_states',
           len(three) == PAPER_N5_TREE_RANK_3, '%d labelled objects' % len(three))
        orb5, closed5 = orbits_under_relabelling(5, three)
        ck('they_form_one_isomorphism_class_and_it_is_the_5_cycle',
           len(orb5) == PAPER_N5_CLASSES and closed5 and C5 in set(orb5[0]),
           '%d class(es), orbit size %d, contains C_5' % (len(orb5), len(orb5[0])))

tr6_le2_masks = [g for g in range(1 << 15) if tree_rank_small(6, g, L6) <= 2]
ck('the_number_of_6x6_objects_of_tree_rank_at_most_2_is_as_the_paper_states',
   len(tr6_le2_masks) == PAPER_TR6_LE_2, '%d of 2^15 = 32768' % len(tr6_le2_masks))

# =====================================================================================
# STEP 9.  the uniform theorem: the same failure at every order
# =====================================================================================

head('Step 9: the uniform family, re-run for n = %d..%d' % (PAPER_UNIFORM_FROM, PAPER_UNIFORM_TO))

tr6_le2 = bytearray(1 << 15)
for g in tr6_le2_masks:
    tr6_le2[g] = 1
ONES_M = parse_edge_words(PAPER_NONEDGES)

for n in range(PAPER_UNIFORM_FROM, PAPER_UNIFORM_TO + 1):
    Gn_edges = [(i, j) for (i, j) in pair_list(n) if (i, j) not in ONES_M]
    Gn = mask_of(n, Gn_edges)
    ck('the_matrix_of_order_%d_has_the_same_six_one_entries_and_no_isolated_vertex' % n,
       sorted(set(pair_list(n)) - set(Gn_edges)) == sorted(ONES_M)
       and not isolated_vertices(n, Gn),
       '%d edges, %d one-entries, so tree rank = r' % (popcount(Gn), len(ONES_M)))
    three_parts = [[[1, 2, 3], [4, 5, 6, 7]] + [[u] for u in range(8, n + 1)],
                   [[1], [2], [4], [7]], [[5], [6]]]
    three = [parts_to_mask(n, p) for p in three_parts]
    ck('the_three_element_cover_of_the_uniform_theorem_is_valid_at_n_%d' % n,
       all(is_complete_multipartite(n, m) and not (m & ~Gn) for m in three)
       and three[0] | three[1] | three[2] == Gn, 'so r <= 3 and tree rank <= 3')
    induced = mask_of(7, [e for e in Gn_edges if e[0] <= 7 and e[1] <= 7])
    ck('the_first_seven_vertices_induce_exactly_G_at_n_%d' % n, induced == G,
       'with r(G) = 3 and monotonicity this forces r >= 3, hence tree rank exactly 3')
    ixn = pair_index(n)
    offender = None
    for S in itertools.combinations(range(1, n + 1), 6):
        rel = {v: k + 1 for k, v in enumerate(S)}
        h = mask_of(6, [(rel[a], rel[b]) for a, b in itertools.combinations(S, 2)
                        if (Gn >> ixn[(a, b)]) & 1])
        if not tr6_le2[h]:
            offender = S
            break
    ck('every_principal_6x6_submatrix_has_tree_rank_at_most_2_at_n_%d' % n,
       offender is None, 'all %d six-element subsets of [%d] checked' % (binom(n, 6), n))

G8 = mask_of(8, [(i, j) for (i, j) in pair_list(8) if (i, j) not in ONES_M])
sub8 = [c for c in CM[8] if c and not (c & ~G8)]
ck('the_number_of_complete_multipartite_subgraphs_at_n_8_is_as_the_paper_states',
   len(sub8) == PAPER_CM_SUBGRAPHS_AT_8, '%d nonempty ones' % len(sub8))
ck('an_exhaustive_pair_search_confirms_r_at_least_3_at_n_8_without_using_monotonicity',
   not any((a | b) == G8 for a in sub8 for b in sub8),
   '%d complete multipartite subgraphs, %d ordered pairs, 0 covers' % (len(sub8), len(sub8) ** 2))

# =====================================================================================
# STEP 10.  the exhaustive census of the whole 0/1 cell at n = 7
# =====================================================================================

head('Step 10: the exhaustive census of the 0/1 cell at n = 7')

N = 1 << 21
ck('the_cell_size_is_2_to_the_21', N == PAPER_CELL_SIZE, '2^21 = %d' % N)
r2set = set()
for a in cm7:
    r2set.update(a | b for b in cm7)
ck('the_number_of_graphs_with_r_at_most_1_is_as_the_paper_states',
   len(cm7) == PAPER_R_LE_1, '%d, the empty graph included' % len(cm7))
ck('the_number_of_graphs_with_r_at_most_2_is_as_the_paper_states',
   len(r2set) == PAPER_R_LE_2, '%d' % len(r2set))

ALL = int.from_bytes(b'\x01' * N, 'little')
R1 = bitset(N, cm7)
R2 = bitset(N, r2set)
inc7 = incidence(7)
ISO = []
for v in range(1, 8):
    acc = [0]
    for k in range(21):
        if not (inc7[v] >> k) & 1:
            acc += [x | (1 << k) for x in acc]
    ISO.append(bitset(N, acc))
iso_ge1 = 0
for x in ISO:
    iso_ge1 |= x
iso_ge2 = 0
for u in range(7):
    for v in range(u + 1, 7):
        iso_ge2 |= ISO[u] & ISO[v]
TR_LE2 = ((ALL ^ iso_ge2) & R2) | (iso_ge2 & R1)
ck('the_number_of_7x7_objects_of_tree_rank_at_most_2_is_as_the_paper_states',
   popcount(TR_LE2) == PAPER_TR7_LE_2, '%d' % popcount(TR_LE2))
ck('the_number_of_7x7_objects_of_tree_rank_at_least_3_is_as_the_paper_states',
   N - popcount(TR_LE2) == PAPER_TR7_GE_3, '%d' % (N - popcount(TR_LE2)))
ck('those_two_counts_exhaust_the_cell',
   popcount(TR_LE2) + (N - popcount(TR_LE2)) == PAPER_CELL_SIZE,
   '%d + %d = 2^21' % (popcount(TR_LE2), N - popcount(TR_LE2)))
ck('the_whole_cell_bitset_agrees_with_the_by_hand_verdict_on_M_and_its_seven_deletions',
   not ((TR_LE2 >> (8 * G)) & 1) and all((TR_LE2 >> (8 * Gdel[i][0])) & 1 for i in range(1, 8)),
   'M excluded, all seven principal 6x6 included')

ix7 = pair_index(7)
S_ALL = ALL
for d in range(1, 8):
    rest = [v for v in range(1, 8) if v != d]
    mp = [0] * 15
    for (a, b), k in pair_index(6).items():
        u, v = rest[a - 1], rest[b - 1]
        mp[k] = 1 << ix7[(min(u, v), max(u, v))]
    lo = [0] * 256
    for h in range(256):
        m = 0
        for k in range(8):
            if (h >> k) & 1:
                m |= mp[k]
        lo[h] = m
    hi = [0] * 128
    for h in range(128):
        m = 0
        for k in range(7):
            if (h >> k) & 1:
                m |= mp[8 + k]
        hi[h] = m
    fv = [0]
    for v in rest:
        f = 1 << ix7[(min(d, v), max(d, v))]
        fv += [x | f for x in fv]
    ba = bytearray(N)
    for h in tr6_le2_masks:
        base = lo[h & 255] | hi[h >> 8]
        for f in fv:
            ba[base | f] = 1
    S_ALL &= int.from_bytes(bytes(ba), 'little')

ck('the_number_of_objects_whose_seven_principal_6x6_all_have_tree_rank_at_most_2_is_as_stated',
   popcount(S_ALL) == PAPER_ALL_SEVEN_OK, '%d' % popcount(S_ALL))
CEX = S_ALL & (ALL ^ TR_LE2)
ck('the_number_of_counterexamples_in_the_cell_is_as_the_paper_states',
   popcount(CEX) == PAPER_COUNTEREXAMPLES, '%d' % popcount(CEX))
ck('the_arithmetic_of_the_two_census_totals_reproduces_that_count',
   PAPER_ALL_SEVEN_OK - PAPER_TR7_LE_2 == PAPER_COUNTEREXAMPLES,
   '%d - %d = %d' % (PAPER_ALL_SEVEN_OK, PAPER_TR7_LE_2, PAPER_COUNTEREXAMPLES))
ck('the_only_if_half_of_the_conjecture_is_silent_over_the_whole_cell',
   popcount(TR_LE2 & (ALL ^ S_ALL)) == PAPER_MONOTONICITY_VIOLATIONS,
   'objects of tree rank <= 2 with a principal 6x6 of tree rank >= 3: %d, and that half is a '
   'theorem, so it must be 0' % popcount(TR_LE2 & (ALL ^ S_ALL)))
ck('the_paper_witness_is_one_of_the_counterexamples', bool((CEX >> (8 * G)) & 1),
   'mask %d' % G)
ck('no_counterexample_in_the_cell_has_an_isolated_vertex', (CEX & iso_ge1) == 0,
   'so at n=7 the conjecture survives on every matrix with an all-ones off-diagonal row')

cex = [i for i, b in enumerate(CEX.to_bytes(N, 'little')) if b]
ck('the_extracted_list_of_counterexamples_has_the_counted_length',
   len(cex) == PAPER_COUNTEREXAMPLES, '%d masks' % len(cex))
note('one-entry counts among the counterexamples: %s'
     % dict(sorted(Counter(21 - popcount(g) for g in cex).items())))

# =====================================================================================
# STEP 11.  the isomorphism classes, and the control the discovery run never ran
# =====================================================================================

head('Step 11: isomorphism classes of the counterexamples')

orbs, closed = orbits_under_relabelling(7, cex)
ck('the_counterexample_set_is_closed_under_relabelling', closed,
   'every image of a counterexample under S_7 is one; %d masks' % len(cex))
ck('the_number_of_isomorphism_classes_of_counterexamples_is_as_the_paper_states',
   len(orbs) == PAPER_CLASSES, '%d classes' % len(orbs))
ck('the_orbit_sizes_are_as_the_paper_states',
   sorted(len(o) for o in orbs) == sorted(PAPER_ORBIT_SIZES),
   'sizes %s, summing to %d' % (sorted(len(o) for o in orbs), sum(len(o) for o in orbs)))
mine = [o for o in orbs if G in set(o)]
ck('the_paper_witness_lies_in_the_smaller_class',
   len(mine) == 1 and len(mine[0]) == min(PAPER_ORBIT_SIZES),
   'the orbit of M has size %d' % len(mine[0]))
other = [o for o in orbs if G not in set(o)][0]
rep = min(other)
ck('the_other_class_has_the_number_of_one_entries_the_paper_states',
   21 - popcount(rep) == PAPER_OTHER_CLASS_ONES,
   'representative zero-entry graph {%s}, %d one-entries'
   % (label(edges_of(7, rep)), 21 - popcount(rep)))
ck('the_other_class_is_the_class_of_G_with_the_edge_45_added',
   (G | bit(7, *PAPER_OTHER_CLASS_IS_G_PLUS)) in set(other),
   'G + %d%d is a counterexample too, so M does not have the fewest one-entries'
   % PAPER_OTHER_CLASS_IS_G_PLUS)

tot = 0
for perm in itertools.permutations(range(1, 8)):
    tab = [0] * 21
    for (a, b), k in ix7.items():
        u, v = perm[a - 1], perm[b - 1]
        tab[k] = 1 << ix7[(min(u, v), max(u, v))]
    nxt = [tab[k].bit_length() - 1 for k in range(21)]
    seenk, cyc = [False] * 21, 0
    for k in range(21):
        if not seenk[k]:
            cyc += 1
            j = k
            while not seenk[j]:
                seenk[j] = True
                j = nxt[j]
    tot += 1 << cyc
ck('control_Burnside_over_the_same_permutation_tables_returns_A000088_of_7',
   tot // 5040 == PAPER_A000088_7,
   'unlabelled graphs on 7 nodes: %d, and OEIS A000088(7) = %d' % (tot // 5040, PAPER_A000088_7))

# =====================================================================================
# scope, then the verdict
# =====================================================================================

head('Scope of this re-run')
note('NOT RE-RUN: the uniform theorem is proved in the paper for every n >= 7; this program '
     'confirms it only for n = 7,...,14, so n >= 15 rests on the hand proof alone.')
note('NOT RE-RUN: whether some local size k >= 7 does certify tree rank 2. The paper leaves '
     'that open and this family cannot settle it, because for k >= 7 the principal k x k '
     'submatrices of the uniform family include M itself.')
note('NOT RE-RUN: the Pachter--Sturmfels notion of tree rank named in the same bullet of the '
     'source. It is a different notion of mixtures, nothing here touches it, and no computation '
     'could -- their text was not read.')
note('NOT RE-RUN: the wording of the conjecture in the printed Combinatorica text. The statement '
     'is quoted from arXiv:0912.1411v1, and no computation can check a transcription.')
note('NOT RE-RUN: dissimilarity matrices that are not 0/1. The source proposition used here as '
     'the decider applies only to the 0/1 sub-family, and the census of Steps 10-11 is a census '
     'of that sub-family alone; a real counterexample of smaller order is not excluded.')
note('NOT RE-RUN: anything about the tropical Grassmannian directly. Every tree rank above is '
     'obtained through the source paper\'s own combinatorial proposition, which this program '
     'takes as given rather than reproving.')

print()
if _FAILED:
    print('VERDICT: %d of %d CHECKS FAILED' % (_FAILED, _PASSED + _FAILED))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _PASSED)
sys.exit(0)

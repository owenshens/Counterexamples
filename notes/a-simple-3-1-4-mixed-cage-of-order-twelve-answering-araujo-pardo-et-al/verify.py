#!/usr/bin/env python3
"""verify.py -- every computational claim of

    "Under a simplicity hypothesis, the order of a [3,1;4]-mixed cage is twelve"

re-derived from the object PRINTED IN THE PAPER, in exact integer arithmetic.

Python 3.9+, standard library only.  No third-party package, no external data file, no solver,
no randomness, no floating-point decision (the one non-integer constant in the paper, Shen's
2.885, is carried as fractions.Fraction(2885, 1000) and every comparison against it is exact).

WHAT IS AND IS NOT CHECKED HERE
  * The order-12 witness is TRANSCRIBED BY HAND into this file from the 36-arc / 6-edge link
    list printed in Section 2 of the paper, and independently rebuilt here from the circulant
    description; the two internal copies are checked label-equal. NOT CHECKED: that either copy
    agrees with what the paper prints -- this program never reads the paper.
  * The girth of every object is computed TWICE by independent routes: (i) exhaustive
    enumeration of closed walks on distinct vertices, which is the definition, and (ii) a
    breadth-first search over links that forbids reusing the same link twice in a row.
  * Lemma L (the whole lower bound) is confirmed by EXHAUSTIVE enumeration of all loopless,
    digon-free digraphs on m <= 5 vertices -- 3^C(m,2) of them, one state per unordered pair.
    That covers the only instance the paper uses (m = 3, delta = 1) and its tight case (m = 4).
  * The paper's exclusion of the orders 6, 8 and 10 at z = 3 is re-derived from that
    exhaustive Lemma L table, not from the printed proof.
  * NOT re-run here: nothing, of the note's own claims.  What this program cannot do is
    check the literature: the published values n[2,1;4] = 8, n[z,4] = 3z+1 and the Moore-type
    value n_AHM[1,1;4] = 6 are QUOTED constants, and the checks that use them are arithmetic
    consequences of those quotations, flagged as such in their check names.
"""
import sys
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product

# ---------------------------------------------------------------------------
# 0. THE CHECK LEDGER.  The verdict number is the counter, so it cannot disagree with the
#    PASS lines: there is exactly one place that prints either.
# ---------------------------------------------------------------------------
_N = [0]
_BAD = [0]


def check(name, ok, detail=''):
    if ok:
        _N[0] += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _BAD[0] += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))
    sys.stdout.flush()


def section(title):
    print('')
    print('=== %s ===' % title)


# ---------------------------------------------------------------------------
# 1. THE OBJECT, EXACTLY AS PRINTED IN SECTION 2 OF THE PAPER
#    (copied from the verbatim block; 'a u v' = arc u->v, 'e u v' = edge {u,v})
# ---------------------------------------------------------------------------
PRINTED_LINK_LIST = """
a 0 1   a 0 2   a 0 7   a 1 2   a 1 3   a 1 8   a 2 3   a 2 4   a 2 9
a 3 4   a 3 5   a 3 10  a 4 5   a 4 6   a 4 11  a 5 0   a 5 6   a 5 7
a 6 1   a 6 7   a 6 8   a 7 2   a 7 8   a 7 9   a 8 3   a 8 9   a 8 10
a 9 4   a 9 10  a 9 11  a 10 0  a 10 5  a 10 11 a 11 0  a 11 1  a 11 6
e 0 6   e 1 7   e 2 8   e 3 9   e 4 10  e 5 11
"""

# The same section's arithmetic description, used as the SECOND, independent route.
N_WITNESS = 12
ARC_DIFFS = (1, 2, 7)
EDGE_DIFF = 6
EXHIBITED_4_CYCLE = (0, 1, 3, 5)

# Labelling B (the canonical-matching form) and the relabelling phi; re-derived below.
OUT_B = {
    0: [2, 3, 4], 1: [2, 3, 5], 2: [4, 5, 6], 3: [4, 5, 7],
    4: [6, 7, 8], 5: [6, 7, 9], 6: [8, 9, 10], 7: [8, 9, 11],
    8: [1, 10, 11], 9: [0, 10, 11], 10: [0, 1, 3], 11: [0, 1, 2],
}


def parse_link_list(text):
    """-> (n, arcs set of (u,v), edges set of frozenset).  Plain text in, structure out."""
    toks = text.split()
    arcs, edges = set(), set()
    i = 0
    while i < len(toks):
        kind = toks[i]
        u, v = int(toks[i + 1]), int(toks[i + 2])
        i += 3
        if kind == 'a':
            arcs.add((u, v))
        elif kind == 'e':
            edges.add(frozenset((u, v)))
        else:
            raise ValueError('unknown link kind %r' % kind)
    verts = set()
    for (u, v) in arcs:
        verts.add(u)
        verts.add(v)
    for e in edges:
        verts |= set(e)
    return (max(verts) + 1), arcs, edges


def circulant(n, diffs, ediff):
    arcs = set((i, (i + d) % n) for i in range(n) for d in diffs)
    edges = set(frozenset((i, (i + ediff) % n)) for i in range(n))
    return arcs, edges


# ---------------------------------------------------------------------------
# 2. THE MIXED-GRAPH PRIMITIVES
# ---------------------------------------------------------------------------
def out_deg(n, arcs):
    d = [0] * n
    for (u, _v) in arcs:
        d[u] += 1
    return d


def in_deg(n, arcs):
    d = [0] * n
    for (_u, v) in arcs:
        d[v] += 1
    return d


def is_perfect_matching(n, edges):
    seen = [0] * n
    for e in edges:
        if len(e) != 2:
            return False
        for x in e:
            seen[x] += 1
    return all(c == 1 for c in seen)


def partner_map(n, edges):
    p = [None] * n
    for e in edges:
        a, b = sorted(e)
        p[a], p[b] = b, a
    return p


def linked(u, v, arcs, edges):
    """u -> v is traversable: an arc forwards, or an edge in either direction."""
    return (u, v) in arcs or frozenset((u, v)) in edges


def two_cycles(n, arcs, edges):
    """Cycles of length 2 need TWO DISTINCT links between the same pair: a digon, or an arc
    parallel to an edge.  A single edge walked out and back reuses one link and is not a
    cycle."""
    out = []
    for u in range(n):
        for v in range(u + 1, n):
            if (u, v) in arcs and (v, u) in arcs:
                out.append(('digon', u, v))
            if frozenset((u, v)) in edges and ((u, v) in arcs or (v, u) in arcs):
                out.append(('arc-parallel-to-edge', u, v))
    return out


def cycles_of_length(n, arcs, edges, L):
    """Every cycle on L >= 3 distinct vertices, up to rotation (first vertex smallest)."""
    found = []
    for tup in permutations(range(n), L):
        if tup[0] != min(tup):
            continue
        if all(linked(tup[i], tup[(i + 1) % L], arcs, edges) for i in range(L)):
            found.append(tup)
    return found


def girth_by_enumeration(n, arcs, edges, upto=5):
    """The definition: shortest closed walk on distinct vertices.  Returns (girth or None,
    a witness)."""
    tc = two_cycles(n, arcs, edges)
    if tc:
        return 2, tc[0]
    for L in range(3, upto + 1):
        cs = cycles_of_length(n, arcs, edges, L)
        if cs:
            return L, cs[0]
    return None, None


def girth_by_bfs(n, arcs, edges):
    """Independent route: BFS over LINKS, forbidding reuse of the same link twice in a row.
    A shortest such closed walk is a cycle, so this must agree with girth_by_enumeration."""
    out = [[] for _ in range(n)]
    for (u, v) in arcs:
        out[u].append(v)
    p = partner_map(n, edges)

    def succ(v):
        s = [(w, ('a', v, w)) for w in out[v]]
        if p[v] is not None:
            s.append((p[v], ('e', frozenset((v, p[v])))))
        return s

    best = None
    for s in range(n):
        frontier = [(w, lk, 1) for (w, lk) in succ(s)]
        seen = {}
        while frontier:
            nxt = []
            for (v, lk, d) in frontier:
                if best is not None and d + 1 >= best:
                    continue
                for (w, nlk) in succ(v):
                    if nlk == lk:
                        continue
                    if w == s:
                        if best is None or d + 1 < best:
                            best = d + 1
                        continue
                    if seen.get((w, nlk), 10 ** 9) <= d + 1:
                        continue
                    seen[(w, nlk)] = d + 1
                    nxt.append((w, nlk, d + 1))
            frontier = nxt
    return best


def validate_mixed(n, arcs, edges, z, tag):
    """Full validation of a candidate [z,1;4]-mixed graph.  Returns list of problems."""
    bad = []
    if any((u, u) in arcs for u in range(n)):
        bad.append('loop')
    od, idg = out_deg(n, arcs), in_deg(n, arcs)
    if any(x != z for x in od):
        bad.append('out-degrees %s' % od)
    if any(x != z for x in idg):
        bad.append('in-degrees %s' % idg)
    if not is_perfect_matching(n, edges):
        bad.append('edges are not a perfect matching')
    g, _w = girth_by_enumeration(n, arcs, edges)
    if g != 4:
        bad.append('girth %s' % g)
    gb = girth_by_bfs(n, arcs, edges)
    if gb != g:
        bad.append('the two girth routes disagree: %s vs %s' % (g, gb))
    return bad, g, gb, tag


# ---------------------------------------------------------------------------
# 3. LEMMA L, BY EXHAUSTIVE ENUMERATION
# ---------------------------------------------------------------------------
def lemma_L_table(mmax=5):
    """For each m <= mmax: the largest minimum-out-degree delta attainable by a LOOPLESS,
    DIGON-FREE, DIRECTED-TRIANGLE-FREE digraph on m vertices (0 if none has delta >= 1), and
    whether every such digraph satisfies m >= 2*delta + 2.

    Enumeration is exact and complete: a loopless digon-free digraph on m labelled vertices is
    exactly a choice, for each unordered pair, of {no arc, i->j, j->i}, so there are
    3^C(m,2) of them and all are visited."""
    table = {}
    for m in range(1, mmax + 1):
        pairs = list(combinations(range(m), 2))
        best_delta, holds, witness = 0, True, None
        for state in product((0, 1, 2), repeat=len(pairs)):
            arcs = set()
            for (idx, (i, j)) in enumerate(pairs):
                if state[idx] == 1:
                    arcs.add((i, j))
                elif state[idx] == 2:
                    arcs.add((j, i))
            od = [0] * m
            for (u, _v) in arcs:
                od[u] += 1
            delta = min(od) if m else 0
            if delta < 1:
                continue
            tri = False
            for (i, j, k) in permutations(range(m), 3):
                if i != min((i, j, k)):
                    continue
                if (i, j) in arcs and (j, k) in arcs and (k, i) in arcs:
                    tri = True
                    break
            if tri:
                continue
            if m < 2 * delta + 2:
                holds = False
            if delta > best_delta:
                best_delta, witness = delta, sorted(arcs)
        table[m] = {'max_delta': best_delta, 'bound_holds': holds, 'witness': witness}
    return table


# ---------------------------------------------------------------------------
# 4. THE PAPER'S LOWER-BOUND ARGUMENT, EXECUTED
# ---------------------------------------------------------------------------
def order_excluded(n, z, L):
    """Re-derive the paper's exclusion of order n for a [z,1;4]-mixed graph, using ONLY the
    exhaustive Lemma L table L.  Returns (excluded, reason)."""
    if n % 2:
        return True, 'order is odd and r = 1 forces a perfect matching'
    if n < 2 * z + 2:
        return True, 'the disjoint sets {u},{v},N+(u),N-(u) already need %d vertices' % (2 * z + 2)
    S = n - 2 * z - 2
    d = z - S
    if d <= 0:
        return False, 'S absorbs every out-arc of A; this branch gives only n >= 3z+2'
    if z not in L:
        return False, 'Lemma L table has no entry for m = %d' % z
    if L[z]['max_delta'] < d:
        return True, ('the digraph induced on N+(u) has m = %d, is loopless, digon-free and '
                      'directed-triangle-free with min out-degree >= %d, and the exhaustive '
                      'Lemma L table caps that at %d' % (z, d, L[z]['max_delta']))
    return False, 'Lemma L permits m = %d with min out-degree %d' % (z, d)


def ceil_div(a, b):
    return -((-a) // b)


def thm3_raw(z):
    """ceil((5z+6)/2), for z >= 2."""
    return ceil_div(5 * z + 6, 2)


def thm3_even(z):
    """the same lifted to an even order: 2*ceil((5z+6)/4)."""
    return 2 * ceil_div(5 * z + 6, 4)


def source_upper(z):
    """Theorem 3.3 of the source: 3(z+1) for odd z, 3z+2 for even z."""
    return 3 * (z + 1) if z % 2 else 3 * z + 2


def even_up(x):
    """least even integer >= x, for an integer x."""
    return x if x % 2 == 0 else x + 1


# ---------------------------------------------------------------------------
# 5. RUN
# ---------------------------------------------------------------------------
def main():
    print('verification of "Under a simplicity hypothesis, the order of a '
          '[3,1;4]-mixed cage is twelve"')
    print('python %s ; exact integer / Fraction arithmetic only, no floating-point decision'
          % sys.version.split()[0])

    # ---- 5.1 the printed object -------------------------------------------------
    section('part 1: the order-12 witness, transcribed by hand from Section 2 and rebuilt from the circulant')
    n, arcs, edges = parse_link_list(PRINTED_LINK_LIST)
    check('printed_list_has_12_vertices_36_arcs_6_edges',
          (n, len(arcs), len(edges)) == (12, 36, 6),
          'n=%d, |A|=%d, |E|=%d' % (n, len(arcs), len(edges)))

    carcs, cedges = circulant(N_WITNESS, ARC_DIFFS, EDGE_DIFF)
    check('printed_arcs_equal_the_circulant_with_differences_1_2_7', arcs == carcs,
          'two independent constructions, label-equal as sets of ordered pairs')
    check('printed_edges_equal_the_difference_6_perfect_matching', edges == cedges,
          '{i, i+6} for i = 0..5')

    check('no_loop', not any((u, u) in arcs for u in range(n)))
    od, idg = out_deg(n, arcs), in_deg(n, arcs)
    check('out_degree_is_3_at_every_vertex', all(x == 3 for x in od), 'min=max=%d' % od[0])
    check('in_degree_is_3_at_every_vertex', all(x == 3 for x in idg), 'min=max=%d' % idg[0])
    check('edges_form_a_perfect_matching_so_r_equals_1', is_perfect_matching(n, edges))

    in_diffs = sorted((N_WITNESS - d) % N_WITNESS for d in ARC_DIFFS)
    check('paper_in_differences_are_5_10_11', in_diffs == [5, 10, 11], str(in_diffs))
    check('no_in_difference_lies_in_the_arc_set_so_no_digon',
          not (set(in_diffs) & set(ARC_DIFFS)),
          '{5,10,11} cap {1,2,7} = empty')
    check('edge_difference_6_is_not_an_arc_difference_so_no_arc_parallel_to_an_edge',
          EDGE_DIFF not in ARC_DIFFS)

    sums3 = [sum(c) % N_WITNESS for c in combinations_with_replacement(ARC_DIFFS, 3)]
    check('the_ten_three_term_sums_are_exactly_those_printed',
          sums3 == [3, 4, 9, 5, 10, 3, 6, 11, 4, 9], str(sums3))
    check('no_three_term_sum_vanishes_so_no_directed_triangle', 0 not in sums3)
    sums2 = [sum(c) % N_WITNESS for c in combinations_with_replacement(ARC_DIFFS, 2)]
    check('the_six_two_term_sums_are_exactly_those_printed',
          sums2 == [2, 3, 8, 4, 9, 2], str(sums2))
    check('no_two_term_sum_equals_the_edge_difference_so_no_mixed_triangle',
          EDGE_DIFF not in sums2)

    check('there_is_no_cycle_of_length_2', two_cycles(n, arcs, edges) == [])
    check('there_is_no_cycle_of_length_3_exhaustively',
          cycles_of_length(n, arcs, edges, 3) == [],
          'all %d ordered triples of distinct vertices examined' % (12 * 11 * 10))

    c4 = EXHIBITED_4_CYCLE
    check('the_exhibited_4_cycle_0_1_3_5_is_a_directed_4_cycle',
          all((c4[i], c4[(i + 1) % 4]) in arcs for i in range(4)),
          '1+2+2+7 = 12 = 0 mod 12')
    g_enum, wit = girth_by_enumeration(n, arcs, edges)
    check('girth_is_exactly_4_by_definitional_enumeration', g_enum == 4,
          'shortest cycle %s' % (wit,))
    g_bfs = girth_by_bfs(n, arcs, edges)
    check('girth_is_exactly_4_by_an_independent_link_bfs', g_bfs == 4, 'bfs girth = %s' % g_bfs)
    check('the_two_girth_routes_agree', g_enum == g_bfs)
    check('the_witness_is_a_valid_3_1_4_mixed_graph_of_order_12',
          validate_mixed(n, arcs, edges, 3, 'C_12(1,2,7)')[0] == [],
          'degrees, matching and girth all as claimed')

    # ---- 5.2 the two labellings -------------------------------------------------
    section('part 2: the second labelling and the relabelling phi')
    phi = {}
    for i in range(6):
        phi[i] = 2 * i
        phi[i + 6] = 2 * i + 1
    check('phi_is_a_bijection_of_Z12', sorted(phi.values()) == list(range(12)))
    arcs_B = set((u, v) for u in OUT_B for v in OUT_B[u])
    edges_B = set(frozenset((2 * i, 2 * i + 1)) for i in range(6))
    check('labelling_B_has_36_arcs_and_6_edges',
          (len(arcs_B), len(edges_B)) == (36, 6),
          '|A|=%d, |E|=%d' % (len(arcs_B), len(edges_B)))
    check('phi_maps_the_arcs_of_labelling_A_onto_the_arcs_of_labelling_B',
          set((phi[u], phi[v]) for (u, v) in arcs) == arcs_B)
    check('phi_maps_the_edges_of_labelling_A_onto_the_canonical_matching',
          set(frozenset((phi[u], phi[v])) for (u, v) in (tuple(e) for e in edges)) == edges_B,
          '{{0,1},{2,3},{4,5},{6,7},{8,9},{10,11}}')
    bad_B, gB, gbB, _ = validate_mixed(12, arcs_B, edges_B, 3, 'labelling B')
    check('labelling_B_is_independently_a_valid_3_1_4_mixed_graph_of_girth_4',
          bad_B == [], 'girth %s by both routes' % gB)

    # ---- 5.3 Lemma L ------------------------------------------------------------
    section('part 3: Lemma L, by exhaustive enumeration of small digraphs')
    L = lemma_L_table(5)
    for m in range(1, 6):
        npairs = m * (m - 1) // 2
        check('lemma_L_holds_for_every_loopless_digon_free_triangle_free_digraph_on_%d_vertices' % m,
              L[m]['bound_holds'],
              'all 3^%d = %d such digraphs enumerated; max attainable min-out-degree = %d'
              % (npairs, 3 ** npairs, L[m]['max_delta']))
    check('no_such_digraph_on_3_vertices_has_min_out_degree_at_least_1',
          L[3]['max_delta'] == 0,
          'this single fact is the whole lower bound at z = 3')
    check('lemma_L_is_tight_at_delta_1_on_4_vertices',
          L[4]['max_delta'] == 1,
          'witness %s, and 4 = 2*1+2' % (L[4]['witness'],))
    check('on_5_vertices_the_min_out_degree_still_cannot_exceed_1',
          L[5]['max_delta'] == 1, 'consistent with 5 >= 2*delta+2 forcing delta <= 1')

    # ---- 5.4 the lower bound at the target -------------------------------------
    section('part 4: the exclusion of every order below 12 at z = 3')
    for n0 in (2, 4, 6, 8, 10):
        ex, why = order_excluded(n0, 3, L)
        check('order_%d_is_impossible_for_a_3_1_4_mixed_graph' % n0, ex, why)
    ex12, why12 = order_excluded(12, 3, L)
    check('order_12_is_NOT_excluded_by_the_argument', not ex12,
          'the argument must not refute its own witness; %s' % why12)
    check('orders_2_through_11_are_excluded_and_girth_of_the_order_12_object_is_4',
          all(order_excluded(k, 3, L)[0] for k in range(2, 12)) and g_enum == 4,
          'each order k with 2 <= k <= 11 excluded, and the order-12 object has girth 4')

    section('part 5: the general bound of Theorem 3')
    check('theorem_3_at_z_3_gives_11_and_parity_lifts_it_to_12',
          (thm3_raw(3), thm3_even(3)) == (11, 12),
          'ceil(21/2) = 11, 2*ceil(21/4) = 12')
    check('theorem_3_at_z_2_reproduces_the_published_exact_value_8',
          (thm3_raw(2), thm3_even(2)) == (8, 8),
          'forced positive: n[2,1;4] = 8 is proved in the source')
    check('theorem_3_never_exceeds_the_source_upper_bound_for_2_le_z_le_200',
          all(thm3_even(z) <= source_upper(z) for z in range(2, 201)),
          'a lower bound above a construction would be a bug')
    check('the_source_upper_bound_is_even_for_every_z',
          all(source_upper(z) % 2 == 0 for z in range(1, 201)),
          'as it must be, since r = 1 forces even order')
    check('transcribed_brackets_at_z_5_and_z_6_are_16_to_18_and_18_to_20',
          (thm3_even(5), source_upper(5), thm3_even(6), source_upper(6)) == (16, 18, 18, 20),
          'the two transcribed formulas of [APdM] Theorem 3 evaluate to the brackets 16..18 at z = 5 and \n'
          '18..20 at z = 6. IMPORTED, NOT PROVED HERE: that those formulas are correct bounds on '
          'n[5,1;4] and n[6,1;4]; no check below establishes either endpoint')

    # ---- 5.5 controls on the other published witnesses -------------------------
    section('part 6: the source\'s other constructions, as forced positives')
    for (z, nn, diffs, ed, label) in ((1, 6, (1,), 3, 'C_6(1): n[1,1;4] <= 6'),
                                      (2, 8, (1, 5), 4, 'C_8(1,5): n[2,1;4] <= 8'),
                                      (4, 14, (1, 2, 8, 9), 7, 'C_14(1,2,8,9): n[4,1;4] <= 14')):
        a, e = circulant(nn, diffs, ed)
        bad, gg, ggb, _ = validate_mixed(nn, a, e, z, label)
        check('published_witness_z%d_n%d_is_a_valid_mixed_graph_of_girth_4' % (z, nn),
              bad == [], '%s ; girth %s by enumeration and %s by bfs' % (label, gg, ggb))

    # ---- 5.6 the scope limit: the multigraph reading ---------------------------
    section('part 7: the multigraph reading, which is a different and undecided problem')
    # arc i -> i+1 with multiplicity 3, edges {i, i+3}, on Z_6.
    M_N, M_DIFF, M_MULT, M_EDIFF = 6, 1, 3, 3
    m_arcs, m_edges = circulant(M_N, (M_DIFF,), M_EDIFF)      # the SUPPORT
    m_multiset = [(i, (i + M_DIFF) % M_N) for i in range(M_N) for _ in range(M_MULT)]
    m_od, m_id = [0] * M_N, [0] * M_N
    for (u, v) in m_multiset:
        m_od[u] += 1
        m_id[v] += 1
    check('multigraph_object_has_out_and_in_degree_3_with_multiplicity',
          len(m_multiset) == 18 and all(x == 3 for x in m_od) and all(x == 3 for x in m_id),
          'counted from the explicit 18-arc multiset: out %s, in %s' % (m_od, m_id))
    check('multigraph_object_edges_form_a_perfect_matching_of_Z6',
          is_perfect_matching(M_N, m_edges))
    check('parallel_arcs_create_no_cycle_of_length_2',
          (M_N - M_DIFF) % M_N not in (M_DIFF,),
          'parallel arcs point the same way; and -1 = 5 is not an arc difference')
    check('multigraph_object_has_no_arc_parallel_to_an_edge', M_EDIFF != M_DIFF)
    check('multigraph_object_has_no_directed_triangle', (3 * M_DIFF) % M_N != 0, '1+1+1 = 3 != 0 mod 6')
    check('multigraph_object_has_no_mixed_triangle', (2 * M_DIFF) % M_N != M_EDIFF, '1+1 = 2 != 3')
    check('multigraph_object_has_a_cycle_of_length_4', (3 * M_DIFF) % M_N == M_EDIFF,
          '0->1->2->3 closed by the edge {3,0}')
    g_supp, _ = girth_by_enumeration(M_N, m_arcs, m_edges)
    check('so_under_a_multigraph_reading_n_3_1_4_is_at_most_6', g_supp == 4,
          'girth of the support is %s, and multiplicity changes no cycle length' % g_supp)
    bad_simple, _, _, _ = validate_mixed(M_N, m_arcs, m_edges, 3, 'multi-arc Z_6')
    check('the_same_object_is_REJECTED_as_a_simple_z3_object', bad_simple != [],
          'anti-control; problems: %s' % '; '.join(bad_simple))

    # ---- 5.7 the near-miss preprint, arithmetic on quoted constants ------------
    section('part 8: arithmetic on QUOTED literature constants (not derived here)')
    N_AHM_1_1_4 = 2 * (1 + 2)                      # source, line 173: n_0[1,r;4] = 2(r+2)
    check('quoted_n_AHM_1_1_4_equals_6', N_AHM_1_1_4 == 6, 'n_0[1,r;4] = 2(r+2) at r = 1')

    def dircage(z):
        return 3 * z + 1                            # n[z,4] as used in arXiv:2401.14768

    def preprint(z):
        return dircage(z) + N_AHM_1_1_4 - 4

    check('the_preprint_inequality_substitutes_to_exactly_12_at_z_3',
          preprint(3) == 12, '10 + 6 - 4 = 12, the same number this paper proves')
    check('the_preprint_inequality_is_false_at_z_2_given_the_published_exact_value_8',
          preprint(2) == 9 and 9 > 8, '7 + 6 - 4 = 9 > 8 = n[2,1;4]')
    check('the_preprint_inequality_is_false_at_z_4_given_the_verified_order_14_construction',
          preprint(4) == 15 and 15 > 14, '13 + 6 - 4 = 15 > 14, and C_14(1,2,8,9) is verified above')

    section('part 9: the baseline comparison printed in Section 4')
    SHEN = Fraction(2885, 1000)
    def shen_even(z):
        return even_up(ceil_div((SHEN * z).numerator, (SHEN * z).denominator))
    check('shen_baseline_and_theorem_3_at_z_3_are_10_and_12',
          (shen_even(3), thm3_even(3)) == (10, 12),
          'the strict improvement is at the target cell')
    check('the_directed_cage_baseline_at_z_3_is_also_only_10',
          even_up(dircage(3)) == 10, 'n[3,4] = 10, already even')
    check('the_directed_cage_baseline_at_z_4_ties_theorem_3_at_14',
          even_up(dircage(4)) == 14 and thm3_even(4) == 14, 'n[4,4] = 13, lifted to 14')
    got = [shen_even(z) for z in range(7, 13)]
    mine = [thm3_even(z) for z in range(7, 13)]
    check('the_z_7_to_12_comparison_table_is_exactly_as_printed',
          got == [22, 24, 26, 30, 32, 36] and mine == [22, 24, 26, 28, 32, 34],
          'shen %s vs theorem 3 %s -- ties at 7,8,9,11 and losses at 10,12' % (got, mine))

    # ---- verdict ---------------------------------------------------------------
    print('')
    print('NOTE SCOPE. This program checks the note\'s own claims and nothing else. It reads the '
          '36-arc/6-edge link list transcribed by hand from Section 2, rebuilds it independently '
          'from the circulant description, and re-derives every difference, sum, girth and '
          'bound; Lemma L is confirmed by complete enumeration of all loopless digon-free '
          'digraphs on at most 5 vertices, which covers the m = 3 instance the lower bound uses. '
          'GAPS NOT COVERED: this program never reads the paper, so the agreement of the '
          'transcribed link list with the printed one is not checked here; and the published '
          'constants n[2,1;4] = 8, n[z,4] = 3z+1 and '
          'n_AHM[1,1;4] = 6 are QUOTED, not derived here, so part 8 and part 9 are arithmetic '
          'consequences of the literature rather than independent verifications of it; nothing '
          'here decides the multigraph variant beyond exhibiting an order-6 object for it; '
          'nothing here claims the order-12 witness is unique; and nothing here touches girth '
          'g != 4 or edge-degree r != 1.')
    print('')
    if _BAD[0]:
        print('VERDICT: %d OF %d CHECKS FAILED' % (_BAD[0], _BAD[0] + _N[0]))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _N[0])
    return 0


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""Verification program for

    "A union of two non-DDMO copies of K_{3,3} is difference distance magic orientable"

It reads the objects printed in the note -- the oriented K_{3,3} of Section 2 and its two
labelings -- and re-derives the quantities the note asserts.  It also carries auxiliary
checks on order-5 oriented graphs (W_4 and relatives) that the note does not claim.

Python 3.9+, standard library only.  Exact integer / Fraction arithmetic throughout; no
floating point anywhere, no external data file, no third-party package.

Convention (the one displayed in the note):
    wt(v) = (sum of f over the IN-neighbours of v) - (sum of f over the OUT-neighbours of v),
and an arc written (u, v) means u -> v.  f is a DDM labeling of an n-vertex oriented graph iff
f is a bijection onto {1,...,n} and wt(v) = 0 for every vertex v.

One `PASS <name>` line per check; the run closes with
    VERDICT: ALL <n> CHECKS PASS
and exits 0 iff every check passed.
"""

import itertools
import sys
from fractions import Fraction
from math import gcd

# ----------------------------------------------------------------------------------------
# the objects, exactly as printed in the paper
# ----------------------------------------------------------------------------------------

# Section 2: the orientation D of K_{3,3}, parts A = {0,1,2}, B = {3,4,5}.
D_ARCS = [(0, 3), (4, 0), (5, 0),
          (3, 1), (1, 4), (1, 5),
          (3, 2), (2, 4), (2, 5)]
A_PART = (0, 1, 2)
B_PART = (3, 4, 5)

# Section 2: the two labelings of the n = 12 witness  D  u  D.
COPY1 = (6, 1, 5, 10, 2, 8)
COPY2 = (12, 3, 9, 11, 4, 7)

# An orientation W of W_4 = K_1 v C_4, hub 0, and a labeling of it (not printed in the note).
W_ARCS = [(1, 0), (0, 2), (0, 3), (0, 4),
          (3, 1), (4, 1), (2, 3), (2, 4)]
W_LABELS = (10, 11, 1, 3, 7)
# A labeling of the K_{3,3} component of the n = 11 object (not printed in the note).
COPY_A = (8, 2, 6, 9, 4, 5)

# An orientation of W_4 that IS a DDMOG (hub 0 carries label 1); not printed in the note.
W4DDM_ARCS = [(1, 0), (4, 0), (0, 3), (0, 2),
              (4, 1), (1, 2), (2, 3), (3, 4)]
W4DDM_LABELS = (1, 2, 4, 3, 5)

# ----------------------------------------------------------------------------------------
# bookkeeping
# ----------------------------------------------------------------------------------------

PASSED = []
FAILED = []


def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    (PASSED if cond else FAILED).append(name)
    print('%s %s%s' % (tag, name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE ' + text)


# ----------------------------------------------------------------------------------------
# primitives
# ----------------------------------------------------------------------------------------

def weights(n, arcs, f):
    """The vector (wt(0), ..., wt(n-1)) of an oriented graph under the labeling f."""
    w = [0] * n
    for (u, v) in arcs:
        w[v] += f[u]      # u is an IN-neighbour of v
        w[u] -= f[v]      # v is an OUT-neighbour of u
    return w


def wmatrix(n, arcs):
    """The integer matrix M with (wt(0),...,wt(n-1)) = M f.  M is skew-symmetric."""
    M = [[0] * n for _ in range(n)]
    for (u, v) in arcs:
        M[v][u] += 1
        M[u][v] -= 1
    return M


def edges_of(arcs):
    return set(frozenset(a) for a in arcs)


def is_orientation(arcs):
    """True iff no arc is repeated and no pair of opposite arcs occurs."""
    s = set(arcs)
    if len(s) != len(arcs):
        return False
    return all((v, u) not in s for (u, v) in arcs)


def degrees(n, edges):
    d = [0] * n
    for e in edges:
        for x in e:
            d[x] += 1
    return d


def all_orientations(edges):
    """Every one of the 2^|E| orientations of the given undirected edge set."""
    es = [tuple(sorted(e)) for e in sorted(edges, key=lambda e: tuple(sorted(e)))]
    for mask in range(1 << len(es)):
        yield [(u, v) if (mask >> i) & 1 else (v, u) for i, (u, v) in enumerate(es)]


def has_ddm(n, arcs):
    """Brute force over all n! bijections onto {1,...,n}: is this oriented graph a DDMOG?"""
    return ddm_labelings(n, arcs, first_only=True) != []


def ddm_labelings(n, arcs, first_only=False):
    """Every DDM labeling (bijection onto {1..n} with all weights 0), by brute force."""
    inn = [[] for _ in range(n)]
    out = [[] for _ in range(n)]
    for (u, v) in arcs:
        inn[v].append(u)
        out[u].append(v)
    found = []
    for p in itertools.permutations(range(1, n + 1)):
        ok = True
        for v in range(n):
            s = 0
            for u in inn[v]:
                s += p[u]
            for u in out[v]:
                s -= p[u]
            if s:
                ok = False
                break
        if ok:
            found.append(p)
            if first_only:
                return found
    return found


def kernel_param(M, n):
    """Reduced row echelon form of M over Q.

    Returns (pivots, free, rows) where `rows` gives, for each pivot column c,
        f_c = (sum_j ints[j] * f_{free[j]}) / den
    with ints and den integers.  Every solution of M f = 0 is obtained by choosing the free
    coordinates arbitrarily, so enumerating the free coordinates enumerates all solutions.
    """
    A = [[Fraction(x) for x in row] for row in M]
    pivots = []
    r = 0
    for c in range(n):
        p = None
        for i in range(r, n):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                fac = A[i][c]
                A[i] = [a - fac * b for a, b in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
        if r == n:
            break
    free = [c for c in range(n) if c not in pivots]
    rows = []
    for k, c in enumerate(pivots):
        coeffs = [-A[k][j] for j in free]
        den = 1
        for x in coeffs:
            den = den * x.denominator // gcd(den, x.denominator)
        rows.append((c, [int(x * den) for x in coeffs], den))
    return pivots, free, rows


def realisable(n, arcs, lo, hi, want_labelings=False):
    """All label sets in {lo,...,hi} realised by this orientation, via the exact kernel.

    A labeling f with all weights zero is exactly a solution of M f = 0, so the free
    coordinates of the kernel determine f.  Each free coordinate is itself a label, hence lies
    in {lo,...,hi}; enumerating injective assignments to the free coordinates is therefore
    complete.  Returns (set of frozensets of realised label sets, list of labelings).
    """
    pivots, free, rows = kernel_param(wmatrix(n, arcs), n)
    sets = set()
    sols = []
    if not free:                      # kernel is {0}: no positive labeling at all
        return sets, sols
    values = range(lo, hi + 1)
    for vals in itertools.permutations(values, len(free)):
        f = [0] * n
        for j, c in enumerate(free):
            f[c] = vals[j]
        ok = True
        for (c, ints, den) in rows:
            s = 0
            for j, iv in enumerate(ints):
                s += iv * vals[j]
            if s % den:
                ok = False
                break
            x = s // den
            if x < lo or x > hi:
                ok = False
                break
            f[c] = x
        if not ok:
            continue
        if len(set(f)) != n:
            continue
        sets.add(frozenset(f))
        if want_labelings:
            sols.append(tuple(f))
    return sets, sols


def canon(n, edges):
    """Canonical form of a graph on n labelled vertices: min over all n! relabelings."""
    best = None
    for p in itertools.permutations(range(n)):
        t = tuple(sorted(tuple(sorted((p[u], p[v]))) for (u, v) in (tuple(e) for e in edges)))
        if best is None or t < best:
            best = t
    return best


ALL_PAIRS_5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
E_K5 = set(frozenset(e) for e in ALL_PAIRS_5)
E_K5E = set(e for e in E_K5 if e != frozenset((3, 4)))
E_W4 = set(e for e in E_K5 if e not in (frozenset((1, 2)), frozenset((3, 4))))
E_K33 = set(frozenset((a, b)) for a in A_PART for b in B_PART)
E_K4 = set(frozenset(e) for e in [(i, j) for i in range(4) for j in range(i + 1, 4)])

# ========================================================================================
print('verification of the objects printed in: A Union of Two Non-DDMO Copies of K_{3,3}')
print('Is Difference Distance Magic Orientable')
print('python %s, exact integer and Fraction arithmetic only' % sys.version.split()[0])
print('')

# ----------------------------------------------------------------------------------------
print('=== Part 1: the oriented K_{3,3} of the note, and an oriented W_4 used below')
# ----------------------------------------------------------------------------------------

check('D_is_an_orientation_no_repeated_and_no_opposite_arcs', is_orientation(D_ARCS),
      '%d arcs' % len(D_ARCS))
check('D_underlying_graph_is_K33_on_parts_012_and_345',
      edges_of(D_ARCS) == E_K33 and len(edges_of(D_ARCS)) == 9,
      'all 9 A-B pairs occur exactly once')
check('D_is_3_regular', degrees(6, edges_of(D_ARCS)) == [3, 3, 3, 3, 3, 3])

din = [sum(1 for (u, v) in D_ARCS if v == x) for x in range(6)]
dout = [sum(1 for (u, v) in D_ARCS if u == x) for x in range(6)]
check('D_in_and_out_degrees_are_as_printed', din == [2, 1, 1, 1, 2, 2] and dout == [1, 2, 2, 2, 1, 1],
      'in %s out %s' % (tuple(din), tuple(dout)))

# The collapse of the six vertex equations to two, checked exactly on the weight matrix.
E1 = [1, -1, -1, 0, 0, 0]          # f0 - f1 - f2
E2 = [0, 0, 0, 1, -1, -1]          # f3 - f4 - f5
MD = wmatrix(6, D_ARCS)
neg = lambda r: [-x for x in r]
rowclass = [(r == E1, r == neg(E1), r == E2, r == neg(E2)) for r in MD]
check('D_weight_system_collapses_to_exactly_two_equations',
      all(sum(rc) == 1 for rc in rowclass)
      and [MD[v] for v in range(6)] == [neg(E2), E2, E2, E1, neg(E1), neg(E1)],
      'every row of M is +-(f0-f1-f2) or +-(f3-f4-f5), so wt=0 iff f0=f1+f2 and f3=f4+f5')

check('W_is_an_orientation_no_repeated_and_no_opposite_arcs', is_orientation(W_ARCS),
      '%d arcs' % len(W_ARCS))
check('W_underlying_graph_is_W4_equals_K5_minus_two_independent_edges',
      edges_of(W_ARCS) == E_W4 and (E_K5 - edges_of(W_ARCS)) == {frozenset((1, 2)), frozenset((3, 4))},
      'complement inside K_5 is {1,2},{3,4}')
check('W_degree_sequence_is_4_3_3_3_3', sorted(degrees(5, edges_of(W_ARCS)), reverse=True) == [4, 3, 3, 3, 3],
      'hub 0 has degree 4; rim is the 4-cycle 1-3-2-4-1')

# ----------------------------------------------------------------------------------------
print('')
print('=== Part 2: the n = 12 witness, two copies of D')
# ----------------------------------------------------------------------------------------

check('n12_the_two_label_sets_are_disjoint_and_their_union_is_1_to_12',
      set(COPY1).isdisjoint(COPY2) and set(COPY1) | set(COPY2) == set(range(1, 13)),
      'copy 1 %s ; copy 2 %s' % (sorted(COPY1), sorted(COPY2)))
check('n12_copy_1_all_six_weights_zero', weights(6, D_ARCS, COPY1) == [0] * 6,
      'f = %s' % (COPY1,))
check('n12_copy_2_all_six_weights_zero', weights(6, D_ARCS, COPY2) == [0] * 6,
      'f = %s' % (COPY2,))

# the disjoint union itself, on vertices 0..11
G12_ARCS = list(D_ARCS) + [(u + 6, v + 6) for (u, v) in D_ARCS]
F12 = list(COPY1) + list(COPY2)
check('n12_the_disjoint_union_is_an_oriented_graph_on_12_vertices',
      is_orientation(G12_ARCS) and len(G12_ARCS) == 18 and len(F12) == 12)
check('n12_all_twelve_weights_of_the_union_are_zero', weights(12, G12_ARCS, F12) == [0] * 12)
check('n12_f_is_a_bijection_onto_1_to_12', sorted(F12) == list(range(1, 13)))
check('n12_is_a_DDM_labeling_of_the_union',
      sorted(F12) == list(range(1, 13)) and weights(12, G12_ARCS, F12) == [0] * 12,
      'so the union is a DDMOG on 12 vertices')
check('n12_the_four_sum_triples_are_as_printed',
      (COPY1[1] + COPY1[2] == COPY1[0] and COPY1[4] + COPY1[5] == COPY1[3]
       and COPY2[1] + COPY2[2] == COPY2[0] and COPY2[4] + COPY2[5] == COPY2[3]),
      '1+5=6, 2+8=10, 3+9=12, 4+7=11')
check('n12_the_partition_into_four_sum_triples_covers_1_to_12',
      sorted([1, 5, 6] + [2, 8, 10] + [3, 9, 12] + [4, 7, 11]) == list(range(1, 13)))
check('n12_every_part_sum_is_even',
      all(s % 2 == 0 for s in (1 + 5 + 6, 2 + 8 + 10, 3 + 9 + 12, 4 + 7 + 11)),
      '12, 20, 24, 22; copy sums 32 and 46, both even')

# ----------------------------------------------------------------------------------------
print('')
print('=== Part 3: an n = 11 zero-weight labelled union, W u D')
# ----------------------------------------------------------------------------------------

check('W_component_of_the_n11_union_all_five_weights_zero', weights(5, W_ARCS, W_LABELS) == [0] * 5,
      'f = %s on label set %s' % (W_LABELS, sorted(W_LABELS)))
check('n11_D_component_all_six_weights_zero', weights(6, D_ARCS, COPY_A) == [0] * 6,
      'f = %s on label set %s' % (COPY_A, sorted(COPY_A)))
G11_ARCS = list(W_ARCS) + [(u + 5, v + 5) for (u, v) in D_ARCS]
F11 = list(W_LABELS) + list(COPY_A)
check('n11_f_is_a_bijection_onto_1_to_11', sorted(F11) == list(range(1, 12)))
check('n11_all_eleven_weights_of_the_union_are_zero', weights(11, G11_ARCS, F11) == [0] * 11,
      'so the union is a DDMOG on 11 vertices')
check('n11_the_W_component_label_set_is_not_an_interval',
      sorted(W_LABELS) != list(range(1, 6)) and set(W_LABELS) == {1, 3, 7, 10, 11},
      'a component of a disconnected DDMOG receives a SUBSET, not the interval {1,...,5}')

# ----------------------------------------------------------------------------------------
print('')
print('=== Part 4: no component is a DDMOG, and K_{3,3} is not DDMO')
# ----------------------------------------------------------------------------------------

k33_or = list(all_orientations(E_K33))
check('K33_has_512_orientations', len(k33_or) == 512 and all(is_orientation(a) for a in k33_or))

bf_hits = [a for a in k33_or if has_ddm(6, a)]
check('no_orientation_of_K33_is_a_DDMOG_brute_force_over_all_512_times_720', bf_hits == [],
      '0 of the 512 orientations admit a DDM labeling; 512 x 6! = 368640 bijections tested')

kern_hits = [a for a in k33_or if frozenset(range(1, 7)) in realisable(6, a, 1, 6)[0]]
check('no_orientation_of_K33_is_a_DDMOG_independent_exact_kernel_route', kern_hits == [],
      'the Fraction-RREF kernel enumeration agrees with the brute force on all 512 orientations')

odd_coeffs = all(
    all((sum(1 for (u, v) in a if u == b) - sum(1 for (u, v) in a if v == b)) % 2 == 1 for b in B_PART)
    for a in k33_or)
check('K33_parity_obstruction_every_B_vertex_has_odd_outdeg_minus_indeg', odd_coeffs and 21 % 2 == 1,
      'so S_B is even for any zero-weight labeling, while 1+...+6 = 21 is odd '
      '(= Altman et al. 2024, Theorem 5 at m=n=3: (3+3)(3+3+1)=42 is not a multiple of 4)')

check('the_exhibited_W_orientation_is_not_a_DDMOG', ddm_labelings(5, W_ARCS) == [],
      '0 of the 5! = 120 bijections onto {1,...,5}; by hand 15 = 3*f0 + 2*f2 has no admissible solution')

check('W4DDM_is_an_orientation_of_W4',
      is_orientation(W4DDM_ARCS) and len(canon(5, edges_of(W4DDM_ARCS))) == 8
      and canon(5, edges_of(W4DDM_ARCS)) == canon(5, E_W4),
      'underlying graph isomorphic to W_4')
check('W4_is_DDMO_the_exhibited_orientation_and_labeling_is_a_DDM_labeling',
      sorted(W4DDM_LABELS) == list(range(1, 6)) and weights(5, W4DDM_ARCS, W4DDM_LABELS) == [0] * 5,
      'f = %s, all weights 0 -- so W_4 IS DDMO (Altman et al. 2024, Theorem 10)' % (W4DDM_LABELS,))

# ----------------------------------------------------------------------------------------
print('')
print('=== Part 5: small-order enumerations (orders 4, 5 and the 5+5 split of 1..10)')
# ----------------------------------------------------------------------------------------

# order at most 4
check('min_degree_at_least_3_forces_order_at_least_4', all(n - 1 < 3 for n in (1, 2, 3)),
      'a graph on n <= 3 vertices has maximum degree n-1 < 3')

k4_or = list(all_orientations(E_K4))
k4_real = [a for a in k4_or if realisable(4, a, 1, 12)[0]]
check('order_4_no_orientation_of_K4_realises_any_injective_labeling_inside_1_to_12',
      len(k4_or) == 64 and k4_real == [],
      '0 of 64 orientations; a BOUNDED statement -- labels are bounded by 12 here, and no '
      'unbounded claim is tested')

# the order-5 classification
five = []
for mask in range(1 << 10):
    es = set(frozenset(ALL_PAIRS_5[i]) for i in range(10) if (mask >> i) & 1)
    if len(es) and min(degrees(5, es)) >= 3:
        five.append(es)
classes = {}
for es in five:
    classes.setdefault(canon(5, es), []).append(es)
check('order_5_exactly_26_labelled_graphs_have_min_degree_at_least_3', len(five) == 26,
      '%d labelled graphs on 5 vertices with min degree >= 3' % len(five))
check('order_5_they_form_exactly_3_isomorphism_classes_with_10_9_and_8_edges',
      len(classes) == 3 and sorted(len(v[0]) for v in classes.values()) == [8, 9, 10],
      'K_5 (10 edges), K_5-e (9), W_4 (8)')
check('order_5_the_three_representatives_used_below_are_those_three_classes',
      set(classes.keys()) == {canon(5, E_K5), canon(5, E_K5E), canon(5, E_W4)},
      'and they are pairwise non-isomorphic')
check('order_5_W4_representative_is_the_underlying_graph_of_W',
      E_W4 == edges_of(W_ARCS))

k5_or = list(all_orientations(E_K5))
k5e_or = list(all_orientations(E_K5E))
w4_or = list(all_orientations(E_W4))
check('order_5_orientation_counts_are_1024_512_and_256_totalling_1792',
      (len(k5_or), len(k5e_or), len(w4_or)) == (1024, 512, 256)
      and len(k5_or) + len(k5e_or) + len(w4_or) == 1792)

bad = []
for a in k5_or + k5e_or:
    if realisable(5, a, 1, 12)[0]:
        bad.append(a)
check('order_5_no_orientation_of_K5_or_K5_minus_e_realises_ANY_label_set_inside_1_to_12',
      bad == [],
      '0 of the 1536 orientations of the two non-DDMO order-5 graphs admit an injective '
      'positive zero-weight labeling with labels <= 12')

w4_real12 = [a for a in w4_or if realisable(5, a, 1, 12)[0]]
check('order_5_exactly_18_of_the_1792_orientations_realise_anything_inside_1_to_12',
      len(w4_real12) == 18,
      'all 18 are orientations of W_4 (%d of its 256); this is the WEAKER predicate '
      '"realises some label set", not "is a DDMOG"' % len(w4_real12))

w4_ddmog = [a for a in w4_or if ddm_labelings(5, a)]
w4_lab_total = sum(len(ddm_labelings(5, a)) for a in w4_ddmog)
aut = [p for p in itertools.permutations(range(5))
       if set(frozenset((p[a], p[b])) for (a, b) in (tuple(e) for e in E_W4)) == E_W4]
orbits = {}
for a in w4_ddmog:
    orbits.setdefault(min(tuple(sorted((p[u], p[v]) for (u, v) in a)) for p in aut), []).append(a)
check('order_5_Aut_W4_has_order_8', len(aut) == 8,
      'the hub is the unique vertex of degree 4, so Aut(W_4) = Aut(C_4) = D_4')
note('order-5 DDMOG census: %d of the 256 labelled orientations of W_4 are DDMOGs; they carry '
     '%d DDM labelings in total; they fall into %d isomorphism class(es) with orbit size(s) %s '
     'under Aut(W_4), and %s DDM labeling(s) per class representative.'
     % (len(w4_ddmog), w4_lab_total, len(orbits), sorted(len(v) for v in orbits.values()),
        sorted(len(ddm_labelings(5, v[0])) for v in orbits.values())))
check('order_5_every_DDMOG_has_underlying_graph_W4',
      bad == [] and len(w4_ddmog) > 0,
      'no orientation of K_5 or of K_5-e is a DDMOG, while %d orientations of W_4 are -- the '
      'first sentence of Altman et al. 2024, Theorem 10' % len(w4_ddmog))
check('order_5_W4_has_10_DDMOG_orientations_in_exactly_2_isomorphism_classes',
      len(w4_ddmog) == 10 and len(orbits) == 2
      and sorted(len(v) for v in orbits.values()) == [2, 8],
      'orbit sizes 2 and 8 under the order-8 group Aut(W_4); "only two unique orientations" '
      'of Altman et al. 2024, Theorem 10')

# The 48 (orientation, DDM labeling) pairs, reduced by Aut(W_4).  The action is FREE: an
# automorphism fixing a pair fixes an injective labeling, hence is the identity.  So the number
# of classes is 48/8 = 6, which is the "exactly 6 such labelings" of Altman et al. Theorem 10.
pairs = [(a, f) for a in w4_ddmog for f in ddm_labelings(5, a)]
pair_orbits = set()
free_action = True
for (a, f) in pairs:
    best = None
    fixers = 0
    for p in aut:
        g = [0] * 5
        for u in range(5):
            g[p[u]] = f[u]
        key = (tuple(sorted((p[u], p[v]) for (u, v) in a)), tuple(g))
        if key == (tuple(sorted(a)), tuple(f)):
            fixers += 1
        if best is None or key < best:
            best = key
    if fixers != 1:
        free_action = False
    pair_orbits.add(best)
check('order_5_the_48_W4_DDM_labelings_form_exactly_6_classes_up_to_isomorphism',
      len(pairs) == 48 and free_action and len(pair_orbits) == 6,
      'Aut(W_4) acts freely on the 48 (orientation, labeling) pairs, so 48/8 = 6 classes -- '
      '"There are exactly 6 such labelings on only two unique orientations", Altman et al. '
      '2024, Theorem 10, reproduced in full')

# the n = 10 lane, split 5 + 5
all_sets10 = set()
nonddmog_sets10 = set()
for a in w4_or:
    s = realisable(5, a, 1, 10)[0]
    all_sets10 |= s
    if frozenset(range(1, 6)) not in s:        # this orientation is NOT a DDMOG
        nonddmog_sets10 |= s
for a in k5_or + k5e_or:                        # realise nothing, checked above; kept for completeness
    s = realisable(5, a, 1, 10)[0]
    all_sets10 |= s
    if frozenset(range(1, 6)) not in s:
        nonddmog_sets10 |= s
check('n10_five_subsets_of_1_to_10_realisable_by_some_order_5_orientation_number_99',
      len(all_sets10) == 99, '%d of the 252 five-subsets' % len(all_sets10))
check('n10_those_realisable_by_a_NON_DDMOG_orientation_number_29',
      len(nonddmog_sets10) == 29,
      '%d -- the subtotal the 5+5 join uses; reading the 29 as the total understates by 3.4x'
      % len(nonddmog_sets10))
full10 = frozenset(range(1, 11))
comp_pairs = [s for s in nonddmog_sets10 if (full10 - s) in nonddmog_sets10]
check('n10_no_5plus5_split_of_1_to_10_has_both_halves_carried_by_a_non_DDMOG_order_5_orientation',
      comp_pairs == [],
      'no 5+5 split of {1,...,10} has both halves carried by a non-DDMOG orientation')





# ----------------------------------------------------------------------------------------
print('')
print('=== Part 6: the source Theorem 5.1 corollary and a k-copies parity arithmetic check')
# ----------------------------------------------------------------------------------------

check('corollary_source_theorem_5_1_is_not_an_if_and_only_if',
      weights(12, G12_ARCS, F12) == [0] * 12 and sorted(F12) == list(range(1, 13))
      and bf_hits == [],
      'the union D u D is a DDMOG while NO component is a DDMOG, so the hypothesis of that '
      'theorem fails and its conclusion holds')

# k copies of D: total label sum 1+...+6k = 3k(6k+1) must be even, since every admissible
# 6-set of D is a union of two sum-triples {a,b,a+b}, each of even sum.
parity_ok = all(((3 * k * (6 * k + 1)) % 2 == 0) == (k % 2 == 0) for k in range(1, 25))
check('k_copies_of_D_parity_condition_3k(6k+1)_even_iff_k_even_verified_for_k_1_to_24',
      parity_ok,
      'the arithmetic equivalence 3k(6k+1) even iff k even is machine-checked for '
      'k = 1..24 only; no general implication is tested')
check('k_equals_2_is_attained_and_k_equals_1_is_impossible',
      weights(12, G12_ARCS, F12) == [0] * 12 and bf_hits == [],
      'k = 2 is the n = 12 witness; k = 1 fails because K_{3,3} is not DDMO')

# ----------------------------------------------------------------------------------------
print('')
note('SCOPE -- what this program does NOT re-run.  (1) The order-6 census of Altman et al. '
     '(2024): their 22 unique 6-vertex oriented DDMOGs, and the 19 isomorphism classes of '
     '6-vertex graphs of minimum degree >= 3, are NOT enumerated here, and no agreement with '
     'them is asserted.  The only order-6 graph the note uses is K_{3,3}, whose 512 '
     'orientations are swept exhaustively in Part 4.  '
     '(2) Component orders 7 and above are never enumerated.  '
     '(3) No lane containing an isolated (order-1) component is enumerated, and no minimum '
     'order and no witness count is established here.  '
     '(4) The sufficiency half of the k-even remark rests on the classical triple-packing '
     'construction of Guy (1976) and is NOT verified here; only its necessity half, and the '
     'case k = 2, are established above.  (5) The order-4 check above bounds labels by 12 '
     'and establishes nothing for unbounded labels.')
print('')

if FAILED:
    print('FAILED CHECKS: %d' % len(FAILED))
    for x in FAILED:
        print('  FAIL %s' % x)
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(PASSED))
sys.exit(0)

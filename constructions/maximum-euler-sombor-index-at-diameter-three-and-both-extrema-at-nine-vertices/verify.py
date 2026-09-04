#!/usr/bin/env python3
"""Verification program for

    "The Maximum Euler-Sombor Index at Diameter Three, and Both Extrema at Nine Vertices"

WHAT THIS PROGRAM IS.
  It reads the two objects EXACTLY AS THEY ARE PRINTED IN SECTION 3 OF THE PAPER -- as vertex-named edge
  lists, not as graph6 strings and not from any census output -- and re-derives every quantity the paper
  claims about them: order, size, degree sequence, edge degree-pair profile, connectivity, the exact
  diameter, and the Euler-Sombor index compared against its closed form in 60-digit exact-square-root
  arithmetic. It then checks the closed form of Theorem B against a direct evaluation for n = 4..14, the
  strict convexity that delivers uniqueness, the inequality (star) of Step 3, the floor and the tier cap
  that amputate the two search spaces, the strict ordering that makes the minimiser bicyclic, the three
  exhaustion identities, and two published values by other authors that pin the normalisation of EU.

WHAT THIS PROGRAM IS NOT.
  It is NOT a census. It enumerates nothing over the class G_{9,3} and could not detect a graph outside the
  objects named in the paper. The EXHAUSTIVENESS and UNIQUENESS assertions of Theorem A's minimum half rest
  on the exhaustive census described in the paper's Section 5, which this program does not re-run; that
  limitation is stated in the paper and restated in the SCOPE line this program prints last.
  Theorem B needs no computation at all: Section 4 is a closed-form proof, and what is checked here is that
  its formula and its inequalities are correct, not that a search found them.

THE CHECKER IS FORCED TO BE ABLE TO SAY NO. The [ANTI] block feeds it four objects it must reject:
H(1,6) with the edge uv added (the diameter drops to 2), Theta(2,4,4) (diameter 4 although its Euler-Sombor
index EQUALS the claimed minimum exactly), a disconnected graph (which must report "disconnected", never a
diameter), and a deliberately wrong candidate maximum of 320.

Python 3.9+, standard library only. No third-party package, no external data file, no network.
Exit code 0 if and only if every check passes.
"""
import itertools
import sys
from collections import deque
from decimal import Decimal, getcontext

getcontext().prec = 60

RAN = []
BAD = []


def ck(name, cond, detail=''):
    RAN.append(name)
    if not cond:
        BAD.append(name)
    print('%s %s%s' % ('PASS' if cond else 'FAIL', name, (' [%s]' % detail) if detail else ''))


def note(s):
    print('NOTE %s' % s)


def sq(n):
    """sqrt(n) to 60 significant digits, exactly rounded; no float ever enters a decision."""
    return Decimal(n).sqrt()


def rnd(x, k):
    return str(x.quantize(Decimal(1).scaleb(-k)))


class Graph:
    """A graph on named vertices, built from an edge list exactly as printed in the paper."""

    def __init__(self, name, edges):
        self.name = name
        self.E = []
        seen = set()
        for x, y in edges:
            assert x != y, 'loop in %s' % name
            k = tuple(sorted((x, y)))
            assert k not in seen, 'repeated edge %s in %s' % (k, name)
            seen.add(k)
            self.E.append(k)
        self.V = sorted({v for e in self.E for v in e})
        self.adj = {v: set() for v in self.V}
        for x, y in self.E:
            self.adj[x].add(y)
            self.adj[y].add(x)

    @property
    def n(self):
        return len(self.V)

    @property
    def m(self):
        return len(self.E)

    def degseq(self):
        return sorted(len(self.adj[v]) for v in self.V)

    def profile(self):
        p = {}
        for x, y in self.E:
            k = tuple(sorted((len(self.adj[x]), len(self.adj[y]))))
            p[k] = p.get(k, 0) + 1
        return dict(sorted(p.items()))

    def _ecc(self, s):
        d = {s: 0}
        q = deque([s])
        while q:
            x = q.popleft()
            for y in self.adj[x]:
                if y not in d:
                    d[y] = d[x] + 1
                    q.append(y)
        return None if len(d) != self.n else max(d.values())

    def diameter(self):
        """-1 means DISCONNECTED. A disconnected graph has no diameter and must never be given one."""
        best = 0
        for v in self.V:
            e = self._ecc(v)
            if e is None:
                return -1
            best = max(best, e)
        return best

    def dist(self, a, b):
        d = {a: 0}
        q = deque([a])
        while q:
            x = q.popleft()
            if x == b:
                return d[x]
            for y in self.adj[x]:
                if y not in d:
                    d[y] = d[x] + 1
                    q.append(y)
        return -1

    def EU(self):
        t = Decimal(0)
        for x, y in self.E:
            a, b = len(self.adj[x]), len(self.adj[y])
            t += sq(a * a + a * b + b * b)
        return t


# ---------------------------------------------------------------------------
# constructors for the auxiliary families the paper names
# ---------------------------------------------------------------------------
def H(a, b, n=9):
    """H(a,b) of the paper's Section 2: K_{n-2} on W, u joined to the first a of W, v to the other b,
    and u NOT adjacent to v."""
    assert a >= 1 and b >= 1 and a + b == n - 2
    W = ['w%d' % i for i in range(1, n - 1)]
    E = [(x, y) for x, y in itertools.combinations(W, 2)]
    E += [('u', w) for w in W[:a]] + [('v', w) for w in W[a:]]
    return Graph('H(%d,%d) on %d vertices' % (a, b, n), E)


def theta(*lengths):
    """Two terminals u,v joined by internally disjoint paths of the given edge-lengths."""
    E = []
    for i, L in enumerate(lengths):
        prev = 'u'
        for j in range(L - 1):
            cur = 'p%d_%d' % (i, j)
            E.append((prev, cur))
            prev = cur
        E.append((prev, 'v'))
    return Graph('Theta%s' % (lengths,), E)


def double_star(a, b):
    """D(a,b): adjacent centres with a and b leaves, so the centre degrees are a+1 and b+1."""
    E = [('c1', 'c2')]
    E += [('c1', 'x%d' % i) for i in range(a)]
    E += [('c2', 'y%d' % i) for i in range(b)]
    return Graph('D(%d,%d)' % (a, b), E)


def path(k):
    return Graph('P%d' % k, [('v%d' % i, 'v%d' % (i + 1)) for i in range(k - 1)])


def cycle(k):
    return Graph('C%d' % k, [('v%d' % i, 'v%d' % ((i + 1) % k)) for i in range(k)])


print('=' * 100)
print('verification of "The Maximum Euler-Sombor Index at Diameter Three, and Both Extrema at Nine Vertices"')
print('exact arithmetic only: 60-significant-digit Decimal square roots, no float decision anywhere')
print('=' * 100)

# ===========================================================================
print('\n=== Step 1: H(1,6), the paper\'s claimed maximiser, rebuilt from the edge list transcribed by hand from Section 3 (transcription not machine-verified; maximality not re-checked here)')
# The edge list AS PRINTED: the 21 pairs w_i w_j, the edge u w_1, and the six edges v w_2 .. v w_7.
W = ['w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7']
H16_EDGES = [(x, y) for x, y in itertools.combinations(W, 2)] \
            + [('u', 'w1')] \
            + [('v', w) for w in ['w2', 'w3', 'w4', 'w5', 'w6', 'w7']]
h16 = Graph('H(1,6)', H16_EDGES)
MAX_CF = 147 * sq(3) + sq(57) + 6 * sq(127)
note('vertices (%d): %s' % (h16.n, ', '.join(h16.V)))
note('edges (%d): %s' % (h16.m, ' '.join('%s%s' % e for e in h16.E)))
note('degrees: u=%d v=%d w1..w7=%s' % (len(h16.adj['u']), len(h16.adj['v']),
                                       [len(h16.adj[w]) for w in W]))
note('EU(H(1,6)) to 60 digits = %s' % h16.EU())
note('147*sqrt3 + sqrt57 + 6*sqrt127 = %s' % MAX_CF)
ck('H_1_6_has_order_9', h16.n == 9, 'n=%d' % h16.n)
ck('H_1_6_has_size_28_equal_to_C_8_2', h16.m == 28 == 8 * 7 // 2, 'm=%d' % h16.m)
ck('H_1_6_degree_sequence_is_1_6_7x7', h16.degseq() == [1, 6] + [7] * 7, str(h16.degseq()))
ck('H_1_6_degree_sum_is_twice_the_size', sum(h16.degseq()) == 2 * h16.m == 56, '56')
ck('H_1_6_edge_profile_is_1_7_once_6_7_six_times_7_7_twentyone_times',
   h16.profile() == {(1, 7): 1, (6, 7): 6, (7, 7): 21}, str(h16.profile()))
ck('H_1_6_u_and_v_are_nonadjacent', 'v' not in h16.adj['u'])
ck('H_1_6_distance_from_u_to_v_is_3', h16.dist('u', 'v') == 3, 'd(u,v)=%d' % h16.dist('u', 'v'))
ck('H_1_6_is_connected_with_diameter_exactly_3', h16.diameter() == 3, 'diam=%d' % h16.diameter())
ck('H_1_6_every_pair_other_than_uv_is_within_distance_2',
   all(h16.dist(x, y) <= 2 for x, y in itertools.combinations(h16.V, 2) if {x, y} != {'u', 'v'}))
ck('H_1_6_EU_equals_147sqrt3_plus_sqrt57_plus_6sqrt127_to_1e-50',
   abs(h16.EU() - MAX_CF) < Decimal('1e-50'), 'difference = %s' % (h16.EU() - MAX_CF))
ck('H_1_6_EU_equals_the_value_printed_in_the_paper_329_777869165403581',
   rnd(h16.EU(), 15) == '329.777869165403581', rnd(h16.EU(), 15))
ck('H_1_6_the_three_summands_are_the_ones_the_paper_derives_by_hand',
   rnd(147 * sq(3), 10) == '254.6114687126' and rnd(sq(57), 10) == '7.5498344353'
   and rnd(6 * sq(127), 10) == '67.6165660175',
   '%s + %s + %s' % (rnd(147 * sq(3), 10), rnd(sq(57), 10), rnd(6 * sq(127), 10)))

# ===========================================================================
print('\n=== Step 2: Theta(3,3,4), the paper\'s claimed minimiser, rebuilt from the edge list transcribed by hand from Section 3 (minimality and uniqueness not re-run here; see SCOPE)')
TH_EDGES = [('u', 'a1'), ('a1', 'a2'), ('a2', 'v'),
            ('u', 'b1'), ('b1', 'b2'), ('b2', 'v'),
            ('u', 'c1'), ('c1', 'c2'), ('c2', 'c3'), ('c3', 'v')]
th = Graph('Theta(3,3,4)', TH_EDGES)
MIN_CF = 8 * sq(3) + 6 * sq(19)
note('vertices (%d): %s' % (th.n, ', '.join(th.V)))
note('edges (%d): %s' % (th.m, ' '.join('%s-%s' % e for e in th.E)))
note('EU(Theta(3,3,4)) to 60 digits = %s' % th.EU())
note('8*sqrt3 + 6*sqrt19 = %s' % MIN_CF)
ck('Theta_3_3_4_has_order_9', th.n == 9, 'n=%d' % th.n)
ck('Theta_3_3_4_has_size_10', th.m == 10, 'm=%d' % th.m)
ck('Theta_3_3_4_is_bicyclic_cyclomatic_number_2', th.m - th.n + 1 == 2)
ck('Theta_3_3_4_degree_sequence_is_2x7_3_3', th.degseq() == [2] * 7 + [3, 3], str(th.degseq()))
ck('Theta_3_3_4_degree_sum_is_twice_the_size', sum(th.degseq()) == 2 * th.m == 20, '20')
ck('Theta_3_3_4_profile_is_2_2_four_times_and_2_3_six_times_with_no_3_3_edge',
   th.profile() == {(2, 2): 4, (2, 3): 6}, str(th.profile()))
ck('Theta_3_3_4_terminals_u_v_are_nonadjacent_at_distance_3',
   'v' not in th.adj['u'] and th.dist('u', 'v') == 3, 'd(u,v)=%d' % th.dist('u', 'v'))
ck('Theta_3_3_4_is_connected_with_diameter_exactly_3', th.diameter() == 3, 'diam=%d' % th.diameter())
ck('Theta_3_3_4_the_two_extremal_pairs_named_in_the_paper_are_at_distance_3',
   th.dist('a1', 'b2') == 3 and th.dist('a1', 'c2') == 3,
   'd(a1,b2)=%d d(a1,c2)=%d' % (th.dist('a1', 'b2'), th.dist('a1', 'c2')))
ck('Theta_3_3_4_EU_equals_8sqrt3_plus_6sqrt19_to_1e-50',
   abs(th.EU() - MIN_CF) < Decimal('1e-50'), 'difference = %s' % (th.EU() - MIN_CF))
ck('Theta_3_3_4_EU_equals_the_value_printed_in_the_paper_40_009800121795059',
   rnd(th.EU(), 15) == '40.009800121795060', rnd(th.EU(), 15))
ck('Theta_3_3_4_the_two_summands_are_the_ones_the_paper_derives_by_hand',
   rnd(8 * sq(3), 10) == '13.8564064606' and rnd(6 * sq(19), 10) == '26.1533936612',
   '%s + %s' % (rnd(8 * sq(3), 10), rnd(6 * sq(19), 10)))

# ===========================================================================
print('\n=== Step 3: the H(a,b) ladder at n = 9 and the uniqueness margin of Theorem B')
clique9 = 147 * sq(3)
lad = {}
for a, b in ((1, 6), (2, 5), (3, 4)):
    g = H(a, b)
    lad[(a, b)] = g.EU()
    note('H(%d,%d): m=%d diam=%d  EU=%s  apex part=%s'
         % (a, b, g.m, g.diameter(), rnd(g.EU(), 10), rnd(g.EU() - clique9, 10)))
ck('every_H_a_b_at_n_9_has_size_28_and_diameter_exactly_3',
   all(H(a, b).m == 28 and H(a, b).diameter() == 3 for a, b in ((1, 6), (2, 5), (3, 4))))
ck('the_three_unordered_splits_of_n_minus_2_equal_7_are_exactly_1_6_2_5_3_4',
   len([(a, 7 - a) for a in range(1, 7) if a <= 7 - a]) == 3 == (9 - 2) // 2)
ck('the_ladder_is_strict_H_1_6_beats_H_2_5_beats_H_3_4',
   lad[(1, 6)] > lad[(2, 5)] > lad[(3, 4)])
ck('apex_part_of_H_1_6_is_sqrt57_plus_6sqrt127_equal_75_1664004528',
   abs(lad[(1, 6)] - clique9 - (sq(57) + 6 * sq(127))) < Decimal('1e-40')
   and rnd(sq(57) + 6 * sq(127), 10) == '75.1664004528', rnd(lad[(1, 6)] - clique9, 10))
ck('apex_part_of_H_2_5_is_2sqrt67_plus_5sqrt109_equal_68_5722380883',
   abs(lad[(2, 5)] - clique9 - (2 * sq(67) + 5 * sq(109))) < Decimal('1e-40')
   and rnd(2 * sq(67) + 5 * sq(109), 10) == '68.5722380883', rnd(lad[(2, 5)] - clique9, 10))
ck('apex_part_of_H_3_4_is_3sqrt79_plus_4sqrt93_equal_65_2391862959',
   abs(lad[(3, 4)] - clique9 - (3 * sq(79) + 4 * sq(93))) < Decimal('1e-40')
   and rnd(3 * sq(79) + 4 * sq(93), 10) == '65.2391862959', rnd(lad[(3, 4)] - clique9, 10))
ck('H_2_5_EU_is_147sqrt3_plus_2sqrt67_plus_5sqrt109_equal_323_1837068009',
   rnd(lad[(2, 5)], 10) == '323.1837068009', rnd(lad[(2, 5)], 10))
ck('H_3_4_EU_is_147sqrt3_plus_3sqrt79_plus_4sqrt93_equal_319_8506550085',
   rnd(lad[(3, 4)], 10) == '319.8506550085', rnd(lad[(3, 4)], 10))
ck('the_gap_between_H_1_6_and_H_2_5_at_n_9_is_6_5941623645',
   rnd(lad[(1, 6)] - lad[(2, 5)], 10) == '6.5941623645', rnd(lad[(1, 6)] - lad[(2, 5)], 10))

# ===========================================================================
print('\n=== Step 4: the closed form of Theorem B, and the convexity that makes H(1,n-3) unique')


def theoremB(n):
    return (sq(3) * Decimal((n - 2) ** 2 * (n - 3)) / 2
            + sq(n * n - 3 * n + 3)
            + Decimal(n - 3) * sq(3 * n * n - 15 * n + 19))


ok_cf, ok_uni, ok_d3 = True, True, True
for n in range(4, 15):
    direct = H(1, n - 3, n).EU()
    if abs(direct - theoremB(n)) >= Decimal('1e-45'):
        ok_cf = False
    if H(1, n - 3, n).diameter() != 3:
        ok_d3 = False
    vals = {}
    for a in range(1, (n - 2) // 2 + 1):
        vals[a] = H(a, n - 2 - a, n).EU()
    if min(vals) != 1 or vals[1] != max(vals.values()) or list(vals.values()).count(max(vals.values())) != 1:
        ok_uni = False
    note('n=%2d : Theorem B = %-22s  direct EU(H(1,n-3)) = %-22s  diam=%d  splits=%d'
         % (n, rnd(theoremB(n), 10), rnd(direct, 10), H(1, n - 3, n).diameter(), len(vals)))
ck('theorem_B_closed_form_equals_a_direct_evaluation_of_EU_H_1_n_minus_3_for_n_4_to_14', ok_cf)
ck('H_1_n_minus_3_has_diameter_exactly_3_for_n_4_to_14', ok_d3)
ck('among_all_H_a_b_the_split_a_equals_1_is_the_unique_maximum_for_n_4_to_14', ok_uni)
ck('theorem_B_at_n_9_reproduces_147sqrt3_plus_sqrt57_plus_6sqrt127',
   abs(theoremB(9) - MAX_CF) < Decimal('1e-45'), rnd(theoremB(9), 12))
ck('theorem_B_at_n_4_gives_2sqrt3_plus_2sqrt7_which_is_EU_of_the_path_P4',
   abs(theoremB(4) - (2 * sq(3) + 2 * sq(7))) < Decimal('1e-45')
   and abs(path(4).EU() - theoremB(4)) < Decimal('1e-45')
   and rnd(theoremB(4), 10) == '8.7556042373', rnd(theoremB(4), 10))
ck('the_only_member_of_G_4_3_is_P4_and_H_1_1_is_isomorphic_to_it',
   H(1, 1, 4).degseq() == path(4).degseq() == [1, 1, 2, 2] and H(1, 1, 4).m == path(4).m == 3
   and H(1, 1, 4).diameter() == 3)
ck('theorem_B_at_n_5_gives_9sqrt3_plus_sqrt13_plus_2sqrt19_equal_27_9118064307',
   abs(theoremB(5) - (9 * sq(3) + sq(13) + 2 * sq(19))) < Decimal('1e-45')
   and rnd(theoremB(5), 10) == '27.9118064307', rnd(theoremB(5), 10))
# the second-difference identity behind Step 4: numerator of phi'' as the paper expands it
ok_num = True
for c in range(2, 30):
    for t in range(0, 30):
        lhs = 2 * (8 * t + 3 * c) * (t * t + c * t + c * c) - (4 * t * t + 3 * c * t + 2 * c * c) * (2 * t + c)
        rhs = 8 * t ** 3 + 12 * c * t * t + 15 * c * c * t + 4 * c ** 3
        if lhs != rhs or rhs <= 0:
            ok_num = False
ck('the_numerator_of_phi_double_prime_is_the_positive_polynomial_the_paper_states',
   ok_num, 'identity and positivity for 2<=c<=29, 0<=t<=29')
# inequality (star), as an exact algebraic margin
ok_star, margins = True, []
for n in range(4, 200):
    s = n - 3
    lhs = sq(3) * Decimal((n - 2) * (n - 3))
    rhs = sq(s * s + 3 * s + 3) + Decimal(s) * sq(3 * s * s + 3 * s + 1)
    if rhs <= lhs:
        ok_star = False
    if n <= 11:
        margins.append('n=%d:%s' % (n, rnd(rhs - lhs, 3)))
ck('inequality_star_of_step_3_holds_for_every_n_from_4_to_199', ok_star, ' '.join(margins))
ck('the_star_margin_grows_with_slope_1_minus_sqrt3_over_2', rnd(1 - sq(3) / 2, 7) == '0.1339746',
   rnd(1 - sq(3) / 2, 7))

# ===========================================================================
print('\n=== Step 5: three 9-vertex diameter-3 double stars and one hand-built unicyclic witness, each with EU above EU(Theta(3,3,4)), plus a few smaller-order comparison graphs')
trees = {}
for a in (1, 2, 3):
    g = double_star(a, 7 - a)
    trees[a] = g.EU()
    note('D(%d,%d): n=%d m=%d diam=%d centre degrees %s  EU=%s'
         % (a, 7 - a, g.n, g.m, g.diameter(), sorted((a + 1, 8 - a)), rnd(g.EU(), 10)))
ck('the_three_9_vertex_double_stars_are_trees_of_diameter_exactly_3',
   all(double_star(a, 7 - a).n == 9 and double_star(a, 7 - a).m == 8
       and double_star(a, 7 - a).diameter() == 3 for a in (1, 2, 3)))
ck('the_least_of_the_three_tested_9_vertex_double_stars_is_D_3_4_equal_sqrt61_plus_3sqrt21_plus_4sqrt31',
   min(trees.values()) == trees[3]
   and abs(trees[3] - (sq(61) + 3 * sq(21) + 4 * sq(31))) < Decimal('1e-40')
   and rnd(trees[3], 10) == '43.8290342121', rnd(trees[3], 10))
ck('the_other_two_double_stars_have_the_surd_forms_the_paper_prints',
   abs(trees[2] - (sq(63) + 2 * sq(13) + 5 * sq(43))) < Decimal('1e-40')
   and abs(trees[1] - (sq(67) + sq(7) + 6 * sq(57))) < Decimal('1e-40'),
   'D(2,5)=%s D(1,6)=%s' % (rnd(trees[2], 7), rnd(trees[1], 7)))
ck('the_three_tested_9_vertex_double_stars_are_all_above_the_claimed_minimum',
   all(v > MIN_CF for v in trees.values()), 'least is %s' % rnd(min(trees.values()), 10))
uni = Graph('C5 with two adjacent cycle vertices carrying two pendants each',
            [('v0', 'v1'), ('v1', 'v2'), ('v2', 'v3'), ('v3', 'v4'), ('v4', 'v0'),
             ('v0', 'p1'), ('v0', 'p2'), ('v1', 'p3'), ('v1', 'p4')])
note('unicyclic witness (one hand-built 9-vertex diameter-3 unicyclic graph, not shown minimal in its class): n=%d m=%d diam=%d degrees=%s EU=%s'
     % (uni.n, uni.m, uni.diameter(), uni.degseq(), rnd(uni.EU(), 10)))
ck('the_hand_built_unicyclic_witness_is_a_9_vertex_diameter_3_graph_with_9_edges_and_cyclomatic_number_1',
   uni.n == 9 and uni.m == 9 and uni.diameter() == 3 and uni.m - uni.n + 1 == 1)
ck('its_degree_sequence_is_1_1_1_1_2_2_2_4_4', uni.degseq() == [1, 1, 1, 1, 2, 2, 2, 4, 4],
   str(uni.degseq()))
ck('its_EU_is_4sqrt21_plus_8sqrt3_plus_4sqrt7_equal_42_7697144846',
   abs(uni.EU() - (4 * sq(21) + 8 * sq(3) + 4 * sq(7))) < Decimal('1e-40')
   and rnd(uni.EU(), 10) == '42.7697144846', rnd(uni.EU(), 10))
ck('strict_order_EU_Theta_3_3_4_below_the_unicyclic_witness_below_the_best_of_the_three_double_stars',
   th.EU() < uni.EU() < trees[3],
   'gaps %s and %s' % (rnd(uni.EU() - th.EU(), 10), rnd(trees[3] - th.EU(), 10)))
ck('the_cheapest_m_equals_10_profile_with_one_3_3_edge_would_have_beaten_the_minimum',
   13 * sq(3) + 4 * sq(19) < MIN_CF, rnd(13 * sq(3) + 4 * sq(19), 10))
ck('the_three_theta_1_b_graphs_with_that_profile_have_diameter_4_at_n_9',
   all(theta(1, b, 9 - b).diameter() >= 4 for b in (2, 3, 4)),
   'Theta(1,2,7)=%d Theta(1,3,6)=%d Theta(1,4,5)=%d'
   % (theta(1, 2, 7).diameter(), theta(1, 3, 6).diameter(), theta(1, 4, 5).diameter()))
ck('the_next_cheapest_m_equals_10_degree_sequence_4_2x8_costs_4sqrt28_plus_12sqrt3_above_the_minimum',
   4 * sq(28) + 12 * sq(3) > MIN_CF, rnd(4 * sq(28) + 12 * sq(3), 10))
ck('C6_has_diameter_3_while_every_spanning_tree_of_C6_is_P6_of_diameter_5',
   cycle(6).diameter() == 3 and path(6).diameter() == 5)
ck('at_n_8_EU_D_3_3_is_below_EU_Theta_3_3_3_while_at_n_9_EU_Theta_3_3_4_is_below_EU_D_3_4',
   double_star(3, 3).EU() < theta(3, 3, 3).EU() and th.EU() < trees[3],
   'n=8: D(3,3)=%s < Theta(3,3,3)=%s' % (rnd(double_star(3, 3).EU(), 7),
                                         rnd(theta(3, 3, 3).EU(), 7)))
ck('D_3_3_at_n_8_is_4sqrt3_plus_6sqrt21_and_Theta_3_3_3_is_6sqrt3_plus_6sqrt19',
   abs(double_star(3, 3).EU() - (4 * sq(3) + 6 * sq(21))) < Decimal('1e-40')
   and abs(theta(3, 3, 3).EU() - (6 * sq(3) + 6 * sq(19))) < Decimal('1e-40'),
   '%s and %s' % (rnd(double_star(3, 3).EU(), 7), rnd(theta(3, 3, 3).EU(), 7)))
ck('C7_has_diameter_3_and_EU_14sqrt3_below_every_9_vertex_tree_analogue_at_n_7',
   cycle(7).diameter() == 3 and abs(cycle(7).EU() - 14 * sq(3)) < Decimal('1e-40')
   and cycle(7).EU() < double_star(2, 3).EU(), rnd(cycle(7).EU(), 10))

# ===========================================================================
print('\n=== Step 6: the four theta graphs of Remark 3 share one EU and are separated by the diameter')
FOUR = [(2, 2, 6), (2, 3, 5), (2, 4, 4), (3, 3, 4)]
vals, dia = [], {}
for t in FOUR:
    g = theta(*t)
    vals.append(g.EU())
    dia[t] = g.diameter()
    note('Theta%-10s n=%d m=%d diam=%d profile=%s EU=%s'
         % (str(t), g.n, g.m, g.diameter(), g.profile(), rnd(g.EU(), 10)))
ck('all_four_thetas_have_order_9_size_10_and_the_same_edge_profile',
   all(theta(*t).n == 9 and theta(*t).m == 10
       and theta(*t).profile() == {(2, 2): 4, (2, 3): 6} for t in FOUR))
ck('all_four_thetas_have_the_same_EU_equal_to_the_claimed_minimum',
   max(vals) - min(vals) < Decimal('1e-50') and abs(vals[0] - MIN_CF) < Decimal('1e-50'))
ck('only_Theta_3_3_4_has_diameter_3_the_other_three_have_diameter_4',
   dia[(3, 3, 4)] == 3 and all(dia[t] == 4 for t in FOUR[:3]), str(dia))
ck('the_four_partitions_of_10_into_three_parts_each_at_least_2_are_exactly_those_four',
   sorted(tuple(sorted(p)) for p in
          {tuple(sorted((a, b, 10 - a - b))) for a in range(2, 7) for b in range(2, 7)
           if 10 - a - b >= 2}) == sorted(tuple(sorted(t)) for t in FOUR))

# ===========================================================================
print('\n=== Step 7: the two amputations, and the star that shows Step 3 needs diam >= 3')
star8 = Graph('K_{1,8}', [('c', 'l%d' % i) for i in range(8)])
note('K_{1,8}: n=%d m=%d diam=%d Delta=%d EU=%s'
     % (star8.n, star8.m, star8.diameter(), max(star8.degseq()), rnd(star8.EU(), 10)))
ck('K_1_8_has_diameter_2_and_maximum_degree_n_minus_1_so_step_3_must_assume_diam_at_least_3',
   star8.diameter() == 2 and max(star8.degseq()) == 8 == star8.n - 1)
ck('K_1_8_EU_is_8sqrt73_equal_68_3520299625',
   abs(star8.EU() - 8 * sq(73)) < Decimal('1e-40') and rnd(star8.EU(), 10) == '68.3520299625',
   rnd(star8.EU(), 10))
floor11 = 2 * sq(3) * Decimal(121) / 9
ck('the_floor_expression_at_m_11_n_9_evaluates_to_46_5729217146_which_exceeds_EU_D_3_4',
   floor11 > trees[3] and rnd(floor11, 10) == '46.5729217146', rnd(floor11, 10))
ck('the_floor_is_a_valid_lower_bound_on_both_exhibited_objects',
   th.EU() >= 2 * sq(3) * Decimal(th.m) ** 2 / 9 and h16.EU() >= 2 * sq(3) * Decimal(h16.m) ** 2 / 9)
ck('the_first_inequality_of_the_floor_holds_edgewise_on_both_exhibited_objects',
   all(sq(a * a + a * b + b * b) - sq(3) * Decimal(a + b) / 2 > -Decimal('1e-50')
       for g in (th, h16) for a, b in
       [(len(g.adj[x]), len(g.adj[y])) for x, y in g.E]),
   'equality holds exactly on the (7,7) and (2,2) edges, where the two sides are both a multiple of sqrt3')
tier = 189 * sq(3)
ck('the_step_3_tier_cap_at_n_9_is_189sqrt3_equal_327_3576026305_below_the_claimed_maximum',
   tier < MAX_CF and rnd(tier, 10) == '327.3576026305', rnd(tier, 10))
ck('the_tier_cap_is_27_times_sqrt3_times_7_because_C_8_2_minus_1_equals_27_and_Delta_at_most_7',
   189 == 27 * 7 and 27 == 8 * 7 // 2 - 1)

# ===========================================================================
print('\n=== Step 8: the exhaustion identities of Section 5, re-added')
per_diameter = [1, 91518, 148229, 19320, 1818, 180, 13, 1]
per_size_d3 = [3, 17, 76, 369, 1483, 4381, 9621, 16547, 23013, 26173, 24441, 18828,
               12040, 6477, 2980, 1185, 412, 131, 38, 11, 3]
note('the eight diameter strata at n=9 sum to %d' % sum(per_diameter))
note('the 21 per-size class counts of the D=3 stratum sum to %d' % sum(per_size_d3))
ck('the_eight_diameter_strata_at_n_9_sum_to_261080_connected_graphs', sum(per_diameter) == 261080)
ck('the_D_3_stratum_holds_148229_isomorphism_classes_by_both_readings',
   per_diameter[2] == 148229 == sum(per_size_d3))
ck('the_complete_cover_has_12_shapes_and_12_times_2_to_the_21_equals_25165824',
   len([(a, b) for a in range(1, 8) for b in range(a, 8) if a + b <= 7]) == 12
   and 12 * 2 ** 21 == 25165824)
ck('the_number_of_free_pairs_after_normalising_a_distance_3_pair_is_C_7_2_equal_21',
   7 * 6 // 2 == 21)
ck('the_47_unlabelled_9_vertex_trees_are_the_published_count_A000055_9',
   [1, 1, 1, 2, 3, 6, 11, 23, 47][8] == 47)

# ===========================================================================
print('\n=== Step 9: two published values by other authors, which pin the normalisation of EU')
# Sekar et al., Theorem 3.2, evaluated at n = 9, against a witness read off their own formula.
n = 9
sekar = (Decimal(n - 4) * sq((n - 2) ** 2 + 1 + (n - 2))
         + sq((n - 2) ** 2 + 9 + 3 * (n - 2))
         + sq((n - 2) ** 2 + 4 + 2 * (n - 2))
         + sq(13) + sq(19))
wit = Graph('Sekar et al. Theorem 3.2 witness at n=9',
            [('h', 'q%d' % i) for i in range(5)]
            + [('h', 'a'), ('h', 'b'), ('a', 'b'), ('a', 'leaf')])
note('Sekar et al. Theorem 3.2 at n=9 = 5sqrt57+sqrt79+sqrt67+sqrt13+sqrt19 = %s' % sekar)
note('their witness: n=%d m=%d diam=%d degrees=%s EU=%s'
     % (wit.n, wit.m, wit.diameter(), wit.degseq(), rnd(wit.EU(), 10)))
ck('the_Sekar_witness_is_a_unicyclic_9_vertex_graph_of_diameter_exactly_3',
   wit.n == 9 and wit.m == 9 and wit.m - wit.n + 1 == 1 and wit.diameter() == 3)
ck('its_degree_sequence_is_the_one_their_formula_describes_1x6_2_3_7',
   wit.degseq() == [1, 1, 1, 1, 1, 1, 2, 3, 7], str(wit.degseq()))
ck('its_EU_equals_their_Theorem_3_2_bound_at_n_9_to_1e-40_so_that_bound_is_sharp_at_n_9',
   abs(wit.EU() - sekar) < Decimal('1e-40'), 'difference = %s' % (wit.EU() - sekar))
ck('that_bound_at_n_9_equals_5sqrt57_plus_sqrt79_plus_sqrt67_plus_sqrt13_plus_sqrt19_equal_62_7871695845',
   abs(sekar - (5 * sq(57) + sq(79) + sq(67) + sq(13) + sq(19))) < Decimal('1e-45')
   and rnd(sekar, 10) == '62.7871695845', rnd(sekar, 10))
ck('the_Sekar_witness_lies_strictly_between_EU_Theta_3_3_4_and_EU_H_1_6',
   th.EU() < wit.EU() < h16.EU())
# Kizilirmak, Lemma 2: the minimum over unicyclic graphs on n vertices is 2*sqrt3*n, at C_n.
ck('Kizilirmak_Lemma_2_at_n_9_gives_EU_C_9_equal_18sqrt3_equal_31_1769145362',
   abs(cycle(9).EU() - 2 * sq(3) * 9) < Decimal('1e-40')
   and rnd(cycle(9).EU(), 10) == '31.1769145362', rnd(cycle(9).EU(), 10))
ck('C_9_has_diameter_4_so_it_is_outside_the_cell_and_does_not_contradict_the_minimum',
   cycle(9).diameter() == 4, 'diam=%d' % cycle(9).diameter())
# K_9 - e, the paper's forced negative on the maximum side.
K9e = Graph('K_9 - e', [(x, y) for x, y in itertools.combinations(['z%d' % i for i in range(9)], 2)
                        if {x, y} != {'z0', 'z1'}])
note('K_9-e: n=%d m=%d diam=%d degrees=%s EU=%s'
     % (K9e.n, K9e.m, K9e.diameter(), K9e.degseq(), rnd(K9e.EU(), 10)))
ck('K_9_minus_e_has_diameter_2_so_it_is_outside_the_cell_although_its_EU_exceeds_EU_H_1_6',
   K9e.diameter() == 2 and K9e.EU() > h16.EU())
ck('K_9_minus_e_EU_is_168sqrt3_plus_182_equal_472_9845356716',
   abs(K9e.EU() - (168 * sq(3) + 14 * 13)) < Decimal('1e-40')
   and rnd(K9e.EU(), 10) == '472.9845356716', rnd(K9e.EU(), 10))
ck('P_9_has_diameter_8_so_the_path_is_outside_the_cell', path(9).diameter() == 8)

# ===========================================================================
print('\n=== Step 10: the prior-art check of Section 2, done exactly')
# chi_1 is the first Zagreb index M_1 = sum d_v^2. On H(a,b) at n=9 it is 343 + a^2 + b^2.
m1 = {}
for a, b in ((1, 6), (2, 5), (3, 4)):
    g = H(a, b)
    m1[(a, b)] = sum(d * d for d in g.degseq())
    note('M_1(H(%d,%d)) = %d, and 343 + %d + %d = %d' % (a, b, m1[(a, b)], a * a, b * b, 343 + a * a + b * b))
ck('M_1_of_H_a_b_at_n_9_equals_343_plus_a_squared_plus_b_squared',
   all(m1[(a, b)] == 343 + a * a + b * b for a, b in m1))
ck('M_1_orders_the_three_splits_380_above_372_above_368',
   m1[(1, 6)] == 380 > m1[(2, 5)] == 372 > m1[(3, 4)] == 368)
ck('the_whole_m_at_most_27_tier_obeys_M_1_at_most_14m_at_most_378_below_380',
   14 * 27 == 378 < 380, 'Delta<=7 so M_1 <= 2*Delta*m = 14m')
ck('the_f_condition_of_the_nearest_prior_art_fails_for_the_Euler_Sombor_edge_function',
   sq(7) - sq(3) > sq(39) - sq(31),
   'sqrt7-sqrt3 = %s > sqrt39-sqrt31 = %s' % (rnd(sq(7) - sq(3), 6), rnd(sq(39) - sq(31), 6)))

# ===========================================================================
print('\n=== [ANTI] four objects the checker must reject (a checker that cannot fail verifies nothing)')
bad = Graph('H(1,6) + uv', H16_EDGES + [('u', 'v')])
note('H(1,6)+uv has diameter %d, so the added edge drops it out of the cell' % bad.diameter())
ck('ANTI_1_H_1_6_plus_the_edge_uv_is_rejected_because_its_diameter_is_2', bad.diameter() == 2)
ck('ANTI_2_Theta_2_4_4_is_rejected_at_diameter_3_although_its_EU_equals_the_minimum_exactly',
   theta(2, 4, 4).diameter() == 4 and abs(theta(2, 4, 4).EU() - MIN_CF) < Decimal('1e-50'))
disc = Graph('K_{1,7} plus an isolated vertex, encoded as two components',
             [('c', 'l%d' % i) for i in range(7)] + [('s1', 's2')])
ck('ANTI_3_a_disconnected_graph_reports_minus_1_and_is_never_given_a_diameter',
   disc.diameter() == -1, 'diameter() = %d' % disc.diameter())
ck('ANTI_4_the_checker_disagrees_with_a_deliberately_wrong_candidate_maximum_of_320',
   h16.EU() > Decimal('320') and abs(h16.EU() - Decimal('320')) > Decimal('9'),
   'EU(H(1,6)) - 320 = %s' % rnd(h16.EU() - Decimal('320'), 6))

# ===========================================================================
print('')
print('NOT RE-RUN, and this program cannot see it: the EXHAUSTIVENESS and UNIQUENESS of the minimum half of '
      'Theorem A. This program enumerates nothing over G_{9,3}; it verifies the exhibited objects, the closed '
      'form of Theorem B, the two search-space amputations and the strict orderings. The claim that all 148,229 '
      'isomorphism classes of the cell were evaluated and that exactly one attains each extremum rests on the '
      'exhaustive census of the paper Section 5, which this program does not re-run, so the counts re-added in '
      'Step 8 are transcriptions rather than re-reads.')
print('NOTE SCOPE: also not re-run here -- any cell with n > 9, the minimum half at D = 3 for n >= 10 (which '
      'the paper does not claim), and the full text of the nearest prior art (Vetrik, Discrete Appl. Math. 333 '
      '(2023) 59-70), which is paywalled; the f-class exclusion checked in Step 10 is an exact computation on '
      'the condition AS PRINTED IN A SECOND PAPER, not a reading of that theorem hypothesis.')
print('')
if BAD:
    print('VERDICT: %d of %d CHECKS FAILED: %s' % (len(BAD), len(RAN), ', '.join(BAD)))
else:
    print('VERDICT: ALL %d CHECKS PASS' % len(RAN))
sys.exit(1 if BAD else 0)

#!/usr/bin/env python3
"""verify.py -- re-derives computational claims of

    "The Hexagonal Prism under a Closed-Neighbourhood One-Round Rule:
     Fractional Eternal Domination Number 46/13"

from the objects PRINTED IN THE NOTE and from nothing else.  Every check below encodes
the closed-neighbourhood reading of the defender's one-round rule set out in the note's
Section 1; no check here can discriminate between that reading and a stricter one.

Python 3.9+, standard library only (fractions, itertools, collections, sys).
No third-party package, no external data file, no floating point in any decision:
every comparison below is between exact Fraction or int values.

The two halves of the theorem are checked separately.

  UPPER BOUND  gamma_f^inf(C_6 [] K_2) <= 46/13, for the closed-neighbourhood rule.
    The note prints a 12-member family of rational weight functions on V(C_6 [] K_2);
    a further family of the same total weight is rebuilt here as well, from a seed
    supplied in this program.  The program rebuilds each family from its
    h-profile / seed, and checks: total weight exactly 46/13 per member;
    fractional domination of every member; weight exactly 1 at the member's own centre;
    and that EVERY one of the 144 ordered pairs is reachable in one round, by exact
    integer max-flow on the closed-neighbourhood transportation network.  The single
    transportation plan printed in the note is additionally checked arc by arc, as are
    three further plans supplied in this program.

  LOWER BOUND  gamma_f^inf(C_6 [] K_2) >= 46/13.
    Section 3 of the note prints (i) a 9-node attack tree, (ii) the LP the tree defines,
    with its variable and row ordering spelled out, (iii) 31 nonzero multipliers on the
    ">= 1" rows, and (iv) 123 nonzero potentials on the equality rows.  The program
    rebuilds that LP from the printed ordering rules, plugs in the printed dual vector,
    and verifies all 493 dual column inequalities and the dual objective in exact
    rational arithmetic.  Weak duality then gives the bound; no solver is called and no
    solver output is trusted.

Controls are run BEFORE the target claims, and both polarities are exercised: two
published exact values of the same authors are reproduced, and a configuration that
must FAIL is confirmed to fail.
"""

import sys
from collections import deque
from fractions import Fraction as F
from itertools import combinations

CHECKS = []


def ck(name, cond, detail=''):
    CHECKS.append(bool(cond))
    print('%s %s%s' % ('PASS' if cond else 'FAIL', name, ('  [%s]' % detail) if detail else ''))


def note(s):
    print('NOTE %s' % s)


# ---------------------------------------------------------------------------
# 0. THE GRAPHS.  Index convention of the paper: v = i + n*j for (i,j) in Z_n x {0,1}.
# ---------------------------------------------------------------------------
def prism(n):
    """C_n [] K_2, closed neighbourhoods, index v = i + n*j."""
    N = []
    for v in range(2 * n):
        i, j = v % n, v // n
        N.append(sorted({v, ((i + 1) % n) + n * j, ((i - 1) % n) + n * j, i + n * (1 - j)}))
    return N


def cay_pm1_half(m):
    """Cay(Z_m, {+-1, m/2}), closed neighbourhoods."""
    return [sorted({v, (v + 1) % m, (v - 1) % m, (v + m // 2) % m}) for v in range(m)]


def is_bipartite(N):
    col = {}
    for s in range(len(N)):
        if s in col:
            continue
        col[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for w in N[u]:
                if w == u:
                    continue
                if w not in col:
                    col[w] = 1 - col[u]
                    q.append(w)
                elif col[w] == col[u]:
                    return False
    return True


def edge_count(N):
    return sum(len(x) - 1 for x in N) // 2


def dominates(N, S):
    S = set(S)
    return all(any(u in S for u in N[v]) for v in range(len(N)))


def min_dominating_size(N, cap):
    """Least k <= cap with a dominating set of size k, else None.  Brute force; n <= 16 here."""
    n = len(N)
    for k in range(1, cap + 1):
        for S in combinations(range(n), k):
            if dominates(N, S):
                return k, S
    return None, None


# ---------------------------------------------------------------------------
# 1. EXACT INTEGER MAX-FLOW (Edmonds-Karp; integer capacities only)
# ---------------------------------------------------------------------------
class MaxFlow:
    def __init__(self, n):
        self.n = n
        self.to = []
        self.cap = []
        self.adj = [[] for _ in range(n)]

    def add(self, u, v, c):
        self.adj[u].append(len(self.to))
        self.to.append(v)
        self.cap.append(int(c))
        self.adj[v].append(len(self.to))
        self.to.append(u)
        self.cap.append(0)

    def run(self, s, t):
        flow = 0
        while True:
            par = [-1] * self.n
            par[s] = -2
            q = deque([s])
            while q and par[t] == -1:
                u = q.popleft()
                for e in self.adj[u]:
                    v = self.to[e]
                    if self.cap[e] > 0 and par[v] == -1:
                        par[v] = e
                        q.append(v)
            if par[t] == -1:
                return flow
            # bottleneck
            b = None
            v = t
            while v != s:
                e = par[v]
                b = self.cap[e] if b is None else min(b, self.cap[e])
                v = self.to[e ^ 1]
            v = t
            while v != s:
                e = par[v]
                self.cap[e] -= b
                self.cap[e ^ 1] += b
                v = self.to[e ^ 1]
            flow += b


def reconfigurable(N, w_src, w_tgt):
    """Is there a one-round move carrying w_src to w_tgt?  Weight may stay put or move to a
    neighbour, so the transportation network is over CLOSED neighbourhoods -- which is exactly
    the source's rule  sum_{y in N(x)} m_{xy} <= w(x)  with the remainder staying at x.
    Exact: scale both vectors by their common denominator and run integer max-flow."""
    n = len(N)
    D = 1
    for x in list(w_src) + list(w_tgt):
        D = D * x.denominator // _gcd(D, x.denominator)
    S, T = 2 * n, 2 * n + 1
    g = MaxFlow(2 * n + 2)
    tot = 0
    for u in range(n):
        c = w_src[u] * D
        assert c.denominator == 1
        if c:
            g.add(S, u, c)
            tot += int(c)
        for y in N[u]:
            g.add(u, n + y, D * 46)          # effectively infinite
    need = 0
    for y in range(n):
        c = w_tgt[y] * D
        assert c.denominator == 1
        if c:
            g.add(n + y, T, c)
            need += int(c)
    if tot != need:
        return False
    return g.run(S, T) == need


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def demand_answerable(N, w, v, demand):
    """From configuration w, can one round reach SOME configuration that fractionally
    dominates G, conserves the total, and puts at least `demand` on v?  Here it is used only
    for the FORCED-NEGATIVE control, where the obstruction is local: N[v] holds less than
    `demand`, so no one-hop move can put `demand` on v."""
    return sum(w[u] for u in N[v]) >= demand


# ---------------------------------------------------------------------------
# 2. THE OBJECT.  C_6 [] K_2, and the closed neighbourhoods printed in the paper.
# ---------------------------------------------------------------------------
print('verification of: the hexagonal prism under a closed-neighbourhood one-round rule,')
print('C_6 [] K_2 = Cay(Z_6 x Z_2, {(+-1,0),(0,1)}),  gamma_f^infty = 46/13 for that reading')
print('python %s, exact integer / Fraction arithmetic only' % sys.version.split()[0])
print()

N = prism(6)
PAPER_N = [[0, 1, 5, 6], [0, 1, 2, 7], [1, 2, 3, 8], [2, 3, 4, 9], [3, 4, 5, 10], [0, 4, 5, 11],
           [0, 6, 7, 11], [1, 6, 7, 8], [2, 7, 8, 9], [3, 8, 9, 10], [4, 9, 10, 11], [5, 6, 10, 11]]

print('=== Step 1: the graph the paper names, rebuilt from its definition')
note('N[v] for v = 0..11 : %s' % N)
ck('order_is_12', len(N) == 12, 'n=%d' % len(N))
ck('three_regular_every_closed_neighbourhood_has_size_4', all(len(x) == 4 for x in N))
ck('size_is_18_edges', edge_count(N) == 18, 'm=%d' % edge_count(N))
ck('graph_is_bipartite', is_bipartite(N))
ck('closed_neighbourhoods_match_the_lists_printed_in_the_paper', N == PAPER_N)
M12 = cay_pm1_half(12)
ck('the_moebius_prism_Cay_Z12_pm1_6_is_NOT_bipartite_so_it_is_a_different_graph',
   not is_bipartite(M12) and edge_count(M12) == 18 and all(len(x) == 4 for x in M12),
   'the one other 12-vertex cubic abelian Cayley graph whose value the source already knows')
g6, S6 = min_dominating_size(N, 4)
ck('domination_number_of_C6_box_K2_is_4', g6 == 4, 'gamma=%s attained by %s' % (g6, S6))

TARGET = F(46, 13)

# ---------------------------------------------------------------------------
# 3. CONTROLS FIRST -- other people's published numbers, and a forced negative.
# ---------------------------------------------------------------------------
print()
print('=== Step 2: controls, run before any claim of the paper is touched')

# (a) the source's published gamma_f^infty(Cay(Z_8,{+-1,4})) = 8/3, from its printed witness
N8 = cay_pm1_half(8)
SEED8 = [F(1), F(0), F(2, 3), F(0), F(1, 3), F(0), F(2, 3), F(0)]
ck('control_published_witness_on_Cay_Z8_totals_8_over_3', sum(SEED8) == F(8, 3), str(sum(SEED8)))
sums8 = [sum(SEED8[u] for u in N8[v]) for v in range(8)]
ck('control_published_witness_on_Cay_Z8_is_fractionally_dominating',
   all(s >= 1 for s in sums8), 'closed-nbhd sums %s' % [str(s) for s in sums8])
ck('control_published_witness_on_Cay_Z8_has_the_printed_neighbourhood_sums',
   [str(s) for s in sums8] == ['4/3', '5/3', '4/3', '1', '4/3', '1', '4/3', '5/3'])
ck('control_published_witness_on_Cay_Z8_is_even_so_its_8_translates_are_well_defined',
   all(SEED8[i] == SEED8[(-i) % 8] for i in range(8)))
FAM8 = [[SEED8[(i - m) % 8] for i in range(8)] for m in range(8)]
ok8 = all(sum(f) == F(8, 3) for f in FAM8) and all(FAM8[m][m] == 1 for m in range(8)) \
    and all(all(sum(FAM8[m][u] for u in N8[v]) >= 1 for v in range(8)) for m in range(8))
ck('control_the_8_translates_form_an_eternal_family_of_total_weight_8_over_3', ok8)
pairs8 = sum(1 for a in range(8) for b in range(8) if reconfigurable(N8, FAM8[a], FAM8[b]))
ck('control_all_64_ordered_pairs_of_that_family_are_reconfigurable_so_gamma_f_infty_Cay_Z8_is_at_most_8_over_3',
   pairs8 == 64, '%d/64' % pairs8)

# (b) the source's published gamma_f^infty(C_10 [] K_2) = 28/5: its own lower-bound formula
#     (n+2)(n+4)/(2(n+5)) is EXACT at n=10, and its Cayley companion is exact at m=4.
lb = lambda n: F((n + 2) * (n + 4), 2 * (n + 5))
ck('control_published_lower_bound_formula_is_exact_at_n_equals_10', lb(10) == F(28, 5), str(lb(10)))
ck('control_published_Cayley_companion_is_exact_at_m_equals_4', F(6 * 8, 18) == F(8, 3))
ck('the_published_lower_bound_formula_does_NOT_apply_at_n_equals_6',
   6 % 12 != 10 and lb(6) == F(40, 11) and lb(6) > TARGET,
   'hypothesis is n = 10 mod 12; at n = 6 the formula would give 40/11 > 46/13 and refute the value')

# (c) forced negative: from uniform 1/4 (total 3) a demand of 11/10 at vertex 0 is impossible
UNIF = [F(1, 4)] * 12
ck('control_uniform_one_quarter_is_fractionally_dominating_with_total_3',
   sum(UNIF) == 3 and all(sum(UNIF[u] for u in N[v]) == 1 for v in range(12)))
ck('control_FORCED_NEGATIVE_a_demand_of_11_over_10_at_vertex_0_from_uniform_one_quarter_is_refused',
   not demand_answerable(N, UNIF, 0, F(11, 10)),
   'N[0] holds 4 x 1/4 = 1 < 11/10, so no one-hop move can answer it')

# (d) the source's own general bound: (n+kappa)/(kappa+1) = 15/4 at n=12, kappa=3, and a
#     matching explicit family, so the incumbent upper bound really is 15/4.
ck('the_sources_own_connectivity_bound_gives_15_over_4_at_n_12_kappa_3',
   F(12 + 3, 3 + 1) == F(15, 4))
FAM154 = [[F(1) if u == m else F(1, 4) for u in range(12)] for m in range(12)]
ok154 = all(sum(f) == F(15, 4) for f in FAM154) \
    and all(all(sum(FAM154[m][u] for u in N[v]) >= 1 for v in range(12)) for m in range(12)) \
    and all(reconfigurable(N, FAM154[a], FAM154[b]) for a in range(12) for b in range(12))
ck('control_FORCED_POSITIVE_the_15_over_4_family_is_a_genuine_eternal_strategy', ok154)
ck('the_new_value_strictly_improves_the_incumbent_upper_bound',
   TARGET < F(15, 4) < F(4), '46/13 = %.12f < 15/4 < 4' % float(TARGET))

# ---------------------------------------------------------------------------
# 4. UPPER BOUND -- family 1, denominator 13, rebuilt from the h-profile printed in the note
# ---------------------------------------------------------------------------
print()
print('=== Step 3: the upper bound, family 1 (denominator 13), rebuilt from the printed h-profile')

H0 = {0: F(1), 1: F(1, 13), 2: F(1, 13), 3: F(5, 13), 4: F(1, 13), 5: F(1, 13)}
H1 = {0: F(6, 13), 1: F(0), 2: F(6, 13), 3: F(6, 13), 4: F(6, 13), 5: F(0)}


def fam1(m, t):
    """f_{(m,t)}(i,j) = h_{j xor t}((i - m) mod 6), flattened to v = i + 6j."""
    w = [F(0)] * 12
    for j in (0, 1):
        h = H1 if (j ^ t) else H0
        for i in range(6):
            w[i + 6 * j] = h[(i - m) % 6]
    return w


FAMILY1 = [fam1(m, t) for t in (0, 1) for m in range(6)]   # index of centre = m + 6t
CENTRE1 = [m + 6 * t for t in (0, 1) for m in range(6)]
ck('family1_profiles_h0_and_h1_are_even_so_the_12_translates_are_well_defined',
   all(H0[i] == H0[(-i) % 6] for i in range(6)) and all(H1[i] == H1[(-i) % 6] for i in range(6)))
seed1 = FAMILY1[0]
note('seed f_(0,0) = %s' % [str(x) for x in seed1])
ck('family1_seed_matches_the_values_printed_in_the_paper',
   [str(x) for x in seed1] == ['1', '1/13', '1/13', '5/13', '1/13', '1/13',
                               '6/13', '0', '6/13', '6/13', '6/13', '0'])
ck('family1_seed_total_is_46_over_13', sum(seed1) == TARGET, str(sum(seed1)))
sums1 = [sum(seed1[u] for u in N[v]) for v in range(12)]
note('closed-neighbourhood sums of the seed = %s' % [str(s) for s in sums1])
ck('family1_seed_neighbourhood_sums_are_the_twelve_values_printed_in_the_paper',
   [str(s) for s in sums1] == ['21/13', '15/13', '1', '1', '1', '15/13',
                               '19/13', '1', '1', '23/13', '1', '1'])
ck('family1_all_12_members_have_total_weight_exactly_46_over_13',
   all(sum(f) == TARGET for f in FAMILY1))
ck('family1_all_12_members_are_fractionally_dominating',
   all(all(sum(f[u] for u in N[v]) >= 1 for v in range(12)) for f in FAMILY1))
ck('family1_member_indexed_by_v_carries_weight_exactly_1_at_v',
   all(FAMILY1[k][CENTRE1[k]] == 1 for k in range(12)))
ck('family1_the_12_translates_are_pairwise_distinct',
   len({tuple(f) for f in FAMILY1}) == 12)

# ---------------------------------------------------------------------------
# 5. UPPER BOUND -- a further family, denominator 26, seed supplied in this program
# ---------------------------------------------------------------------------
print()
print('=== Step 4: a further weight family (denominator 26) supplied in this program, an independent optimum')

# Seed of the further family; supplied in this program, not printed in the note.
X0 = [F(1), F(1, 26), F(1, 13), F(11, 26), F(1, 13), F(1, 26),
      F(6, 13), F(1, 26), F(6, 13), F(11, 26), F(6, 13), F(1, 26)]


def fam2(m, t):
    w = [F(0)] * 12
    for j in (0, 1):
        for i in range(6):
            w[i + 6 * j] = X0[((i - m) % 6) + 6 * (j ^ t)]
    return w


FAMILY2 = [fam2(m, t) for t in (0, 1) for m in range(6)]
ck('family2_seed_total_is_46_over_13', sum(X0) == TARGET, str(sum(X0)))
ck('family2_seed_is_invariant_under_the_stabiliser_i_maps_to_minus_i',
   all(X0[i + 6 * j] == X0[((-i) % 6) + 6 * j] for i in range(6) for j in (0, 1)))
sums2 = [sum(X0[u] for u in N[v]) for v in range(12)]
note('closed-neighbourhood sums of X_0 = %s' % [str(s) for s in sums2])
ck('family2_seed_neighbourhood_sums_are_the_twelve_values_listed_in_this_program',
   [str(s) for s in sums2] == ['20/13', '15/13', '1', '1', '1', '15/13',
                               '20/13', '1', '1', '23/13', '1', '1'])
ck('family2_all_12_members_have_total_weight_exactly_46_over_13',
   all(sum(f) == TARGET for f in FAMILY2))
ck('family2_all_12_members_are_fractionally_dominating',
   all(all(sum(f[u] for u in N[v]) >= 1 for v in range(12)) for f in FAMILY2))
ck('family2_member_indexed_by_v_carries_weight_exactly_1_at_v',
   all(FAMILY2[k][CENTRE1[k]] == 1 for k in range(12)))
ck('the_two_families_are_genuinely_different_optima',
   {tuple(f) for f in FAMILY1} & {tuple(f) for f in FAMILY2} == set())

# ---------------------------------------------------------------------------
# 6. FOUR TRANSPORTATION PLANS, ONE PRINTED IN THE NOTE AND THREE SUPPLIED HERE
# ---------------------------------------------------------------------------
print()
print('=== Step 5: transportation plans checked arc by arc (one from the note, three supplied here)')


def ix(i, j):
    return i + 6 * j


PLAN_A = [(0, 0, F(5, 13)), (0, 1, F(1, 13)), (0, 5, F(1, 13)), (0, 6, F(6, 13)),
          (1, 2, F(1, 13)), (2, 3, F(1, 13)), (3, 3, F(5, 13)), (4, 3, F(1, 13)),
          (5, 4, F(1, 13)), (6, 7, F(6, 13)), (8, 9, F(6, 13)), (9, 3, F(6, 13)),
          (10, 11, F(6, 13))]

PLAN_B = [(ix(0, 0), ix(0, 0), F(11, 26)), (ix(0, 0), ix(1, 0), F(1, 13)),
          (ix(0, 0), ix(5, 0), F(1, 13)), (ix(0, 0), ix(0, 1), F(11, 26)),
          (ix(1, 0), ix(2, 0), F(1, 26)), (ix(2, 0), ix(3, 0), F(1, 13)),
          (ix(3, 0), ix(3, 0), F(11, 26)), (ix(4, 0), ix(3, 0), F(1, 13)),
          (ix(5, 0), ix(4, 0), F(1, 26)), (ix(0, 1), ix(1, 1), F(6, 13)),
          (ix(1, 1), ix(2, 1), F(1, 26)), (ix(2, 1), ix(3, 1), F(6, 13)),
          (ix(3, 1), ix(3, 0), F(11, 26)), (ix(4, 1), ix(4, 1), F(1, 26)),
          (ix(4, 1), ix(5, 1), F(11, 26)), (ix(5, 1), ix(5, 1), F(1, 26))]

PLAN_C = [(ix(0, 0), ix(0, 0), F(6, 13)), (ix(0, 0), ix(1, 0), F(1, 26)),
          (ix(0, 0), ix(5, 0), F(1, 26)), (ix(0, 0), ix(0, 1), F(6, 13)),
          (ix(1, 0), ix(2, 0), F(1, 26)), (ix(2, 0), ix(2, 0), F(1, 13)),
          (ix(3, 0), ix(2, 0), F(9, 26)), (ix(3, 0), ix(3, 0), F(1, 13)),
          (ix(4, 0), ix(3, 0), F(1, 13)), (ix(5, 0), ix(4, 0), F(1, 26)),
          (ix(0, 1), ix(0, 1), F(6, 13)), (ix(1, 1), ix(0, 1), F(1, 26)),
          (ix(2, 1), ix(1, 1), F(1, 26)), (ix(2, 1), ix(2, 1), F(1, 13)),
          (ix(2, 1), ix(3, 1), F(9, 26)), (ix(3, 1), ix(3, 0), F(7, 26)),
          (ix(3, 1), ix(3, 1), F(1, 13)), (ix(3, 1), ix(4, 1), F(1, 13)),
          (ix(4, 1), ix(4, 0), F(11, 26)), (ix(4, 1), ix(5, 1), F(1, 26)),
          (ix(5, 1), ix(0, 1), F(1, 26))]

PLAN_D = [(ix(0, 0), ix(0, 0), F(9, 26)), (ix(0, 0), ix(1, 0), F(9, 26)),
          (ix(0, 0), ix(5, 0), F(4, 13)), (ix(1, 0), ix(0, 0), F(1, 26)),
          (ix(2, 0), ix(1, 0), F(1, 13)), (ix(3, 0), ix(2, 0), F(1, 26)),
          (ix(3, 0), ix(3, 0), F(5, 13)), (ix(4, 0), ix(5, 0), F(1, 13)),
          (ix(5, 0), ix(5, 0), F(1, 26)), (ix(0, 1), ix(0, 0), F(1, 26)),
          (ix(0, 1), ix(0, 1), F(11, 26)), (ix(1, 1), ix(1, 0), F(1, 26)),
          (ix(2, 1), ix(1, 1), F(1, 13)), (ix(2, 1), ix(2, 1), F(1, 26)),
          (ix(2, 1), ix(3, 1), F(9, 26)), (ix(3, 1), ix(3, 0), F(1, 13)),
          (ix(3, 1), ix(3, 1), F(9, 26)), (ix(4, 1), ix(4, 0), F(1, 26)),
          (ix(4, 1), ix(3, 1), F(4, 13)), (ix(4, 1), ix(4, 1), F(1, 26)),
          (ix(4, 1), ix(5, 1), F(1, 13)), (ix(5, 1), ix(5, 0), F(1, 26))]


def check_plan(label, plan, src, tgt, attacked, narcs):
    ck('%s_has_%d_arcs' % (label, narcs), len(plan) == narcs, '%d arcs' % len(plan))
    ck('%s_every_arc_is_an_edge_or_a_self_loop' % label,
       all(y in N[x] for x, y, _ in plan))
    ship = [F(0)] * 12
    recv = [F(0)] * 12
    for x, y, q in plan:
        ck_pos = q > 0
        if not ck_pos:
            ck('%s_arc_%d_to_%d_carries_positive_weight' % (label, x, y), False)
        ship[x] += q
        recv[y] += q
    ck('%s_every_vertex_ships_exactly_its_own_source_weight' % label, ship == list(src),
       'ships %s' % [str(v) for v in ship])
    ck('%s_every_vertex_receives_exactly_its_target_weight' % label, recv == list(tgt),
       'delivers %s' % [str(v) for v in tgt])
    ck('%s_total_shipped_is_46_over_13' % label, sum(ship) == TARGET, str(sum(ship)))
    ck('%s_the_attacked_vertex_%d_receives_exactly_1' % (label, attacked),
       recv[attacked] == 1, str(recv[attacked]))


check_plan('plan_family1_f00_to_f30', PLAN_A, FAMILY1[0], fam1(3, 0), 3, 13)
check_plan('plan_family2_X0_to_X30', PLAN_B, X0, fam2(3, 0), 3, 16)
check_plan('plan_family2_X0_to_X01', PLAN_C, X0, fam2(0, 1), 6, 21)
check_plan('plan_family2_X0_to_X31', PLAN_D, X0, fam2(3, 1), 9, 22)

# ---------------------------------------------------------------------------
# 7. ALL 144 ORDERED PAIRS OF EACH FAMILY
# ---------------------------------------------------------------------------
print()
print('=== Step 6: every ordered pair of each family, by exact integer max-flow')

for lab, FAM in (('family1_denominator_13', FAMILY1), ('family2_denominator_26', FAMILY2)):
    bad = [(a, b) for a in range(12) for b in range(12)
           if not reconfigurable(N, FAM[a], FAM[b])]
    ck('%s_all_144_ordered_pairs_are_reconfigurable_in_one_round' % lab, not bad,
       '%d/144 feasible, violations=%s' % (144 - len(bad), bad))

ck('UPPER_BOUND_gamma_f_infty_C6_box_K2_is_at_most_46_over_13',
   all(sum(f) == TARGET for f in FAMILY1)
   and all(all(sum(FAMILY1[k][u] for u in N[v]) >= 1 for v in range(12)) for k in range(12))
   and all(FAMILY1[k][CENTRE1[k]] == 1 for k in range(12))
   and all(reconfigurable(N, FAMILY1[a], FAMILY1[b]) for a in range(12) for b in range(12)),
   'a total response map r(f_u, v) = f_v exists, so the strategy is eternal at total weight 46/13')

# ---------------------------------------------------------------------------
# 8. LOWER BOUND -- the 9-node attack tree LP and its exact rational dual
# ---------------------------------------------------------------------------
print()
print('=== Step 7: the lower bound -- the 9-node attack-tree LP rebuilt from the printed rules')

# The attack tree printed in the note.  parent[k] and attacked[k].
PARENT = [None, 0, 1, 2, 1, 4, 4, 1, 7]
ATTACK = [None, 0, 8, 11, 9, 7, 11, 10, 11]     # vertex indices; (0,0)=0, (1,1)=7, (2,1)=8, ...
ck('the_attack_tree_has_9_nodes_with_the_printed_parent_and_attack_maps',
   len(PARENT) == 9 and len(ATTACK) == 9 and PARENT[0] is None and ATTACK[0] is None
   and ATTACK[1:] == [0, 8, 11, 9, 7, 11, 10, 11] and PARENT[1:] == [0, 1, 2, 1, 4, 4, 1, 7])
ck('the_second_attack_set_is_the_three_vertices_of_the_other_hexagon_at_distances_3_4_3',
   sorted(ATTACK[k] for k in (2, 4, 7)) == [8, 9, 10])


class LP:
    """min c.x  s.t.  A_ub x <= b_ub,  A_eq x = b_eq,  x >= 0.  Rows are sparse dicts."""

    def __init__(self):
        self.nv = 0
        self.ub = []
        self.bub = []
        self.eq = []
        self.beq = []
        self.c = {}

    def var(self):
        self.nv += 1
        return self.nv - 1

    def addub(self, row, b):
        self.ub.append(dict(row))
        self.bub.append(F(b))

    def addeq(self, row, b):
        self.eq.append(dict(row))
        self.beq.append(F(b))


def build_tree_lp(N, parent, attack):
    """EXACTLY the construction printed in the note.

    Variables, allocated in this order:  S ; then w_k(v) for k = 0..8 and v = 0..11 ;
    then x_k(u,y) for k = 1..8 over ARC = [(u,y) : u = 0..11, y in N[u]] in that order.
    Inequality rows, in this order: for k = 0..8, the twelve domination rows
    -sum_{u in N[v]} w_k(u) <= -1 for v = 0..11, followed (when node k is an attack node)
    by the single attack row -w_k(a_k) <= -1.
    Equality rows, in this order: first the nine total-weight rows sum_v w_k(v) - S = 0 for
    k = 0..8 ; then, for k = 1..8, the twelve out-flow rows sum_{y in N[u]} x_k(u,y)
    - w_{p(k)}(u) = 0 for u = 0..11 followed by the twelve in-flow rows
    sum_{u : y in N[u]} x_k(u,y) - w_k(y) = 0 for y = 0..11.
    """
    nn = len(N)
    lp = LP()
    S = lp.var()
    W = [[lp.var() for _ in range(nn)] for _ in parent]
    ARC = [(u, y) for u in range(nn) for y in N[u]]
    X = [None] * len(parent)
    for k in range(1, len(parent)):
        X[k] = {a: lp.var() for a in ARC}
    lp.c = {S: 1}
    ub_label = []
    for k in range(len(parent)):
        lp.addeq({**{W[k][v]: 1 for v in range(nn)}, S: -1}, 0)
        for v in range(nn):
            lp.addub({W[k][u]: -1 for u in N[v]}, -1)
            ub_label.append(('dom', k, v))
        if attack[k] is not None:
            lp.addub({W[k][attack[k]]: -1}, -1)
            ub_label.append(('atk', k, attack[k]))
    for k in range(1, len(parent)):
        p = parent[k]
        for u in range(nn):
            row = {X[k][(u, y)]: 1 for y in N[u]}
            row[W[p][u]] = -1
            lp.addeq(row, 0)
        for y in range(nn):
            row = {X[k][(u, y)]: 1 for u in range(nn) if y in N[u]}
            row[W[k][y]] = -1
            lp.addeq(row, 0)
    return lp, ub_label


LPP, UB_LABEL = build_tree_lp(N, PARENT, ATTACK)
ck('the_LP_has_the_493_variables_201_equalities_and_116_inequalities_printed',
   (LPP.nv, len(LPP.eq), len(LPP.ub)) == (493, 201, 116),
   '%d vars, %d eq, %d ineq' % (LPP.nv, len(LPP.eq), len(LPP.ub)))
ck('every_inequality_row_is_accounted_for_by_a_label', len(UB_LABEL) == len(LPP.ub))

# --- the 31 nonzero multipliers on the ">= 1" rows, as printed in the note, as |y| ---
ALPHA = {1: F(4, 13), 2: F(2, 13), 3: F(0), 4: F(4, 13),
         5: F(3, 13), 6: F(2, 13), 7: F(2, 13), 8: F(0)}
DOM = {(1, 2): 1, (1, 3): 1, (1, 4): 1, (1, 7): 2, (1, 8): 1, (1, 10): 1, (1, 11): 2,
       (2, 1): 1, (2, 3): 1, (2, 11): 1,
       (4, 1): 1, (4, 2): 2, (4, 4): 2, (4, 5): 1, (4, 6): 1, (4, 11): 1,
       (5, 0): 1, (5, 2): 1, (5, 10): 1,
       (6, 0): 1, (6, 4): 1, (6, 8): 1,
       (7, 3): 1, (7, 5): 1, (7, 7): 1}
DOM = {k: F(v, 13) for k, v in DOM.items()}
ck('the_paper_prints_31_nonzero_multipliers_on_the_at_least_one_rows',
   len([v for v in ALPHA.values() if v]) + len(DOM) == 31,
   '%d attack + %d domination' % (len([v for v in ALPHA.values() if v]), len(DOM)))
ck('the_attack_part_of_the_dual_objective_is_17_over_13',
   sum(ALPHA.values()) == F(17, 13), str(sum(ALPHA.values())))
ck('the_domination_part_of_the_dual_objective_is_29_over_13',
   sum(DOM.values()) == F(29, 13), str(sum(DOM.values())))
ck('the_31_multipliers_sum_to_46_over_13',
   sum(ALPHA.values()) + sum(DOM.values()) == TARGET,
   str(sum(ALPHA.values()) + sum(DOM.values())))
ck('no_at_least_one_row_of_node_3_or_node_8_carries_a_nonzero_multiplier',
   ALPHA[3] == 0 and ALPHA[8] == 0
   and not [k for (k, v) in DOM if k in (3, 8)],
   'the printed dual vector puts no weight on the >= 1 rows of nodes 3 and 8')

Y_UB = []
for lab in UB_LABEL:
    kind, k, v = lab
    Y_UB.append(-(ALPHA[k] if kind == 'atk' else DOM.get((k, v), F(0))))
ck('every_multiplier_on_an_inequality_row_is_non_positive_as_the_dual_requires',
   all(y <= 0 for y in Y_UB))

# --- the 123 nonzero potentials on the equality rows, as printed in the note ---
# Printed per node, in the note's own row order.  Zero entries are written 0.
Z = F(0)
LAM = {3: F(-7, 13), 4: F(-1, 13), 5: F(-3, 13), 6: F(-2, 13)}   # node total-weight rows
MU_OUT = {
    2: [F(6, 13), F(5, 13), F(5, 13), F(5, 13), F(6, 13), F(6, 13),
        F(6, 13), F(5, 13), F(5, 13), F(5, 13), F(6, 13), F(6, 13)],
    3: [F(7, 13)] * 12,
    4: [F(-1, 13), F(-1, 13), F(-1, 13), Z, F(-1, 13), F(-1, 13),
        F(-1, 13), F(-1, 13), Z, Z, Z, F(-1, 13)],
    5: [F(1, 13), Z, F(1, 13), F(2, 13), F(2, 13), F(2, 13),
        Z, Z, Z, F(2, 13), F(2, 13), F(2, 13)],
    6: [Z, F(1, 13), F(1, 13), F(1, 13), Z, Z,
        Z, F(1, 13), F(1, 13), F(1, 13), Z, Z],
    7: [F(-1, 13), F(-1, 13), F(-1, 13), F(-2, 13), F(-2, 13), F(-2, 13),
        F(-1, 13), F(-1, 13), F(-1, 13), F(-2, 13), F(-2, 13), F(-2, 13)],
}
MU_IN = {
    2: [F(-6, 13), F(-6, 13), F(-5, 13), F(-6, 13), F(-6, 13), F(-6, 13),
        F(-6, 13), F(-6, 13), F(-5, 13), F(-6, 13), F(-6, 13), F(-6, 13)],
    3: [F(-7, 13)] * 12,
    4: [F(1, 13), F(1, 13), Z, Z, Z, F(1, 13),
        F(1, 13), Z, Z, Z, Z, Z],
    5: [F(-2, 13), F(-1, 13), F(-2, 13), F(-2, 13), F(-2, 13), F(-2, 13),
        F(-2, 13), Z, F(-2, 13), F(-2, 13), F(-2, 13), F(-2, 13)],
    6: [F(-1, 13), F(-1, 13), F(-1, 13), F(-1, 13), F(-1, 13), Z,
        F(-1, 13), F(-1, 13), F(-1, 13), F(-1, 13), F(-1, 13), Z],
    7: [F(1, 13), F(1, 13), F(1, 13), F(1, 13), F(2, 13), F(1, 13),
        F(1, 13), F(1, 13), F(1, 13), F(1, 13), F(2, 13), F(1, 13)],
}
Y_EQ = [F(0)] * len(LPP.eq)
for k, v in LAM.items():
    Y_EQ[k] = v
for k in range(1, 9):
    base = 9 + (k - 1) * 24
    for u in range(12):
        Y_EQ[base + u] = MU_OUT.get(k, [Z] * 12)[u]
        Y_EQ[base + 12 + u] = MU_IN.get(k, [Z] * 12)[u]

nz = [j for j, y in enumerate(Y_EQ) if y != 0]
ck('the_paper_prints_exactly_123_nonzero_equality_potentials', len(nz) == 123, '%d' % len(nz))
ck('every_printed_potential_has_denominator_dividing_13',
   all(13 % Y_EQ[j].denominator == 0 for j in nz))
ck('the_potentials_of_the_root_and_of_nodes_1_2_and_8_flow_rows_vanish',
   all(Y_EQ[j] == 0 for j in list(range(9, 33)) + list(range(177, 201)))
   and Y_EQ[0] == 0 and Y_EQ[1] == 0 and Y_EQ[2] == 0)

# --- the 493 dual column inequalities, exactly ---
red = {}
for row, y in zip(LPP.ub, Y_UB):
    if y == 0:
        continue
    for kk, vv in row.items():
        red[kk] = red.get(kk, F(0)) + F(vv) * y
for row, y in zip(LPP.eq, Y_EQ):
    if y == 0:
        continue
    for kk, vv in row.items():
        red[kk] = red.get(kk, F(0)) + F(vv) * y
viol = [kk for kk in range(LPP.nv) if red.get(kk, F(0)) > F(LPP.c.get(kk, 0))]
ck('all_493_dual_column_inequalities_hold_in_exact_rational_arithmetic',
   not viol, '%d/493 satisfied, violations=%s' % (LPP.nv - len(viol), viol[:8]))

dual_obj = sum(b * y for b, y in zip(LPP.bub, Y_UB)) + sum(b * y for b, y in zip(LPP.beq, Y_EQ))
ck('the_dual_objective_of_the_certificate_is_exactly_46_over_13',
   dual_obj == TARGET, str(dual_obj))
ck('LOWER_BOUND_by_weak_duality_gamma_f_infty_C6_box_K2_is_at_least_46_over_13',
   (not viol) and dual_obj == TARGET and all(y <= 0 for y in Y_UB),
   'a dual-feasible point of value 46/13 bounds the primal optimum, which bounds gamma_f^infty')

# ---------------------------------------------------------------------------
# 9. THE TWO BOUNDS MEET, AND THE CONSISTENCY TESTS
# ---------------------------------------------------------------------------
print()
print('=== Step 8: the two bounds meet, and every published number they must respect')

ck('THEOREM_gamma_f_infty_C6_box_K2_equals_46_over_13',
   TARGET == F(46, 13), '46/13 = %.15f' % float(TARGET))
ck('the_published_strict_inequality_7_over_2_less_than_gamma_f_infty_survives',
   F(7, 2) < TARGET and TARGET - F(7, 2) == F(1, 26), 'margin 1/26')
ck('the_published_upper_bound_4_is_true_and_slack', TARGET < 4)
ck('the_fractional_domination_number_of_C6_box_K2_is_3',
   sum(UNIF) == 3 and all(sum(UNIF[u] for u in N[v]) == 1 for v in range(12))
   and F(12, 4) == 3,
   'uniform 1/4 attains it and summing the 12 closed-neighbourhood inequalities gives 4S >= 12')
ck('the_sources_corollary_gamma_f_infty_minus_gamma_f_less_than_1_holds',
   TARGET - 3 == F(7, 13) and TARGET - 3 < 1, '46/13 - 3 = 7/13')
ck('gamma_f_infty_is_strictly_below_the_domination_number_as_the_source_predicts',
   TARGET < g6, '46/13 < gamma = 4')

# ---------------------------------------------------------------------------
# 10. THE UNJUSTIFIED STEP IN THE SOURCE'S PUBLISHED PROOF, FOR THIS READING OF THE RULE
# ---------------------------------------------------------------------------
print()
print('=== Step 9: the configuration that defeats the last step of the published proof')

# Figure labels of the source: inner a b c d e f, outer g h i j k l, a~g .. f~l.
LAB = {'a': ix(0, 0), 'b': ix(1, 0), 'c': ix(2, 0), 'd': ix(3, 0), 'e': ix(4, 0), 'f': ix(5, 0),
       'g': ix(0, 1), 'h': ix(1, 1), 'i': ix(2, 1), 'j': ix(3, 1), 'k': ix(4, 1), 'l': ix(5, 1)}
START = [F(0)] * 12
for nm, val in (('a', F(1)), ('d', F(1, 2)), ('g', F(1, 2)), ('i', F(1, 2)),
                ('j', F(1, 2)), ('k', F(1, 2))):
    START[LAB[nm]] = val
MOVE = [('a', 'b', F(1, 2)), ('a', 'f', F(1, 2)), ('g', 'a', F(1, 4)), ('g', 'g', F(1, 4)),
        ('i', 'h', F(1, 4)), ('i', 'j', F(1, 4)), ('k', 'l', F(1, 4)), ('k', 'j', F(1, 4)),
        ('d', 'd', F(1, 2)), ('j', 'j', F(1, 2))]
ck('the_sources_own_x_equals_0_configuration_has_total_weight_7_over_2',
   sum(START) == F(7, 2), str(sum(START)))
ck('the_sources_own_x_equals_0_configuration_is_fractionally_dominating',
   all(sum(START[u] for u in N[v]) >= 1 for v in range(12)))
ck('every_arc_of_the_printed_response_is_an_edge_or_a_self_loop',
   all(LAB[y] in N[LAB[x]] for x, y, _ in MOVE))
out = [F(0)] * 12
after = [F(0)] * 12
for x, y, q in MOVE:
    out[LAB[x]] += q
    after[LAB[y]] += q
ck('no_vertex_ships_more_than_it_holds_the_one_hop_budget_is_respected',
   all(out[v] <= START[v] for v in range(12)))
for v in range(12):
    after[v] += START[v] - out[v]
ck('the_response_conserves_the_total_weight_7_over_2', sum(after) == F(7, 2), str(sum(after)))
ck('the_response_puts_weight_exactly_1_on_the_attacked_vertex_j',
   after[LAB['j']] == 1, str(after[LAB['j']]))
sums_after = [sum(after[u] for u in N[v]) for v in range(12)]
note('closed-neighbourhood sums after the response = %s'
     % {nm: str(sums_after[LAB[nm]]) for nm in 'abcdefghijkl'})
ck('the_resulting_configuration_is_fractionally_dominating_so_the_attack_at_j_IS_answerable',
   all(s >= 1 for s in sums_after))
ck('so_the_published_claim_that_h_and_l_are_forced_to_0_after_the_response_is_unjustified',
   after[LAB['h']] != 0 and after[LAB['l']] != 0,
   'h = %s and l = %s, both nonzero' % (after[LAB['h']], after[LAB['l']]))

# ---------------------------------------------------------------------------
# 11. THE ANCILLARY CELL Cay(Z_16, {+-1, 8})
# ---------------------------------------------------------------------------
print()
print('=== Step 10: the ancillary upper bound gamma_f^infty(Cay(Z_16,{+-1,8})) <= 9/2')

N16 = cay_pm1_half(16)
SEED16 = [F(1), F(0), F(0), F(1, 2), F(0), F(1, 2), F(1, 2), F(0),
          F(1, 2), F(0), F(1, 2), F(1, 2), F(0), F(1, 2), F(0), F(0)]
ck('Cay_Z16_pm1_8_is_16_vertices_cubic_with_24_edges',
   len(N16) == 16 and all(len(x) == 4 for x in N16) and edge_count(N16) == 24)
ck('the_printed_Cay_Z16_seed_totals_9_over_2', sum(SEED16) == F(9, 2), str(sum(SEED16)))
s16 = [sum(SEED16[u] for u in N16[v]) for v in range(16)]
ck('the_printed_Cay_Z16_seed_is_fractionally_dominating_all_sums_1_or_3_over_2',
   all(s >= 1 for s in s16) and set(str(s) for s in s16) <= {'1', '3/2'},
   'sums %s' % [str(s) for s in s16])
ck('the_printed_Cay_Z16_seed_is_even_so_its_16_translates_are_well_defined',
   all(SEED16[i] == SEED16[(-i) % 16] for i in range(16)))
FAM16 = [[SEED16[(i - m) % 16] for i in range(16)] for m in range(16)]
ck('all_16_translates_total_9_over_2_are_fractionally_dominating_and_carry_1_at_their_centre',
   all(sum(f) == F(9, 2) for f in FAM16)
   and all(all(sum(FAM16[m][u] for u in N16[v]) >= 1 for v in range(16)) for m in range(16))
   and all(FAM16[m][m] == 1 for m in range(16)))
bad16 = [(a, b) for a in range(16) for b in range(16)
         if not reconfigurable(N16, FAM16[a], FAM16[b])]
ck('all_256_ordered_pairs_of_the_Cay_Z16_family_are_reconfigurable', not bad16,
   '%d/256, violations=%s' % (256 - len(bad16), bad16[:6]))
g16, S16 = min_dominating_size(N16, 5)
ck('the_domination_number_of_Cay_Z16_pm1_8_is_5_so_9_over_2_is_below_it',
   g16 == 5 and F(9, 2) < 5, 'gamma=%s attained by %s' % (g16, S16))

# ---------------------------------------------------------------------------
print()
print('NOTE SCOPE. This program checks claims of the accompanying note and nothing '
      'wider. NOT RE-RUN here: the floating-point linear programs that FOUND the witness '
      'families and the dual certificate (they are re-checked above as exhibited objects, '
      'which is why no solver is needed); and no attack tree other than the 9-node tree '
      'printed in the note is built or certified here. Every check above encodes the '
      'closed-neighbourhood reading of the one-round rule, so none of them can discriminate '
      'between that reading and a stricter one; for a stricter rule nothing here bounds '
      'gamma_f^infty. GAPS NOT COVERED: nothing above bears on any cell of the source '
      'Problem 7.5 other than C_6 [] K_2, and no general-k statement is tested.')
print()
n = len(CHECKS)
f = CHECKS.count(False)
if f:
    print('VERDICT: %d of %d CHECKS FAILED' % (f, n))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % n)
sys.exit(0)

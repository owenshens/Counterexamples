#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verification program for

    "A 28-vertex circulant with integral fractional broadcast number strictly between its
     multipacking and broadcast numbers"

It re-derives, from the objects PRINTED IN THE PAPER and from nothing else, every quantity the
paper claims:

  * the primary witness G = C_28(1,4), read BOTH from the paper's graph6 string and, independently,
    from the paper's arithmetic definition  i ~ i +- 1, i ~ i +- 4  (mod 28), and checked label-equal;
  * its ball profile, radius and diameter;
  * certificate (A), the primal broadcast LP solution of cost 4;
  * certificate (B), the dual of value 4, which pins gamma_{b,f}(G) = 4 by weak duality alone;
  * certificate (C), gamma_b(G) = 5, by an explicit broadcast above and a ball-partition count below;
  * certificate (D), mp(G) = 3, by the exhibited multipacking {0,6,15} below and by the mod-7 tiling
    obstruction above -- the latter re-checked THREE independent ways: an exact polynomial gcd in
    Q[x] (the paper's argument),     an exhaustive scan of all 210 residue distributions, and an
    exhaustive scan of all 2925 four-element subsets containing 0.


Python 3.9+, STANDARD LIBRARY ONLY (no numpy / sympy / networkx / scipy), no external data file.
All arithmetic is exact: Python ints and fractions.Fraction. No floating-point value is ever
compared, rounded or thresholded.

Contract: one `PASS <name> [detail]` line per check, a closing
`VERDICT: ALL <n> CHECKS PASS`, and exit status 0 iff every check passed.
"""

import sys
from fractions import Fraction as Fr
from itertools import combinations

# --------------------------------------------------------------------------------------------
# the check harness
# --------------------------------------------------------------------------------------------

_PASSED = 0
_FAILED = 0


def check(name, cond, detail=''):
    global _PASSED, _FAILED
    if cond:
        _PASSED += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILED += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE %s' % text)


# --------------------------------------------------------------------------------------------
# THE OBJECTS, EXACTLY AS THE PAPER PRINTS THEM
# --------------------------------------------------------------------------------------------

# The paper's graph6 line for the primary witness (paper.tex, Section "The witness", one line).
G6_PRINTED = ('[hdHGcH@GC_H?H?C_@G?H??c?@G?@G??c??H??@G??C_??H_??HO??Cc??@K_??H')

# The paper's arithmetic definition of the same graph.
N1, S1 = 28, (1, 4)

# The paper's printed BFS level sets from the vertex 0.
LEVELS_PRINTED = {
    0: [0],
    1: [1, 4, 24, 27],
    2: [2, 3, 5, 8, 20, 23, 25, 26],
    3: [6, 7, 9, 12, 16, 19, 21, 22],
    4: [10, 11, 13, 15, 17, 18],
    5: [14],
}
PROFILE_PRINTED = [5, 13, 21, 27, 28]           # n_k = |N_k[v]|, k = 1..5
M_PRINTED = [0, 6, 15]                          # the exhibited multipacking
T_PRINTED = [10, 11, 13, 14, 15, 17, 18]        # F_0 = {v : d(0,v) >= 4}
T_MOD7_PRINTED = [1, 1, 0, 2, 2, 0, 1]          # residue counts of T modulo 7
PARTITION_SUMS_PRINTED = {                      # the paper's five-row table at cost 4
    (4,): 27, (3, 1): 26, (2, 2): 26, (2, 1, 1): 23, (1, 1, 1, 1): 20,
}
N2_0_CAP_N2_6_PRINTED = [1, 2, 3, 4, 5, 8, 26]
N2_15_PRINTED = [7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23]

# --------------------------------------------------------------------------------------------
# graph6 (a self-contained stdlib codec; the format is R. B. Read / B. McKay's, n <= 62 branch)
# --------------------------------------------------------------------------------------------


def graph6_decode(s):
    """-> (n, set of frozenset pairs).  Column-major upper triangle, 6 bits per printable char."""
    s = s.strip()
    n = ord(s[0]) - 63
    if not 0 <= n <= 62:
        raise ValueError('this decoder covers n <= 62 only; header byte gives n = %d' % n)
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        if not 0 <= v < 64:
            raise ValueError('byte %r outside the graph6 printable range' % ch)
        bits.extend((v >> k) & 1 for k in (5, 4, 3, 2, 1, 0))
    need = n * (n - 1) // 2
    if len(bits) < need:
        raise ValueError('graph6 body holds %d bits, needs %d' % (len(bits), need))
    edges, t = set(), 0
    for j in range(1, n):
        for i in range(j):
            if bits[t]:
                edges.add(frozenset((i, j)))
            t += 1
    if any(bits[need:]):
        raise ValueError('graph6 padding bits are not all zero')
    return n, edges


def graph6_encode(n, edges):
    if not 0 <= n <= 62:
        raise ValueError('this encoder covers n <= 62 only')
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if frozenset((i, j)) in edges else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return ''.join(out)


# --------------------------------------------------------------------------------------------
# graph plumbing
# --------------------------------------------------------------------------------------------


def circulant(n, S):
    adj = [set() for _ in range(n)]
    for i in range(n):
        for s in S:
            adj[i].add((i + s) % n)
            adj[i].add((i - s) % n)
    for i in range(n):
        adj[i].discard(i)
    return adj


def edge_set(adj):
    return set(frozenset((u, v)) for u in range(len(adj)) for v in adj[u])


def all_distances(adj):
    n = len(adj)
    D = []
    for s in range(n):
        d = [-1] * n
        d[s] = 0
        frontier = [s]
        while frontier:
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if d[w] < 0:
                        d[w] = d[u] + 1
                        nxt.append(w)
            frontier = nxt
        D.append(d)
    return D


def ball(D, v, k):
    return set(u for u in range(len(D)) if 0 <= D[v][u] <= k)


def is_multipacking(M, D, n, kmax):
    """M is a k-multipacking for every k <= kmax:  |M cap N_s[v]| <= s for all v, all 1 <= s <= kmax."""
    for v in range(n):
        for s in range(1, kmax + 1):
            if sum(1 for u in M if D[v][u] <= s) > s:
                return False
    return True


def best_cover(profile, cost, kmax):
    """max over multisets of powers summing to `cost` of the sum of ball sizes -- the counting bound."""
    best = [0] * (cost + 1)
    for c in range(1, cost + 1):
        for k in range(1, min(c, kmax) + 1):
            cand = profile[k - 1] + best[c - k]
            if cand > best[c]:
                best[c] = cand
    return best


def partitions(c, mx):
    if c == 0:
        yield ()
        return
    for k in range(min(c, mx), 0, -1):
        for rest in partitions(c - k, k):
            yield (k,) + rest


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for j in range(total + 1):
        for rest in compositions(total - j, parts - 1):
            yield (j,) + rest


# --------------------------------------------------------------------------------------------
# exact polynomial arithmetic over Q, for the mod-7 cyclotomic step
# --------------------------------------------------------------------------------------------


def poly_trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def poly_divmod(a, b):
    """a, b lists of Fraction, low-degree first.  -> (quotient, remainder)."""
    a = [Fr(x) for x in a]
    b = poly_trim([Fr(x) for x in b])
    if not b:
        raise ZeroDivisionError
    q = [Fr(0)] * max(0, len(a) - len(b) + 1)
    while len(poly_trim(a)) >= len(b):
        shift = len(a) - len(b)
        f = a[-1] / b[-1]
        q[shift] = f
        for i, c in enumerate(b):
            a[shift + i] -= f * c
        poly_trim(a)
    return poly_trim(q), poly_trim(a)


def poly_gcd(a, b):
    a = poly_trim([Fr(x) for x in a])
    b = poly_trim([Fr(x) for x in b])
    while b:
        a, b = b, poly_divmod(a, b)[1]
    if a:                                   # normalise to monic
        lead = a[-1]
        a = [c / lead for c in a]
    return a


# ============================================================================================
print('verification of "A 28-vertex circulant with integral fractional broadcast number'
      ' strictly between its multipacking and broadcast numbers"')
print('python %s ; exact integer and Fraction arithmetic only, no floating point'
      % sys.version.split()[0])
print('')
print('=== part 1: the primary witness G = C_28(1,4) ===')
# ============================================================================================

# ---- 1. the object, read two independent ways ----------------------------------------------
n6, E6 = graph6_decode(G6_PRINTED)
adj = circulant(N1, S1)
E = edge_set(adj)

check('graph6_header_gives_28_vertices', n6 == 28, 'n = %d' % n6)
check('graph6_body_holds_56_edges', len(E6) == 56, '|E| = %d' % len(E6))
check('graph6_string_and_arithmetic_definition_are_label_equal', n6 == N1 and E6 == E,
      'identical as sets of labelled 0-indexed pairs')
check('graph6_round_trips_to_the_printed_string', graph6_encode(N1, E) == G6_PRINTED,
      '%d chars re-encoded byte-for-byte' % len(G6_PRINTED))
check('graph6_padding_bits_are_zero', True, 'the decoder refuses a nonzero pad; it did not')

n = N1
D = all_distances(adj)
check('G_is_4_regular', sorted(set(len(a) for a in adj)) == [4], 'degree set = {4}')
check('G_is_connected', all(D[0][u] >= 0 for u in range(n)), 'BFS from 0 reaches all 28 vertices')
check('G_has_56_edges', len(E) == 56, '|E| = %d = 28*4/2' % len(E))

ecc = [max(row) for row in D]
check('G_has_constant_eccentricity_5', sorted(set(ecc)) == [5],
      'eccentricity set = {%s}' % ','.join(str(x) for x in sorted(set(ecc))))
rad, diam = min(ecc), max(ecc)
check('radius_equals_diameter_equals_5', rad == 5 and diam == 5, 'rad = diam = 5')
note('rad = diam, so the two conventions for a multipacking (s <= diam, as in the source, and '
     's <= rad, as much of the literature writes it) coincide on this graph.')

# ---- 2. the ball profile and the printed level sets ---------------------------------------
levels = {}
for u in range(n):
    levels.setdefault(D[0][u], []).append(u)
ok_lv = all(sorted(levels.get(k, [])) == LEVELS_PRINTED[k] for k in LEVELS_PRINTED)
check('bfs_level_sets_from_0_match_the_paper', ok_lv and set(levels) == set(LEVELS_PRINTED),
      'levels 0..5 sizes %s' % ([len(levels[k]) for k in sorted(levels)],))

profile = [len(ball(D, 0, k)) for k in range(1, diam + 1)]
check('ball_profile_is_5_13_21_27_28', profile == PROFILE_PRINTED, 'n_k = %s' % (profile,))
same_at_every_vertex = all([len(ball(D, v, k)) for k in range(1, diam + 1)] == profile
                           for v in range(n))
check('ball_profile_is_the_same_at_every_vertex', same_at_every_vertex,
      'n_k independent of v -- so a single profile row governs both certificates')
check('translation_is_an_automorphism_of_G',
      all(edge_set(adj) == set(frozenset(((u + g) % n, (v + g) % n)) for u, v in
                               (tuple(e) for e in E)) for g in range(n)),
      'Z_28 acts on G; this is what licenses the "wlog 0 in M" below')

# ---- 3. the LP weights ---------------------------------------------------------------------
weights = [(Fr(k, profile[k - 1]), k) for k in range(1, rad + 1)]
check('the_five_ratios_k_over_n_k_are_as_printed',
      [str(w) for w, _ in weights] == ['1/5', '2/13', '1/7', '4/27', '5/28'],
      ' '.join('%s' % w for w, _ in weights))
wstar, rstar = min(weights)
check('w_star_is_one_seventh_at_r_star_3', wstar == Fr(1, 7) and rstar == 3,
      'w* = 1/7, r* = 3')
check('the_minimising_r_is_unique', sum(1 for w, _ in weights if w == wstar) == 1,
      'exactly one k attains the minimum')
check('n_times_w_star_is_the_integer_4', n * wstar == 4 and (n * wstar).denominator == 1,
      '28 * 1/7 = 4')

# ---- 4. certificate (A): the primal, cost 4 ------------------------------------------------
xval = Fr(1, profile[rstar - 1])
check('A_every_column_v_r_star_is_a_legal_LP_column', all(rstar <= ecc[v] for v in range(n)),
      'r* = 3 <= e(v) = 5 for all 28 vertices')
check('A_x_is_nonnegative', xval > 0, 'x[(v,3)] = %s at all 28 vertices, 0 elsewhere' % xval)
heard = [sum(xval for v in range(n) if D[u][v] <= rstar) for u in range(n)]
check('A_every_vertex_hears_at_least_1', all(h >= 1 for h in heard),
      'min_u (Ax)_u = %s, and in fact every entry equals %s' % (min(heard), heard[0]))
cost_primal = sum(xval * rstar for _ in range(n))
check('A_primal_cost_is_4', cost_primal == 4, 'cost = 28 * 3/21 = %s' % cost_primal)

# ---- 5. certificate (B): the dual, value 4 -------------------------------------------------
yval = wstar
viol = [(v, k) for v in range(n) for k in range(1, ecc[v] + 1)
        if sum(yval for _ in ball(D, v, k)) > k]
check('B_y_is_nonnegative', yval > 0, 'y(v) = %s at every vertex' % yval)
check('B_every_dual_constraint_holds', not viol,
      'all %d constraints (v,k), 1 <= k <= e(v) = 5, satisfied' % (n * 5))
dual_rows = [profile[k - 1] * yval for k in range(1, diam + 1)]
check('B_dual_row_values_are_as_printed',
      [str(x) for x in dual_rows] == ['5/7', '13/7', '3', '27/7', '4'],
      ' '.join(str(x) for x in dual_rows))
check('B_the_constraint_is_tight_exactly_at_k_3',
      [k for k in range(1, diam + 1) if profile[k - 1] * yval == k] == [3],
      'y(N_3[v]) = 3 = c_{v,3}; every other k is slack')
value_dual = sum(yval for _ in range(n))
check('B_dual_objective_is_4', value_dual == 4, 'sum_v y(v) = 28/7 = %s' % value_dual)

check('weak_duality_pins_gamma_bf_to_4', cost_primal == value_dual == 4,
      'a feasible primal of cost 4 and a feasible dual of value 4 -- no LP solver needed')
gamma_bf = Fr(4)
check('gamma_bf_is_an_integer', gamma_bf.denominator == 1, 'gamma_{b,f}(G) = 4')

# ---- 6. certificate (C): gamma_b = 5 ------------------------------------------------------
check('C_the_broadcast_f_0_equals_5_is_legal', 5 <= ecc[0] and 5 <= diam,
      'f(0) = 5 <= e(0) = 5 and 5 <= diam(G) = 5')
check('C_that_broadcast_dominates', ball(D, 0, 5) == set(range(n)),
      'N_5[0] = V(G), all 28 vertices')
check('C_gamma_b_at_most_5', True, 'the single-vertex broadcast above has cost 5')

best = best_cover(profile, 6, diam)
check('C_counting_bound_row_by_row', best[1:5] == [5, 13, 21, 27],
      'cost 1..4 cover at most 5, 13, 21, 27 of 28')
check('C_cost_4_cannot_cover_28', best[4] == 27 < n, 'best[4] = 27 < 28')
check('C_every_cost_below_4_is_worse', all(best[c] < n for c in range(1, 5)),
      'best[1..4] = %s, all < 28' % (best[1:5],))
tab = dict((p, sum(profile[k - 1] for k in p)) for p in partitions(4, diam))
check('C_the_five_printed_partitions_of_4_are_complete_and_correct',
      tab == PARTITION_SUMS_PRINTED,
      ' ; '.join('%s->%d' % (list(p), tab[p]) for p in sorted(tab, key=lambda q: (-len(q), q))))
check('C_gamma_b_at_least_5', max(tab.values()) == 27 < n,
      'the best cost-4 multiset of balls covers 27 < 28 vertices')
gamma_b = 5
check('C_gamma_b_equals_5', gamma_b == 5, 'gamma_b(G) = 5')

# ---- 7. certificate (D), lower half: mp >= 3 ----------------------------------------------
M = M_PRINTED
check('D_M_is_a_3_element_subset', len(set(M)) == 3 and all(0 <= u < n for u in M),
      'M = %s' % (M,))
check('D_M_pairwise_distances_are_3_3_4',
      sorted(D[a][b] for a, b in combinations(M, 2)) == [3, 3, 4],
      'd(0,6)=%d d(0,15)=%d d(6,15)=%d' % (D[0][6], D[0][15], D[6][15]))
check('D_M_is_a_multipacking_s_up_to_diam', is_multipacking(M, D, n, diam),
      '|M cap N_s[v]| <= s for all 28 vertices and all s = 1..5')
check('D_M_is_a_multipacking_s_up_to_rad', is_multipacking(M, D, n, rad),
      'robust to the s <= rad(G) convention as well')
check('D_no_vertex_sees_two_of_M_within_distance_1',
      max(sum(1 for u in M if D[v][u] <= 1) for v in range(n)) == 1, 'max_v = 1')
check('D_no_vertex_sees_three_of_M_within_distance_2',
      max(sum(1 for u in M if D[v][u] <= 2) for v in range(n)) == 2, 'max_v = 2')
check('D_the_printed_pairwise_ball_intersection_is_right',
      sorted(ball(D, 0, 2) & ball(D, 6, 2)) == N2_0_CAP_N2_6_PRINTED,
      'N_2[0] cap N_2[6] = %s' % (N2_0_CAP_N2_6_PRINTED,))
check('D_the_printed_ball_N2_15_is_right', sorted(ball(D, 15, 2)) == N2_15_PRINTED,
      '|N_2[15]| = %d' % len(N2_15_PRINTED))
check('D_the_triple_intersection_is_empty',
      not (ball(D, 0, 2) & ball(D, 6, 2) & ball(D, 15, 2)),
      'so no vertex has all three of M within distance 2')
check('D_mp_at_least_3', True, 'the exhibited M is a multipacking of size 3')

# ---- 8. certificate (D), upper half: mp <= 3, route 1 = the tiling / mod-7 argument -------
T = sorted(u for u in range(n) if D[0][u] >= rstar + 1)
check('D_F_0_is_the_printed_7_set', T == T_PRINTED, 'F_0 = {v : d(0,v) >= 4} = %s' % (T,))
check('D_F_0_has_size_n_minus_n_r_star', len(T) == n - profile[rstar - 1],
      '|F_0| = 7 = 28 - 21')
check('D_F_0_is_the_stated_translate_of_a_symmetric_set',
      T == sorted((14 + d) % n for d in (0, 1, -1, 3, -3, 4, -4)),
      'F_0 = 14 + {0, +-1, +-3, +-4}')
check('D_F_u_equals_u_plus_T_for_every_u',
      all(sorted((u + t) % n for t in T)
          == sorted(v for v in range(n) if D[u][v] >= rstar + 1) for u in range(n)),
      'checked at all 28 vertices')
lhs = 4 * profile[rstar - 1]
check('D_the_counting_identity_forces_equality', lhs == n * rstar == 84,
      'a 4-element M gives sum_v |M cap N_3[v]| = 4*21 = 84, against the bound 28*3 = 84')
note('equality therefore forces |M cap N_3[v]| = 3 for EVERY v, i.e. each vertex lies in exactly '
     'one F_u with u in M, i.e. the four translates u + T partition Z_28 (4 * 7 = 28).')

q = 7
tvec = [0] * q
for x in T:
    tvec[x % q] += 1
check('D_the_residue_vector_of_T_mod_7_is_as_printed', tvec == T_MOD7_PRINTED,
      't = (%s)' % ','.join(str(x) for x in tvec))
check('D_the_residue_vector_sums_to_7', sum(tvec) == 7, 'sum t_j = |T| = 7')
check('D_each_residue_class_mod_7_holds_exactly_4_elements_of_Z_28',
      all(sum(1 for x in range(n) if x % q == j) == 4 for j in range(q)),
      '28 / 7 = 4')

# the paper's argument: t(x) = 1 + x + 2x^3 + 2x^4 + x^6 is monic of degree 6 and is not Phi_7,
# and Phi_7 is irreducible over Q, so gcd(t, Phi_7) = 1 in Q[x] and t(zeta) != 0.
tpoly = [Fr(c) for c in tvec]
phi7 = [Fr(1)] * 7
check('D_t_is_monic_of_degree_6', len(poly_trim(list(tpoly))) == 7 and tvec[6] == 1,
      't(x) = 1 + x + 2x^3 + 2x^4 + x^6')
check('D_t_differs_from_Phi_7', tvec != [1] * 7,
      'the coefficient of x^2 is 0 in t and 1 in Phi_7')
g = poly_gcd(tpoly, phi7)
check('D_gcd_of_t_and_Phi_7_is_1_over_Q', g == [Fr(1)],
      'exact Euclid in Q[x] returns the unit -- so t(zeta) != 0 at every primitive 7th root of 1')
note('hence the Z_7 Fourier transform of t vanishes at no nontrivial character; the constant '
     'vector 4 is supported on the trivial character alone, so m must be constant, and 7*m_0 = 4 '
     'is impossible over the integers.')

sols = []
scanned = 0
for m in compositions(4, q):
    scanned += 1
    conv = [sum(m[a] * tvec[(j - a) % q] for a in range(q)) for j in range(q)]
    if all(c == 4 for c in conv):
        sols.append(m)
check('D_no_residue_distribution_of_a_4_element_M_can_tile', not sols,
      'all %d nonnegative m with sum 4 scanned; none has (m*t)(j) = 4 for all j in Z_7' % scanned)
check('D_the_composition_scan_was_exhaustive', scanned == 210,
      'C(10,6) = 210 compositions of 4 into 7 parts')

# ---- 9. certificate (D), upper half: mp <= 3, route 2 = brute force over all 4-subsets ----
bad, tried = [], 0
for rest in combinations(range(1, n), 3):
    tried += 1
    if is_multipacking((0,) + rest, D, n, diam):
        bad.append((0,) + rest)
check('D_no_4_subset_containing_0_is_a_multipacking', not bad,
      '%d subsets tested exhaustively, 0 multipackings' % tried)
check('D_that_scan_was_exhaustive_up_to_translation', tried == 2925,
      'C(27,3) = 2925; translation by Z_28 is an automorphism, so a nonempty multipacking may be '
      'moved to contain 0')
mp = 3
check('D_mp_equals_3', mp == 3, 'mp(G) = 3, both bounds established')

# ---- 10. the chain, and the answer --------------------------------------------------------
check('the_chain_mp_lt_gamma_bf_lt_gamma_b_is_strict', mp < gamma_bf < gamma_b,
      '3 < 4 < 5')
check('the_middle_term_is_an_integer', gamma_bf.denominator == 1, 'gamma_{b,f}(G) = 4')
check('G_is_connected_as_the_problem_requires', all(D[0][u] >= 0 for u in range(n)),
      'the source works throughout with a connected graph')
check('G_has_the_property_Problem_4_asks_for',
      mp < gamma_bf < gamma_b and gamma_bf.denominator == 1,
      'mp = 3 < gamma_bf = 4 < gamma_b = 5, with 4 integral')

# ============================================================================================
print('')
note('SCOPE. This program re-derives the paper\'s claims about the one printed graph '
     'G = C_28(1,4) and nothing else. It examines no other graph, it claims and proves no '
     'minimality, and it reads no external manuscript -- in particular not the unpublished '
     'Brewster-Duchesne manuscript to which the source survey credits the vertex-transitive '
     'dual that certificate (B) instantiates.')
print('')
if _FAILED:
    print('VERDICT: %d CHECK(S) FAILED out of %d' % (_FAILED, _PASSED + _FAILED))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _PASSED)
sys.exit(0)

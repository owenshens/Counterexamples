#!/usr/bin/env python3
"""Verification program for the census arithmetic behind

    A 3-vertex census refuting the k-monotonicity clause of a conjecture of
    D'Angeli and Rodaro on typical total simplicity

What this program establishes is the exact arithmetic: the exhaustive 3-vertex census in
the ordered/with-replacement model, giving P(3,2) = 57/148 > P(3,3) = 2566/6877, i.e. one
non-monotone k-step in that reading of the model.  The wording of the k-monotonicity
clause and its location (arXiv:2607.17335v1, totally_main.tex line 1463) are TRANSCRIBED
BY HAND here and are NOT parsed or verified: no source file is read.  Whether the ordered
model is the source's intended reading is likewise not decided here.

Python 3.9+, STANDARD LIBRARY ONLY (itertools, fractions, math, sys).  No third-party
package, no external data file, no network.  Every decision is made in exact integer or
exact rational arithmetic; no floating-point number decides anything (floats appear only
inside f-strings for the reader's convenience).

WHAT IT READS.  The objects the paper exhibits are (i) the model -- ordered k-out graphs
on the vertex set {0,...,n-1}, each vertex choosing k heads independently and uniformly
WITH REPLACEMENT -- and (ii) the finitely many integers the paper prints.  The integers
the paper claims are transcribed once, in the block CLAIMED below, and every one of them
is then re-derived here from scratch by exhaustive enumeration and, independently, by the
closed forms the paper prints.

Prints one `PASS <name> [detail]` line per check and closes with
    VERDICT: ALL <n> CHECKS PASS
exiting 0 if and only if every check passed.
"""

import itertools
import sys
from fractions import Fraction
from math import factorial

# ----------------------------------------------------------------------------------
# THE PAPER'S CLAIMED VALUES, TRANSCRIBED (nothing below is derived from this block;
# it is only compared against).
# ----------------------------------------------------------------------------------
CLAIMED = {
    # (n, k): (A, B)   A = weight of {strongly connected AND totally simple},
    #                  B = weight of {strongly connected}
    (3, 2): (114, 296),
    (3, 3): (5132, 13754),
    (3, 4): (180890, 458000),
    (4, 2): (6720, 20958),
}
CLAIMED_P = {
    (3, 2): Fraction(57, 148),
    (3, 3): Fraction(2566, 6877),
    (3, 4): Fraction(18089, 45800),
    (4, 2): Fraction(160, 499),
}
CLAIMED_DROP = Fraction(12221, 1017796)
CLAIMED_TRACE = {2: 150, 3: 5798, 4: 189518}
CLAIMED_CORRECTION = {2: 36, 3: 666, 4: 8628}
CLAIMED_RIVAL_P = {2: Fraction(34, 65), 3: Fraction(133, 243), 4: Fraction(87, 148)}
CLAIMED_A027834 = {3: 296, 4: 20958}          # OEIS A027834(n) = B(n,2)
CLAIMED_BELL = [2, 5, 15, 52, 203, 877]       # Bell(2..7)
KMAX_CLOSED = 30                              # the closed forms are exercised on k = 2..30

_n_pass = 0
_n_fail = 0


def ck(name, cond, detail=''):
    global _n_pass, _n_fail
    if cond:
        _n_pass += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _n_fail += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(msg):
    print('NOTE %s' % msg)


# ----------------------------------------------------------------------------------
# combinatorial primitives
# ----------------------------------------------------------------------------------
def multinomial(k, counts):
    """k! / prod(c!) -- the weight of one out-multiset in the ORDERED model."""
    r = factorial(k)
    for c in counts:
        r //= factorial(c)
    return r


def count_vectors(n, k):
    """All (c_0,...,c_{n-1}) of nonnegative integers with sum k, i.e. all out-multisets."""
    out = []

    def rec(pos, rem, cur):
        if pos == n - 1:
            out.append(tuple(cur + [rem]))
            return
        for x in range(rem + 1):
            rec(pos + 1, rem - x, cur + [x])
    rec(0, k, [])
    return out


def ordered_rows(n, k):
    """All k-tuples of heads in V^k, collapsed to their count vector.  Weight 1 each --
    this is the DIRECT enumeration over raw edge tuples, carrying no multinomial weight
    anywhere, and it is an independent route to the same numbers."""
    return [tuple(t.count(w) for w in range(n)) for t in itertools.product(range(n), repeat=k)]


def set_partitions(n):
    """Every partition of {0,...,n-1} as a tuple of frozensets."""
    res = []

    def rec(i, blocks):
        if i == n:
            res.append(tuple(frozenset(b) for b in blocks))
            return
        for j in range(len(blocks)):
            blocks[j].append(i)
            rec(i + 1, blocks)
            blocks[j].pop()
        blocks.append([i])
        rec(i + 1, blocks)
        blocks.pop()
    rec(0, [])
    return res


def nontrivial_partitions(n):
    """All partitions OTHER than the discrete one (Delta, n blocks) and the universal one
    (Nabla, 1 block).  Both are excluded by the source: a congruence different from
    Delta_G and Nabla_G is what the source calls nontrivial, and an automaton is `simple`
    when its only congruences are the discrete and the universal one."""
    return [p for p in set_partitions(n) if len(p) not in (1, n)]


def is_lumpable(mat, part):
    """The source's condition: for every pair of blocks C_i, C_j and all u,u' in C_i,
    |delta_u(C_j)| = |delta_u'(C_j)|."""
    for blk in part:
        if len(blk) < 2:
            continue                              # vacuous
        us = sorted(blk)
        for cj in part:
            base = None
            for u in us:
                s = sum(mat[u][w] for w in cj)
                if base is None:
                    base = s
                elif s != base:
                    return False
    return True


def is_totally_simple(mat, nts):
    """No partition other than Delta and Nabla is lumpable."""
    for part in nts:
        if is_lumpable(mat, part):
            return False
    return True


def is_strongly_connected(mat, n):
    """Reachability from every vertex to every vertex, on the SUPPORT of the out-multisets."""
    for v in range(n):
        seen = {v}
        stack = [v]
        while stack:
            u = stack.pop()
            for w in range(n):
                if mat[u][w] > 0 and w not in seen:
                    seen.add(w)
                    stack.append(w)
        if len(seen) != n:
            return False
    return True


def census(n, k, mode='weighted'):
    """Exhaustive census.

    mode='weighted'  -- states are out-multiset vectors, weight = product of multinomials.
                        This IS the ordered model, summed in closed form over each vertex.
    mode='ordered'   -- states are raw k-tuples of heads per vertex, weight 1 each.  A
                        second, independent route with no multinomial weight anywhere.
    mode='rival'     -- states are out-multiset vectors, weight 1 each.  NOT the paper's
                        model; the uniform-over-labelled-multigraphs reading, computed only
                        so the model-sensitivity disclosure can itself be checked.
    Returns dict with A, B, total_weight, n_states, ts_not_scc.
    """
    if mode == 'ordered':
        rows = ordered_rows(n, k)
        wts = [1] * len(rows)
    else:
        rows = count_vectors(n, k)
        wts = [multinomial(k, r) for r in rows] if mode == 'weighted' else [1] * len(rows)
    nts = nontrivial_partitions(n)
    A = B = tot = 0
    ts_not_scc = 0
    for combo in itertools.product(range(len(rows)), repeat=n):
        mat = [rows[i] for i in combo]
        w = 1
        for i in combo:
            w *= wts[i]
        tot += w
        sc = is_strongly_connected(mat, n)
        ts = is_totally_simple(mat, nts)
        if sc:
            B += w
            if ts:
                A += w
        elif ts:
            ts_not_scc += w
    return {'A': A, 'B': B, 'total_weight': tot, 'n_states': len(rows) ** n,
            'ts_not_scc': ts_not_scc}


def binom(a, b):
    if b < 0 or b > a:
        return 0
    return factorial(a) // (factorial(b) * factorial(a - b))


# ----------------------------------------------------------------------------------
# the paper's closed forms at n = 3
# ----------------------------------------------------------------------------------
def B_closed(k):
    """B(3,k) = (3^k - 1)^3 - 3 (2^k - 1)^2 (3^k - 1)."""
    return (3 ** k - 1) ** 3 - 3 * (2 ** k - 1) ** 2 * (3 ** k - 1)


def W_matrix(k):
    """W_k[s][t] = k!/(s! t! (k-s-t)!), zero when s+t > k."""
    return [[(factorial(k) // (factorial(s) * factorial(t) * factorial(k - s - t))
              if s + t <= k else 0) for t in range(k + 1)] for s in range(k + 1)]


def matmul(X, Y):
    p = len(Y)
    q = len(Y[0])
    return [[sum(X[i][r] * Y[r][j] for r in range(p)) for j in range(q)] for i in range(len(X))]


def trace_route(k):
    """(trace, correction, A) with A(3,k) = tr((W_k D)^3) - correction, D = J - I."""
    W = W_matrix(k)
    D = [[0 if i == j else 1 for j in range(k + 1)] for i in range(k + 1)]
    M = matmul(W, D)
    M3 = matmul(matmul(M, M), M)
    tr = sum(M3[i][i] for i in range(k + 1))
    f = [sum(W[t][q] for q in range(1, k + 1)) for t in range(k + 1)]
    corr = 3 * (sum(f) ** 2 - sum(x * x for x in f))
    return tr, corr, tr - corr


def incexc_route(k):
    """A(3,k) = (3^k-1)^3 - 3(3^k-1) S2 + 3 S3 - S4, the inclusion-exclusion closed form."""
    def Wf(x, y):
        return factorial(k) // (factorial(x) * factorial(y) * factorial(k - x - y)) if x + y <= k else 0

    def g(m):
        return (2 ** k - 1) if m == 0 else binom(k, m) * 2 ** (k - m)

    S2 = sum(g(m) ** 2 for m in range(k + 1))
    S3 = sum(Wf(p, q) * g(p) * g(q)
             for p in range(k + 1) for q in range(k + 1 - p) if (p, q) != (0, 0))
    S4 = 0
    for p in range(k + 1):
        for q in range(k + 1 - p):
            for r in range(k + 1):
                if q + r > k or r + p > k:
                    continue
                if [p, q, r].count(0) > 1:
                    continue
                S4 += Wf(p, q) * Wf(r, q) * Wf(r, p)
    A = (3 ** k - 1) ** 3 - 3 * (3 ** k - 1) * S2 + 3 * S3 - S4
    return S2, S3, S4, A


# ==================================================================================
def main():
    print('verification of: the exact 3-vertex census arithmetic behind a claimed')
    print("counterexample to the k-monotonicity clause of a conjecture attributed to")
    print('D\'Angeli-Rodaro (the clause text and its location at arXiv:2607.17335v1,')
    print('totally_main.tex line 1463 are TRANSCRIBED BY HAND, NOT PARSED HERE).')
    print('What is verified: P(3,2) = 57/148 > P(3,3) = 2566/6877 in the')
    print('ordered/with-replacement model. Whether that model is the source\'s intended')
    print('one, and whether line 1463 states the clause as transcribed, are NOT checked.')
    print('python %s, exact integer and rational arithmetic only' % sys.version.split()[0])
    print()

    # ---------------------------------------------------------------- Step 1
    print('=== Step 1: the predicate -- the partition lattice, and Delta/Nabla excluded')
    bell = [len(set_partitions(m)) for m in range(2, 8)]
    ck('bell_numbers_n_2_to_7_reproduced', bell == CLAIMED_BELL, 'Bell(2..7) = %s' % bell)
    ck('nontrivial_partition_count_n3_is_bell3_minus_2',
       len(nontrivial_partitions(3)) == 3 == bell[1] - 2, '3 = Bell(3) - 2')
    ck('nontrivial_partition_count_n4_is_bell4_minus_2',
       len(nontrivial_partitions(4)) == 13 == bell[2] - 2, '13 = Bell(4) - 2')
    # Delta and Nabla are ALWAYS lumpable, which is why the source must exclude both.
    rows32 = count_vectors(3, 2)
    delta = tuple(frozenset([i]) for i in range(3))
    nabla = (frozenset(range(3)),)
    both = all(is_lumpable([rows32[i] for i in c], delta) and
               is_lumpable([rows32[i] for i in c], nabla)
               for c in itertools.product(range(len(rows32)), repeat=3))
    ck('delta_and_nabla_are_lumpable_on_every_state', both,
       'all %d states at n=3,k=2; so excluding both is forced, not chosen' % len(rows32) ** 3)

    # the paper's hand criterion at n = 3, checked against the generic decider
    nts3 = nontrivial_partitions(3)
    mism = 0
    tested = 0
    for k in (2, 3, 4):
        rows = count_vectors(3, k)
        for c in itertools.product(range(len(rows)), repeat=3):
            mat = [rows[i] for i in c]
            hand = (mat[1][0] != mat[2][0]) and (mat[0][1] != mat[2][1]) and (mat[0][2] != mat[1][2])
            tested += 1
            if hand != is_totally_simple(mat, nts3):
                mism += 1
    ck('column_criterion_equals_generic_decider_at_n3', mism == 0,
       '%d states at k=2,3,4; %d mismatches' % (tested, mism))
    print()

    # ---------------------------------------------------------------- Step 2
    print('=== Step 2: the exhaustive census in the paper\'s ordered model, n = 3')
    cen = {}
    for k in (2, 3, 4):
        cen[k] = census(3, k, 'weighted')
        c = cen[k]
        ck('total_weight_identity_n3_k%d' % k, c['total_weight'] == 3 ** (3 * k),
           'enumerated weight %d = 3^%d' % (c['total_weight'], 3 * k))
        ck('state_count_identity_n3_k%d' % k, c['n_states'] == binom(k + 2, k) ** 3,
           '%d = C(%d,%d)^3' % (c['n_states'], k + 2, k))
    for k in (2, 3, 4):
        A, Bv = CLAIMED[(3, k)]
        ck('B_3_%d_equals_%d' % (k, Bv), cen[k]['B'] == Bv, 'strongly connected weight')
        ck('A_3_%d_equals_%d' % (k, A), cen[k]['A'] == A,
           'strongly connected AND totally simple weight')
    print()

    # ---------------------------------------------------------------- Step 3
    print('=== Step 3: a SECOND enumeration, direct over raw ordered edge tuples (weight 1)')
    for k in (2, 3):
        o = census(3, k, 'ordered')
        ck('ordered_tuple_total_is_3_to_the_%d' % (3 * k), o['total_weight'] == 3 ** (3 * k),
           '%d raw tuples enumerated' % o['total_weight'])
        ck('ordered_tuple_route_agrees_at_n3_k%d' % k,
           (o['A'], o['B']) == (cen[k]['A'], cen[k]['B']),
           'A=%d B=%d, no multinomial weight anywhere' % (o['A'], o['B']))
    print()

    # ---------------------------------------------------------------- Step 4
    print('=== Step 4: the closed forms the paper prints, checked against the census')
    ck('B_closed_form_matches_census_n3_k2_3_4',
       all(B_closed(k) == cen[k]['B'] for k in (2, 3, 4)),
       'B(3,k) = (3^k-1)^3 - 3(2^k-1)^2(3^k-1) = %s' % [B_closed(k) for k in (2, 3, 4)])
    tr_ok = corr_ok = A_ok = True
    for k in (2, 3, 4):
        tr, corr, A = trace_route(k)
        tr_ok &= (tr == CLAIMED_TRACE[k])
        corr_ok &= (corr == CLAIMED_CORRECTION[k])
        A_ok &= (A == cen[k]['A'])
    ck('transfer_matrix_traces_are_150_5798_189518', tr_ok,
       'tr((W_k D)^3) for k=2,3,4')
    ck('scc_corrections_are_36_666_8628', corr_ok, '3[(sum f)^2 - sum f^2]')
    ck('trace_minus_correction_equals_census_A', A_ok,
       'A(3,k) = %s' % [trace_route(k)[2] for k in (2, 3, 4)])
    ie_ok = True
    for k in (2, 3, 4):
        S2, S3, S4, A = incexc_route(k)
        ie_ok &= (A == cen[k]['A'])
    ck('inclusion_exclusion_closed_form_equals_census_A', ie_ok,
       'k=2: S2,S3,S4 = %s ; k=3: %s' % (incexc_route(2)[:3], incexc_route(3)[:3]))
    ck('two_independent_closed_forms_agree_for_k_2_to_%d' % KMAX_CLOSED,
       all(trace_route(k)[2] == incexc_route(k)[3] for k in range(2, KMAX_CLOSED + 1)),
       'transfer-matrix route vs inclusion-exclusion route, %d values of k'
       % (KMAX_CLOSED - 1))
    ck('correction_term_equals_weight_of_totally_simple_but_not_scc',
       all(CLAIMED_CORRECTION[k] == cen[k]['ts_not_scc'] for k in (2, 3, 4)),
       '36, 666, 8628 -- two derivations aimed at different questions agree')
    print()

    # ---------------------------------------------------------------- Step 5
    print('=== Step 5: the exact comparison P(3,2) > P(3,3) in the ordered model')
    P = {k: Fraction(cen[k]['A'], cen[k]['B']) for k in (2, 3, 4)}
    for k in (2, 3, 4):
        ck('P_3_%d_equals_%s' % (k, str(CLAIMED_P[(3, k)]).replace('/', '_over_')),
           P[k] == CLAIMED_P[(3, k)], 'exact rational %s = %.9f' % (P[k], float(P[k])))
    ck('cross_multiplication_gives_a_STRICT_decrease_from_k2_to_k3',
       cen[2]['A'] * cen[3]['B'] > cen[3]['A'] * cen[2]['B'],
       '114*13754 = %d  >  5132*296 = %d' % (cen[2]['A'] * cen[3]['B'],
                                             cen[3]['A'] * cen[2]['B']))
    ck('the_decrease_is_decided_without_any_floating_point',
       P[2] > P[3] and isinstance(P[2] - P[3], Fraction),
       'Fraction comparison only')
    ck('drop_equals_12221_over_1017796', P[2] - P[3] == CLAIMED_DROP,
       '%s = %.8f, about 1.2 percentage points' % (CLAIMED_DROP, float(CLAIMED_DROP)))
    ck('the_dip_is_isolated_P_3_4_exceeds_P_3_3', P[4] > P[3],
       '%s > %s: monotonicity resumes at once' % (P[4], P[3]))
    # k = 1 is OUTSIDE the conjecture's own range and is deliberately not used
    c31 = census(3, 1, 'weighted')
    ck('B_3_1_equals_2_and_P_3_1_equals_1_computed_but_not_used',
       c31['B'] == 2 == factorial(2) and Fraction(c31['A'], c31['B']) == 1,
       'B(3,1) = (3-1)! = 2 cyclic permutations, P(3,1) = 1; k=1 is not used above')
    # the whole n = 3 line by closed form
    Pc = {k: Fraction(trace_route(k)[2], B_closed(k)) for k in range(2, KMAX_CLOSED + 1)}
    ck('closed_forms_reproduce_P_3_k_at_k_2_3_4', all(Pc[k] == P[k] for k in (2, 3, 4)))
    viol = [k for k in range(2, KMAX_CLOSED) if Pc[k] > Pc[k + 1]]
    ck('exactly_one_violation_on_the_whole_n3_line_k_2_to_%d' % KMAX_CLOSED,
       viol == [2], 'violating k-steps: %s' % viol)
    ck('P_3_k_strictly_increases_for_every_k_from_3_to_%d' % KMAX_CLOSED,
       all(Pc[k] < Pc[k + 1] for k in range(3, KMAX_CLOSED)),
       'P(3,3) = %.9f rising to P(3,%d) = %.9f' % (float(Pc[3]), KMAX_CLOSED,
                                                   float(Pc[KMAX_CLOSED])))
    print()

    # ---------------------------------------------------------------- Step 6
    print('=== Step 6: the n = 4, k = 2 census, and two transcribed literals compared')
    c42 = census(4, 2, 'weighted')
    ck('total_weight_identity_n4_k2', c42['total_weight'] == 4 ** 8,
       '%d = 4^8' % c42['total_weight'])
    ck('B_4_2_equals_20958_matching_hand_transcribed_A027834_term',
       c42['B'] == CLAIMED_A027834[4],
       'census value vs one transcribed literal, not a database lookup')
    ck('B_3_2_equals_296_matching_hand_transcribed_A027834_term',
       cen[2]['B'] == CLAIMED_A027834[3])
    ck('A_4_2_equals_6720', c42['A'] == CLAIMED[(4, 2)][0])
    ck('P_4_2_equals_160_over_499', Fraction(c42['A'], c42['B']) == CLAIMED_P[(4, 2)],
       '%.9f' % float(Fraction(c42['A'], c42['B'])))
    print()

    # ---------------------------------------------------------------- Step 7
    print('=== Step 7: the model-sensitivity disclosure, itself checked')
    riv = {k: census(3, k, 'rival') for k in (2, 3, 4)}
    for k in (2, 3, 4):
        ck('rival_model_P_3_%d_equals_%s' % (k, str(CLAIMED_RIVAL_P[k]).replace('/', '_over_')),
           Fraction(riv[k]['A'], riv[k]['B']) == CLAIMED_RIVAL_P[k],
           'A_u=%d B_u=%d, %.6f' % (riv[k]['A'], riv[k]['B'],
                                    float(Fraction(riv[k]['A'], riv[k]['B']))))
    ck('rival_model_shows_NO_violation_at_n3_for_k_2_to_4',
       CLAIMED_RIVAL_P[2] < CLAIMED_RIVAL_P[3] < CLAIMED_RIVAL_P[4],
       '34*486 = %d < 266*65 = %d' % (34 * 486, 266 * 65))
    # the source's own displayed identity separates the two readings
    n, k = 3, 2
    rows = count_vectors(n, k)
    wts = [multinomial(k, r) for r in rows]
    tot_o = sum(wts) ** n
    hit_o = sum(w for r, w in zip(rows, wts) if r[0] == 0) ** n
    ck('indegree_zero_probability_in_the_ordered_model_is_64_over_729',
       Fraction(hit_o, tot_o) == Fraction(2, 3) ** (k * n),
       'Pr(no edge lands on v) = (1-1/n)^{kn} = 64/729 in the ordered model')
    hit_r = len([r for r in rows if r[0] == 0]) ** n
    ck('the_same_probability_in_the_rival_model_is_1_over_8_not_64_over_729',
       Fraction(hit_r, len(rows) ** n) == Fraction(1, 8),
       'the two readings give different values for this same event')
    print()

    # ---------------------------------------------------------------- Step 8
    print('=== Step 8: three exhibited 3-vertex states, and the TS-but-not-SCC weights')
    G_cerny = [(0, 2, 0), (0, 1, 1), (1, 0, 1)]        # 0->{1,1}  1->{1,2}  2->{2,0}
    G_notts = [(0, 1, 1), (1, 1, 0), (1, 0, 1)]        # 0->{1,2}  1->{0,1}  2->{0,2}
    G_tsnoscc = [(2, 0, 0), (1, 0, 1), (0, 2, 0)]      # 0->{0,0}  1->{0,2}  2->{1,1}
    ck('cerny_underlying_2_out_graph_is_scc_and_totally_simple',
       is_strongly_connected(G_cerny, 3) and is_totally_simple(G_cerny, nts3),
       'one graph tested: 0->{1,1}, 1->{1,2}, 2->{2,0}')
    ck('printed_state_is_scc_but_NOT_totally_simple',
       is_strongly_connected(G_notts, 3) and not is_totally_simple(G_notts, nts3),
       'partition {{1,2},{0}} is lumpable since a_10 = a_20 = 1')
    ck('printed_state_is_totally_simple_but_NOT_strongly_connected',
       is_totally_simple(G_tsnoscc, nts3) and not is_strongly_connected(G_tsnoscc, 3),
       'one exhibited 3-vertex state; total simplicity does not imply strong connectivity')
    ck('weight_of_totally_simple_and_not_scc_is_36_at_k2_and_666_at_k3',
       cen[2]['ts_not_scc'] == 36 and cen[3]['ts_not_scc'] == 666,
       'nonzero in every cell computed, so Pr(TS|SCC) is NOT Pr(TS)/Pr(SCC)')
    ck('weight_of_totally_simple_and_not_scc_is_8628_at_n3_k4_and_1584_at_n4_k2',
       cen[4]['ts_not_scc'] == 8628 and c42['ts_not_scc'] == 1584,
       'out of 3^12 = 531441 and 4^8 = 65536 respectively')
    ok2 = True
    for kk in (2, 3):
        r2 = count_vectors(2, kk)
        for c in itertools.product(range(len(r2)), repeat=2):
            if not is_totally_simple([r2[i] for i in c], nontrivial_partitions(2)):
                ok2 = False
    ck('every_2_vertex_k_out_graph_is_totally_simple_at_k_2_and_3_control',
       ok2 and nontrivial_partitions(2) == [],
       'Bell(2)-2 = 0 nontrivial partitions, so P(2,2) = P(2,3) = 1')
    print()

    note('SCOPE -- what this program does NOT do. It computes finitely many exact census '
         'values at n = 2,3,4 and one non-monotone k-step among them. It does not read or '
         'parse the source, so the wording of any clause and its location are transcribed '
         'by hand here. Asymptotic statements (Pr -> 1 as n -> infinity, and any '
         '1 - exp(-c_k n + beta_k) + o(1) expansion with c_k increasing) are untouched: no '
         'finite computation can settle them, and nothing here claims to.')
    note('SCOPE -- cells NOT enumerated here. Exhaustively enumerated above: n=3 at k=1,2,3,4 '
         '(and n=3 at k=2,3 a second time over raw ordered tuples), n=4 at k=2, n=2 at k=2,3. '
         'The closed forms are exercised at n=3 for k=2..%d and are n=3 ONLY. Cells at n>=5, '
         'and n=4 at k>=3, are NOT recomputed here; the proposed repair "for fixed n >= 4" is '
         'therefore CONSISTENT WITH, NOT VERIFIED BY, this program.' % KMAX_CLOSED)
    note('SCOPE -- the model is a READING of the source. Both readings are computed above and '
         'the counterexample exists in the ordered/with-replacement one only. The residual '
         'risk in this result is that reading, not the arithmetic.')
    note('SCOPE -- NOT recomputed or looked up here: any sampled table at larger n, and any '
         'attribution of the numbers B(n,2). The values 296 and 20958 are transcribed by '
         'hand in the CLAIMED block above and only compared against this census.')

    print()
    if _n_fail:
        print('VERDICT: %d CHECK(S) FAILED out of %d' % (_n_fail, _n_pass + _n_fail))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _n_pass)
    return 0


if __name__ == '__main__':
    sys.exit(main())

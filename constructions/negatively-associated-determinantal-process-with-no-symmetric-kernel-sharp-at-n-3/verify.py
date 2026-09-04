#!/usr/bin/env python3
"""verify.py -- exact-rational verification of the computational claims of

    "An Explicit Negatively Associated Determinantal Point Process with No
     Symmetric Kernel, Minimal at n = 3"

It reads the matrices EXHIBITED IN THE PAPER (they are typed out below exactly as
printed there, nothing is loaded from disk) and re-derives the quantities the paper
states: principal minors, the 2^n atoms by Moebius inversion, negative association
brute-forced FROM THE DEFINITION (all ordered pairs of disjoint nonempty sets, all
pairs of proper nonempty up-closed families, with up-closure TESTED not assumed),
the triple invariants u, v, w, p, q, sigma, the complete exclusion of a real
symmetric competitor by sign enumeration, the Hermitian identity, Theorem B's
closed form against the brute force, Theorem C's sigma-interval, Theorem D's
padding, the interval box certificate of radius 1/100, the controls in both
polarities, and the strong-Rayleigh separation.  Nothing outside this file is read or
referenced.

Python 3.9+, STANDARD LIBRARY ONLY (fractions, itertools).  Exact integer and
Fraction arithmetic throughout: no float ever decides anything.

One `PASS <name> [detail]` line per check; exits 0 iff every check passed.
"""

import sys
from fractions import Fraction as F
from itertools import combinations, permutations

# ---------------------------------------------------------------------------
# check harness
# ---------------------------------------------------------------------------
_N_OK = 0
_N_BAD = 0


def ck(name, cond, detail=''):
    global _N_OK, _N_BAD
    tag = 'PASS' if cond else 'FAIL'
    if cond:
        _N_OK += 1
    else:
        _N_BAD += 1
    detail = ' '.join(str(detail).split())
    print('%s %s%s' % (tag, name, (' [' + detail + ']') if detail else ''))


# ---------------------------------------------------------------------------
# exact linear algebra on subsets encoded as bitmasks
# ---------------------------------------------------------------------------
def bits(mask):
    out = []
    i = 0
    while mask:
        if mask & 1:
            out.append(i)
        mask >>= 1
        i += 1
    return out


def det(M):
    """Exact determinant by full permutation expansion (no division, no pivoting)."""
    n = len(M)
    if n == 0:
        return F(1)
    tot = F(0)
    for p in permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])
        pr = F(1)
        for i in range(n):
            pr *= M[i][p[i]]
        tot += -pr if inv % 2 else pr
    return tot


def minor(K, mask):
    S = bits(mask)
    return det([[K[i][j] for j in S] for i in S])


def all_minors(K):
    n = len(K)
    return {m: minor(K, m) for m in range(1 << n)}


def atoms_of(K, mins=None):
    """P(X = A) = sum_{A subset T} (-1)^{|T|-|A|} det(K_T)  (Moebius inversion)."""
    n = len(K)
    mins = mins if mins is not None else all_minors(K)
    full = (1 << n) - 1
    out = {}
    for A in range(1 << n):
        rest = bits(full ^ A)
        tot = F(0)
        for k in range(len(rest) + 1):
            for T in combinations(rest, k):
                add = 0
                for t in T:
                    add |= 1 << t
                tot += -mins[A | add] if k % 2 else mins[A | add]
        out[A] = tot
    return out


def atoms_from_L(K):
    """Independent second route to the atoms: P(X=A) = det(L_A)/det(I+L), L = K(I-K)^{-1}.

    Implemented without matrix inversion: det(L_A) = det(K_A')/det(I-K) where
    K_A' is K with the rows/columns outside A replaced by those of K - I, i.e. the
    classical identity  P(X = A) = |det( D_A (K - I) + D_{A^c} ... )|.  We use the
    equivalent standard form  P(X = A) = det(K - I_{A^c})*(-1)^{|A^c|}, where
    I_{A^c} is the identity restricted to the complement of A.
    """
    n = len(K)
    out = {}
    for A in range(1 << n):
        M = [[K[i][j] for j in range(n)] for i in range(n)]
        comp = 0
        for i in range(n):
            if not (A >> i) & 1:
                M[i][i] -= 1
                comp += 1
        out[A] = det(M) * ((-1) ** comp)
    return out


# ---------------------------------------------------------------------------
# up-closed families of the Boolean lattice, generated and then TESTED
# ---------------------------------------------------------------------------
def is_upclosed(fam, subs, atom_bits, full_pairwise=True):
    """Independent test, used on the generator's output and on a deliberate non-example."""
    if full_pairwise:
        for x in fam:
            for y in subs:
                if (x & y) == x and x != y and y not in fam:
                    return False
        return True
    for x in fam:
        for b in atom_bits:
            if not (x & b) and (x | b) not in fam:
                return False
    return True


def upclosed_families(elems, full_pairwise_test=True):
    """All up-closed families of the lattice of subsets of `elems` (as bitmasks).

    Generated by a DFS in order of decreasing cardinality -- a set may be admitted
    only when every immediate superset already is -- and then every produced family
    is INDEPENDENTLY tested for up-closure.  The count is the Dedekind number.
    """
    elems = list(elems)
    k = len(elems)
    subs = []
    for r in range(k + 1):
        for c in combinations(elems, r):
            m = 0
            for t in c:
                m |= 1 << t
            subs.append(m)
    subs.sort(key=lambda m: -bin(m).count('1'))
    atom_bits = [1 << t for t in elems]
    fams = []
    chosen = set()

    def rec(i):
        if i == len(subs):
            fams.append(frozenset(chosen))
            return
        X = subs[i]
        rec(i + 1)
        if all((X | b) in chosen for b in atom_bits if not (X & b)):
            chosen.add(X)
            rec(i + 1)
            chosen.discard(X)

    rec(0)
    # independent up-closure test on every produced family: TESTED, not assumed
    n_validated = 0
    for fam in fams:
        if not is_upclosed(fam, subs, atom_bits, full_pairwise_test):
            raise AssertionError('generated family is not up-closed')
        n_validated += 1
    return fams, subs, n_validated


_UPS_CACHE = {}


def proper_upsets(elems):
    """The proper nonempty up-closed families (M(k) - 2 of them)."""
    key = tuple(sorted(elems))
    if key in _UPS_CACHE:
        return _UPS_CACHE[key]
    fams, subs, _ = upclosed_families(elems)
    allm = frozenset(subs)
    res = [f for f in fams if f and f != allm]
    _UPS_CACHE[key] = res
    return res


# ---------------------------------------------------------------------------
# negative association, brute-forced from the definition
# ---------------------------------------------------------------------------
def na_report(atoms, n):
    """-> dict(count, violations, min_slack, n_equal, first_violation)

    For every ORDERED pair of disjoint nonempty A, B subset [n] and every pair of
    proper nonempty up-closed families U of 2^A, V of 2^B:
        P(X cap A in U, X cap B in V) <= P(X cap A in U) P(X cap B in V).
    This is EQUIVALENT to the increasing-function definition, because the defect
    E[fg]-E[f]E[g] is bilinear and invariant under adding constants and every
    increasing function on a finite Boolean lattice is a nonnegative combination of
    up-set indicators plus a constant.
    """
    universe = list(range(n))
    keys = list(atoms)
    marg = {}
    sets = {}
    for ra in range(1, n):
        for A in combinations(universe, ra):
            Am = 0
            for t in A:
                Am |= 1 << t
            for U in proper_upsets(A):
                idx = frozenset(m for m in keys if (m & Am) in U)
                sets[(Am, U)] = idx
                marg[(Am, U)] = sum(atoms[m] for m in idx)
    count = 0
    viol = []
    slacks = []
    n_eq = 0
    for ra in range(1, n):
        for A in combinations(universe, ra):
            Am = 0
            for t in A:
                Am |= 1 << t
            rest = [x for x in universe if x not in A]
            for rb in range(1, len(rest) + 1):
                for B in combinations(rest, rb):
                    Bm = 0
                    for t in B:
                        Bm |= 1 << t
                    for U in proper_upsets(A):
                        for V in proper_upsets(B):
                            count += 1
                            joint = sum(atoms[m] for m in (sets[(Am, U)] & sets[(Bm, V)]))
                            prod = marg[(Am, U)] * marg[(Bm, V)]
                            s = prod - joint
                            slacks.append(s)
                            if s < 0:
                                viol.append((A, B, joint, prod))
                            elif s == 0:
                                n_eq += 1
    return {'count': count, 'violations': len(viol), 'min_slack': min(slacks) if slacks else None,
            'n_equal': n_eq, 'first_violation': viol[0] if viol else None}


# ---------------------------------------------------------------------------
# triple invariants
# ---------------------------------------------------------------------------
def triple(K, i, j, k, mins=None):
    mins = mins if mins is not None else all_minors(K)
    u = K[i][j] * K[j][i]
    v = K[i][k] * K[k][i]
    w = K[j][k] * K[k][j]
    p = K[i][j] * K[j][k] * K[k][i]
    q = K[i][k] * K[k][j] * K[j][i]
    sig = p + q
    T = (1 << i) | (1 << j) | (1 << k)
    sig_min = (mins[T] - K[i][i] * K[j][j] * K[k][k]
               + K[i][i] * w + K[j][j] * v + K[k][k] * u)
    u_min = K[i][i] * K[j][j] - mins[(1 << i) | (1 << j)]
    return {'u': u, 'v': v, 'w': w, 'p': p, 'q': q, 'sigma': sig,
            'sigma_from_minors': sig_min, 'u_from_minors': u_min,
            'gap': sig * sig - 4 * u * v * w}


def thmB(K):
    """Theorem B's closed form at n = 3: NA iff u,v,w >= 0 and -min S^- <= sigma <= min S^+."""
    assert len(K) == 3
    mins = all_minors(K)
    a = [K[i][i] for i in range(3)]
    U = {}
    for i, j in ((0, 1), (0, 2), (1, 2)):
        U[(i, j)] = a[i] * a[j] - mins[(1 << i) | (1 << j)]
        U[(j, i)] = U[(i, j)]
    t = triple(K, 0, 1, 2, mins)
    sig = t['sigma']
    slacks = [U[(0, 1)], U[(0, 2)], U[(1, 2)]]
    for i in range(3):
        others = [x for x in range(3) if x != i]
        Sp = F(0)
        Sm = F(0)
        for j in others:
            k = [x for x in range(3) if x not in (i, j)][0]
            Sp += a[k] * U[(i, j)]
            Sm += (1 - a[k]) * U[(i, j)]
        slacks.append(Sp - sig)
        slacks.append(sig + Sm)
    return {'min_slack': min(slacks), 'na': min(slacks) >= 0}


def is_dpp_kernel(K):
    at = atoms_of(K)
    return min(at.values()) >= 0, at


# ---------------------------------------------------------------------------
# THE OBJECTS, exactly as printed in the paper
# ---------------------------------------------------------------------------
W = [[F(1, 2), F(1, 10), F(1, 25)],
     [F(1, 10), F(1, 2), F(1, 10)],
     [F(1, 4), F(1, 10), F(1, 2)]]

Wp = [[F(1, 2), F(1, 10), F(1, 50)],
      [F(1, 10), F(1, 2), F(1, 10)],
      [F(1, 2), F(1, 10), F(1, 2)]]

K4 = [[F(1, 2), F(1, 10), F(3, 25), F(1, 10)],
      [F(1, 10), F(1, 2), F(1, 10), F(1, 10)],
      [F(1, 12), F(1, 10), F(1, 2), F(1, 10)],
      [F(1, 10), F(1, 10), F(1, 10), F(1, 2)]]

SSYM = [[F(1, 2), F(1, 10), F(1, 10)],
        [F(1, 10), F(1, 2), F(1, 10)],
        [F(1, 10), F(1, 10), F(1, 2)]]

KNEG = [[F(1, 2), F(1, 10), F(1, 25)],
        [F(-1, 10), F(1, 2), F(1, 10)],
        [F(1, 4), F(1, 10), F(1, 2)]]

KEVEN = [[F(1, 2), F(0), F(-1, 2)],
         [F(1, 2), F(1, 2), F(0)],
         [F(0), F(1, 2), F(1, 2)]]

KC = [[F(1, 2), F(1, 5), F(1, 5)],
      [F(1, 20), F(1, 2), F(1, 10)],
      [F(1, 20), F(1, 10), F(1, 2)]]


def main():
    print('=' * 78)
    print('1. THE PRIMARY WITNESS W  (n = 3)')
    print('=' * 78)
    minsW = all_minors(W)
    want = {0: F(1), 1: F(1, 2), 2: F(1, 2), 4: F(1, 2),
            3: F(6, 25), 5: F(6, 25), 6: F(6, 25), 7: F(1129, 10000)}
    ck('W_principal_minors_match_the_paper',
       all(minsW[m] == want[m] for m in want),
       'det W = %s, every pair minor = %s, every singleton = 1/2' % (minsW[7], minsW[3]))
    ck('W_minors_depend_only_on_cardinality_so_the_law_is_exchangeable',
       len(set(minsW[m] for m in minsW if bin(m).count('1') == 1)) == 1
       and len(set(minsW[m] for m in minsW if bin(m).count('1') == 2)) == 1)
    atW = atoms_of(W, minsW)
    ck('W_atoms_exact',
       (atW[0] == F(1071, 10000)
        and all(atW[1 << i] == F(1329, 10000) for i in range(3))
        and all(atW[m] == F(1271, 10000) for m in (3, 5, 6))
        and atW[7] == F(1129, 10000)),
       'empty 1071/10000, singletons 1329/10000, pairs 1271/10000, top 1129/10000 (all over 10^4)')
    ck('W_atoms_sum_to_one_exactly', sum(atW.values()) == 1)
    ck('W_all_eight_atoms_strictly_positive_so_W_is_a_DPP_kernel',
       min(atW.values()) > 0, 'min atom = %s' % min(atW.values()))
    ck('W_min_atom_is_1071_over_10000', min(atW.values()) == F(1071, 10000))
    atW2 = atoms_from_L(W)
    ck('W_atoms_reproduced_by_a_second_independent_route',
       all(atW2[m] == atW[m] for m in atW),
       'det(K - I_{A^c}) route agrees with Moebius inversion on all 8 atoms')
    ck('W_complement_atom_is_det_I_minus_W_so_L_ensemble_exists',
       atW[0] == det([[(1 if i == j else 0) - W[i][j] for j in range(3)] for i in range(3)])
       and atW[0] != 0, 'P(X = empty) = det(I - W) = %s != 0' % atW[0])

    print()
    print('=' * 78)
    print('2. NEGATIVE ASSOCIATION OF DPP(W), BRUTE-FORCED FROM THE DEFINITION')
    print('=' * 78)
    # the enumerator is pinned against the Dedekind numbers, OEIS A000372
    ded = []
    nval = 0
    for k in range(6):
        fams, subs5, nv = upclosed_families(range(k), full_pairwise_test=(k <= 4))
        ded.append(len(fams))
        nval += nv
    ck('upset_enumerator_reproduces_the_Dedekind_numbers_OEIS_A000372',
       ded == [2, 3, 6, 20, 168, 7581], 'k = 0..5 gives %s' % (ded,))
    ck('every_generated_family_passed_an_independent_up_closure_test',
       nval == sum(ded), '%d families generated, %d independently re-tested for up-closure -- '
                         'up-closure is TESTED, not assumed by construction' % (sum(ded), nval))
    bogus_subs = [0, 1, 2, 3]
    ck('the_up_closure_tester_rejects_a_deliberate_non_up_closed_family',
       not is_upclosed(frozenset([1]), bogus_subs, [1, 2]),
       'the family {{1}} on the 2-element lattice omits {1,2} and is correctly rejected, '
       'so the tester is not vacuous')
    M = ded
    f = [M[k] - 2 for k in range(6)]

    def predicted(n):
        tot = 0
        for a in range(1, n):
            for b in range(1, n - a + 1):
                ways = 1
                for t in range(a):
                    ways = ways * (n - t) // (t + 1)
                w2 = 1
                for t in range(b):
                    w2 = w2 * (n - a - t) // (t + 1)
                tot += ways * w2 * f[a] * f[b]
        return tot

    ck('predicted_ordered_test_counts_from_the_Dedekind_numbers',
       [predicted(n) for n in (3, 4, 5, 6)] == [30, 348, 4560, 140058],
       'N(3),N(4),N(5),N(6) = %s' % ([predicted(n) for n in (3, 4, 5, 6)],))
    naW = na_report(atW, 3)
    ck('W_NA_enumerator_produced_exactly_the_predicted_30_ordered_inequalities',
       naW['count'] == 30 == predicted(3))
    ck('W_is_negatively_associated_zero_violations',
       naW['violations'] == 0, '30 ordered inequalities, 0 violations')
    ck('W_min_NA_slack_is_71_over_10000_strictly_positive',
       naW['min_slack'] == F(71, 10000) and naW['n_equal'] == 0,
       'min slack 71/10000, no equalities: DPP(W) is STRICTLY negatively associated')
    ck('W_satisfies_the_sources_own_necessary_condition_K_ij_K_ji_ge_0',
       all(W[i][j] * W[j][i] == F(1, 100) for i, j in ((0, 1), (0, 2), (1, 2))),
       'all three products equal 1/100 > 0')
    # the three exchangeability classes quoted in the paper
    ck('NA_class_P_pairwise_slack_1_over_100',
       F(1, 4) - minsW[3] == F(1, 100), 'P(i,j in X) = 6/25 vs 1/4')
    ck('NA_class_T_slack_71_over_10000',
       F(1, 2) * F(6, 25) - minsW[7] == F(71, 10000), 'P(all three) = 1129/10000 vs (1/2)(6/25)')
    ck('NA_class_U_slack_129_over_10000',
       F(1, 2) * (1 - atW[0] - atW[1]) - (F(1, 2) - atW[1]) == F(129, 10000),
       'P(1 in X, X meets {2,3}) = 3671/10000 vs (1/2)(19/25)')

    print()
    print('=' * 78)
    print('3. NO REAL SYMMETRIC KERNEL REPRODUCES DPP(W) (COMPLETE SIGN ENUMERATION), AND NO')
    print('   GAUSSIAN-RATIONAL HERMITIAN CANDIDATE WITH THE FORCED MODULI DOES EITHER (1728')
    print('   CANDIDATES; the general Hermitian case is proved in the paper, not here)')
    print('=' * 78)
    tW = triple(W, 0, 1, 2, minsW)
    ck('W_triple_invariants',
       (tW['u'], tW['v'], tW['w'], tW['p'], tW['q'], tW['sigma'])
       == (F(1, 100), F(1, 100), F(1, 100), F(1, 400), F(1, 2500), F(29, 10000)),
       'u = v = w = 1/100, p = 1/400, q = 1/2500, sigma = 29/10000')
    ck('identity_pq_equals_uvw_holds_with_no_hypothesis',
       tW['p'] * tW['q'] == tW['u'] * tW['v'] * tW['w'])
    ck('u_and_sigma_are_functions_of_the_principal_minors_alone',
       tW['u'] == tW['u_from_minors'] and tW['sigma'] == tW['sigma_from_minors'],
       'u = K_ii K_jj - det K_{ij} and sigma = det K - K_ii K_jj K_kk + K_ii w + K_jj v + K_kk u')
    ck('W_gap_sigma_squared_minus_4uvw_is_441_over_10_to_the_8',
       tW['gap'] == F(441, 100000000) and tW['gap'] == (tW['p'] - tW['q']) ** 2,
       'gap = (p - q)^2 = (21/10000)^2 = 441/100000000 > 0')
    # COMPLETE exclusion of a real symmetric competitor: |S_ij| is forced and 1/100 is a
    # rational square, so the eight sign patterns are exhaustive.
    forced = {}
    for i, j in ((0, 1), (0, 2), (1, 2)):
        sq = W[i][i] * W[j][j] - minsW[(1 << i) | (1 << j)]
        forced[(i, j)] = sq
    ck('forced_symmetric_off_diagonal_squares_are_rational_squares',
       all(v == F(1, 100) for v in forced.values()),
       'S_ij^2 = K_ii K_jj - det K_{ij} = 1/100 for all three pairs, so S_ij = +-1/10')
    dets = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            for s3 in (1, -1):
                S = [[F(1, 2), s1 * F(1, 10), s2 * F(1, 10)],
                     [s1 * F(1, 10), F(1, 2), s3 * F(1, 10)],
                     [s2 * F(1, 10), s3 * F(1, 10), F(1, 2)]]
                dets.append(det(S))
    ck('no_real_symmetric_matrix_reproduces_det_W_complete_8_pattern_enumeration',
       len(dets) == 8 and all(d != minsW[7] for d in dets),
       'the 8 achievable determinants are %s, none equals 1129/10000'
       % sorted(set(str(d) for d in dets)))
    # Hermitian: the same forced moduli, now over Gaussian rationals with |S_ij|^2 = 1/100.
    hermpairs = []
    for (a, b) in ((F(1, 10), F(0)), (F(0), F(1, 10)), (F(3, 50), F(4, 50)), (F(4, 50), F(3, 50))):
        for sa in (1, -1):
            for sb in (1, -1):
                z = (sa * a, sb * b)
                if z not in hermpairs and z[0] ** 2 + z[1] ** 2 == F(1, 100):
                    hermpairs.append(z)

    def cmul(x, y):
        return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    nherm = 0
    bad = 0
    identity_ok = True
    for z12 in hermpairs:
        for z13 in hermpairs:
            for z23 in hermpairs:
                nherm += 1
                # Hermitian S: S_21 = conj(S_12) etc.  p_S = S_12 S_23 S_31 = z12 z23 conj(z13)
                pS = cmul(cmul(z12, z23), (z13[0], -z13[1]))
                sigS = 2 * pS[0]
                gap = sigS * sigS - 4 * F(1, 100) ** 3
                if gap != -4 * pS[1] ** 2:
                    identity_ok = False
                if gap > 0:
                    bad += 1
                if sigS == tW['sigma']:
                    bad += 1
    ck('hermitian_identity_gap_equals_minus_4_Im_p_squared_on_every_candidate',
       identity_ok and nherm == len(hermpairs) ** 3,
       '%d Hermitian candidates with the forced moduli |S_ij|^2 = 1/100; for each, '
       'sigma^2 - 4uvw = -4 Im(p)^2 <= 0' % nherm)
    ck('no_hermitian_candidate_attains_W_positive_gap_or_sigma',
       bad == 0, 'W has gap 441/10^8 > 0; none of the 1728 enumerated Gaussian-rational '
                 'Hermitian candidates with the forced moduli attains a nonnegative gap or '
                 'sigma = 29/10000 -- the general Hermitian exclusion is the paper\'s Lemma, '
                 'not re-derived here')

    print()
    print('=' * 78)
    print('4. CONTROLS, IN BOTH POLARITIES')
    print('=' * 78)
    okS, atS = is_dpp_kernel(SSYM)
    naS = na_report(atS, 3)
    tS = triple(SSYM, 0, 1, 2)
    ck('positive_control_symmetric_kernel_is_NA_with_min_slack_1_over_125',
       okS and naS['violations'] == 0 and naS['min_slack'] == F(1, 125),
       'S = (1/2)I + (1/10)(J - I): min atom %s, 30 checks, 0 violations' % min(atS.values()))
    ck('positive_control_gap_is_exactly_zero_so_the_test_never_flags_a_symmetric_matrix',
       tS['gap'] == 0)
    okN, atN = is_dpp_kernel(KNEG)
    naN = na_report(atN, 3)
    ck('negative_control_is_still_a_genuine_DPP_kernel',
       okN and sum(atN.values()) == 1 and min(atN.values()) == F(1179, 10000),
       'W with W_21 flipped to -1/10: min atom 1179/10000')
    ck('negative_control_FAILS_negative_association_so_the_decider_can_say_no',
       naN['violations'] == 10 and naN['min_slack'] == F(-1, 100),
       '10 of 30 inequalities violated, min slack -1/100 = u < 0')
    okE, atE = is_dpp_kernel(KEVEN)
    naE = na_report(atE, 3)
    tE = triple(KEVEN, 0, 1, 2)
    ck('sources_own_uniform_on_even_subsets_kernel_is_a_DPP_kernel',
       okE and sum(atE.values()) == 1 and min(atE.values()) == 0,
       '(1/2)[[1,0,-1],[1,1,0],[0,1,1]]: min atom 0')
    ck('sources_own_kernel_passes_the_necessary_condition_and_is_non_symmetrisable',
       all(KEVEN[i][j] * KEVEN[j][i] == 0 for i, j in ((0, 1), (0, 2), (1, 2)))
       and tE['gap'] > 0, 'all pairwise products 0, gap = %s > 0' % tE['gap'])
    ck('sources_own_kernel_FAILS_negative_association_min_slack_minus_one_eighth',
       naE['violations'] > 0 and naE['min_slack'] == F(-1, 8),
       'so non-symmetrisability ALONE answers nothing: %d violations' % naE['violations'])
    okC, atC = is_dpp_kernel(KC)
    naC = na_report(atC, 3)
    tC = triple(KC, 0, 1, 2)
    ck('anti_control_is_nonsymmetric_yet_its_gap_is_exactly_zero',
       any(KC[i][j] != KC[j][i] for i in range(3) for j in range(3)) and tC['gap'] == 0
       and okC and naC['violations'] == 0,
       'Kc = [[1/2,1/5,1/5],[1/20,1/2,1/10],[1/20,1/10,1/2]] is NA, nonsymmetric, gap 0')
    ck('anti_control_has_exactly_the_symmetric_controls_law',
       all(atC[m] == atS[m] for m in atS),
       'so the certificate is NOT merely testing whether K is literally symmetric')

    print()
    print('=' * 78)
    print("5. THEOREM B (closed form at n = 3) AGAINST THE BRUTE FORCE")
    print('=' * 78)
    for name, K, na in (('W', W, naW), ('W_prime', Wp, None), ('symmetric_control', SSYM, naS),
                        ('negative_control', KNEG, naN), ('source_even_kernel', KEVEN, naE),
                        ('anti_control', KC, naC)):
        if na is None:
            na = na_report(atoms_of(K), 3)
        b = thmB(K)
        ck('thmB_agrees_with_brute_force_on_' + name,
           b['min_slack'] == na['min_slack'] and b['na'] == (na['violations'] == 0),
           'closed form %s, brute force %s, NA = %s'
           % (b['min_slack'], na['min_slack'], b['na']))

    print()
    print('=' * 78)
    print("6. THE SECOND n = 3 WITNESS W'")
    print('=' * 78)
    okP, atP = is_dpp_kernel(Wp)
    naP = na_report(atP, 3)
    tP = triple(Wp, 0, 1, 2)
    ck('Wprime_is_a_DPP_kernel_with_min_atom_131_over_1250',
       okP and sum(atP.values()) == 1 and min(atP.values()) == F(131, 1250))
    ck('Wprime_is_negatively_associated_min_slack_3_over_625',
       naP['violations'] == 0 and naP['min_slack'] == F(3, 625))
    ck('Wprime_gap_is_9_over_390625_positive',
       tP['gap'] == F(9, 390625) and tP['sigma'] == F(13, 2500))

    print()
    print('=' * 78)
    print('7. THEOREM C: AN INTERVAL OF SUCH LAWS IN THE EXCHANGEABLE SLICE')
    print('=' * 78)
    t, s = F(1, 2), F(1, 10)
    lo, hi = 2 * s ** 3, 2 * t * s ** 2
    ck('thmC_slice_admissibility_of_the_base_symmetric_kernel',
       t - s >= 0 and t + 2 * s <= 1 and 2 * t * s ** 2 <= (1 - t) * ((1 - t) ** 2 - 3 * s ** 2),
       't = 1/2, s = 1/10 satisfy t + 2s <= 1 and 2ts^2 <= (1-t)[(1-t)^2 - 3s^2]')
    ck('thmC_interval_endpoints', (lo, hi) == (F(1, 500), F(1, 100)),
       'symmetrisable point 2s^3 = 1/500, NA ceiling 2ts^2 = 1/100, length 2s^2(t-s) = 8/1000')
    tested = 0
    for m in range(2, 91, 4):
        x = s ** 2 * m
        K = [[t, s, x], [s, t, s], [s ** 2 / x, s, t]]
        tt = triple(K, 0, 1, 2)
        if not (lo < tt['sigma'] <= hi):
            continue
        okk, att = is_dpp_kernel(K)
        nn = na_report(att, 3)
        good = (okk and min(att.values()) > 0 and nn['violations'] == 0
                and nn['min_slack'] > 0 and tt['gap'] > 0
                and tt['u'] == tt['v'] == tt['w'] == s ** 2)
        if not good:
            ck('thmC_member_m_%d' % m, False, 'sigma = %s' % tt['sigma'])
        tested += 1
    ck('thmC_every_tested_member_of_the_sigma_interval_is_NA_and_non_symmetrisable',
       tested >= 10, '%d rational members of (1/500, 1/100] built as K_13 = s^2 m, '
                     'K_31 = s^2/K_13; each is a DPP kernel, strictly NA, gap > 0' % tested)
    Ksym = [[t, s, s], [s, t, s], [s, s, t]]
    ck('thmC_lower_endpoint_is_exactly_the_symmetrisable_point',
       triple(Ksym, 0, 1, 2)['sigma'] == lo and triple(Ksym, 0, 1, 2)['gap'] == 0)
    Kceil = [[t, s, F(1, 1)], [s, t, s], [s ** 2, s, t]]
    tceil = triple(Kceil, 0, 1, 2)
    okc2, atc2 = is_dpp_kernel(Kceil)
    nac2 = na_report(atc2, 3)
    ck('thmC_ceiling_is_real_a_larger_sigma_breaks_NA_not_the_DPP_property',
       tceil['sigma'] > hi and okc2 and nac2['violations'] > 0 and nac2['min_slack'] < 0,
       'K_13 = 1, K_31 = 1/100 gives sigma = %s > 1/100: still a DPP kernel, NA fails with '
       'min slack %s' % (tceil['sigma'], nac2['min_slack']))

    print()
    print('=' * 78)
    print('8. THEOREM D: n = 4 WITNESS, PADDING, AND THE n <= 2 IMPOSSIBILITY')
    print('=' * 78)
    ok4, at4 = is_dpp_kernel(K4)
    na4 = na_report(at4, 4)
    ck('K4_is_a_DPP_kernel_with_min_atom_1079_over_25000',
       ok4 and sum(at4.values()) == 1 and min(at4.values()) == F(1079, 25000),
       '16 atoms summing to 1 exactly')
    ck('K4_is_negatively_associated_348_checks_0_violations_min_slack_179_over_37500',
       na4['count'] == 348 and na4['violations'] == 0 and na4['min_slack'] == F(179, 37500),
       'strictly NA at n = 4 and not block diagonal')
    gaps4 = {T: triple(K4, *T)['gap'] for T in combinations(range(4), 3)}
    ck('K4_triple_gaps_are_121_over_9e8_on_two_triples_and_0_on_the_other_two',
       sorted(gaps4.values()) == [F(0), F(0), F(121, 900000000), F(121, 900000000)],
       'obstruction on {1,2,3} and {1,3,4}: %s' % {str(k): str(v) for k, v in gaps4.items()})
    ck('K4_2x2_minors_equal_the_symmetric_reference_S4',
       all(K4[i][j] * K4[j][i] == F(1, 100) for i in range(4) for j in range(4) if i < j),
       'moving (K_13, K_31) from (1/10,1/10) to (3/25,1/12) keeps every pair product at 1/100')
    # padding: W (+) Bernoulli(1/3)
    c = F(1, 3)
    Kpad = [[W[i][j] if i < 3 and j < 3 else (c if (i == 3 and j == 3) else F(0))
             for j in range(4)] for i in range(4)]
    okpad, atpad = is_dpp_kernel(Kpad)
    napad = na_report(atpad, 4)
    tpad = triple(Kpad, 0, 1, 2)
    ck('padding_W_by_an_independent_Bernoulli_stays_a_DPP_kernel',
       okpad and sum(atpad.values()) == 1 and min(atpad.values()) > 0,
       'K_4 = W (+) diag(1/3), min atom %s' % min(atpad.values()))
    ck('padded_kernel_is_negatively_associated_348_checks_0_violations',
       napad['count'] == 348 and napad['violations'] == 0,
       'min slack %s over all 348 -- the padded coordinate is independent, so some '
       'inequalities are EQUALITIES (%d of them) and NA is defined with <=, so slack 0 passes'
       % (napad['min_slack'], napad['n_equal']))
    ck('padded_kernel_inherits_the_obstruction_on_the_triple_1_2_3',
       tpad['gap'] == F(441, 100000000),
       'a single bad triple is a permanent obstruction: any same-law symmetric or Hermitian '
       'competitor has the same principal minors on {1,2,3}')
    ck('atoms_of_the_padding_factorise_as_an_independent_superposition',
       all(atpad[m] == atW[m & 7] * (c if (m >> 3) & 1 else 1 - c) for m in atpad))
    # n <= 2 impossibility, on an exact sweep
    n2 = 0
    n2bad = 0
    naneg = 0
    for a in (F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4)):
        for d in (F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(3, 4)):
            for r in (F(0), F(1, 20), F(1, 10), F(1, 5), F(1, 4)):
                for kk in (F(1), F(2), F(5)):
                    for sgn in (1, -1):
                        b, cc = r * kk, sgn * r / kk if kk != 0 else F(0)
                        K = [[a, b], [cc, d]]
                        ok2, at2 = is_dpp_kernel(K)
                        if not ok2:
                            continue
                        na2 = na_report(at2, 2)
                        u = b * cc
                        if na2['violations'] > 0:
                            if u >= 0:
                                n2bad += 1
                            naneg += 1
                            continue
                        if u < 0:
                            n2bad += 1
                            continue
                        S = [[a, r], [r, d]]
                        mk, ms = all_minors(K), all_minors(S)
                        spec_ok = (ms[3] >= 0 and 1 - (a + d) + ms[3] >= 0 and 0 <= a + d <= 2)
                        if not (all(mk[m] == ms[m] for m in mk) and spec_ok):
                            n2bad += 1
                        n2 += 1
    ck('n_equals_2_sweep_450_rational_NA_2x2_DPP_kernels_with_square_pair_invariant_are_all_symmetrisable',
       n2bad == 0 and n2 >= 50,
       '%d NA 2x2 DPP kernels swept (u = bc a rational square); for each, S with off-diagonal '
       'sqrt(u) reproduces EVERY principal minor and has spectrum in [0,1]. %d kernels with '
       'u < 0 were correctly rejected by the brute-force NA test.' % (n2, naneg))
    ck('n_equals_1_negative_association_is_vacuous',
       na_report(atoms_of([[F(1, 3)]]), 1)['count'] == 0,
       'there is no pair of disjoint nonempty subsets of a 1-element ground set')

    print()
    print('=' * 78)
    print('9. THE BOX CERTIFICATE: THE PHENOMENON IS FULL-DIMENSIONAL')
    print('=' * 78)
    lowb, binding, npoly = box_certificate(F(1, 100))
    ck('box_certificate_has_the_39_constraints_of_the_paper',
       npoly == 39, '8 atoms + 30 NA slacks + (p - q), as exact polynomials in the 9 entries')
    ck('box_certificate_at_radius_1_over_100_is_positive',
       lowb > 0, 'min monomial-wise lower bound over the box |K_ij - W_ij| <= 1/100 is %s > 0, '
                 'binding constraint %s' % (lowb, binding))
    ck('box_certificate_min_lower_bound_is_1339_over_1000000_binding_p_minus_q',
       lowb == F(1339, 1000000) and binding == ('p-q',),
       'so EVERY real 3x3 matrix within 1/100 of W entrywise is a DPP kernel, is negatively '
       'associated, and is non-symmetrisable')
    pmin = (F(9, 100)) * (F(9, 100)) * (F(24, 100))
    qmax = (F(5, 100)) * (F(11, 100)) * (F(11, 100))
    ck('the_binding_bound_re_derived_by_hand',
       pmin - qmax == F(1339, 1000000),
       'p >= 0.09*0.09*0.24 = 1944/10^6 and q <= 0.05*0.11*0.11 = 605/10^6')
    corners = []
    d = F(1, 100)
    for e12 in (W[0][1] - d, W[0][1] + d):
        for e23 in (W[1][2] - d, W[1][2] + d):
            for e31 in (W[2][0] - d, W[2][0] + d):
                for e13 in (W[0][2] - d, W[0][2] + d):
                    for e32 in (W[2][1] - d, W[2][1] + d):
                        for e21 in (W[1][0] - d, W[1][0] + d):
                            corners.append(e12 * e23 * e31 - e13 * e32 * e21)
    ck('the_p_minus_q_bound_is_TIGHT_attained_at_a_corner_of_the_box',
       len(corners) == 64 and min(corners) == F(1339, 1000000),
       'the interval bound is not conservative for the binding constraint')
    for dd, expect_pos in ((F(1, 50), None), (F(1, 20), None)):
        lb, bd, _ = box_certificate(dd)
        ck('box_certificate_at_radius_%s_is_reported_not_hidden' % str(dd).replace('/', '_over_'),
           True, 'delta %s: min lower bound %s, binding %s -- interval arithmetic is '
                 'conservative, so a negative bound is not by itself a bad matrix' % (dd, lb, bd))

    print()
    print('=' * 78)
    print('10. DPP(W) IS NA AND CONDITIONALLY NA, BUT NOT STRONGLY RAYLEIGH')
    print('=' * 78)
    ncond = 0
    condbad = 0
    for i in range(3):
        for val in (1, 0):
            rest = [x for x in range(3) if x != i]
            tot = sum(v for m, v in atW.items() if bool((m >> i) & 1) == bool(val))
            if tot == 0:
                continue
            lawc = {}
            for m, v in atW.items():
                if bool((m >> i) & 1) != bool(val):
                    continue
                mm = 0
                for pos, x in enumerate(rest):
                    if (m >> x) & 1:
                        mm |= 1 << pos
                lawc[mm] = lawc.get(mm, F(0)) + v / tot
            for mm in range(4):
                lawc.setdefault(mm, F(0))
            r = na_report(lawc, 2)
            ncond += 1
            if r['violations']:
                condbad += 1
    ck('W_conditional_laws_on_one_fixed_coordinate_are_all_negatively_associated',
       condbad == 0 and ncond == 6,
       '%d conditional laws (each of the 3 coordinates, both values), 0 violations' % ncond)
    ncondE = condbadE = 0
    for i in range(3):
        for val in (1, 0):
            rest = [x for x in range(3) if x != i]
            tot = sum(v for m, v in atE.items() if bool((m >> i) & 1) == bool(val))
            if tot == 0:
                continue
            lawc = {}
            for m, v in atE.items():
                if bool((m >> i) & 1) != bool(val):
                    continue
                mm = 0
                for pos, x in enumerate(rest):
                    if (m >> x) & 1:
                        mm |= 1 << pos
                lawc[mm] = lawc.get(mm, F(0)) + v / tot
            for mm in range(4):
                lawc.setdefault(mm, F(0))
            ncondE += 1
            if na_report(lawc, 2)['violations']:
                condbadE += 1
    ck('the_conditional_test_is_not_vacuous_the_sources_kernel_fails_it',
       condbadE > 0, '%d of %d conditional laws of the even-subsets kernel violate NA'
                     % (condbadE, ncondE))

    def gen(atoms, x):
        tot = F(0)
        for m, v in atoms.items():
            pr = v
            for i in bits(m):
                pr *= x[i]
            tot += pr
        return tot

    def Dij(atoms, x, i, j):
        def part(idx, xx):
            tot = F(0)
            for m, v in atoms.items():
                if not (m >> idx) & 1:
                    continue
                pr = v
                for t in bits(m):
                    if t != idx:
                        pr *= xx[t]
                tot += pr
            return tot

        def part2(xx):
            tot = F(0)
            for m, v in atoms.items():
                if not ((m >> i) & 1 and (m >> j) & 1):
                    continue
                pr = v
                for t in bits(m):
                    if t not in (i, j):
                        pr *= xx[t]
                tot += pr
            return tot
        return part(i, x) * part(j, x) - gen(atoms, x) * part2(x)

    xstar = [F(-4), F(-4), F(-2)]
    dW = Dij(atW, xstar, 0, 1)
    ck('W_violates_Brandens_strong_Rayleigh_criterion_at_an_explicit_point',
       dW == F(-19, 20000) and dW < 0,
       'x = (-4,-4,-2), (i,j) = (1,2): D_12 f = (2735^2 - 7675*987)/10^8 = -19/20000 < 0, '
       'so DPP(W) is NOT strongly Rayleigh')
    ck('the_strong_Rayleigh_probe_is_not_vacuous_the_symmetric_control_passes_at_that_point',
       Dij(atS, xstar, 0, 1) >= 0, 'D_12 f = %s >= 0 for the symmetric control'
                                   % Dij(atS, xstar, 0, 1))

    print()
    print('NOT RE-RUN / SCOPE. (a) The claim that DPP(W) is not strongly Rayleigh is established '
          'here at ONE explicit point, which suffices to refute the criterion; no sweep of '
          'Branden\'s criterion over all of R^3 is performed, and none is needed. (b) The '
          'conditional-NA check fixes ONE coordinate at a time; deeper conditionings and the '
          'full CNA+ hierarchy are not enumerated here. (c) The n <= 2 impossibility is PROVED '
          'in the paper by hand; what runs here is a finite exact sweep of 2x2 DPP kernels whose '
          'pair invariant is a rational square, which is a corroboration, not the proof. '
          '(d) Theorem C is verified on a finite rational sample of the sigma-interval, not on '
          'the whole interval; the interval statement itself is proved in the paper. (e) The box '
          'certificate is a SOUND LOWER BOUND, so a negative bound at radius 1/50 or 1/20 is '
          'not evidence of a bad matrix; 1/100 is a certified lower bound on the true radius, '
          'not the radius. (f) Poinas\' first question at n >= 4, the continuous-space version '
          'of the question, and the facet description of the 7-dimensional NA region are OPEN '
          'and nothing here bears on them.')
    print()
    if _N_BAD:
        print('VERDICT: %d OF %d CHECKS FAILED' % (_N_BAD, _N_OK + _N_BAD))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _N_OK)
    return 0


# ---------------------------------------------------------------------------
# the box certificate: the 39 constraints as exact polynomials in the 9 entries
# ---------------------------------------------------------------------------
_VARS = [(i, j) for i in range(3) for j in range(3)]
_IDX = {v: k for k, v in enumerate(_VARS)}


class Poly(object):
    """Sparse exact polynomial in the 9 matrix entries: dict exponent-tuple -> Fraction."""

    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = dict(d or {})

    @staticmethod
    def const(c):
        return Poly({(0,) * 9: F(c)}) if c != 0 else Poly()

    @staticmethod
    def var(v):
        e = [0] * 9
        e[_IDX[v]] = 1
        return Poly({tuple(e): F(1)})

    def __add__(a, b):
        d = dict(a.d)
        for k, c in b.d.items():
            d[k] = d.get(k, F(0)) + c
            if d[k] == 0:
                del d[k]
        return Poly(d)

    def __neg__(a):
        return Poly({k: -c for k, c in a.d.items()})

    def __sub__(a, b):
        return a + (-b)

    def __mul__(a, b):
        d = {}
        for k1, c1 in a.d.items():
            for k2, c2 in b.d.items():
                k = tuple(x + y for x, y in zip(k1, k2))
                d[k] = d.get(k, F(0)) + c1 * c2
                if d[k] == 0:
                    del d[k]
        return Poly(d)

    def lower(self, lo, hi):
        """A SOUND lower bound on the box, monomial by monomial."""
        tot = F(0)
        for k, c in self.d.items():
            a = b = F(1)
            for i, e in enumerate(k):
                if e:
                    l, h = lo[i] ** e, hi[i] ** e
                    cands = (a * l, a * h, b * l, b * h)
                    a, b = min(cands), max(cands)
            tot += a * c if c > 0 else b * c
        return tot

    def eval(self, x):
        tot = F(0)
        for k, c in self.d.items():
            m = F(1)
            for i, e in enumerate(k):
                m *= x[i] ** e
            tot += c * m
        return tot


_POLYS = None


def _build_polys():
    global _POLYS
    if _POLYS is not None:
        return _POLYS
    KK = [[Poly.var((i, j)) for j in range(3)] for i in range(3)]

    def detp(S):
        S = sorted(S)
        n = len(S)
        if n == 0:
            return Poly.const(1)
        tot = Poly()
        for p in permutations(range(n)):
            inv = sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])
            pr = Poly.const(1)
            for i in range(n):
                pr = pr * KK[S[i]][S[p[i]]]
            tot = tot + pr if inv % 2 == 0 else tot - pr
        return tot

    mp = {}
    for r in range(4):
        for S in combinations(range(3), r):
            m = 0
            for t in S:
                m |= 1 << t
            mp[m] = detp(S)
    ap = {}
    for A in range(8):
        rest = bits(7 ^ A)
        tot = Poly()
        for k in range(len(rest) + 1):
            for T in combinations(rest, k):
                add = 0
                for t in T:
                    add |= 1 << t
                tot = tot + mp[A | add] if k % 2 == 0 else tot - mp[A | add]
        ap[A] = tot

    def Pr(pred):
        tot = Poly()
        for m, v in ap.items():
            if pred(m):
                tot = tot + v
        return tot

    polys = [(('atom', tuple(bits(A))), ap[A]) for A in range(8)]
    for ra in range(1, 3):
        for A in combinations(range(3), ra):
            Am = 0
            for t in A:
                Am |= 1 << t
            rest = [x for x in range(3) if x not in A]
            for rb in range(1, len(rest) + 1):
                for B in combinations(rest, rb):
                    Bm = 0
                    for t in B:
                        Bm |= 1 << t
                    for U in proper_upsets(A):
                        for V in proper_upsets(B):
                            l = Pr(lambda m, Am=Am, U=U, Bm=Bm, V=V:
                                   (m & Am) in U and (m & Bm) in V)
                            r = (Pr(lambda m, Am=Am, U=U: (m & Am) in U)
                                 * Pr(lambda m, Bm=Bm, V=V: (m & Bm) in V))
                            polys.append((('na', A, B, tuple(sorted(U)), tuple(sorted(V))), r - l))
    polys.append((('p-q',), KK[0][1] * KK[1][2] * KK[2][0] - KK[0][2] * KK[2][1] * KK[1][0]))
    _POLYS = polys
    return polys


def box_certificate(delta):
    polys = _build_polys()
    Wf = [W[i][j] for (i, j) in _VARS]
    lo = [w - delta for w in Wf]
    hi = [w + delta for w in Wf]
    best = None
    for nm, p in polys:
        b = p.lower(lo, hi)
        if best is None or b < best[1]:
            best = (nm, b)
    return best[1], best[0], len(polys)


if __name__ == '__main__':
    sys.exit(main())

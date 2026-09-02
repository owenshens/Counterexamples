#!/usr/bin/env python3
"""Verification program for

    "Pure powers 3^7 and 3^8 in the integer group determinant image of D_54"

Python 3.9+, STANDARD LIBRARY ONLY, no external data file.  Exact integer
arithmetic throughout; no floating point is used for any decision.

Every input is transcribed from the note: the four coefficient vectors of the
main theorem (W1, W2) and the published control witness of [BP] read at n = 9.
Nothing is read from disk.

The necessity statement of [BP] and its achievability results for b = 0 and
b >= 9 are quoted, not verified here.

One `PASS <name>` line per check; the run closes with
    VERDICT: ALL <n> CHECKS PASS
and exits 0 if and only if every check passed.
"""

import itertools
import sys

# ---------------------------------------------------------------------------
# check bookkeeping
# ---------------------------------------------------------------------------
_PASSED = []
_FAILED = []


def ck(name, cond, detail=''):
    if cond:
        _PASSED.append(name)
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILED.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE %s' % text)


# ---------------------------------------------------------------------------
# exact integer linear algebra: fraction-free (Bareiss) determinant
# ---------------------------------------------------------------------------
def det_int(M):
    """Exact determinant of an integer matrix, fraction-free Bareiss.

    Every division below is exact by the Bareiss identity; it is asserted, so a
    silent truncation cannot pass for an answer.
    """
    n = len(M)
    A = [list(row) for row in M]
    for row in A:
        assert len(row) == n
        for x in row:
            assert isinstance(x, int)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            piv = None
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    piv = i
                    break
            if piv is None:
                return 0
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = A[i][j] * A[k][k] - A[i][k] * A[k][j]
                q, r = divmod(num, prev)
                assert r == 0, 'Bareiss division was not exact'
                A[i][j] = q
            A[i][k] = 0
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


# ---------------------------------------------------------------------------
# the group D_2n = < r, s | r^n = s^2 = 1, s r s = r^{-1} >
# elements are pairs (i, e) meaning r^i s^e
# ---------------------------------------------------------------------------
def elements(n):
    return [(i, 0) for i in range(n)] + [(i, 1) for i in range(n)]


def gmul(n, a, b):
    i, e = a
    j, f = b
    if e == 0:
        return ((i + j) % n, f)
    return ((i - j) % n, 0 if f == 1 else 1)


def ginv(n, a):
    i, e = a
    return ((-i) % n, 0) if e == 0 else (i, 1)


def group_det(n, A, B):
    """The DEFINITION: det(a_{u v^{-1}}) over all 2n x 2n pairs of elements,
    with a_{r^i} = A[i] and a_{r^i s} = B[i]."""
    els = elements(n)

    def coef(a):
        i, e = a
        return A[i] if e == 0 else B[i]

    M = [[coef(gmul(n, u, ginv(n, v))) for v in els] for u in els]
    return det_int(M)


# ---------------------------------------------------------------------------
# polynomials in Z[x]/(x^n - 1) and honest polynomial resultants
# ---------------------------------------------------------------------------
def vec(n, d):
    v = [0] * n
    for i, c in d.items():
        v[i % n] += c
    return v


def star(n, f):
    return [f[(-i) % n] for i in range(n)]


def cmul(n, a, b):
    r = [0] * n
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[(i + j) % n] += x * y
    return r


def hpoly(n, f, g):
    ff = cmul(n, f, star(n, f))
    gg = cmul(n, g, star(n, g))
    return [ff[i] - gg[i] for i in range(n)]


def circulant_det(n, h):
    return det_int([[h[(i - j) % n] for j in range(n)] for i in range(n)])


def poly_trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def poly_div_exact(a, b):
    """a / b in Z[x], asserted exact."""
    a = list(a)
    b = poly_trim(b)
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(q) - 1, -1, -1):
        num = a[i + len(b) - 1]
        c, r = divmod(num, b[-1])
        assert r == 0, 'polynomial division was not exact'
        q[i] = c
        for j, y in enumerate(b):
            a[i + j] -= c * y
    assert all(x == 0 for x in a), 'polynomial division left a remainder'
    return q


def cyclotomic(d):
    """Phi_d, by dividing x^d - 1 by the Phi_e for the proper divisors e."""
    num = [-1] + [0] * (d - 1) + [1]
    for e in range(1, d):
        if d % e == 0:
            num = poly_div_exact(num, cyclotomic(e))
    return num


def resultant(a, b):
    """Res(a, b) as the Sylvester determinant; coefficients low -> high."""
    a = poly_trim(a)
    b = poly_trim(b)
    m = len(a) - 1
    k = len(b) - 1
    N = m + k
    M = [[0] * N for _ in range(N)]
    for i in range(k):
        for j, c in enumerate(a):
            M[i][i + (m - j)] = c
    for i in range(m):
        for j, c in enumerate(b):
            M[k + i][i + (k - j)] = c
    return det_int(M)


def xn_minus_1(n):
    return [-1] + [0] * (n - 1) + [1]


def vp(x, p):
    """v_p(x); None for x = 0, which is a different answer from 0."""
    if x == 0:
        return None
    x = abs(x)
    k = 0
    while x % p == 0:
        x //= p
        k += 1
    return k


def is_square(x):
    if x < 0:
        return False
    r = int(x ** 0.5)
    while r * r > x:
        r -= 1
    while (r + 1) * (r + 1) <= x:
        r += 1
    return r * r == x, r


# ---------------------------------------------------------------------------
# THE OBJECTS, TRANSCRIBED FROM THE PAPER
# ---------------------------------------------------------------------------
# Theorem 1: coefficient vectors of f and g, as printed.
W1_F = {1: 1, 4: 1}
W1_G = {0: -1, 1: 1, 4: 1}
W1_H_PRINTED = {0: -1, 1: 1, 4: 1, 23: 1, 26: 1}
W1_H1_PRINTED = 3
W1_D_PRINTED = 2187

W2_F = {0: 1, 1: 1, 2: 1, 4: 1, 5: 1}
W2_G = {1: 1, 2: 1, 4: 1, 5: 1}
W2_H_PRINTED = {0: 1, 1: 1, 2: 1, 4: 1, 5: 1, 22: 1, 23: 1, 25: 1, 26: 1}
W2_H1_PRINTED = 9
W2_D_PRINTED = 6561

# Sections 5 and 6: the published p^5 witness of [BP], f = 1+x, g = 1+x-x^2.
CTRL_F = {0: 1, 1: 1}
CTRL_G = {0: 1, 1: 1, 2: -1}
CTRL_D_AT_N9 = 243            # the published p^5 measure of D_18
CTRL_D_AT_N27 = 6143283       # = 3^7 * 53^2, quoted in Sections 5 and 6

RES_PRINTED = 9               # Res(Phi_{3^m}, h) for m = 1, 2, 3, both witnesses


def main():
    print('verification of "Pure powers 3^7 and 3^8 in the integer group '
          'determinant image of D_54"')
    print('python %s -- exact integer arithmetic only, standard library only'
          % sys.version.split()[0])
    print('')

    # ==================================================================
    print('=== Step 1: the group D_54, built from its presentation')
    n = 27
    els = elements(n)
    ck('group_has_order_54', len(els) == 54 and len(set(els)) == 54,
       'order = %d' % len(els))
    ck('closure_under_multiplication',
       all(gmul(n, u, v) in set(els) for u in els for v in els))
    e = (0, 0)
    ck('identity_element', all(gmul(n, e, u) == u and gmul(n, u, e) == u for u in els))
    ck('inverses_are_two_sided',
       all(gmul(n, u, ginv(n, u)) == e and gmul(n, ginv(n, u), u) == e for u in els))
    bad = [(u, v, w) for u in els for v in els for w in els
           if gmul(n, gmul(n, u, v), w) != gmul(n, u, gmul(n, v, w))]
    ck('associativity_on_all_157464_triples', not bad,
       '54^3 = %d triples, %d failures' % (54 ** 3, len(bad)))
    r, s = (1, 0), (0, 1)
    rinv = ginv(n, r)
    ck('dihedral_relation_s_r_s_equals_r_inverse',
       gmul(n, gmul(n, s, r), s) == rinv)
    ordr = 1
    x = r
    while x != e:
        x = gmul(n, x, r)
        ordr += 1
    ck('r_has_order_27', ordr == 27, 'order(r) = %d' % ordr)
    ck('s_has_order_2', gmul(n, s, s) == e and s != e)

    # ==================================================================
    print('')
    print('=== Step 2: a wiring control on the group determinant')
    # [BP]-style spot value: a_0 = 0 and every other coefficient 1.
    A = vec(n, {i: 1 for i in range(1, 27)})
    B = [1] * 27
    d0 = group_det(n, A, B)
    ck('control_group_determinant_of_0_then_all_ones_is_minus_53', d0 == -53,
       'D = %d' % d0)
    note('a NEGATIVE value: the 54x54 determinant is not being computed up to sign '
         'or in absolute value.')

    # multiplicativity of D_G on Z[G], on small deterministic elements
    def gring_mul(nn, a1, b1, a2, b2):
        """(f1 + g1 s)(f2 + g2 s) = (f1 f2 + g1 g2*) + (f1 g2 + g1 f2*) s."""
        f = [u + v for u, v in zip(cmul(nn, a1, a2), cmul(nn, b1, star(nn, b2)))]
        g = [u + v for u, v in zip(cmul(nn, a1, b2), cmul(nn, b1, star(nn, a2)))]
        return f, g

    # ⛔ every pair below has a NONZERO product: 0 == 0 would satisfy the identity
    # for a routine that returned 0 unconditionally.
    pairs = [(({1: 1, 4: 1}, {0: -1, 1: 1, 4: 1}), ({0: 1}, {3: 1, 5: -1})),
             (({0: 2}, {1: -1}), ({0: 1, 1: -1}, {0: 1})),
             (({0: 5}, {1: 1, 2: 1}), ({0: 2, 1: 1}, {}))]
    okm, det_m = True, []
    for (p1, q1), (p2, q2) in pairs:
        a1, b1 = vec(n, p1), vec(n, q1)
        a2, b2 = vec(n, p2), vec(n, q2)
        f3, g3 = gring_mul(n, a1, b1, a2, b2)
        lhs = group_det(n, a1, b1) * group_det(n, a2, b2)
        rhs = group_det(n, f3, g3)
        det_m.append((lhs, rhs))
        okm = okm and (lhs == rhs) and lhs != 0
    ck('group_determinant_is_multiplicative_on_three_deterministic_nonzero_pairs',
       okm, '; '.join('%d==%d' % t for t in det_m))

    # ==================================================================
    print('')
    print('=== Step 3: the two witnesses of the main theorem, by four routes')
    divisors_27 = [1, 3, 9, 27]
    results = {}
    for tag, fd, gd, hpr, h1pr, dpr, v3exp in (
            ('W1', W1_F, W1_G, W1_H_PRINTED, W1_H1_PRINTED, W1_D_PRINTED, 7),
            ('W2', W2_F, W2_G, W2_H_PRINTED, W2_H1_PRINTED, W2_D_PRINTED, 8)):
        f = vec(n, fd)
        g = vec(n, gd)
        h = hpoly(n, f, g)

        ck('%s_g_equals_f_minus_1_as_printed' % tag,
           g == [f[0] - 1] + f[1:], 'g - (f-1) = 0')
        ck('%s_h_recomputed_from_f_and_g_matches_the_printed_h' % tag,
           {i: c for i, c in enumerate(h) if c} == hpr,
           'h = %s' % {i: c for i, c in enumerate(h) if c})
        # the shortcut of equation (2): g = f - 1  =>  h = f + f* - 1
        hshort = [f[i] + star(n, f)[i] for i in range(n)]
        hshort[0] -= 1
        ck('%s_h_equals_f_plus_fstar_minus_one_equation_2' % tag, hshort == h)
        ck('%s_h_is_symmetric' % tag, all(h[i] == h[(-i) % n] for i in range(n)))
        ck('%s_h_at_1_matches_the_printed_value' % tag, sum(h) == h1pr,
           'h(1) = %d' % sum(h))
        ck('%s_h_at_1_is_odd_so_h_is_a_genuine_reduced_polynomial' % tag,
           sum(h) % 2 == 1)
        # the realizability converse stated in the paper: a_0=(h_0+1)/2, a_i=h_i
        assert (h[0] + 1) % 2 == 0
        a_re = [0] * n
        a_re[0] = (h[0] + 1) // 2
        for i in range(1, 14):
            a_re[i] = h[i]
        g_re = [a_re[0] - 1] + a_re[1:]
        ck('%s_the_stated_realization_of_h_reproduces_h' % tag,
           hpoly(n, a_re, g_re) == h,
           'a_0 = %d, a_1..a_13 = %s' % (a_re[0], a_re[1:14]))

        d_group = group_det(n, f, g)
        d_circ = circulant_det(n, h)
        d_res = resultant(xn_minus_1(n), h)
        resd = {d: resultant(cyclotomic(d), h) for d in divisors_27}
        prod = 1
        for d in divisors_27:
            prod *= resd[d]

        ck('%s_route1_54x54_group_table_determinant' % tag, d_group == dpr,
           'D = %d' % d_group)
        ck('%s_route2_27x27_circulant_of_h_agrees' % tag, d_circ == d_group,
           'D = %d' % d_circ)
        ck('%s_route3_resultant_of_x27_minus_1_and_h_agrees' % tag,
           d_res == d_group, 'D = %d' % d_res)
        ck('%s_route4_product_of_cyclotomic_resultants_agrees' % tag,
           prod == d_group, 'h(1)*Res(Phi_3)*Res(Phi_9)*Res(Phi_27) = %d' % prod)
        ck('%s_res_at_Phi_1_is_h_at_1' % tag, resd[1] == sum(h))
        ck('%s_each_Res_Phi_3m_equals_9_as_printed' % tag,
           all(resd[3 ** m] == RES_PRINTED for m in (1, 2, 3)),
           'Res = %s' % {3 ** m: resd[3 ** m] for m in (1, 2, 3)})
        sq = [is_square(resd[3 ** m]) for m in (1, 2, 3)]
        ck('%s_each_Res_Phi_3m_is_a_perfect_square_P_m_squared_lemma_2' % tag,
           all(t[0] for t in sq), 'P_m = +-%s' % [t[1] for t in sq])
        ck('%s_v3_of_each_P_m_is_1' % tag,
           all(vp(t[1], 3) == 1 for t in sq),
           'v3(P_m) = %s' % [vp(t[1], 3) for t in sq])
        ck('%s_valuation_formula_of_lemma_2_holds' % tag,
           vp(d_group, 3) == vp(sum(h), 3) + 2 * sum(vp(t[1], 3) for t in sq),
           'v3(D) = v3(h(1)) + 2*sum v3(P_m) = %d + %d'
           % (vp(sum(h), 3), 2 * sum(vp(t[1], 3) for t in sq)))
        ck('%s_v3_of_D_is_%d' % (tag, v3exp), vp(d_group, 3) == v3exp,
           'v3 = %s' % vp(d_group, 3))
        ck('%s_v2_of_D_is_0' % tag, vp(d_group, 2) == 0)
        cof = abs(d_group) // (3 ** vp(d_group, 3))
        ck('%s_cofactor_after_removing_all_2s_and_3s_is_1_a_PURE_power' % tag,
           cof == 1, 'cofactor = %d' % cof)
        ck('%s_D_equals_3_to_the_%d_exactly' % (tag, v3exp),
           d_group == 3 ** v3exp, '%d == 3^%d' % (d_group, v3exp))

        # the hand identity for P_1 stated in the paper's Remark
        H1 = sum(h[i] for i in range(1, 14) if i % 3 in (1, 2))
        P1 = sum(h) - 3 * H1
        ck('%s_hand_identity_P1_equals_h1_minus_3H1' % tag,
           P1 * P1 == resd[3], 'H_1 = %d, P_1 = %d, P_1^2 = %d' % (H1, P1, P1 * P1))

        # the other polarity
        ck('%s_swapping_f_and_g_negates_the_determinant' % tag,
           group_det(n, g, f) == -d_group, 'D = %d' % group_det(n, g, f))
        results[tag] = d_group

    ck('the_two_witnesses_give_3_to_the_7_and_3_to_the_8',
       (results['W1'], results['W2']) == (3 ** 7, 3 ** 8),
       '%d and %d' % (results['W1'], results['W2']))

    # ==================================================================
    print('')
    print('=== Step 4: forced positive control on a PUBLISHED value')
    # [BP]'s p^5 witness for D_18: f = 1+x, g = 1+x-x^2, n = 9  ->  243
    f9 = vec(9, CTRL_F)
    g9 = vec(9, CTRL_G)
    h9 = hpoly(9, f9, g9)
    d9 = [group_det(9, f9, g9), circulant_det(9, h9),
          resultant(xn_minus_1(9), h9)]
    p9 = 1
    for d in (1, 3, 9):
        p9 *= resultant(cyclotomic(d), h9)
    d9.append(p9)
    ck('published_p5_witness_of_D18_returns_243_by_all_four_routes',
       all(v == CTRL_D_AT_N9 for v in d9), 'D = %s' % d9)
    ck('published_p5_witness_of_D18_has_v3_5_and_unit_cofactor',
       vp(CTRL_D_AT_N9, 3) == 5 and CTRL_D_AT_N9 == 3 ** 5)

    # the same witness read at n = 27: existence of b = 7, but NOT a pure power
    f27 = vec(27, CTRL_F)
    g27 = vec(27, CTRL_G)
    d27 = group_det(27, f27, g27)
    ck('the_same_published_witness_at_n_27_gives_3_to_the_7_times_53_squared',
       d27 == CTRL_D_AT_N27 and d27 == 3 ** 7 * 53 ** 2,
       'D = %d = 3^7 * 53^2' % d27)
    ck('and_its_cofactor_is_NOT_1_so_it_is_not_a_pure_power',
       abs(d27) // 3 ** vp(d27, 3) == 53 ** 2,
       'cofactor = %d' % (abs(d27) // 3 ** vp(d27, 3)))

    # ==================================================================
    print('')
    print('=== Step 5: negative control -- the forbidden valuation band 1..6')
    # A deterministic family: every symmetric h with h_0 odd and support in
    # {0,+-1,+-2,+-3,+-4}, coefficients in {-1,0,1}.  Lemma 4.4 of [BP] at
    # n = 27 forbids v_3(D) in 1..6 and v_2(D) = 1.
    hist, hist2 = {}, {}
    trials = 0
    for sup in itertools.product([-1, 0, 1], repeat=5):
        hv = [0] * 27
        hv[0] = 2 * sup[0] - 1        # odd by construction
        for i in range(1, 5):
            hv[i] = sup[i]
            hv[27 - i] = sup[i]
        D = resultant(xn_minus_1(27), hv)
        trials += 1
        hist[vp(D, 3)] = hist.get(vp(D, 3), 0) + 1
        hist2[vp(D, 2)] = hist2.get(vp(D, 2), 0) + 1
    band = [v for v in hist if v is not None and 1 <= v <= 6]
    ck('no_v3_in_the_forbidden_band_1_to_6_over_%d_deterministic_cases' % trials,
       not band, 'v3 histogram = %s'
       % dict(sorted(hist.items(), key=lambda t: (t[0] is None, t[0]))))
    ck('the_band_test_is_not_silent_both_7_and_8_actually_occur',
       hist.get(7, 0) > 0 and hist.get(8, 0) > 0,
       'v3 = 7 in %d cases, v3 = 8 in %d cases' % (hist.get(7, 0), hist.get(8, 0)))
    ck('no_v2_equal_to_1_over_the_same_family', 1 not in hist2,
       'v2 histogram = %s'
       % dict(sorted(hist2.items(), key=lambda t: (t[0] is None, t[0]))))

    # ==================================================================
    print('')
    print('=== Step 6: the exponent semigroup used in the completeness theorem')
    # [BP] supplies b = 0 and every b >= 9; the note adds b = 7 and b = 8.
    BOUND = 60
    gens = set([0, 7, 8]) | set(range(9, BOUND + 1))
    S = set(gens)
    for _ in range(8):
        S |= {a + b for a in S for b in S if a + b <= BOUND}
    ck('attainable_b_set_up_to_%d_is_exactly_0_together_with_7_and_above' % BOUND,
       sorted(v for v in S if v <= BOUND) == [0] + list(range(7, BOUND + 1)),
       'min positive attainable b = %d' % min(v for v in S if v > 0))
    ck('the_two_new_pure_powers_are_what_closes_the_gap_7_and_8',
       7 in gens and 8 in gens
       and sorted(v for v in ({0} | set(range(9, BOUND + 1))) if v <= 8) == [0],
       'without b = 7, 8 the source alone attains only b = 0 or b >= 9')
    ck('the_band_1_to_6_stays_unattainable_matching_Theorem_5_4_necessity',
       not any(1 <= v <= 6 for v in S))

    # ==================================================================
    print('')
    note('SCOPE: every claim re-derived above is a claim the note makes about the '
         'objects printed in it. NOT RE-RUN, and not covered by any check here: '
         '(1) the achievability results of [BP] that the completeness theorem '
         'consumes -- every 2^a m with (m,6)=1 and a = 0 or a >= 2 (the b = 0 '
         'members), and every pure power 3^b with b >= 9 -- which are quoted from '
         'the source and not verified; (2) the necessity statement of [BP] '
         'Theorem 5.4 / Lemma 4.4 in general, including the p = 2 clauses, which '
         'is likewise quoted and not verified: Step 5 is a finite deterministic '
         'CONTROL over 243 symmetric h with support in {0,+-1,+-2,+-3,+-4} and '
         'coefficients in {-1,0,1}, not a proof; '
         '(3) the measure set of any group other than D_54 -- nothing here '
         'determines S(D_162), S(D_250) or S of any p > 3; '
         '(4) minimality, uniqueness or exhaustiveness of W1 and W2, on which '
         'the note makes no claim and this program runs no census.')
    print('')
    n_ok, n_bad = len(_PASSED), len(_FAILED)
    if n_bad:
        print('VERDICT: %d of %d CHECKS FAILED' % (n_bad, n_ok + n_bad))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % n_ok)
    return 0


if __name__ == '__main__':
    sys.exit(main())

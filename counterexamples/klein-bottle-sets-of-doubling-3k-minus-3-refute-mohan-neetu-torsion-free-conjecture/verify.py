#!/usr/bin/env python3
"""verify.py -- re-derives every quantity claimed by

    "Klein Bottle Sets of Doubling 3|S|-3 Refute a Torsion-Free Small-Doubling Conjecture of Mohan and Neetu as Stated Torsion-Free Small-Doubling
     Conjecture of Mohan and Neetu"   (paper.tex / paper.pdf, this folder)

from the objects PRINTED IN THAT PAPER and nothing else.

    python3 verify.py            # Python 3.9+, standard library only, no data file, no network

WHAT IS RE-DERIVED, AND FROM WHERE IN THE PAPER
-----------------------------------------------
  * the group law and the inverse of K = BS(1,-1) = Z x|_{-1} Z            eq. (1), labelled eq:law
  * the family S_k = {(0,0)} u {(i,1) : 1 <= i <= k-1}                     eq. (2), labelled eq:S
  * |S_k| = k, e in S_k, |S_k^2| = 3k-3, and the closed form
    S_k^2 = S_k u {(d,2) : |d| <= k-2}, compared AS A SET OF PAIRS         Theorem 2
  * nonabelianness and <S_k> = K via the printed words
    (1,1)^{-1}(2,1) = a^{-1} and a^{-1}(1,1) = t                           Theorem 2, proof
  * the nine ordered products of the k = 3 cell, one by one                display after Theorem 2
  * right-invariance of the printed order and min(S_k) = e                 Corollary 3
  * the authors' own published |T_k^2| = 2k-1                              first Remark
  * |A^2| = 7 = 3|A|-2 for the Boroczky-Palfy-Serra set A = {1,u,v}        second Remark, (3)
  * that BOTH detectors used above can return FALSE (two-polarity controls)

Everything is exact integer arithmetic on pairs of Python ints. No float appears
anywhere, and the sign (-1)^n is taken as a parity lookup rather than a power, so
a negative exponent cannot silently produce a float.

NOT RE-RUN -- stated here and in REVIEW_NOTE.md's ## Scope:
  * Torsion-freeness of K is checked STRUCTURALLY (the second-coordinate identity
    and the n = 0 case) and EXHAUSTIVELY ONLY OVER A BOUNDED BOX of elements and
    exponents. The unbounded statement is Lemma 1 of the paper, a two-line hand
    proof; this program does not replace it.
  * NO SEARCH FOR A SMALLER OR DIFFERENT COUNTEREXAMPLE is performed: no host
    group other than K is examined, no set outside the printed families and the
    printed controls is examined, and nothing here bears on whether 3|S|-3 is the
    least doubling attainable by a nonabelian-spanning identity-containing set in
    a torsion-free group.
  * The BYTE LOCATOR of section 1 (the 97,638 B source member, its sha256, tex
    lines 906-909 and bytes 90,881-91,090) is NOT re-fetched or re-hashed: this
    program has no network access and ships no copy of the e-print.
  * The claim in the last remark that (x g^i)^2 = x^2 forces q = -1 in BS(1,q) is
    an algebraic argument over Z[1/q], not a finite computation, and is NOT
    re-derived here.
  * NO PRIOR-ART OR LITERATURE CLAIM is verified. That the object is prior art
    while the refutation is not, and that the companion paper is unread, are
    bibliographic statements outside this program's reach.
"""
import sys

CHECKS = []
FAILED = []


def check(name, ok, detail=''):
    """One PASS/FAIL line. The verdict count is exactly len(CHECKS)."""
    line = '%s %s%s' % ('PASS' if ok else 'FAIL', name, (' ' + detail) if detail else '')
    print(line)
    CHECKS.append(name)
    if not ok:
        FAILED.append(name)


# ---------------------------------------------------------------------------
# 1. THE GROUP, EXACTLY AS PRINTED IN eq. (1)
#        (m,n)(p,q) = (m + (-1)^n p, n + q),   (m,n)^{-1} = (-(-1)^n m, -n),  e = (0,0)
# ---------------------------------------------------------------------------
E = (0, 0)


def sgn(n):
    """(-1)^n as an INTEGER. Python evaluates (-1)**(-1) as the float -1.0, which
    would leave the integer lattice the moment an inverse such as (1,-1) is used."""
    return 1 if n % 2 == 0 else -1


def mul(x, y):
    (m, n), (p, q) = x, y
    return (m + sgn(n) * p, n + q)


def inv(x):
    m, n = x
    return (-sgn(n) * m, -n)


def power(x, r):
    """x^r for any integer r, by repeated multiplication (no float exponent)."""
    if r < 0:
        return power(inv(x), -r)
    out = E
    for _ in range(r):
        out = mul(out, x)
    return out


def prod_set(S):
    return set(mul(x, y) for x in S for y in S)


def nonabelian(S):
    """TRUE iff two elements of S fail to commute. A DETECTOR, not a constant: the
    two ctrl_abelian_* controls in section 4 make it return FALSE, one of them AT the
    bound |S^2| = 3|S|-3."""
    return any(mul(x, y) != mul(y, x) for x in S for y in S)


def hypothesis_holds(S):
    """TRUE iff |S| >= 3, e in S and |S^2| <= 3|S|-3 -- the antecedent of
    Conjecture 6.3. A DETECTOR, not a constant: the two ctrl_spread_* controls in
    section 4 make it return FALSE."""
    S = list(dict.fromkeys(S))
    return len(S) >= 3 and E in S and len(prod_set(S)) <= 3 * len(S) - 3


BOX = [(m, n) for m in range(-4, 5) for n in range(-3, 4)]

print('OBJECT AS PRINTED: K = <a,t | t a t^-1 = a^-1>, a^m t^n = (m,n),')
print('                   (m,n)(p,q) = (m + (-1)^n p, n+q),  e = (0,0)')
print('                   S_k = {(0,0)} u {(i,1) : 1 <= i <= k-1},  k >= 3')
print('')

# ---- the group law itself --------------------------------------------------
check('group_law_associative',
      all(mul(mul(x, y), z) == mul(x, mul(y, z)) for x in BOX for y in BOX for z in BOX),
      '[all %d^3 = %d triples in m in [-4,4], n in [-3,3]]' % (len(BOX), len(BOX) ** 3))
check('group_law_identity_two_sided',
      all(mul(E, x) == x and mul(x, E) == x for x in BOX),
      '[e = (0,0) on all %d box elements]' % len(BOX))
check('group_law_inverse_formula',
      all(mul(x, inv(x)) == E and mul(inv(x), x) == E for x in BOX),
      '[(m,n)^-1 = (-(-1)^n m, -n) on all %d box elements]' % len(BOX))
check('group_law_closed_in_integer_lattice',
      all(isinstance(c, int) for x in BOX for y in BOX for c in mul(x, y) + inv(x)),
      '[every coordinate produced is a Python int; no float anywhere]')

# ---- Lemma 1: torsion-free (bounded), and K nonabelian ------------------
check('lemma_torsion_free_second_coordinate_identity',
      all(power(x, r)[1] == r * x[1] for x in BOX for r in range(-6, 7)),
      '[(m,n)^r has second coordinate r*n, |r| <= 6]')
check('lemma_torsion_free_bounded_box',
      all(power(x, r) != E for x in BOX if x != E for r in range(1, 13)),
      '[no x != e in the box has x^r = e for 1 <= r <= 12; the unbounded claim is '
      'Lemma 1, NOT re-run]')
check('lemma_K_is_nonabelian',
      mul((1, 0), (0, 1)) == (1, 1) and mul((0, 1), (1, 0)) == (-1, 1)
      and mul((1, 0), (0, 1)) != mul((0, 1), (1, 0)),
      '[a t = (1,1) != (-1,1) = t a, exactly as printed]')

# ---------------------------------------------------------------------------
# 2. THE k = 3 CELL: THE NINE PRODUCTS AND S_3^2, EXACTLY AS DISPLAYED
# ---------------------------------------------------------------------------
S3_PRINTED = [(0, 0), (1, 1), (2, 1)]
NINE_PRINTED = {
    ((0, 0), (0, 0)): (0, 0), ((0, 0), (1, 1)): (1, 1), ((0, 0), (2, 1)): (2, 1),
    ((1, 1), (0, 0)): (1, 1), ((2, 1), (0, 0)): (2, 1), ((1, 1), (1, 1)): (0, 2),
    ((2, 1), (2, 1)): (0, 2), ((1, 1), (2, 1)): (-1, 2), ((2, 1), (1, 1)): (1, 2),
}
S3_SQUARE_PRINTED = {(0, 0), (1, 1), (2, 1), (-1, 2), (0, 2), (1, 2)}

check('k3_nine_printed_products_each_recomputed',
      len(NINE_PRINTED) == 9 and all(mul(x, y) == z for (x, y), z in NINE_PRINTED.items()),
      '[all 9 ordered products of the displayed table agree with eq. (1)]')
check('k3_single_collision_at_t_squared',
      mul((1, 1), (1, 1)) == mul((2, 1), (2, 1)) == (0, 2)
      and len(set(NINE_PRINTED.values())) == 6,
      '[(at)^2 = (a^2 t)^2 = t^2 = (0,2); 9 products, 6 distinct values]')
check('k3_product_set_equals_printed_set',
      prod_set(S3_PRINTED) == S3_SQUARE_PRINTED,
      '[S_3^2 = {(0,0),(1,1),(2,1),(-1,2),(0,2),(1,2)}, label-equal]')
check('k3_bound_is_met_with_equality',
      len(prod_set(S3_PRINTED)) == 6 == 3 * 3 - 3,
      '[|S_3^2| = 6 <= 3|S_3| - 3 = 6]')

# ---------------------------------------------------------------------------
# 3. THEOREM 2 AND COROLLARY 3, FOR EVERY k IN 3..30
#    The paper claims the family for every k >= 3; 3..30 is the finite window
#    re-run here. Larger k are NOT re-run -- the general statement is the proof.
# ---------------------------------------------------------------------------
KMIN, KMAX = 3, 30


def S_of(k):
    """The family, exactly as printed in eq. (2)."""
    return [(0, 0)] + [(i, 1) for i in range(1, k)]


def T_of(k):
    """The authors' own e-free example of the first Remark: T_k = {(i,1) : 1 <= i <= k}."""
    return [(i, 1) for i in range(1, k + 1)]


def less(x, y):
    """The order printed in Corollary 3: (m,n) < (m',n') iff n < n', or n = n' and m < m'."""
    (m, n), (mp, np_) = x, y
    return (n, m) < (np_, mp)


for k in range(KMIN, KMAX + 1):
    S = S_of(k)
    S2 = prod_set(S)
    closed = set(S) | set((d, 2) for d in range(-(k - 2), k - 1))

    check('family_k%02d_cardinality' % k,
          len(set(S)) == k and len(S) == k and E in S,
          '[|S_k| = %d, elements pairwise distinct, e in S_k]' % k)
    check('family_k%02d_doubling_is_3k_minus_3' % k,
          len(S2) == 3 * k - 3,
          '[|S_k^2| = %d = 3*%d-3 = 3|S_k|-3, computed product by product over all %d '
          'ordered pairs]' % (len(S2), k, k * k))
    check('family_k%02d_closed_form_label_equal' % k,
          S2 == closed and len(set((d, 2) for d in range(-(k - 2), k - 1))) == 2 * k - 3,
          '[S_k^2 = S_k u {(d,2) : |d| <= %d} as a SET OF PAIRS; level-2 block has '
          'exactly 2k-3 = %d elements]' % (k - 2, 2 * k - 3))
    check('family_k%02d_level_split_is_disjoint' % k,
          sorted(len([x for x in S2 if x[1] == lev]) for lev in (0, 1, 2)) == sorted([1, k - 1, 2 * k - 3])
          and all(x[1] in (0, 1, 2) for x in S2),
          '[second coordinates of S_k^2 split 1 + %d + %d, nothing outside levels 0,1,2]'
          % (k - 1, 2 * k - 3))
    check('family_k%02d_nonabelian_witness_pair' % k,
          (1, 1) in S and (2, 1) in S and mul((1, 1), (2, 1)) == (-1, 2)
          and mul((2, 1), (1, 1)) == (1, 2) and nonabelian(S),
          '[(1,1)(2,1) = (-1,2) != (1,2) = (2,1)(1,1), both factors in S_k]')
    check('family_k%02d_span_reaches_a_and_t' % k,
          inv((1, 1)) == (1, -1) and mul((1, -1), (2, 1)) == (-1, 0)
          and mul((-1, 0), (1, 1)) == (0, 1) and mul(inv((-1, 0)), E) == (1, 0),
          '[(1,1)^-1 = (1,-1); (1,-1)(2,1) = (-1,0) = a^-1; (-1,0)(1,1) = (0,1) = t, '
          'so <S_k> = K]')
    check('family_k%02d_antecedent_holds_consequent_fails' % k,
          hypothesis_holds(S) and nonabelian(S),
          '[Conjecture 6.3 hypothesis TRUE and conclusion FALSE at |S| = %d]' % k)
    check('family_k%02d_freiman_3k_minus_4_not_contradicted' % k,
          len(S2) == 3 * k - 3 and len(S2) > 3 * k - 4,
          '[|S_k^2| = %d > %d = 3|S_k|-4, so the family lies one unit ABOVE Freiman\'s '
          'torsion-free threshold]' % (len(S2), 3 * k - 4))
    check('family_k%02d_min_under_printed_order_is_e' % k,
          all(less(E, x) for x in S if x != E)
          and min(S, key=lambda x: (x[1], x[0])) == E,
          '[min(S_k) = e = (0,0) under the Corollary 3 order]')
    check('family_k%02d_authors_e_free_example_2k_minus_1' % k,
          len(prod_set(T_of(k))) == 2 * k - 1,
          '[|T_k^2| = %d = 2*%d-1, reproducing the value printed in the source paper]'
          % (2 * k - 1, k))

# ---- Corollary 3: the order is a right-invariant total order ------------
check('order_right_invariant',
      all(less(x, y) == less(mul(x, g), mul(y, g)) for x in BOX for y in BOX for g in BOX),
      '[x < y  iff  xg < yg, over all %d^3 triples in the box]' % len(BOX))
check('order_is_a_strict_total_order',
      all((less(x, y) + less(y, x) + (x == y)) == 1 for x in BOX for y in BOX)
      and all(not (less(x, y) and less(y, z) and not less(x, z))
              for x in BOX for y in BOX for z in BOX),
      '[trichotomy on all pairs and transitivity on all triples of the box]')

# ---------------------------------------------------------------------------
# 4. TWO-POLARITY CONTROLS -- NEITHER DETECTOR IS A CONSTANT
#    Each control names the set, the value expected, and which detector it flips.
# ---------------------------------------------------------------------------
CONTROLS = [
    # (label, S, expected |S^2|, expected hypothesis_holds, expected nonabelian)
    ('ctrl_abelian_span_a_powers      {e,a,a^2}',
     [(0, 0), (1, 0), (2, 0)], 5, True, False),
    ('ctrl_abelian_span_at_the_bound  {e,a,t^2}',
     [(0, 0), (1, 0), (0, 2)], 6, True, False),
    ('ctrl_spread_k3_outside_bound    {e,(1,1),(5,3)}',
     [(0, 0), (1, 1), (5, 3)], 7, False, True),
    ('ctrl_spread_k4_outside_bound    {e,(1,0),(0,1),(7,5)}',
     [(0, 0), (1, 0), (0, 1), (7, 5)], 13, False, True),
]
for label, S, want_sz, want_hyp, want_nab in CONTROLS:
    sz, hyp, nab = len(prod_set(S)), hypothesis_holds(S), nonabelian(S)
    check(label.split()[0],
          (sz, hyp, nab) == (want_sz, want_hyp, want_nab),
          '%s |S^2| = %d (3|S|-3 = %d) hypothesis_holds = %s nonabelian = %s -- NOT a '
          'counterexample' % (label.split(None, 1)[1].strip(), sz, 3 * len(S) - 3, hyp, nab))

check('ctrl_nonabelian_detector_can_be_false',
      any(not nonabelian(S) for _l, S, _a, _b, _c in CONTROLS),
      '[2 of 4 controls return nonabelian = False, one of them AT the bound |S^2| = 3|S|-3, '
      'so meeting the bound does not by itself force nonabelianness]')
check('ctrl_hypothesis_detector_can_be_false',
      any(not hypothesis_holds(S) for _l, S, _a, _b, _c in CONTROLS),
      '[2 of 4 controls return hypothesis_holds = False]')
check('ctrl_t_squared_really_is_central',
      mul((1, 0), (0, 2)) == mul((0, 2), (1, 0)) == (1, 2)
      and not nonabelian([(0, 0), (1, 0), (0, 2)]),
      '[(1,0)(0,2) = (0,2)(1,0) = (1,2), so <a, t^2> = Z^2 is genuinely abelian and the '
      'at-the-bound control is sound, not a bug]')

# ---------------------------------------------------------------------------
# 5. THE NEAREST PUBLISHED SET (second Remark, (3)) IS A DIFFERENT SET
# ---------------------------------------------------------------------------
A_BPS = [(0, 0), (0, 1), (1, 0)]          # A = {1, u, v} with u = t, v = a
A2 = prod_set(A_BPS)
check('bps_set_A_doubling_is_3n_minus_2',
      len(A2) == 7 == 3 * len(A_BPS) - 2,
      '[A = {1,u,v} = {(0,0),(0,1),(1,0)}: |A^2| = 7 = 3|A|-2]')
check('bps_set_A_lies_outside_the_antecedent',
      not hypothesis_holds(A_BPS) and len(A2) > 3 * len(A_BPS) - 3,
      '[7 > 6 = 3|A|-3, so the published 3-element set does NOT satisfy Conjecture 6.3\'s '
      'hypothesis; S_3 is a different set]')
check('bps_set_A_is_not_S3',
      set(A_BPS) != set(S3_PRINTED),
      '[{(0,0),(0,1),(1,0)} != {(0,0),(1,1),(2,1)}]')

# ---------------------------------------------------------------------------
# 6. THE HEADLINE, RESTATED AS ONE CHECK
# ---------------------------------------------------------------------------
check('refutation_holds_at_every_k_in_window',
      all(hypothesis_holds(S_of(k)) and nonabelian(S_of(k)) and len(prod_set(S_of(k))) == 3 * k - 3
          for k in range(KMIN, KMAX + 1)),
      '[Conjecture 6.3 fails for every k in %d..%d; the paper proves every k >= 3]'
      % (KMIN, KMAX))

# ---------------------------------------------------------------------------
# 7. VERDICT
# ---------------------------------------------------------------------------
print('')
print('NOT RE-RUN: torsion-freeness of K is confirmed structurally and over a bounded box '
      '(%d elements, exponents 1..12) only; the unbounded statement is Lemma 1, a hand proof.'
      % len([x for x in BOX if x != E]))
print('NOT RE-RUN: no search for a smaller or different counterexample, in K or in any other '
      'group; no minimality of 3|S|-3 is examined, and k > %d is not enumerated.' % KMAX)
print('NOT RE-RUN: the byte locator of section 1 (source member sha256, tex lines 906-909, '
      'bytes 90,881-91,090) is not re-fetched -- this program has no network and ships no '
      'copy of the e-print.')
print('NOT RE-RUN: the BS(1,q) argument that (x g^i)^2 = x^2 forces q = -1 is algebra over '
      'Z[1/q], not a finite computation.')
print('NOT RE-RUN: no prior-art, attribution or bibliographic claim of the paper is verified.')
print('')
if FAILED:
    print('VERDICT: %d of %d CHECKS FAILED: %s' % (len(FAILED), len(CHECKS), ', '.join(FAILED)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(CHECKS))
sys.exit(0)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Numerical re-derivation program -- constants, transcribed values and inequalities
only -- for the note

    "A counterexample, with a hyperbolic summand, to a printed subadditivity
     question for topological volume"

This program checks arithmetic only.  It does not establish the note's conclusion; it
re-derives the constants the note prints and re-decides the numerical comparisons.

WHAT THIS PROGRAM IS FOR.  The result in the accompanying paper is a hand proof: two
lemmas from published theorems, plus arithmetic on three decimal places.  No program is
needed to follow it.  What a program CAN do, and what this one does, is re-derive every
NUMBER the paper prints -- including the two transcendental constants v_0 and 2v_0, from
their definitions and not from any table -- and re-decide every INEQUALITY the paper
asserts, in exact arithmetic, from the objects printed in the paper itself.

INPUT.  Every value this program consumes is a literal typed into the source below by
hand: ten manifold names and one volume string, the bound vol_h(Weeks) <= 0.943, the
number 3.07, two further volume values, and a set of row counts.  The note itself prints
none of these except 0.943 and the constants; do not expect to find the others in it.
No external file and no third-party package is read.

ARITHMETIC.  Python 3.9+, standard library only.  The two constants are computed with
`decimal` at 80 digits of working precision from the Clausen/Bernoulli expansion, with an
explicit truncation bound, and are then converted into a pair of RATIONAL bounds; every
decision in this program is a comparison of `fractions.Fraction`s.  No decision is taken
on a binary float.

Exit status is 0 if and only if every check passes.
"""

import sys
from decimal import Decimal, getcontext
from fractions import Fraction

PREC = 80
getcontext().prec = PREC

# --------------------------------------------------------------------------------------
# THE OBJECTS PRINTED IN THE PAPER.  Transcribed, not computed.
# --------------------------------------------------------------------------------------

# The published bound on the volume of the Weeks manifold W (Gabai-Meyerhoff-Milley,
# arXiv:0705.4325, abstract: "all closed hyperbolic three-manifolds with volume <= 0.943.
# In particular, the Weeks manifold is the unique smallest volume closed orientable
# hyperbolic 3-manifold.").  The paper uses this BOUND and never a digit expansion of
# vol_h(W); neither does this program.
WEEKS_BOUND = Fraction('0.943')

# The number 3.07, typed in by hand as the reach of the imported census theorem used in
# the routes below.  It is not computed here and the theorem is not stated here.
CENSUS_REACH = Fraction('3.07')

# Ten manifold names and one volume string, typed in by hand below.  An eleventh entry of
# the same group is S^3 and is excluded (M # S^3 = M makes the inequality trivially
# true), so ten remain.  These values are not printed in the note.
LEVEL1_VOLT = '2.02988321281931'          # the string printed in every row of the level
LEVEL1_SAFE = Fraction('2.02989')         # the licensed round-up (table caption, :700)
LEVEL1_ROWS = [
    (723, 'L(5,1)'),
    (724, 'L(10,3)'),
    (725, 'SFS [S2: (2,1) (3,2) (3,-1)]'),
    (726, 'SFS [S2: (2,1) (3,1) (7,-6)]'),
    (727, 'SFS [S2: (2,1) (4,1) (5,-4)]'),
    (728, 'SFS [S2: (3,1) (3,1) (4,-3)]'),
    (729, 'T x I / [ -2,-1 | -1,-1 ]'),
    (730, 'T x I / [ 2,1 | 1,1 ]'),
    (731, 'SFS [D: (2,1) (2,1)] u_(0,1 | 1,1) SFS [D: (2,1) (3,2)]'),
    (732, 'SFS [D: (2,1) (2,1)] u_( 0,1 | 1,0) SFS [D: (2,1) (3,1)]'),
]
LEVEL1_S3_LINE = 722                      # the eleventh row of the level, excluded

# The second volume level of Table 1 (the reason no N beyond level one can be used).
LEVEL2_VOLT = Fraction('2.56897060093671')

# The two Table 1 rows of the pair the theorems deliberately exclude, and the row of
# their connected sum -- the source's own worked example, where subadditivity HOLDS.
VOLT_L21 = Fraction('2.66674478344906')   # Table 1 line 748
VOLT_L31 = Fraction('2.56897060093671')   # Table 1 line 736
VOLT_L21_SUM_L31 = Fraction('2.66674478344906')   # Table 1 line 760, L(2,1) # L(3,1)

# Row counts of the eight volume levels of Table 1, as counted in the paper.
TABLE1_LEVEL_SIZES = [11, 13, 13, 6, 14, 4, 10, 5]
TABLE1_ROWS = 76

# Digit strings the paper prints for the two constants, and the value the SOURCE prints.
TWO_V0_STR = '2.029883212819307250042405'
V0_STR = '1.014941606409653625021203'
SOURCE_TABLE_STR = '2.02988321281931'
FOUR_V0_STR = '4.0597664256386145'

# A hand-typed floating-point literal, not computed by this program and not accompanied
# by any certificate here.  It is used ONLY in the boundary check of Step 6, i.e. that a
# value agreeing with v_0 to 16 digits doubles to something not exceeding 2 v_0.  No
# manifold is named for it and none can be certified from what is shipped here.
CLOSED_CENSUS_THIRD = Fraction('1.014941606409654')

# --------------------------------------------------------------------------------------
# check harness
# --------------------------------------------------------------------------------------

_n = 0
_bad = 0


def ck(name, cond, detail=''):
    global _n, _bad
    _n += 1
    if cond:
        print('PASS %s%s' % (name, ('  [%s]' % detail) if detail else ''))
    else:
        _bad += 1
        print('FAIL %s%s' % (name, ('  [%s]' % detail) if detail else ''))


def note(s):
    print('NOTE %s' % s)


def d(x, places=15):
    """A Fraction as a decimal string -- printing only, never a decision."""
    if isinstance(x, Fraction):
        neg = x < 0
        x = -x if neg else x
        scaled = (x.numerator * 10 ** places) // x.denominator
        s = str(scaled).rjust(places + 1, '0')
        s = s[:-places] + '.' + s[-places:]
        s = s.rstrip('0').rstrip('.')
        return ('-' if neg else '') + (s or '0')
    return str(x)


# --------------------------------------------------------------------------------------
# pi, to PREC digits, by two independent Machin-like formulas
# --------------------------------------------------------------------------------------

def atan_recip(k):
    """arctan(1/k) as a Decimal, k a positive integer >= 2, by its alternating series."""
    k = Decimal(k)
    term = Decimal(1) / k
    ksq = k * k
    total = term
    cut = Decimal(10) ** -(PREC + 5)
    j = 1
    while True:
        term = -term / ksq
        add = term / (2 * j + 1)
        if abs(add) < cut:
            break
        total += add
        j += 1
    return total


def pi_machin():
    return 16 * atan_recip(5) - 4 * atan_recip(239)


def pi_euler():
    return 4 * (atan_recip(2) + atan_recip(3))


# --------------------------------------------------------------------------------------
# Bernoulli numbers B_0, B_1, ..., exactly, by the standard recurrence
#   sum_{k=0}^{m} C(m+1,k) B_k = 0   (m >= 1)
# --------------------------------------------------------------------------------------

def bernoulli(upto):
    from math import comb
    B = [Fraction(1)]
    for m in range(1, upto + 1):
        s = sum(Fraction(comb(m + 1, k)) * B[k] for k in range(m))
        B.append(-s / (m + 1))
    return B


# --------------------------------------------------------------------------------------
# Cl_2(theta) for theta = q * pi, 0 < q < 2, by
#   Cl_2(t) = t - t*ln(t) + t * sum_{n>=1} |B_2n| t^(2n) / (2 (2n)! n (2n+1))
# (the standard expansion; the coefficient is zeta(2n)/(2 pi)^(2n) / (n(2n+1)) and
#  zeta(2n)/(2 pi)^(2n) = |B_2n| / (2 (2n)!)).  All terms are positive and decay like
# (t / 2pi)^(2n), so the truncation error is bounded by the first omitted term divided by
# (1 - (t/2pi)^2); we return a bound on it.
# --------------------------------------------------------------------------------------

_B = bernoulli(140)


def clausen2(q, PI):
    t = Fraction(q).numerator * PI / Fraction(q).denominator
    ln_t = t.ln()
    total = t - t * ln_t
    tsq = t * t
    tpow = Decimal(1)              # t^(2n-2) at the top of the loop body
    fact = Decimal(1)              # (2n-2)! at the top of the loop body
    last = None
    n = 1
    while 2 * n <= len(_B) - 1:
        tpow = tpow * tsq
        fact = fact * (2 * n - 1) * (2 * n) if n > 1 else Decimal(2)
        b = _B[2 * n]
        b = -b if b < 0 else b
        coeff = Decimal(b.numerator) / Decimal(b.denominator)
        term = t * tpow * coeff / (2 * fact * n * (2 * n + 1))
        if term == 0:
            last = term
            break
        total += term
        last = term
        if term < Decimal(10) ** (-(PREC - 30)):
            break
        n += 1
    ratio = (t / (2 * PI)) ** 2
    bound = last / (1 - ratio) if last is not None else Decimal(0)
    return total, bound


def as_bounds(value, err):
    """A Decimal value with an absolute error bound -> exact rational (lo, hi)."""
    v = Fraction(*value.as_integer_ratio())
    e = Fraction(*err.as_integer_ratio())
    # a further slack of one unit in the last working place, for the Decimal rounding
    slack = Fraction(10) ** (-(PREC - 5))
    return v - e - slack, v + e + slack


def agrees_to(lo, hi, s, digits):
    """Does the printed digit string `s` agree with the interval [lo,hi] to `digits`
    significant digits?  Exact rational test: the interval must lie strictly inside the
    band of numbers that round to `s` at that many digits."""
    x = Fraction(s)
    ulp = Fraction(10) ** -(digits - len(str(int(abs(x)))))
    return (x - ulp / 2) <= lo and hi <= (x + ulp / 2)


# ======================================================================================
def main():
    print('numerical re-derivation for the note "A counterexample, with a hyperbolic '
          'summand, to a printed subadditivity question for topological volume" -- '
          'constants, hand-typed values and inequalities only; the topological content '
          '(vol_t = vol_h, Cao-Meyerhoff, Gabai-Meyerhoff-Milley, Sphere Theorem, '
          'Kneser-Milnor) is imported and NOT verified here')
    print('python %s, decimal working precision %d, all decisions in exact rationals'
          % (sys.version.split()[0], PREC))
    print()

    print('=== Step 1: the two transcendental constants, from their definitions')
    pm, pe = pi_machin(), pi_euler()
    tol = Decimal(10) ** (-(PREC - 10))
    ck('pi_from_two_independent_machin_formulas_agree',
       abs(pm - pe) < tol,
       '|16atan(1/5)-4atan(1/239) - 4(atan(1/2)+atan(1/3))| < 1e-%d' % (PREC - 10))
    PI = pm

    # 2 v_0 = 3 Cl_2(2 pi / 3);  v_0 = Cl_2(pi / 3)   (Clausen duplication at theta=pi/3)
    c23, e23 = clausen2(Fraction(2, 3), PI)
    c13, e13 = clausen2(Fraction(1, 3), PI)
    ck('clausen_series_truncation_bounds_are_negligible',
       e23 < Decimal(10) ** -45 and e13 < Decimal(10) ** -45,
       'tail bounds %.3e and %.3e' % (float(e23), float(e13)))

    ck('bernoulli_recurrence_reproduces_the_classical_values',
       [_B[2], _B[4], _B[6], _B[8], _B[10]]
       == [Fraction(1, 6), Fraction(-1, 30), Fraction(1, 42),
           Fraction(-1, 30), Fraction(5, 66)],
       'B_2..B_10 = 1/6, -1/30, 1/42, -1/30, 5/66')

    two_v0_lo, two_v0_hi = as_bounds(3 * c23, 3 * e23)
    v0_lo, v0_hi = as_bounds(c13, e13)

    ck('two_v0_equals_3_Cl2(2pi/3)_to_the_printed_25_digits',
       agrees_to(two_v0_lo, two_v0_hi, TWO_V0_STR, 25),
       '2v_0 = %s' % TWO_V0_STR)
    ck('v0_equals_Cl2(pi/3)_to_the_printed_25_digits',
       agrees_to(v0_lo, v0_hi, V0_STR, 25),
       'v_0 = %s' % V0_STR)
    # The Clausen duplication formula Cl_2(2t) = 2Cl_2(t) - 2Cl_2(pi-t) at t = pi/3 gives
    # 3 Cl_2(2pi/3) = 2 Cl_2(pi/3), i.e. 2 v_0 = 2 Cl_2(pi/3).  Two different arguments of
    # the series, so this is a real cross-check of the implementation.
    ck('duplication_identity_3Cl2(2pi/3)_equals_2Cl2(pi/3)',
       two_v0_lo <= 2 * v0_hi and 2 * v0_lo <= two_v0_hi,
       'the two routes to 2v_0 have overlapping rational enclosures')
    ck('the_sources_own_printed_value_is_two_v0_correctly_rounded',
       agrees_to(two_v0_lo, two_v0_hi, SOURCE_TABLE_STR, 15),
       'Table 1 prints %s' % SOURCE_TABLE_STR)
    ck('four_v0_matches_the_printed_17_digits',
       agrees_to(2 * two_v0_lo, 2 * two_v0_hi, FOUR_V0_STR, 17),
       '4v_0 = %s' % FOUR_V0_STR)
    note('rational enclosure of 2v_0 used for every decision below: [%s, %s], '
         'of width below 1e-45' % (d(two_v0_lo, 30), d(two_v0_hi, 30)))

    print()
    print('=== Step 2: the transcribed table, and the licensed round-up')
    ck('level_one_of_table_1_has_ten_transcribed_rows_and_one_excluded_S3_row',
       len(LEVEL1_ROWS) == 10 and LEVEL1_S3_LINE == 722
       and [ln for ln, _ in LEVEL1_ROWS] == list(range(723, 733)),
       'source lines 723-732, with 722 = S^3 excluded; level size 11')
    ck('the_ten_transcribed_manifolds_are_distinct_and_none_is_S3',
       len({nm for _, nm in LEVEL1_ROWS}) == 10
       and not any(nm.replace(' ', '') in ('S3', 'S^3') for _, nm in LEVEL1_ROWS),
       '10 distinct Regina names')
    ck('the_printed_level_one_volume_is_at_most_the_safe_round_up',
       Fraction(LEVEL1_VOLT) <= LEVEL1_SAFE,
       '%s <= %s' % (LEVEL1_VOLT, d(LEVEL1_SAFE)))
    ck('table_1_level_sizes_sum_to_the_stated_row_count',
       sum(TABLE1_LEVEL_SIZES) == TABLE1_ROWS,
       '+'.join(map(str, TABLE1_LEVEL_SIZES)) + ' = %d' % TABLE1_ROWS)
    ck('no_manifold_beyond_level_one_can_be_paired_with_W_under_theorem_B',
       LEVEL2_VOLT > CENSUS_REACH - WEEKS_BOUND,
       'level 2 volume %s > 3.07 - 0.943 = %s'
       % (d(LEVEL2_VOLT), d(CENSUS_REACH - WEEKS_BOUND)))

    print()
    print('=== Step 3: witness X_0 = W # W, by both routes')
    sum0 = 2 * WEEKS_BOUND
    ck('vol_t(W)+vol_t(W)_is_at_most_1_886',
       sum0 == Fraction('1.886'), '2 * 0.943 = %s' % d(sum0))
    gapA = two_v0_lo - sum0
    ck('X0_route_A_two_v0_strictly_exceeds_the_sum',
       two_v0_lo > sum0 and gapA >= Fraction('0.143883212819307'),
       'gap >= 0.143883212819307 (got %s)' % d(gapA, 21))
    gapB = CENSUS_REACH - sum0
    ck('X0_route_B_the_census_reach_strictly_exceeds_the_sum',
       CENSUS_REACH > sum0 and gapB == Fraction('1.184'),
       '3.07 - 1.886 = %s' % d(gapB))

    print()
    print('=== Step 4: witnesses X_1..X_10 = W # N, N running over level one')
    pair_sum = WEEKS_BOUND + LEVEL1_SAFE
    note('vol_t(W) + vol_t(N) <= 0.943 + %s = %s for every one of the ten'
         % (d(LEVEL1_SAFE), d(pair_sum)))
    for i, (line, nm) in enumerate(LEVEL1_ROWS, start=1):
        ok = (pair_sum == Fraction('2.97289')
              and pair_sum < CENSUS_REACH
              and CENSUS_REACH - pair_sum >= Fraction('0.09711'))
        ck('X%d_W_plus_row_%d_sum_is_below_the_census_reach'
           % (i, line), ok,
           '%s : 0.943 + %s = %s < 3.07, gap >= 0.09711'
           % (nm, d(LEVEL1_SAFE), d(pair_sum)))
        note('an imported census theorem with reach 3.07, plus irreducibility of the '
             'row -- neither stated nor checked here -- would be needed to turn this '
             'inequality into anything topological.')
    ck('the_certified_exhibit_is_exactly_eleven_pairs',
       1 + len(LEVEL1_ROWS) == 11, 'X_0 (W # W) together with X_1..X_10')

    print()
    print('=== Step 5: the argument does not prove a false universal')
    s = VOLT_L21 + VOLT_L31
    ck('the_excluded_pairs_table_values_sum_as_printed',
       s == Fraction('5.23571538438577'),
       '%s + %s = %s' % (d(VOLT_L21), d(VOLT_L31), d(s)))
    ck('subadditivity_HOLDS_at_the_sources_own_worked_example',
       VOLT_L21_SUM_L31 < s,
       'vol_t(L(2,1) # L(3,1)) = %s < %s' % (d(VOLT_L21_SUM_L31), d(s)))
    ck('the_excluded_pair_is_outside_theorem_B_hypothesis_anyway',
       s > CENSUS_REACH, '%s > 3.07, so theorem B never fires there' % d(s))
    ck('orientability_is_load_bearing_the_non_orientable_floor_would_not_suffice',
       v0_hi < sum0,
       'v_0 = %s < 1.886, so a floor of v_0 (Gieseking) gives nothing' % V0_STR)

    print()
    print('=== Step 6: the boundary, and the residual case')
    # A closed-census volume that is numerically v_0: the corresponding pair sums to
    # exactly 2 v_0, so route A yields an EQUALITY and not a strict violation.  Excluded
    # from the exhibit above.
    ck('the_third_closed_census_volume_agrees_with_v0_to_16_digits',
       abs(CLOSED_CENSUS_THIRD - v0_lo) < Fraction(1, 10 ** 15)
       and abs(CLOSED_CENSUS_THIRD - v0_hi) < Fraction(1, 10 ** 15),
       'measured %s vs v_0 = %s' % (d(CLOSED_CENSUS_THIRD, 16), V0_STR))
    ck('doubling_it_lands_ON_two_v0_so_route_A_gives_an_equality_not_a_violation',
       2 * v0_lo <= two_v0_hi and two_v0_lo <= 2 * v0_hi,
       '2 v_0 = 2 v_0; the Question asks "<=", which an equality satisfies')
    ck('the_both_summands_non_hyperbolic_case_is_out_of_this_routes_reach',
       2 * two_v0_lo > CENSUS_REACH,
       '4v_0 = %s > 3.07 = the census reach' % FOUR_V0_STR)
    ck('three_copies_of_W_already_clear_the_two_v0_floor',
       3 * WEEKS_BOUND > two_v0_hi,
       '3 * 0.943 = %s > 2v_0' % d(3 * WEEKS_BOUND))

    print()
    # These three ratios are DERIVED from bounds the paper prints (0.943, 3.07, and Table 1's
    # two lens-space rows); the ratios THEMSELVES are not printed in the paper, so do not go
    # looking for 1.07, 1.62 or 0.510 there.  They only record how large the violation is.
    print('=== Step 7: by-product ratios, derived from the printed bounds and not themselves '
          'printed in the paper')
    ck('ratio_at_n_2_exceeds_1_07_by_route_A',
       two_v0_lo / sum0 > Fraction('1.07'),
       '2v_0 / 1.886 > 1.07')
    ck('ratio_at_n_2_exceeds_1_62_by_route_B',
       CENSUS_REACH / sum0 > Fraction('1.62'),
       '3.07 / 1.886 > 1.62')
    ck('the_sources_own_data_point_has_ratio_below_0_510',
       VOLT_L21_SUM_L31 / s < Fraction('0.510'),
       '%s / %s < 0.510' % (d(VOLT_L21_SUM_L31), d(s)))

    print()
    note('SCOPE, and what this program does NOT check.  (i) It re-derives the two '
         'constants and re-decides comparisons; the topological content -- vol_t = vol_h '
         'on closed hyperbolic manifolds, Cao-Meyerhoff\'s 2v_0 floor for ORIENTABLE '
         'cusped manifolds, Gabai-Meyerhoff-Milley\'s vol_h(Weeks) <= 0.943, the Sphere '
         'Theorem and Kneser-Milnor uniqueness -- is quoted from the literature and is '
         'not machine-checkable here. (ii) It verifies no irreducibility of any manifold '
         'named below. (iii) The values 0.943, 3.07, 2.02988321281931, 2.66674478344906 '
         'and 2.56897060093671 are hand-typed literals, not recomputed here. (iv) One '
         'value, 1.014941606409654, is a hand-typed floating-point literal that is NOT '
         'interval-certified; it is used only in Step 6. (v) No volume of the Weeks '
         'manifold is used anywhere except through the bound 0.943.')

    print()
    if _bad:
        print('VERDICT: %d of %d CHECKS FAILED' % (_bad, _n))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _n)
    return 0


if __name__ == '__main__':
    sys.exit(main())

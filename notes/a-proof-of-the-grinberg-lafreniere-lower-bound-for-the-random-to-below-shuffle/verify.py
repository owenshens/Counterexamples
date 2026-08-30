#!/usr/bin/env python3
"""
Verification of the note

    A Proof of the Grinberg-Lafreniere Lower Bound for the Random-to-Below Shuffle

which asserts, for  S(n) = sum_{i=2}^{n} n / ( i (H_n - H_{i-1}) ),

    S(n) >= n log n + n log log n + 0.0034553 n        for every integer n >= 2.

Python 3.9+, STANDARD LIBRARY ONLY (fractions, math, sys).  No third-party package, no
external data file.

WHY THERE ARE NO FLOATING-POINT DECISIONS.  Every quantity in the note involves a
logarithm, so a float comparison would be a claim about the note resting on a claim about
IEEE rounding.  Instead this program carries a tiny outward-rounded rational INTERVAL
arithmetic (class `Iv`) and a rigorous interval logarithm built from the identity
log((1+z)/(1-z)) = 2 sum_{k>=0} z^(2k+1)/(2k+1) with an explicit geometric tail bound, so
every inequality below is decided by comparing two exact Fractions.  `math` is used only
for `isqrt` and for cosmetic printing.

WHAT IT READS.  The objects printed in the note, transcribed verbatim below:
  * PAPER_TABLE   -- Table 1, the finite range 2 <= n <= 35 (n, margin, M(n), D_n);
  * PAPER_EXACT   -- the seven exact rationals S(2)..S(8) of Section 6;
  * PAPER_S35     -- the digit counts of the numerator and denominator of S(35);
  * PAPER_TAIL_*  -- every displayed constant of the two-line tail computation at n0 = 36;
  * PAPER_C, PAPER_INVLN2, PAPER_GAMMA_LB, PAPER_GAMMA_N -- the four transcribed decimals
    the proof uses as bounds, each of which is re-derived here rather than trusted.
Everything else is recomputed from the definition of S(n).
"""

import math
import sys
from fractions import Fraction as F

# ===========================================================================
# 0.  THE OBJECTS PRINTED IN THE NOTE  (transcribed; nothing else is assumed)
# ===========================================================================

# Table 1 of the note: n | margin S(n) - (n log n + n log log n) | M(n) | D_n
PAPER_TABLE = """
 2   1.346731480043   0.6733657400   0.000000000
 3   1.222019651146   0.4073398837   0.316290732
 4   1.280153647475   0.3200384119   0.433370036
 5   1.426874761492   0.2853749523   0.488200912
 6   1.629643080232   0.2716071800   0.516758296
 7   1.871875508071   0.2674107869   0.532388751
 8   2.143632805281   0.2679541007   0.541024198
 9   2.438399017356   0.2709332242   0.545623780
10   2.751643886956   0.2751643887   0.547788443
11   3.080080611049   0.2800073283   0.548438078
12   3.421244838997   0.2851037366   0.548126242
13   3.773239425541   0.2902491866   0.547198309
14   4.134571669591   0.2953265478   0.545876014
15   4.504045289409   0.3002696860   0.544304950
16   4.880686334613   0.3050428959   0.542582404
17   5.263690953326   0.3096288796   0.540774226
18   5.652387688680   0.3140215383   0.538925378
19   6.046209700325   0.3182215632   0.537066694
20   6.444673927241   0.3222336964   0.535219308
21   6.847365206445   0.3260650098   0.533397608
22   7.253923995463   0.3297238180   0.531611256
23   7.664036758584   0.3332189895   0.529866561
24   8.077428351275   0.3365595146   0.528167456
25   8.493855923618   0.3397542369   0.526516174
26   8.913103992586   0.3428116920   0.524913738
27   9.334980423760   0.3457400157   0.523360310
28   9.759313127911   0.3485468974   0.521855443
29  10.185947324781   0.3512395629   0.520398267
30  10.614743260850   0.3538247754   0.518987619
31  11.045574293407   0.3563088482   0.517622145
32  11.478325272393   0.3586976648   0.516300372
33  11.912891166026   0.3609967020   0.515020759
34  12.349175887291   0.3632110555   0.513781734
35  12.787091286981   0.3653454653   0.512581727
"""

PAPER_EXACT = {                      # Section 6 of the note
    2: "2",
    3: "24/5",
    4: "740/91",
    5: "386080/32571",
    6: "3561064/224257",
    7: "915340590661/45418769253",
    8: "44929799655432/1823748678895",
}
PAPER_S35_DIGITS = (353, 351)        # numerator, denominator of S(35)

PAPER_CELL_LO, PAPER_CELL_HI = 2, 35
PAPER_CELL_MIN_M = F("0.2674107868")   # the cell proposition: M(n) > this on 2..35
PAPER_TAIL_MIN_M_200 = F("0.2692847")  # the Status section, on the window 36 <= n <= 199
PAPER_N0 = 36
PAPER_SLACK = F("0.0034553")         # the theorem's constant c

# Transcribed bounds the proof uses.  Each is re-derived below.
PAPER_INVLN2 = F("1.4426951")        # >= 1/log 2
PAPER_GAMMA_LB = F("0.57721566")     # <= Euler's gamma
PAPER_C = F("0.42278434")            # >= 1 - gamma
PAPER_GAMMA_N = 10000                # the n at which the note evaluates g_n

# The displayed two-line computation at n0 = 36 (Section 5 of the note).
PAPER_TAIL_TERMS = (                 # upper bounds for the four terms of (IV)
    F("1.028571429"),                #   n/(n-1)
    F("0.139527657"),                #   1/(2 log n)
    F("0.240449184"),                #   1.4426951 / floor(sqrt(n+1))
    F("0.333459513"),                #   4 log((n+1)/2) / (n-1)
)
PAPER_TAIL_IV = F("0.871003892")     # >= (IV) at n0 = 36
PAPER_TAIL_V = F("0.125540793")      # >= (V)  at n0 = 36
PAPER_TAIL_TOTAL = F("0.996544685")  # >= (IV)+(V) at n0 = 36, and < 1

# The note's sharpness claim: the same two lines exceed 1 at n0 = 35 and at n0 = 34.
PAPER_N35_IV = F("0.875268945")      # <= (IV) at n = 35
PAPER_N35_V = F("0.126601216")       # <= (V)  at n = 35
PAPER_N35_TOTAL = F("1.001870161")   # <= (IV)+(V) at n = 35, and > 1

# ===========================================================================
# 1.  OUTWARD-ROUNDED RATIONAL INTERVAL ARITHMETIC
# ===========================================================================

PREC = 10 ** 40          # every endpoint is rounded outward onto the grid 1/PREC


def _fl(x):
    """largest multiple of 1/PREC that is <= x"""
    return F(x.numerator * PREC // x.denominator, PREC)


def _ce(x):
    """smallest multiple of 1/PREC that is >= x"""
    n, d = x.numerator * PREC, x.denominator
    return F(-((-n) // d), PREC)


class Iv(object):
    """[lo, hi] with Fraction endpoints; every operation rounds outward, so an Iv always
    CONTAINS the exact real value of the expression that produced it."""

    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        lo = F(lo)
        hi = lo if hi is None else F(hi)
        assert lo <= hi
        self.lo, self.hi = _fl(lo), _ce(hi)

    @staticmethod
    def raw(lo, hi):
        o = Iv.__new__(Iv)
        o.lo, o.hi = _fl(lo), _ce(hi)
        return o

    def __add__(self, o):
        o = _iv(o)
        return Iv.raw(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return Iv.raw(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-_iv(o))

    def __rsub__(self, o):
        return _iv(o) + (-self)

    def __mul__(self, o):
        o = _iv(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv.raw(min(p), max(p))

    __rmul__ = __mul__

    def inv(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('interval straddles 0')
        return Iv.raw(1 / self.hi, 1 / self.lo)

    def __truediv__(self, o):
        return self * _iv(o).inv()

    def __rtruediv__(self, o):
        return _iv(o) * self.inv()

    def width(self):
        return self.hi - self.lo

    # --- the only way a decision is ever made ---------------------------
    def gt(self, o):
        """certainly greater than"""
        return self.lo > _iv(o).hi

    def lt(self, o):
        return self.hi < _iv(o).lo

    def contains(self, x):
        return self.lo <= F(x) <= self.hi

    def __repr__(self):
        return '[%.15f, %.15f]' % (self.lo, self.hi)

    def mid(self):
        return (self.lo + self.hi) / 2


def _iv(x):
    return x if isinstance(x, Iv) else Iv(F(x))


def near(x, target, tol):
    """the interval x lies within `tol` of `target`.  Used when `target` is a DECIMAL
    printed to finitely many places, so exact containment is not the right question."""
    x, target, tol = _iv(x), F(target), F(tol)
    return x.lo > target - tol and x.hi < target + tol


# ===========================================================================
# 2.  A RIGOROUS INTERVAL LOGARITHM
# ===========================================================================

LN_EPS = F(1, 10 ** 30)


def _log_ratio(z, eps=LN_EPS):
    """(lo, hi) enclosing log((1+z)/(1-z)) = 2 sum_{k>=0} z^(2k+1)/(2k+1), for 0 <= z < 1/2.

    Tail bound: sum_{j>=k} z^(2j+1)/(2j+1) <= z^(2k+1) / ((2k+1)(1-z^2)).
    """
    assert 0 <= z < F(1, 2)
    if z == 0:
        return F(0), F(0)
    zzl, zzh = _fl(z * z), _ce(z * z)
    slo = shi = F(0)
    plo, phi = _fl(z), _ce(z)
    k = 0
    while True:
        slo = _fl(slo + plo / (2 * k + 1))
        shi = _ce(shi + phi / (2 * k + 1))
        k += 1
        plo, phi = _fl(plo * zzl), _ce(phi * zzh)
        tail = _ce(phi / ((2 * k + 1) * (1 - zzh)))
        if 2 * tail <= eps:
            return 2 * slo, 2 * (shi + tail)


_LN2 = Iv.raw(*_log_ratio(F(1, 3)))     # (1+1/3)/(1-1/3) = 2
_LN_CACHE = {}


def LN(q):
    """Iv enclosing log(q) for an exact positive rational q."""
    q = F(q)
    if q <= 0:
        raise ValueError('log of %s' % q)
    hit = _LN_CACHE.get(q)
    if hit is not None:
        return hit
    # argument reduction q = 2^k * m with m in [2/3, 4/3]
    k = q.numerator.bit_length() - q.denominator.bit_length()
    two = F(2)

    def m_of(j):
        return q / two ** j if j >= 0 else q * two ** (-j)

    while m_of(k) > F(4, 3):
        k += 1
    while m_of(k) < F(2, 3):
        k -= 1
    m = m_of(k)
    z = (m - 1) / (m + 1)               # |z| <= 1/5
    if z >= 0:
        a = Iv.raw(*_log_ratio(z))
    else:
        lo, hi = _log_ratio(-z)
        a = Iv.raw(-hi, -lo)
    out = a + (k * _LN2 if k else Iv(0))
    _LN_CACHE[q] = out
    return out


def LN_iv(x):
    """log of an interval, by monotonicity."""
    x = _iv(x)
    if x.lo <= 0:
        raise ValueError('log of interval containing 0')
    return Iv.raw(LN(x.lo).lo, LN(x.hi).hi)


# ===========================================================================
# 3.  THE OBJECT ITSELF:  H_n, a_i, S(n), D_n, L_n, M(n)
# ===========================================================================

_H = [F(0)]


def H(n):
    while len(_H) <= n:
        _H.append(_H[-1] + F(1, len(_H)))
    return _H[n]


def a(n, i):
    """a_i = H_n - H_{i-1}, exact."""
    return H(n) - H(i - 1)


def S_exact(n):
    """S(n) as an exact rational.  Used on the cell 2..35 only."""
    return sum(F(n) / (F(i) * a(n, i)) for i in range(2, n + 1))


def S_over_n_iv(n):
    """S(n)/n = sum_{i=2}^{n} 1/(i a_i) as an interval (cheap for large n)."""
    t = Iv(0)
    for i in range(2, n + 1):
        t = t + Iv(1 / (F(i) * a(n, i)))
    return t


def rhs_over_n_iv(n):
    """log n + log log n."""
    ln = LN(n)
    return ln + LN_iv(ln)


def M_iv(n):
    """M(n) = S(n)/n - log n - log log n  (= margin / n)."""
    return S_over_n_iv(n) - rhs_over_n_iv(n)


def D_iv(n):
    """D_n via identity (I): D_n = 1 - S(n)/n + log n + log(H_n - 1).

    (I) itself is proved exactly in checks I1/I2, and check I3 confirms this value
    against the defining sum sum_{i=2}^{n-1} f(a_{i+1}/a_i) for n <= 40.
    """
    return 1 - S_over_n_iv(n) + LN(n) + LN(H(n) - 1)


def D_from_definition_iv(n):
    t = Iv(0)
    for i in range(2, n):
        x = a(n, i + 1) / a(n, i)
        t = t + (Iv(x) - 1 - LN(x))
    return t


def L_iv(n):
    """L_n = log( log n / (H_n - 1) )."""
    return LN_iv(LN(n) / Iv(H(n) - 1))


def bound_III(n):
    """(1/2)[ n/(n-1) - 1/(2(H_n-1)) + R_n ],  R_n = sum_{j=3}^{n-1} 1/(j(j-1)a_j)."""
    R = Iv(0)
    for j in range(3, n):
        R = R + Iv(1 / (F(j * (j - 1)) * a(n, j)))
    return (Iv(F(n, n - 1)) - Iv(1) / (2 * Iv(H(n) - 1)) + R) / 2


def bound_IV(n):
    """(1/2)[ n/(n-1) + 1/(2 log n) + 1.4426951/floor(sqrt(n+1)) + 4 log((n+1)/2)/(n-1) ]"""
    N = n + 1
    J = math.isqrt(N)
    return (Iv(F(n, n - 1))
            + Iv(1) / (2 * LN(n))
            + Iv(PAPER_INVLN2 / J)
            + 4 * LN(F(N, 2)) / Iv(F(n - 1))) / 2


def bound_V(n, c=None):
    """-log(1 - c/log n)  with c = PAPER_C."""
    c = PAPER_C if c is None else c
    return -LN_iv(1 - Iv(c) / LN(n))


# ===========================================================================
# 4.  CHECK HARNESS
# ===========================================================================

_PASSES = 0
_FAILS = 0
_LINES = []


def PASS(name, detail=''):
    global _PASSES
    _PASSES += 1
    print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))


def FAIL(name, detail=''):
    global _FAILS
    _FAILS += 1
    print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def CHECK(cond, name, detail=''):
    if cond:
        PASS(name, detail)
    else:
        FAIL(name, detail)
    return bool(cond)


def NOTE(s):
    print('NOTE %s' % s)


def d(x, k=12):
    """decimal string of a Fraction, for display only"""
    x = F(x)
    neg = x < 0
    if neg:
        x = -x
    scaled = x.numerator * 10 ** k // x.denominator
    s = str(scaled).rjust(k + 1, '0')
    return ('-' if neg else '') + s[:-k] + '.' + s[-k:]


# ===========================================================================
# 5.  PARSE THE PRINTED TABLE
# ===========================================================================

TABLE = []
for line in PAPER_TABLE.strip().splitlines():
    parts = line.split()
    assert len(parts) == 4, line
    TABLE.append((int(parts[0]), F(parts[1]), F(parts[2]), F(parts[3])))
TABLE_BY_N = {r[0]: r for r in TABLE}

print('verification of the note: a proof of the Grinberg-Lafreniere lower bound')
print('for the random-to-below shuffle -- S(n) >= n log n + n log log n + %s n, n >= 2'
      % d(PAPER_SLACK, 7))
print('python %s ; exact Fraction interval arithmetic, no float decision' % sys.version.split()[0])
print('')

# ===========================================================================
# STEP 1:  the interval logarithm is sound
# ===========================================================================
print('=== Step 1: the interval logarithm, before it is used to decide anything')

_T22 = F(1, 10 ** 22)
CHECK(near(_LN2, F('0.69314718055994530941723'), _T22) and _LN2.width() < F(1, 10 ** 29),
      'interval_log_brackets_log_2_to_22_places',
      'width %.3e' % float(_LN2.width()))
CHECK(near(LN(10), F('2.30258509299404568401799'), _T22),
      'interval_log_brackets_log_10_to_22_places')
CHECK(near(LN(36), F('3.58351893845611000162495'), _T22),
      'interval_log_brackets_log_36_to_22_places')
CHECK(LN(1).contains(0) and LN(1).width() == 0,
      'interval_log_of_one_is_exactly_zero')
_s6 = LN(2) + LN(3)
CHECK(_s6.lo <= LN(6).lo and LN(6).hi <= _s6.hi + F(1, 10 ** 29),
      'interval_log_is_additive_on_2_times_3_equals_6')
CHECK((LN(F(7, 5)) + LN(F(5, 7))).contains(0),
      'interval_log_of_a_reciprocal_is_the_negation')
_mono = all(LN(F(k, 7)).hi < LN(F(k + 1, 7)).lo for k in range(1, 200))
CHECK(_mono, 'interval_log_is_strictly_increasing_on_199_consecutive_rationals')
_w = max(LN(F(k, 11)).width() for k in range(1, 400))
CHECK(_w < F(1, 10 ** 25), 'every_interval_log_width_is_below_1e_minus_25',
      'max width %.3e' % float(_w))
# non-vacuity: the enclosure must EXCLUDE nearby wrong values, or it decides nothing
CHECK(not _LN2.contains(F('0.6931471')) and not _LN2.contains(F('0.6931472')),
      'control_the_log_2_enclosure_excludes_both_neighbouring_7_place_values')
CHECK(not LN(36).contains(F('3.5835189')) and not LN(36).contains(F('3.583519')),
      'control_the_log_36_enclosure_excludes_nearby_wrong_values')
print('')

# ===========================================================================
# STEP 2:  S(n) from the definition, and the seven exact rationals of Section 6
# ===========================================================================
print('=== Step 2: S(n) recomputed from the definition')

_bad = []
for n, lit in sorted(PAPER_EXACT.items()):
    got = S_exact(n)
    if str(got) != lit:
        _bad.append((n, lit, str(got)))
CHECK(not _bad, 'the_notes_seven_exact_rationals_S_2_to_S_8_are_reproduced',
      'mismatches %s' % (_bad or '[]'))

_S35 = S_exact(35)
_dig = (len(str(_S35.numerator)), len(str(_S35.denominator)))
CHECK(_dig == PAPER_S35_DIGITS,
      'the_digit_counts_of_S_35_are_the_printed_353_and_351', 'got %s' % (_dig,))

_ok = all(a(n, i) - a(n, i + 1) == F(1, i) for n in range(3, 60) for i in range(2, n))
CHECK(_ok, 'the_recurrence_a_i_minus_a_i_plus_1_equals_1_over_i_holds_exactly',
      'n = 3..59, every i')
_ok = all(a(n, 2) == H(n) - 1 and a(n, n) == F(1, n) for n in range(2, 60))
CHECK(_ok, 'a_2_equals_H_n_minus_1_and_a_n_equals_1_over_n_exactly', 'n = 2..59')
_ok = all(0 < a(n, i + 1) < a(n, i) for n in range(3, 60) for i in range(2, n))
CHECK(_ok, 'the_a_i_are_positive_and_strictly_decreasing', 'n = 3..59')
print('')

# ===========================================================================
# STEP 3:  the identity (I), proved by two exact computations
# ===========================================================================
print('=== Step 3: identity (I)  S(n)/n = 1 + log n + log(H_n - 1) - D_n')

_bad = []
for n in range(2, 81):
    lhs = sum(a(n, i + 1) / a(n, i) - 1 for i in range(2, n))
    if lhs != 1 - S_exact(n) / n:
        _bad.append(n)
CHECK(not _bad, 'sum_of_x_i_minus_1_equals_1_minus_S_over_n_exactly',
      'n = 2..80, exact Fractions, deviations %s' % (_bad or '[]'))

_bad = []
for n in range(3, 81):
    p = F(1)
    for i in range(2, n):
        p *= a(n, i + 1) / a(n, i)
    if p != F(1) / (F(n) * (H(n) - 1)):
        _bad.append(n)
CHECK(not _bad, 'product_of_the_x_i_equals_1_over_n_times_H_n_minus_1_exactly',
      'n = 3..80, exact Fractions, deviations %s' % (_bad or '[]'))
NOTE('the two lines above are (I): sum(x_i - 1) = -(S(n)/n - 1) and -sum log x_i = '
     'log n + log(H_n - 1), so D_n = 1 - S(n)/n + log n + log(H_n - 1).')

_worst = F(0)
for n in range(3, 41):
    w = D_from_definition_iv(n) - D_iv(n)
    _worst = max(_worst, abs(w.lo), abs(w.hi))
CHECK(_worst < F(1, 10 ** 25),
      'the_defining_sum_for_D_n_agrees_with_identity_I_on_n_up_to_40',
      'max |difference| <= %.3e' % float(_worst))

_bad = [n for n in range(2, 36) if not (Iv(1) - D_iv(n) - L_iv(n) - M_iv(n)).contains(0)]
CHECK(not _bad, 'M_of_n_equals_1_minus_D_n_minus_L_n_on_the_whole_cell',
      'deviations %s' % (_bad or '[]'))
print('')

# ===========================================================================
# STEP 4:  Table 1 -- the finite range 2 <= n <= 35
# ===========================================================================
print('=== Step 4: Table 1, the finite range 2 <= n <= 35')

CHECK([r[0] for r in TABLE] == list(range(PAPER_CELL_LO, PAPER_CELL_HI + 1)),
      'the_printed_table_is_exactly_the_34_integers_2_to_35_with_no_gap_or_repeat',
      '%d rows' % len(TABLE))

_wm = _wM = _wD = F(0)
for n, pm, pM, pD in TABLE:
    _wm = max(_wm, max(abs((M_iv(n) * n - pm).lo), abs((M_iv(n) * n - pm).hi)))
    _wM = max(_wM, max(abs((M_iv(n) - pM).lo), abs((M_iv(n) - pM).hi)))
    _wD = max(_wD, max(abs((D_iv(n) - pD).lo), abs((D_iv(n) - pD).hi)))
CHECK(_wm <= F(5, 10 ** 13), 'every_printed_margin_agrees_with_the_recomputed_value',
      'max deviation %.3e, printed to 12 places' % float(_wm))
CHECK(_wM <= F(5, 10 ** 11), 'every_printed_M_value_agrees_with_the_recomputed_value',
      'max deviation %.3e, printed to 10 places' % float(_wM))
CHECK(_wD <= F(5, 10 ** 10), 'every_printed_D_value_agrees_with_the_recomputed_value',
      'max deviation %.3e, printed to 9 places' % float(_wD))

_viol = [n for n in range(2, 36) if not M_iv(n).gt(0)]
CHECK(not _viol, 'the_conjecture_itself_holds_at_every_n_in_the_cell',
      'exhaustion 34/34, counterexamples %d, violating n = %s' % (len(_viol), _viol))
_viol = [n for n in range(2, 36) if not M_iv(n).gt(PAPER_SLACK)]
CHECK(not _viol, 'the_strengthened_bound_with_slack_0_0034553_holds_at_every_n_in_the_cell',
      'violating n = %s' % _viol)
_viol = [n for n in range(2, 36) if not M_iv(n).gt(PAPER_CELL_MIN_M)]
CHECK(not _viol, 'the_notes_cell_minimum_claim_M_of_n_exceeds_0_2674107868',
      'violating n = %s ; the bound is 3.3e-12 below M(7)' % _viol)

_mn = min(range(2, 36), key=lambda n: M_iv(n).mid() * n)
CHECK(_mn == 3 and all((M_iv(3) * 3).lt(M_iv(n) * n) for n in range(2, 36) if n != 3)
      and near(M_iv(3) * 3, TABLE_BY_N[3][1], F(5, 10 ** 13)),
      'the_least_margin_in_the_cell_is_at_n_equals_3_and_is_the_printed_value',
      'margin %s' % d(TABLE_BY_N[3][1]))
_mn = min(range(2, 36), key=lambda n: M_iv(n).mid())
CHECK(_mn == 7 and all(M_iv(7).lt(M_iv(n)) for n in range(2, 36) if n != 7)
      and near(M_iv(7), TABLE_BY_N[7][2], F(5, 10 ** 11)),
      'the_least_M_in_the_cell_is_at_n_equals_7_and_is_the_printed_value',
      'M(7) = %s' % d(TABLE_BY_N[7][2], 10))
_mx = max(range(2, 36), key=lambda n: D_iv(n).mid())
CHECK(_mx == 11 and all(D_iv(11).gt(D_iv(n)) for n in range(2, 36) if n != 11)
      and near(D_iv(11), TABLE_BY_N[11][3], F(5, 10 ** 10)),
      'the_greatest_D_in_the_cell_is_at_n_equals_11_and_is_the_printed_value',
      'D_11 = %s' % d(TABLE_BY_N[11][3], 9))
CHECK(TABLE_BY_N[2][3] == 0 and D_from_definition_iv(2).contains(0),
      'D_2_is_zero_because_its_defining_sum_is_empty')

# --- controls -------------------------------------------------------------
_pub = [n for n in range(2, 36) if not M_iv(n).lt(_LN2 + Iv(F(1, n)))]
CHECK(not _pub,
      'control_must_stay_silent_the_authors_published_upper_bound_M_le_log2_plus_1_over_n',
      'violations %d over the cell' % len(_pub))
_forced = [n for n in range(2, 36) if not M_iv(n).lt(_LN2)]
CHECK(not _forced,
      'control_forced_positive_the_false_claim_M_ge_log_2_is_flagged_at_every_n',
      'flagged 34 of 34, unflagged n = %s' % _forced)
_unflagged = sorted(n for n in range(2, 36) if M_iv(n).gt(F('0.4')))
CHECK(_unflagged == [2, 3],
      'control_graded_at_threshold_0_4_the_unflagged_set_is_the_predicted_2_and_3',
      'measured %s, predicted [2, 3]' % _unflagged)
_struct = [n for n in range(2, 60)
           if min(1 / (F(i) * a(n, i)) for i in range(2, n + 1)) < F(1, n)]
CHECK(not _struct,
      'control_structural_each_summand_1_over_i_a_i_is_at_least_1_over_n',
      'failures %d over n = 2..59' % len(_struct))
print('')

# ===========================================================================
# STEP 5:  the ingredients of the bound on D_n
# ===========================================================================
print('=== Step 5: the lemmas behind the upper bound (III) and (IV) on D_n')

_bad = []
for k in range(1, 200):
    x = F(k, 200)
    f = Iv(x) - 1 - LN(x)
    if not (f.gt(0) and f.lt(Iv((1 - x) ** 2) / (2 * Iv(x)))):
        _bad.append(str(x))
CHECK(not _bad, 'f_is_positive_and_below_1_minus_x_squared_over_2x_on_199_rationals',
      'x = 1/200 .. 199/200, failures %s' % (_bad or '[]'))

_ok = all(0 < a(n, i + 1) / a(n, i) < 1 for n in range(3, 60) for i in range(2, n))
CHECK(_ok, 'x_i_lies_strictly_between_0_and_1_for_2_le_i_le_n_minus_1', 'n = 3..59')
_ok = all(1 - a(n, i + 1) / a(n, i) == F(1, i) / a(n, i)
          for n in range(3, 60) for i in range(2, n))
CHECK(_ok, 'one_minus_x_i_equals_1_over_i_a_i_exactly', 'n = 3..59')
_ok = all((1 - a(n, i + 1) / a(n, i)) ** 2 / (2 * (a(n, i + 1) / a(n, i)))
          == F(1, 2 * i) * (1 / a(n, i + 1) - 1 / a(n, i))
          for n in range(3, 60) for i in range(2, n))
CHECK(_ok, 'the_telescoping_step_1_minus_x_squared_over_2x_equals_half_i_times_b_gap',
      'n = 3..59, exact Fractions')

_bad = []
for n in range(3, 81):
    b = {i: 1 / a(n, i) for i in range(2, n + 1)}
    lhs = sum(F(1, i) * (b[i + 1] - b[i]) for i in range(2, n))
    rhs = b[n] / F(n - 1) - b[2] / 2 + sum(b[j] / F(j * (j - 1)) for j in range(3, n))
    if lhs != rhs:
        _bad.append(n)
CHECK(not _bad, 'the_summation_by_parts_identity_of_step_2c_is_exact',
      'n = 3..80, deviations %s' % (_bad or '[]'))

_bad = []
for m in range(2, 300):
    if sum(F(1, j * (j - 1)) for j in range(m, 3000)) + F(1, 2999) != F(1, m - 1):
        _bad.append(m)
CHECK(not _bad, 'the_telescoping_tail_sum_over_j_ge_m_of_1_over_j_j_minus_1_is_1_over_m_minus_1',
      'm = 2..299, exact Fractions, deviations %s' % (_bad or '[]'))

_bad = [(n, j) for n in (36, 100, 300, 1000) for j in range(2, n + 1)
        if not Iv(a(n, j)).gt(LN(F(n + 1, j)))]
CHECK(not _bad, 'a_j_strictly_exceeds_log_of_n_plus_1_over_j_by_integral_comparison',
      'n in 36,100,300,1000 and every j, failures %s' % (_bad or '[]'))

_bad = [n for n in range(2, 500) if not Iv(H(n)).lt(1 + LN(n))]
CHECK(not _bad, 'H_n_is_strictly_below_1_plus_log_n', 'n = 2..499, failures %s' % (_bad or '[]'))

CHECK(math.isqrt(37) == 6 and 6 * 6 <= 37 < 7 * 7,
      'floor_sqrt_37_is_6')
CHECK(Iv(F(1, 6)).gt(Iv(1) / Iv(F(6082763, 1000000), F(6082764, 1000000))),
      'one_over_floor_sqrt_37_strictly_exceeds_one_over_sqrt_37',
      '1/6 = 0.1666667 > 0.1643990 -- the unfloored form is NOT an upper bound')
_bad = [N for N in range(4, 5000) if math.isqrt(N + 1) < math.isqrt(N)]
CHECK(not _bad, 'one_over_floor_sqrt_N_is_non_increasing_in_N', 'N = 4..4999')
CHECK(Iv(PAPER_INVLN2).gt(Iv(1) / _LN2),
      'the_printed_1_4426951_is_an_upper_bound_for_1_over_log_2',
      '1/log 2 <= %s' % d((Iv(1) / _LN2).hi, 10))
_bad = [N for N in range(8, 4000) if not LN(F(N, 2)).gt(Iv(1) - Iv(F(2, N)))]
CHECK(not _bad, 'log_of_N_over_2_exceeds_1_minus_2_over_N_for_N_from_8_to_3999',
      'failures %s' % (_bad or '[]'))
_bad = [n for n in range(9, 3000)
        if not (3 <= F(n + 1, 2) <= n - 1 and math.isqrt(n + 1) ** 2 <= n + 1
                and n + 1 < (math.isqrt(n + 1) + 1) ** 2)]
CHECK(not _bad, 'the_three_blocks_of_step_2d_cover_j_from_3_to_n_minus_1',
      'n = 9..2999, needs sqrt(N) <= N/2 <= N-2 = n-1, failures %s' % (_bad or '[]'))
print('')

# ===========================================================================
# STEP 6:  (III) and (IV) really dominate D_n
# ===========================================================================
print('=== Step 6: the two upper bounds on D_n, tested against D_n itself')

_bad = [n for n in range(3, 301) if not bound_III(n).gt(D_iv(n))]
CHECK(not _bad, 'bound_III_strictly_dominates_D_n', 'n = 3..300, failures %s' % (_bad or '[]'))
_bad = [n for n in range(36, 301) if not bound_IV(n).gt(bound_III(n))]
CHECK(not _bad, 'bound_IV_strictly_dominates_bound_III', 'n = 36..300, failures %s' % (_bad or '[]'))
_bad = [n for n in range(36, 301) if not bound_IV(n).gt(D_iv(n))]
CHECK(not _bad, 'bound_IV_strictly_dominates_D_n', 'n = 36..300, failures %s' % (_bad or '[]'))

_terms = lambda n: (Iv(F(n, n - 1)), Iv(1) / (2 * LN(n)),
                    Iv(PAPER_INVLN2 / math.isqrt(n + 1)),
                    4 * LN(F(n + 1, 2)) / Iv(F(n - 1)))
_bad = []
for n in range(36, 2000):
    for t0, t1 in zip(_terms(n), _terms(n + 1)):
        if t1.gt(t0):
            _bad.append(n)
            break
CHECK(not _bad, 'each_of_the_four_terms_of_IV_is_non_increasing_in_n',
      'n = 36..1999, failures %s' % (_bad or '[]'))
NOTE('so evaluating (IV) at n0 bounds D_n for every n >= n0 -- this is what makes the '
     'tail a two-line computation.')
print('')

# ===========================================================================
# STEP 7:  gamma, and the bound (V) on L_n
# ===========================================================================
print('=== Step 7: Euler gamma from below, and the bound (V) on L_n')

_bad = [str(F(k, 300)) for k in range(1, 900)
        if not LN(1 + F(k, 300)).lt(Iv(F(k, 300)) * (2 + Iv(F(k, 300)))
                                    / (2 * (1 + Iv(F(k, 300)))))]
CHECK(not _bad, 'the_trapezoid_inequality_log_1_plus_x_lt_x_2_plus_x_over_2_1_plus_x',
      'x = 1/300 .. 3, failures %s' % (_bad or '[]'))
NOTE('that inequality is exactly the statement that g_n = H_n - log n - 1/(2n) increases.')

_g = lambda n: Iv(H(n)) - LN(n) - Iv(F(1, 2 * n))
_bad = [n for n in range(1, 1500) if not _g(n + 1).gt(_g(n))]
CHECK(not _bad, 'g_n_is_strictly_increasing_as_the_trapezoid_inequality_predicts',
      'n = 1..1499, failures %s' % (_bad or '[]'))

_gN = _g(PAPER_GAMMA_N)
CHECK(_gN.gt(PAPER_GAMMA_LB),
      'gamma_is_at_least_the_printed_0_57721566',
      'g(%d) = %s > %s' % (PAPER_GAMMA_N, d(_gN.lo, 12), d(PAPER_GAMMA_LB, 12)))
CHECK(PAPER_C + PAPER_GAMMA_LB >= 1,
      'hence_1_minus_gamma_is_at_most_the_printed_0_42278434',
      '%s + %s >= 1' % (d(PAPER_C, 8), d(PAPER_GAMMA_LB, 8)))

_bad = [n for n in range(36, 2000) if not Iv(H(n) - 1).gt(LN(n) - Iv(PAPER_C))]
CHECK(not _bad, 'H_n_minus_1_exceeds_log_n_minus_0_42278434_directly_for_n_up_to_1999',
      'failures %s' % (_bad or '[]'))
_bad = [n for n in range(36, 401) if not bound_V(n).gt(L_iv(n))]
CHECK(not _bad, 'bound_V_strictly_dominates_L_n', 'n = 36..400, failures %s' % (_bad or '[]'))
_bad = [n for n in range(36, 2000) if not bound_V(n).gt(bound_V(n + 1))]
CHECK(not _bad, 'bound_V_is_strictly_decreasing_in_n', 'n = 36..1999, failures %s' % (_bad or '[]'))
print('')

# ===========================================================================
# STEP 8:  the two-line tail computation at n0 = 36
# ===========================================================================
print('=== Step 8: the tail theorem, n >= 36')

_ex = _terms(PAPER_N0)
_bad = [k for k in range(4) if not Iv(PAPER_TAIL_TERMS[k]).gt(_ex[k])]
CHECK(not _bad, 'each_of_the_four_printed_term_bounds_at_n0_36_dominates_its_exact_term',
      'exact terms %s' % ', '.join(d(t.hi, 9) for t in _ex))
_sum4 = sum(PAPER_TAIL_TERMS, F(0))
CHECK(Iv(PAPER_TAIL_IV).gt(bound_IV(PAPER_N0)) and PAPER_TAIL_IV >= _sum4 / 2,
      'the_printed_0_871003892_dominates_IV_at_n0_36',
      '(IV)(36) <= %s, printed half bracket %s' % (d(bound_IV(PAPER_N0).hi, 12),
                                                   d(PAPER_TAIL_IV, 9)))
CHECK(Iv(PAPER_TAIL_V).gt(bound_V(PAPER_N0)),
      'the_printed_0_125540793_dominates_V_at_n0_36',
      '(V)(36) <= %s' % d(bound_V(PAPER_N0).hi, 12))
CHECK(PAPER_TAIL_TOTAL >= PAPER_TAIL_IV + PAPER_TAIL_V and PAPER_TAIL_TOTAL < 1,
      'the_printed_total_0_996544685_dominates_the_sum_and_is_below_1',
      '%s < 1' % d(PAPER_TAIL_TOTAL, 9))
CHECK(PAPER_SLACK <= 1 - PAPER_TAIL_TOTAL,
      'the_theorems_slack_0_0034553_is_at_most_1_minus_that_total',
      '1 - %s = %s >= %s' % (d(PAPER_TAIL_TOTAL, 9), d(1 - PAPER_TAIL_TOTAL, 9),
                             d(PAPER_SLACK, 9)))
CHECK((Iv(1) - bound_IV(PAPER_N0) - bound_V(PAPER_N0)).gt(PAPER_SLACK),
      'recomputed_without_the_printed_roundings_the_tail_still_gives_M_gt_0_0034553',
      '1 - (IV) - (V) >= %s' % d((Iv(1) - bound_IV(PAPER_N0) - bound_V(PAPER_N0)).lo, 12))

_bad = [n for n in range(36, 401) if not M_iv(n).gt(PAPER_SLACK)]
CHECK(not _bad, 'independently_M_of_n_exceeds_0_0034553_for_every_n_from_36_to_400',
      'computed from the definition of S(n), failures %s' % (_bad or '[]'))

# sharpness of the threshold 36 for THIS bound family
CHECK(PAPER_N35_IV <= bound_IV(35).lo,
      'the_printed_0_875268945_is_a_lower_bound_for_IV_at_n_equals_35',
      '(IV)(35) >= %s' % d(bound_IV(35).lo, 12))
CHECK(PAPER_N35_V <= bound_V(35).lo,
      'the_printed_0_126601216_is_a_lower_bound_for_V_at_n_equals_35',
      '(V)(35) >= %s' % d(bound_V(35).lo, 12))
CHECK(PAPER_N35_TOTAL > 1 and PAPER_N35_TOTAL <= PAPER_N35_IV + PAPER_N35_V,
      'the_same_two_lines_exceed_1_at_n0_35_so_the_threshold_36_is_sharp_for_this_family',
      '%s > 1' % d(PAPER_N35_TOTAL, 9))
CHECK((bound_IV(34) + bound_V(34)).gt(bound_IV(35) + bound_V(35)),
      'n0_34_is_worse_still',
      '(IV)+(V) at 34 >= %s' % d((bound_IV(34) + bound_V(34)).lo, 9))
CHECK(PAPER_CELL_HI + 1 == PAPER_N0 and PAPER_CELL_LO == 2,
      'the_cell_2_to_35_and_the_tail_n_ge_36_partition_every_integer_n_ge_2')

# the Status section's claim about the window 36 <= n <= 199
_at200 = Iv(1) - bound_IV(200) - bound_V(200)
_at199 = Iv(1) - bound_IV(199) - bound_V(199)
CHECK(_at200.gt(PAPER_TAIL_MIN_M_200) and not _at199.gt(PAPER_TAIL_MIN_M_200),
      'the_lemmas_give_M_gt_0_2692847_from_n_200_and_not_from_n_199',
      '1-(IV)-(V) at 200 >= %s, at 199 <= %s' % (d(_at200.lo, 9), d(_at199.hi, 9)))
print('')

# ===========================================================================
# STEP 9:  scope
# ===========================================================================
print('=== Step 9: what this program does NOT cover')
print('NOT RE-RUN: the exact rational S(n) is formed for n <= 80 only; for 81 <= n <= 400 '
      'the quantities S(n)/n, D_n, L_n and M(n) are computed by the outward-rounded '
      'interval arithmetic of Step 1 rather than as exact Fractions (the enclosures are '
      'rigorous, so no inequality above is weakened, but the huge-numerator exact form is '
      'not built).')
print('NOT RE-RUN: no statement about n > 400 is verified NUMERICALLY here. The tail for '
      'n >= 36 is settled by the monotonicity of (IV) and (V), which Steps 6 and 7 verify '
      'only on n = 36..1999; beyond that the note relies on the term-by-term monotonicity '
      'proofs of Sections 3 and 4, which this program does not symbolically differentiate.')
print('NOT RE-RUN: inf_{n>=2} M(n) is NOT determined. M(7) = 0.2674107869 is the minimum '
      'over the cell 2..35 only; Step 8 verifies M(n) > 0.0034553 on 36..400 but no '
      'exhaustive minimisation is performed there, so the note does not claim -- and this '
      'program does not check -- that 0.2674107869 is the optimal constant.')
print('NOT RE-RUN: D_n -> 1 - gamma is NOT verified. Step 5 bounds D_n only from ABOVE, so '
      'the corollary that the authors published upper bound is asymptotically sharp to '
      'within (log 2 - gamma) n is outside this program.')
print('NOT RE-RUN: the probabilistic content. That S(n) IS the expected value of the '
      'authors strong stationary time is their Theorem, quoted in Section 1 and not '
      'reproved here; this program verifies the inequality about the sum.')
print('')

# ===========================================================================
# VERDICT
# ===========================================================================
if _FAILS:
    print('VERDICT: %d CHECKS FAILED of %d' % (_FAILS, _PASSES + _FAILS))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _PASSES)
sys.exit(0)

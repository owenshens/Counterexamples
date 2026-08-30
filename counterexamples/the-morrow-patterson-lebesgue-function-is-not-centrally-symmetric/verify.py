#!/usr/bin/env python3
"""verify.py -- re-derives every quantity printed in paper.tex for

    "The Morrow-Patterson Lebesgue Function Is Not Centrally Symmetric"

from the objects the paper exhibits, in EXACT arithmetic, using the Python standard library only
(Python 3.9+; no numpy, no sympy, no external data file).

WHAT IS EXACT AND WHY NO FLOAT DECIDES ANYTHING
-----------------------------------------------
Every number that occurs at n = 2 lies in the real field  F = Q(sqrt2, sqrt5), which we represent as

    a + b*sqrt2 + c*sqrt5 + d*sqrt10 ,   a,b,c,d in Q   (fractions.Fraction)

1, sqrt2, sqrt5, sqrt10 are linearly independent over Q, so equality in F is equality of the four
rational coordinates -- decided exactly, never numerically. The only place an ORDER decision is
needed is the absolute value inside the Lebesgue function; there `sign()` brackets the element with
rational lower/upper bounds on sqrt2, sqrt5, sqrt10 obtained from `math.isqrt`, and refines the
precision until the bracket excludes 0. Since a nonzero element of F has nonzero coordinates, that
loop terminates, and the sign it returns is proved by the bracket, not guessed from a float.

The general-n part of the paper (the parity mechanism) is a statement about integer index sets, and
it is re-derived here with integers only, for every even n from 2 to 40.

Input: the block PAPER below is transcribed from paper.tex -- the six nodes, the six weights, the
six kernel entries, the six values on -MP_2, the witness and its value, and the two corner values --
as the STRINGS the paper prints. They are parsed by the small recursive-descent parser in section 2
and then re-derived from the definitions; nothing below is taken on trust from the block.
"""

import math
import re
import sys
from fractions import Fraction

# ---------------------------------------------------------------------------
# 0. THE OBJECT, AS PRINTED IN paper.tex
# ---------------------------------------------------------------------------
PAPER = {
    'n': 2,
    'D': 6,
    # MP_2, in the order printed in the paper's display (2.1).
    'nodes': [
        ('sqrt(2)/2',  '(sqrt(5)-1)/4'),
        ('sqrt(2)/2',  '-(sqrt(5)+1)/4'),
        ('0',          '(sqrt(5)+1)/4'),
        ('0',          '(1-sqrt(5))/4'),
        ('-sqrt(2)/2', '(sqrt(5)-1)/4'),
        ('-sqrt(2)/2', '-(sqrt(5)+1)/4'),
    ],
    # Table 1 of the paper: the six weights, in the same node order.
    'weights': [
        '(2+phi)/20', '(3-phi)/20', '(3-phi)/10',
        '(2+phi)/10', '(2+phi)/20', '(3-phi)/20',
    ],
    # Remark 2.3 of the paper gives the same four numbers in reciprocal form.
    'weights_alt_forms': ['1/(4*(3-phi))', '(5-sqrt(5))/40', '1/(2*(2+phi))', '1/(2*(3-phi))'],
    # The witness Q0 = (0, cos(2 pi/5)) and its antipode, and the kernel column at Q0.
    'witness': ('0', '(sqrt(5)-1)/4'),
    'antipode': ('0', '(1-sqrt(5))/4'),
    'antipode_index_m_k': (2, 2),
    'kernel_at_witness': ['4-2*phi', '-2', '2', '2', '4-2*phi', '-2'],
    'lambda_at_witness': '2-sqrt(5)/5',
    'lambda_at_witness_other_forms': ['2-1/sqrt(5)', '(10-sqrt(5))/5'],
    'lambda_at_antipode': '1',
    'difference': '1-sqrt(5)/5',
    # -MP_2, in the order obtained by negating `nodes`, with the values of Table 2.
    'lambda_on_minus_MP2': [
        '5/2-3*sqrt(5)/10', '2+3*sqrt(5)/5', '2+sqrt(5)/5',
        '2-sqrt(5)/5', '5/2-3*sqrt(5)/10', '2+3*sqrt(5)/5',
    ],
    # Remark 2.6 (corroboration only; attributed in the paper to the corner question).
    'corner_pp': 'sqrt(10)+6*sqrt(5)/5',
    'corner_pp_square': '12*sqrt(2)+86/5',
    'corner_mm': '3',
    # Section 3: the parity mechanism is checked for these even degrees.
    'even_n_range': (2, 40),
}

# ---------------------------------------------------------------------------
# 1. EXACT ARITHMETIC IN F = Q(sqrt2, sqrt5)
# ---------------------------------------------------------------------------
# An element is the 4-tuple (a, b, c, d) meaning a + b*sqrt2 + c*sqrt5 + d*sqrt10.
Z = (Fraction(0),) * 4


def E(a=0, b=0, c=0, d=0):
    return (Fraction(a), Fraction(b), Fraction(c), Fraction(d))


ONE = E(1)
SQ2 = E(0, 1)
SQ5 = E(0, 0, 1)
SQ10 = E(0, 0, 0, 1)


def add(u, v):
    return (u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3])


def sub(u, v):
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3])


def neg(u):
    return (-u[0], -u[1], -u[2], -u[3])


def mul(u, v):
    a1, b1, c1, d1 = u
    a2, b2, c2, d2 = v
    return (a1 * a2 + 2 * b1 * b2 + 5 * c1 * c2 + 10 * d1 * d2,
            a1 * b2 + b1 * a2 + 5 * (c1 * d2 + d1 * c2),
            a1 * c2 + c1 * a2 + 2 * (b1 * d2 + d1 * b2),
            a1 * d2 + d1 * a2 + b1 * c2 + c1 * b2)


def is_zero(u):
    return all(x == 0 for x in u)


def inv(u):
    """Multiplicative inverse, by solving the 4x4 rational system M(u) x = (1,0,0,0)."""
    if is_zero(u):
        raise ZeroDivisionError('inverse of 0 in Q(sqrt2,sqrt5)')
    cols = [mul(u, e) for e in (ONE, SQ2, SQ5, SQ10)]
    M = [[cols[j][i] for j in range(4)] + [Fraction(1) if i == 0 else Fraction(0)] for i in range(4)]
    for col in range(4):
        piv = next(r for r in range(col, 4) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(4):
            if r != col and M[r][col] != 0:
                fac = M[r][col]
                M[r] = [x - fac * y for x, y in zip(M[r], M[col])]
    return tuple(M[i][4] for i in range(4))


def div(u, v):
    return mul(u, inv(v))


def _root_bounds(N, digits):
    """Rational lo <= sqrt(N) <= hi with hi - lo = 10**-digits, from integer square roots only."""
    scale = 10 ** digits
    lo = Fraction(math.isqrt(N * scale * scale), scale)
    return lo, lo + Fraction(1, scale)


def sign(u):
    """-1, 0 or +1, decided exactly. 0 iff every coordinate is 0 (1,sqrt2,sqrt5,sqrt10 are
    Q-linearly independent). Otherwise a rational bracket that excludes 0 is exhibited."""
    if is_zero(u):
        return 0
    a, b, c, d = u
    digits = 12
    while True:
        lo2, hi2 = _root_bounds(2, digits)
        lo5, hi5 = _root_bounds(5, digits)
        lo10, hi10 = _root_bounds(10, digits)
        low = a + (b * lo2 if b >= 0 else b * hi2) + (c * lo5 if c >= 0 else c * hi5) \
                + (d * lo10 if d >= 0 else d * hi10)
        high = a + (b * hi2 if b >= 0 else b * lo2) + (c * hi5 if c >= 0 else c * lo5) \
                 + (d * hi10 if d >= 0 else d * lo10)
        if low > 0:
            return 1
        if high < 0:
            return -1
        digits *= 2
        if digits > 4096:                     # unreachable for a nonzero element; a guard, not a policy
            raise ArithmeticError('sign() failed to separate %r from 0' % (u,))


def absval(u):
    return u if sign(u) >= 0 else neg(u)


def approx(u, places=18):
    """A decimal rendering for the human reader, correct to `places` digits. ⛔ Formatted from the
    exact Fraction, never through `float`: a float carries ~17 significant digits and would silently
    corrupt the tail of an 18-place figure a referee is meant to compare against the paper."""
    a, b, c, d = u
    lo2, hi2 = _root_bounds(2, places + 6)
    lo5, hi5 = _root_bounds(5, places + 6)
    lo10, hi10 = _root_bounds(10, places + 6)
    v = a + b * (lo2 + hi2) / 2 + c * (lo5 + hi5) / 2 + d * (lo10 + hi10) / 2
    sgn = '-' if v < 0 else ''
    digits = str(round(abs(v) * 10 ** places)).rjust(places + 1, '0')
    return sgn + digits[:len(digits) - places] + '.' + digits[len(digits) - places:]


def show(u):
    """A canonical printable form, so a reader can match a PASS detail to the paper."""
    if is_zero(u):
        return '0'
    parts = []
    for coef, sym in zip(u, ('', '*sqrt(2)', '*sqrt(5)', '*sqrt(10)')):
        if coef:
            parts.append('%s%s' % (coef, sym))
    return ' + '.join(parts).replace('+ -', '- ')


# ---------------------------------------------------------------------------
# 2. THE PARSER FOR THE PAPER'S PRINTED STRINGS
# ---------------------------------------------------------------------------
PHI = div(add(ONE, SQ5), E(2))                     # (1+sqrt5)/2
_TOK = re.compile(r'\s*(\d+|sqrt|phi|[()+\-*/])')


def parse(text):
    toks, i = [], 0
    while i < len(text):
        m = _TOK.match(text, i)
        if not m:
            raise ValueError('cannot tokenise %r at %d' % (text, i))
        toks.append(m.group(1))
        i = m.end()
    pos = [0]

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else None

    def take(t=None):
        got = toks[pos[0]]
        if t is not None and got != t:
            raise ValueError('expected %r, got %r in %r' % (t, got, text))
        pos[0] += 1
        return got

    def atom():
        t = peek()
        if t == '(':
            take('(')
            v = expr()
            take(')')
            return v
        if t == 'sqrt':
            take('sqrt')
            take('(')
            v = expr()
            take(')')
            if v == E(2):
                return SQ2
            if v == E(5):
                return SQ5
            if v == E(10):
                return SQ10
            if v[1] == v[2] == v[3] == 0:          # a rational whose root is rational
                num, den = v[0].numerator, v[0].denominator
                rn, rd = math.isqrt(num), math.isqrt(den)
                if rn * rn == num and rd * rd == den:
                    return E(Fraction(rn, rd))
            raise ValueError('sqrt of %r is outside Q(sqrt2,sqrt5)' % (v,))
        if t == 'phi':
            take('phi')
            return PHI
        if t is not None and t.isdigit():
            return E(int(take()))
        raise ValueError('unexpected token %r in %r' % (t, text))

    def unary():
        if peek() == '-':
            take('-')
            return neg(unary())
        if peek() == '+':
            take('+')
        return atom()

    def term():
        v = unary()
        while peek() in ('*', '/'):
            op = take()
            w = unary()
            v = mul(v, w) if op == '*' else div(v, w)
        return v

    def expr():
        v = term()
        while peek() in ('+', '-'):
            op = take()
            w = term()
            v = add(v, w) if op == '+' else sub(v, w)
        return v

    v = expr()
    if pos[0] != len(toks):
        raise ValueError('trailing tokens in %r' % text)
    return v


def parse_pt(pair):
    return (parse(pair[0]), parse(pair[1]))


# ---------------------------------------------------------------------------
# 3. THE CHECK HARNESS
# ---------------------------------------------------------------------------
_STATE = {'n': 0, 'bad': 0}


def check(name, ok, detail=''):
    _STATE['n'] += 1
    if ok:
        print('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        _STATE['bad'] += 1
        print('FAIL %s%s' % (name, (' ' + detail) if detail else ''))


# ---------------------------------------------------------------------------
# 4. THE MORROW-PATTERSON MACHINERY AT n = 2
# ---------------------------------------------------------------------------
N = PAPER['n']
IJ = [(i, j) for i in range(N + 1) for j in range(N - i + 1)]        # i + j <= n
C_N = E(Fraction(8, (N + 2) * (N + 3)))                             # Thm cub_formula: C_2 = 2/5


def U(j, t):
    """Chebyshev polynomial of the second kind, by its recursion."""
    a, b = ONE, mul(E(2), t)
    if j == 0:
        return a
    for _ in range(j - 1):
        a, b = b, sub(mul(mul(E(2), t), b), a)
    return b


def kernel(P, Q):
    tot = Z
    for (i, j) in IJ:
        tot = add(tot, mul(mul(U(i, P[0]), U(j, P[1])), mul(U(i, Q[0]), U(j, Q[1]))))
    return tot


def weight_defining(P):
    """1/omega = sum_{i+j<=n} U_i(x)^2 U_j(y)^2  (the reproducing-kernel normalisation)."""
    tot = Z
    for (i, j) in IJ:
        tot = add(tot, mul(mul(U(i, P[0]), U(i, P[0])), mul(U(j, P[1]), U(j, P[1]))))
    return inv(tot)


def weight_closed(P):
    """omega = C_n (1-x^2)(1-y^2)  -- the paper's own cubature formula."""
    return mul(C_N, mul(sub(ONE, mul(P[0], P[0])), sub(ONE, mul(P[1], P[1]))))


def lam(Q, nodes, weights):
    tot = Z
    for P, w in zip(nodes, weights):
        tot = add(tot, mul(w, absval(kernel(P, Q))))
    return tot


# ---------------------------------------------------------------------------
# 5. THE CHECKS
# ---------------------------------------------------------------------------
def main():
    print('re-derivation of "The Morrow-Patterson Lebesgue Function Is Not Centrally Symmetric"')
    print('exact arithmetic in Q(sqrt2,sqrt5); integer arithmetic for the parity mechanism')
    print('python %d.%d, standard library only' % sys.version_info[:2])
    print('')

    # --- 5.1 the cosine values the node set is built from -----------------
    # cos(k pi/5), k=1..4, are pinned by T_5(t) = (-1)^k together with the ordering: T_5 - 1 has
    # roots cos(2k pi/5) only, T_5 + 1 has roots cos((2k+1) pi/5) only, in (-1,1).
    def T5(t):
        t2 = mul(t, t)
        t3 = mul(t2, t)
        t5 = mul(t2, t3)
        return add(sub(mul(E(16), t5), mul(E(20), t3)), mul(E(5), t))

    c = {1: parse('(sqrt(5)+1)/4'), 2: parse('(sqrt(5)-1)/4'),
         3: parse('(1-sqrt(5))/4'), 4: parse('-(sqrt(5)+1)/4')}
    for k in (1, 2, 3, 4):
        want = ONE if k % 2 == 0 else neg(ONE)
        check('cos-%dpi-over-5-satisfies-T5-eq-%d' % (k, 1 if k % 2 == 0 else -1),
              T5(c[k]) == want, '[T_5(%s) = %s]' % (show(c[k]), show(T5(c[k]))))
    ordered = all(sign(sub(c[k], c[k + 1])) == 1 for k in (1, 2, 3)) \
        and sign(c[2]) == 1 and sign(c[3]) == -1 \
        and all(sign(sub(ONE, absval(c[k]))) == 1 for k in (1, 2, 3, 4))
    check('fifth-cosines-strictly-decreasing-and-interior', ordered,
          '[%s > %s > 0 > %s > %s, all in (-1,1)]'
          % (approx(c[1], 12), approx(c[2], 12), approx(c[3], 12), approx(c[4], 12)))

    a = parse('sqrt(2)/2')
    check('cos-pi-over-4-satisfies-T2-eq-0-and-is-positive',
          sub(mul(E(2), mul(a, a)), ONE) == Z and sign(a) == 1,
          '[2*(%s)^2 - 1 = 0, value %s]' % (show(a), approx(a, 12)))

    # --- 5.2 the node set ------------------------------------------------
    # MP_n = {(cos(m pi/(n+2)), y) : m = 1..n+1}, y = cos(2k pi/(n+3)) for m odd,
    # y = cos((2k-1) pi/(n+3)) for m even, k = 1..n/2+1.
    xs = {1: a, 2: Z, 3: neg(a)}
    built = []
    for m in (1, 2, 3):
        for k in (1, 2):
            jj = 2 * k if m % 2 == 1 else 2 * k - 1
            built.append((xs[m], c[jj], m, k))
    derived = [(p[0], p[1]) for p in built]
    printed = [parse_pt(p) for p in PAPER['nodes']]

    D = (N + 1) * (N + 2) // 2
    check('cardinality-of-MP2-is-6-equals-dim-P2-squared',
          len(derived) == D == PAPER['D'] == len(set(derived)),
          '[|MP_2| = %d = (n+1)(n+2)/2 = dim P_2^2, all six distinct]' % len(derived))
    check('MP2-derived-from-the-definition-equals-the-set-printed-in-the-paper',
          derived == printed,
          '[node 1 = (%s, %s), node 4 = (%s, %s)]'
          % (show(derived[0][0]), show(derived[0][1]), show(derived[3][0]), show(derived[3][1])))

    Q0 = parse_pt(PAPER['witness'])
    A0 = parse_pt(PAPER['antipode'])
    check('witness-Q0-is-not-a-node-of-MP2', Q0 not in derived,
          '[Q0 = (0, %s) = (0, %s); the m=2 column carries y in {%s, %s} only]'
          % (show(Q0[1]), approx(Q0[1], 12), show(c[1]), show(c[3])))
    check('antipode-of-Q0-is-the-node-with-index-m2-k2',
          A0 == (neg(Q0[0]), neg(Q0[1])) and A0 in derived
          and (built[[i for i, p in enumerate(built) if (p[0], p[1]) == A0][0]][2],
               built[[i for i, p in enumerate(built) if (p[0], p[1]) == A0][0]][3])
          == PAPER['antipode_index_m_k'],
          '[-Q0 = (0, %s) is node (m,k) = %s]' % (show(A0[1]), str(PAPER['antipode_index_m_k'])))

    # --- 5.3 the weights, two independent ways ---------------------------
    w_def = [weight_defining(P) for P in derived]
    w_clo = [weight_closed(P) for P in derived]
    w_pap = [parse(s) for s in PAPER['weights']]
    check('weights-defining-sum-and-the-papers-cubature-closed-form-agree-at-all-six-nodes',
          w_def == w_clo,
          '[%s]' % ', '.join(approx(x, 7) for x in w_def))
    check('weights-match-the-table-printed-in-the-paper', w_def == w_pap,
          '[omega_2 = %s = %s]' % (show(w_def[1]), approx(w_def[1], 12)))
    check('the-reciprocal-forms-printed-in-the-papers-remark-are-the-same-four-numbers',
          [parse(s) for s in PAPER['weights_alt_forms']] == w_def[:4],
          '[1/(4(3-phi)), (5-sqrt5)/40, 1/(2(2+phi)), 1/(2(3-phi))]')
    tot = Z
    for w in w_def:
        tot = add(tot, w)
    check('weights-sum-to-one-exactly', tot == ONE, '[sum omega_P = %s]' % show(tot))

    # --- 5.4 unisolvence and the Lagrange property -----------------------
    delta_ok = True
    for i, Pi in enumerate(derived):
        for j, Pj in enumerate(derived):
            want = ONE if i == j else Z
            if mul(w_def[i], kernel(Pi, Pj)) != want:
                delta_ok = False
    check('lagrange-delta-property-omega_i-K(P_i,P_j)-eq-delta_ij-on-all-36-pairs', delta_ok,
          '[36 exact equalities; MP_2 is unisolvent for P_2^2]')
    ones = [lam(P, derived, w_def) for P in derived]
    check('lambda2-equals-1-exactly-at-all-six-nodes-of-MP2', all(v == ONE for v in ones),
          '[max deviation exactly 0]')

    # --- 5.5 the kernel column at the witness ----------------------------
    Kcol = [kernel(P, Q0) for P in derived]
    Kpap = [parse(s) for s in PAPER['kernel_at_witness']]
    check('kernel-column-at-Q0-matches-the-six-values-printed-in-the-paper', Kcol == Kpap,
          '[%s]' % ', '.join(show(x) for x in Kcol))
    collapse = all(
        Kcol[i] == sub(add(ONE, mul(sub(PHI, ONE), sub(U(1, P[1]), U(2, P[1])))), U(2, P[0]))
        for i, P in enumerate(derived))
    check('kernel-collapse-identity-at-Q0-U1(0)=0-U2(0)=-1', collapse,
          '[K(P,Q0) = 1 + (phi-1)[U_1(y_P) - U_2(y_P)] - U_2(x_P) at all six P]')
    rep = Z
    for w, K in zip(w_def, Kcol):
        rep = add(rep, mul(w, K))
    check('reproducing-identity-sum-omega_P-K(P,Q0)-eq-1', rep == ONE, '[sum = %s]' % show(rep))
    negs = [i for i, K in enumerate(Kcol) if sign(K) < 0]
    check('exactly-two-of-the-six-kernel-entries-at-Q0-are-negative',
          negs == [1, 5] and all(Kcol[i] == E(-2) for i in negs),
          '[negative at P_2 = (sqrt(2)/2, -cos(pi/5)) and P_6 = (-sqrt(2)/2, -cos(pi/5)), K = -2 each]')

    # --- 5.6 THE REFUTATION ----------------------------------------------
    lam_Q0 = lam(Q0, derived, w_def)
    lam_A0 = lam(A0, derived, w_def)
    check('lambda2-at-Q0-equals-2-minus-sqrt5-over-5', lam_Q0 == parse(PAPER['lambda_at_witness']),
          '[lambda_2(0, cos(2pi/5)) = %s = %s]' % (show(lam_Q0), approx(lam_Q0, 18)))
    check('the-three-printed-forms-of-that-value-are-the-same-field-element',
          all(parse(s) == lam_Q0 for s in PAPER['lambda_at_witness_other_forms']),
          '[2 - sqrt(5)/5 = 2 - 1/sqrt(5) = (10-sqrt(5))/5]')
    check('lambda2-at-the-antipode-equals-1-exactly',
          lam_A0 == parse(PAPER['lambda_at_antipode']) == ONE,
          '[-Q0 is an interpolation node, so lambda_2(-Q0) = 1 with no arithmetic]')
    diff = sub(lam_Q0, lam_A0)
    check('central-symmetry-FAILS-the-difference-is-1-minus-sqrt5-over-5-and-is-nonzero',
          diff == parse(PAPER['difference']) and sign(diff) == 1,
          '[lambda_2(Q0) - lambda_2(-Q0) = %s = %s > 0]' % (show(diff), approx(diff, 18)))
    check('the-hand-argument-1-plus-8-omega_2-gives-the-same-value',
          add(ONE, mul(E(8), w_def[1])) == lam_Q0,
          '[1 + 2*(negative mass) = 1 + 8*omega_2 = %s]' % show(lam_Q0))

    # --- 5.7 all six antipodal node pairs at n = 2 -----------------------
    minus = [(neg(P[0]), neg(P[1])) for P in derived]
    lam_minus = [lam(Q, derived, w_def) for Q in minus]
    check('minus-MP2-values-match-the-six-printed-in-the-paper',
          lam_minus == [parse(s) for s in PAPER['lambda_on_minus_MP2']],
          '[%s]' % ', '.join(show(v) for v in lam_minus))
    check('none-of-the-six-values-on-minus-MP2-equals-1',
          all(v != ONE and sign(sub(v, ONE)) == 1 for v in lam_minus),
          '[all six exceed 1, so symmetry fails at all 6 of the 6 antipodal node pairs]')
    check('Q0-carries-the-smallest-of-the-six-violations',
          all(sign(sub(v, lam_Q0)) >= 0 for v in lam_minus),
          '[min over -MP_2 is %s at Q0, the most conservative witness]' % approx(lam_Q0, 12))
    Kmat = [[kernel(P, Q) for Q in minus] for P in derived]
    nneg = sum(1 for row in Kmat for x in row if sign(x) < 0)
    check('the-6x6-kernel-matrix-MP2-by-minus-MP2-is-not-entrywise-nonnegative',
          nneg > 0,
          '[%d of its 36 entries are negative; one negative entry already kills the criterion]' % nneg)

    # --- 5.8 the parity mechanism, integers only, every even n in 2..40 --
    lo, hi = PAPER['even_n_range']
    evens = list(range(lo, hi + 1, 2))

    def mp_index(n):
        return set((m, jj) for m in range(1, n + 2) for jj in range(1, n + 3)
                   if (jj % 2 == 0) == (m % 2 == 1))

    def mp_from_paper_parametrisation(n):
        s = set()
        for m in range(1, n + 2):
            for k in range(1, n // 2 + 2):
                s.add((m, 2 * k if m % 2 == 1 else 2 * k - 1))
        return s

    def grid(n):
        return set((m, jj) for m in range(1, n + 2) for jj in range(1, n + 3))

    check('the-papers-y-parametrisation-equals-the-parity-rule-jj-even-iff-m-odd',
          all(mp_from_paper_parametrisation(n) == mp_index(n) for n in evens),
          '[even n = %d..%d]' % (lo, hi))
    check('cardinality-of-MP_n-is-(n+1)(n+2)/2-eq-dim-P_n_squared',
          all(len(mp_index(n)) == (n + 1) * (n + 2) // 2 for n in evens),
          '[even n = %d..%d]' % (lo, hi))
    check('sigma_x-(m,jj)->(n+2-m,jj)-preserves-MP_n-because-p=n+2-is-even',
          all(set((n + 2 - m, jj) for (m, jj) in mp_index(n)) == mp_index(n) for n in evens),
          '[this is the paper-proved reflection lemma, even n = %d..%d]' % (lo, hi))
    check('sigma_y-(m,jj)->(n+3-jj,)-flips-the-jj-parity-because-q=n+3-is-odd',
          all(all((jj % 2) != ((n + 3 - jj) % 2) for (m, jj) in mp_index(n)) for n in evens),
          '[even n = %d..%d]' % (lo, hi))
    check('minus-MP_n-is-DISJOINT-from-MP_n',
          all(not (set((n + 2 - m, n + 3 - jj) for (m, jj) in mp_index(n)) & mp_index(n))
              for n in evens),
          '[even n = %d..%d; so every point of -MP_n is a non-node]' % (lo, hi))
    check('MP_n-together-with-minus-MP_n-PARTITION-the-full-tensor-grid',
          all(set((n + 2 - m, n + 3 - jj) for (m, jj) in mp_index(n)) | mp_index(n) == grid(n)
              and len(grid(n)) == (n + 1) * (n + 2) for n in evens),
          '[even n = %d..%d]' % (lo, hi))

    # --- 5.9 corroboration at the corners (attributed, not claimed here) -
    cpp = lam((ONE, ONE), derived, w_def)
    cmm = lam((neg(ONE), neg(ONE)), derived, w_def)
    check('corner-value-lambda2(1,1)-equals-sqrt10-plus-6sqrt5-over-5',
          cpp == parse(PAPER['corner_pp']), '[%s = %s]' % (show(cpp), approx(cpp, 18)))
    check('that-corner-value-squared-equals-12sqrt2-plus-86-over-5',
          mul(cpp, cpp) == parse(PAPER['corner_pp_square']), '[%s]' % show(mul(cpp, cpp)))
    check('corner-value-lambda2(-1,-1)-equals-3-exactly',
          cmm == parse(PAPER['corner_mm']) == E(3), '[%s]' % show(cmm))
    check('the-two-corner-values-differ-so-the-corner-pair-is-a-second-witness',
          sign(sub(cpp, cmm)) == 1,
          '[difference %s, i.e. 12sqrt2 + 41/5 > 0 after squaring]' % approx(sub(cpp, cmm), 12))

    # --- scope ------------------------------------------------------------
    print('')
    print('NOT RE-RUN: (1) the FAILURE of central symmetry at n >= 4. Only the integer parity '
          'mechanism of section 5.8 is re-derived for even n up to 40; the numerical minima of '
          'lambda_n over -MP_n reported for n = 4..30 need multiprecision transcendental arithmetic '
          'and are outside the standard library, so they are NOT recomputed here.')
    print('NOT RE-RUN: (2) the uniform candidate family and the negativity of its one alternating '
          'sum, which the paper states is unproved and only computed; nothing here bears on it.')
    print('NOT RE-RUN: (3) the 101x101 grid census at n = 30 and the diagnosis of the authors\' '
          'plotting script; both are corroboration reported elsewhere and neither is a claim of '
          'this paper.')
    print('NOT RE-RUN: (4) the corner values of section 5.9 are checked but are NOT claimed as a '
          'result of this paper -- they belong to the separate corner-attainment question and are '
          'carried only as an independent second witness to the asymmetry.')
    print('')
    if _STATE['bad']:
        print('VERDICT: %d of %d CHECKS FAILED' % (_STATE['bad'], _STATE['n']))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _STATE['n'])
    return 0


if __name__ == '__main__':
    sys.exit(main())

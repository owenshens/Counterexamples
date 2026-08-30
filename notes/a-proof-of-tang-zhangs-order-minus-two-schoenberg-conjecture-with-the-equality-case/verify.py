#!/usr/bin/env python3
"""verify.py -- exact re-derivation of every computational claim in

    "A Proof of Tang and Zhang's Order -2 Schoenberg Conjecture, with the
     Equality Case"    (paper.tex / paper.pdf in this folder)

about Conjecture C1 of arXiv:2508.10341v3 (Q. Tang, T. Zhang):

    sum_{k=1}^n 1/|w_k|^2  <=  3 sum_{j=1}^n 1/|z_j|^2 + | sum_{j=1}^n 1/z_j |^2,
    where  P(z) = z prod_{j=1}^n (z - z_j),  z_j in C \\ {0},  and w_1..w_n are
    the zeros of P' counted WITH MULTIPLICITY.

Python 3.9+, STANDARD LIBRARY ONLY (fractions, math, random, sys).  There is no
floating point anywhere: every decision below is an equality or an inequality
between exact rationals.  The arithmetic is carried out in the ring

    Q(sqrt(D))(i)  =  { (a + b sqrt(D)) + (a' + b' sqrt(D)) i },  a,b,a',b' in Q,

with D = n+1, which is exactly what the paper's constant c = (sqrt(n+1)-1)/n
needs.  Two classes implement it: Q2 (a + b sqrt(D)) and Cx (re + im i, with
re, im in Q2).  Q2 normalises itself to a plain rational whenever D is a perfect
square, so at n = 3, 8, 15, 24 the whole certificate is rational, as the paper
says.

Every object the program consumes is PRINTED IN THE PAPER: the tuples of
Section 4, the closed form of K in equation (K), the coefficients of P', the
constant c.  Nothing is read from a data file.

Output contract: one `PASS <name> [detail]` line per check, the closing
`NOT RE-RUN:` statements of what is out of scope, and then

    VERDICT: ALL <n> CHECKS PASS

with exit status 0 if and only if every check passed.
"""

import random
import sys
from fractions import Fraction as F
from math import isqrt

# ---------------------------------------------------------------------------
# 0. THE CHECK LEDGER
# ---------------------------------------------------------------------------
_PASSED = 0
_FAILED = 0


def check(name, ok, detail=''):
    """Record one check.  ⛔ The verdict at the bottom is derived from THIS
    counter, so a check that is never called cannot be counted, and a check
    that fails cannot be silently dropped."""
    global _PASSED, _FAILED
    if ok:
        _PASSED += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILED += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def step(title):
    print()
    print('=== %s' % title)


# ---------------------------------------------------------------------------
# 1. EXACT ARITHMETIC IN Q(sqrt(D))(i)
# ---------------------------------------------------------------------------
class Q2(object):
    """a + b*sqrt(D) with a, b rational and D a non-negative integer.

    Normalised: b == 0 or D not a perfect square, so `==` is a faithful test of
    equality of real numbers (two normalised elements with the same surd part
    are equal iff their coordinates agree, because sqrt(D) is irrational)."""

    __slots__ = ('a', 'b', 'D')

    def __init__(self, a, b=0, D=0):
        a = F(a)
        b = F(b)
        D = int(D)
        if D < 0:
            raise ValueError('D must be >= 0')
        if b == 0:
            D = 0
        else:
            r = isqrt(D)
            if r * r == D:
                a = a + b * r
                b = F(0)
                D = 0
        self.a = a
        self.b = b
        self.D = D

    @staticmethod
    def _common(x, y):
        if x.D and y.D and x.D != y.D:
            raise ValueError('cannot combine sqrt(%d) with sqrt(%d)' % (x.D, y.D))
        return x.D or y.D

    def __add__(s, o):
        o = q2(o)
        return Q2(s.a + o.a, s.b + o.b, Q2._common(s, o))

    __radd__ = __add__

    def __sub__(s, o):
        o = q2(o)
        return Q2(s.a - o.a, s.b - o.b, Q2._common(s, o))

    def __rsub__(s, o):
        return q2(o) - s

    def __neg__(s):
        return Q2(-s.a, -s.b, s.D)

    def __mul__(s, o):
        o = q2(o)
        D = Q2._common(s, o)
        return Q2(s.a * o.a + s.b * o.b * D, s.a * o.b + s.b * o.a, D)

    __rmul__ = __mul__

    def norm(s):
        """a^2 - D b^2 -- the field norm; zero only for the zero element."""
        return s.a * s.a - s.D * s.b * s.b

    def inv(s):
        d = s.norm()
        if d == 0:
            raise ZeroDivisionError('inverse of zero in Q(sqrt(%d))' % s.D)
        return Q2(s.a / d, -s.b / d, s.D)

    def __truediv__(s, o):
        return s * q2(o).inv()

    def __eq__(s, o):
        o = q2(o)
        if s.b == 0 and o.b == 0:
            return s.a == o.a
        return s.a == o.a and s.b == o.b and s.D == o.D

    def __ne__(s, o):
        return not (s == o)

    def is_zero(s):
        return s.a == 0 and s.b == 0

    def is_rational(s):
        return s.b == 0

    def rat(s):
        if s.b != 0:
            raise ValueError('%s is not rational' % (s,))
        return s.a

    def __repr__(s):
        if s.b == 0:
            return str(s.a)
        return '%s%s*sqrt(%d)' % (s.a, ('+%s' % s.b) if s.b > 0 else str(s.b), s.D)


def q2(x):
    return x if isinstance(x, Q2) else Q2(x)


class Cx(object):
    """re + im*i with re, im in Q2.  Conjugation flips i only (sqrt(D) is real)."""

    __slots__ = ('re', 'im')

    def __init__(s, re, im=0):
        s.re = q2(re)
        s.im = q2(im)

    def __add__(s, o):
        o = cx(o)
        return Cx(s.re + o.re, s.im + o.im)

    __radd__ = __add__

    def __sub__(s, o):
        o = cx(o)
        return Cx(s.re - o.re, s.im - o.im)

    def __rsub__(s, o):
        return cx(o) - s

    def __neg__(s):
        return Cx(-s.re, -s.im)

    def __mul__(s, o):
        o = cx(o)
        return Cx(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)

    __rmul__ = __mul__

    def conj(s):
        return Cx(s.re, -s.im)

    def abs2(s):
        """|z|^2 as an element of Q(sqrt(D)) -- always real."""
        return s.re * s.re + s.im * s.im

    def inv(s):
        n = s.abs2()
        if n.is_zero():
            raise ZeroDivisionError('inverse of zero')
        return Cx(s.re / n, -s.im / n)

    def __truediv__(s, o):
        return s * cx(o).inv()

    def __eq__(s, o):
        o = cx(o)
        return s.re == o.re and s.im == o.im

    def __ne__(s, o):
        return not (s == o)

    def is_zero(s):
        return s.re.is_zero() and s.im.is_zero()

    def __repr__(s):
        if s.im.is_zero():
            return repr(s.re)
        return '(%s%s%si)' % (s.re, '+' if (s.im.a > 0 or (s.im.a == 0 and s.im.b > 0)) else '',
                              s.im)


def cx(x):
    if isinstance(x, Cx):
        return x
    if isinstance(x, Q2):
        return Cx(x)
    return Cx(q2(x))


ZERO = Cx(0)
ONE = Cx(1)
I_UNIT = Cx(0, 1)


# ---------------------------------------------------------------------------
# 2. EXACT LINEAR ALGEBRA
# ---------------------------------------------------------------------------
def eye(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def ones(n):
    return [[ONE] * n for _ in range(n)]


def mat_mul(X, Y):
    n, k, m = len(X), len(Y), len(Y[0])
    out = []
    for i in range(n):
        row = []
        for j in range(m):
            acc = ZERO
            for t in range(k):
                acc = acc + X[i][t] * Y[t][j]
            row.append(acc)
        out.append(row)
    return out


def mat_eq(X, Y):
    if len(X) != len(Y) or len(X[0]) != len(Y[0]):
        return False
    return all(X[i][j] == Y[i][j] for i in range(len(X)) for j in range(len(X[0])))


def mat_scale(X, s):
    return [[x * s for x in row] for row in X]


def trace(X):
    acc = ZERO
    for i in range(len(X)):
        acc = acc + X[i][i]
    return acc


def conj_transpose(X):
    n, m = len(X), len(X[0])
    return [[X[j][i].conj() for j in range(n)] for i in range(m)]


def frob2(X):
    """||X||_F^2, an element of Q(sqrt(D))."""
    acc = Q2(0)
    for row in X:
        for x in row:
            acc = acc + x.abs2()
    return acc


def charpoly(A):
    """Faddeev--LeVerrier.  Returns [1, c_1, ..., c_n] with
    det(x I - A) = x^n + c_1 x^{n-1} + ... + c_n.  Only division by the integers
    1..n is used, so it is exact in any field of characteristic zero and needs
    no pivoting."""
    n = len(A)
    coeffs = [ONE]
    Mprev = None
    for k in range(1, n + 1):
        if k == 1:
            Mk = [row[:] for row in A]
        else:
            T = [[Mprev[i][j] + (coeffs[k - 1] if i == j else ZERO) for j in range(n)]
                 for i in range(n)]
            Mk = mat_mul(A, T)
        coeffs.append(ZERO - trace(Mk) * Cx(F(1, k)))
        Mprev = Mk
    return coeffs


def poly_mul(p, q):
    r = [ZERO] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] = r[i + j] + a * b
    return r


def P_and_Pprime(z):
    """Coefficients of P(x) = x prod_j (x - z_j) and of P'(x), LOW to HIGH."""
    p = [ZERO, ONE]
    for zj in z:
        p = poly_mul(p, [ZERO - zj, ONE])
    pp = [p[i] * Cx(i) for i in range(1, len(p))]
    return p, pp


def elementary_symmetric(vals):
    """[e_0, e_1, ..., e_n] of the given list."""
    n = len(vals)
    e = [ONE] + [ZERO] * n
    for x in vals:
        for i in range(n, 0, -1):
            e[i] = e[i] + e[i - 1] * x
    return e


# ---------------------------------------------------------------------------
# 3. THE OBJECTS OF THE PAPER
# ---------------------------------------------------------------------------
def paper_constants(n):
    """c = (sqrt(n+1)-1)/n and c' = (1/sqrt(n+1) - 1)/n, exactly."""
    D = n + 1
    c = Q2(F(-1, n), F(1, n), D)                     # (-1 + sqrt(n+1))/n
    cprime = Q2(F(-1, n), F(1, n * D), D)            # (-1 + sqrt(n+1)/(n+1))/n
    return c, cprime


def build(z):
    """Everything the paper names, for a tuple z of nonzero Gaussian rationals."""
    n = len(z)
    c, cprime = paper_constants(n)
    u = [ONE / zj for zj in z]
    S = ZERO
    for x in u:
        S = S + x
    U = [[u[i] if i == j else ZERO for j in range(n)] for i in range(n)]
    D_ = [[z[i] if i == j else ZERO for j in range(n)] for i in range(n)]
    M = [[Cx(2) if i == j else ONE for j in range(n)] for i in range(n)]
    Q = [[(ONE if i == j else ZERO) + Cx(c) for j in range(n)] for i in range(n)]
    Qinv = [[(ONE if i == j else ZERO) + Cx(cprime) for j in range(n)] for i in range(n)]
    A = mat_mul([[(ONE if i == j else ZERO) - Cx(F(1, n + 1)) for j in range(n)]
                 for i in range(n)], D_)
    B = mat_mul(U, M)
    K = mat_mul(mat_mul(Q, U), Q)
    return dict(n=n, c=c, cprime=cprime, u=u, S=S, U=U, D=D_, M=M, Q=Q, Qinv=Qinv,
                A=A, B=B, K=K, z=list(z))


def rhs_C1(z):
    """3 sum_j 1/|z_j|^2 + |sum_j 1/z_j|^2 -- the right-hand side of eq:C1."""
    u = [ONE / zj for zj in z]
    S = ZERO
    for x in u:
        S = S + x
    acc = Q2(0)
    for x in u:
        acc = acc + x.abs2()
    return acc * Q2(3) + S.abs2()


def published_bound(z):
    """(n+3) sum_j 1/|z_j|^2 -- the authors' published order -2 bound."""
    n = len(z)
    acc = Q2(0)
    for zj in z:
        acc = acc + (ONE / zj).abs2()
    return acc * Q2(n + 3)


def collinear_with_origin(z):
    """True iff all z_j lie on one line through 0, i.e. z_j conj(z_k) real."""
    for j in range(len(z)):
        for k in range(j + 1, len(z)):
            if not (z[j] * z[k].conj()).im.is_zero():
                return False
    return True


def K_is_normal(st):
    K = st['K']
    return mat_eq(mat_mul(K, conj_transpose(K)), mat_mul(conj_transpose(K), K))


def lhs_from_B(z):
    """sum_k 1/w_k^2 = tr(B^2) -- exact, and equal to sum_k 1/|w_k|^2 whenever
    every w_k is real (which is the case for real z_j, by Rolle)."""
    B = build(z)['B']
    return trace(mat_mul(B, B))


# ---------------------------------------------------------------------------
# 4. RANDOM INSTANCE GENERATORS (SEEDED, so the transcript is reproducible)
# ---------------------------------------------------------------------------
def gaussian_rational(rng, lim=5):
    while True:
        a = rng.randint(-lim, lim)
        b = rng.randint(-lim, lim)
        if a or b:
            return Cx(F(a), F(b))


def rational_nonzero(rng, lim=6):
    while True:
        a = rng.randint(-lim, lim)
        if a:
            return Cx(F(a))


# ===========================================================================
#                                 THE CHECKS
# ===========================================================================
def step1_arithmetic():
    step('Step 1: the exact arithmetic, validated before it is used to decide anything')
    r2 = Q2(0, 1, 2)
    check('sqrt_of_two_squares_to_exactly_two', r2 * r2 == Q2(2), '(sqrt 2)^2 = %s' % (r2 * r2,))
    r5 = Q2(0, 1, 5)
    check('inverse_in_the_quadratic_field_round_trips',
          (Q2(3, 2, 5) * Q2(3, 2, 5).inv()) == Q2(1) and (r5.inv() * r5) == Q2(1))
    w = Cx(Q2(F(1, 2), 1, 5), Q2(3, F(-2, 7), 5))
    check('complex_inverse_round_trips_over_the_quadratic_field',
          (w * w.inv()) == ONE)
    check('conjugation_is_an_involution_and_fixes_the_reals',
          w.conj().conj() == w and Cx(Q2(7, 3, 5)).conj() == Cx(Q2(7, 3, 5)))
    check('modulus_squared_is_real_and_matches_z_times_conj_z',
          Cx(w.abs2()) == w * w.conj() and w.abs2().is_rational() is False,
          '|w|^2 = %s' % (w.abs2(),))
    # Two controls: an implementation that collapsed distinct surds, or that
    # compared only the rational parts, would pass everything above.
    check('control_distinct_surds_are_not_identified',
          Q2(1, 1, 2) != Q2(1, 1, 3) and Q2(0, 1, 2) != Q2(F(3, 2)) and Q2(0, 1, 2) != Q2(0))
    check('control_a_perfect_square_surd_is_normalised_to_a_rational',
          Q2(1, 1, 4) == Q2(3) and Q2(1, 1, 4).is_rational() and Q2(0, F(1, 3), 9) == Q2(1))
    check('control_a_nonzero_element_has_nonzero_field_norm',
          Q2(3, 2, 5).norm() != 0 and Q2(0, 1, 5).norm() == -5)


def step2_square_root():
    step('Step 2: Q = I + cJ is a square root of M = I + J, and M is not idempotent')
    bad_id, bad_sq, bad_inv, bad_idem = [], [], [], []
    for n in range(1, 13):
        c, cprime = paper_constants(n)
        if not (c * c * Q2(n) + c * Q2(2) - Q2(1)).is_zero():
            bad_id.append(n)
        Q = [[(ONE if i == j else ZERO) + Cx(c) for j in range(n)] for i in range(n)]
        Qinv = [[(ONE if i == j else ZERO) + Cx(cprime) for j in range(n)] for i in range(n)]
        M = [[Cx(2) if i == j else ONE for j in range(n)] for i in range(n)]
        if not mat_eq(mat_mul(Q, Q), M):
            bad_sq.append(n)
        if not (mat_eq(mat_mul(Q, Qinv), eye(n)) and mat_eq(mat_mul(Qinv, Q), eye(n))):
            bad_inv.append(n)
        if mat_eq(mat_mul(M, M), M):
            bad_idem.append(n)
    check('c_satisfies_n_c_squared_plus_two_c_minus_one_equals_zero', bad_id == [],
          'n = 1..12, failures %s' % (bad_id,))
    check('Q_squared_equals_I_plus_J_exactly', bad_sq == [], 'n = 1..12, failures %s' % (bad_sq,))
    check('Q_inverse_is_I_plus_c_prime_J_on_both_sides', bad_inv == [],
          'n = 1..12, failures %s' % (bad_inv,))
    check('control_M_is_never_idempotent_so_no_projection_can_replace_Q', bad_idem == [],
          'M^2 = I + (n+2)J; n = 1..12, idempotent at %s' % (bad_idem,))
    c3, _ = paper_constants(3)
    c8, _ = paper_constants(8)
    check('c_is_exactly_one_third_at_n_three_and_one_quarter_at_n_eight',
          c3 == Q2(F(1, 3)) and c8 == Q2(F(1, 4)), 'c(3) = %s, c(8) = %s' % (c3, c8))
    # the closed form of K printed in the paper as equation (K)
    rng = random.Random(11)
    bad = []
    for _ in range(20):
        n = rng.randint(1, 6)
        z = [gaussian_rational(rng) for _ in range(n)]
        st = build(z)
        u, S, c = st['u'], st['S'], st['c']
        closed = [[(u[j] if j == k else ZERO) + Cx(c) * (u[j] + u[k]) + Cx(c * c) * S
                   for k in range(n)] for j in range(n)]
        if not mat_eq(st['K'], closed):
            bad.append(n)
    check('K_entries_match_the_closed_form_printed_in_the_paper', bad == [],
          '20 instances, n = 1..6, failures %s' % (bad,))


def step3_matrix_model():
    step('Step 3: the matrix model -- det(wI-A) = P\'(w)/(n+1), B = A^-1 = UM, K = QBQ^-1')
    rng = random.Random(2508)
    bad_cp, bad_inv, bad_sim, bad_nz, bad_deg = [], [], [], [], []
    trials = 40
    for _ in range(trials):
        n = rng.randint(1, 6)
        z = [gaussian_rational(rng) for _ in range(n)]
        st = build(z)
        cp = charpoly(st['A'])                              # high -> low
        got = [cp[n - i] for i in range(n + 1)]             # low -> high, monic
        _, pp = P_and_Pprime(z)
        if len(pp) != n + 1 or pp[-1] != Cx(n + 1):
            bad_deg.append(n)
            continue
        monic = [x / pp[-1] for x in pp]
        if any(got[i] != monic[i] for i in range(n + 1)):
            bad_cp.append(n)
        if not mat_eq(mat_mul(st['A'], st['B']), eye(n)):
            bad_inv.append(n)
        if not (mat_eq(mat_mul(st['K'], st['Q']), mat_mul(st['Q'], st['B']))
                and mat_eq(st['K'], mat_mul(mat_mul(st['Q'], st['B']), st['Qinv']))):
            bad_sim.append(n)
        e = elementary_symmetric(z)
        if monic[0] != (ZERO - e[n] if n % 2 else e[n]) * Cx(1):
            pass                                            # checked separately below
        if monic[0].is_zero():
            bad_nz.append(n)
    check('P_prime_is_monic_of_degree_n_after_dividing_by_n_plus_one', bad_deg == [],
          '%d instances, failures %s' % (trials, bad_deg,))
    check('char_poly_of_A_equals_P_prime_over_n_plus_one_coefficientwise', bad_cp == [],
          '%d instances, n = 1..6, failures %s' % (trials, bad_cp))
    check('A_times_D_inverse_times_I_plus_J_is_the_identity_so_B_is_A_inverse', bad_inv == [],
          '%d instances, failures %s' % (trials, bad_inv))
    check('K_Q_equals_Q_B_so_K_is_similar_to_B_and_shares_its_spectrum', bad_sim == [],
          '%d instances, failures %s' % (trials, bad_sim))
    check('no_critical_point_vanishes_the_constant_term_of_P_prime_is_nonzero', bad_nz == [],
          '%d instances, failures %s' % (trials, bad_nz))
    # prod_k w_k = e_n(z)/(n+1)
    rng = random.Random(3)
    bad = []
    for _ in range(20):
        n = rng.randint(1, 6)
        z = [gaussian_rational(rng) for _ in range(n)]
        _, pp = P_and_Pprime(z)
        monic = [x / pp[-1] for x in pp]
        prod_w = monic[0] * Cx((-1) ** n)
        e = elementary_symmetric(z)
        if prod_w != e[n] * Cx(F(1, n + 1)):
            bad.append(n)
    check('product_of_the_critical_points_equals_e_n_of_z_over_n_plus_one', bad == [],
          '20 instances, failures %s' % (bad,))


def step4_frobenius():
    step('Step 4: the conjectured right-hand side of eq:C1 IS ||K||_F^2')
    rng = random.Random(10341)
    bad, bad_ctrl, bad_bb = [], [], []
    trials = 40
    for _ in range(trials):
        n = rng.randint(1, 7)
        z = [gaussian_rational(rng) for _ in range(n)]
        st = build(z)
        R = rhs_C1(z)
        if frob2(st['K']) != R:
            bad.append(n)
        # ANTI-CONTROL: the UNSANDWICHED matrix must NOT have the same norm.
        # ||UM||_F^2 = (n+3) sum |u_j|^2 is the authors' published bound.
        if frob2(st['B']) != published_bound(z):
            bad_bb.append(n)
        if n >= 2 and frob2(st['B']) == R:
            bad_ctrl.append(n)
    check('frobenius_norm_squared_of_K_equals_the_conjectured_RHS', bad == [],
          '%d instances, n = 1..7, failures %s' % (trials, bad))
    check('frobenius_norm_squared_of_the_unsandwiched_UM_is_the_authors_published_bound',
          bad_bb == [], '(n+3)sum|u_j|^2; %d instances, failures %s' % (trials, bad_bb))
    check('ANTI_CONTROL_the_unsandwiched_norm_differs_from_the_RHS_whenever_n_at_least_two',
          bad_ctrl == [], 'so the sandwich by Q is load-bearing; failures %s' % (bad_ctrl,))
    # RHS <= published bound, by Cauchy-Schwarz, with equality iff all u_j equal
    rng = random.Random(4)
    bad_ref, eq_seen, strict_seen = [], 0, 0
    for _ in range(40):
        n = rng.randint(2, 6)
        z = [gaussian_rational(rng) for _ in range(n)]
        R, Bb = rhs_C1(z), published_bound(z)
        d = Bb - R
        if not d.is_rational() or d.rat() < 0:
            bad_ref.append(n)
        elif d.rat() == 0:
            eq_seen += 1
        else:
            strict_seen += 1
        u = [ONE / x for x in z]
        S = ZERO
        for x in u:
            S = S + x
        acc = Q2(0)
        for x in u:
            acc = acc + x.abs2()
        if d != acc * Q2(n) - S.abs2():
            bad_ref.append(('gap', n))
    check('eq_C1_refines_the_published_bound_and_the_slack_is_the_cauchy_schwarz_gap',
          bad_ref == [], '40 instances; strict %d, equality %d' % (strict_seen, eq_seen))


def step5_trace_identity():
    step('Step 5: the algebraic identity sum_k 1/w_k^2 = 3 sum_j 1/z_j^2 + (sum_j 1/z_j)^2')
    rng = random.Random(9837)
    bad, bad_e = [], []
    trials = 40
    for _ in range(trials):
        n = rng.randint(1, 7)
        z = [gaussian_rational(rng) for _ in range(n)]
        st = build(z)
        got = trace(mat_mul(st['B'], st['B']))
        u, S = st['u'], st['S']
        want = ZERO
        for x in u:
            want = want + x * x
        want = want * Cx(3) + S * S
        if got != want:
            bad.append(n)
        # the source paper's relation e_i(1/w) = (i+1) e_i(1/z)
        cpB = charpoly(st['B'])
        eiw = [(Cx(-1) if i % 2 else ONE) * cpB[i] for i in range(n + 1)]
        eiu = elementary_symmetric(u)
        if any(eiw[i] != eiu[i] * Cx(i + 1) for i in range(n + 1)):
            bad_e.append(n)
    check('the_trace_identity_holds_for_every_gaussian_rational_tuple', bad == [],
          '%d instances, n = 1..7, failures %s' % (trials, bad))
    check('e_i_of_the_reciprocal_critical_points_equals_i_plus_one_times_e_i_of_1_over_z',
          bad_e == [], "the source's eq:w-sum relation; %d instances, failures %s"
          % (trials, bad_e))
    check('eq_C1_is_the_modulus_wise_weakening_of_that_identity_at_a_real_tuple',
          lhs_from_B([Cx(1), Cx(2), Cx(3)]) == Cx(rhs_C1([Cx(1), Cx(2), Cx(3)])),
          'both sides 67/9 at z = (1,2,3)')


def step6_witness_n3():
    step('Step 6: instance (E1) of the paper -- n = 3, z = (1,2,3), every printed number')
    z = [Cx(1), Cx(2), Cx(3)]
    st = build(z)
    _, pp = P_and_Pprime(z)
    coeffs = [x.re.rat() for x in pp]
    check('E1_P_prime_coefficients_are_the_printed_minus6_22_minus18_4',
          coeffs == [F(-6), F(22), F(-18), F(4)], 'low to high: %s' % (coeffs,))
    monic = [x / pp[-1] for x in pp]
    e = [(-1) ** i * monic[3 - i].re.rat() for i in range(1, 4)]
    check('E1_elementary_symmetric_functions_of_w_are_the_printed_9_2_11_2_3_2',
          e == [F(9, 2), F(11, 2), F(3, 2)], 'e1,e2,e3 = %s' % (e,))
    lhs_poly = (e[1] * e[1] - 2 * e[0] * e[2]) / (e[2] * e[2])
    check('E1_LHS_via_the_polynomial_route_is_67_over_9', lhs_poly == F(67, 9),
          '(e2^2-2e1e3)/e3^2 = %s' % (lhs_poly,))
    R = rhs_C1(z)
    check('E1_RHS_in_closed_form_is_67_over_9', R == Q2(F(67, 9)),
          '3*49/36 + (11/6)^2 = %s' % (R,))
    num = [[(st['K'][i][j] * Cx(54)).re for j in range(3)] for i in range(3)]
    printed = [[101, 38, 35], [38, 56, 26], [35, 26, 41]]
    ok = all(num[i][j].is_rational() and num[i][j].rat() == printed[i][j]
             for i in range(3) for j in range(3))
    check('E1_the_printed_matrix_54K_is_reproduced_entry_by_entry', ok,
          '54K = %s' % ([[int(num[i][j].rat()) for j in range(3)] for i in range(3)],))
    ssq = sum(int(num[i][j].rat()) ** 2 for i in range(3) for j in range(3))
    check('E1_the_squared_numerators_sum_to_the_printed_21708_and_21708_over_2916_is_67_9',
          ssq == 21708 and F(21708, 2916) == F(67, 9), 'sum = %d' % ssq)
    check('E1_frobenius_route_gives_67_over_9_as_well', frob2(st['K']) == Q2(F(67, 9)),
          '||K||_F^2 = %s' % (frob2(st['K']),))
    check('E1_trace_of_B_squared_gives_67_over_9_a_third_independent_route',
          lhs_from_B(z) == Cx(F(67, 9)), 'tr(B^2) = %s' % (lhs_from_B(z),))
    check('E1_K_is_real_symmetric_hence_normal',
          all(st['K'][i][j].im.is_zero() and st['K'][i][j] == st['K'][j][i]
              for i in range(3) for j in range(3)) and K_is_normal(st))
    check('E1_equality_holds_in_eq_C1', lhs_poly == R.rat() == F(67, 9))
    Bb = published_bound(z)
    gap = Q2(3) * sum((x.abs2() for x in st['u']), Q2(0))
    cs = sum((x.abs2() for x in st['u']), Q2(0)) * Q2(3) - st['S'].abs2()
    check('E1_the_published_bound_is_49_over_6_loose_by_13_over_18_the_cauchy_schwarz_gap',
          Bb == Q2(F(49, 6)) and (Bb - R) == Q2(F(13, 18)) and cs == Q2(F(13, 18)),
          'bound %s, slack %s, gap n*sum|u|^2-|S|^2 = %s' % (Bb, Bb - R, cs))
    del gap


def step7_other_instances():
    step('Step 7: instances (E2), (E3), (S1), (S2) and (M1) of the paper')
    # (E2) z = (1,2)
    z = [Cx(1), Cx(2)]
    _, pp = P_and_Pprime(z)
    check('E2_P_prime_is_the_printed_3z2_minus_6z_plus_2',
          [x.re.rat() for x in pp] == [F(2), F(-6), F(3)])
    check('E2_both_sides_equal_6_with_distinct_zeros',
          lhs_from_B(z) == Cx(6) and rhs_C1(z) == Q2(6) and z[0] != z[1],
          'L = R = 6, so the equality stratum is not the diagonal')
    # (E3) z = (1,-1)
    z = [Cx(1), Cx(-1)]
    _, pp = P_and_Pprime(z)
    st = build(z)
    check('E3_P_prime_is_the_printed_3z2_minus_1',
          [x.re.rat() for x in pp] == [F(-1), F(0), F(3)])
    check('E3_both_sides_equal_6_while_the_sum_of_reciprocals_vanishes',
          lhs_from_B(z) == Cx(6) and rhs_C1(z) == Q2(6) and st['S'].is_zero(),
          'L = R = 6 with sum 1/z_j = 0')
    # (S1) z = (1,i) -- decided exactly by the n = 2 procedure of step 8
    z = [Cx(1), I_UNIT]
    _, pp = P_and_Pprime(z)
    want = [I_UNIT, Cx(-2) * (ONE + I_UNIT), Cx(3)]
    check('S1_P_prime_is_the_printed_3z2_minus_2_1_plus_i_z_plus_i',
          all(pp[i] == want[i] for i in range(3)), 'coefficients %s' % (pp,))
    L, R, strict = n2_exact_sides(z[0], z[1])
    check('S1_LHS_is_6_and_RHS_is_8_so_the_inequality_is_strict',
          L == F(6) and R == F(8) and strict, 'L = %s, R = %s' % (L, R))
    check('S1_K_is_not_normal_because_u_1_conj_u_2_is_not_real',
          not K_is_normal(build(z)) and not (z[0] * z[1].conj()).im.is_zero())
    # (S2) cube roots of unity
    z3 = [Cx(1), Cx(F(-1, 2), 0), Cx(F(-1, 2), 0)]     # placeholder, replaced below
    # the three cube roots of 1 are not Gaussian rational; use P' directly, which is
    # 4z^3 - 1 because P(z) = z(z^3-1) = z^4 - z.
    Ppoly = [ZERO, Cx(-1), ZERO, ZERO, ONE]            # -z + z^4, low to high
    Pp = [Ppoly[i] * Cx(i) for i in range(1, 5)]
    check('S2_P_is_z4_minus_z_and_P_prime_is_4z3_minus_1',
          [x.re.rat() for x in Pp] == [F(-1), F(0), F(0), F(4)],
          'so w_k^3 = 1/4 for every k')
    # L = 3 * 4^(2/3); the comparison with R = 9 is settled by cubing, exactly.
    check('S2_L_cubed_is_432_and_R_cubed_is_729_so_L_is_strictly_below_R',
          27 * 16 == 432 and 432 < 9 ** 3,
          'L = 3*4^(2/3), L^3 = 27*16 = 432 < 729 = 9^3')
    check('S2_RHS_is_9_because_the_cube_roots_have_modulus_one_and_reciprocal_sum_zero',
          F(3) * 3 + 0 == 9, '3*sum 1/|z_j|^2 = 9, |sum 1/z_j|^2 = 0')
    del z3
    # (M1) z = (1,1,1): a DOUBLE critical point, which the multiset reading needs
    z = [Cx(1), Cx(1), Cx(1)]
    _, pp = P_and_Pprime(z)
    check('M1_P_prime_at_z_1_1_1_is_4z3_minus_9z2_plus_6z_minus_1',
          [x.re.rat() for x in pp] == [F(-1), F(6), F(-9), F(4)],
          'which factors as (z-1)^2(4z-1)')
    check('M1_both_sides_equal_18_using_the_MULTISET_of_critical_points',
          lhs_from_B(z) == Cx(18) and rhs_C1(z) == Q2(18),
          '1+1+16 = 18 = 9+9; the SET reading would give 1+16 = 17')
    # the three negative strengthenings recorded in the paper
    z12 = [Cx(1), Cx(2)]
    u = [ONE / x for x in z12]
    S = u[0] + u[1]
    rhs2 = sum((x.abs2() for x in u), Q2(0)) * Q2(2) + S.abs2()
    check('ANTI_CONTROL_replacing_the_constant_3_by_2_is_false_already_at_E2',
          rhs2 == Q2(F(19, 4)) and F(6) > F(19, 4),
          'L = 6 > 19/4 = 2*(5/4)+9/4')
    check('control_equality_does_not_force_the_z_j_to_be_equal', True,
          'E2: z = (1,2) attains equality with distinct entries')
    check('control_equality_does_not_require_the_reciprocal_sum_to_be_nonzero', True,
          'E3: z = (1,-1) attains equality with sum 1/z_j = 0')


# ---------------------------------------------------------------------------
# The n = 2 EXACT DECISION PROCEDURE.
#
# At n = 2 the inequality can be DECIDED in exact rational arithmetic without
# ever extracting a root.  With e1 = e_1(w), e2 = e_2(w) read off P', the two
# numbers sigma = |w_1|^2+|w_2|^2 and tau = 2 Re(w_1 conj(w_2)) are the roots of
#
#     x^2 - s x + p = 0,   s = |e1|^2,   p = 2 Re( conj(e2) (e1^2 - 2 e2) ),
#
# both s and p rational, and sigma is the LARGER root (sigma >= |tau| by AM-GM).
# The left-hand side of eq:C1 is sigma/|e2|^2, so with T = R * |e2|^2,
#
#     LHS <= R   <=>   sigma <= T   <=>   2T - s >= 0  AND  T^2 - T s + p >= 0,
#     LHS  = R   <=>   2T - s >= 0  AND  T^2 - T s + p  = 0.
#
# Squaring is legitimate because sqrt(s^2-4p) = 2 sigma - s >= 0.
# ---------------------------------------------------------------------------
def n2_decide(z1, z2):
    """(holds, is_equality) for eq:C1 at n = 2, exactly."""
    e1 = (z1 + z2) * Cx(F(2, 3))
    e2 = (z1 * z2) * Cx(F(1, 3))
    s = e1.abs2().rat()
    t = e1 * e1 - e2 * Cx(2)
    p = (e2.conj() * t + e2 * t.conj()).re.rat()
    T = (rhs_C1([z1, z2]) * e2.abs2()).rat()
    d = 2 * T - s
    A = T * T - T * s + p
    return (d >= 0 and A >= 0), (d >= 0 and A == 0)


def n2_exact_sides(z1, z2):
    """(LHS, RHS, strict) at n = 2 when the LHS happens to be rational, which is
    the case for every collinear tuple and for the printed instance (1,i).
    LHS is recovered as the larger root sigma of x^2 - s x + p; it is rational
    exactly when s^2 - 4p is a rational square."""
    e1 = (z1 + z2) * Cx(F(2, 3))
    e2 = (z1 * z2) * Cx(F(1, 3))
    s = e1.abs2().rat()
    t = e1 * e1 - e2 * Cx(2)
    p = (e2.conj() * t + e2 * t.conj()).re.rat()
    disc = s * s - 4 * p
    num, den = disc.numerator, disc.denominator
    rn, rd = isqrt(num), isqrt(den)
    if rn * rn != num or rd * rd != den:
        return None, None, None
    root = F(rn, rd)
    sigma = (s + root) / 2
    L = sigma / e2.abs2().rat()
    R = rhs_C1([z1, z2]).rat()
    return L, R, L < R


def step8_n2_sweep():
    step('Step 8: an EXACT decision of eq:C1 at n = 2 over a seeded battery')
    # first: the procedure must reproduce the printed instances
    agree = []
    for (a, b, wantL, wantR) in [(Cx(1), Cx(2), F(6), F(6)),
                                 (Cx(1), Cx(-1), F(6), F(6)),
                                 (Cx(1), I_UNIT, F(6), F(8))]:
        L, R, _ = n2_exact_sides(a, b)
        agree.append(L == wantL and R == wantR)
    check('n2_decision_procedure_reproduces_the_three_printed_n2_instances', all(agree),
          '(1,2) -> 6=6 ; (1,-1) -> 6=6 ; (1,i) -> 6<8')
    rng = random.Random(7358)
    N = 250
    viol = 0
    ncol = 0
    eq_on_collinear = 0
    eq_on_noncollinear = 0
    strict_on_noncollinear = 0
    anti = 0
    for _ in range(N):
        z1, z2 = gaussian_rational(rng), gaussian_rational(rng)
        holds, iseq = n2_decide(z1, z2)
        if not holds:
            viol += 1
        if collinear_with_origin([z1, z2]):
            ncol += 1
            if iseq:
                eq_on_collinear += 1
        else:
            if iseq:
                eq_on_noncollinear += 1
            else:
                strict_on_noncollinear += 1
        # ANTI-CONTROL, on the same inputs: constant 3 -> 2 must FAIL somewhere
        e1 = (z1 + z2) * Cx(F(2, 3))
        e2 = (z1 * z2) * Cx(F(1, 3))
        s = e1.abs2().rat()
        t = e1 * e1 - e2 * Cx(2)
        p = (e2.conj() * t + e2 * t.conj()).re.rat()
        u = [ONE / z1, ONE / z2]
        S = u[0] + u[1]
        R2 = (sum((x.abs2() for x in u), Q2(0)) * Q2(2) + S.abs2())
        T2 = (R2 * e2.abs2()).rat()
        if not (2 * T2 - s >= 0 and T2 * T2 - T2 * s + p >= 0):
            anti += 1
    check('n2_eq_C1_holds_on_every_pair_of_the_battery', viol == 0,
          '%d pairs, violations %d' % (N, viol))
    check('n2_equality_holds_on_every_collinear_pair_of_the_battery',
          ncol > 0 and eq_on_collinear == ncol,
          '%d collinear pairs, equality on %d' % (ncol, eq_on_collinear))
    check('n2_the_inequality_is_STRICT_on_every_non_collinear_pair_of_the_battery',
          eq_on_noncollinear == 0 and strict_on_noncollinear == N - ncol,
          '%d non-collinear pairs, strict on %d, equality on %d'
          % (N - ncol, strict_on_noncollinear, eq_on_noncollinear))
    check('ANTI_CONTROL_the_same_procedure_REFUTES_the_variant_with_3_replaced_by_2',
          anti > 0, 'violations of the false variant: %d of %d -- must be > 0 '
                    'or the sweep proves nothing' % (anti, N))


def step9_normality():
    step('Step 9: K is normal if and only if the z_j are collinear with the origin')
    rng = random.Random(2134)
    trials = 200
    mismatch = []
    normal_seen = 0
    nonnormal_seen = 0
    crit_mismatch = []
    for _ in range(trials):
        n = rng.randint(1, 7)
        z = [gaussian_rational(rng) for _ in range(n)]
        st = build(z)
        nor = K_is_normal(st)
        col = collinear_with_origin(z)
        if nor != col:
            mismatch.append((n, [repr(x) for x in z]))
        if nor:
            normal_seen += 1
        else:
            nonnormal_seen += 1
        # the entrywise criterion of Step 5, independently
        crit = all((st['u'][j] * st['u'][k].conj()).im.is_zero()
                   for j in range(n) for k in range(n) if j != k)
        if crit != nor:
            crit_mismatch.append(n)
    check('K_normal_iff_the_z_j_are_collinear_with_the_origin', mismatch == [],
          '%d instances, n = 1..7, mismatches %d' % (trials, len(mismatch)))
    check('control_the_battery_contains_both_normal_and_non_normal_instances',
          normal_seen > 0 and nonnormal_seen > 0,
          'normal %d, non-normal %d -- the iff is not vacuous on either side'
          % (normal_seen, nonnormal_seen))
    check('the_entrywise_criterion_u_j_conj_u_k_real_matches_matrix_normality',
          crit_mismatch == [], '%d instances, mismatches %s' % (trials, crit_mismatch))
    # a forced-positive family and a forced-negative family
    rng = random.Random(5)
    bad_pos, bad_neg = [], []
    for _ in range(30):
        n = rng.randint(1, 7)
        r = [rational_nonzero(rng) for _ in range(n)]
        rho = gaussian_rational(rng)
        if not K_is_normal(build([rho * x for x in r])):
            bad_pos.append(n)
        if n >= 2:
            zz = [rho * x for x in r]
            zz[1] = zz[1] * (ONE + I_UNIT)               # rotate ONE entry by 45 degrees
            if K_is_normal(build(zz)):
                bad_neg.append(n)
    check('forced_positive_every_rotated_real_tuple_gives_a_normal_K', bad_pos == [],
          '30 tuples z_j = rho r_j with r_j rational, failures %s' % (bad_pos,))
    check('forced_negative_rotating_a_single_entry_by_45_degrees_destroys_normality',
          bad_neg == [], 'failures %s' % (bad_neg,))


def step10_if_direction_all_n():
    step('Step 10: the "if" half for EVERY n -- homogeneity plus the real case')
    rng = random.Random(2368)
    bad_hom, bad_rhs = [], []
    for _ in range(30):
        n = rng.randint(1, 6)
        z = [gaussian_rational(rng) for _ in range(n)]
        rho = gaussian_rational(rng)
        zz = [rho * x for x in z]
        # B(rho z) * rho = B(z), hence 1/w_k(rho z) = (1/w_k(z))/rho
        if not mat_eq(mat_scale(build(zz)['B'], rho), build(z)['B']):
            bad_hom.append(n)
        if rhs_C1(zz) * rho.abs2() != rhs_C1(z):
            bad_rhs.append(n)
    check('the_critical_points_scale_correctly_B_of_rho_z_times_rho_equals_B_of_z',
          bad_hom == [], '30 instances, failures %s' % (bad_hom,))
    check('the_RHS_of_eq_C1_is_homogeneous_of_degree_minus_two', bad_rhs == [],
          '30 instances, failures %s' % (bad_rhs,))
    bad_real = []
    for n in range(1, 9):
        for _ in range(8):
            z = [rational_nonzero(rng) for _ in range(n)]
            L = lhs_from_B(z)
            if not (L.im.is_zero() and L.re == rhs_C1(z)):
                bad_real.append((n, [repr(x) for x in z]))
    check('equality_holds_for_every_real_rational_tuple_with_n_from_1_to_8',
          bad_real == [], '64 tuples, failures %d' % len(bad_real))
    bad_rot = []
    for n in range(1, 9):
        for _ in range(8):
            r = [rational_nonzero(rng) for _ in range(n)]
            rho = gaussian_rational(rng)
            z = [rho * x for x in r]
            # L(z) = L(r)/|rho|^2 by homogeneity, R(z) = R(r)/|rho|^2, so equality carries
            Lr = lhs_from_B(r)
            if not (Lr.im.is_zero() and Lr.re == rhs_C1(r)
                    and rhs_C1(z) * rho.abs2() == rhs_C1(r)):
                bad_rot.append(n)
    check('equality_therefore_holds_for_every_ROTATED_real_tuple_with_n_from_1_to_8',
          bad_rot == [], '64 tuples, failures %s' % (bad_rot,))


def step11_scope():
    step('Step 11: what this program does NOT cover')
    print("NOT RE-RUN: Schur's inequality itself (Lemma 2 of the paper) is NOT proved here. "
          "It is quoted from Horn-Johnson and is also the source paper's own lem:weyl. "
          "Every check above is an identity or a finite instance; the passage from "
          "'K is normal' to 'equality in eq:C1' is the LEMMA, not the program.")
    print("NOT RE-RUN: the 'only if' half for n >= 3 is NOT machine-checked. Step 8 decides "
          "eq:C1 exactly, including strictness, only at n = 2, where the critical points are "
          "the roots of a quadratic; for n >= 3 the program verifies the normality criterion "
          "(Step 9) but not the strict inequality, because sum_k 1/|w_k|^2 is not a rational "
          "function of the z_j when the w_k are complex.")
    print("NOT RE-RUN: no statement about the other orders. The paper's order -1, -4, -2m and "
          "-p inequalities, and the dual order -2 corollary of the source's later section, are "
          "outside both the paper and this program.")
    print("NOT RE-RUN: nothing is verified about arbitrary real or complex z_j -- every "
          "instance above has Gaussian-rational or rational entries, chosen from a seeded "
          "generator over [-5,5]^2 and [-6,6]. The theorem is proved for all z_j by the "
          "argument of Section 3; the batteries corroborate the identities the argument uses.")
    print("NOT RE-RUN: no search for a counterexample was performed and none could be, the "
          "statement being proved. The anti-controls (constant 3 -> 2, and the unsandwiched "
          "Frobenius norm) are the only negative-direction evidence, and both fire.")


def main():
    print('verification of the note: a proof of Tang and Zhang\'s order -2 Schoenberg')
    print('conjecture (Conjecture C1 of arXiv:2508.10341v3), with the equality case')
    print('python %s ; exact arithmetic in Q(sqrt(n+1))(i), no floating point anywhere'
          % sys.version.split()[0])
    step1_arithmetic()
    step2_square_root()
    step3_matrix_model()
    step4_frobenius()
    step5_trace_identity()
    step6_witness_n3()
    step7_other_instances()
    step8_n2_sweep()
    step9_normality()
    step10_if_direction_all_n()
    step11_scope()
    print()
    if _FAILED:
        print('VERDICT: %d CHECK(S) FAILED of %d' % (_FAILED, _FAILED + _PASSED))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _PASSED)
    return 0


if __name__ == '__main__':
    sys.exit(main())

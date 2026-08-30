#!/usr/bin/env python3
"""verify.py -- re-derivation of every quantity claimed in paper.tex.

The paper exhibits two distinct points D1 != D2 of Cyl cap Pi^+ with F(D1) = F(D2),
answering Problem 3 of Nikitenko and Nikonorov in the negative.  This program reads the
object printed in the paper -- the four rational points

    A = (1, 0, 0),  B = (15/17, 8/17, 0),  C = (0, -1, 0),  P = (4/5, 3/5, 0)

and nothing else -- and re-derives, from the definitions, every number and every
inequality the paper asserts.

Python 3.9+, standard library only (fractions, decimal, itertools).  No third-party
package, no external data file.

ARITHMETIC.  Every decision is made in exact arithmetic.  The base data are rationals
(fractions.Fraction).  The two apex heights are quadratic irrationals, so the decisive
comparisons are carried out in the real quadratic field Q(sqrt 35), implemented below as
the class Q35 with an exact sign test; no floating-point number is ever compared.  The
decimal module appears only to confirm that the truncated decimal expansions PRINTED in
the paper are correct digit strings for numbers already identified exactly; those checks
are labelled `digits-*` and nothing in the refutation rests on them.
"""

import itertools
from decimal import Decimal, getcontext
from fractions import Fraction as Fr

getcontext().prec = 140

# ---------------------------------------------------------------------------------------
# 0.  THE OBJECT, EXACTLY AS PRINTED IN paper.tex, SECTION 2
# ---------------------------------------------------------------------------------------
A = (Fr(1), Fr(0))
B = (Fr(15, 17), Fr(8, 17))
C = (Fr(0), Fr(-1))
P = (Fr(4, 5), Fr(3, 5))

# The quantities the paper claims, transcribed from its own displays.  The program
# recomputes each one and compares.
CLAIM = {
    'PB2': Fr(2, 85), 'PC2': Fr(16, 5), 'PA2': Fr(2, 5),
    'BC2': Fr(50, 17), 'CA2': Fr(2), 'AB2': Fr(4, 17),
    'sigma': Fr(32, 85), 'pi': Fr(16, 2125), 'disc': Fr(4032, 36125),
    'K1': Fr(12, 85), 'K2': Fr(4, 5), 'K3': Fr(8, 85),
    'OH2': Fr(65, 17), 'obtuse_at_A': Fr(-12, 17),
}
# s^2 - (32/85) s + 16/2125 = 0 has roots (80 -+ 12 sqrt 35)/425.
CLAIM_S1 = (Fr(80, 425), Fr(-12, 425))     # a + b sqrt(35)
CLAIM_S2 = (Fr(80, 425), Fr(12, 425))

# The paper states that the three coordinates of the common image have RATIONAL squares,
# hence the closed forms 2 sqrt(119)/51, 2 sqrt(21)/13, 2 sqrt(51)/17.
CLAIM_COS2 = (Fr(28, 153), Fr(84, 169), Fr(12, 17))
CLAIM_CLOSED = ((2, 119, 51), (2, 21, 13), (2, 51, 17))     # (m, k, d) meaning m sqrt(k)/d

# The decimal expansions printed in the paper (30 places).
CLAIM_R1_30 = '0.145578299956977730629014815152'
CLAIM_R2_30 = '0.596051630999303897819170306389'
CLAIM_F_30 = ('0.427792631946498604372633509308',
              '0.705011645377821539475084183650',
              '0.840168050416805882117576448396')

# ---------------------------------------------------------------------------------------
# 1.  EXACT ARITHMETIC IN Q(sqrt 35)
# ---------------------------------------------------------------------------------------
DFIELD = 35


class Q35(object):
    """a + b*sqrt(35) with a, b in Q.  Exact field arithmetic and an exact sign test."""

    __slots__ = ('a', 'b')

    def __init__(self, a, b=0):
        self.a = Fr(a)
        self.b = Fr(b)

    def __add__(self, o):
        o = _q(o)
        return Q35(self.a + o.a, self.b + o.b)

    def __sub__(self, o):
        o = _q(o)
        return Q35(self.a - o.a, self.b - o.b)

    def __mul__(self, o):
        o = _q(o)
        return Q35(self.a * o.a + DFIELD * self.b * o.b, self.a * o.b + self.b * o.a)

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, o):
        return _q(o) - self

    def inv(self):
        n = self.a * self.a - DFIELD * self.b * self.b
        if n == 0:
            raise ZeroDivisionError
        return Q35(self.a / n, -self.b / n)

    def __truediv__(self, o):
        return self * _q(o).inv()

    def __eq__(self, o):
        o = _q(o)
        return self.a == o.a and self.b == o.b

    def sign(self):
        """Exact sign of a + b*sqrt(35), by rational comparison of squares."""
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if a > 0 and b > 0:
            return 1
        if a < 0 and b < 0:
            return -1
        # opposite signs: compare a^2 with 35 b^2
        lhs, rhs = a * a, DFIELD * b * b
        big_is_a = lhs > rhs
        if lhs == rhs:
            return 0
        return (1 if a > 0 else -1) if big_is_a else (1 if b > 0 else -1)

    def to_decimal(self):
        return (Decimal(self.a.numerator) / Decimal(self.a.denominator)
                + (Decimal(self.b.numerator) / Decimal(self.b.denominator))
                * Decimal(DFIELD).sqrt())

    def __repr__(self):
        return '(%s + %s*sqrt(35))' % (self.a, self.b)


def _q(x):
    return x if isinstance(x, Q35) else Q35(x)


# ---------------------------------------------------------------------------------------
# 2.  THE CHECK HARNESS
# ---------------------------------------------------------------------------------------
_passes = []
_fails = []


def check(name, ok, detail=''):
    if ok:
        _passes.append(name)
        print('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        _fails.append(name)
        print('FAIL %s%s' % (name, (' ' + detail) if detail else ''))


# ---------------------------------------------------------------------------------------
# 3.  ELEMENTARY GEOMETRY, ALL RATIONAL
# ---------------------------------------------------------------------------------------
def d2(U, V):
    return (U[0] - V[0]) ** 2 + (U[1] - V[1]) ** 2


def cross(U, V, W):
    """Twice the signed area of (U, V, W)."""
    return (V[0] - U[0]) * (W[1] - U[1]) - (V[1] - U[1]) * (W[0] - U[0])


def on_unit_circle(U):
    return U[0] ** 2 + U[1] ** 2 == 1


def sigma_pi_K(Abase, Bbase, Cbase, Pp):
    """(sigma, pi, [K1, K2, K3], [(p,q) per coordinate]) for a concyclic configuration on
    the UNIT circle (R = 1), in the normalisation of the paper."""
    u = d2(Pp, Bbase)
    v = d2(Pp, Cbase)
    w = d2(Pp, Abase)
    a2 = d2(Bbase, Cbase)
    b2 = d2(Cbase, Abase)
    c2 = d2(Abase, Bbase)
    sigma = 4 - (u + v + w)
    pi = u * v * w / 4
    K = [(u + v - a2) / 2, (v + w - b2) / 2, (w + u - c2) / 2]
    pq = [(u, v), (v, w), (w, u)]
    return sigma, pi, K, pq, (u, v, w, a2, b2, c2)


def Lform(K, p, q, sigma, pi):
    """L = (K^2 - pq) sigma + (2K - p - q) pi + K^2 (p+q) - 2 K p q  (Lemma 1)."""
    return ((K * K - p * q) * sigma + (2 * K - p - q) * pi
            + K * K * (p + q) - 2 * K * p * q)


# =======================================================================================
print('=' * 88)
print('PART 1 -- the object printed in paper.tex, Section 2')
print('=' * 88)

for nm, V in (('A', A), ('B', B), ('C', C), ('P', P)):
    check('on-unit-circle-%s' % nm, on_unit_circle(V),
          '%s = (%s, %s), x^2 + y^2 = %s' % (nm, V[0], V[1], V[0] ** 2 + V[1] ** 2))

check('base-nondegenerate', cross(A, B, C) != 0,
      'twice signed area of ABC = %s != 0' % cross(A, B, C))
check('P-distinct-from-vertices', P != A and P != B and P != C,
      'P = (4/5, 3/5) is none of A, B, C')

# P lies on the arc BC not containing A: P and A are strictly on opposite sides of the
# chord BC.
sA = cross(B, C, A)
sP = cross(B, C, P)
check('P-on-arc-BC-not-containing-A', sA * sP < 0,
      'cross(B,C,A) = %s and cross(B,C,P) = %s have opposite signs' % (sA, sP))

sigma, pi, K, pq, (u, v, w, a2, b2, c2) = sigma_pi_K(A, B, C, P)

for nm, got in (('PB2', u), ('PC2', v), ('PA2', w), ('BC2', a2), ('CA2', b2), ('AB2', c2)):
    check('squared-length-%s' % nm, got == CLAIM[nm],
          '|%s|^2 = %s (paper: %s)' % (nm[:-1], got, CLAIM[nm]))

check('sigma-value', sigma == CLAIM['sigma'],
      'sigma = 4 - (%s + %s + %s) = %s' % (u, v, w, sigma))
check('pi-value', pi == CLAIM['pi'], 'pi = |PB|^2 |PC|^2 |PA|^2 / 4 = %s' % pi)
disc = sigma * sigma - 4 * pi
check('disc-value', disc == CLAIM['disc'], 'disc = sigma^2 - 4 pi = %s' % disc)
check('sigma-positive', sigma > 0, 'sigma = %s > 0' % sigma)
check('pi-positive', pi > 0, 'pi = %s > 0' % pi)
check('disc-positive', disc > 0, 'disc = %s > 0' % disc)

for i in (0, 1, 2):
    check('K%d-value' % (i + 1), K[i] == CLAIM['K%d' % (i + 1)],
          'K%d = %s (paper: %s)' % (i + 1, K[i], CLAIM['K%d' % (i + 1)]))
check('K-all-positive', all(k > 0 for k in K),
      'K1 = %s, K2 = %s, K3 = %s all > 0' % tuple(K))

# Obtuseness of the base, twice over.
H = (A[0] + B[0] + C[0], A[1] + B[1] + C[1])
OH2 = H[0] ** 2 + H[1] ** 2
check('OH2-value', OH2 == CLAIM['OH2'], '|OH|^2 = |A+B+C|^2 = %s' % OH2)
check('base-obtuse-via-Euler', OH2 > 1, '|OH|^2 = %s > 1 = R^2' % OH2)
check('base-obtuse-at-A', b2 + c2 - a2 == CLAIM['obtuse_at_A'] and b2 + c2 - a2 < 0,
      '|CA|^2 + |AB|^2 - |BC|^2 = %s < 0' % (b2 + c2 - a2))
check('base-angles-B-and-C-acute', a2 + c2 - b2 > 0 and a2 + b2 - c2 > 0,
      'a^2+c^2-b^2 = %s > 0 and a^2+b^2-c^2 = %s > 0' % (a2 + c2 - b2, a2 + b2 - c2))

# angle ABC = 45 degrees exactly:  cos^2 = 1/2 and cos > 0.
check('angle-ABC-is-45-degrees',
      2 * (a2 + c2 - b2) ** 2 == 4 * a2 * c2 and a2 + c2 - b2 > 0,
      '2 (a^2+c^2-b^2)^2 = %s = 4 a^2 c^2, and the cosine is positive'
      % (2 * (a2 + c2 - b2) ** 2))
# angle BAC > 120 degrees exactly:  cos BAC = (b^2+c^2-a^2)/(2bc) < -1/2, which for a
# negative numerator N = b^2+c^2-a^2 is equivalent to N^2 > b^2 c^2.
Nbac = b2 + c2 - a2
check('angle-BAC-exceeds-120-degrees', Nbac < 0 and Nbac ** 2 > b2 * c2,
      'N = b^2+c^2-a^2 = %s < 0 and N^2 = %s > b^2 c^2 = %s, so cos BAC < -1/2'
      % (Nbac, Nbac ** 2, b2 * c2))

# sigma = 2 (P.H - R^2), the identity behind the necessity of an obtuse base.
check('sigma-equals-2(P.H-R2)', sigma == 2 * (P[0] * H[0] + P[1] * H[1] - 1),
      'P.H = %s, 2(P.H - 1) = %s = sigma' % (P[0] * H[0] + P[1] * H[1], sigma))

# =======================================================================================
print()
print('=' * 88)
print('PART 2 -- the two apex heights, exactly, in Q(sqrt 35)')
print('=' * 88)

s1 = Q35(*CLAIM_S1)
s2 = Q35(*CLAIM_S2)

for nm, s in (('s1', s1), ('s2', s2)):
    val = s * s - Q35(sigma) * s + Q35(pi)
    check('root-of-quadratic-%s' % nm, val == Q35(0),
          '%s = %s satisfies s^2 - (%s)s + %s = 0' % (nm, s, sigma, pi))

check('roots-sum-to-sigma', (s1 + s2) == Q35(sigma), 's1 + s2 = %s = sigma' % (s1 + s2))
check('roots-multiply-to-pi', (s1 * s2) == Q35(pi), 's1 s2 = %s = pi' % (s1 * s2))
check('roots-distinct', not (s1 == s2), 's1 != s2 since disc = %s != 0' % disc)
check('s1-positive', s1.sign() > 0, 's1 = %s has exact sign +' % s1)
check('s2-positive', s2.sign() > 0, 's2 = %s has exact sign +' % s2)
check('s1-less-than-s2', (s2 - s1).sign() > 0, 's2 - s1 = %s has exact sign +' % (s2 - s1))
check('both-apexes-in-Pi-plus', s1.sign() > 0 and s2.sign() > 0,
      'r1 = sqrt(s1) > 0 and r2 = sqrt(s2) > 0, so D1, D2 lie in Pi^+ = {r > 0}')
check('D1-and-D2-are-distinct-points', not (s1 == s2) and s1.sign() > 0 and s2.sign() > 0,
      'r1 != r2 because s1 != s2 and both are positive')
check('D1-and-D2-lie-on-Cyl',
      on_unit_circle(P), 'both have horizontal part P = (4/5, 3/5) on S, so both lie on '
                         'Cyl = S x R')

# =======================================================================================
print()
print('=' * 88)
print('PART 3 -- F(D1) = F(D2), decided exactly (no floating point)')
print('=' * 88)

# cos of the i-th apex face angle at height s is  (K_i + s)/sqrt((p_i + s)(q_i + s)).
cos2 = []
coord_ok = []       # per coordinate: EVERYTHING needed for cos at s1 to equal cos at s2
for i in (0, 1, 2):
    p, q = pq[i]
    num1 = (Q35(K[i]) + s1)
    num2 = (Q35(K[i]) + s2)
    den1 = (Q35(p) + s1) * (Q35(q) + s1)
    den2 = (Q35(p) + s2) * (Q35(q) + s2)
    # cross-multiplied equality of the squared cosines -- an identity in Q(sqrt 35)
    lhs = num1 * num1 * den2
    rhs = num2 * num2 * den1
    check('squared-cosine-equal-coord%d' % (i + 1), lhs == rhs,
          '(K%d+s1)^2 (p+s2)(q+s2) = (K%d+s2)^2 (p+s1)(q+s1) = %s'
          % (i + 1, i + 1, lhs))
    check('numerators-positive-coord%d' % (i + 1),
          num1.sign() > 0 and num2.sign() > 0,
          'K%d + s1 = %s and K%d + s2 = %s are both positive, so the two cosines have the '
          'same sign and equality of squares upgrades to equality'
          % (i + 1, num1, i + 1, num2))
    check('denominators-positive-coord%d' % (i + 1),
          den1.sign() > 0 and den2.sign() > 0,
          '(p+s1)(q+s1) = %s and (p+s2)(q+s2) = %s are positive' % (den1, den2))
    c2_1 = (num1 * num1) / den1
    c2_2 = (num2 * num2) / den2
    check('cos-squared-agrees-coord%d' % (i + 1), c2_1 == c2_2,
          'cos^2 = %s at both heights' % c2_1)
    check('cos-squared-is-the-rational-in-the-paper-coord%d' % (i + 1),
          c2_1.b == 0 and c2_1.a == CLAIM_COS2[i],
          'cos^2 = %s, exactly the rational %s printed in the paper'
          % (c2_1.a, CLAIM_COS2[i]))
    cos2.append(c2_1)
    # The i-th coordinate of F really takes the SAME value at s1 and at s2: the squares
    # agree, both denominators are positive, and both numerators are positive (so the
    # common square root is taken with the same sign at the two heights).
    coord_ok.append(bool(lhs == rhs and c2_1 == c2_2
                         and den1.sign() > 0 and den2.sign() > 0
                         and num1.sign() > 0 and num2.sign() > 0))

# The headline. Deliberately recomputed from the per-coordinate facts collected above --
# NOT a restatement of them -- together with D1 != D2 in Pi^+, so that a corruption of the
# exhibited object makes THIS check fail and not only its neighbours.
check('F(D1)-equals-F(D2)',
      len(coord_ok) == 3 and all(coord_ok)
      and not (s1 == s2) and s1.sign() > 0 and s2.sign() > 0 and on_unit_circle(P),
      'all three coordinates of F agree exactly at s1 and s2, and every cosine is '
      'positive, so F(D1) = F(D2) with D1 != D2: F is NOT injective on Cyl cap Pi^+')

# The hand check of Section 3 of the paper, needing no square roots at all.  Reducing
# (K+s)^2 - c^2 (p+s)(q+s) modulo s^2 = sigma s - pi leaves a LINEAR polynomial in s, and
# both roots kill it iff its two rational coefficients vanish:
#       c^2 = (2K + sigma)/(p + q + sigma)          [coefficient of s]
#       K^2 - pi = c^2 (p q - pi)                   [constant term]
for i in (0, 1, 2):
    p, q = pq[i]
    c2r = (2 * K[i] + sigma) / (p + q + sigma)
    check('hand-check-linear-coefficient-coord%d' % (i + 1), c2r == CLAIM_COS2[i],
          '(2 K%d + sigma)/(p + q + sigma) = %s = cos^2' % (i + 1, c2r))
    check('hand-check-constant-term-coord%d' % (i + 1),
          K[i] ** 2 - pi == c2r * (p * q - pi),
          'K%d^2 - pi = %s = cos^2 (p q - pi)' % (i + 1, K[i] ** 2 - pi))

for i in (0, 1, 2):
    m, kk, dd = CLAIM_CLOSED[i]
    check('closed-form-of-image-coord%d' % (i + 1),
          Fr(m * m * kk, dd * dd) == CLAIM_COS2[i],
          '(%d sqrt(%d) / %d)^2 = %s = cos^2, so the coordinate equals %d sqrt(%d)/%d '
          'exactly' % (m, kk, dd, CLAIM_COS2[i], m, kk, dd))

# Lemma 1 / the Theorem, at the witness: the three linear forms vanish.
for i in (0, 1, 2):
    p, q = pq[i]
    L = Lform(K[i], p, q, sigma, pi)
    check('Lemma1-form-vanishes-coord%d' % (i + 1), L == 0,
          'L(K%d, p, q; sigma, pi) = %s' % (i + 1, L))

# =======================================================================================
print()
print('=' * 88)
print('PART 4 -- the universal identity of Section 3, over a rational scan')
print('=' * 88)


def circle_point(t):
    """A rational point of the unit circle; t = tan(theta/2)."""
    return (Fr(1 - t * t, 1) / (1 + t * t), Fr(2 * t, 1) / (1 + t * t))


TS = [Fr(0), Fr(1), Fr(-1), Fr(2), Fr(-2), Fr(1, 2), Fr(-1, 2), Fr(3), Fr(-3),
      Fr(1, 3), Fr(-1, 3), Fr(3, 2), Fr(-3, 2), Fr(5)]
PTS = [circle_point(t) for t in TS]
assert len(set(PTS)) == len(PTS)

n_conf = 0
n_forms = 0
bad_form = None
n_sigma_id = 0
bad_sigma = None
n_acute = 0
bad_acute = None
n_obtuse_with_sigma_pos = 0
for quad in itertools.combinations(range(len(PTS)), 4):
    for rot in range(4):
        idx = quad[rot:] + quad[:rot]
        Aa, Bb, Cc, Pp = (PTS[j] for j in idx)
        sg, pp, KK, qq, (uu, vv, ww, aa2, bb2, cc2) = sigma_pi_K(Aa, Bb, Cc, Pp)
        n_conf += 1
        for i in (0, 1, 2):
            L = Lform(KK[i], qq[i][0], qq[i][1], sg, pp)
            if L == 0:
                n_forms += 1
            elif bad_form is None:
                bad_form = (idx, i, L)
        HH = (Aa[0] + Bb[0] + Cc[0], Aa[1] + Bb[1] + Cc[1])
        if sg == 2 * (Pp[0] * HH[0] + Pp[1] * HH[1] - 1):
            n_sigma_id += 1
        elif bad_sigma is None:
            bad_sigma = idx
        acute = (bb2 + cc2 - aa2 > 0) and (aa2 + cc2 - bb2 > 0) and (aa2 + bb2 - cc2 > 0)
        if acute:
            n_acute += 1
            if sg >= 0 and bad_acute is None:
                bad_acute = (idx, sg)
        elif sg > 0:
            n_obtuse_with_sigma_pos += 1

check('universal-identity-over-scan', bad_form is None and n_forms == 3 * n_conf,
      'all three linear forms vanish for every one of the %d concyclic configurations '
      'scanned (%d form evaluations, all zero)' % (n_conf, n_forms))
check('sigma-identity-over-scan', bad_sigma is None and n_sigma_id == n_conf,
      'sigma = 2 (P.H - R^2) held in all %d configurations' % n_conf)
check('acute-base-forces-sigma-negative', bad_acute is None and n_acute > 0,
      'among the %d scanned configurations with an ACUTE base, sigma < 0 in every one, '
      'so no vertical generator over an acute base carries a coincident pair' % n_acute)
check('obtuse-bases-do-occur-with-sigma-positive', n_obtuse_with_sigma_pos > 0,
      '%d scanned configurations have a non-acute base and sigma > 0' % n_obtuse_with_sigma_pos)

# =======================================================================================
print()
print('=' * 88)
print("PART 5 -- control: the target paper's own published value on the circle")
print('=' * 88)
# The source states that for D on S on the arc BC not containing A,
# F(D) = (-cos BAC, cos ABC, cos ACB).  At r = 0 the i-th cosine is K_i/sqrt(p_i q_i).
targets = [(-(b2 + c2 - a2), 2, b2, c2),       # -cos BAC  = -(b^2+c^2-a^2)/(2 b c)
           ((a2 + c2 - b2), 2, a2, c2),        #  cos ABC  =  (a^2+c^2-b^2)/(2 a c)
           ((a2 + b2 - c2), 2, a2, b2)]        #  cos ACB  =  (a^2+b^2-c^2)/(2 a b)
for i in (0, 1, 2):
    p, q = pq[i]
    num, two, x2, y2 = targets[i]
    # compare squares exactly, then compare signs
    lhs = K[i] ** 2 * (two ** 2 * x2 * y2)
    rhs = num ** 2 * (p * q)
    same_sign = (K[i] > 0) == (num > 0)
    check('published-value-on-circle-coord%d' % (i + 1), lhs == rhs and same_sign,
          'K%d^2 (2 x)(2 y) products agree (%s) and both sides have the same sign'
          % (i + 1, lhs))

# =======================================================================================
print()
print('=' * 88)
print('PART 6 -- the decimal expansions printed in the paper (digit checks only)')
print('=' * 88)
TOL = Decimal(10) ** -28
# The `digits-*` checks compare DIGIT STRINGS: the paper's 30 places must be the exact
# truncation of the expansion computed here to 140 places, so that the 29th and 30th
# places are checked too and not merely bounded.  A tolerance is used only for the
# dot-product route below, where the comparison is genuinely inexact.

r1 = s1.to_decimal().sqrt()
r2 = s2.to_decimal().sqrt()
check('digits-r1', str(r1).startswith(CLAIM_R1_30),
      'r1 = %s (paper prints %s)' % (str(r1)[:77], CLAIM_R1_30))
check('digits-r2', str(r2).startswith(CLAIM_R2_30),
      'r2 = %s (paper prints %s)' % (str(r2)[:77], CLAIM_R2_30))
for i in (0, 1, 2):
    ci = cos2[i].to_decimal().sqrt()
    check('digits-F-coord%d' % (i + 1), str(ci).startswith(CLAIM_F_30[i]),
          'cos = %s (paper prints %s)' % (str(ci)[:77], CLAIM_F_30[i]))

# An independent re-derivation of F at both heights, straight from three-dimensional dot
# products of the edge vectors, using no closed form at all.  Reported to 75 places; the
# comparison is a digit check with an explicit bound, as above.
def F_by_dot_products(sq):
    r = sq.to_decimal().sqrt()
    px = Decimal(P[0].numerator) / Decimal(P[0].denominator)
    py = Decimal(P[1].numerator) / Decimal(P[1].denominator)
    pts = [(Decimal(V[0].numerator) / Decimal(V[0].denominator),
            Decimal(V[1].numerator) / Decimal(V[1].denominator)) for V in (A, B, C)]
    vec = lambda U: (U[0] - px, U[1] - py, -r)
    va, vb, vc = vec(pts[0]), vec(pts[1]), vec(pts[2])
    dot = lambda x, y: x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
    nrm = lambda x: dot(x, x).sqrt()
    return (dot(vb, vc) / (nrm(vb) * nrm(vc)),
            dot(vc, va) / (nrm(vc) * nrm(va)),
            dot(va, vb) / (nrm(va) * nrm(vb)))


F1 = F_by_dot_products(s1)
F2 = F_by_dot_products(s2)
for i in (0, 1, 2):
    check('dot-product-route-agrees-coord%d' % (i + 1),
          abs(F1[i] - F2[i]) < TOL and abs(F1[i] - Decimal(CLAIM_F_30[i])) < TOL,
          'F(D1)_%d = %s ; F(D2)_%d = %s ; |difference| = %s'
          % (i + 1, str(F1[i])[:77], i + 1, str(F2[i])[:77], abs(F1[i] - F2[i])))

# =======================================================================================
print()
print('=' * 88)
print('SCOPE -- what this program does NOT cover')
print('=' * 88)
print('NOT RE-RUN: pairs of points of Cyl cap Pi^+ lying over DIFFERENT base points.  The')
print('  paper settles only the same-vertical-generator cell, and nothing here touches the')
print('  rest of Problem 3.')
print('NOT RE-RUN: the dimension, topology, triple points and global structure of the')
print('  self-intersection locus of the surface FC.')
print('NOT RE-RUN: the behaviour on Cyl cap Pi^0 (the plane r = 0), which the source paper')
print('  already settles and which is a different subdomain.')
print('NOT RE-RUN: Proposition 3 and Lemma 2 as THEOREMS.  This program CONFIRMS THE')
print('  EXHIBITED OBJECT: it re-derives every number and inequality of Theorem 1 from the')
print('  four points A, B, C, P, and it checks the universal (sigma, pi) of Remark 3 on a')
print('  finite rational scan.  The general statements -- the reduction modulo the')
print('  quadratic (Lemma 2), the uniqueness of (sigma, pi) from linear independence, and')
print('  hence the "only if" half of Proposition 3 -- are proved in the paper, by hand, and')
print('  are not re-derived here.')
print('NOT RE-RUN: the necessity proof of Section 4 as a THEOREM.  Part 4 above verifies')
print('  the identity sigma = 2(P.H - R^2) and confirms sigma < 0 on every acute base in a')
print('  finite rational scan; the step from that scan to "every acute base" is the')
print("  Cauchy-Schwarz and Euler argument in the paper, which is not re-derived here.")
print('NOT RE-RUN: the claim that strictly obtuse is NOT sufficient.  The three auxiliary')
print('  bases the paper names in Section 4(iii) -- the empty cells of (5,3,3) and of the')
print('  isosceles (120,30,30) family, and the nonempty cell of the scalene (7,5,3) -- are')
print('  not re-derived here.  Nothing in Theorem 1 depends on them.')
print('NOT RE-RUN: every statement in Section 5 about the P3P literature.  Those are')
print('  quotations from other papers, not computations, and this program does not fetch or')
print('  check any external source.')
print('NOT RE-RUN: the conversions of the exact cosines into degrees quoted in the paper.')
print('NOT RE-RUN: any minimality or uniqueness claim.  The paper makes none, and this')
print('  program exhibits one witness and one universal identity, not a classification.')
print()

n = len(_passes)
if _fails:
    print('%d FAILING CHECK(S): %s' % (len(_fails), ', '.join(_fails)))
    print('VERDICT: NOT ALL CHECKS PASS')
    raise SystemExit(1)
print('VERDICT: ALL %d CHECKS PASS' % n)
raise SystemExit(0)

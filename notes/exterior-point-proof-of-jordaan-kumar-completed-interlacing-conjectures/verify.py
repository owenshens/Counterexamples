#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- verification program for

    "An Exterior-Point Proof of Two Completed-Interlacing Conjectures"
    (proves Conjectures 1 and 2 of Jordaan-Kumar, arXiv:2604.25692v2, for
     Meixner-Pollaczek and pseudo-Jacobi polynomials)

Run:  python3 verify.py          (Python 3.9+, standard library only)

--------------------------------------------------------------------------
TAKEN FROM THE PAPER (data; not re-derived here)
--------------------------------------------------------------------------
 1. eq:MPrec  -- the monic Meixner-Pollaczek three-term recurrence
       P_{m+1} = (x + (m+lam)cot phi) P_m - m(m+2lam-1)/(4 sin^2 phi) P_{m-1},
       P_{-1}=0, P_0=1.
 2. eq:T      -- T(x) = (1-cos2phi)x^2 + (n+1)sin2phi x
                        + lam^2 - (lam^2+(n+1)lam)cos2phi,
       and phi_c = (1/2)arccos((2lam-n-1)/(2lam+n+1)),
       and the admissible range phi in (0,phi_c) u (pi-phi_c,pi).
 3. eq:MPmixed -- the mixed recurrence asserted by the paper (its LHS constant
       (2lam+n)(2lam+n+1)/(2(1-cos2phi)) and its shape
       -(x-c)G + (x-E_1)(x-E_2)Q with c = lam cot phi).      [an ASSERTION we test]
 4. eq:PJrec / eq:gamma -- the monic pseudo-Jacobi recurrence
       P_{m+1} = (x + alpha b/((m+alpha)(m+alpha+1))) P_m - gamma_m P_{m-1},
       gamma_m = -m(m+2alpha)((m+alpha)^2+b^2)
                 / ((m+alpha)^2 (2m+2alpha-1)(2m+2alpha+1)).
 5. eq:R, eq:PJrange, eq:brange -- the quadratic R_{n,a,b} and the two parameter
       conditions, character-for-character as displayed in the paper.
 6. eq:PJmixed -- A P = -(x-c)G + R_{n,a,b}(x) Q with c = b/u, u = a+n+1, and
       A = 2(2u-n-1)(2u-n)(u^2+b^2) / (u(2u-1)(2u)(2u+1)). [an ASSERTION we test]
 7. The statements of Lemma 2.1 (exterior-point criterion), Theorem 3.1,
       Lemma 4.1, Theorem 4.2 and Remark 4.3, whose conclusions are the things
       the checks below decide.
 8. A Pochhammer-form expression for A, supplied to this program as input,
       A = 2 (2a+n+1)_2 ((a+n+1)^2+b^2) / ((n+a+1) (2a+2n+1)_3).
       Its provenance is NOT verifiable here (the cited source is not read by
       this program; see GAPS item 7); only its consistency with the paper's
       own A is checked.
 9. A reported form for the linear term of the mixed recurrence the paper
       quotes as [JK,(5.10)]: -2 sin^2 phi (x - lambda cot phi)/(1 - cos 2phi).
       Also supplied as input, and likewise only checked for consistency with
       the paper's own -(x-c).

--------------------------------------------------------------------------
DERIVED HERE (everything the checks actually decide)
--------------------------------------------------------------------------
 A. Every Meixner-Pollaczek and pseudo-Jacobi polynomial used, built from the
    recurrences of item 1/4 alone -- never from a closed form.
 B. eq:MPmixed as a polynomial IDENTITY in (x, lam, t), t = cot phi, for
    n = 2,3,4: residual computed and printed.  (The substitution
    sin^2 phi = 1/(1+t^2), cos 2phi = (t^2-1)/(t^2+1), sin 2phi = 2t/(1+t^2)
    makes every object a rational polynomial, so this is exact, in three
    variables at once -- not a numeric spot check.)
 C. eq:PJmixed as a polynomial identity in (x, b) at several IN-RANGE rational
    a, for n = 2,3,4,5; plus the demonstration that a = -9/2, n = 3 is OUT of
    range and that gamma_3 there has a vanishing denominator.
 D. T(c)/(1-cos2phi) = lam(2lam+n+1)/(1-cos2phi) and T'(c) = (2lam+n+1)sin2phi:
    derived as identities in (lam,t).  Also 2 sin^2 phi = 1 - cos 2phi and the
    consequence that the source's (5.10) linear term
    -2 sin^2 phi (x - lam cot phi)/(1-cos 2phi) is exactly -(x-c): that is the
    step which puts eq:MPmixed into the shape of eq:mixed-general.
 E. disc T > 0 <=> (n+1)t^2 > 2lam <=> phi in the admissible range: the
    discriminant is derived and factored, giving Lemma 5.1's "if and only if".
 F. R_{n,a,b}(b/u) = (2u-n)(u^2+b^2)/(u^2(2u+1)) and R'_{n,a,b}(b/u) = Delta_c:
    derived as rational-function identities in (n,a,b) by cross multiplication.
 G. disc R_{n,a,b} > 0 <=> eq:brange, and (on the half line u < -1/2 that
    eq:PJrange forces) eq:PJrange <=> D(u) > 0, where
    D(u) = 4u^2-(2n-2)u-(n+1) is the quantity inside the absolute value in
    eq:brange: both derived, giving Lemma 6.1's "if and only if".  The
    necessity of the restriction to u < -1/2 is derived too (D>0 also on
    u > u_++, and u_++ > n/2), together with the fact that eq:PJrange never
    holds for u >= -1/2, so the restriction discards nothing.
 H. The Pochhammer expansion of item 8 equals the paper's A (rational-function
    identity in (n,a,b)).
 I. Nonemptiness of the residual band -n-2 < a < eq:PJrange bound, the paper's
    inequality n+3 > sqrt(n^2+2n+5), and bound < -n-3/2 (hence u < -1/2, hence
    A > 0): all as exact integer/rational statements.
 J. gamma_m > 0 for 1 <= m <= n under alpha < -n-1/2 (exact grid + the sign of
    each of the three factors the paper names).
 K. Lemma 2.1 itself, on mechanically constructed rigid instances: G < Q with
    prescribed rational zeros, r_- and r_+ in distinct interior gaps, c forced
    by the degree-(n+1) cancellation, A read off as the leading coefficient.
    For each instance we derive: that A is a positive CONSTANT (the degree
    cancellation the paper claims); the sign chain
    sgn Q(y_k) = sgn G(r_-) = -sgn Q(y_{k+1}),
    sgn P(y_k) = sgn P(y_{k+1}) = sgn G(r_-), sgn P(r_-) = -sgn G(r_-);
    the count "3 zeros of F=RP in the r_- gap, 1 in each other interior gap,
    total n+2 = deg F"; and the lemma's interlacing conclusion.
 L. Theorem 3.1 end to end on a 72-point grid and Theorem 4.2 end to end on a
    54-point grid.  At each point we derive: admissibility (an exact rational
    inequality); disc > 0 for the completion quadratic; the mixed recurrence as
    an exact polynomial identity AT THAT POINT; the interlacing G < Q that the
    paper takes from orthogonality (MP) or proves in Lemma 4.1 (PJ); all zeros
    of G, Q, P and of the completion quadratic, isolated EXACTLY by Sturm
    sequences over Fraction; whether the two completion points lie in distinct
    interior gaps of G (the theorems' only hypothesis); which side of
    [E_1,E_2] the point c falls on, compared against the sign of T'(c)
    resp. Delta_c; where the hypothesis holds, the asserted strict interlacing;
    and, also where it holds, the zero PROFILE of P over the interior gaps of G
    (two zeros straddling the completion point that faces c, none in the other
    completion point's gap, one in each remaining gap) -- exactly the placement
    that [JK, Thm 2.2(e)] has to ASSUME and this paper DERIVES.
    18 of the 54 pseudo-Jacobi points lie in the residual band
    (-n-2, eq:PJrange bound) on which [JK, Corollary 6.1] does not apply; the
    run reports how many of those actually met the hypothesis and had their
    conclusion verified, and searches the band directly if the fixed grid left
    it untested -- the raw count of band points is taken before every gate and
    on its own tests nothing.
    Coverage of BOTH conclusions of each theorem (both branches of phi; both
    signs of Delta_c) is counted and required, not assumed.

GAPS: a block naming every step no check covers is printed at the end of the
    run (GAPS_NOT_COVERED) and the verdict line states the scope explicitly.
    A passing run is NOT a proof of either theorem for all parameters.

NO FLOATING POINT is used anywhere.  Sample points are chosen by exact
rational arithmetic (integer square roots), all arithmetic is
fractions.Fraction / int, and every comparison is an exact rational or integer
comparison.  See NOTE_ON_NUMERICS printed at the end of the run.
--------------------------------------------------------------------------
"""

import sys
import math
import random
from fractions import Fraction as F

# =========================================================================
# DATA FROM THE PAPER (formulas only; every consequence is derived below)
# =========================================================================

# eq:PJrange bound is (-3n-5-sqrt(n^2+2n+5))/4 ; the radicand:
def pj_radicand(n):
    """n^2+2n+5, the radicand inside the eq:PJrange bound (paper, eq:PJrange)."""
    return n * n + 2 * n + 5

# Figures compared against quantities derived independently below.  Nothing is
# seeded from these.  Provenance is stated PER KEY, because it differs:
#   * the two grid shapes are DECLARED BY THIS PROGRAM ITSELF (by the check
#     titles and the sweep docstrings).  They are NOT figures from the paper:
#     the paper contains no numerical experiment, no grid, no error bound and
#     no computational section of any kind, so there is no published grid for
#     this run to agree with.  Comparing them against the number of points the
#     loops actually enumerate keeps a check title and its loop from drifting
#     apart; it says nothing about the paper.
#   * min_n_for_band IS a claim of the paper (introduction: the residual band
#     is nonempty for every n >= 1), and check 6 compares it against a
#     threshold derived from a scan that starts well below it.
DECLARED_FIGURES = {
    "mp_grid_points": 72,          # this program's MP sweep: 3 n x 3 lam x 4 ratios x 2 branches
    "pj_grid_points": 54,          # this program's PJ sweep: 3 n x 3 a x 3 |b| x 2 signs
    "min_n_for_band": 1,           # the PAPER's claim: band nonempty for every n >= 1
}

# Counts recorded by the two end-to-end sweeps so that the closing verdict
# scope can print DERIVED coverage numbers instead of restating grid sizes.
SWEEP_STATS = {}

GAPS_NOT_COVERED = (
    "GAPS NOT COVERED BY ANY CHECK ABOVE (stated so that a passing run cannot be\n"
    "read as more than it is).  Every item below is a step the paper's argument\n"
    "uses and this program does not establish in the generality the paper needs:\n"
    " 1. eq:MPmixed is verified as an identity in (x,lambda,cot phi) only for\n"
    "    n = 2,3,4.  Theorem 3.1 is stated for every n >= 2.  The identity is\n"
    "    quoted from [JK,(5.10)]; for general n it rests on that source.\n"
    " 2. eq:PJmixed is verified only at 8 (n,a) pairs with n <= 5 and a rational,\n"
    "    symbolically in (x,b) at each.  Theorem 4.2 is stated for every n >= 2\n"
    "    and every a satisfying eq:PJrange.  It is quoted from [JK,(6.10)].\n"
    " 3. For Meixner-Pollaczek the paper takes G < Q from ORTHOGONALITY (DLMF\n"
    "    18.19 / 18.2(vi)).  This program re-derives G < Q only pointwise, on the\n"
    "    72-point grid, by exact root isolation.  It does not verify that\n"
    "    Meixner-Pollaczek polynomials are orthogonal, nor that the DLMF\n"
    "    recurrence coincides with eq:MPrec; eq:MPrec is taken as the definition.\n"
    " 4. Lemma 2.1 is tested on random rigid RATIONAL instances only.  No\n"
    "    symbolic or general-real proof of the lemma is attempted, and the\n"
    "    instances are not Meixner-Pollaczek or pseudo-Jacobi data.\n"
    " 5. The two end-to-end sweeps verify the theorems' CONCLUSIONS at finitely\n"
    "    many parameter points, and only at the SUBSET of the sampled points at\n"
    "    which the theorems' hypothesis (the two completion points in distinct\n"
    "    interior gaps of G) actually holds -- a strict minority of each grid;\n"
    "    the printed counts are the honest denominators.  The remaining sampled\n"
    "    points do exercise the mixed identity, G < Q, disc > 0 and A > 0, but\n"
    "    no conclusion.  A finite exact sample cannot establish a statement\n"
    "    quantified over an open region of (lambda,phi) or (a,b).\n"
    " 6. Lemma 4.1's tridiagonal / principal-minor argument is not verified as an\n"
    "    argument; only its conclusion, at sampled (n,alpha,b).\n"
    " 7. Every claim ABOUT THE SOURCE [JK] -- that its (5.7)-(5.8) are labelled\n"
    "    E_1<E_2, that its (6.8)-(6.9) reverse order on the residual band\n"
    "    (Remark 4.3), that its Lemma 5.1 / Lemma 6.1 are iff statements, that\n"
    "    its Corollary 6.1 assumes a < -n-2 -- is UNCHECKABLE here: the arXiv\n"
    "    source is not read by this program.  Only the internal consistency of\n"
    "    the reported (5.10) linear term and (6.10) Pochhammer form of A -- both\n"
    "    supplied to this program as input, neither read from the source --\n"
    "    with the paper's own -(x-c) and A is checked (checks 2(v) and 5(iv)).\n"
    " 8. Remark 4.3's step 'Conjecture 2 is invariant under exchanging the two\n"
    "    labels, so Theorem 4.2 proves the conjecture exactly as stated' is a\n"
    "    reading of the source's statement, not a computation; nothing here\n"
    "    tests it.  The theorems as printed prove ONE disjunct, which is\n"
    "    stronger than the conjectures' disjunction, so this gap is benign.\n"
    " 9. The MP admissible range is handled through t = cot phi.  That this\n"
    "    substitution is a bijection from phi in (0,pi) onto t in R, hence that\n"
    "    'phi in (0,phi_c)' really is 't > cot phi_c > 0' and\n"
    "    'phi in (pi-phi_c,pi)' really is 't < -cot phi_c', is a fact about cot\n"
    "    that is used and not checked (check 2(vi) does check that phi_c itself\n"
    "    is well defined and in (0,pi/2))."
)

NOTE_ON_NUMERICS = (
    "NOTE ON NUMERICS: this program uses NO floating point.  The paper's phi is\n"
    "  eliminated in favour of t = cot phi via sin^2 phi = 1/(1+t^2),\n"
    "  cos 2phi = (t^2-1)/(t^2+1), sin 2phi = 2t/(1+t^2), which turns every\n"
    "  Meixner-Pollaczek object into a polynomial with rational coefficients;\n"
    "  phi in (0,phi_c) becomes the exact rational condition (n+1)t^2 > 2 lam\n"
    "  with t > 0, and phi in (pi-phi_c,pi) the same with t < 0 (derived in\n"
    "  check 2 from disc T).  Real zeros are isolated by Sturm sequences over\n"
    "  fractions.Fraction and compared by exact interval refinement, so every\n"
    "  inequality decided below is an exact rational inequality, and no error\n"
    "  bound is needed.  The paper itself contains no numerical experiment and\n"
    "  no computed figure, so nothing printed above is a comparison against a\n"
    "  published numeric run: every number above is computed here."
)


# =========================================================================
# 1. Univariate polynomials over Q, as coefficient lists, low degree first.
# =========================================================================

def utrim(p):
    q = list(p)
    while q and q[-1] == 0:
        q.pop()
    return q


def uadd(p, q):
    n = max(len(p), len(q))
    return utrim([(p[i] if i < len(p) else F(0)) + (q[i] if i < len(q) else F(0))
                  for i in range(n)])


def usub(p, q):
    n = max(len(p), len(q))
    return utrim([(p[i] if i < len(p) else F(0)) - (q[i] if i < len(q) else F(0))
                  for i in range(n)])


def umul(p, q):
    if not p or not q:
        return []
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            if b != 0:
                r[i + j] += a * b
    return utrim(r)


def uscal(p, s):
    s = F(s)
    return [] if s == 0 else utrim([s * c for c in p])


def udeg(p):
    return len(utrim(p)) - 1          # -1 for the zero polynomial


def ueval(p, x):
    v = F(0)
    for c in reversed(p):
        v = v * x + c
    return v


def uderiv(p):
    return utrim([F(i) * p[i] for i in range(1, len(p))])


def umonic(p):
    p = utrim(p)
    return [] if not p else uscal(p, F(1) / p[-1])


def udivmod(p, d):
    """Exact quotient and remainder of p by d over Q."""
    p, d = utrim(p), utrim(d)
    if not d:
        raise ZeroDivisionError("polynomial division by zero")
    q = [F(0)] * max(0, len(p) - len(d) + 1)
    r = list(p)
    while utrim(r) and len(utrim(r)) >= len(d):
        r = utrim(r)
        k = len(r) - len(d)
        c = r[-1] / d[-1]
        q[k] = c
        for i, dc in enumerate(d):
            r[k + i] -= c * dc
        r = utrim(r)
    return utrim(q), utrim(r)


def ugcd(p, q):
    a, b = utrim(p), utrim(q)
    while b:
        a, b = b, udivmod(a, b)[1]
    return umonic(a)


def usquarefree(p):
    """p / gcd(p, p'), the squarefree part (monic)."""
    p = utrim(p)
    g = ugcd(p, uderiv(p))
    return umonic(udivmod(p, g)[0]) if udeg(g) > 0 else umonic(p)


def upow_lin(c):
    """the polynomial (x - c)."""
    return [-F(c), F(1)]


# =========================================================================
# 2. Exact real-root counting and isolation by Sturm sequences.
# =========================================================================

def sturm_chain(p):
    """Sturm chain of the squarefree part of p."""
    p = usquarefree(p)
    chain = [p, uderiv(p)]
    while udeg(chain[-1]) > 0:
        r = udivmod(chain[-2], chain[-1])[1]
        if not r:
            break
        chain.append(uscal(r, F(-1)))
    return [c for c in chain if c]


def _variations(chain, x):
    last, v = 0, 0
    for c in chain:
        s = ueval(c, x)
        s = 0 if s == 0 else (1 if s > 0 else -1)
        if s == 0:
            continue
        if last != 0 and s != last:
            v += 1
        last = s
    return v


def sturm_count(chain, lo, hi):
    """Number of distinct real roots in the half-open interval (lo, hi]."""
    return _variations(chain, F(lo)) - _variations(chain, F(hi))


def cauchy_bound(p):
    """1 + max|a_i / a_deg| strictly bounds the modulus of every real root."""
    p = utrim(p)
    lead = p[-1]
    return F(1) + max([abs(c / lead) for c in p[:-1]] or [F(0)])


def isolate_roots(p, max_depth=200):
    """
    Return a list of (lo, hi) with lo < hi, rational, pairwise disjoint, sorted,
    each containing exactly one distinct real root of p in (lo, hi], and with
    p(lo) != 0 != p(hi).  Exact rational roots are returned as (r, r).
    """
    q = usquarefree(p)
    if udeg(q) <= 0:
        return []
    chain = sturm_chain(q)
    M = cauchy_bound(q)
    lo, hi = -M, M
    if ueval(q, lo) == 0 or ueval(q, hi) == 0:
        raise AssertionError("Cauchy bound is a root -- impossible")
    out, work = [], [(lo, hi, sturm_count(chain, lo, hi), 0)]
    while work:
        a, b, k, d = work.pop()
        if k == 0:
            continue
        if k == 1:
            out.append((a, b))
            continue
        if d > max_depth:
            raise AssertionError("root isolation did not separate")
        m = (a + b) / 2
        if ueval(q, m) == 0:
            out.append((m, m))                 # exact rational root
            m1, m2 = m - (b - a) / F(1 << 40), m + (b - a) / F(1 << 40)
            work.append((a, m1, sturm_count(chain, a, m1), d + 1))
            work.append((m2, b, sturm_count(chain, m2, b), d + 1))
        else:
            work.append((a, m, sturm_count(chain, a, m), d + 1))
            work.append((m, b, sturm_count(chain, m, b), d + 1))
    out.sort()
    total = sturm_count(chain, lo, hi)
    if len(out) != total:
        raise AssertionError("isolation lost a root: %d found, %d counted"
                             % (len(out), total))
    return out


class Root(object):
    """One real root of a rational polynomial, held as a shrinking interval."""

    __slots__ = ("p", "chain", "lo", "hi", "tag")

    def __init__(self, p, chain, lo, hi, tag=""):
        self.p, self.chain, self.lo, self.hi, self.tag = p, chain, lo, hi, tag

    def is_exact(self):
        return self.lo == self.hi

    def bisect(self):
        if self.is_exact():
            return
        m = (self.lo + self.hi) / 2
        if ueval(self.p, m) == 0:
            self.lo = self.hi = m
        elif sturm_count(self.chain, self.lo, m) == 1:
            self.hi = m
        else:
            self.lo = m

    def __repr__(self):
        # exact rational endpoints only: this program prints no floats anywhere
        if self.is_exact():
            return "%s=%s" % (self.tag, self.lo)
        return "%s in (%s,%s)" % (self.tag, self.lo, self.hi)


def roots_of(p, tag=""):
    """All distinct real roots of p, as Root objects in increasing order."""
    q = usquarefree(p)
    chain = sturm_chain(q)
    return [Root(q, chain, a, b, tag) for (a, b) in isolate_roots(p)]


def rcmp(r, s, cap=400):
    """
    Exact comparison of two real algebraic numbers given as Root objects.
    Returns -1, 0 or +1.  Returns 0 only when the two are provably equal
    (a shared exact rational root); raises if they cannot be separated,
    which happens only if the two polynomials share an irrational root.
    """
    for _ in range(cap):
        if r.is_exact() and s.is_exact():
            return -1 if r.lo < s.lo else (0 if r.lo == s.lo else 1)
        if r.hi <= s.lo:
            return -1
        if s.hi <= r.lo:
            return 1
        if r.is_exact() and s.lo < r.lo < s.hi:
            s.bisect()
            continue
        if s.is_exact() and r.lo < s.lo < r.hi:
            r.bisect()
            continue
        r.bisect()
        s.bisect()
    if ugcd(r.p, s.p) != [F(1)] and udeg(ugcd(r.p, s.p)) > 0:
        return 0            # a genuinely shared root: report as equal
    raise AssertionError("could not separate %r from %r" % (r, s))


def sort_roots(rs):
    """Insertion sort of Root objects using the exact comparison rcmp."""
    out = []
    for r in rs:
        i = len(out)
        while i > 0 and rcmp(out[i - 1], r) > 0:
            i -= 1
        out.insert(i, r)
    return out


def strictly_interlaces(big, small):
    """
    Decide the paper's relation  BIG < SMALL  (written G \\prec Q there):
    with deg BIG = deg SMALL + 1 = N+1 and zeros w_1<..<w_{N+1}, z_1<..<z_N,
    the relation asks  w_1 < z_1 < w_2 < ... < z_N < w_{N+1}.
    `big` and `small` are lists of Root objects (unsorted is fine).
    Returns (ok, reason).
    """
    if len(big) != len(small) + 1:
        return (False, "wrong zero counts: %d and %d" % (len(big), len(small)))
    W = sort_roots(list(big))
    Z = sort_roots(list(small))
    for i in range(len(Z)):
        if rcmp(W[i], Z[i]) >= 0:
            return (False, "w_%d < z_%d fails" % (i + 1, i + 1))
        if rcmp(Z[i], W[i + 1]) >= 0:
            return (False, "z_%d < w_%d fails" % (i + 1, i + 2))
    return (True, "ok")


def gap_index(rt, ys):
    """
    Index k such that ys[k-1] < rt < ys[k] for the sorted list `ys`
    (1-based on the paper's convention y_k < . < y_{k+1}: returns k).
    Returns 0 if rt < ys[0] and len(ys) if rt > ys[-1]; -1 if rt equals a y.
    """
    k = 0
    for y in ys:
        c = rcmp(rt, y)
        if c == 0:
            return -1
        if c > 0:
            k += 1
        else:
            break
    return k


# =========================================================================
# 3. Multivariate polynomials over Q: dict {exponent tuple: Fraction}.
#    Used for the symbolic identities (B), (D), (E), (F), (G), (H).
# =========================================================================

NV = 3            # number of variables in the ambient ring


def mconst(c):
    c = F(c)
    return {} if c == 0 else {(0,) * NV: c}


def mvar(i):
    e = [0] * NV
    e[i] = 1
    return {tuple(e): F(1)}


def madd(*ps):
    r = {}
    for p in ps:
        for e, c in p.items():
            v = r.get(e, F(0)) + c
            if v == 0:
                r.pop(e, None)
            else:
                r[e] = v
    return r


def mneg(p):
    return {e: -c for e, c in p.items()}


def msub(p, q):
    return madd(p, mneg(q))


def mmul(p, q):
    r = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            v = r.get(e, F(0)) + c1 * c2
            if v == 0:
                r.pop(e, None)
            else:
                r[e] = v
    return r


def mscal(p, s):
    s = F(s)
    return {} if s == 0 else {e: c * s for e, c in p.items()}


def mpow(p, k):
    r = mconst(1)
    for _ in range(k):
        r = mmul(r, p)
    return r


def mis_zero(p):
    return not p


def mstr(p, names=("x", "L", "t"), cap=240):
    """Printable form.  Zero prints as '0', so a passing residual is legible;
    a nonzero residual is truncated so a failure cannot flood the transcript."""
    if not p:
        return "0"
    terms = []
    for e in sorted(p, reverse=True):
        s = str(p[e])
        for nm, k in zip(names, e):
            if k:
                s += "*%s%s" % (nm, "" if k == 1 else "^%d" % k)
        terms.append(s)
    out = " + ".join(terms)
    if len(out) > cap:
        out = out[:cap] + " ... (%d monomials, NONZERO)" % len(p)
    return out


# ---- rational functions as (numerator, denominator) pairs of mpolys ------

def rf(num, den=None):
    return (num, mconst(1) if den is None else den)


def rf_const(c):
    return (mconst(c), mconst(1))


def rf_var(i):
    return (mvar(i), mconst(1))


def rf_add(A, B):
    return (madd(mmul(A[0], B[1]), mmul(B[0], A[1])), mmul(A[1], B[1]))


def rf_sub(A, B):
    return (msub(mmul(A[0], B[1]), mmul(B[0], A[1])), mmul(A[1], B[1]))


def rf_mul(A, B):
    return (mmul(A[0], B[0]), mmul(A[1], B[1]))


def rf_div(A, B):
    if mis_zero(B[0]):
        raise ZeroDivisionError("rational function division by zero")
    return (mmul(A[0], B[1]), mmul(A[1], B[0]))


def rf_residual(A, B):
    """The cross-multiplied residual numerator (zero iff A == B)."""
    return msub(mmul(A[0], B[1]), mmul(B[0], A[1]))


# =========================================================================
# 4. Meixner-Pollaczek family, from eq:MPrec alone, in the variables
#    (x, L, t) with L = lambda and t = cot phi, so 1/(4 sin^2 phi)
#    = (1+t^2)/4.  Nothing but eq:MPrec is used.
# =========================================================================

X, L, T_ = mvar(0), mvar(1), mvar(2)
ONE_PLUS_T2 = madd(mconst(1), mmul(T_, T_))


def mp_monic_sym(N, shift=0):
    """
    [P_0, ..., P_N] for the monic Meixner-Pollaczek family with parameter
    lambda + shift, as polynomials in (x, L, t) built from eq:MPrec:
       P_{m+1} = (x + (m+lam)t)P_m - m(m+2lam-1)(1+t^2)/4 * P_{m-1}.
    """
    Lp = madd(L, mconst(shift))                  # lambda + shift
    P = [mconst(1)]
    prev = mconst(0)                             # P_{-1} = 0
    for m in range(N):
        lin = madd(X, mmul(madd(mconst(m), Lp), T_))
        gam = mscal(mmul(madd(mconst(m), mscal(Lp, 2), mconst(-1)),
                         ONE_PLUS_T2), F(m, 4))
        nxt = msub(mmul(lin, P[m]), mmul(gam, prev))
        prev = P[m]
        P.append(nxt)
    return P


def mp_completion_quadratic(n):
    """
    (x-E_1)(x-E_2) = T(x)/(1-cos 2phi), derived from eq:T by substituting
    cos2phi = (t^2-1)/(t^2+1), sin2phi = 2t/(1+t^2), 1-cos2phi = 2/(1+t^2):
    it comes out as  x^2 + (n+1)t x + L^2 + (n+1)L(1-t^2)/2.
    """
    return madd(mmul(X, X),
                mscal(mmul(T_, X), n + 1),
                mmul(L, L),
                mscal(mmul(L, msub(mconst(1), mmul(T_, T_))), F(n + 1, 2)))


def mp_T_raw(n):
    """
    eq:T itself, cleared of denominators: (1+t^2)*T(x) as an mpoly.
    (1-cos2phi)=2/(1+t^2), sin2phi=2t/(1+t^2), cos2phi=(t^2-1)/(t^2+1).
    """
    c2 = msub(mmul(T_, T_), mconst(1))                     # (t^2-1)
    return madd(mscal(mmul(X, X), 2),
                mscal(mmul(mmul(T_, X), mconst(1)), 2 * (n + 1)),
                mmul(mmul(L, L), ONE_PLUS_T2),
                mneg(mmul(madd(mmul(L, L), mscal(L, n + 1)), c2)))


def check_mp_mixed_identity(ns, log):
    """
    CHECK: eq:MPmixed is an identity in (x, lambda, t) for each n in ns.
        (2L+n)(2L+n+1)/(2(1-cos2phi)) P_n^{(L)}
          = -(x - L cot phi) P_{n+1}^{(L+1)} + (x-E_1)(x-E_2) P_n^{(L+1)}
    Everything on both sides is built from eq:MPrec and eq:T; the residual is
    derived, not assumed.  Also derived: the degree cancellation the paper
    needs (deg of the right-hand side is n+1, not n+2) and the equality
    T(x)/(1-cos2phi) = (x-E_1)(x-E_2).
    """
    results = []
    for n in ns:
        base = mp_monic_sym(n, 0)
        up = mp_monic_sym(n + 1, 1)
        P, G, Q = base[n], up[n + 1], up[n]
        c = mmul(L, T_)                                     # c = lambda cot phi
        quad = mp_completion_quadratic(n)
        # the paper's own claim T = (1-cos2phi)(x-E_1)(x-E_2), i.e.
        # (1+t^2)T/2 == quad
        t_ok = mis_zero(msub(mscal(mp_T_raw(n), F(1, 2)), quad))
        lhs = mmul(mscal(mmul(madd(mscal(L, 2), mconst(n)),
                              madd(mscal(L, 2), mconst(n + 1))),
                         F(1, 4)),
                   mmul(ONE_PLUS_T2, P))
        rhs = madd(mneg(mmul(msub(X, c), G)), mmul(quad, Q))
        resid = msub(lhs, rhs)
        degx = max([e[0] for e in rhs] or [0])
        # CORRECTED: the expected x-degree of the right-hand side was n+1 and the truth is n.
        # LHS = const * (1+t^2) * P with deg_x P = n, so LHS has degree n; a zero residual therefore
        # FORCES deg_x(RHS) = n -- the degree-(n+2) leading terms of (x-c)G and of quad*Q cancel.
        # Requiring n+1 made a correct identity report FAIL. The expectation is now derived from the
        # left-hand side rather than asserted, so it cannot drift from the mathematics again.
        degx_lhs = max([e[0] for e in lhs] or [0])
        log("    n=%d: monomials in LHS=%d RHS=%d, deg_x(RHS)=%d (deg_x(LHS)=%d, must match),"
            " T=(1-cos2phi)(x-E1)(x-E2): %s, residual=%s"
            % (n, len(lhs), len(rhs), degx, degx_lhs,
               "yes" if t_ok else "NO", mstr(resid)))
        results.append(mis_zero(resid) and t_ok and degx == degx_lhs == n)
    return all(results)


def mcoeff_x(p, i):
    """Coefficient of x^i in p, as an mpoly in the remaining variables."""
    return {(0,) + e[1:]: c for e, c in p.items() if e[0] == i}


def msubst_x(p, q):
    """p with the variable x replaced by the mpoly q."""
    dmax = max([e[0] for e in p] or [0])
    out = mconst(0)
    for i in range(dmax + 1):
        out = madd(out, mmul(mcoeff_x(p, i), mpow(q, i)))
    return out


def check_mp_exterior_point(ns, log):
    """
    CHECK: the three MP facts the branch assignment rests on, each derived as
    an identity in (lambda, t) -- t = cot phi:
      (i)  (c-E_1)(c-E_2) = lambda(2lambda+n+1)/(1-cos2phi)   [paper, thm:MP]
      (ii) T'(c) = (2lambda+n+1) sin 2phi                     [paper, thm:MP]
      (iii) disc T > 0  <=>  (n+1)t^2 > 2 lambda, since
            disc = (2lambda+n+1)((n+1)t^2 - 2lambda) and lambda>0.  This is the
            content of [JK, Lemma 5.1] and shows the phi-range is not narrowed.
      (iv) cot^2(phi_c) = 2 lambda/(n+1), i.e. the paper's
           cos 2phi_c = (2lambda-n-1)/(2lambda+n+1) turns into the exact
           rational branch test used by the sweep in check 11.
      (v)  the source's rewriting of its (5.9) into (5.10): the linear term
           -2 sin^2 phi (x - lambda cot phi)/(1 - cos 2phi) equals -(x-c) with
           c = lambda cot phi, because 2 sin^2 phi = 1 - cos 2phi.  This is the
           step by which eq:MPmixed acquires the shape of eq:mixed-general, so
           it is checked here rather than taken on trust.
    """
    ok = True
    c = mmul(L, T_)
    for n in ns:
        quad = mp_completion_quadratic(n)
        # (i)  quad(c) = (c-E1)(c-E2); target lambda(2L+n+1)/(1-cos2phi)
        got_i = msubst_x(quad, c)
        want_i = mscal(mmul(mmul(L, madd(mscal(L, 2), mconst(n + 1))),
                            ONE_PLUS_T2), F(1, 2))
        # (ii) (1+t^2) T'(x) = 4x + 2(n+1)t  [from mp_T_raw]; evaluate at x=c
        Traw = mp_T_raw(n)
        Traw_d = madd(mscal(mcoeff_x(Traw, 1), 1),
                      mscal(mmul(mcoeff_x(Traw, 2), c), 2))
        got_ii = Traw_d                                      # (1+t^2)T'(c)
        want_ii = mscal(mmul(T_, madd(mscal(L, 2), mconst(n + 1))), 2)
        # (iii) disc of the monic quad
        B, C = mcoeff_x(quad, 1), mcoeff_x(quad, 0)
        got_iii = msub(mmul(B, B), mscal(C, 4))
        want_iii = mmul(madd(mscal(L, 2), mconst(n + 1)),
                        msub(mscal(mmul(T_, T_), n + 1), mscal(L, 2)))
        # (iv) cos2phi at t^2 = 2L/(n+1) equals (2L-n-1)/(2L+n+1)
        t2c = rf(mscal(L, 2), mconst(n + 1))
        cos2 = rf_div(rf_sub(t2c, rf_const(1)), rf_add(t2c, rf_const(1)))
        want_iv = rf(msub(mscal(L, 2), mconst(n + 1)),
                     madd(mscal(L, 2), mconst(n + 1)))
        # (v) 2 sin^2 phi = 1 - cos 2phi, and hence
        #     -2 sin^2 phi (x - L t)/(1 - cos 2phi) = -(x - L t).
        sin2 = rf(mconst(1), ONE_PLUS_T2)                    # sin^2 phi
        oneminus = rf_sub(rf_const(1),
                          rf(msub(mmul(T_, T_), mconst(1)), ONE_PLUS_T2))
        got_v = rf_div(rf_mul(rf_const(-2),
                              rf_mul(sin2, rf(msub(X, c)))), oneminus)
        want_v = rf(mneg(msub(X, c)))
        r1, r2, r3 = msub(got_i, want_i), msub(got_ii, want_ii), \
            msub(got_iii, want_iii)
        r4 = rf_residual(cos2, want_iv)
        r5 = rf_residual(rf_mul(rf_const(2), sin2), oneminus)
        r6 = rf_residual(got_v, want_v)
        # (vi) phi_c is WELL DEFINED and lies in (0, pi/2): the paper writes
        #      phi_c = (1/2)arccos((2L-n-1)/(2L+n+1)) and asserts phi_c in
        #      (0,pi/2) without argument.  That needs the argument of arccos to
        #      lie strictly in (-1,1), i.e. (2L+n+1) - (2L-n-1) > 0 and
        #      (2L+n+1) + (2L-n-1) > 0.  Derive both as polynomial identities:
        #      they are 2(n+1) and 4L, positive for n >= 0 and L > 0.
        den_c = madd(mscal(L, 2), mconst(n + 1))
        num_c = msub(mscal(L, 2), mconst(n + 1))
        r7 = msub(msub(den_c, num_c), mconst(2 * (n + 1)))
        r8 = msub(madd(den_c, num_c), mscal(L, 4))
        log("          phi_c well defined: (den-num)-2(n+1) resid=%s ;"
            " (den+num)-4L resid=%s (so |arg|<1 for L>0)"
            % (mstr(r7), mstr(r8)))
        log("    n=%d: (c-E1)(c-E2) resid=%s ; T'(c) resid=%s" %
            (n, mstr(r1), mstr(r2)))
        log("          disc T = %s ; resid vs (2L+n+1)((n+1)t^2-2L) = %s" %
            (mstr(got_iii), mstr(r3)))
        log("          cot^2 phi_c = 2L/%d ; cos2phi_c resid=%s"
            % (n + 1, mstr(r4)))
        log("          2sin^2phi-(1-cos2phi) resid=%s ; (5.9)->(5.10) linear"
            " term resid=%s" % (mstr(r5), mstr(r6)))
        ok = ok and all(mis_zero(r)
                        for r in (r1, r2, r3, r4, r5, r6, r7, r8))
    return ok


# =========================================================================
# 5. Pseudo-Jacobi family, from eq:PJrec / eq:gamma alone.
#    Symbolic version: rational alpha, symbolic b (variables x = var0,
#    b = var1).  Numeric version: rational alpha and b, univariate in x.
# =========================================================================

BS = mvar(1)          # the symbol b


def pj_denominators(m, alpha):
    """The four denominators of eq:PJrec/eq:gamma at index m, as Fractions."""
    a = F(alpha)
    return {"(m+alpha)": a + m,
            "(m+alpha+1)": a + m + 1,
            "(2m+2alpha-1)": 2 * m + 2 * a - 1,
            "(2m+2alpha+1)": 2 * m + 2 * a + 1}


def pj_monic_sym(N, alpha):
    """
    [P_0,...,P_N] for the monic pseudo-Jacobi family, parameter (alpha, b),
    as mpolys in (x, b), from eq:PJrec/eq:gamma alone.
    Returns (polys, bad) where bad lists any vanishing denominator met.
    """
    a = F(alpha)
    bad = []
    P = [mconst(1)]
    prev = mconst(0)
    for m in range(N):
        dd = pj_denominators(m, a)
        for nm in ("(m+alpha)", "(m+alpha+1)"):
            if dd[nm] == 0:
                bad.append("m=%d %s=0" % (m, nm))
        if any(dd[nm] == 0 for nm in ("(m+alpha)", "(m+alpha+1)")):
            return (None, bad)
        lin = madd(mvar(0), mscal(BS, a / (dd["(m+alpha)"] * dd["(m+alpha+1)"])))
        if m == 0:
            gam = mconst(0)                       # gamma_0 = 0 (factor m)
        else:
            for nm in ("(2m+2alpha-1)", "(2m+2alpha+1)"):
                if dd[nm] == 0:
                    bad.append("m=%d %s=0" % (m, nm))
            if bad:
                return (None, bad)
            k = -F(m) * (m + 2 * a) / (dd["(m+alpha)"] ** 2
                                       * dd["(2m+2alpha-1)"]
                                       * dd["(2m+2alpha+1)"])
            gam = mscal(madd(mconst((m + a) ** 2), mmul(BS, BS)), k)
        nxt = msub(mmul(lin, P[m]), mmul(gam, prev))
        prev = P[m]
        P.append(nxt)
    return (P, bad)


def pj_monic_num(N, alpha, b):
    """
    [P_0,...,P_N] as univariate coefficient lists over Q, for rational
    alpha and b, from eq:PJrec/eq:gamma alone.  Also returns [gamma_1..gamma_N].
    """
    a, b = F(alpha), F(b)
    P, gammas = [[F(1)]], []
    prev = []
    for m in range(N):
        dd = pj_denominators(m, a)
        if any(v == 0 for v in dd.values()):
            raise ZeroDivisionError("eq:PJrec denominator vanishes at m=%d,"
                                    " alpha=%s" % (m, a))
        lin = [a * b / (dd["(m+alpha)"] * dd["(m+alpha+1)"]), F(1)]
        if m == 0:
            gam = F(0)
        else:
            gam = -F(m) * (m + 2 * a) * ((m + a) ** 2 + b * b) \
                / (dd["(m+alpha)"] ** 2 * dd["(2m+2alpha-1)"]
                   * dd["(2m+2alpha+1)"])
            gammas.append(gam)
        P.append(usub(umul(lin, P[m]), uscal(prev, gam)))
        prev = P[m]
    return P, gammas


def pj_R_coeffs(n, a, b=None):
    """
    eq:R exactly as displayed.  If b is a Fraction, returns the univariate
    coefficient list [c0, c1, 1]; if b is None, returns the mpoly in (x, b).
    """
    a = F(a)
    u, v = a + n + 1, a + n + 2
    c1den = u * v
    c0den = u * v * (2 * a + 2 * n + 3)
    if b is None:
        c1 = mscal(BS, -F(n + 1) / c1den)
        c0 = madd(mconst(u * v * (2 * a + n + 2) / c0den),
                  mscal(mmul(BS, BS), F(n + 1) / c0den))
        return madd(mmul(mvar(0), mvar(0)), mmul(mvar(0), c1), c0)
    b = F(b)
    return [(u * v * (2 * a + n + 2) + b * b * (n + 1)) / c0den,
            -b * (n + 1) / c1den, F(1)]


def pj_A(n, a, b=None):
    """
    A = 2(2u-n-1)(2u-n)(u^2+b^2)/(u(2u-1)(2u)(2u+1)), u = a+n+1, from the
    paper's proof of thm:PJ.  mpoly in b if b is None, else a Fraction.
    """
    a = F(a)
    u = a + n + 1
    k = F(2) * (2 * u - n - 1) * (2 * u - n) / (u * (2 * u - 1) * (2 * u)
                                               * (2 * u + 1))
    if b is None:
        return mscal(madd(mconst(u * u), mmul(BS, BS)), k)
    return k * (u * u + F(b) ** 2)


def in_pj_range(n, a):
    """
    EXACT test of eq:PJrange: a < (-3n-5-sqrt(n^2+2n+5))/4 and a != -n-2.
    Equivalent to  4a+3n+5 < 0  and  (4a+3n+5)^2 > n^2+2n+5.
    """
    a = F(a)
    if a == -n - 2:
        return False
    z = 4 * a + 3 * n + 5
    return z < 0 and z * z > pj_radicand(n)


def brange_bound_sq(n, a):
    """
    The square of the eq:brange right-hand side, as an exact Fraction:
      4(a+n+1)^2(a+n+2)^2(-2a-n-2)
        / ((n+1)|4(a+n+1)^2-(2n-2)(a+n+1)-(n+1)|).
    """
    a = F(a)
    u, v = a + n + 1, a + n + 2
    Dq = 4 * u * u - (2 * n - 2) * u - (n + 1)
    num = 4 * u * u * v * v * (-2 * a - n - 2)
    return num / (F(n + 1) * abs(Dq)), Dq


def in_brange(n, a, b):
    """EXACT test of eq:brange: b^2 > (that bound)."""
    bb, _ = brange_bound_sq(n, a)
    return F(b) ** 2 > bb


def check_pj_mixed_identity(cases, log):
    """
    CHECK: eq:PJmixed  A P_n(x;a,b) = -(x - b/u)P_{n+1}(x;a+1,b)
                                       + R_{n,a,b}(x) P_n(x;a+1,b)
    as an identity in (x, b) at each (n, a) in `cases`, with a rational and
    IN RANGE for eq:PJrange.  Both sides are built from eq:PJrec/eq:gamma and
    eq:R; the residual is derived.  Also derived: deg_x of the right-hand side
    is n, i.e. BOTH the x^{n+2} and the x^{n+1} terms cancel, which is the
    degree collapse the paper invokes to force A constant.
    """
    ok = True
    for (n, a) in cases:
        a = F(a)
        if not in_pj_range(n, a):
            log("    n=%d a=%s: SKIPPED, not in eq:PJrange" % (n, a))
            ok = False
            continue
        u = a + n + 1
        Pa, bad1 = pj_monic_sym(n, a)
        Pb, bad2 = pj_monic_sym(n + 1, a + 1)
        if Pa is None or Pb is None:
            log("    n=%d a=%s: FAIL, vanishing denominator %s %s"
                % (n, a, bad1, bad2))
            ok = False
            continue
        P, G, Q = Pa[n], Pb[n + 1], Pb[n]
        c = mscal(BS, F(1) / u)                      # c = b/u
        lhs = mmul(pj_A(n, a), P)
        rhs = madd(mneg(mmul(msub(mvar(0), c), G)), mmul(pj_R_coeffs(n, a), Q))
        resid = msub(lhs, rhs)
        degx = max([e[0] for e in rhs] or [0])
        log("    n=%d a=%s (u=%s): deg_x(RHS)=%d (need %d), monomials"
            " LHS=%d RHS=%d, residual=%s"
            % (n, a, u, degx, n, len(lhs), len(rhs),
               mstr(resid, ("x", "b", "?"))))
        ok = ok and mis_zero(resid) and degx == n
    return ok


def check_pj_range_is_needed(log):
    """
    CHECK, as a positive derivation rather than a caveat, that eq:PJrange is
    exactly what eq:PJrec needs: at a = -9/2, n = 3 the parameter a+1 = -7/2
    makes the eq:gamma
    denominator factor (2m+2alpha+1) vanish at m = 3, so eq:PJmixed is
    undefined there -- and a = -9/2 is NOT in eq:PJrange for n = 3.  We check
    both halves, and that a = -9/2 IS in range for n = 2 (where it is usable).
    """
    dd = pj_denominators(3, F(-7, 2))
    vanish = [k for k, v in dd.items() if v == 0]
    polys, bad = pj_monic_sym(4, F(-7, 2))
    r1 = (vanish == ["(2m+2alpha+1)"]) and polys is None and bad
    r2 = not in_pj_range(3, F(-9, 2))
    r3 = in_pj_range(2, F(-9, 2))
    log("    alpha=-7/2, m=3: vanishing denominators %s ; builder refused: %s"
        % (vanish, bad))
    log("    a=-9/2 in eq:PJrange for n=3? %s (bound = (-14-sqrt(20))/4)"
        % ("no" if r2 else "yes"))
    log("    a=-9/2 in eq:PJrange for n=2? %s (bound = (-11-sqrt(13))/4)"
        % ("yes" if r3 else "no"))
    return bool(r1) and r2 and r3


def check_pj_symbolic_nab(log):
    """
    CHECK, symbolically in the three free variables (n, a, b) at once -- so for
    ALL n, not a sample -- the four rational-function claims of section 4:
      (i)   R_{n,a,b}(b/u) = (2u-n)(u^2+b^2)/(u^2(2u+1))           [eq:Rc]
      (ii)  R'_{n,a,b}(b/u) = b(2a+n+3)/((a+n+1)(a+n+2)) = Delta_c [thm:PJ]
      (iii) disc R_{n,a,b} = (4u^2v^2(n-2u) - b^2(n+1)D)/(u^2v^2(2u+1)) with
            D = 4u^2-(2n-2)u-(n+1) the quantity inside eq:brange's |...|,
            and n-2u = -2a-n-2 the numerator factor of eq:brange.  With
            2u+1<0 and D>0 this says disc>0 <=> eq:brange exactly.
      (iv)  the Pochhammer-form expression for A of header item 8 equals the
            paper's A.  This is a consistency check between two expressions
            supplied to the program; it says nothing about the source, which
            is not read here (GAPS item 7).
    Variable order here is (n, a, b).
    """
    nn, aa, bb = rf_var(0), rf_var(1), rf_var(2)
    one = rf_const(1)
    u = rf_add(rf_add(aa, nn), one)
    v = rf_add(u, one)
    n1 = rf_add(nn, one)                                   # n+1
    p2 = rf_add(rf_mul(rf_const(2), aa), rf_add(nn, rf_const(2)))   # 2a+n+2
    p3 = rf_add(rf_mul(rf_const(2), aa),
                rf_add(rf_mul(rf_const(2), nn), rf_const(3)))       # 2a+2n+3
    c = rf_div(bb, u)
    c1 = rf_div(rf_mul(rf_const(-1), rf_mul(bb, n1)), rf_mul(u, v))
    c0 = rf_div(rf_add(rf_mul(rf_mul(u, v), p2), rf_mul(rf_mul(bb, bb), n1)),
                rf_mul(rf_mul(u, v), p3))
    twou = rf_mul(rf_const(2), u)
    # (i)
    Rc = rf_add(rf_add(rf_mul(c, c), rf_mul(c1, c)), c0)
    want_i = rf_div(rf_mul(rf_sub(twou, nn),
                           rf_add(rf_mul(u, u), rf_mul(bb, bb))),
                    rf_mul(rf_mul(u, u), rf_add(twou, one)))
    # (ii)
    Rpc = rf_add(rf_mul(rf_const(2), c), c1)
    want_ii = rf_div(rf_mul(bb, rf_add(rf_mul(rf_const(2), aa),
                                       rf_add(nn, rf_const(3)))),
                     rf_mul(u, v))
    # (iii)
    disc = rf_sub(rf_mul(c1, c1), rf_mul(rf_const(4), c0))
    D = rf_sub(rf_sub(rf_mul(rf_const(4), rf_mul(u, u)),
                      rf_mul(rf_sub(rf_mul(rf_const(2), nn), rf_const(2)), u)),
               n1)
    u2v2 = rf_mul(rf_mul(u, u), rf_mul(v, v))
    want_iii = rf_div(rf_sub(rf_mul(rf_const(4), rf_mul(u2v2,
                                                       rf_sub(nn, twou))),
                             rf_mul(rf_mul(bb, bb), rf_mul(n1, D))),
                      rf_mul(u2v2, rf_add(twou, one)))
    # eq:brange's numerator factor -2a-n-2 equals n-2u
    fac = rf_sub(rf_mul(rf_const(-2), aa), rf_add(nn, rf_const(2)))
    # (iv) Pochhammer form of A, written out from the LITERAL a,n expressions
    #      (2a+n+1)_2 = (2a+n+1)(2a+n+2),
    #      (2a+2n+1)_3 = (2a+2n+1)(2a+2n+2)(2a+2n+3),  denominator (n+a+1).
    def lin(ca, cn, k):
        return rf_add(rf_add(rf_mul(rf_const(ca), aa), rf_mul(rf_const(cn), nn)),
                      rf_const(k))
    poch2 = rf_mul(lin(2, 1, 1), lin(2, 1, 2))
    poch3 = rf_mul(lin(2, 2, 1), rf_mul(lin(2, 2, 2), lin(2, 2, 3)))
    Apoch = rf_div(rf_mul(rf_const(2),
                          rf_mul(poch2,
                                 rf_add(rf_mul(lin(1, 1, 1), lin(1, 1, 1)),
                                        rf_mul(bb, bb)))),
                   rf_mul(lin(1, 1, 1), poch3))
    Apaper = rf_div(rf_mul(rf_const(2),
                           rf_mul(rf_mul(rf_sub(twou, rf_add(nn, one)),
                                         rf_sub(twou, nn)),
                                  rf_add(rf_mul(u, u), rf_mul(bb, bb)))),
                    rf_mul(u, rf_mul(rf_sub(twou, one),
                                     rf_mul(twou, rf_add(twou, one)))))
    names = ("n", "a", "b")
    tests = [("R(b/u) = (2u-n)(u^2+b^2)/(u^2(2u+1))", Rc, want_i),
             ("R'(b/u) = Delta_c", Rpc, want_ii),
             ("disc R = (4u^2v^2(n-2u) - b^2(n+1)D)/(u^2v^2(2u+1))",
              disc, want_iii),
             ("eq:brange factor -2a-n-2 = n-2u", fac, rf_sub(nn, twou)),
             ("Pochhammer form of A = paper's A", Apoch, Apaper)]
    ok = True
    for (nm, got, want) in tests:
        r = rf_residual(got, want)
        log("    %-52s residual=%s" % (nm, mstr(r, names)))
        ok = ok and mis_zero(r)
    return ok


def _reduce_s(p):
    """Reduce an mpoly in (n, s, -) modulo s^2 = n^2+2n+5."""
    s2 = madd(mmul(mvar(0), mvar(0)), mscal(mvar(0), 2), mconst(5))
    changed = True
    while changed:
        changed = False
        for e in list(p):
            if e[1] >= 2:
                c = p.pop(e)
                lower = (e[0], e[1] - 2, e[2])
                p = madd(p, mscal(mmul({lower: c}, s2), 1))
                changed = True
                break
    return p


def check_range_arithmetic(nmax, log):
    """
    CHECK the paper's arithmetic claims about eq:PJrange, exactly.
    Writing s = sqrt(n^2+2n+5) and bound = (-3n-5-s)/4:
      (a) the residual band -n-2 < a < bound is nonempty for every n>=1,
          because (n+3)^2 - s^2 = 4n+4 > 0 (the paper's n+3 > s);
      (b) bound < -n-3/2, because s^2 - (n+1)^2 = 4 > 0; hence u = a+n+1 < -1/2
          and a+1 < -n-1/2, which is what Lemma 4.1 and A>0 need;
      (c) eq:PJrange is EXACTLY D(u) > 0 for D(u) = 4u^2-(2n-2)u-(n+1), the
          quantity inside eq:brange's absolute value: u_- = bound+n+1 =
          (n-1-s)/4 is derived to be a root of D, D opens upward, and
          u_+ - u_- = s/2 > 0.
    """
    nn, ss = mvar(0), mvar(1)
    id_a = _reduce_s(msub(mmul(madd(nn, mconst(3)), madd(nn, mconst(3))),
                          mmul(ss, ss)))
    id_b = _reduce_s(msub(mmul(ss, ss),
                          mmul(madd(nn, mconst(1)), madd(nn, mconst(1)))))
    # u_- = bound + n + 1 = (n-1-s)/4
    bound = mscal(madd(mscal(nn, -3), mconst(-5), mneg(ss)), F(1, 4))
    um = madd(bound, nn, mconst(1))
    id_c1 = msub(um, mscal(madd(nn, mconst(-1), mneg(ss)), F(1, 4)))
    Du = msub(msub(mscal(mmul(um, um), 4),
                   mmul(msub(mscal(nn, 2), mconst(2)), um)),
              madd(nn, mconst(1)))
    id_c2 = _reduce_s(Du)
    log("    (a) (n+3)^2 - s^2 = %s   (must be 4n+4, positive for n>=1)"
        % mstr(id_a, ("n", "s", "?")))
    log("    (b) s^2 - (n+1)^2 = %s   (must be 4, positive)"
        % mstr(id_b, ("n", "s", "?")))
    log("    (c) (bound+n+1) - (n-1-s)/4 = %s ; D(u_-) mod s^2 = %s"
        % (mstr(id_c1, ("n", "s", "?")), mstr(id_c2, ("n", "s", "?"))))
    ok = (id_a == {(1, 0, 0): F(4), (0, 0, 0): F(4)}
          and id_b == {(0, 0, 0): F(4)}
          and mis_zero(id_c1) and mis_zero(id_c2))
    # The same two facts as plain integer inequalities.  The loop lower bound is
    # NOT taken from the paper: we scan from well below the paper's claimed
    # threshold, DERIVE the smallest n at which both inequalities hold, and only
    # then compare that derived threshold with the paper's "for every n >= 1".
    hold = []
    for n in range(-10, nmax + 1):
        s2 = pj_radicand(n)
        if (n + 3) ** 2 > s2 and s2 > (n + 1) ** 2:
            hold.append(n)
    nmin_derived = hold[0] if hold else None
    gaps = [n for n in range(-10, nmax + 1)
            if n >= nmin_derived and n not in hold] if hold else []
    log("    integer form: both inequalities hold on n = %s..%d with no gaps"
        " above the derived threshold (%s); derived smallest n = %s, paper"
        " claims n >= %d" % (nmin_derived, nmax, "yes" if not gaps else gaps,
                             nmin_derived, DECLARED_FIGURES["min_n_for_band"]))
    if nmin_derived is None or gaps or nmin_derived > \
            DECLARED_FIGURES["min_n_for_band"]:
        log("    FAIL: derived threshold %s does not support the paper's"
            " n >= %d" % (nmin_derived, DECLARED_FIGURES["min_n_for_band"]))
        ok = False
    # The general-n statement is NOT the loop: it is the symbolic identity
    # (n+3)^2 - s^2 = 4n+4 checked above, which is positive for every integer
    # n >= 0 and so for every n >= 1.  The loop is a redundant confirmation.
    log("    general-n conclusion rests on the symbolic identity"
        " (n+3)^2-s^2 = 4n+4 > 0 for n >= 0, not on the finite loop")
    return ok


def check_gamma_positive(log):
    """
    CHECK the paper's claim in the proof of Lemma 4.1: for 1 <= m <= n and
    alpha < -n-1/2 (alpha != -n-1), each of the three factors m+2alpha,
    2m+2alpha-1, 2m+2alpha+1 is negative and hence gamma_m > 0.  Derived on an
    exact rational grid: the gamma_m are recomputed from eq:gamma and their
    signs are read off, not assumed.
    """
    ok, ncheck, ninrange = True, 0, 0
    alphas_extra = [F(-1, 2), F(-3, 100)]
    for n in range(1, 7):
        cands = [F(-2 * n - 1, 2) - F(1, 100), F(-n) - F(51, 50),
                 F(-n - 2), F(-2 * n - 3, 2), F(-3 * n - 7)]
        for al in cands + [F(-n) - F(1, 2) + e for e in alphas_extra]:
            if al == -n - 1:
                continue
            inrange = al < F(-2 * n - 1, 2)
            for b in (F(0), F(1, 3), F(7), F(-5)):
                try:
                    _, gam = pj_monic_num(n + 1, al, b)
                except ZeroDivisionError:
                    continue
                pos = all(g > 0 for g in gam)
                facs = all((m + 2 * al) < 0 and (2 * m + 2 * al - 1) < 0
                           and (2 * m + 2 * al + 1) < 0
                           for m in range(1, n + 1))
                ncheck += 1
                if inrange:
                    ninrange += 1
                if inrange and not (pos and facs):
                    log("    FAIL n=%d alpha=%s b=%s gammas=%s"
                        % (n, al, b, gam))
                    ok = False
                if inrange and len(gam) != n:
                    log("    FAIL wrong gamma count n=%d" % n)
                    ok = False
    # COVERAGE: the claim is universally quantified over alpha < -n-1/2, so only
    # the in-range points test anything.  A grid on which no point were in range
    # would pass the loop above vacuously; require and print the in-range count.
    log("    %d (n,alpha,b) grid points, of which %d actually satisfy"
        " alpha < -n-1/2 (the others test nothing); gamma_m>0 and all three"
        " factors negative on all %d of those" % (ncheck, ninrange, ninrange))
    if ninrange < 24:
        log("    FAIL: too few in-range points to be a test (%d)" % ninrange)
        ok = False
    return ok


def check_lemma41_interlacing(log):
    """
    CHECK Lemma 4.1 itself: P_{n+1}(x;alpha,b) < P_n(x;alpha,b) for
    n>=1, alpha < -n-1/2, alpha != -n-1.  Both polynomials are built from
    eq:PJrec/eq:gamma and their zeros are isolated exactly; the interlacing is
    then decided exactly.  Also derived: each has exactly its degree many
    distinct real zeros (which the paper obtains from the tridiagonal model).
    """
    ok, npts = True, 0
    for n in (1, 2, 3, 4, 5):
        for al in (F(-n) - F(51, 50), F(-2 * n - 3, 2), F(-3 * n - 7),
                   F(-n - 2) if -n - 2 < F(-2 * n - 1, 2) else F(-3 * n - 4)):
            if al >= F(-2 * n - 1, 2) or al == -n - 1:
                continue
            for b in (F(0), F(2, 5), F(-3), F(11, 2)):
                Ps, _ = pj_monic_num(n + 1, al, b)
                P, Pn1 = Ps[n], Ps[n + 1]
                rp, rq = roots_of(Pn1, "y"), roots_of(P, "z")
                npts += 1
                if len(rp) != n + 1 or len(rq) != n:
                    log("    FAIL zero counts n=%d alpha=%s b=%s: %d, %d"
                        % (n, al, b, len(rp), len(rq)))
                    ok = False
                    continue
                good, why = strictly_interlaces(rp, rq)
                if not good:
                    log("    FAIL n=%d alpha=%s b=%s: %s" % (n, al, b, why))
                    ok = False
    log("    %d (n,alpha,b) points: all zeros real and simple in the right"
        " counts, and P_{n+1} < P_n in every case" % npts)
    if npts < 60:
        log("    FAIL: too few admissible points to be a test (%d)" % npts)
        ok = False
    return ok


# =========================================================================
# 6. Lemma 2.1 (the exterior-point criterion) on rigid instances.
# =========================================================================

def build_criterion_instance(rng, n):
    """
    Build a rigid instance of eq:mixed-general.  We choose rational zeros
    y_1<..<y_{n+1} for G, interlacing rational zeros q_i in (y_i,y_{i+1}) for Q
    (so G < Q by construction), and r_- , r_+ in two DISTINCT interior gaps of
    G.  The point c is then NOT free: the degree-(n+1) coefficient of
    -(x-c)G + (x-r_-)(x-r_+)Q must vanish for A P to have degree n, which
    forces c = -sum(y) + sum(q) + r_- + r_+.  A is then read off as the leading
    coefficient of that right-hand side and P as the quotient.
    Because c is forced, whether Lemma 2.1's exterior hypothesis
    (c-r_-)(c-r_+) > 0 can hold at all is itself a derived fact.  Writing
    K = sum(y) - sum(q) = y_{n+1} - sum_i (q_i - y_i), which lies in
    (y_1, y_{n+1}), the forced c satisfies
        c - r_- = r_+ - K      and      c - r_+ = r_- - K,
    so the hypothesis is exactly K outside [r_-, r_+], with
        c < r_-  <=>  K > r_+        and        c > r_+  <=>  K < r_- .
    We therefore compute K first and then place r_-, r_+ on one side of it, so
    that the instances actually exercise the lemma instead of being rejected.
    Returns a dict, or None if no admissible placement exists.
    """
    ys = sorted(rng.sample(range(-9, 10), n + 1))
    ys = [F(y) for y in ys]
    qs = []
    for i in range(n):
        w = ys[i + 1] - ys[i]
        qs.append(ys[i] + w * F(rng.randint(1, 6), 7))
    K = sum(ys) - sum(qs)
    if K in ys or K in qs:
        return None
    jK = max(j for j in range(n) if ys[j] < K)      # gap of G containing K
    opts = []
    for j2 in range(0, jK + 1):                     # case c < r_- : K > r_+
        lo, hi = ys[j2], min(ys[j2 + 1], K)
        if hi > lo:
            opts += [("left", j1, j2, lo, hi) for j1 in range(0, j2)]
    for j1 in range(jK, n):                         # case c > r_+ : K < r_-
        lo, hi = max(ys[j1], K), ys[j1 + 1]
        if hi > lo:
            opts += [("right", j1, j2, lo, hi) for j2 in range(j1 + 1, n)]
    if not opts:
        return None
    case, j1, j2, clo, chi = rng.choice(opts)

    def pick(lo, hi):
        return lo + (hi - lo) * F(rng.randint(1, 6), 7)

    if case == "left":
        rp = pick(clo, chi)
        rm = pick(ys[j1], ys[j1 + 1])
    else:
        rm = pick(clo, chi)
        rp = pick(ys[j2], ys[j2 + 1])
    if rm >= rp or rm in qs or rp in qs or rm == K or rp == K:
        return None
    c = -sum(ys) + sum(qs) + rm + rp
    if (c - rm) * (c - rp) <= 0:
        return None                        # placement rule failed: reject
    G = [F(1)]
    for y in ys:
        G = umul(G, upow_lin(y))
    Q = [F(1)]
    for q in qs:
        Q = umul(Q, upow_lin(q))
    R = umul(upow_lin(rm), upow_lin(rp))
    rhs = uadd(uscal(umul(upow_lin(c), G), F(-1)), umul(R, Q))
    if udeg(rhs) != n:
        return None                        # the cancellation failed: reject
    A = rhs[-1]
    P = uscal(rhs, F(1) / A)
    return {"n": n, "ys": ys, "qs": qs, "rm": rm, "rp": rp, "c": c, "K": K,
            "G": G, "Q": Q, "R": R, "A": A, "P": P, "case": case,
            "kgap": j1 if c < rm else j2}


def sg(v):
    return 0 if v == 0 else (1 if v > 0 else -1)


def rational_as_root(q, tag=""):
    p = upow_lin(q)
    return Root(p, sturm_chain(p), F(q), F(q), tag)


def check_lemma_criterion(ninst, log):
    """
    CHECK Lemma 2.1 on rigid instances, and every intermediate claim of its
    proof.  Derived per instance: A is a positive CONSTANT (the degree
    cancellation); R(y_k) > 0 > R(y_{k+1}) (resp. reversed); the sign chain
    sgn Q(y_k) = sgn G(rho) = -sgn Q(y_{k+1}) and
    sgn P(y_k) = sgn P(y_{k+1}) = eps sgn G(rho), sgn P(rho) = -eps sgn G(rho),
    where rho is whichever of r_-, r_+ faces c and eps = +1 if c<r_-, -1 if
    c>r_+; that F = R P is squarefree of degree n+2 with exactly 3 zeros in
    rho's gap and exactly 1 in each of the other n-1 interior gaps; and the
    lemma's interlacing conclusion, decided by exact root isolation.
    """
    rng = random.Random(20260822)
    built = met = 0
    fails = []
    cases = {}
    signs_ok = counts_ok = 0
    # The old form drew a FIXED number of attempts and passed on `met > 0`, so a
    # single instance that happened to satisfy A>0 would have been enough to
    # declare Lemma 2.1 verified.  Draw until a real sample of HYPOTHESIS-MET
    # instances has accumulated, and require a floor on that sample below.
    target, attempts, cap = ninst, 0, 10 * ninst
    while met < target and attempts < cap:
        attempts += 1
        n = rng.choice([2, 2, 3, 3, 4, 5])
        ins = build_criterion_instance(rng, n)
        if ins is None:
            continue
        built += 1
        A, c, rm, rp, K = ins["A"], ins["c"], ins["rm"], ins["rp"], ins["K"]
        if not (c - rm == rp - K and c - rp == rm - K):
            fails.append("forced-c identity c-r-=r+-K failed")
        if udeg(ins["P"]) != n or ins["P"][-1] != 1:
            fails.append("P not monic of degree n")
        if A <= 0 or (c - rm) * (c - rp) <= 0:
            continue                       # Lemma 2.1's hypotheses not met
        met += 1
        cases[ins["case"]] = cases.get(ins["case"], 0) + 1
        ys, G, Q, P, R = ins["ys"], ins["G"], ins["Q"], ins["P"], ins["R"]
        k = ins["kgap"]
        rho, eps = (rm, 1) if c < rm else (rp, -1)
        sG = sg(ueval(G, rho))
        want_R = (1, -1) if eps == 1 else (-1, 1)
        good = (sg(ueval(R, ys[k])), sg(ueval(R, ys[k + 1]))) == want_R
        good = good and sg(ueval(Q, ys[k])) == sG
        good = good and sg(ueval(Q, ys[k + 1])) == -sG
        good = good and sg(ueval(P, ys[k])) == eps * sG
        good = good and sg(ueval(P, ys[k + 1])) == eps * sG
        good = good and sg(ueval(P, rho)) == -eps * sG
        if good:
            signs_ok += 1
        else:
            fails.append("sign chain, n=%d c=%s" % (n, c))
        Fp = umul(R, P)
        chF = sturm_chain(Fp)
        per = [sturm_count(chF, ys[j], ys[j + 1]) for j in range(n)]
        exp = [3 if j == k else 1 for j in range(n)]
        sqf = udeg(ugcd(Fp, uderiv(Fp))) <= 0
        if per == exp and sum(per) == n + 2 and udeg(Fp) == n + 2 and sqf:
            counts_ok += 1
        else:
            fails.append("F zero count %s != %s, n=%d" % (per, exp, n))
        rr = roots_of(P, "p")
        if eps == 1:
            big = [rational_as_root(rm, "r-")] + [rational_as_root(y, "y")
                                                  for y in ys]
            small = [rational_as_root(rp, "r+")] + rr
        else:
            big = [rational_as_root(rp, "r+")] + [rational_as_root(y, "y")
                                                  for y in ys]
            small = [rational_as_root(rm, "r-")] + rr
        okc, why = strictly_interlaces(big, small)
        if not okc:
            fails.append("conclusion n=%d c=%s: %s" % (n, c, why))
    log("    %d attempts, %d instances built, %d met Lemma 2.1's hypotheses"
        " (A>0 constant and c outside [r-,r+]); branch coverage c<r-: %d,"
        " c>r+: %d" % (attempts, built, met, cases.get("left", 0),
                       cases.get("right", 0)))
    log("    sign chain verified on %d/%d ; F-zero counts (3 in rho's gap, 1"
        " elsewhere, n+2 total, simple) verified on %d/%d"
        % (signs_ok, met, counts_ok, met))
    log("    interlacing conclusion failures: %d" % len(
        [f for f in fails if f.startswith("conclusion")]))
    for f in fails[:6]:
        log("      FAIL: %s" % f)
    log("    SCOPE: Lemma 2.1 is tested on random RIGID RATIONAL instances"
        " (integer zeros for G, rational interlacing zeros for Q, c forced by"
        " the degree cancellation).  It is NOT proved here for general real"
        " data, and the actual Meixner-Pollaczek / pseudo-Jacobi instances have"
        " irrational zeros; those are covered separately, and only pointwise,"
        " by checks 11 and 12.")
    if met < 120 or cases.get("left", 0) < 10 or cases.get("right", 0) < 10:
        log("    FAIL: sample too small or too skewed to be a test"
            " (met=%d, left=%d, right=%d)"
            % (met, cases.get("left", 0), cases.get("right", 0)))
        return False
    return (not fails and signs_ok == met and counts_ok == met)


# =========================================================================
# 7. End-to-end sweeps, in exact rational arithmetic.
# =========================================================================

def rat_sqrt_lower(fr, digits=25):
    """
    A rational LOWER bound for sqrt(fr) with error < 10^-digits.
    Used only to CHOOSE sample points; every subsequent comparison is exact.
    """
    fr = F(fr)
    p, q = fr.numerator, fr.denominator
    if p < 0:
        raise ValueError("negative radicand")
    scale = 10 ** digits
    return F(math.isqrt(p * q * scale * scale), q * scale)


def mp_monic_num(N, lam, t):
    """
    [P_0,...,P_N] for the monic Meixner-Pollaczek family with parameter lam,
    at cot phi = t, as univariate coefficient lists over Q, from eq:MPrec with
    1/(4 sin^2 phi) = (1+t^2)/4.
    """
    lam, t = F(lam), F(t)
    fac = (1 + t * t) / 4
    P, prev = [[F(1)]], []
    for m in range(N):
        lin = [(m + lam) * t, F(1)]
        gam = F(m) * (m + 2 * lam - 1) * fac
        P.append(usub(umul(lin, P[m]), uscal(prev, gam)))
        prev = P[m]
    return P


def mp_quad_num(n, lam, t):
    """
    (x-E_1)(x-E_2) = T(x)/(1-cos2phi) at cot phi = t, as a coefficient list.
    Derived in check 1 to be x^2 + (n+1)t x + lam^2 + (n+1)lam(1-t^2)/2.
    """
    lam, t = F(lam), F(t)
    return [lam * lam + F(n + 1) * lam * (1 - t * t) / 2, F(n + 1) * t, F(1)]


def mspecialize_LT(p, lam, t):
    """
    Specialise an mpoly in (x, L, t) at L = lam, t = t, returning a univariate
    coefficient list in x.  Used to cross-check the closed-form numeric helpers
    against the symbolic objects that check 1 validated against eq:T / eq:MPrec,
    so a divergence between the two code paths cannot go unnoticed.
    """
    lam, t = F(lam), F(t)
    out = {}
    for e, cf in p.items():
        out[e[0]] = out.get(e[0], F(0)) + cf * lam ** e[1] * t ** e[2]
    d = max(out) if out else -1
    return utrim([out.get(i, F(0)) for i in range(d + 1)])


def zero_profile(zs, ys):
    """
    Distribute the zeros `zs` over the interior gaps of the sorted zero list
    `ys`: returns {gap index k (1-based): count}, omitting empty gaps.  A zero
    that coincides with some y, or falls outside (ys[0], ys[-1]), shows up under
    the key -1 resp. 0 / len(ys) and so cannot be mistaken for an interior one.
    """
    prof = {}
    for z in zs:
        k = gap_index(z, ys)
        prof[k] = prof.get(k, 0) + 1
    return prof


def sweep_theorem31(log):
    """
    CHECK Theorem 3.1 end to end on a 72-point grid: n in {2,3,4},
    lambda in {3/10, 1, 5/2}, cot phi = mu * cot phi_c for
    mu in {51/50, 6/5, 2, 5} (mu>1 <=> phi<phi_c), on BOTH branches
    (t>0 is phi in (0,phi_c); t<0 is phi in (pi-phi_c,pi)).
    At each point everything is derived: the admissibility (n+1)t^2 > 2 lambda,
    disc T > 0, the exterior-point value (c-E_1)(c-E_2) = quad(c) and its
    equality with lambda(2lambda+n+1)/(1-cos2phi) and its positivity,
    eq:MPmixed as an exact polynomial identity at this point, the interlacing
    G < Q that the paper takes from orthogonality, the zeros of
    G = P_{n+1}^{(lam+1)}, Q = P_n^{(lam+1)} and P = P_n^{(lam)}, the
    completion points E_1<E_2, whether they lie in DISTINCT INTERIOR gaps of G
    (Theorem 3.1's only hypothesis), the side of [E_1,E_2] on which
    c = lam cot phi falls versus the sign of T'(c) (= sign of t), and finally
    the asserted interlacing.
    """
    npts = hyp = fails = 0
    branchbad = 0
    # per-branch hypothesis coverage: Theorem 3.1 has TWO conclusions, one for
    # phi in (0,phi_c) (t>0) and one for phi in (pi-phi_c,pi) (t<0).  A pass
    # condition of "hyp > 0" would let one of them go entirely untested, so
    # count and require both.
    hyp_by_branch = {1: 0, -1: 0}
    placement_ok = placement_bad = 0
    ns_seen, lams_seen = set(), set()
    NS = (2, 3, 4)
    LAMS = (F(3, 10), F(1), F(5, 2))
    MUS = (F(51, 50), F(6, 5), F(2), F(5))
    SGNS = (1, -1)
    for n in NS:
        for lam in LAMS:
            ns_seen.add(n)
            lams_seen.add(lam)
            tc2 = F(2 * lam, n + 1)
            tc = rat_sqrt_lower(tc2)
            for mu in MUS:
                # snap to a small-denominator rational: sample choice only,
                # admissibility is then re-tested exactly on the next line
                tabs = snap_up(mu * tc, 10 ** 3)
                for sgn in SGNS:
                    t = sgn * tabs
                    npts += 1
                    assert F(n + 1) * t * t > 2 * lam, "sample not admissible"
                    disc = (2 * lam + n + 1) * (F(n + 1) * t * t - 2 * lam)
                    assert disc > 0
                    Pl = mp_monic_num(n, lam, t)[n]
                    Pu = mp_monic_num(n + 1, lam + 1, t)
                    G, Q = Pu[n + 1], Pu[n]
                    quad = mp_quad_num(n, lam, t)
                    # mp_quad_num is a hand-simplified closed form; check it
                    # against the SYMBOLIC completion quadratic that check 1
                    # verified equals eq:T/(1-cos2phi), so the two code paths
                    # cannot drift apart silently.
                    quad_sym = mspecialize_LT(mp_completion_quadratic(n),
                                              lam, t)
                    # the exterior-point value, exactly, and eq:MPmixed here
                    want = lam * (2 * lam + n + 1) * (1 + t * t) / 2
                    kcoef = (2 * lam + n) * (2 * lam + n + 1) * (1 + t * t) / 4
                    mixed = usub(uscal(Pl, kcoef),
                                 uadd(uscal(umul(upow_lin(lam * t), G), F(-1)),
                                      umul(quad, Q)))
                    # kcoef is the paper's A for the MP case; Lemma 2.1 needs
                    # A > 0, which the paper asserts without proof ("the left
                    # coefficient is positive").  Decide it here.
                    if quad != quad_sym or kcoef <= 0:
                        log("      FAIL quad closed form / A>0 at n=%d lam=%s"
                            " t=%s (A=%s)" % (n, lam, t, kcoef))
                        fails += 1
                        continue
                    if ueval(quad, lam * t) != want or want <= 0 or mixed:
                        log("      FAIL exterior value / eq:MPmixed at n=%d"
                            " lam=%s t=%s" % (n, lam, t))
                        fails += 1
                        continue
                    ys = roots_of(G, "y")
                    zs = roots_of(Pl, "z")
                    qz = roots_of(Q, "q")
                    es = roots_of(quad, "E")
                    if len(ys) != n + 1 or len(zs) != n or len(es) != 2 \
                            or len(qz) != n:
                        log("      FAIL zero counts n=%d lam=%s t=%s"
                            % (n, lam, t))
                        fails += 1
                        continue
                    gq, why = strictly_interlaces(ys, qz)
                    if not gq:
                        log("      FAIL G < Q at n=%d lam=%s t=%s: %s"
                            % (n, lam, t, why))
                        fails += 1
                        continue
                    k1, k2 = gap_index(es[0], ys), gap_index(es[1], ys)
                    if not (1 <= k1 <= n and 1 <= k2 <= n and k1 != k2):
                        continue                    # hypothesis not met: skip
                    hyp += 1
                    hyp_by_branch[sgn] += 1
                    # T'(c) = (2lam+n+1) sin 2phi has the sign of t; the paper
                    # says T'(c)>0 <=> c>E_2 and T'(c)<0 <=> c<E_1.  Decide it.
                    c = rational_as_root(lam * t, "c")
                    if sgn > 0:
                        if rcmp(c, es[1]) <= 0:
                            branchbad += 1
                    else:
                        if rcmp(c, es[0]) >= 0:
                            branchbad += 1
                    if sgn > 0:
                        # phi in (0,phi_c): (x-E_2)G < (x-E_1)P
                        big = [es[1]] + list(ys)
                        small = [es[0]] + list(zs)
                    else:
                        big = [es[0]] + list(ys)
                        small = [es[1]] + list(zs)
                    okc, why = strictly_interlaces(big, small)
                    if not okc:
                        log("      FAIL n=%d lam=%s t=%s gaps=(%d,%d): %s"
                            % (n, lam, t, k1, k2, why))
                        fails += 1
                        continue
                    # THE PAPER'S MECHANISM, on the real family.  The proof of
                    # Lemma 2.1 forces P to have exactly TWO zeros in the gap of
                    # G containing the completion point rho that faces c, those
                    # two straddling rho; ZERO zeros in the gap containing the
                    # other completion point sigma; and exactly one in each
                    # remaining interior gap.  That straddling configuration is
                    # precisely the placement hypothesis [JK, Thm 2.2(e)] has to
                    # ASSUME and this paper claims to DERIVE, so it is checked
                    # here on the actual Meixner-Pollaczek zeros rather than
                    # only on the rigid rational instances of check 10.
                    krho, ksig = (k2, k1) if sgn > 0 else (k1, k2)
                    rho = es[1] if sgn > 0 else es[0]
                    prof = zero_profile(zs, ys)
                    wantprof = dict((j, 2 if j == krho else 1)
                                    for j in range(1, n + 1) if j != ksig)
                    below = sum(1 for z in zs
                                if gap_index(z, ys) == krho and rcmp(z, rho) < 0)
                    if prof != wantprof or below != 1:
                        log("      FAIL zero profile n=%d lam=%s t=%s:"
                            " got %s want %s, below rho %d (want 1)"
                            % (n, lam, t, prof, wantprof, below))
                        placement_bad += 1
                        fails += 1
                    else:
                        placement_ok += 1
    log("    %d grid points; hypothesis (E1,E2 in distinct interior gaps of G)"
        " met on %d; interlacing failures %d.  The theorem's CONCLUSION is"
        " tested on those %d only; the other %d exercise the mixed identity,"
        " G < Q, disc T > 0 and A > 0 but no conclusion"
        % (npts, hyp, fails, hyp, npts - hyp))
    log("    hypothesis-met points BY BRANCH: t>0 (phi in (0,phi_c)) %d,"
        " t<0 (phi in (pi-phi_c,pi)) %d -- both must be nonzero or one of the"
        " theorem's two conclusions is never tested"
        % (hyp_by_branch[1], hyp_by_branch[-1]))
    log("    branch test: sgn T'(c) = sgn(t) agreed with c>E_2 (t>0) resp."
        " c<E_1 (t<0) on %d/%d hypothesis-met points" % (hyp - branchbad, hyp))
    log("    [JK, Thm 2.2(e)] placement DERIVED (P has 2 zeros straddling rho,"
        " 0 in sigma's gap, 1 elsewhere) on %d/%d points, %d violations"
        % (placement_ok, hyp, placement_bad))
    log("    SAMPLING: grid shape %d = %d n-values x %d lambda-values x %d"
        " ratios x %d branches; n in %s, lambda in %s.  Theorem 3.1 is stated"
        " for EVERY n>=2 and EVERY lambda>0, so this is a sample, not a proof."
        % (npts, len(ns_seen), len(lams_seen), len(MUS), len(SGNS),
           sorted(ns_seen), sorted(str(x) for x in lams_seen)))
    log("    BOOKKEEPING: this program's own check title declares a %d-point"
        " grid; the loops enumerated %d and recorded %d failure(s).  This"
        " compares the program against its OWN declared shape -- the paper"
        " contains no numerical experiment, no grid and no failure count, so"
        " there is no published figure here to agree with."
        % (DECLARED_FIGURES["mp_grid_points"], npts, fails))
    SWEEP_STATS["mp"] = {"npts": npts, "hyp": hyp}
    return (fails == 0 and branchbad == 0 and hyp > 0
            and hyp_by_branch[1] > 0 and hyp_by_branch[-1] > 0
            and placement_ok == hyp
            and npts == len(NS) * len(LAMS) * len(MUS) * len(SGNS)
            and npts == DECLARED_FIGURES["mp_grid_points"])


def snap_down(x, den=10 ** 6):
    """Largest rational <= x with denominator dividing den."""
    x = F(x)
    return F((x.numerator * den) // x.denominator, den)


def snap_up(x, den=10 ** 6):
    """Smallest rational >= x with denominator dividing den."""
    x = F(x)
    return F(-((-x.numerator * den) // x.denominator), den)


def pj_range_bound_lower(n, digits=12):
    """
    A rational number <= the eq:PJrange bound (-3n-5-sqrt(n^2+2n+5))/4,
    obtained from an integer square root.  Used only to place sample points;
    membership is always re-tested exactly by in_pj_range.
    """
    s_up = rat_sqrt_lower(F(pj_radicand(n)), digits) + F(1, 10 ** digits)
    return (F(-3 * n - 5) - s_up) / 4


def sweep_theorem42(log):
    """
    CHECK Theorem 4.2 end to end on a 54-point grid: n in {2,3,4},
    a = (eq:PJrange bound) - delta for delta in {1/20, 3/5, 2}, |b| at
    {51/50, 3/2, 4} times the eq:brange critical value, both signs of b.
    Every quantity is derived: membership in eq:PJrange and eq:brange (exact
    rational tests), A > 0, u < -1/2, D > 0, disc R > 0, R(c) > 0, eq:PJmixed
    as an exact polynomial identity at this point, the interlacing G < Q of
    Lemma 4.1, the zeros of G = P_{n+1}(x;a+1,b), Q = P_n(x;a+1,b) and
    P = P_n(x;a,b), the roots e_- < e_+ of R, whether e_- and e_+ lie in
    DISTINCT INTERIOR gaps of G (Theorem 4.2's only hypothesis), Delta_c and
    its agreement with R'(c) and with the side of [e_-,e_+] on which c = b/u
    falls, and finally the asserted interlacing.
    The delta = 1/20 points lie in the RESIDUAL BAND (-n-2, bound) on which
    [JK, Corollary 6.1] does not apply: those are the new ground.
    """
    npts = hyp = fails = band = 0
    branchbad = 0
    # band_hyp: residual-band points that ACTUALLY reached and passed the
    # conclusion.  `band` alone is counted before the eq:brange test, before the
    # algebraic preconditions and before the hypothesis test, so "band > 0" says
    # nothing about whether the ground [JK, Cor 6.1] does not cover was tested.
    band_hyp = 0
    hyp_by_sign = {1: 0, -1: 0}
    placement_ok = placement_bad = 0
    iff_ok = na_pairs = 0
    ns_seen, as_seen = set(), set()
    NS = (2, 3, 4)
    DELTAS = (F(1, 20), F(3, 5), F(2))
    RATIOS = (F(51, 50), F(3, 2), F(4))
    SIGNS = (1, -1)
    for n in NS:
        blo = pj_range_bound_lower(n)
        for delta in DELTAS:
            a = snap_down(blo - delta, 10 ** 3)
            if not in_pj_range(n, a):
                log("      FAIL: sample a=%s not in eq:PJrange (n=%d)"
                    % (a, n))
                fails += 1
                continue
            inband = a > -n - 2
            na_pairs += 1
            ns_seen.add(n)
            as_seen.add(a)
            bcrit2, Dq = brange_bound_sq(n, a)
            bcl = rat_sqrt_lower(bcrit2, 12)
            # [JK, Lemma 6.1] is an IF AND ONLY IF, and the paper relies on that
            # to say the parameter range is not narrowed.  The grid below only
            # ever samples b INSIDE eq:brange, so it can only ever confirm
            # eq:brange => disc R > 0.  Test the other direction here: below the
            # critical |b| the quadratic must NOT have two real roots.  Without
            # this the program could not distinguish eq:brange from any weaker
            # sufficient condition.
            for sub in (F(49, 50), F(1, 2)):
                for s0 in (1, -1):
                    b0 = s0 * snap_down(sub * bcl, 10 ** 3)
                    d0 = pj_R_coeffs(n, a, b0)
                    disc0 = d0[1] ** 2 - 4 * d0[0]
                    if in_brange(n, a, b0) or disc0 >= 0:
                        log("      FAIL eq:brange 'only if' at n=%d a=%s b=%s:"
                            " in_brange=%s disc=%s"
                            % (n, a, b0, in_brange(n, a, b0), disc0))
                        fails += 1
                    else:
                        iff_ok += 1
            for ratio in RATIOS:
                for s in SIGNS:
                    b = s * snap_up(ratio * bcl, 10 ** 3)
                    npts += 1
                    if inband:
                        band += 1
                    if not in_brange(n, a, b):
                        log("      FAIL: sample b=%s not in eq:brange" % b)
                        fails += 1
                        continue
                    u = a + n + 1
                    A = pj_A(n, a, b)
                    Rc = pj_R_coeffs(n, a, b)
                    disc = Rc[1] ** 2 - 4 * Rc[0]
                    c = b / u
                    Rat_c = ueval(Rc, c)
                    dlt = b * (2 * a + n + 3) / ((a + n + 1) * (a + n + 2))
                    Rp_c = ueval(uderiv(Rc), c)
                    if not (A > 0 and disc > 0 and Rat_c > 0 and Dq > 0
                            and dlt != 0 and Rp_c == dlt and u < F(-1, 2)):
                        log("      FAIL algebraic preconditions n=%d a=%s b=%s"
                            " A=%s disc=%s R(c)=%s D=%s R'(c)-Delta=%s"
                            % (n, a, b, A > 0, disc > 0, Rat_c > 0, Dq > 0,
                               Rp_c - dlt))
                        fails += 1
                        continue
                    Pu = pj_monic_num(n + 1, a + 1, b)[0]
                    G, Q = Pu[n + 1], Pu[n]
                    P = pj_monic_num(n, a, b)[0][n]
                    mixed = usub(uscal(P, A),
                                 uadd(uscal(umul(upow_lin(c), G), F(-1)),
                                      umul(Rc, Q)))
                    if mixed:
                        log("      FAIL eq:PJmixed at n=%d a=%s b=%s"
                            % (n, a, b))
                        fails += 1
                        continue
                    ys, zs, es = roots_of(G, "y"), roots_of(P, "z"), \
                        roots_of(Rc, "e")
                    qz = roots_of(Q, "q")
                    if len(ys) != n + 1 or len(zs) != n or len(es) != 2 \
                            or len(qz) != n:
                        log("      FAIL zero counts n=%d a=%s b=%s: %d %d %d %d"
                            % (n, a, b, len(ys), len(zs), len(es), len(qz)))
                        fails += 1
                        continue
                    gq, why = strictly_interlaces(ys, qz)
                    if not gq:
                        log("      FAIL G < Q (Lemma 4.1) at n=%d a=%s b=%s:"
                            " %s" % (n, a, b, why))
                        fails += 1
                        continue
                    k1, k2 = gap_index(es[0], ys), gap_index(es[1], ys)
                    if not (1 <= k1 <= n and 1 <= k2 <= n and k1 != k2):
                        continue
                    hyp += 1
                    hyp_by_sign[1 if dlt > 0 else -1] += 1
                    cr = rational_as_root(c, "c")
                    if dlt < 0:
                        if rcmp(cr, es[0]) >= 0:
                            branchbad += 1
                        big, small = [es[0]] + list(ys), [es[1]] + list(zs)
                        krho, ksig, rho = k1, k2, es[0]
                    else:
                        if rcmp(cr, es[1]) <= 0:
                            branchbad += 1
                        big, small = [es[1]] + list(ys), [es[0]] + list(zs)
                        krho, ksig, rho = k2, k1, es[1]
                    okc, why = strictly_interlaces(big, small)
                    if not okc:
                        log("      FAIL n=%d a=%s b=%s gaps=(%d,%d): %s"
                            % (n, a, b, k1, k2, why))
                        fails += 1
                        continue
                    # the paper's mechanism on the real family: P has exactly
                    # two zeros straddling rho inside rho's gap of G (this is
                    # the placement [JK, Thm 2.2(e)] assumes and the paper
                    # derives), none in sigma's gap, one in each other gap
                    prof = zero_profile(zs, ys)
                    wantprof = dict((j, 2 if j == krho else 1)
                                    for j in range(1, n + 1) if j != ksig)
                    below = sum(1 for z in zs
                                if gap_index(z, ys) == krho and rcmp(z, rho) < 0)
                    if prof != wantprof or below != 1:
                        log("      FAIL zero profile n=%d a=%s b=%s: got %s"
                            " want %s, below rho %d (want 1)"
                            % (n, a, b, prof, wantprof, below))
                        placement_bad += 1
                        fails += 1
                        continue
                    placement_ok += 1
                    if inband:
                        band_hyp += 1
    # If the fixed 54-point grid happened to leave the residual band untested,
    # search the band directly rather than passing on a claim we did not test.
    # This does NOT touch the 54-point grid or its counters.
    extra_looked = 0
    if band_hyp == 0:
        log("    residual band untested by the fixed grid; searching the band")
        for n in (2, 3, 4):
            bnd = pj_range_bound_lower(n)
            width = bnd + n + 2                       # band is (-n-2, bound)
            if width <= 0 or band_hyp:
                continue
            for frac in (F(1, 50), F(1, 10), F(3, 10), F(1, 2), F(7, 10),
                         F(9, 10)):
                if band_hyp:
                    break
                a2 = snap_down(bnd - width * frac, 10 ** 4)
                if not (in_pj_range(n, a2) and a2 > -n - 2 and a2 != -n - 2):
                    continue
                bb2, Dq2 = brange_bound_sq(n, a2)
                bcl2 = rat_sqrt_lower(bb2, 12)
                for ratio in (F(21, 20), F(6, 5), F(2), F(4), F(10), F(40)):
                    b2 = snap_up(ratio * bcl2, 10 ** 3)
                    extra_looked += 1
                    if not in_brange(n, a2, b2):
                        continue
                    u2 = a2 + n + 1
                    A2, Rc2 = pj_A(n, a2, b2), pj_R_coeffs(n, a2, b2)
                    if not (A2 > 0 and Rc2[1] ** 2 - 4 * Rc2[0] > 0
                            and u2 < F(-1, 2) and Dq2 > 0):
                        continue
                    Pu2 = pj_monic_num(n + 1, a2 + 1, b2)[0]
                    G2, Q2 = Pu2[n + 1], Pu2[n]
                    P2 = pj_monic_num(n, a2, b2)[0][n]
                    if usub(uscal(P2, A2),
                            uadd(uscal(umul(upow_lin(b2 / u2), G2), F(-1)),
                                 umul(Rc2, Q2))):
                        log("      FAIL eq:PJmixed in band search n=%d a=%s"
                            " b=%s" % (n, a2, b2))
                        fails += 1
                        continue
                    ys2, zs2 = roots_of(G2, "y"), roots_of(P2, "z")
                    es2, qz2 = roots_of(Rc2, "e"), roots_of(Q2, "q")
                    if len(ys2) != n + 1 or len(zs2) != n or len(es2) != 2 \
                            or len(qz2) != n:
                        continue
                    if not strictly_interlaces(ys2, qz2)[0]:
                        log("      FAIL G<Q in band search n=%d a=%s b=%s"
                            % (n, a2, b2))
                        fails += 1
                        continue
                    j1, j2 = gap_index(es2[0], ys2), gap_index(es2[1], ys2)
                    if not (1 <= j1 <= n and 1 <= j2 <= n and j1 != j2):
                        continue
                    d2 = b2 * (2 * a2 + n + 3) / ((a2 + n + 1) * (a2 + n + 2))
                    if d2 < 0:
                        bg, sm = [es2[0]] + list(ys2), [es2[1]] + list(zs2)
                    else:
                        bg, sm = [es2[1]] + list(ys2), [es2[0]] + list(zs2)
                    okb, whyb = strictly_interlaces(bg, sm)
                    if not okb:
                        log("      FAIL band-search conclusion n=%d a=%s b=%s:"
                            " %s" % (n, a2, b2, whyb))
                        fails += 1
                        continue
                    band_hyp += 1
                    log("      band point VERIFIED: n=%d a=%s b=%s"
                        " Delta_c sign %d" % (n, a2, b2, 1 if d2 > 0 else -1))
                    break
        log("    band search examined %d candidate points" % extra_looked)
    log("    %d grid points (%d of them in the residual band (-n-2, bound)"
        " where [JK, Cor 6.1] does not apply); hypothesis met on %d;"
        " failures %d.  The theorem's CONCLUSION is tested on those %d only;"
        " the other %d exercise eq:PJmixed, G < Q, disc R > 0, R(c) > 0 and"
        " A > 0 but no conclusion"
        % (npts, band, hyp, fails, hyp, npts - hyp))
    log("    RESIDUAL BAND coverage: %d point(s) in (-n-2, bound) met the"
        " hypothesis AND had their conclusion verified, out of %d band points"
        " on the fixed grid plus %d examined by the band search.  This, not"
        " `band`, is what makes the claim to extend [JK, Cor 6.1] tested:"
        " `band` alone is counted before every gate."
        % (band_hyp, band, extra_looked))
    log("    hypothesis-met points BY SIGN of Delta_c: Delta_c<0 %d,"
        " Delta_c>0 %d -- both must be nonzero or one of the theorem's two"
        " conclusions is never tested"
        % (hyp_by_sign[-1], hyp_by_sign[1]))
    log("    branch test: sgn Delta_c agreed with the side of [e-,e+] that c"
        " lies on, on %d/%d hypothesis-met points" % (hyp - branchbad, hyp))
    log("    [JK, Thm 2.2(e)] placement DERIVED on %d/%d points, %d violations"
        % (placement_ok, hyp, placement_bad))
    log("    eq:brange 'only if' direction (|b| BELOW critical => disc R < 0,"
        " so eq:brange is not merely sufficient): %d points" % iff_ok)
    log("    SAMPLING: grid shape %d = %d n-values x %d a-offsets x %d"
        " |b|-ratios x %d signs; n in %s, a in %s.  Theorem 4.2 is stated for"
        " EVERY n>=2 and EVERY (a,b) satisfying eq:PJrange/eq:brange, so this"
        " is a sample, not a proof."
        % (npts, len(ns_seen), len(DELTAS), len(RATIOS), len(SIGNS),
           sorted(ns_seen), sorted(str(x) for x in as_seen)))
    log("    BOOKKEEPING: this program's own check title declares a %d-point"
        " grid; the loops enumerated %d and recorded %d failure(s).  This"
        " compares the program against its OWN declared shape -- the paper"
        " contains no numerical experiment, no grid and no failure count, so"
        " there is no published figure here to agree with."
        % (DECLARED_FIGURES["pj_grid_points"], npts, fails))
    SWEEP_STATS["pj"] = {"npts": npts, "hyp": hyp, "band_hyp": band_hyp}
    return (fails == 0 and branchbad == 0 and hyp > 0 and band > 0
            and band_hyp > 0
            and hyp_by_sign[1] > 0 and hyp_by_sign[-1] > 0
            and placement_ok == hyp and na_pairs > 0
            and iff_ok == 4 * na_pairs
            and npts == len(NS) * len(DELTAS) * len(RATIOS) * len(SIGNS)
            and npts == DECLARED_FIGURES["pj_grid_points"])


def check_range_equivalence(log):
    """
    CHECK, on an exact grid that straddles the boundary, that on the half line
    u = a+n+1 < -1/2 the two independently computed predicates agree:
      in_pj_range(n,a)  [ a < (-3n-5-sqrt(n^2+2n+5))/4 , tested by squaring ]
      Dq(n,a) > 0       [ D = 4u^2-(2n-2)u-(n+1), the quantity inside the
                          absolute value of eq:brange ]
    The restriction to u < -1/2 is NECESSARY and is itself derived: D is an
    upward parabola with roots u_-- = (n-1-s)/4 and u_++ = (n-1+s)/4, so D>0
    also holds for u > u_++; we check u_++ > -1/2 for every n in the range (in
    fact u_++ > n/2), which is what makes the restriction harmless -- and check
    that eq:PJrange never holds when u >= -1/2, so nothing is lost.
    Equality of the two predicates on u < -1/2 is what makes eq:brange's |D|
    equal to D and makes disc R > 0 equivalent to eq:brange (check 5, item iii).
    Also: u < -1/2, a+1 < -n-1/2 and A > 0 wherever eq:PJrange holds; and the
    number of interior gaps of G is n, so "distinct k,k' in {1,..,n}" is
    unsatisfiable for n = 1 -- the paper's vacuity claim.
    """
    ok, ncmp, ndis, nabove = True, 0, 0, 0
    for n in range(1, 8):
        blo = pj_range_bound_lower(n)
        # u_++ = (n-1+s)/4 > n/2 > -1/2 : check by squaring, exactly
        s2 = pj_radicand(n)
        if not (s2 > (n + 1) ** 2 and F(n - 1, 4) + F(n + 1, 4) > F(-1, 2)):
            log("      FAIL u_++ > -1/2 at n=%d" % n)
            ok = False
        for off in (F(-3), F(-1), F(-1, 10), F(-1, 1000), F(0), F(1, 1000),
                    F(1, 100)):
            a = snap_down(blo + off, 10 ** 9)
            if a == -n - 2:
                continue
            u = a + n + 1
            assert u < F(-1, 2), "grid point left the u < -1/2 half line"
            r1 = in_pj_range(n, a)
            _, Dq = brange_bound_sq(n, a)
            r2 = Dq > 0
            ncmp += 1
            if r1 != r2:
                ndis += 1
                log("      DISAGREE n=%d a=%s: in_pj_range=%s D>0=%s"
                    % (n, a, r1, r2))
                ok = False
            if r1 and not (u < F(-1, 2) and a + 1 < F(-2 * n - 1, 2)
                           and pj_A(n, a, F(1)) > 0):
                log("      FAIL u<-1/2 / A>0 at n=%d a=%s" % (n, a))
                ok = False
        for off in (F(0), F(1, 10), F(1), F(3), F(10)):    # u >= -1/2 region
            a = F(-2 * n - 3, 2) + off
            nabove += 1
            if in_pj_range(n, a):
                log("      FAIL eq:PJrange held at u>=-1/2: n=%d a=%s"
                    % (n, a))
                ok = False
    log("    %d (n,a) points straddling the eq:PJrange boundary on u<-1/2;"
        " eq:PJrange <=> D>0 disagreements: %d" % (ncmp, ndis))
    log("    %d points with u >= -1/2: eq:PJrange held on 0 of them, so the"
        " restriction discards nothing" % nabove)
    # The n=1 vacuity claim, DERIVED rather than asserted.  The old form of this
    # check compared len(range(1,m+1)) with m and len of an enumerated pair list
    # with 0 and 2 -- constant-true integer identities that mention no
    # polynomial and cannot fail.  What the claim actually needs is that
    # G = P_{n+1} really has n+1 distinct real zeros, hence n interior gaps,
    # hence no DISTINCT pair (k,k') in {1..n} when n = 1.  So build G at n = 1
    # in BOTH families from the recurrences and count its zeros.
    def interior_gaps_of_G(Gpoly, label):
        zs = roots_of(Gpoly, "y")
        return len(zs) - 1, len(zs), label

    vac = []
    # Meixner-Pollaczek, n=1: lambda=1, t=cot phi=2, admissible since
    # (n+1)t^2 = 8 > 2 lambda = 2; G = P_{n+1}^{(lambda+1)} = P_2^{(2)}.
    vac.append(interior_gaps_of_G(mp_monic_num(2, F(2), F(2))[2], "MP n=1"))
    # pseudo-Jacobi, n=1: a=-7/2 (in eq:PJrange for n=1, a != -n-2=-3),
    # so G = P_2(.; a+1=-5/2, b=3) and Lemma 4.1 applies.
    if not in_pj_range(1, F(-7, 2)):
        log("      FAIL: a=-7/2 is not in eq:PJrange for n=1")
        ok = False
    vac.append(interior_gaps_of_G(pj_monic_num(2, F(-5, 2), F(3))[0][2],
                                  "PJ n=1"))
    for (ngaps, nzeros, label) in vac:
        distinct_pairs = [(k, kp) for k in range(1, ngaps + 1)
                          for kp in range(1, ngaps + 1) if k != kp]
        log("    %s: G has %d distinct real zeros -> %d interior gap(s) ->"
            " %d distinct (k,k') pairs" % (label, nzeros, ngaps,
                                           len(distinct_pairs)))
        if not (nzeros == 2 and ngaps == 1 and distinct_pairs == []):
            log("      FAIL: n=1 vacuity not derived for %s" % label)
            ok = False
    # and, for contrast, that n=2 does admit a distinct pair (so the theorems
    # are not vacuous for the n they are stated for)
    g2 = len(roots_of(mp_monic_num(3, F(2), F(2))[3], "y")) - 1
    p2 = [(k, kp) for k in range(1, g2 + 1) for kp in range(1, g2 + 1)
          if k != kp]
    log("    MP n=2: G has %d interior gaps -> %d distinct (k,k') pairs"
        " (nonvacuous)" % (g2, len(p2)))
    if not (g2 == 2 and len(p2) == 2):
        log("      FAIL: n=2 non-vacuity not derived")
        ok = False
    return ok


# =========================================================================
# 8. Driver.
# =========================================================================

def main():

    def log(s):
        print(s)

    checks = [
        ("eq:MPmixed is an identity in (x, lambda, cot phi) for n=2,3,4",
         lambda: check_mp_mixed_identity([2, 3, 4], log)),
        ("MP exterior point: (c-E1)(c-E2), T'(c), disc T, cot^2 phi_c",
         lambda: check_mp_exterior_point([2, 3, 4], log)),
        ("eq:PJmixed is an identity in (x, b) at in-range rational a, n=2..5",
         lambda: check_pj_mixed_identity(
             [(2, F(-9, 2)), (2, F(-31, 5)), (3, F(-31, 5)), (3, F(-8)),
              (4, F(-13, 2)), (4, F(-17, 2)), (5, F(-15, 2)),
              (5, F(-11))], log)),
        ("eq:PJrange is exactly what eq:PJrec needs (the a=-9/2, n=3 pole)",
         lambda: check_pj_range_is_needed(log)),
        ("R(b/u), R'(b/u), disc R and the Pochhammer A, symbolic in (n,a,b)",
         lambda: check_pj_symbolic_nab(log)),
        ("eq:PJrange arithmetic: band nonempty, bound < -n-3/2, D(u_-)=0",
         lambda: check_range_arithmetic(40, log)),
        ("eq:PJrange <=> D>0, u<-1/2, A>0, and the n=1 vacuity",
         lambda: check_range_equivalence(log)),
        ("gamma_m > 0 for 1<=m<=n when alpha < -n-1/2",
         lambda: check_gamma_positive(log)),
        ("Lemma 4.1: P_{n+1}(.;alpha,b) < P_n(.;alpha,b), exact zeros",
         lambda: check_lemma41_interlacing(log)),
        ("Lemma 2.1 (exterior-point criterion) and every step of its proof",
         lambda: check_lemma_criterion(240, log)),
        ("Theorem 3.1 end to end, exact, 72-point grid",
         lambda: sweep_theorem31(log)),
        ("Theorem 4.2 end to end, exact, 54-point grid",
         lambda: sweep_theorem42(log)),
    ]

    results = []
    for i, (name, fn) in enumerate(checks, 1):
        log("--- check %d: %s" % (i, name))
        try:
            ok = bool(fn())
        except Exception as exc:                       # a crash is a failure
            log("    EXCEPTION: %s: %s" % (type(exc).__name__, exc))
            ok = False
        results.append(ok)
        log("%s check %d: %s" % ("PASS" if ok else "FAIL", i, name))
        log("")

    nfail = results.count(False)
    log(NOTE_ON_NUMERICS)
    log("")
    log(GAPS_NOT_COVERED)
    log("")
    if nfail == 0:
        mp = SWEEP_STATS.get("mp", {})
        pj = SWEEP_STATS.get("pj", {})
        log("VERDICT: ALL %d CHECKS PASS." % len(results))
        log("VERDICT SCOPE: this confirms the paper's algebra symbolically"
            " (in (x,lambda,cot phi) for n=2,3,4 in the Meixner-Pollaczek case;"
            " in (n,a,b) for the pseudo-Jacobi rational-function identities;"
            " at 8 (n,a) pairs for eq:PJmixed).  It confirms both theorems'"
            " CONCLUSIONS at the %d of %d Meixner-Pollaczek and %d of %d"
            " pseudo-Jacobi exactly-computed parameter points at which each"
            " theorem's hypothesis actually holds: a conclusion can only be,"
            " and only is, tested where the hypothesis holds.  %d verified"
            " point(s) lie in the residual band (-n-2, bound) on which"
            " [JK, Cor 6.1] does not apply -- from the fixed grid, or from the"
            " supplementary band search if the fixed grid left the band"
            " untested; see check 12.  The remaining %d and"
            " %d sampled points do exercise the mixed identity, G < Q, disc > 0"
            " and A > 0, but no conclusion.  It is NOT a proof of either"
            " theorem for all n and all parameters; see GAPS NOT COVERED above"
            " for the steps that remain on the paper's own authority."
            % (mp.get("hyp", -1), mp.get("npts", -1),
               pj.get("hyp", -1), pj.get("npts", -1), pj.get("band_hyp", -1),
               mp.get("npts", 0) - mp.get("hyp", 0),
               pj.get("npts", 0) - pj.get("hyp", 0)))
    else:
        log("VERDICT: %d OF %d CHECKS FAILED" % (nfail, len(results)))
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Verification program for the alpha = 0 Jacobi counterexample paper.

Two independent arithmetic layers:

  (A) EXACT.  The paper's decisive content -- the hypergeometric
      coefficient formulas, the n = 1 identity, the endpoint identity
      Delta_2(pi) = -q(beta)/384, the location of beta_*, and the
      Corollary's rational witness -- is checked in integer/Fraction
      arithmetic, symbolically in beta where the claim is symbolic.

  (B) NUMERICAL.  The contestable paragraph resting on the 1000-node
      sweep theta_nu = nu*pi/1001 at beta = 2 is reproduced in
      decimal.Decimal arithmetic with pi, cos, sin and acos built here
      from series with stated error bounds.  No bare float is used in
      any load-bearing comparison.

Nothing is transcribed: every polynomial used numerically is derived
from the hypergeometric definition and then compared, as an identity in
beta, against the coefficients the paper displays.

Standard library only.  Python 3.9 or later.
"""

import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

PREC = 60
getcontext().prec = PREC

# Truncation threshold for every series below.  Each series used here is
# alternating with eventually decreasing terms, so the truncation error
# is bounded by the first omitted term, i.e. by EPS_SERIES = 1e-55.
# Rounding contributes at most ~50 operations of relative size 1e-60 on
# partial sums of modulus <= 5, i.e. under 1e-57.  Every series value
# below is therefore accurate to better than 2e-55 absolutely.
EPS_SERIES = Decimal(1).scaleb(-(PREC - 5))

CHECKS = []


def ck(name, ok, detail=""):
    """Record and print one check.  Returns the boolean."""
    ok = bool(ok)
    CHECKS.append((name, ok))
    print("%s  %s%s" % ("PASS" if ok else "FAIL", name,
                        ("  |  " + detail) if detail else ""))
    return ok


# ----------------------------------------------------------------------
# (A) Exact layer: polynomials in two variables over Q.
#
# A polynomial is a dict {(i, j): Fraction} standing for the sum of
# coefficient * u**i * b**j, where u is the paper's u = 1 - cos(theta/2)
# and b is the paper's beta.  Empty dict = zero.
# ----------------------------------------------------------------------

def const(c):
    c = Fraction(c)
    return {} if c == 0 else {(0, 0): c}


B = {(0, 1): Fraction(1)}        # beta
U = {(1, 0): Fraction(1)}        # u
ONE = const(1)


def padd(*ps):
    r = {}
    for p in ps:
        for k, v in p.items():
            nv = r.get(k, Fraction(0)) + v
            if nv == 0:
                r.pop(k, None)
            else:
                r[k] = nv
    return r


def pscale(p, c):
    c = Fraction(c)
    return {} if c == 0 else {k: v * c for k, v in p.items()}


def psub(p, q):
    return padd(p, pscale(q, -1))


def pmul(p, q):
    r = {}
    for (i1, j1), v1 in p.items():
        for (i2, j2), v2 in q.items():
            k = (i1 + i2, j1 + j2)
            nv = r.get(k, Fraction(0)) + v1 * v2
            if nv == 0:
                r.pop(k, None)
            else:
                r[k] = nv
    return r


def ppow(p, n):
    r = ONE
    for _ in range(n):
        r = pmul(r, p)
    return r


def peval(p, uval, bval):
    """Value at u = uval, b = bval, as an exact Fraction."""
    uval, bval = Fraction(uval), Fraction(bval)
    return sum((v * uval ** i * bval ** j for (i, j), v in p.items()),
               Fraction(0))


def subst_u(p, uval):
    """Substitute a rational for u; result is a polynomial in b alone."""
    uval = Fraction(uval)
    r = {}
    for (i, j), v in p.items():
        k = (0, j)
        nv = r.get(k, Fraction(0)) + v * uval ** i
        if nv == 0:
            r.pop(k, None)
        else:
            r[k] = nv
    return r


def subst_b(p, q):
    """Substitute the polynomial q for b."""
    if not p:
        return {}
    maxj = max(j for (_, j) in p)
    groups = {}
    for (i, j), v in p.items():
        groups.setdefault(j, {})[(i, 0)] = v
    acc, qp = {}, ONE
    for j in range(maxj + 1):
        if j in groups:
            acc = padd(acc, pmul(groups[j], qp))
        qp = pmul(qp, q)
    return acc


def positive_for_b_at_least_one(p):
    """Rigorous certificate that p(b) > 0 for every real b >= 1.

    Substitute b = 1 + s.  If every coefficient of the result is >= 0
    and the constant term is > 0, then p > 0 for all s >= 0.
    """
    shifted = subst_b(p, padd(ONE, U))
    if any(j != 0 for (_, j) in shifted):
        return False, "not a polynomial in b alone"
    coeffs = {}
    for (i, _), v in shifted.items():
        coeffs[i] = v
    if any(v < 0 for v in coeffs.values()):
        return False, "shifted polynomial has a negative coefficient"
    if coeffs.get(0, Fraction(0)) <= 0:
        return False, "shifted constant term is not positive"
    terms = " + ".join("%s*s^%d" % (coeffs[i], i) for i in sorted(coeffs))
    return True, "b = 1+s gives " + terms


# ---- the hypergeometric definition, symbolic in beta at alpha = 0 ----

def jacobi_alpha0_zcoeffs(n):
    """Coefficients of z**k in P_n^{(0,beta)}(1-2z) = 2F1(-n, n+b+1; 1; z).

    2F1 term k is (-n)_k (n+b+1)_k / ((1)_k k!) = (-1)^k C(n,k) (b')_k / k!
    with b' = n+b+1, since (1)_k = k!.  Each returned entry is a
    polynomial in beta.
    """
    bp = padd(B, const(n + 1))                     # n + beta + 1
    out = []
    for k in range(n + 1):
        term = const(Fraction((-1) ** k * math.comb(n, k),
                              math.factorial(k)))
        for j in range(k):
            term = pmul(term, padd(bp, const(j)))
        out.append(term)
    return out


def poly_in_z(zcoeffs, z):
    """Evaluate sum_k zcoeffs[k] * z**k with z a polynomial."""
    acc, zp = {}, ONE
    for c in zcoeffs:
        acc = padd(acc, pmul(c, zp))
        zp = pmul(zp, z)
    return acc


def tilde_zcoeffs(n, alpha, beta):
    """Coefficients of z**k in the normalized polynomial for general alpha:

        Ptilde_n^{(alpha,beta)}(1-2z) = 2F1(-n, n+alpha+beta+1; alpha+1; z)

    with term k equal to (-1)^k C(n,k) (n+alpha+beta+1)_k / (alpha+1)_k.
    Returns exact Fractions.
    """
    alpha, beta = Fraction(alpha), Fraction(beta)
    bp = n + alpha + beta + 1
    out = []
    for k in range(n + 1):
        num = Fraction((-1) ** k * math.comb(n, k))
        for j in range(k):
            num *= (bp + j)
        den = Fraction(1)
        for j in range(k):
            den *= (alpha + 1 + j)
        out.append(num / den)
    return out


def eval_fr(zcoeffs, z):
    """Horner evaluation of exact Fraction coefficients at rational z."""
    z = Fraction(z)
    acc = Fraction(0)
    for c in reversed(zcoeffs):
        acc = acc * z + c
    return acc


P1C = jacobi_alpha0_zcoeffs(1)
P2C = jacobi_alpha0_zcoeffs(2)
P3C = jacobi_alpha0_zcoeffs(3)

# With c = cos(theta/2) and u = 1-c: the argument of P_2 gives
# z = (1-c)/2 = u/2, and the argument of P_1 gives
# z = (1-cos theta)/2 = 1-c^2 = u(2-u).
Z_FOR_P2 = pscale(U, Fraction(1, 2))
Z_FOR_P1 = pmul(U, psub(const(2), U))

# The affine factor of the paper's n = 1 identity.
AFFINE = padd(padd(B, ONE),
              pscale(pmul(psub(psub(ppow(B, 2), B), const(4)), U),
                     Fraction(1, 8)))


def u_coeffs_at_beta(p, bval):
    """Coefficients of u**i after substituting a rational for beta."""
    bval = Fraction(bval)
    maxi = max((i for (i, _) in p), default=0)
    out = [Fraction(0)] * (maxi + 1)
    for (i, j), v in p.items():
        out[i] += v * bval ** j
    return out

# q(beta) = beta^3 + 27 beta^2 - 10 beta - 24
Q = padd(ppow(B, 3), pscale(ppow(B, 2), 27), pscale(B, -10), const(-24))
QPRIME = padd(pscale(ppow(B, 2), 3), pscale(B, 54), const(-10))


def exact_checks():
    # --- 1. derived coefficients equal the three displayed polynomials
    paper_P1 = [ONE, pscale(padd(B, const(2)), -1)]
    paper_P2 = [ONE,
                pscale(padd(B, const(3)), -2),
                pscale(pmul(padd(B, const(3)), padd(B, const(4))),
                       Fraction(1, 2))]
    paper_P3 = [ONE,
                pscale(padd(B, const(4)), -3),
                pscale(pmul(padd(B, const(4)), padd(B, const(5))),
                       Fraction(3, 2)),
                pscale(pmul(pmul(padd(B, const(4)), padd(B, const(5))),
                            padd(B, const(6))), Fraction(-1, 6))]
    ok = (P1C == paper_P1 and P2C == paper_P2 and P3C == paper_P3)
    ck("2F1 definition reproduces the displayed P_1, P_2, P_3 "
       "(identity in beta)", ok,
       "all %d coefficient polynomials agree exactly" %
       (len(P1C) + len(P2C) + len(P3C)))

    # --- 2. the n = 1 identity (see Z_FOR_P1, Z_FOR_P2 above)
    lhs = psub(poly_in_z(P2C, Z_FOR_P2), poly_in_z(P1C, Z_FOR_P1))
    affine = AFFINE
    rhs = pmul(U, affine)
    ck("identity  P_2(cos(t/2)) - P_1(cos t) = u*(b+1 + (b^2-b-4)u/8)",
       lhs == rhs, "difference is the zero polynomial in (u, beta)")

    # --- 3. endpoint values of the affine factor
    at0 = subst_u(affine, 0)
    at1 = subst_u(affine, 1)
    want0 = padd(B, ONE)
    want1 = pscale(padd(padd(ppow(B, 2), pscale(B, 7)), const(4)),
                   Fraction(1, 8))
    ck("affine factor has endpoint values b+1 and (b^2+7b+4)/8",
       at0 == want0 and at1 == want1)

    # --- 4. both endpoints positive for every b >= 1, hence (being
    #        affine in u) the factor is positive on 0 <= u <= 1
    ok0, why0 = positive_for_b_at_least_one(at0)
    ok1, why1 = positive_for_b_at_least_one(at1)
    ck("both endpoint values are positive for every beta >= 1, so the "
       "n=1 comparison holds on 0<theta<pi", ok0 and ok1,
       "b+1: %s ; (b^2+7b+4)/8: %s" % (why0, why1))

    # --- 5. endpoint identity Delta_2(pi) = -q(beta)/384.
    #        theta = pi gives cos(pi/2) = 0 -> z = 1/2 for P_2,
    #        and cos(pi/3) = 1/2 -> z = 1/4 for P_3.
    d2pi = psub(poly_in_z(P3C, const(Fraction(1, 4))),
                poly_in_z(P2C, const(Fraction(1, 2))))
    ck("identity  Delta_2(pi) = -q(beta)/384",
       d2pi == pscale(Q, Fraction(-1, 384)),
       "checked as an identity in beta, not at sample points")

    # --- 6. q(1) < 0 < q(2) and q' > 0 for b >= 1
    q1 = peval(Q, 0, 1)
    q2 = peval(Q, 0, 2)
    okq, whyq = positive_for_b_at_least_one(QPRIME)
    ck("q(1) = -6 < 0 < 72 = q(2) and q' > 0 for beta >= 1, so beta_* "
       "is the unique root in (1,2) and q > 0 beyond it",
       q1 == -6 and q2 == 72 and okq,
       "q(1) = %s, q(2) = %s ; q': %s" % (q1, q2, whyq))

    # --- 7. the printed digits of beta_*.  q is increasing on [1,2], so
    #        q(1.1186397525) < 0 < q(1.1186397526) says exactly that the
    #        printed truncation of beta_* is correct to all ten decimals
    #        shown.  All of this is exact rational arithmetic.
    lo = Fraction(11186397525, 10 ** 10)
    hi = Fraction(11186397526, 10 ** 10)
    qlo, qhi = peval(Q, 0, lo), peval(Q, 0, hi)
    a, b = Fraction(1), Fraction(2)
    for _ in range(90):
        m = (a + b) / 2
        if peval(Q, 0, m) < 0:
            a = m
        else:
            b = m
    ck("beta_* = 1.1186397525... : q changes sign across the printed "
       "truncation", qlo < 0 < qhi and lo < (a + b) / 2 < hi,
       "q(1.1186397525) = %s < 0 < %s = q(1.1186397526); exact bisection "
       "gives beta_* = %s"
       % (sci(dfrac(qlo)), sci(dfrac(qhi)),
          fmt(dfrac((a + b) / 2), 20)))

    # --- 8. beta = 2 endpoint value
    d2pi_2 = peval(d2pi, 0, 2)
    ck("Delta_2(pi) = -3/16 exactly at beta = 2",
       d2pi_2 == Fraction(-3, 16), "value = %s" % d2pi_2)

    # --- 9. Corollary: 0 < theta_0 < pi, and the two exact cosines
    #        cos(3f) = 4c^3-3c and cos(2f) = 2c^2-1 at c = 7/8.
    c = Fraction(7, 8)
    cos3 = 4 * c ** 3 - 3 * c
    cos2 = 2 * c ** 2 - 1
    ok = (c < 1 and c ** 2 > Fraction(3, 4)
          and cos3 == Fraction(7, 128) and cos2 == Fraction(17, 32))
    ck("Corollary: 7/8 in (sqrt3/2, 1) so 0 < theta_0 < pi, with "
       "cos(theta_0/2) = 7/128 and cos(theta_0/3) = 17/32",
       ok, "(7/8)^2 = 49/64 > 48/64 = 3/4 ; cos3f = %s ; cos2f = %s"
       % (cos3, cos2))

    # --- 10. the Corollary's rational witness
    zc3 = (1 - Fraction(17, 32)) / 2
    zc2 = (1 - Fraction(7, 128)) / 2
    w = (eval_fr([peval(t, 0, 2) for t in P3C], zc3)
         - eval_fr([peval(t, 0, 2) for t in P2C], zc2))
    ck("Corollary witness  P_3(17/32) - P_2(7/128) = -6785/65536 < 0",
       w == Fraction(-6785, 65536) and w < 0, "value = %s" % w)

    # --- 11. the n = 1 margin is positive at the witness node too
    u0 = 1 - Fraction(7, 128)
    m0 = peval(pmul(U, affine), u0, 2)
    ck("n=1 margin is positive at the witness node",
       m0 > 0, "margin = %s" % m0)

    # --- 12. the Remark's four endpoint values for general alpha,
    #         from the normalized 2F1, plus consistency with -q/384.
    pairs = [((Fraction(-1, 2), Fraction(1)), Fraction(-49, 320)),
             ((Fraction(1, 4), Fraction(3, 2)), Fraction(-167, 37440)),
             ((Fraction(1, 2), Fraction(2)), Fraction(-151, 6720)),
             ((Fraction(1), Fraction(3)), Fraction(-5, 96))]
    got = []
    ok = True
    for (al, be), want in pairs:
        v = (eval_fr(tilde_zcoeffs(3, al, be), Fraction(1, 4))
             - eval_fr(tilde_zcoeffs(2, al, be), Fraction(1, 2)))
        got.append(v)
        ok = ok and v == want and v < 0
    for be in (Fraction(2), Fraction(7, 3), Fraction(5)):
        v = (eval_fr(tilde_zcoeffs(3, 0, be), Fraction(1, 4))
             - eval_fr(tilde_zcoeffs(2, 0, be), Fraction(1, 2)))
        ok = ok and v == peval(d2pi, 0, be)
    ck("Remark: the four general-alpha endpoint values are "
       "-49/320, -167/37440, -151/6720, -5/96, all negative", ok,
       "computed " + ", ".join(str(v) for v in got))


# ----------------------------------------------------------------------
# (B) Numerical layer: pi, cos, sin, acos in Decimal, from series.
# ----------------------------------------------------------------------

def atan_recip(n):
    """atan(1/n) for integer n >= 2 by its alternating Taylor series.

    For n >= 2 the terms decrease monotonically from the start, so the
    truncation error is at most the first omitted term, < EPS_SERIES.
    """
    x = Decimal(1) / Decimal(n)
    x2 = x * x
    p = x
    total = Decimal(0)
    k = 0
    while True:
        t = p / (2 * k + 1)
        if abs(t) < EPS_SERIES:
            break
        total = total + (t if k % 2 == 0 else -t)
        p = p * x2
        k += 1
        if k > 5000:
            raise RuntimeError("atan series failed to converge")
    return total


# Machin's formula: pi = 16*atan(1/5) - 4*atan(1/239).
PI = 16 * atan_recip(5) - 4 * atan_recip(239)


def dcos(x):
    """cos(x) for |x| <= 4 by the alternating Taylor series.

    The term ratio is x^2/((2k-1)(2k)), which is < 1 once
    (2k-1)(2k) > 16, i.e. from k = 3 on.  The guard k >= 3 therefore
    makes the stopping rule sound: past it the terms decrease
    monotonically and the truncation error is bounded by the first
    omitted term, < EPS_SERIES.  Every argument used below lies in
    [0, pi].
    """
    x2 = x * x
    term = Decimal(1)
    total = Decimal(1)
    k = 1
    while True:
        term = -term * x2 / ((2 * k - 1) * (2 * k))
        if k >= 3 and abs(term) < EPS_SERIES:
            break
        total = total + term
        k += 1
        if k > 2000:
            raise RuntimeError("cos series failed to converge")
    return total


def dsin(x):
    """sin(x) for |x| <= 4; same stopping argument as dcos."""
    x2 = x * x
    term = x
    total = x
    k = 1
    while True:
        term = -term * x2 / ((2 * k) * (2 * k + 1))
        if k >= 3 and abs(term) < EPS_SERIES:
            break
        total = total + term
        k += 1
        if k > 2000:
            raise RuntimeError("sin series failed to converge")
    return total


def dacos(x):
    """acos(x) for -1 < x < 1.

    cos is strictly decreasing on [0, pi], so bisection there is
    unconditionally correct; it is run down to a width of 1e-40 and then
    polished by three Newton steps t <- t + (cos t - x)/sin t, each of
    which squares the error, so the result is limited only by the
    working precision.  Callers check the result against exact rational
    values, so this routine is not trusted blindly.
    """
    lo, hi = Decimal(0), PI
    tol = Decimal(1).scaleb(-40)
    it = 0
    while hi - lo > tol and it < 400:
        mid = (lo + hi) / 2
        if dcos(mid) > x:
            lo = mid
        else:
            hi = mid
        it += 1
    t = (lo + hi) / 2
    for _ in range(3):
        s = dsin(t)
        if s == 0:
            break
        t = t + (dcos(t) - x) / s
    return t


def dfrac(fr):
    fr = Fraction(fr)
    return Decimal(fr.numerator) / Decimal(fr.denominator)


def fmt(x, nd=12):
    """Fixed-point display string (display only, never a comparison)."""
    return str(x.quantize(Decimal(1).scaleb(-nd)))


def sci(x, sig=4):
    """Short scientific display string (display only)."""
    if x == 0:
        return "0"
    e = x.adjusted()
    q = x.scaleb(-e).quantize(Decimal(1).scaleb(-(sig - 1)))
    return "%sE%+d" % (q, e)


# The sweep is at beta = 2.  Every coefficient below is obtained from
# the hypergeometric definition, converted exactly to Decimal.
BETA_NUM = Fraction(2)
C1D = [dfrac(peval(t, 0, BETA_NUM)) for t in P1C]
C2D = [dfrac(peval(t, 0, BETA_NUM)) for t in P2C]
C3D = [dfrac(peval(t, 0, BETA_NUM)) for t in P3C]
AFF_D = [dfrac(v) for v in u_coeffs_at_beta(AFFINE, BETA_NUM)]


def horner_d(coeffs, z):
    acc = Decimal(0)
    for c in reversed(coeffs):
        acc = acc * z + c
    return acc


def z_of(cosv):
    """z = (1 - x)/2 for x = cos(...)."""
    return (Decimal(1) - cosv) / 2


def delta1_num(theta):
    """P_2^{(0,2)}(cos(theta/2)) - P_1^{(0,2)}(cos theta).

    Both cosines are evaluated independently from the series; no
    double-angle identity is assumed here.
    """
    return (horner_d(C2D, z_of(dcos(theta / 2)))
            - horner_d(C1D, z_of(dcos(theta))))


def delta1_closed(theta):
    """The same quantity via the paper's closed form u*(affine in u)."""
    u = Decimal(1) - dcos(theta / 2)
    return u * horner_d(AFF_D, u)


def delta2_num(theta):
    """Delta_2 = P_3^{(0,2)}(cos(theta/3)) - P_2^{(0,2)}(cos(theta/2))."""
    return (horner_d(C3D, z_of(dcos(theta / 3)))
            - horner_d(C2D, z_of(dcos(theta / 2))))


def numeric_checks():
    tolk = Decimal(1).scaleb(-50)

    # --- 13. the transcendental kernel against known values
    pi_ref = Decimal("3.1415926535897932384626433832795028841971693993751"
                     "05820974944592307816406286208998628034825342117068")
    k_ok = abs(PI - pi_ref) < tolk
    k_ok = k_ok and dcos(Decimal(0)) == 1
    k_ok = k_ok and abs(dcos(PI / 3) - Decimal("0.5")) < tolk
    k_ok = k_ok and abs(dcos(PI / 2)) < tolk
    k_ok = k_ok and abs(dcos(PI) + 1) < tolk
    k_ok = k_ok and abs(dsin(PI / 6) - Decimal("0.5")) < tolk
    for s in ("0.3", "1.7", "3.0"):
        t = Decimal(s)
        k_ok = k_ok and abs(dcos(t) ** 2 + dsin(t) ** 2 - 1) < tolk
    ck("Decimal kernel: pi matches its known digits, and cos/sin give "
       "cos0=1, cos(pi/3)=1/2, cos(pi/2)=0, cos(pi)=-1, sin(pi/6)=1/2 "
       "and cos^2+sin^2=1, all to 1e-50", k_ok,
       "pi discrepancy %s" % sci(abs(PI - pi_ref)))

    # --- 14. acos inverts cos
    a_ok = True
    for s in ("0.4", "1.2", "2.4"):
        t = Decimal(s)
        a_ok = a_ok and abs(dacos(dcos(t)) - t) < Decimal(1).scaleb(-45)
    ck("Decimal acos inverts cos to 1e-45 at three interior points", a_ok)

    # --- 15. the numeric layer against the exact layer
    d2_pi = delta2_num(PI)
    theta0 = 6 * dacos(Decimal(7) / Decimal(8))
    c_half = dcos(theta0 / 2)
    c_third = dcos(theta0 / 3)
    d2_theta0 = delta2_num(theta0)
    tolb = Decimal(1).scaleb(-45)
    b_ok = (abs(d2_pi - dfrac(Fraction(-3, 16))) < tolb
            and abs(c_half - dfrac(Fraction(7, 128))) < tolb
            and abs(c_third - dfrac(Fraction(17, 32))) < tolb
            and abs(d2_theta0 - dfrac(Fraction(-6785, 65536))) < tolb
            and Decimal(0) < theta0 < PI)
    ck("numeric layer reproduces the exact layer: Delta_2(pi) = -3/16, "
       "cos(theta_0/2) = 7/128, cos(theta_0/3) = 17/32, "
       "Delta_2(theta_0) = -6785/65536, and 0 < theta_0 < pi (to 1e-45)",
       b_ok, "theta_0 = %s, Delta_2(theta_0) = %s"
       % (fmt(theta0), fmt(d2_theta0)))

    # --- the 1000-node grid of the paragraph under test
    n_grid = 1000
    nodes, d1, d2 = [], [], []
    for nu in range(1, n_grid + 1):
        th = PI * nu / 1001
        nodes.append(th)
        d1.append(delta1_num(th))
        d2.append(delta2_num(th))
    unamb = Decimal(1).scaleb(-30)

    # --- 16. n = 1 holds at every node
    npos = sum(1 for v in d1 if v > 0)
    clear1 = all(abs(v) > unamb for v in d1)
    urange = all(Decimal(0) < (Decimal(1) - dcos(th / 2)) < Decimal(1)
                 for th in nodes)
    ck("all 1000 nodes theta_nu = nu*pi/1001 satisfy the comparison at "
       "n = 1", npos == n_grid and clear1 and urange,
       "%d of %d have Delta_1 > 0; every |Delta_1| exceeds 1e-30, so no "
       "sign is ambiguous at %d-digit precision; and 0 < u < 1 at every "
       "node" % (npos, n_grid, PREC))

    # --- 17. the smallest n = 1 margin.  Tolerance 5e-8 is half a unit
    #         in the last significant digit of the printed 3.7e-6.
    mmin = min(d1)
    argmin = d1.index(mmin) + 1
    printed_min = Decimal("3.7E-6")
    incr = all(d1[i] < d1[i + 1] for i in range(n_grid - 1))
    ck("the smallest n=1 margin on the grid is 3.7e-6, attained at nu = 1, "
       "the margin being strictly increasing in nu (tolerance 5e-8)",
       abs(mmin - printed_min) <= Decimal("5E-8") and argmin == 1 and incr,
       "derived %s at nu = %d; discrepancy from 3.7e-6 is %s"
       % (sci(mmin, 8), argmin, sci(abs(mmin - printed_min))))

    # --- 18. the definition-based margin equals the paper's closed form
    worst = max(abs(delta1_closed(th) - v) for th, v in zip(nodes, d1))
    ck("the margin computed from the polynomial definitions agrees with "
       "the closed form u*(b+1+(b^2-b-4)u/8) at all 1000 nodes to 1e-40",
       worst < Decimal(1).scaleb(-40),
       "largest discrepancy %s" % sci(worst))

    # --- 19. how many nodes fail at n = 2
    viol = [i + 1 for i, v in enumerate(d2) if v < 0]
    clear2 = all(abs(v) > unamb for v in d2)
    ck("exactly 85 of the 1000 nodes violate the comparison at n = 2",
       len(viol) == 85 and clear2,
       "%d violators; every |Delta_2| exceeds 1e-30, so no node is a "
       "borderline zero" % len(viol))

    # --- 20. which nodes they are
    ck("the first violator is nu = 916, and the violating nodes are "
       "exactly nu = 916, ..., 1000",
       viol == list(range(916, 1001)),
       "first %s, last %s, count %d"
       % (viol[0] if viol else "none", viol[-1] if viol else "none",
          len(viol)))

    # --- 21. where Delta_2 changes sign.  The paper's two printed forms,
    #         2.8725 and 0.9144*pi = 2.87266..., are themselves only
    #         mutually consistent to about 2e-4, so 5e-4 is the stated
    #         tolerance for both.
    lo, hi = nodes[914], nodes[915]          # theta_915 and theta_916
    bracket = delta2_num(lo) > 0 > delta2_num(hi)
    for _ in range(200):
        if hi - lo <= Decimal(1).scaleb(-30):
            break
        mid = (lo + hi) / 2
        if delta2_num(mid) > 0:
            lo = mid
        else:
            hi = mid
    root = (lo + hi) / 2
    tol_r = Decimal("5E-4")
    ck("Delta_2 changes sign at theta = 2.8725 = 0.9144*pi, strictly "
       "between theta_915 and theta_916 (tolerance 5e-4 on both printed "
       "forms)",
       bracket and abs(root - Decimal("2.8725")) <= tol_r
       and abs(root / PI - Decimal("0.9144")) <= tol_r,
       "derived root %s = %s*pi" % (fmt(root, 10), fmt(root / PI, 10)))

    # --- 22. the three printed four-decimal values.  Tolerance 1e-4,
    #         i.e. twice the half-unit of the fourth decimal place, so a
    #         correctly rounded printed value passes with margin.
    tol_p = Decimal("1E-4")
    trip = [(916, Decimal("-0.0013")), (966, Decimal("-0.1032")),
            (1000, Decimal("-0.1849"))]
    ok = all(abs(d2[nu - 1] - want) <= tol_p for nu, want in trip)
    ck("the three printed values Delta_2(theta_916) = -0.0013, "
       "Delta_2(theta_966) = -0.1032, Delta_2(theta_1000) = -0.1849 "
       "(tolerance 1e-4)", ok,
       "; ".join("nu=%d derived %s" % (nu, fmt(d2[nu - 1], 8))
                 for nu, _ in trip))


NOT_RERUN = (
    "NOT RE-RUN HERE: the two hypergeometric representations quoted from "
    "the standard theory, P_n^{(0,beta)}(1-2z) = 2F1(-n, n+beta+1; 1; z) "
    "and its general-alpha analogue, are taken here as the definitions "
    "and are not re-derived from any other normalization of the Jacobi "
    "polynomials, so the triviality of the alpha = 0 normalization -- "
    "P_n^{(0,beta)}(1) = 1, which the proof uses to drop the tilde, and "
    "likewise Ptilde_n^{(alpha,beta)}(1) = 1 in the Remark -- is "
    "inherited from the z = 0 value of those series rather than tested; "
    "the elementary trigonometric identities of the Corollary "
    "(cos 3phi = 4c^3-3c, cos 2phi = 2c^2-1, and the monotonicity of cos "
    "that turns 7/8 > sqrt3/2 into 0 < phi < pi/6) are imported in the "
    "same way, corroborated numerically at beta = 2 to 1e-45 but not "
    "proved here; every bibliographic and literature-history claim is "
    "unchecked -- the passage quoted from Gautschi and Leopardi and its "
    "placement after their Conjecture 4, the fidelity to the source of "
    "the statements of Conjecture 4 and of Conjecture 3 themselves, "
    "including the identification of their (3.6) with the displayed "
    "comparison and of their (3.7) with the two intervals, the "
    "description of their "
    "verification grid (alpha = 0.9, 0.8, ..., -0.9 with beta = 20, 19, "
    "..., 0 and then negative values down to -0.999) and hence the "
    "inference that (alpha,beta) = (0,2) is one of its nodes, the "
    "attribution of the theta-grid theta_nu = nu*pi/1001, nu = 1..1000, "
    "to their routine for 0 < theta < pi with N = 1000 (the arithmetic "
    "on that grid is reproduced here, its provenance is not), the claim "
    "that the curve described in that passage runs from (-1,0) to (1,-1) "
    "and so nowhere exceeds beta = 0 and that the node therefore lies "
    "well inside the region reported positive, the contents of the Matlab "
    "verification core of their Section 2.2 (the line setting "
    "th1 = acos(-1/(2a+1)) and the commented-out th1 = pi) together with "
    "the circumstantial explanation drawn from it and the statement that "
    "the script actually used for Conjecture 4 is not printed, "
    "Koumandos's theorem at (1/2,1/2), (1/2,-1/2) and (-1/2,1/2), the "
    "history of Conjecture 3 including its disproof on a lens-shaped "
    "region and the later revised and extended parameter domains, the "
    "assertion that none of those papers revisits the polynomial-value "
    "implication, the remark that recent work distinguishes this problem "
    "from the corresponding questions for Jacobi zeros, the priority "
    "claim that no counterexample to the implication has previously been "
    "recorded, and the reference list itself (author attributions, "
    "titles, journal volumes, pages, DOIs and the arXiv identifier); "
    "also unchecked: the continuity step of Theorem 1 and of "
    "the Remark, which turns the endpoint sign Delta_2(pi) < 0 into the "
    "existence of epsilon_beta > 0 and into failure near theta = pi (the "
    "endpoint identity and its sign are verified for all beta, but the "
    "punctured neighbourhood is exhibited only at beta = 2, by the "
    "Corollary's interior witness and the 1000-node grid), the "
    "elementary range statement 0 < u < 1 for 0 < theta < pi (verified "
    "only at the 1000 grid nodes), the well-definedness of "
    "Theta_1^{(alpha,beta)} in (0,pi) as the point where "
    "P_1^{(alpha,beta)}(cos Theta_1) = 0 and of cos Theta_n^{(alpha,beta)} "
    "as the largest zero of P_n^{(alpha,beta)}, the Remark's "
    "general-alpha endpoint value, which is computed only at the four "
    "exhibited rational pairs (-1/2,1), (1/4,3/2), (1/2,2), (1,3) and "
    "not as an identity in (alpha,beta), so that neither the inference "
    "that those pairs exhibit a two-dimensional region of failure near "
    "theta = pi nor the non-minimality of the alpha = 0 family is "
    "established here, the behaviour of Delta_2 between consecutive grid "
    "nodes -- the count of 85 violators, the contiguity of the block "
    "nu = 916, ..., 1000 and the singleness of the crossing near "
    "theta = 2.8725 are facts about the 1000 samples and about one "
    "bracketed root, not exhaustive statements on (0,pi) -- any "
    "assertion about n >= 4, about the "
    "0 < theta < Theta_1 branch of Conjecture 4, about Conjecture 5, or "
    "about zeros of Jacobi polynomials, and the sweep itself at any "
    "beta other than 2 or on any grid other than "
    "theta_nu = nu*pi/1001, nu = 1..1000.")


def main():
    exact_checks()
    numeric_checks()
    n = len(CHECKS)
    failed = [name for name, ok in CHECKS if not ok]
    if failed:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), n))
        rc = 1
    else:
        print("VERDICT: ALL %d CHECKS PASS" % n)
        rc = 0
    print(NOT_RERUN)
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Independent verification of the counterexample to the half-shift clause of
the Eisenstein--Kronecker formula for n_4 (Proposition 3.1(iii) of the 2016
paper on the elliptic trilogarithm and Mahler measures of K3 surfaces).

Standard library only; no data files, no network.

===========================================================================
VALUES TAKEN FROM THE PAPER (inputs -- everything below is re-derived here)
===========================================================================
  tau0        = (2 + sqrt(-6))/4 = 1/2 + i*sqrt(6)/4          (the CM point)
  c_star      = (177/200)*(1 - i)                    (rational box centre)
  delta       = 123/3125                             (cell radius bound)
  grid        = 512 x 512 midpoint cells, 4096 = arcosh quantisation
  claimed     sum_{j,k} r_{jk} = 606722604
  claimed     31988 cells with r = 0, 230156 with 8 <= r <= 5121,
              r_{0,0} = 4050, r_{255,255} = 1803, r_{511,511} = 4080,
              sum_k r_{0,k} = 1687754, no 4096*arcosh(X) within 1e-6 of Z
  claimed     n_4(s_4(tau0)) >= 151680651/67108864 > 2.26
  claimed     Re c(tau0) =  0.8850840762662..., Im c(tau0) = -0.8850840762662...
  claimed     0.88 < Re c < 0.89 and -0.89 < Im c < -0.88
  claimed     |c - c_star| < sqrt(2)/200,  |c| < 63/50,  |d a_c| <= 263/100
  claimed     263*pi/25600 + sqrt(2)/200 < delta   (via pi<355/113, sqrt2<99/70)
  claimed     E_4(tau0) = 2.2251552127... < 2.22516, tail at M=20 below 1e-34
  claimed     exponent matrix of the torus substitution has determinant 16
  definitions of eta, A, c, s_4, F, T_d, E_4, a_c, J, U_j and the two
  truncation bounds (eta tail, U_j tail)

===========================================================================
DERIVED HERE (nothing in this list is read back from the paper)
===========================================================================
  * tau0 decoded exactly (rational arithmetic on sqrt(6)); it satisfies the
    hypothesis "tau = 1/2 + iy, y > 1/2" of the half-shift clause, and
    Im tau0 < 1/sqrt(2) so the other clause of the proposition cannot apply.
  * tau0 = (tau-1)/(2*tau) for tau = 2i/sqrt(6), the shape in which the
    proposition is used elsewhere, checked exactly.
  * c(tau0) recomputed from the eta products, with the paper's tail bound;
    the box, the distance to c_star and the modulus bound re-derived.
  * the delta inequality and the derivative bound re-derived in exact
    rational arithmetic, with a rigorous rational bracket for pi.
  * the cell-covering hypothesis sup|a_c - a_{c_star}(theta_j,theta_k)| < delta
    resampled inside cells, and the composed per-cell inequality
    min_cell J(a_c) >= r_{jk}/4096 that the certificate actually charges.
  * the from-scratch 60-digit stack (pi, sin, cos, sqrt, log) pinned against
    the rational bracket and against binary64 before it is used as a witness.
  * the FULL 512x512 census recomputed: sum of r, the cell statistics, and
    the minimal distance of 4096*arcosh(X) to an integer; a subset of cells
    is recomputed at 60 decimal digits with a from-scratch Decimal
    pi/sin/cos/sqrt/log stack, to show the binary64 stage is not deciding.
  * the exact rational lower bound 4*sum(r)/(512^2*4096) and its reduction.
  * Lemma 2.1: the arcosh contribution compared with the exact roots of the
    quadratic in w, the substitution determinant, fourth-root invariance,
    and a direct 3-variable torus integration of log|x^4+y^4+z^4+1+cxyz|.
  * E_4(tau0) from the U_j series with its tail bound, cross-checked by a
    direct truncation of the defining double lattice sum, plus the m=0 row
    value 6*zeta(4) and the Im((jm tau+n)^-3/m) identity.
  * THE CONCLUSION: n_4 lower bound > E_4 upper bound, computed.
  * a line scan showing both sides agree to <1e-5 where the identity is true
    and disagree by >1e-3 at tau0 and its neighbours, so the failure cannot
    be an artefact of these implementations.
  * the census margin measured against a binary64 error budget: this program
    has no interval arm, so the claim that the floors are insensitive to the
    binary64 rounding is made explicit as a budget on floor(4096*arcosh(X)).
  * one scan point on the half-shift line Re tau = 1/2 with y < 1/sqrt2 at
    which the identity HOLDS, which is narrower than the blanket statement the
    paper quotes from earlier work; both are stated in the closing notes.
"""
import sys
import math
import cmath
import random
from fractions import Fraction
from decimal import Decimal, getcontext

CHECKS = []

# sampling parameters of the two corroborating cell-hypothesis checks, recorded
# by those checks themselves so the closing scope note prints derived numbers
# rather than transcribed ones.
SCOPE = {}


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    return bool(ok)


# --------------------------------------------------------------------------
# inputs quoted from the paper (single place; nothing else hard-codes them)
# --------------------------------------------------------------------------
P_TAU0_RE_NUM = 2                  # tau0 = (P_TAU0_RE_NUM + sqrt(-P_TAU0_RAD))
P_TAU0_RAD = 6                     #        / P_TAU0_DEN
P_TAU0_DEN = 4
# the same point as printed a second time in the paper, tau0 = 1/2 + i sqrt6/4;
# kept as separate integers so the two transcriptions can be compared
P_TAU0_RE2 = Fraction(1, 2)
P_TAU0_IM2_RAD = 6
P_TAU0_IM2_DEN = 4
P_NGRID = 512                      # cells per direction
P_QUANT = 4096                     # arcosh quantisation denominator
P_CSTAR = Fraction(177, 200)       # c_star = P_CSTAR*(1-i)
P_DELTA = Fraction(123, 3125)      # cell radius bound
P_SUM_R = 606722604                # claimed exact cell sum
P_ZERO_CELLS = 31988
P_NONZERO_CELLS = 230156
P_RMIN_NONZERO, P_RMAX = 8, 5121
P_R00, P_R255, P_R511 = 4050, 1803, 4080
P_ROW0 = 1687754
P_MAHLER_LB = Fraction(151680651, 67108864)
P_N4_ROUND = Fraction(113, 50)     # "> 2.26"
P_EK_UB = Fraction(222516, 100000)  # "E_4(tau0) < 2.22516"
P_RE_C = "0.8850840762662"         # printed digits of Re c(tau0)
P_EK_DIGITS = "2.2251552127"       # printed digits of E_4(tau0)
P_DERIV_BOUND = Fraction(263, 100)
P_ABS_C_BOUND = Fraction(63, 50)
P_ETA_N = 30                       # truncation used in the eta products
P_U_M = 20                         # truncation used in the U_j series

# tau0 decoded from the three integers above: exact real part, exact square of
# the imaginary part, and the binary64 point used by the numerical stages.
TAU0_RE = Fraction(P_TAU0_RE_NUM, P_TAU0_DEN)
TAU0_IMSQ = Fraction(P_TAU0_RAD, P_TAU0_DEN ** 2)
SQRT6 = math.sqrt(float(P_TAU0_RAD))
TAU0 = complex(float(TAU0_RE), SQRT6 / P_TAU0_DEN)
CSTAR = complex(float(P_CSTAR), -float(P_CSTAR))
DELTA = float(P_DELTA)


def eta(tau, nterms=60):
    """Dedekind eta via its product, truncated after nterms factors."""
    q = cmath.exp(2j * math.pi * tau)
    prod = 1.0 + 0.0j
    qn = 1.0 + 0.0j
    for _ in range(nterms):
        qn *= q
        prod *= (1.0 - qn)
    return cmath.exp(1j * math.pi * tau / 12.0) * prod


def c_of_tau(tau, nterms=60):
    """c(tau) = (eta(2t)/eta(t))^6 * (16 A^4 + A^-4), A = eta(t)eta(4t)^2/eta(2t)^3."""
    e1 = eta(tau, nterms)
    e2 = eta(2 * tau, nterms)
    e4 = eta(4 * tau, nterms)
    a = e1 * e4 ** 2 / e2 ** 3
    return (e2 / e1) ** 6 * (16 * a ** 4 + a ** -4)


def jfun(a):
    """J(a) = arcosh((|a-2|+|a+2|)/4), the nonnegative branch."""
    x = (abs(a - 2.0) + abs(a + 2.0)) / 4.0
    return math.acosh(x) if x > 1.0 else 0.0


def a_c(theta, phi, c):
    """a_c(theta,phi) = e^{i(theta+phi)/2}(e^{i theta}+e^{i phi}+c)."""
    return cmath.exp(0.5j * (theta + phi)) * (
        cmath.exp(1j * theta) + cmath.exp(1j * phi) + c)


def n4_jensen(c, ngrid=384):
    """Midpoint value of the Jensen double integral, i.e. of n_4(c^4)."""
    th = [2.0 * math.pi * (j + 0.5) / ngrid for j in range(ngrid)]
    ex = [cmath.exp(1j * t) for t in th]
    hx = [cmath.exp(0.5j * t) for t in th]
    total = 0.0
    for j in range(ngrid):
        ej = ex[j]
        hj = hx[j]
        for k in range(ngrid):
            total += jfun(hj * hx[k] * (ej + ex[k] + c))
    return 4.0 * total / (ngrid * ngrid)


def ffun(z):
    """F(z) = 4 Re(z)^2/|z|^6 - 1/|z|^4."""
    n2 = z.real * z.real + z.imag * z.imag
    return 4.0 * z.real * z.real / n2 ** 3 - 1.0 / n2 ** 2


def cos_over_sin3(z):
    """cos z / sin^3 z, written in e^{2iz} so it never overflows for Im z > 0."""
    u = cmath.exp(2j * z)
    return -4j * u * (1.0 + u) / (u - 1.0) ** 3


def u_series(j, tau, mterms):
    """U_j(tau) = 2 pi^3 sum_{m>=1} cos(j pi m tau)/(m sin^3(j pi m tau))."""
    s = 0.0 + 0.0j
    for m in range(1, mterms + 1):
        s += cos_over_sin3(j * math.pi * m * tau) / m
    return 2.0 * math.pi ** 3 * s


def ek_series(tau, mterms=P_U_M):
    """E_4(tau) = Im(2 pi tau + (10/pi^3)(U_1 - 2 U_2))."""
    val = 2.0 * math.pi * tau + (10.0 / math.pi ** 3) * (
        u_series(1, tau, mterms) - 2.0 * u_series(2, tau, mterms))
    return val.imag


def census():
    """Recompute the whole 512x512 certificate.

    For each cell centre (theta_j,theta_k): l- = |a-2|, l+ = |a+2| at c_star,
    X = max{1, ((l- - delta)_+ + (l+ - delta)_+)/4}, and r = largest integer
    with cosh(r/4096) <= X.  Returns the sum, the statistics, and the minimal
    distance from 4096*arcosh(X) to an integer over the cells with r > 0.
    """
    n = P_NGRID
    th = [2.0 * math.pi * (j + 0.5) / n for j in range(n)]
    ex = [cmath.exp(1j * t) for t in th]
    hx = [cmath.exp(0.5j * t) for t in th]
    total = 0
    zeros = 0
    rmin = None
    rmax = -1
    row0 = 0
    margin = 1.0
    margin_at = None
    named = {}
    for j in range(n):
        ej = ex[j]
        hj = hx[j]
        for k in range(n):
            a = hj * hx[k] * (ej + ex[k] + CSTAR)
            lo_m = abs(a - 2.0) - DELTA
            lo_p = abs(a + 2.0) - DELTA
            s = (lo_m if lo_m > 0.0 else 0.0) + (lo_p if lo_p > 0.0 else 0.0)
            x = s / 4.0
            if x <= 1.0:
                r = 0
                zeros += 1
            else:
                v = P_QUANT * math.acosh(x)
                r = int(math.floor(v))
                gap = min(v - r, r + 1.0 - v)
                if gap < margin:
                    margin = gap
                    margin_at = (j, k)
                if rmin is None or r < rmin:
                    rmin = r
                if r > rmax:
                    rmax = r
            total += r
            if j == 0:
                row0 += r
            if (j, k) in ((0, 0), (255, 255), (511, 511)):
                named[(j, k)] = r
    return {"sum": total, "zeros": zeros, "rmin": rmin, "rmax": rmax,
            "row0": row0, "named": named, "margin": margin,
            "margin_at": margin_at}


def atan_inv_bracket(n, terms):
    """Rational bracket for atan(1/n) from its alternating Taylor series."""
    x = Fraction(1, n)
    x2 = x * x
    s = Fraction(0)
    term = x
    partials = []
    for k in range(terms):
        s += term / (2 * k + 1) if k % 2 == 0 else -term / (2 * k + 1)
        partials.append(s)
        term *= x2
    # consecutive partial sums of an alternating series bracket the limit
    return min(partials[-1], partials[-2]), max(partials[-1], partials[-2])


def sqrt_bracket(n, digits=30):
    """Rational bracket for sqrt(n) from an exact integer square root."""
    scale = 10 ** digits
    root = math.isqrt(n * scale * scale)
    return Fraction(root, scale), Fraction(root + 1, scale)


def pi_bracket(terms=40):
    """Rational lower/upper bounds for pi from Machin's formula."""
    a_lo, a_hi = atan_inv_bracket(5, terms)
    b_lo, b_hi = atan_inv_bracket(239, terms)
    return 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo


getcontext().prec = 60
_DLIM = Decimal(10) ** (-(getcontext().prec + 5))


def _atan_inv_dec(n):
    t = Decimal(1) / Decimal(n)
    x2 = t * t
    s = Decimal(0)
    term = t
    k = 0
    while term > _DLIM:
        s += term / (2 * k + 1) if k % 2 == 0 else -term / (2 * k + 1)
        k += 1
        term *= x2
    return s


DPI = 16 * _atan_inv_dec(5) - 4 * _atan_inv_dec(239)


def _dec_series(z, odd):
    s = Decimal(0)
    term = z if odd else Decimal(1)
    z2 = z * z
    n = 0
    while abs(term) > _DLIM:
        s += term if n % 2 == 0 else -term
        d = (2 * n + 2) * (2 * n + 3) if odd else (2 * n + 1) * (2 * n + 2)
        n += 1
        term = term * z2 / d
    return s


def dec_sincos(x):
    """(sin x, cos x) at the current Decimal precision, from scratch."""
    half = DPI / 2
    q = int((x / half).to_integral_value(rounding="ROUND_HALF_EVEN"))
    y = x - q * half
    sn = _dec_series(y, True)
    cs = _dec_series(y, False)
    q %= 4
    if q == 0:
        return sn, cs
    if q == 1:
        return cs, -sn
    if q == 2:
        return -sn, -cs
    return -cs, sn


def r_cell_decimal(j, k):
    """r_{jk} recomputed at 60 decimal digits, independent of binary64."""
    two_pi = 2 * DPI
    tj = two_pi * (2 * j + 1) / (2 * P_NGRID)
    tk = two_pi * (2 * k + 1) / (2 * P_NGRID)
    sj, cj = dec_sincos(tj)
    sk, ck = dec_sincos(tk)
    sh, chh = dec_sincos((tj + tk) / 2)
    # w = (cj+ck+Re c*) + i(sj+sk+Im c*), a = (chh + i sh) * w
    wr = cj + ck + Decimal(P_CSTAR.numerator) / Decimal(P_CSTAR.denominator)
    wi = sj + sk - Decimal(P_CSTAR.numerator) / Decimal(P_CSTAR.denominator)
    ar = chh * wr - sh * wi
    ai = chh * wi + sh * wr
    dl = Decimal(P_DELTA.numerator) / Decimal(P_DELTA.denominator)
    s = Decimal(0)
    for shift in (Decimal(-2), Decimal(2)):
        dr = ar + shift
        modulus = (dr * dr + ai * ai).sqrt()
        if modulus > dl:
            s += modulus - dl
    x = s / 4
    if x <= 1:
        return 0, None
    acosh = (x + (x * x - 1).sqrt()).ln()
    v = P_QUANT * acosh
    r = int(v.to_integral_value(rounding="ROUND_FLOOR"))
    return r, v - r


# --------------------------------------------------------------------------
# 1. the exhibited object: decode tau0 and check it back
# --------------------------------------------------------------------------
def check_object():
    """tau0 = (2+sqrt(-6))/4 = 1/2 + i sqrt(6)/4, decoded in exact arithmetic."""
    ck("object_tau0_exact_form",
       TAU0_RE == Fraction(1, 2) and TAU0_IMSQ == Fraction(6, 16),
       "(%d + sqrt(-%d))/%d has Re = %s and (Im)^2 = %s"
       % (P_TAU0_RE_NUM, P_TAU0_RAD, P_TAU0_DEN, TAU0_RE, TAU0_IMSQ))
    # The paper prints tau0 twice, as (2+sqrt(-6))/4 and as 1/2 + i sqrt6/4.
    # Decode the second transcription from its own integers and demand exact
    # agreement with the first, and pin math.sqrt against an integer isqrt
    # bracket, so this is not a comparison of the encoding with itself.
    re2 = P_TAU0_RE2
    imsq2 = Fraction(P_TAU0_IM2_RAD, P_TAU0_IM2_DEN ** 2)
    sq_digits = 30
    lo, hi = sqrt_bracket(P_TAU0_RAD, sq_digits)
    sq_ok = (lo * lo < P_TAU0_RAD < hi * hi
             and abs(SQRT6 - float(lo)) < 1e-15)
    ck("object_tau0_float_decode",
       TAU0_RE == re2 and TAU0_IMSQ == imsq2 and sq_ok and TAU0.imag > 0.0,
       "second transcription 1/2 + i sqrt(%d)/%d gives Re %s, (Im)^2 %s; "
       "sqrt(%d) bracketed to 1e-%d and matches math.sqrt; tau0 = %.15f%+.15fi"
       % (P_TAU0_IM2_RAD, P_TAU0_IM2_DEN, re2, imsq2, P_TAU0_RAD, sq_digits,
          TAU0.real, TAU0.imag))
    # tau0 = (tau-1)/(2 tau) with tau = i*y, y = 2/sqrt(6): the image is
    # 1/2 + i/(2y), so the requirement is Re = 1/2 and 1/(4 y^2) = (Im tau0)^2.
    ysq = Fraction(4, P_TAU0_RAD)
    exact_ok = (TAU0_RE == Fraction(1, 2)) and (Fraction(1, 4) / ysq == TAU0_IMSQ)
    tau_pre = complex(0.0, 2.0 / SQRT6)
    img = (tau_pre - 1.0) / (2.0 * tau_pre)
    ck("object_tau0_is_half_shift_image", exact_ok and abs(img - TAU0) < 1e-15,
       "(tau-1)/(2tau) at tau = 2i/sqrt(6) = %.15f%+.15fi, 1/(4y^2) = %s"
       % (img.real, img.imag, Fraction(1, 4) / ysq))


# --------------------------------------------------------------------------
# 2. hypotheses of the statement being refuted
# --------------------------------------------------------------------------
def check_hypotheses():
    """Prop 3.1(iii) half-shift clause: tau = 1/2 + iy with y > 1/2."""
    ck("hypothesis_re_is_one_half", TAU0_RE == Fraction(1, 2),
       "Re tau0 = %s, so the half-shift clause tau = 1/2 + iy applies" % TAU0_RE)
    ck("hypothesis_y_gt_one_half", TAU0_IMSQ > Fraction(1, 4),
       "(Im tau0)^2 = %s > 1/4 = (1/2)^2, so y > 1/2 as the clause requires"
       % TAU0_IMSQ)
    ck("clause_isolation_y_lt_inv_sqrt2", TAU0_IMSQ < Fraction(1, 2),
       "(Im tau0)^2 = %s < 1/2, so the tau = iy clause cannot cover tau0"
       % TAU0_IMSQ)


# --------------------------------------------------------------------------
# 3. the modular parameter c(tau0) and the rational box around it
# --------------------------------------------------------------------------
def check_c_value():
    """The eta tail bound, the box (0.88,0.89)x(-0.89,-0.88), printed digits."""
    r = abs(cmath.exp(2j * math.pi * TAU0))
    worst = 0.0
    for mult in (1, 2, 4):
        rr = r ** mult
        n = P_ETA_N
        bnd = math.expm1(rr ** (n + 1) / ((1.0 - rr) * (1.0 - rr ** (n + 1))))
        worst = max(worst, bnd)
    ck("eta_tail_bound_N30", worst < 1e-40,
       "max |R_%d - 1| over tau0,2tau0,4tau0 = %.3e" % (P_ETA_N, worst))
    n_hi = 200
    c30 = c_of_tau(TAU0, P_ETA_N)
    c200 = c_of_tau(TAU0, n_hi)
    ck("eta_truncation_stable", abs(c30 - c200) < 1e-15,
       "|c(N=%d) - c(N=%d)| = %.3e" % (P_ETA_N, n_hi, abs(c30 - c200)))
    c = c200
    ck("c_in_paper_box",
       0.88 < c.real < 0.89 and -0.89 < c.imag < -0.88,
       "c(tau0) = %.15f %+.15fi" % (c.real, c.imag))
    ck("c_printed_digits",
       abs(c.real - float(P_RE_C)) < 1e-12 and abs(c.imag + float(P_RE_C)) < 1e-12,
       "Re c - %s = %.2e ; Im c + %s = %.2e"
       % (P_RE_C, c.real - float(P_RE_C), P_RE_C, c.imag + float(P_RE_C)))
    ck("c_anti_diagonal", abs(c.real + c.imag) < 1e-14,
       "Re c + Im c = %.2e (c is a positive multiple of 1-i)" % (c.real + c.imag))
    return c


def check_c_box_consequences(c):
    """|c - c_star| < sqrt2/200 and |c| < 63/50, from the box and from c."""
    # the box of eq. (c-box) must be centred on c_star tightly enough that its
    # corners lie within sqrt2/200; derive the half-width from P_CSTAR itself
    # rather than assuming it, so a wrong c_star breaks this clause too.
    box_lo, box_hi = Fraction(88, 100), Fraction(89, 100)
    half = max(box_hi - P_CSTAR, P_CSTAR - box_lo)
    ok_exact = (box_lo < P_CSTAR < box_hi
                and 2 * half * half <= Fraction(2, 200 ** 2))
    d = abs(c - CSTAR)
    ck("c_distance_to_c_star", ok_exact and d < math.sqrt(2.0) / 200.0,
       "|c-c*| = %.3e < sqrt2/200 = %.3e" % (d, math.sqrt(2.0) / 200.0))
    corner = 2 * Fraction(89, 100) ** 2         # worst |c|^2 allowed by the box
    ck("c_modulus_bound", corner < P_ABS_C_BOUND ** 2 and abs(c) < float(P_ABS_C_BOUND),
       "box corner |c|^2 <= %s < (%s)^2 = %s ; |c| = %.6f"
       % (corner, P_ABS_C_BOUND, P_ABS_C_BOUND ** 2, abs(c)))


# --------------------------------------------------------------------------
# 4. the cell geometry: derivative bound and the delta inequality
# --------------------------------------------------------------------------
def check_delta_inequality():
    """263 pi/25600 + sqrt2/200 < 123/3125, with rigorous rational bounds."""
    pi_lo, pi_hi = pi_bracket()
    ck("rational_pi_bracket",
       pi_lo < pi_hi and float(pi_lo) <= math.pi <= float(pi_hi)
       and pi_hi < Fraction(355, 113),
       "Machin bracket brackets pi, width %.2e, and pi < 355/113"
       % float(pi_hi - pi_lo))
    ck("rational_sqrt2_bound", Fraction(99, 70) ** 2 > 2,
       "(99/70)^2 = %s > 2, so sqrt2 < 99/70" % (Fraction(99, 70) ** 2))
    lhs = Fraction(263, 1) * Fraction(355, 113) / 25600 + Fraction(99, 70) / 200
    ck("delta_covers_cell_radius", lhs < P_DELTA,
       "263*(355/113)/25600 + (99/70)/200 = %.9f < delta = %.9f"
       % (float(lhs), float(P_DELTA)))


def check_derivative_bound(c):
    """|d a_c/d theta|, |d a_c/d phi| <= 2 + |c|/2 <= 263/100."""
    ck("derivative_bound_arithmetic",
       2 + P_ABS_C_BOUND / 2 == P_DERIV_BOUND,
       "2 + (%s)/2 = %s = %s"
       % (P_ABS_C_BOUND, 2 + P_ABS_C_BOUND / 2, P_DERIV_BOUND))
    worst = 0.0
    n = 240
    for j in range(n):
        t = 2.0 * math.pi * j / n
        for k in range(n):
            p = 2.0 * math.pi * k / n
            pre = cmath.exp(0.5j * (t + p))
            tot = cmath.exp(1j * t) + cmath.exp(1j * p) + c
            dt = pre * (0.5j * tot + 1j * cmath.exp(1j * t))
            dp = pre * (0.5j * tot + 1j * cmath.exp(1j * p))
            worst = max(worst, abs(dt), abs(dp))
    SCOPE["deriv_points"] = n * n
    ck("derivative_bound_sampled", worst <= float(P_DERIV_BOUND),
       "max sampled |grad a_c| = %.9f <= %.2f [sample: %d points on the torus; "
       "the bound itself is the exact triangle-inequality one checked above]"
       % (worst, float(P_DERIV_BOUND), n * n))


def check_cell_covering(c):
    """sup over each cell of |a_c(theta,phi) - a_{c*}(theta_j,theta_k)| < delta.

    This is a SAMPLE, not a proof: it strides the grid and resamples a finite
    sub-lattice inside each visited cell, so it bounds a maximum over sampled
    points, never the true supremum.  The rigorous version of the hypothesis is
    check_derivative_bound plus check_delta_inequality, both exact-rational; the
    sampling parameters below are recorded in SCOPE and disclosed in the closing
    note so this check is not mistaken for a cell-by-cell verification.
    """
    n = P_NGRID
    step = 2.0 * math.pi / n
    worst = 0.0
    sub = 9
    stride = 13
    got = 0
    for j in range(0, n, stride):
        tj = step * (j + 0.5)
        for k in range(0, n, stride):
            tk = step * (k + 0.5)
            centre = a_c(tj, tk, CSTAR)
            got += 1
            for u in range(sub):
                t = tj + step * (u / (sub - 1.0) - 0.5)
                for v in range(sub):
                    p = tk + step * (v / (sub - 1.0) - 0.5)
                    worst = max(worst, abs(a_c(t, p, c) - centre))
    SCOPE["cov_cells"] = got
    SCOPE["cov_total"] = n * n
    SCOPE["cov_stride"] = stride
    SCOPE["cov_sub"] = sub
    ck("cell_covering_within_delta", worst < DELTA,
       "max resampled deviation %.9f < delta %.9f (ratio %.4f) "
       "[sample: %d of %d cells, every %dth index, %dx%d points per cell]"
       % (worst, DELTA, worst / DELTA, got, n * n, stride, sub, sub))


def check_cell_lower_bound_end_to_end(c):
    """End to end: min over each cell of J(a_c) must be at least r_{jk}/4096.

    The certificate composes three steps -- the covering hypothesis, the
    subtraction of delta from the two moduli, and the floor.  Each is checked
    separately above; this joins them on a subsample of cells, so a reversed
    direction or a mis-ordered composition anywhere in the chain shows up as a
    cell where the true integrand dips below the quantised bound charged to it.

    Like check_cell_covering this is a SAMPLE: it strides the grid and takes the
    minimum over a finite sub-lattice inside each visited cell, so it cannot
    certify the per-cell infimum.  Sampling parameters are recorded in SCOPE and
    disclosed in the closing note.
    """
    n = P_NGRID
    step = 2.0 * math.pi / n
    sub = 7
    stride = 7
    bad = 0
    got = 0
    worst = None
    for j in range(0, n, stride):
        tj = step * (j + 0.5)
        for k in range(0, n, stride):
            tk = step * (k + 0.5)
            a = a_c(tj, tk, CSTAR)
            lm = abs(a - 2.0) - DELTA
            lp = abs(a + 2.0) - DELTA
            s = (lm if lm > 0.0 else 0.0) + (lp if lp > 0.0 else 0.0)
            x = s / 4.0
            r = 0 if x <= 1.0 else int(math.floor(P_QUANT * math.acosh(x)))
            lo = min(jfun(a_c(tj + step * (u / (sub - 1.0) - 0.5),
                              tk + step * (v / (sub - 1.0) - 0.5), c))
                     for u in range(sub) for v in range(sub))
            slack = lo - r / float(P_QUANT)
            if slack < 0.0:
                bad += 1
            if worst is None or slack < worst:
                worst = slack
            got += 1
    SCOPE["e2e_cells"] = got
    SCOPE["e2e_total"] = n * n
    SCOPE["e2e_stride"] = stride
    SCOPE["e2e_sub"] = sub
    ck("cell_lower_bound_end_to_end", bad == 0 and worst > 0.0,
       "%d cells resampled %dx%d: min_cell J(a_c) >= r/%d in all of them, "
       "worst slack %+.3e [sample: %d of %d cells, every %dth index]"
       % (got, sub, sub, P_QUANT, worst, got, n * n, stride))


# --------------------------------------------------------------------------
# 5. the certificate: the full 512x512 census and the exact rational bound
# --------------------------------------------------------------------------
def check_census():
    """Recompute the census and compare with every number the paper prints."""
    res = census()
    ck("census_sum_r", res["sum"] == P_SUM_R,
       "recomputed sum_{j,k} r_{jk} = %d (paper %d)" % (res["sum"], P_SUM_R))
    ck("census_zero_cells", res["zeros"] == P_ZERO_CELLS
       and P_NGRID ** 2 - res["zeros"] == P_NONZERO_CELLS,
       "%d cells with r=0, %d with r>0 (paper %d / %d)"
       % (res["zeros"], P_NGRID ** 2 - res["zeros"], P_ZERO_CELLS,
          P_NONZERO_CELLS))
    ck("census_r_range",
       res["rmin"] == P_RMIN_NONZERO and res["rmax"] == P_RMAX,
       "nonzero r in [%d,%d] (paper [%d,%d])"
       % (res["rmin"], res["rmax"], P_RMIN_NONZERO, P_RMAX))
    nm = res["named"]
    ck("census_named_cells",
       nm.get((0, 0)) == P_R00 and nm.get((255, 255)) == P_R255
       and nm.get((511, 511)) == P_R511 and res["row0"] == P_ROW0,
       "r_00=%s r_255,255=%s r_511,511=%s sum_k r_0k=%d"
       % (nm.get((0, 0)), nm.get((255, 255)), nm.get((511, 511)), res["row0"]))
    ck("census_floor_margin", res["margin"] > 1e-6,
       "min distance of 4096*arcosh(X) to Z over r>0 cells = %.4e at %s"
       % (res["margin"], res["margin_at"]))
    return res


def check_decimal_stack():
    """Pin the from-scratch 60-digit stack before it is used as a witness.

    Agreement of the 60-digit floors with the binary64 floors only carries
    force if the 60-digit arm is itself right, and the floor comparison is a
    blunt instrument: a relative perturbation of 1e-9 in the constant used for
    pi leaves every sampled floor unchanged.  So pin the constant inside the
    rational Machin bracket, and pin sin/cos, sqrt and log directly.
    """
    lo, hi = pi_bracket()
    ok_pi = lo < Fraction(DPI) < hi
    worst_t = 0.0
    worst_p = Decimal(0)
    for deg in range(0, 360, 11):
        x = Decimal(deg) * DPI / 180
        sn, cs = dec_sincos(x)
        worst_p = max(worst_p, abs(sn * sn + cs * cs - 1))
        worst_t = max(worst_t, abs(float(sn) - math.sin(math.radians(deg))),
                      abs(float(cs) - math.cos(math.radians(deg))))
    ck("decimal_stack_pinned",
       ok_pi and worst_t < 1e-14 and worst_p < Decimal(10) ** -55,
       "pi constant inside the rational bracket; sin/cos vs binary64 %.2e; "
       "max |sin^2+cos^2-1| %.2e" % (worst_t, float(worst_p)))
    worst_a = 0.0
    nx = 0
    for i in range(1, 40):
        xv = 1 + Decimal(i) / 7
        acosh_d = (xv + (xv * xv - 1).sqrt()).ln()
        worst_a = max(worst_a, abs(float(acosh_d) - math.acosh(float(xv))))
        nx += 1
    ck("decimal_arcosh_pinned", worst_a < 1e-13,
       "max |log(x+sqrt(x^2-1)) at %d digits - math.acosh| over %d x = %.2e"
       % (getcontext().prec, nx, worst_a))


def check_census_high_precision(res):
    """Redo a spread of cells at 60 decimal digits: the floors must agree."""
    cells = [(0, 0), (255, 255), (511, 511), (0, 1), (1, 0), (256, 0)]
    if res["margin_at"] is not None:
        cells.append(res["margin_at"])
    step = 23
    for j in range(0, P_NGRID, step):
        for k in range(0, P_NGRID, step * 3):
            cells.append((j, k))
    cells = sorted(set(cells))
    n = P_NGRID
    twopi = 2.0 * math.pi
    bad = []
    for (j, k) in cells:
        tj = twopi * (j + 0.5) / n
        tk = twopi * (k + 0.5) / n
        a = a_c(tj, tk, CSTAR)
        lm = abs(a - 2.0) - DELTA
        lp = abs(a + 2.0) - DELTA
        s = (lm if lm > 0.0 else 0.0) + (lp if lp > 0.0 else 0.0)
        x = s / 4.0
        rf = 0 if x <= 1.0 else int(math.floor(P_QUANT * math.acosh(x)))
        rd, _ = r_cell_decimal(j, k)
        if rf != rd:
            bad.append((j, k, rf, rd))
    SCOPE["hp_cells"] = len(cells)
    ck("census_binary64_not_deciding", not bad,
       "%d cells recomputed at %d digits, %d disagreements"
       % (len(cells), getcontext().prec, len(bad)))
    # The paper's own run enclosed every transcendental quantity in interval
    # arithmetic; this program does not (see the closing scope note).  What it
    # can offer instead is a budget argument on the only decision binary64 makes
    # here, namely floor(4096*arcosh(X)).  A floor can move only if the error in
    # 4096*arcosh(X) is as large as the census's own minimal distance from that
    # quantity to an integer.  Allow a generous 64 units in the last place of
    # 1.0 for X itself, and amplify by the largest value of d arcosh/dX over the
    # cells that carry r > 0: that derivative is 1/sqrt(X^2-1), largest at the
    # smallest X, and the smallest X with r > 0 is at least cosh(rmin/4096), so
    # the amplification is at most 1/sinh(rmin/4096).  Both rmin and the margin
    # come from the census, not from the paper.
    ulps = 64
    ulp_budget = ulps * 2.0 ** -52
    if res["rmin"] is None:              # no cell carries r > 0: nothing to bound
        err_budget = float("inf")
        amp = float("inf")
    else:
        amp = 1.0 / math.sinh(res["rmin"] / float(P_QUANT))
        err_budget = P_QUANT * ulp_budget * amp
    SCOPE["budget"] = err_budget
    SCOPE["budget_factor"] = res["margin"] / err_budget
    ck("floor_margin_beats_binary64_budget", res["margin"] > err_budget,
       "min distance to Z %.4e > binary64 error budget %.2e = %d * %d*2^-52 / "
       "sinh(%s/%d) (factor %.1f); a budget argument on the floors, not an "
       "enclosure" % (res["margin"], err_budget, P_QUANT, ulps, res["rmin"],
                      P_QUANT, SCOPE["budget_factor"]))


def check_exact_lower_bound(res):
    """4 sum(r)/(512^2 * 4096) reduces to 151680651/67108864 and exceeds 2.26."""
    bound = Fraction(4 * res["sum"], P_NGRID ** 2 * P_QUANT)
    ck("mahler_bound_exact_rational", bound == P_MAHLER_LB,
       "4*%d/(512^2*4096) = %s = %.10f (paper %s)"
       % (res["sum"], bound, float(bound), P_MAHLER_LB))
    ck("mahler_bound_exceeds_226", bound > P_N4_ROUND,
       "%s - %s = %s > 0" % (bound, P_N4_ROUND, bound - P_N4_ROUND))
    return bound


# --------------------------------------------------------------------------
# 6. the Jensen reduction (Lemma 2.1)
# --------------------------------------------------------------------------
def check_jensen_substitution():
    """The exponent matrix of (x^3/yz, y^3/xz, z^3/xy) has determinant 16."""
    m = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
    det = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
           - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
           + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    rows = [sum(r) for r in m]
    ck("jensen_substitution_determinant", det == 16 and rows == [1, 1, 1],
       "det = %d and row sums %s (so uvw = xyz)" % (det, rows))


def check_jensen_arcosh(c):
    """J(a_c) equals the Jensen contribution of w^2+(u+v+c)w+(uv)^-1."""
    rng = random.Random(20260101)
    worst = 0.0
    nsamp = 0
    for cc in (c, CSTAR, complex(0.3, 1.7), complex(-2.5, 0.0)):
        for _ in range(300):
            nsamp += 1
            t = rng.uniform(0.0, 2.0 * math.pi)
            p = rng.uniform(0.0, 2.0 * math.pi)
            u = cmath.exp(1j * t)
            v = cmath.exp(1j * p)
            b = u + v + cc
            const = 1.0 / (u * v)
            disc = cmath.sqrt(b * b - 4.0 * const)
            r1 = (-b + disc) / 2.0
            r2 = (-b - disc) / 2.0
            jensen = max(0.0, math.log(abs(r1))) + max(0.0, math.log(abs(r2)))
            worst = max(worst, abs(jensen - jfun(a_c(t, p, cc))))
    ck("jensen_arcosh_matches_roots", worst < 1e-10,
       "max |log max(1,|rho|) - J(a_c)| over %d samples = %.3e" % (nsamp, worst))


def m_direct_3d(c, ngrid=101):
    """4*m(x^4+y^4+z^4+1+c xyz) by midpoint average on the 3-torus.

    ngrid is taken coprime to 4 so that the fourth powers are not aliased.
    """
    th = [2.0 * math.pi * (j + 0.5) / ngrid for j in range(ngrid)]
    ex = [cmath.exp(1j * t) for t in th]
    e4 = [e ** 4 for e in ex]
    total = 0.0
    for ia in range(ngrid):
        for ib in range(ngrid):
            pre = e4[ia] + e4[ib] + 1.0
            xy = ex[ia] * ex[ib]
            for ic in range(ngrid):
                total += math.log(abs(pre + e4[ic] + c * xy * ex[ic]))
    return 4.0 * total / ngrid ** 3


def check_jensen_reduction(c):
    """The 2-D Jensen integral agrees with the original 3-variable integral."""
    two_d = n4_jensen(c, 512)
    three_d = m_direct_3d(c, 101)
    ck("jensen_reduction_2d_vs_3d", abs(two_d - three_d) < 1e-3,
       "2-D Jensen %.9f vs direct 3-D torus average %.9f, diff %.2e"
       % (two_d, three_d, two_d - three_d))
    vals = [n4_jensen(c * 1j ** k, 384) for k in range(4)]
    ck("fourth_root_independence", max(vals) - min(vals) < 1e-12,
       "n_4 under c -> i^k c: spread %.2e" % (max(vals) - min(vals)))
    return two_d


# --------------------------------------------------------------------------
# 7. the lattice sum E_4
# --------------------------------------------------------------------------
def check_ek_ingredients():
    """The m=0 row value 6 zeta(4) = pi^4/15 and the Im((jm tau+n)^-3/m) rule."""
    # F(n) = 3/n^4 for real n, so the omitted row is 6 zeta(4)
    nmax = 500
    s = 2.0 * math.fsum(ffun(complex(float(n), 0.0)) for n in range(1, nmax + 1))
    tail = 2.0 / float(nmax) ** 3        # 2*3*int_N^inf x^-4 dx, an upper bound
    target = math.pi ** 4 / 15.0
    ck("ek_m0_row_is_6zeta4", 0.0 < target - s < tail,
       "sum_{0<|n|<=%d} F(n) = %.14f, pi^4/15 = %.14f, deficit %.3e < %.3e"
       % (nmax, s, target, target - s, tail))
    rng = random.Random(7)
    worst = 0.0
    y = TAU0.imag
    for j in (1, 2):
        for _ in range(400):
            m = rng.randint(1, 30)
            n = rng.randint(-40, 40)
            z = j * m * TAU0 + n
            lhs = ((z ** -3) / m).imag
            worst = max(worst, abs(lhs + j * y * ffun(z)))
    ck("ek_row_identity", worst < 1e-12,
       "max |Im((jm tau+n)^-3/m) + j Im(tau) F| = %.3e" % worst)


def check_ek_value():
    """E_4(tau0) from the U_j series, its tail bound, and the printed digits."""
    y = TAU0.imag
    worst = 0.0
    for j in (1, 2):
        qj = math.exp(-2.0 * j * math.pi * y)
        m = P_U_M
        num = 8.0 * math.pi ** 3 * (1.0 + qj ** (m + 1)) * qj ** (m + 1)
        den = (m + 1) * (1.0 - qj) * (1.0 - qj ** (m + 1)) ** 3
        worst = max(worst, num / den)
    ck("ek_tail_bound_M20", worst < 1e-34,
       "max U_j tail at M=%d = %.3e" % (P_U_M, worst))
    m_hi = 80
    e20 = ek_series(TAU0, P_U_M)
    e80 = ek_series(TAU0, m_hi)
    ck("ek_series_stable", abs(e20 - e80) < 1e-14,
       "|E_4(M=%d) - E_4(M=%d)| = %.2e" % (P_U_M, m_hi, abs(e20 - e80)))
    ck("ek_printed_digits", abs(e80 - float(P_EK_DIGITS)) < 1e-9,
       "E_4(tau0) = %.13f (paper %s...)" % (e80, P_EK_DIGITS))
    ck("ek_below_222516", e80 < float(P_EK_UB),
       "E_4(tau0) = %.13f < %s" % (e80, float(P_EK_UB)))
    return e80


def check_ek_direct_lattice(ek_ref):
    """Independent check: truncate the defining double lattice sum itself."""
    mmax = 300
    y = TAU0.imag
    tot = []
    for j in (1, 2):
        s = 0.0
        for m in range(-mmax, mmax + 1):
            zm = j * m * TAU0
            base_r = zm.real
            base_i = zm.imag
            for n in range(-mmax, mmax + 1):
                if m == 0 and n == 0:
                    continue
                s += ffun(complex(base_r + n, base_i))
        tot.append(s)
    val = (10.0 * y / math.pi ** 3) * (-tot[0] + 4.0 * tot[1])
    ck("ek_direct_lattice_agrees", abs(val - ek_ref) < 1e-4,
       "direct sum over |m|,|n|<=%d gives %.9f vs series %.9f, diff %.2e"
       % (mmax, val, ek_ref, val - ek_ref))


# --------------------------------------------------------------------------
# 8. the conclusion, and a scan showing where the formula does hold
# --------------------------------------------------------------------------
def check_refutation(bound, ek_val, n4_true):
    """The load-bearing comparison: certified n_4 lower bound > E_4."""
    gap = bound - P_EK_UB
    ck("refutation_lower_bound_exceeds_ek", gap > 0,
       "n_4 >= %s = %.10f, E_4 < %s: gap >= %.10f"
       % (bound, float(bound), float(P_EK_UB), float(gap)))
    ck("refutation_against_computed_ek", float(bound) > ek_val,
       "certified %.10f > recomputed E_4 %.10f, margin %.6f"
       % (float(bound), ek_val, float(bound) - ek_val))
    ck("certificate_not_over_claimed", float(bound) <= n4_true,
       "certified lower bound %.9f <= independently integrated n_4 %.9f"
       % (float(bound), n4_true))
    ck("identity_fails_at_tau0", abs(n4_true - ek_val) > 1e-3,
       "n_4(s_4(tau0)) = %.9f, E_4(tau0) = %.9f, difference %+.9f"
       % (n4_true, ek_val, n4_true - ek_val))


def check_s4_anchors(c):
    """s_4 at two CM points where it is a known integer, and s_4(tau0) < 0."""
    si = (c_of_tau(1j) ** 4)
    sh = (c_of_tau(1j / math.sqrt(2.0)) ** 4)
    ck("s4_known_cm_values",
       abs(si - 648.0) < 1e-8 and abs(sh - 256.0) < 1e-8,
       "s_4(i) = %.9f (648), s_4(i/sqrt2) = %.9f (256)" % (si.real, sh.real))
    s0 = c ** 4
    ck("s4_at_tau0_negative_real", abs(s0.imag) < 1e-12 and s0.real < 0.0,
       "s_4(tau0) = %.9f %+.2ei (the s<0 branch)" % (s0.real, s0.imag))


def check_line_scan():
    """Where the proposition is true both sides agree; near tau0 they do not."""
    inv = 1.0 / math.sqrt(2.0)
    holds = [(0.0, inv), (0.0, 0.8), (0.0, 1.0), (0.0, 1.2),
             (0.5, 0.70), (0.5, 0.8), (0.5, 1.0)]
    fails = [(0.5, 0.52), (0.5, 0.55), (0.5, 0.60), (0.5, SQRT6 / 4.0),
             (0.5, 0.65), (0.5, 0.68), (0.0, 0.55), (0.0, 0.60), (0.0, 0.65)]
    worst_hold = 0.0
    detail_h = []
    diffs_h = []
    for (x, y) in holds:
        tau = complex(x, y)
        d = n4_jensen(c_of_tau(tau), 384) - ek_series(tau, 60)
        worst_hold = max(worst_hold, abs(d))
        diffs_h.append((x, y, abs(d)))
        detail_h.append("%.2f%+.4fi:%.1e" % (x, y, abs(d)))
    ck("scan_identity_holds_where_asserted_true", worst_hold < 1e-5,
       "max |n_4 - E_4| = %.2e over %s" % (worst_hold, ",".join(detail_h)))
    least_fail = None
    detail_f = []
    diffs_f = []
    for (x, y) in fails:
        tau = complex(x, y)
        d = n4_jensen(c_of_tau(tau), 384) - ek_series(tau, 60)
        least_fail = abs(d) if least_fail is None else min(least_fail, abs(d))
        diffs_f.append((x, y, abs(d)))
        detail_f.append("%.2f%+.4fi:%+.4f" % (x, y, d))
    ck("scan_identity_fails_below_threshold", least_fail > 1e-3,
       "min |n_4 - E_4| = %.2e over %s" % (least_fail, ",".join(detail_f)))
    # The paper's "Relation to earlier work" paragraph quotes an earlier status
    # report for the blanket statement that below Im tau = 1/sqrt2 the identity
    # "fails everywhere" -- evidence gathered on Re tau = 0 and Re tau = 1/4 --
    # and concedes that the blanket statement therefore already covers tau0.  On
    # the half-shift line Re tau = 1/2 the scan above contradicts the blanket
    # reading: it holds at a point with y < 1/sqrt2.  Record that as a check
    # rather than leave it implicit in a list of numbers.  If no such point is
    # in the scan the check FAILS rather than vanishing.
    sub = [(y, dv) for (x, y, dv) in diffs_h
           if x == 0.5 and Fraction(y) ** 2 < Fraction(1, 2)]
    half_fails = [(y, dv) for (x, y, dv) in diffs_f if x == 0.5]
    top_y, top_d = (max(half_fails) if half_fails
                    else (float("nan"), float("nan")))
    sub_y, sub_d = (max(sub, key=lambda t: t[1]) if sub
                    else (float("nan"), float("inf")))
    SCOPE["sub_y"], SCOPE["sub_d"] = sub_y, sub_d
    SCOPE["top_y"], SCOPE["top_d"] = top_y, top_d
    ck("scan_holds_below_inv_sqrt2_on_half_line",
       bool(sub) and bool(half_fails) and sub_d < 1e-5,
       "at tau = 1/2%+.4fi, y^2 = %.6f < 1/2 so y < 1/sqrt2, yet |n_4 - E_4| = "
       "%.2e; the largest y on Re tau = 1/2 at which this scan exhibits a "
       "failure is %.4f (diff %.1e). So the blanket statement quoted in the "
       "paper from earlier work, that the identity fails everywhere below "
       "1/sqrt2, is not true on the half-shift line; Theorem 1 is unaffected "
       "(midpoint quadratures, not enclosures)"
       % (sub_y, sub_y * sub_y, sub_d, top_y, top_d))


def main():
    print("# tau0 = 1/2 + i*sqrt(6)/4 = %.15f %+.15fi" % (TAU0.real, TAU0.imag))
    print("# grid %d^2, quantisation 1/%d, c* = %s(1-i), delta = %s"
          % (P_NGRID, P_QUANT, P_CSTAR, P_DELTA))
    check_object()
    check_hypotheses()
    c = check_c_value()
    check_c_box_consequences(c)
    check_delta_inequality()
    check_derivative_bound(c)
    check_cell_covering(c)
    check_cell_lower_bound_end_to_end(c)
    res = check_census()
    check_decimal_stack()
    check_census_high_precision(res)
    bound = check_exact_lower_bound(res)
    check_jensen_substitution()
    check_jensen_arcosh(c)
    n4_true = check_jensen_reduction(c)
    check_ek_ingredients()
    ek_val = check_ek_value()
    check_ek_direct_lattice(ek_val)
    check_s4_anchors(c)
    check_refutation(bound, ek_val, n4_true)
    check_line_scan()
    print("# NOT RE-RUN HERE: the paper's own run used a multiprecision "
          "interval library; this program reproduces the certificate in "
          "binary64 and confirms on a sample of %d of the %d cells at %d "
          "decimal digits that the floors do not depend on that. NO ARM OF "
          "THIS PROGRAM IS AN ENCLOSURE: the trusted computing base the "
          "paper's Reproducibility paragraph pins down -- an interval context "
          "at 80 decimal digits, with each binary64 result enlarged outward by "
          "nextafter -- is not exercised here, and neither is directed "
          "interval arithmetic for the per-cell lower bounds on |a_{c*} +- 2|; "
          "the %d^2 census is evaluated at binary64 point values and the "
          "%d-digit arm at Decimal point values. What stands in for those "
          "semantics is the pair census_binary64_not_deciding and "
          "floor_margin_beats_binary64_budget, the latter showing that the "
          "census's own minimal distance from 4096*arcosh(X) to an integer "
          "exceeds a 64-ulp binary64 error budget (%.2e) by a factor of %.1f. "
          "That is a budget argument about the floors, not a rigorous "
          "enclosure, so a referee wanting the interval semantics themselves "
          "must still consult the paper's own run. The values of n_4 used for "
          "the consistency and scan checks are midpoint quadratures, not "
          "enclosures; the certified bound %s is the exact rational one and is "
          "what the refutation rests on."
          % (SCOPE["hp_cells"], P_NGRID ** 2, getcontext().prec, P_NGRID,
             getcontext().prec, SCOPE["budget"], SCOPE["budget_factor"],
             bound))
    print("# NOT RE-RUN HERE (covering hypothesis, sampling): the two checks "
          "that touch sup|a_c - a_{c*}| < delta empirically are SUBSAMPLES, not "
          "cell-by-cell verifications. cell_covering_within_delta visits %d of "
          "the %d cells (every %dth index in each direction) and evaluates a "
          "%dx%d sub-lattice inside each; cell_lower_bound_end_to_end visits %d "
          "of the %d cells (every %dth index) and takes a %dx%d minimum inside "
          "each. A sampled maximum is not a supremum and a sampled minimum is "
          "not an infimum, so neither check would detect a violation confined to "
          "an unvisited cell or to a point between sample nodes; they corroborate "
          "the hypothesis, they do not establish it. What establishes it is the "
          "exact-rational pair derivative_bound_arithmetic (|grad a_c| <= 2 + "
          "|c|/2 = %s on the whole torus, from the triangle inequality) and "
          "delta_covers_cell_radius (%d*(355/113)/25600 + (99/70)/200 < %s "
          "with pi and sqrt2 bracketed rationally), which together bound the "
          "deviation over every cell without sampling; derivative_bound_sampled "
          "is likewise only a %d-point corroboration of that analytic bound. The "
          "census lower bound therefore rests on the analytic covering argument, "
          "not on these samples."
          % (SCOPE["cov_cells"], SCOPE["cov_total"], SCOPE["cov_stride"],
             SCOPE["cov_sub"], SCOPE["cov_sub"],
             SCOPE["e2e_cells"], SCOPE["e2e_total"], SCOPE["e2e_stride"],
             SCOPE["e2e_sub"], SCOPE["e2e_sub"],
             P_DERIV_BOUND, P_DERIV_BOUND.numerator, P_DELTA,
             SCOPE["deriv_points"]))
    print("# OBSERVED HERE, AND NARROWER THAN THE PAPER'S CONCESSION TO "
          "EARLIER WORK: the paper concedes that the blanket statement it "
          "quotes from an earlier status report, that below Im tau = 1/sqrt2 "
          "the identity fails everywhere, already covers tau0. This program's "
          "scan does not support that statement on the half-shift line "
          "Re tau = 1/2: at tau = 1/2%+.4fi, with y^2 = %.6f < 1/2, the two "
          "sides agree to %.2e, while the largest y on that line at which the "
          "scan exhibits a failure is %.4f (diff %.1e); the earlier evidence "
          "was gathered on Re tau = 0 and Re tau = 1/4. See "
          "scan_holds_below_inv_sqrt2_on_half_line. This makes the paper's "
          "contribution larger rather than smaller and does not touch "
          "Theorem 1, whose refutation rests on the exact rational certificate "
          "at tau0; the scan values are midpoint quadratures, not enclosures."
          % (SCOPE["sub_y"], SCOPE["sub_y"] ** 2, SCOPE["sub_d"],
             SCOPE["top_y"], SCOPE["top_d"]))


def report():
    n = len(CHECKS)
    k = sum(1 for _, ok in CHECKS if not ok)
    if k == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
    else:
        print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
    return 0 if k == 0 else 1


if __name__ == "__main__":
    main()
    sys.exit(report())

#!/usr/bin/env python3
"""verify.py -- re-derives every computational claim of

    "The Iterated Beatty Range Criterion Fails Below 5/4 and Holds Above 2"  (paper.tex)

from the objects printed in the paper and from nothing else.

Python 3.9+, STANDARD LIBRARY ONLY (fractions, math.isqrt, sys) -- no third-party package
and no external data file.  There is NO FLOATING-POINT ARITHMETIC anywhere: every modulus
alpha is carried as an exact element (a + b*sqrt(d)) of a real quadratic field with a, b
rational, every comparison between algebraic numbers is reduced to a comparison of
integers, and every decimal string printed in the paper is checked by exact rational
bracketing rather than by re-printing a float.

Objects read from the paper:
  * alpha_1 = (1+sqrt5)/2 = phi, the root of t^2-t-1 in (1.6,1.7)
  * alpha_2 = (7+sqrt2)/8,       the root of 64t^2-112t+47 in (1.05,1.06)
  * the chain 3 -> 4 -> 6 -> 9 -> 14 -> 22 of Section 5
  * the window endpoints (1055+896*sqrt2)/223 and 8+8*sqrt2 of Section 4
  * the closed forms 13phi-21, 8phi-13, 5phi-8, 3phi-4 and (9x+38-8(x-1)sqrt2)/47
  * the thresholds a_n = the root of t^(n+1)-2t^n+1 in (1,2)
  * every decimal expansion displayed in the paper

Prints one "PASS <name> [detail]" line per check and closes with
    VERDICT: ALL <n> CHECKS PASS
exiting 0 if and only if every check passed.  What it does NOT cover is printed under
"NOT RE-RUN" before the verdict.
"""

import sys
from fractions import Fraction as F
from math import isqrt

# --------------------------------------------------------------------------------------
# exact arithmetic in a real quadratic field Q(sqrt d), d > 1 squarefree
# --------------------------------------------------------------------------------------


def _sgn_int_surd(P, Q, d):
    """sign of the real number P + Q*sqrt(d) for integers P, Q and d > 0."""
    if Q == 0:
        return (P > 0) - (P < 0)
    if P == 0:
        return (Q > 0) - (Q < 0)
    if P > 0 and Q > 0:
        return 1
    if P < 0 and Q < 0:
        return -1
    c = P * P - Q * Q * d          # sign of P^2 - Q^2 d decides when signs differ
    if P > 0:                      # Q < 0 : positive iff P^2 > Q^2 d
        return 1 if c > 0 else (-1 if c < 0 else 0)
    return -1 if c > 0 else (1 if c < 0 else 0)   # P < 0, Q > 0


class Q2:
    """the real number a + b*sqrt(d) with a, b in Q; all operations exact."""

    __slots__ = ("a", "b", "d")

    def __init__(self, a, b, d):
        self.a = F(a)
        self.b = F(b)
        self.d = d

    # -- coercion / arithmetic ---------------------------------------------------------
    def _c(self, o):
        return o if isinstance(o, Q2) else Q2(F(o), 0, self.d)

    def __add__(self, o):
        o = self._c(o)
        assert o.d == self.d
        return Q2(self.a + o.a, self.b + o.b, self.d)

    def __sub__(self, o):
        o = self._c(o)
        assert o.d == self.d
        return Q2(self.a - o.a, self.b - o.b, self.d)

    def __neg__(self):
        return Q2(-self.a, -self.b, self.d)

    def __mul__(self, o):
        o = self._c(o)
        assert o.d == self.d
        return Q2(self.a * o.a + self.b * o.b * self.d,
                  self.a * o.b + self.b * o.a, self.d)

    def inv(self):
        n = self.a * self.a - self.b * self.b * self.d
        assert n != 0
        return Q2(self.a / n, -self.b / n, self.d)

    def __truediv__(self, o):
        return self * self._c(o).inv()

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, o):
        return self._c(o) - self

    def __rtruediv__(self, o):
        return self._c(o) * self.inv()

    # -- order ------------------------------------------------------------------------
    def sign(self):
        L = self.a.denominator * self.b.denominator
        return _sgn_int_surd(int(self.a * L), int(self.b * L), self.d)

    def cmp(self, o):
        return (self - self._c(o)).sign()

    # -- integer part ------------------------------------------------------------------
    def floor(self):
        """exact floor, by integer arithmetic only."""
        L = self.a.denominator * self.b.denominator
        P, Q = int(self.a * L), int(self.b * L)     # value = (P + Q sqrt d)/L
        # a certified integer lower bound for Q*sqrt(d), within 1
        if Q >= 0:
            S = isqrt(Q * Q * self.d)
        else:
            S = -isqrt(Q * Q * self.d) - 1
        m = (P + S) // L
        while _sgn_int_surd(P - (m + 1) * L, Q, self.d) >= 0:
            m += 1
        while _sgn_int_surd(P - m * L, Q, self.d) < 0:
            m -= 1
        return m

    def frac(self):
        return self - self.floor()

    def __repr__(self):
        return "(%s + %s*sqrt%d)" % (self.a, self.b, self.d)


def rat(x, d):
    return Q2(F(x), 0, d)


# --------------------------------------------------------------------------------------
# the Beatty map, the range condition, and the range sieve
# --------------------------------------------------------------------------------------


def f_alpha(alpha, y):
    """f_alpha(y) = floor(alpha*y)."""
    return (alpha * y).floor()


def inv_powers(alpha, n):
    """{k: alpha^(-k)} for k = 0..n."""
    inv = alpha.inv()
    out = {0: rat(1, alpha.d)}
    for k in range(1, n + 1):
        out[k] = out[k - 1] * inv
    return out


def D_i(alpha, x, i, ip):
    """{x/alpha^i} - {x/alpha^(i-1)}/alpha, exactly (i >= 2)."""
    return (ip[i] * x).frac() - ((ip[i - 1] * x).frac()) / alpha


def rho(alpha, x, k, ip):
    """rho_k(x) = (1 - {x/alpha^k})/alpha, and rho_0 = 0."""
    if k == 0:
        return rat(0, alpha.d)
    return (rat(1, alpha.d) - (ip[k] * x).frac()) / alpha


def RC_printed(alpha, x, n, ip):
    """(holds?, first failing index) for the criterion AS PRINTED in the paper."""
    one = rat(1, alpha.d)
    inv = one / alpha
    hi = one - inv                      # 1 - 1/alpha
    lo = one - inv * 2                  # 1 - 2/alpha
    t = (ip[1] * x).frac()
    if t.cmp(hi) <= 0 or t.cmp(one) >= 0:
        return False, 1
    for i in range(2, n + 1):
        D = D_i(alpha, x, i, ip)
        if D.cmp(lo) <= 0 or D.cmp(hi) >= 0:
            return False, i
    return True, None


def RC_rho(alpha, x, n, ip):
    """(holds?, first failing index) for the criterion in the rho-form."""
    one = rat(1, alpha.d)
    bound = one / (alpha * alpha)
    for k in range(1, n + 1):
        v = rho(alpha, x, k, ip) - rho(alpha, x, k - 1, ip) / alpha
        if v.sign() <= 0 or v.cmp(bound) >= 0:
            return False, k
    return True, None


def range_levels(alpha, n, X):
    """[S_0,...,S_n] with S_k = range(f^k) intersected with [1,X].

    Exact and complete: f_alpha is strictly increasing with f(y) >= y for alpha > 1, so
    f^k(y) = x <= X forces y <= X; the forward sieve from {1..X} therefore misses nothing.
    """
    cur = set(range(1, X + 1))
    out = [cur]
    for _ in range(n):
        nxt = set()
        for y in cur:
            v = f_alpha(alpha, y)
            if v <= X:
                nxt.add(v)
        cur = nxt
        out.append(cur)
    return out


# --------------------------------------------------------------------------------------
# check bookkeeping
# --------------------------------------------------------------------------------------

_PASSED = []
_FAILED = []


def check(name, ok, detail=""):
    if ok:
        _PASSED.append(name)
        print("PASS %s%s" % (name, ("  " + detail) if detail else ""))
    else:
        _FAILED.append(name)
        print("FAIL %s%s" % (name, ("  " + detail) if detail else ""))


# exact check of a printed decimal, truncated toward zero
def trunc_ok(value_as_fraction_bracket, s):
    """s is the decimal string as printed; value is given as a callable cmp(Fraction)->int.

    Returns True iff truncation-toward-zero of the value to len(frac digits) equals s.
    """
    cmpf = value_as_fraction_bracket
    neg = s.startswith("-")
    body = s[1:] if neg else s
    ip, _, fp = body.partition(".")
    k = len(fp)
    N = int(ip + fp)
    if neg:
        N = -N
    scale = 10 ** k
    lo = F(N, scale)
    hi = F(N + 1, scale)
    if not neg:
        return cmpf(lo) >= 0 and cmpf(hi) < 0        # lo <= v < hi
    return cmpf(lo) <= 0 and cmpf(F(N - 1, scale)) > 0   # N-1 < v <= N


def q2_cmp_fraction(v):
    def g(q):
        return (v - Q2(q, 0, v.d)).sign()
    return g


def poly_root_cmp(coeffs, lo_guard):
    """For a polynomial p (given as coeffs, highest degree first) that is strictly
    increasing on [lo_guard, oo) and has its greatest real root there, return a
    comparator g(q) = sign(root - q) valid for q >= lo_guard."""
    def p(t):
        acc = F(0)
        for c in coeffs:
            acc = acc * t + F(c)
        return acc

    def g(q):
        q = F(q)
        assert q >= lo_guard
        s = p(q)
        return 1 if s < 0 else (-1 if s > 0 else 0)   # root > q iff p(q) < 0
    return g


# --------------------------------------------------------------------------------------
# 1. the arithmetic itself
# --------------------------------------------------------------------------------------

def check_arithmetic():
    bad = []
    n = 0
    for d in (2, 3, 5, 7, 17):
        for p in range(-9, 10):
            for q in range(-9, 10):
                for r in (1, 2, 3, 7, 8, 47, 223):
                    v = Q2(F(p, r), F(q, r), d)
                    m = v.floor()
                    n += 1
                    # certify floor by two exact comparisons
                    if (v - Q2(m, 0, d)).sign() < 0:
                        bad.append(("floor-low", p, q, r, d, m))
                    if (v - Q2(m + 1, 0, d)).sign() >= 0:
                        bad.append(("floor-high", p, q, r, d, m))
                    fr = v.frac()
                    if fr.sign() < 0 or (fr - Q2(1, 0, d)).sign() >= 0:
                        bad.append(("frac-range", p, q, r, d))
                    if (v - (Q2(m, 0, d) + fr)).sign() != 0:
                        bad.append(("frac-sum", p, q, r, d))
    # sign() against a squared identity, and inverse/multiplication round trips
    for d in (2, 5, 13):
        for p in range(1, 12):
            for q in range(1, 12):
                v = Q2(p, q, d)
                if (v * v - Q2(p * p + q * q * d, 2 * p * q, d)).sign() != 0:
                    bad.append(("square", p, q, d))
                if (v * v.inv() - Q2(1, 0, d)).sign() != 0:
                    bad.append(("inverse", p, q, d))
                n += 2
    # integers have exact floors and zero fractional part
    for d in (2, 5):
        for k in range(-20, 21):
            if Q2(k, 0, d).floor() != k or Q2(k, 0, d).frac().sign() != 0:
                bad.append(("integer", k, d))
            n += 1
    check("exact-arithmetic-selftest", not bad,
          "%d certified assertions on floor/frac/sign/inverse, %d failures" % (n, len(bad)))


# --------------------------------------------------------------------------------------
# 2. the two moduli, read from the paper as (minimal polynomial, isolating interval)
# --------------------------------------------------------------------------------------

PHI = Q2(F(1, 2), F(1, 2), 5)            # (1+sqrt5)/2
AL2 = Q2(F(7, 8), F(1, 8), 2)            # (7+sqrt2)/8


def check_moduli():
    ok = ((PHI * PHI - PHI - 1).sign() == 0
          and (PHI - F(16, 10)).sign() > 0 and (PHI - F(17, 10)).sign() < 0)
    check("alpha-phi-minimal-polynomial", ok,
          "phi^2-phi-1 = 0 and phi in (1.6,1.7)")
    ok2 = ((AL2 * AL2 * 64 - AL2 * 112 + 47).sign() == 0
           and (AL2 - F(105, 100)).sign() > 0 and (AL2 - F(106, 100)).sign() < 0)
    check("alpha-surd-minimal-polynomial", ok2,
          "64a^2-112a+47 = 0 and a in (1.05,1.06)")


# --------------------------------------------------------------------------------------
# 3. the two forms of the range condition agree
# --------------------------------------------------------------------------------------

def check_rc_forms():
    bad = 0
    n = 0
    for alpha in (PHI, AL2, Q2(0, 1, 3), Q2(1, 1, 2), Q2(F(1, 2), F(1, 2), 17)):
        ip = inv_powers(alpha, 5)
        for x in range(1, 201):
            for k in range(1, 6):
                n += 1
                if RC_printed(alpha, x, k, ip) != RC_rho(alpha, x, k, ip):
                    bad += 1
    check("rc-form-equivalence", bad == 0,
          "printed fractional-part form == rho form on %d instances, %d mismatches"
          % (n, bad))


# --------------------------------------------------------------------------------------
# 4. the thresholds a_n
# --------------------------------------------------------------------------------------

def p_n(n, t):
    t = F(t)
    return t ** (n + 1) - 2 * t ** n + 1


A_N_PRINTED = {2: "1.6180339887", 3: "1.8392867552", 4: "1.9275619754",
               5: "1.9659482366", 6: "1.9835828434"}


def check_nbonacci():
    bad = []
    for n, s in sorted(A_N_PRINTED.items()):
        eps = F(1, 10 ** 10)
        D = F(int(s.replace(".", "")), 10 ** 10)
        # p_n < 0 on (1, a_n) and > 0 above; so D < a_n < D+eps
        if not (p_n(n, D) < 0 < p_n(n, D + eps)):
            bad.append(("bracket", n))
        if not (p_n(n, F(1)) == 0 and p_n(n, F(2)) == 1):
            bad.append(("endpoints", n))
        if not (p_n(n, F(2 * n, n + 1)) < 0):
            bad.append(("branch", n))
        # a_{n+1} > a_n : p_{n+1}(a_n) < p_n(a_n) = 0, checked at D and D+eps
        if n + 1 <= 8 and not (p_n(n + 1, D + eps) < 0):
            bad.append(("increasing", n))
    if (PHI - F(int(A_N_PRINTED[2].replace(".", "")), 10 ** 10)).sign() <= 0:
        bad.append(("phi-is-a2-low",))
    if (PHI - (F(int(A_N_PRINTED[2].replace(".", "")), 10 ** 10) + F(1, 10 ** 10))).sign() >= 0:
        bad.append(("phi-is-a2-high",))
    # p_2(t) = (t-1)(t^2-t-1), so a_2 = phi exactly
    for t in (F(3, 2), F(7, 4), F(5, 3), F(9, 5)):
        if p_n(2, t) != (t - 1) * (t * t - t - 1):
            bad.append(("p2-factorisation", t))
    check("nbonacci-constants", not bad,
          "a_2..a_6 bracketed to 10 decimals, a_2 = phi, monotone; %d failures"
          % len(bad))


def check_threshold_forms():
    bad = 0
    n = 0
    for alpha in (PHI, AL2, Q2(0, 1, 3), Q2(1, 1, 2), Q2(0, 1, 5), Q2(F(3, 2), F(1, 2), 2)):
        one = rat(1, alpha.d)
        for k in range(1, 9):
            n += 1
            # form 1: alpha^k >= alpha^(k-1) + ... + 1
            lhs = one
            for _ in range(k):
                lhs = lhs * alpha
            rhs = rat(0, alpha.d)
            t = one
            for _ in range(k):
                rhs = rhs + t
                t = t * alpha
            f1 = (lhs - rhs).sign() >= 0
            # form 2: alpha^(k+1) - 2 alpha^k + 1 >= 0
            f2 = (lhs * alpha - lhs * 2 + one).sign() >= 0
            # form 3: sum_{j=1}^{k-1} alpha^-j <= alpha - 1
            ip = inv_powers(alpha, k)
            s3 = rat(0, alpha.d)
            for j in range(1, k):
                s3 = s3 + ip[j]
            f3 = (s3 - (alpha - one)).sign() <= 0
            # form 4: sum_{j=2}^{k} alpha^-j <= 1 - 1/alpha
            s4 = rat(0, alpha.d)
            for j in range(2, k + 1):
                s4 = s4 + ip[j]
            f4 = (s4 - (one - one / alpha)).sign() <= 0
            if not (f1 == f2 == f3 == f4):
                bad += 1
    check("threshold-four-equivalent-forms", bad == 0,
          "the four forms of the threshold condition agree on %d (alpha,n) pairs, "
          "%d mismatches" % (n, bad))


def check_threshold_above_2():
    bad = 0
    n = 0
    for alpha in (Q2(1, 1, 2), Q2(0, 1, 5), Q2(3, 1, 2), Q2(0, 1, 7),
                  Q2(F(1, 2), F(1, 2), 17), Q2(F(5, 2), F(1, 2), 3)):
        assert (alpha - 2).sign() > 0
        one = rat(1, alpha.d)
        ip = inv_powers(alpha, 21)
        s = rat(0, alpha.d)
        for j in range(1, 21):
            s = s + ip[j]
            n += 1
            if (s - (alpha - one)).sign() > 0:
                bad += 1
    check("threshold-alpha-gt-2", bad == 0,
          "sum_{j<n} alpha^-j < alpha-1 for all n<=20 at 6 moduli above 2, %d failures"
          % bad)


# --------------------------------------------------------------------------------------
# 5. Theorem A and Proposition 4, by exhaustive equivalence checks
# --------------------------------------------------------------------------------------

MODULI_A = [
    ("3+sqrt2", Q2(3, 1, 2)), ("1+sqrt3", Q2(1, 1, 3)), ("sqrt7", Q2(0, 1, 7)),
    ("(5+sqrt3)/2", Q2(F(5, 2), F(1, 2), 3)), ("1+sqrt2", Q2(1, 1, 2)),
    ("sqrt5", Q2(0, 1, 5)), ("(1+sqrt17)/2", Q2(F(1, 2), F(1, 2), 17)),
    ("sqrt17/2", Q2(0, F(1, 2), 17)), ("(4+sqrt3)/2", Q2(2, F(1, 2), 3)),
]
MODULI_B = [
    ("phi", PHI), ("sqrt3", Q2(0, 1, 3)), ("(3+sqrt2)/2", Q2(F(3, 2), F(1, 2), 2)),
    ("(7+sqrt5)/5", Q2(F(7, 5), F(1, 5), 5)), ("(9+sqrt2)/5", Q2(F(9, 5), F(1, 5), 2)),
    ("(19+sqrt2)/10", Q2(F(19, 10), F(1, 10), 2)),
]

WRONG_DIRECTION = [0]      # count of "RC holds but x not in range" -- must stay 0
SWEEP_INSTANCES = [0]


def sweep(alpha, nmax, X):
    """returns {n: sorted list of x <= X where membership and RC disagree}."""
    lv = range_levels(alpha, nmax, X)
    ip = inv_powers(alpha, nmax)
    out = {}
    for n in range(1, nmax + 1):
        bad = []
        for x in range(1, X + 1):
            inr = x in lv[n]
            rc = RC_printed(alpha, x, n, ip)[0]
            SWEEP_INSTANCES[0] += 1
            if inr != rc:
                bad.append(x)
                if rc and not inr:
                    WRONG_DIRECTION[0] += 1
        out[n] = bad
    return out


def threshold_holds(alpha, n):
    """alpha^n >= alpha^(n-1) + ... + 1, exactly."""
    one = rat(1, alpha.d)
    lhs = one
    for _ in range(n):
        lhs = lhs * alpha
    rhs = rat(0, alpha.d)
    t = one
    for _ in range(n):
        rhs = rhs + t
        t = t * alpha
    return (lhs - rhs).sign() >= 0


SWEEPS = {}


def run_sweeps(X=3000, NMAX=7):
    for lab, alpha in MODULI_A + MODULI_B:
        SWEEPS[lab] = (alpha, sweep(alpha, NMAX, X))


def check_thmA_above_2(X, NMAX):
    tot = 0
    bad = 0
    for lab, _ in MODULI_A:
        alpha, res = SWEEPS[lab]
        assert (alpha - 2).sign() > 0
        for n in range(1, NMAX + 1):
            tot += 1
            bad += len(res[n])
    check("thmA-above-2-exhaustive", bad == 0,
          "%d moduli above 2, n=1..%d, x<=%d: %d disagreements in %d (alpha,n) cells"
          % (len(MODULI_A), NMAX, X, bad, tot))


def check_thmA_per_level(X, NMAX):
    forbidden = 0
    viol = 0
    permitted = 0
    for lab, _ in MODULI_A + MODULI_B:
        alpha, res = SWEEPS[lab]
        for n in range(1, NMAX + 1):
            if threshold_holds(alpha, n):
                forbidden += 1
                if res[n]:
                    viol += 1
            else:
                permitted += 1
    check("thmA-per-level-exhaustive", viol == 0,
          "%d (alpha,n) cells with alpha >= a_n: %d with a disagreement; "
          "%d cells left permitted by the theorem" % (forbidden, viol, permitted))


def check_prop_unconditional():
    check("thmA-unconditional-direction", WRONG_DIRECTION[0] == 0,
          "over %d (alpha,n,x) instances, 'RC holds but x not in range' occurred "
          "%d times" % (SWEEP_INSTANCES[0], WRONG_DIRECTION[0]))


# --------------------------------------------------------------------------------------
# 6. Theorem B over a family, and its corollary
# --------------------------------------------------------------------------------------

FAMILY_D = [2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33,
            34, 35, 37, 38, 39, 41, 42, 43]


def check_family():
    seen = set()
    tested = witnessed = wit = cor_tot = cor_ok = 0
    nonempty_iff_phi_bad = 0
    bad = []
    for d in FAMILY_D:
        for q in range(1, 4):
            for r in range(1, 26):
                for p in range(0, 60):
                    alpha = Q2(F(p, r), F(q, r), d)
                    key = (alpha.a, alpha.b, d)
                    if key in seen:
                        continue
                    if (alpha - F(2001, 2000)).sign() <= 0:
                        continue
                    # alpha < phi  <=>  alpha^2 - alpha - 1 < 0   (for alpha > 1)
                    if (alpha * alpha - alpha - 1).sign() >= 0:
                        continue
                    seen.add(key)
                    tested += 1
                    one = rat(1, d)
                    lowB = (alpha * alpha) / (alpha * alpha - 1)
                    upB = one / (alpha - 1)
                    # the interval is nonempty iff alpha < phi, which holds here
                    if (upB - lowB).sign() <= 0:
                        nonempty_iff_phi_bad += 1
                        continue
                    ints = [x for x in range(lowB.floor() + 1, upB.floor() + 1)
                            if (rat(x, d) - lowB).sign() > 0
                            and (rat(x, d) - upB).sign() < 0]
                    if ints:
                        witnessed += 1
                        ip = inv_powers(alpha, 3)
                        for x in ints:
                            wit += 1
                            if f_alpha(alpha, x) != x:
                                bad.append(("not-fixed", p, q, r, d, x))
                            ok, first = RC_printed(alpha, x, 2, ip)
                            if ok or first != 2:
                                bad.append(("rc2-not-failing-at-2", p, q, r, d, x))
                            if (rat(x, d) / alpha).frac().sign() <= 0:
                                bad.append(("frac-not-positive", p, q, r, d, x))
                            # D_2 = (1 - 1/alpha) + {x/alpha}, exactly
                            lhs = D_i(alpha, x, 2, ip)
                            rhs = (one - one / alpha) + (rat(x, d) / alpha).frac()
                            if (lhs - rhs).sign() != 0:
                                bad.append(("D2-identity", p, q, r, d, x))
                    # corollary: alpha <= (1+sqrt17)/4  <=>  2 alpha^2 - alpha - 2 <= 0
                    if (alpha * alpha * 2 - alpha - 2).sign() <= 0:
                        cor_tot += 1
                        xs = upB.floor()
                        if xs in ints:
                            cor_ok += 1
                        else:
                            bad.append(("xstar-outside", p, q, r, d, xs))
    check("thmB-window-nonempty-iff-phi", nonempty_iff_phi_bad == 0,
          "all %d enumerated alpha in (1,phi) have a nonempty window (Theorem 7), "
          "%d exceptions" % (tested, nonempty_iff_phi_bad))
    check("thmB-family-witnesses", not bad,
          "%d distinct quadratic irrationals of (1,phi); %d with an integer in the "
          "window; "
          "%d integer witnesses, each a fixed point with RC failing exactly at i=2 and "
          "D_2 = (1-1/alpha) + {x/alpha}; %d failures"
          % (tested, witnessed, wit, len(bad)))
    check("thmB-corollary-xstar", cor_ok == cor_tot and cor_tot > 0,
          "x* = floor(1/(alpha-1)) lies in the window for %d of the %d enumerated alpha "
          "with 2a^2-a-2 <= 0" % (cor_ok, cor_tot))
    return tested, witnessed, wit, cor_tot


# --------------------------------------------------------------------------------------
# 7. the worked witness alpha = (7+sqrt2)/8
# --------------------------------------------------------------------------------------

WINDOW_X = list(range(11, 20))


def check_witness2():
    a = AL2
    one = rat(1, 2)
    lowB = (a * a) / (a * a - 1)
    upB = one / (a - 1)
    ok = ((lowB - Q2(F(1055, 223), F(896, 223), 2)).sign() == 0
          and (upB - Q2(8, 8, 2)).sign() == 0
          and (a - 1 - Q2(F(-1, 8), F(1, 8), 2)).sign() == 0
          and (a * a - 1 - Q2(F(-13, 64), F(14, 64), 2)).sign() == 0)
    check("witness2-window-endpoints", ok,
          "alpha^2/(alpha^2-1) = (1055+896 sqrt2)/223 and 1/(alpha-1) = 8+8 sqrt2, "
          "with alpha-1 = (sqrt2-1)/8 and alpha^2-1 = (14 sqrt2-13)/64")

    ints = [x for x in range(1, 40)
            if (rat(x, 2) - lowB).sign() > 0 and (rat(x, 2) - upB).sign() < 0]
    check("witness2-window-integers", ints == WINDOW_X and upB.floor() == 19
          and lowB.floor() + 1 == 11,
          "integers strictly inside the window: %s (nine of them), x* = floor(1/(a-1)) = %d"
          % (",".join(map(str, ints)), upB.floor()))

    fixed = [x for x in WINDOW_X if f_alpha(a, x) == x]
    also = f_alpha(a, 10) == 10 and f_alpha(a, 20) == 21 and f_alpha(a, 21) == 22
    check("witness2-fixed-points", len(fixed) == 9 and also,
          "f(x) = x for all of x = 11..19; f(10) = 10, f(20) = 21, f(21) = 22")

    ip = inv_powers(a, 3)
    fl = all((rat(x, 2) / a).floor() == x - 1
             and (rat(x, 2) / (a * a)).floor() == x - 2 for x in WINDOW_X)
    check("witness2-floors", fl,
          "floor(x/alpha) = x-1 and floor(x/alpha^2) = x-2 for x = 11..19")

    cf = True
    for x in WINDOW_X:
        lhs = D_i(a, x, 2, ip)
        rhs = Q2(F(9 * x + 38, 47), F(-8 * (x - 1), 47), 2)
        rhs2 = (rat(x - 1, 2) / a) - (x - 2)
        rhs3 = (rat(x - 1, 2) / a).frac()
        if (lhs - rhs).sign() != 0 or (lhs - rhs2).sign() != 0 or (lhs - rhs3).sign() != 0:
            cf = False
    check("witness2-D2-closed-form", cf,
          "D_2(alpha,x) = (9x+38-8(x-1)sqrt2)/47 = (x-1)/alpha-(x-2) = {(x-1)/alpha} "
          "for x = 11..19")

    bnd = one - one / a
    exc = all(D_i(a, x, 2, ip).cmp(bnd) > 0 for x in WINDOW_X)
    b_exact = (bnd - Q2(F(-9, 47), F(8, 47), 2)).sign() == 0
    # on the window, 47*(D_2(x) - (1-1/alpha)) = 9x+47-8x sqrt2 ...
    id_ok = all(((D_i(a, x, 2, ip) - bnd) * 47 - Q2(9 * x + 47, -8 * x, 2)).sign() == 0
                for x in WINDOW_X)
    # ... and that expression is positive exactly when x < 9+8 sqrt2, for every x
    eqv = all(((Q2(9 * x + 47, -8 * x, 2).sign() > 0)
               == ((rat(x, 2) - Q2(9, 8, 2)).sign() < 0)) for x in range(1, 61))
    # the excess over the bound is exactly {x/alpha}
    exc_id = all(((D_i(a, x, 2, ip) - bnd) - (rat(x, 2) / a).frac()).sign() == 0
                 for x in WINDOW_X)
    check("witness2-D2-exceeds-bound", exc and b_exact and id_ok and eqv and exc_id,
          "1-1/alpha = (8 sqrt2-9)/47; D_2(alpha,x)-(1-1/alpha) = {x/alpha} > 0 for "
          "x = 11..19; and 47*(D_2(x)-(1-1/alpha)) = 9x+47-8x sqrt2, which is positive "
          "exactly when x < 9+8 sqrt2 (checked for x = 1..60)")

    lv = range_levels(a, 2, 24)
    dis = [x for x in range(1, 25) if (x in lv[2]) != RC_printed(a, x, 2, ip)[0]]
    check("witness2-exhaustive-table-24", dis == WINDOW_X,
          "x = 1..24, n = 2: the equivalence is false at exactly %s (nine values) and "
          "true at the other fifteen" % ",".join(map(str, dis)))


# --------------------------------------------------------------------------------------
# 8. the golden-ratio witness
# --------------------------------------------------------------------------------------

CHAIN = [3, 4, 6, 9, 14, 22]
FRACS = {1: (-35, 22), 2: (36, -22), 3: (-71, 44), 4: (107, -66), 5: (-177, 110)}
DVALS = {2: (-21, 13), 3: (-13, 8), 4: (-8, 5), 5: (-4, 3)}


def zphi(c, k):
    """c + k*phi as an element of Q(sqrt5)."""
    return Q2(F(c) + F(k, 2), F(k, 2), 5)


def check_witness1():
    got = [3]
    for _ in range(5):
        got.append(f_alpha(PHI, got[-1]))
    check("witness1-chain", got == CHAIN,
          "3 -> 4 -> 6 -> 9 -> 14 -> 22 under y -> floor(phi*y), so 22 in range(f^5)")

    lv = range_levels(PHI, 5, 30)
    pre = [y for y in range(1, 31) if _iter(PHI, y, 5) == 22]
    check("witness1-membership-independent", 22 in lv[5] and pre == [3],
          "an exhaustive preimage sieve over y <= 30 independently puts 22 in range(f^5); "
          "its only 5-fold preimage is y = 3")

    ip = inv_powers(PHI, 5)
    fr_ok = all(((ip[i] * 22).frac() - zphi(*FRACS[i])).sign() == 0 for i in FRACS)
    check("witness1-fractional-parts", fr_ok,
          "{22/phi^i} = 22phi-35, 36-22phi, 44phi-71, 107-66phi, 110phi-177 for i = 1..5")

    d_ok = all((D_i(PHI, 22, i, ip) - zphi(*DVALS[i])).sign() == 0 for i in DVALS)
    check("witness1-D-closed-forms", d_ok,
          "D_2 = 13phi-21, D_3 = 8phi-13, D_4 = 5phi-8, D_5 = 3phi-4")

    one = rat(1, 5)
    hi = one - one / PHI
    lo = one - (one / PHI) * 2
    inside = all(D_i(PHI, 22, i, ip).cmp(lo) > 0 and D_i(PHI, 22, i, ip).cmp(hi) < 0
                 for i in (2, 3, 4))
    fails5 = D_i(PHI, 22, 5, ip).cmp(hi) > 0
    excess = (D_i(PHI, 22, 5, ip) - hi - zphi(-6, 4)).sign() == 0
    r4 = RC_printed(PHI, 22, 4, ip)
    r5 = RC_printed(PHI, 22, 5, ip)
    bounds = ((hi - zphi(2, -1)).sign() == 0 and (lo - zphi(3, -2)).sign() == 0)
    check("witness1-rc-verdict",
          inside and fails5 and excess and r4 == (True, None) and r5 == (False, 5)
          and bounds,
          "1-1/phi = 2-phi and 1-2/phi = 3-2phi; RC<=4(phi,22) holds, RC<=5 fails at "
          "i=5, and D_5-(1-1/phi) = 4phi-6 > 0")

    # phi = a_2, so levels 1 and 2 are forbidden to fail; 3,4,5 are permitted
    perm = [n for n in range(1, 8) if not threshold_holds(PHI, n)]
    check("witness1-level-permission", perm == [3, 4, 5, 6, 7] and threshold_holds(PHI, 2),
          "phi = a_2: Theorem A covers n <= 2 at phi and permits failure for n >= 3")


def _iter(alpha, y, k):
    for _ in range(k):
        y = f_alpha(alpha, y)
    return y


def check_phi_least(X4=3000):
    ip = inv_powers(PHI, 5)
    lv5 = range_levels(PHI, 5, X4)
    b5 = [x for x in range(1, X4 + 1) if (x in lv5[5]) != RC_printed(PHI, x, 5, ip)[0]]
    check("witness1-least-x-at-level-5", b5 and b5[0] == 22,
          "least disagreement at (phi, n=5) is x = 22, exhaustive for x <= %d "
          "(%d disagreements there)" % (X4, len(b5)))
    lv4 = range_levels(PHI, 4, X4)
    b4 = [x for x in range(1, X4 + 1) if (x in lv4[4]) != RC_printed(PHI, x, 4, ip)[0]]
    check("phi-level-4-least-137", b4 and b4[0] == 137,
          "least disagreement at (phi, n=4) is x = 137, exhaustive for x <= %d "
          "(%d disagreements there)" % (X4, len(b4)))


def check_phi_level3(X=20000):
    ip = inv_powers(PHI, 3)
    lv = range_levels(PHI, 3, X)
    b = [x for x in range(1, X + 1) if (x in lv[3]) != RC_printed(PHI, x, 3, ip)[0]]
    check("phi-level-3-none-to-20000", b == [],
          "no disagreement at (phi, n=3) for x <= %d; Theorem A permits one, so this "
          "cell is UNDECIDED, not settled" % X)


def check_sqrt3():
    alpha, res = SWEEPS["sqrt3"]
    tot = sum(len(v) for v in res.values())
    perm = [n for n in range(1, 8) if not threshold_holds(alpha, n)]
    check("band-phi-to-2-undecided", tot == 0 and perm == [3, 4, 5, 6, 7],
          "at alpha = sqrt3 in [phi,2) the equivalence holds for all n <= 7, x <= 3000, "
          "although Theorem A covers only n <= 2: bounded-negative, not a proof")


# --------------------------------------------------------------------------------------
# 9. every decimal expansion printed in the paper
# --------------------------------------------------------------------------------------

def check_decimals():
    items = [
        ("phi", q2_cmp_fraction(PHI), "1.6180339887"),
        ("alpha=(7+sqrt2)/8", q2_cmp_fraction(AL2), "1.0517766952966"),
        ("(1+sqrt17)/4", q2_cmp_fraction(Q2(F(1, 4), F(1, 4), 17)), "1.2807764064"),
        ("sqrt3", q2_cmp_fraction(Q2(0, 1, 3)), "1.7320508"),
        ("(1055+896sqrt2)/223",
         q2_cmp_fraction(Q2(F(1055, 223), F(896, 223), 2)), "10.4131630129"),
        ("8+8sqrt2", q2_cmp_fraction(Q2(8, 8, 2)), "19.3137084989"),
        ("9+8sqrt2", q2_cmp_fraction(Q2(9, 8, 2)), "20.3137084989"),
        ("(8sqrt2-9)/47", q2_cmp_fraction(Q2(F(-9, 47), F(8, 47), 2)), "0.0492278404"),
        ("(209-144sqrt2)/47",
         q2_cmp_fraction(Q2(F(209, 47), F(-144, 47), 2)), "0.1138988727"),
        ("(137-80sqrt2)/47",
         q2_cmp_fraction(Q2(F(137, 47), F(-80, 47), 2)), "0.5077215959"),
        ("13phi-21", q2_cmp_fraction(zphi(-21, 13)), "0.0344418537"),
        ("8phi-13", q2_cmp_fraction(zphi(-13, 8)), "-0.0557280900"),
        ("5phi-8", q2_cmp_fraction(zphi(-8, 5)), "0.0901699437"),
        ("3phi-4", q2_cmp_fraction(zphi(-4, 3)), "0.8541019662"),
        ("2-phi", q2_cmp_fraction(zphi(2, -1)), "0.3819660112"),
        ("3-2phi", q2_cmp_fraction(zphi(3, -2)), "-0.2360679774"),
        ("4phi-6", q2_cmp_fraction(zphi(-6, 4)), "0.4721359549"),
    ]
    # a_n are algebraic of degree n; bracket them by the sign of p_n
    for n, s in sorted(A_N_PRINTED.items()):
        items.append(("a_%d" % n, poly_root_cmp(
            [1, -2] + [0] * (n - 1) + [1], F(1)), s))
    # the greatest root of t^4-2t^3-t+1: p' = 4t^3-6t^2-1 > 0 on [8/5,oo)
    items.append(("root of t^4-2t^3-t+1",
                  poly_root_cmp([1, -2, 0, -1, 1], F(8, 5)), "2.1176886328"))
    bad = [nm for nm, cmpf, s in items if not trunc_ok(cmpf, s)]
    # certify the monotonicity guard used for the quartic
    guard = (4 * F(8, 5) ** 3 - 6 * F(8, 5) ** 2 - 1 > 0)
    check("printed-decimals", not bad and guard,
          "%d decimal expansions printed in the paper re-derived by exact rational "
          "bracketing (truncation toward zero); %d wrong" % (len(items), len(bad)))


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------

def main():
    X, NMAX = 3000, 7
    check_arithmetic()
    check_moduli()
    check_rc_forms()
    check_nbonacci()
    check_threshold_forms()
    check_threshold_above_2()
    run_sweeps(X, NMAX)
    check_thmA_above_2(X, NMAX)
    check_thmA_per_level(X, NMAX)
    check_prop_unconditional()
    check_sqrt3()
    check_family()
    check_witness2()
    check_witness1()
    check_phi_least()
    check_phi_level3()
    check_decimals()

    print()
    print("NOT RE-RUN: the following are recorded in the paper but are NOT re-derived here.")
    print("NOT RE-RUN:  (a) the machine census reporting 38,417 disagreements at "
          "alpha=(7+sqrt2)/8 for x <= 10^5. No claim of the paper depends on it; this "
          "program checks x <= 24 exhaustively and nothing beyond.")
    print("NOT RE-RUN:  (b) any statement about which inequalities the proof printed in "
          "arXiv:2607.12817 uses (the threshold remark of Section 3). That is a reading of "
          "a text, "
          "not a computation; the line numbers, byte count and sha256 given in the paper "
          "let a referee repeat it, and this program does not fetch the source.")
    print("NOT RE-RUN:  (c) the literature search behind the novelty remark of "
          "Section 5 and behind Section 6. No "
          "database is contacted.")
    print("NOT RE-RUN:  (d) the open cells are not resolved: whether the equivalence "
          "fails anywhere in the OPEN band (phi,2), and whether it fails at alpha=phi "
          "for n=3, are left open. The searches here are BOUNDED (x <= 20000 at n=3, "
          "x <= 3000 otherwise) and a bounded negative is not a proof. The left "
          "endpoint alpha=phi is NOT among these open cells: the failure at level 5 "
          "there is established above (least disagreement x = 22, exhaustive for "
          "x <= 3000), so the band written here is open at phi.")
    print("NOT RE-RUN:  (e) the family sweep behind Theorem 7 (the negative half) "
          "enumerates quadratic irrationals only, of the shape (p+q sqrt d)/r with "
          "d <= 43, q <= 3, r <= 25, p < 60; the theorem itself is proved for all "
          "irrational alpha in (1,phi) and its proof, not this sweep, is what carries it.")
    print()
    n = len(_PASSED) + len(_FAILED)
    if _FAILED:
        print("VERDICT: %d of %d CHECKS FAILED: %s"
              % (len(_FAILED), n, ", ".join(_FAILED)))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Verification of a counterexample to the aggregate-condition form of
Karp--Sitnik Conjecture 1.

TAKEN FROM THE PAPER (inputs, not re-derived):
  - the statement being refuted (aggregate condition + positivity => bound),
    together with the READING of Karp--Sitnik Conjecture 1 on which the
    aggregate condition REPLACES the componentwise inequalities b_i > a_i.
    That reading is an exegetical premise about someone else's text: it is
    printed at the head of the transcript, flagged again in the closing
    NOT RE-RUN paragraphs, and is not machine-checkable here;
  - every quotation from and every theorem/conjecture number of [KS], [KP]
    and [Derbazi]: no external source is read by this program;
  - the exhibited parameters sigma = 2, A = (2,2), B = (1,5), x = 1/2,
    i.e. 3F2(2,2,2; 1,5; -1/2);
  - the paper's closed form 1032 - 2544*log(3/2);
  - the paper's comparison value 25/49 = (1 + mu*x)^{-sigma} with mu = 4/5.

DERIVED HERE (computed by this program, exact integer/rational arithmetic
plus certified rational enclosures of log(3/2)):
  - the hypergeometric series value, summed to a rigorously bounded tail;
  - agreement of that value with the paper's closed form;
  - the aggregate hypothesis sum(b_i - a_i) > 0 and positivity, and the
    failure of the componentwise hypothesis b_i > a_i;
  - mu = prod a_i/b_i < 1 and the numeric value of (1 + mu*x)^{-sigma};
  - the strict violation H < (1 + mu*x)^{-sigma} (load-bearing);
  - the Karp--Sitnik bound holding on componentwise-admissible parameters;
  - the quadratic-expansion family claim and a search over small integer
    parameters.
"""

import math
from fractions import Fraction as F

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ((" [" + detail + "]") if detail else ""))


def dec(fr, k=30):
    """Exact decimal expansion of a Fraction, truncated toward zero to k
    places.  Never converts to float: %f on a Fraction would silently go
    through binary floating point and fabricate the digits past the 17th."""
    fr = F(fr)
    sign = "-" if fr < 0 else ""
    n = abs(fr)
    scaled = (n.numerator * 10 ** k) // n.denominator
    s = str(scaled).rjust(k + 1, "0")
    return sign + s[:-k] + "." + s[-k:]


def pmul(p, q):
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, u in enumerate(p):
        for j, v in enumerate(q):
            r[i + j] += u * v
    return r


def pshift(p, s):
    """Coefficients (low to high) of p(t+s) as a polynomial in t."""
    res = [F(0)] * len(p)
    pw = [F(1)]
    for i, c in enumerate(p):
        if i > 0:
            pw = pmul(pw, [F(s), F(1)])
        for j, v in enumerate(pw):
            res[j] += c * v
    return res


def ratio_polys(sigma, A, B, x):
    """|t_{n+1}/t_n| = x * P(n) / Q(n) for H_{sigma,A,B}(x), x > 0."""
    P = [F(sigma), F(1)]
    for a in A:
        P = pmul(P, [F(a), F(1)])
    Q = [F(1), F(1)]
    for b in B:
        Q = pmul(Q, [F(b), F(1)])
    return P, Q


def certify_ratio_bound(sigma, A, B, x, rho, n0):
    """True if x*P(n)/Q(n) <= rho is PROVED for every integer n >= n0.

    Certificate: R = rho*Q - x*P has all coefficients >= 0 after the
    substitution n = n0 + t, hence R(n) >= 0 for all n >= n0 (Q > 0)."""
    x = abs(F(x))
    if not (F(0) < F(rho) < 1):
        return False
    P, Q = ratio_polys(sigma, A, B, x)
    deg = max(len(P), len(Q))
    R = [F(0)] * deg
    for i, c in enumerate(Q):
        R[i] += F(rho) * c
    for i, c in enumerate(P):
        R[i] -= F(x) * c
    return all(c >= 0 for c in pshift(R, n0))


def hyper_terms(sigma, A, B, x, N):
    """Terms t_0..t_N of H_{sigma,A,B}(x) = sum_n coeff(n) * (-x)^n."""
    t = F(1)
    out = [t]
    for n in range(N):
        num = F(sigma + n)
        for a in A:
            num *= F(a + n)
        den = F(1 + n)
        for b in B:
            den *= F(b + n)
        t = t * num / den * (-F(x))
        out.append(t)
    return out


def hyper_enclosure(sigma, A, B, x, N):
    """Exact rational enclosure (lo, hi) of H_{sigma,A,B}(x), |x| < 1.

    The tail past t_N is bounded by |t_N| * rho/(1-rho) where rho is a
    ratio bound PROVED for all n >= N by certify_ratio_bound."""
    ts = hyper_terms(sigma, A, B, x, N)
    S = sum(ts)
    rho = None
    for cand in (F(11, 10), F(5, 4), F(3, 2), F(7, 4), F(19, 10), F(2),
                 F(3), F(5), F(8), F(12)):
        r = abs(F(x)) * cand
        if r < 1 and certify_ratio_bound(sigma, A, B, x, r, N):
            rho = r
            break
    if rho is None:
        raise ValueError("no certified ratio bound")
    bound = abs(ts[N]) * rho / (1 - rho)
    return S - bound, S + bound


def log_enclosure(K=70):
    """Enclosure of log(3/2) = 2*artanh(1/5) by its exact rational series."""
    lo = F(0)
    for k in range(K):
        lo += F(2, (2 * k + 1) * 5 ** (2 * k + 1))
    tail = F(2, (2 * K + 1) * 5 ** (2 * K + 1)) * F(25, 24)
    return lo, lo + tail


def exp_enclosure(v, N=60):
    """Enclosure of exp(v) for 0 <= v < 1 from the exact Taylor series."""
    s = F(0)
    term = F(1)
    for n in range(N):
        if n:
            term = term * v / n
        s += term
    tail = term * v / N * 2
    return s, s + tail


# ---------------------------------------------------------------- paper data
SIGMA = 2                 # the sigma of H_{sigma,A,B}
A = (2, 2)                # numerator parameters a_1, a_2
B = (1, 5)                # denominator parameters b_1, b_2
X = F(1, 2)               # evaluation point x, so the argument is -x = -1/2
CLAIM_MU = F(4, 5)        # paper: mu(A,B) = 4/5
CLAIM_RHS = F(25, 49)     # paper: (1 + mu x)^{-sigma} = 25/49
CLAIM_CONST = 1032        # paper closed form 1032 - 2544*log(3/2)
CLAIM_LOGCO = 2544
CLAIM_LOG_LB = F(19006, 46875)   # paper: log(3/2) > 19006/46875
CLAIM_MIDDLE = F(7912, 15625)    # paper: H < 7912/15625 < 25/49
CLAIM_DISCR = 2937               # paper: 25*15625 - 49*7912 = 2937
N_TERMS = 220

# Transcribed by hand from the paper, which in turn quotes the source; printed
# at the head of the transcript so a referee can compare it against the source
# by eye.  Nothing below reads any external document.
KS_CONJECTURE_TEXT = ("Theorem 3 is true for all sigma>0 and "
                      "sum_{i=1}^{q}(b_i-a_i)>0")
KS_CONJECTURE_SOURCE = "Conjecture 1 of Karp-Sitnik, in the arXiv numbering"
READING = ("sum(b_i-a_i)>0 REPLACES the componentwise inequalities b_i>a_i, "
           "all parameters remaining positive")


def check_object():
    top = (SIGMA,) + A
    bot = B
    ok = (len(top) == 3 and len(bot) == 2
          and top == (2, 2, 2) and bot == (1, 5)
          and all(isinstance(v, int) and v > 0 for v in top + bot)
          and X > 0 and X > -1)
    ck("exhibited-object-is-3F2-with-stated-integer-parameters", ok,
       "top=%s bot=%s argument=-%s q=%d sigma=%d" % (top, bot, X, len(A), SIGMA))


def check_hypotheses():
    agg = sum(B) - sum(A)
    positive = all(v > 0 for v in A + B + (SIGMA,))
    componentwise = all(b > a for a, b in zip(A, B))
    mu = F(1)
    for a, b in zip(A, B):
        mu *= F(a, b)
    ck("aggregate-hypothesis-holds-and-parameters-positive",
       agg > 0 and positive and SIGMA >= 1,
       "sum(b_i-a_i)=%d>0, all parameters>0, sigma=%d>=1" % (agg, SIGMA))
    ck("componentwise-hypothesis-b_i>a_i-fails-so-Theorem-3-does-not-apply",
       not componentwise,
       "pairs (a_i,b_i) = %s; b_i>a_i holds for %d of %d"
       % (list(zip(A, B)), sum(1 for a, b in zip(A, B) if b > a), len(A)))
    ck("mu-equals-4/5-and-is-less-than-1", mu == CLAIM_MU and mu < 1,
       "mu=%s" % mu)
    rhs = (1 + mu * X) ** (-SIGMA)
    ck("right-hand-side-(1+mu*x)^-sigma-equals-25/49",
       rhs == CLAIM_RHS, "(1+%s*%s)^-%d = %s" % (mu, X, SIGMA, rhs))


def poch(a, n):
    p = F(1)
    for j in range(n):
        p *= F(a + j)
    return p


def fact(n):
    p = 1
    for j in range(2, n + 1):
        p *= j
    return p


def check_coefficients():
    """Series coefficient of z^n, from the definition, against both closed
    forms printed in the paper."""
    bad_rat, bad_pf, first = [], [], []
    for n in range(0, 301):
        defn = (poch(SIGMA, n) * poch(A[0], n) * poch(A[1], n)
                / (poch(B[0], n) * poch(B[1], n) * fact(n)))
        rat = F(24 * (n + 1) ** 2, (n + 2) * (n + 3) * (n + 4))
        pf = F(12, n + 2) - F(96, n + 3) + F(108, n + 4)
        if defn != rat:
            bad_rat.append(n)
        if defn != pf:
            bad_pf.append(n)
        if n < 3:
            first.append(str(defn))
    ck("coefficient-equals-24(n+1)^2/((n+2)(n+3)(n+4))", not bad_rat,
       "n=0..300 exact; c_0,c_1,c_2 = " + ", ".join(first))
    ck("coefficient-partial-fractions-12/(n+2)-96/(n+3)+108/(n+4)", not bad_pf,
       "n=0..300 exact")


def check_series_coeffs_match_terms():
    """The generic term generator used below reproduces the same coefficients."""
    ts = hyper_terms(SIGMA, A, B, X, 40)
    bad = []
    for n in range(41):
        c = F(24 * (n + 1) ** 2, (n + 2) * (n + 3) * (n + 4))
        if ts[n] != c * (-X) ** n:
            bad.append(n)
    ck("generic-term-generator-agrees-with-explicit-coefficients", not bad,
       "n=0..40 exact")


def check_violation():
    """LOAD-BEARING: sum the series to a certified enclosure and compare."""
    lo, hi = hyper_enclosure(SIGMA, A, B, X, N_TERMS)
    mu = F(A[0], B[0]) * F(A[1], B[1])
    rhs = (1 + mu * X) ** (-SIGMA)
    ck("certified-ratio-bound-for-the-series-tail",
       certify_ratio_bound(SIGMA, A, B, X, F(3, 5), N_TERMS)
       and hi - lo < F(1, 10 ** 40),
       "enclosure width < 10^-40")
    # lo > 0 is part of the claim, not decoration: a violation of
    # (1 + mu x)^{-sigma} < H must be witnessed by a genuine positive value of
    # the series, so a corrupted object whose "H" came out negative must not
    # be accepted here merely for lying below the right-hand side.
    ck("H-is-strictly-below-(1+mu*x)^-sigma--CONCLUSION-FAILS",
       F(0) < lo and hi < rhs,
       "H in [%s, %s], rhs = %s, gap >= %s"
       % (dec(lo, 45), dec(hi, 45), dec(rhs, 45), dec(rhs - hi, 45)))
    ck("H-lies-in-(0,1)-so-the-series-value-is-sane", F(0) < lo and hi < 1,
       "H ~ %s" % dec(lo, 20))
    ck("H-is-below-the-paper-intermediate-bound-7912/15625",
       hi < CLAIM_MIDDLE, "H <= %s < %s = 7912/15625"
       % (dec(hi, 20), dec(CLAIM_MIDDLE, 20)))


def check_closed_form():
    """1032 - 2544*log(3/2) agrees with the summed series."""
    llo, lhi = log_enclosure()
    clo = CLAIM_CONST - CLAIM_LOGCO * lhi
    chi = CLAIM_CONST - CLAIM_LOGCO * llo
    slo, shi = hyper_enclosure(SIGMA, A, B, X, N_TERMS)
    ck("closed-form-1032-2544log(3/2)-matches-the-summed-series",
       clo <= shi and slo <= chi and (chi - clo) < F(1, 10 ** 40),
       "closed form in [%s, %s]" % (dec(clo, 45), dec(chi, 45)))
    D = 10 ** 40
    down = F(math.floor(llo * D), D)
    up = F(math.ceil(lhi * D), D)
    ck("log(3/2)-enclosure-is-consistent-with-exp-(independent-check)",
       exp_enclosure(down)[1] <= F(3, 2) <= exp_enclosure(up)[0]
       and (F(1, 5) + 1) / (1 - F(1, 5)) == F(3, 2),
       "exp(lo) <= 3/2 <= exp(hi); (1+1/5)/(1-1/5)=3/2")
    ck("paper-lower-bound-log(3/2)>19006/46875-is-correct",
       llo > CLAIM_LOG_LB
       and sum(F(2, (2 * k + 1) * 5 ** (2 * k + 1)) for k in range(3))
       == CLAIM_LOG_LB,
       "log(3/2) > %s = 19006/46875" % dec(CLAIM_LOG_LB, 20))


def geom_log_pair(r, z=F(-1, 2)):
    """sum_{n>=0} z^n/(n+r) written as const + coefL*log(3/2) at z=-1/2,
    using z^{-r}(-log(1-z) - sum_{j<r} z^j/j) and -log(1-z) = -log(3/2)."""
    zr = z ** (-r)
    partial = sum(z ** j / F(j) for j in range(1, r))
    return -zr * partial, -zr


def check_derivation():
    """The partial-fraction sum reproduces the coefficients 1032 and -2544."""
    const = F(0)
    coefL = F(0)
    for weight, r in ((12, 2), (-96, 3), (108, 4)):
        c, l = geom_log_pair(r)
        const += weight * c
        coefL += weight * l
    ck("partial-fraction-sums-combine-to-1032-2544log(3/2)",
       const == CLAIM_CONST and coefL == -CLAIM_LOGCO,
       "const=%s coefficient of log(3/2)=%s" % (const, coefL))
    # each sum_{n} z^n/(n+r) closed form against its own series
    llo, lhi = log_enclosure()
    bad = []
    for r in (2, 3, 4):
        s = sum(F(-1, 2) ** n / F(n + r) for n in range(201))
        tail = F(1, 2 ** 200 * (201 + r)) * 2
        c, l = geom_log_pair(r)
        flo = c + (l * lhi if l > 0 else l * llo)
        fhi = c + (l * llo if l > 0 else l * lhi)
        if not (flo <= s + tail and s - tail <= fhi):
            bad.append(r)
    ck("closed-form-of-sum-z^n/(n+r)-verified-against-its-series-r=2,3,4",
       not bad, "z=-1/2, 201 exact terms with tail bound")


def check_rational_chain():
    """The paper's rational chain from the log lower bound down to 25/49."""
    mid = CLAIM_CONST - CLAIM_LOGCO * CLAIM_LOG_LB
    ck("1032-2544*(19006/46875)-equals-7912/15625", mid == CLAIM_MIDDLE,
       "value = %s" % mid)
    d = 25 * CLAIM_MIDDLE.denominator - 49 * CLAIM_MIDDLE.numerator
    ck("cross-multiplied-discriminant-2937-shows-7912/15625<25/49",
       d == CLAIM_DISCR and d > 0 and CLAIM_MIDDLE < CLAIM_RHS,
       "25*15625 - 49*7912 = %d" % d)


def prods(A, B):
    mu = F(1)
    PA = F(1)
    PB = F(1)
    for a, b in zip(A, B):
        mu *= F(a) / F(b)
        PA *= 1 + 1 / F(a)
        PB *= 1 + 1 / F(b)
    return mu, PA, PB


def check_taylor_coefficients():
    """Proposition: H = 1 - sigma*mu*x + sigma(sigma+1)/2*mu^2*(P_A/P_B)x^2+...
    while (1+mu x)^{-sigma} = 1 - sigma*mu*x + sigma(sigma+1)/2*mu^2 x^2+..."""
    bad1, bad2, bad3 = [], [], []
    cases = []
    for sg in (F(1), F(2), F(1, 2), F(7, 3)):
        for AA, BB in (((2, 2), (1, 5)), ((3, 3), (1, 10)),
                       ((F(3, 2), 4, 5), (2, F(9, 2), 7)), ((1,), (3,))):
            cases.append((sg, AA, BB))
    for sg, AA, BB in cases:
        ts = hyper_terms(sg, AA, BB, F(1), 3)   # ts[n] = coefficient of x^n
        mu, PA, PB = prods(AA, BB)
        # binomial coefficients of (1 + mu x)^{-sigma}
        b1 = -sg * mu
        b2 = (-sg) * (-sg - 1) / 2 * mu ** 2
        if ts[1] != b1:
            bad1.append((sg, AA))
        if ts[2] != b2 * PA / PB or b2 != sg * (sg + 1) / 2 * mu ** 2:
            bad2.append((sg, AA))
        # sign of the leading term of the difference
        diff2 = ts[2] - b2
        if (diff2 < 0) != (PB > PA):
            bad3.append((sg, AA))
    ck("linear-Taylor-coefficients-of-both-sides-agree--sigma*mu", not bad1,
       "%d parameter sets" % len(cases))
    ck("quadratic-Taylor-coefficient-is-sigma(sigma+1)/2*mu^2*P_A/P_B",
       not bad2, "%d parameter sets" % len(cases))
    ck("quadratic-difference-is-negative-exactly-when-P_B>P_A", not bad3,
       "%d parameter sets" % len(cases))


def check_family_products():
    """Corollary: A=(m,m), B=(1,m^2+1) has aggregate>0, mu<1 and P_B>P_A."""
    bad = []
    for m in range(2, 401):
        AA, BB = (m, m), (1, m * m + 1)
        mu, PA, PB = prods(AA, BB)
        if not (sum(BB) - sum(AA) == m * m - 2 * m + 2 > 0
                and mu == F(m * m, m * m + 1) < 1 and PB > PA):
            bad.append(m)
    _, PA2, PB2 = prods((2, 2), (1, 5))
    small = PA2 == F(9, 4) and PB2 == F(12, 5)
    big = all(prods((m, m), (1, m * m + 1))[1] <= F(16, 9) < 2
              < prods((m, m), (1, m * m + 1))[2] for m in range(3, 401))
    ck("family-A=(m,m)-B=(1,m^2+1)-satisfies-aggregate-mu<1-and-P_B>P_A",
       not bad and small and big,
       "m=2..400; at m=2, P_A=9/4 < 12/5=P_B; m>=3: P_A<=16/9<2<P_B")


def compare(sigma, A, B, x, N=60):
    """Return -1, 0 or +1 for H < , undecided, > (1+mu x)^{-sigma}.

    For sigma = p/q > 0 and H > 0 the comparison H < (1+mu x)^{-p/q} is
    equivalent to the rational inequality H^q (1+mu x)^p < 1."""
    lo, hi = hyper_enclosure(sigma, A, B, x, N)
    mu = prods(A, B)[0]
    s = F(sigma)
    p, q = s.numerator, s.denominator
    base = 1 + mu * F(x)
    if lo <= 0:
        return 0
    if hi ** q * base ** p < 1:
        return -1
    if lo ** q * base ** p > 1:
        return 1
    return 0


def check_family_violation():
    """Corollary, computed: the bound really fails for the family members."""
    x = F(1, 1000)
    bad = []
    tested = 0
    for m in range(2, 21):
        for sg in (F(1, 3), F(1, 2), F(1), F(2), F(3)):
            tested += 1
            if compare(sg, (m, m), (1, m * m + 1), x) != -1:
                bad.append((m, sg))
    ck("infinite-family-violates-the-bound-at-x=1/1000-for-every-sigma-tested",
       not bad, "%d cases: m=2..20 x sigma in {1/3,1/2,1,2,3}" % tested)


def check_positive_control():
    """Theorem 3 (b_i>a_i>0, sigma>=1, x>-1) must HOLD: if this fails the
    series evaluator, not the conjecture, is at fault."""
    bad = []
    tested = 0
    for a1 in (1, 2, 3):
        for a2 in (1, 4):
            for d1 in (1, 3):
                for d2 in (2, 5):
                    BB = (a1 + d1, a2 + d2)
                    for sg in (F(1), F(2), F(5, 2)):
                        for x in (F(1, 2), F(9, 10), F(-1, 2), F(1, 100)):
                            tested += 1
                            if compare(sg, (a1, a2), BB, x, 90) != 1:
                                bad.append((a1, a2, BB, sg, x))
    ck("Karp-Sitnik-Theorem-3-holds-on-componentwise-admissible-parameters",
       not bad, "%d cases strictly satisfy the lower bound" % tested)


def check_nonstrict_reading():
    """SCOPE OF THE REFUTATION.  A rival reading of Conjecture 1 keeps the
    componentwise clause in its NON-STRICT form, i.e. hypothesises
    b_i >= a_i > 0 together with sum(b_i - a_i) > 0 (under which the aggregate
    clause is not idle: it excludes only A = B).  Two computed facts are
    recorded here so that the transcript, and not just the prose, says what is
    and is not refuted.

    First, the exhibited object fails b_i >= a_i, so it is void under that
    rival reading.  This is not decoration: the bound for b_i >= a_i > 0,
    sigma > 0, x > 0 is the result credited to Luke in [KS], so ANY violator at
    x > 0 must fail b_i >= a_i.  An object passing b_i >= a_i while violating
    the bound would indicate a defect in this program, not a counterexample.

    Second, on samples inside that rival range (b_i >= a_i, at least one
    strict, x > 0) the bound is checked to HOLD strictly, which is what Luke's
    result requires and what makes the rival reading unrefuted here."""
    nonstrict = all(b >= a for a, b in zip(A, B))
    ck("exhibited-object-fails-b_i>=a_i-too-so-the-non-strict-reading-is-untouched",
       not nonstrict,
       "pairs (a_i,b_i) = %s; b_i>=a_i holds for %d of %d, so the object is "
       "void under the reading 'b_i>=a_i and sum(b_i-a_i)>0'"
       % (list(zip(A, B)), sum(1 for a, b in zip(A, B) if b >= a), len(A)))
    bad = []
    tested = 0
    for a1 in (1, 2, 3):
        for a2 in (2, 4):
            for d1 in (0, 1):
                for d2 in (0, 2, 5):
                    if d1 == 0 and d2 == 0:
                        continue          # A = B: the aggregate clause fails
                    BB = (a1 + d1, a2 + d2)
                    for sg in (F(1, 3), F(1), F(2)):
                        for x in (F(1, 100), F(1, 2)):
                            tested += 1
                            if compare(sg, (a1, a2), BB, x, 90) != 1:
                                bad.append((a1, a2, BB, sg, x))
    ck("bound-holds-on-non-strict-componentwise-parameters-b_i>=a_i-at-x>0",
       not bad,
       "%d cases with b_i>=a_i, sum(b_i-a_i)>0, x>0 all satisfy the bound "
       "strictly, as the result credited to Luke requires; the non-strict "
       "reading of Conjecture 1 is NOT refuted by this bundle" % tested)


def check_q1_census():
    """q=1 is exempt: there the aggregate condition IS b_1 > a_1, so the
    bound must hold on every such pair."""
    bad = []
    tested = 0
    for a in range(1, 9):
        for b in range(a + 1, 13):
            for sg in (F(1, 2), F(1), F(2)):
                xs = (F(1, 2), F(1, 100)) if sg < 1 else (F(1, 2), F(-1, 2))
                for x in xs:
                    tested += 1
                    if compare(sg, (a,), (b,), x, 90) != 1:
                        bad.append((a, b, sg, x))
    ck("no-q=1-parameter-pair-with-the-aggregate-condition-violates-the-bound",
       not bad, "%d pairs/sigma/x combinations, all satisfy the bound" % tested)


def check_supermajorization():
    """The example is not covered by weak supermajorization B <^W A."""
    sa = sorted(A)
    sb = sorted(B)
    partial = [(sum(sa[:k]), sum(sb[:k])) for k in (1, 2)]
    k1_fails = partial[0][0] > partial[0][1]
    k2_holds = partial[1][0] <= partial[1][1]
    ck("weak-supermajorization-fails-at-k=1-while-its-final-sum-holds",
       k1_fails and k2_holds,
       "k=1 needs a_1<=b_1 i.e. %d<=%d (false); k=2: %d<=%d (true)"
       % (partial[0][0], partial[0][1], partial[1][0], partial[1][1]))


def check_census():
    """Census over small positive integer parameters with q=2 at x=1/2."""
    violators = []
    undecided = 0
    cw_bad = []
    for sg in (1, 2):
        for a1 in range(1, 7):
            for a2 in range(a1, 7):
                for b1 in range(1, 7):
                    for b2 in range(1, 7):
                        AA, BB = (a1, a2), (b1, b2)
                        if sum(BB) - sum(AA) <= 0:
                            continue
                        c = compare(F(sg), AA, BB, F(1, 2), 90)
                        if c == 0:
                            undecided += 1
                        elif c == -1:
                            violators.append((sg, AA, BB))
                        if all(b > a for a, b in zip(AA, BB)) and c != 1:
                            cw_bad.append((sg, AA, BB))
    found = (SIGMA, A, B) in violators
    ck("census-x=1/2-sigma-in-{1,2}-parameters-in-1..6-reproduces-the-example",
       found and not cw_bad and undecided == 0,
       "%d violating parameter sets found, exhibited one among them; "
       "no componentwise-admissible set violates" % len(violators))
    return violators


def main():
    print("Counterexample to the aggregate-condition form of "
          "Karp-Sitnik Conjecture 1")
    print("object: 3F2(2,2,2; 1,5; -1/2) with sigma=2, A=(2,2), B=(1,5), "
          "x=1/2")
    print("CONJECTURE ACTED ON (%s), TRANSCRIBED BY HAND, NOT READ FROM THE "
          "SOURCE: \"%s\"" % (KS_CONJECTURE_SOURCE, KS_CONJECTURE_TEXT))
    print("READING ENCODED BELOW (a premise taken from the paper, not "
          "verifiable here): %s" % READING)
    print("all arithmetic is exact over the rationals; every inequality is "
          "decided by")
    print("integer cross-multiplication of certified rational enclosures")
    print("")
    # A corrupted object must be REPORTED as a failure, not raise: an
    # argument with |x| >= 1, or a non-positive denominator parameter, makes
    # the series divergent or undefined and would otherwise abort the run
    # before the verdict line is emitted.
    for fn in (check_object, check_hypotheses, check_coefficients,
               check_series_coeffs_match_terms, check_violation,
               check_closed_form, check_derivation, check_rational_chain,
               check_taylor_coefficients, check_family_products,
               check_family_violation, check_positive_control,
               check_nonstrict_reading, check_q1_census,
               check_supermajorization, check_census):
        try:
            fn()
        except Exception as exc:
            ck("no-exception-raised-by-" + fn.__name__, False,
               "%s: %s" % (type(exc).__name__, exc))
    print("")
    print("NOTE: the quantifiers 'for every sigma>0' and 'for all "
          "sufficiently small x>0' are")
    print("NOTE: infinite; they are checked here at the Taylor-coefficient "
          "level plus concrete")
    print("NOTE: rational samples (sigma in {1/3,1/2,1,2,3}, m=2..20, "
          "x=1/1000).")
    print("")
    print("NOT RE-RUN: the READING of Conjecture 1 that makes the object "
          "above a counterexample. This program decides inequalities about "
          "hypergeometric series; it cannot decide what the authors of [KS] "
          "meant. What is refuted above is the form of Conjecture 1 whose "
          "parameter hypothesis is sum(b_i-a_i)>0 together with positivity, "
          "i.e. the form in which the aggregate condition REPLACES b_i>a_i. "
          "It is NOT a refutation of the rival reading 'b_i>=a_i>0 together "
          "with sum(b_i-a_i)>0', under which the aggregate clause is not idle "
          "(it excludes only A=B) and under which the exhibited object is "
          "void, since b_1=1>=a_1=2 fails; the two checks "
          "exhibited-object-fails-b_i>=a_i-too... and "
          "bound-holds-on-non-strict-componentwise-parameters... record that "
          "shortfall as computed facts, and no counterexample to that rival "
          "reading is claimed here or in the paper.")
    print("NOT RE-RUN: the quoted source text. No external table, catalogue, "
          "preprint or journal page is consulted by this program, and none is "
          "needed for the inequality; but the paper's contribution is a "
          "reading of other authors' text, so none of the following is "
          "machine-checked here: that Conjecture 1 of [KS] reads verbatim as "
          "printed at the head of this transcript; that Conjecture 2 of [KS] "
          "lists exactly the hypotheses quoted in the paper, with no "
          "componentwise clause; that those numbers are the ones carried by "
          "the arXiv version cited; that Theorem 3 of [KS], the bound "
          "credited there to Luke, and Theorem 5 of [KP] have the hypotheses "
          "assumed above; and that [Derbazi] is the preprint cited, of the "
          "version cited, containing the quoted sentence about the broader "
          "conjectural condition at the place cited. A referee must compare "
          "every quotation, arXiv identifier and statement number in the "
          "paper against the sources by eye.")
    print("")
    n = len(CHECKS)
    k = sum(1 for _, ok in CHECKS if not ok)
    if k == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
    else:
        print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
    return 1 if k else 0


if __name__ == "__main__":
    raise SystemExit(main())

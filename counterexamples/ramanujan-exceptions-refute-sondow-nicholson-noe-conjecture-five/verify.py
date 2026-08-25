#!/usr/bin/env python3
"""Verification of counterexamples to the m>=20 clause of a conjecture on
Ramanujan primes (Sondow-Nicholson-Noe, Conjecture 5).

The conjecture clause acted on:  rho(mn) <= m*rho(n) for all m >= 20, n >= 2,
where R_n is the nth Ramanujan prime and rho(n) = pi(R_n).  The paper exhibits
ten exceptions (nine of them new) and classifies the exception set.

TAKEN FROM THE PAPER (inputs; nothing here is used to justify itself):
  - the exception set E = {22,23,24,25,26,32,33,36,37,38} at n = 9;
  - the attribution in the paper's introduction that the pair (38,9) is the one
    exception already in the literature and was there asserted to be unique.
    This is a statement ABOUT the cited literature; it is an input here and no
    check below establishes it (see NOT VERIFIED HERE);
  - the table of (m, 9m, R_{9m}, rho(9m), 20m, excess) rows;
  - R_9 = 71, rho(9) = 20; R_1244 = 24169, rho(1244) = 2690;
  - the sieve bound B = p_15672 = 172069 and the cut-offs 1244, 1245, 5225;
  - the counts #F = 21803 and the exact minimum 101/1244;
  - the analytic constants a = 1.472, b = 2.51, c = log 2, and the numerical
    inequalities quoted in the tail argument (t_0 bounds, Q(t_0) > 2.15,
    Q'(t_0) > 7.36, Q''(t_0) > 2.46, H(t_0) < -6.44, W > 1.09, gamma < 0.079).

DERIVED HERE (recomputed from first principles by this program):
  - a sieve of Eratosthenes, hence pi(x) and the kth prime;
  - the Ramanujan primes R_j from their defining property, via the
    last-deficit identity R_j = 1 + max{x : pi(x) - pi(floor(x/2)) < j},
    and rho(j) = pi(R_j);
  - for every exhibited pair: that each hypothesis of the clause holds and
    that the conclusion FAILS, by exact integer comparison;
  - an exhaustive re-search for exceptions over the paper's finite region F
    and over a strictly larger region than F;
  - the exact rational minimum of rho(n)/(2n) - 1 on 2 <= n <= 1244;
  - rigorous rational enclosures of every logarithm above, from an
    interval-arithmetic atanh series, so that no decision uses a float.

NOT VERIFIED HERE (quoted from the cited literature, not reproved):
  - the reduction of the desired inequality to the positivity of W_m(n);
  - the upper estimate rho(k) < 2k(1 + gamma(k)) for k >= 5225;
  - the bound R_j < p_{3j}, which is however confirmed for every j in range;
  - the closed forms of Q', Q'', H, u_1' and u_2', and the identity
    Q(t) = t^3 q(e^t): this program evaluates the paper's expressions in
    certified rational interval arithmetic, it does not re-derive them;
  - the historical claim that the previously known exception (38,9) was
    asserted in the literature to be the ONLY exception.  What is verified is
    that (38,9) is an exception and that nine further pairs are exceptions;
    that the earlier work claimed uniqueness cannot be checked from the files
    shipped beside this program, and no check below asserts it.
Consequently pairs with mn beyond the sieve range are not tested individually;
the program prints exactly which range was and was not re-run.

RANGES: the bounds Q''(t) > 2.46 for t >= t_0, gamma(x) < 101/1244 for
x >= 5225, and the decrease of u_1 and u_2 for t >= t_0, are needed on
unbounded ranges.  Each is established here by a monotone-envelope reduction
to finitely many rational inequalities, all of which are checked, so no such
bound rests on sampling.  Where a grid still appears (gamma's strict decrease,
Q's increase) it is corroboration only and is described as such.
"""
import sys
from fractions import Fraction

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return bool(ok)


def finish():
    n = len(CHECKS)
    k = sum(1 for _, ok in CHECKS if not ok)
    if k == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
        sys.exit(0)
    print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
    sys.exit(1)


# ----------------------------------------------------------------- paper data
PAPER_E = [22, 23, 24, 25, 26, 32, 33, 36, 37, 38]
PAPER_N_EXC = 9                 # the exceptional n
PAPER_M_FLOOR = 20              # the clause's hypothesis m >= 20
PAPER_N_FLOOR = 2               # the clause's hypothesis n >= 2
PAPER_R9 = 71
PAPER_RHO9 = 20
PAPER_PRIOR_PAIR = (38, 9)      # the one exception already in the literature.
# Whether the earlier work asserted this pair to be the ONLY exception is a
# statement about that work, not about arithmetic; it is not checked here.
# m, 9m, R_{9m}, rho(9m), 20m, excess
PAPER_TABLE = [
    (22, 198, 3167, 448, 440, 8),
    (23, 207, 3319, 467, 460, 7),
    (24, 216, 3467, 486, 480, 6),
    (25, 225, 3607, 504, 500, 4),
    (26, 234, 3761, 523, 520, 3),
    (32, 288, 4789, 644, 640, 4),
    (33, 297, 4967, 664, 660, 4),
    (36, 324, 5477, 723, 720, 3),
    (37, 333, 5651, 743, 740, 3),
    (38, 342, 5807, 762, 760, 2),
]
PAPER_TAIL_START = 1245         # analytic tail covers n >= 1245
PAPER_N_TOP = 1244              # finite part covers 2 <= n <= 1244
PAPER_K0 = 5225                 # rho(k) < 2k(1+gamma(k)) for k >= 5225
PAPER_F_SIZE = 21803            # #F
PAPER_F_NMAX = 261              # largest n occurring in F
PAPER_EPS = Fraction(101, 1244)  # min over 2<=n<=1244 of rho(n)/(2n) - 1
PAPER_EPS_ARGMIN = 1244
PAPER_R1244 = 24169
PAPER_RHO1244 = 2690
PAPER_SIEVE_INDEX = 15672       # B = p_15672
PAPER_SIEVE_B = 172069
PAPER_NMIN_IN_E = 10
PAPER_NMIN_OUT_E = 2

SIEVE_LIMIT = 800000            # our own sieve, deliberately larger than B


def sieve(limit):
    """Return (pi, primes): pi[x] = #{p <= x} for 0 <= x <= limit."""
    flag = bytearray([1]) * (limit + 1)
    flag[0] = flag[1] = 0
    i = 2
    while i * i <= limit:
        if flag[i]:
            flag[i * i::i] = bytearray(len(flag[i * i::i]))
        i += 1
    pi = [0] * (limit + 1)
    primes = []
    c = 0
    for x in range(limit + 1):
        if flag[x]:
            c += 1
            primes.append(x)
        pi[x] = c
    return pi, primes


def ramanujan_primes(pi, limit, jmax):
    """R_j for 1 <= j <= jmax, from the definition.

    d(x) = pi(x) - pi(floor(x/2)) is constant on [k, k+1) for integer k, so the
    real-variable condition pi(x) - pi(x/2) >= j for all x >= X is equivalent to
    d(y) >= j for all integers y >= X.  With g(x) = min{d(y) : x <= y < limit},
    which is non-decreasing, R_j is the least x with g(x) >= j; equivalently
    R_j = 1 + max{x < limit : d(x) < j}.
    """
    d = [0] * limit
    for x in range(limit):
        d[x] = pi[x] - pi[x >> 1]
    g = [0] * limit
    cur = d[limit - 1]
    for x in range(limit - 1, -1, -1):
        if d[x] < cur:
            cur = d[x]
        g[x] = cur
    R = [0] * (jmax + 1)
    x = 0
    for j in range(1, jmax + 1):
        while g[x] < j:
            x += 1
            if x >= limit:
                raise RuntimeError("sieve window too small for j=%d" % j)
        R[j] = x
    return R


# --------------------------------------- rigorous rational logarithm enclosures
# Every quantity below is a rational with denominator dividing GRAIN; rounding is
# always outward, so an enclosure never shrinks and no bound is ever optimistic.
GRAIN = 10 ** 40
TERMS = 35


def _fl(x):
    """Largest multiple of 1/GRAIN that is <= x."""
    return Fraction(x.numerator * GRAIN // x.denominator, GRAIN)


def _ce(x):
    """Smallest multiple of 1/GRAIN that is >= x."""
    return -_fl(-x)


def _atanh2_lower(z, terms=TERMS):
    """Certified lower bound for 2*atanh(z), 0 <= z <= 1/3."""
    s = Fraction(0)
    zz = _fl(z * z)
    p = z
    for i in range(terms):
        s = _fl(s + p / (2 * i + 1))
        p = _fl(p * zz)
    return _fl(2 * s)


def _atanh2_upper(z, terms=TERMS):
    """Certified upper bound for 2*atanh(z), 0 <= z <= 1/3."""
    s = Fraction(0)
    zz = _ce(z * z)
    p = z
    for i in range(terms):
        s = _ce(s + p / (2 * i + 1))
        p = _ce(p * zz)
    tail = _ce(2 * p / ((2 * terms + 1) * (1 - zz)))
    return _ce(2 * s + tail)


_LN2_LO = _atanh2_lower(Fraction(1, 3))
_LN2_HI = _atanh2_upper(Fraction(1, 3))


def ln_bounds(q):
    """(lo, hi) rationals with lo <= ln(q) <= hi, for rational q > 0.

    Uses ln q = k*ln 2 + 2*atanh((y-1)/(y+1)) with y = q/2^k in [1,2), so the
    series argument satisfies 0 <= z <= 1/3 and converges geometrically.
    """
    q = Fraction(q)
    if q <= 0:
        raise ValueError("ln of non-positive")
    k = 0
    while q >= 2:
        q /= 2
        k += 1
    while q < 1:
        q *= 2
        k -= 1
    ylo, yhi = _fl(q), _ce(q)
    lo = _atanh2_lower(_fl((ylo - 1) / (ylo + 1)))
    hi = _atanh2_upper(_ce((yhi - 1) / (yhi + 1)))
    if k >= 0:
        return _fl(lo + k * _LN2_LO), _ce(hi + k * _LN2_HI)
    return _fl(lo + k * _LN2_HI), _ce(hi + k * _LN2_LO)


class Iv(object):
    """Rational interval; every arithmetic result encloses the true value."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        a = Fraction(lo)
        b = a if hi is None else Fraction(hi)
        if a > b:
            raise ValueError("empty interval")
        self.lo = _fl(a)          # outward rounding keeps denominators bounded
        self.hi = _ce(b)          # and can only widen the enclosure

    def __add__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    def __radd__(self, o):
        return Iv(o) + self

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Iv) else Iv(o)))

    def __rsub__(self, o):
        return Iv(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv(min(p), max(p))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        if o.lo <= 0 <= o.hi:
            raise ZeroDivisionError("interval spans zero")
        return self * Iv(1 / o.hi, 1 / o.lo)

    def __rtruediv__(self, o):
        return Iv(o) / self

    def __pow__(self, k):
        r = Iv(1)
        for _ in range(k):
            r = r * self
        return r

    def gt(self, v):
        return self.lo > Fraction(v)

    def lt(self, v):
        return self.hi < Fraction(v)

    def __repr__(self):
        return "[%.6f,%.6f]" % (float(self.lo), float(self.hi))


def iv_ln(x):
    """Enclosure of ln on an interval of positive rationals."""
    if x.lo <= 0:
        raise ValueError("ln of non-positive interval")
    return Iv(ln_bounds(x.lo)[0], ln_bounds(x.hi)[1])


KNOWN_PI = [(10, 4), (100, 25), (1000, 168), (10000, 1229), (100000, 9592)]
KNOWN_R = [2, 11, 17, 29, 41, 47, 59, 67, 71, 97]   # R_1 .. R_10, literature


def check_sieve(pi, primes):
    bad = [(x, w, pi[x]) for x, w in KNOWN_PI if pi[x] != w]
    ck("sieve_reproduces_known_pi_values", not bad,
       "pi(10^5)=%d, %d anchors" % (pi[100000], len(KNOWN_PI)))
    okb = (len(primes) >= PAPER_SIEVE_INDEX
           and primes[PAPER_SIEVE_INDEX - 1] == PAPER_SIEVE_B
           and pi[PAPER_SIEVE_B] == PAPER_SIEVE_INDEX)
    ck("paper_sieve_bound_B_equals_p_15672", okb,
       "p_%d=%d, pi(%d)=%d" % (PAPER_SIEVE_INDEX, primes[PAPER_SIEVE_INDEX - 1],
                               PAPER_SIEVE_B, pi[PAPER_SIEVE_B]))


def check_ramanujan(pi, primes, R, jmax):
    got = R[1:11]
    ck("ramanujan_primes_match_known_initial_segment", got == KNOWN_R,
       "R_1..R_10=%s" % (got,))
    isp = set(primes)
    bad_p = [j for j in range(1, jmax + 1) if R[j] not in isp]
    bad_i = [j for j in range(2, jmax + 1) if R[j] <= R[j - 1]]
    ck("every_computed_R_j_is_prime_and_strictly_increasing",
       not bad_p and not bad_i,
       "j<=%d, R_%d=%d" % (jmax, jmax, R[jmax]))
    bad_b = [j for j in range(1, jmax + 1) if R[j] >= primes[3 * j - 1]]
    ck("computed_R_j_below_p_3j_as_the_cited_bound_requires", not bad_b,
       "%d values of j, tightest at j=%d" %
       (jmax, max(range(1, jmax + 1),
                  key=lambda j: Fraction(R[j], primes[3 * j - 1]))))


def check_exhibited_object(R, rho):
    tbl_m = [row[0] for row in PAPER_TABLE]
    ok = (len(PAPER_E) == 10
          and sorted(set(PAPER_E)) == PAPER_E
          and all(isinstance(v, int) for v in PAPER_E)
          and tbl_m == PAPER_E)
    ck("exhibited_exception_set_is_well_formed", ok,
       "|E|=%d, E=%s, n=%d" % (len(set(PAPER_E)), PAPER_E, PAPER_N_EXC))
    ck("R_9_and_rho_9_recomputed_from_the_sieve",
       R[PAPER_N_EXC] == PAPER_R9 and rho[PAPER_N_EXC] == PAPER_RHO9,
       "R_9=%d, rho(9)=%d" % (R[PAPER_N_EXC], rho[PAPER_N_EXC]))
    bad = []
    for mm, prod, rr, rh, mrho, exc in PAPER_TABLE:
        k = PAPER_N_EXC * mm
        if (prod != k or R[k] != rr or rho[k] != rh
                or mrho != mm * rho[PAPER_N_EXC]
                or exc != rho[k] - mm * rho[PAPER_N_EXC]):
            bad.append(mm)
    ck("every_table_entry_reproduced_by_independent_computation", not bad,
       "10 rows; first R_198=%d rho=%d, last R_342=%d rho=%d"
       % (R[198], rho[198], R[342], rho[342]))


def check_hypotheses(rho):
    """Every hypothesis of the clause rho(mn) <= m*rho(n), m>=20, n>=2."""
    fails = []
    for mm in PAPER_E:
        if not (isinstance(mm, int) and mm >= PAPER_M_FLOOR):
            fails.append(("m", mm))
        if not (isinstance(PAPER_N_EXC, int) and PAPER_N_EXC >= PAPER_N_FLOOR):
            fails.append(("n", mm))
        if rho[PAPER_N_EXC * mm] == 0 or rho[PAPER_N_EXC] == 0:
            fails.append(("defined", mm))
    ck("all_pairs_satisfy_every_hypothesis_of_the_clause", not fails,
       "m in [%d,%d] all >= %d; n = %d >= %d"
       % (min(PAPER_E), max(PAPER_E), PAPER_M_FLOOR, PAPER_N_EXC, PAPER_N_FLOOR))


def check_conclusion_violated(rho):
    """Load-bearing: the conclusion is FALSE at each exhibited pair."""
    excess = []
    for mm in PAPER_E:
        excess.append(rho[PAPER_N_EXC * mm] - mm * rho[PAPER_N_EXC])
    ok = all(e > 0 for e in excess)
    ck("conclusion_rho_mn_le_m_rho_n_fails_at_every_exhibited_pair", ok,
       "rho(9m) - m*rho(9) = %s (all > 0)" % (excess,))
    paper_exc = [row[5] for row in PAPER_TABLE]
    ck("computed_excesses_agree_with_the_stated_excess_column",
       excess == paper_exc, "min excess %d at m=%d"
       % (min(excess), PAPER_E[excess.index(min(excess))]))


def check_prior_pair_and_nine_more(rho):
    """The pair already in the literature is an exception, and so are nine more.

    Purely arithmetic content.  This check says nothing about what the earlier
    work claimed: the assertion that (38,9) was there declared unique is an
    input taken from the paper and is disclosed as unverified in the closing
    scope line.
    """
    exc = [(mm, PAPER_N_EXC) for mm in PAPER_E
           if rho[PAPER_N_EXC * mm] > mm * rho[PAPER_N_EXC]]
    extra = [p for p in exc if p != PAPER_PRIOR_PAIR]
    ok = (PAPER_PRIOR_PAIR in exc and len(extra) == 9)
    ck("prior_pair_is_an_exception_and_nine_further_exceptions_exist", ok,
       "verified exceptions %d, additional to (%d,%d): %d; the literature's "
       "uniqueness assertion is an unverified input, not a check"
       % (len(exc), PAPER_PRIOR_PAIR[0], PAPER_PRIOR_PAIR[1], len(extra)))


def region_F_pairs():
    for n in range(PAPER_N_FLOOR, PAPER_N_TOP + 1):
        top = (PAPER_K0 - 1) // n
        for mm in range(PAPER_M_FLOOR, top + 1):
            yield mm, n


def check_region_F(rho):
    pairs = list(region_F_pairs())
    formula = sum((PAPER_K0 - 1) // n - (PAPER_M_FLOOR - 1)
                  for n in range(PAPER_N_FLOOR, PAPER_F_NMAX + 1))
    nmax = max(n for _, n in pairs)
    ck("cardinality_of_the_finite_region_F_is_21803",
       len(pairs) == PAPER_F_SIZE and formula == PAPER_F_SIZE
       and nmax == PAPER_F_NMAX,
       "enumerated %d, closed form %d, max n %d" % (len(pairs), formula, nmax))
    found = sorted(p for p in pairs if rho[p[0] * p[1]] > p[0] * rho[p[1]])
    want = sorted((mm, PAPER_N_EXC) for mm in PAPER_E)
    ck("exhaustive_search_of_F_finds_exactly_the_claimed_exceptions",
       found == want, "%d pairs tested, %d exceptions, all at n=%d"
       % (len(pairs), len(found), PAPER_N_EXC))
    return want


def check_extended_census(rho, jmax, want):
    """Search a region strictly larger than F: every m>=20, n>=2 with mn<=jmax."""
    found = []
    tested = 0
    for n in range(PAPER_N_FLOOR, jmax // PAPER_M_FLOOR + 1):
        rn = rho[n]
        for mm in range(PAPER_M_FLOOR, jmax // n + 1):
            tested += 1
            if rho[mm * n] > mm * rn:
                found.append((mm, n))
    ck("no_further_exception_with_mn_up_to_the_full_sieve_range",
       sorted(found) == sorted(want),
       "%d pairs tested, mn <= %d (F needs only mn < %d), %d exceptions"
       % (tested, jmax, PAPER_K0, len(found)))


def check_epsilon(R, rho):
    best = None
    arg = None
    ties = 0
    for n in range(PAPER_N_FLOOR, PAPER_N_TOP + 1):
        v = Fraction(rho[n], 2 * n) - 1
        if best is None or v < best:
            best, arg, ties = v, n, 1
        elif v == best:
            ties += 1
    ck("exact_minimum_of_rho_n_over_2n_minus_one_is_101_over_1244",
       best == PAPER_EPS and arg == PAPER_EPS_ARGMIN and ties == 1,
       "min = %s at n = %d, %d minimiser(s)" % (best, arg, ties))
    ck("R_1244_and_rho_1244_recomputed",
       R[PAPER_N_TOP] == PAPER_R1244 and rho[PAPER_N_TOP] == PAPER_RHO1244,
       "R_1244=%d, rho(1244)=%d" % (R[PAPER_N_TOP], rho[PAPER_N_TOP]))


def gamma_iv(x):
    """Enclosure of gamma(x) = (c + c/L + 0.565/L^2)/(L + log L - c - c/L)."""
    c = Iv(*ln_bounds(2))
    L = iv_ln(Iv(Fraction(x)))
    num = c + c / L + Iv(Fraction(565, 1000)) / (L ** 2)
    den = L + iv_ln(L) - c - c / L
    return num / den


def check_gamma(rho):
    g = gamma_iv(PAPER_K0)
    ck("gamma_at_5225_is_below_the_stated_bound_and_below_101_over_1244",
       g.lt(Fraction(79, 1000)) and g.lt(PAPER_EPS),
       "gamma(%d) in %s < 0.079 < %s" % (PAPER_K0, g, PAPER_EPS))
    # the paper's own rational chain, with its three logarithm facts checked
    L = Iv(*ln_bounds(PAPER_K0))
    facts = (L.gt(Fraction(856, 100))
             and iv_ln(L).gt(Fraction(214, 100))
             and Iv(*ln_bounds(2)).lt(Fraction(694, 1000)))
    l8 = Fraction(856, 100)
    chain = ((Fraction(694, 1000) + Fraction(694, 1000) / l8
              + Fraction(565, 1000) / l8 ** 2)
             / (l8 + Fraction(214, 100) - Fraction(694, 1000)
                - Fraction(694, 1000) / l8))
    ck("logarithm_facts_backing_the_stated_gamma_chain_hold",
       facts and chain < Fraction(79, 1000) < PAPER_EPS,
       "log 5225 > 8.56, loglog 5225 > 2.14, log 2 < 0.694, chain = %s"
       % float(chain))
    # Monotone envelope over the whole unbounded range x >= 5225, so that the
    # bound below is not a grid claim.  Write L = log x, L0 = log 5225.  For
    # x >= 5225 we have L >= L0 > 1; since c > 0, each numerator term
    # (c, c/L, 0.565/L^2) is non-increasing in L and each denominator term
    # (L, log L, -c/L) is non-decreasing, so with den0 > 0
    #        gamma(x) = num(L)/den(L) <= num(L0)/den(L0).
    # Every rational inequality this reduction needs (c > 0, L0 > 1, num0 > 0,
    # den0 > 0, num0/den0 < 101/1244) is checked below.
    L0 = Iv(*ln_bounds(PAPER_K0))
    num0 = C_IV + C_IV / L0 + Iv(Fraction(565, 1000)) / (L0 ** 2)
    den0 = L0 + iv_ln(L0) - C_IV - C_IV / L0
    env = num0 / den0
    envelope = (C_IV.gt(0) and L0.gt(1) and num0.gt(0) and den0.gt(0)
                and env.lt(PAPER_EPS))
    grid = [PAPER_K0]
    while grid[-1] < 10 ** 18:
        grid.append(grid[-1] * 4)
    vals = [gamma_iv(x) for x in grid]
    mono = all(vals[i + 1].hi < vals[i].lo for i in range(len(vals) - 1))
    small = all(v.lt(PAPER_EPS) for v in vals)
    ck("gamma_below_101_over_1244_for_every_x_at_least_5225_by_envelope",
       envelope and mono and small,
       "envelope sup over all x >= %d is %s < %s (monotone reduction, not "
       "sampling); strict decrease corroborated at %d grid points %d..%d, "
       "gamma there %s..%s"
       % (PAPER_K0, env, PAPER_EPS, len(grid), grid[0], grid[-1],
          vals[0], vals[-1]))
    bound = g.hi
    bad = [n for n in range(PAPER_N_FLOOR, PAPER_N_TOP + 1)
           if not Fraction(rho[n]) > 2 * n * (1 + bound)]
    ck("decisive_inequality_2mn_times_1_plus_gamma_stays_under_m_rho_n",
       not bad, "rho(n) > 2n(1+%s) for all %d values of n; the constant is the "
       "envelope sup of gamma on x >= %d, so this covers every mn >= %d"
       % (float(bound), PAPER_N_TOP - PAPER_N_FLOOR + 1, PAPER_K0, PAPER_K0))


A_IV = Iv(Fraction(1472, 1000))     # a = 1.472
B_IV = Iv(Fraction(251, 100))       # b = 2.51
C_IV = Iv(*ln_bounds(2))            # c = log 2


def phi_iv(t):
    """Phi(t) = a log t + b."""
    return A_IV * iv_ln(t) + B_IV


def q_iv(n):
    """q(n) = 0.03 + c/log n - phi(n)/log^2 n - phi(n)/log^3 n."""
    L = iv_ln(Iv(Fraction(n)))
    ph = A_IV * iv_ln(L) + B_IV
    return Iv(Fraction(3, 100)) + C_IV / L - ph / (L ** 2) - ph / (L ** 3)


def Q_iv(t):
    return (Iv(Fraction(3, 100)) * t ** 3 + C_IV * t ** 2
            - phi_iv(t) * (t + Iv(1)))


def check_tail_constants():
    t0 = Iv(*ln_bounds(PAPER_TAIL_START))
    lt0 = iv_ln(t0)
    ok = (t0.gt(Fraction(7126, 1000)) and t0.lt(Fraction(7127, 1000))
          and C_IV.gt(Fraction(6931, 10000)) and C_IV.lt(Fraction(6932, 10000))
          and lt0.gt(Fraction(1963, 1000)) and lt0.lt(Fraction(1964, 1000)))
    ck("quoted_enclosures_of_t0_log2_and_log_t0_are_correct", ok,
       "t0 in %s, c in %s, log t0 in %s" % (t0, C_IV, lt0))
    return t0


def check_tail_convexity(t0):
    q = Q_iv(t0)
    qp = (Iv(Fraction(9, 100)) * t0 ** 2 + 2 * C_IV * t0
          - A_IV * (t0 + Iv(1)) / t0 - phi_iv(t0))
    qpp = (Iv(Fraction(18, 100)) * t0 + 2 * C_IV
           - A_IV / t0 + A_IV / t0 ** 2)
    # The paper needs Q'' > 2.46 on the whole ray t >= t_0, not just at t_0.
    # For t >= t_0 > 0:  0.18t >= 0.18t_0,  -a/t >= -a/t_0,  a/t^2 > 0, so
    #        Q''(t) > 0.18 t_0 + 2c - a/t_0 =: qpp_ray,
    # a single rational quantity; qpp_ray > 2.46 therefore settles the ray.
    qpp_ray = Iv(Fraction(18, 100)) * t0 + 2 * C_IV - A_IV / t0
    ok = (q.gt(Fraction(215, 100)) and qp.gt(Fraction(736, 100))
          and qpp.gt(Fraction(246, 100)) and qpp_ray.gt(Fraction(246, 100))
          and A_IV.gt(0) and t0.gt(2))
    ck("Q_Q_prime_at_t0_and_Q_double_prime_on_all_t_ge_t0_meet_the_bounds", ok,
       "Q=%s>2.15, Q'=%s>7.36, Q''(t0)=%s>2.46, and inf over t>=t0 of Q'' "
       "is >= %s > 2.46" % (q, qp, qpp, qpp_ray))
    grid = [t0] + [Iv(Fraction(v)) for v in (8, 9, 10, 12, 15, 20, 30, 45, 69)]
    vals = [Q_iv(t) for t in grid]
    ck("Q_is_positive_and_increasing_across_the_tail_grid",
       all(v.gt(0) for v in vals)
       and all(vals[i + 1].lo > vals[i].hi for i in range(len(vals) - 1)),
       "%d points t0..69, Q from %s to %s" % (len(grid), vals[0], vals[-1]))
    qn = q_iv(PAPER_TAIL_START)
    ck("q_is_positive_at_the_tail_threshold_1245", qn.gt(0),
       "q(%d) in %s" % (PAPER_TAIL_START, qn))


def check_tail_W(t0):
    lt0 = iv_ln(t0)
    ph = phi_iv(t0)
    H = A_IV * (lt0 - C_IV) + ph - 2 * ph * (lt0 - C_IV)
    hp_num = A_IV - A_IV * (lt0 - C_IV) - ph          # sign of H'(t0)
    u1 = ph / t0
    u2 = ph * (lt0 - C_IV) / t0 ** 2
    # Both decreases are needed on the whole ray t >= t_0, and both reduce to
    # rational inequalities at t_0.  The numerator of H' rearranges to
    # a(1+c) - b - 2a log t, which is non-increasing in t when a > 0, so it
    # stays <= its value at t_0 (< 0) and H stays <= H(t_0) < -6.44 < 0,
    # giving u_2' = H/t^3 < 0 throughout.  And u_1'(t) = (a - Phi(t))/t^2 < 0
    # throughout because Phi is increasing with Phi(t_0) > a.
    ck("H_negative_and_u1_u2_decreasing_on_all_of_t_ge_t0",
       H.lt(Fraction(-644, 100)) and hp_num.lt(0) and A_IV.gt(0)
       and (ph - A_IV).gt(0) and u1.gt(0) and u2.gt(0),
       "H=%s < -6.44, H' numerator %s < 0 and non-increasing in t, "
       "Phi(t0)=%s > a so u1 decreases" % (H, hp_num, ph))
    l20 = Iv(*ln_bounds(20))
    W = (C_IV - Iv(Fraction(3, 100))) * l20 - u1 - u2
    dropped = ((C_IV * l20 - Iv(Fraction(565, 1000))).gt(0)
               and ((l20 + Iv(Fraction(1, 2))) * C_IV
                    - Iv(Fraction(565, 1000))).gt(0))
    ck("tail_lower_bound_for_W_exceeds_1_09_with_dropped_terms_positive",
       W.gt(Fraction(109, 100)) and dropped,
       "W > %s, u1=%s, u2=%s" % (W, u1, u2))


def check_nmin(rho, jmax):
    """N_min(m) = min{N>=2 : no exception at any n >= N}, over the tested range."""
    bad = []
    mtop = jmax // PAPER_N_FLOOR
    roomy = 0
    for mm in range(PAPER_M_FLOOR, mtop + 1):
        ntop = jmax // mm
        if ntop > PAPER_NMIN_IN_E:
            roomy += 1
        exc = [n for n in range(PAPER_N_FLOOR, ntop + 1)
               if rho[mm * n] > mm * rho[n]]
        want = PAPER_NMIN_IN_E if mm in PAPER_E else PAPER_NMIN_OUT_E
        got = (max(exc) + 1) if exc else PAPER_N_FLOOR
        if got != want:
            bad.append((mm, got, want))
    ck("least_valid_threshold_N_min_is_10_on_E_and_2_off_E", not bad,
       "m from %d to %d, each n <= %d/m; %d values of m reach past n=10; "
       "%d disagreements" % (PAPER_M_FLOOR, mtop, jmax, roomy, len(bad)))


def main():
    pi, primes = sieve(SIEVE_LIMIT)
    jmax = len(primes) // 3          # largest j with p_3j inside the sieve
    R = ramanujan_primes(pi, SIEVE_LIMIT, jmax)
    rho = [0] * (jmax + 1)
    for j in range(1, jmax + 1):
        rho[j] = pi[R[j]]
    print("sieve limit %d, pi = %d, R_j computed for j <= %d (R_%d = %d)"
          % (SIEVE_LIMIT, len(primes), jmax, jmax, R[jmax]))
    print("exhibited object: %d pairs (m,%d) with m in %s"
          % (len(PAPER_E), PAPER_N_EXC, PAPER_E))
    print("decoded rows: " + "; ".join(
        "R_%d=%d rho=%d vs %d*rho(9)=%d"
        % (PAPER_N_EXC * mm, R[PAPER_N_EXC * mm], rho[PAPER_N_EXC * mm],
           mm, mm * rho[PAPER_N_EXC]) for mm in PAPER_E))

    check_sieve(pi, primes)
    check_ramanujan(pi, primes, R, jmax)
    check_exhibited_object(R, rho)
    check_hypotheses(rho)
    check_conclusion_violated(rho)
    check_prior_pair_and_nine_more(rho)
    want = check_region_F(rho)
    check_extended_census(rho, jmax, want)
    check_epsilon(R, rho)
    check_gamma(rho)
    t0 = check_tail_constants()
    check_tail_convexity(t0)
    check_tail_W(t0)
    check_nmin(rho, jmax)

    print("NOT RE-RUN HERE: (i) pairs with mn > %d are not tested one by one; "
          "that range is covered by the two analytic ingredients above, whose "
          "quoted numerical inequalities are verified but whose reduction to "
          "W_m(n) and whose estimate rho(k) < 2k(1+gamma(k)) for k >= %d are "
          "quoted from the cited literature and are not reproved here; nor are "
          "the closed forms of Q', Q'', H, u_1', u_2' or the identity "
          "Q(t) = t^3 q(e^t), which this program evaluates in certified "
          "rational interval arithmetic but does not re-derive. "
          "(ii) Three bounds the argument needs on unbounded ranges -- "
          "Q'' > 2.46 for t >= t_0, gamma(x) < 101/1244 for x >= %d, and the "
          "decrease of u_1 and u_2 for t >= t_0 -- are established here for "
          "those full ranges, by monotone-envelope reductions to finitely many "
          "rational inequalities that are all checked; the grids that remain "
          "(gamma's strict decrease at 25 geometric points up to about "
          "1.4e18, Q's increase at 10 points) corroborate only and carry no "
          "load. (iii) Nothing here verifies the paper's attribution that the "
          "pair (38,9) was asserted in the earlier literature to be the unique "
          "exception: that (38,9) is an exception and that nine further pairs "
          "are exceptions is verified by exact integer arithmetic, but the "
          "content of the cited claim is an input taken from the paper and is "
          "not checkable from the files shipped here."
          % (jmax, PAPER_K0, PAPER_K0))
    finish()


if __name__ == "__main__":
    main()

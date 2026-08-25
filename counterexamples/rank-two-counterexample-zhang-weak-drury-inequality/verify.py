#!/usr/bin/env python3
"""Verification of a rank-two counterexample to the weak permanent inequality
Phi_1(A) <= per(A) for Hermitian positive semidefinite A, where

    Phi_1(A) = sum_{k=2}^{n} |a_{1k}|^2 per A[hat{1,k}]

and A[hat{1,k}] deletes rows and columns 1 and k (per of the empty matrix = 1).

TAKEN FROM THE PAPER (inputs; transcribed literally, nothing else):
  * the twelve Gaussian-integer columns of the 2 x 12 matrix V, whose Gram
    matrix A = V^* V is the exhibited object;
  * the claimed order 12 and rank 2;
  * the claimed values per(A), Phi_1(A) and their difference;
  * the table of the eleven weights |a_{1k}|^2 and the eleven minor permanents
    per A[hat{1,k}], including the grouping of equal entries;
  * the claimed band 4 <= n_min <= 12 for the least order of a counterexample;
  * two rationals the paper reports for a matrix printed elsewhere.

DERIVED HERE (recomputed from V alone; no paper value is used as an input to
any computation, only as a target to compare against):
  * A itself, its Hermitian symmetry, membership in Z[i], positive diagonal;
  * the exact rank by elimination over Q(i), and positive semidefiniteness by
    an exact principal-minor certificate that does not use the factorization;
  * per(A) by Ryser's formula over exact Gaussian integers, cross-checked
    against the defining permutation sum on small orders and against the
    paper's rank-two coefficient identity;
  * all eleven minor permanents and hence Phi_1(A), by two independent routes;
  * the sign and exact size of Phi_1(A) - per(A), the load-bearing quantity;
  * the rank-one identity Phi_1 = per/n, and the order-2 and order-3 gap
    identities, each proved on a grid large enough to force the polynomial
    identity rather than merely sampled;
  * exhaustively, that no Gram matrix of three vectors from a named finite set
    violates the inequality, so the paper's lower bound is reproduced;
  * the diagonal-scaling homogeneity of per and of Phi_1 -- proved on a full
    grid at orders 2 and 3, and forced coordinatewise (three distinct values of
    each d_j, with both sides of degree at most 2 in it) for the exhibited
    order-12 matrix along the two rational chains actually used, the full
    twelve-variable identity being out of reach -- the direct-sum step at the
    orders listed by DIRECT_SUM_PADDINGS, and an explicit positive definite
    perturbation of the order-12 matrix and of every direct sum formed here,
    which together give the consequences;
  * the claim that the order-(n-1)-minor variant of the second term is a
    different, nonhomogeneous quantity: it scales by d_1^2 prod_j d_j^2 under a
    diagonal scaling where the permanent scales by prod_j d_j^2;
  * negative controls that damage the object and require the certificates to
    reject it;
  * a bounded sweep of rank-two configurations of orders 4 to 11.

NOT RE-RUN (each printed again at the point of use):
  * the paper's ATTRIBUTION of the refuted inequality -- Zhang's Eq. (21) at
    p. 315, its status as the weak consequence of Eq. (20) surviving
    Hutchinson's disproof, and Hutchinson's indexing B_{k-1,k-1} =
    A[hat{1,k}] -- is transcribed from the cited articles and cannot be
    checked here: there is no network access and no copy of those articles, so
    what is refuted is the inequality as printed in this paper;
  * no exhaustive census over all positive semidefinite matrices of orders 4 to
    11 is attempted;
  * the external matrix behind the two reported rationals is not restated in
    this paper, so only the printed ordering of those two rationals is checked;
  * the corollary's CORRELATION-matrix form is not exhibited: the normalization
    d_j = a_jj^(-1/2) is irrational for the diagonal of the exhibited A, so no
    unit-diagonal matrix can be built in the exact rational arithmetic used
    here, and the claim "in every order N >= 12" is recomputed only at the
    finitely many orders 12 + DIRECT_SUM_PADDINGS (the list is printed, not
    described, at the point of use), the rest resting on the paper's induction.

STANDARD LIBRARY ONLY.  Every decision uses exact integer or Fraction
arithmetic; no floating point anywhere.
"""

import sys
from fractions import Fraction
from itertools import combinations, permutations

CHECKS = []


class Cx(object):
    """Exact complex number over Z or Q (components int or Fraction)."""

    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re = re
        self.im = im

    def __add__(self, o):
        return Cx(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return Cx(self.re - o.re, self.im - o.im)

    def __neg__(self):
        return Cx(-self.re, -self.im)

    def __mul__(self, o):
        if isinstance(o, Cx):
            return Cx(self.re * o.re - self.im * o.im,
                      self.re * o.im + self.im * o.re)
        return Cx(self.re * o, self.im * o)

    def __truediv__(self, o):
        d = o.re * o.re + o.im * o.im
        num = self * o.conj()
        return Cx(Fraction(num.re, 1) / d, Fraction(num.im, 1) / d)

    def conj(self):
        return Cx(self.re, -self.im)

    def norm2(self):
        return self.re * self.re + self.im * self.im

    def is_zero(self):
        return self.re == 0 and self.im == 0

    def __eq__(self, o):
        return self.re == o.re and self.im == o.im

    def __repr__(self):
        if self.im == 0:
            return str(self.re)
        return "%s%s%si" % (self.re, "+" if self.im > 0 else "-", abs(self.im))


ONE = Cx(1, 0)
ZERO = Cx(0, 0)


def permanent(M):
    """Ryser's formula with Gray-code subset enumeration.  Exact."""
    n = len(M)
    if n == 0:
        return Cx(1, 0)
    s = [Cx(0, 0) for _ in range(n)]
    total = Cx(0, 0)
    prev = 0
    pop = 0
    for k in range(1, 1 << n):
        g = k ^ (k >> 1)
        diff = g ^ prev
        j = diff.bit_length() - 1
        if g & diff:
            pop += 1
            for i in range(n):
                s[i] = s[i] + M[i][j]
        else:
            pop -= 1
            for i in range(n):
                s[i] = s[i] - M[i][j]
        prev = g
        prod = s[0]
        for i in range(1, n):
            if prod.is_zero():
                break
            prod = prod * s[i]
        if not prod.is_zero():
            if pop % 2 == 0:
                total = total + prod
            else:
                total = total - prod
    if n % 2:
        total = -total
    return total


def gram(cols):
    """A = V^* V for V whose columns are the given 2-vectors (as Cx pairs)."""
    m = len(cols)
    return [[cols[j][0].conj() * cols[k][0] + cols[j][1].conj() * cols[k][1]
             for k in range(m)] for j in range(m)]


def drop_rows_cols(M, idxs):
    keep = [i for i in range(len(M)) if i not in idxs]
    return [[M[i][j] for j in keep] for i in keep]


ANOMALIES = []


def phi1(M):
    """Phi_1(A) = sum_{k>=2} |a_{1k}|^2 per A[hat{1,k}], 1-based k.  For
    Hermitian M every minor permanent is real; any exception is recorded and
    reported by a check of its own rather than silently discarded."""
    n = len(M)
    tot = 0
    for k in range(1, n):
        w = M[0][k].norm2()
        if w == 0:
            continue
        p = permanent(drop_rows_cols(M, {0, k}))
        if p.im != 0:
            ANOMALIES.append((n, k, p.im))
        tot += w * p.re
    return tot


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    return bool(ok)


def finish():
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


# ----------------------------------------------------------------------
# VALUES TAKEN FROM THE PAPER
# ----------------------------------------------------------------------
# The two rows of V, as (real, imaginary) integer pairs.
V_ROW_X = [(1, 0), (1, 0), (7, 0), (7, 0), (9, 0), (5, 0),
           (5, 0), (7, 0), (9, 0), (5, 0), (5, 0), (0, 0)]
V_ROW_Y = [(0, 0), (-2, 0), (5, 0), (-4, 3), (2, -6), (8, 6),
           (-3, 9), (-4, -3), (2, 6), (-3, -9), (8, -6), (1, 0)]

PAPER_ORDER = 12
PAPER_RANK = 2
PAPER_PER = 463995226880175960000000
PAPER_PHI1 = 464101796319782967000000
PAPER_GAP = 106569439607007000000

# The table in the paper: k (1-based) -> (|a_{1k}|^2, per A[hat{1,k}]).
PAPER_TABLE = {
    2: (1, 9076601601430046100000),
    3: (49, 1684042529828861700000),
    4: (49, 1687819148108486100000),
    8: (49, 1687819148108486100000),
    5: (81, 1042294059228024900000),
    9: (81, 1042294059228024900000),
    6: (25, 362453951893503780000),
    11: (25, 362453951893503780000),
    7: (25, 402529981049836740000),
    10: (25, 402529981049836740000),
}
PAPER_ZERO_WEIGHT_INDEX = 12  # a_{1,12} = 0

# Paper's stated bounds on the least order of a counterexample.
PAPER_NMIN_LOWER = 4
PAPER_NMIN_UPPER = 12


def build_columns():
    return [(Cx(*V_ROW_X[j]), Cx(*V_ROW_Y[j])) for j in range(len(V_ROW_X))]


def check_exhibited_object(cols):
    """Decode V, count it, print it back."""
    m = len(cols)
    all_int = all(isinstance(c.re, int) and isinstance(c.im, int)
                  for col in cols for c in col)
    e1 = cols[0][0] == ONE and cols[0][1] == ZERO
    e2 = cols[-1][0] == ZERO and cols[-1][1] == ONE
    print("  V row 1: " + " ".join(repr(c[0]) for c in cols))
    print("  V row 2: " + " ".join(repr(c[1]) for c in cols))
    ck("V_is_2_by_%d_gaussian_integer" % PAPER_ORDER,
       m == PAPER_ORDER and all_int, "columns=%d" % m)
    ck("V_first_and_last_columns_are_e1_e2", e1 and e2,
       "col1=%s,%s col%d=%s,%s" % (cols[0][0], cols[0][1], m,
                                   cols[-1][0], cols[-1][1]))


def check_gram_shape(A):
    n = len(A)
    herm = all(A[j][k] == A[k][j].conj() for j in range(n) for k in range(n))
    integral = all(isinstance(A[j][k].re, int) and isinstance(A[j][k].im, int)
                   for j in range(n) for k in range(n))
    diag_pos = all(A[j][j].im == 0 and A[j][j].re > 0 for j in range(n))
    ck("A_is_hermitian_order_%d_over_Z[i]" % PAPER_ORDER,
       n == PAPER_ORDER and herm and integral, "order=%d" % n)
    ck("A_diagonal_strictly_positive", diag_pos,
       "diag=" + ",".join(str(A[j][j].re) for j in range(n)))
    return diag_pos


def det_exact(M):
    """Determinant by fraction-free Bareiss-style elimination over Cx."""
    n = len(M)
    if n == 0:
        return Cx(1, 0)
    W = [[Cx(Fraction(M[i][j].re), Fraction(M[i][j].im)) for j in range(n)]
         for i in range(n)]
    d = Cx(Fraction(1), Fraction(0))
    for c in range(n):
        p = None
        for r in range(c, n):
            if not W[r][c].is_zero():
                p = r
                break
        if p is None:
            return Cx(Fraction(0), Fraction(0))
        if p != c:
            W[c], W[p] = W[p], W[c]
            d = -d
        piv = W[c][c]
        d = d * piv
        for r in range(c + 1, n):
            if W[r][c].is_zero():
                continue
            f = W[r][c] / piv
            for j in range(c, n):
                W[r][j] = W[r][j] - f * W[c][j]
    return d


def rank_exact(M):
    """Exact rank by Gaussian elimination over Q(i)."""
    rows = [[Cx(Fraction(e.re), Fraction(e.im)) for e in row] for row in M]
    n = len(rows)
    m = len(rows[0]) if n else 0
    r = 0
    for c in range(m):
        p = None
        for i in range(r, n):
            if not rows[i][c].is_zero():
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        piv = rows[r][c]
        for i in range(r + 1, n):
            if rows[i][c].is_zero():
                continue
            f = rows[i][c] / piv
            for j in range(c, m):
                rows[i][j] = rows[i][j] - f * rows[r][j]
        r += 1
        if r == n:
            break
    return r


def psd_certificate(A):
    """Hermitian A is PSD iff every principal minor is >= 0.  When the rank is
    r, every principal minor of order > r vanishes, so it suffices to test
    orders 1 through r -- and ONLY through r: stopping short of the rank would
    certify matrices that are not PSD (Hermitian order 4 with diagonal 5 and
    every off-diagonal entry -2 has all principal minors of orders 1, 2 and 3
    positive yet determinant -343).  Returns (is_psd, rank, smallest 2x2
    minor)."""
    n = len(A)
    r = rank_exact(A)
    minors2 = [A[j][j].re * A[k][k].re - A[j][k].norm2()
               for j, k in combinations(range(n), 2)]
    ok1 = all(A[j][j].re >= 0 and A[j][j].im == 0 for j in range(n))
    higher_ok = True
    for order in range(3, r + 1):
        for t in combinations(range(n), order):
            d = det_exact(drop_rows_cols(A, set(range(n)) - set(t)))
            if d.im != 0 or d.re < 0:
                higher_ok = False
                break
        if not higher_ok:
            break
    return (ok1 and all(v >= 0 for v in minors2) and higher_ok,
            r, min(minors2))


def check_psd_and_rank(A):
    n = len(A)
    is_psd, r, worst = psd_certificate(A)
    ck("A_rank_is_exactly_%d" % PAPER_RANK, r == PAPER_RANK, "rank=%d" % r)
    ck("A_positive_semidefinite_exact_minor_certificate",
       is_psd and r == PAPER_RANK and worst >= 0,
       "min 2x2 principal minor=%d over %d minors" % (worst, n * (n - 1) // 2))
    higher = [det_exact(drop_rows_cols(A, set(range(n)) - set(t)))
              for t in combinations(range(n), 3)]
    ck("A_all_order_3_principal_minors_vanish",
       all(d.is_zero() for d in higher) and len(higher) == binom(n, 3),
       "tested=%d of the %d order-3 principal minors" % (len(higher),
                                                        binom(n, 3)))


def factorial(k):
    f = 1
    for t in range(2, k + 1):
        f *= t
    return f


def binom(n, k):
    """Exact binomial coefficient, used only to DERIVE the size a sweep is
    supposed to have so that the printed count is asserted, never assumed."""
    if k < 0 or k > n:
        return 0
    num, den = 1, 1
    for t in range(k):
        num *= n - t
        den *= t + 1
    return num // den


def multiset_count(pool, size):
    """Number of size-multisets from a pool of the given size, i.e. the exact
    number of iterations the combinations(range(pool + size - 1), size) trick
    below must perform."""
    return binom(pool + size - 1, size)


def per_rank2(cols):
    """Permanent of V^*V for a 2-row V, via the paper's coefficient identity
    per G = sum_r r!(m-r)! |c_r|^2 where prod_j (x_j + y_j z) = sum_r c_r z^r."""
    c = [Cx(1, 0)]
    for (x, y) in cols:
        nc = [Cx(0, 0)] * (len(c) + 1)
        for r, cr in enumerate(c):
            nc[r] = nc[r] + cr * x
            nc[r + 1] = nc[r + 1] + cr * y
        c = nc
    m = len(cols)
    return sum(factorial(r) * factorial(m - r) * c[r].norm2()
               for r in range(m + 1))


def check_permanent(A, cols):
    p = permanent(A)
    ck("permanent_A_is_real", p.im == 0, "Im=%s" % p.im)
    ck("permanent_A_matches_paper_value", p.re == PAPER_PER,
       "ryser=%d paper=%d" % (p.re, PAPER_PER))
    q = per_rank2(cols)
    ck("permanent_A_independent_rank_two_route_agrees", q == p.re,
       "coefficient identity=%d" % q)
    return p.re


def check_minor_table(A, cols):
    n = len(A)
    bad_w, bad_m = [], []
    computed = {}
    compared = 0
    for k in range(2, n + 1):
        w = A[0][k - 1].norm2()
        sub = drop_rows_cols(A, {0, k - 1})
        pm = permanent(sub)
        if pm.im != 0:
            bad_m.append(("im", k))
        computed[k] = (w, pm.re)
        if k in PAPER_TABLE:
            compared += 1
            pw, pmv = PAPER_TABLE[k]
            if w != pw:
                bad_w.append((k, w, pw))
            if pm.re != pmv:
                bad_m.append((k, pm.re, pmv))
        elif k == PAPER_ZERO_WEIGHT_INDEX:
            compared += 1
            if w != 0:
                bad_w.append((k, w, 0))
    # Coverage is asserted, not assumed: every k = 2..n must have been matched
    # against a paper entry, so a missing table row cannot pass as "no
    # mismatches".
    full = (compared == n - 1)
    ck("weights_abs_a_1k_squared_match_table", not bad_w and full,
       "compared %d of %d indices k=2..%d; mismatches=%s"
       % (compared, n - 1, n, bad_w if bad_w else "none"))
    ck("all_eleven_minor_permanents_match_table", not bad_m and full,
       "compared %d of %d indices k=2..%d; mismatches=%s"
       % (compared, n - 1, n, bad_m if bad_m else "none"))
    groups = [(4, 8), (5, 9), (6, 11), (7, 10)]
    eq = all(computed[a][1] == computed[b][1] for a, b in groups)
    dist = len(set(computed[k][1] for k in range(2, n)))
    # The expected number of distinct values is derived from the grouping the
    # paper prints, not asserted: (n - 2) indices k = 2..n-1, minus one for
    # each declared coincident pair.
    flat = [i for pair in groups for i in pair]
    expect_dist = (n - 2) - len(groups)
    ck("grouped_indices_share_minor_permanent",
       eq and dist == expect_dist and len(set(flat)) == len(flat)
       and all(2 <= i < n for i in flat),
       "pairs=%s distinct values among k=2..%d: %d (expected %d = %d indices "
       "- %d pairs)" % (groups, n - 1, dist, expect_dist, n - 2, len(groups)))
    agree = [k for k in range(2, n + 1)
             if per_rank2([cols[j] for j in range(n)
                           if j not in (0, k - 1)]) == computed[k][1]]
    ck("minor_permanents_reproduced_by_rank_two_identity",
       len(agree) == n - 1,
       "%d of %d minors (k=2..%d) agree by two independent routes"
       % (len(agree), n - 1, n))
    return computed


def check_violation(per_value, A):
    phi = phi1(A)
    ck("Phi1_A_matches_paper_value", phi == PAPER_PHI1,
       "computed=%d paper=%d" % (phi, PAPER_PHI1))
    gap = phi - per_value
    ck("inequality_is_REVERSED_Phi1_greater_than_per", gap > 0,
       "Phi1-per=%d" % gap)
    ck("gap_matches_paper_value", gap == PAPER_GAP,
       "computed=%d paper=%d" % (gap, PAPER_GAP))
    print("  NOT RE-RUN: what is refuted here is the inequality "
          "Phi_1(A) <= per(A) exactly as printed in this paper, with Phi_1 "
          "defined by the paper's own displayed formula "
          "Phi_1(A) = sum_{k=2}^n |a_{1k}|^2 per A[hat{1,k}].  The paper's "
          "ATTRIBUTION of that inequality -- that it is Zhang's Eq. (21) at "
          "p. 315, that it is the weak consequence of Eq. (20) which survived "
          "Hutchinson's disproof of the stronger conjecture, and that under "
          "Hutchinson's indexing B_{k-1,k-1} = A[hat{1,k}] -- is transcribed "
          "from the cited articles and is NOT checked here: this program has "
          "no network access and no copy of those articles, so it cannot "
          "confirm the equation numbers, the page, or the index convention.  "
          "If that identification were wrong, the arithmetic below would "
          "still hold but it would refute a differently-attributed statement.")
    return gap


class LCG(object):
    """Deterministic generator, so the census below is reproducible."""

    def __init__(self, seed=20240917):
        self.s = seed & 0xFFFFFFFFFFFF

    def nxt(self, k):
        self.s = (25214903917 * self.s + 11) & 0xFFFFFFFFFFFF
        return (self.s >> 16) % k

    def small(self, lo, hi):
        return lo + self.nxt(hi - lo + 1)


def permanent_by_definition(M):
    """per M = sum over all permutations, straight from the definition.
    Used only to cross-check the Ryser routine on small orders."""
    m = len(M)
    tot = Cx(0, 0)
    for s in permutations(range(m)):
        pr = Cx(1, 0)
        for j in range(m):
            pr = pr * M[j][s[j]]
        tot = tot + pr
    return tot


def check_ryser_against_definition():
    """The permanent routine used everywhere else is Ryser's inclusion-
    exclusion; compare it with the defining permutation sum, and confirm the
    empty-matrix convention per() = 1 that Phi_1 needs at order 2."""
    g = LCG(101)
    bad = []
    orders = range(1, 8)
    reps = 4
    tested = 0
    for m in orders:
        for _ in range(reps):
            M = [[Cx(g.small(-4, 4), g.small(-4, 4)) for _ in range(m)]
                 for _ in range(m)]
            tested += 1
            if permanent(M) != permanent_by_definition(M):
                bad.append(m)
    A = gram(build_columns())
    A2 = [[A[0][0], A[0][1]], [A[1][0], A[1][1]]]
    empty_ok = (permanent([]) == ONE and phi1(A2) == A[0][1].norm2())
    ck("ryser_permanent_matches_the_permutation_sum",
       not bad and empty_ok and tested == len(orders) * reps and tested > 0,
       "matrices compared=%d (orders %d..%d, %d each); empty permanent=%s, "
       "Phi1 of the 2x2 leading block=%d"
       % (tested, min(orders), max(orders), reps, permanent([]), phi1(A2)))


def check_lemma_against_ryser():
    """Lemma: per(V^*V) = sum_r r!(m-r)!|c_r|^2.  Test it against Ryser on
    rank<=2 Gram matrices generated independently of the paper's V."""
    g = LCG()
    bad = []
    tested = 0
    orders = range(1, 8)
    reps = 6
    for m in orders:
        for _ in range(reps):
            cols = [(Cx(g.small(-3, 3), g.small(-3, 3)),
                     Cx(g.small(-3, 3), g.small(-3, 3))) for _ in range(m)]
            A = gram(cols)
            p = permanent(A)
            q = per_rank2(cols)
            tested += 1
            if p.im != 0 or p.re != q:
                bad.append((m, p.re, q))
    ck("rank_two_permanent_identity_agrees_with_ryser",
       not bad and tested == len(orders) * reps and tested > 0,
       "cases=%d of %d (orders %d..%d, %d each) mismatches=%s"
       % (tested, len(orders) * reps, min(orders), max(orders), reps,
          bad if bad else "none"))


def check_rank_one_corollary():
    """Rank <= 1: per A = n! P and Phi_1(A) = (n-1)! P = per(A)/n."""
    g = LCG(777)
    bad = []
    orders = range(2, 8)
    tested = 0
    for n in orders:
        tested += 1
        u = [Cx(g.small(-4, 4), g.small(-4, 4)) for _ in range(n)]
        if any(e.is_zero() for e in u):
            u = [e if not e.is_zero() else Cx(1, 1) for e in u]
        A = [[u[j].conj() * u[k] for k in range(n)] for j in range(n)]
        P = 1
        for e in u:
            P *= e.norm2()
        per = permanent(A)
        phi = phi1(A)
        if (per.re != factorial(n) * P or per.im != 0
                or phi != factorial(n - 1) * P or phi * n != per.re
                or not (phi < per.re)):
            bad.append((n, per.re, phi, P))
    zero_gaps = []
    for n in range(1, 6):
        Z = [[Cx(0, 0)] * n for _ in range(n)]
        zero_gaps.append(phi1(Z) - permanent(Z).re)
    ck("rank_one_gives_Phi1_equal_per_over_n",
       not bad and all(v == 0 for v in zero_gaps)
       and tested == len(orders) and tested > 0 and len(zero_gaps) > 0,
       "orders %d..%d tested=%d mismatches=%s; rank-0 gaps=%s"
       % (min(orders), max(orders), tested, bad if bad else "none", zero_gaps))


def check_order_two_and_three_identities():
    """Order 2: per - Phi_1 = a11 a22.  Order 3: per - Phi_1 =
    det A + 2a|e|^2 + d|c|^2 + f|b|^2.  Both sides are polynomials of degree
    at most 1 in each real diagonal entry and at most 2 in each real
    coordinate of each off-diagonal entry, so agreement on the product grid
    {0,1} x ... x {0,1,2} x ... below is a PROOF of the identity, not a
    sample: a nonzero polynomial cannot vanish on a grid exceeding its
    per-variable degree in every variable."""
    diag_vals = (0, 1)
    smalls = (0, 1, 2)
    bad2 = []
    pts2 = 0
    for a in diag_vals:
        for d in diag_vals:
            for br in smalls:
                for bi in smalls:
                    b = Cx(br, bi)
                    A = [[Cx(a, 0), b], [b.conj(), Cx(d, 0)]]
                    pts2 += 1
                    if permanent(A).re - phi1(A) != a * d:
                        bad2.append((a, d, br, bi))
    want2 = len(diag_vals) ** 2 * len(smalls) ** 2
    ck("order_2_gap_identity_equals_a11_a22_on_a_proving_grid",
       not bad2 and pts2 == want2 and len(diag_vals) > 1 and len(smalls) > 2,
       "grid points=%d of %d, diagonal in %s (degree 1), off-diagonal "
       "coordinates in %s (degree 2), mismatches=%s"
       % (pts2, want2, str(diag_vals), str(smalls),
          bad2 if bad2 else "none"))

    bad3 = []
    pts3 = 0
    for a in diag_vals:
        for d in diag_vals:
            for f in diag_vals:
                for br in smalls:
                    for bi in smalls:
                        for cr in smalls:
                            for ci in smalls:
                                for er in smalls:
                                    for ei in smalls:
                                        b, c, e = (Cx(br, bi), Cx(cr, ci),
                                                   Cx(er, ei))
                                        A = [[Cx(a, 0), b, c],
                                             [b.conj(), Cx(d, 0), e],
                                             [c.conj(), e.conj(), Cx(f, 0)]]
                                        det = det_exact(A)
                                        rhs = (det.re + 2 * a * e.norm2()
                                               + d * c.norm2()
                                               + f * b.norm2())
                                        pts3 += 1
                                        if (det.im != 0 or
                                                permanent(A).re - phi1(A)
                                                != rhs):
                                            bad3.append((a, d, f, br, bi))
    want3 = len(diag_vals) ** 3 * len(smalls) ** 6
    ck("order_3_gap_identity_of_the_corollary_on_a_proving_grid",
       not bad3 and pts3 == want3 and len(diag_vals) > 1 and len(smalls) > 2,
       "grid points=%d of %d, diagonal in %s (degree 1 each), six off-diagonal "
       "coordinates in %s (degree 2 each), mismatches=%s"
       % (pts3, want3, str(diag_vals), str(smalls),
          bad3 if bad3 else "none"))
    return (not bad2 and not bad3 and pts2 == want2 and pts3 == want3)


def gram_general(vecs):
    m = len(vecs)
    dim = len(vecs[0])
    return [[sum((vecs[j][t].conj() * vecs[k][t] for t in range(dim)),
                 Cx(0, 0)) for k in range(m)] for j in range(m)]


def check_small_orders_exhaustively(ident_proved):
    """Orders 1, 2, 3 cannot host a counterexample: exhaustively over every
    Gram matrix of three vectors drawn from a named finite set (all ranks up
    to 3), and by the proved identities above."""
    base = [Cx(0, 0), Cx(1, 0), Cx(0, 1), Cx(1, 1)]
    vecs = [(p, q, r) for p in base for q in base for r in base]
    worst = None
    hits = 0
    tried = 0
    for combo in combinations(range(len(vecs) + 2), 3):
        vs = [vecs[combo[t] - t] for t in range(3)]
        A = gram_general(vs)
        g = phi1(A) - permanent(A).re
        tried += 1
        if worst is None or g > worst:
            worst = g
        if g > 0:
            hits += 1
    want = multiset_count(len(vecs), 3)
    ck("exhaustive_order_3_gram_family_has_no_violation",
       hits == 0 and tried == want and tried > 0,
       "Grams tested=%d of the %d 3-multisets from %d vectors, "
       "max(Phi1-per)=%s" % (tried, want, len(vecs), worst))
    order1 = [phi1([[Cx(t, 0)]]) - permanent([[Cx(t, 0)]]).re
              for t in range(0, 5)]
    # The bound is "one more than the largest order settled here", derived from
    # the orders actually covered (1, 2 and 3), not written down as a literal.
    orders_settled = ([1] if all(v <= 0 for v in order1) else []) + (
        [2, 3] if (ident_proved and hits == 0 and tried == want) else [])
    established_lower = (max(orders_settled) + 1
                         if orders_settled == [1, 2, 3] else None)
    ck("least_order_lower_bound_matches_the_paper",
       established_lower == PAPER_NMIN_LOWER,
       "established n_min >= %s, paper states %d; order-1 gaps=%s"
       % (established_lower, PAPER_NMIN_LOWER, order1))


def scale_by_diagonal(M, ds):
    n = len(M)
    return [[M[j][k] * (ds[j] * ds[k]) for k in range(n)] for j in range(n)]


def prod_of_squares(ds):
    f = 1
    for d in ds:
        f *= d * d
    return f


def check_scaling_homogeneity_on_a_proving_grid():
    """The homogeneity per(DAD) = (prod d_j^2) per(A) and the same for Phi_1
    are POLYNOMIAL IDENTITIES, and at orders 2 and 3 they are proved here, not
    sampled.  On the Hermitian family with diagonal entries a_jj and real and
    imaginary parts of the off-diagonal entries as free parameters, both sides
    have degree at most 1 in each a_jj and at most 2 in each off-diagonal
    coordinate (a_jk occurs together with its conjugate), and degree at most 2
    in each d_j.  So agreement on the product grid {0,1} for the diagonal,
    {0,1,2} for each off-diagonal coordinate and {1,2,3} for each d_j forces
    the identity: a nonzero polynomial cannot vanish on a grid that exceeds
    its degree in every variable."""
    dgrid = (1, 2, 3)
    for n, params in ((2, "a11,a22 and Re b,Im b"),
                      (3, "a11,a22,a33 and Re,Im of b,c,e")):
        diag_vals = (0, 1)
        off_vals = (0, 1, 2)
        offs = [(j, k) for j in range(n) for k in range(j + 1, n)]
        bad_per = []
        bad_phi = []
        pts = 0
        diags = [[]]
        for _ in range(n):
            diags = [d + [v] for d in diags for v in diag_vals]
        offsets = [[]]
        for _ in offs:
            offsets = [o + [Cx(r, i)] for o in offsets
                       for r in off_vals for i in off_vals]
        dsets = [[]]
        for _ in range(n):
            dsets = [d + [v] for d in dsets for v in dgrid]
        for dg in diags:
            for ov in offsets:
                M = [[Cx(0, 0)] * n for _ in range(n)]
                for j in range(n):
                    M[j][j] = Cx(dg[j], 0)
                for t, (j, k) in enumerate(offs):
                    M[j][k] = ov[t]
                    M[k][j] = ov[t].conj()
                p0 = permanent(M).re
                f0 = phi1(M)
                for ds in dsets:
                    S = scale_by_diagonal(M, ds)
                    f = prod_of_squares(ds)
                    pts += 1
                    if permanent(S).re != f * p0:
                        bad_per.append((dg, ov, ds))
                    if phi1(S) != f * f0:
                        bad_phi.append((dg, ov, ds))
        expected = len(diags) * len(offsets) * len(dsets)
        covered = (pts == expected and len(diag_vals) > 1
                   and len(off_vals) > 2 and len(dgrid) > 2)
        detail = ("grid points=%d of %d (%d diagonals x %d off-diagonal tuples "
                  "x %d scalings) over %s, diagonal in %s, off-diagonal "
                  "coordinates in %s, d in %s; mismatches=%s"
                  % (pts, expected, len(diags), len(offsets), len(dsets),
                     params, str(diag_vals), str(off_vals), str(dgrid), "%s"))
        ck("scaling_homogeneity_of_per_order_%d_on_a_proving_grid" % n,
           not bad_per and covered,
           detail % (len(bad_per) if bad_per else "none"))
        ck("scaling_homogeneity_of_Phi1_order_%d_on_a_proving_grid" % n,
           not bad_phi and covered,
           detail % (len(bad_phi) if bad_phi else "none"))


def check_diagonal_scaling(A, gap):
    """per(DAD) and Phi_1(DAD) are both homogeneous of degree 2 in each d_j,
    so a positive diagonal scaling cannot change the sign of the gap.  This is
    the step that turns the exhibited A into a correlation matrix."""
    n = len(A)
    for name, ds in (("integers", [1, 2, 1, 3, 1, 2, 1, 1, 2, 1, 1, 3]),
                     ("rationals", [Fraction(1, 2), 1, Fraction(3, 2), 1, 1,
                                    Fraction(1, 3), 2, 1, 1, 1, 1, 3])):
        S = [[A[j][k] * (ds[j] * ds[k]) for k in range(n)] for j in range(n)]
        f = 1
        for d in ds:
            f *= d * d
        p = permanent(S)
        ph = phi1(S)
        ok = (p.re == f * PAPER_PER and ph == f * PAPER_PHI1
              and ph - p.re == f * gap and ph - p.re > 0)
        ck("diagonal_scaling_%s_preserves_the_violation" % name, ok,
           "scale factor prod d_j^2=%s, scaled gap=%s" % (f, ph - p.re))
    check_scaling_homogeneity_coordinatewise(A)
    check_correlation_normalization_is_irrational(A)


def check_scaling_homogeneity_coordinatewise(A):
    """The two scalings above are only two points.  Here the homogeneity is
    forced for the ACTUAL order-12 A, one coordinate at a time: with the other
    eleven entries of d held fixed, per(DAD) and Phi_1(DAD) are polynomials of
    degree at most 2 in the remaining d_j, so agreement with d_j^2 times the
    unscaled value at THREE distinct d_j is an identity in d_j, not a sample.
    Chaining the twelve coordinates walks from d = (1,...,1) to each target
    diagonal, every step proved on its own three-point grid."""
    n = len(A)
    targets = (("integers", [1, 2, 1, 3, 1, 2, 1, 1, 2, 1, 1, 3]),
               ("rationals", [Fraction(1, 2), 1, Fraction(3, 2), 1, 1,
                              Fraction(1, 3), 2, 1, 1, 1, 1, 3]))
    for name, tgt in targets:
        cur = [1] * n
        bad = []
        pts = 0
        fewest_distinct = None
        for j in range(n):
            vals = [tgt[j]]
            for v in (1, 2, 3, 4):
                if len(vals) < 3 and v not in vals:
                    vals.append(v)
            nd = len(set(vals))
            if fewest_distinct is None or nd < fewest_distinct:
                fewest_distinct = nd
            base = scale_by_diagonal(A, cur)
            p0 = permanent(base).re
            f0 = phi1(base)
            for v in vals:
                ds = list(cur)
                ds[j] = v
                S = scale_by_diagonal(A, ds)
                pts += 1
                if (permanent(S).re != v * v * p0 or phi1(S) != v * v * f0):
                    bad.append((j + 1, v))
            cur[j] = tgt[j]
        ck("scaling_homogeneity_order_12_coordinatewise_to_the_%s" % name,
           not bad and cur == list(tgt) and pts == 3 * n
           and fewest_distinct is not None and fewest_distinct >= 3,
           "%d coordinates, %d scaled evaluations, fewest distinct d_j on any "
           "one coordinate=%s (>= 3 is what forces a degree-2 polynomial), "
           "mismatches=%s" % (n, pts, fewest_distinct,
                              bad if bad else "none"))
    print("  NOT RE-RUN: the twelve coordinatewise steps just printed force "
          "the homogeneity in ONE d_j at a time, with the other eleven held at "
          "the values reached so far, so they prove the identity along the two "
          "chains from d=(1,...,1) to the two target diagonals -- they do NOT "
          "prove the full twelve-variable identity at order 12, which would "
          "need a grid of 3^12 = %d points and is not run here.  The full "
          "identity is proved only at orders 2 and 3 above." % (3 ** 12))


def check_correlation_normalization_is_irrational(A):
    """Corollary 4 states the counterexample as a positive definite
    CORRELATION matrix, i.e. after the normalization d_j = a_jj^(-1/2).  That
    scaling vector is checked here to be irrational, which is why no
    unit-diagonal matrix is exhibited anywhere in this program: the exact
    rational engine used throughout cannot represent it.  The correlation form
    therefore rests on the homogeneity proved above (which holds over the
    reals), not on a recomputation here."""
    n = len(A)
    diag = [A[j][j].re for j in range(n)]
    nonsquares = []
    for j, a in enumerate(diag):
        r = 0
        while r * r < a:
            r += 1
        if r * r != a:
            nonsquares.append((j + 1, a))
    ck("correlation_normalization_needs_an_irrational_scaling",
       bool(nonsquares) and all(a > 0 for a in diag),
       "diag=%s; a_jj not a perfect square, so d_j=a_jj^(-1/2) is irrational, "
       "at j=%s" % (diag, [j for j, _ in nonsquares]))
    print("  NOT RE-RUN: no unit-diagonal (correlation) matrix is produced "
          "here.  The normalization d_j = a_jj^(-1/2) is irrational for the "
          "diagonal just printed, so the correlation form of Corollary 4 "
          "cannot be exhibited in this program's exact rational arithmetic; "
          "it follows from the homogeneity identity over the reals, which is "
          "proved above on a full grid at orders 2 and 3 and, at order 12, "
          "only along the two rational chains just printed -- never at the "
          "irrational d of this normalization -- and so the correlation form "
          "rests on the paper's argument, not on any matrix computed here.")


def direct_sum_with_identity(A, m):
    n = len(A)
    N = n + m
    C = [[Cx(0, 0)] * N for _ in range(N)]
    for j in range(n):
        for k in range(n):
            C[j][k] = A[j][k]
    for j in range(n, N):
        C[j][j] = Cx(1, 0)
    return C


DIRECT_SUM_PADDINGS = (1, 2, 3, 4)


def check_direct_sum(gap):
    """C = A (+) I_m has Phi_1(C) - per(C) = Phi_1(A) - per(A), the step
    behind 'counterexamples exist in every order at least 12'."""
    A = gram(build_columns())
    n = len(A)
    for m in DIRECT_SUM_PADDINGS:
        N = n + m
        C = direct_sum_with_identity(A, m)
        p = permanent(C)
        ph = phi1(C)
        rk = rank_exact(C)
        ok = (p.re == PAPER_PER and ph == PAPER_PHI1 and ph - p.re == gap
              and ph - p.re > 0 and rk == PAPER_RANK + m and len(C) == N)
        ck("direct_sum_with_I_%d_order_%d_keeps_the_same_gap" % (m, N), ok,
           "order=%d per=%d Phi1-per=%d rank=%d (paper's 2+m=%d)"
           % (len(C), p.re, ph - p.re, rk, PAPER_RANK + m))
    orders = [n] + [n + m for m in DIRECT_SUM_PADDINGS]
    print("  NOT RE-RUN: the abstract and Corollary 4 assert counterexamples "
          "in EVERY order N >= %d.  Only the finitely many orders %s are "
          "recomputed here (m = 0 and m = %s in C = A (+) I_m); every larger N "
          "rests on the induction in the paper, which is not machine-checked, "
          "and the paper's rank claim rank(C) = 2 + m is checked only at those "
          "same orders.  The sums formed here pad the UNNORMALIZED A, not the "
          "correlation matrix A_0 of the corollary's proof; A_0 (+) I_m equals "
          "D'(A (+) I_m)D' for the positive diagonal D' that normalizes A and "
          "fixes the identity block, so the two have the same gap sign by the "
          "homogeneity identity - which is proved above at orders 2 and 3 and "
          "coordinatewise at order 12, and is NOT re-verified at orders %s."
          % (n, ",".join(str(o) for o in orders),
             ",".join(str(m) for m in DIRECT_SUM_PADDINGS),
             ",".join(str(n + m) for m in DIRECT_SUM_PADDINGS)))


def leading_minors_positive(M):
    n = len(M)
    for r in range(1, n + 1):
        d = det_exact([row[:r] for row in M[:r]])
        if d.im != 0 or d.re <= 0:
            return False, r
    return True, n


def check_positive_definite_perturbation():
    """For A + eps I with eps = 1/T, scaling by T gives the integer matrix
    T*A + I; positive definiteness and the sign of the gap are unchanged.
    Find a T for which the perturbed matrix is positive definite and still
    violates the inequality."""
    A = gram(build_columns())
    n = len(A)
    good = []
    denominators = (1, 10, 100, 10 ** 4, 10 ** 6)
    seen = 0
    for T in denominators:
        seen += 1
        B = [[A[j][k] * T for k in range(n)] for j in range(n)]
        for j in range(n):
            B[j][j] = B[j][j] + ONE
        pd, _ = leading_minors_positive(B)
        p = permanent(B)
        ph = phi1(B)
        if pd and ph - p.re > 0:
            good.append(T)
        print("  eps=1/%d: positive definite=%s, sign(Phi1-per)=%d"
              % (T, pd, (ph - p.re > 0) - (ph - p.re < 0)))
    ck("positive_definite_perturbation_still_violates",
       bool(good) and seen == len(denominators),
       "eps=1/T works for %d of the %d T tried, T in %s"
       % (len(good), seen, good))
    check_positive_definite_perturbation_of_the_direct_sum(
        A, good[0] if good else 10 ** 6)


def check_positive_definite_perturbation_of_the_direct_sum(A, T):
    """Corollary 4 perturbs the DIRECT SUM C = A (+) I_{N-12}, not A itself,
    so the perturbation is repeated on C at EVERY order the direct sum is
    recomputed at (all of DIRECT_SUM_PADDINGS, not a prefix of it).  Scaling
    A + eps I by T = 1/eps again gives an integer matrix T*C + I whose positive
    definiteness and gap sign are the ones asserted."""
    n = len(A)
    ok_all = []
    for m in DIRECT_SUM_PADDINGS:
        C = direct_sum_with_identity(A, m)
        N = len(C)
        B = [[C[j][k] * T for k in range(N)] for j in range(N)]
        for j in range(N):
            B[j][j] = B[j][j] + ONE
        pd, _ = leading_minors_positive(B)
        p = permanent(B)
        ph = phi1(B)
        ok_all.append((N, pd, ph - p.re > 0))
        print("  order %d, eps=1/%d: positive definite=%s, sign(Phi1-per)=%d"
              % (N, T, pd, (ph - p.re > 0) - (ph - p.re < 0)))
    ck("positive_definite_perturbation_of_the_direct_sum_still_violates",
       len(ok_all) == len(DIRECT_SUM_PADDINGS)
       and all(pd and pos for _, pd, pos in ok_all),
       "orders %s (%d of %d paddings) at eps=1/%d, base order %d"
       % ([N for N, _, _ in ok_all], len(ok_all),
          len(DIRECT_SUM_PADDINGS), T, n))


def check_order_and_rank_bounds(A, gap):
    """The exhibited object supplies the paper's upper bounds: rank 2 is
    attained, and order 12 is attained."""
    _, r, _ = psd_certificate(A)
    ck("least_rank_upper_bound_two_is_attained",
       gap > 0 and r == PAPER_RANK, "violating matrix has rank %d" % r)
    ck("least_order_upper_bound_matches_the_paper",
       gap > 0 and len(A) == PAPER_NMIN_UPPER,
       "violating matrix has order %d, paper states n_min <= %d"
       % (len(A), PAPER_NMIN_UPPER))


def check_reported_comparison_of_an_external_matrix():
    """The paper also reports, for a real matrix printed in another article,
    Phi_1 = 38673/25000 < 2879693/800000 = per.  That matrix is not
    reproduced here, so only the printed ordering is checked."""
    ln, ld = 38673, 25000
    rn, rd = 2879693, 800000
    ck("reported_external_values_are_ordered_as_printed", ln * rd < rn * ld,
       "cross products %d < %d" % (ln * rd, rn * ld))
    print("  NOT RE-RUN: the matrix behind those two rationals is printed in "
          "a different article and is not restated in this paper, so its "
          "entries cannot be recomputed from the text verified here.")


def check_order_n_minus_one_variant_is_nonhomogeneous(A):
    """The paper asserts that the variant using order-(n-1) principal minors
    A[hat k] in place of A[hat{1,k}] is NONHOMOGENEOUS, hence a different
    quantity from Phi_1.  That assertion is checkable here: under a positive
    diagonal scaling D, per(DAD) picks up prod_j d_j^2, whereas the variant
    picks up d_1^2 * prod_j d_j^2, because the weight |a_{1k}|^2 contributes
    d_1^2 d_k^2 while the surviving minor supplies prod_{j != k} d_j^2."""
    n = len(A)

    def variant(M):
        tot = 0
        for k in range(1, n):
            w = M[0][k].norm2()
            if w == 0:
                continue
            p = permanent(drop_rows_cols(M, {k}))
            tot += w * p.re
        return tot

    ds = [2] + [1] * (n - 1)
    S = [[A[j][k] * (ds[j] * ds[k]) for k in range(n)] for j in range(n)]
    f = 1
    for d in ds:
        f *= d * d
    v0, v1 = variant(A), variant(S)
    p0, p1 = permanent(A).re, permanent(S).re
    scales_right = (v1 == ds[0] * ds[0] * f * v0 and p1 == f * p0)
    differs = (ds[0] * ds[0] != 1 and v1 * p0 != v0 * p1)
    ck("order_n_minus_1_variant_scales_differently_from_per",
       scales_right and differs and v0 != phi1(A),
       "variant factor=%d, per factor=%d, variant(A)=%d != Phi1(A)"
       % (v1 // v0, p1 // p0, v0))


def gap_rank2(cols):
    """Phi_1(V^*V) - per(V^*V) for a 2-row V, using only the coefficient
    identity, so a whole family can be swept cheaply."""
    m = len(cols)
    x0, y0 = cols[0]
    total = 0
    for k in range(1, m):
        a1k = x0.conj() * cols[k][0] + y0.conj() * cols[k][1]
        w = a1k.norm2()
        if w:
            total += w * per_rank2([cols[j] for j in range(m) if j not in (0, k)])
    return total - per_rank2(cols)


def census_column_set(xmax, ys):
    return [(Cx(x, 0), Cx(yr, yi)) for x in range(xmax + 1) for (yr, yi) in ys]


def check_census():
    """Positive control first: the same sweep evaluator must flag the paper's
    own 12-column configuration.  Then sweep smaller orders."""
    paper_gap = gap_rank2(build_columns())
    ck("sweep_evaluator_reproduces_the_paper_gap", paper_gap == PAPER_GAP,
       "sweep gap=%d" % paper_gap)

    ys4 = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (1, 1), (1, -1),
           (2, 1), (2, -1), (1, 2), (3, 2)]
    S4 = census_column_set(3, ys4)
    hits, tried, worst = [], 0, None
    for combo in combinations(range(len(S4) + 3), 4):
        cols = [S4[combo[t] - t] for t in range(4)]  # multisets of size 4
        g = gap_rank2(cols)
        tried += 1
        if worst is None or g > worst:
            worst = g
        if g > 0:
            hits.append(cols)
    want4 = multiset_count(len(S4), 4)
    ck("exhaustive_order_4_family_has_no_violation",
       not hits and tried == want4 and tried > 0,
       "multisets tested=%d of %d, columns in set=%d max(Phi1-per)=%s"
       % (tried, want4, len(S4), worst))

    ys5 = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (1, -1), (2, 3)]
    S5 = census_column_set(2, ys5)
    hits5, tried5, worst5 = 0, 0, None
    for combo in combinations(range(len(S5) + 4), 5):
        cols = [S5[combo[t] - t] for t in range(5)]
        g = gap_rank2(cols)
        tried5 += 1
        if worst5 is None or g > worst5:
            worst5 = g
        if g > 0:
            hits5 += 1
    want5 = multiset_count(len(S5), 5)
    ck("exhaustive_order_5_family_has_no_violation",
       hits5 == 0 and tried5 == want5 and tried5 > 0,
       "multisets tested=%d of %d, columns in set=%d max(Phi1-per)=%s"
       % (tried5, want5, len(S5), worst5))

    rng = LCG(90210)
    found, total = [], 0
    sweep_orders = range(4, 12)
    per_order = 900
    for m in sweep_orders:
        for t in range(per_order):
            if t % 2:
                cols = [(ONE, ZERO)] + [
                    (Cx(rng.small(0, 9), 0),
                     Cx(rng.small(-9, 9), rng.small(-9, 9)))
                    for _ in range(m - 2)] + [(ZERO, ONE)]
            else:
                cols = [(Cx(rng.small(-6, 6), rng.small(-6, 6)),
                         Cx(rng.small(-6, 6), rng.small(-6, 6)))
                        for _ in range(m)]
            total += 1
            if gap_rank2(cols) > 0:
                found.append(m)
    want_sweep = len(sweep_orders) * per_order
    ck("random_sweep_orders_4_to_11_finds_no_violation",
       not found and total == want_sweep and total > 0,
       "configurations tested=%d of %d (orders %d..%d, %d each) violations=%s"
       % (total, want_sweep, min(sweep_orders), max(sweep_orders), per_order,
          sorted(set(found))))
    print("  NOT RE-RUN: no exhaustive census over all Hermitian positive "
          "semidefinite matrices of orders 4..11 exists here, so the paper's "
          "lower bound n_min >= 4 is verified only for orders 2 and 3 (by the "
          "exact identities above) and the band 4 <= n_min <= 12 is left as "
          "stated; ranks above 2 are not swept.")


def check_negative_controls():
    """The certificates above must be able to say no.  Each control below is
    a deliberately damaged input that they are required to reject."""
    A = gram(build_columns())
    n = len(A)
    B = [row[:] for row in A]
    B[1][1] = B[1][1] - Cx(5, 0)          # diagonal entry 5 -> 0, still Hermitian
    is_psd, _, worst = psd_certificate(B)
    ck("control_non_psd_hermitian_matrix_is_rejected",
       (not is_psd) and worst < 0, "smallest 2x2 principal minor=%d" % worst)

    insensitive = []
    edit_at = (2, 4, 6, 9)
    edited = []
    for j in edit_at:
        cols = build_columns()
        cols[j] = (cols[j][0], cols[j][1] + ONE)
        M = gram(cols)
        p = permanent(M).re
        ph = phi1(M)
        edited.append(j + 1)
        if p == PAPER_PER or ph == PAPER_PHI1 or ph - p == PAPER_GAP:
            insensitive.append(j)
    ck("control_single_column_edits_change_both_sides",
       not insensitive and len(edited) == len(edit_at) and edited,
       "edited columns %s (1-based); values unchanged for=%s"
       % (",".join(str(c) for c in edited),
          insensitive if insensitive else "none"))

    ident_orders = (4, 8, 12)
    holds = []
    for m in ident_orders:
        I = [[ONE if j == k else ZERO for k in range(m)] for j in range(m)]
        holds.append((permanent(I).re, phi1(I)))
    ck("control_identity_matrices_satisfy_the_inequality",
       len(holds) == len(ident_orders) and bool(holds)
       and all(ph < p for p, ph in holds),
       "orders %s: (per,Phi1)=%s"
       % (",".join(str(m) for m in ident_orders), holds))

    g = LCG(5150)
    rank3_orders = (4, 5, 6)
    rank3 = []
    for m in rank3_orders:
        cols3 = [[Cx(g.small(-3, 3), g.small(-3, 3)) for _ in range(3)]
                 for _ in range(m)]
        M = [[sum((cols3[j][t].conj() * cols3[k][t] for t in range(3)),
                  Cx(0, 0)) for k in range(m)] for j in range(m)]
        p = permanent(M)
        rank3.append((m, rank_exact(M), p.im == 0, p.re - phi1(M)))
    ck("control_rank_three_grams_are_measured_not_assumed",
       len(rank3) == len(rank3_orders) and bool(rank3)
       and all(r >= 3 and real for _, r, real, _ in rank3),
       "order,rank,per real,per-Phi1=%s" % rank3)
    print("  the last control only records the gap sign on rank-3 samples; "
          "the inequality is not claimed to be true for them.")


def main():
    cols = build_columns()
    print("== the exhibited object ==")
    check_exhibited_object(cols)
    A = gram(cols)
    check_gram_shape(A)
    for j, row in enumerate(A):
        print("  A row %2d: " % (j + 1) + " ".join(repr(e) for e in row))
    check_psd_and_rank(A)
    print("== the two sides of the inequality ==")
    per_value = check_permanent(A, cols)
    check_minor_table(A, cols)
    gap = check_violation(per_value, A)
    print("== the identities the proofs rest on ==")
    check_ryser_against_definition()
    check_lemma_against_ryser()
    check_rank_one_corollary()
    ident = check_order_two_and_three_identities()
    print("== the consequences ==")
    check_small_orders_exhaustively(ident)
    check_order_and_rank_bounds(A, gap)
    check_reported_comparison_of_an_external_matrix()
    check_order_n_minus_one_variant_is_nonhomogeneous(A)
    check_scaling_homogeneity_on_a_proving_grid()
    check_diagonal_scaling(A, gap)
    check_direct_sum(gap)
    check_positive_definite_perturbation()
    print("== negative controls ==")
    check_negative_controls()
    print("== how far down the order the sweep reaches ==")
    check_census()
    ck("no_complex_minor_permanent_seen_anywhere", not ANOMALIES,
       "hermitian inputs must give real minor permanents; anomalies=%d"
       % len(ANOMALIES))
    return finish()


if __name__ == "__main__":
    sys.exit(main())

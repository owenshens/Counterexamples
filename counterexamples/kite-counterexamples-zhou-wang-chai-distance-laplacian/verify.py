#!/usr/bin/env python3
"""Verification of kite-graph counterexamples to the Zhou-Wang-Chai
distance-Laplacian eigenvalue-sum bound.

TRANSCRIBED AND NOT TESTABLE HERE (this program has no network access, so
the wording of the conjecture cannot be collated with its published source):
  * The bound under attack: U_r(G) <= W(G) + C(r+2,3) for 2 <= r <= n,
    where U_r(G) is the sum of the r largest distance Laplacian
    eigenvalues and W(G) is the Wiener index.  Every check below tests
    this inequality exactly AS QUOTED.  If the original attached a
    hypothesis that was lost in transcription -- a diameter, girth or
    r-range restriction -- the checks would be unaffected but the target
    would be a different statement.  The one piece of offline
    corroboration is that the exhaustive diameter-two census below
    reproduces precisely the classification quoted from earlier work for
    this same inequality, namely K_{1,3} at r=2 and K_{1,4} at r=3.
    This is disclosed again in the closing NOT RE-RUN lines.

TAKEN FROM THE PAPER (inputs, quoted and then tested, never assumed):
  * The exhibited family: H_n is the kite obtained from the path
    v_1 v_2 ... v_n by adding the edge v_1 v_3, i.e. a triangle with a
    pendant path.
  * The claimed violating indices: r = n-2 for every n >= 5, and r = 2m
    along the subfamily n = 3m, m >= 4.
  * The claimed excess lower bounds (3n-14)/2 and
    (5m^3-12m^2-17m+12)/6, the asymptotic constant 5/162, the value
    p(4) = 72 and the recurrence p(m+1)-p(m) = 3(m+1)(5m-8).
  * The closed forms W(H_n) = C(n+1,3)-(n-2), tau(v_3) = (n^2-5n+10)/2,
    tau(v_4) = (n^2-7n+22)/2, the Rayleigh value (n^2-6n+18)/2, and the
    selected-transmission sums (20m^3-2m)/3 and (20m^3-17m+12)/3.
  * The displayed matrix D^L(H_5), the factorisation
    x(x-6)(x-8)(x^2-20x+95), the spectrum {10+-sqrt5, 8, 6, 0},
    W(H_5) = 17 and U_3(H_5) = 28 > 27.
  * The claim that H_5 has minimum order among unicyclic
    counterexamples, the spectra listed for C_3, C_4 and H_4, the
    remark's path bounds (3n-10)/2 and the m >= 3 threshold, and the two
    counterexamples quoted from earlier work, K_{1,3} at r=2 and
    K_{1,4} at r=3.

DERIVED HERE (computed from scratch, with nothing but the definitions):
  * H_n and every quantity attached to it: BFS distances, transmissions,
    Wiener index, distance Laplacian, and the characteristic polynomial
    of that matrix by Faddeev-LeVerrier over the integers.
  * Eigenvalue localisation by exact integer root counting (Descartes'
    rule, which is an equality for the real-rooted polynomials arising
    from symmetric matrices), giving certified rational brackets.  No
    floating-point number decides any check; floats appear only inside
    printed detail strings.
  * The excess U_r - W - C(r+2,3) at the claimed indices, computed and
    compared with the paper's lower bounds.
  * Exhaustive censuses of connected graphs of small order, including a
    proof of equality in the borderline cases via a norm-based
    separation bound for algebraic integers.
"""

import sys
from fractions import Fraction
from itertools import combinations, permutations

CHECKS = []

# Filled in by the exhaustive census so that the closing disclosure prints
# derived counts rather than transcribed ones.
CENSUS_STATS = {}


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (" [" + detail + "]" if detail else ""))
    return bool(ok)


def finish():
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


def binom(a, b):
    if b < 0 or b > a:
        return 0
    num = 1
    for i in range(b):
        num = num * (a - i) // (i + 1)
    return num


def kite_edges(n):
    """H_n: path v_1..v_n plus the chord v_1 v_3 (1-indexed in the paper,
    0-indexed here)."""
    e = set()
    for i in range(n - 1):
        e.add((i, i + 1))
    e.add((0, 2))
    return sorted(e)


def path_edges(n):
    return [(i, i + 1) for i in range(n - 1)]


def star_edges(k):
    """K_{1,k} on k+1 vertices, centre 0."""
    return [(0, i) for i in range(1, k + 1)]


def adjacency(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def all_distances(n, edges):
    """BFS from every vertex; returns None if the graph is disconnected."""
    adj = adjacency(n, edges)
    dist = [[-1] * n for _ in range(n)]
    for s in range(n):
        dist[s][s] = 0
        frontier = [s]
        while frontier:
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if dist[s][w] < 0:
                        dist[s][w] = dist[s][u] + 1
                        nxt.append(w)
            frontier = nxt
        if any(d < 0 for d in dist[s]):
            return None
    return dist


def transmissions(dist):
    return [sum(row) for row in dist]


def wiener(dist):
    return sum(sum(row) for row in dist) // 2


def distance_laplacian(dist):
    n = len(dist)
    tau = transmissions(dist)
    return [[(tau[i] if i == j else 0) - dist[i][j] for j in range(n)] for i in range(n)]


def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def charpoly(A):
    """det(xI-A) for an integer matrix, by Faddeev-LeVerrier.
    Returns integer coefficients, highest degree first."""
    n = len(A)
    coeffs = [1]
    M = [row[:] for row in A]
    tr = sum(M[i][i] for i in range(n))
    coeffs.append(-tr)
    for k in range(2, n + 1):
        shifted = [row[:] for row in M]
        for i in range(n):
            shifted[i][i] += coeffs[k - 1]
        M = mat_mul(A, shifted)
        tr = sum(M[i][i] for i in range(n))
        assert tr % k == 0, "non-integral trace division in Faddeev-LeVerrier"
        coeffs.append(-(tr // k))
    return coeffs


def poly_from_roots(roots):
    """prod (x - lam), highest degree first."""
    prod = [1]
    for lam in roots:
        new = [0] * (len(prod) + 1)
        for i, c in enumerate(prod):
            new[i] += c
            new[i + 1] -= c * lam
        prod = new
    return prod


def poly_eval(coeffs, x):
    v = 0
    for c in coeffs:
        v = v * x + c
    return v


def taylor_shift(C, a):
    """Coefficients (highest degree first) of P(y+a) from those of P(z)."""
    work = list(C)
    out = []
    while work:
        carry = 0
        bs = []
        for c in work:
            carry = c + carry * a
            bs.append(carry)
        out.append(bs[-1])
        work = bs[:-1]
    out.reverse()
    return out


def sign_variations(C):
    s = [1 if c > 0 else -1 for c in C if c != 0]
    return sum(1 for i in range(1, len(s)) if s[i] != s[i - 1])


class RootCounter(object):
    """Counts roots strictly greater than a/DEN for a polynomial all of whose
    roots are real (guaranteed here: the matrices are symmetric).  For a
    real-rooted polynomial Descartes' rule of signs is an equality, so the
    count is exact and uses integer arithmetic only."""

    def __init__(self, coeffs, den):
        n = len(coeffs) - 1
        self.den = den
        self.scaled = [coeffs[j] * den ** j for j in range(n + 1)]

    def greater(self, a_num):
        return sign_variations(taylor_shift(self.scaled, a_num))


def eig_bracket(counter, k, hi_num, n):
    """Bracket the k-th largest eigenvalue (k=1 is the largest):
    returns integers (lo, hi) with lo/den < lam_k <= hi/den."""
    den = counter.den
    lo, hi = -den, hi_num
    assert counter.greater(lo) >= k and counter.greater(hi) < k
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if counter.greater(mid) >= k:
            lo = mid
        else:
            hi = mid
    return lo, hi


class Spec(object):
    """Everything computed from a graph: distances, transmissions, Wiener
    index, distance Laplacian, exact characteristic polynomial, and certified
    rational brackets for the eigenvalues."""

    def __init__(self, n, edges, prec=24):
        self.n = n
        self.edges = sorted(edges)
        self.dist = all_distances(n, edges)
        if self.dist is None:
            raise ValueError("disconnected graph")
        self.tau = transmissions(self.dist)
        self.W = wiener(self.dist)
        self.DL = distance_laplacian(self.dist)
        self.coeffs = charpoly(self.DL)
        self.den = 1 << prec
        self.counter = RootCounter(self.coeffs, self.den)
        self.hi_num = (2 * self.W + 1) * self.den
        self._br = {}

    def bracket(self, k):
        if k not in self._br:
            self._br[k] = eig_bracket(self.counter, k, self.hi_num, self.n)
        return self._br[k]

    def eig_lo(self, k):
        return Fraction(self.bracket(k)[0], self.den)

    def eig_hi(self, k):
        return Fraction(self.bracket(k)[1], self.den)

    def U_lower(self, r):
        return sum((self.eig_lo(k) for k in range(1, r + 1)), Fraction(0))

    def U_upper(self, r):
        return sum((self.eig_hi(k) for k in range(1, r + 1)), Fraction(0))

    def bound(self, r):
        return Fraction(self.W + binom(r + 2, 3))

    def integer_spectrum(self):
        """The full spectrum as integers, or None if some eigenvalue is
        irrational.  Multiplicities come from exact synthetic division."""
        if hasattr(self, "_ints"):
            return self._ints
        rem = list(self.coeffs)
        found = []
        for t in range(0, 2 * self.W + 1):
            while len(rem) > 1 and poly_eval(rem, t) == 0:
                q, carry = [], 0
                for c in rem:
                    carry = c + carry * t
                    q.append(carry)
                rem = q[:-1]
                found.append(t)
        self._ints = sorted(found, reverse=True) if len(found) == self.n else None
        return self._ints

    def U_bounds(self, r):
        """Rational bounds (lo, hi) with lo <= U_r <= hi; equal when U_r is
        known exactly.  For r >= n-1 the value is 2W, because the kernel of the
        distance Laplacian of a connected graph is one dimensional -- which is
        verified here by root counting rather than assumed."""
        if r >= self.n - 1:
            assert (self.counter.greater(0) == self.n - 1
                    and poly_eval(self.coeffs, 0) == 0), "kernel not simple"
            v = Fraction(2 * self.W)
            return v, v
        ints = self.integer_spectrum()
        if ints is not None:
            v = Fraction(sum(ints[:r]))
            return v, v
        return self.U_lower(r), self.U_upper(r)

    def is_symmetric(self):
        return all(self.DL[i][j] == self.DL[j][i]
                   for i in range(self.n) for j in range(self.n))


def check_kite_structure(nmax=26):
    """H_n is a well-formed connected unicyclic kite: n vertices, n edges,
    a triangle on v1v2v3, one vertex of degree 3, one leaf, girth 3."""
    bad = []
    for n in range(4, nmax + 1):
        e = kite_edges(n)
        deg = [0] * n
        for u, v in e:
            deg[u] += 1
            deg[v] += 1
        dist = all_distances(n, e)
        triangle = all(p in e for p in [(0, 1), (0, 2), (1, 2)])
        degs = sorted(deg)
        want = sorted([2, 2, 3] + [2] * (n - 4) + [1]) if n >= 5 else sorted([2, 2, 3, 1])
        diam = max(max(row) for row in dist) if dist else -1
        if not (len(e) == n and dist is not None and triangle and degs == want
                and diam == n - 2):
            bad.append(n)
    ck("kite_H_n_is_connected_unicyclic_kite", not bad,
       "n=4..%d; H_5 edges %s" % (nmax, kite_edges(5)) if not bad else "bad n=%s" % bad)


def check_dl_construction(nmax=14):
    """The built matrix really is the distance Laplacian: symmetric, row sums
    zero, trace 2W, and x^T D^L x = sum_{u<v} d(u,v)(x_u-x_v)^2."""
    bad = []
    for n in range(4, nmax + 1):
        S = Spec(n, kite_edges(n))
        if not S.is_symmetric():
            bad.append((n, "asym"))
        if any(sum(row) != 0 for row in S.DL):
            bad.append((n, "rowsum"))
        if sum(S.DL[i][i] for i in range(n)) != 2 * S.W:
            bad.append((n, "trace"))
        for seed in range(4):
            x = [((i * i * (seed + 3) + 7 * i + seed) % 11) - 5 for i in range(n)]
            lhs = sum(x[i] * S.DL[i][j] * x[j] for i in range(n) for j in range(n))
            rhs = sum(S.dist[i][j] * (x[i] - x[j]) ** 2
                      for i in range(n) for j in range(i + 1, n))
            if lhs != rhs:
                bad.append((n, "quadform", seed))
    ck("distance_laplacian_construction_identities", not bad,
       "n=4..%d: symmetry, 1 in kernel, trace=2W, quadratic form" % nmax
       if not bad else "bad %s" % bad[:4])


PAPER_DL_H5 = [[7, -1, -1, -2, -3],
               [-1, 7, -1, -2, -3],
               [-1, -1, 5, -1, -2],
               [-2, -2, -1, 6, -1],
               [-3, -3, -2, -1, 9]]
PAPER_W_H5 = 17
PAPER_U3_H5 = 28
PAPER_BOUND_H5 = 27


def check_H5_matrix():
    """The displayed D^L(H_5) is reproduced from the distances of H_5."""
    S = Spec(5, kite_edges(5))
    ok = S.DL == PAPER_DL_H5 and S.W == PAPER_W_H5
    print("      D^L(H_5) recomputed from BFS distances:")
    for row in S.DL:
        print("        " + " ".join("%3d" % v for v in row))
    ck("H5_distance_laplacian_matches_displayed_matrix", ok,
       "W(H_5)=%d, tau=%s" % (S.W, S.tau))


def check_H5_charpoly():
    """det(xI-D^L(H_5)) = x(x-6)(x-8)(x^2-20x+95), expanded independently."""
    S = Spec(5, kite_edges(5))
    factors = [[1, 0], [1, -6], [1, -8], [1, -20, 95]]
    prod = [1]
    for f in factors:
        new = [0] * (len(prod) + len(f) - 1)
        for i, a in enumerate(prod):
            for j, b in enumerate(f):
                new[i + j] += a * b
        prod = new
    ck("H5_characteristic_polynomial_matches_paper_factorisation",
       S.coeffs == prod,
       "computed %s vs factored %s" % (S.coeffs, prod))


def check_H5_spectrum():
    """Spec(D^L(H_5)) = {10+sqrt5, 8, 10-sqrt5, 6, 0}: 0,6,8 are exact roots,
    x^2-20x+95 divides the characteristic polynomial exactly, and 10+-sqrt5
    are its roots in Z[sqrt 5].  The ordering is certified by root counting."""
    S = Spec(5, kite_edges(5))
    exact_roots = all(poly_eval(S.coeffs, t) == 0 for t in (0, 6, 8))
    # (10 +- s)^2 - 20(10 +- s) + 95 with s^2 = 5, computed in Z[s]
    def q_at(a, b):                    # value of x^2-20x+95 at a + b*sqrt5
        return (a * a + 5 * b * b - 20 * a + 95, 2 * a * b - 20 * b)
    surd_ok = q_at(10, 1) == (0, 0) and q_at(10, -1) == (0, 0)
    # 2 < sqrt5 < 3, so 10+sqrt5 > 8 and 6 < 10-sqrt5 < 8: integer facts only
    surd_size_ok = 2 ** 2 < 5 < 3 ** 2
    D = S.den
    counts = [S.counter.greater(t * D) for t in (0, 6, 7, 8, 12, 13)]
    order_ok = counts == [4, 3, 3, 1, 1, 0]
    ck("H5_spectrum_is_10pm_sqrt5_8_6_0",
       exact_roots and surd_ok and surd_size_ok and order_ok,
       "0,6,8 are exact roots; 10+-sqrt5 root of x^2-20x+95 in Z[sqrt5]; "
       "#eigs above (0,6,7,8,12,13) = %s" % (counts,))


def check_H5_violation():
    """LOAD-BEARING: U_3(H_5) exceeds W(H_5) + C(5,3).  U_3 is obtained as
    trace - lambda_4 - lambda_5 with lambda_4 = 6, lambda_5 = 0 certified by
    exact root counting, so the value 28 is computed, not quoted."""
    S = Spec(5, kite_edges(5))
    D = S.den
    four_above_zero = S.counter.greater(0) == 4          # kernel simple
    three_above_six = S.counter.greater(6 * D) == 3      # lambda_4 = 6 exactly
    six_is_root = poly_eval(S.coeffs, 6) == 0
    U3 = 2 * S.W - 6 - 0
    bound = S.W + binom(3 + 2, 3)
    interval_ok = S.U_lower(3) > bound and S.U_upper(3) >= U3
    ok = (four_above_zero and three_above_six and six_is_root
          and U3 == PAPER_U3_H5 and bound == PAPER_BOUND_H5
          and U3 > bound and interval_ok)
    ck("H5_violates_the_bound_at_r_equals_3", ok,
       "U_3 = 2W - lambda_4 - lambda_5 = %d > %d = W + C(5,3); excess %d; "
       "independent bracket sum lower bound %s" % (U3, bound, U3 - bound,
                                                   S.U_lower(3)))


def check_family_gap_r_n_minus_2(nmax=24):
    """For every n in [5,nmax]: Delta_{n-2}(H_n) > 0 and >= (3n-14)/2.
    U_{n-2} = 2W - lambda_{n-1} exactly (the kernel is simple), and
    lambda_{n-1} is bracketed by exact root counting."""
    rows, bad = [], []
    for n in range(5, nmax + 1):
        S = Spec(n, kite_edges(n))
        if S.counter.greater(0) != n - 1:
            bad.append((n, "kernel"))
            continue
        lo, hi = S.bracket(n - 1)
        lam_hi = Fraction(hi, S.den)          # certified upper bound
        lam_lo = Fraction(lo, S.den)
        r = n - 2
        gap_lower = 2 * S.W - lam_hi - S.bound(r)
        gap_upper = 2 * S.W - lam_lo - S.bound(r)
        claimed = Fraction(3 * n - 14, 2)
        if not (gap_lower > 0 and gap_lower >= claimed):
            bad.append((n, "gap", gap_lower, claimed))
        rows.append((n, gap_lower, gap_upper, claimed))
    for n, gl, gu, cl in rows[:4] + rows[-2:]:
        print("      n=%d: Delta_{n-2} in (%s, %s], paper's bound (3n-14)/2 = %s"
              % (n, float(gl), float(gu), cl))
    ck("kite_family_violates_at_r_equals_n_minus_2", not bad,
       "n=5..%d all violate; every computed excess >= (3n-14)/2" % nmax
       if not bad else "bad %s" % bad[:3])


def check_wiener_formula(nmax=60):
    """W(H_n) = C(n+1,3) - (n-2), i.e. W(P_n) drops by exactly n-2."""
    bad = []
    for n in range(4, nmax + 1):
        dist = all_distances(n, kite_edges(n))
        W = wiener(dist)
        Wp = wiener(all_distances(n, path_edges(n)))
        if not (W == binom(n + 1, 3) - (n - 2) and Wp == binom(n + 1, 3)
                and Wp - W == n - 2):
            bad.append(n)
    ck("wiener_index_formula_for_kite_and_path", not bad,
       "n=4..%d: W(H_n)=C(n+1,3)-(n-2), W(P_n)=C(n+1,3)" % nmax
       if not bad else "bad n=%s" % bad)


def check_transmission_formulas(nmax=60):
    """tau(v_3) = (n^2-5n+10)/2 and tau(v_4) = (n^2-7n+22)/2 in H_n."""
    bad = []
    for n in range(5, nmax + 1):
        tau = transmissions(all_distances(n, kite_edges(n)))
        t3, t4 = tau[2], tau[3]
        f3 = Fraction(n * n - 5 * n + 10, 2)
        f4 = Fraction(n * n - 7 * n + 22, 2)
        if not (t3 == f3 and t4 == f4
                and t3 == 2 + Fraction((n - 3) * (n - 2), 2)
                and t4 == 5 + Fraction((n - 4) * (n - 3), 2)):
            bad.append((n, t3, f3, t4, f4))
    ck("transmissions_of_v3_and_v4_match_closed_forms", not bad,
       "n=5..%d: tau(v3)=(n^2-5n+10)/2, tau(v4)=(n^2-7n+22)/2" % nmax
       if not bad else "bad %s" % bad[:3])


def check_rayleigh_bound(nmax=24):
    """The test vector x = e_3 - e_4 gives Rayleigh quotient
    (tau(v3)+tau(v4)+2)/2 = (n^2-6n+18)/2, and this really does dominate the
    computed lambda_{n-1}."""
    bad = []
    for n in range(5, nmax + 1):
        S = Spec(n, kite_edges(n))
        x = [0] * n
        x[2], x[3] = 1, -1
        num = sum(x[i] * S.DL[i][j] * x[j] for i in range(n) for j in range(n))
        rq = Fraction(num, 2)
        closed = Fraction(n * n - 6 * n + 18, 2)
        alt = Fraction(S.tau[2] + S.tau[3] + 2, 2)
        # the certified UPPER bracket: rq >= lam_hi >= lambda_{n-1} is a proof
        # that the test vector really does dominate the eigenvalue, whereas
        # comparing with the lower bracket would prove nothing.
        lam_hi = Fraction(S.bracket(n - 1)[1], S.den)
        perp = sum(x) == 0
        if not (rq == closed == alt and perp and rq >= lam_hi):
            bad.append((n, rq, closed, alt))
    ck("rayleigh_quotient_at_e3_minus_e4_equals_paper_value", not bad,
       "n=5..%d: x perp 1, x^T D^L x / x^T x = (n^2-6n+18)/2 >= lambda_{n-1}"
       % nmax if not bad else "bad %s" % bad[:3])


def selected_set(m):
    """S = {v_1,...,v_m} u {v_{2m+1},...,v_{3m}}, zero-indexed."""
    return list(range(m)) + list(range(2 * m, 3 * m))


def check_selected_transmission_sums(mmax=20):
    """sum_{v in S} tau_{P_3m}(v) = (20m^3-2m)/3 and the same sum for H_{3m}
    is that minus (5m-4), i.e. (20m^3-17m+12)/3."""
    bad = []
    for m in range(2, mmax + 1):
        n = 3 * m
        S = selected_set(m)
        tp = transmissions(all_distances(n, path_edges(n)))
        th = transmissions(all_distances(n, kite_edges(n)))
        sp, sh = sum(tp[v] for v in S), sum(th[v] for v in S)
        expand = 2 * sum(binom(i, 2) + binom(3 * m - i + 1, 2)
                         for i in range(1, m + 1))
        if not (sp == Fraction(20 * m ** 3 - 2 * m, 3) and sp == expand
                and sh == Fraction(20 * m ** 3 - 17 * m + 12, 3)
                and sp - sh == 5 * m - 4):
            bad.append((m, sp, sh, expand))
    ck("selected_transmission_sums_match_closed_forms", not bad,
       "m=2..%d: path sum (20m^3-2m)/3, kite sum (20m^3-17m+12)/3, "
       "difference 5m-4" % mmax if not bad else "bad %s" % bad[:3])


def check_cubic_gap_at_r_2m(mmax=8):
    """For n = 3m, m >= 4: Delta_{2m}(H_{3m}) >= (5m^3-12m^2-17m+12)/6 > 0.
    Two independent certificates: the Ky Fan lower bound from the selected
    transmission sum, and a bracket-sum lower bound on the true U_{2m}."""
    bad, rows = [], []
    for m in range(4, mmax + 1):
        n, r = 3 * m, 2 * m
        S = Spec(n, kite_edges(n))
        sel = sum(S.tau[v] for v in selected_set(m))
        claimed = Fraction(5 * m ** 3 - 12 * m ** 2 - 17 * m + 12, 6)
        kyfan_gap = sel - S.bound(r)
        true_lower = S.U_lower(r) - S.bound(r)
        if not (kyfan_gap == claimed and kyfan_gap > 0
                and S.U_lower(r) >= sel and true_lower > 0
                and true_lower >= claimed):
            bad.append((m, kyfan_gap, claimed, true_lower))
        rows.append((m, claimed, true_lower))
    for m, cl, tl in rows[:3] + rows[-2:]:
        print("      m=%d (n=%d): paper's lower bound %s, computed excess > %s"
              % (m, 3 * m, cl, float(tl)))
    ck("cubic_violation_at_r_equals_2m_for_n_equals_3m", not bad,
       "m=4..%d: selected-sum excess equals (5m^3-12m^2-17m+12)/6 exactly and "
       "the computed U_{2m} excess is at least that" % mmax
       if not bad else "bad %s" % bad[:3])


def check_p_recurrence(mmax=400):
    """p(m) = 5m^3-12m^2-17m+12 satisfies p(4) = 72 and
    p(m+1)-p(m) = 3(m+1)(5m-8), hence p(m) > 0 for all m >= 4."""
    def p(m):
        return 5 * m ** 3 - 12 * m ** 2 - 17 * m + 12
    bad = []
    if p(4) != 72:
        bad.append(("p(4)", p(4)))
    for m in range(1, mmax + 1):
        if p(m + 1) - p(m) != 3 * (m + 1) * (5 * m - 8):
            bad.append(("rec", m))
    for m in range(4, mmax + 1):
        if p(m) <= 0:
            bad.append(("sign", m))
    if p(3) > 0 or p(2) > 0:
        bad.append(("m<4 not covered by the claim", p(2), p(3)))
    ck("p_of_m_recurrence_and_positivity_from_m_equals_4", not bad,
       "p(4)=72, p(m+1)-p(m)=3(m+1)(5m-8) for m<=%d, p>0 exactly from m=4 "
       "(p(3)=%d, p(2)=%d)" % (mmax, p(3), p(2)) if not bad else "bad %s" % bad[:3])


def check_cubic_asymptotics(mmax=40):
    """The lower bound is 5n^3/162 + O(n^2): its distance from 5n^3/162 is at
    most n^2/2, and the excess is O(n^3) from above since Delta_r <= 2W."""
    bad = []
    for m in range(4, mmax + 1):
        n = 3 * m
        claimed = Fraction(5 * m ** 3 - 12 * m ** 2 - 17 * m + 12, 6)
        lead = Fraction(5, 162) * n ** 3
        if not (abs(claimed - lead) * 2 <= n * n):
            bad.append((m, claimed, lead))
        if not Fraction(5 * m ** 3, 6) == lead:
            bad.append((m, "leading coefficient"))
    for m in range(4, 9):
        n, r = 3 * m, 2 * m
        S = Spec(n, kite_edges(n))
        if not S.U_upper(r) - S.bound(r) <= 2 * S.W:
            bad.append((m, "upper side"))
        # the substantive half of the O(n^3) statement: the trace itself, which
        # caps the excess, is at most n^3/3.  Computed from the BFS Wiener
        # index, so a wrong distance matrix breaks it.
        if not 6 * S.W <= n ** 3:
            bad.append((m, "trace is not O(n^3)", S.W))
    ck("cubic_asymptotic_constant_is_5_over_162", not bad,
       "m=4..%d: |lower bound - 5n^3/162| <= n^2/2, and the excess never "
       "exceeds trace = 2W" % mmax if not bad else "bad %s" % bad[:3])


def check_path_remark(nmax=24):
    """Remark: Delta_{n-2}(P_n) >= (3n-10)/2 > 0 for n >= 4.  Also the exact
    arithmetic behind it: C(n,2) - (tau(v2)+tau(v3)+2)/2 = (3n-10)/2."""
    bad, rows = [], []
    for n in range(4, nmax + 1):
        S = Spec(n, path_edges(n))
        if S.counter.greater(0) != n - 1:
            bad.append((n, "kernel"))
            continue
        lam_hi = Fraction(S.bracket(n - 1)[1], S.den)
        gap_lower = 2 * S.W - lam_hi - S.bound(n - 2)
        claimed = Fraction(3 * n - 10, 2)
        arith = Fraction(binom(n, 2)) - Fraction(S.tau[1] + S.tau[2] + 2, 2)
        if not (gap_lower > 0 and gap_lower >= claimed and arith == claimed):
            bad.append((n, gap_lower, claimed, arith))
        rows.append((n, gap_lower, claimed))
    ck("path_family_also_violates_at_r_equals_n_minus_2", not bad,
       "n=4..%d: computed excess >= (3n-10)/2 > 0; e.g. n=4 excess > %s"
       % (nmax, float(rows[0][1])) if not bad else "bad %s" % bad[:3])


def check_path_kyfan_remark(mmax=40, true_mmax=6):
    """Remark: the selected transmission sum for P_{3m} already exceeds
    W(P_{3m}) + C(2m+2,3) for every m >= 3 (and not for m = 2)."""
    bad = []
    for m in range(2, mmax + 1):
        n, r = 3 * m, 2 * m
        tp = transmissions(all_distances(n, path_edges(n)))
        sel = sum(tp[v] for v in selected_set(m))
        W = wiener(all_distances(n, path_edges(n)))
        excess = Fraction(sel - W - binom(r + 2, 3))
        if excess != Fraction(5 * m ** 3 - 12 * m ** 2 - 5 * m, 6):
            bad.append((m, "closed form", excess))
        if m >= 3 and excess <= 0:
            bad.append((m, "not positive", excess))
        if m == 2 and excess > 0:
            bad.append((m, "threshold m>=3 is wrong", excess))
    for m in range(3, true_mmax + 1):
        S = Spec(3 * m, path_edges(3 * m))
        if not S.U_lower(2 * m) > S.bound(2 * m):
            bad.append((m, "true U_2m does not exceed the bound"))
    ck("path_selected_sum_exceeds_bound_exactly_from_m_equals_3", not bad,
       "m=3..%d: excess = (5m^3-12m^2-5m)/6 > 0, and m=2 fails; the true "
       "U_{2m}(P_{3m}) exceeds the bound for m=3..%d" % (mmax, true_mmax)
       if not bad else "bad %s" % bad[:3])


def decide_against_bound(S, r, prec_cap=4000):
    """Decide U_r(G) > W + C(r+2,3), U_r == it, or U_r < it, with proof.

    U_r is a sum of r roots of a monic integer polynomial whose roots all lie
    in [0, 2W], so U_r - B is an algebraic integer of degree at most C(n,r)
    all of whose conjugates are bounded by M = 2Wr + B.  If it is nonzero its
    norm is a nonzero integer, so |U_r - B| >= M^{1-C(n,r)}.  Refining the
    eigenvalue brackets below that separation therefore proves equality when
    the strict comparisons both fail."""
    B = S.bound(r)
    n = S.n
    M = 2 * S.W * r + int(B) + 1
    expo = binom(n, r) - 1
    eps = Fraction(1, M ** expo)
    need = (r * M ** expo).bit_length() + 2
    if need > prec_cap:
        return "undecided", None
    T = Spec(n, S.edges, prec=need)
    lo, hi = T.U_lower(r), T.U_upper(r)
    if lo > B:
        return "violates", lo - B
    if hi <= B:
        return "complies", hi - B
    assert hi - lo < eps, "refinement did not reach the separation bound"
    return "equality", Fraction(0)


def violating_indices(S):
    """All r in [2,n] with U_r(G) > W(G) + C(r+2,3), decided by exact rational
    arithmetic: returns (list_of_r, indeterminate_count)."""
    out, indet = [], 0
    for r in range(2, S.n + 1):
        lo, hi = S.U_bounds(r)
        if lo - S.bound(r) > 0:
            out.append(r)
        elif hi - S.bound(r) <= 0:
            pass
        else:
            verdict, _ = decide_against_bound(S, r)
            if verdict == "violates":
                out.append(r)
            elif verdict == "undecided":
                indet += 1
    return out, indet


def canonical(n, edges):
    """Canonical form of a small graph: lexicographic minimum over all
    relabellings.  Used only to name and to group isomorphic graphs."""
    best = None
    verts = list(range(n))
    for perm in permutations(verts):
        img = tuple(sorted(tuple(sorted((perm[u], perm[v]))) for u, v in edges))
        if best is None or img < best:
            best = img
    return best


def enumerate_connected(n, exact_edges=None):
    """All connected labelled graphs on n vertices (optionally with a fixed
    number of edges), as edge lists."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if exact_edges is None:
        subsets = (
            [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
            for mask in range(1 << len(pairs)))
    else:
        subsets = (list(c) for c in combinations(pairs, exact_edges))
    for edges in subsets:
        if len(edges) < n - 1:
            continue
        if all_distances(n, edges) is not None:
            yield edges


def census(n, exact_edges=None):
    """Exhaustive census of connected graphs on n vertices, grouped by exact
    isomorphism class.  Returns (classes, indeterminate) where classes maps a
    canonical form to (representative edge list, violating r values, diam)."""
    classes, indet = {}, 0
    for edges in enumerate_connected(n, exact_edges):
        key = canonical(n, edges)
        if key in classes:
            continue
        S = Spec(n, edges)
        rs, bad = violating_indices(S)
        indet += bad
        diam = max(max(row) for row in S.dist)
        classes[key] = (sorted(edges), rs, diam)
    return classes, indet


def check_minimal_unicyclic_order():
    """Exhaustive: no connected unicyclic graph on 3 or 4 vertices violates the
    bound for any r, and one on 5 vertices does.  Hence the minimum order of a
    unicyclic counterexample is exactly 5, as the paper claims for H_5."""
    detail, bad, indet = [], [], 0
    winners5 = []
    for n in (3, 4, 5, 6):
        classes, ind = census(n, exact_edges=n)
        indet += ind
        viol = {k: v for k, v in classes.items() if v[1]}
        detail.append("order %d: %d unicyclic iso classes, %d violating"
                      % (n, len(classes), len(viol)))
        if n in (3, 4) and viol:
            bad.append((n, "unexpected violator", list(viol.values())[:1]))
        if n == 5:
            winners5 = [(v[0], v[1]) for v in viol.values()]
            if not viol:
                bad.append((5, "no violating unicyclic class on 5 vertices"))
        if n == 6 and not viol:
            bad.append((6, "no violator on 6 vertices"))
    h5 = canonical(5, kite_edges(5))
    ok_h5 = any(canonical(5, e) == h5 and 3 in rs for e, rs in winners5)
    ck("minimum_order_of_a_unicyclic_counterexample_is_five",
       not bad and indet == 0 and ok_h5,
       "; ".join(detail) + "; H_5 is among the order-5 violators, at r=3"
       if not bad else "bad %s" % bad[:2])


def check_small_unicyclic_spectra():
    """The paper's spectra for C_3, C_4 and H_4 = Ki_{4,3}, checked as exact
    characteristic polynomials, and none of them violates the bound."""
    cases = [("C_3", 3, [(0, 1), (1, 2), (0, 2)], [3, 3, 0]),
             ("C_4", 4, [(0, 1), (1, 2), (2, 3), (0, 3)], [6, 6, 4, 0]),
             ("H_4", 4, kite_edges(4), [7, 5, 4, 0])]
    bad, detail = [], []
    for name, n, edges, spec_claimed in cases:
        S = Spec(n, edges)
        prod = poly_from_roots(spec_claimed)
        rs, ind = violating_indices(S)
        if S.coeffs != prod or rs or ind:
            bad.append((name, S.coeffs, prod, rs, ind))
        detail.append("%s spectrum %s" % (name, spec_claimed))
    ck("spectra_of_C3_C4_H4_and_their_compliance", not bad,
       "; ".join(detail) + "; none violates for any 2<=r<=n"
       if not bad else "bad %s" % bad[:2])


def check_prior_star_counterexamples():
    """The two counterexamples attributed to earlier work really do fail:
    K_{1,3} at r=2 and K_{1,4} at r=3, with the excess computed."""
    detail, bad = [], []
    for k, r in ((3, 2), (4, 3)):
        S = Spec(k + 1, star_edges(k))
        lo, hi = S.U_bounds(r)
        rs, ind = violating_indices(S)
        if not (lo == hi and lo > S.bound(r) and r in rs and ind == 0):
            bad.append((k, r, lo, hi, S.bound(r)))
        detail.append("K_{1,%d}: U_%d = %s > %s = W + C(%d,3), violating r = %s"
                      % (k, r, lo, S.bound(r), r + 2, rs))
    ck("prior_star_counterexamples_reproduce", not bad, "; ".join(detail))


def check_diameter_two_census(nmax=6):
    """Exhaustive over all connected graphs of order <= nmax: the only
    diameter-two violating pairs are (K_{1,3}, 2) and (K_{1,4}, 3), which is
    what the classification quoted in the paper asserts."""
    diam2, others, indet, counted = [], [], 0, 0
    for n in range(3, nmax + 1):
        seen = set()
        for edges in enumerate_connected(n):
            S = Spec(n, edges)
            diam = max(max(row) for row in S.dist)
            # Whether the bound fails depends only on the spectrum and W, both
            # of which the characteristic polynomial fixes; the diameter is
            # carried in the key so that no diameter class can be masked.
            key = (tuple(S.coeffs), diam)
            if key in seen:
                continue
            seen.add(key)
            counted += 1
            rs, ind = violating_indices(S)
            indet += ind
            if rs:
                (diam2 if diam == 2 else others).append(
                    (n, canonical(n, edges), tuple(rs)))
    want = set([(4, canonical(4, star_edges(3)), (2,)),
                (5, canonical(5, star_edges(4)), (3,))])
    got = set(diam2)
    ok = indet == 0 and got == want
    CENSUS_STATS["order_max"] = nmax
    CENSUS_STATS["classes_searched"] = counted
    CENSUS_STATS["diam2_violators"] = len(got)
    CENSUS_STATS["nonD2_violators"] = len(others)
    CENSUS_STATS["nonD2_orders"] = sorted(set(n for n, _, _ in others))
    print("      diameter-two violating classes found: %d" % len(got))
    print("      violators of diameter >= 3 found: %d classes, orders %s"
          % (CENSUS_STATS["nonD2_violators"], CENSUS_STATS["nonD2_orders"]))
    ck("diameter_two_violators_are_exactly_the_two_stars", ok,
       "%d spectral classes of connected graphs of order 3..%d searched, "
       "0 undecided; diameter-two violators = K_{1,3} at r=2 and K_{1,4} at r=3"
       % (counted, nmax))


def check_psd_and_simple_kernel(nmax=24):
    """D^L(H_n) is positive semidefinite with a one-dimensional kernel: all
    roots lie in [0, 2W], exactly n-1 of them are positive, and 0 is a simple
    root.  This is what licenses U_{n-2} = 2W - lambda_{n-1}."""
    bad = []
    for n in range(4, nmax + 1):
        S = Spec(n, kite_edges(n))
        D = S.den
        all_above_minus_one = S.counter.greater(-D) == n
        pos = S.counter.greater(0) == n - 1
        zero_root = poly_eval(S.coeffs, 0) == 0
        simple = S.coeffs[-2] != 0        # coefficient of x in det(xI - D^L)
        top = S.counter.greater((2 * S.W + 1) * D) == 0
        if not (all_above_minus_one and pos and zero_root and simple and top):
            bad.append((n, all_above_minus_one, pos, zero_root, simple, top))
    ck("distance_laplacian_is_psd_with_simple_kernel", not bad,
       "n=4..%d: n eigenvalues in [0,2W], exactly n-1 positive, 0 simple" % nmax
       if not bad else "bad %s" % bad[:3])


def check_gap_algebra_identities(nmax=60):
    """The final algebraic steps of the proof are polynomial identities, so
    checking them at more points than their degree verifies them for ALL n:

      [C(n+1,3)-(n-2)] - C(n,3) - (n^2-6n+18)/2  ==  (3n-14)/2         (kite)
      [(20m^3-17m+12)/3] - [C(3m+1,3)-(3m-2)] - C(2m+2,3)
                                     == (5m^3-12m^2-17m+12)/6          (n=3m)
      C(n,2) - (tau(v_2)+tau(v_3)+2)/2           ==  (3n-10)/2          (paths)
    Both sides of each are polynomials of degree at most 3, and they are
    checked at far more than 4 points."""
    bad = []
    for n in range(4, nmax + 1):
        lhs = (Fraction(binom(n + 1, 3) - (n - 2)) - binom(n, 3)
               - Fraction(n * n - 6 * n + 18, 2))
        if lhs != Fraction(3 * n - 14, 2):
            bad.append(("kite", n, lhs))
    for m in range(2, nmax + 1):
        lhs = (Fraction(20 * m ** 3 - 17 * m + 12, 3)
               - (binom(3 * m + 1, 3) - (3 * m - 2)) - binom(2 * m + 2, 3))
        if lhs != Fraction(5 * m ** 3 - 12 * m ** 2 - 17 * m + 12, 6):
            bad.append(("cubic", m, lhs))
    for n in range(4, nmax + 1):
        tau = transmissions(all_distances(n, path_edges(n)))
        lhs = Fraction(binom(n, 2)) - Fraction(tau[1] + tau[2] + 2, 2)
        if lhs != Fraction(3 * n - 10, 2):
            bad.append(("path", n, lhs))
    ck("closing_algebra_holds_as_polynomial_identities", not bad,
       "all three identities hold at n,m = 2..%d, i.e. at more points than "
       "their degree, hence identically" % nmax
       if not bad else "bad %s" % bad[:3])


KITE_NMAX = 40          # largest kite verified spectrally, r = n-2
CUBIC_MMAX = 12         # largest m verified spectrally at r = 2m, n = 3m
CENSUS_NMAX = 6         # largest order for the exhaustive graph census
PATH_KYFAN_MMAX = 40    # largest m for the path selected-sum closed form
PATH_TRUE_MMAX = 6      # largest m with the true U_{2m}(P_{3m}) bracketed


def deflation_context():
    """The paper's Remark 2 is deflationary: it concedes that P_n violates the
    bound for every n >= 4 and that P_{3m} alone already gives a Theta(n^3)
    excess, so unicyclicity is the only one of the three headline properties
    that paths do not already supply.  Two checks above test exactly those
    concessions, and the census above quantifies how far from rare the failure
    is on small orders.  Print that, with derived counts, so no reader has to
    reconstruct it."""
    keys = ("nonD2_violators", "nonD2_orders", "order_max")
    if not all(k in CENSUS_STATS for k in keys):
        ck("census_counts_recorded_for_the_closing_disclosure", False,
           "the exhaustive census did not record its violator counts, so the "
           "deflationary context below could not be derived")
        return
    print("CONTEXT: this program corroborates the paper's deflationary "
          "Remark 2 as well as its Theorem 1.  Checks "
          "path_family_also_violates_at_r_equals_n_minus_2 and "
          "path_selected_sum_exceeds_bound_exactly_from_m_equals_3 confirm "
          "that the path P_n already violates the bound for every n >= 4 and "
          "that P_{3m} alone already exceeds it by Theta(n^3); and the "
          "exhaustive census found %d violating spectral classes of diameter "
          ">= 3 (keyed by characteristic polynomial and diameter) on orders "
          "%s, i.e. on at most %d vertices.  So, as Remark 2 itself states, "
          "the unicyclic restriction -- not infinitude and not the cubic "
          "order of the excess -- is what separates Theorem 1 from the "
          "unrestricted statement."
          % (CENSUS_STATS["nonD2_violators"], CENSUS_STATS["nonD2_orders"],
             CENSUS_STATS["order_max"]))


def not_rerun():
    print("NOT RE-RUN: the theorem is an infinite family; the spectral "
          "verification above covers H_n for n = 5..%d and the cubic index "
          "r = 2m for m = 4..%d only." % (KITE_NMAX, CUBIC_MMAX))
    print("NOT RE-RUN: the remark's path statements are checked on finite "
          "ranges too: the excess of P_n at r = n-2 spectrally for n = 4..%d, "
          "the selected-sum excess of P_{3m} in closed form for m = 2..%d, "
          "and the true U_{2m}(P_{3m}) against the bound only for "
          "m = 3..%d." % (KITE_NMAX, PATH_KYFAN_MMAX, PATH_TRUE_MMAX))
    print("NOT RE-RUN: the inequality under attack is transcribed from its "
          "published source and is never collated with it, because this "
          "program has no network access; every check tests "
          "U_r(G) <= W(G) + C(r+2,3) for 2 <= r <= n exactly as the paper "
          "quotes it, so a hypothesis dropped in transcription would not be "
          "detected here.  The only offline corroboration is the check "
          "diameter_two_violators_are_exactly_the_two_stars, which "
          "reproduces for this same inequality the diameter-two "
          "classification quoted from earlier work.")
    print("NOT RE-RUN: the closing algebra is verified as a polynomial "
          "identity (hence for all n and m), but the two variational "
          "principles used by the proof, Rayleigh-Ritz and Ky Fan, are "
          "classical inputs and are exercised numerically, not proved here.")
    print("NOT RE-RUN: the exhaustive census covers connected graphs of order "
          "at most %d; the diameter-two classification quoted from earlier work "
          "is confirmed only within that range." % CENSUS_NMAX)


def run_all():
    check_kite_structure()
    check_dl_construction()
    check_psd_and_simple_kernel(KITE_NMAX)
    check_H5_matrix()
    check_H5_charpoly()
    check_H5_spectrum()
    check_H5_violation()
    check_family_gap_r_n_minus_2(KITE_NMAX)
    check_wiener_formula()
    check_transmission_formulas()
    check_rayleigh_bound(KITE_NMAX)
    check_selected_transmission_sums()
    check_cubic_gap_at_r_2m(CUBIC_MMAX)
    check_cubic_asymptotics()
    check_p_recurrence()
    check_gap_algebra_identities()
    check_path_remark(KITE_NMAX)
    check_path_kyfan_remark(PATH_KYFAN_MMAX, PATH_TRUE_MMAX)
    check_small_unicyclic_spectra()
    check_prior_star_counterexamples()
    check_minimal_unicyclic_order()
    check_diameter_two_census(CENSUS_NMAX)
    deflation_context()
    not_rerun()


if __name__ == "__main__":
    try:
        run_all()
    except Exception as exc:                 # a graph that fails a structural
        ck("verification_ran_to_completion",  # precondition must still be
           False,                             # reported, not raised
           "aborted: %s: %s" % (type(exc).__name__, exc))
    finish()

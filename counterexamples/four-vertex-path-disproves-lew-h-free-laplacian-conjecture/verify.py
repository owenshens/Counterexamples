#!/usr/bin/env python3
"""Verification program for a refutation of an H-free Laplacian-gap conjecture.

Standard library only; every decision is made in exact integer or rational
arithmetic (or in exact arithmetic in Q(sqrt 2)).  No float ever decides a check.

Notation.  eps_k(G) = (sum of the k largest Laplacian eigenvalues of G) - |E(G)|,
and ex(n;H) is the largest number of edges of an n-vertex H-free graph.  The
conjecture under attack asserts eps_k(G) <= ex(k+1;H) for every k >= 1, every
graph H, and every H-free graph G with |V(G)| >= k.

The paper carries one shared theorem counter, so its statements typeset as
Theorem 1 (the P_4 refutation), Lemma 2 (pendant amplification), Proposition 3,
Corollary 4 (every H with delta(H) >= 2), Corollary 5 (every r >= 3, k >= r-1),
Remark 6 (the size of the violation and the approximate bound), Proposition 7
(minimality) and Remark 8 (H = K_2 + K_1).  Those are the numbers used below.

TAKEN FROM THE PAPER (inputs, transcribed; nothing here is checked against
itself):
  * the verbatim statement of Conjecture 7.1 of arXiv:2601.17575v1, as quoted
    in the paper, and the approximate bound ex(k;H) + (4k-2)sqrt(k) that the
    paper attributes to Section 7 of the same preprint.  This program has no
    access to that preprint and CANNOT corroborate either quotation; a
    transcription error in the conjecture's hypotheses on H would not be
    detected here.  See the closing NOT RE-RUN line;
  * the exhibited graph: the path P_4 on four vertices, edges 12, 23, 34;
  * the parameters of the refutation: k = 2, H = K_3, bound value 2;
  * the claimed characteristic polynomial t(t-2)(t^2-4t+2) and the claimed
    Laplacian spectrum 2+sqrt2, 2, 2-sqrt2, 0;
  * the claimed value eps_2(P_4) = 1 + sqrt2;
  * the pendant-amplification lemma: for F on k vertices with e edges and
    Laplacian eigenvalues mu_i, the graph F^<N> (N fresh pendants on every
    vertex of F) has Laplacian spectrum {1 with multiplicity k(N-1)} together
    with the roots alpha_i < 1 < beta_i of t^2-(N+1+mu_i)t+mu_i, and
    eps_k(F^<N>) = k + e - sum alpha_i > k + e - 2e/(N+1);
  * the claimed minimality: no counterexample with k = 1 or with |V(G)| <= k+1
    when H has an edge and no isolated vertices, so k >= 2 and |V(G)| >= 4;
  * the claimed remark that H = K_2 + K_1, G = K_2, k = 2 already fails.

DERIVED HERE (recomputed from the objects themselves):
  * P_4 decoded, counted and printed back; degrees, connectivity, path shape;
  * triangle-freeness by clique search and, independently, by injection search;
  * det(tI - L(P_4)) by Faddeev-LeVerrier over the integers, compared with the
    expansion of the paper's factorisation;
  * the claimed spectrum verified to annihilate that polynomial in Q(sqrt 2),
    to be strictly decreasing, and to sum to the trace;
  * an independent enclosure of every eigenvalue by exact rational bisection on
    root counts (Descartes' rule applied to shifted polynomials, exact because
    Laplacian characteristic polynomials are real-rooted);
  * ex(3;K_3), ex(2;K_3), ex(3;K_2+K_1), ex(n;K_r) for n <= 6 and r in {3,4,5},
    and ex(4;C_4) by exhaustive enumeration of all labelled graphs;
  * THE LOAD-BEARING STEP: eps_2(P_4) > ex(3;K_3), decided from a rational
    lower bound obtained from P_4's own Laplacian, plus the exact gap sqrt2 - 1;
  * minimality: eps_1(G) <= 1 for all 33867 labelled graphs on <= 6 vertices;
    eps_k(G) = |E(G)| <= ex(k+1;K_3) for all graphs on <= 5 vertices with
    k >= |V(G)|-1; no triangle-free graph on <= 3 vertices violates at k = 2;
    and a full census of 4-vertex triangle-free graphs at k = 2;
  * the pendant lemma as an exact polynomial identity, the multiplicity of the
    eigenvalue 1, the separation of the beta_i, and the strict inequality;
  * the two corollaries on concrete instances, including the excess bound;
  * Remark 6's consistency assertion, eps_k(G) <= ex(k;H) + (4k-2)sqrt(k), for
    every counterexample instance this program builds, decided exactly by
    squaring (no float) from a rational UPPER bound on eps_k(G).
"""

import itertools
from fractions import Fraction

CHECKS = []

# Every counterexample instance this program builds is recorded here as it is
# built, for the Remark 6 test against the paper's transcribed approximate bound.
# Entries are (label, k, eps_hi, ex_k_H): eps_hi is a RATIONAL UPPER BOUND for
# eps_k(G) and ex_k_H is the exact value of ex(k;H), both produced upstream.  The
# consumer checks the expected number of entries, so an instance that stops being
# registered turns into a FAILED check instead of a silently empty loop.
EXPECTED_INSTANCES = 6
APPROX_INSTANCES = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ((" [" + detail + "]") if detail else ""))
    return bool(ok)


# ----------------------------------------------------------------- graph basics
def laplacian(n, edges):
    """Integer Laplacian D - A of the simple graph (n, edges)."""
    M = [[0] * n for _ in range(n)]
    for (u, v) in edges:
        if u == v or not (0 <= u < n and 0 <= v < n):
            raise ValueError("not a simple graph")
        M[u][v] -= 1
        M[v][u] -= 1
        M[u][u] += 1
        M[v][v] += 1
    return M


def degrees(n, edges):
    d = [0] * n
    for (u, v) in edges:
        d[u] += 1
        d[v] += 1
    return d


def is_connected(n, edges):
    if n == 0:
        return True
    adj = [set() for _ in range(n)]
    for (u, v) in edges:
        adj[u].add(v)
        adj[v].add(u)
    seen = {0}
    stack = [0]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == n


def nbr_masks(n, edges):
    m = [0] * n
    for (u, v) in edges:
        m[u] |= 1 << v
        m[v] |= 1 << u
    return m


def has_clique(n, edges, r):
    """True iff the graph contains K_r as a subgraph."""
    if r <= 1:
        return n >= r
    mask = nbr_masks(n, edges)

    def grow(clique, cand):
        if len(clique) == r:
            return True
        c = cand
        while c:
            low = c & -c
            v = low.bit_length() - 1
            c ^= low
            if grow(clique + [v], cand & mask[v] & ~((low << 1) - 1)):
                return True
        return False

    return grow([], (1 << n) - 1)


def contains_copy(n, edges, hn, hedges):
    """True iff (n, edges) contains (hn, hedges) as a (not nec. induced) subgraph."""
    if hn > n:
        return False
    eset = set()
    for (u, v) in edges:
        eset.add((min(u, v), max(u, v)))
    for inj in itertools.permutations(range(n), hn):
        if all((min(inj[a], inj[b]), max(inj[a], inj[b])) in eset for (a, b) in hedges):
            return True
    return False


# ------------------------------------------- exact characteristic polynomial etc.
def charpoly(A):
    """Faddeev-LeVerrier: exact integer coefficients of det(tI - A), low->high."""
    n = len(A)
    c = [0] * (n + 1)
    c[n] = 1
    M = [[0] * n for _ in range(n)]
    for k in range(1, n + 1):
        if k == 1:
            M = [row[:] for row in A]
        else:
            B = [row[:] for row in M]
            for i in range(n):
                B[i][i] += c[n - k + 1]
            M = [[sum(A[i][t] * B[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        tr = sum(M[i][i] for i in range(n))
        num = -tr
        if num % k != 0:
            raise ArithmeticError("non-integral charpoly coefficient")
        c[n - k] = num // k
    return c


def poly_shift(p, x):
    """Coefficients of p(s + x), exact (x an int or Fraction)."""
    res = [0]
    for coef in reversed(p):
        new = [0] * (len(res) + 1)
        for i, a in enumerate(res):
            new[i + 1] += a
            new[i] += a * x
        new[0] += coef
        res = new
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res


def count_gt(p, x):
    """Number of roots of p (all real) strictly greater than x, by Descartes'
    rule of signs applied to p(s+x) -- exact for real-rooted polynomials."""
    q = poly_shift(p, x)
    i = 0
    while i < len(q) and q[i] == 0:
        i += 1
    q = q[i:]
    if not q:
        return 0
    changes = 0
    last = 0
    for coef in q:
        if coef == 0:
            continue
        s = 1 if coef > 0 else -1
        if last and s != last:
            changes += 1
        last = s
    return changes


def eig_bracket(p, i, iters=48):
    """Rational (lo, hi) with lo < mu_i <= hi, mu_i the i-th largest root."""
    deg = len(p) - 1
    lo = Fraction(-1)
    hi = Fraction(deg + 1)
    for _ in range(iters):
        mid = (lo + hi) / 2
        if count_gt(p, mid) >= i:
            lo = mid
        else:
            hi = mid
    return lo, hi


def top_k_sum_bracket(p, k, iters=48):
    lo = Fraction(0)
    hi = Fraction(0)
    for i in range(1, k + 1):
        a, b = eig_bracket(p, i, iters)
        lo += a
        hi += b
    return lo, hi


def eps_k_bracket(n, edges, k, iters=48):
    """Rational enclosure of eps_k(G) = sum of k largest Laplacian eigenvalues - |E|."""
    p = charpoly(laplacian(n, edges))
    lo, hi = top_k_sum_bracket(p, k, iters)
    m = len(edges)
    return lo - m, hi - m


# -------------------------------------------- exact arithmetic in Q(sqrt(2)) and
# -------------------------------------------- exact Turan / extremal numbers
class Q2(object):
    """a + b*sqrt(2) with a, b rational; all comparisons exact."""

    def __init__(self, a, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(self.a + o.a, self.b + o.b)

    def __sub__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(self.a - o.a, self.b - o.b)

    def __mul__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(self.a * o.a + 2 * self.b * o.b, self.a * o.b + self.b * o.a)

    def sign(self):
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        if self.a > 0 and self.b > 0:
            return 1
        if self.a < 0 and self.b < 0:
            return -1
        # opposite signs: compare a^2 with 2 b^2
        d = self.a * self.a - 2 * self.b * self.b
        if d == 0:
            return 0
        return 1 if (d > 0) == (self.a > 0) else -1

    def __eq__(self, o):
        return (self - (o if isinstance(o, Q2) else Q2(o))).sign() == 0

    def __gt__(self, o):
        return (self - (o if isinstance(o, Q2) else Q2(o))).sign() > 0

    def __repr__(self):
        return "%s + %s*sqrt2" % (self.a, self.b)

    def bracket(self):
        """Rational enclosure, from 1.414213562 < sqrt(2) < 1.414213563."""
        r_lo = Fraction(1414213562, 10 ** 9)
        r_hi = Fraction(1414213563, 10 ** 9)
        vals = [self.a + self.b * r_lo, self.a + self.b * r_hi]
        return min(vals), max(vals)


def all_graphs(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(pairs)):
        yield [pairs[t] for t in range(len(pairs)) if mask >> t & 1]


def ex_number_bruteforce(n, hn, hedges, clique_r=None):
    """Max edges of an n-vertex graph containing no copy of H, by full enumeration."""
    best = -1
    for edges in all_graphs(n):
        if clique_r is not None:
            bad = has_clique(n, edges, clique_r)
        else:
            bad = contains_copy(n, edges, hn, hedges)
        if not bad and len(edges) > best:
            best = len(edges)
    return best


def turan_edges(n, r):
    """|E(T_{r-1}(n))|, the balanced complete (r-1)-partite graph on n vertices."""
    parts = r - 1
    sizes = [n // parts + (1 if i < n % parts else 0) for i in range(parts)]
    return (n * n - sum(s * s for s in sizes)) // 2


def turan_graph(n, r):
    parts = r - 1
    lab = []
    sizes = [n // parts + (1 if i < n % parts else 0) for i in range(parts)]
    for i, s in enumerate(sizes):
        lab.extend([i] * s)
    edges = [(u, v) for u in range(n) for v in range(u + 1, n) if lab[u] != lab[v]]
    return n, edges


def pendant_blowup(n, edges, N):
    """F^{<N>}: attach N fresh pendant vertices to every vertex of F."""
    out = list(edges)
    nxt = n
    for v in range(n):
        for _ in range(N):
            out.append((v, nxt))
            nxt += 1
    return nxt, out


# ------------------------------------------------------ the exhibited object
# P_4, the path on four vertices 1-2-3-4 (transcribed from the paper).
P4_N = 4
P4_EDGES = [(0, 1), (1, 2), (2, 3)]
K3 = (3, [(0, 1), (0, 2), (1, 2)])
K_PARAM = 2          # the parameter k of the conjecture
PAPER_BOUND = 2      # ex(3; K_3) as stated in the paper
SQRT2 = Q2(0, 1)


def check_object():
    n, E = P4_N, P4_EDGES
    deg = degrees(n, E)
    L = laplacian(n, E)
    print("object: n=%d, |E|=%d, edges=%s" % (n, len(E), sorted((u + 1, v + 1) for (u, v) in E)))
    print("object: degree sequence %s, Laplacian rows %s" % (sorted(deg), L))
    is_path = (
        len(E) == n - 1
        and is_connected(n, E)
        and sorted(deg) == [1, 1, 2, 2]
        and max(deg) == 2
    )
    ck("exhibited_graph_is_P4", n == 4 and len(E) == 3 and is_path,
       "4 vertices, 3 edges, connected, degrees 1,2,2,1")
    ck("P4_edge_count_matches_paper_3", len(set((min(a, b), max(a, b)) for a, b in E)) == 3,
       "|E(P_4)| = %d" % len(E))
    tri_clique = has_clique(n, E, 3)
    tri_inject = contains_copy(n, E, K3[0], K3[1])
    ck("P4_is_triangle_free", (not tri_clique) and (not tri_inject),
       "no K_3 by clique search and by injection search")


def check_charpoly():
    L = laplacian(P4_N, P4_EDGES)
    got = charpoly(L)
    # the paper's factorisation t(t-2)(t^2-4t+2), expanded here from its factors
    def mul(p, q):
        r = [0] * (len(p) + len(q) - 1)
        for i, a in enumerate(p):
            for j, b in enumerate(q):
                r[i + j] += a * b
        return r
    want = mul(mul([0, 1], [-2, 1]), [2, -4, 1])
    print("charpoly(L(P_4)) low->high = %s ; paper factorisation expands to %s" % (got, want))
    ck("P4_charpoly_equals_paper_factorisation", got == want,
       "t^4-6t^3+10t^2-4t from det(tI-L) and from t(t-2)(t^2-4t+2)")


def poly_eval_q2(p, x):
    acc = Q2(0)
    for coef in reversed(p):
        acc = acc * x + Q2(coef)
    return acc


def check_spectrum():
    p = charpoly(laplacian(P4_N, P4_EDGES))
    claimed = [Q2(2) + SQRT2, Q2(2), Q2(2) - SQRT2, Q2(0)]
    zeros = [poly_eval_q2(p, x).sign() == 0 for x in claimed]
    distinct = all(
        (claimed[i] - claimed[j]).sign() != 0
        for i in range(4) for j in range(i + 1, 4)
    )
    ck("P4_claimed_spectrum_are_the_roots",
       all(zeros) and distinct and len(p) - 1 == 4,
       "2+r2, 2, 2-r2, 0 all annihilate the degree-4 charpoly and are distinct")
    tr = sum(laplacian(P4_N, P4_EDGES)[i][i] for i in range(P4_N))
    s = Q2(0)
    for x in claimed:
        s = s + x
    ck("P4_spectrum_trace_identity", s.sign() > 0 and (s - Q2(tr)).sign() == 0 and tr == 2 * len(P4_EDGES),
       "sum of claimed eigenvalues = trace L = 2|E| = %d" % tr)
    ordered = all((claimed[i] - claimed[i + 1]).sign() > 0 for i in range(3))
    # the ordering must be the ordering of the COMPUTED spectrum, not merely of the
    # transcribed literals: the i-th claimed value must have rank i among the roots.
    ranks_ok = True
    for i, x in enumerate(claimed, start=1):
        xl, xh = x.bracket()
        above = count_gt(p, xh)
        at_least = count_gt(p, xl) + root_multiplicity(p, xl)
        if above != i - 1 or at_least < i:
            ranks_ok = False
    ck("P4_spectrum_strictly_decreasing", ordered and ranks_ok,
       "2+r2 > 2 > 2-r2 > 0, and each is the root of that rank in det(tI-L)")
    # independent route: rational bisection brackets on the same charpoly
    ok = True
    detail = []
    for i, x in enumerate(claimed, start=1):
        lo, hi = eig_bracket(p, i)
        xl, xh = x.bracket()
        ok = ok and (lo <= xh) and (xl <= hi)
        detail.append("mu_%d in (%.9f, %.9f]" % (i, float(lo), float(hi)))
    ck("P4_bisection_brackets_agree_with_closed_forms", ok, "; ".join(detail))


def check_eps2_value():
    p = charpoly(laplacian(P4_N, P4_EDGES))
    m = len(P4_EDGES)
    # exact value from the two largest eigenvalues, in Q(sqrt 2).  The two summands
    # are first shown to BE the top two roots of the computed charpoly, so this is
    # not an identity among transcribed constants.
    top2 = [Q2(2) + SQRT2, Q2(2)]
    roots_ok = all(poly_eval_q2(p, x).sign() == 0 for x in top2)
    tl, th = Q2(2).bracket()
    rank_ok = (count_gt(p, th) == 1
               and count_gt(p, tl) + root_multiplicity(p, tl) == 2)
    exact = top2[0] + top2[1] - Q2(m)
    ck("eps2_P4_equals_1_plus_sqrt2",
       roots_ok and rank_ok and (exact - (Q2(1) + SQRT2)).sign() == 0,
       "top two roots of det(tI-L) are 2+r2 and 2, and (2+r2) + 2 - 3 = 1 + r2 exactly")
    lo, hi = eps_k_bracket(P4_N, P4_EDGES, K_PARAM)
    xl, xh = (Q2(1) + SQRT2).bracket()
    print("eps_2(P_4) enclosure (%.12f, %.12f]; 1+sqrt2 in (%.12f, %.12f)"
          % (float(lo), float(hi), float(xl), float(xh)))
    ck("eps2_P4_enclosure_contains_1_plus_sqrt2",
       lo <= xh and xl <= hi and (hi - lo) < Fraction(1, 10 ** 6),
       "independent bisection enclosure of width < 1e-6 straddles 1+sqrt2")


def check_bound_and_refutation():
    ex3 = ex_number_bruteforce(3, K3[0], K3[1], clique_r=3)
    ck("ex_3_K3_equals_2_by_enumeration", ex3 == PAPER_BOUND,
       "brute force over all 8 graphs on 3 vertices gives ex(3;K_3) = %d" % ex3)
    # hypotheses of the conjecture, as literally stated
    hyp_k = K_PARAM >= 1
    hyp_v = P4_N >= K_PARAM
    hyp_free = not has_clique(P4_N, P4_EDGES, 3)
    ck("conjecture_hypotheses_all_met", hyp_k and hyp_v and hyp_free,
       "k = %d >= 1, |V(G)| = %d >= k, G is K_3-free" % (K_PARAM, P4_N))
    # THE load-bearing computation: eps_2(P_4) > ex(3;K_3), from a rational
    # lower bound produced by root counting on the graph's own charpoly.
    lo, hi = eps_k_bracket(P4_N, P4_EDGES, K_PARAM)
    strict = lo > Fraction(ex3)
    ck("REFUTATION_eps2_P4_strictly_exceeds_ex3_K3", strict,
       "%.12f > %d, margin >= %.12f" % (float(lo), ex3, float(lo - ex3)))
    exact_gap = (Q2(1) + SQRT2) - Q2(ex3)
    ck("REFUTATION_exact_gap_is_sqrt2_minus_1_positive",
       exact_gap.sign() > 0 and (exact_gap - (SQRT2 - Q2(1))).sign() == 0,
       "eps_2(P_4) - ex(3;K_3) = sqrt2 - 1 > 0")
    APPROX_INSTANCES.append(("P_4 (Theorem 1), H=K_3", K_PARAM, hi,
                             ex_number_bruteforce(K_PARAM, K3[0], K3[1], clique_r=3)))


def pmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                r[i + j] += a * b
    return r


def ppow(p, e):
    r = [1]
    for _ in range(e):
        r = pmul(r, p)
    return r


def root_multiplicity(p, x):
    """Multiplicity of x as a root of p."""
    q = poly_shift(p, x)
    m = 0
    while m < len(q) and q[m] == 0:
        m += 1
    return m


def check_minimality_k1():
    ex2 = ex_number_bruteforce(2, K3[0], K3[1], clique_r=3)
    ck("ex_2_K3_equals_1_by_enumeration", ex2 == 1, "ex(2;K_3) = %d" % ex2)
    worst = None
    bad = 0
    total = 0
    for n in range(1, 7):
        for E in all_graphs(n):
            total += 1
            p = charpoly(laplacian(n, E))
            # eps_1(G) <= 1  <=>  mu_1 <= |E| + 1, decided over the integers
            if count_gt(p, len(E) + 1) != 0:
                bad += 1
                worst = (n, E)
    ck("minimality_in_k_no_counterexample_at_k_equals_1", bad == 0,
       "eps_1(G) <= 1 = ex(2;K_3) for every one of the %d labelled graphs on <= 6 vertices%s"
       % (total, "" if worst is None else "; failure at %s" % (worst,)))


def check_minimality_small_order():
    """|V(G)| <= k+1 forces eps_k(G) = |E(G)| <= ex(k+1;H): no counterexample there."""
    ident_bad = 0
    bound_bad = 0
    tested = 0
    ex_tab = dict((m, ex_number_bruteforce(m, K3[0], K3[1], clique_r=3)) for m in range(2, 7))
    for n in range(2, 6):
        for E in all_graphs(n):
            for k in (n - 1, n):
                if k < 1:
                    continue
                tested += 1
                lo, hi = eps_k_bracket(n, E, k, 40)
                m = Fraction(len(E))
                if not (lo - Fraction(1, 10 ** 8) <= m <= hi + Fraction(1, 10 ** 8)):
                    ident_bad += 1
                # the conjecture's own bound at this k, for the H-free graphs:
                # no rigorous violation certificate, and |E(G)| <= ex(k+1;K_3)
                if not has_clique(n, E, 3):
                    if lo > Fraction(ex_tab[k + 1]) or len(E) > ex_tab[k + 1]:
                        bound_bad += 1
    ck("small_order_eps_k_equals_edge_count", ident_bad == 0,
       "eps_k(G) = |E(G)| verified on %d (graph, k) pairs with n-1 <= k <= n, n <= 5" % tested)
    ck("minimality_in_order_no_counterexample_when_nV_le_k_plus_1", bound_bad == 0,
       "eps_k(G) = |E(G)| <= ex(k+1;K_3), no violation certificate, for every triangle-free "
       "G on n <= 5 vertices with k in {n-1, n}")


def check_minimality_order_four():
    ex3 = ex_number_bruteforce(3, K3[0], K3[1], clique_r=3)
    small_bad = 0
    for n in (2, 3):
        for E in all_graphs(n):
            if has_clique(n, E, 3):
                continue
            lo, hi = eps_k_bracket(n, E, 2, 40)
            if hi > Fraction(ex3) + Fraction(1, 10 ** 8):
                small_bad += 1
    ck("minimality_in_order_none_below_four_vertices", small_bad == 0,
       "no triangle-free graph on 2 or 3 vertices violates the bound at k = 2")
    ex_bound = Fraction(ex3)
    violators = []
    for E in all_graphs(4):
        if has_clique(4, E, 3):
            continue
        lo, hi = eps_k_bracket(4, E, 2, 40)
        if lo > ex_bound:
            violators.append(E)
    all_are_p4 = all(
        len(E) == 3 and is_connected(4, E) and sorted(degrees(4, E)) == [1, 1, 2, 2]
        for E in violators
    )
    p4_in = any(sorted((min(a, b), max(a, b)) for a, b in E)
                == sorted((min(a, b), max(a, b)) for a, b in P4_EDGES) for E in violators)
    ck("order_four_census_P4_is_the_only_violator",
       p4_in and all_are_p4 and len(violators) == 12,
       "%d labelled triangle-free 4-vertex violators at k=2, all isomorphic to P_4 (12 = 4!/2)"
       % len(violators))


def check_turan_numbers():
    bad = []
    for r in (3, 4, 5):
        for n in range(2, 7):
            brute = ex_number_bruteforce(n, r, [], clique_r=r)
            if brute != turan_edges(n, r):
                bad.append((n, r, brute, turan_edges(n, r)))
    ck("turan_numbers_match_enumeration", not bad,
       "ex(n;K_r) = |E(T_{r-1}(n))| for all 2 <= n <= 6, r in {3,4,5}%s" % ("" if not bad else " %s" % bad))
    inc_bad = []
    for r in range(3, 8):
        for k in range(r - 1, 41):
            lhs = turan_edges(k + 1, r) - turan_edges(k, r)
            if lhs != k - k // (r - 1):
                inc_bad.append((r, k, lhs))
            if k + turan_edges(k, r) - turan_edges(k + 1, r) != k // (r - 1) or k // (r - 1) < 1:
                inc_bad.append((r, k, "d"))
    ck("turan_increment_and_deficit_formulas", not inc_bad,
       "ex(k+1;K_r)-ex(k;K_r) = k-floor(k/(r-1)) and k+ex(k;K_r)-ex(k+1;K_r) = floor(k/(r-1)) >= 1, "
       "r <= 7, k <= 40%s" % ("" if not inc_bad else " %s" % inc_bad))


def pendant_spectrum_prediction(n, edges, N):
    """The Laplacian charpoly of F^{<N>} as predicted by the pendant lemma:
    (t-1)^{k(N-1)} * sum_j a_j u^j (t-1)^{k-j},  u = t^2-(N+1)t,  a = charpoly(L(F))."""
    a = charpoly(laplacian(n, edges))
    k = n
    u = [0, -(N + 1), 1]
    tm1 = [-1, 1]
    acc = [0]
    for j, aj in enumerate(a):
        if aj == 0:
            continue
        term = pmul(ppow(u, j), ppow(tm1, k - j))
        term = [aj * c for c in term]
        if len(term) > len(acc):
            acc = list(acc) + [0] * (len(term) - len(acc))
        else:
            acc = list(acc)
        for i, c in enumerate(term):
            acc[i] += c
    return pmul(ppow(tm1, k * (N - 1)), acc)


def check_pendant_lemma():
    cases = [("K_2", 2, [(0, 1)], 2), ("P_3", 3, [(0, 1), (1, 2)], 4),
             ("K_3", 3, [(0, 1), (0, 2), (1, 2)], 6), ("P_4", 4, P4_EDGES, 3)]
    poly_bad, mult_bad, top_bad = [], [], []
    for (name, n, E, N) in cases:
        gn, gE = pendant_blowup(n, E, N)
        got = charpoly(laplacian(gn, gE))
        want = pendant_spectrum_prediction(n, E, N)
        want = want + [0] * (len(got) - len(want))
        if got != want[:len(got)] or any(want[len(got):]):
            poly_bad.append(name)
        if root_multiplicity(got, 1) != n * (N - 1):
            mult_bad.append((name, root_multiplicity(got, 1), n * (N - 1)))
        # the k largest eigenvalues are the beta_i: exactly k eigenvalues exceed 1,
        # and exactly k are >= N+1 (equality only for the mu_i = 0 of connected F)
        above_one = count_gt(got, 1)
        at_least_N1 = count_gt(got, N + 1) + root_multiplicity(got, N + 1)
        if above_one != n or at_least_N1 != n:
            top_bad.append((name, above_one, at_least_N1, n))
    ck("pendant_lemma_charpoly_identity", not poly_bad,
       "det(tI-L(F^<N>)) equals the lemma's product formula for F in K_2,P_3,K_3,P_4%s"
       % ("" if not poly_bad else " %s" % poly_bad))
    ck("pendant_lemma_eigenvalue_one_multiplicity", not mult_bad,
       "eigenvalue 1 has multiplicity exactly k(N-1) in every case%s"
       % ("" if not mult_bad else " %s" % mult_bad))
    ck("pendant_lemma_k_largest_are_the_beta_i", not top_bad,
       "exactly k Laplacian eigenvalues of F^<N> are > 1, and exactly k are >= N+1%s"
       % ("" if not top_bad else " %s" % top_bad))


def check_pendant_identity_and_bound():
    cases = [("K_2", 2, [(0, 1)], 2), ("P_3", 3, [(0, 1), (1, 2)], 4),
             ("K_3", 3, [(0, 1), (0, 2), (1, 2)], 6)]
    id_bad, bd_bad, al_bad = [], [], []
    for (name, k, E, N) in cases:
        e = len(E)
        gn, gE = pendant_blowup(k, E, N)
        p = charpoly(laplacian(gn, gE))
        top_lo, top_hi = top_k_sum_bracket(p, k, 40)
        eps_lo, eps_hi = top_lo - len(gE), top_hi - len(gE)
        # sum of the k SMALLEST eigenvalues: the alpha_i of the lemma
        a_lo = Fraction(0)
        a_hi = Fraction(0)
        alpha_his = []
        for i in range(gn - k + 1, gn + 1):
            lo, hi = eig_bracket(p, i, 40)
            a_lo += max(lo, Fraction(0))
            a_hi += hi
            alpha_his.append(hi)
        tol = Fraction(1, 10 ** 8)
        rhs_lo, rhs_hi = k + e - a_hi, k + e - a_lo
        if not (eps_lo - tol <= rhs_hi and rhs_lo - tol <= eps_hi):
            id_bad.append(name)
        if not (eps_lo > k + e - Fraction(2 * e, N + 1)):
            bd_bad.append(name)
        # a_lo is clamped at 0 by PSD-ness, so "a_lo >= 0" would be vacuous; the
        # failable content is the lemma's per-root bound alpha_i < 1.
        if not (a_hi < Fraction(2 * e, N + 1) and all(h < 1 for h in alpha_his)):
            al_bad.append((name, float(a_hi), float(max(alpha_his)),
                           float(Fraction(2 * e, N + 1))))
    ck("pendant_eps_identity_k_plus_e_minus_sum_alpha", not id_bad,
       "eps_k(F^<N>) from the top k eigenvalues agrees with k+e-sum(alpha_i) from the bottom k%s"
       % ("" if not id_bad else " %s" % id_bad))
    ck("pendant_strict_bound_eps_k_gt_k_plus_e_minus_2e_over_N1", not bd_bad,
       "strict inequality of the lemma holds in every case%s" % ("" if not bd_bad else " %s" % bd_bad))
    ck("pendant_sum_alpha_below_2e_over_N_plus_1", not al_bad,
       "0 <= sum(alpha_i) < 2e/(N+1) in every case%s" % ("" if not al_bad else " %s" % al_bad))


def complete_graph(k):
    return k, [(i, j) for i in range(k) for j in range(i + 1, k)]


def check_corollary_all_H():
    """H with delta(H) >= 2: F = K_{h-1}, N = 2e gives a connected H-free counterexample."""
    C4 = (4, [(0, 1), (1, 2), (2, 3), (0, 3)])
    cases = [("K_3", K3, 3), ("C_4", C4, 4)]
    rows = []
    free_bad, viol_bad, d_bad = [], [], []
    for (name, H, h) in cases:
        k = h - 1
        fn, fE = complete_graph(k)
        e = len(fE)
        N = 2 * e
        gn, gE = pendant_blowup(fn, fE, N)
        exh = ex_number_bruteforce(h, H[0], H[1])
        exk1 = ex_number_bruteforce(k + 1, H[0], H[1])
        if not (k + e == h * (h - 1) // 2 and exh <= h * (h - 1) // 2 - 1 and k + e - exk1 >= 1):
            d_bad.append(name)
        if contains_copy(gn, gE, H[0], H[1]) or not is_connected(gn, gE):
            free_bad.append(name)
        lo, hi = eps_k_bracket(gn, gE, k, 40)
        rows.append("%s: |V|=%d, eps_%d>=%.6f vs ex(%d;%s)=%d" % (name, gn, k, float(lo), k + 1, name, exk1))
        if not lo > Fraction(exk1):
            viol_bad.append(name)
        APPROX_INSTANCES.append(("K_%d^<2e> (Corollary 4), H=%s" % (k, name), k, hi,
                                 ex_number_bruteforce(k, H[0], H[1])))
    print("pendant amplification instances: " + "; ".join(rows))
    ck("amplified_graphs_are_connected_and_H_free", not free_bad,
       "K_2^<2> is K_3-free, K_3^<6> is C_4-free, both connected%s"
       % ("" if not free_bad else " %s" % free_bad))
    ck("amplification_deficit_d_at_least_one", not d_bad,
       "k+e = C(h,2) and ex(h;H) <= C(h,2)-1, so d >= 1 for H = K_3, C_4%s"
       % ("" if not d_bad else " %s" % d_bad))
    ck("corollary_every_H_with_min_degree_two_has_counterexample", not viol_bad,
       "Corollary 4: eps_k(K_{h-1}^<2e>) > ex(k+1;H) for H = K_3, C_4%s"
       % ("" if not viol_bad else " %s" % viol_bad))


def check_corollary_complete():
    """r >= 3, k >= r-1: F = T_{r-1}(k), N = 2q gives a K_r-free counterexample."""
    viol_bad, free_bad, excess_bad = [], [], []
    rows = []
    for (r, k) in [(3, 3), (3, 4), (4, 3)]:
        fn, fE = turan_graph(k, r)
        q = len(fE)
        exk = ex_number_bruteforce(k, r, [], clique_r=r)
        exk1 = ex_number_bruteforce(k + 1, r, [], clique_r=r)
        gn, gE = pendant_blowup(fn, fE, 2 * q)
        if has_clique(gn, gE, r) or not is_connected(gn, gE):
            free_bad.append((r, k))
        lo, hi = eps_k_bracket(gn, gE, k, 40)
        rows.append("r=%d,k=%d: |V|=%d, eps_%d in (%.6f,%.6f], ex(%d;K_%d)=%d"
                    % (r, k, gn, k, float(lo), float(hi), k + 1, r, exk1))
        if not (q == exk and lo > Fraction(exk1)):
            viol_bad.append((r, k))
        APPROX_INSTANCES.append(("T_%d(%d)^<2q> (Corollary 5), H=K_%d" % (r - 1, k, r), k, hi, exk))
        fl = k // (r - 1)
        if not (Fraction(0) < lo - exk1 and hi - exk1 < Fraction(fl)):
            excess_bad.append((r, k, float(lo - exk1), fl))
    print("Turan-based instances: " + "; ".join(rows))
    ck("Kr_free_instances_are_connected_and_Kr_free", not free_bad,
       "T_{r-1}(k)^<2q> is connected and K_r-free for (r,k) in {(3,3),(3,4),(4,3)}%s"
       % ("" if not free_bad else " %s" % free_bad))
    ck("corollary_Kr_free_counterexample_at_parameter_k", not viol_bad,
       "Corollary 5: |E(T_{r-1}(k))| = ex(k;K_r) and eps_k > ex(k+1;K_r) in every instance%s"
       % ("" if not viol_bad else " %s" % viol_bad))
    ck("excess_positive_but_below_floor_k_over_r_minus_1", not excess_bad,
       "Remark 6 (excess): 0 < eps_k(G) - ex(k+1;K_r) < floor(k/(r-1)) in every instance%s"
       % ("" if not excess_bad else " %s" % excess_bad))


def check_remark_approximate_bound():
    """Remark 6: the counterexamples are consistent with the approximate bound
    eps_k(G) <= ex(k;H) + (4k-2)sqrt(k) that the paper records from Section 7 of
    the cited preprint.  Decided exactly: with D = eps_hi - ex(k;H) rational and
    c = 4k-2, the assertion D <= c*sqrt(k) is equivalent to D <= 0 or
    D^2 <= c^2 k, so no float and no sqrt approximation enters the decision.
    eps_hi is an upper bound for eps_k(G), so a PASS is a rigorous verification."""
    ck("remark_six_all_counterexample_instances_registered",
       len(APPROX_INSTANCES) == EXPECTED_INSTANCES,
       "%d of %d expected instances registered upstream (P_4, 2 from Corollary 4, "
       "3 from Corollary 5)" % (len(APPROX_INSTANCES), EXPECTED_INSTANCES))
    bad = []
    rows = []
    for (label, k, eps_hi, exk) in APPROX_INSTANCES:
        c = 4 * k - 2
        D = Fraction(eps_hi) - Fraction(exk)
        ok = (D <= 0) or (D * D <= Fraction(c * c * k))
        if not ok:
            bad.append(label)
        rows.append("%s: eps_%d <= %.6f, ex(%d;H) = %d, allowance %d*sqrt(%d)"
                    % (label, k, float(eps_hi), k, exk, c, k))
    print("Remark 6 approximate-bound instances: " + "; ".join(rows))
    ck("remark_six_instances_consistent_with_approximate_bound",
       (not bad) and len(APPROX_INSTANCES) == EXPECTED_INSTANCES,
       "eps_k(G) <= ex(k;H) + (4k-2)sqrt(k) for all %d instances, decided by exact "
       "squaring%s" % (len(APPROX_INSTANCES), "" if not bad else "; violations %s" % bad))


def check_remark_isolated_vertex():
    """H = K_2 + K_1, G = K_2, k = 2: a two-vertex failure of the literal statement."""
    H = (3, [(0, 1)])
    exH3 = ex_number_bruteforce(3, H[0], H[1])
    gn, gE = complete_graph(2)
    hfree = not contains_copy(gn, gE, H[0], H[1])
    lo, hi = eps_k_bracket(gn, gE, 2, 40)
    ck("ex_3_K2_plus_K1_equals_0", exH3 == 0,
       "every 3-vertex graph with an edge contains K_2 + K_1, so ex(3;H) = %d" % exH3)
    ck("remark_K2_violates_at_k_two_with_isolated_vertex_H",
       hfree and gn >= 2 and lo > Fraction(exH3) and hi < Fraction(exH3) + 2,
       "Remark 8: K_2 is (K_2+K_1)-free and eps_2(K_2) in (%.6f, %.6f] > 0"
       % (float(lo), float(hi)))


def main():
    check_object()
    check_charpoly()
    check_spectrum()
    check_eps2_value()
    check_bound_and_refutation()
    check_minimality_k1()
    check_minimality_small_order()
    check_minimality_order_four()
    check_turan_numbers()
    check_pendant_lemma()
    check_pendant_identity_and_bound()
    check_corollary_all_H()
    check_corollary_complete()
    check_remark_approximate_bound()
    check_remark_isolated_vertex()
    print("NOT RE-RUN: no exhaustive census beyond the ranges printed above -- the k=1 census "
          "covers all labelled graphs on at most 6 vertices, the order censuses all graphs on at "
          "most 5 vertices (all triangle-free graphs on exactly 4 vertices at k=2), Turan numbers "
          "are enumerated for n <= 6 and the increment formula is checked for r <= 7, k <= 40; the "
          "general theorems for arbitrary H and k are verified only on the listed instances. "
          "Remark 6 is only partly evaluated: its consistency assertion eps_k <= ex(k;H)+(4k-2)"
          "sqrt(k) is tested on the six instances listed above and not for all r and k, and its "
          "asymptotic comparison (excess O(k) against ex(k+1;K_r) = Theta(k^2)) is verified only "
          "through the per-instance excess bound, not in general. NOT CHECKED AGAINST ANY SOURCE: "
          "the wording and hypotheses of Conjecture 7.1, and the approximate bound, are "
          "transcribed from the paper's own quotation of arXiv:2601.17575v1; this program has no "
          "access to that preprint, so a transcription error in the conjecture -- for instance a "
          "dropped hypothesis on H, which Remark 8 shows would matter -- would not be detected "
          "here.")
    finish()


def finish():
    bad = [n for n, ok in CHECKS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), len(CHECKS)))
        raise SystemExit(1)
    print("VERDICT: ALL %d CHECKS PASS" % len(CHECKS))
    raise SystemExit(0)


if __name__ == "__main__":
    main()

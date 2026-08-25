#!/usr/bin/env python3
"""Verification of a counterexample to a Cartesian-product conjecture for
node resistance curvature.

Setting.  For a finite connected simple graph G with unit conductances, let
L_G be its Laplacian, L_G^+ the Moore-Penrose pseudoinverse, and
omega_G(u,v) = (e_u - e_v)^T L_G^+ (e_u - e_v) the effective resistance.
The node resistance curvature at v is

    p_v(G) = 1 - (1/2) * sum_{u ~ v} omega_G(u,v).

The conjecture acted on:  if p_i(G_1) <= 0 and p_j(G_2) <= 0, then
p_{(i,j)}(G_1 box G_2) < 0.  That statement is TRANSCRIBED BY HAND from the
source it is attributed to (Conjecture 4 of Dawkins et al.,
arXiv:2403.01037v1); this program does not read the source, so the
transcription is an unverified input on the same footing as the edge lists.
It is printed at the head of the transcript and flagged again in the closing
NOT RE-RUN paragraphs so that a referee can compare it against the source by
eye.

VALUES TAKEN FROM THE PAPER (inputs, nothing else is assumed):
  * G_1: V = {0,...,4}, E = {01, 02, 13, 34};
  * G_2: V = {0,...,7}, E = {01, 02, 03, 04, 15, 56, 57};
  * the marked vertices i = j = 1;
  * the asserted values, which are re-derived here and compared:
        p_i(G_1) = p_j(G_2) = 0,
        tau(H) = 723963572556125 for H = G_1 box G_2,
        tau(H/e) = 359244778688125, 359244778688125,
                   354728028191625, 363877814994375
            for e = {x,(0,1)}, {x,(3,1)}, {x,(1,0)}, {x,(1,5)}, x = (1,1),
        2*tau(H) = 1447927145112250,
        2*tau(H) - sum_e tau(H/e) = 10831744550000,
        p_x(H) = 32552200/4351396379 > 0.

DERIVED HERE (computed from the two edge lists alone, in exact integer or
rational arithmetic; no floating point enters any decision):
  * both graphs are simple, connected and acyclic, with the stated vertex
    and edge counts, and the marked vertex has degree 2 in each;
  * the effective resistances and hence p_1(G_1), p_1(G_2), from the
    pseudoinverse definition, so the hypotheses of the conjecture hold;
  * the Cartesian product H, its order, size, connectivity and the
    neighbour set of x = (1,1);
  * tau(H) by the matrix-tree theorem (fraction-free integer determinant)
    and, independently, from the characteristic polynomials of the two
    factor Laplacians via a resultant;
  * the four contracted counts tau(H/e), each from the weighted Laplacian
    of the contracted multigraph;
  * p_x(H) by three independent routes (all-minors matrix-tree, spanning
    tree contraction ratios, and an exact rational solve against the
    grounded Laplacian), the reduced fraction, and its strict positivity,
    which contradicts the conjectured conclusion;
  * a census: all rooted trees on at most 8 vertices are generated up to
    isomorphism and the root curvature law p_root = 1 - deg(root)/2 is
    checked on each; then every pair of rooted trees with root degree 2,
    each factor on at most 8 vertices and at most 13 vertices in total, is
    tested, cross-validating two routes.
"""
import sys
from fractions import Fraction
from itertools import combinations_with_replacement

CHECKS = []

def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)

def det_int(M):
    """Determinant of an integer matrix by the Bareiss fraction-free
    algorithm.  Exact: every intermediate value is an integer."""
    A = [row[:] for row in M]
    n = len(A)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            piv = None
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    piv = r
                    break
            if piv is None:
                return 0
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = A[i][j] * A[k][k] - A[i][k] * A[k][j]
                q, r = divmod(num, prev)
                if r != 0:
                    raise ArithmeticError("Bareiss division was not exact")
                A[i][j] = q
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def laplacian(n, edges):
    """Weighted Laplacian of a multigraph on {0,...,n-1}.  Each entry of
    `edges` is (u, v) or (u, v, w); parallel edges add up."""
    L = [[0] * n for _ in range(n)]
    for e in edges:
        u, v = e[0], e[1]
        w = e[2] if len(e) > 2 else 1
        if u == v:
            continue
        L[u][u] += w
        L[v][v] += w
        L[u][v] -= w
        L[v][u] -= w
    return L


def minor_del(L, drop):
    """L with the rows and columns in `drop` deleted."""
    keep = [i for i in range(len(L)) if i not in drop]
    return [[L[i][j] for j in keep] for i in keep]


def tau(n, edges):
    """Number of spanning trees (weighted, with multiplicity) by the
    matrix-tree theorem: any principal cofactor of the Laplacian."""
    L = laplacian(n, edges)
    return det_int(minor_del(L, {0}))


def resistance_minors(n, edges, u, v):
    """omega(u,v) as a ratio of Laplacian minors (all-minors matrix-tree
    theorem): the number of 2-forests separating u from v, over tau."""
    L = laplacian(n, edges)
    num = det_int(minor_del(L, {u, v}))
    den = det_int(minor_del(L, {u}))
    return Fraction(num, den)


def contract(n, edges, u, v):
    """Contract the edge {u,v}: identify v with u, delete resulting self
    loops, keep parallel edges with multiplicity.  Returns (n-1, edges)."""
    rel = {}
    nxt = 0
    for w in range(n):
        if w == v:
            continue
        rel[w] = nxt
        nxt += 1
    rel[v] = rel[u]
    out = []
    for e in edges:
        a, b = rel[e[0]], rel[e[1]]
        if a == b:
            continue
        out.append((a, b, e[2] if len(e) > 2 else 1))
    return nxt, out


def resistance_contraction(n, edges, u, v):
    """omega(u,v) for an edge {u,v} as tau(G/e)/tau(G), i.e. the
    probability that e lies in a uniformly random spanning tree."""
    m, ce = contract(n, edges, u, v)
    return Fraction(tau(m, ce), tau(n, edges))


def solve_exact(A, b):
    """Solve A y = b over the rationals by Gaussian elimination with
    Fraction entries.  Returns None if A is singular."""
    n = len(A)
    M = [[Fraction(A[i][j]) for j in range(n)] + [Fraction(b[i])]
         for i in range(n)]
    for c in range(n):
        piv = None
        for r in range(c, n):
            if M[r][c] != 0:
                piv = r
                break
        if piv is None:
            return None
        M[c], M[piv] = M[piv], M[c]
        inv = Fraction(1, 1) / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][j] - f * M[c][j] for j in range(n + 1)]
    return [M[i][n] for i in range(n)]


def resistance_pinv(n, edges, u, v):
    """omega(u,v) straight from the definition: solve L x = e_u - e_v on
    the subspace orthogonal to the constants (ground vertex n-1), then
    omega = x_u - x_v.  Independent of any matrix-tree identity."""
    L = laplacian(n, edges)
    r = n - 1
    idx = [i for i in range(n) if i != r]
    A = [[L[i][j] for j in idx] for i in idx]
    b = [(1 if i == u else 0) - (1 if i == v else 0) for i in idx]
    y = solve_exact(A, b)
    if y is None:
        raise ArithmeticError("grounded Laplacian is singular")
    x = {r: Fraction(0)}
    for k, i in enumerate(idx):
        x[i] = y[k]
    return x[u] - x[v]


def neighbours(n, edges, v):
    out = set()
    for e in edges:
        if e[0] == v:
            out.add(e[1])
        elif e[1] == v:
            out.add(e[0])
    return sorted(out)


def curvature(n, edges, v, resistance):
    s = sum(resistance(n, edges, v, u) for u in neighbours(n, edges, v))
    return Fraction(1) - Fraction(1, 2) * s


def is_simple(n, edges):
    seen = set()
    for e in edges:
        u, v = e[0], e[1]
        if u == v or not (0 <= u < n and 0 <= v < n):
            return False
        key = (min(u, v), max(u, v))
        if key in seen:
            return False
        seen.add(key)
    return True


def is_connected(n, edges):
    adj = dict((i, []) for i in range(n))
    for e in edges:
        adj[e[0]].append(e[1])
        adj[e[1]].append(e[0])
    seen = {0}
    stack = [0]
    while stack:
        w = stack.pop()
        for z in adj[w]:
            if z not in seen:
                seen.add(z)
                stack.append(z)
    return len(seen) == n


def is_tree(n, edges):
    return is_simple(n, edges) and is_connected(n, edges) and \
        len(edges) == n - 1


def cartesian(n1, e1, n2, e2):
    """G_1 box G_2 with (a,b) indexed as a*n2 + b."""
    def ix(a, b):
        return a * n2 + b
    edges = []
    for (a, c) in [(e[0], e[1]) for e in e1]:
        for b in range(n2):
            edges.append((ix(a, b), ix(c, b)))
    for a in range(n1):
        for (b, d) in [(e[0], e[1]) for e in e2]:
            edges.append((ix(a, b), ix(a, d)))
    return n1 * n2, edges, ix


def charpoly(M):
    """Coefficients (ascending) of det(x*I - M) for an integer matrix M,
    by exact evaluation at n+1 integer points and Lagrange interpolation."""
    n = len(M)
    pts = list(range(n + 1))
    vals = []
    for t in pts:
        A = [[(t if i == j else 0) - M[i][j] for j in range(n)]
             for i in range(n)]
        vals.append(det_int(A))
    coef = [Fraction(0)] * (n + 1)
    for k, t in enumerate(pts):
        basis = [Fraction(1)]
        den = Fraction(1)
        for s in pts:
            if s == t:
                continue
            den *= (t - s)
            new = [Fraction(0)] * (len(basis) + 1)
            for i, c in enumerate(basis):
                new[i] += c * (-s)
                new[i + 1] += c
            basis = new
        scale = Fraction(vals[k]) / den
        for i, c in enumerate(basis):
            coef[i] += c * scale
    out = []
    for c in coef:
        if c.denominator != 1:
            raise ArithmeticError("interpolated coefficient not an integer")
        out.append(int(c))
    return out


def poly_shift_signs(coef):
    """Given f(x) = sum c_k x^k of degree a, return h(y) = prod (lam_i + y)
    where the lam_i are the roots of the monic f, i.e. (-1)^a f(-y)."""
    a = len(coef) - 1
    return [c * (-1) ** ((a + k) % 2) for k, c in enumerate(coef)]


def companion(coef):
    """Companion matrix of the monic polynomial with ascending coefs."""
    d = len(coef) - 1
    C = [[0] * d for _ in range(d)]
    for i in range(1, d):
        C[i][i - 1] = 1
    for i in range(d):
        C[i][d - 1] = -coef[i]
    return C


def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]


def poly_of_matrix(coef, M):
    """Evaluate the polynomial with ascending coefficients `coef` at the
    square integer matrix M."""
    d = len(M)
    acc = [[0] * d for _ in range(d)]
    power = [[1 if i == j else 0 for j in range(d)] for i in range(d)]
    for c in coef:
        if c:
            for i in range(d):
                for j in range(d):
                    acc[i][j] += c * power[i][j]
        power = mat_mul(power, M)
    return acc


def tau_product_via_spectra(n1, e1, n2, e2):
    """tau(G_1 box G_2) from the factor Laplacian characteristic
    polynomials.  With L(G) having eigenvalues 0 = lam_0 < lam_1 <= ...,

        n1*n2*tau(H) = (prod_{i>0} lam_i)(prod_{j>0} mu_j)
                       * prod_{i>0, j>0} (lam_i + mu_j),

    and the last product is det(h(C)), where h(y) = prod_{i>0}(lam_i + y)
    and C is the companion matrix of the deflated char. polynomial of
    L(G_2).  Uses only integer polynomial and matrix arithmetic."""
    f = charpoly(laplacian(n1, e1))
    g = charpoly(laplacian(n2, e2))
    if f[0] != 0 or g[0] != 0:
        raise ArithmeticError("Laplacian char. polynomial lacks the root 0")
    fd, gd = f[1:], g[1:]
    prod_f = (-1) ** ((n1 - 1) % 2) * fd[0]
    prod_g = (-1) ** ((n2 - 1) % 2) * gd[0]
    h = poly_shift_signs(fd)
    cross = det_int(poly_of_matrix(h, companion(gd)))
    total = prod_f * prod_g * cross
    q, r = divmod(total, n1 * n2)
    if r != 0:
        raise ArithmeticError("spectral product not divisible by n1*n2")
    return q


def partitions_into(total):
    """All multisets of positive integers summing to `total`, each given as
    a non-increasing tuple."""
    if total == 0:
        yield ()
        return
    for first in range(total, 0, -1):
        for rest in partitions_into(total - first):
            if not rest or rest[0] <= first:
                yield (first,) + rest


_ROOTED = {1: [()]}


def rooted_trees(n):
    """All rooted trees on n vertices up to rooted isomorphism, encoded
    canonically: a tree is the sorted tuple of its subtrees' encodings."""
    if n in _ROOTED:
        return _ROOTED[n]
    found = set()
    for part in partitions_into(n - 1):
        sizes = sorted(set(part))
        choices = []
        for s in sizes:
            mult = part.count(s)
            choices.append(list(combinations_with_replacement(
                rooted_trees(s), mult)))

        def walk(k, acc):
            if k == len(choices):
                found.add(tuple(sorted(acc)))
                return
            for pick in choices[k]:
                walk(k + 1, acc + list(pick))
        walk(0, [])
    _ROOTED[n] = sorted(found)
    return _ROOTED[n]


def canon_size(form):
    return 1 + sum(canon_size(c) for c in form)


def tree_from_canon(form):
    """Realise a canonical rooted tree as (n, edges) with root 0."""
    edges = []
    counter = [0]

    def build(f):
        me = counter[0]
        counter[0] += 1
        for child in f:
            kid = build(child)
            edges.append((me, kid))
        return me
    build(form)
    return counter[0], edges


def canon_of(n, edges, root):
    """Canonical encoding of the rooted tree (edges, root), in the same
    format the generator produces."""
    adj = dict((i, []) for i in range(n))
    for e in edges:
        adj[e[0]].append(e[1])
        adj[e[1]].append(e[0])

    def rec(v, parent):
        return tuple(sorted(rec(w, v) for w in adj[v] if w != parent))
    return rec(root, None)


# ---------------------------------------------------------------------------
# Values taken from the paper.
# ---------------------------------------------------------------------------
N1 = 5
E1 = [(0, 1), (0, 2), (1, 3), (3, 4)]
N2 = 8
E2 = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (5, 6), (5, 7)]
MARK1 = 1
MARK2 = 1
PAPER_TAU_H = 723963572556125
PAPER_COFACTORS = {(0, 1): 359244778688125,
                   (3, 1): 359244778688125,
                   (1, 0): 354728028191625,
                   (1, 5): 363877814994375}
PAPER_TWO_TAU = 1447927145112250
PAPER_NUMERATOR = 10831744550000
PAPER_CURVATURE = Fraction(32552200, 4351396379)
PAPER_NEIGHBOURS = [(0, 1), (1, 0), (1, 5), (3, 1)]

# ---------------------------------------------------------------------------
# The conjecture being refuted.  This is a HAND TRANSCRIPTION, not a
# computation: no part of this program fetches or parses the source, so the
# attribution (which numbered statement, in which version of which preprint)
# and the exact logical form below are inputs, not results.  The two
# predicates are the only place the logical form is used, so the statement
# printed in the transcript is exactly the statement the verdict is about.
# ---------------------------------------------------------------------------
CONJECTURE_SOURCE = ("Conjecture 4 of Dawkins et al., "
                     "arXiv:2403.01037v1 (hand transcription)")
CONJECTURE_TEXT = ("for graphs G1, G2 and vertices i in V(G1), j in V(G2): "
                   "if p_i(G1) <= 0 and p_j(G2) <= 0, "
                   "then p_(i,j)(G1 box G2) < 0")


def conj_hypothesis(p_i, p_j):
    """The transcribed hypothesis: both factor curvatures NON-POSITIVE (the
    inequalities are not strict, so equality at 0 is admitted)."""
    return p_i <= 0 and p_j <= 0


def conj_conclusion(p_x):
    """The transcribed conclusion: the product curvature is STRICTLY
    negative."""
    return p_x < 0


def check_factor(tag, n, edges, mark, expect_n, expect_m):
    ok = (n == expect_n and len(edges) == expect_m and is_tree(n, edges))
    degs = [len(neighbours(n, edges, v)) for v in range(n)]
    ck("factor_%s_is_a_tree" % tag, ok,
       "n=%d, m=%d, edges=%s, degrees=%s" %
       (n, len(edges), sorted((min(e), max(e)) for e in edges), degs))
    ck("factor_%s_marked_vertex_degree_2" % tag, degs[mark] == 2,
       "deg(%d) = %d, neighbours %s" %
       (mark, degs[mark], neighbours(n, edges, mark)))
    return ok


def check_factor_curvature(tag, n, edges, mark):
    p_def = curvature(n, edges, mark, resistance_pinv)
    p_mt = curvature(n, edges, mark, resistance_minors)
    res = [resistance_pinv(n, edges, mark, u)
           for u in neighbours(n, edges, mark)]
    ck("factor_%s_root_curvature_is_zero" % tag,
       p_def == 0 and p_mt == p_def,
       "resistances %s, p_%d = %s (pseudoinverse) = %s (matrix-tree)" %
       ([str(r) for r in res], mark, p_def, p_mt))
    ck("factor_%s_spanning_tree_count_is_one" % tag, tau(n, edges) == 1,
       "tau = %d" % tau(n, edges))
    return p_def


def check_product_structure(n, edges, ix):
    expect_n = N1 * N2
    expect_m = len(E1) * N2 + N1 * len(E2)
    ok = (n == expect_n and len(edges) == expect_m and
          is_simple(n, edges) and is_connected(n, edges))
    ck("product_is_a_connected_simple_graph_on_40_vertices", ok,
       "|V| = %d (expected %d = %d*%d), |E| = %d (expected %d = %d*%d + "
       "%d*%d)" %
       (n, expect_n, N1, N2, len(edges), expect_m,
        len(E1), N2, N1, len(E2)))
    x = ix(MARK1, MARK2)
    nb = [(v // N2, v % N2) for v in neighbours(n, edges, x)]
    ck("product_neighbourhood_of_marked_vertex",
       sorted(nb) == sorted(PAPER_NEIGHBOURS) and len(nb) == 4,
       "x = (%d,%d); neighbours %s" % (MARK1, MARK2, sorted(nb)))
    return x


def check_tau(n, edges):
    t_mt = tau(n, edges)
    ck("product_spanning_tree_count_matrix_tree", t_mt == PAPER_TAU_H,
       "tau(H) = %d, paper %d" % (t_mt, PAPER_TAU_H))
    t_sp = tau_product_via_spectra(N1, E1, N2, E2)
    ck("product_spanning_tree_count_spectral_route", t_sp == t_mt == PAPER_TAU_H,
       "n1*n2*tau(H) from the factor characteristic polynomials gives "
       "tau(H) = %d" % t_sp)
    return t_mt


def check_cofactors(n, edges, ix, x):
    total = 0
    allok = True
    detail = []
    inc = set(frozenset((e[0], e[1])) for e in edges)
    for (a, b) in sorted(PAPER_COFACTORS):
        u = ix(a, b)
        if u == x or frozenset((x, u)) not in inc:
            allok = False
            detail.append("(%d,%d): not adjacent to x" % (a, b))
            continue
        m, ce = contract(n, edges, x, u)
        got = tau(m, ce)
        want = PAPER_COFACTORS[(a, b)]
        allok = allok and got == want and m == n - 1
        total += got
        detail.append("(%d,%d):%d" % (a, b, got))
    ck("product_contracted_spanning_tree_counts", allok,
       "tau(H/e) for e = {x,v}: " + ", ".join(detail))
    return total


def check_curvature_at_x(n, edges, x, t_h, cof_sum):
    p_a = curvature(n, edges, x, resistance_minors)
    p_b = curvature(n, edges, x, resistance_contraction)
    p_c = curvature(n, edges, x, resistance_pinv)
    ck("product_curvature_three_routes_agree",
       p_a == p_b == p_c,
       "all-minors %s, contraction ratios %s, rational pseudoinverse %s"
       % (p_a, p_b, p_c))
    ck("product_curvature_equals_paper_value", p_c == PAPER_CURVATURE,
       "p_x(H) = %s, paper %s" % (p_c, PAPER_CURVATURE))
    frac = Fraction(2 * t_h - cof_sum, 2 * t_h)
    ck("paper_intermediate_arithmetic",
       2 * t_h == PAPER_TWO_TAU and 2 * t_h - cof_sum == PAPER_NUMERATOR
       and frac == PAPER_CURVATURE
       and (frac.numerator, frac.denominator) ==
       (PAPER_CURVATURE.numerator, PAPER_CURVATURE.denominator),
       "2*tau = %d, 2*tau - sum tau(H/e) = %d, reduced %s"
       % (2 * t_h, 2 * t_h - cof_sum, frac))
    return p_c


def check_transcribed_predicates():
    """The statement printed at the head of the transcript claims a
    NON-STRICT inequality in each hypothesis and a STRICT inequality in the
    conclusion.  The two predicates above are the only encoding of that
    logical form, and they are a separate object from the printed text, so
    probe their truth table on the boundary cases.  This does NOT check the
    transcription against the source -- nothing in this program can -- but it
    stops the printed statement and the code that decides the verdict from
    disagreeing about which inequalities are strict, which is precisely
    where the refutation lives (the exhibited pair sits at p_i = p_j = 0)."""
    zero, pos, neg = Fraction(0), Fraction(1, 7), Fraction(-1, 7)
    hyp_table = [((zero, zero), True), ((neg, zero), True),
                 ((zero, neg), True), ((neg, neg), True),
                 ((pos, zero), False), ((zero, pos), False),
                 ((pos, pos), False)]
    con_table = [(neg, True), (zero, False), (pos, False)]
    hyp_got = [(a, b, bool(conj_hypothesis(a, b))) for (a, b), _w in hyp_table]
    con_got = [(p, bool(conj_conclusion(p))) for p, _w in con_table]
    ok = (all(g[2] == w for g, (_ab, w) in zip(hyp_got, hyp_table)) and
          all(g[1] == w for g, (_p, w) in zip(con_got, con_table)))
    ck("transcribed_conjecture_predicates_are_the_printed_inequalities", ok,
       "hypothesis(p_i,p_j) -> %s (non-strict: admits 0, rejects any "
       "positive); conclusion(p_x) -> %s (strict: rejects 0)"
       % (", ".join("(%s,%s):%s" % (a, b, g) for a, b, g in hyp_got),
          ", ".join("%s:%s" % (p, g) for p, g in con_got)))


def check_refutation(p_i, p_j, p_x):
    check_transcribed_predicates()
    hyp = conj_hypothesis(p_i, p_j)
    ck("conjecture_hypotheses_are_satisfied", hyp,
       "p_i(G1) = %s <= 0 and p_j(G2) = %s <= 0" % (p_i, p_j))
    sign = ">" if p_x > 0 else ("=" if p_x == 0 else "<")
    reading = ("so the statement is refuted" if hyp and p_x > 0
               else "which does not refute the statement")
    ck("conjecture_conclusion_fails_at_the_marked_vertex",
       hyp and not conj_conclusion(p_x) and p_x > 0,
       "conclusion demands p_x(H) < 0; computed p_x(H) = %s %s 0, %s"
       % (p_x, sign, reading))


ROOTED_COUNTS = [1, 1, 2, 4, 9, 20, 48, 115]   # rooted trees up to 8 vertices


def check_enumeration():
    got = [len(rooted_trees(k)) for k in range(1, 9)]
    ck("rooted_tree_enumeration_counts", got == ROOTED_COUNTS,
       "rooted trees on 1..8 vertices, up to isomorphism: %s" % got)
    c1 = canon_of(N1, E1, MARK1)
    c2 = canon_of(N2, E2, MARK2)
    ok = (canon_size(c1) == N1 and canon_size(c2) == N2 and
          c1 in rooted_trees(N1) and c2 in rooted_trees(N2))
    ck("exhibited_factors_appear_in_enumeration", ok,
       "G1 rooted at %d = %s, G2 rooted at %d = %s" % (MARK1, c1, MARK2, c2))
    return c1, c2


def check_curvature_law():
    bad = []
    total = 0
    maxn = 0
    for n in range(1, 9):
        maxn = max(maxn, n)
        for form in rooted_trees(n):
            nn, ee = tree_from_canon(form)
            deg = len(form)
            if nn == 1:
                p = Fraction(1)
            else:
                p = curvature(nn, ee, 0, resistance_pinv)
            total += 1
            if p != Fraction(1) - Fraction(deg, 2):
                bad.append(form)
    ck("root_curvature_law_on_all_rooted_trees_upto_8",
       not bad and total == sum(ROOTED_COUNTS),
       "p_root = 1 - deg(root)/2 verified on all %d rooted trees with at "
       "most %d vertices (expected %d trees); %d exceptions, so p_root = 0 "
       "exactly for root degree 2"
       % (total, maxn, sum(ROOTED_COUNTS), len(bad)))


CENSUS_MAX = 8          # largest factor searched
CENSUS_SUM = 13         # largest total order searched
CENSUS_PAIRS = 427      # pairs the search must cover
CENSUS_VIOLATORS = 6    # pairs found with p_x > 0


def run_census():
    """All unordered pairs of rooted trees with root curvature zero, both
    factors on at most CENSUS_MAX vertices and CENSUS_SUM vertices in
    total.  Returns (pairs, disagreements, hypothesis failures, violators)."""
    pool = {}
    for k in range(3, CENSUS_MAX + 1):
        pool[k] = [f for f in rooted_trees(k) if len(f) == 2]
    pairs = 0
    disagree = []
    hypfail = []
    violators = []
    for a in range(3, CENSUS_MAX + 1):
        for b in range(a, CENSUS_MAX + 1):
            if a + b > CENSUS_SUM:
                continue
            for fa in pool[a]:
                n1, e1 = tree_from_canon(fa)
                pa = curvature(n1, e1, 0, resistance_minors)
                for fb in pool[b]:
                    if a == b and fb < fa:
                        continue
                    n2, e2 = tree_from_canon(fb)
                    pb = curvature(n2, e2, 0, resistance_minors)
                    if pa != 0 or pb != 0:
                        hypfail.append((fa, fb))
                    n, edges, ix = cartesian(n1, e1, n2, e2)
                    x = ix(0, 0)
                    p = curvature(n, edges, x, resistance_minors)
                    q = curvature(n, edges, x, resistance_contraction)
                    if p != q:
                        disagree.append((fa, fb))
                    pairs += 1
                    if p > 0:
                        violators.append((a, b, fa, fb, p))
    return pairs, disagree, hypfail, violators


def check_census(c1, c2):
    pairs, disagree, hypfail, violators = run_census()
    ck("census_two_routes_agree_on_every_pair",
       not disagree and pairs == CENSUS_PAIRS,
       "%d pairs searched (all factors <= %d vertices, total <= %d); "
       "all-minors and contraction routes disagree on %d"
       % (pairs, CENSUS_MAX, CENSUS_SUM, len(disagree)))
    ck("census_every_pair_satisfies_the_hypothesis", not hypfail,
       "root curvature recomputed as 0 for both factors of all %d pairs; "
       "%d failures" % (pairs, len(hypfail)))
    found = [(a, b, p) for (a, b, fa, fb, p) in violators
             if {(a, fa), (b, fb)} == {(N1, c1), (N2, c2)}]
    ck("census_exhibited_pair_is_a_violator",
       len(found) == 1 and found[0][2] == PAPER_CURVATURE
       and len(violators) == CENSUS_VIOLATORS,
       "%d violating pairs in range; the exhibited pair occurs with "
       "p_x = %s" % (len(violators),
                     found[0][2] if found else "not found"))
    small = [(a, b, str(p)) for (a, b, fa, fb, p) in violators
             if a + b < CENSUS_SUM]
    ck("census_no_violator_of_smaller_total_order", not small,
       "every violator found has %d vertices in total, as does the "
       "exhibited pair (%d + %d); factor orders %s; violators of smaller "
       "total order: %s"
       % (CENSUS_SUM, N1, N2,
          sorted(set((a, b) for (a, b, _f, _g, _p) in violators)), small))
    print("NOT RE-RUN: pairs of rooted trees with more than %d vertices in "
          "total, or with a factor on more than %d vertices, were not "
          "enumerated; nor were factors that are not trees; nor factors of "
          "root curvature strictly negative (root degree at least 3), which "
          "also satisfy the hypothesis, so the census above is a statement "
          "about curvature-zero pairs only; the paper makes no minimality "
          "claim." % (CENSUS_SUM, CENSUS_MAX))


def main():
    print("CONJECTURE ACTED ON (%s): %s"
          % (CONJECTURE_SOURCE, CONJECTURE_TEXT))
    print("G1: V = {0..%d}, E = %s" % (N1 - 1, E1))
    print("G2: V = {0..%d}, E = %s" % (N2 - 1, E2))
    print("marked vertices: i = %d in G1, j = %d in G2" % (MARK1, MARK2))
    check_factor("G1", N1, E1, MARK1, 5, 4)
    p_i = check_factor_curvature("G1", N1, E1, MARK1)
    check_factor("G2", N2, E2, MARK2, 8, 7)
    p_j = check_factor_curvature("G2", N2, E2, MARK2)
    n, edges, ix = cartesian(N1, E1, N2, E2)
    x = check_product_structure(n, edges, ix)
    t_h = check_tau(n, edges)
    cof_sum = check_cofactors(n, edges, ix, x)
    p_x = check_curvature_at_x(n, edges, x, t_h, cof_sum)
    check_refutation(p_i, p_j, p_x)
    c1, c2 = check_enumeration()
    check_curvature_law()
    check_census(c1, c2)
    print("NOT RE-RUN: the conjecture refuted above is TRANSCRIBED BY HAND "
          "from the source it is attributed to (%s), namely \"%s\". This "
          "program never reads that preprint, so none of the following is "
          "machine-checked here: that the statement is numbered 4; that it "
          "appears in version v1; that its authors are the ones cited; and, "
          "above all, that its logical form is the one encoded above -- a "
          "NON-STRICT inequality in each hypothesis (so that the exhibited "
          "p_i = p_j = 0 qualifies) and a STRICT inequality in the "
          "conclusion (so that p_x = 0 would already refute it). The check "
          "transcribed_conjecture_predicates_are_the_printed_inequalities "
          "pins those two inequalities in the CODE to the two inequalities "
          "in the TEXT printed at the head, and nothing more: text to source "
          "is unchecked. Every "
          "check above would still pass if that transcription were wrong, "
          "and the paper's arithmetic would still be right while its "
          "attribution was not; a referee must compare the statement "
          "printed at the head of this transcript against the source by "
          "eye. Everything else established above is derived from the two "
          "edge lists alone." % (CONJECTURE_SOURCE, CONJECTURE_TEXT))
    print("NOT RE-RUN: the curvature p_v itself is likewise a transcribed "
          "definition, p_v(G) = 1 - (1/2) sum_{u ~ v} omega_G(u,v) with "
          "unit conductances and omega the effective resistance; it is "
          "computed here by three independent routes, but that it is the "
          "quantity the conjecture speaks of is not checked against the "
          "source either.")

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:                                  # noqa: BLE001
        ck("run_completed_without_an_internal_error", False,
           "%s: %s" % (type(exc).__name__, exc))
    n = len(CHECKS)
    k = sum(1 for _, ok in CHECKS if not ok)
    if k == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
        sys.exit(0)
    else:
        print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
        sys.exit(1)

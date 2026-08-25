#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Verification program for "Disconnected Counterexamples to a Real-Rootedness
Conjecture for Weighted Bond Posets".

Standard library only.  Every decision is exact integer or Fraction arithmetic;
no floating point is used anywhere, and no external file is read.

--------------------------------------------------------------------------
TAKEN FROM THE PAPER (inputs -- these are NOT checks by themselves)
--------------------------------------------------------------------------
  A(t) = 2 + 5t + 2t^2    the Moebius polynomials of the weighted bond posets
  B(t) = 1 + 3t + t^2     of K_3 and of P_3 respectively.  The paper quotes
                          them, together with the multiplicativity
                          mu_{rX} = mu_X^r over disjoint unions, from
                          Equation (1.5) and Proposition 2.2 of the cited
                          article.  The definition of the weighted bond poset
                          lives in that article, which this folder does not
                          reproduce, so A and B are INPUTS here: no poset is
                          enumerated and the identification of A with
                          mu_{K_3} and B with mu_{P_3} is not checked.  The
                          closing NOT RE-RUN line states this and the rest of
                          what is out of scope.
  The printed expansions and constants: A - B = (1+t)^2;
  A^2 - B^2 = (1+t)^2 (3 + 8t + 3t^2), whose quadratic factor has
  discriminant 28; A^3 - B^3 = (1+t)^2 (7 + 37t + 63t^2 + 37t^3 + 7t^4);
  the reduced quadratic 7u^2 + 37u + 49 and its discriminant -3; the leading
  coefficient 5 - 4c of Q_r; and the gamma-expansion of A^r - B^r.  Each is
  recomputed below from A and B alone.  Where BOTH sides of a comparison are
  quantities the paper prints, the check is named a CONSISTENCY check of the
  printed arithmetic: such a check can fail, and a control below shows one
  failing, but it tests the paper's expansions, not its inputs.

--------------------------------------------------------------------------
DERIVED HERE (the checks)
--------------------------------------------------------------------------
  * the graphs rK_3 and rP_3 are built as edge sets: vertex counts, component
    counts by union-find, the spanning-subgraph relation with exactly one edge
    deleted per copy, and the parity of the sign exponent 3r - r.
  * A - B, the values A(-1) = B(-1) = -1, coprimality of A and B by monic gcd
    and by the Sylvester resultant, and A = 2X + t, B = X + t for X = (1+t)^2.
  * the sum-of-squares identity A^2 - 2cAB + B^2 = (A - cB)^2 + (1-c^2)B^2 as
    an identity in Z[c][t]; its leading coefficient 5 - 4c; and the resultant
    Res_t(A - cB, B) computed as a polynomial in c -- it is the constant 1, so
    A - cB and B have no common zero for ANY c, which together with the
    identity is exactly why Q_c(t) > 0 for every real t whenever |c| < 1.
  * an exact real-zero count, WITH MULTIPLICITY, of A^r - B^r for r = 1..12,
    by squarefree decomposition (Yun) plus Sturm sequences over Q.  For
    r = 1, 2 the count equals the degree, so those members are real-rooted;
    for r = 3..12 it is strictly smaller, which is the theorem for those r.
  * the r = 3 member in full: the expansion of A^3 - B^3, its factorization,
    the identity C_3 = A^2 + AB + B^2 = Q(c) at c = cos(2 pi/3) = -1/2, the
    palindromic reduction to 7u^2 + 37u + 49, and the discriminant -3.
  * the r = 2 member, whose quadratic factor has discriminant 28 > 0 and two
    real zeros, so the threshold r >= 3 is exact within this family.
  * the gamma-expansion identity for r = 1..12, positivity of its
    coefficients, and the failure of the same sum truncated one term early.
  * controls: the real-zero counter and the resultant routine on inputs with
    known values; Q_c at c = +-1, where it does have real zeros, so the
    hypothesis |c| < 1 carries weight; and a one-unit mutation of A that
    breaks the printed factorization.
"""

from fractions import Fraction as F
from itertools import permutations
from math import comb, gcd
import sys

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)


def note(msg):
    print("NOTE " + msg)


def finish():
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


# ----------------------------------------------------------------- polynomials
# A polynomial is a list of coefficients, lowest degree first; entries are ints
# or Fractions and the zero polynomial is [].

def ptrim(p):
    q = list(p)
    while q and q[-1] == 0:
        q.pop()
    return q


def pdeg(p):
    return len(ptrim(p)) - 1


def padd(a, b):
    n = max(len(a), len(b))
    return ptrim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                  for i in range(n)])


def psub(a, b):
    n = max(len(a), len(b))
    return ptrim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
                  for i in range(n)])


def pmul(a, b):
    a, b = ptrim(a), ptrim(b)
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y != 0:
                out[i + j] += x * y
    return ptrim(out)


def pscal(a, s):
    return ptrim([s * x for x in a])


def ppow(a, k):
    out = [1]
    for _ in range(k):
        out = pmul(out, a)
    return out


def pderiv(a):
    return ptrim([i * a[i] for i in range(1, len(a))])


def peval(a, x):
    acc = F(0)
    for coef in reversed(ptrim(a)):
        acc = acc * x + coef
    return acc


def pdivmod(a, b):
    """Quotient and remainder over Q."""
    a = [F(x) for x in ptrim(a)]
    b = [F(x) for x in ptrim(b)]
    if not b:
        raise ZeroDivisionError("division by the zero polynomial")
    q = [F(0)] * max(0, len(a) - len(b) + 1)
    while a and len(a) >= len(b):
        coef = a[-1] / b[-1]
        shift = len(a) - len(b)
        q[shift] = coef
        for i in range(len(b)):
            a[shift + i] -= coef * b[i]
        a = ptrim(a)
    return ptrim(q), a


def pdiv_exact(a, b):
    q, r = pdivmod(a, b)
    if pdeg(r) >= 0:
        raise ValueError("division was not exact")
    return q


def pgcd(a, b):
    """Monic greatest common divisor over Q; [] only if both inputs vanish."""
    a = ptrim([F(x) for x in a])
    b = ptrim([F(x) for x in b])
    while b:
        _, r = pdivmod(a, b)
        a, b = b, r
    if not a:
        return []
    lc = a[-1]
    return [x / lc for x in a]


def pnorm_pos(p):
    """Rescale p by a POSITIVE rational to a primitive integer polynomial.

    Positive scaling changes no sign, so this is free to use inside a Sturm
    sequence, and it is what keeps the coefficients from exploding.
    """
    p = ptrim([F(x) for x in p])
    if not p:
        return []
    den = 1
    for x in p:
        d = x.denominator
        den = den * d // gcd(den, d)
    q = [int(x * den) for x in p]
    g = 0
    for x in q:
        g = gcd(g, abs(x))
    if g > 1:
        q = [x // g for x in q]
    return q


def pstr(p, var="t"):
    p = ptrim(p)
    if not p:
        return "0"
    parts = []
    for i, x in enumerate(p):
        if x == 0:
            continue
        s = str(x)
        if i == 0:
            term = s
        elif i == 1:
            term = s + "*" + var
        else:
            term = "%s*%s^%d" % (s, var, i)
        if parts and not term.startswith("-"):
            term = "+" + term
        parts.append(term)
    return "".join(parts)


# --------------------------------------------------------- real zeros, exactly
def sturm_sequence(f):
    """Sturm sequence of a SQUAREFREE nonconstant polynomial, kept over Z."""
    f = pnorm_pos(f)
    if pdeg(f) < 1:
        raise ValueError("Sturm needs a nonconstant polynomial")
    if pdeg(pgcd(f, pderiv(f))) != 0:
        raise ValueError("Sturm needs a squarefree polynomial")
    seq = [f, pnorm_pos(pderiv(f))]
    while pdeg(seq[-1]) > 0:
        _, rem = pdivmod(seq[-2], seq[-1])
        if pdeg(rem) < 0:
            raise ValueError("zero remainder from a squarefree polynomial")
        seq.append(pnorm_pos([-x for x in rem]))
    return seq


def _variations(seq, at_plus_infinity):
    signs = []
    for p in seq:
        d = pdeg(p)
        if d < 0:
            continue
        s = 1 if p[-1] > 0 else -1
        if not at_plus_infinity and d % 2 == 1:
            s = -s
        signs.append(s)
    return sum(1 for i in range(1, len(signs)) if signs[i] != signs[i - 1])


def distinct_real_roots(f):
    """Number of DISTINCT real zeros of f, via its squarefree part."""
    f = ptrim([F(x) for x in f])
    if pdeg(f) < 1:
        return 0
    g = pgcd(f, pderiv(f))
    if pdeg(g) > 0:
        f = pdiv_exact(f, g)
    if pdeg(f) < 1:
        return 0
    seq = sturm_sequence(f)
    return _variations(seq, False) - _variations(seq, True)


def squarefree_decomposition(p):
    """Yun's algorithm: [(k, f_k), ...] with p = lc(p) * prod f_k^k, each f_k
    monic squarefree of positive degree, the f_k pairwise coprime."""
    p = ptrim([F(x) for x in p])
    if pdeg(p) < 1:
        return []
    a = [x / p[-1] for x in p]
    c = pgcd(a, pderiv(a))
    w = pdiv_exact(a, c)
    out = []
    k = 1
    while pdeg(w) > 0:
        y = pgcd(w, c)
        f = pdiv_exact(w, y)
        if pdeg(f) > 0:
            out.append((k, f))
        w = y
        c = pdiv_exact(c, y)
        k += 1
        if k > pdeg(p) + 2:
            raise ValueError("squarefree decomposition did not terminate")
    return out


def sqfree_reconstructs(p):
    """Does lc(p) * prod f_k^k give p back?"""
    p = ptrim([F(x) for x in p])
    if not p:
        return True
    acc = [p[-1]]
    for k, f in squarefree_decomposition(p):
        acc = pmul(acc, ppow(f, k))
    return acc == p


def real_roots_with_multiplicity(p):
    return sum(k * distinct_real_roots(f) for k, f in squarefree_decomposition(p))


# ---------------------------------------------------- polynomials in c over t
# For the symbolic work in the parameter c, a polynomial in t carries
# coefficients that are themselves polynomials in c.

def bt(P):
    Q = [ptrim(x) for x in P]
    while Q and not Q[-1]:
        Q.pop()
    return Q


def badd(P, Q):
    n = max(len(P), len(Q))
    return bt([padd(P[i] if i < len(P) else [], Q[i] if i < len(Q) else [])
               for i in range(n)])


def bsub(P, Q):
    n = max(len(P), len(Q))
    return bt([psub(P[i] if i < len(P) else [], Q[i] if i < len(Q) else [])
               for i in range(n)])


def bmul(P, Q):
    P, Q = bt(P), bt(Q)
    if not P or not Q:
        return []
    out = [[] for _ in range(len(P) + len(Q) - 1)]
    for i, x in enumerate(P):
        for j, y in enumerate(Q):
            out[i + j] = padd(out[i + j], pmul(x, y))
    return bt(out)


def bcmul(P, cpoly):
    return bt([pmul(x, cpoly) for x in P])


def lift(p):
    """A plain t-polynomial, viewed as one with coefficients in Z[c]."""
    return bt([[x] for x in p])


def det_poly(M):
    """Determinant of a square matrix over Z[c], by permutation expansion."""
    n = len(M)
    total = []
    for perm in permutations(range(n)):
        sign = 1
        for i in range(n):
            for j in range(i + 1, n):
                if perm[i] > perm[j]:
                    sign = -sign
        term = [1]
        for i in range(n):
            term = pmul(term, M[i][perm[i]])
            if not term:
                break
        if term:
            total = padd(total, pscal(term, sign))
    return ptrim(total)


def sylvester_resultant(f, g):
    """Res_t(f, g) for f, g in t with coefficients in Z[c], as a poly in c."""
    f, g = bt(f), bt(g)
    if not f or not g:
        raise ValueError("resultant of the zero polynomial")
    m, n = len(f) - 1, len(g) - 1
    size = m + n
    if size == 0:
        return [1]
    M = [[[] for _ in range(size)] for _ in range(size)]
    for i in range(n):
        for j, coef in enumerate(reversed(f)):
            M[i][i + j] = list(coef)
    for i in range(m):
        for j, coef in enumerate(reversed(g)):
            M[n + i][i + j] = list(coef)
    return det_poly(M)


# ----------------------------------------------------------------- the graphs
def triangles_union(r, drop_one_edge):
    """rK_3, or rP_3 when drop_one_edge is True (one edge deleted per copy)."""
    edges = set()
    for i in range(r):
        a, b, d = 3 * i, 3 * i + 1, 3 * i + 2
        edges.add((a, b))
        edges.add((b, d))
        if not drop_one_edge:
            edges.add((a, d))
    return 3 * r, edges


def component_count(n, edges):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u, v in edges:
        a, b = find(u), find(v)
        if a != b:
            parent[a] = b
    return len({find(x) for x in range(n)})


# -------------------------------------------------------- the paper's numbers
A = [2, 5, 2]                                   # mu_{K_3}, a quoted input
B = [1, 3, 1]                                   # mu_{P_3}, a quoted input
X = [1, 2, 1]                                   # X = (1+t)^2
PAPER_A_MINUS_B = [1, 2, 1]                     # (1+t)^2
PAPER_QUAD_R2 = [3, 8, 3]                       # the r = 2 quadratic factor
PAPER_DISC_R2 = 28
PAPER_QUARTIC_R3 = [7, 37, 63, 37, 7]           # the r = 3 quartic factor
PAPER_REDUCED_R3 = [49, 37, 7]                  # 7u^2 + 37u + 49
PAPER_DISC_R3 = -3
PAPER_LEADING_Q = [5, -4]                       # 5 - 4c
R_MAX = 12                                      # exact zero counts up to here
GRID_DEN = 24                                   # rational grid for c
R_GRAPH = 8                                     # graph families up to here


def Q_of_c(c):
    """Q(t) = A^2 - 2c A B + B^2 at a rational c."""
    return padd(psub(pmul(A, A), pscal(pmul(A, B), 2 * F(c))), pmul(B, B))


# ------------------------------------------------------------------- checks
def check_graph_families():
    ok_v = ok_c = ok_sub = ok_sign = True
    for r in range(1, R_GRAPH + 1):
        nG, EG = triangles_union(r, False)
        nH, EH = triangles_union(r, True)
        ok_v = ok_v and nG == 3 * r and nH == 3 * r
        ok_c = ok_c and component_count(nG, EG) == r and component_count(nH, EH) == r
        ok_sub = ok_sub and EH < EG and len(EG) == 3 * r and len(EH) == 2 * r
        ok_sign = ok_sign and (3 * r - r) % 2 == 0
    note("for r = 1..%d, rK_3 has 3r vertices, 3r edges and r components and "
         "rP_3 has 3r vertices, 2r edges and r components" % R_GRAPH)
    ck("both_graphs_have_3r_vertices", ok_v)
    ck("both_graphs_have_exactly_r_connected_components_by_union_find", ok_c)
    ck("rP3_is_a_proper_spanning_subgraph_with_one_edge_deleted_per_copy", ok_sub,
       "E(rP_3) is a proper subset of E(rK_3) and the sizes differ by r")
    ck("sign_exponent_3r_minus_r_is_even_so_the_conjectures_sign_is_plus_one",
       ok_sign, "3r - r = 2r for r = 1..%d" % R_GRAPH)


def check_A_and_B():
    d = psub(A, B)
    ck("consistency_A_minus_B_equals_one_plus_t_squared", d == PAPER_A_MINUS_B,
       "A - B = %s" % pstr(d))
    ck("A_and_B_both_take_the_value_minus_one_at_t_equals_minus_one",
       peval(A, F(-1)) == -1 and peval(B, F(-1)) == -1,
       "A(-1) = %s, B(-1) = %s" % (peval(A, F(-1)), peval(B, F(-1))))
    g = pgcd(A, B)
    res = sylvester_resultant(lift(A), lift(B))
    ck("A_and_B_are_coprime_so_they_share_no_zero_at_all",
       pdeg(g) == 0 and res == [1],
       "monic gcd = %s, Res_t(A, B) = %s" % (pstr(g), pstr(res, "c")))
    ck("A_equals_2X_plus_t_and_B_equals_X_plus_t_for_X_equal_one_plus_t_squared",
       A == padd(pscal(X, 2), [0, 1]) and B == padd(X, [0, 1]))


def check_resultant_routine():
    r1 = sylvester_resultant(lift([-1, 1]), lift([-2, 1]))     # Res(t-1, t-2)
    r2 = sylvester_resultant(lift(A), lift(A))                 # common factor
    r3 = sylvester_resultant(lift([1, 0, 1]), lift([1, 1]))    # Res(t^2+1, t+1)
    ck("control_resultant_routine_reproduces_known_values",
       r1 == [-1] and r2 == [] and r3 == [2],
       "Res(t-1,t-2) = %s, Res(A,A) = %s, Res(t^2+1,t+1) = %s"
       % (pstr(r1, "c"), pstr(r2, "c"), pstr(r3, "c")))


def check_Q_symbolically():
    bA, bB, c = lift(A), lift(B), [0, 1]
    lhs = badd(bsub(bmul(bA, bA), bcmul(bmul(bA, bB), [0, 2])), bmul(bB, bB))
    AmcB = bsub(bA, bcmul(bB, c))
    rhs = badd(bmul(AmcB, AmcB), bcmul(bmul(bB, bB), [1, 0, -1]))
    ck("sum_of_squares_identity_for_Q_holds_identically_in_Z_c_t", lhs == rhs,
       "A^2 - 2cAB + B^2 = (A - cB)^2 + (1 - c^2)B^2, coefficientwise in c")
    lead = lhs[-1] if lhs else []
    ck("leading_coefficient_of_Q_in_t_is_five_minus_four_c",
       lead == PAPER_LEADING_Q and len(lhs) - 1 == 4,
       "deg_t Q = %d, leading coefficient = %s" % (len(lhs) - 1, pstr(lead, "c")))
    res = sylvester_resultant(AmcB, bB)
    note("Res_t(A - cB, B) = %s, a nonzero CONSTANT in c, so A - cB and B have "
         "no common zero for any c whatever; with the identity above that is "
         "why Q(t) > 0 for every real t as soon as 1 - c^2 > 0" % pstr(res, "c"))
    ck("resultant_of_A_minus_cB_and_B_is_the_constant_one_for_every_c",
       res == [1])


def check_Q_on_the_grid():
    bad = []
    for k in range(-(GRID_DEN - 1), GRID_DEN):
        c = F(k, GRID_DEN)
        q = Q_of_c(c)
        if pdeg(q) != 4 or q[0] <= 0 or q[-1] <= 0:
            bad.append((str(c), "degree or an endpoint coefficient"))
        elif distinct_real_roots(q) != 0:
            bad.append((str(c), "has a real zero"))
    ck("Q_has_no_real_zero_at_every_rational_c_in_the_open_unit_interval",
       not bad,
       "%d values c = k/%d with |k| < %d, each with leading and constant "
       "coefficient 5 - 4c > 0 and no real zero%s"
       % (2 * GRID_DEN - 1, GRID_DEN, GRID_DEN,
          "" if not bad else "; failures %s" % bad[:3]))
    q_plus, q_minus = Q_of_c(F(1)), Q_of_c(F(-1))
    n_plus = real_roots_with_multiplicity(q_plus)
    n_minus = real_roots_with_multiplicity(q_minus)
    ck("control_at_c_equal_plus_or_minus_one_Q_does_have_real_zeros",
       q_plus == ppow(PAPER_A_MINUS_B, 2) and q_minus == ppow(padd(A, B), 2)
       and n_plus == 4 and n_minus == 4,
       "Q at c=1 is (A-B)^2 = (1+t)^4 with %d real zeros counted with "
       "multiplicity, Q at c=-1 is (A+B)^2 with %d" % (n_plus, n_minus))


def check_r_equals_three():
    p3 = psub(ppow(A, 3), ppow(B, 3))
    note("A^3 - B^3 = %s" % pstr(p3))
    ck("consistency_printed_factorization_of_A_cubed_minus_B_cubed",
       p3 == pmul(PAPER_A_MINUS_B, PAPER_QUARTIC_R3),
       "(1+t)^2 * (%s) expands to the same polynomial" % pstr(PAPER_QUARTIC_R3))
    C3 = pdiv_exact(p3, PAPER_A_MINUS_B)
    sym = padd(padd(pmul(A, A), pmul(A, B)), pmul(B, B))
    note("at r = 3 the parameter is c = cos(2 pi/3) = -1/2, so the factor Q is "
         "rational and can be compared exactly")
    ck("the_quartic_cofactor_equals_A2_plus_AB_plus_B2_and_equals_Q_at_c_minus_half",
       C3 == sym and C3 == Q_of_c(F(-1, 2)),
       "cofactor = %s" % pstr(C3))
    ck("leading_coefficient_of_the_quartic_is_five_minus_four_c_equals_seven",
       C3[-1] == 7 and peval(PAPER_LEADING_Q, F(-1, 2)) == 7)
    q2, q3, q4 = C3[2], C3[3], C3[4]
    ck("the_quartic_is_palindromic", C3[0] == q4 and C3[1] == q3,
       "coefficients %s" % pstr(C3))
    red = [q2 - 2 * q4, q3, q4]
    back = padd(padd(pscal([1, 0, 2, 0, 1], red[2]),
                     pscal([0, 1, 0, 1, 0], red[1])),
                pscal([0, 0, 1, 0, 0], red[0]))
    ck("reduction_by_u_equals_t_plus_one_over_t_reproduces_the_quartic",
       red == PAPER_REDUCED_R3 and back == C3,
       "t^2 * (%s) = %s" % (pstr(red, "u"), pstr(back)))
    disc = red[1] * red[1] - 4 * red[2] * red[0]
    ck("discriminant_of_the_reduced_quadratic_is_minus_three_so_it_is_positive",
       disc == PAPER_DISC_R3 and disc < 0 and red[2] > 0,
       "disc = %s^2 - 4*%s*%s = %s" % (red[1], red[2], red[0], disc))
    ck("quartic_has_no_real_zero_by_an_independent_sturm_count",
       distinct_real_roots(C3) == 0 and real_roots_with_multiplicity(C3) == 0)


def check_mutation_control():
    Amut = [2, 6, 2]                       # A with one coefficient moved by 1
    pm = psub(ppow(Amut, 3), ppow(B, 3))
    same = (pm == pmul(PAPER_A_MINUS_B, PAPER_QUARTIC_R3))
    _, rem = pdivmod(pm, PAPER_A_MINUS_B)
    ck("control_a_one_unit_mutation_of_A_breaks_the_printed_factorization",
       (not same) and pdeg(rem) >= 0,
       "with A = %s the printed factorization fails and (1+t)^2 no longer "
       "divides the difference" % pstr(Amut))


def check_root_counter_controls():
    cases = [
        ([1, 0, 1], 0),                              # t^2 + 1
        ([-1, 0, 1], 2),                             # t^2 - 1
        ([1, 2, 1], 2),                              # (t+1)^2
        (ppow([1, 0, 1], 2), 0),                     # (t^2+1)^2
        ([0, -1, 0, 1], 3),                          # t^3 - t
        (pmul([1, 0, 1], ppow([-3, 1], 3)), 3),      # (t^2+1)(t-3)^3
        (pmul([1, 2, 1], [3, 8, 3]), 4),             # (1+t)^2 (3+8t+3t^2)
        ([7], 0),                                    # a nonzero constant
    ]
    bad = [(pstr(p), want, real_roots_with_multiplicity(p))
           for p, want in cases if real_roots_with_multiplicity(p) != want]
    ck("control_real_zero_counter_agrees_with_known_counts", not bad,
       "%d control polynomials%s" % (len(cases), "" if not bad else "; %s" % bad))
    ck("control_squarefree_decomposition_reconstructs_the_control_polynomials",
       all(sqfree_reconstructs(p) for p, _ in cases))


def check_the_theorem():
    counts, recon_ok, pred_ok, real_rooted = {}, True, True, {}
    for r in range(1, R_MAX + 1):
        p = psub(ppow(A, r), ppow(B, r))
        pred_ok = pred_ok and pdeg(p) == 2 * r
        recon_ok = recon_ok and sqfree_reconstructs(p)
        m = real_roots_with_multiplicity(p)
        counts[r] = m
        pred_ok = pred_ok and m == (2 if r % 2 else 4)
        real_rooted[r] = (m == 2 * r)
    note("real zeros of A^r - B^r counted with multiplicity, against the "
         "degree 2r: " + ", ".join("r=%d: %d of %d" % (r, counts[r], 2 * r)
                                   for r in sorted(counts)))
    ck("squarefree_decomposition_reconstructs_every_A_r_minus_B_r", recon_ok,
       "r = 1..%d" % R_MAX)
    ck("real_zero_count_is_two_for_odd_r_and_four_for_even_r", pred_ok,
       "the double zero of (1+t)^2, together with the two zeros of A+B when r "
       "is even, and nothing else")
    ck("REFUTATION_A_r_minus_B_r_is_not_real_rooted_for_every_r_from_three",
       all(not real_rooted[r] for r in range(3, R_MAX + 1)),
       "for r = 3..%d the real zeros number %s against degrees %s"
       % (R_MAX, [counts[r] for r in range(3, R_MAX + 1)],
          [2 * r for r in range(3, R_MAX + 1)]))
    ck("the_members_r_equal_one_and_two_ARE_real_rooted_so_the_threshold_is_exact",
       real_rooted[1] and real_rooted[2],
       "r=1: %d of 2, r=2: %d of 4" % (counts[1], counts[2]))


def check_r_equals_two():
    p2 = psub(ppow(A, 2), ppow(B, 2))
    ck("consistency_printed_factorization_of_A_squared_minus_B_squared",
       p2 == pmul(PAPER_A_MINUS_B, PAPER_QUAD_R2) and PAPER_QUAD_R2 == padd(A, B),
       "A^2 - B^2 = (1+t)^2 * (%s), and A + B = %s"
       % (pstr(PAPER_QUAD_R2), pstr(padd(A, B))))
    disc = PAPER_QUAD_R2[1] ** 2 - 4 * PAPER_QUAD_R2[2] * PAPER_QUAD_R2[0]
    ck("discriminant_of_the_r_equal_two_factor_is_twenty_eight_and_positive",
       disc == PAPER_DISC_R2 and disc > 0
       and distinct_real_roots(PAPER_QUAD_R2) == 2,
       "disc = %d, and Sturm finds two real zeros" % disc)


def check_gamma_expansion():
    ident_ok, pos_ok, trunc_ok = True, True, True
    for r in range(1, R_MAX + 1):
        target = psub(ppow(A, r), ppow(B, r))
        terms = []
        for j in range(r):
            coef = comb(r, j) * (2 ** (r - j) - 1)
            pos_ok = pos_ok and coef > 0
            terms.append(pscal(pmul([0] * j + [1], ppow(X, r - j)), coef))
        acc = []
        for term in terms:
            acc = padd(acc, term)
        ident_ok = ident_ok and acc == target
        short = []
        for term in terms[:-1]:
            short = padd(short, term)
        trunc_ok = trunc_ok and short != target
    note("the j = r term of the gamma-expansion carries the coefficient "
         "C(r,r)(2^0 - 1) = 0, which is why the sum stops at j = r-1")
    ck("gamma_expansion_identity_holds_for_every_r_in_range", ident_ok,
       "A^r - B^r = sum_j C(r,j)(2^{r-j} - 1) t^j (1+t)^{2r-2j}, r = 1..%d" % R_MAX)
    ck("gamma_coefficients_are_positive_and_dropping_the_last_term_breaks_it",
       pos_ok and trunc_ok,
       "so the difference is gamma-positive, and no term of the sum is idle")


SUITE = [
    check_graph_families,
    check_A_and_B,
    check_resultant_routine,
    check_Q_symbolically,
    check_Q_on_the_grid,
    check_r_equals_three,
    check_mutation_control,
    check_root_counter_controls,
    check_the_theorem,
    check_r_equals_two,
    check_gamma_expansion,
]


def main():
    print("verification of the disconnected pairs rK_3 vs rP_3 against a "
          "real-rootedness conjecture for weighted bond posets")
    print("python %s, exact integer and Fraction arithmetic only"
          % sys.version.split()[0])
    print("INPUTS quoted from the paper, not derived here: A = mu_{K_3} = %s, "
          "B = mu_{P_3} = %s, and mu_{rX} = mu_X^r" % (pstr(A), pstr(B)))
    print("checks named 'consistency' compare two quantities the paper prints; "
          "they test its expansions, not those inputs")
    for fn in SUITE:
        marker = len(CHECKS)
        try:
            fn()
        except Exception as exc:                       # noqa: BLE001
            ck(fn.__name__, False, "raised %s: %s" % (type(exc).__name__, exc))
        if len(CHECKS) == marker:
            ck(fn.__name__, False, "check produced no result")
    print("NOT RE-RUN: (i) the Moebius polynomials mu_{K_3} = 2 + 5t + 2t^2 and "
          "mu_{P_3} = 1 + 3t + t^2, and the multiplicativity mu_{rX} = mu_X^r, "
          "are INPUTS here, quoted by the paper from Equation (1.5) and "
          "Proposition 2.2 of the cited article; the definition of the weighted "
          "bond poset is in that article and is not reproduced in this folder, "
          "so no poset is enumerated and the identification of A and B with "
          "mu_{K_3} and mu_{P_3} is not checked. (ii) Also transcribed rather "
          "than checked: the numbering of Conjectures 4.13(2), 4.13(1) and "
          "4.12(2), and the reading that 4.13(2) is stated there for arbitrary "
          "k with no connectivity hypothesis, which is the scope claim on which "
          "admitting disconnected graphs depends. (iii) The exact real-zero "
          "count with multiplicity is carried out for r = 1..%d only; for larger "
          "r the theorem rests on the paper's Q_r argument, whose ingredients "
          "are verified above in the parameter c -- the sum-of-squares identity "
          "in Z[c][t], Res_t(A - cB, B) = 1 for every c, the leading "
          "coefficient 5 - 4c, and no real zero at %d rational c in (-1,1) -- "
          "while the factorization A^r - B^r = prod_j (A - zeta^j B) over C and "
          "the value c = cos(2 pi/r) are taken as read, except at r = 3, where "
          "c = -1/2 is rational and the factor is verified exactly. (iv) Two "
          "further truncations of quantifiers over r: the graph-family checks "
          "-- 3r vertices, exactly r components by union-find, rP_3 a proper "
          "spanning subgraph with one edge deleted per copy, and the even sign "
          "exponent 3r - r -- run for r = 1..%d only, and the gamma-expansion "
          "identity together with the positivity of its coefficients "
          "C(r,j)(2^{r-j} - 1) runs for r = 1..%d only, whereas the paper "
          "asserts the graph-family facts for every r >= 3 and gamma-positivity "
          "for every r >= 1; both constructions are uniform in r, but neither "
          "quantifier is exhausted here and no induction is machine-checked. "
          "(v) No trigonometric or complex arithmetic occurs anywhere in this "
          "program: that zeta = e^{2 pi i/r} is nonreal with "
          "|c| = |cos(2 pi/r)| < 1 for every r >= 3, that the leading "
          "coefficient |2 - zeta|^2 equals 5 - 4c, and that Q_r divides "
          "A^r - B^r are all taken as read -- the last verified only at r = 3 "
          "-- and the c-grid is %d rationals k/%d, which contains a genuine "
          "cos(2 pi/r) only in that same r = 3 case; the passage from the "
          "ingredients checked above to the paper's conclusions that "
          "Q_r(t) > 0 for every real t and hence that A^r - B^r fails to be "
          "real-rooted for every r >= 3 -- beyond the members r = 3..%d "
          "decided outright above -- is its own inference, printed here as a "
          "NOTE and not as a check. (vi) Nothing bibliographic or editorial "
          "is fetched or tested: the existence, authorship and contents of "
          "arXiv:2608.08692v1, the novelty or priority of the refutation, and "
          "the paper's negative scope claims -- that the connected restriction "
          "is not addressed, that no assertion is made about "
          "Conjecture 4.13(1), and that these examples do not refute the "
          "gamma-positivity Conjecture 4.12(2) -- are transcribed, as is the "
          "convention that 'real-rooted' means all 2r zeros real, which is the "
          "standard against which the counts above are compared. (vii) Checks "
          "named 'consistency' compare two printed quantities, as disclosed "
          "above." % (R_MAX, 2 * GRID_DEN - 1, R_GRAPH, R_MAX,
                      2 * GRID_DEN - 1, GRID_DEN, R_MAX))
    finish()


if __name__ == "__main__":
    main()

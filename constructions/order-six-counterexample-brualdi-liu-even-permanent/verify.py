#!/usr/bin/env python3
"""Verification program for the order-six counterexample to the Brualdi-Liu
even-permanent conjecture.

The statement under test (the conjecture the paper refutes): for every n >= 3
and every X in the convex hull of the n x n even permutation matrices,
    per^ev(X) := sum over even sigma of prod_i x_{i,sigma(i)}  >=  n!/(2 n^n),
with equality only at the matrix all of whose entries are 1/n.

TAKEN FROM THE PAPER (transcribed inputs, none of them trusted):
  * n = 6, and the integer matrix A with M = A/10 (the exhibited object).
  * The four permutations, in cycle notation, and the integer coefficients
    3, 2, 3, 2 of the claimed decomposition A = 3P_1 + 2P_2 + 3P_3 + 2P_4.
  * The six-row expansion table (image tuple, inversion number, product).
  * The asserted values per^ev(M) = 661/100000, the bound 5/648, the
    cross-multiplication 661*648 = 428328, and the count |A_6| = 360 used in
    the paper's decomposition of J_6/6.

DERIVED HERE (computed, never asserted):
  * A is decoded and printed back; it is integral, nonnegative and doubly
    stochastic after division by 10, so M satisfies the hypotheses.
  * The four permutations are decoded from cycle notation and shown to be even
    by two independent parity computations; the decomposition is rebuilt
    entrywise and its weights shown to be a convex combination, which places M
    in the even hull -- the second hypothesis.
  * per^ev(M) by exhaustive expansion over all 360 even permutations, and
    again by the independent identity per^ev = (per + det)/2 with per by full
    expansion and det by exact elimination.
  * The load-bearing comparison per^ev(M) < 6!/(2*6^6), as exact rationals,
    with the deficit and ratio; and that 6!/(2*6^6) really is 5/648.
  * The paper's expansion table, row by row, and that no further even
    permutation contributes; the support statistics of A.
  * That the bound equals per^ev(J_n/n) for n = 3..6, so the constant compared
    against is the conjectured equality case.
  * That M is a boundary point (it has zero entries), and that the paper's
    interior family M_t = (1-t)M + t J_6/6 has a strictly positive coefficient
    on each of the 360 vertices, rebuilds M_t exactly, and still violates the
    bound for small t -- including the fact, not stated in the paper, that
    per^ev is non-monotone along the segment, so "sufficiently small" is
    essential.
  * Corroboration at other orders: an exhaustive denominator-10 search over
    4-vertex combinations refutes the inequality at n = 4; a weight sweep over
    an exhibited 4-vertex set refutes it at n = 5; and a denominator-30 grid
    over the whole 3-vertex hull at n = 3 attains its minimum exactly at the
    bound, so the conjecture's first order is consistent.
  * A weight sweep at n = 6 showing the paper's (3,2,3,2)/10 is the strict
    minimiser of its simplex at denominators 10, 20 and 30, and that most
    points of that simplex satisfy the inequality.

NOT RE-RUN, and printed as such at the end: no exhaustive census over subsets
of the 360 even permutation matrices of order 6 (or the 60 of order 5).

STANDARD LIBRARY ONLY. Exact integer/rational arithmetic; floats appear only
inside display strings, never in a decision.
"""

import sys
from fractions import Fraction
from itertools import combinations, permutations

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))


def finish():
    n = len(CHECKS)
    bad = [c for c, ok in CHECKS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


def parity_by_cycles(p):
    """Parity of a permutation given as a tuple p with p[i] = image of i (0-based)."""
    n = len(p)
    seen = [False] * n
    par = 0
    for i in range(n):
        if not seen[i]:
            j = i
            c = 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                c += 1
            par += c - 1
    return par % 2


def inversions(p):
    n = len(p)
    return sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])


def is_permutation(p, n):
    return len(p) == n and sorted(p) == list(range(n))


def decode_cycles(cycles, n):
    """Turn 1-based cycle notation, e.g. [[1,2,3,4],[5,6]], into a 0-based tuple."""
    p = list(range(n))
    for cyc in cycles:
        for k in range(len(cyc)):
            src = cyc[k] - 1
            dst = cyc[(k + 1) % len(cyc)] - 1
            p[src] = dst
    return tuple(p)


def perm_matrix(p, n):
    """Matrix with a 1 in position (i, p(i))."""
    return [[1 if p[i] == j else 0 for j in range(n)] for i in range(n)]


def evens(n):
    return [p for p in permutations(range(n)) if parity_by_cycles(p) == 0]


def odds(n):
    return [p for p in permutations(range(n)) if parity_by_cycles(p) == 1]


def perm_sum(X, perms):
    """sum over the given permutations of prod_i X[i][sigma(i)]  (exact)."""
    n = len(X)
    total = 0
    for p in perms:
        prod = 1
        for i in range(n):
            v = X[i][p[i]]
            if v == 0:
                prod = 0
                break
            prod *= v
        total += prod
    return total


def permanent(X):
    return perm_sum(X, list(permutations(range(len(X)))))


def determinant(X):
    """Exact determinant by fraction-free (Bareiss-style) elimination."""
    n = len(X)
    B = [[Fraction(v) for v in row] for row in X]
    det = Fraction(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if B[r][c] != 0:
                piv = r
                break
        if piv is None:
            return Fraction(0)
        if piv != c:
            B[c], B[piv] = B[piv], B[c]
            det = -det
        det *= B[c][c]
        for r in range(c + 1, n):
            f = B[r][c] / B[c][c]
            if f:
                for k in range(c, n):
                    B[r][k] -= f * B[c][k]
    return det


def combine(weights, verts, n):
    """sum_k weights[k] * P_{verts[k]}, entrywise, exactly."""
    X = [[0] * n for _ in range(n)]
    for w, p in zip(weights, verts):
        for i in range(n):
            X[i][p[i]] += w
    return X


def compositions(total, k):
    """All k-tuples of strictly positive integers summing to total."""
    if k == 1:
        yield (total,)
        return
    for a in range(1, total - k + 2):
        for rest in compositions(total - a, k - 1):
            yield (a,) + rest


# ----------------------------------------------------------------------
# VALUES TAKEN FROM THE PAPER (transcription only; nothing below is trusted)
# ----------------------------------------------------------------------
N = 6

A_PAPER = [
    [3, 2, 5, 0, 0, 0],
    [0, 3, 2, 3, 0, 2],
    [0, 0, 3, 4, 3, 0],
    [2, 3, 0, 3, 2, 0],
    [5, 0, 0, 0, 0, 5],
    [0, 2, 0, 0, 5, 3],
]
SCALE = 10  # M = A / 10

# pi_1 = (2 4)(5 6), pi_2 = (1 2 3 4)(5 6), pi_3 = (1 3 5), pi_4 = (1 3 4 5)(2 6)
CYCLES_PAPER = [
    [[2, 4], [5, 6]],
    [[1, 2, 3, 4], [5, 6]],
    [[1, 3, 5]],
    [[1, 3, 4, 5], [2, 6]],
]
COEFFS_PAPER = [3, 2, 3, 2]  # A = 3P1 + 2P2 + 3P3 + 2P4

# The paper's expansion table: (image tuple, inv(sigma), product of entries of A)
TABLE_PAPER = [
    ((1, 3, 4, 5, 6, 2), 4, 480),
    ((1, 4, 3, 2, 6, 5), 4, 2025),
    ((2, 3, 4, 1, 6, 5), 4, 800),
    ((2, 3, 4, 5, 1, 6), 4, 480),
    ((3, 2, 5, 4, 1, 6), 6, 2025),
    ((3, 6, 4, 5, 1, 2), 10, 800),
]
PEREV_PAPER = Fraction(661, 100000)
BOUND_PAPER = Fraction(5, 648)
CROSSMULT_PAPER = 428328  # the paper's 661 * 648
NUM_EVEN_PAPER = 360  # |A_6|, used in the paper's J_6/6 decomposition

M = [[Fraction(a, SCALE) for a in row] for row in A_PAPER]


def check_matrix_well_formed():
    """Check 1: decode the exhibited object, print it back, confirm it is a
    6x6 doubly stochastic matrix (a hypothesis of the statement)."""
    print("exhibited matrix A (M = A/%d):" % SCALE)
    for row in A_PAPER:
        print("   " + " ".join("%2d" % a for a in row))
    shape_ok = len(A_PAPER) == N and all(len(r) == N for r in A_PAPER)
    ints_ok = all(isinstance(a, int) for r in A_PAPER for a in r)
    nonneg = all(a >= 0 for r in A_PAPER for a in r)
    rows = [sum(r) for r in A_PAPER]
    cols = [sum(A_PAPER[i][j] for i in range(N)) for j in range(N)]
    sums_ok = all(s == SCALE for s in rows) and all(s == SCALE for s in cols)
    zeros = sum(1 for r in A_PAPER for a in r if a == 0)
    ck("shape_integral_nonnegative", shape_ok and ints_ok and nonneg,
       "%dx%d, min entry %d" % (N, N, min(a for r in A_PAPER for a in r)))
    ck("M_doubly_stochastic", sums_ok,
       "row sums %s, col sums %s, all == %d" % (rows, cols, SCALE))
    mrows = [sum(M[i]) for i in range(N)]
    mcols = [sum(M[i][j] for i in range(N)) for j in range(N)]
    ck("M_rational_row_col_sums_one",
       all(s == 1 for s in mrows) and all(s == 1 for s in mcols)
       and all(0 <= M[i][j] <= 1 for i in range(N) for j in range(N)),
       "support %d of %d entries, %d zeros" % (N * N - zeros, N * N, zeros))
    return zeros


def check_permutations_even():
    """Check 2: the four permutations of the paper are permutations of
    {1,...,6} and are even, by two independent parity computations."""
    perms = [decode_cycles(c, N) for c in CYCLES_PAPER]
    ok_perm = all(is_permutation(p, N) for p in perms)
    par_cyc = [parity_by_cycles(p) for p in perms]
    par_inv = [inversions(p) % 2 for p in perms]
    ck("four_permutations_valid", ok_perm,
       "images (1-based): " + "; ".join(
           str(tuple(v + 1 for v in p)) for p in perms))
    ck("four_permutations_even",
       par_cyc == [0, 0, 0, 0] and par_inv == [0, 0, 0, 0],
       "cycle parities %s, inversion parities %s, inv = %s"
       % (par_cyc, par_inv, [inversions(p) for p in perms]))
    return perms


def check_hull_membership(perms):
    """Check 3 (hypothesis of the conjecture): M is in the convex hull of the
    even permutation matrices, via the paper's explicit decomposition."""
    rebuilt = combine(COEFFS_PAPER, perms, N)
    entrywise = rebuilt == A_PAPER
    lam = [Fraction(c, SCALE) for c in COEFFS_PAPER]
    ck("decomposition_reproduces_A", entrywise,
       "3P1+2P2+3P3+2P4 equals A" if entrywise else "mismatch")
    ck("convex_weights", sum(lam) == 1 and all(w > 0 for w in lam),
       "weights " + ", ".join(str(w) for w in lam) + "; sum = %s" % sum(lam))
    conv = [[sum(l * perm_matrix(p, N)[i][j] for l, p in zip(lam, perms))
             for j in range(N)] for i in range(N)]
    ck("M_in_even_hull", conv == M,
       "convex combination of the 4 even permutation matrices equals M: %s"
       % (conv == M))


def check_vertex_set(EV):
    """Check 4: the vertex set used by the paper's interior argument."""
    n_even = len(EV)
    distinct = len(set(EV)) == n_even
    all_even = all(parity_by_cycles(p) == 0 for p in EV)
    ck("even_permutation_count",
       n_even == NUM_EVEN_PAPER == 720 // 2 and distinct and all_even,
       "|A_6| = %d (paper uses %d)" % (n_even, NUM_EVEN_PAPER))
    S = [[0] * N for _ in range(N)]
    for p in EV:
        for i in range(N):
            S[i][p[i]] += 1
    flat = set(S[i][j] for i in range(N) for j in range(N))
    ck("sum_of_even_permutation_matrices_is_flat",
       flat == {n_even // N},
       "sum over A_6 of P_sigma = %d * J_6, so J_6/6 = (1/%d) sum P_sigma"
       % (n_even // N, n_even))


def check_bound_constant():
    """Check 5: the paper's constant n!/(2 n^n) at n = 6 really is 5/648."""
    fact = 1
    for i in range(2, N + 1):
        fact *= i
    nn = N ** N
    bound = Fraction(fact, 2 * nn)
    ck("bound_constant_5_over_648",
       bound == BOUND_PAPER and fact == 720 and nn == 46656,
       "6! = %d, 6^6 = %d, 6!/(2*6^6) = %s" % (fact, nn, bound))
    return bound


def check_even_permanent(EV, OD):
    """Check 6/7: per^ev(M) by brute-force expansion over all 360 even
    permutations, and again by the independent identity
    per^ev = (per + det)/2.  Neither value is taken from the paper."""
    ev_int = perm_sum(A_PAPER, EV)
    od_int = perm_sum(A_PAPER, OD)
    perev = Fraction(ev_int, SCALE ** N)
    ck("even_permanent_bruteforce", perev == PEREV_PAPER,
       "sum over A_6 = %d, per^ev(M) = %d/%d = %s (paper says %s)"
       % (ev_int, ev_int, SCALE ** N, perev, PEREV_PAPER))
    per = permanent(A_PAPER)
    det = determinant(A_PAPER)
    half = Fraction(per + det, 2)
    ck("even_permanent_via_per_plus_det",
       half == ev_int and per == ev_int + od_int and det == ev_int - od_int,
       "per(A) = %d, det(A) = %s, (per+det)/2 = %s, odd part = %d"
       % (per, det, half, od_int))
    return perev, ev_int


def check_expansion_table(EV):
    """Check 8: the paper's six-row table is exactly the set of even
    permutations contributing a nonzero term, with the stated inversion
    numbers and products."""
    found = []
    for p in EV:
        prod = 1
        for i in range(N):
            prod *= A_PAPER[i][p[i]]
        if prod:
            found.append((tuple(v + 1 for v in p), inversions(p), prod))
    found.sort()
    want = sorted(TABLE_PAPER)
    ck("expansion_table_exact", found == want,
       "%d nonzero even terms; %s"
       % (len(found), "table reproduced" if found == want
          else "computed %s" % (found,)))
    ck("expansion_table_sums",
       sum(t[2] for t in found) == 6610
       and sum(t[2] for t in want) == sum(t[2] for t in found),
       "computed terms sum to %d over 10^6 (the paper's six rows sum to %d)"
       % (sum(t[2] for t in found), sum(t[2] for t in want)))
    supp = [[1 if A_PAPER[i][j] else 0 for j in range(N)] for i in range(N)]
    n_supp = permanent(supp)
    n_even_supp = perm_sum(supp, EV)
    zeros = sum(1 for r in A_PAPER for a in r if a == 0)
    ck("support_pattern",
       zeros == 17 and n_even_supp == 6 and n_supp == 20,
       "%d zero entries; %d permutations lie in the support, %d even and "
       "%d odd; the other %d even permutations meet a zero entry"
       % (zeros, n_supp, n_even_supp, n_supp - n_even_supp,
          len(EV) - n_even_supp))


def check_violation(perev, bound):
    """Check 9 (load bearing): the computed per^ev(M) is strictly below the
    conjectured bound, by exact rational comparison."""
    ck("VIOLATES_conjectured_bound", perev < bound,
       "per^ev(M) = %s = %.8f  <  %s = %.8f"
       % (perev, float(perev), bound, float(bound)))
    lhs = perev.numerator * bound.denominator
    rhs = bound.numerator * perev.denominator
    ck("violation_cross_multiplication",
       lhs == CROSSMULT_PAPER and rhs == 500000 and lhs < rhs,
       "cross-multiplying gives %d < %d (the paper's 661*648 = %d against "
       "5*100000); deficit = %s, ratio to the bound = %s"
       % (lhs, rhs, CROSSMULT_PAPER, bound - perev, perev / bound))


def check_barycentre(bound):
    """Check 10: the conjectured bound is exactly the value of per^ev at the
    flat matrix J_n/n, for n = 3,...,6 -- so the constant M is compared with
    is the claimed equality case, not an arbitrary number."""
    detail = []
    ok = True
    for n in range(3, N + 1):
        fact = 1
        for i in range(2, n + 1):
            fact *= i
        EVn = evens(n)
        flat = [[Fraction(1, n)] * n for _ in range(n)]
        val = perm_sum(flat, EVn)
        want = Fraction(fact, 2 * n ** n)
        ok = ok and val == want and len(EVn) == fact // 2
        detail.append("n=%d: %s" % (n, val))
    ck("bound_is_value_at_flat_matrix", ok, "; ".join(detail)
       + " (n=6 equals %s)" % bound)


def check_boundary(zeros):
    """Check 11: M itself has zero entries, hence lies on the facet x_ij >= 0
    of the polytope and is NOT in the relative interior -- which is why the
    paper needs a separate interior claim."""
    zpos = [(i + 1, j + 1) for i in range(N) for j in range(N)
            if A_PAPER[i][j] == 0]
    ck("M_is_a_boundary_point", len(zpos) == zeros and zeros > 0,
       "%d zero entries, e.g. positions %s -> M lies on a facet"
       % (zeros, zpos[:4]))


def interior_point(t, EV, perms):
    """M_t = (1-t)M + t J_6/6 together with its coefficient vector over all
    360 vertices, built from the paper's two decompositions."""
    coeff = dict((p, t / len(EV)) for p in EV)
    for w, p in zip(COEFFS_PAPER, perms):
        if p not in coeff:
            # p is not one of the 360 vertices (it is not even), so its mass
            # cannot be placed on the vertex set at all; leaving it out makes
            # the coefficients sum to less than 1 and the check below fail,
            # rather than raising and skipping the remaining checks.
            continue
        coeff[p] += (1 - t) * Fraction(w, SCALE)
    Mt = [[(1 - t) * M[i][j] + t * Fraction(1, N) for j in range(N)]
          for i in range(N)]
    return Mt, coeff


def check_relative_interior(EV, perms, bound):
    """Check 12/13: for t in (0,1) the point M_t is a convex combination with
    a strictly positive coefficient on every vertex (hence relatively
    interior), and per^ev(M_t) is still strictly below the bound for small t."""
    t = Fraction(1, 100)
    Mt, coeff = interior_point(t, EV, perms)
    rebuilt = [[sum(coeff[p] for p in EV if p[i] == j) for j in range(N)]
               for i in range(N)]
    pos = all(c > 0 for c in coeff.values())
    entries_pos = all(Mt[i][j] > 0 for i in range(N) for j in range(N))
    verts_ok = set(perms) <= set(EV)
    ck("interior_point_positive_coefficients",
       verts_ok and sum(coeff.values()) == 1 and pos and rebuilt == Mt
       and entries_pos,
       "t = %s: %d vertex coefficients, min %s, sum %s, rebuild M_t: %s, "
       "all entries of M_t positive: %s"
       % (t, len(coeff), min(coeff.values()), sum(coeff.values()),
          rebuilt == Mt, entries_pos))
    rows = []
    ok = True
    for tt in (Fraction(1, 32), Fraction(1, 100), Fraction(1, 1000),
               Fraction(1, 10000)):
        val = perm_sum(interior_point(tt, EV, perms)[0], EV)
        ok = ok and val < bound
        rows.append("t=%s: per^ev=%.9f" % (tt, float(val)))
    ck("interior_points_also_violate", ok,
       "; ".join(rows) + "; all exactly below %s = %.9f"
       % (bound, float(bound)))
    dyadic = []
    for e in range(1, 8):
        tt = Fraction(1, 2 ** e)
        val = perm_sum(interior_point(tt, EV, perms)[0], EV)
        dyadic.append((tt, val < bound))
    obey = [str(t) for t, b in dyadic if not b]
    viol = [str(t) for t, b in dyadic if b]
    ck("interior_violation_needs_small_t",
       (not dyadic[0][1]) and dyadic[-1][1]
       and all(b for t, b in dyadic if t <= Fraction(1, 32)),
       "per^ev along the segment is not monotone: dyadic t obeying the bound "
       "%s, violating it %s -- so 'sufficiently small t' is essential"
       % (obey, viol))
    at_zero = perm_sum(interior_point(Fraction(0), EV, perms)[0], EV)
    at_one = perm_sum(interior_point(Fraction(1), EV, perms)[0], EV)
    ck("interior_family_endpoints",
       at_zero == PEREV_PAPER and at_one == bound,
       "t=0 gives per^ev = %s (should be the exhibited M's value %s), t=1 "
       "gives %s (the flat matrix, should be the bound %s)"
       % (at_zero, PEREV_PAPER, at_one, bound))


def check_invariance(EV, OD, ev_int):
    """Check 14: per^ev is unchanged by permuting rows by an even permutation
    and swaps to the odd permanent for an odd one.  This licenses the
    normalisation used by the searches below."""
    od_int = perm_sum(A_PAPER, OD)
    bad_even = 0
    for r in EV:
        Y = [A_PAPER[r[i]] for i in range(N)]
        if perm_sum(Y, EV) != ev_int:
            bad_even += 1
    bad_odd = 0
    for r in OD[:60]:
        Y = [A_PAPER[r[i]] for i in range(N)]
        if perm_sum(Y, EV) != od_int:
            bad_odd += 1
    ck("even_permanent_invariance", bad_even == 0 and bad_odd == 0,
       "per^ev(P_rho A) = %d for all %d even rho; = %d (the odd permanent) "
       "for %d odd rho tested" % (ev_int, len(EV), od_int, len(OD[:60])))


def search_min(n, verts_pool, k, total, EVn, fixed=None):
    """Exhaustive minimum of per^ev over integer/total convex combinations of
    k vertices drawn from verts_pool (with the first one fixed if given)."""
    best = None
    ws = list(compositions(total, k))
    base = () if fixed is None else (fixed,)
    pool = [p for p in verts_pool if p != fixed]
    for sub in combinations(pool, k - len(base)):
        vs = base + sub
        for w in ws:
            v = perm_sum(combine(w, vs, n), EVn)
            if best is None or v < best[0]:
                best = (v, w, vs)
    return best[0], Fraction(best[0], total ** n), best[1], best[2]


def check_weight_sweep(EV, perms, bound, perev):
    """Check 15: sweep the weights on the paper's four vertices.  The paper's
    (3,2,3,2)/10 is the strict minimiser over every positive integer weight
    vector of denominator 10, 20 and 30, and it is a delicate choice: most
    points of that simplex satisfy the conjecture."""
    lines = []
    ok = True
    for D in (10, 20, 30):
        vals = []
        for w in compositions(D, 4):
            vals.append((perm_sum(combine(w, perms, N), EV), w))
        mn = min(vals)
        below = sum(1 for v, w in vals if Fraction(v, D ** N) < bound)
        ok = ok and Fraction(mn[0], D ** N) == perev
        ok = ok and mn[1] == tuple(c * D // 10 for c in COEFFS_PAPER)
        lines.append("D=%d: %d points, min %s at %s, %d below the bound"
                     % (D, len(vals), Fraction(mn[0], D ** N), mn[1], below))
    ck("paper_point_minimises_its_simplex", ok, "; ".join(lines))


def check_small_order_failures():
    """Check 16/17: the conjecture already fails at n = 4 and n = 5, found by
    search here rather than quoted.  Exhaustive at n = 4 over all 4-vertex
    combinations of denominator 10 (one vertex normalised to the identity by
    the invariance above); at n = 5 a verified exhibited combination plus an
    exhaustive weight sweep of its simplex."""
    EV4 = evens(4)
    bound4 = Fraction(24, 2 * 4 ** 4)
    v, frac4, w4, vs4 = search_min(4, EV4, 4, 10, EV4, fixed=tuple(range(4)))
    ck("order_four_search_refutes", frac4 < bound4,
       "min per^ev over C(11,3)*84 combinations = %s < %s = 4!/(2*4^4); "
       "weights %s on %s" % (frac4, bound4, w4,
                             [tuple(x + 1 for x in p) for p in vs4]))
    EV5 = evens(5)
    bound5 = Fraction(120, 2 * 5 ** 5)
    vs5 = [(0, 1, 2, 3, 4), (0, 2, 1, 4, 3), (1, 0, 3, 2, 4), (3, 4, 2, 0, 1)]
    ok5 = all(is_permutation(p, 5) and parity_by_cycles(p) == 0 for p in vs5)
    best = None
    for w in compositions(20, 4):
        val = perm_sum(combine(w, vs5, 5), EV5)
        if best is None or val < best[0]:
            best = (val, w)
    frac5 = Fraction(best[0], 20 ** 5)
    ck("order_five_search_refutes", ok5 and frac5 < bound5,
       "min per^ev over the %d weightings of denominator 20 on 4 even "
       "vertices = %s < %s = 5!/(2*5^5); weights %s"
       % (len(list(compositions(20, 4))), frac5, bound5, best[1]))


def check_order_three_holds():
    """Check 18: the conjecture is not vacuous at its first order.  At n = 3
    the hull has only 3 vertices; a full grid of denominator 30 attains its
    minimum exactly at the flat matrix, with value 3!/(2*3^3) = 1/9."""
    EV3 = evens(3)
    bound3 = Fraction(6, 2 * 27)
    D = 30
    best = None
    for a in range(D + 1):
        for b in range(D - a + 1):
            w = (a, b, D - a - b)
            val = Fraction(perm_sum(combine(w, EV3, 3), EV3), D ** 3)
            if best is None or val < best[0]:
                best = (val, w)
    ck("order_three_grid_minimum_is_the_bound",
       best[0] == bound3 and best[1] == (D // 3, D // 3, D // 3)
       and len(EV3) == 3,
       "grid minimum %s attained at weights %s (the flat matrix); bound is %s"
       % (best[0], best[1], bound3))


def main():
    print("order-six counterexample to the Brualdi-Liu even-permanent "
          "conjecture: verification")
    print("conjecture under test: per^ev(X) >= n!/(2 n^n) for all X in the "
          "convex hull of the even permutation matrices, n >= 3")
    EV = evens(N)
    OD = odds(N)
    zeros = check_matrix_well_formed()
    perms = check_permutations_even()
    check_hull_membership(perms)
    check_vertex_set(EV)
    bound = check_bound_constant()
    perev, ev_int = check_even_permanent(EV, OD)
    check_expansion_table(EV)
    check_violation(perev, bound)
    check_barycentre(bound)
    check_boundary(zeros)
    check_relative_interior(EV, perms, bound)
    check_invariance(EV, OD, ev_int)
    check_weight_sweep(EV, perms, bound, perev)
    check_order_three_holds()
    check_small_order_failures()
    print("NOT RE-RUN: no exhaustive census over all 4-element subsets of the "
          "360 even permutation matrices of order 6 (nor of the 60 of order "
          "5) was attempted; the order-6 sweep is over the weights on the "
          "four exhibited vertices only, and the order-5 refutation sweeps "
          "the weights on one exhibited 4-vertex set. The order-4 refutation "
          "is exhaustive over 4-vertex sets at denominator 10.")


if __name__ == "__main__":
    # Any unexpected exception is reported as a failed check so that the
    # verdict line is always emitted and a crash is never mistaken for a
    # clean run.
    try:
        main()
    except Exception as exc:
        ck("no_unhandled_error", False,
           "%s: %s" % (type(exc).__name__, exc))
    finish()

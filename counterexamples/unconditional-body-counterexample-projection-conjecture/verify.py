#!/usr/bin/env python3
"""
Verification program for the counterexample to the conjecture

    |P_E K|_n <= 2^n |P_E (K cap R_+^m)|_n

for unconditional convex bodies K in R^m and n-dimensional subspaces E.

TAKEN FROM THE PAPER (inputs, copied verbatim, never used as a check target
unless independently recomputed here):
  * the body K_0 = { x in R^4 : 4|x1|+3|x3|+3|x4| <= 1, 3|x2|+4|x3| <= 1 };
  * the subspace generators b1 = (3,-1,1,2), b2 = (-3,1,2,0);
  * the claimed vertex list V_+ of K_+ = K_0 cap R_+^4 (9 points), and the
    claim that of the 15 four-element active sets, 9 are feasible and
    nonsingular, 4 are singular, and 2 are infeasible, the infeasible ones
    being (0,0,1/3,0) and (0,-1/9,1/3,0);
  * the claimed Gram matrix G = B^T B = [[15,-8],[-8,14]] with det G = 146;
  * the claimed vertex lists of Q_+ = Phi(K_+) (5 points) and Q = Phi(K_0)
    (8 points), and the two half-space certificates for them;
  * the claimed areas area(Q_+) = 7/12, area(Q) = 1399/576, and the claimed
    ratio 1399/336 = 4 + 55/336.

DERIVED HERE (recomputed from the definitions with exact integer/rational
arithmetic only; no floating point enters any decision):
  * the 12 linear inequalities used below really do ENCODE the paper's
    absolute-value definition of K_0: they are exactly the union of the two
    sign orbits, and a literal transcription of "4|x1|+3|x3|+3|x4| <= 1 and
    3|x2|+4|x3| <= 1" is compared against them pointwise.  Nothing else in the
    program ever looks at the absolute-value form, so without this check every
    other check could pass for the wrong body;
  * K_0 has 0 in its interior, is invariant under all 16 coordinate sign
    changes, and is NOT invariant under any nonidentity coordinate permutation
    (the paper's scope remark; in particular K_0 is not the unit ball of any
    l_q norm, since those are permutation invariant);
  * x in K_0 iff |x| in K_+, which is the mechanism behind the equality case
    n = m -- proved exactly (row set closed under sign changes and under
    coefficientwise absolute value), with a rational grid as a cross-check;
  * K_0 and K_+ are BOUNDED, decided exactly by showing their recession cones
    are {0} via complete extreme-ray enumeration (boundedness is load-bearing
    twice: a polyhedron equals the hull of its vertices only if bounded, and
    both the sign-hull step and the two half-space certificates below use
    exactly that implication);
  * the six inequalities used for K_+ define exactly K_0 cap R_+^4, proved by
    row containment in one direction and coefficientwise domination on the
    orthant in the other, rather than by sampling a finite grid;
  * the complete vertex set of K_+ by exhaustive enumeration of all 15
    four-element active sets of its 6 defining inequalities, with the
    singular/infeasible census;
  * the complete vertex set of K_0 by exhaustive enumeration of all 495
    four-element active sets of its 12 linear inequalities, and its equality
    with the sign-orbit of V_+;
  * G and det G from b1, b2, and the area-distortion factor of Phi via exact
    Gram determinants of triangles inside E_0 (squared, to stay rational);
  * the boundedness of the two half-space certificates, their vertex sets,
    their agreement with the projected generators, their counterclockwise
    order, and their exact shoelace areas;
  * the ratio area(Q)/area(Q_+) and the strict failure of the conjectured
    bound for (m,n) = (4,2) and for every m >= 4, 2 <= n <= m-2.

NOT COVERED HERE (also printed at the end of the run, so a reader of the output
alone is not misled):
  * the positive cases n = 1 and n = m-1 of the classification, which are
    quoted from the cited work;
  * the case n = m, which is supported only through the orthant reduction, not
    by a 4-dimensional volume computation;
  * the product identities P_E K = P_{E_0}K_0 x [-1,1]^{n-2} and
    P_E(K cap R_+^m) = P_{E_0}K_+ x [0,1]^{n-2}, which are verified as sets
    only for the single pair (m,n) = (5,3) and only on a rational grid; the
    all-(m,n) statement is then the arithmetic consequence 2^{n-2}r > 2^n,
    enumerated for m <= 20;
  * the statement of the refuted conjecture itself, and in particular its
    hypothesis class.  What is proved below is that K_0 is invariant under all
    16 coordinate sign changes (hence unconditional) and is invariant under no
    nonidentity coordinate permutation; that this sign-invariance is the whole
    hypothesis of the conjecture is transcribed from the cited work, which is
    not part of this material and is not read by this program;
  * every bibliographic datum: the preprint identifier, the author names, and
    the pinpoint citations to the two propositions supplying the positive
    cases n = 1 and n = m-1.  The negative half of the classification is
    recomputed here in full; the positive half rests on those citations.
"""

from fractions import Fraction as F
import functools
import itertools
import sys

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
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


def solve_exact(rows, rhs):
    """Exact Gaussian elimination. Returns the unique solution or None if the
    square system is singular."""
    n = len(rhs)
    M = [[F(rows[i][j]) for j in range(n)] + [F(rhs[i])] for i in range(n)]
    for c in range(n):
        p = None
        for r in range(c, n):
            if M[r][c] != 0:
                p = r
                break
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][j] - f * M[c][j] for j in range(n + 1)]
    return tuple(M[i][n] for i in range(n))


def satisfies(x, ineqs):
    """ineqs is a list of (a, b) meaning a.x <= b."""
    for a, b in ineqs:
        if sum(F(ai) * xi for ai, xi in zip(a, x)) > F(b):
            return False
    return True


def vertices_of(ineqs, dim):
    """All vertices of {x : a.x <= b for all (a,b)}, by exhaustive enumeration
    of dim-element active sets.  Also returns the singular / infeasible census
    and the number of ACTIVE SETS that were feasible.  That last number is what
    the paper's census counts, and it is not the same quantity as the number of
    distinct vertices: two different active sets can pin the same point."""
    verts = set()
    singular = 0
    infeasible = []
    feasible_sets = 0
    for sub in itertools.combinations(range(len(ineqs)), dim):
        rows = [ineqs[i][0] for i in sub]
        rhs = [ineqs[i][1] for i in sub]
        x = solve_exact(rows, rhs)
        if x is None:
            singular += 1
        elif satisfies(x, ineqs):
            verts.add(x)
            feasible_sets += 1
        else:
            infeasible.append(x)
    return verts, singular, infeasible, feasible_sets


def nullspace(rows, ncols):
    """Exact basis of {d : r.d = 0 for every row r}."""
    M = [[F(v) for v in r] for r in rows]
    pivots = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [v / pv for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [M[i][j] - f * M[r][j] for j in range(ncols)]
        pivots.append(c)
        r += 1
    basis = []
    for fc in [c for c in range(ncols) if c not in pivots]:
        d = [F(0)] * ncols
        d[fc] = F(1)
        for i, pc in enumerate(pivots):
            d[pc] = -M[i][fc]
        basis.append(tuple(d))
    return basis


def cone_trivial(normals, dim):
    """Exactly decide whether the recession cone C = {d : a.d <= 0 for every
    normal a} is {0}, i.e. whether {x : a.x <= b} is bounded (boundedness does
    not depend on the right-hand sides).

    If rank(normals) < dim then ker(normals) is a nonzero subspace of C.
    Otherwise C is pointed, so C != {0} implies C has an extreme ray, and an
    extreme ray is a 1-dimensional face, hence carries dim-1 linearly
    independent active normals.  Enumerating all (dim-1)-subsets whose
    nullspace is a line therefore finds an extreme ray if one exists, so the
    test below is complete, not merely necessary."""
    normals = list(normals)
    if rank(normals) < dim:
        return False
    for sub in itertools.combinations(range(len(normals)), dim - 1):
        line = nullspace([normals[i] for i in sub], dim)
        if len(line) != 1:
            continue
        d = line[0]
        if all(x == 0 for x in d):
            continue
        for sgn in (F(1), F(-1)):
            e = tuple(sgn * x for x in d)
            if all(dot(a, e) <= 0 for a in normals):
                return False
    return True


def ccw_sorted(points):
    """Sort planar points counterclockwise about their centroid, exactly."""
    k = len(points)
    cx = sum(p[0] for p in points) / k
    cy = sum(p[1] for p in points) / k

    def quad(p):
        dx, dy = p[0] - cx, p[1] - cy
        if dx > 0 and dy >= 0:
            return 0
        if dx <= 0 and dy > 0:
            return 1
        if dx < 0 and dy <= 0:
            return 2
        return 3

    def cmp(p, q):
        if quad(p) != quad(q):
            return -1 if quad(p) < quad(q) else 1
        cr = (p[0] - cx) * (q[1] - cy) - (p[1] - cy) * (q[0] - cx)
        if cr > 0:
            return -1
        if cr < 0:
            return 1
        return 0

    return sorted(points, key=functools.cmp_to_key(cmp))


def shoelace2(poly):
    """Twice the signed area of an ordered polygon, exact."""
    s = F(0)
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]
        s += x1 * y2 - x2 * y1
    return s


# ---------------------------------------------------------------------------
# Data taken from the paper.
# ---------------------------------------------------------------------------

# K_0 = { x : 4|x1| + 3|x3| + 3|x4| <= 1,  3|x2| + 4|x3| <= 1 }.
# The two absolute-value forms of the paper, as (coefficients of |x_j|, rhs).
PAPER_ABS_FORMS = [((4, 0, 3, 3), 1), ((0, 3, 4, 0), 1)]

# Expanded into 8 + 4 = 12 linear inequalities a.x <= 1.
K0_INEQS = []
for s1, s3, s4 in itertools.product((1, -1), repeat=3):
    K0_INEQS.append(((4 * s1, 0, 3 * s3, 3 * s4), 1))
for s2, s3 in itertools.product((1, -1), repeat=2):
    K0_INEQS.append(((0, 3 * s2, 4 * s3, 0), 1))

# K_+ = K_0 cap R_+^4: the six inequalities of the paper.
KPLUS_INEQS = [((4, 0, 3, 3), 1), ((0, 3, 4, 0), 1),
               ((-1, 0, 0, 0), 0), ((0, -1, 0, 0), 0),
               ((0, 0, -1, 0), 0), ((0, 0, 0, -1), 0)]

V_PLUS_PAPER = {
    (F(0), F(0), F(0), F(0)),
    (F(0), F(0), F(0), F(1, 3)),
    (F(0), F(0), F(1, 4), F(0)),
    (F(0), F(0), F(1, 4), F(1, 12)),
    (F(0), F(1, 3), F(0), F(0)),
    (F(0), F(1, 3), F(0), F(1, 3)),
    (F(1, 16), F(0), F(1, 4), F(0)),
    (F(1, 4), F(0), F(0), F(0)),
    (F(1, 4), F(1, 3), F(0), F(0)),
}

INFEASIBLE_PAPER = {
    (F(0), F(0), F(1, 3), F(0)),
    (F(0), F(-1, 9), F(1, 3), F(0)),
}

B1 = (3, -1, 1, 2)
B2 = (-3, 1, 2, 0)
G_PAPER = ((15, -8), (-8, 14))
DETG_PAPER = 146

QPLUS_PAPER = [(F(-1, 3), F(1, 3)), (F(3, 4), F(-3, 4)), (F(2, 3), F(0)),
               (F(5, 12), F(1, 2)), (F(1, 4), F(1, 2))]

Q_PAPER = [(F(-13, 12), F(13, 12)), (F(-1), F(1, 3)),
           (F(-5, 12), F(-1, 2)), (F(-1, 16), F(-11, 16)),
           (F(13, 12), F(-13, 12)), (F(1), F(-1, 3)),
           (F(5, 12), F(1, 2)), (F(1, 16), F(11, 16))]

# Half-space certificate for Q_+ (s+t >= 0 written as -s-t <= 0).
HQPLUS = [((-1, -1), 0), ((9, 1), 6), ((6, 3), 4), ((0, 2), 1), ((-2, 7), 3)]

# Half-space certificate for Q: four symmetric pairs |a.y| <= b.
HQ = []
for (a, b) in [((27, 3), 26), ((30, 21), 23), ((36, 68), 49), ((19, 55), 39)]:
    HQ.append((a, b))
    HQ.append(((-a[0], -a[1]), b))

AREA_QPLUS_PAPER = F(7, 12)
AREA_Q_PAPER = F(1399, 576)
RATIO_PAPER = F(1399, 336)


def phi(x):
    """Phi(x) = B^T x, with B the 4x2 matrix of columns b1, b2."""
    return (sum(F(a) * F(xi) for a, xi in zip(B1, x)),
            sum(F(a) * F(xi) for a, xi in zip(B2, x)))


def paper_membership(x):
    """A LITERAL transcription of the paper's definition of K_0,

        4|x_1| + 3|x_3| + 3|x_4| <= 1   and   3|x_2| + 4|x_3| <= 1,

    written out by hand so that it shares no code with the 12-row expansion it
    is used to test."""
    return (4 * abs(x[0]) + 3 * abs(x[2]) + 3 * abs(x[3]) <= 1
            and 3 * abs(x[1]) + 4 * abs(x[2]) <= 1)


def check_encoding_matches_the_paper():
    """Everything below reads K_0 through the 12 linear rows K0_INEQS, so if the
    expansion of the absolute values were wrong the whole run could pass for the
    wrong body.  Two independent checks pin the encoding.

    (exact) The row set is exactly the union of the two sign orbits of the
    paper's coefficient vectors, and those vectors are nonnegative.  For a
    nonnegative c and any x in R^4, max over eps in {-1,1}^4 of (eps.c).x equals
    sum_j c_j|x_j| (take eps_j = sign(x_j)), so satisfying every row of an orbit
    is EQUIVALENT to sum_j c_j|x_j| <= b -- for all of R^4, not just on a grid.

    (cross-check) The hand-written formula is evaluated against the 12 rows at
    every point of a rational grid."""
    rows = set((tuple(a), b) for a, b in K0_INEQS)
    orbit = set()
    for c, b in PAPER_ABS_FORMS:
        for eps in itertools.product((1, -1), repeat=4):
            orbit.add((tuple(e * cj for e, cj in zip(eps, c)), b))
    ck("K0_rows_are_exactly_the_two_sign_orbits_of_the_paper_forms",
       rows == orbit and len(rows) == 12,
       "%d rows, %d orbit rows, symmetric difference %d"
       % (len(rows), len(orbit), len(rows ^ orbit)))
    ck("paper_form_coefficients_are_nonnegative_so_orbit_max_is_the_abs_form",
       all(all(cj >= 0 for cj in c) for c, _ in PAPER_ABS_FORMS),
       "forms %s" % (PAPER_ABS_FORMS,))
    grid = [F(k, 12) for k in range(-6, 7)]
    bad = 0
    tested = 0
    for x in itertools.product(grid, repeat=4):
        tested += 1
        if satisfies(x, K0_INEQS) != paper_membership(x):
            bad += 1
    ck("K0_rows_agree_with_the_literal_paper_formula_on_a_grid", bad == 0,
       "%d grid points, %d mismatches" % (tested, bad))


def check_body_is_a_body():
    """K_0 is a full-dimensional bounded convex body: the origin is strictly
    interior, and every coordinate direction is blocked by some inequality."""
    origin = (F(0), F(0), F(0), F(0))
    slacks = [F(b) - dot(ai, origin) for ai, b in K0_INEQS]
    interior = all(s > 0 for s in slacks)
    ck("K0_origin_strictly_interior", interior, "min slack %s" % min(slacks))
    blocked = []
    for i in range(4):
        for sgn in (1, -1):
            d = [0, 0, 0, 0]
            d[i] = sgn
            blocked.append(any(sum(F(a) * dj for a, dj in zip(ai, d)) > 0
                               for ai, _ in K0_INEQS))
    # Necessary condition only: blocking the 8 axis rays does NOT imply
    # boundedness (the unbounded strip |x1+x2| <= 1 blocks all four axis rays
    # in the plane).  Kept for information; the decision is made by the
    # complete recession-cone test on the next line.
    ck("K0_axis_rays_blocked_necessary_condition_only", all(blocked),
       "%d/8 blocked" % sum(blocked))
    ck("K0_bounded_recession_cone_is_trivial",
       cone_trivial([a for a, _ in K0_INEQS], 4),
       "complete extreme-ray enumeration over 12 normals")
    ck("Kplus_bounded_recession_cone_is_trivial",
       cone_trivial([a for a, _ in KPLUS_INEQS], 4),
       "complete extreme-ray enumeration over 6 normals")


def check_unconditional():
    """Every coordinate sign change permutes the 12 inequalities of K_0, and no
    nonidentity coordinate PERMUTATION does (the paper's scope remark)."""
    base = set((tuple(a), b) for a, b in K0_INEQS)
    ok = len(base) == 12
    for eps in itertools.product((1, -1), repeat=4):
        img = set((tuple(e * ai for e, ai in zip(eps, a)), b)
                  for a, b in K0_INEQS)
        if img != base:
            ok = False
    ck("K0_invariant_under_all_16_sign_changes", ok, "12 distinct facets")
    # Scope claim of the paper: K_0 is not permutation invariant, hence is not
    # the unit ball of any l_q norm (those are permutation invariant).  Only the
    # identity may fix the facet set.
    fixing = [p for p in itertools.permutations(range(4))
              if set((tuple(a[p[j]] for j in range(4)), b)
                     for a, b in K0_INEQS) == base]
    ck("K0_not_invariant_under_any_nonidentity_coordinate_permutation",
       len(fixing) == 1 and fixing[0] == (0, 1, 2, 3),
       "%d of 24 permutations fix the facet set" % len(fixing))


def check_orthant_reduction():
    """x in K_0 iff |x| in K_+.  This is the partition into 2^m congruent orthant
    pieces underlying the case n = m, so it carries weight and is proved exactly;
    the rational grid is only a cross-check.

    Exact proof, from two properties of the row set R of K_0 checked below:
      (closure) (a,b) in R implies (eps.a, b) in R for every sign pattern eps;
      (abs)     (a,b) in R implies (|a|, b) in R.
    If x in K_0 then for any (a,b) in R pick eps_j = sign(x_j); then
    a.|x| = (eps.a).x <= b by closure.  Conversely if |x| in K_0 then for any
    (a,b) in R, a.x <= |a|.|x| <= b by (abs).  So x in K_0 iff |x| in K_0, and
    |x| >= 0, so this is the same as |x| in K_0 cap R_+^4 = K_+ (the last
    equality is the exact row argument in the next function)."""
    rows = set((tuple(a), b) for a, b in K0_INEQS)
    closure = all((tuple(e * ai for e, ai in zip(eps, a)), b) in rows
                  for a, b in K0_INEQS
                  for eps in itertools.product((1, -1), repeat=4))
    absrow = all((tuple(abs(ai) for ai in a), b) in rows for a, b in K0_INEQS)
    ck("orthant_reduction_exact_row_set_closed_under_signs_and_abs",
       closure and absrow,
       "closure %s, abs-row %s over %d rows" % (closure, absrow, len(rows)))
    grid = [F(k, 16) for k in range(-5, 6)] + [F(1, 3), F(-1, 3)]
    tested = 0
    bad = 0
    for x in itertools.product(grid, repeat=4):
        a = satisfies(x, K0_INEQS)
        b = satisfies(tuple(abs(xi) for xi in x), KPLUS_INEQS)
        tested += 1
        if a != b:
            bad += 1
    ck("orthant_reduction_x_in_K0_iff_abs_x_in_Kplus_grid_cross_check", bad == 0,
       "%d grid points, %d mismatches" % (tested, bad))


def check_kplus_is_exactly_K0_cap_orthant():
    """Exact (not grid-sampled) proof that the six inequalities KPLUS_INEQS
    define K_0 cap R_+^4, so that the polygon Q_+ below really is Phi of the
    positive part of the body in the conjecture.

    (K_0 cap R_+^4 subset K_+) every KPLUS row is either one of the four
    -x_j <= 0 rows or a row of K_0, so any x in K_0 with x >= 0 satisfies all
    six of them.
    (K_+ subset K_0 cap R_+^4) the four orthant rows are present, and every row
    a.x <= b of K_0 is dominated on R_+^4 by some KPLUS row a'.x <= b' with
    a'_j >= a_j for all j and b' <= b; then for x >= 0, a.x <= a'.x <= b' <= b."""
    k0rows = set((tuple(a), b) for a, b in K0_INEQS)
    orthant = set((tuple(-1 if j == i else 0 for j in range(4)), 0)
                  for i in range(4))
    kprows = set((tuple(a), b) for a, b in KPLUS_INEQS)
    ck("Kplus_rows_are_K0_rows_or_orthant_rows",
       kprows <= (k0rows | orthant) and orthant <= kprows,
       "%d rows, %d of them orthant" % (len(kprows), len(kprows & orthant)))
    undominated = []
    for a, b in K0_INEQS:
        if not any(all(F(ap[j]) >= F(a[j]) for j in range(4)) and F(bp) <= F(b)
                   for ap, bp in KPLUS_INEQS):
            undominated.append((a, b))
    ck("K0_rows_dominated_on_the_orthant_by_a_Kplus_row", not undominated,
       "%d of %d K0 rows undominated" % (len(undominated), len(K0_INEQS)))


def check_kplus_vertices():
    """Exhaustive census of the 15 four-element active sets of K_+."""
    verts, singular, infeas, feasible_sets = vertices_of(KPLUS_INEQS, 4)
    total = 15
    ck("Kplus_active_set_census_9_feasible_4_singular_2_infeasible",
       feasible_sets == 9 and len(verts) == 9 and singular == 4
       and len(infeas) == 2 and len(set(infeas)) == 2
       and singular + len(infeas) + feasible_sets == total,
       "feasible active sets %d giving %d distinct vertices, singular %d, "
       "infeasible %d, total %d of %d"
       % (feasible_sets, len(verts), singular, len(infeas),
          feasible_sets + singular + len(infeas), total))
    ck("Kplus_infeasible_points_match_paper", set(infeas) == INFEASIBLE_PAPER,
       "%d distinct" % len(set(infeas)))
    ck("Kplus_vertex_set_equals_V_plus", verts == V_PLUS_PAPER,
       "computed %d, paper %d" % (len(verts), len(V_PLUS_PAPER)))
    return verts


def check_sign_hull(verts):
    """K_0 = conv(sign orbit of V_+), the paper's eq. (signhull).  The orbit has
    25 points, of which only 16 are vertices (the origin and the eight points
    with a single nonzero coordinate are not extreme), so the two inclusions are
    NOT "vertex set = orbit": what is checked is
        vert(K_0) subset orbit      and      orbit subset K_0,
    which with K_0 bounded (checked above) gives
        K_0 = conv(vert K_0) subset conv(orbit) subset K_0."""
    orbit = set()
    for v in verts:
        for eps in itertools.product((1, -1), repeat=4):
            orbit.add(tuple(e * vi for e, vi in zip(eps, v)))
    k0verts, _, _, _ = vertices_of(K0_INEQS, 4)
    ck("K0_vertices_are_a_subset_of_the_sign_orbit_of_V_plus",
       k0verts <= orbit, "%d vertices, orbit size %d" % (len(k0verts), len(orbit)))
    ck("sign_orbit_of_V_plus_lies_inside_K0",
       all(satisfies(p, K0_INEQS) for p in orbit),
       "orbit %d points, all feasible" % len(orbit))
    # Extremality certificate: for each computed vertex p, the sum of its
    # active facet normals is a functional strictly maximised at p alone.
    good = 0
    for p in k0verts:
        active = [a for a, b in K0_INEQS
                  if sum(F(ai) * xi for ai, xi in zip(a, p)) == F(b)]
        if rank(active) != 4:
            continue
        c = [sum(F(a[j]) for a in active) for j in range(4)]
        val = sum(cj * pj for cj, pj in zip(c, p))
        if all(sum(cj * qj for cj, qj in zip(c, q)) < val
               for q in k0verts if q != p):
            good += 1
    ck("all_16_K0_vertices_carry_a_strict_supporting_functional",
       good == len(k0verts) == 16, "%d of %d" % (good, len(k0verts)))
    return k0verts


def dot(u, v):
    return sum(F(a) * F(b) for a, b in zip(u, v))


def check_gram():
    G = ((dot(B1, B1), dot(B1, B2)), (dot(B2, B1), dot(B2, B2)))
    ck("gram_matrix_B_transpose_B", G == tuple(tuple(F(v) for v in r)
                                              for r in G_PAPER),
       "G = [[%s,%s],[%s,%s]]" % (G[0][0], G[0][1], G[1][0], G[1][1]))
    det = G[0][0] * G[1][1] - G[0][1] * G[1][0]
    ck("det_gram_146_so_dim_E0_is_2", det == DETG_PAPER and det != 0,
       "det G = %s" % det)
    return det


def check_chart_factor(det):
    """area(Phi(A)) = sqrt(det G) * |P_{E0} A|_2.  Verified squared, so that
    all arithmetic stays rational: for triangles with vertices B c_i inside
    E_0, the Euclidean area uses the Gram determinant of the genuine
    4-dimensional edge vectors, while the image area uses Phi."""
    triples = [((0, 0), (1, 0), (0, 1)), ((1, 2), (3, -1), (-2, 5)),
               ((-3, 1), (2, 2), (0, -4)), ((5, 0), (0, 7), (-1, -1))]
    bad = 0
    for c0, c1, c2 in triples:
        u = [F(c1[0] - c0[0]) * F(B1[i]) + F(c1[1] - c0[1]) * F(B2[i])
             for i in range(4)]
        v = [F(c2[0] - c0[0]) * F(B1[i]) + F(c2[1] - c0[1]) * F(B2[i])
             for i in range(4)]
        # squared Euclidean area of the triangle inside E_0
        areaE_sq = (dot(u, u) * dot(v, v) - dot(u, v) ** 2) / 4
        w1 = phi(u)
        w2 = phi(v)
        areaP = abs(w1[0] * w2[1] - w1[1] * w2[0]) / 2
        if areaP ** 2 != F(det) * areaE_sq or areaE_sq <= 0:
            bad += 1
    ck("chart_area_factor_sqrt_146_on_4_triangles", bad == 0,
       "%d triangles, %d mismatches" % (len(triples), bad))


def check_phi_kills_orthogonal_complement():
    """Phi = Phi o P_{E0}: Phi vanishes on E_0^perp, so the planar areas below
    really are areas of projections onto E_0.

    (The two vectors are ASSERTED to be orthogonal to b1 and b2 rather than
    filtered on that condition: filtering first and then testing phi(z) = 0
    would be the same two equations twice, a clause that cannot fail.  With
    det G != 0 the rank of phi is 2, so ker phi is 2-dimensional and exhibiting
    a rank-2 subspace of it identifies ker phi = E_0^perp exactly.)"""
    basis = [(1, 3, 0, 0), (0, -4, 2, -3)]
    orthogonal = all(dot(B1, z) == 0 and dot(B2, z) == 0 for z in basis)
    ok = (orthogonal and rank(basis) == 2
          and all(phi(z) == (F(0), F(0)) for z in basis))
    ck("phi_vanishes_on_a_rank_2_orthogonal_complement", ok,
       "%d vectors, orthogonal to b1,b2: %s, rank %d"
       % (len(basis), orthogonal, rank(basis)))


def check_polygon(tag, gens, hineqs, paper_pts, paper_area):
    """gens: the projected generators whose hull is the polygon.  Show
    conv(gens) = {y : hineqs} by two inclusions, then take the area.

    The two inclusions are conv(gens) subset H (every generator satisfies every
    half-space) and H subset conv(gens) (every vertex of H is a generator).
    The second one needs H = conv(vert H), which is FALSE for an unbounded H:
    the half-strip {t >= 0, t <= 2, s >= 0, s >= t-1} has exactly the three
    vertices (0,0), (0,1), (1,2), so a generator such as (10,1) would sit
    inside H and outside the triangle they span, and the shoelace area below
    would be wrong.  Hence H must be certified bounded first."""
    ck(tag + "_halfspace_certificate_is_bounded",
       cone_trivial([a for a, _ in hineqs], 2),
       "%d half-spaces, recession cone trivial" % len(hineqs))
    ck(tag + "_projected_generators_satisfy_halfspaces",
       all(satisfies(g, hineqs) for g in gens),
       "%d generators" % len(gens))
    hv, _, _, _ = vertices_of(hineqs, 2)
    ck(tag + "_halfspace_vertices_equal_paper_vertex_list",
       hv == set(paper_pts),
       "computed %d, paper %d" % (len(hv), len(paper_pts)))
    ck(tag + "_every_halfspace_vertex_is_a_projected_generator",
       hv <= set(gens), "%d vertices" % len(hv))
    if not ck(tag + "_has_at_least_three_vertices_so_the_area_is_defined",
              len(hv) >= 3, "%d vertices" % len(hv)):
        ck(tag + "_paper_vertex_order_is_counterclockwise", False, "degenerate")
        ck(tag + "_shoelace_area", False, "degenerate")
        return F(0), frozenset(hv)
    ordered = ccw_sorted(list(hv))
    twice = shoelace2(ordered)
    ck(tag + "_paper_vertex_order_is_counterclockwise",
       shoelace2(paper_pts) == twice and twice > 0,
       "2*area = %s" % twice)
    ck(tag + "_shoelace_area", twice / 2 == paper_area,
       "area = %s" % (twice / 2))
    return twice / 2, frozenset(hv)


def check_nesting(qplus_verts, q_verts):
    """Q_+ = Phi(K_+) is contained in Q = Phi(K_0), as it must be.  The vertex
    sets used here are the ones DERIVED from the half-space certificates, not
    the paper's lists, so this check does not lean on the paper's data even
    though the two were already shown to agree."""
    ck("Qplus_contained_in_Q", all(satisfies(p, HQ) for p in qplus_verts),
       "%d derived vertices of Q_+ inside Q" % len(qplus_verts))
    strict = [p for p in q_verts if not satisfies(p, HQPLUS)]
    ck("Q_strictly_larger_than_Qplus", len(strict) > 0,
       "%d of %d derived vertices of Q outside Q_+"
       % (len(strict), len(q_verts)))


def check_ratio(a_plus, a_full):
    """The chart factor sqrt(146) cancels in the ratio, so the planar areas
    give the projection ratio exactly."""
    if not ck("Qplus_area_is_strictly_positive", a_plus > 0,
              "area(Q_+) = %s" % a_plus):
        for nm in ("ratio_equals_1399_over_336",
                   "ratio_decomposition_4_plus_55_over_336",
                   "conjecture_fails_for_m4_n2"):
            ck(nm, False, "no ratio: area(Q_+) is not positive")
        return F(0)
    r = a_full / a_plus
    ck("ratio_equals_1399_over_336", r == RATIO_PAPER, "ratio = %s" % r)
    ck("ratio_decomposition_4_plus_55_over_336", r == 4 + F(55, 336),
       "%s" % (r - 4))
    ck("conjecture_fails_for_m4_n2", r > F(2) ** 2,
       "%s > 4, excess %s" % (r, r - 4))
    return r


def check_product_extension(r):
    """K = K_0 x [-1,1]^{m-4}, E = E_0 + R^{n-2}: the ratio scales by 2^{n-2}
    because P_E K gains a factor 2^{n-2} and P_E(K cap R_+^m) gains 1."""
    pairs = [(m, n) for m in range(4, 21) for n in range(2, m - 1)]
    bad = [(m, n) for m, n in pairs
           if not (F(2) ** (n - 2) * r > F(2) ** n)]
    ck("product_extension_beats_2_to_the_n_for_all_pairs", not bad,
       "%d pairs (m<=20), %d failures" % (len(pairs), len(bad)))
    # The margin depends only on n, and "grow" is asserted, not assumed:
    # strict monotonicity in n is tested, not merely positivity.
    margins = {}
    for _, n in pairs:
        margins[n] = F(2) ** (n - 2) * r - F(2) ** n
    ns = sorted(margins)
    grows = all(margins[ns[i + 1]] > margins[ns[i]] for i in range(len(ns) - 1))
    ck("product_extension_margins_are_positive_and_grow",
       all(x > 0 for x in margins.values())
       and min(margins.values()) == F(55, 336) and grows and len(ns) >= 2,
       "smallest excess %s at n=%d, strictly increasing over n=%d..%d: %s"
       % (min(margins.values()), ns[0], ns[0], ns[-1], grows))

def check_product_body_m5(a_plus, a_full, qverts, orbit):
    """The case m = 5, n = 3: K = K_0 x [-1,1] is unconditional,
    K cap R_+^5 = K_+ x [0,1], and in the chart Psi(x) = (Phi(x_1..x_4), x_5) the
    projection onto E = E_0 + R e_5 is Q x [-1,1] with positive part Q_+ x [0,1].

    This is the ONLY (m,n) pair beyond (4,2) at which the product identities of
    the paper's Theorem are checked as sets; the general statement is the
    arithmetic 2^{n-2}r > 2^n above.  Strength of each part:
      * unconditionality and the positive part are checked on a rational grid
        (all 32 sign patterns), not proved;
      * Psi(K) subset Q x [-1,1] is checked on the grid here, but it is ALSO an
        exact consequence of a check already made above: Q_projected_generators_
        satisfy_halfspaces says Phi(v) in H_Q for every v in the sign orbit, and
        K = K_0 x [-1,1] = conv(orbit x {-1,1}) with Psi linear, so
        Psi(K) = conv{(Phi(v),t)} subset H_Q x [-1,1] = Q x [-1,1];
      * the reverse inclusion is exact: all 2*|vert Q| vertices of Q x [-1,1]
        are attained by explicit preimages, and Psi(K) is convex.
    So the name "_is_Q_times_interval" is earned, but only via those two earlier
    exact facts -- the grid alone would give only one inclusion."""
    grid = [F(k, 3) for k in range(-3, 4)]
    in_K = lambda x: satisfies(x[:4], K0_INEQS) and abs(x[4]) <= 1
    bad_unc = bad_orth = 0
    pts = list(itertools.product(grid, repeat=5))
    # The grid is symmetric about 0, so every sign image of a grid point is a
    # grid point; memoise membership and then test ALL 32 sign patterns.  A
    # hand-picked handful of patterns would not test unconditionality, only
    # invariance under those particular reflections.
    member = {}
    for x in pts:
        member[x] = in_K(x)
    signs = list(itertools.product((1, -1), repeat=5))
    for x in pts:
        mx = member[x]
        for eps in signs:
            if mx != member[tuple(e * xi for e, xi in zip(eps, x))]:
                bad_unc += 1
        pos = mx and all(xi >= 0 for xi in x)
        if pos != (satisfies(x[:4], KPLUS_INEQS) and 0 <= x[4] <= 1):
            bad_orth += 1
    ck("product_body_m5_is_unconditional", bad_unc == 0 and len(signs) == 32,
       "%d grid points x %d sign patterns, %d mismatches"
       % (len(pts), len(signs), bad_unc))
    ck("product_body_m5_positive_part_is_Kplus_times_unit_interval",
       bad_orth == 0, "%d grid points, %d mismatches" % (len(pts), bad_orth))
    img_ok = all(satisfies(phi(x[:4]), HQ) and abs(x[4]) <= 1
                 for x in pts if member[x])
    # Every vertex of Q x [-1,1] must be ATTAINED by an actual point of
    # K = K_0 x [-1,1] under Psi; counting 2*len(qverts) alone would merely
    # restate how many vertices Q has, so each one is realised by exhibiting a
    # preimage among the sign-orbit generators of K_0.
    attained = set()
    for v in qverts:
        pre = [x for x in orbit if phi(x) == v]
        if not pre:
            continue
        for last in (F(-1), F(1)):
            x5 = tuple(list(pre[0]) + [last])
            if in_K(x5) and (phi(x5[:4]), x5[4]) == (v, last):
                attained.add((v[0], v[1], last))
    ck("product_body_m5_chart_image_is_Q_times_interval",
       img_ok and len(attained) == 2 * len(qverts) and len(qverts) >= 3,
       "grid images inside Q x [-1,1]: %s (the exact inclusion is the earlier "
       "generator check plus linearity); %d of %d product vertices attained "
       "exactly" % (img_ok, len(attained), 2 * len(qverts)))
    if a_plus <= 0:
        ck("product_body_m5_ratio_exceeds_2_cubed", False,
           "no ratio: area(Q_+) is not positive")
        return
    r3 = (F(2) * a_full) / (F(1) * a_plus)
    ck("product_body_m5_ratio_exceeds_2_cubed", r3 > 8 and r3 == F(1399, 168),
       "ratio = %s > 8" % r3)


def rank(rows):
    M = [[F(v) for v in r] for r in rows]
    r = 0
    ncols = len(M[0]) if M else 0
    for c in range(ncols):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [v / pv for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [M[i][j] - f * M[r][j] for j in range(ncols)]
        r += 1
    return r


if __name__ == "__main__":
    check_encoding_matches_the_paper()
    check_body_is_a_body()
    check_unconditional()
    check_orthant_reduction()
    check_kplus_is_exactly_K0_cap_orthant()
    VP = check_kplus_vertices()
    check_sign_hull(VP)
    detG = check_gram()
    check_chart_factor(detG)
    check_phi_kills_orthogonal_complement()
    gens_plus = [phi(v) for v in sorted(VP)]
    orbit = set()
    for v in VP:
        for eps in itertools.product((1, -1), repeat=4):
            orbit.add(tuple(e * vi for e, vi in zip(eps, v)))
    gens_full = [phi(v) for v in sorted(orbit)]
    a_plus, qplus_verts = check_polygon("Qplus", gens_plus, HQPLUS,
                                        QPLUS_PAPER, AREA_QPLUS_PAPER)
    a_full, q_verts = check_polygon("Q", gens_full, HQ, Q_PAPER, AREA_Q_PAPER)
    check_nesting(qplus_verts, q_verts)
    ratio = check_ratio(a_plus, a_full)
    check_product_extension(ratio)
    check_product_body_m5(a_plus, a_full, q_verts, orbit)
    print("NOT RE-RUN HERE (1): the positive cases n = 1 and n = m-1 of the "
          "classification are quoted from the cited work and are not "
          "recomputed.")
    print("NOT RE-RUN HERE (2): the case n = m is supported only through the "
          "orthant reduction above, not by a 4-dimensional volume "
          "computation.")
    print("NOT RE-RUN HERE (3): the product identities "
          "P_E K = P_{E_0}K_0 x [-1,1]^{n-2} and "
          "P_E(K cap R_+^m) = P_{E_0}K_+ x [0,1]^{n-2} are verified as sets "
          "only for the single pair (m,n) = (5,3), and there only on a "
          "rational grid; the all-(m,n) statement above is the arithmetic "
          "consequence 2^{n-2}*(1399/336) > 2^n, enumerated for m <= 20 only.")
    print("NOT RE-RUN HERE (4): checks whose detail line says 'grid' are "
          "finite samples, not proofs over R^m; the load-bearing orthant "
          "reduction and the K_0 encoding also have exact, sample-free "
          "counterparts above.")
    print("NOT RE-RUN HERE (5): the hypothesis class is taken from the paper. "
          "Proved above: K_0 is invariant under all 16 coordinate sign changes "
          "(hence unconditional) and under no nonidentity coordinate "
          "permutation, so it satisfies the sign-invariance hypothesis and not "
          "the stronger permutation-invariant one. That sign invariance is the "
          "whole hypothesis of the conjecture being refuted is transcribed "
          "from the cited work, which is not part of this material and is not "
          "read by this program.")
    print("NOT RE-RUN HERE (6): no computation here touches bibliographic "
          "data -- the preprint identifier and version, the author names, and "
          "the pinpoint citations to the two propositions that supply the "
          "positive cases n = 1 and n = m-1. The negative half of the "
          "if-and-only-if classification (every m >= 4, 2 <= n <= m-2) is "
          "recomputed above; the positive half rests on those citations "
          "together with the elementary case n = m.")
    finish()

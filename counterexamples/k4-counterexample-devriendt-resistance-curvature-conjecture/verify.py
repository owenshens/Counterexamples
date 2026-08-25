#!/usr/bin/env python3
"""Verification of a K_4 counterexample to a resistance-curvature conjecture.

Setting.  G = K_4 on V = {0,1,2,3}, edges ordered 01,02,03,12,13,23, with
positive conductances c.  L is the weighted Laplacian, L^+ its Moore-Penrose
pseudoinverse, Omega_uv = (e_u - e_v)^T L^+ (e_u - e_v) the effective
resistance matrix and K = Omega^{-1}.  The weighting is NORMALIZED when
1^T K 1 = 1.  With r_uv = c_uv * Omega_uv, the resistance curvature is
p_v = 1 - (1/2) sum_{u : uv in E} r_uv, and the resistance capacity is
tau_{} = 0, tau_{v} = 1/2, tau_U = 1/2 + 1/2 * (1^T Omega[U]^{-1} 1)^{-1}
for |U| >= 2.  The conjecture under test asserts: tau submodular ==> p >= 0.

TAKEN FROM THE PAPER (inputs; transcribed verbatim in the block below,
never used as a substitute for a computation):
  * the integer weight vector w = (1,1,2,1,4,4) and the normalizing scalar
    3485/15876, which together define the exhibited conductance vector;
  * the displayed Laplacian L(w) and resistance matrix (1/126)*Omega(w);
  * the scalar 1^T Omega(w)^{-1} 1 = 15876/3485;
  * the relative resistances (22/63, 22/63, 5/9, 2/7, 46/63, 46/63);
  * the curvature vector (47/126, 20/63, 20/63, -1/126);
  * the count 24 of elementary submodularity inequalities, the 15 tabulated
    orbit representatives, and the minimum slack 1123/69700;
  * the remark's uniform-K_4 curvature (1/4,1/4,1/4,1/4).

DERIVED HERE (computed by exact rational arithmetic from w alone):
  * L(w), L^+ and the Moore-Penrose identities it must satisfy;
  * Omega(w), the normalizing scalar, and hence the exhibited c;
  * verification that c satisfies the HYPOTHESES of the conjecture, namely
    normalization and submodularity of tau (all 24 slacks, strictly);
  * verification that c VIOLATES the CONCLUSION, p_3 < 0 -- computed from
    the weights, not asserted;
  * the 15 tabulated slacks, the minimum slack and its location, and the
    orbit structure under the vertex involution 1 <-> 2;
  * the minimality remark: curvature is nonnegative on every weighted K_2,
    P_3 and K_3 in an integer range, and the triangle closed form;
  * stability under a grid of rational perturbations;
  * censuses over integer weightings of K_4 (see the closing NOT RE-RUN line
    for exactly what is and is not covered).

Standard library only; every decision is an exact rational comparison.
"""

from fractions import Fraction as F
from math import gcd
import itertools
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


def finish():
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


def mat_inv(M):
    """Exact inverse of a square rational matrix by Gauss-Jordan elimination.

    Raises ValueError if the matrix is singular.
    """
    n = len(M)
    A = [[F(M[i][j]) for j in range(n)] + [F(1 if i == j else 0) for j in range(n)]
         for i in range(n)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            raise ValueError("singular matrix")
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return [row[n:] for row in A]


def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)]
            for i in range(n)]


def ones_quad(Minv):
    """1^T Minv 1 for a square rational matrix."""
    return sum(sum(row) for row in Minv)


def laplacian(n, edges, weights):
    L = [[F(0)] * n for _ in range(n)]
    for (u, v), c in zip(edges, weights):
        c = F(c)
        L[u][u] += c
        L[v][v] += c
        L[u][v] -= c
        L[v][u] -= c
    return L


def pinv_laplacian(L):
    """Moore-Penrose pseudoinverse of a connected weighted Laplacian:
    (L + J/n)^{-1} - J/n, with J the all-ones matrix."""
    n = len(L)
    q = F(1, n)
    shifted = [[L[i][j] + q for j in range(n)] for i in range(n)]
    Inv = mat_inv(shifted)
    return [[Inv[i][j] - q for j in range(n)] for i in range(n)]


def resistance_matrix(L):
    n = len(L)
    P = pinv_laplacian(L)
    return [[P[i][i] + P[j][j] - 2 * P[i][j] for j in range(n)] for i in range(n)]


def submatrix(M, U):
    return [[M[i][j] for j in U] for i in U]


def tau(Om, U):
    """Resistance capacity of the vertex subset U (tuple of indices)."""
    if len(U) == 0:
        return F(0)
    if len(U) == 1:
        return F(1, 2)
    s = ones_quad(mat_inv(submatrix(Om, U)))
    return F(1, 2) + F(1, 2) / s


def curvature(n, edges, weights, Om):
    """p_v = 1 - (1/2) sum_{u: uv in E} c_uv * omega_uv."""
    p = [F(1)] * n
    for (u, v), c in zip(edges, weights):
        r = F(c) * Om[u][v]
        p[u] -= r / 2
        p[v] -= r / 2
    return p


def elementary_slacks(n, Om):
    """All Delta_ij(S) = tau(S+i) + tau(S+j) - tau(S) - tau(S+i+j),
    for i<j and S ranging over subsets of V minus {i,j}.
    Returns a dict keyed by ((i,j), S)."""
    V = list(range(n))
    out = {}
    for i, j in itertools.combinations(V, 2):
        rest = [v for v in V if v not in (i, j)]
        for k in range(len(rest) + 1):
            for S in itertools.combinations(rest, k):
                Si = tuple(sorted(S + (i,)))
                Sj = tuple(sorted(S + (j,)))
                Sij = tuple(sorted(S + (i, j)))
                out[((i, j), S)] = (tau(Om, Si) + tau(Om, Sj)
                                    - tau(Om, S) - tau(Om, Sij))
    return out


def normalize(n, edges, weights):
    """Rescale conductances so that 1^T Omega^{-1} 1 = 1, and return
    (scaled weights, Omega of the scaled weighting)."""
    Om = resistance_matrix(laplacian(n, edges, weights))
    s = ones_quad(mat_inv(Om))
    # Omega(lam*w) = Omega(w)/lam, so 1^T Omega(lam w)^{-1} 1 = lam * s.
    lam = F(1) / s
    w2 = [F(x) * lam for x in weights]
    Om2 = [[x / lam for x in row] for row in Om]
    return w2, Om2


##  Values transcribed from the paper (inputs).  ------------------------------
N = 4
EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
W_PAPER = [1, 1, 2, 1, 4, 4]
SCALE_PAPER = F(3485, 15876)
L_PAPER = [[4, -1, -1, -2], [-1, 6, -1, -4], [-1, -1, 6, -4], [-2, -4, -4, 10]]
OMEGA_NUM_PAPER = [[0, 44, 44, 35], [44, 0, 36, 23],
                   [44, 36, 0, 23], [35, 23, 23, 0]]
OMEGA_DEN_PAPER = 126
ONES_QUAD_PAPER = F(15876, 3485)
R_PAPER = [F(22, 63), F(22, 63), F(5, 9), F(2, 7), F(46, 63), F(46, 63)]
P_PAPER = [F(47, 126), F(20, 63), F(20, 63), F(-1, 126)]
MIN_SLACK_PAPER = F(1123, 69700)
MIN_SLACK_ARGMIN_PAPER = ((0, 3), (1, 2))
N_INEQ_PAPER = 24
TABLE_PAPER = {
    ((0, 1), ()): F(713, 6970), ((0, 1), (2,)): F(3888, 17425),
    ((0, 1), (3,)): F(2107, 20910), ((0, 1), (2, 3)): F(2389, 41820),
    ((0, 3), ()): F(128, 697), ((0, 3), (1,)): F(112, 615),
    ((0, 3), (1, 2)): F(1123, 69700),
    ((1, 2), ()): F(1217, 6970), ((1, 2), (0,)): F(5148, 17425),
    ((1, 2), (3,)): F(207, 2788), ((1, 2), (0, 3)): F(64, 2091),
    ((1, 3), ()): F(1018, 3485), ((1, 3), (0,)): F(3038, 10455),
    ((1, 3), (2,)): F(2673, 13940), ((1, 3), (0, 2)): F(1346, 52275),
}
UNIFORM_K4_CURVATURE_PAPER = [F(1, 4)] * 4


def check_object_well_formed():
    ok = (N == 4
          and EDGES == sorted(itertools.combinations(range(N), 2))
          and len(EDGES) == 6
          and len(W_PAPER) == 6
          and all(F(x) > 0 for x in W_PAPER))
    ck("object_is_a_positively_weighted_K4", ok,
       "V=%d, E=%d, w=%s" % (N, len(EDGES), tuple(W_PAPER)))


def check_laplacian():
    L = laplacian(N, EDGES, W_PAPER)
    ok = L == [[F(x) for x in row] for row in L_PAPER]
    rowsums_zero = all(sum(row) == 0 for row in L)
    ck("weighted_laplacian_matches_paper", ok and rowsums_zero,
       "rows=%s" % ([[str(x) for x in row] for row in L],))


def check_pseudoinverse_identities():
    L = laplacian(N, EDGES, W_PAPER)
    P = pinv_laplacian(L)
    sym = all(P[i][j] == P[j][i] for i in range(N) for j in range(N))
    kern = all(sum(P[i]) == 0 for i in range(N))
    LPL = mat_mul(mat_mul(L, P), L)
    PLP = mat_mul(mat_mul(P, L), P)
    ok = (sym and kern and LPL == L and PLP == P)
    ck("pseudoinverse_satisfies_moore_penrose_identities", ok,
       "L*Ld*L=L, Ld*L*Ld=Ld, symmetric, Ld*1=0")


def check_omega():
    Om = resistance_matrix(laplacian(N, EDGES, W_PAPER))
    target = [[F(OMEGA_NUM_PAPER[i][j], OMEGA_DEN_PAPER) for j in range(N)]
              for i in range(N)]
    diag0 = all(Om[i][i] == 0 for i in range(N))
    pos = all(Om[i][j] > 0 for i in range(N) for j in range(N) if i != j)
    ck("resistance_matrix_matches_paper", Om == target and diag0 and pos,
       "126*Omega=%s" % ([[str(126 * x) for x in row] for row in Om],))


def check_normalizing_scalar():
    Om = resistance_matrix(laplacian(N, EDGES, W_PAPER))
    s = ones_quad(mat_inv(Om))
    lam = F(1) / s
    ck("normalizing_scalar_matches_paper",
       s == ONES_QUAD_PAPER and lam == SCALE_PAPER,
       "1^T Omega(w)^-1 1 = %s, lambda = %s" % (s, lam))


def exhibited():
    """The conductance vector the paper exhibits, built from the paper's own
    scalar and weight vector, together with its resistance matrix."""
    c = [SCALE_PAPER * F(x) for x in W_PAPER]
    Om = resistance_matrix(laplacian(N, EDGES, c))
    return c, Om


def check_hypothesis_normalized():
    c, Om = exhibited()
    K = mat_inv(Om)
    s = ones_quad(K)
    tv = tau(Om, tuple(range(N)))
    single = all(tau(Om, (v,)) == F(1, 2) for v in range(N))
    empty = tau(Om, ()) == 0
    ck("hypothesis_normalized_holds_for_exhibited_c",
       s == 1 and tv == 1 and single and empty,
       "1^T K 1 = %s, tau(V) = %s" % (s, tv))


def check_relative_resistances():
    c, Om = exhibited()
    r = [F(ci) * Om[u][v] for ci, (u, v) in zip(c, EDGES)]
    Omw = resistance_matrix(laplacian(N, EDGES, W_PAPER))
    rw = [F(wi) * Omw[u][v] for wi, (u, v) in zip(W_PAPER, EDGES)]
    ck("relative_resistances_match_paper_and_are_scale_invariant",
       r == R_PAPER and rw == R_PAPER,
       "r = %s" % ([str(x) for x in r],))


def check_curvature_vector():
    c, Om = exhibited()
    p = curvature(N, EDGES, c, Om)
    ck("curvature_vector_matches_paper", p == P_PAPER,
       "p = %s" % ([str(x) for x in p],))


def check_conclusion_violated():
    """Load bearing: the conclusion p >= 0 of the conjecture must FAIL."""
    c, Om = exhibited()
    p = curvature(N, EDGES, c, Om)
    neg = [v for v in range(N) if p[v] < 0]
    ck("conclusion_nonnegative_curvature_is_violated",
       len(neg) == 1 and neg == [3] and p[3] < 0 and p[3] == P_PAPER[3],
       "negative coordinates %s, p_3 = %s" % (neg, p[3]))


def check_hypothesis_submodular():
    """Load bearing: the hypothesis of the conjecture must HOLD, strictly."""
    c, Om = exhibited()
    sl = elementary_slacks(N, Om)
    nonpos = sorted(k for k, v in sl.items() if v <= 0)
    ck("hypothesis_all_24_elementary_slacks_strictly_positive",
       len(sl) == N_INEQ_PAPER and not nonpos,
       "%d inequalities, %d nonpositive" % (len(sl), len(nonpos)))


def check_slack_table():
    c, Om = exhibited()
    sl = elementary_slacks(N, Om)
    bad = []
    for key, want in TABLE_PAPER.items():
        got = sl.get(key)
        if got != want:
            bad.append((key, str(want), str(got)))
    ck("published_table_of_15_orbit_representatives_is_exact",
       len(TABLE_PAPER) == 15 and not bad,
       "15 published values reproduced" if not bad else "mismatches %s" % (bad,))


def check_min_slack():
    c, Om = exhibited()
    sl = elementary_slacks(N, Om)
    m = min(sl.values())
    argmins = sorted(k for k, v in sl.items() if v == m)
    detail = "min = %s at Delta_{%d%d}(%s)" % ((m,) + argmins[0][0]
                                               + (argmins[0][1],))
    if m != MIN_SLACK_PAPER or argmins != [MIN_SLACK_ARGMIN_PAPER]:
        detail += "; paper states %s at %s" % (MIN_SLACK_PAPER,
                                               MIN_SLACK_ARGMIN_PAPER)
    ck("minimum_slack_value_and_location_match_paper",
       m == MIN_SLACK_PAPER and argmins == [MIN_SLACK_ARGMIN_PAPER] and m > 0,
       detail)


def apply_perm(perm, weights):
    """Push a vertex permutation forward to the edge-weight vector."""
    d = {}
    for (u, v), x in zip(EDGES, weights):
        a, b = sorted((perm[u], perm[v]))
        d[(a, b)] = x
    return [d[e] for e in EDGES]


def check_involution_orbits():
    c, Om = exhibited()
    sl = elementary_slacks(N, Om)
    sigma = {0: 0, 1: 2, 2: 1, 3: 3}
    fixes_c = apply_perm(sigma, c) == c
    bad = []
    orbits = set()
    for (i, j), S in sl:
        ii, jj = sorted((sigma[i], sigma[j]))
        SS = tuple(sorted(sigma[v] for v in S))
        if sl[((i, j), S)] != sl[((ii, jj), SS)]:
            bad.append(((i, j), S))
        orbits.add(min(((i, j), S), ((ii, jj), SS)))
    reps_agree = orbits == set(TABLE_PAPER.keys())
    ck("involution_1_2_is_a_symmetry_with_15_orbits",
       fixes_c and not bad and len(orbits) == 15 and reps_agree,
       "sigma fixes c, %d orbits, representatives equal the published table"
       % len(orbits))


def check_uniform_k4():
    w = [1] * 6
    c, Om = normalize(N, EDGES, w)
    p = curvature(N, EDGES, c, Om)
    sl = elementary_slacks(N, Om)
    ck("uniform_k4_normalized_curvature_is_one_quarter",
       p == UNIFORM_K4_CURVATURE_PAPER and min(sl.values()) >= 0,
       "p = %s, min slack = %s" % ([str(x) for x in p], min(sl.values())))


SMALL_GRAPHS = [("K2", 2, [(0, 1)]),
                ("P3", 3, [(0, 1), (1, 2)]),
                ("K3", 3, [(0, 1), (0, 2), (1, 2)])]


def check_small_order_minimality():
    """No connected simple graph on fewer than four vertices admits a negative
    resistance curvature, so no counterexample of order < 4 can exist."""
    worst = None
    total = 0
    where = None
    for name, n, edges in SMALL_GRAPHS:
        for w in itertools.product(range(1, 7), repeat=len(edges)):
            total += 1
            Om = resistance_matrix(laplacian(n, edges, w))
            p = curvature(n, edges, w, Om)
            m = min(p)
            if worst is None or m < worst:
                worst, where = m, (name, w)
    ck("no_negative_curvature_on_graphs_of_order_below_four",
       worst is not None and worst >= 0,
       "%d weightings of K2,P3,K3 (conductances 1..6); min curvature %s at %s"
       % (total, worst, where))


def check_triangle_closed_form():
    """Paper's series-parallel formula p_v = (x_e+x_f)/(2X), x_h = 1/c_h."""
    n, edges = 3, [(0, 1), (0, 2), (1, 2)]
    inc = {v: tuple(h for h, e in enumerate(edges) if v in e) for v in range(n)}
    bad = 0
    total = 0
    for w in itertools.product(range(1, 8), repeat=3):
        Om = resistance_matrix(laplacian(n, edges, w))
        p = curvature(n, edges, w, Om)
        x = [F(1, t) for t in w]
        X = sum(x)
        for v in range(3):
            i, j = inc[v]
            total += 1
            if p[v] != (x[i] + x[j]) / (2 * X) or p[v] <= 0:
                bad += 1
    ck("triangle_curvature_matches_series_parallel_closed_form", bad == 0,
       "%d (weighting, vertex) pairs agree and are positive" % total)


def check_tree_relative_resistances():
    bad = 0
    total = 0
    for name, n, edges in SMALL_GRAPHS[:2]:
        for w in itertools.product(range(1, 7), repeat=len(edges)):
            Om = resistance_matrix(laplacian(n, edges, w))
            for wi, (u, v) in zip(w, edges):
                total += 1
                if F(wi) * Om[u][v] != 1:
                    bad += 1
    ck("tree_relative_resistances_equal_one", bad == 0,
       "%d tree edges have c_uv * omega_uv = 1" % total)


def is_counterexample(w):
    """(negative curvature, strictly submodular) for the normalization of w."""
    c, Om = normalize(N, EDGES, w)
    p = curvature(N, EDGES, c, Om)
    if min(p) >= 0:
        return False, min(p), None
    sl = elementary_slacks(N, Om)
    m = min(sl.values())
    return m > 0, min(p), m


def check_perturbation_stability():
    step = F(1, 64)
    bad = []
    worst_p, worst_s, total = None, None, 0
    for d in itertools.product([-step, F(0), step], repeat=6):
        w = [F(a) + b for a, b in zip(W_PAPER, d)]
        total += 1
        c, Om = normalize(N, EDGES, w)
        p = curvature(N, EDGES, c, Om)
        m = min(elementary_slacks(N, Om).values())
        if not (p[3] < 0 and m > 0):
            bad.append(tuple(str(x) for x in w))
        if worst_p is None or p[3] > worst_p:
            worst_p = p[3]
        if worst_s is None or m < worst_s:
            worst_s = m
    ck("counterexample_persists_under_perturbation_by_one_sixtyfourth",
       total == 729 and not bad,
       "%d corner perturbations, %d failing: sup p_3 = %s (need < 0), "
       "inf min-slack = %s (need > 0)" % (total, len(bad), worst_p, worst_s))


def primitive_weights(bound, sum_cap=None):
    for w in itertools.product(range(1, bound + 1), repeat=6):
        if sum_cap is not None and sum(w) > sum_cap:
            continue
        g = 0
        for x in w:
            g = gcd(g, x)
        if g == 1:
            yield w


def check_proved_direction_census():
    """Devriendt's proved implication p >= 0 ==> submodular must never fail;
    a failure would indicate the tau/Delta implementation is wrong."""
    viol = []
    total, ces = 0, 0
    for w in primitive_weights(4):
        total += 1
        c, Om = normalize(N, EDGES, w)
        p = curvature(N, EDGES, c, Om)
        m = min(elementary_slacks(N, Om).values())
        if min(p) >= 0 and m < 0:
            viol.append(w)
        if min(p) < 0 and m > 0:
            ces += 1
    ck("census_primitive_weights_to_four_respects_the_proved_direction",
       total == 4031 and not viol and ces > 0,
       "%d weightings, %d violations of p>=0 ==> submodular, %d counterexamples"
       % (total, len(viol), ces))


def check_minimal_weight_sum():
    sums = []
    total = 0
    for w in primitive_weights(13, sum_cap=13):
        total += 1
        good, _, _ = is_counterexample(w)
        if good:
            sums.append(sum(w))
    ck("no_counterexample_of_smaller_conductance_weight_sum",
       sums and min(sums) == 13 and sum(W_PAPER) == 13,
       "%d primitive weightings with sum<=13; %d counterexamples, all of sum %d"
       % (total, len(sums), min(sums) if sums else -1))


def check_automorphism_orbit():
    orbit = set()
    for perm in itertools.permutations(range(N)):
        orbit.add(tuple(apply_perm(dict(enumerate(perm)), W_PAPER)))
    bad = [w for w in sorted(orbit) if not is_counterexample(w)[0]]
    ck("every_member_of_the_symmetric_group_orbit_is_a_counterexample",
       len(orbit) == 12 and not bad,
       "orbit size %d, all strictly submodular with negative curvature"
       % len(orbit))


def main():
    print("Setting: K_4, edges ordered 01,02,03,12,13,23; "
          "conductances c = (3485/15876)*(1,1,2,1,4,4).")
    suite = [
        check_object_well_formed, check_laplacian,
        check_pseudoinverse_identities, check_omega,
        check_normalizing_scalar, check_hypothesis_normalized,
        check_relative_resistances, check_curvature_vector,
        check_conclusion_violated, check_hypothesis_submodular,
        check_slack_table, check_min_slack, check_involution_orbits,
        check_uniform_k4, check_small_order_minimality,
        check_triangle_closed_form, check_tree_relative_resistances,
        check_perturbation_stability, check_proved_direction_census,
        check_minimal_weight_sum, check_automorphism_orbit,
    ]
    for fn in suite:
        marker = len(CHECKS)
        try:
            fn()
        except Exception as exc:                       # noqa: BLE001
            if len(CHECKS) == marker:
                ck(fn.__name__, False, "raised %s" % (exc,))
        if len(CHECKS) == marker:
            ck(fn.__name__, False, "check produced no result")
    print("NOT RE-RUN: the cited article's own theorems are not reproved here; "
          "only the implication p>=0 ==> submodularity is spot-checked, on the "
          "4031 primitive integer weightings of K_4 with conductances at most "
          "4. Searches over conductances are restricted to the integer ranges "
          "named in each check, and stability is tested on a finite grid of "
          "perturbations rather than an open neighbourhood. The paper's two "
          "remarks are proved there analytically and are only spot-checked "
          "here: minimality of order four on the 258 integer weightings of "
          "K_2, P_3 and K_3 with conductances at most 6, together with 1029 "
          "triangle instances against the series-parallel closed form; and "
          "persistence under perturbation only at the 729 corners of the box "
          "w +/- 1/64, not on any open neighbourhood.")
    finish()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verification of a counterexample to a conjecture on resistance curvature.

The conjecture under attack states: for a normalized weighted graph, submodularity
of the resistance capacity tau implies that the resistance curvature p is
nonnegative.  The paper exhibits one weighted graph satisfying both hypotheses
with p_3 < 0.

TAKEN FROM THE PAPER (transcribed inputs and targets, never used as evidence
for themselves):
  * the three DEFINITIONS that constitute the criterion, which no computation
    here can check against the source: "normalized" means 1^T K 1 = 1 with
    K = Omega^-1; the resistance capacity is tau_empty = 0, tau_{v} = 1/2 and
    tau_U = 1/2 + (1/2)(1^T Omega[U]^-1 1)^-1 for |U| >= 2; the resistance
    curvature is p_v = 1 - (1/2) sum_{e at v} c_e w_e.  The singleton value is
    load-bearing (it enters half the elementary slacks) and is the limit of the
    general formula as Omega[U] shrinks to the 1x1 zero matrix, so the
    convention is at least internally consistent;
  * the graph: the diamond K_4 - e on {0,1,2,3} with nonedge 01, in the edge
    order (02, 03, 12, 13, 23);
  * the exhibited weighting c = (109/300)*(1,2,1,2,2)  -- the ONLY object input;
  * the printed Laplacian and resistance matrices, the curvature vector
    (2/5, 2/5, 7/30, -1/30), the edge quantities (1/15)*(7,11,7,11,9), the
    normalizer 3924 and the 16-entry table of scaled submodularity slacks, the
    unit-conductance curvature (3/8,3/8,1/8,1/8), and the small-graph curvature
    vectors (1/2,1/2) and (1/2,0,1/2): all treated as TARGETS to be hit by an
    independent computation from c alone.

DERIVED HERE (recomputed from c in exact rational arithmetic; no floats):
  * the weighted Laplacian and, via L^+ = (L + J/n)^-1 - J/n, the full
    effective-resistance matrix Omega.  Omega is then recomputed on all six
    entries by two further routes that share no linear-algebra code with the
    first: ratios of Laplacian minors (det by Gaussian elimination, not
    inv), and the all-minors matrix-tree ratio omega_uv = F_uv / T obtained
    by enumerating weighted spanning trees and weighted two-component
    spanning forests (no matrix inversion and no determinant at all).  A
    fourth, deliberately partial cross-check recovers the five EDGE
    quantities c_e omega_e as spanning-tree inclusion probabilities; it
    covers 5 of the 6 entries and says nothing about the nonedge entry
    omega_01 = 200/109, which is why the forest route above exists;
  * K = Omega^-1 and 1^T K 1, deciding the normalization hypothesis;
  * the resistance capacity tau on all 16 subsets, all 24 elementary slacks and
    all 256 subset pairs, deciding the submodularity hypothesis;
  * the curvature vector by its definition p_v = 1 - (1/2) sum_{e at v} c_e w_e,
    and the sign of its minimum entry -- the load-bearing refutation;
  * the normalizing factor of (1,2,1,2,2), scale invariance of p, and the
    diamond's unit-conductance curvature, which bound the scope of the result;
  * a census of connected simple graphs on 2 and 3 vertices over an exact
    rational conductance grid, for the vertex-minimality claim.

Some of the checks are SOLVER SELF-TESTS rather than tests of the paper: they
assert theorems that hold for every positively weighted connected graph, or
compare two internal recomputations of the same quantity, so they validate
the arithmetic and the solver but cannot fail however the exhibited object is
perturbed.  They are listed in SELF_TEST_CHECKS below and named, with their
count derived from that list, in the program's closing NOT RE-RUN line:
  * resistance_matrix_is_a_metric (effective resistance is always a metric),
  * resistance_matrix_via_kirchhoff_cofactors (matrix-tree and all-minors
    matrix-tree theorems, read through determinants),
  * resistance_matrix_via_spanning_forest_ratios (the same two theorems, read
    through subset enumeration),
  * edge_probabilities_sum_to_n_minus_1 (Foster's identity, plus the fact
    that no edge of the diamond is a bridge),
  * edge_probabilities_from_spanning_tree_enumeration (matrix-tree theorem),
  * capacity_defined_on_all_16_subsets,
  * capacity_on_pairs_matches_closed_form (2x2 algebra on a zero-diagonal
    block),
  * elementary_slack_census_is_24 (a count fixed by n = 4).
They are kept because they would catch a corrupted
solver, which is a real failure mode here, but they are not evidence about
the object.  The evidence for the paper is carried by the other checks, each
of which was confirmed to report FAIL when the object, a paper target or the
solver is corrupted.

Standard library only.  Exit status is 0 if and only if every check passes.
"""

from fractions import Fraction as F
from itertools import combinations, product

CHECKS = []

# Checks that assert a theorem true of every positively weighted connected
# graph, or that compare two internal recomputations: they police the solver,
# not the exhibited object, and no perturbation of the object can make them
# fail.  Named in the closing NOT RE-RUN line, with the count derived here.
SELF_TEST_CHECKS = [
    "resistance_matrix_is_a_metric",
    "resistance_matrix_via_kirchhoff_cofactors",
    "resistance_matrix_via_spanning_forest_ratios",
    "edge_probabilities_sum_to_n_minus_1",
    "edge_probabilities_from_spanning_tree_enumeration",
    "capacity_defined_on_all_16_subsets",
    "capacity_on_pairs_matches_closed_form",
    "elementary_slack_census_is_24",
]

# ------------------------------------------------------------------
# VALUES TAKEN FROM THE PAPER (transcribed inputs / targets)
# ------------------------------------------------------------------
N = 4                                     # vertex set {0,1,2,3}
EDGES = [(0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]   # diamond, nonedge 01
NONEDGE = (0, 1)
SCALE = F(109, 300)
COND = [SCALE * F(k) for k in (1, 2, 1, 2, 2)]     # the exhibited weighting

PAPER_L_INT = [[3, 0, -1, -2],
               [0, 3, -1, -2],
               [-1, -1, 4, -2],
               [-2, -2, -2, 6]]
PAPER_OMEGA_NUM = [[0, 200, 140, 110],
                   [200, 0, 140, 110],
                   [140, 140, 0, 90],
                   [110, 110, 90, 0]]
PAPER_P = [F(2, 5), F(2, 5), F(7, 30), F(-1, 30)]
PAPER_CW = [F(7, 15), F(11, 15), F(7, 15), F(11, 15), F(9, 15)]
GAMMA_DEN = 3924
# rows keyed by (i,j); columns are S = {}, {a}, {b}, {a,b} with a<b the
# two remaining vertices
PAPER_GAMMA = {
    (0, 1): [162, 560, 165, 0],
    (0, 2): [702, 1100, 414, 249],
    (0, 3): [972, 975, 684, 124],
    (2, 3): [1152, 864, 864, 13],
}
PAPER_UNIT_P = [F(3, 8), F(3, 8), F(1, 8), F(1, 8)]
PAPER_K2_P = [F(1, 2), F(1, 2)]
PAPER_P3_P = [F(1, 2), F(0), F(1, 2)]
SWAP = {0: 1, 1: 0, 2: 2, 3: 3}           # the automorphism 0 <-> 1


# ------------------------------------------------------------------
# exact linear algebra over Fraction
# ------------------------------------------------------------------
def inv(M):
    """Exact inverse by Gauss-Jordan; returns None if singular."""
    n = len(M)
    A = [[F(M[i][j]) for j in range(n)] + [F(1) if i == j else F(0)
                                           for j in range(n)]
         for i in range(n)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        d = A[col][col]
        A[col] = [x / d for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return [row[n:] for row in A]


def matvec(M, v):
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def is_sym(M):
    n = len(M)
    return all(M[i][j] == M[j][i] for i in range(n) for j in range(n))


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)


# ------------------------------------------------------------------
# graph quantities, all exact
# ------------------------------------------------------------------
def laplacian(n, edges, cond):
    L = [[F(0)] * n for _ in range(n)]
    for (u, v), c in zip(edges, cond):
        L[u][u] += c
        L[v][v] += c
        L[u][v] -= c
        L[v][u] -= c
    return L


def resistance_matrix(n, edges, cond):
    """Omega_uv = (e_u-e_v)^T L^+ (e_u-e_v) via L^+ = (L+J/n)^-1 - J/n."""
    L = laplacian(n, edges, cond)
    M = [[L[i][j] + F(1, n) for j in range(n)] for i in range(n)]
    Mi = inv(M)
    if Mi is None:
        return None
    Lp = [[Mi[i][j] - F(1, n) for j in range(n)] for i in range(n)]
    return [[Lp[u][u] + Lp[v][v] - 2 * Lp[u][v] for v in range(n)]
            for u in range(n)], Lp


def connected(n, edges):
    seen, stack = {0}, [0]
    adj = {v: [] for v in range(n)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    while stack:
        for w in adj[stack.pop()]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def tau(Omega, U):
    """Resistance capacity of U (paper's definition, normalized weighting)."""
    U = sorted(U)
    if not U:
        return F(0)
    if len(U) == 1:
        return F(1, 2)
    sub = [[Omega[a][b] for b in U] for a in U]
    Si = inv(sub)
    if Si is None:
        return None
    s = sum(Si[i][j] for i in range(len(U)) for j in range(len(U)))
    if s == 0:
        return None
    return F(1, 2) + F(1, 2) / s


def curvature(n, edges, cond, Omega):
    p = [F(1)] * n
    for (u, v), c in zip(edges, cond):
        p[u] -= F(1, 2) * c * Omega[u][v]
        p[v] -= F(1, 2) * c * Omega[u][v]
    return p


# ------------------------------------------------------------------
# checks: the exhibited object
# ------------------------------------------------------------------
def check_graph():
    pairs = set(tuple(sorted(e)) for e in EDGES)
    allp = set(combinations(range(N), 2))
    deg = sorted(sum(1 for e in EDGES if v in e) for v in range(N))
    ok = (len(EDGES) == len(pairs) == len(COND) == 5 and N == 4
          and all(u != v for u, v in EDGES)
          and allp - pairs == {NONEDGE} and deg == [2, 2, 3, 3]
          and connected(N, EDGES))
    ck("graph_is_diamond_K4_minus_e", ok,
       "n=%d m=%d degrees=%s nonedge=%s connected=%s"
       % (N, len(pairs), deg, NONEDGE, connected(N, EDGES)))


def check_weighting():
    pos = all(c > 0 for c in COND)
    ratios = [c / SCALE for c in COND]
    d = dict(zip((tuple(sorted(e)) for e in EDGES), COND))
    autom = all(d[tuple(sorted((SWAP[u], SWAP[v])))] == c
                for (u, v), c in zip(EDGES, COND))
    ck("weighting_positive_and_0_1_symmetric", pos and autom,
       "c=%s*%s swap-invariant=%s" % (SCALE, [str(r) for r in ratios], autom))


def check_laplacian(L):
    ok = (is_sym(L) and all(sum(row) == 0 for row in L)
          and all(L[i][j] == SCALE * PAPER_L_INT[i][j]
                  for i in range(N) for j in range(N)))
    ck("laplacian_matches_paper", ok, "L = (109/300)*integer matrix, rows sum 0")


def check_omega(Omega):
    ok = (is_sym(Omega)
          and all(Omega[i][i] == 0 for i in range(N))
          and all(Omega[i][j] > 0 for i in range(N) for j in range(N)
                  if i != j)
          and all(109 * Omega[i][j] == PAPER_OMEGA_NUM[i][j]
                  for i in range(N) for j in range(N)))
    tri = all(Omega[i][j] <= Omega[i][k] + Omega[k][j]
              for i in range(N) for j in range(N) for k in range(N))
    ck("resistance_matrix_matches_paper", ok, "109*Omega = paper integer matrix")
    ck("resistance_matrix_is_a_metric", tri,
       "zero diagonal, positive off-diagonal, triangle inequality on %d triples"
       % (N ** 3))


def check_normalized(Omega):
    K = inv(Omega)
    if K is None:
        ck("weighting_is_normalized", False, "Omega singular")
        return None
    K1 = matvec(K, [F(1)] * N)
    tot = sum(K1)
    ck("weighting_is_normalized", tot == 1, "1^T K 1 = %s (hypothesis of 7.6)"
       % tot)
    ck("K1_matches_paper_vector", K1 == PAPER_P,
       "K1 = (" + ", ".join(str(x) for x in K1) + ")")
    return K1


def check_edge_probabilities(Omega):
    cw = [c * Omega[u][v] for (u, v), c in zip(EDGES, COND)]
    ok = all(a == b for a, b in zip(cw, PAPER_CW))
    ck("edge_inclusion_probabilities_match_paper", ok,
       "(c_e w_e) = (" + ", ".join(str(x) for x in cw) + ")")
    ck("edge_probabilities_sum_to_n_minus_1",
       sum(cw) == N - 1 and all(0 < x < 1 for x in cw),
       "sum = %s = n-1, each strictly in (0,1)" % sum(cw))


def det(M):
    M = [row[:] for row in M]
    n = len(M)
    d = F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return d


def check_matrix_tree(Omega):
    """c_e w_e is the inclusion probability of e in the weighted random
    spanning tree: recomputed by enumerating spanning trees."""
    m = len(EDGES)
    total, per = F(0), [F(0)] * m
    trees = 0
    for T in combinations(range(m), N - 1):
        if connected(N, [EDGES[i] for i in T]):
            w = F(1)
            for i in T:
                w *= COND[i]
            total += w
            trees += 1
            for i in T:
                per[i] += w
    probs = [x / total for x in per] if total else None
    cw = [c * Omega[u][v] for (u, v), c in zip(EDGES, COND)]
    ck("edge_probabilities_from_spanning_tree_enumeration",
       probs is not None and probs == cw,
       "%d spanning trees enumerated; weighted inclusion probabilities equal "
       "(c_e w_e), confirming the interpretation used for minimality" % trees)


def check_kirchhoff(Omega, L):
    """Second, independent route to Omega: ratios of Laplacian minors."""
    def minor(drop):
        keep = [i for i in range(N) if i not in drop]
        return det([[L[i][j] for j in keep] for i in keep])

    m = len(EDGES)
    total = F(0)
    for T in combinations(range(m), N - 1):
        if connected(N, [EDGES[i] for i in T]):
            w = F(1)
            for i in T:
                w *= COND[i]
            total += w
    tree_ok = all(minor({v}) == total for v in range(N))
    cof_ok = all(Omega[u][v] * minor({u}) == minor({u, v})
                 for u, v in combinations(range(N), 2))
    ck("resistance_matrix_via_kirchhoff_cofactors", tree_ok and cof_ok,
       "spanning-tree weight %s equals every first cofactor, and "
       "w_uv = det(L minus u,v)/det(L minus u)" % total)


def vertex_partition(n, edges):
    """Vertex classes of the subgraph spanned by an edge subset."""
    root = list(range(n))

    def find(x):
        while root[x] != x:
            root[x] = root[root[x]]
            x = root[x]
        return x

    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            root[ru] = rv
    classes = {}
    for v in range(n):
        classes.setdefault(find(v), set()).add(v)
    return list(classes.values())


def check_forest_ratios(Omega):
    """Third route to Omega, covering ALL SIX entries, the nonedge included.

    All-minors matrix-tree theorem: omega_uv = F_uv / T, where T is the total
    weight of the spanning trees and F_uv the total weight of the spanning
    forests with exactly two components, one containing u and the other v.
    Weights are products of conductances.  This route uses neither inv() nor
    det(): it is pure subset enumeration, so it is independent of the
    pseudoinverse route and of the cofactor route.  It exists because
    check_matrix_tree covers only the five edge entries and therefore never
    touches omega_01, the largest entry and the one driving the 01 row of the
    Gamma table.
    """
    m = len(EDGES)
    total = F(0)
    for T in combinations(range(m), N - 1):
        if connected(N, [EDGES[i] for i in T]):
            w = F(1)
            for i in T:
                w *= COND[i]
            total += w
    sep = {pair: F(0) for pair in combinations(range(N), 2)}
    forests = 0
    for Fs in combinations(range(m), N - 2):
        es = [EDGES[i] for i in Fs]
        part = vertex_partition(N, es)
        if len(part) != N - len(es):        # acyclic iff n - |F| classes
            continue
        forests += 1
        w = F(1)
        for i in Fs:
            w *= COND[i]
        for u, v in sep:
            if not any(u in cl and v in cl for cl in part):
                sep[(u, v)] += w
    got = {pair: (sep[pair] / total if total else None) for pair in sep}
    ok = (total > 0 and forests > 0
          and all(got[(u, v)] == Omega[u][v]
                  for u, v in combinations(range(N), 2)))
    ck("resistance_matrix_via_spanning_forest_ratios", ok,
       "%d two-component spanning forests over tree weight %s reproduce all "
       "%d entries of Omega including the nonedge, omega_01 = %s"
       % (forests, total, len(sep), got[(0, 1)]))


def check_curvature(Omega, K1):
    p = curvature(N, EDGES, COND, Omega)
    ck("curvature_matches_paper_vector", p == PAPER_P,
       "p = (" + ", ".join(str(x) for x in p) + ")")
    ck("curvature_equals_K1_independent_route", K1 is not None and p == K1,
       "definition p_v = 1 - (1/2) sum c_e w_e agrees with Omega^{-1} 1")
    return p


def check_conclusion_violated(p):
    mn = min(p)
    arg = [v for v in range(N) if p[v] == mn]
    ck("CONCLUSION_VIOLATED_curvature_negative", mn < 0,
       "min_v p_v = %s at vertex %s; conjecture demands p >= 0" % (mn, arg))
    ck("negative_coordinate_equals_minus_1_over_30", mn == F(-1, 30),
       "p_3 = %s" % mn)


# ------------------------------------------------------------------
# checks: the resistance capacity and its submodularity
# ------------------------------------------------------------------
def all_subsets(n):
    out = []
    for r in range(n + 1):
        out.extend(combinations(range(n), r))
    return out


def build_tau(Omega):
    tv = {}
    for U in all_subsets(N):
        tv[U] = tau(Omega, U)
    return tv


def check_tau_wellformed(tv, Om):
    defined = all(v is not None for v in tv.values())
    ck("capacity_defined_on_all_16_subsets", defined and len(tv) == 16,
       "every Omega[U] with |U|>=2 invertible and 1^T Omega[U]^-1 1 nonzero")
    if not defined:
        return
    pair = all(tv[(u, v)] == F(1, 2) + Om[u][v] / 4
               for u, v in combinations(range(N), 2))
    ck("capacity_on_pairs_matches_closed_form", pair,
       "tau_{u,v} = 1/2 + w_uv/4 for all 6 pairs, independent of the "
       "2x2 inversion routine")
    mono = all(tv[U] <= tv[tuple(sorted(set(U) | {v}))]
               for U in tv for v in range(N) if v not in U)
    ck("capacity_monotone_and_tau_V_equals_one",
       mono and tv[tuple(range(N))] == 1,
       "monotone on the subset lattice; tau_V = %s, which equals 1 exactly "
       "when the weighting is normalized" % tv[tuple(range(N))])
    inv_ok = all(tv[tuple(sorted(SWAP[v] for v in U))] == tv[U] for U in tv)
    ck("capacity_invariant_under_0_1_swap", inv_ok,
       "justifies deriving rows 12,13 from rows 02,03")


def slacks(tv):
    """All elementary slacks Delta_ij(S), keyed by (i,j,S)."""
    out = {}
    for i, j in combinations(range(N), 2):
        rest = sorted(set(range(N)) - {i, j})
        for r in range(len(rest) + 1):
            for S in combinations(rest, r):
                s = set(S)
                d = (tv[tuple(sorted(s | {i}))] + tv[tuple(sorted(s | {j}))]
                     - tv[tuple(sorted(s))] - tv[tuple(sorted(s | {i, j}))])
                out[(i, j, S)] = d
    return out


def check_slack_census(sl):
    ck("elementary_slack_census_is_24", len(sl) == 24,
       "C(4,2) pairs * 2^2 subsets = %d slacks enumerated" % len(sl))
    ints = all((d * GAMMA_DEN).denominator == 1 for d in sl.values())
    ck("gamma_values_are_integers", ints,
       "3924 = 36*109 clears every denominator")


def check_gamma_table(sl):
    bad = []
    for (i, j), row in PAPER_GAMMA.items():
        rest = sorted(set(range(N)) - {i, j})
        a, b = rest
        cols = [(), (a,), (b,), (a, b)]
        for want, S in zip(row, cols):
            got = sl[(i, j, S)] * GAMMA_DEN
            if got != want:
                bad.append("%d%d,S=%s: got %s want %d" % (i, j, S, got, want))
    ck("gamma_table_matches_paper", not bad,
       "16 tabulated values reproduced" if not bad else "; ".join(bad[:4]))
    # the two rows the paper obtains by applying the automorphism
    img = []
    for (i, j) in [(0, 2), (0, 3)]:
        si, sj = sorted((SWAP[i], SWAP[j]))
        rest = sorted(set(range(N)) - {i, j})
        rs = sorted(set(range(N)) - {si, sj})
        for S in [(), (rest[0],), (rest[1],), tuple(rest)]:
            T = tuple(sorted(SWAP[v] for v in S))
            img.append(sl[(i, j, S)] == sl[(si, sj, T)])
    ck("automorphism_derived_rows_agree", all(img) and len(img) == 8,
       "rows 12,13 equal rows 02,03 under 0<->1; table covers all 24 slacks")


def check_submodular(tv, sl):
    neg = {k: v for k, v in sl.items() if v < 0}
    ck("HYPOTHESIS_capacity_submodular_elementary", not neg,
       "all 24 elementary slacks >= 0, min = %s" % min(sl.values()))
    worst, pairs = None, 0
    for A in all_subsets(N):
        for B in all_subsets(N):
            sa, sb = set(A), set(B)
            d = (tv[A] + tv[B] - tv[tuple(sorted(sa | sb))]
                 - tv[tuple(sorted(sa & sb))])
            pairs += 1
            if worst is None or d < worst:
                worst = d
    ck("HYPOTHESIS_capacity_submodular_all_pairs", worst >= 0,
       "%d ordered subset pairs, min slack = %s" % (pairs, worst))


NONEDGE_SLACK = (0, 1, (2, 3))     # Delta_{01}(V minus {0,1})
PAPER_MARGIN_GAMMA = 13            # the smallest strictly positive table entry


def normalized_tau(raw):
    """tau of the weighting proportional to raw, rescaled to be normalized."""
    Om, _ = resistance_matrix(N, EDGES, raw)
    Ki = inv(Om)
    if Ki is None:
        return None
    t = 1 / sum(sum(r) for r in Ki)
    return build_tau([[x / t for x in row] for row in Om])


def check_submodularity_margin(sl):
    """The minimum slack reported as 0 is FORCED, not a knife edge.

    Delta_{ij}(V minus {i,j}) vanishes identically whenever i,j are
    non-adjacent, so the single zero in the paper's table is a property of the
    nonedge 01 and not of the exhibited weighting.  Two consequences: the
    hypothesis cannot be dismissed as boundary contact, and no non-complete
    graph has a STRICTLY submodular capacity, so the strict reading of the
    conjecture would be vacuous here.  The margin that actually carries the
    hypothesis is the minimum over the other 23 slacks.
    """
    zero = set(k for k, v in sl.items() if v == 0)
    probes = []
    for raw in ([F(1)] * 5, [F(1), F(3), F(2), F(5), F(7)],
                [F(5), F(1), F(2), F(8), F(3)]):
        tvp = normalized_tau(raw)
        probes.append(tvp is not None and slacks(tvp)[NONEDGE_SLACK] == 0)
    ck("nonedge_slack_vanishes_for_structural_reasons",
       zero == {NONEDGE_SLACK} and all(probes) and len(probes) == 3,
       "the one vanishing slack is Delta_01(V-{0,1}); it is still 0 at three "
       "unrelated weightings of the same graph, so the zero comes from the "
       "nonedge, not from c")
    rest = [v for k, v in sl.items() if k != NONEDGE_SLACK]
    mn = min(rest)
    ck("submodularity_margin_strictly_positive",
       len(rest) == 23 and mn > 0 and mn * GAMMA_DEN == PAPER_MARGIN_GAMMA,
       "min over the 23 slacks that are not identically zero is %s, i.e. "
       "Gamma = %s attained at Delta_23({0,1}); submodularity is therefore "
       "strict wherever it can be" % (mn, mn * GAMMA_DEN))


# ------------------------------------------------------------------
# checks: the paper's scope remarks and its minimality proposition
# ------------------------------------------------------------------
def check_normalizing_scale():
    raw = [F(k) for k in (1, 2, 1, 2, 2)]
    Om, _ = resistance_matrix(N, EDGES, raw)
    Kr = inv(Om)
    tot = sum(sum(r) for r in Kr)
    ck("scale_109_over_300_is_the_normalizer", 1 / tot == SCALE,
       "1^T K 1 = %s for c=(1,2,1,2,2), so the normalizing factor is %s"
       % (tot, 1 / tot))
    p_raw = curvature(N, EDGES, raw, Om)
    p_sc = curvature(N, EDGES, COND,
                     resistance_matrix(N, EDGES, COND)[0])
    ck("curvature_invariant_under_rescaling", p_raw == p_sc,
       "p unchanged by the common factor, as the paper asserts")


def check_diamond_resistance_positive():
    unit = [F(1)] * len(EDGES)
    Om, _ = resistance_matrix(N, EDGES, unit)
    p = curvature(N, EDGES, unit, Om)
    ck("unit_conductance_curvature_matches_paper", p == PAPER_UNIT_P,
       "p = (" + ", ".join(str(x) for x in p) + ") >= 0, so the diamond is "
       "resistance positive and 7.7 is untouched")
    t = 1 / sum(sum(r) for r in inv(Om))
    Om2, _ = resistance_matrix(N, EDGES, [t * c for c in unit])
    p2 = curvature(N, EDGES, [t * c for c in unit], Om2)
    norm = sum(sum(r) for r in inv(Om2))
    ck("rescaled_unit_weighting_normalized_and_nonnegative",
       norm == 1 and p2 == PAPER_UNIT_P and min(p2) >= 0,
       "factor %s gives 1^T K 1 = %s with the same nonnegative p" % (t, norm))


def check_small_graphs():
    """Exhaustive over connected simple graphs on 2 and 3 vertices, with a
    grid of positive rational conductances (exact)."""
    vals = [F(1, 3), F(1, 2), F(1), F(2), F(3), F(5, 2), F(7, 3)]
    fams = [("K2", 2, [(0, 1)]),
            ("P3", 3, [(0, 1), (1, 2)]),
            ("K3", 3, [(0, 1), (1, 2), (0, 2)])]
    tested, worst, detail = 0, None, {}
    tree_ok, k3_ok = True, True
    for name, n, es in fams:
        m = len(es)
        wmin = None
        for combo in product(vals, repeat=m):
            Om, _ = resistance_matrix(n, es, list(combo))
            if Om is None:
                continue
            p = curvature(n, es, list(combo), Om)
            tested += 1
            lo = min(p)
            if wmin is None or lo < wmin:
                wmin = lo
            if worst is None or lo < worst:
                worst = lo
            if name == "K2" and p != PAPER_K2_P:
                tree_ok = False
            if name == "P3" and p != PAPER_P3_P:
                tree_ok = False
            if name == "K3":
                cw = {tuple(sorted(e)): c * Om[e[0]][e[1]]
                      for e, c in zip(es, list(combo))}
                for v in range(3):
                    opp = tuple(sorted(set(range(3)) - {v}))
                    if p[v] != cw[opp] / 2 or p[v] <= 0:
                        k3_ok = False
        detail[name] = wmin
    ck("small_graph_census_no_negative_curvature", worst >= 0,
       "%d weightings over K2,P3,K3; min coordinate = %s (minima %s)"
       % (tested, worst, {k: str(v) for k, v in detail.items()}))
    ck("tree_curvature_vectors_are_weight_independent", tree_ok,
       "K2 gives (1/2,1/2) and P3 gives (1/2,0,1/2) at every grid weighting")
    ck("K3_curvature_is_half_the_opposite_edge_probability", k3_ok,
       "p_v = c_e w_e / 2 > 0 for the edge e opposite v, at every weighting")


def finish():
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


def main():
    print("object: diamond K_4 - e on {0,1,2,3}, nonedge %s" % (NONEDGE,))
    print("edge order %s, conductances c = (109/300)*(1,2,1,2,2)"
          % (EDGES,))
    check_graph()
    check_weighting()
    L = laplacian(N, EDGES, COND)
    check_laplacian(L)
    res = resistance_matrix(N, EDGES, COND)
    if res is None:
        ck("resistance_matrix_matches_paper", False, "Laplacian solve failed")
        return finish()
    Omega, _ = res
    check_omega(Omega)
    check_kirchhoff(Omega, L)
    K1 = check_normalized(Omega)
    check_edge_probabilities(Omega)
    check_matrix_tree(Omega)
    check_forest_ratios(Omega)
    p = check_curvature(Omega, K1)
    tv = build_tau(Omega)
    check_tau_wellformed(tv, Omega)
    if all(v is not None for v in tv.values()):
        sl = slacks(tv)
        check_slack_census(sl)
        check_gamma_table(sl)
        check_submodular(tv, sl)
        check_submodularity_margin(sl)
    check_conclusion_violated(p)
    check_normalizing_scale()
    check_diamond_resistance_positive()
    check_small_graphs()
    ran = [name for name, _ in CHECKS]
    drifted = [name for name in SELF_TEST_CHECKS if name not in ran]
    print("NOT RE-RUN: the vertex-minimality Proposition quantifies over ALL "
          "positive conductances on the order-<=3 graphs, a continuum. The "
          "graph enumeration here is exhaustive (K2, P3, K3 are the only "
          "connected simple graphs on 2 or 3 vertices) but the conductances "
          "are a finite exact rational grid, so this census CORROBORATES the "
          "Proposition and does not establish it; the two structural "
          "identities checked above, which hold at every grid point, are "
          "checked here rather than proved. The paper does not rest on this "
          "census: its proof of the Proposition is a complete argument "
          "(c_e omega_e is the spanning-tree inclusion probability of e, so "
          "it is 1 on a bridge and less than 1 on a nonbridge; K2 and P3 are "
          "trees; every vertex of K3 meets two nonbridges), and the grid is "
          "corroboration of that argument, not a substitute for it. Also not "
          "re-run: that the definitions of normalization, resistance "
          "capacity and resistance curvature used here, and the numbering of "
          "the statement refuted and of the companion statement left "
          "standing, agree with the cited source; those are transcribed. And "
          "a caveat on the count: %d of the %d checks are SOLVER SELF-TESTS "
          "rather than tests of the paper -- %s -- each asserting a theorem "
          "true of every positively weighted connected graph, or comparing "
          "two internal recomputations, or counting something fixed by n=4. "
          "They would catch a corrupted solver but cannot fail however the "
          "exhibited object is perturbed, so the refutation is carried by "
          "the remaining %d checks.%s Everything else in the paper is "
          "reproduced exactly, in exact rational arithmetic."
          % (len(SELF_TEST_CHECKS), len(CHECKS), ", ".join(SELF_TEST_CHECKS),
             len(CHECKS) - len(SELF_TEST_CHECKS),
             "" if not drifted else
             " WARNING: %s did not run under those names, so this list has "
             "drifted from the program." % ", ".join(drifted)))
    return finish()


if __name__ == "__main__":
    try:
        CODE = main()
    except Exception as exc:                  # never exit without a verdict
        ck("program_ran_to_completion", False, "unexpected error: %r" % (exc,))
        CODE = finish()
    raise SystemExit(CODE)

#!/usr/bin/env python3
"""Independent verification of the certified RM(7,3,2) mixed radial Moore graph.

TAKEN FROM THE PAPER (inputs, transcribed verbatim, nothing else assumed)
------------------------------------------------------------------------
  * the 104-row adjacency certificate: for each vertex v in {0,...,103},
    the undirected neighbour list U(v) and the out-arc head list A+(v).
    This is the only combinatorial input; it is transcribed below in
    CERTIFICATE exactly as printed, one row per line, "v: U(v) | A+(v)".
  * the target parameters (r, z, k) = (7, 3, 2).
  * the published status vector of the arc-swapped order-108 graphs,
    204^98 205^8 209^2, quoted from the paper's Section 2 discussion of
    the transposed table cell.
  * the numbers the paper asserts, used ONLY as right-hand sides to
    compare recomputed values against:
      |V| = 104, |E| = 364, |A| = 312, M(7,3,2) = 104,
      eccentricity multiset {2^1, 3^103}, radius 2, diameter 3,
      unique central vertex, central layer profile [1, 10, 93],
      Moore status sigma_M = 196, total status 23762, N_1 = 3378,
      and, for the order-108 graphs, baseline 204 and N_1 = 18.

DERIVED HERE (recomputed from the certificate alone)
----------------------------------------------------
  * well-formedness of the certificate: index set, row sizes, no loops,
    no repeats, undirected symmetry, edge and arc cardinalities.
  * every hypothesis of "totally (7,3)-regular simple mixed graph":
    undirected degree 7 everywhere, arc out-degree and in-degree 3
    everywhere, no directed digons, no edge/arc conflicts.
  * the diameter-two mixed Moore bound M(r,z,2) from its defining count
    1 + (r+z) + r(r-1+z) + z(r+z), the two-sided estimate in terms of
    s = r + z, and the exhaustive determination of which (r,z) realise
    the orders 104 and 108.
  * all-pairs distances by breadth-first search, cross-checked against
    an independent iterated-relation-composition computation; hence
    eccentricities, radius, diameter, the centre, the layer profile of
    the central vertex, every vertex status, the total status and N_1.
  * the load-bearing conclusion: order = M(7,3,2), radius = k = 2,
    diameter = k + 1 = 3, so an RM(7,3,2) graph exists, with N_1 = 3378.
  * the arithmetic of the parameter correction: the order-108 baseline
    204 belongs to (r,z) = (3,7) and not to (7,3), whose baseline is
    196; and the quoted status vector has 108 entries, not 104.
  * the identity 1 + (r+z) + r(r-1+z) + z(r+z) = (r+z)^2 + z + 1 and the
    step z = order - s^2 - 1 that the correction is derived from, plus
    the fact that the order and the baseline alone pin the parameter
    pair.  This part of the argument uses no table layout whatever; what
    a printed table's axes are indexed by is a fact about published
    pages and is disclosed as unchecked in the closing lines.

Standard library only, exact integer arithmetic throughout, no input
files, no network. Exit status 0 if and only if every check passes.
"""

import sys
from collections import deque

_CHECKS = []

# ---------------------------------------------------------------------------
# Values quoted from the paper.  Used only as right-hand sides of comparisons.
# ---------------------------------------------------------------------------
R_PAPER = 7                 # undirected degree
Z_PAPER = 3                 # directed in/out degree
K_PAPER = 2                 # claimed radius (diameter is claimed to be K+1)
ORDER_PAPER = 104
EDGES_PAPER = 364
ARCS_PAPER = 312
MOORE_BOUND_PAPER = 104
CENTRAL_PROFILE_PAPER = [1, 10, 93]
ECC_MULTISET_PAPER = {2: 1, 3: 103}
SIGMA_MOORE_PAPER = 196
TOTAL_STATUS_PAPER = 23762
N1_PAPER = 3378

# The order-108 arc-swapped status vector as printed in the literature,
# together with the baseline and norm the paper attributes to it.
STATUS_108_PAPER = [(204, 98), (205, 8), (209, 2)]
BASELINE_108_PAPER = 204
N1_108_PAPER = 18
ORDER_108 = 108

# ---------------------------------------------------------------------------
# The adjacency certificate, transcribed from the printed table.
# One line per vertex:  v: u1,...,u7 | a1,a2,a3
# ---------------------------------------------------------------------------
CERTIFICATE = """
0: 1,2,3,4,5,6,7 | 8,9,10
1: 0,11,12,13,14,15,16 | 17,18,19
2: 0,20,21,22,23,24,25 | 26,27,28
3: 0,29,30,31,32,33,34 | 35,36,37
4: 0,38,39,40,41,42,43 | 44,45,46
5: 0,47,48,49,50,51,52 | 53,54,55
6: 0,56,57,58,59,60,61 | 62,63,64
7: 0,65,66,67,68,69,70 | 71,72,73
8: 74,75,76,77,78,79,80 | 81,82,83
9: 84,85,86,87,88,89,90 | 91,92,93
10: 94,95,96,97,98,99,100 | 101,102,103
11: 1,26,32,40,56,86,93 | 2,15,87
12: 1,16,18,62,83,90,100 | 36,74,103
13: 1,16,23,64,69,73,79 | 5,92,93
14: 1,17,28,46,67,80,98 | 7,57,66
15: 1,22,39,46,48,52,59 | 6,30,61
16: 1,12,13,28,33,57,95 | 23,40,61
17: 14,32,71,76,80,95,96 | 28,34,53
18: 12,27,33,38,71,73,99 | 57,64,70
19: 26,27,36,77,81,87,91 | 11,25,34
20: 2,44,49,53,69,74,75 | 21,24,99
21: 2,35,60,68,79,88,102 | 29,31,94
22: 2,15,23,37,61,73,84 | 13,86,97
23: 2,13,22,29,42,62,89 | 41,43,78
24: 2,33,34,41,48,76,86 | 54,80,101
25: 2,37,45,60,62,67,76 | 4,35,63
26: 11,19,30,36,67,78,85 | 4,17,24
27: 18,19,45,46,55,57,67 | 24,75,90
28: 14,16,48,62,63,71,82 | 50,60,72
29: 3,23,59,63,79,90,98 | 67,84,95
30: 3,26,32,44,47,84,90 | 16,60,100
31: 3,45,50,59,60,87,96 | 13,33,42
32: 3,11,17,30,45,63,86 | 12,23,57
33: 3,16,18,24,34,38,52 | 25,26,79
34: 3,24,33,60,70,74,82 | 1,26,39
35: 21,39,42,44,71,79,80 | 77,78,83
36: 19,26,54,62,81,91,103 | 7,67,102
37: 22,25,60,70,93,101,103 | 49,89,100
38: 4,18,33,48,55,57,77 | 58,67,73
39: 4,15,35,43,61,64,103 | 3,62,94
40: 4,11,49,50,67,78,86 | 1,47,81
41: 4,24,52,83,85,91,95 | 8,16,102
42: 4,23,35,58,70,76,91 | 33,51,98
43: 4,39,53,70,93,99,100 | 48,54,101
44: 20,30,35,71,78,97,102 | 40,65,68
45: 25,27,31,32,73,82,87 | 16,20,89
46: 14,15,27,57,65,95,102 | 62,69,94
47: 5,30,58,61,77,91,100 | 41,42,52
48: 5,15,24,28,38,82,90 | 69,76,84
49: 5,20,40,54,82,85,99 | 3,31,64
50: 5,31,40,86,92,97,100 | 22,61,89
51: 5,52,56,65,73,76,89 | 27,49,75
52: 5,15,33,41,51,72,91 | 35,43,90
53: 20,43,57,61,74,84,101 | 30,59,98
54: 36,49,55,64,65,66,97 | 25,31,85
55: 27,38,54,72,85,100,102 | 22,80,96
56: 6,11,51,61,66,68,81 | 7,46,88
57: 6,16,27,38,46,53,84 | 2,71,96
58: 6,42,47,59,70,88,94 | 15,71,86
59: 6,15,29,31,58,68,90 | 44,82,91
60: 6,21,25,31,34,37,75 | 14,49,51
61: 6,22,39,47,53,56,75 | 2,55,83
62: 12,23,25,28,36,63,85 | 11,48,74
63: 28,29,32,62,70,98,101 | 15,39,56
64: 13,39,54,72,82,95,101 | 77,79,88
65: 7,46,51,54,72,84,96 | 8,77,97
66: 7,54,56,68,74,75,96 | 6,29,37
67: 7,14,25,26,27,40,77 | 11,87,100
68: 7,21,56,59,66,86,103 | 47,52,65
69: 7,13,20,77,82,93,98 | 3,12,27
70: 7,34,37,42,43,58,63 | 41,92,93
71: 17,18,28,35,44,72,73 | 0,19,46
72: 52,55,64,65,71,88,89 | 10,70,75
73: 13,18,22,45,51,71,83 | 9,21,32
74: 8,20,34,53,66,94,99 | 32,84,96
75: 8,20,60,61,66,92,98 | 9,59,82
76: 8,17,24,25,42,51,90 | 14,69,85
77: 8,19,38,47,67,69,92 | 0,18,32
78: 8,26,40,44,80,92,96 | 10,43,79
79: 8,13,21,29,35,89,93 | 19,37,81
80: 8,14,17,35,78,83,93 | 34,51,99
81: 19,36,56,87,89,97,102 | 17,38,68
82: 28,34,45,48,49,64,69 | 18,30,56
83: 12,41,73,80,94,95,101 | 58,63,76
84: 9,22,30,53,57,65,97 | 21,80,85
85: 9,26,41,49,55,62,102 | 28,38,68
86: 9,11,24,32,40,50,68 | 29,98,99
87: 9,19,31,45,81,89,94 | 20,47,56
88: 9,21,58,72,92,94,103 | 6,38,76
89: 9,23,51,72,79,81,87 | 42,74,97
90: 9,12,29,30,48,59,76 | 20,50,70
91: 19,36,41,42,47,52,99 | 13,90,95
92: 50,75,77,78,88,98,101 | 0,45,58
93: 11,37,43,69,79,80,97 | 39,55,66
94: 10,58,74,83,87,88,103 | 4,23,91
95: 10,16,17,41,46,64,83 | 48,52,59
96: 10,17,31,65,66,78,99 | 5,22,73
97: 10,44,50,54,81,84,93 | 5,14,95
98: 10,14,29,63,69,75,92 | 12,65,66
99: 10,18,43,49,74,91,96 | 40,72,87
100: 10,12,43,47,50,55,103 | 44,53,60
101: 37,53,63,64,83,92,102 | 33,36,103
102: 21,44,46,55,81,85,101 | 45,50,88
103: 36,37,39,68,88,94,100 | 1,78,86
"""


def ck(name, ok, detail=""):
    _CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)


def summarize():
    n = len(_CHECKS)
    bad = [nm for nm, ok in _CHECKS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


def parse_certificate(text):
    """Decode the printed table into rows[v] = (list U(v), list A+(v)).

    Deliberately strict: any malformed line, missing vertex label, repeated
    vertex label or non-integer entry raises, so a corrupted certificate
    cannot be silently reinterpreted.
    """
    rows = {}
    order = []
    for raw in text.strip().split("\n"):
        line = raw.strip()
        if not line:
            continue
        head, rest = line.split(":", 1)
        v = int(head.strip())
        und, arc = rest.split("|", 1)
        u_list = [int(x) for x in und.strip().split(",")]
        a_list = [int(x) for x in arc.strip().split(",")]
        if v in rows:
            raise ValueError("vertex %d listed twice" % v)
        rows[v] = (u_list, a_list)
        order.append(v)
    return rows, order


def check_wellformed(rows, order):
    n = len(rows)
    ok_labels = order == list(range(n))
    ck("certificate_rows_are_vertices_0_to_n_minus_1",
       ok_labels and n == ORDER_PAPER,
       "n=%d, labels 0..%d in order: %s" % (n, n - 1, ok_labels))
    ck("order_is_104", n == ORDER_PAPER, "|V|=%d" % n)

    bad_range = [v for v in rows
                 for x in rows[v][0] + rows[v][1] if not 0 <= x < n]
    ck("all_entries_in_vertex_range", not bad_range,
       "offending rows: %s" % sorted(set(bad_range))[:5])

    loops = [v for v in rows if v in rows[v][0] or v in rows[v][1]]
    ck("no_loops", not loops, "rows with a self-reference: %s" % loops[:5])

    dup = [v for v in rows
           if len(set(rows[v][0])) != len(rows[v][0])
           or len(set(rows[v][1])) != len(rows[v][1])]
    ck("no_repeated_entries_within_a_row", not dup,
       "rows with repeats: %s" % dup[:5])

    sizes = sorted(set((len(rows[v][0]), len(rows[v][1])) for v in rows))
    ck("every_row_lists_7_neighbours_and_3_arc_heads",
       sizes == [(R_PAPER, Z_PAPER)],
       "(|U|,|A+|) values present: %s" % sizes)


def check_mixed_graph_axioms(rows):
    """Simple mixed graph: symmetric edge relation, no digon, no conflict."""
    n = len(rows)
    und = dict((v, set(rows[v][0])) for v in rows)
    arc = dict((v, set(rows[v][1])) for v in rows)

    asym = [(v, u) for v in rows for u in und[v] if v not in und[u]]
    ck("undirected_adjacency_is_symmetric", not asym,
       "unmatched pairs: %d, first: %s" % (len(asym), asym[:3]))

    edges = set()
    for v in rows:
        for u in und[v]:
            edges.add((min(u, v), max(u, v)))
    ck("edge_count_is_364",
       len(edges) == EDGES_PAPER and 2 * len(edges) == n * R_PAPER,
       "|E|=%d, 104*7/2=%d" % (len(edges), n * R_PAPER // 2))

    arc_total = sum(len(arc[v]) for v in rows)
    ck("arc_count_is_312",
       arc_total == ARCS_PAPER and arc_total == n * Z_PAPER,
       "|A|=%d, 104*3=%d" % (arc_total, n * Z_PAPER))

    digons = [(v, u) for v in rows for u in arc[v] if v in arc[u]]
    ck("no_directed_digons", not digons,
       "digons: %d, first: %s" % (len(digons), digons[:3]))

    conflicts = [(v, u) for v in rows for u in arc[v] if u in und[v]]
    ck("no_edge_arc_conflicts", not conflicts,
       "conflicting pairs: %d, first: %s" % (len(conflicts), conflicts[:3]))
    return und, arc


def check_total_regularity(rows, und, arc):
    """Hypothesis of the statement: totally (7,3)-regular."""
    bad_e = sorted(v for v in rows if len(und[v]) != R_PAPER)
    ck("undirected_degree_is_7_at_every_vertex", not bad_e,
       "offenders: %s" % bad_e[:5])

    bad_out = sorted(v for v in rows if len(arc[v]) != Z_PAPER)
    ck("arc_out_degree_is_3_at_every_vertex", not bad_out,
       "offenders: %s" % bad_out[:5])

    indeg = dict((v, 0) for v in rows)
    for v in rows:
        for u in arc[v]:
            indeg[u] += 1
    bad_in = sorted(v for v in rows if indeg[v] != Z_PAPER)
    ck("every_vertex_is_an_arc_head_exactly_3_times", not bad_in,
       "in-degree values present: %s" % sorted(set(indeg.values())))

    out_nb = dict((v, und[v] | arc[v]) for v in rows)
    sizes = sorted(set(len(out_nb[v]) for v in rows))
    ck("out_neighbourhood_size_is_r_plus_z_equals_10",
       sizes == [R_PAPER + Z_PAPER],
       "sizes present: %s" % sizes)
    return out_nb


def bfs_distances(out_nb, n):
    """dist[s][t] by breadth-first search respecting arc orientation."""
    dist = {}
    for s in range(n):
        d = {s: 0}
        queue = deque([s])
        while queue:
            x = queue.popleft()
            dx = d[x] + 1
            for y in out_nb[x]:
                if y not in d:
                    d[y] = dx
                    queue.append(y)
        dist[s] = d
    return dist


def closure_distances(out_nb, n, cap):
    """Independent recomputation: iterated composition of the out-relation.

    reach[i] is the set of pairs (s,t) with a directed walk of length <= i.
    Distances are read off as the first i at which (s,t) appears.  This uses
    no queue and no visited-set logic, so it is a genuinely separate route to
    the same matrix rather than a restatement of the BFS.
    """
    dist = dict((s, {s: 0}) for s in range(n))
    frontier = dict((s, set([s])) for s in range(n))
    for step in range(1, cap + 1):
        for s in range(n):
            nxt = set()
            for x in frontier[s]:
                nxt |= out_nb[x]
            nxt -= set(dist[s])
            for t in nxt:
                dist[s][t] = step
            frontier[s] = nxt
    return dist


def check_distances_agree(dist_bfs, dist_alt, n):
    complete = [s for s in range(n) if len(dist_bfs[s]) != n]
    ck("every_vertex_reaches_all_104_vertices", not complete,
       "vertices with incomplete reach: %s" % complete[:5])

    mismatch = 0
    for s in range(n):
        if dist_bfs[s] != dist_alt[s]:
            mismatch += 1
    ck("bfs_and_relation_closure_distances_coincide", mismatch == 0,
       "sources with differing distance rows: %d" % mismatch)


def check_radius_diameter(dist, n):
    ecc = dict((s, max(dist[s].values())) for s in range(n))
    multiset = {}
    for s in range(n):
        multiset[ecc[s]] = multiset.get(ecc[s], 0) + 1
    ck("eccentricity_multiset_is_2_once_and_3_a_hundred_and_three_times",
       multiset == ECC_MULTISET_PAPER,
       "computed: %s" % sorted(multiset.items()))

    radius = min(ecc.values())
    diameter = max(ecc.values())
    ck("radius_equals_k_equals_2", radius == K_PAPER, "radius=%d" % radius)
    ck("diameter_equals_k_plus_1_equals_3", diameter == K_PAPER + 1,
       "diameter=%d" % diameter)

    centre = sorted(s for s in range(n) if ecc[s] == radius)
    ck("unique_central_vertex", len(centre) == 1,
       "centre=%s" % centre)

    profile = []
    c = centre[0] if centre else 0
    for layer in range(diameter + 1):
        cnt = sum(1 for t in range(n) if dist[c][t] == layer)
        if cnt:
            profile.append(cnt)
    ck("central_layer_profile_is_1_10_93",
       profile == CENTRAL_PROFILE_PAPER,
       "profile at vertex %d: %s" % (c, profile))
    return ecc, radius, diameter, centre


def check_central_row_as_described(rows, centre):
    """The text says vertex 0 is central, 1..7 are its undirected
    neighbours and 8,9,10 the heads of its arcs."""
    ok_centre = centre == [0]
    u0 = sorted(rows[0][0]) if 0 in rows else []
    a0 = sorted(rows[0][1]) if 0 in rows else []
    ck("vertex_0_is_central_with_neighbours_1_to_7_and_arcs_to_8_9_10",
       ok_centre and u0 == [1, 2, 3, 4, 5, 6, 7] and a0 == [8, 9, 10],
       "centre=%s, U(0)=%s, A+(0)=%s" % (centre, u0, a0))


def check_moore_bound(n):
    """The diameter-two mixed Moore bound, recomputed, not quoted."""
    def moore(r, z):
        return 1 + (r + z) + r * (r - 1 + z) + z * (r + z)

    m = moore(R_PAPER, Z_PAPER)
    ck("moore_bound_M_7_3_2_equals_104_equals_order",
       m == MOORE_BOUND_PAPER and m == n,
       "M(7,3,2)=%d, |V|=%d" % (m, n))

    # The paper's argument is that the order pins down s = r + z via
    # s^2 + 1 <= M(r,z,2) <= s^2 + s + 1.  Recompute the set of s values
    # compatible with each of the two orders in question.
    s_for = {}
    for target in (ORDER_PAPER, ORDER_108):
        s_for[target] = sorted(set(
            r + z for r in range(0, 41) for z in range(0, 41)
            if (r + z) ** 2 + 1 <= target <= (r + z) ** 2 + (r + z) + 1))
    ck("orders_104_and_108_both_force_s_equals_10",
       s_for[ORDER_PAPER] == [10] and s_for[ORDER_108] == [10],
       "s admissible for 104: %s; for 108: %s"
       % (s_for[ORDER_PAPER], s_for[ORDER_108]))

    hits104 = [(r, z) for r in range(0, 41) for z in range(0, 41)
               if moore(r, z) == ORDER_PAPER and r + z > 0]
    hits108 = [(r, z) for r in range(0, 41) for z in range(0, 41)
               if moore(r, z) == ORDER_108 and r + z > 0]
    ck("order_104_forces_r_7_z_3", hits104 == [(R_PAPER, Z_PAPER)],
       "solutions of M(r,z,2)=104: %s" % hits104)
    ck("order_108_forces_r_3_z_7", hits108 == [(3, 7)],
       "solutions of M(r,z,2)=108: %s" % hits108)
    return m


def check_ball_sizes(dist, n, moore_value):
    """The Moore bound is attained by exactly one 2-ball, and the deficiency
    of every other 2-ball is exactly its distance-three count."""
    balls = dict((v, sum(1 for t in range(n) if dist[v][t] <= K_PAPER))
                 for v in range(n))
    q = dict((v, sum(1 for t in range(n) if dist[v][t] == 3)) for v in range(n))
    # Fails as soon as any vertex sits at distance 4 or more from some other,
    # since such a vertex is counted in neither the 2-ball nor q.
    bad = sorted(v for v in range(n) if balls[v] + q[v] != n)
    ck("two_ball_deficiency_equals_distance_three_count", not bad,
       "vertices where |B_2(v)| + q(v) != %d: %s" % (n, bad[:5]))

    tight = sorted(v for v in range(n) if balls[v] == moore_value)
    ck("moore_bound_attained_by_exactly_one_2_ball",
       len(tight) == 1,
       "vertices with |B_2| = %d: %s (max |B_2| = %d)"
       % (moore_value, tight, max(balls.values())))


def check_status_norm(dist, n):
    """sigma(v), the Moore status, the total status and N_1, all recomputed."""
    # The Moore status for RM(r,z,2): 1 vertex at distance 0, r+z at
    # distance 1, the rest at distance 2.  Derived, not quoted.
    s1 = R_PAPER + Z_PAPER
    s2 = n - 1 - s1
    sigma_moore = 0 * 1 + 1 * s1 + 2 * s2
    ck("moore_status_sigma_M_equals_196",
       sigma_moore == SIGMA_MOORE_PAPER,
       "1*%d + 2*%d = %d" % (s1, s2, sigma_moore))

    sigma = dict((v, sum(dist[v].values())) for v in range(n))
    q = dict((v, sum(1 for t in range(n) if dist[v][t] == 3)) for v in range(n))

    # Independent route to the same statuses via the paper's identity
    # sigma(v) = 196 + q(v); a discrepancy at any vertex fails the check.
    bad = sorted(v for v in range(n)
                 if sigma[v] != sigma_moore + q[v])
    ck("status_identity_sigma_v_equals_196_plus_q_v", not bad,
       "vertices where the identity fails: %s" % bad[:5])

    d1 = sorted(set(sum(1 for t in range(n) if dist[v][t] == 1)
                    for v in range(n)))
    ck("exactly_ten_vertices_at_distance_one_from_every_vertex",
       d1 == [s1], "distance-1 counts present: %s" % d1)

    total = sum(sigma.values())
    ck("total_status_is_23762", total == TOTAL_STATUS_PAPER,
       "sum of sigma = %d" % total)

    n1 = sum(abs(sigma[v] - sigma_moore) for v in range(n))
    ck("status_1_norm_N1_equals_3378", n1 == N1_PAPER,
       "N_1 = %d (and total - 104*196 = %d)"
       % (n1, total - n * sigma_moore))

    pairs3 = sum(1 for v in range(n) for t in range(n) if dist[v][t] == 3)
    ck("ordered_pairs_at_distance_three_equals_N1",
       pairs3 == N1_PAPER and pairs3 == n1,
       "pairs at distance 3 = %d" % pairs3)
    return sigma, n1


def check_settlement(rows, und, arc, n, moore_value, radius, diameter,
                     centre):
    """THE LOAD-BEARING CHECK.

    An RM(r,z,k) graph is a totally (r,z)-regular mixed graph of order
    M(r,z,k) with radius k and diameter k+1.  Each conjunct below is
    recomputed from the certificate; the conclusion "an RM(7,3,2) graph
    exists" is the conjunction, not an assertion.
    """
    is_regular = (all(len(und[v]) == R_PAPER for v in rows)
                  and all(len(arc[v]) == Z_PAPER for v in rows))
    indeg = dict((v, 0) for v in rows)
    for v in rows:
        for u in arc[v]:
            indeg[u] += 1
    is_regular = is_regular and all(indeg[v] == Z_PAPER for v in rows)

    simple = (all(v not in und[v] and v not in arc[v] for v in rows)
              and not [1 for v in rows for u in und[v] if v not in und[u]]
              and not [1 for v in rows for u in arc[v] if v in arc[u]]
              and not [1 for v in rows for u in arc[v] if u in und[v]])

    conjuncts = [
        ("simple mixed graph", simple),
        ("totally (7,3)-regular", is_regular),
        ("order = M(7,3,2)", n == moore_value),
        ("radius = 2", radius == K_PAPER),
        ("diameter = 3", diameter == K_PAPER + 1),
    ]
    failed = [name for name, good in conjuncts if not good]
    ck("RM_7_3_2_graph_exists_all_defining_conditions_hold",
       not failed,
       "unmet conditions: %s" % (failed if failed else "none"))

    ck("witness_is_a_radial_not_a_full_moore_graph",
       radius == K_PAPER and diameter == K_PAPER + 1 and len(centre) == 1,
       "radius %d < diameter %d, |centre|=%d"
       % (radius, diameter, len(centre)))


def check_bound_closed_form():
    """The two derivation steps Section 2 prints, recomputed.

    Section 2 goes from the order to the parameter pair via the closed form
    M(r,z,2) = (r+z)^2 + z + 1 and the consequent z = order - s^2 - 1.  Both
    steps are recomputed here over the whole search range.  Note that this
    route uses no row or column headers of any printed table: the pair is
    pinned by the order and the Moore status baseline alone.
    """
    span = range(0, 41)

    def moore(r, z):
        return 1 + (r + z) + r * (r - 1 + z) + z * (r + z)

    bad_form = [(r, z) for r in span for z in span
                if moore(r, z) != (r + z) ** 2 + z + 1]
    ck("expanded_and_closed_forms_of_the_moore_bound_agree", not bad_form,
       "pairs where 1+(r+z)+r(r-1+z)+z(r+z) != (r+z)^2+z+1: %s"
       % bad_form[:5])

    # z = order - s^2 - 1 with s = 10, exactly the step printed in Section 2.
    z108 = ORDER_108 - 10 ** 2 - 1
    r108 = 10 - z108
    z104 = ORDER_PAPER - 10 ** 2 - 1
    r104 = 10 - z104
    ck("directed_degree_recovered_from_the_order_is_7_at_108_and_3_at_104",
       (r108, z108) == (3, 7) and (r104, z104) == (R_PAPER, Z_PAPER)
       and moore(r108, z108) == ORDER_108
       and moore(r104, z104) == ORDER_PAPER,
       "order 108 -> (r,z)=(%d,%d); order 104 -> (r,z)=(%d,%d)"
       % (r108, z108, r104, z104))

    def baseline(order, s):
        return 1 * s + 2 * (order - 1 - s)

    fits108 = [(r, z) for r in span for z in span
               if moore(r, z) == ORDER_108
               and baseline(ORDER_108, r + z) == BASELINE_108_PAPER]
    fits104 = [(r, z) for r in span for z in span
               if moore(r, z) == ORDER_PAPER
               and baseline(ORDER_PAPER, r + z) == BASELINE_108_PAPER]
    ck("length_108_with_baseline_204_pins_r_3_z_7_using_no_table_layout",
       fits108 == [(3, 7)] and fits104 == [],
       "(r,z) of order 108 with baseline 204: %s; of order 104 with "
       "baseline 204: %s" % (fits108, fits104))


def check_parameter_correction():
    """Section 2: the entry N_1 = 18 belongs to (r,z) = (3,7), not (7,3)."""
    def moore(r, z):
        return 1 + (r + z) + r * (r - 1 + z) + z * (r + z)

    vec = []
    for value, mult in STATUS_108_PAPER:
        vec.extend([value] * mult)
    ck("quoted_status_vector_has_108_entries_not_104",
       len(vec) == ORDER_108 and len(vec) != ORDER_PAPER,
       "length %d; order of the exhibited graph is %d"
       % (len(vec), ORDER_PAPER))

    # Baseline for an order-108 diameter-two-Moore-bound graph, derived the
    # same way as sigma_M above: 1 + (r+z) + rest, with (r,z) = (3,7).
    s1 = 3 + 7
    base108 = 1 * s1 + 2 * (ORDER_108 - 1 - s1)
    ck("order_108_moore_status_baseline_is_204",
       base108 == BASELINE_108_PAPER and moore(3, 7) == ORDER_108,
       "baseline %d, M(3,7,2)=%d" % (base108, moore(3, 7)))

    n1_108 = sum(abs(x - base108) for x in vec)
    ck("quoted_vector_has_status_norm_18",
       n1_108 == N1_108_PAPER,
       "8*(205-204) + 2*(209-204) = %d" % n1_108)

    # The transposed cell (r,z) = (7,3) has a different order and a
    # different baseline, so the entry cannot have come from it.
    base104 = 1 * (R_PAPER + Z_PAPER) + 2 * (
        ORDER_PAPER - 1 - (R_PAPER + Z_PAPER))
    ck("the_two_cells_have_different_orders_and_baselines",
       moore(3, 7) != moore(R_PAPER, Z_PAPER)
       and base108 != base104
       and (base108, moore(3, 7)) == (204, 108)
       and (base104, moore(R_PAPER, Z_PAPER)) == (196, 104),
       "(3,7): order %d baseline %d; (7,3): order %d baseline %d"
       % (moore(3, 7), base108, moore(R_PAPER, Z_PAPER), base104))


def main():
    print("Verification of an RM(7,3,2) mixed radial Moore graph on 104")
    print("vertices, recomputed from the printed adjacency certificate.")
    print("")

    try:
        rows, order = parse_certificate(CERTIFICATE)
    except Exception as exc:                       # malformed certificate
        ck("certificate_parses", False, "%s: %s" % (type(exc).__name__, exc))
        return summarize()
    ck("certificate_parses", True,
       "%d rows decoded" % len(rows))

    n = len(rows)
    check_wellformed(rows, order)
    if order != list(range(n)):
        print("NOTE: vertex labels are not 0..n-1; later checks are skipped.")
        return summarize()

    und, arc = check_mixed_graph_axioms(rows)
    out_nb = check_total_regularity(rows, und, arc)

    dist = bfs_distances(out_nb, n)
    dist_alt = closure_distances(out_nb, n, n)
    check_distances_agree(dist, dist_alt, n)

    if any(len(dist[s]) != n for s in range(n)):
        print("NOTE: the mixed graph is not strongly connected; the radius,")
        print("diameter and status checks below are therefore skipped.")
        return summarize()

    ecc, radius, diameter, centre = check_radius_diameter(dist, n)
    check_central_row_as_described(rows, centre)
    moore_value = check_moore_bound(n)
    check_ball_sizes(dist, n, moore_value)
    check_status_norm(dist, n)
    check_settlement(rows, und, arc, n, moore_value, radius, diameter, centre)
    check_bound_closed_form()
    check_parameter_correction()

    print("")
    print("Echo of the decoded object: %d vertices, undirected degree %d,"
          % (n, len(und[0])))
    print("arc out-degree %d, radius %d, diameter %d, centre %s."
          % (len(arc[0]), radius, diameter, centre))
    print("NOT RE-RUN HERE, stated explicitly:")
    print(" (a) the bibliographic claim that a particular printed table")
    print("     carries the entry 18 in the cell (r,z) = (7,3) is a fact")
    print("     about published pages and cannot be checked by this")
    print("     program.  Concretely, the row and column headers of the")
    print("     two disputed tables are NOT reproduced anywhere in this")
    print("     program or in the paper, so nothing here excludes a table")
    print("     whose axes are indexed the other way round, i.e. by")
    print("     (z,r); on that reading the 'transposed cell' correction")
    print("     would itself be a misreading of a correctly printed")
    print("     table, and no check above would detect it.  The same")
    print("     holds for the paper's supporting quotation that the")
    print("     source text names the vector s_{3,7,2}.  What IS checked")
    print("     is the parameter arithmetic, which uses no table layout")
    print("     at all: order 108 with baseline 204 forces undirected")
    print("     degree 3 and directed degree 7, whereas (r,z) = (7,3) has")
    print("     order 104 and baseline 196.")
    print(" (b) no census over an external catalogue of mixed graphs is")
    print("     attempted, and no minimality of N_1 = 3378 is tested;")
    print("     the paper makes no optimality claim either.")
    print(" (c) the exhaustive parameter search above ranges over the")
    print("     pairs (r,z) with 0 <= r,z <= 40, which covers every pair")
    print("     with M(r,z,2) <= 1681 and so both orders in question.")
    print("")
    return summarize()


if __name__ == "__main__":
    sys.exit(main())

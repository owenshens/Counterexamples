#!/usr/bin/env python3
"""Verification of an exhibited (5,2,2) mixed radial Moore graph.

TAKEN FROM THE PAPER (inputs, transcribed):
  - the parameters r=5, z=2, k=2 and the Moore-like bound M(r,z,k);
  - the explicit adjacency data of the exhibited mixed graph on 52 vertices.

DERIVED HERE (computed, nothing asserted):
  - the order equals M(5,2,2)=52;
  - total regularity: every vertex has undirected degree 5, out-degree 2, in-degree 2;
  - the digon-free / simple-mixed conditions;
  - radius exactly 2 and diameter exactly 3 from all-pairs directed distances;
  - hence existence of a (5,2,2) mixed radial Moore graph;
  - the layer sizes 1,7,44, the reference status 95, the edge/arc counts
    130/104, the distance distribution 52/364/1487/801, the total status 5741
    and N_1 = 801;
  - a second computation of all 52x52 = 2704 distances by bitmask closures
    instead of breadth-first search, compared to the search PAIR BY PAIR (this
    re-does the traversal, not the decoding of the table: both arms read the
    same decoded adjacency, as the closing NOT RE-RUN paragraph states);
  - that (5,2) has the least Moore bound among the parameter pairs recorded as
    open, given the documented list of already-settled pairs;
  - a falsification self-test in two parts: (a) local corruptions that break a
    structural condition (symmetry, degrees, edge/arc overlap) are rejected, and
    (b) DEGREE-PRESERVING rewirings -- an arc-head swap and an edge 2-swap --
    which leave the graph totally (5,2)-regular with the same 130 edges and 104
    arcs, and are therefore rejected only by the *metric* conditions (radius,
    diameter, center, distance histogram, sum s(v), N_1). Part (b) is what shows
    the metric checks are capable of failing.

All arithmetic is exact integer arithmetic. Standard library only.
"""
import sys
from collections import deque

R, Z, K = 5, 2, 2          # paper input: parameters (r,z,k)
SIGMA = 95                 # paper input: reference status sigma_{5,2,2}
CLAIM_ORDER = 52           # paper input: M(5,2,2)
CLAIM_EDGES = 130          # paper input
CLAIM_ARCS = 104           # paper input
CLAIM_LAYERS = (1, 7, 44)  # paper input: layer sizes from the central vertex
CLAIM_DIST = {0: 52, 1: 364, 2: 1487, 3: 801}   # paper input
CLAIM_STATUS_SUM = 5741    # paper input
CLAIM_N1 = 801             # paper input
CLAIM_CENTER = [0]         # paper input

# Adjacency certificate transcribed from the paper's table:
#   v : U(v) (undirected neighbours) : A^+(v) (heads of outgoing arcs)
CERT = """
0 : 1,2,3,4,5 : 6,7
1 : 0,8,9,10,11 : 12,13
2 : 0,14,15,16,17 : 18,19
3 : 0,20,21,22,23 : 24,25
4 : 0,26,27,28,29 : 30,31
5 : 0,32,33,34,35 : 36,37
6 : 38,39,40,41,42 : 43,44
7 : 45,46,47,48,49 : 50,51
8 : 1,36,43,48,50 : 14,47
9 : 1,25,33,35,43 : 16,50
10 : 1,21,23,31,44 : 4,33
11 : 1,12,15,32,40 : 23,35
12 : 11,14,16,44,46 : 32,51
13 : 16,25,34,46,50 : 41,45
14 : 2,12,30,35,45 : 26,32
15 : 2,11,23,26,27 : 17,40
16 : 2,12,13,34,44 : 15,41
17 : 2,19,20,41,47 : 24,45
18 : 19,21,33,42,48 : 1,5
19 : 17,18,28,34,50 : 14,42
20 : 3,17,21,28,47 : 33,39
21 : 3,10,18,20,30 : 1,7
22 : 3,26,29,33,41 : 11,28
23 : 3,10,15,24,36 : 0,26
24 : 23,27,38,40,51 : 2,37
25 : 9,13,36,43,49 : 38,39
26 : 4,15,22,30,41 : 31,38
27 : 4,15,24,36,40 : 13,22
28 : 4,19,20,39,51 : 5,34
29 : 4,22,32,33,51 : 2,44
30 : 14,21,26,39,51 : 20,42
31 : 10,32,37,43,44 : 8,34
32 : 5,11,29,31,47 : 9,46
33 : 5,9,18,22,29 : 15,19
34 : 5,13,16,19,39 : 17,22
35 : 5,9,14,46,50 : 18,27
36 : 8,23,25,27,49 : 10,12
37 : 31,38,42,45,46 : 16,30
38 : 6,24,37,43,49 : 0,29
39 : 6,28,30,34,51 : 27,48
40 : 6,11,24,27,50 : 9,29
41 : 6,17,22,26,48 : 21,23
42 : 6,18,37,47,49 : 8,25
43 : 8,9,25,31,38 : 4,48
44 : 10,12,16,31,45 : 43,49
45 : 7,14,37,44,48 : 10,11
46 : 7,12,13,35,37 : 3,6
47 : 7,17,20,32,42 : 3,36
48 : 7,8,18,41,45 : 28,35
49 : 7,25,36,38,42 : 20,47
50 : 8,13,19,35,40 : 21,46
51 : 24,28,29,30,39 : 40,49
"""

CHECKS = []

def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)

def decode(text):
    """Parse the certificate into (vertices, U, A) with U[v], A[v] sorted tuples."""
    U, A, order = {}, {}, []
    for raw in text.strip().splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.split(":")
        if len(parts) != 3:
            raise ValueError("malformed certificate row: %r" % raw)
        v = int(parts[0].strip())
        us = tuple(sorted(int(x) for x in parts[1].split(",") if x.strip()))
        aa = tuple(sorted(int(x) for x in parts[2].split(",") if x.strip()))
        if v in U:
            raise ValueError("duplicate row for vertex %d" % v)
        U[v], A[v] = us, aa
        order.append(v)
    return order, U, A


def moore_bound(r, z, k):
    """Mixed Moore-tree count; for k=2 this is the paper's M(r,z,2)."""
    if k != 2:
        raise ValueError("only k=2 implemented")
    return 1 + (r + z) + r * (r - 1 + z) + z * (r + z)


def moore_layers(r, z):
    """Layer sizes 1, r+z, r(r-1+z)+z(r+z) of the k=2 Moore tree."""
    return (1, r + z, r * (r - 1 + z) + z * (r + z))


def out_neighbours(U, A):
    """Out-neighbourhood U(v) u A^+(v); references outside V are dropped here and
    reported separately by the range check, so a truncated table cannot crash."""
    V = set(U)
    return dict((v, (set(U[v]) | set(A.get(v, ()))) & V) for v in U)


def bfs_dist(out, src):
    """Source-to-target distances: edges either way, arcs forward only."""
    dist = {src: 0}
    q = deque([src])
    while q:
        x = q.popleft()
        for y in out[x]:
            if y not in dist:
                dist[y] = dist[x] + 1
                q.append(y)
    return dist


def all_pairs(out):
    return dict((v, bfs_dist(out, v)) for v in out)


def evaluate(order, U, A):
    """Compute every structural and metric fact used below. No claim is assumed."""
    f = {}
    V = set(order)
    n = len(order)
    f["n"] = n
    f["labels_ok"] = (V == set(range(n)) and len(order) == len(V))
    f["in_range"] = all(0 <= x < n and x in V
                        for v in order for x in tuple(U[v]) + tuple(A[v]))
    f["no_repeat_in_row"] = all(len(set(U[v])) == len(U[v])
                                and len(set(A[v])) == len(A[v]) for v in order)
    f["sym_ok"] = all((v in U[u]) == (u in U[v]) for v in order for u in order)
    f["udeg"] = sorted(len(U[v]) for v in order)
    f["outdeg"] = sorted(len(A[v]) for v in order)
    indeg = dict((v, 0) for v in order)
    for v in order:
        for h in A[v]:
            if h in indeg:
                indeg[h] += 1
    f["indeg"] = sorted(indeg.values())
    f["loops"] = sum(1 for v in order if v in U[v] or v in A[v])
    f["digons"] = sum(1 for v in order for h in A[v] if v in A.get(h, ()))
    f["overlaps"] = sum(1 for v in order for h in A[v]
                        if h in U[v] or v in U.get(h, ()))
    f["edges"] = sum(len(U[v]) for v in order) // 2
    f["edges_even"] = (sum(len(U[v]) for v in order) % 2 == 0)
    f["arcs"] = sum(len(A[v]) for v in order)

    out = out_neighbours(U, A)
    dists = all_pairs(out)
    hist, unreach = {}, 0
    for v in order:
        for u in order:
            if u in dists[v]:
                d = dists[v][u]
                hist[d] = hist.get(d, 0) + 1
            else:
                unreach += 1
    f["hist"] = hist
    f["unreachable"] = unreach
    f["strong"] = (unreach == 0)
    ecc = dict((v, max(dists[v].values())) for v in order) if unreach == 0 else {}
    f["ecc"] = ecc
    if ecc:
        f["radius"] = min(ecc.values())
        f["diameter"] = max(ecc.values())
        f["center"] = sorted(v for v in order if ecc[v] == f["radius"])
        f["status"] = dict((v, sum(dists[v].values())) for v in order)
        f["status_sum"] = sum(f["status"].values())
        f["n1"] = sum(abs(f["status"][v] - SIGMA) for v in order)
        f["q"] = dict((v, sum(1 for u in order if dists[v][u] == 3)) for v in order)
    else:
        f["radius"] = f["diameter"] = None
        f["center"] = []
        f["status"] = {}
        f["status_sum"] = f["n1"] = None
        f["q"] = {}
    f["layer_from_center"] = None
    if 0 in dists:
        lay = {}
        for u, d in dists[0].items():
            lay.setdefault(d, set()).add(u)
        f["layer_from_center"] = lay
    return f


def finish():
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)

def reach_masks(order, U, A):
    """Second, search-independent recomputation: bitmask closures 'reachable in
    <= i steps', i = 1,2,3, by boolean composition rather than breadth-first
    search.  It consumes the SAME decoded adjacency as bfs_dist (both go through
    out_neighbours), so it is an independent traversal, not an independent decode
    of the paper's table; a fault in the decode would corrupt both arms alike.
    That limitation is stated in the closing NOT RE-RUN paragraph."""
    out = out_neighbours(U, A)
    a1 = dict((v, (1 << v) | sum(1 << u for u in out[v])) for v in order)
    def compose(x, y):
        z = {}
        for v in order:
            m, acc = x[v], 0
            u = 0
            while m:
                if m & 1:
                    acc |= y[u]
                m >>= 1
                u += 1
            z[v] = acc
        return z
    a2 = compose(a1, a1)
    a3 = compose(a1, a2)
    return a1, a2, a3



# The paper's full claim, split so the self-test can say WHICH kind of condition
# rejected a corruption.  STRUCTURAL conditions look only at the incidence data;
# METRIC conditions look only at the computed all-pairs directed distances.
STRUCTURAL_NAMES = ("labels", "endpoints_in_range", "no_repeat_in_row",
                    "symmetry", "undirected_degree_5", "outdegree_2",
                    "indegree_2", "no_loops", "no_digons",
                    "no_edge_arc_overlap", "order_M(5,2,2)", "edges_130",
                    "arcs_104")
METRIC_NAMES = ("strongly_connected", "radius_2", "diameter_3", "center_{0}",
                "distance_histogram", "status_sum_5741", "N_1_801")


def violations(order, U, A):
    """Return (structural_failures, metric_failures, facts) for this mixed graph
    against the paper's full claim.  Two empty lists means it conforms.  Nothing
    is assumed: every entry is decided from the recomputed facts in evaluate(),
    which are handed back so the caller can PRINT measured quantities rather than
    restate the claimed ones."""
    try:
        f = evaluate(order, U, A)
    except Exception as exc:
        return ["evaluation_raised_%s" % type(exc).__name__], [], {}
    n = f["n"]
    s = []
    if not f["labels_ok"]:
        s.append("labels")
    if not f["in_range"]:
        s.append("endpoints_in_range")
    if not f["no_repeat_in_row"]:
        s.append("no_repeat_in_row")
    if not f["sym_ok"]:
        s.append("symmetry")
    if f["udeg"] != [R] * n:
        s.append("undirected_degree_5")
    if f["outdeg"] != [Z] * n:
        s.append("outdegree_2")
    if f["indeg"] != [Z] * n:
        s.append("indegree_2")
    if f["loops"] != 0:
        s.append("no_loops")
    if f["digons"] != 0:
        s.append("no_digons")
    if f["overlaps"] != 0:
        s.append("no_edge_arc_overlap")
    if n != moore_bound(R, Z, K):
        s.append("order_M(5,2,2)")
    if not (f["edges_even"] and f["edges"] == CLAIM_EDGES):
        s.append("edges_130")
    if f["arcs"] != CLAIM_ARCS:
        s.append("arcs_104")
    m = []
    if not f["strong"]:
        m.append("strongly_connected")
    if f["radius"] != K:
        m.append("radius_2")
    if f["diameter"] != K + 1:
        m.append("diameter_3")
    if f["center"] != CLAIM_CENTER:
        m.append("center_{0}")
    if f["hist"] != CLAIM_DIST:
        m.append("distance_histogram")
    if f["status_sum"] != CLAIM_STATUS_SUM:
        m.append("status_sum_5741")
    if f["n1"] != CLAIM_N1:
        m.append("N_1_801")
    return s, m, f


def _swap_arc_heads(order, A, v, w):
    """Give v one of w's arc heads and w one of v's.  The multiset of heads is
    unchanged, so every in-degree is preserved; out-degrees stay 2 and the whole
    undirected part is untouched.  Only the metric can notice."""
    a = dict((x, list(A[x])) for x in order)
    hv, hw = sorted(A[v]), sorted(A[w])
    a[v] = [hv[0], hw[1]]
    a[w] = [hw[0], hv[1]]
    return a


def _swap_edges(order, U, a, b, c, d):
    """Replace the edges {a,b} and {c,d} by {a,d} and {c,b}.  Every undirected
    degree is preserved, the edge count is preserved, and no arc is touched."""
    u = dict((x, list(U[x])) for x in order)
    for x, y in ((a, b), (c, d)):
        u[x].remove(y)
        u[y].remove(x)
    for x, y in ((a, d), (c, b)):
        u[x].append(y)
        u[y].append(x)
    return u


def _dp_detail(nm, sf, mf, f2):
    """Detail line for a degree-preserving mutant.  Every number here is read out
    of that mutant's own recomputed facts, so the printed 'still (5,2)-regular
    with 130 edges and 104 arcs' is a MEASUREMENT and not a restatement of the
    paper's claim."""
    def ends(seq):
        return "%s..%s" % (seq[0], seq[-1]) if seq else "?"
    return ("%s measures udeg %s, outdeg %s, indeg %s, %s edges, %s arcs, "
            "loops/digons/overlaps %s/%s/%s [structural violations: %s] "
            "yet fails %s"
            % (nm, ends(f2.get("udeg")), ends(f2.get("outdeg")),
               ends(f2.get("indeg")), f2.get("edges", "?"), f2.get("arcs", "?"),
               f2.get("loops", "?"), f2.get("digons", "?"),
               f2.get("overlaps", "?"),
               ",".join(sf) or "none", ",".join(mf) or "NOTHING"))


def mutants(order, U, A):
    """Corruptions of the exhibited object, for the falsification self-test.

    The first five break a structural condition.  The last two are
    DEGREE-PRESERVING rewirings: they keep total (5,2)-regularity, 130 edges and
    104 arcs, no loops/digons/overlaps, so the only conditions left to reject
    them are the metric ones.  Each is tagged with the class of condition it is
    expected to be caught by, and the test below enforces the tag."""
    out = []
    u2 = dict((v, list(U[v])) for v in order)
    u2[0] = [1, 2, 3, 4, 6]          # breaks undirected symmetry / degrees
    out.append(("edge-rewire", "structural", order, u2, A))
    a2 = dict((v, list(A[v])) for v in order)
    a2[0] = [6, 8]                    # arc head moved: in-degrees change
    out.append(("arc-head-moved", "structural", order, U, a2))
    a3 = dict((v, list(A[v])) for v in order)
    a3[0] = [6, 1]                    # arc parallel to an existing edge
    out.append(("edge-arc-overlap", "structural", order, U, a3))
    u4 = dict((v, list(U[v])) for v in order)
    u4[8] = list(U[8]) + [0]
    u4[0] = list(U[0]) + [8]          # degree 6 at two vertices, still symmetric
    out.append(("degree-bumped", "structural", order, u4, A))
    u5 = dict((v, list(U[v])) for v in order)
    u5[13] = [x for x in U[13] if x != 16]
    u5[16] = [x for x in U[16] if x != 13]   # degree 4 at two vertices, 129 edges
    out.append(("edge-deleted", "structural", order, u5, A))
    # Degree-preserving: 0 and 1 exchange one arc head each, so 0 -> {6,13} and
    # 1 -> {12,7}.  All degrees, loop/digon/overlap counts and both totals are
    # untouched; the second out-neighbourhood of 0 loses {47,48,49,51}.
    out.append(("arc-heads-swapped-0-1", "metric", order, U,
                _swap_arc_heads(order, A, 0, 1)))
    # Degree-preserving: the edges {0,1} and {6,39} become {0,39} and {6,1}.
    # Again no degree, loop, digon, overlap or count changes; the block {8..13}
    # is pushed out of the second out-neighbourhood of 0.
    out.append(("edges-2-swapped-01-639", "metric", order,
                _swap_edges(order, U, 0, 1, 6, 39), A))
    return out


def main():
    order, U, A = decode(CERT)
    f = evaluate(order, U, A)
    n = f["n"]
    print("decoded certificate: %d rows, %d undirected endpoints, %d arc heads"
          % (n, sum(len(U[v]) for v in order), sum(len(A[v]) for v in order)))
    lo, hi = order[0], order[-1]
    print("first row v=%d: U=%s A+=%s ; last row v=%d: U=%s A+=%s"
          % (lo, list(U[lo]), list(A[lo]), hi, list(U[hi]), list(A[hi])))

    ck("certificate_well_formed",
       f["labels_ok"] and f["in_range"] and f["no_repeat_in_row"] and n == 52,
       "vertex set = {0,...,%d}, no repeats, all endpoints in range" % (n - 1))

    mb = moore_bound(R, Z, K)
    ck("order_equals_mixed_moore_bound",
       n == mb == CLAIM_ORDER and mb == (R + Z) ** 2 + Z + 1,
       "|V|=%d, M(%d,%d,%d)=%d" % (n, R, Z, K, mb))

    ck("undirected_part_symmetric", f["sym_ok"],
       "u in U(v) iff v in U(u) over all %d ordered pairs" % (n * n))

    ck("totally_5_2_regular",
       f["udeg"] == [R] * n and f["outdeg"] == [Z] * n and f["indeg"] == [Z] * n,
       "undirected degree %d, outdegree %d, indegree %d at every vertex" % (R, Z, Z))

    ck("no_loops_digons_or_edge_arc_overlap",
       f["loops"] == 0 and f["digons"] == 0 and f["overlaps"] == 0,
       "loops=%d digons=%d overlaps=%d" % (f["loops"], f["digons"], f["overlaps"]))

    ck("edge_and_arc_counts",
       f["edges_even"] and f["edges"] == CLAIM_EDGES and f["arcs"] == CLAIM_ARCS,
       "%d edges, %d arcs" % (f["edges"], f["arcs"]))

    lay = f["layer_from_center"] or {}
    sizes = tuple(len(lay.get(i, ())) for i in range(3))
    ck("layers_from_vertex_0",
       sizes == CLAIM_LAYERS and sum(sizes) == n
       and lay.get(1) == set(range(1, 8)),
       "layer sizes %s, first layer = {1,...,7}" % (sizes,))

    blocks = [(1, range(8, 14)), (2, range(14, 20)), (3, range(20, 26)),
              (4, range(26, 32)), (5, range(32, 38)),
              (6, range(38, 45)), (7, range(45, 52))]
    seen, ok_blocks = set(), True
    for v, rng in blocks:
        desc = (set(U.get(v, ())) | set(A.get(v, ()))) - set([0])
        if desc != set(rng) or seen & desc:
            ok_blocks = False
        seen |= desc
    ck("second_layer_block_partition",
       ok_blocks and seen == set(range(8, 52)),
       "the 7 descendant blocks are disjoint and exhaust {8,...,51}")

    hist = f["hist"]
    ck("all_pairs_distance_distribution",
       f["unreachable"] == 0 and hist == CLAIM_DIST and sum(hist.values()) == n * n,
       "d=0:%d d=1:%d d=2:%d d=3:%d, total %d ordered pairs"
       % (hist.get(0, 0), hist.get(1, 0), hist.get(2, 0), hist.get(3, 0),
          sum(hist.values())))

    ecc = f["ecc"]
    ck("out_eccentricities_radius_diameter_center",
       ecc.get(0) == K and all(ecc.get(v) == K + 1 for v in order if v != 0)
       and f["radius"] == K and f["diameter"] == K + 1
       and f["center"] == CLAIM_CENTER,
       "ecc+(0)=%s, ecc+(v)=%d for v!=0, radius %s, diameter %s, center %s"
       % (ecc.get(0), K + 1, f["radius"], f["diameter"], f["center"]))

    ck("is_RM_5_2_2",
       f["udeg"] == [R] * n and f["outdeg"] == [Z] * n and f["indeg"] == [Z] * n
       and f["radius"] == K and f["diameter"] == K + 1 and n == moore_bound(R, Z, K),
       "totally (%d,%d)-regular, radius %d, diameter %d, order %d = M(%d,%d,%d)"
       % (R, Z, f["radius"], f["diameter"], n, R, Z, K))

    ml = moore_layers(R, Z)
    sigma = sum(i * ml[i] for i in range(3))
    ck("moore_layer_sizes_and_sigma",
       sum(ml) == moore_bound(R, Z, K) and ml == CLAIM_LAYERS and sigma == SIGMA,
       "Moore layers %s sum to %d, sigma = %d" % (ml, sum(ml), sigma))

    st, q = f["status"], f["q"]
    ck("status_identity_s_equals_sigma_plus_q",
       bool(q) and all(st.get(v) == SIGMA + q[v] for v in order if v in q)
       and len(q) == n,
       "s(v) = %d + q(v) at all %d vertices, q range [%s,%s]"
       % (SIGMA, n, min(q.values()) if q else "-", max(q.values()) if q else "-"))

    ck("total_status_and_N1",
       f["status_sum"] == CLAIM_STATUS_SUM and f["n1"] == CLAIM_N1
       and f["n1"] == sum(q.values()) and f["status_sum"] == n * SIGMA + f["n1"],
       "sum s(v) = %d, N_1 = %d" % (f["status_sum"], f["n1"]))

    a1, a2, a3 = reach_masks(order, U, A)
    full = (1 << n) - 1
    pop = lambda m: bin(m).count("1")
    c1 = sum(pop(a1[v]) for v in order)
    c2 = sum(pop(a2[v]) for v in order)
    c3 = sum(pop(a3[v]) for v in order)
    cum = [hist.get(0, 0), hist.get(0, 0) + hist.get(1, 0)]
    cum.append(cum[1] + hist.get(2, 0))
    cum.append(cum[2] + hist.get(3, 0))
    # PER-PAIR agreement, not merely agreement of the three cumulative totals:
    # since d(v,u) <= i exactly when u lies in the i-step closure of v, every one
    # of the n*n distances can be read straight off the masks.  Compare that
    # matrix entry by entry with the breadth-first matrix, so two compensating
    # errors cannot cancel inside a total.
    bfs_dists = all_pairs(out_neighbours(U, A))
    pair_mismatch = []
    for v in order:
        for u in order:
            if u == v:
                dm = 0
            elif (a1[v] >> u) & 1:
                dm = 1
            elif (a2[v] >> u) & 1:
                dm = 2
            elif (a3[v] >> u) & 1:
                dm = 3
            else:
                dm = None                      # not reachable within 3 steps
            if dm != bfs_dists[v].get(u):
                pair_mismatch.append((v, u, dm, bfs_dists[v].get(u)))
    pairs_agree = n * n - len(pair_mismatch)
    ck("reachability_closure_matches_bfs",
       (c1, c2, c3) == (cum[1], cum[2], cum[3])
       and a3 == dict((v, full) for v in order)
       and a2[0] == full
       and all(a2[v] != full for v in order if v != 0)
       and not pair_mismatch and pairs_agree == n * n,
       "|<=1|=%d |<=2|=%d |<=3|=%d; step-2 closure is complete only from vertex 0; "
       "mask-derived d(v,u) equals the BFS value on %d/%d ordered pairs "
       "(%d mismatches%s)"
       % (c1, c2, c3, pairs_agree, n * n, len(pair_mismatch),
          "" if not pair_mismatch else ": " + str(pair_mismatch[:5])))

    settled = [(r, 1) for r in range(1, 40)] + [(1, z) for z in range(1, 40)]
    settled += [(2, 2), (2, 3), (3, 2), (3, 3), (4, 2), (2, 4)]
    settled = set(settled)
    openp = [(r, z) for r in range(1, 40) for z in range(1, 40)
             if (r, z) not in settled]
    best = min(openp, key=lambda p: (moore_bound(p[0], p[1], 2), p))
    ties = [p for p in openp if moore_bound(p[0], p[1], 2) == moore_bound(R, Z, 2)]
    ck("least_moore_bound_among_open_pairs",
       best == (R, Z) and ties == [(R, Z)]
       and moore_bound(R, Z, 2) == CLAIM_ORDER,
       "min M(r,z,2) over open (r,z) with r,z<=39 is %d, attained only at %s"
       % (moore_bound(*best, k=2), best))

    base_s, base_m, _base_f = violations(order, U, A)
    caught, dp = [], []
    for name, kind, o2, U2, A2 in mutants(order, U, A):
        sfail, mfail, f2 = violations(o2, U2, A2)
        caught.append((name, kind, bool(sfail or mfail), sfail, mfail))
        if kind == "metric":
            dp.append((name, sfail, mfail, f2))
    vocab_ok = all(x in STRUCTURAL_NAMES for _, _, _, sf, _ in caught for x in sf) \
        and all(x in METRIC_NAMES for _, _, _, _, mf in caught for x in mf)
    # The tag is enforced in BOTH directions: a mutant declared structural must be
    # rejected by a structural condition (not merely by some metric side effect),
    # and the metric ones are held to the converse by the next check.
    tags_ok = all(bool(sf) for _, kd, _, sf, _ in caught if kd == "structural")
    ck("corruptions_are_detected",
       not base_s and not base_m and vocab_ok and tags_ok
       and all(c for _, _, c, _, _ in caught),
       "original conforms (0 structural, 0 metric violations); %d/%d corruptions "
       "rejected (%s)"
       % (sum(1 for _, _, c, _, _ in caught if c), len(caught),
          ", ".join("%s:%s" % (nm, "+".join(sf + mf) or "NOT-REJECTED")
                    for nm, _, _, sf, mf in caught)))

    # Force for the metric block: the two degree-preserving rewirings must pass
    # every structural condition and still be rejected, which can only happen on
    # a distance-derived quantity.  Without this the self-test above would only
    # ever exercise the degree/symmetry/overlap predicates.  Every quantity in the
    # detail below is the mutant's own RECOMPUTED value, not the claimed one.
    ck("metric_checks_can_fail",
       len(dp) == 2 and all(not sf and mf for _, sf, mf, _ in dp),
       "; ".join(_dp_detail(nm, sf, mf, f2) for nm, sf, mf, f2 in dp))
    print("NOT RE-RUN: that the (5,2,2) case is listed as open in the two cited "
          "tables is bibliographic and is not recomputed here; the previously "
          "settled families (z=1 any r, r=1 any z) and the six pairs "
          "(2,2),(2,3),(3,2),(3,3),(4,2),(2,4) are taken as documented input to "
          "the least-Moore-bound comparison; no exhaustive search over all mixed "
          "graphs of order 52 is attempted, and no minimality of N_1 is tested "
          "(the paper makes no such claim). The falsification self-test is seven "
          "named local corruptions, five structural and two degree-preserving; it "
          "is not an exhaustive perturbation search, so it shows each condition "
          "family CAN fail, not that no corruption whatsoever escapes. The second "
          "computation of the distances (the bitmask closures) agrees with the "
          "breadth-first matrix on every one of the 2704 ordered pairs, but it "
          "reads the SAME decoded adjacency, so it cross-checks the traversal and "
          "not the transcription: an error in decoding the paper's table would "
          "corrupt both arms identically and is caught only by the printed first "
          "and last rows and the well-formedness, degree and count checks.")
    return order, U, A, f


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:                     # a damaged table must not crash out
        ck("certificate_parses_and_evaluates", False,
           "%s: %s" % (type(exc).__name__, exc))
    finish()

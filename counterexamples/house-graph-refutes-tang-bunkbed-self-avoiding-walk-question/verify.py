#!/usr/bin/env python3
"""Verification of a refutation of a self-avoiding-walk question on bunkbed graphs.

TAKEN FROM THE PAPER (inputs, not verified here):
  * the statement of the question: for every finite connected simple graph G and
    every edge uv of G that is not a cut edge, is the number of self-avoiding
    walks (paths) in the bunkbed graph G x K_2 from u_0 to v_0 at most the
    number from u_0 to v_1?
  * the exhibited counterexample: the five-vertex "house" graph
    H on {a,b,c,d,e} and the cycle edge ad.
  * the reported counts |(a_0,d_0)| = 60 and |(a_0,d_1)| = 59.

DERIVED HERE (recomputed from scratch by this program):
  * the house graph is decoded, printed back, and checked to be simple,
    connected, five-vertex, six-edge;
  * ad is an edge of H and is not a cut edge (bridge) of H;
  * the bunkbed graph H x K_2 is built from H (10 vertices, 17 edges) and
    checked to be connected;
  * the numbers of self-avoiding walks a_0 -> d_0 and a_0 -> d_1 are counted by
    independent exhaustive enumeration (DFS over paths, plus a second,
    structurally different counter) and compared to the reported 60 and 59;
  * the strict inequality 60 > 59 that refutes the question is COMPUTED;
  * an independent census: every connected simple graph of order at most 5 is
    generated up to isomorphism by brute force over edge subsets plus a
    canonical form, order 6 by adding one vertex to each order-5 class, and
    Delta_G(u,v) is evaluated at every adjacent non-cut edge.  This reproduces
    the paper's order <= 4 table, shows no violation below order five, and
    rediscovers the house/ad pair with counts 60 and 59 without being told
    where to look.  Every violating pair found at orders 5 and 6 is recounted
    by explicit enumeration.

NOT RE-RUN: no census beyond order six is attempted, so the program does not
search for counterexamples on seven or more vertices; the paper's minimality
claim concerns orders at most four and is settled completely here.
"""

import sys
from itertools import combinations, permutations

CHECKS = []

# ------------------------------------------------------------------ #
# Values taken verbatim from the paper.                              #
# ------------------------------------------------------------------ #
HOUSE_VERTEX_NAMES = ["a", "b", "c", "d", "e"]
HOUSE_EDGE_NAMES = ["ab", "ad", "ae", "bc", "be", "cd"]
HOUSE_CYCLE = "abcda"          # the 4-cycle a b c d a
HOUSE_TRIANGLE = "abea"        # the triangle a b e a
PAPER_COUNT_D0 = 60
PAPER_COUNT_D1 = 59
PAPER_PROFILE_D0 = [1, 0, 2, 1, 9, 14, 13, 13, 7]   # N_k(a_0,d_0), k = 1..9
PAPER_PROFILE_D1 = [0, 2, 0, 7, 9, 12, 16, 9, 4]    # N_k(a_0,d_1), k = 1..9
PAPER_SMALL_TABLE = {           # |V| <= 4 : name -> max Delta, or None for "none"
    "K_2": None,
    "P_3": None,
    "K_3": -1,
    "P_4": None,
    "K_{1,3}": None,
    "C_4": -1,
    "paw": 0,
    "K_4-e": -4,
    "K_4": -14,
}
PAPER_MIN_ORDER = 5
PAPER_SEVEN_K4 = [0, 1, 2, 6]
PAPER_SEVEN_K3 = [3, 4, 5]
PAPER_SEVEN_LINKS = [(0, 5), (1, 4), (2, 3)]
PAPER_SEVEN_COUNT_0 = 5031
PAPER_SEVEN_COUNT_1 = 4801
PAPER_SEVEN_CONNECTIVITY = 3


def ck(name, ok, detail=""):
    """Record and print one named check."""
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return bool(ok)


# ------------------------------------------------------------------ #
# Elementary graph machinery (exact integer / bitmask arithmetic).    #
# ------------------------------------------------------------------ #
def ek(x, y):
    """Canonical undirected edge."""
    return (x, y) if x < y else (y, x)


def adj_masks(n, edges):
    """Neighbourhood bitmasks."""
    out = [0] * n
    for (x, y) in edges:
        out[x] |= 1 << y
        out[y] |= 1 << x
    return out


def is_simple(n, edges):
    """No loops, no repeats, endpoints in range, stored canonically."""
    for (x, y) in edges:
        if not (0 <= x < n and 0 <= y < n) or x == y or x > y:
            return False
    return len(set(edges)) == len(edges)


def is_connected(n, edges, drop=0):
    """Connectivity of the subgraph induced on the vertices outside bitmask drop."""
    live = [v for v in range(n) if not (drop >> v) & 1]
    if len(live) <= 1:
        return True
    nbr = adj_masks(n, edges)
    seen = 1 << live[0]
    stack = [live[0]]
    while stack:
        v = stack.pop()
        cand = nbr[v] & ~seen & ~drop
        while cand:
            b = cand & -cand
            seen |= b
            stack.append(b.bit_length() - 1)
            cand ^= b
    return all((seen >> v) & 1 for v in live)


def is_cut_edge(n, edges, e):
    """True when deleting e disconnects the graph."""
    rest = [f for f in edges if f != e]
    return not is_connected(n, rest)


def degrees(n, edges):
    d = [0] * n
    for (x, y) in edges:
        d[x] += 1
        d[y] += 1
    return d


def vertex_connectivity_at_least(n, edges, k):
    """True iff no set of fewer than k vertices disconnects the graph."""
    if not is_connected(n, edges):
        return False
    for r in range(1, k):
        for S in combinations(range(n), r):
            drop = 0
            for v in S:
                drop |= 1 << v
            if not is_connected(n, edges, drop):
                return False
    return True


def bunkbed(n, edges):
    """G x K_2 with x_0 = x and x_1 = x + n."""
    be = []
    for (x, y) in edges:
        be.append(ek(x, y))
        be.append(ek(x + n, y + n))
    for x in range(n):
        be.append(ek(x, x + n))
    return 2 * n, sorted(be)


def paths_by_length(n, edges, s):
    """Last-step recurrence of the paper.

    Returns out[t][k] = number of self-avoiding walks s -> t using k edges,
    computed from F_s(A,x) = number of such walks with vertex set exactly A
    and final vertex x.  A is a bitmask, so |A| - 1 is the walk length and the
    length profile falls out of the same table.
    """
    nbr = adj_masks(n, edges)
    size = 1 << n
    f = [None] * size
    f[1 << s] = [0] * n
    f[1 << s][s] = 1
    out = [[0] * (n + 1) for _ in range(n)]
    for mask in range(size):
        row = f[mask]
        if row is None:
            continue
        k = bin(mask).count("1") - 1
        for x in range(n):
            c = row[x]
            if not c:
                continue
            out[x][k] += c
            cand = nbr[x] & ~mask
            while cand:
                b = cand & -cand
                cand ^= b
                nxt = mask | b
                if f[nxt] is None:
                    f[nxt] = [0] * n
                f[nxt][b.bit_length() - 1] += c
        f[mask] = None       # release; masks are visited in increasing order
    return out


def enumerate_paths(n, edges, s, t):
    """Independent counter: explicit depth-first listing of every s->t path."""
    nbr = adj_masks(n, edges)
    found = []
    path = [s]

    def rec(v, mask):
        if v == t:
            found.append(tuple(path))
            return
        cand = nbr[v] & ~mask
        while cand:
            b = cand & -cand
            cand ^= b
            u = b.bit_length() - 1
            path.append(u)
            rec(u, mask | b)
            path.pop()

    if s == t:
        return [(s,)]
    rec(s, 1 << s)
    return found


def walk_is_valid(n, edges, walk, s, t):
    """A self-avoiding walk: distinct vertices, consecutive pairs are edges."""
    if len(walk) < 1 or walk[0] != s or walk[-1] != t:
        return False
    if len(set(walk)) != len(walk):
        return False
    eset = set(edges)
    for i in range(len(walk) - 1):
        if ek(walk[i], walk[i + 1]) not in eset:
            return False
    return True


def canon(n, edges):
    """Canonical form: lexicographically least relabelled edge list."""
    best = None
    for p in permutations(range(n)):
        cur = tuple(sorted(ek(p[x], p[y]) for (x, y) in edges))
        if best is None or cur < best:
            best = cur
    return best


def connected_classes_brute(n):
    """All connected simple graphs on n labelled vertices, up to isomorphism."""
    pairs = list(combinations(range(n), 2))
    seen = {}
    for r in range(max(n - 1, 0), len(pairs) + 1):
        for sub in combinations(pairs, r):
            if not is_connected(n, sub):
                continue
            c = canon(n, sub)
            if c not in seen:
                seen[c] = tuple(sorted(sub))
    return [seen[c] for c in sorted(seen)]


def connected_classes_grow(n, smaller):
    """Grow order n from order n-1 by adding one vertex joined to a nonempty set.

    Every connected graph has a vertex whose deletion keeps it connected, so
    this reaches every isomorphism class of connected graphs on n vertices.
    """
    seen = {}
    for base in smaller:
        for r in range(1, n):
            for S in combinations(range(n - 1), r):
                edges = tuple(sorted(list(base) + [ek(v, n - 1) for v in S]))
                c = canon(n, edges)
                if c not in seen:
                    seen[c] = edges
    return [seen[c] for c in sorted(seen)]


def delta_table(n, edges):
    """Delta_G(u,v) for every ordered pair (u,v) with uv an adjacent non-cut edge."""
    bn, be = bunkbed(n, edges)
    cache = {}
    res = []
    for e in edges:
        if is_cut_edge(n, edges, e):
            continue
        for (u, v) in (e, (e[1], e[0])):
            if u not in cache:
                cache[u] = paths_by_length(bn, be, u)
            out = cache[u]
            res.append((u, v, sum(out[v]), sum(out[v + n])))
    return res


def max_delta(n, edges):
    """Maximum of Delta over adjacent non-cut edges, or None when there are none."""
    tab = delta_table(n, edges)
    if not tab:
        return None
    return max(d0 - d1 for (_u, _v, d0, d1) in tab)


def note(text):
    print("note: " + text)


def name_of(x):
    return HOUSE_VERTEX_NAMES[x % 5] + str(x // 5)


def decode_house():
    """Turn the paper's edge names into indices; return (n, edges)."""
    idx = {c: i for i, c in enumerate(HOUSE_VERTEX_NAMES)}
    edges = sorted(ek(idx[s[0]], idx[s[1]]) for s in HOUSE_EDGE_NAMES)
    return 5, edges


def cycle_edges(word):
    """Edges of the closed walk spelled by word (first letter repeated at the end)."""
    idx = {c: i for i, c in enumerate(HOUSE_VERTEX_NAMES)}
    return {ek(idx[word[i]], idx[word[i + 1]]) for i in range(len(word) - 1)}


def check_house_wellformed(n, edges):
    """Check 1: the exhibited object decodes, is simple, and is the house graph."""
    deg = degrees(n, edges)
    note("house graph V = " + ",".join(HOUSE_VERTEX_NAMES))
    note("house graph E = " + ",".join(
        HOUSE_VERTEX_NAMES[x] + HOUSE_VERTEX_NAMES[y] for (x, y) in edges))
    note("house degrees " + ", ".join(
        "%s:%d" % (HOUSE_VERTEX_NAMES[v], deg[v]) for v in range(n)))
    union = cycle_edges(HOUSE_CYCLE) | cycle_edges(HOUSE_TRIANGLE)
    ok = (n == 5 and len(edges) == 6 and is_simple(n, edges)
          and set(edges) == union and sorted(deg) == [2, 2, 2, 3, 3])
    return ck("house_graph_wellformed", ok,
              "n=%d m=%d simple=%s equals_cycle_union_triangle=%s"
              % (n, len(edges), is_simple(n, edges), set(edges) == union))


def check_house_two_connected(n, edges):
    """Check 2: hypothesis 'finite connected simple' plus the paper's 2-connectivity."""
    conn = is_connected(n, edges)
    two = vertex_connectivity_at_least(n, edges, 2)
    three = vertex_connectivity_at_least(n, edges, 3)
    survivors = []
    for v in range(n):
        survivors.append((HOUSE_VERTEX_NAMES[v], is_connected(n, edges, 1 << v)))
    ok = conn and two and not three and min(degrees(n, edges)) == 2
    return ck("house_is_2_connected_exactly", ok,
              "connected=%s kappa>=2:%s kappa>=3:%s single-deletions=%s"
              % (conn, two, three,
                 ",".join("%s:%s" % (a, b) for a, b in survivors)))


def check_ad_non_cut(n, edges):
    """Check 3: ad is an edge and is not a cut edge (the statement's hypothesis)."""
    idx = {c: i for i, c in enumerate(HOUSE_VERTEX_NAMES)}
    e = ek(idx["a"], idx["d"])
    present = e in set(edges)
    cut = is_cut_edge(n, edges, e) if present else True
    cuts = [f for f in edges if is_cut_edge(n, edges, f)]
    ok = present and not cut and len(cuts) == 0
    return ck("edge_ad_is_adjacent_and_not_a_cut_edge", ok,
              "ad_present=%s ad_is_cut_edge=%s total_cut_edges=%d"
              % (present, cut, len(cuts)))


def check_bunkbed_wellformed(n, edges, bn, be):
    """Check 4: G x K_2 has the right size, both levels are copies of G, matching is perfect."""
    lower = sorted(e for e in be if e[0] < n and e[1] < n)
    upper = sorted(ek(x - n, y - n) for (x, y) in be if x >= n and y >= n)
    rungs = sorted(e for e in be if e[0] < n <= e[1])
    deg = degrees(bn, be)
    basedeg = degrees(n, edges)
    ok = (bn == 2 * n and len(be) == 2 * len(edges) + n
          and is_simple(bn, be) and is_connected(bn, be)
          and lower == sorted(edges) and upper == sorted(edges)
          and rungs == [(x, x + n) for x in range(n)]
          and all(deg[x] == basedeg[x % n] + 1 for x in range(bn)))
    note("bunkbed vertices " + ",".join(name_of(x) for x in range(bn)))
    return ck("bunkbed_graph_wellformed", ok,
              "n=%d m=%d levels_copy_G=%s rungs=%d connected=%s"
              % (bn, len(be), lower == sorted(edges) == upper, len(rungs),
                 is_connected(bn, be)))


def check_counts(bn, be, s, t, claimed, tag):
    """Checks 5-6: exhaustive listing of the self-avoiding walks s -> t."""
    walks = enumerate_paths(bn, be, s, t)
    ok = len(walks) == claimed and len(set(walks)) == len(walks)
    ck("saw_count_%s_matches_paper" % tag, ok,
       "enumerated=%d distinct=%d paper=%d" % (len(walks), len(set(walks)), claimed))
    return walks


def check_walks_valid(bn, be, walks_by_target):
    """Check 7: every listed object really is a self-avoiding walk with the right ends."""
    bad = 0
    total = 0
    for (s, t, walks) in walks_by_target:
        for w in walks:
            total += 1
            if not walk_is_valid(bn, be, w, s, t):
                bad += 1
    shortest = min((len(w) for (_s, _t, ws) in walks_by_target for w in ws),
                   default=0)
    longest = max((len(w) for (_s, _t, ws) in walks_by_target for w in ws),
                  default=0)
    return ck("enumerated_walks_are_valid", bad == 0 and total > 0,
              "checked=%d invalid=%d vertices_per_walk=%d..%d"
              % (total, bad, shortest, longest))


def check_recurrence_agrees(bn, be, s, targets, walks_by_target):
    """Check 8: the paper's subset recurrence reproduces the enumeration."""
    out = paths_by_length(bn, be, s)
    rows = []
    ok = True
    for (t, walks) in zip(targets, walks_by_target):
        dp = sum(out[t])
        if dp != len(walks[2]):
            ok = False
        rows.append("%s:dp=%d enum=%d" % (name_of(t), dp, len(walks[2])))
    return ck("recurrence_agrees_with_enumeration", ok, " ".join(rows))


def check_profile(bn, be, s, t, walks, claimed, tag):
    """Checks 9-10: the paper's N_k row, from the recurrence AND from the listing."""
    dp = paths_by_length(bn, be, s)[t]
    dp_row = [dp[k] for k in range(1, len(claimed) + 1)]
    enum_row = [0] * len(claimed)
    spill = 0
    for w in walks:
        k = len(w) - 1
        if 1 <= k <= len(claimed):
            enum_row[k - 1] += 1
        else:
            spill += 1
    ok = (dp_row == claimed and enum_row == claimed and spill == 0
          and sum(dp) == sum(claimed))
    return ck("length_profile_%s" % tag, ok,
              "recurrence=%s enumeration=%s paper=%s outside_k=1..%d:%d total=%d"
              % (dp_row, enum_row, claimed, len(claimed), spill, sum(claimed)))


def check_refutation(bn, be, s, t0, t1):
    """Check 11: the load-bearing computation, |SAW(a_0,d_0)| > |SAW(a_0,d_1)|."""
    out = paths_by_length(bn, be, s)
    c0 = sum(out[t0])
    c1 = sum(out[t1])
    delta = c0 - c1
    ok = (c0 == PAPER_COUNT_D0 and c1 == PAPER_COUNT_D1 and c0 > c1
          and delta == 1)
    return ck("questioned_inequality_is_violated", ok,
              "|SAW(%s,%s)|=%d %s |SAW(%s,%s)|=%d Delta=%+d"
              % (name_of(s), name_of(t0), c0, ">" if c0 > c1 else "<=",
                 name_of(s), name_of(t1), c1, delta))


def classes_by_order():
    """All connected simple graphs up to isomorphism, orders 2..6."""
    out = {}
    for n in (2, 3, 4, 5):
        out[n] = connected_classes_brute(n)
    out[6] = connected_classes_grow(6, out[5])
    return out


def check_class_counts(classes):
    """Check 12: the isomorphism-class census is the known one (1,2,6,21,112)."""
    known = {2: 1, 3: 2, 4: 6, 5: 21, 6: 112}
    got = {n: len(classes[n]) for n in sorted(classes)}
    ok = got == known
    return ck("connected_isomorphism_class_counts", ok,
              "computed=%s expected=%s" % (got, known))


def small_name(n, edges):
    """Name the connected graphs of order <= 4 by order, size and degree sequence."""
    key = (n, len(edges), tuple(sorted(degrees(n, edges))))
    table = {
        (2, 1, (1, 1)): "K_2",
        (3, 2, (1, 1, 2)): "P_3",
        (3, 3, (2, 2, 2)): "K_3",
        (4, 3, (1, 1, 2, 2)): "P_4",
        (4, 3, (1, 1, 1, 3)): "K_{1,3}",
        (4, 4, (2, 2, 2, 2)): "C_4",
        (4, 4, (1, 2, 2, 3)): "paw",
        (4, 5, (2, 2, 3, 3)): "K_4-e",
        (4, 6, (3, 3, 3, 3)): "K_4",
    }
    return table.get(key)


def canon_rooted(n, edges, u, v):
    """Canonical form of a graph together with a distinguished ordered pair."""
    best = None
    for p in permutations(range(n)):
        cur = (tuple(sorted(ek(p[x], p[y]) for (x, y) in edges)), p[u], p[v])
        if best is None or cur < best:
            best = cur
    return best


def check_small_table(classes):
    """Check 13: reproduce the paper's order <= 4 table of max Delta values."""
    got = {}
    for n in (2, 3, 4):
        for g in classes[n]:
            nm = small_name(n, g)
            if nm is None:
                got["UNNAMED(n=%d,m=%d)" % (n, len(g))] = "?"
            else:
                got[nm] = max_delta(n, g)
    ok = got == PAPER_SMALL_TABLE
    shown = ", ".join("%s:%s" % (k, "none" if got[k] is None else got[k])
                      for k in sorted(got))
    return ck("order_at_most_four_delta_table", ok, shown)


def census_violations(n, graphs):
    """Every (graph, ordered adjacent non-cut pair) with Delta > 0."""
    hits = []
    pairs = 0
    worst = None
    for g in graphs:
        for (u, v, d0, d1) in delta_table(n, g):
            pairs += 1
            if worst is None or d0 - d1 > worst:
                worst = d0 - d1
            if d0 > d1:
                hits.append((g, u, v, d0, d1))
    return hits, pairs, worst


def check_minimum_order(classes):
    """Check 14: no violation at order <= 4, and violations exist at order 5."""
    small_hits = 0
    small_pairs = 0
    worst = None
    for n in (1, 2, 3, 4):
        graphs = [()] if n == 1 else classes[n]
        h, p, w = census_violations(n, graphs)
        small_hits += len(h)
        small_pairs += p
        if w is not None and (worst is None or w > worst):
            worst = w
    five, five_pairs, _ = census_violations(5, classes[5])
    ok = (small_hits == 0 and len(five) > 0 and PAPER_MIN_ORDER == 5
          and small_pairs > 0 and worst is not None and worst <= 0)
    return ck("minimum_counterexample_order_is_five", ok,
              "orders<=4: %d pairs tested, %d violations, max Delta=%s; "
              "order 5: %d pairs tested, %d violations"
              % (small_pairs, small_hits, worst, five_pairs, len(five))), five


def check_house_recovered(n, edges, five_hits):
    """Check 15: the blind order-5 census independently rediscovers house/ad."""
    idx = {c: i for i, c in enumerate(HOUSE_VERTEX_NAMES)}
    target = canon_rooted(n, edges, idx["a"], idx["d"])
    matches = [(d0, d1) for (g, u, v, d0, d1) in five_hits
               if canon_rooted(5, g, u, v) == target]
    ok = len(matches) > 0 and all(m == (PAPER_COUNT_D0, PAPER_COUNT_D1)
                                  for m in matches)
    return ck("house_ad_recovered_by_blind_census", ok,
              "census hits isomorphic to (house,ad): %d, counts=%s"
              % (len(matches), sorted(set(matches))))


def check_census_engine(classes, five_hits):
    """Check 16: recount every order 5 and 6 census hit by explicit enumeration."""
    six_hits, six_pairs, six_worst = census_violations(6, classes[6])
    bad = 0
    for (n, hits) in ((5, five_hits), (6, six_hits)):
        for (g, u, v, d0, d1) in hits:
            bn, be = bunkbed(n, g)
            if len(enumerate_paths(bn, be, u, v)) != d0:
                bad += 1
            if len(enumerate_paths(bn, be, u, v + n)) != d1:
                bad += 1
    three_conn = sum(1 for (g, _u, _v, _a, _b) in six_hits
                     if vertex_connectivity_at_least(6, g, 3))
    ok = bad == 0 and len(six_hits) > 0
    note("order 6: %d ordered pairs tested, %d violations, max Delta=%s, "
         "%d of them in 3-connected graphs"
         % (six_pairs, len(six_hits), six_worst, three_conn))
    return ck("census_hits_confirmed_by_enumeration", ok,
              "orders 5-6: %d violating pairs recounted, %d disagreements"
              % (len(five_hits) + len(six_hits), bad))


def build_seven():
    """The paper's 3-connected example on seven vertices."""
    es = set()
    for (x, y) in combinations(PAPER_SEVEN_K4, 2):
        es.add(ek(x, y))
    for (x, y) in combinations(PAPER_SEVEN_K3, 2):
        es.add(ek(x, y))
    for (x, y) in PAPER_SEVEN_LINKS:
        es.add(ek(x, y))
    return 7, sorted(es)


def check_seven_wellformed(n, edges):
    """Check 17: the seven-vertex object is exactly as described."""
    eset = set(edges)
    k4 = all(ek(x, y) in eset for (x, y) in combinations(PAPER_SEVEN_K4, 2))
    k3 = all(ek(x, y) in eset for (x, y) in combinations(PAPER_SEVEN_K3, 2))
    links = [ek(x, y) for (x, y) in PAPER_SEVEN_LINKS]
    touched = [v for e in links for v in e]
    disjoint = len(set(touched)) == 6
    inside = set()
    for (x, y) in combinations(PAPER_SEVEN_K4, 2):
        inside.add(ek(x, y))
    for (x, y) in combinations(PAPER_SEVEN_K3, 2):
        inside.add(ek(x, y))
    extra = sorted(eset - inside - set(links))
    deg = degrees(n, edges)
    ok = (n == 7 and is_simple(n, edges) and len(edges) == 12 and k4 and k3
          and all(e in eset for e in links) and disjoint and not extra
          and sorted(deg) == [3, 3, 3, 3, 4, 4, 4] and min(deg) == 3)
    note("seven-vertex E = " + ",".join("%d%d" % e for e in edges))
    return ck("seven_vertex_graph_wellformed", ok,
              "n=%d m=%d K4_present=%s K3_present=%s links_disjoint=%s "
              "unexpected_edges=%s degrees=%s"
              % (n, len(edges), k4, k3, disjoint, extra, sorted(deg)))


def check_seven_connectivity(n, edges):
    """Check 18: kappa = 3 exactly for the seven-vertex example."""
    at3 = vertex_connectivity_at_least(n, edges, 3)
    at4 = vertex_connectivity_at_least(n, edges, 4)
    witness = None
    for S in combinations(range(n), 3):
        drop = 0
        for v in S:
            drop |= 1 << v
        if not is_connected(n, edges, drop):
            witness = S
            break
    kappa = 3 if (at3 and not at4) else None
    ok = (at3 and not at4 and witness is not None
          and kappa == PAPER_SEVEN_CONNECTIVITY)
    return ck("seven_vertex_is_3_connected", ok,
              "kappa>=3:%s kappa>=4:%s separating_triple=%s kappa=%s"
              % (at3, at4, witness, kappa))


def check_seven_counts(n, edges):
    """Checks 19-20: the seven-vertex counts and the violated inequality."""
    u, v = PAPER_SEVEN_LINKS[0]
    cut = is_cut_edge(n, edges, ek(u, v))
    bn, be = bunkbed(n, edges)
    out = paths_by_length(bn, be, u)
    c0, c1 = sum(out[v]), sum(out[v + n])
    e0 = len(enumerate_paths(bn, be, u, v))
    e1 = len(enumerate_paths(bn, be, u, v + n))
    ck("seven_vertex_counts_match_paper",
       c0 == PAPER_SEVEN_COUNT_0 and c1 == PAPER_SEVEN_COUNT_1
       and e0 == c0 and e1 == c1 and not cut and bn == 14 and len(be) == 31,
       "recurrence=(%d,%d) enumeration=(%d,%d) paper=(%d,%d) "
       "edge_%d%d_is_cut=%s bunkbed n=%d m=%d"
       % (c0, c1, e0, e1, PAPER_SEVEN_COUNT_0, PAPER_SEVEN_COUNT_1,
          u, v, cut, bn, len(be)))
    ck("seven_vertex_inequality_is_violated", c0 > c1,
       "|SAW(0_0,5_0)|=%d %s |SAW(0_0,5_1)|=%d Delta=%+d"
       % (c0, ">" if c0 > c1 else "<=", c1, c0 - c1))


def main():
    print("Verification of a five-vertex refutation of a bunkbed "
          "self-avoiding-walk question")
    n, edges = decode_house()
    bn, be = bunkbed(n, edges)
    a0 = 0
    d0, d1 = 3, 3 + n

    check_house_wellformed(n, edges)
    check_house_two_connected(n, edges)
    check_ad_non_cut(n, edges)
    check_bunkbed_wellformed(n, edges, bn, be)
    w0 = check_counts(bn, be, a0, d0, PAPER_COUNT_D0, "a0_d0")
    w1 = check_counts(bn, be, a0, d1, PAPER_COUNT_D1, "a0_d1")
    check_walks_valid(bn, be, [(a0, d0, w0), (a0, d1, w1)])
    check_recurrence_agrees(bn, be, a0, [d0, d1],
                            [(a0, d0, w0), (a0, d1, w1)])
    check_profile(bn, be, a0, d0, w0, PAPER_PROFILE_D0, "a0_d0")
    check_profile(bn, be, a0, d1, w1, PAPER_PROFILE_D1, "a0_d1")
    check_refutation(bn, be, a0, d0, d1)

    classes = classes_by_order()
    check_class_counts(classes)
    check_small_table(classes)
    _res, five_hits = check_minimum_order(classes)
    check_house_recovered(n, edges, five_hits)
    check_census_engine(classes, five_hits)

    sn, sedges = build_seven()
    check_seven_wellformed(sn, sedges)
    check_seven_connectivity(sn, sedges)
    check_seven_counts(sn, sedges)

    note("not re-run here: no exhaustive census beyond order 6 was attempted, "
         "so the search for further counterexamples stops at six vertices; "
         "the paper makes no claim past order five, whose minimality is "
         "settled by the complete order <= 4 census above.")

    failed = [nm for (nm, ok) in CHECKS if not ok]
    total = len(CHECKS)
    if failed:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


if __name__ == "__main__":
    sys.exit(main())

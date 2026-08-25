#!/usr/bin/env python3
"""
Verification program for "Order-Minimal Pendant-Free Counterexamples to
Conjecture 2 of Nurdin et al."

Standard library only, exact integer arithmetic, no data files, no network.

VALUES TAKEN FROM THE PAPER (inputs; transcribed, never re-derived here)
-----------------------------------------------------------------------
  * E(G1) = {01,12,23,04,34,15,35,06,16} on the vertex set {0,...,6};
    E(G2) = {01,02,23,14,34,05,45,06,36}.
  * The graph6 strings "FhdU?" (G1), "FpUK_" (G2), "FreRW" (atlas K_{3,4}).
  * The claimed degree sequence (2,2,2,2,3,3,4) of G1 and of G2.
  * The claimed values B(G) = 2 and tvs(G) = 3 for G1, G2 and K_{3,4},
    together with the individual ceiling terms displayed for B.
  * The three total 3-labeling certificates of the table: vertex label
    strings, edge label strings read in the paper's stated edge orders, and
    the claimed vertex weight lists.
  * The census counts 1, 3, 11, 61, 507 of connected graphs with minimum
    degree at least two and order 3,4,5,6,7; the total 583; and the claim
    that 580 of them admit a vertex-irregular total B(G)-labeling.
  * Theorem 2's hypotheses (t >= 1, d >= max(2t, t+2), H a (d-t-2)-regular
    graph on d vertices, G = H join complement(K_{2t+2})) and its claim
    B(G) = 2 < tvs(G); Corollary 3's witnesses t = 1, d = n-4,
    H = complement(C_d) for n >= 8, with complement(G) = C_d + K_4.

WHAT IS DERIVED HERE (computed; nothing in this list is assumed)
---------------------------------------------------------------
  * The "FhdU?" and "FpUK_" strings are decoded and compared label for label
    with the paper's edge sets for G1 and G2.  "FreRW" decodes to a DIFFERENT
    labelling of K_{3,4} from the one the certificate's edge order uses, so it
    is verified to be complete bipartite with parts of sizes 3 and 4 and is
    matched to the certificate's copy up to isomorphism, by equality of
    canonical keys; both bipartitions are derived and printed, and the failure
    of the label-for-label match is itself asserted as a check.
  * Order, connectivity, minimum degree and degree sequences of all three
    graphs (the hypotheses of the theorem being settled).
  * B(G) recomputed from the Nurdin et al. formula, every ceiling term.
  * tvs(G) > 2 for each of the three graphs by COMPLETE exhaustion of all
    total 2-labelings (exact search over every edge-label vector combined
    with an exact system-of-distinct-representatives test on vertex labels).
  * tvs(G) <= 3 by recomputing every weight of the table's certificates.
  * Non-isomorphism of G1 and G2, and that neither is complete multipartite,
    while K_{3,4} is.
  * The census: every graph of order at most 7 is regenerated from scratch
    up to isomorphism, filtered to connected with delta >= 2, counted, and
    each graph is either given an explicit machine-verified B(G)-labeling or
    proved to have none by exhaustion.
  * Theorem 2 and Corollary 3: the constructed graphs are built, their
    degrees and B(G) are computed, and tvs(G) > 2 is obtained from a
    computed endpoint-weight obstruction whose implementation is
    cross-validated against complete exhaustion on small members.  These two
    statements are infinite families, so what is checked here is a finite
    window of each: for Theorem 2, every (d-t-2)-regular core up to
    isomorphism when d <= 7 (taken from the regenerated census) and a single
    circulant core when 8 <= d <= 12; for Corollary 3, orders 8..16.  The
    closing "NOT RE-RUN HERE" line states the shortfall exactly.
"""

import itertools
import sys

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)


# ----------------------------------------------------------------- inputs

PAPER_G1_EDGES = [(0, 1), (1, 2), (2, 3), (0, 4), (3, 4), (1, 5), (3, 5),
                  (0, 6), (1, 6)]
PAPER_G2_EDGES = [(0, 1), (0, 2), (2, 3), (1, 4), (3, 4), (0, 5), (4, 5),
                  (0, 6), (3, 6)]
PAPER_G6 = {"G1": "FhdU?", "G2": "FpUK_", "K34": "FreRW"}
PAPER_DEGSEQ_G12 = (2, 2, 2, 2, 3, 3, 4)
PAPER_B_VALUE = 2
PAPER_TVS_VALUE = 3
PAPER_CENSUS = {3: 1, 4: 3, 5: 11, 6: 61, 7: 507}
PAPER_CENSUS_TOTAL = 583
PAPER_LABELABLE = 580
K34_EDGE_ORDER = [(0, 3), (0, 4), (0, 5), (0, 6), (1, 3), (1, 4), (1, 5),
                  (1, 6), (2, 3), (2, 4), (2, 5), (2, 6)]
PAPER_CERTS = {
    "G1": ("3313233", "111111123", [7, 9, 3, 6, 4, 5, 8]),
    "G2": ("3123333", "111111132", [9, 3, 4, 7, 6, 5, 8]),
    "K34": ("2221222", "111211231233", [7, 9, 11, 4, 6, 8, 10]),
}


# -------------------------------------------------------- graph utilities


def graph6_decode(s):
    """Decode a graph6 string into (order, sorted edge list)."""
    n = ord(s[0]) - 63
    bits = []
    for c in s[1:]:
        v = ord(c) - 63
        if v < 0 or v > 63:
            raise ValueError("bad graph6 character")
        for i in range(5, -1, -1):
            bits.append((v >> i) & 1)
    need = n * (n - 1) // 2
    if len(bits) < need:
        raise ValueError("graph6 string too short")
    edges = []
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                edges.append((i, j))
            idx += 1
    if any(bits[need:]):
        raise ValueError("graph6 padding bits are not zero")
    return n, sorted(edges)


def normalise(edges):
    return sorted(tuple(sorted(e)) for e in edges)


def adjacency(n, edges):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        if u == v:
            raise ValueError("loop")
        adj[u].add(v)
        adj[v].add(u)
    return adj


def degrees(n, edges):
    adj = adjacency(n, edges)
    return [len(adj[v]) for v in range(n)]


def is_connected(n, edges):
    if n == 0:
        return False
    adj = adjacency(n, edges)
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for u in adj[v]:
            if u not in seen:
                seen.add(u)
                stack.append(u)
    return len(seen) == n


def nurdin_bound_terms(n, edges):
    """Return (B, [(i, ceiling term) ...]) for the Nurdin et al. bound."""
    deg = degrees(n, edges)
    dmin, dmax = min(deg), max(deg)
    cnt = {}
    for d in deg:
        cnt[d] = cnt.get(d, 0) + 1
    terms = []
    running = 0
    for i in range(dmin, dmax + 1):
        running += cnt.get(i, 0)
        num = dmin + running
        terms.append((i, -((-num) // (i + 1))))
    return max(t[1] for t in terms), terms


# ------------------------------------------------- total labelings, exactly


def sdr(items, k):
    """System of distinct representatives for equal-length integer intervals.

    items is a list of (lo, vertex); vertex may take any weight in
    [lo, lo+k-1].  Returns a dict vertex -> weight, or None if impossible.
    Greedy by left endpoint is exact for intervals of equal length.
    """
    used = set()
    out = {}
    for lo, v in sorted(items):
        w = lo
        while w in used:
            w += 1
        if w > lo + k - 1:
            return None
        used.add(w)
        out[v] = w
    return out


def total_labeling_search(n, edges, k):
    """Exhaustively search for a vertex-irregular total k-labeling.

    Returns (labeling, nodes) where labeling is None iff none exists; the
    search is complete, so None is a proof of infeasibility.
    """
    edges = sorted(normalise(edges), key=lambda e: (e[1], e[0]))
    m = len(edges)
    deg = degrees(n, edges)
    last = {}
    for idx, (u, v) in enumerate(edges):
        last[u] = idx
        last[v] = idx
    fin_after = [[] for _ in range(m)]
    free = [v for v in range(n) if v not in last]
    for v, idx in last.items():
        fin_after[idx].append(v)
    b = [0] * n
    offs = [0] * m
    nodes = [0]
    fin_upto = [list(free)]

    def rec(idx):
        nodes[0] += 1
        if idx == m:
            return sdr([(1 + deg[v] + b[v], v) for v in range(n)], k)
        u, v = edges[idx]
        for o in range(k):
            offs[idx] = o
            b[u] += o
            b[v] += o
            ok = True
            if fin_after[idx]:
                fin_upto.append(fin_upto[-1] + fin_after[idx])
                ok = sdr([(1 + deg[w] + b[w], w) for w in fin_upto[-1]],
                         k) is not None
            got = rec(idx + 1) if ok else None
            if fin_after[idx]:
                fin_upto.pop()
            if got:
                return got
            b[u] -= o
            b[v] -= o
        return None

    wts = rec(0)
    if wts is None:
        return None, nodes[0]
    elab = dict((edges[i], offs[i] + 1) for i in range(m))
    vlab = dict((v, wts[v] - deg[v] - sum(elab[e] - 1 for e in edges
                                         if v in e)) for v in range(n))
    return (vlab, elab), nodes[0]


def verify_labeling(n, edges, vlab, elab, k):
    """Recompute weights; return (weights, ok) with ok verified from scratch."""
    edges = normalise(edges)
    adj_e = [[] for _ in range(n)]
    for e in edges:
        adj_e[e[0]].append(e)
        adj_e[e[1]].append(e)
    ok = all(1 <= vlab[v] <= k for v in range(n))
    ok = ok and all(1 <= elab[e] <= k for e in edges)
    ok = ok and len(elab) == len(edges) and len(vlab) == n
    wts = [vlab[v] + sum(elab[e] for e in adj_e[v]) for v in range(n)]
    ok = ok and len(set(wts)) == n
    return wts, ok


# ------------------------------------ isomorphism and the order <= 7 census


def _refine(n, adj):
    ranks = sorted(set(len(adj[v]) for v in range(n)))
    colors = [ranks.index(len(adj[v])) for v in range(n)]
    for _ in range(n):
        sig = [(colors[v], tuple(sorted(colors[u] for u in adj[v])))
               for v in range(n)]
        uniq = sorted(set(sig))
        new = [uniq.index(sig[v]) for v in range(n)]
        if new == colors:
            break
        colors = new
    return colors


def canonical_key(n, edges):
    """Isomorphism-invariant canonical code: max upper-triangle bitmask over
    all relabelings that respect the (invariant) colour refinement order."""
    adj = adjacency(n, edges)
    colors = _refine(n, adj)
    blocks = {}
    for v in range(n):
        blocks.setdefault(colors[v], []).append(v)
    blocks = [blocks[c] for c in sorted(blocks)]
    best = -1
    perm = [0] * n
    for combo in itertools.product(*[itertools.permutations(b)
                                     for b in blocks]):
        pos = 0
        for tup in combo:
            for v in tup:
                perm[v] = pos
                pos += 1
        code = 0
        for u, v in edges:
            a, c = perm[u], perm[v]
            if a > c:
                a, c = c, a
            code |= 1 << (c * (c - 1) // 2 + a)
        if code > best:
            best = code
    return (n, best)


def key_to_edges(key):
    n, code = key
    edges = []
    for j in range(1, n):
        for i in range(j):
            if code >> (j * (j - 1) // 2 + i) & 1:
                edges.append((i, j))
    return n, edges


def census(nmax):
    """All simple graphs of order 1..nmax up to isomorphism, generated by
    one-vertex extension and reduced by canonical_key."""
    levels = {1: {(1, 0)}}
    for n in range(2, nmax + 1):
        got = set()
        for key in levels[n - 1]:
            _, edges = key_to_edges(key)
            for r in range(n):
                for nbrs in itertools.combinations(range(n - 1), r):
                    new = edges + [(i, n - 1) for i in nbrs]
                    got.add(canonical_key(n, new))
        levels[n] = got
    return levels


# ------------------------------------- constructions and the join obstruction


def complement_edges(n, edges):
    es = set(normalise(edges))
    return [(i, j) for j in range(n) for i in range(j) if (i, j) not in es]


def join(n1, e1, n2, e2):
    """Vertices 0..n1-1 carry the first graph, n1..n1+n2-1 the second."""
    e = list(normalise(e1)) + [(u + n1, v + n1) for u, v in normalise(e2)]
    e += [(i, n1 + j) for i in range(n1) for j in range(n2)]
    return n1 + n2, normalise(e)


def cycle_edges(d):
    return normalise([(i, (i + 1) % d) for i in range(d)])


def circulant_regular(d, r):
    """An r-regular circulant graph on d vertices, or None if r*d is odd."""
    if r < 0 or r > d - 1 or (r * d) % 2:
        return None
    e = set()
    for i in range(d):
        for s in range(1, r // 2 + 1):
            e.add(tuple(sorted((i, (i + s) % d))))
        if r % 2:
            e.add(tuple(sorted((i, (i + d // 2) % d))))
    return sorted(e)


def is_complete_multipartite(n, edges):
    """True iff non-adjacency is transitive, i.e. the complement is a
    disjoint union of cliques."""
    comp = complement_edges(n, edges)
    adj = adjacency(n, comp)
    seen = set()
    for s in range(n):
        if s in seen:
            continue
        comp_vs = {s}
        stack = [s]
        while stack:
            v = stack.pop()
            for u in adj[v]:
                if u not in comp_vs:
                    comp_vs.add(u)
                    stack.append(u)
        seen |= comp_vs
        for u, v in itertools.combinations(sorted(comp_vs), 2):
            if v not in adj[u]:
                return False
    return True


def endpoint_obstruction(n, edges):
    """Computed sufficient condition for tvs(G) > 2.

    At k = 2 the weight of v lies in [1+deg v, 2+2 deg v].  If the union of
    these ranges holds exactly n integers, the n distinct weights exhaust it;
    the minimum forces a minimum-degree vertex and its edges to be labelled
    1, the maximum forces a maximum-degree vertex and its edges to be
    labelled 2, so any two such vertices must be non-adjacent and distinct.
    Returns True only when that fails for every admissible pair.
    """
    deg = degrees(n, edges)
    allowed = set()
    for v in range(n):
        allowed |= set(range(1 + deg[v], 3 + 2 * deg[v]))
    if len(allowed) != n:
        return False
    lo, hi = min(allowed), max(allowed)
    mins = [v for v in range(n) if 1 + deg[v] == lo]
    maxs = [v for v in range(n) if 2 + 2 * deg[v] == hi]
    if not mins or not maxs:
        return False
    adj = adjacency(n, edges)
    return all(u != w and w in adj[u] for u in mins for w in maxs)


# ------------------------------------------------------------ the checks


def safe_decode(s):
    try:
        return graph6_decode(s)
    except (ValueError, IndexError):
        return (0, [])


def check_exhibited_objects():
    """Decode, count and print back the three exhibited graphs."""
    n1, e1 = safe_decode(PAPER_G6["G1"])
    n2, e2 = safe_decode(PAPER_G6["G2"])
    nk, ek = safe_decode(PAPER_G6["K34"])
    ck("graph6_FhdU_decodes_to_paper_E_G1",
       (n1, e1) == (7, normalise(PAPER_G1_EDGES)),
       "n=%d, %d edges: %s" % (n1, len(e1),
                               "".join("%d%d " % e for e in e1).strip()))
    ck("graph6_FpUK_decodes_to_paper_E_G2",
       (n2, e2) == (7, normalise(PAPER_G2_EDGES)),
       "n=%d, %d edges: %s" % (n2, len(e2),
                               "".join("%d%d " % e for e in e2).strip()))
    k34 = (7, normalise(K34_EDGE_ORDER))

    def parts_by_degree(n, e):
        """Group the vertices by degree; for a complete bipartite graph with
        parts of distinct sizes this recovers the bipartition."""
        blocks = {}
        deg = degrees(n, e)
        for v in range(n):
            blocks.setdefault(deg[v], []).append(v)
        return sorted((sorted(p) for p in blocks.values()), key=len)

    def show(blocks):
        """Printable bipartition; never indexes, so a failed decode still
        reports through the normal PASS/FAIL line instead of raising."""
        return " | ".join(str(b) for b in blocks) if blocks else "(none)"

    parts = {}
    for v in range(nk):
        parts.setdefault(degrees(nk, ek)[v], []).append(v)
    sizes = sorted(len(p) for p in parts.values())
    bip_ok = (sizes == [3, 4] and len(ek) == 12
              and all((min(u, v), max(u, v)) in set(ek)
                      for p in parts.values() for q in parts.values()
                      if p is not q for u in p for v in q)
              and not any((min(u, v), max(u, v)) in set(ek)
                          for p in parts.values()
                          for u, v in itertools.combinations(p, 2)))
    dec_parts = parts_by_degree(nk, ek)
    cert_parts = parts_by_degree(*k34)
    ck("graph6_FreRW_is_complete_bipartite_3_plus_4",
       bip_ok and canonical_key(nk, ek) == canonical_key(*k34),
       "order %d, %d edges, part sizes %s; decoded bipartition %s, "
       "compared with the certificate's copy up to isomorphism (equal "
       "canonical keys), not label for label"
       % (nk, len(ek), sizes, show(dec_parts)))
    # The decoded labelling is NOT the certificate's labelling.  A referee who
    # attempts a label-for-label match will get a mismatch, so the mismatch is
    # recorded here as a derived fact rather than left to be discovered.
    ck("graph6_FreRW_is_a_relabelling_of_the_certificate_copy_not_label_equal",
       set(ek) != set(normalise(K34_EDGE_ORDER))
       and canonical_key(nk, ek) == canonical_key(*k34),
       "decoded edges %s with bipartition %s; the paper's certificate edge "
       "order gives bipartition %s, so the two edge sets differ while the "
       "canonical keys agree"
       % ("".join("%d%d " % e for e in ek).strip(),
          show(dec_parts), show(cert_parts)))
    return {"G1": (7, normalise(PAPER_G1_EDGES)),
            "G2": (7, normalise(PAPER_G2_EDGES)), "K34": k34}


def check_hypotheses(graphs):
    """Every hypothesis of the theorem: order <= 7, connected, delta >= 2."""
    rows = []
    ok = True
    for name in ("G1", "G2", "K34"):
        n, e = graphs[name]
        deg = degrees(n, e)
        good = (n == 7 and is_connected(n, e) and min(deg) >= 2)
        ok = ok and good
        rows.append("%s: n=%d conn=%s delta=%d" % (name, n, is_connected(n, e),
                                                   min(deg)))
    ck("hypotheses_order_at_most_7_connected_min_degree_2", ok, "; ".join(rows))
    d1 = tuple(sorted(degrees(*graphs["G1"])))
    d2 = tuple(sorted(degrees(*graphs["G2"])))
    dk = tuple(sorted(degrees(*graphs["K34"])))
    ck("degree_sequences_match_paper",
       d1 == PAPER_DEGSEQ_G12 and d2 == PAPER_DEGSEQ_G12
       and dk == (3, 3, 3, 3, 4, 4, 4),
       "G1=%s G2=%s K34=%s" % (d1, d2, dk))


def check_bound(graphs):
    """B(G) = 2 and each displayed ceiling term."""
    want_terms = {"G1": [(2, 2), (3, 2), (4, 2)], "G2": [(2, 2), (3, 2),
                                                         (4, 2)],
                  "K34": [(3, 2), (4, 2)]}
    ok = True
    rows = []
    for name in ("G1", "G2", "K34"):
        B, terms = nurdin_bound_terms(*graphs[name])
        ok = ok and B == PAPER_B_VALUE and terms == want_terms[name]
        rows.append("B(%s)=%d terms=%s" % (name, B,
                                           [t[1] for t in terms]))
    ck("nurdin_bound_B_equals_2_with_the_displayed_terms", ok, "; ".join(rows))


def check_no_two_labeling(graphs):
    """The load-bearing refutation: no total 2-labeling exists.  Complete
    exhaustion of the edge-label vectors, exact matching on vertex labels."""
    out = {}
    for name in ("G1", "G2", "K34"):
        n, e = graphs[name]
        lab, nodes = total_labeling_search(n, e, 2)
        out[name] = lab is None
        ck("no_vertex_irregular_total_2_labeling_of_" + name, lab is None,
           "search complete, %d nodes, 2^%d edge-label vectors covered"
           % (nodes, len(e)))
    return out


def check_naive_enumeration(graphs):
    """Independent confirmation of the decisive claim with no pruning and no
    matching theory: enumerate literally every total 2-labeling."""
    rows = []
    ok = True
    for name in ("G1", "G2", "K34"):
        n, e = graphs[name]
        adj_e = [[] for _ in range(n)]
        for idx, (u, v) in enumerate(e):
            adj_e[u].append(idx)
            adj_e[v].append(idx)
        found = None
        total = 0
        for evec in itertools.product((1, 2), repeat=len(e)):
            base = [sum(evec[i] for i in adj_e[v]) for v in range(n)]
            for vvec in itertools.product((1, 2), repeat=n):
                total += 1
                w = [base[v] + vvec[v] for v in range(n)]
                if len(set(w)) == n:
                    found = (evec, vvec)
                    break
            if found:
                break
        ok = ok and found is None
        rows.append("%s: %d labelings examined, irregular ones found: %d"
                    % (name, total, 0 if found is None else 1))
    ck("naive_full_enumeration_finds_no_total_2_labeling", ok, "; ".join(rows))


def check_sdr_engine():
    """The pruning engine's exactness: the greedy system-of-distinct-
    representatives routine agrees with brute force on all small instances."""
    bad = []
    tried = 0
    for k in (1, 2, 3):
        for size in (1, 2, 3, 4):
            for los in itertools.product(range(1, 6), repeat=size):
                items = [(los[i], i) for i in range(size)]
                got = sdr(items, k)
                brute = None
                for pick in itertools.product(*[range(lo, lo + k)
                                                for lo in los]):
                    if len(set(pick)) == size:
                        brute = pick
                        break
                tried += 1
                if (got is None) != (brute is None):
                    bad.append((k, los))
                elif got is not None:
                    if any(not (los[i] <= got[i] <= los[i] + k - 1)
                           for i in range(size)) or \
                            len(set(got.values())) != size:
                        bad.append((k, los))
    ck("greedy_sdr_engine_agrees_with_brute_force", not bad,
       "%d instances, %d disagreements" % (tried, len(bad)))


def check_certificates(graphs):
    """The table's total 3-labelings, recomputed weight by weight."""
    out = {}
    orders = {"G1": PAPER_G1_EDGES, "G2": PAPER_G2_EDGES,
              "K34": K34_EDGE_ORDER}
    for name in ("G1", "G2", "K34"):
        n, e = graphs[name]
        vstr, estr, want = PAPER_CERTS[name]
        order = [tuple(sorted(x)) for x in orders[name]]
        good = len(vstr) == n and len(estr) == len(order)
        vlab = dict((v, int(vstr[v])) for v in range(n))
        elab = dict((order[i], int(estr[i])) for i in range(len(order)))
        wts, ok = verify_labeling(n, e, vlab, elab, PAPER_TVS_VALUE)
        good = good and ok and wts == want and set(order) == set(e)
        out[name] = good
        ck("certificate_is_a_valid_total_3_labeling_of_" + name, good,
           "weights %s, all distinct=%s, max label %d"
           % (wts, len(set(wts)) == n,
              max(list(vlab.values()) + list(elab.values()))))
    return out


def check_tvs_value(no_two, certs):
    ok = all(no_two[k] for k in no_two) and all(certs[k] for k in certs)
    ck("tvs_equals_3_and_exceeds_B_equals_2_for_all_three", ok,
       "tvs>%d from exhaustion, tvs<=%d from certificates, B=%d computed"
       % (PAPER_B_VALUE, PAPER_TVS_VALUE, PAPER_B_VALUE))


def check_structure(graphs):
    """Triangle content, non-isomorphism, complete multipartiteness."""
    def triangles(n, e):
        adj = adjacency(n, e)
        return sum(1 for a, b, c in itertools.combinations(range(n), 3)
                   if b in adj[a] and c in adj[b] and c in adj[a])
    t1 = triangles(*graphs["G1"])
    t2 = triangles(*graphs["G2"])
    k1 = canonical_key(*graphs["G1"])
    k2 = canonical_key(*graphs["G2"])
    ck("G1_has_a_triangle_G2_is_triangle_free_and_they_differ",
       t1 >= 1 and t2 == 0 and k1 != k2,
       "triangles: G1=%d G2=%d; canonical codes differ=%s" % (t1, t2,
                                                              k1 != k2))
    rows = []
    ok = True
    for name in ("G1", "G2"):
        n, e = graphs[name]
        adj = adjacency(n, e)
        witness = (3 not in adj[0] and 1 not in adj[3] and 1 in adj[0])
        cm = is_complete_multipartite(n, e)
        ok = ok and witness and not cm
        rows.append("%s: 0~1 but 0/~3/~1 witness=%s, complete multipartite=%s"
                    % (name, witness, cm))
    ok = ok and is_complete_multipartite(*graphs["K34"])
    ck("neither_G1_nor_G2_is_complete_multipartite_while_K34_is", ok,
       "; ".join(rows) + "; K34=%s" % is_complete_multipartite(*graphs["K34"]))


def build_census():
    """Regenerate every simple graph of order <= 7 up to isomorphism and keep
    the connected ones with minimum degree at least two."""
    levels = census(7)
    counts = dict((n, len(levels[n])) for n in levels)
    ck("graph_generator_reproduces_the_known_isomorphism_class_counts",
       [counts[n] for n in range(1, 8)] == [1, 2, 4, 11, 34, 156, 1044],
       "orders 1..7: %s (expected 1,2,4,11,34,156,1044)"
       % [counts[n] for n in range(1, 8)])
    pool = {}
    for n in range(3, 8):
        keep = []
        for key in sorted(levels[n]):
            nn, e = key_to_edges(key)
            if min(degrees(nn, e)) >= 2 and is_connected(nn, e):
                keep.append((nn, e, key))
        pool[n] = keep
    got = dict((n, len(pool[n])) for n in pool)
    ck("census_of_connected_min_degree_2_graphs_matches_the_table",
       got == PAPER_CENSUS and sum(got.values()) == PAPER_CENSUS_TOTAL,
       "orders 3..7: %s, total %d" % ([got[n] for n in range(3, 8)],
                                      sum(got.values())))
    return levels, pool


def check_census_search(pool, graphs):
    """Every census graph: an explicit verified B(G)-labeling, or none at all."""
    labelled = 0
    exceptions = []
    for n in sorted(pool):
        for nn, e, key in pool[n]:
            B, _ = nurdin_bound_terms(nn, e)
            lab, _ = total_labeling_search(nn, e, B)
            if lab is None:
                exceptions.append(key)
                continue
            wts, ok = verify_labeling(nn, e, lab[0], lab[1], B)
            if ok:
                labelled += 1
            else:
                exceptions.append(("bad-labeling", key))
    ck("census_580_graphs_admit_a_verified_B_labeling",
       labelled == PAPER_LABELABLE,
       "%d of %d graphs certified at k=B(G)" % (labelled,
                                               PAPER_CENSUS_TOTAL))
    want = sorted(canonical_key(*graphs[k]) for k in ("G1", "G2", "K34"))
    ck("the_only_census_exceptions_are_the_three_exhibited_graphs",
       sorted(exceptions) == want,
       "%d exception(s); canonical codes match the three graphs=%s"
       % (len(exceptions), sorted(exceptions) == want))
    return exceptions


def check_lower_bound_tight(pool):
    """No census graph admits a total (B(G)-1)-labeling: the Nurdin et al.
    lower bound is confirmed graph by graph, so tvs = B for the 580."""
    viol = []
    for n in sorted(pool):
        for nn, e, key in pool[n]:
            B, _ = nurdin_bound_terms(nn, e)
            lab, _ = total_labeling_search(nn, e, B - 1)
            if lab is not None:
                viol.append(key)
    ck("no_census_graph_has_a_total_B_minus_1_labeling", not viol,
       "0 violations over %d graphs" % PAPER_CENSUS_TOTAL if not viol
       else "%d violations" % len(viol))


def check_obstruction_lemma(pool):
    """Validate the endpoint-weight obstruction against the census: whenever
    it fires, complete exhaustion must agree that k=2 is infeasible."""
    fired = 0
    wrong = []
    for n in sorted(pool):
        for nn, e, key in pool[n]:
            if not endpoint_obstruction(nn, e):
                continue
            fired += 1
            lab, _ = total_labeling_search(nn, e, 2)
            if lab is not None:
                wrong.append(key)
    ck("endpoint_obstruction_never_contradicts_exhaustive_search",
       fired > 0 and not wrong,
       "fired on %d census graphs, %d disagreements" % (fired, len(wrong)))


def check_family(tmax=4, dmax=12):
    """Theorem 2 on every admissible (t,d) in a finite window: the degrees,
    the order, B(G)=2, and tvs(G)>2 from the computed obstruction.

    One core H per pair, a circulant; check_family_all_cores widens this to
    every (d-t-2)-regular core up to isomorphism for d <= 7.
    """
    done = []
    bad = []
    skipped = []
    for t in range(1, tmax + 1):
        for d in range(max(2 * t, t + 2), dmax + 1):
            H = circulant_regular(d, d - t - 2)
            if H is None:
                # No (d-t-2)-regular graph on d vertices exists at all here,
                # so the theorem's hypothesis is empty; assert that reason
                # rather than skipping silently.
                if (d - t - 2) * d % 2 == 0:
                    bad.append(("core-missing-but-parity-allows-it", t, d))
                else:
                    skipped.append((t, d))
                continue
            if set(degrees(d, H)) != {d - t - 2}:
                bad.append(("core-not-regular", t, d))
                continue
            n, e = join(d, H, 2 * t + 2, [])
            deg = degrees(n, e)
            B, _ = nurdin_bound_terms(n, e)
            ok = (n == d + 2 * t + 2
                  and sorted(set(deg)) == sorted({d, d + t})
                  and [deg[v] for v in range(d)] == [d + t] * d
                  and [deg[v] for v in range(d, n)] == [d] * (2 * t + 2)
                  and is_connected(n, e) and min(deg) >= 2
                  and B == 2 and endpoint_obstruction(n, e))
            if ok:
                done.append((t, d))
            else:
                bad.append((t, d, B, sorted(set(deg))))
    ck("theorem_2_family_B_equals_2_and_tvs_gt_2_on_a_finite_window",
       len(done) >= 20 and not bad,
       "%d admissible (t,d) pairs verified with t<=%d, d<=%d, one circulant "
       "core each; %d further pairs carry no (d-t-2)-regular graph at all "
       "((d-t-2)d odd) and are vacuous; %d failures"
       % (len(done), tmax, dmax, len(skipped), len(bad)))
    return done


def check_family_all_cores(levels, tmax=4, dmax=7):
    """Theorem 2 for EVERY admissible core, not just a circulant one.

    The theorem allows H to be any (d-t-2)-regular graph on d vertices, while
    check_family builds a single circulant per (t,d).  Since the census has
    already regenerated every graph of order at most 7 up to isomorphism, all
    admissible cores with d <= 7 are available here for nothing, and each one
    is put through the same conclusions.  For each pair it is also derived
    that cores exist exactly when (d-t-2)d is even, which is what licenses
    check_family's skipping of the remaining pairs.
    """
    bad = []
    pairs = 0
    cores_tested = 0
    for t in range(1, tmax + 1):
        for d in range(max(2 * t, t + 2), dmax + 1):
            r = d - t - 2
            pairs += 1
            if d not in levels:
                # The census does not reach this order, so the enumeration
                # would be empty: record it as a failure, never as silence.
                bad.append(("census-missing-order", t, d))
                continue
            cores = []
            for key in sorted(levels[d]):
                nn, e = key_to_edges(key)
                if set(degrees(nn, e)) == {r}:
                    cores.append(e)
            if (len(cores) > 0) != (r * d % 2 == 0):
                bad.append(("core-existence-vs-parity", t, d, len(cores)))
            for e in cores:
                n, ge = join(d, e, 2 * t + 2, [])
                deg = degrees(n, ge)
                B, _ = nurdin_bound_terms(n, ge)
                cores_tested += 1
                ok = (n == d + 2 * t + 2
                      and [deg[v] for v in range(d)] == [d + t] * d
                      and [deg[v] for v in range(d, n)] == [d] * (2 * t + 2)
                      and is_connected(n, ge) and min(deg) >= 2
                      and B == 2 and endpoint_obstruction(n, ge))
                if not ok:
                    bad.append((t, d, canonical_key(d, e), B))
    ck("theorem_2_holds_for_every_regular_core_up_to_isomorphism_when_d_le_7",
       cores_tested > 0 and not bad,
       "%d (t,d) pairs with t<=%d, d<=%d; %d distinct cores enumerated up to "
       "isomorphism from the regenerated census, each giving B=2 and a firing "
       "obstruction; %d failures"
       % (pairs, tmax, dmax, cores_tested, len(bad)))


def check_family_bruteforce():
    """Cross-validate the obstruction against complete exhaustion on the
    smallest members of the family (t,d) = (1,4), (2,4), (1,5)."""
    rows = []
    ok = True
    for t, d in ((1, 4), (2, 4), (1, 5)):
        H = circulant_regular(d, d - t - 2)
        n, e = join(d, H, 2 * t + 2, [])
        lab, nodes = total_labeling_search(n, e, 2)
        lem = endpoint_obstruction(n, e)
        ok = ok and lab is None and lem
        rows.append("(t=%d,d=%d) n=%d |E|=%d exhaustive=%s (%d nodes) "
                    "obstruction=%s"
                    % (t, d, n, len(e), "infeasible" if lab is None
                       else "FEASIBLE", nodes, lem))
    ck("smallest_family_members_infeasible_at_k_2_by_full_exhaustion", ok,
       "; ".join(rows))


def check_corollary(nmax=16):
    """Corollary 3: for every order in [8, nmax] a connected, non-complete-
    multipartite witness with delta >= 2, B = 2 and tvs > 2."""
    rows = []
    ok = True
    for n in range(8, nmax + 1):
        t, d = 1, n - 4
        H = complement_edges(d, cycle_edges(d))
        reg = set(degrees(d, H)) == {d - 3}
        nn, e = join(d, H, 2 * t + 2, [])
        B, _ = nurdin_bound_terms(nn, e)
        comp = set(complement_edges(nn, e))
        want = set(cycle_edges(d)) | set(itertools.combinations(range(d, nn),
                                                               2))
        good = (reg and nn == n and is_connected(nn, e)
                and min(degrees(nn, e)) >= 2 and B == 2
                and endpoint_obstruction(nn, e)
                and not is_complete_multipartite(nn, e)
                and comp == want)
        ok = ok and good
        if not good:
            rows.append("order %d FAILED (reg=%s B=%s)" % (n, reg, B))
    ck("corollary_3_witnesses_of_every_order_from_8_to_%d" % nmax, ok,
       "orders 8..%d each: (d-3)-regular core, B=2, obstruction fires, "
       "complement = C_d + K_4, not complete multipartite" % nmax
       if ok else "; ".join(rows))


def main():
    graphs = check_exhibited_objects()
    check_hypotheses(graphs)
    check_bound(graphs)
    no_two = check_no_two_labeling(graphs)
    check_naive_enumeration(graphs)
    check_sdr_engine()
    certs = check_certificates(graphs)
    check_tvs_value(no_two, certs)
    check_structure(graphs)
    levels, pool = build_census()
    check_census_search(pool, graphs)
    check_lower_bound_tight(pool)
    check_obstruction_lemma(pool)
    check_family()
    check_family_all_cores(levels)
    check_family_bruteforce()
    check_corollary()
    print("NOT RE-RUN HERE: Theorem 2 and Corollary 3 are infinite families, "
          "so only a finite window of each is verified above, and the window "
          "is narrower than the statements in three ways. (i) Theorem 2 is "
          "checked for the admissible (t,d) pairs with t <= 4 and d <= 12 "
          "only. (ii) The theorem asserts its conclusion for EVERY "
          "(d-t-2)-regular graph H on d vertices; every such core up to "
          "isomorphism is tested here only for d <= 7, where the regenerated "
          "census supplies them all, and for 8 <= d <= 12 a single circulant "
          "core per pair is tested instead. (iii) Corollary 3 is checked for "
          "orders 8..16 only, with the paper's own core H = complement(C_d); "
          "order 7 is "
          "covered by the census instead. For all pairs and orders except "
          "(t,d) = (1,4), (2,4), (1,5), tvs(G) > 2 comes from the computed "
          "endpoint-weight obstruction, not from exhaustion; on those three "
          "smallest members, and on every census graph where the obstruction "
          "fires, it is cross-validated against complete exhaustion. "
          "The order <= 7 classification (Theorem 1) is verified in full: the "
          "census is regenerated from scratch, not read from any catalogue.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:                       # keep the output contract
        ck("program_ran_to_completion", False,
           "%s: %s" % (type(exc).__name__, str(exc)[:120]))
    n = len(CHECKS)
    k = sum(1 for _, ok in CHECKS if not ok)
    if k == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
        sys.exit(0)
    print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
    sys.exit(1)

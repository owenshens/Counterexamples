#!/usr/bin/env python3
"""Verification of a claw-free counterexample to zdmg(G) <= 4 dmg(G).

Standard library only; exact integer arithmetic throughout (no floating point
enters any decision).

TAKEN FROM THE PAPER (inputs, transcribed and not otherwise trusted)
  * the construction G_t: a 6-cycle v0..v5, a vertex x adjacent to v0 and v1,
    and a pendant path x-u1-...-ut, so |V(G_t)| = 7 + t;
  * the graph6 strings "HhEM?C@" (G_2) and "IhEM?C@?G" (G_3) under the vertex
    order v0,...,v5,x,u1,...,ut;
  * the five survivor walks T_v0, T_v2, T_v3, T_x, T_u1 and the four
    off-cycle consecutive triples used by the local trace lemma;
  * the stated involution v0<->v1, v2<->v5, v3<->v4 fixing x and the u_j;
  * the recurrence defining the damage value of a state (D, p, s) for an
    ordinary pursuer (move set N[p]) and for a zombie (move set the
    neighbours of p strictly closer to s);
  * the claimed isomorphism-class counts of connected and of claw-free
    connected graphs of orders 1..7 (1,1,2,6,21,112,853 and 1,1,2,5,14,50,191)
    and the claimed 881 claw-free connected classes at order eight.

DERIVED HERE (recomputed from scratch; nothing below is copied)
  * G_2 and G_3 are connected and claw-free, with the claimed orders,
    eccentricity, radius 3 and diameter t+3;
  * the graph6 strings decode to exactly the constructed graphs;
  * dmg(G_2) = dmg(G_3) = 2 and zdmg(G_2) = 9, zdmg(G_3) = 10, by an exact
    least-fixpoint evaluation of the state recurrence over all 2^n damaged
    sets;
  * hence zdmg(G_2) = 9 > 8 = 4 dmg(G_2), refuting the conjecture, and
    zdmg(G_3)/dmg(G_3) = 5 exactly, giving the claimed lower bound on the
    claw-free supremum;
  * every consecutive triple of every exhibited walk is the unique geodesic
    between its endpoints, the survivor portion of each walk covers V, the
    stated map is an involutive automorphism fixing x and the u_j, and the
    walk starts together with their images under it and the u_j for j >= 2
    exhaust V, so every initial zombie vertex is accounted for;
  * isomorphism classes of graphs of orders 1..7 enumerated by orbit
    traversal under transposition generators, matching the paper's counts;
  * no connected claw-free graph of order at most 8 violates the conjecture:
    orders 1..7 by exhausting isomorphism classes, order 8 by extending every
    order-7 claw-free class by one vertex in every way (sound because every
    connected graph has a non-cut vertex and claw-freeness is inherited by
    induced subgraphs), a step itself validated at order 7; the order-8 sweep
    is shown to be complete by counting its isomorphism classes, which must
    equal the 881 of the paper's table;
  * the possible (dmg, zdmg) pairs of an order-nine violator, obtained by
    exhausting the integers in range rather than by transcribing the paper:
    zdmg >= 4 dmg + 1 together with zdmg <= 9 forces dmg <= 2, and dmg = 2
    forces zdmg = 9, so the paper's "(dmg,zdmg) = (2,9)" is derived here for
    every order-nine violator with dmg >= 2 (the dmg <= 1 band is not swept);
  * dmg(G) = 0 exactly for the classes with a universal vertex, and those
    classes also have zdmg(G) = 0, checked on every claw-free connected class
    of order at most 7;
  * positive controls: the sweep predicate fires on the exhibited
    counterexample, the engine returns zero damage on complete graphs, and the
    order-8 sweep provably examined G_2 minus a vertex while the identical
    extend-and-test pipeline, run one order higher on that graph, does report
    violations -- so a null result at order 8 is informative.

NOT RE-RUN (stated by the program at run time)
  * the order-nine census: the 261080 connected and 4494 claw-free class
    counts at order nine, the count of exactly 151 isomorphism classes of
    minimum-order counterexamples, and the paper's statement that all 151 of
    them have (dmg,zdmg) = (2,9).  The order-nine enumeration is credited in
    the paper to an external generator and is not reproduced here; what is
    established below is that no connected claw-free graph of order at most
    eight violates the conjecture, that G_2 does violate it, and that any
    order-nine violator with dmg >= 2 has values exactly (2,9);
  * the paper's "connected classes" entry 11117 at order eight.  The
    order-eight sweep here enumerates only the claw-free connected graphs
    (via extension of the order-seven claw-free classes), so it confirms the
    881 of the claw-free column but never forms the full connected class set
    at that order.  Nothing in the refutation or in the minimality argument
    depends on that entry: orders 1..7 are exhausted over all labeled graphs,
    and completeness of the order-eight sweep is established by the non-cut
    vertex argument plus the 881 count.
"""

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


def summary():
    n = len(CHECKS)
    bad = [nm for nm, ok in CHECKS if not ok]
    if not bad:
        print("VERDICT: ALL %d CHECKS PASS" % n)
        return 0
    print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
    return 1


def build_G(t):
    """G_t: 6-cycle v0..v5, plus x adjacent to v0,v1, plus path x-u1-...-ut."""
    V = ["v%d" % i for i in range(6)] + ["x"] + ["u%d" % j for j in range(1, t + 1)]
    E = set()

    def add(a, b):
        E.add((a, b))
        E.add((b, a))

    for i in range(6):
        add("v%d" % i, "v%d" % ((i + 1) % 6))
    add("x", "v0")
    add("x", "v1")
    prev = "x"
    for j in range(1, t + 1):
        add(prev, "u%d" % j)
        prev = "u%d" % j
    adj = dict((v, frozenset(b for (a, b) in E if a == v)) for v in V)
    return tuple(V), adj


def bfs_dist(V, adj, s):
    dist = {s: 0}
    frontier = [s]
    while frontier:
        nxt = []
        for v in frontier:
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1
                    nxt.append(w)
        frontier = nxt
    return dist


def all_dist(V, adj):
    return dict((v, bfs_dist(V, adj, v)) for v in V)


def is_connected(V, adj):
    return len(bfs_dist(V, adj, V[0])) == len(V)


def is_claw_free(V, adj):
    for c in V:
        nb = sorted(adj[c])
        for tri in itertools.combinations(nb, 3):
            a, b, d = tri
            if b not in adj[a] and d not in adj[a] and d not in adj[b]:
                return False, (c, tri)
    return True, None


def index_graph(V, adj):
    """Return (n, open_nbhd, closed_nbhd, distance matrix) on indices 0..n-1."""
    idx = dict((v, i) for i, v in enumerate(V))
    n = len(V)
    No = [sorted(idx[w] for w in adj[v]) for v in V]
    Nc = [sorted([i] + No[i]) for i in range(n)]
    D = all_dist(V, adj)
    dist = [[D[V[i]].get(V[j]) for j in range(n)] for i in range(n)]
    return n, No, Nc, dist


def zombie_moves(n, No, dist):
    """M_Z(p,s) = {p' in N(p) : d(p',s) = d(p,s) - 1}."""
    return [[[q for q in No[p] if dist[q][s] == dist[p][s] - 1]
             for s in range(n)] for p in range(n)]


def val_ge(n, No, Nc, dist, k, zombie):
    """True iff val_P(G) >= k, by exact least-fixpoint on states (D,p,s).

    A state is winning if for every pursuer move p' in M_P(p,s) we have
    p' != s and either |D u {s}| >= k or some s' in N[s]\\{p'} keeps the
    state (D u {s}, p', s') winning.  Masks are processed in decreasing
    popcount so that every strictly larger damaged set is already final.
    """
    if k <= 0:
        return True
    MZ = zombie_moves(n, No, dist) if zombie else None
    W = set()
    order = sorted(range(1 << n), key=lambda m: -bin(m).count("1"))
    for D in order:
        if bin(D).count("1") >= k:
            continue
        changed = True
        while changed:
            changed = False
            for p in range(n):
                for s in range(n):
                    if s == p or (D, p, s) in W:
                        continue
                    D2 = D | (1 << s)
                    big = bin(D2).count("1") >= k
                    moves = MZ[p][s] if zombie else Nc[p]
                    ok = bool(moves)
                    for p2 in moves:
                        if p2 == s:
                            ok = False
                            break
                        if big:
                            continue
                        if not any((D2, p2, s2) in W
                                   for s2 in Nc[s] if s2 != p2):
                            ok = False
                            break
                    if ok:
                        W.add((D, p, s))
                        changed = True
    return all(any((0, p, s) in W for s in range(n) if s != p)
               for p in range(n))


def game_value(n, No, Nc, dist, zombie):
    k = 0
    while k < n and val_ge(n, No, Nc, dist, k + 1, zombie):
        k += 1
    return k


def graph6_edges(s):
    """Decode a graph6 string into (n, set of frozenset edges, well_formed).

    well_formed is True only if every character is in the printable graph6
    range, the payload has exactly ceil(n(n-1)/2 / 6) characters, and the
    trailing pad bits are zero -- so a string carrying extra or garbage
    payload cannot silently decode to the right edge set.
    """
    if not s or any(not (63 <= ord(ch) <= 126) for ch in s):
        return None, set(), False
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        val = ord(ch) - 63
        for b in range(5, -1, -1):
            bits.append((val >> b) & 1)
    need = n * (n - 1) // 2
    if len(bits) < need:
        return n, set(), False
    well = len(s) - 1 == -(-need // 6) and not any(bits[need:])
    E = set()
    pos = 0
    for j in range(1, n):
        for i in range(j):
            if bits[pos]:
                E.add(frozenset((i, j)))
            pos += 1
    return n, E, well


def edge_set(V, adj):
    return set(frozenset((V.index(v), V.index(w))) for v in V for w in adj[v])


def unique_geodesic_triple(V, adj, a, b, c):
    """True iff a-b-c is the unique geodesic from a to c (so d(a,c)=2 and b
    is the only common neighbour of a and c)."""
    if b not in adj[a] or c not in adj[b] or a == c:
        return False
    if c in adj[a]:
        return False
    common = adj[a] & adj[c]
    return common == frozenset([b])


def walk_is_traceable(V, adj, walk):
    """Every consecutive triple of the walk is a unique geodesic, and every
    consecutive pair is an edge (hypothesis of the local trace lemma)."""
    for i in range(len(walk) - 1):
        if walk[i + 1] not in adj[walk[i]]:
            return False, ("not an edge", walk[i], walk[i + 1])
    for i in range(len(walk) - 2):
        if not unique_geodesic_triple(V, adj, walk[i], walk[i + 1], walk[i + 2]):
            return False, ("triple", walk[i], walk[i + 1], walk[i + 2])
    return True, None


PAIRS_CACHE = {}


def pairs_of(n):
    if n not in PAIRS_CACHE:
        PAIRS_CACHE[n] = [(i, j) for j in range(n) for i in range(j)]
    return PAIRS_CACHE[n]


def _slot_tables(n, perm, chunk=8):
    """Tables applying a vertex permutation to an edge-slot bitmask."""
    prs = pairs_of(n)
    slot = dict((p, k) for k, p in enumerate(prs))
    dest = [slot[(min(perm[i], perm[j]), max(perm[i], perm[j]))]
            for (i, j) in prs]
    m = len(prs)
    tabs = []
    for base in range(0, m, chunk):
        width = min(chunk, m - base)
        tbl = [0] * (1 << width)
        for val in range(1 << width):
            out = 0
            for b in range(width):
                if (val >> b) & 1:
                    out |= 1 << dest[base + b]
            tbl[val] = out
        tabs.append((base, (1 << width) - 1, tbl))
    return tabs


def _apply(tabs, g):
    out = 0
    for base, msk, tbl in tabs:
        out |= tbl[(g >> base) & msk]
    return out


def iso_class_reps(n):
    """One representative per isomorphism class of labeled graphs of order n,
    found by orbit traversal under adjacent-transposition generators."""
    m = n * (n - 1) // 2
    gens = []
    for i in range(n - 1):
        perm = list(range(n))
        perm[i], perm[i + 1] = perm[i + 1], perm[i]
        gens.append(_slot_tables(n, perm))
    seen = bytearray(1 << m)
    reps = []
    for g0 in range(1 << m):
        if seen[g0]:
            continue
        seen[g0] = 1
        reps.append(g0)
        stack = [g0]
        while stack:
            g = stack.pop()
            for tabs in gens:
                h = _apply(tabs, g)
                if not seen[h]:
                    seen[h] = 1
                    stack.append(h)
    return reps


def mask_to_graph(n, g):
    V = tuple(range(n))
    nb = [set() for _ in range(n)]
    for k, (i, j) in enumerate(pairs_of(n)):
        if (g >> k) & 1:
            nb[i].add(j)
            nb[j].add(i)
    adj = dict((i, frozenset(nb[i])) for i in range(n))
    return V, adj


GRAPH6 = {2: "HhEM?C@", 3: "IhEM?C@?G"}

U = {2: ["x", "u1", "u2"], 3: ["x", "u1", "u2", "u3"]}


def walks(t):
    c = lambda *a: list(a)
    tail = U[t]
    return {
        "T_v0": c("v0", "v5", "v4", "v3", "v2", "v1",
                  "v0", "v5", "v4", "v3", "v2", "v1") + tail,
        "T_v2": c("v2", "v1", "v0", "v5", "v4", "v3", "v2", "v1") + tail,
        "T_v3": c("v3", "v2", "v1", "v0", "v5", "v4", "v3", "v2", "v1") + tail,
        "T_x": c("x", "v1", "v2", "v3", "v4", "v5", "v0",
                 "v1", "v2", "v3", "v4", "v5", "v0") + tail,
        "T_u1": c("u1", "x", "v0", "v5", "v4", "v3", "v2", "v1") + tail,
    }


TRIPLES = [("v2", "v1", "x"), ("v5", "v0", "x"),
           ("v1", "x", "u1"), ("v0", "x", "u1")]


def check_structure(t, V, adj):
    ck("order_G%d_is_%d" % (t, 7 + t), len(V) == 7 + t, "|V|=%d" % len(V))
    n6, E6, well6 = graph6_edges(GRAPH6[t])
    ck("graph6_G%d_decodes_to_construction" % t,
       well6 and n6 == len(V) and E6 == edge_set(V, adj),
       "%s -> n=%s, %d edges, well_formed=%s"
       % (GRAPH6[t], n6, len(E6), well6))
    ck("G%d_connected" % t, is_connected(V, adj))
    cf, wit = is_claw_free(V, adj)
    ck("G%d_claw_free" % t, cf, "witness=%s" % (wit,) if wit else "no induced K_1,3")
    deg3 = sorted(v for v in V if len(adj[v]) == 3)
    edge_in_nbhd = all(any(b in adj[a] for a, b in
                           itertools.combinations(sorted(adj[v]), 2))
                       for v in deg3)
    ck("G%d_degree3_vertices_are_v0_v1_x_with_edge_in_neighborhood" % t,
       deg3 == ["v0", "v1", "x"] and max(len(adj[v]) for v in V) == 3
       and edge_in_nbhd, "deg-3 set=%s" % (deg3,))


def check_metrics(t, V, adj):
    D = all_dist(V, adj)
    ecc = dict((v, max(D[v].values())) for v in V)
    rad = min(ecc.values())
    diam = max(ecc.values())
    ck("G%d_ecc_x_is_3_and_radius_3" % t, ecc["x"] == 3 and rad == 3,
       "ecc(x)=%d rad=%d" % (ecc["x"], rad))
    ck("G%d_diameter_is_t_plus_3" % t,
       diam == t + 3 and D["u%d" % t]["v3"] == t + 3,
       "diam=%d d(u%d,v3)=%d" % (diam, t, D["u%d" % t]["v3"]))


def check_trace_argument(t, V, adj):
    ok = True
    bad = None
    for (a, b, c) in TRIPLES:
        if not unique_geodesic_triple(V, adj, a, b, c):
            ok = False
            bad = (a, b, c)
    ck("G%d_four_off_cycle_triples_are_unique_geodesics" % t, ok,
       "bad=%s" % (bad,) if bad else "all 4 verified")
    good = True
    detail = "5 walks"
    for name, w in sorted(walks(t).items()):
        okw, why = walk_is_traceable(V, adj, w)
        if not okw:
            good = False
            detail = "%s: %s" % (name, why)
            break
        if set(w[2:]) != set(V):
            good = False
            detail = "%s misses %s" % (name, sorted(set(V) - set(w[2:])))
            break
    ck("G%d_survivor_walks_are_traceable_and_cover_all_vertices" % t, good,
       detail)
    auto = {"v0": "v1", "v1": "v0", "v2": "v5", "v5": "v2",
            "v3": "v4", "v4": "v3"}
    for v in V:
        auto.setdefault(v, v)
    is_auto = all((auto[b] in adj[auto[a]]) == (b in adj[a])
                  for a in V for b in V)
    involution = all(auto[auto[v]] == v for v in V)
    fixes_tail = all(auto[v] == v for v in V
                     if v == "x" or v.startswith("u"))
    ck("G%d_stated_map_is_an_involutive_automorphism_fixing_x_and_the_u_j" % t,
       is_auto and involution and fixes_tail,
       "auto=%s involution=%s fixes_tail=%s"
       % (is_auto, involution, fixes_tail))
    # Coverage of initial zombie vertices, derived from the walks themselves
    # rather than asserted: a walk T_z starts (as a_{-2}) at the zombie vertex
    # z named in its label, the automorphism transports each such start to its
    # image, and starts u_j with j >= 2 are handled by the separate passing
    # argument.  Those three sources must exhaust V(G_t).
    W = walks(t)
    labels_match = all(name == "T_" + w[0] for name, w in W.items())
    starts = set(w[0] for w in W.values())
    tail_starts = set("u%d" % j for j in range(2, t + 1))
    covered = starts | set(auto[z] for z in starts) | tail_starts
    ck("G%d_walk_starts_and_their_images_cover_every_zombie_start" % t,
       labels_match and covered == set(V),
       "labels_match=%s starts=%s uncovered=%s"
       % (labels_match, sorted(starts), sorted(set(V) - covered)))


def check_game_values(t, V, adj):
    n, No, Nc, dist = index_graph(V, adj)
    d = game_value(n, No, Nc, dist, False)
    ck("dmg_G%d_equals_2" % t, d == 2, "dmg=%d" % d)
    zfull = val_ge(n, No, Nc, dist, n, True)
    ck("zdmg_G%d_equals_%d_every_vertex_damaged" % (t, 7 + t), zfull,
       "zdmg>=%d and zdmg<=|V|=%d" % (n, n))
    z = n if zfull else game_value(n, No, Nc, dist, True)
    rad = min(max(row) for row in dist)
    ck("dmg_G%d_meets_radius_minus_one_lower_bound_with_equality" % t,
       d == rad - 1, "rad=%d dmg=%d" % (rad, d))
    return d, z


PAPER_TABLE = {1: (1, 1), 2: (1, 1), 3: (2, 2), 4: (6, 5),
               5: (21, 14), 6: (112, 50), 7: (853, 191)}


def violates(V, adj):
    """(violation?, dmg).  A violation needs zdmg >= 4*dmg+1, which is
    impossible once 4*dmg+1 exceeds the order, so the zombie fixpoint is
    only evaluated when it can matter."""
    n, No, Nc, dist = index_graph(V, adj)
    d = game_value(n, No, Nc, dist, False)
    if 4 * d + 1 > n:
        return False, d
    return val_ge(n, No, Nc, dist, 4 * d + 1, True), d


def census_small():
    """Exact census by isomorphism class for orders 1..7."""
    counts = {}
    reps = {}
    viol = []
    mono = True
    mono_total = 0
    mono_nontrivial = 0
    zero_char = True
    zero_bad = None
    n_zero = 0
    n_pos = 0
    for n in range(1, 8):
        cc = cf = 0
        keep = []
        for g in iso_class_reps(n):
            V, adj = mask_to_graph(n, g)
            if not is_connected(V, adj):
                continue
            cc += 1
            if not is_claw_free(V, adj)[0]:
                continue
            cf += 1
            keep.append(g)
            bad, d = violates(V, adj)
            if bad:
                viol.append((n, g))
            nn, No, Nc, dist = index_graph(V, adj)
            # Closed form for the bottom of the damage scale: dmg(G) = 0 iff
            # some vertex is universal (the cop starts there and captures on
            # its first move), and such a graph also has zdmg(G) = 0 (the
            # zombie at a universal vertex is one step from the survivor and
            # its only closer neighbour is the survivor itself).  So a graph
            # with dmg = 0 cannot violate the conjecture.
            universal = any(len(adj[v]) == n - 1 for v in V)
            if (d == 0) != universal:
                zero_char = False
                zero_bad = (n, g, d, universal)
            if d == 0:
                n_zero += 1
                if val_ge(nn, No, Nc, dist, 1, True):
                    zero_char = False
                    zero_bad = (n, g, "zdmg >= 1 with dmg = 0")
            else:
                n_pos += 1
            mono_total += 1
            # val_ge(.., k, ..) is True by definition for k <= 0, so only the
            # classes with dmg >= 1 put the zombie engine under any obligation
            # here; the count is reported so the check cannot look stronger
            # than it is.
            if d >= 1:
                mono_nontrivial += 1
                if not val_ge(nn, No, Nc, dist, d, True):
                    mono = False
        counts[n] = (cc, cf)
        reps[n] = keep
    ck("isomorphism_class_counts_orders_1_to_7_match_table",
       counts == PAPER_TABLE, "%s" % sorted(counts.items()))
    ck("no_claw_free_violation_at_orders_1_to_7", not viol,
       "violations=%s" % (viol,))
    ck("engine_monotonicity_zdmg_ge_dmg_on_all_classes_to_order_7", mono,
       "%d of %d claw-free connected classes have dmg >= 1 and so make the "
       "inequality non-vacuous" % (mono_nontrivial, mono_total))
    ck("dmg_zero_iff_universal_vertex_and_then_zdmg_zero_orders_1_to_7",
       zero_char,
       "counterexample=%s" % (zero_bad,) if zero_bad else
       "%d classes with dmg=0 (each has a universal vertex and zdmg=0) and "
       "%d with dmg>=1 (none has a universal vertex)" % (n_zero, n_pos))
    return reps


def extensions(n, base_reps):
    """All graphs of order n+1 obtained by adding one vertex with a nonempty
    neighborhood to a graph of order n given by an edge-slot mask."""
    src = pairs_of(n)
    dst = dict((p, k) for k, p in enumerate(pairs_of(n + 1)))
    out = []
    for g in base_reps:
        base = 0
        for k, (i, j) in enumerate(src):
            if (g >> k) & 1:
                base |= 1 << dst[(i, j)]
        for sub in range(1, 1 << n):
            m = base
            for i in range(n):
                if (sub >> i) & 1:
                    m |= 1 << dst[(i, n)]
            out.append(m)
    return out


def orbit(n, g0):
    gens = []
    for i in range(n - 1):
        perm = list(range(n))
        perm[i], perm[i + 1] = perm[i + 1], perm[i]
        gens.append(_slot_tables(n, perm))
    seen = set([g0])
    stack = [g0]
    while stack:
        g = stack.pop()
        for tabs in gens:
            h = _apply(tabs, g)
            if h not in seen:
                seen.add(h)
                stack.append(h)
    return seen


def _transposition_tables(n):
    gens = []
    for i in range(n - 1):
        perm = list(range(n))
        perm[i], perm[i + 1] = perm[i + 1], perm[i]
        gens.append(_slot_tables(n, perm))
    return gens


def iso_classes_among(n, masks):
    """Number of isomorphism classes represented by a collection of order-n
    edge-slot masks.  Orbits are traversed one at a time and discarded, so
    peak memory is bounded by the largest single orbit rather than by the
    union of all of them."""
    gens = _transposition_tables(n)
    remaining = set(masks)
    classes = 0
    while remaining:
        g0 = next(iter(remaining))
        orb = set([g0])
        stack = [g0]
        while stack:
            g = stack.pop()
            for tabs in gens:
                h = _apply(tabs, g)
                if h not in orb:
                    orb.add(h)
                    stack.append(h)
        remaining -= orb
        classes += 1
    return classes


def graph_mask(n, V, adj):
    """Edge-slot mask of (V, adj) under the given vertex order."""
    slot = dict((p, k) for k, p in enumerate(pairs_of(n)))
    idx = dict((v, i) for i, v in enumerate(V))
    g = 0
    for v in V:
        for w in adj[v]:
            i, j = sorted((idx[v], idx[w]))
            g |= 1 << slot[(i, j)]
    return g


def claw_free_extensions(n, base_reps):
    out = []
    for m in extensions(n, base_reps):
        V, adj = mask_to_graph(n + 1, m)
        if is_claw_free(V, adj)[0] and is_connected(V, adj):
            out.append(m)
    return out


def check_extension_method(reps6, reps7):
    """Every connected claw-free graph of order n+1 has a non-cut vertex whose
    deletion leaves a connected claw-free graph of order n, so extending the
    order-n classes reaches every order-(n+1) class.  Tested at n=6 against
    the exact order-7 classes."""
    cand = set(claw_free_extensions(6, reps6))
    missed = [g for g in reps7 if not (orbit(7, g) & cand)]
    ck("vertex_extension_reaches_every_class_at_order_7",
       not missed and len(reps7) == 191,
       "reached %d/%d order-7 claw-free classes from %d order-6 classes"
       % (len(reps7) - len(missed), len(reps7), len(reps6)))


ORDER8_CLAW_FREE_CLASSES = 881


def check_order_eight(reps7):
    cand = claw_free_extensions(7, reps7)
    viol = []
    for g in cand:
        V, adj = mask_to_graph(8, g)
        bad, d = violates(V, adj)
        if bad:
            viol.append((g, d))
    classes = iso_classes_among(8, cand)
    ck("order_8_claw_free_class_count_matches_table",
       classes == ORDER8_CLAW_FREE_CLASSES,
       "%d isomorphism classes among %d swept labeled graphs (table says %d)"
       % (classes, len(cand), ORDER8_CLAW_FREE_CLASSES))
    ck("no_claw_free_violation_at_order_8",
       not viol and classes == ORDER8_CLAW_FREE_CLASSES,
       "swept %d connected claw-free graphs of order 8 in %d classes, "
       "violations=%d" % (len(cand), classes, len(viol)))
    return cand


def check_sweep_would_find_order_nine(cand8):
    """Negative control for the minimality sweep.

    Deleting the non-cut vertex u_2 from G_2 leaves an order-8 connected
    claw-free graph.  That graph must be among the ones the order-8 sweep
    actually examined, and running the very same extend-and-test pipeline one
    order higher on it must recover violations -- otherwise the order-8 result
    "no violation" would be consistent with a pipeline incapable of ever
    reporting one.
    """
    V, adj = build_G(2)
    V8 = tuple(v for v in V if v != "u2")
    adj8 = dict((v, frozenset(w for w in adj[v] if w != "u2")) for v in V8)
    g8 = graph_mask(8, V8, adj8)
    seen = bool(orbit(8, g8) & set(cand8))
    ck("order_8_sweep_examined_the_counterexample_minus_a_vertex",
       len(V8) == 8 and is_connected(V8, adj8) and is_claw_free(V8, adj8)[0]
       and seen, "G2-u2 connected claw-free and present in sweep=%s" % seen)
    found = []
    for g in claw_free_extensions(8, [g8]):
        W, A = mask_to_graph(9, g)
        bad, d = violates(W, A)
        if bad:
            found.append(d)
    ck("same_pipeline_reports_violations_one_order_higher",
       bool(found) and set(found) == set([2]),
       "%d violating order-9 extensions of G2-u2, dmg values=%s"
       % (len(found), sorted(set(found))))


def check_order_nine_value_constraint(order):
    """Which (dmg, zdmg) pairs an order-nine violator can have, derived.

    The paper reports that all its order-nine counterexamples have
    (dmg, zdmg) = (2, 9).  That census is not re-run here, but the pair itself
    is not a free parameter: a violator satisfies zdmg >= 4 dmg + 1, and zdmg
    counts damaged vertices so zdmg <= |V| = 9.  Exhausting the integers in
    that range leaves (2, 9) as the only pair with dmg >= 2.  Only the dmg <= 1
    band survives unswept, and that shortfall is named in the closing lines.
    """
    feasible = sorted((d, z) for d in range(order + 1)
                      for z in range(order + 1) if z >= 4 * d + 1)
    heavy = sorted(pair for pair in feasible if pair[0] >= 2)
    low = sorted(set(d for d, _ in feasible if d <= 1))
    ck("order_nine_violator_values_forced_to_2_and_9_when_dmg_at_least_2",
       order == 9 and heavy == [(2, 9)],
       "order=%d, feasible (dmg,zdmg) with dmg>=2: %s; unswept residual band "
       "is dmg in %s with zdmg >= 4*dmg+1" % (order, heavy, low))


def check_positive_controls():
    """The sweep predicate must fire on a genuine counterexample, and the game
    engine must reproduce independently known values."""
    bad, d = violates(*build_G(2))
    ck("sweep_predicate_fires_on_the_nine_vertex_counterexample", bad and d == 2,
       "violates(G2)=(%s, dmg=%d)" % (bad, d))
    okk = True
    for n in range(2, 6):
        V = tuple(range(n))
        adj = dict((i, frozenset(j for j in V if j != i)) for i in V)
        nn, No, Nc, dist = index_graph(V, adj)
        if game_value(nn, No, Nc, dist, False) != 0:
            okk = False
        if game_value(nn, No, Nc, dist, True) != 0:
            okk = False
    ck("engine_gives_zero_damage_on_complete_graphs_K2_to_K5", okk,
       "pursuer captures before any damage")


def main():
    check_positive_controls()
    for t in (2, 3):
        V, adj = build_G(t)
        check_structure(t, V, adj)
        check_metrics(t, V, adj)
        check_trace_argument(t, V, adj)
        d, z = check_game_values(t, V, adj)
        if t == 2:
            ck("G2_refutes_zdmg_le_4dmg", z > 4 * d,
               "zdmg=%d > 4*dmg=%d" % (z, 4 * d))
        else:
            from fractions import Fraction
            r = Fraction(z, d)
            ck("G3_ratio_is_exactly_5_so_c3_at_least_5", r == 5,
               "zdmg/dmg = %s" % r)
    reps = census_small()
    check_extension_method(reps[6], reps[7])
    cand8 = check_order_eight(reps[7])
    check_sweep_would_find_order_nine(cand8)
    check_order_nine_value_constraint(len(build_G(2)[0]))
    print("NOT RE-RUN: the order-nine census itself -- the 261080 connected "
          "and 4494 claw-free class counts at order nine, the count of "
          "exactly 151 isomorphism classes of minimum-order counterexamples, "
          "and the paper's statement that all 151 of them have "
          "(dmg,zdmg)=(2,9) -- is beyond a pure-Python budget, and the "
          "order-nine enumeration the paper credits to an external generator "
          "is not reproduced here. Established above instead: no connected "
          "claw-free graph of order at most eight violates the conjecture, "
          "G_2 (order nine) does violate it, and any order-nine violator with "
          "dmg >= 2 has values exactly (2,9). The order-nine violators with "
          "dmg <= 1, if any exist, are not swept here.")
    print("NOT RE-RUN: the paper's 11117 connected classes at order eight. "
          "The order-eight sweep above enumerates only the claw-free "
          "connected graphs, so it checks the 881 of the claw-free column but "
          "not that entry; no claim verified here depends on it.")
    return summary()


if __name__ == "__main__":
    sys.exit(main())

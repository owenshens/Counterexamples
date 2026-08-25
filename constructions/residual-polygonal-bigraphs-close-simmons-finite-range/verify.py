#!/usr/bin/env python3
"""Verification of the two residual polygonal bigraphs and the finite range.

Standard library only, exact integer arithmetic, no input files, no network.

Definitions used (from the paper's statement section, not re-derived):
  P_{m,k} has bipartition X={x_i}, Y={y_i} over Z_m and edges
  a_i = x_i y_i,  b_i = x_i y_{i-1},  c_i = x_i y_{i+k},  1 <= k <= (m-1)/2.
  A graph is 2-edge-Hamiltonian if every two distinct edges lie on a common
  Hamilton cycle.  The residual arithmetic condition after the general
  results quoted in the paper is 2 < gcd(m+1,k+1) < k+1.

TAKEN FROM THE PAPER (inputs, transcribed and then tested):
  * table (3.1): for P_{34,9}, the five pairs (B_r, C_r), with
    A_r = Z_34 minus (B_r union C_r);
  * table (3.2): for P_{32,14}, the five pairs (A_r, C_r), with
    B_r = Z_32 minus (A_r union C_r);
  * the six X-index traces of the appendix, for r = 2,3,4 in each graph;
  * the isomorphism phi(x_i) = x'_{19i}, phi(y_i) = y'_{19i-1};
  * the residual table for m <= 34 with its reduction arrows k -> l (u,d);
  * the paper's stated X-steps -9 (for P_{34,9}) and 15 (for P_{32,14}).

DERIVED HERE (nothing below is taken on the paper's word):
  * each graph is built from the definition and shown to be simple, cubic,
    bipartite, connected, and equal to the 2m-cycle plus offset-(2k+1) chords;
  * each table row is shown to give a partition of Z_m whose induced edge set
    M_r is a perfect matching, and each complement H_r is walked and shown to
    be a single cycle through all 2m vertices;
  * the appendix traces are compared with the walked X-index sequences;
  * H_0 and H_1 are identified with the stated edge classes, and the constant
    X-step of H_1 is computed and compared with the paper's value;
  * LOAD-BEARING: for each of P_{34,9}, P_{32,14}, P_{34,14}, every one of the
    3m(3m-1)/2 unordered pairs of distinct edges is exhibited inside one of
    five 2-factors that are re-verified here to be Hamilton cycles.  This is
    the conclusion "2-edge-Hamiltonian", computed pair by pair;
  * phi is shown to be a graph isomorphism onto P_{34,14}, and the five cycles
    are transported and re-verified inside P_{34,14} itself;
  * the residual set for m <= 34 is recomputed from the gcd condition, every
    reduction arrow is checked to be an isomorphism of graphs, and an
    exhaustive search over all units u and all shifts d determines which
    residual parameters reduce to no non-residual parameter;
  * a from-scratch 2-edge-Hamiltonicity decision for all 272 parameter pairs
    with m <= 34, that is for every polygonal bigraph on at most 68 vertices,
    by Hamilton-cycle search with two forced edges over rotation orbits of
    edge pairs.  This re-proves the corollary without citing anything.
"""

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


def edge_labels(m):
    """The 3m edge labels of P_{m,k}: a_i, b_i, c_i for i in Z_m."""
    return [(t, i) for i in range(m) for t in ("a", "b", "c")]


def ends(m, k, e):
    """Endpoints of an edge label, as ('x',i)/('y',j) vertices."""
    t, i = e
    i %= m
    if t == "a":
        return ("x", i), ("y", i)
    if t == "b":
        return ("x", i), ("y", (i - 1) % m)
    if t == "c":
        return ("x", i), ("y", (i + k) % m)
    raise ValueError("bad edge type")


def vertices(m):
    return [("x", i) for i in range(m)] + [("y", i) for i in range(m)]


def incidence(m, k, edge_set):
    # The edges are consumed in sorted label order, so every incidence list --
    # and hence every walk taken below, and the direction it takes -- is
    # independent of set-iteration order.  Without this, hash randomisation of
    # the 'a'/'b'/'c' component makes the printed X-step of H_1 flip to its
    # negative from one process to the next, and the transcript stops being
    # byte-reproducible.
    inc = dict((v, []) for v in vertices(m))
    for e in sorted(edge_set):
        u, v = ends(m, k, e)
        inc[u].append(e)
        inc[v].append(e)
    return inc


def is_perfect_matching(m, k, edge_set):
    inc = incidence(m, k, edge_set)
    return all(len(inc[v]) == 1 for v in inc)


def cycle_vertex_walk(m, k, edge_set):
    """If edge_set is a single cycle through all 2m vertices, return the
    vertex sequence; otherwise return None."""
    inc = incidence(m, k, edge_set)
    if any(len(inc[v]) != 2 for v in inc):
        return None
    start = ("x", 0)
    seq = [start]
    prev_edge = None
    cur = start
    for _ in range(2 * m):
        cand = [e for e in inc[cur] if e != prev_edge]
        if not cand:
            return None
        e = cand[0]
        u, v = ends(m, k, e)
        nxt = v if u == cur else u
        prev_edge = e
        cur = nxt
        if cur == start:
            break
        if cur in seq:
            return None
        seq.append(cur)
    if cur != start or len(seq) != 2 * m:
        return None
    return seq


def x_trace(seq):
    return [i for (t, i) in seq if t == "x"]


def cyclic_same(u, v):
    """True if the cyclic sequence u equals v up to rotation or reflection."""
    if len(u) != len(v):
        return False
    n = len(u)
    for w in (v, v[::-1]):
        for s in range(n):
            if all(u[t] == w[(s + t) % n] for t in range(n)):
                return True
    return False


def shift(S, t, m):
    return set((s + t) % m for s in S)


def partitions_zm(m, parts):
    """True iff the given sets are pairwise disjoint with union Z_m."""
    tot = 0
    seen = set()
    for P in parts:
        Q = set(x % m for x in P)
        if len(Q) != len(P):
            return False
        tot += len(Q)
        seen |= Q
    return tot == m and seen == set(range(m))


def matching_from_rows(m, A, B, C):
    return set([("a", i) for i in A] + [("b", i) for i in B]
               + [("c", i) for i in C])


# ---------------------------------------------------------------- paper data

# P_{34,9}: the table gives B_r and C_r, with A_r = Z_34 \ (B_r u C_r).
BC_34_9 = [
    (set(), set(range(34))),
    (set(range(34)), set()),
    ({25, 26}, {8, 17, 24, 33}),
    ({16, 24}, {7, 15, 23, 32}),
    ({15, 23}, {6, 14, 22, 31}),
]

# P_{32,14}: the table gives A_r and C_r, with B_r = Z_32 \ (A_r u C_r).
AC_32_14 = [
    (set(), set(range(32))),
    (set(range(32)), set()),
    ({0, 16}, {1, 17}),
    ({15, 31}, {0, 16}),
    ({14, 30}, {15, 31}),
]

TRACES_34_9 = {
    2: [0, 10, 20, 30, 6, 16, 25, 1, 11, 21, 31, 7, 17, 18, 28, 4, 14,
        24, 15, 5, 29, 19, 9, 8, 32, 22, 12, 2, 26, 27, 3, 13, 23, 33],
    3: [0, 10, 20, 30, 6, 15, 5, 29, 19, 9, 33, 32, 22, 12, 2, 26, 16,
        17, 27, 3, 13, 23, 14, 4, 28, 18, 8, 7, 31, 21, 11, 1, 25, 24],
    4: [0, 10, 20, 30, 6, 7, 17, 27, 3, 13, 22, 12, 2, 26, 16, 15, 25,
        1, 11, 21, 31, 32, 8, 18, 28, 4, 14, 5, 29, 19, 9, 33, 23, 24],
}

TRACES_32_14 = {
    2: [0, 14, 28, 10, 24, 6, 20, 2, 17, 3, 21, 7, 25, 11, 29, 15,
        16, 30, 12, 26, 8, 22, 4, 18, 1, 19, 5, 23, 9, 27, 13, 31],
    3: [0, 18, 4, 22, 8, 26, 12, 30, 31, 13, 27, 9, 23, 5, 19, 1,
        16, 2, 20, 6, 24, 10, 28, 14, 15, 29, 11, 25, 7, 21, 3, 17],
    4: [0, 18, 4, 22, 8, 26, 12, 30, 29, 11, 25, 7, 21, 3, 17, 31,
        16, 2, 20, 6, 24, 10, 28, 14, 13, 27, 9, 23, 5, 19, 1, 15],
}


def rows_34_9():
    out = []
    for B, C in BC_34_9:
        out.append((set(range(34)) - (B | C), set(B), set(C)))
    return out


def rows_32_14():
    out = []
    for A, C in AC_32_14:
        out.append((set(A), set(range(32)) - (A | C), set(C)))
    return out


def check_wellformed(m, k, name):
    """Decode P_{m,k} and confirm it is what the paper says it is."""
    E = edge_labels(m)
    pairs = [frozenset(ends(m, k, e)) for e in E]
    inc = incidence(m, k, E)
    simple = len(set(pairs)) == len(pairs) and all(len(p) == 2 for p in pairs)
    cubic = all(len(inc[v]) == 3 for v in inc)
    # bipartition computed by breadth-first 2-colouring, not assumed
    col = {("x", 0): 0}
    stack = [("x", 0)]
    bad = False
    while stack:
        u = stack.pop()
        for e in inc[u]:
            p, q = ends(m, k, e)
            w = q if p == u else p
            if w not in col:
                col[w] = 1 - col[u]
                stack.append(w)
            elif col[w] == col[u]:
                bad = True
    connected = len(col) == 2 * m
    krange = 1 <= k <= (m - 1) // 2
    ok = (len(inc) == 2 * m and len(E) == 3 * m and simple and cubic
          and not bad and connected and krange)
    ck("wellformed_%s" % name, ok,
       "|V|=%d |E|=%d cubic=%s bipartite=%s connected=%s 1<=k<=%d"
       % (len(inc), len(E), cubic, not bad, connected, (m - 1) // 2))
    # the paper's alternative description: the 2m-cycle plus offset-(2k+1)
    # chords, under x_i -> 2i, y_i -> 2i+1 in Z_{2m}
    def emb(v):
        return (2 * v[1]) % (2 * m) if v[0] == "x" else (2 * v[1] + 1) % (2 * m)
    ring = set(frozenset(((j, (j + 1) % (2 * m)))) for j in range(2 * m))
    got_ring = set(frozenset((emb(p), emb(q)))
                   for e in E if e[0] in ("a", "b")
                   for p, q in [ends(m, k, e)])
    chords = set(frozenset((emb(p), emb(q)))
                 for e in E if e[0] == "c" for p, q in [ends(m, k, e)])
    want_chords = set(frozenset((2 * i, (2 * i + 2 * k + 1) % (2 * m)))
                      for i in range(m))
    ck("cycle_plus_offset_chords_%s" % name,
       got_ring == ring and chords == want_chords and len(chords) == m,
       "2m-cycle from a,b edges; %d chords of offset %d" % (len(chords), 2 * k + 1))


def build_cycles(m, k, rows, name, traces):
    """Turn the paper's table into the five 2-factors H_r and check each is a
    perfect-matching complement that is a single Hamilton cycle."""
    E = set(edge_labels(m))
    cycles = []
    part_ok, pm_ok, ham_ok, tr_ok = True, True, True, True
    sizes = []
    for r, (A, B, C) in enumerate(rows):
        if not partitions_zm(m, [A, B, C]):
            part_ok = False
        if not partitions_zm(m, [A, shift(B, -1, m), shift(C, k, m)]):
            part_ok = False
        M = matching_from_rows(m, A, B, C)
        if len(M) != m or not is_perfect_matching(m, k, M):
            pm_ok = False
        H = E - M
        w = cycle_vertex_walk(m, k, H)
        if w is None:
            ham_ok = False
            cycles.append(None)
            if r in traces:
                tr_ok = False   # a trace cannot certify a non-cycle
            continue
        sizes.append(len(w))
        cycles.append(H)
        if r in traces:
            T = traces[r]
            if sorted(T) != list(range(m)) or not cyclic_same(x_trace(w), T):
                tr_ok = False
    nrows = len(rows)
    ck("row_sets_partition_%s" % name, part_ok,
       "A_r,B_r,C_r and A_r,B_r-1,C_r+%d each partition Z_%d for r=0..%d"
       % (k, m, nrows - 1))
    ck("perfect_matchings_%s" % name, pm_ok,
       "each M_r has %d edges and covers every vertex once" % m)
    ck("hamilton_cycles_%s" % name, ham_ok and len(sizes) == 5
       and all(s == 2 * m for s in sizes),
       "H_0..H_%d are single cycles on all %d vertices" % (nrows - 1, 2 * m))
    ck("appendix_traces_%s" % name, tr_ok,
       "traces for r=%s are permutations of Z_%d and match the walked "
       "X-index sequences up to rotation/reflection"
       % (sorted(traces), m))
    return cycles


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def check_h0_h1(m, k, cycles, name, h0_types, h1_types, step):
    """H_0 is the outer 2m-cycle; H_1 is the union of two edge classes whose
    successive X-indices advance by a fixed step coprime to m."""
    want0 = set(e for e in edge_labels(m) if e[0] in h0_types)
    want1 = set(e for e in edge_labels(m) if e[0] in h1_types)
    ck("H0_is_outer_cycle_%s" % name,
       cycles[0] == want0 and cycle_vertex_walk(m, k, want0) is not None,
       "H_0 = all %s-edges, a single %d-cycle" % ("/".join(h0_types), 2 * m))
    w = cycle_vertex_walk(m, k, cycles[1]) if cycles[1] else None
    if w is None:
        ck("H1_fixed_step_%s" % name, False, "H_1 is not a Hamilton cycle")
        return
    tr = x_trace(w)
    diffs = set((tr[(t + 1) % m] - tr[t]) % m for t in range(m))
    d = diffs.pop() if len(diffs) == 1 else None
    ok = (cycles[1] == want1 and d is not None
          and step % m in (d, (-d) % m) and gcd(m, step % m) == 1)
    ck("H1_fixed_step_%s" % name, ok,
       "H_1 = all %s-edges, X-step %s (reverse orientation %s; paper: %+d = %d "
       "mod %d), gcd(%d,%d)=%d"
       % ("/".join(h1_types), d, None if d is None else (-d) % m, step,
          step % m, m, m, abs(step), gcd(m, abs(step) % m)))


def check_disjoint_columns(m, rows, name, cols):
    """The paper's three-cycle argument needs two of the three columns to be
    pairwise disjoint across r = 2,3,4."""
    idx = {"A": 0, "B": 1, "C": 2}
    ok = True
    detail = []
    for cname in cols:
        S = [rows[r][idx[cname]] for r in (2, 3, 4)]
        pd = all(not (S[i] & S[j]) for i in range(3) for j in range(i + 1, 3))
        detail.append("%s_2,%s_3,%s_4 disjoint=%s" % (cname, cname, cname, pd))
        ok = ok and pd
    ck("disjoint_columns_%s" % name, ok, "; ".join(detail))


def check_two_edge_hamiltonian(m, k, cycles, name, remaining_types):
    """LOAD-BEARING: every unordered pair of distinct edges lies on a common
    Hamilton cycle.  Only 2-factors that are re-verified here as Hamilton
    cycles may be used as witnesses."""
    good = []
    for H in cycles:
        if H is not None and cycle_vertex_walk(m, k, H) is not None:
            good.append(H)
    E = edge_labels(m)
    total = 0
    uncovered = []
    witness = {}
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            e, f = E[i], E[j]
            total += 1
            hit = None
            for r, H in enumerate(good):
                if e in H and f in H:
                    hit = r
                    break
            if hit is None:
                uncovered.append((e, f))
            else:
                witness[(e, f)] = hit
    ck("two_edge_hamiltonian_%s" % name,
       len(good) == 5 and not uncovered and total == 3 * m * (3 * m - 1) // 2,
       "%d/%d edge pairs lie on one of %d verified Hamilton cycles"
       % (total - len(uncovered), total, len(good)))
    # the paper's reduction: H_0 and H_1 already cover everything except the
    # pairs of the two named types
    # H_0 and H_1 are addressed by position, not by "the first two that
    # happened to verify", so the claim is about the cycles the paper names.
    h0 = cycles[0] if len(cycles) > 0 else None
    h1 = cycles[1] if len(cycles) > 1 else None
    ok01 = (h0 is not None and h1 is not None
            and cycle_vertex_walk(m, k, h0) is not None
            and cycle_vertex_walk(m, k, h1) is not None)
    left = set()
    if ok01:
        for i in range(len(E)):
            for j in range(i + 1, len(E)):
                e, f = E[i], E[j]
                if not ((e in h0 and f in h0) or (e in h1 and f in h1)):
                    left.add(tuple(sorted((e[0], f[0]))))
    ck("residual_pairs_are_%s%s_%s" % (remaining_types[0],
                                       remaining_types[1], name),
       ok01 and left == set([tuple(sorted(remaining_types))]),
       "pairs missed by H_0 and H_1 have types %s"
       % (sorted("".join(t) for t in left) if ok01
          else "<H_0 and H_1 are not both verified Hamilton cycles>"))
    return witness


def iso_map(m, u, d):
    """The map x_i -> x'_{ui}, y_j -> y'_{uj+d}."""
    def f(v):
        t, i = v
        return (t, (u * i + (d if t == "y" else 0)) % m)
    return f


def is_iso(m, k, l, u, d):
    """True iff x_i -> x'_{ui}, y_j -> y'_{uj+d} is an isomorphism
    P_{m,k} -> P_{m,l}; both directions of the edge-set image are checked."""
    if gcd(u % m, m) != 1:
        return False
    f = iso_map(m, u, d)
    if len(set(f(v) for v in vertices(m))) != 2 * m:
        return False
    src = set(frozenset(map(f, ends(m, k, e))) for e in edge_labels(m))
    dst = set(frozenset(ends(m, l, e)) for e in edge_labels(m))
    return src == dst and len(src) == 3 * m


def check_phi_34():
    """The paper's explicit isomorphism P_{34,9} = P_{34,14}."""
    m, u, d = 34, 19, -1
    inv = [t for t in range(m) if (u * t) % m == 1]
    ck("19_is_a_unit_mod_34", len(inv) == 1,
       "19 * %d = 1 mod 34" % (inv[0] if inv else -1))
    ck("phi_is_isomorphism_P34_9_to_P34_14", is_iso(m, 9, 14, u, d),
       "x_i -> x'_{19i}, y_i -> y'_{19i-1} carries E(P_34,9) onto E(P_34,14)")
    f = iso_map(m, u, d)
    img = {"a": "b", "b": "c", "c": "a"}
    ok = True
    for t in ("a", "b", "c"):
        for i in range(m):
            got = frozenset(map(f, ends(m, 9, (t, i))))
            want = frozenset(ends(m, 14, (img[t], (u * i) % m)))
            if got != want:
                ok = False
    ck("phi_edge_class_images", ok,
       "phi(a_i)=b'_{19i}, phi(b_i)=c'_{19i}, phi(c_i)=a'_{19i} for all i")
    cong = all((u * (i - 1) - 1) % m == (u * i + 14) % m
               and (u * (i + 9) - 1) % m == (u * i) % m for i in range(m))
    ck("phi_congruences", cong,
       "19(i-1)-1 = 19i+14 and 19(i+9)-1 = 19i mod 34 for all i")
    ck("S_k_shift_identity_9_to_14",
       set((u * s + d) % m for s in (0, -1, 9))
       == set(s % m for s in (0, -1, 14)),
       "19*S_9 - 1 = S_14 in Z_34")


def transported_cycles_34_14(cycles):
    """Push the five P_{34,9} cycles through phi and relabel them as edges of
    P_{34,14}; the images are re-verified inside P_{34,14} itself."""
    m, u, d = 34, 19, -1
    f = iso_map(m, u, d)
    lookup = dict((frozenset(ends(m, 14, e)), e) for e in edge_labels(m))
    out = []
    for H in cycles:
        if H is None:
            out.append(None)
            continue
        img = set()
        bad = False
        for e in H:
            key = frozenset(map(f, ends(m, 9, e)))
            if key not in lookup:
                bad = True
                break
            img.add(lookup[key])
        out.append(None if bad or len(img) != len(H) else img)
    return out


def int_model(m, k):
    """P_{m,k} on integer vertices: x_i -> i, y_j -> m+j."""
    n = 2 * m
    adj = [[] for _ in range(n)]
    lab = {}
    for i in range(m):
        for t in ("a", "b", "c"):
            u, v = ends(m, k, (t, i))
            a, b = u[1], m + v[1]
            adj[a].append(b)
            adj[b].append(a)
            lab[(a, b)] = (t, i)
            lab[(b, a)] = (t, i)
    return adj, lab


def ham_with_forced(m, adj, forced):
    """Exhaustive search for a Hamilton cycle through both edges of `forced`
    (a list of two vertex pairs).  Returns the vertex cycle, or None if the
    search space is exhausted without success."""
    n = 2 * m
    fnb = [[] for _ in range(n)]
    for (u, v) in forced:
        fnb[u].append(v)
        fnb[v].append(u)
    if any(len(f) > 2 for f in fnb):
        return None
    start, second = forced[0]
    visited = [False] * n
    visited[start] = True
    visited[second] = True
    path = [start, second]

    def deg_ok():
        cur = path[-1]
        for w in range(n):
            if visited[w]:
                continue
            c = 0
            for z in adj[w]:
                if not visited[z] or z == start or z == cur:
                    c += 1
            if c < 2:
                return False
        return True

    def rec(cur, prev):
        if len(path) == n:
            if start not in adj[cur]:
                return False
            if any(z != prev and z != start for z in fnb[cur]):
                return False
            if any(z != path[1] and z != cur for z in fnb[start]):
                return False
            return True
        if len(fnb[cur]) == 2 and prev not in fnb[cur]:
            return False
        f = [z for z in fnb[cur] if z != prev]
        cands = [w for w in adj[cur] if not visited[w]]
        if f:
            cands = [w for w in cands if w in f]
        for w in cands:
            visited[w] = True
            path.append(w)
            if deg_ok() and rec(w, cur):
                return True
            path.pop()
            visited[w] = False
        return False

    return list(path) if rec(second, start) else None


def cover_offsets(m, H):
    """For a Hamilton cycle H, the offsets d such that some rotation of H
    contains both (s,0) and (t,d)."""
    by = {"a": [], "b": [], "c": []}
    for (t, i) in H:
        by[t].append(i)
    D = {}
    for s in ("a", "b", "c"):
        for t in ("a", "b", "c"):
            D[(s, t)] = set((j - i) % m for i in by[s] for j in by[t])
    return D


def pair_orbit_reps(m):
    return [(s, t, d) for s in ("a", "b", "c") for t in ("a", "b", "c")
            for d in range(m) if not (s == t and d == 0)]


def census_two_eh(m, k):
    """Decide 2-edge-Hamiltonicity of P_{m,k} from scratch, using the
    rotation x_i->x_{i+1}, y_i->y_{i+1} to reduce edge pairs to orbits.
    Returns (True, ncycles) or (False, reason)."""
    adj, lab = int_model(m, k)
    if len(lab) != 6 * m:
        return (False, "not simple")
    off = {"a": 0, "b": -1, "c": k}
    pool = []
    for (s, t, d) in pair_orbit_reps(m):
        if any(d in D[(s, t)] for D in pool):
            continue
        forced = [(0, m + off[s] % m), (d % m, m + (d + off[t]) % m)]
        cyc = ham_with_forced(m, adj, forced)
        if cyc is None:
            return (False, "no Hamilton cycle through %s_0 and %s_%d"
                    % (s, t, d))
        H = set(lab[(cyc[i], cyc[(i + 1) % len(cyc)])]
                for i in range(len(cyc)))
        if cycle_vertex_walk(m, k, H) is None:
            return (False, "search returned a non-cycle")
        if (s, 0) not in H or (t, d % m) not in H:
            return (False, "search dropped a forced edge")
        pool.append(cover_offsets(m, H))
    return (True, len(pool))


# The residual table of the paper: m -> list of (k, ell, u, d), with
# ell = None for the entries the table leaves unreduced.
RESIDUAL_TABLE = {
    14: [(5, 3, 11, 0)],
    19: [(7, None, None, None)],
    20: [(5, None, None, None), (8, None, None, None)],
    23: [(8, 3, 20, 0)],
    24: [(9, None, None, None)],
    26: [(5, 4, 21, -1), (11, 6, 19, -1)],
    27: [(7, 4, 23, 0), (11, 5, 22, 0)],
    29: [(8, 12, 13, 12), (11, 12, 12, 12)],
    31: [(11, 12, 13, 12)],
    32: [(5, 13, 19, 0), (8, 7, 7, 7), (11, 3, 29, 0),
         (14, None, None, None)],
    34: [(13, 12, 21, -1), (9, 14, 19, -1)],
}

MMAX = 34
PROVED_HERE = set([(32, 14), (34, 9), (34, 14)])
CITED_SETTLED = set([(19, 7), (20, 5), (20, 8), (24, 9)])


def residual_census(mmax):
    """{(m,k) : 2 < gcd(m+1,k+1) < k+1}, the paper's residual condition."""
    out = set()
    for m in range(3, mmax + 1):
        for k in range(1, (m - 1) // 2 + 1):
            g = gcd(m + 1, k + 1)
            if 2 < g < k + 1:
                out.add((m, k))
    return out


def reduction_targets(m, k):
    """All ell in the canonical range reachable by u*S_k + d = S_ell."""
    want = dict((frozenset(x % m for x in (0, -1, l)), l)
                for l in range(1, (m - 1) // 2 + 1))
    out = set()
    for u in range(1, m):
        if gcd(u, m) != 1:
            continue
        for d in range(m):
            key = frozenset((u * x + d) % m for x in (0, -1, k))
            if key in want:
                out.add(want[key])
    return out


def check_residual_table():
    """Re-derive the paper's residual table for m <= 34 and its reductions."""
    cen = residual_census(MMAX)
    src = set((m, e[0]) for m in RESIDUAL_TABLE for e in RESIDUAL_TABLE[m])
    extra = cen - src
    ck("residual_census_matches_table",
       src <= cen and extra == set([(34, 14)]),
       "%d residual parameters with m<=34; every table entry is residual; "
       "the only residual parameter absent from the table's left column is "
       "P_34,14, the target of the arrow 9->14" % len(cen))
    ok_arrow, ok_target = True, True
    detail, tdetail = [], []
    n_arrows = 0
    for m in sorted(RESIDUAL_TABLE):
        for (k, l, u, d) in RESIDUAL_TABLE[m]:
            if l is None:
                continue
            n_arrows += 1
            shifted = set((u * x + d) % m for x in (0, -1, k))
            if (shifted != set(x % m for x in (0, -1, l))
                    or not (1 <= l <= (m - 1) // 2)
                    or not is_iso(m, k, l, u, d)):
                ok_arrow = False
                detail.append("arrow %d:%d->%d fails" % (m, k, l))
            if ((m, l) in cen) != ((m, k, l) == (34, 9, 14)):
                ok_target = False
                tdetail.append("target %d:%d->%d residual" % (m, k, l))
    ck("reduction_arrows_are_isomorphisms", ok_arrow,
       "each arrow k->l (u,d) satisfies u*S_k+d=S_l and x_i->x'_{ui}, "
       "y_j->y'_{uj+d} maps E(P_m,k) onto E(P_m,l); %s"
       % ("; ".join(detail) if detail else
          "all %d arrows of the table verified" % n_arrows))
    ck("arrow_targets_leave_the_residual_set", ok_target,
       "every arrow target is non-residual except 9->14, whose target "
       "P_34,14 is residual and is settled by the theorem; %s"
       % ("; ".join(tdetail) if tdetail else
          "all %d arrow targets classified" % n_arrows))
    irred = set()
    for (m, k) in cen:
        if not any((m, l) not in cen for l in reduction_targets(m, k)):
            irred.add((m, k))
    ck("irreducible_residual_parameters",
       irred == PROVED_HERE | CITED_SETTLED,
       "exhaustive search over all units u and shifts d: the residual "
       "parameters with m<=34 that reduce to no non-residual parameter are "
       "exactly %s" % sorted(irred))
    ck("P32_14_admits_no_reduction",
       reduction_targets(32, 14) == set([14])
       and 12 in reduction_targets(31, 11)
       and (31, 12) not in cen and reduction_targets(32, 8) == set([7, 8])
       and (32, 7) not in cen,
       "P_32,14 reduces only to itself, while P_31,11 -> P_31,12 and "
       "P_32,8 -> P_32,7 land outside the residual set")
    ck("P34_9_class_is_exactly_9_and_14",
       reduction_targets(34, 9) == set([9, 14])
       and reduction_targets(34, 14) == set([9, 14]),
       "the isomorphism class of P_34,9 inside the m=34 parameters is "
       "{9,14}, and both members are residual")


def check_rotation_and_orbits(cases):
    """The rotation used to reduce edge pairs to orbits really is an
    automorphism, and the orbit representatives really do exhaust the pairs."""
    ok_rot, ok_orb = True, True
    for (m, k) in cases:
        for i in range(m):
            for t in ("a", "b", "c"):
                p, q = ends(m, k, (t, i))
                img = frozenset([(p[0], (p[1] + 1) % m), (q[0], (q[1] + 1) % m)])
                if img != frozenset(ends(m, k, (t, (i + 1) % m))):
                    ok_rot = False
        reps = set(pair_orbit_reps(m))
        E = edge_labels(m)
        for i in range(len(E)):
            for j in range(i + 1, len(E)):
                (s, p), (t, q) = E[i], E[j]
                r = (s, t, (q - p) % m)
                if r not in reps:
                    ok_orb = False
                    continue
                # rotating the representative pair by p must return the pair
                rot = lambda v: (v[0], (v[1] + p) % m)
                got = set([frozenset(map(rot, ends(m, k, (s, 0)))),
                           frozenset(map(rot, ends(m, k, (t, r[2]))))])
                want = set([frozenset(ends(m, k, (s, p))),
                            frozenset(ends(m, k, (t, q)))])
                if got != want:
                    ok_orb = False
    ck("rotation_is_an_automorphism", ok_rot,
       "x_i->x_{i+1}, y_i->y_{i+1} sends a_i,b_i,c_i to a_{i+1},b_{i+1},"
       "c_{i+1} in all %d graphs P_m,k with m<=%d"
       % (len(cases), max(c[0] for c in cases)))
    ck("pair_orbits_are_complete", ok_orb,
       "in all %d graphs, every unordered pair of distinct edges is the "
       "rotation by p of a representative {(s,0),(t,d)}; this is the "
       "hypothesis the census below uses to skip covered orbits"
       % len(cases))


def check_range_arithmetic():
    # vertex counts read off the constructed graphs, not off the paper
    n32 = len(incidence(32, 14, edge_labels(32)))
    n34 = len(incidence(34, 9, edge_labels(34)))
    ck("vertex_counts", n32 == 64 and n34 == 68 and n32 <= 66 < n34,
       "P_32,14 has %d vertices (inside the asserted 66-vertex range, so the "
       "omission is a genuine gap); P_34,9 and P_34,14 have %d > 66"
       % (n32, n34))
    ms = [m for m in range(3, 200) if 2 * m <= 68]
    count = sum((m - 1) // 2 for m in ms)
    ck("range_m_le_34_is_68_vertices",
       ms == list(range(3, MMAX + 1)) and count == 272,
       "2m <= 68 iff 3 <= m <= 34, giving %d parameter pairs (m,k)" % count)
    s_hits = [s for s in range(1, 100) if (8 * s, 4 * s - 2) == (32, 14)]
    ck("family_8s_4s_minus_2_at_s_4", s_hits == [4] and 34 % 8 != 0,
       "P_8s,4s-2 equals P_32,14 exactly at s=4, and no s gives m=34, so the "
       "class P_34,9 = P_34,14 is not a member of that family")


def check_census(mmax):
    """As much of the corollary as pure Python can re-run: an independent
    2-edge-Hamiltonicity decision for every parameter pair with m <= mmax."""
    bad = []
    ncyc = 0
    total = 0
    for m in range(3, mmax + 1):
        for k in range(1, (m - 1) // 2 + 1):
            total += 1
            ok, info = census_two_eh(m, k)
            if ok:
                ncyc += info
            else:
                bad.append("P_%d,%d: %s" % (m, k, info))
    ck("census_all_polygonal_bigraphs_m_le_%d" % mmax,
       not bad and total == 272,
       "%d/%d graphs verified 2-edge-Hamiltonian from scratch (%d Hamilton "
       "cycles found and independently re-verified)%s"
       % (total - len(bad), total, ncyc,
          "" if not bad else "; failures: " + "; ".join(bad[:3])))


def main():
    print("-- the three graphs, decoded from the definition")
    check_wellformed(34, 9, "P34_9")
    check_wellformed(32, 14, "P32_14")
    check_wellformed(34, 14, "P34_14")
    # the orbit reduction is a hypothesis of the from-scratch census at the
    # end, so it is verified for every parameter pair the census will touch
    check_rotation_and_orbits([(m, k) for m in range(3, MMAX + 1)
                               for k in range(1, (m - 1) // 2 + 1)])

    print("-- certificate for P_34,9 (table 3.1 and its traces)")
    r9 = rows_34_9()
    c9 = build_cycles(34, 9, r9, "P34_9", TRACES_34_9)
    check_h0_h1(34, 9, c9, "P34_9", ("a", "b"), ("a", "c"), -9)
    check_disjoint_columns(34, r9, "P34_9", ("B", "C"))
    check_two_edge_hamiltonian(34, 9, c9, "P34_9", ("b", "c"))

    print("-- certificate for P_32,14 (table 3.2 and its traces)")
    r14 = rows_32_14()
    c14 = build_cycles(32, 14, r14, "P32_14", TRACES_32_14)
    check_h0_h1(32, 14, c14, "P32_14", ("a", "b"), ("b", "c"), 15)
    check_disjoint_columns(32, r14, "P32_14", ("A", "C"))
    check_two_edge_hamiltonian(32, 14, c14, "P32_14", ("a", "c"))

    print("-- the isomorphism P_34,9 = P_34,14 and the transported cycles")
    check_phi_34()
    t14 = transported_cycles_34_14(c9)
    check_two_edge_hamiltonian(34, 14, t14, "P34_14", ("a", "c"))

    print("-- the residual table and its reductions, for m <= 34")
    check_residual_table()
    check_range_arithmetic()

    print("-- independent census (this re-proves the corollary from scratch)")
    check_census(MMAX)

    print("NOT RE-RUN: what the originally printed residual table did or did "
          "not display cannot be re-derived from arithmetic; the checks above "
          "establish which parameters are residual and which of them reduce, "
          "not which lines a printed table contained.")
    print("NOT RE-RUN: no claim is made or tested here for m > 34; the "
          "statement verified is the finite one.")


if __name__ == "__main__":
    main()
    n = len(CHECKS)
    k = sum(1 for _, o in CHECKS if not o)
    if k:
        print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)

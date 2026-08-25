#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- verification program for the paper

    "Counterexamples to the Xu--He--Lu Conjecture on Partite Cycle Saturation"
    (15-edge C_8-saturated subgraph of K_{6,6} disproves Xu-He-Lu Conjecture 6.2)

Claim under test (Theorem 1, and the abstract):
    For every integer l >= 4 there is a C_{2l}-saturated spanning subgraph
    H_l of K_{2l-2,2l-2} with e(H_l) = l^2 - 1.  Xu-He-Lu Conjecture 6.2 asserts
    sat(K_{n1,n2},C_{2l}) = n1+n2+l^2-3l+1 for l >= 4 and n1,n2 >= l+2; at
    n1 = n2 = 2l-2 that value is l^2+l-3, which exceeds l^2-1 by l-2.  So the
    conjecture is false for every l >= 4.  Smallest case l = 4: 15 edges in
    K_{6,6} where the conjectured value is 17.

--------------------------------------------------------------------------------
TAKEN FROM THE PAPER (data; not re-derived here)
--------------------------------------------------------------------------------
 T1. The construction recipe of H_l (paper, Section 2, path Q of eq. (3)):
       r = l-1;  A_0 = {a_0..a_{r-1}}, B_0 = {b_0..b_{r-1}}, X = {x_1..x_r},
       Y = {y_1..y_r};  parts A = A_0 u X, B = B_0 u Y;
       R = K_{r,r} - a_0 b_0 on A_0 u B_0;  plus the path
       Q: a_0 - y_1 - x_1 - y_2 - x_2 - ... - y_r - x_r - b_0.
     This is the exhibited object.  Its edge set is read off the definition; every
     property of it is computed below.
 T2. The paper's asserted edge count formula  e(H_l) = l^2 - 1.
 T3. The Xu-He-Lu formula quoted by the paper: n1 + n2 + l^2 - 3l + 1
     (their Theorem 1.1 upper bound, and the equality of their Conjecture 6.2),
     with hypothesis l >= 4 and n1,n2 >= l+2.
 T4. The paper's claimed deficit l-2 at n1 = n2 = 2l-2.
 T5. The candidate automorphism sigma of the proof:
       a_0<->b_0, a_k<->b_k (1<=k<=r-1), x_i<->y_{r-i+1} (1<=i<=r).
 T6. Xu-He-Lu's Definition 2.2 construction G^l_{n1,n2}, transcribed from the
     e-print source (arXiv:2410.11194v2) for calibration control (d):
       |V_i| = n_i, A_i c V_i with |A_i| = n_i - l, B_i = V_i \ A_i (so |B_i| = l),
       distinct x_1,x_2 in B_1, y_1 in B_2; edges are exactly
         (1) u y_1 for all u in A_1,   (2) v x_1 for all v in A_2,
         (3) y_1 x_2,                 (4) all of B_1 x (B_2 \ {y_1}).
 T7. Calibration values fixed in advance here, used only as comparison
     targets: cycle-length spectrum of H_4 = {4,6,10,12}, of H_5 =
     {4,6,8,12,14,16}; K_{4,4} has a C_8 but no C_10; H_4 has 21 host non-edges;
     e(G^4_{6,6}) = 17, e(G^4_{7,6}) = 18, e(G^5_{7,7}) = 25.

--------------------------------------------------------------------------------
DERIVED HERE (this is what the checks actually decide)
--------------------------------------------------------------------------------
 D1. The vertex sets and the full edge list of H_l for l = 4,5,6,7, built from T1,
     and the sizes of the two parts (checked = 2l-2 each).
 D2. e(H_l) by counting the edges of the built graph, compared against T2.
 D3. The complete cycle-length spectrum of H_l by exhaustive DFS over all simple
     cycles -- hence in particular whether a C_{2l} exists (claim: it does not).
     The spectrum is also compared with the structural prediction of the proof
     (R-cycles occupy 4,6,...,2r; Q-cycles occupy 2r+4,...,4r; nothing else),
     and with T7 for l = 4,5.
 D4. The 2*(2l-2)^2 - ... host non-edges of K_{2l-2,2l-2}, enumerated, and for
     each one an explicit simple path of length 2l-1 joining its endpoints in H_l
     (found by exhaustive DFS, so a failure to find one is a real refutation).
 D5. Independently of D4, and by a different code path: for every host non-edge e,
     an exhaustive search for a simple cycle of length exactly 2l in H_l + e.
     (Control (c): saturation checked by cycle search, not via the path lemma.)
 D6. The conjectured value l^2+l-3 obtained by evaluating T3 at n1 = n2 = 2l-2,
     the observed deficit conj - e(H_l), compared against T4, and the hypothesis
     test 2l-2 >= l+2.
 D7. That sigma (T5) is an automorphism of H_l and reverses Q.
 D8. Calibration (a): the spectra above, which prove the cycle finder does find
     cycles and that only 2l is missing.
 D9. Calibration (b): K_{4,4} has a C_8 and has no C_10 (positive control).
 D10. Calibration (d): G^l_{n1,n2} of T6 built for (l,n1,n2) =
     (4,6,6), (4,7,6), (5,7,7); its edge count derived by counting and compared
     with the Xu-He-Lu formula T3; verified C_{2l}-free and C_{2l}-saturated by
     the same routines.  At the identical instance K_{6,6}/C_8 the program thus
     certifies both their 17-edge graph and the paper's 15-edge graph as
     C_8-saturated, with 15 < 17.
 D11. Lemma 2 of the paper (the engine of the general-l argument), decided inside
     R for every index in its two ranges and in both blocks, together with the
     TIGHTNESS of those ranges (no a_0-b_0 path of length 2r+1, no a-a_0 path of
     length 2r).  The tightness half is the only place where find_path_exact is
     required to answer "no such path", i.e. the only calibration of its
     refuting direction.
 D12. That the proof's four-case partition of the missing edges really exhausts
     the derived host non-edges, with the class sizes 1, r^2-1, r^2-1, (r-1)^2
     that the structure forces (sum 3r^2-2r = |E(host)| - e(H_l)).
 D13. Calibration (e), a NEGATIVE control on the saturation verdict: H_4 with the
     Q-edge x_r b_0 deleted is still C_8-free but is NOT C_8-saturated, and both
     saturation routes must report a non-empty and IDENTICAL failure set.
     Without it, "zero unwitnessed non-edges" is never shown to be a falsifiable
     outcome, because every other graph fed to the saturation test is expected to
     be saturated.  Everywhere else the two routes' failure sets are compared
     non-edge by non-edge as well, and not merely both tested for emptiness,
     which would compare them only at 0 = 0.

Three of the recorded checks are arithmetic identities that cannot fail and are
labelled as such in their own output lines: the evaluation of n1+n2+l^2-3l+1 at
n1=n2=2l-2 to l^2+l-3, and the hypothesis test n_i >= l+2 (equivalent to l >= 4).
They audit the paper's algebra rather than the construction.  eq. (4) is
NOT checked as e == (r^2-1)+(2r+1), which is the identical condition to
e == l^2-1 when r = l-1; its two summands are counted separately instead.

All arithmetic is over Python ints.  There is no floating point anywhere in this
program and no numeric tolerance: every check is an exact combinatorial decision
(set equality, integer equality, existence/non-existence of a subgraph found by
exhaustive search).

Target: Python 3.9, standard library only.
Exit status 0 iff every check passes.
"""

import sys

# ----------------------------------------------------------------------------
# Data taken from the paper, plus the comparison targets fixed here (T1-T7 above).
# ----------------------------------------------------------------------------

# Range of l the paper's theorem is claimed for; we test the first four cases.
L_VALUES = (4, 5, 6, 7)

# T7: reference expected cycle-length spectra (comparison targets only).
EXPECTED_SPECTRUM = {
    4: {4, 6, 10, 12},
    5: {4, 6, 8, 12, 14, 16},
}

# T7 / T6: Xu-He-Lu Definition 2.2 instances for calibration control (d),
# with their edge counts, fixed in advance here, as comparison targets.
XHL_INSTANCES = (
    (4, 6, 6, 17),
    (4, 7, 6, 18),
    (5, 7, 7, 25),
)

# T7: number of non-edges of the host K_{6,6} relative to H_4.
EXPECTED_NONEDGES_L4 = 21


class Report(object):
    """Accumulates PASS/FAIL lines in the fixed output format."""

    def __init__(self):
        self.results = []

    def check(self, ok, label):
        self.results.append((bool(ok), label))
        print(("PASS " if ok else "FAIL ") + label)
        return bool(ok)

    def note(self, text):
        print("      " + text)

    def verdict(self):
        n = len(self.results)
        bad = sum(1 for ok, _ in self.results if not ok)
        print("")
        if bad == 0:
            print("VERDICT: ALL %d CHECKS PASS" % n)
        else:
            print("VERDICT: %d OF %d CHECKS FAILED" % (bad, n))
        return 0 if bad == 0 else 1


# ----------------------------------------------------------------------------
# Graph primitives.  A graph is a dict  vertex -> set(neighbours).  Vertices are
# strings, so nothing depends on an accidental numbering.
# ----------------------------------------------------------------------------

def new_graph(vertices):
    """Edgeless graph on the given vertices."""
    return dict((v, set()) for v in vertices)


def add_edge(g, u, v):
    """Add the undirected edge uv; both endpoints must already be vertices."""
    if u not in g or v not in g:
        raise KeyError("endpoint outside vertex set: %r %r" % (u, v))
    if u == v:
        raise ValueError("loop requested: %r" % (u,))
    g[u].add(v)
    g[v].add(u)


def edge_list(g):
    """Sorted list of edges as (u, v) with u < v.  Used to count and to print."""
    out = []
    for u in g:
        for v in g[u]:
            if u < v:
                out.append((u, v))
    out.sort()
    return out


def edge_count(g):
    return len(edge_list(g))


def complete_bipartite(part_a, part_b):
    """K_{|A|,|B|} on the given (disjoint) parts."""
    g = new_graph(list(part_a) + list(part_b))
    for u in part_a:
        for v in part_b:
            add_edge(g, u, v)
    return g


def host_non_edges(g, part_a, part_b):
    """Edges of the host complete bipartite graph on (A,B) that are absent from g.

    Returned sorted, as (u, v) with u in part_a, v in part_b.
    """
    out = []
    for u in sorted(part_a):
        for v in sorted(part_b):
            if v not in g[u]:
                out.append((u, v))
    return out


# ----------------------------------------------------------------------------
# The paper's construction H_l  (T1).  Nothing about it is asserted here; the
# function only transcribes the definition and returns the object.
# ----------------------------------------------------------------------------

def build_H(l):
    """Build H_l exactly as the paper's Section 2 defines it.

    Returns a dict with keys:
       'g'   the graph, 'A'/'B' the two partite classes (lists),
       'A0','B0','X','Y' the four blocks, 'Q' the vertex sequence of the path Q,
       'R_edges' the edge list of R = K_{r,r} - a_0 b_0, 'r' = l-1.
    """
    r = l - 1
    A0 = ["a%d" % k for k in range(r)]          # a_0 .. a_{r-1}
    B0 = ["b%d" % k for k in range(r)]          # b_0 .. b_{r-1}
    X = ["x%d" % i for i in range(1, r + 1)]    # x_1 .. x_r
    Y = ["y%d" % j for j in range(1, r + 1)]    # y_1 .. y_r
    A = A0 + X
    B = B0 + Y
    g = new_graph(A + B)

    # R = K_{r,r} - a_0 b_0  on A_0 u B_0.
    R_edges = []
    for a in A0:
        for b in B0:
            if a == "a0" and b == "b0":
                continue
            add_edge(g, a, b)
            R_edges.append((a, b))

    # Q: a_0 - y_1 - x_1 - y_2 - x_2 - ... - y_r - x_r - b_0
    Q = ["a0"]
    for i in range(1, r + 1):
        Q.append("y%d" % i)
        Q.append("x%d" % i)
    Q.append("b0")
    for s in range(len(Q) - 1):
        add_edge(g, Q[s], Q[s + 1])

    return {"g": g, "A": A, "B": B, "A0": A0, "B0": B0, "X": X, "Y": Y,
            "Q": Q, "R_edges": R_edges, "r": r, "l": l}


# ----------------------------------------------------------------------------
# Xu-He-Lu's own Definition 2.2 construction (T6): the graph that realises their
# proven Theorem 1.1 upper bound.  Used as calibration control (d).
# ----------------------------------------------------------------------------

def build_XHL(l, n1, n2):
    """Build G^l_{n1,n2} of Xu-He-Lu Definition 2.2.

    Requires l >= 3 and n1,n2 >= l+2 (their hypothesis).  Returns a dict with
    'g', 'V1', 'V2' and the blocks A1,B1,A2,B2 and the special vertices.
    """
    if l < 3 or n1 < l + 2 or n2 < l + 2:
        raise ValueError("outside Xu-He-Lu hypothesis: l=%d n=(%d,%d)" % (l, n1, n2))
    V1 = ["p%02d" % i for i in range(n1)]
    V2 = ["q%02d" % j for j in range(n2)]
    A1, B1 = V1[:n1 - l], V1[n1 - l:]           # |A1| = n1-l, |B1| = l
    A2, B2 = V2[:n2 - l], V2[n2 - l:]           # |A2| = n2-l, |B2| = l
    x1, x2 = B1[0], B1[1]                       # distinct, in B1
    y1 = B2[0]                                  # in B2
    B2p = B2[1:]                                # B_2' = B_2 \ {y_1}, size l-1

    g = new_graph(V1 + V2)
    for u in A1:                                # (1) u y_1 for u in A_1
        add_edge(g, u, y1)
    for v in A2:                                # (2) v x_1 for v in A_2
        add_edge(g, v, x1)
    add_edge(g, y1, x2)                         # (3) y_1 x_2
    for u in B1:                                # (4) complete B_1 x B_2'
        for v in B2p:
            add_edge(g, u, v)

    return {"g": g, "V1": V1, "V2": V2, "A1": A1, "B1": B1, "A2": A2,
            "B2": B2, "B2p": B2p, "x1": x1, "x2": x2, "y1": y1,
            "l": l, "n1": n1, "n2": n2}


# ----------------------------------------------------------------------------
# Breadth-first distances, used only as a SOUND lower bound for pruning: a
# distance measured in the whole graph is never larger than the distance in the
# subgraph a partial DFS is still allowed to use, so pruning on it can never
# discard a real solution.  Unreachable vertices get a distance of len(g)+1.
# ----------------------------------------------------------------------------

def search_order(g):
    """Adjacency lists sorted by (degree, name).

    This is a pure search-order heuristic and changes no answer: the DFS below is
    exhaustive either way.  Trying low-degree neighbours first makes the searches
    walk into the degree-two path Q early, which is where the witnesses live.
    """
    return dict((v, sorted(g[v], key=lambda w: (len(g[w]), w))) for v in g)


def bfs_dist(g, src):
    inf = len(g) + 1
    dist = dict((v, inf) for v in g)
    dist[src] = 0
    frontier = [src]
    while frontier:
        nxt = []
        for v in frontier:
            dv = dist[v] + 1
            for w in g[v]:
                if dist[w] > dv:
                    dist[w] = dv
                    nxt.append(w)
        frontier = nxt
    return dist


# ----------------------------------------------------------------------------
# Exhaustive search for a simple path of an exact length (in edges).
# Exhaustive means: if it returns None, no such path exists.
# ----------------------------------------------------------------------------

def find_path_exact(g, u, v, length):
    """Return a simple u--v path with exactly `length` edges, or None.

    `length` counts edges, matching the paper's convention ("All path lengths
    below count edges").  The search is a complete DFS with two sound prunings:
    a vertex is never revisited, and a branch dies as soon as the edge budget is
    exhausted or the target is reached too early (it could not be re-entered).
    """
    if u == v:
        raise ValueError("path endpoints coincide: %r" % (u,))
    if length < 1:
        return None
    path = [u]
    seen = set([u])
    dv = bfs_dist(g, v)          # lower bound on edges still needed to reach v
    adj = search_order(g)

    def dfs(cur, used):
        if cur == v:
            return used == length
        if used == length:
            return False
        if dv[cur] > length - used:
            return False
        for w in adj[cur]:
            if w in seen:
                continue
            # w == v is only admissible on the very last edge
            if w == v and used + 1 != length:
                continue
            seen.add(w)
            path.append(w)
            if dfs(w, used + 1):
                return True
            path.pop()
            seen.discard(w)
        return False

    if dfs(u, 0):
        return list(path)
    return None


# ----------------------------------------------------------------------------
# Exhaustive enumeration of ALL simple cycles: gives the complete cycle-length
# spectrum.  Each cycle is generated from its smallest vertex (in the fixed
# ordering), so no cycle is missed and none is attributed to a wrong length.
# ----------------------------------------------------------------------------

def cycle_spectrum(g):
    """Return (spectrum, counts, examples).

    spectrum : set of the lengths (in edges = vertices) of the simple cycles of g
    counts   : dict length -> number of simple cycles of that length
    examples : dict length -> one witnessing cycle, as a vertex list
    """
    order = sorted(g)
    rank = dict((v, i) for i, v in enumerate(order))
    adj = dict((v, sorted(g[v])) for v in g)
    counts = {}
    examples = {}

    for s in order:
        rs = rank[s]
        path = [s]
        seen = set([s])

        def dfs(cur):
            for w in adj[cur]:
                if rank[w] < rs:
                    continue                     # cycle owned by a smaller start
                if w == s:
                    if len(path) >= 3:           # a genuine cycle
                        k = len(path)
                        counts[k] = counts.get(k, 0) + 1
                        if k not in examples:
                            examples[k] = list(path)
                    continue
                if w in seen:
                    continue
                seen.add(w)
                path.append(w)
                dfs(w)
                path.pop()
                seen.discard(w)

        dfs(s)

    # every cycle was traversed once in each direction
    counts = dict((k, n // 2) for k, n in counts.items())
    return set(counts), counts, examples


# ----------------------------------------------------------------------------
# Exhaustive search for a simple cycle of one exact length.  This is a second,
# independent implementation (depth-bounded, no enumeration of other lengths);
# it is the routine used for the direct saturation test, control (c).
# ----------------------------------------------------------------------------

def find_cycle_exact(g, k, must_contain=None):
    """Return a simple cycle on exactly k vertices (= k edges), or None.

    Complete search.  With must_contain=None every cycle is attempted from its
    smallest vertex, so None means no C_k exists in g at all.  With
    must_contain=w the single start w is used and no vertex is excluded, so the
    search is complete for the question "is there a C_k through w?" -- None then
    means no C_k of g passes through w.
    """
    if k < 3:
        return None
    order = sorted(g)
    rank = dict((v, i) for i, v in enumerate(order))
    adj = search_order(g)
    if must_contain is None:
        starts, restrict = order, True
    else:
        starts, restrict = [must_contain], False

    for s in starts:
        rs = rank[s]
        path = [s]
        seen = set([s])
        ds = bfs_dist(g, s)      # lower bound on edges still needed to close up

        def dfs(cur):
            if len(path) == k:
                return s in g[cur]
            if ds[cur] > k - len(path) + 1:
                return False
            for w in adj[cur]:
                if w in seen or (restrict and rank[w] <= rs):
                    continue
                seen.add(w)
                path.append(w)
                if dfs(w):
                    return True
                path.pop()
                seen.discard(w)
            return False

        if dfs(s):
            return list(path)
    return None


# ----------------------------------------------------------------------------
# Independent validators for the objects the searches return, so that a bug in a
# search routine cannot manufacture a pass.
# ----------------------------------------------------------------------------

def validate_path(g, path, u, v, length):
    """True iff `path` really is a simple u--v path of g with `length` edges."""
    if path is None:
        return False
    if len(path) != length + 1:
        return False
    if path[0] != u or path[-1] != v:
        return False
    if len(set(path)) != len(path):
        return False
    for i in range(len(path) - 1):
        if path[i + 1] not in g.get(path[i], ()):
            return False
    return True


def validate_cycle(g, cyc, k, must_use=None):
    """True iff `cyc` really is a simple cycle of g on k vertices.

    If `must_use` is a pair (u,v), the cycle must also traverse that edge.
    """
    if cyc is None:
        return False
    if len(cyc) != k or len(set(cyc)) != k:
        return False
    for i in range(k):
        a, b = cyc[i], cyc[(i + 1) % k]
        if b not in g.get(a, ()):
            return False
    if must_use is not None:
        u, v = must_use
        pairs = set()
        for i in range(k):
            a, b = cyc[i], cyc[(i + 1) % k]
            pairs.add((a, b))
            pairs.add((b, a))
        if (u, v) not in pairs:
            return False
    return True


# ----------------------------------------------------------------------------
# Saturation, route 1: the paper's own route.  Every host non-edge must be joined
# in H by a simple path of length 2l-1, which the missing edge then closes to a
# C_{2l}.
# ----------------------------------------------------------------------------

def saturation_via_paths(g, part_a, part_b, l):
    """Return (non_edges, n_witnessed, failures, sample).

    failures lists the non-edges for which no simple path of length 2l-1 was
    found, or for which the path returned failed independent validation.
    sample is one (non_edge, path) pair for the transcript.
    """
    want = 2 * l - 1
    nes = host_non_edges(g, part_a, part_b)
    failures = []
    sample = None
    n_ok = 0
    for (u, v) in nes:
        p = find_path_exact(g, u, v, want)
        if validate_path(g, p, u, v, want):
            n_ok += 1
            if sample is None:
                sample = ((u, v), p)
        else:
            failures.append((u, v))
    return nes, n_ok, failures, sample


# ----------------------------------------------------------------------------
# Saturation, route 2 (control (c)): forget the path lemma.  Add the missing edge
# and search directly for a cycle of length exactly 2l in H + e.
# ----------------------------------------------------------------------------

def saturation_via_cycles(g, part_a, part_b, l):
    """Return (n_non_edges, n_witnessed, failures, sample, witnesses).

    `witnesses` maps each non-edge to the C_{2l} that was found AND passed
    independent validation (including that it traverses the new edge), so the
    printed transcript is certified rather than merely printed.
    """
    k = 2 * l
    nes = host_non_edges(g, part_a, part_b)
    failures = []
    sample = None
    witnesses = {}
    n_ok = 0
    for (u, v) in nes:
        g2 = dict((w, set(nb)) for w, nb in g.items())
        add_edge(g2, u, v)
        # First look for a C_k through u (the cycle must contain the new edge, so
        # it must contain u); if that fails, fall back to the unrestricted search
        # so that no genuine witness can be missed.
        cyc = find_cycle_exact(g2, k, must_contain=u)
        if not validate_cycle(g2, cyc, k, must_use=(u, v)):
            cyc = find_cycle_exact(g2, k)
        if validate_cycle(g2, cyc, k, must_use=(u, v)):
            n_ok += 1
            witnesses[(u, v)] = cyc
            if sample is None:
                sample = ((u, v), cyc)
        else:
            failures.append((u, v))
    return len(nes), n_ok, failures, sample, witnesses


# ----------------------------------------------------------------------------
# The symmetry the proof uses to halve its case analysis (T5).  We do not take
# its automorphism property on trust: we build the map and test it.
# ----------------------------------------------------------------------------

def sigma_of(H):
    """The map  a_0<->b_0, a_k<->b_k, x_i<->y_{r-i+1}  as a dict."""
    r = H["r"]
    s = {}
    for k in range(r):
        s["a%d" % k] = "b%d" % k
        s["b%d" % k] = "a%d" % k
    for i in range(1, r + 1):
        s["x%d" % i] = "y%d" % (r - i + 1)
        s["y%d" % i] = "x%d" % (r - i + 1)
    return s


def sigma_facts(H):
    """Return (is_bijection, is_involution, preserves_edges, reverses_Q)."""
    g, Q = H["g"], H["Q"]
    s = sigma_of(H)
    bij = (set(s) == set(g)) and (set(s.values()) == set(g))
    inv = all(s[s[v]] == v for v in s) if bij else False
    pres = True
    if bij:
        for (u, v) in edge_list(g):
            if s[v] not in g[s[u]]:
                pres = False
                break
        # and the image edge set must have the same size (it does, s injective)
    rev = bij and all(s[Q[t]] == Q[len(Q) - 1 - t] for t in range(len(Q)))
    return bij, inv, pres, rev


# ----------------------------------------------------------------------------
# What the proof's structural argument predicts the cycle spectrum to be.
# Cycles inside R = K_{r,r} - a_0b_0 realise every even length 4..2r (r >= 3).
# A cycle meeting an internal vertex of Q must traverse all of Q (2r+1 edges)
# and come back through an a_0--b_0 path of R of odd length 2s+1, 1 <= s <= r-1,
# giving every even length 2r+4 .. 4r.  Nothing else can occur.  So the only even
# length in [4,4r] that must be absent is 2r+2 = 2l.
# ----------------------------------------------------------------------------

def bipartite_ok(g, part_a, part_b):
    """True iff (part_a, part_b) is a partition of V(g) and every edge crosses it.

    This is what makes g a spanning subgraph of the host K_{|A|,|B|}: without it
    every later count would be about some other graph.
    """
    A, B = set(part_a), set(part_b)
    if len(A) != len(part_a) or len(B) != len(part_b):
        return False                      # repeated vertex names
    if A & B:
        return False
    if A | B != set(g):
        return False
    for (u, v) in edge_list(g):
        if not ((u in A and v in B) or (u in B and v in A)):
            return False
    return True


def transcription_facts(H):
    """Check that the built H_l really is R = K_{r,r}-a_0b_0 plus the path Q.

    Returns (R_is_right, Q_is_a_path, edges_partition).
    """
    g, r, Q = H["g"], H["r"], H["Q"]
    A0, B0 = H["A0"], H["B0"]
    R = set()
    for a in A0:
        for b in B0:
            if not (a == "a0" and b == "b0"):
                R.add((min(a, b), max(a, b)))
    R_ok = (len(R) == r * r - 1) and ("b0" not in g["a0"])
    Qe = set()
    for t in range(len(Q) - 1):
        Qe.add((min(Q[t], Q[t + 1]), max(Q[t], Q[t + 1])))
    Q_ok = (len(set(Q)) == 2 * r + 2 and len(Qe) == 2 * r + 1
            and all(Q[t + 1] in g[Q[t]] for t in range(len(Q) - 1)))
    part_ok = (not (R & Qe)) and (R | Qe) == set(edge_list(g))
    return R_ok, Q_ok, part_ok


def predicted_spectrum(l):
    r = l - 1
    inside_R = set(range(4, 2 * r + 1, 2))
    through_Q = set(range(2 * r + 4, 4 * r + 1, 2))
    return inside_R, through_Q, inside_R | through_Q


def core_graph(H):
    """R = K_{r,r} - a_0b_0 on its own, as a standalone graph."""
    R = new_graph(H["A0"] + H["B0"])
    for (a, b) in H["R_edges"]:
        add_edge(R, a, b)
    return R


def lemma_core_facts(H):
    """Decide Lemma 2 of the paper -- the engine of the general-l argument.

    (1) R has an a_0--b_0 path of every odd length 2s+1, 1 <= s <= r-1;
    (2) for every a in A_0\\{a_0}, R has an a--a_0 path of every even length 2s,
        1 <= s <= r-1;  and the mirror statements in B_0.
    Also the TIGHTNESS of the two ranges (2r+1 and 2r are unreachable), which is
    the only place in this program where find_path_exact is required to answer
    "no such path" -- i.e. the only calibration of its refuting direction.

    Returns (n_wanted, missing, tight_ok, sample) with `missing` the list of
    (kind, endpoints, length) triples for which no validated path was found.
    """
    R = core_graph(H)
    r = H["r"]
    missing = []
    n_wanted = 0
    sample = None
    for s in range(1, r):
        for (u, v) in (("a0", "b0"), ("b0", "a0")):
            n_wanted += 1
            p = find_path_exact(R, u, v, 2 * s + 1)
            if validate_path(R, p, u, v, 2 * s + 1):
                if sample is None:
                    sample = ("L1(1)", u, v, 2 * s + 1, p)
            else:
                missing.append(("L1(1)", (u, v), 2 * s + 1))
    for s in range(1, r):
        for (block, root) in ((H["A0"], "a0"), (H["B0"], "b0")):
            for a in block:
                if a == root:
                    continue
                n_wanted += 1
                p = find_path_exact(R, a, root, 2 * s)
                if not validate_path(R, p, a, root, 2 * s):
                    missing.append(("L1(2)", (a, root), 2 * s))
    # tightness: only 2r vertices are available in R
    tight_ok = (find_path_exact(R, "a0", "b0", 2 * r + 1) is None
                and find_path_exact(R, "a1", "a0", 2 * r) is None)
    return n_wanted, missing, tight_ok, sample


def nonedge_classes(H):
    """Split the derived host non-edges into the proof's four cases.

    The proof handles: a_0b_0; then a y_j with a in A_0; then X-B_0 by the sigma
    symmetry; then x_i y_j.  It claims "these cases exhaust all missing edges".
    We classify the DERIVED non-edge list and return the observed sizes together
    with the closed forms the structure forces:
        A_0-B_0 : 1                 (only a_0b_0 is absent)
        A_0-Y   : r^2 - 1           (only a_0y_1 is present)
        X-B_0   : r^2 - 1           (only x_r b_0 is present)
        X-Y     : r^2 - (2r-1)      (Q uses x_iy_i and x_iy_{i+1})
    Their sum must be |E(K_{2r,2r})| - e(H_l) = 3r^2 - 2r.
    """
    r = H["r"]
    A0, B0, X, Y = set(H["A0"]), set(H["B0"]), set(H["X"]), set(H["Y"])
    nes = host_non_edges(H["g"], H["A"], H["B"])
    obs = {"A0-B0": 0, "A0-Y": 0, "X-B0": 0, "X-Y": 0}
    unclassified = []
    for (u, v) in nes:
        if u in A0 and v in B0:
            obs["A0-B0"] += 1
        elif u in A0 and v in Y:
            obs["A0-Y"] += 1
        elif u in X and v in B0:
            obs["X-B0"] += 1
        elif u in X and v in Y:
            obs["X-Y"] += 1
        else:
            unclassified.append((u, v))
    want = {"A0-B0": 1, "A0-Y": r * r - 1, "X-B0": r * r - 1,
            "X-Y": r * r - (2 * r - 1)}
    total_ok = (sum(obs.values()) == len(nes) == 3 * r * r - 2 * r)
    return obs, want, unclassified, total_ok


# ----------------------------------------------------------------------------
# Per-l block 1: sizes, edge count, the conjecture arithmetic, the symmetry.
# ----------------------------------------------------------------------------

def check_structure(rep, H):
    l, r, g = H["l"], H["r"], H["g"]
    tag = "[l=%d]" % l
    nA, nB = len(H["A"]), len(H["B"])
    e = edge_count(g)

    print("")
    print("=== paper's H_%d  (r = l-1 = %d) ===" % (l, r))
    print("  A = %s" % (", ".join(H["A"])))
    print("  B = %s" % (", ".join(H["B"])))
    print("  Q : %s" % (" - ".join(H["Q"])))
    print("  derived |A| = %d, |B| = %d, e(H_%d) = %d  (edges: %s)"
          % (nA, nB, l, e,
             ", ".join("%s-%s" % (u, v) for (u, v) in edge_list(g))))

    rep.check(nA == 2 * l - 2 and nB == 2 * l - 2,
              "%s both parts of H_l have size 2l-2 = %d (derived %d, %d)"
              % (tag, 2 * l - 2, nA, nB))
    rep.check(bipartite_ok(g, H["A"], H["B"]),
              "%s A = A_0 u X and B = B_0 u Y partition V(H_l) and every edge "
              "crosses them, so H_l is a spanning subgraph of K_{%d,%d}"
              % (tag, 2 * l - 2, 2 * l - 2))
    R_ok, Q_ok, part_ok = transcription_facts(H)
    rep.check(R_ok and Q_ok and part_ok,
              "%s the built graph is exactly R = K_{r,r} - a_0b_0 (r^2-1 = %d "
              "edges, a_0b_0 absent) disjointly united with the %d-edge path Q "
              "[R=%s Q=%s partition=%s]"
              % (tag, r * r - 1, 2 * r + 1, R_ok, Q_ok, part_ok))
    rep.check(e == l * l - 1,
              "%s e(H_l) counted = %d equals the paper's l^2-1 = %d"
              % (tag, e, l * l - 1))
    # eq. (4): the two SUMMANDS are counted separately from the derived
    # edge list.  (Checking only e == (r^2-1)+(2r+1) would be the identical
    # condition to e == l^2-1, since r = l-1 makes r^2+2r = (r+1)^2-1; that
    # duplicate cannot fail, so it is not used as the test.)
    core = set(H["A0"]) | set(H["B0"])
    nR = sum(1 for (u, v) in edge_list(g) if u in core and v in core)
    nQ = e - nR
    rep.check(nR == r * r - 1 and nQ == 2 * r + 1 and nR + nQ == e,
              "%s eq. (4) split verified summand-by-summand on the "
              "derived edge list: %d edges inside A_0uB_0 (= r^2-1 = %d) plus %d "
              "edges meeting XuY (= 2r+1 = %d), total %d"
              % (tag, nR, r * r - 1, nQ, 2 * r + 1, e))

    # n1, n2 are the DERIVED part sizes, not the intended 2l-2.
    n1, n2 = nA, nB
    conj = n1 + n2 + l * l - 3 * l + 1
    rep.check(conj == l * l + l - 3,
              "%s [arithmetic identity, cannot fail: it audits the paper's algebra"
              ", not the construction] Xu-He-Lu formula n1+n2+l^2-3l+1 at the "
              "derived n1=%d, n2=%d evaluates to %d = the paper's l^2+l-3 = %d"
              % (tag, n1, n2, conj, l * l + l - 3))
    rep.check(conj - e == l - 2 and e < conj,
              "%s deficit conjectured-minus-constructed = %d - %d = %d equals "
              "l-2 = %d and is positive" % (tag, conj, e, conj - e, l - 2))
    rep.check(n1 >= l + 2 and n2 >= l + 2,
              "%s [equivalent to l>=4, so it cannot fail for any l this program "
              "runs; recorded because the refutation is void outside the "
              "hypothesis] derived host K_{%d,%d} satisfies n_i >= l+2 = %d"
              % (tag, n1, n2, l + 2))

    bij, inv, pres, rev = sigma_facts(H)
    rep.check(bij and inv and pres and rev,
              "%s sigma (a_k<->b_k, x_i<->y_{r-i+1}) is an involutive "
              "automorphism of H_l and reverses Q  [bijection=%s involution=%s "
              "edge-preserving=%s reverses-Q=%s]" % (tag, bij, inv, pres, rev))

    # Lemma 2 of the paper -- the engine that makes the argument work for EVERY
    # l.  Both index ranges, both blocks, plus the tightness of the ranges.
    n_want, miss, tight, samp = lemma_core_facts(H)
    rep.check(n_want > 0 and not miss,
              "%s Lemma 2 decided in R for all %d (kind, endpoints, length) "
              "instances -- odd a_0-b_0 lengths 2s+1 and even a-a_0 lengths 2s "
              "for 1<=s<=r-1=%d, in A_0 and in B_0; unrealised: %d %s"
              % (tag, n_want, r - 1, len(miss), miss[:4]))
    if samp is not None:
        print("    Lemma 2(1) sample: %s-%s of length %d : %s"
              % (samp[1], samp[2], samp[3], "-".join(samp[4])))
    rep.check(tight,
              "%s Lemma 2's index ranges are tight: R has NO a_0-b_0 path of "
              "length 2r+1 = %d and NO a_1-a_0 path of length 2r = %d (only 2r "
              "vertices exist).  This is the only test of find_path_exact in its "
              "refuting direction." % (tag, 2 * r + 1, 2 * r))

    # The proof's case partition of the missing edges: does it exhaust them?
    obs, want, unclass, tot = nonedge_classes(H)
    keys = ("A0-B0", "A0-Y", "X-B0", "X-Y")
    rep.check(all(obs[k] == want[k] for k in keys) and not unclass and tot,
              "%s the proof's four non-edge cases exhaust the derived host "
              "non-edges, with the sizes the structure forces: %s vs expected %s"
              " (unclassified %d, total %d = 3r^2-2r = %d)"
              % (tag, ", ".join("%s:%d" % (k, obs[k]) for k in keys),
                 ", ".join("%s:%d" % (k, want[k]) for k in keys),
                 len(unclass), sum(obs.values()), 3 * r * r - 2 * r))
    return e, conj


# ----------------------------------------------------------------------------
# Per-l block 2: the cycle spectrum -- this is where C_{2l}-freeness is decided,
# and where calibration (a) lives (the finder demonstrably DOES find cycles).
# ----------------------------------------------------------------------------

def check_cycles(rep, H):
    l, r, g = H["l"], H["r"], H["g"]
    tag = "[l=%d]" % l
    spec, counts, examples = cycle_spectrum(g)
    inR, thruQ, pred = predicted_spectrum(l)

    print("  derived cycle-length spectrum of H_%d : %s" % (l, sorted(spec)))
    print("  cycle counts by length : %s"
          % ", ".join("C_%d:%d" % (k, counts[k]) for k in sorted(counts)))
    for k in sorted(examples):
        print("    example C_%-2d : %s" % (k, "-".join(examples[k])))
    bad_ex = [k for k in sorted(examples) if not validate_cycle(g, examples[k], k)]
    rep.check(len(examples) > 0 and not bad_ex,
              "%s all %d exhibited example cycles (one per realised length) pass "
              "independent validation; bad: %s" % (tag, len(examples), bad_ex))

    rep.check(2 * l not in spec,
              "%s H_l is C_{2l}-free: length %d absent from the derived "
              "spectrum %s" % (tag, 2 * l, sorted(spec)))
    rep.check(find_cycle_exact(g, 2 * l) is None,
              "%s independent depth-bounded search finds no C_%d in H_l"
              % (tag, 2 * l))
    rep.check(len(spec) > 0 and max(spec) >= 4 * r,
              "%s calibration (a): the cycle finder does find cycles -- %d "
              "cycles in %d distinct lengths, longest %d = 4r"
              % (tag, sum(counts.values()), len(spec), max(spec) if spec else 0))
    rep.check(spec == pred,
              "%s spectrum equals the proof's prediction: R-cycles %s (cap 2r=%d)"
              " union Q-cycles %s (start 2r+4=%d)"
              % (tag, sorted(inR), 2 * r, sorted(thruQ), 2 * r + 4))
    missing = [k for k in range(4, 4 * r + 1, 2) if k not in spec]
    rep.check(missing == [2 * l],
              "%s the only even length in [4,4r]=[4,%d] absent from H_l is "
              "2l = %d (derived missing set %s)" % (tag, 4 * r, 2 * l, missing))
    if l in EXPECTED_SPECTRUM:
        rep.check(spec == EXPECTED_SPECTRUM[l],
                  "%s derived spectrum %s equals the comparison target %s fixed "
                  "in advance here (the paper does not print it)"
                  % (tag, sorted(spec), sorted(EXPECTED_SPECTRUM[l])))
    return spec


# ----------------------------------------------------------------------------
# Per-l block 3: saturation, by two independent routes.
# ----------------------------------------------------------------------------

def check_saturation(rep, H, e):
    l, g = H["l"], H["g"]
    tag = "[l=%d]" % l
    host = len(H["A"]) * len(H["B"])     # derived host edge count, not 2l-2 twice

    nes, n_path, fail_path, samp_p = saturation_via_paths(g, H["A"], H["B"], l)
    print("  host K_{%d,%d} has %d edges; H_l has %d; derived non-edges: %d"
          % (2 * l - 2, 2 * l - 2, host, e, len(nes)))
    rep.check(len(nes) == host - e,
              "%s non-edge count %d equals |E(host)| - e(H_l) = %d - %d"
              % (tag, len(nes), host, e))
    if l == 4:
        rep.check(len(nes) == EXPECTED_NONEDGES_L4,
                  "%s derived non-edge count %d equals the comparison target 21 "
                  "fixed in advance here" % (tag, len(nes)))
    rep.check(not fail_path and n_path == len(nes),
              "%s route 1 (path lemma): all %d non-edges are joined in H_l by a "
              "validated simple path of length 2l-1 = %d; unwitnessed: %d %s"
              % (tag, len(nes), 2 * l - 1, len(fail_path), fail_path[:6]))
    if samp_p is not None:
        print("    e.g. non-edge %s-%s witnessed by %s"
              % (samp_p[0][0], samp_p[0][1], "-".join(samp_p[1])))

    n_ne, n_cyc, fail_cyc, samp_c, wit = saturation_via_cycles(g, H["A"], H["B"], l)
    rep.check(not fail_cyc and n_cyc == n_ne,
              "%s route 2 (control (c), direct cycle search): H_l + e contains a "
              "validated C_%d through e for all %d non-edges e; failures: %d %s"
              % (tag, 2 * l, n_ne, len(fail_cyc), fail_cyc[:6]))
    if samp_c is not None:
        print("    e.g. adding %s-%s creates the C_%d %s"
              % (samp_c[0][0], samp_c[0][1], 2 * l, "-".join(samp_c[1])))
    # The two routes decide logically equivalent predicates, so they must agree
    # NON-EDGE BY NON-EDGE, not merely both report "no failures".  On a saturated
    # graph both sets are empty; the same comparison is exercised on a genuinely
    # non-saturated graph in calibration control (e).
    rep.check(set(fail_path) == set(fail_cyc),
              "%s routes 1 and 2 agree non-edge by non-edge: identical failure "
              "sets (sizes %d and %d)"
              % (tag, len(fail_path), len(fail_cyc)))

    if l == 4:
        print("    full witness list for l=4 (non-edge : validated closing C_8, "
              "each certified to traverse the added edge):")
        for (u, v) in nes:
            print("      %-8s : %s" % ("%s-%s" % (u, v),
                                       "-".join(wit[(u, v)])
                                       if (u, v) in wit else "NONE"))
        rep.check(len(wit) == len(nes),
                  "%s every one of the %d printed l=4 witness lines is a "
                  "validated C_8 through the added edge (%d of %d)"
                  % (tag, len(nes), len(wit), len(nes)))
    return len(nes)


# ----------------------------------------------------------------------------
# Calibration (b): a positive control on the cycle finder.  K_{4,4} must be found
# to contain a C_8 (it is Hamiltonian) and not to contain a C_10 (too few
# vertices).  A finder that always answers "free" fails the first of these.
# ----------------------------------------------------------------------------

def check_positive_control(rep):
    A = ["u%d" % i for i in range(1, 5)]
    B = ["w%d" % i for i in range(1, 5)]
    k44 = complete_bipartite(A, B)
    spec, counts, examples = cycle_spectrum(k44)
    print("")
    print("=== calibration (b): positive control on K_{4,4} ===")
    print("  derived e(K_{4,4}) = %d, spectrum %s, counts %s"
          % (edge_count(k44), sorted(spec),
             ", ".join("C_%d:%d" % (k, counts[k]) for k in sorted(counts))))
    c8 = find_cycle_exact(k44, 8)
    c10 = find_cycle_exact(k44, 10)
    rep.check(edge_count(k44) == 16,
              "[control b] e(K_{4,4}) derived = %d = 16" % edge_count(k44))
    rep.check(validate_cycle(k44, c8, 8) and 8 in spec,
              "[control b] K_{4,4} HAS a C_8, exhibited and validated: %s"
              % ("-".join(c8) if c8 else "NONE"))
    rep.check(c10 is None and 10 not in spec,
              "[control b] K_{4,4} has NO C_10 (both routines agree)")
    rep.check(spec == set([4, 6, 8]),
              "[control b] full K_{4,4} spectrum derived as %s = {4,6,8}"
              % sorted(spec))
    # counts are classical: C_4: 36, C_6: 96, C_8: 72 in K_{4,4}
    rep.check(counts.get(4) == 36 and counts.get(6) == 96 and counts.get(8) == 72,
              "[control b] cycle counts in K_{4,4} derived as C_4:%s C_6:%s "
              "C_8:%s, matching C(4,k)^2 k!(k-1)!/2 = 36, 96, 72"
              % (counts.get(4), counts.get(6), counts.get(8)))


# ----------------------------------------------------------------------------
# Calibration (e): a NEGATIVE control on the
# SATURATION verdict.  Controls (a)-(d) all calibrate the cycle finder; without
# (e) nothing in this program ever shows that the saturation test is capable of
# reporting an unwitnessed non-edge -- every graph it is given is expected to be
# saturated, so "failures == []" is never observed to be a falsifiable outcome.
#
# The witness: H_4 with the last edge x_r b_0 of Q deleted.  It is still C_8-free
# (X u Y now hangs off a_0 only, so no cycle can meet it, and R-cycles cap at 6),
# but it is NOT C_8-saturated: a_0-b_0 has no path of length 7 left, because an
# a_0--b_0 path must stay inside R, where the longest one has length 2r-1 = 5.
# Both routes must therefore report a NON-EMPTY failure set, and the two failure
# sets must be equal -- which is where the route-agreement check acquires content.
# ----------------------------------------------------------------------------

def check_negative_control(rep):
    H = build_H(4)
    g, l, r = H["g"], 4, H["r"]
    e_before = edge_count(g)                 # derived, not seeded
    dead = ("b0", "x%d" % r)                 # the edge x_r b_0 of Q
    g[dead[0]].discard(dead[1])
    g[dead[1]].discard(dead[0])
    print("")
    print("=== calibration (e): NEGATIVE control -- H_4 minus the Q edge %s-%s ==="
          % (dead[1], dead[0]))
    spec, counts, examples = cycle_spectrum(g)
    e = edge_count(g)
    print("  derived e = %d (was %d), spectrum %s" % (e, e_before, sorted(spec)))
    rep.check(e == e_before - 1 and 2 * l not in spec
              and find_cycle_exact(g, 2 * l) is None,
              "[control e] exactly one edge was removed (%d -> %d) and the "
              "mutilated graph is still C_8-free (spectrum %s), so a failure "
              "below cannot be blamed on a stray C_8"
              % (e_before, e, sorted(spec)))

    nes, n_path, fail_path, samp = saturation_via_paths(g, H["A"], H["B"], l)
    n_ne, n_cyc, fail_cyc, samp_c, _w = saturation_via_cycles(
        g, H["A"], H["B"], l)
    print("  non-edges %d ; route-1 witnessed %d, unwitnessed %d ; route-2 "
          "witnessed %d, unwitnessed %d"
          % (len(nes), n_path, len(fail_path), n_cyc, len(fail_cyc)))
    print("  unwitnessed non-edges: %s"
          % ", ".join("%s-%s" % uv for uv in fail_path))
    rep.check(len(fail_path) > 0,
              "[control e] the saturation test CAN fail: route 1 reports %d of %d "
              "non-edges unwitnessed on a deliberately non-saturated graph, so "
              "'zero unwitnessed' elsewhere is a falsifiable outcome"
              % (len(fail_path), len(nes)))
    rep.check(("a0", "b0") in set(fail_path),
              "[control e] the specific hand-checked non-edge a_0b_0 is among "
              "them (no a_0-b_0 path of length 7 survives: R's longest is %d)"
              % (2 * r - 1))
    rep.check(len(fail_cyc) > 0 and set(fail_path) == set(fail_cyc),
              "[control e] routes 1 and 2 agree non-edge by non-edge on a "
              "NON-empty failure set (%d = %d), so their agreement above is not "
              "the trivial 0 = 0" % (len(fail_path), len(fail_cyc)))
    rep.check(n_path + len(fail_path) == len(nes)
              and n_cyc + len(fail_cyc) == n_ne == len(nes),
              "[control e] witnessed + unwitnessed = %d non-edges on both routes "
              "(no non-edge is silently dropped by either loop)" % len(nes))


# ----------------------------------------------------------------------------
# Calibration (d), the strongest control: rebuild Xu-He-Lu's OWN Definition 2.2
# graph and confirm that these same routines certify it as C_{2l}-saturated with
# exactly n1+n2+l^2-3l+1 edges -- i.e. that the program reproduces their PROVEN
# Theorem 1.1 bound.  At the identical instance K_{6,6}/C_8 the program therefore
# certifies both their 17-edge graph and the paper's 15-edge graph, which is what
# rules out a verifier that merely flatters the paper.
# ----------------------------------------------------------------------------

def check_xhl_instance(rep, l, n1, n2, expected_e):
    G = build_XHL(l, n1, n2)
    g = G["g"]
    tag = "[control d, l=%d, K_{%d,%d}]" % (l, n1, n2)
    e = edge_count(g)
    formula = n1 + n2 + l * l - 3 * l + 1
    print("")
    print("=== calibration (d): Xu-He-Lu Definition 2.2, l=%d, n=(%d,%d) ==="
          % (l, n1, n2))
    print("  A1=%s B1=%s (x1=%s x2=%s)" % (G["A1"], G["B1"], G["x1"], G["x2"]))
    print("  A2=%s B2=%s (y1=%s, B2'=%s)" % (G["A2"], G["B2"], G["y1"], G["B2p"]))
    print("  derived e = %d ; edges: %s"
          % (e, ", ".join("%s-%s" % uv for uv in edge_list(g))))

    rep.check(e == formula,
              "%s counted e = %d equals their formula n1+n2+l^2-3l+1 = %d"
              % (tag, e, formula))
    rep.check(e == expected_e,
              "%s counted e = %d equals the comparison target %d fixed in "
              "advance here" % (tag, e, expected_e))
    rep.check(bipartite_ok(g, G["V1"], G["V2"]),
              "%s (V1,V2) partitions their vertex set and every edge crosses it, "
              "so it is a spanning subgraph of K_{%d,%d}" % (tag, n1, n2))
    rep.check(len(G["A1"]) == n1 - l and len(G["B1"]) == l
              and len(G["A2"]) == n2 - l and len(G["B2"]) == l
              and len(G["B2p"]) == l - 1,
              "%s block sizes as in their Definition 2.2: |A1|=%d=n1-l, |B1|=%d=l,"
              " |A2|=%d=n2-l, |B2|=%d=l, |B2'|=%d=l-1"
              % (tag, len(G["A1"]), len(G["B1"]), len(G["A2"]), len(G["B2"]),
                 len(G["B2p"])))

    spec, counts, examples = cycle_spectrum(g)
    print("  derived spectrum %s" % sorted(spec))
    rep.check(2 * l not in spec and find_cycle_exact(g, 2 * l) is None,
              "%s their graph is C_%d-free (derived spectrum %s)"
              % (tag, 2 * l, sorted(spec)))
    rep.check(len(spec) > 0,
              "%s the finder still sees cycles here (spectrum non-empty: %s)"
              % (tag, sorted(spec)))

    nes, n_path, fail_path, samp = saturation_via_paths(g, G["V1"], G["V2"], l)
    rep.check(len(nes) == n1 * n2 - e,
              "%s non-edge count %d = n1*n2 - e = %d - %d"
              % (tag, len(nes), n1 * n2, e))
    rep.check(not fail_path and n_path == len(nes),
              "%s every one of the %d non-edges is joined by a validated simple "
              "path of length 2l-1 = %d; unwitnessed: %d %s"
              % (tag, len(nes), 2 * l - 1, len(fail_path), fail_path[:6]))
    n_ne, n_cyc, fail_cyc, samp_c, _wit = saturation_via_cycles(
        g, G["V1"], G["V2"], l)
    rep.check(not fail_cyc and n_cyc == n_ne,
              "%s direct cycle search: their graph + e has a validated C_%d for "
              "all %d non-edges; failures: %d %s"
              % (tag, 2 * l, n_ne, len(fail_cyc), fail_cyc[:6]))
    rep.check(set(fail_path) == set(fail_cyc),
              "%s routes 1 and 2 agree non-edge by non-edge on their graph too "
              "(failure sets of sizes %d and %d)"
              % (tag, len(fail_path), len(fail_cyc)))
    if samp_c is not None:
        print("    e.g. adding %s-%s creates the C_%d %s"
              % (samp_c[0][0], samp_c[0][1], 2 * l, "-".join(samp_c[1])))
    return e


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    rep = Report()
    print("verify.py -- Counterexamples to the Xu-He-Lu Conjecture on")
    print("                 Partite Cycle Saturation")
    print("")
    print("Testing the paper's H_l for l = %s.  Every quantity below is computed"
          % (", ".join(str(l) for l in L_VALUES)))
    print("from the construction; the paper's stated values are used only as")
    print("comparison targets.  All arithmetic is exact integer arithmetic; no")
    print("floating point and no tolerance is used anywhere in this program.")

    # Several calibration comparisons are guarded by "if l == 4" / "if l in
    # EXPECTED_SPECTRUM".  If L_VALUES ever loses those l they would DISAPPEAR
    # silently instead of failing, so their presence is itself a recorded check.
    rep.check(4 in L_VALUES and set(EXPECTED_SPECTRUM) <= set(L_VALUES)
              and all(l >= 4 for l in L_VALUES),
              "[config] L_VALUES = %s contains every l whose calibration targets "
              "are configured (%s) and nothing below l = 4, so no guarded "
              "calibration check can be silently skipped"
              % (list(L_VALUES), sorted(EXPECTED_SPECTRUM)))

    table = []
    for l in L_VALUES:
        H = build_H(l)
        e, conj = check_structure(rep, H)
        check_cycles(rep, H)
        n_ne = check_saturation(rep, H, e)
        table.append((l, e, conj, conj - e, n_ne))

    check_positive_control(rep)
    check_negative_control(rep)

    xhl_e = {}
    for (l, n1, n2, expected_e) in XHL_INSTANCES:
        xhl_e[(l, n1, n2)] = check_xhl_instance(rep, l, n1, n2, expected_e)

    print("")
    print("=== summary ===")
    print("   l   e(H_l)   conjectured l^2+l-3   deficit   host non-edges")
    for (l, e, conj, d, n_ne) in table:
        print("  %2d   %6d   %19d   %7d   %14d" % (l, e, conj, d, n_ne))

    # The refutation itself, stated as a check over everything derived above.
    rep.check(all(e < conj for (_, e, conj, _, _) in table),
              "[refutation] for every l tested, the constructed C_{2l}-saturated "
              "spanning subgraph of K_{2l-2,2l-2} has fewer edges than Xu-He-Lu "
              "Conjecture 6.2 asserts: %s"
              % ", ".join("l=%d: %d<%d" % (l, e, c)
                          for (l, e, c, _, _) in table))
    e4 = dict((l, e) for (l, e, _, _, _) in table)[4]
    rep.check(e4 == 15 and xhl_e.get((4, 6, 6)) == 17 and e4 < 17,
              "[head-to-head] at the single instance K_{6,6}/C_8 this program "
              "certifies BOTH graphs as C_8-saturated: Xu-He-Lu's own with %s "
              "edges and the paper's with %d, so sat(K_{6,6},C_8) <= %d < 17"
              % (xhl_e.get((4, 6, 6)), e4, e4))

    print("")
    print("SCOPE NOTE -- exactly what is and is not machine-decided here.")
    print("  ESTABLISHED IN FULL: Conjecture 6.2 asserts an EQUALITY, so a single")
    print("  saturated subgraph with fewer edges refutes it.  The l = 4 instance")
    print("  alone does that, and it is decided here exhaustively: 15 < 17 with")
    print("  C_8-freeness and all 21 non-edges witnessed twice over.  l = 5,6,7")
    print("  are three further independent refutations, not supporting evidence.")
    print("  NOT DECIDED HERE (1): the 'for every l >= 4' half of Theorem 1, i.e.")
    print("  the cases l >= 8.  What is checked at each tested l is every general")
    print("  ingredient of the proof separately -- Lemma 2 and the tightness of")
    print("  its two index ranges, the sigma symmetry, the exhaustiveness of the")
    print("  four-case partition of the non-edges, the summand-by-summand edge")
    print("  count, and the structural cycle prediction (R-cycles cap at 2r,")
    print("  Q-cycles start at 2r+4) -- so the general argument is corroborated")
    print("  at each instance, but induction on l is not a computation.")
    print("  NOT DECIDED HERE (2), and not computable: whether the printed")
    print("  Discrete Math. 349 (2026) 114802 text agrees with arXiv:2410.11194v2")
    print("  in the quoted wording and in the labels 'Theorem 1.1'/'Conjecture")
    print("  6.2'.  Control (d) transcribes Definition 2.2 from the v2 source, so")
    print("  it inherits the same v2 scope the paper itself discloses.  What (d)")
    print("  DOES corroborate is that the formula n1+n2+l^2-3l+1 attributed to")
    print("  them is the edge count of their own construction, at three instances.")

    return rep.verdict()


if __name__ == "__main__":
    sys.exit(main())

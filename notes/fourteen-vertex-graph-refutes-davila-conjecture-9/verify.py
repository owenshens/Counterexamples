#!/usr/bin/env python3
"""Independent verification of the 14-vertex refutation of TxGraffiti Conjecture 9
(Z(G) <= gamma(G) + 2 for connected cubic diamond-free graphs) and of the
chain-family excess claims that accompany it.

Standard library only, exact integer arithmetic, no input files, no network.

------------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs, transcribed and then tested)
------------------------------------------------------------------------------
 P1  The graph G on A u B u C u D u {x,y}, A={a1,a2,a3}, ..., D={d1,d2,d3}, with
     E(G) = {a_i b_j : (i,j) != (3,1)} u {c_i d_j : (i,j) != (3,1)}
            u {a3 x, x b1, c3 y, y d1, x y}.
 P2  The statement acted on: Z(G) <= gamma(G) + 2 for connected cubic
     diamond-free graphs; "diamond" = K4 minus an edge, and the paper notes the
     subgraph/induced-subgraph ambiguity in the hypothesis.
 P3  The asserted values gamma(G) = 4 and Z(G) = 7.
 P4  The dominating witness {a3, b1, c3, d1}.
 P5  The zero forcing set S0 = {y, a1, a2, b2, c1, c2, d2} and the 7-step
     forcing sequence b2->a3, d2->c3, c3->d3, c1->d1, y->x, a3->b3, a1->b1.
 P6  The fort families used for Z >= 7: {a1,a2}, {b2,b3}, {a_i,a3,b1,b_j} for
     i in {1,2}, j in {2,3}; and {x,y} u p u q for p in P, q in Q with
     P = {{a1,a3},{a2,a3},{b1,b2},{b1,b3}}, Q the mirror in C u D.
 P7  The lower-bound counting argument gamma >= ceil(14/4) = 4.
 P8  The chain family: k copies of K_{3,3}, one edge of each end block and two
     edges of each middle block subdivided, consecutive subdivision vertices
     joined by bridges; connected, cubic, triangle-free, on 8k-2 vertices.
 P9  The excess Z - gamma equals 4 at n = 30 and 5 at n = 38, so no bound
     Z <= gamma + c with c <= 4 holds in the class.
 P10 The remarks that xy is a bridge and that x, y are centres of induced claws.

Reading fixed here: the paper does not say WHICH two edges of a middle block are
subdivided.  Up to the automorphism group of K_{3,3} there are two readings, an
independent pair and a pair sharing an endpoint.  This program uses the
independent pair; it also computes the other reading and prints it, because the
two readings give different domination numbers, and only the independent one
reproduces P9.

------------------------------------------------------------------------------
DERIVED HERE (nothing below is copied from the paper: every number is recomputed
from the edge rule P1 / construction P8 and then compared with the transcript)
------------------------------------------------------------------------------
 D1  G decoded: order, size, simplicity, degree sequence, connectivity.
 D2  Triangle-freeness and diamond-freeness under BOTH readings of F-freeness.
 D3  gamma(G) by exhaustive minimisation over all vertex subsets, and the
     counting bound of P7 recomputed from the actual closed neighbourhoods.
 D4  Z(G) by exhaustive minimisation over all vertex subsets (this is the
     load-bearing computation: no set of 6 vertices forces).
 D5  The refutation itself: Z(G) - (gamma(G) + 2) computed from D3 and D4.
 D6  Validity of the exhibited witnesses P4, P5 step by step.
 D7  The complete census of minimal forts of G (all 2^14 subsets tested), the
     exact minimum fort cover, and the two-vertex claim of P6.
 D8  The chain family for 2 <= k <= 8: order, cubicity, connectivity,
     triangle-freeness, and isomorphism of the k=2 member with G.
 D9  gamma exactly for 2 <= k <= 8 (branch and bound, cross-checked against
     exhaustive search where exhaustive search is affordable).
 D10 Z exactly for k = 2, 3, 4, 5 (orders 14, 22, 30, 38) by an exhibited
     forcing set together with a fort-cover lower bound: every zero forcing set
     meets every fort, so if no vertex set of size m that meets all enumerated
     forts forces, then Z > m; supersets of forcing sets force, so ruling out
     size exactly m rules out every size <= m at once.  The fort enumeration
     is complete for forts of size <= 10; at n = 14, where a census over all
     2^n subsets is affordable, that enumeration is required to agree exactly
     with the census, and the forcing witnesses are re-tested by a second,
     independently written propagator (round-parallel rather than one force at
     a time), the dominating witnesses by a second domination test.
 D11 The excess Z - gamma at n = 30 and n = 38, compared with P9.
 D12 Internal agreement: the domination numbers of the k = 2..8 table and the
     ones used for the excesses must come out equal, so the program cannot
     report two different chain families in a single run.
"""

import sys
import time
from itertools import combinations

sys.setrecursionlimit(100000)

CHECKS = []
T0 = time.time()


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    sys.stdout.flush()
    return bool(ok)


def pc(m):
    return bin(m).count("1")


def bitlist(m):
    out = []
    while m:
        low = m & -m
        out.append(low.bit_length() - 1)
        m ^= low
    return out


def mask(vs):
    m = 0
    for v in vs:
        m |= 1 << v
    return m


NAMES14 = ["a1", "a2", "a3", "b1", "b2", "b3",
           "c1", "c2", "c3", "d1", "d2", "d3", "x", "y"]
IDX14 = dict((s, i) for i, s in enumerate(NAMES14))


def build_g14():
    """Decode P1 into an edge list (nothing about the result is assumed)."""
    edges = []
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            if (i, j) != (3, 1):
                edges.append(("a%d" % i, "b%d" % j))
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            if (i, j) != (3, 1):
                edges.append(("c%d" % i, "d%d" % j))
    edges += [("a3", "x"), ("x", "b1"), ("c3", "y"), ("y", "d1"), ("x", "y")]
    pairs = [(IDX14[u], IDX14[v]) for u, v in edges]
    return len(NAMES14), pairs


def adjacency(n, pairs):
    adj = [0] * n
    for u, v in pairs:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def is_connected(adj, n):
    seen = 1
    stack = [0]
    while stack:
        v = stack.pop()
        nb = adj[v] & ~seen
        for w in bitlist(nb):
            seen |= 1 << w
            stack.append(w)
    return seen == (1 << n) - 1


def triangles(adj, n):
    """number of triangles, computed directly"""
    t = 0
    for u in range(n):
        for v in bitlist(adj[u]):
            if v > u:
                t += pc(adj[u] & adj[v] & ~((1 << (v + 1)) - 1))
    return t


def max_induced_edges_on_4_sets(adj, n):
    """max number of edges induced by a 4-subset; a diamond needs 5, K4 needs 6"""
    best = 0
    for q in combinations(range(n), 4):
        e = 0
        for a, b in combinations(q, 2):
            if adj[a] >> b & 1:
                e += 1
        if e > best:
            best = e
    return best


def closure(S, adj, full):
    """zero forcing closure: a blue vertex with exactly one white neighbour
       forces it blue, repeated to a fixed point"""
    blue = S
    while blue != full:
        moved = False
        for v in bitlist(blue):
            w = adj[v] & ~blue
            if w and (w & (w - 1)) == 0:
                blue |= w
                moved = True
                break
        if not moved:
            break
    return blue


def is_zfs(S, adj, full):
    return closure(S, adj, full) == full


def dominates(S, adj, full):
    cov = 0
    for v in bitlist(S):
        cov |= adj[v] | (1 << v)
    return cov == full


def closure_rounds(S, adj, full):
    """A second, independently written propagator: collect every legal force of
       a round and apply them together, repeating to a fixed point.  The zero
       forcing closure does not depend on the order of the forces, so this must
       agree with closure() on every input; a disagreement is a defect in one of
       the two, which is what the chain witness check tests for."""
    blue = S
    while True:
        add = 0
        for v in bitlist(blue):
            w = adj[v] & ~blue
            if w and (w & (w - 1)) == 0:
                add |= w
        if add & ~blue == 0:
            return blue
        blue |= add


def dominates_setwise(S, adj, n):
    """A second, independently written domination test: scan the vertices and
       demand that each is in S or has a neighbour in S."""
    for v in range(n):
        if not (S >> v & 1) and not (adj[v] & S):
            return False
    return True


def min_dominating_exhaustive(adj, n):
    full = (1 << n) - 1
    for k in range(0, n + 1):
        for c in combinations(range(n), k):
            if dominates(mask(c), adj, full):
                return k, mask(c)
    return None, 0


def min_zfs_exhaustive(adj, n):
    full = (1 << n) - 1
    for k in range(0, n + 1):
        for c in combinations(range(n), k):
            if is_zfs(mask(c), adj, full):
                return k, mask(c)
    return None, 0


def all_zfs_of_size(adj, n, k):
    full = (1 << n) - 1
    out = []
    for c in combinations(range(n), k):
        m = mask(c)
        if is_zfs(m, adj, full):
            out.append(m)
    return out


def is_fort(F, adj, full):
    """F nonempty and every vertex outside F has 0 or >= 2 neighbours in F"""
    if F == 0:
        return False
    for v in bitlist(full & ~F):
        c = adj[v] & F
        if c and (c & (c - 1)) == 0:
            return False
    return True


def minimal_forts_exhaustive(adj, n):
    """complete census of minimal forts, by testing every subset (n small)"""
    full = (1 << n) - 1
    forts = [F for F in range(1, 1 << n) if is_fort(F, adj, full)]
    forts.sort(key=pc)
    minimal = []
    for F in forts:
        if not any((g & F) == g for g in minimal):
            minimal.append(F)
    return forts, minimal


def minimize_fort(F, adj, full):
    changed = True
    while changed:
        changed = False
        for v in bitlist(F):
            g = F & ~(1 << v)
            if g and is_fort(g, adj, full):
                F = g
                changed = True
                break
    return F


def min_fort_cover(fam, n, ub):
    """exact minimum set of vertices meeting every member of fam"""
    fam = sorted(fam, key=pc)
    best = [ub, None]

    def bound(chosen):
        used = 0
        c = 0
        for f in fam:
            if (f & chosen) == 0 and (f & used) == 0:
                used |= f
                c += 1
        return c

    def rec(chosen, forb, cnt):
        un = None
        for f in fam:
            if not (f & chosen):
                un = f
                break
        if un is None:
            if cnt < best[0]:
                best[0] = cnt
                best[1] = chosen
            return
        if cnt + bound(chosen) >= best[0]:
            return
        av = un & ~forb
        seen = 0
        while av:
            low = av & -av
            av ^= low
            rec(chosen | low, forb | seen, cnt + 1)
            seen |= low

    rec(0, 0, 0)
    return best[0], best[1]


def chain(k, reading="independent"):
    """P8: k copies of K_{3,3}; one edge of each end block and two edges of each
       middle block subdivided; consecutive subdivision vertices bridged.
       reading fixes which two edges of a middle block are subdivided:
         'independent' -> a3b1 and a1b2 (no common endpoint)
         'sharing'     -> a3b1 and a3b2 (common endpoint a3)"""
    names = []
    edges = []

    def new(nm):
        names.append(nm)
        return nm

    blocks = []
    for t in range(k):
        aa = [new("a%d_%d" % (t, i)) for i in (1, 2, 3)]
        bb = [new("b%d_%d" % (t, i)) for i in (1, 2, 3)]
        blocks.append((aa, bb))
    subs = []
    for t in range(k):
        aa, bb = blocks[t]
        gone = set()
        left = right = None
        if t == 0:
            right = new("s%d_R" % t)
            edges += [(aa[2], right), (right, bb[0])]
            gone.add((2, 0))
        elif t == k - 1:
            left = new("s%d_L" % t)
            edges += [(aa[2], left), (left, bb[0])]
            gone.add((2, 0))
        else:
            e1 = (2, 0)
            e2 = (0, 1) if reading == "independent" else (2, 1)
            left = new("s%d_L" % t)
            right = new("s%d_R" % t)
            edges += [(aa[e1[0]], left), (left, bb[e1[1]])]
            edges += [(aa[e2[0]], right), (right, bb[e2[1]])]
            gone.add(e1)
            gone.add(e2)
        for i in range(3):
            for j in range(3):
                if (i, j) not in gone:
                    edges.append((aa[i], bb[j]))
        subs.append((left, right))
    for t in range(k - 1):
        edges.append((subs[t][1], subs[t + 1][0]))
    idx = dict((s, i) for i, s in enumerate(names))
    n = len(names)
    return n, [(idx[u], idx[v]) for u, v in edges], names


def isomorphic(a1, a2, n):
    """backtracking isomorphism test (n is small and both graphs are cubic)"""
    if sorted(pc(m) for m in a1) != sorted(pc(m) for m in a2):
        return False
    order = []
    seen = 0
    stack = [0]
    while stack:
        v = stack.pop()
        if seen >> v & 1:
            continue
        seen |= 1 << v
        order.append(v)
        for w in bitlist(a1[v] & ~seen):
            stack.append(w)
    for v in range(n):
        if not (seen >> v & 1):
            order.append(v)
    phi = [-1] * n

    def rec(i, used):
        if i == n:
            return True
        v = order[i]
        for w in range(n):
            if used >> w & 1 or pc(a1[v]) != pc(a2[w]):
                continue
            good = True
            for u in order[:i]:
                if bool(a1[v] >> u & 1) != bool(a2[w] >> phi[u] & 1):
                    good = False
                    break
            if good:
                phi[v] = w
                if rec(i + 1, used | (1 << w)):
                    return True
                phi[v] = -1
        return False

    return rec(0, 0)


def gamma_bb(adj, n):
    """exact domination number: branch on the closed neighbourhood of an
       uncovered vertex, with the ceil(uncovered / (max closed nbhd)) bound"""
    full = (1 << n) - 1
    cl = [adj[v] | (1 << v) for v in range(n)]
    wid = max(pc(c) for c in cl)
    best = [n + 1, None]

    def rec(cov, chosen, cnt):
        if cov == full:
            if cnt < best[0]:
                best[0] = cnt
                best[1] = chosen
            return
        rem = pc(full & ~cov)
        if cnt + (rem + wid - 1) // wid >= best[0]:
            return
        u = full & ~cov
        v = (u & -u).bit_length() - 1
        for w in bitlist(cl[v]):
            rec(cov | cl[w], chosen | (1 << w), cnt + 1)

    rec(0, 0, 0)
    return best[0], best[1]


def forts_upto(adj, n, maxsize, percap=400000):
    """every fort of size <= maxsize, by repairing violations: if a vertex u
       outside F has exactly one neighbour in F then any fort containing F must
       contain u or a second neighbour of u, so branch on those choices.
       Complete for the stated size bound unless the node cap trips."""
    full = (1 << n) - 1
    out = set()
    truncated = False
    for s in range(n):
        seen = set()
        stack = [1 << s]
        nodes = 0
        while stack:
            F = stack.pop()
            if F in seen:
                continue
            seen.add(F)
            nodes += 1
            if nodes > percap:
                truncated = True
                break
            if pc(F) > maxsize:
                continue
            bad = -1
            for v in bitlist(full & ~F):
                c = adj[v] & F
                if c and (c & (c - 1)) == 0:
                    bad = v
                    break
            if bad < 0:
                out.add(F)
                continue
            stack.append(F | (1 << bad))
            for w in bitlist(adj[bad] & ~F):
                stack.append(F | (1 << w))
    lst = sorted(out, key=pc)
    keep = []
    for f in lst:
        if not any((g & f) == g for g in keep):
            keep.append(f)
    return keep, truncated


def greedy_zfs(adj, n, start=0):
    """deterministic greedy upper bound on Z: repeatedly add the vertex whose
       addition enlarges the closure most"""
    full = (1 << n) - 1
    S = start
    while True:
        cl = closure(S, adj, full)
        if cl == full:
            return S
        bv, bs = -1, -1
        for v in bitlist(full & ~cl):
            sz = pc(closure(S | (1 << v), adj, full))
            if sz > bs:
                bs, bv = sz, v
        S |= 1 << bv


def best_zfs(adj, n):
    """smallest forcing set found by greedy from every single-vertex seed"""
    best = greedy_zfs(adj, n, 0)
    for v in range(n):
        S = greedy_zfs(adj, n, 1 << v)
        if pc(S) < pc(best):
            best = S
    return best


def forcing_cover_of_size(fam, n, m, adj, full, smallcap=80):
    """Search every vertex set of size exactly m that meets all forts in fam
       (each generated once, by forbidding earlier branch vertices) and return
       one that is zero forcing, else None.  Every zero forcing set meets every
       fort, and every superset of a forcing set forces, so None means
       Z > m."""
    fam = sorted(fam, key=pc)
    small = fam[:smallcap]
    leaves = [0]

    def bound(chosen):
        used = 0
        c = 0
        for f in small:
            if (f & chosen) == 0 and (f & used) == 0:
                used |= f
                c += 1
        return c

    def rec(chosen, forb, cnt):
        un = None
        for f in fam:
            if not (f & chosen):
                un = f
                break
        if un is None:
            if cnt == m:
                leaves[0] += 1
                return chosen if closure(chosen, adj, full) == full else None
            rest = [v for v in range(n) if not (chosen >> v & 1)]
            for extra in combinations(rest, m - cnt):
                S = chosen | mask(extra)
                leaves[0] += 1
                if closure(S, adj, full) == full:
                    return S
            return None
        if cnt >= m:
            return None
        av = un & ~forb
        if av == 0 or cnt + bound(chosen) > m:
            return None
        seen = 0
        while av:
            low = av & -av
            av ^= low
            r = rec(chosen | low, forb | seen, cnt + 1)
            if r is not None:
                return r
            seen |= low
        return None

    return rec(0, 0, 0), leaves[0]


PAPER_GAMMA = 4          # P3
PAPER_Z = 7              # P3
PAPER_DOM = ["a3", "b1", "c3", "d1"]                            # P4
PAPER_S0 = ["y", "a1", "a2", "b2", "c1", "c2", "d2"]            # P5
PAPER_SEQ = [("b2", "a3"), ("d2", "c3"), ("c3", "d3"),
             ("c1", "d1"), ("y", "x"), ("a3", "b3"),
             ("a1", "b1")]                                       # P5
PAPER_EXCESS = {30: 4, 38: 5}                                    # P9


def part_a(adj, n, pairs, full):
    """the exhibited object, its hypotheses, and the two invariants"""
    norm = set()
    loops = 0
    for u, v in pairs:
        if u == v:
            loops += 1
        norm.add((min(u, v), max(u, v)))
    print("decoded graph: order %d, size %d" % (n, len(norm)))
    for v in range(n):
        print("  %-3s : %s" % (NAMES14[v],
                               " ".join(NAMES14[w] for w in bitlist(adj[v]))))
    ck("g14_decodes_to_14_vertices_21_simple_edges",
       n == 14 and len(norm) == 21 and len(pairs) == 21 and loops == 0,
       "order=%d size=%d loops=%d" % (n, len(norm), loops))
    degs = [pc(adj[v]) for v in range(n)]
    ck("g14_is_cubic", set(degs) == set([3]),
       "degree multiset %s" % sorted(set(degs)))
    ck("g14_is_connected", is_connected(adj, n))
    tri = triangles(adj, n)
    ck("g14_is_triangle_free", tri == 0, "triangles=%d" % tri)
    me = max_induced_edges_on_4_sets(adj, n)
    ck("g14_is_diamond_free_subgraph_and_induced", me <= 4,
       "max edges on a 4-set = %d, a diamond needs 5" % me)

    domw = mask(IDX14[s] for s in PAPER_DOM)
    ck("g14_paper_dominating_witness_dominates",
       dominates(domw, adj, full) and pc(domw) == 4,
       "|D|=%d" % pc(domw))
    gmin, gw = min_dominating_exhaustive(adj, n)
    ck("g14_gamma_by_exhaustive_search_is_4", gmin == PAPER_GAMMA,
       "gamma=%d witness %s" % (gmin, [NAMES14[i] for i in bitlist(gw)]))
    wid = max(pc(adj[v] | (1 << v)) for v in range(n))
    need = -(-n // wid)
    ck("g14_counting_lower_bound_ceil_n_over_4_is_4",
       wid == 4 and need == 4 and gmin >= need,
       "max closed neighbourhood %d, ceil(%d/%d)=%d" % (wid, n, wid, need))

    s0 = mask(IDX14[s] for s in PAPER_S0)
    ck("g14_paper_forcing_set_S0_forces_all_of_G",
       pc(s0) == 7 and is_zfs(s0, adj, full),
       "|S0|=%d, closure size %d" % (pc(s0), pc(closure(s0, adj, full))))
    blue = s0
    seq_ok = True
    detail = ""
    for u, w in PAPER_SEQ:
        iu, iw = IDX14[u], IDX14[w]
        white = adj[iu] & ~blue
        if not (blue >> iu & 1) or (blue >> iw & 1) or white != (1 << iw):
            seq_ok = False
            detail = "step %s->%s invalid" % (u, w)
            break
        blue |= 1 << iw
    if seq_ok and blue != full:
        seq_ok = False
        detail = "sequence ends with %d blue" % pc(blue)
    ck("g14_paper_forcing_sequence_is_step_by_step_valid", seq_ok,
       detail or "7 forces, all 14 vertices blue")

    six = all_zfs_of_size(adj, n, 6)
    ck("g14_no_six_vertex_set_is_zero_forcing", len(six) == 0,
       "forcing 6-sets found: %d of %d" % (len(six), 3003))
    zmin, zw = min_zfs_exhaustive(adj, n)
    ck("g14_Z_by_exhaustive_search_is_7", zmin == PAPER_Z,
       "Z=%d witness %s" % (zmin, [NAMES14[i] for i in bitlist(zw)]))
    exc = zmin - (gmin + 2)
    ck("refutes_conjecture_9_Z_exceeds_gamma_plus_2", exc > 0,
       "Z=%d, gamma+2=%d, violation by %d" % (zmin, gmin + 2, exc))
    return gmin, zmin


def part_b(adj, n, full, zmin):
    """the fort route of the paper, and its two structural remarks"""
    u_forts = [["a1", "a2"], ["b2", "b3"]]
    for i in ("a1", "a2"):
        for j in ("b2", "b3"):
            u_forts.append([i, "a3", "b1", j])
    ok = all(is_fort(mask(IDX14[s] for s in F), adj, full) for F in u_forts)
    ck("g14_the_six_forts_inside_A_union_B_are_forts", ok and len(u_forts) == 6,
       "%d sets tested against the definition" % len(u_forts))

    P = [["a1", "a3"], ["a2", "a3"], ["b1", "b2"], ["b1", "b3"]]
    Q = [["c1", "c3"], ["c2", "c3"], ["d1", "d2"], ["d1", "d3"]]
    big = []
    for p in P:
        for q in Q:
            big.append(["x", "y"] + p + q)
    okb = all(is_fort(mask(IDX14[s] for s in F), adj, full) for F in big)
    ck("g14_the_sixteen_forts_xy_union_p_union_q_are_forts",
       okb and len(big) == 16, "%d sets tested" % len(big))

    hits = 0
    tgt = [mask(IDX14[s] for s in F) for F in u_forts]
    for pair in combinations(range(n), 2):
        m = mask(pair)
        if all(m & t for t in tgt):
            hits += 1
    cov3, _ = min_fort_cover(tgt, n, n + 1)
    ck("g14_no_two_vertices_meet_all_six_A_B_forts",
       hits == 0 and cov3 == 3,
       "2-sets meeting all six: %d; exact minimum cover of the six: %d"
       % (hits, cov3))

    allf, minf = minimal_forts_exhaustive(adj, n)
    wmirror = []
    for F in u_forts:
        wmirror.append([{"a": "c", "b": "d"}[s[0]] + s[1] for s in F])
    predicted = set()
    for F in u_forts + wmirror + big:
        predicted.add(mask(IDX14[s] for s in F))
    ck("g14_minimal_fort_census_is_exactly_the_paper_family",
       set(minf) == predicted and len(minf) == 28,
       "forts %d, minimal forts %d, predicted %d, sizes %s"
       % (len(allf), len(minf), len(predicted),
          sorted(set(pc(f) for f in minf))))

    cov, cw = min_fort_cover(minf, n, n + 1)
    ck("g14_minimum_fort_cover_equals_Z", cov == zmin == 7,
       "minimum cover %d, Z %d, cover %s"
       % (cov, zmin, [NAMES14[i] for i in bitlist(cw)]))

    mins = all_zfs_of_size(adj, n, zmin)
    lemma = all(all(S & F for S in mins) for F in minf)
    uhalf = mask(range(0, 6))
    whalf = mask(range(6, 12))
    halves = all(pc(S & uhalf) >= 3 and pc(S & whalf) >= 3 for S in mins)
    ck("g14_every_minimum_forcing_set_meets_every_fort_and_each_half_thrice",
       lemma and halves and len(mins) > 0,
       "%d forcing %d-sets, each against %d minimal forts; min |S n (A u B)| "
       "= %d, min |S n (C u D)| = %d"
       % (len(mins), zmin, len(minf),
          min([pc(S & uhalf) for S in mins] or [0]),
          min([pc(S & whalf) for S in mins] or [0])))

    ix, iy = IDX14["x"], IDX14["y"]
    cut = [a for a in adj]
    cut[ix] &= ~(1 << iy)
    cut[iy] &= ~(1 << ix)
    ck("g14_edge_xy_is_a_bridge",
       is_connected(adj, n) and not is_connected(cut, n))
    claw = True
    for v in (ix, iy):
        nb = bitlist(adj[v])
        if len(nb) != 3 or any(adj[a] >> b & 1 for a, b in combinations(nb, 2)):
            claw = False
    ck("g14_x_and_y_are_centres_of_induced_claws", claw,
       "N(x)=%s N(y)=%s, both independent"
       % ([NAMES14[w] for w in bitlist(adj[ix])],
          [NAMES14[w] for w in bitlist(adj[iy])]))


def chain_gamma_lower_by_exhaustion(adj, n, k):
    """no vertex set of size k-1 dominates (used to cross-check the branch and
       bound where the binomial count is affordable)"""
    full = (1 << n) - 1
    for c in combinations(range(n), k - 1):
        if dominates(mask(c), adj, full):
            return False
    return True


def chain_structure(kmax=8):
    orders = []
    struct = []
    for k in range(2, kmax + 1):
        nk, pk, _ = chain(k)
        ak = adjacency(nk, pk)
        orders.append((k, nk, 8 * k - 2))
        struct.append((k, set(pc(m) for m in ak) == set([3]),
                       is_connected(ak, nk), triangles(ak, nk) == 0,
                       len(set((min(u, v), max(u, v)) for u, v in pk))))
    ck("chain_family_has_order_8k_minus_2_for_k_2_to_8",
       all(nk == exp for _, nk, exp in orders),
       " ".join("k=%d:n=%d" % (k, nk) for k, nk, _ in orders))
    ck("chain_family_is_cubic_connected_triangle_free_for_k_2_to_8",
       all(c and cn and tf for _, c, cn, tf, _ in struct),
       " ".join("k=%d:m=%d" % (k, m) for k, _, _, _, m in struct))
    return struct


def chain_invariants(k, reading="independent", exact_z=True):
    nk, pk, nmk = chain(k, reading)
    ak = adjacency(nk, pk)
    full = (1 << nk) - 1
    g, gw = gamma_bb(ak, nk)
    gw_ok = (dominates(gw, ak, full) and pc(gw) == g
             and dominates_setwise(gw, ak, nk))
    S = best_zfs(ak, nk)
    ub = pc(S)
    forcing_ok = is_zfs(S, ak, full) and closure_rounds(S, ak, full) == full
    lb = None
    leaves = 0
    fam = []
    fam_ok = True
    fam_complete = None
    if exact_z:
        fam, trunc = forts_upto(ak, nk, 10)
        fam_ok = all(is_fort(f, ak, full) for f in fam) and not trunc
        if nk <= 16:
            # the lower bound engine is only as good as this enumeration, so at
            # the one order where a complete census over all 2^n subsets is
            # affordable, demand that the two agree exactly
            _, minf = minimal_forts_exhaustive(ak, nk)
            fam_complete = set(fam) == set(f for f in minf if pc(f) <= 10)
        r, leaves = forcing_cover_of_size(fam, nk, ub - 1, ak, full)
        lb = ub if r is None else pc(r)
    return dict(k=k, n=nk, adj=ak, names=nmk, gamma=g, gw=gw, gw_ok=gw_ok,
                S=S, ub=ub, lb=lb, forcing_ok=forcing_ok, fam=fam,
                fam_ok=fam_ok, fam_complete=fam_complete, leaves=leaves)


def part_c():
    chain_structure(8)
    n14, p14 = build_g14()
    a14 = adjacency(n14, p14)
    n2, p2, _ = chain(2)
    a2 = adjacency(n2, p2)
    ck("chain_member_k2_is_isomorphic_to_the_exhibited_graph",
       n2 == n14 and isomorphic(a2, a14, n14),
       "both cubic on %d vertices" % n14)

    gam = {}
    for k in range(2, 9):
        nk, pk, _ = chain(k)
        ak = adjacency(nk, pk)
        g, gw = gamma_bb(ak, nk)
        gam[k] = (g, dominates(gw, ak, (1 << nk) - 1))
    ck("chain_domination_numbers_computed_with_verified_witnesses_k_2_to_8",
       all(okw for _, okw in gam.values()) and gam[2][0] == 4,
       " ".join("k=%d:gamma=%d" % (k, gam[k][0]) for k in sorted(gam)))
    n3, p3, _ = chain(3)
    a3 = adjacency(n3, p3)
    ck("chain_gamma_lower_bounds_cross_checked_by_exhaustive_search",
       chain_gamma_lower_by_exhaustion(a3, n3, gam[3][0])
       and chain_gamma_lower_by_exhaustion(a2, n2, gam[2][0]),
       "no %d-set dominates at n=22 and no %d-set at n=14"
       % (gam[3][0] - 1, gam[2][0] - 1))

    res = {}
    for k in (2, 3, 4, 5):
        r = chain_invariants(k)
        res[k] = r
        print("  chain k=%d n=%d: gamma=%d  Z in [%s,%d]  forts=%d leaves=%d"
              % (k, r["n"], r["gamma"], r["lb"], r["ub"], len(r["fam"]),
                 r["leaves"]))
    ck("chain_enumerated_sets_are_genuine_forts_and_census_not_truncated",
       all(r["fam_ok"] for r in res.values())
       and res[2]["fam_complete"] is True,
       "families of sizes %s; at n=14 the size <= 10 enumeration equals the "
       "exhaustive minimal-fort census over all 2^14 subsets"
       % [len(res[k]["fam"]) for k in sorted(res)])
    ck("chain_exhibited_dominating_and_forcing_witnesses_are_valid",
       all(r["forcing_ok"] and r["gw_ok"] for r in res.values())
       and all(res[k]["gamma"] == gam[k][0] for k in res),
       "%s ; two independent propagators agree and the k=2..8 domination "
       "table is the same graph family as the one used for the excesses"
       % " ".join("k=%d:|D|=%d,|S|=%d" % (k, res[k]["gamma"], res[k]["ub"])
                  for k in sorted(res)))
    ck("chain_Z_determined_exactly_for_k_2_to_5",
       all(res[k]["lb"] == res[k]["ub"] for k in res),
       " ".join("k=%d:Z=%d" % (k, res[k]["ub"]) for k in sorted(res)))
    ck("chain_k2_member_reproduces_gamma_4_Z_7",
       res[2]["gamma"] == 4 and res[2]["ub"] == 7)
    for k, nn in ((4, 30), (5, 38)):
        e = res[k]["ub"] - res[k]["gamma"]
        ck("chain_excess_Z_minus_gamma_at_n_%d_is_%d" % (nn, PAPER_EXCESS[nn]),
           res[k]["n"] == nn and e == PAPER_EXCESS[nn],
           "n=%d gamma=%d Z=%d excess=%d, paper says %d"
           % (res[k]["n"], res[k]["gamma"], res[k]["ub"], e,
              PAPER_EXCESS[nn]))
    beaten = []
    for c in range(0, 5):
        beaten.append(any(res[k]["ub"] > res[k]["gamma"] + c for k in res))
    ck("no_bound_Z_le_gamma_plus_c_with_c_at_most_4_holds_in_the_class",
       all(beaten),
       "c=0..4 each refuted inside the family; largest excess %d"
       % max(res[k]["ub"] - res[k]["gamma"] for k in res))

    alt = chain_invariants(4, "sharing", exact_z=False)
    ck("middle_block_reading_matters_sharing_pair_gives_excess_at_most_3",
       alt["n"] == 30 and alt["gamma"] == 10 and alt["ub"] - alt["gamma"] <= 3,
       "sharing reading at n=30: gamma=%d, Z<=%d, excess<=%d (independent "
       "reading is the one that gives %d)"
       % (alt["gamma"], alt["ub"], alt["ub"] - alt["gamma"],
          PAPER_EXCESS[30]))
    print("SCOPE DISCLOSURE: the chain family as described in the paper is "
          "under-determined.  \"two edges of each middle block\" does not say "
          "WHICH two, and up to Aut(K_{3,3}) there are two readings: an "
          "independent pair and a pair sharing an endpoint.  Both give "
          "connected cubic triangle-free graphs on 8k-2 vertices, and the two "
          "readings agree at k=2, which has no middle block, so the "
          "14-vertex graph and its gamma=%d, Z=%d are unaffected.  They can "
          "differ once a middle block exists (k >= 3), and at k=4 they do: "
          "the excesses %d at n=30 and %d at n=38 checked above hold for the "
          "independent-pair reading, while the sharing-pair reading gives "
          "gamma=%d and Z<=%d at n=30, i.e. excess at most %d.  A reader who "
          "reconstructs the family the other way will not reproduce the "
          "paper's %d at n=30."
          % (res[2]["gamma"], res[2]["ub"],
             res[4]["ub"] - res[4]["gamma"], res[5]["ub"] - res[5]["gamma"],
             alt["gamma"], alt["ub"], alt["ub"] - alt["gamma"],
             res[4]["ub"] - res[4]["gamma"]))
    sys.stdout.flush()
    return res, alt


def main():
    n, pairs = build_g14()
    adj = adjacency(n, pairs)
    full = (1 << n) - 1
    gmin, zmin = part_a(adj, n, pairs, full)
    part_b(adj, n, full, zmin)
    ck("g14_computed_values_match_the_paper_gamma_4_and_Z_7",
       gmin == PAPER_GAMMA and zmin == PAPER_Z,
       "computed gamma=%d Z=%d, paper gamma=%d Z=%d"
       % (gmin, zmin, PAPER_GAMMA, PAPER_Z))
    res, alt = part_c()
    print("NOT RE-RUN: (i) none of the attributions the paper's framing rests "
          "on can be checked by an offline program, and none of them is "
          "checked here -- that the bound Z <= gamma + 2 for connected cubic "
          "diamond-free graphs is Conjecture 9, under that number, of "
          "arXiv:2406.19231v2; that Zenodo DOI 10.5281/zenodo.21269439 "
          "resolves at version 3 of July 8, 2026; and that the sub-path cited "
          "in the bibliography exists inside that record: a reader must check "
          "all three at the sources, and this program's PASS lines say "
          "nothing about them; (ii) the two independent implementations, the "
          "mutation tests and the independent re-verification reported in the "
          "cited artifact are not re-executed -- every value above is instead "
          "recomputed here from the edge rule; (iii) the chain family is "
          "described in the paper only up to an ambiguity: it does not say "
          "which two edges of a middle block are subdivided, so the excesses "
          "%d at n=30 and %d at n=38 are verified here for the "
          "independent-pair reading alone, and under the other reading (a "
          "pair sharing an endpoint) this program computes gamma=%d and "
          "Z<=%d at n=30, an excess of at most %d, so the paper's reported %d "
          "at n=30 holds under one of the two readings and fails under the "
          "other; (iv) Z is determined here for 2 <= k <= 5 of the chain "
          "family and only gamma for 6 <= k <= 8, which is the range the "
          "paper itself reports."
          % (res[4]["ub"] - res[4]["gamma"], res[5]["ub"] - res[5]["gamma"],
             alt["gamma"], alt["ub"], alt["ub"] - alt["gamma"],
             res[4]["ub"] - res[4]["gamma"]))
    print("elapsed %.1f s" % (time.time() - T0))


if __name__ == "__main__":
    main()
    nbad = sum(1 for _, ok in CHECKS if not ok)
    n = len(CHECKS)
    if nbad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (nbad, n))
        raise SystemExit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    raise SystemExit(0)

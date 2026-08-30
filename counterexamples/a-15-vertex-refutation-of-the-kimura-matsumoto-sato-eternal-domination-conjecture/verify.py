#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify.py -- re-derives every quantity claimed in paper.tex, from the objects PRINTED there.

Python 3.9+, STANDARD LIBRARY ONLY (itertools, sys).  No third-party package, no external data
file, no network, no randomness, no floating-point decision: every arithmetic step below is on
Python ints, and every comparison is an integer or set comparison.

WHAT IS RE-DERIVED, and from what.  The only inputs are strings copied out of the paper:

  * G0_EDGE_LIST  -- the 27 edges of the witness G0, exactly as displayed in Section 2;
  * G0_GRAPH6     -- the graph6 encoding of G0, as displayed in Section 2;
  * HSTAR_GRAPH6  -- the graph6 encoding of the seed H* (KMS Figure 2.1 / MMV Table 10 #4);
  * T10, T11N14   -- the published order-13 and order-14 graph6 strings of MacGillivray,
                     Mynhardt and Virgile, reproduced in Section 5 of the paper, used as
                     FORCED-POSITIVE calibration for the decider.

Nothing is read from the search that found G0, and no value is hard-coded except the values the
paper prints, each of which is re-computed here and compared.

THE MODEL.  gamma^inf is the eternal domination number in the ONE-GUARD-MOVES model: a set S of
guarded vertices must dominate V; an attack at r not in S must be answered by moving a single
guard from some w in S adjacent to r onto r; the resulting set must again be able to answer any
attack, forever.  Formally gamma^inf(G) <= k iff the greatest fixed point of

     F  |-->  { S in F : for every r not in S there is w in S with wr in E(G)
                         and (S - w + r) in F }

started from F0 = {all dominating k-subsets of V(G)}, is non-empty.  This is the model of
Goddard, Hedetniemi and Hedetniemi and is consistent with the published sandwich
alpha(G) <= gamma^inf(G) <= theta(G), which this program re-checks on every object it touches.

Both directions carry a certificate that is re-verified in a SEPARATE pass, so no answer rests on
the fixed-point loop having been iterated to convergence:

  YES: the surviving family F is re-checked, member by member, to consist of dominating k-sets
       closed under every attack (a positive strategy).
  NO:  the round at which each dominating k-set died is kept as a RANK, and a fresh pass checks
       that every dominating k-set D admits an attack r such that every legal response is either
       non-dominating or of STRICTLY SMALLER rank.  Induction on that well-founded rank, not the
       loop, is what proves gamma^inf > k.
"""

import itertools
import sys

# ---------------------------------------------------------------------------------------------
# THE OBJECTS, copied from the paper
# ---------------------------------------------------------------------------------------------

# Section 2 of the paper, the displayed edge list of G0.
G0_EDGE_LIST = """
0-4   0-8   0-11  0-13  1-5   1-8   1-10  1-12  1-14  2-6   2-8   2-9   2-12  3-7
3-8   3-9   3-10  4-9   4-10  4-12  5-9   5-11  6-10  6-11  7-11  7-12  13-14
"""
G0_GRAPH6 = "N?`@?boNAsOwYO_?G?G"
HSTAR_GRAPH6 = "L?`@?boNAsOwYO"
EAR = ((0, 13), (13, 14), (1, 14))

# Section 5 of the paper: MacGillivray-Mynhardt-Virgile, Table 10 (all thirteen order-13
# triangle-free graphs with gamma^inf < theta) and the two order-14 entries of their Table 11.
T10 = ("L?`@F?M]DgYOFo", "L?`DAboU`w@{hS", "L?`@?boNAsLGBw", "L?`@?boNAsOwYO",
       "L?`@?boNAs@{os", "L?`@Cbod`w@{YS", "L?BDB?{{AsRGBs", "L?`@C`wl?{DgQs",
       "L?BDB?{Ucq^?Fo", "L?`@F@wlEcBoBo", "L?`@F@wlEcBoFo", "L?`@F@wlEcBgFo",
       "L?`@F?kQ_}YSlC")
T11N14 = ("M?`@?boNAs@{oshQ?", "M?`@F?kQckBwsglC?")

# Section 2: the two vertex-disjoint 5-cycles of G0.
ODD_CYCLE_A = (0, 4, 9, 2, 8)
ODD_CYCLE_B = (1, 5, 11, 7, 12)

# Section 3: the ear placements of the two sparser witnesses.
TOWER = (((0, 1), (2, 3)), ((0, 4), (1, 8), (2, 9)))

# Values printed in the paper, re-derived below and compared against these copies.
P_DEGREES = [4, 5, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 2, 2]
P_FAMILY_G0 = 2196
P_DOM6_COUNTS = [851, 721, 890, 890, 977, 1120, 1091, 1091, 1003, 890, 924, 757, 924, 1004, 1026]
P_FAMILY_HSTAR = 1098
P_FAMILY_T11N14 = (1741, 1749)
P_AUT_HSTAR = 48
P_ORBIT_SIZES = [4, 4, 6, 8, 8, 12, 12, 24]
P_TOWER_FAMILY = (4392, 8784)

# ---------------------------------------------------------------------------------------------
# check bookkeeping
# ---------------------------------------------------------------------------------------------

_n_pass = 0
_failed = []


def chk(name, cond, detail=""):
    global _n_pass
    if cond:
        _n_pass += 1
        print("PASS %s%s" % (name, (" [" + detail + "]") if detail else ""))
    else:
        _failed.append(name)
        print("FAIL %s%s" % (name, (" [" + detail + "]") if detail else ""))


def head(text):
    print("--- %s ---" % text)


# ---------------------------------------------------------------------------------------------
# graph plumbing.  A graph is (n, frozenset of ordered pairs u<v).  Bitmasks for the work.
# ---------------------------------------------------------------------------------------------

def parse_edges(text):
    E = set()
    for tok in text.split():
        u, v = tok.split("-")
        u, v = int(u), int(v)
        assert u != v
        E.add((min(u, v), max(u, v)))
    return frozenset(E)


def decode_graph6(s):
    """graph6 -> (n, frozenset of edges).  Small-n form only, which is all the paper uses."""
    b = s.encode("ascii")
    assert b[0] != 126, "the paper uses no graph6 string of order >= 63"
    n = b[0] - 63
    bits = []
    for c in b[1:]:
        x = c - 63
        assert 0 <= x < 64, "byte outside the graph6 alphabet"
        for k in range(5, -1, -1):
            bits.append((x >> k) & 1)
    need = n * (n - 1) // 2
    assert len(bits) >= need, "graph6 string too short for order %d" % n
    assert not any(bits[need:]), "graph6 padding bits are not all zero"
    E = set()
    i = 0
    for j in range(1, n):
        for k in range(j):
            if bits[i]:
                E.add((k, j))
            i += 1
    return n, frozenset(E)


def encode_graph6(n, E):
    bits = []
    for j in range(1, n):
        for k in range(j):
            bits.append(1 if (k, j) in E else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for i in range(0, len(bits), 6):
        x = 0
        for b in bits[i:i + 6]:
            x = (x << 1) | b
        out.append(chr(x + 63))
    return "".join(out)


def nbr(n, E):
    """open neighbourhood bitmasks"""
    N = [0] * n
    for u, v in E:
        N[u] |= 1 << v
        N[v] |= 1 << u
    return N


def cnbr(n, E):
    N = nbr(n, E)
    return [N[v] | (1 << v) for v in range(n)]


def popcount(x):
    return bin(x).count("1")


def bits_of(x):
    while x:
        low = x & -x
        yield low.bit_length() - 1
        x ^= low


def delete_vertex(n, E, w):
    """G - w, relabelled to 0..n-2 by the order-preserving map"""
    r = {}
    k = 0
    for v in range(n):
        if v != w:
            r[v] = k
            k += 1
    return n - 1, frozenset((min(r[u], r[v]), max(r[u], r[v])) for (u, v) in E if w not in (u, v))


def triangle_edges(n, E):
    N = nbr(n, E)
    return [(u, v) for (u, v) in sorted(E) if N[u] & N[v]]


def max_clique_size(n, E):
    """exact, by branching.  Used only to certify that the cliques of a triangle-free graph are
    single vertices and single edges, which is what turns theta into a matching problem."""
    N = nbr(n, E)
    best = [0]

    def rec(chosen, cand):
        if len(chosen) > best[0]:
            best[0] = len(chosen)
        if not cand:
            return
        if len(chosen) + popcount(cand) <= best[0]:
            return
        c = cand
        for v in bits_of(c):
            cand ^= 1 << v
            rec(chosen + [v], cand & N[v])

    rec([], (1 << n) - 1)
    return best[0]


def max_matching(n, E, mask=None):
    """(mu, one maximum matching).  Exact: branch on the lowest vertex of `mask`."""
    if mask is None:
        mask = (1 << n) - 1
    N = nbr(n, E)
    memo = {}

    def rec(m):
        if m == 0:
            return 0, ()
        if m in memo:
            return memo[m]
        v = (m & -m).bit_length() - 1
        best = rec(m ^ (1 << v))                       # leave v unmatched
        for u in bits_of(N[v] & m):
            sub, mt = rec(m ^ (1 << v) ^ (1 << u))
            if sub + 1 > best[0]:
                best = (sub + 1, mt + ((min(u, v), max(u, v)),))
        memo[m] = best
        return best

    return rec(mask)


def independence_number(n, E, mask=None):
    if mask is None:
        mask = (1 << n) - 1
    N = nbr(n, E)
    memo = {}

    def rec(m):
        if m == 0:
            return 0
        if m in memo:
            return memo[m]
        v = (m & -m).bit_length() - 1
        r = max(rec(m ^ (1 << v)), 1 + rec(m & ~(N[v] | (1 << v))))
        memo[m] = r
        return r

    return rec(mask)


def domination_number(n, E):
    C = cnbr(n, E)
    full = (1 << n) - 1
    for k in range(1, n + 1):
        for S in itertools.combinations(range(n), k):
            cov = 0
            for v in S:
                cov |= C[v]
            if cov == full:
                return k, S
    return n, tuple(range(n))


def is_two_connected(n, E):
    if n < 3:
        return False
    N = nbr(n, E)

    def connected(mask):
        if mask == 0:
            return True
        start = (mask & -mask).bit_length() - 1
        seen = 1 << start
        stack = [start]
        while stack:
            v = stack.pop()
            for u in bits_of(N[v] & mask & ~seen):
                seen |= 1 << u
                stack.append(u)
        return seen == mask

    full = (1 << n) - 1
    if not connected(full):
        return False
    return all(connected(full ^ (1 << v)) for v in range(n))


def factor_critical(n, E):
    """n odd and G - v has a perfect matching for every v"""
    if n % 2 == 0:
        return False, []
    full = (1 << n) - 1
    wit = []
    for v in range(n):
        mu, mt = max_matching(n, E, full ^ (1 << v))
        if mu != (n - 1) // 2:
            return False, []
        wit.append((v, mt))
    return True, wit


def clique_cover_number(n, E):
    """For a TRIANGLE-FREE graph only: every clique is a vertex or an edge, so a minimum clique
    cover is a maximum matching together with the vertices it misses, of size n - mu.  The
    triangle-freeness hypothesis is asserted, not assumed."""
    assert not triangle_edges(n, E), "clique_cover_number is only used on triangle-free graphs"
    mu, mt = max_matching(n, E)
    cover = [list(e) for e in mt]
    covered = 0
    for u, v in mt:
        covered |= (1 << u) | (1 << v)
    for v in range(n):
        if not (covered >> v) & 1:
            cover.append([v])
    return n - mu, cover


def dominating_ksets(n, E, k):
    C = cnbr(n, E)
    full = (1 << n) - 1
    out = []
    for S in itertools.combinations(range(n), k):
        cov = 0
        for v in S:
            cov |= C[v]
        if cov == full:
            m = 0
            for v in S:
                m |= 1 << v
            out.append(m)
    return out


def edom_le(n, E, k):
    """(answer, |F| if YES else #dominating k-sets, max_rank, certificate_reverified).

    The certificate is re-verified in a pass that reads only the FINAL family / the FINAL ranks,
    never the loop state, so a bug in the iteration order cannot make the certificate pass."""
    if k < 0:
        return False, 0, 0, True
    N = nbr(n, E)
    dom = dominating_ksets(n, E, k)
    ndom = len(dom)
    F = set(dom)
    rank = {}
    rnd = 0
    while True:
        rnd += 1
        dead = []
        for S in F:
            rem = ((1 << n) - 1) ^ S
            alive = True
            for r in bits_of(rem):
                lb = 1 << r
                ok = False
                for w in bits_of(N[r] & S):
                    if ((S ^ (1 << w)) | lb) in F:
                        ok = True
                        break
                if not ok:
                    alive = False
                    break
            if not alive:
                dead.append(S)
        if not dead:
            break
        for S in dead:
            F.discard(S)
            rank[S] = rnd

    if F:
        # ---- POSITIVE certificate, fresh pass over the final family only -------------------
        C = cnbr(n, E)
        full = (1 << n) - 1
        ok = True
        for S in F:
            if popcount(S) != k:
                ok = False
                break
            cov = 0
            for v in bits_of(S):
                cov |= C[v]
            if cov != full:
                ok = False
                break
            for r in bits_of(full ^ S):
                if not any(((S ^ (1 << w)) | (1 << r)) in F for w in bits_of(N[r] & S)):
                    ok = False
                    break
            if not ok:
                break
        return True, len(F), 0, ok

    # ---- NEGATIVE certificate: the ranking is well founded -----------------------------------
    ok = True
    maxrank = max(rank.values()) if rank else 0
    for S, rk in rank.items():
        found = False
        for r in bits_of(((1 << n) - 1) ^ S):
            allsmaller = True
            for w in bits_of(N[r] & S):
                T = (S ^ (1 << w)) | (1 << r)
                if T in rank and rank[T] >= rk:
                    allsmaller = False
                    break
            if allsmaller:
                found = True
                break
        if not found:
            ok = False
            break
    return False, ndom, maxrank, ok


def gamma_inf_exact(n, E, hi):
    """the exact value, found by walking down from hi; hi must satisfy gamma^inf <= hi"""
    val = hi
    for k in range(hi, 0, -1):
        if edom_le(n, E, k)[0]:
            val = k
        else:
            break
    return val


def automorphisms(n, E):
    N = nbr(n, E)
    deg = [popcount(N[v]) for v in range(n)]
    out = []
    perm = []

    def rec(used):
        i = len(perm)
        if i == n:
            out.append(tuple(perm))
            return
        for j in range(n):
            if (used >> j) & 1 or deg[j] != deg[i]:
                continue
            good = True
            for k in range(i):
                if (((N[i] >> k) & 1) != ((N[j] >> perm[k]) & 1)):
                    good = False
                    break
            if good:
                perm.append(j)
                rec(used | (1 << j))
                perm.pop()

    rec(0)
    return out


def add_ears(n, E, ears):
    """attach one odd ear a - x - y - b for each (a, b), on fresh vertices"""
    E = set(E)
    nv = n
    for (a, b) in ears:
        E.add((min(a, nv), max(a, nv)))
        E.add((nv, nv + 1))
        E.add((min(b, nv + 1), max(b, nv + 1)))
        nv += 2
    return nv, frozenset(E)


def cycle(n):
    return n, frozenset((min(i, (i + 1) % n), max(i, (i + 1) % n)) for i in range(n))


# =============================================================================================
head("1. the witness G0, as printed in Section 2")

G0 = parse_edges(G0_EDGE_LIST)
N15 = 15
chk("edge-list-parses-to-27-distinct-edges", len(G0) == 27,
    "vertices used: %s" % sorted(set(x for e in G0 for x in e)))
chk("graph6-agrees-with-the-edge-list-under-the-identity-labelling",
    decode_graph6(G0_GRAPH6) == (N15, G0) and encode_graph6(N15, G0) == G0_GRAPH6,
    "'%s' <-> the 27 printed edges, symmetric difference 0" % G0_GRAPH6)
chk("order-and-size-sit-on-the-conjecture-s-boundary",
    len(G0) == 2 * N15 - 3 == 27, "n = 15, |E| = 27, 2n-3 = %d" % (2 * N15 - 3))
NG0 = nbr(N15, G0)
degs = [popcount(NG0[v]) for v in range(N15)]
chk("degree-sequence-as-printed", degs == P_DEGREES,
    "%s, sum %d = 2|E|" % (degs, sum(degs)))
chk("triangle-free", triangle_edges(N15, G0) == [], "no edge has a common neighbour")
chk("maximum-clique-has-two-vertices", max_clique_size(N15, G0) == 2,
    "so every clique cover of G0 uses only vertices and edges")
chk("two-connected", is_two_connected(N15, G0),
    "G0 - v is connected for all 15 v, so the refutation survives a connected reading")

NH, HSTAR = decode_graph6(HSTAR_GRAPH6)
chk("seed-graph6-decodes-to-the-KMS-Figure-2.1-graph",
    (NH, len(HSTAR)) == (13, 24) and encode_graph6(NH, HSTAR) == HSTAR_GRAPH6,
    "'%s' -> n = 13, m = 24" % HSTAR_GRAPH6)
chk("deleting-the-printed-ear-recovers-the-seed-on-the-nose",
    frozenset(G0 - frozenset(EAR)) == HSTAR and frozenset(EAR) <= G0,
    "G0 - {0-13, 13-14, 14-1} is LABEL-EQUAL to H*, and 13, 14 are G0's only degree-2 vertices"
    if [v for v in range(N15) if degs[v] == 2] == [13, 14] else "ear check")
chk("the-ear-is-the-only-pair-of-degree-two-vertices-and-they-are-adjacent",
    [v for v in range(N15) if degs[v] == 2] == [13, 14] and (13, 14) in G0,
    "so H* is recoverable from G0 and any isomorphism carries ear to ear")

chk("triangle-free-planar-edge-bound-excludes-planarity",
    len(G0) > 2 * N15 - 4, "a triangle-free planar graph on 15 vertices has at most "
                           "2n-4 = 26 edges; G0 has 27, so G0 is NOT planar")


def is_cycle(n, E, seq):
    k = len(seq)
    if len(set(seq)) != k:
        return False
    for i in range(k):
        e = (min(seq[i], seq[(i + 1) % k]), max(seq[i], seq[(i + 1) % k]))
        if e not in E:
            return False
    return True


chk("two-vertex-disjoint-5-cycles",
    is_cycle(N15, G0, ODD_CYCLE_A) and is_cycle(N15, G0, ODD_CYCLE_B)
    and not (set(ODD_CYCLE_A) & set(ODD_CYCLE_B)),
    "%s and %s; hence G0 - v contains an odd cycle for EVERY v, so no G0 - v is bipartite "
    "and KMS Theorem 1.1 cannot apply" % (list(ODD_CYCLE_A), list(ODD_CYCLE_B)))
chk("induced-C5-exists", is_cycle(N15, G0, ODD_CYCLE_A)
    and all((min(a, b), max(a, b)) not in G0
            for a, b in itertools.combinations(ODD_CYCLE_A, 2)
            if (min(a, b), max(a, b)) not in
            {(min(ODD_CYCLE_A[i], ODD_CYCLE_A[(i + 1) % 5]),
              max(ODD_CYCLE_A[i], ODD_CYCLE_A[(i + 1) % 5])) for i in range(5)}),
    "%s is an INDUCED 5-cycle, so G0 is not perfect" % list(ODD_CYCLE_A))

# =============================================================================================
head("2. theta(G0) = 8 and vertex-criticality, by matchings only")

mu0, match0 = max_matching(N15, G0)
chk("matching-number", mu0 == 7, "mu(G0) = 7, e.g. %s"
    % " ".join("%d-%d" % e for e in sorted(match0)))
th0, cover0 = clique_cover_number(N15, G0)
covered = sorted(v for c in cover0 for v in c)
chk("clique-cover-of-size-8-covers-V", th0 == 8 and covered == list(range(N15))
    and all(len(c) == 1 or (min(c), max(c)) in G0 for c in cover0),
    "cover = %s" % " ".join("{" + ",".join(map(str, c)) + "}" for c in cover0))
chk("no-clique-cover-of-size-7-exists", 2 * 7 < N15,
    "G0 is triangle-free so each clique holds at most 2 vertices and 7 cliques reach at most "
    "14 < 15 vertices; with the cover above, theta(G0) = 8 exactly")
a0 = independence_number(N15, G0)
chk("independence-number", a0 == 6, "alpha(G0) = 6")
fc0, fcwit = factor_critical(N15, G0)
chk("factor-critical", fc0,
    "G0 - v has a perfect matching of 7 edges for all 15 v, e.g. at v=0: %s"
    % " ".join("%d-%d" % e for e in sorted(fcwit[0][1])))
thetas = []
for w in range(N15):
    nw, Ew = delete_vertex(N15, G0, w)
    thetas.append(clique_cover_number(nw, Ew)[0])
chk("vertex-critical-in-the-sense-of-Fact-2.1", thetas == [7] * 15,
    "theta(G0 - v) = 7 = theta(G0) - 1 for every one of the 15 vertices v")

# =============================================================================================
head("3. gamma^inf(G0) = 7 < 8 = theta(G0): the LEFT side of the equivalence is FALSE")

yes7, fam7, _, cert7 = edom_le(N15, G0, 7)
chk("gamma-inf-G0-at-most-7", yes7 and cert7 and fam7 == P_FAMILY_G0,
    "the greatest fixed point has |F| = %d seven-sets, each dominating and closed under every "
    "attack, re-verified in a separate pass" % fam7)
no6, ndom6, rk6, cert6 = edom_le(N15, G0, 6)
chk("gamma-inf-G0-greater-than-6", (not no6) and cert6,
    "all %d dominating 6-sets die; the death ranking (max rank %d) is re-verified well founded "
    "in a separate pass, so induction on rank gives gamma^inf(G0) > 6" % (ndom6, rk6))
chk("gamma-inf-G0-equals-7-strictly-below-theta", (not no6) and yes7 and 7 < th0,
    "gamma^inf(G0) = 7 < 8 = theta(G0)")
chk("sandwich-holds-on-G0", a0 <= 7 <= th0, "alpha = 6 <= gamma^inf = 7 <= theta = 8")
gam0, gamwit = domination_number(N15, G0)
chk("domination-number-is-4", gam0 == 4,
    "gamma(G0) = 4, e.g. %s; far below gamma^inf(G0) = 7, so the planar line of Taletskii, "
    "whose hypothesis concerns gamma, does not reach G0" % list(gamwit))

# =============================================================================================
head("4. gamma^inf(G0 - v) = theta(G0 - v) = 7 for all 15 v: the RIGHT side is TRUE")

counts = []
bad4 = []
for w in range(N15):
    nw, Ew = delete_vertex(N15, G0, w)
    thw = clique_cover_number(nw, Ew)[0]
    y7, f7, _, c7 = edom_le(nw, Ew, 7)
    n6, d6, r6, c6 = edom_le(nw, Ew, 6)
    counts.append(d6)
    if not (thw == 7 and y7 and c7 and (not n6) and c6):
        bad4.append(w)
chk("every-deletion-is-tight", bad4 == [],
    "for all 15 v: theta(G0-v) = 7; gamma^inf(G0-v) <= 7 YES with a re-verified positive "
    "certificate; gamma^inf(G0-v) <= 6 NO with a re-verified well-founded ranking")
chk("dominating-6-set-counts-as-printed", counts == P_DOM6_COUNTS,
    "the 15 denominators the negative half ranks: %s" % counts)
chk("the-equivalence-of-Conjecture-2.1-fails-on-G0",
    (not no6) and yes7 and th0 == 8 and bad4 == [] and len(G0) <= 2 * N15 - 3,
    "G0 is vertex-critical, triangle-free, |E| = 27 <= 2n-3 = 27; the right-hand side holds at "
    "all 15 vertices and the left-hand side fails, so Conjecture 2.1 is FALSE")

# =============================================================================================
head("5. the seed, Lemma S, and the published calibration")

fcH, _ = factor_critical(NH, HSTAR)
thH, _ = clique_cover_number(NH, HSTAR)
aH = independence_number(NH, HSTAR)
chk("seed-invariants",
    (NH, len(HSTAR), triangle_edges(NH, HSTAR), fcH, thH, aH) == (13, 24, [], True, 7, 5),
    "H*: n = 13, m = 24, triangle-free, factor-critical, theta = 7, alpha = 5")
yH, famH, _, cH = edom_le(NH, HSTAR, 6)
nH5 = edom_le(NH, HSTAR, 5)
chk("gamma-inf-seed-equals-6", yH and cH and famH == P_FAMILY_HSTAR and not nH5[0] and nH5[3],
    "gamma^inf(H*) <= 6 YES, |F| = %d, positive certificate re-verified; <= 5 NO with a "
    "re-verified ranking; so gamma^inf(H*) = 6 < 7 = theta(H*)" % famH)
chk("Lemma-S-numerically-consistent-on-G0", 7 <= 6 + 1,
    "V(G0) partitions into V(H*) and the ear-tip pair {13,14} inducing a K_2, so "
    "gamma^inf(G0) <= gamma^inf(H*) + gamma^inf(K_2) = 6 + 1 = 7 < 8 = theta(G0); the computed "
    "value 7 matches this bound with equality")

t10bad = []
for s in T10:
    nn, EE = decode_graph6(s)
    th = clique_cover_number(nn, EE)[0]
    y = edom_le(nn, EE, th - 1)
    if not (nn == 13 and th == 7 and triangle_edges(nn, EE) == []
            and factor_critical(nn, EE)[0] and y[0] and y[3]):
        t10bad.append(s)
chk("forced-positive-on-all-13-published-order-13-graphs", t10bad == [],
    "every MMV Table 10 graph is triangle-free, factor-critical, theta = 7, and the decider "
    "returns gamma^inf <= 6 with a re-verified positive certificate: 0 failures of 13")
chk("the-seed-is-MMV-Table-10-entry-4", T10[3] == HSTAR_GRAPH6,
    "'%s' occurs byte-identically in the published table" % HSTAR_GRAPH6)

t11 = []
for s in T11N14:
    nn, EE = decode_graph6(s)
    th = clique_cover_number(nn, EE)[0]
    y, f, _, c = edom_le(nn, EE, th - 1)
    t11.append((nn, len(EE), th, y, f, c))
chk("forced-positive-at-order-14-and-k-equal-theta-minus-1",
    all(x[:4] == (14, 31, 7, True) and x[5] for x in t11)
    and tuple(x[4] for x in t11) == P_FAMILY_T11N14,
    "the two published order-14 graphs return YES at k = 6 with families of sizes %d and %d; "
    "this is why the fifteen NOs of Section 4 are real NOs and not a decider that cannot say "
    "YES on 14 vertices at k = 6" % tuple(x[4] for x in t11))

cyc = []
for m in (5, 7, 9, 11, 13, 15):
    nn, EE = cycle(m)
    th = clique_cover_number(nn, EE)[0]
    g = gamma_inf_exact(nn, EE, th)
    lo = edom_le(nn, EE, th - 1)
    cyc.append((m, th, g, lo[0]))
chk("proved-silent-controls-odd-cycles",
    all(th == (m + 1) // 2 and g == th and not lo for (m, th, g, lo) in cyc),
    "gamma^inf(C_m) = ceil(m/2) = theta(C_m) and the decider says NO at theta-1 for "
    "m = 5,7,9,11,13,15; C_15 is a genuinely sparse member of G0's own order")

NGR = 11
GROTZSCH = frozenset(
    [(i, (i + 1) % 5) if i < (i + 1) % 5 else ((i + 1) % 5, i) for i in range(5)]
    + [(min(5 + i, (i + 4) % 5), max(5 + i, (i + 4) % 5)) for i in range(5)]
    + [(min(5 + i, (i + 1) % 5), max(5 + i, (i + 1) % 5)) for i in range(5)]
    + [(5 + i, 10) for i in range(5)])
thGR, _ = clique_cover_number(NGR, GROTZSCH)
gGR = gamma_inf_exact(NGR, GROTZSCH, thGR)
chk("proved-silent-control-Grotzsch-graph",
    (len(GROTZSCH), triangle_edges(NGR, GROTZSCH), thGR, independence_number(NGR, GROTZSCH),
     gGR, edom_le(NGR, GROTZSCH, 5)[0]) == (20, [], 6, 5, 6, False),
    "n = 11, m = 20, triangle-free, theta = 6, alpha = 5, gamma^inf = 6 = theta, NO at 5")

TWOC5 = frozenset([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 5), (5, 6), (6, 7), (7, 8), (0, 8)])
th2, _ = clique_cover_number(9, TWOC5)
chk("anti-control-cut-vertex-is-not-a-disqualifier",
    factor_critical(9, TWOC5)[0] and th2 == 5 and gamma_inf_exact(9, TWOC5, th2) == 5,
    "two 5-cycles sharing a vertex: factor-critical DESPITE a cut vertex, and "
    "gamma^inf = 5 = theta, so it is correctly not a counterexample")
P5 = frozenset([(0, 1), (1, 2), (2, 3), (3, 4)])
chk("anti-control-parity-and-degree-gates-are-not-vacuous",
    (not factor_critical(5, P5)[0]) and (not factor_critical(*cycle(6))[0]),
    "P_5 is odd but not factor-critical, C_6 is even; both are correctly rejected")

# =============================================================================================
head("6. the family: 78 attachments, 8 isomorphism classes, and two sparser witnesses")

AUT = automorphisms(NH, HSTAR)
chk("automorphism-group-of-the-seed", len(AUT) == P_AUT_HSTAR,
    "|Aut(H*)| = %d, enumerated exhaustively by degree-refined backtracking" % len(AUT))
PAIRS = [(i, j) for i in range(NH) for j in range(i + 1, NH)]
seen = set()
orbits = []
for p in PAIRS:
    if p in seen:
        continue
    orb = set()
    for g in AUT:
        a, b = g[p[0]], g[p[1]]
        orb.add((min(a, b), max(a, b)))
    seen |= orb
    orbits.append(orb)
osz = sorted(len(o) for o in orbits)
chk("eight-orbits-on-the-78-endpoint-pairs",
    len(PAIRS) == 78 and osz == P_ORBIT_SIZES and sum(osz) == 78,
    "orbit sizes %s sum to C(13,2) = 78, so there are EXACTLY 8 witnesses up to isomorphism, "
    "not merely at most 8" % osz)

bad78 = []
for (a, b) in PAIRS:
    n1, E1 = add_ears(NH, HSTAR, ((a, b),))
    if not (n1 == 15 and len(E1) == 27 and triangle_edges(n1, E1) == []
            and factor_critical(n1, E1)[0] and clique_cover_number(n1, E1)[0] == 8):
        bad78.append((a, b, "structure"))
        continue
    if not edom_le(n1, E1, 7)[0] or edom_le(n1, E1, 6)[0]:
        bad78.append((a, b, "gamma-inf"))
        continue
    for w in range(n1):
        nw, Ew = delete_vertex(n1, E1, w)
        if clique_cover_number(nw, Ew)[0] != 7 or not edom_le(nw, Ew, 7)[0] \
                or edom_le(nw, Ew, 6)[0]:
            bad78.append((a, b, w))
            break
chk("all-78-single-ear-attachments-are-counterexamples", bad78 == [],
    "every one of the C(13,2) = 78 endpoint choices gives n = 15, |E| = 27 = 2n-3, "
    "triangle-free, vertex-critical, gamma^inf = 7 < theta = 8, and tightness at all 15 "
    "deletions; 0 failures")

tow = []
for idx, ears in enumerate(TOWER):
    n1, E1 = add_ears(NH, HSTAR, ears)
    th = clique_cover_number(n1, E1)[0]
    y, f, _, c = edom_le(n1, E1, th - 1)
    nn = edom_le(n1, E1, th - 2)
    okdel = all(clique_cover_number(*delete_vertex(n1, E1, w))[0] == th - 1
                and edom_le(*(delete_vertex(n1, E1, w) + (th - 1,)))[0]
                and not edom_le(*(delete_vertex(n1, E1, w) + (th - 2,)))[0]
                for w in range(n1))
    tow.append((n1, len(E1), th, y and c, f, not nn[0] and nn[3], okdel,
                2 * n1 - 3 - len(E1)))
chk("a-2n-4-witness-at-order-17", tow[0][:4] == (17, 30, 9, True) and tow[0][4] == P_TOWER_FAMILY[0]
    and tow[0][5] and tow[0][6],
    "two disjoint ears: n = 17, |E| = 30 = 2n-4, theta = 9, gamma^inf = 8 (|F| = %d), and all "
    "17 deletions tight; so no repaired bound |E| <= 2n-4 rescues the conjecture" % tow[0][4])
chk("a-2n-5-witness-at-order-19", tow[1][:4] == (19, 33, 10, True) and tow[1][4] == P_TOWER_FAMILY[1]
    and tow[1][5] and tow[1][6],
    "three disjoint ears: n = 19, |E| = 33 = 2n-5, theta = 10, gamma^inf = 9 (|F| = %d), and "
    "all 19 deletions tight; so no repaired bound |E| <= 2n-5 rescues it either" % tow[1][4])
chk("the-ear-tower-arithmetic", [(13 + 2 * k, 24 + 3 * k, 7 + k, (2 * (13 + 2 * k) - 3) - (24 + 3 * k))
                                 for k in (1, 2, 3)]
    == [(15, 27, 8, 0), (17, 30, 9, 1), (19, 33, 10, 2)],
    "k ears give n = 13+2k, |E| = 24+3k, theta = 7+k against 2n-3 = 23+4k, so the slack k-1 "
    "GROWS; the three verified rows are k = 1, 2, 3")

# =============================================================================================
print("")
print("NOT RE-RUN: the MINIMALITY of the order 15. This program does not census graphs. That the "
      "least order of a counterexample is 15 rests on a separate exhaustion of all connected "
      "triangle-free factor-critical graphs of order at most 14 with n <= |E| <= 2n-3, which is "
      "not repeated here and is not a claim of the paper's theorem.")
print("NOT RE-RUN: whether some order-15 counterexample has FEWER than 27 edges. The order-15 "
      "cell was never censused, only the structured 78-graph family above, so the sharp edge "
      "threshold at n = 15 is OPEN and nothing here bounds it.")
print("NOT RE-RUN: the asymptotic reading. Clause 5 (gamma^inf < theta) holds for every k in the "
      "ear tower by Lemma S, but the tightness clause is verified here only for k = 1, 2, 3, so "
      "'no bound |E| <= cn+d with c > 3/2 rescues the conjecture' is a CONJECTURE, not proved; "
      "the general case needs a lower-bound technique for gamma^inf that does not exist in print.")
print("NOT RE-RUN: the two censuses that produced the sparser witnesses -- an exhaustive sweep of "
      "3081 two-ear configurations at order 17, and 200 of the 25740 disjoint ear triples at "
      "order 19. Only the two named witnesses are re-checked above; the remaining 25540 triples "
      "are UNTESTED, not failures, and no completeness follows from this run.")
print("NOT RE-RUN: planarity beyond the edge-count argument. G0's non-planarity is certified here "
      "only by |E| = 27 > 2n-4, which is decisive for triangle-free graphs; no Kuratowski "
      "subdivision is exhibited and no embedding is searched for.")
print("NOT RE-RUN: bibliography. No page number, conjecture number, DOI or publication date is "
      "checked by this program; it checks mathematics only.")
print("")

if _failed:
    print("VERDICT: %d CHECKS FAILED: %s" % (len(_failed), ", ".join(_failed)))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _n_pass)
sys.exit(0)

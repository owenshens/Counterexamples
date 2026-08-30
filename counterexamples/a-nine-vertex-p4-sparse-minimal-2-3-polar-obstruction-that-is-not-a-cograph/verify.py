#!/usr/bin/env python3
"""verify.py -- an independent re-derivation of every computational claim in

    "A nine-vertex P_4-sparse minimal (2,3)-polar obstruction that is not a cograph"

Python 3.9+, STANDARD LIBRARY ONLY (itertools, sys, hashlib are all that is used); no external
data file.  Every decision is made in exact integer / bitmask arithmetic: there is no float
anywhere in this program, and no randomness -- every search below is an exhaustion.

WHAT IT READS.  Only the two encodings of the witness that are PRINTED IN THE PAPER: the graph6
string of Section 1 and the labelled 15-edge list of Section 1 (both inside Theorem 1), plus the
nine deletion certificates printed in the table of Step 3 in Section 2.  It re-derives every
quantity the paper asserts
from those printed strings; nothing is imported from the discovery code.

Contract: one `PASS <name> [detail]` line per check, a closing
    VERDICT: ALL <n> CHECKS PASS
and exit status 0 if and only if every check passed.
"""

import itertools
import sys

# ---------------------------------------------------------------------------
# THE OBJECT, EXACTLY AS THE PAPER PRINTS IT
# ---------------------------------------------------------------------------
PAPER_GRAPH6 = "HRKx`co"

# Theorem 1, Section 1 of the paper; vertex order 0..8 = s1 s2 k1 k2 a1 a2 b1 b2 x.
PAPER_LABELS = ["s1", "s2", "k1", "k2", "a1", "a2", "b1", "b2", "x"]
PAPER_EDGES = [
    ("s1", "k1"), ("s2", "k2"), ("k1", "k2"),
    ("k1", "a1"), ("k1", "a2"), ("k1", "b1"), ("k1", "b2"), ("k1", "x"),
    ("k2", "a1"), ("k2", "a2"), ("k2", "b1"), ("k2", "b2"), ("k2", "x"),
    ("a1", "a2"), ("b1", "b2"),
]
PAPER_DEGREES = {"s1": 1, "s2": 1, "k1": 7, "k2": 7,
                 "a1": 3, "a2": 3, "b1": 3, "b2": 3, "x": 2}
PAPER_P4 = ("s1", "k1", "k2", "s2")

# The nine (2,3)-polar partitions of the nine one-vertex deletions, as printed in the paper:
#   deleted vertex -> (A as a set of labels, claimed number of parts of G[A],
#                      B as a set of labels, claimed number of clique components of G[B])
PAPER_CERTIFICATES = {
    "s1": (["k1", "k2", "s2"], 2, ["a1", "a2", "b1", "b2", "x"], 3),
    "s2": (["k1", "k2", "s1"], 2, ["a1", "a2", "b1", "b2", "x"], 3),
    "k1": (["k2", "s2", "x"], 2, ["s1", "a1", "a2", "b1", "b2"], 3),
    "k2": (["k1", "s1", "x"], 2, ["s2", "a1", "a2", "b1", "b2"], 3),
    "x":  (["k1", "k2", "s1"], 2, ["s2", "a1", "a2", "b1", "b2"], 3),
    "a1": (["s1", "s2", "a2", "x"], 1, ["k1", "k2", "b1", "b2"], 1),
    "a2": (["s1", "s2", "a1", "x"], 1, ["k1", "k2", "b1", "b2"], 1),
    "b1": (["s1", "s2", "b2", "x"], 1, ["k1", "k2", "a1", "a2"], 1),
    "b2": (["s1", "s2", "b1", "x"], 1, ["k1", "k2", "a1", "a2"], 1),
}

S_PARAM, K_PARAM = 2, 3          # the cell (s,k) of Problem 1 under attack

# ---------------------------------------------------------------------------
# CHECK HARNESS
# ---------------------------------------------------------------------------
_n_pass = 0
_n_fail = 0


def check(name, cond, detail=""):
    global _n_pass, _n_fail
    if cond:
        _n_pass += 1
        print("PASS %s%s" % (name, (" " + detail) if detail else ""))
    else:
        _n_fail += 1
        print("FAIL %s%s" % (name, (" " + detail) if detail else ""))


# ---------------------------------------------------------------------------
# GRAPHS AS BITMASKS.  A graph is (n, adj) with adj[v] the neighbourhood bitmask of v.
# ---------------------------------------------------------------------------
def make(n, edges):
    adj = [0] * n
    for u, v in edges:
        assert u != v, "no loops"
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return (n, adj)


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def popcount(mask):
    return bin(mask).count("1")


def edge_set(g, mask=None):
    n, adj = g
    if mask is None:
        mask = (1 << n) - 1
    return frozenset((u, v) for u in bits(mask) for v in bits(mask) if u < v and adj[u] >> v & 1)


def degree(g, v, mask=None):
    n, adj = g
    if mask is None:
        mask = (1 << n) - 1
    return popcount(adj[v] & mask & ~(1 << v))


def complement(g):
    n, adj = g
    full = (1 << n) - 1
    return (n, [full & ~adj[v] & ~(1 << v) for v in range(n)])


def components(nbr, mask):
    """Connected components of the graph on `mask` whose neighbourhoods are nbr[]."""
    out = []
    rem = mask
    while rem:
        seed = (rem & -rem).bit_length() - 1
        comp = 0
        stack = [seed]
        while stack:
            u = stack.pop()
            if comp >> u & 1:
                continue
            comp |= 1 << u
            todo = nbr[u] & mask & ~comp
            for w in bits(todo):
                stack.append(w)
        out.append(comp)
        rem &= ~comp
    return out


def is_union_of_cliques(g, mask, kmax):
    """(bool, number of components).  A `k`-cluster in the sense of the paper: a disjoint union
    of at most kmax cliques.  False when the components are not all cliques, whatever kmax."""
    n, adj = g
    comps = components(adj, mask)
    for c in comps:
        for v in bits(c):
            if adj[v] & c != (c & ~(1 << v)):
                return False, len(comps)
    return len(comps) <= kmax, len(comps)


def is_complete_multipartite(g, mask, smax):
    """(bool, number of parts).  Complete multipartite with at most smax parts: the parts are
    the components of the COMPLEMENT restricted to mask; each must be independent in g and
    complete to everything else in mask."""
    n, adj = g
    cnbr = [(~adj[v]) & mask & ~(1 << v) for v in range(n)]
    parts = components(cnbr, mask)
    for p in parts:
        rest = mask & ~p
        for v in bits(p):
            if adj[v] & p:
                return False, len(parts)
            if adj[v] & rest != rest:
                return False, len(parts)
    return len(parts) <= smax, len(parts)


def polar(g, mask, s, k):
    """(bool, A, parts, cliques).  EXHAUSTS all 2^|mask| bipartitions of G[mask]."""
    A = mask
    while True:
        okA, pa = is_complete_multipartite(g, A, s)
        if okA:
            okB, kb = is_union_of_cliques(g, mask & ~A, k)
            if okB:
                return True, A, pa, kb
        if A == 0:
            return False, None, None, None
        A = (A - 1) & mask


def minimal_obstruction(g, mask, s, k):
    """(bool, certificates).  Not (s,k)-polar, but every one-vertex deletion is."""
    if polar(g, mask, s, k)[0]:
        return False, {}
    certs = {}
    for v in bits(mask):
        ok, A, pa, kb = polar(g, mask & ~(1 << v), s, k)
        if not ok:
            return False, {}
        certs[v] = (A, pa, kb)
    return True, certs


def p4s(g, mask=None):
    """Every induced P_4 of G[mask], as a list of 4-subsets (vertex tuples)."""
    n, adj = g
    if mask is None:
        mask = (1 << n) - 1
    vs = list(bits(mask))
    out = []
    for q in itertools.combinations(vs, 4):
        qm = sum(1 << i for i in q)
        deg = sorted(popcount(adj[i] & qm) for i in q)
        m = sum(popcount(adj[i] & qm) for i in q) // 2
        if m != 3 or deg != [1, 1, 2, 2]:
            continue
        if len(components(adj, qm)) != 1:          # connected: rules out anything but a path
            continue
        out.append(q)
    return out


def p4_sparse(g, mask=None):
    """No five vertices induce more than one P_4 (the paper's definition)."""
    n, adj = g
    if mask is None:
        mask = (1 << n) - 1
    vs = list(bits(mask))
    if len(vs) < 5:
        return len(p4s(g, mask)) <= 1
    for q in itertools.combinations(vs, 5):
        if len(p4s(g, sum(1 << i for i in q))) > 1:
            return False
    return True


def is_chordal(g):
    """No induced cycle on four or more vertices."""
    n, adj = g
    for r in range(4, n + 1):
        for q in itertools.combinations(range(n), r):
            qm = sum(1 << i for i in q)
            if all(popcount(adj[i] & qm) == 2 for i in q) and len(components(adj, qm)) == 1:
                return False, q
    return True, None


# ---------------------------------------------------------------------------
# GRAPH6, WRITTEN FROM THE FORMAT DESCRIPTION
# ---------------------------------------------------------------------------
def graph6_decode(s):
    data = [ord(c) - 63 for c in s]
    n = data[0]
    assert 0 <= n <= 62, "this decoder handles n <= 62 only"
    body = data[1:]
    need = (n * (n - 1) // 2 + 5) // 6
    assert len(body) == need, "graph6 body is %d sextets, expected %d" % (len(body), need)
    stream = []
    for x in body:
        assert 0 <= x < 64, "sextet out of range"
        for shift in (5, 4, 3, 2, 1, 0):
            stream.append(x >> shift & 1)
    edges = []
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if stream[idx]:
                edges.append((i, j))
            idx += 1
    return make(n, edges)


def graph6_encode(g):
    n, adj = g
    stream = []
    for j in range(1, n):
        for i in range(j):
            stream.append(1 if adj[i] >> j & 1 else 0)
    while len(stream) % 6:
        stream.append(0)
    out = [chr(63 + n)]
    for b in range(0, len(stream), 6):
        x = 0
        for bit in stream[b:b + 6]:
            x = x * 2 + bit
        out.append(chr(63 + x))
    return "".join(out)


# ---------------------------------------------------------------------------
# NAMED SMALL GRAPHS
# ---------------------------------------------------------------------------
def kn(n):
    return make(n, [(i, j) for i in range(n) for j in range(i + 1, n)])


def path(n):
    return make(n, [(i, i + 1) for i in range(n - 1)])


def cycle(n):
    return make(n, [(i, (i + 1) % n) for i in range(n)])


def matchings(m, extra=0):
    """mK_2 plus `extra` isolated vertices."""
    return make(2 * m + extra, [(2 * i, 2 * i + 1) for i in range(m)])


def disjoint_union(*gs):
    off = 0
    edges = []
    for n, adj in gs:
        for u in range(n):
            for v in range(u + 1, n):
                if adj[u] >> v & 1:
                    edges.append((u + off, v + off))
        off += n
    return make(off, edges)


def thin_spider(head_sizes):
    """The paper's construction: S = {s1,s2} independent, K = {k1,k2} a clique, the two legs
    s1-k1 and s2-k2, K complete and S anticomplete to a head that is a disjoint union of
    cliques of the given sizes.  Vertex order 0=s1 1=s2 2=k1 3=k2 then the head."""
    n = 4 + sum(head_sizes)
    edges = [(0, 2), (1, 3), (2, 3)]
    head = []
    at = 4
    for sz in head_sizes:
        blk = list(range(at, at + sz))
        head.append(blk)
        at += sz
        for i in range(len(blk)):
            for j in range(i + 1, len(blk)):
                edges.append((blk[i], blk[j]))
    for blk in head:
        for v in blk:
            edges.append((2, v))
            edges.append((3, v))
    return make(n, sorted(set(edges))), head


# ===========================================================================
# GROUP A -- THE PREDICATES ARE TESTED BEFORE ANY NUMBER BELOW IS BELIEVED.
# The one error that could fake minimality is a too-permissive cluster test.
# ===========================================================================
print("=== A. predicate self-tests (a too-permissive predicate would fake minimality) ===")
FULL = lambda g: (1 << g[0]) - 1

ok, k = is_union_of_cliques(matchings(2), FULL(matchings(2)), 1)
check("pred-uc-2K2-is-not-one-cluster", (ok, k) == (False, 2), "is_uc(2K2,1) = (False, 2 components)")
ok, k = is_union_of_cliques(kn(4), FULL(kn(4)), 1)
check("pred-uc-K4-is-one-clique", (ok, k) == (True, 1), "is_uc(K4,1) = (True, 1)")
ok, k = is_union_of_cliques(path(3), FULL(path(3)), 3)
check("pred-uc-P3-is-no-union-of-cliques", ok is False, "is_uc(P3,3) = False even with k=3")
ok, k = is_union_of_cliques(matchings(2), FULL(matchings(2)), 2)
check("pred-uc-2K2-is-two-cliques", (ok, k) == (True, 2), "is_uc(2K2,2) = (True, 2)")

ok, p = is_complete_multipartite(matchings(2), FULL(matchings(2)), 2)
check("pred-cm-2K2-is-not-complete-multipartite", ok is False, "is_cm(2K2,2) = False")
ok, p = is_complete_multipartite(cycle(4), FULL(cycle(4)), 2)
check("pred-cm-C4-is-K22", (ok, p) == (True, 2), "is_cm(C4,2) = (True, 2)")
ok, p = is_complete_multipartite(kn(3), FULL(kn(3)), 2)
check("pred-cm-K3-needs-three-parts", (ok, p) == (False, 3), "is_cm(K3,2) = (False, 3 parts)")
ok, p = is_complete_multipartite(kn(3), FULL(kn(3)), 3)
check("pred-cm-K3-is-three-parts", (ok, p) == (True, 3), "is_cm(K3,3) = (True, 3)")
ok, p = is_complete_multipartite(make(3, []), 0b111, 1)
check("pred-cm-empty-graph-is-one-part", (ok, p) == (True, 1), "is_cm(3K1,1) = (True, 1)")
ok, p = is_complete_multipartite(kn(1), 0, 1)
check("pred-cm-empty-set-is-zero-parts", is_complete_multipartite(kn(1), 0, 1) == (True, 0),
      "is_cm(emptyset,1) = (True, 0)")

check("pred-p4-count-on-P4", len(p4s(path(4))) == 1, "P_4 induces exactly 1 P_4")
check("pred-p4-count-on-C4", len(p4s(cycle(4))) == 0, "C_4 is P_4-free")
check("pred-p4-count-on-C5", len(p4s(cycle(5))) == 5, "C_5 induces 5 P_4s (one per deleted vertex)")
check("pred-p4-sparse-rejects-C5", p4_sparse(cycle(5)) is False,
      "C_5 is not P_4-sparse: its one 5-set carries 5 > 1 induced P_4s")
check("pred-p4-sparse-accepts-P4", p4_sparse(path(4)) is True, "P_4 itself is P_4-sparse")

# ===========================================================================
# GROUP B -- THE WITNESS, READ FROM THE PAPER'S TWO PRINTED ENCODINGS
# ===========================================================================
print("=== B. the witness, decoded from the strings printed in the paper ===")
IDX = {name: i for i, name in enumerate(PAPER_LABELS)}
H_from_labels = make(9, [(IDX[u], IDX[v]) for u, v in PAPER_EDGES])
H_from_g6 = graph6_decode(PAPER_GRAPH6)
H = H_from_g6                      # everything below is decided on the graph6 graph

check("g6-order", H_from_g6[0] == 9, "graph6 %r decodes to n = 9" % PAPER_GRAPH6)
check("g6-label-equality-with-printed-edge-list",
      edge_set(H_from_g6) == edge_set(H_from_labels),
      "the 15 edges decoded from graph6 are EXACTLY the 15 labelled edges printed in Theorem 1 "
      "(edge-set equality, not a count)")
check("g6-round-trip", graph6_encode(H_from_g6) == PAPER_GRAPH6,
      "re-encoding the decoded graph reproduces %r" % PAPER_GRAPH6)
check("edge-count", len(edge_set(H)) == 15, "m = 15")
check("non-edge-count", 9 * 8 // 2 - len(edge_set(H)) == 21,
      "21 of the C(9,2) = 36 pairs are non-edges")
check("degrees", {PAPER_LABELS[v]: degree(H, v) for v in range(9)} == PAPER_DEGREES,
      "degree sequence matches the paper: " + repr(PAPER_DEGREES))
check("degree-sum", sum(degree(H, v) for v in range(9)) == 30, "sum of degrees = 30 = 2*15")

S = [IDX["s1"], IDX["s2"]]
Km = [IDX["k1"], IDX["k2"]]
R = [IDX[x] for x in ("a1", "a2", "b1", "b2", "x")]
Smask = sum(1 << v for v in S)
Kmask = sum(1 << v for v in Km)
Rmask = sum(1 << v for v in R)

check("spider-partition-covers-V", Smask | Kmask | Rmask == (1 << 9) - 1
      and Smask & Kmask == 0 and Smask & Rmask == 0 and Kmask & Rmask == 0,
      "S, K, R partition the 9 vertices")
check("spider-S-independent", all(H[1][v] & Smask == 0 for v in S), "S = {s1,s2} is independent")
check("spider-K-clique", H[1][Km[0]] >> Km[1] & 1 == 1, "K = {k1,k2} is a clique")
check("spider-legs-are-a-matching",
      all(popcount(H[1][S[i]] & Kmask) == 1 for i in (0, 1))
      and (H[1][S[0]] & Kmask) != (H[1][S[1]] & Kmask),
      "each s_i has exactly one neighbour in K and they differ: a thin spider (|K| = 2, so the "
      "paper's thin d(s)=1 and thick d(s)=|K|-1 coincide)")
check("spider-K-complete-to-head", all(H[1][v] & Rmask == Rmask for v in Km),
      "K is complete to the head R")
check("spider-S-anticomplete-to-head", all(H[1][v] & Rmask == 0 for v in S),
      "S is anticomplete to the head R")
okR, kR = is_union_of_cliques(H, Rmask, 3)
compsR = sorted(popcount(c) for c in components(H[1], Rmask))
check("head-is-2K2-plus-K1", okR and kR == 3 and compsR == [1, 2, 2],
      "the head R is a disjoint union of cliques of sizes %s, i.e. 2K_2 + K_1" % compsR)
chord, cyc = is_chordal(H)
check("chordal", chord, "no induced cycle on 4..9 vertices (all 2^9 vertex subsets examined)")

Imask = sum(1 << IDX[v] for v in ("s1", "s2", "x"))
Wmask = sum(1 << IDX[v] for v in ("a1", "a2", "b1", "b2"))
okW, kW = is_union_of_cliques(H, Wmask, 2)
check("is-2K2-split",
      all(H[1][v] & Imask == 0 for v in bits(Imask))
      and H[1][Km[0]] >> Km[1] & 1 == 1
      and okW and kW == 2
      and sorted(popcount(c) for c in components(H[1], Wmask)) == [2, 2]
      and Imask | Kmask | Wmask == (1 << 9) - 1,
      "V splits as the clique K = {k1,k2}, the independent set {s1,s2,x} and {a1,a2,b1,b2} "
      "inducing 2K_2: the shape of a 2K_2-split graph")

# ===========================================================================
# GROUP C -- THE FIVE CLAUSES OF THE THEOREM
# ===========================================================================
print("=== C. the theorem's five clauses ===")
allp4 = p4s(H)
check("p4-count-is-one", len(allp4) == 1,
      "exactly 1 induced P_4 among all C(9,4) = %d four-subsets" % len(list(itertools.combinations(range(9), 4))))
check("p4-is-the-printed-one", set(PAPER_LABELS[i] for i in allp4[0]) == set(PAPER_P4),
      "the unique induced P_4 is s1-k1-k2-s2, as printed")
check("p4-sparse", p4_sparse(H) is True,
      "no one of the C(9,5) = %d five-subsets induces two P_4s"
      % len(list(itertools.combinations(range(9), 5))))
check("not-a-cograph", len(allp4) > 0, "H* has an induced P_4, so it is not a cograph")
check("head-not-split", polar(H, Rmask, 1, 1)[0] is False,
      "the head 2K_2 + K_1 is not split, i.e. not (1,1)-polar (all 2^5 = 32 bipartitions tried)")

pol, _A, _p, _k = polar(H, FULL(H), S_PARAM, K_PARAM)
check("not-2-3-polar", pol is False,
      "H* is not (2,3)-polar: all 2^9 = 512 bipartitions of V exhausted, none has G[A] complete "
      "multipartite with <= 2 parts and G[B] a union of <= 3 cliques")
check("not-1-1-polar", polar(H, FULL(H), 1, 1)[0] is False, "H* is not split")

# the nine printed certificates, checked one at a time against the DEFINITION
for name in PAPER_LABELS:
    Alab, want_parts, Blab, want_cliques = PAPER_CERTIFICATES[name]
    v = IDX[name]
    rest = FULL(H) & ~(1 << v)
    Am = sum(1 << IDX[x] for x in Alab)
    Bm = sum(1 << IDX[x] for x in Blab)
    okA, pa = is_complete_multipartite(H, Am, S_PARAM)
    okB, kb = is_union_of_cliques(H, Bm, K_PARAM)
    good = (Am | Bm == rest and Am & Bm == 0 and okA and okB
            and pa == want_parts and kb == want_cliques)
    check("certificate-del-%s" % name, good,
          "A = {%s} complete multipartite, %d part(s) <= 2 | B = {%s} a %d-cluster <= 3 | "
          "A and B partition V - %s" % (",".join(Alab), pa, ",".join(Blab), kb, name))

# and, independently of the printed certificates, by exhaustion
mo, certs = minimal_obstruction(H, FULL(H), S_PARAM, K_PARAM)
check("all-nine-deletions-polar-by-exhaustion", len(certs) == 9,
      "an independent exhaustive search finds a (2,3)-polar partition of each of the nine "
      "one-vertex deletions, without reading the printed certificates")
check("minimal-2-3-polar-obstruction", mo is True,
      "H* is a MINIMAL (2,3)-polar obstruction")
check("problem-1-answered-NO",
      mo and p4_sparse(H) and len(allp4) > 0,
      "P_4-sparse + not a cograph + minimal (2,3)-polar obstruction: the answer to Problem 1 is NO")

# not one of the 50 published P_4-sparse minimal 2-polar obstructions
x = IDX["x"]
check("not-2-2-polar", polar(H, FULL(H), 2, 2)[0] is False, "H* is not (2,2)-polar either")
check("not-a-minimal-2-2-obstruction",
      polar(H, FULL(H) & ~(1 << x), 2, 2)[0] is False,
      "H* - x is not (2,2)-polar, so H* is NOT a minimal (2,2)-polar obstruction and therefore "
      "not one of the 50 published P_4-sparse minimal 2-polar obstructions")

# the complement, giving the cell (3,2) for free
Hc = complement(H)
moc, certsc = minimal_obstruction(Hc, FULL(Hc), K_PARAM, S_PARAM)
check("complement-p4-sparse", p4_sparse(Hc) is True, "the complement of H* is P_4-sparse")
check("complement-not-a-cograph", len(p4s(Hc)) > 0, "the complement of H* is not a cograph")
check("complement-minimal-3-2-obstruction", moc is True,
      "the complement is a minimal (3,2)-polar obstruction, so Problem 1 fails at (3,2) too")

# the two published order bounds
check("order-under-hannebauer-cap", 9 <= (S_PARAM + 1) * (K_PARAM + 1),
      "|V(H*)| = 9 <= (s+1)(k+1) = 12, the bound quoted as Theorem 57 of the target")
check("order-under-2K2-split-bound", 9 <= S_PARAM + 2 * K_PARAM + 2,
      "|V(H*)| = 9 <= s + 2k + 2 = 10, the bound of the authors' 2K_2-split paper")

# ===========================================================================
# GROUP D -- CONTROLS.  A run that could not distinguish a refutation from a
# confirmation would be worthless, so both directions are forced here.
# ===========================================================================
print("=== D. controls ===")
G8 = (FULL(H) & ~(1 << x))
okp, A8, p8, k8 = polar(H, G8, S_PARAM, K_PARAM)
check("ctrl-forced-positive-H-minus-x", okp is True,
      "FORCED POSITIVE at the very cell under test: H* - x IS (2,3)-polar, certificate A = {%s} "
      "(%d parts), B = {%s} (%d cliques)"
      % (",".join(PAPER_LABELS[i] for i in bits(A8)), p8,
         ",".join(PAPER_LABELS[i] for i in bits(G8 & ~A8)), k8))
check("ctrl-4K2-polar", polar(matchings(4), FULL(matchings(4)), S_PARAM, K_PARAM)[0] is True,
      "4K_2 is (2,3)-polar (must be True)")
check("ctrl-5K2-not-polar", polar(matchings(5), FULL(matchings(5)), S_PARAM, K_PARAM)[0] is False,
      "5K_2 is not (2,3)-polar (must be False)")
M9 = matchings(4, 1)
mo9, _c9 = minimal_obstruction(M9, FULL(M9), S_PARAM, K_PARAM)
check("ctrl-detector-can-say-YES", mo9 is True and len(p4s(M9)) == 0,
      "4K_2 + K_1 is a minimal (2,3)-polar obstruction AND a cograph, so a conjecture-CONFIRMING "
      "object and a refuting one are distinguishable by this program")
P4_2K2 = disjoint_union(path(4), matchings(2))
P4_3K2 = disjoint_union(path(4), matchings(3))
check("ctrl-P4-plus-2K2-polar", polar(P4_2K2, FULL(P4_2K2), S_PARAM, K_PARAM)[0] is True,
      "P_4 + 2K_2 (8 vertices) is (2,3)-polar (a reproduced dead end)")
check("ctrl-P4-plus-3K2-not-polar-not-minimal",
      polar(P4_3K2, FULL(P4_3K2), S_PARAM, K_PARAM)[0] is False
      and minimal_obstruction(P4_3K2, FULL(P4_3K2), S_PARAM, K_PARAM)[0] is False,
      "P_4 + 3K_2 (10 vertices) is neither (2,3)-polar nor minimal (a reproduced dead end)")

# someone else's published integer: Foldes and Hammer (1977)
sigs = set()
ngraphs = 0
for n in range(1, 6):
    pairs = list(itertools.combinations(range(n), 2))
    for m in range(1 << len(pairs)):
        ngraphs += 1
        g = make(n, [pairs[i] for i in range(len(pairs)) if m >> i & 1])
        if minimal_obstruction(g, FULL(g), 1, 1)[0]:
            sigs.add((n, len(edge_set(g)), tuple(sorted(degree(g, v) for v in range(n)))))
expect = {(4, 2, (1, 1, 1, 1)), (4, 4, (2, 2, 2, 2)), (5, 5, (2, 2, 2, 2, 2))}
check("published-integer-foldes-hammer-1977", sigs == expect and ngraphs == 1099,
      "over all %d labelled graphs on n <= 5 the minimal (1,1)-polar (= split) obstructions are "
      "exactly 2K_2, C_4, C_5, reproducing Foldes-Hammer" % ngraphs)
check("anti-control-only-C5-separates",
      len(p4s(matchings(2))) == 0 and len(p4s(cycle(4))) == 0 and p4_sparse(cycle(5)) is False,
      "ANTI-CONTROL: 2K_2 and C_4 are P_4-free, so they cannot separate a P_4-sparse test from a "
      "P_4-free one; C_5 is the only separator among the three split obstructions and it is "
      "correctly rejected as not P_4-sparse")

# ===========================================================================
# GROUP E -- THE 16-SHAPE SUB-FAMILY EXHAUSTION OF SECTION 3
# ===========================================================================
print("=== E. the sub-family exhaustion that located H* (Section 3) ===")
shapes = [(a, b, c) for a in range(1, 9) for b in range(1, a + 1) for c in range(1, b + 1)
          if a + b + c <= 8]
check("sixteen-head-shapes", len(shapes) == 16,
      "the partitions of 3..8 into exactly three positive parts number 1+1+2+3+4+5 = 16: %s"
      % (sorted(shapes),))
check("sixteen-shapes-obey-the-cap", all(4 + sum(s) <= 12 for s in shapes),
      "every shape has 4 + a+b+c <= 12 = (s+1)(k+1), so the family is the whole search space "
      "allowed by the cap for |S| = |K| = 2 and a head of exactly three cliques")
nonpolar = []
winners = []
for shp in shapes:
    g, _head = thin_spider(list(shp))
    if not polar(g, FULL(g), S_PARAM, K_PARAM)[0]:
        nonpolar.append(shp)
    if minimal_obstruction(g, FULL(g), S_PARAM, K_PARAM)[0]:
        winners.append(shp)
check("nonpolar-shapes-are-exactly-those-with-two-big-cliques",
      sorted(nonpolar) == sorted(s for s in shapes if s[1] >= 2),
      "the spider is NOT (2,3)-polar for exactly the %d shapes with at least two head cliques of "
      "size >= 2, and IS (2,3)-polar for the %d shapes with b = c = 1 (whose head K_a + 2K_1 is "
      "split): %s" % (len(nonpolar), 16 - len(nonpolar), sorted(nonpolar)))
reasons = {}
for shp in nonpolar:
    if shp == (2, 2, 1):
        continue
    a, b, c = shp
    reduced = tuple(sorted([a - 1, b, c], reverse=True)) if c == 1 else \
        tuple(sorted([a, b, c - 1], reverse=True))
    g, _h = thin_spider(list(shp))
    # the deleted vertex is the last vertex of the a-clique (c == 1) or of the c-clique
    v = (4 + a - 1) if c == 1 else (4 + a + b + c - 1)
    reasons[shp] = (reduced, reduced[1] >= 2,
                    not polar(g, FULL(g) & ~(1 << v), S_PARAM, K_PARAM)[0])
check("every-other-nonpolar-shape-has-a-nonpolar-deletion",
      len(reasons) == 9 and all(ok1 and ok2 for _r, ok1, ok2 in reasons.values()),
      "each of the other %d non-polar shapes loses one head vertex (from the largest clique when "
      "c = 1, from the smallest otherwise) to a shape that again has two cliques of size >= 2 and "
      "is again not (2,3)-polar, so it is not minimal" % len(reasons))
check("exactly-one-winning-shape", winners == [(2, 2, 1)],
      "of the 16 shapes exactly one yields a minimal (2,3)-polar obstruction, namely the head "
      "2K_2 + K_1: %s" % (winners,))
gw, _h = thin_spider([2, 2, 1])
check("winning-shape-is-the-witness", edge_set(gw) == edge_set(H),
      "the surviving shape's graph is edge-set identical to the H* printed in the paper under the "
      "labelling s1,s2,k1,k2,a1,a2,b1,b2,x")

# ---------------------------------------------------------------------------
# SCOPE -- WHAT THIS PROGRAM DOES NOT COVER
# ---------------------------------------------------------------------------
print("=== SCOPE ===")
print("NOT RE-RUN: the census of ALL P_4-sparse graphs on at most 12 vertices at the cell (2,3). "
      "Group E exhausts ONLY the thin spiders with |S| = |K| = 2 whose head is a disjoint union "
      "of exactly three cliques (16 shapes). Nothing here shows that 9 is the LEAST order of a "
      "P_4-sparse non-cograph minimal (2,3)-polar obstruction, and the paper claims no such thing.")
print("NOT RE-RUN: any cell other than (2,3) and its complement-dual (3,2). No claim is made or "
      "checked for general s, k, for k >= 4, for thick spiders with |K| >= 3, or for the "
      "finiteness or classification of the family of such obstructions.")
print("NOT RE-RUN: the list of the 50 published P_4-sparse minimal 2-polar obstructions is not "
      "reproduced here; what is checked is the weaker, self-contained fact that H* is not a "
      "minimal (2,2)-polar obstruction at all, which suffices to place it outside that list.")
print("NOT RE-RUN: nothing in this program is evidence about the literature. Priority, and in "
      "particular whether any earlier source exhibits this graph, is not a computation.")

total = _n_pass + _n_fail
if _n_fail:
    print("VERDICT: %d of %d CHECKS FAILED" % (_n_fail, total))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _n_pass)
sys.exit(0)

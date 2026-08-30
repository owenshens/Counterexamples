#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- re-derives every quantity claimed in

    "A 13-Vertex C_4-Free Graph Whose Worst Orientation Beats the Matching Bound"

from the objects PRINTED IN THAT PAPER and nothing else: the graph6 string
`L??E@agXAp@wdC`, the 22-edge list of Section 2, the matching M and the
independent set of Section 3, and the published `mold` values of the source
paper reproduced in Section 5.

Python 3.9+.  STANDARD LIBRARY ONLY -- no third-party package, no external data
file, no network, no randomness, no floating point in any decision.  All
arithmetic is on Python integers; the exhaustive orientation sweeps are done
with exact integer bitsets.

Definitions used (those of Bousquet, Deschamps, Lehtila and Parreau).  For an
orientation D of a graph G and S a subset of V(G), put I_D(u) = N^-_D(u) cap S.
S is LOCATING-DOMINATING in D when, for every u not in S, I_D(u) is non-empty
and the sets I_D(u), u not in S, are pairwise distinct.  gamma_LD(D) is the
least size of such an S (S = V always works, so it exists), and

    mold(G) = max { gamma_LD(D) : D an orientation of G }.

HOW mold <= 6 IS PROVED EXHAUSTIVELY IN PURE PYTHON.  A naive sweep would test
1716 six-sets against each of the 2^22 = 4194304 orientations, about 10^9
predicate evaluations -- out of reach for CPython.  The program instead builds,
for each candidate six-set S, the EXACT SET of orientations for which S is
locating-dominating, as one Python integer used as a 2^22-bit bitset, and
unions those sets.  The construction is exact, not a heuristic:

  * whether S is locating-dominating depends on the orientation only through
    the traces I(v), v not in S, so only the edges with exactly one endpoint in
    S matter, and those edges are partitioned by their endpoint outside S;
  * hence "I(v) is empty" is a CYLINDER: the set of orientations agreeing with
    one fixed pattern on the |N(v) cap S| edges from v to S, free elsewhere;
  * and "I(u) = I(v)" is a union of cylinders, one for each non-empty
    c contained in N(u) cap N(v) cap S, each fixing the edges from u to S and
    from v to S;
  * a cylinder with a fixed pattern of value `val` on a set F of edge-bits is
    { val + t : t a sum of distinct powers of two outside F }, which is built
    by starting from the single bit `val` and doubling once per free bit:
    X |= X << (1 << p).  Exact, and O(m) big-integer operations.

GOOD(S) is then the all-ones bitset with the empty-trace cylinders and the
collision cylinders removed.  mold(G) is the least k for which the union of
GOOD(S) over all k-sets S is the all-ones bitset.  Section 6 CROSS-CHECKS this
machinery in both polarities against the naive per-orientation predicate: on
every control graph with at most 10 edges the two agree on the value of mold,
and on the 13-vertex witness they agree on a deterministic sample of
orientations at k = 6 (where the answer is YES) and at k = 5 (where it is NO).
"""

import sys
from itertools import combinations

# ----------------------------------------------------------------------------
# check bookkeeping
# ----------------------------------------------------------------------------
_N_PASS = 0
_FAILED = []


def check(name, cond, detail=''):
    global _N_PASS
    if cond:
        _N_PASS += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILED.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def head(t):
    print('--- %s ---' % t)


# ----------------------------------------------------------------------------
# the objects printed in the paper
# ----------------------------------------------------------------------------
GRAPH6 = 'L??E@agXAp@wdC'

# The 22 edges of Section 2, IN THE ORDER PRINTED THERE.  That order is lexicographic,
# and it is what fixes the orientation-index convention: bit i of the index is 1 exactly
# when the i-th edge {a,b}, a < b, is oriented a -> b.
EDGE_TEXT = """0-6 0-8 0-12 1-6 1-9 1-10 2-7 2-8 2-9 3-7 3-10 3-12
              4-8 4-10 4-11 5-9 5-11 5-12 6-11 7-11 8-10 9-12"""

ADJ_TEXT = {          # the adjacency table printed in Section 2 of the paper
    0: [6, 8, 12], 1: [6, 9, 10], 2: [7, 8, 9], 3: [7, 10, 12],
    4: [8, 10, 11], 5: [9, 11, 12], 6: [0, 1, 11], 7: [2, 3, 11],
    8: [0, 2, 4, 10], 9: [1, 2, 5, 12], 10: [1, 3, 4, 8],
    11: [4, 5, 6, 7], 12: [0, 3, 5, 9],
}

MATCHING = [(0, 6), (1, 9), (2, 7), (3, 10), (4, 11), (5, 12)]
INDEP = [0, 1, 2, 3, 4, 5]
TRACES = {6: {0, 1}, 7: {2, 3}, 8: {0, 2, 4}, 9: {1, 2, 5},
          10: {1, 3, 4}, 11: {4, 5}, 12: {0, 3, 5}}
TRIANGLES = [(4, 8, 10), (5, 9, 12)]

N_CLAIM, M_CLAIM = 13, 22
ALPHA_CLAIM, MATCH_CLAIM, MOLD_CLAIM = 6, 6, 6
BAD5_CLAIM = 3192596          # orientations with no locating-dominating 5-set


# ----------------------------------------------------------------------------
# graphs
# ----------------------------------------------------------------------------
def norm(pairs):
    out = []
    for a, b in pairs:
        if a == b:
            raise ValueError('loop')
        out.append((min(a, b), max(a, b)))
    return sorted(set(out))


class Graph(object):
    """A simple graph with a FIXED edge indexing; edge i carries orientation bit i."""

    def __init__(self, n, edges):
        self.n = n
        self.edges = norm(edges)
        self.m = len(self.edges)
        self.adj = [set() for _ in range(n)]
        self.eidx = {}
        for i, (a, b) in enumerate(self.edges):
            self.adj[a].add(b)
            self.adj[b].add(a)
            self.eidx[(a, b)] = i

    def deg(self, v):
        return len(self.adj[v])

    # bit value on the edge {v,w} that puts w INTO I(v), i.e. that orients w -> v.
    # Convention: bit i = 1 means a -> b for edge i = {a,b} with a < b.
    def inbit(self, v, w):
        a, b = (v, w) if v < w else (w, v)
        return self.eidx[(a, b)], (1 if a == w else 0)

    def inmasks(self, o):
        """In-neighbourhood bitmasks of the orientation whose index is o."""
        im = [0] * self.n
        for i, (a, b) in enumerate(self.edges):
            if (o >> i) & 1:
                im[b] |= 1 << a          # a -> b
            else:
                im[a] |= 1 << b          # b -> a
        return im


def parse_edges(text):
    out = []
    for tok in text.split():
        a, b = tok.split('-')
        out.append((int(a), int(b)))
    return out


def graph6_decode(s):
    """Decode a graph6 string.  Returns (n, edges, bits_consumed, padding_bits)."""
    data = [ord(c) - 63 for c in s.strip()]
    for d in data:
        if d < 0 or d > 63:
            raise ValueError('character out of graph6 range')
    n = data[0]
    if n > 62:
        raise ValueError('this decoder handles n <= 62 only')
    body = data[1:]
    need = n * (n - 1) // 2
    bits = []
    for d in body:
        for k in (5, 4, 3, 2, 1, 0):
            bits.append((d >> k) & 1)
    if len(bits) < need:
        raise ValueError('graph6 body too short')
    edges = []
    p = 0
    for j in range(1, n):                # column-major upper triangle
        for i in range(j):
            if bits[p]:
                edges.append((i, j))
            p += 1
    pad = bits[need:]
    return n, norm(edges), need, pad


# ----------------------------------------------------------------------------
# the exact orientation bitsets
# ----------------------------------------------------------------------------
def cylinder(fixed, m):
    """Bitset of { orientations o : (o >> e) & 1 == b for every (e, b) in fixed }.

    `fixed` is a sorted tuple of (edge index, required bit).  Exact: the set is
    val + span of the free bits, and doubling once per free bit enumerates that
    span with no repetition, because distinct subsets of powers of two have
    distinct sums.
    """
    val = 0
    fx = set()
    for e, b in fixed:
        fx.add(e)
        if b:
            val |= 1 << e
    X = 1 << val
    for p in range(m):
        if p not in fx:
            X |= X << (1 << p)
    return X


class Cover(object):
    """Exact `mold` by unioning, over k-sets S, the orientations where S is LD."""

    def __init__(self, g):
        self.g = g
        self.ALL = (1 << (1 << g.m)) - 1
        self._cyl = {}

    def _empty_trace(self, v, nb):
        key = tuple(sorted((self.g.inbit(v, w)[0], 1 - self.g.inbit(v, w)[1]) for w in nb))
        X = self._cyl.get(key)
        if X is None:
            X = cylinder(key, self.g.m)
            self._cyl[key] = X
        return X

    def good(self, S):
        """Exact bitset of the orientations for which S is locating-dominating."""
        g, ALL = self.g, self.ALL
        Ss = set(S)
        out = [v for v in range(g.n) if v not in Ss]
        if not out:
            return ALL
        grp = {}
        for v in out:
            nb = sorted(g.adj[v] & Ss)
            if not nb:                    # I(v) is empty in every orientation
                return 0
            grp[v] = nb
        R = ALL
        for v in out:
            R &= ALL ^ self._empty_trace(v, grp[v])
            if not R:
                return 0
        for u, v in combinations(out, 2):
            com = sorted(set(grp[u]) & set(grp[v]))
            if not com:
                continue                  # I(u) = I(v) would force both empty
            for r in range(1, len(com) + 1):
                for c in combinations(com, r):
                    cs = set(c)
                    it = []
                    for w in grp[u]:
                        e, b = g.inbit(u, w)
                        it.append((e, b if w in cs else 1 - b))
                    for w in grp[v]:
                        e, b = g.inbit(v, w)
                        it.append((e, b if w in cs else 1 - b))
                    R &= ALL ^ cylinder(tuple(sorted(it)), g.m)
                    if not R:
                        return 0
        return R

    def _est(self, S):
        """A deterministic ordering weight: prod over v outside S of 2^|N(v) cap S| - 1."""
        g = self.g
        Ss = set(S)
        p = 1
        for v in range(g.n):
            if v in Ss:
                continue
            d = len(g.adj[v] & Ss)
            if d == 0:
                return 0
            p *= (1 << d) - 1
        return p

    def cover_at(self, k):
        """(union bitset, #k-sets examined, #k-sets that contributed)."""
        g = self.g
        if k >= g.n:
            return self.ALL, 1, 1
        cands = [S for S in combinations(range(g.n), k) if self._est(S) > 0]
        cands.sort(key=lambda S: (-self._est(S), S))
        cov = 0
        used = 0
        for i, S in enumerate(cands):
            gg = self.good(S)
            if gg:
                cov |= gg
                used += 1
            if cov == self.ALL:
                return cov, i + 1, used
        return cov, len(cands), used

    def mold(self):
        for k in range(1, self.g.n + 1):
            cov, _, _ = self.cover_at(k)
            if cov == self.ALL:
                return k
        return self.g.n


# ----------------------------------------------------------------------------
# the naive predicate -- the independent implementation used for cross-checks
# ----------------------------------------------------------------------------
def is_ld(inm, Smask, n):
    seen = []
    for u in range(n):
        if (Smask >> u) & 1:
            continue
        t = inm[u] & Smask
        if t == 0 or t in seen:
            return False
        seen.append(t)
    return True


def has_ld_of_size(g, o, k):
    inm = g.inmasks(o)
    for S in combinations(range(g.n), k):
        Sm = 0
        for v in S:
            Sm |= 1 << v
        if is_ld(inm, Sm, g.n):
            return True, S
    return False, None


def gamma_ld_naive(g, o):
    for k in range(1, g.n + 1):
        ok, _ = has_ld_of_size(g, o, k)
        if ok:
            return k
    return g.n


def mold_naive(g):
    return max(gamma_ld_naive(g, o) for o in range(1 << g.m))


# ----------------------------------------------------------------------------
# small combinatorial quantities, exhaustively
# ----------------------------------------------------------------------------
def independence_number(g):
    best, arg = 0, ()
    for mask in range(1 << g.n):
        vs = [v for v in range(g.n) if (mask >> v) & 1]
        if len(vs) <= best:
            continue
        if all(b not in g.adj[a] for a, b in combinations(vs, 2)):
            best, arg = len(vs), tuple(vs)
    return best, arg


def matching_number(g):
    best = [0, ()]

    def rec(i, used, chosen):
        if len(chosen) > best[0]:
            best[0], best[1] = len(chosen), tuple(chosen)
        if i == g.m:
            return
        if len(chosen) + (g.m - i) <= best[0]:
            return
        a, b = g.edges[i]
        if not ((used >> a) & 1) and not ((used >> b) & 1):
            chosen.append((a, b))
            rec(i + 1, used | (1 << a) | (1 << b), chosen)
            chosen.pop()
        rec(i + 1, used, chosen)

    rec(0, 0, [])
    return best[0], best[1]


def count_c4_subgraphs(g):
    """Count 4-cycles by looking at every 4-subset of vertices (no codegree shortcut)."""
    tot = 0
    for q in combinations(range(g.n), 4):
        for perm in ((q[0], q[1], q[2], q[3]), (q[0], q[1], q[3], q[2]), (q[0], q[2], q[1], q[3])):
            a, b, c, d = perm
            if b in g.adj[a] and c in g.adj[b] and d in g.adj[c] and a in g.adj[d]:
                tot += 1
    return tot


def girth(g):
    best = None
    for s in range(g.n):
        dist = {s: 0}
        par = {s: None}
        q = [s]
        while q:
            nq = []
            for v in q:
                for w in g.adj[v]:
                    if w not in dist:
                        dist[w] = dist[v] + 1
                        par[w] = v
                        nq.append(w)
                    elif par[v] != w:
                        c = dist[v] + dist[w] + 1
                        if best is None or c < best:
                            best = c
            q = nq
    return best


def named(kind, n):
    if kind == 'C':
        return Graph(n, [(i, (i + 1) % n) for i in range(n)])
    if kind == 'P':
        return Graph(n, [(i, i + 1) for i in range(n - 1)])
    if kind == 'K':
        return Graph(n, [(i, j) for i in range(n) for j in range(i + 1, n)])
    raise ValueError(kind)


PETERSEN = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
            (5, 7), (7, 9), (9, 6), (6, 8), (8, 5)]
BOWTIE = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]
STAR14 = [(0, 1), (0, 2), (0, 3), (0, 4)]


# ============================================================================
def main():
    edges = parse_edges(EDGE_TEXT)
    G = Graph(N_CLAIM, edges)

    head('1. the graph G, exactly as printed in Section 2')
    check('edge-list-parses-to-22-distinct-edges-on-13-vertices',
          G.m == M_CLAIM and G.n == N_CLAIM and len(set(norm(edges))) == len(edges)
          and sorted(set([v for e in edges for v in e])) == list(range(13)),
          'n = %d, m = %d' % (G.n, G.m))
    check('the-printed-edge-order-is-the-one-the-orientation-index-uses',
          [tuple(e) for e in norm(edges)] == [tuple(e) for e in edges] == G.edges,
          'the list of Section 2 is already lexicographic, so bit i of an orientation index '
          'refers to edge %s for i = 0 and to edge %s for i = 21'
          % (str(G.edges[0]), str(G.edges[21])))

    n6, e6, bits, pad = graph6_decode(GRAPH6)
    check('graph6-decodes-LABEL-EQUAL-to-the-printed-edge-list',
          n6 == G.n and e6 == G.edges,
          "'%s' <-> the 22 printed edges, symmetric difference %d"
          % (GRAPH6, len(set(e6) ^ set(G.edges))))
    check('graph6-consumes-78-bits-with-zero-padding-bits',
          bits == 78 and all(b == 0 for b in pad),
          'bits = %d, padding = %s' % (bits, ''.join(str(b) for b in pad) or 'none'))

    check('adjacency-table-of-Section-2-is-reproduced',
          all(sorted(G.adj[v]) == sorted(ADJ_TEXT[v]) for v in range(G.n)),
          'all 13 rows agree')
    degs = [G.deg(v) for v in range(G.n)]
    check('degree-sequence-is-eight-cubic-and-five-quartic-vertices',
          sorted(degs) == [3] * 8 + [4] * 5 and sum(degs) == 2 * G.m,
          'sum of degrees %d = 2m, degrees %s' % (sum(degs), degs))
    check('G-is-connected',
          len(_reach(G, 0)) == G.n, 'one component')

    head('2. G has no C_4 SUBGRAPH (induced or not)')
    codeg = {}
    for u, v in combinations(range(G.n), 2):
        codeg[(u, v)] = len(G.adj[u] & G.adj[v])
    worst = max(codeg.values())
    check('all-78-unordered-vertex-pairs-have-at-most-one-common-neighbour',
          len(codeg) == 78 and worst <= 1,
          '%d pairs tested (not just the %d adjacent ones), maximum codegree %d'
          % (len(codeg), G.m, worst))
    check('no-4-subset-of-vertices-carries-a-4-cycle',
          count_c4_subgraphs(G) == 0,
          'independent route: all C(13,4) = %d quadruples, 4-cycles found 0'
          % len(list(combinations(range(13), 4))))
    p3 = sum(d * (d - 1) // 2 for d in degs)
    check('the-two-C4-routes-are-numerically-consistent',
          p3 == sum(1 for c in codeg.values() if c == 1)
          and sum(c * (c - 1) // 2 for c in codeg.values()) == 0,
          'paths of length two %d = pairs of codegree one; sum of C(codeg,2) = 0' % p3)
    tri = [t for t in combinations(range(G.n), 3)
           if t[1] in G.adj[t[0]] and t[2] in G.adj[t[0]] and t[2] in G.adj[t[1]]]
    check('exactly-the-two-printed-triangles',
          sorted(tri) == sorted(tuple(sorted(t)) for t in TRIANGLES),
          '%s -- a remark, NOT the reason G has no C_4' % (sorted(tri),))
    check('girth-is-three-so-G-lies-outside-the-girth-at-least-five-literature',
          girth(G) == 3, 'girth = 3')

    head("3. n - alpha'(G) = 7, and alpha(G) = 6")
    ok_m = (len(MATCHING) == 6 and all(tuple(sorted(e)) in G.eidx for e in MATCHING)
            and len(set([v for e in MATCHING for v in e])) == 12)
    check('the-printed-6-edge-matching-is-a-matching-of-G',
          ok_m, '%s -- twelve distinct endpoints' % (MATCHING,))
    am, arg = matching_number(G)
    check('matching-number-is-exactly-six',
          am == MATCH_CLAIM and G.n % 2 == 1 and am <= G.n // 2,
          "alpha' = %d by exhaustive search, and n = 13 is odd so alpha' <= 6; witness %s"
          % (am, list(arg)))
    check('n-minus-the-matching-number-is-seven',
          G.n - am == 7, "13 - 6 = 7 is the Lemma-34 upper bound at issue")
    check('the-printed-6-set-is-independent',
          all(b not in G.adj[a] for a, b in combinations(INDEP, 2)),
          '%s -- all fifteen pairs are non-edges' % (INDEP,))
    al, alarg = independence_number(G)
    check('independence-number-is-exactly-six',
          al == ALPHA_CLAIM,
          'alpha = %d by exhaustive search over all 2^13 = 8192 subsets; witness %s'
          % (al, list(alarg)))
    check('G-is-NOT-a-Koenig-Egervary-graph',
          al + am != G.n,
          "alpha + alpha' = %d != 13 = n, so mold >= alpha does not force attainment" % (al + am))

    head('4. mold(G) >= 6, from one explicit orientation')
    cross = [i for i, (a, b) in enumerate(G.edges) if (a in INDEP) != (b in INDEP)]
    inner = [i for i, (a, b) in enumerate(G.edges) if a not in INDEP and b not in INDEP]
    check('the-independent-6-set-meets-18-edges-and-4-edges-avoid-it',
          len(cross) == 18 and len(inner) == 4 and len(cross) + len(inner) == G.m,
          'cross edges %d, edges inside {6,...,12} %s'
          % (len(cross), [G.edges[i] for i in inner]))

    # D* : every cross edge directed AWAY from INDEP.  Edge i = (a,b), a<b, bit 1 means a->b.
    base = 0
    for i in cross:
        a, b = G.edges[i]
        if a in INDEP:
            base |= 1 << i               # a -> b, away from INDEP
        # else b in INDEP and bit 0 means b -> a, already away
    Sm = 0
    for v in INDEP:
        Sm |= 1 << v

    all_source = True
    all_six = True
    traces_ok = True
    for extra in range(1 << len(inner)):
        o = base
        for j, i in enumerate(inner):
            if (extra >> j) & 1:
                o |= 1 << i
        inm = G.inmasks(o)
        if any(inm[v] != 0 for v in INDEP):
            all_source = False
        if not is_ld(inm, Sm, G.n):
            all_six = False
        tr = {v: set(w for w in range(G.n) if (inm[v] >> w) & 1 and w in set(INDEP))
              for v in range(6, 13)}
        if tr != TRACES:
            traces_ok = False
        for k in range(0, 6):
            ok, S = has_ld_of_size(G, o, k)
            if ok:
                all_six = False
    check('in-D-star-every-vertex-of-the-independent-6-set-is-a-source',
          all_source,
          'checked for all %d orientations of the four inner edges' % (1 << len(inner)))
    check('the-printed-traces-are-the-traces-of-D-star',
          traces_ok, '%s' % ({v: sorted(TRACES[v]) for v in sorted(TRACES)},))
    check('gamma_LD-of-D-star-is-exactly-six',
          all_six,
          'the independent 6-set is locating-dominating, and NO set of size 0..5 is, '
          'exhaustively over all C(13,0)+...+C(13,5) = %d subsets, for each of the %d '
          'orientations of the inner edges'
          % (sum(len(list(combinations(range(13), k))) for k in range(6)), 1 << len(inner)))
    check('hence-mold-is-at-least-six',
          all_six and all_source, 'mold(G) >= gamma_LD(D*) = 6')

    head('5. mold(G) <= 6: ALL 2^22 = 4194304 orientations, exhaustively')
    cov6 = Cover(G)
    bits6, seen6, used6 = cov6.cover_at(6)
    total = 1 << G.m
    pc6 = bin(bits6).count('1')
    check('every-orientation-of-G-admits-a-locating-dominating-set-of-size-six',
          pc6 == total,
          'union over six-sets covers %d of %d = 2^22 orientations; %d of the 1716 six-sets '
          'examined, %d of them contributed' % (pc6, total, seen6, used6))
    check('so-gamma_LD-is-at-most-six-for-every-orientation',
          pc6 == total and total - pc6 == 0,
          'orientations with NO locating-dominating 6-set: 0')

    cov5 = Cover(G)
    bits5, seen5, used5 = cov5.cover_at(5)
    pc5 = bin(bits5).count('1')
    bad5 = total - pc5
    check('exactly-3192596-orientations-admit-no-locating-dominating-5-set',
          bad5 == BAD5_CLAIM,
          '%d of %d orientations admit a 5-set, so bad(5) = %d; this is the informative '
          'non-zero on the same code path that returns 0 at k = 6' % (pc5, total, bad5))
    nb = (~bits5) & cov5.ALL
    first_bad = (nb & -nb).bit_length() - 1
    okd = any(has_ld_of_size(G, first_bad, k)[0] for k in range(6))
    check('the-least-such-orientation-is-confirmed-by-the-naive-predicate',
          (not okd),
          'orientation index %d, in the convention fixed above (bit i = 1 means a -> b for '
          'the i-th edge of the Section 2 list): the naive scan over all 2380 sets of size '
          'AT MOST five finds none, so gamma_LD of that orientation is 6 without appealing '
          'to monotonicity' % first_bad)
    mono_pairs = 0
    mono_bad = 0
    for o in range(0, total, 9973):
        inm = G.inmasks(o)
        for S in combinations(range(G.n), 5):
            Sm = 0
            for v in S:
                Sm |= 1 << v
            if not is_ld(inm, Sm, G.n):
                continue
            for w in range(G.n):
                if (Sm >> w) & 1:
                    continue
                mono_pairs += 1
                if not is_ld(inm, Sm | (1 << w), G.n):
                    mono_bad += 1
    check('monotonicity-of-locating-domination-holds-on-every-instance-sampled',
          mono_bad == 0 and mono_pairs > 1000,
          '%d instances (S of size 5 locating-dominating, S union {w} of size 6) over %d '
          'orientations at a fixed stride: 0 violations of the lemma that a superset of a '
          'locating-dominating set is locating-dominating'
          % (mono_pairs, len(range(0, total, 9973))))
    check('mold-of-G-is-exactly-six',
          pc6 == total and bad5 > 0 and cov6.mold() == MOLD_CLAIM,
          'mold(G) = %d' % cov6.mold())

    head('6. the bitset machinery cross-checked against the naive predicate')
    sample = sorted(set(list(range(0, 1024)) + list(range(0, total, 4096)) + [total - 1]))
    bad = 0
    for o in sample:
        ok, _ = has_ld_of_size(G, o, 6)
        if not (ok and ((bits6 >> o) & 1)):
            bad += 1
    check('positive-polarity-sample-agrees-with-the-exhaustive-bitset',
          bad == 0,
          '%d orientations (0..1023, every 4096th, and the last): the naive scan finds a '
          'locating-dominating 6-set in each, exactly as the bitset says' % len(sample))
    negs = []
    o = 0
    while len(negs) < 256 and o < total:
        if not ((bits5 >> o) & 1):
            negs.append(o)
        o += 1319                        # a fixed stride, not a random choice
    bad = 0
    for o in negs:
        ok, _ = has_ld_of_size(G, o, 5)
        if ok:
            bad += 1
    check('negative-polarity-sample-agrees-with-the-exhaustive-bitset',
          bad == 0 and len(negs) == 256,
          '%d orientations the bitset calls bad at k = 5: the naive scan over all 1287 '
          'five-sets confirms NO locating-dominating 5-set in any of them' % len(negs))

    head('7. published mold values reproduced, both polarities')
    controls = []
    for n in (3, 4, 5, 6, 7, 8, 9, 10):
        controls.append(('C_%d' % n, named('C', n), 3 if n == 4 else (n + 1) // 2,
                         'Corollary 36 (n = 3 or n >= 5); C_4 printed separately as 3'))
    for n in (2, 3, 4, 5, 7, 9):
        controls.append(('P_%d' % n, named('P', n), (n + 1) // 2, 'paths'))
    for n in (3, 4, 5, 6, 7):
        controls.append(('K_%d' % n, named('K', n), (n + 1) // 2, 'complete graphs'))
    controls.append(('K_{1,4}', Graph(5, STAR14), 4, 'stars: mold = n - 1'))
    controls.append(('bowtie', Graph(5, BOWTIE), 3, 'C_4-free, attains n - alpha^prime'))
    controls.append(('Petersen', Graph(10, PETERSEN), 5,
                     'C_4-free and non-bipartite: the negative branch inside the class'))
    agree = []
    for name, g, expect, _why in controls:
        got = Cover(g).mold()
        nat = mold_naive(g) if g.m <= 10 else None
        check('mold-of-%s-is-%d-as-published' % (name.replace('_', '').replace('{', '').replace('}', '').replace(',', '-'), expect),
              got == expect and (nat is None or nat == expect),
              'n = %d, m = %d, mold = %d%s' % (g.n, g.m, got,
                                               '' if nat is None else
                                               ', naive per-orientation sweep over all 2^%d '
                                               'orientations agrees' % g.m))
        if nat is not None:
            agree.append(name)
    check('the-two-independent-implementations-agree-on-every-small-control',
          len(agree) >= 14,
          '%d control graphs with at most 10 edges were decided twice, by the bitset cover and '
          'by the naive per-orientation sweep, with no disagreement' % len(agree))

    c4 = named('C', 4)
    check('mold-of-C4-exceeds-n-minus-its-matching-number',
          Cover(c4).mold() == 3 and c4.n - matching_number(c4)[0] == 2,
          "mold(C_4) = 3 > 2 = n - alpha'(C_4): this alone rules out a "
          'minimum-over-orientations misreading of mold')
    pet = Graph(10, PETERSEN)
    pcod = max(len(pet.adj[u] & pet.adj[v]) for u, v in combinations(range(10), 2))
    check('Petersen-is-in-the-same-class-and-ATTAINS-the-bound',
          pcod <= 1 and Cover(pet).mold() == 5 and pet.n - matching_number(pet)[0] == 5,
          "no C_4 subgraph (max codegree %d over all 45 pairs) and mold = 5 = n - alpha', so "
          'the decider says NO where it must' % pcod)

    head('8. the answer to the problem')
    check('G-has-no-C4-subgraph-and-mold-is-strictly-below-n-minus-the-matching-number',
          worst <= 1 and count_c4_subgraphs(G) == 0 and cov6.mold() == 6 and G.n - am == 7,
          "mold(G) = 6 < 7 = n - alpha'(G) on a graph with no 4-cycle")

    print('')
    print('NOT RE-RUN: minimality of any kind. This program says nothing about whether 13 is the '
          'least order, or 22 the least size, of a graph with these properties, and it does not '
          'reproduce the census that suggested so; that census found 32 such graphs at n = 13, '
          'fourteen of them at m = 22, and its per-graph decisions were never re-run by a second '
          'implementation. The paper claims existence only, which needs none of it.')
    print('NOT RE-RUN: the other thirteen graphs at (n, m) = (13, 22), and the two at m = 24. Only '
          'the single witness printed in the paper is decided here.')
    print('NOT RE-RUN: n >= 14. No graph of order 14 or more is examined, so nothing here bounds '
          'how common such graphs are.')
    print('NOT RE-RUN: the bibliographic locators. The line numbers, the page, the problem number '
          '"Open problem 37" and the DOI printed in the paper are not checked by this program; it '
          'checks mathematics only.')
    print('NOT RE-RUN: prior art. Nothing here bears on whether the answer was known.')
    print("NOT RE-RUN: ld(G) and the proof of the source's Lemma 34. The chain ld <= mold <= "
          "n - alpha' is not re-derived; only the two ends that the strict inequality needs, "
          'mold(G) = 6 and n - alpha\'(G) = 7, are computed.')
    print('')
    if _FAILED:
        print('VERDICT: %d OF %d CHECKS FAILED: %s'
              % (len(_FAILED), _N_PASS + len(_FAILED), ', '.join(_FAILED)))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _N_PASS)
    return 0


def _reach(g, s):
    seen = {s}
    st = [s]
    while st:
        v = st.pop()
        for w in g.adj[v]:
            if w not in seen:
                seen.add(w)
                st.append(w)
    return seen


if __name__ == '__main__':
    sys.exit(main())

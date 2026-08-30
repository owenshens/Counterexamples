#!/usr/bin/env python3
"""verify.py -- checks every computational claim of

    "An Affirmative Answer to the Levit-Kadrawi Multi-Break Unicyclic Closure Problem"

Python 3.9+, STANDARD LIBRARY ONLY (itertools, math.comb, sys), no external data file.
All arithmetic is exact integer arithmetic; no float ever decides anything.

WHAT IT READS FROM THE PAPER.  Every object is rebuilt from the parameters printed in
section 2 of the paper -- the pattern-graph calculus (P_t, Z_k, (G:H)_k^(m), S_{2,t},
T_{1,l}, T_{G,m,t}) together with the labelling rule stated there -- and every integer
the paper prints is repeated below as a literal and compared against the recomputation.
Nothing is read from a data file, and the program does not consult the run that produced
the result.

HOW IT COMPUTES.  Independence polynomials of the six unicyclic witnesses are computed
TWICE, by the two independent identities of Lemma 3 of the paper:

    (A)  I(H) = I(H - uv) - x^2 * I((H - uv) - N[u] - N[v])          uv a cycle edge
    (B)  I(H) = I(H - w)  + x   * I(H - N[w])                        w a cycle vertex

and both are required to agree.  Forest polynomials come from the rooted recursion, and
that recursion is itself checked against brute-force subset enumeration on every forest
and every unicyclic graph of order at most 6.

Each check prints one `PASS <name> [detail]` line; the program closes with
`VERDICT: ALL <n> CHECKS PASS` and exits 0 iff every check passed.  It also prints, in
its own output, the list of claims it does NOT re-derive (`NOT RE-RUN: ...`).
"""

import itertools
import sys
from math import comb

# ---------------------------------------------------------------------------
# exact polynomial arithmetic over Z, coefficient lists, index = degree
# ---------------------------------------------------------------------------


def padd(a, b):
    if len(a) < len(b):
        a, b = b, a
    r = list(a)
    for i, v in enumerate(b):
        r[i] += v
    return r


def psub(a, b):
    return padd(a, [-c for c in b])


def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    r[i + j] += ai * bj
    return r


def pmul_all(ps):
    """Balanced product, so intermediate degrees stay small."""
    ps = list(ps)
    if not ps:
        return [1]
    while len(ps) > 1:
        nxt = [pmul(ps[i], ps[i + 1]) for i in range(0, len(ps) - 1, 2)]
        if len(ps) % 2:
            nxt.append(ps[-1])
        ps = nxt
    return ps[0]


def pstrip(p):
    """Drop trailing zero coefficients; [0] stays [0]."""
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


# ---------------------------------------------------------------------------
# the pattern-graph calculus of section 2, with the labelling rule of section 2
#
# A rooted graph is the triple (edge list, order n, root), vertices 0..n-1.
# Vertices are created in the order the recipe creates them, which is exactly the
# labelling the paper fixes.
# ---------------------------------------------------------------------------


def P(t):
    """Path on t vertices, rooted at a leaf; labelled 0..t-1 along the path."""
    return ([(i, i + 1) for i in range(t - 1)], t, 0)


def Z(H, k):
    """New vertex v_0 = 0, then k disjoint copies of H, v_0 joined to each root."""
    E, n, r = H
    edges = []
    N = 1
    for _ in range(k):
        edges += [(a + N, b + N) for a, b in E]
        edges.append((0, r + N))
        N += n
    return (edges, N, 0)


def patt(G, H, k):
    """(G:H)_k -- disjoint union of G and Z_k(H) plus the edge r(G) v_0; root stays r(G)."""
    GE, gn, gr = G
    ZE, zn, _zr = Z(H, k)
    edges = list(GE) + [(a + gn, b + gn) for a, b in ZE] + [(gr, gn)]
    return (edges, gn + zn, gr)


def patt_m(G, H, k, m):
    """(G:H)_k^(m) -- apply (.:H)_k to G m times, root unchanged."""
    cur = G
    for _ in range(m):
        cur = patt(cur, H, k)
    return cur


def S2(t):
    """S_{2,t} = (P_1:P_1)_1^(t): a centre with t pendant paths of length two."""
    return patt_m(P(1), P(1), 1, t)


def T1(l):
    """T_{1,l} = (P_1:P_2)_l^(1)."""
    return patt_m(P(1), P(2), l, 1)


def TG(m, t):
    """T_{G,m,t} = GT_{m,t} = (P_2:S_{2,t})_3^(m)."""
    return patt_m(P(2), S2(t), 3, m)


def KL(m, n):
    """Kadrawi-Levit T_{3,m,n}: a root with three neighbours carrying 3, m, n pendant
    paths of length two.  Order 10 + 2m + 2n."""
    edges = []
    nid = 1
    for a in (3, m, n):
        c = nid
        nid += 1
        edges.append((0, c))
        for _ in range(a):
            edges.append((c, nid))
            edges.append((nid, nid + 1))
            nid += 2
    return (edges, nid, 0)


def cycle(n):
    return ([(i, (i + 1) % n) for i in range(n)], n, 0)


def adjacency(E, n):
    a = [set() for _ in range(n)]
    for u, v in E:
        a[u].add(v)
        a[v].add(u)
    return a


def add_edge(G, u, v):
    E, n, r = G
    return (list(E) + [(u, v)], n, r)


# ---------------------------------------------------------------------------
# independence polynomials
# ---------------------------------------------------------------------------


def forest_poly(adj, verts):
    """I(F;x) for the subgraph of `adj` induced on `verts`, which MUST be a forest.
    Rooted recursion of Lemma 3(i)."""
    vs = set(verts)
    res = [1]
    seen = set()
    for s in verts:
        if s in seen:
            continue
        par = {s: -1}
        order = [s]
        stack = [s]
        seen.add(s)
        while stack:
            v = stack.pop()
            for u in adj[v]:
                if u in vs and u not in par:
                    par[u] = v
                    order.append(u)
                    stack.append(u)
                    seen.add(u)
        A = {}
        B = {}
        for v in reversed(order):
            kids = [u for u in adj[v] if u in vs and par.get(u) == v]
            A[v] = pmul([0, 1], pmul_all([B[u] for u in kids]))
            B[v] = pmul_all([padd(A[u], B[u]) for u in kids])
        res = pmul(res, padd(A[s], B[s]))
    return res


def is_forest(E, n):
    adj = adjacency(E, n)
    seen = set()
    for s in range(n):
        if s in seen:
            continue
        par = {s: -1}
        order = [s]
        stack = [s]
        seen.add(s)
        while stack:
            v = stack.pop()
            for u in adj[v]:
                if u not in par:
                    par[u] = v
                    order.append(u)
                    stack.append(u)
                    seen.add(u)
                elif par[v] != u:
                    return False
        # component is a tree iff #edges = #vertices - 1
        m = sum(len(adj[v]) for v in order) // 2
        if m != len(order) - 1:
            return False
    return True


def connected(E, n):
    adj = adjacency(E, n)
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for u in adj[v]:
            if u not in seen:
                seen.add(u)
                stack.append(u)
    return len(seen) == n


def unique_cycle(E, n):
    """The vertex set of the unique cycle of a unicyclic graph, by peeling leaves."""
    adj = adjacency(E, n)
    deg = [len(adj[v]) for v in range(n)]
    alive = set(range(n))
    changed = True
    while changed:
        changed = False
        for v in list(alive):
            if deg[v] <= 1:
                alive.discard(v)
                changed = True
                for u in adj[v]:
                    if u in alive:
                        deg[u] -= 1
    return alive


def poly_route_A(E, n, u, v):
    """Identity (A): delete the cycle edge uv."""
    TE = [e for e in E if set(e) != {u, v}]
    adjT = adjacency(TE, n)
    whole = forest_poly(adjT, list(range(n)))
    drop = {u, v} | adjT[u] | adjT[v]
    rest = [w for w in range(n) if w not in drop]
    inner = forest_poly(adjT, rest) if rest else [1]
    return pstrip(psub(whole, pmul([0, 0, 1], inner)))


def poly_route_B(E, n, w):
    """Identity (B): split on the cycle vertex w."""
    adj = adjacency(E, n)
    rest = [z for z in range(n) if z != w]
    closed = {w} | adj[w]
    rest2 = [z for z in range(n) if z not in closed]
    p1 = forest_poly(adj, rest) if rest else [1]
    p2 = forest_poly(adj, rest2) if rest2 else [1]
    return pstrip(padd(p1, pmul([0, 1], p2)))


def brute_poly(E, n):
    """I(G;x) by enumerating all 2^n subsets.  Only used on n <= 6 graphs."""
    adj = adjacency(E, n)
    c = [0] * (n + 1)
    for r in range(n + 1):
        for S in itertools.combinations(range(n), r):
            s = set(S)
            if all(not (adj[x] & s) for x in S):
                c[r] += 1
    return pstrip(c)


# ---------------------------------------------------------------------------
# the predicates of the paper
# ---------------------------------------------------------------------------


def breaks(c):
    """{k : i_k^2 < i_{k-1} i_{k+1}}, 1 <= k <= alpha-1.  Exact integers."""
    return [k for k in range(1, len(c) - 1) if c[k] * c[k] < c[k - 1] * c[k + 1]]


def is_unimodal(c):
    i = 0
    while i + 1 < len(c) and c[i + 1] >= c[i]:
        i += 1
    while i + 1 < len(c) and c[i + 1] <= c[i]:
        i += 1
    return i == len(c) - 1


def unique_mode(c):
    """(index, True) if the maximum is attained exactly once."""
    m = max(c)
    return c.index(m), c.count(m) == 1


def longest_run(bs):
    if not bs:
        return 0
    best = cur = 1
    for a, b in zip(bs, bs[1:]):
        cur = cur + 1 if b == a + 1 else 1
        best = max(best, cur)
    return best


def nonedges(E, n):
    adj = adjacency(E, n)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if j not in adj[i]]


# ---------------------------------------------------------------------------
# graph6
# ---------------------------------------------------------------------------


def g6_decode(s):
    """-> (n, set of edges).  Handles both the n<=62 and the extended header."""
    if ord(s[0]) == 126:
        n = 0
        for ch in s[1:4]:
            n = (n << 6) | (ord(ch) - 63)
        rest = s[4:]
    else:
        n = ord(s[0]) - 63
        rest = s[1:]
    bits = []
    for ch in rest:
        v = ord(ch) - 63
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    E = set()
    p = 0
    for j in range(1, n):
        for i in range(j):
            if bits[p]:
                E.add((i, j))
            p += 1
    return n, E


def edgeset(E):
    return set((min(u, v), max(u, v)) for u, v in E)


# ---------------------------------------------------------------------------
# check bookkeeping
# ---------------------------------------------------------------------------

_PASSES = 0
_FAILS = 0


def check(name, ok, detail=''):
    global _PASSES, _FAILS
    if ok:
        _PASSES += 1
        print('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        _FAILS += 1
        print('FAIL %s%s' % (name, (' ' + detail) if detail else ''))


# ===========================================================================
# THE OBJECTS OF THE PAPER, section 2
# ===========================================================================

# (tree builder, closure edge (u,v), cycle vertices as printed, name)
WITNESSES = [
    ('H_1',       lambda: TG(2, 6),                  (0, 3),  [0, 2, 3]),
    ('H_star',    lambda: patt_m(P(1), S2(4), 2, 3), (0, 3),  [0, 1, 2, 3]),
    ('H_2star',   lambda: patt_m(P(3), S2(4), 2, 8), (0, 5),  [0, 3, 4, 5]),
    ('H_A',       lambda: patt_m(T1(3), S2(4), 2, 16), (0, 10), [0, 8, 9, 10]),
    ('H_B',       lambda: patt_m(T1(7), S2(5), 2, 13), (1, 72), [1, 0, 62, 63, 72]),
    ('H_7',       lambda: TG(7, 6),                  (0, 3),  [0, 2, 3]),
]

# Theorem 1's table, verbatim from the paper.
TABLE = {
    'H_1':     dict(n=82,  alpha=43,  mode=25, B=[40, 42]),
    'H_star':  dict(n=58,  alpha=31,  mode=18, B=[28, 29]),
    'H_2star': dict(n=155, alpha=82,  mode=47, B=[78, 79, 80]),
    'H_A':     dict(n=312, alpha=164, mode=95, B=[162, 163]),
    'H_B':     dict(n=315, alpha=164, mode=97, B=[161, 162, 163]),
    'H_7':     dict(n=282, alpha=148, mode=88,
                    B=[135, 137, 139, 141, 143, 145, 147]),
}

# The coefficients the paper prints, keyed by index.
PRINTED_COEFFS = {
    'H_1': {39: 78555205917, 40: 125984042, 41: 622042, 42: 485, 43: 1},
    'H_2star': {77: 1135175576719, 78: 3247342473, 79: 9297742, 80: 51245,
                81: 315, 82: 1},
    'H_A': {160: 14399685489601, 161: 15059609304, 162: 12154773, 163: 9815,
            164: 9},
    'H_B': {159: 1836675801123505133, 160: 946417052360718, 161: 367642616716,
            162: 146769284, 163: 137465, 164: 129},
    'H_7': {134: 176034211309656718484117918890005878151,
            135: 966322839077538885713606551341161815,
            136: 5610322511779111055524711146550627,
            137: 19499704262596590695297758206411,
            138: 69973942869517833994033432381,
            139: 174438847810507838032800383,
            140: 444508352733368700464454,
            141: 822213872739649589233,
            142: 1559903682459718229,
            143: 2122169611386407,
            144: 3022643280657,
            145: 2853892767,
            146: 2959617,
            147: 1535,
            148: 1},
}

# The full coefficient list the paper prints for H*, i_0..i_31.
HSTAR_COEFFS = [
    1, 58, 1595, 27709, 341818, 3190686, 23448526, 139359013, 682663978,
    2794877442, 9661712907, 28414927226, 71476627753, 154326951878,
    286560418561, 457792570649, 628535097664, 739676512678, 742815300664,
    632430324122, 452340586723, 268386668712, 129808551984, 49924871616,
    14717851602, 3137062048, 435476016, 31145353, 415590, 6304, 113, 1,
]

# The graph6 string the paper prints for H*, the five printed lines concatenated.
HSTAR_G6 = (
    'ylC_H?@G??g??@??_?_?@?@???G?G???E??????G???_??@???C????_'
    '??A????@???C?????G??@??????C????@?????_?????_????O?????@'
    '?????_??????K????????????C??????G??????G??????O??????@??'
    '????A????????_?????@????????@??????C?????????G???????@??'
    '??????O????????G???????A?????????C???????@??????????G'
)

# The products the paper prints, as (name, k, i_k^2, i_{k-1} i_{k+1}, is_break).
PRINTED_PRODUCTS = [
    ('H_1',      40, 15871978838657764, 48864637399022514, True),
    ('H_1',      41, 386936249764, 61102260370, False),
    ('H_1',      42, 235225, 622042, True),
    ('H_star',   27, 970033013494609, 180979477489440, False),
    ('H_star',   28, 172715048100, 196340305312, True),
    ('H_star',   29, 39740416, 46961670, True),
    ('H_star',   30, 12769, 6304, False),
    ('H_2star',  78, 10545233136949755729, 10554569637034468498, True),
    ('H_2star',  79, 86448006298564, 166410065028885, True),
    ('H_2star',  80, 2626050025, 2928788730, True),
    ('H_2star',  81, 99225, 51245, False),
    ('H_A',     161, 226791832389123364416, 175024908397494015573, False),
    ('H_A',     162, 147738506681529, 147810065318760, True),
    ('H_A',     163, 96334225, 109392957, True),
    ('H_B',     161, 135161093625787682624656, 138904953140373090585912, True),
    ('H_B',     162, 21541222725872656, 50537992306864940, True),
    ('H_B',     163, 18896626225, 18933237636, True),
]


def main():
    print('=== verify.py -- An Affirmative Answer to the Levit-Kadrawi '
          'Multi-Break Unicyclic Closure Problem ===')
    print('exact integer arithmetic throughout; no floating point decides anything')
    print()

    # -------------------------------------------------------------------
    # 0.  THE MACHINERY, against brute force and against closed forms
    # -------------------------------------------------------------------
    print('--- 0. the machinery ---')

    # 0a. every forest and every unicyclic graph of order <= 6, against brute force.
    pairs = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    nf = nu = 0
    bad = None
    for n in range(2, 7):
        idx = [(i, j) for i, j in pairs if j < n]
        for mask in range(1 << len(idx)):
            E = [idx[b] for b in range(len(idx)) if mask >> b & 1]
            m = len(E)
            if m == n - 1 and is_forest(E, n) and connected(E, n):
                got = pstrip(forest_poly(adjacency(E, n), list(range(n))))
                if got != brute_poly(E, n):
                    bad = ('tree', n, E)
                nf += 1
            elif m == n and connected(E, n):
                cyc = unique_cycle(E, n)
                u, v = next((a, b) for a, b in E if a in cyc and b in cyc)
                a = poly_route_A(E, n, u, v)
                b = poly_route_B(E, n, min(cyc))
                if a != brute_poly(E, n) or b != brute_poly(E, n):
                    bad = ('unicyclic', n, E)
                nu += 1
    check('machinery-vs-brute-force', bad is None,
          '%d trees and %d unicyclic graphs of order <= 6; forest recursion and both '
          'identities (A) and (B) agree with 2^n subset enumeration on every one'
          % (nf, nu))

    # 0b. paths: i_k(P_n) = C(n-k+1, k); no breaks (P_n is claw-free, hence real-rooted).
    okp = True
    for n in range(3, 41):
        E, nn, _ = P(n)
        c = pstrip(forest_poly(adjacency(E, nn), list(range(nn))))
        # alpha(P_n) = ceil(n/2)
        want = [comb(n - k + 1, k) for k in range(0, (n + 1) // 2 + 1)]
        if c != want or breaks(c) or not is_unimodal(c):
            okp = False
    check('control-paths', okp,
          'i_k(P_n) = C(n-k+1,k) reproduced coefficient by coefficient for n = 3..40; '
          '0 breaks on all 38')

    # 0c. cycles, computed as unicyclic graphs by both routes:
    #     i_k(C_n) = (n/(n-k)) C(n-k, k); no breaks.
    okc = True
    for n in range(3, 41):
        E, nn, _ = cycle(n)
        a = poly_route_A(E, nn, 0, n - 1)
        b = poly_route_B(E, nn, 0)
        want = [1] + [(n * comb(n - k, k)) // (n - k) for k in range(1, n // 2 + 1)]
        for k in range(1, n // 2 + 1):
            if (n * comb(n - k, k)) % (n - k):
                okc = False
        if a != b or a != want or breaks(a) or not is_unimodal(a):
            okc = False
    check('control-cycles', okc,
          'i_k(C_n) = (n/(n-k))C(n-k,k) reproduced for n = 3..40 by BOTH identities; '
          '0 breaks on all 38.  With the paths above: the break detector is silent on '
          'all 76 log-concave controls')

    # 0d. the break detector on hand-checkable synthetic sequences.
    synth = [([1, 4, 3, 1], []), ([1, 1, 10, 1, 1], [1, 3]),
             ([9, 1, 9, 1, 9, 1, 9], [1, 3, 5]), ([1, 2, 4, 8, 16], []),
             ([1, 5, 3, 7, 4], [2])]
    check('control-detector', all(breaks(c) == b for c, b in synth),
          '5 synthetic sequences with 0, 1, 2, 3 and 1 breaks classified exactly')

    # 0e. construction sizes, i.e. equation (2) of the paper.
    oks = True
    for t in range(1, 9):
        if S2(t)[1] != 1 + 2 * t:
            oks = False
    for l in range(1, 9):
        if T1(l)[1] != 2 + 2 * l:
            oks = False
    for m in range(1, 8):
        for t in range(4, 9):
            if TG(m, t)[1] != m * (4 + 6 * t) + 2:
                oks = False
    for m, n in ((3, 3), (4, 4), (5, 8)):
        if KL(m, n)[1] != 10 + 2 * m + 2 * n:
            oks = False
    sized = [(patt_m(P(1), S2(4), 2, 3), 1 + 3 * (1 + 2 * 9)),
             (patt_m(P(3), S2(4), 2, 8), 3 + 8 * 19),
             (patt_m(T1(3), S2(4), 2, 16), 8 + 16 * 19),
             (patt_m(T1(7), S2(5), 2, 13), 16 + 13 * 23)]
    for G, want in sized:
        if G[1] != want:
            oks = False
    check('construction-sizes', oks,
          '|S_{2,t}| = 1+2t, |T_{1,l}| = 2+2l, |T_{G,m,t}| = m(4+6t)+2 and equation '
          '(2) |(G:H)_k^(m)| = |G| + m(1+k|H|) hold on every parameter the paper uses')

    # -------------------------------------------------------------------
    # 1.  PUBLISHED INTEGERS OF OTHER AUTHORS (positive controls)
    # -------------------------------------------------------------------
    print()
    print('--- 1. positive controls on OTHER authors\' published integers ---')

    # Bautista-Ramos, Fig. 1 caption: GT_{2,5} has degree 37 and breaks at 34 and 36.
    E, n, _ = TG(2, 5)
    c = pstrip(forest_poly(adjacency(E, n), list(range(n))))
    check('control-BR-GT25', len(c) - 1 == 37 and breaks(c) == [34, 36],
          'GT_{2,5}: degree %d, breaks %s (published: degree 37, indices 34 and 36)'
          % (len(c) - 1, breaks(c)))

    okm = True
    got = []
    for m in (4, 5):
        E, n, _ = TG(m, 6)
        c = pstrip(forest_poly(adjacency(E, n), list(range(n))))
        got.append((m, breaks(c)))
    if got != [(4, [78, 80, 82, 84]), (5, [97, 99, 101, 103, 105])]:
        okm = False
    check('control-BR-multiplicity', okm,
          'GT_{4,6} breaks %s and GT_{5,6} breaks %s -- m breaks for m = 4, 5, as '
          'Bautista-Ramos reports for t large enough' % (got[0][1], got[1][1]))

    # BRGG cor:infiniteTwo -- (T_{1,3}:S_{2,4})_2^(m), m >= 16: alpha = 10m+4,
    # two consecutive breaks at 10m+2, 10m+3.
    ok2 = True
    det2 = []
    for m in (16, 17, 18):
        E, n, _ = patt_m(T1(3), S2(4), 2, m)
        c = pstrip(forest_poly(adjacency(E, n), list(range(n))))
        if len(c) - 1 != 10 * m + 4 or breaks(c) != [10 * m + 2, 10 * m + 3]:
            ok2 = False
        det2.append('m=%d:alpha=%d,B=%s' % (m, len(c) - 1, breaks(c)))
    check('control-BRGG-two-consecutive', ok2,
          'alpha = 10m+4 with breaks {10m+2, 10m+3} at %s' % ' '.join(det2))

    # BRGG cor:infiniteThree -- (T_{1,7}:S_{2,5})_2^(m), m >= 13: alpha = 12m+8,
    # three consecutive breaks at 12m+5, 12m+6, 12m+7.
    ok3 = True
    det3 = []
    for m in (13, 14):
        E, n, _ = patt_m(T1(7), S2(5), 2, m)
        c = pstrip(forest_poly(adjacency(E, n), list(range(n))))
        if len(c) - 1 != 12 * m + 8 or breaks(c) != [12 * m + 5, 12 * m + 6, 12 * m + 7]:
            ok3 = False
        det3.append('m=%d:alpha=%d,B=%s' % (m, len(c) - 1, breaks(c)))
    check('control-BRGG-three-consecutive', ok3,
          'alpha = 12m+8 with breaks {12m+5, 12m+6, 12m+7} at %s' % ' '.join(det3))

    # Kadrawi-Levit: trees of order 26 already fail log-concavity.  Over the family
    # T_{3,m,n}, 1 <= m <= n <= 8, the least order carrying a break is 26 (T_{3,4,4}).
    least = None
    who = None
    for m in range(1, 9):
        for nn in range(m, 9):
            E, n, _ = KL(m, nn)
            c = pstrip(forest_poly(adjacency(E, n), list(range(n))))
            if breaks(c) and (least is None or n < least):
                least = n
                who = (m, nn, breaks(c))
    check('control-KL-order-26', least == 26 and who[:2] == (4, 4),
          'least order in T_{3,m,n}, 1 <= m,n <= 8, carrying a break is %d, at '
          'T_{3,%d,%d} with breaks %s -- and nothing smaller breaks'
          % (least, who[0], who[1], who[2]))

    # -------------------------------------------------------------------
    # 2.  THE SIX WITNESSES OF THEOREM 1
    # -------------------------------------------------------------------
    print()
    print('--- 2. the six witnesses of Theorem 1 ---')
    polys = {}
    trees = {}
    for name, build, (u, v), cyc_printed in WITNESSES:
        T = build()
        H = add_edge(T, u, v)
        E, n, _ = H
        want = TABLE[name]

        # structure
        ok = (n == want['n'] and len(E) == n and connected(E, n)
              and is_forest(list(T[0]), n) and connected(list(T[0]), n))
        cyc = unique_cycle(E, n)
        ok = ok and cyc == set(cyc_printed)
        check('%s-structure' % name, ok,
              'n = |E| = %d, connected, tree + one edge, unique cycle = %s '
              '(length %d), as printed'
              % (n, sorted(cyc), len(cyc)))

        # the polynomial, twice
        pa = poly_route_A(E, n, u, v)
        pb = poly_route_B(E, n, sorted(cyc)[0])
        check('%s-two-routes' % name, pa == pb,
              'identity (A) (delete cycle edge %s) and identity (B) (split at cycle '
              'vertex %d) give the same %d coefficients'
              % ((u, v), sorted(cyc)[0], len(pa)))
        c = pa
        polys[name] = c

        mode, uniq = unique_mode(c)
        check('%s-alpha-unimodal' % name,
              len(c) - 1 == want['alpha'] and is_unimodal(c) and uniq
              and mode == want['mode'],
              'alpha = %d, unimodal, unique mode at i_%d' % (len(c) - 1, mode))

        bs = breaks(c)
        check('%s-breaks-exactly' % name, bs == want['B'],
              'B(%s) = %s -- every one of the %d indices 1 <= k <= alpha-1 was tested, '
              'so this is exact, not a subset' % (name, bs, len(c) - 2))

        # closure preserves or amplifies
        TE, tn, _ = T
        ct = pstrip(forest_poly(adjacency(TE, tn), list(range(tn))))
        trees[name] = ct
        bt = breaks(ct)
        rel = 'PRESERVES' if bt == bs else ('AMPLIFIES %d -> %d' % (len(bt), len(bs)))
        check('%s-tree-comparison' % name, is_unimodal(ct) and len(ct) == len(c),
              'the tree has alpha = %d and breaks %s; the closure %s'
              % (len(ct) - 1, bt, rel))

    # the paper's printed coefficients
    for name, want in PRINTED_COEFFS.items():
        c = polys[name]
        bad = [k for k, v in want.items() if c[k] != v]
        check('%s-printed-coefficients' % name, not bad,
              '%d coefficients printed in the paper (i_%d..i_%d) all agree'
              % (len(want), min(want), max(want)))

    check('H_star-full-coefficient-list', polys['H_star'] == HSTAR_COEFFS,
          'all 32 coefficients i_0..i_31 printed in the paper agree with the '
          'recomputation')

    # the paper's printed products
    badp = []
    for name, k, sq, cross, isbreak in PRINTED_PRODUCTS:
        c = polys[name]
        if c[k] * c[k] != sq or c[k - 1] * c[k + 1] != cross:
            badp.append((name, k, 'value'))
        elif (sq < cross) != isbreak:
            badp.append((name, k, 'direction'))
    check('printed-products', not badp,
          '%d printed products i_k^2 and i_{k-1}i_{k+1} recomputed exactly, and each '
          'inequality points the way the paper says' % len(PRINTED_PRODUCTS))

    # H_7's seven breaks and six holds, spelled out
    c7 = polys['H_7']
    ok7 = (all(c7[k] ** 2 < c7[k - 1] * c7[k + 1] for k in range(135, 148, 2))
           and all(c7[k] ** 2 >= c7[k - 1] * c7[k + 1] for k in range(136, 147, 2)))
    check('H_7-alternating', ok7,
          'the 7 inequalities at k = 135,137,...,147 all fail and the 6 at '
          'k = 136,138,...,146 all hold')

    # the graph6 string printed for H*
    nn, Eg = g6_decode(HSTAR_G6)
    Hs = add_edge(patt_m(P(1), S2(4), 2, 3), 0, 3)
    check('H_star-graph6', len(HSTAR_G6) == 277 and nn == 58
          and Eg == edgeset(Hs[0]),
          'the 277-character graph6 string printed in the paper decodes to n = %d and '
          'to exactly the %d edges of the constructed H*, label for label'
          % (nn, len(Eg)))

    # the one small structural claim made in prose: i_2 drops by one
    Ts = patt_m(P(1), S2(4), 2, 3)
    cts = trees['H_star']
    check('H_star-i2-drop', cts[2] == 1596 and polys['H_star'][2] == 1595
          and cts[2] == comb(58, 2) - 57 and polys['H_star'][2] == comb(58, 2) - 58,
          'i_2 goes 1596 -> 1595 from T* to H*, so the polynomial really does change '
          'even though the break set does not (Ts order %d)' % Ts[1])

    # -------------------------------------------------------------------
    # 3.  THE TWO EXHAUSTIVE CENSUSES OF PROPOSITION 4
    # -------------------------------------------------------------------
    print()
    print('--- 3. the two exhaustive censuses of Proposition 4 ---')
    for label, build, want_ne, want_dist, want_exact, want_run in (
            ('T_{G,2,6}', lambda: TG(2, 6), 3240, {1: 864, 2: 2376}, (2, [40, 42]), 1),
            ('T*', lambda: patt_m(P(1), S2(4), 2, 3), 1596,
             {0: 417, 1: 207, 2: 972}, (2, [28, 29]), 2)):
        T = build()
        E, n, _ = T
        NE = nonedges(E, n)
        check('%s-nonedge-count' % label,
              len(NE) == want_ne == comb(n, 2) - (n - 1),
              'enumerated %d nonedges; C(%d,2) - (%d-1) = %d -- the enumeration and the '
              'arithmetic agree' % (len(NE), n, n, comb(n, 2) - (n - 1)))
        dist = {}
        runs = {}
        nonuni = 0
        exact = 0
        for (u, v) in NE:
            c = poly_route_A(list(E) + [(u, v)], n, u, v)
            if not is_unimodal(c):
                nonuni += 1
            bs = breaks(c)
            dist[len(bs)] = dist.get(len(bs), 0) + 1
            runs[longest_run(bs)] = runs.get(longest_run(bs), 0) + 1
            if bs == want_exact[1]:
                exact += 1
        check('%s-census' % label,
              nonuni == 0 and dist == want_dist
              and exact == want_dist[want_exact[0]]
              and max(runs) == want_run,
              'complete census of all %d closures: 0 non-unimodal, break-count '
              'distribution %s, all %d of the %d-break closures have break set exactly '
              '%s, longest consecutive run over the cell = %d'
              % (len(NE), sorted(dist.items()), exact, want_exact[0],
                 want_exact[1], max(runs)))

    # -------------------------------------------------------------------
    # 4.  THE REMARK'S TABLE: preserve versus amplify along T_{G,m,6}
    # -------------------------------------------------------------------
    print()
    print('--- 4. the Remark: preserving versus amplifying ---')
    ROW = {1: ([21], [21]),
           2: ([40, 42], [40, 42]),
           3: ([59, 61, 63], [59, 61, 63]),
           4: ([78, 80, 82, 84], [78, 80, 82, 84]),
           5: ([97, 99, 101, 103, 105], [97, 99, 101, 103, 105]),
           6: ([116, 118, 120, 122, 124], [116, 118, 120, 122, 124, 126]),
           7: ([135, 137, 139, 141, 143, 145],
               [135, 137, 139, 141, 143, 145, 147])}
    badr = []
    for m in range(1, 8):
        E, n, _ = TG(m, 6)
        ct = pstrip(forest_poly(adjacency(E, n), list(range(n))))
        c = poly_route_A(list(E) + [(0, 3)], n, 0, 3)
        if (breaks(ct), breaks(c)) != ROW[m] or n != m * 40 + 2 \
                or len(c) - 1 != len(ct) - 1 or not is_unimodal(c) \
                or not is_unimodal(ct):
            badr.append(m)
    check('amplification-table', not badr,
          'all 7 rows m = 1..7 of the table reproduce: orders %s, every tree and every '
          'closure unimodal, the closure preserving the break set for m <= 5 and '
          'adding one break at index alpha-1 for m = 6, 7'
          % [m * 40 + 2 for m in range(1, 8)])

    # the non-monotonicity the paper records in section 6
    E, n, _ = TG(5, 5)
    c55 = pstrip(forest_poly(adjacency(E, n), list(range(n))))
    okmono = not breaks(c55)
    for m in (6, 7):
        E, n, _ = TG(m, 5)
        if breaks(pstrip(forest_poly(adjacency(E, n), list(range(n))))):
            okmono = False
    check('non-monotonicity', okmono,
          'T_{G,m,5} has NO break for m = 5, 6, 7, so the break count is not a '
          'function of m alone -- the caveat section 6 records')

    # -------------------------------------------------------------------
    # scope
    # -------------------------------------------------------------------
    print()
    print('NOT RE-RUN: this program checks the paper, and the paper only.  It does NOT')
    print('NOT RE-RUN: enumerate the nonedges of the four cells other than T_{G,2,6}')
    print('NOT RE-RUN: and T* -- in particular no census of (T_{1,3}:S_{2,4})_2^(16) or')
    print('NOT RE-RUN: of (T_{1,7}:S_{2,5})_2^(13) is performed, and none is claimed;')
    print('NOT RE-RUN: H_A and H_B are checked as the two individual graphs they are.')
    print('NOT RE-RUN: It does NOT search for smaller witnesses, so no minimality claim')
    print('NOT RE-RUN: is supported.  It does NOT verify any statement about arbitrarily')
    print('NOT RE-RUN: many breaks: the largest number of breaks it exhibits is 7, and')
    print('NOT RE-RUN: the m = 1..7 table is a finite computation, not a proof of a law.')
    print('NOT RE-RUN: It does NOT re-read the source e-print arXiv:2603.17114, so the')
    print('NOT RE-RUN: quotations, line numbers and byte counts in section 1 of the paper')
    print('NOT RE-RUN: are outside its reach and must be checked against the e-print.')
    print()

    n = _PASSES + _FAILS
    if _FAILS:
        print('VERDICT: %d of %d CHECKS FAILED' % (_FAILS, n))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _PASSES)
    return 0


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""verify.py -- checks every computational claim of

    "Corner-free 4-cop-win graphs on 20 vertices, and a retraction question of
     Turcotte and Yvon"   (paper.tex / paper.pdf in this folder)

Python 3.9 or later, STANDARD LIBRARY ONLY: no third-party package, no external data file, no
network, no randomness, integer arithmetic throughout.  Every object it consumes is printed in
the paper: the two graph6 strings of Section 2, the two adjacency tables of Section 2, the two
3-element cop sets and the dominating set of Section 3, and the two control graphs (Petersen,
Robertson) of Section 5.  They are transcribed here as literals and nothing else is read.

WHAT IT DOES

  1. Decodes each graph6 string with its own bit-unpacker, re-encodes it, and checks that the
     decoded edge set equals the adjacency table printed beside it in the paper.  So the two
     encodings in the paper are cross-checked against each other, not merely restated.
  2. Re-derives order, size, regularity, connectivity, triangle count, maximum number of common
     neighbours, and the girth (by breadth-first search from every vertex).
  3. Counts CORNERS -- ordered pairs (x, u), x != u, with N(x) subset of N[u] -- over all
     n(n-1) pairs.  Zero corners is what kills every order-(n-1) retract.
  4. Re-checks the forest certificate of Lemma 3.2 for the printed 3-set, and then EXHAUSTIVELY
     over all C(20,3) = 1140 three-element subsets.
  5. Re-checks the perfect-code certificate for W1 and computes the domination number of each
     graph exhaustively.
  6. Decides the cop number of each graph OUTRIGHT, with an independent bitmask least-fixed-point
     solver for the k-cop game (cops move first, either player may pass).  This shares no code
     with the numpy/scipy solver used when the result was found; it is a second implementation.
  7. Counts the automorphism group of each graph by backtracking, which is how the paper's claim
     that W1 and W2 ARE Meringer's published census graphs (|Aut| = 96 and 20) is checked.
  8. Runs both polarities of every detector on published integers that are not ours: the cop
     solver on P_5 (1), C_4 (2), C_5 (2), Petersen (3) and Robertson (4); the corner detector on
     C_4, where it must find corners; the forest lemma on Robertson, where it must succeed, and
     on Petersen with one cop, where it must refuse a bound the literature contradicts.

Output: one `PASS <name>` line per check, the closing verdict, and the NOT RE-RUN list.  Exit 0
if and only if every check passed.
"""
import sys
from itertools import combinations

# ---------------------------------------------------------------------------
# 0. THE OBJECTS, TRANSCRIBED FROM THE PAPER
# ---------------------------------------------------------------------------
# Section 2, the two graph6 strings.
G6 = {
    'W1': 'S????B?K?W?WIg`cIP?j?QSA`GCoGBGG?',
    'W2': 'S??CA?_c?o@_Q_iOAcU@OK`@cO?MO?Bg?',
}

# Section 2, the two adjacency tables, copied line for line from the paper.
ADJ_TABLE = {
    'W1': """
     0:8,13,17,18   1:8,14,16,19   2:9,12,17,19   3:9,14,15,18   4:10,12,16,18
     5:10,13,15,19  6:11,12,13,14  7:11,15,16,17  8:0,1,12,15    9:2,3,13,16
    10:4,5,14,17   11:6,7,18,19   12:2,4,6,8     13:0,5,6,9     14:1,3,6,10
    15:3,5,7,8     16:1,4,7,9     17:0,2,7,10    18:0,3,4,11    19:1,2,5,11
    """,
    'W2': """
     0:6,9,13,15    1:7,12,15,17   2:8,13,16,17   3:9,10,14,16   4:10,11,12,13
     5:11,14,17,18  6:0,12,16,18   7:1,13,18,19   8:2,14,15,19   9:0,3,17,19
    10:3,4,15,18   11:4,5,16,19   12:1,4,6,14    13:0,2,4,7     14:3,5,8,12
    15:0,1,8,10    16:2,3,6,11    17:1,2,5,9     18:5,6,7,10    19:7,8,9,11
    """,
}

# Section 3, the three-element cop sets whose complements induce forests.
COP_SET = {'W1': (0, 14, 16), 'W2': (0, 2, 5)}
# Section 3, the perfect dominating set of W1.
PERFECT_CODE = (8, 9, 10, 11)
# Section 2, the 5-cycles exhibited to pin the girth from above.
FIVE_CYCLE = {'W1': (0, 8, 15, 3, 18), 'W2': (0, 6, 12, 1, 15)}

# Section 5, the two published control graphs.
PETERSEN_G6 = 'IheA@GUAo'
ROBERTSON_G6 = 'R???C@_E?KI_I`bGCi?f?ocAoGAQO?'

# Published integers the controls are held to.  None of these is ours.
PUBLISHED_COP_NUMBER = {'P5': 1, 'C4': 2, 'C5': 2, 'Petersen': 3, 'Robertson': 4}
PUBLISHED_AUT_ORDER = {'Petersen': 120, 'Robertson': 24, 'W1': 96, 'W2': 20}

# ---------------------------------------------------------------------------
# 1. THE CHECK LEDGER
# ---------------------------------------------------------------------------
_passes = []
_fails = []


def check(name, ok, detail=''):
    if ok:
        _passes.append(name)
        print('PASS %-34s %s' % (name, detail))
    else:
        _fails.append(name)
        print('FAIL %-34s %s' % (name, detail))
    return bool(ok)


# ---------------------------------------------------------------------------
# 2. GRAPHS
# ---------------------------------------------------------------------------
def decode_graph6(s):
    """graph6 -> (n, frozenset of (i, j) with i < j).  Own bit-unpacker."""
    b = [ord(c) - 63 for c in s]
    n, bits = b[0], []
    for x in b[1:]:
        for i in range(5, -1, -1):
            bits.append((x >> i) & 1)
    edges, k = set(), 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]:
                edges.add((i, j))
            k += 1
    return n, frozenset(edges)


def encode_graph6(n, edges):
    """The inverse map, written independently of the decoder, for a round trip."""
    assert 1 <= n <= 62, 'this encoder covers only 1 <= n <= 62'
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if (i, j) in edges else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for p in range(0, len(bits), 6):
        v = 0
        for b in bits[p:p + 6]:
            v = v * 2 + b
        out.append(chr(v + 63))
    return ''.join(out)


def parse_adj_table(text):
    """`v:a,b,c,d` tokens -> (n, frozenset of (i, j) with i < j)."""
    edges, rows = set(), set()
    for tok in text.split():
        head, tail = tok.split(':')
        v = int(head)
        rows.add(v)
        for w in tail.split(','):
            w = int(w)
            edges.add((min(v, w), max(v, w)))
    return len(rows), frozenset(edges)


def cycle_graph(n):
    return n, frozenset((min(i, (i + 1) % n), max(i, (i + 1) % n)) for i in range(n))


def path_graph(n):
    return n, frozenset((i, i + 1) for i in range(n - 1))


def adjacency(n, edges):
    """open neighbourhoods as sets, closed neighbourhoods as bitmasks"""
    nbr = [set() for _ in range(n)]
    for a, b in edges:
        nbr[a].add(b)
        nbr[b].add(a)
    closed = [(1 << v) | sum(1 << w for w in nbr[v]) for v in range(n)]
    return nbr, closed


def connected(n, nbr):
    seen, stack = {0}, [0]
    while stack:
        v = stack.pop()
        for w in nbr[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def girth(n, nbr):
    """Shortest cycle length by BFS from every vertex; None if acyclic."""
    best = None
    for root in range(n):
        dist = {root: 0}
        parent = {root: None}
        queue = [root]
        head = 0
        while head < len(queue):
            v = queue[head]
            head += 1
            for w in nbr[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1
                    parent[w] = v
                    queue.append(w)
                elif parent[v] != w:
                    cand = dist[v] + dist[w] + 1
                    if best is None or cand < best:
                        best = cand
    return best


def triangles(edges, nbr):
    return sum(1 for a, b in edges if nbr[a] & nbr[b])


def max_common_neighbours(n, nbr):
    return max(len(nbr[a] & nbr[b]) for a in range(n) for b in range(a + 1, n))


def corners(n, nbr):
    """Ordered pairs (x, u), x != u, with N(x) subset of N[u]: the retraction obstruction."""
    return [(x, u) for x in range(n) for u in range(n)
            if x != u and nbr[x] <= (nbr[u] | {u})]


def is_forest(vertices, edges):
    """union-find over the induced subgraph"""
    par = {v: v for v in vertices}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        par[ra] = rb
    return True


def forest_remainder(n, edges, nbr, S):
    """(R, induced edges of R) for R = V \\ N[S]"""
    closed_S = set(S)
    for v in S:
        closed_S |= nbr[v]
    R = sorted(set(range(n)) - closed_S)
    Rset = set(R)
    RE = sorted((a, b) for a, b in edges if a in Rset and b in Rset)
    return R, RE


def forest_certificate(n, edges, nbr, S):
    R, RE = forest_remainder(n, edges, nbr, S)
    return R, RE, is_forest(R, RE)


def count_forest_sets(n, edges, nbr, size):
    """how many size-element subsets S have V \\ N[S] inducing a forest"""
    return sum(1 for S in combinations(range(n), size)
               if forest_certificate(n, edges, nbr, S)[2])


def domination_number(n, closed, upto=None):
    """(gamma, number of minimum dominating sets) by exhaustive search over increasing sizes"""
    full = (1 << n) - 1
    limit = n if upto is None else upto
    for size in range(1, limit + 1):
        hits = 0
        for S in combinations(range(n), size):
            m = 0
            for v in S:
                m |= closed[v]
            if m == full:
                hits += 1
        if hits:
            return size, hits
    return None, 0


# ---------------------------------------------------------------------------
# 3. THE COP NUMBER, DECIDED FROM SCRATCH
# ---------------------------------------------------------------------------
def cops_win(n, edges, k):
    """True iff k cops catch the robber on (n, edges), under the convention of the target paper:
    the cops choose their k vertices (repeats allowed), then the robber chooses his, then the
    cops move first, and on a move each cop and the robber may stay put or step to a neighbour.

    Least fixed point.  A cop configuration is an ORDERED k-tuple, indexed in base n; for each
    configuration we hold the BITMASK of robber vertices from which the cops win with the cops to
    move.  One iteration is
        A[c]  = bits(c)  |  { r : N[r] subset of R[c] }        (catch now, or every robber reply
                                                                already lost)
        R'[c] = R[c]  |  OR over c' in prod_i N[c_i] of A[c']  (the cops pick their move)
    and the OR over the product is done one coordinate at a time, which is what makes this cheap.
    """
    nbr, closed = adjacency(n, edges)
    step_nbrs = [sorted(nbr[v] | {v}) for v in range(n)]
    full = (1 << n) - 1
    size = n ** k
    place = [n ** (k - 1 - i) for i in range(k)]
    occupied = [0] * size
    for idx in range(size):
        m, rest = 0, idx
        for p in place:
            m |= 1 << (rest // p)
            rest %= p
        occupied[idx] = m
    R = list(occupied)
    iters = 0
    while True:
        iters += 1
        A = [0] * size
        for idx in range(size):
            reach = R[idx]
            a = occupied[idx]
            if reach:
                for r in range(n):
                    if closed[r] & ~reach == 0:
                        a |= 1 << r
            A[idx] = a
        B = A
        for i in range(k):
            stride = place[i]
            block = stride * n
            C = [0] * size
            for base in range(0, size, block):
                for off in range(stride):
                    cells = [B[base + off + j * stride] for j in range(n)]
                    for v in range(n):
                        acc = 0
                        for w in step_nbrs[v]:
                            acc |= cells[w]
                        C[base + off + v * stride] = acc
            B = C
        newR = [R[i] | B[i] for i in range(size)]
        if newR == R:
            break
        R = newR
    return any(m == full for m in R), iters


def cop_number(n, edges, cap):
    """the least k <= cap with k cops winning, or None"""
    for k in range(1, cap + 1):
        win, _ = cops_win(n, edges, k)
        if win:
            return k
    return None


# ---------------------------------------------------------------------------
# 4. AUTOMORPHISMS
# ---------------------------------------------------------------------------
def automorphism_count(n, edges):
    """|Aut(G)| by backtracking along a breadth-first vertex order."""
    A = [[False] * n for _ in range(n)]
    for a, b in edges:
        A[a][b] = A[b][a] = True
    order, seen = [0], {0}
    head = 0
    while len(order) < n:
        if head < len(order):
            v = order[head]
            head += 1
            for w in range(n):
                if A[v][w] and w not in seen:
                    seen.add(w)
                    order.append(w)
        else:
            for w in range(n):
                if w not in seen:
                    seen.add(w)
                    order.append(w)
                    break
    prefix = [list(order[:i]) for i in range(n)]
    total = 0
    image, used = {}, set()

    def extend(i):
        nonlocal total
        if i == n:
            total += 1
            return
        v = order[i]
        for w in range(n):
            if w in used:
                continue
            ok = True
            for u in prefix[i]:
                if A[v][u] != A[w][image[u]]:
                    ok = False
                    break
            if ok:
                image[v] = w
                used.add(w)
                extend(i + 1)
                used.discard(w)
                del image[v]

    extend(0)
    return total


# ---------------------------------------------------------------------------
# 5. THE CHECKS
# ---------------------------------------------------------------------------
print('=== the two witnesses, as printed in Section 2 of the paper ===')

WIT = {}
for name in ('W1', 'W2'):
    s = G6[name]
    n, edges = decode_graph6(s)
    nbr, closed = adjacency(n, edges)
    WIT[name] = (n, edges, nbr, closed)

    check('g6-length-%s' % name, len(s) == 33,
          'graph6 line is %d bytes = 1 + ceil(C(20,2)/6) = 1 + 32' % len(s))
    check('g6-roundtrip-%s' % name, encode_graph6(n, edges) == s,
          're-encoding the decoded graph reproduces the printed string exactly')
    tn, tedges = parse_adj_table(ADJ_TABLE[name])
    check('encodings-agree-%s' % name, (tn, tedges) == (n, edges),
          'graph6 and the printed adjacency table describe the SAME labelled graph '
          '(%d rows, %d edges)' % (tn, len(tedges)))
    check('order-%s' % name, n == 20,
          'n = %d = M_4 + 1, since M_4 = 19 (target paper, its own corollary)' % n)
    check('size-%s' % name, len(edges) == 40, 'm = %d' % len(edges))
    degs = sorted({len(nbr[v]) for v in range(n)})
    check('regular-%s' % name, degs == [4],
          'every degree is 4, so delta = Delta = 4 (degree multiset {%s})'
          % ','.join(str(d) for d in degs))
    check('connected-%s' % name, connected(n, nbr), 'one component')
    tri = triangles(edges, nbr)
    check('triangle-free-%s' % name, tri == 0, '%d triangles over the %d edges' % (tri, len(edges)))
    mx = max_common_neighbours(n, nbr)
    check('c4-free-%s' % name, mx <= 1,
          'max |N(a) cap N(b)| = %d over all C(20,2) = %d pairs, so no C_4'
          % (mx, n * (n - 1) // 2))
    g = girth(n, nbr)
    cyc = FIVE_CYCLE[name]
    cyc_ok = all(cyc[(i + 1) % 5] in nbr[cyc[i]] for i in range(5)) and len(set(cyc)) == 5
    check('girth-%s' % name, g == 5 and cyc_ok,
          'BFS girth = %d, and the exhibited cycle %s is adjacency-checked'
          % (g, '-'.join(str(v) for v in cyc) + '-%d' % cyc[0]))

    cs = corners(n, nbr)
    check('corner-free-%s' % name, not cs,
          '%d corners over all %d ordered pairs (x,u): no induced subgraph on 19 vertices is a '
          'retract, and none on fewer is reached by a chain' % (len(cs), n * (n - 1)))

    S = COP_SET[name]
    R, RE, ok = forest_certificate(n, edges, nbr, S)
    check('forest-certificate-%s' % name, ok,
          'S = %s, |N[S]| = %d, R = %s induces exactly %s -- %d edges on %d vertices and acyclic, '
          'so c <= 4 by Lemma 3.2'
          % (list(S), n - len(R), R, ' '.join('%d-%d' % e for e in RE), len(RE), len(R)))
    nsets = count_forest_sets(n, edges, nbr, 3)
    check('forest-census-%s' % name, nsets > 0,
          '%d of the C(20,3) = %d three-element subsets have a forest remainder (exhaustive)'
          % (nsets, 20 * 19 * 18 // 6))

    lower_ok = (g >= 5 and degs == [4])
    check('aigner-fromme-hypotheses-%s' % name, lower_ok,
          'girth %d >= 5 and delta = 4, which is exactly what Theorem 3.1 (Aigner-Fromme) needs '
          'for c >= delta = 4' % g)

print()
print('=== the domination route, and the domination numbers ===')

n1, e1, nbr1, cl1 = WIT['W1']
cover = {}
for d in PERFECT_CODE:
    for v in nbr1[d] | {d}:
        cover[v] = cover.get(v, 0) + 1
check('perfect-code-W1', len(cover) == 20 and set(cover.values()) == {1},
      'D = %s: the four closed neighbourhoods are disjoint and cover all %d vertices once each, '
      'so gamma(W1) <= 4 and four cops finish on move one'
      % ('{' + ','.join(str(d) for d in PERFECT_CODE) + '}', len(cover)))
g1, m1 = domination_number(n1, cl1, upto=4)
check('gamma-W1', g1 == 4 and m1 == 1,
      'gamma(W1) = %d exhaustively -- no subset of size 3 or less dominates, and D is the UNIQUE '
      'minimum dominating set (%d of the C(20,4) = 4845 four-subsets dominate)' % (g1, m1))
n2, e2, nbr2, cl2 = WIT['W2']
g2, m2 = domination_number(n2, cl2, upto=5)
check('gamma-W2', g2 == 5,
      'gamma(W2) = %d exhaustively -- 0 of the C(20,4) = 4845 four-subsets dominate and %d of the '
      'C(20,5) = 15504 five-subsets do, so the domination route does NOT reach W2 and W2 needs '
      'the forest certificate' % (g2, m2))
check('witnesses-not-isomorphic', g1 != g2,
      'gamma(W1) = %d differs from gamma(W2) = %d, so W1 and W2 are non-isomorphic: the '
      'refutation rests on two distinct graphs' % (g1, g2))

print()
print('=== the cop number, decided outright by a second independent solver ===')

for name in ('W1', 'W2'):
    n, edges, nbr, closed = WIT[name]
    lose3, it3 = cops_win(n, edges, 3)
    check('three-cops-lose-%s' % name, not lose3,
          'the 3-cop game is a robber win (fixed point in %d iterations), so c(%s) >= 4 without '
          'appeal to Theorem 3.1' % (it3, name))
    win4, it4 = cops_win(n, edges, 4)
    check('four-cops-win-%s' % name, win4,
          'the 4-cop game is a cop win (fixed point in %d iterations), so c(%s) <= 4 without '
          'appeal to Lemma 3.2' % (it4, name))
    check('cop-number-%s' % name, (not lose3) and win4,
          'c(%s) = 4 exactly, on 20 = M_4 + 1 vertices' % name)

print()
print('=== provenance: these are Meringer census graphs, not constructions of ours ===')

for name in ('W1', 'W2'):
    n, edges, nbr, closed = WIT[name]
    a = automorphism_count(n, edges)
    check('automorphisms-%s' % name, a == PUBLISHED_AUT_ORDER[name],
          '|Aut(%s)| = %d, matching the published census entry (%s of 20_4_5.asc)'
          % (name, a, 'Graph 2' if name == 'W1' else 'Graph 1'))

print()
print('=== controls, both polarities, on published integers that are not ours ===')

CONTROLS = {
    'P5': path_graph(5),
    'C4': cycle_graph(4),
    'C5': cycle_graph(5),
    'Petersen': decode_graph6(PETERSEN_G6),
    'Robertson': decode_graph6(ROBERTSON_G6),
}

nP, eP = CONTROLS['Petersen']
nbrP, clP = adjacency(nP, eP)
check('control-petersen-shape', (nP, len(eP)) == (10, 15) and
      sorted({len(nbrP[v]) for v in range(nP)}) == [3] and girth(nP, nbrP) == 5,
      'the transcribed graph6 really is the Petersen graph: 10 vertices, 15 edges, cubic, girth 5')
nR, eR = CONTROLS['Robertson']
nbrR, clR = adjacency(nR, eR)
check('control-robertson-shape', (nR, len(eR)) == (19, 38) and
      sorted({len(nbrR[v]) for v in range(nR)}) == [4] and girth(nR, nbrR) == 5,
      'the transcribed graph6 really is the (4,5)-cage: 19 vertices, 38 edges, 4-regular, girth 5')

for name in ('P5', 'C4', 'C5', 'Petersen', 'Robertson'):
    n, edges = CONTROLS[name]
    want = PUBLISHED_COP_NUMBER[name]
    got = cop_number(n, edges, want + 1)
    check('control-cop-number-%s' % name, got == want,
          'solver returns c = %s, published value %d%s'
          % (got, want, '  [the crucial negative: it must not say 4]'
             if name == 'Petersen' else
             '  [the positive control AT the value in question, on the paper\'s own M_4 witness]'
             if name == 'Robertson' else
             '  [the authors\' own k = 2 refuter]' if name == 'C5' else
             '  [M_2 = 4, the minimum 2-cop-win graph]' if name == 'C4' else ''))

nC4, eC4 = CONTROLS['C4']
nbrC4, _ = adjacency(nC4, eC4)
csC4 = corners(nC4, nbrC4)
check('corner-detector-positive-control', len(csC4) == 4,
      'the corner detector finds %d corners in C_4 (girth 4), so its zero on W1 and W2 is a '
      'measurement and not a constant' % len(csC4))
nC5, eC5 = CONTROLS['C5']
nbrC5, _ = adjacency(nC5, eC5)
check('corner-free-C5', not corners(nC5, nbrC5),
      'C_5 is corner-free too, which is why the authors\' own k = 2 counterexample works and why '
      'the mechanism used here is theirs')

fR = count_forest_sets(nR, eR, nbrR, 3)
check('forest-lemma-forced-positive', fR == 555,
      '%d of the C(19,3) = %d three-element subsets of the Robertson graph have a forest '
      'remainder, so Lemma 3.2 delivers c(Robertson) <= 4, and with Theorem 3.1 the published 4'
      % (fR, 19 * 18 * 17 // 6))
fP = count_forest_sets(nP, eP, nbrP, 1)
check('forest-lemma-anti-control', fP == 0,
      '%d of the 10 one-element subsets of the Petersen graph have a forest remainder, so the '
      'same lemma REFUSES the false bound c(Petersen) <= 2; the published value is 3' % fP)

print()
print('NOT RE-RUN: the completeness of the census lane.  This program generates no graphs, so it '
      'does not re-derive that the connected 4-regular graphs of girth at least 5 on 20 vertices '
      'are exactly two.  That count is Meringer\'s published census (file 20_4_5.asc) and OEIS '
      'A058343; what is checked here is that the two graphs named in the paper have the claimed '
      'automorphism orders, 96 and 20.')
print('NOT RE-RUN: M_4 = 19.  It is the target paper\'s own corollary and is used, not verified; '
      'the paper\'s claim is conditional on it in exactly the way the source states it.')
print('NOT RE-RUN: girth 3 and girth 4 graphs on 20 vertices.  Nothing here bounds how many '
      '20-vertex graphs with cop number 4 and no corner exist, so no minimality or uniqueness is '
      'claimed.')
print('NOT RE-RUN: k = 2, k = 3 and k >= 5 of the question.  k = 2 is refuted in the source '
      'itself, k = 3 is TRUE there by its own classification, and k >= 5 is untouched because '
      'M_5 is unknown.  Only k = 4 is settled here.')
print('NOT RE-RUN: the authors\' published table of 25148 cop numbers.  The solver above is '
      'controlled on five published values, not regression-tested against that table.')

print()
total = len(_passes) + len(_fails)
if _fails:
    print('VERDICT: %d of %d CHECKS FAILED -- %s' % (len(_fails), total, ', '.join(_fails)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(_passes))
sys.exit(0)

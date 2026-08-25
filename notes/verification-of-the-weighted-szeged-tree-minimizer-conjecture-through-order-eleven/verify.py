#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Referee verification program for
"Verification of the Weighted Szeged Tree-Minimizer Conjecture Through Order Eleven"
(settles affirmatively, for 3 <= n <= 11, Conjecture 2 of Bok, Furtula,
Jedlickova and Skrekovski).

Standard library only.  Exact integer / Fraction arithmetic throughout; no
floating point is used in any decision.

------------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs; never used to prove themselves)
------------------------------------------------------------------------------
  * all eleven graph6 witness strings of the witness table and their claimed
    wSz values:
        n=3..9 trees   Bo 12, Cq 34, DkG 72, Eh_O 130, FCOf? 204,
                       GiE@?O 306, HiQ@?_G 432
        n=10 tree      I?AA@_gw?     578
        n=10 non-tree  I?AA@_gwO     600
        n=11 tree      J?AA@AOEB_?   762
        n=11 non-tree  J?AACGoIF??   778
  * the paper's structural description of the tree witnesses: the spiders
    S(2,2,2), S(2,2,3), S(2,2,2,2) at n = 7,8,9, and two adjacent degree-3
    vertices carrying legs of length 2, with one leg of length 3 at n = 11,
    at n = 10,11;
  * Table 1: connected-class counts, tree minima, non-tree minima, gaps,
    for n = 3..11, and the grand total 1,018,690,327 classes;
  * the claim "the order-10 census was run as 37 disjoint edge-count shards".
Every one of these is an INPUT that the program recomputes and compares.

------------------------------------------------------------------------------
WHAT IS DERIVED HERE (the checks)
------------------------------------------------------------------------------
  * graph6 decoding, well-formedness and byte-exact re-encoding of all eleven
    witnesses; their order, edge count, connectivity, tree/non-tree status;
  * wSz of each witness from the definition, by TWO independent distance
    kernels (bitmask BFS and Floyd-Warshall), plus a third edge-split kernel
    for trees;
  * that each tree witness attains the exhaustive tree minimum of its order
    and is one of the generated isomorphism classes;
  * that each tree witness at n = 7..11 is isomorphic (AHU canonical form) to
    the shape the paper names for it, rebuilt from that description alone;
  * that the per-edge-count minima of the full n = 7 and n = 8 census are not
    monotone in m -- the fact the SCOPE report cites when it calls the
    truncated n = 10, 11 searches a scope limit rather than a proof;
  * wSz(P_3) = wSz(K_3) = 12 and the n=3 tie;
  * the complete set of free trees on n vertices, 3 <= n <= 11 (AHU canonical
    form), hence the EXACT tree minimum for every n in 3..11;
  * the complete set of isomorphism classes of connected graphs for
    3 <= n <= 8 (refinement-pruned canonical form), hence exact class counts
    and exact tree / non-tree minima, and the fact that for 4 <= n <= 8 every
    minimizer is a tree;
  * a covering enumeration of every isomorphism class of connected graph on
    9 vertices, hence the exact n=9 minima;
  * for n = 10 and n = 11, an exhaustive spanning-tree-plus-k-extra-edges
    search over every connected graph with at most 13 edges;
  * the connected-class counts for n = 3..11 and their total, from first
    principles: Burnside/cycle-index count of all graphs, then the inverse
    Euler transform (multiset-of-connected-components) -- independent of OEIS;
  * the free-tree counts, cross-checked against Otter's formula;
  * the "37 edge-count shards at order 10" arithmetic, validated against the
    exhaustive edge-count sets of the generated n = 3..8 census (the count
    45 - 9 + 1 = 37 on its own is arithmetic, not evidence);
  * the PRECONDITIONS of the two covering arguments, which are what makes them
    exhaustive: the 8-vertex class list feeding the order-9 covering
    enumeration is verified connected, pairwise non-isomorphic and of the full
    cardinality 11117, and the augmentation is verified to reach all 47
    nine-vertex free-tree classes; the spanning-tree seed list of the
    order-10/11 sparse search is verified duplicate-free and of Otter
    cardinality, and every per-edge-count minimum in that search is confirmed
    by the second distance kernel.  (An enumeration size that matches its own
    binomial formula by construction certifies nothing and is not on its own
    the content of any check here.)

NOT re-derived: the >= (n+2)-edge part of the census for n = 10, 11, except for
the dense tail, where wSz(G) >= sum_v deg(v)^2 closes m >= 39 at n = 10 and
m >= 47 at n = 11 by proof.  See the SCOPE lines printed at the end.
"""

import sys
import time
import itertools
from fractions import Fraction
from math import factorial, gcd

_T0 = time.time()
_CHECKS = []


def check(name, ok, detail=""):
    """Record one check.  Prints 'PASS <name>' or 'FAIL <name>'."""
    ok = bool(ok)
    _CHECKS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    sys.stdout.flush()


def note(msg):
    print("-- " + str(msg))
    sys.stdout.flush()


# ---------------------------------------------------------------- paper inputs
P_ORDERS = list(range(3, 12))
P_CLASSES = {3: 2, 4: 6, 5: 21, 6: 112, 7: 853, 8: 11117,
             9: 261080, 10: 11716571, 11: 1006700565}
P_TREEMIN = {3: 12, 4: 34, 5: 72, 6: 130, 7: 204, 8: 306,
             9: 432, 10: 578, 11: 762}
P_NONTREEMIN = {3: 12, 4: 36, 5: 76, 6: 132, 7: 212, 8: 316,
                9: 444, 10: 600, 11: 778}
P_GAP = {3: 0, 4: 2, 5: 4, 6: 2, 7: 8, 8: 10, 9: 12, 10: 22, 11: 16}
P_TOTAL = 1018690327
P_WITNESS = [(10, "tree", "I?AA@_gw?", 578),
             (10, "nontree", "I?AA@_gwO", 600),
             (11, "tree", "J?AA@AOEB_?", 762),
             (11, "nontree", "J?AACGoIF??", 778)]
# The remaining seven rows of the witness table, n = 3..9 (all trees).
P_WITNESS_SMALL = [(3, "Bo", 12),
                   (4, "Cq", 34),
                   (5, "DkG", 72),
                   (6, "Eh_O", 130),
                   (7, "FCOf?", 204),
                   (8, "GiE@?O", 306),
                   (9, "HiQ@?_G", 432)]
# The paper's structural description of the tree witnesses: spiders (one centre
# with legs of the stated lengths) at n = 7,8,9, and two adjacent degree-3
# centres carrying legs at n = 10,11 (one leg of length 3 at n = 11).
P_SPIDER = {7: [2, 2, 2], 8: [2, 2, 3], 9: [2, 2, 2, 2]}
P_TWO_CENTRE = {10: ([2, 2], [2, 2]), 11: ([2, 2], [2, 3])}
P_SHARDS_N10 = 37
P_TREE_TIE_N3 = 12


# ----------------------------------------------------------------- graph6 I/O
def g6_decode(s):
    """Decode a graph6 string.  Returns (n, sorted edge list).

    Raises ValueError on any malformed input: bad character range, wrong
    length for the declared order, or a non-zero padding bit.
    """
    if not s:
        raise ValueError("empty string")
    for ch in s:
        if not (63 <= ord(ch) <= 126):
            raise ValueError("byte out of graph6 range: %r" % ch)
    vals = [ord(c) - 63 for c in s]
    n = vals[0]
    if n >= 63:
        raise ValueError("this program handles orders < 63 only")
    need = n * (n - 1) // 2
    nbytes = (need + 5) // 6
    if len(vals) - 1 != nbytes:
        raise ValueError("expected %d data bytes for n=%d, found %d"
                         % (nbytes, n, len(vals) - 1))
    bits = []
    for v in vals[1:]:
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    if any(bits[need:]):
        raise ValueError("non-zero padding bit")
    edges = []
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                edges.append((i, j))
            idx += 1
    return n, edges


def g6_encode(n, edges):
    """Re-encode (n, edges) as graph6; inverse of g6_decode."""
    es = set()
    for (u, v) in edges:
        es.add((min(u, v), max(u, v)))
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if (i, j) in es else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for i in range(0, len(bits), 6):
        v = 0
        for k in range(6):
            v = (v << 1) | bits[i + k]
        out.append(chr(v + 63))
    return "".join(out)


# --------------------------------------------------- graph basics and kernel 1
def adjacency(n, edges):
    adj = [0] * n
    for (u, v) in edges:
        if u == v:
            raise ValueError("loop")
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def is_connected(n, edges):
    adj = adjacency(n, edges)
    seen = 1
    front = 1
    while front:
        nxt = 0
        f = front
        while f:
            b = f & -f
            nxt |= adj[b.bit_length() - 1]
            f ^= b
        nxt &= ~seen
        seen |= nxt
        front = nxt
    return seen == (1 << n) - 1


def bfs_dist(n, adj):
    """All-pairs distances by repeated bitmask BFS.  None if disconnected."""
    full = (1 << n) - 1
    dist = []
    for s in range(n):
        d = [0] * n
        seen = 1 << s
        front = 1 << s
        lev = 0
        while front:
            lev += 1
            nxt = 0
            f = front
            while f:
                b = f & -f
                nxt |= adj[b.bit_length() - 1]
                f ^= b
            nxt &= ~seen
            if not nxt:
                break
            m = nxt
            while m:
                b = m & -m
                d[b.bit_length() - 1] = lev
                m ^= b
            seen |= nxt
            front = nxt
        if seen != full:
            return None
        dist.append(d)
    return dist


def wsz_bfs(n, edges):
    """KERNEL 1.  wSz(G) = sum_{uv in E} (deg u + deg v) n_u(uv) n_v(uv),
    distances by repeated BFS.  Returns None if G is disconnected."""
    adj = adjacency(n, edges)
    dist = bfs_dist(n, adj)
    if dist is None:
        return None
    deg = [bin(a).count("1") for a in adj]
    tot = 0
    for (u, v) in edges:
        du = dist[u]
        dv = dist[v]
        nu = 0
        nv = 0
        for a, b in zip(du, dv):
            if a < b:
                nu += 1
            elif b < a:
                nv += 1
        tot += (deg[u] + deg[v]) * nu * nv
    return tot


def wsz_fw(n, edges):
    """KERNEL 2.  Same invariant, distances by Floyd-Warshall on an explicit
    matrix with a large finite sentinel.  Returns None if disconnected."""
    INF = 10 ** 6
    D = [[0 if i == j else INF for j in range(n)] for i in range(n)]
    deg = [0] * n
    for (u, v) in edges:
        D[u][v] = 1
        D[v][u] = 1
        deg[u] += 1
        deg[v] += 1
    for k in range(n):
        Dk = D[k]
        for i in range(n):
            Di = D[i]
            dik = Di[k]
            if dik >= INF:
                continue
            for j in range(n):
                alt = dik + Dk[j]
                if alt < Di[j]:
                    Di[j] = alt
    for i in range(n):
        for j in range(n):
            if D[i][j] >= INF:
                return None
    tot = 0
    for (u, v) in edges:
        nu = 0
        nv = 0
        for x in range(n):
            a = D[x][u]
            b = D[x][v]
            if a < b:
                nu += 1
            elif b < a:
                nv += 1
        tot += (deg[u] + deg[v]) * nu * nv
    return tot


def split_counts(n, edges):
    """For every edge, return (n_u, n_v, equidistant) using BFS distances."""
    adj = adjacency(n, edges)
    dist = bfs_dist(n, adj)
    if dist is None:
        return None
    out = []
    for (u, v) in edges:
        du = dist[u]
        dv = dist[v]
        nu = nv = eq = 0
        for a, b in zip(du, dv):
            if a < b:
                nu += 1
            elif b < a:
                nv += 1
            else:
                eq += 1
        out.append((nu, nv, eq))
    return out


def wsz_treesplit(n, edges):
    """KERNEL 3, trees only.  Uses NO distance matrix: deleting edge uv from a
    tree leaves components of sizes a and n-a, and n_u = a, n_v = n-a."""
    if len(edges) != n - 1 or not is_connected(n, edges):
        raise ValueError("not a tree")
    adj = [[] for _ in range(n)]
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    deg = [len(a) for a in adj]
    tot = 0
    for (u, v) in edges:
        stack = [u]
        seen = set([u, v])
        cnt = 1
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    cnt += 1
                    stack.append(y)
        tot += (deg[u] + deg[v]) * cnt * (n - cnt)
    return tot


def is_tree(n, edges):
    return len(edges) == n - 1 and is_connected(n, edges)


# ------------------------------------------- free trees: AHU canonical form
def ahu_canon(n, edges):
    """Canonical string of a free tree: AHU encoding rooted at its centre.
    Two trees are isomorphic iff their canonical strings are equal.
    Raises ValueError on input that is not a tree."""
    if n == 1:
        return "()"
    if not is_tree(n, edges):
        raise ValueError("ahu_canon needs a tree (n=%d, m=%d)"
                         % (n, len(edges)))
    adj = [[] for _ in range(n)]
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    deg = [len(a) for a in adj]
    alive = [True] * n
    rem = n
    layer = [v for v in range(n) if deg[v] <= 1]
    while rem > 2:
        if not layer:
            raise ValueError("leaf stripping stalled: input is not a tree")
        nxt = []
        for v in layer:
            if not alive[v]:
                continue
            alive[v] = False
            rem -= 1
            for w in adj[v]:
                if alive[w]:
                    deg[w] -= 1
                    if deg[w] == 1:
                        nxt.append(w)
        layer = [v for v in nxt if alive[v]]
    centres = [v for v in range(n) if alive[v]]

    def enc(root, banned):
        # iterative post-order encoding, so recursion depth is never an issue
        order = []
        parent = {root: banned}
        stack = [root]
        while stack:
            x = stack.pop()
            order.append(x)
            for y in adj[x]:
                if y != parent[x]:
                    parent[y] = x
                    stack.append(y)
        code = {}
        for x in reversed(order):
            kids = sorted(code[y] for y in adj[x] if y != parent[x])
            code[x] = "(" + "".join(kids) + ")"
        return code[root]

    if len(centres) == 1:
        return enc(centres[0], -1)
    a, b = centres
    return "[" + "".join(sorted([enc(a, b), enc(b, a)])) + "]"


def all_free_trees(nmax):
    """trees[n] = one edge list per isomorphism class of tree on n vertices.
    Built by attaching a new leaf to every vertex of every (n-1)-vertex tree,
    which reaches every tree, then deduplicated by AHU canonical form."""
    trees = {1: [[]]}
    for n in range(2, nmax + 1):
        seen = {}
        for rep in trees[n - 1]:
            for v in range(n - 1):
                e = list(rep) + [(v, n - 1)]
                k = ahu_canon(n, e)
                if k not in seen:
                    seen[k] = e
        trees[n] = list(seen.values())
    return trees


def spider_tree(legs):
    """(n, edge list) of the spider S(legs): one centre joined to disjoint paths
    of the stated lengths.  Built from the paper's description only."""
    edges = []
    nxt = 1
    for L in legs:
        prev = 0
        for _ in range(L):
            edges.append((min(prev, nxt), max(prev, nxt)))
            prev = nxt
            nxt += 1
    return nxt, edges


def two_centre_tree(legs_a, legs_b):
    """(n, edge list) of two adjacent centres carrying the stated legs (paths).
    With two legs on each centre both centres have degree 3, which is the shape
    the paper ascribes to the n = 10 and n = 11 tree witnesses."""
    edges = [(0, 1)]
    nxt = 2
    for (centre, legs) in ((0, legs_a), (1, legs_b)):
        for L in legs:
            prev = centre
            for _ in range(L):
                edges.append((min(prev, nxt), max(prev, nxt)))
                prev = nxt
                nxt += 1
    return nxt, edges


# ------------------------- general graphs: canonical form and class generation
def pair_index(n):
    """T[i][j] = bit position of the pair {i,j} (graph6 column-major order)."""
    T = [[0] * n for _ in range(n)]
    for j in range(1, n):
        for i in range(j):
            T[i][j] = T[j][i] = j * (j - 1) // 2 + i
    return T


def refine_colours(n, adj):
    """1-dimensional Weisfeiler-Leman refinement: an isomorphism-invariant
    vertex colouring, used only to prune the permutation search."""
    col = [bin(adj[v]).count("1") for v in range(n)]
    for _ in range(n):
        keys = []
        for v in range(n):
            nb = []
            m = adj[v]
            while m:
                b = m & -m
                nb.append(col[b.bit_length() - 1])
                m ^= b
            nb.sort()
            keys.append((col[v], tuple(nb)))
        order = sorted(set(keys))
        rank = dict((k, i) for i, k in enumerate(order))
        new = [rank[k] for k in keys]
        if new == col:
            break
        col = new
    return col


def graph_canon(n, edges, T):
    """Canonical key of a graph: the minimum, over all relabellings that
    respect the refined colouring, of the packed upper-triangle bitmask.
    Exact (the colouring is isomorphism-invariant, so no class is split)."""
    adj = adjacency(n, edges)
    col = refine_colours(n, adj)
    groups = {}
    for v in range(n):
        groups.setdefault(col[v], []).append(v)
    blocks = []
    start = 0
    for k in sorted(groups):
        g = groups[k]
        blocks.append((g, list(range(start, start + len(g)))))
        start += len(g)
    perms = [list(itertools.permutations(pos)) for (g, pos) in blocks]
    best = None
    p = [0] * n
    for combo in itertools.product(*perms):
        for (g, pos), asn in zip(blocks, combo):
            for v, q in zip(g, asn):
                p[v] = q
        m = 0
        for (u, v) in edges:
            m |= 1 << T[p[u]][p[v]]
        if best is None or m < best:
            best = m
    return best


def edges_from_mask(n, mask, T):
    out = []
    for j in range(1, n):
        for i in range(j):
            if (mask >> T[i][j]) & 1:
                out.append((i, j))
    return out


def connected_classes(nmax):
    """classes[n] = one edge list per isomorphism class of CONNECTED graph on
    n vertices, for n <= nmax.  Every connected graph has a non-cut vertex, so
    deleting it leaves a connected graph on n-1 vertices; hence attaching a new
    vertex by every non-empty neighbourhood to every (n-1)-class reaches every
    n-class.  Deduplicated by graph_canon."""
    classes = {1: [[]]}
    for n in range(2, nmax + 1):
        T = pair_index(n)
        seen = {}
        for rep in classes[n - 1]:
            for S in range(1, 1 << (n - 1)):
                e = list(rep)
                m = S
                while m:
                    b = m & -m
                    e.append((b.bit_length() - 1, n - 1))
                    m ^= b
                k = graph_canon(n, e, T)
                if k not in seen:
                    seen[k] = edges_from_mask(n, k, T)
        classes[n] = list(seen.values())
    return classes


# ------------------------------------------- independent counts (no OEIS used)
def binom(a, b):
    if b < 0 or b > a:
        return 0
    r = 1
    for i in range(b):
        r = r * (a - i) // (i + 1)
    return r


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield []
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest


def count_all_graphs(n):
    """Number of isomorphism classes of simple graphs on n vertices, by
    Burnside applied to the induced action of S_n on the C(n,2) pairs.  A
    permutation with cycle lengths L_1..L_k has sum floor(L_i/2) +
    sum_{i<j} gcd(L_i,L_j) orbits on pairs."""
    if n == 0:
        return 1
    total = 0
    nfac = factorial(n)
    for lam in partitions(n):
        mult = {}
        for p in lam:
            mult[p] = mult.get(p, 0) + 1
        cls = nfac
        for j, mj in mult.items():
            cls //= (j ** mj) * factorial(mj)
        orb = sum(p // 2 for p in lam)
        for i in range(len(lam)):
            for j in range(i + 1, len(lam)):
                orb += gcd(lam[i], lam[j])
        total += cls * (2 ** orb)
    if total % nfac:
        raise ArithmeticError("Burnside sum not divisible by n!")
    return total // nfac


def mobius(k):
    r = 1
    x = k
    p = 2
    while p * p <= x:
        if x % p == 0:
            x //= p
            if x % p == 0:
                return 0
            r = -r
        p += 1
    if x > 1:
        r = -r
    return r


def inverse_euler(b, N):
    """Given b[0..N] with b[0]=1 counting multisets of the sought objects
    (an arbitrary graph is a multiset of its connected components), recover the
    component counts a[1..N] from  1+sum b_m x^m = prod (1-x^n)^{-a_n}.
    Exact rational arithmetic; a non-integral answer is an error."""
    B = [Fraction(x) for x in b]
    L = [Fraction(0)] * (N + 1)
    for m in range(1, N + 1):
        s = B[m]
        for k in range(1, m):
            s -= Fraction(k, m) * L[k] * B[m - k]
        L[m] = s
    d = [Fraction(0)] * (N + 1)
    for m in range(1, N + 1):
        d[m] = m * L[m]
    a = [0] * (N + 1)
    for m in range(1, N + 1):
        s = Fraction(0)
        for k in range(1, m + 1):
            if m % k == 0:
                s += mobius(m // k) * d[k]
        s = s / m
        if s.denominator != 1:
            raise ArithmeticError("inverse Euler transform not integral")
        a[m] = int(s)
    return a


def rooted_tree_counts(N):
    """r[n] = number of rooted trees on n nodes, from the classical recurrence
    n*r[n+1] = sum_{k=1..n} (sum_{d|k} d*r[d]) * r[n+1-k]."""
    r = [0] * (N + 1)
    if N >= 1:
        r[1] = 1
    for n in range(1, N):
        s = [0] * (n + 1)
        for k in range(1, n + 1):
            tot = 0
            for d in range(1, k + 1):
                if k % d == 0:
                    tot += d * r[d]
            s[k] = tot
        acc = 0
        for k in range(1, n + 1):
            acc += s[k] * r[n + 1 - k]
        if acc % n:
            raise ArithmeticError("rooted tree recurrence not integral")
        r[n + 1] = acc // n
    return r


def otter_free_tree_counts(N):
    """t[n] = number of free trees on n nodes, by Otter's dissimilarity
    formula T(x) = R(x) - (R(x)^2 - R(x^2))/2."""
    r = rooted_tree_counts(N)
    t = [0] * (N + 1)
    for n in range(1, N + 1):
        conv = 0
        for i in range(1, n):
            conv += r[i] * r[n - i]
        if n % 2 == 0:
            conv -= r[n // 2]
        if conv % 2:
            raise ArithmeticError("Otter formula not integral")
        t[n] = r[n] - conv // 2
    return t


# =============================================================================
# CHECKS
# =============================================================================
COMPUTED = {"treemin": {}, "nontreemin": {}, "classes": {},
            "nontree_exhaustive": {}, "edge_minima": {}}


def edge_minimum_descents(orders):
    """orders -> list of (m, min at m, m+1, min at m+1) for every STRICT descent
    of the per-edge-count minimum of the full census of that order, or None when
    that order's census was not generated in this run.  Measured, not asserted:
    the data comes from COMPUTED["edge_minima"], filled in by the census."""
    out = {}
    for n in orders:
        by_m = COMPUTED["edge_minima"].get(n)
        if not by_m:
            out[n] = None
            continue
        out[n] = [(m, by_m[m], m + 1, by_m[m + 1]) for m in sorted(by_m)
                  if m + 1 in by_m and by_m[m + 1] < by_m[m]]
    return out


def describe_descents(descents):
    """Human-readable, and honest when there is nothing to report."""
    parts = []
    for n in sorted(descents):
        d = descents[n]
        if d is None:
            parts.append("n=%d: census not generated in this run" % n)
        elif not d:
            parts.append("n=%d: minima non-decreasing in m, no descent found"
                         % n)
        else:
            parts.append(("n=%d: " % n) + ", ".join(
                "m %d->%d min %d->%d" % (m, m2, w, w2)
                for (m, w, m2, w2) in d))
    return "; ".join(parts)


def check_witnesses():
    """Decode, print back, and evaluate the four exhibited graph6 records."""
    out = {}
    for (n, cls, s, val) in P_WITNESS:
        tag = "n%d_%s" % (n, cls)
        try:
            dn, edges = g6_decode(s)
            back = g6_encode(dn, edges)
            wf = (dn == n and back == s)
            det = "order=%d edges=%d re-encode=%s" % (dn, len(edges), back)
        except ValueError as exc:
            wf = False
            dn, edges = None, None
            det = "decode error: %s" % exc
        check("witness_wellformed_" + tag, wf, det)
        if not wf:
            continue
        conn = is_connected(dn, edges)
        tree = is_tree(dn, edges)
        want_tree = (cls == "tree")
        cyclic = len(edges) >= dn
        ok_cls = conn and (tree == want_tree) and (want_tree or cyclic)
        check("witness_class_" + tag, ok_cls,
              "connected=%s edges=%d n-1=%d tree=%s" %
              (conn, len(edges), dn - 1, tree))
        w1 = wsz_bfs(dn, edges)
        w2 = wsz_fw(dn, edges)
        w3 = wsz_treesplit(dn, edges) if tree else None
        agree = (w1 == val and w2 == val and (w3 is None or w3 == val))
        check("witness_wsz_" + tag, agree,
              "bfs=%s floyd=%s split=%s paper=%d" % (w1, w2, w3, val))
        sc = split_counts(dn, edges)
        if sc is None:
            check("witness_structure_" + tag, False, "graph is disconnected")
            continue
        if want_tree:
            struct = all(a + b == dn and e == 0 for (a, b, e) in sc)
            sdet = "every edge splits V into n_u+n_v=n (tree property)"
        else:
            struct = any(e > 0 for (a, b, e) in sc)
            sdet = "max equidistant on an edge = %d > 0" % max(
                e for (a, b, e) in sc)
        check("witness_structure_" + tag, struct, sdet)
        edge_list = ",".join("%d-%d" % (u, v) for (u, v) in edges)
        note("witness %s %s decodes to n=%d, m=%d, edges %s, wSz=%d"
             % (tag, s, dn, len(edges), edge_list, w1))
        out[(n, cls)] = (dn, edges, w1)
    return out


def tree_witness_strings():
    """order -> graph6 string of the TREE row of the witness table, all eleven
    rows of it (seven small orders plus the two large ones)."""
    out = dict((n, s) for (n, s, _v) in P_WITNESS_SMALL)
    for (n, cls, s, _v) in P_WITNESS:
        if cls == "tree":
            out[n] = s
    return out


def check_small_witnesses(trees):
    """The witness-table rows for n = 3..9, which the large-order witness code
    above does not touch: decode, byte-exact re-encode, evaluate by all three
    kernels, and confirm each attains the exhaustive tree minimum of its order
    and is one of the generated isomorphism classes."""
    for (n, s, val) in P_WITNESS_SMALL:
        tag = "small_witness_n%d" % n
        try:
            dn, edges = g6_decode(s)
            back = g6_encode(dn, edges)
            tree = is_tree(dn, edges)
            det = ("order=%d edges=%d n-1=%d tree=%s re-encode=%s"
                   % (dn, len(edges), dn - 1, tree, back))
        except ValueError as exc:
            dn, edges, tree = None, None, False
            back = None
            det = "decode error: %s" % exc
        ok_form = (dn == n and back == s and tree)
        check(tag + "_wellformed_tree", ok_form, det)
        if not ok_form:
            # An unusable witness FAILS the checks; it never skips them.
            check(tag + "_wsz", False, "witness record unusable")
            check(tag + "_attains_tree_minimum", False,
                  "witness record unusable")
            continue
        w1 = wsz_bfs(dn, edges)
        w2 = wsz_fw(dn, edges)
        w3 = wsz_treesplit(dn, edges)
        check(tag + "_wsz", w1 == val and w2 == val and w3 == val,
              "bfs=%d floyd=%d split=%d paper=%d" % (w1, w2, w3, val))
        target = COMPUTED["treemin"].get(n)
        reps = trees.get(n)
        if target is None or not reps:
            check(tag + "_attains_tree_minimum", False,
                  "no exhaustive tree census for n=%d in this run" % n)
        else:
            member = ahu_canon(dn, edges) in set(
                ahu_canon(n, e) for e in reps)
            check(tag + "_attains_tree_minimum", w1 == target and member,
                  "witness wSz=%d vs exhaustive tree minimum %d (equal: %s); "
                  "the class occurs in the generated set of %d trees: %s"
                  % (w1, target, w1 == target, len(reps), member))
        note("witness n%d_tree %s decodes to n=%d, m=%d, edges %s, wSz=%d"
             % (n, s, dn, len(edges),
                ",".join("%d-%d" % (u, v) for (u, v) in edges), w1))


def check_witness_structures():
    """The paper's structural description of the tree witnesses, checked by
    isomorphism (AHU canonical form) against the shape rebuilt from that
    description alone: spiders at n = 7,8,9 and two adjacent degree-3 centres
    carrying legs at n = 10,11."""
    strings = tree_witness_strings()
    claims = []
    for n in sorted(P_SPIDER):
        legs = P_SPIDER[n]
        claims.append((n, "spider S(%s)" % ",".join(str(L) for L in legs),
                       spider_tree(legs)))
    for n in sorted(P_TWO_CENTRE):
        a, b = P_TWO_CENTRE[n]
        claims.append((n, "two adjacent degree-3 centres with legs %s and %s"
                       % (a, b), two_centre_tree(a, b)))
    for (n, label, (cn, cedges)) in claims:
        name = "witness_is_the_structure_the_paper_names_n%d" % n
        s = strings.get(n)
        if s is None:
            check(name, False, "no tree witness listed for n=%d" % n)
            continue
        try:
            dn, edges = g6_decode(s)
            got = ahu_canon(dn, edges)
        except ValueError as exc:
            check(name, False, "witness unusable: %s" % exc)
            continue
        want = ahu_canon(cn, cedges)
        check(name, dn == n and cn == n and got == want,
              "%s (order %d) is isomorphic to the %s built from the paper's "
              "description (order %d): %s"
              % (s, dn, label, cn, dn == n and cn == n and got == want))


def check_p3_k3():
    """The proof's displayed identity wSz(P_3) = wSz(K_3) = 12."""
    p3 = [(0, 1), (1, 2)]
    k3 = [(0, 1), (1, 2), (0, 2)]
    wp = wsz_bfs(3, p3)
    wk = wsz_bfs(3, k3)
    wp2 = wsz_fw(3, p3)
    wk2 = wsz_fw(3, k3)
    check("P3_K3_both_equal_12",
          wp == wk == wp2 == wk2 == P_TREE_TIE_N3,
          "wSz(P3)=%d wSz(K3)=%d (floyd %d, %d) paper=%d"
          % (wp, wk, wp2, wk2, P_TREE_TIE_N3))
    check("P3_is_tree_K3_is_not", is_tree(3, p3) and not is_tree(3, k3),
          "P3 has 2 edges, K3 has 3")


def check_tree_census(nmax=11):
    """Exhaustive over ALL free trees, 3 <= n <= nmax: the tree-minimum column
    of Table 1 in full, plus exhaustiveness evidence for the generator."""
    t0 = time.time()
    trees = all_free_trees(nmax)
    gen = [len(trees[n]) for n in range(1, nmax + 1)]
    otter = otter_free_tree_counts(nmax)[1:]
    check("free_tree_counts_match_otter_formula", gen == otter,
          "generated %s == Otter %s" % (gen[2:], otter[2:]))
    bad = []
    for n in range(3, nmax + 1):
        mins = []
        for e in trees[n]:
            w1 = wsz_bfs(n, e)
            w3 = wsz_treesplit(n, e)
            if w1 != w3:
                bad.append((n, e, w1, w3))
            mins.append(w1)
        COMPUTED["treemin"][n] = min(mins)
    check("tree_kernels_agree_on_every_tree", not bad,
          "%d trees, distance kernel == edge-split kernel on all"
          % sum(len(trees[n]) for n in range(3, nmax + 1)))
    got = [COMPUTED["treemin"][n] for n in range(3, nmax + 1)]
    want = [P_TREEMIN[n] for n in range(3, nmax + 1)]
    check("tree_minima_column_n3_to_n11", got == want,
          "computed %s vs paper %s" % (got, want))
    note("free trees per order n=3..%d: %s (%.1fs)"
         % (nmax, gen[2:], time.time() - t0))
    return trees


def check_full_census_to_8(nmax=8):
    """Every isomorphism class of connected graph, 3 <= n <= 8, generated and
    evaluated: class counts, both minima, and the minimizer structure."""
    t0 = time.time()
    classes = connected_classes(nmax)
    gen = [len(classes[n]) for n in range(3, nmax + 1)]
    want = [P_CLASSES[n] for n in range(3, nmax + 1)]
    check("connected_class_counts_by_generation_n3_to_n8", gen == want,
          "generated %s vs paper %s" % (gen, want))
    disagree = 0
    n3_min_edges = []
    census_tmin = {}
    for n in range(3, nmax + 1):
        tmin = None
        nmin = None
        gmin = None
        vals = []
        for e in classes[n]:
            w = wsz_bfs(n, e)
            if w != wsz_fw(n, e):
                disagree += 1
            vals.append((w, len(e)))
            if len(e) == n - 1:
                if tmin is None or w < tmin:
                    tmin = w
            else:
                if nmin is None or w < nmin:
                    nmin = w
            if gmin is None or w < gmin:
                gmin = w
        census_tmin[n] = tmin
        COMPUTED["nontreemin"][n] = nmin
        COMPUTED["nontree_exhaustive"][n] = True
        by_m = {}
        for (w, m) in vals:
            if m not in by_m or w < by_m[m]:
                by_m[m] = w
        COMPUTED["edge_minima"][n] = by_m
        mins = [m for (w, m) in vals if w == gmin]
        if n == 3:
            n3_min_edges = sorted(mins)
        else:
            check("every_minimizer_is_a_tree_n%d" % n,
                  gmin == tmin and all(m == n - 1 for m in mins),
                  "global min %d, %d minimizer(s), edge counts %s, n-1=%d"
                  % (gmin, len(mins), sorted(set(mins)), n - 1))
    a = [census_tmin[n] for n in range(3, nmax + 1)]
    b = [COMPUTED["treemin"].get(n) for n in range(3, nmax + 1)]
    check("census_and_tree_enumerator_agree_n3_to_n8", a == b,
          "graph census %s vs independent tree enumerator %s" % (a, b))
    check("kernels_agree_on_all_connected_classes_n3_to_n8", disagree == 0,
          "%d classes evaluated by both kernels, %d disagreements"
          % (sum(gen), disagree))
    check("n3_minimum_attained_by_exactly_P3_and_K3",
          n3_min_edges == [2, 3],
          "minimizing classes at n=3 have edge counts %s (a tree and K_3)"
          % n3_min_edges)
    check("first_order_with_a_connected_non_tree_is_3",
          all(len(e) == n - 1 for n in (1, 2) for e in classes[n])
          and any(len(e) > 2 for e in classes[3]),
          "n=1,2 connected graphs are all trees; n=3 contains K_3")
    got_t = [census_tmin[n] for n in range(3, nmax + 1)]
    got_n = [COMPUTED["nontreemin"][n] for n in range(3, nmax + 1)]
    check("full_census_minima_n3_to_n8",
          got_t == [P_TREEMIN[n] for n in range(3, nmax + 1)]
          and got_n == [P_NONTREEMIN[n] for n in range(3, nmax + 1)],
          "tree %s / non-tree %s" % (got_t, got_n))
    # This is the fact the SCOPE report appeals to when it says the truncated
    # n = 10, 11 searches are a scope limit and not a proof: a minimum over the
    # edge counts BELOW the truncation cannot bound the ones above it, because
    # the per-edge-count minima are not monotone in m even at n = 7, 8.  The
    # check fails if the data is absent, so the SCOPE sentence can never be
    # printed on the strength of a measurement that did not happen.
    descents = edge_minimum_descents((7, 8))
    check("per_edge_count_minima_not_monotone_in_m_n7_and_n8",
          all(descents.get(n) for n in (7, 8)),
          "strict descents of min wSz as m grows: %s" % describe_descents(
              descents))
    note("full census n<=%d done in %.1fs" % (nmax, time.time() - t0))
    return classes


def check_covering_census_9(classes8, trees):
    """Order 9, exhaustive.  Every connected graph on 9 vertices has a non-cut
    vertex, so it arises (up to isomorphism) from some connected 8-vertex class
    by adding a ninth vertex with a non-empty neighbourhood.  Enumerating all
    11117 * 255 such graphs therefore covers every one of the 261080 classes
    (with repetitions, which cannot change a minimum)."""
    t0 = time.time()
    n = 9
    tmin = None
    nmin = None
    cnt = 0
    tree_keys = set()
    for rep in classes8:
        for S in range(1, 1 << 8):
            e = list(rep)
            m = S
            while m:
                b = m & -m
                e.append((b.bit_length() - 1, 8))
                m ^= b
            w = wsz_bfs(n, e)
            cnt += 1
            if len(e) == n - 1:
                tree_keys.add(ahu_canon(n, e))
                if tmin is None or w < tmin:
                    tmin = w
            elif nmin is None or w < nmin:
                nmin = w
    COMPUTED["nontreemin"][n] = nmin
    COMPUTED["nontree_exhaustive"][n] = True
    check("covering_census_n9_minima",
          tmin == P_TREEMIN[9] and nmin == P_NONTREEMIN[9],
          "%d graphs covering all %d classes: tree min %s (paper %d), "
          "non-tree min %s (paper %d)"
          % (cnt, P_CLASSES[9], tmin, P_TREEMIN[9], nmin, P_NONTREEMIN[9]))
    check("every_minimizer_is_a_tree_n9", nmin > tmin,
          "non-tree min %s > tree min %s" % (nmin, tmin))
    # The covering argument is only as good as its INPUT: it certifies order 9
    # exactly when classes8 really is the complete, duplicate-free list of
    # 8-vertex connected classes.  "cnt == len(classes8) * 255" alone is a
    # tautology of the loop and cannot fail, so verify the input as well, and
    # verify the augmentation against the one ground truth available at order 9
    # -- the 47 free trees, enumerated independently by AHU canonical form.
    T8 = pair_index(8)
    keys8 = set()
    shape8 = True
    for rep in classes8:
        if not (7 <= len(rep) <= 28) or not is_connected(8, rep):
            shape8 = False
        keys8.add(graph_canon(8, rep, T8))
    want_tree_keys = set(ahu_canon(n, e) for e in trees[n])
    check("covering_census_n9_size_is_right",
          cnt == len(classes8) * 255
          and shape8
          and len(keys8) == len(classes8) == P_CLASSES[8]
          and tree_keys == want_tree_keys and tree_keys,
          "%d = %d classes x 255 non-empty neighbourhoods; the input classes "
          "are connected 8-vertex graphs, pairwise non-isomorphic (%d distinct "
          "canonical forms) and complete (%d expected); the augmentation "
          "reached all %d nine-vertex tree classes (independent enumerator: "
          "%d)"
          % (cnt, len(classes8), len(keys8), P_CLASSES[8],
             len(tree_keys), len(want_tree_keys)))
    note("order 9 covering census done in %.1fs" % (time.time() - t0))


def check_sparse_census(trees, n, K):
    """Orders 10 and 11, sparse range.  Every connected graph has a spanning
    tree, so up to isomorphism it is one of our tree representatives plus k
    extra edges.  Enumerating k = 0..K therefore covers EVERY connected graph
    on n vertices with at most n-1+K edges."""
    t0 = time.time()
    best = {}
    arg = {}
    cnt = 0
    expect = 0
    for T in trees[n]:
        es = set(T)
        non = [(i, j) for i in range(n) for j in range(i + 1, n)
               if (i, j) not in es]
        for k in range(0, K + 1):
            expect += binom(len(non), k)
        for k in range(0, K + 1):
            for combo in itertools.combinations(non, k):
                e = T + list(combo)
                w = wsz_bfs(n, e)
                m = len(e)
                cnt += 1
                if m not in best or w < best[m]:
                    best[m] = w
                    arg[m] = e
    tmin = best[n - 1]
    nmin = min(v for (m, v) in best.items() if m >= n)
    COMPUTED["nontreemin"][n] = nmin
    COMPUTED["nontree_exhaustive"][n] = False
    check("sparse_census_n%d_tree_min" % n, tmin == P_TREEMIN[n],
          "min over all %d spanning trees = %d (paper %d)"
          % (len(trees[n]), tmin, P_TREEMIN[n]))
    check("sparse_census_n%d_nontree_min_over_m_le_%d"
          % (n, n - 1 + K), nmin == P_NONTREEMIN[n],
          "min over every connected graph with %d..%d edges = %d (paper %d)"
          % (n, n - 1 + K, nmin, P_NONTREEMIN[n]))
    check("minimizers_are_trees_n%d_within_searched_range" % n, nmin > tmin,
          "non-tree min %d > tree min %d" % (nmin, tmin))
    # "cnt == expect" is a tautology of itertools.combinations and cannot fail;
    # on its own it certifies nothing.  Completeness of this search has exactly
    # one precondition -- trees[n] must be the complete, duplicate-free list of
    # free trees on n vertices, since a search seeded with a defective tree list
    # silently returns a WRONG minimum -- and one postcondition, that every
    # edge count in the covered range was actually reached.  Both are verified
    # here, together with a second-kernel confirmation of each per-edge-count
    # minimum (elsewhere in this range only the BFS kernel is used).
    keys = set(ahu_canon(n, e) for e in trees[n])
    otter_n = otter_free_tree_counts(n)[n]
    fw_ok = all(wsz_fw(n, arg[mm]) == best[mm] for mm in sorted(best))
    check("sparse_census_n%d_enumeration_complete" % n,
          cnt == expect
          and len(keys) == len(trees[n]) == otter_n
          and set(best) == set(range(n - 1, n + K))
          and fw_ok,
          "%d graphs evaluated, expected %d; seed tree list: %d entries, %d "
          "distinct isomorphism classes, Otter count %d; edge counts reached "
          "%s (want %s); Floyd-Warshall confirms every per-edge-count "
          "minimum: %s"
          % (cnt, expect, len(trees[n]), len(keys), otter_n, sorted(best),
             list(range(n - 1, n + K)), fw_ok))
    note("order %d, m -> min wSz over ALL connected graphs with m edges: %s"
         % (n, sorted(best.items())))
    note("order %d sparse census: %d graphs (with repeats), %.1fs"
         % (n, cnt, time.time() - t0))


def check_class_counts(classes):
    """The 'connected classes' column and the grand total, from first
    principles: Burnside on the pair action gives all graphs; the inverse Euler
    transform (a graph is a multiset of connected components) gives connected
    graphs.  No OEIS table is consulted."""
    t0 = time.time()
    N = 11
    g = [count_all_graphs(n) for n in range(N + 1)]
    # validate the cycle-index count against brute force for n <= 6
    brute = []
    for n in range(1, 7):
        T = pair_index(n)
        seen = set()
        for mask in range(1 << (n * (n - 1) // 2)):
            seen.add(graph_canon(n, edges_from_mask(n, mask, T), T))
        brute.append(len(seen))
    check("cycle_index_matches_brute_force_orbit_count_n1_to_n6",
          brute == g[1:7],
          "brute force %s == Burnside %s" % (brute, g[1:7]))
    c = inverse_euler(g, N)
    got = [c[n] for n in P_ORDERS]
    want = [P_CLASSES[n] for n in P_ORDERS]
    check("connected_class_counts_n3_to_n11", got == want,
          "computed %s vs paper %s" % (got, want))
    gen = [len(classes[n]) for n in range(3, 9)]
    check("counting_chain_agrees_with_direct_generation_n3_to_n8",
          gen == got[:6], "generated %s == counted %s" % (gen, got[:6]))
    tot = sum(got)
    check("total_census_size", tot == P_TOTAL,
          "sum of the column = %d, paper says %d" % (tot, P_TOTAL))
    for n in P_ORDERS:
        COMPUTED["classes"][n] = c[n]
    note("class counts (%.1fs): %s" % (time.time() - t0, got))


def realizable_edge_counts(n):
    """The set of edge counts m attained by SOME connected graph of order n,
    by explicit construction (spanning star plus extra edges).  Note that the
    star makes every candidate connected, so the is_connected() test below can
    never fire; the construction is therefore validated against exhaustive
    ground truth (the generated census) in check_shard_arithmetic."""
    out = set()
    star = [(0, i) for i in range(1, n)]
    non = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    for k in range(len(non) + 1):
        e = star + non[:k]
        if is_connected(n, e):
            out.add(len(e))
    return out


def check_shard_arithmetic(classes):
    """The paper says order 10 was run as 37 disjoint edge-count shards."""
    n = 10
    realizable = realizable_edge_counts(n)
    # Validate the construction against the exhaustive census: for every order
    # 3..8 the set of edge counts actually realized by a connected class must be
    # exactly {n-1, ..., C(n,2)} and must be exactly what the construction
    # predicts.  Without this the "37" is pure arithmetic on 45 - 9 + 1.
    census_sets_ok = True
    sizes = []
    for k in range(3, 9):
        truth = set(len(e) for e in classes[k])
        want = set(range(k - 1, k * (k - 1) // 2 + 1))
        if truth != want or realizable_edge_counts(k) != truth:
            census_sets_ok = False
        sizes.append(len(truth))
    check("order10_has_exactly_37_edge_count_shards",
          census_sets_ok
          and realizable == set(range(n - 1, n * (n - 1) // 2 + 1))
          and len(realizable) == P_SHARDS_N10,
          "edge counts %d..%d are all realizable, %d values, paper says %d; "
          "the same rule reproduces the exhaustive edge-count sets of the "
          "n=3..8 census (%s values)"
          % (min(realizable), max(realizable), len(realizable),
             P_SHARDS_N10, sizes))
    # Per ORDER, not aggregated: a min/max taken across orders passes as soon as
    # ONE order is tight, so a census that lost every tree at order 8 would slip
    # through.
    bounds_ok = True
    slack = []
    for k in range(3, 9):
        ms = [len(e) for e in classes[k]]
        lo = min(ms) - (k - 1)
        hi = max(ms) - k * (k - 1) // 2
        slack.append((lo, hi))
        if lo != 0 or hi != 0 or any(
                mm < k - 1 or mm > k * (k - 1) // 2 for mm in ms):
            bounds_ok = False
    check("connected_edge_count_bounds_n3_to_n8", bounds_ok,
          "every connected class has n-1 <= m <= C(n,2); per-order slack "
          "(min m -(n-1), max m -C(n,2)) for n=3..8 = %s" % slack)


def check_witness_optimality(W, trees):
    """The witnesses must not merely have the stated values: they must attain
    the computed minima of their classes."""
    for (n, cls, s, val) in P_WITNESS:
        if (n, cls) not in W:
            check("witness_optimal_n%d_%s" % (n, cls), False, "not decoded")
            continue
        dn, edges, w = W[(n, cls)]
        if cls == "tree":
            target = COMPUTED["treemin"].get(n)
            if is_tree(dn, edges):
                keys = set(ahu_canon(n, e) for e in trees[n])
                member = ahu_canon(dn, edges) in keys
            else:
                member = False
            check("witness_optimal_n%d_tree" % n,
                  w == target and member,
                  "witness wSz=%d vs exhaustive tree minimum %d (equal: %s); "
                  "the witness class occurs in the generated set of %d trees: "
                  "%s"
                  % (w, target, w == target, len(trees[n]), member))
        else:
            target = COMPUTED["nontreemin"].get(n)
            # The searched range is m <= n-1+K_SPARSE[n], which is 13 at BOTH
            # n=10 and n=11; the literal "n+3" used here before was right at
            # n=10 and admitted m=14 at n=11, outside the search.
            covered = (COMPUTED["nontree_exhaustive"].get(n, False)
                       or len(edges) <= n - 1 + K_SPARSE.get(n, 0))
            check("witness_optimal_n%d_nontree" % n,
                  w == target and len(edges) >= dn and covered,
                  "witness wSz=%d vs searched non-tree minimum %d (equal: %s); "
                  "m=%d, non-tree: %s, inside the searched range m <= %d: %s"
                  % (w, target, w == target, len(edges), len(edges) >= dn,
                     n - 1 + K_SPARSE.get(n, 0), covered))


def check_table_arithmetic():
    """The gap column, and the theorem's conclusion, from computed minima."""
    gaps = []
    ok_conj = True
    strict = []
    missing = []
    for n in P_ORDERS:
        t = COMPUTED["treemin"].get(n)
        nt = COMPUTED["nontreemin"].get(n)
        if t is None or nt is None:
            missing.append(n)
            gaps.append(None)
            ok_conj = False
            continue
        gaps.append(nt - t)
        if min(t, nt) != t:
            ok_conj = False
        if nt > t:
            strict.append(n)
    if missing:
        note("no computed minimum for orders %s" % missing)
    want = [P_GAP[n] for n in P_ORDERS]
    check("gap_column_n3_to_n11", gaps == want,
          "computed non-tree minus tree = %s, paper = %s" % (gaps, want))
    check("conjecture_2_conclusion_min_over_G_equals_min_over_T", ok_conj,
          "for every n in 3..11 the global minimum over the verified range "
          "equals the tree minimum")
    check("strict_gap_for_n4_to_n11", strict == list(range(4, 12)),
          "orders with a strictly larger non-tree minimum: %s; n=3 ties at %s"
          % (strict, COMPUTED["treemin"].get(3)))
    hdr = "%3s %16s %10s %10s %5s" % ("n", "classes", "tree", "non-tree",
                                      "gap")
    note("recomputed Table 1")
    note(hdr)
    for n, gap in zip(P_ORDERS, gaps):
        note("%3s %16s %10s %10s %5s"
             % (n, COMPUTED["classes"].get(n), COMPUTED["treemin"].get(n),
                COMPUTED["nontreemin"].get(n), gap))


K_SPARSE = {10: 4, 11: 3}


def dense_tail_threshold(n, target):
    """Smallest M with the property that EVERY connected n-vertex graph having
    m >= M edges satisfies wSz(G) > target.

    Proof of the bound used: for every edge uv, u itself satisfies
    d(u,u)=0 < 1 = d(u,v), so n_u >= 1, and symmetrically n_v >= 1; hence
    wSz(G) >= sum_{uv in E} (deg u + deg v) = sum_v deg(v)^2.  Over degree
    sequences with sum 2m the right side is minimised by the balanced sequence,
    giving the exact bound below.  All integer arithmetic."""
    hi = n * (n - 1) // 2
    lb = {}
    for m in range(n - 1, hi + 1):
        q, r = divmod(2 * m, n)
        lb[m] = (n - r) * q * q + r * (q + 1) * (q + 1)
    M = hi + 1
    for m0 in range(hi, n - 2, -1):
        if min(lb[m] for m in range(m0, hi + 1)) > target:
            M = m0
        else:
            break
    return M, lb


def scope_report():
    print("")
    print("SCOPE: what this program re-ran and what it did not.")
    print("SCOPE: exhaustive over EVERY isomorphism class of connected graph "
          "for n = 3,4,5,6,7,8 (generated here) and n = 9 (covering "
          "enumeration); both minima are therefore exact for n <= 9.")
    print("SCOPE: exhaustive over EVERY free tree for n = 3..11, so the whole "
          "tree-minimum column of Table 1 is exact.")
    for n in (10, 11):
        print("SCOPE: for n = %d the non-tree search was exhaustive only over "
              "connected graphs with at most %d of the %d possible edges "
              "(spanning tree plus k <= %d extra edges)."
              % (n, n - 1 + K_SPARSE[n], n * (n - 1) // 2, K_SPARSE[n]))
    print("SCOPE: NOT re-run: the part of the census with at least %d edges at "
          "n = 10 and at least %d edges at n = 11, which is the bulk of the "
          "%d classes and needs nauty geng; the class counts themselves are "
          "instead re-derived exactly by cycle index plus inverse Euler "
          "transform."
          % (10 + K_SPARSE[10], 11 + K_SPARSE[11], P_TOTAL))
    print("SCOPE: the dense tail is partly closed by proof rather than by "
          "search.  Since n_u >= 1 and n_v >= 1 on every edge, "
          "wSz(G) >= sum_v deg(v)^2, which over degree sequences of sum 2m is "
          "minimised by the balanced sequence:")
    for n in (10, 11):
        target = COMPUTED["nontreemin"].get(n)
        if target is None:
            continue
        M, lb = dense_tail_threshold(n, target)
        top = n * (n - 1) // 2
        print("SCOPE:   n = %d: every connected graph with m >= %d edges has "
              "wSz >= %d > %d, so the only edge counts neither searched nor "
              "bounded are %d <= m <= %d (of the full range %d..%d)."
              % (n, M, min(lb[m] for m in range(M, top + 1)), target,
                 n + K_SPARSE[n], M - 1, n - 1, top))
    descents = edge_minimum_descents((7, 8))
    if all(descents.get(n) for n in (7, 8)):
        print("SCOPE: the per-edge-count minima of the full n = 7 and n = 8 "
              "census are NOT monotone in m -- measured in this run, the "
              "strict descents are %s -- so the truncation above is a scope "
              "limit, not a proof; the residual windows just named are the "
              "exact gap between this program and the paper's claim."
              % describe_descents(descents))
    else:
        print("SCOPE: whether the per-edge-count minima of the full n = 7 and "
              "n = 8 census are monotone in m was NOT established in this run "
              "(%s), and the corresponding check above is FAILED; the "
              "truncation above is in any case a scope limit and not a proof, "
              "and the residual windows just named are the exact gap between "
              "this program and the paper's claim."
              % describe_descents(descents))
    print("SCOPE: all eleven rows of the paper's witness table are re-decoded, "
          "re-encoded byte for byte and re-evaluated here, and the structural "
          "descriptions of the tree witnesses at n = 7..11 (spiders, and two "
          "adjacent degree-3 centres with legs) are checked by isomorphism.")
    print("SCOPE: the paper's shard replications (mod 64, 256, 240, 512) are "
          "not reproduced; two distance kernels and, for trees, a third "
          "edge-split kernel are.  The per-edge-count shard SIZES at n = 10 "
          "(the OEIS A054924 comparison) are not recomputed: only the total "
          "class count and the set of realizable edge counts are.")


def verdict():
    """Print the verdict line and return the process exit status."""
    total = len(_CHECKS)
    failed = [nm for (nm, ok) in _CHECKS if not ok]
    if failed:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


def main():
    note("verification program for the weighted Szeged tree-minimizer census")
    note("python %s" % sys.version.split()[0])
    W = check_witnesses()
    check_p3_k3()
    trees = check_tree_census(11)
    check_small_witnesses(trees)
    check_witness_structures()
    classes = check_full_census_to_8(8)
    check_class_counts(classes)
    check_shard_arithmetic(classes)
    check_covering_census_9(classes[8], trees)
    for n in (10, 11):
        check_sparse_census(trees, n, K_SPARSE[n])
    check_witness_optimality(W, trees)
    check_table_arithmetic()
    scope_report()
    print("")
    note("elapsed %.1f s" % (time.time() - _T0))
    return verdict()


if __name__ == "__main__":
    sys.exit(main())

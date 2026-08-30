#!/usr/bin/env python3
"""
verify.py -- re-derives every quantity claimed in paper.tex for

    "The Clow--Zaguia cop bound is tight at independence number four: f(4) = 3"

from the objects PRINTED IN THAT PAPER: the graph6 string of the witness G, the
adjacency rule G is defined by, the four-part clique cover, the three-cell
dominating set, the six-case robber strategy, and the permutation phi that
identifies G with the Latin square graph of the order-4 Klein-group table.

Python 3.9+, STANDARD LIBRARY ONLY (itertools and sys are all that is imported),
no external data file, and no floating point anywhere: every decision below is
made on Python integers or on bitmask set operations.

Prints one `PASS <name> <detail>` line per check, then
    VERDICT: ALL <n> CHECKS PASS
and exits 0 if and only if every check passed. What it does NOT cover is printed
as `NOT RE-RUN:` lines just above the verdict, and repeated in REVIEW_NOTE.md's
`## Scope`.
"""

import sys
from itertools import combinations

# ---------------------------------------------------------------------------
# THE OBJECTS, COPIED FROM THE PAPER
# ---------------------------------------------------------------------------
# paper.tex, Section 2: the graph6 string of G (backslashes escaped for Python).
GRAPH6 = "O?]uf@vmuy\\o\\vmzZmf\\o"

# paper.tex, Section 2: the four broken diagonals Q_t = {(i, i+t mod 4)},
# written with 1-based cell coordinates exactly as in the paper.
CLIQUE_COVER = [
    [(1, 1), (2, 2), (3, 3), (4, 4)],
    [(1, 2), (2, 3), (3, 4), (4, 1)],
    [(1, 3), (2, 4), (3, 1), (4, 2)],
    [(1, 4), (2, 1), (3, 2), (4, 3)],
]

# paper.tex, Section 2: the dominating set.
DOMINATING_SET = [(1, 1), (2, 2), (3, 3)]

# paper.tex, Section 4: the isomorphism phi from the Latin square graph of the
# order-4 Klein-group table L(i,j) = i XOR j onto G. Vertex k of either graph is
# the cell (1 + k // 4, 1 + k % 4).
PHI = [0, 5, 10, 15, 11, 14, 1, 4, 13, 8, 7, 2, 6, 3, 12, 9]

# ---------------------------------------------------------------------------
# CHECK BOOKKEEPING
# ---------------------------------------------------------------------------
_lines = []
_failed = 0


def ck(name, ok, detail=""):
    global _failed
    if ok:
        _lines.append("PASS %s%s" % (name, (" " + detail) if detail else ""))
    else:
        _failed += 1
        _lines.append("FAIL %s%s" % (name, (" " + detail) if detail else ""))
    return bool(ok)


def note(text):
    _lines.append(text)


# ---------------------------------------------------------------------------
# GRAPHS AS BITMASKS
# ---------------------------------------------------------------------------
def pc(x):
    return bin(x).count("1")


def mk(n, edges):
    adj = [0] * n
    for a, b in edges:
        if a == b:
            raise ValueError("loop")
        adj[a] |= 1 << b
        adj[b] |= 1 << a
    return adj


def edges_of(n, adj):
    return [(a, b) for a in range(n) for b in range(a + 1, n) if adj[a] >> b & 1]


def complement(n, adj):
    full = (1 << n) - 1
    return [(full ^ adj[v]) & ~(1 << v) for v in range(n)]


def g6_encode(n, adj):
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if adj[i] >> j & 1 else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = v * 2 + b
        out.append(chr(v + 63))
    return "".join(out)


def g6_decode(s):
    n = ord(s[0]) - 63
    bits = []
    for chx in s[1:]:
        v = ord(chx) - 63
        for k in (5, 4, 3, 2, 1, 0):
            bits.append(v >> k & 1)
    adj = [0] * n
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            idx += 1
    return n, adj


def bfs_ecc(n, adj, s):
    dist = {s: 0}
    frontier = [s]
    while frontier:
        nxt = []
        for x in frontier:
            m = adj[x]
            while m:
                y = (m & -m).bit_length() - 1
                m &= m - 1
                if y not in dist:
                    dist[y] = dist[x] + 1
                    nxt.append(y)
        frontier = nxt
    return dist


def diameter(n, adj):
    """-1 if disconnected."""
    d = 0
    for s in range(n):
        dist = bfs_ecc(n, adj, s)
        if len(dist) < n:
            return -1
        d = max(d, max(dist.values()))
    return d


def girth(n, adj):
    """Length of a shortest cycle, or None if acyclic. Exact BFS from every vertex."""
    best = None
    for s in range(n):
        dist = {s: 0}
        parent = {s: None}
        frontier = [s]
        while frontier:
            nxt = []
            for x in frontier:
                m = adj[x]
                while m:
                    y = (m & -m).bit_length() - 1
                    m &= m - 1
                    if y == parent[x]:
                        continue
                    if y in dist:
                        c = dist[x] + dist[y] + 1
                        if best is None or c < best:
                            best = c
                    else:
                        dist[y] = dist[x] + 1
                        parent[y] = x
                        nxt.append(y)
            frontier = nxt
    return best


def alpha(n, adj):
    """Exact independence number by branch and bound on bitmasks."""
    best = [0]

    def rec(cand, size):
        if size + pc(cand) <= best[0]:
            return
        if cand == 0:
            best[0] = max(best[0], size)
            return
        v = max((u for u in range(n) if cand >> u & 1),
                key=lambda u: pc(adj[u] & cand))
        rec(cand & ~adj[v] & ~(1 << v), size + 1)
        rec(cand & ~(1 << v), size)

    rec((1 << n) - 1, 0)
    return best[0]


def omega(n, adj):
    return alpha(n, complement(n, adj))


def chi(n, adj):
    """Exact chromatic number, DSATUR branch and bound."""
    best = [n]
    color = [-1] * n

    def rec(ncol):
        if ncol >= best[0]:
            return
        u, bs = -1, -1
        for v in range(n):
            if color[v] == -1:
                s = len({color[w] for w in range(n)
                         if adj[v] >> w & 1 and color[w] != -1})
                if s > bs:
                    bs, u = s, v
        if u == -1:
            best[0] = min(best[0], ncol)
            return
        forb = {color[w] for w in range(n) if adj[u] >> w & 1 and color[w] != -1}
        for c in range(min(ncol + 1, best[0] - 1) + 1):
            if c in forb:
                continue
            color[u] = c
            rec(max(ncol, c + 1))
            color[u] = -1

    rec(0)
    return best[0]


def gamma_exact(n, adj, cap=None):
    """Exact domination number by increasing exhaustive search."""
    clo = [adj[v] | 1 << v for v in range(n)]
    full = (1 << n) - 1
    top = n if cap is None else cap
    for k in range(1, top + 1):
        for S in combinations(range(n), k):
            m = 0
            for v in S:
                m |= clo[v]
            if m == full:
                return k
    return None


def dominates(n, adj, cells):
    clo = 0
    for v in cells:
        clo |= adj[v] | 1 << v
    return clo == (1 << n) - 1


def is_hole(n, adj, mask):
    """True iff the induced subgraph on `mask` is a chordless cycle of length >= 5."""
    vs = [v for v in range(n) if mask >> v & 1]
    if len(vs) < 5:
        return False
    for v in vs:
        if pc(adj[v] & mask) != 2:
            return False
    seen = {vs[0]}
    st = [vs[0]]
    while st:
        x = st.pop()
        m = adj[x] & mask
        while m:
            y = (m & -m).bit_length() - 1
            m &= m - 1
            if y not in seen:
                seen.add(y)
                st.append(y)
    return len(seen) == len(vs)


def perfect_by_spgt_scan(n, adj):
    """Strong Perfect Graph Theorem criterion, brute force over ALL odd subsets of
    size >= 5: G is perfect iff it has no induced odd hole and no induced odd
    antihole. Returns (bool, witness_or_None, scanned_count)."""
    co = complement(n, adj)
    scanned = 0
    for mask in range(1 << n):
        k = pc(mask)
        if k < 5 or k % 2 == 0:
            continue
        scanned += 1
        if is_hole(n, adj, mask):
            return False, ("odd hole", k, mask), scanned
        if is_hole(n, co, mask):
            return False, ("odd antihole", k, mask), scanned
    return True, None, scanned


def two_cops_win(n, adj):
    """Exact backward induction (least fixpoint / attractor) over the full state
    space cop-pair x robber-vertex x turn-parity, with passing legal for both
    sides (closed neighbourhoods). Returns True iff two cops have a winning
    strategy, i.e. iff there is a starting pair from which every robber choice
    is already a cop-win state."""
    clo = [adj[v] | 1 << v for v in range(n)]
    nb = [[u for u in range(n) if clo[v] >> u & 1] for v in range(n)]

    def IDX(a, b, r):
        return (a * n + b) * n + r

    sz = n * n * n
    winC = bytearray(sz)          # cops to move
    winR = bytearray(sz)          # robber to move
    cnt = [0] * sz                # unresolved robber successors
    stack = []
    for a in range(n):
        for b in range(n):
            for r in range(n):
                i = IDX(a, b, r)
                cnt[i] = len(nb[r])
                if r == a or r == b:
                    winC[i] = 1
                    winR[i] = 1
                    stack.append((i, 0))
                    stack.append((i, 1))
    while stack:
        i, par = stack.pop()
        r = i % n
        b = (i // n) % n
        a = i // (n * n)
        if par:                    # robber-to-move state is a cop win
            for ap in nb[a]:       # predecessors: cops moved a->ap, b->bp
                for bp in nb[b]:
                    j = IDX(ap, bp, r)
                    if not winC[j]:
                        winC[j] = 1
                        stack.append((j, 0))
        else:                      # cop-to-move state is a cop win
            for rp in nb[r]:
                j = IDX(a, b, rp)
                if not winR[j]:
                    cnt[j] -= 1
                    if cnt[j] == 0:
                        winR[j] = 1
                        stack.append((j, 1))
    for a in range(n):
        for b in range(a, n):
            if all(winC[IDX(a, b, r)] for r in range(n)):
                return True
    return False


def is_isomorphism(n, A, B, phi):
    if sorted(phi) != list(range(n)):
        return False
    for a in range(n):
        for b in range(n):
            if bool(A[a] >> b & 1) != bool(B[phi[a]] >> phi[b] & 1):
                return False
    return True


def srg_parameters(n, adj):
    """(k, lambda, mu) if strongly regular, else None."""
    ks = {pc(adj[v]) for v in range(n)}
    if len(ks) != 1:
        return None
    lam, mu = set(), set()
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            c = pc(adj[a] & adj[b])
            (lam if adj[a] >> b & 1 else mu).add(c)
    if len(lam) != 1 or len(mu) != 1:
        return None
    return (ks.pop(), lam.pop(), mu.pop())


def local_graph_profile(n, adj, v):
    """(order, size, connected, degree multiset) of the subgraph induced on N(v)."""
    nbrs = [u for u in range(n) if adj[v] >> u & 1]
    mask = 0
    for u in nbrs:
        mask |= 1 << u
    size = sum(1 for x, y in combinations(nbrs, 2) if adj[x] >> y & 1)
    seen = {nbrs[0]}
    st = [nbrs[0]]
    while st:
        x = st.pop()
        for y in nbrs:
            if y not in seen and adj[x] >> y & 1:
                seen.add(y)
                st.append(y)
    degs = tuple(sorted(pc(adj[x] & mask) for x in nbrs))
    return (len(nbrs), size, len(seen) == len(nbrs), degs)


# ---------------------------------------------------------------------------
# THE FAMILY K_m x K_m = complement of the m x m rook's graph
# ---------------------------------------------------------------------------
def tensor_complement_of_rook(rows, cols):
    """Cells (i,j), 1-based, indexed cols*(i-1)+(j-1); (i,j) ~ (k,l) iff i != k and j != l."""
    n = rows * cols
    E = [(a, b) for a in range(n) for b in range(a + 1, n)
         if (a // cols) != (b // cols) and (a % cols) != (b % cols)]
    return n, mk(n, E)


def idx(m, cell):
    i, j = cell
    return m * (i - 1) + (j - 1)


def cell(m, k):
    return (1 + k // m, 1 + k % m)


# --- the six-case robber strategy, exactly as printed in the paper ----------
def _least_outside(m, banned):
    for v in range(1, m + 1):
        if v not in banned:
            return v
    return None


def strategy_reply(m, r, d1, d2):
    """The robber's reply, by the paper's named cases. `r` is the robber's cell,
    `d1`, `d2` are the cops' cells AFTER their move. Returns (case, cell)."""
    x, y = r
    p, q = d1
    s, t = d2
    thr1 = (p != x and q != y)          # cop 1 is adjacent to the robber
    thr2 = (s != x and t != y)
    if not thr1 and not thr2:
        return ("0", r)                  # pass
    if thr1 and thr2:
        if p == s:
            v = _least_outside(m, {y, q, t})
            return ("1a", None if v is None else (p, v))
        if q != t:
            return ("1b", (p, t))
        u = _least_outside(m, {x, p, s})
        return ("1c", None if u is None else (u, q))
    if thr2 and not thr1:                # relabel so that cop 1 is the threat
        p, q, s, t = s, t, p, q
    if s == x and t != y:
        if t != q:
            return ("2a", (p, t))
        u = _least_outside(m, {x, p})
        return ("2a", None if u is None else (u, t))
    if t == y and s != x:
        if s != p:
            return ("2b", (s, q))
        v = _least_outside(m, {y, q})
        return ("2b", None if v is None else (s, v))
    return ("impossible", None)          # s == x and t == y means cop 2 sits on r


def strategy_initial(m, c1, c2):
    """The robber's opening choice, by the paper's rule."""
    a, b = c1
    c, d = c2
    if a == c:
        e = _least_outside(m, {b, d})
        return None if e is None else (a, e)
    if b == d:
        e = _least_outside(m, {a, c})
        return None if e is None else (e, b)
    return (a, d)


def shares_line(u, v):
    return u != v and (u[0] == v[0] or u[1] == v[1])


def safe(r, c1, c2):
    """The paper's invariant: the robber shares a line with every cop (hence is
    adjacent to none) and sits on no cop."""
    return shares_line(r, c1) and shares_line(r, c2)


def audit_strategy(m, named=True):
    """Replays the paper's Lemma over EVERY safe state and EVERY cop move.

    Returns a dict with the counts and the first violation found (or None).
    `named=True` checks the specific cell the paper's case analysis names;
    `named=False` only asks that some safe reply exists."""
    cells = [cell(m, k) for k in range(m * m)]
    nbrs = {}
    for u in cells:
        nbrs[u] = [u] + [v for v in cells if v[0] != u[0] and v[1] != u[1]]
    lines = {u: [v for v in cells if shares_line(u, v)] for u in cells}
    res = {"pairs": 0, "safe_states": 0, "cop_moves": 0,
           "init_bad": None, "bad": None, "cases": {}}
    for c1 in cells:
        for c2 in cells:
            res["pairs"] += 1
            r0 = strategy_initial(m, c1, c2) if named else None
            if named:
                if r0 is None or not safe(r0, c1, c2):
                    if res["init_bad"] is None:
                        res["init_bad"] = (c1, c2, r0)
            elif not any(safe(r, c1, c2) for r in cells):
                if res["init_bad"] is None:
                    res["init_bad"] = (c1, c2, None)
    for c1 in cells:
        for c2 in cells:
            for r in lines[c1]:
                if not shares_line(r, c2):
                    continue
                res["safe_states"] += 1
                for d1 in nbrs[c1]:
                    for d2 in nbrs[c2]:
                        res["cop_moves"] += 1
                        if r == d1 or r == d2:
                            if res["bad"] is None:
                                res["bad"] = ("captured", c1, c2, r, d1, d2)
                            continue
                        if named:
                            case, rp = strategy_reply(m, r, d1, d2)
                            res["cases"][case] = res["cases"].get(case, 0) + 1
                            ok = (rp is not None and rp in nbrs[r]
                                  and safe(rp, d1, d2))
                        else:
                            ok = any(safe(rp, d1, d2) for rp in nbrs[r])
                        if not ok:
                            if res["bad"] is None:
                                res["bad"] = ("no reply", c1, c2, r, d1, d2)
    return res


# ---------------------------------------------------------------------------
# CONTROL OBJECTS
# ---------------------------------------------------------------------------
def petersen():
    E = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
         (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
         (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]
    return 10, mk(10, E)


def heawood():
    """Incidence graph of the Fano plane: point i ~ line j iff i in {j, j+1, j+3} mod 7."""
    E = []
    for i in range(7):
        for j in range(7):
            if i in ((j) % 7, (j + 1) % 7, (j + 3) % 7):
                E.append((i, 7 + j))
    return 14, mk(14, E)


def shrikhande():
    """Cayley graph on Z_4 x Z_4 with connection set +-(1,0), +-(0,1), +-(1,1)."""
    S = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)]
    E = set()
    for a in range(16):
        i, j = divmod(a, 4)
        for di, dj in S:
            b = ((i + di) % 4) * 4 + ((j + dj) % 4)
            E.add((min(a, b), max(a, b)))
    return 16, mk(16, sorted(E))


def latin_square_graph(table):
    """Vertices = the m^2 cells; two distinct cells adjacent iff they share a row,
    share a column, or carry the same symbol."""
    m = len(table)
    n = m * m
    E = []
    for a in range(n):
        for b in range(a + 1, n):
            i, j = divmod(a, m)
            k, l = divmod(b, m)
            if i == k or j == l or table[i][j] == table[k][l]:
                E.append((a, b))
    return n, mk(n, E)


def cycle(n):
    return n, mk(n, [(i, (i + 1) % n) for i in range(n)])


def is_bipartite(n, adj):
    col = {}
    for s in range(n):
        if s in col:
            continue
        col[s] = 0
        st = [s]
        while st:
            x = st.pop()
            m = adj[x]
            while m:
                y = (m & -m).bit_length() - 1
                m &= m - 1
                if y not in col:
                    col[y] = 1 - col[x]
                    st.append(y)
                elif col[y] == col[x]:
                    return False
    return True


# ===========================================================================
# 1. THE WITNESS, FROM THE PRINTED ADJACENCY RULE AND THE PRINTED graph6
# ===========================================================================
N, G = tensor_complement_of_rook(4, 4)
E = edges_of(N, G)

ck("graph6-encode",
   g6_encode(N, G) == GRAPH6,
   "the rule-generated graph encodes to the printed string %s" % GRAPH6)

dn, dG = g6_decode(GRAPH6)
ck("graph6-decode",
   dn == N and set(edges_of(dn, dG)) == set(E),
   "the printed string decodes to the same %d labelled edges (both directions checked)"
   % len(E))

ck("order-size-regular",
   N == 16 and len(E) == 72 and {pc(G[v]) for v in range(N)} == {9},
   "n=16, m=72, 9-regular")

nonedges = [(a, b) for a in range(N) for b in range(a + 1, N) if not G[a] >> b & 1]
lines_only = all((a // 4 == b // 4) or (a % 4 == b % 4) for a, b in nonedges)
ck("nonedges-are-exactly-the-lines",
   len(nonedges) == 48 and lines_only,
   "48 = 8 lines x C(4,2) non-edges, every one of them a same-row or same-column pair")

diam = diameter(N, G)
ck("connected-diameter-2", diam == 2, "diameter = 2, so G is connected")

# ===========================================================================
# 2. alpha = 4, omega = 4, theta = 4
# ===========================================================================
a_G = alpha(N, G)
ck("alpha-equals-4", a_G == 4, "exact maximum independent set size = 4")

max_ind = []
for S in combinations(range(N), 4):
    if all(not G[a] >> b & 1 for a, b in combinations(S, 2)):
        max_ind.append(S)
rows = [tuple(range(4 * i, 4 * i + 4)) for i in range(4)]
cols = [tuple(4 * i + j for i in range(4)) for j in range(4)]
ck("maximum-independent-sets-are-the-8-lines",
   sorted(max_ind) == sorted(rows + cols),
   "all %d independent 4-sets are exactly the 4 rows and 4 columns" % len(max_ind))

w_G = omega(N, G)
ck("omega-equals-4", w_G == 4, "exact maximum clique size = 4")

cover_idx = [[idx(4, c) for c in part] for part in CLIQUE_COVER]
parts_are_cliques = all(all(G[a] >> b & 1 for a, b in combinations(part, 2))
                        for part in cover_idx)
ck("clique-cover-parts-are-cliques",
   parts_are_cliques,
   "each of the 4 printed broken diagonals induces a K_4")

flat = [v for part in cover_idx for v in part]
ck("clique-cover-is-a-partition",
   sorted(flat) == list(range(16)),
   "the 4 parts partition all 16 cells, so theta <= 4")

theta_G = chi(N, complement(N, G))
ck("theta-equals-4",
   theta_G == 4 and theta_G >= a_G,
   "exact chi(complement of G) = 4 = theta(G), and theta >= alpha forces equality")

ck("alpha-equals-theta",
   a_G == theta_G == 4,
   "alpha = theta = 4, the identity the weak perfect graph theorem forces on a perfect graph")

# ===========================================================================
# 3. PERFECTION, BY TWO INDEPENDENT ROUTES
# ===========================================================================
ok_p, why_p, scanned = perfect_by_spgt_scan(N, G)
ck("perfect-strong-pgt-scan",
   ok_p,
   "no induced odd hole and no induced odd antihole among all %d odd vertex "
   "subsets of size >= 5 (all 2^16 subsets enumerated)" % scanned)

# Route 2: rook(4,4) IS the line graph of K_{4,4}, hence perfect by Koenig, hence
# its complement G is perfect by the weak perfect graph theorem.
rook = complement(N, G)
k44_edges = [(i, j) for i in range(4) for j in range(4)]        # edge (row i, col j)
line_adj = [0] * 16
for a in range(16):
    for b in range(16):
        if a != b and (k44_edges[a][0] == k44_edges[b][0]
                       or k44_edges[a][1] == k44_edges[b][1]):
            line_adj[a] |= 1 << b
ck("rook-is-the-line-graph-of-K44",
   line_adj == rook,
   "cell (i,j) <-> the edge of K_{4,4} joining row-vertex i to column-vertex j; "
   "sharing a line <-> sharing an endpoint, so complement(G) = L(K_{4,4}) and G "
   "is perfect by Koenig plus the weak PGT")

# ===========================================================================
# 4. c(G) <= 3
# ===========================================================================
dom_idx = [idx(4, c) for c in DOMINATING_SET]
ck("printed-dominating-set-dominates",
   dominates(N, G, dom_idx),
   "{(1,1),(2,2),(3,3)} leaves 0 cells undominated, so c(G) <= gamma(G) <= 3")

g_G = gamma_exact(N, G, cap=3)
ck("gamma-equals-3",
   g_G == 3,
   "no 1-cell or 2-cell dominating set exists (all 136 subsets checked), so gamma = 3 exactly")

# ===========================================================================
# 5. c(G) >= 3: THE SIX-CASE ROBBER STRATEGY, REPLAYED EXHAUSTIVELY
# ===========================================================================
aud = audit_strategy(4, named=True)
ck("opening-safe-cell-exists-for-every-cop-pair",
   aud["init_bad"] is None,
   "the printed opening rule returns a safe cell for all %d ordered cop pairs, "
   "equal pairs included" % aud["pairs"])

ck("capture-impossible-from-a-safe-cell",
   aud["bad"] is None or aud["bad"][0] != "captured",
   "over all %d safe states and all %d cop moves out of them, no cop move ever "
   "lands on the robber" % (aud["safe_states"], aud["cop_moves"]))

ck("named-case-restores-the-invariant",
   aud["bad"] is None,
   "for every safe state and every cop move the case the paper names returns a "
   "legal robber move to a cell that is safe again; case counts %s"
   % ({k: aud["cases"][k] for k in sorted(aud["cases"])}))

ck("robber-strategy-survives-forever",
   aud["init_bad"] is None and aud["bad"] is None,
   "initialisation + invariance + no-capture, so two cops never catch the robber: c(G) >= 3")

tcw = two_cops_win(N, G)
ck("two-cops-lose-exact-solver",
   tcw is False,
   "independent exact backward induction over all 16^3 x 2 states returns "
   "twoCopsWin = 0")

ck("cop-number-equals-3",
   (not tcw) and g_G == 3,
   "3 <= c(G) <= 3, so c(G) = 3 = alpha(G) - 1")

# ===========================================================================
# 6. f(4) = 3
# ===========================================================================
ck("f4-lower-bound-3",
   diam == 2 and ok_p and a_G == 4 and (not tcw),
   "G is connected, perfect, alpha(G) = 4 and c(G) = 3, so any valid f(4) is >= 3")

ck("f4-equals-3",
   a_G == 4 and (not tcw) and g_G == 3,
   "with the target paper's own theorem (connected perfect, alpha >= 4 => c < alpha) "
   "giving f(4) <= 3, the witness forces f(4) = 3")

# ===========================================================================
# 7. CONTROLS -- BOTH POLARITIES, EACH LABELLED WITH ITS OWN OBJECT
# ===========================================================================
pn, P = petersen()
p_perf, p_why, _ = perfect_by_spgt_scan(pn, P)
ck("control-petersen-forced-positive",
   alpha(pn, P) == 4 and gamma_exact(pn, P, cap=3) == 3
   and two_cops_win(pn, P) is False and p_perf is False
   and p_why[0] == "odd hole" and p_why[1] == 5,
   "Petersen (n=10): alpha=4, gamma=3, twoCopsWin=0 -- so the solver CAN return "
   "'two cops lose' -- and perfection is correctly refused on an induced C_5")

hn, H = heawood()
ck("control-heawood-forced-positive",
   is_bipartite(hn, H) and alpha(hn, H) == 7 and girth(hn, H) == 6
   and {pc(H[v]) for v in range(hn)} == {3}
   and gamma_exact(hn, H, cap=4) == 4 and two_cops_win(hn, H) is False,
   "Heawood (n=14): bipartite hence perfect, alpha=7, girth 6, 3-regular, "
   "gamma=4, twoCopsWin=0 so c >= 3 (the published value is c = 3). This is the "
   "q=2 case of the projective-plane family in the paper's scope section, and it "
   "checks that family's two Aigner-Fromme hypotheses (girth >= 5, min degree "
   "q+1) on a concrete member")

sn, SH = shrikhande()
CSH = complement(sn, SH)
csh_perf, csh_why, _ = perfect_by_spgt_scan(sn, CSH)
ck("control-complement-of-shrikhande-published-triple",
   alpha(sn, CSH) == 3 and chi(sn, SH) == 4
   and gamma_exact(sn, CSH, cap=3) == 3 and two_cops_win(sn, CSH) is False
   and csh_perf is False,
   "complement of Shrikhande (n=16): alpha=3, theta=4, gamma=3, twoCopsWin=0 "
   "(so c=3) and NOT perfect -- exactly the triple Char, Maniya and Pradhan "
   "publish, reproduced here")

c4n, C4 = cycle(4)
c4_perf, _, _ = perfect_by_spgt_scan(c4n, C4)
ck("anticontrol-C4",
   c4_perf is True and alpha(c4n, C4) == 2 and two_cops_win(c4n, C4) is True,
   "C_4 is perfect with alpha=2 and two cops win, i.e. c = alpha, NOT alpha-1; "
   "labelled an ANTI-control so it is never used as the forced positive for the "
   "c = alpha-1 shape")

n9, G9 = tensor_complement_of_rook(3, 3)
aud3 = audit_strategy(3, named=True)
ck("silence-control-K3xK3",
   two_cops_win(n9, G9) is True and (aud3["init_bad"] is not None or aud3["bad"] is not None),
   "K_3 x K_3 (n=9): two cops WIN, and the same six-case strategy FAILS there -- "
   "cases 1a and 1c need a fourth row and column. It must fail, because m_3 = 10 "
   "means no 9-vertex graph has cop number 3")

n12, G12 = tensor_complement_of_rook(4, 3)
p12, _, _ = perfect_by_spgt_scan(n12, G12)
ck("silence-control-12-vertex-family-member",
   n12 == 12 and len(edges_of(n12, G12)) == 36 and p12 is True
   and alpha(n12, G12) == 4 and two_cops_win(n12, G12) is True,
   "the complement of L(K_{4,3}) (n=12, m=36) is perfect with alpha=4 but is "
   "2-cop-win, so it does NOT lower the order 16 and the minimum order stays "
   "bracketed 10 <= min <= 16")

# ===========================================================================
# 8. IDENTITY GUARDS AGAINST THE COSPECTRAL TWIN
# ===========================================================================
srg_G = srg_parameters(N, G)
srg_CSH = srg_parameters(sn, CSH)
ck("identity-guard-cospectral-parameters",
   srg_G == (9, 4, 6) and srg_CSH == (9, 4, 6) and alpha(N, G) != alpha(sn, CSH),
   "G and the complement of Shrikhande are BOTH srg(16,9,4,6), hence cospectral, "
   "yet alpha is 4 and 3 respectively, so they are not isomorphic and confusing "
   "them would invalidate everything")

ck("identity-guard-local-graphs",
   srg_parameters(N, rook) == (6, 2, 2) and srg_parameters(sn, SH) == (6, 2, 2)
   and all(local_graph_profile(N, rook, v)[1:] == (6, False, (2,) * 6) for v in range(N))
   and all(local_graph_profile(sn, SH, v)[1:] == (6, True, (2,) * 6) for v in range(sn)),
   "both srg(16,6,2,2); the local graph of complement(G) = rook(4,4) is 2K_3 "
   "(2-regular on 6, DISCONNECTED) at every vertex while Shrikhande's is C_6 "
   "(2-regular on 6, CONNECTED) -- the invariant that tells the twins apart")

# ===========================================================================
# 9. THE PRIOR-ART ATTRIBUTION IS ITSELF CHECKED
# ===========================================================================
klein = [[i ^ j for j in range(4)] for i in range(4)]
z4 = [[(i + j) % 4 for j in range(4)] for i in range(4)]
ln_k, LK = latin_square_graph(klein)
ln_z, LZ = latin_square_graph(z4)

ck("latin-square-graph-order-4-degrees",
   ln_k == 16 and len(edges_of(ln_k, LK)) == 72
   and {pc(LK[v]) for v in range(16)} == {9}
   and ln_z == 16 and {pc(LZ[v]) for v in range(16)} == {9},
   "both order-4 Latin square graphs are 9-regular on 16 vertices with 72 edges, "
   "as srg(16,9,4,6) requires")

ck("witness-is-a-latin-square-graph-of-order-4",
   is_isomorphism(16, LK, G, PHI),
   "the permutation phi printed in the paper is an isomorphism from the Latin "
   "square graph of the Klein-group table L(i,j) = i XOR j onto G, so the object "
   "and its cop number 3 belong to Ahirwar, Bonato, Gittins, Huang, Marbach and "
   "Zaidman and must be cited")

ck("the-other-order-4-latin-square-is-the-twin",
   all(local_graph_profile(16, complement(16, LZ), v)[2] is True for v in range(16))
   and alpha(16, LZ) == 3,
   "the Z_4 table instead gives the complement of Shrikhande (local graph C_6, "
   "connected) with alpha = 3, so it is NOT our witness -- the polarity control "
   "on the identification")

# ===========================================================================
# 10. THE FAMILY: A LOWER BOUND ONLY, AND SAID SO
# ===========================================================================
for m in (5, 6):
    nm, Gm = tensor_complement_of_rook(m, m)
    covm = [[idx(m, (i, 1 + (i - 1 + t) % m)) for i in range(1, m + 1)]
            for t in range(m)]
    cov_ok = (sorted(v for part in covm for v in part) == list(range(m * m))
              and all(all(Gm[a] >> b & 1 for a, b in combinations(part, 2))
                      for part in covm))
    audm = audit_strategy(m, named=True)
    ck("family-K%dxK%d-has-cop-number-3" % (m, m),
       diameter(nm, Gm) == 2 and alpha(nm, Gm) == m and cov_ok
       and dominates(nm, Gm, [idx(m, c) for c in DOMINATING_SET])
       and audm["init_bad"] is None and audm["bad"] is None,
       "n=%d: connected, alpha = theta = %d (broken-diagonal cover verified), the "
       "SAME 3-cell dominating set works, and the SAME six-case strategy closes, "
       "so c = 3 -- which gives 3 <= f(%d) <= %d, a LOWER bound only, never f(%d) = 3"
       % (nm, m, m, m - 1, m))

# ===========================================================================
# SCOPE
# ===========================================================================
note("")
note("NOT RE-RUN: no 3-cop solver was run. c(G) <= 3 is established from the "
     "printed dominating set (c <= gamma) and, independently, from the target "
     "paper's own theorem; the exact solver above decides only the 2-cop game.")
note("NOT RE-RUN: perfection of K_m x K_m for m >= 5 was NOT brute-forced -- the "
     "strong-PGT scan is exponential in the order, so only the 9-, 12- and "
     "16-vertex graphs were scanned. For m >= 5 perfection is by theorem "
     "(line graph of a bipartite graph, plus the weak PGT), not by test.")
note("NOT RE-RUN: no census of connected perfect graphs was performed, so the "
     "MINIMUM ORDER of a connected perfect graph with alpha = 4 and cop number 3 "
     "is not determined here; it is only bracketed 10 <= min <= 16, the lower "
     "bound being a published value and the upper bound this witness.")
note("NOT RE-RUN: nothing here computes f(k) for any k >= 5. The family gives "
     "3 <= f(k) <= k-1 for k >= 4 and this program checks only the lower half at "
     "k = 5 and k = 6. Problem 5.2 itself is NOT solved.")
note("NOT RE-RUN: the family check above stops at m = 6; m = 7 and beyond rest on "
     "the case analysis in the paper, which uses only m >= 4, and were not replayed.")
note("NOT RE-RUN: no literature was fetched. Whether the value c(K_n x K_n) = 3 "
     "is already proved in print -- in particular in Neufeld and Nowakowski, "
     "Discrete Math. 186 (1998) 253-268, which we could not read -- is a "
     "bibliographic question no program can settle.")
note("")

for ln in _lines:
    print(ln)

npass = sum(1 for ln in _lines if ln.startswith("PASS "))
if _failed:
    print("VERDICT: %d CHECK(S) FAILED" % _failed)
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % npass)
sys.exit(0)

#!/usr/bin/env python3
"""Verification program for "The Saturation Number of K_3 u 2K_2" (Problem 3).

Standard library only; no data file is read and nothing is downloaded.  Run as
    python3 verify.py
Exit status 0 iff every check passes.  Progress goes to stderr, checks to stdout.

--------------------------------------------------------------------------------
TAKEN FROM THE PAPER (inputs to be tested, never used as evidence for themselves)
--------------------------------------------------------------------------------
  * the exhibited 13-vertex witnesses, as graph6 strings:
        B_13        = "L~`@?a?_C?O?_?"
        K_6 u I_7   = "L~~w??????????"
    together with the paper's descriptions of them (B_13 = a K_4 core u_1..u_4
    with one pendant vertex at each u_i and five further pendants at u_1;
    degree sequence 9,4,4,4,1^9) and the claim that each has 15 edges;
  * the exhibited 7-vertex extremal graph "FkOow" (claimed 9 edges) and the
    paper's description of G_7;
  * the claimed values sat(n,F) = C(n,2) for 1<=n<=6, n+2 for 7<=n<=12,
    15 for n>=13, and the table of Lemma 5 (9,10,11,12,13,14 for v=7..12);
  * the statement of Lemma 3 (the dichotomy specialised from Proposition 3.4
    of Lin-He-Xu), whose hypotheses and conclusion are recomputed on every
    census representative below;
  * the three enumeration figures quoted in the paper: 1,008,629 isolate-free
    graphs with 1..14 edges, 740,226 of them with exactly 14 edges, and the
    25,937,431 extension candidates of order at least 14 obtained from those;
  * the saturation criterion stated inside the proof of Lemma 5.

--------------------------------------------------------------------------------
DERIVED HERE (every line below is recomputed from scratch)
--------------------------------------------------------------------------------
  1. graph6 decoding of each exhibited string, its order, size, degree sequence,
     and a graph6 re-encoding that must reproduce the published string verbatim;
  2. an independent reconstruction of each witness from the paper's prose
     description, and a canonical-form test that the decoded and reconstructed
     graphs are isomorphic;
  3. F-freeness of each witness by an exhaustive 7-subset search, and
     F-saturation by adding every single nonedge and searching for K_3 u 2K_2 in
     the result (the load-bearing check: 63 nonedges for each 13-vertex witness);
  4. non-isomorphism of the two 13-vertex witnesses (canonical forms), which is
     the paper's sharpness claim at n = 13;
  5. negative controls: deleting any one edge of either witness destroys
     saturation, so the exhibited edge counts are locally minimal;
  6. the constructions of Lemma 2 (B_n for 8<=n<=20, K_6 u I_{n-6} for
     7<=n<=20, G_7) are built from the prose and tested for saturation, giving
     the upper bound sat(n) <= min{n+2,15};
  7. self-tests of the machinery: the F-predicate on F itself and on F minus an
     edge; agreement of the fast F-predicate with a brute-force oracle and of
     the paper's saturation criterion with the definition, over a corpus of
     23,289 graphs containing both verdicts; canonical forms validated against
     an independent Polya count (156 classes on 6 vertices, per-edge-count
     vector) and under random relabellings;
  8. sat(n) for 1<=n<=6 by exhausting all labelled graphs on n<=6 vertices:
     the F-saturated graphs are exactly the complete graphs;
  9. a self-generated census, with isomorphism rejection, of every isolate-free
     F-free graph of order <= 13 with at most 14 edges (no catalogue file);
     its unpruned level counts are checked against Polya counts for the first
     10 of the 14 levels (see the closing scope line), and every
     representative, padded with isolated vertices to each order n <= 13, is
     tested for saturation.  This yields sat(n) for 7<=n<=12 (the table of
     Lemma 5), the complete list of extremal graphs for those orders, and the
     lower bound sat(13) > 14 which with the two witnesses gives sat(13) = 15;
 10. the catalogue figures 1,008,629 and 740,226 recomputed by Polya counting,
     and with them the 25,937,431 extension candidates, obtained by summing the
     three extension multiplicities over the by-order breakdown of the 14-edge
     layer (that breakdown is printed, so the sum is reproducible by hand);
     these are counts of the layers, not saturation tests over them;
 11. Lemma 3 tested on every census representative and, exhaustively, on every
     labelled graph of order at most 6 (which is where the lemma's conclusion
     forces H to live, and which includes K_6 itself): both hypotheses are
     recomputed and wherever they hold the conclusion (two disjoint triangles,
     or a disjoint K_2 and K_4) must hold too, with a non-vacuity gate on each
     of the two families;
 12. an integrity check on the deduplication key that the whole lower bound
     rests on: for every census representative the canonical certificate is
     decoded and shown to be that graph relabelled by the ordering canon()
     selected, so equal certificates force isomorphism and the census cannot
     silently lose an isomorphism class;
 13. cross-validation of the census pipeline against the exhaustive labelled
     search on orders 2..5, where both are available.

The line beginning "NOT RE-RUN:" states exactly which part of the paper's
census is outside this program's budget and is therefore not reverified.
"""
import sys, itertools

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    sys.stdout.flush()
    return bool(ok)


def bits(m):
    while m:
        b = m & -m
        yield b.bit_length() - 1
        m ^= b


def g6_decode(s):
    """graph6 string -> (n, adj) with adj a list of neighbour bitmasks."""
    d = [ord(c) - 63 for c in s]
    for x in d:
        if x < 0 or x > 63:
            raise ValueError("bad graph6 character")
    if d[0] <= 62:
        n = d[0]
        pos = 1
    elif d[1] <= 62:
        n = (d[1] << 12) | (d[2] << 6) | d[3]
        pos = 4
    else:
        raise ValueError("graph6 order too large")
    need = (n * (n - 1) // 2 + 5) // 6
    body = d[pos:]
    if len(body) != need:
        raise ValueError("graph6 body length %d, expected %d" % (len(body), need))
    stream = []
    for x in body:
        for k in (5, 4, 3, 2, 1, 0):
            stream.append((x >> k) & 1)
    adj = [0] * n
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if stream[idx]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            idx += 1
    if any(stream[n * (n - 1) // 2:]):
        raise ValueError("graph6 padding bits nonzero")
    return n, adj


def g6_encode(n, adj):
    """(n, adj) -> graph6 string; inverse of g6_decode."""
    if n > 62:
        raise ValueError("order too large for this encoder")
    stream = []
    for j in range(1, n):
        for i in range(j):
            stream.append(1 if (adj[i] >> j) & 1 else 0)
    while len(stream) % 6:
        stream.append(0)
    out = [chr(63 + n)]
    for k in range(0, len(stream), 6):
        v = 0
        for b in stream[k:k + 6]:
            v = (v << 1) | b
        out.append(chr(63 + v))
    return "".join(out)


def edges_of(n, adj):
    return [(i, j) for i in range(n) for j in bits(adj[i]) if j > i]


def has_edge_in(adj, mask):
    for u in bits(mask):
        if adj[u] & mask:
            return True
    return False


def nu_ge_2(adj, mask):
    """True iff the induced subgraph on `mask` has two vertex-disjoint edges."""
    for u in bits(mask):
        for v in bits(adj[u] & mask):
            if v > u and has_edge_in(adj, mask & ~((1 << u) | (1 << v))):
                return True
    return False


def triangles(n, adj):
    for i in range(n):
        for j in bits(adj[i]):
            if j <= i:
                continue
            for k in bits(adj[i] & adj[j]):
                if k > j:
                    yield (i, j, k)


def contains_F(n, adj):
    """True iff G contains K_3 u 2K_2 as a (not necessarily induced) subgraph.

    Definitional: a triangle plus two further edges, all five parts pairwise
    vertex-disjoint.  Equivalent to: some triangle T has nu(G-T) >= 2.
    """
    full = (1 << n) - 1
    for (i, j, k) in triangles(n, adj):
        if nu_ge_2(adj, full & ~((1 << i) | (1 << j) | (1 << k))):
            return True
    return False


def contains_F_bruteforce(n, adj):
    """Independent slow oracle: search all 7-subsets and all splits of them."""
    full = list(range(n))
    for S in itertools.combinations(full, 7):
        for T in itertools.combinations(S, 3):
            a, b, c = T
            if not (adj[a] >> b & 1 and adj[a] >> c & 1 and adj[b] >> c & 1):
                continue
            rest = [v for v in S if v not in T]
            for pair in itertools.combinations(rest, 2):
                other = [v for v in rest if v not in pair]
                if adj[pair[0]] >> pair[1] & 1 and adj[other[0]] >> other[1] & 1:
                    return True
    return False


def has_disjoint_K3_K2(n, adj, mask):
    """True iff the induced subgraph on `mask` contains vertex-disjoint copies
    of K_3 and K_2."""
    for (i, j, k) in triangles(n, adj):
        t = (1 << i) | (1 << j) | (1 << k)
        if t & ~mask:
            continue
        if has_edge_in(adj, mask & ~t):
            return True
    return False


def dichotomy_hypothesis(n, adj):
    """Hypothesis (2) of Lemma 3: for every vertex z, H - z contains
    vertex-disjoint copies of K_2 and K_3.  (Hypothesis (1), that H contains no
    disjoint K_2,K_2,K_3, is exactly F-freeness and is checked separately.)"""
    full = (1 << n) - 1
    for z in range(n):
        if not has_disjoint_K3_K2(n, adj, full & ~(1 << z)):
            return False
    return True


def two_disjoint_triangles(n, adj):
    ts = [((1 << i) | (1 << j) | (1 << k)) for (i, j, k) in triangles(n, adj)]
    for a in range(len(ts)):
        for b in range(a + 1, len(ts)):
            if not (ts[a] & ts[b]):
                return True
    return False


def disjoint_K2_and_K4(n, adj):
    full = (1 << n) - 1
    for (i, j, k) in triangles(n, adj):
        for l in bits(adj[i] & adj[j] & adj[k]):
            if l <= k:
                continue
            q = (1 << i) | (1 << j) | (1 << k) | (1 << l)
            if has_edge_in(adj, full & ~q):
                return True
    return False


def dichotomy_conclusion(n, adj):
    """The conclusion of Lemma 3: two vertex-disjoint triangles, or
    vertex-disjoint copies of K_2 and K_4."""
    return two_disjoint_triangles(n, adj) or disjoint_K2_and_K4(n, adj)


def is_saturated(n, adj):
    """Definitional test: G is F-free and G+e contains F for every nonedge e."""
    if contains_F(n, adj):
        return False
    for x in range(n):
        for y in range(x + 1, n):
            if adj[x] >> y & 1:
                continue
            b = list(adj)
            b[x] |= 1 << y
            b[y] |= 1 << x
            if not contains_F(n, b):
                return False
    return True


def saturating_nonedge_paper(n, adj, x, y):
    """The paper's criterion for the nonedge xy being saturating, stated inside
    the proof of Lemma 5:
    either G-{x,y} has disjoint K_3 and K_2, or some common neighbour z of
    x and y has nu(G-{x,y,z}) >= 2."""
    full = (1 << n) - 1
    rest = full & ~((1 << x) | (1 << y))
    for (i, j, k) in triangles(n, adj):
        t = (1 << i) | (1 << j) | (1 << k)
        if t & ~rest:
            continue
        if has_edge_in(adj, rest & ~t):
            return True
    for z in bits(adj[x] & adj[y]):
        if nu_ge_2(adj, full & ~((1 << x) | (1 << y) | (1 << z))):
            return True
    return False


def is_saturated_paper(n, adj):
    """Same decision, computed through the paper's stated criterion."""
    if contains_F(n, adj):
        return False
    for x in range(n):
        for y in range(x + 1, n):
            if adj[x] >> y & 1:
                continue
            if not saturating_nonedge_paper(n, adj, x, y):
                return False
    return True


def refine_colors(n, adj):
    """Isomorphism-invariant equitable colouring; returns list of colours."""
    col = [bin(adj[v]).count("1") for v in range(n)]
    while True:
        keys = []
        for v in range(n):
            keys.append((col[v], tuple(sorted(col[u] for u in bits(adj[v])))))
        order = sorted(set(keys))
        newcol = [order.index(keys[v]) for v in range(n)]
        if newcol == col:
            return col
        col = newcol


def canon(n, adj, want_perm=False):
    """Canonical certificate: the lexicographically least column sequence over
    all vertex orderings consistent with the invariant colouring.

    With want_perm, also return the vertex ordering that attains the certificate;
    check_certificate_pins_class uses it to verify that the certificate really is
    this graph relabelled, hence that canon() cannot merge two nonisomorphic
    graphs (which is what makes the census incapable of losing a class).

    Search pruning: at each position only candidates from the current colour
    cell are tried, only those attaining the least column value can start a
    lexicographically least completion, and a candidate whose twin has already
    been tried at this node is skipped (transposing twins is an automorphism
    fixing every other vertex, so it cannot change the certificate)."""
    col = refine_colors(n, adj)
    cells = {}
    for v in range(n):
        cells.setdefault(col[v], []).append(v)
    cellseq = [list(cells[c]) for c in sorted(cells)]
    best = [None]
    bestperm = [None]

    def rec(pos, perm, ci, remaining, cur):
        if pos == n:
            t = tuple(cur)
            if best[0] is None or t < best[0]:
                best[0] = t
                if want_perm:
                    bestperm[0] = list(perm)
            return
        while not remaining[ci]:
            ci += 1
        scored = []
        tried = []
        for v in remaining[ci]:
            if any(adj[v] & ~(1 << w) == adj[w] & ~(1 << v) for w in tried):
                continue          # a twin of v was already tried at this node
            tried.append(v)
            c = 0
            for i, u in enumerate(perm):
                if adj[v] >> u & 1:
                    c |= 1 << i
            scored.append((c, v))
        m = min(s[0] for s in scored)
        cur.append(m)
        if best[0] is None or tuple(cur) <= best[0][:len(cur)]:
            for c, v in scored:
                if c != m:
                    continue
                remaining[ci].remove(v)
                perm.append(v)
                rec(pos + 1, perm, ci, remaining, cur)
                perm.pop()
                remaining[ci].append(v)
        cur.pop()

    rec(0, [], 0, cellseq, [])
    if want_perm:
        return (n, best[0]), bestperm[0]
    return (n, best[0])


def decode_certificate(cert):
    """Rebuild the graph named by a canonical certificate.  Column j of the
    certificate is the neighbourhood of the j-th vertex among the earlier ones,
    so the certificate determines a graph outright."""
    n, cols = cert
    adj = [0] * n
    for j in range(n):
        for i in bits(cols[j]):
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return n, adj


def certificate_pins_class(n, adj):
    """True iff the certificate returned by canon(), decoded, is exactly this
    graph relabelled by the ordering canon() found.  If this holds then
    canon(G) == canon(H) forces G and H isomorphic (both are isomorphic to the
    decoded graph), so deduplicating by canon can only ever merge isomorphic
    graphs -- the census cannot silently lose an isomorphism class."""
    cert, perm = canon(n, adj, True)
    m, dadj = decode_certificate(cert)
    if m != n or perm is None or sorted(perm) != list(range(n)):
        return False
    pos = [0] * n
    for i, v in enumerate(perm):
        pos[v] = i
    b = [0] * n
    for v in range(n):
        for w in bits(adj[v]):
            b[pos[v]] |= 1 << pos[w]
    return b == dadj


def add_edges(n, pairs):
    adj = [0] * n
    for (i, j) in pairs:
        if i == j:
            raise ValueError("loop")
        adj[i] |= 1 << j
        adj[j] |= 1 << i
    return n, adj


def build_B(n):
    """B_n of Lemma 2 (n >= 8): core K_4 on u1..u4, one leaf on each u_i,
    the remaining n-8 leaves on u_1."""
    if n < 8:
        raise ValueError("B_n defined for n >= 8")
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    pairs += [(i, 4 + i) for i in range(4)]
    pairs += [(0, v) for v in range(8, n)]
    return add_edges(n, pairs)


def build_K6_plus_isolates(n):
    """K_6 u I_{n-6}."""
    if n < 6:
        raise ValueError("need n >= 6")
    return add_edges(n, [(i, j) for i in range(6) for j in range(i + 1, 6)])


def build_G7():
    """The 7-vertex graph of Lemma 2: A={0,1,2} a triangle, x_i=3+i with
    a_i x_i, and s=6 joined to every x_i."""
    pairs = [(0, 1), (0, 2), (1, 2)]
    pairs += [(i, 3 + i) for i in range(3)]
    pairs += [(6, 3 + i) for i in range(3)]
    return add_edges(7, pairs)


def isomorphic(g, h):
    return canon(*g) == canon(*h)


def degseq(n, adj):
    return sorted((bin(a).count("1") for a in adj), reverse=True)


def extensions(n, adj, cap):
    """The three ways to grow an isolate-free graph by one edge (the reverse
    step in the proof of Lemma 5): join two nonadjacent existing vertices,
    attach a new leaf,
    or adjoin a disjoint K_2.  Order is capped at `cap`."""
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            if not (adj[i] >> j & 1):
                b = list(adj)
                b[i] |= 1 << j
                b[j] |= 1 << i
                out.append((n, b))
    if n + 1 <= cap:
        for i in range(n):
            b = list(adj) + [0]
            b[i] |= 1 << n
            b[n] |= 1 << i
            out.append((n + 1, b))
    if n + 2 <= cap:
        out.append((n + 2, list(adj) + [1 << (n + 1), 1 << n]))
    return out


def pad(n, adj, k):
    return n + k, list(adj) + [0] * k


def levels_iter(cap, emax, prune_F):
    """Yield (e, {canon: (n, adj)}) for every isolate-free graph of order <= cap
    with e edges, up to isomorphism; if prune_F, drop graphs containing F (valid
    because containing F is preserved by adding edges, so no F-saturated graph
    has an F-containing subgraph)."""
    cur = {canon(2, [2, 1]): (2, [2, 1])}
    yield 1, cur
    for e in range(2, emax + 1):
        nxt = {}
        for (n, adj) in cur.values():
            for (m, b) in extensions(n, adj, cap):
                if prune_F and contains_F(m, b):
                    continue
                c = canon(m, b)
                if c not in nxt:
                    nxt[c] = (m, b)
        cur = nxt
        yield e, cur


def census(cap, emax, report=None):
    """Generate every isolate-free F-free graph of order <= cap with at most
    emax edges, up to isomorphism, and record for each order n <= cap the least
    number of edges of an n-vertex F-saturated graph (with isolated vertices
    allowed, obtained by padding a smaller core).

    Returns (levels, best) where levels[e] is the number of representatives with
    e edges and best[n] = (min_edges, graph6) or None.
    """
    best = {}
    levels = {}
    per_n = {}
    iso_tested = [0]
    iso_sat = []
    integrity = {"pin_tested": 0, "pin_bad": [],
                 "dich_hyp": 0, "dich_bad": []}
    for e, cur in levels_iter(cap, emax, True):
        levels[e] = len(cur)
        for (m, adj) in cur.values():
            # Integrity of the deduplication key: the certificate must pin the
            # isomorphism class, or the census could drop a class silently.
            integrity["pin_tested"] += 1
            if not certificate_pins_class(m, adj):
                integrity["pin_bad"].append(g6_encode(m, adj))
            # Lemma 3 on this representative: it is F-free by construction
            # (hypothesis 1); if hypothesis 2 also holds, the conclusion must.
            if dichotomy_hypothesis(m, adj):
                integrity["dich_hyp"] += 1
                if not dichotomy_conclusion(m, adj):
                    integrity["dich_bad"].append(g6_encode(m, adj))
            for k in range(0, cap - m + 1):
                n, b = pad(m, adj, k)
                per_n[n] = per_n.get(n, 0) + 1
                if k > 0:
                    iso_tested[0] += 1
                if not is_saturated_paper(n, b):
                    continue
                if k > 0:
                    iso_sat.append((n, e, m, g6_encode(n, b)))
                if n not in best:
                    best[n] = (e, [])
                if best[n][0] == e:
                    best[n][1].append(g6_encode(n, b))
        if report:
            report(e, levels[e], best)
    return levels, best, per_n, (iso_tested[0], iso_sat), integrity


def _partitions(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield []
        return
    for p in range(min(n, maxp), 0, -1):
        for rest in _partitions(n - p, p):
            yield [p] + rest


def classes_by_edges(k, emax):
    """Burnside/Polya count: number of isomorphism classes of graphs on exactly
    k vertices with e edges, for e = 0..emax.  Independent of the generator."""
    from math import factorial, gcd
    tot = [0] * (emax + 1)
    fk = factorial(k)
    for lam in _partitions(k):
        mult = {}
        for p in lam:
            mult[p] = mult.get(p, 0) + 1
        size = fk
        for p, m in mult.items():
            size //= (p ** m) * factorial(m)
        cyc = []
        for a in range(len(lam)):
            for b in range(a + 1, len(lam)):
                g = gcd(lam[a], lam[b])
                cyc += [lam[a] * lam[b] // g] * g
        for p in lam:
            if p % 2 == 0:
                cyc.append(p // 2)
                cyc += [p] * ((p - 2) // 2)
            else:
                cyc += [p] * ((p - 1) // 2)
        poly = [0] * (emax + 1)
        poly[0] = 1
        for l in cyc:
            if l > emax:
                continue
            for i in range(emax - l, -1, -1):
                if poly[i]:
                    poly[i + l] += poly[i]
        for e in range(emax + 1):
            tot[e] += size * poly[e]
    return [t // fk for t in tot]


def isolate_free_classes(kmax, emax):
    """B[k][e] = classes of isolate-free graphs on exactly k vertices with e
    edges, via B(k,e) = A(k,e) - A(k-1,e)."""
    A = {0: [1] + [0] * emax}
    for k in range(1, kmax + 1):
        A[k] = classes_by_edges(k, emax)
    return {k: [A[k][e] - A[k - 1][e] for e in range(emax + 1)]
            for k in range(1, kmax + 1)}


def finish():
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


G6_B13 = "L~`@?a?_C?O?_?"
G6_K6I7 = "L~~w??????????"
G6_G7 = "FkOow"
PAPER_DEGSEQ_B13 = [9, 4, 4, 4] + [1] * 9
PAPER_SAT = dict([(n, n * (n - 1) // 2) for n in range(1, 7)] +
                 [(n, n + 2) for n in range(7, 13)] +
                 [(n, 15) for n in range(13, 21)])
PAPER_TABLE = {7: 9, 8: 10, 9: 11, 10: 12, 11: 13, 12: 14}
PAPER_CATALOGUE_TOTAL = 1008629
PAPER_14EDGE_REPS = 740226
PAPER_EXTENSION_CANDIDATES = 25937431
CAP = 13
EMAX = 14
POLYA_EMAX = 10                # range of the independent Polya cross-check of
                               # the generator's own unpruned level counts


def check_F_predicate():
    """The subgraph predicate itself: F = K_3 u 2K_2 must be recognised, and
    must not be recognised in F minus any edge."""
    nF, aF = add_edges(7, [(0, 1), (0, 2), (1, 2), (3, 4), (5, 6)])
    # Not just "7 vertices and 5 edges": the literal edge list is checked to be
    # K_3 u 2K_2 up to isomorphism -- three components of sizes 3, 2, 2 and
    # exactly one triangle, which together with the degrees pins the graph.
    seen, comps = 0, []
    for v in range(nF):
        if seen >> v & 1:
            continue
        stack, c = [v], 0
        seen |= 1 << v
        while stack:
            u = stack.pop()
            c += 1
            for w in bits(aF[u] & ~seen):
                seen |= 1 << w
                stack.append(w)
        comps.append(c)
    ntri = len(list(triangles(nF, aF)))
    ck("F_definition_7_vertices_5_edges",
       nF == 7 and len(edges_of(nF, aF)) == 5 and
       degseq(nF, aF) == [2, 2, 2, 1, 1, 1, 1] and
       sorted(comps, reverse=True) == [3, 2, 2] and ntri == 1,
       "n=%d e=%d degrees=%s components=%s triangles=%d"
       % (nF, len(edges_of(nF, aF)), degseq(nF, aF),
          sorted(comps, reverse=True), ntri))
    ck("F_contains_itself", contains_F(nF, aF) and contains_F_bruteforce(nF, aF))
    survivors = []
    for (i, j) in edges_of(nF, aF):
        b = list(aF)
        b[i] &= ~(1 << j)
        b[j] &= ~(1 << i)
        if contains_F(7, b) or contains_F_bruteforce(7, b):
            survivors.append((i, j))
    ck("F_minus_any_edge_is_F_free", not survivors,
       "%d deletions, %d of them still contain F" % (5, len(survivors)))
    n8, a8 = add_edges(8, [(0, 1), (0, 2), (1, 2), (3, 4), (4, 5)])
    ck("K3_plus_P3_is_F_free", not contains_F(n8, a8) and
       not contains_F_bruteforce(n8, a8), "P_3 has no 2-matching")


def check_witness(name, g6, order, nedges, builder, dseq=None):
    n, adj = g6_decode(g6)
    ck(name + "_graph6_roundtrip", g6_encode(n, adj) == g6, g6)
    E = edges_of(n, adj)
    ck(name + "_order_and_size", n == order and len(E) == nedges,
       "n=%d e=%d" % (n, len(E)))
    if dseq is not None:
        ck(name + "_degree_sequence", degseq(n, adj) == dseq,
           "computed %s" % (degseq(n, adj),))
    cd, cb = canon(n, adj), canon(*builder)
    ck(name + "_matches_paper_construction", cd == cb,
       "decoded certificate %s %s built-from-prose certificate %s"
       % (cd[1], "==" if cd == cb else "!=", cb[1]))
    freeF = not contains_F_bruteforce(n, adj)
    ck(name + "_is_F_free", freeF,
       "7-subset exhaustive search: %s"
       % ("no K_3 u 2K_2 found" if freeF else "a K_3 u 2K_2 WAS found"))
    nonedges = n * (n - 1) // 2 - len(E)
    bad = []
    for x in range(n):
        for y in range(x + 1, n):
            if adj[x] >> y & 1:
                continue
            b = list(adj)
            b[x] |= 1 << y
            b[y] |= 1 << x
            if not contains_F_bruteforce(n, b):
                bad.append((x, y))
    ck(name + "_every_nonedge_creates_F", not bad,
       "%d nonedges, %d failures" % (nonedges, len(bad)))
    ck(name + "_is_F_saturated", is_saturated(n, adj))
    return n, adj, E


def check_sharpness(B, K):
    """The two exhibited 13-vertex extremal graphs, and the negative controls."""
    nB, aB = B
    nK, aK = K
    ck("n13_two_witnesses_same_size",
       len(edges_of(nB, aB)) == len(edges_of(nK, aK)) == 15,
       "computed sizes %d and %d, paper claims 15 for both"
       % (len(edges_of(nB, aB)), len(edges_of(nK, aK))))
    ck("n13_witnesses_nonisomorphic", canon(nB, aB) != canon(nK, aK),
       "distinct canonical forms; degree sequences %s vs %s"
       % (degseq(nB, aB), degseq(nK, aK)))
    ck("n13_witnesses_isolated_vertex_counts_differ",
       degseq(nB, aB).count(0) == 0 and degseq(nK, aK).count(0) == 7,
       "isolated vertices: B_13 has %d (paper: 0), K_6 u I_7 has %d (paper: 7)"
       % (degseq(nB, aB).count(0), degseq(nK, aK).count(0)))
    for (label, nn, aa) in (("B13", nB, aB), ("K6I7", nK, aK)):
        survivors = []
        for (i, j) in edges_of(nn, aa):
            b = list(aa)
            b[i] &= ~(1 << j)
            b[j] &= ~(1 << i)
            if is_saturated(nn, b):
                survivors.append((i, j))
        ck(label + "_minus_any_edge_not_saturated", not survivors,
           "negative control over all %d single-edge deletions: %d still "
           "F-saturated%s" % (len(edges_of(nn, aa)), len(survivors),
                              (" " + str(survivors)) if survivors else ""))


def check_constructions():
    """Lemma 2: sat(n,F) <= min{n+2, 15} for n >= 7, from the two families."""
    e7, s7 = len(edges_of(*build_G7())), is_saturated(*build_G7())
    ck("G7_construction_saturated_9_edges", s7 and e7 == 9,
       "e(G_7)=%d, F-saturated=%s, so sat(7) <= %d" % (e7, s7, e7))
    rows = []
    for n in range(8, 21):
        g = build_B(n)
        e = len(edges_of(*g))
        rows.append((n, e, is_saturated(*g)))
    ck("B_n_saturated_with_n_plus_2_edges",
       all(e == n + 2 and s for (n, e, s) in rows),
       "n = 8..20: %d/%d have e(B_n)=n+2, %d/%d are F-saturated%s"
       % (sum(1 for (n, e, s) in rows if e == n + 2), len(rows),
          sum(1 for r in rows if r[2]), len(rows),
          "" if all(e == n + 2 and s for (n, e, s) in rows) else
          "; offenders " + str([r for r in rows if not (r[1] == r[0] + 2
                                                        and r[2])])))
    rows = []
    for n in range(7, 21):
        g = build_K6_plus_isolates(n)
        e = len(edges_of(*g))
        rows.append((n, e, is_saturated(*g)))
    ck("K6_plus_isolates_saturated_15_edges",
       all(e == 15 and s for (n, e, s) in rows),
       "n = 7..20: %d/%d have 15 edges, %d/%d are F-saturated%s"
       % (sum(1 for r in rows if r[1] == 15), len(rows),
          sum(1 for r in rows if r[2]), len(rows),
          "" if all(e == 15 and s for (n, e, s) in rows) else
          "; offenders " + str([r for r in rows if not (r[1] == 15
                                                        and r[2])])))
    # For each n, take the smallest F-saturated graph actually produced by the
    # two families and compare its edge count with the paper's min{n+2,15}.
    ub, ok = {}, True
    for n in range(7, 21):
        fam = [build_K6_plus_isolates(n)] + ([build_B(n)] if n >= 8 else
                                             [build_G7()])
        sizes = [len(edges_of(*g)) for g in fam if is_saturated(*g)]
        if not sizes:
            ok = False
            continue
        ub[n] = min(sizes)
        if ub[n] != min(n + 2, 15):
            ok = False
    ck("constructed_upper_bound_equals_min_n_plus_2_and_15", ok,
       "n=7..20: " + ",".join("%d:%d" % (n, ub[n]) for n in sorted(ub))[:120])


def all_labelled(n):
    """Every labelled graph on n vertices, as (n, adj)."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for m in range(1 << len(pairs)):
        adj = [0] * n
        for k, (i, j) in enumerate(pairs):
            if m >> k & 1:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
        yield n, adj


def check_canon():
    """Validate the canonical form against the Polya count and against random
    relabellings."""
    import random
    rnd = random.Random(20260824)
    ok = True
    for g6 in (G6_B13, G6_K6I7, G6_G7):
        n, adj = g6_decode(g6)
        c = canon(n, adj)
        for _ in range(25):
            p = list(range(n))
            rnd.shuffle(p)
            b = [0] * n
            for i in range(n):
                for j in bits(adj[i]):
                    b[p[i]] |= 1 << p[j]
            if canon(n, b) != c:
                ok = False
        if canon(n, adj) != c:
            ok = False
    ck("canon_invariant_under_relabelling", ok, "3 graphs x 25 permutations")
    byedges = {}
    for (n, adj) in all_labelled(6):
        byedges.setdefault(len(edges_of(n, adj)), set()).add(canon(n, adj))
    got = [len(byedges.get(e, ())) for e in range(16)]
    want = classes_by_edges(6, 15)
    ck("canon_class_counts_match_polya", got == want,
       "order 6, e=0..15: %s" % (got,))
    ck("canon_total_classes_order6", sum(got) == 156,
       "%d isomorphism classes on 6 vertices" % sum(got))


def check_small_orders():
    """n <= 6 exhaustively over all labelled graphs: the only F-saturated graph
    is K_n, so sat(n,F) = C(n,2)."""
    ok = True
    computed = {}
    counts = {}
    for n in range(1, 7):
        sat = [(m, adj) for (m, adj) in all_labelled(n) if is_saturated(m, adj)]
        forms = set(canon(m, adj) for (m, adj) in sat)
        kn = canon(*add_edges(n, [(i, j) for i in range(n)
                                  for j in range(i + 1, n)]))
        counts[n] = len(forms)
        if not sat or forms != {kn}:
            ok = False
            continue
        computed[n] = min(len(edges_of(m, adj)) for (m, adj) in sat)
    ck("small_orders_only_K_n_is_saturated", ok,
       "isomorphism classes of F-saturated graphs found, n=1..6: %s "
       "(each must be exactly 1, namely K_n)"
       % (",".join("%d:%d" % (n, counts[n]) for n in range(1, 7)),))
    agree = (len(computed) == 6 and
             all(computed[n] == PAPER_SAT[n] for n in range(1, 7)))
    ck("small_orders_min_edges_match_paper", agree,
       "computed " + ",".join("%d:%d" % (n, computed.get(n, -1))
                              for n in range(1, 7)) +
       " vs paper " + ",".join("%d:%d" % (n, PAPER_SAT[n])
                               for n in range(1, 7)))
    return computed


def corpus():
    """A test corpus: every isolate-free graph of order <= 9 with <= 10 edges
    (F-containing ones included), each padded out to order 9, plus seeded
    random graphs."""
    import random
    for e, cur in levels_iter(10, 11, False):
        for (m, adj) in cur.values():
            for k in range(0, 11 - m):
                yield pad(m, adj, k)
    seeds = [g6_decode(s) for s in (G6_B13, G6_K6I7, G6_G7)]
    seeds += [build_B(n) for n in range(8, 14)]
    seeds += [build_K6_plus_isolates(n) for n in range(7, 14)]
    seeds += [build_G7()]
    for (n, adj) in seeds:
        yield n, adj
        for i in range(n):
            for j in range(i + 1, n):
                b = list(adj)
                b[i] ^= 1 << j
                b[j] ^= 1 << i
                yield n, b
    rnd = random.Random(11235813)
    for _ in range(1500):
        n = rnd.randint(7, 11)
        p = rnd.choice([0.15, 0.25, 0.35, 0.5])
        adj = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if rnd.random() < p:
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
        yield n, adj


def check_verifier_logic():
    """Cross-validate the two F-containment routines and the paper's saturation
    criterion against the definition, on a corpus with both verdicts present."""
    tot = withF = satd = 0
    bad_c = bad_s = 0
    for (n, adj) in corpus():
        tot += 1
        c1 = contains_F(n, adj)
        c2 = contains_F_bruteforce(n, adj)
        if c1 != c2:
            bad_c += 1
        if c1:
            withF += 1
        s1 = is_saturated(n, adj)
        s2 = is_saturated_paper(n, adj)
        if s1 != s2:
            bad_s += 1
        if s1:
            satd += 1
    ck("corpus_is_nonvacuous", tot > 20000 and withF > 1000 and satd >= 20,
       "%d graphs, %d contain F, %d are F-saturated" % (tot, withF, satd))
    ck("contains_F_agrees_with_brute_force", bad_c == 0,
       "%d disagreements over %d graphs" % (bad_c, tot))
    ck("paper_saturation_criterion_agrees_with_definition", bad_s == 0,
       "%d disagreements over %d graphs" % (bad_s, tot))


def check_generator_exhaustive(pruned_levels, emax=POLYA_EMAX):
    """The unpruned generator must reproduce the Polya counts exactly, and the
    F-free pruning must account for the difference."""
    B = isolate_free_classes(CAP, emax)
    want = [None] + [sum(B[k][e] for k in range(1, CAP + 1))
                     for e in range(1, emax + 1)]
    got = {}
    withF = {}
    for e, cur in levels_iter(CAP, emax, False):
        got[e] = len(cur)
        withF[e] = sum(1 for (n, adj) in cur.values() if contains_F(n, adj))
    ok = all(got[e] == want[e] for e in range(1, emax + 1))
    ck("generator_matches_polya_counts", ok,
       "order <= %d, e=1..%d: %s" % (CAP, emax,
                                     [got[e] for e in range(1, emax + 1)]))
    ok2 = all(got[e] - withF[e] == pruned_levels[e]
              for e in range(1, emax + 1))
    ck("F_free_pruning_is_exact", ok2,
       "unpruned - F-containing vs pruned, e=1..%d: %s"
       % (emax, [(got[e] - withF[e], pruned_levels[e])
                 for e in range(1, emax + 1)]))


def check_catalogue_counts():
    """The two enumeration figures quoted in the paper, recomputed by Polya
    counting (no catalogue file is read)."""
    cat = 14                      # the paper's census runs over e = 1..14
    B = isolate_free_classes(2 * cat, cat)
    per_e = [sum(B[k][e] for k in range(1, 2 * cat + 1))
             for e in range(cat + 1)]
    total = sum(per_e[1:])
    ck("catalogue_total_1_008_629", total == PAPER_CATALOGUE_TOTAL,
       "computed %d isolate-free classes with 1..14 edges, paper says %d"
       % (total, PAPER_CATALOGUE_TOTAL))
    ck("catalogue_14_edge_reps_740_226", per_e[14] == PAPER_14EDGE_REPS,
       "computed %d, paper says %d" % (per_e[14], PAPER_14EDGE_REPS))
    # Independent confirmation of the Polya counter itself: for e <= 7 every
    # isolate-free graph has order <= 14, so explicit generation with cap 14
    # must reproduce the same all-orders counts.
    gen = {}
    for e, cur in levels_iter(2 * 7, 7, False):
        gen[e] = len(cur)
    ck("polya_counter_confirmed_by_explicit_generation",
       all(gen[e] == per_e[e] for e in range(1, 8)),
       "e=1..7 generated %s, Polya %s"
       % ([gen[e] for e in range(1, 8)], per_e[1:8]))
    print("INFO Polya counts of isolate-free classes, e=1..14: %s"
          % (per_e[1:],))
    # The paper's third enumeration figure, derived from the same table.  The
    # 15-edge layer is generated from the 14-edge representatives by the three
    # extensions of the proof of Lemma 5 -- join two nonadjacent existing
    # vertices (order unchanged, C(k,2)-14 ways), attach a new leaf (order k+1,
    # k ways), adjoin a disjoint K_2 (order k+2, one way) -- keeping only
    # candidates whose resulting order is at least 14, repetitions allowed.
    # This is the by-ORDER breakdown of the 14-edge layer, so it is printed too:
    # the summed row above does not determine it.
    by_order_14 = [(k, B[k][cat]) for k in range(1, 2 * cat + 1) if B[k][cat]]
    ext_total = sum(c * ((k * (k - 1) // 2 - cat if k >= 14 else 0)
                         + (k if k + 1 >= 14 else 0)
                         + (1 if k + 2 >= 14 else 0))
                    for (k, c) in by_order_14)
    print("INFO isolate-free 14-edge classes by order: %s" % (by_order_14,))
    ck("extension_candidates_25_937_431",
       ext_total == PAPER_EXTENSION_CANDIDATES,
       "computed %d one-edge / one-leaf / one-K_2 extensions of the %d "
       "fourteen-edge classes with resulting order >= 14, paper says %d"
       % (ext_total, per_e[cat], PAPER_EXTENSION_CANDIDATES))


def check_isolate_lemma(iso):
    """Lemma 4: an F-saturated graph of order >= 7 with an isolated vertex is
    K_6 u I_{n-6}, which has 15 edges.  Inside the census range (at most 14
    edges) there must therefore be no F-saturated graph whatsoever with an
    isolated vertex.  The complementary half, that K_6 u I_{n-6} really is
    saturated, is checked in check_constructions."""
    tested, sat = iso
    ck("lemma_isolate_none_saturated_below_15_edges",
       tested > 100000 and not sat,
       "%d graphs with at least one isolated vertex and at most %d edges tested,"
       " %d F-saturated (must be 0 since K_6 u I_m needs 15 edges)"
       % (tested, EMAX, len(sat)))


def check_dichotomy(integrity):
    """Lemma 3, which the paper deduces from Proposition 3.4 of Lin-He-Xu and
    then uses to prove Lemma 4.  Every census representative is F-free, i.e.
    satisfies hypothesis (1); hypothesis (2) is recomputed here, and wherever
    both hold the conclusion (two disjoint triangles, or a disjoint K_2 and
    K_4) must hold as well."""
    hyp = integrity["dich_hyp"]
    bad = list(integrity["dich_bad"])
    # Widen the test beyond the census: every labelled graph on at most six
    # vertices, F-containing ones included so that hypothesis (1) is really
    # tested, which brings in K_6 itself -- the graph Lemma 4 uses the
    # dichotomy to identify, and which has 15 > EMAX edges so the census
    # cannot see it.
    small_hyp, small_tot = 0, 0
    for k in range(1, 7):
        for (nn, aa) in all_labelled(k):
            small_tot += 1
            if not any(aa) or contains_F(nn, aa):
                continue                       # nonempty, hypothesis (1)
            if dichotomy_hypothesis(nn, aa):
                small_hyp += 1
                if not dichotomy_conclusion(nn, aa):
                    bad.append(g6_encode(nn, aa))
    ck("lemma_dichotomy_hypotheses_are_satisfiable", hyp > 0 and small_hyp > 0,
       "%d of the %d census representatives and %d of the %d labelled graphs of "
       "order <= 6 satisfy both hypotheses"
       % (hyp, integrity["pin_tested"], small_hyp, small_tot))
    ck("lemma_dichotomy_conclusion_holds_wherever_hypotheses_hold", not bad,
       "%d instances tested, %d satisfying both hypotheses but not the "
       "conclusion%s" % (hyp + small_hyp, len(bad),
                         (": " + " ".join(bad[:5])) if bad else ""))


def check_certificate_pins_class(integrity):
    """The census deduplicates by canonical certificate, so the lower bound
    sat(13) > 14 would be void if two nonisomorphic graphs could share a
    certificate.  For every representative the certificate is decoded and shown
    to be that very graph relabelled by the ordering canon() selected; hence
    equal certificates force isomorphism and no class can be lost."""
    bad = integrity["pin_bad"]
    ck("canon_certificate_pins_the_isomorphism_class", not bad,
       "%d census representatives; %d whose decoded certificate is not the "
       "graph itself relabelled%s"
       % (integrity["pin_tested"], len(bad),
          (": " + " ".join(bad[:5])) if bad else ""))


def check_census(levels, best, small, wit13, per_n):
    """Lemma 5's table and the values of the theorem for n <= 13."""
    table = dict((n, best[n][0]) for n in range(7, 13) if n in best)
    ck("census_table_orders_7_to_12", table == PAPER_TABLE,
       "computed %s, paper %s" % (sorted(table.items()),
                                  sorted(PAPER_TABLE.items())))
    # The lower bound below is only as good as the range the census actually
    # swept, so the range is read back off the run itself: every edge count
    # 1..EMAX must have been generated and every order 2..CAP tested.  (An
    # earlier version of this check ended in "and CAP >= 13 and EMAX >= 14",
    # a comparison between two module constants that no input can falsify.)
    spanned = sorted(levels) == list(range(1, EMAX + 1))
    orders = sorted(n for n in per_n if per_n[n] > 0)
    ck("census_covers_every_order_up_to_13",
       spanned and orders == list(range(2, CAP + 1)),
       "edge counts generated: %s (expected 1..%d); graphs tested per order: %s"
       % (sorted(levels), EMAX,
          ",".join("%d:%d" % (n, per_n.get(n, 0)) for n in range(2, CAP + 1))))
    # Independent cross-validation of the whole census pipeline (extension
    # generator, canonical dedup, padding, is_saturated_paper) against the
    # exhaustive labelled search of check_small_orders, on the orders where both
    # are available: 2..5.  Order 6 must be absent from the census because K_6
    # has 15 > EMAX edges.
    agree = all(n in best and best[n][0] == small[n] for n in range(2, 6))
    ck("census_agrees_with_exhaustive_search_on_small_orders",
       agree and 6 not in best,
       "census %s vs exhaustive %s; order 6 absent from census (K_6 has 15 > "
       "%d edges): %s"
       % (",".join("%d:%s" % (n, best[n][0] if n in best else None)
                   for n in range(2, 6)),
          ",".join("%d:%d" % (n, small[n]) for n in range(2, 6)),
          EMAX, 6 not in best))
    ck("census_no_saturated_graph_of_order_13_with_at_most_14_edges",
       13 not in best and per_n.get(13, 0) > 0,
       "%d isolate-free F-free representatives of order <= %d with <= %d edges; "
       "%d of them pad to order 13; F-saturated among those: %s"
       % (sum(levels.values()), CAP, EMAX, per_n.get(13, 0),
          "none" if 13 not in best else
          "%d with %d edges (%s)" % (len(best[13][1]), best[13][0],
                                     " ".join(best[13][1][:5]))))
    reverified = []
    for n in sorted(best):
        for g6 in best[n][1]:
            m, adj = g6_decode(g6)
            reverified.append(m == n and is_saturated(m, adj) and
                              len(edges_of(m, adj)) == best[n][0])
    ck("census_witnesses_reverified_by_definition", all(reverified),
       "%d extremal representatives re-tested with the definitional test"
       % len(reverified))
    forms = dict((n, set(canon(*g6_decode(g)) for g in best[n][1]))
                 for n in best)
    ck("G7_is_the_unique_extremal_graph_of_order_7",
       7 in forms and forms[7] == {canon(*g6_decode(G6_G7))},
       "%d extremal graph(s) of order 7 with 9 edges" % len(forms.get(7, ())))
    inset = all(canon(*build_B(n)) in forms.get(n, set())
                for n in range(8, 13))
    ck("B_n_is_extremal_for_orders_8_to_12", inset,
       "counts of extremal graphs: " +
       ",".join("%d:%d" % (n, len(forms.get(n, ()))) for n in range(7, 13)))
    computed = dict(small)
    for n in range(7, 13):
        if n in best:
            computed[n] = best[n][0]
    if 13 not in best:
        # lower bound sat(13) > EMAX from the census; upper bound = the edge
        # count actually measured on the two decoded, verified witnesses.
        computed[13] = min(wit13)
    ok = all(computed.get(n) == PAPER_SAT[n] for n in range(1, 14))
    ck("theorem_values_sat_n_for_1_to_13", ok,
       "computed " + ",".join("%d:%s" % (n, computed.get(n))
                              for n in range(1, 14)))
    return computed


def not_rerun(levels, polya_emax):
    """The closing scope line.  Every number in it that describes THIS run is
    taken from the run: the size of the top census level and the range of the
    independent Polya cross-check are passed in, never spelled out."""
    return (
    "NOT RE-RUN: this program is not the program that produced the paper's "
    "census; it is an independent and deliberately narrower re-verification "
    "written from the paper's text, and it reads no catalogue file.  The "
    "paper's census is exhaustive over McKay's catalogues of all 1,008,629 "
    "isolate-free graphs with 1..14 edges (orders up to 28) and over the "
    "25,937,431 one-edge / one-leaf / one-K_2 extensions of the 740,226 "
    "fourteen-edge representatives.  This program re-runs that census only for "
    "orders <= %d and <= %d edges (self-generated), which suffices for sat(n) "
    "with n <= 13; it does NOT re-run orders 14..28 of the 14-edge layer nor "
    "any part of the 15-edge layer, so the values sat(n)=15 for n >= 14 and the "
    "uniqueness of K_6 u I_{n-6} for n >= 14 are NOT reverified here.  All "
    "three enumeration figures (1,008,629, 740,226 and 25,937,431) ARE "
    "recomputed here by Polya counting, but that fixes only the SIZE of those "
    "layers, not the F-saturation test over them.  Finally, the independent "
    "Polya cross-check of this program's own unpruned generator covers "
    "e=1..%d of the %d levels; levels %d..%d, including the top level of %d "
    "representatives on which the bound sat(13) > 14 rests, are supported by "
    "the generator's extension argument together with the certificate-pinning "
    "and small-order cross-validation checks, not by a second independent "
    "count." % (CAP, EMAX, polya_emax, EMAX, polya_emax + 1, EMAX,
                levels[EMAX]))


def main():
    print("INFO verifying: The Saturation Number of K_3 u 2K_2 (Problem 3)")
    print("INFO F = K_3 u 2K_2; census range: order <= %d, edges <= %d"
          % (CAP, EMAX))
    check_F_predicate()
    check_canon()
    B = check_witness("B13", G6_B13, 13, 15, build_B(13), PAPER_DEGSEQ_B13)
    K = check_witness("K6I7", G6_K6I7, 13, 15, build_K6_plus_isolates(13),
                      [5] * 6 + [0] * 7)
    check_witness("G7", G6_G7, 7, 9, build_G7(), [3, 3, 3, 3, 2, 2, 2])
    check_sharpness((B[0], B[1]), (K[0], K[1]))
    check_constructions()
    check_verifier_logic()
    small = check_small_orders()
    sys.stderr.write("running census (order <= %d, edges <= %d)...\n"
                     % (CAP, EMAX))

    def report(e, cnt, best):
        sys.stderr.write("  e=%2d representatives=%7d best=%s\n"
                         % (e, cnt, sorted((n, v[0]) for n, v in best.items())))
        sys.stderr.flush()

    levels, best, per_n, iso, integrity = census(CAP, EMAX, report)
    print("INFO census: %d isolate-free F-free representatives of order <= %d "
          "with at most %d edges; %d padded graphs examined"
          % (sum(levels.values()), CAP, EMAX, sum(per_n.values())))
    print("INFO census levels by edge count: %s"
          % ([levels[e] for e in sorted(levels)],))
    check_generator_exhaustive(levels, POLYA_EMAX)
    check_certificate_pins_class(integrity)
    check_dichotomy(integrity)
    check_isolate_lemma(iso)
    check_catalogue_counts()
    wit13 = [len(edges_of(*g6_decode(g))) for g in (G6_B13, G6_K6I7)]
    computed = check_census(levels, best, small, wit13, per_n)
    for n in range(7, 13):
        if n in best:
            print("INFO all %d extremal graphs of order %d (%d edges), in this "
                  "program's own labelling: %s"
                  % (len(best[n][1]), n, best[n][0], " ".join(best[n][1])))
    print("INFO computed sat(n,F) for n=1..13: %s"
          % ([computed.get(n) for n in range(1, 14)],))
    print(not_rerun(levels, POLYA_EMAX))
    finish()


if __name__ == "__main__":
    main()

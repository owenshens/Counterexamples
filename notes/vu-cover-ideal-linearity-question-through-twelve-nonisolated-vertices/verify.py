#!/usr/bin/env python3
"""
verify.py -- referee verification for
  "Vu's Cover-Ideal Linearity Question Through Twelve Nonisolated Vertices"
  (settles Vu, Question 3.2, affirmatively for every bipartite graph with at
  most twelve nonisolated vertices).

Standard library only.  Exact integer and bit-mask arithmetic throughout; no
floating point, no randomness, no external data file, no tolerance parameter.

=====================================================================
VALUES TAKEN FROM THE PAPER (inputs; each is compared against an
independently recomputed value, never against itself)
=====================================================================
  * PAPER_TABLE      : all eleven rows of Table 1 (n, generated, retained),
                       n = 2..12, and its two column totals
                       635,786,845 generated / 465,766,775 retained;
  * PAPER_RETAINED_ABSTRACT = 465,766,775, the count quoted in the abstract
                       and in the proof of Proposition 5;
  * the encoding of Section 3 (a definition): for n <= 12 and
                       1 <= q <= floor(n/2), p = n-q, a nondecreasing tuple
                       (A_1,...,A_q) of nonempty subsets of [p] with union [p];
  * E(p,q) = C(2^p+q-2, q) and
    K(p,q) = sum_r (-1)^r C(p,r) C(2^(p-r)+q-2, q)      (Section 3 formulas);
  * equations (4) and (5) for alpha(G) and alpha(G-N[e]) (Section 3);
  * Vu's local formula (1)-(3): v_{P_e}(J^t) = c_e+(t-1)g_e,
    c_e = deg u + deg v - 2 + tau(G-N[e]), g_e = min(h(u),h(v)),
    h(u) = deg u + tau(G-N[u])                       -- CITED, not re-proved;
  * the statement of Proposition 5 (every edge of M(G) has an endpoint in a
    maximum independent set) and of Proposition 4 (the edge criterion).

=====================================================================
WHAT THIS PROGRAM DERIVES (the checks)
=====================================================================
 1. Table 1, three times and independently of the paper's inclusion-exclusion:
    (a) by the closed forms E,K; (b) by an unbounded-knapsack dynamic program
    over the 2^p-1 nonempty subsets that uses no inclusion-exclusion and no
    binomial coefficient, for every cell AND every transposed cell; (c) by
    prefix enumeration of every cell of every row, n = 2..12 inclusive.
 2. That the encoding means what the paper says: every retained encoding with
    n <= 9 is built as an explicit graph and shown bipartite (BFS 2-colouring),
    with no isolated vertex and with the declared sides really a bipartition;
    every discarded encoding is shown to have an isolated vertex AND to
    reappear in a smaller census cell once those vertices are deleted, so
    discarding loses no graph; every one of the 77,340 distinct labelled
    bipartite graphs without isolated vertices on n <= 7 vertices (reached
    through 169,656 (bipartition, edge set) pairs) is built from scratch and
    located in the census, and conversely every retained encoding is hit,
    which pins the labelling convention down.
 3. Equations (4) and (5), Lemma 3 (c_e = n-2-alpha(G-N[e]), h(u) >= tau with
    equality iff u lies in a maximum independent set, min_e g_e = tau) and the
    maximum-independent-set support rule, each against brute force over all
    2^n vertex subsets, for all 59,846 retained encodings with n <= 9; tau(G)
    and every tau(G-N[e]) are also recomputed by Koenig maximum matching, and
    equation (2) is evaluated on that alpha-free route before being compared
    with Lemma 3's closed form (otherwise the comparison is an arithmetic
    tautology that holds for every graph).
 4. THE LOAD-BEARING CHECK: zero violations of Proposition 5.  It is computed
    for every one of the 955,674 retained encodings with n <= 10 by the
    literal route of equations (4)-(5), and for the whole census up to n = 12
    by a transposed scan of 110,570,825 left-relabelling-orbit representatives
    (see the NOT-RE-RUN note printed at the end).  The two routes are compared
    graph by graph on 27 cells, including two n = 12 cells, plus two off-census
    shapes added only to exercise the transposed scanner at the widths q = 5
    and q = 6 that carry 98% of the census; and they are tied to a common
    closed-form count of labelled incidence matrices through their orbit
    sizes.  The short circuit inside the scanner that decides 104 of those 110
    million representatives is additionally compared, leaf by leaf, with the
    unabridged computation ON A FAMILY WHOSE VERDICT IS FALSE, since agreeing
    only where nothing violates would not test it at all.
 5. The theorem itself for n <= 9: min_e (c_e+(t-1)g_e) = min_e c_e +
    (t-1)tau(G) for t = 1..20, from Vu's local formula.
 6. Non-vacuity and sensitivity.  The condition of Proposition 5 is shown to
    FAIL for 4,381 general graphs on at most 6 vertices (so the violation
    counter is live and bipartiteness is doing work); 60 of them have no
    argmax edge at all meeting a maximum independent set, so Proposition 4's
    equivalence is tested with both sides false; 26,847 census graphs own an
    edge with neither endpoint in a maximum independent set, so the conclusion
    is a statement about WHICH edges maximise alpha(G-N[e]); and three
    deliberate corruptions -- argmin for argmax, "in every maximiser" for "in
    some maximiser" in the literal scan, and the same misreading inside the
    transposed scanner that carries the census -- are each shown to produce
    violations.

Exit status is 0 if and only if every check passes.
"""

import sys
from itertools import combinations, combinations_with_replacement, permutations

_RESULTS = []


def check(name, ok, detail=""):
    ok = bool(ok)
    _RESULTS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    sys.stdout.flush()
    return ok


def note(msg):
    print("NOTE " + str(msg))
    sys.stdout.flush()


def verdict():
    n = len(_RESULTS)
    bad = [nm for nm, ok in _RESULTS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


# =====================================================================
# INPUTS TAKEN FROM THE PAPER
# =====================================================================

PAPER_NMAX = 12

# Table 1: n -> (generated, retained).
PAPER_TABLE = {
    2: (1, 1),
    3: (3, 1),
    4: (13, 5),
    5: (43, 14),
    6: (235, 98),
    7: (1239, 522),
    8: (10659, 5472),
    9: (98439, 53733),
    10: (1428007, 895828),
    11: (23944527, 16033665),
    12: (610303679, 448777436),
}
PAPER_TOTAL_GENERATED = 635786845
PAPER_TOTAL_RETAINED = 465766775
# Abstract and proof of Proposition 5 quote the retained total again.
PAPER_RETAINED_ABSTRACT = 465766775

# Range of shapes (p,q) used by Section 3.
def cells(n):
    return [(n - q, q) for q in range(1, n // 2 + 1)]


def nCk(a, b):
    if b < 0 or a < b:
        return 0
    r = 1
    for i in range(b):
        r = r * (a - i) // (i + 1)
    return r


def E_closed(p, q):
    """Paper's E(p,q): nondecreasing q-tuples of nonempty subsets of [p]."""
    return nCk(2 ** p + q - 2, q)


def K_closed(p, q):
    """Paper's K(p,q): those whose union is all of [p]."""
    return sum((-1) ** r * nCk(p, r) * nCk(2 ** (p - r) + q - 2, q)
               for r in range(p + 1))


def dp_counts(p, q):
    """Count multisets of size q of nonempty subsets of [p], and those whose
    union is [p], by an unbounded knapsack over the 2^p-1 items.  Uses no
    inclusion-exclusion and no binomial coefficient, so it is independent of
    E_closed / K_closed.  dp[c][U] = number of multisets of size c, drawn from
    the items processed so far, with union U."""
    full = (1 << p) - 1
    dp = [[0] * (1 << p) for _ in range(q + 1)]
    dp[0][0] = 1
    for A in range(1, 1 << p):
        for c in range(1, q + 1):
            prev = dp[c - 1]
            cur = dp[c]
            for U in range(1 << p):
                v = prev[U]
                if v:
                    cur[U | A] += v
    return sum(dp[q]), dp[q][full]


def direct_multiset_count(p, q):
    """Brute-force count using itertools, for tiny cells only."""
    tot = 0
    cov = 0
    full = (1 << p) - 1
    for tup in combinations_with_replacement(range(1, 1 << p), q):
        tot += 1
        u = 0
        for a in tup:
            u |= a
        if u == full:
            cov += 1
    return tot, cov


def check_table_from_closed_forms():
    """Table 1 and its totals against the paper's own closed forms."""
    badg, badr = [], []
    tg = tr = 0
    for n in range(2, PAPER_NMAX + 1):
        g = sum(E_closed(p, q) for p, q in cells(n))
        r = sum(K_closed(p, q) for p, q in cells(n))
        tg += g
        tr += r
        if g != PAPER_TABLE[n][0]:
            badg.append((n, g, PAPER_TABLE[n][0]))
        if r != PAPER_TABLE[n][1]:
            badr.append((n, r, PAPER_TABLE[n][1]))
    check("table1_generated_column_all_11_rows", not badg,
          "mismatches=%s" % badg if badg else "rows n=2..12 reproduced")
    check("table1_retained_column_all_11_rows", not badr,
          "mismatches=%s" % badr if badr else "rows n=2..12 reproduced")
    check("table1_total_generated", tg == PAPER_TOTAL_GENERATED,
          "computed=%d paper=%d" % (tg, PAPER_TOTAL_GENERATED))
    check("table1_total_retained", tr == PAPER_TOTAL_RETAINED,
          "computed=%d paper=%d" % (tr, PAPER_TOTAL_RETAINED))
    check("abstract_retained_count_equals_table_total",
          tr == PAPER_RETAINED_ABSTRACT,
          "computed=%d abstract=%d" % (tr, PAPER_RETAINED_ABSTRACT))
    return tg, tr


def check_counts_by_dp():
    """Two count-independent recomputations of every cell of the census."""
    bad_dp, bad_direct = [], []
    ncell = ndirect = 0
    for n in range(2, PAPER_NMAX + 1):
        for p, q in cells(n):
            ncell += 1
            e, k = dp_counts(p, q)
            if (e, k) != (E_closed(p, q), K_closed(p, q)):
                bad_dp.append((p, q, e, k, E_closed(p, q), K_closed(p, q)))
            # the transposed shape, which the orbit scan enumerates
            ncell += 1
            et, kt = dp_counts(q, p)
            if (et, kt) != (E_closed(q, p), K_closed(q, p)):
                bad_dp.append((q, p, et, kt, E_closed(q, p), K_closed(q, p)))
            if E_closed(p, q) <= 300000:
                ndirect += 1
                if direct_multiset_count(p, q) != (e, k):
                    bad_direct.append((p, q))
    check("knapsack_dp_reproduces_E_and_K_every_cell_n_le_12", not bad_dp,
          "%d cells checked; mismatches=%s" % (ncell, bad_dp))
    check("itertools_multiset_count_agrees_on_small_cells", not bad_direct,
          "%d cells checked; mismatches=%s" % (ndirect, bad_direct))


def primal_cell(p, q, want_keys=False, controls=False):
    """The literal route.  Enumerates EVERY encoding of shape (p,q) exactly as
    Section 3 describes (nondecreasing tuples of nonempty subsets of [p]),
    keeps those with union [p], and evaluates equations (4) and (5) and the
    maximum-independent-set support rule on each one.

    stats = [generated, retained, violations of Proposition 5,
             violations under the argmin corruption,
             violations under the "in every maximiser" corruption,
             encodings whose data disagrees with an earlier encoding of the
             same left-relabelling orbit,
             encodings owning an edge with neither endpoint in any maximum
             independent set, encodings with |M(G)| >= 2]
    """
    fq = (1 << q) - 1
    fp = (1 << p) - 1
    pc = [bin(i).count("1") for i in range(1 << p)]
    bc = [bin(i).count("1") for i in range(1 << q)]
    A = [0] * q
    stats = [0, 0, 0, 0, 0, 0, 0, 0]
    keys = {} if want_keys else None

    def leaf(U):
        stats[0] += 1
        if U[fq] != fp:
            return
        stats[1] += 1
        alpha = -1
        ms = []
        for B in range(1 << q):
            v = bc[B] + p - pc[U[B]]
            if v > alpha:
                alpha = v
                ms = [B]
            elif v == alpha:
                ms.append(B)
        misR = misL = 0
        wR, wL = fq, fp
        for B in ms:
            misR |= B
            misL |= fp & ~U[B]
            wR &= B
            wL &= fp & ~U[B]
        sig = [0] * p
        for i in range(q):
            a = A[i]
            for x in range(p):
                if (a >> x) & 1:
                    sig[x] |= 1 << i
        vals = []
        for x in range(p):
            F = fq & ~sig[x]
            subs = []
            B = F
            while True:
                subs.append((bc[B], U[B]))
                if B == 0:
                    break
                B = (B - 1) & F
            for i in range(q):
                if (sig[x] >> i) & 1:
                    Ai = A[i]
                    b = -1
                    for c, u in subs:
                        v = c + p - pc[Ai | u]
                        if v > b:
                            b = v
                    vals.append((b, x, i))
        top = max(v[0] for v in vals)
        arg = [(x, i) for v, x, i in vals if v == top]
        ok = all(((misL >> x) & 1) or ((misR >> i) & 1) for x, i in arg)
        if not ok:
            stats[2] += 1
        if not all(((misL >> x) & 1) or ((misR >> i) & 1)
                   for v, x, i in vals):
            stats[6] += 1
        if len(arg) > 1:
            stats[7] += 1
        if controls:
            lo = min(v[0] for v in vals)
            amin = [(x, i) for v, x, i in vals if v == lo]
            if not all(((misL >> x) & 1) or ((misR >> i) & 1) for x, i in amin):
                stats[3] += 1
            if not all(((wL >> x) & 1) or ((wR >> i) & 1) for x, i in arg):
                stats[4] += 1
        if keys is not None:
            k = tuple(sorted(sig))
            d = (alpha, top, len(arg), ok, pc[misL], bc[misR], len(vals))
            if k in keys:
                if keys[k] != d:
                    stats[5] += 1
            else:
                keys[k] = d

    def rec(d, start, U):
        if d == q:
            leaf(U)
            return
        for a in range(start, 1 << p):
            A[d] = a
            rec(d + 1, a, U + [u | a for u in U])

    rec(0, 1, [0])
    return stats, keys


def primal_count_cell(p, q):
    """Count the encodings of shape (p,q) without evaluating anything on them.
    Every one of the first q-1 entries of the tuple is enumerated explicitly;
    the admissible last entries are counted exactly (they are the subsets that
    contain the still-uncovered part of [p] and are >= the previous entry), so
    the result is independent of the closed forms E and K."""
    from bisect import bisect_left
    fp = (1 << p) - 1
    sups = []
    for miss in range(1 << p):
        sups.append([m for m in range(1 << p) if m & miss == miss])
    tot = [0, 0]

    def rec(d, start, cover):
        if d == q - 1:
            tot[0] += (1 << p) - start
            s = sups[fp & ~cover]
            tot[1] += len(s) - bisect_left(s, start)
            return
        for a in range(start, 1 << p):
            rec(d + 1, a, cover | a)

    rec(0, 1, 0)
    return tot[0], tot[1]


def check_enumeration_counts():
    """Direct enumeration of every cell, compared with Table 1 row by row."""
    bad_cell, bad_row = [], []
    for n in range(2, PAPER_NMAX + 1):
        g = r = 0
        for p, q in cells(n):
            cg, cr = primal_count_cell(p, q)
            if (cg, cr) != (E_closed(p, q), K_closed(p, q)):
                bad_cell.append((p, q, cg, cr))
            g += cg
            r += cr
        if (g, r) != PAPER_TABLE[n]:
            bad_row.append((n, g, r, PAPER_TABLE[n]))
    check("prefix_enumeration_reproduces_E_and_K_every_cell", not bad_cell,
          "mismatches=%s" % bad_cell)
    check("prefix_enumeration_reproduces_table1_rows_n2_to_n12", not bad_row,
          "mismatches=%s" % bad_row)


def build_graph(A, p, q):
    """Explicit graph of the encoding: left vertices 0..p-1, right vertices
    p..p+q-1, edge x--(p+i) iff x in A_i.  Returns (n, adj, edges)."""
    n = p + q
    adj = [0] * n
    edges = []
    for i in range(q):
        for x in range(p):
            if (A[i] >> x) & 1:
                adj[x] |= 1 << (p + i)
                adj[p + i] |= 1 << x
                edges.append((x, p + i))
    return n, adj, edges


def independent_sets(adj, n):
    """All independent sets of G, sorted by decreasing size (brute force over
    all 2^n vertex subsets, straight from the definition)."""
    out = []
    for S in range(1 << n):
        T = S
        ok = True
        while T:
            b = T & -T
            u = b.bit_length() - 1
            if adj[u] & S:
                ok = False
                break
            T ^= b
        if ok:
            out.append(S)
    out.sort(key=lambda S: -bin(S).count("1"))
    return out


def alpha_brute(ind, allowed):
    """max |independent set| inside the induced subgraph on `allowed`."""
    for S in ind:
        if S & ~allowed == 0:
            return bin(S).count("1")
    return 0


def mis_support_brute(ind, n):
    """Union of all maximum independent sets of G."""
    a = bin(ind[0]).count("1")
    m = 0
    for S in ind:
        if bin(S).count("1") != a:
            break
        m |= S
    return a, m


def is_bipartite_bfs(adj, n):
    """2-colourability by BFS, from the definition."""
    col = [-1] * n
    for s in range(n):
        if col[s] >= 0:
            continue
        col[s] = 0
        stack = [s]
        while stack:
            u = stack.pop()
            m = adj[u]
            while m:
                b = m & -m
                v = b.bit_length() - 1
                m ^= b
                if col[v] < 0:
                    col[v] = 1 - col[u]
                    stack.append(v)
                elif col[v] == col[u]:
                    return False
    return True


def matching_size(A, p, q):
    """Maximum matching by augmenting paths; Koenig gives tau = |matching|."""
    matchR = [-1] * q
    total = 0
    for x in range(p):
        seen = [False] * q

        def try_aug(v):
            for i in range(q):
                if ((A[i] >> v) & 1) and not seen[i]:
                    seen[i] = True
                    if matchR[i] == -1 or try_aug(matchR[i]):
                        matchR[i] = v
                        return True
            return False

        if try_aug(x):
            total += 1
    return total


def tau_of_G_minus_Ne(A, p, q, x, i):
    """tau(G - N[e]) for e = x--i, computed by Koenig maximum matching on the
    induced bipartite subgraph.  This route never mentions alpha, so using it
    inside equation (2) turns Lemma 3's first identity into a real cross-check
    instead of the arithmetic tautology
        deg u + deg v - 2 + (|V| - |N[e]|) - alpha == n - 2 - alpha,
    which holds for every graph as soon as |N[e]| = deg u + deg v."""
    keep = [y for y in range(p) if not ((A[i] >> y) & 1)]
    AA = []
    for j in range(q):
        if (A[j] >> x) & 1:          # right vertex j lies in N(x), so deleted
            continue
        m = 0
        for t, y in enumerate(keep):
            if (A[j] >> y) & 1:
                m |= 1 << t
        AA.append(m)
    return matching_size(AA, len(keep), len(AA))


def iter_encodings(p, q, retained_only=True):
    """Every encoding of shape (p,q), as a list of q masks."""
    fp = (1 << p) - 1
    A = [0] * q

    def rec(d, start):
        if d == q:
            u = 0
            for a in A:
                u |= a
            if (u == fp) or not retained_only:
                yield list(A), u == fp
            return
        for a in range(start, 1 << p):
            A[d] = a
            for z in rec(d + 1, a):
                yield z

    return rec(0, 1)


def graph_report(A, p, q):
    """Everything the paper asserts about one encoded graph, computed twice:
    by brute force over all 2^n vertex subsets (the definitions) and by the
    bit-mask formulas (4), (5) of Section 3."""
    n, adj, edges = build_graph(A, p, q)
    fp, fq, full = (1 << p) - 1, (1 << q) - 1, (1 << (p + q)) - 1
    pop = lambda z: bin(z).count("1")
    U = [0] * (1 << q)
    for B in range(1, 1 << q):
        b = B & -B
        U[B] = U[B ^ b] | A[b.bit_length() - 1]
    ind = independent_sets(adj, n)
    alpha_b, mis_b = mis_support_brute(ind, n)
    tau = n - alpha_b
    alpha_bit = max(pop(B) + p - pop(U[B]) for B in range(1 << q))
    ms = [B for B in range(1 << q) if pop(B) + p - pop(U[B]) == alpha_bit]
    misR = misL = 0
    for B in ms:
        misR |= B
        misL |= fp & ~U[B]
    rep = {}
    rep["bipartite"] = is_bipartite_bfs(adj, n)
    # the third conjunct is the content: the declared sides really are a
    # bipartition of the built graph (n == p+q alone cannot fail, n was set to
    # p+q by build_graph).
    rep["noiso"] = (all(adj[u] for u in range(n)) and n == p + q
                    and all(adj[x] & fp == 0 for x in range(p))
                    and all(adj[p + i] >> p == 0 for i in range(q)))
    # len(edges) alone against sum(pop(A_i)) cannot see a MISPLACED edge, so
    # also require: no repeats, symmetry, and the handshake identity.
    rep["nedges"] = (len(edges) == sum(pop(a) for a in A)
                     and len(set(edges)) == len(edges)
                     and sum(pop(adj[u]) for u in range(n)) == 2 * len(edges)
                     and all((adj[u] >> v) & 1 and (adj[v] >> u) & 1
                             for u, v in edges))
    rep["eq4"] = alpha_bit == alpha_b
    rep["mis"] = (misL | (misR << p)) == mis_b
    rep["koenig"] = matching_size(A, p, q) == tau
    hval = []
    for u in range(n):
        allowed = full & ~((1 << u) | adj[u])
        hval.append(pop(adj[u]) + pop(allowed) - alpha_brute(ind, allowed))
    rep["h_ge_tau_iff_mis"] = all(
        hval[u] >= tau and ((hval[u] == tau) == bool((mis_b >> u) & 1))
        for u in range(n))
    lines, Ne_ok, ce_ok, eq5_ok, tau_ok = [], True, True, True, True
    for (x, v) in edges:
        i = v - p
        Ne = (1 << x) | (1 << v) | adj[x] | adj[v]
        if pop(Ne) != pop(adj[x]) + pop(adj[v]):
            Ne_ok = False
        allowed = full & ~Ne
        a_b = alpha_brute(ind, allowed)
        F = fq & ~(adj[x] >> p)
        a_bit = max(pop(B) + p - pop(A[i] | U[B])
                    for B in range(1 << q) if B & ~F == 0)
        if a_bit != a_b:
            eq5_ok = False
        # tau(G-N[e]) by maximum matching -- an alpha-free route.
        tau_e = tau_of_G_minus_Ne(A, p, q, x, i)
        if tau_e != pop(allowed) - a_b:
            tau_ok = False
        # equation (2) evaluated on the matching route, then compared with
        # Lemma 3's closed form n-2-alpha(G-N[e]).
        c_e = pop(adj[x]) + pop(adj[v]) - 2 + tau_e
        if c_e != n - 2 - a_b:
            ce_ok = False
        lines.append((c_e, min(hval[x], hval[v]), a_b, x, v))
    rep["Ne_size"] = Ne_ok
    rep["eq5"] = eq5_ok
    rep["tau_Ne_matching"] = tau_ok
    rep["ce_identity"] = ce_ok
    rep["min_g_eq_tau"] = min(g for c, g, a, x, v in lines) == tau
    top = max(a for c, g, a, x, v in lines)
    arg = [(x, v) for c, g, a, x, v in lines if a == top]
    rep["prop5"] = all(((mis_b >> x) & 1) or ((mis_b >> v) & 1)
                       for x, v in arg)
    rep["cond3"] = any(((mis_b >> x) & 1) or ((mis_b >> v) & 1)
                       for x, v in arg)
    cmin = min(c for c, g, a, x, v in lines)
    rep["vJ"] = cmin == n - 2 - top
    env = [min(c + (t - 1) * g for c, g, a, x, v in lines)
           for t in range(1, 23)]
    rep["affine"] = all(env[t - 1] == cmin + (t - 1) * tau
                        for t in range(1, 21))
    # Proposition 4 asserts a THREE-way equivalence (1)<=>(2)<=>(3); the
    # program used to test only (2)<=>(3).  (1) is the weaker "the envelope is
    # SOME affine function of t", which is what Vu's question actually asks;
    # here it is computed separately from its differences.  Every breakpoint of
    # a lower envelope of lines with integer intercepts in [0,n] and pairwise
    # distinct integer slopes lies at t-1 <= n, so 22 points settle all t >= 1.
    rep["affine_is_some_line"] = (
        all(env[t] - env[t - 1] == env[1] - env[0]
            for t in range(1, len(env))) == rep["affine"])
    return rep


DEF_KEYS = [
    ("encoded_graphs_are_bipartite_2colourable", "bipartite"),
    ("encoded_graphs_have_no_isolated_vertex_and_n_eq_p_plus_q", "noiso"),
    ("edge_count_equals_sum_of_neighbourhood_sizes", "nedges"),
    ("eq4_alpha_bitset_equals_bruteforce_alpha", "eq4"),
    ("mis_support_rule_equals_bruteforce_union_of_all_mis", "mis"),
    ("koenig_maximum_matching_equals_tau", "koenig"),
    ("lemma3_card_Ne_equals_deg_u_plus_deg_v", "Ne_size"),
    ("eq5_alpha_of_G_minus_Ne_equals_bruteforce", "eq5"),
    ("tau_of_G_minus_Ne_by_matching_equals_size_minus_alpha", "tau_Ne_matching"),
    ("lemma3_ce_equals_n_minus_2_minus_alpha_G_minus_Ne", "ce_identity"),
    ("lemma3_h_ge_tau_with_equality_iff_vertex_in_some_mis",
     "h_ge_tau_iff_mis"),
    ("lemma3_min_over_edges_of_g_e_equals_tau", "min_g_eq_tau"),
    ("prop5_every_argmax_edge_meets_a_mis_bruteforce", "prop5"),
    ("vJ_at_t1_equals_n_minus_2_minus_max_alpha", "vJ"),
    ("theorem_envelope_affine_for_t_1_to_20", "affine"),
    ("prop4_1_iff_2_envelope_is_affine_iff_its_slope_is_tau",
     "affine_is_some_line"),
]


def check_definitions(nmax):
    """Full brute-force pass over every retained encoding with n <= nmax."""
    fails = dict((k, 0) for _, k in DEF_KEYS)
    first = {}
    tot = 0
    equiv_bad = 0
    discard_bad = 0
    discard_seen = 0
    for n in range(2, nmax + 1):
        for p, q in cells(n):
            for A, retained in iter_encodings(p, q, retained_only=False):
                nn, adj, edges = build_graph(A, p, q)
                if not retained:
                    discard_seen += 1
                    # "all(adj[u])" alone cannot fail: discarded MEANS union
                    # != [p], and a left vertex outside the union has degree 0
                    # by construction.  Require the SHARP statement instead --
                    # the isolated vertices of the built graph are exactly the
                    # left vertices missed by the union, and there is at least
                    # one -- which does compare build_graph with the retention
                    # predicate rather than restating it.
                    u = 0
                    for a in A:
                        u |= a
                    iso = set(w for w in range(nn) if adj[w] == 0)
                    want = set(y for y in range(p) if not ((u >> y) & 1))
                    if iso != want or not iso:
                        discard_bad += 1
                    continue
                tot += 1
                rep = graph_report(A, p, q)
                for _, k in DEF_KEYS:
                    if not rep[k]:
                        fails[k] += 1
                        first.setdefault(k, (n, list(A)))
                if rep["affine"] != rep["cond3"]:
                    equiv_bad += 1
    for name, k in DEF_KEYS:
        check(name + "_n_le_%d" % nmax, fails[k] == 0,
              "%d encodings; failures=%d%s"
              % (tot, fails[k],
                 "" if not fails[k] else " first=%s" % (first[k],)))
    check("prop4_equivalence_affine_iff_some_argmax_edge_meets_mis_n_le_%d"
          % nmax, equiv_bad == 0,
          "%d encodings; mismatches=%d" % (tot, equiv_bad))
    check("discarded_encodings_all_have_an_isolated_vertex",
          discard_bad == 0 and discard_seen > 0,
          "%d discarded encodings; the isolated-vertex set of each equals the "
          "left vertices outside the union, and is nonempty; "
          "counterexamples=%d" % (discard_seen, discard_bad))
    return tot


def check_discards_lose_nothing(nmax=9):
    """The previous check cannot fail: "discarded" MEANS union != [p], and a
    left vertex outside the union is exactly an isolated vertex, so it merely
    unfolds a definition.  The statement with content is that discarding loses
    no GRAPH: delete the isolated vertices of a discarded encoding and the
    result must itself be a retained encoding of a smaller census cell (with
    the sides swapped when the surviving left side is the smaller one)."""
    RET = {}
    for n in range(2, nmax + 1):
        for p, q in cells(n):
            RET[(p, q)] = set(tuple(A) for A, _ in iter_encodings(p, q))
    missing, seen = [], 0
    for n in range(2, nmax + 1):
        for p, q in cells(n):
            for A, retained in iter_encodings(p, q, retained_only=False):
                if retained:
                    continue
                seen += 1
                u = 0
                for a in A:
                    u |= a
                keep = [y for y in range(p) if (u >> y) & 1]
                pp = len(keep)
                AA = [sum(1 << t for t, y in enumerate(keep) if (a >> y) & 1)
                      for a in A]
                if pp >= q:
                    cell, key = (pp, q), tuple(sorted(AA))
                else:
                    BB = [sum(1 << i for i in range(q) if (AA[i] >> t) & 1)
                          for t in range(pp)]
                    cell, key = (q, pp), tuple(sorted(BB))
                if key not in RET.get(cell, ()):
                    if len(missing) < 4:
                        missing.append((p, q, tuple(A), cell, key))
    check("discarded_encodings_reappear_in_the_census_once_isolated_vertices_"
          "are_deleted", not missing and seen > 0,
          "%d discarded encodings relocated; not found=%s" % (seen, missing))


_DTAB = {}


def dual_tables(q):
    """bc = |B|; SUBS[C] = all subsets of C; UPD[S] = the pairs used to update
    Z and vb when a left vertex of signature S is appended; SUP[M] = all
    supersets of M."""
    if q in _DTAB:
        return _DTAB[q]
    fq = (1 << q) - 1
    bc = [bin(i).count("1") for i in range(1 << q)]
    SUBS = []
    for C in range(1 << q):
        s, T = [], C
        while True:
            s.append(T)
            if T == 0:
                break
            T = (T - 1) & C
        SUBS.append(tuple(s))
    UPD = [tuple((S | T, T) for T in SUBS[fq & ~S]) for S in range(1 << q)]
    SUP = [tuple(M for M in range(1 << q) if M & Mi == Mi)
           for Mi in range(1 << q)]
    _DTAB[q] = (fq, bc, SUBS, UPD, SUP)
    return _DTAB[q]


def dual_cell(p, q, mode="fast", keydata=None, corrupt=False):
    """Transposed scan of the cell (p,q): one representative per orbit of the
    retained encodings under relabelling of the p-side, obtained by enumerating
    the multiset of the p left signatures S_x subset of [q] (nondecreasing,
    nonempty, covering [q]).  Z[M] = #{x : S_x subset of M} and
    vb[B] = |B| + Z[complement B] = the size of the independent set B u
    (L \\ U(B)), both maintained incrementally down the recursion, so
    alpha(G) = max vb by equation (4).

    stats = [representatives, short-circuited leaves, fully computed leaves,
             violations of Proposition 5, fast/full verdict disagreements]
    """
    fq, bc, SUBS, UPD, SUP = dual_tables(q)
    Z = bytearray(1 << q)
    vb = bytearray(bc)
    stack = [0] * p
    stats = [0, 0, 0, 0, 0]

    def edgeval(S, i):
        """alpha(G - N[e]) for e = x i with signature S_x = S, by (5)."""
        bi = 1 << i
        best = -1
        for B in SUBS[fq & ~S]:
            v = bc[B] + Z[fq ^ (B | bi)]
            if v > best:
                best = v
        return best

    def classify(ms, badR):
        """Split the edge classes (S,i) into those with an endpoint in a
        maximum independent set and those without."""
        bad, good = [], []
        for S in set(stack):
            if corrupt:
                # deliberate misreading: "in EVERY maximiser" (self-test only)
                okL = all(B & S == 0 for B in ms)
            else:
                okL = False
                for B in ms:
                    if B & S == 0:
                        okL = True
                        break
            t = S
            while t:
                b = t & -t
                t ^= b
                i = b.bit_length() - 1
                if okL or not ((badR >> i) & 1):
                    good.append((S, i))
                else:
                    bad.append((S, i))
        return bad, good

    def maximisers():
        alpha = max(vb)
        ms, k = [], 0
        while True:
            pos = vb.find(alpha, k)
            if pos < 0:
                break
            ms.append(pos)
            k = pos + 1
        if corrupt:
            misR = fq
            for B in ms:
                misR &= B
            return alpha, ms, misR
        misR = 0
        for B in ms:
            misR |= B
        return alpha, ms, misR

    def leaf_fast():
        stats[0] += 1
        alpha, ms, misR = maximisers()
        badR = fq & ~misR
        if badR == 0:
            stats[1] += 1
            return
        bad, good = classify(ms, badR)
        if not bad:
            stats[1] += 1
            return
        stats[2] += 1
        vbad = -1
        for (S, i) in bad:
            v = edgeval(S, i)
            if v > vbad:
                vbad = v
        for (S, i) in good:
            if edgeval(S, i) > vbad:
                return
        stats[3] += 1

    def leaf_full():
        """Same leaf without any short circuit: every edge value is computed,
        the argmax set is formed, and each of its edges is tested."""
        stats[0] += 1
        alpha, ms, misR = maximisers()
        badR = fq & ~misR
        bad, good = classify(ms, badR)
        mult = {}
        for S in stack:
            mult[S] = mult.get(S, 0) + 1
        vals = [(edgeval(S, i), S, i, True) for (S, i) in good]
        vals += [(edgeval(S, i), S, i, False) for (S, i) in bad]
        top = max(v[0] for v in vals)
        argmax = [(S, i, g) for v, S, i, g in vals if v == top]
        ok_full = all(g for S, i, g in argmax)
        vbad = max([v for v, S, i, g in vals if not g] + [-1])
        vgood = max([v for v, S, i, g in vals if g] + [-1])
        ok_fast = (badR == 0) or (not bad) or (vgood > vbad)
        if ok_full != ok_fast:
            stats[4] += 1
        if not ok_full:
            stats[3] += 1
        if badR == 0 or not bad:
            stats[1] += 1
        else:
            stats[2] += 1
        if keydata is not None:
            nleft = sum(c for S, c in mult.items()
                        if any(B & S == 0 for B in ms))
            keydata[tuple(stack)] = (
                alpha, top, sum(mult[S] for S, i, g in argmax),
                ok_full, nleft, bc[misR],
                sum(bc[S] for S in stack))

    leaf = leaf_fast if mode == "fast" else leaf_full

    def rec(d, start, cover):
        if d == p - 1:
            for S in SUP[fq & ~cover]:
                if S < start:
                    continue
                stack[d] = S
                for m, t in UPD[S]:
                    Z[m] += 1
                    vb[t] += 1
                leaf()
                for m, t in UPD[S]:
                    Z[m] -= 1
                    vb[t] -= 1
            return
        for S in range(start, 1 << q):
            stack[d] = S
            for m, t in UPD[S]:
                Z[m] += 1
                vb[t] += 1
            rec(d + 1, S, cover | S)
            for m, t in UPD[S]:
                Z[m] -= 1
                vb[t] -= 1

    rec(0, 1, 0)
    return stats


_PERM = {}


def perm_tables(q):
    """For every permutation of the q-side, the induced map on masks."""
    if q not in _PERM:
        tabs = []
        for sigma in permutations(range(q)):
            t = [0] * (1 << q)
            for m in range(1 << q):
                r = 0
                for i in range(q):
                    if (m >> i) & 1:
                        r |= 1 << sigma[i]
                t[m] = r
            tabs.append(t)
        _PERM[q] = tabs
    return _PERM[q]


def canon_key(key, q):
    """Lexicographically least relabelling of the q-side.  Identifies the
    left-orbit invariant of the literal scan with that of the transposed scan.
    """
    best = None
    for t in perm_tables(q):
        c = tuple(sorted(t[s] for s in key))
        if best is None or c < best:
            best = c
    return best


def labeled_matrix_count(p, q):
    """The number of p x q incidence matrices with no zero row and no zero
    column, i.e. of bipartite graphs on labelled sides [p] and [q] with no
    isolated vertex.  Independent of E and K."""
    return sum((-1) ** r * nCk(p, r) * (2 ** (p - r) - 1) ** q
               for r in range(p + 1))


def orderings(seq, m):
    """m! / prod(multiplicity!) -- the size of the orbit of the multiset."""
    from math import factorial
    cnt = {}
    for s in seq:
        cnt[s] = cnt.get(s, 0) + 1
    r = factorial(m)
    for c in cnt.values():
        r //= factorial(c)
    return r


def selftest_transposed_scanner(nmax=8):
    """Mutation test of the scanner that carries the census: with the single
    change "in EVERY maximiser" for "in SOME maximiser" it must report
    violations, otherwise its violation counter is dead code."""
    viol = reps = 0
    for n in range(2, nmax + 1):
        for p, q in cells(n):
            st = dual_cell(p, q, "fast", corrupt=True)
            reps += st[0]
            viol += st[3]
    check("selftest_transposed_scanner_reports_violations_when_corrupted",
          viol > 0, "%d of %d representatives with n <= %d are reported as "
          "violations once the maximum-independent-set rule is corrupted"
          % (viol, reps, nmax))


def check_leaf_paths_agree(nmax=8):
    """The census up to n = 12 is carried by leaf_fast, whose short circuit
    ("stop as soon as some edge meeting a maximum independent set beats every
    edge that does not") decides 104 million of the 110 million
    representatives.  leaf_full recomputes the same verdict with no short
    circuit at all.  On the TRUE census both verdicts are always "no
    violation", so agreeing there says nothing about the short circuit: any
    error that only ever turns a violation into a pass -- for instance the
    single character `>` weakened to `>=`, which is exactly the case where a
    non-meeting edge TIES the maximum and so genuinely violates Proposition 5
    -- would be invisible.  So the two paths are also compared on the corrupted
    family, where the verdict is false for most representatives.  AUX_SHAPES is
    included so that the comparison also runs at q = 5 and q = 6."""
    shapes = [pq for n in range(2, nmax + 1) for pq in cells(n)] + AUX_SHAPES
    for corrupt in (False, True):
        f = [0] * 5
        u = [0] * 5
        for p, q in shapes:
            a = dual_cell(p, q, "fast", corrupt=corrupt)
            b = dual_cell(p, q, "full", {}, corrupt=corrupt)
            for j in range(4):
                f[j] += a[j]
                u[j] += b[j]
            u[4] += b[4]
        if not corrupt:
            check("census_scanner_fast_and_full_agree_on_the_true_family_n_le_%d"
                  % nmax, f[:4] == u[:4] and u[4] == 0 and f[3] == 0,
                  "fast=%s full=%s over %d shapes including the widths %s "
                  "(violations must be 0 on both)"
                  % (f[:4], u[:4], len(shapes),
                     sorted(set(q for p, q in shapes))))
        else:
            check("census_scanner_short_circuit_agrees_with_full_where_the_"
                  "verdict_is_FALSE_n_le_%d" % nmax,
                  f[:4] == u[:4] and u[4] == 0 and f[3] > 0,
                  "%d of %d representatives violate under the corrupted "
                  "maximum-independent-set rule, and the short-circuited and "
                  "fully computed leaves report the same %d" % (f[3], f[0], u[3]))


AUX_SHAPES = [(4, 5), (3, 6)]


def check_two_routes_agree(max_retained=150000):
    """The literal encoding scan and the transposed orbit scan must see the
    same graphs and compute the same invariants on them.

    The affordable census cells all have q <= 4, yet the transposed scan that
    carries n = 11 and n = 12 runs at q = 5 and q = 6 (those two widths cover
    98% of the census).  AUX_SHAPES therefore adds two off-census shapes whose
    BOTH routes are cheap, purely so that the dual_tables/edgeval code is
    cross-checked at the widths it is actually trusted with."""
    bad_set, bad_data, bad_inv, cells_done = [], [], [], []
    bad_orbit = []
    todo = []
    for n in range(2, PAPER_NMAX + 1):
        for p, q in cells(n):
            if K_closed(p, q) <= max_retained:
                todo.append((n, p, q))
    todo += [(None, p, q) for p, q in AUX_SHAPES]
    if True:
        for n, p, q in todo:
            cells_done.append((n, p, q))
            ps, pk = primal_cell(p, q, want_keys=True)
            kd = {}
            st = dual_cell(p, q, "full", kd)
            if ps[5] or st[4]:
                bad_inv.append((p, q, ps[5], st[4]))
            xp = sum(orderings(A, q) for A, _ in iter_encodings(p, q))
            xd = sum(orderings(k, p) for k in kd)
            xc = labeled_matrix_count(p, q)
            if xp != xc or xd != xc:
                bad_orbit.append((p, q, xp, xd, xc))
            cp, cd = {}, {}
            for src, dst in ((pk, cp), (kd, cd)):
                for k, d in src.items():
                    ck = canon_key(k, q)
                    if ck in dst and dst[ck] != d:
                        bad_data.append((p, q, ck, dst[ck], d))
                    dst[ck] = d
            if set(cp) != set(cd):
                bad_set.append((p, q, len(cp), len(cd)))
            else:
                for ck in cp:
                    if cp[ck] != cd[ck]:
                        bad_data.append((p, q, ck, cp[ck], cd[ck]))
    check("literal_and_transposed_scans_see_the_same_graphs", not bad_set,
          "%d cells; mismatched cells=%s" % (len(cells_done), bad_set[:4]))
    check("literal_and_transposed_scans_compute_the_same_invariants",
          not bad_data, "differences=%s" % (bad_data[:4],))
    check("invariants_are_stable_under_relabelling_both_sides", not bad_inv,
          "cells with an unstable value=%s" % (bad_inv[:4],))
    check("both_enumerations_weighted_by_orbit_size_count_the_same_graphs",
          not bad_orbit, "mismatches=%s" % (bad_orbit[:4],))
    check("two_route_comparison_covered_an_n12_cell",
          any(c[0] == 12 for c in cells_done),
          "n=12 cells compared=%s of %d cells in all"
          % ([c for c in cells_done if c[0] == 12], len(cells_done)))
    widths = sorted(set(c[2] for c in cells_done))
    check("two_route_comparison_covered_the_widths_q5_and_q6_used_at_n11_n12",
          5 in widths and 6 in widths,
          "widths q compared=%s; the census at n = 11 and n = 12 uses q up to "
          "%d, and 98%% of its representatives sit at q = 5 and q = 6"
          % (widths, max(q for p, q in cells(PAPER_NMAX))))


def check_literal_census(nmax):
    """Proposition 5 on every single retained encoding with n <= nmax, by the
    literal route of equations (4) and (5), together with the two deliberate
    corruptions of the test."""
    tot = [0] * 8
    rows_bad = []
    for n in range(2, nmax + 1):
        g = r = 0
        for p, q in cells(n):
            st, _ = primal_cell(p, q, controls=True)
            for j in range(8):
                tot[j] += st[j]
            g += st[0]
            r += st[1]
        if (g, r) != PAPER_TABLE[n]:
            rows_bad.append((n, g, r, PAPER_TABLE[n]))
    check("literal_scan_reproduces_table1_rows_n_le_%d" % nmax, not rows_bad,
          "mismatches=%s" % rows_bad)
    check("prop5_zero_violations_literal_scan_n_le_%d" % nmax, tot[2] == 0,
          "%d encodings scanned, %d violations" % (tot[1], tot[2]))
    check("control_argmin_instead_of_argmax_does_violate", tot[3] > 0,
          "%d encodings fail Proposition 5 when M(G) is taken to be the "
          "argMIN of alpha(G-N[e])" % tot[3])
    check("control_in_every_maximiser_rule_does_violate", tot[4] > 0,
          "%d encodings fail when membership in a maximum independent set is "
          "read off EVERY maximiser instead of SOME maximiser" % tot[4])
    check("prop5_is_not_a_statement_about_all_edges", tot[6] > 0,
          "%d of %d encodings own an edge with neither endpoint in any maximum "
          "independent set, so the conclusion constrains WHICH edges maximise "
          "alpha(G-N[e])" % (tot[6], tot[1]))
    check("argmax_set_is_often_larger_than_one_edge", tot[7] > 0,
          "%d of %d encodings have |M(G)| >= 2, so the universal quantifier in "
          "Proposition 5 is stronger than the existential one needed by "
          "Proposition 4(3)" % (tot[7], tot[1]))
    return tot


def check_full_census():
    """Proposition 5 over the whole census up to n = 12, by the transposed
    scan: every bipartite graph with at most twelve nonisolated vertices is
    represented by at least one of these orbit representatives."""
    tot = [0] * 5
    per_n = {}
    bad_count = []
    for n in range(2, PAPER_NMAX + 1):
        s = [0] * 5
        for p, q in cells(n):
            st = dual_cell(p, q, "fast")
            if st[0] != K_closed(q, p):
                bad_count.append((p, q, st[0], K_closed(q, p)))
            for j in range(5):
                s[j] += st[j]
        per_n[n] = s[0]
        for j in range(5):
            tot[j] += s[j]
    check("transposed_scan_representative_count_matches_closed_form",
          not bad_count, "%d representatives over %d cells; mismatches=%s"
          % (tot[0], sum(len(cells(n)) for n in range(2, PAPER_NMAX + 1)),
             bad_count))
    check("PROP5_ZERO_VIOLATIONS_WHOLE_CENSUS_N_LE_12", tot[3] == 0,
          "%d orbit representatives covering all %d retained encodings; "
          "violations=%d" % (tot[0], PAPER_TOTAL_RETAINED, tot[3]))
    check("both_leaf_paths_exercised_in_the_census",
          tot[1] > 0 and tot[2] > 0,
          "short-circuited=%d fully-computed=%d" % (tot[1], tot[2]))
    note("representatives by n: %s"
         % ", ".join("n=%d:%d" % (k, v) for k, v in sorted(per_n.items())))
    return tot, per_n


def check_general_graphs(nmax=6):
    """The same test applied to EVERY graph without isolated vertices on at
    most nmax vertices, bipartite or not.  This is where the condition of
    Proposition 5 is false, so the test is not vacuous; and it is where
    Proposition 4's equivalence has both sides false, so that is not vacuous
    either.  Nothing here is taken from the paper except the two statements."""
    viol = c3false = equiv_bad = 0
    bip_seen = bip_viol = 0
    lin_bad = lin_seen = 0
    total = 0
    ex = None
    for n in range(3, nmax + 1):
        pairs = list(combinations(range(n), 2))
        full = (1 << n) - 1
        for mask in range(1 << len(pairs)):
            adj = [0] * n
            edges = []
            for k, (u, v) in enumerate(pairs):
                if (mask >> k) & 1:
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
                    edges.append((u, v))
            if not edges or any(adj[u] == 0 for u in range(n)):
                continue
            total += 1
            ind = independent_sets(adj, n)
            alpha, mis = mis_support_brute(ind, n)
            tau = n - alpha
            lines = []
            for (u, v) in edges:
                allowed = full & ~((1 << u) | (1 << v) | adj[u] | adj[v])
                a = alpha_brute(ind, allowed)
                h = []
                for w in (u, v):
                    al = full & ~((1 << w) | adj[w])
                    h.append(bin(adj[w]).count("1") + bin(al).count("1")
                             - alpha_brute(ind, al))
                lines.append((n - 2 - a, min(h), a, u, v))
            top = max(x[2] for x in lines)
            arg = [(u, v) for c, g, a, u, v in lines if a == top]
            every = all(((mis >> u) & 1) or ((mis >> v) & 1) for u, v in arg)
            some = any(((mis >> u) & 1) or ((mis >> v) & 1) for u, v in arg)
            cmin = min(x[0] for x in lines)
            env = [min(c + (t - 1) * g for c, g, a, u, v in lines)
                   for t in range(1, 27)]
            affine = all(env[t - 1] == cmin + (t - 1) * tau
                         for t in range(1, 26))
            linear = all(env[t] - env[t - 1] == env[1] - env[0]
                         for t in range(1, len(env)))
            if linear != affine:
                lin_bad += 1
            if not affine:
                lin_seen += 1
            if not every:
                viol += 1
                if ex is None:
                    ex = (n, edges)
            if not some:
                c3false += 1
            if affine != some:
                equiv_bad += 1
            if is_bipartite_bfs(adj, n):
                bip_seen += 1
                if not every:
                    bip_viol += 1
    check("control_general_graphs_do_violate_prop5", viol > 0,
          "%d of %d graphs on <=%d vertices violate it; smallest witness=%s"
          % (viol, total, nmax, ex))
    check("control_general_graphs_with_no_good_argmax_edge_exist",
          c3false > 0, "%d graphs have NO argmax edge meeting a maximum "
          "independent set" % c3false)
    check("prop4_equivalence_holds_on_all_graphs_n_le_%d" % nmax,
          equiv_bad == 0,
          "%d graphs, %d mismatches between affineness and criterion (3)"
          % (total, equiv_bad))
    check("bipartite_graphs_in_that_family_never_violate_prop5",
          bip_viol == 0 and bip_seen > 0,
          "%d bipartite graphs, %d violations" % (bip_seen, bip_viol))
    check("prop4_1_iff_2_on_all_graphs_n_le_%d" % nmax,
          lin_bad == 0 and lin_seen > 0,
          "%d graphs, %d mismatches between 'the envelope is SOME affine "
          "function' (Proposition 4(1), i.e. Vu's question) and 'it is "
          "v(J)+(t-1)tau' (Proposition 4(2)); the two are compared on a family "
          "in which %d graphs make both false" % (total, lin_bad, lin_seen))


def check_representation(nmax=7):
    """Completeness of the encoding.  Every LABELLED bipartite graph without
    isolated vertices on n <= nmax vertices is generated here from scratch (a
    2-colouring S and an arbitrary edge set between S and its complement),
    turned into the tuple of Section 3 with the smaller side as the q-side, and
    looked up in the census.  Conversely every retained encoding must be hit,
    which pins down the labelling convention.

    NOTE ON THE HEADLINE NUMBER: the loop below runs over (bipartition, edge
    set) PAIRS, and a disconnected bipartite graph has several valid
    bipartitions, so the pair count (169,656 for n <= 7) is about 2.2x the
    number of DISTINCT labelled graphs (77,340).  Both are reported; the
    completeness statement is about the distinct graphs."""
    RET, HIT = {}, {}
    for n in range(2, nmax + 1):
        for p, q in cells(n):
            RET[(p, q)] = set(tuple(A) for A, _ in iter_encodings(p, q))
            HIT[(p, q)] = set()
    missing = []
    total = 0
    distinct = set()
    for n in range(2, nmax + 1):
        for k in range(1, n):
            for S in combinations(range(n), k):
                Sl = list(S)
                Tl = [v for v in range(n) if v not in S]
                m = n - k
                for em in range(1 << (k * m)):
                    adj = [0] * n
                    for a in range(k):
                        for b in range(m):
                            if (em >> (a * m + b)) & 1:
                                adj[Sl[a]] |= 1 << Tl[b]
                                adj[Tl[b]] |= 1 << Sl[a]
                    if any(adj[u] == 0 for u in range(n)):
                        continue
                    total += 1
                    distinct.add((n, tuple(adj)))
                    R, L = (Sl, Tl) if k <= m else (Tl, Sl)
                    p, q = len(L), len(R)
                    pos = dict((v, j) for j, v in enumerate(L))
                    A = []
                    for v in R:
                        w = 0
                        for u in pos:
                            if (adj[v] >> u) & 1:
                                w |= 1 << pos[u]
                        A.append(w)
                    key = tuple(sorted(A))
                    if key not in RET.get((p, q), ()):
                        if len(missing) < 4:
                            missing.append((n, (p, q), key))
                    else:
                        HIT[(p, q)].add(key)
    unhit = [(pq, len(RET[pq]) - len(HIT[pq]))
             for pq in RET if HIT[pq] != RET[pq]]
    check("every_labelled_bipartite_graph_n_le_%d_is_in_the_census" % nmax,
          not missing, "%d DISTINCT labelled bipartite graphs without isolated "
          "vertices, reached through %d (bipartition, edge set) pairs; not "
          "found=%s" % (len(distinct), total, missing))
    check("every_retained_encoding_n_le_%d_comes_from_such_a_graph" % nmax,
          not unhit, "cells with unreached encodings=%s" % unhit)
    return len(distinct), total


# =====================================================================
# MAIN
# =====================================================================

LITERAL_NMAX = 10
BRUTE_NMAX = 9
# Completeness of the Section 3 encoding is machine-checked only up to here:
# the outer loop is 2^(k*(n-k)) edge sets per bipartition, so n = 8 is already
# out of reach.  The shortfall is named in the closing NOT RE-RUN note.
REPR_NMAX = 7


def main():
    print("verify.py -- Vu's cover-ideal linearity question through twelve "
          "nonisolated vertices")
    note("exact integer and bit-mask arithmetic only; no floating point, no "
         "randomness, no external data")
    note("INPUTS from the paper: Table 1 (eleven rows and two totals), the "
         "encoding of Section 3, the closed forms E and K, equations (4)-(5), "
         "and the statements of Propositions 4 and 5")
    check_table_from_closed_forms()
    check_counts_by_dp()
    check_enumeration_counts()
    check_definitions(BRUTE_NMAX)
    check_discards_lose_nothing(BRUTE_NMAX)
    repr_distinct, repr_pairs = check_representation(REPR_NMAX)
    check_general_graphs(6)
    selftest_transposed_scanner(8)
    check_leaf_paths_agree(8)
    check_two_routes_agree()
    lit = check_literal_census(LITERAL_NMAX)
    tot, per_n = check_full_census()
    note("SCOPE: what is re-derived for EVERY bipartite graph with at most 12 "
         "nonisolated vertices is Proposition 5 -- every edge maximising "
         "alpha(G-N[e]) has an endpoint in a maximum independent set. Theorem "
         "1's own identity v(J^t) = v(J) + (t-1)tau is exhibited by direct "
         "computation only for n <= %d (t = 1..20); for n = 10, 11, 12 it is "
         "not computed here, but follows from Proposition 5 through "
         "Proposition 4, which is proved analytically in the paper and "
         "machine-checked here only for n <= %d, and through Vu's cited local "
         "formula (1)-(3)." % (BRUTE_NMAX, BRUTE_NMAX))
    note("The n <= %d layer is scanned encoding by encoding: all %d retained "
         "encodings of Table 1, one at a time, by the literal route of "
         "equations (4)-(5)." % (LITERAL_NMAX, lit[1]))
    note("NOT RE-RUN one encoding at a time: the %d retained encodings with "
         "n = 11 and n = 12. They are covered by the %d orbit representatives "
         "the transposed scan has at those two values of n (%d of them at "
         "n = 12; the whole census up to n = 12 takes %d representatives): "
         "every retained "
         "encoding of a cell (p,q) is a relabelling of the p-side of exactly "
         "one representative, and alpha(G), the union of the maximum "
         "independent sets and M(G) are isomorphism invariants, so no graph is "
         "left untested. The equality of the two routes is checked cell by "
         "cell wherever both are affordable, including two n = 12 cells."
         % (PAPER_TABLE[11][1] + PAPER_TABLE[12][1],
            per_n[11] + per_n[12], per_n[12], tot[0]))
    note("Table 1's own counts, including the n = 12 row (%d generated, %d "
         "retained), are nonetheless recomputed three independent ways: the "
         "paper's closed forms, an unbounded-knapsack dynamic program, and "
         "prefix enumeration." % PAPER_TABLE[12])
    note("NOT RE-DERIVED (cited algebra, not computation): Lemma 2 on "
         "isolated variables, the standard gradedness J^t = J^(t) of "
         "Herzog-Hibi-Trung, and Vu's local formula (1)-(3). Everything "
         "downstream of them is recomputed here, and the local formula is "
         "used only in the affineness check.")
    note("NOT RE-RUN: the program that produced Table 1, which the paper "
         "keeps in an auxiliary archive available on request and does not "
         "reproduce -- this program is an independent reimplementation in "
         "Python, written from the paper's printed data alone. NOT VERIFIED BY "
         "MACHINE ANYWHERE ABOVE: that the encoding of Section 3 reaches every "
         "labelled bipartite graph without isolated vertices on 8 to 12 "
         "nonisolated vertices. That completeness is checked by exhaustive "
         "construction only for n <= %d (%d distinct labelled graphs, reached "
         "through %d (bipartition, edge set) pairs); for 8 <= n <= 12 it is "
         "argued from the encoding, not computed."
         % (REPR_NMAX, repr_distinct, repr_pairs))
    return verdict()


if __name__ == "__main__":
    sys.exit(main())

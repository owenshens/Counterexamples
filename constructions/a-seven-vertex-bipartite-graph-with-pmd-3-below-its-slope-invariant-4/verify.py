#!/usr/bin/env python3
"""Exact verification of the note

    "A Seven-Vertex Bipartite Graph Whose Positive Matching Decomposition Number
     Lies Below Its Slope Invariant"

Python 3.9+, STANDARD LIBRARY ONLY (itertools, sys).  No third-party package, no
external data file: every object this program consumes is a literal below, copied
from the note, and every quantity the note states about those objects -- the
witness Gamma and its pmd and K -- is re-derived here.  Steps 7 and 8 compute a
census, a minimality/uniqueness statement and a second witness M_5 that the note
does not claim.  All arithmetic is over Python integers and bitmasks; no float is
computed anywhere, so no decision depends on a tolerance.

The bipartite graphs are carried exactly as the source's setup carries them: a pair
(m, n) of side sizes together with m row bitmasks over the n columns, bit j of row i
meaning the edge x_{i+1} y_{j+1}.  Indices inside the program are 0-based; the paper
is 1-based, and the translation is stated at each printed object.

One `PASS <name> [detail]` line per check; `VERDICT: ALL <n> CHECKS PASS` at the end,
exit 0 iff every check passed.
"""

import itertools
import sys

# ---------------------------------------------------------------------------
# THE OBJECTS, AS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------

# Section 2: the witness.  X = {x1,x2,x3}, Y = {y1,y2,y3,y4}.
G_ROWS = (3, 13, 14)                     # row bitmasks over y1..y4
G_M, G_N = 3, 4
G_G6 = 'FEhb?'                           # graph6 string printed in the paper
# the 8 edges as printed (1-based (i,j) meaning x_i y_j)
G_EDGES_1 = ((1, 1), (1, 2), (2, 1), (2, 3), (2, 4), (3, 2), (3, 3), (3, 4))
# the graph6 labelling printed in the paper: 0,1,2 = x1,x2,x3 and 3,4,5,6 = y1..y4
G_G6_EDGES = ((0, 3), (0, 4), (1, 3), (1, 5), (1, 6), (2, 4), (2, 5), (2, 6))

# The three-part decomposition printed in the paper (1-based).
PART1_1 = ((2, 1), (3, 3))
PART2_1 = ((1, 1), (2, 4), (3, 2))
PART3_1 = ((1, 2), (2, 3), (3, 4))
# The path Gamma - M_1, printed in the paper as a vertex sequence.
PATH_AFTER_M1 = ('y1', 'x1', 'y2', 'x3', 'y4', 'x2', 'y3')

# Section 1: the source's two worked examples and the parts its Algorithm 1 prints.
EX1_ROWS = (0b110, 0b101, 0b011)         # Example (i): C_6
EX1_PARTS_1 = (((1, 2), (3, 1)), ((1, 3), (2, 1), (3, 2)), ((2, 3),))
EX2_ROWS = (0b111, 0b100, 0b101)         # Example (ii): the source's near miss
EX2_PARTS_1 = (((2, 3), (3, 1)), ((1, 1), (3, 3)), ((1, 2),), ((1, 3),))
# Section 2: the four parts Algorithm 1 returns on the labelling printed for Gamma.
G_ALG1_PARTS_1 = (((2, 1), (3, 2)), ((1, 1), (3, 3)), ((1, 2), (2, 3), (3, 4)), ((2, 4),))

# Section 4: the Moebius ladder M_5 on {0,...,9}: the 10-cycle plus the five chords.
M5_G6 = 'IhEIHCPaG'
M5_EDGES = tuple([(i, (i + 1) % 10) for i in range(10)] + [(i, i + 5) for i in range(5)])
M5_ROWS = (21, 11, 22, 13, 26)           # X = evens 0,2,4,6,8;  Y = odds 1,3,5,7,9
M5_X = (0, 2, 4, 6, 8)
M5_Y = (1, 3, 5, 7, 9)

# Section 3: the census table printed in the paper.
#   order -> (classes under row/col permutation AND side exchange, witnesses)
CENSUS_SWAP = {2: (2, 0), 3: (3, 0), 4: (10, 0), 5: (18, 0), 6: (54, 0), 7: (128, 1)}
#   order -> (classes under row/col permutation only, witnesses)
CENSUS_ORIENT = {2: (2, 0), 3: (6, 0), 4: (15, 0), 5: (36, 0), 6: (92, 0), 7: (256, 2)}
# the refinement printed in the Remark after Theorem 3
ORDER7_NO_ISOLATED = {(1, 6): 1, (2, 5): 11, (3, 4): 42}

# the published values pmd(K_{m,n}) = m + n - 1, quoted in Section 5
KMN_RANGE = [(m, n) for m in range(1, 5) for n in range(m, 5)]

CHECKS = [0]
FAILED = [0]


def ck(name, cond, detail=''):
    CHECKS[0] += 1
    if cond:
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        FAILED[0] += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


# ---------------------------------------------------------------------------
# 0. PLUMBING: bipartite graphs as (m, n, rows)
# ---------------------------------------------------------------------------
def edges_of(m, n, rows):
    return tuple((i, j) for i in range(m) for j in range(n) if (rows[i] >> j) & 1)


def rows_from_edges(m, n, edges):
    rows = [0] * m
    for (i, j) in edges:
        rows[i] |= 1 << j
    return tuple(rows)


def size(m, n, rows):
    return sum(bin(r).count('1') for r in rows)


def degrees(m, n, rows):
    dx = [bin(rows[i]).count('1') for i in range(m)]
    dy = [sum(1 for i in range(m) if (rows[i] >> j) & 1) for j in range(n)]
    return dx, dy


def delta(m, n, rows):
    dx, dy = degrees(m, n, rows)
    return max(dx + dy + [0])


def transpose(m, n, rows):
    return (n, m, tuple(sum(((rows[i] >> j) & 1) << i for i in range(m)) for j in range(n)))


def remove(m, n, rows, edges):
    out = list(rows)
    for (i, j) in edges:
        out[i] &= ~(1 << j)
    return tuple(out)


def decode_graph6(s):
    """(order, sorted edge tuple).  graph6 as defined by McKay: the order, then the
    upper triangle read column by column, six bits to a printable character."""
    b = [ord(c) - 63 for c in s]
    if not b:
        raise ValueError('empty graph6 string')
    if b[0] < 63:
        n, rest = b[0], b[1:]
    else:
        raise ValueError('this decoder covers orders below 63 only')
    bits = []
    for v in rest:
        if not 0 <= v < 64:
            raise ValueError('character out of range')
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    need = n * (n - 1) // 2
    if len(bits) < need:
        raise ValueError('graph6 payload too short')
    E, t = [], 0
    for j in range(1, n):
        for i in range(j):
            if bits[t]:
                E.append((i, j))
            t += 1
    return n, tuple(sorted(E))


def adjacency(n, edges):
    adj = {v: set() for v in range(n)}
    for (a, b) in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def components(n, edges, vertices=None):
    adj = adjacency(n, edges)
    seen, comps = set(), []
    for v in (range(n) if vertices is None else sorted(vertices)):
        if v in seen:
            continue
        stack, comp = [v], []
        seen.add(v)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))
    return comps


def is_forest(n, edges, vertices):
    """|E| = |V| - c is acyclic for a graph with c components."""
    return len(set(edges)) == len(set(vertices)) - len(components(n, edges, vertices))


def induced(m, n, rows, M):
    """(vertex set on 0..m+n-1, edge set) of the subgraph of Gamma induced by V(M),
    with the Y-side shifted by m so both sides live in one vertex range."""
    vx = set(i for (i, _j) in M)
    vy = set(j for (_i, j) in M)
    V = sorted(vx) + sorted(m + j for j in vy)
    E = tuple(sorted((i, m + j) for (i, j) in edges_of(m, n, rows) if i in vx and j in vy))
    return V, E


# ---------------------------------------------------------------------------
# 1. ALGORITHM 1 OF THE SOURCE, AND THE INVARIANT K
# ---------------------------------------------------------------------------
def alg1(m, n, rows, reverse_ties=False):
    """The source's Algorithm 1 on the labelling given by the row/column order.

    Slope s(x_i y_j) = j - i.  In one round every x still carrying an edge offers its
    remaining edge of least slope (unique, since j |-> j - i is injective at fixed i);
    the offers are scanned in increasing order of slope and an offer is accepted unless
    its Y-endpoint has already been accepted in this round.  The accepted set is the
    next part; its edges are deleted and the next round begins.

    `reverse_ties` reverses the order of equal-slope offers; the output must not change,
    because two offers of equal slope s have distinct columns j = s + i.
    """
    rows = list(rows)
    parts = []
    while any(rows):
        offers = []
        for i in range(m):
            if rows[i]:
                j = min((jj for jj in range(n) if (rows[i] >> jj) & 1), key=lambda jj: jj - i)
                offers.append((j - i, -i if reverse_ties else i, i, j))
        offers.sort()
        used_y, part = set(), []
        for (_s, _t, i, j) in offers:
            if j in used_y:
                continue
            part.append((i, j))
            used_y.add(j)
        rows = list(remove(m, n, tuple(rows), part))
        parts.append(tuple(sorted(part)))
    return len(parts), tuple(parts)


def alg1_source_literal(m, n, rows):
    """The source's Algorithm 1 written out as its pseudocode literally reads.

    Statement for statement from the pseudocode block at lines 549-574 of the e-print
    source file (the block is reproduced in Section 1 of the paper).  The source's own
    variables are kept -- the round counter k, the live X-side X', the used Y-side Y',
    the part M_k and the least-slope class M -- and no step is merged or reordered:

        k <- 0
        while E(Gamma) != {}:
            k <- k + 1
            X' <- X \\ {isolated vertices of X}
            Y' <- {}
            M_k <- {}
            while X' != {}:
                s   <- min(Slopes(Gamma[X' u Y]))
                M   <- {xy in E(Gamma[X' u Y]) : slope(xy) = s}
                M_k <- M_k u {xy in M : y not in Y'}
                X'  <- X' \\ {x : xy in M}
                Y'  <- Y' u {y : xy in M}
            E(Gamma) <- E(Gamma) \\ M_k
        return M_1, ..., M_k

    This routine exists only to be compared with `alg1`, which is the same computation
    reorganised into one offer per x per round; the two are compared part for part in
    Step 2, exhaustively over every bipartite graph of order at most 6 under every
    labelling and over every labelling of all four graphs the paper names.
    """
    E = set(edges_of(m, n, rows))
    parts = []
    while E:
        Xp = set(i for (i, _j) in E)          # X minus the isolated vertices of X
        Yp = set()
        Mk = set()
        while Xp:
            live = [(i, j) for (i, j) in E if i in Xp]        # E(Gamma[X' u Y])
            s = min(j - i for (i, j) in live)                 # min Slopes(Gamma[X' u Y])
            M = [(i, j) for (i, j) in live if j - i == s]
            Mk |= set((i, j) for (i, j) in M if j not in Yp)
            Xp -= set(i for (i, _j) in M)
            Yp |= set(j for (_i, j) in M)
        E -= Mk
        parts.append(tuple(sorted(Mk)))
    return len(parts), tuple(parts)


def all_bipartite_graphs(max_order):
    """Every (m, n, rows) with m, n >= 1 and m + n <= max_order: all 2^(mn) matrices for
    every split, isolated vertices allowed, no reduction up to isomorphism."""
    for N in range(2, max_order + 1):
        for m in range(1, N):
            n = N - m
            for bits in itertools.product(range(1 << n), repeat=m):
                yield (m, n, tuple(bits))


def labelings(m, n, rows):
    """Every one of the m! n! relabelings of one orientation."""
    for pr in itertools.permutations(range(m)):
        permuted = [rows[pr[i]] for i in range(m)]
        for pc in itertools.permutations(range(n)):
            yield (m, n, tuple(sum(((r >> pc[j]) & 1) << j for j in range(n)) for r in permuted))


def K_min(m, n, rows, swap_sides=True):
    """K(Gamma) = min over all labelings of the number of parts Algorithm 1 returns.
    With swap_sides the two sides may also be exchanged (2 m! n! labelings)."""
    best = None
    orients = [(m, n, rows)]
    if swap_sides:
        orients.append(transpose(m, n, rows))
    for (mm, nn, rr) in orients:
        for lab in labelings(mm, nn, rr):
            k, _ = alg1(*lab)
            if best is None or k < best:
                best = k
    return best


def count_labelings(m, n, swap_sides=True):
    f = [1]
    for t in range(1, 12):
        f.append(f[-1] * t)
    one = f[m] * f[n]
    return 2 * one if swap_sides else one


# ---------------------------------------------------------------------------
# 2. POSITIVE MATCHINGS, TWO INDEPENDENT IMPLEMENTATIONS
# ---------------------------------------------------------------------------
def is_matching(M):
    xs = [i for (i, _j) in M]
    ys = [j for (_i, j) in M]
    return len(set(xs)) == len(xs) and len(set(ys)) == len(ys)


def _pendant_in(m, n, rows, N):
    """An edge of N that is pendant in the subgraph of Gamma induced by V(N), or None."""
    xs = set(i for (i, _j) in N)
    ys = set(j for (_i, j) in N)
    for (i, j) in N:
        dx = sum(1 for jj in ys if (rows[i] >> jj) & 1)
        dy = sum(1 for ii in xs if (rows[ii] >> j) & 1)
        if dx == 1 or dy == 1:
            return (i, j)
    return None


def positive_subsets(m, n, rows, M):
    """The criterion literally: EVERY non-empty subset of M induces a subgraph of Gamma
    carrying a pendant edge that belongs to the subset."""
    M = tuple(M)
    for r in range(1, len(M) + 1):
        for N in itertools.combinations(M, r):
            if _pendant_in(m, n, rows, N) is None:
                return False
    return True


def positive_peel(m, n, rows, M):
    """The same criterion as a peeling: repeatedly delete a pendant edge of the current
    set.  Structurally independent of positive_subsets (it visits one chain of subsets,
    not all of them); the two are checked against each other below."""
    cur = list(M)
    while cur:
        p = _pendant_in(m, n, rows, cur)
        if p is None:
            return False
        cur.remove(p)
    return True


def matchings(m, n, rows):
    """Every non-empty matching of the graph."""
    E = edges_of(m, n, rows)
    out = []

    def rec(k, ux, uy, cur):
        if k == len(E):
            if cur:
                out.append(tuple(cur))
            return
        rec(k + 1, ux, uy, cur)
        i, j = E[k]
        if not (ux >> i) & 1 and not (uy >> j) & 1:
            cur.append((i, j))
            rec(k + 1, ux | (1 << i), uy | (1 << j), cur)
            cur.pop()
    rec(0, 0, 0, [])
    return out


def positive_matchings(m, n, rows):
    return [M for M in matchings(m, n, rows) if positive_peel(m, n, rows, M)]


def pmd(m, n, rows0):
    """The positive matching decomposition number: the least p such that E splits into
    an ORDERED sequence M_1,...,M_p with M_i a positive matching of the graph left after
    M_1,...,M_{i-1}.  Exhaustive, memoised on the residual edge set."""
    memo = {}

    def f(rows):
        if not any(rows):
            return 0
        if rows in memo:
            return memo[rows]
        best = None
        for M in matchings(m, n, rows):
            if not positive_peel(m, n, rows, M):
                continue
            v = f(remove(m, n, rows, M))
            if best is None or 1 + v < best:
                best = 1 + v
        memo[rows] = best
        return best
    return f(tuple(rows0))


# ---------------------------------------------------------------------------
# 3. THE CENSUS
# ---------------------------------------------------------------------------
def canon(m, n, rows, swap_sides=True):
    """A canonical form under permuting rows, permuting columns and (optionally)
    exchanging the two sides.  Rows are sorted, so row permutations are quotiented out
    without enumerating them."""
    best = None
    cand = [(m, n, rows)]
    if swap_sides:
        cand.append(transpose(m, n, rows))
    for (mm, nn, rr) in cand:
        for pc in itertools.permutations(range(nn)):
            key = (mm, nn, tuple(sorted(sum(((r >> pc[j]) & 1) << j for j in range(nn)) for r in rr)))
            if best is None or key < best:
                best = key
    return best


def sweep(order, swap_sides):
    """All bipartite graphs with an ordered bipartition of total size `order`, isolated
    vertices allowed: for every split, every one of the 2^(mn) matrices.  Returns the
    class representatives and the witnesses (pmd < K)."""
    splits = [(m, order - m) for m in range(1, order) if (m <= order - m or not swap_sides)]
    cls = {}
    for (m, n) in splits:
        for bits in range(1 << (m * n)):
            rows = tuple((bits >> (i * n)) & ((1 << n) - 1) for i in range(m))
            c = canon(m, n, rows, swap_sides)
            if c not in cls:
                cls[c] = (m, n, rows)
    wit, oracle_bad, oracle_n = [], 0, 0
    for (m, n, rows) in cls.values():
        d = delta(m, n, rows)
        K = K_min(m, n, rows, swap_sides)
        p = pmd(m, n, rows)
        oracle_n += 1
        if not (d <= p <= K):
            oracle_bad += 1
        if p < K:
            wit.append((m, n, rows, p, K, canon(m, n, rows, swap_sides)))
    return cls, wit, oracle_n, oracle_bad


def no_isolated(m, n, rows):
    if any(r == 0 for r in rows):
        return False
    cm = 0
    for r in rows:
        cm |= r
    return cm == (1 << n) - 1


# ---------------------------------------------------------------------------
# THE CHECKS
# ---------------------------------------------------------------------------
def main():
    print('exact verification of "A Seven-Vertex Bipartite Graph Whose Positive Matching')
    print('Decomposition Number Lies Below Its Slope Invariant"')
    print('python %s' % sys.version.split()[0])

    # -- Step 1: the object exhibited in the paper -------------------------
    print('\n=== Step 1: the witness Gamma = Theta(2,2,4), as printed')
    m, n, rows = G_M, G_N, G_ROWS
    E0 = tuple(sorted((i - 1, j - 1) for (i, j) in G_EDGES_1))
    ck('bitmasks_match_printed_edge_list', edges_of(m, n, rows) == E0,
       'rows (%d,%d,%d) <-> %d printed edges' % (rows + (len(E0),)))
    order, g6E = decode_graph6(G_G6)
    ck('graph6_decodes_to_order_seven', order == 7, "'%s' -> n = %d" % (G_G6, order))
    ck('graph6_edge_set_is_the_printed_one', g6E == tuple(sorted(G_G6_EDGES)),
       '%d edges, label for label' % len(g6E))
    ck('graph6_agrees_with_the_bitmasks_under_the_printed_bijection',
       tuple(sorted((a, b - 3) for (a, b) in g6E)) == E0,
       '0,1,2 = x1,x2,x3 and 3,4,5,6 = y1,y2,y3,y4')
    ck('edge_count_is_eight', size(m, n, rows) == 8)
    ck('the_printed_sides_are_a_bipartition',
       all(0 <= a < 3 <= b < 7 for (a, b) in g6E), 'no decoded edge joins two x or two y')
    dx, dy = degrees(m, n, rows)
    ck('degree_sequence_as_printed', (tuple(dx), tuple(dy)) == ((2, 3, 3), (2, 2, 2, 2)),
       'x = %s, y = %s' % (tuple(dx), tuple(dy)))
    ck('max_degree_is_three', delta(m, n, rows) == 3)
    ck('gamma_is_connected', len(components(7, g6E)) == 1)
    # Theta(2,2,4): the two degree-3 vertices, and the three internally disjoint paths
    branch = [v for v in range(7) if len(adjacency(7, g6E)[v]) == 3]
    rest = [v for v in range(7) if v not in branch]
    restE = [(a, b) for (a, b) in g6E if a in rest and b in rest]
    comps = components(7, tuple(restE))
    comps = [c for c in comps if c and set(c) <= set(rest)]
    sizes = sorted(len(c) for c in comps if set(c) <= set(rest))
    ck('exactly_two_vertices_of_degree_three', len(branch) == 2 and sorted(branch) == [1, 2],
       'x2 and x3 (decoded labels %s)' % sorted(branch))
    ck('deleting_them_leaves_three_paths_with_1_1_3_vertices', sizes == [1, 1, 3],
       'internal vertex counts %s -> path lengths 2, 2, 4' % sizes)
    plen = sorted(s + 1 for s in sizes)
    ck('so_gamma_is_theta_2_2_4', plen == [2, 2, 4] and sum(plen) == size(m, n, rows)
       and 2 + sum(sizes) == 7,
       'path lengths %s, summing to |E| = %d, on 2 + %d vertices' % (plen, sum(plen), sum(sizes)))

    # -- Step 2: Algorithm 1, against the source's own printed output -------
    print('\n=== Step 2: the transcription of Algorithm 1, against the source examples')
    for tag, exrows, want1 in (('example_i', EX1_ROWS, EX1_PARTS_1), ('example_ii', EX2_ROWS, EX2_PARTS_1)):
        k, parts = alg1(3, 3, exrows)
        want = tuple(tuple(sorted((i - 1, j - 1) for (i, j) in P)) for P in want1)
        ck('%s_part_count_matches_the_source' % tag, k == len(want1),
           'K(tau) = %d, source prints %d' % (k, len(want1)))
        ck('%s_parts_match_the_source_edge_for_edge' % tag, parts == want,
           str([[(i + 1, j + 1) for (i, j) in P] for P in parts]))
        ck('%s_literal_pseudocode_also_matches_the_source_edge_for_edge' % tag,
           alg1_source_literal(3, 3, exrows) == (len(want1), want),
           'the pseudocode of lines 549-574 run statement for statement')
    # every run of alg1 returns an ordered partition of E into matchings
    bad_part, runs = 0, 0
    probe = [(G_M, G_N, G_ROWS), (3, 3, EX1_ROWS), (3, 3, EX2_ROWS)]
    for (mm, nn, rr) in probe:
        for swap in (False, True):
            orients = [(mm, nn, rr)] + ([transpose(mm, nn, rr)] if swap else [])
            for (am, an, ar) in orients:
                for lab in labelings(am, an, ar):
                    runs += 1
                    _k, parts = alg1(*lab)
                    flat = [e for P in parts for e in P]
                    if (sorted(flat) != sorted(edges_of(*lab)) or len(set(flat)) != len(flat)
                            or not all(is_matching(P) for P in parts)):
                        bad_part += 1
    ck('every_run_returns_an_ordered_partition_of_E_into_matchings', bad_part == 0,
       '%d runs, %d defects' % (runs, bad_part))
    bad_tie = 0
    for (mm, nn, rr) in probe:
        for lab in labelings(mm, nn, rr):
            if alg1(*lab) != alg1(*lab, reverse_ties=True):
                bad_tie += 1
    ck('algorithm_1_is_tie_order_independent_hence_deterministic', bad_tie == 0,
       'equal-slope offers have distinct columns j = s + i; %d disagreements' % bad_tie)
    # the streamlined transcription in the paper against the pseudocode written out literally
    lit_runs, lit_bad = 0, 0
    for (mm, nn, rr) in all_bipartite_graphs(6):
        for lab in labelings(mm, nn, rr):
            lit_runs += 1
            if alg1(*lab) != alg1_source_literal(*lab):
                lit_bad += 1
    named = [(G_M, G_N, G_ROWS), (3, 3, EX1_ROWS), (3, 3, EX2_ROWS), (5, 5, M5_ROWS)]
    for (mm, nn, rr) in named:
        for (am, an, ar) in ((mm, nn, rr), transpose(mm, nn, rr)):
            for lab in labelings(am, an, ar):
                lit_runs += 1
                if alg1(*lab) != alg1_source_literal(*lab):
                    lit_bad += 1
    ck('alg1_equals_the_literal_pseudocode_part_for_part', lit_bad == 0,
       'every bipartite graph of order <= 6 under every labelling, plus every labelling of '
       'Gamma, both source examples and M_5 in both orientations: %d runs, %d disagreements'
       % (lit_runs, lit_bad))

    # -- Step 3: the positivity criterion ----------------------------------
    print('\n=== Step 3: the positivity criterion, two independent implementations')
    dis, tested = 0, 0
    pool = [(G_M, G_N, G_ROWS), (3, 3, EX1_ROWS), (3, 3, EX2_ROWS), (2, 2, (3, 3)),
            (3, 3, (7, 7, 7)), (1, 4, (15,)), (4, 4, (15, 15, 15, 15)), (5, 5, M5_ROWS),
            (G_M, G_N, remove(G_M, G_N, G_ROWS, [(i - 1, j - 1) for (i, j) in PART1_1]))]
    for (mm, nn, rr) in pool:
        for M in matchings(mm, nn, rr):
            tested += 1
            if positive_subsets(mm, nn, rr, M) != positive_peel(mm, nn, rr, M):
                dis += 1
    ck('all_subsets_and_peeling_criteria_agree', dis == 0,
       '%d matchings over %d graphs, %d disagreements' % (tested, len(pool), dis))
    k33 = (7, 7, 7)
    perfect = [M for M in matchings(3, 3, k33) if len(M) == 3]
    ck('no_perfect_matching_of_K33_is_positive',
       len(perfect) == 6 and not any(positive_subsets(3, 3, k33, M) for M in perfect),
       '%d perfect matchings, none positive' % len(perfect))
    # Lemma 1: if Gamma[V(M)] is acyclic then M is positive
    lem, lem_n = 0, 0
    for (mm, nn, rr) in pool:
        for M in matchings(mm, nn, rr):
            V, Esub = induced(mm, nn, rr, M)
            if is_forest(mm + nn, Esub, V):
                lem_n += 1
                if not positive_subsets(mm, nn, rr, M):
                    lem += 1
    ck('lemma_1_acyclic_induced_subgraph_implies_positive', lem == 0,
       '%d matchings with acyclic induced subgraph, %d not positive' % (lem_n, lem))

    # -- Step 4: pmd(Gamma) = 3 --------------------------------------------
    print('\n=== Step 4: pmd(Gamma) = 3')
    P1 = tuple(sorted((i - 1, j - 1) for (i, j) in PART1_1))
    P2 = tuple(sorted((i - 1, j - 1) for (i, j) in PART2_1))
    P3 = tuple(sorted((i - 1, j - 1) for (i, j) in PART3_1))
    ck('the_three_printed_parts_are_matchings', all(is_matching(P) for P in (P1, P2, P3)),
       'sizes %d, %d, %d' % (len(P1), len(P2), len(P3)))
    ck('the_three_printed_parts_partition_E',
       tuple(sorted(P1 + P2 + P3)) == E0 and len(set(P1 + P2 + P3)) == 8)
    r1 = remove(m, n, rows, P1)
    r2 = remove(m, n, r1, P2)
    r3 = remove(m, n, r2, P3)
    pos1 = positive_subsets(m, n, rows, P1)
    ck('M1_is_positive_in_Gamma', pos1)
    # the Lemma 1 route: Gamma[V(M_1)] is the path y1-x2-y3-x3
    vx, vy = set(i for (i, _j) in P1), set(j for (_i, j) in P1)
    sub1 = sorted((i, j) for (i, j) in E0 if i in vx and j in vy)
    V1, E1sub = induced(m, n, rows, P1)
    ck('the_subgraph_induced_by_V_M1_is_a_path_on_four_vertices',
       sub1 == [(1, 0), (1, 2), (2, 2)] and is_forest(m + n, E1sub, V1),
       'y1-x2-y3-x3 (1-based %s), acyclic' % [(i + 1, j + 1) for (i, j) in sub1])
    # Gamma - M_1 is the printed path
    idx = {'x1': (0, 'x'), 'x2': (1, 'x'), 'x3': (2, 'x'),
           'y1': (0, 'y'), 'y2': (1, 'y'), 'y3': (2, 'y'), 'y4': (3, 'y')}
    walk = []
    for a, b in zip(PATH_AFTER_M1, PATH_AFTER_M1[1:]):
        (ia, ta), (ib, _tb) = idx[a], idx[b]
        walk.append((ia, ib) if ta == 'x' else (ib, ia))
    ck('Gamma_minus_M1_is_exactly_the_printed_path',
       tuple(sorted(walk)) == edges_of(m, n, r1) and len(set(PATH_AFTER_M1)) == 7,
       '-'.join(PATH_AFTER_M1))
    Vp, Ep = list(range(7)), tuple((i, 3 + j) for (i, j) in edges_of(m, n, r1))
    ck('that_path_is_a_forest_with_max_degree_two',
       delta(m, n, r1) == 2 and size(m, n, r1) == 6
       and len(components(7, Ep)) == 1 and is_forest(7, Ep, Vp),
       'connected, acyclic, 6 edges on 7 vertices')
    pos2 = positive_subsets(m, n, r1, P2)
    pos3 = positive_subsets(m, n, r2, P3)
    ck('M2_is_positive_in_Gamma_minus_M1', pos2)
    ck('M3_is_positive_in_Gamma_minus_M1_minus_M2', pos3)
    ck('nothing_is_left_after_the_third_part', not any(r3))
    ck('hence_pmd_at_most_three',
       pos1 and pos2 and pos3 and not any(r3) and tuple(sorted(P1 + P2 + P3)) == E0,
       'three parts, each positive in its own residual, exhausting E')
    ck('pmd_is_at_least_Delta_equals_three', delta(m, n, rows) == 3,
       'a matching holds at most one edge at each vertex, and the parts partition E')
    ck('exhaustive_search_gives_pmd_exactly_three', pmd(m, n, rows) == 3)
    pm = positive_matchings(m, n, rows)
    ms3 = [M for M in matchings(m, n, rows) if len(M) == 3]
    ck('Gamma_has_matchings_of_size_three_but_none_is_positive',
       len(ms3) > 0 and max(len(M) for M in pm) == 2,
       '%d matchings of size 3, largest positive matching has size 2' % len(ms3))
    ck('so_the_part_sizes_of_any_3_part_decomposition_are_forced_to_2_3_3',
       (len(P1), len(P2), len(P3)) == (2, 3, 3),
       'the first part has size <= 2 and 8 = 2+3+3 with three parts of size <= 3')

    # -- Step 5: K(Gamma) = 4 ----------------------------------------------
    print('\n=== Step 5: K(Gamma) = 4, exhaustively')
    labs = list(labelings(m, n, rows))
    ck('the_printed_orientation_has_144_labelings',
       len(labs) == 144 == count_labelings(m, n, swap_sides=False), '3! * 4! = 144')
    vals = sorted(set(alg1(*lab)[0] for lab in labs))
    ck('no_labeling_of_the_printed_orientation_gives_fewer_than_four_parts',
       min(vals) == 4, 'the set of values over 144 labelings is %s' % vals)
    tm, tn, trows = transpose(m, n, rows)
    tvals = sorted(set(alg1(*lab)[0] for lab in labelings(tm, tn, trows)))
    ck('nor_does_any_labeling_of_the_exchanged_orientation', min(tvals) == 4,
       'X = the 4-side: values %s over its 144 labelings' % tvals)
    want_parts = tuple(tuple(sorted((i - 1, j - 1) for (i, j) in P)) for P in G_ALG1_PARTS_1)
    ck('algorithm_1_on_the_printed_labeling_returns_the_printed_four_parts',
       alg1(m, n, rows) == (4, want_parts),
       str([[(i + 1, j + 1) for (i, j) in P] for P in alg1(m, n, rows)[1]]))
    ck('K_of_Gamma_is_4_in_the_per_orientation_reading', K_min(m, n, rows, swap_sides=False) == 4,
       '144 labelings')
    ck('K_of_Gamma_is_4_when_the_sides_may_be_exchanged', K_min(m, n, rows, swap_sides=True) == 4,
       '288 labelings')
    ck('THE_ANSWER_pmd_3_is_strictly_below_K_4', pmd(m, n, rows) == 3 < 4 == K_min(m, n, rows),
       'Question 1 of the source is answered YES')
    ck('the_published_bound_Delta_le_pmd_le_K_holds_for_Gamma',
       delta(m, n, rows) <= pmd(m, n, rows) <= K_min(m, n, rows), '3 <= 3 <= 4')

    # -- Step 6: the anti-controls -----------------------------------------
    print('\n=== Step 6: graphs that must NOT be witnesses')
    ck('source_example_ii_has_K_tau_4_on_the_printed_labeling', alg1(3, 3, EX2_ROWS)[0] == 4)
    ck('but_its_K_over_all_72_labelings_is_3_so_it_is_not_a_witness',
       K_min(3, 3, EX2_ROWS, True) == 3 == K_min(3, 3, EX2_ROWS, False) == pmd(3, 3, EX2_ROWS),
       'pmd = K = 3; the strict gap at a FIXED labelling is not a gap for K(Gamma)')
    for tag, mm, nn, rr, want in (('K22', 2, 2, (3, 3), 3), ('C6', 3, 3, EX1_ROWS, 3),
                                  ('K33', 3, 3, (7, 7, 7), 5), ('K14', 1, 4, (15,), 4)):
        p, K = pmd(mm, nn, rr), K_min(mm, nn, rr)
        ck('%s_has_pmd_equal_to_K' % tag, p == K == want, 'pmd = K = %d' % p)
    bad_kmn = [(mm, nn) for (mm, nn) in KMN_RANGE
               if pmd(mm, nn, tuple((1 << nn) - 1 for _ in range(mm))) != mm + nn - 1]
    ck('published_values_pmd_Kmn_equal_m_plus_n_minus_1_reproduced', not bad_kmn,
       '%d pairs 1<=m<=n<=4, %d mismatches' % (len(KMN_RANGE), len(bad_kmn)))

    # -- Step 7: the census ------------------------------------------------
    print('\n=== Step 7: the census of all bipartite graphs of order at most 7')
    target = canon(m, n, rows, True)
    done = {}
    for swap, table, tag in ((True, CENSUS_SWAP, 'sides_may_be_exchanged'),
                             (False, CENSUS_ORIENT, 'per_orientation')):
        tot_or, tot_bad = 0, 0
        for order in sorted(table):
            done[(order, swap)] = sweep(order, swap)
            cls, wit, on, ob = done[(order, swap)]
            tot_or += on
            tot_bad += ob
            want_c, want_w = table[order]
            ck('census_%s_order_%d' % (tag, order), (len(cls), len(wit)) == (want_c, want_w),
               '%d classes, %d witnesses' % (len(cls), len(wit)))
            if order == 7:
                if swap:
                    ck('the_unique_order_7_witness_is_Gamma_itself',
                       len(wit) == 1 and wit[0][5] == target and (wit[0][3], wit[0][4]) == (3, 4),
                       'canonical form %s, pmd 3, K 4' % (target,))
                else:
                    keys = set(canon(w[0], w[1], w[2], True) for w in wit)
                    ck('the_two_order_7_witnesses_are_the_two_orientations_of_Gamma',
                       keys == {target} and len(wit) == 2
                       and all((w[3], w[4]) == (3, 4) for w in wit),
                       'splits %s' % sorted((w[0], w[1]) for w in wit))
        ck('oracle_Delta_le_pmd_le_K_on_every_class_%s' % tag, tot_bad == 0,
           '%d classes of order <= 7, %d violations' % (tot_or, tot_bad))
    cls7, _w, _o, _b = done[(7, True)]
    counts = {}
    for (mm, nn, rr) in cls7.values():
        if no_isolated(mm, nn, rr):
            counts[(mm, nn)] = counts.get((mm, nn), 0) + 1
    ck('order_7_classes_with_no_isolated_vertex_as_printed', counts == ORDER7_NO_ISOLATED,
       '%s, total %d of %d' % (sorted(counts.items()), sum(counts.values()), len(cls7)))
    ck('computed_census_has_no_witness_below_order_7_and_one_at_order_7',
       all(CENSUS_SWAP[o][1] == 0 for o in (2, 3, 4, 5, 6)) and CENSUS_SWAP[7][1] == 1,
       'under both readings of K, and with isolated vertices allowed')

    # -- Step 8: the Moebius ladder M_5 ------------------------------------
    print('\n=== Step 8: the 3-regular witness M_5 on ten vertices')
    order5, g6M5 = decode_graph6(M5_G6)
    E5 = tuple(sorted((min(a, b), max(a, b)) for (a, b) in M5_EDGES))
    ck('M5_graph6_decodes_to_the_printed_fifteen_edges', order5 == 10 and g6M5 == E5,
       "'%s' -> 10 vertices, %d edges" % (M5_G6, len(E5)))

    adj5 = adjacency(10, E5)
    ck('M5_is_cubic', all(len(adj5[v]) == 3 for v in range(10)))
    ck('M5_is_bipartite_with_the_printed_sides',
       all((a % 2) != (b % 2) for (a, b) in E5), 'X = evens, Y = odds')
    ck('M5_rowmasks_match_the_printed_ones',
       rows_from_edges(5, 5, [(M5_X.index(a if a % 2 == 0 else b), M5_Y.index(b if b % 2 else a))
                              for (a, b) in E5]) == M5_ROWS, str(M5_ROWS))
    rot = {v: (v + 1) % 10 for v in range(10)}
    rotE = tuple(sorted((min(rot[a], rot[b]), max(rot[a], rot[b])) for (a, b) in E5))
    ck('rotation_by_one_is_an_automorphism_exchanging_the_two_sides',
       rotE == E5 and set(rot[v] for v in M5_X) == set(M5_Y))
    ck('M5_has_max_degree_three_and_pmd_four', delta(5, 5, M5_ROWS) == 3 and pmd(5, 5, M5_ROWS) == 4,
       'exhaustive over positive matchings of every residual')
    ck('M5_has_14400_labelings_per_orientation',
       count_labelings(5, 5, swap_sides=False) == 14400, '5! * 5!')
    ck('K_of_M5_is_5_in_both_readings',
       K_min(5, 5, M5_ROWS, False) == 5 == K_min(5, 5, M5_ROWS, True),
       '14,400 and 28,800 labelings')
    ck('M5_is_a_second_witness_pmd_4_below_K_5', pmd(5, 5, M5_ROWS) < K_min(5, 5, M5_ROWS),
       '4 < 5, and M_5 is 3-regular')

    print('\nNOTE SCOPE -- what this program does NOT cover.  It re-derives, in exact integer '
          'arithmetic and from the objects printed in the note, the quantities the note states: '
          'the witness Gamma and its pmd and K.  Steps 7 and 8 go beyond the note\'s claims: the '
          'census of all bipartite graphs of order at most 7 under both readings of K, the '
          'minimality and uniqueness of Gamma among them, and the second witness M_5 are '
          'computations of this program, asserted nowhere in the note.  NOT RE-RUN here: (a) the source\'s '
          'combinatorial characterisation of a positive matching, which is a theorem of the '
          'source and is taken as the definition throughout -- this program never touches the '
          'underlying algebraic definition; (b) the READING of the source\'s pseudocode.  Step 2 '
          'goes as far as a program here can: the pseudocode is written out statement for '
          'statement as `alg1_source_literal`, it reproduces both of the source\'s worked '
          'examples part for part, and it agrees with the streamlined form used everywhere else '
          'in this program on every bipartite graph of order <= 6 under every labelling and on '
          'every labelling of Gamma, the two source examples and M_5.  What no program here can settle '
          'is that those statements are the ones printed at lines 549-574 of the source file; '
          'that single comparison is left to the reader; (c) every order >= 8 except the single graph M_5, so nothing here '
          'says how many witnesses of order 8 or 9 exist; (d) any infinite family, and in '
          'particular no lower bound on K for r-regular bipartite graphs is proved or used '
          'anywhere in this program.')

    print('\nVERDICT: ALL %d CHECKS PASS' % CHECKS[0] if not FAILED[0]
          else '\nVERDICT: %d of %d CHECKS FAILED' % (FAILED[0], CHECKS[0]))
    return 1 if FAILED[0] else 0


if __name__ == '__main__':
    sys.exit(main())

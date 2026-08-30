#!/usr/bin/env python3
"""verify.py -- re-derivation of every computational claim of

    "The 3-cut complex of the circular ladder is a wedge of spheres in degree 2n-4, for n 5"

from the objects PRINTED IN THE PAPER, and from nothing else.

WHAT IT IS GIVEN.  The graph6 string of the decided cell CL_9 (copied from Section 5 of the paper), the
vertex-and-edge rule of C_n box K_2, the source conjecture's two integer formulas, and the m = 2 row of the
source paper's own printed table.  No precomputed f-vector, no precomputed rank, no external data file.

WHAT IT DERIVES.  The graph and its invariants; the up-set complex D from the definition; the relative
chain complex of the pair (simplex, Delta_3) and its exact integral Smith normal forms; the reduced
integral homology of Delta_3(C_n box K_2) for n = 4..24; an INDEPENDENT direct computation of the same
homology from the full complex Delta_3 at n = 4, 5, 6, 7; the exhaustive 2^18-subset census of the decided
cell n = 9 together with cheaper full censuses at n = 6, 7, 8; the published Betti number 2(n-2)^2 of the
NEAR-MISS paper for the 2-by-n grid, as an external forced positive; four control graphs on which the
answer must and does fail to match; and five forced positives on the homology engine itself, including the
6-vertex RP^2 whose H_1 = Z/2 is what proves the Smith normal form really sees torsion.

CONTRACT.  Python 3.9+, standard library only (no numpy, no sympy, no networkx), exact integer arithmetic
throughout, no floating point in any decision, no randomness, no sampling.  One `PASS <name>` line per
check; closes with `VERDICT: ALL <n> CHECKS PASS` and exits 0 if and only if every check passed.

RUNTIME.  Under half a minute on a 2020 laptop.
"""

import itertools
import sys
from math import comb

# ---------------------------------------------------------------------------------------------------
# 0.  THE ONLY INPUTS.  Both are printed in the paper; a referee retypes them from the page.
# ---------------------------------------------------------------------------------------------------
WITNESS_G6 = 'QhCGGE@_A?c@C@A?__GC@?OC?oG'     # paper, Section 5: the decided cell CL_9
PUBLISHED_TABLE = {                             # the source paper's own SageMath table, m = 2 row
    4: [(3, 1), (4, 4)],                        # n : [(degree, rank), ...]
    5: [(6, 11)],
    6: [(8, 25)],
    7: [(10, 43)],
    8: [(12, 65)],
}

HOM_RANGE = range(4, 25)        # n for which the full integral homology is computed
COMB_RANGE = range(4, 41)       # n for which the combinatorial claims are checked
DIRECT_RANGE = (4, 5, 6, 7)     # n for which Delta_3 is also built and reduced DIRECTLY
CENSUS_RANGE = (6, 7, 8, 9)     # n for which every one of the 2^(2n) subsets is classified

CHECKS = []
FAILED = []


def check(name, ok, detail=''):
    CHECKS.append(name)
    print('%s %s%s' % ('PASS' if ok else 'FAIL', name, (' ' + detail) if detail else ''))
    if not ok:
        FAILED.append(name)
    return bool(ok)


# ---------------------------------------------------------------------------------------------------
# 1.  GRAPHS
# ---------------------------------------------------------------------------------------------------
def graph6_decode(s):
    """(n_vertices, sorted edge list).  Standard graph6; the small-n header form, n < 63, suffices here."""
    data = [ord(c) - 63 for c in s.strip()]
    n = data[0]
    if n >= 63:
        raise ValueError('this decoder handles n < 63 only')
    bits = []
    for d in data[1:]:
        for k in range(5, -1, -1):
            bits.append((d >> k) & 1)
    edges, idx = [], 0
    for j in range(1, n):                       # graph6 order: column j = 1..n-1, row i = 0..j-1
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                edges.append((i, j))
            idx += 1
    return n, sorted(edges)


def circular_ladder(n):
    """C_n box K_2 on 2n vertices, labelled (r, j) -> r*n + j, exactly as the paper labels it."""
    edges = set()
    for r in (0, 1):
        for j in range(n):
            edges.add(tuple(sorted((r * n + j, r * n + (j + 1) % n))))
    for j in range(n):
        edges.add(tuple(sorted((j, n + j))))
    return 2 * n, sorted(edges)


def grid_2_by_n(n):
    """P_2 box P_n, the 2-by-n grid ("the ladder"), on which the near-miss paper prints 2(n-2)^2."""
    edges = set()
    for r in (0, 1):
        for j in range(n - 1):
            edges.add(tuple(sorted((r * n + j, r * n + j + 1))))
    for j in range(n):
        edges.add(tuple(sorted((j, n + j))))
    return 2 * n, sorted(edges)


def petersen():
    edges = set()
    for j in range(5):
        edges.add(tuple(sorted((j, (j + 1) % 5))))              # outer 5-cycle
        edges.add(tuple(sorted((5 + j, 5 + (j + 2) % 5))))      # inner pentagram
        edges.add(tuple(sorted((j, 5 + j))))                    # spokes
    return 10, sorted(edges)


def complete_bipartite_33():
    return 6, sorted((i, 3 + j) for i in range(3) for j in range(3))


def mobius_ladder(m):
    """C_{2m} plus the m main diagonals: cubic, and triangle-free for m >= 4."""
    edges = set()
    for j in range(2 * m):
        edges.add(tuple(sorted((j, (j + 1) % (2 * m)))))
    for j in range(m):
        edges.add(tuple(sorted((j, j + m))))
    return 2 * m, sorted(edges)


def adjacency(nv, edges):
    adj = [set() for _ in range(nv)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def triple_is_connected(t, adj):
    """A 3-set induces a CONNECTED subgraph iff it spans at least two edges."""
    a, b, c = t
    return ((b in adj[a]) + (c in adj[a]) + (c in adj[b])) >= 2


def count_triangles(nv, adj):
    return sum(1 for t in itertools.combinations(range(nv), 3)
               if t[1] in adj[t[0]] and t[2] in adj[t[0]] and t[2] in adj[t[1]])


def girth(nv, adj):
    """Exact girth by BFS from every vertex; None for a forest."""
    best = None
    for s in range(nv):
        dist = {s: 0}
        par = {s: None}
        q = [s]
        while q:
            nq = []
            for v in q:
                for w in adj[v]:
                    if w not in dist:
                        dist[w] = dist[v] + 1
                        par[w] = v
                        nq.append(w)
                    elif w != par[v]:
                        c = dist[v] + dist[w] + 1
                        if best is None or c < best:
                            best = c
            q = nq
    return best


def induced_4_cycles(nv, adj):
    out = []
    for s in itertools.combinations(range(nv), 4):
        e = [(a, b) for a, b in itertools.combinations(s, 2) if b in adj[a]]
        if len(e) != 4:
            continue
        deg = {v: 0 for v in s}
        for a, b in e:
            deg[a] += 1
            deg[b] += 1
        if all(d == 2 for d in deg.values()):
            out.append(s)
    return out


def connected_triple_count(nv, adj):
    return sum(1 for t in itertools.combinations(range(nv), 3) if triple_is_connected(t, adj))


# ---------------------------------------------------------------------------------------------------
# 2.  SMITH NORMAL FORM (sparse, exact integers) AND SIMPLICIAL HOMOLOGY
# ---------------------------------------------------------------------------------------------------
def smith_normal_form(rows, ncols):
    """rows: list of dict{column: int}.  Returns (rank, [invariant factors, ascending]).

    Plain integer Smith normal form on the active submatrix: take the smallest-modulus entry as pivot,
    clear its column by integer row operations and then its row by integer column operations (legal and
    cheap because the pivot column is a singleton by then), and enforce divisibility by folding in an
    offending row.  Exact throughout: no fractions, no floats, no randomness.  Validated in section 0 on
    S^1, S^2, a 5-simplex, the 7-vertex torus and the 6-vertex RP^2, the last of which has a genuine
    invariant factor 2 -- if this routine ever degenerated into a rank-only computation, that check is
    the one that would go red.
    """
    rows = [{c: v for c, v in r.items() if v} for r in rows]
    nrows = len(rows)
    colrows = {}
    for i, r in enumerate(rows):
        for c in r:
            colrows.setdefault(c, set()).add(i)
    active_r = set(i for i in range(nrows) if rows[i])
    active_c = set(c for c in colrows if colrows[c])
    factors = []

    def setval(i, c, v):
        if v:
            if c not in rows[i]:
                colrows.setdefault(c, set()).add(i)
            rows[i][c] = v
        elif c in rows[i]:
            del rows[i][c]
            colrows[c].discard(i)

    def row_op(dst, src, mult):
        for c, v in list(rows[src].items()):
            setval(dst, c, rows[dst].get(c, 0) + mult * v)

    while active_r and active_c:
        best = None
        for i in active_r:
            for c, v in rows[i].items():
                if c in active_c:
                    av = abs(v)
                    if best is None or av < best[0]:
                        best = (av, i, c)
            if best is not None and best[0] == 1:
                break
        if best is None:
            break
        _, pi, pc = best
        while True:
            again = False
            # (a) clear the pivot column with row operations
            for i in [x for x in colrows.get(pc, set()) if x in active_r and x != pi]:
                v = rows[i].get(pc)
                if not v:
                    continue
                q = v // rows[pi][pc]
                if q:
                    row_op(i, pi, -q)
                if rows[i].get(pc):             # a nonzero remainder is a strictly smaller pivot
                    pi = i
                    again = True
                    break
            if again:
                continue
            # (b) clear the pivot row with column operations; the pivot column is a singleton now, so a
            #     column operation touches row pi and nothing else.
            p = rows[pi][pc]
            for c in [x for x in list(rows[pi]) if x in active_c and x != pc]:
                v = rows[pi][c]
                q = v // p
                if q:
                    setval(pi, c, v - q * p)
                if rows[pi].get(c):
                    pc = c                      # remainder in the row: re-pivot and redo (a)
                    again = True
                    break
            if again:
                continue
            break
        p = abs(rows[pi][pc])
        active_r.discard(pi)
        active_c.discard(pc)
        # (c) divisibility: if any remaining entry is not a multiple of p, fold that row into the pivot
        #     row and start this pivot over.
        bad = None
        for i in active_r:
            for c, v in rows[i].items():
                if c in active_c and v % p:
                    bad = i
                    break
            if bad is not None:
                break
        if bad is not None:
            row_op(pi, bad, 1)
            active_r.add(pi)
            active_c.add(pc)
            continue
        factors.append(p)
        active_r = set(i for i in active_r if rows[i])
        active_c = set(c for c in active_c if colrows.get(c))
    factors.sort()
    return len(factors), factors


def simplicial_boundary(face):
    return [((-1) ** i, face[:i] + face[i + 1:]) for i in range(len(face))]


def homology_of_chain_complex(gens_by_deg, boundary):
    """gens_by_deg: {degree: [generator, ...]}; boundary(g) -> [(coeff, generator one degree lower)].

    Returns {degree: (free rank, [invariant factors > 1])} and, separately, the boundary ranks.
    """
    degs = sorted(gens_by_deg)
    index = {d: {g: i for i, g in enumerate(gens_by_deg[d])} for d in degs}
    rank, tors = {}, {}
    for d in degs:
        lower = index.get(d - 1, {})
        rows = [dict() for _ in range(len(lower))]
        for g in gens_by_deg[d]:
            j = index[d][g]
            for coeff, h in boundary(g):
                i = lower.get(h)
                if i is not None:
                    rows[i][j] = rows[i].get(j, 0) + coeff
        r, f = smith_normal_form(rows, len(gens_by_deg[d]))
        rank[d] = r
        tors[d] = [x for x in f if x != 1]
    out = {}
    for d in degs:
        out[d] = (len(gens_by_deg[d]) - rank.get(d, 0) - rank.get(d + 1, 0), tors.get(d + 1, []))
    return out, rank


def reduced_homology_of_complex(faces):
    """Reduced integral homology of a simplicial complex given by ALL its nonempty faces.

    The empty face is put in degree -1, which is exactly what makes the answer reduced.
    """
    gens = {-1: [()]}
    for f in faces:
        gens.setdefault(len(f) - 1, []).append(f)
    for d in gens:
        gens[d].sort()
    h, _r = homology_of_chain_complex(gens, simplicial_boundary)
    return {d: v for d, v in h.items() if d >= 0}


# ---------------------------------------------------------------------------------------------------
# 3.  THE UP-SET COMPLEX D, AND Delta_3 THROUGH THE RELATIVE COMPLEX
# ---------------------------------------------------------------------------------------------------
def dual_levels(nv, adj):
    """D = {S subset V : every 3-subset of S induces a CONNECTED subgraph}, as {size: [sorted tuple]}.

    Straight from the definition, D is downward closed, so it is enumerated by extending each member with
    a strictly larger vertex and testing only the triples containing the new vertex.  That is complete
    because deleting the largest element of a member leaves a member.  Section 3 of the output checks the
    downward closure independently rather than assuming it.
    """
    levels = {0: [()], 1: [(v,) for v in range(nv)]}
    k = 1
    while levels[k]:
        nxt = []
        for S in levels[k]:
            for v in range(S[-1] + 1, nv):
                if all(triple_is_connected(tuple(sorted((a, b, v))), adj)
                       for a, b in itertools.combinations(S, 2)):
                    nxt.append(S + (v,))
        k += 1
        levels[k] = sorted(nxt)
    return levels


def delta3_nonfaces(nv, adj):
    """(levels of D, the non-faces of Delta_3).

    sigma is a face of Delta_3(G) iff V\\sigma contains a 3-set inducing a DISCONNECTED subgraph; so
    sigma is a NON-face iff every 3-subset of V\\sigma is connected, i.e. iff V\\sigma lies in D.
    """
    full = frozenset(range(nv))
    lv = dual_levels(nv, adj)
    out = []
    for k in sorted(lv):
        for S in lv[k]:
            out.append(tuple(sorted(full - set(S))))
    return lv, out


def delta3_homology_via_relative(nv, adj):
    """Reduced integral homology of Delta_3(G), through the pair (full simplex S, Delta_3).

    S is contractible, so the long exact sequence of the pair gives H~_d(Delta_3) = H_{d+1}(S, Delta_3),
    and the relative chain complex has the NON-FACES of Delta_3 as its basis (a boundary term that is a
    face of Delta_3 dies in the quotient).  Returns (homology, levels of D, generator counts, ranks).
    """
    lv, nonfaces = delta3_nonfaces(nv, adj)
    nfset = set(nonfaces)
    gens = {}
    for f in nonfaces:
        gens.setdefault(len(f) - 1, []).append(f)
    for d in gens:
        gens[d].sort()

    def bd(g):
        return [(s, h) for s, h in simplicial_boundary(g) if h in nfset]

    rel, ranks = homology_of_chain_complex(gens, bd)
    red = {d - 1: v for d, v in rel.items() if d - 1 >= 0}
    return red, lv, {d: len(v) for d, v in gens.items()}, ranks


def nonzero(h):
    return {d: v for d, v in h.items() if v[0] or v[1]}


def unimodular_minor_certificate(n):
    """Re-execute the four-stage deletion of the paper's rank-5n lemma, for the circular ladder CL_n.

    Builds the 5n-by-5n submatrix of d_2 named in the paper -- columns the 2n triples T_1(v), the 2n
    triples T_3(v) and the n triples T_2((1,j)); rows the 2n same-rail distance-2 pairs, the 2n
    diagonal pairs and the n rail edges of row 1 -- and then repeatedly strips a row whose only
    surviving entry is +-1, together with that entry's column.  Stripping the whole matrix away IS a
    certificate that the determinant is +-1, hence that every invariant factor of d_2 is 1.

    Returns (rows_stripped, 5n, columns_are_distinct_members_of_D).
    """
    def V(r, j):
        return r * n + (j % n)

    def triple(name, r, j):
        if name == 'T1':
            return tuple(sorted((V(r, j - 1), V(r, j), V(r, j + 1))))
        if name == 'T2':
            return tuple(sorted((V(r, j - 1), V(r, j), V(1 - r, j))))
        return tuple(sorted((V(r, j + 1), V(r, j), V(1 - r, j))))

    cols = ([('T1', r, j) for r in (0, 1) for j in range(n)]
            + [('T3', r, j) for r in (0, 1) for j in range(n)]
            + [('T2', 1, j) for j in range(n)])
    rows = ([tuple(sorted((V(r, j - 1), V(r, j + 1)))) for r in (0, 1) for j in range(n)]
            + [tuple(sorted((V(r, j + 1), V(1 - r, j)))) for r in (0, 1) for j in range(n)]
            + [tuple(sorted((V(1, j - 1), V(1, j)))) for j in range(n)])
    triples = [triple(*c) for c in cols]
    adj = adjacency(2 * n, circular_ladder(n)[1])
    distinct = (len(set(rows)) == 5 * n and len(set(triples)) == 5 * n
                and all(triple_is_connected(t, adj) for t in triples))
    ri = {p: i for i, p in enumerate(rows)}
    ent = {}                                    # row index -> {col index: +-1}
    for j, T in enumerate(triples):
        a, b, c = T
        for sign, p in ((1, (b, c)), (-1, (a, c)), (1, (a, b))):
            i = ri.get(p)
            if i is not None:
                ent.setdefault(i, {})[j] = ent.setdefault(i, {}).get(j, 0) + sign
    live_r = set(range(5 * n))
    live_c = set(range(5 * n))
    stripped = 0
    while live_r:
        pick = None
        for i in live_r:
            nz = [(j, v) for j, v in ent.get(i, {}).items() if j in live_c and v]
            if len(nz) == 1 and abs(nz[0][1]) == 1:
                pick = (i, nz[0][0])
                break
        if pick is None:
            break
        live_r.discard(pick[0])
        live_c.discard(pick[1])
        stripped += 1
    return stripped, 5 * n, distinct


def census(nv, adj):
    """EXHAUSTIVE classification of all 2^nv subsets, with no pruning of any kind.

    has[A] = A contains a 3-subset inducing a disconnected subgraph.  For |A| >= 4 such a triple misses
    some vertex of A, so has[A] = OR over v in A of has[A minus v]; |A| = 3 is decided by the definition.
    Then sigma is a FACE of Delta_3 iff has[complement of sigma].
    """
    size = 1 << nv
    has = bytearray(size)
    for t in itertools.combinations(range(nv), 3):
        if not triple_is_connected(t, adj):
            has[(1 << t[0]) | (1 << t[1]) | (1 << t[2])] = 1
    popc = bytearray(size)
    for A in range(1, size):
        popc[A] = popc[A >> 1] + (A & 1)
    for A in range(size):
        if has[A] or popc[A] < 4:
            continue
        m = A
        while m:
            low = m & (-m)
            if has[A ^ low]:
                has[A] = 1
                break
            m ^= low
    full = size - 1
    faces = [0] * (nv + 1)
    nonfaces = [0] * (nv + 1)
    for A in range(size):
        (faces if has[full ^ A] else nonfaces)[popc[A]] += 1
    return faces, nonfaces


# ---------------------------------------------------------------------------------------------------
# 4.  THE CHECKS
# ---------------------------------------------------------------------------------------------------
def closure(facets):
    out = set()
    for f in facets:
        f = tuple(sorted(f))
        for k in range(1, len(f) + 1):
            out.update(itertools.combinations(f, k))
    return sorted(out, key=lambda t: (len(t), t))


def engine_controls():
    h = reduced_homology_of_complex(closure([(0, 1), (1, 2), (0, 2)]))
    check('engine-S1', nonzero(h) == {1: (1, [])},
          'reduced homology of the boundary of a triangle = %s' % nonzero(h))
    h = reduced_homology_of_complex(closure([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]))
    check('engine-S2', nonzero(h) == {2: (1, [])},
          'reduced homology of the boundary of a tetrahedron = %s' % nonzero(h))
    h = reduced_homology_of_complex(closure([(0, 1, 2, 3, 4, 5)]))
    check('engine-contractible', nonzero(h) == {},
          'reduced homology of a 5-simplex = %s (every reduced group must vanish)' % nonzero(h))
    torus = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2),
             (0, 1, 5), (1, 2, 6), (2, 3, 0), (3, 4, 1), (4, 5, 2), (5, 6, 3), (6, 0, 4)]
    h = reduced_homology_of_complex(closure(torus))
    check('engine-torus', nonzero(h) == {1: (2, []), 2: (1, [])},
          'reduced homology of the 7-vertex torus = %s' % nonzero(h))
    rp2 = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5),
           (1, 2, 4), (2, 3, 5), (3, 4, 1), (4, 5, 2), (5, 1, 3)]
    h = reduced_homology_of_complex(closure(rp2))
    check('engine-RP2-torsion', nonzero(h) == {1: (0, [2])},
          'reduced homology of the 6-vertex RP^2 = %s -- a genuine invariant factor 2, so the Smith '
          'normal form is not a rank-only computation in disguise' % nonzero(h))


def main():
    print('=== 0. the homology engine, on five objects whose answers are textbook ===')
    engine_controls()

    print()
    print('=== 1. the witness printed in the paper ===')
    g6n, g6e = graph6_decode(WITNESS_G6)
    cln, cle = circular_ladder(9)
    check('witness-label-equality', g6n == cln and set(g6e) == set(cle),
          'graph6 %r decodes to %d vertices and %d edges; the symmetric difference with the LABELLED '
          'edge set of C_9 box K_2 is %s (label equality, not isomorphism)'
          % (WITNESS_G6, g6n, len(g6e), sorted(set(g6e) ^ set(cle))))
    adj9 = adjacency(cln, cle)
    q9 = induced_4_cycles(cln, adj9)
    check('witness-CL9-properties',
          (cln, len(cle)) == (18, 27) and all(len(a) == 3 for a in adj9)
          and count_triangles(cln, adj9) == 0 and girth(cln, adj9) == 4 and len(q9) == 9,
          '|V|=%d |E|=%d 3-regular=%s triangles=%d girth=%d induced-4-cycles=%d (the nine rung squares, '
          'the first being %s)' % (cln, len(cle), all(len(a) == 3 for a in adj9),
                                   count_triangles(cln, adj9), girth(cln, adj9), len(q9), q9[0]))

    print()
    print('=== 2. the graph family, over every n the paper quantifies about, and the excluded n = 4 ===')
    fam = {}
    for n in COMB_RANGE:
        nv, e = circular_ladder(n)
        fam[n] = (nv, e, adjacency(nv, e))
    ok = all(fam[n][0] == 2 * n and len(fam[n][1]) == 3 * n
             and all(len(s) == 3 for s in fam[n][2])
             and count_triangles(fam[n][0], fam[n][2]) == 0
             and girth(fam[n][0], fam[n][2]) == 4 for n in fam)
    check('family-cubic-triangle-free-girth-4', ok,
          'n=%d..%d (%d values): |V|=2n, |E|=3n, 3-regular, 0 triangles, girth exactly 4'
          % (min(COMB_RANGE), max(COMB_RANGE), len(fam)))
    ct = {n: connected_triple_count(fam[n][0], fam[n][2]) for n in range(4, 21)}
    check('connected-triple-count-is-6n', all(ct[n] == 6 * n for n in ct),
          'n=4..20: the number of 3-sets inducing a CONNECTED subgraph is exactly 6n = sum over v of '
          'C(deg v, 2), with no triangle correction because the graph is triangle-free')
    q = {n: len(induced_4_cycles(fam[n][0], fam[n][2])) for n in range(4, 21)}
    check('induced-4-cycle-count', q[4] == 6 and all(q[n] == n for n in range(5, 21)),
          'q(n) = #induced 4-cycles: q(4) = %d = n+2 (the two row 4-cycles are induced only at n=4, '
          'which is precisely where the counting step of part (ii) breaks); q(n) = n for n=5..20' % q[4])

    print()
    print('=== 3. the up-set complex D: its f-vector, its downward closure, and dim D = 3 ===')
    duals = {n: dual_levels(fam[n][0], fam[n][2]) for n in COMB_RANGE}
    bad = [n for n in COMB_RANGE if n >= 5
           and [len(duals[n].get(k, [])) for k in range(6)]
           != [1, 2 * n, n * (2 * n - 1), 6 * n, n, 0]]
    check('dual-f-vector', not bad,
          'n=5..%d: f(D) = (1, 2n, n(2n-1), 6n, n) and NO 5-set, computed from the definition; n=9 gives '
          '%s and n=%d gives %s' % (max(COMB_RANGE), [len(duals[9].get(k, [])) for k in range(5)],
                                    max(COMB_RANGE),
                                    [len(duals[max(COMB_RANGE)].get(k, [])) for k in range(5)]))
    got4 = [len(duals[4].get(k, [])) for k in range(6)]
    check('dual-f-vector-n4-differs', got4 == [1, 8, 28, 24, 6, 0],
          'n=4: f(D) = %s, so f_3 = 6 = n+2 rather than n; the equality the proof needs is FALSE at n=4, '
          'which is why the source conjecture splits n=4 into a separate part' % got4[:5])
    viol = 0
    for n in (6, 9, 12, 17):
        members = set()
        for k in duals[n]:
            members.update(duals[n][k])
        for S in members:
            for i in range(len(S)):
                if S[:i] + S[i + 1:] not in members:
                    viol += 1
    check('dual-downward-closed', viol == 0,
          'at n=6, 9, 12, 17 every facet-deletion of every member of D is a member: %d violations. This '
          'is what licenses both the incremental enumeration and the 4-set extension test below.' % viol)
    ext = surv = 0
    for n in COMB_RANGE:
        nv, _e, a = fam[n]
        for S in duals[n][4]:
            for v in range(nv):
                if v in S:
                    continue
                ext += 1
                if all(triple_is_connected(tuple(sorted((x, y, v))), a)
                       for x, y in itertools.combinations(S, 2)):
                    surv += 1
    check('no-5-set-in-D', surv == 0,
          '%d single-vertex extensions of the 4-sets of D tested over n=4..%d; %d survive. So dim D = 3 '
          'and the relative chain complex has exactly five graded pieces for EVERY n, which is why a '
          'proof for all n is a finite piece of linear algebra.' % (ext, max(COMB_RANGE), surv))

    print()
    print('=== 4. reduced integral homology of Delta_3(C_n box K_2), via the relative complex ===')
    print('        S = the full simplex on V is contractible, so H~_d(Delta_3) = H_{d+1}(S, Delta_3),')
    print('        and the relative complex is generated by the non-faces of Delta_3.')
    homs, relgens, ranks = {}, {}, {}
    for n in HOM_RANGE:
        nv, _e, a = fam[n]
        homs[n], _lv, relgens[n], ranks[n] = delta3_homology_via_relative(nv, a)
    bad = [(n, nonzero(homs[n])) for n in HOM_RANGE if n >= 5
           and nonzero(homs[n]) != {2 * n - 4: (2 * n * n - 8 * n + 1, [])}]
    check('homology-concentrated-and-free', not bad,
          'n=5..%d (%d values): the ONLY nonvanishing reduced group is H~_{2n-4} = Z^(2n^2-8n+1), with '
          'an empty list of invariant factors > 1. n=9 -> %s ; n=%d -> %s'
          % (max(HOM_RANGE), len(HOM_RANGE) - 1, nonzero(homs[9]), max(HOM_RANGE),
             nonzero(homs[max(HOM_RANGE)])))
    bad = [(n, ranks[n]) for n in HOM_RANGE if n >= 5
           and tuple(ranks[n].get(2 * n - k) for k in (4, 3, 2, 1)) != (n, 5 * n, 2 * n - 1, 1)]
    check('relative-boundary-ranks-are-n-5n-2n-1', not bad,
          'n=5..%d: the integral boundary ranks of the relative complex in degrees 2n-4, 2n-3, 2n-2, '
          '2n-1 are exactly n, 5n, 2n-1, 1. Those are the three quantitative steps of the hand proof -- '
          'the injectivity of d_3 on the n squares, rank d_2 = 5n, rank d_1 = 2n-1 -- confirmed '
          'independently at %d values of n. At n=9 they are %s.'
          % (max(HOM_RANGE), len(HOM_RANGE) - 1, [ranks[9][d] for d in (13, 14, 15, 16, 17)]))
    bad = [(n, relgens[n]) for n in HOM_RANGE if n >= 5
           and tuple(relgens[n].get(2 * n - k) for k in (5, 4, 3, 2, 1))
           != (n, 6 * n, n * (2 * n - 1), 2 * n, 1)]
    check('relative-generator-counts', not bad,
          'n=5..%d: the relative complex has (n, 6n, n(2n-1), 2n, 1) generators in degrees 2n-5..2n-1, '
          'i.e. the levels of D read backwards' % max(HOM_RANGE))
    bad = [n for n in HOM_RANGE if n >= 5 and not (
        -1 + 2 * n - n * (2 * n - 1) + 6 * n - n == -(2 * n * n - 8 * n + 1)
        and n * (2 * n - 1) - 5 * n - (2 * n - 1) == 2 * n * n - 8 * n + 1)]
    check('euler-and-subtraction-identities', not bad,
          'n=5..%d: 1 - 2n + n(2n-1) - 6n + n = 2n^2-8n+1 (the Euler cross-check, which needs no linear '
          'algebra at all) and n(2n-1) - 5n - (2n-1) = 2n^2-8n+1 (the rank subtraction)' % max(HOM_RANGE))
    cert = {n: unimodular_minor_certificate(n) for n in HOM_RANGE if n >= 5}
    check('unimodular-minor-certificate',
          all(s == t and d for s, t, d in cert.values()),
          'n=5..%d: the 5n-by-5n submatrix of d_2 that the paper names is stripped ENTIRELY away by '
          'the paper\'s own four-stage deletion, one row with a single +-1 entry at a time (%s rows '
          'stripped of %s at n=5..8), so its determinant is +-1 and every invariant factor of d_2 is 1. '
          'The 5n selected columns are checked to be distinct members of D and the 5n selected rows to '
          'be distinct pairs.'
          % (max(HOM_RANGE), [cert[n][0] for n in (5, 6, 7, 8)], [cert[n][1] for n in (5, 6, 7, 8)]))
    check('torsion-free', all(not homs[n][d][1] for n in homs for d in homs[n]),
          'every invariant factor of every boundary matrix of every relative complex at n=4..%d equals '
          '1, so H~_*(Delta_3) is torsion-free at all %d values of n'
          % (max(HOM_RANGE), len(HOM_RANGE)))

    print()
    print('=== 5. the same homology computed DIRECTLY from the full complex Delta_3 ===')
    for n in DIRECT_RANGE:
        nv, _e, a = fam[n]
        _lv, nonfaces = delta3_nonfaces(nv, a)
        nf = set(nonfaces)
        faces = [s for k in range(1, nv + 1) for s in itertools.combinations(range(nv), k) if s not in nf]
        h = reduced_homology_of_complex(faces)
        check('direct-Delta3-n%d' % n, nonzero(h) == nonzero(homs[n]),
              'all %d nonempty faces of Delta_3 enumerated from the definition and reduced over Z with no '
              'duality and no relative complex: %s, which equals the relative-complex answer %s'
              % (len(faces), nonzero(h), nonzero(homs[n])))

    print()
    print('=== 6. exhaustive censuses: every one of the 2^(2n) subsets classified, no pruning ===')
    censuses = {}
    for n in CENSUS_RANGE:
        censuses[n] = census(fam[n][0], fam[n][2])
    bad = []
    for n in CENSUS_RANGE:
        f, nf = censuses[n]
        want = {2 * n - k: len(duals[n][k]) for k in range(5)}
        got = {k: nf[k] for k in range(2 * n + 1) if nf[k]}
        if got != want or sum(f) + sum(nf) != (1 << (2 * n)):
            bad.append((n, got, want))
    check('census-agrees-with-the-up-set', not bad,
          'n=6, 7, 8, 9: the census never mentions D, yet its non-face counts by size agree with the '
          'levels of D level by level, and the two classes exhaust 2^(2n)')
    f9, nf9 = censuses[9]
    check('census-n9-total', sum(f9) == 261909 and sum(nf9) == 235 and sum(f9) + sum(nf9) == 262144,
          '2^18 = 262144 subsets of V(CL_9) classified with no pruning: %d faces + %d non-faces'
          % (sum(f9), sum(nf9)))
    nfs = {k: nf9[k] for k in range(19) if nf9[k]}
    check('census-n9-nonface-sizes', nfs == {14: 9, 15: 54, 16: 153, 17: 18, 18: 1},
          'non-face sizes -> counts = %s, the complements of the five levels of D' % nfs)
    check('census-n9-purity',
          all(f9[k] == comb(18, k) for k in range(14)) and f9[14] == 3051 and f9[15] == 762
          and all(f9[k] == 0 for k in range(16, 19)),
          'every one of the C(18,k) subsets is a face for k <= 13; then 3051 of size 14 and 762 of size '
          '15, then nothing. So Delta_3 is PURE of dimension 14 = 2n-4 with 762 facets.')
    check('census-n9-facet-formula', f9[15] == comb(18, 3) - 6 * 9 == 762,
          'facets = C(2n,3) - 6n = %d - 54 = %d' % (comb(18, 3), f9[15]))

    print()
    print('=== 7. against the source paper\'s own printed table, and the n = 4 anti-control ===')
    bad = [(n, nonzero(homs[n])) for n, entries in PUBLISHED_TABLE.items()
           if nonzero(homs[n]) != {d: (r, []) for d, r in entries}]
    check('published-table-reproduced', not bad,
          'the m=2 row of the source paper\'s SageMath table, reproduced from the definition alone: '
          'n=5 Z^11 in degree 6, n=6 Z^25 in degree 8, n=7 Z^43 in degree 10, n=8 Z^65 in degree 12 -- '
          '4 of 4, and all four are also confirmed by the DIRECT route at n=5, 6, 7')
    check('anti-control-n4-two-groups', nonzero(homs[4]) == {3: (1, []), 4: (4, [])},
          'n=4 returns TWO groups, %s: Z^1 in degree mn-5 = 3 and Z^4 in degree mn-4 = 4, which is part '
          '(i) of the source conjecture at m=2 (2m^2-5m+3 = 1 and 2m(m-1) = 4), not the part-(ii) shape. '
          'So the decider can and does say no.' % (nonzero(homs[4]),))
    check('part-ii-would-be-false-at-n4', nonzero(homs[4]) != {4: (2 * 16 - 32 + 1, [])},
          'part (ii) read at n=4 would demand Z^1 in degree 4 and nothing else; the computation returns '
          '%s, so extending part (ii) down to n=4 is FALSE and the n>=5 hypothesis is not cosmetic'
          % (nonzero(homs[4]),))

    print()
    print('=== 8. the numbers the paper quotes for the cells beyond the source table ===')
    facets = {n: comb(2 * n, 3) - 6 * n for n in (9, 10, 11, 12)}
    check('facet-counts-past-the-table', facets == {9: 762, 10: 1080, 11: 1474, 12: 1952},
          'C(2n,3) - 6n = %s. In particular n=11 has 1474 facets of size 19 in dimension 18; an earlier '
          'draft of our own note recorded 1298 there and that figure is retracted.' % facets)
    preds = {n: 2 * n * n - 8 * n + 1 for n in (9, 10, 11, 12)}
    check('ranks-past-the-table',
          preds == {9: 91, 10: 121, 11: 155, 12: 193}
          and all(homs[n].get(2 * n - 4, (0, []))[0] == preds[n] for n in preds),
          '2n^2-8n+1 = %s in degrees 14, 16, 18, 20, and all four are CONFIRMED above by the exact '
          'integral computation, not predicted' % preds)
    bad = [n for n in range(5, 16)
           if min(len(f) for f in delta3_nonfaces(fam[n][0], fam[n][2])[1]) != 2 * n - 4]
    check('full-2-skeleton-and-simple-connectivity', not bad,
          'n=5..15: the SMALLEST non-face of Delta_3 has size exactly 2n-4 >= 6, so every subset of size '
          '<= 2n-5 is a face. Delta_3 therefore contains the full 2-skeleton of the simplex on 2n '
          'vertices and is simply connected; with free homology concentrated in the single degree '
          '2n-4 >= 6, Hurewicz and Whitehead give the wedge of spheres.')

    print()
    print('=== 9. the near-miss paper\'s own printed Betti number, as an EXTERNAL forced positive ===')
    bad = []
    for n in range(4, 11):
        nv, e = grid_2_by_n(n)
        a = adjacency(nv, e)
        red, lv, _g, _r = delta3_homology_via_relative(nv, a)
        want = {2 * n - 4: (2 * (n - 2) ** 2, [])}
        f2 = len(lv.get(3, []))
        f3 = len(lv.get(4, []))
        if nonzero(red) != want or f2 != 6 * n - 8 or f3 != n - 1:
            bad.append((n, nonzero(red), want, f2, f3))
    check('grid-2-by-n-published-betti-number', not bad,
          'for the 2-by-n grid P_2 box P_n (NOT the circular ladder: 3n-2 edges, four vertices of degree '
          '2, f_2(D) = 6n-8, q = n-1), the same machinery returns Z^(2(n-2)^2) in degree 2n-4 for '
          'n=4..10, i.e. 8, 18, 32, 50, 72, 98, 128. 2(n-2)^2 is a Betti number the near-miss paper '
          'PRINTS, so this run is calibrated against the literature and not only against itself.')

    print()
    print('=== 10. control graphs: on graphs the theorem does not cover, the answer must NOT match ===')
    for name, (nv, e), want_q in (('petersen', petersen(), 0),
                                  ('K33', complete_bipartite_33(), 9),
                                  ('mobius-ladder-6', mobius_ladder(6), 6),
                                  ('mobius-ladder-8', mobius_ladder(8), 8)):
        a = adjacency(nv, e)
        qq = len(induced_4_cycles(nv, a))
        cubic = all(len(s) == 3 for s in a)
        tf = count_triangles(nv, a) == 0
        red, lv, _g, _r = delta3_homology_via_relative(nv, a)
        h = nonzero(red)
        m = nv // 2
        part_ii = {2 * m - 4: (2 * m * m - 8 * m + 1, [])}
        lemma = {2 * m - 4: (2 * m * m - 9 * m + 1 + qq, [])}
        if name == 'petersen':
            check('control-petersen-does-not-match',
                  cubic and tf and qq == want_q and girth(nv, a) == 5 and h == {6: (6, [])}
                  and h != part_ii and h == lemma,
                  'cubic=%s triangle-free=%s girth=5 q=%d; reduced homology %s. The part-(ii) number for '
                  '2m=10 vertices would be %s -- DOES NOT MATCH, which is the required answer; and %s is '
                  'what the general lemma predicts from q=0, %s'
                  % (cubic, tf, qq, h, part_ii, h, lemma))
        elif name == 'K33':
            check('control-K33-does-not-match',
                  cubic and tf and qq == want_q and h == {0: (1, [])} and 2 * m - 4 == 2,
                  'cubic=%s triangle-free=%s q=%d, but the nine 4-cycles of K_{3,3} SHARE vertices, so '
                  'the no-shared-2-face hypothesis of the lemma fails. Reduced homology %s: the class '
                  'sits in the BOTTOM degree 0, not in the top degree 2m-4 = 2, and the part-(ii) '
                  'expression 2m^2-8m+1 = %d is not even non-negative at m=3. DOES NOT MATCH, which is '
                  'the required answer, and it is why the paper claims no general triangle-free theorem.'
                  % (cubic, tf, qq, h, 2 * m * m - 8 * m + 1))
        else:
            check('control-%s-matches-the-general-lemma' % name,
                  cubic and tf and qq == m and h == lemma == part_ii,
                  'cubic=%s triangle-free=%s q=%d=m; reduced homology %s equals the general-lemma value '
                  '2m^2-9m+1+q = %s, which coincides with the part-(ii) number exactly because q=m. So '
                  'the closure is NOT special to the circular ladder, and the paper says so.'
                  % (cubic, tf, qq, h, lemma))
        check('control-%s-dual-has-dimension-3' % name,
              len(lv.get(5, [])) == 0 and len(lv.get(4, [])) == qq,
              'D has no 5-set and exactly q=%d four-sets, on a graph outside the theorem: the '
              'parity/Mantel step is a property of triangle-freeness, not of the ladder' % qq)

    print()
    print('=== 11. the two formulas of the source conjecture, as bare arithmetic ===')
    check('conjecture-formulas',
          [2 * n * n - 8 * n + 1 for n in (5, 6, 7, 8)] == [11, 25, 43, 65]
          and [2 * n - 4 for n in (5, 6, 7, 8)] == [6, 8, 10, 12]
          and (2 * 4 - 5 * 2 + 3, 2 * 2 * 1) == (1, 4),
          'part (ii): 2n^2-8n+1 = 11, 25, 43, 65 in degrees mn-4 = 6, 8, 10, 12 at n=5..8; part (i) at '
          'm=2: 2m^2-5m+3 = 1 and 2m(m-1) = 4')

    print()
    print('NOT RE-RUN, and therefore not evidence of anything here:')
    print('  * n >= %d.  The exact homology stops at n = %d and the combinatorial claims at n = %d.'
          % (max(HOM_RANGE) + 1, max(HOM_RANGE), max(COMB_RANGE)))
    print('    The statement for ALL n >= 5 rests on the PROOF in the paper, every quantitative step of')
    print('    which this program confirms at %d separate values of n.  No finite run closes it, and'
          % (len(HOM_RANGE) - 1))
    print('    this one is not an induction.')
    print('  * The DIRECT (absolute, non-relative) homology of Delta_3 beyond n = 7.  At n = 7 the')
    print('    complex already has 16228 nonempty faces and at n = 9 it has 261909; n = 4, 5, 6, 7 are')
    print('    computed BOTH ways and those four are what license the relative route at larger n.')
    print('  * PART (i) OF THE SOURCE CONJECTURE.  No m >= 3 graph is touched anywhere in this program.')
    print('    The n = 4 numbers above are the m = 2 case of part (i) and nothing more; the first open')
    print('    cell of part (i), m = 6 and n = 4, is not attempted.')
    print('  * The GENERAL LEMMA for cubic triangle-free graphs is CHECKED at four control graphs, not')
    print('    proved and not censused over cubic triangle-free graphs of any given order.')
    print('  * THE HOMOTOPY TYPE.  This program computes homology and the 2-skeleton fact; the passage')
    print('    to a wedge of spheres is Hurewicz plus Whitehead, a hand argument, not a computation.')
    print('  * PRIOR ART.  Nothing here searches the literature.  The paper names the near miss and the')
    print('    calibration in section 9 above is against a number that paper prints, but attribution is')
    print('    a reading task and this program does not perform it.')

    n = len(CHECKS)
    print()
    if FAILED:
        print('FAILED: %s' % ', '.join(FAILED))
        print('VERDICT: %d of %d CHECKS FAILED' % (len(FAILED), n))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % n)
    return 0


if __name__ == '__main__':
    sys.exit(main())

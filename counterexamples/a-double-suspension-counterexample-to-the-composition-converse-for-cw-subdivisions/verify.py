#!/usr/bin/env python3
"""verify.py -- checks the computational claims of

    "A Double-Suspension Counterexample to the Composition Converse
     for Strong CW-Regular Subdivisions"

against the two objects PRINTED IN THAT PAPER and nothing else.  Table 1 of the paper
(the 90 tetrahedra of the 16-vertex triangulation P of the Poincare homology 3-sphere)
and Table 2 of the paper (the 53 non-identity labels of an A_5-valued flat connection on
the 2-skeleton of P) are transcribed below as the two data blocks TABLE_1 and TABLE_2:
the same entries in the same order, with the paper's column separators dropped and the
`u v: w` of Table 2 written `u v w`, so that `.split()` reads them.  Everything else in
this file is derived from those two blocks.

Python 3.9+, standard library only: no third-party package, no external data file, no
network.  All arithmetic is exact integer arithmetic; no floating-point value is computed
and no decision is taken on one.

The program prints one `PASS <name>` line per check and closes with

    VERDICT: ALL <n> CHECKS PASS

exiting 0 if and only if every check passed.  What it does NOT cover is printed as
`NOT RE-RUN:` lines immediately before the verdict, and repeated in REVIEW_NOTE.md
under `## Scope`.
"""

import sys
from itertools import combinations, permutations

# ---------------------------------------------------------------------------
# 0.  THE TWO OBJECTS, AS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------

# Table 1 of the paper: the 90 tetrahedra of P, in the paper's order and layout.
TABLE_1 = """
1 2 4 9    1 2 4 15    1 2 6 14    1 2 6 15    1 2 9 14
1 3 4 12   1 3 4 15    1 3 7 10    1 3 7 12    1 3 10 15
1 4 9 12   1 5 6 13    1 5 6 14    1 5 8 11    1 5 8 13
1 5 11 14  1 6 13 15   1 7 8 10    1 7 8 11    1 7 11 12
1 8 10 13  1 9 11 12   1 9 11 14   1 10 13 15  2 3 5 10
2 3 5 11   2 3 7 10    2 3 7 13    2 3 11 13   2 4 9 13
2 4 11 13  2 4 11 15   2 5 8 11    2 5 8 12    2 5 10 12
2 6 10 12  2 6 10 14   2 6 12 15   2 7 9 13    2 7 9 14
2 7 10 14  2 8 11 15   2 8 12 15   3 4 5 14    3 4 5 15
3 4 12 14  3 5 10 15   3 5 11 14   3 7 12 13   3 11 13 14
3 12 13 14 4 5 6 7     4 5 6 14    4 5 7 15    4 6 7 11
4 6 10 11  4 6 10 14   4 7 11 15   4 8 9 12    4 8 9 13
4 8 10 13  4 8 10 14   4 8 12 14   4 10 11 13  5 6 7 13
5 7 9 13   5 7 9 15    5 8 9 12    5 8 9 13    5 9 10 12
5 9 10 15  6 7 11 12   6 7 12 13   6 10 11 12  6 12 13 15
7 8 10 14  7 8 11 15   7 8 14 15   7 9 14 15   8 12 14 15
9 10 11 12 9 10 11 16  9 10 15 16  9 11 14 16  9 14 15 16
10 11 13 16  10 13 15 16  11 13 14 16  12 13 14 15  13 14 15 16
"""

# Table 2 of the paper: `u v w` for each edge uv of P whose label is not the identity;
# w is the one-line word of a permutation of {1,...,5}, i.e. w[i] is the image of i+1.
# EVERY EDGE OF P NOT LISTED HERE CARRIES THE IDENTITY, exactly as the paper says.
TABLE_2 = """
2 3 24153    2 5 35214    2 7 24153    2 8 35214
2 10 24153   2 11 35214   2 12 43521   2 13 41532
3 5 23451    3 11 23451   3 13 54213   3 14 23451
4 5 23451    4 6 23451    4 7 35214    4 8 41532
4 10 41532   4 11 35214   4 13 41532   4 14 23451
5 7 43521    5 9 25413    5 10 51234   5 12 25413
5 15 51234   6 7 43521    6 10 24153   6 11 43521
6 12 43521   7 9 31524    7 13 54213   7 14 31524
7 15 43152   8 9 25413    8 12 25413   8 14 31524
8 15 43152   9 10 35421   9 13 41532   9 15 35421
10 11 54132  10 12 54132  10 14 31524  10 16 54132
11 13 35421  11 15 43152  12 13 54213  12 14 23451
12 15 54213  13 14 54132  13 16 54132  14 15 35421
15 16 54132
"""

# The claims the paper makes, quoted here so that a disagreement is a FAIL and not a
# silent recomputation.
CLAIM_F_P = (16, 106, 180, 90)
CLAIM_F_KX = (19, 157, 546, 948, 810, 270)
CLAIM_CARD_Y = 393
CLAIM_CARD_X = 2751
CLAIM_RANKS = {'sigma': 2, 'tau': 4, 'tau_sigma': 6}
CLAIM_RHO_Y_BOTTOM = 2
CLAIM_RHO_Z_BOTTOM = 6
CLAIM_IMAGE_ORDER = 60


def parse_table_1(text):
    tok = text.split()
    assert len(tok) % 4 == 0, 'Table 1 does not consist of 4-element rows'
    return [tuple(sorted(int(t) for t in tok[i:i + 4])) for i in range(0, len(tok), 4)]


def parse_table_2(text):
    out = {}
    tok = text.split()
    assert len(tok) % 3 == 0, 'Table 2 does not consist of `u v w` triples'
    for i in range(0, len(tok), 3):
        u, v, w = int(tok[i]), int(tok[i + 1]), tok[i + 2]
        key = (min(u, v), max(u, v))
        assert key not in out, 'edge %s labelled twice in Table 2' % (key,)
        out[key] = tuple(int(c) - 1 for c in w)
    return out


# ---------------------------------------------------------------------------
# 1.  THE CHECK HARNESS
# ---------------------------------------------------------------------------
STATE = {'pass': 0, 'fail': 0}


def check(name, ok, detail=''):
    if ok:
        STATE['pass'] += 1
        print('PASS %s%s' % (name, ('  ' + detail) if detail else ''))
    else:
        STATE['fail'] += 1
        print('FAIL %s%s' % (name, ('  ' + detail) if detail else ''))
    return bool(ok)


# ---------------------------------------------------------------------------
# 2.  SIMPLICIAL TOOLKIT (exact, standard library only)
# ---------------------------------------------------------------------------
def faces_of(facets, include_empty=False):
    """All faces of the complex generated by `facets`, as sorted tuples."""
    S = set()
    for f in facets:
        sf = tuple(sorted(f))
        for k in range(1, len(sf) + 1):
            for c in combinations(sf, k):
                S.add(c)
    if include_empty:
        S.add(())
    return S


def f_vector(facets):
    S = faces_of(facets)
    d = max(len(f) for f in S)
    return tuple(sum(1 for f in S if len(f) == k) for k in range(1, d + 1))


def euler(facets):
    return sum((-1) ** k * n for k, n in enumerate(f_vector(facets)))


def reduced_euler_paper(facets):
    """The paper's reduced Euler characteristic (source line 3921 of the e-print):
    even-dimensional cells minus odd-dimensional cells, the empty cell included."""
    fv = f_vector(facets)
    tot = 1                                    # the empty cell, of dimension -1 (odd side? see sign)
    # dimension of a k-element face is k-1; the empty cell has dimension -1.
    even = sum(n for k, n in enumerate(fv, start=1) if (k - 1) % 2 == 0)
    odd = sum(n for k, n in enumerate(fv, start=1) if (k - 1) % 2 == 1) + tot
    return even - odd


def codim1_multiplicities(facets):
    """{codim-1 face: number of facets containing it} for a pure complex."""
    mult = {}
    for f in facets:
        sf = tuple(sorted(f))
        for i in range(len(sf)):
            g = sf[:i] + sf[i + 1:]
            mult[g] = mult.get(g, 0) + 1
    return mult


def dual_graph_connected(facets):
    mult = codim1_multiplicities(facets)
    nbr = {}
    for g, _ in mult.items():
        owners = [i for i, f in enumerate(facets) if set(g) <= set(f)]
        if len(owners) == 2:
            nbr.setdefault(owners[0], set()).add(owners[1])
            nbr.setdefault(owners[1], set()).add(owners[0])
    seen = {0}
    stack = [0]
    while stack:
        x = stack.pop()
        for y in nbr.get(x, ()):
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == len(facets)


def graph_connected(verts, edges):
    if not verts:
        return True
    nbr = {v: set() for v in verts}
    for u, v in edges:
        nbr[u].add(v)
        nbr[v].add(u)
    root = next(iter(verts))
    seen = {root}
    stack = [root]
    while stack:
        x = stack.pop()
        for y in nbr[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == len(verts)


def link(facets, face):
    """The link of `face` in the complex generated by `facets`, as a facet list."""
    s = set(face)
    out = []
    for f in facets:
        if s <= set(f):
            out.append(tuple(sorted(set(f) - s)))
    return out


def is_triangulated_2_sphere(facets):
    """Combinatorial 2-sphere: pure 2-dimensional, every edge in exactly two triangles,
    dual graph connected, Euler characteristic 2."""
    if not facets or any(len(f) != 3 for f in facets):
        return False
    if set(codim1_multiplicities(facets).values()) != {2}:
        return False
    if not dual_graph_connected(facets):
        return False
    return euler(facets) == 2


def is_single_cycle(facets):
    """Pure 1-dimensional, every vertex of degree 2, connected: a triangulated S^1."""
    if not facets or any(len(f) != 2 for f in facets):
        return False
    deg = {}
    for u, v in facets:
        deg[u] = deg.get(u, 0) + 1
        deg[v] = deg.get(v, 0) + 1
    if set(deg.values()) != {2}:
        return False
    return graph_connected(set(deg), facets) and len(facets) == len(deg)


# --- integer Smith normal form, used only for homology ----------------------
def smith_invariants(M):
    """Invariant factors of an integer matrix, by Smith normal form.  Pivots of absolute
    value 1 are preferred, which is why no entry ever grows on the matrices used here."""
    A = [row[:] for row in M]
    m = len(A)
    n = len(A[0]) if m else 0
    inv = []
    r = c = 0
    while r < m and c < n:
        piv = None
        for i in range(r, m):
            for j in range(c, n):
                if A[i][j]:
                    if abs(A[i][j]) == 1:
                        piv = (i, j)
                        break
                    if piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]]):
                        piv = (i, j)
            if piv is not None and abs(A[piv[0]][piv[1]]) == 1:
                break
        if piv is None:
            break
        pi, pj = piv
        A[r], A[pi] = A[pi], A[r]
        for row in A:
            row[c], row[pj] = row[pj], row[c]
        while True:
            p = A[r][c]
            for i in range(r + 1, m):
                if A[i][c]:
                    q = A[i][c] // p
                    if q:
                        for j in range(c, n):
                            A[i][j] -= q * A[r][j]
                    if A[i][c]:
                        A[r], A[i] = A[i], A[r]
                        p = A[r][c]
            for j in range(c + 1, n):
                if A[r][j]:
                    q = A[r][j] // p
                    if q:
                        for i in range(r, m):
                            A[i][j] -= q * A[i][c]
                    if A[r][j]:
                        for i in range(r, m):
                            A[i][c], A[i][j] = A[i][j], A[i][c]
                        p = A[r][c]
            if (all(A[i][c] == 0 for i in range(r + 1, m))
                    and all(A[r][j] == 0 for j in range(c + 1, n))):
                break
        inv.append(abs(A[r][c]))
        r += 1
        c += 1
    return inv


def integral_homology(facets):
    """{k: (betti_k, [torsion invariant factors])} of the simplicial chain complex."""
    S = faces_of(facets)
    d = max(len(f) for f in S) - 1
    F = [sorted(f for f in S if len(f) == k + 1) for k in range(d + 1)]
    rank = {0: 0, d + 1: 0}
    tors = {d + 1: []}
    for k in range(1, d + 1):
        idx = {f: i for i, f in enumerate(F[k - 1])}
        M = [[0] * len(F[k]) for _ in range(len(F[k - 1]))]
        for j, f in enumerate(F[k]):
            for t in range(len(f)):
                M[idx[f[:t] + f[t + 1:]]][j] += (-1) ** t
        iv = smith_invariants(M)
        rank[k] = len(iv)
        tors[k] = [x for x in iv if x > 1]
    return {k: (len(F[k]) - rank.get(k, 0) - rank.get(k + 1, 0), tors.get(k + 1, []))
            for k in range(d + 1)}


# ---------------------------------------------------------------------------
# 3.  A_5, WITH LEFT-TO-RIGHT COMPOSITION
# ---------------------------------------------------------------------------
# A permutation is a 5-tuple p with p[i] = image of i.  The product a*b is "a then b",
# i.e. (a*b)(i) = b(a(i)); with this convention the product of the labels along an edge
# path is a homomorphism on the edge-path group, read left to right.
ID = (0, 1, 2, 3, 4)


def mul(a, b):
    return (b[a[0]], b[a[1]], b[a[2]], b[a[3]], b[a[4]])


def inv(a):
    r = [0] * 5
    for i, x in enumerate(a):
        r[x] = i
    return tuple(r)


def is_even(p):
    return sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j]) % 2 == 0


def cycle_word(p):
    seen = set()
    out = ''
    for i in range(5):
        if i in seen:
            continue
        cyc = []
        j = i
        while j not in seen:
            seen.add(j)
            cyc.append(j + 1)
            j = p[j]
        if len(cyc) > 1:
            out += '(' + ''.join(str(x) for x in cyc) + ')'
    return out or 'e'


def generated_subgroup(gens):
    G = {ID}
    stack = [ID]
    while stack:
        x = stack.pop()
        for g in gens:
            y = mul(x, g)
            if y not in G:
                G.add(y)
                stack.append(y)
    return G


def spanning_tree(verts, edges):
    """A spanning tree of the graph, as a set of edges, or None if disconnected."""
    par = {v: v for v in verts}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    tree = set()
    for e in sorted(edges):
        a, b = find(e[0]), find(e[1])
        if a != b:
            par[a] = b
            tree.add(e)
    return tree if len(tree) == len(verts) - 1 else None


def propagate_from_tree(triangles, edges, verts):
    """Set every edge of one spanning tree to the identity and then apply, repeatedly,
    the deduction "two labels of a triangle determine the third".  Returns
    (number of labels forced, number of free generators left, all_forced_are_identity).
    Every step is an equation in the group, so the conclusion holds over ANY group."""
    tree = spanning_tree(verts, edges)
    lab = {e: (ID if e in tree else None) for e in edges}
    changed = True
    while changed:
        changed = False
        for (u, v, w) in triangles:
            a, b, c = (u, v), (v, w), (u, w)
            known = [k for k in (a, b, c) if lab[k] is not None]
            if len(known) != 2:
                continue
            if lab[a] is None:
                lab[a] = mul(lab[c], inv(lab[b]))
            elif lab[b] is None:
                lab[b] = mul(inv(lab[a]), lab[c])
            else:
                lab[c] = mul(lab[a], lab[b])
            changed = True
    forced = [e for e in edges if lab[e] is not None]
    return len(forced), len(edges) - len(forced), all(lab[e] == ID for e in forced)


# ---------------------------------------------------------------------------
# 4.  RUN THE CHECKS
# ---------------------------------------------------------------------------
def main():
    P = parse_table_1(TABLE_1)
    LAB2 = parse_table_2(TABLE_2)

    # ---------------- A.  the exhibited complex P -------------------------
    check('facet-table-shape',
          len(P) == 90 and len(set(P)) == 90 and all(len(f) == 4 for f in P),
          '90 distinct 4-element facets read from Table 1')

    verts = sorted({v for f in P for v in f})
    check('vertex-set', verts == list(range(1, 17)),
          'vertices are exactly 1..16')

    fP = f_vector(P)
    check('f-vector', fP == CLAIM_F_P, 'f(P) = %s, as claimed' % (fP,))

    check('euler-characteristic', euler(P) == 0,
          '16 - 106 + 180 - 90 = 0')

    check('reduced-euler-characteristic', reduced_euler_paper(P) == -1,
          "(16 + 180) - (1 + 106 + 90) = -1 = (-1)^3, the paper's convention (line 3921)")

    check('f2-equals-2f3', fP[2] == 2 * fP[3], '180 = 2 * 90')

    check('dehn-sommerville', fP[1] == fP[0] + fP[3], '106 = 16 + 90')

    mult = codim1_multiplicities(P)
    check('triangles-in-two-tetrahedra',
          len(mult) == 180 and set(mult.values()) == {2},
          'all 180 triangles lie in exactly 2 of the 90 tetrahedra')

    check('dual-graph-connected', dual_graph_connected(P),
          'the 90 facets form one connected dual graph')

    E = sorted({e for f in P for e in combinations(f, 2)})
    T = sorted({t for f in P for t in combinations(f, 3)})
    check('one-skeleton-connected', graph_connected(set(verts), E),
          '16 vertices, 106 edges, one component')

    bad = [v for v in verts if not is_triangulated_2_sphere(link(P, (v,)))]
    check('vertex-links-are-2-spheres', not bad,
          'all 16 vertex links are combinatorial 2-spheres, so P is a closed 3-manifold')

    bad = [e for e in E if not is_single_cycle(link(P, e))]
    check('edge-links-are-circles', not bad,
          'all 106 edge links are single cycles')

    HP = integral_homology(P)
    want = {0: (1, []), 1: (0, []), 2: (0, []), 3: (1, [])}
    check('homology-of-P', HP == want,
          'H_*(P) = (Z, 0, 0, Z) with no torsion, by Smith normal form over Z')

    betti = tuple(HP[k][0] for k in range(4))
    check('P-is-a-rational-homology-3-sphere',
          betti == (1, 0, 0, 1) and all(HP[k][1] == [] for k in range(4)),
          'Betti numbers (1,0,0,1) and no torsion, so P is a homology 3-sphere over Q, '
          "which is the hypothesis of the paper's ex:CWB0counter")

    # ---------------- B.  pi_1(P) is not trivial --------------------------
    A5 = [p for p in permutations(range(5)) if is_even(p)]
    A5set = set(A5)
    okgrp = (len(A5) == 60
             and all(mul(a, b) in A5set for a in A5 for b in A5)
             and all(mul(a, inv(a)) == ID and mul(inv(a), a) == ID for a in A5)
             and all(mul(ID, a) == a and mul(a, ID) == a for a in A5)
             and all(mul(mul(a, b), A5[7]) == mul(a, mul(b, A5[7])) for a in A5 for b in A5))
    check('A5-group-table', okgrp,
          '60 even permutations of {1..5}; the 3600 products all land inside it, the '
          'identity and all inverses are present and associativity is spot-checked')

    check('label-table-shape',
          len(LAB2) == 53 and all(e in set(E) for e in LAB2) and all(v != ID for v in LAB2.values()),
          '53 non-identity labels from Table 2, every one on an edge of P')

    check('labels-are-even-permutations',
          all(is_even(v) for v in LAB2.values()),
          'all 53 printed labels lie in A_5 (each is a 5-cycle)')

    lab = {e: LAB2.get(e, ID) for e in E}
    check('labelling-is-total',
          len(lab) == 106 and sum(1 for e in E if lab[e] == ID) == 53,
          'all 106 edges labelled; the 53 edges absent from Table 2 carry the identity')

    bad = [t for t in T if mul(lab[(t[0], t[1])], lab[(t[1], t[2])]) != lab[(t[0], t[2])]]
    check('flat-connection-on-180-triangles', not bad,
          'g_uv * g_vw = g_uw at every one of the 180 triangles, so holonomy is a '
          'homomorphism pi_1(|P|) -> A_5')

    idedges = {e for e in E if lab[e] == ID}
    tree = spanning_tree(set(verts), idedges)
    check('identity-edges-contain-spanning-tree',
          tree is not None and len(tree) == 15,
          'the identity-labelled edges span all 16 vertices, so the image of the '
          'holonomy homomorphism is generated by the 106 labels')

    G = generated_subgroup(set(lab.values()))
    check('holonomy-image-is-A5',
          len(G) == CLAIM_IMAGE_ORDER and set(G) == set(A5),
          'the labels generate a group of order %d, i.e. all of A_5' % len(G))

    loop = [(1, 2), (2, 3), (1, 3)]
    hol = mul(mul(lab[(1, 2)], lab[(2, 3)]), inv(lab[(1, 3)]))
    check('explicit-nontrivial-loop',
          all(e in set(E) for e in loop) and (1, 2, 3) not in set(T) and hol != ID,
          'the closed edge path 1-2-3-1 (whose vertex set is NOT a face of P) has '
          'holonomy %s != e' % cycle_word(hol))

    P_not_sphere = (hol != ID) and len(G) > 1
    check('P-is-not-a-sphere',
          P_not_sphere and HP == want,
          'pi_1(|P|) surjects onto A_5, so |P| is not simply connected and is not '
          'homeomorphic to S^3, although H_*(P) = H_*(S^3)')

    # ---------------- controls on the sphere side -------------------------
    BD4 = [tuple(sorted(c)) for c in combinations(range(1, 6), 4)]
    check('control-boundary-delta4-is-indistinguishable-homologically',
          f_vector(BD4) == (5, 10, 10, 5)
          and integral_homology(BD4) == HP
          and reduced_euler_paper(BD4) == reduced_euler_paper(P),
          'f(bd Delta^4) = (5,10,10,5); its integral homology is exactly that of P and '
          'both reduced Euler characteristics are -1')

    E4 = sorted({e for f in BD4 for e in combinations(f, 2)})
    T4 = sorted({t for f in BD4 for t in combinations(f, 3)})
    forced4, free4, allid4 = propagate_from_tree(T4, E4, set(range(1, 6)))
    forcedP, freeP, _ = propagate_from_tree(T, E, set(verts))
    check('control-boundary-delta4-admits-no-nontrivial-labelling',
          forced4 == 10 and free4 == 0 and allid4 and freeP == 53,
          'for bd Delta^4 the triangle relations force all 10 labels to the identity '
          '(0 free generators), so its edge-path group is trivial; for P, 53 generators '
          'survive the same deduction and Table 2 realises a surjection onto A_5')

    # ---------------- C.  the poset triple -------------------------------
    Y = sorted(faces_of(P, include_empty=True), key=lambda t: (len(t), t))
    check('Y-cardinality', len(Y) == CLAIM_CARD_Y,
          '|Y| = |face(P)| = 1 + 16 + 106 + 180 + 90 = %d' % len(Y))

    BD2 = sorted(faces_of([('a', 'b'), ('b', 'c'), ('a', 'c')], include_empty=True),
                 key=lambda t: (len(t), t))
    check('face-poset-of-boundary-delta2', len(BD2) == 7,
          'face(bd Delta^2) = bd B_3 has 7 elements: the empty face, 3 vertices, 3 edges')

    X = [(F, Gc) for F in Y for Gc in BD2]
    check('X-cardinality', len(X) == CLAIM_CARD_X == len(Y) * len(BD2),
          '|X| = 393 * 7 = %d' % len(X))

    rho_X = lambda x: len(x[0]) + len(x[1])
    rho_Y = lambda y: len(y) + 2
    RHO_Z = CLAIM_RHO_Z_BOTTOM

    cnt = {}
    for x in X:
        cnt[rho_X(x)] = cnt.get(rho_X(x), 0) + 1
    fKX = tuple(cnt[k] for k in range(1, 7))
    check('f-vector-of-KX-by-enumeration', fKX == CLAIM_F_KX,
          'f(K_X) = %s, counted over the 2751 pairs (F,G)' % (fKX,))

    conv = [0] * 7
    for i, a in enumerate((1, 16, 106, 180, 90)):
        for j, b in enumerate((1, 3, 3)):
            conv[i + j] += a * b
    check('f-vector-of-KX-by-convolution',
          tuple(conv[1:]) == CLAIM_F_KX and conv[0] == 1 and sum(conv) == CLAIM_CARD_X,
          '(1,16,106,180,90) * (1,3,3) = %s, sum 2751 = 393 * 7 -- an independent route '
          'to the same f-vector' % (tuple(conv),))

    check('dim-KX', max(rho_X(x) for x in X) - 1 == 5,
          'dim K_X = 3 + 1 + 1 = 5')

    check('euler-characteristic-KX',
          sum((-1) ** k * fKX[k] for k in range(6)) == 0,
          '19 - 157 + 546 - 948 + 810 - 270 = 0 = chi(S^5)')

    KXfacets = [tuple(sorted([str(v) for v in f] + list(g))) for f in P
                for g in (('a', 'b'), ('b', 'c'), ('a', 'c'))]
    m5 = codim1_multiplicities(KXfacets)
    check('KX-is-a-closed-5-pseudomanifold',
          len(KXfacets) == 270 and len(m5) == 810 and set(m5.values()) == {2},
          'each of the 810 4-faces lies in exactly two of the 270 5-simplices')

    check('rank-function-values',
          rho_X(((), ())) == 0 and rho_Y(()) == CLAIM_RHO_Y_BOTTOM and RHO_Z == 6,
          'rho_X(0hat) = 0, rho_Y(0hat) = 2 (natural rank shifted by +2), rho_Z(0hat) = 6')

    r_sigma = rho_Y(()) - rho_X(((), ()))
    r_tau = RHO_Z - rho_Y(())
    r_comp = RHO_Z - rho_X(((), ()))
    check('rank-arithmetic',
          (r_sigma, r_tau, r_comp) == (CLAIM_RANKS['sigma'], CLAIM_RANKS['tau'],
                                       CLAIM_RANKS['tau_sigma'])
          and r_sigma + r_tau == r_comp and r_sigma > 0,
          'rank(sigma) = 2, rank(tau) = 4, rank(tau.sigma) = 6 = 2 + 4; rank(sigma) > 0 '
          'places the triple outside the proved half of lem:composeCW')

    # covers of the product poset, which generate its order
    def covers_in(poset):
        s = set(poset)
        out = {}
        for a in poset:
            out[a] = [b for b in s if len(b) == len(a) + 1 and set(a) <= set(b)]
        return out

    covY = covers_in(Y)
    covB = covers_in(BD2)
    sig = lambda x: x[0]
    ncov = 0
    bad = 0
    for x in X:
        F, Gc = x
        for F2 in covY[F]:
            ncov += 1
            if not set(sig(x)) <= set(sig((F2, Gc))):
                bad += 1
        for G2 in covB[Gc]:
            ncov += 1
            if not set(sig(x)) <= set(sig((F, G2))):
                bad += 1
    check('sigma-order-preserving', bad == 0 and ncov == 11433,
          'sigma(F,G) = F is order-preserving on all %d covering pairs of X, and the '
          'covers generate the order' % ncov)

    check('sigma-rank-increasing',
          all(rho_X(x) <= rho_Y(x[0]) for x in X),
          'rho_X(x) <= rho_Y(sigma(x)) at all 2751 elements of X')

    check('sigma-surjective',
          {x[0] for x in X} == set(Y),
          'sigma hits all 393 elements of Y')

    sup = {F: [y for y in Y if set(F) <= set(y)] for F in Y}
    supB = {Gc: [H for H in BD2 if set(Gc) <= set(H)] for Gc in BD2}
    npair = sum(len(v) for v in sup.values())

    inst = 0
    bad = 0
    for F in Y:
        for y in sup[F]:
            for Gc in BD2:
                inst += 1
                if not any(len(y) + len(H) == rho_Y(y) for H in supB[Gc]):
                    bad += 1
    check('sigma-strongly-surjective', bad == 0 and inst == 23359,
          'for each of %d pairs (x, y >= sigma(x)) some x\' >= x has sigma(x\') = y and '
          'rho_X(x\') = rho_Y(y)' % inst)

    bad = [F for F in Y if not any(len(Fp) == 4 for Fp in sup[F])]
    check('tau-strongly-surjective', not bad and npair == 3337,
          'every one of the 393 faces of P lies in a tetrahedron (P is pure), which is '
          'the strong surjectivity of tau : Y -> B_0')

    sums = set()
    inst = 0
    for F in Y:
        for y in sup[F]:
            for Gc in BD2:
                s = sum((-1) ** (rho_Y(y) - (len(y) + len(H))) for H in supB[Gc])
                sums.add(s)
                inst += 1
    check('sigma-strong-formal', sums == {1} and inst == 23359,
          'all %d alternating sums of def:sfs for sigma equal 1' % inst)

    sums = set()
    for F in Y:
        sums.add(sum((-1) ** (RHO_Z - rho_Y(Fp)) for Fp in sup[F]))
    check('tau-strong-formal', sums == {1},
          'all 393 alternating sums of def:sfs for tau equal 1')

    quoted = [1 - 3 + 3, -1 + 2, 1, 1 - 16 + 106 - 180 + 90]
    quoted += [(-1) ** k * (-1) ** (4 - k) for k in range(1, 5)]
    check('paper-quoted-sfs-sums', set(quoted) == {1} and len(quoted) == 8,
          "the five sums displayed in the paper -- 1-3+3, -1+2, +1, 1-16+106-180+90 and "
          "(-1)^k(-1)^(4-k) for k = 1..4 -- all equal 1")

    # clause (ii): sigma is a strong CW-regular subdivision
    below = {}
    for y in Y:
        below[y] = [(F, Gc) for F in Y if set(F) <= set(y) for Gc in BD2]
    bad = []
    for y in Y:
        dim = max(rho_X(x) for x in below[y]) - 1
        if dim != rho_Y(y) - rho_X(((), ())) - 1:
            bad.append(y)
    sigma_dim_ok = not bad
    check('clause-ii-dimension', sigma_dim_ok,
          'dim K_{X,<=y} = rho_Y(y) - rho_X(0hat_X) - 1 at every one of the 393 y in Y')

    bad = []
    for y in Y:
        if y == ():
            continue
        strict = [(F, Gc) for F in Y if set(F) < set(y) for Gc in BD2]
        dimle = max(rho_X(x) for x in below[y]) - 1
        dimlt = max(rho_X(x) for x in strict) - 1
        interior = [x for x in below[y] if x not in set(strict)]
        fac = [tuple(sorted(map(str, set(F) | set(Gc)))) for (F, Gc) in strict
               if rho_X((F, Gc)) - 1 == dimlt]
        m = codim1_multiplicities(fac)
        chi = sum((-1) ** (rho_X(x) - 1) for x in strict if x != ((), ()))
        if not (dimlt == dimle - 1
                and len(interior) == 7
                and set(m.values()) == {2}
                and chi == (0 if dimlt % 2 else 2)):
            bad.append(y)
    sigma_bdry_ok = not bad
    check('clause-ii-boundary-identity', sigma_bdry_ok,
          'for all 392 y != 0hat_Y: K_{X,<y} has dimension dim K_{X,<=y} - 1, is a '
          'closed pseudomanifold with the Euler characteristic of S^{dim}, and '
          'K_{X,<=y} adds exactly the 7 faces (y,G) -- the boundary bookkeeping of '
          'ex:joinballs')

    bottom = [x for x in below[()]]
    bfac = [tuple(sorted(Gc)) for (F, Gc) in bottom if len(Gc) == 2]
    sigma_bottom_ok = (len(bottom) == 7 and is_single_cycle(bfac)
                       and max(rho_X(x) for x in bottom) - 1 == rho_Y(()) - 0 - 1)
    check('clause-ii-sphere-at-bottom', sigma_bottom_ok,
          'K_{X,<=0hat_Y} = bd Delta^2, a single 3-cycle, hence a regular CW-sphere of '
          'dimension 1 = rho_Y(0hat_Y) - rho_X(0hat_X) - 1')

    comp_dim_ok = (RHO_Z - rho_X(((), ())) - 1 == max(rho_X(x) for x in X) - 1 == 5)
    check('clause-iii-dimension', comp_dim_ok,
          'the single element 0hat_Z requires dim K_X = 6 - 0 - 1 = 5, and dim K_X = 5')

    tau_dim_ok = (RHO_Z - rho_Y(()) - 1 == len(max(Y, key=len)) - 1 == 3)
    check('clause-iv-dimension-passes', tau_dim_ok,
          'the dimension clause for tau READS 6 - 2 - 1 = 3 = dim P and PASSES, so the '
          'failure of tau is not a bookkeeping artefact')

    check('clause-iv-sphere-clause-fails', P_not_sphere,
          'K_{Y,<=0hat_Z} = P and |P| is not homeomorphic to S^3, so tau is not a '
          'strong CW-regular subdivision')

    # ⛔ THE ONE INPUT THIS PROGRAM DOES NOT COMPUTE, named as a constant so that it
    # cannot be mistaken for a check: |K_X| is homeomorphic to S^5 by Cannon's double
    # suspension theorem.  Everything else below is derived from the checks above.
    SPHERE_KX_CITED = True
    sigma_verdict = sigma_dim_ok and sigma_bdry_ok and sigma_bottom_ok
    comp_verdict = comp_dim_ok and SPHERE_KX_CITED
    tau_verdict = tau_dim_ok and not P_not_sphere
    check('decider-discriminates',
          (sigma_verdict, comp_verdict, tau_verdict) == (True, True, False),
          'one and the same clause-by-clause test, fed the three maps, returns YES on '
          'sigma (rank 2), YES on tau.sigma (rank 6, using the cited double suspension '
          'theorem) and NO on tau (rank 4): neither a constantly-YES nor a '
          'constantly-NO test reproduces this pattern')

    # r = 0 member of the family: the proved case, which must be silent
    X0 = [(F, ()) for F in Y]
    comp_verdict_r0 = tau_dim_ok and not P_not_sphere      # K_X = P when r = 0
    check('control-r0-is-silent',
          len(X0) == 393 and comp_verdict_r0 is False,
          'at r = 0 one has bd Delta^0 = empty, K_X = P, X = Y x B_0 of size 393, '
          'sigma = id and tau.sigma = tau; the SAME decider then returns NO on '
          'tau.sigma, so hypothesis (iii) fails, the family produces no counterexample '
          'at r = 0 and nothing here contradicts lem:composeCW, whose converse is '
          'proved exactly there')

    def binom(n, k):
        c = 1
        for i in range(k):
            c = c * (n - i) // (i + 1)
        return c

    ok = True
    detail = []
    for r in range(2, 7):
        bnd = [binom(r + 1, k) for k in range(r + 1)]          # f of face(bd Delta^r)
        card = sum(bnd)
        cv = [0] * (5 + len(bnd) - 1)
        for i, a in enumerate((1, 16, 106, 180, 90)):
            for j, b in enumerate(bnd):
                cv[i + j] += a * b
        rho_z = r + 4
        ok = ok and (card == 2 ** (r + 1) - 1
                     and sum(cv) == CLAIM_CARD_Y * card
                     and len(cv) - 1 == 3 + r + 1
                     and rho_z - 0 - 1 == 3 + r
                     and rho_z - r - 1 == 3
                     and {sum((-1) ** (r - k) * binom(r + 1 - j, k - j)
                              for k in range(j, r + 1)) for j in range(r + 1)} == {1})
        detail.append('r=%d:|X|=%d' % (r, sum(cv)))
    check('family-r-2-to-6', ok,
          'the same recipe with bd Delta^r, r = 2..6: |X| = 393(2^(r+1)-1) [%s], '
          'dim K_X = 3+r = rho_Z(0hat)-1, the clause-(iv) dimension is 3 = dim P for '
          'every r, and every def:sfs sum is 1 -- the witness is not tuned to r = 2'
          % ', '.join(detail))

    # ---------------- what this program does not cover --------------------
    print()
    print('NOT RE-RUN: the double suspension theorem, |P| * S^1 = |P| * S^0 * S^0 '
          'homeomorphic to S^5 (J. W. Cannon, Ann. of Math. 110 (1979) 83-112). This is '
          'the ONE input taken from the literature and the only support for the sphere '
          'clause of hypothesis (iii); nothing here reproves it, and no homeomorphism '
          'type is computed anywhere in this program.')
    print('NOT RE-RUN: ex:joinballs (that Delta^n * bd Delta^{n\'} is a ball with '
          'boundary bd Delta^n * bd Delta^{n\'}) is quoted from the target paper; the '
          'program checks its combinatorial consequences -- dimensions, boundary face '
          'sets, pseudomanifold conditions and Euler characteristics at all 393 y -- '
          'not the ball structure itself.')
    print('NOT RE-RUN: the classical identification of pi_1(|P|) with the binary '
          'icosahedral group of order 120. What is proved here is strictly weaker and '
          'strictly sufficient: pi_1(|P|) surjects onto A_5, hence is not trivial.')
    print('NOT RE-RUN: the case rank(sigma) = 1 of the Question. No computation in this '
          'program bears on it; it is OPEN, not empty.')
    print('NOT RE-RUN: any minimality. No census of smaller complexes, fewer vertices or '
          'other homology spheres was performed, and none is claimed; Bjoerner and Lutz '
          'already give non-PL 5-spheres on 18 vertices against the 19 of K_X.')
    print('NOT RE-RUN: the integral homology of K_X, and the exact-iff statement on the '
          'Z = B_0 slice. K_X is checked only to be a closed 5-pseudomanifold with '
          'chi = 0.')
    print('NOT RE-RUN: the literature. This program performs no search and no citation '
          'check.')
    print()
    if STATE['fail']:
        print('VERDICT: %d CHECK(S) FAILED of %d'
              % (STATE['fail'], STATE['fail'] + STATE['pass']))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % STATE['pass'])
    return 0


if __name__ == '__main__':
    sys.exit(main())

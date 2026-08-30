#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify.py -- machine check of the paper
    "Integer height-two rectangle visibility graphs are outerplanar"

WHAT THIS PROGRAM CHECKS.  Every object it consumes is PRINTED IN THE PAPER (Sections 1, 3 and 4);
it reads nothing from disk and imports nothing outside the Python standard library.  It

  (1) re-derives the visibility graph of each printed layout from the DEFINITION (a positive-width
      axis-parallel band meeting both interiors and no other closed rectangle), in exact integer /
      Fraction arithmetic, and compares it with what the paper says that graph is;
  (2) enumerates EVERY layout of integer rectangles whose bounding box is exactly [0,w] x [0,2],
      for w = 1..7 (711,202 layouts), and on each one checks
          - THEOREM A reproduces the definition edge-for-edge,
          - the proof of the main theorem goes through as written: split the FULL rectangles, the
            split layout's own definition-derived graph is drawn with pairwise non-crossing chords
            in the two-layer order of Step 2, and contracting the split pairs returns the original
            graph exactly,
          - no vertex with an independent neighbourhood has degree > 4;
  (3) runs the corollary's census: K_{1,5} occurs in no height-2 layout, K_{1,4} occurs and always
      with a FULL centre, and -- a refinement of the certificate's lemma that this run found --
      K_{1,3} occurs with BOTTOM and TOP centres too, so the "FULL centre" clause is needed only
      at m = 4;
  (4) runs four ANTI-CONTROLS, so that a vacuous pass is distinguishable from a real one: the
      outerplanarity routine must REJECT K_4 and K_{2,3}; a height-THREE layout must produce a K_4
      (so the height-2 hypothesis is load-bearing and the machinery can see a failure); and the
      y-rescaling of that K_4 layout must have bounding-box height exactly 2 with the SAME
      visibility graph, which is why integrality is load-bearing in the theorem statement.

Exit status is 0 if and only if every check passed.
"""

import sys
from fractions import Fraction
from itertools import combinations, permutations

CHECKS = []
FAILED = []


def check(name, ok, detail=''):
    CHECKS.append(name)
    line = '%s %s%s' % ('PASS' if ok else 'FAIL', name, ('  [%s]' % detail) if detail else '')
    print(line)
    if not ok:
        FAILED.append(name)


# ---------------------------------------------------------------------------
# 1.  THE DEFINITION, TWICE
# ---------------------------------------------------------------------------
# A rectangle is a 4-tuple (x1, x2, y1, y2) with x1 < x2 and y1 < y2, meaning the closed set
# [x1,x2] x [y1,y2].  A layout is a list of rectangles with pairwise disjoint INTERIORS (shared
# edges are allowed -- source line 96).  Two rectangles are adjacent in the visibility graph iff
# some axis-parallel band of POSITIVE WIDTH joins them and meets no other CLOSED rectangle.
#
# Both implementations below are EXACT and COMPLETE, not samples.  The set of abscissae at which a
# vertical band joining A and B is unobstructed is
#         (x1A, x2A) n (x1B, x2B)  \  U_C [x1C, x2C],
# an open set minus finitely many closed intervals whose endpoints all lie in the finite set of
# rectangle coordinates.  Such a set is nonempty iff it contains the midpoint of two consecutive
# coordinates.  So testing midpoints decides positive-width visibility outright.

def vis_by_columns(rects, w, h):
    """The definition, specialised to INTEGER coordinates in a w x h box: test the unit columns and
    unit rows, i.e. the midpoints (c+1/2) and (r+1/2).  Pure integer arithmetic."""
    E = set()
    for c in range(w):                                   # vertical bands
        col = [i for i, (x1, x2, y1, y2) in enumerate(rects) if x1 <= c and x2 >= c + 1]
        col.sort(key=lambda i: rects[i][2])              # y-interiors are disjoint here
        for u, v in zip(col, col[1:]):
            E.add((min(u, v), max(u, v)))
    for r in range(h):                                   # horizontal bands
        row = [i for i, (x1, x2, y1, y2) in enumerate(rects) if y1 <= r and y2 >= r + 1]
        row.sort(key=lambda i: rects[i][0])
        for u, v in zip(row, row[1:]):
            E.add((min(u, v), max(u, v)))
    return E


def vis_general(rects):
    """The definition again, for ARBITRARY rational coordinates and with no box assumed: test the
    midpoint of every pair of consecutive distinct x-coordinates, and likewise in y.  Exact
    Fraction arithmetic; no float ever decides anything."""
    R = [tuple(Fraction(v) for v in r) for r in rects]
    xs = sorted(set([r[0] for r in R] + [r[1] for r in R]))
    ys = sorted(set([r[2] for r in R] + [r[3] for r in R]))
    E = set()
    for a, b in zip(xs, xs[1:]):
        m = (a + b) / 2
        col = [i for i, (x1, x2, y1, y2) in enumerate(R) if x1 < m < x2]
        col.sort(key=lambda i: R[i][2])
        for u, v in zip(col, col[1:]):
            E.add((min(u, v), max(u, v)))
    for a, b in zip(ys, ys[1:]):
        m = (a + b) / 2
        row = [i for i, (x1, x2, y1, y2) in enumerate(R) if y1 < m < y2]
        row.sort(key=lambda i: R[i][0])
        for u, v in zip(row, row[1:]):
            E.add((min(u, v), max(u, v)))
    return E


def interiors_disjoint(rects):
    for (a1, a2, b1, b2), (c1, c2, d1, d2) in combinations(
            [tuple(Fraction(v) for v in r) for r in rects], 2):
        if min(a2, c2) > max(a1, c1) and min(b2, d2) > max(b1, d1):
            return False
    return True


# ---------------------------------------------------------------------------
# 2.  THEOREM A
# ---------------------------------------------------------------------------
BOTTOM, TOP, FULL = (0, 1), (1, 2), (0, 2)


def theorem_a(rects):
    """E = path(B* in x-order) u path(T* in x-order)
           u { bt : b BOTTOM, t TOP, open x-intervals overlap positively }."""
    bstar = sorted([i for i, r in enumerate(rects) if r[2] == 0], key=lambda i: rects[i][0])
    tstar = sorted([i for i, r in enumerate(rects) if r[3] == 2], key=lambda i: rects[i][0])
    E = set()
    for u, v in zip(bstar, bstar[1:]):
        E.add((min(u, v), max(u, v)))
    for u, v in zip(tstar, tstar[1:]):
        E.add((min(u, v), max(u, v)))
    bo = [i for i, r in enumerate(rects) if (r[2], r[3]) == BOTTOM]
    to = [i for i, r in enumerate(rects) if (r[2], r[3]) == TOP]
    for i in bo:
        for j in to:
            if min(rects[i][1], rects[j][1]) > max(rects[i][0], rects[j][0]):
                E.add((min(i, j), max(i, j)))
    return E


# ---------------------------------------------------------------------------
# 3.  OUTERPLANARITY
# ---------------------------------------------------------------------------
# A graph is outerplanar iff its vertices admit a CYCLIC ORDER in which no two edges, drawn as
# chords of the circle, cross properly.  Two routes are used and cross-checked:
#   * the PROOF's route -- the two-layer order of Step 2 on the FULL-split layout, then contraction;
#   * an INDEPENDENT route -- brute-force search over all cyclic orders (used on the printed objects
#     and on every layout of width <= 4), which knows nothing about rectangles.

def chords_non_crossing(order, E):
    pos = {v: k for k, v in enumerate(order)}
    ch = []
    for a, b in E:
        p, q = pos[a], pos[b]
        ch.append((min(p, q), max(p, q)))
    for i in range(len(ch)):
        a, b = ch[i]
        for j in range(i + 1, len(ch)):
            c, d = ch[j]
            if len({a, b, c, d}) < 4:
                continue
            if (a < c < b) != (a < d < b):
                return False
    return True


def outerplanar_by_search(n, E):
    """Independent of everything else in this file: is there ANY cyclic order with no crossing?"""
    if n <= 3:
        return True
    for p in permutations(range(1, n)):
        if p[0] > p[-1]:
            continue                                     # each cyclic order once, up to reflection
        if chords_non_crossing((0,) + p, E):
            return True
    return False


def split_fulls(rects):
    """Step 1: replace each FULL by its bottom and top halves.  Returns (layout, owner-map)."""
    out, owner = [], []
    for i, (x1, x2, y1, y2) in enumerate(rects):
        if (y1, y2) == FULL:
            out.append((x1, x2, 0, 1)); owner.append(i)
            out.append((x1, x2, 1, 2)); owner.append(i)
        else:
            out.append((x1, x2, y1, y2)); owner.append(i)
    return out, owner


def two_layer_order(split):
    """Step 2's drawing: bottoms left-to-right, then tops right-to-left."""
    b = sorted([k for k, r in enumerate(split) if (r[2], r[3]) == BOTTOM], key=lambda k: split[k][0])
    t = sorted([k for k, r in enumerate(split) if (r[2], r[3]) == TOP], key=lambda k: split[k][0])
    return b + t[::-1]


def proof_goes_through(rects, w, E):
    """Steps 1-3 of the paper's proof, re-derived on this layout.  None on success, else the reason."""
    split, owner = split_fulls(rects)
    esp = vis_by_columns(split, w, 2)                    # the DEFINITION again, on the split layout
    for i in set(owner):
        ks = [k for k, o in enumerate(owner) if o == i]
        if len(ks) == 2 and (min(ks), max(ks)) not in esp:
            return 'split pair %d is not an edge of G\', so the contraction is illegal' % i
    order = two_layer_order(split)
    if sorted(order) != list(range(len(split))):
        return 'the two-layer order is not a permutation of V(G\')'
    if not chords_non_crossing(order, esp):
        return 'the two-layer drawing of G\' has a crossing'
    contracted = set()
    for a, b in esp:
        u, v = owner[a], owner[b]
        if u != v:
            contracted.add((min(u, v), max(u, v)))
    if contracted != E:
        return 'contracting G\' does not return G'
    return None


# ---------------------------------------------------------------------------
# 4.  LAYOUT ENUMERATION
# ---------------------------------------------------------------------------
def layouts(w, h):
    """Every layout of integer rectangles with interiors pairwise disjoint whose bounding box is
    EXACTLY [0,w] x [0,h].  Rectangles are placed in the order of their minimal cell in column-major
    order, so each layout is generated exactly once; a cell may also be left uncovered."""
    cells = [(c, r) for c in range(w) for r in range(h)]
    idx = {cl: i for i, cl in enumerate(cells)}
    n = len(cells)
    acc = []

    def rec(i, covered, rects):
        while i < n and (covered >> i) & 1:
            i += 1
        if i == n:
            acc.append(tuple(rects))
            return
        c, r = cells[i]
        rec(i + 1, covered | (1 << i), rects)            # leave this cell uncovered
        for y2 in range(r + 1, h + 1):
            for x2 in range(c + 1, w + 1):
                mask, ok = 0, True
                for cc in range(c, x2):
                    for rr in range(r, y2):
                        bit = 1 << idx[(cc, rr)]
                        if covered & bit:
                            ok = False
                            break
                        mask |= bit
                    if not ok:
                        break
                if not ok:
                    break                                # wider is hopeless too
                rects.append((c, x2, r, y2))
                rec(i + 1, covered | mask, rects)
                rects.pop()

    rec(0, 0, [])
    out = []
    for L in acc:
        if not L:
            continue
        if min(r[0] for r in L) or max(r[1] for r in L) != w:
            continue
        if min(r[2] for r in L) or max(r[3] for r in L) != h:
            continue
        out.append(L)
    return out


def adjacency(n, E):
    adj = [set() for _ in range(n)]
    for a, b in E:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def star_centre(n, E):
    """If the graph is a star K_{1,m} with m >= 1, return (centre, m); else None."""
    if n < 2 or len(E) != n - 1:
        return None
    adj = adjacency(n, E)
    deg = [len(a) for a in adj]
    c = deg.index(max(deg))
    return (c, n - 1) if deg[c] == n - 1 else None


def max_independent_neighbourhood_degree(n, E):
    adj = adjacency(n, E)
    best = 0
    for v in range(n):
        N = sorted(adj[v])
        if all((min(a, b), max(a, b)) not in E for a, b in combinations(N, 2)):
            best = max(best, len(N))
    return best


def has_k4(n, E):
    for q in combinations(range(n), 4):
        if all((min(a, b), max(a, b)) in E for a, b in combinations(q, 2)):
            return q
    return None


def box(rects):
    x1 = min(r[0] for r in rects); x2 = max(r[1] for r in rects)
    y1 = min(r[2] for r in rects); y2 = max(r[3] for r in rects)
    w, h = x2 - x1, y2 - y1
    return (len(rects), h, w, w * h, 2 * (w + h))


# ---------------------------------------------------------------------------
# 5.  THE OBJECTS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
# paper.tex Section 4, display (2), "the published height-2 layout of K_{1,4}" (source line 370
# gives the row  K_{1,4} (S_1)  5 & 2 & 5 & 10 & 14).
K14_NAMES = ['a', 'b', 'c', 'd', 'e']
K14 = [(1, 2, 1, 2),        # a  TOP
       (3, 4, 1, 2),        # b  TOP
       (2, 3, 0, 2),        # c  FULL
       (0, 1, 0, 1),        # d  BOTTOM
       (4, 5, 0, 1)]        # e  BOTTOM

# paper.tex Section 4, display (4), "the published height-2 layout of G_5" (source line 417, row
# 7 & 2 & 5 & 10 & 14).
G5_NAMES = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
G5 = [(0, 1, 1, 2),         # a  TOP
      (0, 1, 0, 1),         # b  BOTTOM
      (1, 2, 0, 2),         # c  FULL
      (2, 3, 0, 2),         # d  FULL
      (3, 4, 0, 2),         # e  FULL
      (4, 5, 0, 1),         # f  BOTTOM
      (4, 5, 1, 2)]         # g  TOP

# paper.tex Section 3, Remark "one proof, not two": the layout that refutes the struck lemma.
TRI_NAMES = ['F1', 'b', 'F2']
TRI = [(0, 1, 0, 2),        # F1 FULL
       (1, 2, 0, 1),        # b  BOTTOM
       (2, 3, 0, 2)]        # F2 FULL

# paper.tex Section 4, Remark "the FULL-centre clause is needed only at m=4": K_{1,3} with a
# BOTTOM centre.
K13_NAMES = ['v', 'L', 'R', 't']
K13 = [(1, 2, 0, 1),        # v  BOTTOM, the centre
       (0, 1, 0, 1),        # L  BOTTOM
       (2, 3, 0, 1),        # R  BOTTOM
       (1, 2, 1, 2)]        # t  TOP

# paper.tex Section 4, display (5), last paragraph: a 7-vertex height-2 layout with TEN edges sharing the
# published bounding-box quintuple (7,2,5,10,14) of G_5, so a quintuple is not a graph identifier.
DECOY = [(0, 1, 0, 1), (0, 1, 1, 2), (1, 2, 0, 1), (1, 3, 1, 2),
         (2, 4, 0, 1), (3, 5, 1, 2), (4, 5, 0, 1)]

# paper.tex Section 1, Remark on integrality: a height-THREE layout whose graph is K_4.
K4_H3 = [(0, 2, 0, 1), (0, 2, 2, 3), (1, 2, 1, 2), (2, 3, 0, 3)]

# paper.tex Section 4, display (3): K_{1,5} realised in a 5 x 3 box, one row higher.
K15_NAMES = ['C', 'A', 'B', 'D', 'E', 'G']
K15_H3 = [(2, 3, 0, 2), (1, 2, 0, 1), (0, 1, 1, 2), (3, 4, 0, 1), (4, 5, 1, 2), (2, 3, 2, 3)]


def named(E, names):
    return sorted(tuple(sorted((names[a], names[b]))) for a, b in E)


def main():
    print('verify.py -- integer height-two rectangle visibility graphs are outerplanar')
    print('exact arithmetic only (int / fractions.Fraction); standard library only')
    print('')
    print('--- Part 1: the objects printed in the paper -------------------------------')

    # ---- K_{1,4} ----
    ok = interiors_disjoint(K14)
    e_def = vis_general(K14)
    e_col = vis_by_columns(K14, 5, 2)
    e_thm = theorem_a(K14)
    check('K14_layout_is_legal', ok, 'interiors pairwise disjoint')
    check('K14_definition_equals_theoremA', e_def == e_col == e_thm,
          'definition (Fraction midpoints) = definition (unit columns) = THEOREM A = %d edges'
          % len(e_thm))
    want = [('a', 'c'), ('b', 'c'), ('c', 'd'), ('c', 'e')]
    check('K14_graph_is_the_star_centred_at_the_FULL_rectangle',
          named(e_def, K14_NAMES) == want and K14[2][2:] == FULL,
          'E = %s, centre c = [2,3]x[0,2] is FULL' % (named(e_def, K14_NAMES),))
    check('K14_bounding_box_5_tuple_matches_source_line_370', box(K14) == (5, 2, 5, 10, 14),
          '(n,height,width,area,perimeter) = %s' % (box(K14),))
    check('K14_outerplanar_by_independent_cyclic_order_search',
          outerplanar_by_search(len(K14), e_def), 'brute force over all cyclic orders')

    # ---- G_5 ----
    e_def = vis_general(G5)
    e_thm = theorem_a(G5)
    check('G5_layout_is_legal', interiors_disjoint(G5), 'interiors pairwise disjoint')
    check('G5_definition_equals_theoremA',
          e_def == vis_by_columns(G5, 5, 2) == e_thm, '%d edges by all three routes' % len(e_thm))
    want = [('a', 'b'), ('a', 'c'), ('b', 'c'), ('c', 'd'), ('d', 'e'), ('e', 'f'), ('e', 'g'),
            ('f', 'g')]
    check('G5_is_two_triangles_joined_by_the_path_c_d_e', named(e_def, G5_NAMES) == want,
          '7 vertices, 8 edges: triangles abc and efg, path c-d-e')
    check('G5_bounding_box_5_tuple_matches_source_line_417', box(G5) == (7, 2, 5, 10, 14),
          '(n,height,width,area,perimeter) = %s' % (box(G5),))
    check('G5_outerplanar_by_independent_cyclic_order_search',
          outerplanar_by_search(len(G5), e_def), 'brute force over all cyclic orders')

    # ---- the struck lemma ----
    e_def = vis_general(TRI)
    bstar = [i for i, r in enumerate(TRI) if r[2] == 0]
    tri = has_k4(3, e_def) is None and len(e_def) == 3
    check('struck_lemma_premise_is_false_B_star_contains_a_triangle',
          tri and len(bstar) == 3 and e_def == theorem_a(TRI),
          'B* = {F1,b,F2} (all three), E = %s is a triangle, so a clique inside B* has size 3, '
          'not <= 2' % (named(e_def, TRI_NAMES),))

    # ---- K_{1,3} with a BOTTOM centre ----
    e_def = vis_general(K13)
    sc = star_centre(4, e_def)
    check('K13_is_realised_with_a_BOTTOM_centre',
          interiors_disjoint(K13) and e_def == theorem_a(K13) and sc == (0, 3)
          and K13[0][2:] == BOTTOM and box(K13)[1] == 2,
          'layout %s in a 3x2 box gives the star E = %s centred at the BOTTOM rectangle v, so the '
          '"FULL centre" clause of the star lemma holds at m=4 and not at m=3'
          % (K13, named(e_def, K13_NAMES)))

    # ---- a quintuple is not a graph identifier ----
    e_def = vis_general(DECOY)
    deg = [len(a) for a in adjacency(7, e_def)]
    check('bounding_box_quintuple_is_not_a_graph_identifier',
          interiors_disjoint(DECOY) and e_def == theorem_a(DECOY) and len(e_def) == 10
          and max(deg) == 4 and box(DECOY) == (7, 2, 5, 10, 14) == box(G5),
          'layout %s has the same quintuple (7,2,5,10,14) as G_5 but 10 edges, not 8, and a vertex '
          'of degree 4' % (DECOY,))

    print('')
    print('--- Part 2: exhaustive enumeration, integer layouts of bounding box [0,w]x[0,2] ----')
    W = 7
    total = 0
    global_max_indep = 0
    stars = {}                                           # m -> set of centre y-classes
    for w in range(1, W + 1):
        Ls = layouts(w, 2)
        total += len(Ls)
        bad_thm = bad_proof = 0
        worst_indep = 0
        for L in Ls:
            E = vis_by_columns(L, w, 2)
            if E != theorem_a(L):
                bad_thm += 1
                continue
            why = proof_goes_through(L, w, E)
            if why is not None:
                bad_proof += 1
                continue
            d = max_independent_neighbourhood_degree(len(L), E)
            if d > worst_indep:
                worst_indep = d
            sc = star_centre(len(L), E)
            if sc is not None:
                stars.setdefault(sc[1], set()).add(L[sc[0]][2:])
        global_max_indep = max(global_max_indep, worst_indep)
        check('structure_theorem_w%d' % w, bad_thm == 0,
              '%d layouts, THEOREM A = the definition on all of them' % len(Ls))
        check('outerplanarity_proof_w%d' % w, bad_proof == 0,
              '%d layouts: split, non-crossing two-layer drawing, contraction returns G' % len(Ls))
        check('independent_neighbourhood_degree_bound_w%d' % w, worst_indep <= 4,
              'max degree over vertices with an independent neighbourhood = %d' % worst_indep)

    print('')
    print('--- Part 3: the raw definition, cross-checked, and the corollary -------------------')
    dis = 0
    small = 0
    for w in range(1, 6):
        for L in layouts(w, 2):
            small += 1
            if vis_by_columns(L, w, 2) != vis_general(L):
                dis += 1
    check('integer_and_rational_definitions_agree', dis == 0,
          '%d layouts (w <= 5): unit-column test = Fraction-midpoint test' % small)

    check('K_1_5_is_realised_by_no_height_2_layout', 5 not in stars,
          'over all %d layouts (w <= %d) the star sizes realised are m in %s'
          % (total, W, sorted(stars)))
    check('K_1_4_is_realised_and_its_centre_is_always_FULL',
          4 in stars and stars[4] == {FULL},
          'centre y-classes for m=4: %s' % sorted(stars.get(4, ())))
    check('K_1_3_centre_need_not_be_FULL', stars.get(3) == {BOTTOM, TOP, FULL},
          'centre y-classes for m=3: %s -- so the certificate\'s "FULL centre" clause is needed '
          'only at m=4' % sorted(stars.get(3, ())))
    check('independent_neighbourhood_bound_is_exactly_4', global_max_indep == 4,
          'attained, and never exceeded, over all %d layouts' % total)

    print('')
    print('--- Part 4: anti-controls (a vacuous pass must be distinguishable) -----------------')
    k4 = set(combinations(range(4), 2))
    check('anticontrol_K4_is_rejected_by_the_outerplanarity_search',
          not outerplanar_by_search(4, k4), 'no cyclic order of K_4 avoids a crossing')
    k23 = set()
    for a in (0, 1):
        for b in (2, 3, 4):
            k23.add((a, b))
    check('anticontrol_K_2_3_is_rejected_by_the_outerplanarity_search',
          not outerplanar_by_search(5, k23), 'no cyclic order of K_{2,3} avoids a crossing')

    e3 = vis_general(K4_H3)
    q = has_k4(4, e3)
    check('anticontrol_height_3_admits_a_K4', interiors_disjoint(K4_H3) and q is not None
          and len(e3) == 6 and box(K4_H3)[1] == 3,
          'layout %s in a 3x3 box has visibility graph K_4, so the height-2 hypothesis is '
          'load-bearing and this machinery can see a non-outerplanar graph' % (K4_H3,))

    scale = Fraction(2, 3)
    resc = [(Fraction(x1), Fraction(x2), scale * y1, scale * y2) for (x1, x2, y1, y2) in K4_H3]
    h = max(r[3] for r in resc) - min(r[2] for r in resc)
    check('anticontrol_y_rescaling_gives_K4_bounding_box_height_2_in_the_real_model',
          h == 2 and vis_general(resc) == e3,
          'y -> (2/3)y maps that layout to bounding-box height exactly 2 with the SAME graph K_4, '
          'so without INTEGRALITY the theorem is false')

    e5 = vis_general(K15_H3)
    sc = star_centre(6, e5)
    check('anticontrol_height_3_admits_K_1_5', interiors_disjoint(K15_H3) and sc == (0, 5)
          and box(K15_H3)[1] == 3,
          'layout %s in a 5x3 box realises K_{1,5} centred at C, E = %s'
          % (K15_H3, named(e5, K15_NAMES)))

    bf = 0
    bfbad = 0
    for w in range(1, 5):
        for L in layouts(w, 2):
            bf += 1
            if not outerplanar_by_search(len(L), vis_by_columns(L, w, 2)):
                bfbad += 1
    check('independent_outerplanarity_route_w_le_4', bfbad == 0,
          '%d layouts: a non-crossing cyclic order was FOUND by brute force, without using the '
          'split/contract argument' % bf)

    print('')
    print('NOT RE-RUN: the enumeration is exhaustive for bounding boxes [0,w]x[0,2] with w <= 7 '
          '(711,202 layouts) and NOT beyond; the proof in the paper is what covers all w and all n, '
          'and this program checks that proof only on those widths.')
    print('NOT RE-RUN: the brute-force cyclic-order route runs on w <= 4 only (1,638 layouts); for '
          '5 <= w <= 7 outerplanarity is certified by the proof\'s own construction, not by an '
          'independent search.')
    print('NOT RE-RUN: nothing here touches height(G) for h >= 3 beyond the two exhibited '
          'anti-control layouts, and nothing here verifies that height(K_{1,5}) equals 3 (only '
          'that K_{1,5} is realisable at height 3 and at no height 2).')
    print('NOT RE-RUN: no literature claim is checked by this program. The prior-art status of the '
          'result -- including the reading of the 2014 GD poster abstract cited in the paper, which '
          'was obtained and searched by hand and not by this program -- is outside its scope.')
    print('')
    if FAILED:
        print('FAILURES: %s' % ', '.join(FAILED))
        print('VERDICT: %d of %d CHECKS FAILED' % (len(FAILED), len(CHECKS)))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % len(CHECKS))
    return 0


if __name__ == '__main__':
    sys.exit(main())

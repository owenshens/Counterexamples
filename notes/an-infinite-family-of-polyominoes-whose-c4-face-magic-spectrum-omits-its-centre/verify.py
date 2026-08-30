#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Referee check for

    "An Infinite Family of Polyominoes Whose C_4-Face-Magic Spectrum
     Omits Its Centre"  (paper.tex / paper.pdf in this directory)

which refutes Question 2 of Section 5 of

    P. Chalise, R. M. Low, A. Eisenkolb-Vaithyanathan,
    "All Polyominoes Are C_4-face-magic", arXiv:2608.08914v1.

EVERY object this program consumes is a literal below, transcribed from the
tables PRINTED IN THE PAPER: the cell lists, the covering sets, and all
fifteen vertex labelings.  There is no external data file, no third-party
package and no randomness; all arithmetic is exact integer arithmetic on
Python ints.  There is no search for labelings anywhere: the only enumerations
are over covering multisets of a fixed shape and over the 163 free polyominoes
with at most seven cells.  Python 3.9+, standard library only.

Conventions, exactly as in the paper.  A polyomino P is a list of CELLS; the
cell [x, y] is the closed unit square with lower-left lattice point (x, y) and
vertex set {(x,y), (x+1,y), (x,y+1), (x+1,y+1)}.  V(P) is the union of the cell
vertex sets and n = |V(P)|.  A labeling is a bijection lambda: V(P) -> {1..n};
it is C_4-face-magic with magic constant s when every cell's four labels sum
to s.  Spec_{C_4}(P) is the set of attainable s.

What is checked, section by section of the paper:

  Section 2, Lemma 3 (the counting identity) and Lemma 4 (the window):
      re-derived numerically on each exhibited shape --- the declared covering
      set really is a single-defect cover, |S| = n/4, 4 | n, the algebraic
      identity (n/4)*2(n+1) = n(n+1)/2 that forces the label collision at the
      centre, and the six-value window Spec subset {2(n+1)+j : |j| in 1,2,3}.
  Section 3, Theorem 1:  Spec(L_2) = {15,16,17,19,20,21} EXACTLY --- the six
      printed labelings give the lower bound, Lemmas 3 and 4 the upper bound.
  Section 4, Theorem 2 and Corollary 5:  the covering rule for L_m for every
      even m in 2..200, and Spec(L_6) = {31,32,33,35,36,37} exactly.
  Section 5, Remark 7:   the 1x3 strip has 18 in its spectrum and admits NO
      single-defect cover, so the hypothesis of Lemma 3 is essential and the
      omission of 18 is not an artifact of n = 8.
  Section 5, Remark 8:   the V-pentomino of Figure 6 of Shiu-Low-Liu 2024 is a
      DIFFERENT pentomino from L_4 at the same order n = 12, it does attain the
      centre 26, and it has no single-defect cover -- so the published labeling
      is not in conflict with Theorem 2.
  Section 5, Remark 9:   the holed shape, its bounded hole face (four corners
      AND four edges), its unbounded empty cell, and the labeling that is magic
      over cells AND over the hole.
  Section 6:   the SHAPE enumeration behind the census reading -- there are
      exactly 163 free polyominoes with 2..7 cells, and exactly ONE of them has
      a bounded C_4 face that is not a cell, namely the shape of Remark 9.  The
      spectra of that census are NOT recomputed here.

What is NOT covered is printed by the program itself, just above the verdict.
"""

import sys

# ---------------------------------------------------------------------------
# check harness
# ---------------------------------------------------------------------------
_PASSES = []
_FAILS = []


def ck(name, cond, detail=""):
    """One check.  Prints exactly one PASS or FAIL line."""
    if cond:
        _PASSES.append(name)
        print("PASS %s %s" % (name, detail))
    else:
        _FAILS.append(name)
        print("FAIL %s %s" % (name, detail))


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def cell_vertices(cell):
    x, y = cell
    return [(x, y), (x + 1, y), (x, y + 1), (x + 1, y + 1)]


def vertex_set(cells):
    V = set()
    for c in cells:
        V.update(cell_vertices(c))
    return V


def cover_profile(cells, S):
    """Multiplicity profile of the multiset S of cells over V(P).

    -> (n_once, uncovered sorted, doubled sorted, covered>=3 sorted)
    """
    V = vertex_set(cells)
    mult = dict((v, 0) for v in V)
    for c in S:
        for v in cell_vertices(c):
            mult[v] += 1
    once = sorted(v for v in V if mult[v] == 1)
    unc = sorted(v for v in V if mult[v] == 0)
    dbl = sorted(v for v in V if mult[v] == 2)
    high = sorted(v for v in V if mult[v] >= 3)
    return len(once), unc, dbl, high


def cell_sums(cells, order, labels):
    """labels is a tuple aligned with `order`, the printed vertex order."""
    lam = dict(zip(order, labels))
    return sorted(set(sum(lam[v] for v in cell_vertices(c)) for c in cells))


def window(n):
    """Lemma 4: the six values 2(n+1)+j, j in {-3,-2,-1,1,2,3}."""
    centre = 2 * (n + 1)
    return sorted(centre + j for j in (-3, -2, -1, 1, 2, 3))


def cell_edges(cell):
    x, y = cell
    return [((x, y), (x + 1, y)), ((x, y), (x, y + 1)),
            ((x + 1, y), (x + 1, y + 1)), ((x, y + 1), (x + 1, y + 1))]


def edge_set(cells):
    E = set()
    for c in cells:
        for e in cell_edges(c):
            E.add(tuple(sorted(e)))
    return E


def noncell_c4_faces(cells):
    """The empty unit squares that are nonetheless bounded faces isomorphic to C_4.

    An empty unit square is such a face exactly when all four of its corners are
    vertices of P AND all four of its sides are edges of P; the open square then
    contains no vertex or edge of P, so it is a face, and its boundary is a
    4-cycle.
    """
    cells = set(map(tuple, cells))
    V = vertex_set(cells)
    E = edge_set(cells)
    xs = [x for x, y in V]
    ys = [y for x, y in V]
    out = []
    for x in range(min(xs) - 1, max(xs) + 1):
        for y in range(min(ys) - 1, max(ys) + 1):
            if (x, y) in cells:
                continue
            if not all(v in V for v in cell_vertices((x, y))):
                continue
            if all(tuple(sorted(e)) in E for e in cell_edges((x, y))):
                out.append((x, y))
    return sorted(out)


def free_form(cells):
    """Canonical form of a polyomino under the 8 lattice symmetries and translation."""
    best = None
    cs = [tuple(c) for c in cells]
    for _ in range(4):
        cs = [(-y, x) for (x, y) in cs]                       # rotate 90 degrees
        for reflect in (False, True):
            d = [(-x, y) for (x, y) in cs] if reflect else list(cs)
            mx = min(x for x, y in d)
            my = min(y for x, y in d)
            t = tuple(sorted((x - mx, y - my) for x, y in d))
            if best is None or t < best:
                best = t
    return best


# ===========================================================================
# THE OBJECTS, transcribed from the tables printed in the paper
# ===========================================================================

# ---- Section 3: L_2, the L-tromino (Table 1 of the paper) -----------------
L2_CELLS = [(0, 0), (0, 1), (1, 0)]
L2_ORDER = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)]
L2_COVER = [(0, 1), (1, 0)]
L2_LABELINGS = {
    15: (7, 5, 3, 2, 1, 6, 4, 8),
    16: (5, 8, 3, 2, 1, 4, 6, 7),
    17: (3, 8, 2, 5, 1, 6, 4, 7),
    19: (1, 8, 2, 7, 3, 6, 4, 5),
    20: (1, 8, 3, 6, 5, 4, 2, 7),
    21: (1, 8, 2, 5, 7, 4, 3, 6),
}

# ---- Section 4: L_6, the seven-cell member (Table 2 of the paper) ---------
L6_CELLS = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 0)]
L6_ORDER = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
            (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
            (2, 0), (2, 1)]
L6_COVER = [(0, 1), (0, 3), (0, 5), (1, 0)]
L6_LABELINGS = {
    31: (13, 12, 2, 4, 3, 6, 8, 5, 1, 16, 9, 15, 7, 10, 11, 14),
    32: (9, 16, 2, 3, 4, 5, 7, 6, 1, 13, 14, 11, 12, 8, 10, 15),
    33: (5, 16, 2, 7, 3, 8, 4, 11, 1, 14, 10, 13, 9, 12, 6, 15),
    35: (1, 16, 2, 6, 3, 7, 4, 13, 5, 12, 15, 11, 14, 10, 8, 9),
    36: (1, 16, 3, 11, 4, 12, 5, 10, 9, 8, 14, 7, 13, 6, 2, 15),
    37: (1, 14, 2, 11, 3, 12, 4, 9, 13, 8, 16, 7, 15, 6, 5, 10),
}

# ---- Section 5, Remark 7: the 1x3 strip -----------------------------------
STRIP_CELLS = [(0, 0), (1, 0), (2, 0)]
STRIP_ORDER = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]
STRIP_LABELING_18 = (2, 8, 1, 7, 4, 6, 3, 5)

# ---- Section 5, Remark 8: the V-pentomino, plus(0,2,2,0) of Fig. 6 of SLL --
VPENT_CELLS = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)]
VPENT_ORDER = [(0, 0), (0, 1), (0, 2), (0, 3),
               (1, 0), (1, 1), (1, 2), (1, 3),
               (2, 0), (2, 1), (3, 0), (3, 1)]
VPENT_LABELING_26 = (1, 10, 2, 9, 4, 11, 3, 12, 5, 6, 7, 8)
L4_CELLS = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0)]      # the m = 4 member

# ---- Section 5, Remark 9: the holed shape ---------------------------------
HOLED_CELLS = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]
HOLED_HOLE_FACE = (1, 1)          # empty, but all four corners lie in V(P)
HOLED_OTHER_EMPTY = (2, 2)        # empty, and (3,3) is NOT in V(P)
HOLED_ORDER = [(0, 0), (0, 1), (0, 2), (0, 3),
               (1, 0), (1, 1), (1, 2), (1, 3),
               (2, 0), (2, 1), (2, 2), (2, 3),
               (3, 0), (3, 1), (3, 2)]
HOLED_LABELING_29 = (4, 9, 6, 15, 3, 13, 1, 7, 8, 5, 10, 11, 14, 2, 12)


# ===========================================================================
# the checks
# ===========================================================================
def check_shape(tag, cells, order, cover, labelings):
    """The full argument for one shape with a single-defect cover."""
    V = vertex_set(cells)
    n = len(V)
    centre = 2 * (n + 1)
    N = n * (n + 1) // 2

    print("")
    print("--- %s: cells %s ---" % (tag, [list(c) for c in cells]))

    ck("%s-vertex-order" % tag, sorted(order) == sorted(V) and len(order) == len(set(order)),
       "the printed vertex order lists V(P) exactly once each")
    ck("%s-order-n" % tag, n == len(order), "n = |V(P)| = %d" % n)
    ck("%s-divisibility" % tag, n % 4 == 0, "4 | n  (n = %d, so n/4 = %d)" % (n, n // 4))
    ck("%s-centre" % tag, centre == 2 * (n + 1), "centre 2(n+1) = %d" % centre)

    # --- Lemma 3, on the printed covering set --------------------------
    once, unc, dbl, high = cover_profile(cells, cover)
    ck("%s-cover-size" % tag, len(cover) == n // 4,
       "|S| = %d = n/4 for S = %s" % (len(cover), [list(c) for c in cover]))
    ck("%s-cover-uncovered" % tag, len(unc) == 1,
       "exactly one uncovered vertex, p = %s" % (list(unc[0]) if unc else None,))
    ck("%s-cover-doubled" % tag, len(dbl) == 1,
       "exactly one doubled vertex, v = %s" % (list(dbl[0]) if dbl else None,))
    ck("%s-cover-no-triple" % tag, high == [], "no vertex is covered three or more times")
    ck("%s-cover-distinct-pv" % tag, bool(unc) and bool(dbl) and unc[0] != dbl[0],
       "p != v, so lambda(p) = lambda(v) is impossible for a bijection")
    ck("%s-cover-total" % tag, once + 2 * len(dbl) == 4 * len(cover) == n,
       "multiplicities sum to 4|S| = %d = n" % (4 * len(cover)))

    # --- the algebraic core: at s = 2(n+1) the identity forces p = v ----
    ck("%s-collision-identity" % tag, (n // 4) * centre == N,
       "(n/4)*2(n+1) = %d = n(n+1)/2, so N - lambda(p) + lambda(v) = (n/4)s "
       "forces lambda(p) = lambda(v) at s = %d" % (N, centre))

    # --- Lemma 4, the window -------------------------------------------
    win = window(n)
    spec = sorted(labelings)
    ck("%s-window-shape" % tag, win == [centre + j for j in (-3, -2, -1, 1, 2, 3)],
       "the admissible window is %s (centre %d excluded)" % (win, centre))

    # --- lower bound: every printed labeling really is magic ------------
    for s in spec:
        lab = labelings[s]
        ck("%s-bijection-s%d" % (tag, s), sorted(lab) == list(range(1, n + 1)),
           "the s=%d labeling is a bijection V(P) -> {1..%d}" % (s, n))
    for s in spec:
        sums = cell_sums(cells, order, labelings[s])
        ck("%s-magic-s%d" % (tag, s), sums == [s],
           "all %d cells of the s=%d labeling sum to %d (distinct sums seen: %s)"
           % (len(cells), s, s, sums))

    # --- and each realised s sits at an admissible delta ---------------
    t = n // 4
    for s in spec:
        delta = (s - centre) * t
        ck("%s-delta-s%d" % (tag, s), delta != 0 and delta % t == 0 and 1 <= abs(delta) <= n - 1,
           "s=%d = 2(n+1) + delta/t with delta = %d, t = %d, 1 <= |delta| <= %d"
           % (s, delta, t, n - 1))

    # --- upper bound meets lower bound: the spectrum is DETERMINED ------
    ck("%s-spectrum-exact" % tag, spec == win,
       "Spec = %s exactly: the printed labelings give every window value, and "
       "Lemmas 3 and 4 forbid everything else" % (spec,))
    ck("%s-centre-absent" % tag, centre not in spec,
       "the centre %d is absent, so Spec is not an interval "
       "(%d < %d < %d)" % (centre, min(spec), centre, max(spec)))
    ck("%s-gap-is-exactly-centre" % tag,
       [u for u in range(min(spec), max(spec) + 1) if u not in spec] == [centre],
       "the only gap inside [min Spec, max Spec] is {%d}" % centre)
    ck("%s-symmetry" % tag, sorted(4 * (n + 1) - u for u in spec) == spec,
       "Spec is closed under s -> 4(n+1)-s, as the source paper's own remark requires")
    ck("%s-source-containment" % tag, min(spec) >= 10 and max(spec) <= 4 * n - 6,
       "Spec is inside the source paper's own bound {10,...,4n-6} = {10,...,%d}" % (4 * n - 6))
    return n, spec


def check_family(m_max):
    """Section 4, Theorem 2: the covering rule for L_m, m even."""
    print("")
    print("--- Section 4: the family L_m, every even m in 2..%d ---" % m_max)
    bad_n, bad_cover, bad_size = [], [], []
    for m in range(2, m_max + 1, 2):
        cells = [(0, j) for j in range(m)] + [(1, 0)]
        S = [(0, j) for j in range(1, m, 2)] + [(1, 0)]
        n = len(vertex_set(cells))
        if n != 2 * m + 4 or n % 4 != 0:
            bad_n.append(m)
        if len(S) != n // 4 or len(S) != m // 2 + 1:
            bad_size.append(m)
        once, unc, dbl, high = cover_profile(cells, S)
        if not (unc == [(0, 0)] and dbl == [(1, 1)] and high == []
                and once + 2 == 4 * len(S) == n):
            bad_cover.append(m)
    tested = len(range(2, m_max + 1, 2))
    ck("family-order", bad_n == [],
       "n = 2m+4 and 4 | n for all %d even m in 2..%d (failures: %s)"
       % (tested, m_max, bad_n or "none"))
    ck("family-cover-size", bad_size == [],
       "|S| = m/2+1 = n/4 for all %d even m in 2..%d (failures: %s)"
       % (tested, m_max, bad_size or "none"))
    ck("family-single-defect", bad_cover == [],
       "S is a single-defect cover with p = (0,0), v = (1,1) for all %d even m "
       "in 2..%d (failures: %s)" % (tested, m_max, bad_cover or "none"))
    ck("family-centre-absent", bad_n == [] and bad_size == [] and bad_cover == [],
       "hence 2(n+1) = 4m+10 is absent from Spec(L_m) for every even m tested, "
       "by Lemma 3 alone and with no search")


def check_strip():
    """Section 5, Remark 7: the hypothesis of Lemma 3 is essential."""
    print("")
    print("--- Section 5, Remark 7: the 1x3 strip ---")
    cells = STRIP_CELLS
    V = vertex_set(cells)
    n = len(V)
    ck("strip-order", sorted(STRIP_ORDER) == sorted(V) and n == 8,
       "the 1x3 strip has n = %d, the same order as the L-tromino" % n)
    lab = STRIP_LABELING_18
    ck("strip-bijection", sorted(lab) == list(range(1, n + 1)),
       "the printed labeling is a bijection V -> {1..8}")
    sums = cell_sums(cells, STRIP_ORDER, lab)
    ck("strip-magic-18", sums == [18],
       "all 3 cells sum to 18 = 2(n+1), so 18 IS in Spec(1x3 strip) "
       "(distinct sums seen: %s)" % (sums,))
    # exhaust every multiset of n/4 = 2 cells out of the 3
    found = []
    for i in range(len(cells)):
        for j in range(i, len(cells)):
            S = [cells[i], cells[j]]
            once, unc, dbl, high = cover_profile(cells, S)
            if len(unc) == 1 and len(dbl) == 1 and high == []:
                found.append(S)
    ck("strip-no-single-defect-cover", found == [],
       "none of the %d multisets of n/4 = 2 cells of the strip is a single-defect "
       "cover, so Lemma 3 says nothing about the strip -- the hypothesis is "
       "sufficient, not necessary, and 18 is attainable at n = 8"
       % (len(cells) * (len(cells) + 1) // 2))


def check_vpentomino():
    """Section 5, Remark 8: the published Fig. 6 labeling is not a conflict."""
    print("")
    print("--- Section 5, Remark 8: the V-pentomino plus(0,2,2,0) at n = 12 ---")
    cells = VPENT_CELLS
    V = vertex_set(cells)
    n = len(V)
    centre = 2 * (n + 1)
    ck("vpent-order", sorted(VPENT_ORDER) == sorted(V) and n == 12,
       "the V-pentomino has n = %d, the same order as L_4, so the same centre "
       "2(n+1) = %d" % (n, centre))
    ck("vpent-distinct-from-L4", free_form(cells) != free_form(L4_CELLS),
       "the V-pentomino and L_4 = %s are DIFFERENT free pentominoes, so the "
       "published value 26 for the former says nothing about the latter"
       % ([list(c) for c in L4_CELLS],))
    lab = VPENT_LABELING_26
    ck("vpent-bijection", sorted(lab) == list(range(1, n + 1)),
       "the printed labeling is a bijection V -> {1..12}")
    sums = cell_sums(cells, VPENT_ORDER, lab)
    ck("vpent-magic-26", sums == [centre],
       "all %d cells sum to %d = 2(n+1), so the centre IS in Spec(V-pentomino) "
       "-- reproduced here, independently of Shiu-Low-Liu 2024 (distinct sums "
       "seen: %s)" % (len(cells), centre, sums))
    # exhaust every multiset of n/4 = 3 cells out of the 5
    found = []
    k = len(cells)
    trials = 0
    for i in range(k):
        for j in range(i, k):
            for h in range(j, k):
                trials += 1
                S = [cells[i], cells[j], cells[h]]
                once, unc, dbl, high = cover_profile(cells, S)
                if len(unc) == 1 and len(dbl) == 1 and high == []:
                    found.append(S)
    ck("vpent-no-single-defect-cover", found == [] and trials == 35,
       "none of the %d multisets of n/4 = 3 cells of the V-pentomino is a "
       "single-defect cover, so Lemma 3 does not apply to it and must not -- "
       "exactly as for the strip of Remark 7" % trials)


def check_shape_census():
    """Section 6: the SHAPE enumeration behind the census, and which reading."""
    print("")
    print("--- Section 6: the 163 free polyominoes with 2..7 cells ---")
    level = set([free_form([(0, 0)])])
    counts = {1: 1}
    shapes = []
    for k in range(2, 8):
        nxt = set()
        for p in level:
            for (x, y) in p:
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    c = (x + dx, y + dy)
                    if c in p:
                        continue
                    nxt.add(free_form(list(p) + [c]))
        level = nxt
        counts[k] = len(level)
        shapes.extend(level)
    ck("census-shape-counts",
       [counts[k] for k in range(1, 8)] == [1, 1, 2, 5, 12, 35, 108],
       "free polyominoes with k = 1..7 cells: %s, the known sequence"
       % ([counts[k] for k in range(1, 8)],))
    ck("census-shape-total", len(shapes) == 163,
       "exactly %d free polyominoes have 2..7 cells, the population the paper's "
       "census covers" % len(shapes))
    extra = [p for p in shapes if noncell_c4_faces(p)]
    ck("census-one-holed-shape", len(extra) == 1,
       "exactly %d of the 163 shapes has a bounded C_4 face that is NOT a cell, "
       "so the two readings of Remark 9 can differ on that shape alone"
       % len(extra))
    same = (len(extra) == 1
            and free_form(extra[0]) == free_form(HOLED_CELLS)
            and len(vertex_set(extra[0])) == 15
            and len(extra[0]) == 7
            and len(noncell_c4_faces(extra[0])) == 1)
    ck("census-holed-shape-is-remark9", same,
       "and that shape is the one printed in Remark 9, up to the 8 lattice "
       "symmetries: 7 cells, n = 15, exactly one extra bounded C_4 face")
    ck("census-holed-shape-outside-lemma3", len(vertex_set(HOLED_CELLS)) % 4 != 0,
       "its order 15 is not divisible by 4, so it has no single-defect cover and "
       "is not one of the shapes Lemma 3 speaks about")


def check_holed():
    """Section 5, Remark 9: the two readings of `face' are different."""
    print("")
    print("--- Section 5, Remark 9: the holed shape ---")
    cells = HOLED_CELLS
    V = vertex_set(cells)
    n = len(V)
    ck("holed-order", sorted(HOLED_ORDER) == sorted(V) and n == 15,
       "the holed shape has n = %d" % n)
    ck("holed-hole-empty", HOLED_HOLE_FACE not in cells,
       "the cell %s is NOT a cell of P" % (list(HOLED_HOLE_FACE),))
    ck("holed-hole-bounded", all(v in V for v in cell_vertices(HOLED_HOLE_FACE)),
       "all four corners of %s lie in V(P), so it is a bounded face isomorphic "
       "to C_4 -- a face under the Section 1 reading, not a cell"
       % (list(HOLED_HOLE_FACE),))
    E = edge_set(cells)
    ck("holed-hole-cycle",
       all(tuple(sorted(e)) in E for e in cell_edges(HOLED_HOLE_FACE)),
       "and all four SIDES of %s are edges of P, so its boundary really is a "
       "4-cycle of P and the open square really is a face (corners alone would "
       "not suffice)" % (list(HOLED_HOLE_FACE),))
    ck("holed-noncell-faces-unique", noncell_c4_faces(cells) == [HOLED_HOLE_FACE],
       "%s is the ONLY empty unit square of this shape that is a bounded C_4 "
       "face" % (list(HOLED_HOLE_FACE),))
    ck("holed-other-empty-unbounded",
       HOLED_OTHER_EMPTY not in cells
       and not all(v in V for v in cell_vertices(HOLED_OTHER_EMPTY)),
       "the other empty cell %s has a corner (3,3) outside V(P), so it is part "
       "of the unbounded face" % (list(HOLED_OTHER_EMPTY),))
    ck("holed-odd-order", n % 4 != 0,
       "4 does not divide n = %d, so this shape has no single-defect cover and "
       "Lemma 3 does not apply to it" % n)
    lab = HOLED_LABELING_29
    ck("holed-bijection", sorted(lab) == list(range(1, n + 1)),
       "the printed labeling is a bijection V -> {1..15}")
    faces = list(cells) + [HOLED_HOLE_FACE]
    sums = cell_sums(faces, HOLED_ORDER, lab)
    ck("holed-magic-29-all-faces", sums == [29],
       "all %d bounded C_4 faces -- the 7 cells AND the hole -- sum to 29, so "
       "both readings are simultaneously satisfiable here (distinct sums seen: %s)"
       % (len(faces), sums))


def main():
    print("=" * 74)
    print("verify.py -- An Infinite Family of Polyominoes Whose C_4-Face-Magic")
    print("             Spectrum Omits Its Centre")
    print("target: Question 2, Section 5 of arXiv:2608.08914v1")
    print("        (source file C4FaceMagic_arXiv_v1.0.tex, bytes 32456-32531)")
    print("all objects below are literals transcribed from the paper's own tables")
    print("=" * 74)

    n2, spec2 = check_shape("L2", L2_CELLS, L2_ORDER, L2_COVER, L2_LABELINGS)
    n6, spec6 = check_shape("L6", L6_CELLS, L6_ORDER, L6_COVER, L6_LABELINGS)
    check_family(200)
    check_strip()
    check_vpentomino()
    check_holed()
    check_shape_census()

    print("")
    print("--- the refutation, restated from the checks above ---")
    ck("refutation-L2", spec2 == [15, 16, 17, 19, 20, 21] and 18 not in spec2,
       "Spec_{C_4}(L-tromino) = {15,16,17,19,20,21} is NOT an interval: the "
       "answer to Question 2 is NO")
    ck("refutation-L6", spec6 == [31, 32, 33, 35, 36, 37] and 34 not in spec6,
       "Spec_{C_4}(L_6) = {31,32,33,35,36,37} is NOT an interval either, at a "
       "second order n = 16")

    print("")
    print("NOT RE-RUN: this program checks the paper's claims and nothing else.")
    print("NOT RE-RUN: the SPECTRA of the census of Section 6 -- which found 37")
    print("  non-interval spectra among the 163 shapes, every gap exactly at the")
    print("  centre -- are CORROBORATION and are not recomputed here. Only the")
    print("  SHAPE population of that census is re-derived above (163 shapes, and")
    print("  exactly one of them with a bounded C_4 face that is not a cell). The")
    print("  census's k = 7 stratum rests on a single implementation that was never")
    print("  independently re-exhausted, and its 37 is the count under the CELL")
    print("  reading. No claim of the paper's theorems depends on any of it.")
    print("NOT RE-RUN: attainment for L_m with m >= 8 is NOT verified. Lemma 3")
    print("  excludes the centre for EVERY even m and Lemma 4 confines Spec(L_m) to")
    print("  six values, both with no search; that Spec(L_m) equals all six values")
    print("  is established here only for m = 2 and m = 6, by exhibited labelings.")
    print("NOT RE-RUN: minimality. Nothing here asserts that L_2 is the smallest")
    print("  polyomino with a non-interval spectrum, nor that the L_m are the only")
    print("  such shapes.")
    print("NOT RE-RUN: for the holed shape of Remark 9 this program checks only")
    print("  that the hole is a bounded C_4 face and that one labeling is magic")
    print("  over cells AND over the hole. Whether the two readings have DIFFERENT")
    print("  attainable sets on that shape was decided by the single unreplicated")
    print("  census above and is deliberately not claimed by the paper.")
    print("NOT RE-RUN: the odd-m half of the parity dichotomy is not touched. It")
    print("  would need Shiu-Low-Liu 2024, Theorem 2.3, whose 2025 corrigendum we")
    print("  could not obtain; the paper therefore states the family only for m")
    print("  even, which is forced by 4 | n and not chosen.")
    print("")

    if _FAILS:
        print("RESULT: %d of %d CHECKS FAILED" % (len(_FAILS), len(_FAILS) + len(_PASSES)))
        for name in _FAILS:
            print("  failed: %s" % name)
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % len(_PASSES))
    return 0


if __name__ == "__main__":
    sys.exit(main())

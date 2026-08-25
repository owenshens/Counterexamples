#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification program for "A z-Visibility Representation of K_{6,6} Minus a
Perfect Matching", which refutes Conjecture 6.5 of Bose et al.
(1994), restated as Sec. 5 of Bose et al. (1998).

Standard library only.  All decisions use exact integer or Fraction
arithmetic; no floating point is used anywhere.

=============================================================================
VALUES TAKEN FROM THE PAPER (INPUTS -- never checked against themselves)
=============================================================================
(I1) The exhibited object: the 12 rectangles of the paper's first table, as
     (label, z(R), x0, x1, y0, y1), the rectangle being [x0,x1]x[y0,y1]x{z}:
       h0/0/[6,9]x[0,7]   h3/1/[2,5]x[0,7]   v4/2/[4,8]x[1,3]
       v0/3/[2,4]x[1,6]   v5/4/[4,8]x[3,6]   h1/5/[1,8]x[2,5]
       v3/6/[2,10]x[2,5]  h2/7/[3,7]x[1,6]   v1/8/[3,7]x[1,7]
       h5/9/[2,8]x[1,3]   h4/10/[2,8]x[3,6]  v2/11/[0,8]x[0,4]
(I2) The paper's second table: for each of the 30 pairs h_i,v_j with i != j,
     the exhibited witness cell claimed to lie in F(h_i,v_j).
(I3) The paper's third table: for each of the remaining 36 pairs, the list of
     blocking rectangles claimed to lie in I(A,B) and to cover
     C(A) cap C(B), with the symbol "empty" for the cases where that
     intersection is itself empty.
(I4) The definitions of C(R), I(A,B), F(A,B) (paper, Sec. 2) and the
     definition of z-visibility (paper, Sec. 1, after Bose et al.).
(I5) The target graph G = K_{6,6} minus the perfect matching {h_i v_i}, which
     the program builds from its combinatorial description, not from a table.

Nothing else is taken from the paper.  In particular the numbers 30, 36, 66,
the visibility graph itself, all cell sets, all areas and all witness
cylinders are recomputed below and compared against the paper's assertions.

=============================================================================
WHAT IS DERIVED HERE (THE CHECKS)
=============================================================================
(D1) Well-formedness of the exhibited object: 12 rectangles, labels exactly
     {h0..h5, v0..v5}, integral coordinates, positive area, distinct heights
     0..11, pairwise disjointness in R^3 recomputed geometrically, and the
     bounding box of the configuration.
(D2) Cell sets C(R) computed twice (a brute-force scan over every candidate
     integer cell of the bounding box against the containment predicate
     Q_{a,b} subset pi(R), and a closed product formula) and compared.  This
     comparison is a consistency test of two routines that agree by a
     one-line argument whenever the boundaries are grid-aligned; it is NOT
     advertised as independent evidence.  What is genuinely tested, and can
     fail, is the hypothesis the paper's Lemma actually uses: that no unit
     cell is ever only PARTIALLY inside a projection.  That trichotomy is
     evaluated directly, on rational coordinates if the object supplies them,
     rather than inferred from the integrality of the table.
(D3) The target graph really is K_{6,6} minus a perfect matching -- and this
     is derived from the labels the EXHIBITED OBJECT carries, not from the
     module's constants: the two parts are read off the object's own label
     prefixes, K_{6,6} is built over those two parts, the deleted set is
     formed by pairing labels with a common numeric suffix and is verified to
     be a perfect matching, G is checked to have 30 edges, degree 5
     everywhere, and to be bipartite; only then is it compared with the
     constant description (I5).  Phrased over the constants these gates could
     not fail for any input at all.
(D4) LOAD-BEARING, method 1 (the paper's unit-cell criterion): for each of
     the 66 unordered pairs, F(A,B) is computed from scratch and the
     resulting visibility graph is compared to G.  Independently of that
     comparison, the graph the rectangles actually realise is shown to BE
     K_{6,6} minus a perfect matching using only itself -- bipartite between
     the two parts, 5-regular, 30 edges, and its complement inside K_{6,6}
     recomputed and shown to be a perfect matching -- so the refutation does
     not depend on the paper's choice of labelling.
(D5) LOAD-BEARING, method 2 (definition-level, independent of the paper's
     Lemma): for every pair the exact area of
     (pi(A) cap pi(B)) minus union of pi(R), R in I(A,B), is computed by
     coordinate compression over the true arrangement of all 24 rectangle
     coordinates in exact integer arithmetic.  Area 0 forbids any positive-
     radius disk, hence proves non-visibility; and for every claimed edge an
     explicit closed disk of radius 1/4 is verified in exact Fraction
     arithmetic to lie inside both endpoint projections and to be disjoint
     from every intermediate projection, which exhibits the cylinder the
     definition demands.  The graph so obtained is compared to G.
(D5') LOAD-BEARING, method 3: the same 66 exact areas recomputed by
     INCLUSION-EXCLUSION over the subsets of I(A,B), a routine that
     enumerates no cells and builds no arrangement and so shares no
     load-bearing logic with methods 1 or 2.  All three are required to agree.
     Two honest caveats are recorded rather than glossed: (i) for grid-aligned
     boundaries "F(A,B) nonempty" and "free area positive" are provably the
     same predicate, so the boolean agreement of methods 1 and 2 cannot fail
     and is not evidence -- what carries content is the numeric agreement,
     |F(A,B)| = area(A,B) exactly on all 66 pairs, and the agreement of the
     structurally unrelated method 3; (ii) the three methods are re-compared
     on all 162 mutant configurations of the census (D8), so an error in any
     one of them has 162 further chances to surface.
(D6) The paper's own two certificate tables: that they partition all 66
     pairs, that each of the 30 exhibited witness cells really lies in
     F(h_i,v_j), that every listed blocker really lies strictly between its
     pair in z, that the listed blockers really cover C(A) cap C(B), and
     that the "empty" symbol marks exactly the pairs with empty
     intersection.
(D7) The paper's arithmetic 30 + 36 = C(12,2) = 66, with the summand 30
     supplied by the unit-cell method and the summand 36 by the
     inclusion-exclusion method, so the identity is a cross-method check and
     not a tautology (30 + 30 or 36 + 36 would fail it).
(D8) Falsifiability self-tests: the object is deliberately corrupted in six
     ways (colliding heights, a deleted rectangle, a height moved above the
     stack, a projection shrunk by one unit, a bogus witness cell, and every
     single-blocker deletion from the paper's cover lists) and the
     corresponding check is required to report failure; plus an exhaustive
     mutation census over all 96 admissible single-endpoint +-1 perturbations
     of the twelve projections and all 66 transpositions of the height order.
     A MAJORITY of each family is required to destroy the representation --
     "at least one breaks it" is a gate almost nothing could fail -- and all
     three visibility methods are required to agree on each of the 162
     mutants.
(D9) The claim that the CONTAINMENT reading of "joins" is load-bearing, which
     the paper's Scope paragraph states as a fact.  It is computed here, not
     asserted: the graph the same twelve rectangles realise under the weaker
     reading is built from scratch and must strictly contain G, each
     additional pair must come with an explicit cylinder of exact positive
     radius, and the additional pairs are printed by name and counted from the
     computation.  The weakening evaluated is "the end disks need only MEET
     the two rectangles, the axis still on both", which is exactly
     nonemptiness -- rather than positive area -- of the CLOSED free set, and
     is therefore decidable exactly (see `touching_free_point`).

READING OF THE DEFINITION, and it is load-bearing.  "A closed vertical
cylinder joins A and B" is taken, as in Bose et al., to mean that the
cylinder's two end disks are CONTAINED in A and in B.  Under the weaker
reading in which the end disks need only MEET A and B, these twelve
rectangles realise a strictly larger graph and do NOT represent G, so this is
not a harmless normalisation -- and that is (D9), computed rather than taken
on trust.  The verification is, by contrast, insensitive to whether "meets
another rectangle" is read as closed or as interior intersection: the
exhibited witnesses keep clearance 1/4 from every intermediate projection,
and a free region of area 0 admits no positive-radius disk under either
reading.

NOT MACHINE-CHECKABLE HERE: that Conjecture 6.5 of Bose et al. (1994) is the
statement quoted by the paper (a bibliographic fact), and the general proof
of the paper's Lemma (a mathematical argument).  The Lemma's conclusion is
however independently confirmed on this instance by (D5)/(D5'), which never
use it.  In (D9) the step from the one weakening computed there to EVERY
weaker reading of "joins" is a one-line monotonicity remark -- weakening
"joins" only admits more cylinders -- and not a search over cylinders whose
axis touches neither rectangle.
"""

import itertools
import sys
from fractions import Fraction

_RESULTS = []


def check(name, ok, detail=""):
    """Record one falsifiable check and print it immediately."""
    ok = bool(ok)
    _RESULTS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    return ok


def info(msg):
    """Print a non-check informational line (never affects the verdict)."""
    print("INFO " + str(msg))


# ---------------------------------------------------------------------------
# (I1) INPUT: the exhibited object, transcribed from the paper's first table.
# Fields: label, z(R), x0, x1, y0, y1  ->  rectangle [x0,x1] x [y0,y1] x {z}.
# ---------------------------------------------------------------------------
PAPER_RECTS = [
    ("h0", 0, 6, 9, 0, 7),
    ("h3", 1, 2, 5, 0, 7),
    ("v4", 2, 4, 8, 1, 3),
    ("v0", 3, 2, 4, 1, 6),
    ("v5", 4, 4, 8, 3, 6),
    ("h1", 5, 1, 8, 2, 5),
    ("v3", 6, 2, 10, 2, 5),
    ("h2", 7, 3, 7, 1, 6),
    ("v1", 8, 3, 7, 1, 7),
    ("h5", 9, 2, 8, 1, 3),
    ("h4", 10, 2, 8, 3, 6),
    ("v2", 11, 0, 8, 0, 4),
]

# (I2) INPUT: the paper's witness table, one cell per edge h_i v_j (i != j).
PAPER_WITNESS = {
    ("h0", "v1"): (6, 6), ("h0", "v2"): (6, 0), ("h0", "v3"): (8, 2),
    ("h0", "v4"): (6, 1), ("h0", "v5"): (6, 3),
    ("h1", "v0"): (2, 2), ("h1", "v2"): (1, 2), ("h1", "v3"): (2, 2),
    ("h1", "v4"): (4, 2), ("h1", "v5"): (4, 3),
    ("h2", "v0"): (3, 1), ("h2", "v1"): (3, 1), ("h2", "v3"): (3, 2),
    ("h2", "v4"): (4, 1), ("h2", "v5"): (4, 5),
    ("h3", "v0"): (2, 1), ("h3", "v1"): (3, 6), ("h3", "v2"): (2, 0),
    ("h3", "v4"): (4, 1), ("h3", "v5"): (4, 3),
    ("h4", "v0"): (2, 5), ("h4", "v1"): (3, 3), ("h4", "v2"): (2, 3),
    ("h4", "v3"): (2, 3), ("h4", "v5"): (7, 5),
    ("h5", "v0"): (2, 1), ("h5", "v1"): (3, 1), ("h5", "v2"): (2, 1),
    ("h5", "v3"): (2, 2), ("h5", "v4"): (7, 1),
}

# (I3) INPUT: the paper's blocker table; [] transcribes the symbol "empty".
PAPER_BLOCKERS = {
    ("h0", "h1"): ["v4", "v5"], ("h0", "h2"): ["v4", "v5"], ("h0", "h3"): [],
    ("h0", "h4"): ["v5"], ("h0", "h5"): ["v4"], ("h1", "h2"): ["v3"],
    ("h1", "h3"): ["v0", "v4", "v5"], ("h1", "h4"): ["v3"],
    ("h1", "h5"): ["v3"],
    ("h2", "h3"): ["v0", "v4", "v5"], ("h2", "h4"): ["v1"],
    ("h2", "h5"): ["v1"],
    ("h3", "h4"): ["h2", "v0"], ("h3", "h5"): ["h2", "v0"],
    ("h4", "h5"): [],
    ("v0", "v1"): ["h2"], ("v0", "v2"): ["h1", "h5"], ("v0", "v3"): ["h1"],
    ("v0", "v4"): [], ("v0", "v5"): [], ("v1", "v2"): ["h4", "h5"],
    ("v1", "v3"): ["h2"], ("v1", "v4"): ["h2"], ("v1", "v5"): ["h2"],
    ("v2", "v3"): ["h4", "h5"], ("v2", "v4"): ["h5"], ("v2", "v5"): ["h1"],
    ("v3", "v4"): ["h1"], ("v3", "v5"): ["h1"], ("v4", "v5"): [],
    ("h0", "v0"): [], ("h1", "v1"): ["h2"], ("h2", "v2"): ["v1"],
    ("h3", "v3"): ["h1"], ("h4", "v4"): [], ("h5", "v5"): [],
}


# ---------------------------------------------------------------------------
# Configuration model.  A configuration is a dict so that the self-tests can
# corrupt it without touching the transcribed input.
# ---------------------------------------------------------------------------
def make_config(rects):
    cfg = {"labels": [r[0] for r in rects],
           "z": dict((r[0], r[1]) for r in rects),
           "box": dict((r[0], (r[2], r[3], r[4], r[5])) for r in rects)}
    return cfg


def clone(cfg):
    return {"labels": list(cfg["labels"]),
            "z": dict(cfg["z"]),
            "box": dict(cfg["box"])}


def bounding_box(cfg):
    xs0 = [cfg["box"][l][0] for l in cfg["labels"]]
    xs1 = [cfg["box"][l][1] for l in cfg["labels"]]
    ys0 = [cfg["box"][l][2] for l in cfg["labels"]]
    ys1 = [cfg["box"][l][3] for l in cfg["labels"]]
    return (min(xs0), max(xs1), min(ys0), max(ys1))


def _floor(v):
    """Exact floor for int or Fraction (no floating point)."""
    if isinstance(v, int):
        return v
    return v.numerator // v.denominator


def _ceil(v):
    """Exact ceiling for int or Fraction (no floating point)."""
    if isinstance(v, int):
        return v
    return -((-v.numerator) // v.denominator)


def partial_cell_violations(cfg, bb):
    """The hypothesis the paper's Lemma actually uses is not "the coordinates
    are integers" but its consequence: for every rectangle R and every integer
    cell Q_{a,b}, either Q_{a,b} subset pi(R) or Q_{a,b} cap pi(R) = empty --
    no cell is ever only PARTIALLY inside a projection.  That trichotomy is
    tested directly here, on rational coordinates if need be, so the test is
    falsifiable rather than inferred from integrality."""
    bx0, bx1, by0, by1 = bb
    bad = []
    for l in cfg["labels"]:
        x0, x1, y0, y1 = cfg["box"][l]
        for a in range(_floor(bx0), _ceil(bx1)):
            for b in range(_floor(by0), _ceil(by1)):
                inside = (x0 <= a and a + 1 <= x1 and y0 <= b and b + 1 <= y1)
                # (a,a+1)x(b,b+1) misses the CLOSED rectangle [x0,x1]x[y0,y1]
                disjoint = (a + 1 <= x0 or x1 <= a or b + 1 <= y0 or y1 <= b)
                if not inside and not disjoint:
                    bad.append("%s@(%s,%s)" % (l, a, b))
    return bad


def cells_by_scan(box, bb):
    """C(R) by testing the predicate Q_{a,b} subset pi(R) on every candidate
    integer cell of the whole configuration's bounding box."""
    x0, x1, y0, y1 = box
    bx0, bx1, by0, by1 = bb
    out = set()
    for a in range(_floor(bx0), _ceil(bx1)):
        for b in range(_floor(by0), _ceil(by1)):
            # Q_{a,b} = (a,a+1) x (b,b+1); its closure is contained in
            # [x0,x1] x [y0,y1] exactly when the four inequalities hold, and
            # for integral boundaries containment of the open cell is
            # equivalent to containment of its closure.
            if x0 <= a and a + 1 <= x1 and y0 <= b and b + 1 <= y1:
                out.add((a, b))
    return out


def cells_by_formula(box):
    """C(R) by the closed product formula, a second code path.  Total on
    rational as well as integral boxes: a runs over the integers with
    x0 <= a and a + 1 <= x1, i.e. ceil(x0) <= a <= floor(x1 - 1)."""
    x0, x1, y0, y1 = box
    if x1 <= x0 or y1 <= y0:
        return set()
    return set(itertools.product(
        range(_ceil(x0), _floor(x1 - 1) + 1),
        range(_ceil(y0), _floor(y1 - 1) + 1)))


def all_cell_sets(cfg):
    bb = bounding_box(cfg)
    return dict((l, cells_by_scan(cfg["box"][l], bb)) for l in cfg["labels"])


def strictly_between(cfg, a, b):
    """I(A,B) = {R : min(z(A),z(B)) < z(R) < max(z(A),z(B))}."""
    lo = min(cfg["z"][a], cfg["z"][b])
    hi = max(cfg["z"][a], cfg["z"][b])
    return [l for l in cfg["labels"] if lo < cfg["z"][l] < hi]


def free_cells(cfg, cells, a, b):
    """F(A,B) = (C(A) cap C(B)) minus union of C(R) over R in I(A,B)."""
    s = set(cells[a]) & set(cells[b])
    for r in strictly_between(cfg, a, b):
        s -= cells[r]
    return s


def all_pairs(cfg):
    return list(itertools.combinations(cfg["labels"], 2))


def visibility_graph_by_cells(cfg):
    """Method 1: the paper's unit-cell criterion."""
    cells = all_cell_sets(cfg)
    g = set()
    for a, b in all_pairs(cfg):
        if free_cells(cfg, cells, a, b):
            g.add(frozenset((a, b)))
    return g


def free_area(cfg, a, b):
    """Exact area of (pi(A) cap pi(B)) minus union of pi(R), R in I(A,B),
    by coordinate compression over the arrangement of ALL rectangle
    coordinates.  Integer arithmetic; makes no integrality assumption beyond
    using the coordinates that actually occur."""
    labels = cfg["labels"]
    box = cfg["box"]
    xs = sorted(set([box[l][0] for l in labels] + [box[l][1] for l in labels]))
    ys = sorted(set([box[l][2] for l in labels] + [box[l][3] for l in labels]))
    mid = strictly_between(cfg, a, b)

    def covers(l, cx0, cx1, cy0, cy1):
        p = box[l]
        return p[0] <= cx0 and cx1 <= p[1] and p[2] <= cy0 and cy1 <= p[3]

    total = 0
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            cx0, cx1, cy0, cy1 = xs[i], xs[i + 1], ys[j], ys[j + 1]
            if not (covers(a, cx0, cx1, cy0, cy1)
                    and covers(b, cx0, cx1, cy0, cy1)):
                continue
            if any(covers(m, cx0, cx1, cy0, cy1) for m in mid):
                continue
            total += (cx1 - cx0) * (cy1 - cy0)
    return total


def visibility_graph_by_area(cfg):
    """Method 2: an edge exactly when the exact free area is positive."""
    g = set()
    for a, b in all_pairs(cfg):
        if free_area(cfg, a, b) > 0:
            g.add(frozenset((a, b)))
    return g


def _rect_inter(r, s):
    if r is None or s is None:
        return None
    x0, x1 = max(r[0], s[0]), min(r[1], s[1])
    y0, y1 = max(r[2], s[2]), min(r[3], s[3])
    if x1 <= x0 or y1 <= y0:
        return None
    return (x0, x1, y0, y1)


def _rect_area(r):
    return 0 if r is None else (r[1] - r[0]) * (r[3] - r[2])


def free_area_incl_excl(cfg, a, b):
    """METHOD 3, structurally independent of methods 1 and 2: the same exact
    area, computed by INCLUSION-EXCLUSION over the subsets of I(A,B),

        |base \\ union pi(R)| = sum_{T subset I} (-1)^{|T|} |base cap cap_T pi|,

    where base = pi(A) cap pi(B).  This routine never enumerates unit cells and
    never builds an arrangement, so it shares no load-bearing logic with either
    of the other two methods; only the raw rectangle-intersection primitive."""
    base = _rect_inter(cfg["box"][a], cfg["box"][b])
    if base is None:
        return 0
    mid = strictly_between(cfg, a, b)
    total = _rect_area(base)
    for k in range(1, len(mid) + 1):
        s = 0
        for t in itertools.combinations(mid, k):
            r = base
            for x in t:
                r = _rect_inter(r, cfg["box"][x])
                if r is None:
                    break
            s += _rect_area(r)
        total += (-1) ** k * s
    return total


def visibility_graph_by_incl_excl(cfg):
    g = set()
    for a, b in all_pairs(cfg):
        if free_area_incl_excl(cfg, a, b) > 0:
            g.add(frozenset((a, b)))
    return g


def disk_witness_ok(cfg, a, b, cell):
    """Definition-level witness.  Verify in exact Fraction arithmetic that the
    closed disk of radius 1/4 about the centre of unit cell `cell` lies inside
    pi(A) and pi(B) and misses pi(R) for every other rectangle R whose height
    lies in the closed interval [min z, max z] spanned by the cylinder.  Such
    a disk, extruded from z(A) to z(B), is a closed vertical cylinder of
    positive radius and positive length joining A and B and meeting no other
    rectangle -- exactly the object the definition of z-visibility asks for."""
    cx = Fraction(2 * cell[0] + 1, 2)
    cy = Fraction(2 * cell[1] + 1, 2)
    rad = Fraction(1, 4)
    za, zb = cfg["z"][a], cfg["z"][b]
    if za == zb:
        return False, "cylinder length is zero (z(%s) == z(%s))" % (a, b)
    lo, hi = min(za, zb), max(za, zb)
    for l in cfg["labels"]:
        x0, x1, y0, y1 = cfg["box"][l]
        if l == a or l == b:
            if not (x0 <= cx - rad and cx + rad <= x1
                    and y0 <= cy - rad and cy + rad <= y1):
                return False, "disk not inside pi(%s)" % l
            continue
        if not (lo <= cfg["z"][l] <= hi):
            continue
        dx = max(Fraction(x0) - cx, Fraction(0), cx - Fraction(x1))
        dy = max(Fraction(y0) - cy, Fraction(0), cy - Fraction(y1))
        if dx * dx + dy * dy <= rad * rad:
            return False, "disk meets pi(%s)" % l
    return True, ""


def touching_free_point(cfg, a, b):
    """THE WEAKER READING, made computable.  Under the containment reading a
    cylinder joining A and B needs a whole disk of positive radius inside
    pi(A) and inside pi(B), so the admissible centres form a set of POSITIVE
    AREA.  Weaken "the end disks are contained in A and B" to "the end disks
    MEET A and B" while keeping the axis on the two rectangles, and the
    requirement collapses to nonemptiness of the CLOSED free set

        Fbar(A,B) = (pi(A) cap pi(B)) minus union of pi(R), R in I(A,B),

    the closed intersection being used and no area being demanded: from any
    point of Fbar(A,B) a disk of small enough positive radius meets pi(A) and
    pi(B) (it contains that point, which lies in both) and misses every
    intermediate projection, those being finitely many closed sets that avoid
    the point.  Conversely a cylinder whose axis meets both rectangles and
    misses the intermediate ones has its axis in Fbar(A,B).

    Decided exactly, with no search over radii: every rectangle boundary
    coordinate that occurs anywhere is a cut, membership in each closed
    rectangle is constant on each face of the resulting arrangement, and every
    face meeting the closed base contains one of the representative points
    below -- a cut value, or the exact rational midpoint of two consecutive cut
    values.  Returns such a point of Fbar(A,B) as a pair of Fractions, or None
    when Fbar(A,B) is empty."""
    box = cfg["box"]
    ax0, ax1, ay0, ay1 = box[a]
    bx0, bx1, by0, by1 = box[b]
    x_lo, x_hi = max(ax0, bx0), min(ax1, bx1)
    y_lo, y_hi = max(ay0, by0), min(ay1, by1)
    if x_hi < x_lo or y_hi < y_lo:
        return None
    mid = strictly_between(cfg, a, b)

    def reps(lo, hi, cuts):
        vals = sorted(set([lo, hi] + [c for c in cuts if lo < c < hi]))
        out = []
        for i, v in enumerate(vals):
            out.append(Fraction(v))
            if i + 1 < len(vals):
                out.append(Fraction(v + vals[i + 1], 2))
        return out

    xcuts = [box[l][i] for l in cfg["labels"] for i in (0, 1)]
    ycuts = [box[l][i] for l in cfg["labels"] for i in (2, 3)]
    for px in reps(x_lo, x_hi, xcuts):
        for py in reps(y_lo, y_hi, ycuts):
            blocked = False
            for m in mid:
                mx0, mx1, my0, my1 = box[m]
                if mx0 <= px <= mx1 and my0 <= py <= my1:
                    blocked = True
                    break
            if not blocked:
                return (px, py)
    return None


def touching_visibility_graph(cfg):
    """The graph the rectangles realise under the weaker reading above."""
    g = set()
    for a, b in all_pairs(cfg):
        if touching_free_point(cfg, a, b) is not None:
            g.add(frozenset((a, b)))
    return g


def touching_disk_witness(cfg, a, b, pt):
    """Exhibit, in exact Fraction arithmetic, a positive radius r for which the
    closed disk of radius r about `pt` meets pi(A) and pi(B) and is disjoint
    from pi(R) for every R in I(A,B) -- i.e. an actual cylinder of positive
    radius and positive length joining A and B under the weaker reading.
    Returns (ok, r squared, reason)."""
    cx, cy = pt
    if cfg["z"][a] == cfg["z"][b]:
        return False, None, "cylinder length is zero (z(%s) == z(%s))" % (a, b)
    for l in (a, b):
        x0, x1, y0, y1 = cfg["box"][l]
        if not (x0 <= cx <= x1 and y0 <= cy <= y1):
            return False, None, "axis point not on pi(%s)" % l
    d2s = []
    for l in strictly_between(cfg, a, b):
        x0, x1, y0, y1 = cfg["box"][l]
        dx = max(Fraction(x0) - cx, Fraction(0), cx - Fraction(x1))
        dy = max(Fraction(y0) - cy, Fraction(0), cy - Fraction(y1))
        d2 = dx * dx + dy * dy
        if d2 == 0:
            return False, None, "axis point lies on pi(%s)" % l
        d2s.append((d2, l))
    # r = (distance to the nearest intermediate projection)/2, so r > 0 and
    # r^2 < d^2 for every intermediate; with none, any positive radius serves.
    r2 = Fraction(1, 4) if not d2s else min(d for d, _ in d2s) / 4
    if r2 <= 0:
        return False, None, "radius not positive"
    for d2, l in d2s:
        if d2 <= r2:
            return False, None, "disk meets pi(%s)" % l
    return True, r2, ""


def rectangles_disjoint(cfg, a, b):
    """Exact disjointness of two closed rectangles in R^3."""
    if cfg["z"][a] != cfg["z"][b]:
        return True
    p, q = cfg["box"][a], cfg["box"][b]
    return p[1] < q[0] or q[1] < p[0] or p[3] < q[2] or q[3] < p[2]


H_PART = ["h%d" % i for i in range(6)]
V_PART = ["v%d" % i for i in range(6)]
MATCHING = [frozenset(("h%d" % i, "v%d" % i)) for i in range(6)]


def complete_bipartite():
    return set(frozenset((h, v)) for h in H_PART for v in V_PART)


def target_graph():
    """(I5) G = K_{6,6} minus the perfect matching {h_i v_i}, built from the
    combinatorial description in the paper's Sec. 1."""
    return complete_bipartite() - set(MATCHING)


def show_pair(e):
    # Total, including on the one-element "pair" a duplicated label produces:
    # that case must be reported by the check that detects it, not raised as an
    # unpacking error before the verdict is printed.
    return "-".join(sorted(e))


def check_well_formed(cfg):
    labels = cfg["labels"]
    ok_n = check("object_decoded_12_rectangles", len(labels) == 12,
                 "n=%d" % len(labels))
    ok_lab = check("labels_are_exactly_h0_h5_v0_v5",
                   sorted(labels) == sorted(H_PART + V_PART),
                   "labels=%s" % ",".join(sorted(labels)))
    ok_d = check("labels_distinct", len(set(labels)) == len(labels))
    ints = all(isinstance(v, int)
               for l in labels
               for v in (cfg["z"][l],) + tuple(cfg["box"][l]))
    ok_int = check("all_coordinates_integral_so_boundaries_lie_on_grid_lines",
                   ints)
    nondeg = [l for l in labels
              if cfg["box"][l][1] <= cfg["box"][l][0]
              or cfg["box"][l][3] <= cfg["box"][l][2]]
    ok_area = check("every_projection_has_positive_area", not nondeg,
                    "degenerate=%s" % (",".join(nondeg) if nondeg else "none"))
    # The hypothesis the Lemma really needs, tested rather than inferred.
    pv = partial_cell_violations(cfg, bounding_box(cfg)) if not nondeg else []
    check("no_unit_cell_is_only_partially_inside_any_projection", not pv,
          "partial=%s" % (",".join(pv[:6]) + ("..." if len(pv) > 6 else "")
                          if pv else "none"))
    zs = [cfg["z"][l] for l in labels]
    check("heights_pairwise_distinct", len(set(zs)) == len(zs),
          "n_distinct=%d" % len(set(zs)))
    check("heights_are_the_consecutive_integers_0_to_11",
          sorted(zs) == list(range(12)), "z=%s" % sorted(zs))
    bad = [show_pair(frozenset((a, b))) for a, b in all_pairs(cfg)
           if not rectangles_disjoint(cfg, a, b)]
    check("rectangles_pairwise_disjoint_in_R3", not bad,
          "overlapping=%s" % (",".join(bad) if bad else "none"))
    bb = bounding_box(cfg)
    check("bounding_box_is_0_10_x_0_7", bb == (0, 10, 0, 7),
          "bbox=[%s,%s]x[%s,%s]" % bb)
    # The ONLY structural precondition that makes the rest of the program
    # undefined rather than merely failing is the label set: every downstream
    # routine indexes the configuration by the twelve names h0..h5, v0..v5.
    # Non-integral or degenerate boxes are deliberately NOT treated as fatal --
    # the cell routines are total on rationals and on empty boxes, so those
    # objects are carried all the way through and made to FAIL the geometric
    # checks out loud instead of short-circuiting them.
    del ok_int, ok_area
    usable = bool(ok_n and ok_lab and ok_d)
    return bb, usable


def print_object(cfg, cells):
    """Print the decoded object back, in increasing z, as the paper lists it."""
    info("decoded object, in increasing z(R):")
    for l in sorted(cfg["labels"], key=lambda t: cfg["z"][t]):
        x0, x1, y0, y1 = cfg["box"][l]
        info("  %-2s z=%-2d pi=[%d,%d]x[%d,%d]  |C(R)|=%d"
             % (l, cfg["z"][l], x0, x1, y0, y1, len(cells[l])))


def check_cell_sets(cfg, bb):
    """(D2) C(R) computed two independent ways, then compared."""
    scan = dict((l, cells_by_scan(cfg["box"][l], bb)) for l in cfg["labels"])
    form = dict((l, cells_by_formula(cfg["box"][l])) for l in cfg["labels"])
    mismatch = [l for l in cfg["labels"] if scan[l] != form[l]]
    check("cell_sets_agree_between_scan_and_product_formula", not mismatch,
          "mismatch=%s" % (",".join(mismatch) if mismatch else "none"))
    empty = [l for l in cfg["labels"] if not scan[l]]
    check("every_cell_set_C_R_is_nonempty", not empty,
          "empty=%s" % (",".join(empty) if empty else "none"))
    total = sum(len(scan[l]) for l in cfg["labels"])
    info("total cells over the 12 rectangles: %d" % total)
    return scan


def check_target_graph(cfg):
    """(D3) The graph being represented really is K_{6,6} minus a perfect
    matching.

    Everything here is built from the labels the EXHIBITED OBJECT carries, not
    from the module-level constants: the two parts are read off the object's
    label prefixes, K_{6,6} is built over those two parts, and the deleted set
    is formed by pairing labels that share a numeric suffix.  Only at the end
    is the result compared with the constant description (I5).  Were these
    gates phrased over the hardcoded constants instead they could not fail for
    any input, which is precisely what this rewrite avoids."""
    hs = sorted(set(l for l in cfg["labels"]
                    if len(l) == 2 and l[0] == "h" and l[1].isdigit()))
    vs = sorted(set(l for l in cfg["labels"]
                    if len(l) == 2 and l[0] == "v" and l[1].isdigit()))
    kb = set(frozenset((h, v)) for h in hs for v in vs)
    mset = set(frozenset(("h" + l[1], "v" + l[1])) for l in hs
               if ("v" + l[1]) in vs)
    g = kb - mset
    check("the_exhibited_labels_split_into_two_parts_of_size_6_and_6",
          len(hs) == 6 and len(vs) == 6
          and len(hs) + len(vs) == len(set(cfg["labels"]))
          and len(cfg["labels"]) == 12,
          "|H|=%d |V|=%d of %d labels" % (len(hs), len(vs),
                                          len(cfg["labels"])))
    check("K_6_6_has_36_edges", len(kb) == 36, "|E|=%d" % len(kb))
    covered = set()
    for e in mset:
        covered |= set(e)
    check("deleted_set_is_a_perfect_matching_of_K_6_6",
          len(mset) == 6 and mset <= kb
          and covered == set(hs) | set(vs) and len(covered) == 12,
          "|M|=%d covers %d vertices" % (len(mset), len(covered)))
    check("target_graph_has_30_edges", len(g) == 30, "|E(G)|=%d" % len(g))
    deg = {}
    for e in g:
        for x in e:
            deg[x] = deg.get(x, 0) + 1
    check("every_vertex_of_G_has_degree_5",
          len(deg) == 12 and set(deg.values()) == set([5]),
          "degrees=%s" % sorted(set(deg.values())))
    check("G_is_bipartite_between_H_and_V",
          bool(g) and all(len(set(e) & set(hs)) == 1 for e in g))
    check("target_graph_from_the_objects_labels_equals_the_stated_description",
          g == target_graph() and kb == complete_bipartite()
          and mset == set(MATCHING),
          "|E|=%d vs %d" % (len(g), len(target_graph())))
    return g


def check_represented_graph_shape(cfg, vis):
    """The refutation only needs the represented graph to BE K_{6,6} minus a
    perfect matching.  Everything here is derived from the computed visibility
    graph `vis` alone; the paper's labelling of the deleted matching is not
    used, so this survives any relabelling."""
    hs = set(l for l in cfg["labels"] if l.startswith("h"))
    vs = set(l for l in cfg["labels"] if l.startswith("v"))
    cross = [show_pair(e) for e in vis
             if len(set(e) & hs) != 1 or len(set(e) & vs) != 1]
    check("represented_graph_is_bipartite_between_the_two_parts", not cross,
          "non_crossing_edges=%s" % (",".join(sorted(cross)) if cross
                                     else "none"))
    deg = dict((l, 0) for l in cfg["labels"])
    for e in vis:
        for x in e:
            deg[x] = deg.get(x, 0) + 1
    check("represented_graph_is_5_regular_on_all_12_vertices",
          set(deg.values()) == set([5]) and len(deg) == 12,
          "degrees=%s" % sorted(set(deg.values())))
    comp = set(frozenset((h, v)) for h in hs for v in vs) - vis
    covered = set()
    for e in comp:
        covered |= set(e)
    check("complement_of_represented_graph_inside_K_6_6_is_a_perfect_matching",
          len(comp) == 6 and covered == hs | vs and len(covered) == 12,
          "|complement|=%d covers %d vertices" % (len(comp), len(covered)))
    info("deleted matching recovered from the computed graph: %s"
         % ", ".join(sorted(show_pair(e) for e in comp)))
    check("represented_graph_has_30_edges_so_it_is_K66_minus_a_matching",
          len(vis) == 30, "|E|=%d" % len(vis))


def check_load_bearing_cells(cfg, cells, g):
    """(D4) LOAD-BEARING method 1: recompute F(A,B) for all 66 pairs and
    compare the resulting visibility graph with G."""
    pairs = all_pairs(cfg)
    n = len(cfg["labels"])
    check("enumerated_pair_count_equals_binom_n_2_and_equals_66",
          len(pairs) == n * (n - 1) // 2 and len(pairs) == 66,
          "n=%d pairs=%d binom=%d" % (n, len(pairs), n * (n - 1) // 2))
    vis = set()
    fsize = {}
    for a, b in pairs:
        f = free_cells(cfg, cells, a, b)
        fsize[frozenset((a, b))] = len(f)
        if f:
            vis.add(frozenset((a, b)))
    extra = sorted(show_pair(e) for e in vis - g)
    missing = sorted(show_pair(e) for e in g - vis)
    check("visibility_graph_by_unit_cell_criterion_equals_K66_minus_matching",
          vis == g,
          "extra=%s missing=%s"
          % (",".join(extra) if extra else "none",
             ",".join(missing) if missing else "none"))
    check("recomputed_edge_count_is_30", len(vis) == 30, "|E|=%d" % len(vis))
    nonedges = sum(1 for e in fsize if fsize[e] == 0)
    check("recomputed_nonedge_count_is_36", nonedges == 36,
          "nonedges=%d" % nonedges)
    nz = [show_pair(e) for e in g if fsize[e] == 0]
    check("every_edge_of_G_has_a_nonempty_free_cell_set", not nz,
          "empty_F=%s" % (",".join(sorted(nz)) if nz else "none"))
    nzz = [show_pair(e) for e in
           (set(frozenset(p) for p in pairs) - g) if fsize[e] > 0]
    check("every_nonedge_of_G_has_an_empty_free_cell_set", not nzz,
          "nonempty_F=%s" % (",".join(sorted(nzz)) if nzz else "none"))
    info("min |F(A,B)| over the 30 edges: %d"
         % min(fsize[e] for e in g))
    return vis, fsize


def check_load_bearing_geometry(cfg, cells, g, fsize):
    """(D5) LOAD-BEARING method 2, independent of the paper's Lemma: exact
    areas, and an explicit rational cylinder for every claimed edge."""
    pairs = all_pairs(cfg)
    areas = dict((frozenset(p), free_area(cfg, p[0], p[1])) for p in pairs)
    vis_area = set(e for e in areas if areas[e] > 0)

    extra = sorted(show_pair(e) for e in vis_area - g)
    missing = sorted(show_pair(e) for e in g - vis_area)
    check("visibility_graph_by_exact_free_area_equals_K66_minus_matching",
          vis_area == g,
          "extra=%s missing=%s"
          % (",".join(extra) if extra else "none",
             ",".join(missing) if missing else "none"))

    disagree = sorted(show_pair(e) for e in areas
                      if (areas[e] > 0) != (fsize[e] > 0))
    check("unit_cell_criterion_and_exact_area_agree_on_all_66_pairs",
          not disagree,
          "disagree=%s" % (",".join(disagree) if disagree else "none"))

    # A far tighter comparison than the boolean one above: with the projection
    # boundaries on grid lines the free region is a union of unit cells, so the
    # COUNT of free cells must equal the AREA exactly, pair by pair.  Any
    # off-by-one in the arrangement, in `covers`, or in the cell predicate
    # shows up here even when the sign of the area is unaffected.
    numdis = sorted("%s:%d!=%d" % (show_pair(e), fsize[e], areas[e])
                    for e in areas if fsize[e] != areas[e])
    check("free_cell_count_equals_free_area_exactly_on_all_66_pairs",
          not numdis,
          "mismatch=%s" % (",".join(numdis) if numdis else "none"))

    # METHOD 3: the same areas by inclusion-exclusion over subsets of I(A,B).
    # This shares no logic with methods 1 or 2 beyond rectangle intersection,
    # so agreement here is an implementation cross-check with real content,
    # unlike the comparison of two routines that both enumerate a grid.
    ie = dict((frozenset(p), free_area_incl_excl(cfg, p[0], p[1]))
              for p in pairs)
    iedis = sorted("%s:%d!=%d" % (show_pair(e), ie[e], areas[e])
                   for e in areas if ie[e] != areas[e])
    check("inclusion_exclusion_area_matches_coordinate_compression_on_66_pairs",
          not iedis,
          "mismatch=%s" % (",".join(iedis) if iedis else "none"))
    vis_ie = set(e for e in ie if ie[e] > 0)
    check("visibility_graph_by_inclusion_exclusion_equals_K66_minus_matching",
          vis_ie == g,
          "extra=%s missing=%s"
          % (",".join(sorted(show_pair(e) for e in vis_ie - g)) or "none",
             ",".join(sorted(show_pair(e) for e in g - vis_ie)) or "none"))

    zbad = [show_pair(frozenset(p)) for p in pairs
            if cfg["z"][p[0]] == cfg["z"][p[1]]]
    check("every_pair_admits_a_cylinder_of_positive_length", not zbad,
          "equal_heights=%s" % (",".join(zbad) if zbad else "none"))

    # Sufficiency, at the level of the definition: exhibit the cylinder.
    bad = []
    for e in sorted(g, key=show_pair):
        a, b = sorted(e)
        f = free_cells(cfg, cells, a, b)
        if not f:
            bad.append(show_pair(e) + ":no free cell")
            continue
        cell = sorted(f)[0]
        ok, why = disk_witness_ok(cfg, a, b, cell)
        if not ok:
            bad.append(show_pair(e) + ":" + why)
    check("explicit_radius_one_quarter_cylinder_verified_for_all_30_edges",
          not bad and len(g) == 30,
          "failures=%s" % (";".join(bad) if bad else "none"))

    # Necessity, at the level of the definition: zero area forbids any disk.
    nonzero = sorted(show_pair(e) for e in areas if e not in g
                     and areas[e] != 0)
    check("all_36_nonadjacent_pairs_have_exactly_zero_free_area",
          not nonzero and len([e for e in areas if e not in g]) == 36,
          "positive_area=%s" % (",".join(nonzero) if nonzero else "none"))
    # (D7) the paper's arithmetic 30 + 36 = C(12,2), with the two summands
    # produced by the two independent methods so the identity is not circular:
    # 30 from the unit-cell criterion, 36 from the exact-area computation.
    edges_by_cells = sum(1 for e in fsize if fsize[e] > 0)
    nonedges_by_area = sum(1 for e in ie if ie[e] == 0)
    n = len(cfg["labels"])
    check("paper_arithmetic_30_by_cells_plus_36_by_area_equals_binom_12_2",
          edges_by_cells == 30 and nonedges_by_area == 36
          and sum(1 for e in areas if areas[e] == 0) == 36
          and edges_by_cells + nonedges_by_area == n * (n - 1) // 2 == 66,
          "%d+%d=%d, binom(%d,2)=%d"
          % (edges_by_cells, nonedges_by_area,
             edges_by_cells + nonedges_by_area, n, n * (n - 1) // 2))
    info("free area over the 30 edges: min=%d max=%d"
         % (min(areas[e] for e in g), max(areas[e] for e in g)))
    return areas


def check_witness_table(cfg, cells, g, witness):
    """(D6a) Every cell in the paper's second table really lies in F."""
    keys = set(frozenset(k) for k in witness)
    check("witness_table_has_30_entries", len(witness) == 30
          and len(keys) == 30, "entries=%d distinct=%d"
          % (len(witness), len(keys)))
    check("witness_table_keys_are_exactly_the_edges_of_G", keys == g,
          "extra=%s missing=%s"
          % (",".join(sorted(show_pair(e) for e in keys - g)) or "none",
             ",".join(sorted(show_pair(e) for e in g - keys)) or "none"))
    bad = []
    for (a, b), cell in sorted(witness.items()):
        if a not in cfg["z"] or b not in cfg["z"]:
            bad.append("%s-%s[not a rectangle of the object]" % (a, b))
            continue
        f = free_cells(cfg, cells, a, b)
        if cell not in f:
            reason = []
            if cell not in cells.get(a, set()):
                reason.append("not in C(%s)" % a)
            if cell not in cells.get(b, set()):
                reason.append("not in C(%s)" % b)
            for r in strictly_between(cfg, a, b):
                if cell in cells[r]:
                    reason.append("blocked by %s" % r)
            bad.append("%s-%s(%d,%d)%s"
                       % (a, b, cell[0], cell[1],
                          "[" + ",".join(reason) + "]" if reason else ""))
    check("every_paper_witness_cell_lies_in_F_of_its_pair", not bad,
          "bad=%s" % (";".join(bad) if bad else "none"))
    bad2 = []
    for (a, b), cell in sorted(witness.items()):
        if a not in cfg["z"] or b not in cfg["z"]:
            bad2.append("%s-%s:not a rectangle of the object" % (a, b))
            continue
        ok, why = disk_witness_ok(cfg, a, b, cell)
        if not ok:
            bad2.append("%s-%s:%s" % (a, b, why))
    check("every_paper_witness_cell_yields_a_valid_cylinder", not bad2,
          "bad=%s" % (";".join(bad2) if bad2 else "none"))
    return keys


def check_blocker_table(cfg, cells, g, blockers, wkeys):
    """(D6b) The paper's third table: blockers really are intermediate, really
    cover the intersection, and the "empty" symbol marks exactly the pairs
    whose intersection of cell sets is empty."""
    keys = set(frozenset(k) for k in blockers)
    allp = set(frozenset(p) for p in all_pairs(cfg))
    check("blocker_table_has_36_entries",
          len(blockers) == 36 and len(keys) == 36,
          "entries=%d distinct=%d" % (len(blockers), len(keys)))
    check("blocker_table_keys_are_exactly_the_nonedges_of_G",
          keys == allp - g,
          "extra=%s missing=%s"
          % (",".join(sorted(show_pair(e) for e in keys - (allp - g)))
             or "none",
             ",".join(sorted(show_pair(e) for e in (allp - g) - keys))
             or "none"))
    check("the_two_certificate_tables_partition_all_66_pairs",
          (wkeys | keys) == allp and not (wkeys & keys)
          and len(allp) == 66,
          "union=%d overlap=%d" % (len(wkeys | keys), len(wkeys & keys)))

    notbetween, notcover, wrong_empty, wrong_nonempty, selfref = [], [], [], [], []
    unknown = []
    for (a, b), bl in sorted(blockers.items()):
        if a not in cfg["z"] or b not in cfg["z"]:
            unknown.append("%s-%s:pair not in the object" % (a, b))
            continue
        inter = set(cells[a]) & set(cells[b])
        mid = set(strictly_between(cfg, a, b))
        for r in bl:
            if r not in cfg["z"] or r not in cfg["labels"]:
                unknown.append("%s-%s:%s" % (a, b, r))
            elif r == a or r == b:
                selfref.append("%s-%s:%s" % (a, b, r))
            elif r not in mid:
                notbetween.append("%s-%s:%s(z=%d not in (%d,%d))"
                                  % (a, b, r, cfg["z"][r],
                                     min(cfg["z"][a], cfg["z"][b]),
                                     max(cfg["z"][a], cfg["z"][b])))
        cov = set()
        for r in bl:
            cov |= cells.get(r, set())
        if bl and not (inter <= cov):
            notcover.append("%s-%s:%s" % (a, b, sorted(inter - cov)))
        if not bl and inter:
            wrong_empty.append("%s-%s:%s" % (a, b, sorted(inter)))
        if bl and not inter:
            wrong_nonempty.append("%s-%s" % (a, b))
    check("every_listed_blocker_names_a_rectangle_of_the_object", not unknown,
          "unknown=%s" % (";".join(unknown) if unknown else "none"))
    check("no_listed_blocker_is_one_of_its_own_pair", not selfref,
          "bad=%s" % (";".join(selfref) if selfref else "none"))
    check("every_listed_blocker_lies_strictly_between_its_pair_in_z",
          not notbetween,
          "bad=%s" % (";".join(notbetween) if notbetween else "none"))
    check("listed_blockers_cover_C_A_cap_C_B_for_every_nonempty_entry",
          not notcover,
          "uncovered=%s" % (";".join(notcover) if notcover else "none"))
    check("the_empty_symbol_entries_really_have_C_A_cap_C_B_empty",
          not wrong_empty,
          "bad=%s" % (";".join(wrong_empty) if wrong_empty else "none"))
    check("the_nonempty_entries_really_have_C_A_cap_C_B_nonempty",
          not wrong_nonempty,
          "bad=%s" % (";".join(wrong_nonempty) if wrong_nonempty else "none"))
    return keys


def check_self_tests(cfg, cells, g, blockers):
    """(D8) Deliberate corruptions; each must be caught, so none of the
    checks above can be vacuous."""
    # 1. two rectangles pushed into the same plane, projections overlapping.
    c1 = clone(cfg)
    c1["z"]["h2"] = c1["z"]["h1"]
    check("selftest_equal_heights_break_pairwise_disjointness",
          not rectangles_disjoint(c1, "h1", "h2")
          and len(set(c1["z"].values())) == 11)

    # 2. one rectangle deleted.
    c2 = clone(cfg)
    c2["labels"].remove("v2")
    check("selftest_deleting_a_rectangle_is_detected",
          len(c2["labels"]) != 12
          and visibility_graph_by_cells(c2) != g)

    # 3. a height moved to the top of the stack.
    c3 = clone(cfg)
    c3["z"]["h0"] = 12
    check("selftest_moving_h0_above_the_stack_breaks_the_representation",
          visibility_graph_by_cells(c3) != g
          and visibility_graph_by_area(c3) != g)

    # 4. one projection shrunk by a single unit.
    c4 = clone(cfg)
    c4["box"]["h0"] = (7, 9, 0, 7)
    check("selftest_shrinking_h0_to_x_in_7_9_breaks_the_representation",
          visibility_graph_by_cells(c4) != g
          and visibility_graph_by_area(c4) != g)

    # 5. a bogus witness cell must not be accepted.
    ok_bogus, _ = disk_witness_ok(cfg, "h0", "v1", (0, 0))
    check("selftest_a_bogus_witness_cell_is_rejected",
          (0, 0) not in free_cells(cfg, cells, "h0", "v1")
          and not ok_bogus)

    # 6. the covering test is tight: for every non-empty blocker entry,
    #    dropping any single listed blocker leaves the intersection uncovered.
    slack = []
    for (a, b), bl in sorted(blockers.items()):
        if not bl:
            continue
        if a not in cells or b not in cells or any(r not in cells for r in bl):
            slack.append("%s-%s:names a rectangle absent from the object"
                         % (a, b))
            continue
        inter = set(cells[a]) & set(cells[b])
        for drop in bl:
            cov = set()
            for r in bl:
                if r != drop:
                    cov |= cells[r]
            if inter <= cov:
                slack.append("%s-%s:%s" % (a, b, drop))
    check("selftest_every_listed_blocker_is_necessary_so_the_cover_is_tight",
          not slack,
          "redundant=%s" % (";".join(slack) if slack else "none"))


def check_mutation_census(cfg, g):
    """(D8, continued) Exhaustive census over all single-endpoint +-1
    perturbations of the 12 projections and all 66 transpositions of the
    height order.  This is a rigidity measurement of the exhibited object and
    a demonstration that the load-bearing test really discriminates."""
    # Every mutant is also an independent test bed for the equivalence of the
    # three methods.  On the true object that equivalence is a single data
    # point; here it is re-tested on 162 further configurations, so a coding
    # error in any one method has 162 chances to surface.
    crossbad = []

    def graphs_agree(c, tag):
        g1 = visibility_graph_by_cells(c)
        g2 = visibility_graph_by_area(c)
        g3 = visibility_graph_by_incl_excl(c)
        if not (g1 == g2 == g3):
            crossbad.append(tag)
        return g1

    tried = broke = degenerate = 0
    survivors = []
    for l in cfg["labels"]:
        for k in range(4):
            for d in (-1, 1):
                v = list(cfg["box"][l])
                v[k] += d
                if v[0] >= v[1] or v[2] >= v[3]:
                    degenerate += 1
                    continue
                c = clone(cfg)
                c["box"][l] = tuple(v)
                tried += 1
                if graphs_agree(c, "box:%s%+d@%d" % (l, d, k)) != g:
                    broke += 1
                else:
                    survivors.append("%s->[%d,%d]x[%d,%d]"
                                     % (l, v[0], v[1], v[2], v[3]))
    info("endpoint perturbation census: %d admissible, %d break the "
         "representation, %d preserve it, %d rejected as degenerate"
         % (tried, broke, tried - broke, degenerate))
    check("mutation_census_endpoint_perturbations_are_discriminating",
          tried == 96 and 2 * broke > tried,
          "tried=%d broke=%d (a majority must break)" % (tried, broke))

    ttried = tbroke = 0
    for a, b in all_pairs(cfg):
        c = clone(cfg)
        c["z"][a], c["z"][b] = cfg["z"][b], cfg["z"][a]
        ttried += 1
        if graphs_agree(c, "z:%s<->%s" % (a, b)) != g:
            tbroke += 1
    info("height transposition census: %d transpositions, %d break the "
         "representation, %d preserve it"
         % (ttried, tbroke, ttried - tbroke))
    check("mutation_census_height_transpositions_are_discriminating",
          ttried == 66 and 2 * tbroke > ttried,
          "tried=%d broke=%d (a majority must break)" % (ttried, tbroke))
    check("all_three_methods_agree_on_every_mutant_of_the_census",
          not crossbad and tried + ttried == 162,
          "disagreements=%s over %d mutants"
          % (",".join(crossbad[:6]) if crossbad else "none", tried + ttried))
    if survivors:
        info("perturbations that still represent G (%d): %s"
             % (len(survivors), ", ".join(survivors)))


def check_weaker_reading(cfg, g):
    """(D9) The Scope paragraph's one quantitative side claim, COMPUTED rather
    than asserted: that under the weaker reading of "joins", in which the end
    disks need only meet the two rectangles instead of being contained in them,
    the same twelve rectangles realise a STRICTLY LARGER graph and therefore do
    not represent G -- so the containment reading is not a normalisation of
    convenience.  The weakening evaluated here keeps the cylinder axis on both
    rectangles and drops the containment of the disks; see
    `touching_free_point` for why that is exactly nonemptiness of the closed
    free set, decided in exact arithmetic."""
    tg = touching_visibility_graph(cfg)
    extra = sorted(show_pair(e) for e in tg - g)
    missing = sorted(show_pair(e) for e in g - tg)
    check("weaker_reading_graph_contains_every_edge_of_G", not missing,
          "missing=%s" % (",".join(missing) if missing else "none"))
    check("weaker_reading_graph_is_STRICTLY_larger_than_G_so_the_containment"
          "_reading_is_load_bearing",
          tg > g,
          "|E(weak)|=%d vs |E(G)|=%d extra=%s"
          % (len(tg), len(g), ",".join(extra) if extra else "none"))

    # Not merely a set inequality: each additional pair is realised by an
    # explicit cylinder of exact positive radius under the weaker reading.
    bad = []
    shown = []
    for e in sorted(tg - g, key=show_pair):
        a, b = sorted(e)
        pt = touching_free_point(cfg, a, b)
        if pt is None:
            bad.append(show_pair(e) + ":no axis point")
            continue
        ok, r2, why = touching_disk_witness(cfg, a, b, pt)
        if not ok:
            bad.append(show_pair(e) + ":" + why)
        else:
            shown.append("%s axis=(%s,%s) r^2=%s"
                         % (show_pair(e), pt[0], pt[1], r2))
    check("every_extra_pair_of_the_weaker_reading_carries_an_explicit_positive"
          "_radius_disk",
          bool(shown) and not bad and len(shown) == len(tg - g),
          "witnessed=%d of %d failures=%s"
          % (len(shown), len(tg - g), ";".join(bad) if bad else "none"))

    hs = set(l for l in cfg["labels"] if l.startswith("h"))
    crossing = sorted(show_pair(e) for e in tg - g if len(set(e) & hs) == 1)
    check("the_extra_pairs_all_join_two_rectangles_of_the_SAME_part",
          bool(tg - g) and not crossing,
          "crossing=%s" % (",".join(crossing) if crossing else "none"))
    vismatch = sorted(show_pair(e) for e in MATCHING if e in tg)
    check("the_six_deleted_matching_pairs_stay_nonvisible_under_the_weaker"
          "_reading",
          len(MATCHING) == 6 and not vismatch,
          "visible_matching_pairs=%s"
          % (",".join(vismatch) if vismatch else "none"))

    # Falsifiability of the routine itself: the extra pair h4-h5 exists only
    # because those two projections share the line y=3.  Lift h4 off that line
    # and the extra pair must disappear.
    mut = clone(cfg)
    if "h4" in mut["box"] and "h5" in mut["box"]:
        x0, x1, y0, y1 = mut["box"]["h4"]
        mut["box"]["h4"] = (x0, x1, y0 + 1, y1)
    mtg = touching_visibility_graph(mut)
    hh = frozenset(("h4", "h5"))
    check("selftest_lifting_h4_off_the_line_it_shares_with_h5_removes_that"
          "_extra_pair",
          hh in tg and hh not in mtg,
          "before=%s after=%s" % ("visible" if hh in tg else "not visible",
                                  "visible" if hh in mtg else "not visible"))

    info("weaker-reading graph: %d pairs visible against %d under the "
         "containment reading; the %d additional pair(s): %s"
         % (len(tg), len(g), len(tg - g),
            ", ".join(extra) if extra else "none"))
    if shown:
        info("cylinders realising them: " + "; ".join(shown))
    info("The blocking condition is left at its most restrictive: an "
         "intermediate projection blocks as a CLOSED set, so the additional "
         "pairs above are not an artefact of relaxing 'meets'.  And weakening "
         "'joins' only ever ADMITS cylinders, so under any reading at least as "
         "weak as the one computed here -- including one that lets the end "
         "disks meet the two rectangles without their axis touching either -- "
         "the realised graph contains the graph above and is therefore "
         "strictly larger than G as well.")
    return tg


def finish():
    n = len(_RESULTS)
    bad = [nm for nm, ok in _RESULTS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


def main():
    print("verify.py -- z-visibility representation of K_{6,6} minus a "
          "perfect matching")
    print("all arithmetic exact (int / Fraction); standard library only")
    print("")

    cfg = make_config(PAPER_RECTS)

    print("--- 1. the exhibited object, decoded and checked ---")
    bb, usable = check_well_formed(cfg)
    cells = None
    if usable:
        cells = check_cell_sets(cfg, bb)
        print_object(cfg, cells)
    print("")

    print("--- 2. hypotheses of the statement being refuted ---")
    g = check_target_graph(cfg)
    info("the conjecture asserts NO z-visibility representation exists for "
         "this graph; a single representation refutes it")
    print("")

    if not usable:
        info("the exhibited object does not carry the twelve labels h0..h5, "
             "v0..v5, so every routine below is undefined on it; the run "
             "stops here with a verdict rather than an interpreter traceback.")
        print("")
        return finish()

    print("--- 3. LOAD-BEARING: the visibility graph, unit-cell method ---")
    vis, fsize = check_load_bearing_cells(cfg, cells, g)
    check_represented_graph_shape(cfg, vis)
    print("")

    print("--- 4. LOAD-BEARING: the visibility graph, exact geometry, "
          "independent of the paper's Lemma ---")
    check_load_bearing_geometry(cfg, cells, g, fsize)
    print("")

    print("--- 5. the paper's own finite certificate, both tables ---")
    wkeys = check_witness_table(cfg, cells, g, PAPER_WITNESS)
    check_blocker_table(cfg, cells, g, PAPER_BLOCKERS, wkeys)
    print("")

    print("--- 6. falsifiability self-tests and mutation census ---")
    check_self_tests(cfg, cells, g, PAPER_BLOCKERS)
    check_mutation_census(cfg, g)
    print("")

    print("--- 7. the weaker reading of 'joins', computed not asserted ---")
    check_weaker_reading(cfg, g)
    print("")

    print("--- 8. scope ---")
    info("The paper's claim is a single exhibited object plus the finite "
         "66-pair verification of it; both are reproduced here in full, "
         "with no reduction of range.")
    info("READING OF THE DEFINITION, and it is load-bearing: 'a closed "
         "vertical cylinder joins A and B' is taken, as in Bose et al., to "
         "mean that the cylinder's two end disks are CONTAINED in A and in B "
         "(D subset pi(A) and D subset pi(B)).  Both directions above use that "
         "reading: the 30 witnesses exhibit such a disk, and the 36 "
         "non-visibility proofs bound the area of the set of admissible "
         "centres.  That the weaker reading, in which the end disks need only "
         "MEET A and B, makes the twelve rectangles realise a strictly larger "
         "graph -- so that the containment reading is not a harmless "
         "normalisation but the definition being refuted -- is not merely "
         "asserted here: section 7 computes that graph, exhibits a cylinder "
         "for each additional pair, and reports the additional pairs by name.")
    info("Insensitive, by contrast, to whether 'meets another rectangle' is "
         "read as closed or as interior intersection: the witnesses keep "
         "clearance 1/4 from every intermediate projection, and a zero-area "
         "free region admits no positive-radius disk under either reading.")
    info("Not reproduced, because they are not computations: the "
         "bibliographic claim that Conjecture 6.5 of Bose et al. (1994) is "
         "the statement quoted, and the general proof of the paper's "
         "unit-cell Lemma (whose conclusion is nevertheless confirmed on "
         "this instance by the independent exact-geometry method in "
         "section 4); and, in section 7, the step from the ONE weakening of "
         "'joins' computed there -- cylinder axis on both rectangles, end "
         "disks not contained -- to every weaker reading, which is the "
         "one-line monotonicity remark printed there and not a search over "
         "cylinders whose axis touches neither rectangle.")
    print("")
    return finish()


if __name__ == "__main__":
    sys.exit(main())

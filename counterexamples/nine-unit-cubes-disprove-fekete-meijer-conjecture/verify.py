#!/usr/bin/env python3
"""Exact verification of a unit-cube box-visibility representation of K_9.

Setting.  A box-visibility representation assigns pairwise disjoint
axis-parallel boxes in R^3 to the vertices of a graph; two vertices are
adjacent exactly when their boxes are joined by an axis-parallel cylindrical
channel of positive radius whose end disks lie in the two box interiors and
which meets no other box.  The refuted conjecture asserts that no such
representation of K_9 exists using nine unit cubes.

TAKEN FROM THE PAPER (inputs, transcribed literally, nothing else):
  * scale 16 and the nine lower corners (x_i, y_i, z_i), i = 0..8, each cube
    being [x_i, x_i+16] x [y_i, y_i+16] x [z_i, z_i+16];
  * the table of 36 visibility axes, one per unordered pair, transcribed as
    the strings "ij:d(a,b)" meaning the line l_d(a,b) with the paper's
    coordinate order (l_x(a,b) = {(t,a,b)}, l_y(a,b) = {(a,t,b)},
    l_z(a,b) = {(a,b,t)});
  * the channel radius 1/4 at scale 16 and the end extension 1/4.

DERIVED HERE (computed, never assumed):
  * the nine cubes are decoded, counted and printed back; every side length
    is 16, so after dividing by 16 they are congruent axis-parallel unit
    cubes;
  * the cubes are pairwise disjoint, by exhibiting a separating direction
    for each of the 36 pairs;
  * the paper's structural remarks: exactly four pairs have overlapping
    z-intervals and which they are; all 18 face coordinates are distinct in
    each direction;
  * for every one of the 36 pairs, the cylinder named by the table is built
    from the coordinates and checked to have both end disks strictly inside
    the two cube interiors and to be disjoint from all seven other cubes
    (exact rational arithmetic, squared distances, no floats);
  * INDEPENDENTLY of the table, the full box-visibility graph of the nine
    cubes is computed by an exact decision procedure over the integer-grid
    cell arrangement of the projections, and is found to be complete on 9
    vertices (36 edges);
  * that decision procedure is re-run on the configuration blown up by 4, so
    that it samples a lattice four times finer, and returns the same graph;
  * the same certificate re-verified in unit-scale rational coordinates with
    radius 1/64;
  * five negative controls, one per clause of the definition, showing that
    the disjointness test, the certificate verifier (its transverse
    containment branch, its end-disk branch and its blocking branch, each
    demanded by name so that switching one off is noticed) and the
    independent search all do report failure on corrupted input.

NOT RE-RUN (external, cited by the paper, outside this program's scope):
  the cited theorem that K_10 admits no unit-cube box-visibility
  representation, and the cited K_8 construction.  Only the K_9 lower bound
  is re-established here; the matching upper bound is quoted.
"""

from fractions import Fraction as F
import sys

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)


def stage(fn, *args):
    """Run one group of checks; an unexpected exception is itself a failure."""
    try:
        return fn(*args)
    except Exception as exc:                       # corrupted input, etc.
        ck("stage_" + fn.__name__ + "_completed", False,
           "%s: %s" % (type(exc).__name__, exc))
        return None


def finish():
    n = len(CHECKS)
    bad = [c for c, ok in CHECKS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


# ---------------------------------------------------------------- paper data
SCALE = 16
# lower corners, exactly as displayed in the paper
X = [0, 9, -8, 10, 6, 2, 7, -5, 5]
Y = [0, 2, 3, -1, -7, 11, -6, 8, 1]
Z = [0, 17, 19, 34, 51, 53, 68, 85, 102]
# upper corners, transcribed separately so that "side 16" is a real test
XU = [16, 25, 8, 26, 22, 18, 23, 11, 21]
YU = [16, 18, 19, 15, 9, 27, 10, 24, 17]
ZU = [16, 33, 35, 50, 67, 69, 84, 101, 118]

def build_cubes():
    """A cube is a triple of closed intervals; no width is assumed anywhere."""
    n = min(len(X), len(Y), len(Z), len(XU), len(YU), len(ZU))
    return [((X[i], XU[i]), (Y[i], YU[i]), (Z[i], ZU[i])) for i in range(n)]


CUBES = build_cubes()
RADIUS = F(1, 4)          # channel radius at scale 16
EXTEND = F(1, 4)          # how far each channel reaches into the two cubes
PAPER_MARGIN = F(1, 2)    # the paper's asserted uniform sup-norm margin

TABLE_RAW = [
    "01:z(19/2,5/2)", "02:z(1/2,7/2)", "03:z(21/2,1/2)",
    "04:z(13/2,1/2)", "05:z(17/2,23/2)", "06:z(17/2,19/2)",
    "07:z(17/2,21/2)", "08:z(11/2,3/2)", "12:x(7/2,39/2)",
    "13:z(21/2,5/2)", "14:z(19/2,5/2)", "15:z(19/2,23/2)",
    "16:z(19/2,19/2)", "17:z(19/2,21/2)", "18:z(37/2,31/2)",
    "23:x(7/2,69/2)", "24:z(13/2,7/2)", "25:z(5/2,23/2)",
    "26:z(15/2,19/2)", "27:z(-9/2,17/2)", "28:z(11/2,7/2)",
    "34:z(21/2,-1/2)", "35:z(21/2,23/2)", "36:z(21/2,19/2)",
    "37:z(21/2,21/2)", "38:z(23/2,21/2)", "45:y(13/2,107/2)",
    "46:z(15/2,-11/2)", "47:z(13/2,17/2)", "48:z(13/2,3/2)",
    "56:y(15/2,137/2)", "57:z(5/2,23/2)", "58:z(23/2,23/2)",
    "67:z(15/2,17/2)", "68:z(15/2,3/2)", "78:z(11/2,17/2)",
]

AXIS = {"x": 0, "y": 1, "z": 2}
TRANSVERSE = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}


def parse_table(raw):
    """Decode the paper's labels into {(i,j): (direction, a, b)}."""
    out = {}
    for s in raw:
        head, rest = s.split(":", 1)
        i, j = int(head[0]), int(head[1])
        d = rest[0]
        a, b = rest[2:-1].split(",")
        if d not in AXIS or not (0 <= i < 9 and 0 <= j < 9) or i >= j:
            raise ValueError(s)
        out[(i, j)] = (d, F(a), F(b))
    return out


try:
    TABLE, TABLE_ERR = parse_table(TABLE_RAW), ""
except Exception as _exc:                          # malformed label
    TABLE, TABLE_ERR = {}, "%s: %s" % (type(_exc).__name__, _exc)
PAIRS = [(i, j) for i in range(9) for j in range(i + 1, 9)]


def iv(cubes, k, ax):
    return cubes[k][ax]


def rect(cubes, k, d):
    return tuple(iv(cubes, k, ax) for ax in TRANSVERSE[d])


def inner_margin(q, r):
    """sup-norm distance from q to the complement of the open rectangle r
    (negative if q is outside)."""
    return min(min(q[t] - r[t][0], r[t][1] - q[t]) for t in (0, 1))


def outer_gaps(q, r):
    """per-axis sup-norm gaps from q to the closed rectangle r."""
    return [max(r[t][0] - q[t], q[t] - r[t][1], 0) for t in (0, 1)]


def verify_channel(cubes, i, j, d, q, radius, extend):
    """Build the named cylinder and test it. Returns (problems, margins)."""
    ax = AXIS[d]
    lo_i, hi_i = iv(cubes, i, ax)
    lo_j, hi_j = iv(cubes, j, ax)
    bad, margins = [], []
    if hi_i < lo_j:
        low, high = i, j
    elif hi_j < lo_i:
        low, high = j, i
    else:
        return (["%d%d: %s-intervals of the two cubes are not disjoint"
                 % (i, j, d)], margins)
    u = iv(cubes, low, ax)[1]           # top face of the lower cube
    v = iv(cubes, high, ax)[0]          # bottom face of the upper cube
    p0, p1 = u - extend, v + extend     # the cylinder's span along d
    if p1 < p0:
        bad.append("%d%d: empty cylinder span" % (i, j))
    if not (iv(cubes, low, ax)[0] < p0 < u):
        bad.append("%d%d: lower end disk not inside the lower cube" % (i, j))
    if not (v < p1 < iv(cubes, high, ax)[1]):
        bad.append("%d%d: upper end disk not inside the upper cube" % (i, j))
    for h in (i, j):
        m = inner_margin(q, rect(cubes, h, d))
        margins.append(m)
        if not m > radius:               # closed disk inside the open face
            bad.append("%d%d: axis not interior to cube %d (margin %s)"
                       % (i, j, h, m))
    for k in range(len(cubes)):
        if k in (i, j):
            continue
        k_lo, k_hi = iv(cubes, k, ax)
        if k_hi < p0 or k_lo > p1:
            continue                     # cube k is not at these d-heights
        gx, gy = outer_gaps(q, rect(cubes, k, d))
        margins.append(max(gx, gy))
        if gx * gx + gy * gy <= radius * radius:
            bad.append("%d%d: cylinder meets cube %d" % (i, j, k))
    return (bad, margins)


def verify_certificate(cubes, table, radius, extend):
    """Run verify_channel over every table entry."""
    bad, margins = [], []
    for (i, j) in sorted(table):
        d, a, b = table[(i, j)]
        b1, m1 = verify_channel(cubes, i, j, d, (a, b), radius, extend)
        bad.extend(b1)
        margins.extend(m1)
    return bad, margins


def separating_directions(cubes, i, j):
    out = []
    for d, ax in sorted(AXIS.items()):
        lo_i, hi_i = iv(cubes, i, ax)
        lo_j, hi_j = iv(cubes, j, ax)
        if hi_i < lo_j or hi_j < lo_i:
            out.append(d)
    return out


def visible_in_direction(cubes, i, j, d):
    """Exact decision: is there an axis-parallel channel along d from C_i to
    C_j of some positive radius meeting no other cube?

    All faces are at integer coordinates, so the transverse projections form
    an arrangement of the unit integer grid.  Any nonempty open set left
    after deleting the blocking rectangles from the open overlap contains a
    full open grid cell, hence its half-integer centre; and a cell meeting a
    closed integer-cornered rectangle lies inside it.  So testing the
    half-integer centres of the overlap decides the question exactly.
    """
    ax = AXIS[d]
    lo_i, hi_i = iv(cubes, i, ax)
    lo_j, hi_j = iv(cubes, j, ax)
    if hi_i < lo_j:
        low, high = i, j
    elif hi_j < lo_i:
        low, high = j, i
    else:
        return None
    u, v = iv(cubes, low, ax)[1], iv(cubes, high, ax)[0]
    blockers = []
    for k in range(len(cubes)):
        if k in (i, j):
            continue
        k_lo, k_hi = iv(cubes, k, ax)
        if not (k_hi < u or k_lo > v):
            blockers.append(rect(cubes, k, d))
    ri, rj = rect(cubes, i, d), rect(cubes, j, d)
    span = [(max(ri[t][0], rj[t][0]), min(ri[t][1], rj[t][1])) for t in (0, 1)]
    for m in range(span[0][0], span[0][1]):
        for n in range(span[1][0], span[1][1]):
            q = (F(2 * m + 1, 2), F(2 * n + 1, 2))
            if any(all(r[t][0] <= q[t] <= r[t][1] for t in (0, 1))
                   for r in blockers):
                continue
            return (d, q)
    return None


def visibility_graph(cubes):
    """The full box-visibility edge set, computed from the cubes alone."""
    edges = {}
    n = len(cubes)
    for i in range(n):
        for j in range(i + 1, n):
            for d in ("x", "y", "z"):
                w = visible_in_direction(cubes, i, j, d)
                if w is not None:
                    edges[(i, j)] = w
                    break
    return edges


def check_decoding():
    print("the exhibited object, decoded and printed back (scale %d):" % SCALE)
    for i, c in enumerate(CUBES):
        print("  C_%d = [%4d,%4d] x [%4d,%4d] x [%4d,%4d]"
              % (i, c[0][0], c[0][1], c[1][0], c[1][1], c[2][0], c[2][1]))
    sizes = {len(v) for v in (X, Y, Z, XU, YU, ZU)}
    ck("cube_count_is_nine",
       len(CUBES) == 9 and len(set(CUBES)) == 9 and sizes == {9},
       "%d boxes, %d distinct, coordinate list lengths %s"
       % (len(CUBES), len(set(CUBES)), sorted(sizes)))
    ck("boxes_are_nondegenerate_products_of_closed_intervals",
       all(lo < hi for c in CUBES for (lo, hi) in c) and
       all(len(c) == 3 for c in CUBES),
       "each box is an interval product in x, y, z with lo < hi, hence "
       "axis-parallel with nonempty interior")
    widths = sorted({hi - lo for c in CUBES for (lo, hi) in c})
    ck("every_box_is_a_cube_of_side_16", widths == [SCALE],
       "the 27 edge lengths take the values %s" % widths)


def check_disjointness():
    bad, witness = [], {}
    for (i, j) in PAIRS:
        ds = separating_directions(CUBES, i, j)
        if not ds:
            bad.append((i, j))
        else:
            witness[(i, j)] = ds
    ck("cubes_pairwise_disjoint", not bad and len(witness) == 36,
       "36 pairs, each separated in >=1 direction"
       if not bad else "overlapping pairs: %s" % bad)
    return witness


def check_paper_structure(witness):
    overlap_z = sorted(p for p in PAIRS
                       if "z" not in separating_directions(CUBES, *p))
    ck("z_overlapping_pairs_are_exactly_the_four_listed",
       overlap_z == [(1, 2), (2, 3), (4, 5), (5, 6)],
       "pairs with meeting z-intervals: %s" % (overlap_z,))
    sep = {p: witness.get(p, []) for p in overlap_z}
    ok = (sep.get((1, 2)) == ["x"] and sep.get((2, 3)) == ["x"]
          and sep.get((4, 5)) == ["y"] and sep.get((5, 6)) == ["y"])
    ck("those_four_pairs_are_separated_in_x_x_y_y", ok, "%s" % (sep,))
    dups = []
    for ax, nm in ((0, "x"), (1, "y"), (2, "z")):
        vals = [v for c in CUBES for v in c[ax]]
        if len(set(vals)) != 18:
            dups.append(nm)
    ck("all_18_face_coordinates_distinct_in_each_direction", not dups,
       "repeats in: %s" % dups if dups else "18 distinct per direction")


def check_table_shape():
    ck("table_lists_one_axis_for_each_of_the_36_pairs",
       len(TABLE_RAW) == 36 and sorted(TABLE) == PAIRS and not TABLE_ERR,
       "%d labels, %d distinct pairs, C(9,2)=%d%s"
       % (len(TABLE_RAW), len(TABLE), len(PAIRS),
          "; decode error " + TABLE_ERR if TABLE_ERR else ""))
    halves = all(F(a).denominator == 2 and F(b).denominator == 2
                 for (_, a, b) in TABLE.values())
    dirs = sorted({d for (d, _, _) in TABLE.values()})
    ck("table_axes_are_half_integer_lines_in_x_y_or_z",
       halves and dirs == ["x", "y", "z"],
       "directions used: %s; every listed coordinate is a half-integer: %s"
       % (dirs, halves))


def check_certificate():
    bad, margins = verify_certificate(CUBES, TABLE, RADIUS, EXTEND)
    ck("all_36_exhibited_channels_are_valid", not bad,
       "36/36 cylinders of radius 1/4 have both end disks interior and miss "
       "all other cubes" if not bad else "; ".join(bad[:4]))
    worst = min(margins) if margins else None
    ck("uniform_sup_norm_margin_is_at_least_one_half",
       margins and worst >= PAPER_MARGIN,
       "smallest of the %d sup-norm margins = %s (paper claims >= 1/2), "
       "channel radius %s" % (len(margins), worst, RADIUS))


def check_visibility_graph():
    edges = visibility_graph(CUBES)
    ck("independent_search_finds_all_36_pairs_visible",
       sorted(edges) == PAIRS,
       "edges found %d of 36; missing %s"
       % (len(edges), [p for p in PAIRS if p not in edges]))
    deg = [0] * len(CUBES)
    for (i, j) in edges:
        deg[i] += 1
        deg[j] += 1
    complete = (len(CUBES) == 9 and len(edges) == 9 * 8 // 2
                and deg == [8] * 9)
    ck("box_visibility_graph_equals_K_9", complete,
       "%d vertices, %d edges, degree sequence %s -> K_9"
       % (len(CUBES), len(edges), deg))
    agree = sum(1 for p in PAIRS
                if p in edges and p in TABLE
                and edges[p][0] == TABLE[p][0])
    ck("independent_witnesses_use_the_paper_s_directions", agree == 36,
       "%d of 36 pairs first become visible along the direction the paper "
       "names" % agree)
    return edges


def blown_up(cubes, fac):
    return [tuple((lo * fac, hi * fac) for (lo, hi) in c) for c in cubes]


def check_grid_refinement(edges):
    """The decision procedure tests one point per unit cell.  Blowing the
    configuration up by 4 makes it test a lattice four times finer, which
    must neither miss nor invent edges."""
    fine = visibility_graph(blown_up(CUBES, 4))
    ck("search_agrees_under_fourfold_grid_refinement",
       sorted(fine) == sorted(edges) and len(fine) == 36,
       "coarse %d edges, refined %d edges, identical: %s, both complete: %s"
       % (len(edges), len(fine), sorted(fine) == sorted(edges),
          len(fine) == 36 and len(edges) == 36))
    stack = [((0, 16), (0, 16), (17 * i, 17 * i + 16)) for i in range(9)]
    fine_stack = visibility_graph(blown_up(stack, 4))
    ck("refinement_invents_no_edge_in_the_blocked_stack",
       len(fine_stack) == 8,
       "the column still has %d edges on the finer lattice"
       % len(fine_stack))


def scaled(cubes, den):
    return [tuple((F(lo, den), F(hi, den)) for (lo, hi) in c) for c in cubes]


def check_unit_scale():
    """Divide everything by 16: unit cubes, channels of radius 1/64."""
    unit = scaled(CUBES, SCALE)
    widths = sorted({hi - lo for c in unit for (lo, hi) in c})
    tbl = {p: (d, a / SCALE, b / SCALE) for p, (d, a, b) in TABLE.items()}
    r, e = RADIUS / SCALE, EXTEND / SCALE
    bad, margins = verify_certificate(unit, tbl, r, e)
    ck("unit_cubes_after_dividing_by_16", widths == [F(1)],
       "the 27 edge lengths all equal %s" % widths[0])
    ck("certificate_still_valid_at_unit_scale_radius_1_64",
       not bad and r == F(1, 64) and min(margins) == F(1, 32),
       "36/36 channels of radius %s, smallest margin %s"
       % (r, min(margins)) if not bad else "; ".join(bad[:3]))
    dis = all(separating_directions(unit, i, j) for (i, j) in PAIRS)
    ck("unit_cubes_still_pairwise_disjoint", dis,
       "all 36 pairs separated in some direction at unit scale")


def check_negative_controls():
    # (1) two cubes made to coincide: the disjointness test must object.
    bad_cubes = list(CUBES)
    bad_cubes[2] = CUBES[1]
    ok1 = not separating_directions(bad_cubes, 1, 2)
    ck("control_overlapping_cubes_are_detected", ok1,
       "cloning C_1 onto C_2 leaves no separating direction")
    # (2) one cube shifted by 1 in x: the certificate must be rejected.
    moved = list(CUBES)
    moved[5] = ((X[5] + 1, XU[5] + 1), CUBES[5][1], CUBES[5][2])
    bad2, _ = verify_certificate(moved, TABLE, RADIUS, EXTEND)
    # both clauses of the definition must be seen to reject: the axis leaving a
    # cube's projected face, and the cylinder running into a third cube.  A
    # bare "some test failed" would still pass with the containment clause
    # switched off entirely, so the two failure modes are demanded by name.
    hit_containment = any("axis not interior to cube 5" in m for m in bad2)
    hit_blocking = any("meets cube 5" in m for m in bad2)
    ck("control_perturbed_coordinates_are_rejected",
       len(bad2) > 0 and hit_containment and hit_blocking,
       "shifting C_5 by +1 in x breaks %d of the 36 channel tests "
       "(containment clause fired: %s, blocking clause fired: %s)"
       % (len(bad2), hit_containment, hit_blocking))
    # (3) an axis moved onto a position screened by C_2: the blocking
    #     branch (not the containment branch) must fire.
    tampered = dict(TABLE)
    tampered[(0, 8)] = ("z", F(11, 2), F(7, 2))
    bad3, _ = verify_certificate(CUBES, tampered, RADIUS, EXTEND)
    ck("control_blocked_channel_is_rejected",
       any("meets cube 2" in m for m in bad3),
       "moving the 0-8 axis to (11/2,7/2) is reported as blocked: %s"
       % ("; ".join(bad3) if bad3 else "NOT REPORTED"))
    # (4) the end-disk clause.  For cubes of side 16 the end planes u - 1/4 and
    #     v + 1/4 always fall inside, so that clause never binds on the
    #     exhibited object; it is exercised here on a deliberately thin slab
    #     (thickness 1/8 < 1/4) where the end disk must be reported outside.
    fat, thin = (F(1), F(17)), (F(0), F(1, 8))
    sq = ((F(0), F(4)), (F(0), F(4)))
    low_thin, _ = verify_channel([sq + (thin,), sq + (fat,)], 0, 1, "z",
                                 (F(2), F(2)), RADIUS, EXTEND)
    high_thin, _ = verify_channel([sq + ((F(0), F(16)),), sq + ((F(17), F(137, 8)),)],
                                  0, 1, "z", (F(2), F(2)), RADIUS, EXTEND)
    msgs = low_thin + high_thin
    ck("control_end_disk_outside_a_thin_box_is_rejected",
       any("lower end disk not inside" in m for m in low_thin)
       and any("upper end disk not inside" in m for m in high_thin)
       and not any("axis not interior" in m for m in msgs),
       "a slab of thickness 1/8 cannot hold an end disk 1/4 deep: %s"
       % ("; ".join(msgs) if msgs else "NOT REPORTED"))
    # (5) nine cubes stacked in a column: only consecutive pairs are visible.
    stack = [((0, 16), (0, 16), (17 * i, 17 * i + 16)) for i in range(9)]
    edges = visibility_graph(stack)
    ck("control_blocked_stack_is_not_complete",
       len(edges) == 8 and sorted(edges) == [(i, i + 1) for i in range(8)],
       "a straight column of nine unit cubes yields %d edges, not 36"
       % len(edges))


def main():
    stage(check_decoding)
    witness = stage(check_disjointness) or {}
    stage(check_paper_structure, witness)
    stage(check_table_shape)
    stage(check_certificate)
    edges = stage(check_visibility_graph) or {}
    stage(check_grid_refinement, edges)
    stage(check_unit_scale)
    stage(check_negative_controls)
    if all(ok for _, ok in CHECKS):
        print("nine pairwise disjoint axis-parallel unit cubes whose "
              "box-visibility graph is complete: the K_9 conjecture is "
              "refuted.")
    else:
        print("the exhibited configuration did NOT verify; see the failures "
              "above.")
    print("NOT RE-RUN: the cited nonexistence result for K_10 and the cited "
          "K_8 construction are external and are not recomputed here; only "
          "the K_9 representation is verified.")
    finish()


if __name__ == "__main__":
    main()

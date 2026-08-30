#!/usr/bin/env python3
"""verify.py -- re-derives every computational claim of

    "An Exhaustive Census Confirms a(6) = 1814 for Topologically Distinct Pizza Slicings"

Python 3.9 or later, STANDARD LIBRARY ONLY (fractions, itertools, collections, sys, time).
No third-party package, no external data file, no floating point in any decision: every
coordinate, every crossing parameter and every comparison is exact integer or Fraction
arithmetic.

WHAT IT READS.  Everything it consumes is printed in the paper:
  * PAPER_TABLE1  -- the 47 twelve-tuples of Table 1 of the paper;
  * PAPER_WITNESS_* -- the single witness exhibited in full in Section 4 (Example 1),
    including its 12 circle points, its 30 crossing parameters, its crossing orders and
    its region-graph degree sequence;
  * PAPER_* -- every integer and every table the paper asserts.
Nothing is read from disk.

WHAT IT DOES.  Two halves, matching the paper's two bounds.

  UPPER BOUND.  It re-runs the topological census of chord arrangements in a disk by
  incremental insertion (Section 3), obtains the arrangement types level by level for
  n = 1..6, builds each type's region graph, and counts the graphs up to isomorphism by
  TWO independent routes (a refinement-plus-backtracking canonical form, and a
  Weisfeiler-Leman/distance/triangle invariant followed by an explicit exact isomorphism
  search).  It checks the published confirmed prefix a(2..5) = 2, 5, 19, 130 falls out of
  the same code.

  LOWER BOUND.  It sweeps every one of the 47 x 10395 = 488,565 (point set, perfect
  matching) pairs of Table 1, rejects a pair exactly when two crossings coincide on a
  chord (that IS the no-interior-triple-point clause), computes the topological type of
  each surviving arrangement, and checks (a) that no arrangement falls outside the census
  -- a live falsification channel, one exception would prove a(6) >= 1815 -- and (b) that
  all 1824 census types, hence all 1814 region-graph classes, are realized by straight
  chords with rational endpoints.

  It then verifies the exhibited witness of Example 1 in complete detail, and runs four
  controls in both polarities (a must-fire degeneracy detector, a must-not-change mirror
  symmetry, an independent Prufer enumeration of the trees anchoring the X = 0 shard, and
  a deliberately wrong adjacency rule that must give 18 rather than 19 at n = 4).

CREDIT.  The insertion enumerator, the half-edge map and the two canonical forms follow
the design of the independent enumerator written for this result; they are restated here
in full so that this file stands alone.

Exit status is 0 if and only if every check passed.
"""

import itertools
import sys
import time
from collections import defaultdict, deque
from fractions import Fraction as F

# ===========================================================================
# 0.  WHAT THE PAPER CLAIMS.  Every number below is quoted from the paper and is
#     compared against a quantity this program derives.
# ===========================================================================

PAPER_A272906_PREFIX = {2: 2, 3: 5, 4: 19, 5: 130}      # OEIS A272906, confirmed terms
PAPER_A272906_6 = 1814                                   # the term the paper confirms
PAPER_TYPES_6 = 1824                                     # topological types at n = 6
PAPER_TYPES_PERX = [12, 30, 51, 87, 108, 142, 165, 179,
                    184, 187, 182, 148, 129, 101, 76, 43]
PAPER_CLASSES_PERX = [11, 28, 47, 85, 107, 142, 165, 179,
                      184, 187, 182, 148, 129, 101, 76, 43]
PAPER_COLLAPSES = 10                                     # types minus classes
PAPER_COLLAPSE_SHARDS = {0: 1, 1: 2, 2: 4, 3: 2, 4: 1}   # X -> types lost to a collision
PAPER_A000055_7 = 11                                     # trees on 7 unlabelled nodes
PAPER_A090338_6 = 43                                     # simple arrangements of 6 lines
PAPER_WRONG_ADJACENCY_N4 = 18                            # the vertex rule undercounts 19
PAPER_MATCHINGS_12 = 10395                               # (2k-1)!! for k = 6
PAPER_SWEEP_PAIRS = 488565                               # 47 * 10395
PAPER_SWEEP_ADMISSIBLE = 487831                          # pairs surviving the triple-point test
PAPER_SWEEP_REJECTED = 734                               # pairs with an interior multiple point
PAPER_SWEEP_OUTSIDE_CENSUS = 0                           # the falsification channel

# Table 1 of the paper: 47 strictly increasing twelve-tuples of integers.  Row 1 is the
# witness exhibited in full in Example 1.
PAPER_TABLE1 = [
    [-7, -5, -4, -3, -2, -1, 2, 3, 4, 7, 8, 9],
    [-60, -50, -39, -37, -36, -20, -3, 21, 31, 39, 54, 60],
    [-59, -51, -18, 1, 5, 6, 11, 14, 19, 20, 43, 47],
    [-58, -51, -12, 0, 3, 4, 23, 31, 33, 51, 58, 59],
    [-55, -47, -40, -30, -16, -4, 0, 6, 50, 51, 52, 57],
    [-55, -40, -39, -36, -31, -28, -26, -18, 15, 24, 57, 60],
    [-53, -47, -39, -37, -34, -22, -20, 4, 11, 16, 27, 45],
    [-53, -45, -33, -28, -25, -5, 17, 20, 44, 53, 57, 60],
    [-51, -45, -44, -36, -35, 0, 5, 27, 36, 39, 50, 54],
    [-49, -44, -11, -8, 13, 25, 26, 27, 28, 35, 46, 59],
    [-49, -29, -27, -24, -23, -9, 21, 35, 46, 54, 57, 60],
    [-40, -38, -18, -9, -1, 7, 9, 14, 16, 17, 24, 34],
    [-40, -37, -34, -28, -27, -21, 3, 6, 28, 32, 36, 38],
    [-40, -36, -35, -27, -22, -14, -13, 3, 4, 7, 36, 39],
    [-40, -36, -32, -23, -21, -15, -14, -12, 2, 4, 21, 23],
    [-40, -34, -32, -30, -2, 8, 15, 17, 19, 21, 28, 40],
    [-40, -33, -25, -22, -21, -18, 1, 20, 26, 31, 37, 39],
    [-40, -30, -17, -15, -13, -9, -5, 2, 5, 8, 20, 24],
    [-40, -25, -22, -2, 0, 1, 2, 3, 4, 8, 10, 13],
    [-39, -35, -26, -24, -20, -10, -6, 3, 13, 27, 30, 39],
    [-39, -34, -31, -24, -23, -15, -8, 3, 4, 11, 15, 40],
    [-38, -37, -35, -30, -27, -23, -11, 6, 8, 17, 31, 34],
    [-38, -32, -3, 2, 3, 9, 16, 18, 25, 26, 30, 39],
    [-38, -29, -27, -23, -19, -17, -12, -9, -8, -1, 28, 39],
    [-38, -29, -20, -15, 0, 3, 5, 6, 11, 13, 15, 25],
    [-38, -28, -14, 3, 7, 12, 18, 20, 25, 26, 29, 33],
    [-38, -24, -21, -17, 13, 16, 24, 25, 26, 27, 28, 34],
    [-37, -34, -32, -30, -26, -23, -16, 0, 8, 24, 25, 38],
    [-37, -32, -15, -5, 1, 16, 17, 24, 25, 31, 37, 38],
    [-36, -35, -28, -24, -22, -13, 2, 13, 21, 22, 25, 40],
    [-34, -33, -31, -28, -21, -13, 1, 6, 10, 24, 28, 34],
    [-34, -29, -25, -14, -1, 4, 5, 8, 13, 15, 20, 31],
    [-33, -32, -28, -16, -14, 7, 23, 28, 30, 32, 33, 39],
    [-33, -29, -21, -17, -5, 1, 11, 14, 23, 25, 28, 33],
    [-33, -25, -17, -14, -8, -5, 5, 8, 10, 18, 25, 38],
    [-32, -31, -29, -14, -11, -10, -8, -7, 18, 20, 27, 28],
    [-32, -31, -28, -22, -19, -11, 5, 9, 11, 24, 25, 35],
    [-32, -24, -19, -16, -13, -12, -2, 10, 11, 19, 20, 27],
    [-32, -20, -15, -14, -9, -8, -3, -2, 7, 30, 35, 36],
    [-31, -30, -28, -27, -20, -14, -7, 13, 14, 17, 23, 39],
    [-31, -25, -19, -4, -2, 3, 13, 17, 23, 25, 27, 37],
    [-30, -29, -26, -15, -13, -2, 16, 23, 27, 30, 37, 38],
    [-30, -29, -24, -20, -19, -18, 2, 10, 11, 15, 19, 21],
    [-30, -26, -5, -4, 6, 7, 28, 31, 32, 37, 38, 39],
    [-28, -26, -16, -12, -11, -8, 18, 22, 27, 29, 30, 33],
    [-27, -25, -23, -22, -17, -3, -1, 13, 29, 31, 33, 34],
    [-27, -13, -11, -9, -1, 10, 15, 17, 20, 23, 24, 29],
]

# --- Example 1 of the paper, exhibited in full ------------------------------
PAPER_WITNESS_T = [-7, -5, -4, -3, -2, -1, 2, 3, 4, 7, 8, 9]
PAPER_WITNESS_MATCHING = [(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)]
PAPER_WITNESS_POINTS = [
    (F(-24, 25), F(-7, 25)), (F(-12, 13), F(-5, 13)), (F(-15, 17), F(-8, 17)),
    (F(-4, 5), F(-3, 5)), (F(-3, 5), F(-4, 5)), (F(0), F(-1)),
    (F(-3, 5), F(4, 5)), (F(-4, 5), F(3, 5)), (F(-15, 17), F(8, 17)),
    (F(-24, 25), F(7, 25)), (F(-63, 65), F(16, 65)), (F(-40, 41), F(9, 41)),
]
# the crossing parameter of chord d along chord c, measured from c's lower-indexed
# endpoint, listed for each chord c in increasing order -- so the second component of
# each pair is the crossing ORDER along c.
PAPER_WITNESS_S = [
    [(F(28, 153), 3), (F(11, 51), 2), (F(2, 9), 1), (F(5, 21), 4), (F(16, 51), 5)],
    [(F(5, 18), 3), (F(45, 136), 2), (F(7, 20), 0), (F(3, 8), 4), (F(35, 74), 5)],
    [(F(11, 32), 3), (F(7, 16), 1), (F(9, 20), 0), (F(1, 2), 4), (F(39, 64), 5)],
    [(F(35, 68), 2), (F(5, 9), 1), (F(10, 17), 0), (F(55, 64), 4), (F(15, 17), 5)],
    [(F(26, 35), 0), (F(3, 4), 1), (F(13, 17), 2), (F(117, 128), 3), (F(143, 152), 5)],
    [(F(369, 425), 0), (F(164, 185), 1), (F(123, 136), 2), (F(82, 85), 3), (F(369, 380), 4)],
]
PAPER_WITNESS_X = 15
PAPER_WITNESS_PIECES = 22
PAPER_WITNESS_SUBSEGMENTS = 36
PAPER_WITNESS_V = 27
PAPER_WITNESS_E = 48
PAPER_WITNESS_FACES = 23
PAPER_WITNESS_GRAPH_EDGES = 36
PAPER_WITNESS_DEGREES = [2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                         4, 4, 4, 4, 4, 4, 4, 5]
# The 22 pieces of Example 1 as printed in the paper, each named by the word
# e_0 e_1 ... e_5 in which e_c = 1 when the piece lies on the positive side of the
# directed chord c and e_c = 0 otherwise, in increasing lexicographic order.
PAPER_WITNESS_SIGNS = [
    '000000', '000001', '000011', '000111', '001011', '001111', '010011',
    '011011', '011111', '100000', '100001', '100011', '110000', '110001',
    '110011', '111000', '111001', '111011', '111100', '111101', '111110',
    '111111',
]
PAPER_WITNESS_ONEFLIP_EDGES = 36

# --- the degeneracy control of the paper (Section 6) ------------------------
PAPER_CONTROL_DIAMETER_T = [F(-5), F(-3), F(-2), F(1, 5), F(1, 3), F(1, 2)]
PAPER_CONTROL_DIAMETER_REJECTS = 1
PAPER_CONTROL_DIAMETER_MATCHING = ((0, 3), (1, 4), (2, 5))
PAPER_MATCHINGS_6 = 15

# ===========================================================================
# 1.  THE CHECK LEDGER
# ===========================================================================

_CHECKS = []
_T0 = time.time()


def check(name, ok, detail=''):
    _CHECKS.append((name, bool(ok), detail))
    print('%s %s%s' % ('PASS' if ok else 'FAIL', name, (' ' + detail) if detail else ''),
          flush=True)
    return bool(ok)


# ===========================================================================
# 2.  EXACT GEOMETRY.  The rational parametrisation of the unit circle.
# ===========================================================================

def circle_point(t):
    """t rational -> the exact rational point ((1-t^2)/(1+t^2), 2t/(1+t^2)).

    It satisfies x^2 + y^2 = 1 identically, and it is strictly monotone
    counterclockwise in t, so a strictly increasing t-vector lists boundary
    positions in counterclockwise order.
    """
    t = F(t)
    d = 1 + t * t
    return ((1 - t * t) / d, 2 * t / d)


def det2(ux, uy, vx, vy):
    return ux * vy - uy * vx


def interleaves(c, d):
    """Combinatorial crossing test: chords c=(a,b), d=(p,q) on the circle cross in the
    interior iff their four endpoints are distinct and exactly one of p, q lies strictly
    between a and b in the boundary order."""
    a, b = c
    p, q = d
    if len({a, b, p, q}) != 4:
        return False
    return (a < p < b) != (a < q < b)


def ccw_convex_position(pts):
    """True iff pts, in the given order, are the vertices of a convex polygon traversed
    counterclockwise: every triple i < j < k makes a left turn.  For points on a circle
    that is exactly the statement that the given order IS the counterclockwise boundary
    order, and it is decided by exact integer/rational determinants alone.
    """
    n = len(pts)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if det2(pts[j][0] - pts[i][0], pts[j][1] - pts[i][1],
                        pts[k][0] - pts[i][0], pts[k][1] - pts[i][1]) <= 0:
                    return False
    return True


def segments_cross(A, B, C, D):
    """Exact proper-intersection test for the open segments AB and CD."""
    def s(o, u, v):
        x = det2(u[0] - o[0], u[1] - o[1], v[0] - o[0], v[1] - o[1])
        return (x > 0) - (x < 0)
    return s(A, B, C) * s(A, B, D) < 0 and s(C, D, A) * s(C, D, B) < 0


def crossing_table(ts):
    """For a strictly increasing t-vector of length m, compute for EVERY chord (a,b)
    (a < b, boundary positions) and every chord d crossing it the exact parameter of the
    crossing along it, measured from its lower-indexed endpoint.

    -> (rank, sval) where sval[c][d] is that Fraction and rank[(c, d)] is its integer
    rank along c, EQUAL ranks meaning two crossings coincide, i.e. an interior multiple
    point.  Returns None on a degenerate point set (which cannot happen for distinct
    points on a circle, and is guarded rather than assumed).
    """
    pts = [circle_point(t) for t in ts]
    m = len(pts)
    chords = [(a, b) for a in range(m) for b in range(a + 1, m)]
    sval = defaultdict(dict)
    for i in range(len(chords)):
        for j in range(i + 1, len(chords)):
            c, d = chords[i], chords[j]
            if not interleaves(c, d):
                continue
            A, B = pts[c[0]], pts[c[1]]
            C, D = pts[d[0]], pts[d[1]]
            bx, by = B[0] - A[0], B[1] - A[1]
            dx, dy = D[0] - C[0], D[1] - C[1]
            den = det2(bx, by, dx, dy)
            if den == 0:
                return None
            s = det2(C[0] - A[0], C[1] - A[1], dx, dy) / den
            den2 = det2(dx, dy, bx, by)
            u = det2(A[0] - C[0], A[1] - C[1], bx, by) / den2
            if not (0 < s < 1 and 0 < u < 1):
                return None
            sval[c][d] = s
            sval[d][c] = u
    rank = {}
    for c in chords:
        items = sorted(sval[c].items(), key=lambda kv: kv[1])
        r = 0
        prev = None
        for (d, s) in items:
            if prev is not None and s != prev:
                r += 1
            rank[(c, d)] = r
            prev = s
    return rank, sval


def arrangement(match, rank):
    """(bw, orders) of the straight arrangement cut out by the chords `match`, or None if
    two of its crossings coincide -- which is exactly the interior multiple point the
    definition of A272906 forbids.

    bw[i]     = id of the chord with an endpoint at boundary position i;
    orders[c] = ids of the chords crossing c, in order along c from c's lower-indexed
                endpoint.  Chord ids are assigned by first appearance in bw.
    """
    m = 2 * len(match)
    ch = sorted(match)
    bw = [None] * m
    for cid, c in enumerate(ch):
        bw[c[0]] = cid
        bw[c[1]] = cid
    orders = []
    for cid, c in enumerate(ch):
        parts = []
        for did, d in enumerate(ch):
            if did == cid:
                continue
            if (c, d) in rank:
                parts.append((rank[(c, d)], did))
        rs = [r for (r, _) in parts]
        if len(set(rs)) != len(rs):
            return None                     # interior multiple point: REJECT
        parts.sort()
        orders.append(tuple(did for (_, did) in parts))
    return (tuple(bw), tuple(orders))


def perfect_matchings(k):
    """All (2k-1)!! perfect matchings of {0, ..., 2k-1}, each a sorted tuple of pairs."""
    def rec(rem):
        if not rem:
            yield ()
            return
        a = rem[0]
        for i in range(1, len(rem)):
            rest = rem[1:i] + rem[i + 1:]
            for tail in rec(rest):
                yield ((a, rem[i]),) + tail
    return [tuple(sorted(mm)) for mm in rec(tuple(range(2 * k)))]


def crossings_of(match):
    return sum(1 for i in range(len(match)) for j in range(i + 1, len(match))
               if interleaves(match[i], match[j]))


# ===========================================================================
# 3.  THE PLANAR MAP OF A TOPOLOGICAL ARRANGEMENT
# ===========================================================================

def build_map(bw, orders):
    """The half-edge map of the arrangement (bw, orders) inside the disk, or None if the
    data are not realizable as a topological chord arrangement.

    Vertices: ('b', i) for boundary position i, ('x', c, d) for the crossing of c and d.
    Edges: the m boundary arcs and the sub-segments of the chords.  A rotation system is
    read off the combinatorics alone, faces are traced, and Euler's formula is CHECKED.
    """
    m = len(bw)
    k = m // 2
    pos = defaultdict(list)
    for i, c in enumerate(bw):
        pos[c].append(i)
    for c in range(k):
        if len(pos[c]) != 2:
            return None
    for c in range(k):
        for d in orders[c]:
            if c not in orders[d]:
                return None
    seq = {}
    for c in range(k):
        s = [('b', pos[c][0])]
        for d in orders[c]:
            s.append(('x', min(c, d), max(c, d)))
        s.append(('b', pos[c][1]))
        if len(set(s)) != len(s):
            return None
        seq[c] = s

    edges = []

    def add_edge(u, v, kind, info):
        edges.append((u, v, kind, info))
        return len(edges) - 1

    arc_eid = [add_edge(('b', i), ('b', (i + 1) % m), 'arc', i) for i in range(m)]
    seg_eid = {}
    for c in range(k):
        s = seq[c]
        for t in range(len(s) - 1):
            seg_eid[(c, t)] = add_edge(s[t], s[t + 1], 'seg', (c, t))

    def rev(d):
        return (d[0], 1 - d[1])

    def tail(d):
        u, v, _, _ = edges[d[0]]
        return u if d[1] == 0 else v

    rot_next = {}

    def set_rot(v, cyc):
        if len(set(cyc)) != len(cyc):
            return False
        for j, d in enumerate(cyc):
            if tail(d) != v:
                return False
            rot_next[d] = cyc[(j + 1) % len(cyc)]
        return True

    for i in range(m):
        c = bw[i]
        s = seq[c]
        chord_dart = (seg_eid[(c, 0)], 0) if s[0] == ('b', i) \
            else (seg_eid[(c, len(s) - 2)], 1)
        if not set_rot(('b', i), [(arc_eid[i], 0), chord_dart,
                                 (arc_eid[(i - 1) % m], 1)]):
            return None

    for c in range(k):
        for d in orders[c]:
            if d < c:
                continue
            v = ('x', c, d)
            tc = seq[c].index(v)
            td = seq[d].index(v)
            c1 = (seg_eid[(c, tc - 1)], 1)
            c2 = (seg_eid[(c, tc)], 0)
            d1 = (seg_eid[(d, td - 1)], 1)
            d2 = (seg_eid[(d, td)], 0)
            a, b = pos[c]
            p, q = pos[d]
            if a < p < b:
                cyc = [c1, d1, c2, d2]
            elif a < q < b:
                cyc = [c1, d2, c2, d1]
            else:
                return None
            if not set_rot(v, cyc):
                return None

    darts = [(e, s) for e in range(len(edges)) for s in (0, 1)]
    seen = set()
    faces = []
    for d0 in darts:
        if d0 in seen:
            continue
        f = []
        d = d0
        while d not in seen:
            seen.add(d)
            f.append(d)
            nd = rot_next.get(rev(d))
            if nd is None:
                return None
            d = nd
        if d != d0:
            return None
        faces.append(f)
    X = sum(len(o) for o in orders) // 2
    V = m + X
    E = len(edges)
    if V - E + len(faces) != 2:
        return None
    outer, inner = None, []
    for f in faces:
        if len(f) == m and all(edges[e][2] == 'arc' for (e, _) in f):
            if outer is not None:
                return None
            outer = f
        else:
            inner.append(f)
    if outer is None or len(inner) != 1 + k + X:
        return None
    dartface = {}
    for i, f in enumerate(inner):
        for d in f:
            dartface[d] = i
    return dict(k=k, X=X, m=m, V=V, E=E, nfaces=len(faces), edges=edges, seq=seq,
                pos=pos, arc_eid=arc_eid, seg_eid=seg_eid, inner=inner,
                dartface=dartface)


def region_graph(mp):
    """A272906's graph: vertices = pieces, edges = pairs of pieces separated by a segment.
    Returns None if two pieces are separated by more than one segment or a segment fails
    to separate two distinct pieces (neither happens; both are guarded)."""
    n = len(mp['inner'])
    adj = [set() for _ in range(n)]
    seen = set()
    for (c, t), e in mp['seg_eid'].items():
        a = mp['dartface'].get((e, 0))
        b = mp['dartface'].get((e, 1))
        if a is None or b is None or a == b:
            return None
        key = (min(a, b), max(a, b))
        if key in seen:
            return None
        seen.add(key)
        adj[a].add(b)
        adj[b].add(a)
    return adj


def region_graph_vertex_rule(mp):
    """THE DELIBERATELY WRONG ADJACENCY RULE of the paper's control: two pieces adjacent
    when they share any VERTEX of the arrangement, rather than a segment."""
    n = len(mp['inner'])
    verts = [set() for _ in range(n)]
    for i, f in enumerate(mp['inner']):
        for (e, s) in f:
            u, v, _, _ = mp['edges'][e]
            verts[i].add(u)
            verts[i].add(v)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if verts[i] & verts[j]:
                adj[i].add(j)
                adj[j].add(i)
    return adj


# ===========================================================================
# 4.  CANONICAL FORMS
# ===========================================================================

def canon_arrangement(bw, orders):
    """Canonical form of a topological arrangement under the dihedral action on the
    boundary circle (m rotations x 2 reflections) together with chord relabelling."""
    m = len(bw)
    k = m // 2
    pos = defaultdict(list)
    for i, c in enumerate(bw):
        pos[c].append(i)
    best = None
    for refl in (1, -1):
        for r in range(m):
            oldpos = [(r + refl * j) % m for j in range(m)]
            lab = {}
            nb = []
            for op in oldpos:
                c = bw[op]
                if c not in lab:
                    lab[c] = len(lab)
                nb.append(lab[c])
            newpos = {}
            for op in oldpos:
                newpos.setdefault(lab[bw[op]], []).append(op)
            no = []
            for nc in range(k):
                op0 = newpos[nc][0]
                oc = bw[op0]
                o = list(orders[oc])
                if op0 != pos[oc][0]:
                    o = o[::-1]
                no.append(tuple(lab[x] for x in o))
            cand = (tuple(nb), tuple(no))
            if best is None or cand < best:
                best = cand
    return best


def canon_graph(adj):
    """Canonical adjacency-bitmask certificate of an undirected graph: equitable
    refinement, then a backtracking search over the individualisations consistent with
    the lexicographically least cell, keeping the least certificate.  Route 1 of two."""
    n = len(adj)
    col = [len(adj[i]) for i in range(n)]
    while True:
        sig = [(col[i], tuple(sorted(col[j] for j in adj[i]))) for i in range(n)]
        vals = {s: t for t, s in enumerate(sorted(set(sig)))}
        new = [vals[s] for s in sig]
        if new == col:
            break
        col = new
    best = [None]

    def rec(fixed, remaining):
        if not remaining:
            inv = {v: p for p, v in enumerate(fixed)}
            rows = []
            for v in fixed:
                r = 0
                for j in adj[v]:
                    r |= 1 << inv[j]
                rows.append(r)
            cert = tuple(rows)
            if best[0] is None or cert < best[0]:
                best[0] = cert
            return
        sig = {v: (col[v], tuple(sorted(p for p, u in enumerate(fixed) if u in adj[v])))
               for v in remaining}
        mn = min(sig.values())
        for v in remaining:
            if sig[v] == mn:
                rec(fixed + [v], [u for u in remaining if u != v])

    rec([], list(range(n)))
    return best[0]


def wl_colours(adj, rounds=6):
    n = len(adj)
    col = [len(adj[i]) for i in range(n)]
    for _ in range(rounds):
        sig = [(col[i], tuple(sorted(col[j] for j in adj[i]))) for i in range(n)]
        m = {s: t for t, s in enumerate(sorted(set(sig)))}
        new = [m[s] for s in sig]
        if new == col:
            break
        col = new
    return col


def graph_invariant(adj):
    """A cheap isomorphism-invariant: WL colours, degrees, triangle counts and the sorted
    distance profile of every vertex.  Route 2 of two, used to bucket before an exact
    isomorphism search."""
    n = len(adj)
    col = wl_colours(adj)
    prof = []
    for s in range(n):
        d = [-1] * n
        d[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if d[v] < 0:
                    d[v] = d[u] + 1
                    q.append(v)
        tri = sum(1 for a in adj[s] for b in adj[s] if a < b and b in adj[a])
        prof.append((col[s], len(adj[s]), tri, tuple(sorted(d)),
                     tuple(sorted(col[j] for j in adj[s]))))
    return (n, sum(len(a) for a in adj) // 2, tuple(sorted(col)), tuple(sorted(prof)))


def isomorphic(a1, a2):
    """Exact isomorphism search, WL-pruned; no heuristic acceptance."""
    n = len(a1)
    if n != len(a2) or sorted(len(s) for s in a1) != sorted(len(s) for s in a2):
        return False
    c1, c2 = wl_colours(a1), wl_colours(a2)
    if sorted(c1) != sorted(c2):
        return False
    order = sorted(range(n), key=lambda i: (-len(a1[i]), c1[i]))
    cand = {i: [j for j in range(n) if c2[j] == c1[i] and len(a2[j]) == len(a1[i])]
            for i in range(n)}
    mp = {}
    used = set()

    def rec(kk):
        if kk == n:
            return True
        u = order[kk]
        for v in cand[u]:
            if v in used:
                continue
            ok = True
            for w in a1[u]:
                if w in mp and mp[w] not in a2[v]:
                    ok = False
                    break
            if ok:
                for w in mp:
                    if w not in a1[u] and mp[w] in a2[v]:
                        ok = False
                        break
            if ok:
                mp[u] = v
                used.add(v)
                if rec(kk + 1):
                    return True
                del mp[u]
                used.discard(v)
        return False

    return rec(0)


# ===========================================================================
# 5.  THE INSERTION ENUMERATION (the census of Section 3)
# ===========================================================================

def insertions(bw, orders):
    """Every way to add one more topological chord: pick a starting boundary arc, walk a
    simple path of faces crossing each chord at most once, and end on a boundary arc.

    The only geometric facts used are that two straight chords meet at most once and that
    every face of a straight-chord arrangement in a convex disk is convex, so a straight
    chord meets each face in at most one segment.  Hence the family generated is a
    SUPERSET of the straight-realizable arrangements.
    """
    mp = build_map(bw, orders)
    if mp is None:
        return
    m, k = mp['m'], mp['k']
    dartface, edges = mp['dartface'], mp['edges']
    face_arcs = defaultdict(list)
    arc_face = {}
    for i in range(m):
        e = mp['arc_eid'][i]
        f = dartface.get((e, 0))
        if f is None:
            f = dartface.get((e, 1))
        arc_face[i] = f
        face_arcs[f].append(i)
    face_segs = defaultdict(list)
    for (c, t), e in mp['seg_eid'].items():
        a = dartface.get((e, 0))
        b = dartface.get((e, 1))
        face_segs[a].append((c, t, b))
        face_segs[b].append((c, t, a))

    out = []

    def walk(f, visitedf, crossed, path, start_arc):
        for j in face_arcs[f]:
            out.append((start_arc, j, tuple(path)))
        for (c, t, g) in face_segs[f]:
            if c in crossed or g in visitedf:
                continue
            walk(g, visitedf | {g}, crossed | {c}, path + [(c, t)], start_arc)

    for i in range(m):
        f = arc_face[i]
        walk(f, {f}, frozenset(), [], i)

    seen = set()
    for (i, j, path) in out:
        newid = k
        if i == j:
            nb = []
            for p in range(m):
                nb.append(bw[p])
                if p == i:
                    nb.append(newid)
                    nb.append(newid)
            first_is_start = True
        else:
            nb = []
            for p in range(m):
                nb.append(bw[p])
                if p == i or p == j:
                    nb.append(newid)
            first_is_start = (i < j)
        no = [list(orders[c]) for c in range(k)]
        crossed_in_order = [c for (c, _) in path]
        for (c, t) in path:
            no[c].insert(t, newid)
        no.append(crossed_in_order if first_is_start else crossed_in_order[::-1])
        key = (tuple(nb), tuple(tuple(x) for x in no))
        if key in seen:
            continue
        seen.add(key)
        yield key


def census(nmax):
    """-> {n: {canonical form: (bw, orders)}} for n = 1 .. nmax."""
    start = ((0, 0), ((),))
    levels = {1: {canon_arrangement(*start): start}}
    level = levels[1]
    for n in range(2, nmax + 1):
        nxt = {}
        for (bw, orders) in level.values():
            for (bw2, o2) in insertions(bw, orders):
                if build_map(bw2, o2) is None:
                    continue
                cf = canon_arrangement(bw2, o2)
                if cf not in nxt:
                    nxt[cf] = (bw2, o2)
        level = nxt
        levels[n] = level
    return levels


# ===========================================================================
# 6.  INDEPENDENT ENUMERATION OF THE TREES ON 7 NODES (the X = 0 anchor)
# ===========================================================================

def unlabelled_trees(n):
    """Canonical forms of all unlabelled trees on n nodes, obtained by decoding every one
    of the n^(n-2) Prufer sequences.  Cayley's formula says that enumerates every LABELLED
    tree exactly once, so canonicalising the results enumerates every unlabelled tree.
    Used to check the X = 0 shard against A000055(n) without trusting the census."""
    if n == 1:
        return {canon_graph([set()])}
    if n == 2:
        return {canon_graph([{1}, {0}])}
    out = set()
    nlab = 0
    for seq in itertools.product(range(n), repeat=n - 2):
        deg = [1] * n
        for x in seq:
            deg[x] += 1
        adj = [set() for _ in range(n)]
        d = deg[:]
        for x in seq:
            leaf = min(i for i in range(n) if d[i] == 1)
            adj[leaf].add(x)
            adj[x].add(leaf)
            d[leaf] = 0
            d[x] -= 1
        rest = [i for i in range(n) if d[i] == 1]
        assert len(rest) == 2
        adj[rest[0]].add(rest[1])
        adj[rest[1]].add(rest[0])
        assert sum(len(a) for a in adj) // 2 == n - 1
        nlab += 1
        out.add(canon_graph(adj))
    assert nlab == n ** (n - 2), (nlab, n)
    return out


# ===========================================================================
# 7.  THE PROGRAM
# ===========================================================================

def main():
    print('verify.py -- "An Exhaustive Census Confirms a(6) = 1814 for Topologically '
          'Distinct Pizza Slicings"')
    print('exact integer / Fraction arithmetic only; standard library only; '
          'no external data file')
    print()

    # ------------------------------------------------------------------
    print('--- PART 1: the census (the upper bound a(6) <= 1814) ---')
    lv = census(6)
    print('    [%.1fs] insertion enumeration to n=6 finished' % (time.time() - _T0),
          flush=True)

    # types and region graphs, level by level
    graphs_by_n = {}
    maps_by_n = {}
    for n in range(2, 7):
        gs = []
        ms = []
        for (bw, orders) in lv[n].values():
            mp = build_map(bw, orders)
            g = region_graph(mp)
            gs.append((mp['X'], g, bw, orders))
            ms.append(mp)
        graphs_by_n[n] = gs
        maps_by_n[n] = ms

    check('census-region-graph-well-defined',
          all(g is not None for n in range(2, 7) for (_, g, _, _) in graphs_by_n[n]),
          '[every segment separates two distinct pieces, and no two pieces are '
          'separated by two segments, at every n<=6]')

    check('census-euler',
          all(mp['V'] - mp['E'] + mp['nfaces'] == 2
              for n in range(2, 7) for mp in maps_by_n[n]),
          '[V-E+F=2 on all %d traced maps]'
          % sum(len(maps_by_n[n]) for n in range(2, 7)))

    check('census-pieces-equal-7-plus-X',
          all(len(mp['inner']) == 1 + mp['k'] + mp['X'] for mp in maps_by_n[6]),
          '[all 1824 types at n=6 have 1+6+X = 7+X pieces]')

    # the published confirmed prefix, from the same code
    small = {}
    for n in range(2, 6):
        cls = {}
        for (X, g, _, _) in graphs_by_n[n]:
            cls[(X, canon_graph(g))] = 1
        small[n] = len(cls)
    for n in (2, 3, 4, 5):
        check('census-a%d' % n, small[n] == PAPER_A272906_PREFIX[n],
              '[%d distinct region graphs; A272906(%d) = %d, a confirmed term]'
              % (small[n], n, PAPER_A272906_PREFIX[n]))

    types6 = graphs_by_n[6]
    check('census-types-n6', len(types6) == PAPER_TYPES_6,
          '[%d topological arrangement types at n=6; paper says %d]'
          % (len(types6), PAPER_TYPES_6))

    tperx = [sum(1 for (X, _, _, _) in types6 if X == x) for x in range(16)]
    check('census-types-perX', tperx == PAPER_TYPES_PERX, '[%s]' % tperx)

    # ---- route 1: a canonical form for every region graph
    canon_of = [canon_graph(g) for (_, g, _, _) in types6]
    canon = {}
    for i, (X, _, _, _) in enumerate(types6):
        canon.setdefault((X, canon_of[i]), []).append(i)
    print('    [%.1fs] route 1 (canonical form) finished' % (time.time() - _T0),
          flush=True)
    check('census-classes-n6-canonical-form', len(canon) == PAPER_A272906_6,
          '[%d distinct region graphs among the %d types; paper says %d]'
          % (len(canon), len(types6), PAPER_A272906_6))

    # ---- route 2: invariant buckets + exact isomorphism search
    buckets = defaultdict(list)
    for i, (X, g, _, _) in enumerate(types6):
        buckets[(X, graph_invariant(g))].append(i)
    reps = []
    merged = 0
    for key, idxs in sorted(buckets.items(), key=lambda kv: kv[1][0]):
        local = []
        for i in idxs:
            if not any(isomorphic(types6[i][1], types6[r][1]) for r in local):
                local.append(i)
            else:
                merged += 1
        reps.extend(local)
    print('    [%.1fs] route 2 (invariant + exact isomorphism) finished'
          % (time.time() - _T0), flush=True)
    check('census-classes-n6-invariant-plus-iso', len(reps) == PAPER_A272906_6,
          '[%d buckets, %d non-singleton, %d exact isomorphisms found; %d classes]'
          % (len(buckets), sum(1 for b in buckets.values() if len(b) > 1),
             merged, len(reps)))
    check('census-two-routes-agree', len(reps) == len(canon),
          '[%d = %d; the two isomorphism routes give the same class count]'
          % (len(reps), len(canon)))

    cperx = [sum(1 for i in reps if types6[i][0] == x) for x in range(16)]
    check('census-classes-perX', cperx == PAPER_CLASSES_PERX, '[%s]' % cperx)
    check('census-classes-sum', sum(cperx) == PAPER_A272906_6,
          '[%d = %d]' % (sum(cperx), PAPER_A272906_6))

    check('census-collapse-count', len(types6) - len(reps) == PAPER_COLLAPSES,
          '[%d types collapse onto an earlier one; paper says %d]'
          % (len(types6) - len(reps), PAPER_COLLAPSES))
    shards = {x: tperx[x] - cperx[x] for x in range(16) if tperx[x] != cperx[x]}
    check('census-collapse-shards', shards == PAPER_COLLAPSE_SHARDS,
          '[%s; every collapse has X<=4 and there is none for X>=5]' % shards)

    check('census-upper-bound', len(reps) == PAPER_A272906_6,
          '[a(6) <= %d, since the census family contains every straight-realizable '
          'type] ' % PAPER_A272906_6)

    # ------------------------------------------------------------------
    print()
    print('--- PART 2: controls on the census ---')

    # X=0 shard against an independently enumerated A000055(7)
    zero = [i for i in reps if types6[i][0] == 0]
    zero_are_trees = all(sum(len(s) for s in types6[i][1]) // 2 == 6
                         and len(types6[i][1]) == 7 for i in zero)
    trees7 = unlabelled_trees(7)
    zero_canon = set(canon_of[i] for i in zero)
    check('control-X0-shard-is-all-trees-on-7-nodes',
          zero_are_trees and len(trees7) == PAPER_A000055_7 and zero_canon == trees7,
          '[%d classes with X=0; an independent Prufer enumeration gives %d unlabelled '
          'trees on 7 nodes = A000055(7) = %d, and the two SETS of canonical forms are '
          'identical]' % (len(zero), len(trees7), PAPER_A000055_7))

    check('control-X15-shard-equals-A090338-6', cperx[15] == PAPER_A090338_6,
          '[%d classes with all 15 pairs crossing; A090338(6) = %d simple arrangements '
          'of 6 lines]' % (cperx[15], PAPER_A090338_6))

    wrong = set()
    for mp in maps_by_n[4]:
        wrong.add(canon_graph(region_graph_vertex_rule(mp)))
    check('control-wrong-adjacency-rule-undercounts',
          len(wrong) == PAPER_WRONG_ADJACENCY_N4 != PAPER_A272906_PREFIX[4],
          '[the vertex-sharing rule gives %d graphs at n=4, not %d: reproducing 19 is '
          'evidence the segment rule is the one implemented]'
          % (len(wrong), PAPER_A272906_PREFIX[4]))

    # the degeneracy detector must FIRE, exactly once, on three diameters
    ts = PAPER_CONTROL_DIAMETER_T
    pts = [circle_point(t) for t in ts]
    diam_ok = all(pts[i][0] + pts[i + 3][0] == 0 and pts[i][1] + pts[i + 3][1] == 0
                  for i in range(3))
    ms6 = perfect_matchings(3)
    rk = crossing_table(ts)
    rejected = [M for M in ms6 if rk is not None and arrangement(M, rk[0]) is None]
    check('control-degeneracy-detector-fires-exactly-once',
          diam_ok and len(ms6) == PAPER_MATCHINGS_6
          and len(rejected) == PAPER_CONTROL_DIAMETER_REJECTS
          and rejected[0] == PAPER_CONTROL_DIAMETER_MATCHING,
          '[t = -5,-3,-2,1/5,1/3,1/2 gives P(t)+P(-1/t) = (0,0) exactly, so %s is three '
          'diameters through the origin; of the %d matchings exactly %d is rejected and '
          'it is that one]'
          % (str(PAPER_CONTROL_DIAMETER_MATCHING), PAPER_MATCHINGS_6, len(rejected)))

    # ------------------------------------------------------------------
    print()
    print('--- PART 3: the witness exhibited in Example 1 ---')
    wt = [F(x) for x in PAPER_WITNESS_T]
    check('witness-t-strictly-increasing',
          all(wt[i] < wt[i + 1] for i in range(11)) and len(set(wt)) == 12,
          '[12 distinct rationals in increasing order]')
    wp = [circle_point(t) for t in wt]
    check('witness-points-printed-correctly',
          wp == PAPER_WITNESS_POINTS,
          '[all 12 points agree with the paper, exactly]')
    check('witness-points-on-unit-circle',
          all(p[0] * p[0] + p[1] * p[1] == 1 for p in wp),
          '[x^2+y^2 = 1 exactly at all 12 boundary points]')
    check('witness-boundary-order-is-t-order', ccw_convex_position(wp),
          '[all C(12,3) = 220 triples i<j<k are positively oriented, so the 12 points are '
          'in convex position counterclockwise and the sorted t-order IS the boundary '
          'order]')

    wrk, wsv = crossing_table(wt)
    WM = tuple(PAPER_WITNESS_MATCHING)
    check('witness-crossing-number', crossings_of(WM) == PAPER_WITNESS_X,
          '[%d interleaving pairs among the 6 chords; X = %d, the maximum C(6,2)]'
          % (crossings_of(WM), PAPER_WITNESS_X))
    # the printed s-values, and the crossing orders they induce
    derived = []
    for cid, c in enumerate(sorted(WM)):
        row = []
        for did, d in enumerate(sorted(WM)):
            if did != cid and d in wsv[c]:
                row.append((wsv[c][d], did))
        row.sort()
        derived.append(row)
    check('witness-crossing-parameters-printed-correctly',
          derived == PAPER_WITNESS_S,
          '[all 30 crossing parameters, and the order they induce along each chord, '
          'agree with the paper exactly]')
    check('witness-crossings-strictly-interior',
          all(0 < s < 1 for row in derived for (s, _) in row),
          '[all 30 parameters lie strictly in (0,1), so every crossing is strictly '
          'inside the disk and off the boundary]')
    check('witness-no-interior-multiple-point',
          all(len(set(s for (s, _) in row)) == len(row) for row in derived),
          '[on each chord the 5 crossing parameters are pairwise distinct; a triple point '
          'would put two of them at the same value]')

    wa = arrangement(WM, wrk)
    check('witness-arrangement-accepted', wa is not None,
          '[the witness survives the triple-point test]')
    wmp = build_map(*wa)
    check('witness-map-traced', wmp is not None and wmp['X'] == PAPER_WITNESS_X,
          '[the planar map is consistent and has X = %d]' % PAPER_WITNESS_X)
    nseg = len(wmp['seg_eid'])
    check('witness-subsegment-count', nseg == PAPER_WITNESS_SUBSEGMENTS
          and nseg == 6 + 2 * PAPER_WITNESS_X,
          '[%d chord sub-segments = 6 + 2*%d]' % (nseg, PAPER_WITNESS_X))
    check('witness-euler-bookkeeping',
          (wmp['V'], wmp['E'], wmp['nfaces']) ==
          (PAPER_WITNESS_V, PAPER_WITNESS_E, PAPER_WITNESS_FACES),
          '[V = %d = 12+X, E = %d = 12 arcs + %d sub-segments, F = %d = 2-V+E]'
          % (wmp['V'], wmp['E'], nseg, wmp['nfaces']))
    check('witness-piece-count',
          len(wmp['inner']) == PAPER_WITNESS_PIECES == 7 + PAPER_WITNESS_X,
          '[%d pieces = F - 1 = 7 + X]' % len(wmp['inner']))
    wg = region_graph(wmp)
    check('witness-region-graph-simple-and-connected',
          wg is not None and sum(len(s) for s in wg) // 2 == PAPER_WITNESS_GRAPH_EDGES
          and _connected(wg),
          '[%d vertices, %d edges, connected, no two pieces separated by two segments]'
          % (len(wg), sum(len(s) for s in wg) // 2))
    check('witness-degree-sequence',
          sorted(len(s) for s in wg) == PAPER_WITNESS_DEGREES,
          '[%s]' % sorted(len(s) for s in wg))
    # --- a THIRD, geometric route to the same region graph: name each piece by its sign
    # --- vector with respect to the six directed chords, and read the adjacency off as a
    # --- single sign flip.  Nothing here consults the map's rotation system.
    wch = sorted(WM)

    def vertex_coord(vert):
        if vert[0] == 'b':
            return wp[vert[1]]
        c, d = vert[1], vert[2]
        A, B = wp[wch[c][0]], wp[wch[c][1]]
        s = wsv[wch[c]][wch[d]]
        return (A[0] + s * (B[0] - A[0]), A[1] + s * (B[1] - A[1]))

    def side(p, c):
        A, B = wp[wch[c][0]], wp[wch[c][1]]
        x = det2(B[0] - A[0], B[1] - A[1], p[0] - A[0], p[1] - A[1])
        return (x > 0) - (x < 0)

    face_sign = []
    for f in wmp['inner']:
        vs = set()
        for (e, _) in f:
            u, w, _, _ = wmp['edges'][e]
            vs.add(u)
            vs.add(w)
        cs = [vertex_coord(x) for x in sorted(vs)]
        cx = sum(p[0] for p in cs) / len(cs)
        cy = sum(p[1] for p in cs) / len(cs)
        face_sign.append(tuple(side((cx, cy), c) for c in range(6)))
    as_str = [''.join('1' if x > 0 else ('0' if x < 0 else '?') for x in s)
              for s in sorted(face_sign)]
    check('witness-sign-vectors-name-the-pieces',
          len(set(face_sign)) == PAPER_WITNESS_PIECES
          and all(0 not in s for s in face_sign) and as_str == PAPER_WITNESS_SIGNS,
          '[the centroid of each piece\'s corners is strictly off all six chord lines, '
          'the %d sign words are pairwise distinct, and they are exactly the %d printed '
          'in the paper]' % (len(face_sign), len(PAPER_WITNESS_SIGNS)))
    order = sorted(range(len(face_sign)), key=lambda i: face_sign[i])
    newid = {old: k for k, old in enumerate(order)}
    got = set()
    for i in range(len(wg)):
        for j in wg[i]:
            got.add((min(newid[i], newid[j]), max(newid[i], newid[j])))
    S = [face_sign[order[k]] for k in range(len(order))]
    flips = set((i, j) for i in range(len(S)) for j in range(i + 1, len(S))
                if sum(1 for k in range(6) if S[i][k] != S[j][k]) == 1)
    check('witness-region-graph-is-the-one-sign-flip-graph',
          got == flips and len(flips) == PAPER_WITNESS_ONEFLIP_EDGES,
          '[the %d edges traced in the planar map are exactly the %d pairs of sign vectors '
          'differing in one coordinate -- a third, purely geometric derivation of the same '
          'region graph, so a reader can read Table 2 off by hand]'
          % (len(got), len(flips)))

    census_types = {}
    for i, (X, g, bw, orders) in enumerate(types6):
        census_types[canon_arrangement(bw, orders)] = i
    wcf = canon_arrangement(*wa)
    check('witness-type-is-in-the-census', wcf in census_types,
          '[the witness realizes census type %d, in the X = %d shard]'
          % (census_types.get(wcf, -1), PAPER_WITNESS_X))

    # ------------------------------------------------------------------
    print()
    print('--- PART 4: Table 1 and the sweep (the lower bound a(6) >= 1814) ---')
    check('table-shape', len(PAPER_TABLE1) == 47
          and all(len(r) == 12 and all(isinstance(x, int) for x in r)
                  and all(r[i] < r[i + 1] for i in range(11)) for r in PAPER_TABLE1),
          '[47 rows, each 12 strictly increasing integers]')
    allpts = [[circle_point(F(x)) for x in r] for r in PAPER_TABLE1]
    check('table-points-on-unit-circle',
          all(p[0] * p[0] + p[1] * p[1] == 1 for row in allpts for p in row),
          '[all 47*12 = %d boundary points satisfy x^2+y^2 = 1 exactly]'
          % (47 * 12))
    check('table-boundary-order-is-t-order',
          all(ccw_convex_position(row) for row in allpts),
          '[all 47 * 220 = %d oriented triples are positive, so every row is in convex '
          'position counterclockwise and its t-order is its boundary order]' % (47 * 220))

    # combinatorial crossing test == exact geometric crossing test
    nprs = 0
    agree = True
    for row in allpts:
        chords = [(a, b) for a in range(12) for b in range(a + 1, 12)]
        for i in range(len(chords)):
            for j in range(i + 1, len(chords)):
                c, d = chords[i], chords[j]
                if len({c[0], c[1], d[0], d[1]}) != 4:
                    continue
                nprs += 1
                if interleaves(c, d) != segments_cross(row[c[0]], row[c[1]],
                                                       row[d[0]], row[d[1]]):
                    agree = False
    check('table-interleaving-equals-geometry', agree,
          '[on all %d chord pairs with four distinct endpoints, the combinatorial '
          'interleaving test and the exact segment-intersection test agree]' % nprs)
    print('    [%.1fs] geometry cross-check finished' % (time.time() - _T0), flush=True)

    MS = perfect_matchings(6)
    check('matchings-count', len(MS) == PAPER_MATCHINGS_12,
          '[%d perfect matchings of 12 boundary positions = 11!! ]' % len(MS))

    realized = {}
    outside = 0
    admissible = 0
    pairs = 0
    for k, r in enumerate(PAPER_TABLE1):
        rk = crossing_table([F(x) for x in r])
        if rk is None:
            outside += 1
            continue
        for M in MS:
            pairs += 1
            a = arrangement(M, rk[0])
            if a is None:
                continue
            admissible += 1
            cf = canon_arrangement(*a)
            if cf not in census_types:
                outside += 1
                continue
            realized.setdefault(cf, (k, M))
        if (k + 1) % 12 == 0:
            print('    [%.1fs] swept %d of 47 rows; %d types realized'
                  % (time.time() - _T0, k + 1, len(realized)), flush=True)

    check('sweep-pair-count', pairs == PAPER_SWEEP_PAIRS,
          '[47 * %d = %d (point set, matching) pairs examined]'
          % (PAPER_MATCHINGS_12, pairs))
    check('sweep-admissible-count', admissible == PAPER_SWEEP_ADMISSIBLE
          and pairs - admissible == PAPER_SWEEP_REJECTED,
          '[%d pairs survive the triple-point test; %d are rejected for an interior '
          'multiple point]' % (admissible, pairs - admissible))
    check('sweep-nothing-outside-the-census', outside == PAPER_SWEEP_OUTSIDE_CENSUS,
          '[%d of the %d admissible straight arrangements has a type outside the census '
          '-- ONE would have proved a(6) >= 1815]' % (outside, admissible))
    check('sweep-realizes-every-census-type', len(realized) == len(types6),
          '[%d of %d topological types are realized by straight chords with rational '
          'endpoints]' % (len(realized), len(types6)))
    rperx = [sum(1 for cf in realized if types6[census_types[cf]][0] == x)
             for x in range(16)]
    check('sweep-realized-perX', rperx == PAPER_TYPES_PERX, '[%s]' % rperx)
    realized_classes = set()
    for cf in realized:
        i = census_types[cf]
        realized_classes.add((types6[i][0], canon_of[i]))
    check('sweep-realizes-every-region-graph-class',
          len(realized_classes) == PAPER_A272906_6,
          '[%d region-graph classes realized; a(6) >= %d]'
          % (len(realized_classes), PAPER_A272906_6))
    check('lower-bound-and-upper-bound-meet',
          len(realized_classes) == len(reps) == PAPER_A272906_6,
          '[a(6) >= %d by exhibited exact objects and a(6) <= %d by the census, so '
          'a(6) = %d]' % (len(realized_classes), len(reps), PAPER_A272906_6))
    check('corollary-no-non-stretchable-type-at-n6',
          len(realized) == len(types6) == PAPER_TYPES_6,
          '[all %d topological types, not merely all %d graph classes, are realized by '
          'straight rational chords]' % (PAPER_TYPES_6, PAPER_A272906_6))

    # mirror control: t -> -t must not change the set of types a row realizes
    mirror_ok = True
    for r in PAPER_TABLE1:
        rk = crossing_table([F(x) for x in r])[0]
        mk = crossing_table([F(-x) for x in reversed(r)])[0]
        for M in (tuple((i, i + 6) for i in range(6)),):
            a = arrangement(M, rk)
            b = arrangement(M, mk)
            if (a is None) != (b is None):
                mirror_ok = False
            elif a is not None and canon_arrangement(*a) != canon_arrangement(*b):
                mirror_ok = False
    row1 = crossing_table([F(x) for x in PAPER_TABLE1[0]])[0]
    row1m = crossing_table([F(-x) for x in reversed(PAPER_TABLE1[0])])[0]
    s1 = set()
    s2 = set()
    for M in MS:
        a = arrangement(M, row1)
        if a is not None:
            s1.add(canon_arrangement(*a))
        b = arrangement(M, row1m)
        if b is not None:
            s2.add(canon_arrangement(*b))
    check('control-mirror-symmetry', mirror_ok and s1 == s2 and len(s1) > 1000,
          '[reflecting t -> -t with i -> 11-i is an orientation-reversing symmetry of '
          'the disk: it fixes the canonical type of all 47 all-crossing arrangements, '
          'and row 1 realizes the same set of %d types either way]' % len(s1))

    # ------------------------------------------------------------------
    print()
    print('NOT RE-RUN -- what this program does NOT establish:')
    print('  * NOT RE-RUN: the SUPERSET LEMMA itself.  That two straight chords meet at '
          'most once and that a straight chord meets each convex face in at most one '
          'segment are proved in Section 3 of the paper by hand; this program re-runs '
          'the enumeration those facts license, not the facts.')
    print('  * NOT RE-RUN: any n >= 7.  A272906(7) and beyond are untouched, as are '
          'A241600 (known for n <= 7) and A090338 (known for n <= 9).')
    print('  * NOT RE-RUN: the sister sequence A273280 (chords of a SQUARE; 1, 1, 2, 5, '
          '19, 129, 1806, Giovanni Resta, May 2016), which carries the identical '
          '"unconfirmed" caveat and is not attempted here even though the method '
          'applies.')
    print('  * NOT RE-RUN: the identification of the 43 classes of the X = 15 shard with '
          "the 43 chirotopes of Christ's database of simple arrangements of 6 lines.  "
          'Only the two COUNTS are compared, 43 = A090338(6).')
    print('  * NOT RE-RUN: Jon Hart\'s May 2016 guided random trials, which produced the '
          'value 1814 in the first place.  The value is his; what is re-derived here is '
          'the exhaustive upper bound that makes it a theorem.')
    print('  * NOT RE-RUN: any claim that Table 1 is minimal.  47 point sets suffice; '
          'nothing here says fewer do not.')
    print()

    fails = [n for (n, ok, _) in _CHECKS if not ok]
    print('checks: %d run, %d passed, %d failed, %.1fs'
          % (len(_CHECKS), len(_CHECKS) - len(fails), len(fails), time.time() - _T0))
    if fails:
        for n in fails:
            print('FAILED: %s' % n)
        print('VERDICT: %d OF %d CHECKS FAILED' % (len(fails), len(_CHECKS)))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % len(_CHECKS))
    return 0


def _connected(adj):
    n = len(adj)
    seen = {0}
    q = deque([0])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                q.append(v)
    return len(seen) == n


if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    sys.exit(main())

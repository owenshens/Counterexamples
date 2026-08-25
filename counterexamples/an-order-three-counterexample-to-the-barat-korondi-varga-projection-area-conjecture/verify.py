#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- referee's verification program for

    "An Order-Three Counterexample to the Barat--Korondi--Varga
     Projection-Area Conjecture"

Standard library only.  Exact integer arithmetic throughout; no floating point
is used in any decision.  Single process, well under one minute.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER  (inputs -- these are transcribed, not checked)
--------------------------------------------------------------------------
  P_B          the exhibited base position B = {111,113,122,131,223,311,
               322,331,333} in [3]^3                       (proof, line 1)
  P_PROJ       the three projection sets pi_1(B), pi_2(B), pi_3(B) as
               printed in the displayed align* environment
  P_TABLE      the 18-row build-up certificate: for each k the added vertex
               v_k and the printed neighbourhood N_{S_{k-1}}(v_k)
  P_PROJ_SIZES (6,7,7)      P_AREA 20      P_BOUND 23
  P_FORMULA    the conjectured lower bound n^2 + 6n - 4, transcribed from the
               displayed inequality of the note.  The note quotes that
               formulation at second hand, as the one printed as Conjecture 1
               of Barat-Wanless, and states that it cannot confirm what
               number the inequality carries in Barat-Korondi-Varga; this
               program therefore attaches no number to it either.
  P_RULE       "solution" = independent B with |B| = n^2 from which [n]^3 is
               reachable by adding vertices one at a time, each added vertex
               having EXACTLY three neighbours already present
  P_NSOL       116 order-three solutions      (proof, last sentence; and the
               Exact-verification paragraph)
  P_NCLASS     7 isometry classes                 (paper, from Barat-Wanless)
  P_CLASS_AREA the seven class totals 27,27,25,23,21,20,21
  P_NIND       4224 independent 9-subsets of [3]^3   (proof, last sentence;
               and the Exact-verification paragraph)

--------------------------------------------------------------------------
WHAT IS DERIVED HERE  (the checks; nothing below is copied from the paper)
--------------------------------------------------------------------------
  * the grid graph on [n]^3 and its adjacency, built from the definition,
    and its degree profile (8 corners, 12 edge midpoints, 6 face centres,
    1 centre, 54 edges) for n = 3;
  * that the addition rule really is EXACTLY three and not at least three:
    [3]^3 minus its centre is shown unreachable, minus a corner reachable,
    minus an edge midpoint unreachable, and minus its centre reachable
    under the relaxed rule -- so the two readings are distinguished by the
    code, not merely by a comment;
  * that P_B is well formed, has 9 = 3^2 distinct elements of [3]^3, and is
    an independent set -- every one of the 36 pairs is tested;
  * pi_1(B), pi_2(B), pi_3(B) recomputed from P_B and compared elementwise
    with P_PROJ; their sizes and A(B) = sum |pi_i(B)| recomputed;
  * the conjectured bound n^2+6n-4 evaluated at n=3, and the strict
    inequality A(B) < that value -- the load-bearing refutation step;
  * every row of P_TABLE re-derived: the true neighbourhood of v_k inside
    S_{k-1} is computed from the adjacency and compared with the printed
    set, its cardinality is required to be exactly 3, the v_k are required
    to be distinct and to exhaust [3]^3 \ B, and S_18 is required to be
    all of [3]^3;
  * an INDEPENDENT search (backtracking over addition orders, not using
    P_TABLE at all) that B is a solution;
  * the complete order-three census: all independent 9-subsets of [3]^3 are
    enumerated, each is tested for solution-hood by exact backtracking
    search, the 48-element cube isometry group is built and verified, the
    solutions are partitioned into isometry classes, A is verified invariant
    on those classes, and s_3(3) = min A over solutions is computed;
  * the count of independent 9-subsets recomputed by a SECOND, unrelated
    algorithm (a layer transfer-matrix DP), and the same DP used to size
    the n = 4 census honestly in the printed SCOPE line;
  * the same census for n = 1 and n = 2, where the conjecture is confirmed
    to hold with equality -- a control showing the census machinery does
    not systematically undershoot.

Every P_* value printed in the note is confronted with a quantity computed
from the definitions.  That includes P_NIND = 4224 and P_NSOL = 116: both are
printed in the note -- in the closing sentence of the proof ("of the 4224
independent 9-subsets of [3]^3, exactly 116 are solutions") and again in its
Exact-verification paragraph -- so the two checks that compare them with the
census are tests of the paper like the rest.  They are corroborated further
by a second counting algorithm for 4224 (a transfer-matrix DP) and by the
paper's own class count for 116.

The mutation results quoted in the next two paragraphs were obtained while
this program was being written, by editing copies of it and of the data above
and re-running them.  They are NOT part of the transcript this program
produces: no mutation harness, mutation list or mutation log accompanies it,
and nothing below re-runs them.  The closing SCOPE line says so as well.

The program was mutation-tested against 31 deliberate corruptions of both
the paper's data and its own code (a corrupted vertex of B, a corrupted
table row, a swapped table order, an illegal table order, B replaced by a
non-solution, B replaced by a still-independent 9-set, each quoted constant
perturbed, and sabotages of the geometry, the addition rule, the
projections, the enumeration and the isometry group).  Every one is reported
as a FAIL with exit status 1 except one that is provably equivalent on this
input: interchanging the first two table rows, whose printed neighbourhoods
both lie inside B, so their relative order is immaterial.

An adversarial review then found three mutations that this program used to
report as a PASS, and all three are now caught:
  (i)   deleting the call to any whole section -- section_refutation,
        section_machinery, section_table or section_census3 -- left every
        remaining check passing and printed "ALL 30 CHECKS PASS" with exit
        status 0.  The structural guard now asserts the check count against
        EXPECTED_CHECKS instead of merely reporting it.
  (ii)  a repeated element in a printed projection set was invisible to the
        set-equality comparison, so a display whose entries do not match its
        own stated cardinality passed.  The printed lists are now required to
        be repetition-free and as long as the recomputed images.
  (iii) weakening the load-bearing comparison A(B) < bound to A(B) != bound
        still passed on this input, though it would falsely certify a
        solution of area 25.  The comparison is now a named predicate whose
        direction and strictness are exercised on bound + 2 and on bound.
"""

import sys
from itertools import combinations, permutations, product

# =====================================================================
# BLOCK 1.  VERBATIM PAPER INPUTS.  Corrupt anything here and the checks
# below must report FAIL.
# =====================================================================

# B, exactly as printed in the proof (strings "ijk").
P_B = ["111", "113", "122", "131", "223", "311", "322", "331", "333"]

# The three projection sets as printed in the align* display.
P_PROJ = {
    1: [(1, 1), (1, 3), (2, 2), (2, 3), (3, 1), (3, 3)],
    2: [(1, 1), (1, 2), (1, 3), (2, 3), (3, 1), (3, 2), (3, 3)],
    3: [(1, 1), (1, 2), (1, 3), (2, 2), (3, 1), (3, 2), (3, 3)],
}
P_PROJ_SIZES = (6, 7, 7)     # "projection sizes 6,7,7"
P_AREA = 20                  # "A(B) = 6+7+7 = 20"
P_BOUND = 23                 # "20 < 23 = 3^2 + 6*3 - 4"
P_N = 3                      # order three

# The 18-row build-up table: (k, v_k, printed N_{S_{k-1}}(v_k)).
P_TABLE = [
    (1,  "112", ["111", "113", "122"]),
    (2,  "121", ["111", "122", "131"]),
    (3,  "123", ["113", "122", "223"]),
    (4,  "222", ["122", "223", "322"]),
    (5,  "321", ["311", "322", "331"]),
    (6,  "323", ["223", "322", "333"]),
    (7,  "332", ["322", "331", "333"]),
    (8,  "221", ["121", "222", "321"]),
    (9,  "231", ["131", "221", "331"]),
    (10, "232", ["222", "231", "332"]),
    (11, "233", ["223", "232", "333"]),
    (12, "132", ["122", "131", "232"]),
    (13, "133", ["123", "132", "233"]),
    (14, "211", ["111", "221", "311"]),
    (15, "212", ["112", "211", "222"]),
    (16, "213", ["113", "212", "223"]),
    (17, "312", ["212", "311", "322"]),
    (18, "313", ["213", "312", "323"]),
]

# Census numbers for order three.  NOTE ON PROVENANCE: every constant here is
# printed in the paper.  P_NCLASS (seven isometry classes), P_CLASS_AREA
# (their totals) and P_S3_OF_3 appear in the reverse-inequality paragraph of
# the proof; P_NIND and P_NSOL appear in its closing sentence ("of the 4224
# independent 9-subsets of [3]^3, exactly 116 are solutions") and again in the
# Exact-verification paragraph.  Comparing any of them with the census is
# therefore a test of the paper, and a disagreement would refute it.  The two
# census counts are corroborated twice over besides: 4224 is recomputed by the
# transfer-matrix DP (a wholly different algorithm), and 116 is cross-checked
# against the paper's class count, since the 116 solutions must fall into
# exactly P_NCLASS = 7 isometry orbits.
P_NIND = 4224                                  # independent 9-subsets (paper)
P_NSOL = 116                                   # solutions among them (paper)
P_NCLASS = 7                                   # isometry classes (paper)
P_CLASS_AREA = [27, 27, 25, 23, 21, 20, 21]    # their totals (paper)
P_S3_OF_3 = 20                                 # s_3(3) = 20 (paper)

# This program prints exactly this many checks on the true object.  The
# structural guard in finish() enforces it, so deleting or short-circuiting a
# whole section cannot leave a green verdict behind (it used to: removing the
# call to section_refutation printed "ALL 30 CHECKS PASS" and exited 0).
EXPECTED_CHECKS = 32

_RESULTS = []


def check(name, ok, detail=""):
    """Record one check.  ok must be a genuine boolean computation."""
    _RESULTS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    return bool(ok)


def info(text):
    print("INFO " + str(text))


def finish():
    # Structural guard.  Two independent failure modes must be caught here:
    #   (a) a check that silently ran twice under one name;
    #   (b) a section that never ran at all, which leaves its checks simply
    #       ABSENT -- nothing fails, so the verdict would read
    #       "ALL 30 CHECKS PASS" and the exit status would be 0.  Counting the
    #       names is the only thing that sees (b), so the count is asserted
    #       against EXPECTED_CHECKS, not merely reported.
    names = [nm for nm, _ in _RESULTS]
    dups = sorted(set(nm for nm in names if names.count(nm) > 1))
    want = EXPECTED_CHECKS - 1          # every check except this one
    check("check_names_are_distinct_and_every_section_ran_once",
          len(dups) == 0 and len(names) == len(set(names))
          and len(names) == want,
          "%d checks ran, %d expected, %d distinct names%s"
          % (len(names), want, len(set(names)),
             "" if not dups else ", duplicated: " + str(dups)))
    n = len(_RESULTS)
    bad = [nm for nm, ok in _RESULTS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


# =====================================================================
# BLOCK 2.  The grid graph on [n]^3, straight from the definition.
# =====================================================================

def cube_vertices(n):
    """All of [n]^3, in lexicographic order."""
    return [(i, j, k) for i in range(1, n + 1)
            for j in range(1, n + 1) for k in range(1, n + 1)]


def is_grid_adjacent(u, v):
    """Adjacent iff they differ by 1 in exactly one coordinate.

    Implemented literally: count coordinates where they differ, and require
    that count to be 1 and the difference there to have absolute value 1.
    """
    diffs = [a - b for a, b in zip(u, v)]
    nz = [d for d in diffs if d != 0]
    return len(nz) == 1 and (nz[0] == 1 or nz[0] == -1)


class Cube(object):
    """Index/bitmask machinery for [n]^3 plus neighbour masks."""

    def __init__(self, n):
        self.n = n
        self.verts = cube_vertices(n)
        self.N = len(self.verts)
        self.idx = dict((v, t) for t, v in enumerate(self.verts))
        self.nbr = []          # nbr[t] = list of neighbour indices
        self.nbrmask = []      # nbrmask[t] = bitmask of neighbours
        for t, u in enumerate(self.verts):
            lst = [s for s, w in enumerate(self.verts) if is_grid_adjacent(u, w)]
            self.nbr.append(lst)
            m = 0
            for s in lst:
                m |= 1 << s
            self.nbrmask.append(m)
        self.full = (1 << self.N) - 1

    def mask_of(self, vs):
        m = 0
        for v in vs:
            m |= 1 << self.idx[v]
        return m

    def verts_of(self, mask):
        return [v for t, v in enumerate(self.verts) if (mask >> t) & 1]


def popcount(x):
    return bin(x).count("1")


def parse_vertex(s):
    """Decode the paper's compact notation 'ijk' into a tuple."""
    if not isinstance(s, str) or len(s) != 3 or not s.isdigit():
        return None
    return (int(s[0]), int(s[1]), int(s[2]))


def projections(vs):
    """pi_i deletes the i-th coordinate; return the three image sets."""
    p1 = set((v[1], v[2]) for v in vs)
    p2 = set((v[0], v[2]) for v in vs)
    p3 = set((v[0], v[1]) for v in vs)
    return {1: p1, 2: p2, 3: p3}


def area(vs):
    """A(B) = |pi_1(B)| + |pi_2(B)| + |pi_3(B)|, exact integer."""
    p = projections(vs)
    return len(p[1]) + len(p[2]) + len(p[3])


def conjectured_bound(n):
    """The Barat--Korondi--Varga lower bound n^2 + 6n - 4."""
    return n * n + 6 * n - 4


# =====================================================================
# BLOCK 3.  Independence, and the exact "exactly three neighbours" search.
# =====================================================================

def independent_pairs(vs):
    """Return the adjacent pairs inside vs (empty iff vs is independent)."""
    bad = []
    for u, w in combinations(vs, 2):
        if is_grid_adjacent(u, w):
            bad.append((u, w))
    return bad


def is_solution(cube, base, want_order=False):
    """Is [n]^3 reachable from `base` under the EXACTLY-three rule?

    The rule is not monotone (an addition can push a vertex from three
    present neighbours to four, destroying its eligibility for ever), so a
    greedy closure would be unsound.  We therefore do a complete
    backtracking search over addition orders, memoising states from which
    the full cube is provably unreachable.
    """
    target = cube.full
    dead = set()
    order = []
    nbrmask = cube.nbrmask
    Nv = cube.N

    def dfs(state):
        if state == target:
            return True
        if state in dead:
            return False
        for t in range(Nv):
            bit = 1 << t
            if state & bit:
                continue
            if popcount(state & nbrmask[t]) == 3:
                order.append(cube.verts[t])
                if dfs(state | bit):
                    return True
                order.pop()
        dead.add(state)
        return False

    ok = dfs(base)
    if want_order:
        return ok, list(order)
    return ok


def bootstrap_closure(cube, base):
    """Fixpoint of the relaxed AT-LEAST-three rule (monotone, so confluent).

    Used only for a printed robustness note: it is an upper bound on what
    the exactly-three rule can reach.
    """
    state = base
    changed = True
    while changed:
        changed = False
        for t in range(cube.N):
            bit = 1 << t
            if state & bit:
                continue
            if popcount(state & cube.nbrmask[t]) >= 3:
                state |= bit
                changed = True
    return state


def enumerate_independent(cube, size):
    """All independent subsets of [n]^3 of the given size, as bitmasks."""
    out = []
    Nv = cube.N
    nbrmask = cube.nbrmask

    def rec(start, chosen, count, blocked):
        if count == size:
            out.append(chosen)
            return
        for t in range(start, Nv):
            if Nv - t < size - count:
                break
            if (blocked >> t) & 1:
                continue
            rec(t + 1, chosen | (1 << t), count + 1, blocked | nbrmask[t])

    rec(0, 0, 0, 0)
    return out


# =====================================================================
# BLOCK 4.  The cube isometry group, as permutations of [n]^3.
# =====================================================================

def isometries(cube):
    """The signed coordinate permutations of [n]^3, as index permutations."""
    n = cube.n
    out = set()
    for perm in permutations(range(3)):
        for flip in product((False, True), repeat=3):
            f = []
            for v in cube.verts:
                w = []
                for i in range(3):
                    c = v[perm[i]]
                    w.append(n + 1 - c if flip[i] else c)
                f.append(cube.idx[tuple(w)])
            out.add(tuple(f))
    return sorted(out)


def apply_iso(cube, mask, f):
    r = 0
    m = mask
    t = 0
    while m:
        if m & 1:
            r |= 1 << f[t]
        m >>= 1
        t += 1
    return r


def compose(f, g):
    """(f after g)."""
    return tuple(f[g[t]] for t in range(len(g)))


def group_is_closed(maps):
    """Identity present, closed under composition, closed under inverses."""
    S = set(maps)
    ident = tuple(range(len(maps[0])))
    if ident not in S:
        return False
    for f in maps:
        inv = [0] * len(f)
        for t, s in enumerate(f):
            inv[s] = t
        if tuple(inv) not in S:
            return False
        for g in maps:
            if compose(f, g) not in S:
                return False
    return True


def orbits(cube, masks, maps):
    """Partition `masks` into isometry orbits.

    Returns triples (representative, orbit, orbit_stays_inside_masks).  The
    representative is always a member of `masks`, so downstream code stays
    well defined even when `maps` is not in fact a symmetry group.
    """
    pool = set(masks)
    seen = set()
    out = []
    for m in masks:
        if m in seen:
            continue
        orb = set(apply_iso(cube, m, f) for f in maps)
        seen |= orb
        out.append((m, orb, orb <= pool))
    return out


# =====================================================================
# BLOCK 4b.  Unit tests of the machinery the verdict rests on.
# =====================================================================

def section_machinery(cube):
    """Pin down the graph and, crucially, the EXACTLY-three addition rule.

    The 3x3x3 grid must come out with 8 corners of degree 3, 12 edge
    midpoints of degree 4, 6 face centres of degree 5, one centre of degree
    6, and 54 edges.  And the rule must be exactly-three, not at-least-
    three: [3]^3 minus its centre is a state in which the only missing
    vertex has SIX present neighbours, so it is unreachable under the
    paper's rule though reachable under the relaxed one, while [3]^3 minus
    a corner (three present neighbours) is reachable under both.
    """
    degs = {}
    for lst in cube.nbr:
        degs[len(lst)] = degs.get(len(lst), 0) + 1
    edges = sum(len(l) for l in cube.nbr) // 2
    check("grid_graph_on_[3]^3_has_the_right_degrees_and_54_edges",
          degs == {3: 8, 4: 12, 5: 6, 6: 1} and edges == 54,
          "degree histogram %s, %d edges" % (str(sorted(degs.items())), edges))

    centre = cube.full ^ (1 << cube.idx[(2, 2, 2)])
    corner = cube.full ^ (1 << cube.idx[(1, 1, 1)])
    edgemid = cube.full ^ (1 << cube.idx[(1, 1, 2)])
    r1 = is_solution(cube, centre)
    r2 = is_solution(cube, corner)
    r3 = is_solution(cube, edgemid)
    r4 = bootstrap_closure(cube, centre) == cube.full
    check("addition_rule_is_EXACTLY_three_not_at_least_three",
          (r1 is False) and (r2 is True) and (r3 is False) and (r4 is True),
          "cube minus centre (6 nbrs) reachable? %s; minus corner (3 nbrs)? "
          "%s; minus edge midpoint (4 nbrs)? %s; minus centre under the "
          "relaxed >=3 rule? %s" % (r1, r2, r3, r4))


# =====================================================================
# BLOCK 5.  Checks on the exhibited object itself.
# =====================================================================

def section_object(cube):
    """Decode B, print it back, and test every hypothesis of the definition.

    Returns the decoded B (list of tuples) or None if it does not decode.
    """
    n = cube.n
    decoded = [parse_vertex(s) for s in P_B]
    wellformed = (None not in decoded
                  and all(all(1 <= c <= n for c in v) for v in decoded if v)
                  and len(set(decoded)) == len(decoded))
    check("B_decodes_to_distinct_points_of_[3]^3", wellformed,
          "B = " + ", ".join("".join(map(str, v)) for v in decoded
                             if v is not None))
    if not wellformed:
        return None
    B = decoded
    check("B_has_size_n^2_=_9", len(B) == n * n,
          "|B| = %d, n^2 = %d" % (len(B), n * n))

    bad = independent_pairs(B)
    check("B_is_an_independent_set", len(bad) == 0,
          "%d of %d pairs tested are adjacent%s"
          % (len(bad), len(B) * (len(B) - 1) // 2,
             "" if not bad else ": " + str(bad[:3])))
    return B


def section_projections(B):
    """Recompute the projections from B and confront the printed display."""
    got = projections(B)
    same = all(got[i] == set(P_PROJ[i]) for i in (1, 2, 3))
    # Set equality alone does not see a repeated element in the printed
    # display: {(1,1),(1,1),...} and {(1,1),...} are the same set, so a
    # display listing six entries of which two coincide -- which would make
    # the printed |pi_1| wrong -- used to pass.  The table rows already carry
    # this guard (len(set(printedN)) == 3); the display needs it too.  So the
    # printed lists are required to be repetition-free and to have exactly as
    # many entries as the recomputed image has elements.
    nodup = all(len(P_PROJ[i]) == len(set(P_PROJ[i])) for i in (1, 2, 3))
    samelen = all(len(P_PROJ[i]) == len(got[i]) for i in (1, 2, 3))
    detail = "; ".join("pi_%d=%s (display lists %d entries)"
                       % (i, sorted(got[i]), len(P_PROJ[i])) for i in (1, 2, 3))
    check("projections_recomputed_equal_printed_display",
          same and nodup and samelen, detail)

    sizes = (len(got[1]), len(got[2]), len(got[3]))
    check("projection_sizes_are_(6,7,7)", sizes == P_PROJ_SIZES,
          "computed %s, paper %s" % (str(sizes), str(P_PROJ_SIZES)))

    a = area(B)
    check("total_projection_area_A(B)_=_20", a == P_AREA,
          "computed A(B) = %d, paper %d" % (a, P_AREA))
    return a


def refutes(area_value, bound_value):
    """Does an area of `area_value` refute the lower bound `bound_value`?

    The conjecture asserts s_3(n) >= n^2+6n-4, so a solution refutes it iff
    its area is STRICTLY LESS than the bound.  Isolated as a named predicate
    so that its direction and its strictness can be exercised on inputs where
    a wrong criterion would visibly differ (see section_refutation).
    """
    return area_value < bound_value


def section_refutation(a, n):
    """The load-bearing step: A(B) is strictly below the conjectured bound."""
    bnd = conjectured_bound(n)
    check("conjectured_bound_n^2+6n-4_at_n=3_is_23", bnd == P_BOUND,
          "n^2+6n-4 = %d, paper %d" % (bnd, P_BOUND))
    # The refutation itself, plus a directional self-test of the criterion.
    # A criterion of the wrong shape would still say "yes" on this input:
    # "a != bnd" is true for A(B) = 20 and bnd = 23, and so is "a <= bnd".
    # Both are wrong -- the census below contains real solutions of area 25
    # and 27, which must NOT count as refutations, and a solution of area
    # exactly 23 meets the bound.  So the predicate is required to reject
    # bnd + 2 and to reject bnd, as well as to accept A(B).
    verdict = refutes(a, bnd)
    directional = (refutes(bnd + 2, bnd) is False
                   and refutes(bnd, bnd) is False
                   and refutes(bnd - 1, bnd) is True)
    check("REFUTATION_A(B)_<_n^2+6n-4", verdict is True and directional,
          "%d < %d is %s; criterion rejects area %d and area %d as required: %s"
          % (a, bnd, verdict, bnd + 2, bnd, directional))
    return bnd


# =====================================================================
# BLOCK 6.  The printed build-up certificate, row by row.
# =====================================================================

def section_table(cube, B):
    """Replay the 18-row table; every cell is recomputed from adjacency."""
    n = cube.n
    rows = [(k, parse_vertex(v), [parse_vertex(w) for w in N])
            for (k, v, N) in P_TABLE]
    ks = [k for (k, _, _) in rows]
    vs = [v for (_, v, _) in rows]
    ok_shape = (len(rows) == 27 - n * n
                and ks == list(range(1, len(rows) + 1))
                and None not in vs
                and all(v in cube.idx for v in vs if v)
                and len(set(vs)) == len(vs))
    check("table_has_18_rows_with_distinct_vertices_of_[3]^3", ok_shape,
          "%d rows, %d distinct v_k" % (len(rows), len(set(vs))))
    if not ok_shape:
        return False

    complement = set(cube.verts) - set(B)
    check("table_vertices_are_exactly_[3]^3_minus_B",
          set(vs) == complement,
          "|{v_k}| = %d, |[3]^3 \\ B| = %d, symmetric difference %d"
          % (len(set(vs)), len(complement),
             len(set(vs) ^ complement)))

    S = set(B)
    row_fail = []
    for (k, v, printedN) in rows:
        trueN = set(w for w in S if is_grid_adjacent(v, w))
        if (trueN != set(printedN) or len(trueN) != 3
                or len(printedN) != 3 or len(set(printedN)) != 3):
            row_fail.append((k, "".join(map(str, v)),
                             sorted("".join(map(str, w)) for w in trueN)))
        S.add(v)
    check("table_neighbourhoods_recomputed_and_all_of_size_exactly_3",
          len(row_fail) == 0,
          "18 rows replayed, %d disagree%s" % (len(row_fail),
              "" if not row_fail else ": " + str(row_fail[:3])))
    check("table_terminates_with_S_18_=_[3]^3",
          S == set(cube.verts),
          "|S_18| = %d of %d" % (len(S), cube.N))
    return len(row_fail) == 0 and S == set(cube.verts)


# =====================================================================
# BLOCK 7.  B is a solution -- proved again without the paper's table.
# =====================================================================

def section_search(cube, B, a):
    """Independent backtracking search, ignoring P_TABLE entirely."""
    bm = cube.mask_of(B)
    ok, order = is_solution(cube, bm, want_order=True)
    check("independent_search_certifies_B_is_a_solution", ok,
          ("own order: " + " ".join("".join(map(str, v)) for v in order))
          if ok else "exhaustive search found NO legal build-up from B")
    if ok:
        # Replay the search's own order under the exactly-three rule.
        S = bm
        good = len(order) == cube.N - len(B)
        for v in order:
            t = cube.idx[v]
            if popcount(S & cube.nbrmask[t]) != 3 or (S >> t) & 1:
                good = False
                break
            S |= 1 << t
        check("search_order_replays_legally_to_the_full_cube",
              good and S == cube.full,
              "%d additions, |S| = %d" % (len(order), popcount(S)))
    check("upper_bound_s_3(3)_<=_20", ok and a <= P_S3_OF_3,
          "B is a solution: %s, and A(B) = %d <= %d (paper's s_3(3)): %s"
          % (ok, a, P_S3_OF_3, a <= P_S3_OF_3))
    return ok


# =====================================================================
# BLOCK 8.  The complete order-three census.
# =====================================================================

def census(n, verbose_name):
    """Full census for [n]^3: independent n^2-sets, solutions, min area."""
    cube = Cube(n)
    ind = enumerate_independent(cube, n * n)
    sols = [m for m in ind if is_solution(cube, m)]
    areas = dict((m, area(cube.verts_of(m))) for m in sols)
    best = min(areas.values()) if sols else None
    info("census n=%d (%s): %d independent %d-subsets, %d solutions, "
         "min A = %s, bound n^2+6n-4 = %d"
         % (n, verbose_name, len(ind), n * n, len(sols), str(best),
            conjectured_bound(n)))
    return cube, ind, sols, areas, best


def section_census3(cube, B, ind, sols, areas, best):
    """Counts, isometry classes, class areas, and s_3(3) = min A."""
    check("census_counts_4224_independent_9_subsets_of_[3]^3",
          len(ind) == P_NIND,
          "enumerated %d, paper %d" % (len(ind), P_NIND))
    check("census_counts_116_order_three_solutions",
          len(sols) == P_NSOL,
          "found %d, paper %d" % (len(sols), P_NSOL))
    check("B_is_one_of_the_enumerated_solutions",
          cube.mask_of(B) in set(sols),
          "B %s in the %d solutions found by the independent enumeration"
          % ("IS" if cube.mask_of(B) in set(sols) else "is NOT", len(sols)))

    maps = isometries(cube)
    closed = group_is_closed(maps)
    check("cube_isometry_group_built_has_order_48_and_is_a_group",
          len(maps) == 48 and closed,
          "|G| = %d, closed under composition and inverses: %s"
          % (len(maps), closed))

    orbs = orbits(cube, sols, maps)
    nclosed = sum(1 for (_, _, c) in orbs if c)
    check("isometry_group_permutes_the_solution_set",
          nclosed == len(orbs),
          "%d of %d orbits stay inside the solution set"
          % (nclosed, len(orbs)))

    inv_bad = 0
    for (_, orb, _) in orbs:
        vals = set(areas.get(m, area(cube.verts_of(m))) for m in orb)
        if len(vals) != 1:
            inv_bad += 1
    check("total_projection_area_is_isometry_invariant", inv_bad == 0,
          "%d of %d orbits carry more than one value of A"
          % (inv_bad, len(orbs)))

    check("census_gives_7_isometry_classes", len(orbs) == P_NCLASS,
          "computed %d classes (orbit sizes %s summing to %d), paper %d"
          % (len(orbs), sorted(len(o) for (_, o, _) in orbs),
             sum(len(o) for (_, o, _) in orbs), P_NCLASS))

    got = sorted(areas[rep] for (rep, _, _) in orbs)
    check("class_area_multiset_matches_the_paper's_seven_totals",
          got == sorted(P_CLASS_AREA),
          "computed %s, paper (sorted) %s" % (got, sorted(P_CLASS_AREA)))

    check("s_3(3)_=_min_A_over_all_solutions_=_20", best == P_S3_OF_3,
          "exhaustive minimum %s, paper %d" % (str(best), P_S3_OF_3))

    twenty = [orb for (rep, orb, _) in orbs if areas[rep] == P_S3_OF_3]
    check("B_lies_in_the_unique_class_attaining_the_minimum",
          len(twenty) == 1 and cube.mask_of(B) in twenty[0],
          "%d class(es) attain A = %d" % (len(twenty), P_S3_OF_3))
    return maps


# =====================================================================
# BLOCK 9.  Controls and a reading-robustness test.
# =====================================================================

def section_controls():
    """n = 1 and n = 2, where the conjecture should hold with EQUALITY.

    If the census machinery undershot -- a too-permissive addition rule, a
    wrong projection, a leaky enumeration -- these would come out below the
    bound too, and the n=3 violation would be an artefact.
    """
    for n in (1, 2):
        cube, ind, sols, areas, best = census(n, "control")
        bnd = conjectured_bound(n)
        check("control_n=%d_conjecture_holds_with_equality_s_3(%d)_=_%d"
              % (n, n, bnd), sols and best == bnd,
              "%d solutions, min A = %s, bound %d" % (len(sols), str(best), bnd))


def section_robustness(cube, ind, sols, areas):
    """Does the census depend on reading "exactly 3" rather than ">= 3"?

    The relaxed rule is monotone, so its reachable set is the bootstrap
    closure.  Recomputing the census under that reading is an independent
    route to the same numbers.
    """
    relaxed = [m for m in ind if bootstrap_closure(cube, m) == cube.full]
    rbest = min(area(cube.verts_of(m)) for m in relaxed) if relaxed else None
    check("census_is_the_same_under_the_relaxed_>=3_reading",
          set(relaxed) == set(sols) and rbest == P_S3_OF_3,
          "%d relaxed solutions vs %d exact, min A = %s"
          % (len(relaxed), len(sols), str(rbest)))
    below = [m for m in sols if areas[m] < conjectured_bound(cube.n)]
    info("%d of the %d order-three solutions have A < %d; the smallest A "
         "over ALL %d independent 9-subsets (solution or not) is %d"
         % (len(below), len(sols), conjectured_bound(cube.n), len(ind),
            min(area(cube.verts_of(m)) for m in ind)))


def count_independent_layerdp(n, size):
    """Count independent `size`-subsets of [n]^3 by a transfer-matrix DP.

    Wholly different algorithm from enumerate_independent: it slices [n]^3
    into the n layers i = const, precomputes the independent patterns of the
    n x n grid inside a layer, and sweeps layer to layer forbidding equal
    (j,k) in consecutive layers.  Used (i) to confirm the enumeration count
    at n = 3 by a second route, (ii) to size the n = 4 problem honestly.
    """
    cells = [(j, k) for j in range(n) for k in range(n)]
    pos = dict((c, t) for t, c in enumerate(cells))

    def flat_ok(mask):
        for (j, k) in cells:
            if (mask >> pos[(j, k)]) & 1:
                for (dj, dk) in ((0, 1), (1, 0)):
                    nb = (j + dj, k + dk)
                    if nb in pos and (mask >> pos[nb]) & 1:
                        return False
        return True

    states = [m for m in range(1 << len(cells)) if flat_ok(m)]
    pcs = [popcount(s) for s in states]
    dp = {}
    for t, s in enumerate(states):
        if pcs[t] <= size:
            dp[(s, pcs[t])] = dp.get((s, pcs[t]), 0) + 1
    for _ in range(1, n):
        nd = {}
        for (prev, c), v in dp.items():
            for t, s in enumerate(states):
                if s & prev:
                    continue
                c2 = c + pcs[t]
                if c2 > size:
                    continue
                key = (s, c2)
                nd[key] = nd.get(key, 0) + v
        dp = nd
    return sum(v for (s, c), v in dp.items() if c == size)


# =====================================================================
# BLOCK 10.  Driver.
# =====================================================================

def main():
    n = P_N
    cube = Cube(n)
    info("grid graph on [%d]^3: %d vertices, %d edges, degrees %s"
         % (n, cube.N, sum(len(l) for l in cube.nbr) // 2,
            sorted(set(len(l) for l in cube.nbr))))

    section_machinery(cube)
    B = section_object(cube)
    if B is None:
        finish()
    a = section_projections(B)
    section_refutation(a, n)
    section_table(cube, B)
    section_search(cube, B, a)

    cube3, ind, sols, areas, best = census(n, "full order-three census")
    section_census3(cube3, B, ind, sols, areas, best)
    section_robustness(cube3, ind, sols, areas)
    section_controls()

    dpc = count_independent_layerdp(n, n * n)
    check("independent_9_subset_count_confirmed_by_a_second_algorithm",
          dpc == len(ind) == P_NIND,
          "transfer-matrix DP %d, backtracking enumeration %d, paper %d"
          % (dpc, len(ind), P_NIND))

    big = count_independent_layerdp(4, 16)
    info("SCOPE: the order-three census above is COMPLETE -- all %d "
         "independent 9-subsets of [3]^3 were enumerated and each tested "
         "for solution-hood, so the theorem s_3(3) = 20 < 23 is fully "
         "re-derived here, independently of the Barat-Wanless catalogue. "
         "NOT RE-RUN: n >= 4. The conjecture is stated for all n and the "
         "paper claims nothing beyond n = 3; a census at n = 4 would have "
         "to test %d independent 16-subsets of [4]^3 (counted here by the "
         "transfer-matrix DP), about %d times the order-three workload, "
         "which is many orders of magnitude past a 25-minute budget. ALSO "
         "NOT RE-RUN: the mutation suite described in this program's header. "
         "Those 31 corruptions were run while the program was being written, "
         "by editing copies of it; no mutation harness, list or log ships "
         "with it and nothing above re-executes them, so the note's mutation "
         "claim is NOT evidenced by this transcript. ALSO NOT REPRODUCED: "
         "the note's attribution claims, which are about the "
         "literature and not about mathematics -- that the inequality "
         "transcribed here is the one posed by Barat-Korondi-Varga (the note "
         "quotes it at second hand, as the formulation printed as Conjecture "
         "1 of Barat-Wanless, and states that it cannot confirm what number "
         "it carries in the original; this program asserts no number for it), "
         "that the "
         "seven totals appear in Figure 6 of Barat-Wanless in the printed "
         "left-to-right order (only their MULTISET is checked above, which is "
         "all the theorem needs), and that B is isometric to the "
         "representative there of total projection area 20. Those require "
         "reading the cited papers."
         % (len(ind), big, big // len(ind)))
    finish()


if __name__ == "__main__":
    main()

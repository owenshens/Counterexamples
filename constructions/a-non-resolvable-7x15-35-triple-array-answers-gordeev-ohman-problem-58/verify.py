#!/usr/bin/env python3
"""verify.py -- re-derives every computational claim of

    "A non-resolvable (7 x 15, 35)-triple array", answering Problem 58 of
    Gordeev and Ohman, Electron. J. Combin. 33(2) (2026) #P2.35, Section 9.2.

Python 3.9+, STANDARD LIBRARY ONLY: no third-party package, no external data file, no
network.  Exact integer / Fraction arithmetic throughout; no floating point anywhere, so
no comparison in this program depends on a rounding mode.

THE INPUT IS THE OBJECT PRINTED IN THE PAPER.  The 7 x 15 array below is transcribed from
the display in Section 2 of paper.tex and nothing else is read in.  Everything the paper
asserts about it -- the Definition 6 parameters, the row-set multiset, non-resolvability,
the column design's automorphism group, its parallel classes and resolutions, and the
failure of the quad-array condition -- is recomputed here from that array alone.

Each check prints one `PASS <name> [detail]` line.  The program closes with
`VERDICT: ALL <n> CHECKS PASS` and exits 0 if and only if every check passed.
"""

import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations

# ---------------------------------------------------------------------------
# 0.  THE OBJECT, AS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
ARRAY = [
    [1, 10, 18, 20, 19, 27, 15, 12, 5, 17, 22, 4, 35, 33, 34],
    [2, 12, 8, 22, 24, 23, 20, 30, 31, 3, 34, 26, 13, 7, 18],
    [3, 8, 16, 19, 15, 10, 31, 5, 28, 32, 6, 34, 29, 24, 9],
    [4, 1, 14, 3, 25, 28, 30, 16, 13, 33, 23, 21, 18, 10, 32],
    [5, 13, 15, 21, 26, 20, 11, 27, 25, 29, 33, 8, 6, 14, 7],
    [6, 9, 17, 14, 1, 2, 32, 24, 22, 26, 11, 28, 30, 35, 27],
    [7, 11, 2, 9, 23, 29, 4, 21, 17, 12, 16, 35, 19, 31, 25],
]

# The parameters, as the paper states them.  Nothing below reads these as data: each one is
# recomputed from ARRAY and compared.
R_CLAIM, C_CLAIM, V_CLAIM = 7, 15, 35
E_CLAIM, LRC_CLAIM, LRR_CLAIM, LCC_CLAIM = 3, 3, 5, 1
K_CLAIM = 5                       # group size of a resolution: k = c/e
DISTINCT_ROWSETS_CLAIM = 31
MULT_MULTISET_CLAIM = [1] * 27 + [2] * 4
COINCIDENCES_CLAIM = {
    (1, 3, 5): [5, 15],
    (4, 5, 7): [21, 25],
    (3, 4, 6): [28, 32],
    (1, 6, 7): [17, 35],
}
PARALLEL_CLASSES_CLAIM = 3
RESOLUTIONS_CLAIM = 0
AUT_COLUMN_DESIGN_CLAIM = 6
AUT_PG32 = 20160
TRIPLE_INTERSECTION_CLAIM = {0: 45, 1: 226, 2: 43, 3: 1}

# The Fano plane, used only as the model in the lemma's linear algebra and as a positive
# control for the resolvability decider.
FANO = [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 5, 6)]

_PASSED = 0
_FAILED = 0


def check(name, ok, detail=''):
    global _PASSED, _FAILED
    if ok:
        _PASSED += 1
        print('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        _FAILED += 1
        print('FAIL %s%s' % (name, (' ' + detail) if detail else ''))


# ---------------------------------------------------------------------------
# 1.  DERIVED SETS
# ---------------------------------------------------------------------------
NROW = len(ARRAY)
NCOL = len(ARRAY[0]) if ARRAY else 0
ROWSET = [set(row) for row in ARRAY]
COLSET = [set(ARRAY[i][j] for i in range(NROW)) for j in range(NCOL)]
REP = Counter(s for row in ARRAY for s in row)
SYMBOLS = sorted(REP)

# rs(s) = the triple of rows (1-based) in which symbol s occurs.
RS = {s: tuple(sorted(i + 1 for i in range(NROW) if s in ROWSET[i])) for s in SYMBOLS}
RS_CLASSES = {}
for _s, _t in RS.items():
    RS_CLASSES.setdefault(_t, []).append(_s)

# The column design: one block per symbol, namely the set of columns it occupies.
BLOCKS = sorted(tuple(sorted(j + 1 for j in range(NCOL) if s in COLSET[j])) for s in SYMBOLS)
POINTS = list(range(1, NCOL + 1))


# ---------------------------------------------------------------------------
# 2.  THE DECIDER FOR DEFINITION 21 (resolvability), AND ITS CONTROLS
# ---------------------------------------------------------------------------
def resolution_row_side(rs_map, r, k):
    """Definition 21, row side.  -> (resolvable, info).

    A resolvable array partitions its v symbols into r groups of size k = c/e such that all
    symbols of a group occupy the SAME rows.  So a group is a k-subset of one row-set class,
    and a partition exists iff every class size is a multiple of k.  By the paper's lemma
    the r group row-sets must in addition form a 2-(r, 3, 1) design; both conditions are
    tested here, so the decider can answer NO for either reason and is not a machine that
    only ever prints one word.
    """
    classes = {}
    for s, t in rs_map.items():
        classes.setdefault(t, []).append(s)
    sizes = sorted(len(v) for v in classes.values())
    bad = [t for t, v in classes.items() if len(v) % k]
    if bad:
        return False, {'reason': 'a row-set class has size not a multiple of k=%d' % k,
                       'class_sizes': sizes, 'max_class': max(sizes)}
    groups = sum(len(v) // k for v in classes.values())
    if groups != r:
        return False, {'reason': 'the classes yield %d groups, not r=%d' % (groups, r),
                       'class_sizes': sizes}
    triples = []
    for t, v in classes.items():
        triples.extend([t] * (len(v) // k))
    pair_count = Counter()
    for t in triples:
        if len(t) != 3:
            return False, {'reason': 'a group row-set is not a triple', 'class_sizes': sizes}
        for x, y in combinations(sorted(t), 2):
            pair_count[(x, y)] += 1
    want_pairs = list(combinations(range(1, r + 1), 2))
    if sorted(pair_count) != sorted(want_pairs) or set(pair_count.values()) != {1}:
        return False, {'reason': 'the group row-sets do not form a 2-(%d,3,1) design' % r,
                       'class_sizes': sizes}
    return True, {'reason': 'row side of Definition 21 is satisfiable',
                  'class_sizes': sizes}


def parallel_classes(blocks, points):
    """Every set of |points|/3 pairwise disjoint blocks covering all points."""
    bs = [frozenset(b) for b in blocks]
    universe = frozenset(points)
    out = []

    def dfs(chosen, covered):
        if covered == universe:
            out.append(tuple(sorted(tuple(sorted(b)) for b in chosen)))
            return
        p = min(universe - covered)
        for b in bs:
            if p in b and not (b & covered):
                dfs(chosen + [b], covered | b)

    dfs([], frozenset())
    return sorted(set(out))


def steiner_third(blocks):
    """{(x, y): z} for every block {x, y, z} of a Steiner triple system."""
    third = {}
    for b in blocks:
        x, y, z = b
        for a, c, d in ((x, y, z), (x, z, y), (y, z, x)):
            third[(a, c)] = d
            third[(c, a)] = d
    return third


def sts_automorphism_count(blocks, points):
    """|Aut| of a Steiner triple system, by backtracking with full propagation.

    An automorphism is a permutation f of the points with f(third(x,y)) = third(f(x),f(y)),
    so fixing two images forces a third; the search assigns images in point order and
    propagates every forced value before recursing.  Exact and exhaustive.
    """
    third = steiner_third(blocks)
    n = len(points)
    f, used = {}, set()
    total = 0

    def propagate(seed):
        """-> (ok, [newly forced points]).  Extends f by everything the assignment forces."""
        forced = []
        stack = [seed]
        while stack:
            x = stack.pop()
            for y in list(f):
                if y == x:
                    continue
                z = third[(x, y)]
                zi = third[(f[x], f[y])]
                if z in f:
                    if f[z] != zi:
                        return False, forced
                else:
                    if zi in used:
                        return False, forced
                    f[z] = zi
                    used.add(zi)
                    forced.append(z)
                    stack.append(z)
        return True, forced

    def undo(forced):
        for z in forced:
            used.discard(f[z])
            del f[z]

    def dfs(k):
        nonlocal total
        if k > n:
            total += 1
            return
        if k in f:
            dfs(k + 1)
            return
        for img in points:
            if img in used:
                continue
            f[k] = img
            used.add(img)
            ok, forced = propagate(k)
            if ok:
                dfs(k + 1)
            undo(forced)
            used.discard(img)
            del f[k]

    dfs(1)
    return total


def solve_exact(mat, rhs):
    """Exact Gaussian elimination over Q.  -> (rank, solution or None)."""
    n = len(mat)
    a = [[Fraction(x) for x in row] + [Fraction(rhs[i])] for i, row in enumerate(mat)]
    rank = 0
    where = [-1] * n
    for col in range(n):
        piv = next((r for r in range(rank, n) if a[r][col] != 0), None)
        if piv is None:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        inv = Fraction(1) / a[rank][col]
        a[rank] = [x * inv for x in a[rank]]
        for r in range(n):
            if r != rank and a[r][col] != 0:
                fac = a[r][col]
                a[r] = [x - fac * y for x, y in zip(a[r], a[rank])]
        where[col] = rank
        rank += 1
    if rank < n:
        return rank, None
    return rank, [a[where[c]][n] for c in range(n)]


# ---------------------------------------------------------------------------
# 3.  THE CHECKS
# ---------------------------------------------------------------------------
def main():
    print('object: the 7 x 15 array printed in Section 2 of paper.tex, read from this file'
          ' and from nowhere else')
    print('arithmetic: exact integers and Fractions only; no floating point is used')
    print()

    # --- Definition 6: it really is a (7 x 15, 35)-triple array ----------------------
    print('-- Definition 6: the array is a (7 x 15, 35)-triple array --')
    check('shape', NROW == R_CLAIM and all(len(row) == C_CLAIM for row in ARRAY)
          and sum(len(row) for row in ARRAY) == R_CLAIM * C_CLAIM,
          'r=%d rows, c=%d columns, %d cells' % (NROW, NCOL, NROW * NCOL))
    check('symbol-set', SYMBOLS == list(range(1, V_CLAIM + 1)),
          'the entries are exactly the v=%d symbols 1..%d' % (V_CLAIM, V_CLAIM))
    reps = sorted(set(REP.values()))
    check('equireplicate', reps == [E_CLAIM] and E_CLAIM * V_CLAIM == NROW * NCOL,
          'every symbol occurs e=%d times, and e = rc/v = %d/%d = %d'
          % (E_CLAIM, NROW * NCOL, V_CLAIM, E_CLAIM))
    check('row-binary', sorted(set(len(x) for x in ROWSET)) == [C_CLAIM],
          'each of the %d rows holds %d DISTINCT symbols, so no symbol repeats in a row'
          % (NROW, C_CLAIM))
    check('column-binary', sorted(set(len(x) for x in COLSET)) == [R_CLAIM],
          'each of the %d columns holds %d DISTINCT symbols, so no symbol repeats in a column'
          % (NCOL, R_CLAIM))
    lrc = sorted(set(len(ROWSET[i] & COLSET[j]) for i in range(NROW) for j in range(NCOL)))
    check('lambda_rc', lrc == [LRC_CLAIM],
          '|R_i cap C_j| = %d on all %d row-column pairs' % (LRC_CLAIM, NROW * NCOL))
    lrr = sorted(set(len(ROWSET[i] & ROWSET[j]) for i, j in combinations(range(NROW), 2)))
    check('lambda_rr', lrr == [LRR_CLAIM]
          and LRR_CLAIM * (NROW - 1) == NCOL * (E_CLAIM - 1),
          '|R_i cap R_s| = %d on all %d row pairs, and c(e-1)/(r-1) = %d*%d/%d = %d'
          % (LRR_CLAIM, NROW * (NROW - 1) // 2, NCOL, E_CLAIM - 1, NROW - 1, LRR_CLAIM))
    lcc = sorted(set(len(COLSET[i] & COLSET[j]) for i, j in combinations(range(NCOL), 2)))
    check('lambda_cc', lcc == [LCC_CLAIM]
          and LCC_CLAIM * (NCOL - 1) == NROW * (E_CLAIM - 1),
          '|C_j cap C_t| = %d on all %d column pairs, and r(e-1)/(c-1) = %d*%d/%d = %d'
          % (LCC_CLAIM, NCOL * (NCOL - 1) // 2, NROW, E_CLAIM - 1, NCOL - 1, LCC_CLAIM))
    check('non-trivial', max(NROW, NCOL) < V_CLAIM < NROW * NCOL,
          'max(r,c) = %d < v = %d < rc = %d' % (max(NROW, NCOL), V_CLAIM, NROW * NCOL))
    check('non-extremal', V_CLAIM > NROW + NCOL - 1,
          'v = %d > r + c - 1 = %d, so this is OUTSIDE the extremal case v = r+c-1'
          % (V_CLAIM, NROW + NCOL - 1))
    check('cell-admits-resolutions', NCOL % E_CLAIM == 0 and NCOL // E_CLAIM == K_CLAIM
          and K_CLAIM * NROW == V_CLAIM,
          'k = c/e = %d/%d = %d is an integer and rk = %d*%d = v = %d, so the PARAMETER SET'
          ' does admit resolvable arrays and non-resolvability here is substantive, not an'
          ' integrality failure' % (NCOL, E_CLAIM, K_CLAIM, NROW, K_CLAIM, V_CLAIM))

    # --- the row-set multiset -------------------------------------------------------
    print()
    print('-- the row-set multiset, which is the whole proof --')
    check('row-sets-are-triples', all(len(t) == 3 for t in RS.values()),
          'every symbol occupies exactly %d distinct rows' % E_CLAIM)
    mults = sorted(len(v) for v in RS_CLASSES.values())
    check('distinct-row-sets', len(RS_CLASSES) == DISTINCT_ROWSETS_CLAIM,
          '%d distinct row-sets among the %d symbols' % (len(RS_CLASSES), V_CLAIM))
    check('multiplicity-multiset', mults == MULT_MULTISET_CLAIM,
          'multiplicities are [1 x %d, 2 x %d], summing to %d'
          % (mults.count(1), mults.count(2), sum(mults)))
    coinc = dict((t, sorted(v)) for t, v in RS_CLASSES.items() if len(v) > 1)
    check('the-four-coincidences',
          coinc == dict((t, sorted(v)) for t, v in COINCIDENCES_CLAIM.items()),
          'the only repeated row-sets are ' + ', '.join(
              '%s=%s' % ('='.join(str(s) for s in sorted(v)), str(t))
              for t, v in sorted(coinc.items())))
    check('max-multiplicity-below-k', max(mults) == 2 < K_CLAIM,
          'maximum row-set multiplicity is %d, and a resolution needs k = %d symbols with a'
          ' COMMON row-set' % (max(mults), K_CLAIM))

    # An independent brute force of the same fact: no reasoning, just every 5-subset.
    quints = sum(1 for q in combinations(SYMBOLS, K_CLAIM)
                 if len(set(RS[s] for s in q)) == 1)
    check('brute-force-no-common-row-set-quintuple', quints == 0,
          'of all C(%d,%d) = %d five-element sets of symbols, %d have a common row-set'
          % (V_CLAIM, K_CLAIM, len(list(combinations(range(V_CLAIM), K_CLAIM))), quints))

    resolvable, info = resolution_row_side(RS, R_CLAIM, K_CLAIM)
    check('not-resolvable', resolvable is False,
          'the Definition 21 decider answers NO: %s (class sizes %s)'
          % (info['reason'], info['class_sizes']))

    # --- controls on the decider ----------------------------------------------------
    print()
    print('-- controls: the decider is not a machine that only ever says NO --')
    fano_rs = {}
    s = 1
    for t in FANO:
        for _ in range(K_CLAIM):
            fano_rs[s] = t
            s += 1
    ok_pos, info_pos = resolution_row_side(fano_rs, R_CLAIM, K_CLAIM)
    check('control-positive-fano', ok_pos is True and info_pos['class_sizes'] == [5] * 7,
          'on a synthetic row-set map that is five copies of each Fano triple the SAME'
          ' decider answers YES, class sizes %s' % (info_pos['class_sizes'],))
    # A size-legal but non-Fano input must still be rejected, by the lemma's second clause.
    non_fano = [(1, 2, 3), (1, 2, 4), (1, 2, 5), (1, 2, 6), (1, 2, 7), (3, 4, 5), (5, 6, 7)]
    bad_rs = {}
    s = 1
    for t in non_fano:
        for _ in range(K_CLAIM):
            bad_rs[s] = t
            s += 1
    ok_neg, info_neg = resolution_row_side(bad_rs, R_CLAIM, K_CLAIM)
    check('control-negative-non-fano', ok_neg is False
          and '2-(7,3,1)' in info_neg['reason'],
          'a synthetic map with class sizes %s whose 7 row-sets are NOT a 2-(7,3,1) design'
          ' is rejected for that reason: %s' % (info_neg['class_sizes'], info_neg['reason']))

    # --- the lemma's linear algebra -------------------------------------------------
    print()
    print('-- the lemma: five copies of a Fano plane is the only resolvable pattern --')
    N = [[1 if p in t else 0 for t in FANO] for p in range(1, R_CLAIM + 1)]
    NNT = [[sum(N[i][k] * N[j][k] for k in range(R_CLAIM)) for j in range(R_CLAIM)]
           for i in range(R_CLAIM)]
    want = [[2 + 1 if i == j else 1 for j in range(R_CLAIM)] for i in range(R_CLAIM)]
    check('fano-gram-matrix', NNT == want,
          'N N^T = 2I + J for the 7 x 7 Fano incidence matrix N')
    rank, sol = solve_exact(N, [E_CLAIM] * R_CLAIM)
    check('fano-system-unique-solution',
          rank == R_CLAIM and sol == [Fraction(1)] * R_CLAIM,
          'N has rank %d over Q, so N m = %d*1 has the UNIQUE solution m = (1,...,1):'
          ' each group meets each column-set exactly once, and the column clause of'
          ' Definition 21 is automatic' % (rank, E_CLAIM))

    # --- the second, independent proof: the column design ---------------------------
    print()
    print('-- the independent second proof, via the column design --')
    check('column-design-blocks', len(BLOCKS) == V_CLAIM
          and all(len(b) == 3 for b in BLOCKS) and len(set(BLOCKS)) == V_CLAIM,
          '%d distinct blocks of size 3 on %d points' % (len(BLOCKS), NCOL))
    pair_count = Counter()
    for b in BLOCKS:
        for x, y in combinations(b, 2):
            pair_count[(x, y)] += 1
    check('column-design-is-sts15',
          sorted(pair_count) == sorted(combinations(POINTS, 2))
          and set(pair_count.values()) == {1},
          'every one of the %d point pairs lies in exactly one block: a 2-(%d,3,1) design'
          % (NCOL * (NCOL - 1) // 2, NCOL))
    pcs = parallel_classes(BLOCKS, POINTS)
    check('parallel-classes', len(pcs) == PARALLEL_CLASSES_CLAIM,
          'exhaustive search finds exactly %d parallel classes: %s'
          % (len(pcs), ' | '.join(' '.join('{%d,%d,%d}' % b for b in pc) for pc in pcs)))
    check('no-resolutions', len(pcs) < NROW and RESOLUTIONS_CLAIM == 0,
          'a resolution needs %d pairwise disjoint parallel classes and only %d classes'
          ' exist, so the design has %d resolutions and is one of the non-resolvable'
          ' STS(15)' % (NROW, len(pcs), RESOLUTIONS_CLAIM))
    aut = sts_automorphism_count(BLOCKS, POINTS)
    check('column-design-automorphisms', aut == AUT_COLUMN_DESIGN_CLAIM,
          '|Aut| = %d by exhaustive backtracking' % aut)
    check('column-design-is-not-pg32', aut != AUT_PG32,
          '|Aut| = %d differs from |Aut(PG(3,2))| = %d, so the column design here is NOT'
          ' the projective design carried by the previously known (7 x 15, 35)-triple'
          ' array' % (aut, AUT_PG32))

    # --- not a quad array -----------------------------------------------------------
    print()
    print('-- the array is not a quad array, so Question 60 is untouched --')
    tri = Counter()
    for i, j in combinations(range(NROW), 2):
        for k in range(NCOL):
            tri[len(ROWSET[i] & ROWSET[j] & COLSET[k])] += 1
    check('triple-intersection-multiset', dict(tri) == TRIPLE_INTERSECTION_CLAIM
          and sum(tri.values()) == NROW * (NROW - 1) // 2 * NCOL,
          'over all %d (row-pair, column) triples the multiset of |R_i cap R_s cap C_j| is'
          ' %s, which is NOT constant'
          % (sum(tri.values()), dict(sorted(tri.items()))))

    # --- scope ----------------------------------------------------------------------
    print()
    print('NOT RE-RUN: the CP-SAT search that PRODUCED this array is not re-run and is not'
          ' needed -- the array above is the input, and every claim of the paper is'
          ' re-derived from it by the checks above.')
    print('NOT RE-RUN: no minimality, uniqueness or counting claim. This program does not'
          ' enumerate the (7 x 15, 35)-triple arrays, does not decide how many'
          ' non-resolvable ones exist, and does not identify the column design inside any'
          ' published catalogue of the 80 STS(15) isomorphism classes -- it only computes'
          ' that design\'s own invariants.')
    print('NOT RE-RUN: nothing about the PREVIOUSLY known (7 x 15, 35)-triple array. Its'
          ' array is not printed here, so the statements that its column design is PG(3,2)'
          ' and that it is resolvable are quoted from the literature, not verified.')
    print('NOT RE-RUN: no bibliographic claim. Whether this is the first non-resolvable'
          ' triple array on any non-extremal parameter set is a literature question and is'
          ' outside this program.')
    print('NOT RE-RUN: Problems 56 and 57 and Questions 59 and 60 of the same section, and'
          ' the resolvable-array enumeration of the source paper.')

    print()
    total = _PASSED + _FAILED
    if _FAILED:
        print('VERDICT: %d of %d CHECKS FAILED' % (_FAILED, total))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _PASSED)
    return 0


if __name__ == '__main__':
    sys.exit(main())

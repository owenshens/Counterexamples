#!/usr/bin/env python3
"""Class 69 length-two mesh patterns are NOT equidistributed on involutions.

Standard library only, exact integer arithmetic, no data files, no network.

TAKEN FROM THE PAPER (inputs, transcribed below as PAPER_* constants):
  * the four shading sets of the Class 69 patterns p_r = (12, R_r), a mesh cell
    being indexed by its lower-left corner:
        R_1 = {(0,0),(1,1),(1,2),(2,1)}   R_2 = {(0,1),(1,0),(1,1),(2,2)}
        R_3 = {(0,2),(1,0),(1,1),(2,1)}   R_4 = {(0,1),(1,1),(1,2),(2,0)}
  * the occurrence convention: for i<j with s_i<s_j, put x=(0,i,j,n+1) and
    y=(0,s_i,s_j,n+1); the pair is an occurrence iff for every (a,b) in R there
    is no k with x_a<k<x_{a+1} and y_b<s_k<y_{b+1};
  * the exhibited object: the involutions of S_3, listed as 123,132,213,321;
  * the paper's claimed cell of the unused plot point for each increasing pair
    at n=3, its occurrence table, its polynomials 2+q+q^2 and 1+2q+q^2, the
    avoidance counts 2,2,1,1, the totals 3,3,4,4, the two cells (0,2),(2,0)
    said to be unoccupied, the blocked counts 4,4,3,3, and F_2 = 1+q.

DERIVED HERE (computed, never asserted):
  * the involutions of S_3 are recomputed from s^2 = id and compared with the
    exhibited list;
  * the occurrence counter is validated against three independent closed forms
    (empty shading, fully shaded grid, fully shaded left column) and against a
    second, literal transcription of the definition;
  * the cell of every unused plot point, the occurrence table, the four
    distribution polynomials, the avoidance counts and the totals at n=3;
  * the load-bearing inequality F^I_{3,p_1}(q) != F^I_{3,p_3}(q);
  * reverse-complement: that the cell map (a,b)->(2-a,2-b) carries R_1 to R_2
    and R_3 to R_4, and that occurrence counts are equivariant under it;
  * minimality: n=0,1 admit no occurrence and n=2 gives 1+q for all four;
  * that the four shadings form a single orbit closed under every symmetry of
    the square acting on the cell grid, so that no uniform re-reading of the
    grid diagrams (reflection, rotation, or reordering) can change the set;
  * a census over all 512 shadings of the 3x3 grid, the persistence of the
    involution split for n up to 11, and the competing equidistribution over
    the full symmetric groups for n up to 9.
"""

import itertools

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)


# ----------------------------------------------------------------- paper input
PAPER_R = {
    1: frozenset([(0, 0), (1, 1), (1, 2), (2, 1)]),
    2: frozenset([(0, 1), (1, 0), (1, 1), (2, 2)]),
    3: frozenset([(0, 2), (1, 0), (1, 1), (2, 1)]),
    4: frozenset([(0, 1), (1, 1), (1, 2), (2, 0)]),
}
PAPER_I3 = ["123", "132", "213", "321"]
PAPER_CELLS_N3 = {
    "123": {(1, 2): (2, 2), (1, 3): (1, 1), (2, 3): (0, 0)},
    "132": {(1, 2): (2, 1), (1, 3): (1, 2)},
    "213": {(1, 3): (1, 0), (2, 3): (0, 1)},
    "321": {},
}
PAPER_OCC_N3 = {
    1: {"123": 1, "132": 0, "213": 2, "321": 0},
    2: {"123": 1, "132": 2, "213": 0, "321": 0},
    3: {"123": 2, "132": 1, "213": 1, "321": 0},
    4: {"123": 2, "132": 1, "213": 1, "321": 0},
}
PAPER_POLY_N3 = {1: [2, 1, 1], 2: [2, 1, 1], 3: [1, 2, 1], 4: [1, 2, 1]}
PAPER_AVOIDERS_N3 = {1: 2, 2: 2, 3: 1, 4: 1}
PAPER_TOTALS_N3 = {1: 3, 2: 3, 3: 4, 4: 4}
PAPER_UNOCCUPIED_N3 = frozenset([(0, 2), (2, 0)])
PAPER_BLOCKED_N3 = {1: 4, 2: 4, 3: 3, 4: 3}
PAPER_POLY_N2 = [1, 1]
GRID = [(a, b) for a in range(3) for b in range(3)]

# Sweep bounds, named once so that the printed scope disclosure in main() cannot
# drift from the ranges the checks actually cover.
SG_TOP = 9                        # full symmetric groups S_1..S_SG_TOP
SPLIT_SIZES = list(range(3, 12))  # involution sizes for the persistence check


# ---------------------------------------------------------------------- engine
def occupied_cells(sigma, i, j):
    """Cells (a,b) of the 3x3 grid holding at least one plot point of sigma
    other than the two points of the increasing pair (i,j).  Positions and
    values are 1-based; sigma[k-1] = sigma_k."""
    n = len(sigma)
    if not (1 <= i < j <= n) or sigma[i - 1] >= sigma[j - 1]:
        raise ValueError("not an increasing pair")
    xs = (0, i, j, n + 1)
    ys = (0, sigma[i - 1], sigma[j - 1], n + 1)
    cells = set()
    for k in range(1, n + 1):
        if k == i or k == j:
            continue
        v = sigma[k - 1]
        a = [t for t in range(3) if xs[t] < k < xs[t + 1]]
        b = [t for t in range(3) if ys[t] < v < ys[t + 1]]
        if len(a) == 1 and len(b) == 1:
            cells.add((a[0], b[0]))
    return cells


def increasing_pairs(sigma):
    n = len(sigma)
    return [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)
            if sigma[i - 1] < sigma[j - 1]]


def occ(sigma, R):
    """occ_{(12,R)}(sigma): pairs whose shaded cells are all empty."""
    return sum(1 for (i, j) in increasing_pairs(sigma)
               if not (set(R) & occupied_cells(sigma, i, j)))


def occ_literal(sigma, R):
    """Second, literal transcription of the paper's occurrence definition."""
    n = len(sigma)
    total = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            if sigma[i - 1] >= sigma[j - 1]:
                continue
            xs = (0, i, j, n + 1)
            ys = (0, sigma[i - 1], sigma[j - 1], n + 1)
            ok = True
            for (a, b) in R:
                for k in range(1, n + 1):
                    if xs[a] < k < xs[a + 1] and ys[b] < sigma[k - 1] < ys[b + 1]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                total += 1
    return total


def dist_poly(perms, R, engine=occ):
    """Coefficient list of sum_{sigma} q^{occ(sigma)} over the given set."""
    counts = {}
    for s in perms:
        e = engine(s, R)
        counts[e] = counts.get(e, 0) + 1
    top = max(counts) if counts else 0
    return [counts.get(d, 0) for d in range(top + 1)]


def poly_str(coeffs):
    if not coeffs:
        return "0"
    out = []
    for d, c in enumerate(coeffs):
        if c == 0:
            continue
        if d == 0:
            out.append(str(c))
        elif d == 1:
            out.append(("q" if c == 1 else "%dq" % c))
        else:
            out.append(("q^%d" % d if c == 1 else "%dq^%d" % (c, d)))
    return "+".join(out) if out else "0"


def involutions(n):
    """All sigma in S_n with sigma^2 = id, generated directly."""
    out = []

    def rec(free, image):
        if not free:
            out.append(tuple(image[p] for p in range(1, n + 1)))
            return
        m = free[0]
        image[m] = m
        rec(free[1:], image)
        for idx in range(1, len(free)):
            t = free[idx]
            image[m], image[t] = t, m
            rec(free[1:idx] + free[idx + 1:], image)
        image.pop(m, None)

    if n == 0:
        return [()]
    rec(tuple(range(1, n + 1)), {})
    return out


def rc(sigma):
    """Reverse-complement of a permutation."""
    n = len(sigma)
    return tuple(n + 1 - sigma[n - i] for i in range(1, n + 1))


def rc_cells(R):
    """Reverse-complement acting on shaded cells of a length-two pattern."""
    return frozenset((2 - a, 2 - b) for (a, b) in R)


def word(sigma):
    return "".join(str(v) for v in sigma)


def cell_bit(cell):
    return 1 << (3 * cell[0] + cell[1])


def shading_mask(R):
    m = 0
    for c in R:
        m |= cell_bit(c)
    return m


def pair_masks(sigma):
    """One 9-bit occupied-cell mask per increasing pair of sigma."""
    return [shading_mask(occupied_cells(sigma, i, j))
            for (i, j) in increasing_pairs(sigma)]


def occ_from_masks(masks, rmask):
    return sum(1 for m in masks if not (m & rmask))


def poly_from_counts(counts):
    top = max(counts) if counts else 0
    return [counts.get(d, 0) for d in range(top + 1)]


def dist_poly_masks(perms, rmasks):
    """Distribution polynomials for several shadings in one pass."""
    tallies = [{} for _ in rmasks]
    for s in perms:
        masks = pair_masks(s)
        for idx, rm in enumerate(rmasks):
            e = occ_from_masks(masks, rm)
            tallies[idx][e] = tallies[idx].get(e, 0) + 1
    return [poly_from_counts(t) for t in tallies]


def check_shadings_wellformed():
    ok = True
    detail = []
    for r in (1, 2, 3, 4):
        R = PAPER_R[r]
        inside = all(c in GRID for c in R)
        ok = ok and inside and len(R) == 4
        detail.append("R%d=%s" % (r, sorted(R)))
    distinct = len(set(PAPER_R.values())) == 4
    ok = ok and distinct
    print("  shadings: " + " ; ".join(detail))
    ck("shading_sets_wellformed", ok,
       "four distinct 4-cell subsets of the 3x3 grid: %s" % distinct)


def check_engine_empty_shading():
    """occ_{(12,empty)}(sigma) must equal the number of non-inversions."""
    bad = []
    tested = 0
    for n in range(1, 7):
        for s in itertools.permutations(range(1, n + 1)):
            tested += 1
            want = sum(1 for i in range(len(s)) for j in range(i + 1, len(s))
                       if s[i] < s[j])
            if occ(s, frozenset()) != want or occ_literal(s, frozenset()) != want:
                bad.append(word(s))
    ck("engine_empty_shading_counts_noninversions", not bad,
       "n<=6, all %d permutations; mismatches=%d" % (tested, len(bad)))


def check_engine_full_shading():
    """With all nine cells shaded, any point off the pair kills the pair, so
    occ is 0 for n>=3 and equals the non-inversion count for n<=2."""
    full = frozenset(GRID)
    bad = []
    tested = 0
    for n in range(1, 7):
        for s in itertools.permutations(range(1, n + 1)):
            tested += 1
            got = occ(s, full)
            want = 0 if n >= 3 else (1 if s == (1, 2) else 0)
            if got != want or occ_literal(s, full) != want:
                bad.append((word(s), got, want))
    ck("engine_full_shading_kills_all_pairs", not bad,
       "n<=6, all %d permutations; mismatches=%d" % (tested, len(bad)))


def check_engine_left_column():
    """Shading the whole left column forces i=1, so occ = #{j>1: s_j>s_1}."""
    R = frozenset([(0, 0), (0, 1), (0, 2)])
    bad = []
    tested = 0
    for n in range(1, 7):
        for s in itertools.permutations(range(1, n + 1)):
            tested += 1
            want = sum(1 for j in range(1, len(s)) if s[j] > s[0])
            if occ(s, R) != want or occ_literal(s, R) != want:
                bad.append(word(s))
    ck("engine_left_column_closed_form", not bad,
       "n<=6, all %d permutations; mismatches=%d" % (tested, len(bad)))


def check_engine_three_implementations_agree():
    """Set engine, literal transcription of the definition, and bitmask engine
    must agree for every one of the 512 shadings, n<=5."""
    bad = 0
    total = 0
    for n in range(1, 6):
        for s in itertools.permutations(range(1, n + 1)):
            masks = pair_masks(s)
            for k in range(10):
                for S in itertools.combinations(GRID, k):
                    R = frozenset(S)
                    a = occ(s, R)
                    b = occ_literal(s, R)
                    c = occ_from_masks(masks, shading_mask(R))
                    total += 1
                    if not (a == b == c):
                        bad += 1
    ck("engine_three_implementations_agree", bad == 0,
       "%d (permutation, shading) pairs, n<=5; mismatches=%d" % (total, bad))


def check_exhibited_involutions():
    """Recompute I_3 from sigma^2 = id and print it back."""
    brute = [s for s in itertools.permutations(range(1, 4))
             if all(s[s[i] - 1] == i + 1 for i in range(3))]
    gen = involutions(3)
    words = sorted(word(s) for s in gen)
    print("  I_3 recomputed = {%s}" % ", ".join(words))
    ok = (words == sorted(PAPER_I3)
          and sorted(word(s) for s in brute) == words
          and len(gen) == 4)
    ck("exhibited_object_is_the_involutions_of_S3", ok,
       "|I_3|=%d, matches the four words listed" % len(gen))


def check_unused_point_cells():
    """Each increasing pair at n=3 leaves exactly one unused plot point, and it
    lies in the cell the paper names."""
    bad = []
    seen = []
    for s in involutions(3):
        w = word(s)
        got = {}
        for (i, j) in increasing_pairs(s):
            cells = occupied_cells(s, i, j)
            if len(cells) != 1:
                bad.append((w, (i, j), "not a single cell"))
                continue
            c = next(iter(cells))
            got[(i, j)] = c
            seen.append(c)
        if got != PAPER_CELLS_N3.get(w):
            bad.append((w, got, PAPER_CELLS_N3.get(w)))
        print("  %s -> %s" % (w, sorted(got.items())))
    ck("unused_point_cells_match_paper", not bad,
       "%d increasing pairs over I_3; discrepancies=%d" % (len(seen), len(bad)))
    distinct = len(set(seen)) == len(seen) == 7
    missing = frozenset(GRID) - set(seen)
    ck("seven_pairs_occupy_seven_distinct_cells",
       distinct and missing == PAPER_UNOCCUPIED_N3,
       "distinct=%s, unoccupied cells=%s" % (distinct, sorted(missing)))


def check_blocked_counts():
    """Coarser invariant: how many of the seven occupied cells each R_r hits."""
    seven = set()
    for s in involutions(3):
        for (i, j) in increasing_pairs(s):
            seven |= occupied_cells(s, i, j)
    got = {r: len(PAPER_R[r] & seven) for r in (1, 2, 3, 4)}
    meets = {r: sorted(PAPER_R[r] & PAPER_UNOCCUPIED_N3) for r in (1, 2, 3, 4)}
    ok = got == PAPER_BLOCKED_N3
    ok = ok and meets[1] == [] and meets[2] == []
    ok = ok and meets[3] == [(0, 2)] and meets[4] == [(2, 0)]
    print("  blocked pairs per pattern: %s" % [got[r] for r in (1, 2, 3, 4)])
    ck("blocked_pair_counts_are_4_4_3_3", ok,
       "computed %s ; R_3 meets %s, R_4 meets %s"
       % ([got[r] for r in (1, 2, 3, 4)], meets[3], meets[4]))


def check_occurrence_table():
    """The 4x4 occurrence table of the proof, recomputed."""
    bad = []
    entries = 0
    for r in (1, 2, 3, 4):
        row = {}
        for s in involutions(3):
            a = occ(s, PAPER_R[r])
            b = occ_literal(s, PAPER_R[r])
            if a != b:
                bad.append(("engine disagreement", r, word(s)))
            row[word(s)] = a
            entries += 1
        if row != PAPER_OCC_N3[r]:
            bad.append((r, row, PAPER_OCC_N3[r]))
        print("  occ_p%d over %s = %s"
              % (r, PAPER_I3, [row[w] for w in PAPER_I3]))
    ck("occurrence_table_n3_matches_paper", not bad,
       "%d entries recomputed; discrepancies=%d" % (entries, len(bad)))


def polys_n3():
    return {r: dist_poly(involutions(3), PAPER_R[r], occ_literal)
            for r in (1, 2, 3, 4)}


def check_polynomials_n3():
    got = polys_n3()
    for r in (1, 2, 3, 4):
        print("  F^I_{3,p%d}(q) = %s" % (r, poly_str(got[r])))
    ok = all(got[r] == PAPER_POLY_N3[r] for r in (1, 2, 3, 4))
    ck("distribution_polynomials_n3_match_paper", ok,
       "p1,p2 -> %s ; p3,p4 -> %s" % (poly_str(got[1]), poly_str(got[3])))


def check_refutation_n3():
    """LOAD BEARING.  The conjecture asserts one polynomial for all four
    patterns; compute all four and exhibit an unequal pair."""
    got = polys_n3()
    pairs_equal = got[1] == got[2] and got[3] == got[4]
    witnesses = [(a, b) for a in (1, 2, 3, 4) for b in (1, 2, 3, 4)
                 if a < b and got[a] != got[b]]
    classes = sorted(set(tuple(got[r]) for r in (1, 2, 3, 4)))
    refuted = len(classes) > 1
    print("  distinct involution distributions at n=3: %s"
          % [poly_str(list(c)) for c in classes])
    ck("equidistribution_on_involutions_fails_at_n3",
       refuted and pairs_equal and len(witnesses) == 4,
       "%d distribution classes; unequal pattern pairs %s"
       % (len(classes), witnesses))


def check_avoidance_counts_n3():
    got = {}
    for r in (1, 2, 3, 4):
        got[r] = sum(1 for s in involutions(3)
                     if occ_literal(s, PAPER_R[r]) == 0)
    ok = got == PAPER_AVOIDERS_N3 and got[1] != got[3]
    print("  involutions of S_3 avoiding p_r: %s" % [got[r] for r in (1, 2, 3, 4)])
    ck("avoidance_counts_n3_are_2_2_1_1", ok,
       "computed %s, and 2 != 1 separates the pairs"
       % [got[r] for r in (1, 2, 3, 4)])


def check_totals_n3():
    """Total occurrence counts, i.e. F'(1)."""
    got = {}
    for r in (1, 2, 3, 4):
        got[r] = sum(occ_literal(s, PAPER_R[r]) for s in involutions(3))
    derivs = {r: sum(d * c for d, c in enumerate(PAPER_POLY_N3[r]))
              for r in (1, 2, 3, 4)}
    ok = got == PAPER_TOTALS_N3 and derivs == got
    print("  total occurrences (F'(1)): %s" % [got[r] for r in (1, 2, 3, 4)])
    ck("total_occurrence_counts_are_3_3_4_4", ok,
       "computed %s" % [got[r] for r in (1, 2, 3, 4)])


def check_minimality():
    """n=0,1 admit no occurrence; n=2 gives 1+q for every pattern; so n=3 is
    the smallest size at which the four distributions can differ."""
    ok = True
    for n in (0, 1):
        for r in (1, 2, 3, 4):
            if any(occ_literal(s, PAPER_R[r]) != 0 for s in involutions(n)):
                ok = False
    p2 = {r: dist_poly(involutions(2), PAPER_R[r], occ_literal)
          for r in (1, 2, 3, 4)}
    ok = ok and all(p2[r] == PAPER_POLY_N2 for r in (1, 2, 3, 4))
    print("  F^I_{2,p_r}(q) = %s for r=1..4" % poly_str(p2[1]))
    ck("counterexample_is_minimal_n_at_most_2_agree", ok,
       "no occurrences at n=0,1; all four give %s at n=2" % poly_str(p2[1]))


def check_rc_pattern_pairing():
    """The quoted parenthetical says p_1,p_2 and p_3,p_4 are reverse-complement
    images of one another; the cell action is (a,b) -> (2-a,2-b)."""
    ok = (rc_cells(PAPER_R[1]) == PAPER_R[2]
          and rc_cells(PAPER_R[2]) == PAPER_R[1]
          and rc_cells(PAPER_R[3]) == PAPER_R[4]
          and rc_cells(PAPER_R[4]) == PAPER_R[3])
    cross = rc_cells(PAPER_R[1]) != PAPER_R[3] and rc_cells(PAPER_R[1]) != PAPER_R[4]
    ck("reverse_complement_pairs_the_shadings", ok and cross,
       "rc(R_1)=R_2 and rc(R_3)=R_4, and rc(R_1) is neither R_3 nor R_4")


def check_shading_set_is_symmetry_closed():
    """The transcription of the four grid diagrams is the one input taken on
    trust.  Its exposure is bounded: the four shadings form a single orbit
    closed under every symmetry of the square acting on the cell grid, so any
    uniform re-reading of the diagrams (reflection, rotation, or a different
    listing order) reproduces the same set, and hence the same two
    polynomials.  Only the labelling within the set can move."""
    S = frozenset(PAPER_R.values())
    rev = lambda c: (2 - c[0], c[1])          # reflection in the vertical axis
    tra = lambda c: (c[1], c[0])              # reflection in the main diagonal
    # A group element is recorded as the tuple of images of GRID, so that
    # elements can be compared and de-duplicated.  Close {rev, tra} to a
    # FIXPOINT: a fixed number of doubling passes reaches only words of bounded
    # length in the generators (three passes stop at length 3 and so miss the
    # 180-degree rotation rev.tra.rev.tra), whereas iterating until the element
    # set stops growing is guaranteed to realise the whole group, and the size
    # is asserted below.
    idx = {c: k for k, c in enumerate(GRID)}
    group = {tuple(GRID)}
    while True:
        grown = set(group)
        for g in group:
            for h in (rev, tra):
                grown.add(tuple(h(c) for c in g))
        if grown == group:
            break
        group = grown
    apply_g = lambda g, X: frozenset(g[idx[c]] for c in X)
    rot180 = tuple((2 - a, 2 - b) for (a, b) in GRID)
    images = set(frozenset(apply_g(g, R) for R in S) for g in group)
    orbit = set(apply_g(g, PAPER_R[1]) for g in group)
    print("  symmetries realised = %d (180-degree rotation present: %s) ; "
          "|D_4 orbit of R_1| = %d ; symmetry images of the set = %d"
          % (len(group), rot180 in group, len(orbit), len(images)))
    ck("shading_set_closed_under_all_square_symmetries",
       len(group) == 8 and rot180 in group
       and images == {S} and orbit == S and len(S) == 4,
       "all 8 symmetries of the square are applied, the 180-degree rotation "
       "among them; the 4 shadings are one D_4 orbit, so a uniformly re-read "
       "diagram gives the same set")


def check_rc_equivariance():
    """Hence occ_{p_1}(sigma) = occ_{p_2}(rc sigma), and rc preserves I_n: the
    two trivial equidistributions the conjecture relies on."""
    bad = 0
    tested = 0
    inv_bad = 0
    for n in range(1, 7):
        invs = set(involutions(n))
        for s in invs:
            if rc(s) not in invs:
                inv_bad += 1
        for s in itertools.permutations(range(1, n + 1)):
            t = rc(s)
            tested += 1
            if occ_literal(s, PAPER_R[1]) != occ_literal(t, PAPER_R[2]):
                bad += 1
            if occ_literal(s, PAPER_R[3]) != occ_literal(t, PAPER_R[4]):
                bad += 1
    ck("rc_equivariance_and_involution_closure", bad == 0 and inv_bad == 0,
       "%d permutations, n<=6; occ mismatches=%d, rc leaves I_n %d times"
       % (tested, bad, inv_bad))


def check_split_persists():
    """The two involution distributions stay pairwise equal and mutually
    distinct well beyond the minimal counterexample."""
    rmasks = [shading_mask(PAPER_R[r]) for r in (1, 2, 3, 4)]
    bad = []
    sizes = SPLIT_SIZES
    for n in sizes:
        p = dist_poly_masks(involutions(n), rmasks)
        if not (p[0] == p[1] and p[2] == p[3] and p[0] != p[2]):
            bad.append(n)
        if n <= 5:
            print("  n=%d: p1,p2 -> %s ; p3,p4 -> %s"
                  % (n, poly_str(p[0]), poly_str(p[2])))
    ck("involution_split_persists_for_n_3_to_11", not bad,
       "%d sizes checked (n=%d..%d), |I_%d|=%d; sizes where the split failed=%s"
       % (len(sizes), sizes[0], sizes[-1], sizes[-1],
          len(involutions(sizes[-1])), bad))


def check_symmetric_group_equidistribution():
    """Contrast claim: on the full symmetric groups the four patterns ARE
    equidistributed, so the failure above is specific to involutions."""
    rmasks = [shading_mask(PAPER_R[r]) for r in (1, 2, 3, 4)]
    bad = []
    top = SG_TOP
    tested = 0
    fact = 1
    for n in range(1, top + 1):
        fact *= n
        tested += fact
        p = dist_poly_masks(itertools.permutations(range(1, n + 1)), rmasks)
        if not (p[0] == p[1] == p[2] == p[3]):
            bad.append(n)
        if n <= 4:
            print("  S_%d: %s" % (n, poly_str(p[0])))
    ck("equidistribution_holds_on_full_symmetric_groups_n_to_9", not bad,
       "n<=%d, %d permutations in total (|S_%d|=%d); sizes with a mismatch=%s"
       % (top, tested, top, fact, bad))


def check_census_all_shadings():
    """Census over all 2^9 shadings of a length-two mesh pattern: group them by
    their involution distributions for n<=6 and locate the four patterns."""
    shadings = [frozenset(S) for k in range(10)
                for S in itertools.combinations(GRID, k)]
    rmasks = [shading_mask(R) for R in shadings]
    sig = {}
    for n in range(1, 7):
        polys = dist_poly_masks(involutions(n), rmasks)
        for R, p in zip(shadings, polys):
            sig.setdefault(R, []).append(tuple(p))
    keys = {R: tuple(v) for R, v in sig.items()}
    rc_bad = sum(1 for R in shadings if keys[rc_cells(R)] != keys[R])
    same = keys[PAPER_R[1]] == keys[PAPER_R[2]]
    same = same and keys[PAPER_R[3]] == keys[PAPER_R[4]]
    diff = keys[PAPER_R[1]] != keys[PAPER_R[3]]
    nclasses = len(set(keys.values()))
    print("  %d shadings fall into %d involution classes for n<=6"
          % (len(shadings), nclasses))
    ck("census_of_512_shadings_separates_the_two_pairs",
       len(shadings) == 512 and same and diff and rc_bad == 0,
       "{p1,p2} share a class, {p3,p4} share a different one; "
       "rc-mismatched shadings=%d" % rc_bad)


def guarded(fn):
    """A corrupted input must produce a FAIL line, never a traceback."""
    try:
        fn()
    except Exception as exc:
        ck(fn.__name__, False, "raised %s: %s" % (type(exc).__name__, exc))


def main():
    print("Class 69 length-two mesh patterns on involutions: "
          "recomputation of the counterexample at n=3")
    for fn in (check_shadings_wellformed,
               check_exhibited_involutions,
               check_engine_empty_shading,
               check_engine_full_shading,
               check_engine_left_column,
               check_engine_three_implementations_agree,
               check_unused_point_cells,
               check_blocked_counts,
               check_occurrence_table,
               check_polynomials_n3,
               check_refutation_n3,
               check_avoidance_counts_n3,
               check_totals_n3,
               check_minimality,
               check_rc_pattern_pairing,
               check_shading_set_is_symmetry_closed,
               check_rc_equivariance,
               check_split_persists,
               check_symmetric_group_equidistribution,
               check_census_all_shadings):
        guarded(fn)
    print(("NOT RE-RUN: the equidistribution over the full symmetric groups is "
           "a result quoted from the literature for all n; it is confirmed here "
           "only for n <= %d, and the involution split only for n <= %d. That "
           "the four shadings above are the ones the conjecture displays, and "
           "that Class 69 is the last open length-two class, are bibliographic "
           "facts about the source and are not machine-checkable.")
          % (SG_TOP, SPLIT_SIZES[-1]))


if __name__ == "__main__":
    main()
    nbad = sum(1 for _, ok in CHECKS if not ok)
    if nbad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (nbad, len(CHECKS)))
        raise SystemExit(1)
    print("VERDICT: ALL %d CHECKS PASS" % len(CHECKS))

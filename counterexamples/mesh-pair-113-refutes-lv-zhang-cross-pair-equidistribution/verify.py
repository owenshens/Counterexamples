#!/usr/bin/env python3
"""Verification of a counterexample to an unnumbered cross-pair
equidistribution assertion for mesh patterns.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs, transcribed, never used as evidence
for themselves):

  * the shading of the pair labelled 113:
        R113 = {(0,1),(0,2),(1,1),(1,2),(2,2),(3,0),(3,1),(3,2),(3,3)}
    with both components being (123, R113) and (132, R113);
  * the shading of the pair X_9^(1):
        R9   = {(0,0),(0,1),(0,2),(0,3),(1,1),(1,2),(3,1),(3,2)}
    with both components being (123, R9) and (132, R9);
  * the asserted conclusion (the statement being refuted): for every pair
    P in {113,...,116}, every pair Q in {X_i^(1), Y_i^(1) : 9<=i<=18},
    every n and every (k,l),  T^P_{n,k,l} = T^Q_{n,k,l};
  * the closed form quoted for the Q-family (equation (1) of the source
    paper), for (k,l) != (0,0):
        T^Q_{n,k,l} = (n-1)! * C(k+l,k) * sum_{j=2}^{n-1} c(j-1,k+l)/j!
    where c is the unsigned Stirling number of the first kind, and the
    branch value 2(n-1)! at (k,l) = (0,0);
  * the numbers the paper asserts: T^P_{n,2,0} = 0; T^Q_{4,2,0} = 1;
    T^Q_{n,2,0} >= (n-1)!/6 > 0 for n >= 4; the two S_4 distributions
    (0,0)->12, (0,1)->5, (1,0)->5, (1,1)->2  for the 113 pair and
    (0,0)->12, (0,1)->4, (0,2)->1, (1,0)->4, (1,1)->2, (2,0)->1 for
    X_9^(1); T^{113}_{n,0,0} = 2(n-1)! for 2 <= n <= 7; the S_3 profile
    one (1,0), one (0,1), four (0,0); n = 4 minimal.

--------------------------------------------------------------------------
DERIVED HERE (computed from the definition of a mesh-pattern occurrence,
by exhaustive enumeration of S_n and exact integer/rational arithmetic --
no value below is read off the paper):

  * the two shadings are well formed, of the stated sizes, and distinct;
  * the unsigned Stirling triangle, from its own recurrence;
  * the full bivariate occurrence distributions of both exhibited pairs
    on S_n for every n <= 8, by exhaustive enumeration;
  * that every component of the 113 pair occurs at most once in any
    permutation of S_n, n <= 8, hence its support lies in {0,1}^2 and
    T^{113}_{n,2,0} = 0 -- the hypothesis side of the refutation;
  * that X_9^(1) does attain the cell (2,0), with the exhibited witness
    1234 and its two occurrence position triples, and the identity of the
    shaded cell that kills each of the other two increasing triples;
  * that the closed form quoted for the Q-family reproduces the directly
    enumerated distribution of X_9^(1) cell by cell, so the exhibited Q
    really is a member of the asserted family;
  * the inequality T^Q_{n,2,0} >= (n-1)!/6 as exact rationals;
  * that the (0,0) cells agree, so an avoidance-only test misses it;
  * minimality: the distributions coincide for n <= 3 and differ at n = 4;
  * sensitivity: perturbing single cells of the shading changes the
    computed occurrence counts in the direction the paper's argument
    requires, so the checks above are driven by the exhibited data;
  * the reading of a shaded cell (a,b) as (column, row), calibrated on a
    hand-computed example that uses none of the paper's data;
  * the inversion mechanism the paper invokes to pass from the pair 113
    to a second member of the family, and the same empty cell (2,0) for
    that inverted pair.
--------------------------------------------------------------------------
"""
import sys
from fractions import Fraction
from itertools import combinations, permutations

CHECKS = []

# ---------------------------------------------------------------- inputs
# Transcribed from the paper.  k = 3 for both components, so shadings are
# subsets of {0,1,2,3} x {0,1,2,3}; (a,b) means column a, row b.
R113 = frozenset([(0, 1), (0, 2), (1, 1), (1, 2), (2, 2),
                  (3, 0), (3, 1), (3, 2), (3, 3)])
R9 = frozenset([(0, 0), (0, 1), (0, 2), (0, 3),
                (1, 1), (1, 2), (3, 1), (3, 2)])
PAT_INC = (1, 2, 3)   # the "123" component
PAT_132 = (1, 3, 2)   # the "132" component
NMAX = 8              # exhaustive enumeration range (8! = 40320)


def flatten(vals):
    """Order-isomorphic reduction of a sequence of distinct values."""
    order = sorted(vals)
    return tuple(order.index(v) + 1 for v in vals)


def mesh_occurrences(sigma, shading, pi):
    """Position tuples of the occurrences of the mesh pattern (sigma, R).

    Definition used: positions i_1<...<i_k carry values order-isomorphic
    to sigma; with i_0 = 0, i_{k+1} = n+1 and j_0 = 0 < j_1 < ... < j_k
    the occurrence values in increasing order, j_{k+1} = n+1, the box
    C_{a,b} = {x : i_a < x < i_{a+1}, j_b < pi(x) < j_{b+1}} must be
    empty for every (a,b) in R.  Positions here are 1-indexed.
    """
    n = len(pi)
    k = len(sigma)
    out = []
    for idx in combinations(range(n), k):
        vals = [pi[i] for i in idx]
        if flatten(vals) != sigma:
            continue
        ii = [0] + [i + 1 for i in idx] + [n + 1]
        jj = [0] + sorted(vals) + [n + 1]
        good = True
        for (a, b) in shading:
            lo_v, hi_v = jj[b], jj[b + 1]
            for x in range(ii[a] + 1, ii[a + 1]):
                if lo_v < pi[x - 1] < hi_v:
                    good = False
                    break
            if not good:
                break
        if good:
            out.append(tuple(i + 1 for i in idx))
    return out


def blocking_cells(sigma, pi, idx0):
    """Cells (a,b) occupied by points outside a classical occurrence."""
    n = len(pi)
    k = len(idx0)
    ii = [0] + [i + 1 for i in idx0] + [n + 1]
    jj = [0] + sorted(pi[i] for i in idx0) + [n + 1]
    cells = set()
    for x in range(1, n + 1):
        if x - 1 in idx0:
            continue
        a = max(t for t in range(k + 1) if ii[t] < x)
        b = max(t for t in range(k + 1) if jj[t] < pi[x - 1])
        cells.add((a, b))
    return cells


def stirling_first_unsigned(rmax):
    """c(r,s), 0 <= r,s <= rmax, from c(r,s) = c(r-1,s-1)+(r-1)c(r-1,s)."""
    c = [[0] * (rmax + 1) for _ in range(rmax + 1)]
    c[0][0] = 1
    for r in range(1, rmax + 1):
        for s in range(0, rmax + 1):
            left = c[r - 1][s - 1] if s >= 1 else 0
            c[r][s] = left + (r - 1) * c[r - 1][s]
    return c


def factorial(m):
    out = 1
    for t in range(2, m + 1):
        out *= t
    return out


def binom(m, r):
    if r < 0 or r > m:
        return 0
    num, den = 1, 1
    for t in range(r):
        num *= m - t
        den *= t + 1
    return num // den


def quoted_formula(n, k, l, c):
    """Equation (1) of the source paper, as an exact Fraction.

    Returns (n-1)! * C(k+l,k) * sum_{j=2}^{n-1} c(j-1,k+l)/j!  for
    (k,l) != (0,0), and the quoted branch 2*(n-1)! at (k,l) = (0,0).
    """
    if (k, l) == (0, 0):
        return Fraction(2 * factorial(n - 1), 1)
    tot = Fraction(0)
    for j in range(2, n):
        # c(r,s) = 0 whenever s > r, so an index past the tabulated width
        # is exactly zero rather than an omission.
        s = k + l
        cval = c[j - 1][s] if s < len(c[j - 1]) else 0
        tot += Fraction(cval, factorial(j))
    return Fraction(factorial(n - 1) * binom(k + l, k), 1) * tot


def joint_dist(shading, n):
    """T^{(shading)}_{n,k,l} as a dict keyed by the vector (k,l)."""
    dist = {}
    for pi in permutations(range(1, n + 1)):
        key = (len(mesh_occurrences(PAT_INC, shading, pi)),
               len(mesh_occurrences(PAT_132, shading, pi)))
        dist[key] = dist.get(key, 0) + 1
    return dist

def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return bool(ok)

def main():
    print("Exhibited objects, decoded:")
    check_shadings_well_formed()
    c = stirling_first_unsigned(2 * NMAX + 8)
    check_stirling(c)
    check_pair_shapes()
    check_witness_1234()
    dist113, dist9 = {}, {}
    for n in range(1, NMAX + 1):
        dist113[n] = joint_dist(R113, n)
        dist9[n] = joint_dist(R9, n)
    check_s4_distributions(dist113[4], dist9[4])
    check_refutation(dist113, dist9)
    check_support_of_113(dist113)
    check_quoted_formula(dist9, c)
    check_two_zero_closed_form(dist9, c)
    check_minimality(dist113, dist9)
    check_zero_cell_agreement(dist113, dist9)
    check_shading_sensitivity()
    check_convention_calibration()
    check_inversion_extends_to_second_pair()
    print("NOT RE-RUN: the enumeration above covers S_n for every n <= %d, "
          "for the two pairs the paper exhibits in full -- the pair "
          "labelled 113 and the pair X_9^(1) -- and, over S_n for n <= 6, "
          "for the inverse of the 113 pair, which the paper's own argument "
          "identifies with a second member of the four-pair family. The "
          "two other members of that four-pair family, labelled 115 and "
          "116, and the remaining nineteen members of the X^(1)/Y^(1) "
          "family, are named only by label in the source paper, whose "
          "shading tables are not reproduced here, so their distributions "
          "were not recomputed: two of the four P-pairs are not "
          "enumerated here at all and rest entirely on the citation to Lv "
          "and Zhang's Theorems 6.1 and 6.2, and for the nineteen "
          "unenumerated Q-pairs the refutation rests on the closed form "
          "quoted for that family, which is checked here cell by cell "
          "against the directly enumerated distribution of X_9^(1)."
          % NMAX)


def finish():
    n = len(CHECKS)
    k = sum(1 for _, ok, _ in CHECKS if not ok)
    if k == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
        sys.exit(0)
    print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
    sys.exit(1)

def check_shadings_well_formed():
    """Decode the exhibited shadings, count them, print them back."""
    grid = set((a, b) for a in range(4) for b in range(4))
    ok = True
    for name, R, size in (("R113", R113, 9), ("R9", R9, 8)):
        inside = R <= grid
        print("    %s (|R| = %d): %s" % (name, len(R), sorted(R)))
        ok = ok and inside and len(R) == size
    ck("shading_well_formed",
       ok,
       "both shadings are subsets of {0..3}^2; |R113|=%d (9), |R9|=%d (8)"
       % (len(R113), len(R9)))
    ck("shadings_distinct",
       R113 != R9,
       "symmetric difference has %d cells" % len(R113 ^ R9))
    # The paper's verbal description of R113: the whole rightmost column
    # and the whole row b = 2 are shaded.
    col3 = set((3, b) for b in range(4))
    row2 = set((a, 2) for a in range(4))
    ck("R113_matches_verbal_description",
       col3 <= R113 and row2 <= R113 and R113 == col3 | row2 | {(0, 1), (1, 1)},
       "column a=3 shaded, row b=2 shaded, remainder {(0,1),(1,1)}")


def check_stirling(c):
    """The Stirling triangle used by the quoted closed form."""
    # sum_s c(r,s) = r!  and  c(r,1) = (r-1)!  are independent identities
    # of the recurrence; a wrong recurrence breaks them.
    rows_ok = all(sum(c[r]) == factorial(r) for r in range(1, 9))
    col1_ok = all(c[r][1] == factorial(r - 1) for r in range(1, 9))
    ck("stirling_recurrence",
       rows_ok and col1_ok,
       "row sums equal r! and c(r,1) = (r-1)! for r <= 8")
    ck("stirling_values_used_in_proof",
       c[1][2] == 0 and c[2][2] == 1 and c[3][2] == 3 and c[4][2] == 11
       and c[5][2] == 50 and c[6][2] == 274,
       "c(1,2)=0, c(2,2)=1, c(3,2)=3, c(4,2)=11, c(5,2)=50, c(6,2)=274")


def check_pair_shapes():
    """Each pair is ((123,R),(132,R)) for a single shading R."""
    ok = (flatten((1, 2, 3)) == PAT_INC and flatten((1, 3, 2)) == PAT_132
          and PAT_INC != PAT_132 and len(PAT_INC) == len(PAT_132) == 3)
    # Both members of a pair carry the SAME shading, and the two classical
    # patterns are the two the source paper studies.
    ck("pair_shape",
       ok,
       "components (123,R) and (132,R): classical patterns 123 and 132, "
       "one shading per pair, length 3")
    # A length-3 pattern cannot occur in S_1 or S_2, for either pair.
    small = []
    for R in (R113, R9):
        for n in (1, 2):
            for pi in permutations(range(1, n + 1)):
                small.append(len(mesh_occurrences(PAT_INC, R, pi))
                             + len(mesh_occurrences(PAT_132, R, pi)))
    ck("no_occurrences_below_order_three",
       set(small) == {0},
       "all occurrence counts vanish on S_1 and S_2")


def check_witness_1234():
    """The exhibited witness: 1234 has occurrence vector (2,0) for X_9."""
    pi = (1, 2, 3, 4)
    inc = mesh_occurrences(PAT_INC, R9, pi)
    dec = mesh_occurrences(PAT_132, R9, pi)
    ck("witness_1234_vector_two_zero",
       len(inc) == 2 and len(dec) == 0,
       "1234: occ_(123,R9) = %d, occ_(132,R9) = %d, vector (%d,%d)"
       % (len(inc), len(dec), len(inc), len(dec)))
    ck("witness_1234_occurrence_positions",
       sorted(inc) == [(1, 2, 3), (1, 2, 4)],
       "surviving position triples %s (paper: (1,2,3) and (1,2,4))"
       % (sorted(inc),))
    # The paper names the cell holding the unused point in each case.
    named = {(1, 2, 3): (3, 3), (1, 2, 4): (2, 2),
             (1, 3, 4): (1, 1), (2, 3, 4): (0, 0)}
    computed = {}
    for idx in named:
        computed[idx] = sorted(blocking_cells(PAT_INC, pi,
                                             tuple(i - 1 for i in idx)))
    cells_ok = all(computed[idx] == [named[idx]] for idx in named)
    survive_ok = all((named[idx] in R9) == (idx not in inc) for idx in named)
    ck("witness_1234_cell_bookkeeping",
       cells_ok and survive_ok,
       "unused point lies in cell "
       + ", ".join("%s->%s" % (idx, computed[idx][0]) for idx in sorted(named))
       + "; exactly the two shaded ones are excluded")
    # 1234 has no classical 132 occurrence at all, so the shading is not
    # what makes the second coordinate zero.
    classical_132 = [idx for idx in combinations(range(4), 3)
                     if flatten([pi[i] for i in idx]) == PAT_132]
    ck("witness_1234_no_classical_132",
       classical_132 == [],
       "1234 contains no classical 132 pattern")


def check_s4_distributions(d113, d9):
    """The two bivariate distributions on S_4, exactly as the paper
    displays them."""
    paper113 = {(0, 0): 12, (0, 1): 5, (1, 0): 5, (1, 1): 2}
    paper9 = {(0, 0): 12, (0, 1): 4, (0, 2): 1,
              (1, 0): 4, (1, 1): 2, (2, 0): 1}
    print("    computed T^113_{4,.,.} = %s" % (sorted(d113.items()),))
    print("    computed T^X9_{4,.,.}  = %s" % (sorted(d9.items()),))
    ck("s4_distribution_pair_113", d113 == paper113,
       "matches the printed distribution, total %d" % sum(d113.values()))
    ck("s4_distribution_X9", d9 == paper9,
       "matches the printed distribution, total %d" % sum(d9.values()))
    ck("s4_distributions_total_24",
       sum(d113.values()) == 24 and sum(d9.values()) == 24,
       "both distributions sum to 4! = 24")
    ck("s4_cells_01_and_10_also_differ",
       d113.get((0, 1)) == 5 and d9.get((0, 1)) == 4
       and d113.get((1, 0)) == 5 and d9.get((1, 0)) == 4,
       "cell (0,1): 5 vs 4; cell (1,0): 5 vs 4")


def check_refutation(dist113, dist9):
    """THE load-bearing check: the asserted equality fails, computed."""
    bad = []
    for n in sorted(dist9):
        if n < 4:
            continue
        p = dist113[n].get((2, 0), 0)
        q = dist9[n].get((2, 0), 0)
        if not (p == 0 and q > 0 and p != q):
            bad.append(n)
    detail = ", ".join("n=%d: T^113=%d vs T^X9=%d"
                       % (n, dist113[n].get((2, 0), 0), dist9[n].get((2, 0), 0))
                       for n in sorted(dist9) if n >= 4)
    ck("conclusion_violated_at_cell_2_0", bad == [], detail)
    ck("distributions_differ_somewhere",
       all(dist113[n] != dist9[n] for n in sorted(dist9) if n >= 4),
       "the two bivariate distributions are unequal for every n >= 4 tested")


def check_support_of_113(dist113):
    """The hypothesis the refutation rests on: each component of the
    113 pair occurs at most once, so its support lies in {0,1}^2."""
    worst = (0, 0)
    for n in sorted(dist113):
        for (k, l) in dist113[n]:
            worst = (max(worst[0], k), max(worst[1], l))
    ck("pair_113_support_in_zero_one_square",
       worst == (1, 1),
       "max (occ_123, occ_132) over S_n, n <= %d, is %s"
       % (max(dist113), worst))
    ck("pair_113_cell_2_0_empty",
       all(dist113[n].get((2, 0), 0) == 0 for n in dist113),
       "T^113_{n,2,0} = 0 for every n <= %d" % max(dist113))


def check_quoted_formula(dist9, c):
    """The quoted closed form reproduces X_9's distribution cell by
    cell, so the exhibited Q really satisfies the asserted equation."""
    mismatch = []
    tested = 0
    for n in sorted(dist9):
        kmax = max([k for (k, l) in dist9[n]] + [0])
        lmax = max([l for (k, l) in dist9[n]] + [0])
        for k in range(kmax + 2):
            for l in range(lmax + 2):
                if (k, l) == (0, 0) and n < 2:
                    continue
                want = quoted_formula(n, k, l, c)
                got = Fraction(dist9[n].get((k, l), 0), 1)
                tested += 1
                if want != got:
                    mismatch.append((n, k, l, str(want), got))
    ck("quoted_formula_reproduces_X9",
       mismatch == [],
       "%d (n,k,l) cells agree with (n-1)!C(k+l,k)sum c(j-1,k+l)/j!; "
       "first mismatch %s" % (tested, mismatch[:1] if mismatch else "none"))


def check_two_zero_closed_form(dist9, c):
    """T^Q_{n,2,0} = (n-1)! sum_{j=2}^{n-1} c(j-1,2)/j! >= (n-1)!/6."""
    rows, bad = [], []
    for n in sorted(dist9):
        if n < 4:
            continue
        val = Fraction(0)
        for j in range(2, n):
            val += Fraction(c[j - 1][2], factorial(j))
        val *= factorial(n - 1)
        got = dist9[n].get((2, 0), 0)
        bound = Fraction(factorial(n - 1), 6)
        rows.append("n=%d: %s" % (n, got))
        if not (val == got and val.denominator == 1 and got >= bound
                and got > 0):
            bad.append(n)
    ck("cell_2_0_closed_form_and_bound",
       bad == [],
       "counts " + "; ".join(rows) + "; each an integer, each >= (n-1)!/6 > 0")
    ck("cell_2_0_equals_one_at_order_four",
       dist9[4].get((2, 0), 0) == 1 and c[1][2] == 0,
       "T^X9_{4,2,0} = 1, the only surviving term being j = 3 since c(1,2)=0")


def check_minimality(dist113, dist9):
    """n = 4 is the first order at which the assertion fails."""
    agree = [n for n in (1, 2, 3) if dist113[n] == dist9[n]]
    ck("orders_below_four_agree",
       agree == [1, 2, 3] and dist113[4] != dist9[4],
       "distributions coincide for n = 1,2,3 and differ at n = 4")
    prof = {(0, 0): 4, (0, 1): 1, (1, 0): 1}
    ck("order_three_profile",
       dist113[3] == prof and dist9[3] == prof,
       "on S_3 both pairs give one (1,0), one (0,1), four (0,0)")


def check_zero_cell_agreement(dist113, dist9):
    """The (0,0) cells agree: an avoidance-only test misses the failure."""
    rows, bad = [], []
    for n in sorted(dist113):
        if n < 2:
            continue
        a = dist113[n].get((0, 0), 0)
        b = dist9[n].get((0, 0), 0)
        rows.append("n=%d: %d" % (n, a))
        if not (a == 2 * factorial(n - 1) and a == b):
            bad.append(n)
    ck("zero_cell_is_twice_factorial_and_agrees",
       bad == [],
       "T^113_{n,0,0} = T^X9_{n,0,0} = 2(n-1)!: " + "; ".join(rows))


def check_shading_sensitivity():
    """The checks above are driven by the exhibited shadings: perturbing a
    single cell moves the computed counts as the paper's argument says."""
    pi = (1, 2, 3, 4)
    base = len(mesh_occurrences(PAT_INC, R9, pi))
    shade22 = len(mesh_occurrences(PAT_INC, R9 | {(2, 2)}, pi))
    shade33 = len(mesh_occurrences(PAT_INC, R9 | {(3, 3)}, pi))
    unshade00 = len(mesh_occurrences(PAT_INC, R9 - {(0, 0)}, pi))
    ck("shading_sensitivity_on_witness",
       base == 2 and shade22 == 1 and shade33 == 1 and unshade00 == 3,
       "occ_123(1234): R9 = %d, +(2,2) = %d, +(3,3) = %d, -(0,0) = %d"
       % (base, shade22, shade33, unshade00))
    # Unshading the rightmost column of R113 lets a second occurrence in,
    # so the "at most one occurrence" hypothesis is a property of R113.
    loose = R113 - {(3, 0), (3, 1), (3, 2), (3, 3)}
    hits = max(len(mesh_occurrences(PAT_INC, loose, p))
               for p in permutations(range(1, 6)))
    ck("shading_sensitivity_on_pair_113",
       hits >= 2,
       "unshading column a=3 of R113 admits %d occurrences on some "
       "permutation of S_5, so the bound is not vacuous" % hits)


def check_convention_calibration():
    """Pin the reading of a shaded cell (a,b) = (column, row) against a
    hand-computed example independent of the paper's data.

    Pattern 12 in the permutation 213.  The two classical occurrences sit
    at positions (1,3) and (2,3).  For (1,3) the leftover point is the
    value 1 at position 2, which lies between the two occurrence
    positions and below both occurrence values, i.e. in the cell
    (column 1, row 0); for (2,3) the leftover point is the value 2 at
    position 1, left of both positions and between the two values, i.e.
    in the cell (column 0, row 1).  So shading {(1,0)} must kill exactly
    the first and shading {(0,1)} exactly the second.
    """
    pi = (2, 1, 3)
    a = mesh_occurrences((1, 2), frozenset([(1, 0)]), pi)
    b = mesh_occurrences((1, 2), frozenset([(0, 1)]), pi)
    ck("cell_index_convention",
       a == [(2, 3)] and b == [(1, 3)],
       "in 213: shading {(1,0)} leaves %s, shading {(0,1)} leaves %s"
       % (a, b))


def inverse_perm(pi):
    out = [0] * len(pi)
    for i, v in enumerate(pi):
        out[v - 1] = i + 1
    return tuple(out)


def transpose(shading):
    return frozenset((b, a) for (a, b) in shading)


def check_inversion_extends_to_second_pair():
    """The paper obtains the pairs 114 and 116 from 113 and 115 by
    inversion.  Verify the mechanism, then verify the conclusion for the
    inverted pair, which is a second member of the four-pair family."""
    bad = []
    for pi in permutations(range(1, 7)):
        for sigma in (PAT_INC, PAT_132):
            if len(mesh_occurrences(sigma, transpose(R113), pi)) != \
               len(mesh_occurrences(sigma, R113, inverse_perm(pi))):
                bad.append(pi)
    ck("inversion_transposes_the_shading",
       bad == [],
       "occ_{(sigma,R^T)}(pi) = occ_{(sigma,R)}(pi^{-1}) on all 720 "
       "permutations of S_6, for both components of the 113 pair")
    dists = {n: joint_dist(transpose(R113), n) for n in range(1, 7)}
    sup = max(max(k, l) for n in dists for (k, l) in dists[n])
    ck("inverted_pair_also_avoids_cell_2_0",
       sup <= 1 and all(dists[n].get((2, 0), 0) == 0 for n in dists),
       "the inverse of the 113 pair has support in {0,1}^2 for n <= 6, so "
       "its cell (2,0) is empty as well")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:              # record the failure as a check
        ck("program_ran_to_completion", False, "%r" % (exc,))
    finish()

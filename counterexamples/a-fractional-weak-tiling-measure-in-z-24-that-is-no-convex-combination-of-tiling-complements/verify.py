#!/usr/bin/env python3
"""Verification program for the note

    A fractional weak tiling measure in Z_24 that is no convex combination of tiling complements

It reads the objects PRINTED IN THE PAPER -- the modulus M = 24, the set E = {0,4,6,10}, and the
measure g -- and re-derives every quantity the paper claims about them.  Python 3.9+, standard
library only (fractions, itertools), no external data file.  All arithmetic is exact integer or
Fraction arithmetic; no floating point value is ever formed, so no decision here depends on a
rounding mode.

One `PASS <name> [detail]` line per check; a closing `VERDICT: ALL <n> CHECKS PASS`; exit 0 if and
only if every check passed.
"""

import sys
from fractions import Fraction
from itertools import combinations

FAILURES = []
NPASS = 0


def check(name, ok, detail=''):
    global NPASS
    if ok:
        NPASS += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        FAILURES.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE %s' % text)


# ---------------------------------------------------------------------------
# exact linear algebra over Q -- Gaussian elimination on Fractions
# ---------------------------------------------------------------------------
def rank(rows):
    """Rank of a matrix given as a list of row lists (ints or Fractions)."""
    A = [[Fraction(v) for v in r] for r in rows]
    nr, nc = len(A), (len(A[0]) if A else 0)
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [v / pv for v in A[r]]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == nr:
            break
    return r


def det(rows):
    """Determinant of a square matrix, exactly."""
    A = [[Fraction(v) for v in r] for r in rows]
    n = len(A)
    d = Fraction(1)
    for c in range(n):
        piv = None
        for i in range(c, n):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            return Fraction(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d *= A[c][c]
        inv = A[c][c]
        A[c] = [v / inv for v in A[c]]
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return d


def solve(rows, rhs):
    """-> (rank, solution or None, unique)   for A x = b over Q.

    `solution` is one solution when the system is consistent, None when it is inconsistent;
    `unique` is True iff rank equals the number of columns.
    """
    A = [[Fraction(v) for v in r] + [Fraction(b)] for r, b in zip(rows, rhs)]
    nr = len(A)
    nc = len(rows[0]) if rows else 0
    where = []
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [v / pv for v in A[r]]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        where.append(c)
        r += 1
        if r == nr:
            break
    for i in range(r, nr):
        if all(A[i][c] == 0 for c in range(nc)) and A[i][nc] != 0:
            return r, None, r == nc
    x = [Fraction(0)] * nc
    for i, c in enumerate(where):
        x[c] = A[i][nc]
    return r, x, r == nc


# ---------------------------------------------------------------------------
# the ambient definitions of the source paper
# ---------------------------------------------------------------------------
def convolution(M, E, g):
    """(1_E * g)(y) = sum_{e in E} g(y - e), as a list of Fractions indexed by y."""
    return [sum((g[(y - e) % M] for e in E), Fraction(0)) for y in range(M)]


def circulant(M, E):
    """A[y][t] = 1 iff y - t lies in E."""
    Es = set(E)
    return [[1 if (y - t) % M in Es else 0 for t in range(M)] for y in range(M)]


def cols(A, idx):
    return [[row[t] for t in idx] for row in A]


def difference_set(M, E):
    return set((a - b) % M for a in E for b in E)


def tiles(M, E, T):
    """True iff E (+) T = Z_M, i.e. the |E||T| sums e+t are distinct and exhaust Z_M."""
    seen = 0
    for t in T:
        m = 0
        for e in E:
            m |= 1 << ((e + t) % M)
        if seen & m:
            return False
        seen |= m
    return seen == (1 << M) - 1


def all_tiling_complements(M, E):
    """Every T subseteq Z_M with |T| = M/|E| and E (+) T = Z_M, by exhaustive enumeration."""
    k = M // len(E)
    masks = []
    for t in range(M):
        m = 0
        for e in E:
            m |= 1 << ((e + t) % M)
        masks.append(m)
    full = (1 << M) - 1
    out = []
    for T in combinations(range(M), k):
        seen = 0
        ok = True
        for t in T:
            if seen & masks[t]:
                ok = False
                break
            seen |= masks[t]
        if ok and seen == full:
            out.append(T)
    return out


def pinned_reduced_system(M, E):
    """The pinned polytope P_E = {g >= 0, 1_E * g = 1, g(0) = 1}, reduced by the pin lemma.

    Setting y = e in E in 1_E * g = 1 and using g(0) = 1 with g >= 0 forces g(e - e') = 0 for all
    distinct e, e' in E.  So supp(g) is contained in {0} union (Z_M minus (E - E)); the free
    variables are the residues of C = Z_M minus (E - E), and the equations become A x = b with
    b(y) = 1 - [y in E].  Returns (C, A, b).
    """
    C = sorted(set(range(M)) - difference_set(M, E))
    Es = set(E)
    A = [[1 if (y - t) % M in Es else 0 for t in C] for y in range(M)]
    b = [1 - (1 if y in Es else 0) for y in range(M)]
    return C, A, b


# ===========================================================================
print('verification of the note: a fractional weak tiling measure in Z_24 that is no convex')
print('combination of tiling complements -- M = 24, E = {0,4,6,10}, item 3 of Kadir-Fan')
print('python %s, exact integer and Fraction arithmetic only' % sys.version.split()[0])
print('')

# ---------------------------------------------------------------------------
print('=== Part A: the exhibited object, exactly as printed in the paper')
M = 24
E = [0, 4, 6, 10]
ONES = [0, 8, 16]
HALVES = [7, 9, 11, 19, 21, 23]
g = [Fraction(0)] * M
for t in ONES:
    g[t] = Fraction(1)
for t in HALVES:
    g[t] = Fraction(1, 2)
note('E = %s,  g = 1 on %s, 1/2 on %s, 0 elsewhere' % (E, ONES, HALVES))

check('A1_g_nonnegative', all(v >= 0 for v in g), 'min = %s' % min(g))
check('A2_g_is_pinned_at_zero', g[0] == 1, 'g(0) = %s (the source Definition requires f(0) = 1)' % g[0])
check('A3_g_is_not_zero_one', any(v not in (0, 1) for v in g),
      'nine support values: %s' % ' '.join(str(g[t]) for t in sorted(ONES + HALVES)))

conv = convolution(M, E, g)
bad = [(y, str(conv[y])) for y in range(M) if conv[y] != 1]
check('A4_convolution_is_one_at_all_24_residues', not bad, 'violations: %s' % bad)

mass = sum(g, Fraction(0))
check('A5_mass_equals_M_over_absE', mass == Fraction(M, len(E)) == 6,
      'sum g = %s = 24/4' % mass)

DD = difference_set(M, E)
check('A6_difference_set_is_all_even', sorted(DD) == [0, 2, 4, 6, 10, 14, 18, 20, 22]
      and all(d % 2 == 0 for d in DD), 'E - E = %s' % sorted(DD))
allowed = set([0]) | (set(range(M)) - DD)
check('A7_support_obeys_the_pin_lemma',
      set(t for t in range(M) if g[t] != 0) <= allowed,
      'supp(g) = %s is inside {0} u (Z_24 \\ (E-E)) = %s'
      % (sorted(t for t in range(M) if g[t] != 0), sorted(allowed)))

# the parity split the paper's proof of A4 uses, checked as a statement in its own right
even_translates = [set((e + s) % M for e in E) for s in ONES]
check('A8_three_even_translates_partition_the_twelve_even_residues',
      len(set().union(*even_translates)) == 12
      and all(x % 2 == 0 for x in set().union(*even_translates))
      and sum(len(s) for s in even_translates) == 12,
      'E+0, E+8, E+16 are disjoint and cover {0,2,...,22}')
odd_cover = {}
for s in HALVES:
    for e in E:
        odd_cover[(e + s) % M] = odd_cover.get((e + s) % M, 0) + 1
asym = [t for t in range(M) if g[t] != g[(-t) % M]]
check('A10_g_is_not_positive_definite', asym != [],
      'g is real but not symmetric -- g(%d) = %s while g(%d) = %s -- and a real positive-definite '
      'function on an abelian group satisfies f(-t) = f(t); so the positive-definite variants asked '
      'in the same list are untouched by this object'
      % (asym[0], g[asym[0]], (-asym[0]) % M, g[(-asym[0]) % M]))
check('A9_six_odd_translates_cover_each_odd_residue_exactly_twice',
      sorted(odd_cover) == list(range(1, M, 2)) and set(odd_cover.values()) == {2},
      'multiplicity 2 on all twelve odd residues, and 2 * (1/2) = 1')
print('')

# ---------------------------------------------------------------------------
print('=== Part B: E is a genuine tile, and ALL of its tiling complements')
comps = all_tiling_complements(M, E)
check('B1_E_has_exactly_sixteen_tiling_complements', len(comps) == 16,
      'exhaustive over all C(24,6) = 134596 six-subsets of Z_24: %d found' % len(comps))
sub8 = [0, 8, 16]
structural = set()
for a in range(4):
    for b in range(4):
        structural.add(tuple(sorted(set((2 * a + h) % M for h in sub8)
                                    | set((2 * b + 1 + h) % M for h in sub8))))
check('B2_the_sixteen_are_exactly_the_structural_family',
      set(comps) == structural,
      'T = (2a + <8>) u (2b+1 + <8>), a,b in {0,1,2,3}')
through0 = sorted(T for T in comps if 0 in T)
check('B3_exactly_four_complements_contain_zero',
      len(through0) == 4 and [list(T) for T in through0] == [[0, 1, 8, 9, 16, 17],
                                                             [0, 3, 8, 11, 16, 19],
                                                             [0, 5, 8, 13, 16, 21],
                                                             [0, 7, 8, 15, 16, 23]],
      '%s' % [list(T) for T in through0])
check('B4_only_complements_through_zero_can_carry_positive_weight',
      len([T for T in comps if 0 not in T]) == 12 and g[0] == 1,
      'the other twelve complements have 1_T(0) = 0 while g(0) = 1, so in any convex combination '
      'summing to g the coefficients on them must all vanish')
in_polytope = []
for T in through0:
    iT = [Fraction(1) if t in T else Fraction(0) for t in range(M)]
    in_polytope.append(all(v == 1 for v in convolution(M, E, iT)) and iT[0] == 1)
check('B5_the_four_complements_through_zero_are_points_of_the_pinned_polytope',
      len(in_polytope) == 4 and all(in_polytope),
      'each 1_T is nonnegative, satisfies 1_E * 1_T = 1 and has 1_T(0) = 1, so the hull of item 3 '
      'is a nonempty subset of P_E -- which is what makes the refutation non-trivial')
print('')

# ---------------------------------------------------------------------------
print('=== Part C: the separating functional -- refutation certificate 1, no linear algebra')
L = lambda x: sum((x[h] for h in HALVES), Fraction(0))
Lg = L(g)
check('C1_L_of_g_equals_three', Lg == 3, 'L(x) = x(7)+x(9)+x(11)+x(19)+x(21)+x(23), L(g) = %s' % Lg)
LT = sorted(set(L([Fraction(1) if t in T else Fraction(0) for t in range(M)]) for T in comps))
check('C2_L_is_at_most_two_on_every_tiling_complement', max(LT) == 2,
      'multiset of values of L(1_T) over all sixteen complements = %s, max = %s'
      % ([str(v) for v in LT], max(LT)))
check('C3_separation_two_is_less_than_three', max(LT) < Lg,
      'L is linear, so L is at most 2 on the whole convex hull while L(g) = 3: g is NOT a convex '
      'combination of tiling complements')
print('')

# ---------------------------------------------------------------------------
print('=== Part D: the vertex certificate -- refutation certificate 2, independent of Part C')
A = circulant(M, E)
supp_full = [0, 8, 16] + HALVES
r_full = rank(cols(A, supp_full))
check('D1_full_circulant_rank_on_the_nine_support_columns',
      r_full == 9 == len(supp_full),
      'rows: all 24 of Z_24; columns %s; rank %d = full column rank' % (supp_full, r_full))
reduced_rows = [y for y in range(M) if y not in set(E)]
supp_red = [8, 16] + HALVES
r_red = rank([[A[y][t] for t in supp_red] for y in reduced_rows])
check('D2_reduced_pinned_system_rank_on_its_eight_free_columns',
      r_red == 8 == len(supp_red),
      'rows: the 20 residues of Z_24 \\ E; columns %s (column 0 is the PINNED coordinate, not a '
      'variable, and vanishes identically on these rows); rank %d' % (supp_red, r_red))
check('D3_the_two_rank_statements_are_about_different_systems',
      r_full != r_red and len(supp_full) != len(supp_red),
      'rank 9 on 24 rows / 9 columns and rank 8 on 20 rows / 8 columns are both correct; '
      'conflating them is not')
odd_rows = [13, 15, 17, 1, 3, 5]
block = [[A[y][t] for t in HALVES] for y in odd_rows]
dblock = det(block)
check('D4_odd_block_determinant_is_four', dblock == 4 or dblock == -4,
      'rows %s x columns %s, block diagonal with two copies of [[1,1,0],[0,1,1],[1,0,1]] of '
      'determinant 2; det = %s' % (odd_rows, HALVES, dblock))
check('D5_g_is_a_vertex_of_the_pinned_polytope', r_full == len(supp_full) and g[0] == 1,
      'the support columns of the equality system (circulant rows plus the pin row) are linearly '
      'independent, so g is a basic feasible solution, hence an extreme point of P_E')
check('D6_an_extreme_point_that_is_not_zero_one_is_no_convex_combination',
      r_full == len(supp_full) and any(v not in (0, 1) for v in g)
      and all(g != [Fraction(1) if t in T else Fraction(0) for t in range(M)] for T in comps),
      'g differs from all sixteen 0/1 points 1_T, and an extreme point of P_E is no convex '
      'combination of points of P_E other than itself')
print('')

# ---------------------------------------------------------------------------
print('=== Part E: anti-controls -- the four literal-reading pseudo-refuters are NOT refuters')
for label, (m, e) in [('E1', (3, [0, 1])), ('E2', (7, [0, 1, 3])),
                      ('E3', (4, [0, 1, 2])), ('E4', (6, [0, 1, 3]))]:
    C, Am, bm = pinned_reduced_system(m, e)
    # C empty means the pin forces g = 1_{0}; the system is then feasible iff b is identically 0
    witness = [y for y in range(m) if bm[y] != 0]
    check('%s_pinned_polytope_of_M%d_E%s_is_infeasible' % (label, m, ''.join(map(str, e))),
          C == [] and witness != [],
          'E - E = Z_%d so the pin forces supp(g) = {0}; then (1_E * g)(y) = 0 at y in %s, not 1'
          % (m, witness))
print('')

# ---------------------------------------------------------------------------
print('=== Part F: the UNPINNED reading of item 3, disclosed rather than hidden')
m7, e7 = 7, [0, 1, 3]
A7 = circulant(m7, e7)
d7 = det(A7)
check('F1_the_Z7_circulant_of_the_Singer_set_has_determinant_plus_or_minus_24',
      abs(d7) == 24, 'det = %s; the incidence matrix of the Fano plane PG(2,2), and the symmetric '
      '2-(7,3,1) identity gives det^2 = (k-l)^(v-1)(k-l+lv) = 2^6 * 9 = 576 = 24^2' % d7)
r7, x7, u7 = solve(A7, [1] * m7)
check('F2_its_unique_solution_is_one_third_everywhere',
      u7 and x7 is not None and all(v == Fraction(1, 3) for v in x7),
      'nonsingular, so mu = 1/3 on all seven residues is the unique solution: nonnegative and '
      'FRACTIONAL')
check('F3_no_tiling_complement_exists_because_three_does_not_divide_seven',
      m7 % len(e7) != 0 and not any(tiles(m7, e7, T) for k in range(1, m7 + 1)
                                    for T in combinations(range(m7), k)),
      '|E| = 3 does not divide M = 7, and an exhaustive search over all 127 nonempty subsets of '
      'Z_7 finds no T with E (+) T = Z_7: the hull of item 3 is EMPTY, so the unpinned reading of '
      'item 3 is false trivially, already at M = 7')
print('')

# ---------------------------------------------------------------------------
print('=== Part G: the Z_12 seed, and why the failure is not visible at M = 12')
Ep = [0, 2, 3, 5]
Cp, Ap, bp = pinned_reduced_system(12, Ep)
check('G1_the_pin_leaves_three_free_variables_at_M12', Cp == [4, 6, 8],
      'C = Z_12 \\ (E\' - E\') = %s' % Cp)
rp, xp, up = solve(Ap, bp)
check('G2_the_pinned_polytope_of_Z12_is_a_single_integral_point',
      up and xp is not None and xp == [Fraction(1), Fraction(0), Fraction(1)],
      'unique solution x_4 = 1, x_6 = 0, x_8 = 1, i.e. g = 1_{{0,4,8}}: integral, so M = 12 '
      'carries no counterexample of this shape')
print('')

# ---------------------------------------------------------------------------
print('=== Part H: the infinite family -- item 3 is CHECKED to fail at M = 12d for d = 2..8 '
      '(M = 24, 36, 48, 60, 72, 84, 96); the step from these members to EVERY M divisible by 12 '
      'with M >= 24 rests on the block-diagonal / rank argument stated in the note below, which is '
      'NOT verified for general d by any check here')
M12 = circulant(12, Ep)
U = [0, 4, 8]
V = [3, 4, 5, 9, 10, 11]
check('H1_Eprime_plus_U_is_all_of_Z12',
      sorted(set((a + b) % 12 for a in Ep for b in U)) == list(range(12))
      and len(Ep) * len(U) == 12,
      "{0,2,3,5} (+) {0,4,8} = Z_12 is what this check computes; the dilation step to "
      'd.{0,2,3,5} (+) (d.{0,4,8} + {0..d-1}) = Z_{12d} for every d -- E always has a tiling '
      'complement through 0, so the hull is never empty -- is asserted, not computed here, and is '
      'checked directly only for d = 2..8, in H4')
check('H2_block_rank_on_the_integral_coset_columns', rank(cols(M12, U)) == 3,
      'the 12 x 3 submatrix of the Z_12 circulant on columns %s has rank 3' % U)
rV = rank(cols(M12, V))
dV = det([[M12[j][t] for t in V] for j in [6, 7, 8, 0, 1, 2]])
check('H3_block_rank_on_the_half_integral_coset_columns', rV == 6 and abs(dV) == 4,
      'the 12 x 6 submatrix on columns %s has rank 6, witnessed by the minor on rows [6,7,8,0,1,2] '
      'of determinant %s' % (V, dV))
note('the ARGUMENT of the family, asserted here and NOT executed for general d: the circulant of '
     'd.E\' in Z_{12d} is claimed to be BLOCK DIAGONAL over the d cosets of d.Z_{12d} with every '
     'block equal to the Z_12 circulant above, and the support of g_d is claimed to meet those '
     'blocks in exactly one U block and d-1 V blocks; GRANTING both of those claims, the checked '
     'ranks H2 and H3 give full column rank 3 + 6(d-1) for every d, which is the whole proof of '
     'the family. No check below proves the block decomposition or the support pattern for '
     'arbitrary d -- the rank equality itself is executed only for d = 2..8, in H4')

for d in range(2, 9):
    Md = 12 * d
    Ed = [d * e for e in Ep]
    gd = [Fraction(0)] * Md
    for j in U:
        gd[d * j] = Fraction(1)
    for c in range(1, d):
        for j in V:
            gd[d * j + c] = Fraction(1, 2)
    convd = convolution(Md, Ed, gd)
    Td = sorted(set((d * t + c) % Md for t in U for c in range(d)))
    Ad = circulant(Md, Ed)
    suppd = [t for t in range(Md) if gd[t] != 0]
    rd = rank(cols(Ad, suppd))
    ok = (gd[0] == 1 and all(v >= 0 for v in gd) and all(v == 1 for v in convd)
          and sum(gd, Fraction(0)) == 3 * d and any(v not in (0, 1) for v in gd)
          and rd == len(suppd) == 6 * d - 3
          and 0 in Td and tiles(Md, Ed, Td))
    check('H4_family_member_d%d_M%d' % (d, Md), ok,
          'E = %s; g_d in P_E, mass %s = 3d, half-integral, %d support columns of rank %d (a '
          'vertex), and T = %s is a tiling complement through 0'
          % (Ed, sum(gd, Fraction(0)), len(suppd), rd, Td if Md <= 36 else '(d.{0,4,8} + {0..d-1})'))
    if d == 2:
        check('H5_family_member_d2_is_exactly_the_object_of_Part_A', gd == g and Ed == E,
              'the d = 2 member of the family IS M = 24, E = {0,4,6,10} with the printed g')
print('')

# ---------------------------------------------------------------------------
print('=== Scope')
note('SCOPE: what is verified above is item 3 of the source\'s open problems under the source\'s '
     'OWN Definition -- a nonnegative f with f(0) = 1 and 1_E * f = 1 -- in a CYCLIC group. '
     'NOT RE-RUN and NOT CLAIMED here: (i) that M = 24 is the least modulus carrying such an '
     'object -- an exhaustive census over 8 <= M <= 24 was run elsewhere in this project and its '
     'raw transcript is not part of this folder, so minimality is left open here; (ii) the '
     'positive-definite variants asked in the same list, which this object cannot reach at all '
     '(check A10: g is not even symmetric, hence not positive definite); (iii) the two continuous '
     'Kolountzakis-Lev-'
     'Matolcsi versions, for convex bodies in R^d and for finite unions of intervals in R, which a '
     'finite cyclic object says nothing about; (iv) further seed orbits at moduli 18, 20 and 24 '
     'recorded in this project\'s artifacts but not reproduced here, so the family addressed above '
     'is exactly {12d : d >= 2} and no more -- and within it only the members d = 2..8 are '
     'checked, the rest resting on the block-diagonal argument granted in Part H\'s note; (v) whether a cyclic group contains a LONELY weak '
     'tile, i.e. a set with a weak tiling measure and no tiling complement at all, which is '
     'untouched by everything above.')
print('')

if FAILURES:
    print('VERDICT: %d CHECK(S) FAILED: %s' % (len(FAILURES), ', '.join(FAILURES)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % NPASS)
sys.exit(0)

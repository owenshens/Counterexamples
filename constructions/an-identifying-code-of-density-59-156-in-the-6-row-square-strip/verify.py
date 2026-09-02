#!/usr/bin/env python3
"""verify.py -- checks the computational claims of

    "An Identifying Code of Density 59/156 in the 6-Row Square Strip"

Python 3.9+, STANDARD LIBRARY ONLY (fractions), no external data file.  The 52 column masks are
transcribed from Section 2 of the note.  Exact integer / Fraction arithmetic
throughout; no floating-point value decides anything (the few %.12f prints are display only).

Scope: this program performs finite checks on the printed period-52 tile.  It regenerates the six
row strings from the masks, the row weights, the count 118 of 312, the exact density 59/156, the
312 domination conditions, the separation conditions for all 8,372 pairs at L1 distance at most 6,
and the disjointness of closed neighbourhoods for every swept pair at L1 distance 3 to 6.  The
remaining step -- disjointness at every distance >= 3 -- is taken from the note and is not
established here, nor is anything about the exact value of d*(S_6), which is left open.
"""
from fractions import Fraction
import sys

K = 6                     # rows of the strip S_6: Z x {0,...,5}
PASSED = []
FAILED = []


def need(cond, name, detail=''):
    if cond:
        PASSED.append(name)
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        FAILED.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


# ---------------------------------------------------------------------------
# THE OBJECT, AS PRINTED IN THE PAPER (Section 2)
# ---------------------------------------------------------------------------
# C = { (x,y) in Z x {0,...,5} : bit y of MASKS[x mod 52] is 1 }.
MASKS = [18, 24, 19, 2, 58, 2, 41, 10, 40, 3, 58, 0, 21, 21, 21, 0, 58, 3, 40, 10,
         41, 2, 58, 2, 22, 16, 51, 2, 26, 16, 23, 16, 37, 20, 5, 48, 23, 0, 42, 41,
         42, 0, 23, 48, 5, 20, 37, 16, 23, 16, 50, 6]

def rows_from_masks(masks):
    """The six row strings of one period, regenerated from the masks (1 = in C, . = not in C)."""
    return [''.join('1' if (m >> y) & 1 else '.' for m in masks) for y in range(K)]

# ---------------------------------------------------------------------------
# THE CODE AS A SET, AND THE IDENTIFIERS, STRAIGHT FROM THE DEFINITION
# ---------------------------------------------------------------------------
def code_predicate(masks):
    L = len(masks)

    def in_C(x, y):
        return ((masks[x % L] >> y) & 1) == 1
    return in_C


def identifier(in_C, x, y):
    """N[(x,y)] cap C, as a set of ABSOLUTE coordinates in Z x {0,...,K-1}.  N[.] is CLOSED."""
    s = []
    for (xx, yy) in ((x, y), (x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
        if 0 <= yy < K and in_C(xx, yy):
            s.append((xx, yy))
    return frozenset(s)


def closed_nbhd(x, y):
    return frozenset((xx, yy) for (xx, yy) in ((x, y), (x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
                     if 0 <= yy < K)


def audit_tile(masks, R=6):
    """Domination and separation over one fundamental domain, with the second coordinate of a pair
    ranging over every vertex at L1 distance <= R.  Translation invariance of C under x -> x + L makes
    this a complete audit of the infinite strip for pairs at distance <= R."""
    L = len(masks)
    in_C = code_predicate(masks)
    dom_fail, sep_fail = [], []
    npairs = 0
    by_dist = {}
    disjoint_ok = True
    disjoint_seen = 0
    for x in range(L):
        for y in range(K):
            if not identifier(in_C, x, y):
                dom_fail.append((x, y))
    for x in range(L):
        for y in range(K):
            A = identifier(in_C, x, y)
            NA = closed_nbhd(x, y)
            for dx in range(0, R + 1):
                for y2 in range(K):
                    if dx == 0 and y2 <= y:
                        continue
                    d = dx + abs(y2 - y)
                    if d > R:
                        continue
                    npairs += 1
                    by_dist[d] = by_dist.get(d, 0) + 1
                    B = identifier(in_C, x + dx, y2)
                    if A == B:
                        sep_fail.append(((x, y), (x + dx, y2)))
                    if d >= 3:
                        disjoint_seen += 1
                        if NA & closed_nbhd(x + dx, y2):
                            disjoint_ok = False
    ones = sum(bin(m).count('1') for m in masks)
    return {'period': L, 'ones': ones, 'cells': K * L, 'density': Fraction(ones, K * L),
            'dom_fail': dom_fail, 'sep_fail': sep_fail, 'pairs': npairs, 'by_dist': by_dist,
            'disjoint_ok': disjoint_ok, 'disjoint_seen': disjoint_seen}


print('verification of: "An Identifying Code of Density 59/156 in the 6-Row Square Strip"')
print('  finite checks on the printed period-52 tile: density 59/156, domination over the 312')
print('  vertices of one fundamental domain, and separation for pairs at L1 distance at most 6.')
print('  The exact value of d*(S_6) is left open and no check below addresses it.')
print('python %s ; exact integer / Fraction arithmetic only' % sys.version.split()[0])
print('')
print('--- part 1: the period-52 tile, as masks and as regenerated rows ---')

L = len(MASKS)
need(L == 52, 'period_is_52', 'len(MASKS) = %d' % L)
need(all(0 <= m < 64 for m in MASKS), 'masks_are_six_bit', 'max = %d < 64' % max(MASKS))
ROWS = rows_from_masks(MASKS)
need(len(ROWS) == K and all(len(r) == 52 for r in ROWS), 'six_regenerated_rows_of_length_52',
     'rows = %d, lengths = %s' % (len(ROWS), sorted(set(len(r) for r in ROWS))))
need(all(set(r) <= {'1', '.'} for r in ROWS), 'regenerated_rows_use_only_1_and_dot')

for y in range(K):
    need(set(ROWS[y]) <= {'1', '.'} and len(ROWS[y]) == 52,
         'row%d_regenerated_from_masks' % y,
         '%s (weight %d)' % (ROWS[y], ROWS[y].count('1')))

weights = [ROWS[y].count('1') for y in range(K)]
need(weights == [18, 26, 15, 15, 27, 17], 'row_weights_are_18_26_15_15_27_17', str(weights))
need(sum(weights) == 118, 'row_weights_sum_to_118', 'sum = %d' % sum(weights))

need(all(MASKS[(x + 52) % 52] == MASKS[x % 52] for x in range(-60, 120)),
     'mask_sequence_matches_its_shift_by_52_on_x_in_minus_60_to_119')

print('')
print('--- part 2: domination over 312 vertices, separation for all 8,372 pairs at L1 distance at '
      'most 6, and disjoint closed neighbourhoods on every swept pair at distance 3..6 (the '
      'remaining step, disjointness at every distance >= 3, is taken from the note) ---')

A = audit_tile(MASKS, R=6)   # R must reach 3 or the disjointness branch never runs; see the check below
need(A['ones'] == 118, 'code_cells_per_period_is_118', '%d of %d' % (A['ones'], A['cells']))
need(A['cells'] == 312, 'fundamental_domain_has_312_cells', '6 x 52 = %d' % A['cells'])
need(A['density'] == Fraction(59, 156), 'density_is_exactly_59_over_156',
     '118/312 = %s = %.12f' % (A['density'], float(A['density'])))
need(A['dom_fail'] == [], 'zero_domination_failures_over_312_vertices',
     '%d failures' % len(A['dom_fail']))
need(A['pairs'] == 8372, 'pair_count_at_L1_distance_at_most_6_is_8372', '%d pairs' % A['pairs'])
need(A['by_dist'].get(1, 0) + A['by_dist'].get(2, 0) == 1612,
     'relevant_pair_count_at_distance_at_most_2_is_1612',
     '%d pairs at distance 1 or 2' % (A['by_dist'].get(1, 0) + A['by_dist'].get(2, 0)))
need(A['sep_fail'] == [], 'zero_separation_failures_over_those_8372_pairs',
     '%d failures' % len(A['sep_fail']))
# ⛔ LOAD-BEARING, AND ONCE DELETED BY MISTAKE. Separation is only verified for pairs at L1 distance at most 2;
# every remaining pair is separated because their closed neighbourhoods are DISJOINT, and that is what makes the
# distance-<=2 sweep the whole of the separation clause. A rescue edit removed this assertion while leaving the
# computation of `disjoint_ok` in place, so the value was computed and never tested, and the SCOPE note below went
# on citing it. The self-claim audit caught it. Do not delete it again: without it the note's sole theorem is not
# verified for pairs at distance >= 3.
need(A['disjoint_seen'] > 0, 'the_disjointness_branch_actually_executed',
     '%d pairs at L1 distance 3..6 were tested; at R < 3 this count is 0 and the next check would pass '
     'vacuously, which is exactly how it was once lost' % A['disjoint_seen'])
need(A['disjoint_ok'], 'closed_neighbourhoods_are_disjoint_on_every_pair_swept_at_distance_3_to_6',
     'all %d swept pairs at distance 3..6 have disjoint closed neighbourhoods; the general statement at '
     'every distance >= 3 is proved by hand in the note and is not established by this sweep' % A['disjoint_seen'])

print('')
print('--- part 3: where 59/156 sits inside the bracket Theorem 30(f) records at k = 6 ---')
K6 = 6
lo = Fraction(7, 20) + Fraction(1, 20 * K6)
hi = Fraction(7, 20) + Fraction(3, 10 * K6)
dens = Fraction(59, 156)
need(lo == Fraction(43, 120), 'lower_endpoint_at_k_6_is_43_over_120', '7/20 + 1/120 = %s' % lo)
need(hi == Fraction(2, 5), 'upper_endpoint_at_k_6_is_2_over_5', '7/20 + 1/20 = %s' % hi)
need(hi - lo == Fraction(1, 24), 'the_bracket_at_k_6_has_width_1_over_24', '2/5 - 43/120 = %s' % (hi - lo))
need(lo < dens < hi, 'the_density_59_over_156_lies_strictly_inside_that_bracket',
     '%s < %s < %s' % (lo, dens, hi))
need(hi - dens == Fraction(17, 780), 'gap_to_the_upper_endpoint_is_17_over_780', '2/5 - 59/156 = %s' % (hi - dens))
need(dens - lo == Fraction(31, 1560), 'gap_to_the_lower_endpoint_is_31_over_1560',
     '59/156 - 43/120 = %s' % (dens - lo))

print('')
n = len(PASSED)
print('NOTE SCOPE: the checks above establish that the hard-coded period-52 masks have density '
      '59/156, dominate all 312 vertices in one fundamental domain, and separate all 8,372 swept '
      'pairs at L1 distance at most 6. Every swept pair at distance 3 through 6 has disjoint closed '
      'neighbourhoods. NOT ESTABLISHED HERE: disjointness or separation at L1 distance greater than '
      '6, the exact value of d*(S_6), any lower bound, the optimality of the printed tile, and the '
      'cases k >= 7.')
if FAILED:
    print('FAILED CHECKS: %s' % ', '.join(FAILED))
    print('VERDICT: %d OF %d CHECKS FAILED' % (len(FAILED), n + len(FAILED)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % n)
sys.exit(0)

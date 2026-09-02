#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verification program for

    "A Cyclic Completely Quasi-Uniform Nested SQS(22)"

WHAT IT READS.  Nothing but the 19 lines printed in Section 2 of the paper, pasted verbatim into
BASE_BLOCKS_AS_PRINTED below, and the integers the paper itself prints.  No external data file, no
third-party module.  Python 3.9+, standard library only.  All arithmetic is over the integers and
Fraction; no floating-point value is computed, so no decision depends on one.

WHAT IT CHECKS.  From the printed table, the program checks that the entries are valid nested
blocks over Z_22, computes the stated orbit lengths and their union of 385 nested blocks, verifies
that the 385 underlying quadruples are distinct and cover every one of the 1540 triples exactly
once, checks that all 231 pairs occur as halves with multiplicity profile {3: 154, 4: 77}, and
verifies invariance under x -> x+1.

Output contract: one `PASS <name> [detail]` line per check, a closing
`VERDICT: ALL <n> CHECKS PASS`, and exit status 0 if and only if every check passed.
"""

import itertools
import re
import sys
from collections import Counter

# ----------------------------------------------------------------------------------------------
# THE OBJECT, EXACTLY AS PRINTED IN SECTION 2 OF THE PAPER
# ----------------------------------------------------------------------------------------------
BASE_BLOCKS_AS_PRINTED = """
B01  { 0,  1 |  2,  4}   22        B11  { 0,  2 | 12, 19}   22
B02  { 0,  1 |  5,  6}   22        B12  { 0, 18 |  2,  6}   22
B03  { 0,  7 |  1,  8}   22        B13  { 0, 18 |  3,  9}   22
B04  { 0,  1 |  9, 19}   22        B14  { 0, 11 |  3, 14}   11
B05  { 0, 10 |  1, 11}   22        B15  { 0,  8 |  4, 15}   22
B06  { 0, 14 |  1, 20}   22        B16  { 0,  9 |  2,  7}   22
B07  { 0,  2 |  8, 14}   22        B17  { 0, 13 |  2, 11}   11
B08  { 0,  3 |  7, 17}   22        B18  { 0, 17 |  4,  9}   22
B09  { 0, 13 |  3,  6}   22        B19  { 0, 16 |  5, 11}   11
B10  { 0,  5 |  2, 10}   22
"""

V = 22

# The diameter decompositions printed in the proof of Lemma 2.
PRINTED_SHORT_ORBIT_DECOMPOSITIONS = {
    'B14': ((0, 11), (3, 14)),
    'B17': ((0, 11), (2, 13)),
    'B19': ((0, 11), (5, 16)),
}

# Integers printed in the paper, listed here so each comparison can actually fail.
PAPER = {
    'base_blocks': 19,
    'long_orbits': 16,
    'short_orbits': 3,
    'blocks': 385,
    'triples': 1540,
    'pairs': 231,
    'mult_lo': 3,
    'mult_hi': 4,
    'n_lo': 154,
    'n_hi': 77,
    'mult_sum': 770,
}

_CHECKS = []
_FAILED = []


def check(name, ok, detail=''):
    """Record one check. ⛔ A PASS line is printed only when ok is true; a FAIL line otherwise, and
    the process exits nonzero at the end. The verdict count is len(_CHECKS), i.e. the number of
    PASS lines actually emitted -- it cannot drift from them."""
    if ok:
        _CHECKS.append(name)
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILED.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


# ----------------------------------------------------------------------------------------------
# Step 1 -- parse the printed table
# ----------------------------------------------------------------------------------------------
print('verification of "A Cyclic Completely Quasi-Uniform Nested SQS(22)"')
print('python %s, standard library only, exact integer arithmetic only'
      % sys.version.split()[0])
print('the object is read from the 19 lines printed in Section 2 of the paper')
print('')
print('--- Step 1: the printed table of base nested blocks ---')

ROW_RE = re.compile(r'(B\d\d)\s+\{\s*(\d+)\s*,\s*(\d+)\s*\|\s*(\d+)\s*,\s*(\d+)\s*\}\s+(\d+)')

parsed = []
for m in ROW_RE.finditer(BASE_BLOCKS_AS_PRINTED):
    label = m.group(1)
    a, b, c, d = (int(m.group(i)) for i in (2, 3, 4, 5))
    parsed.append((label, ((a, b), (c, d)), int(m.group(6))))

check('parsed_19_base_nested_blocks_from_the_printed_table',
      len(parsed) == PAPER['base_blocks'],
      'parsed %d rows, paper prints %d' % (len(parsed), PAPER['base_blocks']))

labels = [p[0] for p in parsed]
check('labels_are_B01_through_B19_each_exactly_once',
      sorted(labels) == ['B%02d' % i for i in range(1, 20)],
      '19 distinct labels B01..B19 (the table is printed in two columns, so the reading order is '
      'B01,B11,B02,...)')

ok_points = all(0 <= x < V for _l, halves, _L in parsed for h in halves for x in h)
check('every_printed_point_lies_in_Z_22', ok_points, 'all 76 entries in 0..21')

ok_4set = all(len({halves[0][0], halves[0][1], halves[1][0], halves[1][1]}) == 4
              for _l, halves, _L in parsed)
check('every_base_block_is_a_4_set_split_into_two_disjoint_pairs', ok_4set,
      '19 blocks, 4 distinct points each, halves disjoint')

check('printed_orbit_lengths_are_16_twenty_twos_and_3_elevens',
      sorted(Counter(L for _l, _h, L in parsed).items()) == [(11, 3), (22, 16)],
      'multiset of printed lengths = {22: 16, 11: 3}')

check('printed_orbit_lengths_sum_to_385',
      sum(L for _l, _h, L in parsed) == PAPER['blocks'],
      '16*22 + 3*11 = %d' % sum(L for _l, _h, L in parsed))


# ----------------------------------------------------------------------------------------------
# Step 2 -- Lemma 2, the orbits
# ----------------------------------------------------------------------------------------------
def shift_nb(nb, t):
    (a, b), (c, d) = nb
    h1 = tuple(sorted(((a + t) % V, (b + t) % V)))
    h2 = tuple(sorted(((c + t) % V, (d + t) % V)))
    return tuple(sorted((h1, h2)))


def norm(nb):
    return shift_nb(nb, 0)


print('')
print('--- Step 2: Lemma 2, the orbits under x -> x+1 ---')

orbits = []
for label, halves, printed_len in parsed:
    orb = set(shift_nb(halves, t) for t in range(V))
    orbits.append((label, norm(halves), orb, printed_len))

check('computed_orbit_length_equals_the_printed_length_for_all_19_blocks',
      all(len(orb) == printed_len for _l, _n, orb, printed_len in orbits),
      'orbit lengths ' + ','.join(str(len(o)) for _l, _n, o, _p in orbits))

check('every_orbit_length_divides_22_and_is_22_or_11',
      all(len(orb) in (11, 22) for _l, _n, orb, _p in orbits),
      'stabilisers are trivial or {0,11}')

short = sorted(l for l, _n, orb, _p in orbits if len(orb) == 11)
check('the_three_short_orbits_are_exactly_B14_B17_B19', short == ['B14', 'B17', 'B19'],
      'short-orbit labels = %s' % ','.join(short))


def is_diameter_union(halves):
    q = set(halves[0]) | set(halves[1])
    return all(((x + 11) % V) in q for x in q)


check('a_base_block_has_a_short_orbit_iff_its_quadruple_is_a_union_of_two_diameters',
      all((len(orb) == 11) == is_diameter_union(halves)
          for (label, halves, _L), (_l2, _n, orb, _p) in zip(parsed, orbits)),
      'the biconditional of Lemma 2 holds on all 19 rows')

dec_ok = True
for label, (d1, d2) in PRINTED_SHORT_ORBIT_DECOMPOSITIONS.items():
    halves = dict((l, h) for l, h, _L in parsed)[label]
    q = set(halves[0]) | set(halves[1])
    dec_ok &= (set(d1) | set(d2) == q
               and (d1[0] + 11) % V == d1[1] and (d2[0] + 11) % V == d2[1])
check('the_diameter_decompositions_printed_in_the_proof_of_Lemma_2_are_correct', dec_ok,
      'B14 = {0,11}+{3,14}, B17 = {0,11}+{2,13}, B19 = {0,11}+{5,16}')

check('all_three_splittings_of_a_diameter_union_quadruple_are_sigma_invariant',
      all(shift_nb(sp, 11) == norm(sp)
          for a in range(V) for b in range(V)
          if a < b and b != (a + 11) % V
          for sp in (((a, (a + 11) % V), (b, (b + 11) % V)),
                     ((a, b), ((a + 11) % V, (b + 11) % V)),
                     ((a, (b + 11) % V), ((a + 11) % V, b)))),
      'checked for every unordered pair of distinct diameters, all 3 splittings')

B = set()
for _l, _n, orb, _p in orbits:
    B |= orb
check('the_19_orbits_are_pairwise_disjoint_and_give_385_nested_blocks',
      len(B) == PAPER['blocks'] and len(B) == sum(len(o) for _l, _n, o, _p in orbits),
      '|B| = %d' % len(B))

check('the_nested_block_set_is_invariant_under_x_to_x_plus_1',
      all(shift_nb(nb, 1) in B for nb in B), 'B is a union of Z_22-orbits, verified elementwise')


# ----------------------------------------------------------------------------------------------
# Step 3 -- Theorem 1: the underlying quadruples form an SQS(22)
# ----------------------------------------------------------------------------------------------
print('')
print('--- Step 3: Theorem 1, the underlying quadruples form an SQS(22) ---')

quads = [tuple(sorted(set(h1) | set(h2))) for (h1, h2) in B]
check('the_385_nested_blocks_have_385_distinct_underlying_quadruples',
      len(set(quads)) == len(quads) == PAPER['blocks'],
      '%d nested blocks, %d distinct quadruples' % (len(quads), len(set(quads))))

Q = set(quads)
tri = Counter()
for q in Q:
    for T in itertools.combinations(q, 3):
        tri[T] += 1
check('the_number_of_triples_of_a_22_set_is_1540',
      len(list(itertools.combinations(range(V), 3))) == PAPER['triples'],
      'C(22,3) = %d' % PAPER['triples'])
check('every_one_of_the_1540_triples_is_covered_exactly_once',
      len(tri) == PAPER['triples'] and set(tri.values()) == {1},
      '%d triples covered, multiplicity set %s' % (len(tri), sorted(set(tri.values()))))


# ----------------------------------------------------------------------------------------------
# Step 4 -- Theorem 1: complete, quasi-uniform, not uniform
# ----------------------------------------------------------------------------------------------
print('')
print('--- Step 4: Theorem 1, completeness and the quasi-uniform profile ---')

mu = Counter()
for (h1, h2) in B:
    mu[h1] += 1
    mu[h2] += 1

check('the_number_of_pairs_of_a_22_set_is_231',
      len(list(itertools.combinations(range(V), 2))) == PAPER['pairs'], 'C(22,2) = 231')
check('all_231_pairs_occur_as_a_half_so_the_nesting_is_complete',
      len(mu) == PAPER['pairs'] and min(mu.values()) >= 1,
      '%d distinct pairs occur, minimum multiplicity %d' % (len(mu), min(mu.values())))

profile = Counter(mu.values())
check('the_multiplicity_profile_is_154_pairs_at_3_and_77_pairs_at_4',
      dict(profile) == {PAPER['mult_lo']: PAPER['n_lo'], PAPER['mult_hi']: PAPER['n_hi']},
      'profile = %s' % dict(sorted(profile.items())))
check('the_nesting_is_quasi_uniform_and_not_uniform',
      len(profile) == 2 and max(mu.values()) - min(mu.values()) == 1,
      'two values, %d and %d, differing by 1' % (min(mu.values()), max(mu.values())))
check('the_multiplicity_sum_is_770_equals_twice_385',
      sum(mu.values()) == PAPER['mult_sum'] == 2 * PAPER['blocks'],
      '154*3 + 77*4 = %d = 2*385' % sum(mu.values()))

# ----------------------------------------------------------------------------------------------
# scope and verdict
# ----------------------------------------------------------------------------------------------
print('')
print('NOTE SCOPE: this program checks the single completely quasi-uniform nested SQS(22) '
      'printed in the note: its orbit lengths, 385 distinct underlying quadruples, unique '
      'coverage of all 1540 triples, completeness, multiplicity profile {3: 154, 4: 77}, and '
      'invariance under x -> x+1. It does not test isomorphism, enumerate other systems, or '
      'determine the full automorphism group.')
print('')
if _FAILED:
    print('VERDICT: %d CHECK(S) FAILED: %s' % (len(_FAILED), ', '.join(_FAILED)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(_CHECKS))
sys.exit(0)

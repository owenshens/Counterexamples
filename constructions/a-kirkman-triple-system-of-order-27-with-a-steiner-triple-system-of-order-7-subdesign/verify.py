#!/usr/bin/env python3
"""Verification of

    "A Kirkman Triple System of Order 27 with a Steiner Triple System of Order 7
     as a Subdesign"

Python 3.9+, STANDARD LIBRARY ONLY, exact integer arithmetic only: no
third-party package, no external data file, no floating-point value anywhere, so
no decision below depends on one.

The only input is the text pasted into TABLE below -- copied by hand from Table 1
of the paper -- and re-parsed from that text by this program.  The quantities
listed in the closing SCOPE note are then re-derived and compared, so every
comparison can fail.  That list was inventoried by hand from the paper when this
program was written: no check below opens paper.tex, so neither the transcription
of Table 1 nor the completeness of the list is established here.

One PASS line per check; a closing VERDICT line; exit 0 iff every check passed.
"""

import itertools
import sys

# ----------------------------------------------------------------------------
# THE OBJECT: Table 1 of the paper, copied by hand into the string below.
# Nothing here re-reads paper.tex, so the transcription itself is assumed, not
# checked.  A line beginning `Cnn:` opens a class; a following line without a
# label continues it.  Triples are separated by `|`, and a trailing `|` marks
# the continuation.
# ----------------------------------------------------------------------------
TABLE = """
C00: 0 1 3 | 2 17 18 | 4 11 21 | 5 22 23 | 6 12 19 |
     7 9 14 | 8 15 20 | 10 13 24 | 16 25 26
C01: 0 7 8 | 1 2 4 | 3 20 24 | 5 19 26 | 6 18 25 |
     9 10 11 | 12 15 22 | 13 16 23 | 14 17 21
C02: 0 11 16 | 1 19 25 | 2 3 5 | 4 9 20 | 6 7 24 |
     8 21 26 | 10 14 22 | 12 13 18 | 15 17 23
C03: 0 9 17 | 1 10 18 | 2 20 26 | 3 4 6 | 5 7 25 |
     8 22 24 | 11 12 23 | 13 14 19 | 15 16 21
C04: 0 4 5 | 1 16 20 | 2 10 23 | 3 14 18 | 6 21 22 |
     7 11 13 | 8 17 19 | 9 12 26 | 15 24 25
C05: 0 10 15 | 1 5 6 | 2 11 19 | 3 7 26 | 4 18 24 |
     8 23 25 | 9 13 21 | 12 14 20 | 16 17 22
C06: 0 2 6 | 1 9 22 | 3 21 23 | 4 15 19 | 5 13 20 |
     7 10 12 | 8 16 18 | 11 14 25 | 17 24 26
C07: 0 19 23 | 1 15 26 | 2 7 22 | 3 13 17 | 4 8 12 |
     5 14 24 | 6 9 16 | 10 21 25 | 11 18 20
C08: 0 20 21 | 1 8 13 | 2 16 24 | 3 12 25 | 4 7 23 |
     5 10 17 | 6 14 15 | 9 18 19 | 11 22 26
C09: 0 18 22 | 1 7 21 | 2 8 14 | 3 11 15 | 4 17 25 |
     5 12 16 | 6 13 26 | 9 23 24 | 10 19 20
C10: 0 14 26 | 1 12 17 | 2 9 25 | 3 10 16 | 4 13 22 |
     5 8 11 | 6 20 23 | 7 15 18 | 19 21 24
C11: 0 12 24 | 1 14 23 | 2 13 15 | 3 8 9 | 4 10 26 |
     5 18 21 | 6 11 17 | 7 16 19 | 20 22 25
C12: 0 13 25 | 1 11 24 | 2 12 21 | 3 19 22 | 4 14 16 |
     5 9 15 | 6 8 10 | 7 17 20 | 18 23 26
"""

V = 27                                   # the order, as claimed
U = frozenset(range(7))                  # the subdesign's point set, as printed

# The seven blocks the paper's proof lists as the ones inside U, in the order it
# prints them, and the classes it says they lie in.
INNER_PRINTED = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6),
                 (0, 4, 5), (1, 5, 6), (0, 2, 6)]
INNER_CLASS_PRINTED = [0, 1, 2, 3, 4, 5, 6]

# sigma, as printed in Section 4, in cycle form; every point not named is fixed.
SIGMA_CYCLES = [(1, 2, 4), (3, 6, 5), (9, 10, 11), (12, 13, 14), (15, 16, 17),
                (18, 19, 20), (21, 22, 23), (24, 25, 26)]
SIGMA_FIXED_PRINTED = (0, 7, 8)

# pi, the induced permutation of the classes, as printed in Section 4.
PI_CYCLES = [(0, 6, 4), (2, 3, 5), (7, 8, 9), (10, 11, 12)]
PI_FIXED_PRINTED = (1,)

# the five classes Section 4 says suffice
BASE_CLASSES = (0, 1, 2, 7, 10)

# admissible orders of a Steiner triple system (1 or 3 mod 6), up to 27
ADMISSIBLE = [u for u in range(1, V + 1) if u % 6 in (1, 3)]


# ----------------------------------------------------------------------------
# check bookkeeping
# ----------------------------------------------------------------------------
_n_pass = 0
_n_fail = 0


def ok(name, detail=''):
    global _n_pass
    _n_pass += 1
    print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))


def bad(name, detail=''):
    global _n_fail
    _n_fail += 1
    print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def check(name, cond, detail=''):
    if cond:
        ok(name, detail)
    else:
        bad(name, detail)
    return bool(cond)


def note(s):
    print('NOTE %s' % s)


def perm_from_cycles(cycles, n):
    p = list(range(n))
    for c in cycles:
        for i, x in enumerate(c):
            p[x] = c[(i + 1) % len(c)]
    return p


# ----------------------------------------------------------------------------
print('verification of "A Kirkman Triple System of Order 27 with a Steiner '
      'Triple System of Order 7 as a Subdesign"')
print('python %d.%d.%d, standard library only, exact integer arithmetic only'
      % sys.version_info[:3])
print('the only input is the text pasted into TABLE below, copied by hand from '
      'Table 1 of the paper and re-parsed here; no check reads paper.tex')
print('')

# ============================================================================
print('=== Step 1: parsing the printed table')
# ============================================================================
lines = [ln for ln in TABLE.splitlines() if ln.strip()]
check('table_has_26_lines_two_per_class', len(lines) == 26, '%d lines' % len(lines))

labels = []
raw = []          # accumulated text per class
for ln in lines:
    head = ln.split(':', 1)
    if len(head) == 2 and head[0].strip().startswith('C') and head[0].strip()[1:].isdigit():
        labels.append(head[0].strip())
        raw.append(head[1])
    else:
        if not raw:
            bad('continuation_line_before_any_label', ln)
            print('VERDICT: ABORTED')
            sys.exit(1)
        raw[-1] += ' ' + ln

check('labels_are_C00_through_C12_in_order',
      labels == ['C%02d' % i for i in range(13)], ' '.join(labels))

classes = []
for text in raw:
    trips = []
    for piece in text.split('|'):
        piece = piece.strip()
        if not piece:
            continue
        trips.append(tuple(int(t) for t in piece.split()))
    classes.append(trips)

check('every_class_parses_to_nine_triples',
      all(len(c) == 9 for c in classes),
      'sizes %s' % sorted(set(len(c) for c in classes)))

flat = [b for c in classes for b in c]
check('every_triple_has_three_distinct_points_in_0_26',
      all(len(b) == 3 and len(set(b)) == 3 and all(0 <= x < V for x in b)
          for b in flat),
      '%d triples' % len(flat))
check('every_triple_is_printed_in_increasing_order',
      all(list(b) == sorted(b) for b in flat))
check('total_block_count_is_117', len(flat) == 117, '13*9 = %d' % (13 * 9))
check('the_117_blocks_are_pairwise_distinct',
      len(set(flat)) == len(flat), '%d distinct' % len(set(flat)))

blocks = [tuple(sorted(b)) for b in flat]
blockset = set(blocks)
cls_of = {}
for ci, c in enumerate(classes):
    for b in c:
        cls_of[tuple(sorted(b))] = ci
print('')

# ============================================================================
print('=== Step 2: Theorem 1(a), the Kirkman triple system')
# ============================================================================
check('each_of_the_13_classes_is_a_partition_of_the_27_points',
      all(sorted(x for b in c for x in b) == list(range(V)) for c in classes),
      '13 classes of 9 disjoint triples')

pairs = {}
for b in blocks:
    for p in itertools.combinations(b, 2):
        pairs[p] = pairs.get(p, 0) + 1
check('pair_slots_number_351', sum(pairs.values()) == 351,
      '117 blocks * 3 pairs = %d' % (117 * 3))
check('arithmetic_binom_27_2_is_351', V * (V - 1) // 2 == 351,
      '27*26/2 = %d' % (V * (V - 1) // 2))
check('every_pair_of_points_is_covered_exactly_once',
      len(pairs) == 351 and set(pairs.values()) == {1},
      '%d distinct pairs, multiplicities %s'
      % (len(pairs), sorted(set(pairs.values()))))
check('the_design_is_an_STS_27',
      len(pairs) == V * (V - 1) // 2 and set(pairs.values()) == {1}
      and all(len(b) == 3 for b in blocks))
check('the_13_classes_form_a_resolution_so_it_is_a_KTS_27',
      sum(len(c) for c in classes) == len(blockset)
      and all(sorted(x for b in c for x in b) == list(range(V)) for c in classes))
check('number_of_classes_equals_v_minus_1_over_2', len(classes) == (V - 1) // 2,
      '(27-1)/2 = %d' % ((V - 1) // 2))

rep = {}
for b in blocks:
    for x in b:
        rep[x] = rep.get(x, 0) + 1
check('every_point_lies_in_exactly_13_blocks',
      set(rep) == set(range(V)) and set(rep.values()) == {13},
      'replication %s' % sorted(set(rep.values())))
check('every_point_lies_in_exactly_one_block_of_every_class',
      all(sum(1 for b in c if x in b) == 1 for c in classes for x in range(V)))
print('')

# ============================================================================
print('=== Step 3: Theorem 1(b), the STS(7) subdesign on U = {0,...,6}')
# ============================================================================
inner = sorted(b for b in blocks if set(b) <= U)
check('exactly_seven_blocks_lie_inside_U', len(inner) == 7,
      '%d blocks' % len(inner))
note('blocks inside U = %s' % (inner,))
check('the_seven_inner_blocks_are_the_ones_printed_in_the_theorem',
      inner == sorted(tuple(sorted(b)) for b in INNER_PRINTED))

ip = {}
for b in inner:
    for p in itertools.combinations(b, 2):
        ip[p] = ip.get(p, 0) + 1
check('inner_blocks_cover_all_21_pairs_of_U_exactly_once',
      len(ip) == 21 and set(ip.values()) == {1},
      '7*6/2 = %d pairs, multiplicities %s' % (7 * 6 // 2, sorted(set(ip.values()))))
check('U_with_its_seven_blocks_is_an_STS_7',
      len(U) == 7 and len(inner) == 7 and len(ip) == len(U) * (len(U) - 1) // 2
      and set(ip.values()) == {1})
check('every_point_of_U_lies_on_three_inner_blocks',
      all(sum(1 for b in inner if x in b) == 3 for x in sorted(U)),
      'r = (7-1)/2 = 3')
check('any_two_inner_blocks_meet_in_exactly_one_point',
      all(len(set(a) & set(b)) == 1
          for a, b in itertools.combinations(inner, 2)),
      'the Fano plane, 21 pairs of lines')
check('the_seven_inner_blocks_lie_one_each_in_C00_to_C06_in_the_printed_order',
      [cls_of[tuple(sorted(b))] for b in INNER_PRINTED] == INNER_CLASS_PRINTED,
      'classes %s' % [cls_of[tuple(sorted(b))] for b in INNER_PRINTED])
check('no_block_meets_U_in_exactly_two_points',
      not any(len(set(b) & U) == 2 for b in blocks),
      'forced: such a pair is already covered inside U')
print('')

# ============================================================================
print('=== Step 4: the 7 + 70 + 40 census of the Remark')
# ============================================================================
prof = []
for c in classes:
    a = sum(1 for b in c if set(b) <= U)
    cr = sum(1 for b in c if len(set(b) & U) == 1)
    outr = sum(1 for b in c if not (set(b) & U))
    prof.append((a, cr, outr))
note('per-class (inner, one-point, disjoint) profiles = %s' % (prof,))
check('classes_C00_to_C06_have_profile_one_inner_four_cross_four_outer',
      all(p == (1, 4, 4) for p in prof[:7]))
check('classes_C07_to_C12_have_profile_no_inner_seven_cross_two_outer',
      all(p == (0, 7, 2) for p in prof[7:]))
tot = tuple(sum(p[i] for p in prof) for i in range(3))
check('census_splits_117_as_7_plus_70_plus_40', tot == (7, 70, 40),
      'totals %s, sum %d' % (tot, sum(tot)))
check('cross_block_count_is_seven_points_times_ten_blocks_each',
      all(sum(1 for b in blocks if set(b) & U == {x}) == 10 for x in sorted(U))
      and 7 * 10 == tot[1], '13 blocks per point, 3 of them inner')
print('')

# ============================================================================
print('=== Step 5: sigma, the prescribed automorphism of Section 4')
# ============================================================================
sigma = perm_from_cycles(SIGMA_CYCLES, V)
check('sigma_is_a_permutation_of_the_27_points',
      sorted(sigma) == list(range(V)))
check('sigma_has_order_3',
      all(sigma[sigma[sigma[x]]] == x for x in range(V))
      and any(sigma[x] != x for x in range(V)))
fix = tuple(x for x in range(V) if sigma[x] == x)
check('sigma_fixes_exactly_the_points_0_7_8', fix == SIGMA_FIXED_PRINTED,
      'fixed points %s' % (fix,))
cyc_lens = []
seen = set()
for x in range(V):
    if x in seen:
        continue
    orb = [x]
    seen.add(x)
    y = sigma[x]
    while y != x:
        orb.append(y)
        seen.add(y)
        y = sigma[y]
    cyc_lens.append(len(orb))
check('sigma_orbit_shape_is_three_fixed_points_and_eight_3_cycles',
      sorted(cyc_lens) == [1, 1, 1] + [3] * 8,
      'orbit lengths %s' % sorted(cyc_lens))
check('sigma_preserves_the_block_set',
      all(tuple(sorted(sigma[x] for x in b)) in blockset for b in blockset))
check('the_fixed_set_0_7_8_is_a_block_of_the_design',
      tuple(SIGMA_FIXED_PRINTED) in blockset)
check('that_block_lies_in_class_C01',
      cls_of[tuple(SIGMA_FIXED_PRINTED)] == 1,
      'C%02d' % cls_of[tuple(SIGMA_FIXED_PRINTED)])
check('sigma_preserves_U_setwise', set(sigma[x] for x in U) == set(U))
check('sigma_on_U_is_x_maps_to_2x_mod_7',
      all(sigma[x] == (2 * x) % 7 for x in sorted(U)),
      'a collineation of the Fano plane')

pi = {}
consistent = True
for ci, c in enumerate(classes):
    img = set(cls_of[tuple(sorted(sigma[x] for x in b))] for b in c)
    if len(img) != 1:
        consistent = False
        break
    pi[ci] = img.pop()
check('sigma_maps_each_class_onto_a_single_class', consistent)
pi_printed = perm_from_cycles(PI_CYCLES, 13)
check('the_induced_class_permutation_is_the_printed_pi',
      consistent and [pi[i] for i in range(13)] == pi_printed,
      'pi = %s' % ([pi[i] for i in range(13)] if consistent else '?'))
check('the_printed_pi_fixes_only_C01',
      tuple(i for i in range(13) if pi_printed[i] == i) == PI_FIXED_PRINTED)
check('the_induced_class_permutation_has_order_3',
      all(pi_printed[pi_printed[pi_printed[i]]] == i for i in range(13))
      and any(pi_printed[i] != i for i in range(13)))
print('')

# ============================================================================
print('=== Step 6: the compact form -- sigma plus five base classes')
# ============================================================================
orbits = []
seen = set()
for i in range(13):
    if i in seen:
        continue
    o = [i]
    seen.add(i)
    j = pi_printed[i]
    while j != i:
        o.append(j)
        seen.add(j)
        j = pi_printed[j]
    orbits.append(tuple(o))
note('orbits of pi on the 13 classes = %s' % (orbits,))
check('the_five_named_classes_meet_every_orbit_of_pi_exactly_once',
      len(orbits) == 5
      and all(len(set(o) & set(BASE_CLASSES)) == 1 for o in orbits),
      'base classes %s' % (tuple('C%02d' % b for b in BASE_CLASSES),))

rebuilt = {}
for b in BASE_CLASSES:
    cur = [tuple(sorted(t)) for t in classes[b]]
    idx = b
    for _ in range(3):
        if idx in rebuilt:
            break
        rebuilt[idx] = sorted(cur)
        cur = [tuple(sorted(sigma[x] for x in t)) for t in cur]
        idx = pi_printed[idx]
check('sigma_and_the_five_base_classes_rebuild_all_13_classes',
      sorted(rebuilt) == list(range(13)),
      '%d classes rebuilt' % len(rebuilt))
check('the_rebuilt_design_equals_the_printed_design_class_by_class',
      all(rebuilt.get(i) == sorted(tuple(sorted(t)) for t in classes[i])
          for i in range(13)))
print('')

# ============================================================================
print('=== Step 7: Proposition 3, the subdesign census')
# ============================================================================
third = {}
for b in blocks:
    x, y, z = b
    third[(x, y)] = z
    third[(x, z)] = y
    third[(y, z)] = x


def closure(seed):
    s = set(seed)
    growing = True
    while growing:
        growing = False
        for p, q in itertools.combinations(sorted(s), 2):
            t = third[(p, q)] if p < q else third[(q, p)]
            if t not in s:
                s.add(t)
                growing = True
                break
    return frozenset(s)


check('arithmetic_binom_27_3_is_2925',
      len(list(itertools.combinations(range(V), 3))) == 2925,
      '27*26*25/6 = %d' % (V * (V - 1) * (V - 2) // 6))

closed = set()
for tri in itertools.combinations(range(V), 3):
    closed.add(closure(tri))
check('closure_of_every_point_triple_is_a_closed_set',
      all(all((third[(p, q)] in s) for p, q in itertools.combinations(sorted(s), 2))
          for s in closed),
      '%d distinct closed sets from 2925 triples' % len(closed))
sizes = {}
for s in closed:
    sizes[len(s)] = sizes.get(len(s), 0) + 1
note('closed-set sizes found (size: count) = %s' % sorted(sizes.items()))
check('every_closed_set_of_size_at_least_three_is_a_subdesign',
      all(sum(1 for b in blocks if set(b) <= s) == len(s) * (len(s) - 1) // 6
          for s in closed if len(s) >= 3))

s7 = [s for s in closed if len(s) == 7]
check('closed_sets_of_size_seven_number_exactly_one', len(s7) == 1,
      '%d' % len(s7))
check('the_unique_order_7_subdesign_is_U', len(s7) == 1 and s7[0] == U,
      'U = %s' % sorted(s7[0]) if s7 else 'none')
check('no_subdesign_of_order_nine', not any(len(s) == 9 for s in closed))
check('no_subdesign_of_order_thirteen', not any(len(s) == 13 for s in closed))
s3 = set(s for s in closed if len(s) == 3)
check('the_closed_sets_of_size_three_are_exactly_the_117_blocks',
      s3 == set(frozenset(b) for b in blocks) and len(s3) == 117,
      '%d' % len(s3))

fano_through = {x: sum(1 for s in s7 if x in s) for x in range(V)}
check('fano_subdesigns_through_a_point_is_1_on_U_and_0_off_U',
      all(fano_through[x] == (1 if x in U else 0) for x in range(V)),
      'values %s' % sorted(set(fano_through.values())))
check('the_point_invariant_is_not_constant_so_Aut_is_not_point_transitive',
      len(set(fano_through.values())) > 1,
      '7 points with 1, 20 points with 0')
check('lemma_2w_plus_1_excludes_a_proper_subdesign_of_order_7_9_or_13_inside_one',
      all(u < 2 * w + 1
          for u in (7, 9, 13) for w in ADMISSIBLE if 3 < w < u),
      'so every such subdesign is the closure of three of its own points')
print('')

# ============================================================================
print('=== Step 8: admissibility arithmetic -- of these four the paper states '
      'only that 3 does not divide 7 (Section 1, on the vacuity of a sub-KTS(7)) '
      'and the inequality v >= 2u+1 (Lemma 2); the two congruences correspond to '
      'nothing in the paper and are recorded here as background')
# ============================================================================
check('27_is_admissible_for_a_Kirkman_triple_system', V % 6 == 3,
      '27 mod 6 = %d' % (V % 6))
check('7_is_admissible_for_a_Steiner_triple_system', 7 % 6 == 1,
      '7 mod 6 = %d' % (7 % 6))
check('the_necessary_condition_v_at_least_2u_plus_1_holds', V >= 2 * 7 + 1,
      '27 >= 15')
check('a_sub_KTS_7_is_vacuous_since_3_does_not_divide_7', 7 % 3 != 0,
      '7 points admit no partition into triples')
print('')

note('SCOPE: this program certifies exactly ONE object -- the KTS(27) of '
     'Table 1 of the paper, copied by hand into TABLE above -- and the '
     'following quantities, which were inventoried by hand from the paper when '
     'this program was written: the resolution, the 351 pairs, the seven inner '
     'blocks and their Fano structure, the 7+70+40 census, sigma and its '
     'induced class permutation, the reconstruction from five base classes, '
     'and the subdesign census behind Proposition 3. NOT ESTABLISHED HERE: no '
     'check below opens paper.tex, so two things are asserted rather than '
     'verified -- that the table above is the paper\'s Table 1, and that the '
     'list just given is COMPLETE for what the paper states about the object; '
     'a quantity the paper states and that list omits would go unnoticed. NOT '
     'COVERED: (a) any COUNT of KTS(27) with an '
     'STS(7) subdesign -- one object is exhibited and no enumeration or '
     'isomorph rejection is performed, so nothing here bounds how many exist '
     'or claims minimality or canonicity; (b) the ORDER of Aut(D) -- only that '
     'it contains sigma and, by the invariant above, that it is not transitive '
     'on points; (c) the orders v = 33 and 39 of the quoted sentence, and Open '
     'Problem 4 itself, which are untouched; (d) the BIBLIOGRAPHIC claims -- '
     'the byte locator of the quoted passage in arXiv:2110.07874v1, the '
     'sentence quoted from Colbourn-Magliveras-Mathon, and the question of '
     'priority, which the paper\'s Scope section expressly declines to claim -- '
     'none of which this program fetches or re-reads; (e) a second '
     'witness recorded in the underlying run but NOT printed in this paper, '
     'which no check here touches. NOT RE-RUN: the constraint search that '
     'produced the design -- this program checks the object, not the search, '
     'and the paper needs only the object.')

if _n_fail:
    print('VERDICT: %d CHECK(S) FAILED of %d' % (_n_fail, _n_pass + _n_fail))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _n_pass)
sys.exit(0)

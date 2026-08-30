#!/usr/bin/env python3
"""verify.py -- re-derivation of every computational claim of the accompanying note

    "An explicit 32-point set with a free half-turn and no convex 7-gon"

The program READS THE OBJECT AS PRINTED IN THE PAPER: the sixteen half-turn orbit
representatives of Section 2 (list REPS below, transcribed from the first display of
Section 2) and, independently, the thirty-two expanded points of the second display
(list PRINTED_32).  Nothing is loaded from a file; the paper is the input.

Python 3.9+, standard library only.  Every decision is made on exact integers: the only
divisions in the program are the exact ones inside `Fraction`-free integer comparisons
(angular order is decided by cross-multiplication, never by a quotient), and no floating
point value is ever compared, rounded or summed.

Output convention: one `PASS <name> [detail]` line per check, and a closing

    VERDICT: ALL <n> CHECKS PASS

with exit status 0 if and only if every check passed.  The scope statement the note
repeats under `## Scope` is printed verbatim, as a `NOT RE-RUN:` line, before the verdict.
"""

import sys
from itertools import combinations
from functools import cmp_to_key

# ---------------------------------------------------------------------------
# THE OBJECT, TRANSCRIBED FROM THE PAPER
# ---------------------------------------------------------------------------
# Section 2, first display: the sixteen orbit representatives P.
REPS = [(18, -4), (6, -18), (-17, 4), (-3, -18),
        (-5, 20), (12, 21), (15, 23), (-13, 5),
        (20, -25), (-22, 2), (23, 27), (-22, -26),
        (23, -2), (-30, -31), (27, -1), (4, 19)]

# Section 2, second display: the thirty-two points of W* = P u (-P), in the printed order.
PRINTED_32 = [
    (18, -4), (-18, 4), (6, -18), (-6, 18), (-17, 4), (17, -4), (-3, -18), (3, 18),
    (-5, 20), (5, -20), (12, 21), (-12, -21), (15, 23), (-15, -23), (-13, 5), (13, -5),
    (20, -25), (-20, 25), (-22, 2), (22, -2), (23, 27), (-23, -27), (-22, -26), (22, 26),
    (23, -2), (-23, 2), (-30, -31), (30, 31), (27, -1), (-27, 1), (4, 19), (-4, -19),
]

# Section 2, the hexagon exhibited after the table of convex-position counts.
PRINTED_HEXAGON = [(18, -4), (-18, 4), (6, -18), (-6, 18), (-3, -18), (3, 18)]

# Every number the paper asserts about W*, in one place, so that a reader can see what is
# being compared against what.  Keys are quoted in the PASS detail fields below.
CLAIMED = {
    'n': 32,
    'max_abs_coordinate': 31,
    'orbits': 16,
    'C(32,3)': 4960,
    'collinear_triples': 0,
    'min_abs_determinant': 1,
    'max_abs_determinant': 2892,
    'convex_counts': {3: 4960, 4: 24322, 5: 59002, 6: 56261},
    'max_convex_position_size': 6,
    'C(32,7)': 3365856,
    'convex_7_subsets': 0,
    'layer_signature': [6, 6, 6, 6, 6, 2],
    'admissible_rotation_orders': [1, 2, 4],
}

# ---------------------------------------------------------------------------
# CHECK PLUMBING
# ---------------------------------------------------------------------------
_n_pass = 0
_n_fail = 0


def ok(name, cond, detail=''):
    """Record one check.  A check is a check whatever its outcome, so a failure prints FAIL
    on the same counter and the verdict is withheld by a nonzero exit status."""
    global _n_pass, _n_fail
    if cond:
        _n_pass += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _n_fail += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE %s' % text)


def head(text):
    print()
    print('=== %s' % text)


# ---------------------------------------------------------------------------
# EXACT PRIMITIVES
# ---------------------------------------------------------------------------
def cross(o, a, b):
    """Twice the signed area of (o, a, b).  Integer in, integer out."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def in_convex_position(pts):
    """Caratheodory criterion, exactly as the source paper states it (main.tex line 251): a
    finite set in general position is in convex position iff no member lies strictly inside
    the triangle spanned by three others.  Pure integer arithmetic."""
    m = len(pts)
    for q in range(m):
        others = [t for t in range(m) if t != q]
        Q = pts[q]
        for a, b, c in combinations(others, 3):
            d1 = cross(pts[a], pts[b], Q)
            d2 = cross(pts[b], pts[c], Q)
            d3 = cross(pts[c], pts[a], Q)
            if (d1 > 0 and d2 > 0 and d3 > 0) or (d1 < 0 and d2 < 0 and d3 < 0):
                return False
    return True


def convex_hull(pts):
    """Monotone chain, exact integers, strict turns only (so with a set in general position
    the result is exactly the set of hull vertices, counterclockwise)."""
    pts = sorted(set(pts))
    if len(pts) <= 2:
        return list(pts)
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def convex_layers(pts):
    """The onion decomposition: peel hull vertices until nothing is left."""
    rest = list(pts)
    layers = []
    while rest:
        h = convex_hull(rest)
        if not h:
            break
        layers.append(h)
        hs = set(h)
        rest = [p for p in rest if p not in hs]
    return layers


def convex_position_counts(pts):
    """EXHAUSTIVE enumeration of every subset in convex position, AT EVERY SIZE AT ONCE, by
    growing angularly ordered left-turning chains from each vertex in turn.

    This is a different code path from `in_convex_position`: it never asks the Caratheodory
    question and never enumerates k-subsets for any fixed k.  Each convex polygon is
    generated exactly once, from its lexicographically smallest vertex v0, with the
    remaining vertices in increasing angular order about v0 -- all candidates are
    lexicographically greater than v0 and therefore lie in the open right half plane
    together with the upward vertical ray, so those angles lie in (-90, 90] and the turn at
    v0 itself is automatically left.  Because the enumeration is not truncated at any size,
    the largest size it reports IS the maximum size of a subset in convex position; nothing
    above that size needs a separate pass.

    Angular order is decided by integer cross-multiplication (dy_i * dx_j vs dy_j * dx_i
    with both dx positive), so no quotient and no float is formed.

    Returns (counts_by_size, max_size, a witness of max size, chains_visited).
    """
    S = list(pts)
    n = len(S)
    order = sorted(range(n), key=lambda i: S[i])
    counts = {}
    state = {'max': 0, 'wit': None, 'chains': 0}

    for si in range(n):
        v0 = S[order[si]]
        cand = [S[order[j]] for j in range(si + 1, n)]

        def acmp(p, q):
            dxp, dyp = p[0] - v0[0], p[1] - v0[1]
            dxq, dyq = q[0] - v0[0], q[1] - v0[1]
            # lexicographic order guarantees dx > 0, or dx == 0 with dy > 0
            if dxp == 0 and dxq == 0:
                raise AssertionError('two candidates vertically aligned with v0: collinear')
            if dxp == 0:
                return 1        # +90 degrees is the largest angle
            if dxq == 0:
                return -1
            lhs = dyp * dxq
            rhs = dyq * dxp
            if lhs == rhs:
                raise AssertionError('angular tie at %r: %r and %r are collinear with it'
                                     % (v0, p, q))
            return -1 if lhs < rhs else 1

        cand.sort(key=cmp_to_key(acmp))
        m = len(cand)

        # iterative DFS: a stack of chains, each a list of indices into `cand`
        stack = [[]]
        while stack:
            chain = stack.pop()
            state['chains'] += 1
            k = len(chain)
            if k >= 2:
                a = cand[chain[-2]]
                b = cand[chain[-1]]
                if cross(a, b, v0) > 0:          # closing turn back to v0
                    sz = k + 1
                    counts[sz] = counts.get(sz, 0) + 1
                    if sz > state['max']:
                        state['max'] = sz
                        state['wit'] = [v0] + [cand[t] for t in chain]
            start = chain[-1] + 1 if chain else 0
            for nxt in range(start, m):
                if k == 0:
                    stack.append([nxt])
                else:
                    a = v0 if k == 1 else cand[chain[-2]]
                    b = cand[chain[-1]]
                    if cross(a, b, cand[nxt]) > 0:
                        stack.append(chain + [nxt])
    return counts, state['max'], state['wit'], state['chains']


def binom(n, k):
    if k < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r


# ---------------------------------------------------------------------------
print('verification of the note: an explicit 32-point set with a free half-turn')
print('and no convex 7-gon -- W* = P u (-P), n = 32, k = 7, s = 2')
print('python %s, exact integer arithmetic only' % sys.version.split()[0])

# ===========================================================================
head('Step 1: the object, as printed in the paper')
# ===========================================================================
ok('reps_are_16_integer_points',
   len(REPS) == 16 and all(isinstance(c, int) for p in REPS for c in p),
   '16 orbit representatives, all coordinates Python ints')

EXPANDED = []
for (x, y) in REPS:
    EXPANDED.append((x, y))
    EXPANDED.append((-x, -y))

ok('expansion_of_the_16_reps_reproduces_the_printed_32_points',
   EXPANDED == PRINTED_32,
   'p -> (p, -p) over the first display equals the second display, element by element')

W = list(PRINTED_32)
n = len(W)
Wset = set(W)

ok('the_printed_points_are_32_and_pairwise_distinct',
   n == CLAIMED['n'] and len(Wset) == CLAIMED['n'],
   'len(list) = %d, len(set) = %d, claimed %d' % (n, len(Wset), CLAIMED['n']))

ok('the_origin_is_not_a_point_of_the_set', (0, 0) not in Wset,
   'required for the half-turn to act freely')

maxc = max(max(abs(x), abs(y)) for (x, y) in W)
ok('max_absolute_coordinate_is_31', maxc == CLAIMED['max_abs_coordinate'],
   'max |coordinate| = %d, claimed %d (attained by (-30,-31))' % (maxc, CLAIMED['max_abs_coordinate']))

ok('the_hexagon_exhibited_in_section_2_is_a_subset_of_the_set',
   all(p in Wset for p in PRINTED_HEXAGON),
   'all 6 printed hexagon vertices occur among the 32')

# ===========================================================================
head('Step 2: the half-turn h(z) = -z acts freely, with 16 orbits of size 2')
# ===========================================================================
h = lambda p: (-p[0], -p[1])

ok('h_maps_the_set_onto_itself', set(map(h, W)) == Wset,
   'h(W*) = W* as sets')

ok('h_has_order_exactly_two',
   all(h(h(p)) == p for p in W) and any(h(p) != p for p in W),
   'h o h = id on W*, and h is not the identity on W*')

ok('h_has_no_fixed_point_in_the_set', all(h(p) != p for p in W),
   'the unique fixed point of h in the plane is the origin, which is absent')

orbits = set()
for p in W:
    orbits.add(min(p, h(p)))
ok('the_action_has_16_free_orbits_of_size_two',
   len(orbits) == CLAIMED['orbits'] and all(len({p, h(p)}) == 2 for p in W),
   '%d orbits, every one of size 2, 16 * 2 = 32' % len(orbits))

ok('the_orbit_representatives_printed_are_a_transversal',
   {min(p, h(p)) for p in REPS} == orbits and len(REPS) == 16,
   'the 16 printed representatives meet each orbit exactly once')

# h is linear, so it preserves orientation determinants exactly.  Checked, not asserted.
same = all(cross(h(W[a]), h(W[b]), h(W[c])) == cross(W[a], W[b], W[c])
           for a, b, c in combinations(range(n), 3))
ok('h_preserves_every_orientation_determinant',
   same,
   'all C(32,3) = %d determinants are unchanged by h, so h preserves collinearity and '
   'convex position' % binom(32, 3))

# ===========================================================================
head('Step 3: general position -- all C(32,3) orientation determinants')
# ===========================================================================
dets = [abs(cross(W[a], W[b], W[c])) for a, b, c in combinations(range(n), 3)]
ok('all_C_32_3_triples_were_examined', len(dets) == CLAIMED['C(32,3)'],
   'triples examined = %d = C(32,3) = %d' % (len(dets), binom(32, 3)))
ok('no_three_points_are_collinear',
   sum(1 for d in dets if d == 0) == CLAIMED['collinear_triples'],
   'collinear triples = %d, claimed %d' % (sum(1 for d in dets if d == 0),
                                           CLAIMED['collinear_triples']))
ok('smallest_absolute_determinant_is_1', min(dets) == CLAIMED['min_abs_determinant'],
   'min |det| = %d, claimed %d' % (min(dets), CLAIMED['min_abs_determinant']))
ok('largest_absolute_determinant_is_2892', max(dets) == CLAIMED['max_abs_determinant'],
   'max |det| = %d, claimed %d' % (max(dets), CLAIMED['max_abs_determinant']))
ok('every_determinant_fits_in_12_unsigned_bits', max(dets) < 4096,
   'max |det| = %d < 4096; %d signed bits suffice' % (max(dets), max(dets).bit_length() + 1))

# ===========================================================================
head('Step 4: convex position -- exhaustive chain enumeration over ALL sizes at once')
# ===========================================================================
counts, maxsize, witness, chains = convex_position_counts(W)
note('chains visited = %d' % chains)
for sz in sorted(counts):
    note('subsets in convex position of size %d = %d' % (sz, counts[sz]))

for sz in sorted(CLAIMED['convex_counts']):
    ok('convex_position_subsets_of_size_%d_number_%d' % (sz, CLAIMED['convex_counts'][sz]),
       counts.get(sz, 0) == CLAIMED['convex_counts'][sz],
       'got %d, claimed %d, out of C(32,%d) = %d'
       % (counts.get(sz, 0), CLAIMED['convex_counts'][sz], sz, binom(32, sz)))

ok('the_triangle_count_equals_C_32_3_as_general_position_forces',
   counts.get(3, 0) == binom(32, 3),
   'every triple of a set in general position is in convex position: %d = %d'
   % (counts.get(3, 0), binom(32, 3)))

ok('no_subset_of_size_7_is_in_convex_position', counts.get(7, 0) == 0,
   'the enumeration reports %d subsets of size 7 in convex position' % counts.get(7, 0))

ok('the_maximum_size_of_a_subset_in_convex_position_is_exactly_6',
   maxsize == CLAIMED['max_convex_position_size'],
   'max = %d, claimed %d; witness %r' % (maxsize, CLAIMED['max_convex_position_size'], witness))

ok('the_enumeration_is_untruncated_so_no_size_above_6_occurs',
   max(counts) == 6 and all(counts.get(sz, 0) == 0 for sz in range(7, 33)),
   'sizes present: %s; sizes 7..32 all absent, so no convex m-gon exists for any m >= 7'
   % sorted(counts))

ok('the_max_size_witness_is_in_convex_position_by_the_caratheodory_test',
   in_convex_position(witness) and len(witness) == 6,
   'the chain enumeration and the Caratheodory test agree on the extremal witness')

ok('the_hexagon_printed_in_section_2_is_in_convex_position',
   in_convex_position(PRINTED_HEXAGON),
   'the printed hexagon %r passes the Caratheodory test' % (PRINTED_HEXAGON,))

# ===========================================================================
head('Step 5: independent cross-check -- Caratheodory over every 4- and 7-subset')
# ===========================================================================
# 5a. The two code paths are made to agree on a size where BOTH are cheap.
c4 = sum(1 for s in combinations(W, 4) if in_convex_position(list(s)))
ok('caratheodory_and_chain_enumeration_agree_on_all_C_32_4_subsets',
   c4 == counts.get(4, 0) == CLAIMED['convex_counts'][4],
   'Caratheodory over all C(32,4) = %d subsets gives %d; the chain enumeration gives %d'
   % (binom(32, 4), c4, counts.get(4, 0)))

# 5b. The decisive claim, by the second path, at full size and with no pruning.
# For every triple, the bitmask of the other points strictly inside its (oriented) triangle.
INS = {}
for (a, b, c) in combinations(range(n), 3):
    A, B, C = W[a], W[b], W[c]
    if cross(A, B, C) < 0:
        B, C = C, B
    mask = 0
    for p in range(n):
        if p in (a, b, c):
            continue
        P = W[p]
        if cross(A, B, P) > 0 and cross(B, C, P) > 0 and cross(C, A, P) > 0:
            mask |= 1 << p
    INS[(a, b, c)] = mask

BIT = [1 << i for i in range(n)]
seen = 0
convex7 = 0
for sub in combinations(range(n), 7):
    seen += 1
    M = 0
    for i in sub:
        M |= BIT[i]
    good = True
    for tri in combinations(sub, 3):
        if INS[tri] & M:
            good = False
            break
    if good:
        convex7 += 1

ok('every_one_of_the_C_32_7_seven_subsets_was_examined',
   seen == CLAIMED['C(32,7)'] == binom(32, 7),
   '7-subsets examined = %d = C(32,7) = %d, no pruning' % (seen, binom(32, 7)))
ok('brute_force_finds_zero_seven_subsets_in_convex_position',
   convex7 == CLAIMED['convex_7_subsets'],
   'independent Caratheodory path over all %d subsets: %d in convex position, claimed %d'
   % (seen, convex7, CLAIMED['convex_7_subsets']))

# 5c. Heredity, which is why the maximum settles every larger size at once.
hered = True
for s in combinations(PRINTED_HEXAGON, 5):
    if not in_convex_position(list(s)):
        hered = False
ok('convex_position_is_inherited_by_subsets_on_the_exhibited_hexagon',
   hered,
   'all C(6,5) = 6 five-subsets of the printed hexagon are in convex position; heredity is '
   'why "no 7 in convex position" already implies "no m in convex position" for every m > 7')

# ===========================================================================
head('Step 6: controls, both polarities, on the criterion and on the enumeration')
# ===========================================================================
# POSITIVE control 1: points on a parabola -- EVERY subset is in convex position, so the
# enumeration must return the full binomial profile.  Forces the counter to be able to
# say yes, at every size, including 7.
PARA = [(i, i * i) for i in range(12)]
pcounts, pmax, _, _ = convex_position_counts(PARA)
ok('positive_control_parabola_enumeration_returns_the_full_binomial_profile',
   all(pcounts.get(k, 0) == binom(12, k) for k in range(3, 13)) and pmax == 12,
   '12 points on y = x^2: sizes 3..12 counted %s, expected %s'
   % ([pcounts.get(k, 0) for k in range(3, 13)], [binom(12, k) for k in range(3, 13)]))
ok('positive_control_parabola_caratheodory_accepts_all_C_12_7_seven_subsets',
   sum(1 for s in combinations(PARA, 7) if in_convex_position(list(s))) == binom(12, 7),
   'all C(12,7) = %d seven-subsets of the parabola read "in convex position"' % binom(12, 7))

# POSITIVE control 2: a rounded regular heptagon -- 7 points that must read YES.
HEPT = [(1000, 0), (623, 782), (-223, 975), (-901, 434),
        (-901, -434), (-223, -975), (623, -782)]
ok('positive_control_regular_heptagon_is_in_general_position',
   all(cross(*t) != 0 for t in combinations(HEPT, 3)),
   'all C(7,3) = %d determinants of the rounded regular heptagon are nonzero' % binom(7, 3))
ok('positive_control_regular_heptagon_reads_in_convex_position',
   in_convex_position(HEPT),
   'the criterion CAN return True on 7 points, so a False on W* is not vacuous')

# NEGATIVE control: a hexagon plus one interior point -- 7 points that must read NO.
HEX_PLUS = [(2, 0), (1, 2), (-1, 2), (-2, 0), (-1, -2), (1, -2), (0, 1)]
ok('negative_control_hexagon_plus_interior_point_is_in_general_position',
   all(cross(*t) != 0 for t in combinations(HEX_PLUS, 3)),
   'all C(7,3) = %d determinants nonzero, so the control is a legal input' % binom(7, 3))
ok('negative_control_hexagon_plus_interior_point_reads_not_in_convex_position',
   not in_convex_position(HEX_PLUS),
   'the criterion CAN return False on 7 points')
hcounts, hmax, _, _ = convex_position_counts(HEX_PLUS)
ok('negative_control_chain_enumeration_caps_the_hexagon_plus_point_at_6',
   hmax == 6 and hcounts.get(7, 0) == 0,
   'the SECOND code path also reports max = %d and 0 subsets of size 7' % hmax)

# ===========================================================================
head('Step 7: the onion decomposition and the lemma it satisfies')
# ===========================================================================
layers = convex_layers(W)
sig = [len(L) for L in layers]
note('convex layer sizes = %s' % sig)
ok('the_convex_layer_signature_is_6_6_6_6_6_2', sig == CLAIMED['layer_signature'],
   'got %s, claimed %s' % (sig, CLAIMED['layer_signature']))
ok('the_layer_sizes_sum_to_32', sum(sig) == 32, 'sum(%s) = %d' % (sig, sum(sig)))
ok('every_convex_layer_is_carried_to_itself_by_h',
   all(set(map(h, L)) == set(L) for L in layers),
   'negation is an affine bijection, so it permutes the onion layers; here it fixes each')
ok('every_convex_layer_has_even_size', all(len(L) % 2 == 0 for L in layers),
   'each layer is h-invariant and the origin is absent, so h acts freely on it: sizes %s' % sig)
ok('every_convex_layer_has_size_at_most_6', max(sig) <= 6,
   'a layer is in convex position, so its size is at most k - 1 = 6; max = %d' % max(sig))
ok('the_number_of_layers_is_at_least_ceil_32_over_6',
   len(sig) >= -(-32 // 6),
   '%d layers, and the lemma forces at least ceil(32/6) = %d' % (len(sig), -(-32 // 6)))
ok('every_layer_is_in_convex_position_by_the_caratheodory_test',
   all(in_convex_position(L) for L in layers),
   'all %d layers pass the independent test' % len(layers))

# ===========================================================================
head('Step 8: which rotation orders can act at all on 32 points without a 7-gon')
# ===========================================================================
# A rotation of order s permutes the set with all orbits of size s except at most one fixed
# point (the centre), so 32 = 0 or 1 (mod s).  And an orbit of size s >= 7 is the vertex set
# of a regular s-gon, hence contains 7 points in convex position.
by_congruence = [s for s in range(1, 33) if 32 % s in (0, 1)]
ok('the_congruence_32_mod_s_in_0_1_admits_exactly_1_2_4_8_16_31_32',
   by_congruence == [1, 2, 4, 8, 16, 31, 32],
   'orders passing 32 = 0 or 1 (mod s): %s' % by_congruence)
ok('orders_3_5_and_6_are_excluded_by_the_congruence_alone',
   all(32 % s not in (0, 1) for s in (3, 5, 6)),
   '32 mod 3 = %d, 32 mod 5 = %d, 32 mod 6 = %d, and none is 0 or 1'
   % (32 % 3, 32 % 5, 32 % 6))
survivors = [s for s in by_congruence if s < 7]
ok('the_admissible_rotation_orders_are_exactly_1_2_and_4',
   survivors == CLAIMED['admissible_rotation_orders'],
   'after deleting every s >= 7 (an orbit of size s >= 7 is a regular s-gon and so contains '
   '7 points in convex position): %s' % survivors)
# The exclusion of s >= 7 is itself checked on a concrete regular 7-gon and 8-gon.
OCT = [(1000, 0), (707, 707), (0, 1000), (-707, 707),
       (-1000, 0), (-707, -707), (0, -1000), (707, -707)]
ok('a_regular_8_gon_orbit_does_contain_7_points_in_convex_position',
   any(in_convex_position(list(s)) for s in combinations(OCT, 7)),
   'the s >= 7 exclusion is exhibited, not merely asserted: a rounded regular octagon has a '
   '7-subset in convex position')

# ===========================================================================
head('Scope')
# ===========================================================================
print('NOT RE-RUN: (1) the SEARCH that produced W* is not re-run here and is not '
      'reproducible from a seed -- each restart of the original grid descent was cut by a '
      '90-second wall-clock budget, so W* is pinned as the integer list printed in the '
      'paper and not as a random seed; nothing above depends on how it was found. '
      '(2) The claim that 31 is MINIMAL for the largest absolute coordinate of such a set '
      'is NOT made and NOT tested: 31 is where a descending ladder of grid bounds stopped '
      'on budget, and no exhaustive search of any smaller grid was performed. '
      '(3) NO statement about g(7) is re-derived, because none is made: g(7) is unknown, '
      'g(7) >= 33 is classical and independent of this object, and the 4-fold negative '
      'result at n = 32 is quoted from the source paper and NOT recomputed here. '
      '(4) The source paper\'s own published 16-point 4-fold configuration is NOT re-counted '
      'here (its coordinates are not reproduced in this note); the calibration of the '
      'convex-position criterion against that published object was done elsewhere, and the '
      'criterion is instead pinned here by the four controls of Step 6. '
      '(5) The nine further 32-point witnesses at larger coordinate bounds mentioned in the '
      'note are NOT verified here; only W* is. '
      '(6) The REFLECTION (mirror) symmetry of order 2 at n = 32 is untouched: this program '
      'says nothing about it. '
      '(7) In Step 8 the exclusion of every rotation order s >= 7 is checked on one concrete '
      'rounded regular octagon, not on all such s; the general argument (an orbit of size s '
      'is a regular s-gon, and any 7 of its vertices are in convex position) is a hand proof '
      'and is given in the note.')

print()
if _n_fail:
    print('VERDICT: %d CHECKS FAILED of %d' % (_n_fail, _n_fail + _n_pass))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _n_pass)
sys.exit(0)

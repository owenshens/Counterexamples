#!/usr/bin/env python3
"""verify.py -- re-derivation of every computational claim in

    "An order-49 counterexample to the heavy-chopped-rectangle interval
     conjecture of Gottlieb, Krnc and Mursic"

Python 3.9+, STANDARD LIBRARY ONLY: no numpy, no sympy, no networkx, no
external data file.  All arithmetic is exact integer arithmetic; no float
value decides anything.

The program READS THE OBJECTS PRINTED IN THE PAPER -- the partition
lambda = [11,11,11,8,8], its cell (a,b) = (10,4), the 45-option value set,
the second minimal witness, and both of the source paper's printed
appendix tables -- and re-derives every quantity the paper asserts about
them from the move rule alone.

Two independent Sprague-Grundy engines are implemented:

  ENGINE A  works on the parts tuple.  A row move keeps a proper
            sub-multiset of the parts (enumerated by run multiplicity); a
            column move is the same operation carried out on the
            conjugate and conjugated back.

  ENGINE B  works on the Young diagram as a set of cells.  A row move
            keeps a proper subset of the ROW INDEX SET and re-sorts; a
            column move keeps a proper subset of the COLUMN INDEX SET and
            recounts each row.  It never calls a conjugation routine and
            it never collapses equal rows, so it is an independent check
            both of the conjugate implementation and of the
            run-collapsing shortcut Engine A takes.

Output: one `PASS <name> [detail]` line per check, then

    VERDICT: ALL <n> CHECKS PASS

and exit status 0 if and only if every check passed.
"""

import itertools
import sys

# ---------------------------------------------------------------------------
# 0.  THE OBJECTS, EXACTLY AS THE PAPER PRINTS THEM
# ---------------------------------------------------------------------------

LAMBDA = (11, 11, 11, 8, 8)          # the counterexample
CELL = (10, 4)                       # (a, b)
LAMBDA_ORDER = 49
LAMBDA_SG = 10                       # the paper's computed value
LAMBDA_REQUIRED = 15                 # a + b + 1, what the conjecture demands
LAMBDA_CONJUGATE = (5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3)
OPTION_COUNT = 45
OPTION_VALUE_SET = frozenset([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14])

# The second partition of order 49 that the paper exhibits, in a different cell.
SECOND = (9, 9, 7, 7, 6, 5, 4, 2)
SECOND_CELL = (8, 7)
SECOND_SG = 3
SECOND_REQUIRED = 16

# The two non-heavy members of cell (11,4) the paper lists.
CELL_11_4_NONHEAVY = ((12, 12, 12, 11, 10), (12, 12, 12, 11, 11))

# arXiv:2506.04991v1, Appendix A ("Small partitions with sg value 2"): all
# partitions of order at most 26 with sg = 2, up to conjugation, transcribed
# from the \youngdiagram{...} arguments in the source of arXiv1.tex.
APPENDIX_A_SG2 = (
    (1, 1),
    (3, 3, 3, 1),
    (3, 3, 3, 3),
    (5, 4, 4, 3, 1),
    (5, 5, 5, 3, 3, 1),
    (5, 5, 5, 4, 3, 1),
    (5, 5, 5, 5, 3, 1),
    (5, 5, 5, 4, 4, 1),
    (5, 5, 5, 3, 3, 3),
    (5, 5, 5, 5, 4, 1),
    (5, 5, 5, 5, 5, 1),
    (5, 5, 5, 5, 3, 3),
    (5, 5, 5, 4, 4, 3),
)

# arXiv:2506.04991v1, Appendix B ("Small heavy partitions under 1-PNim"): all
# heavy partitions of order at most 8, up to conjugation, transcribed the same
# way.
APPENDIX_B_HEAVY = (
    (1,), (1, 1), (2, 1), (1, 1, 1), (2, 1, 1), (1, 1, 1, 1),
    (3, 1, 1), (2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1), (3, 2, 1), (3, 1, 1, 1),
    (2, 2, 2), (2, 2, 1, 1), (2, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1), (4, 1, 1, 1), (3, 2, 2),
    (3, 2, 1, 1), (3, 1, 1, 1, 1), (2, 2, 2, 1), (2, 2, 1, 1, 1), (2, 1, 1, 1, 1, 1),
    (1, 1, 1, 1, 1, 1, 1),
    (4, 2, 1, 1), (4, 1, 1, 1, 1), (3, 3, 1, 1), (3, 2, 2, 1), (3, 2, 1, 1, 1),
    (3, 1, 1, 1, 1, 1),
    (2, 2, 2, 1, 1), (2, 2, 1, 1, 1, 1), (2, 1, 1, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1, 1, 1),
)

MIN_ORDER = 49          # the paper's minimality claim
MIN_ORDER_ATTAINED_BY = (LAMBDA, SECOND)

# ---------------------------------------------------------------------------
# 1.  THE CHECK HARNESS
# ---------------------------------------------------------------------------

_passed = 0
_failed = 0


def say(*a):
    print(*a, flush=True)


def check(name, cond, detail=''):
    global _passed, _failed
    if cond:
        _passed += 1
        say('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _failed += 1
        say('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(*a):
    say('NOTE', *a)


# ---------------------------------------------------------------------------
# 2.  ENGINE A -- parts tuple, run-collapsed moves, conjugate for columns
# ---------------------------------------------------------------------------

def norm(p):
    """A partition as a weakly decreasing tuple of positive parts."""
    return tuple(x for x in sorted(p, reverse=True) if x > 0)


def conjugate(p):
    """Exchange rows and columns of the Young diagram."""
    if not p:
        return ()
    return tuple(sum(1 for x in p if x >= j) for j in range(1, p[0] + 1))


def _runs(seq):
    out = []
    for x in seq:
        if out and out[-1][0] == x:
            out[-1][1] += 1
        else:
            out.append([x, 1])
    return out


def proper_sub_multisets(seq):
    """Every kept sub-multiset of a weakly decreasing `seq` except `seq` itself.

    The empty selection IS proper: it is the move to the empty partition.
    """
    rr = _runs(seq)
    full = tuple(c for _, c in rr)
    for choice in itertools.product(*[range(c + 1) for _, c in rr]):
        if choice == full:
            continue
        out = []
        for (v, _), k in zip(rr, choice):
            out.extend([v] * k)
        yield tuple(out)


def options_A(p):
    """The 1-PNim option set of `p`: row moves (1) and column moves (2) of the
    source Definition."""
    out = set()
    for kept in proper_sub_multisets(p):
        out.add(norm(kept))
    c = conjugate(p)
    for kept in proper_sub_multisets(c):
        out.add(norm(conjugate(norm(kept))))
    return out


MEMO_A = {(): 0}


def sg_A(p):
    """Sprague-Grundy value, iterative so no recursion limit is involved."""
    p = norm(p)
    stack = [p]
    while stack:
        cur = stack[-1]
        if cur in MEMO_A:
            stack.pop()
            continue
        opts = options_A(cur)
        pending = [o for o in opts if o not in MEMO_A]
        if pending:
            stack.extend(pending)
            continue
        seen = {MEMO_A[o] for o in opts}
        m = 0
        while m in seen:
            m += 1
        MEMO_A[cur] = m
        stack.pop()
    return MEMO_A[p]


def longest_play(p):
    """The source's Proposition: the longest 1-PNim play from a nonempty
    partition has length lambda_1 + r - 1."""
    return p[0] + len(p) - 1 if p else 0


def heavy(p):
    """The source's definition (line 484): heavy iff sg equals the length of
    the longest play."""
    return bool(p) and sg_A(p) == longest_play(p)


# ---------------------------------------------------------------------------
# 3.  ENGINE B -- Young diagram as a cell set, raw index subsets, no conjugate
# ---------------------------------------------------------------------------

def options_B(p):
    """Independent implementation.  A row move keeps a proper subset of the row
    INDICES; a column move keeps a proper subset of the column INDICES and
    recounts every row.  Nothing here collapses equal rows and nothing here
    calls `conjugate`.  Subsets are enumerated as bit masks over the raw index
    ranges, so every one of the 2^r - 1 and 2^c - 1 proper index subsets is
    visited separately even when rows or columns repeat."""
    p = norm(p)
    r = len(p)
    c = p[0] if p else 0
    out = set()
    for mask in range((1 << r) - 1):            # mask == all ones is improper
        out.add(norm([p[i] for i in range(r) if (mask >> i) & 1]))
    for mask in range((1 << c) - 1):
        # pref[x] = how many kept columns lie strictly left of column x, so row
        # i of length p[i] keeps pref[p[i]] cells.
        pref = [0] * (c + 1)
        for j in range(c):
            pref[j + 1] = pref[j] + ((mask >> j) & 1)
        out.add(norm([pref[p[i]] for i in range(r)]))
    return out


MEMO_B = {(): 0}


def sg_B(p):
    p = norm(p)
    stack = [p]
    while stack:
        cur = stack[-1]
        if cur in MEMO_B:
            stack.pop()
            continue
        opts = options_B(cur)
        pending = [o for o in opts if o not in MEMO_B]
        if pending:
            stack.extend(pending)
            continue
        seen = {MEMO_B[o] for o in opts}
        m = 0
        while m in seen:
            m += 1
        MEMO_B[cur] = m
        stack.pop()
    return MEMO_B[p]


# ---------------------------------------------------------------------------
# 4.  PARTITION UTILITIES
# ---------------------------------------------------------------------------

def partitions_of(n):
    def rec(n, mx):
        if n == 0:
            yield ()
            return
        for k in range(min(n, mx), 0, -1):
            for rest in rec(n - k, k):
                yield (k,) + rest
    yield from rec(n, n)


def partitions_upto(n):
    return [p for m in range(1, n + 1) for p in partitions_of(m)]


def young_le(lam, mu):
    """The source's Young's lattice order: lam <= mu iff ell(lam) <= ell(mu)
    and lam_j <= mu_j for every j = 1..ell(lam)."""
    if len(lam) > len(mu):
        return False
    return all(lam[j] <= mu[j] for j in range(len(lam)))


def conj_class(p):
    return frozenset((norm(p), conjugate(norm(p))))


def binom(n, k):
    if k < 0 or k > n:
        return 0
    num = 1
    for i in range(k):
        num = num * (n - i) // (i + 1)
    return num


def catalan(n):
    return binom(2 * n, n) // (n + 1)


def staircase(c, r):
    """str_{c,r} = [c, c-1, ..., c-r+1] of the source, for c >= r >= 1."""
    return tuple(range(c, c - r, -1))


def interval_members(a, b):
    """Every partition lambda with [a+1,a,...,a-b+1] <= lambda <= [(a+1)^(b+1)],
    built by the parameterisation: exactly b+1 rows, row j between the
    staircase floor and a+1, weakly decreasing."""
    lo = staircase(a + 1, b + 1)
    res = []

    def rec(i, prev, acc):
        if i == b + 1:
            res.append(tuple(acc))
            return
        for x in range(lo[i], min(a + 1, prev) + 1):
            rec(i + 1, x, acc + [x])
    rec(0, a + 1, [])
    return res


def interval_members_by_filter(a, b):
    """The same set, built the slow honest way: enumerate EVERY partition in
    the (b+1) x (a+1) Young box and keep those the coded Young order places
    between the two endpoints."""
    lo = staircase(a + 1, b + 1)
    hi = tuple([a + 1] * (b + 1))
    box = []

    def rec(i, prev, acc):
        if i == b + 1:
            box.append(norm(acc))
            return
        for x in range(0, prev + 1):
            rec(i + 1, x, acc + [x])
    rec(0, a + 1, [])
    box = set(box)
    return sorted(p for p in box if young_le(lo, p) and young_le(p, hi))


def lower_endpoint_order(a, b):
    """(b+1)(2a+2-b)/2 -- the order of the staircase [a+1,a,...,a-b+1], which
    is the least order of any member of cell (a,b)."""
    return (b + 1) * (2 * a + 2 - b) // 2


def admissible(a, b):
    """The conjecture's hypothesis: the rectangle [(a+1)^(b+1)] is heavy.  By
    the source's obs:lucas that is binom(a+b,b) odd, which Kummer's theorem
    makes a & b == 0.  The commented clause a >= b is carried separately."""
    return (a & b) == 0


# ===========================================================================
say('verification of: an order-49 counterexample to the heavy-chopped-rectangle')
say('interval conjecture of Gottlieb, Krnc and Mursic (arXiv:2506.04991v1,')
say('conj:heavychoppedrect).  Exact integer arithmetic; standard library only.')
say('python %s' % sys.version.split()[0])
say('')

a, b = CELL

# ---------------------------------------------------------------------------
say('=== Step 1: the exhibited object, read as the paper prints it')
# ---------------------------------------------------------------------------
note('lambda =', list(LAMBDA), ' cell (a,b) =', CELL)
check('lambda_is_a_partition_weakly_decreasing_and_positive',
      LAMBDA == norm(LAMBDA) and all(x > 0 for x in LAMBDA),
      'norm(lambda) = %s' % (list(norm(LAMBDA)),))
check('lambda_order_is_49', sum(LAMBDA) == LAMBDA_ORDER == 49,
      'sum = %d' % sum(LAMBDA))
check('lambda_first_part_is_a_plus_one', LAMBDA[0] == a + 1,
      '%d = %d + 1' % (LAMBDA[0], a))
check('lambda_length_is_b_plus_one', len(LAMBDA) == b + 1,
      '%d = %d + 1' % (len(LAMBDA), b))
check('longest_play_from_lambda_is_15', longest_play(LAMBDA) == 15,
      'lambda_1 + r - 1 = 11 + 5 - 1 = %d' % longest_play(LAMBDA))
check('the_value_the_conjecture_demands_is_a_plus_b_plus_one_equals_the_longest_play',
      a + b + 1 == longest_play(LAMBDA) == LAMBDA_REQUIRED,
      'a+b+1 = %d' % (a + b + 1))
note('conjugate(lambda) =', list(conjugate(LAMBDA)))
check('conjugate_of_lambda_is_5_to_the_8_then_3_to_the_3',
      conjugate(LAMBDA) == LAMBDA_CONJUGATE,
      '[5^8,3^3] = %s' % (list(LAMBDA_CONJUGATE),))
check('conjugation_is_an_involution_and_preserves_order',
      all(conjugate(conjugate(p)) == p and sum(conjugate(p)) == sum(p)
          for p in partitions_upto(9)),
      'checked on all %d nonempty partitions of order <= 9' % len(partitions_upto(9)))
say('')

# ---------------------------------------------------------------------------
say('=== Step 2: two independent engines for the same move rule')
# ---------------------------------------------------------------------------
oa = options_A(LAMBDA)
ob = options_B(LAMBDA)
check('the_two_engines_produce_the_same_option_set_for_lambda', oa == ob,
      '%d options each, symmetric difference %d' % (len(oa), len(oa ^ ob)))

# the reachable set of lambda: the closure of lambda under the move relation
reach = set()
frontier = [norm(LAMBDA)]
while frontier:
    cur = frontier.pop()
    if cur in reach:
        continue
    reach.add(cur)
    frontier.extend(options_A(cur))
note('reachable set of lambda under the move relation:', len(reach), 'states')
check('the_two_engines_produce_the_same_option_set_at_every_reachable_state',
      all(options_A(p) == options_B(p) for p in reach),
      '%d states' % len(reach))
check('the_two_engines_agree_on_sg_at_every_reachable_state',
      all(sg_A(p) == sg_B(p) for p in reach),
      '%d states, 0 disagreements' % len(reach))
small = partitions_upto(9)
check('the_two_engines_agree_on_sg_for_every_partition_of_order_at_most_9',
      all(sg_A(p) == sg_B(p) for p in small),
      '%d partitions' % len(small))
check('sg_is_conjugation_invariant_as_obs_conjinvariance_states',
      all(sg_A(p) == sg_A(conjugate(p)) for p in partitions_upto(14)),
      'all %d nonempty partitions of order <= 14' % len(partitions_upto(14)))
say('')

# ---------------------------------------------------------------------------
say('=== Step 3: the hypothesis of the conjecture holds at (a,b) = (10,4)')
# ---------------------------------------------------------------------------
rect = tuple([a + 1] * (b + 1))
note('rectangle [(a+1)^(b+1)] = [11^5]; binom(a+b,b) = binom(14,4) =', binom(14, 4))
check('a_and_b_are_disjoint_in_binary', (a & b) == 0, '10 & 4 = %d' % (a & b))
check('a_is_at_least_b_so_the_commented_clause_also_holds', a >= b,
      '%d >= %d' % (a, b))
check('binom_14_4_is_1001_and_odd', binom(14, 4) == 1001 and binom(14, 4) % 2 == 1,
      'binom(14,4) = %d' % binom(14, 4))
check('rectangle_11_5_has_sg_15_by_engine_A', sg_A(rect) == 15, 'sg = %d' % sg_A(rect))
check('rectangle_11_5_has_sg_15_by_engine_B', sg_B(rect) == 15, 'sg = %d' % sg_B(rect))
check('rectangle_11_5_sg_matches_the_closed_form_of_thm_rect',
      sg_A(rect) == ((len(rect) - 1) ^ (rect[0] - 1)) + 1,
      '((5-1) XOR (11-1)) + 1 = %d' % (((5 - 1) ^ (11 - 1)) + 1))
check('rectangle_11_5_is_heavy_so_the_hypothesis_is_satisfied',
      heavy(rect) and sg_A(rect) == longest_play(rect) == 15,
      'sg = L = %d' % longest_play(rect))
bad_lucas = [(x, y) for x in range(65) for y in range(65)
             if (binom(x + y, y) % 2 == 1) != ((x & y) == 0)]
check('obs_lucas_hypothesis_is_exactly_a_and_b_equals_zero',
      not bad_lucas, 'binom(x+y,y) odd <=> x & y == 0 on all %d pairs x,y <= 64'
      % (65 * 65))
bad_rect = [(r, c) for r in range(1, 10) for c in range(1, 10)
            if sg_A(tuple([c] * r)) != ((r - 1) ^ (c - 1)) + 1]
check('thm_rect_closed_form_holds_on_all_81_rectangles_r_c_at_most_9',
      not bad_rect, '81 rectangles, %d mismatches' % len(bad_rect))
bad_obs = [(r, c) for r in range(1, 10) for c in range(1, 10)
           if heavy(tuple([c] * r)) != (binom(c + r - 2, r - 1) % 2 == 1)]
check('obs_lucas_agrees_with_the_engine_on_all_81_rectangles',
      not bad_obs, '81 rectangles, %d mismatches' % len(bad_obs))
say('')

# ---------------------------------------------------------------------------
say('=== Step 4: lambda lies strictly inside the conjectured interval')
# ---------------------------------------------------------------------------
lo = staircase(a + 1, b + 1)
hi = rect
note('lower endpoint [a+1,a,...,a-b+1] =', list(lo), ' upper endpoint =', list(hi))
check('lower_endpoint_is_the_staircase_str_11_5_of_the_source',
      lo == (11, 10, 9, 8, 7), 'str_{11,5} = %s' % (list(lo),))
check('lower_endpoint_is_heavy_by_prop_padded_stair', heavy(lo),
      'sg = L = %d' % sg_A(lo))
check('lower_endpoint_is_at_most_lambda_in_youngs_lattice', young_le(lo, LAMBDA),
      '11<=11, 10<=11, 9<=11, 8<=8, 7<=8 on 5 rows')
check('lambda_is_at_most_the_upper_endpoint_in_youngs_lattice',
      young_le(LAMBDA, hi), 'every part <= 11 and 5 <= 5 rows')
check('lambda_is_neither_endpoint_so_it_is_an_interior_member',
      LAMBDA != lo and LAMBDA != hi)
mem_param = sorted(interval_members(a, b))
mem_filter = sorted(interval_members_by_filter(a, b))
check('the_parameterised_member_set_equals_the_set_got_by_filtering_the_whole_5x11_box',
      mem_param == mem_filter,
      '%d members both ways, box size binom(16,5) = %d'
      % (len(mem_param), binom(16, 5)))
check('cell_10_4_has_catalan_5_equals_42_members',
      len(mem_param) == catalan(b + 1) == 42, '%d members' % len(mem_param))
check('lambda_is_one_of_them', norm(LAMBDA) in set(mem_param))
say('')

# ---------------------------------------------------------------------------
say('=== Step 5: the option set of lambda, exactly as the paper describes it')
# ---------------------------------------------------------------------------
row_opts = set()
for i in range(4):
    for j in range(3):
        if (i, j) == (3, 2):
            continue
        row_opts.add(norm([11] * i + [8] * j))
col_opts = set()
for p in range(9):
    for q in range(4):
        if (p, q) == (8, 3):
            continue
        col_opts.add(norm([p + q, p + q, p + q, p, p]))
check('the_paper_s_row_option_family_11_to_the_i_8_to_the_j_has_11_members',
      len(row_opts) == 11, '0<=i<=3, 0<=j<=2, (i,j) != (3,2)')
check('the_paper_s_column_option_family_p_q_p_q_p_q_p_p_has_35_members',
      len(col_opts) == 35, '0<=p<=8, 0<=q<=3, (p,q) != (8,3)')
check('the_two_families_overlap_in_exactly_one_place_the_empty_partition',
      row_opts & col_opts == {()},
      'overlap = %s' % sorted(map(list, row_opts & col_opts)))
check('their_union_is_exactly_the_engine_s_option_set_of_45',
      (row_opts | col_opts) == oa and len(oa) == OPTION_COUNT,
      '%d = 11 + 35 - 1 options' % len(oa))
check('lambda_itself_is_the_excluded_column_tuple_p_q_equals_8_3',
      norm([8 + 3, 8 + 3, 8 + 3, 8, 8]) == norm(LAMBDA),
      'so the column family is closed and lambda is not its own option')
say('')

# ---------------------------------------------------------------------------
say('=== Step 6: the mex table, and the refutation')
# ---------------------------------------------------------------------------
opts = sorted(oa, key=lambda p: (-sum(p), p))
vals = {}
for o in opts:
    vals[o] = sg_A(o)
    say('NOTE OPT %-24s order %3d sg %2d' % (str(list(o)), sum(o), vals[o]))
check('prop_bound_holds_at_every_reachable_state',
      all(sg_A(p) <= longest_play(p) for p in reach if p),
      'sg <= lambda_1 + r - 1 on all %d nonempty reachable states' % len(reach - {()}))
check('every_option_loses_at_least_one_row_or_column_the_authors_f_inequality',
      all(longest_play(o) <= longest_play(LAMBDA) - 1 for o in opts if o),
      'mu_1 + ell(mu) <= lambda_1 + ell(lambda) - 1 for all 44 nonempty options')
check('hence_no_option_can_carry_the_value_15',
      all(vals[o] <= 14 for o in opts),
      'max option value = %d' % max(vals.values()))
vset = set(vals.values())
check('the_45_option_values_are_exactly_0_to_9_and_11_to_14',
      vset == OPTION_VALUE_SET,
      'value set = %s' % sorted(vset))
check('the_value_10_is_realised_by_no_option', 10 not in vset,
      '10 is the unique missing value below 15')
mex = 0
while mex in vset:
    mex += 1
check('the_mex_of_the_option_values_is_10', mex == 10, 'mex = %d' % mex)
check('sg_of_lambda_is_10_by_engine_A', sg_A(LAMBDA) == LAMBDA_SG == 10,
      'sg = %d' % sg_A(LAMBDA))
check('sg_of_lambda_is_10_by_engine_B_independently',
      sg_B(LAMBDA) == LAMBDA_SG == 10, 'sg = %d' % sg_B(LAMBDA))
check('lambda_is_not_heavy_its_sg_falls_5_short_of_the_longest_play',
      not heavy(LAMBDA) and longest_play(LAMBDA) - sg_A(LAMBDA) == 5,
      'sg = 10, L = 15, deficit 5')
check('REFUTATION_a_member_of_a_heavy_rectangle_interval_is_not_heavy',
      (admissible(a, b) and a >= b and heavy(rect)
       and young_le(lo, LAMBDA) and young_le(LAMBDA, hi)
       and not heavy(LAMBDA)),
      'hypothesis holds at (10,4), lambda = [11,11,11,8,8] is in the interval, '
      'sg(lambda) = 10 != 15 = a+b+1')
say('')

# ---------------------------------------------------------------------------
say('=== Step 7: how much of the certificate needs no solver')
# ---------------------------------------------------------------------------
# The options that are rectangles: their value is forced by thm:rect, a proved
# theorem of the source; the empty partition has value 0 by definition.
forced = {}
for o in opts:
    if not o:
        forced[o] = 0
    elif len(set(o)) == 1:
        forced[o] = ((len(o) - 1) ^ (o[0] - 1)) + 1
note('options whose value a published closed form already forces:', len(forced))
check('the_forced_values_agree_with_the_engine_on_every_such_option',
      all(forced[o] == vals[o] for o in forced),
      '%d options: %d rectangles by thm:rect plus the empty partition'
      % (len(forced), len(forced) - 1))
check('16_of_the_45_options_are_rectangles_and_one_more_is_the_empty_partition',
      len(forced) == 17 and sum(1 for o in forced if o) == 16,
      '17 options carry a value no solver is needed for')
check('the_forced_values_already_cover_0_through_9_so_sg_at_least_10_needs_no_solver',
      set(range(10)) <= set(forced.values()),
      'forced value set = %s' % sorted(set(forced.values())))
# prop:2rows of the source covers the three two-row options.


def prop_2rows(c1, c2):
    return c1 - 1 if (c1 == c2 and c1 % 2 == 0) else c1 + 1


two_row = [o for o in opts if len(o) == 2]
check('prop_2rows_of_the_source_forces_the_three_two_row_options',
      len(two_row) == 3 and all(prop_2rows(o[0], o[1]) == vals[o] for o in two_row),
      'options %s' % [list(o) for o in two_row])
solver_free = set(forced) | set(two_row)
check('so_18_of_the_45_values_are_solver_free_and_the_other_27_are_not',
      len(solver_free) == 18 and OPTION_COUNT - len(solver_free) == 27,
      '17 by thm:rect and the empty partition, [11,8] added by prop:2rows; the '
      'REFUTING half -- that no option carries 10 -- rests on the other 27')
bad2 = [(c1, c2) for c1 in range(1, 15) for c2 in range(1, c1 + 1)
        if sg_A((c1, c2)) != prop_2rows(c1, c2)]
check('prop_2rows_closed_form_holds_for_all_two_part_partitions_with_c1_at_most_14',
      not bad2, '%d pairs, %d mismatches' % (14 * 15 // 2, len(bad2)))
say('')

# ---------------------------------------------------------------------------
say('=== Step 8: the source paper\'s two printed tables, reproduced')
# ---------------------------------------------------------------------------
p8 = partitions_upto(8)
h8 = [p for p in p8 if heavy(p)]
cl8 = set(conj_class(p) for p in h8)
note('order <= 8: nonempty partitions %d, heavy %d, heavy conjugation classes %d'
     % (len(p8), len(h8), len(cl8)))
check('there_are_66_nonempty_partitions_of_order_at_most_8',
      len(p8) == 66, 'p(1)+...+p(8) = 1+2+3+5+7+11+15+22 = %d' % len(p8))
check('62_of_them_are_heavy', len(h8) == 62)
check('they_fall_into_34_conjugation_classes', len(cl8) == 34)
printed_b = set(conj_class(p) for p in APPENDIX_B_HEAVY)
check('the_34_printed_appendix_B_diagrams_are_pairwise_inequivalent',
      len(printed_b) == len(APPENDIX_B_HEAVY) == 34,
      '%d transcribed entries, %d classes' % (len(APPENDIX_B_HEAVY), len(printed_b)))
check('every_printed_appendix_B_partition_is_heavy_and_of_order_at_most_8',
      all(heavy(p) and sum(p) <= 8 for p in APPENDIX_B_HEAVY))
check('census_classes_equal_the_printed_appendix_B_symmetric_difference_empty',
      cl8 == printed_b, 'both %d classes, symmetric difference %d'
      % (len(cl8), len(cl8 ^ printed_b)))
p9 = partitions_upto(9)
cl9 = set(conj_class(p) for p in p9 if heavy(p))
note('order <= 9: nonempty partitions %d, heavy conjugation classes %d'
     % (len(p9), len(cl9)))
check('a_census_to_order_9_could_not_have_returned_34_classes',
      len(p9) == 96 and len(cl9) == 49 and len(cl9 - cl8) == 15,
      '96 nonempty, 49 classes, 15 new classes at order exactly 9')
p26 = partitions_upto(26)
sg2_cl = set(conj_class(p) for p in p26 if sg_A(p) == 2)
printed_a = set(conj_class(p) for p in APPENDIX_A_SG2)
note('order <= 26: nonempty partitions %d, classes with sg = 2: %d'
     % (len(p26), len(sg2_cl)))
check('there_are_11731_nonempty_partitions_of_order_at_most_26',
      len(p26) == 11731, '%d' % len(p26))
check('exactly_13_conjugation_classes_of_order_at_most_26_have_sg_2',
      len(sg2_cl) == 13 and len(printed_a) == 13)
check('census_classes_equal_the_printed_appendix_A_symmetric_difference_empty',
      sg2_cl == printed_a, 'symmetric difference %d' % len(sg2_cl ^ printed_a))
check('the_printed_appendix_A_orders_are_2_10_12_17_22_23_24_24_24_25_26_26_26',
      sorted(sum(p) for p in APPENDIX_A_SG2)
      == [2, 10, 12, 17, 22, 23, 24, 24, 24, 25, 26, 26, 26])
check('every_printed_appendix_A_partition_really_has_sg_2',
      all(sg_A(p) == 2 for p in APPENDIX_A_SG2))
say('')

# ---------------------------------------------------------------------------
say('=== Step 9: the source paper\'s proved heaviness results still hold')
# ---------------------------------------------------------------------------
check('prop_hook_c_1_to_the_r_minus_1_is_heavy_for_all_r_c_at_most_12',
      all(heavy(norm([c] + [1] * (r - 1)))
          for r in range(1, 13) for c in range(1, 13)),
      '144 hooks')
check('prop_padded_stair_str_c_r_is_heavy_for_all_10_ge_c_ge_r_ge_1',
      all(heavy(staircase(c, r)) for c in range(1, 11) for r in range(1, c + 1)),
      '55 staircases; str_{10,10} alone has order 55')


def interval_thm_pi(p):
    """The right-hand side of thm:partition_interval."""
    if norm(p) == (1,):
        return True
    for r in range(2, len(p) + p[0] + 2):
        lo_r = norm([r, r] + list(range(r - 1, 1, -1)))
        hi_r = tuple([r] * r)
        if young_le(lo_r, p) and young_le(p, hi_r):
            return True
    return False


bad_pi = [p for p in p26 if (sg_A(p) == 1) != interval_thm_pi(p)]
check('thm_partition_interval_holds_as_an_iff_on_all_11731_partitions_of_order_at_most_26',
      not bad_pi, '%d counterexamples to the iff' % len(bad_pi))
check('the_empty_partition_is_the_unique_terminal_position_with_sg_0',
      sg_A(()) == 0 and options_A(()) == set()
      and all(sg_A(p) > 0 for p in partitions_upto(12)),
      'every nonempty partition of order <= 12 has sg >= 1')
say('')

# ---------------------------------------------------------------------------
say('=== Step 10: controls, both polarities, on objects named in advance')
# ---------------------------------------------------------------------------
neg = [(2, 2), (4, 4), (6, 6)]
for p in neg:
    note('forced-NEGATIVE', list(p), 'sg', sg_A(p), 'L', longest_play(p),
         'heavy', heavy(p))
check('forced_negative_controls_are_all_not_heavy_so_the_detector_can_say_no',
      all(not heavy(p) for p in neg),
      '[2,2] sg 1 vs 3, [4,4] sg 3 vs 5, [6,6] sg 5 vs 7')
pos = [(5, 4, 3), (5, 5, 5), (1,), (3, 1, 1)]
for p in pos:
    note('forced-POSITIVE', list(p), 'sg', sg_A(p), 'L', longest_play(p),
         'heavy', heavy(p))
check('forced_positive_controls_are_all_heavy_so_the_detector_can_say_yes',
      all(heavy(p) for p in pos),
      'a staircase, a rectangle, a point and a hook, each proved heavy in the source')
check('the_maximal_length_option_of_lambda_is_itself_heavy',
      heavy((10, 10, 10, 8, 8)) and sg_A((10, 10, 10, 8, 8)) == 14,
      '[10,10,10,8,8] sg 14 = L: the hole sits far below the top of the interval')
check('all_41_other_members_of_cell_10_4_are_heavy',
      all(heavy(p) for p in mem_param if p != norm(LAMBDA)),
      'so the failure is not an artefact of the cell being broken everywhere')
say('')

# ---------------------------------------------------------------------------
say('=== Step 11: the inventory of non-heavy interval members')
# ---------------------------------------------------------------------------


def sweep(x, y, order_cap=None, print_all=False):
    rct = tuple([x + 1] * (y + 1))
    if not admissible(x, y):
        return None
    if not heavy(rct):
        return None
    mem = interval_members(x, y)
    if order_cap is not None:
        mem = [p for p in mem if sum(p) <= order_cap]
    bad = sorted((p for p in mem if not heavy(p)), key=lambda p: (sum(p), p))
    if print_all:
        for p in bad:
            say('NOTE NONHEAVY %-30s cell %s order %3d sg %2d required %2d'
                % (str(list(p)), (x, y), sum(p), sg_A(p), longest_play(p)))
    return len(mem), bad


n104, bad104 = sweep(10, 4, print_all=True)
check('cell_10_4_has_exactly_one_non_heavy_member_and_it_is_lambda',
      n104 == 42 and bad104 == [norm(LAMBDA)],
      '1 of 42')
n114, bad114 = sweep(11, 4, print_all=True)
check('cell_11_4_is_admissible_and_has_exactly_the_two_printed_non_heavy_members',
      n114 == 42 and sorted(bad114) == sorted(CELL_11_4_NONHEAVY),
      '2 of 42, both sg 9 against required 16')
check('cell_8_7_is_admissible_and_its_rectangle_9_8_is_heavy',
      admissible(8, 7) and heavy(tuple([9] * 8)) and sg_A(tuple([9] * 8)) == 16,
      'sg([9^8]) = ((7 XOR 8) + 1) = %d' % (((8 - 1) ^ (9 - 1)) + 1))
n87, bad87 = sweep(8, 7, print_all=True)
check('cell_8_7_has_exactly_40_non_heavy_members_of_its_1430',
      n87 == 1430 == catalan(8) and len(bad87) == 40,
      '40 of 1430; all 40 are printed above, none is left as a bare count')
check('the_second_witness_is_the_least_order_non_heavy_member_of_cell_8_7',
      bad87[0] == norm(SECOND) and sg_A(SECOND) == SECOND_SG
      and longest_play(SECOND) == SECOND_REQUIRED,
      '%s sg %d against required %d'
      % (list(SECOND), sg_A(SECOND), longest_play(SECOND)))
check('the_full_inventory_over_these_three_cells_is_43_objects',
      len(bad104) + len(bad114) + len(bad87) == 43, '1 + 2 + 40 = 43')
# the a = 3 mod 8 family: an observation, NOT a theorem
fam = [(11, 19), (19, 27)]
fam_vals = []
for x in (11, 19, 27):
    p = tuple([x + 1] * 3 + [x] * 2)
    fam_vals.append((x, sg_A(p), longest_play(p)))
    note('family member a =', x, list(p), 'sg', sg_A(p), 'L', longest_play(p))
check('the_family_a_plus_1_cubed_a_squared_is_non_heavy_with_deficit_7_at_a_11_19_27',
      all(admissible(x, 4) and L - s == 7 for x, s, L in fam_vals),
      'a = 11, 19, 27, each a == 3 mod 8; an OBSERVATION, not a proved family')
say('')

# ---------------------------------------------------------------------------
say('=== Step 12: order 49 is minimal, and it is attained twice')
# ---------------------------------------------------------------------------
bad_form = [(x, y) for x in range(0, 30) for y in range(0, x + 1)
            if sum(staircase(x + 1, y + 1)) != lower_endpoint_order(x, y)]
check('the_lower_endpoint_order_formula_b_plus_1_times_2a_plus_2_minus_b_over_2_is_right',
      not bad_form, 'checked against the enumerated staircase on %d cells'
      % (30 * 31 // 2))
cells = [(x, y) for y in range(0, 60) for x in range(y, 400)
         if admissible(x, y) and lower_endpoint_order(x, y) <= MIN_ORDER]
note('admissible cells (a >= b, a & b = 0) whose least member has order <= 49:',
     len(cells))
by_b = {}
for x, y in cells:
    by_b.setdefault(y, []).append(x)
for y in sorted(by_b):
    note('  b = %d : a in %s' % (y, by_b[y]))
check('no_cell_with_b_at_least_8_can_hold_a_member_of_order_at_most_49',
      max(by_b) == 7,
      'least order over admissible cells with b >= 8 is %d, attained at (a,b) = %s'
      % min(((lower_endpoint_order(x, y), (x, y))
             for y in range(8, 60) for x in range(y, 400) if admissible(x, y))))
check('b_equals_0_cells_need_no_sweep_because_the_interval_is_the_rectangle_itself',
      all(interval_members(x, 0) == [(x + 1,)] for x in by_b.get(0, [])),
      '%d cells with b = 0, each interval a single point, heavy by hypothesis'
      % len(by_b.get(0, [])))
found = []
swept = 0
for x, y in cells:
    if y == 0:
        continue
    res = sweep(x, y, order_cap=MIN_ORDER)
    if res is None:
        continue
    swept += 1
    _n, bad = res
    for p in bad:
        found.append((sum(p), p, (x, y), sg_A(p), longest_play(p)))
for o, p, c, s, Lv in sorted(found):
    note('MINIMAL-RANGE NONHEAVY', list(p), 'cell', c, 'order', o, 'sg', s,
         'required', Lv)
check('every_such_cell_with_b_at_least_1_was_swept_over_its_members_of_order_at_most_49',
      swept == len([1 for x, y in cells if y >= 1]),
      '%d cells swept' % swept)
check('exactly_two_non_heavy_members_of_order_at_most_49_exist_over_all_of_them',
      len(found) == 2, '%d found' % len(found))
check('both_have_order_exactly_49_so_49_is_the_minimum',
      all(o == MIN_ORDER for o, _p, _c, _s, _L in found),
      'orders = %s' % [o for o, _p, _c, _s, _L in found])
check('the_two_minimisers_are_the_two_partitions_the_paper_names',
      sorted(p for _o, p, _c, _s, _L in found) == sorted(map(norm, MIN_ORDER_ATTAINED_BY)),
      '%s and %s' % (list(LAMBDA), list(SECOND)))
check('so_the_banked_witness_TIES_for_least_order_and_is_not_the_unique_minimiser',
      len(found) == 2 and norm(LAMBDA) in [p for _o, p, _c, _s, _L in found])
say('')

# ---------------------------------------------------------------------------
say('=== Scope')
# ---------------------------------------------------------------------------
say('NOT RE-RUN: (1) HEAVINESS OF WHOLE CELLS OUTSIDE THE ORDER-49 WINDOW. Step 12')
say('  sweeps every admissible cell whose least member has order <= 49, but only over')
say('  that cell\'s members OF ORDER <= 49. Cells with b >= 8, and the parts of cells')
say('  with b <= 7 lying above order 49, are NOT claimed clean -- they are only shown')
say('  unable to contain a counterexample of order <= 49. Only cells (10,4), (11,4)')
say('  and (8,7) are swept in full (Step 11), so the inventory of 43 is a lower bound')
say('  on the number of non-heavy interval members, not a total.')
say('NOT RE-RUN: (2) THE INFINITE FAMILY. The last check of Step 11 evaluates')
say('  [(a+1)^3, a^2] at a = 11, 19, 27 only. That every a == 3 (mod 8) gives deficit')
say('  exactly 7 is an OBSERVATION suggested by three data points and is NOT proved')
say('  here; the paper states it as a conjecture of ours and nothing else in this')
say('  program depends on it.')
say('NOT RE-RUN: (3) THE SOURCE TEXT AND THE PRIOR-ART SEARCH. The two appendix tables')
say('  and the closed forms of thm:rect, prop:bound, prop:hook, prop:padded_stair,')
say('  prop:2rows and thm:partition_interval are TRANSCRIBED from arXiv1.tex of')
say('  arXiv:2506.04991v1 into the constants at the head of this file; this program')
say('  cannot fetch that file and does not verify the transcription. Nor does it')
say('  check any literature claim: that no published theorem covers a 5-row order-49')
say('  partition is an editorial finding, not a computation.')
say('NOT RE-RUN: (4) THE FINITE RANGES ARE FINITE. thm:rect is checked for r,c <= 9')
say('  (and its closed form used, not reproved, elsewhere), prop:hook for r,c <= 12,')
say('  prop:padded_stair for c,r <= 10, prop:2rows for c1 <= 14, obs:lucas as a')
say('  binomial identity for a,b <= 64, thm:partition_interval as an iff to order 26,')
say('  and the two-engine agreement over the reachable set of lambda plus every')
say('  partition of order <= 9. Nothing outside those ranges is asserted.')
say('')

n = _passed + _failed
if _failed:
    say('VERDICT: %d of %d CHECKS FAILED' % (_failed, n))
    sys.exit(1)
say('VERDICT: ALL %d CHECKS PASS' % n)
sys.exit(0)

#!/usr/bin/env python3
"""Independent verification of three small counterexamples to a trace
conjecture of Frankl and Wang.

Setting.  For a family F of subsets of [n] and Y subset of [n], the trace is
F|_Y = {S cap Y : S in F}.  The arrow relation (n,m) -> (a,b) says every
F with |F| >= m admits an a-set Y with |F|_Y| >= b, and
m(n,a,b) = min{m : (n,m) -> (a,b)}.  Since (n,m) -> (a,b) fails exactly when
some family of size >= m has all a-traces of size <= b-1, one has
m(n,a,b) = 1 + max{|F| : every a-trace of F has size <= b-1}.

VALUES TAKEN FROM THE PAPER (inputs; nothing here is assumed correct)
  * the conjectured formula  m(n,l+1, 3*2^(l-1)+1) = 1 + prod_{0<=i<l}
    floor((n+l+i)/l)  for all n > l > 0;
  * the six triples D = {136,145,156,234,236,245} and the three 4-sets
    A = {1235,1246,3456} on [6], and F6 = C([6],<=2) u (C([6],3)\\D) u A;
  * the seven triples L = {123,145,167,246,257,347,356} on [7] and
    F4 = C([7],<=3) \\ L;
  * C = {1,4}, M = {123,456,147} on [7], Q = {4-sets meeting C and
    containing no member of M}, F5 = C([7],<=2) u (C([7],3)\\M) u Q;
  * the claimed values m(6,5,25) = 40, m(7,5,25) >= 58, m(7,6,49) >= 80,
    against formula values 37, 55, 73;
  * one external ingredient, the squashing lemma of the cited literature
    (an arbitrary family can be replaced by a down-set of the same size
    with no larger traces).  It is not proved here; it is exercised on
    random instances, and the direct random test of the upper bound does
    not use it at all.  It is relied on TWICE below: for the upper bound
    m(6,5,25) <= 40, and for the auxiliary equalities m(5,3,7) = 13 and
    m(5,4,13) = 19, whose searches range over down-sets of 2^[5] only.

DERIVED HERE (every number below is recomputed from the definitions)
  * the three exhibited families are decoded, counted, and tested for the
    down-set property and for every structural claim made about D, A, L,
    M and Q (regularity, non-containment, Fano axioms, disjoint cones);
  * every trace size is computed from the definition {S cap Y : S in F},
    never from a formula: all six 5-traces of F6 equal 24, all twenty-one
    5-traces of F4 equal 24, all seven 6-traces of F5 equal 48;
  * the resulting lower bounds 40, 58, 80 and the comparison with the
    formula values 37, 55, 73 -- the refutation itself;
  * the matching upper bound m(6,5,25) <= 40 by EXHAUSTIVE census of every
    down-set of size 40 in 2^[6] (none has all 5-traces <= 24), with the
    same search at size 39 returning 60 families, one of which is F6;
  * the paper's proof steps, independently: the Kruskal-Katona step is
    replaced by an exhaustive minimisation of the upper 4-shadow over all
    125970 eight-element families of triples of [6], the upper cone sizes
    64/32/16/8 that exclude small sets from the complementary up-set, and
    the double counting step by an exhaustive scan of the admissible (t,q);
  * the tautology of the formula at n = l+1; the auxiliary values
    m(3,3,7) = 7 and m(4,3,7) = 10 by unrestricted brute force over
    every one of the 2^(2^n) families; and m(5,3,7) = 13 and
    m(5,4,13) = 19 as maxima over all 7581 down-sets of 2^[5] -- there
    the ">=" direction is unconditional but the equality needs the
    quoted squashing lemma.  All four agree with the formula, so the
    formula's failure begins exactly where the paper says it does.

NOT RE-RUN: the squashing lemma is not proved (only random-tested at n=6
and confirmed by exhaustive comparison at n=4); the auxiliary values
m(5,3,7) and m(5,4,13) are maxima over down-sets only, so those two rest
on the same lemma; and no census is run for n=7, where the paper claims
only lower bounds.
"""

import random
import sys
from itertools import combinations

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return bool(ok)


def bits(*labels):
    """Bitmask of a set given by 1-based point labels."""
    m = 0
    for p in labels:
        m |= 1 << (p - 1)
    return m


def size(m):
    return bin(m).count("1")


def levels(n, k):
    """All subsets of [n] of size k, as bitmasks."""
    return [sum(1 << (p - 1) for p in c) for c in combinations(range(1, n + 1), k)]


def upto(n, k):
    out = []
    for j in range(k + 1):
        out.extend(levels(n, j))
    return out


def trace(fam, ymask):
    """The trace family {S cap Y : S in fam}, from the definition."""
    return set(s & ymask for s in fam)


def is_downset(fam):
    fam = set(fam)
    for s in fam:
        b = s
        while b:
            low = b & -b
            if (s ^ low) not in fam:
                return False
            b ^= low
    return True


def fw_rhs(n, ell):
    """Right-hand side of the conjectured formula."""
    p = 1
    for i in range(ell):
        p *= (n + ell + i) // ell
    return 1 + p


def fw_b(ell):
    return 3 * 2 ** (ell - 1) + 1


def max_trace(fam, n, a):
    return max(len(trace(fam, y)) for y in levels(n, a))


D_TRIPLES = [(1, 3, 6), (1, 4, 5), (1, 5, 6), (2, 3, 4), (2, 3, 6), (2, 4, 5)]
A_QUADS = [(1, 2, 3, 5), (1, 2, 4, 6), (3, 4, 5, 6)]


def build_f6():
    dfam = [bits(*t) for t in D_TRIPLES]
    afam = [bits(*q) for q in A_QUADS]
    fam = set(upto(6, 2))
    fam |= set(t for t in levels(6, 3) if t not in dfam)
    fam |= set(afam)
    return fam, dfam, afam


def profile(fam, n):
    return [sum(1 for s in fam if size(s) == k) for k in range(n + 1)]


def show(fam, n=7):
    """Print a family back in the paper's notation, e.g. 136,145,..."""
    out = []
    for s in sorted(fam, key=lambda m: (size(m), m)):
        out.append("".join(str(p) for p in range(1, n + 1) if (s >> (p - 1)) & 1))
    return ",".join(out)


def check_f6_structure(fam, dfam, afam):
    ck("f6_members_are_subsets_of_[6]", all(0 <= s < 64 for s in fam))
    ck("f6_cardinality_39", len(fam) == 39, "|F6|=%d" % len(fam))
    prof = profile(fam, 6)
    ck("f6_level_profile_1_6_15_14_3_0_0", prof == [1, 6, 15, 14, 3, 0, 0],
       "by size: " + ",".join(str(v) for v in prof))
    ck("f6_is_a_down_set", is_downset(fam))
    deg_d = [sum(1 for t in dfam if (t >> x) & 1) for x in range(6)]
    ck("D_is_six_triples_3_regular",
       len(set(dfam)) == 6 and all(size(t) == 3 for t in dfam)
       and deg_d == [3] * 6, "degrees " + str(deg_d))
    deg_a = [sum(1 for q in afam if (q >> x) & 1) for x in range(6)]
    ck("A_is_three_4_sets_2_regular",
       len(set(afam)) == 3 and all(size(q) == 4 for q in afam)
       and deg_a == [2] * 6, "degrees " + str(deg_a))
    ck("no_member_of_D_lies_in_a_member_of_A",
       all((t & q) != t for t in dfam for q in afam))
    miss = set(levels(6, 3)) - fam
    tops = set(s for s in fam if size(s) == 4)
    ck("f6_decodes_back_to_missing_triples_D_and_top_layer_A",
       miss == set(dfam) and tops == set(afam),
       "omitted triples " + show(miss, 6) + " ; 4-sets " + show(tops, 6))
    ok = all(trace(fam, y) == set(s for s in fam if s & ~y == 0)
             for y in levels(6, 5))
    ck("down_set_trace_identity_for_F6", ok)


def check_f6_traces(fam):
    sizes = {}
    decomp = set()
    for y in levels(6, 5):
        tr = trace(fam, y)
        sizes[y] = len(tr)
        decomp.add(tuple(sum(1 for s in tr if size(s) == k) for k in (0, 1, 2, 3, 4)))
    vals = sorted(set(sizes.values()))
    ck("f6_every_5_trace_has_size_24", vals == [24], "trace sizes " + str(vals))
    ck("f6_trace_decomposition_16_7_1",
       decomp == {(1, 5, 10, 7, 1)}, str(sorted(decomp)))
    top = max(sizes.values())
    ck("m_6_5_25_lower_bound_is_40", top <= 24 and len(fam) + 1 == 40,
       "max 5-trace %d < 25 on %d sets" % (top, len(fam)))
    # The bound is released only if the exhibited family really is admissible;
    # a cardinality alone proves nothing, so a bad trace returns 0.
    return len(fam) + 1 if top <= 24 else 0


def census_downsets_6(fsize, cap=24):
    """Every down-set of 2^[6] of size fsize whose six 5-traces are <= cap.

    Enumerates the complementary up-set G (size 64-fsize) top down: a set may
    join G only once all its immediate supersets are in G, which is exactly
    up-closure.  For a down-set F, |F|_{Y}| = |F cap 2^Y|, so the trace
    condition reads |G cap 2^{[6]\\{x}}| >= 32 - cap for every point x.
    """
    need = 32 - cap
    gsize = 64 - fsize
    order = sorted(range(64), key=lambda m: -size(m))
    pos = dict((m, i) for i, m in enumerate(order))
    sup = [[pos[m | (1 << b)] for b in range(6) if not (m >> b) & 1] for m in order]
    rest = [[0] * 6 for _ in range(65)]
    for i in range(63, -1, -1):
        for x in range(6):
            rest[i][x] = rest[i + 1][x] + (0 if (order[i] >> x) & 1 else 1)
    ing = [False] * 64
    cnt = [0] * 6
    found = []

    def rec(i, count):
        if count == gsize:
            if min(cnt) >= need:
                found.append(frozenset(m for m in range(64) if not ing[m]))
            return
        if i == 64 or count + 64 - i < gsize:
            return
        for x in range(6):
            if cnt[x] + rest[i][x] < need:
                return
        m = order[i]
        if all(ing[order[j]] for j in sup[i]):
            ing[m] = True
            for x in range(6):
                if not (m >> x) & 1:
                    cnt[x] += 1
            rec(i + 1, count + 1)
            ing[m] = False
            for x in range(6):
                if not (m >> x) & 1:
                    cnt[x] -= 1
        rec(i + 1, count)

    rec(0, 0)
    return found


def check_census():
    fam6 = build_f6()[0]
    at39 = census_downsets_6(39)
    good = all(len(f) == 39 and is_downset(f) and max_trace(f, 6, 5) == 24
               for f in at39)
    ck("census_at_39_recovers_the_exhibited_family",
       len(at39) >= 1 and frozenset(fam6) in at39 and good,
       "%d extremal down-sets of size 39, all re-validated" % len(at39))
    at40 = census_downsets_6(40)
    ck("exhaustive_census_no_40_member_down_set_has_all_5_traces_le_24",
       len(at40) == 0, "%d survivors" % len(at40))
    control = census_downsets_6(40, cap=25)
    ok = all(len(f) == 40 and is_downset(f) and max_trace(f, 6, 5) == 25
             for f in control)
    ck("positive_control_the_same_search_at_size_40_is_not_vacuous",
       len(control) > 0 and ok,
       "%d down-sets of size 40 with all 5-traces <= 25, all re-validated"
       % len(control))
    ck("m_6_5_25_equals_40",
       len(at40) == 0 and len(at39) > 0,
       "max admissible down-set size 39, so m=40 (squashing lemma quoted)")


def binom(n, k):
    if k < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r


def check_kruskal_katona_step():
    """Exhaustive replacement for the Kruskal-Katona appeal in the paper."""
    triples = levels(6, 3)
    quads = levels(6, 4)
    up = dict((t, [q for q in quads if q & t == t]) for t in triples)
    down = dict((t, [p for p in levels(6, 2) if p & t == p]) for t in triples)
    best = None
    mismatch = 0
    total = 0
    for fam in combinations(triples, 8):
        total += 1
        shadow = set()
        for t in fam:
            shadow.update(up[t])
        lower = set()
        for t in fam:
            lower.update(down[63 ^ t])
        if len(lower) != len(shadow):
            mismatch += 1
        if best is None or len(shadow) < best:
            best = len(shadow)
    ck("kk_cascade_identities_8_and_10",
       binom(4, 3) + binom(3, 2) + binom(1, 1) == 8
       and binom(4, 2) + binom(3, 1) + binom(1, 0) == 10)
    ck("upper_4_shadow_of_any_8_triples_of_[6]_is_at_least_10",
       best == 10 and total == binom(20, 8),
       "minimum %s over all %d eight-element families" % (best, total))
    ck("complementation_matches_upper_and_lower_shadows",
       mismatch == 0, "%d discrepancies" % mismatch)


def check_upper_cone_sizes():
    """The two cone counts the upper-bound argument quotes, recomputed.

    G is an up-set of size 64-40 = 24.  A set of size <= 1 cannot lie in G
    because its upper cone is already larger than 24, and if a pair lies in G
    then exactly 24 - 16 = 8 members of G lie outside its cone.
    """
    cone = dict((k, sorted(set(sum(1 for u in range(64) if u & s == s)
                               for s in levels(6, k)))) for k in range(4))
    gsize = 64 - 40
    ok = (cone[0] == [64] and cone[1] == [32] and cone[2] == [16]
          and cone[3] == [8]
          and min(cone[1]) > gsize
          and gsize - min(cone[2]) == 8)
    ck("upper_cone_of_a_point_exceeds_|G|_and_a_pair_leaves_exactly_8", ok,
       "cone sizes by base size " + str([cone[k] for k in range(4)])
       + ", |G|=%d" % gsize)


def check_double_counting_step():
    """G is an up-set of size 24, all members of size >= 3, all 5-sets and
    [6] present; t triples, q quadruples.  The paper's two inequalities are
    re-derived by scanning every admissible (t,q)."""
    bad = []
    for t in range(0, 21):
        for q in range(0, 16):
            if t + q != 17:
                continue
            if 3 * t + 2 * q + 6 < 48:
                continue
            if not (t >= 8 and q <= 9):
                bad.append((t, q))
    ck("counting_forces_t_ge_8_and_q_le_9", not bad, "violations " + str(bad))
    fails = 0
    rng = random.Random(20240501)
    triples = levels(6, 3)
    quads = levels(6, 4)
    for _ in range(400):
        tt = rng.sample(triples, rng.randint(0, 12))
        gfam = set(tt) | set(levels(6, 5)) | set([63])
        for t in tt:
            gfam.update(q for q in quads if q & t == t)
        extra = [q for q in quads if q not in gfam]
        gfam.update(rng.sample(extra, rng.randint(0, len(extra))))
        lhs = sum(len(set(s for s in gfam if (s >> x) & 1 == 0)) for x in range(6))
        tq = profile(gfam, 6)
        # The paper's premise -- no member of size at most two, all six 5-sets
        # and [6] present -- must be asserted, and the trailing constant must
        # be the literal 6, or the equation degenerates into a tautology that
        # holds for every family whatever its level profile.
        premise = (tq[0] == tq[1] == tq[2] == 0 and tq[5] == 6 and 63 in gfam)
        if not premise or lhs != 3 * tq[3] + 2 * tq[4] + 6:
            fails += 1
    ck("double_counting_identity_sum_x_|G_cap_2^Yx|_eq_3t+2q+6",
       fails == 0, "%d failures over 400 random up-sets" % fails)


FANO = [(1, 2, 3), (1, 4, 5), (1, 6, 7), (2, 4, 6),
        (2, 5, 7), (3, 4, 7), (3, 5, 6)]


def build_f4():
    lines = [bits(*t) for t in FANO]
    fam = set(upto(7, 3)) - set(lines)
    return lines, fam


def check_fano(lines):
    ck("L_is_seven_triples_of_[7]",
       len(set(lines)) == 7 and all(size(t) == 3 and t < 128 for t in lines))
    pair_hits = [sum(1 for t in lines if t & p == p) for p in levels(7, 2)]
    ck("every_pair_of_points_lies_on_exactly_one_line",
       set(pair_hits) == {1}, "counts " + str(sorted(set(pair_hits))))
    meets = set(size(a & b) for a, b in combinations(lines, 2))
    ck("two_distinct_lines_meet_in_exactly_one_point",
       meets == {1}, "intersection sizes " + str(sorted(meets)))
    inside, hitting = set(), set()
    for p in levels(7, 2):
        y = 127 ^ p
        inside.add(sum(1 for t in lines if t & y == t))
        hitting.add(sum(1 for t in lines if t & p))
    ck("each_pair_meets_5_lines_leaving_2_lines_inside_its_complement",
       hitting == {5} and inside == {2},
       "meeting " + str(sorted(hitting)) + ", inside " + str(sorted(inside)))


def check_f4(fam, lines):
    ck("f4_cardinality_57", len(fam) == 64 - 7 == 57, "|F4|=%d" % len(fam))
    prof = profile(fam, 7)
    ck("f4_level_profile_1_7_21_28", prof == [1, 7, 21, 28, 0, 0, 0, 0],
       "by size: " + ",".join(str(v) for v in prof))
    ck("f4_is_a_down_set", is_downset(fam))
    miss = set(levels(7, 3)) - fam
    ck("f4_omits_exactly_the_seven_Fano_lines", miss == set(lines),
       "omitted triples " + show(miss))
    sizes = sorted(set(len(trace(fam, y)) for y in levels(7, 5)))
    ck("f4_every_one_of_21_5_traces_has_size_24",
       sizes == [24] and len(levels(7, 5)) == 21, "trace sizes " + str(sizes))
    top = max(len(trace(fam, y)) for y in levels(7, 5))
    ck("m_7_5_25_lower_bound_is_58", top <= 24 and len(fam) + 1 == 58,
       "max 5-trace %d < 25 on %d sets" % (top, len(fam)))
    return len(fam) + 1 if top <= 24 else 0


CSET = (1, 4)
M_TRIPLES = [(1, 2, 3), (4, 5, 6), (1, 4, 7)]


def build_f5():
    cmask = bits(*CSET)
    mfam = [bits(*t) for t in M_TRIPLES]
    qfam = set(q for q in levels(7, 4)
               if (q & cmask) and all(q & m != m for m in mfam))
    fam = set(upto(7, 2))
    fam |= set(t for t in levels(7, 3) if t not in mfam)
    fam |= qfam
    return qfam, fam


def check_q_structure(qfam):
    cmask = bits(*CSET)
    mfam = [bits(*t) for t in M_TRIPLES]
    quads = levels(7, 4)
    meeting = [q for q in quads if q & cmask]
    ck("four_sets_meeting_C_number_30",
       len(meeting) == 35 - binom(5, 4) == 30, "%d of %d" % (len(meeting), len(quads)))
    cones = [set(q for q in quads if q & m == m) for m in mfam]
    ck("each_member_of_M_has_4_four_set_supersets_all_meeting_C",
       all(len(c) == 4 for c in cones)
       and all(q & cmask for c in cones for q in c),
       "cone sizes " + str([len(c) for c in cones]))
    ck("the_three_cones_are_pairwise_disjoint_since_unions_are_big",
       all(not (a & b) for a, b in combinations(cones, 2))
       and all(size(a | b) >= 5 for a, b in combinations(mfam, 2)))
    ck("q_cardinality_18", len(qfam) == 30 - 12 == 18, "|Q|=%d" % len(qfam))


def check_f5(fam, qfam):
    cmask = bits(*CSET)
    mfam = [bits(*t) for t in M_TRIPLES]
    ck("f5_cardinality_79", len(fam) == 29 + 32 + 18 == 79, "|F5|=%d" % len(fam))
    prof = profile(fam, 7)
    ck("f5_level_profile_1_7_21_32_18", prof == [1, 7, 21, 32, 18, 0, 0, 0],
       "by size: " + ",".join(str(v) for v in prof))
    ck("f5_is_a_down_set", is_downset(fam))
    miss = set(levels(7, 3)) - fam
    tops = set(s for s in fam if size(s) == 4)
    ck("f5_decodes_back_to_missing_triples_M_and_top_layer_Q",
       miss == set(mfam) and tops == set(qfam),
       "omitted triples " + show(miss) + " ; Q = " + show(qfam))
    detail = []
    ok = True
    for x in range(1, 8):
        y = 127 ^ bits(x)
        inm = sum(1 for m in mfam if m & y == m)
        inq = sum(1 for q in qfam if q & y == q)
        want = (1, 7) if x in CSET else (2, 8)
        detail.append("x=%d:(%d,%d)" % (x, inm, inq))
        if (inm, inq) != want:
            ok = False
    ck("per_point_counts_of_M_and_Q_inside_Y_x", ok, " ".join(detail))
    obs = []
    for x in range(1, 8):
        y = 127 ^ bits(x)
        inside = [q for q in levels(7, 4) if q & y == q]
        hits = [q for q in inside if q & cmask]
        cones = [[q for q in hits if q & m == m] for m in mfam if m & y == m]
        obs.append((len(inside), len(hits), tuple(sorted(len(c) for c in cones)),
                    max(len(set.intersection(*[set(c) for c in cones])) if len(cones) > 1 else 0, 0)))
    want_in = [(15, 10, (3,), 0) if x in CSET else (15, 14, (3, 3), 0)
               for x in range(1, 8)]
    ck("counts_of_four_sets_in_Y_x_meeting_C_and_over_each_M",
       obs == want_in, "observed " + str(obs))
    sizes = sorted(set(len(trace(fam, y)) for y in levels(7, 6)))
    ck("f5_every_one_of_7_6_traces_has_size_48",
       sizes == [48] and len(levels(7, 6)) == 7, "trace sizes " + str(sizes))
    top = max(len(trace(fam, y)) for y in levels(7, 6))
    ck("m_7_6_49_lower_bound_is_80", top <= 48 and len(fam) + 1 == 80,
       "max 6-trace %d < 49 on %d sets" % (top, len(fam)))
    return len(fam) + 1 if top <= 48 else 0


def check_conjecture(bound6, bound7a, bound7b):
    """The refutation: compare the computed bounds with the formula."""
    params = [(6, 4, bound6), (7, 4, bound7a), (7, 5, bound7b)]
    shape = [(ell + 1, fw_b(ell)) for _, ell, _ in params]
    ck("conjecture_instances_have_a_and_b_equal_to_5_25__5_25__6_49",
       shape == [(5, 25), (5, 25), (6, 49)], str(shape))
    rhs = [fw_rhs(n, ell) for n, ell, _ in params]
    ck("formula_values_are_37_55_73", rhs == [37, 55, 73], str(rhs))
    verdict = []
    for (n, ell, low), r in zip(params, rhs):
        verdict.append((n, ell, low, r, low > r))
    # Each low was released by its own trace test; an inadmissible family
    # arrives here as 0 and cannot manufacture a refutation from cardinality.
    ck("computed_bounds_exceed_the_formula_so_it_is_false",
       all(v[4] for v in verdict),
       "; ".join("m(%d,%d,%d)>=%d vs formula %d" %
                 (n, ell + 1, fw_b(ell), low, r) for n, ell, low, r, _ in verdict))
    # Which instances are actually refuted must be READ OFF the computed
    # bounds, not off the literal parameter list, or the check is a tautology.
    refuted = sorted((n, ell) for (n, ell, low), r in zip(params, rhs) if low > r)
    firstfail = sorted(set(ell for n, ell in refuted if n == ell + 2))
    ck("failure_reaches_n_eq_l_plus_2_for_l_eq_4_and_5",
       refuted == [(6, 4), (7, 4), (7, 5)] and firstfail == [4, 5]
       and all(fw_rhs(ell + 1, ell) == fw_b(ell) for ell in (4, 5)),
       "refuted instances (n,l) " + str(refuted)
       + "; smallest n-l attained " + str(min([n - ell for n, ell in refuted])
                                          if refuted else None))
    taut = all(fw_rhs(ell + 1, ell) == fw_b(ell) for ell in range(1, 15))
    ck("formula_is_tautological_at_n_eq_l_plus_1", taut,
       "checked l=1..14: rhs equals the target trace size")


def all_families_max(n, a, cap):
    """Largest family of subsets of [n] all of whose a-traces are <= cap.

    Brute force over every one of the 2^(2^n) families; only for n <= 4.
    """
    ys = levels(n, a)
    best = -1
    for mask in range(1 << (1 << n)):
        if bin(mask).count("1") <= best:
            continue
        fam = [s for s in range(1 << n) if (mask >> s) & 1]
        if all(len(trace(fam, y)) <= cap for y in ys):
            best = len(fam)
    return best


def all_downsets(n):
    order = sorted(range(1 << n), key=size)
    low = [[m ^ (1 << b) for b in range(n) if (m >> b) & 1] for m in order]
    inf = [False] * (1 << n)
    out = []

    def rec(i, cur):
        if i == len(order):
            out.append(tuple(cur))
            return
        rec(i + 1, cur)
        m = order[i]
        if all(inf[s] for s in low[i]):
            inf[m] = True
            cur.append(m)
            rec(i + 1, cur)
            cur.pop()
            inf[m] = False

    rec(0, [])
    return out


def downset_max(n, a, cap, dsets):
    ys = levels(n, a)
    best = -1
    for fam in dsets:
        if len(fam) <= best:
            continue
        if all(len(set(s & y for s in fam)) <= cap for y in ys):
            best = len(fam)
    return best


def check_small_cases():
    ds4 = all_downsets(4)
    ds5 = all_downsets(5)
    ck("down_set_counts_of_2^[4]_and_2^[5]_are_168_and_7581",
       len(ds4) == 168 and len(ds5) == 7581,
       "%d and %d" % (len(ds4), len(ds5)))
    full4 = all_families_max(4, 3, 6)
    ds_only4 = downset_max(4, 3, 6, ds4)
    ck("squashing_reduction_is_exact_at_n_eq_4",
       full4 == ds_only4 == 9,
       "unrestricted max %d, down-set max %d" % (full4, ds_only4))
    m437 = full4 + 1
    m537 = downset_max(5, 3, 6, ds5) + 1
    m5413 = downset_max(5, 4, 12, ds5) + 1
    got = [m437, m537, m5413]
    want = [fw_rhs(4, 2), fw_rhs(5, 2), fw_rhs(5, 3)]
    # Scope disclosure.  Only m(4,3,7) is an unrestricted brute force (every
    # one of the 2^16 families of 2^[4]).  The two n=5 values search the
    # down-sets of 2^[5] only, so what they establish unconditionally is
    # m(5,3,7) >= 13 and m(5,4,13) >= 19; the equalities asserted here are
    # contingent on the same quoted squashing lemma as the m(6,5,25) upper
    # bound, which is why the check name and the detail say so.
    scope = ("; m(4,3,7) unrestricted over all %d families of 2^[4], "
             "m(5,3,7) and m(5,4,13) maxima over the %d down-sets of 2^[5] "
             "so equality there is contingent on the quoted squashing lemma"
             % (1 << (1 << 4), len(ds5)))
    ck("formula_holds_for_(n,l)_eq_(4,2)_(5,2)_(5,3)_n_eq_5_over_down_sets_only",
       got == want,
       "computed " + str(got) + " vs formula " + str(want) + scope)
    m337 = all_families_max(3, 3, 6) + 1
    ck("tautological_case_m_3_3_7_equals_7", m337 == 7 == fw_rhs(3, 2),
       "m(3,3,7)=%d" % m337)


def compress(fam, n):
    """Down-compress a family to a down-set of the same cardinality."""
    cur = set(fam)
    changed = True
    while changed:
        changed = False
        for b in range(n):
            nxt = set()
            for s in cur:
                t = s ^ (1 << b)
                if (s >> b) & 1 and t not in cur:
                    nxt.add(t)
                    changed = True
                else:
                    nxt.add(s)
            cur = nxt
    return cur


def check_squashing_random():
    rng = random.Random(88)
    pool = list(range(64))
    ys5 = levels(6, 5)
    bad_size = bad_down = bad_trace = 0
    small = 64
    for _ in range(300):
        fam = rng.sample(pool, 40)
        comp = compress(fam, 6)
        if len(comp) != 40:
            bad_size += 1
        if not is_downset(comp):
            bad_down += 1
        for y in ys5:
            if len(trace(comp, y)) > len(trace(fam, y)):
                bad_trace += 1
        small = min(small, max(len(trace(comp, y)) for y in ys5))
    ck("compression_keeps_size_reaches_a_down_set_and_never_grows_a_trace",
       bad_size == 0 and bad_down == 0 and bad_trace == 0,
       "size %d, non-down-set %d, trace increases %d"
       % (bad_size, bad_down, bad_trace))
    ck("compressed_40_member_down_sets_agree_with_the_census",
       small >= 25, "smallest maximal 5-trace after compression: %d" % small)
    worst = 64
    for _ in range(20000):
        fam = rng.sample(pool, 40)
        worst = min(worst, max(len(trace(fam, y)) for y in ys5))
    ck("random_40_member_families_always_have_a_5_trace_of_at_least_25",
       worst >= 25,
       "smallest maximal 5-trace seen over 20000 random families: %d" % worst)


def main():
    fam6, dfam, afam = build_f6()
    check_f6_structure(fam6, dfam, afam)
    bound6 = check_f6_traces(fam6)
    check_census()
    check_kruskal_katona_step()
    check_upper_cone_sizes()
    check_double_counting_step()
    lines, fam4 = build_f4()
    check_fano(lines)
    bound7a = check_f4(fam4, lines)
    qfam, fam5 = build_f5()
    check_q_structure(qfam)
    bound7b = check_f5(fam5, qfam)
    check_conjecture(bound6, bound7a, bound7b)
    check_small_cases()
    check_squashing_random()
    print("NOTE not re-run: the squashing lemma (down-set reduction) is "
          "quoted from the cited literature, not proved; it is exercised on "
          "300 random 40-member families at n=6 and confirmed exhaustively "
          "at n=4. Two results above rest on it. First, the upper bound "
          "m(6,5,25)<=40, which uses that one quoted lemma plus the "
          "exhaustive down-set census run above. Second, the auxiliary "
          "equalities m(5,3,7)=13 and m(5,4,13)=19, whose searches range "
          "over the 7581 down-sets of 2^[5] only, so without the lemma they "
          "are only the lower bounds m(5,3,7)>=13 and m(5,4,13)>=19; of the "
          "auxiliary values, just m(3,3,7)=7 and m(4,3,7)=10 come from "
          "unrestricted brute force over every family. The three lower "
          "bounds 40, 58, 80, which are what refute the formula, use no "
          "external input. No census is attempted at n=7.")


def report():
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    main()
    sys.exit(report())

#!/usr/bin/env python3
"""Verification of "A Proof that s_2(3) = 13, a Value Asserted without Proof by Sun".

The only inputs are the two objects PRINTED IN THE PAPER, transcribed below as text blocks
and parsed:

  * the extremal multiset W of Section 5 (twelve elements of (Z/9)^2), and
  * Table 1 of Section 3: the four antipodal pairs, the eight origin-missing lines of
    AG(2,3), and the eight linear equations E1..E8 in the paper's own naming
    A, A', B, B', C, C', D, D' of the frozen classes.

Nothing else is read: no witness file, no solver, no network, no third-party module.
Python 3.9+, standard library only. All arithmetic is on integers mod 3 and mod 9; there is
no floating point anywhere in this program and no decision depends on one.

What is re-derived here:

  (a) the combinatorics of F_3^2 that Fact 1 asserts (72 admissible classes, 8 origin-missing
      lines, incidence counts, no antipodal pair on a line), and that Table 1 as printed is
      exactly those 8 lines;
  (b) every claim the paper makes about W -- admissible, BAD, size 12, class multiplicity 2,
      support of size 4, no full line, its four candidates and their sums, and MAXIMALITY
      against all 72 admissible one-element extensions -- which is the lower half s_2(3) >= 13;
  (c) two directed negatives, so a "bad = False" answer means something: the must-fire triple
      of Fact 2(a), and the |I| <= n misreading of the definition, under which W is NOT bad;
  (d) the exhaustion of Section 4: the master bound 3|T| - |A(T)| over ALL 2^8 = 256 supports,
      which of them it settles, and which it does not;
  (e) the two eliminations of Section 4 as displayed, as identities over all 9^3 lifts; and
  (f) STRONGER THAN THE PAPER: for EVERY support with |T| >= 7 -- all eight of size 7, not
      only the one the paper relabels to, plus the one of size 8 -- the full-line system mod 9
      is shown to have NO solution at all, by reducing it to a linear system over F_3 and
      solving that exhaustively. That is the upper half s_2(3) <= 13, given Facts 2 and 3.

What is NOT covered is printed by the program itself as a closing SCOPE note.
"""
import itertools
import re
import sys

N, MOD = 3, 9

# ---------------------------------------------------------------------------
# THE OBJECTS, EXACTLY AS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------

# Section 5, the displayed multiset W.
W_PRINTED = """
W = { (1,0), (4,0), (4,0), (8,0), (5,0), (5,0),
      (0,1), (0,4), (0,4), (0,8), (0,5), (0,5) }
"""

# Section 3, Table 1: the four antipodal pairs of F_3^2 \ {0}.
PAIRS_PRINTED = """
P1 = { (1,0), (2,0) }
P2 = { (0,1), (0,2) }
P3 = { (1,1), (2,2) }
P4 = { (1,2), (2,1) }
"""

# Section 3, Table 1: the eight origin-missing lines with the paper's names for the frozen
# classes, A = x_(1,0), A' = x_(2,0), B = x_(0,1), B' = x_(0,2), C = x_(1,1), C' = x_(2,2),
# D = x_(1,2), D' = x_(2,1).  Each row is  <name> : <three points> : <three class symbols>.
LINES_PRINTED = """
E1 : (1,0) (0,1) (2,2) : A  B  C'
E2 : (1,0) (0,2) (2,1) : A  B' D'
E3 : (1,0) (1,1) (1,2) : A  C  D
E4 : (2,0) (0,2) (1,1) : A' B' C
E5 : (2,0) (0,1) (1,2) : A' B  D
E6 : (2,0) (2,2) (2,1) : A' C' D'
E7 : (0,1) (1,1) (2,1) : B  C  D'
E8 : (0,2) (2,2) (1,2) : B' C' D
"""

# Section 3, Table 1: which point of F_3^2 each class symbol names.
SYMBOL_POINT = {"A": (1, 0), "A'": (2, 0), "B": (0, 1), "B'": (0, 2),
                "C": (1, 1), "C'": (2, 2), "D": (1, 2), "D'": (2, 1)}


def parse_vectors(block):
    """Every '(x,y)' in a printed block, in order of appearance, as a list of pairs."""
    return [(int(a), int(b)) for a, b in re.findall(r"\((\d+)\s*,\s*(\d+)\)", block)]


def parse_lines_table(block):
    """-> [(name, [3 points], [3 symbols])] from the printed table."""
    out = []
    for row in block.strip().splitlines():
        name, pts, syms = [p.strip() for p in row.split(":")]
        out.append((name, parse_vectors(pts), syms.split()))
    return out


# ---------------------------------------------------------------------------
# THE DEFINITION, CODED ONCE
# ---------------------------------------------------------------------------

def admissible(a):
    """Sun's definition: a is NOT congruent to 0 mod n (some coordinate is nonzero mod 3)."""
    return any(x % N for x in a)


def viol(vs):
    """The condition Sun asks for: the sum is 0 mod n but NOT 0 mod n^2."""
    s = [sum(v[j] for v in vs) for j in range(2)]
    if any(x % N for x in s):
        return False
    return any(x % MOD for x in s)


def is_bad(S, size=N):
    """No `size` DISTINCT INDICES of S give a violating sum (repeated VALUES are allowed)."""
    return all(not viol(t) for t in itertools.combinations(S, size))


def is_bad_le(S, size=N):
    """The |I| <= n MISREADING of the definition, used only as an anti-control."""
    return all(not viol(t) for k in range(1, size + 1)
               for t in itertools.combinations(S, k))


def red(v):
    return (v[0] % N, v[1] % N)


ALPH = [t for t in itertools.product(range(MOD), repeat=2) if admissible(t)]
NONZERO = [v for v in itertools.product(range(N), repeat=2) if v != (0, 0)]
LINES = [frozenset(t) for t in itertools.combinations(NONZERO, 3)
         if all(sum(v[j] for v in t) % N == 0 for j in range(2))]

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append(bool(ok))
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, (" [%s]" % detail) if detail else ""))
    return ok


# ---------------------------------------------------------------------------
# STEP 1.  The ambient combinatorics: Fact 1 and the printed Table 1
# ---------------------------------------------------------------------------

def step1(pairs, table):
    print("=== Step 1: the ambient combinatorics of F_3^2, and Table 1 as printed")
    ck("alphabet_size_is_72", len(ALPH) == 72,
       "|{a in (Z/9)^2 : a not= 0 mod 3}| = %d = 81 - 9" % len(ALPH))

    # Fact 1: the reduction-triples that can be candidates.
    trips = [t for t in itertools.combinations_with_replacement(NONZERO, 3)
             if all(sum(v[j] for v in t) % N == 0 for j in range(2))]
    eq = [t for t in trips if len(set(t)) == 1]
    di = [t for t in trips if len(set(t)) == 3]
    mixed = [t for t in trips if len(set(t)) == 2]
    ck("fact1_no_mixed_candidate_shape", len(mixed) == 0,
       "candidate reduction-triples = %d: all-equal %d, all-distinct %d, MIXED %d (must be 0)"
       % (len(trips), len(eq), len(di), len(mixed)))
    ck("fact1_distinct_candidates_are_the_origin_missing_lines",
       {frozenset(t) for t in di} == set(LINES) and len(di) == len(LINES),
       "%d all-distinct candidates, %d origin-missing lines, same set" % (len(di), len(LINES)))
    ck("origin_missing_lines_number_eight", len(LINES) == 8, "8 of the 12 lines of AG(2,3)")
    inc = sorted(sum(1 for L in LINES if v in L) for v in NONZERO)
    ck("every_nonzero_point_lies_on_three_of_them", inc == [3] * 8,
       "incidences 8 x 3 = 24 = 8 lines x 3 points")

    pr = [set(x) for x in pairs]
    ck("printed_antipodal_pairs_are_the_four_of_F_3_2",
       len(pr) == 4 and all(len(p) == 2 for p in pr)
       and set().union(*pr) == set(NONZERO)
       and all({tuple((N - x) % N for x in v) for v in p} == p for p in pr),
       "P1..P4 partition the 8 nonzero points into v,-v pairs")
    ck("no_origin_missing_line_holds_an_antipodal_pair",
       all(not any(p <= set(L) for p in pr) for L in LINES),
       "so each line meets 3 of the 4 pairs in one point each")
    ck("each_origin_missing_line_meets_exactly_three_pairs",
       all(sum(1 for p in pr if p & set(L)) == 3 for L in LINES))

    # The printed table IS the eight lines, and its class symbols agree with the points.
    named = {}
    for name, pts, syms in table:
        named[name] = (frozenset(pts), syms)
    ck("printed_table_has_eight_rows", len(table) == 8, " ".join(n for n, _, _ in table))
    ck("printed_table_rows_are_exactly_the_eight_origin_missing_lines",
       {v[0] for v in named.values()} == set(LINES) and len(named) == len(LINES),
       "Table 1 transcription is correct: 8 distinct rows, and as a set of triples they are "
       "the 8 computed lines")
    ck("printed_table_symbols_name_the_printed_points",
       all(sorted(SYMBOL_POINT[s] for s in syms) == sorted(pts)
           for _, pts, syms in table),
       "A,A',B,B',C,C',D,D' <-> (1,0),(2,0),(0,1),(0,2),(1,1),(2,2),(1,2),(2,1)")
    ck("printed_table_uses_each_symbol_three_times",
       all(sum(1 for _, _, syms in table if s in syms) == 3 for s in SYMBOL_POINT))

    # Fact 2(a): three copies of one admissible class always violate.
    ck("fact2a_triple_of_one_class_is_always_violating",
       all(viol([a, a, a]) for a in ALPH),
       "for all 72 admissible a: 3a = 0 mod 3 and 3a not= 0 mod 9, so no class occurs 3 times")
    return named


# ---------------------------------------------------------------------------
# STEP 2.  The printed witness W: the lower half s_2(3) >= 13
# ---------------------------------------------------------------------------

def fibres(S):
    T = {}
    for v in S:
        T.setdefault(red(v), []).append(v)
    return T


def full_lines(T):
    return [L for L in LINES if set(L) <= set(T)]


def step2(W):
    print()
    print("=== Step 2: the multiset W printed in Section 5 -- the lower half")
    ck("W_has_twelve_elements", len(W) == 12, "|W| = %d" % len(W))
    ck("W_is_admissible", all(admissible(v) for v in W),
       "every element is nonzero mod 3")
    ck("W_is_BAD", is_bad(W),
       "no 3 distinct indices sum to 0 mod 3 while not 0 mod 9")
    mult = max(W.count(v) for v in set(W))
    ck("W_class_multiplicity_is_two", mult == 2,
       "max multiplicity %d, the cap of Fact 2(a)" % mult)
    F = fibres(W)
    ck("W_support_has_size_four", len(F) == 4,
       "T(W) = %s" % sorted(F))
    ck("W_support_is_two_whole_antipodal_pairs",
       sorted(F) == [(0, 1), (0, 2), (1, 0), (2, 0)],
       "P1 u P2")
    ck("W_fibres_all_have_size_three", sorted(len(x) for x in F.values()) == [3, 3, 3, 3],
       "|S| = 3|T| = 12 with equality in the master bound")
    fl = full_lines(F)
    ck("W_contains_no_full_line", len(fl) == 0,
       "full lines 0, so |A| = 0 and 3|T| - |A| = 12")
    cands = [t for t in itertools.combinations(range(12), 3)
             if all(sum(W[i][j] for i in t) % N == 0 for j in range(2))]
    sums = sorted(tuple(sum(W[i][j] for i in t) for j in range(2)) for t in cands)
    ck("W_has_exactly_four_candidates", len(cands) == 4,
       "the four same-fibre triples")
    ck("W_candidate_sums_are_the_four_printed", sums == [(0, 9), (0, 18), (9, 0), (18, 0)],
       "sums %s, every one 0 mod 9" % (sums,))
    ext = [e for e in ALPH if is_bad(list(W) + [e])]
    ck("W_is_MAXIMAL_over_all_72_admissible_extensions", len(ext) == 0,
       "bad_extensions = %d of 72" % len(ext))
    ck("lower_half_s2_3_at_least_13", len(W) == 12 and is_bad(W),
       "BAD is downward closed and s_2(3) = 1 + max|BAD|, so s_2(3) >= 13")


# ---------------------------------------------------------------------------
# STEP 3.  Directed negatives, so a "not bad" answer means something
# ---------------------------------------------------------------------------

def step3(W):
    print()
    print("=== Step 3: directed negatives -- the decider can fire")
    three = [(1, 0)] * 3
    ck("must_fire_three_copies_of_1_0_are_not_BAD", not is_bad(three),
       "sum (3,0) is 0 mod 3 and not 0 mod 9, so is_bad returns False as it must")
    ck("anti_control_W_is_not_BAD_under_the_size_at_most_3_misreading", not is_bad_le(W),
       "e.g. (1,0)+(5,0) = (6,0); |I| = n EXACTLY is load-bearing")
    ck("empty_and_singleton_multisets_are_BAD", is_bad([]) and is_bad([(1, 0)]),
       "downward closure has the expected base")


# ---------------------------------------------------------------------------
# STEP 4.  The master bound over all 256 supports
# ---------------------------------------------------------------------------

def step4():
    print()
    print("=== Step 4: the master bound 3|T| - |A(T)| over ALL 2^8 = 256 supports")
    best, unsettled = {}, []
    for k in range(9):
        best[k] = 0
    for r in range(9):
        for T in itertools.combinations(NONZERO, r):
            fl = full_lines(T)
            A = set().union(*[set(L) for L in fl]) if fl else set()
            b = 3 * len(T) - len(A)
            best[r] = max(best[r], b)
            if b > 12:
                unsettled.append(T)
    table = [best[k] for k in range(9)]
    ck("master_bound_table_by_support_size", table == [0, 3, 6, 9, 12, 12, 12, 14, 16],
       "max(3|T| - |A|) for |T| = 0..8 is %s" % (table,))
    ck("master_bound_settles_247_of_256_supports", 256 - len(unsettled) == 247,
       "settled %d, not settled %d" % (256 - len(unsettled), len(unsettled)))
    ck("the_unsettled_supports_are_exactly_the_nine_with_T_at_least_7",
       len(unsettled) == 9 and sorted(len(T) for T in unsettled) == [7] * 8 + [8],
       "C(8,7) + C(8,8) = 8 + 1 = 9; 247 + 9 = 256, fraction exhausted 1")
    ck("supports_of_size_five_and_six_are_tight_at_twelve",
       best[5] == 12 and best[6] == 12,
       "so |T| = 5 and |T| = 6 always force enough full lines")
    return unsettled


# ---------------------------------------------------------------------------
# STEP 5.  The two eliminations, as identities over all lifts
# ---------------------------------------------------------------------------

def add(*vs):
    return (sum(v[0] for v in vs) % MOD, sum(v[1] for v in vs) % MOD)


def neg(v):
    return ((-v[0]) % MOD, (-v[1]) % MOD)


def lifts(v):
    """The nine elements of (Z/9)^2 reducing to v mod 3."""
    return [((v[0] + 3 * i) % MOD, (v[1] + 3 * j) % MOD) for i in range(3) for j in range(3)]


def step5():
    print()
    print("=== Step 5: the eliminations of Section 4, as identities over all 9^3 = 729 lifts")
    trip = [(a, b, c) for a in lifts((1, 0)) for b in lifts((0, 1)) for c in lifts((1, 1))]
    ck("there_are_729_lift_triples_for_A_B_C", len(trip) == 729)

    # |T| = 8.  E1 -> C', E3 -> D, E5 -> A', E4 -> B', E2 -> D'; then E7 must give 3C.
    ok8 = True
    for A, B, C in trip:
        Cp = neg(add(A, B))                 # E1: A + B + C' = 0
        D = neg(add(A, C))                  # E3: A + C + D  = 0
        Ap = neg(add(B, D))                 # E5: A' + B + D = 0
        Bp = neg(add(Ap, C))                # E4: A' + B' + C = 0
        Dp = neg(add(A, Bp))                # E2: A + B' + D' = 0
        if red(Cp) != (2, 2) or red(D) != (1, 2) or red(Ap) != (2, 0) \
           or red(Bp) != (0, 2) or red(Dp) != (2, 1):
            ok8 = False
        if add(B, C, Dp) != ((3 * C[0]) % MOD, (3 * C[1]) % MOD):
            ok8 = False
    ck("T8_elimination_gives_E7_equal_to_3C", ok8,
       "for all 729 lifts: B + C + D' = 3C mod 9, and every derived class has the right "
       "reduction")
    ck("T8_contradiction_pi_C_is_nonzero", red((1, 1)) != (0, 0),
       "3C = 0 mod 9 forces C = 0 mod 3, but pi(C) = (1,1)")

    # |T| = 7, missing (2,1).  E1 -> C', E3 -> D, E5 -> A', E4 -> B'; then E8 gives -3(A+C).
    ok7 = True
    for A, B, C in trip:
        Cp = neg(add(A, B))
        D = neg(add(A, C))
        Ap = neg(add(B, D))
        Bp = neg(add(Ap, C))
        want = neg(((3 * (A[0] + C[0])) % MOD, (3 * (A[1] + C[1])) % MOD))
        if add(Bp, Cp, D) != want:
            ok7 = False
    ck("T7_elimination_missing_2_1_gives_E8_equal_to_minus_3_A_plus_C", ok7,
       "for all 729 lifts: B' + C' + D = -3(A + C) mod 9")
    ck("T7_contradiction_pi_A_plus_pi_C_is_nonzero",
       add((1, 0), (1, 1))[0] % N or add((1, 0), (1, 1))[1] % N,
       "pi(A) + pi(C) = (1,0) + (1,1) = (2,1) not= 0")

    # |T| = 7, missing (1,2), the paper's independent second route.
    ok7b = True
    for A, B, C in trip:
        Cp = neg(add(A, B))                 # E1
        Dp = neg(add(B, C))                 # E7: B + C + D' = 0
        Bp = neg(add(A, Dp))                # E2: A + B' + D' = 0
        Ap = neg(add(Bp, C))                # E4: A' + B' + C = 0
        want = neg(((3 * (B[0] + C[0])) % MOD, (3 * (B[1] + C[1])) % MOD))
        if add(Ap, Cp, Dp) != want:
            ok7b = False
    ck("T7_elimination_missing_1_2_gives_E6_equal_to_minus_3_B_plus_C", ok7b,
       "for all 729 lifts: A' + C' + D' = -3(B + C) mod 9, a second and different route")
    ck("T7_second_contradiction_pi_B_plus_pi_C_is_nonzero",
       add((0, 1), (1, 1))[0] % N or add((0, 1), (1, 1))[1] % N,
       "pi(B) + pi(C) = (0,1) + (1,1) = (1,2) not= 0")


# ---------------------------------------------------------------------------
# STEP 6.  STRONGER THAN THE PAPER: every support with |T| >= 7, solved exhaustively
# ---------------------------------------------------------------------------

def system_has_solution(T):
    """Is there an assignment of frozen classes x_v in (Z/9)^2, pi(x_v) = v for v in A(T),
    satisfying x_p + x_q + x_r = 0 mod 9 for every full line {p,q,r} of T?

    Writing x_v = v + 3 u_v with u_v in F_3^2 and v its canonical lift in {0,1,2}^2, a line
    equation becomes  w_L + u_p + u_q + u_r = 0 in F_3^2,  where the canonical lifts of L sum
    to 3 w_L.  The two coordinates decouple, so this is two independent linear systems over
    F_3 in |A(T)| unknowns each, and both are solved by exhaustive substitution.
    """
    fl = full_lines(T)
    if not fl:
        return True, [], set()
    A = sorted(set().union(*[set(L) for L in fl]))
    idx = {v: i for i, v in enumerate(A)}
    eqs = []
    for L in fl:
        pts = sorted(L)
        s = [sum(p[j] for p in pts) for j in range(2)]
        assert all(x % N == 0 for x in s)
        eqs.append(([idx[p] for p in pts], [(s[j] // N) % N for j in range(2)]))
    for j in range(2):
        sat = False
        for u in itertools.product(range(N), repeat=len(A)):
            if all((w[j] + sum(u[i] for i in ii)) % N == 0 for ii, w in eqs):
                sat = True
                break
        if not sat:
            return False, fl, set(A)
    return True, fl, set(A)


def step6(unsettled):
    print()
    print("=== Step 6: every support the master bound leaves open, solved mod 9 in full")
    ck("nine_supports_reach_step_6", len(unsettled) == 9,
       "the 8 supports with |T| = 7 and the 1 with |T| = 8")
    every_one_frozen_and_inconsistent = True
    for T in sorted(unsettled, key=lambda t: (-len(t), t)):
        missing = sorted(set(NONZERO) - set(T))
        tag = "T8_all_eight_points" if not missing else \
              "T7_missing_%d_%d" % missing[0]
        solvable, fl, A = system_has_solution(T)
        ck("%s_covered_by_full_lines" % tag, A == set(T),
           "full lines %d, A(T) = T (%d points), so every fibre is frozen" % (len(fl), len(A)))
        ck("%s_has_NO_solution_mod_9" % tag, not solvable,
           "the %d full-line equations are inconsistent over Z/9" % len(fl))
        if solvable or A != set(T):
            every_one_frozen_and_inconsistent = False
    ck("upper_half_s2_3_at_most_13",
       len(unsettled) == 9 and every_one_frozen_and_inconsistent,
       "|T| <= 6 gives |S| <= 12 by the master bound; and each of the %d supports the bound "
       "leaves open was just shown frozen and inconsistent mod 9; so, GRANTING Facts 2 and 3, "
       "which are proved by hand and not here, max|BAD| = 12" % len(unsettled))


# ---------------------------------------------------------------------------
# STEP 7.  The arithmetic of the source's own statements
# ---------------------------------------------------------------------------

def step7():
    print()
    print("=== Step 7: the source's own numbers")
    ck("conjectured_value_4n_plus_1_at_n_3_is_13", 4 * 3 + 1 == 13)
    ck("published_bracket_at_r_2_is_13_to_17",
       3 * 2 ** 2 + 1 == 13 and 2 * 3 ** 2 - 1 == 17,
       "3*2^r+1 <= s_r(3) <= 2*3^r-1 at r = 2")
    ck("remark_as_printed_gives_19_which_its_own_paper_refutes",
       2 * 3 ** 2 + 1 == 19 and 19 > 17,
       "2*3^r+1 at r = 2 is 19 > 17, so the printed formula is a typo for 3*2^r+1 = 13")
    ck("remark_alternative_agrees_at_r_1", 3 * 2 ** 1 + 1 == 7 and 2 * 3 ** 1 + 1 == 7,
       "both formulas give 7 at r = 1, so r = 2 is the only place they can be told apart")
    ck("sun_double_count_at_r_2_is_16", 2 * (3 ** 2 - 1) == 16,
       "2(3^r - 1) = 16, whence the published s_2(3) <= 17 that this paper improves to 13")


def main():
    print('verification of "A Proof that s_2(3) = 13, a Value Asserted without Proof by Sun"')
    print("-- inputs are the multiset W and Table 1, both printed")
    print("in the paper and transcribed into this program as text.  Exact integers mod 3 and")
    print("mod 9 only; no solver, no randomness, no network, no third-party module.")
    print()
    W = parse_vectors(W_PRINTED)
    pairs = [parse_vectors(r) for r in PAIRS_PRINTED.strip().splitlines()]
    table = parse_lines_table(LINES_PRINTED)
    step1(pairs, table)
    step2(W)
    step3(W)
    unsettled = step4()
    step5()
    step6(unsettled)
    step7()
    print()
    print("NOTE SCOPE. This program re-derives every quantity the paper claims about the "
          "objects printed in it, and it settles the two cases the master bound leaves open "
          "for EVERY support of size 7 or 8, not only for the relabelled representative the "
          "paper argues. NOT COVERED, and stated as such in the paper: (i) Sun's Lemma 2.1 "
          "(Fact 2) and Sun's collapsing claim (Fact 3) are quoted from the source and proved "
          "there and in Section 3 by hand -- only Fact 2(a)'s must-fire mechanism is machine "
          "checked here, the cancellation arguments of Fact 2(b) and Fact 3 are not; (ii) the "
          "master bound |S| <= 3|T| - |A| itself is a counting step, taken from Fact 2(b) and "
          "Fact 3 rather than re-derived; (iii) NO independent exhaustive search over the "
          "2.25 x 10^34 multisets in {0,1,2}^72 is performed here, and none is shipped in "
          "this folder, so the upper bound rests on the hand proof plus the checks above; "
          "(iv) no minimality or uniqueness of W -- the "
          "paper claims none, and other extremal multisets exist at |T| = 5 and 6; (v) n >= 4 "
          "and r >= 3, which are untouched; (vi) the line and byte locators quoted from the "
          "e-print source, which are transcribed and not fetched by this program.")
    n = len(CHECKS)
    if all(CHECKS):
        print("VERDICT: ALL %d CHECKS PASS" % n)
        return 0
    print("VERDICT: %d of %d CHECKS FAILED" % (CHECKS.count(False), n))
    return 1


if __name__ == "__main__":
    sys.exit(main())

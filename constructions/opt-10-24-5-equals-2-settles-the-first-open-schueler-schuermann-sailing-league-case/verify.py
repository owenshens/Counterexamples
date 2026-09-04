#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification program for the note

    opt(10, 24, 5) = 2:  the sailing league problem for N_teams = 10, N_inrace = 5,
    N_flights = 24, i.e. the n = 3 member of the family Schueler-Schuermann left open.

Python 3.9+, STANDARD LIBRARY ONLY (no numpy / sympy / networkx), no external data file.
All arithmetic is exact: integers and fractions.Fraction only; no float ever decides anything.

The program reads the objects PRINTED IN THE PAPER -- the 24 x 10 tournament plan (in both the
human-readable grid form and the canonical row-string form), the paper's own lambda table, and its list
of the fifteen lambda = 12 pairs -- and re-derives every quantity the paper claims about them, together
with a Schueler-Schuermann 16-flight table transcribed HERE (not printed in the paper) from the arXiv
e-print LaTeX source of arXiv:2212.02865 and used only as an outside calibration control, and with the
finite exhaustive
censuses behind Lemmas B and C and the utility-2 profile census.

Prints one `PASS <name> [detail]` line per check and closes with
    VERDICT: ALL <n> CHECKS PASS
exiting 0 iff every check passed.
"""

import sys
from fractions import Fraction
from itertools import combinations

# ----------------------------------------------------------------------------------------------------
# THE OBJECT, AS PRINTED IN THE PAPER -- both renderings, so a transcription slip cannot pass silently.
# ----------------------------------------------------------------------------------------------------

# (a) the human-readable grid, copied character-for-character from the paper's Table 1.
GRID = """
    1    |   1  1  1  1  2  2  2  2  1  2
    2    |   1  1  1  2  1  2  1  2  2  2
    3    |   1  1  1  2  2  1  2  1  2  2
    4    |   1  1  2  1  1  1  2  2  2  2
    5    |   1  2  2  2  2  2  1  1  1  1
    6    |   1  2  1  1  1  1  2  2  2  2
    7    |   1  1  2  1  2  2  1  2  2  1
    8    |   1  2  2  1  2  1  1  2  1  2
    9    |   1  2  1  2  2  1  2  2  1  1
   10    |   1  1  2  1  2  2  1  2  2  1
   11    |   1  2  2  1  2  1  1  2  1  2
   12    |   1  2  1  2  2  1  2  2  1  1
   13    |   1  1  2  1  2  2  2  1  1  2
   14    |   1  2  2  1  2  1  2  1  2  1
   15    |   1  2  1  2  2  1  1  1  2  2
   16    |   1  1  2  2  1  1  1  2  2  2
   17    |   1  2  2  2  1  2  1  2  1  1
   18    |   1  2  1  1  1  2  2  2  1  2
   19    |   1  1  2  2  1  2  2  1  1  2
   20    |   1  2  2  1  1  2  2  1  2  1
   21    |   1  2  1  1  2  2  1  1  2  2
   22    |   1  1  2  2  2  1  2  1  2  1
   23    |   1  2  2  2  1  1  1  1  2  2
   24    |   1  2  1  2  1  2  1  2  2  1
"""

# (b) the canonical machine-readable form: the 24 row strings, flight 1 first.
PLAN = [
    "1111222212", "1112121222", "1112212122", "1121112222", "1222221111", "1211112222",
    "1121221221", "1221211212", "1212212211", "1121221221", "1221211212", "1212212211",
    "1121222112", "1221212121", "1212211122", "1122111222", "1222121211", "1211122212",
    "1122122112", "1221122121", "1211221122", "1122212121", "1222111122", "1212121221",
]

# The complete lambda table exactly as printed in the paper (teams 1..10; 0 on the diagonal).
LAMBDA_TABLE = [
    [0, 10, 10, 12, 10, 12, 12, 10, 10, 10],
    [10, 0, 10, 12, 12, 10, 10, 12, 10, 10],
    [10, 10, 0, 10, 12, 12, 10, 10, 12, 10],
    [12, 12, 10, 0, 10, 10, 10, 10, 12, 10],
    [10, 12, 12, 10, 0, 10, 12, 10, 10, 10],
    [12, 10, 12, 10, 10, 0, 10, 12, 10, 10],
    [12, 10, 10, 10, 12, 10, 0, 10, 10, 12],
    [10, 12, 10, 10, 10, 12, 10, 0, 10, 12],
    [10, 10, 12, 12, 10, 10, 10, 10, 0, 12],
    [10, 10, 10, 10, 10, 10, 12, 12, 12, 0],
]

# The fifteen lambda = 12 pairs, as printed in the paper (1-based team labels).
LAMBDA12_PAIRS = [(1, 4), (1, 6), (1, 7), (2, 4), (2, 5), (2, 8), (3, 5), (3, 6), (3, 9),
                  (4, 9), (5, 7), (6, 8), (7, 10), (8, 10), (9, 10)]

# The authors' own published 16-flight plan.  NOT printed in the paper: it is transcribed HERE, from
# the arXiv e-print LaTeX source of arXiv:2212.02865 (Table tab:tournament-plan-10teams-16flights-
# 5inrace), so a referee auditing it must fetch that source.  Used ONLY as an outside calibration
# control (steps 11 and 12); no claim of the paper depends on it.  Its first 8 rows are the authors'
# (10, 8, 5) optimum.
PAPER16 = [
    "1212112122", "1122121122", "1111222212", "1121122221", "1211222121", "1221211122",
    "1122212112", "1122211221", "1221112212", "1212121212", "1222122111", "1121212212",
    "1212212211", "1221221112", "1111221222", "1122212121",
]

N_TEAMS = 10
N_INRACE = 5
PAIRS = list(combinations(range(N_TEAMS), 2))          # 45 unordered pairs, 0-based
TRIPLES = list(combinations(range(N_TEAMS), 3))        # 120 triples

# ----------------------------------------------------------------------------------------------------
# check bookkeeping
# ----------------------------------------------------------------------------------------------------

_passed = []
_failed = []


def note(msg):
    print("NOTE " + msg)


def chk(name, cond, detail=""):
    if cond:
        _passed.append(name)
        print("PASS %s%s" % (name, (" [%s]" % detail) if detail else ""))
    else:
        _failed.append(name)
        print("FAIL %s%s" % (name, (" [%s]" % detail) if detail else ""))


# ----------------------------------------------------------------------------------------------------
# basic machinery
# ----------------------------------------------------------------------------------------------------

def lam(rows):
    """lambda(t,u) = number of flights in which teams t and u sail the same race."""
    out = {}
    for (t, u) in PAIRS:
        out[(t, u)] = sum(1 for r in rows if r[t] == r[u])
    return out


def profile(lm):
    """sorted [(value, number of pairs)]."""
    h = {}
    for v in lm.values():
        h[v] = h.get(v, 0) + 1
    return sorted(h.items())


def degrees(edges):
    d = [0] * N_TEAMS
    for (a, b) in edges:
        d[a] += 1
        d[b] += 1
    return d


def is_cut(edge_set):
    """Exhibit a bipartition and VERIFY it against all 45 pairs.  -> (bool, |U|, k*(10-k))."""
    # the only candidate: U = {0} union {v : 0v is an edge}'s complement side.
    f = [0] * N_TEAMS
    for v in range(1, N_TEAMS):
        f[v] = 1 if (0, v) in edge_set else 0
    for (a, b) in PAIRS:
        if ((f[a] != f[b]) != ((a, b) in edge_set)):
            return False, None, None
    k = sum(f)
    k = min(k, N_TEAMS - k)
    return True, k, k * (N_TEAMS - k)


def is_two_cliques(edge_set):
    """The F-odd shape: the complement of a cut, i.e. a disjoint union of two cliques."""
    comp = set(p for p in PAIRS if p not in edge_set)
    return is_cut(comp)


# ====================================================================================================
print("verification of the note: opt(10,24,5) = 2 -- the sailing league problem for")
print("N_teams = 10, N_inrace = 5, N_flights = 24 (the n = 3 member left open by")
print("Schueler and Schuermann, J. Combin. Des. 32 (2024), no. 4, 171-189)")
print("python %s, exact integer / Fraction arithmetic only" % sys.version.split()[0])
print("")
print("=== Step 1: the exhibited object, read from the paper in both of its printed forms")
# ====================================================================================================

grid_rows = []
for line in GRID.strip().splitlines():
    left, right = line.split("|")
    grid_rows.append((int(left.strip()), "".join(right.split())))

chk("grid_has_24_flight_lines_numbered_1_to_24",
    [n for n, _ in grid_rows] == list(range(1, 25)),
    "flights %d..%d" % (grid_rows[0][0], grid_rows[-1][0]))
chk("printed_grid_and_printed_row_string_list_agree",
    [s for _, s in grid_rows] == PLAN,
    "24 of 24 rows identical")
chk("plan_has_24_rows_of_length_10_over_the_alphabet_1_2",
    len(PLAN) == 24 and all(len(r) == 10 and set(r) <= {"1", "2"} for r in PLAN))
chk("every_flight_is_a_five_five_split",
    all(r.count("1") == 5 and r.count("2") == 5 for r in PLAN),
    "all 24 rows have five 1s and five 2s")

note("row strings: " + " ".join(PLAN))

# ====================================================================================================
print("")
print("=== Step 2: the lambda table, the profile and the utility of the exhibited plan")
# ====================================================================================================

lm = lam(PLAN)
recomputed = [[0] * N_TEAMS for _ in range(N_TEAMS)]
for (t, u), v in lm.items():
    recomputed[t][u] = recomputed[u][t] = v

chk("all_45_entries_of_the_printed_lambda_table_are_reproduced",
    recomputed == LAMBDA_TABLE,
    "45 of 45 entries agree")

lmin, lmax = min(lm.values()), max(lm.values())
chk("lambda_min_is_10", lmin == 10, "lambda_min = %d" % lmin)
chk("lambda_max_is_12", lmax == 12, "lambda_max = %d" % lmax)
chk("utility_of_the_exhibited_plan_is_2", lmax - lmin == 2, "lambda_max - lambda_min = %d" % (lmax - lmin))

prof = profile(lm)
chk("profile_is_30_pairs_at_10_and_15_pairs_at_12", prof == [(10, 30), (12, 15)], "profile %s" % (prof,))
chk("no_pair_sits_at_lambda_11", all(v != 11 for v in lm.values()), "0 pairs at 11")

total = sum(lm.values())
chk("grand_total_is_480_equals_24_times_2_times_C_5_2",
    total == 480 == 24 * 2 * (N_INRACE * (N_INRACE - 1) // 2),
    "sum over the 45 pairs = %d" % total)
chk("profile_arithmetic_closes", 30 * 10 + 15 * 12 == total, "30*10 + 15*12 = %d" % (30 * 10 + 15 * 12))

per_team = [sum(recomputed[t][u] for u in range(N_TEAMS) if u != t) for t in range(N_TEAMS)]
chk("every_team_sum_is_96_equals_24_times_4",
    per_team == [96] * N_TEAMS, "per-team sums %s" % (sorted(set(per_team)),))

# ====================================================================================================
print("")
print("=== Step 3: the lambda = 12 graph -- cubic, and the Petersen graph")
# ====================================================================================================

E12 = set(p for p, v in lm.items() if v == 12)
printed12 = set((a - 1, b - 1) for (a, b) in LAMBDA12_PAIRS)
chk("the_fifteen_lambda_12_pairs_printed_in_the_paper_are_exactly_the_computed_ones",
    E12 == printed12, "%d edges, set equality" % len(E12))
d12 = degrees(sorted(E12))
chk("the_lambda_12_graph_is_cubic", d12 == [3] * N_TEAMS, "degree sequence %s" % (d12,))

adj12 = [set() for _ in range(N_TEAMS)]
for (a, b) in E12:
    adj12[a].add(b)
    adj12[b].add(a)
srg_lambda = set()
srg_mu = set()
for (a, b) in PAIRS:
    common = len(adj12[a] & adj12[b])
    (srg_lambda if b in adj12[a] else srg_mu).add(common)
chk("the_lambda_12_graph_is_strongly_regular_with_parameters_10_3_0_1",
    srg_lambda == {0} and srg_mu == {1},
    "adjacent pairs share %s common neighbours, non-adjacent pairs share %s"
    % (sorted(srg_lambda), sorted(srg_mu)))

# explicit isomorphism to the Kneser graph K(5,2)
KV = [frozenset(s) for s in combinations(range(1, 6), 2)]
kadj = [set() for _ in range(10)]
for i in range(10):
    for j in range(10):
        if i != j and not (KV[i] & KV[j]):
            kadj[i].add(j)


def find_iso(pos, used, phi):
    if pos == N_TEAMS:
        return list(phi)
    for cand in range(10):
        if cand in used:
            continue
        ok = True
        for q in range(pos):
            if (q in adj12[pos]) != (phi[q] in kadj[cand]):
                ok = False
                break
        if ok:
            phi.append(cand)
            used.add(cand)
            got = find_iso(pos + 1, used, phi)
            if got is not None:
                return got
            used.discard(cand)
            phi.pop()
    return None


iso = find_iso(0, set(), [])
chk("an_explicit_bijection_to_the_kneser_graph_K_5_2_petersen_exists",
    iso is not None and all(((b in adj12[a]) == (iso[b] in kadj[iso[a]])) for (a, b) in PAIRS),
    "team t -> 2-subset: " + ", ".join("%d:%s" % (t + 1, sorted(KV[iso[t]])) for t in range(N_TEAMS))
    if iso else "no bijection found")

chk("design_parameters_of_the_exhibited_plan_are_v_10_k_5_r_24_b_48",
    (N_TEAMS, N_INRACE, len(PLAN), 2 * len(PLAN)) == (10, 5, 24, 48),
    "v=10 k=5 r=24 b=48; lambda_1 = 10 on intersecting 2-subsets, lambda_2 = 12 on disjoint ones")

# ====================================================================================================
print("")
print("=== Step 4: the repeated flights, disclosed in the paper, and the periodicity checks")
# ====================================================================================================

distinct = sorted(set(PLAN))
chk("the_plan_uses_21_distinct_flights_of_24", len(distinct) == 21,
    "%d distinct rows" % len(distinct))
repeats = sorted(tuple(sorted(i + 1 for i, r in enumerate(PLAN) if r == s))
                 for s in distinct if PLAN.count(s) > 1)
chk("the_repeated_flights_are_exactly_7_10_and_8_11_and_9_12",
    repeats == [(7, 10), (8, 11), (9, 12)], "repeat classes %s" % (repeats,))
chk("the_plan_is_not_12_periodic", PLAN[:12] != PLAN[12:],
    "rows 1-12 differ from rows 13-24, so the plan is not a doubled 12-flight plan")
chk("the_plan_is_not_8_periodic", not (PLAN[:8] == PLAN[8:16] == PLAN[16:]),
    "rows 1-8, 9-16, 17-24 are not all equal, so the plan is not three copies of one base schedule")

# ====================================================================================================
print("")
print("=== Step 5: Lemma A -- the authors' average at this cell, and utility 0")
# ====================================================================================================

F = 24
lam_avg = Fraction(F * (N_INRACE - 1), N_TEAMS - 1)
chk("lambda_average_is_32_over_3_by_the_authors_lemma",
    lam_avg == Fraction(96, 9) == Fraction(32, 3),
    "lambda = 24*4/9 = %s" % lam_avg)
chk("lambda_average_is_not_an_integer_so_utility_0_is_impossible",
    lam_avg.denominator != 1, "denominator %d, so no perfect pairing list exists at 24 flights"
    % lam_avg.denominator)
chk("the_per_team_sum_and_the_grand_total_agree_with_the_average",
    Fraction(sum(per_team), N_TEAMS * (N_TEAMS - 1)) == lam_avg,
    "(sum of the ten per-team sums) / (10*9) = %d/90 = %s"
    % (sum(per_team), Fraction(sum(per_team), 90)))

# ====================================================================================================
print("")
print("=== Step 6: Lemma B exhaustively -- every triple, in every five-five split")
# ====================================================================================================

splits = [frozenset((0,) + c) for c in combinations(range(1, N_TEAMS), N_INRACE - 1)]
chk("there_are_126_five_five_splits_of_10_teams", len(splits) == 126, "%d splits" % len(splits))
viol = 0
nchecks = 0
for S in splits:
    for (a, b, c) in TRIPLES:
        same = ((a in S) == (b in S)) + ((a in S) == (c in S)) + ((b in S) == (c in S))
        nchecks += 1
        if same % 2 == 0:
            viol += 1
chk("lemma_B_holds_on_every_split_and_every_triple",
    viol == 0 and nchecks == 15120,
    "%d checks (126 splits x 120 triples), %d violations" % (nchecks, viol))
chk("hence_every_triangle_sum_is_congruent_to_N_flights_mod_2",
    all((lm[(a, b)] + lm[(a, c)] + lm[(b, c)]) % 2 == F % 2 for (a, b, c) in TRIPLES),
    "verified on the exhibited plan too: all 120 triangle sums even (N_flights = 24)")

# ====================================================================================================
print("")
print("=== Step 7: Lemma C exhaustively -- every graph on 10 labelled vertices with all triples even")
# ====================================================================================================
# COMPLETENESS OF THE ENUMERATION: if every triple spans an even number of edges then, for each triple
# (0,u,v), e(u,v) = e(0,u) + e(0,v) mod 2.  So such a graph is DETERMINED by the 9 bits e(0,v).  The 512
# candidates below therefore exhaust the possibilities; each is then tested on all 120 triples.

cuts = []
for mask in range(512):
    e = {}
    for v in range(1, N_TEAMS):
        e[(0, v)] = (mask >> (v - 1)) & 1
    for u in range(1, N_TEAMS):
        for v in range(u + 1, N_TEAMS):
            e[(u, v)] = (e[(0, u)] + e[(0, v)]) % 2
    if all((e[(a, b)] + e[(a, c)] + e[(b, c)]) % 2 == 0 for (a, b, c) in TRIPLES):
        cuts.append(frozenset(p for p in PAIRS if e[p]))

chk("the_number_of_graphs_on_10_vertices_with_all_120_triples_even_is_512",
    len(cuts) == 512 and len(set(cuts)) == 512, "%d graphs, all distinct" % len(cuts))

not_a_cut = 0
edge_counts = set()
degseqs = set()
six_regular = 0
side_sizes = set()
for g in cuts:
    ok, k, kk = is_cut(g)
    if not ok:
        not_a_cut += 1
        continue
    side_sizes.add(k)
    edge_counts.add(len(g))
    d = tuple(sorted(degrees(sorted(g))))
    degseqs.add(d)
    if len(set(d)) == 1 and d[0] == 6:
        six_regular += 1
chk("every_one_of_them_is_a_complete_bipartite_cut",
    not_a_cut == 0, "0 of 512 fail; each bipartition exhibited and verified against all 45 pairs")
chk("the_realised_edge_counts_are_exactly_0_9_16_21_24_25",
    sorted(edge_counts) == [0, 9, 16, 21, 24, 25], "k(10-k) for k in %s -> %s"
    % (sorted(side_sizes), sorted(edge_counts)))
chk("the_realised_degree_sequences_are_exactly_the_six_listed_in_the_paper",
    sorted(degseqs) == sorted([tuple([0] * 10), tuple([1] * 9 + [9]), tuple([2] * 8 + [8] * 2),
                               tuple([3] * 7 + [7] * 3), tuple([4] * 6 + [6] * 4), tuple([5] * 10)]),
    "0^10, 1^9 9, 2^8 8^2, 3^7 7^3, 4^6 6^4, 5^10")
chk("none_of_the_512_is_6_regular", six_regular == 0,
    "%d six-regular graphs among all graphs with all triples even" % six_regular)

REG_DEGREES = set(d[0] for d in degseqs if len(set(d)) == 1)
chk("the_only_regular_all_even_triple_graphs_have_degree_0_or_5",
    REG_DEGREES == {0, 5}, "regular degrees realised: %s" % sorted(REG_DEGREES))

# ====================================================================================================
print("")
print("=== Step 8: Theorem 1 -- the lower bound opt(10,24,5) >= 2, with no assumption left implicit")
# ====================================================================================================


def utility_at_most_one_windows(flights):
    """All windows {m, m+1} compatible with the per-team sum 4*flights over 9 opponents."""
    S = 4 * flights
    return [m for m in range(0, S + 1) if 9 * m <= S <= 9 * m + 9]


def utility_one_forced_odd_graph(flights, m):
    """-> (odd value, forced regular degree of the odd-lambda graph) for the window {m, m+1}."""
    S = 4 * flights
    a = S - 9 * m                    # number of neighbours at value m+1, per team
    assert 0 <= a <= 9
    if m % 2 == 1:                   # m odd, m+1 even -> the odd graph is the m-graph
        return m, 9 - a
    return m + 1, a                  # m even, m+1 odd -> the odd graph is the (m+1)-graph


wins24 = utility_at_most_one_windows(24)
chk("at_24_flights_utility_at_most_1_forces_the_window_10_11",
    wins24 == [10], "9m <= 96 <= 9m+9 has the unique solution m = %d" % wins24[0])
oddval, deg = utility_one_forced_odd_graph(24, 10)
chk("that_window_forces_the_odd_lambda_graph_to_be_6_regular",
    (oddval, deg) == (11, 6),
    "11 is the only odd value in {10,11}; 90 + a(t) = 96 gives a(t) = 6 for every team")
chk("no_6_regular_graph_has_all_triples_even_so_utility_1_is_impossible_at_24_flights",
    deg not in REG_DEGREES and (deg * 10 // 2) not in edge_counts,
    "6-regular needs 30 edges; 30 is not of the form k(10-k) and 6 is not a realised regular degree")
chk("theorem_1_opt_10_24_5_is_at_least_2",
    lam_avg.denominator != 1 and deg not in REG_DEGREES,
    "utility 0 killed by Lemma A, utility 1 killed by Lemmas B and C")
chk("combined_with_the_exhibited_plan_opt_10_24_5_equals_2",
    (lmax - lmin) == 2 and lam_avg.denominator != 1 and deg not in REG_DEGREES,
    "upper bound 2 from Step 2, lower bound 2 from Theorem 1")

# ====================================================================================================
print("")
print("=== Step 9: the utility-2 profile census at 24 flights")
# ====================================================================================================
# For a window {m, m+1, m+2}: 9m <= 96 <= 9m+18, so m in {9, 10}.  For each candidate ODD graph (which
# by Lemmas B and C must be one of the 512 cuts) and each team, ask whether nonnegative counts
# (x0, x1, x2) of neighbours at m, m+1, m+2 exist with x0+x1+x2 = 9, the per-team sum equal to 96, and
# the odd-lambda degree equal to that team's degree in the candidate graph.

wins2 = [m for m in range(0, 97) if 9 * m <= 96 <= 9 * m + 18]
chk("at_24_flights_utility_2_allows_only_the_windows_9_10_11_and_10_11_12",
    wins2 == [9, 10], "9m <= 96 <= 9m+18 gives m in %s" % wins2)


_OPT_CACHE = {}


def team_options(m, d):
    """All (x0, x1, x2) for one team, given its degree d in the candidate odd-lambda graph."""
    if (m, d) in _OPT_CACHE:
        return _OPT_CACHE[(m, d)]
    out = []
    for x2 in range(10):
        for x1 in range(10 - x2):
            x0 = 9 - x1 - x2
            if x0 < 0:
                continue
            if m * x0 + (m + 1) * x1 + (m + 2) * x2 != 96:
                continue
            odd_deg = x1 if m % 2 == 0 else x0 + x2
            if odd_deg == d:
                out.append((x0, x1, x2))
    _OPT_CACHE[(m, d)] = out
    return out


survivors = {9: [], 10: []}
for m in (9, 10):
    for g in cuts:
        dd = degrees(sorted(g))
        opts = [team_options(m, dd[t]) for t in range(N_TEAMS)]
        if all(opts):
            survivors[m].append((g, opts))

chk("the_window_9_10_11_is_impossible", len(survivors[9]) == 0,
    "0 of 512 candidate odd graphs admit per-team counts summing to 96")

found = {}
for (g, opts) in survivors[10]:
    counts = {10: 0, 11: 0, 12: 0}
    consistent = all(len(o) == 1 for o in opts)
    if not consistent:
        continue
    for t in range(N_TEAMS):
        x0, x1, x2 = opts[t][0]
        counts[10] += x0
        counts[11] += x1
        counts[12] += x2
    key = tuple(sorted((v, c // 2) for v, c in counts.items() if c))
    found[key] = found.get(key, 0) + 1

chk("every_surviving_case_pins_its_per_team_counts_uniquely",
    len(found) and sum(found.values()) == len(survivors[10]),
    "%d surviving odd graphs, all with a unique per-team solution" % len(survivors[10]))
P1 = ((10, 30), (12, 15))
P2 = ((10, 18), (11, 24), (12, 3))
chk("exactly_two_utility_2_profiles_survive_up_to_relabelling",
    set(found) == {P1, P2}, "profiles %s" % sorted(found))
chk("profile_P1_arithmetic_closes_at_480", 30 * 10 + 15 * 12 == 480, "30*10 + 15*12 = 480")
chk("profile_P2_arithmetic_closes_at_480", 18 * 10 + 24 * 11 + 3 * 12 == 480,
    "18*10 + 24*11 + 3*12 = 480")
chk("the_surviving_odd_graphs_are_the_empty_graph_and_the_210_four_six_cuts",
    sorted(len(g) for (g, _) in survivors[10]) == [0] + [24] * 210,
    "1 empty graph (profile P1) and 210 cuts K_{4,6} (profile P2)")
chk("the_exhibited_plan_realises_profile_P1_with_a_cubic_12_graph",
    tuple(prof) == P1 and d12 == [3] * N_TEAMS, "profile %s, 12-graph cubic" % (prof,))
chk("the_window_5_at_which_the_odd_graph_would_be_K_5_5_is_excluded",
    all(len(g) != 25 for (g, _) in survivors[10]),
    "k = 5 (odd graph K_{5,5}, degree sequence 5^10) survives the crude per-team bound but not "
    "d(t) = 6 - 2c(t), which forces d(t) even")

# ====================================================================================================
print("")
print("=== Step 10: by-products of the same route at other even flight counts")
# ====================================================================================================


def lower_bound_two(flights):
    """-> (verdict string, detail).  Decides utility 0 and utility 1 at an EVEN flight count."""
    avg = Fraction(flights * (N_INRACE - 1), N_TEAMS - 1)
    zero_dead = avg.denominator != 1
    ws = utility_at_most_one_windows(flights)
    details = []
    one_dead = True
    for m in ws:
        ov, dg = utility_one_forced_odd_graph(flights, m)
        alive = (dg in REG_DEGREES) and (dg * 10 // 2 in edge_counts)
        details.append("window {%d,%d}: odd graph is the %d-graph, forced %d-regular (%d edges) -> %s"
                       % (m, m + 1, ov, dg, dg * 10 // 2, "POSSIBLE" if alive else "impossible"))
        if alive:
            one_dead = False
    return (zero_dead and one_dead), "lambda = %s; %s" % (avg, "; ".join(details))


for flights, label in ((12, "opt_10_12_5"), (16, "opt_10_16_5"), (32, "opt_10_32_5")):
    ok, detail = lower_bound_two(flights)
    chk("byproduct_%s_is_at_least_2" % label, ok, detail)

ok8, detail8 = lower_bound_two(8)
chk("control_the_same_route_at_8_flights_gives_only_at_least_2_and_never_contradicts_the_published_3",
    ok8, detail8 + " -- the authors prove opt(10,8,5) = 3, so >= 2 is weaker and consistent")

# ====================================================================================================
print("")
print("=== Step 11: calibration controls on the authors' OWN published 16-flight table")
# ====================================================================================================

chk("the_published_16_flight_table_is_a_legal_pairing_list",
    all(len(r) == 10 and r.count("1") == 5 and r.count("2") == 5 for r in PAPER16),
    "16 rows, each a five-five split")

lm16 = lam(PAPER16)
p16 = profile(lm16)
u16 = max(lm16.values()) - min(lm16.values())
pt16 = [sum(lm16[tuple(sorted((t, u)))] for u in range(N_TEAMS) if u != t) for t in range(N_TEAMS)]
chk("control_A_the_authors_16_flight_theorem_is_reproduced",
    (min(lm16.values()), max(lm16.values()), u16) == (6, 8, 2) and p16 == [(6, 20), (8, 25)]
    and sum(lm16.values()) == 320 and pt16 == [64] * 10,
    "lambda_min=6 lambda_max=8 UTILITY=2, profile %s, total 320, per-team sums 64 -- their theorem "
    "states utility 2" % (p16,))

lm8 = lam(PAPER16[:8])
p8 = profile(lm8)
u8 = max(lm8.values()) - min(lm8.values())
pt8 = [sum(lm8[tuple(sorted((t, u)))] for u in range(N_TEAMS) if u != t) for t in range(N_TEAMS)]
chk("control_B_the_authors_8_flight_theorem_is_reproduced_including_its_lambda_min_and_lambda_max",
    (min(lm8.values()), max(lm8.values()), u8) == (2, 5, 3) and sum(lm8.values()) == 160
    and pt8 == [32] * 10,
    "lambda_min=2 lambda_max=5 UTILITY=3, profile %s, total 160, per-team sums 32 -- their "
    "Theorem thm:asian_pacific states lambda_min=2, lambda_max=5, utility 3" % (p8,))

# ====================================================================================================
print("")
print("=== Step 12: the NON-VACUOUS parity test -- truncations of the authors' table that carry ODD")
print("             lambda values, so that 0 violations is not guaranteed in advance")
# ====================================================================================================
# On an object all of whose lambda values are even, 'every triangle sum is even' is true before any
# computation and is therefore no evidence at all.  Truncating the published table to F rows produces
# objects with odd lambda values, and Lemmas B and C then make a FALSIFIABLE prediction: the odd-lambda
# graph is a complete bipartite cut when F is even, and the complement of one (two cliques) when F is
# odd.

for Ftrunc in (7, 8, 9, 11, 13, 16):
    lmt = lam(PAPER16[:Ftrunc])
    odd = set(p for p, v in lmt.items() if v % 2 == 1)
    if Ftrunc % 2 == 0:
        ok, k, kk = is_cut(odd)
        shape = "a complete bipartite cut"
    else:
        ok, k, kk = is_two_cliques(odd)
        shape = "the complement of a cut (two cliques)"
    chk("parity_prediction_holds_on_the_published_table_truncated_to_%d_flights" % Ftrunc,
        ok, "%d pairs have odd lambda; predicted shape = %s; realised with |U| = %s and k(10-k) = %s"
        % (len(odd), shape, k, kk))

# ====================================================================================================
print("")
print("=== Scope")
# ====================================================================================================
note("SCOPE: this program verifies the single cell (N_teams, N_flights, N_inrace) = (10, 24, 5), i.e. "
     "the n = 3 member of the family the source leaves open. The family N_flights = 8n with n >= 3 is "
     "NOT closed: n = 4 (32 flights), n = 5 (40 flights) and every larger n remain open, and only the "
     "bound opt(10,32,5) >= 2 is established above -- no 32-flight witness exists in this record.")
note("SCOPE: the exhibited plan uses 21 distinct flights, three of them twice (7 = 10, 8 = 11, "
     "9 = 12). That is legal because the source defines a pairing LIST and says explicitly that races "
     "and even flights may repeat. The distinct-flights variant of this cell is a DIFFERENT question "
     "and is not answered here, nor is its optimum.")
note("NOT RE-RUN: the search that FOUND the plan. The plan was produced by a CP-SAT run whose captured "
     "stdout is not part of the record; this program does not re-run any solver, needs none, and takes "
     "the plan purely as a printed object. Nothing above depends on the search.")
note("NOT RE-RUN: profile P2 is shown to be the only other arithmetically admissible utility-2 "
     "profile; whether any (10, 24, 5) plan REALISES P2 is unknown and irrelevant to opt(10,24,5) = 2. "
     "Likewise the exact values of opt(10,12,5) and opt(10,16,5) are not claimed here -- only the "
     "bound >= 2 in each case (opt(10,16,5) = 2 is the authors' own theorem).")
note("NOT RE-RUN: the bibliographic and prior-art record (the journal reference, the DOI, the citer "
     "registries and the design catalogues). Those are reads, not computations, and are not checked here.")

print("")
if _failed:
    print("VERDICT: %d CHECK(S) FAILED of %d: %s" % (len(_failed), len(_failed) + len(_passed),
                                                     ", ".join(_failed)))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % len(_passed))
sys.exit(0)

#!/usr/bin/env python3
"""Verification of a ten-bamboo counterexample to the Reduce-Max backlog-2 conjecture.

Setting.  A bamboo instance has positive growth rates r_1..r_n with sum 1 and all
heights initially 0.  Each round every bamboo grows by its rate, then one bamboo of
maximum height is cut to 0.  The backlog is the supremum over rounds of the maximum
pre-cut height.  The conjecture under attack asserts backlog <= 2 after this
normalization.

VALUES TAKEN FROM THE PAPER (inputs, transcribed and never used as evidence for
themselves; each is re-derived below from the rate vector alone):
  * the rate vector r = w/937 with w = (313,261,131,118,19,19,19,19,19,19);
  * the claimed backlog 1878/937 = 2 + 4/937, attained in round 60;
  * the certificate table: eight round intervals with their cut words and their
    m / M / Delta columns;
  * the nine canonical post-cut states at the interval endpoints;
  * the claim that the only ties for the global maximum are rounds 43,44,51,52,55;
  * the five decisive comparisons for rounds 56-60;
  * seven-periodicity of the canonical trajectory from the post-round-137 state,
    with cycle word 214321S and cycle peak 1252/937;
  * the corollary construction for n = 10+m: rates scaled by 999/1000 plus m
    bamboos of rate 1/(1000m), with round-60 height 999/1000 * 1878/937
    = 2 + 2122/937000;
  * the remark that the search supporting the conjecture enumerated integer
    partitions of the total growth H in {5,10,...,35}.

DERIVED HERE (the checks; all integer or exact-rational, no floats):
  1. the exhibited instance is well formed: 10 positive weights, sum 937, gcd 1,
     and the rates sum to exactly 1;
  2. the greedy cut word over rounds 1..144 recomputed from w equals the
     concatenation of the table's eight words, with the stated lengths;
  3. the m / M / Delta columns recomputed per interval;
  4. the nine canonical endpoint states recomputed;
  5. the load-bearing violation: the maximum pre-cut height over rounds 1..144 is
     computed to be 1878/937, attained only in round 60, and 1878/937 > 2;
  6. the five decisive comparisons, including that the multipliers really are the
     rounds since the most recent cut and that the two sides really are the cut
     height and the best competing rate class;
  7. tie structure: over the whole simulated horizon the rounds with a tied global
     maximum are exactly 43,44,51,52,55 -- so the periodic tail is tie-free -- and
     every such tie is between equal-rate, equal-height slow bamboos;
  8. tie-independence by EXHAUSTIVE enumeration of every tie-breaking rule for
     rounds 1..144 (branching on all maximizers), showing one canonical state per
     round and peak 1878 in round 60 under all of them;
  9. seven-periodicity from round 137, cycle word, cycle peak 1252, and that the
     supremum over a long horizon is still 1878;
 10. the corollary: exact rational simulation of the (10+m)-bamboo instances for
     m = 1..12, verifying positivity, sum 1, that no added bamboo is cut through
     round 60, the growth bound used in the proof, and a round-60 height above 2;
 11. a re-run of the conjecture-supporting census: all integer partitions of
     H in {5,10,...,35}, whose worst backlog under lowest-index tie-breaking is
     below 2, and that H = 937 is outside that range.

NOT RE-RUN (stated for the record, and printed by the program): the 2000-bamboo
and 2702-bamboo instances of the prior disproof, and their reported lower bounds
2.0004 and 2.076, are not reconstructed here -- those instances are not exhibited
in the paper being checked; the corollary is simulated exactly only for
n = 10..22; nothing is claimed or tested about fewer than ten bamboos; and in the
H <= 35 census the backlog is exact only where the state recurred inside the
4000-round horizon and is computed under lowest-index tie-breaking only.  That
census merely corroborates a remark about an earlier search; no part of the main
theorem rests on it, and the theorem's tie-independence is instead established by
exhaustive expansion of every tie-breaking.
"""
from fractions import Fraction
from math import gcd

# ---------------------------------------------------------------- paper inputs
W = (313, 261, 131, 118, 19, 19, 19, 19, 19, 19)
P_SUM = 937
P_PEAK = 1878
P_PEAK_ROUND = 60
P_HORIZON = 144
P_CYCLE_START = 137
P_PERIOD = 7
P_CYCLE_WORD = "214321S"
P_CYCLE_PEAK = 1252
P_TABLE = (
    (1, 6, 6, "121321", 313, 939, 2),
    (7, 42, 36, "421321" * 6, 708, 939, 4),
    (43, 60, 18, "SS214321SS21SS4321", 817, 1878, 1),
    (61, 72, 12, "214321" * 2, 522, 1252, 3),
    (73, 100, 28, "S214321" * 4, 570, 1252, 23),
    (101, 108, 8, "SS214321", 874, 1252, 27),
    (109, 137, 29, "S214321" * 4 + "S", 684, 1252, 40),
    (138, 144, 7, "214321S", 783, 1252, 40),
)
P_ENDPOINTS = {
    0: (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    6: (0, 261, 262, 708, 114, 114, 114, 114, 114, 114),
    42: (0, 261, 262, 590, 798, 798, 798, 798, 798, 798),
    60: (0, 261, 262, 354, 76, 95, 152, 171, 304, 323),
    72: (0, 261, 262, 354, 304, 323, 380, 399, 532, 551),
    100: (0, 261, 262, 354, 114, 247, 380, 513, 836, 855),
    108: (0, 261, 262, 354, 114, 133, 266, 399, 532, 665),
    137: (313, 522, 393, 472, 0, 133, 266, 399, 532, 665),
    144: (313, 522, 393, 472, 0, 133, 266, 399, 532, 665),
}
P_TIE_ROUNDS = (43, 44, 51, 52, 55)
# (round, multiplier_cut, rate_cut, multiplier_rival, rate_rival)
P_DECISIVE = (
    (56, 56, 19, 9, 118),
    (57, 10, 118, 9, 131),
    (58, 10, 131, 5, 261),
    (59, 6, 261, 5, 313),
)
P_ROUND60 = (60, 6, 313)
P_DDN_H = tuple(range(5, 36, 5))
P_COR_ROUND60 = Fraction(2) + Fraction(2122, 937000)

CHECKS = []


def canon(state):
    """Keep the four fast coordinates, sort the six equal-rate slow ones."""
    return tuple(state[:4]) + tuple(sorted(state[4:]))


def symbol(j):
    """Cut-word symbol: digit for a singleton rate class, S for the slow class."""
    return str(j + 1) if j < 4 else "S"


def simulate(w, rounds, zero=0):
    """Reduce-Max with lowest-index tie-breaking.

    Returns a list of rounds 1..rounds; entry t-1 is a dict with the pre-cut
    height vector, the index cut, its pre-cut height, the number of maximizers,
    and the post-cut state.  Arithmetic is whatever type `zero` and w are.
    """
    n = len(w)
    h = [zero] * n
    out = []
    for _ in range(rounds):
        pre = [h[i] + w[i] for i in range(n)]
        top = max(pre)
        maximizers = [i for i in range(n) if pre[i] == top]
        j = maximizers[0]
        h = list(pre)
        h[j] = zero
        out.append({
            "pre": tuple(pre),
            "cut": j,
            "height": top,
            "nmax": len(maximizers),
            "post": tuple(h),
        })
    return out


def cross_rate_gap(w, pre, j):
    """pre-cut height of the cut bamboo minus the best pre-cut height in any
    other rate class (the paper's delta_t)."""
    rivals = [pre[i] for i in range(len(w)) if w[i] != w[j]]
    return pre[j] - max(rivals)


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)


def check_instance():
    """Hypotheses of the statement: 10 bamboos, positive rates summing to 1."""
    n = len(W)
    total = sum(W)
    g = 0
    for x in W:
        g = gcd(g, x)
    rates = [Fraction(x, total) for x in W]
    ok = (n == 10 and total == P_SUM and all(x > 0 for x in W)
          and sum(rates) == 1 and all(x > 0 for x in rates) and g == 1)
    ck("instance_well_formed", ok,
       "n=%d, sum(w)=%d, gcd=%d, sum(rates)=%s, min rate=%s"
       % (n, total, g, sum(rates), min(rates)))
    ck("rate_classes_as_described",
       len(set(W[:4])) == 4 and len(set(W[4:])) == 1 and len(W[4:]) == 6
       and all(W[i] not in W[4:] for i in range(4)),
       "four singleton classes 313,261,131,118 and one class of six 19s")


def check_cut_word(run):
    """The recomputed greedy cut word equals the certificate's eight words."""
    word = "".join(symbol(r["cut"]) for r in run[:P_HORIZON])
    lengths_ok = all(b - a + 1 == ln == len(wd) for a, b, ln, wd, _, _, _ in P_TABLE)
    covers = [(a, b) for a, b, _, _, _, _, _ in P_TABLE]
    contiguous = covers[0][0] == 1 and covers[-1][1] == P_HORIZON and all(
        covers[i + 1][0] == covers[i][1] + 1 for i in range(len(covers) - 1))
    joined = "".join(wd for _, _, _, wd, _, _, _ in P_TABLE)
    ck("certificate_cut_word_matches_greedy", word == joined and lengths_ok
       and contiguous and len(joined) == P_HORIZON,
       "144 symbols agree; first mismatch: %s"
       % next((str(i + 1) for i in range(min(len(word), len(joined)))
               if word[i] != joined[i]), "none"))


def check_table_columns(run):
    """The m, M and Delta columns recomputed interval by interval."""
    bad = []
    for a, b, _, _, m, M, delta in P_TABLE:
        seg = run[a - 1:b]
        sel = [r["height"] for r in seg]
        gaps = [cross_rate_gap(W, r["pre"], r["cut"]) for r in seg]
        if (min(sel), max(sel), min(gaps)) != (m, M, delta):
            bad.append("%d-%d got (%d,%d,%d)"
                       % (a, b, min(sel), max(sel), min(gaps)))
    ck("certificate_m_M_delta_columns", not bad, "; ".join(bad) or "8/8 rows exact")
    # The recurrence is replayed here from scratch, so that "the cut attains the
    # global maximum" is tested against an independently rebuilt pre-cut vector
    # rather than against the value the simulator happened to store.
    replay_ok, prev, worst_gap = True, (0,) * len(W), None
    for r in run[:P_HORIZON]:
        pre = tuple(prev[i] + W[i] for i in range(len(W)))
        top = max(pre)
        j = r["cut"]
        post = tuple(0 if i == j else pre[i] for i in range(len(W)))
        gap = cross_rate_gap(W, pre, j)
        if worst_gap is None or gap < worst_gap:
            worst_gap = gap
        if not (pre == r["pre"] and post == r["post"] and pre[j] == top
                and r["height"] == top and gap > 0):
            replay_ok = False
        prev = post
    ck("all_cuts_globally_greedy", replay_ok and worst_gap > 0,
       "recurrence replayed for rounds 1..%d; every cut agrees with the "
       "independently rebuilt pre-cut vector and its maximum: %s; smallest "
       "cross-rate gap: %s" % (P_HORIZON, replay_ok, worst_gap))


def check_endpoints(run):
    """The nine canonical post-cut states at the interval endpoints."""
    bad = []
    for t, claimed in sorted(P_ENDPOINTS.items()):
        raw = (0,) * 10 if t == 0 else run[t - 1]["post"]
        # the paper's states are already canonical: last six nondecreasing
        sorted_tail = list(claimed[4:]) == sorted(claimed[4:])
        # fast coordinates compared in place, slow class as a multiset
        if not (sorted_tail and raw[:4] == claimed[:4]
                and sorted(raw[4:]) == sorted(claimed[4:])
                and canon(raw) == claimed):
            bad.append("t=%d got %s" % (t, canon(raw)))
    ck("canonical_endpoint_states", not bad, "; ".join(bad) or "9/9 states exact")


def check_violation(run):
    """LOAD-BEARING: the greedy backlog is computed, then compared with 2."""
    seg = run[:P_HORIZON]
    peak = max(r["height"] for r in seg)
    where = [i + 1 for i, r in enumerate(seg) if r["height"] == peak]
    backlog = Fraction(peak, P_SUM)
    ck("peak_pre_cut_height_is_1878_in_round_60",
       peak == P_PEAK and where == [P_PEAK_ROUND],
       "computed peak %d attained in round(s) %s" % (peak, where))
    ck("backlog_exceeds_two",
       backlog > 2 and backlog == Fraction(2) + Fraction(4, P_SUM)
       and backlog == Fraction(P_PEAK, P_SUM),
       "computed backlog %s = 2 %s %s; compared with 2 in exact rationals"
       % (backlog, "+" if backlog >= 2 else "-", abs(backlog - 2)))
    r60 = run[P_PEAK_ROUND - 1]
    ck("round_60_cuts_the_fastest_bamboo",
       r60["cut"] == 0 and r60["height"] == P_PEAK and r60["nmax"] == 1,
       "round 60 cuts bamboo 1 at height %d, uniquely maximal" % r60["height"])


def rounds_since_cut(run, t, i):
    """How many rounds bamboo i has grown when round t's pre-cut height forms."""
    for s in range(t - 1, 0, -1):
        if run[s - 1]["cut"] == i:
            return t - s
    return t


def check_decisive(run):
    """Rounds 56-60: each side of the paper's comparison is recomputed."""
    bad = []
    for t, mult, rate, rmult, rrate in P_DECISIVE:
        r = run[t - 1]
        j = r["cut"]
        lhs, rhs = mult * rate, rmult * rrate
        rivals = [(r["pre"][i], i) for i in range(len(W)) if W[i] != W[j]]
        best, bi = max(rivals)
        ok = (W[j] == rate and r["height"] == lhs and lhs > rhs
              and rounds_since_cut(run, t, j) == mult and best == rhs
              and W[bi] == rrate and rounds_since_cut(run, t, bi) == rmult)
        if not ok:
            bad.append("round %d: cut w=%d h=%d age=%d, rival h=%d w=%d"
                       % (t, W[j], r["height"], rounds_since_cut(run, t, j),
                          best, W[bi]))
    t, mult, rate = P_ROUND60
    r = run[t - 1]
    if not (W[r["cut"]] == rate and rounds_since_cut(run, t, r["cut"]) == mult
            and r["height"] == mult * rate == 2 * P_SUM + 4):
        bad.append("round 60 arithmetic")
    ck("decisive_comparisons_rounds_56_to_60", not bad,
       "; ".join(bad) or "5/5 comparisons exact, multipliers = rounds since last cut")


def check_ties(run):
    """The only ties for the global maximum, and their nature."""
    # Scanned over the whole simulated horizon, not merely rounds 1..144: the
    # claim "every tie-breaking" needs the later rounds to be tie-free too.
    tied = tuple(i + 1 for i, r in enumerate(run) if r["nmax"] > 1)
    slow_only = True
    for t in tied:
        r = run[t - 1]
        top = r["height"]
        idx = [i for i in range(len(W)) if r["pre"][i] == top]
        if not all(i >= 4 for i in idx) or len({W[i] for i in idx}) != 1:
            slow_only = False
    ck("tie_rounds_are_exactly_43_44_51_52_55", tied == P_TIE_ROUNDS,
       "computed tie rounds %s over rounds 1..%d" % (list(tied), len(run)))
    ck("every_tie_is_within_the_slow_equal_rate_class", slow_only and bool(tied),
       "all tied maximizers have rate 19/937 and equal height")


def check_all_tie_breakings():
    """"Every tie-breaking": branch on ALL maximizers, rounds 1..144.

    A layer is the set of exact reachable post-cut states after round t.  The
    canonical class of each layer must be a single state (that is the paper's
    tie-independence), and the maximum pre-cut height over the whole tree must be
    1878, reached only in round 60.
    """
    layer = {(0,) * 10}
    peak, peak_rounds, widths, visited = 0, set(), [], 0
    for t in range(1, P_HORIZON + 1):
        nxt = set()
        for state in layer:
            pre = tuple(state[i] + W[i] for i in range(10))
            top = max(pre)
            if top > peak:
                peak, peak_rounds = top, {t}
            elif top == peak:
                peak_rounds.add(t)
            for j in range(10):
                if pre[j] == top:
                    child = list(pre)
                    child[j] = 0
                    nxt.add(tuple(child))
        layer = nxt
        visited += len(layer)
        # width measured two ways: the paper's C, and independently the fast
        # coordinates in place together with the slow multiset
        widths.append(max(len({canon(s) for s in layer}),
                          len({s[:4] for s in layer}),
                          len({tuple(sorted(s[4:])) for s in layer})))
    ck("tie_independence_exhaustive",
       max(widths) == 1 and len(layer) >= 1,
       "%d exact states explored over all tie-breakings; canonical layer width "
       "max %d" % (visited, max(widths)))
    ck("peak_1878_under_every_tie_breaking",
       peak == P_PEAK and sorted(peak_rounds) == [P_PEAK_ROUND],
       "peak %d over the whole tie-breaking tree, only in round(s) %s"
       % (peak, sorted(peak_rounds)))


def check_period(run):
    """Seven-periodicity from round 137, its word, and its peak."""
    start = canon(run[P_CYCLE_START - 1]["post"])
    ret = canon(run[P_CYCLE_START + P_PERIOD - 1]["post"])
    cyc = run[P_CYCLE_START:P_CYCLE_START + P_PERIOD]
    word = "".join(symbol(r["cut"]) for r in cyc)
    long_ok = all(canon(run[t - 1]["post"]) == canon(run[t + P_PERIOD - 1]["post"])
                  for t in range(P_CYCLE_START, len(run) - P_PERIOD + 1))
    ck("canonical_trajectory_is_seven_periodic",
       start == ret and word == P_CYCLE_WORD and long_ok,
       "C(H(137)) = C(H(144)), word %s, and C(H(t+7)) = C(H(t)) for "
       "137 <= t <= %d" % (word, len(run) - P_PERIOD))
    cyc_peak = max(r["height"] for r in cyc)
    ck("cycle_peak_is_1252",
       cyc_peak == P_CYCLE_PEAK and cyc_peak < P_PEAK,
       "computed cycle peak %d, claimed %d; %d < %d so no round of the cycle "
       "beats round 60" % (cyc_peak, P_CYCLE_PEAK, cyc_peak, P_PEAK))
    tail = max(r["height"] for r in run[P_PEAK_ROUND:])
    overall = max(r["height"] for r in run)
    ck("supremum_over_long_horizon_is_1878",
       overall == P_PEAK and tail == P_CYCLE_PEAK,
       "over %d rounds the computed max is %d (claimed %d); after round 60 it "
       "never exceeds %d" % (len(run), overall, P_PEAK, tail))


def check_corollary(run):
    """n = 10+m for m = 1..12, simulated in exact rationals."""
    ten_word = "".join(symbol(r["cut"]) for r in run[:P_PEAK_ROUND])
    bad = []
    new_max = Fraction(0)
    for m in range(1, 13):
        rates = [Fraction(999, 1000) * Fraction(x, P_SUM) for x in W] \
            + [Fraction(1, 1000 * m)] * m
        big = simulate(rates, P_PEAK_ROUND, zero=Fraction(0))
        h60 = big[P_PEAK_ROUND - 1]["height"]
        word = "".join(symbol(r["cut"]) for r in big)
        new_max = max([new_max] + [r["pre"][i] for r in big
                                   for i in range(10, 10 + m)])
        if not (sum(rates) == 1 and all(x > 0 for x in rates)
                and all(r["cut"] < 10 for r in big) and word == ten_word
                and h60 == P_COR_ROUND60 and h60 > 2
                and h60 == Fraction(999, 1000) * Fraction(P_PEAK, P_SUM)
                and max(r["height"] for r in big) == h60):
            bad.append("m=%d h60=%s" % (m, h60))
    ck("corollary_positive_rate_instances_for_n_10_to_22", not bad,
       "; ".join(bad) or "12/12: sum 1, only original bamboos cut through "
       "round 60, round-60 height %s > 2" % P_COR_ROUND60)
    growth = Fraction(P_PEAK_ROUND, 1000)
    floor_sel = min(r["height"] for r in run[:P_PEAK_ROUND])
    barrier = Fraction(999, 1000) * Fraction(floor_sel, P_SUM)
    # new_max is the largest height actually reached by any added bamboo in the
    # simulated instances, so this compares two measured quantities; growth is
    # the closed-form ceiling t/(1000m) <= 60/1000 used by the proof for all m.
    ck("corollary_growth_bound_is_strict",
       growth < barrier and new_max <= growth and new_max < barrier
       and floor_sel == 313,
       "largest height reached by any added bamboo through round 60: %s; "
       "closed-form ceiling: %s; smallest selected original height: %s = "
       "%d/937 (the ceiling must be the smaller of the last two)"
       % (new_max, growth, barrier, floor_sel))


def greedy_peak(w, cap=4000):
    """Peak pre-cut height of lowest-index Reduce-Max, with cycle detection."""
    n = len(w)
    h = [0] * n
    seen = set()
    best = 0
    for _ in range(cap):
        pre = [h[i] + w[i] for i in range(n)]
        top = max(pre)
        if top > best:
            best = top
        h = pre
        h[h.index(top)] = 0
        key = tuple(h)
        if key in seen:
            return best, True
        seen.add(key)
    return best, False


def partitions(n, cap=None):
    """All integer partitions of n into positive parts, nonincreasing."""
    if cap is None:
        cap = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, cap), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def check_census():
    """Re-run of the enumeration that supported the conjecture: every integer
    partition of H in {5,10,...,35}.  None reaches backlog 2, and the exhibited
    instance has H = 937, outside that range."""
    worst = (Fraction(0), None)
    total = closed = 0
    for h in P_DDN_H:
        for w in partitions(h):
            total += 1
            peak, cycled = greedy_peak(w)
            closed += 1 if cycled else 0
            ratio = Fraction(peak, h)
            if ratio > worst[0]:
                worst = (ratio, w)
    ck("census_of_H_up_to_35_backlog_below_2_lowest_index_tie_breaking",
       total == 23297 and worst[0] < 2 and worst[0] > 0,
       "%d partitions enumerated; %d closed by a proven repeated state, %d only "
       "to a %d-round horizon; worst backlog %s at w=%s"
       % (total, closed, total - closed, 4000, worst[0], worst[1]))
    ck("exhibited_instance_outside_that_enumeration",
       P_SUM not in P_DDN_H and P_SUM > max(P_DDN_H)
       and sum(W) == P_SUM and len(W) == 10,
       "integer total growth H = %d exceeds the enumerated maximum %d, while "
       "n = 10" % (P_SUM, max(P_DDN_H)))


def main():
    run = simulate(W, 3000)
    check_instance()
    check_cut_word(run)
    check_table_columns(run)
    check_endpoints(run)
    check_violation(run)
    check_decisive(run)
    check_ties(run)
    check_all_tie_breakings()
    check_period(run)
    check_corollary(run)
    check_census()
    print("NOT RE-RUN: (a) the 2000- and 2702-bamboo instances of the earlier "
          "disproof, and their reported bounds 2.0004 and 2.076, are not "
          "reconstructed, since those rate vectors are not exhibited in the "
          "text being checked; (b) the corollary is simulated exactly for "
          "n = 10..22 only, its general growth bound t/(1000m) <= 60/1000 "
          "being verified in closed form for all m >= 1; (c) nothing about "
          "fewer than ten bamboos is claimed or tested; (d) in the census of "
          "H <= 35 the backlog is exact only for the partitions whose state "
          "recurred within the 4000-round horizon, is a horizon-limited lower "
          "bound for the remaining 95, and is computed under lowest-index "
          "tie-breaking alone -- that census only corroborates the remark "
          "about the earlier search and carries no part of the main theorem, "
          "whose tie-independence is established by exhaustive expansion above.")
    return run


if __name__ == "__main__":
    main()
    n = len(CHECKS)
    k = sum(1 for _, ok in CHECKS if not ok)
    if k == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
    else:
        print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
    raise SystemExit(0 if k == 0 else 1)

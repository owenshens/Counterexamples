#!/usr/bin/env python3
"""Verification of small-parameter counterexamples to a conjecture on
subtraction Nim with a one-time pass (three parameterized families).

Standard library only; all arithmetic is exact integer arithmetic on
Grundy values.  No floats are used in any decision.

One caveat on the inputs below.  The wording of the conjecture itself --
the three identities and the three families -- is a transcription that
nothing here can check against its source, since no external text is
read.  What is checked is that the transcription is a contingent
statement rather than a tautology or an absurdity: the final control
confirms that the very same identities DO hold, over thousands of
evaluations, at parameters the paper does not single out.  Every other
transcribed value (table entries, thresholds, mex evaluations, lemma
closed forms, census, scope) is compared against a value computed here.

TAKEN FROM THE PAPER (inputs, transcribed and then tested):
  * the two-argument Grundy recursion for a subtraction game with a
    one-time pass:
        G(x,0) = mex{ G(x-s,0) : s in S, s <= x },
        G(x,1) = mex( { G(x-s,1) : s in S, s <= x } u { G(x,0) } ),
    with G(0,1) = 0;
  * the three conjectured identities, printed there as (35), (36), (37):
        (35) G(x,0) = mex{ G(x+s,0) : s in S },
        (36) G(x,1) = mex{ G(x-s,1) : s in S },
        (37) G(x,1) = mex{ G(x+s,1) : s in S };
  * the three families and thresholds of the conjecture's subparts:
        (a)   S = {a, 2an, 2an+a},      threshold 6an+(4a+1),
        (b)   S = {a, (2n+1)a, (2n+3)a}, threshold a(2n+3)+1,
        (iii) S = {a, (2n+1)a, (2n+5)a}, threshold a(2n+5)+1,
    each quantified over a, n in N, each a conjunction that also asserts
    (35) for all x and, in (b) and (iii), G(x,1) < 4 for all x;
  * the exhibited witnesses: (a) n=2, x=18a; (b) n=1, x=6a;
    (iii) n=2, x=10a;
  * the printed Grundy table for the subtraction set {1,4,5}, x = 0..23
    (both rows), and the printed mex evaluations mex{2,3,1}=0,
    mex{2,1,0}=3, mex{2,2,2}=0 with G(18,1)=4, G(6,1)=1, G(10,1)=1;
  * the two closed forms of the paper's lemmas (scaling, and the odd
    subtraction sets containing 1);
  * the scope claims: which members of the three families have the form
    {2,4m,4m+2}, and the census "for a=1 and n=1..12 there are no
    exceptional n other than 2, 1 and 2 respectively".

DERIVED HERE (recomputed from the recursion; nothing below is copied):
  * the Grundy engine is cross-validated against an independent
    win/loss (boolean) recursion on both layers of the pass game;
  * the printed table is regenerated from the recursion and compared
    entry by entry;
  * both lemmas are re-derived over a range;
  * for each subpart the hypotheses (family membership, distinctness,
    threshold, witness above threshold) are checked, and then the
    conclusion is COMPUTED at the witness and found violated, for
    a = 1..20 and under both possible conventions for a pass at a
    positive position with no legal subtraction;
  * the conjuncts the paper does NOT contradict are confirmed to hold
    at the witnesses ((35), and the bound G(x,1) < 4);
  * the paper's threshold specializations are recomputed, and so is its
    census, over n <= CENSUS_NMAX and threshold <= x <= CENSUS_XMAX --
    the paper's census sentence is unbounded in x, so for the
    non-exceptional n this program certifies the two identities only up
    to that cap; the closing note says so.

HARNESS INTEGRITY: several stages register one check per subpart inside a
loop, so a crash on the second or third subpart could leave the first
check standing and the rest missing.  The driver therefore (i) records a
FAILED check, with the traceback, for any stage that raises, (ii) records
a FAILED check for any stage whose number of registered checks differs
from the number declared in the stage table, and (iii) reconciles the
total against the count derived from that table before printing the
verdict, refusing to exit 0 on any mismatch.  The same reconciliation
check pins SUBPART_KEYS -- which drives both the loops and the declared
counts -- to the transcribed subpart dictionaries, so coverage cannot
shrink while declaration and execution still agree.  Every sweep bound
that appears in a printed line is a module constant used by the loop that
does the work, so no reported range can drift from the executed one.
"""

import sys
import traceback

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)


# ---------------------------------------------------------------- inputs
# Values transcribed from the paper.  Every one of them is compared below
# against a value this program computes from the recursion.

PAPER_TABLE_145 = {
    0: [0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3],
    1: [0, 2, 1, 0, 1, 4, 0, 2, 3, 0, 1, 3, 0, 1, 3, 0, 1, 2, 4, 2, 3, 0, 1, 0],
}

# subpart label -> (base subtraction set at scale 1, value of n, witness
# multiplier x/a, family generator, threshold generator, does the paper
# claim (37) fails here too?)
SUBPARTS = {
    "a": dict(base=(1, 4, 5), n=2, xmul=18,
              family=lambda a, n: (a, 2 * a * n, 2 * a * n + a),
              thresh=lambda a, n: 6 * a * n + (4 * a + 1),
              forward_fails=True),
    "b": dict(base=(1, 3, 5), n=1, xmul=6,
              family=lambda a, n: (a, (2 * n + 1) * a, (2 * n + 3) * a),
              thresh=lambda a, n: a * (2 * n + 3) + 1,
              forward_fails=False),
    "iii": dict(base=(1, 5, 9), n=2, xmul=10,
                family=lambda a, n: (a, (2 * n + 1) * a, (2 * n + 5) * a),
                thresh=lambda a, n: a * (2 * n + 5) + 1,
                forward_fails=False),
}

# The three subparts, in the order the paper prints them.  The stages that
# loop over these register one check per key; the driver derives its
# expected check count from this tuple.
SUBPART_KEYS = ("a", "b", "iii")

PAPER_MEX_CLAIMS = {
    # subpart -> (multiset feeding the backward mex at the witness,
    #             its mex, the Grundy value at the witness)
    "a": ([2, 3, 1], 0, 4),
    "b": ([2, 2, 2], 0, 1),
    "iii": ([2, 2, 2], 0, 1),
}
PAPER_FORWARD_CLAIM_A = ([2, 1, 0], 3, 4)   # part (a), forward mex at 18
PAPER_CENSUS = {"a": [2], "b": [1], "iii": [2]}   # a=1, n=1..12
CENSUS_NMAX = 12
CENSUS_XMAX = 400
AMAX = 20

# Sweep bounds.  Each of these is used BOTH in the loop that performs the
# work and in the text that reports it, so no reported bound can drift
# away from the bound actually executed.
ENGINE_XMAX = 120        # cross-check of the Grundy engine, x = 0..this
SCALE_AMAX = 8           # scaling lemma, a = 1..this
SCALE_KMAX = 60          # scaling lemma, k = 0..this
SCOPE_MAX = 60           # family-shape sweep, a and n = 1..this
M_SWEEP = 400            # m range used to build the shape {2,4m,4m+2}
M_MIN_PROVED = 3         # the cited theorems treat that shape only from here
BOUND_LT = 4             # the bound conjunct of parts (b) and (iii)


def mex(values):
    """Least non-negative integer not in the finite set `values`."""
    seen = set(values)
    m = 0
    while m in seen:
        m += 1
    return m


def grundy(S, N, stuck_pass=True):
    """Grundy values G(x,0) and G(x,1) for x = 0..N.

    `stuck_pass` selects the convention at a positive position with no
    legal subtraction: True lets the mover spend the pass (so the only
    option is G(x,0)), False forbids it (no options at all).  Every
    position used by the counterexamples admits the subtraction min(S),
    so the checks below are run under both settings.
    """
    S = sorted(set(S))
    g0 = [0] * (N + 1)
    g1 = [0] * (N + 1)
    for x in range(1, N + 1):
        opts = [x - s for s in S if s <= x]
        if not opts:
            g0[x] = 0
            g1[x] = mex([g0[x]]) if stuck_pass else 0
            continue
        g0[x] = mex(g0[y] for y in opts)
        g1[x] = mex([g0[x]] + [g1[y] for y in opts])
    return g0, g1


def winning(S, N, stuck_pass=True):
    """Independent boolean recursion: is (x, p) a win for the mover?

    p = 1 means the one-time pass is still available.  A position is a
    win iff some option is a loss.  This never mentions mex, so it is a
    genuine cross-check of the Grundy engine through the Sprague-Grundy
    criterion "value 0 iff loss".
    """
    S = sorted(set(S))
    w0 = [False] * (N + 1)
    w1 = [False] * (N + 1)
    for x in range(1, N + 1):
        opts = [x - s for s in S if s <= x]
        w0[x] = any(not w0[y] for y in opts)
        if not opts:
            w1[x] = (not w0[x]) if stuck_pass else False
        else:
            w1[x] = (not w0[x]) or any(not w1[y] for y in opts)
    return w0, w1


def check_engine():
    """Grundy value 0 must coincide with "the mover loses", on both layers."""
    bad = []
    tested = 0
    N = ENGINE_XMAX
    sets = [(1,), (1, 2), (1, 4, 5), (1, 3, 5), (1, 5, 9), (2, 8, 10),
            (3, 12, 15), (2, 3, 7), (4, 5)]
    for S in sets:
        for stuck in (True, False):
            g0, g1 = grundy(S, N, stuck)
            w0, w1 = winning(S, N, stuck)
            for x in range(N + 1):
                tested += 2
                if (g0[x] == 0) != (not w0[x]):
                    bad.append((S, stuck, x, 0))
                if (g1[x] == 0) != (not w1[x]):
                    bad.append((S, stuck, x, 1))
    ck("engine_matches_independent_win_loss_recursion", not bad,
       "%d position-layers, %d mismatches" % (tested, len(bad)))


def check_table_wellformed():
    """The exhibited object: decode the printed table and print it back."""
    ncols = len(PAPER_TABLE_145[0])
    ok = sorted(PAPER_TABLE_145) == [0, 1] and ncols == 24
    for row in (0, 1):
        vals = PAPER_TABLE_145[row]
        ok = ok and len(vals) == ncols and all(
            isinstance(v, int) and v >= 0 for v in vals)
    S = SUBPARTS["a"]["base"]
    ok = ok and S == (1, 4, 5) and len(set(S)) == 3
    print("  exhibited table, subtraction set {%s}:"
          % ",".join(map(str, sorted(S))))
    print("    x        " + " ".join("%2d" % x for x in range(ncols)))
    print("    G(x,0)   " + " ".join("%2d" % v for v in PAPER_TABLE_145[0]))
    print("    G(x,1)   " + " ".join("%2d" % v for v in PAPER_TABLE_145[1]))
    ck("printed_table_wellformed_24_columns_two_rows", ok,
       "rows=%d entries=%d+%d" % (len(PAPER_TABLE_145),
                                  len(PAPER_TABLE_145[0]),
                                  len(PAPER_TABLE_145[1])))


def check_table_recomputed():
    """Regenerate the table from the recursion and compare entry by entry."""
    S = tuple(sorted(SUBPARTS["a"]["base"]))
    xmax = max(len(PAPER_TABLE_145[r]) for r in (0, 1)) - 1
    g0, g1 = grundy(S, xmax)
    rows = {0: g0, 1: g1}
    entries = sum(len(PAPER_TABLE_145[r]) for r in (0, 1))
    diffs = [(x, r) for r in (0, 1) for x in range(len(PAPER_TABLE_145[r]))
             if PAPER_TABLE_145[r][x] != rows[r][x]]
    ck("printed_table_reproduced_by_recursion", not diffs and entries > 0,
       "%d entries, %d discrepancies" % (entries, len(diffs)))


def check_scaling_lemma():
    """G_{aS}(ak,p) = G_S(k,p) whenever 1 in S, for both conventions."""
    bad = []
    tested = 0
    kmax = SCALE_KMAX
    for S in [(1,), (1, 2, 3), (1, 4, 5), (1, 3, 5), (1, 5, 9), (1, 2, 7)]:
        base0, base1 = grundy(S, kmax)
        for a in range(1, SCALE_AMAX + 1):
            for stuck in (True, False):
                aS = tuple(a * s for s in S)
                g0, g1 = grundy(aS, a * kmax, stuck)
                for k in range(kmax + 1):
                    tested += 2
                    if g0[a * k] != base0[k]:
                        bad.append((S, a, k, 0))
                    if g1[a * k] != base1[k]:
                        bad.append((S, a, k, 1))
    ck("scaling_lemma_holds_on_multiples_of_a", not bad,
       "%d comparisons, %d failures" % (tested, len(bad)))


def odd_closed_form(S, x):
    """The paper's closed form for odd subtraction sets containing 1."""
    g0 = 0 if x % 2 == 0 else 1
    if x == 0:
        g1 = 0
    elif x % 2 == 0:
        g1 = 1
    elif x in S:
        g1 = 2
    else:
        g1 = 0
    return g0, g1


def check_odd_lemma():
    """Re-derive the closed form for every 3-element odd set containing 1."""
    xmax = 200
    sets = []
    for b in range(3, 20, 2):
        for c in range(b + 2, 22, 2):
            sets.append((1, b, c))
    bad = []
    tested = 0
    for S in sets:
        g0, g1 = grundy(S, xmax)
        for x in range(xmax + 1):
            tested += 2
            f0, f1 = odd_closed_form(set(S), x)
            if g0[x] != f0:
                bad.append((S, x, 0))
            if g1[x] != f1:
                bad.append((S, x, 1))
    ck("odd_set_lemma_closed_form_correct", not bad and len(sets) >= 40,
       "%d odd sets, %d comparisons, %d failures" % (len(sets), tested,
                                                     len(bad)))


def subpart_data(key, a, stuck=True):
    """Family set, threshold, witness and the three quantities of (36),(37)."""
    d = SUBPARTS[key]
    n = d["n"]
    S = tuple(sorted(d["family"](a, n)))
    thr = d["thresh"](a, n)
    x = d["xmul"] * a
    g0, g1 = grundy(S, x + max(S), stuck)
    lhs1 = g1[x]
    back = mex(g1[x - s] for s in S)
    fwd = mex(g1[x + s] for s in S)
    return dict(S=S, n=n, thr=thr, x=x, g0=g0, g1=g1,
                lhs1=lhs1, back=back, fwd=fwd)


def check_hypotheses():
    """Each witness must satisfy every hypothesis of the subpart it attacks."""
    for key in SUBPART_KEYS:
        d = SUBPARTS[key]
        bad = []
        for a in range(1, AMAX + 1):
            r = subpart_data(key, a)
            S, x, thr = r["S"], r["x"], r["thr"]
            scaled = tuple(sorted(a * s for s in d["base"]))
            conds = [
                a >= 1,                       # a in N
                r["n"] >= 1,                  # n in N, inside the quantifier
                len(set(S)) == 3,             # a genuine 3-element set
                S == scaled,                  # family member = a * base set
                all(s >= 1 for s in S),       # subtraction set in N
                x >= thr,                     # witness above the threshold
                x >= max(S),                  # (36) uses only x-s >= 0
                min(S) <= x,                  # witness admits a subtraction
            ]
            if not all(conds):
                bad.append((a, [i for i, c in enumerate(conds) if not c]))
        r1 = subpart_data(key, 1)
        ck("hypotheses_satisfied_subpart_%s" % key, not bad,
           "a=1..%d; at a=1: S=%s n=%d threshold=%d witness x=%d; %d failures"
           % (AMAX, "{%s}" % ",".join(map(str, r1["S"])), r1["n"], r1["thr"],
              r1["x"], len(bad)))


def check_refutation_backward():
    """LOAD-BEARING: (36) is computed at each witness and found violated."""
    for key in SUBPART_KEYS:
        bad = []
        sample = None
        for a in range(1, AMAX + 1):
            for stuck in (True, False):
                r = subpart_data(key, a, stuck)
                if r["back"] == r["lhs1"]:
                    bad.append((a, stuck))
                if a == 1 and stuck:
                    sample = r
        ck("equation36_fails_at_witness_subpart_%s" % key, not bad,
           "a=1..%d x2 conventions; at a=1: mex{G(x-s,1)}=%d vs G(x,1)=%d "
           "at x=%d; %d parameters where (36) did NOT fail"
           % (AMAX, sample["back"], sample["lhs1"], sample["x"], len(bad)))


def check_refutation_forward():
    """(37) at the witnesses: must fail for (a), and hold for (b) and (iii)."""
    for key in SUBPART_KEYS:
        want_fail = SUBPARTS[key]["forward_fails"]
        bad = []
        sample = None
        for a in range(1, AMAX + 1):
            for stuck in (True, False):
                r = subpart_data(key, a, stuck)
                differs = (r["fwd"] != r["lhs1"])
                if differs != want_fail:
                    bad.append((a, stuck))
                if a == 1 and stuck:
                    sample = r
        name = ("equation37_fails_at_witness_subpart_%s" % key if want_fail
                else "equation37_still_holds_at_witness_subpart_%s" % key)
        ck(name, not bad,
           "a=1..%d x2 conventions; at a=1: mex{G(x+s,1)}=%d vs G(x,1)=%d "
           "(paper claims %s); %d deviations"
           % (AMAX, sample["fwd"], sample["lhs1"],
              "failure" if want_fail else "no failure", len(bad)))


def check_printed_mex_evaluations():
    """Every displayed Grundy value and mex evaluation, recomputed at a=1."""
    bad = []
    detail = []
    for key in SUBPART_KEYS:
        r = subpart_data(key, 1)
        got = [r["g1"][r["x"] - s] for s in r["S"]]
        want_list, want_mex, want_lhs = PAPER_MEX_CLAIMS[key]
        if got != want_list:
            bad.append((key, "options", got, want_list))
        if mex(got) != want_mex:
            bad.append((key, "mex", mex(got), want_mex))
        if r["lhs1"] != want_lhs:
            bad.append((key, "lhs", r["lhs1"], want_lhs))
        detail.append("%s:mex%s=%d vs %d" % (key, tuple(got), mex(got),
                                            r["lhs1"]))
    ra = subpart_data("a", 1)
    gotf = [ra["g1"][ra["x"] + s] for s in ra["S"]]
    wl, wm, wlhs = PAPER_FORWARD_CLAIM_A
    if gotf != wl or mex(gotf) != wm or ra["lhs1"] != wlhs:
        bad.append(("a", "forward", gotf, wl))
    ck("printed_mex_evaluations_recomputed", not bad,
       "; ".join(detail) + "; forward a:mex%s=%d; %d mismatches"
       % (tuple(gotf), mex(gotf), len(bad)))


def check_surviving_conjunct_35():
    """(35) is not contradicted: it holds at the witnesses and on a range."""
    bad_wit = []
    for key in SUBPART_KEYS:
        for a in range(1, AMAX + 1):
            for stuck in (True, False):
                r = subpart_data(key, a, stuck)
                if r["g0"][r["x"]] != mex(r["g0"][r["x"] + s] for s in r["S"]):
                    bad_wit.append((key, a, stuck))
    bad_rng = []
    xmax = 300
    for key in SUBPART_KEYS:
        S = tuple(sorted(SUBPARTS[key]["base"]))
        g0, _ = grundy(S, xmax + max(S))
        for x in range(xmax + 1):
            if g0[x] != mex(g0[x + s] for s in S):
                bad_rng.append((key, x))
    ck("equation35_conjunct_holds_at_witnesses_and_on_range",
       not bad_wit and not bad_rng,
       "witnesses a=1..%d x2 conventions: %d failures; a=1 and x<=%d: "
       "%d failures" % (AMAX, len(bad_wit), xmax, len(bad_rng)))


def check_surviving_conjunct_bound():
    """(b) and (iii) also assert G(x,1) < 4; at their witnesses it holds."""
    xmax = 400
    bad = []
    worst = 0
    shown = []
    for key in ("b", "iii"):
        S = tuple(sorted(SUBPARTS[key]["base"]))
        shown.append("{%s}" % ",".join(map(str, S)))
        _, g1 = grundy(S, xmax)
        worst = max(worst, max(g1))
        bad += [(key, x) for x in range(xmax + 1) if g1[x] >= BOUND_LT]
    ck("bound_conjunct_G_lt_%d_holds_for_odd_witness_sets" % BOUND_LT, not bad,
       "a=1, x<=%d, sets %s: max G(x,1)=%d, %d violations"
       % (xmax, " and ".join(shown), worst, len(bad)))


def check_mixed_reading():
    """The mixed reading of (37) also fails at the part (a) witness."""
    r = subpart_data("a", 1)
    rhs = mex(r["g1"][r["x"] + s] for s in r["S"])
    lhs0 = r["g0"][r["x"]]
    # The two pinned values are not magic numbers: the forward mex is the one
    # the paper prints for part (a), and G(x,0) is that column of the paper's
    # own table (itself reproduced from the recursion above).
    want_rhs = PAPER_FORWARD_CLAIM_A[1]
    want_lhs0 = PAPER_TABLE_145[0][r["x"]]
    ck("mixed_reading_of_equation37_also_fails",
       rhs != lhs0 and rhs == want_rhs and lhs0 == want_lhs0,
       "x=%d: mex{G(x+s,1)}=%d, G(x,0)=%d" % (r["x"], rhs, lhs0))


def check_convention_independence():
    """The witnesses are insensitive to the pass convention at stuck positions."""
    bad = []
    for key in SUBPART_KEYS:
        for a in range(1, AMAX + 1):
            r0 = subpart_data(key, a, True)
            r1 = subpart_data(key, a, False)
            if (r0["lhs1"], r0["back"], r0["fwd"]) != \
               (r1["lhs1"], r1["back"], r1["fwd"]):
                bad.append((key, a))
            S, x = r0["S"], r0["x"]
            if not all(min(S) <= y for y in [x] + [x - s for s in S if s < x]
                       + [x + s for s in S]):
                bad.append((key, a, "stuck position used"))
    ck("witness_values_independent_of_pass_convention", not bad,
       "%d subparts x a=1..%d: %d disagreements"
       % (len(SUBPART_KEYS), AMAX, len(bad)))


def check_threshold_specializations():
    """The paper's specialized thresholds and the witness margins."""
    bad = []
    for a in range(1, AMAX + 1):
        if SUBPARTS["a"]["thresh"](a, 2) != 16 * a + 1:
            bad.append(("a", a))
        if SUBPARTS["b"]["thresh"](a, 1) != 5 * a + 1:
            bad.append(("b", a))
        if SUBPARTS["iii"]["thresh"](a, 2) != 9 * a + 1:
            bad.append(("iii", a))
        for key in ("b", "iii"):
            r = subpart_data(key, a)
            if r["x"] != min(r["S"]) + max(r["S"]):
                bad.append((key, a, "not min+max"))
            if r["x"] - r["thr"] != a - 1:
                bad.append((key, a, "margin"))
        if subpart_data("a", a)["x"] - SUBPARTS["a"]["thresh"](a, 2) != 2 * a - 1:
            bad.append(("a", a, "margin"))
    eq_at_1 = all(subpart_data(k, 1)["x"] == subpart_data(k, 1)["thr"]
                  for k in ("b", "iii"))
    ck("threshold_specializations_and_margins", not bad and eq_at_1,
       "a=1..%d: 6an+(4a+1)|n=2 = 16a+1, a(2n+3)+1|n=1 = 5a+1, "
       "a(2n+5)+1|n=2 = 9a+1; margins 2a-1, a-1, a-1; %d failures"
       % (AMAX, len(bad)))


def exceptional_n(key, nmax, xmax):
    """Which n in 1..nmax break (36) or (37) at a=1, for some x in
    [threshold, xmax].  A break at x > xmax is not detected, so "no
    exceptional n other than these" is a statement about that window."""
    out = []
    d = SUBPARTS[key]
    for n in range(1, nmax + 1):
        S = tuple(sorted(set(d["family"](1, n))))
        if len(S) < 3:
            continue
        thr = d["thresh"](1, n)
        _, g1 = grundy(S, xmax + max(S))
        broke = False
        for x in range(thr, xmax + 1):
            if max(S) <= x and g1[x] != mex(g1[x - s] for s in S):
                broke = True
                break
            if g1[x] != mex(g1[x + s] for s in S):
                broke = True
                break
        if broke:
            out.append(n)
    return out


def check_census():
    """The paper's census: at a=1 the only exceptional n are 2, 1 and 2."""
    got = {}
    bad = []
    for key in SUBPART_KEYS:
        got[key] = exceptional_n(key, CENSUS_NMAX, CENSUS_XMAX)
        if got[key] != PAPER_CENSUS[key]:
            bad.append((key, got[key], PAPER_CENSUS[key]))
    ck("census_of_exceptional_n_matches_paper", not bad,
       "a=1, n=1..%d, x in [threshold,%d]: exceptional n = %s; the paper's "
       "census is unbounded in x, so for the other n this certifies (36) and "
       "(37) only up to x=%d"
       % (CENSUS_NMAX, CENSUS_XMAX,
          " ".join("%s:%s" % (k, got[k]) for k in SUBPART_KEYS), CENSUS_XMAX))


def check_scope_family_overlap():
    """Which members of the three families have the form {2,4m,4m+2}.

    The cited theorems treat that family only from m = 3 onwards, so
    "no witness is covered by what is proved" is NOT settled by showing
    that a witness merely has the shape {2,4m,4m+2}: it must miss the
    shape for every m >= 3.  Both sets are therefore built, and the
    load-bearing conjunct is that no witness lands in the second.
    """
    form = set((2, 4 * m, 4 * m + 2) for m in range(1, M_SWEEP))
    treated = set((2, 4 * m, 4 * m + 2) for m in range(M_MIN_PROVED, M_SWEEP))
    hits = dict((key, set()) for key in SUBPART_KEYS)
    covered = []
    for key in SUBPART_KEYS:
        n_wit = SUBPARTS[key]["n"]
        for a in range(1, SCOPE_MAX + 1):
            for n in range(1, SCOPE_MAX + 1):
                S = tuple(sorted(set(SUBPARTS[key]["family"](a, n))))
                if len(S) == 3 and S in form:
                    hits[key].add((a, n))
            W = tuple(sorted(set(SUBPARTS[key]["family"](a, n_wit))))
            if W in treated:
                covered.append((key, a, W))
    a_values = sorted({a for a, _ in hits["a"]})
    witness_a2 = tuple(sorted(SUBPARTS["a"]["family"](2, 2)))
    ok = (a_values == [2] and not hits["b"] and not hits["iii"]
          and witness_a2 == (2, 8, 10) and witness_a2 in form
          and witness_a2 not in treated and not covered)
    ck("scope_overlap_with_proved_family", ok,
       "a=1..%d,n=1..%d: part (a) meets the shape only at a in %s (a=2,n=2 "
       "gives {%s}, i.e. m=%d, below the m>=%d range proved); parts "
       "(b),(iii) shape matches: %d,%d; witnesses inside the m>=%d range: %d"
       % (SCOPE_MAX, SCOPE_MAX, a_values,
          ",".join(map(str, witness_a2)), witness_a2[1] // 4, M_MIN_PROVED,
          len(hits["b"]), len(hits["iii"]), M_MIN_PROVED, len(covered)))


def check_identity_control():
    """Control: at non-exceptional parameters (36) and (37) do hold.

    Without this, a mis-transcribed identity would make the refutation
    checks pass vacuously for every parameter.
    """
    bad = []
    tested = 0
    controls = [("a", 3), ("a", 5), ("b", 2), ("b", 4), ("iii", 1), ("iii", 3)]
    xmax = 400
    for key, n in controls:
        d = SUBPARTS[key]
        S = tuple(sorted(set(d["family"](1, n))))
        thr = d["thresh"](1, n)
        _, g1 = grundy(S, xmax + max(S))
        for x in range(thr, xmax + 1):
            tested += 2
            if max(S) <= x and g1[x] != mex(g1[x - s] for s in S):
                bad.append((key, n, x, 36))
            if g1[x] != mex(g1[x + s] for s in S):
                bad.append((key, n, x, 37))
    ck("control_identities_hold_at_non_exceptional_parameters",
       not bad and tested > 1000,
       "%d identity evaluations at %d control parameters, %d failures"
       % (tested, len(controls), len(bad)))


def main():
    print("Subtraction Nim with a one-time pass: verification of the")
    print("counterexamples to the three subparts of the conjecture.")
    print("")
    # (label, stage, number of checks the stage is declared to register).
    # The three subpart stages register one check per subpart, so their
    # declared count is derived from SUBPART_KEYS rather than written out.
    nsub = len(SUBPART_KEYS)
    stages = [
        ("engine_matches_independent_win_loss_recursion", check_engine, 1),
        ("printed_table_wellformed_24_columns_two_rows",
         check_table_wellformed, 1),
        ("printed_table_reproduced_by_recursion", check_table_recomputed, 1),
        ("scaling_lemma_holds_on_multiples_of_a", check_scaling_lemma, 1),
        ("odd_set_lemma_closed_form_correct", check_odd_lemma, 1),
        ("hypotheses_satisfied", check_hypotheses, nsub),
        ("equation36_fails_at_witness", check_refutation_backward, nsub),
        ("equation37_at_witness", check_refutation_forward, nsub),
        ("printed_mex_evaluations_recomputed",
         check_printed_mex_evaluations, 1),
        ("equation35_conjunct_holds", check_surviving_conjunct_35, 1),
        ("bound_conjunct_G_lt_%d" % BOUND_LT,
         check_surviving_conjunct_bound, 1),
        ("mixed_reading_of_equation37_also_fails", check_mixed_reading, 1),
        ("witness_values_independent_of_pass_convention",
         check_convention_independence, 1),
        ("threshold_specializations_and_margins",
         check_threshold_specializations, 1),
        ("census_of_exceptional_n_matches_paper", check_census, 1),
        ("scope_overlap_with_proved_family", check_scope_family_overlap, 1),
        ("control_identities_hold_at_non_exceptional_parameters",
         check_identity_control, 1),
    ]
    # The "+ 1" is the harness reconciliation check registered after the
    # loop, which counts itself.
    expected = sum(want for _, _, want in stages) + 1
    trouble = []
    for label, fn, want in stages:
        before = len(CHECKS)
        try:
            fn()
        except Exception as exc:      # a crash is ALWAYS a failed check, even
            trouble.append(label)     # if the stage already logged some checks
            print("  TRACEBACK for stage %s:" % label)
            for line in traceback.format_exc().rstrip().splitlines():
                print("    " + line)
            ck("stage_%s_raised_an_exception" % label, False,
               "aborted after %d of its %d declared checks: %s: %s"
               % (len(CHECKS) - before, want, type(exc).__name__, exc))
            continue
        got = len(CHECKS) - before
        if got != want:
            trouble.append(label)
            ck("stage_%s_registered_unexpected_check_count" % label, False,
               "declared %d checks, registered %d" % (want, got))
    # The per-subpart stages are driven by SUBPART_KEYS, and their declared
    # count is len(SUBPART_KEYS); so if SUBPART_KEYS ever failed to list every
    # transcribed subpart, coverage would silently shrink with declaration and
    # execution still agreeing.  Pin it to the transcribed dictionaries.
    keys = tuple(sorted(SUBPART_KEYS))
    key_coverage = (len(SUBPART_KEYS) == len(set(SUBPART_KEYS))
                    and keys == tuple(sorted(SUBPARTS))
                    and keys == tuple(sorted(PAPER_MEX_CLAIMS))
                    and keys == tuple(sorted(PAPER_CENSUS)))
    ck("harness_registered_every_declared_check",
       len(CHECKS) + 1 == expected and not trouble and key_coverage,
       "%d stages declaring %d checks, plus this reconciliation = %d expected;"
       " %d registered; %d stage(s) aborted or miscounted; %d subpart keys "
       "%s every transcribed subpart"
       % (len(stages), expected - 1, expected, len(CHECKS) + 1, len(trouble),
          len(SUBPART_KEYS), "cover" if key_coverage else "DO NOT cover"))
    print("")
    print("NOT RE-RUN HERE: the conjecture quantifies over all a and n and "
          "all x above the threshold.")
    print("  Refutation is computed exactly for a = 1..%d (both pass "
          "conventions); the reduction" % AMAX)
    print("  to a = 1 for all larger a rests on the scaling lemma, which is "
          "verified numerically")
    print("  for a <= %d and k <= %d but not proved here.  The census over n "
          "is finite (n <= %d) and" % (SCALE_AMAX, SCALE_KMAX, CENSUS_NMAX))
    print("  is also bounded in x, where the paper's census sentence is not: "
          "that sentence says that")
    print("  for every other n the two identities hold above the stated "
          "threshold, i.e. at every")
    print("  x >= threshold, whereas each non-exceptional n is certified here "
          "only for threshold <= x")
    print("  <= %d, so a break at larger x would not be seen.  The three "
          "exceptional n the theorem" % CENSUS_XMAX)
    print("  exhibits are settled outright at their witnesses and do not rest "
          "on the census.  The")
    print("  eventual-periodicity conjuncts are not tested, and the paper "
          "asserts nothing about them.")
    print("  No external catalogue or table is read.")
    print("")
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if n != expected:
        print("HARNESS ERROR: %d checks executed, %d expected from the stage "
              "table; the verdict below is not a clean pass." % (n, expected))
    if bad or n != expected:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())

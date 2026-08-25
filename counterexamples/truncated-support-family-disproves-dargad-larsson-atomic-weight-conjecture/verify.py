#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- verification program for

    "An Infinite Family of Counterexamples to an Atomic-Weight Conjecture
     of Dargad and Larsson"

Claim under test (paper, Theorem 1 and Corollary):  for every integer t >= 1,
writing H_n = TS_t(n; 1, 2t+1) for the truncated-support partizan subtraction
game,

        H_{2t+2} = 0,   H_{2t+3} = *,   H_{2t+4} = down,

so aw(H_{2t+4}) = -1, whereas Conjecture 7.1 of Dargad-Larsson (arXiv
2607.27989v1) predicts -floor((n - R_betahat)/a) = 0 at n = R_betahat = 2t+4.

--------------------------------------------------------------------------
TAKEN FROM THE PAPER (data; not re-derived here)
--------------------------------------------------------------------------
 P1. The rules of TS_tau(n; a, b)  [paper section 1, citing DL Definition 5.1]:
       TS(0) = 0 (empty option sets);
       Left options  = { TS(n-k) : 1       <= k <= a, n-k >= 0 };
       Right options = { TS(n-k) : tau+1   <= k <= b, n-k >= 0 };
     normal play (a player with no legal move loses).
 P2. The literal definitions  * = {0|0},  up = {0|*},  down = {*|0}.
 P3. The conjectured formula and its ingredients  [paper eq. (1), DL (5),(13)]:
       aw(TS_tau(n;a,b)) = -floor((n - R_betahat)/a),
       betahat = floor((b-a)/(b-a-tau)),   R_betahat = betahat*(a+tau+1),
     under the standing hypotheses 0 < a <= b, a < b, tau < b-a, n >= R_betahat.
 P4. The paper's asserted values, used ONLY as the right-hand side of the
     comparisons: H_{2t+2}=0, H_{2t+3}=*, H_{2t+4}=down; betahat=2;
     R_betahat=2t+4; predicted atomic weight 0; aw(down) = -1.
 P5. Additional instances beyond the paper's statements, again only as
     right-hand sides: TS_t(2t+5;1,2t+1) = down+down+*  for t = 1..4
     (predicted aw -1); TS_tau(R_betahat;1,b) = down for (b,tau) in
     {(4,2),(5,3),(7,4),(9,6)}.
 P6. Two standard atomic-weight facts, cited by the paper as DL Definition
     A.17 and Lemma A.19 (= Siegel, Prop. 7.12):
       (A) aw(up) = 1, aw(-G) = -aw(G), aw(*) = 0;
       (B) aw is additive: aw(G+H) = aw(G) + aw(H).
     Together with the standard characterisation
       (C) for an all-small G,  aw(G) = 0  <=>  G is infinitesimal with
           respect to up, i.e. for EVERY integer j >= 1 the game j*G is
           less than or confused with up, and greater than or confused with
           down.  The strict form "down < j*G < up" is NOT the criterion and
           is false: aw(*) = 0 while * is confused with up and with down.
     these are the ONLY imported statements about atomic weight; every use of
     them is printed explicitly where it happens.

--------------------------------------------------------------------------
DERIVED HERE (this is what the checks actually decide)
--------------------------------------------------------------------------
 D1. The game trees of TS_tau(n;a,b) themselves, built from P1 alone.
 D2. The partial order on short games, from the standard recursive test
       G <= H  iff  no G^L >= H and no H^R <= G,
     with equality tested as <= in both directions (no canonicalisation).
 D3. Disjunctive sums G+H and negatives -G, built structurally.
 D4. An INDEPENDENT second decision procedure: outcome classes by normal-play
     backward induction (who wins moving first).  "G = X" is then re-decided as
     "the second player wins G + (-X)", which is the route the paper's own
     proof takes.  Both engines must agree on every value check.
 D5. betahat, R_betahat and the conjecture's predicted atomic weight, computed
     from the integer formulas of P3 for every instance used.
 D6. Verification that each instance really does satisfy the conjecture's
     hypotheses (so it is in scope) and really does violate the hypothesis
     tau < (b-a)/2 of DL Theorem 1.5 (so that theorem is untouched).
 D7. Refutation of the predicted atomic weight WITHOUT computing aw: from
     P6(B,C), aw(G) = P implies that for every j >= 1 the game j*(G - P*up)
     is less than or confused with up and greater than or confused with
     down; we exhibit a j with  up <= j*(G - P*up)  or  j*(G - P*up) <= down,
     which refutes aw(G) = P outright.
 D8. aw(down) = -1 and aw(down+down+*) = -2 by integer arithmetic on P6(A,B).
     Only aw(up) = 1 and aw(*) = 0 are imported numbers; aw(down) is derived
     from aw(-G) = -aw(G), and the premise of that step (that the game the
     paper calls down really is -up) is decided here by the order engine.
     The -2 is obtained twice, by additivity and from the engine-decided
     identity (down+down+*) + 2*up = *, and the two routes must agree.
 D9. A self-test of both engines against 17 textbook CGT relations.
 D10. RULE ANCHORS: the built trees are tested against move-by-move facts the
     paper's own proof states in prose (Right cannot move from a heap <= tau;
     from heap 2t+2 his options are exactly heaps 1..t+1; he can empty a heap
     of 2t+1; Left's only option from heap m is m-1).  Both decision engines
     read the same ts(), so their agreement cannot catch a mis-transcribed
     subtraction set; these anchors can.
 D11. The paper's REMARK, computed rather than cited: for (a,b,tau)=(3,7,2)
     betahat and R_betahat are derived and compared with the paper's 2 and 12,
     and TS_2(n;3,7) is built for 12 <= n <= 34 to confirm that the order test
     refutes the prediction at no such n.  This replaces the Remark's appeal to
     DL Table 4 with a computation.
 D12. Theorem 1 is additionally checked for t = 8..12.
 D13. A sensitivity self-test of the refutation probe of D7: it must refute
     atomic weights that are wrong (aw(down) = 0, aw(down+down+*) = -1) and
     must NOT refute ones that are right (aw(*) = 0, aw(down) = -1).  Without
     this, a "not refuted" outcome anywhere in the transcript would carry no
     information at all.

NOT COVERED (declared, not stubbed): a source-fidelity audit of
arXiv:2607.27989v1 -- checking this program's transcription of Conjecture 7.1,
of Definition 5.1 and of Table 4 against the preprint itself -- needs a network
fetch of a third-party document and is out of scope for a stdlib offline
program.  The program prints a NOTE saying so; it is not counted as a check.  main() also
prints a GAPS section listing every step between the checks and the paper's
claims that no check covers -- foremost that if the transcription of
Conjecture 7.1's scope in P3 is wrong (for instance if a = 1 is inadmissible),
every check here still passes while the paper is wrong.

Exact arithmetic only: the whole computation is on integers and finite game
trees.  There is no floating point anywhere in this file.

Python 3.9+, standard library only.
"""

import sys

# ---------------------------------------------------------------------------
# P4/P5: the paper's own figures.  These appear ONLY on the right-hand side of
# comparisons; nothing below is seeded from them.
# ---------------------------------------------------------------------------

# Theorem 1 of the paper: claimed values of H_n = TS_t(n;1,2t+1).
PAPER_THEOREM1 = {2: "0", 3: "*", 4: "down"}      # keyed by n - (2t)
PAPER_T_RANGE = range(1, 8)                        # t = 1..7

PAPER_BETAHAT = 2                                  # paper: betahat = 2
PAPER_R_FORMULA = lambda t: 2 * t + 4              # paper: R_betahat = 2t+4
PAPER_PREDICTED_AW_AT_R = 0                        # paper: formula predicts 0
PAPER_AW_DOWN = -1                                 # paper: aw(down) = -1

# Additional instances beyond the paper's statements.
PAPER_STRENGTHEN_T_RANGE = range(1, 5)             # t = 1..4
PAPER_STRENGTHEN_VALUE = "down+down+*"             # TS_t(2t+5;1,2t+1)
PAPER_STRENGTHEN_AW = -2                           # aw of that value
PAPER_STRICT_INSIDE = [(4, 2), (5, 3), (7, 4), (9, 6)]   # (b, tau), a = 1
PAPER_STRICT_INSIDE_VALUE = "down"

sys.setrecursionlimit(200000)


# ===========================================================================
# D1/D3: short partizan games as interned option-set pairs.
# ===========================================================================

class Game(object):
    """A short partizan game: a pair of finite sets of Games.

    Instances are interned on their structure, so structurally identical trees
    share one object and one integer id.  That makes the memo tables below
    exact dictionaries keyed by (id, id) with no hashing of deep trees.
    Structural identity is only ever used as a *sufficient* condition; game
    equality is always decided by the order test in leq().
    """

    __slots__ = ("L", "R", "gid", "key")

    def __init__(self, L, R, gid, key):
        self.L = L          # tuple of Games, Left's options
        self.R = R          # tuple of Games, Right's options
        self.gid = gid
        self.key = key

    def __repr__(self):
        return "{%s|%s}" % (
            ",".join(repr(x) for x in self.L),
            ",".join(repr(x) for x in self.R),
        )


_INTERN = {}
_ALL_GAMES = []


def game(L, R):
    """Intern the game with Left options L and Right options R (iterables)."""
    lk = tuple(sorted(set(g.gid for g in L)))
    rk = tuple(sorted(set(g.gid for g in R)))
    key = (lk, rk)
    got = _INTERN.get(key)
    if got is not None:
        return got
    byid = {}
    for g in list(L) + list(R):
        byid[g.gid] = g
    g = Game(tuple(byid[i] for i in lk), tuple(byid[i] for i in rk),
             len(_ALL_GAMES), key)
    _INTERN[key] = g
    _ALL_GAMES.append(g)
    return g


# ===========================================================================
# D2: the standard recursive order test.  G <= H iff no G^L >= H and no H^R <= G.
# No canonicalisation: equality is "<= in both directions".
# ===========================================================================

_LEQ = {}


def leq(g, h):
    """Decide G <= H for short partizan games, exactly."""
    k = (g.gid, h.gid)
    got = _LEQ.get(k)
    if got is not None:
        return got
    # Recursion terminates because every recursive call strictly decreases the
    # sum of the two arguments' formal birthdays; no memo cell is ever read
    # while it is being computed, so no placeholder value is needed.
    res = True
    for gl in g.L:
        if leq(h, gl):          # some G^L >= H
            res = False
            break
    if res:
        for hr in h.R:
            if leq(hr, g):      # some H^R <= G
                res = False
                break
    _LEQ[k] = res
    return res


def equal(g, h):
    """Game equality: G = H iff G <= H and H <= G."""
    return leq(g, h) and leq(h, g)


def less(g, h):
    return leq(g, h) and not leq(h, g)


def greater(g, h):
    return less(h, g)


def confused(g, h):
    """G || H : neither G <= H nor H <= G."""
    return (not leq(g, h)) and (not leq(h, g))


# ===========================================================================
# D3: negative and disjunctive sum, built structurally (memoised, so the sum
# of two DAGs stays a DAG instead of blowing up into a tree).
# ===========================================================================

_NEG = {}
_ADD = {}


def neg(g):
    """-G : swap the roles of the players everywhere."""
    got = _NEG.get(g.gid)
    if got is not None:
        return got
    out = game([neg(r) for r in g.R], [neg(l) for l in g.L])
    _NEG[g.gid] = out
    return out


def add(g, h):
    """G + H : move in exactly one component."""
    k = (g.gid, h.gid)
    got = _ADD.get(k)
    if got is not None:
        return got
    L = [add(gl, h) for gl in g.L] + [add(g, hl) for hl in h.L]
    R = [add(gr, h) for gr in g.R] + [add(g, hr) for hr in h.R]
    out = game(L, R)
    _ADD[k] = out
    return out


def sub(g, h):
    return add(g, neg(h))


# ===========================================================================
# P2: the reference values, written down exactly as the paper defines them:
#   0 = {|},  * = {0|0},  up = {0|*},  down = {*|0}.
# down is built from its own definition, NOT as -up; that identity is then
# checked in the engine self-test.
# ===========================================================================

ZERO = game([], [])
STAR = game([ZERO], [ZERO])
UP = game([ZERO], [STAR])
DOWN = game([STAR], [ZERO])


def integer(n):
    """The game n: {n-1 | } for n > 0, { | n+1} for n < 0."""
    if n == 0:
        return ZERO
    if n > 0:
        return game([integer(n - 1)], [])
    return game([], [integer(n + 1)])


def up_multiple(k):
    """k * up, for any integer k (negative k gives |k| copies of down)."""
    if k == 0:
        return ZERO
    base = UP if k > 0 else DOWN
    out = base
    for _ in range(abs(k) - 1):
        out = add(out, base)
    return out


def reference_table():
    """Named values used to *report* derived game values.  Reporting only:
    every PASS/FAIL below compares against a value the paper names."""
    tab = [("0", ZERO), ("*", STAR), ("up", UP), ("down", DOWN),
           ("up+*", add(UP, STAR)), ("down+*", add(DOWN, STAR))]
    for k in (2, 3, 4):
        tab.append(("%d*up" % k, up_multiple(k)))
        tab.append(("%d*up+*" % k, add(up_multiple(k), STAR)))
        tab.append(("%d*down" % k, up_multiple(-k)))
        tab.append(("%d*down+*" % k, add(up_multiple(-k), STAR)))
    for k in (1, 2, -1, -2):
        tab.append(("%d" % k, integer(k)))
    return tab


def name_of(g, tab=None):
    """Report a derived game value by matching it against the table, or say
    that it is outside the table.  Never used to decide a check."""
    if tab is None:
        tab = reference_table()
    for nm, ref in tab:
        if equal(g, ref):
            return nm
    return "<not in reference table>"


# ===========================================================================
# D1: the game trees of TS_tau(n; a, b), straight from P1.
#     TS(0) = 0;  Left  moves n -> n-k for 1     <= k <= a with n-k >= 0;
#                 Right moves n -> n-k for tau+1 <= k <= b with n-k >= 0.
# Built bottom-up so the construction needs no deep recursion.
# ===========================================================================

_TS = {}


def ts(n, tau, a, b):
    """TS_tau(n; a, b) as an interned game."""
    if n < 0:
        raise ValueError("heap size must be non-negative")
    if not (0 < a <= b):
        raise ValueError("rules require 0 < a <= b")
    if tau < 0:
        raise ValueError("truncation must be non-negative")
    for m in range(0, n + 1):
        key = (m, tau, a, b)
        if key in _TS:
            continue
        L = [_TS[(m - k, tau, a, b)] for k in range(1, a + 1) if m - k >= 0]
        R = [_TS[(m - k, tau, a, b)] for k in range(tau + 1, b + 1)
             if m - k >= 0]
        _TS[key] = game(L, R)
    return _TS[(n, tau, a, b)]


def ts_option_report(n, tau, a, b):
    """Human-readable statement of the option sets actually generated, so a
    reader can check the tree against Definition 5.1 by eye."""
    lo = [n - k for k in range(1, a + 1) if n - k >= 0]
    ro = [n - k for k in range(tau + 1, b + 1) if n - k >= 0]
    return "TS_%d(%d;%d,%d): L-options -> heaps %s ; R-options -> heaps %s" % (
        tau, n, a, b, lo if lo else "none", ro if ro else "none")


# ===========================================================================
# D4: a SECOND, independent decision procedure -- normal-play backward
# induction.  It never consults leq().  The paper's own proof works this way:
# it shows H_{2t+3} + * and H_{2t+4} + up are second-player wins.  We use the
# standard fact quoted in the paper: a short game equals 0 exactly when the
# second player wins, whichever player that is.
# ===========================================================================

_WINS = {}


def wins_moving_first(g, player):
    """True iff `player` ('L' or 'R'), moving first in g, wins under normal
    play (a player with no legal move loses)."""
    k = (g.gid, player)
    got = _WINS.get(k)
    if got is not None:
        return got
    opts = g.L if player == "L" else g.R
    other = "R" if player == "L" else "L"
    res = False
    for o in opts:
        if not wins_moving_first(o, other):
            res = True
            break
    _WINS[k] = res
    return res


def second_player_wins(g):
    """g = 0, decided by outcome class alone."""
    return (not wins_moving_first(g, "L")) and (not wins_moving_first(g, "R"))


def outcome_class(g):
    """P (second player), N (first player), L (Left) or R (Right) wins."""
    l = wins_moving_first(g, "L")
    r = wins_moving_first(g, "R")
    if l and r:
        return "N"
    if l and not r:
        return "L"
    if r and not l:
        return "R"
    return "P"


def equal_via_outcome(g, h):
    """g = h re-decided independently: second player wins g + (-h)."""
    return second_player_wins(sub(g, h))


# ===========================================================================
# Check harness.
# ===========================================================================

RESULTS = []


def check(ok, label, detail=""):
    ok = bool(ok)
    RESULTS.append((ok, label))
    line = ("PASS " if ok else "FAIL ") + label
    if detail:
        line += "   [" + detail + "]"
    print(line)
    return ok


def note(text):
    print("      " + text)


# ===========================================================================
# D9: self-test of both engines against textbook CGT relations.  If any of
# these fails, nothing else in the transcript means anything.
# ===========================================================================

def engine_selftest():
    two_up = add(UP, UP)
    rel = [
        ("0 < up", less(ZERO, UP)),
        ("down < 0", less(DOWN, ZERO)),
        ("* is confused with 0", confused(STAR, ZERO)),
        ("up + down = 0", equal(add(UP, DOWN), ZERO)),
        ("-{0|*} equals {*|0}, i.e. -up = down", equal(neg(UP), DOWN)),
        ("* + * = 0", equal(add(STAR, STAR), ZERO)),
        # The next three relations are CONFUSIONS, not inequalities: up+* (= up-star) is confused
        # with 0, * is incomparable with up, and down is incomparable with *. Any strict-inequality
        # reading of these three is false, which is exactly why the atomic-weight criterion in
        # refute_aw_equals() below is stated as "less than or confused with" rather than "<".
        ("up + * is confused with 0",
         not greater(add(UP, STAR), ZERO) and not less(add(UP, STAR), ZERO)
         and not equal(add(UP, STAR), ZERO)),
        ("* is confused with up", not less(STAR, UP) and not greater(STAR, UP)
         and not equal(STAR, UP)),
        ("down is confused with *", not less(DOWN, STAR) and not greater(DOWN, STAR)
         and not equal(DOWN, STAR)),
        ("up + up > up", greater(two_up, UP)),
        ("up + up > *", greater(two_up, STAR)),
        ("1 > up  (integer dominates the infinitesimal)",
         greater(integer(1), UP)),
        ("1 + (-1) = 0", equal(add(integer(1), integer(-1)), ZERO)),
        ("-(up+up+*) = down+down+*",
         equal(neg(add(two_up, STAR)), add(up_multiple(-2), STAR))),
        ("outcome classes: 0 is P, * is N, up is L, down is R",
         (outcome_class(ZERO), outcome_class(STAR), outcome_class(UP),
          outcome_class(DOWN)) == ("P", "N", "L", "R")),
        ("both engines agree that up + down = 0",
         equal_via_outcome(add(UP, DOWN), ZERO)),
        ("both engines agree that up + * is NOT 0",
         not equal_via_outcome(add(UP, STAR), ZERO)),
    ]
    bad = [lab for lab, ok in rel if not ok]
    for lab, ok in rel:
        note(("ok   " if ok else "BAD  ") + lab)
    check(not bad, "engine self-test: %d/%d textbook relations hold"
          % (len(rel) - len(bad), len(rel)),
          "failures: " + "; ".join(bad) if bad else "no floating point used")
    return not bad


# ===========================================================================
# D5/D6: the conjecture's own integer arithmetic, from P3.  All exact integer
# arithmetic; Python's // is floor division, which is the floor() of the paper.
# ===========================================================================

def conjecture_data(a, b, tau, n):
    """betahat, R_betahat and the predicted atomic weight at heap n, computed
    from the formulas of P3 -- never seeded."""
    if b - a - tau <= 0:
        raise ValueError("betahat undefined unless tau < b-a")
    betahat = (b - a) // (b - a - tau)
    R = betahat * (a + tau + 1)
    pred = None
    if n >= R:
        pred = -((n - R) // a)
    d = {
        "a": a, "b": b, "tau": tau, "n": n,
        "betahat": betahat, "R": R, "pred": pred,
        "hyp_0_lt_a_le_b": (0 < a <= b),
        "hyp_a_lt_b": (a < b),
        "hyp_tau_lt_b_minus_a": (tau < b - a),
        "hyp_n_ge_R": (n >= R),
        # DL Theorem 1.5 assumes 1 <= tau < (b-a)/2; kept in integers as
        # 2*tau < b-a.  The counterexamples must FAIL this.
        "thm15_applies": (1 <= tau and 2 * tau < b - a),
        "tau_eq_half": (2 * tau == b - a),
        "strictly_inside": (b - a < 2 * tau < 2 * (b - a)),
    }
    return d


def describe_instance(d):
    return ("a=%d b=%d tau=%d n=%d -> betahat=floor(%d/%d)=%d, "
            "R_betahat=%d*(%d+%d+1)=%d, predicted aw=%s"
            % (d["a"], d["b"], d["tau"], d["n"], d["b"] - d["a"],
               d["b"] - d["a"] - d["tau"], d["betahat"], d["betahat"],
               d["a"], d["tau"], d["R"],
               "n/a" if d["pred"] is None else str(d["pred"])))


def in_scope(d):
    """Every standing hypothesis of the conjecture, as derived booleans."""
    return (d["hyp_0_lt_a_le_b"] and d["hyp_a_lt_b"]
            and d["hyp_tau_lt_b_minus_a"] and d["hyp_n_ge_R"])


# ===========================================================================
# D7: refute "aw(G) = P" using only the order relation.
#
# From P6(B) additivity, aw(G) = P iff aw(G - P*up) = 0; from P6(C), aw(H) = 0
# requires that for EVERY integer j >= 1 the game j*H be less than or confused
# with up and greater than or confused with down.  So exhibiting a single j with
#     up <= j*H   or   j*H <= down
# refutes aw(G) = P outright.  This is a decidable order computation: no atomic
# weight is ever computed.
# ===========================================================================

# The sweep runs j = 1, 2. Each j is an independent sufficient condition for refutation, so a bound
# on j only ever weakens a refutation attempt; it can never turn a failure into a pass. The node cap
# bounds the sums formed for j > 1 within the memory budget of this program, and the detail line
# names every j that was skipped rather than passing over it silently.
AW_J_MAX = 2
AW_J_NODE_CAP = 20000


def refute_aw_equals(g, P):
    """Return (refuted, detail, probed) for the assertion aw(g) = P, where
    probed is the tuple of j actually swept.

    The decision uses only leq(); the name of H is cosmetic, so it is only
    computed for small H (matching a big sum against the whole reference table
    is the one superlinear step in this program).
    """
    # The criterion for aw(H) = 0 is that H be infinitesimal with respect to up: for every integer
    # j >= 1, j*H is less than OR CONFUSED WITH up, and greater than or confused with down. The
    # refutable half at a given j is therefore
    #     refuted  iff  up <= j*H  or  j*H <= down
    # and NOT the strict pair down < j*H < up. Many all-small games are CONFUSED with up rather than
    # below it -- aw(*) = 0 while * is confused with up and with down -- so a strict form would
    # reject correct atomic weights. Concretely, the strict form rejects the weights PROVEN in
    # Dargad-Larsson's Theorem 1.5 at (a,b,tau) = (3,7,1), at exactly the heap sizes where the proven
    # weight is odd, so it is unsound as a refutation criterion.
    # The criterion quantifies over ALL positive j, and testing only j = 1 would be weaker than the
    # criterion: for the assertion aw(down+down+*) = -1 the game H = down+* satisfies both halves at
    # j = 1, while 2H = down+down (the stars cancel) is <= down and refutes. Hence the sweep over
    # small j. Each j is a sound refutation on its own; more j can only strengthen.
    H = sub(g, up_multiple(P))
    nH = count_nodes(H)
    nm = name_of(H) if nH <= 64 else "<%d positions, not named>" % nH
    refuted, why, jH = False, [], H
    probed = []
    for j in range(1, AW_J_MAX + 1):
        if j > 1:
            jH = add(jH, H)
            if count_nodes(jH) > AW_J_NODE_CAP:      # keep the sweep bounded; report the cut
                why.append("j=%d skipped (too large)" % j)
                break
        up_le, le_dn = leq(UP, jH), leq(jH, DOWN)
        probed.append(j)
        why.append("j=%d: up<=jH %s, jH<=down %s" % (j, up_le, le_dn))
        if up_le or le_dn:
            refuted = True
            break
    detail = "H := G - (%d)*up = %s ; %s" % (P, nm, " ; ".join(why))
    return refuted, detail, tuple(probed)


# ===========================================================================
# D13: sensitivity of the probe just defined.  Several checks below report that
# the probe found NO refutation -- most importantly the check on the paper's
# Remark, where the derived games are large and both order comparisons come back
# False at every j.  Such a line means nothing unless the probe is known to fire
# when it should.  So run it on four assertions whose truth is settled by the
# imported facts P6(A),(B) plus the engine self-test: two false weights, which it
# must refute, and two true ones, which it must leave standing.
# ===========================================================================

def check_probe_sensitivity():
    print("")
    print("--- Sensitivity of the atomic-weight refutation probe ---")
    down_down_star = add(up_multiple(-2), STAR)
    cases = [
        # (assertion, game, claimed weight, must the probe refute it?)
        ("aw(down) = 0", DOWN, 0, True),
        ("aw(down+down+*) = -1", down_down_star, -1, True),
        ("aw(*) = 0", STAR, 0, False),
        ("aw(down) = -1", DOWN, -1, False),
    ]
    ok = []
    for lab, g, P, want in cases:
        rf, detail, js = refute_aw_equals(g, P)
        ok.append(rf == want)
        note("%-22s : this weight is %s, so the probe must%s refute it ; probe "
             "says %-5s ; j swept %s  [%s]"
             % (lab, "FALSE" if want else "TRUE", "" if want else " NOT",
                rf, list(js), detail))
    check(all(ok),
          "sensitivity of the refutation probe: it refutes the two FALSE "
          "assertions aw(down)=0 and aw(down+down+*)=-1 (the second only at "
          "j=2, so the sweep over j is load-bearing) and does NOT refute the "
          "two TRUE ones aw(*)=0 and aw(down)=-1, so a 'not refuted' line "
          "elsewhere in this transcript is the verdict of a probe that can "
          "actually fire",
          "%d/%d cases as expected" % (sum(1 for x in ok if x), len(ok)))
    return all(ok)


# ===========================================================================
# D8: atomic weight as integer arithmetic on the cited facts P6(A) and P6(B).
# The ONLY imported numbers are the two cited primitives aw(up) = 1 and
# aw(*) = 0.  aw(down) is NOT written down: it is derived from the cited rule
# aw(-G) = -aw(G), and the premise of that rule -- that the game the paper
# calls down really is -up -- is decided here by the order engine on the
# literal definitions of P2, not assumed.  aw(0) is derived the same way from
# additivity and the engine-verified identity up + down = 0.
# ===========================================================================

AW_UP = 1          # P6(A), cited: aw(up) = 1
AW_STAR = 0        # P6(A), cited: aw(*) = 0


def derive_aw_atoms():
    """Return (atoms, working_lines, engine_ok).

    atoms maps an atom name to its atomic weight.  Only AW_UP and AW_STAR are
    imported; the rest are derived, and every derivation names the order-engine
    identity that licenses it.  engine_ok is False if any of those identities
    fails, in which case the derivation is void.
    """
    lines, ok = [], []
    id_down = equal(DOWN, neg(UP))            # {*|0} = -{0|*} ?
    aw_down = -AW_UP                          # P6(A): aw(-G) = -aw(G)
    ok.append(id_down)
    lines.append("engine decides  down = -up  : %s   =>  aw(down) = -aw(up) = %d"
                 % (id_down, aw_down))
    id_zero = equal(add(UP, DOWN), ZERO)      # up + down = 0 ?
    aw_zero = AW_UP + aw_down                 # P6(B): additivity
    ok.append(id_zero)
    lines.append("engine decides  up + down = 0 : %s   =>  aw(0) = aw(up) + "
                 "aw(down) = %d" % (id_zero, aw_zero))
    atoms = {"up": AW_UP, "*": AW_STAR, "down": aw_down, "0": aw_zero}
    return atoms, lines, all(ok)


def aw_of_decomposition(atoms, names):
    """aw of a formal sum of atoms, by additivity.  Returns (value, working)."""
    total = 0
    parts = []
    for at in names:
        if at not in atoms:
            raise KeyError("no cited atomic weight for atom %r" % (at,))
        total += atoms[at]
        parts.append("aw(%s)=%d" % (at, atoms[at]))
    working = " + ".join(parts) + " = %d" % total
    return total, working


# ===========================================================================
# Reporting one instance: derive the value two independent ways and print the
# intermediate quantities.  Returns (order_engine_says_equal,
#                                    outcome_engine_says_equal).
# ===========================================================================

def report_value(n, tau, a, b, ref, refname, tab):
    g = ts(n, tau, a, b)
    e_order = equal(g, ref)
    e_outcome = equal_via_outcome(g, ref)
    derived = name_of(g, tab)
    note("  " + ts_option_report(n, tau, a, b))
    note("    tree nodes reachable: %d ; outcome class: %s ; derived value: %s"
         % (count_nodes(g), outcome_class(g), derived))
    note("    = %s ? order test: %s ; independent outcome test on TS-(%s): %s"
         % (refname, e_order, refname, e_outcome))
    return e_order, e_outcome


def count_nodes(g):
    """Number of distinct positions reachable from g (a derived statistic)."""
    seen = set()
    stack = [g]
    while stack:
        x = stack.pop()
        if x.gid in seen:
            continue
        seen.add(x.gid)
        stack.extend(x.L)
        stack.extend(x.R)
    return len(seen)


# ===========================================================================
# Rule anchors.  Both decision engines share ts(), so a mis-transcribed
# subtraction set would be agreed on by both and the "engines agree" checks
# would give no protection.  These anchors compare the BUILT TREES against
# facts the paper's own proof states in prose:
#   "Left removes one token, while Right removes any number from t+1 through
#    2t+1";  "If Left starts, she leaves 2t+1 tokens and Right removes the
#    entire heap";  "If Right starts, he leaves a heap of size m in
#    {1,...,t+1}; Left moves to m-1 <= t, after which Right has no legal move".
# A wrong endpoint in Right's set (tau instead of tau+1, or b-1 instead of b)
# fails one of these.
# ===========================================================================

def check_definition_anchors():
    print("")
    print("--- Rule anchors: the built trees against the paper's own prose ---")
    small, leftopt, r_low, r_set, distinct = [], [], [], [], []
    for t in PAPER_T_RANGE:
        a, b, tau = 1, 2 * t + 1, t
        ts(2 * t + 5, tau, a, b)          # build the whole segment once
        small.append(all((not ts(m, tau, a, b).R)
                         and equal(ts(m, tau, a, b), integer(m))
                         for m in range(0, t + 1)))
        leftopt.append(all(tuple(x.gid for x in ts(m, tau, a, b).L)
                           == (ts(m - 1, tau, a, b).gid,)
                           for m in range(1, 2 * t + 6)))
        r_low.append(ts(0, tau, a, b).gid
                     in set(x.gid for x in ts(2 * t + 1, tau, a, b).R))
        got = set(x.gid for x in ts(2 * t + 2, tau, a, b).R)
        want = set(ts(m, tau, a, b).gid for m in range(1, t + 2))
        r_set.append(got == want)
        distinct.append(len(set(ts(m, tau, a, b).gid
                                for m in range(0, t + 3))) == t + 3)
        note("t=%d : heaps 0..%d have no Right option and equal the integer: "
             "%s ; R-options of heap %d are exactly heaps 1..%d: %s ; Right "
             "can empty a heap of %d: %s" % (t, t, small[-1], 2 * t + 2,
                                             t + 1, r_set[-1], 2 * t + 1,
                                             r_low[-1]))
    check(all(small), "rule anchor: Right has no move from a heap of size "
                      "<= tau and the value there is the integer n, as the "
                      "paper's proof asserts (pins Right's lower endpoint "
                      "tau+1)")
    check(all(r_set) and all(r_low) and all(distinct),
          "rule anchor: Right's options from heap 2t+2 are exactly heaps "
          "1..t+1 and Right can empty a heap of 2t+1 (pins both endpoints of "
          "{tau+1,...,b}; the heaps 0..tau+2 are pairwise structurally "
          "distinct, so the set comparison is not an interning artefact)")
    check(all(leftopt), "rule anchor: Left's only option from heap m >= 1 is "
                        "heap m-1, i.e. a = 1 (pins Left's subtraction set)")
    return all(small) and all(r_set) and all(r_low) and all(leftopt)


# ===========================================================================
# The paper's Theorem 1.
#   TS_t(2t+2;1,2t+1) = 0, TS_t(2t+3;1,2t+1) = *, TS_t(2t+4;1,2t+1) = down,
#   for t = 1..7.
# ===========================================================================

def check_theorem1(tab):
    refs = {2: (ZERO, "0"), 3: (STAR, "*"), 4: (DOWN, "down")}
    hits = {2: [], 3: [], 4: []}
    cross = []
    print("")
    print("--- Theorem 1: TS_t(n;1,2t+1) for n = 2t+2, 2t+3, 2t+4 ---")
    for t in PAPER_T_RANGE:
        a, b, tau = 1, 2 * t + 1, t
        note("t=%d : a=%d, b=2t+1=%d, tau=t=%d" % (t, a, b, tau))
        for off in (2, 3, 4):
            ref, refname = refs[off]
            if PAPER_THEOREM1[off] != refname:
                raise AssertionError("reference table disagrees with P4")
            eo, ec = report_value(2 * t + off, tau, a, b, ref, refname, tab)
            hits[off].append(eo)
            cross.append(eo == ec)
    ts_all = list(PAPER_T_RANGE)
    for off in (2, 3, 4):
        nm = PAPER_THEOREM1[off]
        bad = [t for t, ok in zip(ts_all, hits[off]) if not ok]
        check(not bad,
              "Theorem 1: TS_t(2t+%d;1,2t+1) = %s for all t in %d..%d"
              % (off, nm, ts_all[0], ts_all[-1]),
              "%d/%d instances agree%s" % (len(ts_all) - len(bad), len(ts_all),
                                           "" if not bad
                                           else "; failed at t=" + str(bad)))
    check(all(cross),
          "the two independent engines (order test, outcome test) agree on "
          "all %d Theorem 1 value decisions" % len(cross),
          "%d/%d agree" % (sum(1 for c in cross if c), len(cross)))
    return all(hits[2]) and all(hits[3]) and all(hits[4]) and all(cross)


# ===========================================================================
# The integer arithmetic of the Corollary, derived per t.
# ===========================================================================

def check_corollary_arithmetic():
    print("")
    print("--- betahat, R_betahat and the predicted atomic weight ---")
    bh, rr, pr, scope, t15, half = [], [], [], [], [], []
    for t in PAPER_T_RANGE:
        d = conjecture_data(1, 2 * t + 1, t, 2 * t + 4)
        note("t=%d : %s" % (t, describe_instance(d)))
        note("    paper says betahat=%d, R_betahat=%d, prediction=%d"
             % (PAPER_BETAHAT, PAPER_R_FORMULA(t), PAPER_PREDICTED_AW_AT_R))
        note("    in scope: 0<a<=b %s, a<b %s, tau<b-a %s, n>=R %s ; "
             "n == R_betahat: %s"
             % (d["hyp_0_lt_a_le_b"], d["hyp_a_lt_b"],
                d["hyp_tau_lt_b_minus_a"], d["hyp_n_ge_R"],
                d["n"] == d["R"]))
        note("    DL Thm 1.5 hypothesis 2*tau < b-a : %s (2*tau=%d, b-a=%d) ; "
             "tau = (b-a)/2 : %s"
             % (d["thm15_applies"], 2 * d["tau"], d["b"] - d["a"],
                d["tau_eq_half"]))
        bh.append(d["betahat"] == PAPER_BETAHAT)
        rr.append(d["R"] == PAPER_R_FORMULA(t) == d["n"])
        pr.append(d["pred"] == PAPER_PREDICTED_AW_AT_R)
        scope.append(in_scope(d))
        t15.append(not d["thm15_applies"])
        half.append(d["tau_eq_half"])
    span = "t=%d..%d" % (min(PAPER_T_RANGE), max(PAPER_T_RANGE))
    check(all(bh), "derived betahat = floor((b-a)/(b-a-tau)) equals the "
                   "paper's 2 for " + span)
    check(all(rr), "derived R_betahat = betahat*(a+tau+1) equals the paper's "
                   "2t+4, and equals n, for " + span)
    check(all(pr), "derived prediction -floor((n-R_betahat)/a) equals the "
                   "paper's 0 for " + span)
    check(all(scope), "every instance satisfies the conjecture's standing "
                      "hypotheses (0<a<=b, a<b, tau<b-a, n>=R_betahat)")
    check(all(t15), "no instance satisfies DL Theorem 1.5's hypothesis "
                    "tau < (b-a)/2, so that theorem is untouched")
    check(all(half), "every instance sits exactly at tau = (b-a)/2, as the "
                     "paper's Remark states")
    return all(bh) and all(rr) and all(pr) and all(scope) and all(t15)


# ===========================================================================
# The refutation.  Two routes, both reported:
#   (R1) order-theoretic, decidable here: for some j >= 1 the game
#        j*(G - pred*up) satisfies up <= it or it <= down, so aw(G) != pred by
#        P6(B)+P6(C);
#   (R2) the paper's own route: aw(G) = aw(down) = -1 by P6(A), != pred.
# ===========================================================================

def check_refutation(tab):
    print("")
    print("--- The refutation at n = R_betahat ---")
    ref_ok, val_ok, preds = [], [], []
    for t in PAPER_T_RANGE:
        a, b, tau, n = 1, 2 * t + 1, t, 2 * t + 4
        d = conjecture_data(a, b, tau, n)
        g = ts(n, tau, a, b)
        is_down = equal(g, DOWN) and equal_via_outcome(g, DOWN)
        val_ok.append(is_down)
        preds.append(d["pred"])
        refuted, detail, _js = refute_aw_equals(g, d["pred"])
        ref_ok.append(refuted)
        note("t=%d : TS_%d(%d;%d,%d) derived value %s ; conjecture predicts "
             "aw = %d" % (t, tau, n, a, b, name_of(g, tab), d["pred"]))
        note("    %s  -> aw != %d : %s" % (detail, d["pred"], refuted))
    atoms, awlines, atoms_ok = derive_aw_atoms()
    for ln in awlines:
        note("aw derivation: " + ln)
    aw_down, working = aw_of_decomposition(atoms, ["down"])
    note("aw of the derived value, by additivity on the cited primitives: "
         + working)
    check(all(val_ok),
          "the game at n = R_betahat is down (both engines, %s)"
          % ("t=%d..%d" % (min(PAPER_T_RANGE), max(PAPER_T_RANGE))))
    check(all(ref_ok),
          "aw(TS_t(R_betahat;1,2t+1)) != 0 by the order test alone "
          "(no atomic weight computed), refuting Conjecture 7.1")
    check(atoms_ok,
          "the identities licensing aw(down) = -aw(up) are decided by the "
          "order engine, not assumed (down = -up and up + down = 0)")
    check(aw_down == PAPER_AW_DOWN,
          "aw(down), DERIVED as -aw(up) from the engine-verified identity "
          "down = -up, equals the paper's -1",
          "derived %d, paper %d" % (aw_down, PAPER_AW_DOWN))
    # The Corollary itself: for EVERY t, the derived value's atomic weight
    # differs from the prediction DERIVED from the conjecture's own formula.
    # Both operands vary with the computation; neither is a literal read back.
    contra = [(v and (aw_down != p)) for v, p in zip(val_ok, preds)]
    check(all(contra) and atoms_ok,
          "CONTRADICTION (the Corollary): for every t the derived value is "
          "down with derived aw = %d, while the prediction derived from "
          "-floor((n-R_betahat)/a) is %s" % (aw_down, sorted(set(preds))),
          "%d/%d instances contradict the formula" % (sum(1 for c in contra
                                                          if c), len(contra)))
    return all(val_ok) and all(ref_ok) and all(contra)


# ===========================================================================
# One step past R_betahat, beyond the paper's claim:
#   TS_t(2t+5;1,2t+1) = down+down+*  for t = 1..4, so aw = -2 where the
#   conjecture predicts -1.  This closes the "R_betahat off by one" door.
# ===========================================================================

def check_strengthening_next_heap(tab):
    print("")
    print("--- One heap past R_betahat, n = 2t+5 ---")
    target = add(up_multiple(-2), STAR)
    note("reference target built from P2: down+down+* = %s (nodes %d)"
         % (name_of(target, tab), count_nodes(target)))
    val, cross, preds, ref_ok = [], [], [], []
    for t in PAPER_STRENGTHEN_T_RANGE:
        a, b, tau, n = 1, 2 * t + 1, t, 2 * t + 5
        d = conjecture_data(a, b, tau, n)
        note("t=%d : %s (n = R_betahat + %d)" % (t, describe_instance(d),
                                                 n - d["R"]))
        eo, ec = report_value(n, tau, a, b, target, PAPER_STRENGTHEN_VALUE, tab)
        val.append(eo)
        cross.append(eo == ec)
        preds.append(d["pred"] == -1)
        refuted, detail, _js = refute_aw_equals(ts(n, tau, a, b), d["pred"])
        ref_ok.append(refuted)
        note("    %s  -> aw != %d : %s" % (detail, d["pred"], refuted))
    atoms, _awlines, atoms_ok = derive_aw_atoms()
    aw_val, working = aw_of_decomposition(atoms, ["down", "down", "*"])
    note("aw(down+down+*) by additivity: " + working)
    # Second, engine-backed route to the same number: the order engine decides
    # (down+down+*) + 2*up = *, and additivity then forces
    # aw(down+down+*) = aw(*) - 2*aw(up).
    id_ok = equal(add(target, up_multiple(2)), STAR)
    aw_val2 = AW_STAR - 2 * AW_UP
    note("engine decides  (down+down+*) + up + up = *  : %s   =>  "
         "aw(down+down+*) = aw(*) - 2*aw(up) = %d" % (id_ok, aw_val2))
    span = "t=%d..%d" % (min(PAPER_STRENGTHEN_T_RANGE),
                         max(PAPER_STRENGTHEN_T_RANGE))
    check(all(val), "strengthening (beyond the paper's claim): "
                    "TS_t(2t+5;1,2t+1) = down+down+* for " + span,
          "%d/%d instances agree" % (sum(1 for x in val if x), len(val)))
    check(all(cross), "both engines agree on all %d strengthening decisions"
          % len(cross))
    check(all(preds), "derived prediction at n = R_betahat+1 is -1 for " + span)
    check(all(ref_ok), "strengthening (beyond the paper's claim): "
                       "aw(TS_t(2t+5;1,2t+1)) != -1 by the order test alone, "
                       "so the failure persists one heap past R_betahat")
    check(aw_val == PAPER_STRENGTHEN_AW and aw_val2 == PAPER_STRENGTHEN_AW
          and id_ok and atoms_ok,
          "aw(down+down+*) = -2 by TWO routes that agree: additivity on the "
          "derived atoms, and the engine-decided identity (down+down+*) + "
          "2*up = *", "additivity %d, engine route %d, reference %d"
          % (aw_val, aw_val2, PAPER_STRENGTHEN_AW))
    return all(val) and all(cross) and all(preds) and all(ref_ok)


# ===========================================================================
# Away from the boundary, beyond the paper's claim:
#   TS_tau(R_betahat;1,b) = down for (b,tau) in {(4,2),(5,3),(7,4),(9,6)},
#   each strictly inside (b-a)/2 < tau < b-a.  These show the failure is not
#   an artefact of the equality tau = (b-a)/2.
# ===========================================================================

def check_strictly_inside(tab):
    print("")
    print("--- a=1, (b,tau) strictly inside (b-a)/2 < tau < b-a ---")
    val, cross, inside, preds, ref_ok = [], [], [], [], []
    for (b, tau) in PAPER_STRICT_INSIDE:
        a = 1
        # R_betahat does not depend on n, so derive it first (n=0 gives no
        # prediction), then re-derive the full instance at n = R_betahat.
        R = conjecture_data(a, b, tau, 0)["R"]
        d = conjecture_data(a, b, tau, R)
        note("(b,tau)=(%d,%d) : %s" % (b, tau, describe_instance(d)))
        note("    strictly inside: b-a=%d < 2*tau=%d < 2(b-a)=%d : %s"
             % (b - a, 2 * tau, 2 * (b - a), d["strictly_inside"]))
        eo, ec = report_value(R, tau, a, b, DOWN, PAPER_STRICT_INSIDE_VALUE,
                              tab)
        val.append(eo)
        cross.append(eo == ec)
        inside.append(d["strictly_inside"])
        preds.append(d["pred"] == 0)
        refuted, detail, _js = refute_aw_equals(ts(R, tau, a, b), d["pred"])
        ref_ok.append(refuted)
        note("    %s  -> aw != %d : %s" % (detail, d["pred"], refuted))
    n_i = len(PAPER_STRICT_INSIDE)
    check(all(inside), "all %d instances satisfy (b-a)/2 < tau < b-a strictly"
          % n_i)
    check(all(val), "strengthening (beyond the paper's claim): "
                    "TS_tau(R_betahat;1,b) = down for all %d instances" % n_i,
          "%d/%d agree" % (sum(1 for x in val if x), n_i))
    check(all(cross), "both engines agree on all %d strictly-inside decisions"
          % len(cross))
    check(all(preds), "derived prediction at n = R_betahat is 0 for all %d"
          % n_i)
    check(all(ref_ok), "strengthening (beyond the paper's claim): aw != 0 by "
                       "the order test alone for all %d instances, so the "
                       "failure is not confined to tau = (b-a)/2" % n_i)
    return all(val) and all(cross) and all(inside) and all(ref_ok)


# ===========================================================================
# Further values of t.  Theorem 1 is stated for EVERY t >= 1 and no
# induction on t is verified anywhere in this program; the only thing a
# computation can do about a universal quantifier is widen the finite window.
# This is cheap (heap 2t+4 has at most 2t+5 positions), so we widen it, and we
# say plainly in the label that the window is still finite.
# ===========================================================================

EXTENDED_T_RANGE = range(8, 13)


def check_extended_range():
    print("")
    print("--- Further values of t: Theorem 1 for t = 8..12 ---")
    ok = []
    for t in EXTENDED_T_RANGE:
        a, b, tau = 1, 2 * t + 1, t
        g2 = ts(2 * t + 2, tau, a, b)
        g3 = ts(2 * t + 3, tau, a, b)
        g4 = ts(2 * t + 4, tau, a, b)
        this = (equal(g2, ZERO) and equal_via_outcome(g2, ZERO)
                and equal(g3, STAR) and equal_via_outcome(g3, STAR)
                and equal(g4, DOWN) and equal_via_outcome(g4, DOWN))
        ok.append(this)
        note("t=%2d : H_%d = 0, H_%d = *, H_%d = down, both engines : %s"
             % (t, 2 * t + 2, 2 * t + 3, 2 * t + 4, this))
    check(all(ok),
          "Theorem 1 also holds for t = %d..%d "
          "(the window is still FINITE: no induction on t is verified anywhere "
          "in this program)"
          % (min(EXTENDED_T_RANGE), max(EXTENDED_T_RANGE)),
          "%d/%d instances agree" % (sum(1 for x in ok if x), len(ok)))
    return all(ok)


# ===========================================================================
# The paper's Remark.  It claims that tau = (b-a)/2 is not itself the
# obstruction, because for (a,b) = (3,7) and tau = 2 = (b-a)/2 one has
# betahat = 2, R_betahat = 12, and the atomic weights tabulated in DL Table 4
# agree with (1) for every tabulated n with 12 <= n <= 34.  That claim carries
# the abstract's assertion "the equality tau=(b-a)/2 ... is not itself the
# obstruction", so it is a claim of the paper, not a citation we may skip.
# The first half is computable offline: we DERIVE betahat and R_betahat.  The
# second half is only PARTLY computable: we DERIVE the games TS_2(n;3,7) and ask
# whether the order test refutes the prediction, which can contradict the Remark
# but can never confirm it, and which does not read Table 4 at all.  A
# refutation anywhere in 12 <= n <= 34 would contradict the Remark; the absence
# of one is weaker evidence, and the check below says so in its own label.
# ===========================================================================

def check_remark_boundary(tab):
    print("")
    print("--- The paper's Remark: boundary instance (a,b,tau) = (3,7,2) ---")
    a, b, tau = 3, 7, 2
    heaps = range(12, 35)
    d0 = conjecture_data(a, b, tau, 12)
    note("derived: %s" % describe_instance(d0))
    note("paper's Remark asserts betahat=2 and R_betahat=12 ; "
         "2*tau == b-a (the boundary): %s" % (2 * tau == b - a))
    refuted_at, examined, cut = [], [], []
    for n in heaps:
        g = ts(n, tau, a, b)
        pred = conjecture_data(a, b, tau, n)["pred"]
        rf, detail, js = refute_aw_equals(g, pred)
        examined.append(n)
        if len(js) < AW_J_MAX:
            cut.append(n)
        if rf:
            refuted_at.append(n)
        note("  n=%2d  value %-24s predicted aw=%d  refuted: %-5s  [%s]"
             % (n, name_of(g, tab), pred, rf, detail))
    note("probe coverage here: %d of %d heap sizes were swept over the whole "
         "range j = 1..%d ; %d were swept at j = 1 only, the larger j being cut "
         "by the %d-node cap (n = %s)"
         % (len(examined) - len(cut), len(examined), AW_J_MAX, len(cut),
            AW_J_NODE_CAP, cut if cut else "none"))
    note("what this does NOT establish: DL Table 4 is never read by this "
         "program, and the order test can only ever REFUTE an atomic weight, so "
         "the absence of a refutation is consistent both with the tabulated "
         "weights agreeing with the formula and with their disagreeing in a way "
         "this probe cannot see.  The probe's ability to fire at all is what "
         "the sensitivity self-test above establishes.")
    check(d0["betahat"] == PAPER_BETAHAT and d0["R"] == 12
          and 2 * tau == b - a,
          "the Remark's boundary instance: derived betahat=%d and "
          "R_betahat=%d equal the paper's 2 and 12, and tau = (b-a)/2 exactly"
          % (d0["betahat"], d0["R"]))
    check((not refuted_at) and len(examined) == len(list(heaps)),
          "the Remark's boundary instance is examined by computation instead of "
          "by reading DL Table 4: at (3,7,2) the order test refutes the "
          "conjecture's prediction at NO n in 12..34, which is consistent with "
          "the Remark that tau = (b-a)/2 is not by itself the obstruction.  "
          "This is NOT a confirmation of the Remark: the order test can refute "
          "but never confirm an atomic weight, %d of the %d heap sizes were "
          "swept at j = 1 only, and DL Table 4's tabulated values are never "
          "read here -- a refutation, however, WOULD have contradicted the "
          "paper" % (len(cut), len(examined)),
          "%d heap sizes examined, %d refuted, %d swept at j=1 only"
          % (len(examined), len(refuted_at), len(cut))
          if not refuted_at
          else "REFUTED at n = %s: this CONTRADICTS the paper's Remark and "
               "its abstract" % (refuted_at,))
    return not refuted_at


# ===========================================================================
# Survey (printed derivation, NOT a check): the whole initial segment of each
# family, so a reader can see where R_betahat sits and what happens on either
# side of it.  Nothing here is asserted by the paper.
# ===========================================================================

def survey(tab):
    print("")
    print("--- Survey (no checks): derived values of TS_t(n;1,2t+1) ---")
    refuted_n, inconclusive_n = 0, 0
    for t in (1, 2, 3):
        a, b, tau = 1, 2 * t + 1, t
        R = conjecture_data(a, b, tau, 0)["R"]
        note("t=%d (a=1, b=%d, tau=%d, R_betahat=%d):" % (t, b, tau, R))
        for n in range(0, 2 * t + 9):
            g = ts(n, tau, a, b)
            line = "  n=%2d  value %-22s outcome %s" % (
                n, name_of(g, tab), outcome_class(g))
            if n >= R:
                pred = -((n - R) // a)
                rf, _det, _js = refute_aw_equals(g, pred)
                line += "   prediction aw=%d : %s" % (
                    pred, "REFUTED" if rf else "not refuted by order test")
                if rf:
                    refuted_n += 1
                else:
                    inconclusive_n += 1
            note(line)
    note("surveyed heaps with n >= R_betahat: %d refuted, %d inconclusive "
         "(the order test can refute but never confirm an atomic weight)"
         % (refuted_n, inconclusive_n))


# ===========================================================================
# The gaps between what is checked and what the paper claims.  Printed, not
# hidden: each item is a step no check above covers.
# ===========================================================================

def print_gaps():
    print("")
    print("--- GAPS between the checks above and the paper's claims ---")
    note("G1  SCOPE, and the one that matters most.  Every check is "
         "conditioned on P3, this program's transcription of Conjecture 7.1 "
         "and of its standing hypotheses -- above all that a = 1 is "
         "admissible.  Dargad-Larsson's preprint also carries a '3 <= a < b' "
         "hypothesis elsewhere, in an unrelated section on domination.  If the "
         "conjecture's real scope excludes a = 1, or if its hypotheses are not "
         "0<a<=b, a<b, tau<b-a, n>=R_betahat, then EVERY CHECK HERE STILL "
         "PASSES and the paper is nevertheless wrong.  No offline computation "
         "can close this; it needs a fidelity audit against the preprint "
         "itself, which this program does not attempt because it would need a "
         "network fetch of a third-party document.")
    note("G1a Two concrete ways G1 could bite, both invisible to every check "
         "here: (i) if DL's equations (5)/(13) give an R_betahat larger than "
         "2t+4, then n = 2t+4 is below the smallest governed heap and the "
         "instances are out of scope; (ii) if DL's Definition 5.1 gives Right "
         "the subtraction set {tau,...,b} or {tau+1,...,b-1} rather than "
         "{tau+1,...,b}, the games computed here are the wrong games.  The "
         "rule anchors below test the built trees against the prose of THIS "
         "paper's own proof, which establishes internal consistency, not "
         "fidelity to Dargad-Larsson's Definition 5.1.")
    note("G2  FINITE vs UNIVERSAL.  Theorem 1 is stated for every integer "
         "t >= 1.  It is verified here for t = 1..12 only, by direct "
         "computation; no induction on t is verified.  What IS fully settled "
         "is the Corollary, because falsifying a conjecture needs only ONE "
         "instance: any single t above suffices, so the refutation of "
         "Conjecture 7.1 does not depend on the universal quantifier.")
    note("G3  ATOMIC WEIGHT IS NEVER COMPUTED.  This program implements the "
         "order relation, not the atomic-weight calculus, so it can REFUTE a "
         "predicted atomic weight but can never independently CONFIRM one.  "
         "The imported facts are aw(up)=1, aw(*)=0, additivity, and the "
         "characterisation P6(C): for an all-small G, aw(G)=0 exactly when G is "
         "infinitesimal with respect to up, i.e. when for EVERY integer j >= 1 "
         "the game j*G is less than or confused with up and greater than or "
         "confused with down.  The refutation used here is the contrapositive: "
         "if some j >= 1 gives up <= j*(G - P*up) or j*(G - P*up) <= down, then "
         "aw(G) != P.  Note that the STRICT form 'aw(G)=0 => down < G < up' is "
         "NOT the criterion and is false -- aw(*)=0 while * is confused with "
         "up and with down, as the engine self-test above records -- so it is "
         "not used or relied on anywhere in this program.  aw(down) = -1 and "
         "aw(down+down+*) = -2 are then derived, each from an identity the "
         "order engine decides here.")
    note("G4  BOTH ENGINES SHARE ts().  The order test and the outcome "
         "induction are independent decision procedures, but they read the "
         "same game trees, so their agreement cannot detect a mis-transcribed "
         "subtraction set.  The rule-anchor checks above address that "
         "separately, by testing the built trees against move-by-move facts "
         "the paper's own proof states in prose.")
    note("G5  THE REMARK IS THE WEAKEST CHECK HERE.  The Remark's appeal to DL "
         "Table 4 is answered by a direct computation of TS_2(n;3,7) for "
         "12 <= n <= 34, but that computation can only fail to refute: it does "
         "not read Table 4, it does not confirm the tabulated weights, and at "
         "the largest heap sizes the j-sweep was cut to j = 1 by the node cap "
         "(the check's own detail line prints how many).  A refutation there "
         "would have contradicted the paper's Remark and its abstract, so the "
         "check is not empty -- but its PASS establishes consistency with the "
         "Remark, not the Remark.  The abstract's clause that tau = (b-a)/2 is "
         "not itself the obstruction therefore still rests on Dargad-Larsson's "
         "Table 4, which nothing in this program has read.")


# ===========================================================================
# main
# ===========================================================================

def main():
    print("verify.py -- counterexamples to Dargad-Larsson Conjecture 7.1")
    print("Exact integer / game-tree arithmetic only; no floating point.")
    print("")
    print("--- Engine self-test ---")
    engine_selftest()
    tab = reference_table()

    check_probe_sensitivity()
    check_definition_anchors()
    check_theorem1(tab)
    check_corollary_arithmetic()
    check_refutation(tab)
    check_extended_range()
    check_strengthening_next_heap(tab)
    check_strictly_inside(tab)
    check_remark_boundary(tab)
    survey(tab)
    print_gaps()

    print("")
    print("--- NOT COVERED (declared, not checked) ---")
    note("Fidelity to arXiv:2607.27989v1 itself -- that Conjecture 7.1, its "
         "standing hypotheses, Definition 5.1 and Table 4 are as transcribed in "
         "P1/P3 above -- is NOT checked: it would need a network fetch of a "
         "third-party preprint, and this program is offline and stdlib-only.  "
         "It is NOT counted as a check.")
    note("The two atomic-weight facts P6(A)/(B) and the characterisation "
         "P6(C) are imported from the cited literature, not proved here; "
         "every check that uses them says so on its own line.")
    note("The values of DL Table 4 are never read.  The paper's Remark, and "
         "with it the abstract's clause that tau = (b-a)/2 is not itself the "
         "obstruction, is therefore only tested to the extent that the order "
         "test found no refutation at the boundary instance (3,7,2) for "
         "12 <= n <= 34 -- and at the largest of those heap sizes the sweep was "
         "cut to j = 1 by the node cap, as the check's detail line records.")
    note("Theorem 1 is verified for t = 1..12 only, by direct computation over "
         "finite game trees; no induction on t is verified anywhere, and no "
         "second implementation of these values is shipped alongside this one.")
    note("Total distinct game positions constructed: %d ; order-test memo "
         "entries: %d" % (len(_ALL_GAMES), len(_LEQ)))

    note("Checks whose label says \"strengthening (beyond the paper's claim)\" "
         "test additional instances, not statements of the paper; the paper's "
         "own claims are the Theorem 1 and refutation checks.")

    print("")
    total = len(RESULTS)
    failed = [lab for ok, lab in RESULTS if not ok]
    for lab in failed:
        note("recap of failure: " + lab)
    if failed:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


if __name__ == "__main__":
    sys.exit(main())

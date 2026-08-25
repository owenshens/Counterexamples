#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- verification program for

    "A Counterexample to the Unrestricted Dargad--Larsson Domination Problem"

The paper studies the knife-edge truncated-support partizan subtraction game

    G_n = TS_{b-a}(n; a, b),      S_L = {1,...,a},   S_R = {q,...,b},   q = b-a+1,

    G_n = { G_{n-i} : 1 <= i <= min(a,n)  |  G_{n-j} : q <= j <= min(b,n) },

proves the order  G_r <= G_s  <=>  floor(r/q) = floor(s/q) and r <= s  (0 <= r,s <= b),
and concludes that the unrestricted wording of Dargad--Larsson Problem 7.1 is false
unless q | (a-1); the exhibited counterexample is (a,b) = (3,5), where the Left option
G_4 of G_6 is dominated by G_5, a comparison at index b-1 = 4 lying outside the printed
range {b-a+1,...,b-2} = {3} of their Proposition 7.2.

--------------------------------------------------------------------------------
TAKEN FROM THE PAPER (data / definitions; not re-derived here)
--------------------------------------------------------------------------------
 T1. The ruleset: S_L = {1,...,a}, S_R = {q,...,b} with q = b-a+1, 0 < a < b, and the
     recursion (eq:recursion) for G_n.  This is the object under study.
 T2. The standard short-game comparison criterion (Siegel): X <= Y iff no Left option
     X^L has Y <= X^L and no Right option Y^R has Y^R <= X.  A definition, used as the
     *only* engine; nothing about specific games is assumed.
 T3. The witness parameters: (a,b) = (3,5), heap size 6, and the claim that its Left
     options are the heaps 5, 4, 3.  (The heap sizes are read off the paper; the option
     set itself is recomputed from T1 and compared.)
 T4. Quoted from Dargad--Larsson via the paper: Proposition 7.2's printed comparison
     range  {b-a+1,...,b-3,b-2} \ Q  with  Q = {(m+1)q - 1 : m >= 1};
     Observation 5.2's assertion  TS_tau(n;a,b) = n  for 1 <= n <= tau  (tau = b-a);
     the periodicity  TS_kappa(n) = TS_kappa(n mod (b+1))  with claimed period p = b+1;
     and that Problem 7.1 as printed carries no restriction on the heap size n.
     These three transcriptions are INPUT DATA: nothing here can check them against the
     source, and gap (1) at the end of the run says so.
 T5. The statements to be tested: Theorem 1's closed form for <=, its block
     decomposition B_m = {mq,...,min((m+1)q-1, b)}, Corollary 2's domination rule,
     and the boundary dichotomy "G_{b-1}, G_b incomparable  <=>  q | (a-1)".
 T6. The paper's own convention (its section 3, first paragraph): "given by
     Proposition 7.2" means "derivable from Observation 5.2 and Proposition 7.2 by
     transitivity and the periodicity", and BOTH of the paper's verdicts are taken
     under that reading.  CHECK 7 tests the verdicts under exactly this convention;
     without it, "b-1 lies outside the printed range" would be a strictly weaker
     statement than the one the paper makes.

--------------------------------------------------------------------------------
DERIVED HERE (everything the checks actually decide)
--------------------------------------------------------------------------------
 D1. Every game comparison, from T1+T2 alone, by the memoised predicate leq(n,m)
     ("G_n <= G_m").  No closed form is ever used to produce a comparison.
 D2. The witness: G_3 < G_4 < G_5 for (a,b) = (3,5); the actual Left-option target set
     of heap 6; G_4 <= G_5 with G_5 not <= G_4, hence strict domination; the actual
     Proposition-7.2 comparison index set for (3,5) and the fact that b-1 = 4 misses it.
 D3. The truth value of Theorem 1's iff, tested pointwise on all (r,s) in [0,b]^2 for
     every 0 < a < b <= 20, and the pairwise distinctness of G_0,...,G_b.
 D4. Periodicity G_n = G_{n mod (b+1)} for n <= 4(b+1), and the *least* positive period,
     with an explicit refuting witness n for each smaller candidate period.
 D5. The boundary dichotomy at (b-1, b) for every 0 < a < b < 20, and the equivalence
     of "q | (a-1)" with "b is the first element of its block".
 D6. a = 1: the per-player option counts at every heap size (hence: nothing is dominated).
     a = 2: that Proposition 7.2's range is empty, that B_0 = {0,...,b-2}, that
     Observation 5.2's value G_n = n holds there (checked as a *game* equality against
     the integer game n built from {n-1 | }), and that the residue-b witness survives.
 D7. Corollary 2's domination rule, tested against the dominated sets actually computed
     from D1 at every heap size in a range, via residues mod p = b+1.
 D8. That distinct options of one player always land on distinct residues mod b+1.
 D9. The transitive closure of the comparisons the source supplies (Observation 5.2's
     whole first-block chain plus Proposition 7.2's printed range minus Q), and from
     it BOTH halves of the paper's if-and-only-if:
       - if q does not divide a-1, the needed comparison G_{b-1} <= G_b is not in the
         closure, and no supplied pair touches residue b at all;
       - if q divides a-1, every dominated option at every heap size n <= 3(b+1) has
         a witnessing comparison inside the closure;
       - for every (a,b), every dominated option of every representative position
         0..b is inside the closure (the affirmative contextual verdict);
       - for every a >= 2, some domination is witnessed only by a comparison outside
         Proposition 7.2's supply, i.e. by Observation 5.2.
D10. The premise of the mod-(b+1) reduction, re-derived from the engine at each
     target actually reduced, because CHECK 2b's periodicity sweep stops at b = 10
     while the domination sweeps run further.

HONEST ACCOUNTING.  CHECK 4b (Corollary 2) is a consequence of CHECK 2a, not
independent evidence; and six checks state definitional identities that cannot fail.
main() prints both facts at the end rather than leaving them implicit.

All arithmetic is integer (// and %).  No floating point occurs anywhere, so no check
is numeric or tolerance-based: every check is an exact combinatorial identity.

Python 3.9, standard library only.  Usage: python3 verify.py
"""

import sys

# ----------------------------------------------------------------------------
# Data taken from the paper (see T3, T4, T5 above).
# ----------------------------------------------------------------------------

# T3: the exhibited counterexample.
W_A, W_B = 3, 5              # (a,b) of the counterexample
W_HEAP = 6                   # the heap size at which the domination occurs (= b+1)
W_LEFT_TARGETS = (5, 4, 3)   # paper: "Left may move to G_5, G_4, or G_3"
W_CHAIN = (3, 4, 5)          # paper: "G_3 < G_4 < G_5"
W_DOMINATED, W_DOMINATOR = 4, 5   # paper: "G_4 is dominated by G_5"
W_PROP72_INDEX_SET = (3,)    # paper: Prop 7.2 supplies only G_3 <= G_4 here
W_COMPARISON_INDEX = 4       # paper: the needed comparison sits at index b-1 = 4

# T5 / spec ranges for the exhaustive sweeps.
THM1_BMAX = 20               # Theorem 1 order test:      0 < a < b <= 20
PERIOD_BMAX = 10             # periodicity + least period: 0 < a < b <= 10
PERIOD_REACH = 4             # n up to PERIOD_REACH * (b+1)
BOUNDARY_BMAX = 19           # boundary dichotomy:        0 < a < b < 20
OBS52_BMAX = 14              # Observation 5.2 game-equality test
DOM_BMAX = 12                # Corollary 2 domination-rule test
DOM_REACH = 3                # n up to DOM_REACH * (b+1)
DERIV_BMAX = 16              # CHECK 7: derivability closure, end-to-end iff
DERIV_REACH = 3              # n up to DERIV_REACH * (b+1) in CHECK 7

sys.setrecursionlimit(1000000)

_RESULTS = []


def record(ok, label, detail=""):
    """Print one PASS/FAIL line and remember the outcome."""
    _RESULTS.append((bool(ok), label))
    line = ("PASS " if ok else "FAIL ") + label
    if detail:
        line += "  |  " + detail
    print(line)
    return bool(ok)


def say(text):
    """Print a derived intermediate quantity (indented, never a verdict line)."""
    print("    " + text)


# ----------------------------------------------------------------------------
# The engine: short-game comparison built from the recursion (T1) and the
# comparison criterion (T2) only.
#
# Nodes are ('G', n)  = the heap-n position of TS_{b-a}(.;a,b), and
#           ('I', k)  = the integer game k, defined by 0 = { | } and k = {k-1 | },
# the latter used only to certify Observation 5.2's claim G_n = n as a genuine
# equality of games.  Every option strictly decreases the node's index, so the
# measure index(x) + index(y) strictly decreases in every recursive call and the
# recursion terminates.
# ----------------------------------------------------------------------------

class TSGame(object):
    def __init__(self, a, b):
        if not (0 < a < b):
            raise ValueError("need 0 < a < b, got a=%d b=%d" % (a, b))
        self.a = a
        self.b = b
        self.q = b - a + 1          # q = b - a + 1 = kappa + 1
        self.tau = b - a            # knife-edge truncation kappa = b - a
        self.p = b + 1              # claimed period
        self._memo = {}

    # --- the move rules (T1) -------------------------------------------------
    def left_targets(self, n):
        """Heap sizes Left can move to from heap n: n-i, 1 <= i <= min(a,n)."""
        return [n - i for i in range(1, min(self.a, n) + 1)]

    def right_targets(self, n):
        """Heap sizes Right can move to from heap n: n-j, q <= j <= min(b,n)."""
        return [n - j for j in range(self.q, min(self.b, n) + 1)]

    def _options(self, x):
        kind, n = x
        if kind == 'G':
            return ([('G', t) for t in self.left_targets(n)],
                    [('G', t) for t in self.right_targets(n)])
        if n <= 0:
            return ([], [])                      # integer 0 = { | }
        return ([('I', n - 1)], [])               # integer k = {k-1 | }

    # --- the comparison criterion (T2) --------------------------------------
    def leq_node(self, x, y):
        key = (x, y)
        hit = self._memo.get(key)
        if hit is not None:
            return hit
        xl = self._options(x)[0]
        yr = self._options(y)[1]
        res = True
        for xo in xl:
            if self.leq_node(y, xo):
                res = False
                break
        if res:
            for yo in yr:
                if self.leq_node(yo, x):
                    res = False
                    break
        self._memo[key] = res
        return res

    # --- convenience wrappers on heap indices -------------------------------
    def leq(self, n, m):
        return self.leq_node(('G', n), ('G', m))

    def eq(self, n, m):
        return self.leq(n, m) and self.leq(m, n)

    def lt(self, n, m):
        return self.leq(n, m) and not self.leq(m, n)

    def incomparable(self, n, m):
        return (not self.leq(n, m)) and (not self.leq(m, n))

    def eq_integer(self, n, k):
        """True iff the heap-n game equals the integer game k."""
        return (self.leq_node(('G', n), ('I', k)) and
                self.leq_node(('I', k), ('G', n)))


_GAMES = {}


def game(a, b):
    key = (a, b)
    g = _GAMES.get(key)
    if g is None:
        g = TSGame(a, b)
        _GAMES[key] = g
    return g


# ----------------------------------------------------------------------------
# The paper's / the source's combinatorial predicates, transcribed literally so
# that the derived game facts can be compared against them.  These functions
# make *claims*; the TSGame engine above makes the *findings*.
# ----------------------------------------------------------------------------

def paper_predicts_leq(r, s, q):
    """Theorem 1: G_r <= G_s  <=>  floor(r/q) = floor(s/q) and r <= s."""
    return (r // q == s // q) and (r <= s)


def block_index(n, q):
    """The q-block containing n, i.e. m with n in B_m."""
    return n // q


def block(m, q, b):
    """B_m = {mq, ..., min((m+1)q - 1, b)}  (Theorem 1); [] if it misses [0,b]."""
    lo = m * q
    hi = min((m + 1) * q - 1, b)
    if lo > b:
        return []
    return list(range(lo, hi + 1))


def prop72_range(a, b):
    """Prop 7.2's printed comparison indices before removing Q:
    {b-a+1, ..., b-3, b-2} = {q, ..., b-2}, empty when q > b-2."""
    q = b - a + 1
    return list(range(q, b - 1))


def in_Q(n, q):
    """Q = {(m+1)q - 1 : m >= 1} = {2q-1, 3q-1, ...}; membership test, no float."""
    if (n + 1) % q != 0:
        return False
    m = (n + 1) // q - 1          # n = (m+1)q - 1
    return m >= 1


def Q_upto(limit, q):
    """The elements of Q that are <= limit, generated from the definition."""
    out = []
    m = 1
    while (m + 1) * q - 1 <= limit:
        out.append((m + 1) * q - 1)
        m += 1
    return out


def prop72_supplied_indices(a, b):
    """The comparison indices Proposition 7.2 actually supplies: printed range
    minus the exceptional set Q."""
    q = b - a + 1
    return [n for n in prop72_range(a, b) if not in_Q(n, q)]


def divides(d, n):
    """Exact integer divisibility d | n (d > 0)."""
    return n % d == 0


# ----------------------------------------------------------------------------
# CHECK 1 -- THE WITNESS.  (a,b) = (3,5), q = 3.
# ----------------------------------------------------------------------------

def check_witness():
    a, b = W_A, W_B
    g = game(a, b)
    print("[CHECK 1] the witness (a,b) = (%d,%d)" % (a, b))
    say("derived q = b-a+1 = %d ; kappa = b-a = %d ; claimed period p = b+1 = %d"
        % (g.q, g.tau, g.p))
    say("derived S_L = %s ; S_R = %s"
        % (sorted(range(1, a + 1)), sorted(range(g.q, b + 1))))
    for n in range(0, b + 2):
        say("heap %d: Left -> %s   Right -> %s"
            % (n, g.left_targets(n), g.right_targets(n)))
    rel = []
    for r in range(0, b + 1):
        rel.append("".join("<" if g.lt(r, s) else ("=" if g.eq(r, s) else
                   (">" if g.lt(s, r) else "|")) for s in range(0, b + 1)))
    for r in range(0, b + 1):
        say("order row r=%d vs s=0..%d : %s" % (r, b, rel[r]))

    r0, r1, r2 = W_CHAIN
    chain = g.lt(r0, r1) and g.lt(r1, r2)
    record(chain, "witness chain G_%d < G_%d < G_%d for (a,b)=(%d,%d)"
           % (r0, r1, r2, a, b),
           "G_%d<G_%d=%s, G_%d<G_%d=%s" % (r0, r1, g.lt(r0, r1),
                                           r1, r2, g.lt(r1, r2)))

    lt6 = g.left_targets(W_HEAP)
    record(sorted(lt6, reverse=True) == sorted(W_LEFT_TARGETS, reverse=True),
           "heap %d has Left options exactly to heaps %s"
           % (W_HEAP, ",".join(str(t) for t in W_LEFT_TARGETS)),
           "derived Left targets of heap %d = %s" % (W_HEAP, lt6))

    dom = W_DOMINATED
    by = W_DOMINATOR
    both = (dom in lt6) and (by in lt6)
    strict = g.leq(dom, by) and not g.leq(by, dom)
    record(both and strict,
           "Left option to heap %d is strictly dominated by the option to heap %d"
           % (dom, by),
           "both are Left options of heap %d: %s ; G_%d<=G_%d=%s ; G_%d<=G_%d=%s"
           % (W_HEAP, both, dom, by, g.leq(dom, by), by, dom, g.leq(by, dom)))

    printed = prop72_range(a, b)
    supplied = prop72_supplied_indices(a, b)
    qset = Q_upto(b, g.q)
    say("derived Prop 7.2 printed index range {q,...,b-2} = %s ; Q cap [0,%d] = %s ; "
        "supplied indices = %s" % (printed, b, qset, supplied))
    say("the comparison the counterexample needs is G_%d <= G_%d, index b-1 = %d"
        % (dom, by, W_COMPARISON_INDEX))
    ok = (printed == list(W_PROP72_INDEX_SET) and
          supplied == list(W_PROP72_INDEX_SET) and
          W_COMPARISON_INDEX not in printed and
          W_COMPARISON_INDEX == b - 1)
    record(ok,
           "needed index b-1=%d lies outside Prop 7.2's range, which is exactly {%s}"
           % (W_COMPARISON_INDEX, ",".join(str(x) for x in W_PROP72_INDEX_SET)),
           "printed=%s supplied=%s needed=%d" % (printed, supplied,
                                                 W_COMPARISON_INDEX))


# ----------------------------------------------------------------------------
# CHECK 2a -- THEOREM 1: the order on the b+1 periodic values.
# ----------------------------------------------------------------------------

def check_theorem1_order():
    print("[CHECK 2a] Theorem 1's closed form on [0,b]^2 for all 0 < a < b <= %d"
          % THM1_BMAX)
    pairs = 0
    params = 0
    true_rel = 0
    viol = []
    dist_viol = []
    for b in range(2, THM1_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            params += 1
            for r in range(0, b + 1):
                for s in range(0, b + 1):
                    pairs += 1
                    got = g.leq(r, s)
                    want = paper_predicts_leq(r, s, g.q)
                    if got:
                        true_rel += 1
                    if got != want:
                        if len(viol) < 6:
                            viol.append((a, b, r, s, got, want))
                    if r != s and g.eq(r, s):
                        if len(dist_viol) < 6:
                            dist_viol.append((a, b, r, s))
    say("(a,b) pairs tested = %d ; (r,s) comparisons decided = %d ; "
        "of these G_r<=G_s held %d times" % (params, pairs, true_rel))
    for v in viol:
        say("VIOLATION a=%d b=%d r=%d s=%d derived_leq=%s paper_predicts=%s" % v)
    record(not viol,
           "G_r <= G_s  <=>  floor(r/q)=floor(s/q) and r<=s, for all (r,s) in [0,b]^2, "
           "all 0<a<b<=%d" % THM1_BMAX,
           "%d/%d comparisons agree" % (pairs - len(viol), pairs))
    for v in dist_viol:
        say("EQUAL PAIR a=%d b=%d r=%d s=%d (should be distinct)" % v)
    record(not dist_viol,
           "G_0,...,G_b are pairwise distinct for all 0<a<b<=%d" % THM1_BMAX,
           "%d parameter pairs, no two of the b+1 values equal" % params)


# ----------------------------------------------------------------------------
# CHECK 2b -- PERIODICITY G_n = G_{n mod (b+1)}, and minimality of the period.
# ----------------------------------------------------------------------------

def check_periodicity():
    print("[CHECK 2b] periodicity to n = %d(b+1) and least period, 0 < a < b <= %d"
          % (PERIOD_REACH, PERIOD_BMAX))
    tested = 0
    viol = []
    minfail = []
    shown = 0
    for b in range(2, PERIOD_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            nmax = PERIOD_REACH * g.p
            for n in range(0, nmax + 1):
                tested += 1
                if not g.eq(n, n % g.p):
                    if len(viol) < 6:
                        viol.append((a, b, n, n % g.p))
            # minimality: every 0 < d < p must fail to be a period, with a witness
            wits = []
            for d in range(1, g.p):
                w = None
                for n in range(0, g.p):
                    if n + d <= nmax and not g.eq(n, n + d):
                        w = n
                        break
                if w is None:
                    minfail.append((a, b, d))
                else:
                    wits.append((d, w))
            if shown < 3:
                say("a=%d b=%d p=b+1=%d : smaller candidate periods refuted by "
                    "witnesses (d, n with G_n != G_{n+d}) = %s"
                    % (a, b, g.p, wits))
                shown += 1
    say("heap sizes n checked against their residue: %d" % tested)
    for v in viol:
        say("VIOLATION a=%d b=%d : G_%d != G_%d" % v)
    record(not viol,
           "G_n = G_{n mod (b+1)} for all n <= %d(b+1), all 0<a<b<=%d"
           % (PERIOD_REACH, PERIOD_BMAX),
           "%d/%d heap sizes agree with their residue" % (tested - len(viol), tested))
    for v in minfail:
        say("NO WITNESS a=%d b=%d : d=%d was not refuted as a period" % v)
    record(not minfail,
           "no 0 < d < b+1 is a period, so b+1 is the least positive period, "
           "all 0<a<b<=%d" % PERIOD_BMAX,
           "every smaller candidate period has an explicit refuting heap size")


# ----------------------------------------------------------------------------
# CHECK 3a -- THE BOUNDARY: G_{b-1} vs G_b, and the role of q | (a-1).
# ----------------------------------------------------------------------------

def check_boundary():
    print("[CHECK 3a] boundary pair (b-1, b) for all 0 < a < b < %d"
          % (BOUNDARY_BMAX + 1))
    inc_bad = []
    lt_bad = []
    blk_bad = []
    domi_bad = []
    n_div = 0
    n_ndiv = 0
    shown = 0
    for b in range(2, BOUNDARY_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            div = divides(g.q, a - 1)
            inc = g.incomparable(b - 1, b)
            lt = g.lt(b - 1, b)
            if div:
                n_div += 1
            else:
                n_ndiv += 1
            if inc != div:
                if len(inc_bad) < 6:
                    inc_bad.append((a, b, g.q, div, inc))
            if lt != (not div):
                if len(lt_bad) < 6:
                    lt_bad.append((a, b, g.q, div, lt))
            # "the final block is the singleton {b} exactly when q | (a-1)"
            fin = block(block_index(b, g.q), g.q, b)
            if (fin == [b]) != div:
                if len(blk_bad) < 6:
                    blk_bad.append((a, b, g.q, div, fin))
            # when q does not divide a-1: a>=2 and heaps b-1, b are both Left
            # options of heap b+1, so G_{b-1} really is a dominated option there
            if not div:
                lts = g.left_targets(b + 1)
                if not (a >= 2 and (b in lts) and (b - 1 in lts) and
                        g.lt(b - 1, b)):
                    if len(domi_bad) < 6:
                        domi_bad.append((a, b, g.q, lts))
            if shown < 6:
                say("a=%d b=%d q=%d : q|(a-1)=%s -> derived incomparable=%s, "
                    "G_{b-1}<G_b=%s, block of b = %s"
                    % (a, b, g.q, div, inc, lt, fin))
                shown += 1
    say("parameter pairs with q|(a-1): %d ; with q not dividing a-1: %d"
        % (n_div, n_ndiv))
    for v in inc_bad:
        say("VIOLATION a=%d b=%d q=%d div=%s incomparable=%s" % v)
    record(not inc_bad,
           "G_{b-1} and G_b are incomparable exactly when q | (a-1), all 0<a<b<%d"
           % (BOUNDARY_BMAX + 1),
           "%d parameter pairs" % (n_div + n_ndiv))
    for v in lt_bad:
        say("VIOLATION a=%d b=%d q=%d div=%s lt=%s" % v)
    record(not lt_bad,
           "G_{b-1} < G_b exactly when q does not divide a-1, all 0<a<b<%d"
           % (BOUNDARY_BMAX + 1),
           "%d parameter pairs" % (n_div + n_ndiv))
    for v in blk_bad:
        say("VIOLATION a=%d b=%d q=%d div=%s final block=%s" % v)
    record(not blk_bad,
           "the block of b is the singleton {b} exactly when q | (a-1)",
           "%d parameter pairs" % (n_div + n_ndiv))
    for v in domi_bad:
        say("VIOLATION a=%d b=%d q=%d Left targets of b+1 = %s" % v)
    record(not domi_bad,
           "when q does not divide a-1: a>=2 and heaps b-1,b are both Left options "
           "of heap b+1, so G_{b-1} is a dominated option at index b-1",
           "%d such parameter pairs, all exhibiting the residue-b domination" % n_ndiv)


# ----------------------------------------------------------------------------
# Domination, computed from the derived order alone (no closed form used).
# A Left option u is dominated if some other Left option v has G_u <= G_v;
# a Right option u is dominated if some other Right option v has G_v <= G_u.
# Targets are reduced mod p = b+1 first, which is legitimate only because the
# periodicity is itself checked in CHECK 2b.
# ----------------------------------------------------------------------------

def option_residues(g, n, side):
    ts = g.left_targets(n) if side == 'L' else g.right_targets(n)
    return [t % g.p for t in ts]


def reduction_violations(g, n):
    """The reduction of option targets mod p = b+1 is only licensed by the
    periodicity.  CHECK 2b proves periodicity for b <= PERIOD_BMAX only, while
    the domination sweeps run to larger b, so the premise is re-established
    HERE, from the engine, for every target actually reduced.  Returns the list
    of targets t > b for which G_t = G_{t mod p} could not be derived."""
    bad = []
    for side in ('L', 'R'):
        ts = g.left_targets(n) if side == 'L' else g.right_targets(n)
        for t in ts:
            if t > g.b and not g.eq(t, t % g.p):
                bad.append((n, side, t, t % g.p))
    return bad


def dominated_residues(g, n, side):
    """The residues of the dominated options of that player at heap n."""
    res = option_residues(g, n, side)
    out = []
    for u in res:
        for v in res:
            if u == v:
                continue
            if (g.leq(u, v) if side == 'L' else g.leq(v, u)):
                out.append(u)
                break
    return sorted(set(out))


def domination_witnesses(g, n, side):
    """[(u, (x,y))] : u is a dominated option residue at heap n for that player,
    and "G_x <= G_y" is a comparison that witnesses it (x=u,y=v for Left, since
    a Left option is dominated when G_u <= G_v; x=v,y=u for Right).  All
    witnesses are listed, so a later "is it derivable" test may pick any."""
    res = sorted(set(option_residues(g, n, side)))
    out = []
    for u in res:
        for v in res:
            if u == v:
                continue
            if side == 'L':
                if g.leq(u, v):
                    out.append((u, (u, v)))
            else:
                if g.leq(v, u):
                    out.append((u, (v, u)))
    return out


def corollary_predicts_dominated(g, n, side):
    """Corollary 2: within each q-block of residues, all options but the largest
    (Left) / smallest (Right) residue are dominated; nothing across blocks."""
    res = sorted(set(option_residues(g, n, side)))
    keep = {}
    for u in res:
        m = block_index(u, g.q)
        if m not in keep:
            keep[m] = u
        elif side == 'L':
            keep[m] = max(keep[m], u)
        else:
            keep[m] = min(keep[m], u)
    survivors = set(keep.values())
    return sorted(u for u in res if u not in survivors)


# ----------------------------------------------------------------------------
# CHECK 3b -- a = 1: one legal move size per player, so nothing is dominated.
# ----------------------------------------------------------------------------

def check_a1():
    print("[CHECK 3b] a = 1 : single option per player, no domination anywhere")
    a = 1
    setsize_bad = []
    count_bad = []
    dom_bad = []
    div_bad = []
    heaps = 0
    for b in range(2, BOUNDARY_BMAX + 1):
        g = game(a, b)
        sl = list(range(1, a + 1))
        sr = list(range(g.q, b + 1))
        if not (len(sl) == 1 and len(sr) == 1):
            setsize_bad.append((b, sl, sr))
        if not divides(g.q, a - 1):
            div_bad.append((b, g.q))
        for n in range(0, DOM_REACH * g.p + 1):
            heaps += 1
            nl = len(g.left_targets(n))
            nr = len(g.right_targets(n))
            if nl > 1 or nr > 1:
                if len(count_bad) < 6:
                    count_bad.append((b, n, nl, nr))
            dl = dominated_residues(g, n, 'L')
            dr = dominated_residues(g, n, 'R')
            if dl or dr:
                if len(dom_bad) < 6:
                    dom_bad.append((b, n, dl, dr))
        if b <= 4:
            say("a=1 b=%d : q=%d, S_L=%s, S_R=%s, max |Left opts|=%d, "
                "max |Right opts|=%d over n<=%d"
                % (b, g.q, sl, sr,
                   max(len(g.left_targets(n)) for n in range(0, DOM_REACH * g.p + 1)),
                   max(len(g.right_targets(n)) for n in range(0, DOM_REACH * g.p + 1)),
                   DOM_REACH * g.p))
    say("heap sizes examined for a=1: %d (b = 2..%d)" % (heaps, BOUNDARY_BMAX))
    for v in setsize_bad:
        say("VIOLATION b=%d S_L=%s S_R=%s" % v)
    for v in count_bad:
        say("VIOLATION b=%d n=%d |Left|=%d |Right|=%d" % v)
    for v in dom_bad:
        say("VIOLATION b=%d n=%d dominated Left=%s dominated Right=%s" % v)
    record(not setsize_bad and not count_bad and not dom_bad,
           "a=1: |S_L|=|S_R|=1, every heap has at most one option per player, and no "
           "option is dominated at any heap size <= %d(b+1)" % DOM_REACH,
           "%d heap sizes, %d values of b, zero dominated options"
           % (heaps, BOUNDARY_BMAX - 1))
    for v in div_bad:
        say("VIOLATION b=%d q=%d does not divide a-1=0" % v)
    record(not div_bad,
           "a=1: q=b divides a-1=0, so the boundary dichotomy also calls the "
           "unrestricted statement true here",
           "checked for b = 2..%d" % BOUNDARY_BMAX)


# ----------------------------------------------------------------------------
# CHECK 3c -- a = 2: Prop 7.2's range is empty, Observation 5.2 covers exactly
# B_0 = {0,...,b-2}, and the residue-b witness survives.
# ----------------------------------------------------------------------------

def check_a2():
    print("[CHECK 3c] a = 2 : empty Prop 7.2 range, B_0 = {0,...,b-2}, "
          "residue-b witness present")
    a = 2
    rng_bad = []
    b0_bad = []
    obs_bad = []
    wit_bad = []
    n_b = 0
    for b in range(3, BOUNDARY_BMAX + 1):
        g = game(a, b)
        n_b += 1
        rng = prop72_range(a, b)
        if rng:
            rng_bad.append((b, g.q, rng))
        b0 = block(0, g.q, b)
        if b0 != list(range(0, b - 1)):
            b0_bad.append((b, b0, list(range(0, b - 1))))
        # Observation 5.2 on B_0: G_n = n for 0 <= n <= tau = b-a, and Right is stuck
        for n in b0:
            if g.right_targets(n) or not g.eq_integer(n, n):
                if len(obs_bad) < 6:
                    obs_bad.append((b, n, g.right_targets(n), g.eq_integer(n, n)))
        # the residue-b witness
        div = divides(g.q, a - 1)
        lts = g.left_targets(b + 1)
        good = ((not div) and g.lt(b - 1, b) and (b in lts) and (b - 1 in lts) and
                (b - 1) in dominated_residues(g, b + 1, 'L'))
        if not good:
            if len(wit_bad) < 6:
                wit_bad.append((b, g.q, div, g.lt(b - 1, b), lts))
        if b <= 6:
            say("a=2 b=%d : q=%d, tau=%d, Prop 7.2 range %s, B_0=%s, "
                "q|(a-1)=%s, G_%d<G_%d=%s, dominated Left residues at heap %d = %s"
                % (b, g.q, g.tau, rng, b0, div, b - 1, b, g.lt(b - 1, b),
                   b + 1, dominated_residues(g, b + 1, 'L')))
    for v in rng_bad:
        say("VIOLATION b=%d q=%d Prop 7.2 range not empty: %s" % v)
    record(not rng_bad,
           "a=2: Proposition 7.2's range {q,...,b-2} is empty for every b, so it "
           "supplies no comparison",
           "%d values of b (3..%d)" % (n_b, BOUNDARY_BMAX))
    for v in b0_bad:
        say("VIOLATION b=%d derived B_0=%s expected %s" % v)
    for v in obs_bad:
        say("VIOLATION b=%d n=%d Right targets=%s G_n equals integer n: %s" % v)
    record(not b0_bad and not obs_bad,
           "a=2: B_0 = {0,...,b-2}, Right has no move there, and G_n equals the "
           "integer game n on B_0 (Observation 5.2)",
           "%d values of b, all of B_0 certified as an integer game" % n_b)
    for v in wit_bad:
        say("VIOLATION b=%d q=%d div=%s lt=%s Left targets of b+1=%s" % v)
    record(not wit_bad,
           "a=2: q=b-1 never divides a-1=1, G_{b-1}<G_b, and heap b+1 has heap b-1 "
           "as a dominated Left option -- the residue-b witness survives",
           "%d values of b (3..%d)" % (n_b, BOUNDARY_BMAX))


# ----------------------------------------------------------------------------
# CHECK 4 -- Corollary 2 (the complete domination rule) and the claim that
# distinct options of one player have distinct residues mod p = b+1.
# ----------------------------------------------------------------------------

def check_domination_corollary():
    print("[CHECK 4] Corollary 2's domination rule and residue distinctness, "
          "0 < a < b <= %d, n <= %d(b+1)" % (DOM_BMAX, DOM_REACH))
    bad = []
    res_bad = []
    spread_bad = []
    red_bad = []
    positions = 0
    dominated_seen = 0
    for b in range(2, DOM_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            # Move-set spreads, recovered from the game object's OWN move
            # generator at a heap large enough that no truncation applies, so
            # this is not a restatement of range(1,a+1).  NOTE the honest
            # accounting: "spread == a-1" is definitional (S_L={1..a} and
            # S_R={q..b} with q=b-a+1), so only "spread < p" carries content;
            # the substantive claim of this check is res_bad below.
            nbig = b + g.p
            move_L = sorted(nbig - t for t in g.left_targets(nbig))
            move_R = sorted(nbig - t for t in g.right_targets(nbig))
            spread_L = max(move_L) - min(move_L)
            spread_R = max(move_R) - min(move_R)
            if not (spread_L == a - 1 and spread_R == a - 1 and
                    spread_L < g.p and spread_R < g.p):
                spread_bad.append((a, b, spread_L, spread_R, g.p))
            for n in range(0, DOM_REACH * g.p + 1):
                positions += 1
                for v in reduction_violations(g, n):
                    if len(red_bad) < 6:
                        red_bad.append((a, b) + v)
                for side in ('L', 'R'):
                    rs = option_residues(g, n, side)
                    if len(set(rs)) != len(rs):
                        if len(res_bad) < 6:
                            res_bad.append((a, b, n, side, rs))
                    got = dominated_residues(g, n, side)
                    want = corollary_predicts_dominated(g, n, side)
                    dominated_seen += len(got)
                    if got != want:
                        if len(bad) < 6:
                            bad.append((a, b, n, side, got, want))
    say("positions examined = %d ; dominated options found = %d"
        % (positions, dominated_seen))
    for v in spread_bad:
        say("VIOLATION a=%d b=%d spreads %d,%d vs p=%d" % v)
    for v in res_bad:
        say("VIOLATION a=%d b=%d n=%d side=%s repeated residues %s" % v)
    record(not spread_bad and not res_bad,
           "distinct options of one player always have distinct residues mod b+1 "
           "(spreads read back off the move generator; note 'spread = a-1' is "
           "definitional, the content is the absence of residue collisions)",
           "%d positions, no residue collision" % positions)
    for v in red_bad:
        say("VIOLATION a=%d b=%d n=%d side=%s target %d not equal to residue %d" % v)
    record(not red_bad,
           "every option target reduced mod b+1 in this check was independently "
           "shown equal to its residue (the reduction's premise, re-derived here "
           "because CHECK 2b only reaches b <= %d)" % PERIOD_BMAX,
           "%d positions" % positions)
    for v in bad:
        say("VIOLATION a=%d b=%d n=%d side=%s derived dominated=%s corollary=%s" % v)
    record(not bad,
           "the dominated options at every heap size are exactly those Corollary 2 "
           "names (blockwise non-maximal for Left, non-minimal for Right)",
           "%d positions, %d dominated options, all accounted for"
           % (positions, dominated_seen))
    # ANTI-VACUITY.  Both sides of the comparison above are empty at a position
    # with fewer than two options, so an empty universe would print PASS.  The
    # check is only meaningful if dominations were actually found.
    record(positions > 0 and dominated_seen > 0,
           "anti-vacuity: the domination sweep was non-empty and did find dominated "
           "options (so the agreement above is not agreement over an empty set)",
           "positions=%d dominated options found=%d" % (positions, dominated_seen))


# ----------------------------------------------------------------------------
# The paper's OWN convention (tex lines 206-213):
#
#   "Throughout, 'given by Proposition 7.2' means 'derivable from Observation 5.2
#    and Proposition 7.2 by transitivity and the periodicity'.  ... both verdicts
#    below are taken under it."
#
# So neither verdict is about the printed index range alone: the negative verdict
# needs G_{b-1} <= G_b to be UNREACHABLE in that closure, and the positive verdict
# needs every domination to be REACHABLE in it.  The two functions below build the
# supplied relation and close it; CHECK 7 then decides both verdicts.
# ----------------------------------------------------------------------------

def supplied_comparisons(g):
    """Ordered pairs (u,v), u,v in [0,b], for which "G_u <= G_v" is supplied by
    the source.  Deliberately as GENEROUS as the source allows, so that a
    non-derivability verdict below is as strong as possible:

      * reflexivity, (u,u);
      * Observation 5.2: G_n = n for 1 <= n <= tau and G_0 = 0, so the whole
        chain 0 <= 1 <= ... <= tau is supplied, not merely consecutive steps;
      * Proposition 7.2: (n, n+1) for every n in its printed range {q,...,b-2}
        minus the exceptional set Q.

    Periodicity contributes nothing further because every index here is already
    a residue mod p = b+1: that is exactly how periodicity is used."""
    pairs = set()
    for u in range(0, g.b + 1):
        pairs.add((u, u))
    hi = min(g.tau, g.b)
    for u in range(0, hi + 1):
        for v in range(u, hi + 1):
            pairs.add((u, v))
    for n in prop72_supplied_indices(g.a, g.b):
        if 0 <= n and n + 1 <= g.b:
            pairs.add((n, n + 1))
    return pairs


def transitive_closure(pairs, b):
    """Reflexive-transitive closure of a relation on {0,...,b} (Floyd-Warshall,
    pure integer/boolean; no arithmetic on the relation itself)."""
    reach = [[False] * (b + 1) for _ in range(b + 1)]
    for (u, v) in pairs:
        reach[u][v] = True
    for u in range(0, b + 1):
        reach[u][u] = True
    for k in range(0, b + 1):
        rk = reach[k]
        for i in range(0, b + 1):
            if reach[i][k]:
                ri = reach[i]
                for j in range(0, b + 1):
                    if rk[j]:
                        ri[j] = True
    return reach


# ----------------------------------------------------------------------------
# CHECK 5 -- the source's own comparisons: for q <= n <= b-1,
#   G_n < G_{n+1}  <=>  n not in Q,  and the two are incomparable when n in Q;
# and Observation 5.2 + Proposition 7.2 together supply exactly the consecutive
# comparisons that hold among the option values G_0,...,G_{b-1}.
# ----------------------------------------------------------------------------

def check_prop72_consecutive():
    print("[CHECK 5] the Q-dichotomy on consecutive pairs and coverage of "
          "G_0..G_{b-1}, 0 < a < b <= %d" % THM1_BMAX)
    qbad = []
    covbad = []
    tested = 0
    shown = 0
    for b in range(2, THM1_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            for n in range(g.q, b):
                tested += 1
                inq = in_Q(n, g.q)
                if g.lt(n, n + 1) != (not inq) or g.incomparable(n, n + 1) != inq:
                    if len(qbad) < 6:
                        qbad.append((a, b, n, inq, g.lt(n, n + 1),
                                     g.incomparable(n, n + 1)))
            first_block = [n for n in range(0, b - 1) if n <= g.q - 2]
            supplied = set(first_block) | set(prop72_supplied_indices(a, b))
            holds = set(n for n in range(0, b - 1) if g.lt(n, n + 1))
            if holds != supplied:
                if len(covbad) < 6:
                    covbad.append((a, b, sorted(holds), sorted(supplied)))
            if shown < 5 and a >= 3:
                say("a=%d b=%d q=%d : Q cap [0,b] = %s ; consecutive comparisons "
                    "holding on n in [0,b-2] = %s ; Obs 5.2 gives %s, Prop 7.2 gives %s"
                    % (a, b, g.q, Q_upto(b, g.q), sorted(holds), first_block,
                       prop72_supplied_indices(a, b)))
                shown += 1
    say("consecutive pairs (n, n+1) with q <= n <= b-1 tested: %d" % tested)
    for v in qbad:
        say("VIOLATION a=%d b=%d n=%d in_Q=%s lt=%s incomparable=%s" % v)
    record(not qbad,
           "for q <= n <= b-1: G_n < G_{n+1} iff n not in Q, incomparable iff n in Q, "
           "all 0<a<b<=%d" % THM1_BMAX,
           "%d consecutive pairs" % tested)
    for v in covbad:
        say("VIOLATION a=%d b=%d holding=%s supplied=%s" % v)
    record(not covbad,
           "Observation 5.2 (first block) plus Proposition 7.2 (range minus Q) supply "
           "exactly the consecutive comparisons that hold among G_0..G_{b-1}",
           "checked for every 0<a<b<=%d" % THM1_BMAX)


# ----------------------------------------------------------------------------
# CHECK 5c -- the paper prints Q twice, as {(m+1)(b-a)+m : m>=1} and as
# {(m+1)q-1 : m>=1}, and asserts the two agree.  Only the second form is used
# elsewhere in this program, so the first is checked against it here; a
# transcription slip between them would silently move every Q-exception.
# ----------------------------------------------------------------------------

def check_Q_two_forms():
    print("[CHECK 5c] the paper's two printed descriptions of Q agree")
    bad = []
    tested = 0
    for b in range(2, THM1_BMAX + 1):
        for a in range(1, b):
            q = b - a + 1
            tau = b - a
            for m in range(1, 2 * b + 3):
                tested += 1
                form1 = (m + 1) * tau + m        # {(m+1)(b-a)+m}
                form2 = (m + 1) * q - 1          # {(m+1)q-1}
                if form1 != form2:
                    if len(bad) < 6:
                        bad.append((a, b, m, form1, form2))
            # in_Q must accept exactly the generated elements and nothing else.
            # This is falsifiable: it pins the m>=1 floor, so q-1 (which is
            # (m+1)q-1 at m=0) must be REJECTED.  Dropping that floor would add
            # a spurious exception at the end of the first block and change every
            # coverage verdict below.
            gen = set(Q_upto(4 * b, q))
            for n in range(0, 4 * b + 1):
                if in_Q(n, q) != (n in gen):
                    if len(bad) < 6:
                        bad.append((a, b, n, int(in_Q(n, q)), int(n in gen)))
            if in_Q(q - 1, q):
                if len(bad) < 6:
                    bad.append((a, b, q - 1, -1, -1))
    say("(a,b,m) triples tested = %d" % tested)
    for v in bad:
        say("VIOLATION a=%d b=%d m=%d (m+1)(b-a)+m=%d (m+1)q-1=%d" % v)
    record(tested > 0 and not bad,
           "Q's two printed forms {(m+1)(b-a)+m} and {(m+1)q-1} coincide term by term "
           "(a transcription cross-check: it fails only if the paper's two "
           "descriptions disagree), and in_Q accepts exactly the generated set, "
           "rejecting q-1",
           "%d triples" % tested)


# ----------------------------------------------------------------------------
# CHECK 6 -- Observation 5.2 in general: TS_tau(n;a,b) = n for 1 <= n <= tau,
# certified as an equality with the integer game n, and B_0 = {0,...,tau}.
# ----------------------------------------------------------------------------

def check_obs52():
    print("[CHECK 6] Observation 5.2: G_n = n for 1 <= n <= tau = b-a, and "
          "B_0 = {0,...,tau}, 0 < a < b <= %d" % OBS52_BMAX)
    intbad = []
    rightbad = []
    b0bad = []
    tested = 0
    shown = 0
    for b in range(2, OBS52_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            if block(0, g.q, b) != list(range(0, min(g.tau, b) + 1)):
                b0bad.append((a, b, block(0, g.q, b), g.tau))
            for n in range(1, g.tau + 1):
                tested += 1
                if g.right_targets(n):
                    if len(rightbad) < 6:
                        rightbad.append((a, b, n, g.right_targets(n)))
                if not g.eq_integer(n, n):
                    if len(intbad) < 6:
                        intbad.append((a, b, n))
            if shown < 4:
                say("a=%d b=%d : tau=%d, B_0=%s, heaps 1..tau certified equal to the "
                    "integer games 1..%d" % (a, b, g.tau, block(0, g.q, b), g.tau))
                shown += 1
    say("heap sizes n with 1 <= n <= tau certified as integers: %d" % tested)
    for v in rightbad:
        say("VIOLATION a=%d b=%d n=%d Right can move to %s" % v)
    for v in intbad:
        say("VIOLATION a=%d b=%d : G_%d does not equal the integer game %d"
            % (v[0], v[1], v[2], v[2]))
    record(not rightbad and not intbad,
           "for 1 <= n <= tau Right has no move and G_n equals the integer game n, "
           "all 0<a<b<=%d" % OBS52_BMAX,
           "%d heap sizes certified" % tested)
    for v in b0bad:
        say("VIOLATION a=%d b=%d B_0=%s tau=%d" % v)
    record(not b0bad,
           "the first block B_0 is exactly {0,...,tau}, so Observation 5.2 covers it "
           "entirely",
           "all 0<a<b<=%d" % OBS52_BMAX)


# ----------------------------------------------------------------------------
# CHECK 7a -- THE NEGATIVE VERDICT, under the paper's own convention.
# Not "b-1 is outside the printed range" (that is weaker than the paper's claim)
# but "G_{b-1} <= G_b is not in the transitive closure of what the source supplies".
# ----------------------------------------------------------------------------

def check_negative_verdict():
    print("[CHECK 7a] the counterexample under the paper's stated convention: "
          "G_{b-1} <= G_b is NOT derivable from Obs 5.2 + Prop 7.2 by "
          "transitivity and periodicity, 0 < a < b <= %d" % DERIV_BMAX)
    bad = []
    n_case = 0
    shown = 0
    for b in range(2, DERIV_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            if divides(g.q, a - 1):
                continue                      # positive regime, handled in 7b
            n_case += 1
            sup = supplied_comparisons(g)
            reach = transitive_closure(sup, b)
            derivable = reach[b - 1][b]
            holds = g.lt(b - 1, b)
            lts = g.left_targets(b + 1)
            real = (b in lts) and (b - 1 in lts)
            # the residue b never appears on either side of any supplied pair,
            # which is the structural reason; recorded as a derived fact
            touch_b = sorted(set([u for (u, v) in sup if v == b and u != b] +
                                 [v for (u, v) in sup if u == b and v != b]))
            if (not holds) or (not real) or derivable or touch_b:
                if len(bad) < 6:
                    bad.append((a, b, g.q, holds, real, derivable, touch_b))
            if shown < 5:
                say("a=%d b=%d q=%d : derived G_%d<G_%d=%s ; both are Left options "
                    "of heap %d: %s ; supplied closure contains G_%d<=G_%d: %s ; "
                    "supplied pairs touching residue b other than (b,b): %s"
                    % (a, b, g.q, b - 1, b, holds, b + 1, real, b - 1, b,
                       derivable, touch_b))
                shown += 1
    say("parameter pairs in the q-does-not-divide-(a-1) regime: %d" % n_case)
    for v in bad:
        say("VIOLATION a=%d b=%d q=%d lt=%s both_options=%s derivable=%s "
            "pairs_touching_b=%s" % v)
    record(n_case > 0 and not bad,
           "whenever q does not divide a-1 the domination at heap b+1 is real "
           "(G_{b-1}<G_b, both Left options) and the comparison it needs is NOT in "
           "the transitive closure of Obs 5.2 + Prop 7.2 -- the refutation survives "
           "the paper's own transitivity-and-periodicity latitude",
           "%d parameter pairs, closure computed on {0..b} for each" % n_case)


# ----------------------------------------------------------------------------
# CHECK 7b -- THE POSITIVE VERDICT (the other half of the paper's iff, tex line
# 234-237): when q | (a-1), Obs 5.2 and Prop 7.2 together with transitivity and
# periodicity "account for every dominated option at every heap size".
# ----------------------------------------------------------------------------

def check_positive_verdict():
    print("[CHECK 7b] when q | (a-1): every dominated option at every heap "
          "n <= %d(b+1) is derivable from Obs 5.2 + Prop 7.2 by transitivity and "
          "periodicity, 0 < a < b <= %d" % (DERIV_REACH, DERIV_BMAX))
    bad = []
    red_bad = []
    n_case = 0
    n_nontrivial = 0
    n_dom = 0
    n_just = 0
    shown = 0
    for b in range(2, DERIV_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            if not divides(g.q, a - 1):
                continue
            n_case += 1
            if a >= 2:
                n_nontrivial += 1
            reach = transitive_closure(supplied_comparisons(g), b)
            seen_here = 0
            for n in range(0, DERIV_REACH * g.p + 1):
                for v in reduction_violations(g, n):
                    if len(red_bad) < 6:
                        red_bad.append((a, b) + v)
                for side in ('L', 'R'):
                    wit = domination_witnesses(g, n, side)
                    dominated = sorted(set(u for (u, _) in wit))
                    for u in dominated:
                        n_dom += 1
                        seen_here += 1
                        ok = any(reach[x][y] for (w, (x, y)) in wit if w == u)
                        if ok:
                            n_just += 1
                        else:
                            if len(bad) < 6:
                                bad.append((a, b, n, side, u,
                                            [xy for (w, xy) in wit if w == u]))
            if shown < 5:
                say("a=%d b=%d q=%d : q|(a-1)=True, final block = %s, dominated "
                    "options found over n<=%d = %d, all derivable so far: %s"
                    % (a, b, g.q, block(block_index(b, g.q), g.q, b),
                       DERIV_REACH * g.p, seen_here, not bad))
                shown += 1
    say("parameter pairs with q|(a-1): %d (of which a>=2: %d) ; dominated options "
        "encountered: %d ; derivable from the source: %d"
        % (n_case, n_nontrivial, n_dom, n_just))
    for v in red_bad:
        say("VIOLATION a=%d b=%d n=%d side=%s target %d not equal to residue %d" % v)
    for v in bad:
        say("VIOLATION a=%d b=%d n=%d side=%s dominated residue %d has no derivable "
            "witness; witnesses were %s" % v)
    record(not red_bad, "CHECK 7b's mod-(b+1) reductions all verified against the "
           "engine", "%d parameter pairs" % n_case)
    record(n_case > 0 and n_nontrivial > 0 and n_dom > 0 and not bad,
           "positive half of the iff: when q | (a-1) every dominated option at every "
           "heap size is accounted for by Obs 5.2 + Prop 7.2 under transitivity and "
           "periodicity",
           "%d parameter pairs (%d with a>=2), %d/%d dominations derivable"
           % (n_case, n_nontrivial, n_just, n_dom))


# ----------------------------------------------------------------------------
# CHECK 7c -- THE CONTEXTUAL VERDICT (abstract, sentence 4): for the b+1
# representative positions the source's comparisons ARE complete.  Stated as
# domination coverage, not merely as agreement on consecutive comparisons.
# Also certifies the paper's remark at tex 211-213: without the transitivity
# latitude the answer is negative for every a >= 2, because some domination is
# witnessed only inside the first block, i.e. by Obs 5.2 and never by Prop 7.2.
# ----------------------------------------------------------------------------

def check_contextual_verdict():
    print("[CHECK 7c] representative positions 0 <= n <= b: every dominated option "
          "is derivable from Obs 5.2 + Prop 7.2, 0 < a < b <= %d" % DERIV_BMAX)
    bad = []
    trivial_bad = []
    n_dom = 0
    n_pairs = 0
    n_a2 = 0
    shown = 0
    for b in range(2, DERIV_BMAX + 1):
        for a in range(1, b):
            g = game(a, b)
            n_pairs += 1
            sup = supplied_comparisons(g)
            reach = transitive_closure(sup, b)
            prop_pairs = set((n, n + 1) for n in prop72_supplied_indices(a, b))
            needs_obs = False
            for n in range(0, b + 1):
                for side in ('L', 'R'):
                    wit = domination_witnesses(g, n, side)
                    for u in sorted(set(w for (w, _) in wit)):
                        n_dom += 1
                        xys = [xy for (w, xy) in wit if w == u]
                        if not any(reach[x][y] for (x, y) in xys):
                            if len(bad) < 6:
                                bad.append((a, b, n, side, u, xys))
                        if all(xy not in prop_pairs for xy in xys):
                            needs_obs = True
            if a >= 2:
                n_a2 += 1
                if not needs_obs:
                    if len(trivial_bad) < 6:
                        trivial_bad.append((a, b, g.q))
            if shown < 4:
                say("a=%d b=%d q=%d : Prop 7.2 supplies pairs %s ; representative "
                    "heaps 0..%d examined ; some domination needs Obs 5.2 rather "
                    "than Prop 7.2: %s"
                    % (a, b, g.q, sorted(prop_pairs), b, needs_obs))
                shown += 1
    say("parameter pairs = %d ; dominations at representative heaps = %d"
        % (n_pairs, n_dom))
    for v in bad:
        say("VIOLATION a=%d b=%d n=%d side=%s dominated residue %d not derivable; "
            "witnesses %s" % v)
    record(n_pairs > 0 and n_dom > 0 and not bad,
           "contextual reading of Problem 7.1 is affirmative: every dominated option "
           "of every representative position 0..b is derivable from Obs 5.2 + Prop 7.2 "
           "by transitivity and periodicity",
           "%d parameter pairs, %d dominations, all derivable" % (n_pairs, n_dom))
    for v in trivial_bad:
        say("VIOLATION a=%d b=%d q=%d : no domination required Obs 5.2" % v)
    record(n_a2 > 0 and not trivial_bad,
           "for every a >= 2 some domination is witnessed only by a comparison "
           "Prop 7.2 does not supply (it lies in the first block, i.e. comes from "
           "Obs 5.2), so without the transitivity latitude the answer is negative -- "
           "the paper's remark at the head of section 3",
           "%d parameter pairs with a>=2" % n_a2)


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    print("verify.py -- verifier for 'A Counterexample to the Unrestricted "
          "Dargad--Larsson Domination Problem'")
    print("game: G_n = TS_{b-a}(n;a,b), S_L={1..a}, S_R={q..b}, q=b-a+1, "
          "p=b+1; every comparison is computed from the recursion and the")
    print("short-game criterion X<=Y iff (no X^L with Y<=X^L) and "
          "(no Y^R with Y^R<=X).  Integer arithmetic only; no floating point,")
    print("so no check is numeric or tolerance-based.")
    print("")
    check_witness()
    print("")
    check_theorem1_order()
    print("")
    check_periodicity()
    print("")
    check_boundary()
    print("")
    check_a1()
    print("")
    check_a2()
    print("")
    check_domination_corollary()
    print("")
    check_prop72_consecutive()
    print("")
    check_Q_two_forms()
    print("")
    check_obs52()
    print("")
    check_negative_verdict()
    print("")
    check_positive_verdict()
    print("")
    check_contextual_verdict()
    print("")
    print("scope note: the sweeps above are exhaustive over the stated finite ranges "
          "of (a,b) (b <= %d for the order," % THM1_BMAX)
    print("b <= %d for periodicity and least period, b < %d for the boundary "
          "dichotomy, b <= %d for the derivability closure of CHECK 7)."
          % (PERIOD_BMAX, BOUNDARY_BMAX + 1, DERIV_BMAX))
    print("The paper's claim for all 0<a<b rests on its Theorem 1 hand proof; this "
          "program certifies the computational")
    print("content, not the induction beyond those ranges.")
    print("")
    print("gaps this program does NOT close, stated rather than hidden:")
    print("  (1) QUOTATION FIDELITY.  That Problem 7.1, Proposition 7.2 and "
          "Observation 5.2 are numbered and worded in")
    print("      arXiv:2607.27989v1 as the paper quotes them -- in particular that "
          "Prop 7.2's printed comparison range")
    print("      really ends at b-2, that Q really is {(m+1)q-1 : m>=1}, and that "
          "Problem 7.1 as printed carries no")
    print("      restriction on the heap size n -- is input data here, not a verified "
          "fact; the identifier itself is not")
    print("      resolved by this program either.  Every verdict of CHECK 1d, "
          "CHECK 5b and CHECK 7 is conditional on")
    print("      those three transcriptions being right; a single typo in any of them "
          "would invert the negative verdict.")
    print("  (2) UNRESTRICTED (a,b).  Nothing above is an induction; the finite "
          "sweeps cannot certify all 0<a<b.")
    print("  (3) DECORATION, declared.  These checks state true but definitional "
          "facts and could not fail for any input:")
    print("      CHECK 3a-iii (q|b <=> q|(a-1) is forced by b=q+a-1); CHECK 3a-iv "
          "(implied by CHECK 3a-ii on the same")
    print("      parameters); CHECK 3b-i's empty dominated sets (a=1 gives at most "
          "one option per player, so emptiness is")
    print("      structural); CHECK 3b-ii (q divides 0 always); CHECK 6b (B_0 ends "
          "at q-1 and tau=q-1 by definition);")
    print("      the 'spread = a-1' clause of CHECK 4a.  They are kept for the "
          "record, but the load-bearing checks are")
    print("      1a-1c, 2a, 2b, 3a-i/ii, 3c-ii/iii, 4a's residue-collision clause, "
          "4b, 5a, 5b, 6a and 7a-7c.")
    print("  (4) CHECK 4b is a consequence of CHECK 2a, not independent evidence: "
          "once the derived order is known to be")
    print("      'same block and r<=s', the blockwise-extremum domination rule "
          "follows algebraically.")
    print("")
    total = len(_RESULTS)
    failed = [lab for ok, lab in _RESULTS if not ok]
    if failed:
        for lab in failed:
            print("failed check: " + lab)
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


if __name__ == "__main__":
    sys.exit(main())

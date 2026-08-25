#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- independent verifier for

    "A Counterexample to He and Huang's Runge-Kutta Necessity Conjecture"
    (the embedded order-five formula of Verner's eight-stage pair RK(8-6:5)a
     violates QO(m) for every admissible (m,n), refuting Conjecture 2.1 of
     arXiv:2605.16995v3)

Python 3.9+, standard library only.  All arithmetic is exact: every number in
this program is a Python int or a fractions.Fraction.  There is NO floating
point anywhere, so no error bound is required and none is printed.

--------------------------------------------------------------------------
TAKEN FROM THE PAPER (input data -- transcribed, not derived)
--------------------------------------------------------------------------
  T1. The 8x8 coefficient matrix A of Verner's RK(8-6:5)a, as printed in the
      tableau of Theorem 2 (56 rational entries, strictly lower triangular).
  T2. The abscissa column c printed to the left of that tableau.
  T3. The embedded weight vector b_hat printed below the rule of that tableau.
  T4. The definitions quoted from He-Huang: q_r = A c^{or} - c^{o(r+1)}/(r+1),
      Q_1 = span{q_0}, Q_2 = span{q_1, A q_0, q_0 o c, q_0}, Q_1 <= Q_2 <= ...,
      QO(m): b o Q_m = {0} (Hadamard), and the admissibility constraints
      m >= n-1, m+n+1 >= p on nonnegative integers m, n.
  T5. Butcher's tree criterion as quoted in Theorem 2: Phi(.) = 1,
      Phi([t_1..t_k])_i = prod_j (A Phi(t_j))_i, gamma(.) = 1,
      gamma(t) = |t| prod_j gamma(t_j), order p  <=>  b^T Phi(t) = 1/gamma(t)
      for all rooted trees with |t| <= p.
  T6. The paper's *claimed numerical values*, held only as comparison targets
      (never used as inputs to a derivation): c itself as printed, Ac,
      q_1, b_hat o q_1, b_hat^T A^4 c = 1/672, the tree counts 1,1,2,4,9,20,
      "eight violations at order six", "min m = 2", the stage index set
      {1,3,4,5,7,8}, b_hat_2 = -1/5, c_2 = 1/5.

--------------------------------------------------------------------------
DERIVED HERE (what the checks actually decide)
--------------------------------------------------------------------------
  D1. That A is strictly lower triangular and hence c_1 = (A1)_1 = 0.
  D2. c = A*1 recomputed row by row from A alone, compared with the printed c.
  D3. sum(b_hat), compared with 1.
  D4. q_0 = A*1 - c, hence Q_1 = span{0} = {0}: the obstruction genuinely
      needs m >= 2, it is not available at m = 1.
  D5. Ac by matrix-vector product; q_1 = Ac - c^{o2}/2 componentwise.
  D6. The exact set of stages at which (Ac)_i = c_i^2/2, and the values of
      c_i and (Ac)_i at the remaining stages.
  D7. b_hat o q_1 componentwise, and whether it is the zero vector.
  D8. b_hat o g for every printed generator g of Q_2, hence whether QO(2)
      holds; and by Q_2 <= Q_m, whether QO(m) can hold for any m >= 2.
  D9. b_hat_2 * c_2^2, the quantity the paper's Proposition forces to vanish.
 D10. A^k c for k = 1..4 by repeated exact matrix-vector products, and the
      four chain conditions b_hat^T A^{k} c against 1/6, 1/24, 1/120, 1/720.
 D11. The complete set of rooted trees with at most 6 nodes, generated from
      scratch by canonical multiset recursion (no count is assumed), their
      densities gamma(t) and elementary weights Phi(t) over the Fractions.
 D12. Every order condition residual b_hat^T Phi(t) - 1/gamma(t) for |t| <= 6,
      the number of violations at each order, and thus that the order is
      exactly five.
 D13. b_hat^T q_1 (the Scope paragraph) and the general identity
      b^T q_1 = b^T Ac - (1/2) b^T c^{o2} instantiated on this method.
 D14. The Proposition's key algebraic step (q_1)_2 = -c_2^2/2, re-derived on a
      deterministic family of unrelated internally consistent explicit
      tableaux, to test the lemma and not just the counterexample.
 D15. min{m : exists n >= 0, m >= n-1, m+n+1 >= p} by integer enumeration at the
      DERIVED order p, cross-checked against the closed form 2m+2 >= p (which is
      what removes the finite-box truncation), plus the half-integer relaxation
      showing that integrality of m is load-bearing.
 D16. A falsification arm for D14: one generator, two classes.  Inside the
      Proposition's hypothesis the lemma never fails; outside it, all three
      lemma statements fail on explicitly exhibited tableaux.  Without this,
      D14 is 74799 iterations of an identity that cannot fail.
 D17. The quadrature conditions B(k): b_hat^T c^{o(k-1)} = 1/k for every
      k = 1..p at the DERIVED order p, computed directly from b_hat and c
      (the bush trees of D12 impose the same conditions by a different route).
      B(p) is therefore the one member of the Theorem 1.2 package besides
      QO(m) that this program decides, and the GAPS block's list of members it
      does NOT evaluate is derived from these checks rather than asserted.
 D18. An explicit GAPS block (printed last) naming the seven steps between the
      arithmetic above and the paper's conclusion that NO check here decides.

WHAT THIS PROGRAM CANNOT DECIDE -- see the GAPS block at the end of the run.
The definitions of q_r, of Q_2's generators, of the chain Q_1 <= Q_2 <= ...,
the Hadamard reading of QO(m), and the admissibility inequalities are all
transcribed from the paper's quotation of He-Huang.  Every check below passes
if any of them is misquoted.
"""

import sys
import textwrap
from fractions import Fraction as F
from itertools import product as _iproduct

# ---------------------------------------------------------------------------
# T1/T2/T3: the tableau of Theorem 2, transcribed verbatim.
# A[i][j] is a_{i+1,j+1}; the method is explicit so A is strictly lower
# triangular.  Trailing zeros are written out in full so the transcription can
# be diffed against the printed array shape.
# ---------------------------------------------------------------------------
Z = F(0)

A = [
    [Z,               Z,             Z,               Z,              Z,                Z,           Z, Z],
    [F(1, 5),         Z,             Z,               Z,              Z,                Z,           Z, Z],
    [F(3, 32),        F(5, 32),      Z,               Z,              Z,                Z,           Z, Z],
    [F(60, 343),      F(-440, 343),  F(576, 343),     Z,              Z,                Z,           Z, Z],
    [F(847, 17496),   F(770, 729),   F(-17024, 19683), F(84721, 157464), Z,             Z,           Z, Z],
    [F(523, 2240),    Z,             F(-5, 57),       F(245, 2496),   F(-1215, 27664),  Z,           Z, Z],
    [F(15, 2212),     F(314, 395),   F(3448, 1501),   F(-2695, 4108), F(203391, 273182), F(-864, 395), Z, Z],
    [F(1185, 10976),  F(-2543, 1715), F(36262, 19551), F(-245, 1248), F(1215, 13832),   F(1, 5),     Z, Z],
]

# The abscissae as printed in the left column of the tableau.
C_PAPER = [Z, F(1, 5), F(1, 4), F(4, 7), F(7, 9), F(1, 5), F(1), F(4, 7)]

# The embedded (order five) weight vector printed below the rule.
BHAT = [F(43, 560), F(-1, 5), F(2816, 7695), F(-41, 84240),
        F(19683, 69160), F(1, 5), F(79, 1080), F(1, 5)]

N = 8

# ---------------------------------------------------------------------------
# T6: values the paper claims.  These are comparison targets ONLY.  Nothing
# below ever reads one of these to build the quantity it is compared against.
# ---------------------------------------------------------------------------
PAPER_AC = [Z, Z, F(1, 32), F(8, 49), F(49, 162), Z, F(1, 2), F(8, 49)]
PAPER_Q1 = [Z, F(-1, 50), Z, Z, Z, F(-1, 50), Z, Z]
PAPER_B_HAD_Q1 = [Z, F(1, 250), Z, Z, Z, F(-1, 250), Z, Z]
PAPER_EQ_STAGES = (1, 3, 4, 5, 7, 8)      # 1-based stages with (Ac)_i = c_i^2/2
PAPER_ODD_C = F(1, 5)                     # c_i at the two remaining stages
PAPER_B2 = F(-1, 5)                       # abstract: b_hat_2
PAPER_C2 = F(1, 5)                        # abstract: c_2
PAPER_CHAIN = {3: F(1, 6), 4: F(1, 24), 5: F(1, 120), 6: F(1, 720)}
PAPER_B_A4_C = F(1, 672)                  # the failing order-six chain value
PAPER_TREE_COUNTS = {1: 1, 2: 1, 3: 2, 4: 4, 5: 9, 6: 20}
PAPER_ORDER6_VIOLATIONS = 8
PAPER_EXACT_ORDER = 5
PAPER_MIN_M = 2
PAPER_P = 5

_RESULTS = []      # list of (bool ok, str message)


def check(ok, label, detail=""):
    """Record and immediately print one check line in the fixed format."""
    ok = bool(ok)
    _RESULTS.append((ok, label))
    line = ("PASS " if ok else "FAIL ") + label
    if detail:
        line += "  [" + detail + "]"
    print(line)
    return ok


def say(text):
    """Print a derived intermediate quantity (not itself a check)."""
    print("      " + text)


def vstr(v):
    """Render a vector of Fractions compactly and exactly."""
    return "(" + ", ".join(str(x) for x in v) + ")"


def matvec(M, v):
    """Exact M*v for a square list-of-lists M and list v."""
    n = len(M)
    out = []
    for i in range(n):
        s = F(0)
        row = M[i]
        for j in range(n):
            if row[j]:
                s += row[j] * v[j]
        out.append(s)
    return out


def dot(u, v):
    """Exact u^T v."""
    s = F(0)
    for a, b in zip(u, v):
        s += a * b
    return s


def hadamard(u, v):
    """Exact componentwise product u o v."""
    return [a * b for a, b in zip(u, v)]


def is_zero(v):
    return all(x == 0 for x in v)


def check_tableau_structure():
    """D1: A strictly lower triangular (explicit method), hence c_1 = 0.

    Also guards the shapes: dot() and hadamard() use zip(), which SILENTLY
    truncates to the shorter operand, so a dropped entry in b_hat or in a row
    of A would otherwise produce a wrong scalar with no complaint.
    """
    shapes = ([len(A)] + [len(r) for r in A] + [len(C_PAPER), len(BHAT),
              len(PAPER_AC), len(PAPER_Q1), len(PAPER_B_HAD_Q1)])
    check(all(s == N for s in shapes),
          "every transcribed array has exactly %d entries (zip cannot truncate)" % N,
          "shapes = %s" % (shapes,))
    bad = [(i + 1, j + 1) for i in range(N) for j in range(i, N) if A[i][j] != 0]
    check(not bad, "A is strictly lower triangular (method is explicit)",
          "offending entries: %s" % (bad,) if bad else "all 36 on/above-diagonal entries vanish")
    row1 = A[0]
    check(is_zero(row1), "row 1 of A vanishes, so (A*1)_1 = 0 and c_1 = 0",
          "A[1,:] = " + vstr(row1))
    return not bad


def check_internal_consistency():
    """D2/D3: recompute c = A*1 from A, and sum(b_hat), and compare to print."""
    ones = [F(1)] * N
    c_derived = matvec(A, ones)
    say("c = A*1 derived from A : " + vstr(c_derived))
    say("c as printed in tableau: " + vstr(C_PAPER))
    rows_ok = [c_derived[i] == C_PAPER[i] for i in range(N)]
    for i in range(N):
        if not rows_ok[i]:
            say("  row %d MISMATCH: (A*1)_%d = %s but printed c_%d = %s"
                % (i + 1, i + 1, c_derived[i], i + 1, C_PAPER[i]))
    check(all(rows_ok),
          "internal consistency c = A*1 holds on all %d rows" % N,
          "%d/%d rows agree" % (sum(1 for r in rows_ok if r), N))
    s = sum(BHAT, F(0))
    say("sum(b_hat) derived = " + str(s))
    check(s == 1, "sum(b_hat) = 1 (quadrature/order-one condition B(1))",
          "derived %s, expected 1" % s)
    return c_derived


def check_q0(c):
    """D4: q_0 = A c^{o0} - c^{o1}/1, evaluated against the PRINTED abscissae.

    HONESTY NOTE (was a tautology before).  Computing q_0 as
    matvec(A, ones) - c with c itself defined as matvec(A, ones) is zero for
    ANY matrix whatsoever, so such a check cannot fail and decides nothing.
    Here q_0 is formed against C_PAPER, the abscissa column as PRINTED, so a
    mistyped abscissa really does break this line.

    CONVENTION NOTE.  c^{o0} is read as the all-ones vector, i.e. 0^0 = 1 in
    component 1.  Under the opposite convention q_0 = -A[:,1] != 0 and Q_1
    would be nonzero.  Nothing downstream depends on the choice: the smallest
    admissible m at p = 5 is 2 (section 12), so the witness only ever has to
    live in Q_2, never in Q_1.  This item is therefore informative, not
    load-bearing.
    """
    ones = [F(1)] * N
    q0 = [x - y for x, y in zip(matvec(A, ones), C_PAPER)]
    say("q_0 = A c^{o0} - c^{o1}/1 = A*1 - c(printed) = " + vstr(q0))
    ok = check(is_zero(q0),
               "q_0 = 0 against the PRINTED c, hence Q_1 = span{q_0} = {0}",
               "0^0 = 1 convention; not load-bearing, since min admissible m = 2")
    # keep the derived-c copy for the Q_2 generator sweep
    q0d = [x - y for x, y in zip(matvec(A, ones), c)]
    return q0d, ok


def check_q1(c):
    """D5/D6: Ac, then q_1 = Ac - c^{o2}/2, and the stage set (Ac)_i = c_i^2/2."""
    ac = matvec(A, c)
    say("Ac derived                : " + vstr(ac))
    say("Ac as printed in the paper: " + vstr(PAPER_AC))
    check(ac == PAPER_AC, "derived Ac equals the vector printed in the paper",
          "derived " + vstr(ac))
    csq = hadamard(c, c)
    say("c^{o2} derived             : " + vstr(csq))
    q1 = [ac[i] - csq[i] / 2 for i in range(N)]
    say("q_1 = Ac - c^{o2}/2 derived: " + vstr(q1))
    say("q_1 as printed in the paper: " + vstr(PAPER_Q1))
    check(q1 == PAPER_Q1,
          "derived q_1 equals the printed (0,-1/50,0,0,0,-1/50,0,0)",
          "derived " + vstr(q1))
    # D6: which stages satisfy (Ac)_i = c_i^2 / 2, i.e. (q_1)_i = 0?
    eq = tuple(i + 1 for i in range(N) if ac[i] == csq[i] / 2)
    ne = tuple(i + 1 for i in range(N) if ac[i] != csq[i] / 2)
    say("stages with (Ac)_i = c_i^2/2 : %s   (remaining: %s)" % (eq, ne))
    check(eq == PAPER_EQ_STAGES,
          "the stage set where (Ac)_i = c_i^2/2 is exactly {1,3,4,5,7,8}",
          "derived %s, paper says %s" % (eq, PAPER_EQ_STAGES))
    detail = "; ".join("stage %d: c=%s, (Ac)=%s" % (i, c[i - 1], ac[i - 1]) for i in ne)
    check(all(c[i - 1] == PAPER_ODD_C and ac[i - 1] == 0 for i in ne) and len(ne) == 2,
          "the two remaining stages both have c_i = 1/5 and (Ac)_i = 0", detail)
    return q1, ac, csq


def check_hadamard_violation(q1):
    """D7: b_hat o q_1, the actual violation of QO(m)."""
    h = hadamard(BHAT, q1)
    say("b_hat o q_1 derived                : " + vstr(h))
    say("b_hat o q_1 as printed in the paper: " + vstr(PAPER_B_HAD_Q1))
    check(h == PAPER_B_HAD_Q1,
          "derived b_hat o q_1 equals the printed (0,1/250,0,0,0,-1/250,0,0)",
          "derived " + vstr(h))
    nz = [(i + 1, h[i]) for i in range(N) if h[i] != 0]
    check(not is_zero(h),
          "b_hat o q_1 != 0, so q_1 witnesses the failure of QO(m)",
          "nonzero components: " + ", ".join("(b o q_1)_%d = %s" % t for t in nz))
    return h


def check_QO2_generators(q0, q1, c):
    """D8: b_hat o g over the printed spanning set of Q_2, hence QO(m), m >= 2."""
    gens = [("q_1", q1), ("A q_0", matvec(A, q0)),
            ("q_0 o c", hadamard(q0, c)), ("q_0", q0)]
    nonzero_gens = []
    for name, g in gens:
        h = hadamard(BHAT, g)
        say("generator %-8s = %s" % (name, vstr(g)))
        say("   b_hat o %-8s = %s%s" % (name, vstr(h), "" if is_zero(h) else "   <-- NONZERO"))
        if not is_zero(h):
            nonzero_gens.append(name)
    # b o span{g_1..g_k} = {0}  <=>  b o g_i = 0 for every i, since b o (.) is linear.
    check(nonzero_gens == ["q_1"],
          "QO(2) fails, and q_1 is the only printed Q_2 generator that breaks it",
          "generators with b_hat o g != 0: %s" % (nonzero_gens,))
    # What the arithmetic decides is QO(2) alone.  The step to "every m >= 2"
    # uses the quoted chain Q_2 <= Q_m, which is an ASSUMED input (gap G2), not
    # something this program can evaluate, so the label says only that.
    check(len(nonzero_gens) > 0,
          "QO(2) is decided FALSE by arithmetic; the extension to every m >= 2"
          " rests on the quoted chain Q_2 <= Q_m (see gap G2)",
          "witness q_1 in Q_2 with b_hat o q_1 != 0")


def check_proposition_quantity(c, q1):
    """D9: the quantity the Proposition forces to vanish is nonzero here."""
    b2, c2 = BHAT[1], c[1]
    say("b_hat_2 derived from the printed weight vector = " + str(b2))
    say("c_2 derived from A*1                           = " + str(c2))
    check(b2 == PAPER_B2 and c2 == PAPER_C2,
          "abstract's values b_hat_2 = -1/5 and c_2 = 1/5 are the tableau's",
          "derived b_hat_2 = %s, c_2 = %s" % (b2, c2))
    prod = b2 * c2 * c2
    say("b_hat_2 * c_2^2 = " + str(prod))
    check(prod != 0,
          "b_2 c_2^2 != 0, contradicting the Proposition's conclusion under QO(m)",
          "derived b_hat_2 c_2^2 = %s" % prod)
    check(q1[1] == -c2 * c2 / 2,
          "the Proposition's step (q_1)_2 = -c_2^2/2 holds on this tableau",
          "(q_1)_2 = %s and -c_2^2/2 = %s" % (q1[1], -c2 * c2 / 2))


def check_chain_conditions(c):
    """D10: b_hat^T A^k c for k = 1..4 against 1/6, 1/24, 1/120, 1/720.

    The rooted chain (path) on k+2 nodes has Phi = A^k c and gamma = (k+2)!,
    so these are the order-(k+2) chain conditions of the paper.
    """
    v = c[:]                       # A^0 c
    vals = {}
    for k in range(1, 5):
        v = matvec(A, v)           # now A^k c
        order = k + 2
        val = dot(BHAT, v)
        vals[order] = val
        say("A^%d c = %s" % (k, vstr(v)))
        say("   b_hat^T A^%d c = %s   (order-%d chain, target 1/gamma = %s)"
            % (k, val, order, PAPER_CHAIN[order]))
    for order in (3, 4, 5):
        check(vals[order] == PAPER_CHAIN[order],
              "order-%d chain condition b_hat^T A^%d c = %s holds"
              % (order, order - 2, PAPER_CHAIN[order]),
              "derived %s" % vals[order])
    check(vals[6] == PAPER_B_A4_C,
          "derived b_hat^T A^4 c = 1/672, the paper's printed value",
          "derived %s" % vals[6])
    check(vals[6] != PAPER_CHAIN[6],
          "order-six chain condition FAILS: b_hat^T A^4 c != 1/720",
          "residual = %s - %s = %s"
          % (vals[6], PAPER_CHAIN[6], vals[6] - PAPER_CHAIN[6]))
    return vals[6] - PAPER_CHAIN[6]


# ---------------------------------------------------------------------------
# D11: rooted trees.  A tree is the tuple of its child subtrees, kept sorted so
# that the tuple is a canonical form for the unordered rooted tree; () is the
# single node "bullet".  Isomorphic trees therefore collapse to one object in a
# set, which is what makes the derived counts trustworthy: nothing about the
# number of trees is assumed anywhere.
# ---------------------------------------------------------------------------

def _forests(total, trees_by_order, memo):
    """All canonical multisets of trees whose node counts sum to `total`."""
    if total in memo:
        return memo[total]
    if total == 0:
        memo[0] = [()]
        return memo[0]
    res = set()
    for size in range(1, total + 1):
        for t in trees_by_order[size]:
            for rest in _forests(total - size, trees_by_order, memo):
                res.add(tuple(sorted((t,) + rest)))
    memo[total] = sorted(res)
    return memo[total]


def gen_trees(max_order):
    """trees_by_order[k] = every rooted tree on k nodes, up to isomorphism."""
    trees_by_order = {1: [()]}
    memo = {}
    for n in range(2, max_order + 1):
        # a tree on n nodes is a root plus a forest of n-1 nodes
        trees_by_order[n] = sorted(set(_forests(n - 1, trees_by_order, memo)))
    return trees_by_order


def tree_order(t):
    """|t|, the number of nodes."""
    return 1 + sum(tree_order(s) for s in t)


def tree_str(t):
    """Bracket notation: '.' for the single node, '[..]' for a root's children."""
    if not t:
        return "."
    return "[" + "".join(tree_str(s) for s in t) + "]"


def gamma(t):
    """Density: gamma(.) = 1, gamma(t) = |t| * prod_j gamma(t_j).  Exact int."""
    g = tree_order(t)
    for s in t:
        g *= gamma(s)
    return g


def phi(t):
    """Elementary weight vector: Phi(.) = 1, Phi([t_1..t_k])_i = prod_j (A Phi(t_j))_i."""
    if not t:
        return [F(1)] * N
    out = [F(1)] * N
    for s in t:
        w = matvec(A, phi(s))
        out = [out[i] * w[i] for i in range(N)]
    return out


def check_tree_counts(trees_by_order, max_order):
    """D11: the derived tree census against the paper's 1,1,2,4,9,20."""
    derived = {k: len(v) for k, v in trees_by_order.items()}
    say("rooted-tree census derived: " +
        ", ".join("order %d: %d" % (k, derived[k]) for k in sorted(derived)))
    say("paper's census            : " +
        ", ".join("order %d: %d" % (k, PAPER_TREE_COUNTS[k])
                  for k in sorted(PAPER_TREE_COUNTS)))
    check(derived == PAPER_TREE_COUNTS,
          "derived rooted-tree counts are 1,1,2,4,9,20 at orders 1..6",
          "derived %s" % ([derived[k] for k in sorted(derived)],))
    n_low = sum(derived[k] for k in range(1, max_order))
    check(n_low == 17,
          "the orders 1..5 carry exactly 17 order conditions",
          "derived %d = 1+1+2+4+9" % n_low)
    return derived


def sweep_order_conditions(trees_by_order, max_order):
    """D12: residual b_hat^T Phi(t) - 1/gamma(t) for every tree, |t| <= max_order."""
    viol = {}
    print("      --- order conditions, one line per rooted tree ---")
    for k in range(1, max_order + 1):
        viol[k] = []
        for t in trees_by_order[k]:
            g = gamma(t)
            lhs = dot(BHAT, phi(t))
            res = lhs - F(1, g)
            if res != 0:
                viol[k].append((t, lhs, g, res))
            say("order %d  t = %-14s gamma = %-4d b^T Phi = %-14s 1/gamma = %-10s %s"
                % (k, tree_str(t), g, str(lhs), str(F(1, g)),
                   "OK" if res == 0 else "VIOLATED by " + str(res)))
    print("      --- end of order conditions ---")
    for k in range(1, max_order + 1):
        say("order %d: %d trees, %d violations" % (k, len(trees_by_order[k]), len(viol[k])))
    return viol


def check_exact_order(trees_by_order, viol, max_order, chain6_residual):
    """D12: zero violations up to order five, exactly eight at order six."""
    low = sum(len(viol[k]) for k in range(1, max_order))
    check(low == 0,
          "all 17 order conditions of orders 1..5 hold exactly (0 violations)",
          "violations found at orders 1..5: %d" % low)
    n6 = len(viol[max_order])
    check(n6 == PAPER_ORDER6_VIOLATIONS,
          "exactly 8 of the 20 order-six conditions fail",
          "derived %d violations out of %d trees"
          % (n6, len(trees_by_order[max_order])))
    check(low == 0 and n6 > 0,
          "the embedded formula has order exactly %d" % PAPER_EXACT_ORDER,
          "orders 1..5 all satisfied, order 6 not satisfied")
    # cross-link: the order-six chain must be one of the eight failures, with
    # the residual already computed independently as b^T A^4 c - 1/720.
    chain = ()
    for _ in range(max_order - 1):
        chain = (chain,)
    hit = [v for v in viol[max_order] if v[0] == chain]
    check(len(hit) == 1 and hit[0][3] == chain6_residual,
          "the order-six chain %s is among the failures, with the residual"
          " already derived from b_hat^T A^4 c" % tree_str(chain),
          "tree residual %s vs chain residual %s"
          % (hit[0][3] if hit else None, chain6_residual))
    return n6


def check_scope(c, q1, ac, csq):
    """D13: b_hat^T q_1 = 0, and the general third-order argument behind it."""
    s = dot(BHAT, q1)
    say("b_hat^T q_1 derived = " + str(s))
    check(s == 0,
          "b_hat^T q_1 = 0: scalar orthogonality is NOT violated (Scope)",
          "derived %s, so the refutation is specific to the Hadamard reading" % s)
    b_ac = dot(BHAT, ac)          # third-order condition b^T Ac
    b_csq = dot(BHAT, csq)        # third-order condition b^T c^{o2}
    say("b_hat^T Ac = %s   b_hat^T c^{o2} = %s" % (b_ac, b_csq))
    check(b_ac == F(1, 6) and b_csq == F(1, 3),
          "the two third-order conditions give b^T Ac = 1/6 and b^T c^{o2} = 1/3",
          "derived %s and %s" % (b_ac, b_csq))
    lhs = b_ac - b_csq / 2
    # HONESTY NOTE: the "lhs == s" conjunct is an ALGEBRAIC IDENTITY, not a
    # test.  Since q_1 = Ac - c^{o2}/2 and the dot product is bilinear,
    # b^T Ac - (1/2) b^T c^{o2} equals b^T q_1 for any b and any A over the
    # rationals, so that conjunct cannot fail and decides nothing.  The content
    # of this line is the "lhs == 0" conjunct, and that in turn is the same
    # statement as the b_hat^T q_1 = 0 check immediately above; what is new is
    # only that the value 0 is exhibited as 1/6 - (1/2)(1/3), i.e. as a
    # consequence of the two third-order conditions, which is the paper's
    # explanation for why scalar orthogonality survives.
    check(lhs == 0 and lhs == s,
          "the Scope identity b^T q_1 = b^T Ac - (1/2) b^T c^{o2} = 1/6 - 1/6 = 0"
          " (the second conjunct is bilinearity, not a test)",
          "derived %s, matching b_hat^T q_1 = %s" % (lhs, s))


def check_admissible_m(p, bound=40):
    """D15: min m over nonnegative integers with m >= n-1 and m+n+1 >= p.

    `p` is the DERIVED exact order, not a literal, so the integer argument is
    driven by the tree sweep rather than by a constant.

    The enumeration alone would only be a statement about the box [0,bound]^2;
    an unbounded n could in principle rescue a small m.  It cannot, and the
    reason is elementary: m >= n-1 caps n at m+1, so the largest attainable
    value of m+n+1 for a fixed m is 2m+2.  Hence

        m is admissible  <=>  2m+2 >= p  <=>  m >= ceil((p-2)/2),

    which is checked in closed form against the enumeration below.  That
    closed form is monotone in m, so no m outside the box can be smaller.
    """
    check(p == PAPER_P, "the order used in the integer argument is the DERIVED"
          " order, p = %d" % p, "derived p = %d, paper's p = %d" % (p, PAPER_P))
    pairs = [(m, n) for m in range(bound + 1) for n in range(bound + 1)
             if m >= n - 1 and m + n + 1 >= p]
    check(bool(pairs), "the admissible set for p = %d is nonempty" % p,
          "%d pairs (m,n) found with 0 <= m,n <= %d" % (len(pairs), bound))
    min_m = min(m for m, _ in pairs)
    witnesses = sorted(n for m, n in pairs if m == min_m)
    say("smallest admissible m derived = %d, attained at n in %s" % (min_m, witnesses))
    check(min_m == PAPER_MIN_M,
          "min m over admissible (m,n) at p = %d is 2" % p,
          "derived min m = %d (paper says %d)" % (min_m, PAPER_MIN_M))
    # closed form, which is what removes the box truncation
    enum_adm = set(m for m, _ in pairs)
    closed_adm = set(m for m in range(bound + 1) if 2 * m + 2 >= p)
    ceil_m = -((-(p - 2)) // 2)
    say("closed form: m admissible <=> 2m+2 >= %d <=> m >= ceil((p-2)/2) = %d"
        % (p, ceil_m))
    check(enum_adm == closed_adm and ceil_m == min_m,
          "the closed form 2m+2 >= p reproduces the enumeration exactly, so no"
          " m outside the box can be admissible and min m = 2 is unconditional",
          "enumeration and closed form agree on %d values of m; ceil = %d"
          % (len(closed_adm), ceil_m))
    # INTEGRALITY IS LOAD-BEARING -- demonstrate it instead of asserting it.
    # (The paper attributes "nonnegative integers m,n" to Theorem 1.2; the
    #  source writes only "for some m, n".  See the GAPS block.)
    grid = [F(k, 2) for k in range(2 * bound + 1)]
    hpairs = [(m, n) for m in grid for n in grid
              if m >= n - 1 and m + n + 1 >= p]
    min_half = min(m for m, _ in hpairs)
    say("min m over the HALF-INTEGER grid = %s (real optimum (p-2)/2 = %s)"
        % (min_half, F(p - 2, 2)))
    check(min_half == F(p - 2, 2) and min_half < PAPER_MIN_M,
          "integrality is LOAD-BEARING: dropping it lowers min m to 3/2 < 2,"
          " at which Q_m is not even defined",
          "half-integer min m = %s < %d" % (min_half, PAPER_MIN_M))
    return min_m


def check_proposition_lemma():
    """D14: test the Proposition's algebra, not just this one counterexample.

    For EVERY strictly lower triangular A with c = A*1 the paper claims
    c_1 = 0 and (q_1)_2 = (Ac)_2 - c_2^2/2 = a_21 c_1 - c_2^2/2 = -c_2^2/2,
    whence (b o q_1)_2 = -b_2 c_2^2 / 2 for any b, so QO(m) with m >= 2 forces
    b_2 c_2^2 = 0.  We re-derive this on an exhaustive deterministic family of
    small tableaux built from fixed rational value sets -- no floating point,
    no randomness, so the family is identical on every run.
    """
    families = [(3, (F(-2), F(-1, 3), F(0), F(1, 2), F(3))),
                (4, (F(-2), F(-1, 3), F(0), F(1, 2), F(3))),
                (5, (F(-5, 7), F(0), F(4, 3)))]
    total = 0
    bad_c1 = bad_q1 = bad_had = 0
    half = F(1, 2)
    for n, vals in families:
        # positions of the strictly lower triangular entries, row by row
        rows_idx = [list(range(i)) for i in range(n)]
        nfree = sum(len(r) for r in rows_idx)
        b_probe = [F(i + 2, 3) for i in range(n)]
        b2 = b_probe[1]
        for assign in _iproduct(vals, repeat=nfree):
            # unpack the flat assignment back into rows of A
            rows = []
            k = 0
            for i in range(n):
                w = len(rows_idx[i])
                rows.append(assign[k:k + w])
                k += w
            cc = [sum(r) if r else F(0) for r in rows]
            acc = []
            for i in range(n):
                s = F(0)
                for j, v in enumerate(rows[i]):
                    if v:
                        s += v * cc[j]
                acc.append(s)
            qq = [acc[i] - half * cc[i] * cc[i] for i in range(n)]
            total += 1
            if cc[0] != 0:
                bad_c1 += 1
            if qq[1] != -half * cc[1] * cc[1]:
                bad_q1 += 1
            # (b o q_1)_2 = -b_2 c_2^2 / 2 for a fixed arbitrary b
            if b2 * qq[1] != -b2 * half * cc[1] * cc[1]:
                bad_had += 1
    say("independent tableaux tested (sizes 3,4,5): %d" % total)
    say("HONESTY NOTE: this family is generated strictly lower triangular with")
    say("  c := A*1, so the three checks below hold IDENTICALLY on it -- they")
    say("  confirm the Proposition's algebra but have no discriminating power")
    say("  by themselves.  The discriminating power is supplied by the")
    say("  falsification arm in the next block, which shows the same three")
    say("  statements FAIL as soon as the hypothesis is dropped.")
    check(bad_c1 == 0,
          "c_1 = 0 on every internally consistent explicit tableau tested",
          "%d/%d failures" % (bad_c1, total))
    check(bad_q1 == 0,
          "(q_1)_2 = -c_2^2/2 on every such tableau (the Proposition's step)",
          "%d/%d failures" % (bad_q1, total))
    check(bad_had == 0,
          "(b o q_1)_2 = -b_2 c_2^2/2, so QO(m>=2) forces b_2 c_2^2 = 0",
          "%d/%d failures" % (bad_had, total))
    return total


def check_lemma_falsification():
    """FALSIFICATION ARM for D14: show the Proposition's hypothesis is what does
    the work, so the positive sweep above is not an empty exercise.

    One generator produces ALL 3x3 matrices with entries drawn from a fixed
    three-element rational set; c := A*1 throughout, so internal consistency is
    kept in both classes and the ONLY thing that varies is whether A is
    strictly lower triangular.  Inside that class the three lemma statements
    must hold with zero exceptions; outside it each of them must fail at least
    once.  A check that could not fail would show up here as a zero count.
    """
    vals = (F(0), F(1), F(-1, 2))
    b2 = F(1)                      # any nonzero b_2 exposes the third statement
    half = F(1, 2)
    n = 3
    slt_tot = slt_bad = 0
    oth_tot = 0
    oth_c1 = oth_q1 = oth_had = 0
    ex = {}
    for flat in _iproduct(vals, repeat=n * n):
        M = [list(flat[i * n:(i + 1) * n]) for i in range(n)]
        is_slt = all(M[i][j] == 0 for i in range(n) for j in range(i, n))
        cc = [sum(M[i], F(0)) for i in range(n)]
        acc = [sum((M[i][j] * cc[j] for j in range(n)), F(0)) for i in range(n)]
        qq = [acc[i] - half * cc[i] * cc[i] for i in range(n)]
        f_c1 = (cc[0] != 0)
        f_q1 = (qq[1] != -half * cc[1] * cc[1])
        f_had = (b2 * qq[1] != -b2 * half * cc[1] * cc[1])
        if is_slt:
            slt_tot += 1
            if f_c1 or f_q1 or f_had:
                slt_bad += 1
        else:
            oth_tot += 1
            oth_c1 += 1 if f_c1 else 0
            oth_q1 += 1 if f_q1 else 0
            oth_had += 1 if f_had else 0
            if f_q1 and "q1" not in ex:
                ex["q1"] = (M, cc, qq[1], -half * cc[1] * cc[1])
    say("generator produced %d matrices: %d strictly lower triangular, %d not"
        % (slt_tot + oth_tot, slt_tot, oth_tot))
    say("outside the hypothesis: c_1 != 0 on %d, (q_1)_2 != -c_2^2/2 on %d,"
        " (b o q_1)_2 mismatched on %d" % (oth_c1, oth_q1, oth_had))
    if "q1" in ex:
        M, cc, got, want = ex["q1"]
        say("witness A = [%s]" % "; ".join(vstr(r) for r in M))
        say("   -> c = %s, (q_1)_2 = %s but -c_2^2/2 = %s" % (vstr(cc), got, want))
    check(slt_tot > 0 and oth_tot > 0,
          "the falsification generator populates BOTH classes (no empty universe)",
          "%d in-hypothesis, %d out-of-hypothesis" % (slt_tot, oth_tot))
    check(slt_bad == 0,
          "inside the hypothesis the lemma never fails, on the same generator",
          "%d/%d failures" % (slt_bad, slt_tot))
    check(oth_c1 > 0 and oth_q1 > 0 and oth_had > 0,
          "DROPPING strict lower triangularity breaks all three lemma statements,"
          " so the positive sweep is discriminating and not a tautology",
          "counterexample counts c_1: %d, (q_1)_2: %d, (b o q_1)_2: %d"
          % (oth_c1, oth_q1, oth_had))
    return slt_tot + oth_tot


def derive_exact_order(viol, max_order):
    """Largest k <= max_order with zero violations at every order 1..k.

    This is read off the sweep, so the value of p fed to the integer argument
    is DERIVED and not the literal 5.
    """
    k = 0
    for order in range(1, max_order + 1):
        if viol[order]:
            break
        k = order
    return k


def check_quadrature_B(c, p):
    """D17: the quadrature conditions B(k), b_hat^T c^{o(k-1)} = 1/k, k = 1..p.

    B(p) is the first member of the Theorem 1.2 package.  It is evaluated here
    directly from b_hat and the DERIVED c, at the DERIVED order p, so that the
    GAPS block can name exactly which members of the package go unevaluated
    instead of asserting a list.  The section-10 sweep already imposes the same
    conditions through the bush trees ., [.], [..], ... by a different route
    (elementary weights Phi(t)); agreement of the two routes is why the sweep's
    'OK' lines and these checks must stand or fall together.

    The powers are built by repeated componentwise multiplication starting from
    the all-ones vector, so c^{o0} = 1 holds by construction and no 0^0
    convention is needed even though c_1 = 0.
    """
    verified = []
    cpow = [F(1)] * N                      # c^{o0}
    for k in range(1, p + 1):
        lhs = dot(BHAT, cpow)
        rhs = F(1, k)
        say("b_hat^T c^{o%d} = %-8s (B(%d) component, target 1/%d = %s)"
            % (k - 1, str(lhs), k, k, rhs))
        ok = check(lhs == rhs,
                   "B(%d) component b_hat^T c^{o%d} = 1/%d holds" % (k, k - 1, k),
                   "derived %s, target %s" % (lhs, rhs))
        if ok:
            verified.append(k)
        cpow = hadamard(cpow, c)
    return verified


def _b_gap_sentence(p, b_verified):
    """The B(p) clause of gap G6, DERIVED from what check_quadrature_B decided.

    B(p) is deliberately NOT in G6's list of unevaluated conditions, because the
    program does evaluate it.  This sentence is generated from the returned list
    of verified components so that it cannot drift away from the checks.
    """
    full = list(range(1, p + 1))
    if b_verified == full:
        body = ("B(p) is NOT among them: it IS evaluated, at the derived"
                " p = %d, both directly (b_hat^T c^{o(k-1)} = 1/k for"
                " k = 1..%d, all %d components verified above) and again by the"
                " bush trees of the order-1..%d sweep, and it HOLDS on this"
                " tableau.  A method may of course satisfy B(p) and still fail"
                " QO(m), which is exactly what happens here."
                % (p, p, len(full), p))
    else:
        missing = [k for k in full if k not in b_verified]
        body = ("B(p) at the derived p = %d is evaluated above and does NOT"
                " hold: the components k = %s fail, so the arithmetic here is"
                " inconsistent with the paper's order-%d claim and the whole"
                " run must be treated as failed."
                % (p, missing, p))
    return textwrap.fill(body, width=74,
                         initial_indent="     ", subsequent_indent="     ")


def print_gaps(p, b_verified):
    """THE GAP BLOCK: what this program does NOT decide.

    Everything above is exact arithmetic on the transcribed tableau.  The steps
    below sit between those facts and the paper's conclusion and are NOT
    machine-checkable from the tableau; they are quoted from the paper (which in
    turn quotes arXiv:2605.16995).  A reader who rejects any one of them should
    not accept the refutation, however many checks pass above.
    """
    print("\n" + "=" * 78)
    print("GAPS -- steps NOT decided by this program (read before trusting it)")
    print("=" * 78)
    gaps = [
        ("G1 SOURCE FIDELITY.  q_r = A c^{or} - c^{o(r+1)}/(r+1), the generator"
         " list\n     Q_2 = span{q_1, A q_0, q_0 o c, q_0}, the chain"
         " Q_1 <= Q_2 <= ..., the reading\n     of QO(m) as the Hadamard"
         " b o Q_m = {0}, and the admissibility constraints\n     m >= n-1,"
         " m+n+1 >= p are all TRANSCRIBED from the paper's quotation of\n"
         "     He-Huang.  No line of this program can detect a misquotation."
         "  If any of\n     them is misquoted, every check above still passes"
         " and the paper is wrong."),
        ("G2 Q_2 <= Q_m FOR m >= 2.  The witness q_1 is shown to lie in the"
         " printed\n     spanning set of Q_2 and to satisfy b o q_1 != 0."
         "  That it also lies in\n     Q_m for every m >= 2 is taken from the"
         " quoted chain; this program never\n     constructs Q_m for m >= 3,"
         " so non-monotonicity of the He-Huang spaces\n     would defeat the"
         " 'for every admissible (m,n)' step while leaving the\n     QO(2)"
         " computation intact."),
        ("G3 INTEGRALITY OF (m,n).  min m = 2 needs m to be an integer;"
         " section 12\n     shows the half-integer relaxation gives 3/2."
         "  The paper attributes\n     'nonnegative integers m,n' to"
         " Theorem 1.2, while the source is reported\n     to write only"
         " 'for some m, n'.  This program verifies the claim UNDER the\n"
         "     paper's stated hypothesis; it cannot verify the hypothesis."),
        ("G4 HADAMARD VERSUS SCALAR.  The refutation is specific to the"
         " componentwise\n     reading.  Section 11 shows b_hat^T q_1 = 0"
         " EXACTLY, so under a scalar\n     reading of QO(m) nothing here"
         " refutes anything.  The paper says so; a\n     reader who thinks"
         " He-Huang meant the scalar product should reject the\n     paper,"
         " not this program."),
        ("G5 ATTRIBUTION TO VERNER.  That this tableau is Table 4 of Verner"
         " (2014),\n     RK(8-6:5)a, is NOT checked and cannot be checked"
         " here.  It is also not\n     needed: the order-five property is"
         " derived from the tableau itself, so a\n     misattribution would"
         " be a citation defect, not a mathematical one."),
        ("G6 THE OTHER FOUR CONDITIONS.  DO(n), QD_weak(m,n), PR(n) and QR(m)"
         "\n     are never evaluated.  This is sound -- refuting one"
         " member of a\n     conjunction refutes the conjunction -- but it"
         " means the program says\n     nothing about whether those four"
         " members of the Theorem 1.2 package\n     hold.\n"
         + _b_gap_sentence(p, b_verified)),
        ("G7 NO PRIOR-ART CHECK.  Novelty of the counterexample is a literature"
         "\n     question and is outside the scope of any arithmetic."),
    ]
    for g in gaps:
        print("  " + g)
    print("=" * 78)


def main():
    print("=" * 78)
    print("verify.py -- Verner RK(8-6:5)a embedded formula vs He-Huang")
    print("                 Conjecture 2.1 (arXiv:2605.16995v3)")
    print("All arithmetic is exact (int / fractions.Fraction); no floating point,")
    print("hence no numeric tolerance and no error bound is needed.")
    print("=" * 78)

    print("\n-- 1. structure of the transcribed tableau -----------------------")
    check_tableau_structure()

    print("\n-- 2. internal consistency c = A*1 and sum(b_hat) = 1 ------------")
    c = check_internal_consistency()

    print("\n-- 3. q_0, and why the obstruction needs m >= 2 ------------------")
    q0, _ = check_q0(c)

    print("\n-- 4. Ac and q_1 = Ac - c^{o2}/2 ---------------------------------")
    q1, ac, csq = check_q1(c)

    print("\n-- 5. the violation b_hat o q_1 != 0 -----------------------------")
    check_hadamard_violation(q1)

    print("\n-- 6. QO(2) over the printed generators of Q_2 -------------------")
    check_QO2_generators(q0, q1, c)

    print("\n-- 7. the Proposition's forced quantity b_2 c_2^2 ----------------")
    check_proposition_quantity(c, q1)

    print("\n-- 8. chain conditions b_hat^T A^k c -----------------------------")
    chain6_residual = check_chain_conditions(c)

    print("\n-- 9. rooted trees to order six, generated from scratch ----------")
    max_order = 6
    trees_by_order = gen_trees(max_order)
    check_tree_counts(trees_by_order, max_order)

    print("\n-- 10. every order condition through order six -------------------")
    viol = sweep_order_conditions(trees_by_order, max_order)
    check_exact_order(trees_by_order, viol, max_order, chain6_residual)
    p_derived = derive_exact_order(viol, max_order)
    say("exact order DERIVED from the sweep (not read from a constant): p = %d"
        % p_derived)

    print("\n-- 10b. the quadrature conditions B(k) at the DERIVED p ----------")
    b_verified = check_quadrature_B(c, p_derived)

    print("\n-- 11. Scope: scalar orthogonality is not violated ---------------")
    check_scope(c, q1, ac, csq)

    print("\n-- 12. admissible (m,n) at the DERIVED order p -------------------")
    check_admissible_m(p_derived)

    print("\n-- 13. the Proposition's lemma on independent tableaux -----------")
    check_proposition_lemma()

    print("\n-- 14. falsification arm: is that lemma sweep discriminating? ----")
    check_lemma_falsification()

    print("\n" + "=" * 78)
    print("SUMMARY OF THE REFUTATION (all quantities above were derived here)")
    print("  the method is internally consistent and explicit;")
    print("  it satisfies all 17 order conditions of orders 1..5 and fails 8 of")
    print("  the 20 at order 6, so its order is exactly 5 >= 5;")
    print("  q_1 in Q_2 <= Q_m has b_hat o q_1 != 0, so QO(m) fails, and every")
    print("  admissible (m,n) at p = 5 has m >= 2, so QO(m) fails for all of them;")
    print("  a method of order 5 therefore violates a member of the Theorem 1.2")
    print("  package, refuting the necessity asserted by Conjecture 2.1")
    print("  -- CONDITIONAL on the seven items in the GAPS block below, which")
    print("  are quoted from the paper and are not decided by any check here.")
    print("=" * 78)

    print_gaps(p_derived, b_verified)

    n = len(_RESULTS)
    failed = [lab for ok, lab in _RESULTS if not ok]
    for lab in failed:
        print("FAILED CHECK: " + lab)
    if failed:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())

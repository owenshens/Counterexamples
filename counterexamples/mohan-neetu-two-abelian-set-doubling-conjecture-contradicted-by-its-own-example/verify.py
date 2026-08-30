#!/usr/bin/env python3
"""verify.py -- re-derives every quantity claimed in paper.tex / paper.pdf of

    "A Four-Element Counterexample to Conjecture~6.1 of Mohan and Neetu on Two Disjoint Abelian Setsnt Abelian Sets Is Contradicted
     by the Example Printed Below It"

from the objects PRINTED IN THAT PAPER and from nothing else.

Python 3.9+, standard library only.  No third-party package, no external data file, no
network.  All arithmetic is on pairs of Python integers; no float appears anywhere, and the
sign (-1)^n of the group law is taken as a PARITY LOOKUP rather than a power, so a negative
exponent cannot silently produce a float.

Every input is transcribed literally from the paper:

  * the group law and the inverse, equation (1);
  * the family A_m, B_m, S_m of equation (2);
  * the closed form for S_m^2 in Theorem 3(2);
  * the sixteen products of S_2 displayed after the proof of Theorem 3;
  * the set S = {(0,2),(0,3),(1,3),(0,4)} of Theorem 4 with its A, B, its sixteen products
    and its nine-element S^2;
  * the authors' own published integers quoted in the Remarks: |S^2| = 3|S|-3 for the
    three-set Example 6.1 at |S| = 4, 8, 12, 16, 20, and |S^2| = 2k-1 for Example 6.2.

The program prints one `PASS <name> <detail>` line per check, closes with
`VERDICT: ALL <n> CHECKS PASS`, and exits 0 if and only if every check passed.
"""

import sys

# ---------------------------------------------------------------------------
# THE GROUP.  Equation (1) of the paper, which is [MN, section 5] at q = -1:
#     (r,n)(s,m) = (r + (-1)^n s, n + m),      (r,n)^-1 = (-(-1)^n r, -n).
# ---------------------------------------------------------------------------
E = (0, 0)


def sgn(n):
    """(-1)**n as a parity LOOKUP.  Never a power: n may be negative."""
    return 1 if n % 2 == 0 else -1


def mul(u, v):
    (r, n), (s, m) = u, v
    return (r + sgn(n) * s, n + m)


def inv(u):
    r, n = u
    return (-sgn(n) * r, -n)


def power(x, k):
    """x^k for k >= 1, by repeated multiplication in the group law."""
    y = x
    for _ in range(k - 1):
        y = mul(y, x)
    return y


def has_finite_order_upto(x, kmax):
    """TRUE iff x^k = e for some 1 <= k <= kmax.  Incremental, so no recomputation."""
    y = x
    for _ in range(kmax):
        if y == E:
            return True
        y = mul(y, x)
    return False


def prodset(S):
    """S^2 over ORDERED pairs, so s^2 in S^2 for every s in S."""
    return frozenset(mul(a, b) for a in S for b in S)


def pairwise_commuting(X):
    return all(mul(a, b) == mul(b, a) for a in X for b in X)


def nonabelian_witness(X):
    """A non-commuting pair of X, or None.  `pairwise_commuting` is the abelian-set test
    justified by Lemma 2 of the paper; a failure of it is a proof that <X> is nonabelian."""
    for a in X:
        for b in X:
            if mul(a, b) != mul(b, a):
                return (a, b)
    return None


# ---------------------------------------------------------------------------
# THE ONE DECIDER, used in BOTH polarities below.  A predicate that can only say YES is
# worthless, so the same function is applied to the witnesses and to the controls.
# ---------------------------------------------------------------------------
def refutes_conjecture_6_1(S, A, B):
    """(bool, reason).  TRUE iff (S, A, B) meets every clause of the antecedent of
    Conjecture 6.1 -- S = A u B with A, B disjoint abelian sets and |S^2| <= 3|S|-3 -- and
    fails its conclusion, i.e. <S> is nonabelian.  The host group K is torsion-free
    (Lemma 1), which is the remaining hypothesis and does not depend on S."""
    S, A, B = frozenset(S), frozenset(A), frozenset(B)
    if A | B != S:
        return False, 'A u B is not S'
    if A & B:
        return False, 'A and B are not disjoint'
    if not pairwise_commuting(A):
        return False, 'A is not an abelian set'
    if not pairwise_commuting(B):
        return False, 'B is not an abelian set'
    d = len(prodset(S))
    if d > 3 * len(S) - 3:
        return False, '|S^2| = %d exceeds 3|S|-3 = %d' % (d, 3 * len(S) - 3)
    if pairwise_commuting(S):
        return False, '<S> is abelian, so the conclusion holds'
    return True, '|S| = %d, |S^2| = %d <= 3|S|-3 = %d, <S> nonabelian' % (len(S), d, 3 * len(S) - 3)


# ---------------------------------------------------------------------------
# THE CHECK HARNESS
# ---------------------------------------------------------------------------
_N = [0]
_FAILED = []


def check(name, cond, detail=''):
    _N[0] += 1
    if cond:
        print('PASS %s %s' % (name, detail))
    else:
        print('FAIL %s %s' % (name, detail))
        _FAILED.append(name)


def head(title):
    print()
    print('=== %s ===' % title)


# ---------------------------------------------------------------------------
# THE OBJECTS, transcribed from the paper
# ---------------------------------------------------------------------------
def A_m(m):
    return [(-1, 2 * i - 1) for i in range(1, m + 1)]


def B_m(m):
    return [(1, 2 * i - 1) for i in range(1, m + 1)]


def S_m(m):
    return A_m(m) + B_m(m)


def C_m(m):
    """The third part of the authors' Example 6.1, restored."""
    return [(c, 2 * i) for i in range(1, m + 1) for c in (-1, 1)]


def closed_form_S_m_sq(m):
    """Theorem 3(2): S_m^2 = {(d,2j) : d in {-2,0,2}, 1 <= j <= 2m-1}."""
    return frozenset((d, 2 * j) for d in (-2, 0, 2) for j in range(1, 2 * m))


# Theorem 4's object, from the %-commented example at source line 724.
T4_S = [(0, 2), (0, 3), (1, 3), (0, 4)]
T4_A = [(0, 2), (0, 3)]
T4_B = [(1, 3), (0, 4)]
T4_SQ_PRINTED = frozenset([(0, 4), (0, 5), (1, 5), (-1, 6), (0, 6), (1, 6), (0, 7), (1, 7), (0, 8)])

# The sixteen products of T4_S, exactly as the paper's table prints them.
T4_TABLE_PRINTED = {
    ((0, 2), (0, 2)): (0, 4), ((0, 2), (0, 3)): (0, 5), ((0, 2), (1, 3)): (1, 5),
    ((0, 2), (0, 4)): (0, 6), ((0, 3), (0, 2)): (0, 5), ((0, 3), (0, 3)): (0, 6),
    ((0, 3), (1, 3)): (-1, 6), ((0, 3), (0, 4)): (0, 7), ((1, 3), (0, 2)): (1, 5),
    ((1, 3), (0, 3)): (1, 6), ((1, 3), (1, 3)): (0, 6), ((1, 3), (0, 4)): (1, 7),
    ((0, 4), (0, 2)): (0, 6), ((0, 4), (0, 3)): (0, 7), ((0, 4), (1, 3)): (1, 7),
    ((0, 4), (0, 4)): (0, 8),
}

# The sixteen products of S_2, exactly as the paper's table prints them.
S2_TABLE_PRINTED = {
    ((-1, 1), (-1, 1)): (0, 2), ((-1, 1), (1, 1)): (-2, 2),
    ((1, 1), (-1, 1)): (2, 2), ((1, 1), (1, 1)): (0, 2),
    ((-1, 1), (-1, 3)): (0, 4), ((-1, 1), (1, 3)): (-2, 4),
    ((1, 1), (-1, 3)): (2, 4), ((1, 1), (1, 3)): (0, 4),
    ((-1, 3), (-1, 1)): (0, 4), ((-1, 3), (1, 1)): (-2, 4),
    ((1, 3), (-1, 1)): (2, 4), ((1, 3), (1, 1)): (0, 4),
    ((-1, 3), (-1, 3)): (0, 6), ((-1, 3), (1, 3)): (-2, 6),
    ((1, 3), (-1, 3)): (2, 6), ((1, 3), (1, 3)): (0, 6),
}
S2_SQ_PRINTED = frozenset([(-2, 2), (0, 2), (2, 2), (-2, 4), (0, 4), (2, 4), (-2, 6), (0, 6), (2, 6)])

# The authors' own published integers, Remark "the authors' three-set claim is correct".
EX61_PRINTED = {1: (4, 9), 2: (8, 21), 3: (12, 33), 4: (16, 45), 5: (20, 57)}

MMAX = 40           # the paper says the family is machine-checked for m <= 40
BOX = range(-4, 5)  # coordinate box for the group-axiom and torsion checks
TORSION_EXP = 200   # bounded torsion search


def main():
    # -----------------------------------------------------------------
    head('1. THE GROUP LAW, equation (1), and Lemma 1')
    box = [(r, n) for r in BOX for n in BOX]

    check('identity',
          all(mul(E, x) == x and mul(x, E) == x for x in box),
          '(0,0) is a two-sided identity on the %d-element box r,n in [%d,%d]'
          % (len(box), BOX[0], BOX[-1]))

    small = [(r, n) for r in range(-2, 3) for n in range(-2, 3)]
    check('associativity',
          all(mul(mul(a, b), c) == mul(a, mul(b, c)) for a in small for b in small for c in small),
          'over all %d ordered triples from the 25-element box' % (len(small) ** 3))

    check('inverse-formula',
          all(mul(x, inv(x)) == E and mul(inv(x), x) == E for x in box),
          '(r,n)^-1 = (-(-1)^n r, -n) is two-sided on the %d-element box' % len(box))

    check('parity-not-power',
          sgn(-3) == -1 and sgn(-4) == 1 and isinstance(sgn(-3), int),
          '(-1)^n is a parity lookup, so negative n cannot produce a float')

    ab = nonabelian_witness([(-1, 1), (1, 1)])
    check('lemma1-nonabelian',
          mul((-1, 1), (1, 1)) == (-2, 2) and mul((1, 1), (-1, 1)) == (2, 2) and ab is not None,
          '(-1,1)(1,1) = (-2,2) != (2,2) = (1,1)(-1,1), as Lemma 1 prints')

    bad = [x for x in box if x != E and has_finite_order_upto(x, TORSION_EXP)]
    check('lemma1-torsion-box',
          not bad,
          'no non-identity element of the %d-element box has order <= %d (BOUNDED; the '
          'unbounded statement is the hand proof of Lemma 1)' % (len(box), TORSION_EXP))

    # The authors' own two printed identities in Example 6.2 -- these are what pin q = -1.
    rng = range(-6, 7)
    check('ex62-square-identity',
          all(mul((i, 1), (i, 1)) == (0, 2) for i in rng),
          "(i,1)^2 = (0,2) for i in [-6,6], the identity Example 6.2 of [MN] prints")
    check('ex62-product-identity',
          all(mul((i, 1), (j, 1)) == (i - j, 2) for i in rng for j in rng),
          "(i,1)(j,1) = (i-j,2) for i,j in [-6,6]; this forces q = -1")

    # -----------------------------------------------------------------
    head('2. THEOREM 3 -- the authors\' Example 6.1 with its third part deleted')
    for m in range(1, MMAX + 1):
        A, B, S = A_m(m), B_m(m), S_m(m)
        Sset = frozenset(S)
        sq = prodset(S)
        cf = closed_form_S_m_sq(m)
        check('thm3-structure-m%02d' % m,
              len(Sset) == 2 * m and frozenset(A) | frozenset(B) == Sset
              and not (frozenset(A) & frozenset(B))
              and pairwise_commuting(A) and pairwise_commuting(B),
              '|S_m| = %d, A_m and B_m disjoint, both pairwise commuting (abelian sets by '
              'Lemma 2)' % len(Sset))
        check('thm3-productset-m%02d' % m,
              sq == cf,
              'S_m^2 equals the closed form {(d,2j): d in {-2,0,2}, 1<=j<=%d} AS A SET '
              '(%d elements), not merely in cardinality' % (2 * m - 1, len(sq)))
        check('thm3-count-m%02d' % m,
              len(sq) == 6 * m - 3 == 3 * len(Sset) - 3,
              '|S_m^2| = %d = 6m-3 = 3|S_m|-3, with EQUALITY, so |S_m^2| <= 3|S_m|-3'
              % len(sq))
        w = nonabelian_witness(S)
        ok, why = refutes_conjecture_6_1(S, A, B)
        check('thm3-refutes-m%02d' % m,
              ok and w is not None,
              'decider says YES: %s; non-commuting pair %s' % (why, w))

    # the m = 2 table, printed in full in the paper
    got = {(a, b): mul(a, b) for a in S_m(2) for b in S_m(2)}
    check('thm3-m2-table',
          got == S2_TABLE_PRINTED,
          'all 16 ordered products of S_2 agree with the table printed in the paper')
    check('thm3-m2-productset',
          prodset(S_m(2)) == S2_SQ_PRINTED and len(S2_SQ_PRINTED) == 9,
          'S_2^2 is the nine elements printed, and 9 = 3*4-3 = 3|S_2|-3')

    # -----------------------------------------------------------------
    head('3. THEOREM 4 -- the %-commented example at source line 724')
    S, A, B = T4_S, T4_A, T4_B
    Sset, sq = frozenset(S), prodset(S)
    check('thm4-size', len(Sset) == 4, '|S| = 4 for S = {(0,2),(0,3),(1,3),(0,4)}')
    check('thm4-partition',
          frozenset(A) | frozenset(B) == Sset and not (frozenset(A) & frozenset(B)),
          'A = {(0,2),(0,3)} and B = {(1,3),(0,4)} are disjoint with union S')
    check('thm4-A-abelian', pairwise_commuting(A),
          'A is pairwise commuting, hence an abelian set by Lemma 2')
    check('thm4-B-abelian', pairwise_commuting(B),
          'B is pairwise commuting: (1,3)(0,4) = %s = (0,4)(1,3)' % (mul((1, 3), (0, 4)),))
    got = {(a, b): mul(a, b) for a in S for b in S}
    check('thm4-table', got == T4_TABLE_PRINTED,
          'all 16 ordered products of S agree with the table printed in the paper')
    check('thm4-productset', sq == T4_SQ_PRINTED,
          'S^2 is exactly the nine elements printed in Theorem 4')
    check('thm4-count', len(sq) == 9 == 3 * len(Sset) - 3,
          "|S^2| = 9 = 3|S|-3 with EQUALITY -- the authors' own printed integer 9")
    check('thm4-nonabelian',
          mul((0, 3), (1, 3)) == (-1, 6) and mul((1, 3), (0, 3)) == (1, 6),
          '(0,3)(1,3) = (-1,6) != (1,6) = (1,3)(0,3), so <S> is nonabelian')
    ok, why = refutes_conjecture_6_1(S, A, B)
    check('thm4-refutes', ok, 'decider says YES: %s' % why)

    # -----------------------------------------------------------------
    head("4. POSITIVE CONTROLS -- the authors' own published integers must reproduce")
    for m in sorted(EX61_PRINTED):
        A, B, C = A_m(m), B_m(m), C_m(m)
        S = A + B + C
        Sset, sq = frozenset(S), prodset(S)
        want_n, want_d = EX61_PRINTED[m]
        check('ex61-threeset-m%d' % m,
              len(Sset) == want_n and len(sq) == want_d == 3 * len(Sset) - 3
              and pairwise_commuting(A) and pairwise_commuting(B) and pairwise_commuting(C)
              and nonabelian_witness(S) is not None,
              "|S| = %d and |S^2| = %d = 3|S|-3, all three parts abelian, <S> nonabelian "
              "-- THEIR integers, not ours" % (len(Sset), len(sq)))
    for k in range(2, 9):
        S = [(i, 1) for i in range(1, k + 1)]
        sq = prodset(S)
        check('ex62-doubling-k%d' % k,
              len(sq) == 2 * k - 1 and sq == frozenset((d, 2) for d in range(1 - k, k)),
              "|S^2| = %d = 2k-1 and S^2 = {(1-k,2),...,(k-1,2)}, both as Example 6.2 "
              "prints them" % len(sq))

    # -----------------------------------------------------------------
    head('5. NEGATIVE CONTROLS -- the decider must be able to say NO')
    ab4 = [(0, 0), (0, 1), (0, 2), (0, 3)]
    sq = prodset(ab4)
    ok, why = refutes_conjecture_6_1(ab4, ab4[:2], ab4[2:])
    check('control-abelian-inside-antecedent',
          (not ok) and len(sq) == 7 <= 3 * len(ab4) - 3 and pairwise_commuting(ab4),
          'S = {(0,0),(0,1),(0,2),(0,3)} has |S^2| = 7 <= 9 = 3|S|-3 AND abelian span, so '
          'the decider correctly returns NO (%s)' % why)

    ok, why = refutes_conjecture_6_1([(0, 1), (1, 1)], [(0, 1)], [(1, 1)])
    check('control-antecedent-at-the-bound',
          ok,
          'and the decider is not constant-NO on small sets: {(0,1),(1,1)} split into '
          'singletons does refute (%s)' % why)

    check('anticontrol-not-a-commuting-pair',
          not pairwise_commuting([(0, 1), (1, 1)]),
          '{(0,1),(1,1)} is NOT a commuting pair -- (0,1)(1,1) = %s, (1,1)(0,1) = %s -- so '
          'it may not serve as an abelian PART; a checker that swallowed it would be broken'
          % (mul((0, 1), (1, 1)), mul((1, 1), (0, 1))))

    bps = [(0, 0), (0, 1), (1, 0)]
    sq = prodset(bps)
    ok, why = refutes_conjecture_6_1(bps, bps[:2], bps[2:])
    check('control-outside-antecedent',
          (not ok) and len(sq) == 7 == 3 * len(bps) - 2,
          'the Boroczky-Palfy-Serra set {e,(0,1),(1,0)} has |A^2| = 7 = 3|A|-2, one ABOVE '
          'the bound, so it lies outside the antecedent and the decider returns NO (%s)'
          % why)

    check('control-nonabelian-part-rejected',
          not refutes_conjecture_6_1([(0, 1), (1, 1), (0, 0)],
                                     [(0, 1), (1, 1)], [(0, 0)])[0],
          'a candidate whose first part is the non-commuting {(0,1),(1,1)} is rejected, '
          'even though its span is nonabelian')

    # -----------------------------------------------------------------
    head('6. THE SHARPNESS REMARK -- A_m^2 = B_m^2, so one hypothesis of Theorem 2.10 fails')
    for m in range(1, MMAX + 1):
        A, B = A_m(m), B_m(m)
        sqA, sqB = prodset(A), prodset(B)
        want = frozenset((0, 2 * j) for j in range(1, 2 * m))
        check('sharp-squares-m%02d' % m,
              sqA == sqB == want and len(sqA) == 2 * m - 1,
              'A_m^2 = B_m^2 = {(0,2),...,(0,%d)} with |A_m^2| = %d' % (4 * m - 2, len(sqA)))
    for m in range(4, 11):
        A, B = A_m(m), B_m(m)
        sqA, sqB = prodset(A), prodset(B)
        check('sharp-thm210-m%02d' % m,
              len(A) >= 4 and len(B) >= 3 and sqA == sqB and len(sqA) > 1,
              '|A| = %d >= 4 and |B| = %d >= 3 hold, while max(A^2) <= min(B^2) cannot, '
              'since A^2 = B^2 has %d > 1 elements -- so it fails for EVERY right-order'
              % (len(A), len(B), len(sqA)))
    check('sharp-thm27-untouched',
          all(len(prodset(S_m(m))) > 3 * len(frozenset(S_m(m))) - 4 for m in range(1, MMAX + 1)),
          'for every m <= %d, |S_m^2| = 3|S_m|-3 > 3|S_m|-4, so Theorem 2.7 of [MN] is not '
          'contradicted -- the family misses its hypothesis by exactly one' % MMAX)

    # -----------------------------------------------------------------
    head('7. WHAT THIS PROGRAM DOES NOT COVER')
    for line in (
        'NOT RE-RUN: torsion-freeness of K is confirmed only on the %d-element coordinate '
        'box r,n in [%d,%d] and only for exponents 1..%d; the unbounded statement is the '
        'hand proof of Lemma 1.' % (len(box), BOX[0], BOX[-1], TORSION_EXP),
        'NOT RE-RUN: the family of equation (2) is machine-checked for 1 <= m <= %d only; '
        'the statement for all m is Theorem 3, proved by hand.' % MMAX,
        'NOT RE-RUN: no search for a witness of smaller cardinality (in particular nothing '
        'is decided about |S| = 3), for a smaller doubling, or in any host group other '
        'than K = BS(1,-1).',
        'NOT RE-RUN: this program does not fetch the e-print, so the byte locator, the '
        'line numbers and the SHA-256 of section 1 of the paper are NOT re-checked here; '
        'they were read off the source file directly.',
        'NOT RE-RUN: nothing here bears on Conjecture 6.2 or Conjecture 6.3 of [MN], nor on '
        'whether any repaired form of Conjecture 6.1 is true.',
        'NOT RE-RUN: the provenance question of whether the line-724 object is inherited '
        'from the closed-access companion paper is not a computation and is not addressed.',
    ):
        print(line)

    print()
    if _FAILED:
        print('CHECKS FAILED: %d of %d -- %s' % (len(_FAILED), _N[0], ', '.join(_FAILED)))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _N[0])
    return 0


if __name__ == '__main__':
    sys.exit(main())

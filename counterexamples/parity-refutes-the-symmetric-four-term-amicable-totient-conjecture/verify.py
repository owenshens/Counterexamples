#!/usr/bin/env python3
"""Verification program for

    "Parity Refutes the Printed Four-Term Amicable-Totient Equations"

Everything this program consumes is TRANSCRIBED FROM THE PAPER, in the literal
blocks PAPER_WITNESS, PAPER_CELL_TABLE and PAPER_CENSUS below.  No external data
file, no network, no third-party module: Python 3.9+ standard library only.  All
arithmetic is exact integer arithmetic; there is no floating point anywhere and
no decision is taken on a float.

The program re-derives, from the primes printed in the paper:

  A. the witness object of Theorem 3 -- factorisations, primality, gcd, sigma,
     phi, S, x, y, and the paper's own hand-multiplication chains;
  B. the six achievable values of each bracket over each member;
  C. all SIXTEEN admissible equations of (5) over all 576 labellings;
  D. Proposition 4 (the corrected identity) and every step of its proof;
  E. Lemma 1 over a family of 4845 squarefree odd integers that contains no
     amicable pair at all -- the lemma is re-tested as an OUTPUT, not assumed;
  F. consistency with the three identities the target paper PROVES;
  G. the seven further cell members tabulated in Section 5, from scratch;
  H. the internal arithmetic of the census figures quoted in Section 5.

It prints one `PASS <name> [detail]` line per check, a SCOPE block naming what it
does NOT cover, and closes with `VERDICT: ALL <n> CHECKS PASS`, exiting 0 if and
only if every check passed.
"""

import sys
from itertools import combinations, permutations
from math import gcd, prod

# ---------------------------------------------------------------------------
# THE PAPER'S OWN NUMBERS, TRANSCRIBED
# ---------------------------------------------------------------------------
# Theorem 3 and its proof, verbatim as numbers.
PAPER_WITNESS = """
n                       2
A                       32642324
B                       35095276
primes_A                11 13 149 383
primes_B                17 47 79 139
gcd                     4
chain_A                 143 57067 8160581 32642324
chain_B                 799 10981 8773819 35095276
sigma_A_factors         7 12 14 150 384
sigma_B_factors         7 18 48 80 140
sigma                   67737600
phi_A_factors           2 10 12 148 382
phi_B_factors           2 16 46 78 138
phi_A                   13568640
phi_B                   15844608
S                       29413248
A_minus_S               3229076
B_minus_S               5682028
x                       807269
y                       1420507
Sigma_A                 556
Sigma_B                 282
brackets_A              133143 654023 684443 767391 798279 1369751
brackets_B              185163 256331 301451 581591 628511 703583
spot_A_low              143 532 57067 133143
spot_A_high             57067 24 143 1369751
spot_B_high             10981 64 799 703583
required_offset_free    807551 807825 1420789 1421063
gap_to_798279           9272
overshoot_over_703583   103968
E1                      556
E2                      69978
E3                      1445684
F1                      282
F2                      25732
F3                      876966
"""

# Section 5, the eight tabulated members of the cell.  Columns:  n | primes of A
# | primes of B | A.
PAPER_CELL_TABLE = """
2   11 13 149 383            17 47 79 139               32642324
3   19 71 97 503             41 83 89 223               526552472
4   59 107 223 271           101 127 167 179            6104216464
5   113 191 587 719          167 251 359 607            291493436768
6   197 499 719 2311         449 479 577 1319           10453833569728
7   487 1307 1433 2339       653 701 1951 2389          273080651722624
8   1187 2239 2549 3527      1439 1847 2099 4283        6116740737238784
5   47 9437 59929 117503     103 2111 4861 2987279      99946542990038176
"""

# Section 5, the census figures.  n = 2..8.
PAPER_CENSUS = """
bfile_terms             415523
cell_members            14484
labellings_per_pair     576
labellings_total        8342784
n_distribution          111 1368 4085 4413 3077 1276 154
type_4_2_members        1073
"""

# The paper's own (2,2) example, quoted from the target paper's Theorem 2.1.
PAPER_22_EXAMPLE = """
A                       2620
B                       2924
primes_A                5 131
primes_B                17 43
sigma                   5544
S                       2384
x                       59
y                       135
"""


def parse_block(block):
    """key -> list of ints, for the paper blocks above."""
    out = {}
    for line in block.strip().splitlines():
        parts = line.split()
        out[parts[0]] = [int(v) for v in parts[1:]]
    return out


def parse_table(block):
    """-> [(n, primes_A, primes_B, A)]"""
    rows = []
    for line in block.strip().splitlines():
        v = [int(t) for t in line.split()]
        assert len(v) == 10, line
        rows.append((v[0], v[1:5], v[5:9], v[9]))
    return rows


# ---------------------------------------------------------------------------
# CHECK HARNESS
# ---------------------------------------------------------------------------
_STATE = {'pass': 0, 'fail': 0}


def chk(name, cond, detail=''):
    if cond:
        _STATE['pass'] += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _STATE['fail'] += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


# ---------------------------------------------------------------------------
# EXACT ARITHMETIC PRIMITIVES (nothing imported does any of this)
# ---------------------------------------------------------------------------
def is_prime(m):
    """Deterministic trial division. The largest integer tested here is 2987279,
    whose square root is under 1729, so this is both exact and instant."""
    if m < 2:
        return False
    if m % 2 == 0:
        return m == 2
    d = 3
    while d * d <= m:
        if m % d == 0:
            return False
        d += 2
    return True


def sigma_from(n, odd_primes):
    """sigma(2^n * p1 p2 p3 p4) for DISTINCT odd primes, from the factorisation."""
    return (2 ** (n + 1) - 1) * prod(p + 1 for p in odd_primes)


def phi_from(n, odd_primes):
    """phi(2^n * p1 p2 p3 p4) for DISTINCT odd primes, n >= 1."""
    return 2 ** (n - 1) * prod(p - 1 for p in odd_primes)


def esym(P):
    """[e_0, e_1, ..., e_k] of a list of integers."""
    return [sum(prod(c) for c in combinations(P, k)) if k else 1 for k in range(len(P) + 1)]


def bracket(p, q, r, s):
    """The paper's printed bracket, [pq(r+s)+rs]."""
    return p * q * (r + s) + r * s


def xy_of(n, PA, PB):
    """(x, y, S, A, B) from the two prime lists and n, by the paper's definitions."""
    A = 2 ** n * prod(PA)
    B = 2 ** n * prod(PB)
    S = phi_from(n, PA) + phi_from(n, PB)
    assert (A - S) % 2 ** n == 0 and (B - S) % 2 ** n == 0
    return (A - S) // 2 ** n, (B - S) // 2 ** n, S, A, B


def sixteen_equations(x, y):
    """The sixteen (LHS, bracket-side, sum-side) selectors of (5).

    LHS in {x-1, y-1, x, y}; bracket in {B_A, B_B}; sum in {Sigma_A, Sigma_B}.
    Returned as (label, lhs_name, L, which_bracket, which_sum), the two `which_`
    fields in {'A','B'}. ⛔ The left-hand side is identified by NAME, never by its
    value: two of the four values could coincide on some pair and a value-based
    filter would then silently merge two different equations.
    """
    eqs = []
    for lname, L in (('x-1', x - 1), ('y-1', y - 1), ('x', x), ('y', y)):
        for bname in ('A', 'B'):
            for sname in ('A', 'B'):
                eqs.append(('%s = B_%s - Sigma_%s' % (lname, bname, sname),
                            lname, L, bname, sname))
    return eqs


def sweep_576(n, PA, PB):
    """Every one of the 4! x 4! = 576 labellings against every one of the sixteen
    admissible equations, with NO pruning: 9216 integer comparisons per pair.

    -> (satisfied_per_equation dict, printed_x_count, printed_y_count,
        printed_both_count, exchanged_count, corrected_identity_count)
    """
    x, y, _S, _A, _B = xy_of(n, PA, PB)
    eqs = sixteen_equations(x, y)
    hits = dict((lab, 0) for lab, _, _, _, _ in eqs)
    cx = cy = cboth = cswap = ctrue = 0
    E = esym(PA)
    F = esym(PB)
    true_x = F[3] + F[1] - E[2] - 1
    true_y = E[3] + E[1] - F[2] - 1
    for pa in permutations(PA):
        BA = bracket(*pa)
        SA = sum(pa)
        for pb in permutations(PB):
            BB = bracket(*pb)
            SB = sum(pb)
            B_ = {'A': BA, 'B': BB}
            S_ = {'A': SA, 'B': SB}
            for lab, _ln, L, bn, sn in eqs:
                if L == B_[bn] - S_[sn]:
                    hits[lab] += 1
            okx = (x - 1 == BA - SB)          # the printed x-equation (3.3)
            oky = (y - 1 == BB - SA)          # the printed y-equation (3.2)
            cx += okx
            cy += oky
            cboth += (okx and oky)
            cswap += ((x - 1 == BB - SB) and (y - 1 == BA - SA))
            ctrue += (x == true_x and y == true_y)
    return hits, cx, cy, cboth, cswap, ctrue


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    W = parse_block(PAPER_WITNESS)
    TAB = parse_table(PAPER_CELL_TABLE)
    CEN = parse_block(PAPER_CENSUS)
    EX22 = parse_block(PAPER_22_EXAMPLE)

    n = W['n'][0]
    PA = W['primes_A']
    PB = W['primes_B']
    Ap = W['A'][0]
    Bp = W['B'][0]

    # =====================================================================
    print()
    print('=== A. THE WITNESS OBJECT AS PRINTED IN THEOREM 3 ===')
    # =====================================================================
    A = 2 ** n * prod(PA)
    B = 2 ** n * prod(PB)
    chk('A1-factorisation-of-A', A == Ap,
        '2^%d * %s = %d, paper prints %d' % (n, '*'.join(map(str, PA)), A, Ap))
    chk('A2-factorisation-of-B', B == Bp,
        '2^%d * %s = %d, paper prints %d' % (n, '*'.join(map(str, PB)), B, Bp))
    ch = W['chain_A']
    chk('A3-hand-multiplication-chain-for-A',
        [PA[0] * PA[1], PA[2] * PA[3], PA[0] * PA[1] * PA[2] * PA[3], A] == ch,
        'paper: %s' % ' -> '.join(map(str, ch)))
    ch = W['chain_B']
    chk('A4-hand-multiplication-chain-for-B',
        [PB[0] * PB[1], PB[2] * PB[3], PB[0] * PB[1] * PB[2] * PB[3], B] == ch,
        'paper: %s' % ' -> '.join(map(str, ch)))
    chk('A5-all-eight-odd-parts-are-prime', all(is_prime(p) for p in PA + PB),
        'trial division to the square root, largest %d' % max(PA + PB))
    chk('A6-eight-primes-distinct-so-type-is-(4,4)',
        len(set(PA + PB)) == 8 and len(PA) == len(PB) == 4,
        'both odd parts squarefree with exactly four prime factors')
    chk('A7-gcd-is-exactly-2^n', gcd(A, B) == 2 ** n == W['gcd'][0],
        'gcd(A,B) = %d = 2^%d, paper prints %d' % (gcd(A, B), n, W['gcd'][0]))
    sA = sigma_from(n, PA)
    sB = sigma_from(n, PB)
    chk('A8-sigma(A)-equals-the-printed-product',
        prod(W['sigma_A_factors']) == sA == W['sigma'][0],
        '%s = %d' % ('*'.join(map(str, W['sigma_A_factors'])), sA))
    chk('A9-sigma(B)-equals-the-printed-product',
        prod(W['sigma_B_factors']) == sB == W['sigma'][0],
        '%s = %d' % ('*'.join(map(str, W['sigma_B_factors'])), sB))
    chk('A10-the-pair-is-amicable', sA == sB == A + B and A != B,
        'sigma(A) = sigma(B) = A+B = %d' % (A + B))
    phA = phi_from(n, PA)
    phB = phi_from(n, PB)
    chk('A11-phi(A)-equals-the-printed-product',
        prod(W['phi_A_factors']) == phA == W['phi_A'][0],
        '%s = %d' % ('*'.join(map(str, W['phi_A_factors'])), phA))
    chk('A12-phi(B)-equals-the-printed-product',
        prod(W['phi_B_factors']) == phB == W['phi_B'][0],
        '%s = %d' % ('*'.join(map(str, W['phi_B_factors'])), phB))
    S = phA + phB
    x, y, S2, _, _ = xy_of(n, PA, PB)
    chk('A13-S-and-the-two-deficits',
        S == S2 == W['S'][0] and A - S == W['A_minus_S'][0] == 2 ** n * W['x'][0]
        and B - S == W['B_minus_S'][0] == 2 ** n * W['y'][0],
        'S = %d, A-S = %d = 4*%d, B-S = %d = 4*%d' % (S, A - S, W['x'][0], B - S, W['y'][0]))
    chk('A14-x-and-y-match-the-paper-and-are-odd',
        x == W['x'][0] and y == W['y'][0] and x % 2 == 1 and y % 2 == 1,
        'x = %d (odd), y = %d (odd)' % (x, y))
    chk('A15-the-two-subtracted-sums-are-even',
        sum(PA) == W['Sigma_A'][0] and sum(PB) == W['Sigma_B'][0]
        and sum(PA) % 2 == 0 and sum(PB) % 2 == 0,
        'Sigma_A = %d, Sigma_B = %d' % (sum(PA), sum(PB)))

    # =====================================================================
    print()
    print('=== B. THE SIX ACHIEVABLE VALUES OF EACH BRACKET ===')
    # =====================================================================
    def bracket_multiset(P):
        vals = {}
        for pp in permutations(P):
            vals[bracket(*pp)] = vals.get(bracket(*pp), 0) + 1
        return vals

    mA = bracket_multiset(PA)
    mB = bracket_multiset(PB)
    chk('B1-bracket-over-A-takes-exactly-6-values-each-4-times',
        len(mA) == 6 and set(mA.values()) == {4} and sum(mA.values()) == 24,
        'C(4,2) = 6 values, 4 of the 24 orderings each')
    chk('B2-bracket-over-B-takes-exactly-6-values-each-4-times',
        len(mB) == 6 and set(mB.values()) == {4} and sum(mB.values()) == 24)
    chk('B3-the-six-values-over-A-match-the-paper',
        sorted(mA) == W['brackets_A'], '%s' % sorted(mA))
    chk('B4-the-six-values-over-B-match-the-paper',
        sorted(mB) == W['brackets_B'], '%s' % sorted(mB))
    chk('B5-all-twelve-achievable-bracket-values-are-odd',
        all(v % 2 == 1 for v in list(mA) + list(mB)),
        'each is a sum of three products of odd primes')
    sp = W['spot_A_low']
    sp2 = W['spot_A_high']
    chk('B6-the-two-spot-computations-over-A',
        sp[0] * sp[1] + sp[2] == sp[3] and sp2[0] * sp2[1] + sp2[2] == sp2[3]
        and sp[3] == min(mA) and sp2[3] == max(mA),
        '%d*%d+%d = %d (min) and %d*%d+%d = %d (max)' % tuple(sp + sp2))
    sp3 = W['spot_B_high']
    chk('B7-the-spot-computation-for-the-maximum-over-B',
        sp3[0] * sp3[1] + sp3[2] == sp3[3] == max(mB),
        '%d*%d+%d = %d = max over B' % tuple(sp3))

    # =====================================================================
    print()
    print('=== C. THE SIXTEEN ADMISSIBLE EQUATIONS ON THE WITNESS ===')
    # =====================================================================
    eqs = sixteen_equations(x, y)
    printed_x = ('x-1 = B_A - Sigma_B', 'x-1', x - 1, 'A', 'B')
    printed_y = ('y-1 = B_B - Sigma_A', 'y-1', y - 1, 'B', 'A')
    chk('C1-there-are-sixteen-admissible-equations-and-both-printed-ones-are-among-them',
        len(eqs) == 16 and len(set(l for l, _, _, _, _ in eqs)) == 16
        and printed_x in eqs and printed_y in eqs,
        '4 left-hand sides x 2 brackets x 2 subtracted sums')
    hits, cx, cy, cboth, cswap, ctrue = sweep_576(n, PA, PB)
    chk('C2-printed-x-equation-fails-on-all-576-labellings', cx == 0,
        'x-1 = [ab(c+d)+cd] - (e+f+g+h): 0 of 576')
    chk('C3-printed-y-equation-fails-on-all-576-labellings', cy == 0,
        'y-1 = [ef(g+h)+gh] - (a+b+c+d): 0 of 576')
    chk('C4-the-two-printed-equations-fail-simultaneously', cboth == 0, '0 of 576')
    chk('C5-the-bracket-exchanged-variant-also-fails', cswap == 0,
        'x-1 = B_B - Sigma_B and y-1 = B_A - Sigma_A: 0 of 576')
    req_even = [(L + (sum(PA) if sn == 'A' else sum(PB))) % 2
                for _lab, ln, L, _bn, sn in eqs if ln in ('x-1', 'y-1')]
    chk('C6-the-eight-offset-carrying-equations-die-by-parity',
        len(req_even) == 8 and set(req_even) == {0},
        'each demands an EVEN bracket value while all twelve achievable values are odd')
    req_free = sorted(set(L + (sum(PA) if sn == 'A' else sum(PB))
                          for _lab, ln, L, _bn, sn in eqs if ln in ('x', 'y')))
    chk('C7-the-eight-offset-free-equations-demand-unachievable-values',
        req_free == W['required_offset_free']
        and all(v % 2 == 1 for v in req_free)
        and not (set(req_free) & (set(mA) | set(mB))),
        'required %s, all ODD hence parity-admissible, none achievable' % req_free)
    chk('C8-the-gap-at-the-closest-near-miss',
        min(req_free) - max(v for v in mA if v < min(req_free)) == W['gap_to_798279'][0]
        and max(v for v in mA if v < min(req_free)) == 798279,
        '807551 - 798279 = %d' % W['gap_to_798279'][0])
    chk('C9-the-smallest-overshoot-of-the-maximum-over-B',
        min(req_free) - max(mB) == W['overshoot_over_703583'][0]
        and all(v > max(mB) for v in req_free),
        '807551 - 703583 = %d; all four required values exceed max over B'
        % W['overshoot_over_703583'][0])
    chk('C10-brute-force-all-sixteen-equations-x-all-576-labellings',
        set(hits.values()) == {0} and len(hits) == 16,
        '%d integer comparisons, 0 satisfied' % (16 * 576))

    # =====================================================================
    print()
    print('=== D. PROPOSITION 4, THE CORRECTED IDENTITY ===')
    # =====================================================================
    E = esym(PA)
    F = esym(PB)
    chk('D1-elementary-symmetric-functions-of-A-primes-match-the-paper',
        [E[1], E[2], E[3]] == [W['E1'][0], W['E2'][0], W['E3'][0]],
        'E1 = %d, E2 = %d, E3 = %d, E4 = %d' % (E[1], E[2], E[3], E[4]))
    chk('D2-elementary-symmetric-functions-of-B-primes-match-the-paper',
        [F[1], F[2], F[3]] == [W['F1'][0], W['F2'][0], W['F3'][0]],
        'F1 = %d, F2 = %d, F3 = %d, F4 = %d' % (F[1], F[2], F[3], F[4]))
    chk('D3-x-equals-F3+F1-E2-1', x == F[3] + F[1] - E[2] - 1,
        '%d + %d - %d - 1 = %d = x' % (F[3], F[1], E[2], x))
    chk('D4-y-equals-E3+E1-F2-1', y == E[3] + E[1] - F[2] - 1,
        '%d + %d - %d - 1 = %d = y' % (E[3], E[1], F[2], y))
    chk('D5-FORCED-POSITIVE-the-corrected-identity-holds-on-all-576-labellings',
        ctrue == 576,
        'the same evaluator that returns 0 for the printed equations returns 576/576 here')
    T_E = prod(p + 1 for p in PA)
    T_F = prod(p + 1 for p in PB)
    chk('D6-the-sigma-forced-relation-prod(p+1)=prod(q+1)', T_E == T_F,
        'both equal %d; this is the only non-arithmetic input to Prop. 4' % T_E)
    chk('D7-phi(M)=T-2(E3+E1)-and-phi(N)=T-2(F3+F1)',
        prod(p - 1 for p in PA) == T_E - 2 * (E[3] + E[1])
        and prod(p - 1 for p in PB) == T_F - 2 * (F[3] + F[1]),
        'the expansion step of the proof of Prop. 4')
    chk('D8-E4-T-equals-minus(E3+E2+E1+1)', E[4] - T_E == -(E[3] + E[2] + E[1] + 1),
        'the last step of the proof of Prop. 4')

    # =====================================================================
    print()
    print('=== E. LEMMA 1 RE-TESTED AS AN OUTPUT, ON A FAMILY WITH NO AMICABLE PAIR ===')
    # =====================================================================
    odd_primes_40 = [p for p in range(3, 1000) if is_prime(p)][:40]
    odd_primes_20 = odd_primes_40[:20]
    second_20 = odd_primes_40[20:]
    quads = list(combinations(odd_primes_20, 4))
    quads_b = list(combinations(second_20, 4))
    chk('E1-16-divides-phi(M)-for-every-squarefree-M-with-four-odd-prime-factors',
        len(quads) == 4845 and all(prod(p - 1 for p in q) % 16 == 0 for q in quads),
        'C(20,4) = %d subsets of {%d,...,%d}, no amicability used'
        % (len(quads), odd_primes_20[0], odd_primes_20[-1]))
    chk('E2-phi(M)=E4-E3+E2-E1+1-for-every-such-M',
        all(prod(p - 1 for p in q) ==
            esym(list(q))[4] - esym(list(q))[3] + esym(list(q))[2] - esym(list(q))[1] + 1
            for q in quads),
        'the expansion used in Lemma 1 and Prop. 4, on %d subsets' % len(quads))
    # ⭐ Drawn from two DISJOINT prime pools, so every pair is coprime by construction and
    # the family size is deterministic -- no rejection sampling, no dependence on luck.
    fam = [(quads[(7 * i + 3) % len(quads)], quads_b[(263 * i + 101) % len(quads_b)])
           for i in range(2000)]
    bad_int = bad_odd = 0
    for qa, qb in fam:
        for nn in range(1, 9):
            Aa = 2 ** nn * prod(qa)
            Bb = 2 ** nn * prod(qb)
            Ss = phi_from(nn, qa) + phi_from(nn, qb)
            if (Aa - Ss) % 2 ** nn or (Bb - Ss) % 2 ** nn:
                bad_int += 1
                continue
            if ((Aa - Ss) // 2 ** nn) % 2 != 1 or ((Bb - Ss) // 2 ** nn) % 2 != 1:
                bad_odd += 1
    chk('E3-x-and-y-are-integers-and-both-odd-across-the-whole-family',
        len(fam) == 2000 and bad_int == 0 and bad_odd == 0,
        '%d coprime prime-set pairs x n = 1..8 = %d instances, 0 non-integral, 0 even; '
        'not one is an amicable pair' % (len(fam), 8 * len(fam)))
    bad_par = 0
    for qa, qb in fam[:400]:
        for pa in permutations(qa):
            if bracket(*pa) % 2 != 1 or sum(pa) % 2 != 0:
                bad_par += 1
    chk('E4-the-printed-bracket-is-always-odd-and-the-subtracted-sum-always-even',
        bad_par == 0,
        '400 subsets x 24 orderings = 9600 instances; so B - Sigma is ODD and forces x EVEN')
    pairs2 = list(combinations(odd_primes_20, 2))
    chk('E5-two-distinct-odd-prime-factors-already-suffice-for-the-lemma',
        all(prod(p - 1 for p in q) % 4 == 0 for q in pairs2)
        and all((2 ** 1 * prod(a) - (phi_from(1, a) + phi_from(1, b))) // 2 % 2 == 1
                for a, b in zip(pairs2[:50], pairs2[50:100])),
        'C(20,2) = %d subsets: 4 | phi(M), and x odd on 50 sample pairs at n = 1' % len(pairs2))
    m1, n1 = [3], [5]
    x1 = prod(m1) - (prod(p - 1 for p in m1) + prod(p - 1 for p in n1)) // 2
    chk('E6-the-hypothesis-is-not-vacuous',
        x1 % 2 == 0,
        'with ONE odd prime factor each (M=3, N=5) the same formula gives x = %d, EVEN' % x1)

    # =====================================================================
    print()
    print('=== F. CONSISTENCY WITH THE THREE IDENTITIES THE TARGET PAPER PROVES ===')
    # =====================================================================
    e = EX22
    pA2, pB2 = e['primes_A'], e['primes_B']
    A2 = 4 * prod(pA2)
    B2 = 4 * prod(pB2)
    S2 = phi_from(2, pA2) + phi_from(2, pB2)
    x2 = (A2 - S2) // 4
    y2 = (B2 - S2) // 4
    chk('F1-the-target-papers-own-(2,2)-example-reproduces-its-Theorem-2.1',
        A2 == e['A'][0] and B2 == e['B'][0]
        and sigma_from(2, pA2) == sigma_from(2, pB2) == A2 + B2 == e['sigma'][0]
        and S2 == e['S'][0] and x2 == e['x'][0] == sum(pB2) - 1
        and y2 == e['y'][0] == sum(pA2) - 1 and x2 % 2 == 1 and y2 % 2 == 1,
        'A = %d, B = %d: x = %d = c+d-1, y = %d = a+b-1, both odd' % (A2, B2, x2, y2))
    chk('F2-Theorem-2.1-shape-c+d-1-is-always-odd',
        all((a + b - 1) % 2 == 1 for a, b in pairs2),
        'odd + odd - 1 on all %d two-subsets' % len(pairs2))
    tri = list(combinations(odd_primes_20, 3))
    chk('F3-Theorem-2.3-shape-(d+e)-(a+b+c)-is-always-odd',
        all((d + ee - sum(t)) % 2 == 1 for t in tri[:200] for d, ee in pairs2[:20]),
        'even - odd on %d instances' % (200 * 20))
    chk('F4-Theorem-2.5-shape-x-1=e2(d,e,f)-(a+b+c)-is-always-even',
        all((esym(list(u))[2] - sum(t)) % 2 == 0 for t in tri[:200] for u in tri[:20]),
        'odd - odd on %d instances, so x is ODD there too' % (200 * 20))
    chk('F5-the-monomial-counts-that-locate-the-error',
        len(list(combinations('abc', 2))) == 3 and len(list(combinations('abcd', 3))) == 4
        and len(list(combinations('abcd', 2))) == 6,
        'e2 of 3 letters: 3 monomials (ODD); e3 of 4: 4 (EVEN); e2 of 4: 6. '
        '[ab(c+d)+cd] has 3 monomials, so it carries the 3-letter parity to 4 letters')
    chk('F6-the-corrected-form-restores-the-parity',
        all((esym(list(qb))[3] + esym(list(qb))[1] - esym(list(qa))[2] - 2) % 2 == 0
            for qa, qb in fam[:500]),
        'x-1 = F3+F1-E2-2 is EVEN on 500 instances, as Lemma 1 demands')

    # =====================================================================
    print()
    print('=== G. THE EIGHT TABULATED MEMBERS OF THE CELL, RE-DERIVED FROM SCRATCH ===')
    # =====================================================================
    for idx, (nn, qa, qb, Aprinted) in enumerate(TAB, 1):
        Aa = 2 ** nn * prod(qa)
        Bb = 2 ** nn * prod(qb)
        xx, yy, _Ss, _, _ = xy_of(nn, qa, qb)
        h, c1, c2, c3, c4, c5 = sweep_576(nn, qa, qb)
        Ea, Fb = esym(qa), esym(qb)
        ok = (Aa == Aprinted
              and all(is_prime(p) for p in qa + qb) and len(set(qa + qb)) == 8
              and sigma_from(nn, qa) == sigma_from(nn, qb) == Aa + Bb and Aa != Bb
              and gcd(Aa, Bb) == 2 ** nn
              and xx % 2 == 1 and yy % 2 == 1
              and set(h.values()) == {0} and c1 == c2 == c3 == c4 == 0
              and c5 == 576
              and xx == Fb[3] + Fb[1] - Ea[2] - 1 and yy == Ea[3] + Ea[1] - Fb[2] - 1)
        chk('G%d-cell-member-n=%d-A=%d' % (idx, nn, Aa), ok,
            'amicable, gcd = 2^%d exactly, type (4,4), x = %d and y = %d both odd, '
            '0 of 576 labellings satisfy any of the 16 equations, corrected identity 576/576'
            % (nn, xx, yy))
    chk('G9-every-tabulated-A-equals-2^n-times-its-prime-list',
        all(2 ** nn * prod(qa) == Ap_ for nn, qa, _qb, Ap_ in TAB),
        '8 of 8 rows transcribe consistently')
    chk('G10-the-first-tabulated-row-is-the-witness-of-Theorem-3',
        TAB[0][0] == n and TAB[0][1] == PA and TAB[0][2] == PB and TAB[0][3] == Ap)
    chk('G11-the-table-spans-every-n-that-occurs-in-the-census',
        sorted(set(r[0] for r in TAB)) == list(range(2, 9)),
        'n = 2,...,8, matching the census distribution')

    # =====================================================================
    print()
    print('=== H. THE CENSUS FIGURES QUOTED IN SECTION 5 (INTERNAL ARITHMETIC ONLY) ===')
    # =====================================================================
    chk('H1-the-n-distribution-sums-to-the-cell-size',
        sum(CEN['n_distribution']) == CEN['cell_members'][0]
        and len(CEN['n_distribution']) == 7,
        '%s = %d' % (' + '.join(map(str, CEN['n_distribution'])), CEN['cell_members'][0]))
    chk('H2-labellings-per-pair-times-cell-size',
        CEN['labellings_per_pair'][0] == 24 * 24
        and CEN['labellings_per_pair'][0] * CEN['cell_members'][0] == CEN['labellings_total'][0],
        '4! x 4! = 576, and 576 * %d = %d' % (CEN['cell_members'][0], CEN['labellings_total'][0]))
    chk('H3-the-witness-is-the-smallest-A-among-the-tabulated-members',
        Ap == min(r[3] for r in TAB),
        'consistent with the census claim that it is the least member of the cell; '
        'the census itself is NOT re-run here')

    # =====================================================================
    print()
    print('=== SCOPE ===')
    print('NOT RE-RUN: the 415,523-term census of OEIS A002025. It needs the b-file, which is an '
          'external data file this program is forbidden to read, so the figures 14,484 cell members, '
          '8,342,784 labellings, the n-distribution and the 1,073 type-(4,2) members are checked here '
          'only for internal arithmetic consistency (section H). No statement of the paper rests on '
          'them: Theorem 2 is a parity argument and Theorem 3 is one explicit pair.')
    print('NOT RE-RUN: the two minimality claims that do rest on the census -- that A = 32642324 is the '
          'SMALLEST type-(4,4) amicable pair with gcd a power of two, and that no member of the cell has '
          'n = 1 or n > 8 below 10^17. Section H checks only that the witness is smallest among the eight '
          'pairs printed in the paper.')
    print('NOT RE-RUN: the provenance of the target. The wording of the conjecture, the line numbers '
          '576-582, 565, 222-230, 567-573, 585, 21, 33, the file size 25,122 B and the 664 lines of '
          'Amicable27122025.tex were established by fetching and reading the arXiv:2512.22319v1 e-print '
          'source; nothing here re-fetches it. This program checks mathematics, not bibliography.')
    print('NOT RE-RUN: the three theorems of the target paper are checked here only for PARITY SHAPE and '
          'on the one (2,2) example the paper prints (section F). Their proofs are not reproduced, and '
          'the paper\'s (3,3) example is not recomputed.')
    print('NOT RE-RUN: the offset-free equations for type-(4,4) pairs OTHER than the eight tabulated '
          'here. Parity does not reach them and no general argument is offered; they are open, as the '
          'paper says in its Status list.')
    print('NOT RE-RUN: the prior-art question. Whether Proposition 4 is already in print was not '
          'settled -- the full texts of Borho-Hoffmann, Costello and Garcia-Pedersen-te Riele were not '
          'obtained, and MathSciNet was never consulted.')
    print('NOT COVERED: the assertion at line 585 of the target, that x is in general a function of the '
          '(k-1)-th and (k-2)-th elementary symmetric polynomials of its OWN prime factors. Nothing '
          'here refutes it.')

    print()
    if _STATE['fail']:
        print('VERDICT: %d of %d CHECKS FAILED' % (_STATE['fail'], _STATE['fail'] + _STATE['pass']))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _STATE['pass'])
    return 0


if __name__ == '__main__':
    sys.exit(main())

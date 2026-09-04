#!/usr/bin/env python3
"""Verification program for

    "A mod 31^2 counterexample to a generalized cubic partition congruence of
     Das, Maity and Saikia"

Python 3.9+, STANDARD LIBRARY ONLY (no numpy, no sympy, no external data file).
Exact integer arithmetic throughout; no floating point value is ever compared or branched on.

WHAT IT CHECKS.  Every quantity printed in the paper is re-derived here from the definition, and the
integers printed in the paper are pasted below and compared against the fresh computation as decimal
STRINGS, so a transcription error in the paper fails this program.

    statement under test (Das-Maity-Saikia, arXiv:2503.19399v2, Conjecture 7.1, first line):

            for all n >= 0,     a_25(31^2 n + 644) == 0   (mod 31^2)

    a_c(n) is defined by the source's own generating function

            sum_{n>=0} a_c(n) q^n  =  1 / (f_1 f_2^{c-1}),      f_k = (q^k; q^k)_infinity

    i.e. odd parts monochromatic, even parts in c colours.

TWO DISJOINT EXACT ROUTES to the witness, sharing no arithmetic:

  R-A  q-PRODUCT ROUTE.  D = f_1 f_2^{24} truncated, built by sparse multiplications indexed by the
       generalized pentagonal numbers; then A = 1/D from A[0] = 1, A[n] = -sum_{j=1..n} D[j] A[n-j].
       D[0] = 1, so no division and no modular inverse is ever needed.
  R-B  LOGARITHMIC-DERIVATIVE ROUTE.  From q d/dq log(1/(f_1 f_2^{c-1})),
            n a_c(n) = sum_{j=1..n} B(j) a_c(n-j),   B(j) = sigma(j) + 2(c-1) sigma(j/2) [2 | j],
       exact integers, every division by n asserted exact.  Forms no q-product, inverts nothing and
       never touches a pentagonal number.
  R-C  a third, purely COMBINATORIAL route at tiny n: a part of even size used k times may be
       coloured in C(c+k-1, k) ways, an odd part in one way.  Used to pin a_25(0..6) with no series
       machinery at all.

CONTROLS, BOTH POLARITIES.  Published third-party integers (OEIS A000041 and A002513) pin the
definition on numbers nobody here minted; the eight congruences the source PROVES must stay silent
(its Theorem 1.3, all seven lines, and its Theorem 1.4, which is a mod p^2 statement built by the
same Radu machinery); every residue class mod 31^2 must be live on a_25 while exactly the class
644 mod 31 = 24 is dead, the latter being a PUBLISHED theorem of Guadalupe; and the source's PROSE
reading 1/f_1^c (which contradicts its own generating function) is run as an anti-control.

Exit code 0 iff every check passed.
"""
import sys

FAILS = []
NCHECK = 0


def out(*a):
    print(*a)
    sys.stdout.flush()


def ck(name, cond, detail=''):
    global NCHECK
    NCHECK += 1
    if not cond:
        FAILS.append(name)
    out('%s %s%s' % ('PASS' if cond else 'FAIL', name, ('  | ' + detail) if detail else ''))


# ---------------------------------------------------------------------------
# generalized pentagonal exponents of f_step = (q^step; q^step)_inf, as (exponent, sign)
# ---------------------------------------------------------------------------
def pent(step, L):
    terms = [(0, 1)]
    k = 1
    while True:
        e1 = step * (k * (3 * k - 1) // 2)
        e2 = step * (k * (3 * k + 1) // 2)
        if e1 >= L and e2 >= L:
            break
        s = -1 if k % 2 else 1
        if e1 < L:
            terms.append((e1, s))
        if e2 < L:
            terms.append((e2, s))
        k += 1
    return terms


# ---------------------------------------------------------------------------
# R-A, exact over Z
# ---------------------------------------------------------------------------
def sparse_mul_exact(dense, terms, L):
    res = [0] * L
    for (e, s) in terms:
        if e >= L:
            continue
        if s == 1:
            for i in range(L - e):
                res[e + i] += dense[i]
        else:
            for i in range(L - e):
                res[e + i] -= dense[i]
    return res


def series_exact_qproduct(c, L):
    D = [0] * L
    D[0] = 1
    D = sparse_mul_exact(D, pent(1, L), L)
    t2 = pent(2, L)
    for _ in range(c - 1):
        D = sparse_mul_exact(D, t2, L)
    assert D[0] == 1
    A = [0] * L
    A[0] = 1
    for n in range(1, L):
        acc = 0
        for j in range(1, n + 1):
            dj = D[j]
            if dj:
                acc += dj * A[n - j]
        A[n] = -acc
    return A


# ---------------------------------------------------------------------------
# R-B, exact over Z
# ---------------------------------------------------------------------------
def series_exact_divisorsum(c, L):
    sig = [0] * L
    for d in range(1, L):
        for m in range(d, L, d):
            sig[m] += d
    B = [0] * L
    for j in range(1, L):
        B[j] = sig[j] + (2 * (c - 1) * sig[j // 2] if j % 2 == 0 else 0)
    A = [0] * L
    A[0] = 1
    for n in range(1, L):
        acc = 0
        for j in range(1, n + 1):
            bj = B[j]
            if bj:
                acc += bj * A[n - j]
        q, r = divmod(acc, n)
        assert r == 0, 'division by n=%d was not exact -- the recurrence is wrong' % n
        A[n] = q
    return A


# ---------------------------------------------------------------------------
# R-C, combinatorial, tiny n
# ---------------------------------------------------------------------------
def binom(n, k):
    num = 1
    for i in range(k):
        num = num * (n - i) // (i + 1)
    return num


def series_combinatorial(c, L):
    """a_c(0..L-1) by a DP over part sizes.  An even part of size k used j times contributes
    C(c+j-1, j) colourings; an odd part contributes 1.  No q-product, no recurrence, no inversion."""
    A = [0] * L
    A[0] = 1
    for k in range(1, L):
        new = [0] * L
        for n in range(L):
            if A[n] == 0:
                continue
            j = 0
            while n + j * k < L:
                w = 1 if k % 2 else binom(c + j - 1, j)
                new[n + j * k] += A[n] * w
                j += 1
        A = new
    return A


# ---------------------------------------------------------------------------
# THE MODULAR ENGINE: never forms f_1 f_2^{c-1} at all.  1/(f_1 f_2^{c-1}) is obtained by dividing 1
# by f_1 and then by f_2, (c-1) times, each division by a SPARSE series costing O(L * #terms).
# ---------------------------------------------------------------------------
def divide_by_sparse(H, terms, M):
    """G with G * f == H (mod M), where f = 1 + sum_{e>0} s_e q^e is given by `terms`."""
    L = len(H)
    plus = sorted(e for (e, s) in terms if s == 1 and e > 0)
    minus = sorted(e for (e, s) in terms if s == -1)
    G = list(H)
    for n in range(1, L):
        acc = G[n]
        for e in plus:
            if e > n:
                break
            acc -= G[n - e]
        for e in minus:
            if e > n:
                break
            acc += G[n - e]
        G[n] = acc % M
    return G


def series_mod(c, L, M, prose=False):
    """a_c(0..L-1) mod M.  prose=True computes the source's PROSE object 1/f_1^c instead."""
    A = [0] * L
    A[0] = 1
    t1 = pent(1, L)
    if prose:
        for _ in range(c):
            A = divide_by_sparse(A, t1, M)
    else:
        A = divide_by_sparse(A, t1, M)
        t2 = pent(2, L)
        for _ in range(c - 1):
            A = divide_by_sparse(A, t2, M)
    return A


def v_p(x, p):
    if x == 0:
        return None
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


out('verification of: a mod 31^2 counterexample to a generalized cubic partition congruence')
out('of Das, Maity and Saikia (arXiv:2503.19399v2, Conjecture 7.1, first displayed congruence)')
out('python %s' % sys.version.split()[0])
out('exact integer arithmetic only; standard library only')
out('')

# ===========================================================================
# Section 1.  The integers PRINTED IN THE PAPER, and their internal arithmetic
# ===========================================================================
out('=== Section 1: the integers printed in the paper')
PAPER_WITNESS = (
    '14714658623812815578026685166319028447102593034654589643689045403133888112437860270187'
    '17043793415486827990833622853525010333330703396635284385120674666433997413001085330998'
    '1128071099'
)
PAPER_QUOTIENT = (
    '47466640721976824445247371504254930474524493660176095624803372268173832620767291194152'
    '16270301340280090293011686624274226881711946440758981887486047311077411009680920422574'
    '55744229'
)
W = int(PAPER_WITNESS)
Q = int(PAPER_QUOTIENT)
ck('paper_witness_has_182_digits', len(PAPER_WITNESS) == 182, '%d digits' % len(PAPER_WITNESS))
ck('paper_quotient_has_180_digits', len(PAPER_QUOTIENT) == 180, '%d digits' % len(PAPER_QUOTIENT))
ck('paper_integers_are_arithmetically_consistent',
   31 * Q == W and W % 961 == 682 and W % 31 == 0 and Q % 31 == 22,
   '31*quotient == witness; witness mod 961 = %d; witness mod 31 = %d; quotient mod 31 = %d'
   % (W % 961, W % 31, Q % 31))
ck('the_index_is_2566', 31 ** 2 * 2 + 644 == 2566, 'N = 31^2 * 2 + 644 = 2566')
ck('682_is_22_times_31', 682 == 22 * 31 and 682 % 961 != 0,
   '682 = 22*31, nonzero modulo 961')

# ===========================================================================
# Section 2.  The witness, exactly, by two disjoint exact routes
# ===========================================================================
out('')
out('=== Section 2: the witness a_25(2566), exactly, by two disjoint exact routes')
NW = 2566
LW = NW + 1
EA = series_exact_qproduct(25, LW)
EB = series_exact_divisorsum(25, LW)
ck('route_A_equals_route_B_at_every_coefficient', EA == EB,
   'exact agreement over all %d coefficients a_25(0..%d)' % (LW, NW))
a2566 = EA[NW]
ck('fresh_computation_equals_the_witness_printed_in_the_paper', str(a2566) == PAPER_WITNESS,
   '%d digits, compared as a decimal string' % len(str(a2566)))
ck('fresh_quotient_equals_the_quotient_printed_in_the_paper',
   a2566 % 31 == 0 and str(a2566 // 31) == PAPER_QUOTIENT, '')
ck('a_25_2566_is_682_mod_961', a2566 % 961 == 682, 'a_25(2566) mod 31^2 = %d' % (a2566 % 961))
ck('v_31_of_a_25_2566_is_exactly_1', v_p(a2566, 31) == 1,
   '31 divides a_25(2566) and 31^2 does not')
ck('CONJECTURE_7_1_LINE_1_IS_FALSE_AT_n_equals_2', a2566 % 961 != 0,
   'a_25(31^2*2 + 644) = a_25(2566) is NOT congruent to 0 modulo 31^2')

out('NOTE a_25(644)  has %d digits, v_31 = %s, residue mod 961 = %d'
    % (len(str(EA[644])), v_p(EA[644], 31), EA[644] % 961))
out('NOTE a_25(1605) has %d digits, v_31 = %s, residue mod 961 = %d'
    % (len(str(EA[1605])), v_p(EA[1605], 31), EA[1605] % 961))
ck('the_cells_n_0_and_n_1_DO_vanish_mod_961_with_v_31_exactly_2',
   EA[644] % 961 == 0 and EA[1605] % 961 == 0
   and v_p(EA[644], 31) == 2 and v_p(EA[1605], 31) == 2
   and len(str(EA[644])) == 84 and len(str(EA[1605])) == 141,
   'so v_31 is NOT identically 1 on the progression, and a spot check at n = 0, 1 would have '
   'supported the conjecture')

# ===========================================================================
# Section 3.  The definition, pinned three ways
# ===========================================================================
out('')
out('=== Section 3: the definition of a_c, pinned by a third route and by published integers')
SMALL = [1, 1, 26, 27, 377, 403, 3979, 4355, 33956, 37908, 247690, 281242, 1597895]
out('NOTE a_25(0..12) = %s' % SMALL)
ck('the_small_values_printed_in_the_paper_are_reproduced', EA[:13] == SMALL, '')
CO = series_combinatorial(25, 7)
ck('a_third_purely_combinatorial_route_agrees_at_n_0_to_6', CO == SMALL[:7],
   'colourings counted directly: a(2) = 1 + 25, a(3) = 1 + 25 + 1, '
   'a(4) = 1 + 25 + C(26,2) + 1 + 25 = 377')
P1 = series_exact_divisorsum(1, 201)
ck('c_equals_1_reproduces_OEIS_A000041',
   [P1[50], P1[100], P1[200]] == [204226, 190569292, 3972999029388],
   'p(50) = %d, p(100) = %d, p(200) = %d' % (P1[50], P1[100], P1[200]))
A002513 = [1, 1, 3, 4, 9, 12, 23, 31, 54, 73, 118, 159, 246, 329, 489, 651, 940, 1242, 1751,
           2298, 3177]
P2 = series_exact_divisorsum(2, 21)
ck('c_equals_2_reproduces_OEIS_A002513', P2 == A002513, 'a_2(0..20) = %s' % P2)
PR = series_mod(25, LW, 961, prose=True)
ck('the_sources_PROSE_reading_1_over_f1_to_the_c_is_a_different_object',
   PR[644] % 961 == 550 and PR[2566] % 961 == 492,
   'under 1/f_1^25 the target would already fail at n = 0: residues 550 and 492 mod 961')

# ===========================================================================
# Section 4.  The window printed in the paper, modulo 31^2
# ===========================================================================
out('')
out('=== Section 4: a_25(961 n + 644) mod 961 for n = 0..12, and the residue-class sweeps')
NMAX = 12
L961 = 961 * NMAX + 644 + 1
A961 = series_mod(25, L961, 961)
ck('the_modular_engine_agrees_with_the_exact_routes_on_0_to_2566',
   all(A961[i] == EA[i] % 961 for i in range(LW)),
   'a third independent engine, on the same coefficients')
row = [A961[961 * n + 644] for n in range(NMAX + 1)]
out('NOTE a_25(961n+644) mod 961, n = 0..12: ' + ' '.join('%d:%d' % (n, v) for n, v in enumerate(row)))
ck('the_window_printed_in_the_paper_is_reproduced',
   row == [0, 0, 682, 217, 837, 124, 31, 403, 868, 465, 0, 124, 775], '')
ck('the_vanishing_cells_of_the_window_are_exactly_n_0_1_10',
   [n for n, v in enumerate(row) if v == 0] == [0, 1, 10], '')
ck('every_nonvanishing_cell_of_the_window_is_a_nonzero_multiple_of_31',
   all(v % 31 == 0 for v in row) and len([v for v in row if v]) == 10,
   '10 nonzero cells, all divisible by 31: Guadalupe Theorem 3 at p = 31 is satisfied while the '
   'mod 31^2 claim fails')
ck('the_neighbour_residues_printed_in_the_paper_are_reproduced',
   [A961[r] for r in (641, 642, 643, 644, 645, 646, 647)] == [561, 336, 174, 0, 431, 817, 159],
   'a_25(641..647) mod 961')
dead961 = [r for r in range(961) if not any(A961[r::961])]
ck('no_residue_class_mod_961_is_identically_zero_on_a_25', dead961 == [],
   'the decider is not stuck on "zero": 961 of 961 classes carry a nonzero value over the swept '
   'range, and it says YES on r = 644 at n = 2')
dead31 = [r for r in range(31) if not any(x % 31 for x in A961[r::31])]
ck('exactly_the_class_24_mod_31_is_identically_zero_on_a_25', dead31 == [24],
   '644 mod 31 = %d, and 13(31^2-1)/24 = %d = 31*16 + 24. THIS IS A PUBLISHED THEOREM '
   '(Guadalupe, Integers 25 (2025) #A20, Thm 3), reproduced here, not discovered here'
   % (644 % 31, 13 * (31 ** 2 - 1) // 24))

# ===========================================================================
# Section 5.  Must-stay-silent controls: the congruences the source PROVES
# ===========================================================================
out('')
out('=== Section 5: must-stay-silent controls -- the eight congruences the source proves')
L49 = 49 * 130 + 39 + 1
A3 = series_mod(3, L49, 49)
viol = [n for n in range(131) if A3[49 * n + 39] % 49]
dead49 = [r for r in range(49) if not any(A3[r::49])]
ck('Theorem_1_4_a_3_49n_plus_39_is_0_mod_49_stays_silent', viol == [] and dead49 == [39],
   'n = 0..130, 0 violations; the unique identically-zero class of the 49 is [39]. This is the '
   'source own mod p^2 theorem, built by the SAME Radu machinery as the refuted conjecture')
THM1 = [(37, 43, 12), (41, 47, 21), (53, 59, 56), (61, 67, 19), (65, 71, 32), (73, 79, 62),
        (77, 83, 79)]
for (c, p, r) in THM1:
    NM = 40
    Ac = series_mod(c, p * NM + r + 1, p)
    v = [n for n in range(NM + 1) if Ac[p * n + r] % p]
    dead = [t for t in range(p) if not any(Ac[t::p])]
    ck('Theorem_1_3_a_%d_%dn_plus_%d_is_0_mod_%d_stays_silent' % (c, p, r, p),
       v == [] and dead == [r] and (13 * (p * p - 1) // 24) % p == r,
       'n = 0..%d, 0 violations; unique identically-zero class [%d] of the %d; and '
       '13(p^2-1)/24 mod p = %d, the published residue'
       % (NM, r, p, (13 * (p * p - 1) // 24) % p))

# ===========================================================================
# Section 6.  The two sibling congruences of the same conjecture environment
# ===========================================================================
out('')
out('=== Section 6: the other two displayed congruences of Conjecture 7.1')
for (c, p, r, want) in ((41, 47, 256, 1927), (61, 67, 555, 2747)):
    m = p * p
    N = m + r
    Ax = series_mod(c, N + 1, m)
    val = Ax[N]
    ck('sibling_a_%d_%d_squared_plus_%d_is_%d_mod_%d_squared' % (c, p, r, want, p),
       val == want and val % p == 0 and val != 0 and Ax[r] % m == 0 and want == (want // p) * p,
       'a_%d(%d) = %d = %d*%d mod %d: nonzero mod %d and zero mod %d, so v_%d = 1 exactly; and '
       'the n = 0 cell a_%d(%d) DOES vanish mod %d. This line is false at n = 1.'
       % (c, N, val, want // p, p, m, m, p, p, c, r, m))

# ===========================================================================
# Section 7.  The certificate window, from the source's own printed formula
# ===========================================================================
out('')
out('=== Section 7: the size of Radu certificate window for this progression')
PT = sorted((644 + 31 * u) % 961 for u in range(31) if legendre(1 + 24 * u, 31) == 1)
ck('P_of_644_has_fifteen_residues_not_one',
   PT == [55, 148, 303, 334, 365, 427, 489, 613, 644, 675, 706, 768, 799, 861, 892]
   and (24 * 644 - 49) % 961 == 31,
   'P(644) = %s; the source prints {644}. So Radu Lemma 4.5 asks for 15*3814 = %d cells, not '
   '3814 -- see the paper Scope section' % (PT, 15 * 3814))

# ===========================================================================
out('')
out('NOTE SCOPE: this program re-derives every integer printed in the paper and nothing beyond it. '
    'It does NOT prove the conjecture on any other progression, it does NOT discharge Radu '
    'certificate hypothesis for the mod 31 companion congruence (that congruence is a published '
    'theorem of Guadalupe and is quoted, not proved, here), and it does NOT verify the source '
    'printed Radu parameter table against Radu original paper -- Section 7 re-derives P(644) from '
    'the formula AS PRINTED IN THE SOURCE only. The residue-class sweeps of Section 4 are '
    'bounded-range statements over the swept window, not uniqueness theorems; Sections 5 and 6 run '
    'to n = 40 and n = 1 respectively.')
out('')
if FAILS:
    out('*** %d CHECK(S) FAILED: %s' % (len(FAILS), ', '.join(FAILS)))
    sys.exit(1)
out('VERDICT: ALL %d CHECKS PASS' % NCHECK)
sys.exit(0)

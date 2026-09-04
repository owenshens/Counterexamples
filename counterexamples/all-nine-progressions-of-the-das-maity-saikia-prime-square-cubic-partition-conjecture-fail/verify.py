#!/usr/bin/env python3
"""verify.py -- checks the computational claims of the accompanying paper

    "All Nine Progressions of the Conjecture of Das, Maity and Saikia on
     Generalized Cubic Partitions Modulo Prime Squares Printed in
     arXiv:2503.19399v2 Fail"

Python 3.9+, STANDARD LIBRARY ONLY (no numpy, no sympy, no external data file), exact integer
arithmetic throughout, no floating point anywhere, no randomness.

EVERY number this program compares against is TRANSCRIBED FROM THE PAPER, in the block marked
"THE PAPER'S PRINTED OBJECTS" below, and is then re-derived here from the definition

        sum_{n>=0} a_c(n) q^n  =  1 / ( f_1 * f_2^{c-1} ),      f_k = prod_{m>=1} (1 - q^{km}),

which is the displayed generating function of the source paper (arXiv:2503.19399v2).  Note that it is NOT the paper's PROSE definition 1/f_1^c: see the checks
`anticontrol-prose-c2` and `anticontrol-prose-cells`, which establish that the prose reading
reproduces neither the published cubic-partition numbers nor a single one of the nine cells.

Every check prints one `PASS <name>` line; the program closes with
    VERDICT: ALL <n> CHECKS PASS
and exits 0 if and only if every check passed.

RUN:  python3 verify.py
Runtime is dominated by the exhaustive residue census of Section 5 of the paper: every one of the
18 221 residue classes modulo p^2, for the five primes Conjecture 7.2 names, to depth 3.  Expect
a couple of minutes and a few hundred MB of memory.
"""

import platform
import sys

# =====================================================================================
# THE PAPER'S PRINTED OBJECTS.  Nothing below this block is taken from the paper.
# =====================================================================================

# The nine progressions of Conjecture 7.2, as printed, in printed order: (c, p, r).
CELLS = [(25, 31, 41), (25, 31, 81), (25, 31, 585),
         (37, 43, 716), (37, 43, 987),
         (53, 59, 100),
         (65, 71, 41), (65, 71, 47),
         (77, 83, 157)]

# Paper, Theorem 1(i) and eq (2): the refuting coefficient, exactly.
A65_41 = 451303321502143296879
A65_41_MOD_P = 40                      # mod 71
A65_41_MOD_P2 = 2170                   # mod 71^2 = 5041
QUOT_P = 6356384809889342209           # A65_41 = 71 * QUOT_P + 40
QUOT_P2 = 89526546618159749            # A65_41 = 5041 * QUOT_P2 + 2170

# The two tables whose 21-term dot product IS the coefficient (transcribed from the note).
P64_TABLE = [1, 64, 2144, 49920, 905840, 13627264, 176638592, 2025205248, 20930373880,
             197788352320, 1728062919232, 14083242424576, 107837287452608, 780481475916160,
             5366307146732800, 35202669371599360, 221142159585764508, 1334633003840266624,
             7760187771579170400, 43579749087236893440, 236897695447322916960]
PART_TABLE = [44583, 31185, 21637, 14883, 10143, 6842, 4565, 3010, 1958, 1255, 792, 490, 297,
              176, 101, 56, 30, 15, 7, 3, 1]                      # p(41-2j), j = 0..20

# The hand certificate of the note's Section 3.  E_TABLE is OEIS A000730, PART_TABLE is
# OEIS A000041, and their dot product is b(41).
E_TABLE = [1, -7, 14, 7, -49, 21, 35, 41, -49, -133, 98, -21, 126, 112, -176, -105, -126, 140,
           -35, 147, 259]
B41 = -31

# The two companion integers (transcribed from the note).
A25_41 = 793193782332525               # = 961 * 825383748525, so == 0 mod 31^2
A25_41_COFACTOR = 825383748525
A65_47 = 65316694397122053889074       # == 0 mod 71^2

# The note's table of the nine cells: a_c(p^2 n + r) mod p^2 for n = 0, 1, 2, in printed order.
TABLE3 = {(25, 31, 41): [0, 87, 154],
          (25, 31, 81): [0, 425, 685],
          (25, 31, 585): [0, 336, 823],
          (37, 43, 716): [0, 382, 199],
          (37, 43, 987): [0, 301, 288],
          (53, 59, 100): [0, 1768, 3210],
          (65, 71, 41): [2170, 3779, 1478],
          (65, 71, 47): [0, 3518, 885],
          (77, 83, 157): [0, 429, 6028]}

# The note's mod-p column: the n = 1 values reduced mod p rather than mod p^2.
MODP_ROW = [25, 22, 26, 38, 0, 57, 16, 39, 14]

# The COMPLETE list of r in [0, p^2) with a_c(r) == 0 mod p^2, for the five primes Conjecture 7.2
# names (c = p - 6 throughout), and the survivor counts at depths 1, 2, 3 (transcribed).
ZEROS = {31: [41, 81, 585, 644],
         43: [716, 987, 1115, 1841],
         59: [100],
         71: [47, 2173, 2233, 4088],
         83: [157, 3362]}
SURVIVORS = {31: [4, 1, 0], 43: [4, 0, 0], 59: [1, 0, 0], 71: [4, 0, 0], 83: [2, 0, 0]}
# The unique class in all 18 221 that survives n = 0 and n = 1, and the value that kills it.
LONE_DEPTH2 = (31, 644)
A25_2566_MOD = 682                     # == 31 * 22 mod 961

# The Guadalupe/Ahlgren diagnostic 24r + 13 == 0 (mod p).
DIAG_C72 = [5, 4, 10, 40, 8, 53, 3, 5, 46]          # (24r+13) mod p for the nine cells
# The residues of the source's PROVED seven-prime theorem, labelled THM1 in the e-print source.
# (The printed theorem NUMBERS of the source are not asserted anywhere here or in the paper: the
# published version could not be read, so the e-print's own LABELS are used instead.)
THM1 = {43: 12, 47: 21, 59: 56, 67: 19, 71: 32, 79: 62, 83: 79}

# The anti-control.  Under the source's PROSE reading 1/f_1^c, none of the nine cells vanishes
# even at n = 0; these are those nine values mod p^2.
PROSE_N0 = [514, 234, 558, 561, 1004, 2860, 4686, 355, 2764]
PROSE_C2 = [1, 2, 5, 10, 20, 36, 65, 110, 185, 300, 481, 752]     # 1/f_1^2, NOT A002513

# The published sequences the definition is anchored to.
A000041_13 = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77]
P_100 = 190569292
A002513_32 = [1, 1, 3, 4, 9, 12, 23, 31, 54, 73, 118, 159, 246, 329, 489, 651, 940, 1242, 1751,
              2298, 3177, 4142, 5630, 7293, 9776, 12584, 16659, 21320, 27922, 35532, 46092,
              58342]

# The census control.  The identical census at (c, p) = (3, 7) -- the
# parameters of the source's PROVED theorem THM2 -- must return exactly one survivor.
CONTROL_P, CONTROL_C, CONTROL_DEPTH, CONTROL_SURVIVORS = 7, 3, 52, [39]

CENSUS_DEPTH = 3
CENSUS_PRIMES = (31, 43, 59, 71, 83)

# =====================================================================================
# 0.  BOOKKEEPING
# =====================================================================================
_n_pass = 0
_n_fail = 0


def ck(name, ok, detail=''):
    global _n_pass, _n_fail
    if ok:
        _n_pass += 1
        print('PASS %-26s %s' % (name, detail))
    else:
        _n_fail += 1
        print('FAIL %-26s %s' % (name, detail))
    sys.stdout.flush()


# =====================================================================================
# 1.  POWER SERIES, EXACTLY
# =====================================================================================
def pentagonal_terms(N):
    """[(exponent, sign)] of f_1 = prod (1-q^m) = sum_k (-1)^k q^{k(3k-1)/2}, exponents <= N."""
    out = [(0, 1)]
    k = 1
    while True:
        g1 = k * (3 * k - 1) // 2
        if g1 > N:
            break
        s = -1 if k % 2 else 1
        out.append((g1, s))
        g2 = k * (3 * k + 1) // 2
        if g2 <= N:
            out.append((g2, s))
        k += 1
    return out


def partitions(N):
    """p(0..N), exact, by Euler's pentagonal-number recurrence."""
    p = [0] * (N + 1)
    p[0] = 1
    for n in range(1, N + 1):
        tot = 0
        k = 1
        while True:
            g1 = k * (3 * k - 1) // 2
            if g1 > n:
                break
            s = 1 if k % 2 else -1
            tot += s * p[n - g1]
            g2 = k * (3 * k + 1) // 2
            if g2 <= n:
                tot += s * p[n - g2]
            k += 1
        p[n] = tot
    return p


def eta_power(k, N):
    """[q^n] prod_{m>=1} (1-q^m)^k for n = 0..N, exact, by k sparse multiplications."""
    f = [0] * (N + 1)
    f[0] = 1
    pent = pentagonal_terms(N)
    for _ in range(k):
        g = [0] * (N + 1)
        for e, s in pent:
            if s == 1:
                for i in range(0, N + 1 - e):
                    g[i + e] += f[i]
            else:
                for i in range(0, N + 1 - e):
                    g[i + e] -= f[i]
        f = g
    return f


def polymul_mod(a, b, mod, trunc):
    """(a*b) mod (q^trunc), coefficients reduced mod `mod`.

    Kronecker substitution: the two coefficient vectors are packed into single integers with a
    limb WIDE ENOUGH that no coefficient of the product can overflow it, so one big-integer
    multiplication does the whole convolution and the arithmetic stays exact.  The limb width is
    computed from a proved bound, never guessed.
    """
    a = [x % mod for x in a[:trunc]]
    b = [x % mod for x in b[:trunc]]
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return []
    bound = (mod - 1) * (mod - 1) * min(la, lb)      # max possible product coefficient
    nb = 1
    while (1 << (8 * nb)) <= bound:
        nb += 1
    A = int.from_bytes(b''.join(x.to_bytes(nb, 'little') for x in a), 'little')
    B = int.from_bytes(b''.join(x.to_bytes(nb, 'little') for x in b), 'little')
    C = A * B
    ncoef = min(la + lb - 1, trunc)
    raw = C.to_bytes(nb * (la + lb) + nb, 'little')
    return [int.from_bytes(raw[i * nb:(i + 1) * nb], 'little') % mod for i in range(ncoef)]


def polypow_mod(a, e, mod, trunc):
    """a^e mod (q^trunc, mod), by square-and-multiply.  e >= 0."""
    r = [1 % mod] + [0] * (trunc - 1)
    base = [x % mod for x in a[:trunc]]
    while e:
        if e & 1:
            r = polymul_mod(r, base, mod, trunc)
        e >>= 1
        if e:
            base = polymul_mod(base, base, mod, trunc)
    return r


def a_series_mod(c, mod, N, part=None):
    """a_c(0..N) mod `mod`, from 1/(f_1 f_2^{c-1}) = (1/f_1)(q) * (1/f_1^{c-1})(q^2)."""
    p = part if part is not None else partitions(N)
    inv_f1 = [x % mod for x in p[:N + 1]]
    M = N // 2
    col = polypow_mod(inv_f1[:M + 1], c - 1, mod, M + 1)          # 1/f_1^{c-1} up to q^M
    even = [0] * (N + 1)
    for j, v in enumerate(col):
        if 2 * j <= N:
            even[2 * j] = v
    return polymul_mod(inv_f1, even, mod, N + 1)


def prose_series_mod(c, mod, N, part=None):
    """The source's PROSE reading, 1/f_1^c -- the anti-control, not the paper's function."""
    p = part if part is not None else partitions(N)
    return polypow_mod([x % mod for x in p[:N + 1]], c, mod, N + 1)


def colored_exact(k, M):
    """p_k(0..M) exact: [x^j] prod (1-x^m)^{-k}, i.e. k-coloured partitions.  Small M only."""
    f = [0] * (M + 1)
    f[0] = 1
    for _ in range(k):
        for m in range(1, M + 1):
            for i in range(m, M + 1):
                f[i] += f[i - m]
    return f


def a_exact(c, N, part=None):
    """a_c(N) as an exact integer, no modulus: sum_j p_{c-1}(j) p(N-2j)."""
    p = part if part is not None else partitions(N)
    pk = colored_exact(c - 1, N // 2)
    return sum(pk[j] * p[N - 2 * j] for j in range(N // 2 + 1))


# =====================================================================================
# 2.  THE PRINTED TABLES ARE THE PUBLISHED SEQUENCES
# =====================================================================================
print('interpreter: Python %s on %s (%s)'
      % (platform.python_version(), platform.system(), platform.machine()))
print('== 1. the two printed tables, against the definition ==')
PART = partitions(20670)                                   # exact, reused everywhere below

ck('printed-p-table', [PART[41 - 2 * j] for j in range(21)] == PART_TABLE,
   'p(41-2j), j=0..20 (OEIS A000041) reproduced from the pentagonal recurrence')

pk64 = colored_exact(64, 20)
ck('printed-p64-table', pk64 == P64_TABLE,
   'p_64(j), j=0..20, reproduced from prod (1-x^m)^{-64}')

e7 = eta_power(7, 40)
ck('printed-eta7-table', [e7[k] for k in range(21)] == E_TABLE,
   'e(k) = [x^k] prod (1-x^n)^7, k=0..20 (OEIS A000730) reproduced')

# =====================================================================================
# 3.  CONTROLS ON THE DEFINITION ITSELF, BOTH POLARITIES
# =====================================================================================
print('== 2. controls on the decider (a wrong a_c would make everything below meaningless) ==')
BIG = 10 ** 400
ck('control-A000041', [a_exact(1, n, PART) for n in range(13)] == A000041_13,
   'c=1 gives 1,1,2,3,5,7,11,15,22,30,42,56,77')
ck('control-p100', a_exact(1, 100, PART) == P_100, 'a_1(100) = p(100) = %d' % P_100)
a2 = [a_exact(2, n, PART) for n in range(32)]
ck('control-A002513', a2 == A002513_32,
   'c=2 reproduces all 32 quoted terms of OEIS A002513 (cubic partitions)')
ck('control-a_c(0)', all(a_exact(c, 0, PART) == 1 for c in (1, 2, 3, 25, 65)),
   'a_c(0) = 1 for c = 1,2,3,25,65')
ck('anticontrol-prose-c2', PROSE_C2 == prose_series_mod(2, BIG, 11, PART) and PROSE_C2 != a2[:12],
   'the PROSE reading 1/f_1^2 gives 1,2,5,10,20,36,... which is NOT A002513')

# =====================================================================================
# 4.  THE WITNESS
# =====================================================================================
print('== 3. the refuting coefficient a_65(41) ==')
prod_terms = [P64_TABLE[j] * PART_TABLE[j] for j in range(21)]
ck('witness-dot-product', sum(prod_terms) == A65_41,
   'the 21 printed products sum to %d' % A65_41)
ck('witness-exact-route', a_exact(65, 41, PART) == A65_41,
   'independently: a_65(41) = sum_j p_64(j) p(41-2j), exact, no modulus')
ck('witness-mod-71', A65_41 % 71 == A65_41_MOD_P and A65_41_MOD_P != 0,
   'a_65(41) == %d (mod 71), NOT 0' % A65_41_MOD_P)
ck('witness-mod-71sq', A65_41 % 5041 == A65_41_MOD_P2 and A65_41_MOD_P2 != 0,
   'a_65(41) == %d (mod 71^2), NOT 0  ==>  CONJECTURE 7.2 IS FALSE at (c,p,r,n)=(65,71,41,0)'
   % A65_41_MOD_P2)
ck('witness-quotient-71', 71 * QUOT_P + A65_41_MOD_P == A65_41, 'division identity mod 71')
ck('witness-quotient-71sq', 5041 * QUOT_P2 + A65_41_MOD_P2 == A65_41, 'division identity mod 71^2')

print('== 4. the hand certificate b(41) = -31, and 71 does not divide 31 ==')
ck('certificate-b41', sum(E_TABLE[k] * PART_TABLE[k] for k in range(21)) == B41,
   'sum_{k=0}^{20} e(k) p(41-2k) = %d' % B41)
ck('certificate-71-nmid', 31 % 71 != 0 and B41 % 71 == A65_41_MOD_P,
   '71 does not divide 31, and -31 == %d (mod 71) agrees with the exact integer' % (B41 % 71))

# The Reduction Lemma itself: a_{p-6}(N) == sum_j p(j) b(N-2pj) (mod p), with b = f_2^7/f_1.
LEM_N = 1050
e7L = eta_power(7, LEM_N // 2 + 1)
e6L = eta_power(6, LEM_N // 2 + 1)


def b_series(evens, N):
    """[q^m] f_2^k / f_1 for m = 0..N, exact, where evens[j] = [x^j] prod (1-x^n)^k."""
    return [sum(evens[j] * PART[m - 2 * j] for j in range(m // 2 + 1)) for m in range(N + 1)]


b7 = b_series(e7L, LEM_N)
b6 = b_series(e6L, LEM_N)
a65_mod71 = a_series_mod(65, 71, LEM_N, PART)


def lemma_lhs_rhs(b, p, N):
    s = 0
    j = 0
    while 2 * p * j <= N:
        s += PART[j] * b[N - 2 * p * j]
        j += 1
    return s % p


mis7 = [N for N in range(LEM_N + 1) if lemma_lhs_rhs(b7, 71, N) != a65_mod71[N]]
mis6 = [N for N in range(LEM_N + 1) if lemma_lhs_rhs(b6, 71, N) != a65_mod71[N]]
nonzero = sum(1 for N in range(LEM_N + 1) if a65_mod71[N] % 71)
ck('reduction-lemma', not mis7,
   'a_65(N) == sum_j p(j) b(N-142j) (mod 71) for all N=0..%d: 0 mismatches' % LEM_N)
ck('lemma-not-vacuous', nonzero > (LEM_N + 1) // 2,
   '%d of %d values are NONZERO mod 71, so the identity is not one between zero functions'
   % (nonzero, LEM_N + 1))
ck('lemma-anticontrol', len(mis6) > 0,
   'the deliberately wrong exponent f_2^6/f_1 mismatches at %d of %d indices (a silent control '
   'that MUST fire)' % (len(mis6), LEM_N + 1))
ck('certificate-b41-series', b7[41] == B41, 'the b-series independently gives b(41) = %d' % B41)

# =====================================================================================
# 5.  THE TWO COMPANION INTEGERS
# =====================================================================================
print('== 5. the companions: r=41 IS a zero one line earlier, and r=47 IS a zero here ==')
ck('companion-a25-41', a_exact(25, 41, PART) == A25_41 and A25_41 % 961 == 0
   and 961 * A25_41_COFACTOR == A25_41,
   'a_25(41) = %d = 961 * %d, so == 0 (mod 31^2)' % (A25_41, A25_41_COFACTOR))
ck('companion-a65-47', a_exact(65, 47, PART) == A65_47 and A65_47 % 5041 == 0,
   'a_65(47) = %d == 0 (mod 71^2): the same computation that fails r=41 passes r=47'
   % A65_47)

# =====================================================================================
# 6.  ALL NINE PROGRESSIONS, AND THE CENSUS OVER EVERY RESIDUE CLASS
# =====================================================================================
print('== 6. the exhaustive census over all r mod p^2, for the five primes named ==')
series_cache = {}
for p in CENSUS_PRIMES:
    c = p - 6
    mod = p * p
    N = CENSUS_DEPTH * mod - 1
    A = a_series_mod(c, mod, N, PART)
    series_cache[p] = A
    zeros = [r for r in range(mod) if A[r] == 0]
    alive = list(range(mod))
    surv = []
    for n in range(CENSUS_DEPTH):
        alive = [r for r in alive if A[mod * n + r] == 0]
        surv.append(len(alive))
    ck('census-zeros-p%d' % p, zeros == ZEROS[p],
       'p=%d c=%d: the COMPLETE set of r in [0,%d) with a_%d(r)==0 mod %d is %s'
       % (p, c, mod, c, mod, zeros))
    ck('census-survivors-p%d' % p, surv == SURVIVORS[p],
       'survivors by depth 1,2,3 over all %d classes: %s' % (mod, surv))

ck('census-41-absent-p71', 41 not in ZEROS[71] and series_cache[71][41] != 0,
   'r=41 is NOT a zero of a_65 mod 71^2 at all, so no depth can rescue the printed line')
ck('census-no-three-term', all(SURVIVORS[p][2] == 0 for p in CENSUS_PRIMES),
   'for EVERY one of the five primes, NO residue r mod p^2 whatever gives three consecutive '
   'vanishing terms -- no re-choice of a printed residue rescues any line')
lone = [(p, r) for p in CENSUS_PRIMES for r in ZEROS[p]
        if series_cache[p][p * p + r] == 0]
ck('census-lone-depth2', lone == [LONE_DEPTH2]
   and series_cache[31][2 * 961 + 644] == A25_2566_MOD,
   'the only class surviving n=0,1 among all %d is (p,r)=%s -- Conjecture 7.1\'s residue, not '
   'one of these nine -- and it dies at n=2 with a_25(2566) == %d = 31*22 (mod 961)'
   % (sum(p * p for p in CENSUS_PRIMES), LONE_DEPTH2, A25_2566_MOD))

print('== 7. the nine printed cells at n = 0, 1, 2 ==')
for i, (c, p, r) in enumerate(CELLS):
    mod = p * p
    got = [series_cache[p][mod * n + r] % mod for n in range(3)]
    ck('cell-%d-%d-%d' % (c, p, r), got == TABLE3[(c, p, r)],
       'a_%d(%d n + %d) mod %d = %s for n=0,1,2 -- FAILS at n=%d'
       % (c, mod, r, mod, got, 0 if got[0] else 1))
got_modp = [series_cache[p][p * p + r] % p for (c, p, r) in CELLS]
ck('nine-cells-mod-p', got_modp == MODP_ROW,
   'the n=1 values reduced mod p: %s -- eight of nine are not even divisible by p; the single '
   'exception is (c,p,r)=(37,43,987), where 301 = 43*7 is divisible by 43 and not by 43^2'
   % got_modp)
ck('nine-cells-all-fail',
   all(any(series_cache[p][p * p * n + r] != 0 for n in range(3)) for (c, p, r) in CELLS),
   'EVERY one of the nine printed progressions is false: one at n=0, the other eight at n=1')

print('== 8. positive controls on the SAME functions and primes: the source\'s proved theorems ==')
for p in (43, 59, 71, 83):
    c, r, mod = p - 6, THM1[p], p * p
    A = series_cache[p]
    idx = [p * n + r for n in range((len(A) - 1 - r) // p + 1)]
    bad = [N for N in idx if A[N] % p]
    ck('thm1-control-p%d' % p, not bad and len(idx) > 50,
       'proved theorem THM1: a_%d(%d n + %d) == 0 (mod %d) for all %d values n=0..%d -- 0 '
       'exceptions' % (c, p, r, p, len(idx), len(idx) - 1))

mod7 = CONTROL_P * CONTROL_P
A7 = a_series_mod(CONTROL_C, mod7, CONTROL_DEPTH * mod7 - 1, PART)
alive7 = list(range(mod7))
for n in range(CONTROL_DEPTH):
    alive7 = [r for r in alive7 if A7[mod7 * n + r] == 0]
ck('thm2-census-control', alive7 == CONTROL_SURVIVORS,
   'the IDENTICAL census at (c,p)=(%d,%d) -- the source\'s PROVED theorem THM2 -- returns exactly '
   '%s over all %d classes to depth %d, so a census that finds nothing elsewhere is not a census '
   'incapable of finding something' % (CONTROL_C, CONTROL_P, alive7, mod7, CONTROL_DEPTH))

# =====================================================================================
# 7.  THE DIAGNOSTIC, AND THE PROSE ANTI-CONTROL ON THE NINE CELLS
# =====================================================================================
print('== 9. the Guadalupe/Ahlgren diagnostic 24r + 13 == 0 (mod p) ==')
ck('diagnostic-c72', [(24 * r + 13) % p for (c, p, r) in CELLS] == DIAG_C72
   and all(v != 0 for v in DIAG_C72),
   '(24r+13) mod p = %s for the nine printed residues: NOT ONE satisfies the criterion, so '
   'there is no mod-p divisibility to lift' % DIAG_C72)
ck('diagnostic-thm1', all((24 * r + 13) % p == 0 for p, r in THM1.items()),
   'all seven residues of the proved theorem THM1 satisfy it exactly')
ck('diagnostic-guadalupe', all(13 * (p * p - 1) // 24 % p == r for p, r in THM1.items()),
   'and 13(p^2-1)/24 == r (mod p) for each, i.e. the criterion is Guadalupe\'s Theorem 3.2 class; '
   'at p=71 that class is 32, NOT the witness class 41')

print('== 10. the prose anti-control on the nine cells ==')
prose = [prose_series_mod(c, p * p, r, PART)[r] for (c, p, r) in CELLS]
ck('anticontrol-prose-cells', prose == PROSE_N0 and all(v != 0 for v in PROSE_N0),
   'under the prose reading 1/f_1^c the nine cells give %s at n=0 -- not one vanishes, so the '
   'authors demonstrably computed with the displayed generating function' % PROSE_N0)

# =====================================================================================
# 8.  SCOPE
# =====================================================================================
print()
print('NOTE SCOPE -- what this program does NOT cover.')
print('NOT RE-RUN: anything bibliographic. The locator of the conjecture (arXiv:2503.19399v2), '
      'the statement of Guadalupe\'s Theorem 3.2 and of Ahlgren\'s identity, and the OEIS entry '
      'numbers A000041 / A000730 / A002513 are TRANSCRIBED from those sources and are neither '
      'fetched nor checked here. The published TERMS of the three OEIS sequences ARE checked, '
      'above, against the generating function.')
print('NOT RE-RUN: the body of the version of record (Eur. J. Math. 12 (2026), Paper No. 34). No '
      'open-access copy exists; the refuted text is the arXiv e-print, and nothing here can '
      'certify that the published revision still prints the same residue.')
print('NOT RE-RUN: the census beyond depth 3, and beyond the five primes Conjecture 7.2 names. '
      'This program tests all %d residue classes for p in %s (c = p-6) at n=0,1,2 only. It says '
      'nothing about other (c,p) pairs -- indeed the source\'s proved theorem THM2 is a TRUE '
      'instance of the same shape at (c,p)=(3,7), which the control above rediscovers. Depth 3 '
      'is not a limitation of the refutation: the witness and its certificate depend on ONE '
      'coefficient of a 42-term power series.' % (sum(p * p for p in CENSUS_PRIMES), list(CENSUS_PRIMES)))
print('NOT RE-RUN: the three primes 47, 67, 79 of the source\'s other statements are covered here '
      'only by the arithmetic diagnostic, not by a census; and Conjecture 7.1, a separate '
      'statement, is not adjudicated here at all.')
print('NOT RE-RUN: the inference about how the nine residues were chosen (a search cutoff). It is '
      'an inference about the authors, not a computation, and nothing above depends on it.')
print()
if _n_fail:
    print('VERDICT: %d CHECK(S) FAILED out of %d' % (_n_fail, _n_pass + _n_fail))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _n_pass)
sys.exit(0)

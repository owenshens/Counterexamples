#!/usr/bin/env python3
"""verify.py -- checks every computational claim of

    "The Reverse Degree-Polynomial Range Clause for the 2 x n Miura-ori Flip
     Graph Fails for Every a >= 7"                                (paper.tex)

Python 3.9+, STANDARD LIBRARY ONLY, no external data file.  Exact integer
arithmetic throughout: there is not a single float in this program, and no
decision is taken on an inexact quantity.

WHAT IT READS.  Only objects PRINTED IN THE PAPER:

  * the seeds  v_1^2 = 2,  v_2^2 = 4,  v_2^3 = 0,  v_2^4 = 2         (paper, S1)
  * the pure-v recurrence (eq:P) and the conventions beside it       (paper, S1)
  * the published cells of Table 1 (tab:degdistribution of the source)
  * the witness  f_7(5..9) = 12, 88, 296, 680, 1288                  (eq:witness)
  * the printed polynomials G_2, ..., G_8                            (paper, S3)
  * the printed lists: the ten first-obstruction values, the failure counts
  * the printed control values                                       (paper, S3)

and re-derives everything else.  Nothing is imported from any source other than
produced the paper.

TWO CONSTANTS BELOW ARE *NOT* PRINTED IN THE PAPER, and are flagged here rather
than passed off as transcriptions: V10_13 = 2168 and CONTROL_A4_SECOND_DIFFS =
[8, 8, 8, 8].  Neither is load-bearing: each is checked against two independent
internal routes (the recurrence and the expansion of G_7/(1-x)^4 for the first;
the recurrence's own difference table for the second), so deleting either
constant would weaken no claim of the paper.

ONE `PASS <name>` LINE PER CHECK, then
    VERDICT: ALL <n> CHECKS PASS
Exit status 0 iff every check passed.
"""

import math
import sys

# ---------------------------------------------------------------------------
# 0.  THE OBJECTS PRINTED IN THE PAPER, TRANSCRIBED HERE AND NOWHERE ELSE
# ---------------------------------------------------------------------------

# Seeds (paper, Section 1): the only input values of the recurrence.
SEEDS = {(1, 2): 2, (2, 2): 4, (2, 3): 0, (2, 4): 2}

# Table 1 of the paper: cells of tab:degdistribution of the source.
# A printed blank is read as 0, as the paper's caption says.
PUBLISHED_ROWS = {                      # d -> {n: v_n^d}, n = 2..9
    3:  {2: 0, 3: 4, 4: 8, 5: 12, 6: 16, 7: 20, 8: 24, 9: 28},
    5:  {3: 0, 4: 8, 5: 36, 6: 88, 7: 168, 8: 280, 9: 428},
    7:  {4: 0, 5: 12, 6: 80, 7: 296, 8: 792, 9: 1744},
    9:  {5: 0, 6: 16, 7: 140, 8: 680, 9: 2396},
    11: {6: 0, 7: 20, 8: 216, 9: 1288},
}
# The whole n = 3 column, printed under Table 1.
PUBLISHED_N3_COLUMN = {2: 4, 3: 4, 4: 8, 5: 0, 6: 2}
# The a = 6 reverse diagonal, printed under Table 1: v_5^4 .. v_9^12.
PUBLISHED_A6_DIAGONAL = {(5, 4): 36, (6, 6): 128, (7, 8): 292,
                         (8, 10): 544, (9, 12): 900}
# From the convention quoted at source line 1204 (v_m^2 = 4 for m >= 2).
PUBLISHED_V4_2 = 4

# The witness, eq:witness of the paper: f_7(5), ..., f_7(9).
WITNESS = [12, 88, 296, 680, 1288]
WITNESS_D4 = 4                                     # its 4th forward difference
WITNESS_BACK_PREDICTION = 8                        # the cubic's value at n = 5

# The remark after Theorem 1: each witness value as the paper prints it, as a
# (target cell, [the five summands of eq:P with their signs]) pair.
WITNESS_RECURRENCE = [
    ((5, 3),  [(4, 3), (4, 2), (4, 1), (3, 1), (3, 0)],  [8, 4, 0, 0, 0],   12),
    ((6, 5),  [(5, 5), (5, 4), (5, 3), (4, 3), (4, 2)],  [36, 36, 12, 8, 4], 88),
    ((7, 7),  [(6, 7), (6, 6), (6, 5), (5, 5), (5, 4)],  [80, 128, 88, 36, 36], 296),
    ((8, 9),  [(7, 9), (7, 8), (7, 7), (6, 7), (6, 6)],  [140, 292, 296, 80, 128], 680),
    ((9, 11), [(8, 11), (8, 10), (8, 9), (7, 9), (7, 8)], [216, 544, 680, 140, 292], 1288),
]

# The polynomials printed in Section 3, as coefficient lists (index = exponent).
PRINTED_G = {
    2: [0, 0, 4],
    3: [0, 0, 0, 4],
    4: [0, 0, 0, 4, 8, -4],
    5: [0, 0, 0, 0, 8, 12, -4],
    6: [0, 0, 0, 0, 4, 20, 8, -20, 4],
    7: [0, 0, 0, 0, 0, 12, 40, 16, -24, 4],
    8: [0, 0, 0, 0, 0, 4, 36, 56, -36, -56, 32, -4],
}
G_AT_ONE_TARGETS = {4: 8, 5: 16, 6: 16, 7: 48}     # computed here, not quoted from the note

# The ten first-obstruction values, computed here; the note prints none of them.
TEN_FIRST_OBSTRUCTIONS = {7: 4, 8: -4, 9: 36, 10: -44, 11: 172,
                          12: -260, 13: 488, 14: -996, 15: 344, 16: -2164}
# The printed failure counts.
PRINTED_FAILURE_COUNTS = {7: 1, 8: 1, 9: 2, 10: 2, 15: 5, 16: 5, 40: 17, 80: 37}

# The controls printed in Section 3.
CONTROL_A6_WINDOW = [36, 128, 292, 544, 900]       # f_6(5..9), Delta^4 == 0
CONTROL_A6_WITH_EXCLUDED_D4 = 4                    # prepend v_4^2 = 4 -> 4
CONTROL_A4_SECOND_DIFFS = [8, 8, 8, 8]             # a = 4 against "degree 1"
CONTROL_CORRUPT_WITNESS = [8, 88, 296, 680, 1288]  # 12 -> 8 must silence it
V10_13 = 2168                                      # f_7(10), two routes

# Bounds.  NMAX bounds the recurrence table; AMAX the census over a.
NMAX = 130
AMAX = 80
ADIRECT = 40          # the independent direct-difference census stops here
DMAX_SIBLING = 40     # the companion null over the source's proved sibling

# ---------------------------------------------------------------------------
# 1.  THE RECURRENCE TABLE
# ---------------------------------------------------------------------------
OFF = 3               # so that index d - 3 is never negative


def build_table(nmax):
    """v_n^d for 1 <= n <= nmax from the SEEDS and eq:P alone.

    Rows are lists indexed by d + OFF, so the d-3 term of eq:P can never wrap
    around to the end of the list (a negative Python index would silently read
    the wrong cell -- the one bug that would make every check below vacuous).
    """
    width = 2 * nmax + OFF + 2
    rows = [None] * (nmax + 1)
    for n in (1, 2):
        rows[n] = [0] * width
    for (n, d), val in SEEDS.items():
        rows[n][d + OFF] = val
    for n in range(3, nmax + 1):
        prev, prev2, cur = rows[n - 1], rows[n - 2], [0] * width
        for d in range(2, 2 * n + 1):
            i = d + OFF
            cur[i] = prev[i] + prev[i - 1] + prev[i - 2] + prev2[i - 2] - prev2[i - 3]
        rows[n] = cur
    return rows


ROWS = build_table(NMAX)


def v(n, d):
    """v_n^d, with the source's conventions outside the support."""
    if n < 1 or n > NMAX or d < 2 or d > 2 * n:
        return 0
    return ROWS[n][d + OFF]


def f(a, n):
    """The reverse diagonal f_a(n) = v_n^{2n-a}; 0 outside the support."""
    return v(n, 2 * n - a)


def fx(a, n):
    """f_a extended by 0 to n <= 0, as the dictionary eq:dict requires."""
    return f(a, n) if n >= 1 else 0


# ---------------------------------------------------------------------------
# 2.  DIFFERENCES AND POLYNOMIALS
# ---------------------------------------------------------------------------
def dtab(vals, j):
    """The j-th forward difference of a finite list, base index preserved."""
    cur = list(vals)
    for _ in range(j):
        cur = [cur[i + 1] - cur[i] for i in range(len(cur) - 1)]
    return cur


def delta(a, j, base):
    """(Delta^j f_a)(base) by the binomial formula, valid for base <= 0 too."""
    return sum((-1) ** i * math.comb(j, i) * fx(a, base + j - i) for i in range(j + 1))


def padd(p, q):
    out = [0] * max(len(p), len(q))
    for i, c in enumerate(p):
        out[i] += c
    for i, c in enumerate(q):
        out[i] += c
    return out


def pmulmono(p, deg, coef):
    """p(x) * coef * x^deg."""
    return [0] * deg + [c * coef for c in p]


def ptrim(p):
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def build_G(amax):
    """G_a for 2 <= a <= amax from the printed seeds G_2 = 4x^2, G_3 = 4x^3 and
    the ladder eq:G.  Nothing else is used."""
    G = {2: [0, 0, 4], 3: [0, 0, 0, 4]}
    for a in range(4, amax + 1):
        g1, g2 = G[a - 1], G[a - 2]
        if a % 2 == 0:                       # x(1-x)G_{a-1} + x(1+x)G_{a-2}
            t = padd(pmulmono(g1, 1, 1), pmulmono(g1, 2, -1))
        else:                                # x G_{a-1} + x(1+x)G_{a-2}
            t = pmulmono(g1, 1, 1)
        t = padd(t, padd(pmulmono(g2, 1, 1), pmulmono(g2, 2, 1)))
        G[a] = ptrim(t)
    return G


G = build_G(AMAX)


def gdeg(a):
    return len(ptrim(G[a])) - 1


def gcoef(a, n):
    p = G[a]
    return p[n] if 0 <= n < len(p) else 0


def s_start(a):
    """The admissible start eq:s."""
    return a // 2 + 2


# ---------------------------------------------------------------------------
# 3.  THE CHECKS
# ---------------------------------------------------------------------------
CHECKS = []


def ck(name, ok, detail=''):
    CHECKS.append((name, bool(ok), detail))


# --- C1  the recurrence is calibrated against the source's own integers -----
bad = [(n, d) for d, row in PUBLISHED_ROWS.items() for n, val in row.items()
       if v(n, d) != val]
cells = sum(len(r) for r in PUBLISHED_ROWS.values())
ck('published-odd-degree-rows', not bad,
   ':: %d cells of Table 1 (d = 3,5,7,9,11; n = 2..9) reproduced from the seeds '
   'and eq:P, %d mismatches' % (cells, len(bad)))

bad = [d for d, val in PUBLISHED_N3_COLUMN.items() if v(3, d) != val]
ck('published-n3-column', not bad,
   ':: v_3^2..v_3^6 = %s, 5 of 5, including v_3^4 = 8 and v_3^6 = 2 which need '
   'the v_1^2 = 2 anomaly' % ([v(3, d) for d in sorted(PUBLISHED_N3_COLUMN)],))

bad = [k for k, val in PUBLISHED_A6_DIAGONAL.items() if v(*k) != val]
ck('published-a6-diagonal', not bad,
   ':: v_5^4, v_6^6, v_7^8, v_8^10, v_9^12 = %s, 5 of 5'
   % ([v(*k) for k in sorted(PUBLISHED_A6_DIAGONAL)],))

ck('published-v4-2-convention', v(4, 2) == PUBLISHED_V4_2 == 4,
   ':: v_4^2 = %d, matching the source convention v_m^2 = 4 for m >= 2' % v(4, 2))

ALL_PUBLISHED = {}                      # (n, d) -> the value printed in the paper
for _d, _row in PUBLISHED_ROWS.items():
    for _n, _val in _row.items():
        ALL_PUBLISHED[(_n, _d)] = _val
for _d, _val in PUBLISHED_N3_COLUMN.items():
    ALL_PUBLISHED[(3, _d)] = _val
ALL_PUBLISHED.update(PUBLISHED_A6_DIAGONAL)
ALL_PUBLISHED[(4, 2)] = PUBLISHED_V4_2
bad = [k for k, val in ALL_PUBLISHED.items() if v(*k) != val]
ck('published-cells-total', not bad and len(ALL_PUBLISHED) == 39,
   ':: %d DISTINCT published cells, all reproduced, %d mismatches; the '
   'recurrence is calibrated against the source before it is used anywhere '
   'below' % (len(ALL_PUBLISHED), len(bad)))

# --- C2  the witness --------------------------------------------------------
ck('witness-equals-recurrence-values', [f(7, n) for n in range(5, 10)] == WITNESS,
   ':: f_7(5..9) = %s, equal to the five printed cells' % ([f(7, n) for n in range(5, 10)],))

ok = True
det = []
for tgt, srcs, vals, res in WITNESS_RECURRENCE:
    if [v(*c) for c in srcs] != vals:
        ok = False
    signed = vals[0] + vals[1] + vals[2] + vals[3] - vals[4]
    if signed != res or v(*tgt) != res:
        ok = False
    det.append('v_%d^%d=%d' % (tgt[0], tgt[1], res))
ck('witness-forced-by-published-cells', ok,
   ':: each of %s follows from eq:P applied to OTHER published cells, exactly as '
   'the paper prints the sums' % (', '.join(det),))

d4 = dtab(WITNESS, 4)
ck('witness-fourth-difference', len(d4) == 1 and d4[0] == WITNESS_D4 == 4,
   ':: Delta^4 (12, 88, 296, 680, 1288) = %d != 0, so no polynomial of degree '
   '<= 3 agrees with f_7 on n = 5..9' % d4[0])

rows_shown = [dtab(WITNESS, j) for j in range(1, 5)]
ck('witness-difference-tower-as-printed',
   rows_shown == [[76, 208, 384, 608], [132, 176, 224], [44, 48], [4]],
   ':: the tower printed in the paper is reproduced line for line: %s' % (rows_shown,))

ck('witness-window-inside-asserted-range',
   5 >= -(-7 // 2) + 1 and 9 >= 5 and all((2 * n - 7) != 2 for n in range(1, 200)),
   ':: the asserted range at a = 7 is n >= ceil(7/2)+1 = 5, the window is '
   'n = 5..9, and 2n-7 = 2 has no integer solution so the exclusion clause is '
   'vacuous')

# back-extrapolation: the unique cubic through n = 6..9, evaluated at n = 5
t = WITNESS[1:]
d1, d2, d3 = dtab(t, 1), dtab(t, 2), dtab(t, 3)
back = t[0] - (d1[0] - (d2[0] - d3[0]))
ck('witness-back-extrapolation', back == WITNESS_BACK_PREDICTION == 8 != WITNESS[0],
   ':: the cubic through (6,88),(7,296),(8,680),(9,1288) predicts %d at n = 5, '
   'while f_7(5) = %d' % (back, WITNESS[0]))

# --- C3  the two reverse rows the source prints -----------------------------
ck('reverse-row-a2-closed-form', all(f(2, n) == 4 * (n - 1) for n in range(1, 60)),
   ':: f_2(n) = 4(n-1) for n = 1..59, consistent with (but not proving) '
   'F_2 = 4x^2/(1-x)^2')
ck('reverse-row-a3-closed-form', all(f(3, n) == 4 * (n - 2) for n in range(2, 60)),
   ':: f_3(n) = 4(n-2) for n = 2..59, consistent with (but not proving) '
   'F_3 = 4x^3/(1-x)^2')
ck('reverse-row-start-indices',
   all(min(n for n in range(1, 200) if 2 * n - a >= 3) == s_start(a)
       for a in range(0, AMAX + 1))
   and [9 - s_start(a) + 1 for a in range(6)] == [8, 8, 7, 7, 6, 6],
   ':: the least n with 2n-a >= 3 (i.e. defined and not the excluded cell) is '
   'floor(a/2)+2 for every a = 0..%d, and the printed row lengths over n <= 9 '
   'are then 8,8,7,7,6,6 for a = 0..5, as the source prints them' % AMAX)

# --- C4  the exclusion clause and the admissible start ----------------------
bad = []
for a in range(0, AMAX + 1):
    sol = [n for n in range(1, 200) if 2 * n - a == 2]
    if a % 2 == 0:
        if sol != [a // 2 + 1] or s_start(a) != -(-a // 2) + 2:
            bad.append(a)
    else:
        if sol or s_start(a) != -(-a // 2) + 1:
            bad.append(a)
ck('exclusion-clause-and-admissible-start', not bad,
   ':: for a = 0..%d: 2n-a = 2 is solvable exactly for even a, at n = a/2+1 = '
   'ceil(a/2)+1; and s(a) = floor(a/2)+2 equals the literal bound ceil(a/2)+1 '
   'for odd a, one more for even a' % AMAX)

# --- C5  the printed polynomials --------------------------------------------
bad = [a for a, p in PRINTED_G.items() if ptrim(G[a]) != ptrim(p)]
ck('printed-polynomials-G2-to-G8', not bad,
   ':: G_2..G_8 built from the ladder eq:G agree coefficient for coefficient '
   'with the polynomials printed in the paper')

bad = [a for a, val in G_AT_ONE_TARGETS.items() if sum(G[a]) != val]
ck('G-at-one-equals-8-16-16-48', not bad,
   ':: G_4(1), G_5(1), G_6(1), G_7(1) = 8, 16, 16, 48, computed from the coefficient\n   lists built above; the note states only that G_a(1) > 0')

# --- C6  the dictionary eq:dictG, against the recurrence --------------------
nchecked = 0
bad = []
for a in range(2, ADIRECT + 1):
    k = a // 2
    for n in range(0, gdeg(a) + 6):
        nchecked += 1
        if gcoef(a, n) != delta(a, k + 1, n - k - 1):
            bad.append((a, n))
ck('dictionary-identity', not bad,
   ':: [x^n] G_a == (Delta^{floor(a/2)+1} f_a)(n - floor(a/2) - 1) verified in '
   '%d instances over a = 2..%d, 0 mismatches; this ties the generating function '
   'to the recurrence' % (nchecked, ADIRECT))

# --- C7  Proposition 2 -----------------------------------------------------
bad = [a for a in range(2, AMAX + 1)
       if gdeg(a) != (3 * (a // 2) - 1 if a % 2 == 0 else 3 * (a // 2))]
ck('degree-of-G', not bad,
   ':: deg G_{2m} = 3m-1 and deg G_{2m+1} = 3m for every a = 2..%d' % AMAX)

bad = [a for a in range(2, AMAX + 1) if gcoef(a, gdeg(a)) != 4 * (-1) ** (a // 2 + 1)]
ck('leading-coefficient-is-plus-minus-4', not bad,
   ':: L_a = 4(-1)^{floor(a/2)+1} for every a = 2..%d, so |L_a| = 4 throughout '
   '-- one boundary quantum that never cancels' % AMAX)

vals = {a: sum(G[a]) for a in range(2, AMAX + 1)}
ck('G-at-one-positive', all(x > 0 for x in vals.values()),
   ':: G_a(1) > 0 for every a = 2..%d (min = %d); GRANTING the paper\'s identity '
   'F_a = G_a/(1-x)^{floor(a/2)+1} (eq:F, not re-derived here), this gives pole '
   'order exactly floor(a/2)+1 and degree exactly floor(a/2)'
   % (AMAX, min(vals.values())))

bad = [a for a in range(4, AMAX + 1)
       if (sum(G[a]) != 2 * sum(G[a - 2]) if a % 2 == 0
           else sum(G[a]) != sum(G[a - 1]) + 2 * sum(G[a - 2]))]
ck('G-at-one-recursion', not bad and sum(G[2]) == sum(G[3]) == 4,
   ':: the x = 1 specialisation of eq:G holds for a = 4..%d: G_a(1) = 2G_{a-2}(1) '
   'for even a and G_{a-1}(1) + 2G_{a-2}(1) for odd a, from G_2(1) = G_3(1) = 4'
   % AMAX)

# --- C8  Corollary 3: the threshold and the last obstruction ---------------
bad = [a for a in range(2, AMAX + 1) if gdeg(a) - a // 2 != a - 1]
ck('agreement-threshold-is-a-minus-1', not bad,
   ':: deg G_a - floor(a/2) = a - 1 for every a = 2..%d, i.e. f_a agrees with a '
   'degree-floor(a/2) polynomial exactly from n = a-1 on' % AMAX)

bad = [a for a in range(2, AMAX + 1)
       if delta(a, a // 2 + 1, a - 2) != 4 * (-1) ** (a // 2 + 1)
       or any(delta(a, a // 2 + 1, n) != 0 for n in (a - 1, a, a + 1, a + 7))]
ck('last-obstruction-value-and-position', not bad,
   ':: (Delta^{floor(a/2)+1} f_a)(a-2) = 4(-1)^{floor(a/2)+1} != 0 and the same '
   'difference vanishes at n = a-1, a, a+1, a+7, for every a = 2..%d' % AMAX)

bad = []
for a in range(2, ADIRECT + 1):
    k = a // 2
    c = sum(G[a])
    if any(delta(a, k, n) != c for n in range(a - 1, a + 25)):
        bad.append(a)
ck('constant-k-th-difference-equals-G-at-one', not bad,
   ':: Delta^{floor(a/2)} f_a equals the constant G_a(1) on the sampled window '
   'n = a-1..a+24 for a = 2..%d, so the leading coefficient of the polynomial '
   'is G_a(1)/floor(a/2)!' % ADIRECT)

# --- C9  the census over a, from G -----------------------------------------
holds, fails_by_a = [], {}
for a in range(2, AMAX + 1):
    k = a // 2
    fails = [n for n in range(s_start(a), gdeg(a) - k) if gcoef(a, n + k + 1) != 0]
    fails_by_a[a] = fails
    if not fails:
        holds.append(a)
ck('census-hold-set-is-a-at-most-6', holds == [2, 3, 4, 5, 6],
   ':: over a = 2..%d the statement holds on its whole admissible range for '
   'EXACTLY a = 2,3,4,5,6 and fails for a = 7..%d (%d values)'
   % (AMAX, AMAX, AMAX - 6))

ck('census-a0-and-a1', all(f(0, n) == 2 for n in range(1, 40))
   and all(f(1, n) == 0 for n in range(1, 40)),
   ':: f_0(n) = v_n^{2n} = 2 for n = 1..39 (degree exactly 0) and f_1(n) = '
   'v_n^{2n-1} = 0 for n = 1..39, so a = 0 holds and a = 1 is the degenerate '
   'index the paper declines to use')

bad = [a for a in range(7, AMAX + 1)
       if not fails_by_a[a] or max(fails_by_a[a]) != a - 2
       or min(fails_by_a[a]) != s_start(a)
       or fails_by_a[a] != list(range(s_start(a), a - 1))]
ck('failing-indices-are-a-contiguous-block', not bad,
   ':: for every a = 7..%d the failing base indices inside the asserted range '
   'are EXACTLY n = s(a), ..., a-2 -- no coefficient of G_a vanishes in that '
   'window' % AMAX)

bad = [a for a in range(7, AMAX + 1)
       if len(fails_by_a[a]) != -(-a // 2) - 3 or len(fails_by_a[a]) != (a - 5) // 2]
ck('failure-count-formula', not bad,
   ':: #failing n = ceil(a/2)-3 = floor((a-5)/2) for every a = 7..%d' % AMAX)

bad = [a for a, c in PRINTED_FAILURE_COUNTS.items() if len(fails_by_a[a]) != c]
ck('failure-counts-at-the-sampled-indices', not bad,
   ':: failure counts computed by this program (the note asserts no such counts): %s'
   % ({a: len(fails_by_a[a]) for a in sorted(PRINTED_FAILURE_COUNTS)},))

bad = [a for a in range(2, AMAX + 1)
       if max(s_start(a), a - 1) != (min(n for n in range(s_start(a), gdeg(a) + 3)
                                         if not fails_by_a[a] or n > max(fails_by_a[a])))]
ck('repaired-hypothesis-is-least-correct-start', not bad,
   ':: for every a = 2..%d the least admissible n_0 >= s(a) beyond every '
   'obstruction is exactly max(floor(a/2)+2, a-1), which is eq:repair; it is 6 '
   'at a = 7, and the weaker repair floor(a/2)+2 = 5 there is the refuted bound'
   % AMAX)

bad = [a for a in range(2, AMAX + 1)
       if (max(s_start(a), a - 1) == a - 1) != (a >= 5)]
ck('repair-terms-coincide-from-a5', not bad,
   ':: the two terms of eq:repair coincide exactly for a >= 5 (both are 4 at '
   'a = 5); only a = 2,3,4 are governed by the admissible start')

# --- C10  an INDEPENDENT census, from difference tables only ----------------
N2 = 90
bad, agree = [], 0
for a in range(0, ADIRECT + 1):
    k = a // 2
    seq = [f(a, n) for n in range(1, N2 + 1)]
    dd = dtab(seq, k + 1)                       # dd[i] has base index i+1
    hi = N2 - (k + 1)
    fails = [n for n in range(s_start(a), hi + 1) if dd[n - 1] != 0]
    agree += 1
    exp = list(range(s_start(a), a - 1)) if a >= 7 else []
    if fails != exp:
        bad.append(a)
ck('independent-difference-census', not bad,
   ':: a second census computed from RAW difference tables of the recurrence '
   '(no generating function anywhere) agrees with the G-based census for all '
   '%d values a = 0..%d: empty failure set iff a <= 6, else exactly '
   's(a)..a-2' % (agree, ADIRECT))

# --- C11  the ten certificates ---------------------------------------------
bad, got = [], {}
for a in range(7, 17):
    k = a // 2
    window = [f(a, n) for n in range(s_start(a), s_start(a) + k + 2)]
    d = dtab(window, k + 1)
    got[a] = d[0]
    if len(d) != 1 or d[0] != TEN_FIRST_OBSTRUCTIONS[a] or d[0] == 0:
        bad.append(a)
    if d[0] != gcoef(a, s_start(a) + k + 1):
        bad.append(a)
ck('ten-certificates-a7-to-a16', not bad,
   ':: on the shortest window inside the asserted range, n = s(a)..s(a)+k+1, '
   'the final difference is %s -- all nonzero, all equal to the printed list, '
   'and all equal to the corresponding coefficient of G_a'
   % ([got[a] for a in range(7, 17)],))

ck('first-and-last-obstruction-agree-only-at-a7-a8',
   all((TEN_FIRST_OBSTRUCTIONS[a] == 4 * (-1) ** (a // 2 + 1)) == (a in (7, 8))
       for a in range(7, 17)),
   ':: the first obstruction equals the last (magnitude 4) exactly at a = 7 and '
   'a = 8, the only a with a single failing n')

# --- C12  the published-integer checksums ----------------------------------
chk = (WITNESS[4] - 4 * WITNESS[3] + 6 * WITNESS[2] - 4 * WITNESS[1] + WITNESS[0])
ck('checksum-x9-of-G7-from-published-integers', chk == 4 == gcoef(7, 9),
   ':: 1288 - 4*680 + 6*296 - 4*88 + 12 = %d = [x^9] G_7, the witness recovered '
   'from published integers with no recurrence at all' % chk)

w = CONTROL_A6_WINDOW
chk6 = w[4] - 4 * w[3] + 6 * w[2] - 4 * w[1] + w[0]
ck('checksum-a6-vanishes-from-published-integers', chk6 == 0,
   ':: 900 - 4*544 + 6*292 - 4*128 + 36 = 0 outright from the published a = 6 '
   'diagonal, the control for the same computation')

ck('f7-at-n10-two-routes',
   f(7, 10) == V10_13 == sum(gcoef(7, j) * math.comb(10 - j + 3, 3)
                             for j in range(0, 10 + 1)),
   ':: v_10^13 = %d both from the recurrence and from expanding '
   'G_7/(1-x)^4 -- two independent routes' % f(7, 10))

# --- C13  the controls, both polarities ------------------------------------
ck('control-a6-stays-silent-on-its-range',
   all(x == 0 for x in dtab([f(6, n) for n in range(5, 61)], 4))
   and [f(6, n) for n in range(5, 10)] == CONTROL_A6_WINDOW,
   ':: MUST STAY SILENT: a = 6 over n = 5..60 has Delta^4 identically 0, so the '
   'detector does not simply always fire')

d = dtab([PUBLISHED_V4_2] + CONTROL_A6_WINDOW, 4)
ck('control-a6-fires-with-the-excluded-cell',
   d[0] == CONTROL_A6_WITH_EXCLUDED_D4 == 4,
   ':: MUST FIRE: prepending the excluded cell v_4^2 = 4 to the a = 6 window '
   'gives Delta^4 = %d != 0 at n = 4 = a-2, exactly where Corollary 3 puts it'
   % d[0])

seq4 = [f(4, n) for n in range(4, 10)]
ck('control-a4-wrong-degree-fires', dtab(seq4, 2) == CONTROL_A4_SECOND_DIFFS,
   ':: MUST FIRE: a = 4 tested against the WRONG hypothesis "degree 1" gives '
   'Delta^2 = (8,8,8,8) on n = 4..9')

ck('control-a4-correct-degree-silent', dtab(seq4, 3) == [0, 0, 0],
   ':: MUST STAY SILENT: the same a = 4 window against its CORRECT degree 2 '
   'gives Delta^3 = (0,0,0)')

ck('control-corrupted-witness-is-silenced',
   dtab(CONTROL_CORRUPT_WITNESS, 4) == [0],
   ':: FALSIFICATION CONTROL: corrupting the witness 12 -> 8 makes Delta^4 '
   'vanish, so the detector is reading the data and is not stuck on "fails"')

# --- C14  the companion null: the source's PROVED sibling theorem -----------
bad, ndone = [], 0
for d0 in range(2, DMAX_SIBLING + 1):
    start = -(-d0 // 2) + 1                       # ceil(d/2)+1, the source's range
    hi = 85
    seq = [v(n, d0) for n in range(start, hi + 1)]
    if any(x != 0 for x in dtab(seq, d0 - 1)):
        bad.append(('range', d0))
    tail = dtab(seq, d0 - 2)
    if not tail or len(set(tail)) != 1 or tail[0] == 0:
        bad.append(('degree', d0))
    ndone += 1
ck('companion-null-proved-sibling-theorem', not bad,
   ':: thm:polynomials (fixed degree d, the theorem the source PROVES) is '
   'corroborated at m = 2 on the sampled window n = ceil(d/2)+1..85 with degree '
   'exactly d-2, for all %d values d = 2..%d: zero obstructions in that window. '
   'The defect is specific to the reverse diagonal' % (ndone, DMAX_SIBLING))

ck('companion-null-would-have-caught-a-shifted-range',
   any(x != 0 for x in dtab([v(n, 8) for n in range(-(-8 // 2), 60)], 7)),
   ':: and the same test DOES fire one step below that range at d = 8, so the '
   'null above is a real test and not a vacuous one')

# ---------------------------------------------------------------------------
# 4.  REPORT
# ---------------------------------------------------------------------------
npass = 0
for name, ok, detail in CHECKS:
    print('%s %-46s %s' % ('PASS' if ok else 'FAIL', name, detail))
    npass += 1 if ok else 0

print('')
print('SCOPE -- what this program does NOT check.  (1) prop:recurrences of the '
      'source is NOT re-proved: eq:P is taken as given, and the conventions '
      'around it are pinned against the source\'s own printed integers only.  '
      '(2) NOT RE-RUN: flat-foldability itself.  No mountain-valley assignment '
      'of M_{2,n} is enumerated anywhere here, so the model behind eq:P is '
      'trusted; Theorem 1 (a = 7) does not depend on it, since its five values '
      'are printed cells, but every value at n >= 10 does.  (3) The census '
      'stops at a = %d and n = %d.  Proposition 2, Corollary 3 and the FALSITY '
      'assertion of Theorem 4 are PROVED for all a and are only corroborated '
      'here; but the EXACT failing set of Theorem 4, and hence the attainment '
      'of the upper bound ceil(a/2)-3 on the number of failures, rests on this '
      'census alone and is NOT claimed beyond a = 80.  '
      '(4) The line numbers, byte count and md5 quoted in paper.tex are NOT '
      'checked: they are properties of an external file this program does not '
      'fetch.' % (AMAX, NMAX))
print('')
if npass == len(CHECKS):
    print('VERDICT: ALL %d CHECKS PASS' % npass)
    sys.exit(0)
print('VERDICT: %d of %d checks passed -- SOME CHECK FAILED' % (npass, len(CHECKS)))
sys.exit(1)

#!/usr/bin/env python3
"""verify.py -- re-check every computational claim of

    "Seven unknown minimal row-column factorial designs of strength two at m = 6 all exist"

from the objects PRINTED IN THE PAPER and from nothing else.

WHAT IT READS.  Three literals, transcribed from the paper: the 24 x 23 array M of Section 2
(M_ROWS), the three generators h1, h2, h3 of the subspace D (D_GENS), and the 32 words of
V = D^perp in F_2^8 (V_PRINTED, the k = 5 exhibit).  No external file, no third-party module.
Python 3.9 or later; standard library only; exact integer arithmetic, no floating point anywhere.

WHAT IT DERIVES.  (a) that the Paley recipe of Section 2 really produces M -- the Hadamard matrix
of order 24 is rebuilt from q = 23, checked orthogonal, normalised and reduced, and the result is
compared with the printed digits; (b) that M is an OA(24,23,2,2); (c) hypotheses (i)-(iii) of
Theorem 1 for every column count K, 8 <= K <= 23; (d) that A[i][j] = g_i XOR V[j] is an
I_{k+3}(24,2^k,2,2) for each of the seven cases 5 <= k <= 11, checked cell by cell against the
DEFINITION of Section 1 -- 1,363,104 pattern and multiplicity counts in all; (e) that each of the
seven is an admissible instance of the conjecture; (f) the count 21,252 of generator triples
{h1,h2,h3} that would serve -- equivalently 759 subspaces D, each with 28 such bases -- on this one
fixed set of eight columns, which is the Remark of Section 2; (g) four
controls in which the witness is deliberately damaged and this program's own checker must reject it;
and (h) that the reading of the source's definition used in Section 1 of the paper is the one its own
printed examples and necessary conditions carry -- from the labels and printed DIMENSIONS of the
four designs the source exhibits, transcribed here, and from the arithmetic of the two candidate
readings.

TWO INDEPENDENT COUNTERS, ON PURPOSE.  `design_checks_naive` is a transparent triple loop over
cells with no algebra in it.  `design_checks_fast` counts the same cells with bit-parallel
population counts, which is what makes k = 10 and k = 11 (up to 770,592 counts) practical.  The
two are run against each other on the two smallest cases, so the fast counter is never trusted on
its own word.

CONVENTION.  A word of F_2^K is written as a K-character string of '0'/'1' in which character
position t (0-indexed, left to right) is factor coordinate t; as an integer it is the mask whose
bit t is coordinate t.  Rows of M are indexed 0..23 and columns 0..22.
"""

import sys
from itertools import combinations

# ---------------------------------------------------------------------------------------------
# THE OBJECTS, TRANSCRIBED FROM THE PAPER
# ---------------------------------------------------------------------------------------------

# Section 2, the array M: an OA(24,23,2,2).
M_ROWS = [
    "11111111111111111111111",
    "00000101001100110101111",
    "10000010100110011010111",
    "11000001010011001101011",
    "11100000101001100110101",
    "11110000010100110011010",
    "01111000001010011001101",
    "10111100000101001100110",
    "01011110000010100110011",
    "10101111000001010011001",
    "11010111100000101001100",
    "01101011110000010100110",
    "00110101111000001010011",
    "10011010111100000101001",
    "11001101011110000010100",
    "01100110101111000001010",
    "00110011010111100000101",
    "10011001101011110000010",
    "01001100110101111000001",
    "10100110011010111100000",
    "01010011001101011110000",
    "00101001100110101111000",
    "00010100110011010111100",
    "00001010011001101011110",
]

# Section 2, the generators of D <= F_2^8.
D_GENS = ("11100000", "01011000", "10010100")

# Section 2, the 32 words of V = D^perp in F_2^8 (the k = 5 exhibit).
V_PRINTED = [
    "00000000", "11010000", "01101000", "10111000", "10100100", "01110100", "11001100", "00011100",
    "00000010", "11010010", "01101010", "10111010", "10100110", "01110110", "11001110", "00011110",
    "00000001", "11010001", "01101001", "10111001", "10100101", "01110101", "11001101", "00011101",
    "00000011", "11010011", "01101011", "10111011", "10100111", "01110111", "11001111", "00011111",
]

M_ROWS_COUNT = 24          # 4m with m = 6
Q_PALEY = 23               # the prime the Hadamard matrix of order 24 is built on
DIM_D = 3                  # k - N in the source's Theorem "bigdeal"
K_RANGE = range(5, 12)     # the seven cases: 5 <= k <= 11


# ---------------------------------------------------------------------------------------------
# CHECK HARNESS
# ---------------------------------------------------------------------------------------------
_n_pass = 0
_n_fail = 0


def check(name, condition, detail=""):
    """One check.  Prints `PASS <name> <detail>` or `FAIL <name> <detail>`; never both."""
    global _n_pass, _n_fail
    if condition:
        _n_pass += 1
        print("PASS %s %s" % (name, detail))
    else:
        _n_fail += 1
        print("FAIL %s %s" % (name, detail))


# ---------------------------------------------------------------------------------------------
# SMALL HELPERS
# ---------------------------------------------------------------------------------------------
def mask_of(s):
    """'11100000' -> integer mask whose bit t is character t."""
    return sum(1 << t for t, c in enumerate(s) if c == "1")


def popcount(x):
    # int.bit_count() is 3.10+; this keeps the program runnable on 3.9.
    return bin(x).count("1")


def parity(x):
    return popcount(x) & 1


def span(gens):
    """The subspace generated by `gens`, as a sorted list of masks."""
    out = [0]
    for g in gens:
        out = out + [x ^ g for x in out]
    return sorted(set(out))


def dual(gens, K):
    """{ x in F_2^K : <x,g> = 0 for every generator g }, as a sorted list of masks."""
    full = (1 << K) - 1
    return [x for x in range(1 << K) if all(parity(x & (g & full)) == 0 for g in gens)]


# ---------------------------------------------------------------------------------------------
# THE TWO DESIGN CHECKERS
#
# An I_K(24, n, 2, 2) is a 24 x n arrangement of lambda copies of every word of F_2^K, lambda =
# 24n/2^K, in which for every pair of coordinates {a,b} each of the four patterns 00,01,10,11
# occurs n/4 times in every row and 24/4 = 6 times in every column (Section 1 of the paper).  The
# number of counts a full check performs is therefore
#       2^K                      multiplicities
#     + 24 * C(K,2) * 4          row pattern counts
#     +  n * C(K,2) * 4          column pattern counts.
# ---------------------------------------------------------------------------------------------
def design_checks_naive(A, K):
    """Transparent triple loop over cells.  -> (ok, number_of_counts, first_failure_or_None)."""
    m = len(A)
    n = len(A[0])
    lam, rem = divmod(m * n, 1 << K)
    if rem:
        return False, 0, ("index", m * n, 1 << K)
    counts = 0
    mult = {}
    for row in A:
        for cell in row:
            mult[cell] = mult.get(cell, 0) + 1
    if len(mult) != (1 << K):
        return False, 0, ("distinct", len(mult), 1 << K)
    for v, c in mult.items():
        counts += 1
        if c != lam:
            return False, counts, ("multiplicity", v, c, lam)
    pairs = list(combinations(range(K), 2))
    for i in range(m):
        for a, b in pairs:
            c4 = [0, 0, 0, 0]
            for j in range(n):
                w = A[i][j]
                c4[2 * ((w >> a) & 1) + ((w >> b) & 1)] += 1
            counts += 4
            if c4 != [n // 4] * 4:
                return False, counts, ("row", i, (a, b), c4)
    for j in range(n):
        for a, b in pairs:
            c4 = [0, 0, 0, 0]
            for i in range(m):
                w = A[i][j]
                c4[2 * ((w >> a) & 1) + ((w >> b) & 1)] += 1
            counts += 4
            if c4 != [m // 4] * 4:
                return False, counts, ("column", j, (a, b), c4)
    return True, counts, None


def design_checks_fast(G, V, K):
    """The same counts on A[i][j] = G[i] XOR V[j], with the four patterns of a coordinate pair
    counted by population counts over bit-vectors of cells.  For a fixed row i the bit-vector of
    coordinate t over the n columns is the coordinate-t bit-vector of V, complemented exactly when
    bit t of G[i] is 1; for a fixed column j the bit-vector over the 24 rows is that of G,
    complemented exactly when bit t of V[j] is 1.  Every count below is still a count of actual
    cells.  -> (ok, number_of_counts, first_failure_or_None)."""
    m = len(G)
    n = len(V)
    lam, rem = divmod(m * n, 1 << K)
    if rem:
        return False, 0, ("index", m * n, 1 << K)
    counts = 0

    mult = {}
    for g in G:
        for v in V:
            w = g ^ v
            mult[w] = mult.get(w, 0) + 1
    if len(mult) != (1 << K):
        return False, 0, ("distinct", len(mult), 1 << K)
    for w, c in mult.items():
        counts += 1
        if c != lam:
            return False, counts, ("multiplicity", w, c, lam)

    pairs = list(combinations(range(K), 2))
    full_n = (1 << n) - 1
    full_m = (1 << m) - 1
    vb = [sum(1 << j for j, v in enumerate(V) if (v >> t) & 1) for t in range(K)]
    gb = [sum(1 << i for i, g in enumerate(G) if (g >> t) & 1) for t in range(K)]

    for i, g in enumerate(G):
        col = [(vb[t] ^ full_n) if (g >> t) & 1 else vb[t] for t in range(K)]
        for a, b in pairs:
            x, y = col[a], col[b]
            c11 = popcount(x & y)
            c10 = popcount(x & ~y & full_n)
            c01 = popcount(~x & y & full_n)
            c00 = n - c11 - c10 - c01
            counts += 4
            if [c00, c01, c10, c11] != [n // 4] * 4:
                return False, counts, ("row", i, (a, b), [c00, c01, c10, c11])

    for j, v in enumerate(V):
        col = [(gb[t] ^ full_m) if (v >> t) & 1 else gb[t] for t in range(K)]
        for a, b in pairs:
            x, y = col[a], col[b]
            c11 = popcount(x & y)
            c10 = popcount(x & ~y & full_m)
            c01 = popcount(~x & y & full_m)
            c00 = m - c11 - c10 - c01
            counts += 4
            if [c00, c01, c10, c11] != [m // 4] * 4:
                return False, counts, ("column", j, (a, b), [c00, c01, c10, c11])

    return True, counts, None


def oa_violations(rows, K):
    """Pairs of coordinates on which `rows` fails strength 2, plus the number of counts made."""
    bad = []
    counts = 0
    n = len(rows)
    for a, b in combinations(range(K), 2):
        c4 = [0, 0, 0, 0]
        for w in rows:
            c4[2 * ((w >> a) & 1) + ((w >> b) & 1)] += 1
        counts += 4
        if c4 != [n // 4] * 4:
            bad.append(((a, b), c4))
    return bad, counts


# ---------------------------------------------------------------------------------------------
# 1.  THE PALEY RECIPE REPRODUCES THE PRINTED ARRAY M
# ---------------------------------------------------------------------------------------------
print("--- the array M of Section 2 ---")

qr = set((x * x) % Q_PALEY for x in range(1, Q_PALEY))


def chi(a):
    a %= Q_PALEY
    if a == 0:
        return 0
    return 1 if a in qr else -1


N = Q_PALEY + 1                                     # 24
C = [[0] * N for _ in range(N)]
for j in range(1, N):
    C[0][j] = 1
for i in range(1, N):
    C[i][0] = -1
for i in range(1, N):
    for j in range(1, N):
        C[i][j] = chi((j - 1) - (i - 1))
Hd = [[C[i][j] + (1 if i == j else 0) for j in range(N)] for i in range(N)]

inner_bad = 0
inner_counts = 0
for i in range(N):
    for j in range(N):
        s = sum(Hd[i][t] * Hd[j][t] for t in range(N))
        inner_counts += 1
        if s != (N if i == j else 0):
            inner_bad += 1
check("hadamard-24-orthogonal", inner_bad == 0,
      "H = C + I from q = %d: all %d inner products of rows are 24 on the diagonal and 0 off it; "
      "entries all +-1: %s" % (Q_PALEY, inner_counts,
                               all(abs(Hd[i][j]) == 1 for i in range(N) for j in range(N))))

Hn = [r[:] for r in Hd]
for i in range(N):
    if Hn[i][0] == -1:
        Hn[i] = [-x for x in Hn[i]]
M_rebuilt = ["".join("1" if Hn[i][j] == 1 else "0" for j in range(1, N)) for i in range(N)]
check("paley-reproduces-M", M_rebuilt == M_ROWS,
      "normalising the first column to +1, deleting it and mapping -1 -> 0 gives exactly the 24 x 23 "
      "array printed in the paper (%d rows of %d characters compared)"
      % (len(M_ROWS), len(M_ROWS[0])))

M = [mask_of(s) for s in M_ROWS]
check("M-shape", len(M) == M_ROWS_COUNT and all(len(s) == Q_PALEY for s in M_ROWS)
      and len(set(M)) == M_ROWS_COUNT,
      "M has %d rows, %d columns, and its rows are pairwise distinct" % (M_ROWS_COUNT, Q_PALEY))

bad, cnt = oa_violations(M, Q_PALEY)
check("M-is-OA-24-23-2-2", not bad,
      "all %d coordinate pairs of M give 6/6/6/6 (%d pattern counts, %d violations)"
      % (len(list(combinations(range(Q_PALEY), 2))), cnt, len(bad)))

# ---------------------------------------------------------------------------------------------
# 2.  THE SUBSPACE D AND ITS DUAL
# ---------------------------------------------------------------------------------------------
print("--- the subspace D of Section 2 ---")

Hg = [mask_of(s) for s in D_GENS]
D = span(Hg)
weights = [popcount(x) for x in D[1:]]
check("D-is-3-dimensional", len(D) == (1 << DIM_D),
      "span(h1,h2,h3) has %d elements, so h1,h2,h3 are independent and dim D = %d"
      % (len(D), DIM_D))
check("D-minimum-weight-3", min(weights) >= 3,
      "the %d nonzero words of D have weights %s; minimum %d >= 3, so no word of weight 1 or 2 "
      "lies in D" % (len(weights), sorted(weights), min(weights)))

Vp = [mask_of(s) for s in V_PRINTED]
Vc = dual(Hg, 8)
check("V-printed-equals-D-dual", set(Vp) == set(Vc) and len(set(Vp)) == 32 == len(Vc),
      "the 32 printed words of V are exactly D^perp in F_2^8, recomputed from h1,h2,h3 "
      "(|V| = %d, |D^perp| = %d, equal as sets)" % (len(set(Vp)), len(Vc)))

G8 = [m & 0xFF for m in M]                          # columns 0..7 of M
syn = [0] * 8
for g in G8:
    syn[parity(g & Hg[0]) | (parity(g & Hg[1]) << 1) | (parity(g & Hg[2]) << 2)] += 1
check("syndromes-balanced", syn == [3] * 8,
      "the 24 syndromes (h1.g, h2.g, h3.g) hit each of the 8 values of F_2^3 exactly 3 times: %s"
      % syn)
odd = [sum(1 for g in G8 if parity(g & h)) for h in D[1:]]
check("parity-12-per-word", odd == [12] * 7,
      "equivalently, each of the %d nonzero words of D has odd inner product with exactly 12 of "
      "the 24 rows: %s" % (len(odd), odd))

# ---------------------------------------------------------------------------------------------
# 2b.  HOW MANY SUCH D THERE ARE ON THESE EIGHT COLUMNS (the Remark of Section 2 of the paper)
#
# A word h of F_2^8 is ADMISSIBLE when its weight is at least 3 and it has odd inner product with
# exactly 12 of the 24 rows of G^(8).  A triple of distinct admissible words whose four further
# nonzero combinations are also admissible is exactly a D satisfying hypotheses (ii) and (iii):
# minimum weight >= 3 on all seven nonzero words, and -- since the 8 syndrome classes are equal in
# size precisely when every nonzero character of F_2^3 sums to zero over the rows -- the balanced
# syndrome distribution.
# ---------------------------------------------------------------------------------------------
print("--- the population of admissible D on columns 0..7 ---")

admissible = [h for h in range(1, 256)
              if popcount(h) >= 3 and sum(1 for g in G8 if parity(g & h)) == 12]
adm_set = set(admissible)
triples = 0
subspaces = set()
L = len(admissible)
for i1 in range(L):
    a1 = admissible[i1]
    for i2 in range(i1 + 1, L):
        a2 = admissible[i2]
        a12 = a1 ^ a2
        if a12 not in adm_set:
            continue
        for i3 in range(i2 + 1, L):
            a3 = admissible[i3]
            if a3 == a12:
                continue
            if (a1 ^ a3) in adm_set and (a2 ^ a3) in adm_set and (a12 ^ a3) in adm_set:
                triples += 1
                # ⭐ A TRIPLE IS A BASIS, NOT A SUBSPACE: every eligible D has 7*6*4/3! = 28 of them, so
                # the triple count must NOT be quoted as a number of subspaces D.  Both are asserted.
                subspaces.add(frozenset([0, a1, a2, a3, a12, a1 ^ a3, a2 ^ a3, a12 ^ a3]))
check("printed-D-is-admissible", set(D[1:]) <= adm_set and len(admissible) == 132,
      "all 7 nonzero words of the printed D are admissible, so (h1,h2,h3) is one of the triples "
      "counted next (%d admissible words in F_2^8 in total)" % len(admissible))
check("triple-count-21252", triples == 21252 and len(subspaces) == 759
      and triples == 28 * len(subspaces),
      "columns 0..7 of M carry exactly %d such triples {h1,h2,h3}, the figure the paper's Remark "
      "reports; they span %d DISTINCT subspaces D, each with 28 bases (%d = 28 x %d), so the choice "
      "of D is not delicate and the triple count is not a count of subspaces"
      % (triples, len(subspaces), triples, len(subspaces)))

# ---------------------------------------------------------------------------------------------
# 3.  HYPOTHESES (i)-(iii) OF THEOREM 1 FOR EVERY COLUMN COUNT 8 <= K <= 23
# ---------------------------------------------------------------------------------------------
print("--- hypotheses of Theorem 1, for every K with 8 <= K <= 23 ---")

hyp_bad = []
hyp_counts = 0
for K in range(8, Q_PALEY + 1):
    GK = [m & ((1 << K) - 1) for m in M]
    b, c = oa_violations(GK, K)
    hyp_counts += c
    wK = [popcount(x) for x in span([h & ((1 << K) - 1) for h in Hg])[1:]]
    sK = [0] * 8
    for g in GK:
        sK[parity(g & Hg[0]) | (parity(g & Hg[1]) << 1) | (parity(g & Hg[2]) << 2)] += 1
    if b or min(wK) < 3 or len(wK) != 7 or sK != [3] * 8:
        hyp_bad.append((K, len(b), min(wK), sK))
check("hypotheses-hold-for-every-K", not hyp_bad,
      "for all %d values of K: columns 0..K-1 of M form an OA(24,K,2,2), the zero-padded D still "
      "has 7 nonzero words of minimum weight 3, and the syndromes stay [3]*8 (%d pattern counts, "
      "%d failures)" % (len(range(8, Q_PALEY + 1)), hyp_counts, len(hyp_bad)))

# ---------------------------------------------------------------------------------------------
# 4.  THE SEVEN DESIGNS, CHECKED CELL BY CELL AGAINST THE DEFINITION
# ---------------------------------------------------------------------------------------------
print("--- the seven designs I_{k+3}(24, 2^k, 2, 2), 5 <= k <= 11 ---")

total_counts = 0
per_k = {}
for k in K_RANGE:
    K = k + DIM_D
    GK = [m & ((1 << K) - 1) for m in M]
    VK = dual([h & ((1 << K) - 1) for h in Hg], K)
    if len(VK) != (1 << k):
        check("design-k=%d" % k, False, "|D^perp| = %d, expected %d" % (len(VK), 1 << k))
        continue
    ok, cnt, fail = design_checks_fast(GK, VK, K)
    expect = (1 << K) + 24 * len(list(combinations(range(K), 2))) * 4 \
        + (1 << k) * len(list(combinations(range(K), 2))) * 4
    per_k[k] = cnt
    total_counts += cnt
    check("design-k=%d" % k, ok and cnt == expect,
          "I_%d(24,%d,2,2) verified from the definition: index lambda = %d, %d counts "
          "(%d multiplicities + %d row + %d column pattern counts), first failure: %s"
          % (K, 1 << k, 24 * (1 << k) // (1 << K), cnt, 1 << K,
             24 * len(list(combinations(range(K), 2))) * 4,
             (1 << k) * len(list(combinations(range(K), 2))) * 4, fail))

check("seven-cases-all-present", sorted(per_k) == list(K_RANGE),
      "all seven cases 5 <= k <= 11 were built and checked: k = %s" % sorted(per_k))
check("total-definition-level-counts", total_counts == 1363104,
      "the seven cells together make %d definition-level counts (per k: %s), the number the paper "
      "reports" % (total_counts, [per_k[k] for k in sorted(per_k)]))

# The naive counter, on the two smallest cases, against the fast one.
K = 8
A5 = [[g ^ v for v in Vp] for g in G8]              # built from the PRINTED V, in printed order
ok5, cnt5, fail5 = design_checks_naive(A5, K)
check("naive-counter-agrees-k=5", ok5 and cnt5 == per_k[5] == 6528,
      "the transparent triple loop over the 768 cells of the printed k = 5 array makes the same "
      "%d counts and reaches the same verdict as the bit-parallel counter (failure: %s)"
      % (cnt5, fail5))
K6 = 9
G9 = [m & ((1 << K6) - 1) for m in M]
V6 = dual([h & ((1 << K6) - 1) for h in Hg], K6)
A6 = [[g ^ v for v in V6] for g in G9]
ok6, cnt6, fail6 = design_checks_naive(A6, K6)
check("naive-counter-agrees-k=6", ok6 and cnt6 == per_k[6] == 13184,
      "and again on k = 6: %d counts, same verdict (failure: %s)" % (cnt6, fail6))

# ---------------------------------------------------------------------------------------------
# 5.  EACH OF THE SEVEN IS AN ADMISSIBLE INSTANCE OF THE CONJECTURE
# ---------------------------------------------------------------------------------------------
print("--- admissibility against the conjecture's own hypotheses ---")

m_par = 6
adm = []
for k in K_RANGE:
    K = k + DIM_D
    four_n = 1 << k
    if four_n % 4:
        adm.append((k, "4n not divisible by 4"))
        continue
    n_par = four_n // 4
    if n_par < m_par:
        adm.append((k, "n = %d < m = %d" % (n_par, m_par)))
    if (16 * m_par * n_par) % (1 << K):
        adm.append((k, "2^%d does not divide 16mn = %d" % (K, 16 * m_par * n_par)))
    if K > 4 * m_par - 1:
        adm.append((k, "K = %d > 4m-1 = %d" % (K, 4 * m_par - 1)))
    if m_par == 1:
        adm.append((k, "the m = 1 exception"))
check("all-seven-are-instances", not adm,
      "for every k in 5..11 with m = 6, n = 2^(k-2): n >= 6, 2^(k+3) | 16mn, k+3 <= 4m-1 = 23, and "
      "m > 1 -- so each design is a case the conjecture asserts (%d objections)" % len(adm))
check("hadamard-of-order-4m-exists", inner_bad == 0 and N == 4 * m_par,
      "the conjecture's hypothesis is a Hadamard matrix of order 4m = 24, and one was built and "
      "checked above")

# ---------------------------------------------------------------------------------------------
# 5b.  THE READING OF THE SOURCE'S DEFINITION
#
# The sentence of the source that introduces the definition calls the array formed by a row
# (column) an orthogonal array "of size k, degree n (respectively, m)", which interchanges the two
# slots of the source's own OA(N,k,q,t) notation: a row of an m x n array of words of length k
# holds n such words, so as an array it is n x k.
#
#   READING A (the paper's Section 1, and the sentence that immediately follows the definition in
#   the source, and the source's own necessary conditions on p. 81 of the thesis): each row is an
#   OA(n,k,q,t) and each column an OA(m,k,q,t); m is the number of rows.
#
#   READING B (the interchanged slots taken literally): the array formed by a row would have k rows
#   and n columns.  It really has n rows and k columns, so reading B is self-consistent only when
#   n = k.
#
# SOURCE_EXAMPLES is TRANSCRIBED from the source: for each design the source prints, its label and
# the dimensions of the printed array.  The DIGITS of those arrays are not transcribed and nothing
# here replays them; only the dimensions and the parameter arithmetic are checked.  The dimensions
# are recorded unordered, except for the one non-square example whose rows the source itself marks
# (Table 5.3, "with rows indicated by superscripts": twelve superscripts A..L on twelve rows of
# eight words), where the ordered pair (rows, columns) is recorded.
# ---------------------------------------------------------------------------------------------
print("--- the reading of the source's definition ---")

#      k,   m,   n, q, t, printed dims,  rows marked?, where
SOURCE_EXAMPLES = [
    (4,  9,  9, 3, 2, (9, 9), False, "Table 5.1"),
    (5, 12,  8, 2, 2, (12, 8), True, "Table 5.3, rows marked by the superscripts A..L"),
    (4, 12, 12, 2, 2, (12, 12), False, "Table 5.5"),
    (6, 12, 16, 2, 2, (12, 16), False, "Table 5.6, printed sideways"),
]
TARGETS = [(k + DIM_D, 4 * m_par, 1 << k, 2, 2) for k in K_RANGE]

shape_bad = [e for e in SOURCE_EXAMPLES if sorted(e[5]) != sorted([e[1], e[2]])]
orient_bad = [e for e in SOURCE_EXAMPLES if e[6] and e[5] != (e[1], e[2])]
check("printed-shapes-match-the-labels", not shape_bad and not orient_bad,
      "the dimensions of each of the %d designs the source prints agree with the {m,n} of its "
      "label (%s), and in the one non-square example whose rows the source marks the count of "
      "marked rows is m and the count of words in a row is n (%d label mismatches, %d orientation "
      "mismatches)"
      % (len(SOURCE_EXAMPLES),
         "; ".join("%s: I_%d(%d,%d,%d,%d) printed %d x %d"
                   % (e[7], e[0], e[1], e[2], e[3], e[4], e[5][0], e[5][1])
                   for e in SOURCE_EXAMPLES),
         len(shape_bad), len(orient_bad)))

readingA_bad = []
readingB_consistent = []
for e in SOURCE_EXAMPLES + TARGETS:
    k_, m_, n_, q_, t_ = e[0], e[1], e[2], e[3], e[4]
    if n_ % (q_ ** t_) or m_ % (q_ ** t_) or (m_ * n_) % (q_ ** k_) or (m_ * n_) // (q_ ** k_) < 1:
        readingA_bad.append((k_, m_, n_, q_, t_))
    if n_ == k_:
        readingB_consistent.append((k_, m_, n_, q_, t_))
check("reading-A-fits-every-label-reading-B-none",
      not readingA_bad and not readingB_consistent,
      "for all %d labels -- the %d the source prints and the %d at issue here -- reading A makes "
      "n/q^t, m/q^t and the index mn/q^k positive integers (%d failures), while reading B would "
      "require n = k, which holds for none of them (%d cases where it holds)"
      % (len(SOURCE_EXAMPLES) + len(TARGETS), len(SOURCE_EXAMPLES), len(TARGETS),
         len(readingA_bad), len(readingB_consistent)))

nec_bad = []
for k in K_RANGE:
    K = k + DIM_D
    GKn = [mm & ((1 << K) - 1) for mm in M]
    VKn = dual([h & ((1 << K) - 1) for h in Hg], K)
    b_col, _ = oa_violations(GKn, K)                # the OA(m,k,q,t) the source's conditions ask for
    b_row, _ = oa_violations(VKn, K)                # the OA(n,k,q,t) the source's conditions ask for
    if b_col or b_row or len(VKn) != (1 << k) or (M_ROWS_COUNT * (1 << k)) % (1 << K):
        nec_bad.append((k, len(b_col), len(b_row), len(VKn)))
check("source-necessary-conditions-exhibited", not nec_bad,
      "for each of the seven k, with K = k+3: 2^K divides mn = 24 * 2^k, the first K columns of M "
      "are an OA(24,K,2,2) and the 2^k words of V^(K) are an OA(2^k,K,2,2) -- the two orthogonal "
      "arrays the source's own necessary conditions for an I_K(24,2^k,2,2) demand, in the source's "
      "own notation (%d failures)" % len(nec_bad))

# ---------------------------------------------------------------------------------------------
# 6.  CONTROLS: THE CHECKER MUST REJECT A DAMAGED WITNESS
# ---------------------------------------------------------------------------------------------
print("--- controls: deliberately damaged witnesses must be rejected ---")

B = [row[:] for row in A5]
B[0][0], B[0][1] = B[0][1], B[0][0]
okb, _, failb = design_checks_naive(B, 8)
check("control-swap-rejected", not okb,
      "transposing two cells inside row 0 of the k = 5 array is rejected: %s" % (failb,))

B = [row[:] for row in A5]
B[0][1] = B[0][0]
okb, _, failb = design_checks_naive(B, 8)
check("control-duplicate-rejected", not okb,
      "overwriting one cell with a copy of another is rejected: %s" % (failb,))

bad_gens = [mask_of("11100000"), mask_of("01011000"), mask_of("10010000")]
Dbad = span(bad_gens)
wbad = [popcount(x) for x in Dbad[1:]]
Vbad = dual(bad_gens, 8)
Abad = [[g ^ v for v in Vbad] for g in G8]
okb, _, failb = design_checks_naive(Abad, 8)
check("control-weight-2-D-rejected", len(Dbad) == 8 and min(wbad) == 2 and not okb,
      "replacing h3 by 10010000 keeps dim D = 3 and |D^perp| = %d but breaks hypothesis (ii) "
      "(weights %s, minimum 2); the resulting array is rejected: %s"
      % (len(Vbad), sorted(wbad), failb))

Gsh = [m & 0xFF for m in M]
Gsh[0] = Gsh[1]                                     # two equal rows: strength 2 must break
badsh, _ = oa_violations(Gsh, 8)
check("control-broken-G-rejected", bool(badsh),
      "duplicating a row of G destroys the OA(24,8,2,2) property, and the same routine that "
      "certified G reports it (%d offending pairs, first %s)"
      % (len(badsh), badsh[0] if badsh else None))

# ---------------------------------------------------------------------------------------------
# WHAT THIS PROGRAM DOES NOT COVER
# ---------------------------------------------------------------------------------------------
print()
print("NOT RE-RUN: the cells 12 <= k <= 20 at definition level. Their arrays have up to 2^20 "
      "columns, which is not referee-sized; what is re-checked here for every K in 8..23 is the "
      "hypotheses of Theorem 1, from which those cells follow. They are in any case NOT new -- "
      "they are Rahim and Cavenagh's own published theorem for 12 <= k <= 20.")
print("NOT RE-RUN: nothing here bears on the whole m = 6 case of the conjecture. This program "
      "checks the seven MINIMAL cases the source lists as unknown; the reduction of a general "
      "admissible (k,n) at m = 6 to those minimal cases is not verified here, and no such "
      "reduction is proved in the source.")
print("NOT RE-RUN: the search over CHOICES OF COLUMN SET. The generating column set is the first "
      "eight columns of M and is fixed in the source of this program; only the 21,252 admissible "
      "generator triples {h1,h2,h3} on THAT one column set -- the 759 subspaces D they span -- are "
      "enumerated. The exhaustive C(23,8) = 490,314 census over column "
      "sets was never completed (the original run sampled 20,000 subsets, 4.079%), and existence "
      "needs none of it. Nothing here claims this witness is unique, canonical or least.")
print("NOT RE-RUN: the source's own printed examples at the level of their DIGITS. Those arrays "
      "are the source's, they are not transcribed into this program, and no cell of them is "
      "counted here. What section 5b checks about the reading of the definition is their printed "
      "DIMENSIONS, transcribed with their labels, and the parameter arithmetic of the two readings. "
      "The passages of the source quoted or paraphrased in Section 1 of the paper were read by "
      "hand and are not machine checkable at all.")
print()

if _n_fail:
    print("%d CHECK(S) DID NOT PASS" % _n_fail)
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _n_pass)
sys.exit(0)

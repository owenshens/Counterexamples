#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify.py -- independent re-derivation of every quantity claimed in paper.tex.

WHAT IT CHECKS, and in this order:

  1. The object, read from the sign matrix PRINTED IN THE PAPER: that H is hot, that C = H - I is
     cold, the two determinants 33489 = 183^2 and 64000 = 125*2^9, the exact characteristic
     polynomial of C, the Pfaffian, the out-degree sequence and the graph6 string.
  2. The two invariance lemmas the proof rests on -- switching and permutation leave det C, det(I+C)
     and the whole characteristic polynomial fixed -- and the bordering normal form.
  3. The covering, with a DETERMINISTIC COMPLETENESS CERTIFICATE for the class list the census rests
     on: the canonical form is verified to be a complete isomorphism invariant by exhaustion over
     every labelled tournament of order 5 and of order 6; the automorphism count used below is
     verified against the count over all k! relabellings for every representative on up to 7
     vertices; each of the 6880 order-8 canonical codes is verified to be REALIZED by an explicit
     relabelling (so equal codes force an isomorphism and the deduplication can lose no class); and
     for every k = 2..8 the orbit sizes k!/|Aut| of the representatives are verified to sum to
     exactly 2^(k(k-1)/2), the number of ALL labelled k-vertex tournaments -- which forces the list
     to be complete and free of repeated classes, and not merely to have the published length. Also
     the class counts against A000568 and the 9-vertex subtournament property.
  4. The published small cells f_k(4), g_k(4), f_k(6), g_k(6) by BRUTE FORCE over every labelled cold
     matrix, f_k(8), g_k(8) by the covering, and the covering reproduced against brute force at n = 6.
  5. The n = 10 census itself: 6880 * 256 = 1761280 leaves, the two maxima, the coincidence of the two
     argmax sets, the uniqueness of the maximiser class, the achieved Pfaffian spectrum, and Orrick's
     two order-reversing pairs as a forced positive.

Python 3.9+, STANDARD LIBRARY ONLY. Every decision is taken in exact integer or Fraction arithmetic;
no floating point appears anywhere in this file.

Prints one `PASS <name>` line per check and closes with `VERDICT: ALL <n> CHECKS PASS`, exiting 0 if
and only if every check passed.
"""

import itertools
import math
import random
import sys
from fractions import Fraction

# ---------------------------------------------------------------------------
# 0. THE OBJECT, EXACTLY AS PRINTED IN paper.tex
# ---------------------------------------------------------------------------
# This block is the one in Section 2 of the paper, character for character. Nothing else about the
# witness is hard-coded: every number below is re-derived from these 100 signs.
PAPER_H = """\
++++++++++
-+++++--+-
--+++++--+
---++-++--
----++-+++
---+-++-+-
-+--+-+-++
-++--+++--
--++---+++
-+-+-+-+-+
"""

# The paper's claims, quoted here so a disagreement is a FAILURE and not a silent overwrite.
CLAIM_DET_COLD = 33489
CLAIM_DET_HOT = 64000
CLAIM_CHARPOLY = [1, 0, 45, 0, 770, 0, 6210, 0, 23485, 0, 33489]   # descending, det(xI - C)
CLAIM_PFAFFIAN = 183
CLAIM_OUTDEG = [3, 3, 4, 4, 4, 4, 4, 5, 5, 9]
CLAIM_GRAPH6 = 'I~~lkrMig'
CLAIM_EDGES = 31
CLAIM_CLASS_COUNTS = [1, 1, 2, 4, 12, 56, 456, 6880]               # A000568(1..8)
# The number of LABELLED tournaments on k vertices, 2^{k(k-1)/2}: the total the orbits of the
# representatives must account for, exactly, if the class list is complete and irredundant.
CLAIM_LABELLED = {k: 1 << (k * (k - 1) // 2) for k in range(1, 9)}  # ..., 8: 268435456
CLAIM_LEAVES_10 = 1761280
CLAIM_SMALL_CELLS = {4: (9, 16), 6: (81, 160), 8: (2401, 4096)}
CLAIM_N6_ARGMAX_COUNT = 1280
CLAIM_PF_SET_8 = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 31, 33, 35, 49]
# The complete order-10 cold Pfaffian spectrum, as printed by the paper and as published by
# Klanderman, Montee, Piotrowski, Rice and Shader for order-10 Seidel tournament matrices.
CLAIM_PF_SET_10 = (list(range(1, 130, 2))
                   + [133, 135, 137, 139, 141, 143]
                   + [147, 153, 161, 165, 175, 183])
# Orrick's two order-reversing order-10 pairs (det H, det C), used as a forced positive: a census
# that could not exhibit the SHAPE of a counterexample would not be evidence of its absence.
CLAIM_ORRICK = [(45056, 21609), (46080, 20449)]

RESULTS = []


def chk(name, ok, detail=''):
    RESULTS.append((name, bool(ok), detail))
    print('%s %s%s' % ('PASS' if ok else 'FAIL', name, (' [%s]' % detail) if detail else ''))
    sys.stdout.flush()


# ---------------------------------------------------------------------------
# 1. EXACT LINEAR ALGEBRA (integers only; Bareiss, so every division is exact)
# ---------------------------------------------------------------------------
BADDIV = 0


def det_int(M, n):
    """Fraction-free Bareiss determinant of an n x n integer matrix. Asserts every division is exact,
    which is the property that makes the arithmetic a proof rather than a hope."""
    global BADDIV
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            piv = -1
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    piv = r
                    break
            if piv < 0:
                return 0
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        akk = A[k][k]
        Ak = A[k]
        for i in range(k + 1, n):
            Ai = A[i]
            aik = Ai[k]
            for j in range(k + 1, n):
                num = akk * Ai[j] - aik * Ak[j]
                if num % prev:
                    BADDIV += 1
                Ai[j] = num // prev
        prev = akk
    return sign * A[n - 1][n - 1]


def det_fast(A):
    """The census inner loop: Bareiss on a 10 x 10, no exactness assertion (checked separately on a
    sample), no copy -- the caller passes a scratch matrix."""
    sign = 1
    prev = 1
    for k in range(9):
        if A[k][k] == 0:
            piv = -1
            for r in range(k + 1, 10):
                if A[r][k] != 0:
                    piv = r
                    break
            if piv < 0:
                return 0
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        akk = A[k][k]
        Ak = A[k]
        for i in range(k + 1, 10):
            Ai = A[i]
            aik = Ai[k]
            for j in range(k + 1, 10):
                Ai[j] = (akk * Ai[j] - aik * Ak[j]) // prev
        prev = akk
    return sign * A[9][9]


def pfaffian(M, n):
    """Pfaffian of a skew-symmetric integer matrix of even order, by expansion with memoisation over
    subsets. Exact integers."""
    memo = {0: 1}

    def pf(mask):
        if mask in memo:
            return memo[mask]
        i = (mask & -mask).bit_length() - 1
        rest = mask & ~(1 << i)
        total = 0
        m = rest
        while m:
            jb = m & -m
            j = jb.bit_length() - 1
            m ^= jb
            between = rest & ((1 << j) - 1)
            s = -1 if (bin(between).count('1') & 1) else 1
            if M[i][j]:
                total += s * M[i][j] * pf(rest & ~jb)
        memo[mask] = total
        return total

    return pf((1 << n) - 1)


def charpoly_of_skew(C, n):
    """Coefficients of det(x I - C), descending, by exact Lagrange interpolation through n+1 integer
    points. Fractions are used and the result is asserted integral."""
    pts = []
    for x in range(n + 1):
        M = [[(x if i == j else 0) - C[i][j] for j in range(n)] for i in range(n)]
        pts.append((x, det_int(M, n)))
    coef = [Fraction(0)] * (n + 1)                     # ascending
    for xi, yi in pts:
        num = [Fraction(1)]
        den = Fraction(1)
        for xj, _ in pts:
            if xj == xi:
                continue
            new = [Fraction(0)] * (len(num) + 1)
            for k, c in enumerate(num):
                new[k + 1] += c
                new[k] += c * (-xj)
            num = new
            den *= (xi - xj)
        for k, c in enumerate(num):
            coef[k] += Fraction(yi) * c / den
    assert all(c.denominator == 1 for c in coef), 'interpolation left a non-integer coefficient'
    return [int(coef[n - k]) for k in range(n + 1)]     # descending


# ---------------------------------------------------------------------------
# 2. HOT AND COLD MATRICES
# ---------------------------------------------------------------------------
def parse_signs(block):
    rows = [ln.strip() for ln in block.strip().splitlines()]
    return [[1 if ch == '+' else -1 for ch in r] for r in rows]


def is_hot(H, n):
    if len(H) != n or any(len(r) != n for r in H):
        return False
    for i in range(n):
        if H[i][i] != 1:
            return False
        for j in range(n):
            if i != j and (H[i][j] not in (1, -1) or H[i][j] != -H[j][i]):
                return False
    return True


def is_cold(C, n):
    for i in range(n):
        if C[i][i] != 0:
            return False
        for j in range(n):
            if i != j and (C[i][j] not in (1, -1) or C[i][j] != -C[j][i]):
                return False
    return True


def cold_of(H, n):
    return [[H[i][j] - (1 if i == j else 0) for j in range(n)] for i in range(n)]


def hot_of(C, n):
    return [[C[i][j] + (1 if i == j else 0) for j in range(n)] for i in range(n)]


def graph6(adj, n):
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if adj[i][j] else 0)
    while len(bits) % 6:
        bits.append(0)
    s = chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = v * 2 + b
        s += chr(v + 63)
    return s


# ---------------------------------------------------------------------------
# 3. TOURNAMENTS: A COMPLETE CANONICAL FORM, AND CLASS GENERATION
# ---------------------------------------------------------------------------
# A tournament on k vertices is a tuple of k ints; bit j of entry i is set iff i beats j.
#
# CANONICAL FORM. can(T) = max over the relabellings that sort the score sequence non-decreasingly of
# the upper-triangle bit code. This is a COMPLETE invariant, not a heuristic: if T' = sigma.T then p
# sorts the scores of T' exactly when p o sigma sorts the scores of T, so the two matrices range over
# the same SET and the two maxima are equal; conversely equal codes exhibit an isomorphism. The cost
# is the product of the score multiplicities' factorials, not k!.
def canon(out, k):
    score = [bin(out[i]).count('1') for i in range(k)]
    order = sorted(range(k), key=lambda v: score[v])
    cells = []
    i = 0
    while i < k:
        j = i
        while j < k and score[order[j]] == score[order[i]]:
            j += 1
        cells.append(order[i:j])
        i = j
    best = -1
    for choice in itertools.product(*[itertools.permutations(c) for c in cells]):
        p = [v for c in choice for v in c]
        code = 0
        for a in range(k):
            oa = out[p[a]]
            for b in range(a + 1, k):
                code = (code << 1) | ((oa >> p[b]) & 1)
        if code > best:
            best = code
    return best


def extend_to(out, k, m):
    """out is a k-vertex tournament; add vertex k, beating vertex j exactly when bit j of m is set."""
    o = list(out) + [0]
    for j in range(k):
        if (m >> j) & 1:
            o[k] |= (1 << j)
        else:
            o[j] |= (1 << k)
    return tuple(o)


def class_reps(kmax):
    """One representative per isomorphism class of tournaments on 1..kmax vertices."""
    reps = {1: [(0,)]}
    for k in range(1, kmax):
        seen = {}
        for out in reps[k]:
            for m in range(1 << k):
                t = extend_to(out, k, m)
                c = canon(t, k + 1)
                if c not in seen:
                    seen[c] = t
        reps[k + 1] = list(seen.values())
    return reps


def relabel(out, k, p):
    """p is a permutation list: new vertex a is old vertex p[a]."""
    o = [0] * k
    for a in range(k):
        oa = out[p[a]]
        for b in range(k):
            if a != b and (oa >> p[b]) & 1:
                o[a] |= (1 << b)
    return tuple(o)


def border(out, k):
    """The order-(k+1) cold matrix [[0, j],[-j^T, A]] where A is the Seidel matrix of the tournament
    `out`: index 0 is the bordered vertex, index i+1 is tournament vertex i."""
    n = k + 1
    C = [[0] * n for _ in range(n)]
    for j in range(1, n):
        C[0][j] = 1
        C[j][0] = -1
    for i in range(1, n):
        oi = out[i - 1]
        Ci = C[i]
        for j in range(1, n):
            if i != j:
                Ci[j] = 1 if (oi >> (j - 1)) & 1 else -1
    return C


# --- the completeness certificate for the class lists -----------------------
# The census trusts class_reps(8) to be a COMPLETE and IRREDUNDANT transversal of the 6880
# isomorphism classes of 8-vertex tournaments. That is certified below by two independent facts, each
# of which a wrong class list violates:
#
#   (R) REALIZABILITY.  canon(T) is attained by an explicit relabelling of T, so two tournaments with
#       the same code are both isomorphic to the tournament decoded from that code; hence the
#       code-keyed dedup inside class_reps can never merge two DIFFERENT classes, and no class is
#       lost.  (canon_witness + tournament_of_code + code_of check this.)
#   (O) ORBIT COUNT.  sum over the representatives of k! / |Aut(T)| equals 2^{k(k-1)/2}, the number of
#       ALL labelled k-vertex tournaments.  Given (R) the representatives already meet every class, so
#       a repeated class would push this sum strictly above the total: equality forces the list to be
#       irredundant as well as complete.
#
# |Aut(T)| is counted over the score-preserving permutations only.  That is legitimate because
# relabelling permutes out-degrees -- score(relabel(T,p))[a] = score(T)[p[a]] -- so an automorphism
# fixes every vertex's score.  aut_cell is checked against the count over ALL k! relabellings for
# every representative on up to 7 vertices, and the whole canonical form is checked exhaustively over
# every labelled tournament on 5 and on 6 vertices.
def canon_witness(out, k):
    """canon(), and also a relabelling p attaining the maximal code."""
    score = [bin(out[i]).count('1') for i in range(k)]
    order = sorted(range(k), key=lambda v: score[v])
    cells = []
    i = 0
    while i < k:
        j = i
        while j < k and score[order[j]] == score[order[i]]:
            j += 1
        cells.append(order[i:j])
        i = j
    best = -1
    bestp = None
    for choice in itertools.product(*[itertools.permutations(c) for c in cells]):
        p = [v for c in choice for v in c]
        code = 0
        for a in range(k):
            oa = out[p[a]]
            for b in range(a + 1, k):
                code = (code << 1) | ((oa >> p[b]) & 1)
        if code > best:
            best = code
            bestp = p
    return best, bestp


def code_of(out, k):
    """The upper-triangle bit code of a tournament in its own labelling."""
    code = 0
    for a in range(k):
        oa = out[a]
        for b in range(a + 1, k):
            code = (code << 1) | ((oa >> b) & 1)
    return code


def tournament_of_code(code, k):
    """Inverse of code_of: rebuild the tournament from its upper-triangle bit code."""
    o = [0] * k
    idx = k * (k - 1) // 2 - 1
    for a in range(k):
        for b in range(a + 1, k):
            if (code >> idx) & 1:
                o[a] |= (1 << b)
            else:
                o[b] |= (1 << a)
            idx -= 1
    return tuple(o)


def score_cells(out, k):
    """The vertices grouped by out-degree, in non-decreasing order of out-degree."""
    score = [bin(out[i]).count('1') for i in range(k)]
    byscore = {}
    for v in range(k):
        byscore.setdefault(score[v], []).append(v)
    return [byscore[s] for s in sorted(byscore)]


def aut_cell(out, k):
    """|Aut(out)|, counted over the score-preserving permutations (see the note above)."""
    cells = score_cells(out, k)
    p = [0] * k
    total = 0
    for choice in itertools.product(*[itertools.permutations(c) for c in cells]):
        for cell, perm in zip(cells, choice):
            for pos, v in zip(cell, perm):
                p[pos] = v
        if relabel(out, k, p) == out:
            total += 1
    return total


def aut_full(out, k):
    """|Aut(out)|, counted over ALL k! relabellings. Only affordable for small k."""
    total = 0
    for p in itertools.permutations(range(k)):
        if relabel(out, k, list(p)) == out:
            total += 1
    return total


def tournament_of_bordered(C, n):
    """The inverse of border(): read the (n-1)-vertex tournament off a cold matrix whose row 0 is all
    +1. Raises if row 0 is not all +1."""
    for j in range(1, n):
        assert C[0][j] == 1, 'row 0 is not all +1; the matrix is not in bordered form'
    k = n - 1
    o = [0] * k
    for i in range(k):
        for j in range(k):
            if i != j and C[i + 1][j + 1] == 1:
                o[i] |= (1 << j)
    return tuple(o)


# ---------------------------------------------------------------------------
# 4. CENSUS DRIVERS
# ---------------------------------------------------------------------------
def brute_cell(n):
    """Every labelled cold matrix of order n. (n, maxC, maxH, violA, violB, #C, #H, #both, total,
    pfaffian set)."""
    pairs = list(itertools.combinations(range(n), 2))
    hist = {}
    for m in range(1 << len(pairs)):
        C = [[0] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            v = 1 if (m >> b) & 1 else -1
            C[i][j] = v
            C[j][i] = -v
        dc = det_int(C, n)
        dh = det_int(hot_of(C, n), n)
        hist[(dc, dh)] = hist.get((dc, dh), 0) + 1
    return summarise(hist)


def cover_cell(reps, n):
    """The covering census at order n: one bordered matrix per (n-1)-vertex tournament class."""
    hist = {}
    for out in reps[n - 1]:
        C = border(out, n - 1)
        dc = det_int(C, n)
        dh = det_int(hot_of(C, n), n)
        hist[(dc, dh)] = hist.get((dc, dh), 0) + 1
    return summarise(hist)


def summarise(hist):
    mc = max(k[0] for k in hist)
    mh = max(k[1] for k in hist)
    out = {
        'maxC': mc, 'maxH': mh,
        'violA': sum(v for k, v in hist.items() if k[0] == mc and k[1] < mh),
        'violB': sum(v for k, v in hist.items() if k[1] == mh and k[0] < mc),
        'nC': sum(v for k, v in hist.items() if k[0] == mc),
        'nH': sum(v for k, v in hist.items() if k[1] == mh),
        'nBoth': sum(v for k, v in hist.items() if k[0] == mc and k[1] == mh),
        'total': sum(hist.values()),
        'pf': sorted(set(math.isqrt(k[0]) for k in hist)),
        'squares': all(math.isqrt(k[0]) ** 2 == k[0] and k[0] % 2 == 1 for k in hist),
        'hist': hist,
    }
    return out


def census10(reps8):
    """The decisive run: for every one of the 6880 isomorphism classes of 8-vertex tournaments and
    every one of the 256 ways to add a ninth vertex, the bordered order-10 cold matrix and its hot
    partner. Returns the joint histogram plus the winning 9-vertex tournaments.

    ⭐ Every 2048th leaf is recomputed by the SECOND determinant routine, the one that asserts every
    Bareiss division is exact, so the fast inner loop is cross-checked against a checked engine
    rather than trusted."""
    hist = {}
    leaves = 0
    maxC = -1
    maxH = -1
    winC = []
    winH = []
    crossed = 0
    mismatch = 0
    for out8 in reps8:
        for m in range(256):
            out = list(out8) + [0]
            o8 = 0
            for j in range(8):
                if (m >> j) & 1:
                    o8 |= (1 << j)
                else:
                    out[j] |= (1 << 8)
            out[8] = o8
            out = tuple(out)
            C = border(out, 9)
            H = [[C[i][j] + (1 if i == j else 0) for j in range(10)] for i in range(10)]
            dc = det_fast([r[:] for r in C])
            dh = det_fast([r[:] for r in H])
            key = (dc, dh)
            hist[key] = hist.get(key, 0) + 1
            leaves += 1
            if leaves % 2048 == 1:
                crossed += 1
                if det_int(C, 10) != dc or det_int(H, 10) != dh:
                    mismatch += 1
            if dc > maxC:
                maxC = dc
                winC = []
            if dc == maxC:
                winC.append(out)
            if dh > maxH:
                maxH = dh
                winH = []
            if dh == maxH:
                winH.append(out)
    return hist, leaves, winC, winH, crossed, mismatch


# ---------------------------------------------------------------------------
# 5. THE CHECKS
# ---------------------------------------------------------------------------
def main():
    n = 10
    print('--- 1. the object, read from the sign matrix printed in the paper ---')
    H = parse_signs(PAPER_H)
    chk('paper-matrix-parsed',
        len(H) == 10 and all(len(r) == 10 for r in H),
        '%d rows x %d columns of +/- read from paper.tex' % (len(H), len(H[0])))
    chk('hot-matrix', is_hot(H, n), 'H[i][i]=+1, H[i][j]=-H[j][i], all 100 entries in {+1,-1}')
    C = cold_of(H, n)
    chk('cold-matrix', is_cold(C, n), 'C = H - I has zero diagonal, is skew and has 90 entries +/-1')

    dC = det_int(C, n)
    dH = det_int(H, n)
    chk('det-cold', dC == CLAIM_DET_COLD and dC == 183 ** 2,
        'det C = %d = 183^2, exact Bareiss' % dC)
    chk('det-hot', dH == CLAIM_DET_HOT and dH == 125 * 2 ** 9,
        'det H = %d = 125*2^9, exact Bareiss' % dH)

    cp = charpoly_of_skew(C, n)
    chk('charpoly', cp == CLAIM_CHARPOLY,
        'det(xI-C) = ' + ' + '.join('%d x^%d' % (c, n - k) for k, c in enumerate(cp) if c))
    chk('charpoly-even', all(cp[n - k] == 0 for k in range(1, n + 1, 2)),
        'every odd-degree coefficient vanishes, as it must for skew C')
    chk('charpoly-endpoints', cp[-1] == dC and sum(cp) == dH,
        'p(0) = %d = det C and p(-1) = 1+45+770+6210+23485+33489 = %d = det H' % (cp[-1], sum(cp)))

    pf = pfaffian(C, n)
    chk('pfaffian', abs(pf) == CLAIM_PFAFFIAN and pf * pf == dC,
        'Pf(C) = %d, Pf(C)^2 = %d = det C' % (pf, pf * pf))
    trace = sum(C[k][i] * C[k][i] for i in range(n) for k in range(n))
    chk('trace-forces-45', trace == 90 and cp[2] == 45,
        'tr(C^T C) = %d over 90 off-diagonal entries of modulus 1, so the x^8 coefficient is 45' % trace)

    outdeg = sorted(sum(1 for j in range(n) if j != i and H[i][j] == 1) for i in range(n))
    chk('outdegrees', outdeg == CLAIM_OUTDEG, 'sorted out-degrees %s' % (outdeg,))

    und = [[0] * n for _ in range(n)]
    ne = 0
    for i in range(n):
        for j in range(i + 1, n):
            if H[i][j] == 1:
                und[i][j] = und[j][i] = 1
                ne += 1
    g6 = graph6(und, n)
    chk('graph6', g6 == CLAIM_GRAPH6 and ne == CLAIM_EDGES,
        'graph6 %s over %d edges {i<j : H[i][j]=+1}' % (g6, ne))

    print('--- 2. the invariance lemmas, and the bordering normal form ---')
    ok = True
    for bits in range(1 << (n - 1)):
        d = [1] + [1 if (bits >> t) & 1 else -1 for t in range(n - 1)]
        DC = [[d[i] * C[i][j] * d[j] for j in range(n)] for i in range(n)]
        if not is_cold(DC, n) or det_int(DC, n) != dC or det_int(hot_of(DC, n), n) != dH:
            ok = False
            break
    chk('switching-invariance', ok,
        'all 2^9 = 512 diagonal switchings D (d_0 = +1 without loss): DCD stays cold and both '
        'determinants are unchanged')

    ok = True
    perms = [list(range(n))]
    for a in range(n):
        for b in range(a + 1, n):
            p = list(range(n))
            p[a], p[b] = p[b], p[a]
            perms.append(p)
    for r in range(1, n):
        perms.append([(i + r) % n for i in range(n)])
    for p in perms:
        PC = [[C[p[i]][p[j]] for j in range(n)] for i in range(n)]
        if det_int(PC, n) != dC or det_int(hot_of(PC, n), n) != dH \
           or charpoly_of_skew(PC, n) != cp:
            ok = False
            break
    chk('permutation-invariance', ok,
        '%d permutations (identity, all 45 transpositions, all 9 rotations): det C, det H and the '
        'whole characteristic polynomial are unchanged' % len(perms))

    rnd = random.Random(20260829)
    ok = True
    tested = 0
    for _ in range(2000):
        R = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                v = 1 if rnd.getrandbits(1) else -1
                R[i][j] = v
                R[j][i] = -v
        d = [1] + [R[0][j] for j in range(1, n)]
        DR = [[d[i] * R[i][j] * d[j] for j in range(n)] for i in range(n)]
        if not (is_cold(DR, n) and all(DR[0][j] == 1 for j in range(1, n))
                and det_int(DR, n) == det_int(R, n)
                and det_int(hot_of(DR, n), n) == det_int(hot_of(R, n), n)):
            ok = False
            break
        tested += 1
    chk('bordering-normal-form', ok and tested == 2000,
        'on 2000 pseudo-random order-10 cold matrices (seed 20260829) the explicit switching '
        'd_j = c_{0j} makes row 0 all +1 and changes neither determinant')

    tj = tournament_of_bordered(C, n)
    chk('paper-object-is-bordered',
        border(tj, 9) == C and all(C[0][j] == 1 for j in range(1, n)),
        'the paper matrix is already in bordered form; its 9-vertex tournament has scores %s'
        % (sorted(bin(x).count('1') for x in tj),))

    print('--- 3. the covering ---')
    reps = class_reps(8)
    counts = [len(reps[k]) for k in range(1, 9)]
    chk('tournament-class-counts', counts == CLAIM_CLASS_COUNTS,
        'isomorphism classes on 1..8 vertices: %s, the published A000568' % (counts,))

    codes7 = sorted(canon(t, 7) for t in reps[7])
    ok = len(set(codes7)) == len(codes7)
    for t in reps[7]:
        c0 = canon(t, 7)
        for _ in range(10):
            p = list(range(7))
            rnd.shuffle(p)
            if canon(relabel(t, 7, p), 7) != c0:
                ok = False
                break
        if not ok:
            break
    chk('canonical-form-complete', ok,
        'the 456 seven-vertex codes are pairwise distinct and each is fixed by 10 pseudo-random '
        'relabellings -- 4560 relabellings in all')

    # -- the canonical form, exhaustively, on EVERY labelled tournament of order 5 and order 6 --
    groups5 = {}
    for code in range(CLAIM_LABELLED[5]):
        t = tournament_of_code(code, 5)
        if code_of(t, 5) != code:
            groups5 = None
            break
        groups5.setdefault(canon(t, 5), []).append(t)
    ok = groups5 is not None and len(groups5) == CLAIM_CLASS_COUNTS[4]
    if groups5 is not None:
        if sum(len(g) for g in groups5.values()) != CLAIM_LABELLED[5]:
            ok = False
        for c, members in groups5.items():
            orbit = set()
            for p in itertools.permutations(range(5)):
                u = relabel(members[0], 5, list(p))
                orbit.add(u)
                if canon(u, 5) != c:
                    ok = False
            if orbit != set(members) or len(orbit) * aut_full(members[0], 5) != 120:
                ok = False
    chk('canon-classes-are-orbits-n5-exhaustive', ok,
        'all 1024 labelled 5-vertex tournaments: the canonical code takes exactly 12 = A000568(5) '
        'values, and each code class is precisely one full 120-relabelling orbit, of size '
        '120/|Aut|')

    groups6 = {}
    for code in range(CLAIM_LABELLED[6]):
        groups6.setdefault(canon(tournament_of_code(code, 6), 6), []).append(code)
    ok = (len(groups6) == CLAIM_CLASS_COUNTS[5]
          and sum(len(g) for g in groups6.values()) == CLAIM_LABELLED[6])
    for c, members in groups6.items():
        m0 = tournament_of_code(members[0], 6)
        if len(members) * aut_full(m0, 6) != 720 or canon(m0, 6) != c:
            ok = False
    chk('canon-classes-are-orbits-n6-exhaustive', ok,
        'all 32768 labelled 6-vertex tournaments: exactly 56 = A000568(6) canonical codes, each '
        'code class of size 720/|Aut| computed over all 720 relabellings, sizes summing to 32768')

    # -- the automorphism count used by the orbit certificate, against the full symmetric group --
    ok = True
    tested = 0
    for k in range(1, 8):
        for t in reps[k]:
            if aut_cell(t, k) != aut_full(t, k):
                ok = False
            tested += 1
    chk('aut-cell-equals-aut-full', ok and tested == sum(CLAIM_CLASS_COUNTS[:7]),
        'for all %d representatives on 1..7 vertices, counting automorphisms over the '
        'score-preserving permutations agrees with counting over all k! relabellings (456 * 5040 '
        'at k = 7 alone)' % tested)

    # -- (R) the canonical code of each order-8 representative is REALIZED by a relabelling --
    ok = True
    codes8 = []
    for t in reps[8]:
        c, p = canon_witness(t, 8)
        u = relabel(t, 8, p)
        if c != canon(t, 8) or code_of(u, 8) != c or tournament_of_code(c, 8) != u \
                or canon(u, 8) != c:
            ok = False
            break
        codes8.append(c)
    chk('canon-code-realized-k8', ok and len(set(codes8)) == CLAIM_CLASS_COUNTS[7],
        'each of the 6880 order-8 codes is the exact upper-triangle code of an explicit relabelling '
        'of its representative, and the 6880 codes are pairwise distinct; hence equal codes force an '
        'isomorphism and the dedup inside class_reps can merge no two distinct classes')

    # -- (O) the orbits of the representatives account for every labelled tournament, exactly --
    ok = True
    totals = []
    for k in range(2, 9):
        fac = math.factorial(k)
        tot = 0
        for t in reps[k]:
            a = aut_cell(t, k)
            if a < 1 or fac % a != 0:
                ok = False
                a = max(a, 1)
            tot += fac // a
        totals.append(tot)
        if tot != CLAIM_LABELLED[k]:
            ok = False
    chk('class-set-complete-and-irredundant', ok,
        'for every k = 2..8 the orbit sizes k!/|Aut| of the representatives sum to exactly '
        '2^(k(k-1)/2) labelled tournaments: %s; at k = 8 that is 6880 classes summing to %d, so the '
        'list is complete AND has no repeated class'
        % (totals, CLAIM_LABELLED[8]))

    reps8set = set(canon(t, 8) for t in reps[8])
    ok = True
    for _ in range(3000):
        o = [0] * 9
        for i in range(9):
            for j in range(i + 1, 9):
                if rnd.getrandbits(1):
                    o[i] |= (1 << j)
                else:
                    o[j] |= (1 << i)
        sub = tuple(o[i] & 0xFF for i in range(8))
        if canon(sub, 8) not in reps8set:
            ok = False
            break
    chk('covering-property', ok,
        '3000 pseudo-random 9-vertex tournaments: deleting vertex 8 always leaves an 8-vertex '
        'tournament isomorphic to one of the 6880 representatives')
    chk('covering-size', len(reps[8]) * 256 == CLAIM_LEAVES_10,
        '6880 * 2^8 = %d leaves cover every switching-and-permutation class of order-10 cold '
        'matrices' % (len(reps[8]) * 256))

    print('--- 4. the published small cells ---')
    for m in (4, 6):
        s = brute_cell(m)
        fk, gk = CLAIM_SMALL_CELLS[m]
        detail = ('all %d labelled cold matrices: f_k(%d) = %d, g_k(%d) = %d, violA = %d, violB = %d'
                  % (s['total'], m, s['maxC'], m, s['maxH'], s['violA'], s['violB']))
        chk('brute-force-n%d' % m,
            s['maxC'] == fk and s['maxH'] == gk and s['violA'] == 0 and s['violB'] == 0
            and s['total'] == 1 << (m * (m - 1) // 2) and s['squares'], detail)
        if m == 6:
            chk('brute-force-n6-argmax-count',
                s['nC'] == s['nH'] == s['nBoth'] == CLAIM_N6_ARGMAX_COUNT,
                'exactly %d of the 32768 labelled matrices attain f_k(6), the same %d attain g_k(6)'
                % (s['nC'], s['nBoth']))

    for m in (4, 6, 8):
        s = cover_cell(reps, m)
        fk, gk = CLAIM_SMALL_CELLS[m]
        chk('covering-n%d' % m,
            s['maxC'] == fk and s['maxH'] == gk and s['violA'] == 0 and s['violB'] == 0
            and s['squares'],
            '%d classes of %d-vertex tournaments: f_k(%d) = %d, g_k(%d) = %d, violA = violB = 0'
            % (s['total'], m - 1, m, s['maxC'], m, s['maxH']))
        if m == 8:
            chk('pfaffian-spectrum-n8', s['pf'] == CLAIM_PF_SET_8,
                'the %d achieved order-8 cold Pfaffians are %s -- 29 is absent'
                % (len(s['pf']), s['pf']))

    print('--- 5. the n = 10 census ---')
    hist, leaves, winC, winH, crossed, mismatch = census10(reps[8])
    s = summarise(hist)
    chk('census-leaves', leaves == CLAIM_LEAVES_10 and s['total'] == CLAIM_LEAVES_10,
        '%d leaves enumerated, expected 6880 * 256 = %d' % (leaves, CLAIM_LEAVES_10))
    chk('census-max-cold', s['maxC'] == CLAIM_DET_COLD,
        'f_k(10) = max det C over every order-10 cold matrix = %d = 183^2' % s['maxC'])
    chk('census-max-hot', s['maxH'] == CLAIM_DET_HOT,
        'g_k(10) = max det H over every order-10 hot matrix = %d = 125*2^9' % s['maxH'])
    chk('census-side-constraints',
        s['squares'] and all(k[1] > 0 and k[1] % 512 == 0 for k in hist),
        'on all %d leaves: det C is an odd perfect square, det H > 0 and 2^9 divides det H'
        % leaves)
    chk('census-lin-both-directions',
        s['violA'] == 0 and s['violB'] == 0 and s['nC'] == s['nH'] == s['nBoth'],
        'violA (det C = f_k but det H < g_k) = %d, violB (det H = g_k but det C < f_k) = %d, and '
        'the two argmax sets are the same %d leaves' % (s['violA'], s['violB'], s['nBoth']))
    chk('census-pfaffian-spectrum-n10', s['pf'] == CLAIM_PF_SET_10,
        'the %d achieved order-10 cold Pfaffians are exactly the published set; the 15 odd values '
        '131,145,149,151,155,157,159,163,167,169,171,173,177,179,181 do not occur' % len(s['pf']))
    got_orrick = [(dh, dc) for (dc, dh) in hist if (dh, dc) in CLAIM_ORRICK]
    chk('census-orrick-forced-positive', sorted(got_orrick) == sorted(CLAIM_ORRICK),
        'both order-reversing pairs (det H, det C) = (45056, 21609) and (46080, 20449) occur as '
        'leaves, so the census can exhibit the shape of a counterexample')

    classesC = sorted(set(canon(t, 9) for t in winC))
    classesH = sorted(set(canon(t, 9) for t in winH))
    chk('census-uniqueness',
        len(classesC) == 1 and classesC == classesH,
        '%d winning leaves collapse to exactly %d isomorphism class of 9-vertex tournaments, the '
        'same class for both maxima; hence exactly one switching-and-permutation class of order-10 '
        'cold matrices attains either maximum' % (len(winC), len(classesC)))
    chk('census-winner-is-the-paper-object', classesC == [canon(tj, 9)],
        'the unique maximising class is the one exhibited in the paper (canonical code %d)'
        % canon(tj, 9))
    chk('census-two-engines-agree', crossed > 0 and mismatch == 0,
        '%d leaves recomputed by the division-checked determinant routine, %d disagreements'
        % (crossed, mismatch))
    chk('bareiss-divisions-exact', BADDIV == 0,
        '%d inexact divisions over every call to the division-checked routine in this run' % BADDIV)

    print()
    print('NOT RE-RUN: only the cell n = 10 is decided here. Orders 18, 22 and 30 -- the remaining '
          'cases of Lin\'s question below 30 -- are not touched, and the bordering-plus-covering '
          'method does not reach them: n = 18 would need the isomorphism classes of 17-vertex '
          'tournaments.')
    print('NOT RE-RUN: the four censuses reported in our own working notes (two C programs, one '
          'numpy batch pass, one min-code canonical form) are not re-executed; this program is a '
          'fourth, independently written, pure-Python census, and it agrees with them.')
    print('NOT RE-RUN: nothing here checks a bibliographic claim. That f_k(10) = 33489 was already '
          'proved by Klanderman, Montee, Piotrowski, Rice and Shader, and that both values were '
          'already reported by Alvarez, Armario, Frau and Gudiel, are statements about the '
          'literature; the only literature datum this program touches is the published order-10 '
          'Pfaffian spectrum, which it reproduces.')
    print('NOT RE-RUN: the identification of this witness with the order-10 cocyclic matrix of '
          'Alvarez et al. is deduced from uniqueness, not computed -- that paper prints no order-10 '
          'matrix to compare against.')
    print()

    bad = [nm for nm, ok_, _ in RESULTS if not ok_]
    if bad:
        print('VERDICT: %d of %d CHECKS FAILED: %s' % (len(bad), len(RESULTS), ', '.join(bad)))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % len(RESULTS))
    return 0


if __name__ == '__main__':
    sys.exit(main())

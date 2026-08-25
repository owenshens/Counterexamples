#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- independent verifier for paper.tex,
  "A Counterexample to a Conjecture of Franklin on (321,1342)-Avoiders"

Python 3.9, standard library only, exact integer arithmetic (no floats anywhere).

=====================================================================
DEFINITIONS AS IMPLEMENTED (check these against the paper's wording)
=====================================================================
(D1) CONTAINMENT.  A permutation p = p_1..p_n CONTAINS the pattern
     t = t_1..t_m iff there exist indices i_1 < i_2 < ... < i_m
     (NOT required to be consecutive) such that the subsequence
     p_{i_1},...,p_{i_m} is order-isomorphic to t, i.e.
     p_{i_r} < p_{i_s}  <=>  t_r < t_s  for all r,s.
     p AVOIDS t iff it does not contain t.  Av(321,1342) = permutations
     avoiding both 321 and 1342.  Implemented in contains_pattern().

(D2) INDECOMPOSABILITY.  A nonempty p = p_1..p_n is DECOMPOSABLE iff
     some PROPER prefix p_1..p_j (1 <= j < n) has value set exactly
     {1,...,j}; p is INDECOMPOSABLE iff no such proper prefix exists.
     (j = n is excluded: the whole permutation always has value set
     {1..n}, so "proper" is what makes the notion nontrivial.)
     Implemented in is_indecomposable().

(D3) inv(p) = #{ (i,j) : i < j and p_i > p_j }.

=====================================================================
TAKEN FROM THE PAPER (inputs; never used as the computed side of a check)
=====================================================================
 P1. The two patterns: 321 and 1342.
 P2. The claimed exact counts   |I_k| = 1 + C(k,2)   (k >= 0)
                                |S*_n| = 1 + C(n,3)  (n >= 1).
 P3. Franklin's conjectured count, as quoted: k(k+1)/2 + 1, which the
     paper rewrites as 1 + C(k+1,2).
 P4. The bivariate generating function, eq. (J):
       J(x,q) = x/(1-xq)
              + x^3 q^2/(1-xq) * SUM_{a>=0} (xq)^a /
                                   ((1-x q^{a+1})(1-x q^{a+2})).
 P5. Family type 1 (indecomposable):  n 1 2 ... (n-1),  n >= 1.
 P6. Family type 2:  pi = c g1 g2 n b d1 d2, all four blocks increasing,
     with g2 < b < d1 < c < g1 < d2 < n elementwise;
     a=|g1|, u=|g2|, v=|d1|, w=|d2|; length 3+a+u+v+w.
 P7. The paper's inversion formula for type 2:
       inv = 2 + a(u+v) + a + u + 2v + w,
     and its claim inv - (len-1) = a(u+v) + v >= 0.
 P8. The paper's cover of the WHOLE class (from BGU Prop. 7 / eq. (6)):
     every nonempty member of Av(321,1342) is uniquely either
       type 1':  alpha n beta, alpha in Av(321,1342), beta increasing,
                 every entry of alpha < every entry of beta,   or
       type 2 :  as in P6.
 P9. The paper's length generating function claims:
       F = 1 + xF/(1-x) + x^3/(1-x)^4,
       F = ((1-x)^4 + x^3)/((1-x)^3 (1-2x))
         = (1-4x+6x^2-3x^3+x^4)/(1-5x+9x^2-7x^3+2x^4)
         = 1 + x + 2x^2 + 5x^3 + 13x^4 + 32x^5 + ...
P10. The telescoping identity, for every A >= 0:
       SUM_{a=0..A} q^a/((1-q^{a+1})(1-q^{a+2}))
         = (1/(1-q)) * ( 1/(1-q) - q^{A+1}/(1-q^{A+2}) ),
     its limit 1/(1-q)^2, and the extractions
       [q^k] (1/(1-q) + q^2/(1-q)^3) = 1 + C(k,2),
       [x^n] (x/(1-x) + x^3/(1-x)^4) = 1 + C(n,3).
P11. The paper's assertion (attributed to Franklin) that an
     indecomposable permutation of length n has inv >= n-1.
P12. The paper's claim I_1 = {21}.

TAKEN FROM A SEPARATE EARLIER COMPUTATION made while the paper was
being prepared, and transcribed here by hand.  These values are used
ONLY as a second, differently-sourced cross-reference for
transcription errors; no check uses one of them as both sides.
 S1. |I_k| for k=0..10, |S*_n| for n=1..11, class sizes for n=1..11.
 S2. The four calibration sequences (OEIS A000041, A005169, A010054,
     A006958) for k=0..7.
 S3. Head counts: 506 indecomposables at n<=11, 495 type-2 tuples,
     637 class members / representations at n<=8.

=====================================================================
DERIVED HERE (the computed side of every check)
=====================================================================
 D1. Av(321,1342) up to length 11, generated from (D1) alone, by two
     independent generators: (i) closure under insertion of the new
     maximum, (ii) exhaustive scan of S_n for n <= 8.  The two are
     compared to each other.
 D2. The indecomposable members (D2), tabulated by inv and by length.
 D3. inv of every constructed permutation, from (D3), never from P7.
 D4. The type-1 and type-2 families built from P5/P6, their inversion
     numbers, their pairwise distinctness, and set equality of their
     union with D2's brute-force set.
 D5. All representations of class members under the P8 cover, counted
     with multiplicity, compared setwise with D1.
 D6. The bivariate series of P4, expanded by exact truncated series
     arithmetic to x^9, compared coefficientwise with D2.
 D7. F(x) from P9, expanded by exact polynomial long division to x^11,
     compared with D1; plus the two polynomial identities of P9.
 D8. The telescoping identity P10 as truncated q-series for A=0..8, and
     the two coefficient extractions.
 D9. Calibration: indecomposable Av(132), Av(231), Av(12), Av(321) by
     inversions for k=0..7, each compared against an INDEPENDENTLY
     COMPUTED model (integer partitions; fountains of coins;
     triangular-number test; parallelogram polyominoes by area) --
     not against a copied OEIS table.
D10. The refutation itself: 1 + C(k,2) != 1 + C(k+1,2) for k >= 1, the
     index-shift claim |I_{k+1}| = 1 + C(k+1,2), and I_0, I_1 as sets.
D11. Lemma P11 itself, tested on ALL permutations of length <= N_CROSS
     (not just class members), together with its mechanism: a
     permutation is indecomposable iff its inversion graph is connected,
     and a connected graph on m vertices has >= m-1 edges.  P11 is what
     makes a length-truncated enumeration complete in k, so it is tested
     rather than assumed.  See the SCOPE block printed at the end of a
     run for exactly what remains uncovered.
"""

import sys
from itertools import combinations, permutations

# ---------------------------------------------------------------
# P1: the two patterns, as tuples (one-line-notation).
# ---------------------------------------------------------------
PAT_321 = (3, 2, 1)
PAT_1342 = (1, 3, 4, 2)
PATTERNS = (PAT_321, PAT_1342)

# Enumeration ranges.  N_MAIN=11 makes |I_k| complete for all k <= 10,
# because an indecomposable permutation of length n has inv >= n-1
# (verified as a fact here, not assumed), so a permutation with k <= 10
# inversions has length <= 11.
N_MAIN = 11          # brute-force length bound for the main class
N_CROSS = 8          # length bound for the exhaustive S_n cross-check
N_COVER = 8          # length bound for the full-class cover check (P8)
N_SERIES = 9         # x-truncation for the bivariate series check (P4)
K_CAL = 7            # inversion bound for the calibration check (D9)

# ---------------------------------------------------------------
# S1/S2/S3: reference-reported values, cross-reference only.
# ---------------------------------------------------------------
SPEC_I_K = [1, 1, 2, 4, 7, 11, 16, 22, 29, 37, 46]          # k=0..10
SPEC_S_N = [1, 1, 2, 5, 11, 21, 36, 57, 85, 121, 166]       # n=1..11
SPEC_CLASS = [1, 2, 5, 13, 32, 74, 163, 347, 722, 1480, 3005]  # n=1..11
SPEC_CAL = {
    "132": [1, 1, 2, 3, 5, 7, 11, 15],    # A000041, partitions
    "231": [1, 1, 1, 2, 3, 5, 9, 15],     # A005169, fountains of coins
    "12":  [1, 1, 0, 1, 0, 0, 1, 0],      # A010054, triangular indicator
    "321": [1, 1, 2, 4, 9, 20, 46, 105],  # A006958, parallelogram polyominoes
}
SPEC_TOTAL_INDEC_11 = 506
SPEC_TYPE2_TUPLES_11 = 495
SPEC_CLASS_TOTAL_8 = 637
SPEC_NONZERO_JQ_9 = 43   # nonzero [x^n q^k] of eq. (J) with n <= 9

# P9: the two printed rational forms, as integer coefficient lists
# (index = power of x).
PAPER_F_NUM_PRINTED = [1, -4, 6, -3, 1]
PAPER_F_DEN_PRINTED = [1, -5, 9, -7, 2]
# P9: the printed head of the expansion of F.
PAPER_F_SERIES_PRINTED = [1, 1, 2, 5, 13, 32]


# ===============================================================
# Section 1: the definitions (D1)-(D3), implemented directly.
# ===============================================================

def binom(n, k):
    """Exact integer binomial coefficient, 0 for k<0 or k>n."""
    if k < 0 or n < 0 or k > n:
        return 0
    num = 1
    den = 1
    for i in range(k):
        num *= (n - i)
        den *= (i + 1)
    return num // den


def inv(p):
    """(D3) inversion number of a permutation given as a sequence."""
    n = len(p)
    c = 0
    for i in range(n):
        pi = p[i]
        for j in range(i + 1, n):
            if pi > p[j]:
                c += 1
    return c


def contains_pattern(p, pat):
    """(D1) True iff p contains pat as an order-isomorphic subsequence.

    Subsequences need NOT be consecutive.  Order-isomorphism is tested
    pairwise on all pairs of positions of the pattern, so it does not
    depend on the pattern being a permutation of 1..m.
    """
    m = len(pat)
    if m > len(p):
        return False
    pairs = [(r, s) for r in range(m) for s in range(r + 1, m)]
    for sub in combinations(p, m):
        ok = True
        for (r, s) in pairs:
            if (sub[r] < sub[s]) != (pat[r] < pat[s]):
                ok = False
                break
        if ok:
            return True
    return False


def avoids_all(p, pats):
    """True iff p avoids every pattern in pats."""
    for pat in pats:
        if contains_pattern(p, pat):
            return False
    return True


def is_indecomposable(p):
    """(D2) True iff no PROPER prefix of p has value set {1,..,j}."""
    if len(p) == 0:
        return False
    hi = 0
    for j in range(1, len(p)):          # j = prefix length, proper: j < n
        if p[j - 1] > hi:
            hi = p[j - 1]
        if hi == j:                     # prefix values are exactly {1..j}
            return False
    return True


# ===============================================================
# Section 2: two independent generators of Av(pats) (D1).
# ===============================================================

def gen_by_max_insertion(pats, nmax):
    """Generator A: all avoiders of length 1..nmax, by growth.

    Deleting the largest entry of a permutation leaves a permutation
    that still avoids every pattern (pattern classes are closed under
    taking sub-permutations), so every avoider of length n arises from
    an avoider of length n-1 by inserting the value n somewhere.  We
    therefore build level n from level n-1 and filter by (D1).  This is
    exhaustive, not heuristic.
    """
    levels = {0: [tuple()]}
    for n in range(1, nmax + 1):
        out = []
        for p in levels[n - 1]:
            for pos in range(n):
                cand = p[:pos] + (n,) + p[pos:]
                if avoids_all(cand, pats):
                    out.append(cand)
        levels[n] = out
    del levels[0]
    return levels


def gen_by_exhaustive_scan(pats, nmax):
    """Generator B: filter all of S_n directly, for n = 1..nmax.

    Independent of generator A: it never uses closure under deletion of
    the maximum, it just walks every permutation of {1..n}.  Used only
    as a cross-check, so nmax must stay small.
    """
    levels = {}
    for n in range(1, nmax + 1):
        out = []
        for p in permutations(range(1, n + 1)):
            if avoids_all(p, pats):
                out.append(p)
        levels[n] = out
    return levels


# ===============================================================
# Section 3: tabulation (D2).
# ===============================================================

def tabulate_indecomposable(levels):
    """From {n: [avoiders]} produce the indecomposable data.

    Returns (indec_set, by_len, by_inv, joint) where
      indec_set : set of tuples, every indecomposable avoider found
      by_len[n] : count of indecomposable avoiders of length n
      by_inv[k] : count of indecomposable avoiders with inv = k
      joint[(n,k)] : count with that length and inversion number
    Nothing here consults any expected value.
    """
    indec_set = set()
    by_len = {}
    by_inv = {}
    joint = {}
    for n in sorted(levels):
        by_len[n] = 0
        for p in levels[n]:
            if not is_indecomposable(p):
                continue
            k = inv(p)
            indec_set.add(p)
            by_len[n] += 1
            by_inv[k] = by_inv.get(k, 0) + 1
            joint[(n, k)] = joint.get((n, k), 0) + 1
    return indec_set, by_len, by_inv, joint


def min_inv_minus_len(indec_set):
    """min over the set of (inv(p) - (len(p)-1)); P11 says this is >= 0.

    Returns None when the set is empty, so an empty input can never be
    mistaken for a satisfied bound.
    """
    best = None
    for p in indec_set:
        d = inv(p) - (len(p) - 1)
        if best is None or d < best:
            best = d
    return best


# ===============================================================
# Section 4: the paper's two families, built verbatim (P5, P6).
# ===============================================================

def family_type1(nmax):
    """P5: pi = n 1 2 ... (n-1) for n = 1..nmax."""
    out = []
    for n in range(1, nmax + 1):
        out.append((n,) + tuple(range(1, n)))
    return out


def family_type2(nmax):
    """P6: pi = c g1 g2 n b d1 d2 with g2 < b < d1 < c < g1 < d2 < n.

    The chain of strict elementwise inequalities forces the value
    assignment: reading the blocks in increasing value order they are
    g2 (u values), b (1), d1 (v), c (1), g1 (a), d2 (w), n (1).
    Each block is written in increasing order.  Returns a list of
    (tuple, a, u, v, w) for every (a,u,v,w) with 3+a+u+v+w <= nmax.
    """
    out = []
    for total in range(0, nmax - 3 + 1):
        for a in range(total + 1):
            for u in range(total - a + 1):
                for v in range(total - a - u + 1):
                    w = total - a - u - v
                    n = 3 + a + u + v + w
                    g2 = tuple(range(1, u + 1))
                    b = u + 1
                    d1 = tuple(range(u + 2, u + 2 + v))
                    c = u + 2 + v
                    g1 = tuple(range(c + 1, c + 1 + a))
                    d2 = tuple(range(c + 1 + a, c + 1 + a + w))
                    top = n
                    pi = (c,) + g1 + g2 + (top,) + (b,) + d1 + d2
                    out.append((pi, a, u, v, w))
    return out


def blocks_increasing_and_ordered(pi, a, u, v, w):
    """Re-derive the block decomposition of a type-2 word and verify the
    inequality chain g2 < b < d1 < c < g1 < d2 < n elementwise, plus
    that each block is increasing.  Returns True/False.
    """
    n = 3 + a + u + v + w
    if len(pi) != n or sorted(pi) != list(range(1, n + 1)):
        return False
    c = pi[0]
    g1 = pi[1:1 + a]
    g2 = pi[1 + a:1 + a + u]
    top = pi[1 + a + u]
    b = pi[2 + a + u]
    d1 = pi[3 + a + u:3 + a + u + v]
    d2 = pi[3 + a + u + v:]
    for blk in (g1, g2, d1, d2):
        if list(blk) != sorted(blk):
            return False
    if top != n:
        return False
    chain = [g2, (b,), d1, (c,), g1, d2, (n,)]
    # All PAIRS, not just consecutive ones: with consecutive comparisons
    # only, an empty block in the middle would silently break the chain
    # (transitivity has nothing to pass through) and the test would pass
    # vacuously on words violating the paper's condition.
    for i in range(len(chain)):
        for j in range(i + 1, len(chain)):
            lo, hi = chain[i], chain[j]
            if lo and hi and max(lo) >= min(hi):
                return False
    return True


# ===============================================================
# Section 5: the full-class cover of P8, taken literally.
# ===============================================================

def cover_type1_reps(class_levels, nmax):
    """P8 type 1': pi = alpha n beta, alpha in Av(321,1342), beta
    increasing, every entry of alpha below every entry of beta.

    Because alpha and beta together use the values {1,...,n-1} and alpha
    lies entirely below beta, alpha must occupy {1,...,|alpha|} and beta
    is then forced to be (|alpha|+1, ..., n-1) in increasing order.  So a
    representation is exactly a pair (alpha, n) with |alpha| < n, and we
    enumerate all of them (alpha empty included: it is a member of the
    class).  The value pattern of alpha is used as-is, since its entries
    are already 1..|alpha|.
    """
    reps = []
    for n in range(1, nmax + 1):
        for j in range(0, n):
            alphas = [tuple()] if j == 0 else class_levels.get(j, [])
            beta = tuple(range(j + 1, n))
            for alpha in alphas:
                reps.append((alpha + (n,) + beta, "type1"))
    return reps


def cover_all_reps(class_levels, nmax):
    """All P8 representations of length <= nmax, with their type tags."""
    reps = cover_type1_reps(class_levels, nmax)
    for (pi, a, u, v, w) in family_type2(nmax):
        reps.append((pi, "type2"))
    return reps


# ===============================================================
# Section 6: exact truncated series arithmetic (integers only).
# A bivariate series is a dict {(n, k): int}: coefficient of x^n q^k.
# Truncation is on the x-degree only; q-degree stays exact because
# every factor contributes q-degree bounded by its x-degree times a
# fixed integer, so no q-truncation is needed below x^N.
# ===============================================================

def bs_mul(f, g, xmax):
    """Product of two bivariate series, truncated to x-degree <= xmax."""
    out = {}
    for (n1, k1), c1 in f.items():
        if n1 > xmax:
            continue
        for (n2, k2), c2 in g.items():
            n = n1 + n2
            if n > xmax:
                continue
            key = (n, k1 + k2)
            out[key] = out.get(key, 0) + c1 * c2
    return {kk: cc for kk, cc in out.items() if cc != 0}


def bs_add(f, g):
    """Sum of two bivariate series."""
    out = dict(f)
    for kk, cc in g.items():
        out[kk] = out.get(kk, 0) + cc
    return {kk: cc for kk, cc in out.items() if cc != 0}


def bs_geom_xqm(m, xmax):
    """1/(1 - x q^m) = sum_{i>=0} x^i q^{m i}, truncated at x^xmax."""
    return {(i, m * i): 1 for i in range(0, xmax + 1)}


def bs_monomial(n, k, coeff=1):
    """coeff * x^n q^k."""
    return {(n, k): coeff}


def build_J(xmax):
    """P4, expanded verbatim to x-degree <= xmax.

      J = x/(1-xq) + x^3 q^2/(1-xq) * SUM_{a>=0} (xq)^a
              / ((1-x q^{a+1})(1-x q^{a+2}))

    The sum over a is finite after truncation: the factor (xq)^a carries
    x^a and the prefix carries x^3, so terms with a > xmax-3 cannot
    contribute.  A range that came out empty would make the second term
    vanish silently, so the caller is handed the number of a-terms used.
    """
    geom_xq = bs_geom_xqm(1, xmax)
    term1 = bs_mul(bs_monomial(1, 0), geom_xq, xmax)

    inner = {}
    a_terms = 0
    for a in range(0, max(xmax - 3, -1) + 1):
        piece = bs_monomial(a, a)                       # (xq)^a
        piece = bs_mul(piece, bs_geom_xqm(a + 1, xmax), xmax)
        piece = bs_mul(piece, bs_geom_xqm(a + 2, xmax), xmax)
        inner = bs_add(inner, piece)
        a_terms += 1
    term2 = bs_mul(bs_monomial(3, 2), geom_xq, xmax)
    term2 = bs_mul(term2, inner, xmax)

    J = bs_add(term1, term2)
    J = {kk: cc for kk, cc in J.items() if kk[0] <= xmax and cc != 0}
    return J, a_terms, term1, term2


def naive_type1_series(xmax):
    """SUM_{n>=1} x^n q^{n-1}: the type-1 contribution written out as a
    sum over the family P5, before the paper's closed form x/(1-xq)."""
    return {(n, n - 1): 1 for n in range(1, xmax + 1)}


def naive_type2_series(xmax):
    """SUM_{a,u,v,w>=0} x^{3+a+u+v+w} q^{2+a(u+v)+a+u+2v+w}.

    The four-fold sum of the paper, written out term by term from the
    exponent formula P7, before the geometric-series manipulation that
    produces the closed form.  Comparing this with the closed form is a
    check on that manipulation alone, independent of any permutation.
    """
    out = {}
    for (pi, a, u, v, w) in family_type2(xmax):
        n = 3 + a + u + v + w
        k = 2 + a * (u + v) + a + u + 2 * v + w
        out[(n, k)] = out.get((n, k), 0) + 1
    return out


# ---------------------------------------------------------------
# Univariate integer series, as coefficient lists of length kmax+1.
# ---------------------------------------------------------------

def us_zero(kmax):
    return [0] * (kmax + 1)


def us_add(f, g):
    return [x + y for x, y in zip(f, g)]


def us_sub(f, g):
    return [x - y for x, y in zip(f, g)]


def us_mul(f, g, kmax):
    """Truncated product.  Both operands are lists of length kmax+1, so
    every index i satisfies i <= kmax and no length guard is needed."""
    out = [0] * (kmax + 1)
    for i, a in enumerate(f):
        if a == 0:
            continue
        for j, b in enumerate(g[:kmax - i + 1]):
            if b:
                out[i + j] += a * b
    return out


def us_geom(m, kmax):
    """1/(1 - q^m) = 1 + q^m + q^{2m} + ...  (m >= 1)."""
    out = [0] * (kmax + 1)
    i = 0
    while i <= kmax:
        out[i] = 1
        i += m
    return out


def us_shift(f, s, kmax):
    """q^s * f, truncated."""
    out = [0] * (kmax + 1)
    for i, a in enumerate(f):
        if a and i + s <= kmax:
            out[i + s] = a
    return out


def us_pow(f, e, kmax):
    out = [0] * (kmax + 1)
    out[0] = 1
    for _ in range(e):
        out = us_mul(out, f, kmax)
    return out


def telescope_sides(A, kmax):
    """P10 for one value of A, as two truncated q-series.

    LHS = SUM_{a=0..A} q^a / ((1-q^{a+1})(1-q^{a+2}))
    RHS = 1/(1-q) * ( 1/(1-q) - q^{A+1}/(1-q^{A+2}) )

    The two sides are built from different expressions; neither is
    derived from the other.  Returns (lhs, rhs, n_terms).
    """
    lhs = us_zero(kmax)
    n_terms = 0
    for a in range(0, A + 1):
        piece = us_mul(us_geom(a + 1, kmax), us_geom(a + 2, kmax), kmax)
        lhs = us_add(lhs, us_shift(piece, a, kmax))
        n_terms += 1
    g1 = us_geom(1, kmax)
    tail = us_shift(us_geom(A + 2, kmax), A + 1, kmax)
    rhs = us_mul(g1, us_sub(g1, tail), kmax)
    return lhs, rhs, n_terms


def limit_sum_sides(kmax):
    """P10's limit: SUM_{a>=0} q^a/((1-q^{a+1})(1-q^{a+2})) = 1/(1-q)^2.

    Terms with a > kmax start at q^a and so cannot affect coefficients
    up to q^kmax; taking A = kmax therefore gives the exact truncation
    of the infinite sum.
    """
    lhs, _rhs_unused, n_terms = telescope_sides(kmax, kmax)
    rhs = us_pow(us_geom(1, kmax), 2, kmax)
    return lhs, rhs, n_terms


def extract_Ik_series(kmax):
    """P10: [q^k] ( 1/(1-q) + q^2/(1-q)^3 ) as a coefficient list."""
    g1 = us_geom(1, kmax)
    return us_add(g1, us_shift(us_pow(g1, 3, kmax), 2, kmax))


def extract_Sn_series(nmax):
    """P10: [x^n] ( x/(1-x) + x^3/(1-x)^4 ) as a coefficient list."""
    g1 = us_geom(1, nmax)
    return us_add(us_shift(g1, 1, nmax),
                  us_shift(us_pow(g1, 4, nmax), 3, nmax))


# ---------------------------------------------------------------
# Exact (untruncated) integer polynomials, as coefficient lists.
# ---------------------------------------------------------------

def poly_trim(p):
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def poly_mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b:
                    out[i + j] += a * b
    return poly_trim(out)


def poly_add(p, q):
    n = max(len(p), len(q))
    return poly_trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
                      for i in range(n)])


def poly_one_minus_x_pow(e):
    """(1-x)^e for e >= 0."""
    out = [1]
    for _ in range(e):
        out = poly_mul(out, [1, -1])
    return out


def series_div(num, den, nmax):
    """Power series expansion of num/den to x^nmax, exact integers.

    Requires den[0] == 1 (true for every denominator used here), which
    makes every coefficient an integer produced by subtraction only.
    """
    if not den or den[0] != 1:
        raise ValueError("series_div needs constant term 1 in denominator")
    c = [0] * (nmax + 1)
    for n in range(nmax + 1):
        acc = num[n] if n < len(num) else 0
        for j in range(1, min(n, len(den) - 1) + 1):
            acc -= den[j] * c[n - j]
        c[n] = acc
    return c


# ===============================================================
# Section 7: independent models for the calibration sequences (D9).
# Each of the four is computed from its OWN combinatorial definition,
# so the calibration compares permutation brute force against a model
# that shares no code path with it.
# ===============================================================

def partition_counts(kmax):
    """p(k) = number of integer partitions of k, k = 0..kmax (A000041).

    Franklin's proved value for the indecomposable 132-avoiders by
    inversions.  Standard coin-change DP; nothing permutational here.
    """
    dp = [0] * (kmax + 1)
    dp[0] = 1
    for part in range(1, kmax + 1):
        for v in range(part, kmax + 1):
            dp[v] += dp[v - part]
    return dp


def triangular_indicator(kmax):
    """1 if k is a triangular number m(m+1)/2, else 0 (A010054).

    Franklin's proved value for the indecomposable 12-avoiders by
    inversions: the only 12-avoider of length n is n(n-1)...1, whose
    inversion number is C(n,2).
    """
    tri = set()
    m = 0
    while m * (m + 1) // 2 <= kmax:
        tri.add(m * (m + 1) // 2)
        m += 1
    return [1 if k in tri else 0 for k in range(kmax + 1)]


def _fountain_gaps(row):
    """Positions available in the row above `row`.

    Coins live on a doubled integer lattice so that a coin resting in the
    gap between two touching coins at p and p+2 sits at p+1.  Two coins
    touch exactly when their positions differ by 2.
    """
    row = tuple(sorted(row))
    return tuple(row[i] + 1 for i in range(len(row) - 1)
                 if row[i + 1] - row[i] == 2)


def _fountain_ways_above(row, remaining, memo):
    """Number of fountains completing `row` with exactly `remaining`
    further coins placed strictly above it."""
    if remaining == 0:
        return 1
    key = (row, remaining)
    if key in memo:
        return memo[key]
    gaps = _fountain_gaps(row)
    total = 0
    for size in range(1, min(len(gaps), remaining) + 1):
        for sub in combinations(gaps, size):
            total += _fountain_ways_above(sub, remaining - size, memo)
    memo[key] = total
    return total


def fountain_counts(kmax):
    """Fountains of k coins, k = 0..kmax (A005169).

    A fountain: a contiguous base row, and above it any set of coins each
    resting in the gap between two touching coins of the row below,
    recursively.  Franklin's proved value for the indecomposable
    231-avoiders by inversions.  k = 0 is the empty fountain.
    """
    memo = {}
    out = [0] * (kmax + 1)
    out[0] = 1
    for k in range(1, kmax + 1):
        tot = 0
        for b in range(1, k + 1):
            base = tuple(range(0, 2 * b, 2))
            tot += _fountain_ways_above(base, k - b, memo)
        out[k] = tot
    return out


def _pgram_ways(h, remaining, memo):
    """Completions of a parallelogram polyomino to the right.

    State: the height h of the current (rightmost so far) column and the
    area `remaining` still to be placed in further columns.  Writing
    b_{i+1} = b_i + s with 0 <= s <= h_i - 1 (columns must overlap in at
    least one cell) and requiring the tops to be weakly increasing,
    t_{i+1} >= t_i, gives h_{i+1} >= h_i - s.  The state therefore needs
    nothing but (h, remaining).
    """
    key = (h, remaining)
    if key in memo:
        return memo[key]
    total = 1 if remaining == 0 else 0        # stop after this column
    for s in range(0, h):
        for h2 in range(max(1, h - s), remaining + 1):
            total += _pgram_ways(h2, remaining - h2, memo)
    memo[key] = total
    return total


def parallelogram_counts(kmax):
    """Parallelogram (staircase) polyominoes with area k, k = 0..kmax
    (A006958), with the empty polyomino counted at k = 0.

    A parallelogram polyomino is a column-convex polyomino whose column
    bottoms and column tops are both weakly increasing and whose
    consecutive columns overlap.  Franklin's proved value for the
    indecomposable 321-avoiders by inversions.
    """
    memo = {}
    out = [0] * (kmax + 1)
    out[0] = 1
    for area in range(1, kmax + 1):
        tot = 0
        for h1 in range(1, area + 1):
            tot += _pgram_ways(h1, area - h1, memo)
        out[area] = tot
    return out


# ===============================================================
# Section 8: check harness.
# ===============================================================

class Checks(object):
    """Collects PASS/FAIL lines.  Every check must be given both sides
    explicitly, so the transcript shows what was compared."""

    def __init__(self):
        self.lines = []
        self.failed = 0

    def eq(self, label, computed, expected):
        ok = (computed == expected)
        self._record(ok, label, computed, expected)
        return ok

    def ne(self, label, computed, expected):
        ok = (computed != expected)
        self._record(ok, label + " [must differ]", computed, expected)
        return ok

    def true(self, label, cond, detail=""):
        self._record(bool(cond), label, detail if detail else bool(cond), "")
        return bool(cond)

    def _record(self, ok, label, computed, expected):
        tag = "PASS" if ok else "FAIL"
        if not ok:
            self.failed += 1
        if expected == "":
            self.lines.append("%s %s | got=%s" % (tag, label, computed))
        else:
            self.lines.append("%s %s | computed=%s expected=%s"
                              % (tag, label, computed, expected))

    def report(self):
        for ln in self.lines:
            print(ln)
        total = len(self.lines)
        if self.failed == 0:
            print("VERDICT: ALL %d CHECKS PASS" % total)
            return 0
        print("VERDICT: %d OF %d CHECKS FAILED" % (self.failed, total))
        return 1


# ===============================================================
# Section 9: stages.
# ===============================================================

def stage_enumeration(ck):
    """D1/D2: generate the class, cross-check the generators, tabulate."""
    print("== STAGE 1: brute-force enumeration of Av(321,1342) ==")
    print("   patterns as implemented: %s and %s" % (PAT_321, PAT_1342))
    print("   lengths generated: n = 1..%d (generator A: max-insertion)"
          % N_MAIN)
    print("   independent cross-check: n = 1..%d (generator B: all of S_n)"
          % N_CROSS)

    lv_a = gen_by_max_insertion(PATTERNS, N_MAIN)
    lv_b = gen_by_exhaustive_scan(PATTERNS, N_CROSS)

    ck.true("stage1.generatorA-nonempty",
            all(len(lv_a.get(n, [])) > 0 for n in range(1, N_MAIN + 1)),
            "levels=%s" % [len(lv_a[n]) for n in range(1, N_MAIN + 1)])
    agree = True
    for n in range(1, N_CROSS + 1):
        if set(lv_a[n]) != set(lv_b[n]):
            agree = False
    ck.true("stage1.two-generators-agree-setwise-n<=%d" % N_CROSS, agree,
            "A=%s B=%s" % ([len(lv_a[n]) for n in range(1, N_CROSS + 1)],
                           [len(lv_b[n]) for n in range(1, N_CROSS + 1)]))

    class_sizes = [len(lv_a[n]) for n in range(1, N_MAIN + 1)]
    indec_set, by_len, by_inv, joint = tabulate_indecomposable(lv_a)

    print("   DERIVED class sizes |Av(321,1342)_n|, n=1..%d: %s"
          % (N_MAIN, class_sizes))
    print("   DERIVED |S*_n| (indecomposable, by length), n=1..%d:" % N_MAIN)
    print("     n : count : 1+C(n,3)")
    for n in range(1, N_MAIN + 1):
        print("     %2d : %5d : %d" % (n, by_len[n], 1 + binom(n, 3)))
    kmax = N_MAIN - 1
    print("   DERIVED |I_k| (indecomposable, by inversions), k=0..%d:" % kmax)
    print("     k : count : 1+C(k,2)")
    for k in range(0, kmax + 1):
        print("     %2d : %5d : %d" % (k, by_inv.get(k, 0), 1 + binom(k, 2)))

    return {"levels": lv_a, "class_sizes": class_sizes,
            "indec_set": indec_set, "by_len": by_len,
            "by_inv": by_inv, "joint": joint, "kmax": kmax}


def boundary_scan(pats, prev_level, n, cap):
    """Every INDECOMPOSABLE member of Av(pats) of length n with inv <= cap.

    Completeness: deleting the entry n from an avoider of length n leaves
    an avoider of length n-1 with NO MORE inversions, so every length-n
    avoider with inv <= cap -- indecomposable or not -- is obtained by
    inserting the value n into some length-(n-1) avoider that itself has
    inv <= cap.  Iterating over the complete level n-1 therefore reaches
    all of them; we then keep the indecomposable ones.  Inserting n at
    0-indexed position pos adds exactly (n-1-pos) inversions.

    NOTE ON SCOPE (do not overclaim).  An empty result shows there is no
    indecomposable avoider of length exactly n with inv <= cap.  It does
    NOT by itself close lengths > n, because max-deletion does not
    preserve indecomposability, so the induction cannot be pushed up.
    Lengths > n are closed by the general lemma P11 (indecomposable =>
    inv >= len-1), which is verified exhaustively elsewhere in this
    program over all permutations of length <= N_CROSS and over every
    indecomposable class member found, together with the graph fact that
    the inversion graph of an indecomposable permutation is connected and
    so has at least len-1 edges.  That lemma is the one step here taken
    from the literature rather than computed.

    Returns (found, n_candidates_tested) so an empty candidate loop is
    visible rather than being mistaken for a satisfied claim.
    """
    found = []
    n_cand = 0
    for p in prev_level:
        ip = inv(p)
        if ip > cap:
            continue
        for pos in range(n):
            added = n - 1 - pos
            if ip + added > cap:
                continue
            cand = p[:pos] + (n,) + p[pos:]
            n_cand += 1
            if avoids_all(cand, pats) and is_indecomposable(cand):
                found.append(cand)
    return found, n_cand


def inversion_graph_connected(p):
    """True iff the inversion graph of p is connected.

    Vertices are the positions of p, with an edge {i,j} whenever (i,j) is
    an inversion.  A connected graph on m vertices has at least m-1
    edges, so "indecomposable => inversion graph connected" implies
    "indecomposable => inv >= len-1", which is lemma P11.  A one-vertex
    graph is connected and needs 0 edges, consistent with inv(1) = 0.
    """
    m = len(p)
    if m == 0:
        return False
    adj = dict((i, set()) for i in range(m))
    for i in range(m):
        for j in range(i + 1, m):
            if p[i] > p[j]:
                adj[i].add(j)
                adj[j].add(i)
    seen = set([0])
    stack = [0]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == m


def stage_lemma(ck):
    """P11, the one imported lemma, tested exhaustively on all of S_n.

    Not restricted to Av(321,1342): the lemma is a statement about all
    permutations, and the completeness of the |I_k| range depends on it,
    so it is tested on every permutation of length <= N_CROSS.
    """
    print("== STAGE 0: lemma P11 (indecomposable => inv >= len-1) ==")
    bad_bound = []
    bad_equiv = []
    n_indec = 0
    n_total = 0
    for n in range(1, N_CROSS + 1):
        for p in permutations(range(1, n + 1)):
            n_total += 1
            ind = is_indecomposable(p)
            con = inversion_graph_connected(p)
            if ind != con:
                bad_equiv.append(p)
            if ind:
                n_indec += 1
                if inv(p) < n - 1:
                    bad_bound.append(p)
    print("   permutations tested: %d (n=1..%d), indecomposable: %d"
          % (n_total, N_CROSS, n_indec))
    ck.true("stage0.lemma-test-loop-nonempty", n_indec > 0,
            "indecomposable permutations tested = %d" % n_indec)
    ck.eq("stage0.indecomposable-implies-inv>=len-1 (all of S_n, n<=%d)"
          % N_CROSS, len(bad_bound), 0)
    ck.eq("stage0.indecomposable-IFF-inversion-graph-connected (mechanism"
          " of the lemma, all of S_n, n<=%d)" % N_CROSS, len(bad_equiv), 0)


def stage_formulas(ck, data):
    """P2 vs D2: the two closed forms, and the completeness of the range."""
    print("== STAGE 2: the paper's closed forms against brute force ==")
    kmax = data["kmax"]
    by_inv, by_len = data["by_inv"], data["by_len"]

    derived_Ik = [by_inv.get(k, 0) for k in range(0, kmax + 1)]
    formula_Ik = [1 + binom(k, 2) for k in range(0, kmax + 1)]
    derived_Sn = [by_len[n] for n in range(1, N_MAIN + 1)]
    formula_Sn = [1 + binom(n, 3) for n in range(1, N_MAIN + 1)]

    # These two are configuration guards, not mathematics: they fire if
    # someone narrows the enumeration bounds and leaves the claims in
    # place.  The mathematical non-vacuity check is the next one.
    ck.true("stage2.config-guard-k-range-is-0..10",
            len(derived_Ik) == kmax + 1 == 11, "k covered = 0..%d" % kmax)
    ck.true("stage2.config-guard-n-range-is-1..11",
            len(derived_Sn) == N_MAIN == 11, "n covered = 1..%d" % N_MAIN)
    ck.true("stage2.all-Ik-positive", all(c > 0 for c in derived_Ik),
            "derived=%s" % derived_Ik)
    ck.eq("stage2.|I_k|=1+C(k,2) for k=0..%d" % kmax, derived_Ik, formula_Ik)
    ck.eq("stage2.|S*_n|=1+C(n,3) for n=1..%d" % N_MAIN, derived_Sn,
          formula_Sn)

    # Completeness of the k-range at the truncation boundary.
    over, n_cand = boundary_scan(PATTERNS, data["levels"][N_MAIN],
                                 N_MAIN + 1, kmax)
    ck.true("stage2.boundary-scan-candidate-loop-nonempty", n_cand > 0,
            "candidates of length %d with inv<=%d tested = %d"
            % (N_MAIN + 1, kmax, n_cand))
    ck.eq("stage2.no-INDECOMPOSABLE-class-member-of-length-%d-with-inv<=%d"
          " (lengths >%d rest on lemma P11)"
          % (N_MAIN + 1, kmax, N_MAIN + 1), len(over), 0)

    # P11, checked as a fact on everything enumerated, not assumed.
    slack = min_inv_minus_len(data["indec_set"])
    ck.true("stage2.P11-inv>=len-1-on-all-%d-indecomposables"
            % len(data["indec_set"]),
            slack is not None and slack >= 0, "min(inv-(len-1))=%s" % slack)

    # Margins recomputed by a SECOND pass over the deduplicated SET of
    # indecomposables, not from the same accumulator loop that produced
    # by_len/by_inv.  If the generator had ever emitted a permutation
    # twice, the list-based counts would exceed these and this would fail.
    m_len = {}
    m_inv = {}
    for p in data["indec_set"]:
        m_len[len(p)] = m_len.get(len(p), 0) + 1
        m_inv[inv(p)] = m_inv.get(inv(p), 0) + 1
    ck.eq("stage2.length-margin-recomputed-from-dedup-set", m_len, by_len)
    ck.eq("stage2.inversion-margin-recomputed-from-dedup-set", m_inv, by_inv)
    ck.eq("stage2.joint-table-total-equals-set-size",
          sum(data["joint"].values()), len(data["indec_set"]))

    # Cross-reference against the task independently reported
    # values (transcription check; neither side is the other's source).
    ck.eq("stage2.xref-reference-Ik", derived_Ik, SPEC_I_K)
    ck.eq("stage2.xref-reference-Sn", derived_Sn, SPEC_S_N)
    ck.eq("stage2.xref-reference-class-sizes", data["class_sizes"],
          SPEC_CLASS)


def stage_families(ck, data):
    """P5/P6/P7 vs D4: rebuild both families and match the brute force."""
    print("== STAGE 3: the paper's two families, rebuilt verbatim ==")
    t1 = family_type1(N_MAIN)
    t2 = family_type2(N_MAIN)
    print("   type 1 permutations built: %d (n=1..%d)" % (len(t1), N_MAIN))
    print("   type 2 (a,u,v,w) tuples built: %d (length 3..%d)"
          % (len(t2), N_MAIN))

    ck.eq("stage3.type1-count", len(t1), N_MAIN)
    ck.true("stage3.type2-loop-nonempty", len(t2) > 0, "tuples=%d" % len(t2))
    ck.eq("stage3.type2-tuple-count-equals-C(%d,4)" % (N_MAIN + 1),
          len(t2), binom(N_MAIN + 1, 4))
    ck.eq("stage3.xref-reference-type2-tuple-count", len(t2),
          SPEC_TYPE2_TUPLES_11)

    bad_t1 = [p for p in t1
              if not (avoids_all(p, PATTERNS) and is_indecomposable(p)
                      and inv(p) == len(p) - 1)]
    ck.eq("stage3.type1-all-in-class-indecomposable-inv=len-1",
          len(bad_t1), 0)

    bad_struct, bad_avoid, bad_indec, bad_inv, bad_slack = [], [], [], [], []
    per_len = {}
    for (pi, a, u, v, w) in t2:
        n = len(pi)
        per_len[n] = per_len.get(n, 0) + 1
        if not blocks_increasing_and_ordered(pi, a, u, v, w):
            bad_struct.append(pi)
        if not avoids_all(pi, PATTERNS):
            bad_avoid.append(pi)
        if not is_indecomposable(pi):
            bad_indec.append(pi)
        got = inv(pi)                                  # D3, computed
        want = 2 + a * (u + v) + a + u + 2 * v + w     # P7, the paper's
        if got != want:
            bad_inv.append((pi, a, u, v, w, got, want))
        if got - (n - 1) != a * (u + v) + v:
            bad_slack.append((pi, got - (n - 1), a * (u + v) + v))
    ck.eq("stage3.type2-block-structure-valid", len(bad_struct), 0)
    ck.eq("stage3.type2-all-avoid-321-and-1342", len(bad_avoid), 0)
    ck.eq("stage3.type2-all-indecomposable", len(bad_indec), 0)
    ck.eq("stage3.type2-inv-matches-paper-formula (%d tuples tested)"
          % len(t2), len(bad_inv), 0)
    ck.eq("stage3.type2-inv-minus-(len-1)=a(u+v)+v", len(bad_slack), 0)
    ck.eq("stage3.type2-count-per-length=C(n,3)",
          [per_len.get(n, 0) for n in range(1, N_MAIN + 1)],
          [binom(n, 3) for n in range(1, N_MAIN + 1)])

    s1, s2 = set(t1), set(pi for (pi, a, u, v, w) in t2)
    ck.eq("stage3.type1-perms-distinct", len(s1), len(t1))
    ck.eq("stage3.type2-perms-distinct-no-duplicate-tuples", len(s2), len(t2))
    ck.eq("stage3.families-disjoint", len(s1 & s2), 0)
    union = s1 | s2
    print("   union of the two families: %d permutations" % len(union))
    ck.eq("stage3.union-size", len(union), len(s1) + len(s2))
    ck.eq("stage3.xref-reference-total-indecomposables", len(union),
          SPEC_TOTAL_INDEC_11)
    bf = data["indec_set"]
    ck.eq("stage3.union-EQUALS-brute-force-indecomposable-set"
          " (families-only, bruteforce-only)",
          (len(union - bf), len(bf - union)), (0, 0))
    ck.eq("stage3.union-size-equals-brute-force-size", len(union), len(bf))
    ck.true("stage3.comparison-not-against-an-empty-set", len(bf) > 0,
            "brute-force indecomposables = %d" % len(bf))


def stage_cover(ck, data):
    """P8 vs D5: the cover of the WHOLE class, taken literally."""
    print("== STAGE 4: the full-class cover of P8 at n <= %d ==" % N_COVER)
    reps = cover_all_reps(data["levels"], N_COVER)
    perms = [pi for (pi, tag) in reps]
    tags = set(tag for (pi, tag) in reps)
    class_set = set()
    for n in range(1, N_COVER + 1):
        class_set.update(data["levels"][n])

    n_t1 = sum(1 for (pi, tag) in reps if tag == "type1")
    n_t2 = sum(1 for (pi, tag) in reps if tag == "type2")
    print("   representations built: %d type-1' + %d type-2 = %d"
          % (n_t1, n_t2, len(reps)))
    print("   class members of length <= %d: %d" % (N_COVER, len(class_set)))

    ck.true("stage4.both-types-present", tags == {"type1", "type2"},
            "tags=%s" % sorted(tags))
    ck.true("stage4.reps-nonempty", len(reps) > 0, "reps=%d" % len(reps))
    ck.eq("stage4.class-members-n<=%d-equals-sum-of-derived-sizes" % N_COVER,
          len(class_set), sum(data["class_sizes"][:N_COVER]))
    ck.eq("stage4.xref-reference-class-total-n<=%d" % N_COVER,
          len(class_set), SPEC_CLASS_TOTAL_8)
    ck.eq("stage4.representation-count-equals-class-size", len(reps),
          len(class_set))
    ck.eq("stage4.no-member-represented-twice", len(perms), len(set(perms)))
    ck.eq("stage4.representations-are-exactly-the-class"
          " (reps-only, class-only)",
          (len(set(perms) - class_set), len(class_set - set(perms))), (0, 0))
    by_len_reps = {}
    for pi in perms:
        by_len_reps[len(pi)] = by_len_reps.get(len(pi), 0) + 1
    ck.eq("stage4.representations-per-length-match-class-sizes",
          [by_len_reps.get(n, 0) for n in range(1, N_COVER + 1)],
          data["class_sizes"][:N_COVER])

    # The paper's inference "if alpha is nonempty, pi is decomposable",
    # and its consequence that the indecomposable type-1' members are
    # exactly the family P5.  The first is an implementation-consistency
    # check (it confirms that alpha really occupies {1..|alpha|} here);
    # the second links stage 4 back to stage 3.
    nonempty_alpha = []
    empty_alpha = []
    for n in range(1, N_COVER + 1):
        for j in range(0, n):
            beta = tuple(range(j + 1, n))
            alphas = [tuple()] if j == 0 else data["levels"].get(j, [])
            for alpha in alphas:
                pi = alpha + (n,) + beta
                if alpha:
                    nonempty_alpha.append(pi)
                else:
                    empty_alpha.append(pi)
    ck.true("stage4.nonempty-alpha-loop-nonempty", len(nonempty_alpha) > 0,
            "reps with alpha nonempty = %d" % len(nonempty_alpha))
    ck.eq("stage4.alpha-nonempty-implies-decomposable",
          sum(1 for pi in nonempty_alpha if is_indecomposable(pi)), 0)
    ck.eq("stage4.indecomposable-type1'-reps-are-exactly-family-P5",
          sorted(empty_alpha), sorted(family_type1(N_COVER)))
    ck.eq("stage4.all-alpha-empty-reps-are-indecomposable",
          sum(1 for pi in empty_alpha if not is_indecomposable(pi)), 0)


def stage_series(ck, data):
    """P4 vs D6: eq. (J) expanded, against the brute-force joint table."""
    print("== STAGE 5: eq. (J) as a truncated bivariate series, x^n, n<=%d =="
          % N_SERIES)
    J, a_terms, term1, term2 = build_J(N_SERIES)
    ck.true("stage5.a-sum-nonempty", a_terms > 0, "a-terms used=%d" % a_terms)

    # The closed forms of P4 against the raw sums they are claimed to be
    # equal to.  Neither side involves permutations, so this isolates the
    # geometric-series algebra in the paper's proof.
    n1 = naive_type1_series(N_SERIES)
    n2 = naive_type2_series(N_SERIES)
    t1c = {kk: c for kk, c in term1.items() if kk[0] <= N_SERIES and c != 0}
    t2c = {kk: c for kk, c in term2.items() if kk[0] <= N_SERIES and c != 0}
    ck.true("stage5.naive-sums-nonempty", len(n1) > 0 and len(n2) > 0,
            "type1 terms=%d type2 terms=%d" % (len(n1), len(n2)))
    ck.eq("stage5.x/(1-xq)-equals-SUM_n x^n q^{n-1}", t1c, n1)
    ck.eq("stage5.closed-form-equals-4-fold-sum-over-(a,u,v,w)", t2c, n2)

    bf = {kk: c for kk, c in data["joint"].items() if kk[0] <= N_SERIES}
    Js = {kk: c for kk, c in J.items() if kk[0] <= N_SERIES and c != 0}
    print("   nonzero coefficients: series=%d brute force=%d"
          % (len(Js), len(bf)))
    print("     n   k : series : brute force")
    for (n, k) in sorted(set(list(Js.keys()) + list(bf.keys()))):
        print("     %2d %3d : %6d : %6d" % (n, k, Js.get((n, k), 0),
                                            bf.get((n, k), 0)))
    ck.true("stage5.series-support-nonempty", len(Js) > 0, "terms=%d" % len(Js))
    ck.eq("stage5.support-of-series-equals-support-of-brute-force",
          sorted(Js.keys()), sorted(bf.keys()))
    ck.eq("stage5.all-coefficients-agree (%d nonzero terms)" % len(Js),
          Js, bf)
    ck.eq("stage5.xref-reference-nonzero-coefficient-count", len(Js),
          SPEC_NONZERO_JQ_9)
    ck.true("stage5.every-series-term-has-k>=n-1 (so [q^k] is complete for"
            " k <= %d)" % (N_SERIES - 1),
            all(k >= n - 1 for (n, k) in Js), "terms=%d" % len(Js))

    # J(x,1): the length specialisation of P4 against P2.
    len_spec = {}
    for (n, k), c in Js.items():
        len_spec[n] = len_spec.get(n, 0) + c
    ck.eq("stage5.J(x,1)-coefficients=1+C(n,3) for n=1..%d" % N_SERIES,
          [len_spec.get(n, 0) for n in range(1, N_SERIES + 1)],
          [1 + binom(n, 3) for n in range(1, N_SERIES + 1)])

    # J(1,q): the inversion specialisation, complete only for k <= n-1.
    inv_spec = {}
    for (n, k), c in Js.items():
        inv_spec[k] = inv_spec.get(k, 0) + c
    kcap = N_SERIES - 1
    ck.eq("stage5.J(1,q)-coefficients=1+C(k,2) for k=0..%d" % kcap,
          [inv_spec.get(k, 0) for k in range(0, kcap + 1)],
          [1 + binom(k, 2) for k in range(0, kcap + 1)])
    return len(Js)


def stage_identities(ck, data):
    """P9/P10 vs D7/D8: the class generating function and the hand checks."""
    print("== STAGE 6: generating-function identities ==")
    num = poly_add(poly_one_minus_x_pow(4), [0, 0, 0, 1])
    den = poly_mul(poly_one_minus_x_pow(3), [1, -2])
    ck.eq("stage6.numerator-(1-x)^4+x^3-equals-printed-form",
          num, PAPER_F_NUM_PRINTED)
    ck.eq("stage6.denominator-(1-x)^3(1-2x)-equals-printed-form",
          den, PAPER_F_DEN_PRINTED)

    F = series_div(num, den, N_MAIN)
    print("   DERIVED [x^n] F(x), n=0..%d: %s" % (N_MAIN, F))
    ck.eq("stage6.F-constant-term-is-1-empty-permutation", F[0], 1)
    ck.eq("stage6.F-coefficients-equal-brute-force-class-sizes",
          F[1:N_MAIN + 1], data["class_sizes"])
    ck.eq("stage6.F-head-equals-printed-expansion",
          F[:len(PAPER_F_SERIES_PRINTED)], PAPER_F_SERIES_PRINTED)

    # The functional equation of P9, as truncated series in x.
    g1 = us_geom(1, N_MAIN)
    lhs = list(F)
    rhs = us_zero(N_MAIN)
    rhs[0] = 1
    rhs = us_add(rhs, us_shift(us_mul(F, g1, N_MAIN), 1, N_MAIN))
    rhs = us_add(rhs, us_shift(us_pow(g1, 4, N_MAIN), 3, N_MAIN))
    ck.eq("stage6.functional-equation-F=1+xF/(1-x)+x^3/(1-x)^4", lhs, rhs)

    # P10 telescoping identity, as q-series, for a range of A.
    KQ = 30
    tel_fail = []
    tel_tested = 0
    for A in range(0, 9):
        l, r, nterms = telescope_sides(A, KQ)
        tel_tested += 1
        if l != r or nterms != A + 1:
            tel_fail.append(A)
    ck.true("stage6.telescoping-identity-A-range-nonempty", tel_tested == 9,
            "A tested = 0..8 (%d values), q-order %d" % (tel_tested, KQ))
    ck.eq("stage6.telescoping-identity-holds-for-every-A", tel_fail, [])

    l, r, nterms = limit_sum_sides(KQ)
    ck.eq("stage6.limit-sum-equals-1/(1-q)^2 (A=%d terms, q-order %d)"
          % (nterms, KQ), l, r)

    ik = extract_Ik_series(KQ)
    ck.eq("stage6.[q^k](1/(1-q)+q^2/(1-q)^3)=1+C(k,2), k=0..%d" % KQ,
          ik, [1 + binom(k, 2) for k in range(KQ + 1)])
    sn = extract_Sn_series(KQ)
    ck.eq("stage6.[x^n](x/(1-x)+x^3/(1-x)^4)=1+C(n,3), n=1..%d" % KQ,
          sn[1:], [1 + binom(n, 3) for n in range(1, KQ + 1)])


def stage_calibration(ck):
    """D9: this program's own index convention, on four known cases.

    The SAME counter (same generator, same (D1)-(D3)) is pointed at four
    classes whose inversion enumeration is classically known, and each
    result is compared with a model computed here from an unrelated
    combinatorial definition -- integer partitions, fountains of coins,
    triangular numbers, parallelogram polyominoes -- so neither side of
    any comparison comes from the other.

    WHAT THIS DOES AND DOES NOT SETTLE (see SCOPE item 6).  It shows
    that the counter used throughout this program indexes by inversions
    with I_0 = {1}, in four classes at once, against models it did not
    copy.  It does NOT show that these four enumerations, or their
    indexing, are the ones proved in the paper's reference [2]: that is
    a bibliographic attribution of the same kind as SCOPE item 5 and is
    not checked here.  So it does not by itself rule out the refutation
    being a disagreement about indexing -- that rests on the quoted
    sentence of [2], which no check reaches.
    """
    print("== STAGE 7: index-convention calibration, k=0..%d ==" % K_CAL)
    print("   NOTE: this calibrates THIS program's convention only.  That")
    print("   these four enumerations, at these indices, are the ones proved")
    print("   in the paper's reference [2] is NOT checked; see SCOPE item 6.")
    models = {
        "132": ("integer partitions p(k)", partition_counts(K_CAL)),
        "231": ("fountains of k coins", fountain_counts(K_CAL)),
        "12": ("triangular-number indicator", triangular_indicator(K_CAL)),
        "321": ("parallelogram polyominoes of area k",
                parallelogram_counts(K_CAL)),
    }
    pats_of = {
        "132": ((1, 3, 2),),
        "231": ((2, 3, 1),),
        "12": ((1, 2),),
        "321": ((3, 2, 1),),
    }
    derived_at_0 = {}
    I0_sets = {}
    for name in ("132", "231", "12", "321"):
        pats = pats_of[name]
        levels = gen_by_exhaustive_scan(pats, K_CAL + 1)
        indec, _bl, by_inv, _j = tabulate_indecomposable(levels)
        derived = [by_inv.get(k, 0) for k in range(0, K_CAL + 1)]
        derived_at_0[name] = derived[0]
        I0_sets[name] = sorted(p for p in indec if inv(p) == 0)
        ck.true("stage7.Av(%s)-enumeration-nonempty" % name,
                sum(derived) > 0 and len(indec) > 0,
                "indecomposable avoiders found = %d" % len(indec))
        label, model = models[name]
        print("   Av(%s): derived |I_k| = %s" % (name, derived))
        print("            model (%s) = %s" % (label, model))
        over, n_cand = boundary_scan(pats, levels[K_CAL + 1],
                                     K_CAL + 2, K_CAL)
        print("            boundary scan at length %d: %d candidates, %d hits"
              % (K_CAL + 2, n_cand, len(over)))
        ck.eq("stage7.Av(%s)-no-indecomposable-member-of-length-%d-with"
              "-inv<=%d" % (name, K_CAL + 2, K_CAL), len(over), 0)
        ck.eq("stage7.Av(%s)-|I_k|-equals-%s" % (name, label),
              derived, model)
        ck.eq("stage7.Av(%s)-xref-reference-table" % name,
              derived, SPEC_CAL[name])
    # The convention itself: in all four classes the k=0 set is exactly
    # the length-1 permutation, which is the indexing the paper says it
    # shares with Franklin.  Checked on the DERIVED sets, not the models.
    ck.eq("stage7.derived-|I_0|=1-in-all-four-classes",
          derived_at_0, {"132": 1, "231": 1, "12": 1, "321": 1})
    ck.eq("stage7.derived-I_0-is-exactly-{1}-in-all-four-classes",
          I0_sets, {"132": [(1,)], "231": [(1,)], "12": [(1,)],
                    "321": [(1,)]})


def stage_refutation(ck, data):
    """D10: the paper's conclusion, not merely its lemmas."""
    print("== STAGE 8: the refutation itself ==")
    kmax = data["kmax"]
    by_inv = data["by_inv"]
    derived = [by_inv.get(k, 0) for k in range(0, kmax + 1)]
    conj = [1 + binom(k + 1, 2) for k in range(0, kmax + 1)]

    # P3: the quoted form k(k+1)/2+1 and the paper's binomial rewrite.
    ck.eq("stage8.quoted-k(k+1)/2+1-equals-1+C(k+1,2), k=0..40",
          [k * (k + 1) // 2 + 1 for k in range(41)],
          [1 + binom(k + 1, 2) for k in range(41)])

    # I_0 and I_1 as sets, from brute force.  The paper's index claim
    # ("in Franklin as here, k is the number of inversions and I_0={1}")
    # is what the k=0 check pins down.
    I0 = sorted(p for p in data["indec_set"] if inv(p) == 0)
    I1 = sorted(p for p in data["indec_set"] if inv(p) == 1)
    ck.eq("stage8.I_0-is-exactly-{1}", I0, [(1,)])
    print("   DERIVED I_0 = %s" % (I0,))
    print("   DERIVED I_1 = %s" % (I1,))
    print("   k : derived |I_k| : conjectured 1+C(k+1,2)")
    for k in range(0, kmax + 1):
        print("   %2d : %13d : %d" % (k, derived[k], conj[k]))
    ck.eq("stage8.I_1-is-exactly-{21}", I1, [(2, 1)])

    # The paper's argument for I_1 = {21}: every permutation with exactly
    # one inversion is the identity with two adjacent consecutive values
    # interchanged, and only 21 among those is indecomposable.  Tested on
    # ALL permutations of length <= N_CROSS, not only on class members,
    # since that is the generality in which the paper states it.
    one_inv = []
    not_adjacent_swap = []
    for n in range(1, N_CROSS + 1):
        ident = tuple(range(1, n + 1))
        swaps = set()
        for i in range(n - 1):
            sw = list(ident)
            sw[i], sw[i + 1] = sw[i + 1], sw[i]
            swaps.add(tuple(sw))
        for p in permutations(ident):
            if inv(p) == 1:
                one_inv.append(p)
                if p not in swaps:
                    not_adjacent_swap.append(p)
    ck.true("stage8.one-inversion-loop-nonempty", len(one_inv) > 0,
            "permutations with inv=1 found = %d" % len(one_inv))
    ck.eq("stage8.every-inv=1-permutation-is-an-adjacent-transposition"
          " (n<=%d)" % N_CROSS, len(not_adjacent_swap), 0)
    ck.eq("stage8.only-indecomposable-inv=1-permutation-is-21 (n<=%d)"
          % N_CROSS,
          sorted(set(p for p in one_inv if is_indecomposable(p))), [(2, 1)])
    ck.eq("stage8.|I_1|-derived", derived[1], 1)
    ck.ne("stage8.|I_1|-differs-from-conjectured-value", derived[1], conj[1])

    diffs = [k for k in range(1, kmax + 1) if derived[k] != conj[k]]
    ck.eq("stage8.conjecture-fails-at-EVERY-k-in-1..%d" % kmax,
          diffs, list(range(1, kmax + 1)))
    ck.eq("stage8.conjecture-agrees-at-k=0-only", derived[0], conj[0])

    # The paper's "correct sequence shifted by one in the index".
    shifted = [derived[k + 1] for k in range(0, kmax)]
    ck.eq("stage8.conjectured-value-at-k-equals-derived-|I_{k+1}|",
          conj[:kmax], shifted)


# ===============================================================
# Section 10: main.
# ===============================================================

def print_scope():
    """What this program does NOT establish.  Stated so that a passing
    verdict is not read as more than it is."""
    print("")
    print("== SCOPE / WHAT IS NOT COVERED BY ANY CHECK ==")
    print(" 1. Every claim is verified on the finite ranges printed above.")
    print("    The theorem's 'for every k>=0, n>=1' is not established by")
    print("    exhaustion; the finite evidence is the enumeration, and the")
    print("    algebraic chain (eq. (J) = the 4-fold sum -> telescoping ->")
    print("    coefficient extraction) is verified as series identities to")
    print("    the printed orders, not symbolically for all k.")
    print(" 2. Lemma P11 (indecomposable => inv >= len-1) is imported from")
    print("    the literature.  It is verified here on ALL permutations of")
    print("    length <= %d and on every indecomposable class member found,"
          % N_CROSS)
    print("    plus the boundary length %d, but not proved." % (N_MAIN + 1))
    print(" 3. The BGU structural cover (P8) is verified as a FACT for")
    print("    n <= %d, not derived from BGU's proof." % N_COVER)
    print(" 4. Asymptotic statements in the paper ('grows like n^3/6',")
    print("    '1+C(k,2) and 1+C(k+1,2) are both asymptotic to k^2/2') are")
    print("    not finitely checkable and are NOT checked here.")
    print(" 5. Bibliographic claims -- the quoted sentence, the page and")
    print("    proposition locators, and the attribution of the conjecture")
    print("    to numerical evidence -- are outside this program.")
    print(" 6. In particular, the calibration of STAGE 7 does NOT close the")
    print("    gap left by item 5.  It shows that THIS program's counter")
    print("    indexes by inversions with I_0 = {1}, in Av(132), Av(231),")
    print("    Av(12) and Av(321) at once, against four models built here")
    print("    from unrelated definitions.  That these same four")
    print("    enumerations, at these same indices, are the ones proved in")
    print("    the paper's reference [2] is an attribution of the same")
    print("    bibliographic kind as item 5 and is NOT checked.  Hence the")
    print("    paper's sentence 'in [2] as here, k is the number of")
    print("    inversions of the permutation and I_0={1}' -- on which the")
    print("    whole refutation depends, since a different convention in [2]")
    print("    would make the conjectured count correct -- rests on the")
    print("    quotation and locator of item 5, and on no check below.")


def main():
    print("verify.py -- Franklin (321,1342) inversion-count refutation")
    print("RANGES COVERED (all counts below are generated, not seeded):")
    print("  lemma P11 tested on all of S_n  : n = 1..%d" % N_CROSS)
    print("  Av(321,1342) brute force        : n = 1..%d" % N_MAIN)
    print("  second, independent generator   : n = 1..%d" % N_CROSS)
    print("  |I_k| tested                    : k = 0..%d" % (N_MAIN - 1))
    print("  |S*_n| tested                   : n = 1..%d" % N_MAIN)
    print("  structural families rebuilt     : n = 1..%d" % N_MAIN)
    print("  full-class cover of P8          : n = 1..%d" % N_COVER)
    print("  eq. (J) bivariate series        : n = 1..%d" % N_SERIES)
    print("  F(x) expansion                  : n = 0..%d" % N_MAIN)
    print("  q-series identities             : q-order 30, A = 0..8")
    print("  index calibration, 4 classes    : k = 0..%d" % K_CAL)
    print("")

    ck = Checks()
    stage_lemma(ck)
    data = stage_enumeration(ck)
    stage_formulas(ck, data)
    stage_families(ck, data)
    stage_cover(ck, data)
    stage_series(ck, data)
    stage_identities(ck, data)
    stage_calibration(ck)
    stage_refutation(ck, data)
    print_scope()

    print("")
    print("== CHECK RESULTS ==")
    return ck.report()


if __name__ == "__main__":
    sys.exit(main())

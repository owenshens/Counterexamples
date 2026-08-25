#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- referee's verification program for

    "A Counterexample to the Etzion-Vardy-Yaakobi Conference-Matrix
     Conjecture"

Standard library only.  Exact integer / finite-field arithmetic throughout;
no floating point is used in any decision.  Single process.

-------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER  (inputs -- transcribed, NOT checked)
-------------------------------------------------------------------------
  P_ROWS      the 30 printed rows of W over the alphabet {0,+,-}, exactly
              as they appear in the verbatim display in the paper
  P_X         x = (25,5,26,12,0,6,19,16,0,22,5,1,26,23,14,0^15) in F_29^30
  P_XW        the printed product xW = (1,0^14,15,6,0,8,27,0,18,14,0,15,
              4,11,20,17,19)
  P_SUPPORT   the printed support {1,16,17,19,20,22,23,25,26,27,28,29,30}
  P_WEIGHT    13            ("Hamming weight 13", "d(C) <= 13")
  P_DET       -8584878      (leading 15x15 principal minor of W)
  P_DET_MOD   21            (-8584878 = 21 mod 29)
  P_P         29            the prime p
  P_N, P_K    30, 15        the claimed code parameters [30,15]_29
  P_MDS_D     16            "an MDS [30,15]_29 code has d = 30-15+1 = 16"

-------------------------------------------------------------------------
WHAT IS DERIVED HERE  (the checks; nothing below is copied from the paper)
-------------------------------------------------------------------------
  * P_ROWS is decoded to a 30x30 integer matrix and its shape, alphabet,
    zero diagonal, +-1 off-diagonal entries and symmetry are recomputed;
  * all 900 integer inner products of rows of W, i.e. W W^T = 29 I_30
    (30 diagonal entries = 29, all 435 unordered off-diagonal pairs = 0),
    so that W is verified to BE a conference matrix of order 30 -- the
    hypothesis of Conjecture 5 -- rather than assumed to be one; and,
    independently of those inner products, |det W| = 29^15, the exact
    30x30 fraction-free integer determinant forced by W W^T = 29 I_30;
    (the 30 diagonal inner products are NOT independent evidence: with a
    zero diagonal and off-diagonal entries in {+1,-1} each row must have
    self-product n-1, so that line cannot fail unless one of the two
    alphabet checks does -- the off-diagonal orthogonality is the content)
  * the exact leading 15x15 principal minor of W, by fraction-free
    Bareiss elimination and, independently, by Gaussian elimination over
    the rationals; and its residue mod 29;
  * rank_{F_29}(W) by Gauss-Jordan elimination over F_29, its reduced row
    echelon form, and the pivot columns;
  * self-duality: W W^T = 0 in F_29 (C is contained in C-perp), the
    F_29 null space of W (= C-perp) with its dimension AND its computed
    F_29 rank (so dim C-perp = 15 is measured, not inherited from the way
    the basis was constructed; with the nonzero minor this pins
    rank_{F_29}(W) = 15 without trusting the pivot count), and the fact
    that every null-space basis vector lies in the row space of W (C-perp
    is contained in C), giving C = C-perp by computation, not by citation.
    Both containment lines require the basis to have the full size n-k:
    an empty or short basis satisfies "every basis vector lies in C"
    vacuously, and a full-rank corrupted W would otherwise PASS them;
  * xW recomputed in F_29 from the decoded W and compared entry by entry
    with P_XW; its support and Hamming weight;
  * the same weight-13 codeword re-derived WITHOUT using x at all, as the
    first row of the systematic generator matrix [I_15 | A] obtained from
    the reduced row echelon form of W;
  * the Singleton bound d = n-k+1 for the DERIVED n and k, the conjectured
    parameters [p+1,(p+1)/2,(p+3)/2] evaluated at the derived p, and the
    strict inequality (weight found) < (MDS distance) -- the load-bearing
    refutation step;
  * a partial minimum-distance census: the exact minimum Hamming weight
    over ALL nonzero codewords of C whose information vector (in the
    systematic coordinates) has weight at most the census depth (default
    5, i.e. 5.25 x 10^10 codewords), enumerated exhaustively up to scalar
    multiples, with its witness re-verified to lie in the F_29 row span of
    W.  This confirms d(C) <= 13 independently of the paper's x, and in
    fact yields d(C) <= 12; the paper's bound is correct but not tight.
    What is NOT re-run is printed as a NOTE.

The closing notes also disclose the make-up of the headline count, with the
figures derived from the check log rather than written out by hand: which
lines are machinery self-tests, which are logical restatements of an earlier
line and so cannot fail on their own, and how many lines are left carrying
force independent of every other line.  The switching orbit under which every
check is insensitive is likewise computed from x and xW, not asserted.

Usage:  python3 verify.py              (default, about 200 s)
        python3 verify.py --fast       (census depth 4, about 3 s)
        python3 verify.py --quick      (census depth 2, instant)
Everything except the census finishes in well under a second.
Exit status is 0 if and only if every check passed.
"""

import sys
import time

# =====================================================================
# CHECK HARNESS
# =====================================================================

_RESULTS = []


def check(name, condition, detail=""):
    """Record and print one check.  Every call must be able to print FAIL."""
    ok = bool(condition)
    _RESULTS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return ok


def note(text):
    print("NOTE " + text)


def verdict():
    n = len(_RESULTS)
    bad = [nm for nm, ok in _RESULTS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


# =====================================================================
# BLOCK 1.  VERBATIM PAPER INPUTS.  Corrupt any character here and the
# checks below must report FAIL.
# =====================================================================

# The 30 printed rows of W, alphabet 0/+/- meaning 0/+1/-1.
P_ROWS = [
    "0+++++++++++++++++++++++++++++",
    "+0--------------++++++++++++++",
    "+-0------+++++++-------+++++++",
    "+--0----+-++++++-++++++------+",
    "+---0+++----+++++---+++---+++-",
    "+---+0++-+-+--++--++-+++++----",
    "+---++0++-++-+--+-++-----+-+++",
    "+---+++0++--++---+-++--++-+--+",
    "+--+--++0+++---++++-+---+-+-+-",
    "+-+--+-++0+-+-+-+++---++-++---",
    "+-++--+-++0--++-+--+++-+-+--+-",
    "+-++-++-+--0+--++-+--+-++--+-+",
    "+-+++--+-+-+0-+-++-+-+----++-+",
    "+-+++-++--+--0-+-+--+-++-+-+-+",
    "+-++++---++-+-0---+++-+-+--++-",
    "+-++++--+--+-+-0-+---++-+++-+-",
    "++--+-+-+++++---0+---+++---++-",
    "++-+---+++--++-++0-+--+-++-+--",
    "++-+-++-++-+--+---0-+-+--+++-+",
    "++-+-+++--+-+-+--+-0-+--++--++",
    "++-++--++-+--++---+-0+-++-++--",
    "++-+++----+++--++--++0-+-++---",
    "++-+++---+---++++++---0+----++",
    "+++--+-+-+++-+--+---+++0+----+",
    "+++--+-++--+--++-+-++--+0--++-",
    "+++--++--++--+-+-+++-+---0++--",
    "+++-+--+++--+--+--+-++---+0-++",
    "+++-+-+----++++-+++-+---++-0--",
    "+++-+-+-+-+---+++--+--+-+-+-0+",
    "++++--++---+++----++--++--+-+0",
]

P_P = 29                      # the prime p
P_X = [25, 5, 26, 12, 0, 6, 19, 16, 0, 22, 5, 1, 26, 23, 14] + [0] * 15
P_XW = [1] + [0] * 14 + [15, 6, 0, 8, 27, 0, 18, 14, 0, 15, 4, 11, 20, 17, 19]
P_SUPPORT = [1, 16, 17, 19, 20, 22, 23, 25, 26, 27, 28, 29, 30]
P_WEIGHT = 13                 # "Hamming weight 13"; "d(C) <= 13"
P_DET = -8584878              # leading 15x15 principal minor
P_DET_MOD = 21                # its residue mod 29
P_N = 30                      # claimed length
P_K = 15                      # claimed dimension
P_MDS_D = 16                  # claimed MDS distance 30-15+1


# =====================================================================
# BLOCK 2.  EXACT ARITHMETIC HELPERS (no floating point anywhere)
# =====================================================================

SYMBOL = {"0": 0, "+": 1, "-": -1}


def decode(rows):
    """Decode the 0/+/- display into a matrix of integers, or None."""
    out = []
    for r in rows:
        if any(c not in SYMBOL for c in r):
            return None
        out.append([SYMBOL[c] for c in r])
    return out


def is_prime(m):
    if m < 2:
        return False
    d = 2
    while d * d <= m:
        if m % d == 0:
            return False
        d += 1
    return True


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def det_bareiss(A):
    """Exact determinant by fraction-free Bareiss elimination."""
    A = [row[:] for row in A]
    m = len(A)
    sign = 1
    prev = 1
    for k in range(m - 1):
        if A[k][k] == 0:
            for r in range(k + 1, m):
                if A[r][k] != 0:
                    A[k], A[r] = A[r], A[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, m):
            for j in range(k + 1, m):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[m - 1][m - 1]


def det_rational(A):
    """Exact determinant, second independent route: Gaussian elimination
    over the rationals.  Returns a Fraction (integral for integer input)."""
    from fractions import Fraction
    A = [[Fraction(x) for x in row] for row in A]
    m = len(A)
    det = Fraction(1)
    for k in range(m):
        piv = None
        for r in range(k, m):
            if A[r][k] != 0:
                piv = r
                break
        if piv is None:
            return Fraction(0)
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            det = -det
        det *= A[k][k]
        for r in range(k + 1, m):
            f = A[r][k] / A[k][k]
            if f:
                for c in range(k, m):
                    A[r][c] -= f * A[k][c]
    return det


def rref_mod(A, p):
    """Gauss-Jordan over F_p.  Returns (reduced rows, pivot columns)."""
    A = [[x % p for x in row] for row in A]
    m = len(A)
    ncol = len(A[0]) if m else 0
    piv = []
    r = 0
    for c in range(ncol):
        pr = None
        for i in range(r, m):
            if A[i][c]:
                pr = i
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(ncol)]
        piv.append(c)
        r += 1
        if r == m:
            break
    return A, piv


def null_space_mod(A, p):
    """Basis of the right null space {y : A y^T = 0} over F_p."""
    R, piv = rref_mod(A, p)
    ncol = len(A[0])
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for f in free:
        y = [0] * ncol
        y[f] = 1
        for i, c in enumerate(piv):
            y[c] = (-R[i][f]) % p
        basis.append(y)
    return basis


def in_row_space(v, R, piv, p):
    """True iff v is an F_p combination of the reduced rows R (pivots piv)."""
    w = [x % p for x in v]
    for i, c in enumerate(piv):
        f = w[c]
        if f:
            w = [(w[j] - f * R[i][j]) % p for j in range(len(w))]
    return all(x == 0 for x in w)


def fit(v, n):
    """Pad/truncate a paper vector to length n, so that a corrupted input
    of the wrong length still produces FAIL lines instead of a traceback."""
    return (list(v) + [0] * n)[:n]


def weight(v):
    return sum(1 for x in v if x % P_P != 0)


def support(v):
    return [j + 1 for j in range(len(v)) if v[j] % P_P != 0]


def census_min_weight(A, k, p, depth):
    """EXACT minimum Hamming weight of a codeword of the systematic code
    [I_k | A] whose information vector has weight in [1, depth].

    Every such codeword is enumerated (one representative per scalar class,
    leading coefficient fixed to 1; scaling does not change the weight).
    Returns (best_weight, best_info_pairs, n_representatives_examined).
    """
    m2 = len(A[0])
    inv = [0] + [pow(a, p - 2, p) for a in range(1, p)]
    mul = [[tuple((c * a) % p for a in row) for c in range(p)] for row in A]
    best = [10 ** 9, None]
    count = [0]

    def finalize(v, i, mlen, pre):
        """Weights of v + c*A[i] for every c in 1..p-1, without building
        p-1 vectors: v_j + c a_j = 0 iff a_j != 0 and c = -v_j/a_j."""
        a = A[i]
        cnt = [0] * p
        z0 = 0
        for j in range(m2):
            aj = a[j]
            if aj == 0:
                if v[j] == 0:
                    z0 += 1
            else:
                cnt[((-v[j]) * inv[aj]) % p] += 1
        base = mlen + 1 + m2 - z0
        count[0] += p - 1
        for c in range(1, p):
            w = base - cnt[c]
            if w < best[0]:
                best[0] = w
                best[1] = pre + [(i, c)]

    def go(start, v, mlen, budget, pre):
        for i in range(start, k):
            finalize(v, i, mlen, pre)
            if budget > 1:
                mu = mul[i]
                for c in range(1, p):
                    mc = mu[c]
                    nv = tuple((v[j] + mc[j]) % p for j in range(m2))
                    go(i + 1, nv, mlen + 1, budget - 1, pre + [(i, c)])

    for i1 in range(k):
        w = 1 + sum(1 for a in A[i1] if a % p)
        count[0] += 1
        if w < best[0]:
            best[0] = w
            best[1] = [(i1, 1)]
        if depth >= 2:
            go(i1 + 1, tuple(A[i1]), 1, depth - 1, [(i1, 1)])
    return best[0], best[1], count[0]


def vec_mat_mod(x, M, p):
    """x M over F_p, x a row vector, M a list of rows."""
    ncol = len(M[0])
    out = [0] * ncol
    for i, xi in enumerate(x):
        if xi % p:
            row = M[i]
            for j in range(ncol):
                out[j] = (out[j] + xi * row[j]) % p
    return out


# =====================================================================
# BLOCK 3.  THE CHECKS
# =====================================================================

def stage_well_formed():
    """Decode the display, print it back, and verify W IS a conference
    matrix of order 30 -- the hypothesis of Conjecture 5."""
    p = P_P
    W = decode(P_ROWS)
    check("display_decodes_over_alphabet_0_plus_minus", W is not None,
          "%d rows" % len(P_ROWS))
    if W is None:
        sys.exit(verdict())
    n = len(W)
    square = (n > 0 and all(len(r) == n for r in W))
    check("display_is_square_with_rows_of_equal_length", square,
          "rows=%d, row lengths=%s" % (n, sorted(set(len(r) for r in W))))
    if not square:
        note("aborting: a non-square array admits none of the checks below")
        sys.exit(verdict())
    check("order_of_the_matrix_is_30", n == P_N,
          "order %d, paper %d" % (n, P_N))

    note("decoded matrix, re-encoded (compare with the paper's display):")
    inv_sym = {0: "0", 1: "+", -1: "-"}
    for i, row in enumerate(W):
        note("  W[%2d] %s" % (i + 1, "".join(inv_sym[e] for e in row)))
    npos = sum(1 for row in W for e in row if e == 1)
    nneg = sum(1 for row in W for e in row if e == -1)
    nzero = sum(1 for row in W for e in row if e == 0)
    note("entry census: %d(+1) %d(-1) %d(0), total %d"
         % (npos, nneg, nzero, npos + nneg + nzero))

    diag_ok = all(W[i][i] == 0 for i in range(n))
    off_ok = all(W[i][j] in (1, -1) for i in range(n)
                 for j in range(n) if i != j)
    check("zero_diagonal", diag_ok,
          "%d diagonal entries" % n)
    check("off_diagonal_entries_all_plus_or_minus_one", off_ok,
          "%d(+1) %d(-1) %d(0)" % (npos, nneg, nzero))
    check("matrix_is_symmetric", all(W[i][j] == W[j][i]
                                    for i in range(n) for j in range(i)),
          "%d unordered pairs" % (n * (n - 1) // 2))
    check("p_is_prime_and_order_equals_p_plus_1",
          is_prime(p) and n == p + 1, "p=%d, order=%d" % (p, n))
    return W


def stage_conference(W):
    """W W^T = 29 I_30 by exact integer arithmetic: 30 diagonal entries
    and all 435 unordered off-diagonal pairs."""
    n = len(W)
    p = P_P
    dbad = [(i, dot(W[i], W[i])) for i in range(n)
            if dot(W[i], W[i]) != n - 1]
    check("gram_diagonal_all_equal_n_minus_1", not dbad,
          "%d rows tested, target %d, offenders %s"
          % (n, n - 1, dbad[:3]))
    obad = []
    npairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            npairs += 1
            s = dot(W[i], W[j])
            if s != 0:
                obad.append((i, j, s))
    check("gram_off_diagonal_all_zero", not obad,
          "%d unordered pairs tested, offenders %s" % (npairs, obad[:3]))
    check("number_of_off_diagonal_pairs_is_435", npairs == n * (n - 1) // 2
          and npairs == 435, "%d" % npairs)
    check("conference_constant_n_minus_1_equals_p", n - 1 == p,
          "n-1=%d, p=%d" % (n - 1, p))
    # A bare restatement of the two Gram checks above could never fail on its
    # own, so this line also verifies an INDEPENDENT
    # consequence of W W^T = (n-1) I_n by completely different arithmetic:
    # det(W)^2 = (n-1)^n, i.e. |det W| = 29^15, from an exact 30x30
    # fraction-free integer determinant that touches none of the 900 inner
    # products.  A bug in the Gram enumeration therefore also reports FAIL.
    detW = det_bareiss(W)
    det_target = (n - 1) ** (n // 2)
    det_ok = (n % 2 == 0 and abs(detW) == det_target)
    check("W_is_a_conference_matrix_of_order_30",
          not dbad and not obad and n == 30 and det_ok,
          "W W^T = %d I_%d, det W = %d, |det W| = %d^%d = %d ? %s"
          % (n - 1, n, detW, n - 1, n // 2, det_target, det_ok))
    # C is contained in C-perp: the same Gram matrix over F_29.
    zbad = [(i, j) for i in range(n) for j in range(n)
            if dot(W[i], W[j]) % p != 0]
    check("gram_matrix_vanishes_mod_p_so_C_is_self_orthogonal", not zbad,
          "%d ordered pairs tested, offenders %s" % (n * n, zbad[:3]))


def stage_minor(W):
    """The exact leading 15x15 principal minor, two independent ways."""
    m = min(P_K, len(W))
    sub = [row[:m] for row in W[:m]]
    d1 = det_bareiss(sub)
    d2 = det_rational(sub)
    check("leading_15x15_minor_equals_minus_8584878_bareiss", d1 == P_DET,
          "computed %d, paper %d" % (d1, P_DET))
    check("leading_15x15_minor_equals_minus_8584878_rational",
          d2 == P_DET, "computed %s, paper %d" % (d2, P_DET))
    check("minor_residue_mod_p_equals_21", d1 % P_P == P_DET_MOD,
          "%d mod %d = %d, paper %d" % (d1, P_P, d1 % P_P, P_DET_MOD))
    check("minor_is_nonzero_mod_p_so_rank_at_least_15", d1 % P_P != 0,
          "residue %d" % (d1 % P_P))
    return d1


def stage_code(W):
    """rank, systematic form, and C = C-perp, all by computation."""
    n = len(W)
    p = P_P
    R, piv = rref_mod(W, p)
    k = len(piv)
    check("rank_over_F_29_equals_15", k == P_K,
          "computed rank %d, paper %d" % (k, P_K))
    check("pivot_columns_are_the_first_15", piv == list(range(P_K)),
          "pivots %s" % (piv,))
    check("dimension_is_half_the_length", 2 * k == n,
          "2*%d vs n=%d" % (k, n))
    G = [R[i] for i in range(k)]
    ns = null_space_mod(W, p)
    annih = all(dot(W[i], y) % p == 0 for y in ns for i in range(n))
    # The exhibited null-space basis is verified INDEPENDENT over F_29 rather
    # than left as an unstated invariant of null_space_mod (a machinery
    # cross-check: it cannot fail unless that routine is wrong).  Its point is
    # that "dim C-perp = 15" then rests on a computed rank, so with the nonzero
    # minor (rank >= 15) rank_{F_29}(W) = 15 does not depend on the pivot
    # count alone: 15 independent annihilators force rank <= 15.
    ns_indep = bool(ns) and len(rref_mod(ns, p)[1]) == len(ns)
    check("C_perp_has_dimension_15_and_annihilates_W",
          len(ns) == P_N - P_K and annih and ns_indep,
          "dim C-perp = %d, paper n-k = %d, basis rank %d, "
          "W y = 0 for all basis y: %s"
          % (len(ns), P_N - P_K, len(rref_mod(ns, p)[1]) if ns else 0, annih))
    # NOTE ON VACUITY.  "every basis vector lies in C" is satisfied by an EMPTY
    # basis, so a corrupted W of full F_29 rank (C-perp = 0) would PASS both
    # containment lines while C != C-perp; and dim C-perp = n - k is a
    # rank-nullity IDENTITY of this code, so binding to it would be no test at
    # all.  The condition that actually makes containment equivalent to
    # equality is dim C-perp = dim C, i.e. len(ns) = k, and that is what is
    # required here: it fails for every rank other than n/2.
    ns_dim_ok = (len(ns) == k and len(ns) > 0)
    inside = ns_dim_ok and all(in_row_space(y, R, piv, p) for y in ns)
    check("every_C_perp_basis_vector_lies_in_C", inside,
          "%d basis vectors (dim C = %d must match, else containment is not "
          "equality) reduced against the row space" % (len(ns), k))
    # Same containment again, but solved against the RAW rows of W, so it
    # does not rely on the row reduction having preserved the row space.
    raw = ns_dim_ok and all(solve_left(W, y, p) is not None for y in ns)
    check("C_perp_basis_solves_over_the_raw_rows_of_W", raw,
          "%d systems y = z W solved exactly over F_29 (dim C = %d required)"
          % (len(ns), k))
    check("code_is_euclidean_self_dual", inside and len(ns) == k,
          "dim C = dim C-perp = %d and C-perp is contained in C" % k)
    gbad = [(i, j) for i in range(k) for j in range(k)
            if dot(G[i], G[j]) % p != 0]
    check("systematic_generator_gram_vanishes_mod_p", not gbad,
          "%d ordered pairs, offenders %s" % (k * k, gbad[:3]))
    return R, piv, G, k


def stage_codeword(W, G, k):
    """Recompute xW in F_29, and re-derive the same word without x."""
    n = len(W)
    p = P_P
    nz = [i for i, e in enumerate(P_X) if e]
    xok = (len(P_X) == n and all(0 <= e < p for e in P_X)
           and any(e % p for e in P_X)
           and all(e == 0 for e in P_X[k:]))
    check("x_is_a_well_formed_nonzero_vector_of_F_29_30", xok,
          "len=%d, entries in [0,%d), %d trailing zeros"
          % (len(P_X), p, (len(P_X) - 1 - nz[-1]) if nz else len(P_X)))
    xt = fit(P_X, n)
    xW = vec_mat_mod(xt, W, p)
    xWp = fit(P_XW, n)
    mism = [(j + 1, xW[j], xWp[j]) for j in range(n) if xW[j] != xWp[j]]
    check("xW_recomputed_matches_the_printed_product", not mism,
          "computed %s, mismatches %s" % (xW, mism[:3]))
    w = weight(xW)
    check("xW_has_hamming_weight_13", w == P_WEIGHT,
          "computed weight %d, paper %d" % (w, P_WEIGHT))
    check("xW_support_matches_the_printed_support",
          support(xW) == P_SUPPORT, "computed %s" % (support(xW),))
    check("xW_is_a_nonzero_codeword", w > 0, "weight %d" % w)
    # Independent re-derivation: xW must be the first row of [I_15 | A],
    # which is computed from W alone and never touches x.
    g0 = G[0] if G else []
    check("codeword_re_derived_from_rref_without_using_x", g0 == xW,
          "rref row 1 = %s" % (g0,))
    check("xW_information_part_is_the_first_unit_vector",
          xW[:k] == [1] + [0] * (k - 1), "%s" % (xW[:k],))
    return xW, w


def stage_refutation(n, k, w):
    """The load-bearing step: the conclusion of Conjecture 5 fails."""
    p = P_P
    d_singleton = n - k + 1
    check("singleton_bound_for_derived_n_and_k_is_16",
          d_singleton == P_MDS_D,
          "n-k+1 = %d-%d+1 = %d, paper %d" % (n, k, d_singleton, P_MDS_D))
    conj = (p + 1, (p + 1) // 2, (p + 3) // 2)
    check("conjectured_parameters_p1_p1half_p3half_are_30_15_16",
          conj == (n, k, d_singleton) and (p + 1) % 2 == 0,
          "conjecture predicts %s, derived (n,k,MDS d) = %s"
          % (conj, (n, k, d_singleton)))
    check("derived_dimension_agrees_with_the_conjecture", k == conj[1],
          "k=%d" % k)
    check("minimum_distance_upper_bound_beats_the_conjectured_distance",
          w < conj[2], "d(C) <= %d < %d" % (w, conj[2]))
    check("code_is_not_MDS", w < d_singleton,
          "d(C) <= %d, MDS would need %d" % (w, d_singleton))
    return d_singleton


def solve_left(rows, target, p):
    """Coefficients y with sum_i y_i rows[i] = target over F_p, or None.
    All len(target) equations are used, so consistency is verified."""
    nr = len(rows)
    nc = len(target)
    aug = [[rows[i][j] % p for i in range(nr)] + [target[j] % p]
           for j in range(nc)]
    R, piv = rref_mod(aug, p)
    if nr in piv:                      # pivot in the augmented column
        return None
    y = [0] * nr
    for i, c in enumerate(piv):
        y[c] = R[i][nr]
    for j in range(nc):
        if sum(y[i] * rows[i][j] for i in range(nr)) % p != target[j] % p:
            return None
    return y


def binom(a, b):
    r = 1
    for i in range(b):
        r = r * (a - i) // (i + 1)
    return r


# Minimum weight over information weight <= 2, DERIVED BY THIS PROGRAM
# (it is not a value from the paper; the paper only claims d(C) <= 13).
D_CENSUS_MIN = 12


def stage_census(W, R, piv, G, k, depth):
    """Partial minimum-distance census: the EXACT minimum weight over all
    codewords whose systematic information vector has weight <= depth."""
    n = len(W)
    p = P_P
    pre = (piv == list(range(k)) and k == P_K)
    check("census_prerequisite_generator_is_systematic", pre,
          "pivots %s" % (piv,))
    if not pre:
        note("census SKIPPED: W has no [I_k | A] form, so the information "
             "coordinates are not the first k")
        return None
    A = [g[k:] for g in G]
    t0 = time.time()
    bw, pairs, cnt = census_min_weight(A, k, p, depth)
    secs = time.time() - t0
    total = sum(binom(k, m) * (p - 1) ** (m - 1) for m in range(1, depth + 1))
    check("census_enumeration_is_complete_for_its_range", cnt == total,
          "%d scalar classes examined, %d expected" % (cnt, total))
    u = [0] * k
    for i, c in pairs:
        u[i] = c
    cw = vec_mat_mod(u, G, p)
    check("census_witness_weight_equals_the_census_minimum",
          weight(cw) == bw, "witness weight %d, census min %d"
          % (weight(cw), bw))
    # ON VACUITY AND TAUTOLOGY.  cw was built as u*G and the rows of G are
    # elementary row operations applied to W, so "cw is an F_29 combination of
    # the rows of W" holds for EVERY input on which the census runs: the raw
    # solve alone can never report FAIL.  The membership is therefore also
    # decided by the independent DUAL criterion -- cw lies in C = (C-perp)-perp
    # iff cw is orthogonal to every basis vector of the computed null space --
    # and the witness is required to be a NONZERO word with a full-size dual
    # basis, so a degenerate census cannot pass this line by default.
    y = solve_left([row[:] for row in W], cw, p)
    ns2 = null_space_mod(W, p)
    dual_ok = (len(ns2) == k and len(ns2) > 0
               and all(dot(cw, yv) % p == 0 for yv in ns2))
    cw_nonzero = any(v % p for v in cw)
    check("census_witness_lies_in_the_F_29_row_span_of_W",
          y is not None and cw_nonzero and dual_ok,
          "witness support %s; raw solve %s, nonzero %s, orthogonal to all "
          "%d C-perp basis vectors %s"
          % (support(cw), y is not None, cw_nonzero, len(ns2), dual_ok))
    check("census_minimum_is_at_most_the_papers_bound_13", bw <= P_WEIGHT,
          "census min %d, paper's claimed bound %d" % (bw, P_WEIGHT))
    if depth >= 2:
        check("census_minimum_equals_12_derived_here", bw == D_CENSUS_MIN,
              "census min %d, value derived by this program %d"
              % (bw, D_CENSUS_MIN))
    note("census depth %d: min weight %d, %d scalar classes "
         "(= %d codewords), %.1f s"
         % (depth, bw, cnt, cnt * (p - 1), secs))
    note("census witness: information support %s, codeword %s"
         % ([(i + 1, c) for i, c in pairs], cw))
    return bw


# =====================================================================
# BLOCK 4.  DRIVER
# =====================================================================

CENSUS_DEPTH_DEFAULT = 5      # ~200 s here; depth 6 would need hours


# ---------------------------------------------------------------------
# CHECK INVENTORY.  A headline count is only honest if the lines that
# cannot fail on their own are named, so both groups are listed here and
# printed, with counts derived from _RESULTS, at the end of every run.
#
#   MACHINERY_SELF_TESTS   test this program rather than the paper: a count
#                          against a closed form, or a second route that must
#                          agree with a first.
#   DEPENDENT_RESTATEMENTS (name, the earlier check that already implies it):
#                          each is printed because it names a step of the
#                          argument, but it adds no independent failure mode.
# Nothing here changes what is checked; it changes only what is disclosed.
# ---------------------------------------------------------------------

MACHINERY_SELF_TESTS = [
    "number_of_off_diagonal_pairs_is_435",
    "census_enumeration_is_complete_for_its_range",
    "census_witness_weight_equals_the_census_minimum",
    "census_witness_lies_in_the_F_29_row_span_of_W",
]

DEPENDENT_RESTATEMENTS = [
    ("conference_constant_n_minus_1_equals_p",
     "p_is_prime_and_order_equals_p_plus_1, which already forces n = p+1"),
    ("minor_is_nonzero_mod_p_so_rank_at_least_15",
     "minor_residue_mod_p_equals_21, and 21 is not 0"),
    ("dimension_is_half_the_length",
     "rank_over_F_29_equals_15 with order_of_the_matrix_is_30"),
    ("code_is_euclidean_self_dual",
     "every_C_perp_basis_vector_lies_in_C, whose condition it repeats"),
    ("xW_is_a_nonzero_codeword", "xW_has_hamming_weight_13"),
    ("xW_information_part_is_the_first_unit_vector",
     "codeword_re_derived_from_rref_without_using_x with "
     "pivot_columns_are_the_first_15"),
    ("derived_dimension_agrees_with_the_conjecture",
     "conjectured_parameters_p1_p1half_p3half_are_30_15_16"),
    ("code_is_not_MDS",
     "minimum_distance_upper_bound_beats_the_conjectured_distance, the "
     "Singleton value having been checked equal to the conjectured distance"),
]

# Only dependent when the depth >= 2 census line below it is present.
DEPENDENT_RESTATEMENTS_DEEP = [
    ("census_minimum_is_at_most_the_papers_bound_13",
     "census_minimum_equals_12_derived_here, and 12 <= 13"),
]


def main(argv):
    depth = CENSUS_DEPTH_DEFAULT
    for a in argv[1:]:
        if a.startswith("--census-depth="):
            depth = int(a.split("=", 1)[1])
        elif a == "--quick":
            depth = 2
        elif a == "--fast":
            depth = 4
        else:
            print("usage: verify.py [--census-depth=N] [--fast] [--quick]")
            print("  default depth 5 (~200 s here); --fast depth 4 (~3 s);")
            print("  --quick depth 2 (instant).  Everything except the census")
            print("  runs in well under a second at any depth.")
            return 2
    if not 1 <= depth <= P_K:
        print("census depth must lie in 1..%d" % P_K)
        return 2
    t0 = time.time()
    note("verify.py -- Etzion-Vardy-Yaakobi Conjecture 5 at p = 29")
    note("exact integer and F_29 arithmetic only; census depth %d" % depth)

    W = stage_well_formed()
    stage_conference(W)
    stage_minor(W)
    R, piv, G, k = stage_code(W)
    xW, w = stage_codeword(W, G, k)
    stage_refutation(len(W), k, w)
    bw = stage_census(W, R, piv, G, k, depth)

    total_words = P_P ** P_K - 1
    covered = sum(binom(P_K, m) * (P_P - 1) ** m for m in range(1, depth + 1))
    note("NOT RE-RUN: the exact minimum distance of C. The paper claims only "
         "d(C) <= 13 and that claim is fully verified above. The census "
         "settles the minimum weight over %d of the %d nonzero codewords "
         "(information weight <= %d); codewords of information weight >= %d "
         "are NOT enumerated, because depth %d costs about %d times depth %d."
         % (covered, total_words, depth, depth + 1, depth + 1,
            (P_K - depth) * (P_P - 1) // (depth + 1), depth))
    note("no exhaustive census over conference matrices of order 30 is "
         "attempted: the paper exhibits one matrix and claims nothing about "
         "the others.")
    note("NOT VERIFIED: the provenance sentence (representative 3 of Spence's "
         "order-30 conference two-graph file). This program is standard "
         "library only and fetches nothing; the paper states that the "
         "provenance is not used in the proof, and the checks above are "
         "self-contained.")
    # ---- Disclosure of the check inventory, with DERIVED counts. ----------
    # The headline count is what a referee reads first, so the lines that
    # carry no force of their own are named here rather than left to be
    # rediscovered: the machinery self-tests, and the lines that are logical
    # consequences of an earlier check and so cannot fail independently of it.
    ck_total = len(_RESULTS)
    recorded = set(nm for nm, ok in _RESULTS)
    machinery = [nm for nm in MACHINERY_SELF_TESTS if nm in recorded]
    restated_spec = list(DEPENDENT_RESTATEMENTS)
    if depth >= 2:
        restated_spec += DEPENDENT_RESTATEMENTS_DEEP
    restated = [(a, b) for a, b in restated_spec if a in recorded]
    independent = ck_total - len(machinery) - len(restated)
    note("MACHINERY SELF-TESTS (%d of the %d checks; they test this program, "
         "not the paper, and cannot fail for a well-formed input): %s. Each "
         "is a count checked against a closed form, or a second route that "
         "must agree with a first. A loop-bound or elimination error here "
         "reports FAIL, but a wrong matrix is caught by the Gram, minor, rank "
         "and codeword checks, not by these."
         % (len(machinery), ck_total, "; ".join(machinery)))
    note("DEPENDENT RESTATEMENTS (%d of the %d checks): each of the following "
         "is a logical consequence of an earlier PASS line, so, given that "
         "line, it cannot fail on its own -- %s. They are printed because "
         "each names the step of the argument it licenses, but the number of "
         "lines carrying force independent of every other line is %d, not "
         "%d: a referee should read the VERDICT count below as %d independent "
         "checks."
         % (len(restated), ck_total,
            "; ".join("%s (follows from %s)" % (a, b) for a, b in restated),
            independent, ck_total, independent))
    stale = ([nm for nm in MACHINERY_SELF_TESTS if nm not in recorded]
             + [a for a, b in restated_spec if a not in recorded])
    if stale:
        note("DISCLOSURE INVENTORY is out of date: it names %d check(s) that "
             "this run did not record, so the two counts above may be wrong: "
             "%s" % (len(stale), "; ".join(stale)))
    # ---- Switching insensitivity, DERIVED from P_X and P_XW. --------------
    # D W D is again a conference matrix, and every check above returns the
    # same verdict on it, exactly when the -1 entries of D are confined to the
    # coordinates where x and xW BOTH vanish: there x D = x and (x W) D = x W,
    # hence x (D W D) = x W entry for entry.  The coordinate set and the orbit
    # size are computed from the two printed vectors, not asserted.
    xf = fit(P_X, P_N)
    xWf = fit(P_XW, P_N)
    both_zero = [j + 1 for j in range(P_N)
                 if xf[j] % P_P == 0 and xWf[j] % P_P == 0]
    if len(both_zero) > 1:
        coords = (", ".join(str(c) for c in both_zero[:-1])
                  + " and " + str(both_zero[-1]))
    else:
        coords = ", ".join(str(c) for c in both_zero) or "(none)"
    note("KNOWN INSENSITIVITY: the checks are invariant under the switchings "
         "D W D, D = diag(+-1), whose -1 entries are confined to the "
         "coordinates where x and xW both vanish -- computed from the two "
         "vectors printed above, those are the %d coordinates %s, an orbit of "
         "2^%d = %d displays. Each orbit member is itself a conference matrix "
         "whose F_29 row span is a self-dual [%d,%d] code containing the same "
         "weight-%d word, so the theorem holds verbatim for every one of them; "
         "only the provenance sentence, which the paper excludes from the "
         "proof, distinguishes them. The decoded matrix is printed back above "
         "so the display itself can be compared character by character."
         % (len(both_zero), coords, len(both_zero), 2 ** len(both_zero),
            len(W), k, w))
    note("elapsed %.1f s" % (time.time() - t0))
    return verdict()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

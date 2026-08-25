#!/usr/bin/env python3
"""
Exact verifier for "A Nine-Element Counterexample to Sinclair's Front-Loading
Conjecture".  Standard library only; every decision is made
with exact integer / Fraction arithmetic.  No floating point anywhere.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs; these are NOT checks)
--------------------------------------------------------------------------
  (I1) The 5x9 integer matrix A of Theorem 1 (hard-coded below as A_PAPER).
       This is the exhibited object.  Corrupt it and this program must fail.
  (I2) The transcription of Sinclair's operator formula, eq. (2) of the paper:
           L_a^beta g_Y = sum over covers X < Y with a in Y\\X of g_X
           D_beta   g_Y = rk(Y) * sum over covers X < Y of (|Y\\X|/(n-|X|)) g_X
       on the rescaled Moebius-algebra basis g_X = beta_X e_X, together with
       the chain  C_0 = Q g_E,  C_{k+1} = D_beta C_k + sum_a L_a^beta C_k
       (and C_k^0 obtained by omitting D_beta).  A verifier cannot check a
       definition against Sinclair's paper; it is an input here.
  (I3) The statement of Conjecture C:  delta_k >= delta_{r-k} for k <= r/2,
       where delta_k = h_k^beta - h_k^0.
  (I4) The numeric assertions of the paper, used only as comparison targets:
       9 elements, rank 5, simple, 67 bases, 119 flats,
       W = (1,9,32,49,27,1), h^0 = (1,9,27,27,9,1),
       h^beta = (1,9,31,32,9,1), delta = (0,0,4,5,0,0),
       and the two printed Hilbert polynomials of the abstract.
  (I5) THE INDEX ORIENTATION, and it is load-bearing.  C_0 = Q g_E lives in
       the top degree r of the graded Moebius algebra and every operator
       lowers the degree by one, so C_k sits in H_{r-k}.  The paper indexes
       h_k by the NUMBER OF OPERATOR APPLICATIONS (h_k = dim C_k): its table
       prints W_{r-k} against h_k and its abstract prints
       1+9q+31q^2+32q^3+9q^4+q^5.  This verifier uses that convention.  The
       opposite reading, indexing by the DEGREE of the graded piece
       (h_j = dim C_{r-j}), reverses delta to (0,0,5,4,0,0), which SATISFIES
       Conjecture C.  Nothing internal pins the orientation: both readings fit
       inside the ambient dimensions, and h^0 is palindromic so it cannot
       distinguish them.  The program prints the reversed reading and its
       verdict, and shows that under the reversed reading Conjecture C is
       already false for simple rank-5 matroids on nine elements OTHER than
       M, so the reversed reading cannot be the one in which the conjecture
       was posed.  That evidence is no longer one hand-picked matrix.  Under
       the reversed reading the conjecture asserts that no defect is ever
       STRICTLY front-loaded (delta_k > delta_{r-k} for some k <= r/2), so
       every strictly front-loaded matroid falsifies it; section_controls
       therefore runs a DETERMINISTIC scan written out in full here
       (ORIENT_SCAN_CELLS x ORIENT_SCAN_VALUES applied to A_CTRL_FRONT -- no
       random source, no seed, no external search code) and reports how many
       of its members are simple of rank 5 on nine elements, and how many of
       those are strictly front-loaded.  A referee can regenerate the scan
       line for line and enlarge it by widening either list.  It is EVIDENCE
       about the orientation, not a derivation of it: the orientation must
       still be confirmed against Sinclair's preprint, and the closing
       "NOT RE-RUN:" lines say so.

--------------------------------------------------------------------------
DERIVED HERE (the actual checks; computed from A_PAPER alone)
--------------------------------------------------------------------------
  * rank of every one of the 2^9 column subsets, by two independent
    implementations (rational Gauss-Jordan and rank over a 31-bit prime
    field), plus rank submodularity over all 2^9 x 2^9 subset pairs;
  * loops / parallel pairs, hence simplicity; rank 5; 9 elements;
  * closures, the 119 flats, the Whitney numbers W_k, the covering relations,
    and the paper's proof step "X v a = Y iff a in Y\\X for a cover X < Y";
  * the characteristic polynomial by two independent routes (Whitney rank
    formula over 2^9 subsets vs. Moebius function of the flat lattice);
  * the number of bases, twice (rank oracle, and 126 exact 5x5 Bareiss
    determinants);  beta_X for every flat, and beta_E = 67*5!;
  * the operator matrices, the two chains, and h_k^0, h_k^beta with a
    two-sided CERTIFICATE for each: a nonzero h_k x h_k minor (lower bound,
    by fraction-free Bareiss) and explicitly verified annihilating
    functionals (upper bound), the generating sets for k <= 3 being all
    words of length k applied to g_E so the bound needs no induction.  The
    certificates do not avoid elimination: the rows of the minor, its
    columns and the annihilators are all SELECTED by the same incremental
    pivoting used by span_basis.  What is independent is the WITNESS, not
    the search that found it (see the note in section_certificates);
  * delta_k, and whether Conjecture C's inequality holds at each k <= r/2;
  * log-concavity of (1,9,31,32,9,1) (the paper's Remark);
  * two CONTROL matroids run through the identical pipeline, for which the
    conjecture test must report HOLDS -- so the refutation test is two-sided
    and not a gate that always fires.

Nothing in the paper's claim is left out: there is no census or minimality
claim to reproduce ("No claim of minimality is made", Remark 2).  This file is
the ancillary verifier referred to in the paper's "Exact verification" section.

FALSIFIABILITY.  Mutations were run against single entries of A, whole columns,
forced loops and parallel pairs, rank drops, shape changes, each claimed
constant, the transcription of D_beta, and each internal routine in turn.  The
honest summary, after an adversarial review:

  * a single-entry corruption of A that CHANGES THE MATROID is caught by the
    load-bearing lines: the paper's row 3, column 8 entry 2 -> 1 fails 17 of the
    69 checks and 2 -> 3 fails 13, in both cases including
    beta_hilbert_function, differential_defect_delta, certificate_hbeta_k2 and
    k3, CONJECTURE_C_VIOLATED_AT_k2 and CONJECTURE_C_IS_REFUTED_BY_M;  row 1,
    column 6 entry 1 -> 0 fails 24, starting at no_parallel_pairs;
  * a corruption of A that does NOT change the column matroid -- rescaling a
    column, permuting two columns, A[5][8] 3 -> 4 or A[5][9] 1 -> 2 (all
    verified to leave the rank function on all 2^9 subsets identical) -- leaves
    every mathematical check passing, correctly, because Theorem 1 is a
    statement about the matroid;  what pins the PRINTED matrix in that case is
    only the transcription line last_four_columns_match_paper, which compares
    A_PAPER against a second, independent hand transcription of the same
    printed matrix (PAPER_EXTRA_COLS) and so is a double-entry guard, not a
    mathematical check;
  * dropping the rk(Y) factor from D_beta, or replacing rk(Y) by rk(X), cannot
    matter, because rk is constant on each graded piece, so D_beta|H_k is only
    rescaled and the span is unchanged;  the program recomputes h^beta without
    the factor and prints whether it changed.  This matroid's chain is also
    unchanged if one L_a, or D_beta at the single step k=3, is dropped -- but
    h^0 does change, so the mutation is still caught by delta;
  * two checks were found by review to be STRUCTURALLY unfalsifiable and have
    been strengthened: "every_h_k_at_most_dim_H_(r-k)" (the computed half is
    forced by span_basis stopping at `ambient` pivots; it now also requires the
    PAPER's claimed h vectors to fit inside the COMPUTED Whitney numbers) and
    "beta_X_positive_..." (every flat of every matroid has a basis, so
    positivity alone can never fail; beta_X is now also recounted with the
    second, independent rank implementation).  Several further checks are
    theorems about any matrix rather than facts about this one --
    rank_submodular_on_all_subset_pairs, rank_unit_increase_and_bounded,
    closure_is_idempotent, flats_closed_under_intersection,
    flat_lattice_is_graded, delta_is_nonnegative, C_k^0_is_a_subspace_of
    C_k^beta, word_generators_and_iterative_chain_agree and
    the_two_forms_of_conjecture_C_agree_at_every_k.  They are retained as
    implementation self-tests (each does fail under an injected bug in the
    routine it exercises), but they are not evidence for the paper's claim.

WHAT THE PRINTED TOTAL CONTAINS.  Every line printed as PASS is counted in the
total, self-tests included, so the total is NOT a count of evidence about M.
The closing "ACCOUNTING:" lines split it four ways, and the four numbers are
COUNTED from the check names actually executed rather than hardcoded:
  * evidence about M itself;
  * evidence about the two control matroids and the orientation scan (these
    bear on input I5 and on the two-sidedness of the predicate, not on M);
  * implementation self-tests: the nine theorems-about-any-matrix listed above,
    plus the accounting check itself;
  * lines that only restate the input matrix and so carry no information about
    the matroid: matrix_shape_5x9, matrix_entries_are_integers,
    first_five_columns_are_identity, ground_set_has_9_elements, and
    last_four_columns_match_paper, which is a double-entry transcription guard
    comparing one hand transcription of the printed matrix against another.
A name declared non-evidence that fails to run makes the check
declared_non_evidence_check_names_all_ran FAIL, so the split cannot drift out
of step with the checks it describes.  The classification is the one argued
above and is deliberately conservative: a referee may judge further lines
(bottom_flat_is_empty and ground_set_is_the_unique_top_flat, for instance) to
be theorems about any matroid as well.

Usage:  python3 verify.py        (exit 0 iff every check passes)
"""

import sys
from fractions import Fraction as Fr
from itertools import combinations

_FAILED = []
_NAMES = []
_N = 0


def check(name, cond, detail=""):
    """Record one falsifiable check."""
    global _N
    _N += 1
    _NAMES.append(name)
    cond = bool(cond)
    if not cond:
        _FAILED.append(name)
    line = ("PASS " if cond else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    return cond


def info(*parts):
    print("#", *parts)


# --------------------------------------------------------------------------
# Exact linear algebra.  Two independent rank implementations on purpose.
# --------------------------------------------------------------------------

P_MOD = (1 << 31) - 1          # 2147483647, prime


def rref(rows):
    """Reduced row echelon form over Q.  Returns (basis rows, pivot columns)."""
    m = [list(map(Fr, rw)) for rw in rows]
    if not m:
        return [], []
    ncols = len(m[0])
    piv = []
    top = 0
    for c in range(ncols):
        sel = None
        for i in range(top, len(m)):
            if m[i][c] != 0:
                sel = i
                break
        if sel is None:
            continue
        m[top], m[sel] = m[sel], m[top]
        pv = m[top][c]
        if pv != 1:
            m[top] = [x / pv for x in m[top]]
        for i in range(len(m)):
            if i != top and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[top])]
        piv.append(c)
        top += 1
        if top == len(m):
            break
    return m[:top], piv


def rank_q(rows):
    """Rank over Q by rational Gauss-Jordan (implementation 1)."""
    return len(rref(rows)[0])


def rank_modp(rows, p=P_MOD):
    """Rank of an integer matrix over GF(p) (implementation 2, independent)."""
    m = [[int(x) % p for x in rw] for rw in rows]
    if not m:
        return 0
    ncols = len(m[0])
    top = 0
    for c in range(ncols):
        sel = None
        for i in range(top, len(m)):
            if m[i][c]:
                sel = i
                break
        if sel is None:
            continue
        m[top], m[sel] = m[sel], m[top]
        inv = pow(m[top][c], p - 2, p)
        m[top] = [(x * inv) % p for x in m[top]]
        for i in range(len(m)):
            if i != top and m[i][c]:
                f = m[i][c]
                m[i] = [(a - f * b) % p for a, b in zip(m[i], m[top])]
        top += 1
        if top == len(m):
            break
    return top


def det_bareiss(mat):
    """Exact determinant by fraction-free Bareiss elimination.

    Used as an INDEPENDENT witness of nonsingularity (it never divides by a
    pivot the way rref does).  Accepts Fractions; scales rows to clear
    denominators first, then works over Z, and divides the scaling out.
    """
    n = len(mat)
    if n == 0:
        return Fr(1)
    scale = Fr(1)
    m = []
    for rw in mat:
        assert len(rw) == n
        d = 1
        for x in rw:
            d = d * Fr(x).denominator // _gcd(d, Fr(x).denominator)
        scale *= Fr(1, d)
        m.append([int(Fr(x) * d) for x in rw])
    sign = 1
    prev = 1
    for k in range(n - 1):
        if m[k][k] == 0:
            sel = None
            for i in range(k + 1, n):
                if m[i][k] != 0:
                    sel = i
                    break
            if sel is None:
                return Fr(0)
            m[k], m[sel] = m[sel], m[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) // prev
            m[i][k] = 0
        prev = m[k][k]
    return Fr(sign * m[n - 1][n - 1]) * scale


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def kernel_basis(rows, ncols):
    """Basis of {x in Q^ncols : R x = 0} for R with the given rows."""
    R, piv = rref(rows) if rows else ([], [])
    free = [c for c in range(ncols) if c not in set(piv)]
    ker = []
    for f in free:
        v = [Fr(0)] * ncols
        v[f] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -R[i][f]
        ker.append(v)
    return ker


# --------------------------------------------------------------------------
# The matroid of a matrix: subset ranks, closures, flats, covers.
# --------------------------------------------------------------------------

def popcount(x):
    return bin(x).count("1")


class Matroid(object):
    """Column matroid of an integer matrix, built from scratch by elimination."""

    def __init__(self, A):
        self.A = [list(rw) for rw in A]
        self.nrows = len(self.A)
        self.n = len(self.A[0])
        self.cols = [tuple(self.A[i][j] for i in range(self.nrows))
                     for j in range(self.n)]
        self.RK = [0] * (1 << self.n)
        for mask in range(1 << self.n):
            self.RK[mask] = rank_q([self.cols[j] for j in range(self.n)
                                    if mask >> j & 1])
        self.r = self.RK[(1 << self.n) - 1]
        self.CL = [0] * (1 << self.n)
        for mask in range(1 << self.n):
            c = mask
            rk = self.RK[mask]
            for a in range(self.n):
                if not (mask >> a & 1) and self.RK[mask | (1 << a)] == rk:
                    c |= 1 << a
            self.CL[mask] = c
        self.flats = sorted(set(self.CL))
        self.byrank = dict((k, []) for k in range(self.r + 1))
        for X in self.flats:
            self.byrank[self.RK[X]].append(X)
        self.W = [len(self.byrank[k]) for k in range(self.r + 1)]
        self.idx = dict((k, dict((X, i) for i, X in enumerate(self.byrank[k])))
                        for k in self.byrank)
        # covers: up[X] = [(Y, Y\X) : X covered by Y]
        self.up = dict((X, []) for X in self.flats)
        for k in range(self.r):
            for X in self.byrank[k]:
                for Y in self.byrank[k + 1]:
                    if X & Y == X and X != Y:
                        self.up[X].append((Y, Y & ~X))

    def n_bases(self):
        full = self.r
        return sum(1 for S in combinations(range(self.n), full)
                   if self.RK[sum(1 << j for j in S)] == full)

    def beta(self, X):
        """Number of ORDERED bases of the restriction to the flat X."""
        k = self.RK[X]
        els = [j for j in range(self.n) if X >> j & 1]
        cnt = sum(1 for S in combinations(els, k)
                  if self.RK[sum(1 << j for j in S)] == k)
        f = 1
        for i in range(2, k + 1):
            f *= i
        return cnt * f


# --------------------------------------------------------------------------
# Sinclair's operators, eq. (2), in the g-basis.  Both lower degree by one.
# op = integer a  ->  L_a^beta ;  op = 'D'  ->  D_beta.
# --------------------------------------------------------------------------

def apply_op(M, v, k, op):
    """Apply one operator to v in H_k (coords indexed by rank-k flats)."""
    out = [Fr(0)] * M.W[k - 1]
    ik = M.idx[k]
    for X in M.byrank[k - 1]:
        s = Fr(0)
        if op == 'D':
            denom = M.n - popcount(X)
            for (Y, d) in M.up[X]:
                c = v[ik[Y]]
                if c:
                    s += Fr(M.RK[Y] * popcount(d), denom) * c
        else:
            bit = 1 << op
            for (Y, d) in M.up[X]:
                if d & bit:
                    s += v[ik[Y]]
        out[M.idx[k - 1][X]] = s
    return out


def span_basis(vecs, ambient):
    """Row-reduced basis of the span, built incrementally with early exit.

    Every stored pivot row is kept zero at all other pivot columns, so the
    single reduction pass below is exact whatever order the pivots arose in.
    """
    pivots = {}
    order = []
    for v in vecs:
        w = list(v)
        for c in list(pivots):
            if w[c]:
                f = w[c]
                w = [x - f * y for x, y in zip(w, pivots[c])]
        c = None
        for i, x in enumerate(w):
            if x:
                c = i
                break
        if c is None:
            continue
        pv = w[c]
        if pv != 1:
            w = [x / pv for x in w]
        for c2 in list(pivots):
            b = pivots[c2]
            if b[c]:
                f = b[c]
                pivots[c2] = [x - f * y for x, y in zip(b, w)]
        pivots[c] = w
        order.append(c)
        if len(order) == ambient:
            break
    return [pivots[c] for c in order]


def top_vector(M):
    """g_E, the generator of C_0, as a vector in H_r (which is 1-dimensional)."""
    return [Fr(1)] * M.W[M.r]


def chain_iterative(M, use_D):
    """h_k = dim C_k, C_{k+1} = D C_k + sum_a L_a C_k, C_0 = Q g_E.

    Returns (h, gens) where gens[k] is the generating set actually used for
    C_k (images of a basis of C_{k-1}); gens[0] = [g_E].
    """
    ops = list(range(M.n)) + (['D'] if use_D else [])
    C = [top_vector(M)]
    h = [len(C)]
    gens = [[top_vector(M)]]
    for k in range(M.r, 0, -1):
        g = []
        for v in C:
            for op in ops:
                g.append(apply_op(M, v, k, op))
        C = span_basis(g, M.W[k - 1])
        gens.append(g)
        h.append(len(C))
    return h, gens


def words_generators(M, use_D, kmax):
    """All words of length k in the operators applied to g_E, for k<=kmax.

    span{words of length k} = C_k exactly (the operators are linear), so this
    is an induction-free generating set for C_k.
    """
    ops = list(range(M.n)) + (['D'] if use_D else [])
    cur = [top_vector(M)]
    out = {0: [top_vector(M)]}
    for k in range(M.r, M.r - kmax, -1):
        nxt = []
        for v in cur:
            for op in ops:
                nxt.append(apply_op(M, v, k, op))
        out[M.r - k + 1] = nxt
        cur = nxt
    return out


# --------------------------------------------------------------------------
# Two-sided dimension certificates.  Elimination is used to SEARCH for the
# witnesses and is not trusted to establish them: independent_indices picks the
# rows, rref picks the columns and kernel_basis proposes the annihilators, but
# the lower bound is then established by a fraction-free Bareiss determinant of
# the selected minor, and the upper bound by dotting every proposed annihilator
# against every generator and by a second Bareiss minor witnessing that the
# annihilators are independent.  A bug in the elimination can therefore lose a
# witness (and fail the check) but cannot manufacture one.
# --------------------------------------------------------------------------

def independent_indices(vecs, ambient):
    """Greedy: indices of the ORIGINAL vectors forming a basis of span(vecs).

    Same pivot maintenance as span_basis, so no column is ever reintroduced.
    """
    pivots = {}
    keep = []
    for i, v in enumerate(vecs):
        w = list(v)
        for c in list(pivots):
            if w[c]:
                f = w[c]
                w = [x - f * y for x, y in zip(w, pivots[c])]
        c = next((j for j, x in enumerate(w) if x), None)
        if c is None:
            continue
        pv = w[c]
        w = [x / pv for x in w]
        for c2 in list(pivots):
            b = pivots[c2]
            if b[c]:
                f = b[c]
                pivots[c2] = [x - f * y for x, y in zip(b, w)]
        pivots[c] = w
        keep.append(i)
        if len(keep) == ambient:
            break
    return keep


def minor_is_nonsingular(vecs, rows_idx):
    """det of a square submatrix of the chosen rows, by Bareiss; True if != 0."""
    rows = [vecs[i] for i in rows_idx]
    _, piv = rref(rows)
    if len(piv) != len(rows):
        return False, None
    sub = [[Fr(rows[i][c]) for c in piv] for i in range(len(rows))]
    d = det_bareiss(sub)
    return d != 0, d


def certify_dim(gens, ambient):
    """Return (lower, upper, det, nker) with lower <= dim span(gens) <= upper.

    lower: size of an explicitly nonsingular minor among the generators.
    upper: ambient minus the number of verified independent annihilators.
    """
    keep = independent_indices(gens, ambient)
    ok, det = minor_is_nonsingular(gens, keep)
    lower = len(keep) if ok else -1
    ker = kernel_basis(gens, ambient)
    good = 0
    for x in ker:
        if all(sum(a * b for a, b in zip(g, x) if a) == 0 for g in gens):
            good += 1
    # The upper bound ambient - |ker| is only licensed when EVERY annihilator
    # was verified against EVERY generator and the annihilators are provably
    # independent; otherwise report no usable bound (ambient + 1) so the
    # calling check fails instead of accepting a bound it cannot justify.
    indep = True
    if ker:
        indep, _ = minor_is_nonsingular(ker, list(range(len(ker))))
    if good == len(ker) and indep:
        upper = ambient - len(ker)
    else:
        upper = ambient + 1
    return lower, upper, det, good


# --------------------------------------------------------------------------
# INPUT (I1): the exhibited object, Theorem 1.
# --------------------------------------------------------------------------

A_PAPER = [
    [1, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 0, 1, 1, 2, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1, 0, 1, 3, 1],
]

# INPUT (I4): the paper's assertions, used only as comparison targets.
PAPER_EXTRA_COLS = [(1, 0, 1, 0, 0), (0, 1, 1, 0, 1),
                    (0, 1, 2, 0, 3), (0, 0, 0, 1, 1)]
PAPER_N = 9
PAPER_RANK = 5
PAPER_BASES = 67
PAPER_FLATS = 119
PAPER_W = [1, 9, 32, 49, 27, 1]
PAPER_H0 = [1, 9, 27, 27, 9, 1]
PAPER_HB = [1, 9, 31, 32, 9, 1]
PAPER_DELTA = [0, 0, 4, 5, 0, 0]
PAPER_H0_POLY = "1+9q+27q^2+27q^3+9q^4+q^5"
PAPER_HB_POLY = "1+9q+31q^2+32q^3+9q^4+q^5"


def poly_string(h):
    """Render a Hilbert vector the way the paper's abstract prints it."""
    out = []
    for k, c in enumerate(h):
        if c == 0:
            continue
        if k == 0:
            out.append(str(c))
        else:
            coef = "" if c == 1 else str(c)
            out.append(coef + ("q" if k == 1 else "q^%d" % k))
    return "+".join(out)


def mask_str(X, n):
    return "{" + ",".join(str(j + 1) for j in range(n) if X >> j & 1) + "}"


def conjecture_C_violations(delta, r):
    """Indices k <= r/2 where Conjecture C's inequality delta_k >= delta_{r-k}
    FAILS.  Returns a list; empty means the conjecture holds for this matroid."""
    return [k for k in range(r + 1) if 2 * k <= r and delta[k] < delta[r - k]]


def char_poly_subsets(M):
    """chi_M(t) = sum_{S subset E} (-1)^{|S|} t^{r-rk(S)}   (Whitney rank formula).
    Returned as a list of coefficients indexed by the power of t."""
    c = [0] * (M.r + 1)
    for mask in range(1 << M.n):
        c[M.r - M.RK[mask]] += (-1) ** popcount(mask)
    return c


def char_poly_moebius(M):
    """chi_M(t) = sum_{X flat} mu(bottom, X) t^{r-rk(X)} via the Moebius
    function of the lattice of flats -- an independent second route."""
    mu = {}
    for X in M.flats:            # M.flats is sorted, so subsets come first
        if M.RK[X] == 0:
            mu[X] = 1
            continue
        s = 0
        for Y in M.flats:
            if Y != X and Y & X == Y:
                s += mu[Y]
        mu[X] = -s
    c = [0] * (M.r + 1)
    for X in M.flats:
        c[M.r - M.RK[X]] += mu[X]
    return c, mu


def row_space_basis_indices(M):
    """Indices of a maximal independent set of ROWS of A.  Deleting the other
    rows changes neither the row space nor the column matroid, and makes the
    r x r minors square, so the determinant route works at any rank."""
    return independent_indices([[Fr(x) for x in rw] for rw in M.A], M.r)


def moebius_from_subsets(M):
    """mu(bottom, X) = sum over S subset X with cl(S)=X of (-1)^|S|.

    Computed from the closure map on all 2^n subsets, i.e. without ever using
    the lattice recursion -- an independent value to compare against.
    """
    g = dict((X, 0) for X in M.flats)
    for mask in range(1 << M.n):
        g[M.CL[mask]] += (-1) ** popcount(mask)
    return g


def cover_lemma_ok(M):
    """The paper's proof step: for a cover X < Y and an atom a,
    cl(X u a) = Y  iff  a in Y \\ X.  Checked over every (cover, atom) pair."""
    bad = 0
    ncov = 0
    for X in M.flats:
        for (Y, d) in M.up[X]:
            ncov += 1
            for a in range(M.n):
                join_is_Y = (M.CL[X | (1 << a)] == Y)
                in_diff = bool(d >> a & 1)
                if join_is_Y != in_diff:
                    bad += 1
    return bad, ncov


def connected_separator(M):
    """A 1-separation of M, or None if M is connected.

    M = M|S (+) M|(E-S) iff rk(S) + rk(E-S) = rk(E) for a proper nonempty S.
    A counterexample that were a direct sum would be a weaker object (it would
    reduce to a statement about a smaller matroid), so this is reported.
    """
    full = (1 << M.n) - 1
    for S in range(1, full):
        if M.RK[S] + M.RK[full ^ S] == M.r:
            return S
    return None


def lattice_is_graded(M):
    """Every flat pair X < Y with rk(Y)-rk(X) >= 2 has an intermediate flat
    of rank rk(X)+1 inside Y."""
    bad = 0
    for X in M.flats:
        for Y in M.flats:
            if X != Y and X & Y == X and M.RK[Y] - M.RK[X] >= 2:
                if not any(Z & Y == Z for (Z, _) in M.up[X]):
                    bad += 1
    return bad


def section_object(M):
    """(1) The exhibited object is well formed and is what the paper says."""
    info("--- the exhibited object, decoded and printed back ---")
    for i, rw in enumerate(M.A):
        info("A row %d:" % (i + 1), " ".join("%2d" % x for x in rw))
    for j, c in enumerate(M.cols):
        info("column %d = %s" % (j + 1, str(c)))

    check("matrix_shape_5x9",
          len(M.A) == 5 and all(len(rw) == 9 for rw in M.A),
          "%dx%d" % (len(M.A), len(M.A[0])))
    check("matrix_entries_are_integers",
          all(isinstance(x, int) for rw in M.A for x in rw))
    ident = [tuple(1 if i == j else 0 for i in range(5)) for j in range(5)]
    check("first_five_columns_are_identity", M.cols[:5] == ident)
    check("last_four_columns_match_paper", list(M.cols[5:]) == PAPER_EXTRA_COLS,
          str(list(M.cols[5:])))
    check("ground_set_has_9_elements", M.n == PAPER_N, "n=%d" % M.n)

    rfull = rank_q(M.cols)
    check("rank_is_5", rfull == PAPER_RANK, "rk(A)=%d" % rfull)
    check("rank_agrees_over_GF(2^31-1)",
          rank_modp([[M.A[i][j] for j in range(M.n)] for i in range(M.nrows)])
          == rfull)

    loops = [j for j in range(M.n) if M.RK[1 << j] == 0]
    check("no_loops", not loops, "zero columns: %s" % loops)
    par = []
    for a in range(M.n):
        for b in range(a + 1, M.n):
            if M.RK[(1 << a) | (1 << b)] < 2:
                par.append((a + 1, b + 1))
    check("no_parallel_pairs", not par, "parallel: %s" % par)
    check("M_is_simple_rank_5_on_9_elements",
          (not loops) and (not par) and M.r == 5 and M.n == 9)
    check("every_singleton_is_a_rank_1_flat",
          all(M.CL[1 << j] == 1 << j for j in range(M.n)) and M.W[1] == M.n,
          "W_1=%d" % M.W[1])
    sep = connected_separator(M)
    check("M_is_connected_not_a_direct_sum", sep is None,
          "no 1-separation" if sep is None else "separator %s"
          % mask_str(sep, M.n))


def section_rank_oracle(M):
    """The rank function everything else is built on, validated two ways."""
    full = (1 << M.n) - 1
    disagree = []
    for mask in range(1 << M.n):
        rows = [[M.A[i][j] for j in range(M.n) if mask >> j & 1]
                for i in range(M.nrows)]
        rp = rank_modp(rows) if mask else 0
        if rp != M.RK[mask]:
            disagree.append(mask)
    check("subset_ranks_agree_two_implementations", not disagree,
          "512 subsets, %d disagreements" % len(disagree))

    bad = 0
    for mask in range(1 << M.n):
        if not (0 <= M.RK[mask] <= min(popcount(mask), M.r)):
            bad += 1
        for a in range(M.n):
            if not (mask >> a & 1):
                d = M.RK[mask | (1 << a)] - M.RK[mask]
                if d not in (0, 1):
                    bad += 1
    check("rank_unit_increase_and_bounded", bad == 0, "%d violations" % bad)

    bad = 0
    for S in range(1 << M.n):
        rs = M.RK[S]
        for T in range(S, 1 << M.n):
            if M.RK[S | T] + M.RK[S & T] > rs + M.RK[T]:
                bad += 1
    check("rank_submodular_on_all_subset_pairs", bad == 0,
          "%d ordered-pair violations among 2^9 x 2^9" % bad)

    check("closure_is_idempotent",
          all(M.CL[M.CL[mask]] == M.CL[mask] for mask in range(1 << M.n)))
    check("closure_contains_and_preserves_rank",
          all(M.CL[m] & m == m and M.RK[M.CL[m]] == M.RK[m]
              for m in range(1 << M.n)))
    fl = set(M.flats)
    check("flats_closed_under_intersection",
          all((X & Y) in fl for X in M.flats for Y in M.flats))
    check("ground_set_is_the_unique_top_flat",
          M.W[M.r] == 1 and M.byrank[M.r] == [full])
    check("bottom_flat_is_empty", M.byrank[0] == [0])


def section_lattice(M):
    """Flats, Whitney numbers, covers, characteristic polynomial."""
    info("--- flat lattice ---")
    info("W =", M.W, " total flats =", len(M.flats))
    for k in range(M.r + 1):
        if M.W[k] <= 12:
            info("  rank %d flats:" % k,
                 " ".join(mask_str(X, M.n) for X in M.byrank[k]))
    check("whitney_numbers_second_kind", M.W == PAPER_W, str(M.W))
    check("total_number_of_flats_is_119", len(M.flats) == PAPER_FLATS,
          str(len(M.flats)))

    c1 = char_poly_subsets(M)
    c2, mu = char_poly_moebius(M)
    info("chi_M(t) =", " ".join("%+d t^%d" % (c1[k], k)
                               for k in range(M.r, -1, -1) if c1[k]))
    check("characteristic_polynomial_two_routes_agree", c1 == c2,
          "Whitney rank formula %s vs Moebius %s" % (c1, c2))
    check("moebius_lattice_recursion_matches_subset_formula",
          moebius_from_subsets(M) == mu,
          "%d flats, mu from cl() on 2^%d subsets vs mu from the lattice"
          % (len(M.flats), M.n))

    # The covering relation itself, cross-checked against the independent
    # characterisation "the flats covering X are exactly the cl(X u {a})".
    # Without this, a MISSING cover would simply never be tested below.
    miss = 0
    for X in M.flats:
        mine = set(Y for (Y, _) in M.up[X])
        theirs = set(M.CL[X | (1 << a)] for a in range(M.n)
                     if not (X >> a & 1))
        theirs.discard(X)
        if mine != theirs:
            miss += 1
    check("cover_set_equals_the_closures_cl(X_union_a)", miss == 0,
          "%d of %d flats with a mismatched cover set" % (miss, len(M.flats)))

    bad, ncov = cover_lemma_ok(M)
    info("covering pairs in the lattice:", ncov)
    check("cover_join_lemma_X_join_a_equals_Y_iff_a_in_Y_minus_X",
          bad == 0 and ncov > 0,
          "%d covers x 9 atoms, %d mismatches" % (ncov, bad))
    check("flat_lattice_is_graded", lattice_is_graded(M) == 0)

    nb = M.n_bases()
    # Independent route: an r-subset of columns is a basis iff the r x r minor
    # on any basis of the row space is nonsingular.  Bareiss determinants.
    rowbasis = row_space_basis_indices(M)
    nb2 = 0
    nsub = 0
    for S in combinations(range(M.n), M.r):
        nsub += 1
        sub = [[M.A[i][j] for j in S] for i in rowbasis]
        if det_bareiss(sub) != 0:
            nb2 += 1
    check("number_of_bases_is_67", nb == PAPER_BASES, "%d" % nb)
    check("number_of_bases_confirmed_by_exact_determinants",
          nb2 == nb and len(rowbasis) == M.r,
          "%d of %d %d-subsets nonsingular on rows %s"
          % (nb2, nsub, M.r, [i + 1 for i in rowbasis]))

    betas = dict((X, M.beta(X)) for X in M.flats)
    # POSITIVITY ALONE IS VACUOUS: every flat of a matroid has a basis, so
    # beta_X >= 1 holds for every matrix that could ever be substituted here.
    # The value is therefore also recounted by the second, independent rank
    # implementation (GF(2^31-1) elimination on the columns of the flat rather
    # than the rational rank oracle that built M.RK), so that a wrong beta --
    # or a wrong flat -- makes this line fail instead of certifying itself.
    bad_beta = []
    for X in M.flats:
        k = M.RK[X]
        els = [j for j in range(M.n) if X >> j & 1]
        cnt = 0
        for S in combinations(els, k):
            rows = [[M.A[i][j] for j in S] for i in range(M.nrows)]
            if (rank_modp(rows) if k else 0) == k:
                cnt += 1
        f = 1
        for i in range(2, k + 1):
            f *= i
        if betas[X] <= 0 or betas[X] != cnt * f:
            bad_beta.append((mask_str(X, M.n), betas[X], cnt * f))
    check("beta_X_positive_and_reconfirmed_for_every_flat",
          not bad_beta,
          "min beta_X = %d over %d flats, %d mismatches %s"
          % (min(betas.values()), len(M.flats), len(bad_beta), bad_beta[:3]))
    fact = 1
    for i in range(2, M.r + 1):
        fact *= i
    check("beta_E_equals_bases_times_r_factorial",
          betas[(1 << M.n) - 1] == nb * fact,
          "beta_E=%d = %d*%d" % (betas[(1 << M.n) - 1], nb, fact))


def section_hilbert(M):
    """The two chains, with a two-sided certificate for every dimension."""
    info("--- Sinclair's chains on the graded Moebius algebra ---")
    h0, g0 = chain_iterative(M, False)
    hb, gb = chain_iterative(M, True)
    info("h^0    =", h0)
    info("h^beta =", hb)
    check("apolar_hilbert_function_h0", h0 == PAPER_H0, str(h0))
    check("beta_hilbert_function", hb == PAPER_HB, str(hb))
    check("h0_polynomial_string_matches_abstract",
          poly_string(h0) == PAPER_H0_POLY, poly_string(h0))
    check("hbeta_polynomial_string_matches_abstract",
          poly_string(hb) == PAPER_HB_POLY, poly_string(hb))
    check("h0_is_palindromic",
          all(h0[k] == h0[M.r - k] for k in range(M.r + 1)), str(h0))
    check("both_chains_reach_degree_0_of_the_Moebius_algebra",
          h0[M.r] == 1 and hb[M.r] == 1,
          "h_r^0=%d, h_r^beta=%d" % (h0[M.r], hb[M.r]))
    # NOTE ON FALSIFIABILITY.  The first half of this check (the COMPUTED h_k
    # against the computed ambient dimension) cannot fail for any input:
    # span_basis stops at `ambient` pivots by construction, so h_k <= W_{r-k}
    # is forced by the code, not by the object.  It is retained only as an
    # assertion about the data structure.  The second half is the falsifiable
    # one: the PAPER's claimed Hilbert vectors must fit inside the ambient
    # graded pieces of the COMPUTED flat lattice, which fails as soon as a
    # corruption of A shrinks any Whitney number (e.g. W_2 = 30 < 32).
    check("every_h_k_at_most_dim_H_(r-k)",
          all(hb[k] <= M.W[M.r - k] and h0[k] <= M.W[M.r - k]
              and PAPER_HB[k] <= M.W[M.r - k] and PAPER_H0[k] <= M.W[M.r - k]
              for k in range(M.r + 1)),
          "W_{r-k} = %s" % [M.W[M.r - k] for k in range(M.r + 1)])

    # C_k^0 must sit inside C_k^beta (the beta chain has one extra operator).
    contained = True
    for k in range(1, M.r + 1):
        b0 = span_basis(g0[k], M.W[M.r - k])
        both = span_basis(list(b0) + list(gb[k]), M.W[M.r - k])
        if len(both) != hb[k]:
            contained = False
    check("C_k^0_is_a_subspace_of_C_k^beta", contained)

    # Honest disclosure of what the answer is NOT sensitive to: rk(Y) is
    # constant on H_k, so D_beta|H_k is a nonzero scalar multiple of the same
    # operator without that factor and the spans are unchanged.  Recomputed
    # here rather than argued.
    hb2, _ = chain_iterative_no_rank_factor(M)
    info("h^beta recomputed with the rk(Y) factor dropped from D_beta:", hb2,
         "-- unchanged" if hb2 == hb else "-- CHANGED")
    return h0, hb, g0, gb


def chain_iterative_no_rank_factor(M):
    """Same chain with rk(Y) omitted from D_beta (see the note above)."""
    C = [top_vector(M)]
    h = [1]
    gens = [[top_vector(M)]]
    for k in range(M.r, 0, -1):
        g = []
        for v in C:
            for a in range(M.n):
                g.append(apply_op(M, v, k, a))
            out = [Fr(0)] * M.W[k - 1]
            for X in M.byrank[k - 1]:
                s = Fr(0)
                for (Y, d) in M.up[X]:
                    s += Fr(popcount(d), M.n - popcount(X)) * v[M.idx[k][Y]]
                out[M.idx[k - 1][X]] = s
            g.append(out)
        C = span_basis(g, M.W[k - 1])
        gens.append(g)
        h.append(len(C))
    return h, gens


def section_certificates(M, h0, hb, g0, gb, kmax_words=3):
    """Re-derive each h_k by explicit certificates.

    For k <= kmax_words the generating set is ALL words of length k applied to
    g_E, so span(gens) = C_k with no reliance on the iterative construction;
    the annihilators then bound h_k from above unconditionally.  For larger k
    the claimed value equals dim H_{r-k}, so the ambient dimension is the
    upper bound and only the nonsingular minor is needed.
    """
    info("--- dimension certificates ---")
    info("CERTIFICATE SCOPE: the load-bearing dimensions 31 (k=2) and 32 (k=3)")
    info("CERTIFICATE SCOPE: do not rely on the ITERATIVE chain -- their")
    info("CERTIFICATE SCOPE: generators are all words of length k applied to")
    info("CERTIFICATE SCOPE: g_E -- but they DO use the same elimination as")
    info("CERTIFICATE SCOPE: span_basis to SEARCH for their witnesses:")
    info("CERTIFICATE SCOPE: independent_indices selects the rows of the minor,")
    info("CERTIFICATE SCOPE: rref selects its columns, kernel_basis proposes")
    info("CERTIFICATE SCOPE: the annihilators.  Only the WITNESSES are verified")
    info("CERTIFICATE SCOPE: independently: a nonzero Bareiss determinant for")
    info("CERTIFICATE SCOPE: the lower bound, and every annihilator dotted")
    info("CERTIFICATE SCOPE: against every generator plus a Bareiss minor for")
    info("CERTIFICATE SCOPE: their independence for the upper bound.  A bug in")
    info("CERTIFICATE SCOPE: the elimination can lose a witness, not invent one.")
    kmax_words = min(kmax_words, M.r)
    w0 = words_generators(M, False, kmax_words)
    wb = words_generators(M, True, kmax_words)
    agree = True
    for k in range(1, kmax_words + 1):
        amb = M.W[M.r - k]
        if len(span_basis(w0[k], amb)) != h0[k]:
            agree = False
        if len(span_basis(wb[k], amb)) != hb[k]:
            agree = False
    check("word_generators_and_iterative_chain_agree", agree,
          "k <= %d, %d/%d words" % (kmax_words, len(w0[kmax_words]),
                                    len(wb[kmax_words])))
    # The certificates are compared with the value the PAPER claims, so that a
    # corrupted object makes them fail rather than silently certify some other
    # number.  They independently re-derive that value from both sides.
    for tag, claim, words, gens in (("h0", PAPER_H0, w0, g0),
                                    ("hbeta", PAPER_HB, wb, gb)):
        for k in range(1, M.r + 1):
            amb = M.W[M.r - k]
            gset = words[k] if k <= kmax_words else gens[k]
            lo, up, det, nann = certify_dim(gset, amb)
            if k > kmax_words:
                up = amb           # dim C_k <= dim H_(r-k); tight iff equal
            want = claim[k] if k < len(claim) else None
            ok = (want is not None and lo == want and up == want)
            check("certificate_%s_k%d_equals_%s" % (tag, k, want), ok,
                  "%d <= dim <= %d, minor det=%s, %d annihilators, ambient %d"
                  % (lo, up, det, nann, amb))


def section_refutation(M, h0, hb):
    """(3) THE LOAD-BEARING CHECK: the object violates Conjecture C."""
    info("--- the differential defect and Conjecture C ---")
    delta = [b - a for a, b in zip(h0, hb)]
    info("k        :", " ".join("%3d" % k for k in range(M.r + 1)))
    info("W_(r-k)  :", " ".join("%3d" % M.W[M.r - k] for k in range(M.r + 1)))
    info("h_k^0    :", " ".join("%3d" % x for x in h0))
    info("h_k^beta :", " ".join("%3d" % x for x in hb))
    info("delta_k  :", " ".join("%3d" % x for x in delta))
    check("differential_defect_delta", delta == PAPER_DELTA, str(delta))
    check("delta_is_nonnegative", all(d >= 0 for d in delta))
    check("conjecture_C_index_k2_is_in_range_k_le_r_over_2", 2 * 2 <= M.r,
          "2*k=4 <= r=%d" % M.r)
    if M.r != 5:
        # The paper's claim is about a rank-5 matroid; the indices below are
        # meaningless otherwise.  Fail loudly instead of raising.
        check("rank_is_5_so_the_k=2_vs_k=3_comparison_is_defined", False,
              "rk(M)=%d" % M.r)
        return delta

    viol = conjecture_C_violations(delta, M.r)
    for k in range(M.r + 1):
        if 2 * k <= M.r:
            holds = delta[k] >= delta[M.r - k]
            info("Conjecture C at k=%d: delta_%d=%d vs delta_%d=%d -> %s"
                 % (k, k, delta[k], M.r - k, delta[M.r - k],
                    "holds" if holds else "VIOLATED"))
    check("conjecture_C_holds_at_k0", delta[0] >= delta[M.r])
    check("conjecture_C_holds_at_k1", delta[1] >= delta[M.r - 1])
    check("CONJECTURE_C_VIOLATED_AT_k2", delta[2] < delta[3],
          "delta_2=%d < delta_3=%d" % (delta[2], delta[3]))
    check("equivalent_form_h2beta_less_than_h3beta", hb[2] < hb[3],
          "h_2^beta=%d < h_3^beta=%d" % (hb[2], hb[3]))
    check("CONJECTURE_C_IS_REFUTED_BY_M", viol == [2],
          "violating indices k <= r/2: %s" % viol)

    # ----------------------------------------------------------------------
    # DISCLOSURE (I5): THE INDEX ORIENTATION IS LOAD-BEARING AND IS AN INPUT.
    # The chain starts at g_E, which has Moebius-algebra degree r, and every
    # operator lowers the degree by one, so C_k sits in H_{r-k}.  The paper
    # indexes h_k by the NUMBER OF OPERATOR APPLICATIONS, h_k = dim C_k: its
    # table prints W_{r-k} against h_k, and its abstract prints
    # 1+9q+31q^2+32q^3+9q^4+q^5.  This verifier uses the paper's convention.
    # It is not derivable from the paper alone, and it matters: indexing
    # instead by the DEGREE of the graded piece (h_j = dim C_{r-j}) reverses
    # delta, and the reversed sequence SATISFIES Conjecture C.  Neither
    # reading is excluded by the ambient dimensions, and h^0 is palindromic,
    # so h^0 cannot pin the orientation either.  The one piece of evidence a
    # verifier can compute is in section_controls: under the reversed reading
    # Conjecture C is already false for a randomly found simple rank-5 matroid
    # on nine elements, so the reversed reading cannot be the one in which the
    # conjecture was posed.  A referee should still confirm the orientation
    # against Sinclair's preprint.
    rev = delta[::-1]
    info("ORIENTATION (input I5): h_k = dim C_k is the paper's indexing;")
    info("ORIENTATION: under the reversed (degree) indexing delta would read",
         rev, "with violating indices", conjecture_C_violations(rev, M.r),
         "-- i.e. the verdict is orientation-dependent and the orientation is")
    info("ORIENTATION: taken from the paper, not derived (see control 2).")
    # The paper's claim that the two forms of the conjecture are equivalent
    # (it follows from the palindromicity of h^0, itself checked above).
    check("the_two_forms_of_conjecture_C_agree_at_every_k",
          all((delta[k] >= delta[M.r - k]) == (hb[k] >= hb[M.r - k])
              for k in range(M.r + 1) if 2 * k <= M.r))
    return delta


def section_remark(hb):
    """(4) Remark 2: the beta-sequence is strictly log-concave."""
    info("--- Remark 2 ---")
    check("hbeta_has_no_internal_zeros", all(x > 0 for x in hb), str(hb))
    lc = [(k, hb[k] * hb[k], hb[k - 1] * hb[k + 1])
          for k in range(1, len(hb) - 1)]
    for k, sq, pr in lc:
        info("log-concavity at k=%d: %d^2 = %d > %d" % (k, hb[k], sq, pr))
    check("hbeta_is_strictly_log_concave", all(sq > pr for _, sq, pr in lc))


# --------------------------------------------------------------------------
# CONTROLS.  The refutation test above must not be a gate that always fires.
# These two simple rank-5 matroids on 9 elements go through the identical
# pipeline and the identical predicate, which must report "holds" for them.
# --------------------------------------------------------------------------

# Vandermonde columns (1,t,t^2,t^3,t^4), t = 1..9: the uniform matroid U(5,9).
A_CTRL_UNIFORM = [[t ** i for t in range(1, 10)] for i in range(5)]

# A second simple rank-5 matroid on 9 elements whose defect is nonzero and
# genuinely front-loaded (delta_2 > delta_3).  It is the seed of the
# deterministic scan below, so its provenance no longer matters: whatever search
# first produced it, the scan around it is written out in full here and can be
# rerun and enlarged by a referee with no random source and no seed.
A_CTRL_FRONT = [
    [1, 0, 0, 0, 0, 1, 0, 3, 0],
    [0, 1, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 1, 0, 1, 1],
    [0, 0, 0, 1, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 3, 0, 2, 1],
]

# THE ORIENTATION SCAN (widening the evidence for input I5).  Under the reversed
# (degree) index orientation, Conjecture C reads delta_{r-k} >= delta_k for
# k <= r/2, i.e. it asserts that NO simple matroid has a strictly front-loaded
# defect.  Hence every simple matroid with delta_k > delta_{r-k} for some
# k <= r/2 refutes the reversed reading outright.  The family scanned is the
# seed A_CTRL_FRONT together with every single-entry replacement of it given by
# the cross product of the two lists below (0-based (row, column), all inside
# the four non-identity columns, so the rank stays 5 and the identity block is
# untouched).  Widen either list to enlarge the scan; nothing else changes.
ORIENT_SCAN_CELLS = [(0, 7), (2, 7), (4, 5), (4, 8)]
ORIENT_SCAN_VALUES = [2, 4]


def is_simple(M):
    if any(M.RK[1 << j] == 0 for j in range(M.n)):
        return False
    for a in range(M.n):
        for b in range(a + 1, M.n):
            if M.RK[(1 << a) | (1 << b)] < 2:
                return False
    return True


def binom(n, k):
    c = 1
    for i in range(k):
        c = c * (n - i) // (i + 1)
    return c


def orientation_scan(seed):
    """Deterministically enumerate the family described at ORIENT_SCAN_CELLS.

    Returns (rows, members, excluders, skipped), where rows carries one entry
    (label, delta, reversed-violation indices) per matrix scanned.  A matrix
    that is not simple of rank 5 on nine elements has delta None and is counted
    in `skipped`, so it is REPORTED rather than dropped: the scan always prints
    its own denominator, and a scan that produced no members at all fails the
    check below instead of vanishing.
    """
    mats = [("seed", [list(rw) for rw in seed])]
    for (i, j) in ORIENT_SCAN_CELLS:
        for val in ORIENT_SCAN_VALUES:
            B = [list(rw) for rw in seed]
            B[i][j] = val
            if all(B != C for _, C in mats):
                mats.append(("A[%d][%d]=%d" % (i + 1, j + 1, val), B))
    rows = []
    members = 0
    excluders = 0
    skipped = 0
    for label, B in mats:
        MB = Matroid(B)
        if not (MB.r == 5 and MB.n == 9 and is_simple(MB)):
            skipped += 1
            rows.append((label, None, None))
            continue
        members += 1
        h0b, _ = chain_iterative(MB, False)
        hbb, _ = chain_iterative(MB, True)
        d = [b - a for a, b in zip(h0b, hbb)]
        # Front-loaded at k (delta_k > delta_{r-k}, k <= r/2) is exactly a
        # violation of Conjecture C read on the REVERSED sequence.
        rev = conjecture_C_violations(d[::-1], MB.r)
        if rev:
            excluders += 1
        rows.append((label, d, rev))
    return rows, members, excluders, skipped


def section_controls():
    info("--- controls: the same predicate applied to other matroids ---")
    U = Matroid(A_CTRL_UNIFORM)
    h0u, _ = chain_iterative(U, False)
    hbu, _ = chain_iterative(U, True)
    du = [b - a for a, b in zip(h0u, hbu)]
    info("control 1 = U(5,9) Vandermonde: W =", U.W, " h^0 =", h0u,
         " h^beta =", hbu, " delta =", du)
    check("control1_is_simple_rank5_on_9_elements",
          U.r == 5 and U.n == 9 and is_simple(U))
    check("control1_is_the_uniform_matroid_U_5_9",
          U.W == [binom(9, k) for k in range(5)] + [1], str(U.W))
    check("control1_conjecture_C_holds",
          conjecture_C_violations(du, U.r) == [], "delta = %s" % du)

    F = Matroid(A_CTRL_FRONT)
    h0f, _ = chain_iterative(F, False)
    hbf, _ = chain_iterative(F, True)
    df = [b - a for a, b in zip(h0f, hbf)]
    info("control 2: W =", F.W, " h^0 =", h0f, " h^beta =", hbf,
         " delta =", df)
    check("control2_is_simple_rank5_on_9_elements",
          F.r == 5 and F.n == 9 and is_simple(F))
    check("control2_has_a_nonzero_differential_defect", any(d != 0 for d in df),
          "delta = %s" % df)
    check("control2_defect_is_strictly_front_loaded_at_k2", df[2] > df[3],
          "delta_2=%d > delta_3=%d" % (df[2], df[3]))
    check("control2_conjecture_C_holds",
          conjecture_C_violations(df, F.r) == [], "delta = %s" % df)
    # The one computable piece of evidence that the paper's index orientation
    # (h_k = dim C_k) is the one in which Conjecture C was posed: under the
    # reversed (degree) indexing this randomly found matroid ALSO violates the
    # conjecture, so the reversed reading would make the conjecture false for
    # generic simple rank-5 matroids on 9 elements -- not a statement anyone
    # would conjecture, and not a refutation anyone would need 9 elements for.
    # If this line ever fails, both orientations are tenable and the paper's
    # refutation is ambiguous rather than established.
    dfrev = df[::-1]
    info("control 2 under the REVERSED index orientation: delta =", dfrev,
         "violating indices", conjecture_C_violations(dfrev, F.r))
    check("reversed_index_orientation_is_excluded_by_control2",
          conjecture_C_violations(dfrev, F.r) != [],
          "reversed delta = %s violates at %s, so the reversed reading is not "
          "the conjecture's" % (dfrev, conjecture_C_violations(dfrev, F.r)))

    # WIDENING the same argument off a single matrix.  The reversed reading
    # asserts that no simple matroid is strictly front-loaded, so every
    # strictly front-loaded member of the deterministic scan below falsifies it
    # on its own.  The scan is written out in full (ORIENT_SCAN_CELLS x
    # ORIENT_SCAN_VALUES around A_CTRL_FRONT): no random source, no seed, no
    # external search code, and enlarging it means widening one of two lists.
    rows, members, excluders, skipped = orientation_scan(A_CTRL_FRONT)
    for label, d, rev in rows:
        if d is None:
            info("orientation scan %-12s: not simple of rank 5 on 9 elements,"
                 " skipped (counted in the tally below)" % label)
        else:
            info("orientation scan %-12s: delta = %s  strictly front-loaded at"
                 " %s" % (label, d, rev if rev else "no index"))
    info("orientation scan: %d matrices scanned, %d simple rank-5 members on 9"
         " elements, %d skipped, %d members strictly front-loaded"
         % (len(rows), members, skipped, excluders))
    check("orientation_scan_reproduces_control2_delta",
          bool(rows) and rows[0][1] == df,
          "seed member delta = %s vs control 2 delta = %s"
          % (rows[0][1] if rows else None, df))
    check("reversed_index_orientation_excluded_by_a_regenerable_scan",
          members >= 1 and excluders >= 1,
          "%d of %d scanned members are strictly front-loaded; the reversed "
          "reading asserts that none is, so it is not the reading in which "
          "Conjecture C was posed" % (excluders, members))


# --------------------------------------------------------------------------
# ACCOUNTING.  The printed total counts every PASS, self-tests included, so it
# is not a count of evidence about M.  The classification below is the one
# argued in the module docstring; the four numbers reported from it are COUNTED
# from the check names that actually ran.
# --------------------------------------------------------------------------

ACCOUNTING_CHECK = "declared_non_evidence_check_names_all_ran"

# Theorems about any matrix / any matroid, retained as implementation
# self-tests (each fails under an injected bug in the routine it exercises).
NON_EVIDENCE_SELF_TESTS = frozenset([
    "rank_unit_increase_and_bounded",
    "rank_submodular_on_all_subset_pairs",
    "closure_is_idempotent",
    "flats_closed_under_intersection",
    "flat_lattice_is_graded",
    "delta_is_nonnegative",
    "C_k^0_is_a_subspace_of_C_k^beta",
    "word_generators_and_iterative_chain_agree",
    "the_two_forms_of_conjecture_C_agree_at_every_k",
    ACCOUNTING_CHECK,
])

# Lines that only restate the input matrix, plus the double-entry transcription
# guard, which compares one hand transcription of the printed matrix against a
# second one -- input against input, not a mathematical check.
NON_EVIDENCE_INPUT_RESTATEMENTS = frozenset([
    "matrix_shape_5x9",
    "matrix_entries_are_integers",
    "first_five_columns_are_identity",
    "ground_set_has_9_elements",
    "last_four_columns_match_paper",
])

# Checks about the control matroids and the orientation scan: real facts, but
# about other objects, bearing on input I5 and on the two-sidedness of the
# predicate rather than on M.
CONTROL_PREFIXES = ("control1_", "control2_", "reversed_", "orientation_")


def section_accounting():
    """Split the printed total four ways, from the names that actually ran."""
    declared = NON_EVIDENCE_SELF_TESTS | NON_EVIDENCE_INPUT_RESTATEMENTS
    # ACCOUNTING_CHECK names the line being emitted right now, so its own name
    # is not yet in _NAMES while this condition is being evaluated.
    missing = sorted(declared - set(_NAMES) - set([ACCOUNTING_CHECK]))
    check(ACCOUNTING_CHECK, not missing,
          "%d names declared non-evidence, %d of them never ran %s"
          % (len(declared), len(missing), missing))
    selft = sum(1 for nm in _NAMES if nm in NON_EVIDENCE_SELF_TESTS)
    restate = sum(1 for nm in _NAMES if nm in NON_EVIDENCE_INPUT_RESTATEMENTS)
    ctrl = sum(1 for nm in _NAMES
               if nm.startswith(CONTROL_PREFIXES)
               and nm not in NON_EVIDENCE_SELF_TESTS
               and nm not in NON_EVIDENCE_INPUT_RESTATEMENTS)
    evid = _N - selft - restate - ctrl
    info("ACCOUNTING: the total below counts EVERY printed PASS, so it is not")
    info("ACCOUNTING: a count of evidence about M.  Of %d checks: %d are"
         % (_N, evid))
    info("ACCOUNTING: evidence about M itself; %d are about the two control"
         % ctrl)
    info("ACCOUNTING: matroids and the orientation scan, which bear on input I5")
    info("ACCOUNTING: and on the two-sidedness of the predicate, not on M; %d"
         % selft)
    info("ACCOUNTING: are implementation self-tests, i.e. theorems true of any")
    info("ACCOUNTING: matrix, retained because each fails under an injected bug")
    info("ACCOUNTING: in the routine it exercises; and %d only restate the"
         % restate)
    info("ACCOUNTING: input matrix -- its shape, its integrality, the identity")
    info("ACCOUNTING: block, n=9, and last_four_columns_match_paper, a")
    info("ACCOUNTING: double-entry guard comparing one transcription against")
    info("ACCOUNTING: another.  All four numbers are counted from the names")
    info("ACCOUNTING: that ran.  The classification is conservative: a referee")
    info("ACCOUNTING: may judge further lines (bottom_flat_is_empty,")
    info("ACCOUNTING: ground_set_is_the_unique_top_flat) to be theorems about")
    info("ACCOUNTING: any matroid as well.")


def main():
    t0 = _now()
    info("verifier for: A Nine-Element Counterexample to Sinclair's")
    info("              Front-Loading Conjecture (Theorem 1, Remark 2)")
    info("exact arithmetic only; python", sys.version.split()[0])
    M = Matroid(A_PAPER)
    section_object(M)
    section_rank_oracle(M)
    section_lattice(M)
    h0, hb, g0, gb = section_hilbert(M)
    section_certificates(M, h0, hb, g0, gb)
    section_refutation(M, h0, hb)
    section_remark(hb)
    section_controls()
    section_accounting()
    info("SCOPE: the paper's full claim is re-derived from the matrix A alone.")
    info("SCOPE: the paper makes no minimality or census claim (Remark 2),")
    info("SCOPE: so there is no exhaustive search to reproduce; the two")
    info("SCOPE: control matroids are extra, not part of the paper's claim.")
    info("SCOPE: two values are TAKEN ON TRUST, both flagged above: the")
    info("SCOPE: transcription of Sinclair's operator formula eq. (2) (I2),")
    info("SCOPE: and the index orientation h_k = dim C_k (I5).  The verdict")
    info("SCOPE: below is conditional on both; the orientation reverses delta")
    info("SCOPE: and hence the verdict, and what argues for the paper's reading")
    info("SCOPE: is control 2 together with the deterministic orientation scan.")
    info("NOT RE-RUN: the transcription of Sinclair's operator formula eq. (2)")
    info("NOT RE-RUN: (input I2).  It is a reading of the source rather than a")
    info("NOT RE-RUN: computation, and nothing in this program can check it.")
    info("NOT RE-RUN: the index orientation h_k = dim C_k (input I5).  The")
    info("NOT RE-RUN: verdict is orientation-dependent: under the reversed")
    info("NOT RE-RUN: (degree) reading, delta for M is the reverse of the")
    info("NOT RE-RUN: sequence printed above and Conjecture C HOLDS for M --")
    info("NOT RE-RUN: both values are printed on the ORIENTATION lines.  What")
    info("NOT RE-RUN: is computed here is EVIDENCE against that reading (it")
    info("NOT RE-RUN: would already be false for control 2 and for every")
    info("NOT RE-RUN: strictly front-loaded member of the orientation scan),")
    info("NOT RE-RUN: not a derivation of the orientation, which must still be")
    info("NOT RE-RUN: confirmed against Sinclair's preprint.")
    info("NOT RE-RUN: nothing else.  The paper makes no census and no")
    info("NOT RE-RUN: minimality claim (Remark 2), so no exhaustive search has")
    info("NOT RE-RUN: to be reproduced, and every number of Theorem 1 and")
    info("NOT RE-RUN: Remark 2 is re-derived above from the matrix A alone.")
    info("elapsed %.1f s" % (_now() - t0))
    if _FAILED:
        info("failed:", ", ".join(_FAILED))
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(_FAILED), _N))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % _N)
    return 0


def _now():
    import time
    return time.time()


if __name__ == "__main__":
    sys.exit(main())


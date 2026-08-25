#!/usr/bin/env python3
"""Verification of an order-14 counterexample to Conjecture 4 of Pioge et al.

The statement acted on: for A in H_n (positive semidefinite Hermitian) put
(F_A)_{ij} = a_{ij} per(A(i,j)), where A(i,j) deletes row i and column j.
Conjecture 4 asserts lambda_max(Re F_A) <= per(A).  The paper exhibits an
order-14 A for which this fails; the published counterexample had order 16.

TAKEN FROM THE PAPER (inputs, transcribed and nothing else):
  * the 2 x 14 Gaussian-integer matrix X = X_RE + i X_IM whose columns
    generate the rank-two Gram matrix A = X^* X;
  * the integer witness vector v;
  * the numbers the paper asserts, kept apart from the computation so that
    they are compared against, never used by it: det of the first two
    columns of X, v^T v, the 90-digit gap Delta, the printed decimal digits
    of the Rayleigh ratio and of the bunching coefficient, the diagonal
    weights s, a, b, c and the scalar eta.

DERIVED HERE (computed, never asserted):
  * A = X^* X, then Hermiticity, positive semidefiniteness and rank exactly
    two from principal minors in exact integer arithmetic;
  * per(A) and all 196 cofactor permanents per(A(i,j)) by the rank-two
    expansion, each one re-run independently by Ryser's formula, which is
    itself first calibrated on closed-form permanents;
  * F_A, its Hermiticity, and the Laplace identity that every row and column
    of F_A sums to per(A);
  * the load-bearing inequality: v^T (Re F_A) v - per(A) v^T v > 0, hence
    lambda_max(Re F_A) >= v^T (Re F_A) v / v^T v > per(A);
  * calibration that the test is not vacuous: twelve unrelated rank-two
    Gaussian Gram matrices show no violation at the same v, and the all-ones
    Gram matrix attains equality;
  * an explicit entrywise positive correlation matrix B with
    per(A o B) > per(A) in exact rational arithmetic, refuting Conjecture 3
    in order 14 rather than only its second-order expansion;
  * the diagonal-congruence lemma, XWX^* = eta I_2, AWA = eta A (so the
    rescaled matrix is a rank-two orthogonal projection of trace 2), and the
    order-15 padding, whose permanent and whole F matrix are recomputed from
    scratch by Ryser's formula.

NOT RE-RUN (the closing line printed by the program states the same list):
the original authors' unsuccessful numerical search in order 15; orders at most
13, which the paper does not claim; the "every order n >= 14" corollary beyond
k = 1, only the order-15 instance of the block identity being recomputed; the
explicit unitary completion U in C^{14x14} of the two orthonormal rows of
Y = eta^{-1/2} X W^{1/2}, which is not exactly representable here; and the
transcendental Gaussian time-delay matrix S(d) itself.
"""

import sys
from fractions import Fraction

N = 14

# ---- values TAKEN FROM THE PAPER -------------------------------------------
X_RE = [
    [1178, -223, 145, 278, -388, -78, -9, 163, -347, 602, -344, 529, -168, 77],
    [-30, -50, 68, -14, 17, 369, 866, -686, -467, -88, -133, 1199, -794, 247],
]
X_IM = [
    [317, -59, 327, 310, 426, -109, 942, 521, -1308, 672, -52, 117, -497, -368],
    [838, 744, -92, -307, 438, 270, -403, -274, -107, 907, -728, 386, -41, 166],
]
V = [-24, 3, 24, 17, 14, -5, 13, -16, 22, -17, -12, -4, 8, -23]

CLAIM_DET_X12 = (-350880, 1045686)
CLAIM_VTV = 3622
CLAIM_DELTA = (188558958243265516176136227271677528258899031 * 10 ** 45
               + 880644911529452618269391965199062697646489600)
CLAIM_RATIO_DIGITS = "1.001531003891338"
CLAIM_BUNCH_DIGITS = "0.003062007782676"
CLAIM_S = 5975547235929
CLAIM_A = 91619936424609
CLAIM_B = 170572929732992
CLAIM_C = 119336587939769
CLAIM_ETA = 144241568073276041742
# ---------------------------------------------------------------------------

CHECKS = []


def gadd(z, w):
    return (z[0] + w[0], z[1] + w[1])


def gsub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gscale(z, c):
    return (z[0] * c, z[1] * c)


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return bool(ok)


def finish():
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


def build_X():
    """Columns of X as pairs of Gaussian integers (x_{1j}, x_{2j})."""
    cols = []
    for j in range(N):
        cols.append(((X_RE[0][j], X_IM[0][j]), (X_RE[1][j], X_IM[1][j])))
    return cols


def build_A(cols):
    """A = X^* X, so a_{ij} = <col_i, col_j> = sum_r conj(x_{ri}) x_{rj}."""
    A = []
    for i in range(N):
        row = []
        for j in range(N):
            s = (0, 0)
            for r in range(2):
                s = gadd(s, gmul(gconj(cols[i][r]), cols[j][r]))
            row.append(s)
        A.append(row)
    return A


def check_shapes():
    """Runs before anything indexes the transcribed data."""
    ok = (len(X_RE) == 2 and len(X_IM) == 2
          and all(len(r) == N for r in X_RE + X_IM) and len(V) == N)
    ok = ok and all(isinstance(e, int) for r in X_RE + X_IM for e in r)
    ok = ok and all(isinstance(e, int) for e in V)
    return ck("X_shape_2_by_14_gaussian_integers", ok,
              "rows %s, witness length %d"
              % ([len(r) for r in X_RE + X_IM], len(V)))


def check_wellformed(cols, A):
    nz = sum(1 for c in cols if c != ((0, 0), (0, 0)))
    ck("X_has_14_nonzero_columns", nz == N,
       "%d of %d columns nonzero, so the degree-14 product is nonzero"
       % (nz, N))
    ck("A_is_14_by_14", len(A) == N and all(len(r) == N for r in A),
       "order %d, %d entries" % (len(A), sum(len(r) for r in A)))
    print("  A[0][0..3] = " + " ".join(str(A[0][j]) for j in range(4)))
    print("  diag(A)    = " + " ".join(str(A[i][i][0]) for i in range(N)))


def check_hermitian(A):
    bad = [(i, j) for i in range(N) for j in range(N)
           if A[i][j] != gconj(A[j][i])]
    ck("A_is_hermitian", not bad, "%d violations" % len(bad))
    ck("A_diagonal_real_positive",
       all(A[i][i][1] == 0 and A[i][i][0] > 0 for i in range(N)),
       "min diag %d" % min(A[i][i][0] for i in range(N)))


def gdet3(M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, k = M[2]
    t1 = gmul(a, gsub(gmul(e, k), gmul(f, h)))
    t2 = gmul(b, gsub(gmul(d, k), gmul(f, g)))
    t3 = gmul(c, gsub(gmul(d, h), gmul(e, g)))
    return gadd(gsub(t1, t2), t3)


def check_psd_and_rank(cols, A):
    """PSD: for a Hermitian matrix, PSD <=> every principal minor >= 0.
    Rank 2 makes all principal minors of size >= 3 vanish, so sizes 1 and 2
    decide positive semidefiniteness."""
    neg1 = [i for i in range(N) if A[i][i][0] < 0]
    neg2 = []
    for i in range(N):
        for j in range(i + 1, N):
            m = A[i][i][0] * A[j][j][0] - (A[i][j][0] ** 2 + A[i][j][1] ** 2)
            if m < 0:
                neg2.append((i, j, m))
    bad3 = []
    for i in range(N):
        for j in range(i + 1, N):
            for k in range(j + 1, N):
                idx = (i, j, k)
                sub = [[A[r][c] for c in idx] for r in idx]
                if gdet3(sub) != (0, 0):
                    bad3.append(idx)
    ck("A_all_3x3_principal_minors_vanish", not bad3,
       "%d of 364 nonzero" % len(bad3))
    big2 = max(A[i][i][0] * A[j][j][0] - (A[i][j][0] ** 2 + A[i][j][1] ** 2)
               for i in range(N) for j in range(i + 1, N))
    ck("A_rank_exactly_two", not bad3 and big2 > 0,
       "largest 2x2 principal minor = %d > 0" % big2)
    ck("A_is_positive_semidefinite", not neg1 and not neg2 and not bad3,
       "%d negative 1x1, %d negative 2x2 principal minors"
       % (len(neg1), len(neg2)))
    d = gsub(gmul(cols[0][0], cols[1][1]), gmul(cols[1][0], cols[0][1]))
    ck("det_of_first_two_columns_matches_paper", d == CLAIM_DET_X12,
       "computed %d%+di" % d)


FACT = [1]
for _t in range(1, 40):
    FACT.append(FACT[-1] * _t)


def poly_of(cols, idx):
    """Coefficients q_k of prod_{j in idx} (x_{1j} t + x_{2j})."""
    q = [(1, 0)]
    for j in idx:
        a, b = cols[j][0], cols[j][1]
        new = [(0, 0)] * (len(q) + 1)
        for k, c in enumerate(q):
            new[k] = gadd(new[k], gmul(c, b))
            new[k + 1] = gadd(new[k + 1], gmul(c, a))
        q = new
    return q


def per_rank_two(qI, qJ):
    """Lemma: per(X_I^* X_J) = sum_k k!(m-k)! conj(q_{I,k}) q_{J,k}."""
    m = len(qI) - 1
    assert len(qJ) - 1 == m
    s = (0, 0)
    for k in range(m + 1):
        s = gadd(s, gscale(gmul(gconj(qI[k]), qJ[k]), FACT[k] * FACT[m - k]))
    return s


def ryser(M):
    """Independent permanent by Ryser's formula with Gray-code subsets."""
    n = len(M)
    if n == 0:
        return (1, 0)
    rowsum = [(0, 0)] * n
    total = (0, 0)
    prev = 0
    for g in range(1, 1 << n):
        gray = g ^ (g >> 1)
        diff = gray ^ prev
        col = diff.bit_length() - 1
        if gray & diff:
            for i in range(n):
                rowsum[i] = gadd(rowsum[i], M[i][col])
        else:
            for i in range(n):
                rowsum[i] = gsub(rowsum[i], M[i][col])
        prev = gray
        prod = (1, 0)
        for i in range(n):
            prod = gmul(prod, rowsum[i])
        if bin(gray).count("1") & 1:
            total = gsub(total, prod)
        else:
            total = gadd(total, prod)
    return total if n % 2 == 0 else (-total[0], -total[1])


def gram(cols, I, J):
    out = []
    for i in I:
        row = []
        for j in J:
            s = (0, 0)
            for r in range(2):
                s = gadd(s, gmul(gconj(cols[i][r]), cols[j][r]))
            row.append(s)
        out.append(row)
    return out


def check_lemma_identity(cols):
    """The rank-two permanent formula is validated against Ryser's formula on
    genuine submatrices X_I^* X_J of the exhibited X."""
    cases = [((0, 1, 2), (3, 4, 5)), ((0, 2, 4, 6), (1, 3, 5, 7)),
             ((0, 1, 2, 3, 4), (0, 1, 2, 3, 4)),
             ((2, 5, 7, 9, 11, 13), (0, 1, 4, 6, 8, 12)),
             ((1, 3, 5, 7, 9, 11, 13), (1, 3, 5, 7, 9, 11, 13))]
    agree = 0
    for I, J in cases:
        lhs = ryser(gram(cols, I, J))
        rhs = per_rank_two(poly_of(cols, I), poly_of(cols, J))
        if lhs == rhs:
            agree += 1
        else:
            print("  mismatch on |I|=%d: %s vs %s" % (len(I), lhs, rhs))
    ck("rank_two_permanent_formula_matches_ryser", agree == len(cases),
       "%d of %d submatrix pairs agree" % (agree, len(cases)))


def permanent_data(cols):
    full = poly_of(cols, range(N))
    p = per_rank_two(full, full)
    leave = [poly_of(cols, [j for j in range(N) if j != i]) for i in range(N)]
    cof = [[per_rank_two(leave[i], leave[j]) for j in range(N)]
           for i in range(N)]
    return p, cof


def check_permanent(cols, A, p, cof):
    ck("permanent_of_A_is_a_positive_integer", p[1] == 0 and p[0] > 0,
       "per(A) = %d (%d digits), Im = %d" % (p[0], len(str(abs(p[0]))), p[1]))
    r = ryser(A)
    ck("permanent_of_A_agrees_with_ryser", r == p, "ryser gives %s" % (r,))
    pairs = [(i, j) for i in range(N) for j in range(N)]
    agree = 0
    for i, j in pairs:
        I = [r_ for r_ in range(N) if r_ != i]
        J = [c_ for c_ in range(N) if c_ != j]
        if ryser(gram(cols, I, J)) == cof[i][j]:
            agree += 1
    ck("sampled_cofactor_permanents_agree_with_ryser", agree == len(pairs),
       "%d of %d order-13 cofactors re-run by ryser" % (agree, len(pairs)))
    hb = [(i, j) for i in range(N) for j in range(N)
          if cof[i][j] != gconj(cof[j][i])]
    ck("cofactor_permanents_are_hermitian_conjugates", not hb,
       "%d violations" % len(hb))


def build_F(A, cof):
    return [[gmul(A[i][j], cof[i][j]) for j in range(N)] for i in range(N)]


def check_F(F, p):
    bad = [(i, j) for i in range(N) for j in range(N)
           if F[i][j] != gconj(F[j][i])]
    ck("F_A_is_hermitian_so_Re_F_A_is_symmetric", not bad,
       "%d violations" % len(bad))
    rs, cs = [], []
    for i in range(N):
        s = (0, 0)
        for j in range(N):
            s = gadd(s, F[i][j])
        rs.append(s)
    for j in range(N):
        s = (0, 0)
        for i in range(N):
            s = gadd(s, F[i][j])
        cs.append(s)
    ck("every_row_and_column_sum_of_F_A_equals_per_A",
       all(s == p for s in rs) and all(s == p for s in cs),
       "%d rows and %d columns equal per(A)"
       % (sum(1 for s in rs if s == p), sum(1 for s in cs if s == p)))


def trunc_digits(fr, ndec):
    """Decimal truncation of a positive Fraction to ndec places, as a string."""
    scaled = fr.numerator * 10 ** ndec // fr.denominator
    s = str(scaled).rjust(ndec + 1, "0")
    return s[:-ndec] + "." + s[-ndec:]


def check_witness_vector():
    ck("witness_vector_sums_to_zero_and_has_norm_3622",
       sum(V) == 0 and sum(x * x for x in V) == CLAIM_VTV,
       "sum = %d, v^T v = %d" % (sum(V), sum(x * x for x in V)))


def check_violation(F, p):
    quad = 0
    for i in range(N):
        for j in range(N):
            quad += V[i] * V[j] * F[i][j][0]
    imag = sum(V[i] * V[j] * F[i][j][1] for i in range(N) for j in range(N))
    vtv = sum(x * x for x in V)
    delta = quad - p[0] * vtv
    ck("quadratic_form_of_imaginary_part_vanishes", imag == 0,
       "v^T Im(F_A) v = %d" % imag)
    ck("delta_matches_the_published_integer", delta == CLAIM_DELTA,
       "computed %d digits, %s"
       % (len(str(abs(delta))), "equal" if delta == CLAIM_DELTA
          else "DIFFERS by %d" % (delta - CLAIM_DELTA)))
    ck("conclusion_of_conjecture_4_is_violated", delta > 0,
       "lambda_max(Re F_A) >= v^T(Re F_A)v / v^T v > per(A) since "
       "v^T(Re F_A)v - per(A)*v^T v = %d > 0" % delta)
    ratio = Fraction(quad, p[0] * vtv)
    got = trunc_digits(ratio, len(CLAIM_RATIO_DIGITS.split(".")[1]))
    ck("rayleigh_ratio_digits_match_paper", got == CLAIM_RATIO_DIGITS,
       "computed %s, paper %s" % (got, CLAIM_RATIO_DIGITS))
    bunch = 2 * (ratio - 1)
    gotb = trunc_digits(bunch, len(CLAIM_BUNCH_DIGITS.split(".")[1]))
    ck("bunching_coefficient_digits_match_paper", gotb == CLAIM_BUNCH_DIGITS,
       "computed %s, paper %s" % (gotb, CLAIM_BUNCH_DIGITS))
    ck("conjecture_3_quadratic_coefficient_is_positive",
       Fraction(quad, vtv) - p[0] > 0,
       "tau^T(Re F_A)tau - per(A) = %s > 0"
       % trunc_digits(Fraction(delta, vtv), 3))
    return delta, ratio - 1


def minor(M, i, j):
    return [[M[r][c] for c in range(len(M)) if c != j]
            for r in range(len(M)) if r != i]


def F_by_ryser(M):
    n = len(M)
    return [[gmul(M[i][j], ryser(minor(M, i, j))) for j in range(n)]
            for i in range(n)]


def check_congruence_lemma(cols):
    """per(DAD) = rho per(A) and F_{DAD} = rho F_A with rho = prod d_i^2,
    tested on a genuine rank-two Gram submatrix, not asserted."""
    I = (0, 3, 6, 9, 12)
    B = gram(cols, I, I)
    d = [2, 1, 3, 1, 5]
    rho = 1
    for x in d:
        rho *= x * x
    C = [[gscale(B[i][j], d[i] * d[j]) for j in range(5)] for i in range(5)]
    okp = ryser(C) == gscale(ryser(B), rho)
    FB, FC = F_by_ryser(B), F_by_ryser(C)
    okF = all(FC[i][j] == gscale(FB[i][j], rho)
              for i in range(5) for j in range(5))
    ck("diagonal_congruence_lemma_holds_on_a_rank_two_gram_block", okp and okF,
       "per scaling %s, F scaling %s" % (okp, okF))


def check_projection(cols, A, p):
    W = [CLAIM_S] * N
    W[3], W[4], W[7] = CLAIM_A, CLAIM_B, CLAIM_C
    ck("scaling_weights_are_positive_integers",
       all(isinstance(w, int) and w > 0 for w in W), "min weight %d" % min(W))
    M = [[(0, 0), (0, 0)], [(0, 0), (0, 0)]]
    for r in range(2):
        for s in range(2):
            acc = (0, 0)
            for j in range(N):
                acc = gadd(acc, gscale(gmul(cols[j][r], gconj(cols[j][s])),
                                       W[j]))
            M[r][s] = acc
    want = [[(CLAIM_ETA, 0), (0, 0)], [(0, 0), (CLAIM_ETA, 0)]]
    ck("X_W_Xstar_equals_eta_times_identity", M == want,
       "computed [[%s,%s],[%s,%s]]" % (M[0][0], M[0][1], M[1][0], M[1][1]))
    awa_ok = True
    for i in range(N):
        for j in range(N):
            acc = (0, 0)
            for k in range(N):
                acc = gadd(acc, gscale(gmul(A[i][k], A[k][j]), W[k]))
            if acc != gscale(A[i][j], CLAIM_ETA):
                awa_ok = False
    ck("P_is_idempotent_via_exact_identity_AWA_equals_eta_A", awa_ok,
       "A W A = eta A entrywise: %s" % awa_ok)
    tr = sum(W[i] * A[i][i][0] for i in range(N))
    ck("projection_has_trace_two_and_diagonal_in_unit_interval",
       tr == 2 * CLAIM_ETA
       and all(0 < W[i] * A[i][i][0] <= CLAIM_ETA for i in range(N)),
       "trace = %s, max diagonal entry = %s"
       % (Fraction(tr, CLAIM_ETA),
          trunc_digits(Fraction(max(W[i] * A[i][i][0] for i in range(N)),
                                CLAIM_ETA), 6)))


def check_order_fifteen(A, F, p, delta):
    """Corollary: A_1 = A (+) I_1 is a PSD order-15 witness. Its permanent and
    its whole F matrix are recomputed from scratch by Ryser's formula."""
    n = N + 1
    A1 = [[(0, 0)] * n for _ in range(n)]
    for i in range(N):
        for j in range(N):
            A1[i][j] = A[i][j]
    A1[N][N] = (1, 0)
    p1 = ryser(A1)
    ck("padded_order_15_matrix_has_the_same_permanent", p1 == p,
       "per(A (+) I_1) = per(A): %s" % (p1 == p))
    F1 = F_by_ryser(A1)
    ok = True
    for i in range(n):
        for j in range(n):
            exp = (p if (i == N and j == N) else
                   ((0, 0) if (i == N or j == N) else F[i][j]))
            if F1[i][j] != exp:
                ok = False
    ck("padded_F_matrix_is_the_direct_sum_F_A_plus_per_A", ok,
       "F_{A_1} = F_A (+) per(A) I_1: %s" % ok)
    v1 = V + [0]
    quad1 = sum(v1[i] * v1[j] * F1[i][j][0] for i in range(n)
                for j in range(n))
    d1 = quad1 - p1[0] * sum(x * x for x in v1)
    ck("conjecture_4_also_fails_in_order_15", d1 > 0 and d1 == delta,
       "order-15 gap equals the order-14 gap and is positive: %s" % (d1 > 0))


def check_ryser_sanity():
    """The independent oracle is itself validated on closed-form permanents."""
    ok = True
    for n in range(1, 8):
        J = [[(1, 0)] * n for _ in range(n)]
        I = [[(1, 0) if i == j else (0, 0) for j in range(n)]
             for i in range(n)]
        if ryser(J) != (FACT[n], 0) or ryser(I) != (1, 0):
            ok = False
    small = [[(1, 0), (2, 0)], [(3, 0), (4, 0)]]
    ok = ok and ryser(small) == (10, 0)
    tri = [[(0, 1), (1, 0), (0, 0)], [(1, 0), (0, 0), (1, 0)],
           [(0, 0), (1, 0), (0, 1)]]
    ok = ok and ryser(tri) == (0, 2)
    ck("ryser_oracle_reproduces_closed_form_permanents", ok,
       "per(J_n)=n!, per(I_n)=1 for n<=7, plus two hand-checked cases")


def check_A_matches_lemma(cols, A):
    """The order-1 instance of the rank-two formula must return A itself."""
    bad = 0
    for i in range(N):
        for j in range(N):
            if per_rank_two(poly_of(cols, [i]), poly_of(cols, [j])) != A[i][j]:
                bad += 1
    ck("gram_entries_agree_with_the_rank_two_formula_at_m_equals_one",
       bad == 0, "%d of %d entries disagree" % (bad, N * N))


def _delta_for(cols, vec=None):
    if vec is None:
        vec = V
    A = build_A(cols)
    full = poly_of(cols, range(N))
    p = per_rank_two(full, full)
    leave = [poly_of(cols, [j for j in range(N) if j != i]) for i in range(N)]
    quad = 0
    for i in range(N):
        for j in range(N):
            quad += vec[i] * vec[j] * gmul(A[i][j],
                                           per_rank_two(leave[i],
                                                        leave[j]))[0]
    return quad - p[0] * sum(x * x for x in vec)


def check_controls():
    """Non-vacuity: the same witness vector applied to unrelated rank-two
    Gaussian Gram matrices must NOT produce a violation, so the positive gap
    is a property of the exhibited matrix and not of the test itself."""
    state = 12345
    viol = []
    for trial in range(12):
        cols = []
        for _ in range(N):
            vals = []
            for _ in range(4):
                state = (1103515245 * state + 12345) % (1 << 31)
                vals.append(-1300 + state % 2601)
            cols.append(((vals[0], vals[1]), (vals[2], vals[3])))
        if _delta_for(cols) > 0:
            viol.append(trial)
    ck("control_rank_two_gram_matrices_show_no_violation", not viol,
       "%d of 12 pseudorandom controls violate the inequality" % len(viol))
    ones = [((1, 0), (0, 0))] * N
    A1 = build_A(ones)
    q = poly_of(ones, range(N))
    eq = (per_rank_two(q, q) == (FACT[N], 0)
          and _delta_for(ones, [1] * N) == 0
          and _delta_for(ones) < 0
          and all(A1[i][j] == (1, 0) for i in range(N) for j in range(N)))
    ck("all_ones_gram_matrix_attains_equality_in_the_conjecture", eq,
       "per(J_14) = 14! and the gap at the all-ones vector is exactly 0, "
       "so the bound is tight and the test is calibrated")


def hadamard_gap(A, p, k):
    """With tau = v/sqrt(3622), eps^2 = 1/k and
        B_ij = (1+eps^2 tau_i tau_j)/sqrt((1+eps^2 tau_i^2)(1+eps^2 tau_j^2)),
    B = D C D with C_ij = K + v_i v_j, K = 3622k and D positive diagonal, so
        per(A o B) = per(C o A)/prod_i (K + v_i^2)
    is an exact rational.  Returns (per(A o B) - p, per(A o B)/p - 1, C)."""
    K = CLAIM_VTV * k
    C = [[K + V[i] * V[j] for j in range(N)] for i in range(N)]
    M = [[gscale(A[i][j], C[i][j]) for j in range(N)] for i in range(N)]
    num = ryser(M)
    den = 1
    for i in range(N):
        den *= K + V[i] * V[i]
    per_ab = Fraction(num[0], den)
    return per_ab - p[0], per_ab / p[0] - 1, C, num[1]


def check_conjecture_three(A, p, ratio_minus_one):
    """Corollary: an explicit entrywise positive correlation matrix B with
    per(A o B) > per(A), certified in exact rational arithmetic."""
    gap, rel, C, im = hadamard_gap(A, p, 1)
    ck("hadamard_witness_C_is_entrywise_positive_and_psd_of_rank_two",
       all(C[i][j] > 0 for i in range(N) for j in range(N))
       and all(C[i][i] * C[j][j] - C[i][j] ** 2 >= 0
               for i in range(N) for j in range(N))
       and all(gdet3([[(C[r][c], 0) for c in idx] for r in idx]) == (0, 0)
               for idx in [(i, j, l) for i in range(N)
                           for j in range(i + 1, N)
                           for l in range(j + 1, N)]),
       "min entry %d, all 3x3 minors vanish" % min(min(r) for r in C))
    ck("permanent_of_hadamard_product_is_real", im == 0, "Im = %d" % im)
    ck("conjecture_3_fails_in_order_14_exactly", gap > 0,
       "per(A o B) - per(A) > 0 with relative excess %s" % trunc_digits(rel, 9))
    r10 = hadamard_gap(A, p, 10)[1] * 10
    r100 = hadamard_gap(A, p, 100)[1] * 100
    ck("hadamard_expansion_coefficient_approached_from_above",
       ratio_minus_one < r100 < r10,
       "c = %s < %s < %s" % (trunc_digits(ratio_minus_one, 12),
                             trunc_digits(r100, 12), trunc_digits(r10, 12)))


def check_time_delay_identity(F, p, delta):
    """The d^2 coefficient of per(P o S(d))/per(P) equals
    -sum_{ij}(tau_i-tau_j)^2 (F_P)_ij / per(P); recomputing it by that route
    must return the published bunching coefficient."""
    s = 0
    for i in range(N):
        for j in range(N):
            s += (V[i] - V[j]) ** 2 * F[i][j][0]
    ck("time_delay_quadratic_coefficient_recomputed_independently",
       s == -2 * delta and s < 0,
       "sum (v_i-v_j)^2 (F_A)_ij = -2*gap, coefficient %s"
       % trunc_digits(Fraction(-s, p[0] * CLAIM_VTV), 15))


if __name__ == "__main__":
    if not check_shapes():
        finish()
    COLS = build_X()
    AA = build_A(COLS)
    check_wellformed(COLS, AA)
    check_hermitian(AA)
    check_psd_and_rank(COLS, AA)
    check_ryser_sanity()
    check_A_matches_lemma(COLS, AA)
    check_lemma_identity(COLS)
    P14, COF = permanent_data(COLS)
    check_permanent(COLS, AA, P14, COF)
    if P14[0] <= 0 or P14[1] != 0:
        print("HALTED: per(A) is not a positive integer, so the ratios that "
              "follow are undefined.")
        finish()
    FA = build_F(AA, COF)
    check_F(FA, P14)
    check_witness_vector()
    DELTA, RM1 = check_violation(FA, P14)
    check_controls()
    check_conjecture_three(AA, P14, RM1)
    check_time_delay_identity(FA, P14, DELTA)
    check_congruence_lemma(COLS)
    check_projection(COLS, AA, P14)
    check_order_fifteen(AA, FA, P14, DELTA)
    print("NOT RE-RUN: the original authors' unsuccessful numerical search "
          "in order 15; orders at most 13, about which the paper claims "
          "nothing; the padding corollary beyond k = 1, i.e. orders 16 and "
          "above are not each recomputed, only the order-15 instance of the "
          "same block identity; the explicit unitary completion "
          "U in C^{14x14} of the two orthonormal rows of "
          "Y = eta^{-1/2} X W^{1/2}, which is not constructed here and is "
          "not exactly representable, eta^{-1/2} being irrational -- what is "
          "checked instead is X W X^* = eta I_2, i.e. that those two rows are "
          "orthonormal, together with A W A = eta A and the trace-2 identity, "
          "so that P = eta^{-1} D A D is an exact Hermitian rank-two "
          "idempotent and the completion follows by Gram-Schmidt; and the "
          "transcendental Gaussian time-delay "
          "matrix S(d) itself, whose d^2 coefficient is instead obtained "
          "exactly from F_A by the two independent routes checked above.")
    finish()

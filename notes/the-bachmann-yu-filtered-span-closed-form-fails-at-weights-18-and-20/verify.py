#!/usr/bin/env python3
# verify.py -- re-derives every computational claim of
#   "The Bachmann-Yu Closed Form for the Filtered Integral Span Fails at Weights 18 and 20".
#
# Python 3.9+, STANDARD LIBRARY ONLY (sys, time, fractions).  No third-party package, no external
# data file, no floating point, no randomness.  The ONLY inputs are the definitions printed in the
# paper:
#
#   x_m       = q^m/(1-q^m)^2 = sum_{j>=1} j q^{mj}                          (paper, Sec. 2)
#   A_r       = e_r(x_1, x_2, ...)                                           (paper, Sec. 2)
#   Lambda_n  = Z-span of the A-monomials A_{mu_1}...A_{mu_s} with |mu| <= n
#             = sum_{|lambda| <= n} Z g(lambda)                              (paper, Lemma 2.1)
#   F_{<=2m}M = Q-span of the A-monomials of degree <= m,
#               dim F_{<=2m}M = #{ (a,b,c) >= 0 : a + 2b + 3c <= m }         (paper, Sec. 2)
#   L_N(k)    = { f in F_{<=k}M : [q^i] f in Z for 0 <= i <= N }             (paper, Sec. 3)
#   i(k,n,N)  = [ L_N(k) : L_N(k) cap Lambda_n ]                             (paper, Sec. 3)
#
# The object the paper exhibits is the pair of finite integer chains
#   i(18, n, 300) = 4782969, 243, 9, 3, 1     (n = 9, ..., 13)
#   i(20, n, 360) = 5811307335, 32805, 27, 3, 1 (n = 10, ..., 14)
# each closed by an ATTAINED 1, and this program recomputes both from scratch, together with all
# fourteen integers of the Bachmann-Yu table for k <= 16 that the paper uses as its control.
#
# All linear algebra is exact: mod-p pivot selection (verified against the predicted rank on every
# call, so an unlucky pivot aborts instead of biasing), Bareiss integer determinants, and a modular
# Hermite normal form.  Indices are determinant RATIOS of integer lattices, hence exact integers,
# and every division is asserted exact before it is taken.

import sys
import time
from fractions import Fraction

T0 = time.time()
_NPASS = 0
_NFAIL = 0


def log(*a):
    print("[%6.1fs]" % (time.time() - T0), *a)
    sys.stdout.flush()


def check(name, got, want):
    """One PASS/FAIL line.  `got == want` is the whole contract; both are printed on failure."""
    global _NPASS, _NFAIL
    if got == want:
        _NPASS += 1
        print("PASS %-30s %s" % (name, got))
    else:
        _NFAIL += 1
        print("FAIL %-30s got %s want %s" % (name, got, want))
    sys.stdout.flush()


# ---------------------------------------------------------------------------
# 1.  q-series over Z, truncated at q^NMAX
# ---------------------------------------------------------------------------
NMAX = 360          # the k = 20 chain is computed at N = 360, the k = 18 chain at N = 300
RMAX = 14           # A_1 ... A_14 ;  Lambda_n for n <= 14


def smul(a, b, prec):
    r = [0] * (prec + 1)
    lb = [(j, bj) for j, bj in enumerate(b[:prec + 1]) if bj]
    for i in range(prec + 1):
        ai = a[i]
        if not ai:
            continue
        lim = prec - i
        for j, bj in lb:
            if j > lim:
                break
            r[i + j] += ai * bj
    return r


# A_r = e_r(x_1, x_2, ...) by the elementary-symmetric recursion
#   e_r <- e_r + x_m * e_{r-1}   (one m at a time),   x_m = sum_{j>=1} j q^{mj}
A = [[0] * (NMAX + 1) for _ in range(RMAX + 1)]
A[0][0] = 1
for m in range(1, NMAX + 1):
    for r in range(RMAX, 0, -1):
        src, tgt = A[r - 1], A[r]
        for j in range(1, NMAX // m + 1):
            p = m * j
            for i in range(p, NMAX + 1):
                s = src[i - p]
                if s:
                    tgt[i] += j * s
log("A_1 ... A_%d built to q^%d" % (RMAX, NMAX))

check("A1-is-sigma1", A[1][1:8], [1, 3, 4, 7, 6, 12, 8])
check("A2-head", [A[2][i] for i in range(3, 9)], [1, 3, 9, 15, 30, 45])

# --- an INDEPENDENT second route to A_r: Newton's identities on the power sums of the x_m.
# p_j = sum_{m>=1} x_m^j ;   r e_r = sum_{j=1}^{r} (-1)^{j-1} e_{r-j} p_j .
# This shares no code with the recursion above and its integrality (r | the Newton sum) is itself a
# check on both.
# x_m(q) = X(q^m) with X(t) = sum_{j>=1} j t^j, so x_m^j = (X^j)(q^m) and
#   p_j = sum_{m>=1} x_m^j  has  [q^n] p_j = sum_{d | n} [t^{n/d}] X^j .
NCTL = 190
X = [0] * (NCTL + 1)
for j in range(1, NCTL + 1):
    X[j] = j
Xp = [None] * (RMAX + 1)
Xp[1] = X
for j in range(2, RMAX + 1):
    Xp[j] = smul(Xp[j - 1], X, NCTL)
psum = [None] * (RMAX + 1)
for j in range(1, RMAX + 1):
    tot = [0] * (NCTL + 1)
    Y = Xp[j]
    for d in range(1, NCTL + 1):
        for n in range(d, NCTL + 1, d):
            c = Y[n // d]
            if c:
                tot[n] += c
    psum[j] = tot
enew = [[0] * (NCTL + 1) for _ in range(RMAX + 1)]
enew[0][0] = 1
newton_exact = True
for r in range(1, RMAX + 1):
    acc = [0] * (NCTL + 1)
    for j in range(1, r + 1):
        sgn = 1 if (j % 2 == 1) else -1
        t = smul(enew[r - j], psum[j], NCTL)
        for i in range(NCTL + 1):
            if t[i]:
                acc[i] += sgn * t[i]
    for i in range(NCTL + 1):
        if acc[i] % r:
            newton_exact = False
        enew[r][i] = acc[i] // r
check("newton-divisibility-exact", newton_exact, True)
check("Ar-vs-newton-r1..%d" % RMAX,
      all(enew[r][:NCTL + 1] == A[r][:NCTL + 1] for r in range(1, RMAX + 1)), True)
log("Newton cross-check of A_1..A_%d to q^%d done" % (RMAX, NCTL))


# ---------------------------------------------------------------------------
# 2.  the A-monomials, i.e. generators of Lambda_n
# ---------------------------------------------------------------------------
def parts(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield ()
        return
    for f in range(min(n, mx), 0, -1):
        for rest in parts(n - f, f):
            yield (f,) + rest


PART = {d: list(parts(d)) for d in range(RMAX + 1)}
MON = {(): A[0]}
for d in range(1, RMAX + 1):
    for p in PART[d]:
        MON[p] = smul(A[p[0]], MON[p[1:]], NMAX)
log("A-monomials of degree <= %d built: %d" % (RMAX, len(MON)))
check("monomials-deg-le-13", sum(len(PART[d]) for d in range(14)), 373)
check("monomials-deg-le-14", sum(len(PART[d]) for d in range(15)), 508)


def gens(n, N):
    """all A-monomials of degree <= n, truncated at q^N -- a generating set of Lambda_n"""
    out = []
    for d in range(n + 1):
        for p in PART[d]:
            out.append(MON[p][:N + 1])
    return out


def dimF(m):
    """dim_Q F_{<=2m}M = #{(a,b,c) >= 0 : a + 2b + 3c <= m}"""
    return sum(1 for c in range(m // 3 + 1)
               for b in range((m - 3 * c) // 2 + 1)
               for a in range(m - 3 * c - 2 * b + 1))


# the same count as partitions with all parts <= 3 and size <= m -- two spellings of one dimension
def dim_partitions(m):
    tot = 0
    for d in range(m + 1):
        tot += sum(1 for p in PART[d] if not p or p[0] <= 3)
    return tot


for m, want in ((9, 53), (10, 67), (12, 102), (13, 123), (14, 147)):
    check("dimF-le-%d" % (2 * m), dimF(m), want)
check("dim-two-spellings-agree", all(dimF(m) == dim_partitions(m) for m in range(RMAX + 1)), True)


# ---------------------------------------------------------------------------
# 3.  exact integer linear algebra
# ---------------------------------------------------------------------------
P0 = (1 << 61) - 1
P1 = (1 << 31) - 1


def indep_rows_mod(rows, ncol, target=None, p=P0, order=None):
    """Greedily pick Z-independent rows (independence certified mod p) and their pivot columns.

    `order` permutes the COLUMN scan, so a second call with a different order yields a different
    injective projection of the same lattice -- used below to confirm that no index depends on the
    choice of coordinates.
    """
    cols = list(range(ncol)) if order is None else order
    M, sel, pivcols = [], [], []
    for idx, rw in enumerate(rows):
        if target is not None and len(sel) >= target:
            break
        v = [x % p for x in rw]
        for pc, pr in zip(pivcols, M):
            if v[pc]:
                f = v[pc]
                v = [(v[t] - f * pr[t]) % p for t in range(ncol)]
        nz = None
        for c in cols:
            if v[c]:
                nz = c
                break
        if nz is None:
            continue
        inv = pow(v[nz], p - 2, p)
        v = [(x * inv) % p for x in v]
        M.append(v)
        pivcols.append(nz)
        sel.append(idx)
    return sel, pivcols


def bareiss(Mat):
    n = len(Mat)
    M = [row[:] for row in Mat]
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    sw = i
                    break
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
            M[i][k] = 0
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def xgcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def lat_hnf(rows, D):
    """Modular HNF (upper triangular, rank D) of the lattice generated by integer rows in Z^D."""
    sel, _ = indep_rows_mod(rows, D, target=D)
    assert len(sel) == D, "generators do not have full rank %d (got %d)" % (D, len(sel))
    d0 = abs(bareiss([rows[i][:] for i in sel]))
    assert d0 != 0
    mod = d0
    H = [[0] * D for _ in range(D)]
    for i in range(D):
        H[i][i] = d0
    for g0 in rows:
        g = [x % mod for x in g0]
        if not any(g):
            continue
        for j in range(D):
            if g[j] == 0:
                continue
            a, b = H[j][j], g[j]
            d, x, y = xgcd(a, b)
            Hj = H[j]
            nH = [(x * Hj[t] + y * g[t]) % mod for t in range(D)]
            ag, bg = a // d, b // d
            ng = [(ag * g[t] - bg * Hj[t]) % mod for t in range(D)]
            if nH[j] == 0:
                nH[j] = mod
            H[j] = nH
            g = ng
        nd = 1
        for j in range(D):
            nd *= H[j][j]
        if nd < mod:
            mod = nd
            for j in range(D):
                H[j] = [x % mod for x in H[j]]
                if H[j][j] == 0:
                    H[j][j] = mod
    return H


def lat_det(rows, D):
    H = lat_hnf(rows, D)
    r = 1
    for j in range(D):
        r *= H[j][j]
    return r


def tri_inverse(H, D):
    Inv = [[Fraction(0)] * D for _ in range(D)]
    for i in range(D - 1, -1, -1):
        Inv[i][i] = Fraction(1, H[i][i])
        for j in range(i + 1, D):
            s = Fraction(0)
            for t in range(i + 1, j + 1):
                s += H[i][t] * Inv[t][j]
            Inv[i][j] = -s / H[i][i]
    return Inv


# ---------------------------------------------------------------------------
# 4.  L_N(k) and the index i(k,n,N)
# ---------------------------------------------------------------------------
def LN_basis(k, N):
    """A Z-basis of L_N(k), as integral q-series, plus the base index [L_N(k) : Lambda_{k/2}].

    Write V = F_{<=k}M and pick dm = dim V independent A-monomials F_1..F_dm of degree <= k/2 as a
    Q-basis of V.  For c in Q^dm the q^i coefficient of sum_t c_t F_t is <c, col_i> with
    col_i = (F_t[i])_t in Z^dm, so L_N(k) is exactly the DUAL lattice of the lattice U spanned by
    col_0..col_N.  A Z-basis of U^* is read off the columns of the inverse of an HNF basis of U.
    """
    m = k // 2
    dm = dimF(m)
    Gm = gens(m, N)
    sel, pivc = indep_rows_mod(Gm, N + 1, target=dm)
    assert len(sel) == dm, "rank %d != dim %d: truncation N=%d too short" % (len(sel), dm, N)
    F = [Gm[i] for i in sel]
    U = [[F[t][i] for t in range(dm)] for i in range(N + 1)]
    HU = lat_hnf(U, dm)
    Inv = tri_inverse(HU, dm)
    LNq = []
    for i in range(dm):
        c = [Inv[t][i] for t in range(dm)]
        ser = [Fraction(0)] * (N + 1)
        for t in range(dm):
            if c[t]:
                Ft = F[t]
                for j in range(N + 1):
                    if Ft[j]:
                        ser[j] += c[t] * Ft[j]
        assert all(v.denominator == 1 for v in ser), "an L_N basis vector is not integral"
        LNq.append([int(v) for v in ser])
    dLam = lat_det([[g[i] for i in pivc] for g in Gm], dm)
    dLN = lat_det([[g[i] for i in pivc] for g in LNq], dm)
    assert dLam % dLN == 0, "Lambda_{k/2} is not contained in L_N"
    return LNq, dLam // dLN


def index_at(k, n, N, LNq, order=None, p=P0):
    """i(k,n,N) = [L_N(k) : L_N(k) cap Lambda_n] = det(Lambda_n) / det(Lambda_n + L_N(k)),
    both determinants taken in a set of Dn = dim F_{<=2n}M coordinates on which Lambda_n injects."""
    Dn = dimF(n)
    Gn = gens(n, N)
    seln, pivn = indep_rows_mod(Gn, N + 1, target=Dn, p=p, order=order)
    assert len(seln) == Dn, "rank %d != dim %d at n=%d, N=%d" % (len(seln), Dn, n, N)
    rowsG = [[g[i] for i in pivn] for g in Gn]
    rowsL = [[g[i] for i in pivn] for g in LNq]
    dG = lat_det(rowsG, Dn)
    dGL = lat_det(rowsG + rowsL, Dn)
    assert dG % dGL == 0, "index is not an integer"
    return dG // dGL


def chain(k, N, nmax):
    """i(k,n,N) for n = k/2 upward, stopping at the first attained 1.  Returns (dict, B, base)."""
    LNq, base = LN_basis(k, N)
    out, B = {}, None
    for n in range(k // 2, nmax + 1):
        ix = index_at(k, n, N, LNq)
        out[n] = ix
        log("  k=%2d n=%2d  i(k,n,N=%d) = %s" % (k, n, N, ix))
        if ix == 1:
            B = n
            break
    return out, B, base, LNq


def factor(x):
    f, d = {}, 2
    while d * d <= x:
        while x % d == 0:
            f[d] = f.get(d, 0) + 1
            x //= d
        d += 1
    if x > 1:
        f[x] = f.get(x, 0) + 1
    return f


def guess(k):
    return k // 2 + k // 6


# ---------------------------------------------------------------------------
# 5.  rank controls
# ---------------------------------------------------------------------------
for n, N, want in ((9, 300, 53), (10, 360, 67), (13, 300, 123)):
    sel, _ = indep_rows_mod(gens(n, N), N + 1, target=want + 1)
    check("rank-Lambda%d-at-N%d" % (n, N), len(sel), want)

# the paper's own remark: 72 coefficients suffice to pin the 67-dimensional space F_{<=20}M ...
sel, _ = indep_rows_mod(gens(10, 72), 73, target=68)
check("rank-F20-at-72-coeffs", len(sel), 67)
# ... and 170 coefficients already pin the 123-dimensional F_{<=26}M, so no index below uses a
# truncation short enough for the projection to lose information.
sel, _ = indep_rows_mod(gens(13, 170), 171, target=124)
check("rank-F26-at-170-coeffs", len(sel), 123)
log("rank controls done")


# ---------------------------------------------------------------------------
# 6.  CONTROL: the fourteen published integers of Bachmann-Yu for k <= 16
# ---------------------------------------------------------------------------
PUB_INDEX = {4: 1, 6: 3, 8: 3, 10: 9, 12: 81, 14: 729, 16: 6561}
PUB_B = {4: 2, 6: 4, 8: 5, 10: 6, 12: 8, 14: 9, 16: 10}
NPUB = 300
pub_chains = {}
for k in (4, 6, 8, 10, 12, 14, 16):
    ch, B, base, LN16 = chain(k, NPUB, 13)
    pub_chains[k] = ch
    check("published-index-k%d" % k, base, PUB_INDEX[k])
    check("published-B-k%d" % k, B, PUB_B[k])
    if k == 16:
        LN_k16 = LN16
check("published-B-matches-guess-k4..16",
      all(PUB_B[k] == guess(k) for k in PUB_B), True)
check("published-indices-are-3-powers",
      sorted(set(sum((list(factor(v)) for v in PUB_INDEX.values() if v > 1), []))), [3])
log("all 14 published integers reproduced")

# FORCED NEGATIVE, inside the very range where the k = 18 growth appears: past its own closure the
# k = 16 chain must stay flat at 1.  An engine that manufactured growth would show it here.
flat16 = {}
for n in (11, 12, 13):
    flat16[n] = index_at(16, n, NPUB, LN_k16)
check("forced-negative-k16-flat", [flat16[n] for n in (11, 12, 13)], [1, 1, 1])

# FORCED NEGATIVE / POSITIVE on a hand-sized certificate at k = 6, where B(6) = 4:
#   h0 = (A_1^3 - A_2)/3 - 3 A_1 A_2 + 7 A_3   lies in Lambda_4 but not in Lambda_3.
A1 = MON[(1,)]
A2 = MON[(2,)]
A3 = MON[(3,)]
A111 = MON[(1, 1, 1)]
A21 = MON[(2, 1)]
h0f = [Fraction(A111[i] - A2[i], 3) - 3 * A21[i] + 7 * A3[i] for i in range(NMAX + 1)]
check("k6-certificate-integral", all(v.denominator == 1 for v in h0f), True)
h0 = [int(v) for v in h0f]


def in_lattice(h, n, N):
    Dn = dimF(n)
    Gn = gens(n, N)
    seln, pivn = indep_rows_mod(Gn, N + 1, target=Dn)
    assert len(seln) == Dn
    rowsG = [[g[i] for i in pivn] for g in Gn]
    hp = [h[i] for i in pivn]
    return lat_det(rowsG, Dn) == lat_det(rowsG + [hp], Dn)


check("k6-certificate-in-Lambda4", in_lattice(h0, 4, NPUB), True)
check("k6-certificate-not-in-Lambda3", in_lattice(h0, 3, NPUB), False)
log("controls done")


# ---------------------------------------------------------------------------
# 7.  THE RESULT: the two chains
# ---------------------------------------------------------------------------
CH18_WANT = {9: 4782969, 10: 243, 11: 9, 12: 3, 13: 1}
CH20_WANT = {10: 5811307335, 11: 32805, 12: 27, 13: 3, 14: 1}

ch18, B18, base18, LN18 = chain(18, 300, 13)
for n in range(9, 14):
    check("index-18-%d" % n, ch18[n], CH18_WANT[n])
check("base-index-18", base18, 4782969)
check("B18", B18, 13)

ch20, B20, base20, LN20 = chain(20, 360, 14)
for n in range(10, 15):
    check("index-20-%d" % n, ch20[n], CH20_WANT[n])
check("base-index-20", base20, 5811307335)
check("B20", B20, 14)

# the closed form, and the size of the failure
check("guess-18", guess(18), 12)
check("guess-20", guess(20), 13)
check("B18-exceeds-guess-by-1", B18 - guess(18), 1)
check("B20-exceeds-guess-by-1", B20 - guess(20), 1)
check("squeeze-closes-18", (ch18[13], ch18[12] > 1), (1, True))
check("squeeze-closes-20", (ch20[14], ch20[13] > 1), (1, True))
# the chain must be a divisibility chain: i(k,n+1) | i(k,n).
check("chain-18-divides",
      all(ch18[n] % ch18[n + 1] == 0 for n in range(9, 13)), True)
check("chain-20-divides",
      all(ch20[n] % ch20[n + 1] == 0 for n in range(10, 14)), True)

# COORDINATE INDEPENDENCE of the two decisive entries per weight: recomputed with the column scan
# rotated by 47 and with a different pivot modulus, i.e. over a different injective projection of the
# same lattice.  An index is a lattice invariant, so it must not move.
def rot(ncol, off=47):
    return list(range(off, ncol)) + list(range(off))


check("index-18-12-other-coords", index_at(18, 12, 300, LN18, order=rot(301), p=P1), 3)
check("index-18-13-other-coords", index_at(18, 13, 300, LN18, order=rot(301), p=P1), 1)
check("index-20-13-other-coords", index_at(20, 13, 360, LN20, order=rot(361), p=P1), 3)
check("index-20-14-other-coords", index_at(20, 14, 360, LN20, order=rot(361), p=P1), 1)

# the factorisations printed in the paper
check("factor-4782969", factor(4782969), {3: 14})
check("factor-5811307335", factor(5811307335), {3: 19, 5: 1})
check("factor-32805", factor(32805), {3: 8, 5: 1})
check("k18-chain-is-3-primary",
      sorted(set(sum((list(factor(v)) for v in ch18.values() if v > 1), []))), [3])
check("k20-chain-not-3-primary",
      sorted(set(sum((list(factor(v)) for v in ch20.values() if v > 1), []))), [3, 5])
check("prime-5-enters-first-at-k20",
      (base18 % 5, base20 % 5 == 0, ch20[11] % 5 == 0), (4, True, True))

# the paper's excess sequence e(n) = B(2n) - n
e_seq = [PUB_B[2 * n] - n for n in range(2, 9)] + [B18 - 9, B20 - 10]
check("excess-sequence", e_seq, [0, 1, 1, 1, 2, 2, 2, 4, 4])
check("excess-skips-3", 3 in e_seq, False)
check("B-values-k4..20", [PUB_B[k] for k in (4, 6, 8, 10, 12, 14, 16)] + [B18, B20],
      [2, 4, 5, 6, 8, 9, 10, 13, 14])


# ---------------------------------------------------------------------------
# 8.  TRUNCATION CONTROL: the base index as a function of N
# ---------------------------------------------------------------------------
NLIST = (72, 80, 100, 120, 150, 200, 250, 300)
for N in NLIST:
    _, b = LN_basis(18, N)
    check("base18-N%d" % N, b, 4782969)
for N in NLIST:
    _, b = LN_basis(20, N)
    check("base20-N%d" % N, b, 87169610025 if N == 72 else 5811307335)
check("factor-87169610025", factor(87169610025), {3: 20, 5: 2})
# The point of the sweep: the authors' own 72 coefficients are enough for the RANK claim (checked
# above) and NOT enough for the integral lattice at k = 20 -- the base index at N = 72 is a proper
# multiple of its stable value, so a finite coefficient check is demonstrably insufficient, not
# merely formally so.  The squeeze, not the truncation length, is what makes the values exact.
check("N72-insufficient-at-k20", 87169610025 % 5811307335 == 0 and 87169610025 != 5811307335, True)


# ---------------------------------------------------------------------------
# 9.  verdict
# ---------------------------------------------------------------------------
print("")
print("NOT RE-RUN: the following are OUTSIDE this program and are not asserted by it.")
print("NOT RE-RUN: (a) B(k) for k = 22 and beyond, and for odd k -- only k <= 20 is computed here;")
print("NOT RE-RUN: (b) any replacement closed form for B(k) -- nine data points are printed, no")
print("NOT RE-RUN:     formula is fitted or tested;")
print("NOT RE-RUN: (c) Bachmann-Yu mainconj:integralspan (R_g = M_Z) in general -- what is verified")
print("NOT RE-RUN:     is only its weight-18 and weight-20 instances, which follow from the attained")
print("NOT RE-RUN:     index 1 at n = 13 and n = 14;")
print("NOT RE-RUN: (d) any explicit element h of F_{<=18}M_Z \\ Lambda_12 -- no witness element is")
print("NOT RE-RUN:     claimed or checked, because 31 or 259 printed coefficients cannot pin a")
print("NOT RE-RUN:     53-dimensional space; the deciding object is the index chain itself;")
print("NOT RE-RUN: (e) the literature: novelty, priority and the Craig near-miss are not machine")
print("NOT RE-RUN:     checkable and nothing here speaks to them;")
print("NOT RE-RUN: (f) the byte offsets quoted from the arXiv source -- this program reads no")
print("NOT RE-RUN:     external file and cannot confirm a line number in someone else's TeX.")
print("")
if _NFAIL:
    print("VERDICT: %d CHECKS FAILED" % _NFAIL)
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _NPASS)
sys.exit(0)

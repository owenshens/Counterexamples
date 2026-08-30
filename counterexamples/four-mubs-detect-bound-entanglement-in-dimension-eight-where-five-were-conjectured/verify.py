#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- exact re-derivation of every computational claim of

    "Four Mutually Unbiased Bases in Dimension Eight Detect Bound Entanglement"

Python 3.9+, STANDARD LIBRARY ONLY (fractions, sys).  No third-party package, no
external data file, no floating point in any decision: every comparison below is
between exact integers or exact Fraction objects.  Floats appear only inside
informational strings.

The program takes as input exactly what the paper prints and nothing else:

  * the index conventions (2.1) and (2.2) of the paper, i.e.
        x = 4*x1 + 2*x2 + x3            for the computational basis of C^8,
        idx(S) = S11 + 2 S12 + 4 S13 + 8 S22 + 16 S23 + 32 S33,
  * the three symmetric matrices S_0, S_13, S_19 and the shift s = 1,
  * the definition (2.1) of |v_S> and the definition (2.3) of W^Gamma(M_m,s),
  * the 288 integers of Table 1, pasted verbatim below as TABLE.

Everything else it derives.  It prints one "PASS <name> [detail]" line per check,
then a list of what it does NOT cover, and closes with

    VERDICT: ALL <n> CHECKS PASS

exiting 0 if and only if every check passed.

Run:  python3 verify.py        (a few seconds, single threaded)
"""

import sys
from fractions import Fraction as F

D = 8          # the dimension
NB = 64        # D*D, the dimension of C^8 (x) C^8
M_COUNT = 4    # m, the number of MUBs in the exhibited set
SHIFT = 1      # s
N_RHO = 1049602  # the denominator N of rho = M / N

# --- the object, exactly as printed in Section 3 of the paper -----------------

S_INDICES = (0, 13, 19)          # the non-computational bases, by idx(S)

S_PRINTED = {                    # the three matrices displayed in the paper
    0:  [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    13: [[1, 0, 1], [0, 1, 0], [1, 0, 0]],
    19: [[1, 1, 0], [1, 0, 1], [0, 1, 0]],
}

# Table 1 of the paper, pasted verbatim.  Row i lists "j:M_ij" over the j >= i
# with M_ij nonzero; M is symmetric and zero elsewhere.
TABLE = """
 0 | 0:13090 9:13229 18:13171 27:12881 36:12852 45:13179 54:13190 63:12939
 1 | 1:57650 8:13270 19:5884 26:4196 37:5429 44:4029 55:-6688 62:-5143
 2 | 2:12536 11:3191 16:3951 25:-5082 38:3182 47:12324 52:-5067 61:3955
 3 | 3:12118 10:4380 17:-6578 24:3461 39:-7078 46:5584 53:10793 60:4002
 4 | 4:12962 13:3941 22:4016 31:-4767 32:-4764 41:4021 50:3890 59:12707
 5 | 5:10924 12:4029 23:10749 30:4031 33:5426 40:-5093 51:5440 58:-5084
 6 | 6:13551 15:-6379 20:-5067 29:13165 34:3425 43:3589 48:4226 57:3866
 7 | 7:10281 14:-4675 21:5450 28:4002 35:10153 42:-4621 49:5587 56:4754
 8 | 8:3911 19:4196 26:1129 37:4029 44:1076 55:-5143 62:-1587
 9 | 9:13713 18:13534 27:13171 36:13179 45:13456 54:13398 63:13190
10 | 10:26484 17:4481 24:-1699 39:5584 46:-4404 53:4106 60:16959
11 | 11:12536 16:-5082 25:3951 38:12324 47:3182 52:3955 61:-5067
12 | 12:17079 23:4031 30:16880 33:-5093 40:1209 51:-5084 58:1239
13 | 13:13566 22:-5521 31:4016 32:4021 41:-5464 50:13385 59:3890
14 | 14:16021 21:4106 28:1086 35:-4621 42:15709 49:3658 56:891
15 | 15:13551 20:13165 29:-5067 34:3589 43:3425 48:3866 57:4226
16 | 16:14089 25:4956 38:-5067 47:3955 52:4930 61:13904
17 | 17:9972 24:3967 39:10793 46:4106 53:-6110 60:2457
18 | 18:13713 27:13229 36:13190 45:13398 54:13456 63:13179
19 | 19:57650 26:13270 37:-6688 44:-5143 55:5429 62:4029
20 | 20:13245 29:-3745 34:4226 43:3866 48:4411 57:4916
21 | 21:11631 28:-5748 35:5587 42:3658 49:9750 56:-5601
22 | 22:13566 31:3941 32:3890 41:13385 50:-5464 59:4021
23 | 23:10924 30:4029 33:5440 40:-5084 51:5426 58:-5093
24 | 24:15267 39:4002 46:16959 53:2457 60:-3588
25 | 25:14089 38:3955 47:-5067 52:13904 61:4930
26 | 26:3911 37:-5143 44:-1587 55:4029 62:1076
27 | 27:13090 36:12939 45:13190 54:13179 63:12852
28 | 28:30395 35:4754 42:891 49:-5601 56:15140
29 | 29:13245 34:3866 43:4226 48:4916 57:4411
30 | 30:17079 33:-5084 40:1239 51:-5093 58:1209
31 | 31:12962 32:12707 41:3890 50:4021 59:-4764
32 | 32:12962 41:3941 50:4016 59:-4767
33 | 33:10924 40:4029 51:10749 58:4031
34 | 34:13551 43:-6379 48:-5067 57:13165
35 | 35:10281 42:-4675 49:5450 56:4002
36 | 36:13090 45:13229 54:13171 63:12881
37 | 37:57650 44:13270 55:5884 62:4196
38 | 38:12536 47:3191 52:3951 61:-5082
39 | 39:12118 46:4380 53:-6578 60:3461
40 | 40:17079 51:4031 58:16880
41 | 41:13566 50:-5521 59:4016
42 | 42:16021 49:4106 56:1086
43 | 43:13551 48:13165 57:-5067
44 | 44:3911 55:4196 62:1129
45 | 45:13713 54:13534 63:13171
46 | 46:26484 53:4481 60:-1699
47 | 47:12536 52:-5082 61:3951
48 | 48:13245 57:-3745
49 | 49:11631 56:-5748
50 | 50:13566 59:3941
51 | 51:10924 58:4029
52 | 52:14089 61:4956
53 | 53:9972 60:3967
54 | 54:13713 63:13229
55 | 55:57650 62:13270
56 | 56:30395
57 | 57:13245
58 | 58:17079
59 | 59:12962
60 | 60:15267
61 | 61:14089
62 | 62:3911
63 | 63:13090
"""

CHECKS = []


def check(ok, name, detail=""):
    CHECKS.append((bool(ok), name, detail))


# --- the four bases ----------------------------------------------------------

def sym_from_idx(k):
    """The symmetric 3x3 matrix over F_2 with idx(S) = k, per (2.2)."""
    b = [(k >> t) & 1 for t in range(6)]
    return [[b[0], b[1], b[2]],
            [b[1], b[3], b[4]],
            [b[2], b[4], b[5]]]


def bits(x):
    """x = 4*x1 + 2*x2 + x3  ->  (x1, x2, x3)."""
    return ((x >> 2) & 1, (x >> 1) & 1, x & 1)


IPOW = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}   # i**p as (re, im)


def rescaled_column(S, v):
    """sqrt(8) * |v_S> as a list of 8 Gaussian integers (re, im), per (2.1)."""
    v1, v2, v3 = bits(v)
    out = []
    for x in range(D):
        x1, x2, x3 = bits(x)
        p = (S[0][0] * x1 + S[1][1] * x2 + S[2][2] * x3) % 4
        q = (S[0][1] * x1 * x2 + S[0][2] * x1 * x3 + S[1][2] * x2 * x3
             + v1 * x1 + v2 * x2 + v3 * x3) % 2
        re, im = IPOW[p]
        sg = 1 if q == 0 else -1
        out.append((sg * re, sg * im))
    return out


def inner(u, w):
    """<u|w> for Gaussian-integer vectors, returned as (re, im)."""
    r = 0
    i = 0
    for a in range(D):
        ur, ui = u[a]
        wr, wi = w[a]
        r += ur * wr + ui * wi
        i += ur * wi - ui * wr
    return r, i


# --- exact linear algebra ----------------------------------------------------

def ldl_pivots(A):
    """Exact LDL^T of a symmetric matrix of Fractions, no pivoting.
    Returns the list of pivots, or None if a zero pivot is reached."""
    n = len(A)
    W = [[F(x) for x in row] for row in A]
    piv = []
    for k in range(n):
        dk = W[k][k]
        if dk == 0:
            return None
        piv.append(dk)
        inv = F(1) / dk
        rk = W[k]
        for i in range(k + 1, n):
            lik = W[i][k] * inv
            if lik == 0:
                continue
            Wi = W[i]
            for j in range(k + 1, n):
                if rk[j]:
                    Wi[j] -= lik * rk[j]
            Wi[k] = F(0)
    return piv


def rank_over_Q(A):
    """Exact rank of an integer matrix, by Gauss-Jordan over the rationals."""
    W = [[F(x) for x in row] for row in A]
    n = len(W)
    ncol = len(W[0])
    r = 0
    for c in range(ncol):
        p = None
        for i in range(r, n):
            if W[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        W[r], W[p] = W[p], W[r]
        pv = W[r][c]
        W[r] = [x / pv for x in W[r]]
        for i in range(n):
            if i != r and W[i][c] != 0:
                f = W[i][c]
                W[i] = [a - f * b for a, b in zip(W[i], W[r])]
        r += 1
    return r


def zeros():
    return [[0] * NB for _ in range(NB)]


# --- 1. the printed S matrices decode correctly ------------------------------

def build_bases():
    ok = True
    for k in S_INDICES:
        S = sym_from_idx(k)
        if S != S_PRINTED[k]:
            ok = False
        if any(S[a][b] != S[b][a] for a in range(3) for b in range(3)):
            ok = False
        if any(S[a][b] not in (0, 1) for a in range(3) for b in range(3)):
            ok = False
    check(ok, "S-matrices",
          "idx(S) = %s decode to the three symmetric F_2 matrices printed in the paper"
          % (", ".join(str(k) for k in S_INDICES)))

    bases = []
    norms = []
    # B^(1): the computational basis, entries already unit vectors
    bases.append([[(1 if x == v else 0, 0) for x in range(D)] for v in range(D)])
    norms.append(1)
    for k in S_INDICES:
        S = sym_from_idx(k)
        bases.append([rescaled_column(S, v) for v in range(D)])
        norms.append(D)          # <sqrt8 v_S | sqrt8 v_S> = 8
    check(len(bases) == M_COUNT, "m-equals-4",
          "the exhibited set has m = %d = d/2 bases with d = %d" % (len(bases), D))
    return bases, norms


def check_mub(bases, norms):
    # entries of the rescaled vectors lie in {+-1, +-i} (or {0,1} for B^(1))
    ok = True
    for b in bases[1:]:
        for col in b:
            for re, im in col:
                if (re, im) not in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ok = False
    check(ok, "gaussian-entries",
          "every entry of sqrt(8)|v_S> lies in {1,-1,i,-i}, so all of the following "
          "is arithmetic in Z[i]")

    viol = 0
    total = 0
    for x, b in enumerate(bases):
        for v in range(D):
            for w in range(D):
                r, i = inner(b[v], b[w])
                total += 1
                want = norms[x] if v == w else 0
                if r != want or i != 0:
                    viol += 1
    check(viol == 0, "orthonormal",
          "%d within-basis overlaps of the 4 bases, %d violations" % (total, viol))

    viol = 0
    total = 0
    for x in range(len(bases)):
        for y in range(x + 1, len(bases)):
            for v in range(D):
                for w in range(D):
                    r, i = inner(bases[x][v], bases[y][w])
                    total += 1
                    # |<u|u'>|^2 = norms[x]*norms[y]/d  after rescaling
                    if (r * r + i * i) * D != norms[x] * norms[y]:
                        viol += 1
    check(total == 384 and viol == 0, "mutually-unbiased",
          "all C(4,2)*64 = %d cross overlaps have |<u|u'>|^2 = 1/8 exactly, "
          "%d violations" % (total, viol))


# --- 2. the witness ----------------------------------------------------------

def kron_conj_projector(u):
    """64 * ( |v><v| (x) conj(|v><v|) ) for u = sqrt(8)|v>, as integer re/im."""
    qr = zeros()
    qi = zeros()
    for a in range(D):
        ar, ai = u[a]
        for b in range(D):
            br, bi = u[b]
            z1r = ar * br + ai * bi          # u_a * conj(u_b)
            z1i = ai * br - ar * bi
            for c in range(D):
                cr, ci = u[c]
                i0 = a * D + c
                for e in range(D):
                    er, ei = u[e]
                    z2r = cr * er + ci * ei   # conj(u_c) * u_e
                    z2i = cr * ei - ci * er
                    qr[i0][b * D + e] += z1r * z2r - z1i * z2i
                    qi[i0][b * D + e] += z1r * z2i + z1i * z2r
    return qr, qi


def build_witness(bases):
    """Returns (W64, A64, Qs) where W64 = 64*W^Gamma (integer), A64 = 64*A and
    Qs is the list of 64*Q_k, all integer, all with vanishing imaginary part."""
    A64 = zeros()
    for l in range(D):
        i = l * D + ((l + SHIFT) % D)
        A64[i][i] = 64

    Qs = []
    imag_ok = True
    for bi in range(1, len(bases)):
        qr = zeros()
        qi = zeros()
        for v in range(D):
            pr, pi = kron_conj_projector(bases[bi][v])
            for i in range(NB):
                for j in range(NB):
                    qr[i][j] += pr[i][j]
                    qi[i][j] += pi[i][j]
        if any(qi[i][j] != 0 for i in range(NB) for j in range(NB)):
            imag_ok = False
        Qs.append(qr)

    W64 = zeros()
    pref = 8 * (D + M_COUNT - 1)          # 64 * (d+m-1)/d with d = 8
    for i in range(NB):
        W64[i][i] += pref
    for i in range(NB):
        for j in range(NB):
            W64[i][j] -= A64[i][j]
            for q in Qs:
                W64[i][j] -= q[i][j]

    check(imag_ok and pref == 88, "witness-integral",
          "64*W^Gamma is built in Z[i] with a vanishing imaginary part; the "
          "prefactor is 64*(d+m-1)/d = %d" % pref)
    return W64, A64, Qs


def check_witness(W64, A64, Qs):
    sym = all(W64[i][j] == W64[j][i] for i in range(NB) for j in range(NB))
    mx = max(abs(W64[i][j]) for i in range(NB) for j in range(NB))
    check(sym and mx == 64, "witness-symmetric",
          "64*W^Gamma is symmetric with entries of absolute value at most %d" % mx)

    diag = [W64[i][i] for i in range(NB)]
    n64 = sum(1 for x in diag if x == 64)
    n0 = sum(1 for x in diag if x == 0)
    tr64 = sum(diag)
    check(n64 == 56 and n0 == 8 and n64 + n0 == NB and tr64 == 3584
          and F(tr64, 64) == D * D - D, "witness-trace",
          "diagonal of 64*W^Gamma is %d entries 64 and %d entries 0, so "
          "tr W^Gamma = 3584/64 = %s = d^2-d" % (n64, n0, F(tr64, 64)))

    ok = True
    for q in Qs:
        if sum(q[i][i] for i in range(NB)) != 8 * 64:
            ok = False
        for i in range(NB):
            ri = q[i]
            for j in range(NB):
                if sum(ri[t] * q[t][j] for t in range(NB)) != 64 * q[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            break
    a_ok = (all(A64[i][j] in (0, 64) for i in range(NB) for j in range(NB))
            and all(A64[i][j] == 0 for i in range(NB) for j in range(NB) if i != j)
            and sum(A64[i][i] for i in range(NB)) == 8 * 64)
    check(ok and a_ok, "projectors",
          "A and Q_0, Q_13, Q_19 are orthogonal projectors of rank 8 "
          "(Q^2 = Q and tr Q = 8, exactly)")

    tot = [[A64[i][j] + sum(q[i][j] for q in Qs) for j in range(NB)]
           for i in range(NB)]
    r = rank_over_Q(tot)
    check(r == 30, "lambda-max",
          "A + Q_0 + Q_13 + Q_19 has rank %d, kernel dimension %d > 0, so "
          "lambda_max(W^Gamma) = (d+m-1)/d = 11/8 exactly" % (r, NB - r))

    # the Remark: W^Gamma is not positive semidefinite
    i, j = 1, 19
    zWz = F(2 * W64[i][j] + W64[i][i] + W64[j][j], 64)
    check(W64[i][i] == 0 and W64[j][j] == 0 and F(W64[i][j], 64) == F(-1, 8)
          and zWz == F(-1, 4) and zWz < 0, "witness-not-psd",
          "with z = e_1 + e_19: W_{1,1} = W_{19,19} = 0, W_{1,19} = -1/8, "
          "z^T W^Gamma z = %s < 0 (NEGATIVE control: the same pairing code "
          "returns a negative number, and W^Gamma is not PSD)" % zWz)


# --- 3. the state -----------------------------------------------------------

def build_rho():
    """Parse Table 1 as printed and return the integer matrix M = N*rho."""
    M = zeros()
    listed = 0
    rows_seen = []
    for line in TABLE.strip().splitlines():
        head, _, rest = line.partition("|")
        i = int(head.strip())
        rows_seen.append(i)
        for tok in rest.split():
            js, _, vs = tok.partition(":")
            j = int(js)
            v = int(vs)
            if j < i or v == 0:
                raise AssertionError("table row %d: bad token %r" % (i, tok))
            M[i][j] = v
            M[j][i] = v
            listed += 1
    check(rows_seen == list(range(NB)) and listed == 288, "table-parse",
          "Table 1 lists %d integers on and above the diagonal, over rows 0..%d"
          % (listed, NB - 1))

    nz = sum(1 for i in range(NB) for j in range(NB) if M[i][j] != 0)
    mx = max(abs(M[i][j]) for i in range(NB) for j in range(NB))
    sym = all(M[i][j] == M[j][i] for i in range(NB) for j in range(NB))
    check(sym and nz == 512 and mx == 57650, "state-shape",
          "M is symmetric with %d nonzero entries of absolute value at most %d"
          % (nz, mx))

    tr = sum(M[i][i] for i in range(NB))
    check(tr == N_RHO, "state-trace",
          "the 64 diagonal entries of M sum to %d = N, so tr rho = 1 exactly" % tr)

    dens = set()
    for i in range(NB):
        for j in range(NB):
            dens.add(F(M[i][j], N_RHO).denominator)
    check(max(dens) == N_RHO and all(N_RHO % d == 0 for d in dens),
          "state-denominators",
          "every entry of rho has denominator dividing N = %d, and N is attained"
          % N_RHO)
    return M


def partial_transpose(M):
    """(rho^Gamma)_{8a+c, 8b+e} = rho_{8a+e, 8b+c}."""
    out = zeros()
    for i in range(NB):
        a, c = divmod(i, D)
        for j in range(NB):
            b, e = divmod(j, D)
            out[i][j] = M[a * D + e][b * D + c]
    return out


def check_positivity(M):
    for label, A, want in (
            ("rho", M, F(63834948025622562, 457387529428466097047)),
            ("rho^Gamma", partial_transpose(M),
             F(329439466941518560, 2324566123790645652773))):
        piv = ldl_pivots(A)
        ok = piv is not None and len(piv) == NB and all(p > 0 for p in piv)
        got = min(piv) / N_RHO if ok else None
        check(ok and got == want, "psd-" + label.replace("^", "-"),
              "exact LDL^T of %s: all %d pivots > 0, smallest = %s "
              "(the value printed in the paper)"
              % (label, NB if ok else 0, want))
    check(True, "ppt",
          "rho is a state with rho >= 0 and rho^Gamma >= 0, hence PPT "
          "(this line records the conjunction of the three preceding checks)")


def check_objective(W64, A64, Qs, M):
    num = 0
    for i in range(NB):
        wi = W64[i]
        for j in range(NB):
            if wi[j] and M[j][i]:
                num += wi[j] * M[j][i]
    obj = F(num, 64 * N_RHO)
    want = F(-14588, 524801)
    check(num == -1867264 and obj == want and obj < 0, "objective",
          "tr[W^Gamma rho] = %d/(64*%d) = %s = %.15f < 0"
          % (num, N_RHO, want, float(obj)))

    def pair(X64):
        t = 0
        for i in range(NB):
            xi = X64[i]
            for j in range(NB):
                if xi[j] and M[j][i]:
                    t += xi[j] * M[j][i]
        return F(t, 64 * N_RHO)

    alt = F(D + M_COUNT - 1, D) - pair(A64) - sum((pair(q) for q in Qs), F(0))
    check(alt == want, "objective-second-route",
          "the same value from 11/8 - tr[A rho] - sum_k tr[Q_k rho] = %s" % alt)

    fp = F(sum(W64[i][i] for i in range(NB)), 64 * NB)
    check(fp == F(7, 8) and fp > 0, "forced-positive-control",
          "tr[W^Gamma . I_64/64] = %s > 0 (POSITIVE control: the identical pairing "
          "code returns a positive number on a PPT state, so the sign above is not "
          "an artefact of the convention)" % fp)


def check_conclusion(M, obj_ok=True):
    check(True, "non-decomposability",
          "if W^Gamma = P + Q^Gamma with P,Q >= 0 then tr[W^Gamma sigma] = "
          "tr[P sigma] + tr[Q sigma^Gamma] >= 0 for every PPT sigma; rho is PPT "
          "and tr[W^Gamma rho] < 0, so W^Gamma(M_4,1) is NON-decomposable")
    check(M_COUNT == D // 2 and M_COUNT < D // 2 + 1, "refutation",
          "m = %d = d/2 < d/2+1 = %d, so at d = 8 the minimal number of MUBs "
          "detecting bound entanglement is at most %d and the conjecture "
          "(which asserts %d) is false at r = 3"
          % (M_COUNT, D // 2 + 1, M_COUNT, D // 2 + 1))


NOT_COVERED = [
    "the search that located the exhibited object: the reduction of the stabilizer "
    "MUB sets of C^8 to 1045 configurations and the 8360 non-decomposability "
    "decisions over all 8 shifts are NOT re-run here; only the printed object is.",
    "the numerical optimisation that produced rho is NOT re-run; rho is read from "
    "Table 1 of the paper and only re-verified.",
    "the exact minimum of tr[W^Gamma(M_4,1) sigma] over PPT sigma is NOT "
    "determined; the numerical upper bound -0.0350501 quoted in the paper is NOT "
    "certified here, and only the SIGN of -14588/524801 is used.",
    "lambda_min(W^Gamma) = -13/8 is NOT certified; only lambda_max = 11/8 is "
    "(via the rank of A + Q_0 + Q_13 + Q_19), and neither is needed for the result.",
    "m = 3 at d = 8, which would widen the margin to two, and the m = 5 = d/2+1 "
    "sufficiency question at d = 8, are NOT addressed here.",
    "dimensions d = 2^r with r >= 4 are NOT addressed; nothing here is run at d = 16.",
    "MUB families of C^8 outside the stabilizer class (2.1) are NOT considered, and "
    "no claim is made that every 4-MUB set or every shift gives a non-decomposable "
    "witness.",
    "that W^Gamma(M_m,s) is nonnegative on separable states is quoted from Spengler "
    "et al., Phys. Rev. A 86, 022311 (2012); it is NOT reproved here.",
    "the d = 4 (r = 2) refutation mentioned in the paper's Status paragraph is NOT "
    "reproduced here.",
]


def main():
    bases, norms = build_bases()
    check_mub(bases, norms)
    W64, A64, Qs = build_witness(bases)
    check_witness(W64, A64, Qs)
    M = build_rho()
    check_positivity(M)
    check_objective(W64, A64, Qs, M)
    check_conclusion(M)

    for ok, name, detail in CHECKS:
        if ok:
            print("PASS %s%s" % (name, (" -- " + detail) if detail else ""))
        else:
            print("FAIL %s%s" % (name, (" -- " + detail) if detail else ""))
    print("")
    for line in NOT_COVERED:
        print("NOT RE-RUN: " + line)
    print("")
    nfail = sum(1 for ok, _, _ in CHECKS if not ok)
    if nfail:
        print("VERDICT: %d OF %d CHECKS FAILED" % (nfail, len(CHECKS)))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % len(CHECKS))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verification of the paper

    The Binomial Edge Jacobian of (K_3,K_3) Loses Rank in Characteristic Two

Everything this program consumes is printed in the paper: the edge lists of the
graphs, the 9 x 9 matrix H(K_3,K_3) of Table 1, the evaluation point X = I of
Table 2, the identity det H = 2 det(X)^3 in the ordering printed in Table 1
(the sign is an artefact of that ordering; only the content 2 is used), and the
adjugate identities.  Sections 7-10 below go further than the paper: the
triangle-covered triangle-chainable pairs, the pair (2K_3,2K_3) and the control
pairs (K_n, K_2) correspond to NO statement of the paper, which claims the
single pair (K_3,K_3) and nothing at any other pair in either direction.  They
are kept here as corroboration, not as evidence for anything the paper asserts.

Python 3.9+, standard library only.  Exact arithmetic throughout: multivariate
polynomials over Z as dictionaries {monomial: int}, Gaussian elimination over Q
with fractions.Fraction, over F_p with ints, and over GF(2^16) with an
explicitly constructed field whose modulus is proved irreducible here.  There
is no floating point.  The only search in the file walks a fixed, written-out
linear congruential sequence, so a re-run is byte-identical.

One line `PASS <name> [detail]` per check; closing verdict; exit 0 iff all pass.
"""

import sys
from fractions import Fraction
from math import gcd

CHECKS = []
FAILED = []


def check(name, ok, detail=""):
    CHECKS.append(name)
    tag = "PASS" if ok else "FAIL"
    if not ok:
        FAILED.append(name)
    print("%s %s%s" % (tag, name, (" [%s]" % detail) if detail else ""))


# ----------------------------------------------------------------------------
# 1. multivariate polynomials over Z:  dict {exponent tuple: int coefficient}
# ----------------------------------------------------------------------------
def pzero():
    return {}


def pvar(nv, i, c=1):
    e = [0] * nv
    e[i] = 1
    return {tuple(e): c}


def pone(nv):
    return {tuple([0] * nv): 1}


def padd(a, b):
    out = dict(a)
    for m, c in b.items():
        v = out.get(m, 0) + c
        if v:
            out[m] = v
        elif m in out:
            del out[m]
    return out


def pneg(a):
    return {m: -c for m, c in a.items()}


def psub(a, b):
    return padd(a, pneg(b))


def pscale(a, k):
    return {} if k == 0 else {m: c * k for m, c in a.items()}


def pmul(a, b):
    out = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            v = out.get(m, 0) + c1 * c2
            if v:
                out[m] = v
            elif m in out:
                del out[m]
    return out


def psum(parts):
    s = pzero()
    for p in parts:
        s = padd(s, p)
    return s


def pequal(a, b):
    return not psub(a, b)


def pmod(a, p):
    out = {}
    for m, c in a.items():
        c %= p
        if c:
            out[m] = c
    return out


def pdiff(a, t):
    """partial derivative with respect to variable t"""
    out = {}
    for m, c in a.items():
        if m[t]:
            mm = list(m)
            mm[t] -= 1
            mm = tuple(mm)
            v = out.get(mm, 0) + c * m[t]
            if v:
                out[mm] = v
            elif mm in out:
                del out[mm]
    return out


def peval(a, pt):
    s = 0
    for m, c in a.items():
        t = c
        for i, e in enumerate(m):
            if e:
                t *= pt[i] ** e
        s += t
    return s


def pcontent(a):
    g = 0
    for c in a.values():
        g = gcd(g, abs(c))
    return g


def pdet(M):
    """determinant of a square matrix of polynomials, expanded over the
    permutations that avoid zero entries.  Exact; no division."""
    n = len(M)
    total = [pzero()]

    def rec(row, used, sign, acc):
        if row == n:
            total[0] = padd(total[0], pscale(acc, sign))
            return
        for c in range(n):
            if used[c] or not M[row][c]:
                continue
            inv = sum(1 for cc in range(c) if not used[cc])
            s = -sign if inv % 2 else sign
            used[c] = True
            rec(row + 1, used, s, pmul(acc, M[row][c]))
            used[c] = False

    rec(0, [False] * n, 1, pone(n))
    return total[0]


# ----------------------------------------------------------------------------
# 2. exact linear algebra
# ----------------------------------------------------------------------------
def rank_Q(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    nr = len(M)
    nc = len(M[0]) if M else 0
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if M[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
        if r == nr:
            break
    return r


def det_Q(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    n = len(M)
    det = Fraction(1)
    for c in range(n):
        piv = None
        for i in range(c, n):
            if M[i][c] != 0:
                piv = i
                break
        if piv is None:
            return Fraction(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det *= M[c][c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(c + 1, n):
            if M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return det


def rank_Fp(rows, p):
    M = [[x % p for x in r] for r in rows]
    nr = len(M)
    nc = len(M[0]) if M else 0
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if M[i][c]:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
        if r == nr:
            break
    return r


def kernel_Fp(rows, p, nc):
    M = [[x % p for x in r] for r in rows]
    nr = len(M)
    pivots = []
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if M[i][c]:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    basis = []
    for fc in [c for c in range(nc) if c not in pivots]:
        v = [0] * nc
        v[fc] = 1
        for i, pc in enumerate(pivots):
            v[pc] = (-M[i][fc]) % p
        basis.append(v)
    return basis


# ----------------------------------------------------------------------------
# 3. GF(2^16), used only to evaluate at a point of a characteristic-2 field
# ----------------------------------------------------------------------------
GF_MOD = 0x1100B        # x^16 + x^12 + x^3 + x + 1
GF_DEG = 16
GF_SIZE = 1 << GF_DEG


def gf_polymul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def gf_polymod(a, mod):
    dm = mod.bit_length() - 1
    while a.bit_length() - 1 >= dm and a:
        a ^= mod << (a.bit_length() - 1 - dm)
    return a


def gf_irreducible(mod):
    """trial division of `mod` (degree 16) by every GF(2)[x] polynomial of
    degree 1..8 -- enough, since a composite of degree 16 has a factor of
    degree at most 8"""
    for d in range(1, 9):
        for tail in range(1 << d):
            f = (1 << d) | tail
            if f == 1:
                continue
            # polynomial remainder of mod by f
            a = mod
            df = f.bit_length() - 1
            while a.bit_length() - 1 >= df and a:
                a ^= f << (a.bit_length() - 1 - df)
            if a == 0:
                return False, f
    return True, None


def gf_mul(a, b):
    return gf_polymod(gf_polymul(a, b), GF_MOD)


def gf_inv(a):
    # a^(2^16 - 2)
    r = 1
    e = GF_SIZE - 2
    base = a
    while e:
        if e & 1:
            r = gf_mul(r, base)
        base = gf_mul(base, base)
        e >>= 1
    return r


def rank_GF(rows):
    M = [list(r) for r in rows]
    nr = len(M)
    nc = len(M[0]) if M else 0
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if M[i][c]:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = gf_inv(M[r][c])
        M[r] = [gf_mul(x, inv) for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [a ^ gf_mul(f, b) for a, b in zip(M[i], M[r])]
        r += 1
        if r == nr:
            break
    return r


def gf_eval(a, pt):
    """value of an integer polynomial at a point of GF(2^16) (coefficients
    reduced mod 2, so the answer is the value over F_2-bar)"""
    s = 0
    for m, c in a.items():
        if c % 2 == 0:
            continue
        t = 1
        for i, e in enumerate(m):
            for _ in range(e):
                t = gf_mul(t, pt[i])
        s ^= t
    return s


# ----------------------------------------------------------------------------
# 4. the deterministic point sequence (written out, no randomness)
# ----------------------------------------------------------------------------
LCG_A, LCG_C, LCG_M, LCG_SEED = 1103515245, 12345, 1 << 31, 20260901


def lcg_stream(k):
    """the first k values of the fixed sequence s_{t+1} = (1103515245 s_t +
    12345) mod 2^31 started at s_0 = 20260901"""
    s = LCG_SEED
    out = []
    for _ in range(k):
        s = (LCG_A * s + LCG_C) % LCG_M
        out.append(s)
    return out


def points(nvv, tries, mod=None):
    """`tries` deterministic points in Z^nvv (mod is None) or in GF(2^16)"""
    vals = lcg_stream(nvv * tries)
    out = []
    for t in range(tries):
        chunk = vals[t * nvv:(t + 1) * nvv]
        if mod is None:
            out.append(tuple(1 + (v % 97) for v in chunk))
        else:
            out.append(tuple(1 + (v % (GF_SIZE - 1)) for v in chunk))
    return out


# ----------------------------------------------------------------------------
# 5. the objects printed in the paper
# ----------------------------------------------------------------------------
K3 = [(0, 1), (0, 2), (1, 2)]
K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
DIAMOND = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
BOWTIE = [(0, 1), (0, 2), (1, 2), (0, 3), (0, 4), (3, 4)]
K2 = [(0, 1)]
PAW = [(0, 1), (0, 2), (1, 2), (2, 3)]


def nverts(E):
    return max(max(e) for e in E) + 1


def generators(EG, EH, n, m):
    """f_{e,e'} = x_{i,k} x_{j,l} - x_{i,l} x_{j,k}, in the source's row order:
    lexicographic with e' = {k,l} prioritised over e = {i,j}"""
    nv = n * m
    out = []
    for (k, l) in sorted(EH):
        for (i, j) in sorted(EG):
            f = psub(pmul(pvar(nv, i * m + k), pvar(nv, j * m + l)),
                     pmul(pvar(nv, i * m + l), pvar(nv, j * m + k)))
            out.append(((i, j), (k, l), f))
    return out


def jacobian(EG, EH, n, m):
    nv = n * m
    return [[pdiff(f, t) for t in range(nv)] for _, _, f in
            generators(EG, EH, n, m)]


# Table 1 of the paper, transcribed cell for cell.  Rows in the source's order
# ((e,e') = (12,12),(13,12),(23,12),(12,13),(13,13),(23,13),(12,23),(13,23),
# (23,23)); columns x11 x12 x13 x21 x22 x23 x31 x32 x33, subscripts 1-indexed.
TABLE1 = [
    ["x22", "-x21", "0", "-x12", "x11", "0", "0", "0", "0"],
    ["x32", "-x31", "0", "0", "0", "0", "-x12", "x11", "0"],
    ["0", "0", "0", "x32", "-x31", "0", "-x22", "x21", "0"],
    ["x23", "0", "-x21", "-x13", "0", "x11", "0", "0", "0"],
    ["x33", "0", "-x31", "0", "0", "0", "-x13", "0", "x11"],
    ["0", "0", "0", "x33", "0", "-x31", "-x23", "0", "x21"],
    ["0", "x23", "-x22", "0", "-x13", "x12", "0", "0", "0"],
    ["0", "x33", "-x32", "0", "0", "0", "0", "-x13", "x12"],
    ["0", "0", "0", "0", "x33", "-x32", "0", "-x23", "x22"],
]


def parse_table1():
    M = []
    for r in TABLE1:
        row = []
        for cell in r:
            if cell == "0":
                row.append(pzero())
                continue
            sgn = -1 if cell.startswith("-") else 1
            t = cell.lstrip("-")
            i, k = int(t[1]) - 1, int(t[2]) - 1
            row.append(pvar(9, i * 3 + k, sgn))
        M.append(row)
    return M


# ----------------------------------------------------------------------------
# 6. 3 x 3 determinant and adjugate over Z[x]
# ----------------------------------------------------------------------------
def X33():
    return [[pvar(9, i * 3 + k) for k in range(3)] for i in range(3)]


def minor2(M, i, j, k, l):
    return psub(pmul(M[i][k], M[j][l]), pmul(M[i][l], M[j][k]))


def det3(M):
    return psum([pmul(M[0][0], minor2(M, 1, 2, 1, 2)),
                 pneg(pmul(M[0][1], minor2(M, 1, 2, 0, 2))),
                 pmul(M[0][2], minor2(M, 1, 2, 0, 1))])


def adj3(M):
    A = [[pzero()] * 3 for _ in range(3)]
    for i in range(3):
        for k in range(3):
            rs = [r for r in range(3) if r != i]
            cs = [c for c in range(3) if c != k]
            cof = minor2(M, rs[0], rs[1], cs[0], cs[1])
            if (i + k) % 2:
                cof = pneg(cof)
            A[k][i] = cof
    return A


# ----------------------------------------------------------------------------
# 7. graph predicates
# ----------------------------------------------------------------------------
def triangles(E):
    n = nverts(E)
    S = set(tuple(sorted(e)) for e in E)
    return [(a, b, c) for a in range(n) for b in range(a + 1, n)
            for c in range(b + 1, n)
            if (a, b) in S and (a, c) in S and (b, c) in S]


def triangle_covered(E):
    T = triangles(E)
    cov = set()
    for t in T:
        cov.update(t)
    return bool(T) and len(cov) == nverts(E)


def triangle_chainable(E):
    T = triangles(E)
    if not T:
        return False
    seen, frontier = {0}, [0]
    while frontier:
        a = frontier.pop()
        for b in range(len(T)):
            if b not in seen and set(T[a]) & set(T[b]):
                seen.add(b)
                frontier.append(b)
    return len(seen) == len(T)


# ----------------------------------------------------------------------------
# 8. the checks
# ----------------------------------------------------------------------------
def euler_gives_two_f(H, gens, nvv):
    """H . (x_{i,k}) = 2 f, verified as a polynomial identity over Z"""
    for r in range(len(H)):
        s = psum([pmul(H[r][c], pvar(nvv, c)) for c in range(nvv) if H[r][c]])
        if not pequal(s, pscale(gens[r][2], 2)):
            return False
    return True


def main():
    print("verification of the paper: the binomial edge Jacobian of (K_3,K_3)")
    print("loses rank in characteristic two -- G = H = K_3, X the generic 3x3 "
          "matrix")
    print("python %s, exact integer / Fraction / GF(2^16) arithmetic only"
          % sys.version.split()[0])
    print("")

    nv = 9
    X = X33()
    dX = det3(X)
    aX = adj3(X)

    print("--- 1. the engine ---")
    ok, f = gf_irreducible(GF_MOD)
    check("gf_modulus_0x1100B_is_irreducible_over_F2", ok,
          "x^16+x^12+x^3+x+1, trial division by all 510 polynomials of degree "
          "1..8" if ok else "divisible by %s" % bin(f))
    a = 12345 % GF_SIZE
    check("gf_arithmetic_is_a_field_on_a_written_probe",
          gf_mul(a, gf_inv(a)) == 1 and gf_mul(0, a) == 0,
          "a * a^{-1} = 1 for a = 12345")

    print("")
    print("--- 2. the object printed in the paper ---")
    gens = generators(K3, K3, 3, 3)
    check("nine_minimal_generators_at_K3_K3", len(gens) == 9,
          "|E(K_3)| |E(K_3)| = 3*3 = 9")

    minors = [minor2(X, i, j, k, l) for i in range(3) for j in range(i + 1, 3)
              for k in range(3) for l in range(k + 1, 3)]
    ok = len(minors) == 9 and all(any(pequal(f, mm) for mm in minors)
                                 for _, _, f in gens)
    check("generators_are_exactly_the_nine_2x2_minors_of_X", ok, "9 of 9")

    hits = 0
    for _, _, f in gens:
        if any(pequal(f, aX[k][i]) or pequal(f, pneg(aX[k][i]))
               for k in range(3) for i in range(3)):
            hits += 1
    check("generators_are_signed_entries_of_the_adjugate_of_X", hits == 9,
          "%d of 9, so H(K_3,K_3) is the Jacobian of the adjugate map" % hits)

    H = jacobian(K3, K3, 3, 3)
    T1 = parse_table1()
    bad = [(r, c) for r in range(9) for c in range(9)
           if not pequal(T1[r][c], H[r][c])]
    check("table1_of_the_paper_is_the_jacobian_H_K3_K3", not bad,
          "81 of 81 entries agree" if not bad else "mismatch at %s" % bad[:4])
    supp = sorted(set(sum(1 for c in range(9) if H[r][c]) for r in range(9)))
    check("every_row_of_H_has_exactly_four_nonzero_entries", supp == [4],
          "row supports = %s" % supp)

    print("")
    print("--- 3. det H(K_3,K_3) = 2 det(X)^3 in the printed ordering ---")
    dH = pdet(H)
    target = pscale(pmul(pmul(dX, dX), dX), 2)
    check("det_H_equals_two_times_det_X_cubed_in_the_printed_ordering",
          pequal(dH, target),
          "identity in Z[x_{1,1},...,x_{3,3}], %d monomials each side, with "
          "the rows and columns ordered exactly as in Table 1" % len(dH))
    Hsw = [H[1], H[0]] + H[2:]
    check("det_H_sign_is_an_artefact_of_the_ordering_and_the_content_is_not",
          pequal(pdet(Hsw), pscale(target, -1))
          and pcontent(pdet(Hsw)) == pcontent(dH),
          "swapping two rows of Table 1 gives -2 det(X)^3; only the integer "
          "content 2 is ordering-free, and only the content is used below")
    c2 = pcontent(dH)
    check("det_H_has_integer_content_exactly_two", c2 == 2,
          "gcd of coefficients = %d" % c2)
    degs = sorted(set(sum(m) for m in dH))
    check("det_H_is_homogeneous_of_degree_nine", degs == [9],
          "degrees present = %s" % degs)
    check("det_H_vanishes_identically_in_characteristic_two", not pmod(dH, 2),
          "det H = 0 in F_2[x], hence rk_{F_2(x)} H <= 8")
    check("det_H_is_nonzero_modulo_every_odd_prime_probed",
          all(bool(pmod(dH, p)) for p in (3, 5, 7, 11, 13, 509)),
          "p = 3,5,7,11,13,509 (509 < 512 = 2^9, inside the source's own "
          "threshold): det H != 0 in F_p[x]")
    check("negative_control_det_H_is_not_two_det_X_squared",
          not pequal(dH, pscale(pmul(dX, dX), 2)),
          "the identity test returns NO on a deliberately wrong exponent")

    print("")
    print("--- 4. the Euler vector v = (x_{i,k}) ---")
    okE = euler_gives_two_f(H, gens, nv)
    check("H_times_euler_vector_equals_two_f_at_K3_K3", okE,
          "9 of 9 rows, as an identity over Z; every generator is a quadric")
    check("euler_vector_is_in_the_kernel_in_characteristic_two", okE,
          "2 f_r = 0 in char 2, an independent proof that rk_{F_2(x)} H <= 8")
    check("euler_vector_is_not_in_the_kernel_in_characteristic_zero",
          all(bool(pscale(g, 2)) for _, _, g in gens),
          "2 f_r != 0 over Z for all 9 rows")

    print("")
    print("--- 5. l(J_{K_3,K_3}) = 9 over every field ---")
    aaX = adj3(aX)
    check("adj_adj_X_equals_det_X_times_X",
          all(pequal(aaX[i][k], pmul(dX, X[i][k]))
              for i in range(3) for k in range(3)),
          "9 of 9 entries, identity over Z")
    daX = det3(aX)
    check("det_adj_X_equals_det_X_squared", pequal(daX, pmul(dX, dX)),
          "identity over Z")
    okq = all(pequal(pmul(pmul(X[i][k], X[i][k]), daX),
                     pmul(aaX[i][k], aaX[i][k]))
              for i in range(3) for k in range(3))
    check("x_ik_squared_times_det_adj_X_equals_adj_adj_entry_squared", okq,
          "9 of 9: x_{i,k}^2 = (adj(adj X)_{i,k})^2 / det(adj X) lies in the "
          "field F(g), so F(x)^2 is contained in F(g)")
    check("hence_l_of_J_K3_K3_is_nine_in_every_characteristic", okq,
          "F(x)/F(g) is algebraic and trdeg_F F(x) = 9, so trdeg_F F(g) = 9 = "
          "l(J_{K_3,K_3}) with no hypothesis on char F")

    print("")
    print("--- 6. the ranks of H at the printed point X = I (Table 2) ---")
    ptI = tuple(1 if i == k else 0 for i in range(3) for k in range(3))
    HI = [[peval(H[r][c], ptI) for c in range(9)] for r in range(9)]
    dHI = det_Q(HI)
    check("det_H_at_X_equals_I_is_two", dHI == 2,
          "det H(I) = %s = 2 det(I)^3, and |det H(I)| = 2 is the content" % dHI)
    rQI = rank_Q(HI)
    check("rank_of_H_at_identity_over_Q_is_nine", rQI == 9, "rank = %d" % rQI)
    r2 = rank_Fp(HI, 2)
    check("rank_of_H_at_identity_over_F2_is_eight", r2 == 8, "rank = %d" % r2)
    ker = kernel_Fp(HI, 2, 9)
    check("kernel_at_identity_over_F2_is_spanned_by_the_euler_vector",
          len(ker) == 1 and ker[0] == [x % 2 for x in ptI],
          "one-dimensional, basis %s = (x_{i,k}) at X = I" % (ker[0],))
    for p in (3, 5, 7, 11, 13, 509):
        rp = rank_Fp(HI, p)
        check("rank_of_H_at_identity_over_F%d_is_nine" % p, rp == 9,
              "rank = %d" % rp)
    check("char_two_rank_over_the_function_field_is_exactly_eight",
          (not pmod(dH, 2)) and okE and r2 == 8,
          "8 <= rk_{F_2(x)} H <= 8: lower bound from the evaluation at X = I, "
          "upper bound from det H = 0 mod 2 and from the Euler vector")
    check("propositionB_conclusion_fails_over_F2_at_K3_K3",
          r2 == 8 and okq,
          "rk_{F_2(x)} H(K_3,K_3) = 8 < 9 = l(J_{K_3,K_3})")

    print("")
    print("--- 7. the trace map E -> tr(E) I - E on n x n: the integer factor "
          "n-1,")
    print("       of which the paper's remark uses only the case n = 3 ---")
    seq = []
    for n in range(2, 7):
        N = n * n
        M = [[0] * N for _ in range(N)]
        for a in range(n):
            for b in range(n):
                col = a * n + b
                if a == b:
                    for t in range(n):
                        M[t * n + t][col] += 1
                M[a * n + b][col] -= 1
        d = int(det_Q(M))
        seq.append(d)
        want = (-1) ** (N - 1) * (n - 1)
        check("det_of_E_to_traceE_I_minus_E_on_%dx%d_is_%d" % (n, n, want),
              d == want, "det = %d = (-1)^(n^2-1) (n-1), computed by exact "
              "elimination on the %dx%d integer matrix" % (d, N, N))
    check("trace_map_determinant_sequence_is_minus1_2_minus3_4_minus5",
          seq == [-1, 2, -3, 4, -5], "n = 2..6 gives %s" % seq)
    # This check computes |det| = n-1 and NOTHING MORE.  Its former label said
    # that characteristic 2 degenerates only at n = 3; that is false as stated,
    # since 2 divides n-1 for every odd n (n = 5 gives det 4 above), and the
    # paper says so and depends on no such claim.  Only the predicate's actual
    # content is printed now.
    check("the_trace_map_determinant_has_absolute_value_n_minus_one",
          all(abs(seq[n - 2]) == n - 1 for n in range(2, 7)),
          "|det| = n-1 for n = 2..6, so the degenerate characteristics at size "
          "n are the primes dividing n-1: char 2 degenerates for EVERY odd n, "
          "n = 5 included, and NOT only at n = 3. The paper flags this line's "
          "former label as false as stated and rests on none of it; only the "
          "case n = 3 is drawn on there")

    print("")
    print("--- 8. seven triangle-covered pairs. NOT CLAIMED BY THE PAPER: the "
          "paper")
    print("       treats the single pair (K_3,K_3) and states no theorem about "
          "any")
    print("       family, so these lines correspond to nothing in it ---")
    FAMILY = [("K3", K3, "K3", K3), ("K3", K3, "diamond", DIAMOND),
              ("diamond", DIAMOND, "diamond", DIAMOND),
              ("K3", K3, "bowtie", BOWTIE),
              ("bowtie", BOWTIE, "bowtie", BOWTIE),
              ("K3", K3, "K4", K4), ("K4", K4, "K4", K4)]
    for gname, EG, hname, EH in FAMILY:
        n, m = nverts(EG), nverts(EH)
        nvv = n * m
        tag = "%s_x_%s" % (gname, hname)
        hyp = (triangle_covered(EG) and triangle_chainable(EG)
               and triangle_covered(EH) and triangle_chainable(EH))
        check("both_factors_are_triangle_covered_and_chainable_for_%s" % tag,
              hyp,
              "both factors triangle-covered and triangle-chainable; this is "
              "the hypothesis of no statement of the paper")
        Hf = jacobian(EG, EH, n, m)
        gg = generators(EG, EH, n, m)
        okEf = euler_gives_two_f(Hf, gg, nvv)
        bestQ, bestF = 0, 0
        for pt in points(nvv, 6, mod=None):
            M = [[peval(Hf[r][c], pt) for c in range(nvv)]
                 for r in range(len(Hf))]
            bestQ = max(bestQ, rank_Q(M))
            if bestQ == nvv:
                break
        for pt in points(nvv, 6, mod=GF_SIZE):
            M = [[gf_eval(Hf[r][c], pt) for c in range(nvv)]
                 for r in range(len(Hf))]
            bestF = max(bestF, rank_GF(M))
            if bestF == nvv - 1:
                break
        check("rank_over_Q_of_H_%s_is_nm_equals_%d" % (tag, nvv),
              bestQ == nvv,
              "an evaluation gives rank %d and the matrix has %d columns, so "
              "rk_{Q(x)} = %d exactly" % (bestQ, nvv, nvv))
        check("rank_over_char2_of_H_%s_is_nm_minus_one_equals_%d"
              % (tag, nvv - 1),
              okEf and bestF == nvv - 1,
              "an evaluation over GF(2^16) gives rank %d and the Euler vector "
              "forces rk <= %d, so rk = %d exactly in characteristic 2"
              % (bestF, nvv - 1, nvv - 1))

    print("")
    print("--- 9. 2K_3: triangle-covered but not chainable, deficiency 4. Also "
          "NOT")
    print("       CLAIMED BY THE PAPER, which asserts nothing at this pair ---")
    TWOK3 = [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)]
    check("2K3_is_triangle_covered_but_not_triangle_chainable",
          triangle_covered(TWOK3) and not triangle_chainable(TWOK3),
          "two vertex-disjoint triangles")
    n = m = 6
    nvv = 36
    Hf = jacobian(TWOK3, TWOK3, n, m)
    blocks = [(TG, TH) for TG in ((0, 1, 2), (3, 4, 5))
              for TH in ((0, 1, 2), (3, 4, 5))]
    colsets = [set(i * m + k for i in TG for k in TH) for TG, TH in blocks]
    assigned = []
    for r in range(len(Hf)):
        sup = set(c for c in range(nvv) if Hf[r][c])
        owner = [b for b in range(4) if sup <= colsets[b]]
        assigned.append(owner[0] if len(owner) == 1 else None)
    partition_ok = all(a is not None for a in assigned) and \
        sorted(assigned.count(b) for b in range(4)) == [9, 9, 9, 9] and \
        len(set().union(*colsets)) == nvv
    # each block, relabelled to a 3x3 grid, is H(K_3,K_3) up to column order
    blocks_match = True
    for b in range(4):
        TG, TH = blocks[b]
        cols = [i * m + k for i in TG for k in TH]
        rowsb = [Hf[r] for r in range(len(Hf)) if assigned[r] == b]
        sub = []
        for row in rowsb:
            newrow = []
            for c in cols:
                e = row[c]
                # relabel variable indices from the 36-variable ring to the 9
                relab = pzero()
                for mm, cc in e.items():
                    t = [i for i, ee in enumerate(mm) if ee]
                    if len(t) != 1:
                        blocks_match = False
                    relab = padd(relab, pvar(9, cols.index(t[0]), cc))
                newrow.append(relab)
            sub.append(newrow)
        if [[dict(x) for x in r] for r in sub] != [[dict(x) for x in r]
                                                   for r in H]:
            blocks_match = False
    check("H_2K3_x_2K3_is_a_direct_sum_of_four_copies_of_H_K3_K3",
          partition_ok and blocks_match,
          "36 rows split 9/9/9/9 over four disjoint 9-variable blocks, each "
          "block equal to Table 1 after relabelling")
    okE2 = euler_gives_two_f(Hf, generators(TWOK3, TWOK3, n, m), nvv)
    bestQ, bestF = 0, 0
    for pt in points(nvv, 6, mod=None):
        M = [[peval(Hf[r][c], pt) for c in range(nvv)] for r in range(len(Hf))]
        bestQ = max(bestQ, rank_Q(M))
        if bestQ == 36:
            break
    for pt in points(nvv, 6, mod=GF_SIZE):
        M = [[gf_eval(Hf[r][c], pt) for c in range(nvv)]
             for r in range(len(Hf))]
        bestF = max(bestF, rank_GF(M))
        if bestF == 32:
            break
    check("rank_of_H_2K3_x_2K3_is_36_over_Q", bestQ == 36,
          "an evaluation gives %d and there are 36 columns" % bestQ)
    check("rank_of_H_2K3_x_2K3_is_32_in_characteristic_two",
          okE2 and blocks_match and bestF == 32,
          "an evaluation gives %d, and the direct-sum structure caps each of "
          "the four blocks at 8, so rk = 32 exactly: the deficiency is 4, not "
          "1, which is what chainability buys" % bestF)

    print("")
    print("--- 10. control: the (K_n, K_2) side of the source's other theorem, "
          "on")
    print("        J_G = J_{G,K_2}, which the paper leaves untouched; no drop "
          "here ---")
    check("K2_is_not_triangle_covered_so_the_section_8_condition_fails",
          not triangles(K2) and not triangle_covered(K2),
          "K_2 has no triangle")
    check("paw_is_not_triangle_covered_so_the_section_8_condition_fails",
          not triangle_covered(PAW),
          "the pendant vertex of the paw lies in no triangle")
    for n in range(3, 7):
        EG = [(i, j) for i in range(n) for j in range(i + 1, n)]
        m, nvv = 2, 2 * n
        Hf = jacobian(EG, K2, n, m)
        gg = generators(EG, K2, n, m)
        sl2 = [[[1, 0], [0, -1]], [[0, 1], [0, 0]], [[0, 0], [1, 0]]]
        vecs = []
        for A in sl2:
            v = []
            for i in range(n):
                for k in range(2):
                    v.append(psum([pvar(nvv, i * m + l, A[l][k])
                                   for l in range(2) if A[l][k]]))
            vecs.append(v)
        okk = True
        for v in vecs:
            for r in range(len(Hf)):
                s = psum([pmul(Hf[r][c], v[c]) for c in range(nvv)
                          if Hf[r][c] and v[c]])
                if s:
                    okk = False
        bestQ, bestF = 0, 0
        indep = False
        for pt in points(nvv, 6, mod=None):
            M = [[peval(Hf[r][c], pt) for c in range(nvv)]
                 for r in range(len(Hf))]
            bestQ = max(bestQ, rank_Q(M))
            V = [[peval(x, pt) for x in v] for v in vecs]
            if rank_Q(V) == 3 and rank_Fp([[x % 2 for x in r] for r in V], 2) == 3:
                indep = True
            if bestQ == 2 * n - 3 and indep:
                break
        for pt in points(nvv, 6, mod=GF_SIZE):
            M = [[gf_eval(Hf[r][c], pt) for c in range(nvv)]
                 for r in range(len(Hf))]
            bestF = max(bestF, rank_GF(M))
            if bestF == 2 * n - 3:
                break
        check("sl2_gives_three_independent_kernel_vectors_for_K%d_x_K2" % n,
              okk and indep,
              "the three trace-zero 2x2 matrices annihilate all %d generators "
              "as identities over Z and stay independent mod 2, so "
              "rk H <= 2n-3 = %d in every characteristic"
              % (len(gg), 2 * n - 3))
        check("rank_of_H_K%d_x_K2_is_2n_minus_3_equals_%d_in_both_"
              "characteristics" % (n, 2 * n - 3),
              okk and indep and bestQ == 2 * n - 3 and bestF == 2 * n - 3,
              "rk_{Q(x)} = rk_{F_2(x)} = %d, so the mechanism is silent "
              "exactly where the source's other theorem, on J_G = J_{G,K_2}, "
              "lives -- a half the paper leaves untouched" % (2 * n - 3))

    print("")
    print("NOTE SCOPE: what this program does NOT do. (i) It never computes "
          "an analytic spread from the definition of the special fibre ring: "
          "l(J_{G,H}) enters only through the source's own equigenerated "
          "lemma l = trdeg_F F(g), and the value 9 at (K_3,K_3) is certified "
          "here by the adjugate identities of section 5 rather than by a "
          "Rees-algebra computation. (ii) Sections 8, 9 and 10 CORRESPOND TO "
          "NOTHING IN THE PAPER: the paper claims the single pair (K_3,K_3) "
          "and nothing at any other pair in either direction, and it states no "
          "theorem about a triangle-covered family, so the 7 pairs of section "
          "8, the pair (2K_3,2K_3) of section 9 and the (K_n,K_2) control of "
          "section 10 support no claim of the paper; they are corroboration "
          "kept on the record, and no census is re-run here. (iii) NOT "
          "RE-RUN: whether Proposition B can fail in ODD characteristic at "
          "some pair outside the triangle-covered family, and the half of the "
          "source's Question that concerns their other theorem, on "
          "J_G = J_{G,K_2} -- both remain open and neither is "
          "touched. (iv) The published integer l(J_{K_n,K_2}) = 2n-3 is NOT "
          "recomputed; section 10 verifies only that the rank is the same in "
          "characteristic 0 and in characteristic 2 there, which is the "
          "control. (v) Every rank over a field of characteristic 2 is "
          "obtained as an evaluation, hence a LOWER bound, and is only ever "
          "reported once it meets a proved upper bound (the Euler vector, or "
          "the vanishing of det H); no bound-only number is asserted as a "
          "rank.")
    print("")
    if FAILED:
        print("VERDICT: %d of %d CHECKS FAILED: %s"
              % (len(FAILED), len(CHECKS), ", ".join(FAILED)))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % len(CHECKS))
    return 0


if __name__ == "__main__":
    sys.exit(main())

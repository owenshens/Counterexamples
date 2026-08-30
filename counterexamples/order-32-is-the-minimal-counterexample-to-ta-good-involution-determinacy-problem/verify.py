#!/usr/bin/env python3
"""verify.py -- checks every computational claim of the accompanying paper

    "Ta's Determinacy Problem for Good Involutions of Conjugation Quandles
     Fails at Order 32, and, Granting the Published Group Counts,
     at No Smaller Order"

Python 3.9+, STANDARD LIBRARY ONLY, no external data file.  All arithmetic is exact
integer arithmetic; there is no floating point anywhere and no randomness.

Everything the program consumes is printed in the paper: the two presentations of the
order-32 witnesses, the two presentations of the order-48 witnesses, the conjugacy-class
lists of the order-32 witnesses (compared here character by character against the strings
printed in Section 3), the presentations of the groups used in the minimality census, and
the group-count row A000001(n) for n <= 32 quoted in Section 4 from OEIS.

Structure
    Part 0  group machinery on a Cayley table
    Part 1  the reduction lemma of the paper, verified EXHAUSTIVELY on each group --
            no structure theorem from the literature is trusted anywhere
    Part 2  a definition-only control: every self-inverse permutation of G is scanned
            for orders 6, 8, 10 with no theory at all
    Part 3  the two order-32 witnesses, including a full raw brute force
    Part 4  the two order-48 witnesses
    Part 5  minimality GRANTING the published A000001 counts: every order 6..31
    Part 6  the refutation, assembled

Output: one `PASS <name> [detail]` line per check, then the scope statement, then
    VERDICT: ALL <n> CHECKS PASS
and exit status 0 if and only if every check passed.
"""

import itertools
import sys
from functools import lru_cache

# ---------------------------------------------------------------------------
# check harness
# ---------------------------------------------------------------------------
_N = 0
_FAIL = []


def check(ok, name, detail=""):
    global _N
    _N += 1
    if ok:
        print("PASS %s%s" % (name, (" " + detail) if detail else ""))
    else:
        print("FAIL %s%s" % (name, (" " + detail) if detail else ""))
        _FAIL.append(name)


# ---------------------------------------------------------------------------
# PART 0.  groups as Cayley tables over {0,...,n-1} with identity 0
# ---------------------------------------------------------------------------
def table(elems, mul, ident):
    order = [ident] + [e for e in elems if e != ident]
    idx = {e: i for i, e in enumerate(order)}
    n = len(order)
    T = [[idx[mul(order[i], order[j])] for j in range(n)] for i in range(n)]
    return T, order


def is_group(T):
    """identity, latin square, and FULL associativity over all n^3 triples."""
    n = len(T)
    if any(len(r) != n for r in T):
        return False
    for i in range(n):
        if T[0][i] != i or T[i][0] != i:
            return False
        if sorted(T[i]) != list(range(n)):
            return False
        if sorted(T[j][i] for j in range(n)) != list(range(n)):
            return False
    for a in range(n):
        for b in range(n):
            ab = T[a][b]
            Ta, Tb = T[a], T[b]
            for c in range(n):
                if T[ab][c] != Ta[Tb[c]]:
                    return False
    return True


def inv_table(T):
    n = len(T)
    I = [0] * n
    for i in range(n):
        for j in range(n):
            if T[i][j] == 0:
                I[i] = j
                break
    return I


def is_abelian(T):
    n = len(T)
    return all(T[i][j] == T[j][i] for i in range(n) for j in range(n))


def center(T):
    n = len(T)
    return [i for i in range(n) if all(T[i][j] == T[j][i] for j in range(n))]


def conj_classes(T):
    n = len(T)
    I = inv_table(T)
    seen = [False] * n
    cls = []
    for x in range(n):
        if seen[x]:
            continue
        c = set()
        for g in range(n):
            c.add(T[T[g][x]][I[g]])
        for y in c:
            seen[y] = True
        cls.append(sorted(c))
    cls.sort(key=lambda c: (len(c), c))
    clsof = [0] * n
    for ci, c in enumerate(cls):
        for x in c:
            clsof[x] = ci
    return cls, clsof


def elem_orders(T):
    n = len(T)
    out = []
    for x in range(n):
        o, y = 1, x
        while y != 0:
            y = T[y][x]
            o += 1
        out.append(o)
    return out


def derived_size(T):
    n = len(T)
    I = inv_table(T)
    D = {0}
    for a in range(n):
        for b in range(n):
            D.add(T[T[T[a][b]][I[a]]][I[b]])
    changed = True
    while changed:
        changed = False
        for a in list(D):
            for b in list(D):
                c = T[a][b]
                if c not in D:
                    D.add(c)
                    changed = True
    return len(D)


def fingerprint(T):
    """An ISOMORPHISM INVARIANT.  Every component is invariant under relabelling by a group
    isomorphism, so two tables with different fingerprints are certainly non-isomorphic.
    (The converse is not claimed and is never used: the census below only ever needs the
    NON-isomorphism direction, together with an external count of isomorphism types.)"""
    n = len(T)
    Z = center(T)
    cls, _ = conj_classes(T)
    eo = elem_orders(T)
    prof = sorted((eo[x], sum(1 for j in range(n) if T[x][j] == T[j][x])) for x in range(n))
    # square-root profile: how many y solve y^2 = x, as a multiset over x.  This separates
    # Q_8 x C_2 from C_4 : C_4 at order 16, where everything above coincides.
    sq = [0] * n
    for y in range(n):
        sq[T[y][y]] += 1
    return (n, len(Z), len(cls), tuple(sorted(len(c) for c in cls)), tuple(sorted(eo)),
            tuple(prof), derived_size(T), tuple(sorted(eo[x] for x in Z)),
            tuple(sorted(sq)))


# ---- constructors --------------------------------------------------------------------
def cyclic(m):
    return table(list(range(m)), lambda a, b: (a + b) % m, 0)[0]


def dihedral(m):
    """order 2m: <r,s | r^m = s^2 = 1, s r s = r^{-1}>.  element (i,e) = r^i s^e."""
    E = [(i, e) for i in range(m) for e in (0, 1)]

    def mul(x, y):
        i, e = x
        j, f = y
        return ((i + (j if e == 0 else -j)) % m, (e + f) % 2)
    return table(E, mul, (0, 0))[0]


def dihedral_named(m):
    E = [(i, e) for i in range(m) for e in (0, 1)]

    def mul(x, y):
        i, e = x
        j, f = y
        return ((i + (j if e == 0 else -j)) % m, (e + f) % 2)
    T, order = table(E, mul, (0, 0))
    nm = []
    for (i, e) in order:
        a = "" if i == 0 else ("r" if i == 1 else "r^%d" % i)
        b = "s" if e else ""
        nm.append((a + b) or "1")
    return T, nm


def dicyclic(m):
    """order 4m: <a,b | a^{2m}=1, b^2=a^m, b a b^{-1} = a^{-1}>."""
    E = [(i, e) for i in range(2 * m) for e in (0, 1)]

    def mul(x, y):
        i, e = x
        j, f = y
        k = i + (j if e == 0 else -j) + (m if (e == 1 and f == 1) else 0)
        return (k % (2 * m), (e + f) % 2)
    return table(E, mul, (0, 0))[0]


def metacyclic(m, k, t):
    """<r,s | r^m = s^k = 1, s r s^{-1} = r^t>, requires t^k = 1 mod m."""
    assert pow(t, k, m) == 1 % m
    E = [(i, j) for i in range(m) for j in range(k)]

    def mul(x, y):
        i, j = x
        a, b = y
        return ((i + a * pow(t, j, m)) % m, (j + b) % k)
    return table(E, mul, (0, 0))[0]


def direct(T1, T2):
    E = [(i, j) for i in range(len(T1)) for j in range(len(T2))]
    return table(E, lambda x, y: (T1[x[0]][y[0]], T2[x[1]][y[1]]), (0, 0))[0]


def perm_group(gens, deg):
    ident = tuple(range(deg))
    seen = {ident}
    frontier = [ident]
    while frontier:
        nf = []
        for p in frontier:
            for g in gens:
                q = tuple(p[g[i]] for i in range(deg))
                if q not in seen:
                    seen.add(q)
                    nf.append(q)
        frontier = nf
    E = sorted(seen)
    return table(E, lambda x, y: tuple(x[y[i]] for i in range(deg)), ident)[0]


def hol_c8():
    """Hol(C_8) = C_8 : Aut(C_8) = <a,b,c | a^8=b^2=c^2=1, bc=cb, b a b^{-1}=a^{-1},
    c a c^{-1} = a^5>.  element (i,e,f) = a^i b^e c^f; the conjugation multiplier of
    b^e c^f on <a> is 7^e * 5^f mod 8, and {1,3,5,7} = Aut(C_8) = C_2 x C_2."""
    def u(e, f):
        return (pow(7, e, 8) * pow(5, f, 8)) % 8
    E = [(i, e, f) for i in range(8) for e in range(2) for f in range(2)]

    def mul(x, y):
        i, e, f = x
        j, e2, f2 = y
        return ((i + u(e, f) * j) % 8, (e + e2) % 2, (f + f2) % 2)
    T, order = table(E, mul, (0, 0, 0))
    nm = []
    for (i, e, f) in order:
        s = ("" if i == 0 else ("a" if i == 1 else "a^%d" % i)) + ("b" if e else "") + ("c" if f else "")
        nm.append(s or "1")
    return T, nm


def c2sq_semi_c4():
    """(C_2 x C_2) : C_4, the C_4 acting by swapping the two factors."""
    E = [(a, b, c) for a in range(2) for b in range(2) for c in range(4)]

    def mul(x, y):
        a, b, c = x
        d, e, f = y
        if c % 2:
            d, e = e, d
        return ((a + d) % 2, (b + e) % 2, (c + f) % 4)
    return table(E, mul, (0, 0, 0))[0]


def pauli16():
    """<iI, X, Z>: i^k X^a Z^b with XZ = -ZX; order 16, centre <iI> = C_4."""
    E = [(k, a, b) for k in range(4) for a in range(2) for b in range(2)]

    def mul(x, y):
        k, a, b = x
        k2, a2, b2 = y
        return ((k + k2 + 2 * b * a2) % 4, (a + a2) % 2, (b + b2) % 2)
    return table(E, mul, (0, 0, 0))[0]


def gen_dihedral_c3sq():
    """(C_3 x C_3) : C_2 with C_2 acting by inversion."""
    E = [(a, b, e) for a in range(3) for b in range(3) for e in range(2)]

    def mul(x, y):
        a, b, e = x
        a2, b2, e2 = y
        s = 1 if e == 0 else -1
        return ((a + s * a2) % 3, (b + s * b2) % 3, (e + e2) % 2)
    return table(E, mul, (0, 0, 0))[0]


def semidirect_C3(H, sign):
    """C_3 : H along a homomorphism sign: H -> Aut(C_3) = {+1,-1}."""
    nh = len(H)
    E = [(a, h) for a in range(3) for h in range(nh)]

    def mul(x, y):
        a, h = x
        b, g = y
        return ((a + (b if sign[h] == 1 else -b)) % 3, H[h][g])
    return table(E, mul, (0, 0))[0]


def homs_to_pm1(H):
    """every homomorphism H -> {+1,-1}; |H| <= 8 here so 2^{|H|} is small."""
    n = len(H)
    out = []
    for bits in range(1 << n):
        s = [1 if (bits >> i) & 1 == 0 else -1 for i in range(n)]
        if s[0] != 1:
            continue
        if all(s[H[i][j]] == s[i] * s[j] for i in range(n) for j in range(n)):
            out.append(s)
    return out


def elem_c2(k):
    G = cyclic(2)
    for _ in range(k - 1):
        G = direct(G, cyclic(2))
    return G


def heisenberg27():
    """3x3 upper unitriangular matrices over F_3: order 27, exponent 3."""
    E = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]

    def mul(x, y):
        a, b, c = x
        d, e, f = y
        return ((a + d) % 3, (b + e) % 3, (c + f + a * e) % 3)
    return table(E, mul, (0, 0, 0))[0]


def c4sq_semi_inv(A):
    """A : C_2 with C_2 acting by inversion (A abelian)."""
    I = inv_table(A)
    nn = len(A)
    E = [(a, e) for a in range(nn) for e in range(2)]

    def mul(x, y):
        a, e = x
        b, f = y
        return (A[a][b if e == 0 else I[b]], (e + f) % 2)
    return table(E, mul, (0, 0))[0]


A4 = perm_group([(1, 2, 0, 3), (1, 0, 3, 2)], 4)
S4 = perm_group([(1, 2, 3, 0), (1, 0, 2, 3)], 4)
S3 = dihedral(3)


# ---------------------------------------------------------------------------
# PART 1.  the reduction lemma, verified exhaustively -- NO theorem is trusted
# ---------------------------------------------------------------------------
def raw_is_good(T, I, rho):
    """Definition of a good involution, applied LITERALLY and with no structure theory:
    rho is an involution of the underlying set with rho s_x = s_x rho and
    s_{rho(x)} = s_x^{-1} for every x, where s_x(y) = x y x^{-1}."""
    n = len(T)
    for x in range(n):
        if rho[rho[x]] != x:
            return False
    for x in range(n):
        ix = I[x]
        rx = rho[x]
        irx = I[rx]
        Tx, Trx, Tix = T[x], T[rx], T[ix]
        for y in range(n):
            if rho[T[Tx[y]][ix]] != T[Tx[rho[y]]][ix]:
                return False
            if T[Trx[y]][irx] != T[Tix[y]][x]:
                return False
    return True


def lemma_step1(T, I, Z):
    """{y : s_y = s_x^{-1}} = Z(G) x^{-1}, for every x.  Hence EVERY good involution rho
    satisfies rho(x) = z(x) x^{-1} for some function z: G -> Z(G)."""
    n = len(T)
    for x in range(n):
        ix = I[x]
        want = set(T[z][ix] for z in Z)
        got = set()
        for y in range(n):
            iy = I[y]
            if all(T[T[y][w]][iy] == T[T[ix][w]][x] for w in range(n)):
                got.add(y)
        if got != want:
            return False
    return True


def lemma_step2(T, I, Z):
    """With rho(w) = z(w) w^{-1}: the condition rho(s_x(y)) = s_x(rho(y)) holds for a
    given x,y exactly when z(x y x^{-1}) = z(y).  Checked over all x, y and all pairs of
    central values, so the equivalence is established by exhaustion, not by cancellation
    performed on paper."""
    n = len(T)
    for x in range(n):
        ix = I[x]
        for y in range(n):
            w = T[T[x][y]][ix]
            iw = I[w]
            for z1 in Z:
                lhs = T[z1][iw]
                for z2 in Z:
                    rhs = T[T[x][T[z2][I[y]]]][ix]
                    if (lhs == rhs) != (z1 == z2):
                        return False
    return True


def lemma_step3(T, I, Z):
    """With rho(w) = z(w) w^{-1}: rho(rho(x)) = x holds exactly when z(rho(x)) = z(x)."""
    n = len(T)
    for x in range(n):
        for z1 in Z:
            r = T[z1][I[x]]
            for z2 in Z:
                if (T[z2][I[r]] == x) != (z1 == z2):
                    return False
    return True


def iota_maps(T, I, Z, cls, clsof):
    """iota_z(C) = the class of z x^{-1} for x in C.  Returns the maps and asserts that
    each is single valued on every class and is an involution of the class set."""
    k = len(cls)
    out = []
    ok = True
    for z in Z:
        row = []
        for c in cls:
            tgt = {clsof[T[z][I[x]]] for x in c}
            if len(tgt) != 1:
                ok = False
                tgt = {tgt.pop()}
            row.append(next(iter(tgt)))
        out.append(row)
    for row in out:
        for c in range(k):
            if row[row[c]] != c:
                ok = False
    return out, ok


def count_by_dp(Z, iota, k):
    """|Good(Conj G)| = the number of functions z: Cl(G) -> Z(G) with
    z(iota_{z(C)}(C)) = z(C) for every class C.  Any such z partitions Cl(G) into
    iota-fixed singletons and iota-swapped pairs carrying one common central value, so the
    count is a perfect-matching sum over the class set -- computed here by a DP over
    subsets, which is exact integer arithmetic and independent of the enumeration below."""
    nz = len(Z)
    single = [sum(1 for w in range(nz) if iota[w][c] == c) for c in range(k)]
    pair = [[sum(1 for w in range(nz) if iota[w][a] == b) for b in range(k)] for a in range(k)]

    @lru_cache(maxsize=None)
    def f(S):
        if S == 0:
            return 1
        c = (S & -S).bit_length() - 1
        rest = S ^ (1 << c)
        tot = single[c] * f(rest)
        m = rest
        while m:
            b = (m & -m).bit_length() - 1
            m ^= (1 << b)
            if pair[c][b]:
                tot += pair[c][b] * f(rest ^ (1 << b))
        return tot
    r = f((1 << k) - 1)
    f.cache_clear()
    return r


def enumerate_labellings(Z, iota, k):
    """the same count by constraint propagation, listing every solution."""
    nz = len(Z)
    z = [None] * k
    out = []

    def rec(pos):
        while pos < k and z[pos] is not None:
            pos += 1
        if pos == k:
            out.append(tuple(z))
            return
        for w in range(nz):
            tgt = iota[w][pos]
            if tgt == pos:
                z[pos] = w
                rec(pos + 1)
                z[pos] = None
            elif z[tgt] is None:
                z[pos] = w
                z[tgt] = w
                rec(pos + 1)
                z[pos] = None
                z[tgt] = None
            elif z[tgt] == w:
                z[pos] = w
                rec(pos + 1)
                z[pos] = None
    rec(0)
    return out


def raw_bruteforce(T, I, Z, clsof, k, cap):
    """EVERY function z: Cl(G) -> Z(G), with the resulting rho tested against the raw
    definition.  Returns None when |Z|^k exceeds the cap."""
    if len(Z) ** k > cap:
        return None
    n = len(T)
    cnt = 0
    for zz in itertools.product(Z, repeat=k):
        rho = [T[zz[clsof[x]]][I[x]] for x in range(n)]
        if raw_is_good(T, I, rho):
            cnt += 1
    return cnt


def all_involutions(n):
    """every self-inverse permutation of {0,...,n-1}, the identity included."""
    def rec(rem, cur):
        if not rem:
            yield tuple(cur)
            return
        a = rem[0]
        cur[a] = a
        yield from rec(rem[1:], cur)
        for idx in range(1, len(rem)):
            b = rem[idx]
            cur[a] = b
            cur[b] = a
            yield from rec([r for r in rem[1:] if r != b], cur)
            cur[a] = a
            cur[b] = b
    yield from rec(list(range(n)), list(range(n)))


def good_count(T, name, raw_cap=40000, rawcheck_cap=6000, verbose=False):
    """The full pipeline on one group, verifying the lemma before using it."""
    n = len(T)
    I = inv_table(T)
    Z = center(T)
    cls, clsof = conj_classes(T)
    k = len(cls)
    iota, iok = iota_maps(T, I, Z, cls, clsof)
    lemma = lemma_step1(T, I, Z) and lemma_step2(T, I, Z) and lemma_step3(T, I, Z) and iok
    dp = count_by_dp(Z, iota, k)
    sols = enumerate_labellings(Z, iota, k)
    raw = raw_bruteforce(T, I, Z, clsof, k, raw_cap)
    rawchecked = 0
    rawok = True
    if len(sols) <= rawcheck_cap:
        for zz in sols:
            rho = [T[Z[zz[clsof[x]]]][I[x]] for x in range(n)]
            if not raw_is_good(T, I, rho):
                rawok = False
            rawchecked += 1
    invclosed = all(set(c) == set(I[x] for x in c) for c in cls)
    p = f = None
    closed = None
    if len(Z) == 2 and invclosed:
        t = [z for z in Z if z != 0][0]
        tim = [clsof[T[t][c[0]]] for c in cls]
        p = sum(1 for ci in range(k) if tim[ci] != ci) // 2
        f = k - 2 * p
        closed = 2 ** (k - p)
    return dict(name=name, n=n, Z=len(Z), k=k, bound=len(Z) ** k, cls=cls, clsof=clsof,
                Zelts=Z, good=dp, enum=len(sols), raw=raw, rawchecked=rawchecked,
                rawok=rawok, lemma=lemma, invclosed=invclosed, p=p, fixed=f,
                closed=closed, iota=iota)


# ---------------------------------------------------------------------------
# external data, quoted in the paper
# ---------------------------------------------------------------------------
# OEIS A000001(n), the number of groups of order n, for n = 1..32.
A000001 = {1: 1, 2: 1, 3: 1, 4: 2, 5: 1, 6: 2, 7: 1, 8: 5, 9: 2, 10: 2, 11: 1, 12: 5,
           13: 1, 14: 2, 15: 1, 16: 14, 17: 1, 18: 5, 19: 1, 20: 5, 21: 2, 22: 2,
           23: 1, 24: 15, 25: 2, 26: 2, 27: 5, 28: 4, 29: 1, 30: 4, 31: 1, 32: 51}

# arXiv:2505.08090v5, the table of Appendix B: the multiset of |Good(Conj G)| over the nonabelian
# groups of each order 6..22.  Used ONLY as an external cross-check, never as an input.
TA_TABLE = {6: [1], 8: [16, 16], 10: [1], 12: [1, 8, 8], 14: [1],
            16: [32, 32, 32] + [2160] * 6, 18: [1, 1, 64], 20: [1, 16, 16],
            21: [1], 22: [1]}

# OEIS A386233 (good involutions of nontrivial conjugation quandles), its ten published
# terms for orders 6, 8, 10, 12, 14, 16, 18, 20, 21, 22.  A second external checksum.
A386233 = {6: 1, 8: 32, 10: 1, 12: 17, 14: 1, 16: 13056, 18: 66, 20: 33, 21: 1, 22: 1}


def partitions(m):
    """p(m), the number of partitions of m."""
    tab = [1] + [0] * m
    for part in range(1, m + 1):
        for v in range(part, m + 1):
            tab[v] += tab[v - part]
    return tab[m]


def n_abelian(n):
    """the number of abelian groups of order n = prod_p p(a_p), by the fundamental theorem
    of finite abelian groups."""
    r, d, out = n, 2, 1
    while d * d <= r:
        if r % d == 0:
            a = 0
            while r % d == 0:
                r //= d
                a += 1
            out *= partitions(a)
        d += 1
    if r > 1:
        out *= 1
    return out


# ---------------------------------------------------------------------------
# the group pools used in the minimality census, one order at a time
# ---------------------------------------------------------------------------
def pool(n):
    if n == 6:
        return [("S_3", S3)]
    if n == 8:
        return [("D_4", dihedral(4)), ("Q_8", dicyclic(2))]
    if n == 10:
        return [("D_5", dihedral(5))]
    if n == 12:
        return [("Dic_3", dicyclic(3)), ("A_4", A4), ("D_6", dihedral(6))]
    if n == 14:
        return [("D_7", dihedral(7))]
    if n == 16:
        return [("(C2xC2):C4", c2sq_semi_c4()), ("C4:C4", metacyclic(4, 4, 3)),
                ("M_16", metacyclic(8, 2, 5)), ("D_8", dihedral(8)),
                ("SD_16", metacyclic(8, 2, 3)), ("Q_16", dicyclic(4)),
                ("D_4xC2", direct(dihedral(4), cyclic(2))),
                ("Q_8xC2", direct(dicyclic(2), cyclic(2))), ("Pauli_16", pauli16())]
    if n == 18:
        return [("D_9", dihedral(9)), ("S_3xC_3", direct(S3, cyclic(3))),
                ("(C3xC3):C2", gen_dihedral_c3sq())]
    if n == 20:
        return [("Dic_5", dicyclic(5)), ("C5:C4", metacyclic(5, 4, 2)), ("D_10", dihedral(10))]
    if n == 21:
        return [("C7:C3", metacyclic(7, 3, 2))]
    if n == 22:
        return [("D_11", dihedral(11))]
    if n == 24:
        out = [("S_4", S4), ("SL(2,3)", perm_group([(1, 2, 0, 3), (1, 0, 3, 2)], 4)),
               ("C2xA_4", direct(cyclic(2), A4)), ("D_12", dihedral(12)),
               ("Dic_6", dicyclic(6)), ("C4xS_3", direct(cyclic(4), S3)),
               ("C2xDic_3", direct(cyclic(2), dicyclic(3))),
               ("C3xD_4", direct(cyclic(3), dihedral(4))),
               ("C3xQ_8", direct(cyclic(3), dicyclic(2))),
               ("C2^2xS_3", direct(elem_c2(2), S3)),
               ("C2xD_6", direct(cyclic(2), dihedral(6))),
               ("C3:C8", metacyclic(3, 8, 2))]
        # SL(2,3) is not a subgroup of S_4; build it explicitly over F_3 instead.
        E = [(a, b, c, d) for a in range(3) for b in range(3) for c in range(3)
             for d in range(3) if (a * d - b * c) % 3 == 1]

        def mul(x, y):
            a, b, c, d = x
            e, ff, g, h = y
            return ((a * e + b * g) % 3, (a * ff + b * h) % 3,
                    (c * e + d * g) % 3, (c * ff + d * h) % 3)
        out[1] = ("SL(2,3)", table(E, mul, (1, 0, 0, 1))[0])
        # every C_3 : H with H of order 8 -- these supply the remaining types
        Hs = [("C8", cyclic(8)), ("C4xC2", direct(cyclic(4), cyclic(2))),
              ("C2^3", elem_c2(3)), ("D_4", dihedral(4)), ("Q_8", dicyclic(2))]
        for hn, H in Hs:
            for si, s in enumerate(homs_to_pm1(H)):
                out.append(("C3:%s[%d]" % (hn, si), semidirect_C3(H, s)))
        return out
    if n == 26:
        return [("D_13", dihedral(13))]
    if n == 27:
        return [("C9:C3", metacyclic(9, 3, 4)), ("Heis_27", heisenberg27())]
    if n == 28:
        return [("D_14", dihedral(14)), ("Dic_7", dicyclic(7))]
    if n == 30:
        return [("D_15", dihedral(15)), ("C3xD_5", direct(cyclic(3), dihedral(5))),
                ("C5xS_3", direct(cyclic(5), S3))]
    return []


# ---------------------------------------------------------------------------
def fmt_class(cls, nm):
    return "{" + ",".join(nm[x] for x in cls) + "}"


def main():
    print("=" * 78)
    print("verify.py -- Ta's determinacy problem for good involutions of conjugation")
    print("            quandles fails at order 32, and, GRANTING the published A000001")
    print("            group counts, at no smaller order")
    print("=" * 78)

    # -----------------------------------------------------------------------
    print("\n--- PART 2: definition-only control, no theory used at all ---")
    for nm, T, want in (("S_3", S3, 1), ("D_4", dihedral(4), 16),
                        ("Q_8", dicyclic(2), 16), ("D_5", dihedral(5), 1)):
        I = inv_table(T)
        c = sum(1 for rho in all_involutions(len(T)) if raw_is_good(T, I, rho))
        check(c == want, "definition-only census %s" % nm,
              "scanned every self-inverse permutation of the %d-element set: %d good"
              % (len(T), c))
    ident_good = raw_is_good(dihedral(4), inv_table(dihedral(4)), tuple(range(8)))
    check(ident_good, "identity counted",
          "the identity permutation of D_4 satisfies the definition, so |Good(Conj D_4)| = 16 "
          "includes it; this pins the convention against the value 16 in Appendix B of "
          "arXiv:2505.08090v5")

    # -----------------------------------------------------------------------
    print("\n--- PART 3: the two order-32 witnesses ---")
    G1, nm1 = dihedral_named(16)
    G2, nm2 = hol_c8()

    PAPER_G1 = ["{1}", "{r^8}", "{r,r^15}", "{r^2,r^14}", "{r^3,r^13}", "{r^4,r^12}",
                "{r^5,r^11}", "{r^6,r^10}", "{r^7,r^9}",
                "{s,r^2s,r^4s,r^6s,r^8s,r^10s,r^12s,r^14s}",
                "{rs,r^3s,r^5s,r^7s,r^9s,r^11s,r^13s,r^15s}"]
    PAPER_G2 = ["{1}", "{a^4}", "{c,a^4c}", "{a^2,a^6}", "{a^2c,a^6c}",
                "{b,a^2b,a^4b,a^6b}", "{bc,a^2bc,a^4bc,a^6bc}", "{a,a^3,a^5,a^7}",
                "{ac,a^3c,a^5c,a^7c}", "{ab,a^3b,a^5b,a^7b}", "{abc,a^3bc,a^5bc,a^7bc}"]

    R32 = []
    for label, T, nm, paper, want_p, want_good in (
            ("G_1 = D_32 = <r,s | r^16, s^2, srs=r^-1>", G1, nm1, PAPER_G1, 4, 128),
            ("G_2 = Hol(C_8) = <a,b,c | a^8,b^2,c^2, bab^-1=a^-1, cac^-1=a^5>", G2, nm2,
             PAPER_G2, 1, 1024)):
        short = label.split("=")[0].strip()
        check(is_group(T), "group axioms %s" % short,
              "identity, all inverses and FULL associativity on all %d^3 = %d triples"
              % (len(T), len(T) ** 3))
        check(not is_abelian(T), "nonabelian %s" % short,
              "so Conj G is a nontrivial quandle and the Problem applies")
        r = good_count(T, label)
        R32.append(r)
        check((r["n"], r["Z"], r["k"]) == (32, 2, 11), "invariant triple %s" % short,
              "(n, |Z(G)|, k(G)) = (%d, %d, %d)" % (r["n"], r["Z"], r["k"]))
        got = [fmt_class(c, nm) for c in r["cls"]]
        check(sorted(got) == sorted(paper), "class list %s" % short,
              "all %d conjugacy classes agree, element by element, with the list printed "
              "in Section 3" % r["k"])
        check(r["Zelts"] == sorted(r["Zelts"]) and len(r["Zelts"]) == 2
              and set(nm[z] for z in r["Zelts"]) == ({"1", "r^8"} if short == "G_1"
                                                     else {"1", "a^4"}),
              "centre %s" % short, "Z(G) = {%s}" % ", ".join(nm[z] for z in r["Zelts"]))
        check(r["invclosed"], "classes inverse-closed %s" % short,
              "every one of the %d classes satisfies C = C^{-1}" % r["k"])
        check(r["lemma"], "reduction lemma %s" % short,
              "the three set-level equivalences of Lemma 2.1 hold by exhaustion over all "
              "x, y in G and all pairs of central values; no published theorem is assumed")
        check(r["p"] == want_p and r["fixed"] == r["k"] - 2 * want_p,
              "central-involution action %s" % short,
              "multiplication by t on Cl(G): %d fixed classes and %d two-cycles, "
              "%d + 2*%d = %d = k" % (r["fixed"], r["p"], r["fixed"], r["p"], r["k"]))
        check(r["closed"] == want_good, "closed form %s" % short,
              "2^(k-p) = 2^(%d-%d) = %d" % (r["k"], r["p"], r["closed"]))
        check(r["good"] == want_good, "matching DP %s" % short,
              "|Good(Conj G)| = %d" % r["good"])
        check(r["enum"] == want_good, "constraint enumeration %s" % short,
              "%d labellings z: Cl(G) -> Z(G) satisfy the Lemma 2.1 constraint" % r["enum"])
        check(r["raw"] == want_good, "FULL raw brute force %s" % short,
              "all %d functions z: Cl(G) -> Z(G) tested by building rho and applying the "
              "definition literally: %s good" % (r["bound"], r["raw"]))
        check(r["rawok"] and r["rawchecked"] == want_good,
              "every counted rho re-checked %s" % short,
              "%d of %d counted involutions satisfy the raw definition"
              % (r["rawchecked"], r["good"]))

    a, b = R32
    check((a["n"], a["Z"], a["k"]) == (b["n"], b["Z"], b["k"]),
          "order-32 pair shares the triple",
          "both at (n, |Z(G)|, k(G)) = (32, 2, 11)")
    check(a["good"] != b["good"], "order-32 pair separates |Good|",
          "%d vs %d" % (a["good"], b["good"]))
    check(fingerprint(G1) != fingerprint(G2), "order-32 pair non-isomorphic",
          "isomorphism invariants already differ, independently of the counts")

    # -----------------------------------------------------------------------
    print("\n--- PART 4: the two order-48 witnesses ---")
    H1 = dihedral(24)
    H2 = direct(S3, dihedral(4))
    R48 = []
    for label, T, want_p, want_good in (
            ("H_1 = D_48 = <r,s | r^24, s^2, srs=r^-1>", H1, 6, 512),
            ("H_2 = S_3 x D_4 (D_4 of order 8)", H2, 3, 4096)):
        short = label.split("=")[0].strip()
        check(is_group(T), "group axioms %s" % short,
              "identity, all inverses and FULL associativity on all %d^3 = %d triples"
              % (len(T), len(T) ** 3))
        r = good_count(T, label, raw_cap=40000, rawcheck_cap=5000)
        R48.append(r)
        check((r["n"], r["Z"], r["k"]) == (48, 2, 15) and not is_abelian(T),
              "invariant triple %s" % short,
              "nonabelian, (n, |Z(G)|, k(G)) = (%d, %d, %d)" % (r["n"], r["Z"], r["k"]))
        check(r["invclosed"] and r["lemma"], "structure %s" % short,
              "all %d classes inverse-closed and Lemma 2.1 verified by exhaustion" % r["k"])
        check(r["p"] == want_p, "central-involution action %s" % short,
              "%d fixed classes, %d two-cycles" % (r["fixed"], r["p"]))
        check(r["closed"] == want_good == r["good"] == r["enum"],
              "three counts agree %s" % short,
              "closed form 2^(%d-%d) = %d, matching DP = %d, constraint enumeration = %d"
              % (r["k"], r["p"], r["closed"], r["good"], r["enum"]))
        check(r["raw"] == want_good, "FULL raw brute force %s" % short,
              "all %d functions z tested against the literal definition: %s"
              % (r["bound"], r["raw"]))
    check((R48[0]["n"], R48[0]["Z"], R48[0]["k"]) == (R48[1]["n"], R48[1]["Z"], R48[1]["k"])
          and R48[0]["good"] != R48[1]["good"], "order-48 pair refutes too",
          "(48,2,15): %d vs %d" % (R48[0]["good"], R48[1]["good"]))

    # -----------------------------------------------------------------------
    print("\n--- PART 5: minimality GRANTING the published A000001 counts"
          " -- every order 6..31 ---")
    census = {}
    settled = {}
    for n in range(6, 32):
        want = A000001[n] - n_abelian(n)
        if want <= 1:
            settled[n] = True
            check(True, "order %d vacuous" % n,
                  "A000001(%d) = %d groups, %d abelian, so %d nonabelian type%s and no "
                  "pair to compare" % (n, A000001[n], n_abelian(n), want,
                                       "" if want == 1 else "s"))
            continue
        P = pool(n)
        bad = [nmx for nmx, T in P if not is_group(T) or len(T) != n]
        check(not bad, "order %d pool is groups of order %d" % (n, n),
              "%d constructions, each with full associativity checked on all %d triples"
              % (len(P), n ** 3))
        seen = {}
        for nmx, T in P:
            if is_abelian(T):
                continue
            fp = fingerprint(T)
            if fp not in seen:
                seen[fp] = (nmx, T)
        check(len(seen) == want, "order %d census exhaustive GRANTING A000001" % n,
              "%d pairwise non-isomorphic nonabelian groups exhibited (distinct "
              "isomorphism invariants) and A000001(%d) - (abelian types) = %d - %d = %d, "
              "so, GRANTING that published count (not verified here), the census is complete"
              % (len(seen), n, A000001[n], n_abelian(n), want))
        rows = []
        for nmx, T in sorted(seen.values(), key=lambda kv: kv[0]):
            r = good_count(T, nmx, raw_cap=3000, rawcheck_cap=1200)
            check(r["lemma"], "order %d lemma %s" % (n, nmx),
                  "Lemma 2.1 verified by exhaustion on this group")
            check(r["good"] == r["enum"], "order %d two counts agree %s" % (n, nmx),
                  "matching DP = %d, constraint enumeration = %d" % (r["good"], r["enum"]))
            if r["raw"] is not None:
                check(r["raw"] == r["good"], "order %d raw brute force %s" % (n, nmx),
                      "all %d functions z tested literally: %d" % (r["bound"], r["raw"]))
            if r["rawchecked"]:
                check(r["rawok"], "order %d counted rho literal %s" % (n, nmx),
                      "%d of %d counted involutions satisfy the raw definition"
                      % (r["rawchecked"], r["good"]))
            rows.append(r)
        census[n] = rows
        buckets = {}
        for r in rows:
            buckets.setdefault((r["Z"], r["k"]), []).append(r)
        live = {key: v for key, v in buckets.items() if len(v) > 1}
        detail = "; ".join("(%d,%d,%d): %s all at %d" % (n, key[0], key[1],
                                                         "+".join(x["name"] for x in v),
                                                         v[0]["good"])
                           for key, v in sorted(live.items()))
        settled[n] = all(len({x["good"] for x in v}) == 1 for v in live.values())
        check(settled[n],
              "order %d determinacy" % n,
              "%d isomorphism types, %d bucket%s of the invariant triple hold two or more "
              "groups%s" % (len(rows), len(live), "" if len(live) == 1 else "s",
                            (" -- " + detail) if detail else ""))
    TA_SUM = {}
    for n in sorted(TA_TABLE):
        if n in census:
            got = sorted(r["good"] for r in census[n])
        else:
            P = [(nmx, T) for nmx, T in pool(n) if not is_abelian(T)]
            got = sorted(good_count(T, nmx, raw_cap=3000, rawcheck_cap=1200)["good"]
                         for nmx, T in P)
        check(got == sorted(TA_TABLE[n]), "external cross-check order %d" % n,
              "our |Good(Conj G)| multiset %s equals the table of Appendix B of "
              "arXiv:2505.08090v5" % (got,))
        TA_SUM[n] = sum(got)
    check(TA_SUM == A386233, "external cross-check OEIS A386233",
          "per-order sums %s equal the ten published terms of A386233"
          % ([TA_SUM[n] for n in sorted(TA_SUM)],))

    # the Remark after Corollary 2.2: no |Z|^{orbits} rule at order 27
    r27 = good_count(metacyclic(9, 3, 4), "C9:C3", raw_cap=0, rawcheck_cap=0)
    k27 = r27["k"]
    seen27 = [False] * k27
    orb = 0
    for c in range(k27):
        if seen27[c]:
            continue
        orb += 1
        stack = [c]
        seen27[c] = True
        while stack:
            u = stack.pop()
            for row in r27["iota"]:
                v = row[u]
                if not seen27[v]:
                    seen27[v] = True
                    stack.append(v)
    check(orb == 5 and r27["good"] == 324 and 3 ** orb == 243,
          "no |Z|^orbits rule beyond |Z|=2",
          "at order 27 the group generated by inversion and by multiplication by Z(G) has "
          "%d orbits on the %d classes, so an |Z|^orbits rule would give 3^%d = %d, whereas "
          "|Good(Conj G)| = %d" % (orb, k27, orb, 3 ** orb, r27["good"]))

    # -----------------------------------------------------------------------
    print("\n--- PART 6: the refutation, assembled ---")
    zero = [n for n in (23, 25, 29, 31) if A000001[n] - n_abelian(n) == 0]
    check(zero == [23, 25, 29, 31] and A000001[26] - n_abelian(26) == 1,
          "orders 23, 25, 26, 29, 31 carry no pair",
          "A000001 minus the abelian types gives 0 nonabelian groups at 23, 25, 29, 31 "
          "(three primes and 5^2) and exactly 1 at 26, so no two groups can share a triple")
    check(a["good"] == 128 and b["good"] == 1024,
          "Problem 11.6, first question: NO",
          "D_32 and Hol(C_8) both have (n,|Z(G)|,k(G)) = (32,2,11) but |Good(Conj G)| = "
          "128 and 1024")
    check(True, "Problem 11.6, second question: MOOT",
          "it asks for a formula in the triple CONDITIONAL on the first answer being yes; "
          "with the antecedent refuted the conditional is vacuous, not false")
    check(sorted(settled) == list(range(6, 32)) and all(settled.values()),
          "GRANTING the published A000001 counts, 32 is minimal",
          "all %d orders 6..31 settled: each is either vacuous (at most one nonabelian "
          "type) or carries a census, complete only by the quoted A000001 term, in which "
          "every invariant-triple bucket is "
          "constant; and orders 1..5 have no nonabelian group at all"
          % len(settled))

    # -----------------------------------------------------------------------
    print("\n--- SUMMARY ---")
    print("%-52s %4s %4s %4s %10s" % ("group", "n", "|Z|", "k", "|Good|"))
    for r in R32 + R48:
        print("%-52s %4d %4d %4d %10d" % (r["name"][:52], r["n"], r["Z"], r["k"], r["good"]))
    print("")
    for n in sorted(census):
        for r in census[n]:
            print("  order %2d  %-14s |Z|=%d k=%2d |Good|=%d"
                  % (n, r["name"], r["Z"], r["k"], r["good"]))

    print("\nNOT RE-RUN: the isomorphism-type counts A000001(n) are taken from OEIS as "
          "external input and are not proved here; the program only exhibits enough "
          "pairwise non-isomorphic groups to meet them.")
    print("NOT RE-RUN: the third bucket reported in the source row, (32,4,14) with 21600 "
          "against 552960, is not re-derived here and nothing above depends on it.")
    print("NOT RE-RUN: the search that FOUND the witnesses -- a sweep of 2450 nonabelian "
          "groups of orders 6..48 -- is not reproduced; this program verifies exhibited "
          "objects only.")
    print("NOT RE-RUN: orders 33..47 and every order above 48 are not examined, so nothing "
          "here claims 32 is the ONLY order at which the triple fails.")
    print("NOT RE-RUN: the definition-only scan over ALL self-inverse permutations of G is "
          "done at orders 6, 8, 10 only; above that the count relies on Lemma 2.1, which "
          "is itself verified by exhaustion on every group used.")
    print("NOT RE-RUN: the table of Appendix B of arXiv:2505.08090v5 is used as an "
          "external cross-check for orders 6..22; it is never an input to any claim "
          "above.")

    if _FAIL:
        print("\nVERDICT: %d OF %d CHECKS FAILED: %s" % (len(_FAIL), _N, ", ".join(_FAIL)))
        return 1
    print("\nVERDICT: ALL %d CHECKS PASS" % _N)
    return 0


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    sys.exit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- re-derives every computational claim of

    "An interior rank-one coefficient refutes the generalized Kaplansky-Lvov
     conjecture of Panja, Saini and Singh"   (paper.tex, this folder)

from the object PRINTED IN THE PAPER and from nothing else.

The object is five 2x2 matrices with entries in {0,1}.  They are not hard-coded as Python
literals: they are parsed below out of PAPER_DISPLAY, which is the verbatim text of the
display in Section 2 of paper.tex.  Everything else the paper asserts -- the exact image of
x A y, the finite-field cardinalities 10 / 50 / 33 / 339, the rank-r family, and all five
positive controls -- is recomputed here.

CONTRACT
    Python 3.9+, STANDARD LIBRARY ONLY (no numpy, no sympy, no external data file).
    Exact arithmetic throughout: integers, fractions.Fraction over Q, and Z/pZ for the
    finite-field work.  NO FLOAT EVER DECIDES ANYTHING.  No randomness anywhere.
    One `PASS <name> -- <detail>` line per check; exits 0 iff every check passed.

WHAT THIS PROGRAM DOES NOT COVER is printed as `NOT RE-RUN:` lines at the end of the run,
and repeated in REVIEW_NOTE.md's `## Scope`.
"""

import re
import sys
from fractions import Fraction

# ----------------------------------------------------------------------------------------
# THE OBJECT, verbatim from the display in Section 2 of paper.tex.
# ----------------------------------------------------------------------------------------
PAPER_DISPLAY = """
A   = [[1, 0], [0, 0]]
X1  = [[1, 0], [0, 1]]
Y1  = [[1, 0], [0, 1]]
X2  = [[0, 0], [1, 0]]
Y2  = [[0, 1], [0, 0]]
"""

CHECKS = []
FAILED = []


def chk(name, ok, detail=""):
    """Record one check.  `name` is a single token so the transcript stays machine-readable."""
    CHECKS.append((name, bool(ok), detail))
    if ok:
        print("PASS %s -- %s" % (name, detail))
    else:
        FAILED.append(name)
        print("FAIL %s -- %s" % (name, detail))


# ----------------------------------------------------------------------------------------
# parsing
# ----------------------------------------------------------------------------------------
def parse_display(text):
    """Read `NAME = [[a, b], [c, d]]` lines into {name: tuple-of-tuples}.  No eval()."""
    out = {}
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        m = re.match(r"^([A-Za-z]\w*)\s*=\s*\[\s*(\[.*\])\s*\]$", line)
        if not m:
            raise ValueError("unparseable display line: %r" % line)
        rows = re.findall(r"\[([^\[\]]*)\]", m.group(2))
        mat = tuple(tuple(int(e.strip()) for e in r.split(",")) for r in rows)
        widths = set(len(r) for r in mat)
        if len(widths) != 1 or len(mat) != widths.pop():
            raise ValueError("display matrix %s is not square" % m.group(1))
        out[m.group(1)] = mat
    return out


# ----------------------------------------------------------------------------------------
# exact linear algebra
# ----------------------------------------------------------------------------------------
def mmulZ(P, Q):
    """matrix product over Z, hence over any commutative ring, hence over any field."""
    inner, wide = len(Q), len(Q[0])
    return tuple(tuple(sum(P[i][t] * Q[t][j] for t in range(inner)) for j in range(wide))
                 for i in range(len(P)))


def maddZ(P, Q):
    return tuple(tuple(P[i][j] + Q[i][j] for j in range(len(P[0]))) for i in range(len(P)))


def rank_Q(M):
    """exact rank over the rationals: Gaussian elimination in Fraction, no floats."""
    A = [[Fraction(e) for e in row] for row in M]
    rows, cols, r = len(A), len(A[0]), 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = Fraction(1) / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def rank_Fp(M, p):
    """exact rank over Z/pZ, p prime."""
    A = [[e % p for e in row] for row in M]
    rows, cols, r = len(A), len(A[0]), 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def _tuples(k, m, alphabet=None):
    """all k-tuples over range(m), or over `alphabet` when given (len(alphabet) == m)."""
    alpha = list(range(m)) if alphabet is None else list(alphabet)
    out = [()]
    for _ in range(k):
        out = [t + (a,) for t in out for a in alpha]
    return out


def outer(u, v, q=None):
    """the outer product u v^T of a column u and a row v."""
    if q is None:
        return tuple(tuple(u[i] * v[j] for j in range(len(v))) for i in range(len(u)))
    return tuple(tuple((u[i] * v[j]) % q for j in range(len(v))) for i in range(len(u)))


def basis_mats(n):
    """the n^2 matrix units e_ij over Z, indexed i*n + j."""
    return [tuple(tuple(1 if (a, b) == (i, j) else 0 for b in range(n)) for a in range(n))
            for i in range(n) for j in range(n)]


def diag(n, ones):
    return tuple(tuple(1 if (i == j and i < ones) else 0 for j in range(n)) for i in range(n))


def rank1_count(n, q):
    """the paper's closed form for |{M in M_n(F_q) : rank M <= 1}|."""
    return 1 + (q ** n - 1) ** 2 // (q - 1)


# ----------------------------------------------------------------------------------------
# M_n(F_q) for prime q.  Matrices are tuples of rows; rows are tuples of ints in [0, q).
# ----------------------------------------------------------------------------------------
class Alg(object):
    def __init__(self, n, q):
        self.n, self.q = n, q
        self.vecs = _tuples(n, q)
        self.mats = _tuples(n, q ** n, alphabet=self.vecs)
        self.matset = set(self.mats)
        self.eye = diag(n, n)
        self._rc = {}

    def vmul(self, v, Y):
        """row vector times matrix, memoised on (v, Y).  This is what lets the all-pairs
        enumerations below finish in seconds with no numpy: there are only q^n * q^(n^2)
        distinct (row, matrix) pairs, and each is computed once."""
        key = (v, Y)
        r = self._rc.get(key)
        if r is None:
            q, n = self.q, self.n
            r = tuple(sum(v[t] * Y[t][j] for t in range(n)) % q for j in range(n))
            self._rc[key] = r
        return r

    def mul(self, X, Y):
        return tuple(self.vmul(row, Y) for row in X)

    def add(self, X, Y):
        q, n = self.q, self.n
        return tuple(tuple((X[i][j] + Y[i][j]) % q for j in range(n)) for i in range(n))

    def sub(self, X, Y):
        q, n = self.q, self.n
        return tuple(tuple((X[i][j] - Y[i][j]) % q for j in range(n)) for i in range(n))

    def smul(self, c, X):
        q, n = self.q, self.n
        return tuple(tuple((c * X[i][j]) % q for j in range(n)) for i in range(n))

    def red(self, M):
        q, n = self.q, self.n
        return tuple(tuple(M[i][j] % q for j in range(n)) for i in range(n))

    def zero(self):
        return tuple((0,) * self.n for _ in range(self.n))

    def rank_le_locus(self, r):
        return set(M for M in self.mats if rank_Fp(M, self.q) <= r)

    def is_subspace(self, S):
        """S closed under addition and under every scalar multiple, and containing 0."""
        if self.zero() not in S:
            return False
        for X in S:
            for c in range(self.q):
                if self.smul(c, X) not in S:
                    return False
        for X in S:
            for Y in S:
                if self.add(X, Y) not in S:
                    return False
        return True


# ========================================================================================
print("=" * 88)
print("PART 1 -- THE WITNESS PRINTED IN THE PAPER (parsed out of PAPER_DISPLAY, not")
print("          hard-coded), over Z and hence over every field")
print("=" * 88)

D = parse_display(PAPER_DISPLAY)
A = D["A"]
I2 = ((1, 0), (0, 1))
E11 = ((1, 0), (0, 0))
E22 = ((0, 0), (0, 1))
E21 = ((0, 0), (1, 0))
E12 = ((0, 1), (0, 0))

chk("display-parsed", sorted(D) == ["A", "X1", "X2", "Y1", "Y2"],
    "PAPER_DISPLAY yields exactly A, X1, Y1, X2, Y2, each 2x2")
chk("A-is-e11", A == E11,
    "A = %s = e_11, the INTERIOR coefficient of omega(x_1,x_2) = x_1 A x_2" % (A,))
chk("A-nonzero", any(e for row in A for e in row),
    "A != 0, so the pair (omega, A) is inside the conjecture's quantifier")
chk("entries-0-1", all(e in (0, 1) for M in D.values() for row in M for e in row),
    "every printed entry is 0 or 1, so the witness transports to EVERY field along Z -> K")
chk("X1-is-I", D["X1"] == I2 and D["Y1"] == I2, "first evaluation point is (X,Y) = (I_2, I_2)")
chk("X2-is-e21", D["X2"] == E21 and D["Y2"] == E12,
    "second evaluation point is (X,Y) = (e_21, e_12)")

w1 = mmulZ(mmulZ(D["X1"], A), D["Y1"])
w2 = mmulZ(mmulZ(D["X2"], A), D["Y2"])
s = maddZ(w1, w2)

chk("omega-I-I", w1 == E11, "omega(I_2, I_2) = I A I = %s = e_11" % (w1,))
chk("omega-e21-e12", w2 == E22, "omega(e_21, e_12) = e_21 e_11 e_12 = %s = e_22" % (w2,))
chk("sum-is-identity", s == I2, "e_11 + e_22 = %s = I_2" % (s,))
chk("rankQ-A", rank_Q(A) == 1, "rank A = 1 over Q")
chk("rankQ-w1", rank_Q(w1) == 1, "rank omega(I,I) = 1 over Q")
chk("rankQ-w2", rank_Q(w2) == 1, "rank omega(e_21,e_12) = 1 over Q")
chk("rankQ-sum", rank_Q(s) == 2, "rank(e_11 + e_22) = rank I_2 = 2 over Q")
chk("det-sum", s[0][0] * s[1][1] - s[0][1] * s[1][0] == 1,
    "det(e_11 + e_22) = 1 is nonzero in every field, so the sum has rank 2 in every field")
for p in (2, 3, 5, 7, 11):
    chk("ranks-mod-%d" % p,
        (rank_Fp(A, p), rank_Fp(w1, p), rank_Fp(w2, p), rank_Fp(s, p)) == (1, 1, 1, 2),
        "over GF(%d): rank A = rank omega(I,I) = rank omega(e_21,e_12) = 1, rank(sum) = 2" % p)
chk("witness-strict", rank_Q(s) > rank_Q(A),
    "rank(sum) = %d > %d = rank A -- this inequality is the entire refutation"
    % (rank_Q(s), rank_Q(A)))

# ========================================================================================
print()
print("=" * 88)
print("PART 2 -- THE FACTORIZATION IDENTITY  X e_11 Y = (col_1 X)(row_1 Y)")
print("=" * 88)

for n in (2, 3, 4):
    B = basis_mats(n)
    e11n = B[0]
    bad = sum(1 for X in B for Y in B
              if mmulZ(mmulZ(X, e11n), Y) != outer(tuple(X[i][0] for i in range(n)), Y[0]))
    chk("identity-basis-n%d" % n, bad == 0,
        "X e_11 Y == (col_1 X)(row_1 Y) on all %d ordered pairs of matrix units of M_%d, "
        "over Z" % (len(B) ** 2, n))

for (n, q) in ((2, 2), (2, 3), (3, 2)):
    alg = Alg(n, q)
    e11n = basis_mats(n)[0]
    bad = 0
    for X in alg.mats:
        XA = alg.mul(X, e11n)
        col1 = tuple(X[i][0] for i in range(n))
        for Y in alg.mats:
            if alg.mul(XA, Y) != outer(col1, Y[0], q):
                bad += 1
    chk("identity-allpairs-n%d-q%d" % (n, q), bad == 0,
        "identity verified on ALL %d ordered pairs (X,Y) in M_%d(GF(%d))^2, 0 violations"
        % (len(alg.mats) ** 2, n, q))

# additivity and homogeneity in each variable: "multilinear" checked, not assumed
for (n, q, probe_only) in ((2, 2, False), (2, 3, True)):
    alg = Alg(n, q)
    e11n = basis_mats(n)[0]
    probe = [alg.red(m) for m in basis_mats(n)] + [alg.eye] if probe_only else alg.mats
    badL = badR = badS = 0
    for Y in probe:
        for X in alg.mats:
            for Xp in alg.mats:
                if alg.mul(alg.mul(alg.add(X, Xp), e11n), Y) != alg.add(
                        alg.mul(alg.mul(X, e11n), Y), alg.mul(alg.mul(Xp, e11n), Y)):
                    badL += 1
            for c in range(q):
                if alg.mul(alg.mul(alg.smul(c, X), e11n), Y) != alg.smul(
                        c, alg.mul(alg.mul(X, e11n), Y)):
                    badS += 1
    for X in probe:
        for Y in alg.mats:
            for Yp in alg.mats:
                if alg.mul(alg.mul(X, e11n), alg.add(Y, Yp)) != alg.add(
                        alg.mul(alg.mul(X, e11n), Y), alg.mul(alg.mul(X, e11n), Yp)):
                    badR += 1
    chk("multilinear-n%d-q%d" % (n, q), badL == 0 and badR == 0 and badS == 0,
        "omega = x e_11 y is additive in x (%d violations), additive in y (%d) and "
        "homogeneous in x (%d); the free variable ran over all %d of M_%d(GF(%d)), the "
        "spectator over %s" % (badL, badR, badS, len(alg.mats), n, q,
                               "all of it" if not probe_only
                               else "the %d matrix units and I" % (n * n)))

# ========================================================================================
print()
print("=" * 88)
print("PART 3 -- THE IMAGE IS EXACTLY THE RANK-<=1 LOCUS, AND IT IS NOT A SUBSPACE")
print("=" * 88)

for (n, q, allpairs) in ((2, 2, True), (2, 3, True), (3, 2, False), (3, 3, False)):
    alg = Alg(n, q)
    e11n = basis_mats(n)[0]
    img = set(outer(u, v, q) for u in alg.vecs for v in alg.vecs)
    if allpairs:
        img2 = set(alg.mul(alg.mul(X, e11n), Y) for X in alg.mats for Y in alg.mats)
        chk("image-allpairs-agrees-n%d-q%d" % (n, q), img == img2,
            "the image over ALL %d pairs (X,Y) equals the image over the %d outer products"
            % (len(alg.mats) ** 2, len(alg.vecs) ** 2))
    locus = alg.rank_le_locus(1)
    chk("image-is-rank1-locus-n%d-q%d" % (n, q), img == locus,
        "image(x e_11 y) = {M in M_%d(GF(%d)) : rank M <= 1} exactly; both sides have %d "
        "elements out of %d" % (n, q, len(img), len(alg.mats)))
    chk("image-size-n%d-q%d" % (n, q), len(img) == rank1_count(n, q),
        "|image| = %d = 1 + (q^n - 1)^2/(q - 1) at n = %d, q = %d" % (len(img), n, q))
    a11, a22 = alg.red(basis_mats(n)[0]), alg.red(basis_mats(n)[n + 1])
    chk("summands-in-image-n%d-q%d" % (n, q), a11 in img and a22 in img,
        "e_11 = omega(I,I) and e_22 = omega(e_21,e_12) both lie in the image over GF(%d)" % q)
    chk("sum-outside-image-n%d-q%d" % (n, q), alg.add(a11, a22) not in img,
        "e_11 + e_22 has rank 2 and is NOT in the image: the image is not closed under "
        "addition over GF(%d)" % q)
    chk("image-scalar-closed-n%d-q%d" % (n, q),
        alg.zero() in img and all(alg.smul(c, M) in img for M in img for c in range(q)),
        "the image DOES contain 0 and IS closed under scalar multiplication -- it is a cone, "
        "and only additivity fails")
    chk("image-not-subspace-n%d-q%d" % (n, q), not alg.is_subspace(img),
        "VERDICT for M_%d(GF(%d)): the image of the multilinear omega(x,y) = x e_11 y is NOT "
        "a vector space" % (n, q))

# ========================================================================================
print()
print("=" * 88)
print("PART 4 -- THE FAMILY: AN INTERIOR A OF RANK r WITH 0 < r < n GIVES THE RANK-<=r LOCUS")
print("=" * 88)

for (n, q, r) in ((3, 2, 1), (3, 2, 2), (3, 3, 1), (4, 2, 1), (4, 2, 2)):
    alg = Alg(n, q)
    Dr = diag(n, r)
    chk("rankA-n%d-q%d-r%d" % (n, q, r), rank_Fp(Dr, q) == r,
        "A = diag(1^%d, 0^%d) has rank %d over GF(%d), and 0 < %d < %d"
        % (r, n - r, r, q, r, n))
    Us = _tuples(n, q ** r, alphabet=_tuples(r, q))
    Vs = _tuples(r, q ** n, alphabet=_tuples(n, q))
    img = set()
    for U in Us:
        for V in Vs:
            img.add(tuple(tuple(sum(U[i][t] * V[t][j] for t in range(r)) % q
                                for j in range(n)) for i in range(n)))
    chk("family-image-n%d-q%d-r%d" % (n, q, r), img == alg.rank_le_locus(r),
        "image(x A y) = {M : rank M <= %d} in M_%d(GF(%d)), %d of %d matrices, computed from "
        "all %d factorizations X A Y = (X's first %d columns)(Y's first %d rows)"
        % (r, n, q, len(img), len(alg.mats), len(Us) * len(Vs), r, r))
    M1, M2 = Dr, alg.red(basis_mats(n)[r * n + r])
    chk("family-sum-escapes-n%d-q%d-r%d" % (n, q, r),
        M1 in img and M2 in img and alg.add(M1, M2) not in img
        and rank_Fp(alg.add(M1, M2), q) == r + 1,
        "diag(1^%d,0^%d) and e_%d%d lie in the image; their sum has rank %d > %d and does not"
        % (r, n - r, r + 1, r + 1, r + 1, r))

# ========================================================================================
print()
print("=" * 88)
print("PART 5 -- POSITIVE CONTROLS: THE NEIGHBOURING CASES WHERE THE IMAGE *IS* A VECTOR")
print("          SPACE, so the attack is not breaking everything it touches")
print("=" * 88)

# (a) A invertible
for (n, q) in ((2, 2), (2, 3), (3, 2), (3, 3)):
    alg = Alg(n, q)
    img = set(alg.mul(alg.mul(alg.eye, alg.eye), M) for M in alg.mats)
    ok = (img == alg.matset)
    if len(alg.mats) <= 600:
        ok = ok and alg.is_subspace(img)
    chk("control-A-invertible-n%d-q%d" % (n, q), ok,
        "A = I_%d invertible: image(x A y) is all %d of M_%d(GF(%d)), a vector space%s"
        % (n, len(alg.mats), n, q,
           "" if len(alg.mats) > 600 else " (closure re-checked elementwise)"))

# (b) m = 1
for (n, q) in ((2, 2), (2, 3), (3, 2)):
    alg = Alg(n, q)
    e11n = alg.red(basis_mats(n)[0])
    img = set(alg.mul(e11n, M) for M in alg.mats)
    chk("control-m1-n%d-q%d" % (n, q), alg.is_subspace(img) and len(img) == q ** n,
        "omega(x) = e_11 x with m = 1: image = {M : rows 2..%d vanish}, %d elements, IS a "
        "subspace" % (n, len(img)))

# (c) coefficient at an END
for (n, q) in ((2, 2), (2, 3), (3, 2)):
    alg = Alg(n, q)
    e11n = alg.red(basis_mats(n)[0])
    reduced = set(alg.mul(e11n, M) for M in alg.mats)
    if n == 2:
        pairs = set(alg.mul(alg.mul(e11n, X), Y) for X in alg.mats for Y in alg.mats)
        same = (pairs == reduced)
    else:
        same = all(alg.mul(alg.mul(e11n, X), alg.eye) in reduced for X in alg.mats)
    chk("control-end-coefficient-n%d-q%d" % (n, q), same and alg.is_subspace(reduced),
        "omega(x,y) = e_11 x y with the coefficient at an END: image = {e_11 M} = "
        "{M : col M contained in col A}, %d elements, IS a subspace" % len(reduced))

# (d) the generalized-commutator shape the original row's key named: Y = I collapses it
for (n, q) in ((2, 2), (2, 3), (3, 2), (3, 3)):
    alg = Alg(n, q)
    B = [alg.red(m) for m in basis_mats(n)]
    tested = [(B[0], B[n + 1]), (B[0], B[0]), (alg.eye, B[0]), (B[1], B[n])]
    ok = True
    for (Am, Bm) in tested:
        target = alg.sub(Am, Bm)
        for X in alg.mats:
            lhs = alg.sub(alg.mul(alg.mul(Am, X), alg.eye), alg.mul(alg.mul(Bm, alg.eye), X))
            if lhs != alg.mul(target, X):
                ok = False
                break
    chk("control-YeqI-collapse-n%d-q%d" % (n, q), ok,
        "for omega(x,y) = A x y - B y x: omega(X, I) = (A - B)X for all %d X and all %d "
        "tested (A,B) over GF(%d), so that image always CONTAINS the subspace (A-B)M_%d -- "
        "which is why no rank obstruction can reach that shape"
        % (len(alg.mats), len(tested), q, n))

# (e)+(f) the commutator set, computed once per (n,q) and used twice
for (n, q, sl_size) in ((2, 2, 8), (2, 3, 27), (3, 2, 256)):
    alg = Alg(n, q)
    comm = set()
    for X in alg.mats:
        for Y in alg.mats:
            comm.add(alg.sub(alg.mul(X, Y), alg.mul(Y, X)))
    trace0 = set(M for M in alg.mats if sum(M[i][i] for i in range(n)) % q == 0)
    chk("control-commutators-sl%d-q%d" % (n, q),
        len(comm) == sl_size and comm == trace0 and alg.is_subspace(comm),
        "image of xy - yx on M_%d(GF(%d)) over all %d pairs is exactly {tr = 0}: %d = "
        "%d^(%d^2-1) elements, a subspace, matching the published Albert-Muckenhoupt count"
        % (n, q, len(alg.mats) ** 2, len(comm), q, n))
    e11n = alg.red(basis_mats(n)[0])
    capped = set(alg.mul(e11n, C) for C in comm)
    chk("control-e11-commutator-n%d-q%d" % (n, q),
        alg.is_subspace(capped) and all(rank_Fp(M, q) <= 1 for M in capped),
        "A = B = e_11 on M_%d(GF(%d)): the image of e_11(xy - yx) has %d elements, EVERY one "
        "of rank <= 1, and it IS a subspace -- a rank CEILING alone refutes nothing; cell 1 "
        "fails because its confinement is a rank condition, not a linear one"
        % (n, q, len(capped)))

# ========================================================================================
print()
print("=" * 88)
print("SCOPE -- WHAT THIS PROGRAM DOES NOT COVER")
print("=" * 88)
print("NOT RE-RUN: the PUBLISHED wording of the conjecture.  Every locator is the e-print "
      "arXiv:2312.13865v1 (conc.tex, lines 21-23).  The Wiley text of Mathematika 71 (2025) "
      "e70031 returned HTTP 403 and was never read, so whether the published sentence "
      "restricts the coefficient positions is settled neither here nor in the paper.")
print("NOT RE-RUN: the generalized-commutator cell of the original row -- omega(x,y) = "
      "A x y - B y x over an algebraically closed field of characteristic 0.  Parts 5(d) and "
      "5(e) show only why a rank obstruction cannot reach that shape; nothing here decides "
      "it, and the paper claims nothing about it.")
print("NOT RE-RUN: any infinite, characteristic-0 or algebraically closed field.  Every "
      "enumeration runs over GF(2) or GF(3).  The refutation does not depend on them: Part 1 "
      "is an exact computation valid over an arbitrary field, and Parts 2-4 verify the "
      "structural identities that carry it to every field and every n >= 2.")
print("NOT RE-RUN: n >= 5; at n = 4 only q = 2, and only r in {1, 2}.  The rank-r family "
      "claim is PROVED in the paper by rank factorization -- Part 4 is a finite "
      "confirmation, not a proof by exhaustion.")
print("NOT RE-RUN: the two AWS finite-field censuses filed with the original row "
      "(census_kl_job1.py, census_kl_job2.py, numpy-based, 78 s and 1982 s).  Their stdout "
      "was never captured to a file and the claim they supported was struck in review as "
      "off-cell.  Their decisive negative-polarity integers -- |image| = 10, 50, 33, 339 -- "
      "are recomputed from scratch in Part 3 here, with no numpy and no floats.")

print()
n_pass = sum(1 for _, ok, _ in CHECKS if ok)
if FAILED:
    print("FAILURES: %s" % ", ".join(FAILED))
    print("VERDICT: %d of %d CHECKS FAILED" % (len(FAILED), len(CHECKS)))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % n_pass)
sys.exit(0)

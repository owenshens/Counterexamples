#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- referee verification program for

    "A Length-14 Negative Answer to Alderson's q=4 Optimality Question"
    (Theorem 1 and its Corollary; refutes optimality in Alderson's Problem 9.2)

Standard library only, exact finite-field arithmetic, no floating point.

======================================================================
VALUES TAKEN FROM THE PAPER  (inputs -- these are NOT checks)
======================================================================
  * The ambient objects: E = F_4^4, U = F_2^4 <= E, Pi = PG(E) = PG(3,4),
    and the standard Baer subgeometry B = { <u>_{F_4} : 0 != u in U }.
  * The matrix
        T = [[0,0,1,0],[1,0,1,0],[0,1,0,0],[0,0,0,1]]
    asserted by the paper to lie in GL(4,2) and to have order 7.
  * The four vectors
        u  = (1,0,0,w),    v  = (0,1,w^2,0),
        u' = (1,0,0,w^2),  v' = (0,1,w,0),
    with L = <u,v>_{F_4} and L' = <u',v'>_{F_4}, w^2+w+1 = 0.
  * The indexing of the fourteen lines
        l_{j+1} = PG(T^j L),  l_{j+8} = PG(T^j L'),  0 <= j <= 6.
  * The definition C = { (rho_1(a),...,rho_14(a)) : a in E } for F_4-linear
    epimorphisms rho_i : E -> F_16 with ker rho_i = V_i, l_i = PG(V_i).
  * The paper's asserted determinant tuples and weight enumerator, quoted
    verbatim so that the recomputed values can be compared against them:
        (det[u,v,T^k u,T^k v])_{k=1..6}      = (w,1,w^2,w^2,1,w)
        (det[u',v',T^k u',T^k v'])_{k=1..6}  = (w^2,1,w,w,1,w^2)
        (det[u,v,T^k u',T^k v'])_{k=0..6}    = (1,w,w,w^2,w,w^2,w^2)
        W_C(z) = 1 + 210 z^13 + 45 z^14
  * The comparison length 112 from Alderson's Theorem 6.6 / Problem 9.2.

======================================================================
DERIVED HERE  (what the program actually computes and therefore checks)
======================================================================
  * F_4 and F_16 are built from scratch (polynomial arithmetic) and their
    field axioms are verified; F_16 = F_4[x]/(x^2+x+w) is checked to be a
    field of order 16 containing F_4.
  * T's invertibility, its exact multiplicative order, and the minimal
    polynomial X^3+X+1 of its upper-left 3x3 block.
  * All 19 determinants are recomputed over F_4 from T, u, v, u', v' and
    compared with the paper's tuples; nonvanishing is checked.
  * The fourteen lines are built as explicit 5-point subsets of the 85
    points of PG(3,4): pairwise disjointness (91 pairs), disjointness from
    B, the 70-point union, and the mixed partition Pi = B + l_1 + ... + l_14.
  * All 357 lines of PG(3,4) are enumerated; none lies inside B.
  * rho_1..rho_14 are constructed by linear algebra (not quoted), C is
    built, |C| = 256, C is an F_4-subspace of F_4-dimension 4, C is
    faithful, the zero-coordinate lemma holds for all 255 nonzero a, and
    the weight enumerator, minimum distance and MDS/Singleton equality are
    computed.
  * The 16-ary extension C^+ is constructed and d(C^+) = 14 is computed.
  * MAXIMALITY (load-bearing, computed not asserted): all 4^8 = 65536
    F_4-linear maps rho : E -> F_16 are enumerated and for each one the
    minimum distance of { (c(a), rho(a)) } is COMPUTED from the recomputed
    weights of C (no distance is assigned from a literal); none reaches 14.
    Completeness of that enumeration is itself certified rather than counted:
    the 65536 basis-image 4-tuples (rho(e_1),...,rho(e_4)) are checked to be
    all of F_16^4, which is exactly the statement that the enumeration is the
    whole of Hom_{F_4}(E, F_16).
  * C is checked NOT to be F_16-linear.  This is an independent consistency
    test of the refutation: an F_16-linear [14,2,13] MDS code extends to a
    [15,2,14] linear -- hence F_4-additive -- code.
  * B is checked to have no 2-secant, the paper's other assertion about B.
  * Choice-independence: the construction is repeated with alternative
    complements defining rho_i, and with four genuinely different F_2-linear
    isomorphisms tau : E/U -> F_16 given by invertible 4x4 F_2 matrices (each
    verified F_2-linear, onto and with kernel exactly U), the minimum distance
    of the resulting C^+ being recomputed over all pairs each time.  Merely
    permuting the coset coordinates would NOT be an independent choice of tau:
    the only feature of tau used downstream is its zero set, and every
    permutation leaves that equal to U.
  * Mutation sensitivity: the whole battery is re-run on seven deliberately
    corrupted versions of the exhibited object and is required to report
    failures for each, and the maximality search is run on a positive control
    (one 2-dimensional subspace artificially removed from the weight-13 set)
    where it must find exactly |GL(2,4)| = 180 distance-14 extensions.  A
    search that could only ever answer "no extension" would fail that check.
  * The paper's own positive control, on a hole set that DOES contain lines:
    a regular spread of PG(3,4) is built here (not quoted) and checked to be
    seventeen pairwise disjoint lines partitioning Pi; fourteen of them replace
    l_1,...,l_14, leaving the three remaining lines as holes; the resulting
    code is checked to have the same enumerator 1 + 210 z^13 + 45 z^14; and the
    SAME exhaustive search over all 65536 maps rho must then ACCEPT exactly
    3|GL(2,4)| = 540 of them -- the expected count being derived by enumerating
    ordered bases of each hole's annihilator, not quoted -- and the 540
    accepted maps are exhibited as, for each hole line PG(V), the 180 maps with
    kernel exactly V.
  * The primary battery is run inside an exception barrier, so a corrupt
    exhibited object is reported as FAILures with a verdict line rather than
    aborting the program with a traceback and no verdict.

Caveats a referee should keep in mind.  (i) Faithfulness is checked in the
sense the paper itself uses -- every coordinate projection is onto F_16 -- and
that property is automatic for the epimorphisms this program constructs; if
Alderson's notion of "faithful" is stronger, that is not tested here.  (ii)
The statements of Alderson's Problem 9.2 and of "extendable" / "additively
maximal", and the length 112, are taken from the cited paper and are not
re-derived; what is verified is that this code has all seven properties and
has length 14 < 112.  (iii) The Remark's identification of Alderson's
scattered F_2-linear set with a Baer subgeometry projectively equivalent to B,
and the classifications of van Dam and Mellinger, are not verified here.
Every one of these caveats is also PRINTED by the program, on the NOT RE-RUN
lines just above the verdict, so that the transcript alone discloses them.

Usage:  python3 verify.py    (exit status 0 iff every check passes)
Runtime: about 10 seconds single-process; no arguments, no input files.
"""

import itertools
import sys

PAPER_DET_L = ("w", "1", "w2", "w2", "1", "w")
PAPER_DET_LP = ("w2", "1", "w", "w", "1", "w2")
PAPER_DET_MIX = ("1", "w", "w", "w2", "w", "w2", "w2")
PAPER_ENUM = {0: 1, 13: 210, 14: 45}
ALDERSON_LENGTH = 112


class Battery(object):
    """Collects (name, ok, detail) triples."""

    def __init__(self):
        self.rows = []

    def add(self, name, ok, detail=""):
        self.rows.append((name, bool(ok), str(detail)))
        return bool(ok)

    def failures(self):
        return [r for r in self.rows if not r[1]]

    def result(self, name):
        """True only if a check of this name ran and passed."""
        for n, ok, _ in self.rows:
            if n == name:
                return ok
        return False


# ----------------------------------------------------------------------
# F_4 = F_2[w]/(w^2+w+1).  Element a = a1*w + a0 is stored as 2*a1 + a0,
# so 0,1,2=w,3=w^2=w+1 and addition is XOR.
# ----------------------------------------------------------------------
F4 = (0, 1, 2, 3)
OMEGA = 2
OMEGA2 = 3
F4NAME = {0: "0", 1: "1", 2: "w", 3: "w2"}


def _build_f4_mul():
    tab = [[0] * 4 for _ in range(4)]
    for a in F4:
        for b in F4:
            a1, a0 = a >> 1, a & 1
            b1, b0 = b >> 1, b & 1
            c2 = a1 & b1              # coefficient of w^2
            c1 = (a1 & b0) ^ (a0 & b1)
            c0 = a0 & b0
            c1 ^= c2                  # w^2 = w + 1
            c0 ^= c2
            tab[a][b] = (c1 << 1) | c0
    return tab


F4MUL = _build_f4_mul()


def m4(a, b):
    return F4MUL[a][b]


def inv_f4(a):
    for b in (1, 2, 3):
        if F4MUL[a][b] == 1:
            return b
    raise ZeroDivisionError("0 has no inverse in F_4")


def frob4(a):
    """Frobenius x -> x^2 on F_4."""
    return F4MUL[a][a]


# ----------------------------------------------------------------------
# F_16 = F_4[x]/(x^2 + x + w), element (hi, lo) = hi*x + lo with hi,lo in F_4.
# The F_4-subfield is { (0, c) }.  Addition is componentwise in F_4.
# ----------------------------------------------------------------------
F16 = tuple((hi, lo) for hi in F4 for lo in F4)
ZERO16 = (0, 0)
ONE16 = (0, 1)
X16 = (1, 0)


def a16(p, q):
    return (p[0] ^ q[0], p[1] ^ q[1])


def m16(p, q):
    """(p1 x + p0)(q1 x + q0) with x^2 = x + w (char 2)."""
    p1, p0 = p
    q1, q0 = q
    t = m4(p1, q1)
    hi = t ^ m4(p1, q0) ^ m4(p0, q1)
    lo = m4(OMEGA, t) ^ m4(p0, q0)
    return (hi, lo)


def scal16(c, p):
    """Multiply an F_16 element by the F_4-scalar c."""
    return (m4(c, p[0]), m4(c, p[1]))


# ----------------------------------------------------------------------
# Linear algebra over F_4 (4x4 matrices, column vectors).
# ----------------------------------------------------------------------
I4 = [[1 if i == j else 0 for j in range(4)] for i in range(4)]


def matmul(A, B):
    out = []
    for i in range(4):
        row = []
        for j in range(4):
            s = 0
            for k in range(4):
                s ^= m4(A[i][k], B[k][j])
            row.append(s)
        out.append(row)
    return out


def matvec(A, x):
    out = []
    for i in range(4):
        s = 0
        for k in range(4):
            s ^= m4(A[i][k], x[k])
        out.append(s)
    return out


def vadd(x, y):
    return [x[i] ^ y[i] for i in range(4)]


def vscal(c, x):
    return [m4(c, x[i]) for i in range(4)]


def det4(cols):
    """Determinant over F_4 of the matrix whose columns are cols[0..3]."""
    tot = 0
    for p in itertools.permutations(range(4)):
        t = 1
        for i in range(4):
            t = m4(t, cols[p[i]][i])
        tot ^= t
    return tot


def rref(vectors):
    """Row-reduce a list of length-4 F_4 vectors; return (basis, pivots)."""
    rows = [list(v) for v in vectors]
    piv = []
    r = 0
    for c in range(4):
        pr = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                pr = i
                break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        iv = inv_f4(rows[r][c])
        rows[r] = vscal(iv, rows[r])
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = vadd(rows[i], vscal(rows[i][c], rows[r]))
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    return [row for row in rows[:r]], piv


def rank4(vectors):
    return len(rref(vectors)[0])


def _rank_general(vectors):
    """Rank over F_4 of a list of vectors of arbitrary common length."""
    rows = [list(v) for v in vectors]
    n = len(rows[0]) if rows else 0
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                pr = i
                break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        iv = inv_f4(rows[r][c])
        rows[r] = [m4(iv, t) for t in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [rows[i][j] ^ m4(f, rows[r][j]) for j in range(n)]
        r += 1
        if r == len(rows):
            break
    return r


def inv4(A):
    """Inverse of a 4x4 F_4 matrix, or None if singular."""
    M = [A[i][:] + [1 if i == j else 0 for j in range(4)] for i in range(4)]
    r = 0
    for c in range(4):
        pr = None
        for i in range(r, 4):
            if M[i][c]:
                pr = i
                break
        if pr is None:
            return None
        M[r], M[pr] = M[pr], M[r]
        iv = inv_f4(M[r][c])
        M[r] = [m4(iv, t) for t in M[r]]
        for i in range(4):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [M[i][j] ^ m4(f, M[r][j]) for j in range(8)]
        r += 1
    return [row[4:] for row in M]


# ----------------------------------------------------------------------
# PG(3,4): points are F_4^4 \ {0} normalized at the first nonzero coordinate.
# ----------------------------------------------------------------------
ALLVEC = [list(t) for t in itertools.product(F4, repeat=4)]
NONZERO = [x for x in ALLVEC if any(x)]


def pt(x):
    """Projective normalization of a nonzero vector."""
    for i in range(4):
        if x[i]:
            iv = inv_f4(x[i])
            return tuple(m4(iv, c) for c in x)
    raise ValueError("zero vector has no projective point")


POINTS = sorted(set(pt(x) for x in NONZERO))
BAER = sorted(set(pt(list(t)) for t in itertools.product((0, 1), repeat=4)
                  if any(t)))
BAERSET = set(BAER)


def span2(a, b):
    """All 16 vectors of the F_4-span of a, b."""
    out = set()
    for s in F4:
        for t in F4:
            out.add(tuple(vadd(vscal(s, a), vscal(t, b))))
    return out


def line_of(a, b):
    """The projective point set of <a,b>_{F_4}; 5 points iff a,b independent."""
    return frozenset(pt(list(w)) for w in span2(a, b) if any(w))


def all_lines_of_pg34():
    """Every line of PG(3,4), built from pairs of distinct points."""
    out = set()
    for i in range(len(POINTS)):
        for j in range(i + 1, len(POINTS)):
            out.add(line_of(list(POINTS[i]), list(POINTS[j])))
    return sorted(out, key=lambda s: sorted(s))


ALL_LINES = all_lines_of_pg34()


# ----------------------------------------------------------------------
# The exhibited object of the paper.
# ----------------------------------------------------------------------
PAPER_T = [[0, 0, 1, 0], [1, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
PAPER_U = [1, 0, 0, OMEGA]
PAPER_V = [0, 1, OMEGA2, 0]
PAPER_UP = [1, 0, 0, OMEGA2]
PAPER_VP = [0, 1, OMEGA, 0]


class Spread(object):
    """Decoded exhibited object: T, its powers, the two seed planes, the
    fourteen 2-dimensional subspaces V_i and the fourteen lines l_i."""

    def __init__(self, T, u, v, up, vp):
        self.T = [row[:] for row in T]
        self.u, self.v, self.up, self.vp = u[:], v[:], up[:], vp[:]
        self.pw = [I4]
        for _ in range(7):
            self.pw.append(matmul(self.pw[-1], self.T))
        self.gens = []
        for j in range(7):
            self.gens.append((matvec(self.pw[j], self.u),
                              matvec(self.pw[j], self.v)))
        for j in range(7):
            self.gens.append((matvec(self.pw[j], self.up),
                              matvec(self.pw[j], self.vp)))
        self.V = []          # spanning pairs, reduced
        self.lines = []      # 5-point sets (or fewer, if the input is corrupt)
        for (a, b) in self.gens:
            self.V.append(rref([a, b])[0])
            self.lines.append(line_of(a, b))

    def line_index_of_point(self, p):
        return [i for i in range(len(self.lines)) if p in self.lines[i]]


def build_rho(V, complements):
    """An F_4-linear epimorphism E -> F_16 with kernel exactly span(V).

    V is a list of vectors spanning a 2-dimensional subspace; `complements`
    is an ordered list of candidate vectors used to extend V to a basis of E.
    Returns a callable, or None when span(V) is not 2-dimensional.
    """
    bas = rref(V)[0]
    if len(bas) != 2:
        return None
    full = [r[:] for r in bas]
    for e in complements:
        if len(full) == 4:
            break
        if rank4(full + [list(e)]) == len(full) + 1:
            full.append(list(e))
    if len(full) != 4:
        return None
    A = [[full[k][i] for k in range(4)] for i in range(4)]   # columns = basis
    Ainv = inv4(A)
    if Ainv is None:
        return None
    imgs = (ZERO16, ZERO16, ONE16, X16)

    def rho(x):
        y = matvec(Ainv, x)
        out = ZERO16
        for k in range(4):
            if y[k]:
                out = a16(out, scal16(y[k], imgs[k]))
        return out

    return rho


UNIT_VECTORS = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


class AdditiveCode(object):
    """C = { (rho_1(a),...,rho_14(a)) : a in E } together with its extension."""

    def __init__(self, spread, complements=None):
        comp = UNIT_VECTORS if complements is None else complements
        self.rhos = [build_rho(V, comp) for V in spread.V]
        self.ok = all(r is not None for r in self.rhos)
        self.words = {}
        if not self.ok:
            return
        for a in ALLVEC:
            self.words[tuple(a)] = tuple(r(a) for r in self.rhos)

    def word(self, a):
        return self.words[tuple(a)]

    def weight(self, a):
        return sum(1 for s in self.words[tuple(a)] if s != ZERO16)


# ----------------------------------------------------------------------
# The 16-ary extension.  Over F_2, E = U + wU, so a = a0 + w*a1 with
# a0,a1 in F_2^4 and a + U <-> a1.  tau : F_2^4 -> F_16 is the F_2-linear
# isomorphism a1 |-> (a1[2] + w a1[3]) x + (a1[0] + w a1[1]).
# ----------------------------------------------------------------------
def coset_coords(a):
    """The F_2^4 coordinate vector of a + U in E/U."""
    return [c >> 1 for c in a]


def tau(a, perm=(0, 1, 2, 3)):
    c = coset_coords(a)
    c = [c[perm[0]], c[perm[1]], c[perm[2]], c[perm[3]]]
    return (c[2] | (c[3] << 1), c[0] | (c[1] << 1))


# Genuinely different F_2-linear isomorphisms E/U -> F_16: apply an invertible
# 4x4 matrix over F_2 to the coset coordinates.  (A mere permutation of the
# coset coordinates is NOT an independent choice of tau, because the only
# feature of tau used downstream is its zero set, and every permutation
# leaves that equal to U.)
ALT_TAU_MATRICES = (
    ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
    ((1, 1, 0, 0), (0, 1, 1, 0), (0, 0, 1, 1), (0, 0, 0, 1)),
    ((0, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0)),
    ((0, 0, 0, 1), (1, 0, 0, 1), (0, 1, 0, 0), (0, 0, 1, 0)),
)


def tau_mat(a, M):
    """The F_16 element whose F_2-coordinate vector is M * (a + U)."""
    c = coset_coords(a)
    cp = [0, 0, 0, 0]
    for i in range(4):
        s = 0
        for j in range(4):
            s ^= M[i][j] & c[j]
        cp[i] = s
    return (cp[2] | (cp[3] << 1), cp[0] | (cp[1] << 1))


def in_U(a):
    return all(c >> 1 == 0 for c in a)


def dotf4(m, a):
    s = 0
    for i in range(4):
        s ^= m4(m[i], a[i])
    return s


def exhaustive_maximality(weight_groups):
    """Enumerate ALL 4^8 = 65536 F_4-linear maps rho : F_4^4 -> F_16 and, for
    each, COMPUTE the minimum distance of the length-15 F_4-additive code
    { (c(a), rho(a)) : a in E },

        d(rho) = min over nonzero a of  wt(c(a)) + [rho(a) != 0].

    A map is parametrised by a pair (m1, m2) of F_4-linear functionals, rho(a)
    = (m1.a) x + (m2.a) in the F_4-basis {x, 1} of F_16; every F_4-linear map
    arises exactly once, which the caller certifies from the returned set of
    basis images.

    `weight_groups` is a list of (weight, vectors) pairs, sorted by increasing
    weight, covering every nonzero a of E.  Nothing is assigned: the returned
    dmax is a maximum of computed minima.

    Returns (count_ge_14, examined, dmax, basis_images).
    """
    rows = [list(t) for t in itertools.product(F4, repeat=4)]
    count14 = 0
    examined = 0
    dmax = 0
    basis_images = set()
    for m1 in rows:
        pre = [[dotf4(m1, a) for a in vecs] for _, vecs in weight_groups]
        for m2 in rows:
            examined += 1
            # rho(e_k) = (m1[k], m2[k]) in the (hi, lo) coordinates of F_16.
            basis_images.add(((m1[0], m2[0]), (m1[1], m2[1]),
                              (m1[2], m2[2]), (m1[3], m2[3])))
            best = None
            for gi in range(len(weight_groups)):
                wt = weight_groups[gi][0]
                if best is not None and best <= wt:
                    break                      # later groups are heavier
                vecs = weight_groups[gi][1]
                d1 = pre[gi]
                in_ker = False
                out_ker = False
                for idx in range(len(vecs)):
                    if d1[idx] == 0 and dotf4(m2, vecs[idx]) == 0:
                        in_ker = True
                    else:
                        out_ker = True
                    if in_ker and out_ker:
                        break
                cand = wt if in_ker else (wt + 1 if out_ker else None)
                if cand is not None and (best is None or cand < best):
                    best = cand
            if best is None:
                continue
            if best > dmax:
                dmax = best
            if best >= 14:
                count14 += 1
    return count14, examined, dmax, basis_images


# ----------------------------------------------------------------------
# Battery A -- the arithmetic and geometric models this program builds.
# ----------------------------------------------------------------------
def checks_model(bat):
    ok = True
    for a in F4:
        for b in F4:
            if m4(a, b) != m4(b, a):
                ok = False
            for c in F4:
                if m4(m4(a, b), c) != m4(a, m4(b, c)):
                    ok = False
                if m4(a, b ^ c) != (m4(a, b) ^ m4(a, c)):
                    ok = False
    units = all(any(m4(a, b) == 1 for b in F4) for a in (1, 2, 3))
    bat.add("f4_model_is_a_field_of_order_4", ok and units and m4(0, 1) == 0,
            "|F_4|=%d" % len(set(F4)))
    bat.add("f4_omega_satisfies_X2_plus_X_plus_1",
            (m4(OMEGA, OMEGA) ^ OMEGA ^ 1) == 0 and OMEGA2 == m4(OMEGA, OMEGA),
            "w^2=%s" % F4NAME[m4(OMEGA, OMEGA)])

    ok = True
    for p in F16:
        for q in F16:
            if m16(p, q) != m16(q, p):
                ok = False
            for r in F16:
                if m16(m16(p, q), r) != m16(p, m16(q, r)):
                    ok = False
                if m16(p, a16(q, r)) != a16(m16(p, q), m16(p, r)):
                    ok = False
    inverses = all(any(m16(p, q) == ONE16 for q in F16)
                   for p in F16 if p != ZERO16)
    e = ONE16
    order_x = None
    for t in range(1, 17):
        e = m16(e, X16)
        if e == ONE16:
            order_x = t
            break
    bat.add("f16_model_is_a_field_of_order_16", ok and inverses
            and len(set(F16)) == 16, "|F_16|=%d" % len(set(F16)))
    bat.add("f16_multiplicative_order_of_x_is_15", order_x == 15,
            "ord(x)=%s" % order_x)
    sub = set((0, c) for c in F4)
    closed = all(m16(p, q) in sub and a16(p, q) in sub for p in sub
                 for q in sub)
    bat.add("f16_contains_f4_as_the_subfield_of_scalars", closed and len(sub) == 4,
            "|subfield|=%d" % len(sub))
    bat.add("pg_3_4_has_85_points", len(POINTS) == 85, "%d points" % len(POINTS))
    bat.add("pg_3_4_has_357_lines", len(ALL_LINES) == 357,
            "%d lines" % len(ALL_LINES))
    bat.add("baer_subgeometry_B_has_15_points", len(BAER) == 15,
            "|B|=%d" % len(BAER))


# ----------------------------------------------------------------------
# Battery B -- the exhibited object is well formed and is what the paper says.
# ----------------------------------------------------------------------
def checks_object(bat, sp):
    entries_f2 = all(c in (0, 1) for row in sp.T for c in row)
    dT = det4([[sp.T[i][j] for i in range(4)] for j in range(4)])
    bat.add("T_has_F2_entries_and_is_invertible", entries_f2 and dT != 0,
            "det T = %s" % F4NAME[dT])
    ordT = None
    for k in range(1, 8):
        if sp.pw[k] == I4:
            ordT = k
            break
    bat.add("T_has_multiplicative_order_exactly_7", ordT == 7,
            "order = %s" % ordT)

    M = [[sp.T[i][j] for j in range(3)] for i in range(3)]
    M2 = [[0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            s = 0
            for k in range(3):
                s ^= m4(M[i][k], M[k][j])
            M2[i][j] = s
    M3 = [[0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            s = 0
            for k in range(3):
                s ^= m4(M2[i][k], M[k][j])
            M3[i][j] = s
    idn = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
    vanishes = all(M3[i][j] ^ M[i][j] ^ idn[i][j] == 0
                   for i in range(3) for j in range(3))
    flat = [[idn[i][j] for i in range(3) for j in range(3)],
            [M[i][j] for i in range(3) for j in range(3)],
            [M2[i][j] for i in range(3) for j in range(3)]]
    indep = _rank_general(flat) == 3
    bat.add("T_upper_left_block_has_min_poly_X3_plus_X_plus_1",
            vanishes and indep,
            "M^3+M+I=0 : %s, deg 3 : %s" % (vanishes, indep))

    bat.add("L_and_Lprime_are_2dimensional_F4_subspaces",
            rank4([sp.u, sp.v]) == 2 and rank4([sp.up, sp.vp]) == 2,
            "dim L = %d, dim L' = %d" % (rank4([sp.u, sp.v]),
                                         rank4([sp.up, sp.vp])))
    fro = set(tuple(frob4(c) for c in w) for w in span2(sp.u, sp.v))
    bat.add("Lprime_is_the_Frobenius_conjugate_of_L",
            fro == set(tuple(w) for w in span2(sp.up, sp.vp)),
            "|phi(L)|=%d" % len(fro))

    for tag, x, y in (("L", sp.u, sp.v), ("Lprime", sp.up, sp.vp)):
        reps = [x, y, vadd(x, y), vadd(x, vscal(OMEGA, y)),
                vadd(x, vscal(OMEGA2, y))]
        if rank4([x, y]) != 2 or not all(any(c for c in r) for r in reps):
            # A corrupt object can make one of the five representatives zero;
            # that is a FAILURE, not a reason to abort the report.
            for nm in ("PG_%s_five_points_are_the_paper_representatives" % tag,
                       "PG_%s_representatives_each_have_a_coordinate_outside_F2"
                       % tag):
                bat.add(nm, False,
                        "the span of the two generators is not 2-dimensional, "
                        "or a representative is the zero vector")
            continue
        normalized = all(any(c for c in r) and r[[i for i in range(4)
                         if r[i]][0]] == 1 for r in reps)
        pts = set(pt(r) for r in reps)
        bat.add("PG_%s_five_points_are_the_paper_representatives" % tag,
                len(pts) == 5 and pts == set(line_of(x, y)) and normalized,
                "%d distinct, normalized=%s" % (len(pts), normalized))
        outside = all(any(c not in (0, 1) for c in r) for r in reps)
        bat.add("PG_%s_representatives_each_have_a_coordinate_outside_F2" % tag,
                outside and not (pts & BAERSET),
                "outside_F2=%s, meets B in %d" % (outside, len(pts & BAERSET)))


# ----------------------------------------------------------------------
# Battery C -- the 19 determinants of the paper's "direct calculation".
# ----------------------------------------------------------------------
def checks_determinants(bat, sp):
    d1 = [det4([sp.u, sp.v, matvec(sp.pw[k], sp.u), matvec(sp.pw[k], sp.v)])
          for k in range(1, 7)]
    d2 = [det4([sp.up, sp.vp, matvec(sp.pw[k], sp.up), matvec(sp.pw[k], sp.vp)])
          for k in range(1, 7)]
    d3 = [det4([sp.u, sp.v, matvec(sp.pw[k], sp.up), matvec(sp.pw[k], sp.vp)])
          for k in range(0, 7)]
    n1 = tuple(F4NAME[t] for t in d1)
    n2 = tuple(F4NAME[t] for t in d2)
    n3 = tuple(F4NAME[t] for t in d3)
    bat.add("det_u_v_Tku_Tkv_k1to6_matches_paper", n1 == PAPER_DET_L,
            "computed %s, paper %s" % (",".join(n1), ",".join(PAPER_DET_L)))
    bat.add("det_up_vp_Tkup_Tkvp_k1to6_matches_paper", n2 == PAPER_DET_LP,
            "computed %s, paper %s" % (",".join(n2), ",".join(PAPER_DET_LP)))
    bat.add("det_u_v_Tkup_Tkvp_k0to6_matches_paper", n3 == PAPER_DET_MIX,
            "computed %s, paper %s" % (",".join(n3), ",".join(PAPER_DET_MIX)))
    alld = d1 + d2 + d3
    bat.add("all_19_determinants_are_nonzero",
            len(alld) == 19 and all(t != 0 for t in alld),
            "%d values, %d zero" % (len(alld), sum(1 for t in alld if t == 0)))


# ----------------------------------------------------------------------
# Battery D -- the fourteen lines and the mixed partition of PG(3,4).
# ----------------------------------------------------------------------
def checks_geometry(bat, sp):
    nl = len(sp.lines)
    sizes = [len(l) for l in sp.lines]
    bat.add("fourteen_lines_each_carry_exactly_5_points",
            nl == 14 and sizes == [5] * nl,
            "%d lines, sizes %s" % (nl, sizes))
    real = nl == 14 and all(l in set(ALL_LINES) for l in sp.lines)
    bat.add("fourteen_lines_are_genuine_lines_of_PG_3_4", real,
            "all %d occur in the 357-line enumeration: %s" % (nl, real))
    meets = [(i + 1, j + 1) for i in range(nl) for j in range(i + 1, nl)
             if sp.lines[i] & sp.lines[j]]
    npairs = nl * (nl - 1) // 2
    bat.add("fourteen_lines_are_pairwise_disjoint", nl == 14 and not meets,
            "%d of %d pairs meet" % (len(meets), npairs))
    hits = [i + 1 for i in range(nl) if sp.lines[i] & BAERSET]
    bat.add("fourteen_lines_are_disjoint_from_the_Baer_subgeometry", not hits,
            "%d lines meet B" % len(hits))
    union = set()
    for l in sp.lines:
        union |= set(l)
    bat.add("union_of_the_fourteen_lines_has_70_points",
            len(union) == 70 and nl == 14,
            "|union| = %d, %d*5 = %d" % (len(union), nl, nl * 5))
    partition = (union | BAERSET) == set(POINTS) and not (union & BAERSET)
    bat.add("mixed_partition_B_plus_14_lines_equals_PG_3_4",
            partition and len(union) + len(BAER) == len(POINTS),
            "70+15=%d, |PG(3,4)|=%d, covered=%s"
            % (len(union) + len(BAER), len(POINTS), partition))
    counts = [len(sp.line_index_of_point(p)) for p in POINTS]
    every = all(counts[i] == (0 if POINTS[i] in BAERSET else 1)
                for i in range(len(POINTS)))
    bat.add("every_point_outside_B_lies_on_exactly_one_of_the_14_lines", every,
            "multiplicities over the %d points: %s"
            % (len(POINTS), sorted(set(counts))))
    inter = max((len(l & BAERSET) for l in ALL_LINES), default=-1)
    three = sum(1 for l in ALL_LINES if len(l & BAERSET) == 3)
    bat.add("no_line_of_PG_3_4_lies_inside_B",
            inter == 3 and three == 35,
            "max |line cap B| = %d over %d lines; %d are 3-secants"
            % (inter, len(ALL_LINES), three))
    # The paper's other assertion about B: the F_4-line through two distinct
    # points of B meets B in the three points of the corresponding F_2-subline,
    # i.e. B has no 2-secant.
    tally = dict()
    for l in ALL_LINES:
        s = len(l & BAERSET)
        tally[s] = tally.get(s, 0) + 1
    bat.add("B_has_no_2_secant_so_two_points_of_B_span_an_F2_subline",
            tally.get(2, 0) == 0 and tally.get(1, 0) == 210
            and tally.get(0, 0) == 112 and sum(tally.values()) == 357,
            "line-B intersection sizes %s" % sorted(tally.items()))
    img = []
    for i in range(nl):
        tgt = frozenset(pt(matvec(sp.T, list(p))) for p in sp.lines[i])
        cand = [k for k in range(nl) if sp.lines[k] == tgt]
        img.append(cand[0] if len(cand) == 1 else None)
    cyc = (img[:7] == [1, 2, 3, 4, 5, 6, 0]
           and img[7:] == [8, 9, 10, 11, 12, 13, 7])
    bat.add("T_permutes_the_14_lines_as_two_7_cycles", cyc,
            "image indices %s" % [None if t is None else t + 1 for t in img])


# ----------------------------------------------------------------------
# Battery E -- the additive code C.
# ----------------------------------------------------------------------
def checks_code(bat, sp, code):
    if not bat.add("rho_1_to_rho_14_exist_with_kernels_V_1_to_V_14", code.ok):
        return None
    kern_ok = len(code.rhos) == 14
    for i in range(len(code.rhos)):
        ker = [a for a in ALLVEC if code.rhos[i](a) == ZERO16]
        if len(ker) != 16 or set(map(tuple, ker)) != set(
                tuple(w) for w in span2(*sp.gens[i])):
            kern_ok = False
    onto = [len(set(r(a) for a in ALLVEC)) for r in code.rhos]
    bat.add("each_ker_rho_i_equals_V_i_and_rho_i_is_onto_F16",
            kern_ok and onto == [16] * len(code.rhos),
            "kernels agree with V_1..V_14 : %s ; image sizes %s"
            % (kern_ok, sorted(set(onto))))
    cw = set(code.words.values())
    bat.add("C_has_256_words_and_a_to_c_a_is_injective",
            len(cw) == 256 and len(code.words) == 256,
            "|C| = %d = 16^2" % len(cw))
    lin = True
    for a in ALLVEC:
        for c in F4:
            want = tuple(scal16(c, s) for s in code.word(a))
            if code.word(vscal(c, a)) != want:
                lin = False
    add_ok = True
    for a in ALLVEC:
        wa = code.word(a)
        for b in ALLVEC:
            if code.word(vadd(a, b)) != tuple(
                    a16(s, t) for s, t in zip(wa, code.word(b))):
                add_ok = False
                break
        if not add_ok:
            break
    flat = [[c for s in code.word(a) for c in s] for a in ALLVEC]
    dim = _rank_general(flat)
    bat.add("C_is_an_F4_subspace_of_F16_14_of_F4_dimension_4",
            lin and add_ok and dim == 4 and len(cw) == 4 ** dim,
            "F_4-dim = %d, |C| = %d, F_4-linear = %s"
            % (dim, len(cw), lin and add_ok))
    nn = len(next(iter(cw)))
    surj = [len(set(w[i] for w in cw)) for i in range(nn)]
    bat.add("C_is_faithful_every_coordinate_projection_is_onto_F16",
            nn == 14 and surj == [16] * nn,
            "%d coordinates, image sizes %s" % (nn, sorted(set(surj))))
    zc = True
    for a in NONZERO:
        zeros = set(i for i in range(nn) if code.word(a)[i] == ZERO16)
        if zeros != set(sp.line_index_of_point(pt(a))):
            zc = False
    bat.add("zero_coordinates_of_c_a_are_exactly_the_lines_through_a", zc,
            "checked all %d nonzero a" % len(NONZERO))
    return cw


def hdist(x, y):
    return sum(1 for s, t in zip(x, y) if s != t)


def min_distance(words):
    """Minimum Hamming distance over all unordered pairs (no linearity used)."""
    ws = list(words)
    best = len(ws[0])
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            d = hdist(ws[i], ws[j])
            if d < best:
                best = d
    return best


# ----------------------------------------------------------------------
# Battery F -- weight enumerator, minimum distance, MDS.
# ----------------------------------------------------------------------
def checks_parameters(bat, sp, code, cw):
    enum = {}
    for w in cw:
        k = sum(1 for s in w if s != ZERO16)
        enum[k] = enum.get(k, 0) + 1
    bat.add("weight_enumerator_is_1_plus_210z13_plus_45z14",
            enum == PAPER_ENUM,
            "computed %s, paper %s" % (sorted(enum.items()),
                                       sorted(PAPER_ENUM.items())))
    union = set()
    for l in sp.lines:
        union |= set(l)
    bat.add("weight_counts_equal_3x_the_geometric_point_counts",
            enum.get(14, 0) == 3 * len(BAER)
            and enum.get(13, 0) == 3 * len(union),
            "45 = 3*%d ; 210 = 3*%d" % (len(BAER), len(union)))
    d = min_distance(cw)
    bat.add("minimum_distance_of_C_is_13",
            d == 13, "d = %d (over all %d unordered pairs)"
            % (d, len(cw) * (len(cw) - 1) // 2))
    n = len(next(iter(cw)))
    k = None
    for t in range(0, 6):
        if 16 ** t == len(cw):
            k = t
    bat.add("C_is_MDS_singleton_bound_is_attained",
            k == 2 and n == 14 and d == n - k + 1,
            "n=%d, k=log_16|C|=%s, d=%d, n-k+1=%s"
            % (n, k, d, None if k is None else n - k + 1))
    bat.add("length_14_is_forced_by_the_partition_85_minus_15_over_5",
            n == (len(POINTS) - len(BAER)) // 5 == 14,
            "(%d-%d)/5 = %d coordinates" % (len(POINTS), len(BAER),
                                            (len(POINTS) - len(BAER)) // 5))
    bat.add("C_is_a_14_2_13_code_over_F16_additive_over_F4",
            len(cw) == 16 ** 2 and n == 14 and d == 13,
            "(n,k,d) = (%d,%s,%d)" % (n, k, d))
    # C must NOT be F_16-linear: an F_16-linear [14,2,13] MDS code is an arc of
    # 14 points of PG(1,16) and extends to a [15,2,14] linear -- hence
    # F_4-additive -- code, which would contradict additive maximality.  So
    # this is an independent consistency test of the load-bearing refutation.
    f16lin = all(tuple(m16(X16, s) for s in w) in cw for w in cw)
    bat.add("C_is_genuinely_additive_and_not_F16_linear", not f16lin,
            "closed under multiplication by the primitive x of F_16 : %s"
            % f16lin)


# ----------------------------------------------------------------------
# Battery G -- the length-15 extension C^+ (16-ary, F_2-linear).
# ----------------------------------------------------------------------
def checks_extension(bat, sp, code, perm=(0, 1, 2, 3)):
    tvals = [tau(a, perm) for a in ALLVEC]
    tmap = dict((tuple(a), tvals[i]) for i, a in enumerate(ALLVEC))
    f2lin = all(tmap[tuple(vadd(a, b))] == a16(tmap[tuple(a)], tmap[tuple(b)])
                for a in ALLVEC for b in ALLVEC)
    kern = [a for a in ALLVEC if tau(a, perm) == ZERO16]
    tau_ok = (len(set(tvals)) == 16 and f2lin
              and set(map(tuple, kern)) == set(tuple(a) for a in ALLVEC
                                               if in_U(a)))
    bat.add("tau_is_an_F2_linear_isomorphism_E_mod_U_onto_F16", tau_ok,
            "|image| = %d, kernel = U (%d vectors), F_2-linear = %s"
            % (len(set(tvals)), len(kern), f2lin))
    plus = {}
    for a in ALLVEC:
        plus[tuple(a)] = code.word(a) + (tau(a, perm),)
    pw = set(plus.values())
    punct = set(w[:14] for w in pw)
    bat.add("C_plus_has_256_words_of_length_15_and_punctures_back_to_C",
            len(pw) == 256 and all(len(w) == 15 for w in pw)
            and punct == set(code.words.values()),
            "|C^+| = %d, punctured = C : %s"
            % (len(pw), punct == set(code.words.values())))
    dplus = min_distance(pw)
    bat.add("minimum_distance_of_C_plus_is_14_so_C_is_extendable",
            dplus == 14, "d(C^+) = %d (all %d unordered pairs)"
            % (dplus, len(pw) * (len(pw) - 1) // 2))
    add2 = all(tuple(a16(s, t) for s, t in zip(plus[tuple(a)], plus[tuple(b)]))
               in pw for a in ALLVEC for b in ALLVEC)
    w_closed = all(tuple(scal16(OMEGA, s) for s in plus[tuple(a)]) in pw
                   for a in ALLVEC)
    bat.add("C_plus_is_F2_linear_but_not_F4_additive",
            add2 and not w_closed,
            "F_2-closed = %s, closed under w = %s" % (add2, w_closed))
    return pw


# ----------------------------------------------------------------------
# Battery H -- additive maximality.  This is the load-bearing refutation:
# C is extendable over F_16 but admits no F_4-additive extension.
# ----------------------------------------------------------------------
_SUBS_CACHE = []


def two_dim_subspaces():
    """Every 2-dimensional F_4-subspace of E, as a frozenset of 16 vectors."""
    if _SUBS_CACHE:
        return _SUBS_CACHE[0]
    out = set()
    for i in range(len(NONZERO)):
        for j in range(i + 1, len(NONZERO)):
            a, b = NONZERO[i], NONZERO[j]
            if rank4([a, b]) == 2:
                out.add(frozenset(span2(a, b)))
    _SUBS_CACHE.append(sorted(out, key=lambda s: sorted(s)))
    return _SUBS_CACHE[0]


def checks_maximality(bat, sp, code):
    w13 = [a for a in NONZERO if code.weight(a) == 13]
    w14 = [a for a in NONZERO if code.weight(a) == 14]
    cone = set(tuple(vscal(c, a)) for c in (1, 2, 3) for a in NONZERO
               if in_U(a))
    bat.add("weight14_vectors_are_exactly_the_F4_cone_over_B",
            set(map(tuple, w14)) == cone and len(w14) == 45,
            "%d weight-14 vectors, |cone over B| = %d" % (len(w14), len(cone)))
    subs = two_dim_subspaces()
    w13set = set(map(tuple, w13))
    bad = [W for W in subs if not (set(W) & w13set)]
    bat.add("all_357_two_dim_subspaces_contain_a_weight13_vector",
            len(subs) == 357 and not bad,
            "%d subspaces enumerated, %d lie inside the cone over B"
            % (len(subs), len(bad)))
    groups = [(13, w13), (14, w14)]
    count14, examined, dmax, images = exhaustive_maximality(groups)
    bat.add("no_F4_linear_rho_extends_C_to_distance_14", count14 == 0,
            "%d of %d F_4-linear maps E->F_16 reach d >= 14"
            % (count14, examined))
    # Completeness of the enumeration is certified, not merely counted: an
    # F_4-linear map E -> F_16 is determined by, and free on, the images of
    # e_1..e_4, so the enumeration is the whole of Hom_{F_4}(E, F_16) exactly
    # when the set of basis-image 4-tuples is all of F_16^4.
    allimg = set(itertools.product(F16, repeat=4))
    bat.add("all_65536_F4_linear_maps_were_examined",
            examined == 4 ** 8 and len(images) == 4 ** 8 and images == allimg,
            "%d pairs examined, %d distinct maps, = every element of F_16^4 : %s"
            % (examined, len(images), images == allimg))
    bat.add("best_F4_additive_extension_has_distance_13_not_14", dmax == 13,
            "max over all rho of the computed d({(c(a),rho(a))}) = %d" % dmax)
    # Positive control: if one 2-dimensional subspace W were free of
    # weight-13 vectors, the same search must FIND the extensions with
    # ker rho = W -- exactly |GL(2,4)| = 180 of them.  This shows the search
    # is not hard-wired to answer zero.
    if subs:
        W = set(subs[0])
        ctrl = [a for a in w13 if tuple(a) not in W]
        c14, _, cdmax, _ = exhaustive_maximality([(13, ctrl), (14, w14)])
        bat.add("maximality_search_finds_d14_on_a_positive_control",
                c14 == 180 and cdmax == 14,
                "control with one 2-dim subspace removed from S_13 yields "
                "%d maps (|GL(2,4)| = 180) and d = %d" % (c14, cdmax))
    return count14, examined, dmax


# ----------------------------------------------------------------------
# Battery H2 -- the SPREAD positive control the paper describes: a hole set
# that DOES contain lines, on which the same exhaustive search must ACCEPT
# exactly 3|GL(2,4)| = 540 of the 65536 maps rho, namely, for each hole line
# PG(V), the 180 maps with kernel V.
#
# E = F_4^4 is identified with F_16^2 by the F_4-linear bijection
#       a = (a_0,a_1,a_2,a_3)  <->  (X, Y),   X = a_0 x + a_1, Y = a_2 x + a_3.
# The seventeen one-dimensional F_16-subspaces of F_16^2 are F_4-subspaces of
# dimension 2 and partition E: they are a (regular) spread of PG(3,4).  Nothing
# about the spread is quoted; it is built here and its defining properties are
# checked before it is used as a control.
# ----------------------------------------------------------------------
def _e16(X, Y):
    """The vector of E corresponding to the pair (X, Y) in F_16^2."""
    return [X[0], X[1], Y[0], Y[1]]


def regular_spread():
    """The seventeen 2-dimensional F_4-subspaces of a regular spread, each as
    a (basis pair, frozenset of its 16 vectors) couple, in a fixed order."""
    out = []
    for m in F16:
        b = [_e16(ONE16, m), _e16(X16, m16(m, X16))]
        out.append((b, frozenset(span2(b[0], b[1]))))
    b = [_e16(ZERO16, ONE16), _e16(ZERO16, X16)]
    out.append((b, frozenset(span2(b[0], b[1]))))
    return sorted(out, key=lambda t: sorted(t[1]))


class _SubspaceList(object):
    """Stand-in for Spread carrying only the field AdditiveCode reads, .V."""

    def __init__(self, bases):
        self.V = [rref(b)[0] for b in bases]


_CTRL_DEPENDENTS = (
    "control_code_from_14_spread_lines_has_enumerator_1_210_45",
    "spread_positive_control_accepts_540_of_the_65536_maps",
    "the_540_accepted_maps_are_the_180_with_each_hole_as_kernel",
)


def checks_spread_positive_control(bat):
    """Run the paper's spread control.  Every check below is registered
    unconditionally: if the control code cannot be built, its dependents are
    recorded as FAILures rather than silently skipped."""
    sprd = regular_spread()
    dims_ok = all(rank4(b) == 2 for b, _ in sprd)
    disjoint = True
    for i in range(len(sprd)):
        for j in range(i + 1, len(sprd)):
            if len(sprd[i][1] & sprd[j][1]) != 1:      # only the zero vector
                disjoint = False
    covered = set()
    for _, S in sprd:
        covered |= S
    slines = [line_of(b[0], b[1]) for b, _ in sprd]
    allset = set(ALL_LINES)
    genuine = all(len(L) == 5 and L in allset for L in slines)
    ptunion = set()
    for L in slines:
        ptunion |= set(L)
    bat.add("regular_spread_of_PG_3_4_is_17_disjoint_lines_partitioning_Pi",
            len(sprd) == 17 and dims_ok and disjoint and genuine
            and len(covered) == len(NONZERO) + 1
            and len(ptunion) == len(POINTS)
            and sum(len(L) for L in slines) == len(POINTS),
            "%d subspaces, all of F_4-dimension 2: %s, pairwise intersection "
            "trivial: %s, all genuine lines of PG(3,4): %s, %d of %d points "
            "covered with 17*5 = %d incidences"
            % (len(sprd), dims_ok, disjoint, genuine, len(ptunion),
               len(POINTS), sum(len(L) for L in slines)))

    part, holes = sprd[:14], sprd[14:]
    holepts = [line_of(b[0], b[1]) for b, _ in holes]
    hole_disjoint = all(not (holepts[i] & holepts[j])
                        for i in range(len(holepts))
                        for j in range(i + 1, len(holepts)))
    holevecs = set()
    for _, S in holes:
        holevecs |= S
    subs = two_dim_subspaces()
    inside = [W for W in subs if set(W) <= holevecs]
    bat.add("control_hole_set_is_three_disjoint_lines_and_does_contain_lines",
            len(holes) == 3 and hole_disjoint
            and sum(len(L) for L in holepts) == 15
            and len(inside) == 3
            and set(frozenset(W) for W in inside) == set(S for _, S in holes),
            "%d hole lines, pairwise disjoint: %s, %d hole points; %d of the "
            "%d two-dimensional subspaces lie inside the hole set, and they "
            "are exactly the three holes (the paper's Baer hole set admits "
            "none, which is the whole difference)"
            % (len(holes), hole_disjoint, sum(len(L) for L in holepts),
               len(inside), len(subs)))

    ctrl = AdditiveCode(_SubspaceList([b for b, _ in part]))
    if not ctrl.ok:
        for nm in _CTRL_DEPENDENTS:
            bat.add(nm, False,
                    "the control code over the fourteen spread lines could "
                    "not be built, so this check did not run; an absent check "
                    "is recorded here as a FAILURE, not as a silent skip")
        return
    cenum = {}
    for w in set(ctrl.words.values()):
        k = sum(1 for s in w if s != ZERO16)
        cenum[k] = cenum.get(k, 0) + 1
    cw13 = [a for a in NONZERO if ctrl.weight(a) == 13]
    cw14 = [a for a in NONZERO if ctrl.weight(a) == 14]
    holenz = set(t for t in holevecs if any(t))
    w14_is_holes = set(map(tuple, cw14)) == holenz
    bat.add(_CTRL_DEPENDENTS[0],
            cenum == {0: 1, 13: 210, 14: 45} and w14_is_holes,
            "computed enumerator %s; its %d weight-14 vectors are exactly the "
            "%d nonzero vectors of the three hole subspaces: %s"
            % (sorted(cenum.items()), len(cw14), len(holenz), w14_is_holes))

    c540, cexam, cdmax, _ = exhaustive_maximality([(13, cw13), (14, cw14)])
    # The expected count is DERIVED, not quoted: rho(a) = (m1.a)x + (m2.a) has
    # kernel exactly V iff (m1, m2) is an ordered basis of the annihilator
    # V^perp, so the maps with kernel V are counted by enumerating those.
    per_hole = []
    for b, _ in holes:
        V = rref(b)[0]
        vperp = [m for m in ALLVEC if all(dotf4(m, a) == 0 for a in V)]
        n = sum(1 for m1 in vperp for m2 in vperp if rank4([m1, m2]) == 2)
        per_hole.append((V, vperp, n))
    expected = sum(n for _, _, n in per_hole)
    bat.add(_CTRL_DEPENDENTS[1],
            c540 == expected and expected == 3 * 180 and cdmax == 14
            and cexam == 4 ** 8,
            "the same exhaustive search accepts %d of %d maps and attains "
            "d = %d; derived count %s = %d = 3|GL(2,4)|"
            % (c540, cexam, cdmax,
               " + ".join(str(n) for _, _, n in per_hole), expected))

    wt = dict((tuple(a), ctrl.weight(a)) for a in NONZERO)
    exhibited = 0
    kernels_ok = True
    for V, vperp, _ in per_hole:
        vnz = set(t for t in span2(V[0], V[1]) if any(t))
        for m1 in vperp:
            for m2 in vperp:
                if rank4([m1, m2]) != 2:
                    continue
                exhibited += 1
                ker = set()
                dd = None
                for a in NONZERO:
                    zero = (dotf4(m1, a) == 0 and dotf4(m2, a) == 0)
                    if zero:
                        ker.add(tuple(a))
                    cand = wt[tuple(a)] + (0 if zero else 1)
                    if dd is None or cand < dd:
                        dd = cand
                if dd != 14 or ker != vnz:
                    kernels_ok = False
    bat.add(_CTRL_DEPENDENTS[2],
            kernels_ok and exhibited == expected and exhibited == c540,
            "%d maps exhibited, %d per hole line, each with kernel exactly "
            "that hole and recomputed d = 14: %s; that equals the %d the "
            "search accepted, so the accepted set is exactly these"
            % (exhibited, exhibited // len(per_hole), kernels_ok, c540))


# ----------------------------------------------------------------------
# Battery I -- the claims must not depend on the (non-canonical) choice of
# the epimorphisms rho_i or of tau.
# ----------------------------------------------------------------------
def _lcg_order(seed):
    """A deterministic permutation of NONZERO, used as complement candidates."""
    idx = list(range(len(NONZERO)))
    s = seed
    for i in range(len(idx) - 1, 0, -1):
        s = (1103515245 * s + 12345) % (1 << 31)
        j = s % (i + 1)
        idx[i], idx[j] = idx[j], idx[i]
    return [NONZERO[t] for t in idx]


def checks_choice_independence(bat, sp, base_code):
    base_enum = {}
    for w in set(base_code.words.values()):
        k = sum(1 for s in w if s != ZERO16)
        base_enum[k] = base_enum.get(k, 0) + 1
    base_s13 = set(tuple(a) for a in NONZERO if base_code.weight(a) == 13)
    agree, tried = 0, 0
    for seed in (1, 7, 2026, 90210):
        tried += 1
        alt = AdditiveCode(sp, _lcg_order(seed))
        if not alt.ok:
            continue
        cwa = set(alt.words.values())
        enum = {}
        for w in cwa:
            k = sum(1 for s in w if s != ZERO16)
            enum[k] = enum.get(k, 0) + 1
        s13 = set(tuple(a) for a in NONZERO if alt.weight(a) == 13)
        faithful = all(len(set(w[i] for w in cwa)) == 16
                       for i in range(len(next(iter(cwa)))))
        dplus = min(alt.weight(a) + (0 if tau(a) == ZERO16 else 1)
                    for a in NONZERO)
        if (len(cwa) == 256 and enum == base_enum and s13 == base_s13
                and faithful and dplus == 14):
            agree += 1
    bat.add("claims_are_independent_of_the_choice_of_rho_i",
            tried > 0 and agree == tried,
            "%d of %d alternative complement choices reproduce W_C, "
            "faithfulness, S_13 and d(C^+)=14" % (agree, tried))
    # Four genuinely different F_2-linear isomorphisms E/U -> F_16.  For each
    # one the extension is rebuilt and its minimum distance is recomputed over
    # all unordered pairs, so neither linearity nor the weight enumerator is
    # assumed; each tau is also checked to be F_2-linear, onto F_16 and to have
    # kernel exactly U.
    dvals = []
    valid = []
    for M in ALT_TAU_MATRICES:
        tv = dict((tuple(a), tau_mat(a, M)) for a in ALLVEC)
        lin = all(tv[tuple(vadd(a, b))] == a16(tv[tuple(a)], tv[tuple(b)])
                  for a in ALLVEC for b in ALLVEC)
        ker_u = set(k for k, val in tv.items() if val == ZERO16) == set(
            tuple(a) for a in ALLVEC if in_U(a))
        valid.append(lin and ker_u and len(set(tv.values())) == 16)
        dvals.append(min_distance(
            set(base_code.word(a) + (tv[tuple(a)],) for a in ALLVEC)))
    bat.add("d_C_plus_is_14_for_every_tested_choice_of_tau",
            all(valid) and dvals == [14] * len(ALT_TAU_MATRICES),
            "%d of %d tau are F_2-linear isomorphisms with kernel U; "
            "recomputed d(C^+) values %s"
            % (sum(1 for t in valid if t), len(valid), dvals))


# ----------------------------------------------------------------------
# Battery J -- the corollary: Problem 9.2's q=4 clause is answered
# negatively at length 14 < 112.
# ----------------------------------------------------------------------
REQUIRED = (
    "mixed_partition_B_plus_14_lines_equals_PG_3_4",
    "C_is_a_14_2_13_code_over_F16_additive_over_F4",
    "C_is_an_F4_subspace_of_F16_14_of_F4_dimension_4",
    "C_is_faithful_every_coordinate_projection_is_onto_F16",
    "C_is_MDS_singleton_bound_is_attained",
    "minimum_distance_of_C_plus_is_14_so_C_is_extendable",
    "no_F4_linear_rho_extends_C_to_distance_14",
)


def checks_corollary(bat, length):
    """`length` is the recomputed length of the exhibited code, not a
    constant from the paper."""
    have = [bat.result(nm) for nm in REQUIRED]
    bat.add("problem_9_2_hypotheses_all_verified_for_this_code", all(have),
            "%d of %d required properties hold: %s"
            % (sum(1 for h in have if h), len(have),
               ", ".join(nm for nm, h in zip(REQUIRED, have) if not h)
               or "none missing"))
    bat.add("negative_answer_length_14_is_below_Alderson_112",
            all(have) and length == 14 and length < ALDERSON_LENGTH,
            "exhibited length %d < %d, so 112 is not optimal"
            % (length, ALDERSON_LENGTH))


# ----------------------------------------------------------------------
# Driver.
# ----------------------------------------------------------------------
def run_battery(T, u, v, up, vp, with_model=True, with_extras=True, bat=None):
    if bat is None:
        bat = Battery()
    if with_model:
        checks_model(bat)
    sp = Spread(T, u, v, up, vp)
    checks_object(bat, sp)
    checks_determinants(bat, sp)
    checks_geometry(bat, sp)
    code = AdditiveCode(sp)
    cw = checks_code(bat, sp, code)
    if cw is None:
        return bat, sp, None
    checks_parameters(bat, sp, code, cw)
    checks_extension(bat, sp, code)
    checks_maximality(bat, sp, code)
    if with_extras:
        checks_choice_independence(bat, sp, code)
    checks_corollary(bat, len(next(iter(cw))))
    return bat, sp, code


# ----------------------------------------------------------------------
# Battery K -- mutation sensitivity.  A corrupted exhibited object must make
# the battery above report failures; otherwise its checks are vacuous.
# ----------------------------------------------------------------------
def mutants():
    T2 = [row[:] for row in PAPER_T]
    T2[0][2] = 0                                   # destroys invertibility
    T3 = [[PAPER_T[j][i] for j in range(4)] for i in range(4)]   # transpose
    T4 = [row[:] for row in PAPER_T]
    T4[3][3] = OMEGA                               # leaves GL(4,2)
    return [
        ("T_entry_flipped", T2, PAPER_U, PAPER_V, PAPER_UP, PAPER_VP),
        ("T_transposed", T3, PAPER_U, PAPER_V, PAPER_UP, PAPER_VP),
        ("T_entry_outside_F2", T4, PAPER_U, PAPER_V, PAPER_UP, PAPER_VP),
        ("u_moved_into_F2_4", PAPER_T, [1, 0, 0, 1], PAPER_V, PAPER_UP,
         PAPER_VP),
        ("u_last_coordinate_w2_not_w", PAPER_T, [1, 0, 0, OMEGA2], PAPER_V,
         PAPER_UP, PAPER_VP),
        ("v_third_coordinate_w_not_w2", PAPER_T, PAPER_U, [0, 1, OMEGA, 0],
         PAPER_UP, PAPER_VP),
        ("Lprime_replaced_by_L", PAPER_T, PAPER_U, PAPER_V, PAPER_U, PAPER_V),
    ]


def checks_mutation_sensitivity(bat):
    for label, T, u, v, up, vp in mutants():
        try:
            mb, _, _ = run_battery(T, u, v, up, vp, with_model=False,
                                   with_extras=False)
            fails = mb.failures()
            detail = "%d of %d checks fail; first = %s" % (
                len(fails), len(mb.rows), fails[0][0] if fails else "-")
            ok = len(fails) > 0
        except Exception as exc:                  # a corrupt object may abort
            ok = True
            detail = "battery aborted: %s: %s" % (type(exc).__name__, exc)
        bat.add("corrupting_%s_is_detected" % label, ok, detail)


def note(text):
    print("NOTE " + text)


def show_object(sp, code):
    """Print the decoded exhibited object back out."""
    print("# EXHIBITED OBJECT (decoded from the paper)")
    print("#   E = F_4^4, w^2+w+1 = 0, B = PG(F_2^4) has %d points, "
          "|PG(3,4)| = %d" % (len(BAER), len(POINTS)))
    for i in range(4):
        print("#   T row %d = %s" % (i + 1, [F4NAME[c] for c in sp.T[i]]))
    print("#   u  = %s   v  = %s" % ([F4NAME[c] for c in sp.u],
                                     [F4NAME[c] for c in sp.v]))
    print("#   u' = %s   v' = %s" % ([F4NAME[c] for c in sp.up],
                                     [F4NAME[c] for c in sp.vp]))
    for i in range(len(sp.lines)):
        pts = sorted(sp.lines[i])
        print("#   l_%-2d = %s" % (i + 1, " ".join(
            "".join(F4NAME[c] for c in p) for p in pts)))
    if code is not None and code.ok:
        enum = {}
        for w in set(code.words.values()):
            k = sum(1 for s in w if s != ZERO16)
            enum[k] = enum.get(k, 0) + 1
        print("# recomputed weight enumerator: " + " + ".join(
            "%dz^%d" % (enum[k], k) for k in sorted(enum)))


def main(argv):
    print("# verify.py -- A Length-14 Negative Answer to Alderson's "
          "q=4 Optimality Question")
    # The primary battery must never abort: a corrupt exhibited object has to
    # come out as reported FAILures, not as a traceback with no verdict line.
    bat = Battery()
    sp = code = None
    aborted = None
    try:
        _, sp, code = run_battery(PAPER_T, PAPER_U, PAPER_V, PAPER_UP,
                                  PAPER_VP, bat=bat)
    except Exception as exc:
        aborted = "%s: %s" % (type(exc).__name__, exc)
    bat.add("primary_battery_ran_to_completion_without_error", aborted is None,
            aborted or "no exception raised while checking the paper's object")
    if sp is not None:
        show_object(sp, code)
    checks_mutation_sensitivity(bat)
    checks_spread_positive_control(bat)
    print("# NOTE: every enumeration below is exhaustive and complete: all 85 "
          "points and all 357 lines of PG(3,4), all 256 codewords of C and of "
          "C^+, all 357 two-dimensional F_4-subspaces of E, and all 4^8 = "
          "65536 F_4-linear maps E -> F_16.  Nothing in the paper's claim was "
          "sampled, truncated or assumed.")
    for name, ok, detail in bat.rows:
        line = ("PASS " if ok else "FAIL ") + name
        if detail:
            line += " [" + detail + "]"
        print(line)
    # The limitations must be readable from the TRANSCRIPT, not only from this
    # file's docstring: a referee who sees the PASS lines and nothing else
    # would otherwise have no hint that anything above rests on an input.
    note("NOT RE-RUN: five things above are inputs taken from the literature "
         "and are recomputed by no check in this transcript. "
         "(1) Alderson's Problem 9.2 itself, together with his Theorem 6.6 "
         "and its parameters (112,2,104)_{16/4}: this program reads no "
         "external paper, so the statement being answered and the comparison "
         "length 112 are transcriptions, and a misreading of either would be "
         "caught by nothing above. "
         "(2) The definitions of 'extendable' and of 'additively maximal': "
         "the 7 properties counted in "
         "problem_9_2_hypotheses_all_verified_for_this_code are this "
         "program's rendering of those definitions, not a derivation of them. "
         "(3) The Remark's projective-equivalence claim. The 112 lines "
         "external to B counted above are external to THIS B; that Alderson's "
         "scattered F_2-linear set for t=2 is a Baer subgeometry projectively "
         "equivalent to it, and hence that these fourteen lines are fourteen "
         "of his own 112 external lines, is not checked here. "
         "(4) The spread classifications of van Dam and of Mellinger. The "
         "fourteen explicit lines are verified above to partition Pi together "
         "with B, which is a self-contained substitute, but no classification "
         "is confirmed and it is not tested that this spread lies in one of "
         "their three orbits. "
         "(5) Faithfulness, tested only in the paper's own sense that every "
         "coordinate projection is onto F_16; a stronger notion is not "
         "tested. "
         "Also NOT attempted: any lower bound, and any search for the global "
         "minimum length -- what is established above is 14 < 112, not that "
         "14 is optimal, and the nonsquare, nonprime clause of Problem 9.2 is "
         "untouched. Apart from the five items above, the only paper values "
         "this program reads are the three determinant tuples and the weight "
         "enumerator W_C, and those are quoted solely as comparison targets "
         "printed beside independently recomputed values.")
    nfail = len(bat.failures())
    ntot = len(bat.rows)
    if nfail == 0:
        print("VERDICT: ALL %d CHECKS PASS" % ntot)
        return 0
    print("VERDICT: %d OF %d CHECKS FAILED" % (nfail, ntot))
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))

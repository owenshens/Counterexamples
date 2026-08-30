#!/usr/bin/env python3
"""verify.py -- re-derives every computational claim of

    "Two Zonotopes with the Same Arithmetic Matroid and Different Graded Ehrhart Series"
    (paper.tex / paper.pdf in this folder)

from the objects PRINTED IN THAT PAPER and from nothing else.  Python 3.9+, standard library
only (fractions, itertools, math) -- no numpy, no sympy, no external data file, so a referee can
run it with nothing installed.  Every decision is made in exact integer or exact rational
arithmetic; there is no floating point anywhere and no randomness, hence no seed.

WHAT IT READS.  The two generator matrices, the two lattice-point sets, the eight graded
coefficient lists, the eight separating cubics with their claimed values, and the control values
quoted from the literature are all typed in below EXACTLY as the paper prints them, in the block
marked "THE PAPER'S CLAIMS".  Everything else the program derives.

WHAT IT DOES NOT COVER is printed by the program itself, as `NOT RE-RUN:` lines just before the
verdict, and repeated in REVIEW_NOTE.md under `## Scope`.
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd

# =====================================================================================
# THE PAPER'S CLAIMS -- typed in from paper.tex, not computed
# =====================================================================================
# Section 2, the witness.  Columns of the two 2x2 integer generator matrices A_1, A_2.
A1 = [(1, 0), (1, 5)]
A2 = [(1, 0), (2, 5)]

# Section 2, the shared arithmetic matroid: mu on the four subsets of {1,2}, identity labelling.
CLAIM_MU = {(): 1, (0,): 1, (1,): 1, (0, 1): 5}
CLAIM_TUTTE_X1 = {0: 4, 1: 5, 2: 8, 3: 13}          # T_M(x,1) = x^2 + 4, evaluated at x = 0,1,2,3

# Section 2, the two lattice-point sets at m = 1, as printed.
CLAIM_S1 = [(0, 0), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 5)]
CLAIM_S2 = [(0, 0), (1, 0), (1, 1), (1, 2), (2, 3), (2, 4), (2, 5), (3, 5)]

# Section 2, Table 1: the graded Ehrhart coefficient lists, ascending from q^0.
CLAIM_SERIES = {
    ('Z1', 1): [1, 2, 2, 1, 1, 1],
    ('Z1', 2): [1, 2, 3, 4, 4, 3, 3, 2, 1, 1, 1],
    ('Z1', 3): [1, 2, 3, 4, 5, 6, 6, 5, 5, 4, 3, 3, 2, 1, 1, 1],
    ('Z1', 4): [1, 2, 3, 4, 5, 6, 7, 8, 8, 7, 7, 6, 5, 5, 4, 3, 3, 2, 1, 1, 1],
    ('Z2', 1): [1, 2, 3, 2],
    ('Z2', 2): [1, 2, 3, 4, 5, 6, 4],
    ('Z2', 3): [1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 1],
    ('Z2', 4): [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 9, 2],
}
CLAIM_COUNTS = {1: 8, 2: 25, 3: 52, 4: 89}          # |mZ cap Z^2| = 5m^2 + 2m + 1

# Section 3, Table 2: the eight separating cubics for S_2 and their values at the omitted point.
# Each entry is (index omitted, printed expression, callable, printed value at the omitted point).
CLAIM_CUBICS = [
    (0, '(x-1)(x-2)(x-3)',      lambda x, y: (x - 1) * (x - 2) * (x - 3),            -6),
    (1, '(y-2x)(x-2)(y-2x+1)',  lambda x, y: (y - 2 * x) * (x - 2) * (y - 2 * x + 1), -2),
    (2, 'y(x-2)(2y-3x-1)',      lambda x, y: y * (x - 2) * (2 * y - 3 * x - 1),        2),
    (3, 'y(x-2)(y-2x+1)',       lambda x, y: y * (x - 2) * (y - 2 * x + 1),           -2),
    (4, '(y-5)(x-1)(y-2x)',     lambda x, y: (y - 5) * (x - 1) * (y - 2 * x),          2),
    (5, '(y-5)(x-1)(3x-2y)',    lambda x, y: (y - 5) * (x - 1) * (3 * x - 2 * y),      2),
    (6, '(y-2x)(x-1)(y-2x+1)',  lambda x, y: (y - 2 * x) * (x - 1) * (y - 2 * x + 1),  2),
    (7, 'x(x-1)(x-2)',          lambda x, y: x * (x - 1) * (x - 2),                    6),
]

# Section 5, the controls quoted from the literature: (name, generators, m, coefficient list).
CLAIM_CONTROLS = [
    ('CP-hexagon-Paper.tex-311',  [(1, 0), (0, 1), (1, 1)], 1, [1, 2, 3, 1]),
    ('CP-Example-7.2-diamond',    [(1, -1), (1, 1)],        1, [1, 2, 2]),
    ('RR24-anc-line-157',         [(1, 0), (1, 3)],         1, [1, 2, 2, 1]),
    ('RR24-anc-line-146',         [(1, 0), (0, 3)],         1, [1, 2, 2, 2, 1]),
]
# Section 5, the forced-positive control at fixed cardinality (point sets, not zonotopes).
CLAIM_FORCED = [(((0, 0), (1, 0), (2, 0)), [1, 1, 1]), (((0, 0), (1, 0), (0, 1)), [1, 2])]
# Section 5, the anti-collapse control: these residues return the OTHER member's series.
CLAIM_ANTICOLLAPSE = [([(1, 0), (3, 5)], 'Z2'), ([(1, 0), (4, 5)], 'Z1')]

# Section 4, the family and the two caveats.
CLAIM_FAMILY_ODD_K = list(range(5, 16, 2))
CLAIM_FAMILY_ALL_K = list(range(5, 21))
CLAIM_NO_WITNESS_K = 6
CLAIM_BAD_CHOICE = (30, 5)                          # k = 30, a = 5: column (5,30) is imprimitive
CLAIM_MINIMALITY_DET = 4                            # no witness in this cell with |det| <= 4
CLAIM_BOX = 4                                       # normal-form validation box: entries in [-4,4]

# =====================================================================================
# 0. THE CHECK HARNESS
# =====================================================================================
_N = [0]
_BAD = [0]


def ok(name, detail=''):
    _N[0] += 1
    print('PASS %s%s' % (name, ('  ' + detail) if detail else ''))


def bad(name, detail=''):
    _BAD[0] += 1
    print('CHECK-FAILED %s%s' % (name, ('  ' + detail) if detail else ''))


def want(cond, name, detail=''):
    if cond:
        ok(name, detail)
    else:
        bad(name, detail)


# =====================================================================================
# 1. EXACT LINEAR ALGEBRA OVER Q
# =====================================================================================
_MONS = {}


def monomials_upto(j):
    """All monomials x^a y^b with a + b <= j, ordered by total degree then by a."""
    if j not in _MONS:
        _MONS[j] = [(a, t - a) for t in range(j + 1) for a in range(t + 1)]
    return _MONS[j]


def monomials_deg(j):
    return [(a, j - a) for a in range(j + 1)]


def _evalmat(S, mons):
    return [[Fraction(p[0] ** m[0] * p[1] ** m[1]) for m in mons] for p in S]


def rank_profile(S, jmax):
    """[r_0, ..., r_jmax] with r_j = rank of the evaluation matrix of monomials of degree <= j.

    Columns are eliminated in monomial order, so the pivot count after the last column of degree
    j IS r_j: one elimination yields the whole profile.  Exact rational arithmetic.
    """
    mons = monomials_upto(jmax)
    M = _evalmat(S, mons)
    nr, nc = len(M), len(mons)
    ends = {len(monomials_upto(j)) - 1: j for j in range(jmax + 1)}
    prof, r = [], 0
    for c in range(nc):
        if r < nr:
            piv = None
            for i in range(r, nr):
                if M[i][c] != 0:
                    piv = i
                    break
            if piv is not None:
                M[r], M[piv] = M[piv], M[r]
                pv, Mr = M[r][c], M[r]
                for i in range(r + 1, nr):
                    if M[i][c] != 0:
                        f, Mi = M[i][c] / pv, M[i]
                        for k in range(c, nc):
                            if Mr[k]:
                                Mi[k] -= f * Mr[k]
                r += 1
        if c in ends:
            prof.append(r)
    return prof


def rank_rows(rows, nc):
    if not rows:
        return 0
    M = [list(map(Fraction, row)) for row in rows]
    nr, r = len(M), 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if M[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv, Mr = M[r][c], M[r]
        for i in range(r + 1, nr):
            if M[i][c] != 0:
                f, Mi = M[i][c] / pv, M[i]
                for k in range(c, nc):
                    if Mr[k]:
                        Mi[k] -= f * Mr[k]
        r += 1
        if r == nr:
            break
    return r


def kernel(rows, nc):
    """A basis of the null space of the given matrix, exact over Q."""
    M = [list(map(Fraction, row)) for row in rows]
    nr, pivots, r = len(M), [], 0
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
        Mr = M[r]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f, Mi = M[i][c], M[i]
                for k in range(c, nc):
                    if Mr[k]:
                        Mi[k] -= f * Mr[k]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    basis = []
    for fc in [c for c in range(nc) if c not in pivots]:
        v = [Fraction(0)] * nc
        v[fc] = Fraction(1)
        for i, pc in enumerate(pivots):
            v[pc] = -M[i][fc]
        basis.append(v)
    return basis


# =====================================================================================
# 2. THE GRADED EHRHART SERIES, TWO INDEPENDENT ROUTES
# =====================================================================================
def hilbert_at_degree(S, deg):
    """The coefficient list of Hilb(Orb(S); q) ASSUMING the paper's claimed top degree `deg`,
    together with the two facts that pin that degree: r_deg = |S| and r_{deg-1} < |S|."""
    prof = rank_profile(S, deg)
    h, prev = [], 0
    for rj in prof:
        h.append(rj - prev)
        prev = rj
    return h, prof


def hilbert(S):
    """The coefficient list with no claimed degree: the profile is extended until it saturates."""
    n, j = len(S), 1
    while True:
        prof = rank_profile(S, j)
        if prof[-1] == n:
            break
        j += 2
        if j > 90:
            raise AssertionError('runaway degree')
    h, prev = [], 0
    for rj in prof:
        h.append(rj - prev)
        prev = rj
        if prev == n:
            break
    return h


def top_form_spaces(S, jmax):
    """[W_0, ..., W_jmax], each a list of coefficient vectors on monomials_deg(j).

    W_j = { top-degree-j homogeneous component of f : f in I(S), deg f <= j }.  This is the
    SECOND route to the graded Ehrhart series and it does not use the rank-differencing identity:
    dim Orb(S)_j = dim C[x,y]_j - dim (gr I(S))_j = (j+1) - dim W_j, valid once (gr I(S))_j = W_j,
    which `closed_under_multiplication` below verifies rather than assumes.
    """
    Ws = []
    for j in range(jmax + 1):
        mons = monomials_upto(j)
        idx = {m: i for i, m in enumerate(mons)}
        K = kernel(_evalmat(S, mons), len(mons))
        top = monomials_deg(j)
        Ws.append([[v[idx[m]] for m in top] for v in K])
    return Ws


def closed_under_multiplication(Ws):
    """x * W_{j-1} + y * W_{j-1} contained in W_j for every j -- the hypothesis that makes the
    ideal generated by the top forms equal to the graded pieces W_j degree by degree."""
    for j in range(1, len(Ws)):
        base, d = Ws[j], rank_rows(Ws[j], j + 1)
        pos = {m: i for i, m in enumerate(monomials_deg(j))}
        extra = []
        for v in Ws[j - 1]:
            xv = [Fraction(0)] * (j + 1)
            yv = [Fraction(0)] * (j + 1)
            for i, (a, b) in enumerate(monomials_deg(j - 1)):
                xv[pos[(a + 1, b)]] = v[i]
                yv[pos[(a, b + 1)]] = v[i]
            extra += [xv, yv]
        if extra and rank_rows(base + extra, j + 1) != d:
            return False
    return True


def orb_from_top_forms(S, jmax):
    Ws = top_form_spaces(S, jmax)
    return [(j + 1) - rank_rows(Ws[j], j + 1) for j in range(jmax + 1)], Ws


# =====================================================================================
# 3. LATTICE POINTS OF A ZONOTOPE, TWO INDEPENDENT ROUTES
# =====================================================================================
def points_support(gens, m):
    """Pure-INTEGER route: p in m*Z(v_1..v_n) iff <u,p> <= m * sum_i max(0, <u,v_i>) for every
    edge normal u = +-(v_i[1], -v_i[0]).  No division anywhere."""
    ns = sorted({u for v in gens for u in ((v[1], -v[0]), (-v[1], v[0])) if u != (0, 0)})
    sup = {u: m * sum(max(0, u[0] * v[0] + u[1] * v[1]) for v in gens) for u in ns}
    lo_x = m * sum(min(0, v[0]) for v in gens)
    hi_x = m * sum(max(0, v[0]) for v in gens)
    lo_y = m * sum(min(0, v[1]) for v in gens)
    hi_y = m * sum(max(0, v[1]) for v in gens)
    out = []
    for x in range(lo_x, hi_x + 1):
        for y in range(lo_y, hi_y + 1):
            if all(u[0] * x + u[1] * y <= sup[u] for u in ns):
                out.append((x, y))
    return sorted(out)


def points_inverse(gens, m):
    """Exact-RATIONAL route, for n = 2 only: p in m*Z(A) iff A^{-1}p/m lies in [0,1]^2."""
    (a, c), (b, d) = gens
    det = a * d - b * c
    if det == 0:
        raise ValueError('singular')
    xs = [0, m * a, m * b, m * (a + b)]
    ys = [0, m * c, m * d, m * (c + d)]
    out = []
    for x in range(min(xs), max(xs) + 1):
        for y in range(min(ys), max(ys) + 1):
            t1 = Fraction(d * x - b * y, det * m)
            t2 = Fraction(a * y - c * x, det * m)
            if 0 <= t1 <= 1 and 0 <= t2 <= 1:
                out.append((x, y))
    return sorted(out)


# =====================================================================================
# 4. THE ARITHMETIC MATROID, BY SMITH NORMAL FORM
# =====================================================================================
def snf_divisors(mat):
    """The nonzero elementary divisors of an integer matrix, by integer row/column reduction."""
    M = [row[:] for row in mat]
    nr = len(M)
    nc = len(M[0]) if nr else 0
    divs, r, c = [], 0, 0
    while r < nr and c < nc:
        piv, best = None, None
        for i in range(r, nr):
            for j in range(c, nc):
                if M[i][j] != 0 and (best is None or abs(M[i][j]) < best):
                    best, piv = abs(M[i][j]), (i, j)
        if piv is None:
            break
        i0, j0 = piv
        M[r], M[i0] = M[i0], M[r]
        for row in M:
            row[c], row[j0] = row[j0], row[c]
        while True:
            changed = False
            for i in range(r + 1, nr):
                if M[i][c] != 0:
                    q = M[i][c] // M[r][c]
                    for k in range(c, nc):
                        M[i][k] -= q * M[r][k]
                    if M[i][c] != 0:
                        M[r], M[i] = M[i], M[r]
                        changed = True
            for j in range(c + 1, nc):
                if M[r][j] != 0:
                    q = M[r][j] // M[r][c]
                    for i in range(r, nr):
                        M[i][j] -= q * M[i][c]
                    if M[r][j] != 0:
                        for i in range(r, nr):
                            M[i][c], M[i][j] = M[i][j], M[i][c]
                        changed = True
            if not changed:
                break
        divs.append(abs(M[r][c]))
        r, c = r + 1, c + 1
    return [d for d in divs if d]


def _cols_matrix(cols):
    return [[cols[j][i] for j in range(len(cols))] for i in range(2)]


def mult(cols):
    """mu(S) = the index of the lattice spanned by S in (real span of S) cap Z^2 = the product of
    the elementary divisors of the matrix whose columns are S.  mu(empty) = 1."""
    if not cols:
        return 1
    p = 1
    for d in snf_divisors(_cols_matrix(cols)):
        p *= d
    return p


def arithmetic_matroid(gens):
    """((subset, rank, mu), ...) over every subset of the ground set, identity labelling."""
    out = []
    for k in range(len(gens) + 1):
        for S in combinations(range(len(gens)), k):
            cols = [gens[i] for i in S]
            out.append((S, len(snf_divisors(_cols_matrix(cols))) if cols else 0, mult(cols)))
    return tuple(out)


def arithmetic_tutte_at(gens, x, y):
    """T_M(x,y) = sum_S mu(S) (x-1)^{r(E)-r(S)} (y-1)^{|S|-r(S)}, exact rational."""
    am = arithmetic_matroid(gens)
    rE = max(r for _S, r, _m in am)
    tot = Fraction(0)
    for S, r, m in am:
        tot += m * Fraction(x - 1) ** (rE - r) * Fraction(y - 1) ** (len(S) - r)
    return tot


# =====================================================================================
# 5. HERMITE NORMAL FORM AND THE CENSUS
# =====================================================================================
def hnf(v1, v2):
    """(a, b, d) with left-GL_2(Z) normal form [[a,b],[0,d]], a > 0, d > 0, 0 <= b < d, for the
    2x2 matrix whose COLUMNS are v1, v2.  Unique, so it is a complete invariant of the pair up to
    the GL_2(Z) action."""
    a, c = v1
    b, d = v2
    while c != 0:
        q = a // c
        a, c = c, a - q * c
        b, d = d, b - q * d
    if a < 0:
        a, b = -a, -b
    if d < 0:
        d = -d
    return (a, b % d, d)


def hnf_gens(h):
    a, b, d = h
    return [(a, 0), (b, d)]


def series1(gens):
    return tuple(hilbert(points_support(gens, 1)))


def hnf_census(k):
    """Every left-GL_2(Z) class of ordered pairs in Z^2 with |det| = k, with its arithmetic
    matroid (up to relabelling the two generators) and its m = 1 graded series."""
    rows = []
    for a in range(1, k + 1):
        if k % a:
            continue
        d = k // a
        for b in range(d):
            g = hnf_gens((a, b, d))
            am = arithmetic_matroid(g)
            key = (tuple(sorted(m for S, _r, m in am if len(S) == 1)), k)
            rows.append(((a, b, d), key, series1(g)))
    return rows


def split_types(k):
    groups = {}
    for _h, key, s in hnf_census(k):
        groups.setdefault(key, set()).add(s)
    return groups, {kk: v for kk, v in groups.items() if len(v) > 1}


# =====================================================================================
# 6. THE CHECKS
# =====================================================================================
def main():
    print('verifying "Two Zonotopes with the Same Arithmetic Matroid and Different Graded '
          'Ehrhart Series"')
    print('objects read from the paper: A_1 columns %s, A_2 columns %s' % (A1, A2))
    print('')

    # --- 6.1 the arithmetic matroid ---------------------------------------------------
    am1, am2 = arithmetic_matroid(A1), arithmetic_matroid(A2)
    want(am1 == am2, 'arithmetic-matroids-equal',
         'identical on all 4 subsets under the identity labelling: %s'
         % [(S, r, m) for S, r, m in am1])
    want({S: m for S, _r, m in am1} == CLAIM_MU and {S: m for S, _r, m in am2} == CLAIM_MU,
         'multiplicity-matches-paper', 'mu = %s, by Smith normal form (not the gcd/det shortcut)'
         % {('{%s}' % ','.join(str(i + 1) for i in S)) or '{}': m for S, _r, m in am1})
    dets = [abs(A1[0][0] * A1[1][1] - A1[1][0] * A1[0][1]),
            abs(A2[0][0] * A2[1][1] - A2[1][0] * A2[0][1])]
    want(dets == [5, 5] and all(gcd(v[0], v[1]) == 1 for v in A1 + A2),
         'all-four-columns-primitive', '|det| = %s and every column gcd is 1' % dets)
    tut = {x: arithmetic_tutte_at(A1, x, 1) for x in CLAIM_TUTTE_X1}
    want(tut == {x: Fraction(v) for x, v in CLAIM_TUTTE_X1.items()}
         and {x: arithmetic_tutte_at(A2, x, 1) for x in CLAIM_TUTTE_X1} == tut,
         'shared-arithmetic-tutte', 'T_M(x,1) = x^2 + 4 for both, checked at x = 0,1,2,3: %s'
         % {x: int(v) for x, v in sorted(tut.items())})

    # --- 6.2 the lattice points ------------------------------------------------------
    pts = {}
    for nm, A in (('Z1', A1), ('Z2', A2)):
        for m in (1, 2, 3, 4):
            ps, pi = points_support(A, m), points_inverse(A, m)
            if ps != pi:
                bad('lattice-points-two-routes', '%s m=%d disagree' % (nm, m))
                return 1
            pts[(nm, m)] = ps
    want(True, 'lattice-points-two-routes',
         'the integer support-function test and the exact-Fraction inverse-matrix test return '
         'identical sets for both zonotopes at m = 1,2,3,4')
    want([list(p) for p in pts[('Z1', 1)]] == [list(p) for p in CLAIM_S1]
         and [list(p) for p in pts[('Z2', 1)]] == [list(p) for p in CLAIM_S2],
         'point-sets-match-paper', 'S_1 and S_2 are exactly the 8+8 points printed in Section 2')
    want(all(len(pts[(nm, m)]) == CLAIM_COUNTS[m] for nm in ('Z1', 'Z2') for m in (1, 2, 3, 4)),
         'ungraded-counts-agree', '|mZ_1| = |mZ_2| = %s for m = 1,2,3,4'
         % [CLAIM_COUNTS[m] for m in (1, 2, 3, 4)])
    want(all(CLAIM_COUNTS[m] == 5 * m * m + 2 * m + 1 for m in (1, 2, 3, 4))
         and all(m * m * arithmetic_tutte_at(A, Fraction(m + 1, m), 1) == CLAIM_COUNTS[m]
                 for A in (A1, A2) for m in (1, 2, 3, 4)),
         'dadderio-moci-identity',
         'm^2 T_M((m+1)/m, 1) = 5m^2 + 2m + 1 = the count, both zonotopes, m = 1,2,3,4')
    pick = []
    for nm, A in (('Z1', A1), ('Z2', A2)):
        for m in (1, 2, 3, 4):
            area, B = 5 * m * m, 4 * m         # primitive edges, so each of the 4 edges has m
            I = area - B // 2 + 1              # Pick: area = I + B/2 - 1
            pick.append(I + B == len(pts[(nm, m)]))
    want(all(pick), 'picks-theorem',
         'area 5m^2, B = 4m (all four edges primitive), I = 5m^2 - 2m + 1, I + B = the count')

    # --- 6.3 the graded series -------------------------------------------------------
    got = {}
    for nm in ('Z1', 'Z2'):
        for m in (1, 2, 3, 4):
            claim = CLAIM_SERIES[(nm, m)]
            S = pts[(nm, m)]
            h, prof = hilbert_at_degree(S, len(claim) - 1)
            got[(nm, m)] = h
            want(h == claim and prof[-1] == len(S) and prof[-2] < len(S),
                 'graded-series-%s-m%d' % (nm, m),
                 '%s; degree is exactly %d (r_%d = %d = |S|, r_%d = %d)'
                 % (h, len(claim) - 1, len(claim) - 1, prof[-1], len(claim) - 2, prof[-2]))
    want(all(sum(got[(nm, m)]) == CLAIM_COUNTS[m] for nm in ('Z1', 'Z2') for m in (1, 2, 3, 4)),
         'no-series-truncated',
         'every coefficient list sums to |mZ cap Z^2|, which is the total-dimension property of '
         'Orb; a truncated series would fail here')
    want(all(got[('Z1', m)] != got[('Z2', m)] for m in (1, 2, 3, 4)),
         'series-differ-at-every-m', 'm = 1,2,3,4 all separate the pair')
    want(got[('Z1', 1)][2] == 2 and got[('Z2', 1)][2] == 3
         and len(got[('Z1', 1)]) - 1 == 5 and len(got[('Z2', 1)]) - 1 == 3,
         'the-decisive-difference',
         'q^2 coefficient 2 vs 3, degree 5 vs 3, at m = 1 -- one m suffices, so the first clause '
         'of Question 7.3 is answered NO')
    want(all(sum(got[('Z1', m)]) == sum(got[('Z2', m)]) for m in (1, 2, 3, 4)),
         'equal-at-q-equals-1',
         'the decider returns EQUAL at q = 1 on the same pair (%s), so it is not a negative-only '
         'test' % [sum(got[('Z1', m)]) for m in (1, 2, 3, 4)])

    # --- 6.4 the graded series a second way, and the gr-ideal hypothesis -------------
    second = []
    for nm in ('Z1', 'Z2'):
        for m in (1, 2):
            S = pts[(nm, m)]
            jm = len(CLAIM_SERIES[(nm, m)]) - 1
            h2, Ws = orb_from_top_forms(S, jm)
            second.append((h2 == got[(nm, m)], closed_under_multiplication(Ws)))
    want(all(a for a, _b in second), 'series-second-route',
         'dim Orb_j = (j+1) - dim W_j, computed from an explicit basis of I(S)_{<=j} and its '
         'degree-j top forms, agrees with rank differencing for both zonotopes at m = 1,2')
    want(all(b for _a, b in second), 'gr-ideal-equals-top-forms',
         'x*W_{j-1} + y*W_{j-1} is contained in W_j at every j, so the ideal generated by the top '
         'forms has degree-j part exactly W_j -- the identity dim Orb_j = r_j - r_{j-1} is '
         'verified, not assumed')

    # --- 6.5 the hand certificate for S_1 -------------------------------------------
    S1 = [tuple(p) for p in CLAIM_S1]
    prof1 = rank_profile(S1, 5)
    Nj = lambda j: (j + 1) * (j + 2) // 2
    want(prof1 == [1, 3, 5, 6, 7, 8] and all(prof1[j] == j + 3 for j in range(2, 6)),
         'certificate-S1-rank-profile', 'r_j = j + 3 for 2 <= j <= 5, r_0 = 1, r_1 = 3: %s'
         % prof1)
    forced = []
    for j in range(2, 6):
        # dim I(S_1)_{<=j} from the rank, against the certificate's (x-1)*{g(0,0)=g(2,5)=0}
        dim_I = Nj(j) - prof1[j]
        mons = monomials_upto(j - 1)
        cond = [[Fraction(p[0] ** m[0] * p[1] ** m[1]) for m in mons] for p in [(0, 0), (2, 5)]]
        dim_cof = len(mons) - rank_rows(cond, len(mons))
        forced.append((dim_I, dim_cof, Nj(j - 1) - 2))
    want(all(a == b == c for a, b, c in forced), 'certificate-S1-forced-line-factor',
         'dim I(S_1)_{<=j} = dim (x-1)*{deg <= j-1, vanishing at (0,0) and (2,5)} = N_{j-1} - 2 '
         'for j = 2..5, so the inclusion of the certificate is an equality: %s' % forced)
    want(sum(1 for p in S1 if p[0] == 1) == 6,
         'certificate-S1-six-collinear', 'six of the eight points of S_1 lie on x = 1, so a form '
         'of degree j <= 5 vanishing on S_1 has 6 >= j+1 roots there and is divisible by (x-1)')

    # --- 6.6 the hand certificate for S_2 -------------------------------------------
    S2 = [tuple(p) for p in CLAIM_S2]
    prof2 = rank_profile(S2, 3)
    want(prof2 == [1, 3, 6, 8], 'certificate-S2-rank-profile',
         'r_0..r_3 = %s; r_2 = 6 says no conic passes through S_2 and r_3 = 8 says the eight '
         'points impose independent conditions on cubics' % prof2)
    col1 = [p for p in S2 if p[0] == 1]
    col2 = [p for p in S2 if p[0] == 2]
    want(len(col1) == 3 and len(col2) == 3
         and (0 - 1) * (0 - 2) == 2,
         'certificate-S2-no-conic-argument',
         'x = 1 carries %s and x = 2 carries %s, so a conic through S_2 is a multiple of '
         '(x-1)(x-2), which takes the value 2 at (0,0)' % (col1, col2))
    cub = []
    for i, expr, f, val in CLAIM_CUBICS:
        vals = [f(*p) for p in S2]
        cub.append((all(v == 0 for j, v in enumerate(vals) if j != i), vals[i] == val,
                    expr, vals[i]))
    want(all(a and b for a, b, _e, _v in cub), 'certificate-S2-eight-separating-cubics',
         'each of the eight printed cubics vanishes at the other seven points and takes its '
         'printed value at the omitted one: %s'
         % [(e, v) for _a, _b, e, v in cub])
    coll = max(sum(1 for p in S2 if (b[0] - a[0]) * (p[1] - a[1]) == (b[1] - a[1]) * (p[0] - a[0]))
               for a, b in combinations(S2, 2))
    want(len(S2) == 8 == 2 * 3 + 2 and coll == 3 < 3 + 2 and prof2[2] == 6,
         'certificate-S2-cayley-bacharach',
         'the Ellia-Peskine / Eisenbud-Green-Harris route with form degree d = 3: |S_2| = 8 = '
         '2d+2, the largest collinear subset has %d < d+2 = 5 points, and S_2 lies on no conic, '
         'so the conditions are independent and r_3 = 8' % coll)
    sig = sorted((3 - x, 5 - y) for x, y in S2)
    want(sig == sorted(S2), 'certificate-S2-central-symmetry',
         'sigma(x,y) = (3-x, 5-y) maps S_2 onto itself, pairing P1<->P8, P2<->P7, P3<->P6, '
         'P4<->P5; each cubic is nevertheless printed explicitly')

    # --- 6.7 lattice inequivalence --------------------------------------------------
    equiv = []
    for order in ((0, 1), (1, 0)):
        for s0 in (1, -1):
            for s1 in (1, -1):
                B = [tuple(s0 * t for t in A2[order[0]]), tuple(s1 * t for t in A2[order[1]])]
                (a, c), (b, d) = A1
                det = a * d - b * c
                U = [[Fraction(B[0][i] * d - B[1][i] * c, det),
                      Fraction(-B[0][i] * b + B[1][i] * a, det)] for i in range(2)]
                integral = all(x.denominator == 1 for row in U for x in row)
                unimod = abs(U[0][0] * U[1][1] - U[0][1] * U[1][0]) == 1
                equiv.append(integral and unimod)
    want(not any(equiv), 'not-lattice-equivalent',
         'none of the 8 signed orderings of A_2 equals U A_1 for an integer U with |det U| = 1, '
         'so Z_1 and Z_2 are genuinely different lattice polygons (the residue a = 1 vs a = 2 '
         'mod 5); the direct argument U(1,0) = (1,0) forces b = 1/5')

    # --- 6.8 the controls -----------------------------------------------------------
    for nm, gens, m, claim in CLAIM_CONTROLS:
        S = points_support(gens, m)
        h = hilbert(S)
        want(h == claim and sum(h) == len(S), 'control-%s' % nm,
             'generators %s at m=%d: %d points, series %s -- matches the published value'
             % (gens, m, len(S), h))
    fp = [(hilbert(list(S)) == claim) for S, claim in CLAIM_FORCED]
    want(all(fp) and len(CLAIM_FORCED[0][0]) == len(CLAIM_FORCED[1][0]) == 3,
         'control-forced-positive-at-fixed-cardinality',
         'three collinear points give 1 + q + q^2 and three non-collinear give 1 + 2q; the '
         'decider separates two sets of the SAME cardinality')
    silent = []
    U8 = [((1, 0), (0, 1)), ((1, 0), (0, -1)), ((-1, 0), (0, 1)), ((-1, 0), (0, -1)),
          ((0, 1), (1, 0)), ((0, 1), (-1, 0)), ((0, -1), (1, 0)), ((0, -1), (-1, 0))]
    for U in U8:
        for s0 in (1, -1):
            for s1 in (1, -1):
                g = []
                for s, v in ((s0, A1[0]), (s1, A1[1])):
                    g.append((s * (U[0][0] * v[0] + U[0][1] * v[1]),
                              s * (U[1][0] * v[0] + U[1][1] * v[1])))
                for m in (1, 2):
                    silent.append(hilbert(points_support(g, m)) == CLAIM_SERIES[('Z1', m)])
    want(len(silent) == 64 and all(silent), 'control-must-stay-silent',
         '8 signed permutation matrices x 4 column sign patterns = 32 lattice equivalences of '
         'Z_1, all 32 reproducing Z_1\'s own series at m = 1 and m = 2 (64 comparisons)')
    anti = [(series1(g) == tuple(CLAIM_SERIES[(who, 1)])) for g, who in CLAIM_ANTICOLLAPSE]
    want(all(anti), 'control-anti-collapse',
         'the residue theory predicts (1,0),(3,5) -> Z_2\'s series (3 = 2^{-1} mod 5) and '
         '(1,0),(4,5) -> Z_1\'s series (4 = -1 mod 5); both hold, so the separation tracks the '
         'residue class and not the numerals')

    # --- 6.9 minimality, the k = 6 gap, and the normal form -------------------------
    nosplit = []
    for k in range(1, CLAIM_MINIMALITY_DET + 1):
        _g, sp = split_types(k)
        nosplit.append((k, len(_g), sp == {}))
    want(all(f for _k, _n, f in nosplit), 'minimality-no-witness-below-five',
         'exhaustive over ALL left-GL_2(Z) classes with |det| = k for k = 1,2,3,4: every '
         'arithmetic-matroid type carries exactly one m = 1 series. %s'
         % [(k, n) for k, n, _f in nosplit])
    groups5, split5 = split_types(5)
    got5 = sorted(list(split5.values())[0]) if len(split5) == 1 else None
    want(len(split5) == 1 and got5 == sorted([tuple(CLAIM_SERIES[('Z1', 1)]),
                                              tuple(CLAIM_SERIES[('Z2', 1)])]),
         'first-split-is-at-determinant-five',
         'at |det| = 5 exactly one of the %d arithmetic-matroid types splits, namely mu = (1,1), '
         'and it splits into exactly the two series of the paper: %s'
         % (len(groups5), got5))
    g6, sp6 = split_types(CLAIM_NO_WITNESS_K)
    want(sp6 == {}, 'no-witness-at-determinant-six',
         'the family misses k = 6 for a real reason: at |det| = 6 no arithmetic-matroid type '
         'splits, because (Z/6)^* = {1,5} collapses to one class under a -> -a')
    mism, seen, cache = 0, 0, {}
    B = CLAIM_BOX
    for a in range(-B, B + 1):
        for b in range(-B, B + 1):
            for c in range(-B, B + 1):
                for d in range(-B, B + 1):
                    det = a * d - b * c
                    if not (1 <= abs(det) <= 5):
                        continue
                    seen += 1
                    h = hnf((a, c), (b, d))
                    if h not in cache:
                        cache[h] = series1(hnf_gens(h))
                    if series1([(a, c), (b, d)]) != cache[h]:
                        mism += 1
    want(mism == 0, 'normal-form-validated-by-brute-force',
         'the Hermite normal form is a HYPOTHESIS about which pairs share a series, so it was '
         'checked: all %d integer 2x2 matrices with entries in [-%d,%d] and 1 <= |det| <= 5 '
         '(falling into %d normal forms) have the m = 1 series of their own normal form, %d '
         'mismatches' % (seen, B, B, len(cache), mism))

    # --- 6.10 the infinite family and its two caveats -------------------------------
    fam = []
    for k in CLAIM_FAMILY_ODD_K:
        g1, g2 = [(1, 0), (1, k)], [(1, 0), (2, k)]
        s1, s2 = series1(g1), series1(g2)
        fam.append((k, arithmetic_matroid(g1) == arithmetic_matroid(g2), s1 != s2,
                    len(s1) - 1 == k, len(s2) - 1 <= (k + 5) // 2))
    want(all(a and b and c and d for _k, a, b, c, d in fam), 'family-odd-k',
         'for every odd k in %s the pair (1,0),(1,k) and (1,0),(2,k) has equal arithmetic '
         'matroid and different series, with deg i_{Z_1(k)} = k exactly and deg i_{Z_2(k)} <= '
         '(k+5)/2' % CLAIM_FAMILY_ODD_K)
    want(all(len(series1([(1, 0), (2, k)])) - 1 < k for k in CLAIM_FAMILY_ODD_K if k >= 7),
         'family-degree-gap-is-strict',
         'for odd k >= 7 the bound (k+5)/2 is strictly below k, so the degrees alone separate '
         'the pair without computing a single coefficient; k = 5 needs the computation')
    mech = []
    for k in CLAIM_FAMILY_ODD_K:
        S = points_support([(1, 0), (2, k)], 1)
        cols = {}
        for x, y in S:
            cols[x] = cols.get(x, 0) + 1
        mech.append((k, sorted(cols), sorted(cols.values()),
                     sorted(cols) == [0, 1, 2, 3]
                     and max(cols.values()) == (k + 1) // 2
                     and len(series1([(1, 0), (2, k)])) - 1 <= 4 + (k + 1) // 2 - 2))
    want(all(f for _k, _c, _v, f in mech), 'family-degree-bound-mechanism',
         'the reason for the bound, checked: Z((1,0),(2,k)) cap Z^2 meets exactly the four '
         'vertical lines x = 0,1,2,3, with at most (k+1)/2 points on each, so line-by-line '
         'interpolation is surjective by degree 4 + (k+1)/2 - 2 = (k+5)/2. Line multiplicities: '
         '%s' % [(k, v) for k, _c, v, _f in mech])
    quant = []
    for k in CLAIM_FAMILY_ALL_K:
        a = next(x for x in range(2, k + 2) if gcd(x, k) == 1)
        g1, g2 = [(1, 0), (1, k)], [(1, 0), (a, k)]
        works = (arithmetic_matroid(g1) == arithmetic_matroid(g2)) and series1(g1) != series1(g2)
        quant.append((k, a, works))
    want(all(w for k, _a, w in quant if k != CLAIM_NO_WITNESS_K)
         and not [w for k, _a, w in quant if k == CLAIM_NO_WITNESS_K][0],
         'family-a-must-be-quantified-not-fixed',
         'taking a = the least integer >= 2 coprime to k gives a witness for every k in %s '
         'EXCEPT k = 6, which has none: %s'
         % (CLAIM_FAMILY_ALL_K, [(k, a, w) for k, a, w in quant]))
    kk, aa = CLAIM_BAD_CHOICE
    gbad = [(1, 0), (aa, kk)]
    mus = {S: m for S, _r, m in arithmetic_matroid(gbad)}
    want(mus[(1,)] == 5 != 1, 'family-fixed-a-would-be-wrong',
         'the caveat is real: at k = %d the choice a = %d makes the column (%d,%d) imprimitive, '
         'mu({2}) = %d, so the arithmetic matroids DIFFER and the pair is not a counterexample '
         'at all' % (kk, aa, aa, kk, mus[(1,)]))

    # --- the scope statement, in the program's own output ---------------------------
    print('')
    print('NOT RE-RUN: the SECOND clause of Question 7.3. Nothing above decides whether some '
          'invariant strictly between the arithmetic matroid and the full lattice class supports '
          'a generalisation of the paper\'s Proposition on Ehrhart polynomials; that stays open.')
    print('NOT RE-RUN: minimality outside the cell d = 2, n = 2. The census is exhaustive over '
          'GL_2(Z) classes of PAIRS only; no search over three or more generators is performed '
          'here, so nothing above excludes a smaller or differently shaped counterexample at '
          'n >= 3.')
    print('NOT RE-RUN: the independent top-form route to dim Orb_j (checks '
          '`series-second-route` and `gr-ideal-equals-top-forms`) covers m = 1 and m = 2 only; '
          'the m = 3 and m = 4 series rest on rank differencing alone.')
    print('NOT RE-RUN: bibliographic facts. The line numbers quoted for arXiv:2603.07873, the '
          'numbering of Question 7.3, the ancillary-file line numbers of arXiv:2407.06511 and '
          'every claim about the literature are read from sources, not computed here; only the '
          'NUMERICAL values taken from those sources are recomputed, as the `control-*` checks.')
    print('NOT RE-RUN: the infinite family is checked for the listed finite ranges of k, not '
          'proved for all k; the proof for all k is the argument in Section 4 of the paper.')
    print('')
    if _BAD[0]:
        print('VERDICT: %d CHECK(S) DID NOT PASS' % _BAD[0])
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _N[0])
    return 0


if __name__ == '__main__':
    sys.exit(main())

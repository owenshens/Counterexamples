#!/usr/bin/env python3
"""verify.py -- re-derives every computational claim of

    "Four distinct triangles force at most seven points, and the regular heptagon"

from the objects printed in paper.tex.  Python 3.9+, STANDARD LIBRARY ONLY (itertools,
fractions, math, sys): no third-party package, no external data file.

Everything decided here is decided in exact integer or Fraction arithmetic.  No floating
point number is computed anywhere in this file, and no comparison is made with a
tolerance.  The two algebraic number fields that occur are handled symbolically:

  * Q[c]/(c^3 + c^2 - 2c - 1),  c = 2cos(2*pi/7)   -- the regular heptagon,
  * Q                                              -- the hexagon-plus-centre, whose
    squared distances printed in the paper are the integers 1, 3, 4.

Layout, with the paper's numbering:

  A. Section 2   the regular heptagon: its squared chords, its four congruence classes.
  B. Section 2   the hexagon-plus-centre: 5 classes counted the paper's way, 4 counted
                 the rival way.  This is the object the Scope paragraph turns on.
  C. Lemma 4.1   the census of edge-colourings of K_7 in the cell.
  D. Lemma 4.2   the Cayley-Menger filter over the 161351 census leaves, and the S_7
                 orbit count of its survivors.
  E. Lemma 4.3   the exact elimination of 29 of the 30 orbits.
  F. Lemma 4.4   the surviving orbit is the gap colouring of Z_7, and Schoenberg's
                 criterion pins its squared distances to (s_1 : s_2 : s_3).
  G. Lemma 4.5   no 7-point planar set has at most three classes; the two-distance lane.

Closing `NOT RE-RUN:` lines say what this program does not cover.
"""
import itertools
import math
import sys
from fractions import Fraction

# ---------------------------------------------------------------------------
# check bookkeeping
# ---------------------------------------------------------------------------
_N_PASS = 0
_N_FAIL = 0


def check(name, ok, detail=''):
    global _N_PASS, _N_FAIL
    if ok:
        _N_PASS += 1
        print('PASS %-34s %s' % (name, detail))
    else:
        _N_FAIL += 1
        print('FAIL %-34s %s' % (name, detail))
    sys.stdout.flush()


# ---------------------------------------------------------------------------
# 0.  multivariate polynomials over Z, in y_2, y_3, y_4  (y_1 = 1 by scale)
#     exponent tuples of length 3; the empty dict is the zero polynomial
# ---------------------------------------------------------------------------
ONE = {(0, 0, 0): 1}


def pmul(a, b):
    r = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = (ea[0] + eb[0], ea[1] + eb[1], ea[2] + eb[2])
            r[e] = r.get(e, 0) + ca * cb
    return {e: c for e, c in r.items() if c}


def padd(*ps):
    r = {}
    for p in ps:
        for e, c in p.items():
            r[e] = r.get(e, 0) + c
    return {e: c for e, c in r.items() if c}


def pneg(a):
    return {e: -c for e, c in a.items()}


def pvar(colour):
    """the squared distance carried by colour `colour`; colour 1 is normalised to 1."""
    if colour == 1:
        return dict(ONE)
    e = [0, 0, 0]
    e[colour - 2] = 1
    return {tuple(e): 1}


def one_signed(p):
    """True iff p is a nonzero polynomial all of whose coefficients have one sign.
    Such a p is >0 at every point of the positive orthant, or <0 at every point, so it
    cannot vanish there: the colouring is not realizable.  (All exponents are >= 0.)"""
    if not p:
        return False
    return len(set(1 if c > 0 else -1 for c in p.values())) == 1


def det3(m):
    return padd(pmul(m[0][0], padd(pmul(m[1][1], m[2][2]), pneg(pmul(m[1][2], m[2][1])))),
                pneg(pmul(m[0][1], padd(pmul(m[1][0], m[2][2]), pneg(pmul(m[1][2], m[2][0]))))),
                pmul(m[0][2], padd(pmul(m[1][0], m[2][1]), pneg(pmul(m[1][1], m[2][0])))))


# ---------------------------------------------------------------------------
# 0b.  univariate rational polynomials, as coefficient lists low -> high
# ---------------------------------------------------------------------------
def _trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def ugcd(a, b):
    a, b = _trim([Fraction(x) for x in a]), _trim([Fraction(x) for x in b])
    while b:
        while a and len(a) >= len(b):
            f = a[-1] / b[-1]
            sh = len(a) - len(b)
            for i, c in enumerate(b):
                a[sh + i] -= f * c
            a = _trim(a)
        a, b = b, a
    if a:
        lc = a[-1]
        a = [x / lc for x in a]
    return a


def positive_rational_roots(coeffs):
    """the positive rational roots of a nonzero rational polynomial, by the rational root
    theorem after clearing denominators.  Used only where the caller also checks the
    degree, so no root is ever silently dropped."""
    h = _trim(coeffs)
    den = 1
    for c in h:
        den = den * c.denominator // math.gcd(den, c.denominator)
    H = [int(c * den) for c in h]
    a0, an = H[0], H[-1]
    cand = set()
    for p in range(1, abs(a0) + 1):
        if a0 % p:
            continue
        for q in range(1, abs(an) + 1):
            if an % q:
                continue
            cand.add(Fraction(p, q))
    return sorted(r for r in cand if sum(c * r ** i for i, c in enumerate(H)) == 0)


def poly_str(coeffs, name):
    """a monic rational polynomial as text, low -> high, for the transcript."""
    out = []
    for k, c in reversed(list(enumerate(_trim(coeffs)))):
        if not c:
            continue
        t = ('%s' % c) if (k == 0 or abs(c) != 1) else ('-' if c < 0 else '')
        if k == 1:
            t += name
        elif k > 1:
            t += '%s^%d' % (name, k)
        out.append(t)
    return ' + '.join(out).replace('+ -', '- ') or '0'


# ---------------------------------------------------------------------------
# A.  THE REGULAR HEPTAGON, in Q[c]/(c^3 + c^2 - 2c - 1),  c = 2cos(2*pi/7)
# ---------------------------------------------------------------------------
MINPOLY = (-1, -2, 1, 1)            # -1 - 2c + c^2 + c^3, low -> high
Z3 = (Fraction(0),) * 3


def rred(p):
    """reduce a coefficient list modulo c^3 = -c^2 + 2c + 1."""
    p = [Fraction(x) for x in p] + [Fraction(0)] * (6 - len(p))
    for k in range(len(p) - 1, 2, -1):
        a = p[k]
        if a:
            p[k] = Fraction(0)
            p[k - 1] -= a
            p[k - 2] += 2 * a
            p[k - 3] += a
    return tuple(p[:3])


def radd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def rsub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def rmul(a, b):
    r = [Fraction(0)] * 5
    for i in range(3):
        for j in range(3):
            r[i + j] += a[i] * b[j]
    return rred(r)


def rsc(a, k):
    return tuple(x * Fraction(k) for x in a)


RONE = (Fraction(1), Fraction(0), Fraction(0))
RC = (Fraction(0), Fraction(1), Fraction(0))


def section_A():
    print('--- A.  the regular heptagon (paper, Section 2) ---')
    # A1  the minimal polynomial of 2cos(2pi/7) is irreducible over Q, so the ring is a field
    bad = [r for r in (1, -1) if sum(c * r ** i for i, c in enumerate(MINPOLY)) == 0]
    check('heptagon-minpoly-irreducible', not bad,
          'c^3+c^2-2c-1 has no rational root (candidates +-1), hence is irreducible over Q '
          'and Q[c]/(c^3+c^2-2c-1) is a field')

    # A2  c_g = 2cos(2 pi g/7) from the Chebyshev recurrence, checked closed at g = 7
    cs = [rsc(RONE, 2), RC]
    for g in range(2, 8):
        cs.append(rsub(rmul(RC, cs[g - 1]), cs[g - 2]))
    check('heptagon-cos-recurrence', cs[7] == cs[0] and all(cs[g] == cs[7 - g] for g in (1, 2, 3)),
          'c_g := 2cos(2*pi*g/7) by c_{g+1}=c*c_g-c_{g-1}: c_7=c_0=2 and c_g=c_{7-g} for g=1,2,3')

    # A3  the squared chords s_g = 2 - c_g, and their minimal polynomial
    s = [None] + [rsub(rsc(RONE, 2), cs[g]) for g in (1, 2, 3)]
    def poly_at(x):
        x2 = rmul(x, x)
        return radd(radd(rmul(x2, x), rsc(x2, -7)), radd(rsc(x, 14), rsc(RONE, -7)))
    check('heptagon-chord-minpoly', all(poly_at(s[g]) == Z3 for g in (1, 2, 3)),
          's_g = |1-zeta^g|^2 = 2-2cos(2*pi*g/7) satisfies x^3-7x^2+14x-7 = 0 for g = 1,2,3')
    check('heptagon-chords-symmetric',
          radd(radd(s[1], s[2]), s[3]) == rsc(RONE, 7)
          and radd(radd(rmul(s[1], s[2]), rmul(s[1], s[3])), rmul(s[2], s[3])) == rsc(RONE, 14)
          and rmul(rmul(s[1], s[2]), s[3]) == rsc(RONE, 7),
          'e_1 = 7, e_2 = 14, e_3 = 7, so {s_1,s_2,s_3} is exactly the root set of that cubic')
    check('heptagon-chords-distinct', len({s[1], s[2], s[3]}) == 3,
          's_1 = 2-c, s_2 = 4-c^2, s_3 = 1+c+c^2 are three distinct elements of the field, '
          'so gap g |-> s_g is injective and a triangle class is a multiset of gaps')
    # positivity of the three roots, from the signs of the coefficients alone: for x <= 0 a
    # term a_k x^k is <= 0 when a_k > 0 and k is odd, or a_k < 0 and k is even.  If every
    # term is <= 0 and the constant term is < 0 the cubic has no root x <= 0.
    coeffs = [-7, 14, -7, 1]                     # x^3 - 7x^2 + 14x - 7, low -> high
    nonpos = all((c > 0) if (k % 2) else (c < 0) for k, c in enumerate(coeffs))
    check('heptagon-chords-positive', nonpos and coeffs[0] < 0,
          'x^3-7x^2+14x-7 has coefficient signs (-,+,-,+) from the constant up, so every '
          'term is <= 0 for x <= 0 and the constant is -7 < 0: no root is <= 0, hence '
          's_1, s_2, s_3 > 0')

    # A4  the 35 triples of the heptagon fall into exactly four congruence classes
    def gap(i, j):
        d = (i - j) % 7
        return min(d, 7 - d)
    classes = {}
    for t in itertools.combinations(range(7), 3):
        key = tuple(sorted(gap(a, b) for a, b in itertools.combinations(t, 2)))
        classes.setdefault(key, []).append(t)
    check('heptagon-four-classes', len(classes) == 4 and sum(len(v) for v in classes.values()) == 35,
          '35 triples, %d congruence classes: %s'
          % (len(classes), ' '.join(''.join(map(str, k)) for k in sorted(classes))))
    check('heptagon-classes-are-partitions',
          sorted(classes) == [(1, 1, 2), (1, 2, 3), (1, 3, 3), (2, 2, 3)],
          'the four classes are the arc-gap partitions 5+1+1, 4+2+1, 3+3+1, 3+2+2 of 7, '
          'carrying {s1,s1,s2}, {s1,s2,s3}, {s1,s3,s3}, {s2,s2,s3}')
    check('heptagon-three-distances', len({s[1], s[2], s[3]}) == 3 and len(set(
              gap(a, b) for a, b in itertools.combinations(range(7), 2))) == 3,
          'the heptagon realises exactly three distinct distances, so D = 3')

    # A5  none of the four classes is degenerate: 16*area^2 = 2(ab+bc+ca)-(a^2+b^2+c^2) != 0
    def sixteen_area_sq(tri):
        a, b, cc = (s[g] for g in tri)
        two = radd(radd(rmul(a, b), rmul(b, cc)), rmul(a, cc))
        sq = radd(radd(rmul(a, a), rmul(b, b)), rmul(cc, cc))
        return rsub(rsc(two, 2), sq)
    areas = {k: sixteen_area_sq(k) for k in sorted(classes)}
    check('heptagon-no-degenerate-class', all(v != Z3 for v in areas.values()),
          '16*area^2 = 2(ab+bc+ca)-(a^2+b^2+c^2) is a nonzero field element for all four '
          'classes, so no triple of heptagon vertices is collinear')
    check('heptagon-count-convention-free', len(classes) == 4,
          'hence the heptagon spans 4 classes whether or not collinear triples are counted')
    return s, cs


# ---------------------------------------------------------------------------
# B.  THE HEXAGON PLUS ITS CENTRE (paper, Section 2 Remark), over Q
# ---------------------------------------------------------------------------
def section_B():
    print('--- B.  the hexagon plus its centre (paper, Remark 2.2) ---')
    # coordinates as printed: unit hexagon (x, y*sqrt(3)) with y in (1/2)Z, plus the origin.
    # every squared distance is an integer, so this needs no field extension at all.
    pts = [(Fraction(1), Fraction(0)), (Fraction(1, 2), Fraction(1, 2)),
           (Fraction(-1, 2), Fraction(1, 2)), (Fraction(-1), Fraction(0)),
           (Fraction(-1, 2), Fraction(-1, 2)), (Fraction(1, 2), Fraction(-1, 2)),
           (Fraction(0), Fraction(0))]
    # second coordinate is in units of sqrt(3): |P-Q|^2 = dx^2 + 3*dy^2
    d2 = {}
    for i, j in itertools.combinations(range(7), 2):
        dx = pts[i][0] - pts[j][0]
        dy = pts[i][1] - pts[j][1]
        d2[(i, j)] = dx * dx + 3 * dy * dy
    vals = sorted(set(d2.values()))
    check('hexcentre-squared-distances', vals == [Fraction(1), Fraction(3), Fraction(4)],
          'the 21 squared distances take exactly the three integer values %s'
          % ' '.join(str(v) for v in vals))
    cls = set()
    collinear_triples = 0
    for t in itertools.combinations(range(7), 3):
        m = tuple(sorted(d2[(a, b)] for a, b in itertools.combinations(t, 2)))
        cls.add(m)
        a, b, c = m
        if 2 * (a * b + b * c + a * c) - (a * a + b * b + c * c) == 0:
            collinear_triples += 1
    deg = set()
    for m in cls:
        a, b, c = m
        if 2 * (a * b + b * c + a * c) - (a * a + b * b + c * c) == 0:
            deg.add(m)
    check('hexcentre-five-classes', len(cls) == 5,
          'counting collinear triples, as the target paper does: %d classes, %s'
          % (len(cls), ' '.join('{%s}' % ','.join(str(x) for x in sorted(m)) for m in sorted(cls))))
    check('hexcentre-one-degenerate', deg == {(Fraction(1), Fraction(1), Fraction(4))},
          'exactly one class is collinear, {1,1,4} = {1,1,2} in distances (a vertex, the '
          'centre, the opposite vertex)')
    check('hexcentre-four-noncollinear', len(cls) - len(deg) == 4,
          'so under the noncollinear convention it spans 4 classes on 7 points: the rival '
          'reading of the conjecture is FALSE, which is what the Scope paragraph records')
    check('hexcentre-not-a-heptagon', collinear_triples == 3,
          'it contains %d collinear triples (the three diameters) while the regular heptagon '
          'contains 0, so the two configurations are not similar' % collinear_triples)


# ---------------------------------------------------------------------------
# C/D/E/G.  the census machinery, for a general n
# ---------------------------------------------------------------------------
def edges_of(n):
    E = [(i, j) for j in range(1, n) for i in range(j)]
    EI = {}
    for k, (i, j) in enumerate(E):
        EI[(i, j)] = k
        EI[(j, i)] = k
    return E, EI


def cm_poly(cols6):
    """the 5x5 Cayley-Menger determinant of a 4-point subset whose six squared distances
    carry colours cols6, ordered (01),(02),(03),(12),(13),(23).  It vanishes iff the four
    points are coplanar, so it must vanish for every 4-subset of a planar set."""
    pairs = list(itertools.combinations(range(4), 2))
    M = [[{} for _ in range(5)] for _ in range(5)]
    for j in range(1, 5):
        M[0][j] = dict(ONE)
        M[j][0] = dict(ONE)
    for p, c in zip(pairs, cols6):
        M[p[0] + 1][p[1] + 1] = pvar(c)
        M[p[1] + 1][p[0] + 1] = pvar(c)
    tot = {}
    for perm in itertools.permutations(range(5)):
        sgn = 1
        for i in range(5):
            for j in range(i + 1, 5):
                if perm[i] > perm[j]:
                    sgn = -sgn
        term = dict(ONE)
        zero = False
        for i in range(5):
            e = M[i][perm[i]]
            if not e:
                zero = True
                break
            term = pmul(term, e)
        if not zero:
            tot = padd(tot, {e: c * sgn for e, c in term.items()})
    return tot


_CM_KIND = {}


def cm_kind(cols6):
    """0 = no information, 1 = a nonzero CONSTANT determinant, 2 = a nonconstant
    determinant of one sign on the positive orthant.  1 and 2 both refute the colouring."""
    v = _CM_KIND.get(cols6)
    if v is None:
        p = cm_poly(cols6)
        if not p:
            v = 0
        elif len(p) == 1 and next(iter(p)) == (0, 0, 0):
            v = 1
        elif one_signed(p):
            v = 2
        else:
            v = 0
        _CM_KIND[cols6] = v
    return v


def gram_minors(word, n, EI):
    """every 3x3 minor of the matrix 2G^{(b)}, G^{(b)}_{xy} = <x-b, y-b>, over all base
    points b.  A planar realization has rank(G) <= 2, so all of them vanish."""
    out = []
    for b in range(n):
        others = [v for v in range(n) if v != b]
        M = {}
        for x in others:
            for y in others:
                M[(x, y)] = padd(pvar(word[EI[(b, x)]]), pvar(word[EI[(b, y)]]),
                                 pneg(pvar(word[EI[(x, y)]])) if x != y else {})
        for R in itertools.combinations(others, 3):
            for C in itertools.combinations(others, 3):
                d = det3([[M[(r, c)] for c in C] for r in R])
                if d:
                    out.append(d)
    return out


def restricted_growth(seq):
    m = {}
    out = []
    for c in seq:
        if c not in m:
            m[c] = len(m) + 1
        out.append(m[c])
    return tuple(out)


def act(word, perm, E, EI):
    out = [0] * len(E)
    for k, (i, j) in enumerate(E):
        out[EI[(perm[i], perm[j])]] = word[k]
    return restricted_growth(out)


def orbits(words, n, E, EI):
    """S_n orbits of a set of restricted-growth words, by union-find over two generators
    of S_n (a transposition and a full cycle), which is exact and avoids enumerating n!."""
    idx = {w: i for i, w in enumerate(sorted(words))}
    par = list(range(len(idx)))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    gens = [[1, 0] + list(range(2, n)), list(range(1, n)) + [0]]
    closed = True
    for w in idx:
        for g in gens:
            im = act(w, g, E, EI)
            if im not in idx:
                closed = False
                continue
            a, b = find(idx[w]), find(idx[im])
            if a != b:
                par[a] = b
    inv = {i: w for w, i in idx.items()}
    groups = {}
    for i in range(len(idx)):
        groups.setdefault(find(i), []).append(inv[i])
    reps = sorted((min(v), len(v)) for v in groups.values())
    return reps, closed


def census_k7():
    """Lemma 4.1: every restricted-growth colouring of the 21 edges of K_7 with at least
    three colours and at most four distinct triple colour-multisets."""
    n = 7
    E, EI = edges_of(n)
    NE = len(E)
    completed = []
    for k, (b, c) in enumerate(E):
        completed.append(tuple((EI[(a, b)], EI[(a, c)], k) for a in range(b)))
    leaves = []
    nodes = 0
    sys.setrecursionlimit(10000)

    def dfs(k, col, ncol, msets):
        nonlocal nodes
        nodes += 1
        if k == NE:
            if ncol >= 3:
                leaves.append((tuple(col), ncol, len(msets)))
            return
        for c in range(1, ncol + 2):
            col[k] = c
            new = msets
            bad = False
            for e1, e2, e3 in completed[k]:
                m = tuple(sorted((col[e1], col[e2], c)))
                if m not in new:
                    if len(new) >= 4:
                        bad = True
                        break
                    if new is msets:
                        new = set(msets)
                    new.add(m)
            if not bad:
                dfs(k + 1, col, max(ncol, c), new)
            col[k] = 0
    dfs(0, [0] * NE, 0, set())
    return leaves, nodes, E, EI


def section_CD(leaves, nodes, E, EI):
    print('--- C.  the census of K_7 colourings in the cell (paper, Lemma 4.1) ---')
    check('census-dfs-nodes', nodes == 2915582, 'the search tree has %d nodes; the census is '
          'complete and uncapped (no node budget, no time limit)' % nodes)
    check('census-leaves', len(leaves) == 161351,
          '%d restricted-growth colourings with D >= 3 colours and at most 4 triple '
          'colour-multisets' % len(leaves))
    byD = {}
    cross = {}
    for w, nc, nm in leaves:
        byD[nc] = byD.get(nc, 0) + 1
        cross[(nc, nm)] = cross.get((nc, nm), 0) + 1
    check('census-max-distances', max(byD) == 4 and sorted(byD) == [3, 4],
          'D <= 4 everywhere in the cell: D=3 in %d leaves, D=4 in %d' % (byD[3], byD[4]))
    check('census-cross-tabulation',
          cross == {(3, 3): 3948, (3, 4): 152762, (4, 4): 4641},
          '(D, #classes): (3,3) %d, (3,4) %d, (4,4) %d; and 3948+152762+4641 = %d'
          % (cross[(3, 3)], cross[(3, 4)], cross[(4, 4)], sum(cross.values())))
    hep = restricted_growth([min((j - i) % 7, (i - j) % 7) for (i, j) in E])
    words = set(w for w, _, _ in leaves)
    check('census-contains-heptagon', hep in words,
          "the heptagon's own gap colouring %s is one of the leaves"
          % ''.join(map(str, hep)))

    print('--- D.  the Cayley-Menger filter and the S_7 orbit count (paper, Lemma 4.2) ---')
    quads = list(itertools.combinations(range(7), 4))
    QE = [tuple(EI[(a, b)] for a, b in itertools.combinations(q, 2)) for q in quads]
    check('census-quad-count', len(quads) == 35, 'C(7,4) = %d four-point subsets per colouring'
          % len(quads))
    kills_const = 0
    kills_signed = 0
    surv = []
    for w, nc, nm in leaves:
        worst = 0
        for qe in QE:
            k = cm_kind((w[qe[0]], w[qe[1]], w[qe[2]], w[qe[3]], w[qe[4]], w[qe[5]]))
            if k > worst:
                worst = k
            if worst == 1:
                break
        if worst == 1:
            kills_const += 1
        elif worst == 2:
            kills_signed += 1
        else:
            surv.append((w, nc, nm))
    check('cm-kill-nonzero-constant', kills_const == 81306,
          '%d colourings carry a 4-subset whose Cayley-Menger determinant is a nonzero '
          'CONSTANT, so it cannot vanish for any squared distances' % kills_const)
    check('cm-kill-one-signed', kills_signed == 70286,
          '%d more carry one whose determinant is nonconstant with all coefficients of one '
          'sign, so it cannot vanish at any POSITIVE squared distances' % kills_signed)
    check('cm-survivors', len(surv) == 9759, '%d colourings survive' % len(surv))
    check('cm-exhaustion-arithmetic', kills_const + kills_signed + len(surv) == len(leaves),
          '%d + %d + %d = %d, the whole census accounted for'
          % (kills_const, kills_signed, len(surv), len(leaves)))
    reps, closed = orbits(set(w for w, _, _ in surv), 7, E, EI)
    check('cm-survivors-S7-closed', closed,
          'the survivor set is closed under the S_7 action, as it must be: the filter reads '
          'the colouring only through its 4-subsets')
    prof = {}
    for _, sz in reps:
        prof[sz] = prof.get(sz, 0) + 1
    check('cm-orbit-count', len(reps) == 30, 'the 9759 survivors form %d orbits under S_7'
          % len(reps))
    check('cm-orbit-size-profile',
          prof == {70: 1, 105: 5, 120: 1, 140: 1, 210: 5, 252: 2, 315: 8, 420: 1, 630: 5,
                   1260: 1},
          'orbit sizes with multiplicities: %s'
          % ' '.join('%d^%d' % (k, v) for k, v in sorted(prof.items())))
    check('cm-orbit-sizes-sum', sum(sz for _, sz in reps) == 9759,
          'the orbit sizes sum to %d' % sum(sz for _, sz in reps))
    return reps, hep


def section_E(reps, E, EI):
    print('--- E.  exact elimination of 29 of the 30 orbits (paper, Lemma 4.3) ---')
    killed_signed = []
    killed_zero = []
    open_orbits = []
    for w, sz in reps:
        nv = max(w) - 1
        ms = [p for p in gram_minors(w, 7, EI) if p]
        if any(one_signed(p) for p in ms):
            killed_signed.append((w, sz))
            continue
        # forced-value loop: a minor that involves a single unknown y is a univariate
        # polynomial that must vanish, so y is a root of the gcd of all such minors.
        forced = {}
        killer = None
        trace = []
        for _ in range(nv):
            progress = False
            for vi in range(nv):
                if vi in forced:
                    continue
                U = [p for p in ms
                     if any(e[vi] > 0 for e in p)
                     and all(all(k == 0 for jj, k in enumerate(e) if jj != vi) for e in p)]
                if not U:
                    continue
                g = None
                for p in U:
                    d = max(e[vi] for e in p)
                    lst = [Fraction(0)] * (d + 1)
                    for e, c in p.items():
                        lst[e[vi]] += c
                    g = lst if g is None else ugcd(g, lst)
                k = 0
                while k < len(g) and g[k] == 0:
                    k += 1
                h = g[k:]
                trace.append('gcd of the %d minors in y_%d alone = %s'
                             % (len(U), vi + 2, poly_str(g, 'y_%d' % (vi + 2))))
                if len(h) == 1:
                    if k == 0:
                        killer = 'that gcd is a nonzero CONSTANT, so those minors have no ' \
                                 'common root at all'
                    else:
                        killer = 'that gcd is a nonzero multiple of y_%d^%d, forcing ' \
                                 'y_%d = 0, which no two distinct points admit' \
                                 % (vi + 2, k, vi + 2)
                    break
                roots = positive_rational_roots(h)
                if len(h) - 1 == 1 and len(roots) == 1:
                    val = roots[0]
                    forced[vi] = val
                    new = []
                    for p in ms:
                        r = {}
                        for e, c in p.items():
                            ne = list(e)
                            k2 = ne[vi]
                            ne[vi] = 0
                            r[tuple(ne)] = r.get(tuple(ne), 0) + Fraction(c) * val ** k2
                        r = {e: c for e, c in r.items() if c}
                        if r:
                            new.append(r)
                    ms = new
                    progress = True
            if killer or not progress:
                break
        if killer:
            killed_zero.append((w, sz, forced, trace, killer))
        else:
            open_orbits.append((w, sz, forced))
    for w, sz, forced, trace, why in killed_zero:
        print('     orbit %s (size %d), D = %d:' % (''.join(map(str, w)), sz, max(w)))
        for line in trace:
            print('        %s' % line)
        print('        => forced %s; %s'
              % (', '.join('y_%d = %s' % (k + 2, v) for k, v in sorted(forced.items())) or '-',
                 why))
    check('elim-one-signed-gram-minor', len(killed_signed) == 27,
          '%d of the 30 orbits carry a 3x3 Gram minor that is a nonzero one-signed '
          'polynomial, hence cannot vanish at positive squared distances' % len(killed_signed))
    check('elim-univariate-gcd', len(killed_zero) == 2,
          '%d further orbits die on the univariate gcd of their Gram minors, as traced above'
          % len(killed_zero))
    check('elim-single-survivor', len(open_orbits) == 1,
          '%d orbit survives both exact tests: %s of size %d'
          % (len(open_orbits), ''.join(map(str, open_orbits[0][0])), open_orbits[0][1]))
    check('elim-accounting', len(killed_signed) + len(killed_zero) + len(open_orbits) == 30,
          '27 + 2 + 1 = 30')
    return open_orbits[0][0]


def section_F(survivor, s, cs, E, EI, hep):
    print('--- F.  the survivor is the gap colouring, and it is the heptagon '
          '(paper, Lemma 4.4) ---')
    # F1  the survivor equals, up to relabelling the vertices, the cyclic gap colouring
    good = [p for p in itertools.permutations(range(7)) if act(survivor, p, E, EI) == hep]
    check('survivor-is-gap-colouring', bool(good),
          'a vertex relabelling carries the surviving word %s to the gap colouring %s; there '
          'are %d such relabellings, so its stabiliser has order %d = |AGL(1,7)|'
          % (''.join(map(str, survivor)), ''.join(map(str, hep)), len(good), len(good)))
    check('survivor-orbit-size', len(good) * 120 == 5040,
          'consistent with the orbit size 120 found in D: 5040/42 = 120')
    # colour == cyclic gap on all 21 edges, and each colour class is a single 7-cycle
    ok = all(hep[EI[(i, j)]] == min((i - j) % 7, (j - i) % 7) for i, j in E)
    cyc = True
    sizes = []
    for g in (1, 2, 3):
        cl = [(i, j) for (i, j) in E if hep[EI[(i, j)]] == g]
        sizes.append(len(cl))
        adj = {v: [] for v in range(7)}
        for i, j in cl:
            adj[i].append(j)
            adj[j].append(i)
        if len(cl) != 7 or any(len(adj[v]) != 2 for v in range(7)):
            cyc = False
            continue
        seen = {0}
        cur, prev = adj[0][0], 0
        while cur != 0:
            seen.add(cur)
            nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
            cur, prev = nxt, cur
        if len(seen) != 7:
            cyc = False
    check('gap-colouring-explicit', ok and cyc and sizes == [7, 7, 7],
          'colour g is exactly cyclic gap g on all 21 edges, and each of the three colour '
          'classes is a single 7-cycle (sizes %s)' % sizes)
    # F2  Schoenberg for a circulant: mu_m = sum_g y_g * c_{gm}, m = 1,2,3
    def idx(k):
        k %= 7
        return min(k, 7 - k)
    def mu(y, m):
        t = Z3
        for g in (1, 2, 3):
            t = radd(t, rmul(y[g], cs[idx(g * m)]))
        return t
    y = {1: s[1], 2: s[2], 3: s[3]}
    vals = [mu(y, m) for m in (1, 2, 3)]
    check('fourier-heptagon-spectrum',
          vals == [rsc(RONE, -7), Z3, Z3],
          'for y = (s_1,s_2,s_3): (mu_1, mu_2, mu_3) = (-7, 0, 0) are the eigenvalues of M '
          'on the sum-zero hyperplane, the value -7 with multiplicity 2, so M is negative '
          'semidefinite there of rank exactly 2 and equivalently -M is positive semidefinite '
          'of rank exactly 2, which is what Schoenberg requires for R^2')
    # F3  rank of the two linear conditions mu_2 = mu_3 = 0 is 2, so the ray is unique
    rows = [[cs[idx(g * m)] for g in (1, 2, 3)] for m in (2, 3)]
    minor = rsub(rmul(rows[0][0], rows[1][1]), rmul(rows[0][1], rows[1][0]))
    check('fourier-rank-two', minor != Z3,
          'the 2x3 coefficient matrix of (mu_2, mu_3) has the nonzero 2x2 minor '
          'c_2*c_1 - c_3*c_3, hence rank 2 over the field, hence a 1-dimensional kernel')
    check('fourier-kernel-is-the-heptagon', vals[1] == Z3 and vals[2] == Z3 and minor != Z3,
          '(s_1,s_2,s_3) lies in that 1-dimensional kernel, so it SPANS it: the only '
          'circulant solution with mu_1 < 0 is the regular heptagon up to scale')
    # F4  the other two sign choices are the same point set, relabelled by a multiplier
    y2 = {1: s[2], 2: s[3], 3: s[1]}
    y3 = {1: s[3], 2: s[1], 3: s[2]}
    v2 = [mu(y2, m) for m in (1, 2, 3)]
    v3 = [mu(y3, m) for m in (1, 2, 3)]
    check('fourier-case-mu2-negative', v2 == [Z3, rsc(RONE, -7), Z3],
          'y = (s_2,s_3,s_1) gives (0,-7,0): this is the heptagon with its distance classes '
          'relabelled by the multiplier g -> 2g on Z_7')
    check('fourier-case-mu3-negative', v3 == [Z3, Z3, rsc(RONE, -7)],
          'y = (s_3,s_1,s_2) gives (0,0,-7): the multiplier g -> 3g.  All three sign choices '
          'therefore give the SAME point set')
    check('fourier-multiplier-orbit',
          [min((2 * g) % 7, 7 - (2 * g) % 7) for g in (1, 2, 3)] == [2, 3, 1]
          and [min((3 * g) % 7, 7 - (3 * g) % 7) for g in (1, 2, 3)] == [3, 1, 2],
          'g -> 2g permutes the gaps (1,2,3) cyclically as (2,3,1) and g -> 3g as (3,1,2), '
          'which is exactly the permutation of the rays above')


def section_G(survivor, leaves):
    print('--- G.  no 7-point planar set spans at most three classes (paper, Lemma 4.5) ---')
    nm = None
    for w, nc, m in leaves:
        if w == survivor:
            nm = m
    check('survivor-has-four-classes', nm == 4,
          'the unique surviving colouring realises 4 distinct triple multisets, not 3, so no '
          '7-point planar set with at least 3 distances spans at most 3 classes')
    print('--- G2.  the lanes D = 1 and D = 2 (paper, Lemma 4.5) ---')
    check('one-distance-lane-empty', cm_kind((1, 1, 1, 1, 1, 1)) == 1,
          'the Cayley-Menger determinant of a monochromatic 4-subset is a nonzero constant, '
          'so no 4 points in the plane are equidistant and D = 1 is impossible for n >= 4')
    print('--- G2.  the two-distance lane at n = 6 (paper, Lemma 4.5, second half) ---')
    n = 6
    E6, EI6 = edges_of(n)
    NE = len(E6)
    quads = list(itertools.combinations(range(n), 4))
    QE = [tuple(EI6[(a, b)] for a, b in itertools.combinations(q, 2)) for q in quads]
    words = []
    for mask in range(1 << (NE - 1)):
        w = tuple([1] + [1 + ((mask >> b) & 1) for b in range(NE - 1)])
        if max(w) == 1:
            continue
        words.append(w)
    check('n6-two-colour-word-count', len(words) == 2 ** (NE - 1) - 1 == 16383,
          'K_6 has 15 edges, so there are 2^14 - 1 = %d restricted-growth colourings using '
          'exactly two colours, and all of them lie in the cell' % len(words))
    alive = [w for w in words
             if not any(cm_kind(tuple(w[x] for x in qe)) for qe in QE)]
    reps6, closed6 = orbits(set(alive), n, E6, EI6)
    check('n6-orbit-count', closed6 and len(reps6) == 42,
          '%d survive the Cayley-Menger filter and form %d orbits under S_6'
          % (len(alive), len(reps6)))
    k1 = k2 = 0
    stuck = []
    for w, sz in reps6:
        ms = [p for p in gram_minors(w, n, EI6) if p]
        if any(one_signed(p) for p in ms):
            k1 += 1
            continue
        U = [p for p in ms if any(e[0] > 0 for e in p)]
        g = None
        for p in U:
            d = max(e[0] for e in p)
            lst = [Fraction(0)] * (d + 1)
            for e, c in p.items():
                lst[e[0]] += c
            g = lst if g is None else ugcd(g, lst)
        if g is None:
            stuck.append(w)
            continue
        k = 0
        while k < len(g) and g[k] == 0:
            k += 1
        if len(g[k:]) == 1:
            k2 += 1
        else:
            stuck.append(w)
    check('n6-two-distance-lane-empty', k1 + k2 == len(reps6) and not stuck,
          '%d orbits die on a one-signed Gram minor and %d are forced to a vanishing '
          'distance: NO 6-point planar two-distance set exists (Erdos-Fishburn g(2) = 5), '
          'so at n = 7 the lane D <= 2 is empty by heredity' % (k1, k2))


def main():
    s, cs = section_A()
    section_B()
    leaves, nodes, E, EI = census_k7()
    reps, hep = section_CD(leaves, nodes, E, EI)
    survivor = section_E(reps, E, EI)
    section_F(survivor, s, cs, E, EI, hep)
    section_G(survivor, leaves)

    print()
    print('NOT RE-RUN: Schoenberg\'s criterion itself (a point set with squared-distance '
          'matrix M embeds in R^d iff -M is positive semidefinite of rank <= d on the '
          'sum-zero hyperplane) is quoted from the literature, not proved here.')
    print('NOT RE-RUN: part (c) of Theorem 1.1, the step from seven points to eight.  It '
          'is a hand argument, Proposition 5.1 of the paper, and contains no computation.')
    print('NOT RE-RUN: any claim about the NONCOLLINEAR convention beyond the single '
          'hexagon-plus-centre configuration checked in B.  In particular this program does '
          'NOT decide the maximality half under that convention, and the paper does not '
          'claim it.')
    print('NOT RE-RUN: minimality or uniqueness of the certificates used above.  The '
          'one-signed test and the forced-zero test are SUFFICIENT conditions for '
          'non-realizability; nothing here claims they are the shortest such certificates, '
          'or that the 30 orbits could not be eliminated by a different route.')
    print('NOT RE-RUN: the published classifications this result is adjacent to '
          '(Shinohara 2004 on planar three-distance sets, Erdos-Fishburn 1996 on g(k)).  '
          'The two-distance lane is re-derived above; the three-distance classification is '
          'NOT used anywhere in this proof and is NOT re-derived.')
    print()
    if _N_FAIL:
        print('VERDICT: %d CHECK(S) FAILED of %d' % (_N_FAIL, _N_FAIL + _N_PASS))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _N_PASS)
    return 0


if __name__ == '__main__':
    sys.exit(main())

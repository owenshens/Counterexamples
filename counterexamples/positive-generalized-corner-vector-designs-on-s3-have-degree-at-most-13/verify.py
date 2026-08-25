#!/usr/bin/env python3
"""
verify.py -- referee verification for
  "Positive Generalized Corner-Vector Designs on S^3 Have Degree at Most 13"
  (refutation of Problem 6.1 of Tanino-Tamaru-Hirao-Sawa, Algebr. Comb. 8 (2025)).

Standard library only.  All decisions use exact integer / Fraction arithmetic;
no floating point is used anywhere in a comparison.

----------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER  (inputs, never checked against themselves)
----------------------------------------------------------------------
  * PAPER_COEF : the exhibited object.  The 11 coefficients of
        p = sum c_{ijk} e_2^i e_3^j e_4^k,   e_j = e_j(x_1^2,...,x_4^2),
    given by eq. (2) of the paper.  This is the ONLY mathematical input.
  * PAPER_MEAN, PAPER_H, PAPER_SUBDIV, PAPER_BMIN, PAPER_SCALE :
    the numbers the paper *asserts* about that object (Lemma 3's spherical
    mean, the three reduced polynomials H_1,H_2,H_3 with their integrality
    scalings 1, 32, 729, the three subdivisions of [0,1], and the nine
    minimum Bernstein coefficients of the table).  These are recomputed
    from PAPER_COEF and compared; each such comparison can fail.

----------------------------------------------------------------------
WHAT IS DERIVED HERE  (the checks)
----------------------------------------------------------------------
  1. p is decoded, expanded into 130 monomials in y_i = x_i^2, printed back,
     shown symmetric in y (hence B_4-invariant) and of degree exactly 14 in x.
  2. The Dirichlet(1/2,1/2,1/2,1/2) moment functional is built and validated
     independently: <e_1^n> = 1 for n = 0..7; every odd-degree moment vanishes,
     established by the COMPUTED sign-subgroup symmetrization of the monomial
     (with an even-degree control where that symmetrization is 16, not 0);
     agreement with the 24-cell (an exact 5-design on S^3) through degree 5 --
     together with a control showing the 24-cell is NOT exact at degree 6; and
     agreement with the iterated-Laplacian formula on all 330 even monomials of
     x-degree <= 14, i.e. over the whole range p occupies, by an
     algorithmically unrelated route.  The Pochhammer / (2a-1)!! pair is also
     compared, but that pair is an algebraic REWRITING of one formula, so it is
     labelled in the output as a transcription check, not as independence.
  3. <p> is recomputed from the expansion and compared with PAPER_MEAN, and
     its sign is computed.
  4. The three orbit profiles (u,c,...,c,0,...,0), c=(1-u)/s, are substituted
     into p; the resulting univariate polynomials, their minimal integer
     scalings, and their nine minimum Bernstein coefficients over the paper's
     subdivisions are recomputed and compared with the table.
  5. Positivity of H_1,H_2,H_3 on [0,1] is certified a SECOND, independent way
     by exact Sturm root counting (with a control on polynomials of known
     root count), and a third way by exact evaluation of p at 609 genuine
     orbit points with rational a^2.
  6. B_4 orbits of representative corner vectors are enumerated explicitly
     (384 signed permutations); orbit sizes are computed, p is verified
     constant on each orbit, and every degree-15 monomial is verified to have
     zero orbit sum and (via its computed sign-subgroup symmetrization, with a
     degree-14 control where that symmetrization is nonzero) zero spherical
     mean -- this is the paper's Remark 5: exactness through degree 14 ==
     exactness through degree 15.
  7. eq. (3) itself is verified rather than trusted: for all nine intervals
     sum_k beta_k C(d,k) z^k (1-z)^{d-k} is checked to equal H(r+(q-r)z)
     coefficient by coefficient, with the affine shift pinned by 9-point
     agreement and a control in which perturbing one beta_k breaks the
     identity.  Only then does the bound H >= min_k beta_k carry any weight.
     The two facts that turn min_k beta_k into a bound are themselves computed
     rather than asserted: (a) each basis element is verified, coefficient by
     coefficient, to be comb(d,k) z^k (1-z)^{d-k} with comb(d,k) > 0, and the
     basis is verified to sum to the constant polynomial 1 (partition of
     unity), and (b) for each of the nine intervals the residual
     H(r+(q-r)z) - min_k beta_k is verified to equal sum_k (beta_k - min beta)
     C(d,k) z^k (1-z)^{d-k} exactly, with every coefficient beta_k - min beta
     computed to be >= 0.  The only step left outside the arithmetic is then
     the sign fact z >= 0 and 1-z >= 0 for z in [0,1], which is the definition
     of the interval.  The bound property is additionally spot-checked at 189
     points spread over the nine intervals, and the resulting global bound is
     checked to lie below every value of p actually evaluated.
     The quantifier over a is discharged on the u side as well: u = a^2/(a^2+s)
     is verified to be strictly inside (0,1) and inside one of the paper's
     subintervals at every sampled a, to be strictly increasing in a^2, and to
     be onto: for 65 rational targets u* spread over (0,1) the preimage
     a^2 = s u*/(1-u*) is formed and p is re-evaluated at that genuine orbit
     point and compared with H_s(u*)/scale_s and with the interval's Bernstein
     minimum.
  8. The load-bearing conclusion is computed, not asserted: a strictly
     positive rational lower bound L for p over EVERY generalized
     corner-vector point v_{a,s} (all a>0, s in {0,1,2,3}) is assembled from
     the certified Bernstein minima and from the COMPUTED value of p at
     (1,0,0,0) for the s = 0 branch (no literal from the paper enters L), and
     L > 0 > <p> is checked.  Any
     normalized positive-weight formula therefore returns >= L for p while
     the true integral is negative.  The contradiction is also exhibited
     numerically on 200 random normalized positive-weight configurations.

There is no exhaustive census in this paper: the "for all a>0" quantifier is a
continuum, and it is discharged exactly by the Bernstein and Sturm
certificates, so the verification is complete rather than sampled.

Falsifiability: no mutation log, corrupted-run transcript or second
implementation accompanies this program, and none should be assumed.  What a
referee can check from this file alone is that the falsifiability is built in
and re-run on every invocation: several checks carry an explicit control whose
value is computed and printed beside the check it protects -- a polynomial with
a negative Bernstein coefficient, a beta_k perturbed by 1 that breaks the
eq. (3) identity, a perturbed residual coefficient that breaks the nonnegative
representation, a weight moved off 1 that destroys the partition of unity,
polynomials of known Sturm root count, the 24-cell failing at degree 6, and
even-degree monomials whose sign-subgroup symmetrization is nonzero.  A
control that stopped discriminating would itself turn its check into a FAIL.
Every comparison against a number the paper asserts is likewise a comparison
with a quantity recomputed here from the object of eq. (2) alone, so a wrong
paper constant fails rather than being absorbed.

Output contract: one "PASS <name>" / "FAIL <name>" line per check, then a NOTE
line and a "NOT RE-RUN:" line stating what this program does not cover, then
"VERDICT: ALL <n> CHECKS PASS" or "VERDICT: <k> OF <n> CHECKS FAILED".
Exit status 0 iff every check passed.
"""

import sys
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, gcd

_RESULTS = []


def check(name, ok, detail=""):
    """Record one check.  ok must be a genuine computed boolean."""
    ok = bool(ok)
    _RESULTS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    return ok


def verdict():
    n = len(_RESULTS)
    bad = [nm for nm, ok in _RESULTS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


# ======================================================================
# INPUTS TAKEN FROM THE PAPER
# ======================================================================

# eq. (2): p = sum_{(i,j,k)} c * e_2^i e_3^j e_4^k.
PAPER_COEF = {
    (0, 0, 0): 2,
    (1, 0, 0): 240,
    (0, 1, 0): 1345,
    (0, 0, 1): 227924,
    (2, 0, 0): -3518,
    (1, 1, 0): -448,
    (0, 2, 0): 283604,
    (1, 0, 1): -1529173,
    (3, 0, 0): 12886,
    (0, 1, 1): 4855151,
    (2, 1, 0): -132972,
}

# The paper's displayed monomial list, in the order printed in eq. (2),
# used only to confirm the decoded tuple is the printed polynomial.
PAPER_TERMS = [
    ("", 2), ("e_2", 240), ("e_3", 1345), ("e_4", 227924),
    ("e_2^2", -3518), ("e_2e_3", -448), ("e_3^2", 283604),
    ("e_2e_4", -1529173), ("e_2^3", 12886), ("e_3e_4", 4855151),
    ("e_2^2e_3", -132972),
]

PAPER_DEGREE_IN_X = 14           # "p ... has degree 14 in x_1,...,x_4"
PAPER_MEAN = F(-60763, 143360)   # Lemma 3
PAPER_P_AT_S0 = 2                # "p(v_{a,0}) = 2"

# Lemma 4: p(v_{a,1}) = H_1(u), 32 p(v_{a,2}) = H_2(u), 729 p(v_{a,3}) = H_3(u).
PAPER_SCALE = {1: 1, 2: 32, 3: 729}
# H_s as coefficient lists, index = power of u.
PAPER_H = {
    1: [2, 240, -3758, 19922, -42176, 38658, -12886],
    2: [1391, -42268, 440357, -1750294, 3540713, -3979288, 2387827, -598374],
    3: [39713, -1819702, 26917938, -130946506, 299081978, -355558770,
        213893327, -51606520],
}
# Table: subdivision 0=r_0<r_1<r_2<r_3=1 and the min Bernstein coefficient
# of H_s on each of the three successive intervals.
PAPER_SUBDIV = {
    1: [F(0), F(1, 6), F(5, 6), F(1)],
    2: [F(0), F(1, 12), F(3, 4), F(1)],
    3: [F(0), F(1, 20), F(3, 5), F(1)],
}
PAPER_BMIN = {
    1: [F(14287, 7776), F(46531, 23328), F(14287, 7776)],
    2: [F(137627845, 2985984), F(414113333, 5971968), F(219941, 5376)],
    3: [F(542097724127, 448000000), F(17317137, 15625), F(4479262, 3125)],
}


# ======================================================================
# univariate polynomials over Q: list of coefficients, index = power
# ======================================================================

def upmul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y == 0:
                continue
            r[i + j] += x * y
    return r


def upadd(a, b):
    r = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        r[i] += x
    for i, x in enumerate(b):
        r[i] += x
    return r


def upscale(a, s):
    return [s * x for x in a]


def uppow(a, n):
    r = [F(1)]
    for _ in range(n):
        r = upmul(r, a)
    return r


def uptrim(a):
    r = list(a)
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    return r


def upeval(a, t):
    s = F(0)
    for c in reversed(a):
        s = s * t + c
    return s


# ======================================================================
# multivariate polynomials in y_1..y_4: dict exponent-4-tuple -> coeff
# ======================================================================

def mpmul(A, B):
    r = {}
    for ea, ca in A.items():
        for eb, cb in B.items():
            k = (ea[0] + eb[0], ea[1] + eb[1], ea[2] + eb[2], ea[3] + eb[3])
            r[k] = r.get(k, 0) + ca * cb
    return {k: v for k, v in r.items() if v != 0}


def mpadd(A, B):
    r = dict(A)
    for k, v in B.items():
        r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v != 0}


def mppow(A, n):
    r = {(0, 0, 0, 0): 1}
    for _ in range(n):
        r = mpmul(r, A)
    return r


def elem_sym_y(j):
    """e_j(y_1,...,y_4) as a multivariate polynomial in y."""
    r = {}
    for cb in combinations(range(4), j):
        e = [0, 0, 0, 0]
        for i in cb:
            e[i] = 1
        r[tuple(e)] = 1
    return r


def mpeval(A, y):
    s = F(0)
    for e, c in A.items():
        t = F(c)
        for yi, ei in zip(y, e):
            if ei:
                t *= yi ** ei
        s += t
    return s


def build_p_in_y(coef):
    """Expand sum c_{ijk} e_2^i e_3^j e_4^k into monomials in y."""
    E = {2: elem_sym_y(2), 3: elem_sym_y(3), 4: elem_sym_y(4)}
    out = {}
    for (i, j, k), c in coef.items():
        t = mpmul(mpmul(mppow(E[2], i), mppow(E[3], j)), mppow(E[4], k))
        out = mpadd(out, {m: c * v for m, v in t.items()})
    return out


# ======================================================================
# the spherical moment functional on S^3, two independent closed forms
# ======================================================================

def poch(z, n):
    r = F(1)
    for i in range(n):
        r *= (z + i)
    return r


def mom_y_pochhammer(alpha):
    """<y_1^a1...y_4^a4> = prod (1/2)_{a_i} / (2)_{sum a_i}   (Dirichlet route)."""
    num = F(1)
    for a in alpha:
        num *= poch(F(1, 2), a)
    return num / poch(F(2), sum(alpha))


def dfact_odd(a):
    """(2a-1)!! ; equals 1 for a = 0."""
    r = 1
    for i in range(1, a + 1):
        r *= (2 * i - 1)
    return r


def mom_y_doublefact(alpha):
    """Same moment via prod (2a_i-1)!! / (4*6*...*(2m+2)),  m = sum a_i.

    This is the classical surface-measure formula for the unit sphere of R^4
    and is coded from a different closed form than mom_y_pochhammer.
    """
    m = sum(alpha)
    num = 1
    for a in alpha:
        num *= dfact_odd(a)
    den = 1
    for k in range(1, m + 1):
        den *= (2 * k + 2)
    return F(num, den)


def mom_x(alpha):
    """<x_1^a1...x_4^a4> on S^3; zero unless every exponent is even."""
    if any(a % 2 for a in alpha):
        return F(0)
    return mom_y_pochhammer(tuple(a // 2 for a in alpha))


def mean_of_y_poly(A):
    return sum(F(c) * mom_y_pochhammer(e) for e, c in A.items())


def cell24_points():
    """Vertices of the 24-cell: +-e_i and (+-1,+-1,+-1,+-1)/2.  Exact rationals.

    These 24 points are a spherical 5-design on S^3 (verified below through
    degree 5, with a control at degree 6 where they are NOT exact).
    """
    pts = []
    for i in range(4):
        for sg in (F(1), F(-1)):
            v = [F(0)] * 4
            v[i] = sg
            pts.append(tuple(v))
    for sg in product((F(1, 2), F(-1, 2)), repeat=4):
        pts.append(sg)
    return pts


def cell24_average(alpha, pts):
    s = F(0)
    for pt in pts:
        t = F(1)
        for c, a in zip(pt, alpha):
            if a:
                t *= c ** a
        s += t
    return s / len(pts)


# ======================================================================
# orbit profiles:  squared coordinates of v_{a,s} are (u, c,...,c, 0,...,0)
# with u = a^2/(a^2+s) and c = (1-u)/s.  Substituting gives H_s(u)/scale.
# ======================================================================

def profile_polys(s):
    """The four squared coordinates of v_{a,s} as univariate polynomials in u."""
    u = [F(0), F(1)]                 # u
    c = [F(1, s), F(-1, s)] if s else None   # (1-u)/s
    vals = [u]
    for _ in range(s):
        vals.append(list(c))
    while len(vals) < 4:
        vals.append([F(0)])
    return vals


def elem_sym_uni(vals, j):
    tot = [F(0)]
    for cb in combinations(range(len(vals)), j):
        t = [F(1)]
        for i in cb:
            t = upmul(t, vals[i])
        tot = upadd(tot, t)
    return tot


def p_on_profile(coef, s):
    """p restricted to the s-profile, as a univariate polynomial in u over Q."""
    vals = profile_polys(s)
    E = {j: elem_sym_uni(vals, j) for j in (2, 3, 4)}
    tot = [F(0)]
    for (i, j, k), cf in coef.items():
        t = upmul(upmul(uppow(E[2], i), uppow(E[3], j)), uppow(E[4], k))
        tot = upadd(tot, upscale(t, F(cf)))
    return uptrim(tot)


def integer_scaling(poly):
    """Least positive integer L with L*poly integral, and the scaled coeffs."""
    L = 1
    for x in poly:
        d = x.denominator
        L = L * d // gcd(L, d)
    return L, [x * L for x in poly]


# ======================================================================
# exact Bernstein coefficients on a subinterval  (eq. (3) of the paper)
# ======================================================================

def affine_shift(coefs, r, q):
    """Coefficients c_j of H(r + (q-r) z) in the monomial basis of z."""
    d = len(coefs) - 1
    out = [F(0)] * (d + 1)
    for j, cj in enumerate(coefs):
        if cj == 0:
            continue
        for t in range(j + 1):
            out[t] += F(cj) * comb(j, t) * (F(r) ** (j - t)) * (F(q - r) ** t)
    return out


def bernstein_coeffs(coefs, r, q):
    """beta_{I,k} = sum_{j<=k} c_j C(k,j)/C(d,j) for I = [r,q]."""
    c = affine_shift(coefs, r, q)
    d = len(coefs) - 1
    return [sum(c[j] * F(comb(k, j), comb(d, j)) for j in range(k + 1))
            for k in range(d + 1)]


def bernstein_reconstruct(beta):
    """sum_k beta_k C(d,k) z^k (1-z)^{d-k}, as coefficients in z.

    Used to VERIFY eq. (3) as an exact polynomial identity rather than trusting
    it: the Bernstein bound H >= min_k beta_k on I is sound only if the beta_k
    really are the Bernstein coefficients of H on I, i.e. only if this
    reconstruction returns H(r + (q-r) z) exactly.
    """
    d = len(beta) - 1
    out = [F(0)] * (d + 1)
    for k, bk in enumerate(beta):
        if bk == 0:
            continue
        w = F(bk) * comb(d, k)
        for i in range(d - k + 1):
            out[k + i] += w * comb(d - k, i) * ((-1) ** i)
    return out


# ======================================================================
# exact Sturm root counting on (a, b]  -- independent positivity certificate
# ======================================================================

def upderiv(a):
    return uptrim([F(i) * a[i] for i in range(1, len(a))] or [F(0)])


def uprem(a, b):
    """Remainder of a mod b over Q."""
    a = list(a)
    db = len(uptrim(b)) - 1
    lb = uptrim(b)[-1]
    while True:
        a = uptrim(a)
        da = len(a) - 1
        if a == [F(0)] or da < db:
            return a
        f = a[da] / lb
        for i in range(db + 1):
            a[da - db + i] -= f * b[i]
        a[da] = F(0)


def sturm_chain(a):
    a = uptrim(a)
    ch = [a, upderiv(a)]
    while uptrim(ch[-1]) != [F(0)]:
        r = uprem(ch[-2], ch[-1])
        if uptrim(r) == [F(0)]:
            break
        ch.append([-x for x in r])
    return [c for c in ch if uptrim(c) != [F(0)]]


def sign_changes(chain, t):
    sg = []
    for c in chain:
        v = upeval(c, F(t))
        if v != 0:
            sg.append(1 if v > 0 else -1)
    return sum(1 for i in range(len(sg) - 1) if sg[i] != sg[i + 1])


def sturm_count(a, lo, hi):
    """Number of distinct real roots of squarefree a in (lo, hi]."""
    ch = sturm_chain(a)
    return sign_changes(ch, lo) - sign_changes(ch, hi)


def upgcd(a, b):
    a, b = uptrim(a), uptrim(b)
    while b != [F(0)]:
        a, b = b, uptrim(uprem(a, b))
    return a


def is_squarefree(a):
    """True iff gcd(a, a') is a nonzero constant (so Sturm counting is valid)."""
    g = upgcd(a, upderiv(a))
    return len(uptrim(g)) == 1 and uptrim(g) != [F(0)]


# ======================================================================
# the hyperoctahedral group B_4 acting on R^4 by signed permutations
# ======================================================================

def b4_elements():
    els = []
    for perm in permutations(range(4)):
        for sg in product((1, -1), repeat=4):
            els.append((perm, sg))
    return els


def b4_apply(el, v):
    perm, sg = el
    return tuple(sg[i] * v[perm[i]] for i in range(4))


def b4_orbit(v, els):
    return sorted({b4_apply(el, v) for el in els})


def corner_vector_int(a_num, a_den, s):
    """Unnormalized integer corner vector (a, 1,...,1, 0,...,0) scaled by a_den."""
    v = [a_num] + [a_den] * s + [0] * (3 - s)
    return tuple(v)


def render_term(e):
    """(i,j,k) -> the paper's printed monomial name, e.g. (2,1,0) -> 'e_2^2e_3'."""
    out = ""
    for base, ex in zip(("e_2", "e_3", "e_4"), e):
        if ex == 1:
            out += base
        elif ex > 1:
            out += base + "^" + str(ex)
    return out


# ======================================================================
# CHECK GROUP 1 -- the exhibited object: decode it, count it, print it back
# ======================================================================

def checks_object(P_y):
    # NOTE: "len(set(PAPER_COEF)) == 11" used to appear here; dict keys are
    # already distinct, so that conjunct could not fail.  What can fail, and is
    # what was meant, is that the 11 rendered monomial names are distinct and
    # that the paper's printed list has no duplicate entry.
    ok = (len(PAPER_COEF) == 11
          and len({render_term(e) for e in PAPER_COEF}) == 11
          and len(set(PAPER_TERMS)) == 11
          and all(c != 0 for c in PAPER_COEF.values())
          and all(len(e) == 3 and all(x >= 0 for x in e) for e in PAPER_COEF))
    check("object_tuple_wellformed", ok,
          "%d distinct exponent triples, all coefficients nonzero"
          % len(PAPER_COEF))

    decoded = sorted((render_term(e), c) for e, c in PAPER_COEF.items())
    printed = sorted(PAPER_TERMS)
    check("object_matches_printed_eq2", decoded == printed,
          "11 terms of eq. (2) reproduced from the coefficient tuple")

    order = [nm for nm, _ in PAPER_TERMS]
    items = sorted(PAPER_COEF.items(),
                   key=lambda t: (order.index(render_term(t[0]))
                                  if render_term(t[0]) in order else 99))
    print("     p = " + " ".join(
        ("+ " if c > 0 else "- ") + str(abs(c)) + render_term(e)
        for e, c in items).lstrip("+ "))
    print("     expanded: %d monomials in y_1..y_4" % len(P_y))

    degs = [sum(e) for e in P_y]
    dy = max(degs)
    top = [e for e in P_y if sum(e) == dy]
    # NOTE: "len(top) > 0" used to be a conjunct here; a maximum is always
    # attained, so it could not fail for any input and has been removed.
    check("degree_in_x_is_14",
          2 * dy == PAPER_DEGREE_IN_X,
          "max y-degree %d -> x-degree %d (paper: %d); %d top monomials, "
          "e.g. y^%s with coefficient %d; no y_i occurs above power %d"
          % (dy, 2 * dy, PAPER_DEGREE_IN_X, len(top), min(top), P_y[min(top)],
             max(max(e) for e in P_y)))

    sym = True
    for perm in permutations(range(4)):
        for e, c in P_y.items():
            if P_y.get(tuple(e[perm[i]] for i in range(4)), 0) != c:
                sym = False
                break
        if not sym:
            break
    check("y_expansion_symmetric_hence_B4_invariant", sym,
          "invariant under all 24 permutations of y = (x_1^2,...,x_4^2)")

    v0 = mpeval(P_y, (F(1), F(0), F(0), F(0)))
    check("p_at_corner_s0_equals_2", v0 == PAPER_P_AT_S0,
          "p(v_{a,0}) = p(1,0,0,0) = %s (paper: %d)" % (v0, PAPER_P_AT_S0))


# ======================================================================
# CHECK GROUP 2 -- the moment functional used in Lemma 3, validated
# ======================================================================

def alphas_up_to(deg, nvars=4):
    for e in product(range(deg + 1), repeat=nvars):
        if sum(e) <= deg:
            yield e


def sign_group_sum(alpha):
    """sum over the 16 sign flips eps in {+-1}^4 of prod_i eps_i^{alpha_i}.

    This is the B_4 sign-subgroup symmetrization of the monomial x^alpha,
    computed by explicit summation (no special case for odd exponents).  It is
    0 exactly when some alpha_i is odd and 16 when all are even.  Since sigma
    is invariant under the sign subgroup, <x^alpha> equals the symmetrized
    mean, so a zero symmetrization *forces* <x^alpha> = 0.
    """
    tot = 0
    for eps in product((1, -1), repeat=4):
        t = 1
        for e, a in zip(eps, alpha):
            t *= e ** a
        tot += t
    return tot


def checks_moments():
    e1 = {(1, 0, 0, 0): 1, (0, 1, 0, 0): 1, (0, 0, 1, 0): 1, (0, 0, 0, 1): 1}
    vals = []
    for n in range(0, 8):
        vals.append(mean_of_y_poly(mppow(e1, n)))
    check("moment_normalization_and_e1_powers_are_one",
          len(vals) == 8 and all(v == 1 for v in vals),
          "<(y_1+...+y_4)^n> = 1 for n = 0..7 (e_1 = 1 on S^3); the n = 7 "
          "case sums 120 monomial moments with multinomial weights; got %s"
          % sorted(set(str(v) for v in vals)))

    # An earlier version of this check only asked whether mom_x returned 0 on
    # odd-degree monomials.  mom_x short-circuits on odd exponents, so that
    # condition could not fail for ANY paper and ANY object: it was vacuous.
    # It is now anchored to the computed sign-subgroup symmetrization, which is
    # the actual reason the moments vanish, together with an even-degree
    # control where the symmetrization is nonzero.
    odd = [a for a in alphas_up_to(15) if sum(a) % 2 == 1]
    bad_sym = [a for a in odd if sign_group_sum(a) != 0]
    ctrl = [(14, 0, 0, 0), (2, 4, 4, 4), (0, 0, 0, 2), (0, 0, 0, 0)]
    ctrl_ok = all(sign_group_sum(a) == 16 for a in ctrl)
    check("odd_moments_vanish_by_computed_sign_symmetrization",
          len(odd) > 2000 and not bad_sym
          and all(mom_x(a) == 0 for a in odd) and ctrl_ok,
          "the sign-subgroup symmetrization of all %d monomials of odd degree "
          "<= 15 is computed to be 0, hence <x^a> = 0, and mom_x agrees on "
          "every one; control: the same summation returns 16 (not 0) on the "
          "%d even monomials %s, so the gate is not vacuous"
          % (len(odd), len(ctrl), ctrl))

    pts = cell24_points()
    check("cell24_is_on_the_sphere",
          len(pts) == 24 and len(set(pts)) == 24
          and all(sum(c * c for c in v) == 1 for v in pts),
          "24 distinct unit vectors (the 24-cell)")

    agree = [a for a in alphas_up_to(5)
             if cell24_average(a, pts) == mom_x(a)]
    tot = list(alphas_up_to(5))
    check("moment_agrees_with_24cell_5design",
          len(agree) == len(tot),
          "all %d monomials of x-degree <= 5 integrate identically over the "
          "24-cell and under the moment functional" % len(tot))

    d6 = [a for a in product(range(7), repeat=4) if sum(a) == 6
          and cell24_average(a, pts) != mom_x(a)]
    check("cell24_control_fails_at_degree_6", len(d6) > 0,
          "%d degree-6 monomials disagree (e.g. %s: 24-cell %s vs sphere %s) "
          "-- the previous check is not vacuous"
          % (len(d6), d6[0], cell24_average(d6[0], pts), mom_x(d6[0])))

    bad = [a for a in alphas_up_to(7)
           if mom_y_pochhammer(a) != mom_y_doublefact(a)]
    check("two_closed_forms_for_the_moment_agree", not bad,
          "Pochhammer and double-factorial forms agree on all %d exponent "
          "vectors with |a| <= 7.  HONEST LABEL: these two forms are "
          "algebraically the SAME formula rewritten -- (1/2)_a = (2a-1)!!/2^a "
          "and (2)_m = (m+1)! turn one into the other -- so this is a "
          "transcription check, not an independent derivation.  Genuine "
          "independence is supplied by the next check and by the Laplacian "
          "route to <p>" % len(list(alphas_up_to(7))))

    # A genuinely independent validation of the moment functional over the FULL
    # degree range that p occupies (x-degree 14).  The 24-cell check above only
    # reaches degree 5 and the double-factorial form is an algebraic rewriting;
    # the iterated-Laplacian formula shares no code and no closed form with
    # mom_y_pochhammer.
    evens = [a for a in alphas_up_to(7)]
    disagree = [a for a in evens
                if mom_x(tuple(2 * ai for ai in a))
                != mean_homogeneous_laplacian(
                    {tuple(2 * ai for ai in a): 1}, 2 * sum(a))]
    check("moment_matches_laplacian_route_through_degree_14",
          not disagree and len(evens) == 330,
          "Dirichlet/Pochhammer moment equals Delta^m x^a / (2^m m! "
          "4*6*...*(2m+2)) on all %d even monomials of x-degree <= 14 -- the "
          "whole range used by p, by an algorithmically unrelated route"
          % len(evens))


# ======================================================================
# CHECK GROUP 3 -- Lemma 3: the spherical mean of p, two independent ways
# ======================================================================

def to_x_poly(P_y):
    """y_i = x_i^2 : double every exponent."""
    return {tuple(2 * a for a in e): c for e, c in P_y.items()}


def laplacian_x(Q):
    """Euclidean Laplacian of a polynomial in x_1..x_4 (dict form)."""
    out = {}
    for e, c in Q.items():
        for i in range(4):
            if e[i] >= 2:
                f = list(e)
                f[i] -= 2
                k = tuple(f)
                out[k] = out.get(k, 0) + c * e[i] * (e[i] - 1)
    return {k: v for k, v in out.items() if v != 0}


def mean_homogeneous_laplacian(Q, n):
    """<Q> on S^3 for Q homogeneous of even degree n, via Delta^{n/2} Q.

    <Q> = Delta^{m} Q / ( 2^m m! * 4*6*...*(2m+2) ),  m = n/2, d = 4.
    Independent of the Dirichlet-moment route.
    """
    m = n // 2
    R = dict(Q)
    for _ in range(m):
        R = laplacian_x(R)
    const = R.get((0, 0, 0, 0), 0)
    if set(R) - {(0, 0, 0, 0)}:
        raise AssertionError("Delta^m of a degree-n form is not constant")
    den = F(2) ** m
    for i in range(1, m + 1):
        den *= i
    for k in range(m):
        den *= (4 + 2 * k)
    return F(const) / den


def mean_via_laplacian(P_y):
    Q = to_x_poly(P_y)
    tot = F(0)
    for n in sorted({sum(e) for e in Q}):
        part = {e: c for e, c in Q.items() if sum(e) == n}
        tot += mean_homogeneous_laplacian(part, n)
    return tot


def checks_mean(P_y):
    m1 = mean_of_y_poly(P_y)
    check("mean_of_p_recomputed_matches_paper", m1 == PAPER_MEAN,
          "<p> = %s ; paper Lemma 3: %s" % (m1, PAPER_MEAN))
    m2 = mean_via_laplacian(P_y)
    check("mean_of_p_second_route_laplacian_agrees", m2 == m1,
          "harmonic/Laplacian route gives %s" % m2)
    check("mean_of_p_is_strictly_negative", m1 < 0,
          "%s < 0 (numerator %d)" % (m1, m1.numerator))


# ======================================================================
# CHECK GROUP 4 -- Lemma 4, step 1: reduction of p on the three profiles
# ======================================================================

def checks_reduction():
    ok = True
    for s in (1, 2, 3):
        tot = [F(0)]
        for v in profile_polys(s):
            tot = upadd(tot, v)
        if uptrim(tot) != [F(1)]:
            ok = False
    check("profile_squared_coordinates_sum_to_one", ok,
          "u + s*(1-u)/s = 1 for s = 1,2,3, so each profile lies on S^3.  "
          "SCOPE: this is a self-test of profile_polys, not a test of the "
          "paper -- it cannot fail for any object, only for a corrupted "
          "profile.  That the profile really is the squared-coordinate vector "
          "of v_{a,s} is tested separately, by agreeing with the independently "
          "built squared_coords() at 609 rational a^2")

    Hs = {}
    for s in (1, 2, 3):
        Pu = p_on_profile(PAPER_COEF, s)
        L, H = integer_scaling(Pu)
        Hint = [int(x) for x in H]
        Hs[s] = (L, [F(x) for x in H])
        good = (L == PAPER_SCALE[s] and Hint == PAPER_H[s])
        check("reduction_s%d_matches_paper_H%d" % (s, s), good,
              "least integer scaling %d (paper %d), deg %d, coeffs %s"
              % (L, PAPER_SCALE[s], len(Hint) - 1,
                 "match" if Hint == PAPER_H[s] else "%s vs %s"
                 % (Hint, PAPER_H[s])))
    return Hs


# ======================================================================
# CHECK GROUP 5 -- Lemma 4, step 2: Bernstein certificate + Sturm control
# ======================================================================

def checks_subdivisions():
    ok = True
    for s in (1, 2, 3):
        r = PAPER_SUBDIV[s]
        if not (len(r) == 4 and r[0] == 0 and r[3] == 1
                and all(r[i] < r[i + 1] for i in range(3))):
            ok = False
    check("subdivisions_are_partitions_of_unit_interval", ok,
          "0 = r_0 < r_1 < r_2 < r_3 = 1 for s = 1,2,3: "
          + "; ".join(str([str(x) for x in PAPER_SUBDIV[s]]) for s in (1, 2, 3)))


def checks_bernstein_selftest(Hs):
    ends = True
    for s in (1, 2, 3):
        H = Hs[s][1]
        for j in range(3):
            r, q = PAPER_SUBDIV[s][j], PAPER_SUBDIV[s][j + 1]
            b = bernstein_coeffs(H, r, q)
            if b[0] != upeval(H, r) or b[-1] != upeval(H, q):
                ends = False
    neg = bernstein_coeffs([F(-1, 2), F(1)], F(0), F(1))   # u - 1/2
    check("bernstein_implementation_selftest",
          ends and min(neg) == F(-1, 2) and max(neg) == F(1, 2),
          "endpoint coefficients reproduce H(r), H(q) on all 9 intervals; "
          "control: min beta of u-1/2 on [0,1] is %s < 0" % min(neg))

    # The endpoint test above pins only beta_0 and beta_d.  The lower bound
    # H >= min_k beta_k -- the sole source of the load-bearing constant L -- is
    # sound only if the whole beta vector really is the Bernstein coefficient
    # vector of H on I.  Verify that as an EXACT polynomial identity, and
    # verify the affine shift itself by 9-point agreement (H has degree <= 7,
    # so 9 points determine it).
    ident = True
    shift_ok = True
    for s in (1, 2, 3):
        H = Hs[s][1]
        for j in range(3):
            r, q = PAPER_SUBDIV[s][j], PAPER_SUBDIV[s][j + 1]
            shifted = affine_shift(H, r, q)
            for i in range(9):
                z = F(i, 8)
                if upeval(shifted, z) != upeval(H, r + (q - r) * z):
                    shift_ok = False
            back = bernstein_reconstruct(bernstein_coeffs(H, r, q))
            if [x for x in back] != [x for x in shifted]:
                ident = False
    bad_beta = bernstein_coeffs(Hs[1][1], F(0), F(1, 6))
    bad_beta[1] += 1
    ctrl = (bernstein_reconstruct(bad_beta)
            != affine_shift(Hs[1][1], F(0), F(1, 6)))
    # The identity alone does not yield the bound: that step also needs the
    # basis to be nonnegative on [0,1] and to sum to 1, and it needs the
    # residual H - min_k beta_k to be a NONNEGATIVE combination of the basis.
    # Both were formerly asserted in prose; compute them and report them here,
    # so that a failure of either fails this check.
    basis_ok, basis_txt = bernstein_basis_facts(Hs)
    resid_ok, resid_txt = bernstein_residual_facts(Hs)
    check("bernstein_coefficients_verified_as_exact_identity",
          ident and shift_ok and ctrl and basis_ok and resid_ok,
          "for all 9 intervals, sum_k beta_k C(d,k) z^k (1-z)^{d-k} equals "
          "H(r+(q-r)z) coefficient by coefficient, and the affine shift agrees "
          "with H at 9 points (> deg H), so eq. (3) is verified rather than "
          "assumed; control: perturbing one beta_k by 1 breaks the identity "
          "(%s).  %s.  %s.  The only step of the bound left outside the exact "
          "arithmetic is then the sign fact z >= 0, 1-z >= 0 for z in [0,1], "
          "which is the definition of the interval"
          % (ctrl, basis_txt, resid_txt))


def basis_poly(d, k):
    """comb(d,k) z^k (1-z)^{d-k}, expanded as coefficients in z."""
    zk = [F(0)] * k + [F(1)]
    return upscale(upmul(zk, uppow([F(1), F(-1)], d - k)), F(comb(d, k)))


def bernstein_basis_facts(Hs):
    """The two facts that turn min_k beta_k into a lower bound, COMPUTED.

    Both were previously left as assertions in prose: that the Bernstein basis
    is nonnegative on [0,1] and that it sums to 1.  Nonnegativity is reduced to
    an exact coefficient-by-coefficient identification of each basis element
    with comb(d,k) z^k (1-z)^{d-k}, comb(d,k) > 0 -- after which the only step
    outside the arithmetic is that z and 1-z are >= 0 on [0,1], which is what
    the interval means.

    Returns (ok, text); it is reported inside the eq. (3) identity check,
    which is the check whose bound it licenses.
    """
    degs = sorted({len(Hs[s][1]) - 1 for s in (1, 2, 3)})
    pos = True
    prod_ok = True
    unity_ok = True
    for d in degs:
        for k in range(d + 1):
            if comb(d, k) < 1:
                pos = False
            e_k = [F(0)] * (d + 1)
            e_k[k] = F(1)
            if bernstein_reconstruct(e_k) != basis_poly(d, k):
                prod_ok = False
        if bernstein_reconstruct([F(1)] * (d + 1)) != [F(1)] + [F(0)] * d:
            unity_ok = False
    # Control: the partition-of-unity test must not pass for a weight vector
    # that is not identically 1.
    dctrl = degs[0]
    wctrl = [F(1)] * (dctrl + 1)
    wctrl[1] = F(2)
    ctrl = bernstein_reconstruct(wctrl) != [F(1)] + [F(0)] * dctrl
    return (pos and prod_ok and unity_ok and ctrl,
            "BASIS, for d in %s: every basis element used by the "
            "reconstruction is identified coefficient by coefficient with "
            "comb(d,k) z^k (1-z)^{d-k} and every comb(d,k) >= 1, so each is "
            ">= 0 on [0,1]; and sum_k comb(d,k) z^k (1-z)^{d-k} is the "
            "constant polynomial 1 exactly -- nonnegativity and partition of "
            "unity are computed here, not asserted; control: changing one "
            "weight from 1 to 2 destroys the partition of unity (%s)"
            % (degs, ctrl))


def bernstein_residual_facts(Hs):
    """H(r+(q-r)z) - min_k beta_k as an exact NONNEGATIVE combination.

    Returns (ok, text); reported inside the eq. (3) identity check.
    """
    ident = True
    coef_nonneg = True
    ctrl = True
    worst_pair = None
    for s in (1, 2, 3):
        H = Hs[s][1]
        for j in range(3):
            r, q = PAPER_SUBDIV[s][j], PAPER_SUBDIV[s][j + 1]
            beta = bernstein_coeffs(H, r, q)
            m = min(beta)
            resid = [b - m for b in beta]
            if any(x < 0 for x in resid):
                coef_nonneg = False
            shifted = affine_shift(H, r, q)
            rhs = list(shifted)
            rhs[0] = rhs[0] - m
            if bernstein_reconstruct(resid) != rhs:
                ident = False
            # Control: the identity must NOT survive a perturbation of the
            # residual vector, i.e. it really pins the representation.
            bad = list(resid)
            bad[1] = bad[1] + 1
            if bernstein_reconstruct(bad) == rhs:
                ctrl = False
            if worst_pair is None or min(resid) < worst_pair:
                worst_pair = min(resid)
    return (ident and coef_nonneg and ctrl,
            "RESIDUAL, for all 9 intervals: H(r+(q-r)z) - min_k beta_k equals "
            "sum_k (beta_k - min beta) comb(d,k) z^k (1-z)^{d-k} coefficient "
            "by coefficient, and all 9 residual coefficient vectors are "
            ">= 0 (smallest entry %s), so with the basis facts above the step "
            "min_k beta_k -> lower bound on H is a computation, not an "
            "assertion; control: perturbing one residual coefficient by 1 "
            "breaks the identity on every interval (%s)" % (worst_pair, ctrl))


def checks_bernstein_table(Hs):
    """Recompute the nine minimum Bernstein coefficients of the table."""
    mins = {}
    for s in (1, 2, 3):
        H = Hs[s][1]
        got = []
        for j in range(3):
            r, q = PAPER_SUBDIV[s][j], PAPER_SUBDIV[s][j + 1]
            got.append(min(bernstein_coeffs(H, r, q)))
        mins[s] = got
        check("bernstein_minima_s%d_match_table" % s, got == PAPER_BMIN[s],
              "computed %s ; paper %s"
              % ([str(x) for x in got], [str(x) for x in PAPER_BMIN[s]]))
    allmins = [x for s in (1, 2, 3) for x in mins[s]]
    check("all_nine_bernstein_minima_are_positive",
          len(allmins) == 9 and all(x > 0 for x in allmins),
          "smallest of the nine is %s > 0" % min(allmins))

    # The table entries are used as LOWER BOUNDS; test that property directly.
    nsamp = 0
    lb_ok = True
    slack = None
    for s in (1, 2, 3):
        H = Hs[s][1]
        for j in range(3):
            r, q = PAPER_SUBDIV[s][j], PAPER_SUBDIV[s][j + 1]
            for i in range(21):
                t = r + (q - r) * F(i, 20)
                d = upeval(H, t) - mins[s][j]
                nsamp += 1
                if d < 0:
                    lb_ok = False
                if slack is None or d < slack:
                    slack = d
    check("table_entries_really_are_lower_bounds_for_H", lb_ok,
          "min_k beta_{I,k} <= H_s(t) at %d sample points spread over the "
          "nine intervals; smallest slack %s >= 0 (the relation is <=, not <: "
          "equality occurs at an interval endpoint where beta_0 = H(r) is the "
          "minimum).  This is a spot-check of a bound established exactly, and "
          "for every u, by bernstein_coefficients_verified_as_exact_identity "
          "above: the eq. (3) identity, the computed nonnegativity and "
          "partition of unity of the basis, and the exact nonnegative "
          "representation of H(r+(q-r)z) - min_k beta_k" % (nsamp, slack))
    return mins


def checks_sturm(Hs):
    for s in (1, 2, 3):
        H = Hs[s][1]
        sf = is_squarefree(H)
        nroots = sturm_count(H, F(0), F(1)) if sf else -1
        h0 = upeval(H, F(0))
        h1 = upeval(H, F(1))
        check("sturm_H%d_has_no_root_in_unit_interval" % s,
              sf and nroots == 0 and h0 > 0 and h1 > 0,
              "squarefree=%s, roots in (0,1] = %s, H_%d(0) = %s, H_%d(1) = %s"
              % (sf, nroots, s, h0, s, h1))
    c1 = sturm_count([F(1, 8), F(-3, 4), F(1)], F(0), F(1))   # (u-1/4)(u-1/2)
    c2 = sturm_count([F(1), F(0), F(1)], F(0), F(1))          # u^2+1
    c3 = sturm_count([F(-1), F(2)], F(0), F(1))               # 2u-1
    c4 = sturm_count(Hs[1][1], F(-10), F(10))
    check("sturm_implementation_control",
          (c1, c2, c3) == (2, 0, 1) and c4 >= 1,
          "root counts in (0,1]: (u-1/4)(u-1/2) -> %d, u^2+1 -> %d, "
          "2u-1 -> %d; H_1 has %d real roots in (-10,10] so the "
          "no-root verdict above is not vacuous" % (c1, c2, c3, c4))


# ======================================================================
# CHECK GROUP 6 -- p at genuine corner-vector points (exact, rational a^2)
# ======================================================================

def squared_coords(a2, s):
    """Squared coordinates of v_{a,s} for rational a^2 = a2 > 0."""
    D = a2 + s
    return tuple([a2 / D] + [F(1, 1) / D] * s + [F(0)] * (3 - s))


def eval_p_y(P_y, y):
    dmax = max(max(e) for e in P_y)
    pw = []
    for yi in y:
        col = [F(1)]
        for _ in range(dmax):
            col.append(col[-1] * yi)
        pw.append(col)
    tot = F(0)
    for e, c in P_y.items():
        t = F(c)
        for i in range(4):
            if e[i]:
                t *= pw[i][e[i]]
        tot += t
    return tot


def a2_grid():
    g = [F(j, 64) for j in range(1, 129)]
    g += [F(j) for j in range(1, 41)]
    g += [F(2) ** k for k in range(6, 25)]
    g += [F(1, 2) ** k for k in range(7, 25)]
    return sorted(set(g))


def checks_orbit_grid(P_y, Hs, mins):
    grid = a2_grid()
    # COMPUTED, not the literal 2 taken from the paper: v_{a,0} = (1,0,0,0) for
    # every a > 0, so the s = 0 entry of the observed minima is p(1,0,0,0).
    observed = {0: eval_p_y(P_y, (F(1), F(0), F(0), F(0)))}
    for s in (1, 2, 3):
        L, H = Hs[s]
        bound = min(mins[s]) / L
        worst = None
        agree = True
        for a2 in grid:
            y = squared_coords(a2, s)
            val = eval_p_y(P_y, y)
            u = a2 / (a2 + s)
            if val != upeval(H, u) / L:
                agree = False
            if worst is None or val < worst:
                worst = val
        observed[s] = worst
        check("p_positive_on_%d_true_orbit_points_s%d" % (len(grid), s),
              agree and worst > 0 and worst >= bound,
              "min over %d rational a^2 in [2^-24, 2^24] is %s > 0; "
              "matches H_%d(u)/%d at every point: %s; certified lower "
              "bound %s" % (len(grid), worst, s, L, agree, bound))
    y0 = (F(1), F(0), F(0), F(0))
    v0 = eval_p_y(P_y, y0)
    check("p_positive_on_orbit_s0",
          v0 > 0 and v0 == PAPER_P_AT_S0,
          "p(v_{a,0}) = p(1,0,0,0) = %s > 0 for every a > 0 (paper: %d)"
          % (v0, PAPER_P_AT_S0))

    H1 = Hs[1][1]
    pal = all(upeval(H1, F(k, 12)) == upeval(H1, 1 - F(k, 12))
              for k in range(13))
    check("H1_symmetric_under_u_to_1_minus_u", pal,
          "H_1(u) = H_1(1-u), as forced by symmetry of p and the 2-point "
          "profile (u,1-u,0,0); this is why the outer two table entries "
          "coincide")
    return observed


def u_targets():
    """Rational targets spread over the open interval (0,1)."""
    ts = [F(k, 64) for k in range(1, 64)]
    ts += [F(1, 1000), F(999, 1000)]
    return sorted(set(ts))


def u_sweep_facts(P_y, Hs, mins):
    """The a-quantifier on the u side: u = a^2/(a^2+s) lies in (0,1) and
    sweeps it.  Returns (ok, text), reported inside the bound check.

    The certificate is stated for u in [0,1]; what makes it apply to every
    generalized corner vector is that u(a) never leaves (0,1) and that every
    u in (0,1) is realised by some a > 0.  Both were previously left implicit.
    Surjectivity is turned into a test with force by inverting the map at each
    target -- a^2 = s u/(1-u) -- and re-evaluating p at that genuine orbit
    point against H_s(u)/scale_s and against the interval's Bernstein minimum.
    """
    targets = u_targets()
    inside = True
    inverse_ok = True
    covered = True
    reduction_ok = True
    bound_ok = True
    npts = 0
    for s in (1, 2, 3):
        scale, H = Hs[s]
        rr = PAPER_SUBDIV[s]
        for u in targets:
            a2 = F(s) * u / (1 - u)
            if a2 <= 0 or a2 / (a2 + s) != u:
                inverse_ok = False
            if not (0 < u < 1):
                inside = False
            idx = [j for j in range(3) if rr[j] <= u <= rr[j + 1]]
            if not idx:
                covered = False
                continue
            val = eval_p_y(P_y, squared_coords(a2, s))
            npts += 1
            if val != upeval(H, u) / scale:
                reduction_ok = False
            if val < mins[s][idx[0]] / scale:
                bound_ok = False

    mono = True
    grid_cov = True
    grid = a2_grid()
    for s in (1, 2, 3):
        rr = PAPER_SUBDIV[s]
        us = [a2 / (a2 + s) for a2 in grid]
        if any(us[i] >= us[i + 1] for i in range(len(us) - 1)):
            mono = False
        for u in us:
            if not (0 < u < 1) or not any(rr[j] <= u <= rr[j + 1]
                                          for j in range(3)):
                grid_cov = False
    return (inside and inverse_ok and covered and reduction_ok and bound_ok
            and mono and grid_cov and npts == 3 * len(targets),
            "u-RANGE: u = a^2/(a^2+s) is strictly increasing in a^2 across %d "
            "sampled a^2 in [2^-24, 2^24] and lies strictly inside (0,1) and "
            "inside one of the paper's three subintervals at every one of "
            "them; and it is onto, so nothing is missed and nothing beyond "
            "[0,1] is needed: for each of %d rational targets u in (0,1) "
            "(down to 1/1000 and up to 999/1000) the preimage "
            "a^2 = s u/(1-u) reproduces u exactly, and p at that genuine "
            "orbit point equals H_s(u)/scale_s and is at least that "
            "subinterval's Bernstein minimum -- %d further exact evaluations "
            "of p" % (len(grid), len(targets), npts))


# ======================================================================
# CHECK GROUP 7 -- B_4 orbits: sizes, p constant, antipodality, degree 15
# ======================================================================

REPS = [(1, 1, 0), (3, 1, 0), (3, 1, 1), (1, 1, 1), (1, 2, 1),
        (3, 1, 2), (1, 1, 2), (3, 1, 3), (1, 1, 3), (1, 2, 3)]


def checks_orbits(P_y, Hs):
    els = b4_elements()
    check("b4_has_384_elements", len(els) == 384 and len(set(els)) == 384,
          "4! * 2^4 signed permutations of R^4")

    sizes = {}
    const_ok = True
    anti_ok = True
    sum_ok = True
    for (an, ad, s) in REPS:
        v = corner_vector_int(an, ad, s)
        orb = b4_orbit(v, els)
        sizes[(an, ad, s)] = len(orb)
        n2 = sum(c * c for c in v)
        vals = set()
        tot = F(0)
        for w in orb:
            y = tuple(F(c * c, n2) for c in w)
            val = eval_p_y(P_y, y)
            vals.add(val)
            tot += val
        if len(vals) != 1:
            const_ok = False
        if tuple(-c for c in v) not in set(orb):
            anti_ok = False
        # expected value from the INDEPENDENT reduced-polynomial route
        a2 = F(an * an, ad * ad)
        if s == 0:
            exp = F(PAPER_P_AT_S0)
        else:
            L, H = Hs[s]
            exp = upeval(H, a2 / (a2 + s)) / L
        if tot != len(orb) * exp:
            sum_ok = False
    check("orbit_sizes_enumerated",
          sorted(sizes.values()) == sorted(
              [8, 8, 48, 24, 48, 96, 32, 64, 16, 64]),
          "; ".join("a=%s,s=%d:%d" % (F(an, ad), s, sizes[(an, ad, s)])
                    for (an, ad, s) in REPS))
    check("p_is_constant_on_every_orbit", const_ok,
          "one value of p per orbit over all %d representatives, so "
          "sum_{x in orbit} p(x) = |orbit| p(v)" % len(REPS))
    check("orbit_sum_equals_size_times_reduced_polynomial", sum_ok,
          "sum_{x in orbit} p(x) = |orbit| * H_s(u)/%s, computed by explicit "
          "summation over the group orbit on one side and by the reduced "
          "univariate polynomial on the other"
          % str(tuple(PAPER_SCALE[s] for s in (1, 2, 3))))
    check("orbits_are_antipodal", anti_ok,
          "-v lies in the B_4 orbit of v for all %d representatives "
          "(since -I in B_4).  SCOPE: a self-test of the group construction -- "
          "-I is in B_4 for every v, so this cannot fail for any object, only "
          "for a corrupted b4_elements()" % len(REPS))
    return els


def monomials_of_degree(d):
    out = []
    for i in range(d + 1):
        for j in range(d - i + 1):
            for k in range(d - i - j + 1):
                out.append((i, j, k, d - i - j - k))
    return out


def orbit_monomial_sum(orb, alpha, dmax):
    pw = []
    for w in orb:
        col = []
        for c in w:
            row = [1]
            for _ in range(dmax):
                row.append(row[-1] * c)
            col.append(row)
        pw.append(col)
    tot = 0
    for col in pw:
        t = 1
        for i in range(4):
            t *= col[i][alpha[i]]
        tot += t
    return tot


def checks_degree15(els):
    m15 = monomials_of_degree(15)
    m14 = monomials_of_degree(14)
    zero_all = True
    for (an, ad, s) in REPS:
        orb = b4_orbit(corner_vector_int(an, ad, s), els)
        for a in m15:
            if orbit_monomial_sum(orb, a, 15) != 0:
                zero_all = False
                break
        if not zero_all:
            break
    check("degree15_orbit_sums_vanish", len(m15) == 816 and zero_all,
          "all %d monomials of degree 15 have zero sum over every "
          "corner-vector orbit" % len(m15))

    # As with the odd-moment gate, "mom_x returns 0 on odd degree" could not
    # fail for any input.  Anchor it to the computed sign-subgroup
    # symmetrization (the reason sigma kills these monomials) and add a
    # degree-14 control where that symmetrization is nonzero.
    bad_sym15 = [a for a in m15 if sign_group_sum(a) != 0]
    ctrl14sym = [a for a in m14 if sign_group_sum(a) != 0]
    check("degree15_spherical_means_vanish",
          not bad_sym15 and all(mom_x(a) == 0 for a in m15)
          and len(ctrl14sym) > 0,
          "the sign-subgroup symmetrization of all %d degree-15 monomials is "
          "computed to be 0, so each has zero spherical mean (mom_x agrees), "
          "hence exactness through degree 14 is equivalent to exactness "
          "through degree 15 (Remark 5); control: %d of the %d degree-14 "
          "monomials have NONZERO symmetrization, so the gate is not vacuous"
          % (len(m15), len(ctrl14sym), len(m14)))

    orb = b4_orbit(corner_vector_int(3, 1, 3), els)
    nz = [a for a in m14 if orbit_monomial_sum(orb, a, 14) != 0]
    check("degree14_orbit_sums_control",
          len(nz) > 0,
          "%d of %d degree-14 monomials have NONZERO orbit sum (e.g. %s -> "
          "%d), so the vanishing above is a property of odd degree, not of "
          "the summation" % (len(nz), len(m14), nz[0],
                             orbit_monomial_sum(orb, nz[0], 14)))


# ======================================================================
# CHECK GROUP 8 -- the conclusion: the cubature formula cannot be exact
# ======================================================================

def certified_lower_bound(Hs, mins, p_at_s0):
    """Strictly positive rational L with p(v_{a,s}) >= L for ALL a>0, s.

    Bernstein basis functions are nonnegative and sum to 1 -- both facts
    computed, not assumed, by bernstein_basis_facts, and combined with the
    residual representation of bernstein_residual_facts -- so on each
    subinterval H_s >= min_k beta_k; the subintervals cover [0,1] and
    u = a^2/(a^2+s) lies in (0,1) and sweeps it (u_sweep_facts);
    p(v_{a,s}) = H_s(u)/scale_s.

    p_at_s0 must be the COMPUTED value of p at (1,0,0,0).  An earlier version
    seeded this with the literal F(2) taken from the paper, which made the
    s = 0 branch of the load-bearing bound asserted rather than derived.
    """
    L = F(p_at_s0)                 # the s = 0 value, computed from the object
    for s in (1, 2, 3):
        L = min(L, min(mins[s]) / Hs[s][0])
    return L


def checks_conclusion(P_y, Hs, mins, els, observed):
    L = certified_lower_bound(Hs, mins,
                              eval_p_y(P_y, (F(1), F(0), F(0), F(0))))
    mean = mean_of_y_poly(P_y)
    # The "for every a > 0" of this check needs the u-side of the quantifier:
    # u = a^2/(a^2+s) never leaves (0,1), and every u in (0,1) is attained.
    sweep_ok, sweep_txt = u_sweep_facts(P_y, Hs, mins)
    check("certified_positive_lower_bound_for_p_on_all_orbits",
          L > 0 and sweep_ok,
          "p(v_{a,s}) >= %s > 0 for every a > 0 and every s in {0,1,2,3}.  %s"
          % (L, sweep_txt))
    check("certified_bound_is_below_every_evaluated_value",
          all(L <= observed[s] for s in (0, 1, 2, 3)),
          "L = %s does not exceed the smallest p value actually observed for "
          "any s (%s)" % (L, min(observed.values())))
    check("refutation_lower_bound_exceeds_the_spherical_mean",
          L > 0 > mean,
          "any positive normalized formula returns >= %s for p, while "
          "int p dsigma = %s < 0 : Theorem 1 holds, Problem 6.1 has a "
          "negative answer" % (L, mean))

    import random
    rnd = random.Random(20250824)
    size_cache = {}
    worst = None
    ok = True
    for _ in range(200):
        m = rnd.randint(1, 5)
        rhs = F(0)
        terms = []
        for _ in range(m):
            s = rnd.randint(0, 3)
            an, ad = rnd.randint(1, 24), rnd.randint(1, 8)
            v = corner_vector_int(an, ad, s)
            key = (s, F(an, ad))
            if key not in size_cache:
                size_cache[key] = len(b4_orbit(v, els))
            n2 = sum(c * c for c in v)
            y = tuple(F(c * c, n2) for c in v)
            terms.append((size_cache[key], eval_p_y(P_y, y),
                          rnd.randint(1, 50)))
        norm = sum(w * sz for sz, _, w in terms)
        for sz, val, w in terms:
            rhs += F(w, 1) / norm * sz * val
        if not (rhs >= L and rhs > 0 and rhs != mean):
            ok = False
        if worst is None or rhs < worst:
            worst = rhs
    check("200_random_positive_normalized_formulas_all_disagree", ok,
          "every sampled sum_i W_i |orbit_i| p(v_i) with W_i > 0 and "
          "sum_i W_i |orbit_i| = 1 lies above %s; smallest seen %s, versus "
          "the true value %s" % (L, worst, mean))


# ======================================================================
def main():
    print("verify.py -- Positive Generalized Corner-Vector Designs on S^3")
    print("     input from the paper: the 11-term coefficient tuple of eq. (2)")
    P_y = build_p_in_y(PAPER_COEF)
    checks_object(P_y)
    checks_moments()
    checks_mean(P_y)
    Hs = checks_reduction()
    checks_subdivisions()
    checks_bernstein_selftest(Hs)
    mins = checks_bernstein_table(Hs)
    checks_sturm(Hs)
    observed = checks_orbit_grid(P_y, Hs, mins)
    els = checks_orbits(P_y, Hs)
    checks_degree15(els)
    checks_conclusion(P_y, Hs, mins, els, observed)
    print("NOTE: the paper contains no exhaustive census.  Its universal "
          "quantifier (all a > 0, all s in {0,1,2,3}) is a continuum and is "
          "discharged exactly here by the Bernstein certificate over a "
          "partition of [0,1], independently re-certified by Sturm root "
          "counting; every computational assertion of the refutation itself "
          "(eq. (2), Lemma 3, Lemma 4 with its table, Remark 5 and Theorem 1) "
          "is recomputed above from the eleven coefficients of eq. (2) alone.")
    print("NOT RE-RUN: (i) the SHARPNESS half of the corollary -- the two "
          "positive-weight 13-designs of this type reported in the cited "
          "paper, whose nodes and weights are given there numerically to six "
          "significant figures -- is not verified here, as the paper itself "
          "states; (ii) claims about the cited literature (its definition of "
          "a weighted design, and its record that signed 15-designs exist on "
          "S^3) are read, not computed; (iii) two elementary sign facts are "
          "used as definitions rather than computed: z >= 0 and 1-z >= 0 for "
          "z in [0,1], which is what the Bernstein interval means, and "
          "a^2 > 0 implies 0 < a^2/(a^2+s) < 1; (iv) this folder ships no "
          "mutation log, no corrupted-run transcript and no second "
          "implementation -- the falsifiability evidence here is the controls "
          "printed beside the checks, which are re-run on every invocation.")
    return verdict()


if __name__ == "__main__":
    sys.exit(main())

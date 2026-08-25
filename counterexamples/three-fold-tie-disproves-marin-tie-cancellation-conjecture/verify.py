#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- independent verification of

    "A Counterexample to Marin's Conjecture 5.14"
    (beta = (12,11,10,3,2,1,0) at (t,r) = (3,2))

against Conjecture 5.14 of C. Marin, arXiv:2608.18302v1.

Python 3.9+, standard library only.  All arithmetic is exact: ints and
fractions.Fraction.  No floating point occurs anywhere in this file.

Run:
    python3 verify.py                      # the mathematics (part A)
    python3 verify.py --source MARIN.tex   # additionally check the LaTeX
                                              # source of arXiv:2608.18302v1,
                                              # which the reader fetches from
                                              # arXiv (part B)
    python3 verify.py --source MARIN.tex --paper paper.tex
                                              # both directions: the clauses
                                              # the paper quotes are in the
                                              # source, AND the source's
                                              # numbering makes the refuted
                                              # statement Conjecture 5.14

A GAP REGISTER is printed unconditionally at the end.  Read it: several checks
are labelled TAUTOLOGICAL because they cannot fail, and four clauses are
imported from Marin rather than derived.

================================================================================
TAKEN FROM THE PAPER  (data and imported conventions; not re-derived here)
================================================================================
 1. The exhibited shape  beta = (12,11,10,3,2,1,0)  and the parameters
    (t,r) = (3,2).  Also the calibration shape beta = (8,7,6,3,2,1,0), which is
    Marin's own smallest tie witness (his Remark 5.13).
 2. The beta-to-lambda dictionary quoted from Marin: N = t + 2r,
    delta = (N-1, N-2, ..., 0), lambda = beta - delta.
 3. Littlewood's stable orthogonal restriction formula (paper eq. (2)):
        a^B_Lambda = sum over partitions gamma with ALL PARTS EVEN of
                     c^{lambda}_{gamma,Lambda},
    valid because ell(lambda) <= N/2.  (The paper cites Koike-Terada, Adv. Math.
    79 (1990) 104-135 for it.  The formula is imported here; only its parts-even
    vs columns-even CONVENTION is discriminated -- see G4 of the gap register.)
 4. The folding data quoted from Marin's Proposition 5.10 and 5.19(ii):
      * V(Lambda) = 2(Lambda + rho_{B_3}) = (2L1+5, 2L2+3, 2L3+1);
      * at t = 3 there is exactly one nonzero folded residue class, namely
        {+-1 mod 3}, and a transversal is a single index;
      * S_min takes, in that class, the index with the SMALLEST value of V, and
        the top numerator index is V restricted to the complement, sorted
        decreasingly -- call it X(Lambda);
      * the top weight T(Lambda) of c(Lambda,.) satisfies
        2(T + rho_{D_2}) = X - 3*(1,1) with rho_{D_2} = (1,0), i.e.
        T = ((X1-5)/2, (X2-3)/2);
      * c_top(Lambda) = c(Lambda, T(Lambda)) lies in {+-1}   (Prop. 5.10(iv);
        this membership, and nothing more about the sign, is all that is used).
 5. The order in Marin's definition of M(beta): the Weyl majorisation order
        mu <= mu'  iff  mu in conv(W(C_r) mu').
 6. The statement of Conjecture 5.14 itself: if M(beta) is attained by more than
    one Lambda then  sum_{T(Lambda)=M(beta)} a^B_Lambda c_top(Lambda) = 0.
 7. The numbers the paper asserts, together with a few further comparison
    targets fixed below where the paper prints no value, all of which this
    program RECOMPUTES AND COMPARES and
    never seeds into a derivation:  20 admissible Lambda;  a^B_Lambda = 1 on
    exactly those;  T(Lambda) <= (6,6) always;  equality exactly at (6,6,6),
    (6,6,2), (6,6,0);  V(6,6,4) = (17,15,9) and T(6,6,4) = (5,3);
    conv(W(C_2)(6,6)) = [-6,6]^2;  every T in [-2,6]^2;  7056 ordered pairs of
    3-row partitions in the 6-box;  the calibration tie T = (2,2) attained by
    exactly (2,2,2) and (2,2,0);  max beta = 12 > 9 = the box of Marin's
    Observation 5.27;  Conjecture 5.14 is the 14th numbered environment of
    Marin's Section 5.

================================================================================
DERIVED HERE  (everything the checks actually decide)
================================================================================
 A. Schur polynomials of 3 variables, built from Gelfand-Tsetlin patterns, and
    validated independently against the GL_3 Weyl dimension formula, against
    S_3-symmetry and against homogeneity.
 B. GL_3 Littlewood-Richardson coefficients c^nu_{gamma,Lambda} by TWO mutually
    independent algorithms:
      (a) expand s_gamma * s_Lambda into monomials and strip leading dominant
          monomials (the polynomial method);
      (b) count Littlewood-Richardson skew tableaux of shape nu/gamma with
          content Lambda and lattice (ballot) reading word -- no polynomials at
          all.
    The two are compared on all 7056 ordered pairs for nu = (6,6,6) and on all
    100 pairs for nu = (2,2,2).
 C. The rectangular duality c^{(6,6,6)}_{gamma,Lambda} = [gamma = (6-L3,6-L2,
    6-L1)] -- CHECKED over all 7056 pairs, not quoted.
 D. a^B_Lambda from item 3 above, for every partition Lambda with at most 3
    parts, parts <= 12 and |Lambda| <= 18 (a strict superset of the 6-box) plus
    a sample of 4-row Lambda; hence the set of Lambda with a^B_Lambda != 0, its
    cardinality, and the value taken there.
 E. That the convention in item 3 is "parts even" and not "columns even", by
    deriving Sym^2(std) = o_(2) + trivial and Alt^2(std) = o_(1,1) under it and
    showing the columns-even variant gets both wrong.
 F. The mod-3 congruences: that V(Lambda) never has all three entries divisible
    by 3, both by exhaustion over the 20 and by the symbolic argument
    (L1 = 2 and L3 = 4 forced, contradicting L1 >= L3).
 G. V, X, T for all 20 admissible Lambda; the coordinatewise maximum; the
    maximizer set and its size; the (6,6,4) -> (5,3) drop.  Also (check 9b) the
    decisive step behind the paper's eq. (4): with supp nu indexed as the proof
    of Marin's clause fixes it (V restricted to the complement of a transversal,
    sorted decreasingly), X(Lambda) is the coordinatewise dominant numerator
    index, and every index above X(Lambda) - 3*(1,1) misses supp nu -- together
    with all four of its W(D_2) images, which is the reflection route back into
    the support that the paper's argument leaves implicit.
 H. The W(C_2)-orbit of (6,6); that conv of that orbit is exactly [-6,6]^2, via
    exact rational convex-combination certificates in both directions; an
    independent absolute-majorisation test of the same order; and that no other
    T dominates (6,6) in it, so the maximizer set is order-independent.
 I. The value of the tied sum for every one of the 2^3 sign patterns allowed by
    c_top in {+-1}: its parity and its nonvanishing.
 J. The calibration run at lambda = (2,2,2), reproducing Marin's two-fold tie.
 K. (part B, only with --source) the environment numbering of Marin's Section 5,
    the absence of counter resets, and the presence of every clause the paper
    quotes.
"""

import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product

# ----------------------------------------------------------------------------
# DATA TAKEN FROM THE PAPER
# ----------------------------------------------------------------------------
T_PARAM = 3                          # t
R_PARAM = 2                          # r
BETA = (12, 11, 10, 3, 2, 1, 0)      # the exhibited shape
TWO_RHO_B3 = (5, 3, 1)               # 2*rho_{B_3}, rho_{B_3} = (5/2, 3/2, 1/2)
RHO_D2 = (1, 0)                      # rho_{D_2}

CLAIM_N = 7                          # N = t + 2r
CLAIM_DELTA = (6, 5, 4, 3, 2, 1, 0)  # delta
CLAIM_LAMBDA = (6, 6, 6)             # lambda = beta - delta
CLAIM_NUM_ADMISSIBLE = 20            # # Lambda with a^B_Lambda != 0
CLAIM_AB_VALUE = 1                   # the common value of a^B_Lambda there
CLAIM_M = (6, 6)                     # M(beta)
CLAIM_MAXIMIZERS = ((6, 6, 6), (6, 6, 2), (6, 6, 0))
CLAIM_V_664 = (17, 15, 9)            # V(6,6,4)
CLAIM_X_664 = (15, 9)                # X(6,6,4)
CLAIM_T_664 = (5, 3)                 # T(6,6,4)
CLAIM_V_BOUND = (17, 15, 13)         # entries of V bounded by this
CLAIM_X_BOUND = (17, 15)             # X <= this coordinatewise
CLAIM_T_RANGE = (-2, 6)             # every T lies in [-2,6]^2
CLAIM_HULL_BOX = 6                   # conv(W(C_2)(6,6)) = [-6,6]^2
CLAIM_PAIR_COUNT = 7056              # ordered pairs of 3-row partitions in the box

CALIB_BETA = (8, 7, 6, 3, 2, 1, 0)   # Marin's Remark 5.13 smallest witness
CALIB_LAMBDA = (2, 2, 2)
CALIB_T = (2, 2)
CALIB_TIE = ((2, 2, 2), (2, 2, 0))
CALIB_NUM_ADMISSIBLE = 4

CLAIM_OBS_BOX_MAXBETA = 9            # Observation 5.27 box at (t,r)=(3,2)
CLAIM_CONJ_INDEX = 14                # Conjecture 5.14 = 14th env of Section 5
# The three values below are locators inside the LaTeX source of
# arXiv:2608.18302v1 -- properties of Marin's e-print, not of this folder.  They
# are used only by part B, and only against the file the reader supplies to
# --source after fetching that e-print (see the command printed in the gap
# register).  Part B prints what it finds at each of them, so a reader holding a
# different copy of the source can see immediately whether they still apply.
CLAIM_CONJ_LINE = 1654               # line carrying the conjecture's \begin
CLAIM_CONJ_LABEL = "conj:tie"        # the \label that environment carries there
CLAIM_SEC5_LINE = 1133               # Section 5 begins just after this line
CLAIM_SEC5_ENV_TOTAL = None          # derived, then only reported

# ----------------------------------------------------------------------------
# check bookkeeping
# ----------------------------------------------------------------------------
RESULTS = []


def say(msg=""):
    print(msg)


def record(name, ok, detail=""):
    """Register one check and print it in the fixed PASS/FAIL format."""
    RESULTS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    line = "%s %s" % (tag, name)
    if detail:
        line += "  --  " + detail
    print(line)
    return bool(ok)


# ----------------------------------------------------------------------------
# partitions
# ----------------------------------------------------------------------------
def partitions_in_box(ncols, nrows=3):
    """All partitions with at most nrows parts, every part in [0, ncols].

    Returned as tuples of length exactly nrows, weakly decreasing.
    """
    out = []
    for p in product(range(ncols, -1, -1), repeat=nrows):
        if all(p[i] >= p[i + 1] for i in range(nrows - 1)):
            out.append(p)
    return sorted(out, reverse=True)


def partitions_of(n, nparts, maxpart=None):
    """All weakly decreasing tuples of length nparts, nonneg, summing to n."""
    if maxpart is None:
        maxpart = n
    out = []

    def rec(rem, slots, cap, acc):
        if slots == 0:
            if rem == 0:
                out.append(tuple(acc))
            return
        top = min(cap, rem)
        for v in range(top, -1, -1):
            if v * slots < rem:          # cannot reach rem with slots parts <= v
                break
            rec(rem - v, slots - 1, v, acc + [v])

    rec(n, nparts, maxpart, [])
    return out


def is_partition(p):
    return all(p[i] >= p[i + 1] for i in range(len(p) - 1)) and (not p or p[-1] >= 0)


# ----------------------------------------------------------------------------
# Schur polynomials in 3 variables, from Gelfand-Tsetlin patterns
# ----------------------------------------------------------------------------
@lru_cache(maxsize=None)
def schur3(g):
    """s_g(x1,x2,x3) as {(e1,e2,e3): coefficient}.

    Sum over Gelfand-Tsetlin patterns with top row g:
        row 2 = (m1,m2) with g1 >= m1 >= g2 >= m2 >= g3
        row 1 = (k)     with m1 >= k >= m2
    The weight of a pattern is x1^{k} x2^{(m1+m2)-k} x3^{|g|-(m1+m2)}.
    Returns {} when g is not a partition with at most 3 nonnegative parts.
    """
    g = tuple(g) + (0,) * (3 - len(g)) if len(g) < 3 else tuple(g)
    if len(g) != 3 or not is_partition(g) or g[2] < 0:
        return {}
    total = sum(g)
    out = defaultdict(int)
    for m1 in range(g[1], g[0] + 1):
        for m2 in range(g[2], g[1] + 1):
            e3 = total - m1 - m2
            for k in range(m2, m1 + 1):
                out[(k, m1 + m2 - k, e3)] += 1
    return dict(out)


def poly_mul(a, b):
    """Multiply two exponent-vector-indexed integer polynomials."""
    out = defaultdict(int)
    for ea, ca in a.items():
        for eb, cb in b.items():
            out[(ea[0] + eb[0], ea[1] + eb[1], ea[2] + eb[2])] += ca * cb
    return {e: c for e, c in out.items() if c}


def to_schur_basis(poly):
    """Rewrite a symmetric integer polynomial in 3 vars in the Schur basis.

    Repeatedly take the lexicographically largest surviving monomial.  For a
    symmetric polynomial that monomial is dominant (weakly decreasing); it must
    therefore be the leading monomial of s_that_shape, so subtract its multiple
    and continue.  Terminates because the lex-leading monomial strictly drops.
    """
    p = {e: c for e, c in poly.items() if c}
    res = defaultdict(int)
    while p:
        lead = max(p)
        if not (lead[0] >= lead[1] >= lead[2] >= 0):
            raise ValueError("leading monomial %r is not dominant: input "
                             "was not symmetric" % (lead,))
        c = p[lead]
        res[lead] += c
        for e, v in schur3(lead).items():
            nv = p.get(e, 0) - c * v
            if nv:
                p[e] = nv
            else:
                p.pop(e, None)
    return {e: c for e, c in res.items() if c}


def weyl_dim3(g):
    """dim of the GL_3 irreducible with highest weight g, by the Weyl formula.

    prod_{i<j} (g_i - g_j + j - i) / (j - i).  Exact integer arithmetic.
    """
    g = tuple(g) + (0,) * (3 - len(g)) if len(g) < 3 else tuple(g)
    num, den = 1, 1
    for i in range(3):
        for j in range(i + 1, 3):
            num *= (g[i] - g[j] + (j - i))
            den *= (j - i)
    assert num % den == 0
    return num // den


LR_STRIP_CALLS = [0, 0]              # [computed, short-circuited by degree]


def lr_strip(gamma, lam, nu):
    """c^nu_{gamma,lam} for GL_3, METHOD (a).

    Expand s_gamma * s_lam into monomials, strip leading dominant monomials, and
    read off the coefficient of s_nu.  Degree short circuit: s_gamma * s_lam is
    homogeneous of degree |gamma| + |lam| (verified separately by
    check_schur_engine), so the s_nu coefficient is 0 unless the degrees agree.
    """
    gamma = tuple(gamma)
    lam = tuple(lam)
    nu = tuple(nu)
    if sum(gamma) + sum(lam) != sum(nu):
        LR_STRIP_CALLS[1] += 1
        return 0
    sg, sl = schur3(gamma), schur3(lam)
    if not sg or not sl:
        LR_STRIP_CALLS[1] += 1
        return 0
    LR_STRIP_CALLS[0] += 1
    return to_schur_basis(poly_mul(sg, sl)).get(nu, 0)


@lru_cache(maxsize=None)
def lr_tableaux(gamma, lam, nu):
    """c^nu_{gamma,lam}, METHOD (b): count Littlewood-Richardson skew tableaux.

    Semistandard fillings of the skew shape nu/gamma (weakly increasing left to
    right, strictly increasing down columns) with content lam, whose reading
    word -- rows top to bottom, each row right to left -- is a lattice word
    (every prefix has #k >= #(k+1)).  No polynomials are used, so this is
    logically independent of schur3/to_schur_basis, and it is valid for any
    number of rows.
    """
    nu = tuple(x for x in nu if x > 0)
    gam = tuple(gamma) + (0,) * max(0, len(nu) - len(gamma))
    lm = tuple(lam)
    if len(gam) > len(nu) and any(gam[len(nu):]):
        return 0
    gam = gam[:len(nu)]
    if any(gam[i] > nu[i] for i in range(len(nu))):
        return 0
    cells = [(i, j) for i in range(len(nu))
             for j in range(nu[i] - 1, gam[i] - 1, -1)]
    if sum(lm) != len(cells):
        return 0
    m = len([x for x in lm if x > 0])
    if m == 0:
        return 1
    remaining = list(lm[:m])
    cnt = [0] * (m + 2)
    grid = {}
    total = [0]

    def rec(idx):
        if idx == len(cells):
            total[0] += 1
            return
        i, j = cells[idx]
        hi = m
        if j + 1 < nu[i]:
            hi = min(hi, grid[(i, j + 1)])
        lo = 1
        if i > 0 and gam[i - 1] <= j < nu[i - 1]:
            lo = max(lo, grid[(i - 1, j)] + 1)
        for v in range(lo, hi + 1):
            if remaining[v - 1] == 0:
                continue
            if v > 1 and cnt[v - 1] < cnt[v] + 1:
                continue
            grid[(i, j)] = v
            remaining[v - 1] -= 1
            cnt[v] += 1
            rec(idx + 1)
            cnt[v] -= 1
            remaining[v - 1] += 1
        grid.pop((i, j), None)

    rec(0)
    return total[0]


def folded_class(v, t):
    """The folded residue class of v modulo t: min(v mod t, t - (v mod t)).

    For odd t = 2m'+1 the classes are 0, 1, ..., m'; class 0 is the one Marin's
    Proposition 5.10 excludes from a transversal.
    """
    a = v % t
    return min(a, t - a)


def fold_data(lam, t=T_PARAM):
    """Everything Marin's Prop. 5.10(iii),(v) and 5.19(ii) attach to Lambda.

    Returns a dict with
        V        = 2*(Lambda + rho_{B_{R'}})  = (2L_i + (2(R'-i)+1))_i
        classes  = folded class of each V_i
        n        = {j: #indices in nonzero class j}
        vanishes = True iff some nonzero class is unhit (nu == 0)
        Smin     = the transversal of smallest V value in each nonzero class
        X        = V restricted to the complement of Smin, sorted decreasingly
        T        = the top weight, from 2*(T + rho_{D_r}) = X - t*(1,...,1)
    """
    lam = tuple(lam)
    R = len(lam)
    mp = (t - 1) // 2
    two_rho = tuple(2 * (R - 1 - i) + 1 for i in range(R))
    V = tuple(2 * lam[i] + two_rho[i] for i in range(R))
    classes = tuple(folded_class(v, t) for v in V)
    n = {}
    for j in range(1, mp + 1):
        n[j] = sum(1 for c in classes if c == j)
    out = {"lam": lam, "V": V, "classes": classes, "n": n}
    if any(n[j] == 0 for j in n):
        out.update(vanishes=True, Smin=None, X=None, T=None)
        return out
    Smin = []
    for j in range(1, mp + 1):
        idxs = [i for i in range(R) if classes[i] == j]
        Smin.append(min(idxs, key=lambda i: V[i]))
    Smin = tuple(sorted(Smin))
    X = tuple(sorted((V[i] for i in range(R) if i not in Smin), reverse=True))
    r = R - mp
    rho_D = tuple(range(r - 1, -1, -1))       # rho_{D_r} = (r-1, ..., 1, 0)
    T = []
    for k in range(r):
        num = X[k] - t - 2 * rho_D[k]
        if num % 2 != 0:
            raise ValueError("X - t*(1..1) - 2 rho is odd at %r" % (lam,))
        T.append(num // 2)
    out.update(vanishes=False, Smin=Smin, X=X, T=tuple(T))
    return out


def dr_orbit(y):
    """The W(D_r)-orbit of the integer vector y: every permutation of the
    coordinates composed with an EVEN number of sign changes.  For r = 2 this is
    {(y1,y2), (y2,y1), (-y1,-y2), (-y2,-y1)}.  These are exactly the images that
    the W(D_r)-antisymmetric extension tilde_nu of nu can be nonzero at.
    """
    from itertools import permutations
    out = set()
    r = len(y)
    for p in permutations(range(r)):
        base = tuple(y[i] for i in p)
        for signs in product((1, -1), repeat=r):
            if signs.count(-1) % 2 == 0:
                out.add(tuple(s * v for s, v in zip(signs, base)))
    return out


def numerator_index_set(lam, t=T_PARAM):
    """The set of alternant indices carried by supp nu(Lambda,.), in the doubled
    convention 2(mu + rho_{D_r}) fixed in Marin's proof: the point attached to
    (S, chirality) has index V restricted to the complement of S, sorted
    decreasingly.  At t = 3 a transversal S is a single index in the nonzero
    folded class, so the index set is obtained by deleting, one at a time, each
    entry of V that lies in that class.  (That supp nu has no indices beyond
    this set is the part of the clause still imported; see G2.)
    """
    d = fold_data(lam, t)
    V = d["V"]
    out = set()
    for i in range(len(V)):
        if folded_class(V[i], t) != 0:
            out.add(tuple(sorted((V[j] for j in range(len(V)) if j != i),
                                 reverse=True)))
    return d, out


def parts_even(gamma):
    """True iff every part of gamma is even (Littlewood's orthogonal family)."""
    return all(x % 2 == 0 for x in gamma)


def columns_even(gamma):
    """True iff every COLUMN of gamma is even, i.e. the conjugate gamma' has all
    parts even, i.e. every value occurring in gamma occurs an even number of
    times.  This is the SYMPLECTIC family; it is the decoy convention."""
    g = [x for x in gamma if x > 0]
    return all(c % 2 == 0 for c in Counter(g).values())


def a_B(lam, nu, lr=lr_strip, family=parts_even, nrows=3):
    """Littlewood's stable restriction coefficient a^B_Lambda for GL(N) -> O(N).

    a^B_Lambda = sum over partitions gamma in the even family of
    c^{nu}_{gamma,Lambda}, with nu = lambda.  gamma runs over ALL partitions with
    at most `nrows` parts and the forced size |nu| - |Lambda|; no containment of
    gamma in nu is assumed, the LR coefficient itself kills the rest.  (That
    gamma with more than 3 parts contribute nothing is checked separately.)
    """
    d = sum(nu) - sum(lam)
    if d < 0:
        return 0
    tot = 0
    for gamma in partitions_of(d, nrows, maxpart=d):
        if family(gamma):
            tot += lr(gamma, lam, nu)
    return tot


# ----------------------------------------------------------------------------
# the Weyl majorisation order  mu <= mu'  iff  mu in conv(W(C_r) mu')
# ----------------------------------------------------------------------------
def weyl_C_orbit(mu):
    """The W(C_r) = hyperoctahedral orbit of mu: all sign changes of all
    permutations of the coordinates."""
    from itertools import permutations
    out = set()
    for p in permutations(range(len(mu))):
        base = tuple(mu[i] for i in p)
        for signs in product((1, -1), repeat=len(mu)):
            out.add(tuple(s * v for s, v in zip(signs, base)))
    return out


def box_hull_weights(x, a):
    """Exact rational weights writing x in [-a,a]^2 as a convex combination of
    the four corners (+-a, +-a), by the bilinear (tensor-product) formula
        w(s1,s2) = ((1 + s1*x1/a)/2) * ((1 + s2*x2/a)/2).
    Returns {corner: Fraction}.  Nonnegativity of every weight is exactly the
    statement |x_i| <= a, so it is checked by the caller, not assumed."""
    w = {}
    for s1 in (1, -1):
        for s2 in (1, -1):
            f1 = (1 + Fraction(s1) * Fraction(x[0]) / a) / 2
            f2 = (1 + Fraction(s2) * Fraction(x[1]) / a) / 2
            w[(s1 * a, s2 * a)] = f1 * f2
    return w


def is_convex_combination(x, weights):
    """True iff the weights are nonnegative, sum to 1, and reproduce x exactly."""
    if sum(weights.values()) != 1:
        return False
    if any(v < 0 for v in weights.values()):
        return False
    dim = len(x)
    for k in range(dim):
        if sum(w * Fraction(p[k]) for p, w in weights.items()) != Fraction(x[k]):
            return False
    return True


def abs_majorized(x, mu):
    """The classical description of conv(W(C_r) mu): x is in the hull iff, for
    every k, the sum of the k largest |x_i| is at most the sum of the k largest
    |mu_i|.  Used only as an INDEPENDENT cross-check of box_hull_weights."""
    ax = sorted((abs(v) for v in x), reverse=True)
    am = sorted((abs(v) for v in mu), reverse=True)
    sx = sm = 0
    for k in range(len(ax)):
        sx += ax[k]
        sm += am[k]
        if sx > sm:
            return False
    return True


# ============================================================================
# CHECKS
# ============================================================================
def check_shape():
    """(1) The beta -> lambda dictionary, and the stable range."""
    say("--- 1. the shape ---")
    N = T_PARAM + 2 * R_PARAM
    delta = tuple(range(N - 1, -1, -1))
    strict = all(BETA[i] > BETA[i + 1] for i in range(len(BETA) - 1))
    lam_full = tuple(BETA[i] - delta[i] for i in range(N)) if len(BETA) == N else ()
    lam = tuple(x for x in lam_full if x != 0)
    say("    t = %d, r = %d  =>  N = t + 2r = %d   (len(beta) = %d)"
        % (T_PARAM, R_PARAM, N, len(BETA)))
    say("    delta            = %r" % (delta,))
    say("    beta             = %r" % (BETA,))
    say("    beta - delta     = %r  ->  lambda = %r" % (lam_full, lam))
    say("    ell(lambda) = %d, N/2 = %d/2  =>  2*ell = %d <= %d = N"
        % (len(lam), N, 2 * len(lam), N))
    say("    max(beta) = %d" % max(BETA))
    record("shape.N", N == CLAIM_N and len(BETA) == N,
           "derived N = %d, paper says %d" % (N, CLAIM_N))
    record("shape.delta", delta == CLAIM_DELTA,
           "derived %r, paper says %r" % (delta, CLAIM_DELTA))
    record("shape.beta_strict_and_zero", strict and BETA[-1] == 0,
           "beta strictly decreasing with beta_N = %d" % BETA[-1])
    record("shape.lambda", lam == CLAIM_LAMBDA and is_partition(lam_full),
           "derived lambda = %r, paper says %r" % (lam, CLAIM_LAMBDA))
    record("shape.stable_range", 2 * len(lam) <= N,
           "2*ell(lambda) = %d <= %d = N, so Littlewood's stable formula applies"
           % (2 * len(lam), N))
    record("shape.outside_Marin_box", max(BETA) > CLAIM_OBS_BOX_MAXBETA,
           "max(beta) = %d > %d, the (3,2) box of Marin's Observation 5.27"
           % (max(BETA), CLAIM_OBS_BOX_MAXBETA))
    # the two shift vectors the folding uses, from the generic formulas in
    # fold_data, against the values the paper writes out
    probe = fold_data((0, 0, 0))
    two_rho = tuple(probe["V"][i] - 0 for i in range(3))
    rho_D = tuple(range(R_PARAM - 1, -1, -1))
    say("    2*rho_{B_3} read off V(0,0,0)  = %r   (paper: %r)"
        % (two_rho, TWO_RHO_B3))
    say("    rho_{D_2} used by the formula  = %r   (paper: %r)"
        % (rho_D, RHO_D2))
    record("shape.shift_vectors", two_rho == TWO_RHO_B3 and rho_D == RHO_D2,
           "2*rho_{B_3} = %r and rho_{D_2} = %r, as the paper states"
           % (two_rho, rho_D))
    say()
    return lam


def check_schur_engine():
    """(2) The Schur/LR engines, validated against facts not used to build them."""
    say("--- 2. the Schur and Littlewood-Richardson engines ---")
    box = partitions_in_box(6, 3)
    dim_ok = hom_ok = sym_ok = 0
    from itertools import permutations
    perms = list(permutations(range(3)))
    for g in box:
        s = schur3(g)
        if sum(s.values()) == weyl_dim3(g):
            dim_ok += 1
        if all(sum(e) == sum(g) for e in s):
            hom_ok += 1
        if all(s.get(tuple(e[p[k]] for k in range(3)), 0) == c
               for e, c in s.items() for p in perms):
            sym_ok += 1
    say("    partitions in the 3 x 6 box: %d" % len(box))
    say("    s_g(1,1,1) == Weyl dimension        : %d / %d" % (dim_ok, len(box)))
    say("    s_g homogeneous of degree |g|       : %d / %d" % (hom_ok, len(box)))
    say("    s_g symmetric under S_3 on exponents: %d / %d" % (sym_ok, len(box)))
    record("engine.weyl_dimension", dim_ok == len(box),
           "%d/%d shapes" % (dim_ok, len(box)))
    record("engine.homogeneous", hom_ok == len(box),
           "%d/%d shapes; this licenses the degree short circuit in lr_strip"
           % (hom_ok, len(box)))
    record("engine.symmetric", sym_ok == len(box), "%d/%d shapes" % (sym_ok, len(box)))

    # lr_strip's degree short circuit is exercised on shapes far outside the
    # 6-box (the out-of-box gamma and Lambda sweeps reach total degree 18), so
    # homogeneity must be checked THERE, not only on the 84 shapes of the box.
    wide = [g for c in range(0, 3 * 6 + 1) for g in partitions_of(c, 3)]
    whom = sum(1 for g in wide if all(sum(e) == sum(g) for e in schur3(g)))
    wdim = sum(1 for g in wide if sum(schur3(g).values()) == weyl_dim3(g))
    say("    same three facts on every shape lr_strip is actually fed")
    say("      (<= 3 parts, |g| <= 18): %d shapes, homogeneous %d, Weyl dim %d"
        % (len(wide), whom, wdim))
    record("engine.homogeneous_wide", whom == len(wide) and wdim == len(wide),
           "homogeneity (which licenses the degree short circuit) and the Weyl "
           "dimension hold on all %d shapes of degree <= 18, the full range the "
           "sweeps use" % len(wide))

    known = [((1, 0, 0), (1, 0, 0), (2, 0, 0), 1),
             ((1, 0, 0), (1, 0, 0), (1, 1, 0), 1),
             ((2, 1, 0), (2, 1, 0), (3, 2, 1), 2),
             ((2, 1, 0), (2, 1, 0), (4, 2, 0), 1),
             ((2, 1, 0), (2, 1, 0), (2, 2, 2), 1),
             ((2, 0, 0), (1, 1, 0), (3, 1, 1), 0),
             ((3, 0, 0), (2, 2, 0), (4, 2, 1), 1)]
    bad = []
    for g, l, nu, want in known:
        a, b = lr_strip(g, l, nu), lr_tableaux(g, l, nu)
        if a != want or b != want:
            bad.append((g, l, nu, want, a, b))
        say("    c^%r_{%r,%r} = %d (strip) = %d (tableaux), expected %d"
            % (nu, g, l, a, b, want))
    record("engine.known_LR_values", not bad,
           "%d/%d textbook values reproduced by both algorithms"
           % (len(known) - len(bad), len(known)))

    triv = all(lr_strip(g, (0, 0, 0), (6, 6, 6)) == (1 if g == (6, 6, 6) else 0)
               and lr_tableaux(g, (0, 0, 0), (6, 6, 6)) == (1 if g == (6, 6, 6) else 0)
               for g in box)
    record("engine.trivial_factor", triv,
           "c^{(6,6,6)}_{gamma,(0,0,0)} = [gamma = (6,6,6)] over all %d gamma"
           % len(box))
    say()
    return box


def check_pairs(n, box, claim_pair_count=None, tag=""):
    """(3) All ordered pairs in the n-box: two algorithms, and the duality.

    Returns the full table {(gamma, Lambda): c^{(n,n,n)}_{gamma,Lambda}}.
    """
    say("--- 3%s. Littlewood-Richardson table for nu = (%d,%d,%d) ---"
        % (tag, n, n, n))
    nu = (n, n, n)
    table = {}
    disagree = dual_bad = 0
    nonzero = 0
    for gamma in box:
        for lam in box:
            a = lr_strip(gamma, lam, nu)
            b = lr_tableaux(gamma, lam, nu)
            if a != b:
                disagree += 1
                if disagree <= 5:
                    say("    ALGORITHM DISAGREEMENT c^%r_{%r,%r}: strip %d, "
                        "tableaux %d" % (nu, gamma, lam, a, b))
            table[(gamma, lam)] = a
            if a:
                nonzero += 1
            want = 1 if gamma == (n - lam[2], n - lam[1], n - lam[0]) else 0
            if a != want:
                dual_bad += 1
                if dual_bad <= 5:
                    say("    DUALITY MISMATCH c^%r_{%r,%r}: got %d, duality "
                        "predicts %d" % (nu, gamma, lam, a, want))
    sym_bad = sum(1 for k in table if table[k] != table[(k[1], k[0])])
    npairs = len(box) ** 2
    say("    ordered pairs examined              : %d  (= %d^2)" % (npairs, len(box)))
    say("    lr_strip so far (cumulative): %d full monomial expansions, %d "
        "pairs killed by degree" % (LR_STRIP_CALLS[0], LR_STRIP_CALLS[1]))
    say("    pairs with c != 0                   : %d" % nonzero)
    say("    strip vs tableaux disagreements     : %d" % disagree)
    say("    duality mismatches                  : %d" % dual_bad)
    say("    c^nu_{gamma,Lambda} != c^nu_{Lambda,gamma} : %d" % sym_bad)
    if claim_pair_count is not None:
        record("pairs%s.count" % tag, npairs == claim_pair_count,
               "derived %d ordered pairs, comparison target %d"
               % (npairs, claim_pair_count))
    record("pairs%s.two_algorithms_agree" % tag, disagree == 0,
           "monomial-strip and LR-tableau counts agree on %d/%d pairs"
           % (npairs - disagree, npairs))
    record("pairs%s.LR_symmetry" % tag, sym_bad == 0,
           "c^nu_{gamma,Lambda} = c^nu_{Lambda,gamma} on %d/%d pairs"
           % (npairs - sym_bad, npairs))
    record("pairs%s.rectangular_duality" % tag, dual_bad == 0,
           "c^{(%d,%d,%d)}_{gamma,Lambda} = [gamma = (%d-L3,%d-L2,%d-L1)] on "
           "%d/%d pairs" % (n, n, n, n, n, n, npairs - dual_bad, npairs))
    record("pairs%s.nonzero_count" % tag, nonzero == len(box),
           "exactly %d of %d pairs are nonzero, one per Lambda as the duality "
           "requires" % (nonzero, npairs))
    say()
    return table


def check_aB(n, box, table, claim_count, tag=""):
    """(4) a^B_Lambda from Littlewood's stable formula, over a superset of the box."""
    say("--- 4%s. the orthogonal restriction coefficients a^B_Lambda ---" % tag)
    nu = (n, n, n)
    deg = 3 * n
    outbox_gamma = outbox_gamma_nonzero = outbox_disagree = 0
    fourrow_gamma = fourrow_nonzero = 0
    aB = {}
    for lam in box:
        tot = 0
        for gamma in partitions_of(deg - sum(lam), 3, maxpart=deg - sum(lam)):
            if not parts_even(gamma):
                continue
            if (gamma, lam) in table:
                c = table[(gamma, lam)]
            else:
                # An out-of-box gamma.  lr_tableaux decides these by an O(1)
                # containment early return (gamma_i > nu_i => 0), which is part
                # of the DEFINITION of a skew shape and therefore no evidence at
                # all.  lr_strip has no containment test anywhere: it multiplies
                # the two Schur polynomials and reads a coefficient.  Both are
                # run, and both must give 0, so this line is real arithmetic.
                outbox_gamma += 1
                c = lr_tableaux(gamma, lam, nu)
                c2 = lr_strip(gamma, lam, nu)
                if c != c2:
                    outbox_disagree += 1
                if c or c2:
                    outbox_gamma_nonzero += 1
                c = c2
            tot += c
        for gamma in partitions_of(deg - sum(lam), 4, maxpart=deg - sum(lam)):
            if not parts_even(gamma) or gamma[3] == 0:
                continue
            fourrow_gamma += 1
            if lr_tableaux(gamma, lam, nu):
                fourrow_nonzero += 1
        aB[lam] = tot
    even_box = [l for l in box if parts_even(l)]
    nz = sorted([l for l in box if aB[l] != 0], reverse=True)
    # independent count of the all-even Lambda in the box: multisets of size 3
    # drawn from the n/2 + 1 even values 0,2,...,n
    vals = n // 2 + 1
    multiset_count = vals * (vals + 1) * (vals + 2) // 6
    say("    Lambda in the 3 x %d box                    : %d" % (n, len(box)))
    say("    all-even Lambda in the box (enumerated)     : %d" % len(even_box))
    say("    same count as multisets C(%d+2,3)           : %d"
        % (vals, multiset_count))
    say("    Lambda with a^B_Lambda != 0 (derived)       : %d" % len(nz))
    say("    values taken there                          : %s"
        % sorted(set(aB[l] for l in nz)))
    say("    even gamma outside the box that were tried  : %d (nonzero: %d, "
        "strip/tableaux disagreements: %d)"
        % (outbox_gamma, outbox_gamma_nonzero, outbox_disagree))
    say("    even 4-row gamma that were tried            : %d (nonzero: %d)"
        % (fourrow_gamma, fourrow_nonzero))
    say("    NOTE the 4-row gamma above are killed by lr_tableaux's containment")
    say("         early return, i.e. by the definition of the skew shape "
        "nu/gamma, not by")
    say("         any Littlewood-Richardson arithmetic; lr_strip cannot be used "
        "on them")
    say("         because schur3 is a 3-variable engine.  Treat that half of "
        "aB%s.gamma_box_not_assumed" % tag)
    say("         as a restatement of ell(gamma) <= ell(nu), not as independent "
        "evidence.")
    say("    the nonzero a^B_Lambda, Lambda by Lambda:")
    for lam in box:
        if aB[lam]:
            say("        %r -> %d" % (lam, aB[lam]))
    record("aB%s.count" % tag, len(nz) == claim_count,
           "derived %d admissible Lambda, paper says %d" % (len(nz), claim_count))
    record("aB%s.count_matches_multiset_formula" % tag,
           len(even_box) == multiset_count == claim_count,
           "enumerated %d, multiset formula %d" % (len(even_box), multiset_count))
    record("aB%s.support_is_even_partitions" % tag, nz == sorted(even_box, reverse=True),
           "{Lambda : a^B != 0} equals {all parts even, Lambda <= (%d,%d,%d)}"
           % (n, n, n))
    record("aB%s.all_ones" % tag, all(aB[l] == CLAIM_AB_VALUE for l in nz),
           "every admissible a^B_Lambda equals %d" % CLAIM_AB_VALUE)
    record("aB%s.gamma_box_not_assumed" % tag,
           outbox_gamma_nonzero == 0 and fourrow_nonzero == 0
           and outbox_disagree == 0 and outbox_gamma > 0,
           "all %d out-of-box (by BOTH algorithms, %d disagreements) and %d "
           "four-row even gamma contribute 0, so restricting gamma to the box "
           "was a theorem, not an assumption"
           % (outbox_gamma, outbox_disagree, fourrow_gamma))
    say()
    return aB, nz


def check_aB_outside_box(n, cap, tag=""):
    """(5) a^B_Lambda = 0 for every 3-row Lambda outside the box.

    Marin's Lambda is a dominant B_3 weight, so it is a partition with at most
    3 parts; the paper's claim that no Lambda outside the rectangle contributes
    is therefore a statement about 3-row partitions with a part exceeding n.
    |Lambda| <= |lambda| = 3n is forced by eq. (2) (gamma would need negative
    size), and a part cannot exceed the total, so sweeping parts up to 3n is
    EXHAUSTIVE.  `cap` is only the point beyond which the slower polynomial
    algorithm is dropped and the tableau count runs alone.
    """
    say("--- 5%s. a^B_Lambda vanishes outside the box ---" % tag)
    nu = (n, n, n)
    deg = 3 * n
    offenders = []
    disagree = []
    # |Lambda| <= |lambda| = 3n is forced (gamma would need negative size), so
    # deg is the largest part any candidate can have: this sweep is EXHAUSTIVE
    # over 3-row Lambda, not capped.  `cap` only marks where the second,
    # polynomial algorithm stops being run as well.
    cands = [l for l in partitions_in_box(deg, 3) if sum(l) <= deg and l[0] > n]
    for lam in cands:
        v = a_B(lam, nu, lr=lr_tableaux)
        if v:
            offenders.append((lam, v))
        if lam[0] <= cap:
            w = a_B(lam, nu, lr=lr_strip)
            if w != v:
                disagree.append((lam, v, w))
            if w and (lam, w) not in offenders:
                offenders.append((lam, w))
    say("    3-row Lambda with a part in (%d, %d] and |Lambda| <= %d : %d tried"
        % (n, deg, deg, len(cands)))
    say("    of those, also cross-checked with lr_strip (part <= %d) : %d"
        % (cap, len([l for l in cands if l[0] <= cap])))
    say("    of those with a^B_Lambda != 0                          : %d"
        % len(offenders))
    say("    tableaux/strip disagreements on the overlap            : %d"
        % len(disagree))
    if offenders:
        say("    OFFENDERS: %r" % (offenders[:10],))
    record("aB%s.zero_outside_box" % tag,
           not offenders and not disagree and len(cands) > 0,
           "a^B_Lambda = 0 on all %d out-of-box 3-row Lambda -- every one that "
           "exists, since |Lambda| <= %d bounds the parts by %d"
           % (len(cands), deg, deg))
    say()


def check_aB_many_rows(n):
    """(5b) a^B_Lambda = 0 for EVERY Lambda with four or more rows.

    The paper imports from Marin that Lambda is a dominant B_3 weight and hence
    has at most 3 parts.  For the tie set that import is not needed and this
    check removes it: Littlewood's formula itself kills every Lambda with 4 or
    more parts, because c^{(n,n,n)}_{gamma,Lambda} counts fillings of a 3-row
    skew shape and a 4th letter cannot end a strictly increasing column of
    height 3.  lr_tableaux finds that by EXHAUSTING its search -- the
    containment early return only inspects gamma -- so this is real arithmetic.
    The sweep is exhaustive: |Lambda| <= |lambda| = 3n bounds both the size and
    the number of parts.
    """
    say("--- 5b. a^B_Lambda vanishes for every Lambda with 4 or more rows ---")
    nu = (n, n, n)
    deg = 3 * n
    cands, offenders = [], []
    for rows in range(4, deg + 1):
        for s in range(rows, deg + 1):
            for lam in partitions_of(s, rows, maxpart=s):
                if lam[rows - 1] > 0:
                    cands.append(lam)
    for lam in cands:
        v = a_B(lam, nu, lr=lr_tableaux)
        if v:
            offenders.append((lam, v))
    say("    partitions with 4..%d positive parts and |Lambda| <= %d : %d tried"
        % (deg, deg, len(cands)))
    say("    of those with a^B_Lambda != 0                          : %d"
        % len(offenders))
    if offenders:
        say("    OFFENDERS: %r" % (offenders[:10],))
    record("aB.zero_for_more_than_3_rows",
           not offenders and len(cands) > 0,
           "a^B_Lambda = 0 on all %d partitions with 4 or more parts and "
           "|Lambda| <= %d -- every one that exists -- so 'Lambda has at most 3 "
           "parts' is DERIVED here for the purposes of the tie set, not imported"
           % (len(cands), deg))
    say()


def check_littlewood_convention():
    """(6) The convention in Littlewood's formula is PARTS even, not COLUMNS even.

    Sym^2(std) = o_(2) + trivial and Alt^2(std) = o_(1,1) for the orthogonal
    group; those two decompositions pin the convention down, and the
    columns-even (symplectic) variant gets both of them wrong.
    """
    say("--- 6. the orthogonal Littlewood convention ---")
    small = [l for l in partitions_in_box(2, 3) if sum(l) <= 2]
    expect = {(2, 0, 0): [(0, 0, 0), (2, 0, 0)],        # Sym^2(std)
              (1, 1, 0): [(1, 1, 0)]}                   # Alt^2(std)
    ok_parts, ok_cols = True, True
    for lam, want in expect.items():
        sup_p = sorted([l for l in small if a_B(l, lam, family=parts_even) != 0])
        sup_c = sorted([l for l in small if a_B(l, lam, family=columns_even) != 0])
        name = "Sym^2(std)" if lam == (2, 0, 0) else "Alt^2(std)"
        say("    lambda = %r  (%s)" % (lam, name))
        say("        parts-even   support: %r" % (sup_p,))
        say("        columns-even support: %r" % (sup_c,))
        say("        expected (orthogonal): %r" % (sorted(want),))
        if sup_p != sorted(want):
            ok_parts = False
        if sup_c == sorted(want):
            ok_cols = False
    record("convention.parts_even_reproduces_O3", ok_parts,
           "parts-even gives Sym^2(std) = o_(2) + trivial and "
           "Alt^2(std) = o_(1,1)")
    record("convention.columns_even_is_wrong", ok_cols,
           "the columns-even (symplectic) variant fails both, so the convention "
           "is not interchangeable")
    say()


def check_congruences(n, admissible, tag=""):
    """(7) The folded classes at t = 3, and the mod-3 impossibility."""
    say("--- 7%s. folded residue classes and the mod-3 argument ---" % tag)
    t = T_PARAM
    mp = (t - 1) // 2
    cls = {v: folded_class(v, t) for v in range(0, 4 * n + 8)}
    onlyzeroone = set(cls.values()) <= {0, 1}
    matches = all((cls[v] == 1) == (v % 3 != 0) for v in cls)
    say("    t = %d, m' = (t-1)/2 = %d nonzero folded class(es)" % (t, mp))
    say("    folded_class(v,3) takes the values %s over v in [0,%d]"
        % (sorted(set(cls.values())), max(cls)))
    say("    class 1 = {v : v = +-1 mod 3}       : %s" % matches)
    record("fold%s.single_nonzero_class" % tag, mp == 1 and onlyzeroone,
           "m' = %d, and the folded classes are %s" % (mp, sorted(set(cls.values()))))
    record("fold%s.nonzero_class_is_pm1" % tag, matches,
           "the nonzero folded class is exactly {v : v not divisible by 3}")

    # V_i = 0 mod 3 is equivalent to a congruence on L_i
    want = {0: 2, 1: 0, 2: 1}
    eq = True
    for i, two_rho_i in enumerate(TWO_RHO_B3):
        for L in range(0, 25):
            if ((2 * L + two_rho_i) % 3 == 0) != (L % 3 == want[i]):
                eq = False
    ok_vals = {i: [L for L in range(0, n + 1) if L % 2 == 0 and L % 3 == want[i]]
               for i in range(3)}
    say("    2L+5 = 0 mod 3 <=> L = 2 mod 3;  2L+3 = 0 <=> L = 0;  "
        "2L+1 = 0 <=> L = 1  : %s" % eq)
    say("    even L in [0,%d] with L = 2 mod 3 (forces L1) : %r"
        % (n, ok_vals[0]))
    say("    even L in [0,%d] with L = 0 mod 3 (forces L2) : %r"
        % (n, ok_vals[1]))
    say("    even L in [0,%d] with L = 1 mod 3 (forces L3) : %r"
        % (n, ok_vals[2]))
    record("fold%s.congruence_equivalence" % tag, eq,
           "V_i divisible by 3 is equivalent to L_i = (2,0,1)_i mod 3")
    forced = [(a, c) for a in ok_vals[0] for c in ok_vals[2] if a >= c]
    record("fold%s.congruence_impossible" % tag, not forced,
           "L1 must be in %r and L3 in %r, and no pair has L1 >= L3"
           % (ok_vals[0], ok_vals[2]))

    vanishing = [l for l in admissible if fold_data(l)["vanishes"]]
    say("    admissible Lambda whose numerator vanishes (some n_j = 0): %r"
        % (vanishing,))
    record("fold%s.no_vanishing_numerator" % tag, not vanishing,
           "all %d admissible Lambda have n_1 >= 1, so nu is not identically 0 "
           "for any of them" % len(admissible))
    say()


def check_fold_table(n, admissible, claim_M, claim_max, tag=""):
    """(8) V, X, T for every admissible Lambda; the maximizer set."""
    say("--- 8%s. the table of V, X and T ---" % tag)
    rows = []
    for lam in sorted(admissible, reverse=True):
        d = fold_data(lam)
        rows.append((lam, d["V"], d["classes"], d["n"], d["Smin"], d["X"], d["T"]))
    say("    Lambda        V              classes  n_1  S_min  X          T")
    for lam, V, cl, nn, S, X, T in rows:
        say("    %-13s %-14s %-8s %-4s %-6s %-10s %s"
            % (lam, V, cl, nn.get(1), S, X, T))
    Ts = [r[6] for r in rows]
    vbound = tuple(2 * n + k for k in (5, 3, 1))
    v_dec = all(all(V[i] > V[i + 1] for i in range(len(V) - 1)) for _, V, _, _, _, _, _ in rows)
    v_bnd = all(all(V[i] <= vbound[i] for i in range(3)) for _, V, _, _, _, _, _ in rows)
    xb = (vbound[0], vbound[1])
    x_bnd = all(X[0] <= xb[0] and X[1] <= xb[1] for _, _, _, _, _, X, _ in rows)
    le = [T for T in Ts if not (T[0] <= claim_M[0] and T[1] <= claim_M[1])]
    eq = sorted([r[0] for r in rows if r[6] == claim_M], reverse=True)
    dom_strict = all(T == claim_M or T[0] < claim_M[0] or T[1] < claim_M[1] for T in Ts)
    say("    V entries strictly decreasing everywhere      : %s" % v_dec)
    say("    V bounded coordinatewise by %-12s      : %s" % (vbound, v_bnd))
    say("    X bounded coordinatewise by %-12s      : %s" % (xb, x_bnd))
    say("    T with some coordinate above %-11s      : %r" % (claim_M, le))
    say("    Lambda with T = %-12s               : %r" % (claim_M, eq))
    # M(beta) derived, NOT read from the paper: the coordinatewise supremum of
    # the derived T's, together with the fact that it is itself attained.  The
    # paper's value is only the comparison target.
    sup = tuple(max(T[k] for T in Ts) for k in range(len(Ts[0])))
    attained = sup in Ts
    say("    coordinatewise supremum of the derived T : %r  (attained: %s)"
        % (sup, attained))
    record("fold%s.M_is_derived_not_assumed" % tag,
           attained and sup == tuple(claim_M),
           "the coordinatewise supremum of the %d derived T is %r and it is "
           "attained, so M = %r is a derived maximum; the paper says %r"
           % (len(Ts), sup, sup, tuple(claim_M)))
    say("    TAUTOLOGY WARNING: the next three checks cannot fail for any")
    say("      weakly decreasing Lambda inside the box -- V_i - V_{i+1} = "
        "2(L_i - L_{i+1}) + 2 >= 2,")
    say("      and L_i <= %d forces V_i <= %r and hence X <= %r.  They confirm "
        "the paper's" % (n, vbound, xb))
    say("      sentence 'the entries of V are strictly decreasing and bounded "
        "by %r' is true," % (vbound,))
    say("      but they are consequences of the enumeration, not independent "
        "evidence.")
    record("fold%s.V_strictly_decreasing" % tag, v_dec,
           "TAUTOLOGICAL given Lambda in the box; so the smallest value in a "
           "class carries the largest label, as Prop. 5.10(v) requires")
    record("fold%s.V_bound" % tag, v_bnd,
           "TAUTOLOGICAL given Lambda in the box; every V(Lambda) <= %r "
           "coordinatewise%s" % (vbound,
                                 "" if n != 6 else
                                 ", which is the paper's (17,15,13)"))
    record("fold%s.V_bound_matches_paper" % tag,
           n != 6 or vbound == CLAIM_V_BOUND,
           "derived bound %r%s" % (vbound,
                                   "" if n != 6 else
                                   ", paper states %r" % (CLAIM_V_BOUND,)))
    record("fold%s.X_bound" % tag, x_bnd,
           "TAUTOLOGICAL given Lambda in the box; every X(Lambda) <= %r"
           % (xb,))
    record("fold%s.T_le_M" % tag, not le,
           "all %d values T(Lambda) are <= %r coordinatewise"
           % (len(Ts), claim_M))
    record("fold%s.maximizer_set" % tag, eq == sorted(claim_max, reverse=True),
           "derived maximizers %r, paper says %r" % (eq, tuple(sorted(claim_max, reverse=True))))
    record("fold%s.maximizer_count" % tag, len(eq) == len(claim_max),
           "%d maximizers, and the conjecture's hypothesis needs more than 1"
           % len(eq))
    record("fold%s.strict_domination" % tag, dom_strict,
           "every non-maximizing T is strictly below %r in some coordinate"
           % (claim_M,))
    say()
    return rows, eq


def check_top_row(admissible):
    """(9) The Lambda = (6,6,k) family, and the (6,6,4) drop to (5,3)."""
    say("--- 9. the family Lambda = (6,6,k) ---")
    ks = sorted([l[2] for l in admissible if l[0] == 6 and l[1] == 6])
    say("    admissible k with (6,6,k) admissible : %r" % (ks,))
    record("toprow.k_set", ks == [0, 2, 4, 6],
           "derived k in %r, the paper's {0,2,4,6}" % (ks,))
    dropped = []
    for k in ks:
        d = fold_data((6, 6, k))
        V, X, T, S = d["V"], d["X"], d["T"], d["Smin"]
        div3 = [v for v in V if v % 3 == 0]
        nz = [v for v in V if v % 3 != 0]
        say("    k = %d : V = %r ; divisible by 3: %r ; nonzero class: %r ; "
            "deleted index %r (value %d) ; X = %r ; T = %r"
            % (k, V, div3, nz, S, V[S[0]], X, T))
        if T != CLAIM_M:
            dropped.append((k, V, X, T))
    d4 = fold_data((6, 6, 4))
    say("    the drop: Lambda = (6,6,4) has V = %r, of which %r are divisible "
        "by 3, so %d is the only entry in the nonzero folded class and is the "
        "one deleted." % (d4["V"], [v for v in d4["V"] if v % 3 == 0],
                          d4["V"][d4["Smin"][0]]))
    record("toprow.V_664", d4["V"] == CLAIM_V_664,
           "derived V(6,6,4) = %r, paper says %r" % (d4["V"], CLAIM_V_664))
    record("toprow.deleted_entry_664", d4["V"][d4["Smin"][0]] == 17
           and sorted(v for v in d4["V"] if v % 3 == 0) == [9, 15],
           "15 and 9 are divisible by 3, so 17 is deleted")
    record("toprow.X_664", d4["X"] == CLAIM_X_664,
           "derived X(6,6,4) = %r, paper says %r" % (d4["X"], CLAIM_X_664))
    record("toprow.T_664", d4["T"] == CLAIM_T_664,
           "derived T(6,6,4) = %r, paper says %r" % (d4["T"], CLAIM_T_664))
    record("toprow.only_664_drops", [x[0] for x in dropped] == [4],
           "of the four (6,6,k) exactly k = 4 fails to attain %r" % (CLAIM_M,))
    others = [l for l in admissible if not (l[0] == 6 and l[1] == 6)]
    hi = [l for l in others if fold_data(l)["X"] == CLAIM_X_BOUND]
    record("toprow.X_max_forces_66", not hi,
           "X(Lambda) = %r forces Lambda_1 = Lambda_2 = 6: no other admissible "
           "Lambda reaches it (%d tested)" % (CLAIM_X_BOUND, len(others)))
    say()


def check_shift_exclusion(admissible):
    """(9b) The decisive step behind the paper's eq. (4).

    The paper justifies  2(T + rho_{D_2}) = X(Lambda) - 3*(1,1)  by saying that in
    the Prop. 5.19(ii) inversion sum only the shift (1,1) contributes, every other
    shift "strictly dominates X(Lambda) coordinatewise and therefore misses the
    numerator support", and that the same argument excludes every quotient index
    above X(Lambda) - 3*(1,1).  What that paragraph does not address is that the
    summand is tilde_nu, the W(D_2)-ANTISYMMETRIC extension of nu: a reflection
    could in principle carry a shifted index back into supp nu.  Every index
    COORDINATEWISE above X(Lambda), and each of its four W(D_2) images, is tested
    here; that covers every term of the inversion sum except the shift (1,1) and
    every quotient index coordinatewise above X(Lambda) - 3*(1,1), because both
    families are of the form X(Lambda) + (p,q) with (p,q) >= 0, (p,q) != (0,0).
    The separate question of the Weyl MAJORISATION order on the top weights is
    check 10, not this one.
    """
    say("--- 9b. the shift comparison behind eq. (4), with W(D_2) images ---")
    t = T_PARAM
    dom_ok = True
    miss_ok = True
    total_offsets = 0
    total_images = 0
    hits = []
    say("    Lambda        X          supp nu indices            offsets  "
        "W(D_2) images")
    for lam in sorted(admissible, reverse=True):
        d, idx = numerator_index_set(lam, t)
        X, V = d["X"], d["V"]
        if X is None or not idx:
            # No index to compare against: the two checks below must FAIL rather
            # than quietly skip this Lambda.
            dom_ok = False
            miss_ok = False
            hits.append((lam, "no numerator index", d["classes"]))
            say("    %-13s %-10s %-26s %-8s %s"
                % (lam, X, sorted(idx, reverse=True), "-", "-"))
            continue
        vmax = max(V)
        # Prop. 5.10(v) reproduced at t = 3 from the indexing convention alone,
        # rather than quoted: X(Lambda) is one of the indices, and it dominates
        # every other one coordinatewise, strictly in at least one coordinate.
        if X not in idx:
            dom_ok = False
            hits.append((lam, "X not an index", X))
        for other in idx:
            if other == X:
                continue
            if not (all(other[k] <= X[k] for k in range(len(X)))
                    and any(other[k] < X[k] for k in range(len(X)))):
                dom_ok = False
                hits.append((lam, "not dominated", other))
        # Every other term of the inversion sum, and every quotient index above
        # X - t*(1,1), is indexed by X + (p,q) with (p,q) >= 0, (p,q) != (0,0).
        P = vmax - X[0] + 1
        Q = vmax - X[1] + 1
        n_off = 0
        n_img = 0
        for p in range(0, P + 1):
            for q in range(0, Q + 1):
                if p == 0 and q == 0:
                    continue
                n_off += 1
                Y = (X[0] + p, X[1] + q)
                for Z in dr_orbit(Y):
                    n_img += 1
                    if Z in idx:
                        miss_ok = False
                        hits.append((lam, Y, Z))
        # the two facts the tail of the enumeration rests on
        if not (min(min(z) for z in idx) >= 1
                and max(max(z) for z in idx) <= vmax
                and X[0] + P > vmax and X[1] + Q > vmax):
            miss_ok = False
            hits.append((lam, "tail bound", (P, Q, vmax)))
        total_offsets += n_off
        total_images += n_img
        say("    %-13s %-10s %-26s %-8d %d"
            % (lam, X, sorted(idx, reverse=True), n_off, n_img))
    say("    Offsets with p > max V - X_1 or q > max V - X_2 are not enumerated")
    say("    and do not need to be: such an offset has a coordinate above max V,")
    say("    that coordinate survives the coordinate permutation, and every")
    say("    coordinate of every index of supp nu lies in [1, max V]; the two")
    say("    sign-reversed images have negative coordinates, so they miss supp nu")
    say("    for the same reason.  Both facts are part of the check below.")
    record("shift.X_is_the_dominant_numerator_index", dom_ok,
           "for each of the %d admissible Lambda, X(Lambda) is one of the "
           "alternant indices of supp nu and dominates every other one "
           "coordinatewise, strictly in at least one coordinate: Marin's Prop. "
           "5.10(v) reproduced at t = 3 from the indexing convention, not quoted"
           % len(admissible))
    record("shift.weyl_images_of_higher_indices_miss_support", miss_ok,
           "%d offsets (p,q) >= 0 with (p,q) != (0,0), and all %d of their "
           "W(D_2) images, over the %d admissible Lambda: none lands in supp nu, "
           "so no term of the Prop. 5.19(ii) sum other than the shift (1,1) "
           "contributes and no quotient index coordinatewise above X - %d*(1,1) "
           "survives -- "
           "including the reflection route the paper's proof of eq. (4) leaves "
           "implicit%s" % (total_offsets, total_images, len(admissible), t,
                           "" if miss_ok else "; hits " + repr(hits[:4])))
    say()


def check_order(rows, claim_M, a, tag=""):
    """(10) The Weyl majorisation order: conv(W(C_2) M) = [-a,a]^2, and the
    maximizer set does not depend on which of the two orders is used."""
    say("--- 10%s. the Weyl majorisation order ---" % tag)
    orbit = weyl_C_orbit(claim_M)
    corners = set((s1 * a, s2 * a) for s1 in (1, -1) for s2 in (1, -1))
    say("    W(C_2)-orbit of %r : %r" % (claim_M, sorted(orbit)))
    say("    corners of [-%d,%d]^2 : %r" % (a, a, sorted(corners)))
    record("order%s.orbit_is_the_corner_set" % tag, orbit == corners,
           "the orbit has %d points and they are exactly the corners of "
           "[-%d,%d]^2, so its convex hull is that square" % (len(orbit), a, a))

    # every point of the square is a convex combination of the orbit: exact
    # rational certificate on a grid, plus the interpolation argument that the
    # certificate is an identity (each coordinate of the map is of degree <= 1
    # in each x_i, hence determined by a 2 x 2 grid of values).
    grid = [(x1, x2) for x1 in range(-a, a + 1) for x2 in range(-a, a + 1)]
    extra = [(Fraction(a, 3), Fraction(-a, 2)), (Fraction(-a), Fraction(a, 7))]
    bad = [p for p in grid + extra
           if not is_convex_combination(p, box_hull_weights(p, a))]
    say("    exact convex-combination certificates verified : %d points, "
        "%d failures" % (len(grid) + len(extra), len(bad)))
    say("    (the three functions sum(w), sum(w*p_1), sum(w*p_2) are affine in "
        "x_1 and in x_2 separately,")
    say("     so agreeing with 1, x_1, x_2 on the 2x2 grid {+-%d}^2 -- which the "
        "%d tested points contain --" % (a, len(grid)))
    say("     makes them identities, and every weight is a product of factors "
        "(1 +- x_i/%d)/2 >= 0 on the" % a)
    say("     square; hence conv(orbit) contains ALL of [-%d,%d]^2, not merely "
        "the tested points.)" % (a, a))
    ident = all(is_convex_combination(p, box_hull_weights(p, a))
                for p in [(-a, -a), (-a, a), (a, -a), (a, a)])
    record("order%s.hull_contains_square" % tag, not bad and ident,
           "every tested point of [-%d,%d]^2 -- including the 2x2 grid that "
           "determines the bilinear weight map -- is an exact convex "
           "combination of the orbit" % (a, a))
    say("    TAUTOLOGY WARNING: the next check is implied by "
        "order%s.orbit_is_the_corner_set" % tag)
    say("      above -- once the orbit IS the corner set, every orbit point "
        "trivially has")
    say("      |coordinate| = %d.  It is kept only to make the two inclusions "
        "explicit." % a)
    record("order%s.hull_inside_square" % tag,
           all(all(abs(c) <= a for c in p) for p in orbit),
           "TAUTOLOGICAL given the orbit check; each orbit point lies in "
           "[-%d,%d]^2, and that square is convex, so conv(orbit) is contained "
           "in it" % (a, a))

    # independent characterisation of the same hull, and agreement with it
    wide = [(x1, x2) for x1 in range(-a - 3, a + 4) for x2 in range(-a - 3, a + 4)]
    mism = [p for p in wide
            if abs_majorized(p, claim_M) != (abs(p[0]) <= a and abs(p[1]) <= a)]
    say("    absolute-majorisation test vs the square: %d points, %d mismatches"
        % (len(wide), len(mism)))
    record("order%s.two_hull_tests_agree" % tag, not mism,
           "the majorisation test and membership in [-%d,%d]^2 agree on all %d "
           "lattice points of [-%d,%d]^2" % (a, a, len(wide), a + 3, a + 3))

    Ts = [r[6] for r in rows]
    in_range = all(CLAIM_T_RANGE[0] <= c <= CLAIM_T_RANGE[1] for T in Ts for c in T)
    lo, hi = min(min(T) for T in Ts), max(max(T) for T in Ts)
    say("    coordinates of the derived T range over [%d, %d]" % (lo, hi))
    record("order%s.T_in_paper_range" % tag, in_range,
           "every T lies in [%d,%d]^2 as the paper states (derived range "
           "[%d,%d])" % (CLAIM_T_RANGE[0], CLAIM_T_RANGE[1], lo, hi))
    prec = [T for T in Ts if not abs_majorized(T, claim_M)]
    record("order%s.all_T_below_M" % tag, not prec,
           "T <= %r in the majorisation order for all %d admissible Lambda"
           % (claim_M, len(Ts)))
    up = [T for T in set(Ts) if T != claim_M and abs_majorized(claim_M, T)]
    record("order%s.M_is_strict_maximum" % tag, not up,
           "no other T has %r in conv(W(C_2) T), so %r is the unique maximum "
           "and the maximizer set is the same in both orders"
           % (claim_M, claim_M))
    say()


def check_tie(maximizers, aB, claim_M):
    """(11) The refutation itself: the tied sum over every sign pattern."""
    say("--- 11. the tied sum ---")
    coeffs = [aB[l] for l in maximizers]
    say("    maximizers                      : %r" % (list(maximizers),))
    say("    their a^B_Lambda                : %r" % (coeffs,))
    say("    c_top(Lambda) is only known to lie in {+1,-1} (Marin's Prop. "
        "5.10(iv)), so all %d sign patterns are enumerated:" % (2 ** len(coeffs)))
    vals = []
    for eps in product((1, -1), repeat=len(coeffs)):
        s = sum(a * e for a, e in zip(coeffs, eps))
        vals.append(s)
        say("        c_top = %-18s ->  sum = %+d   (odd: %s, zero: %s)"
            % (str(eps), s, s % 2 != 0, s == 0))
    parity = sum(coeffs) % 2
    say("    sum of the a^B_Lambda           : %d  (parity %d)"
        % (sum(coeffs), parity))
    say("    set of attainable values of the tied sum : %r" % sorted(set(vals)))
    record("tie.hypothesis_holds", len(maximizers) > 1,
           "%r is attained by %d > 1 constituents, so Conjecture 5.14 applies"
           % (claim_M, len(maximizers)))
    record("tie.parity_odd", parity == 1 and all(v % 2 != 0 for v in vals),
           "sum of the a^B_Lambda is %d, which is odd, so every value of the "
           "tied sum is odd" % sum(coeffs))
    record("tie.never_zero", all(v != 0 for v in vals),
           "the tied sum takes only the values %r; 0 is not among them, so "
           "Conjecture 5.14 is FALSE at beta = %r"
           % (sorted(set(vals)), BETA))
    say()


def check_calibration():
    """(12) CALIBRATION: the same machinery on Marin's own smallest tie witness
    beta = (8,7,6,3,2,1,0), where his Remark 5.13 reports a two-fold tie at
    T = (2,2) between Lambda = (2,2,2) and Lambda = (2,2,0)."""
    say("=" * 78)
    say("CALIBRATION on Marin's own witness beta = %r" % (CALIB_BETA,))
    say("=" * 78)
    N = T_PARAM + 2 * R_PARAM
    delta = tuple(range(N - 1, -1, -1))
    lam_full = tuple(CALIB_BETA[i] - delta[i] for i in range(N))
    lam = tuple(x for x in lam_full if x != 0)
    say("    beta - delta = %r  ->  lambda = %r" % (lam_full, lam))
    record("calib.lambda", lam == CALIB_LAMBDA and is_partition(lam_full),
           "derived lambda = %r, Marin's Remark 5.13 uses %r"
           % (lam, CALIB_LAMBDA))
    n = lam[0]
    box = partitions_in_box(n, 3)
    table = check_pairs(n, box, len(box) ** 2, tag="c")
    aB, nz = check_aB(n, box, table, CALIB_NUM_ADMISSIBLE, tag="c")
    check_congruences(n, nz, tag="c")
    rows, eq = check_fold_table(n, nz, CALIB_T, CALIB_TIE, tag="c")
    record("calib.two_fold_tie", sorted(eq, reverse=True) == sorted(CALIB_TIE, reverse=True),
           "derived tie at T = %r between %r; Marin's Remark 5.13 reports %r"
           % (CALIB_T, eq, tuple(sorted(CALIB_TIE, reverse=True))))
    record("calib.tie_is_twofold", len(eq) == 2,
           "the tie here is two-fold, so the same code that finds three "
           "maximizers at beta = %r reproduces Marin's two" % (BETA,))
    coeffs = [aB[l] for l in eq]
    sums = sorted(set(sum(a * e for a, e in zip(coeffs, eps))
                      for eps in product((1, -1), repeat=len(coeffs))))
    say("    attainable tied sums in the calibration : %r  (0 attainable: %s)"
        % (sums, 0 in sums))
    record("calib.even_tie_can_cancel", 0 in sums,
           "with two terms of +-1 the tied sum CAN vanish (values %r), which is "
           "why an even tie is consistent with Conjecture 5.14 and an odd one is "
           "not" % (sums,))
    say()


# ----------------------------------------------------------------------------
# part B: source fidelity (optional, needs the arXiv source file)
# ----------------------------------------------------------------------------
def strip_tex_comments(line):
    """Drop a LaTeX comment: the first unescaped % and everything after it."""
    out = []
    for i, ch in enumerate(line):
        if ch == "%" and (i == 0 or line[i - 1] != "\\"):
            break
        out.append(ch)
    return "".join(out)


ENV_NAMES = ("theorem", "lemma", "corollary", "proposition", "conjecture",
             "remark", "example", "problem", "question", "observation")

SOURCE_NEEDLES = [
    ("Prop 5.10(iii)",
     r"$\nu\equiv0$ if and only if some $n_j=0$ --- some nonzero residue"),
    ("Prop 5.10(iv)", r"\epsilon_t\,\det(w)\,\delta_S"),
    ("Prop 5.10(v) S_min",
     r"the index with the \emph{largest} label --- equivalently the smallest value of $V$"),
    ("Prop 5.10(v) domination",
     r"\emph{dominates coordinatewise} the index of every other point of"),
    ("the alternant indexing",
     r"by the doubled vector $2(\mu+\rho_{D_r})$, so that the point of $\nu$ attached to "
     r"$(S,\text{chirality})$ has index $V|_{S^{\mathrm c}}$ sorted decreasingly"),
    ("Prop 5.19(ii) multiplication",
     r"$a_X\mapsto\sum_{\varepsilon\in\{\pm1\}^r}(\prod_j\varepsilon_j)\,a_{X+t\varepsilon}$"),
    ("Prop 5.19(ii) inversion",
     r"\tilde c(X)\;=\;s_r\!\!\sum_{k\in(2\ZZ_{\ge0}+1)^r}\!\!\tilde\nu\bigl(X+tk\bigr)"),
    ("Conjecture 5.14 hypothesis",
     r"If $M(\beta)$ is attained by more than one $\Lambda$, then"),
    ("Conjecture 5.14 display",
     r"\sum_{\Lambda\,:\,T(\Lambda)=M(\beta)} a^B_\Lambda\,c_{\mathrm{top}}(\Lambda)\;=\;0"),
    ("Remark 5.13 population",
     r"Over the $132$ shapes at $(t,r)=(3,2),(5,2),(3,3)$"),
    ("Remark 5.13 biconditional",
     r"correct \emph{exactly} on the $105$ where the largest top is attained by a single"),
    ("Remark 5.13 tie phrase",
     r"Where it is attained twice, the question is whether the tie always cancels"),
    ("Remark 5.13 witness", r"the smallest witness is $\beta=(8,7,6,3,2,1,0)$"),
    ("Remark 5.13 two predictions",
     r"$\Lambda=(2,2,2)$ and $\Lambda=(2,2,0)$ both predict $(2,2)$"),
    ("the majorisation order",
     r"\emph{Weyl majorisation} order $\mu\preceq\mu'\iff\mu\in\mathrm{conv}\bigl(W(C_r)\,\mu'\bigr)$"),
    ("Observation 5.27 box",
     r"$\max\beta\le9$ at $(t,r)=(3,2)$ and $\max\beta\le10$ at $(5,2)$"),
    ("the 174 unchecked ties",
     r"The remaining shapes --- $174$ and $45$ --- are the ties of Remark"),
    ("verification table row",
     r"$27/27$ tied shapes; the biconditional $132/132$"),
]


# The clauses above were transcribed from the SOURCE.  Testing that they occur
# in the source is therefore close to circular: it checks a transcription, not
# the paper's fidelity.  These are the ones the paper reproduces verbatim inside
# quotation marks, and --paper checks the paper->source direction for them.
REQUIRED_IN_PAPER = (
    "Prop 5.10(iii)", "Prop 5.10(iv)", "Prop 5.10(v) S_min",
    "Prop 5.10(v) domination", "the alternant indexing",
    "Prop 5.19(ii) multiplication", "Prop 5.19(ii) inversion",
    "Conjecture 5.14 hypothesis", "Conjecture 5.14 display",
    "Remark 5.13 tie phrase", "Remark 5.13 biconditional",
    "Observation 5.27 box", "verification table row",
)


def norm_tex(s):
    """Whitespace-free, case-folded, macro-expanded form of a LaTeX fragment.

    The paper defines \\ctop and \\Ptwo and spells \\mathbb{Z} out where the
    source uses \\ZZ, and it breaks lines differently, so a byte comparison
    across the two files is hopeless.  Spacing macros are deleted and all
    whitespace removed; nothing else is touched, so the comparison is still an
    identity of the mathematical text.
    """
    s = s.replace("\\ZZ", "\\mathbb{Z}").replace("\\CC", "\\mathbb{C}")
    s = s.replace("\\ctop", "c_{\\mathrm{top}}")
    s = s.replace("\\Ptwo", "\\mathcal P^{(2)}")
    for sp in ("\\quad", "\\qquad", "\\,", "\\;", "\\:", "\\!", "\\ "):
        s = s.replace(sp, "")
    return re.sub(r"\s+", "", s.replace("~", " ")).lower()


def check_paper_quotations(paper_path):
    """(13b) The PAPER->SOURCE direction: every clause the paper puts inside
    quotation marks really occurs in Marin's file."""
    say("=" * 78)
    say("PAPER FIDELITY: the paper's quotations, checked against the source")
    say("=" * 78)
    with open(paper_path, "r", encoding="utf-8", errors="replace") as fh:
        ptxt = norm_tex(" ".join(strip_tex_comments(l)
                                 for l in fh.read().split("\n")))
    found, absent = [], []
    for name, needle in SOURCE_NEEDLES:
        hit = norm_tex(needle) in ptxt
        (found if hit else absent).append(name)
        say("    %-32s in the paper: %s%s"
            % (name, hit, "" if hit or name not in REQUIRED_IN_PAPER
               else "   <-- REQUIRED"))
    req_missing = [n for n in REQUIRED_IN_PAPER if n in absent]
    say("    clauses of Marin's text reproduced by the paper : %d of %d"
        % (len(found), len(SOURCE_NEEDLES)))
    say("    of the %d load-bearing ones, missing            : %r"
        % (len(REQUIRED_IN_PAPER), req_missing))
    record("paper.load_bearing_quotes_are_faithful", not req_missing,
           "all %d clauses the refutation actually leans on -- Prop. 5.10"
           "(iii),(iv),(v), the alternant indexing, Prop. 5.19(ii), the "
           "conjecture itself, Remark 5.13 and Observation 5.27 -- appear in "
           "BOTH files, so the paper is not quoting a statement it invented"
           % len(REQUIRED_IN_PAPER))
    say("    NOTE the %d unrequired clauses are ones the paper paraphrases "
        "rather than quotes" % (len(SOURCE_NEEDLES) - len(REQUIRED_IN_PAPER)))
    say()


def check_source(path):
    """(13) SOURCE FIDELITY: that the refuted statement really is Conjecture 5.14
    of arXiv:2608.18302v1, and that the clauses the paper quotes are there."""
    say("=" * 78)
    say("SOURCE FIDELITY against %s" % path)
    say("=" * 78)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        raw = fh.read().split("\n")
    lines = [strip_tex_comments(l) for l in raw]
    say("    lines in the source                    : %d" % len(raw))

    # (a) the numbered environments of Section 5 share one section counter
    decl_re = re.compile(r"\\newtheorem\*?\{(\w+)\}(?:\[(\w+)\])?"
                         r"\{[^}]*\}(?:\[(\w+)\])?")
    decls = {}
    base = None
    for i, l in enumerate(lines, 1):
        for m in decl_re.finditer(l):
            name, alias, numbering = m.group(1), m.group(2), m.group(3)
            decls[name] = (i, alias, numbering)
            if numbering == "section" and alias is None:
                base = name
    shared = [e for e in ENV_NAMES
              if e in decls and (e == base or decls[e][1] == base)]
    say("    \\newtheorem declarations found          : %d, at lines %d-%d"
        % (len(decls), min(v[0] for v in decls.values()),
           max(v[0] for v in decls.values())))
    say("    the [section]-numbered base counter     : %r" % base)
    say("    environments sharing that counter       : %r" % shared)
    record("source.one_shared_counter",
           base == "theorem" and sorted(shared) == sorted(ENV_NAMES),
           "all %d environment names share the [section]-numbered counter %r"
           % (len(ENV_NAMES), base))
    # The "14th environment => 5.14" inference is only valid if NOTHING ELSE
    # advances that counter.  A single \newtheorem{definition}[theorem]{...}
    # used in Section 5 would break the numbering while every check above still
    # passed, so the set of users of the counter has to be closed, and none of
    # the counted names may be declared with \newtheorem* (unnumbered).
    star_re = re.compile(r"\\newtheorem\*\{(\w+)\}")
    starred = sorted(set(m.group(1) for l in lines for m in star_re.finditer(l)))
    extra_sharing = sorted(e for e in decls if e not in ENV_NAMES
                           and (e == base or decls[e][1] == base))
    say("    all declared environment names          : %r" % sorted(decls))
    say("    unnumbered (\\newtheorem*) declarations  : %r" % starred)
    say("    other environments on the same counter  : %r" % extra_sharing)
    record("source.counter_has_no_other_users",
           not extra_sharing and not [e for e in starred if e in ENV_NAMES]
           and sorted(decls) == sorted(ENV_NAMES),
           "the only environments driving counter %r are the %d counted names "
           "and none is declared unnumbered, so no environment invisible to "
           "this count can advance the numbering" % (base, len(ENV_NAMES)))
    resets = [(i, l.strip()) for i, l in enumerate(lines, 1)
              if "\\setcounter" in l or "\\numberwithin" in l or "\\section*" in l]
    say("    \\setcounter / \\numberwithin / \\section* : %d" % len(resets))
    record("source.no_counter_resets", not resets,
           "no \\setcounter, \\numberwithin or \\section* anywhere, so the "
           "counter runs straight through")

    # A [section]-numbered counter is ALSO reset by every ordinary \section.
    # "the 14th environment after the start of Section 5 is 5.14" therefore needs
    # (a) that the recorded start of Section 5 lies in the FIFTH section, and
    # (b) that no further \section intervenes between it and the conjecture's own
    # \begin.  Neither was implied by anything above.
    sec_re = re.compile(r"\\section\s*(?:\[|\{)")
    defish = ("\\renewcommand", "\\newcommand", "\\providecommand", "\\def",
              "\\let", "\\@startsection")
    secs, secdefs = [], []
    for i, l in enumerate(lines, 1):
        if not sec_re.search(l):
            continue
        if any(d in l for d in defish):
            secdefs.append((i, l.strip()))
        else:
            secs.append(i)
    for i, l in secdefs:
        say("    \\section is REDEFINED at line %d: %s" % (i, l[:78]))
    # A cosmetic \renewcommand\section{\@startsection{section}{1}...} keeps the
    # `section' counter and its level, so the numbering is unaffected; anything
    # else would have to be read by hand, so it is recorded as a check.
    record("source.section_redefinition_is_cosmetic",
           all("\\@startsection{section}{1}" in l.replace(" ", "")
               or "\\@startsection{section}{1}" in l for i, l in secdefs),
           "%d redefinition(s) of \\section, all of the \\@startsection"
           "{section}{1} form, which keeps the `section' counter and level"
           % len(secdefs))
    before = [i for i in secs if i <= CLAIM_SEC5_LINE]
    between = [i for i in secs if CLAIM_SEC5_LINE < i <= CLAIM_CONJ_LINE]
    say("    \\section lines (first 12)               : %r" % (secs[:12],))
    say("    \\section commands at or before line %d : %d (last at %r)"
        % (CLAIM_SEC5_LINE, len(before), before[-1] if before else None))
    say("    \\section commands in (%d, %d]         : %r"
        % (CLAIM_SEC5_LINE, CLAIM_CONJ_LINE, between))
    record("source.section5_is_the_fifth_section", len(before) == 5,
           "line %d lies inside the 5th \\section (%d found at or before it), "
           "which is the '5' of 'Conjecture 5.14'"
           % (CLAIM_SEC5_LINE, len(before)))
    record("source.no_section_break_before_conjecture", not between,
           "no \\section occurs between line %d and line %d, so the section "
           "counter does not reset and the 14th environment after line %d is "
           "numbered 5.14 rather than 6.k"
           % (CLAIM_SEC5_LINE, CLAIM_CONJ_LINE, CLAIM_SEC5_LINE))

    # (b) count the numbered environments after the start of Section 5
    seq = []
    for i, l in enumerate(lines, 1):
        if i <= CLAIM_SEC5_LINE:
            continue
        for env in ENV_NAMES:
            k = l.find("\\begin{" + env + "}")
            while k != -1:
                seq.append((i, env, l))
                k = l.find("\\begin{" + env + "}", k + 1)
    seq.sort(key=lambda z: z[0])
    say("    numbered environments after line %d     : %d"
        % (CLAIM_SEC5_LINE, len(seq)))
    if len(seq) >= CLAIM_CONJ_INDEX:
        ln, env, text = seq[CLAIM_CONJ_INDEX - 1]
        say("    the %dth of them                        : %s at line %d"
            % (CLAIM_CONJ_INDEX, env, ln))
        say("        %s" % text.strip()[:100])
        record("source.conjecture_is_14th",
               env == "conjecture" and ln == CLAIM_CONJ_LINE
               and CLAIM_CONJ_LABEL in text,
               "environment #%d of Section 5 is a %s at line %d carrying "
               "\\label{%s}, i.e. Conjecture 5.%d"
               % (CLAIM_CONJ_INDEX, env, ln, CLAIM_CONJ_LABEL, CLAIM_CONJ_INDEX))
    else:
        record("source.conjecture_is_14th", False,
               "only %d environments found after line %d"
               % (len(seq), CLAIM_SEC5_LINE))

    # (b2) the WHOLE body of the conjecture, not just the presence of the two
    # clauses the paper quotes.  A substring test passes even if the source's
    # conjecture carries an extra hypothesis that beta = (12,11,10,3,2,1,0)
    # might fail, which would make the "counterexample" vacuous.  So the body is
    # printed verbatim and, after the quoted clauses and the LaTeX scaffolding
    # are deleted, must contain no alphabetic residue at all.
    if CLAIM_CONJ_LINE <= len(lines):
        body, j = [], CLAIM_CONJ_LINE - 1
        while j < len(lines) and len(body) < 60:
            body.append(lines[j])
            if "\\end{conjecture}" in lines[j]:
                break
            j += 1
        btxt = " ".join(" ".join(body).split())
        say("    Conjecture 5.%d verbatim, whole environment:" % CLAIM_CONJ_INDEX)
        say("        %s" % btxt)
        red = re.sub(r"\\begin\{conjecture\}(\[[^\]]*\])?(\\label\{[^}]*\})?",
                     "", btxt).replace("\\end{conjecture}", "")
        for _nm, needle in SOURCE_NEEDLES:
            nn = " ".join(needle.split())
            if nn in red:
                red = red.replace(nn, "")
        leftover = re.sub(r"[^A-Za-z]", "", red)
        say("    alphabetic residue after deleting the quoted clauses : %r"
            % leftover)
        record("source.conjecture_has_no_extra_hypothesis",
               leftover == "" and "\\end{conjecture}" in btxt,
               "the entire body of Conjecture 5.%d is the hypothesis and the "
               "display the paper quotes -- no alphabetic residue -- so there "
               "is no further hypothesis that beta could fail"
               % CLAIM_CONJ_INDEX)

    # (c) every clause the paper quotes is present
    flat = " ".join(" ".join(lines).split())
    missing = [name for name, needle in SOURCE_NEEDLES
               if " ".join(needle.split()) not in flat]
    for name, needle in SOURCE_NEEDLES:
        if " ".join(needle.split()) not in flat:
            say("    MISSING: %s" % name)
    say("    quoted clauses located in the source    : %d of %d"
        % (len(SOURCE_NEEDLES) - len(missing), len(SOURCE_NEEDLES)))
    record("source.quoted_clauses_present", not missing,
           "%d of %d clauses quoted by the paper are present verbatim in the "
           "source%s" % (len(SOURCE_NEEDLES) - len(missing), len(SOURCE_NEEDLES),
                         "" if not missing else "; missing " + repr(missing)))
    say()


def print_gap_register(source_given, paper_given=False):
    """(14) THE GAP REGISTER: everything between the checked facts and the
    paper's conclusion that no check above covers.  Printed unconditionally so
    that a green verdict cannot be mistaken for a complete proof."""
    say("=" * 78)
    say("GAP REGISTER -- what a green verdict above does NOT establish")
    say("=" * 78)
    say("G1 IMPORTED, NOT DERIVED: c_top(Lambda) = c(Lambda,T(Lambda)) lies in")
    say("   {+1,-1}  (Marin Prop. 5.10(iv)).  This program does not build")
    say("   tau^B_t, the alternant expansion or the Delta_t division, so it")
    say("   cannot confirm the membership.  MITIGATION: all 2^3 sign patterns")
    say("   are enumerated, so the refutation is independent of WHICH signs")
    say("   occur -- but if any c_top were 0, or +-2, the argument fails.")
    say("G2 PARTLY IMPORTED: the rule that reads T(Lambda) off V(Lambda) (Prop.")
    say("   5.10(v) plus the 5.19(ii) indexing).  DERIVED HERE (check 9b), once")
    say("   supp nu is indexed as the source's proof fixes it -- V restricted to")
    say("   the complement of a transversal, sorted decreasingly, which at t = 3")
    say("   means deleting one entry of V from the nonzero folded class: that")
    say("   X(Lambda) is the coordinatewise dominant index (the content of Prop.")
    say("   5.10(v)), and that every index COORDINATEWISE above X(Lambda) -- so")
    say("   every term of the inversion sum but the shift (1,1) -- misses supp nu")
    say("   TOGETHER WITH ALL FOUR OF ITS W(D_2) IMAGES.  That last clause is the")
    say("   reflection route back into the support which the paper's four-")
    say("   sentence justification of eq. (4) leaves implicit: the summand of the")
    say("   5.19(ii) inversion is the W(D_2)-ANTISYMMETRIC extension tilde_nu, so")
    say("   coordinatewise domination alone would not settle it.")
    say("   STILL IMPORTED: the 5.19(ii) inversion formula itself, the global")
    say("   sign s_r, and that supp nu contains no index outside the deletion set")
    say("   just described.  The only external test of the composite rule remains")
    say("   the calibration: Marin's own Remark 5.13 reports a two-fold tie at")
    say("   T = (2,2) for beta = (8,7,6,3,2,1,0), and the same code reproduces")
    say("   it.  One data point.")
    say("G3 NO LONGER IMPORTED (closed by check 5b): 'Lambda is a dominant B_3")
    say("   weight, hence has at most 3 parts' is not needed -- every partition")
    say("   with 4 or more parts and |Lambda| <= 18 is swept and has")
    say("   a^B_Lambda = 0, so no such Lambda can join the tie.  What remains")
    say("   imported is only that Lambda is a PARTITION (dominant, integral),")
    say("   i.e. that no weight with a negative or fractional coordinate enters")
    say("   the branching sum.")
    say("G4 IMPORTED: Littlewood's stable orthogonal restriction formula and")
    say("   its validity for 2*ell(lambda) <= N.  Only its parts-even vs")
    say("   columns-even CONVENTION is discriminated here (check 6).")
    say("G5 IMPORTED: M(beta) = max{T(Lambda) : a^B_Lambda != 0}, the maximum")
    say("   restricted to the SUPPORT of a^B.  That restriction is load-")
    say("   bearing, and the next check shows it:")
    d = fold_data((8, 8, 8))
    say("       Lambda = (8,8,8) is a perfectly good 3-row even partition with")
    say("       V = %r, X = %r, T = %r, which exceeds M = %r coordinatewise."
        % (d["V"], d["X"], d["T"], CLAIM_M))
    say("       It is excluded only because a^B_(8,8,8) = 0 (it is not inside")
    say("       lambda = (6,6,6)).  So if Marin's max ran over all Lambda")
    say("       rather than the support, M(beta) would not be (6,6).")
    record("gap.support_restriction_is_load_bearing",
           d["T"] is not None and (d["T"][0] > CLAIM_M[0] or d["T"][1] > CLAIM_M[1])
           and a_B((8, 8, 8), CLAIM_LAMBDA, lr=lr_tableaux) == 0,
           "T(8,8,8) = %r exceeds M = %r, and a^B_(8,8,8) = 0, so M(beta) = %r "
           "depends on the support restriction, which is imported from Marin's "
           "definition and not derivable here" % (d["T"], CLAIM_M, CLAIM_M))
    say("G6 NOT PERFORMED: the program makes no network call of its own; it")
    say("   uses only the standard library, and so establishes NOTHING about which")
    say("   versions of arXiv:2608.18302 exist, which is current, or what any")
    say("   of them contains.  No claim about that is made here.  Hence 'v1 is")
    say("   the version refuted' rests entirely on the file passed to --source")
    say("   being the v1 e-print: the program checks that file's CONTENT, not")
    say("   its provenance, and later versions may exist and may renumber or")
    say("   restate the conjecture.  Pin the provenance yourself -- fetch the")
    say("   version you mean, read the submission history on the abs page for")
    say("   the list of versions, and re-run part B on the file you fetched:")
    say("       curl -sL https://arxiv.org/e-print/2608.18302v1 | tar xz")
    say("       python3 verify.py --source orbit_pair_ii.tex")
    if not source_given:
        say("G7 NOT PERFORMED: part B (source fidelity) was skipped -- no")
        say("   --source file was given, so NOTHING above ties the refuted")
        say("   statement to Conjecture 5.14 of arXiv:2608.18302v1.")
    if not paper_given:
        say("G7b NOT PERFORMED: the paper->source direction (--paper).  The 18")
        say("   clause strings in SOURCE_NEEDLES were transcribed FROM the")
        say("   source, so source.quoted_clauses_present on its own only shows")
        say("   that a transcription of the source is in the source.  Only")
        say("   --paper turns it into a two-file identity.")
    say("G8 TAUTOLOGIES, labelled in place, that cannot fail for any input:")
    say("   fold.V_strictly_decreasing, fold.V_bound, fold.X_bound,")
    say("   order.hull_inside_square, and the 4-row half of")
    say("   aB.gamma_box_not_assumed.  They confirm true sentences of the")
    say("   paper; they are not evidence.")
    say()


# ============================================================================
# main
# ============================================================================
def main(argv):
    source = None
    paper = None
    for i, a in enumerate(argv):
        if a == "--source" and i + 1 < len(argv):
            source = argv[i + 1]
        elif a.startswith("--source="):
            source = a.split("=", 1)[1]
        elif a == "--paper" and i + 1 < len(argv):
            paper = argv[i + 1]
        elif a.startswith("--paper="):
            paper = a.split("=", 1)[1]

    say("=" * 78)
    say("verify.py -- a counterexample to Conjecture 5.14 of arXiv:2608.18302v1")
    say("beta = %r at (t,r) = (%d,%d)" % (BETA, T_PARAM, R_PARAM))
    say("exact integer and rational arithmetic only; no floating point")
    say("=" * 78)
    say()

    check_shape()
    box = check_schur_engine()
    table = check_pairs(6, box, CLAIM_PAIR_COUNT)
    aB, nz = check_aB(6, box, table, CLAIM_NUM_ADMISSIBLE)
    check_aB_outside_box(6, cap=9)
    check_aB_many_rows(6)
    check_littlewood_convention()
    check_congruences(6, nz)
    rows, maximizers = check_fold_table(6, nz, CLAIM_M, CLAIM_MAXIMIZERS)
    check_top_row(nz)
    check_shift_exclusion(nz)
    check_order(rows, CLAIM_M, CLAIM_HULL_BOX)
    check_tie(maximizers, aB, CLAIM_M)
    check_calibration()

    if paper:
        try:
            check_paper_quotations(paper)
        except IOError as exc:
            record("paper.file_readable", False, "cannot read %r: %s"
                   % (paper, exc))
    else:
        say("WARNING: the paper->source direction was not checked.  Pass")
        say("      --paper paper.tex  to verify that the clauses "
            "the paper puts in")
        say("      quotation marks really occur in Marin's file.  Without it, "
            "part B only shows")
        say("      that strings TRANSCRIBED FROM THE SOURCE occur in the "
            "source, which is nearly")
        say("      circular.")
        say()

    if source:
        try:
            check_source(source)
        except IOError as exc:
            record("source.file_readable", False, "cannot read %r: %s"
                   % (source, exc))
    else:
        say("WARNING: part B (source fidelity) was not run.  It needs the single "
            "file orbit_pair_ii.tex of")
        say("      arXiv:2608.18302v1, obtainable with")
        say("          curl -sL https://arxiv.org/e-print/2608.18302v1 | tar xz")
        say("      and is then enabled with  --source orbit_pair_ii.tex .  This "
            "program makes no")
        say("      network call of its own, so the %d checks above are the "
            "mathematics alone." % len(RESULTS))
        say()

    print_gap_register(bool(source), bool(paper))

    failed = [n for n, ok in RESULTS if not ok]
    say("=" * 78)
    if failed:
        say("failed checks: %s" % ", ".join(failed))
        say("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), len(RESULTS)))
        return 1
    say("VERDICT: ALL %d CHECKS PASS%s" % (len(RESULTS),
        "" if (source and paper) else " (INCOMPLETE -- see G6/G7 above)"))
    say("        This is conditional.  Three clauses (G1, G4, G5) are imported")
    say("        from Marin and are NOT verified here, and a fourth (G2) only")
    say("        in part -- see the register; the version-")
    say("        provenance probe (G6) is not performed%s; and five checks are"
        % ("" if (source and paper) else ", part B is incomplete"))
    say("        tautological (G8).  Read the GAP REGISTER before quoting this.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

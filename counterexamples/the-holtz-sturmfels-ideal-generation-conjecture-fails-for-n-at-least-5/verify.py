#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
verify.py -- referee's verification program for

    "The Holtz-Sturmfels Ideal-Generation Conjecture Fails for n >= 5"

The paper refutes Holtz-Sturmfels [HS, Conjecture 14] (J_n = P_n) by exhibiting
four nonzero SL_2^5-invariant quartics F_2,...,F_5 in (P_5)_4 that cannot lie in
HD(5) = (J_5)_4, because HD(5) carries no nonzero SL_2^5-invariant.

======================================================================
VALUES TAKEN FROM THE PAPER (inputs; these are NOT checks)
======================================================================
 I1  n = 5; coordinates A_alpha, alpha in {0,1}^5, identified with A_I, I <= [5].
 I2  The alternating form  eps_01 = 1, eps_10 = -1, eps_00 = eps_11 = 0.
 I3  The definition q^{(i)}_{ab}(A) = sum_{alpha_i=a, beta_i=b} A_alpha A_beta
     prod_{j != i} eps_{alpha_j beta_j},  Q_i = (q^{(i)}_{ab}),  B_i = det Q_i.
 I4  The homogeneous principal-minor map phi_n[X:t] = [ t^{n-|I|} det X_I ]_I.
 I5  The 12 signed monomials of the pentad Pi(X) as printed in Section 3.
 I6  F_i = B_1 - B_i for i = 2,3,4,5.
 I7  The 4 x 4 coefficient matrix printed in Lemma 3 (rows = the four monomials
     m_1..m_4 named there, columns = F_2,F_3,F_4,F_5) and its determinant -512.
 I8  Oeding's G-module decomposition of HD(n) [Oeding, Prop. 5.2]:
     HD(n) = sum over T <= [n], |T| = 3, of
        (tensor_{j in T} S_{(2,2)} V_j^*) tensor (tensor_{j notin T} S_{(4)} V_j^*).
 I9  Oeding's set-theoretic theorem sqrt(J_n) = P_n [Oeding, Thm 1.3], and
     HD(n) <= (P_n)_4 [Oeding, Prop. 5.5]; J_n is generated in degree four.
 I10 Cayley's 2x2x2 hyperdeterminant, in its standard 12-term expansion, as the
     generator whose GL_2^n-orbit spans HD(n).

======================================================================
WHAT THIS PROGRAM DERIVES (the checks)
======================================================================
 * Q_i is symmetric, B_i is a nonzero quartic of subset weight 10, the five B_i
   are pairwise distinct -- all recomputed from I2/I3 alone.
 * The pentad printed in the paper (I5) is re-derived independently as the signed
   sum over the 12 Hamiltonian cycles of K_5, sign = parity of the cycle word.
 * THE LOAD-BEARING IDENTITY: B_i(phi_5(X,t)) - 12 t^10 Pi(X)^2 = 0
   coefficientwise in Z[x_11..x_55, t] (16 variables), for every i = 1..5, by
   exact expansion; hence F_i(phi_5(X,t)) == 0 and F_i in (P_5)_4.
 * B_i(g.A) = (prod_j det g_j)^2 B_i(A) for integer g in GL_2^5 -- so the F_i are
   SL_2^5-invariant, and the test is shown to discriminate (det != 1 rescales).
 * The 4 x 4 coefficient matrix of I7 is recomputed from the actual F_i, its
   determinant is recomputed (-512), and the rank of the full 4 x 280 integer
   coefficient matrix is recomputed as 4.
 * dim S_(2,2) C^2 = 1 and dim S_(4) C^2 = 5 by counting semistandard tableaux;
   S_lambda C^2 restricted to SL_2 is Sym^{lam1-lam2} by comparing characters;
   dim Sym^m(C^2)^{SL_2} = [m = 0] by exact kernels of two unipotents. Hence,
   from I8: dim HD(n) = binom(n,3) 5^{n-3} and HD(n)^{H_S} = 0.
 * HD(5) is not merely cited.  Because W is SL_2^5-trivial, W can only meet the
   TORUS-WEIGHT-ZERO part of HD(5); that part is built explicitly, as the span
   of 130 elements got by pushing all 40 hyperdeterminant slices to weight zero
   with the sl_2^5 operators and with transvection words.  Its exact rank is 10,
   the dimension predicted from I8; adjoining F_2..F_5 raises the rank to 14
   (so W cap HD(5) = 0), and the joint kernel of the ten sl_2^5 operators on
   that 10-dimensional space is 0 (so HD(5)^{G_5} = 0), both computed.
 * The matching rank 10 alone would only be a coincidence test against the cited
   decomposition, so the UPPER bound is computed as well: for each of the ten
   3-subsets T, the weight-zero projection of the sampled GL_2^5-orbit of the
   four slices D_{T,fix} is found to be exactly ONE dimensional and saturated,
   whence dim HD(5)^{weight 0} <= 10 with no appeal to I8; adjoining all sampled
   orbit elements to the 130 constructed ones does not raise the rank above 10.
 * The refutation is then restated in a form that quotes no number from the
   literature: W <= (P_5)_4 (computed), dim W = 4 (computed), W cap HD(5) = 0
   (computed) and (J_5)_4 = HD(5) (input I9) already give J_5 != P_5, whatever
   dim HD(5) is; the values 250 and 254 are transcription checks, not premises.
 * The F_i are annihilated by all ten generators of sl_2^5 -- infinitesimal
   invariance, a second and independent route to Lemma 1.
 * dim(P_5)_4 >= 254 > 250 = dim HD(5), so J_5 != P_5.
 * Census: the n >= 5 reduction is re-run exactly for n = 6..12, over EVERY
   5-subset S and every i, at several exact integer points; the HD(n) module
   facts are re-run for n = 3..14 over every (3-subset, 5-subset) pair.
 * Cayley's hyperdeterminant of each of the forty 2x2x2 slices of phi_5(X,t) is
   shown to vanish EXACTLY, as the zero polynomial of Z[x_11..x_55, t] -- the
   computational content of "HD(5) <= (P_5)_4" (I9) at the orbit generators, not
   a sample.  The numeric GL_2^5-twisted slices and a generic-array control that
   Cayley is not identically zero are retained beside it.
 * Negative controls: thirteen deliberate corruptions of the exhibited object are
   applied and each is confirmed to make a load-bearing check fail.  Controls of
   the form "corrupting X breaks a check" are demonstrations that a test
   discriminates; they are not, and are not counted as, verifications.

Standard library only; exact integer arithmetic throughout; no floating point.
Every load-bearing algebraic identity is verified coefficientwise over Z, not at
sampled points.  Sampling is used in exactly three places, and in each the
sampled statement is also established exactly by another check:
  * the finite-group invariance of B_i and of W (integer GL_2^5 / SL_2(Z) tuples)
    -- also proved exactly, and infinitesimally, by the sl_2^5 annihilation check;
  * the numeric evaluations of the pullback identity and of the n = 6..12 census
    -- the n = 5 identity itself is symbolic and the census is exhaustive over
    the discrete data (n, S, i);
  * the SATURATION bound dim HD(5)^{weight 0} <= 10 of check 27b, which samples
    the GL_2^5-orbit of each slice.  A saturated sample is evidence, not proof,
    of an upper bound; the same value 10 is derived independently from Oeding's
    decomposition (input I8), and the two agree.  Any nonzero orbit element the
    sample missed would only enlarge HD(5)^{weight 0}, so the check is a genuine
    guard against the cited decomposition understating HD(5): if it ever raised
    the rank above 10, the separation verdict would fail rather than pass.

THE REAL FORM.  The paper works over C, while Conjecture 14 is posed over R
(ring R[A_.], group SL_2(R)^n), and the paper's Status paragraph descends the
refutation to that form.  Everything this program contributes to that descent is
arithmetic and not a separate argument: every polynomial built here -- the B_i,
the F_i, the pentad, Cayley's slices, the F_i^S of the census, and every
constructed or sampled generator of HD(5)^{weight 0} -- has integer
coefficients, and every rank and kernel above is computed over Q by exact
elimination, so each of those ranks is also the rank of the same integer matrix
over R.  The step of the descent that is taken from the paper and NOT
re-derived here is the one about the module itself: that HD(n), the span of the
orbit of an integral quartic under a group defined over Q, has a basis of
quartics with integer coefficients.  The closing NOT RE-RUN block records this
shortfall, with the coefficient scan reported as a derived count.
"""

import itertools
import math
import random
import sys
from fractions import Fraction

RESULTS = []


def record(name, ok, detail=""):
    """Register one check.  Printing happens here so order is stable."""
    RESULTS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    sys.stdout.flush()
    return bool(ok)


def section(title):
    print("")
    print("--- %s ---" % title)


# ======================================================================
# Exact sparse multivariate polynomials over Z.
# A monomial is an exponent vector packed into one Python int, six bits per
# variable (exponents stay well under 63), so monomial multiplication is a
# single integer addition.  Coefficients are Python ints: exact, unbounded.
# ======================================================================

SH = 6                                   # bits per exponent field


def pk(pairs):
    """Pack an iterable of (variable index, exponent) into a monomial key."""
    m = 0
    for k, e in pairs:
        m += e << (SH * k)
    return m


def unpk(m, nv):
    """Unpack a monomial key into a tuple of nv exponents."""
    mask = (1 << SH) - 1
    return tuple((m >> (SH * k)) & mask for k in range(nv))


def pmul(P, Q):
    R = {}
    for m1, c1 in P.items():
        for m2, c2 in Q.items():
            m = m1 + m2
            v = R.get(m)
            R[m] = c1 * c2 if v is None else v + c1 * c2
    return dict((m, c) for m, c in R.items() if c)


def pscale(P, s):
    if s == 0:
        return {}
    return dict((m, c * s) for m, c in P.items())


def padd(P, Q, s=1):
    R = dict(P)
    for m, c in Q.items():
        v = R.get(m, 0) + s * c
        if v:
            R[m] = v
        elif m in R:
            del R[m]
    return R


def pshift(P, k, e):
    """Multiply P by (variable k)^e."""
    if e == 0:
        return dict(P)
    t = 1 << (SH * k)
    t *= e
    return dict((m + t, c) for m, c in P.items())


def permsign(w):
    """Sign of the word w read as a permutation (parity of inversions)."""
    s = 1
    n = len(w)
    for a in range(n):
        for b in range(a + 1, n):
            if w[a] > w[b]:
                s = -s
    return s


def det_int(M):
    """Exact determinant of a square integer matrix, by Leibniz (small n)."""
    n = len(M)
    if n == 0:
        return 1
    tot = 0
    for p in itertools.permutations(range(n)):
        pr = permsign(p)
        for k in range(n):
            pr *= M[k][p[k]]
        tot += pr
    return tot


def rank_exact(rows):
    """Exact rank over Q of an integer matrix given as a list of rows."""
    A = [[Fraction(x) for x in r] for r in rows]
    if not A:
        return 0
    nc = len(A[0])
    r = 0
    for c in range(nc):
        piv = None
        for k in range(r, len(A)):
            if A[k][c] != 0:
                piv = k
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        for k in range(len(A)):
            if k != r and A[k][c] != 0:
                f = A[k][c] / pv
                A[k] = [a - f * b for a, b in zip(A[k], A[r])]
        r += 1
        if r == len(A):
            break
    return r


# ======================================================================
# INPUTS TAKEN FROM THE PAPER
# ======================================================================

N = 5                                            # I1

EPS = {(0, 1): 1, (1, 0): -1, (0, 0): 0, (1, 1): 0}          # I2

# I5: the twelve signed monomials of the pentad, transcribed from Section 3.
# Each entry is (coefficient, [five 1-based index pairs (i,j), i < j]).
PENTAD_PAPER = [
    (+1, [(1, 2), (1, 3), (2, 4), (3, 5), (4, 5)]),
    (-1, [(1, 2), (1, 3), (2, 5), (3, 4), (4, 5)]),
    (-1, [(1, 2), (1, 4), (2, 3), (3, 5), (4, 5)]),
    (+1, [(1, 2), (1, 4), (2, 5), (3, 4), (3, 5)]),
    (+1, [(1, 2), (1, 5), (2, 3), (3, 4), (4, 5)]),
    (-1, [(1, 2), (1, 5), (2, 4), (3, 4), (3, 5)]),
    (+1, [(1, 3), (1, 4), (2, 3), (2, 5), (4, 5)]),
    (-1, [(1, 3), (1, 4), (2, 4), (2, 5), (3, 5)]),
    (-1, [(1, 3), (1, 5), (2, 3), (2, 4), (4, 5)]),
    (+1, [(1, 3), (1, 5), (2, 4), (2, 5), (3, 4)]),
    (+1, [(1, 4), (1, 5), (2, 3), (2, 4), (3, 5)]),
    (-1, [(1, 4), (1, 5), (2, 3), (2, 5), (3, 4)]),
]

PULLBACK_CONST_PAPER = 12                        # the 12 in 12 t^10 Pi^2

# I7: Lemma 3.  Rows are m_1..m_4; columns are (F_2, F_3, F_4, F_5).
LEMMA3_MONOMIALS = [
    ("A_134 A_234 A_15 A_25", [(1, 3, 4), (2, 3, 4), (1, 5), (2, 5)]),
    ("A_124 A_234 A_15 A_35", [(1, 2, 4), (2, 3, 4), (1, 5), (3, 5)]),
    ("A_124 A_134 A_25 A_35", [(1, 2, 4), (1, 3, 4), (2, 5), (3, 5)]),
    ("A_123 A_234 A_15 A_45", [(1, 2, 3), (2, 3, 4), (1, 5), (4, 5)]),
]
LEMMA3_MATRIX_PAPER = [[0, 4, 4, 4],
                       [4, 0, 4, 4],
                       [-4, -4, 0, 0],
                       [4, 4, 0, 4]]
LEMMA3_DET_PAPER = -512

DIM_HD5_PAPER = 250                              # I8 consequence printed in (4)
DIM_P54_BOUND_PAPER = 254                        # Theorem 1 / Proposition 4

# I8: the two Schur functors occurring in Oeding's decomposition.
LAMBDA_IN_T = (2, 2)
LAMBDA_OUT_T = (4,)

# ---------------------------------------------------------------------
# Coordinate bookkeeping for the 32 coordinates A_alpha (n = 5).
# ---------------------------------------------------------------------
ALPHAS = list(itertools.product((0, 1), repeat=N))
AIDX = dict((a, k) for k, a in enumerate(ALPHAS))
IDX2ALPHA = dict((k, a) for a, k in AIDX.items())


def alpha_of_set(S):
    """Indicator tuple of a 1-based index set S <= [5]."""
    a = [0] * N
    for k in S:
        a[k - 1] = 1
    return tuple(a)


def set_of_alpha(a):
    return tuple(k + 1 for k in range(len(a)) if a[k] == 1)


def avar(S):
    return AIDX[alpha_of_set(S)]


# ======================================================================
# The exhibited object, rebuilt from the paper's DEFINITION (I2, I3).
# A-space monomials are sorted tuples of indicator tuples alpha.
# ======================================================================

def amul(P, Q):
    R = {}
    for m1, c1 in P.items():
        for m2, c2 in Q.items():
            m = tuple(sorted(m1 + m2))
            R[m] = R.get(m, 0) + c1 * c2
    return dict((m, c) for m, c in R.items() if c)


def asub(P, Q):
    R = dict(P)
    for m, c in Q.items():
        v = R.get(m, 0) - c
        if v:
            R[m] = v
        elif m in R:
            del R[m]
    return R


def q_entry(i, a, b, n=N, eps=None):
    """q^{(i)}_{ab} as a quadratic form in the coordinates A_alpha, alpha in
    {0,1}^n.  Exactly formula (1) of the paper: only pairs with
    eps_{alpha_j beta_j} != 0 for all j != i survive, i.e. beta_j = 1 - alpha_j."""
    if eps is None:
        eps = EPS
    others = [j for j in range(n) if j != i]
    R = {}
    for al_all in itertools.product((0, 1), repeat=n):
        if al_all[i] != a:
            continue
        for be_all in itertools.product((0, 1), repeat=n):
            if be_all[i] != b:
                continue
            sg = 1
            for j in others:
                sg *= eps[(al_all[j], be_all[j])]
                if sg == 0:
                    break
            if sg == 0:
                continue
            key = tuple(sorted((al_all, be_all)))
            R[key] = R.get(key, 0) + sg
    return dict((m, c) for m, c in R.items() if c)


def build_Q_and_B(i, n=N, eps=None):
    """Returns (Q, B) with Q a dict (a,b) -> quadratic form and B = det Q."""
    Q = {}
    for a in (0, 1):
        for b in (0, 1):
            Q[(a, b)] = q_entry(i, a, b, n, eps)
    B = asub(amul(Q[(0, 0)], Q[(1, 1)]), amul(Q[(0, 1)], Q[(1, 0)]))
    return Q, B


def subset_weight(mon):
    """Total subset weight sum |I_k| of an A-space monomial."""
    return sum(sum(al) for al in mon)


# ======================================================================
# The symmetric-matrix side: Z[x_ij (i <= j), t] and the principal minors.
# ======================================================================

class SymRing(object):
    """Polynomial ring Z[x_ij : 1 <= i <= j <= n][t] with packed monomials."""

    def __init__(self, n):
        self.n = n
        self.pairs = [(i, j) for i in range(n) for j in range(i, n)]
        self.xidx = dict((p, k) for k, p in enumerate(self.pairs))
        self.tvar = len(self.pairs)
        self.nv = self.tvar + 1
        self.minor = {}
        for r in range(n + 1):
            for I in itertools.combinations(range(n), r):
                self.minor[I] = self._detpoly(I)

    def x(self, i, j):
        """Monomial key of x_ij for 0-based i, j (symmetric)."""
        return 1 << (SH * self.xidx[(min(i, j), max(i, j))])

    def _detpoly(self, I):
        """det X_I as an exact polynomial, Leibniz expansion (|I| <= 5)."""
        k = len(I)
        if k == 0:
            return {0: 1}
        R = {}
        for p in itertools.permutations(range(k)):
            sg = permsign(p)
            m = 0
            for a in range(k):
                m += self.x(I[a], I[p[a]])
            R[m] = R.get(m, 0) + sg
        return dict((mm, c) for mm, c in R.items() if c)

    def t_power(self, e):
        return {pk([(self.tvar, e)]): 1}


def pentad_from_paper(R5):
    """The pentad as printed in the paper (input I5), as a polynomial."""
    P = {}
    for c, edges in PENTAD_PAPER:
        m = 0
        for (i, j) in edges:
            m += R5.x(i - 1, j - 1)
        P[m] = P.get(m, 0) + c
    return dict((m, c) for m, c in P.items() if c)


def hamiltonian_cycles_K5():
    """The distinct Hamiltonian cycles of K_5, each as (frozenset of edges,
    list of all 10 cyclic words).  Derived, nothing taken from the paper."""
    cyc = {}
    for p in itertools.permutations(range(1, 6)):
        edges = frozenset(frozenset((p[k], p[(k + 1) % 5])) for k in range(5))
        cyc.setdefault(edges, []).append(p)
    return cyc


def pentad_derived(R5):
    """Independent reconstruction of the pentad: sum over the 12 Hamiltonian
    cycles C of K_5 of sgn(word(C)) * prod_{e in C} x_e.  The sign is an
    invariant of C: rotating a 5-cycle word is an even permutation and so is
    reversing positions 2..5, so every one of the 10 words has the same parity.
    Returns (polynomial, well_defined_flag, number_of_cycles)."""
    cyc = hamiltonian_cycles_K5()
    P = {}
    welldef = True
    for edges, words in cyc.items():
        signs = set(permsign(w) for w in words)
        if len(signs) != 1:
            welldef = False
            continue
        s = signs.pop()
        m = 0
        for e in edges:
            i, j = sorted(e)
            m += R5.x(i - 1, j - 1)
        P[m] = P.get(m, 0) + s
    return dict((m, c) for m, c in P.items() if c), welldef, len(cyc)


# ======================================================================
# The pullback along phi_n:  A_alpha  |-->  t^{n-|alpha|} det X_{I(alpha)}
# ======================================================================

def phi_images(R):
    """dict alpha -> polynomial t^{n-|alpha|} det X_I, for all alpha in {0,1}^n."""
    n = R.n
    out = {}
    for al in itertools.product((0, 1), repeat=n):
        I = tuple(k for k in range(n) if al[k] == 1)
        out[al] = pshift(R.minor[I], R.tvar, n - len(I))
    return out


def subst_A_poly(R, P, images):
    """Apply the substitution `images` to an A-space polynomial P."""
    out = {}
    for mon, c in P.items():
        acc = {0: c}
        for al in mon:
            acc = pmul(acc, images[al])
        out = padd(out, acc)
    return out


def rhs_pullback(R, pentad, const, texp):
    """const * t^texp * pentad^2."""
    return pshift(pscale(pmul(pentad, pentad), const), R.tvar, texp)


# ======================================================================
# The GL_2^n action on the coordinates A, and Cayley's hyperdeterminant.
# ======================================================================

def act_on_A(g, Aval, n):
    """(g.A)_alpha = sum_beta prod_j (g_j)_{alpha_j beta_j} A_beta, on integer
    values.  g is a list of n integer 2x2 matrices."""
    betas = list(itertools.product((0, 1), repeat=n))
    out = {}
    for al in betas:
        s = 0
        for be in betas:
            p = 1
            for j in range(n):
                p *= g[j][al[j]][be[j]]
                if p == 0:
                    break
            if p:
                s += p * Aval[be]
        out[al] = s
    return out


def eval_A_poly(P, Aval):
    """Evaluate an A-space polynomial at integer values Aval."""
    tot = 0
    for mon, c in P.items():
        pr = c
        for al in mon:
            pr *= Aval[al]
        tot += pr
    return tot


def eval_q(Aval, i, n, eps=None):
    """Q_i and B_i = det Q_i at integer values, straight from formula (1).
    Independent of the polynomial machinery: an evaluation code path."""
    if eps is None:
        eps = EPS
    others = [j for j in range(n) if j != i]
    q = {}
    for a in (0, 1):
        for b in (0, 1):
            s = 0
            for rest in itertools.product((0, 1), repeat=n - 1):
                al = [0] * n
                be = [0] * n
                al[i] = a
                be[i] = b
                sg = 1
                for k, j in enumerate(others):
                    al[j] = rest[k]
                    be[j] = 1 - rest[k]
                    sg *= eps[(al[j], be[j])]
                if sg:
                    s += sg * Aval[tuple(al)] * Aval[tuple(be)]
            q[(a, b)] = s
    return q, q[(0, 0)] * q[(1, 1)] - q[(0, 1)] * q[(1, 0)]


# Input I10: Cayley's 2x2x2 hyperdeterminant, standard 12-term expansion.
# This table is the SINGLE source of truth for it; everything else evaluates it.
#   Det(a) = sum_{4 terms} a_S^2 a_{S'}^2
#          - 2 sum_{6 terms} a a a a  +  4 sum_{2 terms} a a a a
CAYLEY_TERMS = (
    (+1, ('000', '000', '111', '111')),
    (+1, ('001', '001', '110', '110')),
    (+1, ('010', '010', '101', '101')),
    (+1, ('100', '100', '011', '011')),
    (-2, ('000', '001', '110', '111')),
    (-2, ('000', '010', '101', '111')),
    (-2, ('000', '011', '100', '111')),
    (-2, ('001', '010', '101', '110')),
    (-2, ('001', '011', '100', '110')),
    (-2, ('010', '011', '100', '101')),
    (+4, ('000', '011', '101', '110')),
    (+4, ('001', '010', '100', '111')),
)


def cayley_A_poly():
    """Cayley's hyperdeterminant as an A-space polynomial in the eight
    coordinates A_alpha, alpha in {0,1}^3."""
    P = {}
    for c, strs in CAYLEY_TERMS:
        mon = tuple(sorted(tuple(int(ch) for ch in s) for s in strs))
        P[mon] = P.get(mon, 0) + c
    return dict((m, c) for m, c in P.items() if c)


def cayley_2x2x2(a):
    """Cayley's hyperdeterminant of a 2x2x2 integer array a[(i,j,k)]."""
    return eval_A_poly(cayley_A_poly(), a)


# ======================================================================
# Representation theory of SL_2, computed from scratch.
# ======================================================================

def ssyt_list(lam, k=2):
    """All semistandard Young tableaux of shape lam with entries in 1..k,
    returned as content lists.  Rows weakly increase, columns strictly."""
    lam = [l for l in lam if l > 0]
    rows = len(lam)
    out = []

    def rec(r, prev, acc):
        if r == rows:
            out.append(list(acc))
            return
        row = []

        def fill(c, last):
            if c == lam[r]:
                rec(r + 1, list(row), acc + row)
                return
            lo = last
            if r > 0 and c < len(prev):
                lo = max(lo, prev[c] + 1)
            for val in range(lo, k + 1):
                row.append(val)
                fill(c + 1, val)
                row.pop()

        fill(0, 1)

    rec(0, [], [])
    return out


def schur_dim(lam, k=2):
    """dim S_lam(C^k) = number of SSYT of shape lam with entries in 1..k."""
    return len(ssyt_list(lam, k))


def schur_character_2(lam):
    """Character of S_lam(C^2) as a dict (a,b) -> multiplicity of z_1^a z_2^b."""
    d = {}
    for T in ssyt_list(lam, 2):
        a = sum(1 for val in T if val == 1)
        b = len(T) - a
        d[(a, b)] = d.get((a, b), 0) + 1
    return d


def sym_character_2(m, shift_by=0):
    """Character of Sym^m(C^2) tensored with det^{shift_by}."""
    return dict(((k + shift_by, m - k + shift_by), 1) for k in range(m + 1))


def sym_power_matrix(g, m):
    """Matrix of Sym^m(g) on the basis e0^{m-r} e1^r, r = 0..m."""
    M = [[0] * (m + 1) for _ in range(m + 1)]
    for r in range(m + 1):
        coeffs = {0: 1}
        for _ in range(m - r):                    # multiply by g(e0)
            nxt = {}
            for kk, c in coeffs.items():
                nxt[kk] = nxt.get(kk, 0) + c * g[0][0]
                nxt[kk + 1] = nxt.get(kk + 1, 0) + c * g[1][0]
            coeffs = nxt
        for _ in range(r):                        # multiply by g(e1)
            nxt = {}
            for kk, c in coeffs.items():
                nxt[kk] = nxt.get(kk, 0) + c * g[0][1]
                nxt[kk + 1] = nxt.get(kk + 1, 0) + c * g[1][1]
            coeffs = nxt
        for kk, c in coeffs.items():
            M[kk][r] += c
    return M


def sym_invariant_dim(m):
    """dim Sym^m(C^2)^{SL_2}, computed as the joint fixed space of the two
    unipotents u = [[1,1],[0,1]] and l = [[1,0],[1,1]], which generate SL_2."""
    rows = []
    for g in ([[1, 1], [0, 1]], [[1, 0], [1, 1]]):
        M = sym_power_matrix(g, m)
        for a in range(m + 1):
            rows.append([M[a][b] - (1 if a == b else 0) for b in range(m + 1)])
    return (m + 1) - rank_exact(rows)


# ======================================================================
# GROUP 1 -- the exhibited object is well formed and is what the paper says
# ======================================================================

def checks_object(ctx):
    section("GROUP 1: the exhibited object, decoded and printed back")

    # 1. eps is the SL_2-invariant symplectic form:  g^T eps g = (det g) eps.
    E = [[EPS[(0, 0)], EPS[(0, 1)]], [EPS[(1, 0)], EPS[(1, 1)]]]
    rng = random.Random(20260824)
    bad = []
    for _ in range(200):
        g = [[rng.randint(-6, 6) for _ in range(2)] for _ in range(2)]
        dg = g[0][0] * g[1][1] - g[0][1] * g[1][0]
        lhs = [[sum(g[c][a] * E[c][d] * g[d][b] for c in (0, 1) for d in (0, 1))
                for b in (0, 1)] for a in (0, 1)]
        if any(lhs[a][b] != dg * E[a][b] for a in (0, 1) for b in (0, 1)):
            bad.append(g)
    ok = (not bad) and E[0][1] == 1 and E[1][0] == -1 and E[0][0] == 0 and E[1][1] == 0
    record("eps_is_SL2_invariant_alternating_form", ok,
           "g^T eps g = (det g) eps on 200 integer g; failures=%d" % len(bad))

    # 2. 32 coordinates, indicator <-> subset bijection.
    subsets = [tuple(sorted(s)) for r in range(N + 1)
               for s in itertools.combinations(range(1, N + 1), r)]
    round_trip = all(set_of_alpha(alpha_of_set(s)) == s for s in subsets)
    ok = (len(ALPHAS) == 32 and len(set(subsets)) == 32 and round_trip
          and len(set(AIDX.values())) == 32)
    record("coordinate_set_is_32_subsets_of_[5]", ok,
           "|{0,1}^5|=%d, |2^[5]|=%d, bijection=%s"
           % (len(ALPHAS), len(set(subsets)), round_trip))

    # 3. Q_i symmetric (the paper's "the four contractions make Q_i symmetric").
    sym = [ctx['Q'][i][(0, 1)] == ctx['Q'][i][(1, 0)] for i in range(N)]
    nonzero = [len(ctx['Q'][i][(0, 1)]) > 0 for i in range(N)]
    record("Q_i_is_symmetric_for_all_i", all(sym) and all(nonzero),
           "q01 == q10 for i=1..5: %s; q01 has %d monomials"
           % (all(sym), len(ctx['Q'][0][(0, 1)])))

    # 4. B_i nonzero quartics, 192 monomials, pairwise distinct.
    Bs = ctx['B']
    degs = set(len(m) for B in Bs for m in B)
    counts = [len(B) for B in Bs]
    distinct = all(Bs[i] != Bs[j] for i in range(N) for j in range(i + 1, N))
    ok = degs == set([4]) and all(c > 0 for c in counts) and distinct
    record("B_i_are_distinct_nonzero_quartics", ok,
           "monomial counts %s, all of degree 4: %s, pairwise distinct: %s"
           % (counts, degs == set([4]), distinct))

    # 5. Subset-weight grading (the proof's reduction to the chart t = 1).
    wq = True
    for i in range(N):
        for (a, b) in ((0, 0), (0, 1), (1, 0), (1, 1)):
            ws = set(subset_weight(m) for m in ctx['Q'][i][(a, b)])
            if ws != set([4 + a + b]):
                wq = False
    wb = set(subset_weight(m) for B in Bs for m in B)
    record("subset_weight_grading_q=4+a+b_and_B=10", wq and wb == set([10]),
           "q weights match 4+a+b: %s; B weights = %s" % (wq, sorted(wb)))

    # 6. The pentad as transcribed: 12 distinct squarefree quintic monomials,
    #    coefficients +-1 (six of each), off-diagonal variables only.
    PIp = ctx['PIp']
    R5 = ctx['R5']
    exps = [unpk(m, R5.nv) for m in PIp]
    sqfree = all(max(e) <= 1 for e in exps)
    quintic = all(sum(e) == 5 for e in exps)
    offdiag = all(e[R5.xidx[(i, i)]] == 0 for e in exps for i in range(N))
    signs = sorted(PIp.values())
    ok = (len(PIp) == 12 and sqfree and quintic and offdiag
          and signs.count(1) == 6 and signs.count(-1) == 6)
    record("pentad_transcription_12_squarefree_quintics", ok,
           "terms=%d squarefree=%s degree5=%s offdiagonal=%s (+1)x%d (-1)x%d"
           % (len(PIp), sqfree, quintic, offdiag, signs.count(1),
              signs.count(-1)))

    # 7. Independent reconstruction of the pentad from K_5.
    PId, welldef, ncyc = ctx['PId'], ctx['PI_welldef'], ctx['PI_ncyc']
    record("pentad_equals_signed_hamiltonian_cycle_sum_of_K5",
           welldef and ncyc == 12 and PId == PIp,
           "12 Ham cycles of K_5 = %d, sign well defined on all 10 words each: "
           "%s, derived pentad == printed pentad: %s"
           % (ncyc, welldef, PId == PIp))

    # 7b. Pi is alternating under relabelling of the five indices.  This is what
    #     makes Pi^2 -- and hence the whole statement -- independent of the
    #     ordering of a chosen 5-subset S, which the n > 5 argument relies on.
    bad = []
    for sig in itertools.permutations(range(N)):
        Q = relabel_poly(R5, PIp, sig)
        if Q != pscale(PIp, permsign(sig)):
            bad.append(sig)
    sq = all(relabel_poly(R5, pmul(PIp, PIp), sig) == pmul(PIp, PIp)
             for sig in itertools.permutations(range(N)))
    record("pentad_is_alternating_and_pentad_squared_is_symmetric",
           not bad and sq,
           "Pi(sigma X) = sgn(sigma) Pi(X) for all %d permutations "
           "(failures=%d); Pi^2 invariant under all of them: %s"
           % (math.factorial(N), len(bad), sq))


def relabel_poly(R, P, sig):
    """Relabel the matrix indices of a polynomial by the 0-based permutation
    sig: x_{ij} |--> x_{sig(i) sig(j)}.  Leaves t alone."""
    out = {}
    for m, c in P.items():
        e = unpk(m, R.nv)
        mm = e[R.tvar] << (SH * R.tvar)
        for k, (i, j) in enumerate(R.pairs):
            if e[k]:
                a, b = sig[i], sig[j]
                mm += e[k] << (SH * R.xidx[(min(a, b), max(a, b))])
        out[mm] = out.get(mm, 0) + c
    return dict((m, c) for m, c in out.items() if c)


def random_sym(rng, n, lo=-6, hi=6):
    X = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = rng.randint(lo, hi)
            X[i][j] = v
            X[j][i] = v
    return X


def minor_values(X, t, n, only=None):
    """A_alpha = t^{n-|I|} det X_I for alpha in {0,1}^n (all, or only those
    supported inside the index set `only`)."""
    out = {}
    for al in itertools.product((0, 1), repeat=n):
        I = [k for k in range(n) if al[k] == 1]
        if only is not None and any(k not in only for k in I):
            continue
        M = [[X[i][j] for j in I] for i in I]
        out[al] = (t ** (n - len(I))) * det_int(M)
    return out


# ======================================================================
# GROUP 2 -- every hypothesis of the statement being acted on
# ======================================================================

def checks_hypotheses(ctx):
    section("GROUP 2: hypotheses of Holtz-Sturmfels Conjecture 14 and of Oeding")
    R5 = ctx['R5']
    im = ctx['im']

    # 8. phi_5 is the map the paper says it is.
    empty = tuple([0] * N)
    ok_empty = im[empty] == {pk([(R5.tvar, N)]): 1}
    graded = True
    for al, P in im.items():
        w = sum(al)
        for m in P:
            e = unpk(m, R5.nv)
            if sum(e[:R5.tvar]) != w or e[R5.tvar] != N - w:
                graded = False
    distinct = len(set(tuple(sorted(P.items())) for P in im.values())) == 32
    sizes = [len(im[alpha_of_set(tuple(range(1, r + 1)))]) for r in range(N + 1)]
    record("phi5_coordinates_are_the_32_principal_minors",
           ok_empty and graded and distinct,
           "A_empty -> t^5: %s; bidegree (|I|, 5-|I|) on every coordinate: %s; "
           "32 distinct: %s; monomials of det X_[1..r] = %s"
           % (ok_empty, graded, distinct, sizes))

    # 9. Lemma 1: B_i is SL_2^5-invariant.  Verified in the sharper covariant
    #    form, which also makes the test discriminating.
    rng = random.Random(1729)
    n_sl2, n_gen, bad_Q, bad_B, discriminating = 0, 0, 0, 0, 0
    for trial in range(40):
        Aval = dict((al, rng.randint(-9, 9)) for al in ALPHAS)
        if trial % 2 == 0:                        # honest SL_2(Z) elements
            g = []
            for _ in range(N):
                a, b = rng.randint(-3, 3), rng.randint(-3, 3)
                g.append([[1 + a * b, a], [b, 1]])   # det = 1
            n_sl2 += 1
        else:
            g = [[[rng.randint(-3, 3) for _ in range(2)] for _ in range(2)]
                 for _ in range(N)]
            n_gen += 1
        dets = [gg[0][0] * gg[1][1] - gg[0][1] * gg[1][0] for gg in g]
        D = 1
        for d in dets:
            D *= d
        gA = act_on_A(g, Aval, N)
        for i in range(N):
            q0, b0 = eval_q(Aval, i, N)
            q1, b1 = eval_q(gA, i, N)
            pref = 1
            for j in range(N):
                if j != i:
                    pref *= dets[j]
            gi = g[i]
            M = [[sum(gi[a][c] * q0[(c, d)] * gi[b][d]
                      for c in (0, 1) for d in (0, 1)) for b in (0, 1)]
                 for a in (0, 1)]
            if any(q1[(a, b)] != pref * M[a][b] for a in (0, 1) for b in (0, 1)):
                bad_Q += 1
            if b1 != D * D * b0:
                bad_B += 1
            if D * D != 1 and b0 != 0 and b1 != b0:
                discriminating += 1
    ok = (bad_Q == 0 and bad_B == 0 and n_sl2 > 0 and discriminating > 0)
    record("B_i_is_SL2^5_invariant_covariance_law", ok,
           "Q_i(gA) = (prod_{j!=i} det g_j) g_i Q_i g_i^T and B_i(gA) = "
           "(prod det g_j)^2 B_i(A) on %d SL_2 and %d general GL_2 tuples "
           "(x5 values of i): Q failures=%d, B failures=%d; %d trials had "
           "(prod det)^2 != 1 with B_i(gA) != B_i(A), so the test discriminates"
           % (n_sl2, n_gen, bad_Q, bad_B, discriminating))

    # 10/11/12/13. Oeding's decomposition (input I8) evaluated from scratch.
    d22 = schur_dim(LAMBDA_IN_T)
    d4 = schur_dim(LAMBDA_OUT_T)
    record("schur_dims_of_C2_by_counting_SSYT", d22 == 1 and d4 == 5,
           "dim S_(2,2) C^2 = %d, dim S_(4) C^2 = %d" % (d22, d4))

    chi22 = schur_character_2(LAMBDA_IN_T) == sym_character_2(0, 2)
    chi4 = schur_character_2(LAMBDA_OUT_T) == sym_character_2(4, 0)
    record("schur_restricts_to_SL2_as_Sym^(lam1-lam2)", chi22 and chi4,
           "char S_(2,2) = (z1 z2)^2 * char Sym^0: %s; char S_(4) = "
           "char Sym^4: %s" % (chi22, chi4))

    inv = [sym_invariant_dim(m) for m in range(7)]
    ok = inv == [1, 0, 0, 0, 0, 0, 0]
    record("Sym^m(C2)^SL2_is_zero_for_m>0", ok,
           "dim Sym^m(C^2)^{SL_2} for m=0..6 = %s (joint kernel of the two "
           "unipotent generators, exact rational rank)" % (inv,))

    ctx['inv_in_T'] = sym_invariant_dim(LAMBDA_IN_T[0] - LAMBDA_IN_T[1])
    ctx['inv_out_T'] = sym_invariant_dim(LAMBDA_OUT_T[0])
    ctx['d22'], ctx['d4'] = d22, d4

    # dim HD(5): sum over the 10 summands of Oeding's decomposition.
    tot = 0
    for T in itertools.combinations(range(N), 3):
        p = 1
        for j in range(N):
            p *= d22 if j in T else d4
        tot += p
    closed = math.comb(N, 3) * 5 ** (N - 3)
    record("dim_HD5_equals_250_from_Oeding_decomposition",
           tot == DIM_HD5_PAPER and closed == DIM_HD5_PAPER,
           "sum over the %d summands = %d; binom(5,3)*5^2 = %d; paper says %d"
           % (math.comb(N, 3), tot, closed, DIM_HD5_PAPER))
    ctx['dim_HD5'] = tot

    # HD(5)^{G_5} = 0: the invariants of a tensor product of SL_2-modules is the
    # tensor product of the invariants, so each summand contributes the product
    # of the per-factor invariant dimensions.
    per_summand = []
    for T in itertools.combinations(range(N), 3):
        p = 1
        for j in range(N):
            p *= ctx['inv_in_T'] if j in T else ctx['inv_out_T']
        per_summand.append(p)
    record("HD5_has_no_nonzero_G5_invariant", sum(per_summand) == 0
           and len(per_summand) == 10,
           "invariant dimension of each of the %d summands = %s, total = %d"
           % (len(per_summand), sorted(set(per_summand)), sum(per_summand)))

    # 14. HD(n) <= (P_n)_4 at the level of the orbit generator: Cayley's
    #     2x2x2 hyperdeterminant on every 3-subset slice of phi_5(X,t), and on
    #     GL_2^5-twists of it (which is what spans the orbit span HD(5)).
    rng = random.Random(31337)
    slices_checked = 0
    nz_plain, nz_twist, nz_control = 0, 0, 0
    for _ in range(6):
        X = random_sym(rng, N)
        t = rng.randint(1, 5)
        Av = minor_values(X, t, N)
        g = [[[rng.randint(-3, 3) for _ in range(2)] for _ in range(2)]
             for _ in range(N)]
        gA = act_on_A(g, Av, N)
        for T in itertools.combinations(range(N), 3):
            rest = [j for j in range(N) if j not in T]
            for fix in itertools.product((0, 1), repeat=N - 3):
                sub, subg = {}, {}
                for b in itertools.product((0, 1), repeat=3):
                    al = [0] * N
                    for k, j in enumerate(T):
                        al[j] = b[k]
                    for k, j in enumerate(rest):
                        al[j] = fix[k]
                    sub[b] = Av[tuple(al)]
                    subg[b] = gA[tuple(al)]
                slices_checked += 1
                if cayley_2x2x2(sub) != 0:
                    nz_plain += 1
                if cayley_2x2x2(subg) != 0:
                    nz_twist += 1
    for _ in range(50):                            # control: Cayley is not 0
        arr = dict((b, rng.randint(-9, 9))
                   for b in itertools.product((0, 1), repeat=3))
        if cayley_2x2x2(arr) != 0:
            nz_control += 1
    # ... and, so that this does not rest on sampling, the SAME statement as an
    # exact coefficientwise identity: the pullback along phi_5 of Cayley's
    # hyperdeterminant of each of the 40 (3-subset, complement-pattern) slices
    # is the zero polynomial of Z[x_11..x_55, t].  This is the computational
    # content of "HD(5) <= (P_5)_4" (input I9) at the orbit generators.
    sym_slices, sym_nonzero, sym_terms = 0, 0, 0
    for T in itertools.combinations(range(N), 3):
        for fix in itertools.product((0, 1), repeat=N - 3):
            D = embed_cayley_slice(T, fix)
            sym_slices += 1
            sym_terms += len(D)
            if subst_A_poly(ctx['R5'], D, ctx['im']) != {}:
                sym_nonzero += 1
    record("cayley_hyperdet_vanishes_on_phi5_and_its_GL2^5_twists",
           nz_plain == 0 and nz_twist == 0 and nz_control > 0
           and sym_slices == 40 and sym_terms == 40 * len(cayley_A_poly())
           and sym_nonzero == 0,
           "%d (3-subset, complement-pattern) slices: nonvanishing on phi_5 = "
           "%d, on GL_2^5-twisted phi_5 = %d; EXACT symbolic pullback of all %d "
           "slices (%d A-monomials) along phi_5: %d nonzero residuals; control: "
           "Cayley nonzero on %d/50 generic arrays"
           % (slices_checked, nz_plain, nz_twist, sym_slices, sym_terms,
              sym_nonzero, nz_control))

    # 15. Identify the paper's construction with Cayley's hyperdeterminant at
    #     n = 3: det Q_i = -Det_{2x2x2} as polynomials in the eight A_alpha.
    cay = cayley_A_poly()
    negcay = dict((m, -c) for m, c in cay.items())
    same = []
    for i in range(3):
        _, B3 = build_Q_and_B(i, n=3)
        same.append(B3 == negcay)
    record("paper_construction_at_n=3_is_minus_cayley_hyperdet",
           all(same) and len(cay) == 12,
           "det Q_i == -Det_{2x2x2} coefficientwise for i=1,2,3: %s "
           "(Cayley has %d monomials)" % (same, len(cay)))


def phi_images_affine(R):
    """The chart t = 1: A_alpha |--> det X_I."""
    n = R.n
    out = {}
    for al in itertools.product((0, 1), repeat=n):
        I = tuple(k for k in range(n) if al[k] == 1)
        out[al] = dict(R.minor[I])
    return out


def eval_xt(P, R, xval, tval):
    """Evaluate a packed polynomial at integer x_ij = xval[(i,j)], t = tval."""
    tot = 0
    for m, c in P.items():
        e = unpk(m, R.nv)
        pr = c
        for k, (i, j) in enumerate(R.pairs):
            if e[k]:
                pr *= xval[i][j] ** e[k]
        if e[R.tvar]:
            pr *= tval ** e[R.tvar]
        tot += pr
    return tot


# ======================================================================
# GROUP 3 -- the load-bearing computation: the conclusion is violated
# ======================================================================

def checks_loadbearing(ctx):
    section("GROUP 3: the load-bearing identity and the refutation")
    R5, im, PIp = ctx['R5'], ctx['im'], ctx['PIp']

    # 16. Proposition 2, in the fully homogeneous form (t carried along).
    rhs = rhs_pullback(R5, PIp, PULLBACK_CONST_PAPER, 10)
    pulls = []
    for i in range(N):
        pulls.append(subst_A_poly(R5, ctx['B'][i], im))
    eq = [p == rhs for p in pulls]
    residuals = [len(padd(p, rhs, -1)) for p in pulls]
    record("pullback_B_i_equals_12_t^10_pentad^2_for_all_i",
           all(eq) and len(rhs) > 0,
           "coefficientwise in Z[x_11..x_55, t] (16 variables): matches for "
           "i=1..5: %s; |B_i o phi_5| = %s monomials, |12 t^10 Pi^2| = %d; "
           "nonzero residual monomials %s"
           % (all(eq), [len(p) for p in pulls], len(rhs), residuals))
    ctx['pulls'] = pulls

    # 17. The same identity on the affine chart t = 1, dehomogenised, computed
    #     through a different substitution.
    ima = phi_images_affine(R5)
    rhs_a = pscale(pmul(PIp, PIp), PULLBACK_CONST_PAPER)
    eq_a = [subst_A_poly(R5, ctx['B'][i], ima) == rhs_a for i in range(N)]
    record("pullback_identity_on_the_affine_chart_t=1", all(eq_a),
           "B_i((det X_I)_I) == 12 Pi(X)^2 coefficientwise for i=1..5: %s "
           "(%d monomials)" % (all(eq_a), len(rhs_a)))

    # 18. Independent code path: exact integer evaluation, no polynomial
    #     expansion at all (formula (1) evaluated directly).
    rng = random.Random(90210)
    bad, npts, nz = 0, 0, 0
    for _ in range(30):
        X = random_sym(rng, N, -8, 8)
        t = rng.randint(-4, 4)
        Av = minor_values(X, t, N)
        target = PULLBACK_CONST_PAPER * (t ** 10) * eval_xt(PIp, R5, X, t) ** 2
        if target != 0:
            nz += 1
        for i in range(N):
            _, b = eval_q(Av, i, N)
            npts += 1
            if b != target:
                bad += 1
    record("pullback_identity_at_exact_integer_points", bad == 0 and nz >= 20,
           "%d exact evaluations (30 integer (X,t), 5 values of i) via formula "
           "(1) directly: mismatches=%d; %d/30 points had 12 t^10 Pi^2 != 0"
           % (npts, bad, nz))

    # 19. F_i = B_1 - B_i pulls back to the zero polynomial: F_i in (P_5)_4.
    Fs = [asub(ctx['B'][0], ctx['B'][i]) for i in range(1, N)]
    ctx['F'] = Fs
    zero = [subst_A_poly(R5, F, im) == {} for F in Fs]
    zero_a = [subst_A_poly(R5, F, ima) == {} for F in Fs]
    ctx['F_in_P5'] = all(zero) and all(zero_a)
    record("F_i_pull_back_to_zero_hence_F_i_in_P5_degree4",
           all(zero) and all(zero_a),
           "F_i o phi_5 == 0 in Z[x,t] for i=2,3,4,5: %s; and on the chart "
           "t=1: %s" % (all(zero), all(zero_a)))

    # 20. F_i are nonzero homogeneous quartics.
    counts = [len(F) for F in Fs]
    degs = set(len(m) for F in Fs for m in F)
    wts = set(subset_weight(m) for F in Fs for m in F)
    record("F_i_are_nonzero_homogeneous_quartics",
           all(c > 0 for c in counts) and degs == set([4]) and wts == set([10]),
           "monomial counts %s, degrees %s, subset weights %s"
           % (counts, sorted(degs), sorted(wts)))

    # 21. Lemma 3: recompute the printed 4 x 4 coefficient matrix and its det.
    mat = []
    for name, sets in LEMMA3_MONOMIALS:
        mon = tuple(sorted(alpha_of_set(s) for s in sets))
        mat.append([Fs[c].get(mon, 0) for c in range(4)])
    d = det_int(mat)
    agree = mat == LEMMA3_MATRIX_PAPER
    print("    recomputed Lemma 3 matrix (rows m_1..m_4, cols F_2..F_5):")
    for name, row in zip([nm for nm, _ in LEMMA3_MONOMIALS], mat):
        print("      %-22s %s" % (name, row))
    record("lemma3_coefficient_matrix_and_determinant_recomputed",
           agree and d == LEMMA3_DET_PAPER,
           "matrix equals the one printed in Lemma 3: %s; recomputed det = %d, "
           "paper says %d" % (agree, d, LEMMA3_DET_PAPER))

    # 22. Rank over Q of the FULL coefficient matrix of (F_2,...,F_5): this does
    #     not depend on the four monomials the paper happened to choose.
    mons = sorted(set(m for F in Fs for m in F))
    rows = [[F.get(m, 0) for m in mons] for F in Fs]
    rk = rank_exact(rows)
    record("F_2..F_5_span_a_4_dimensional_space", rk == 4,
           "exact rank over Q of the 4 x %d integer coefficient matrix = %d"
           % (len(mons), rk))
    ctx['dim_W'] = rk

    # 23. W is G_5-trivial, hence meets HD(5) only in 0.
    rng = random.Random(6060842)
    bad_inv, ntr = 0, 0
    for _ in range(24):
        Aval = dict((al, rng.randint(-9, 9)) for al in ALPHAS)
        g = []
        for _ in range(N):
            a, b = rng.randint(-3, 3), rng.randint(-3, 3)
            g.append([[1 + a * b, a], [b, 1]])          # det = 1
        gA = act_on_A(g, Aval, N)
        coefs = [rng.randint(-5, 5) for _ in range(4)]
        w0 = sum(c * eval_A_poly(F, Aval) for c, F in zip(coefs, Fs))
        w1 = sum(c * eval_A_poly(F, gA) for c, F in zip(coefs, Fs))
        ntr += 1
        if w0 != w1:
            bad_inv += 1
    ok = (bad_inv == 0 and ctx['dim_W'] == 4 and ctx['dim_HD5'] == DIM_HD5_PAPER
          and ctx['inv_out_T'] == 0)
    record("W_is_G5_trivial_and_W_cap_HD5_is_zero", ok,
           "every element of W = span(F_2..F_5) is SL_2^5-invariant "
           "(%d random elements x random SL_2(Z) tuples, failures=%d); "
           "dim W = %d; HD(5)^{G_5} = 0, so W cap HD(5) = 0"
           % (ntr, bad_inv, ctx['dim_W']))

    # 24. The conclusion of Conjecture 14 is violated.  The criterion names the
    #     two facts the bound actually rests on -- W <= (P_5)_4 (the pullback
    #     computation of check 19) and dim W > 0 -- and not merely the arithmetic
    #     250 + 4 = 254, which would hold even if the membership had failed.
    lower = ctx['dim_HD5'] + ctx['dim_W']
    ok = (ctx['F_in_P5'] and ctx['dim_W'] > 0
          and lower == DIM_P54_BOUND_PAPER and lower > ctx['dim_HD5'])
    record("dim_P5_4_at_least_254_strictly_exceeds_dim_HD5_so_J5_neq_P5", ok,
           "F_i in (P_5)_4 (computed, check 19): %s; dim(P_5)_4 >= dim HD(5) + "
           "dim W = %d + %d = %d > %d = dim HD(5) = dim (J_5)_4  ==>  J_5 is "
           "strictly contained in P_5"
           % (ctx['F_in_P5'], ctx['dim_HD5'], ctx['dim_W'], lower,
              ctx['dim_HD5']))


# ======================================================================
# GROUP 3b -- HD(5) is not merely cited: its weight-zero subspace, the only
# part of it that could ever meet the G_5-trivial space W, is constructed.
#
# The torus diag(z_j, z_j^{-1}) in factor j acts on A_alpha by z_j^{1-2 alpha_j},
# so a quartic monomial A_{a1}..A_{a4} has j-weight 4 - 2 sum_k (a_k)_j.  The
# Lie algebra sl_2 in factor j acts on R_4 by the derivations
#    lie_lo_j : A_alpha -> A_{alpha with bit j cleared}  (if alpha_j = 1)
#    lie_hi_j : A_alpha -> A_{alpha with bit j set}      (if alpha_j = 0)
# (these are d/ds of the two one-parameter transvection subgroups at s = 0).
# ======================================================================

def weight_vec(mon):
    return tuple(4 - 2 * sum(al[j] for al in mon) for j in range(N))


def lie_op(P, j, hi):
    """Apply the derivation lie_hi_j (hi=True) or lie_lo_j (hi=False)."""
    out = {}
    want = 0 if hi else 1
    new = 1 if hi else 0
    for mon, c in P.items():
        for k in range(len(mon)):
            al = mon[k]
            if al[j] == want:
                nb = list(al)
                nb[j] = new
                nm = tuple(sorted(mon[:k] + (tuple(nb),) + mon[k + 1:]))
                out[nm] = out.get(nm, 0) + c
    return dict((m, c) for m, c in out.items() if c)


def transvection(P, j, c, hi):
    """The algebra automorphism rho(I + c X_j) = exp(c * lie op), exact."""
    out = {}
    want = 0 if hi else 1
    new = 1 if hi else 0
    for mon, co in P.items():
        acc = {(): co}
        for al in mon:
            opts = [(al, 1)]
            if al[j] == want:
                nb = list(al)
                nb[j] = new
                opts.append((tuple(nb), c))
            nxt = {}
            for key, cc in acc.items():
                for nal, mult in opts:
                    k2 = key + (nal,)
                    nxt[k2] = nxt.get(k2, 0) + cc * mult
            acc = nxt
        for key, cc in acc.items():
            k = tuple(sorted(key))
            out[k] = out.get(k, 0) + cc
    return dict((m, c) for m, c in out.items() if c)


def embed_cayley_slice(T, fix):
    """Cayley's hyperdeterminant in the coordinates A_alpha with alpha|_T free
    and alpha|_{[5]-T} = fix.  These 40 forms all lie in HD(5)."""
    rest = [j for j in range(N) if j not in T]
    out = {}
    for mon, c in cayley_A_poly().items():
        nm = []
        for b in mon:
            al = [0] * N
            for k, j in enumerate(T):
                al[j] = b[k]
            for k, j in enumerate(rest):
                al[j] = fix[k]
            nm.append(tuple(al))
        key = tuple(sorted(nm))
        out[key] = out.get(key, 0) + c
    return dict((m, c) for m, c in out.items() if c)


def project_weight0(P):
    z = (0,) * N
    return dict((m, c) for m, c in P.items() if weight_vec(m) == z)


def coeff_rows(polys):
    mons = sorted(set(m for P in polys for m in P))
    return [[P.get(m, 0) for m in mons] for P in polys], len(mons)


def checks_hd5_weight_zero(ctx):
    section("GROUP 3b: HD(5) constructed in weight zero, and W tested against it")
    Fs = ctx['F']

    # 25. Infinitesimal invariance of the F_i: annihilated by all ten
    #     generators of sl_2^5, and of torus weight zero.  Independent of the
    #     finite-group covariance test in Group 2.
    z = (0,) * N
    wt_ok = all(weight_vec(m) == z for F in Fs for m in F)
    killed = all(lie_op(F, j, hi) == {} for F in Fs for j in range(N)
                 for hi in (True, False))
    # Control: EACH of the ten operators is separately nonzero -- Cayley's
    # hyperdeterminant slices, which are NOT SL_2^5-invariant, do move.  An
    # aggregate "some operator moved something" control would be satisfied by
    # a single surviving operator while the other nine silently returned {},
    # which is exactly the failure mode that would make `killed` vacuous, so
    # the control is per operator.
    ctrl_tot, ctrl_nz = 0, 0
    ctrl_ops = set()
    for T in itertools.combinations(range(N), 3):
        for fix in itertools.product((0, 1), repeat=2):
            D = embed_cayley_slice(T, fix)
            for j in range(N):
                for hi in (True, False):
                    ctrl_tot += 1
                    if lie_op(D, j, hi) != {}:
                        ctrl_nz += 1
                        ctrl_ops.add((j, hi))
    record("F_i_annihilated_by_sl2^5_and_of_torus_weight_zero",
           wt_ok and killed and len(ctrl_ops) == 2 * N,
           "all monomials of F_2..F_5 have torus weight (0,0,0,0,0): %s; "
           "lie_lo_j F = lie_hi_j F = 0 for j=1..5: %s; control: the (not "
           "invariant) hyperdeterminant slices have nonzero image in %d of %d "
           "(slice, j, op) cases, and every one of the %d operators "
           "lie_lo_j/lie_hi_j is nonzero on at least one slice (%d/%d), so no "
           "operator is silently the zero map"
           % (wt_ok, killed, ctrl_nz, ctrl_tot, 2 * N, len(ctrl_ops), 2 * N))

    # 26. Predicted dimension of the weight-zero part of HD(5), from the
    #     multiplicity of the zero SL_2-weight in each Schur factor.
    def zero_weight_mult(lam):
        return sum(mult for (a, b), mult in schur_character_2(lam).items()
                   if a == b)
    m22, m4 = zero_weight_mult(LAMBDA_IN_T), zero_weight_mult(LAMBDA_OUT_T)
    pred = math.comb(N, 3) * (m22 ** 3) * (m4 ** 2)
    lam_str = lambda lam: "(" + ",".join(str(p) for p in lam) + ")"
    record("predicted_dim_of_weight_zero_part_of_HD5_is_10",
           m22 == 1 and m4 == 1 and pred == 10,
           "zero-weight multiplicity: S_%s -> %d, S_%s -> %d; so "
           "dim HD(5)^{torus} = binom(5,3) * %d^3 * %d^2 = %d"
           % (lam_str(LAMBDA_IN_T), m22, lam_str(LAMBDA_OUT_T), m4,
              m22, m4, pred))
    ctx['pred_w0'] = pred

    # 27. Construct weight-zero elements of HD(5) explicitly: push each of the
    #     40 hyperdeterminant slices to weight zero with the Lie operators, and
    #     add elements obtained through transvection words.  Their rank must be
    #     exactly the predicted 10.
    base = []
    for T in itertools.combinations(range(N), 3):
        rest = [j for j in range(N) if j not in T]
        for fix in itertools.product((0, 1), repeat=2):
            P = embed_cayley_slice(T, fix)
            for k, j in enumerate(rest):
                for _ in range(2):
                    P = lie_op(P, j, hi=(fix[k] == 0))
            if P:
                base.append(P)
    rng = random.Random(24680)
    Tlist = list(itertools.combinations(range(N), 3))
    extra = []
    for _ in range(150):
        T = Tlist[rng.randrange(len(Tlist))]
        rest = [j for j in range(N) if j not in T]
        fix = (rng.randrange(2), rng.randrange(2))
        P = embed_cayley_slice(T, fix)
        for k, j in enumerate(rest):
            for _ in range(2):
                P = transvection(P, j, rng.choice([1, -1, 2]), hi=(fix[k] == 0))
        for _ in range(rng.randint(0, 2)):
            j = T[rng.randrange(3)]
            P = transvection(P, j, rng.choice([1, -1]), hi=False)
            P = transvection(P, j, rng.choice([1, -1]), hi=True)
        Z = project_weight0(P)
        if Z:
            extra.append(Z)
    allP = base + extra
    wz = all(weight_vec(m) == z for P in allP for m in P)
    rows, nmon = coeff_rows(allP)
    rk = rank_exact(rows)
    record("computed_weight_zero_part_of_HD5_has_rank_10",
           rk == ctx['pred_w0'] and len(base) == 40 and wz,
           "%d explicit elements of HD(5) (the %d hyperdeterminant slices "
           "pushed to weight zero, plus %d transvection-word elements), all of "
           "weight zero: %s; exact rank over Q on %d monomials = %d, and the "
           "dimension predicted from Oeding's decomposition is %d"
           % (len(allP), len(base), len(extra), wz, nmon, rk, ctx['pred_w0']))
    ctx['HD0'] = allP
    ctx['HD0_rank'] = rk

    # 27b. SATURATION.  Check 27 shows only that the constructed span has the
    # dimension PREDICTED by Oeding's decomposition (input I8).  If that
    # citation were misquoted -- if the true HD(5) were larger -- the span could
    # be a PROPER subspace of HD(5)^{weight 0} of coincidentally matching rank,
    # and the separation of check 28 would be tested against the wrong space
    # while still printing PASS.  So the upper bound is computed here too.
    #
    # HD(5) is the span of the GL_2^5-orbit of the forty slices, so
    #    HD(5)^{wt 0} = sum over (T, fix) of pi_0( span of the orbit of D_{T,fix} ),
    # pi_0 being the (linear) projection to torus weight zero.  For each of the
    # ten 3-subsets T we sample that orbit by words of transvections and find
    # the weight-zero projection of the whole sampled orbit to be exactly ONE
    # dimensional -- saturated: no further sample enlarges it.  Ten subsets then
    # give at most 10, which together with the rank 10 of check 27 pins
    # dim HD(5)^{wt 0} = 10 without using the cited decomposition, and shows the
    # constructed space already CONTAINS every sampled orbit element.
    rng = random.Random(13579)
    per_T, orbit_samples = [], []
    for T in Tlist:
        rest = [j for j in range(N) if j not in T]
        fam = []
        for fix in itertools.product((0, 1), repeat=N - 3):
            D = embed_cayley_slice(T, fix)
            for _ in range(4):
                # the two moves that actually displace this slice (the others,
                # inside T, fix it: Cayley's hyperdeterminant is SL_2^3
                # invariant), plus random extra transvections in any factor.
                steps = [(rest[0], fix[0] == 0), (rest[1], fix[1] == 0)]
                for _ in range(rng.randint(1, 3)):
                    steps.append((rng.randrange(N), bool(rng.randrange(2))))
                rng.shuffle(steps)
                P = D
                for (j, hi) in steps:
                    P = transvection(P, j, rng.choice([1, -1, 2, -2, 3]), hi=hi)
                Z = project_weight0(P)
                if Z:
                    fam.append(Z)
        r = rank_exact(coeff_rows(fam)[0]) if fam else 0
        per_T.append((len(fam), r))
        orbit_samples.extend(fam)
    rk_orb = rank_exact(coeff_rows(orbit_samples)[0])
    rk_union = rank_exact(coeff_rows(allP + orbit_samples)[0])
    wz_orb = all(weight_vec(m) == z for P in orbit_samples for m in P)
    ok = (len(per_T) == math.comb(N, 3)
          and all(cnt > 0 and r == 1 for cnt, r in per_T)
          and wz_orb and rk_orb == rk and rk_union == rk)
    record("HD5_weight_zero_span_is_saturated_by_the_GL2^5_orbit", ok,
           "for each of the %d 3-subsets T the weight-zero projection of the "
           "sampled GL_2^5-orbit of the four slices D_{T,fix} has rank 1 "
           "(sample sizes %s), so dim HD(5)^{wt 0} <= %d; the %d samples "
           "together have rank %d and adjoining them to the %d constructed "
           "generators leaves the rank at %d -- the constructed space is not a "
           "proper subspace of the sampled orbit's weight-zero part"
           % (len(per_T), [c for c, _ in per_T], math.comb(N, 3),
              len(orbit_samples), rk_orb, len(allP), rk_union))
    ctx['HD0_saturated'] = ok
    ctx['orbit_samples'] = orbit_samples

    # 28. W meets the computed space in 0 -- the refutation, computed directly.
    rows, nmon = coeff_rows(allP + orbit_samples + Fs)
    rk2 = rank_exact(rows)
    record("W_meets_the_computed_HD5_weight_zero_space_in_zero",
           rk2 == rk + 4 and rk == ctx['pred_w0'] and ctx['HD0_saturated']
           and ctx['dim_W'] == 4,
           "rank(HD(5)^{weight 0} generators together with the %d sampled orbit "
           "elements) = %d, rank(same together with F_2..F_5) = %d = %d + 4 on "
           "%d monomials, and that space is saturated under the sampled orbit "
           "(check 27b), so span(F_2..F_5) cap HD(5) = 0"
           % (len(ctx['orbit_samples']), rk, rk2, rk, nmon))
    ctx['W_cap_HD5_zero_computed'] = (rk2 == rk + 4 and ctx['HD0_saturated']
                                      and ctx['dim_W'] == 4)

    # 29. The computed space carries no nonzero SL_2^5-invariant: solve the
    #     linear system lie_op(sum c_k v_k) = 0 for all ten operators.
    basis = []
    seen_rows = []
    for P in allP:                                # extract an independent basis
        cand = seen_rows + [P]
        rws, _ = coeff_rows(cand)
        if rank_exact(rws) == len(cand):
            seen_rows = cand
            basis.append(P)
        if len(basis) == rk:
            break
    imgs = []
    for j in range(N):
        for hi in (True, False):
            imgs.append([lie_op(P, j, hi) for P in basis])
    rowsM = []
    for group in imgs:
        mons = sorted(set(m for P in group for m in P))
        for m in mons:
            rowsM.append([P.get(m, 0) for P in group])
    rk3 = rank_exact(rowsM) if rowsM else 0
    inv_dim = len(basis) - rk3
    record("computed_HD5_weight_zero_space_has_no_G5_invariant",
           len(basis) == ctx['pred_w0'] and inv_dim == 0,
           "on the %d-dimensional computed space, the ten operators lie_lo_j, "
           "lie_hi_j have joint rank %d, so the invariant subspace has "
           "dimension %d" % (len(basis), rk3, inv_dim))
    ctx['HD0_basis'] = basis

    # 30. The conclusion, restated so that it uses NO number taken from the
    #     literature.  Check 24 quotes dim HD(5) = 250 and the printed bound
    #     254, both of which descend from Oeding's decomposition (I8); but the
    #     strict containment does not need either.  It needs exactly:
    #        (a) W <= (P_5)_4                 -- computed, check 19
    #        (b) dim W = 4 > 0                -- computed, check 22
    #        (c) W cap HD(5) = 0              -- computed, checks 27b + 28
    #        (d) (J_5)_4 = HD(5) <= (P_5)_4   -- input I9 (J_5 is generated in
    #            degree four, and Oeding Prop. 5.5), with the generators' side
    #            of it recomputed exactly in check 14.
    #     Then (J_5)_4 = HD(5) is a proper subspace of HD(5) + W <= (P_5)_4,
    #     so (J_5)_4 != (P_5)_4 and a fortiori J_5 != P_5, whatever
    #     dim HD(5) happens to be.
    ok = (ctx['F_in_P5'] and ctx['dim_W'] == 4
          and ctx['W_cap_HD5_zero_computed'])
    record("J5_neq_P5_from_computed_facts_only_no_dimension_cited", ok,
           "W <= (P_5)_4: %s; dim W = %d; W cap HD(5) = 0 by direct computation "
           "in weight zero (saturated: %s)  ==>  (J_5)_4 = HD(5) is strictly "
           "contained in HD(5) + W <= (P_5)_4, hence J_5 != P_5 -- with no "
           "appeal to the value 250 or to the printed bound 254"
           % (ctx['F_in_P5'], ctx['dim_W'], ctx['HD0_saturated']))


# ======================================================================
# GROUP 4 -- the census: the extension of the argument to n >= 5
# ======================================================================

CENSUS_N_MAX_RELATIONS = 12        # n range re-run for the relations F_i^S
CENSUS_N_MAX_MODULE = 14           # n range re-run for the HD(n) module facts


def checks_census(ctx):
    section("GROUP 4: census over n (the paper claims every n >= 5)")
    R5, PIp = ctx['R5'], ctx['PIp']

    # 25. For n = 6..N_MAX: for EVERY 5-subset S of [n] and EVERY i, the copy
    #     F_i^S is a relation, and B_i^S o phi_n = 12 t^{4(n-5)+10} Pi(X_S)^2.
    rng = random.Random(555)
    per_n = []
    bad_rel, bad_scale, tested, nz = 0, 0, 0, 0
    for n in range(6, CENSUS_N_MAX_RELATIONS + 1):
        subs = list(itertools.combinations(range(n), 5))
        cnt = 0
        for rep in range(3):
            X = random_sym(rng, n, -5, 5)
            t = rng.randint(1, 4)
            for S in subs:
                XS = [[X[a][b] for b in S] for a in S]
                piS = eval_xt(PIp, R5, XS, 1)
                target = (PULLBACK_CONST_PAPER * t ** (4 * (n - 5))
                          * t ** 10 * piS * piS)
                if target != 0:
                    nz += 1
                Av = {}
                for al in itertools.product((0, 1), repeat=5):
                    I = [S[k] for k in range(5) if al[k] == 1]
                    M = [[X[a][b] for b in I] for a in I]
                    Av[al] = (t ** (n - len(I))) * det_int(M)
                vals = []
                for i in range(5):
                    _, b = eval_q(Av, i, 5)
                    vals.append(b)
                    tested += 1
                    if b != target:
                        bad_scale += 1
                for i in range(1, 5):
                    if vals[0] - vals[i] != 0:
                        bad_rel += 1
                cnt += 1
        per_n.append((n, len(subs), cnt))
    ok = bad_rel == 0 and bad_scale == 0 and nz > 0
    record("extension_to_n_6..%d_every_5subset_every_i"
           % CENSUS_N_MAX_RELATIONS, ok,
           "(n, #5-subsets, #(S,point) pairs) = %s; %d exact evaluations: "
           "F_i^S(phi_n) != 0 in %d cases, B_i^S(phi_n) != 12 t^{4(n-5)+10} "
           "Pi(X_S)^2 in %d cases; %d evaluations had nonzero target"
           % (per_n, tested, bad_rel, bad_scale, nz))

    # 26. The module side for n = 3..N_MAX: dim HD(n) and HD(n)^{H_S} = 0.
    d22, d4 = ctx['d22'], ctx['d4']
    iv_in, iv_out = ctx['inv_in_T'], ctx['inv_out_T']
    dim_rows, bad_dim, bad_inv, pairs = [], 0, 0, 0
    min_out = None
    for n in range(3, CENSUS_N_MAX_MODULE + 1):
        tot = 0
        for T in itertools.combinations(range(n), 3):
            p = 1
            for j in range(n):
                p *= d22 if j in T else d4
            tot += p
        closed = math.comb(n, 3) * 5 ** (n - 3)
        if tot != closed:
            bad_dim += 1
        dim_rows.append((n, tot))
        if n >= 5:
            for T in itertools.combinations(range(n), 3):
                Ts = set(T)
                for S in itertools.combinations(range(n), 5):
                    out = len([j for j in S if j not in Ts])
                    min_out = out if min_out is None else min(min_out, out)
                    p = 1
                    for j in S:
                        p *= iv_in if j in Ts else iv_out
                    pairs += 1
                    if p != 0:
                        bad_inv += 1
    ok = bad_dim == 0 and bad_inv == 0 and min_out == 2 and pairs > 0
    record("HDn_dimension_and_H_S_invariants_census_n=3..%d"
           % CENSUS_N_MAX_MODULE, ok,
           "dim HD(n) == binom(n,3) 5^{n-3} for all n in 3..%d "
           "(mismatches=%d); over %d (3-subset T, 5-subset S) pairs the "
           "H_S-invariant dimension of the summand was nonzero %d times; "
           "min |S \\ T| = %s; dims %s"
           % (CENSUS_N_MAX_MODULE, bad_dim, pairs, bad_inv, min_out,
              dim_rows[:3] + ['...'] + dim_rows[-1:]))

    print("    NOT RE-RUN: the paper asserts J_n != P_n for EVERY n >= 5, "
          "an infinite family.")
    print("    Re-run here: n = 5 in full, exactly and symbolically (the "
          "exhibited object);")
    print("    the relations F_i^S for n = 6..%d, over every one of the "
          "%d 5-subsets S and"
          % (CENSUS_N_MAX_RELATIONS,
             sum(math.comb(n, 5) for n in range(6, CENSUS_N_MAX_RELATIONS + 1))))
    print("    every i, at exact integer points; and the HD(n) module facts "
          "(dimension and")
    print("    vanishing of H_S-invariants) for n = 3..%d over every "
          "(3-subset, 5-subset) pair." % CENSUS_N_MAX_MODULE)
    print("    NOT re-run: every n > %d, and the relations for n = %d..%d. "
          "Those rest on the"
          % (CENSUS_N_MAX_MODULE, CENSUS_N_MAX_RELATIONS + 1,
             CENSUS_N_MAX_MODULE))
    print("    paper's uniform argument (restriction of phi_n to a 5-subset), "
          "not on computation.")

    # The remaining shortfall against the paper: it also descends the
    # refutation to the REAL form of Conjecture 14 (Status paragraph).  Report
    # what this program does supply towards that descent -- integrality of every
    # coefficient it handles, scanned and counted here, and ranks over Q, which
    # for an integer matrix are unchanged over R -- and name the step it does
    # not supply.
    scanned, non_integral = 0, 0
    for P in ([ctx['PIp'], ctx['PId']] + list(ctx['B']) + list(ctx['F'])
              + list(ctx['HD0']) + list(ctx['orbit_samples'])):
        for c in P.values():
            scanned += 1
            if not isinstance(c, int):
                non_integral += 1
    print("    NOT re-run: the descent of the Status paragraph to the REAL "
          "form of Conjecture 14")
    print("    (ring R[A_.], group SL_2(R)^n).  Supplied here: every "
          "coefficient of the pentad, the")
    print("    B_i, the F_i and the %d constructed and sampled generators of "
          "HD(5)^{weight 0} is an"
          % (len(ctx['HD0']) + len(ctx['orbit_samples'])))
    print("    integer (%d coefficients scanned, %d non-integral), and every "
          "rank and kernel above is"
          % (scanned, non_integral))
    print("    computed over Q by exact elimination, so it is the rank of the "
          "same integer matrix")
    print("    over R.  NOT supplied: that HD(n) itself has a basis of "
          "quartics with integer")
    print("    coefficients.  That is taken from the paper (orbit span of an "
          "integral quartic under")
    print("    a group defined over Q) and not re-derived, so the real form is "
          "not re-verified here.")


# ======================================================================
# GROUP 5 -- negative controls: corrupt the object, demand a failure
# ======================================================================

def checks_negative_controls(ctx):
    section("GROUP 5: negative controls (each corruption must be detected)")
    R5, im, PIp, Bs, Fs = ctx['R5'], ctx['im'], ctx['PIp'], ctx['B'], ctx['F']
    rhs = rhs_pullback(R5, PIp, PULLBACK_CONST_PAPER, 10)
    outcomes = []

    def note(label, detected):
        outcomes.append((label, detected))
        print("    control %-38s detected=%s" % (label, detected))

    # (a) flip the sign of one pentad monomial
    m0 = sorted(PIp)[0]
    bad = dict(PIp)
    bad[m0] = -bad[m0]
    note("pentad_term_sign_flipped",
         subst_A_poly(R5, Bs[0], im)
         != rhs_pullback(R5, bad, PULLBACK_CONST_PAPER, 10))

    # (b) delete one pentad monomial
    bad = dict(PIp)
    del bad[m0]
    note("pentad_term_deleted",
         subst_A_poly(R5, Bs[0], im)
         != rhs_pullback(R5, bad, PULLBACK_CONST_PAPER, 10))

    # (c) wrong constant in 12 t^10 Pi^2
    note("pullback_constant_12_to_13",
         subst_A_poly(R5, Bs[0], im) != rhs_pullback(R5, PIp, 13, 10))

    # (d) wrong power of t
    note("pullback_t_exponent_10_to_9",
         subst_A_poly(R5, Bs[0], im)
         != rhs_pullback(R5, PIp, PULLBACK_CONST_PAPER, 9))

    # (e) make the alternating form symmetric
    eps_bad = {(0, 1): 1, (1, 0): 1, (0, 0): 0, (1, 1): 0}
    _, Bbad = build_Q_and_B(0, n=N, eps=eps_bad)
    note("epsilon_made_symmetric", subst_A_poly(R5, Bbad, im) != rhs)

    # (f) perturb a single coefficient of B_2, so F_2 is no longer a relation
    B2 = dict(Bs[1])
    mm = sorted(B2)[0]
    B2[mm] = B2[mm] + 1
    note("one_B_2_coefficient_perturbed",
         subst_A_poly(R5, asub(Bs[0], B2), im) != {})

    # (g) perturb one entry of the Lemma 3 matrix
    mat = [row[:] for row in LEMMA3_MATRIX_PAPER]
    mat[0][1] += 1
    note("lemma3_matrix_entry_perturbed", det_int(mat) != LEMMA3_DET_PAPER)

    # (h) pretend Sym^4 had an SL_2-invariant: HD(5)^{G_5} would be nonzero
    tot = 0
    for T in itertools.combinations(range(N), 3):
        p = 1
        for j in range(N):
            p *= ctx['inv_in_T'] if j in T else 1
        tot += p
    note("Sym4_invariant_dim_forced_to_1", tot != 0)

    # (i) replace the Hamiltonian-cycle sign rule by all +1
    allplus = {}
    for edges in hamiltonian_cycles_K5():
        m = 0
        for e in edges:
            i, j = sorted(e)
            m += R5.x(i - 1, j - 1)
        allplus[m] = 1
    note("hamiltonian_cycle_signs_dropped", allplus != PIp)

    # (j) structural corruption: in one pentad term, replace one variable by
    #     another (a relabelling of ALL indices would only change the sign of
    #     Pi, hence nothing in Pi^2, so the corruption must be local).
    e0 = unpk(m0, R5.nv)
    kin = min(k for k in range(R5.tvar) if e0[k] == 1)
    kout = min(k for k in range(R5.tvar) if e0[k] == 0)
    bad = dict(PIp)
    c0 = bad.pop(m0)
    newkey = m0 - (1 << (SH * kin)) + (1 << (SH * kout))
    bad[newkey] = bad.get(newkey, 0) + c0
    note("one_pentad_variable_substituted",
         subst_A_poly(R5, Bs[0], im)
         != rhs_pullback(R5, bad, PULLBACK_CONST_PAPER, 10))

    # (k) positive control on the membership test of Group 3b: adding a genuine
    #     element of HD(5)^{weight 0} must NOT raise the rank, whereas adding
    #     F_2 does.  So the "W cap HD(5) = 0" verdict is not an artefact of the
    #     rank computation always increasing.
    HD0, basis = ctx['HD0'], ctx['HD0_basis']
    combo = {}
    for k, P in enumerate(basis):
        for m, c in P.items():
            combo[m] = combo.get(m, 0) + (k + 1) * c
    combo = dict((m, c) for m, c in combo.items() if c)
    r_base = rank_exact(coeff_rows(HD0)[0])
    r_in = rank_exact(coeff_rows(HD0 + [combo])[0])
    r_out = rank_exact(coeff_rows(HD0 + [Fs[0]])[0])
    note("rank_test_admits_HD5_elements_and_rejects_F_2",
         r_base == r_in and r_out == r_base + 1)

    # (l) the SATURATION test of check 27b must be able to fail: adjoining a
    #     non-member (F_2, or any weight-zero element outside HD(5)) to a
    #     per-T orbit family must raise that family's rank above 1.
    samp = ctx['orbit_samples']
    fam0 = [samp[0]]
    note("saturation_rank1_test_rejects_a_non_member",
         rank_exact(coeff_rows(fam0)[0]) == 1
         and rank_exact(coeff_rows(fam0 + [Fs[0]])[0]) == 2)

    # (m) the exact symbolic pullback of check 14 must be able to fail: a single
    #     perturbed coefficient in a Cayley slice must stop it vanishing.
    Dbad = embed_cayley_slice((0, 1, 2), (0, 0))
    mb = sorted(Dbad)[0]
    Dbad = dict(Dbad)
    Dbad[mb] = Dbad[mb] + 1
    note("perturbed_cayley_slice_no_longer_pulls_back_to_zero",
         subst_A_poly(R5, Dbad, im) != {})

    ndet = sum(1 for _, d in outcomes if d)
    record("negative_controls_all_detected", ndet == len(outcomes),
           "%d of %d deliberate corruptions of the exhibited object made a "
           "load-bearing check fail" % (ndet, len(outcomes)))


# ======================================================================
# Driver
# ======================================================================

def build_context():
    ctx = {}
    ctx['R5'] = SymRing(N)
    ctx['im'] = phi_images(ctx['R5'])
    Qs, Bs = [], []
    for i in range(N):
        Q, B = build_Q_and_B(i, n=N)
        Qs.append(Q)
        Bs.append(B)
    ctx['Q'], ctx['B'] = Qs, Bs
    ctx['PIp'] = pentad_from_paper(ctx['R5'])
    PId, wd, nc = pentad_derived(ctx['R5'])
    ctx['PId'], ctx['PI_welldef'], ctx['PI_ncyc'] = PId, wd, nc
    return ctx


def print_object(ctx):
    section("THE EXHIBITED OBJECT, PRINTED BACK")
    print("n = %d;  ring C[A_I : I subset of [5]] has %d coordinates."
          % (N, len(ALPHAS)))
    print("epsilon: eps_01 = %d, eps_10 = %d, eps_00 = %d, eps_11 = %d."
          % (EPS[(0, 1)], EPS[(1, 0)], EPS[(0, 0)], EPS[(1, 1)]))
    print("Q_1 entries (as quadratic forms in A), monomial counts: "
          "q00=%d q01=%d q10=%d q11=%d"
          % tuple(len(ctx['Q'][0][k]) for k in ((0, 0), (0, 1), (1, 0), (1, 1))))
    print("B_i = det Q_i: %s monomials for i = 1..5; coefficients seen: %s"
          % ([len(B) for B in ctx['B']],
             sorted(set(c for B in ctx['B'] for c in B.values()))))
    print("pentad Pi: %d terms, one per Hamiltonian cycle of K_5; for example"
          % len(PENTAD_PAPER))
    for c, edges in PENTAD_PAPER[:2]:
        print("    %+d * %s" % (c, " ".join("x_%d%d" % e for e in edges)))
    print("F_i = B_1 - B_i, i = 2,3,4,5.")


def main():
    print("verify.py -- The Holtz-Sturmfels ideal-generation conjecture "
          "fails for n >= 5")
    print("exact integer arithmetic, standard library only, no external data")
    ctx = build_context()
    print_object(ctx)
    checks_object(ctx)
    checks_hypotheses(ctx)
    checks_loadbearing(ctx)
    checks_hd5_weight_zero(ctx)
    checks_census(ctx)
    checks_negative_controls(ctx)

    total = len(RESULTS)
    failed = [nm for nm, ok in RESULTS if not ok]
    print("")
    if failed:
        print("failed checks: %s" % ", ".join(failed))
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


if __name__ == "__main__":
    sys.exit(main())

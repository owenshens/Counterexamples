#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verifier for "Supermodular Log-Submodular Valuations Obey the Jaccard
Triangle Inequality on Every Lattice", which answers affirmatively the second
open problem of Badica--Badica, arXiv:2608.18194v1, Section 7 (file
main-arxiv-v4.tex, line 631).

Python 3.9+, standard library only, exact integer / Fraction / Z[sqrt 5]
arithmetic.  No float ever decides anything.

=====================================================================
TAKEN FROM THE PAPER (inputs -- these are not checked, they ARE the
objects the paper exhibits, or the numbers the paper claims)
=====================================================================
  * The COVER RELATIONS of the three lattices printed in Section 3:
        M_3 : 0 < a, 0 < b, 0 < c, a < 1, b < 1, c < 1
        N_5 : 0 < a, 0 < b, b < c, a < 1, c < 1
        B_2 : 0 < x, 0 < y, x < 1, y < 1
    Only the covers are hard-coded.  The order, the joins, the meets,
    distributivity, relative complementation are all derived here.
  * The four valuations printed on B_2 in Section 3, as the tuple
    (f(0), f(x), f(y), f(1)):
        f_1 = (2, 3, 3, 20)     f_2 = (1, 3, 3, 4)
        f_3 = (1, 2, 3, 4)      f_4 = (1, 1, 1, 1)
    and the two claimed violation margins 7/30 (at the triple (x,0,y))
    and 1/4 (at the triple (x,1,y)).
  * The near-miss valuation g on B_2 printed in Section 4,
        g(0) = (5 - sqrt 5)/2,  g(x) = g(y) = 2,  g(1) = (5 + sqrt 5)/2,
    and the claim d_g(x,0) + d_g(0,y) = d_g(x,y) = (sqrt 5 - 1)/2.
  * The claimed sweep totals of Section 3 (subsection "A finite control"):
    over the 22 isomorphism classes
    of lattices of order 4, 5 and 6 and the value grid {1,2,3,4,6,9}
    normalised by f(0) = 1, there are 1814 admissible valuations, the
    largest triangle-inequality margin is exactly 0, and neither Lemma 1
    nor Lemma 2 ever fails.

  NAMING, so that every label below can be located in the paper: what this
  program calls Lemma 1 and Lemma 2 are the two inequalities the paper
  establishes in Step 1 and Step 2 of the proof of its Theorem 1, displayed
  there as equations (6)  r >= p + q - b  and (7)  R <= PQ/b.  The paper has
  no numbered Lemma 1 or Lemma 2; the names are this program's shorthand.
  * The enumerator control: the numbers of isomorphism classes of lattices
    of order n = 2,...,6 are 1, 1, 2, 5, 15 (OEIS A006966).

=====================================================================
DERIVED HERE (the checks; every one can fail if the paper is wrong)
=====================================================================
  A. The three printed cover sets really do define lattices; M_3 and N_5
     really are non-distributive of order 5; M_3 really is relatively
     complemented and N_5 really is not (so M_3 lacks exactly one of the
     two traits named in the open problem, N_5 lacks both).
  B. The enumerator reproduces OEIS A006966 for n = 2..6; every lattice of
     order at most 4 is a chain or is distributive AND relatively
     complemented, so |L| <= 4 is already settled by the source's own two
     theorems; and M_3, N_5 are the ONLY non-distributive lattices of
     order 5, i.e. exactly the |L| = 5 instance of the open problem.
  C. The two algebraic identities that carry Steps 4 and 5 of the proof
     are verified as EXACT POLYNOMIAL IDENTITIES in Z[p,P,q,Q,b], by
     sparse expansion -- not by sampling.
  D. The four order-theoretic side conditions Lemma 1 and Lemma 2 use
     ((A^B) v (B^C) <= B, A^B^C <= A^C, (AvB) ^ (BvC) >= B,
     A v B v C >= A v C) hold at every ordered triple of every isomorphism
     class of lattices of order at most 6.
  E. Both hypotheses are necessary for the two lemmas: dropping
     log-submodularity (f_1) or supermodularity (f_2) breaks the triangle
     inequality on B_2 with exactly the margins the paper prints, and
     names exactly the lemma the paper says it breaks; while the two
     admissible controls f_3, f_4 have margin exactly 0.
  F. The sweep of Section 3 is re-run in cleared integer arithmetic and
     its three totals are compared with the paper's.
  G. The near-miss valuation g is re-checked in exact Z[sqrt 5]: monotone,
     supermodular, reciprocal supermodular, NOT log-submodular, and it
     satisfies the triangle inequality at all 64 triples with equality at
     the printed one.  This is why the paper claims sufficiency for a
     strictly smaller class and NOT a converse or a characterisation.
"""
import itertools
import sys
from fractions import Fraction

CHECKS = 0
FAILED = 0


def ok(name, detail=''):
    global CHECKS
    CHECKS += 1
    print('PASS %s [%s]' % (name, detail))


def bad(name, detail=''):
    global FAILED
    FAILED += 1
    print('FAIL %s [%s]' % (name, detail))


def note(text):
    print('NOTE %s' % text)


def claim(name, cond, detail=''):
    (ok if cond else bad)(name, detail)


# =====================================================================
# 1. LATTICES FROM COVER RELATIONS
# =====================================================================
def order_from_covers(n, covers):
    """Reflexive-transitive closure of a cover list.  leq[j] = bitmask of the
    elements <= j.  Returns None if the closure is not antisymmetric."""
    leq = [1 << j for j in range(n)]
    for (i, j) in covers:
        leq[j] |= 1 << i
    changed = True
    while changed:
        changed = False
        for j in range(n):
            m = leq[j]
            for i in range(n):
                if (m >> i) & 1 and (leq[i] & ~m):
                    leq[j] |= leq[i]
                    m = leq[j]
                    changed = True
    for i in range(n):
        for j in range(n):
            if i != j and (leq[j] >> i) & 1 and (leq[i] >> j) & 1:
                return None                      # a 2-cycle: not a partial order
    return leq


def lattice_tables(n, leq):
    """(join, meet) tables, or None when some pair has no unique bound."""
    up = [0] * n
    for i in range(n):
        for j in range(n):
            if (leq[j] >> i) & 1:
                up[i] |= 1 << j
    join = [[0] * n for _ in range(n)]
    meet = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            U = up[i] & up[j]                    # common upper bounds
            cand = [k for k in range(n) if (U >> k) & 1]
            mins = [k for k in cand if not any(k2 != k and ((leq[k] >> k2) & 1) for k2 in cand)]
            if len(mins) != 1:
                return None
            join[i][j] = mins[0]
            L = leq[i] & leq[j]                  # common lower bounds
            cand = [k for k in range(n) if (L >> k) & 1]
            maxs = [k for k in cand if not any(k2 != k and ((leq[k2] >> k) & 1) for k2 in cand)]
            if len(maxs) != 1:
                return None
            meet[i][j] = maxs[0]
    return join, meet


def build(n, covers):
    leq = order_from_covers(n, covers)
    if leq is None:
        return None
    t = lattice_tables(n, leq)
    if t is None:
        return None
    return (leq, t[0], t[1])


def le(leq, i, j):
    return bool((leq[j] >> i) & 1)


def is_chain(n, leq):
    return all(le(leq, i, j) or le(leq, j, i) for i in range(n) for j in range(n))


def is_distributive(n, join, meet):
    for x in range(n):
        for y in range(n):
            for z in range(n):
                if meet[x][join[y][z]] != join[meet[x][y]][meet[x][z]]:
                    return False, (x, y, z)
    return True, None


def is_relatively_complemented(n, leq, join, meet):
    """Every element of every interval [u,w] has a complement in that interval."""
    for u in range(n):
        for w in range(n):
            if not le(leq, u, w):
                continue
            for v in range(n):
                if not (le(leq, u, v) and le(leq, v, w)):
                    continue
                found = False
                for t in range(n):
                    if not (le(leq, u, t) and le(leq, t, w)):
                        continue
                    if meet[v][t] == u and join[v][t] == w:
                        found = True
                        break
                if not found:
                    return False, (u, v, w)
    return True, None


# =====================================================================
# 2. ENUMERATION OF ALL LATTICES OF SMALL ORDER (the control)
# =====================================================================
def naturally_labelled_orders(n):
    """Every finite poset has a linear extension, so up to isomorphism we may
    assume i < j whenever x_i < x_j.  Enumerate the strict down-sets."""
    res = []

    def rec(j, D):
        if j == n:
            res.append(list(D))
            return
        for m in range(1 << j):
            good = True
            for i in range(j):
                if (m >> i) & 1 and (D[i] & ~m):
                    good = False
                    break
            if good:
                rec(j + 1, D + [m])
    rec(0, [])
    return res


def canonical_form(n, leq):
    best = None
    for perm in itertools.permutations(range(n)):
        rel = tuple(sorted((perm[i], perm[j])
                           for i in range(n) for j in range(n) if (leq[j] >> i) & 1))
        if best is None or rel < best:
            best = rel
    return best


def all_lattices(n):
    """(number of naturally-labelled lattices, one representative per iso class)."""
    labelled = 0
    reps = {}
    for D in naturally_labelled_orders(n):
        leq = [D[j] | (1 << j) for j in range(n)]
        t = lattice_tables(n, leq)
        if t is None:
            continue
        labelled += 1
        c = canonical_form(n, leq)
        if c not in reps:
            reps[c] = (leq, t[0], t[1])
    return labelled, [reps[k] for k in sorted(reps)]


def isomorphic(n, A, B):
    return canonical_form(n, A[0]) == canonical_form(n, B[0])


# =====================================================================
# 3. VALUATIONS: the hypotheses and the triangle inequality, cleared of
#    all denominators so that integer inputs stay integer.
# =====================================================================
def hypothesis_failure(n, leq, join, meet, f):
    """The FIRST hypothesis f breaks, or None.  Works over any ordered ring."""
    for i in range(n):
        for j in range(n):
            if le(leq, i, j) and f[i] > f[j]:
                return 'MONOTONE'
    for i in range(n):
        for j in range(i + 1, n):
            J, M = f[join[i][j]], f[meet[i][j]]
            if J + M < f[i] + f[j]:
                return 'SUPERMODULAR'
            if J * M > f[i] * f[j]:
                return 'LOGSUBMODULAR'
    return None


def cleared_margin(join, meet, f, A, B, C):
    """P*Q*R * ( S(A,B) + S(B,C) - 1 - S(A,C) ) with S(X,Y)=f(X^Y)/f(XvY).
    All of P,Q,R are > 0, so the sign of this equals the sign of the margin,
    and the triangle inequality at (A,B,C) is exactly `<= 0`."""
    p, P = f[meet[A][B]], f[join[A][B]]
    q, Q = f[meet[B][C]], f[join[B][C]]
    r, R = f[meet[A][C]], f[join[A][C]]
    return p * Q * R + q * P * R - P * Q * R - r * P * Q


def true_margin(join, meet, f, A, B, C):
    """The margin itself, as an exact Fraction (used only for printing)."""
    p, P = f[meet[A][B]], f[join[A][B]]
    q, Q = f[meet[B][C]], f[join[B][C]]
    r, R = f[meet[A][C]], f[join[A][C]]
    return Fraction(p, P) + Fraction(q, Q) - 1 - Fraction(r, R)


def lemma_failure(n, join, meet, f, A, B, C):
    """Which of the paper's two lemmas fails at this triple, or None."""
    p, P = f[meet[A][B]], f[join[A][B]]
    q, Q = f[meet[B][C]], f[join[B][C]]
    r, R = f[meet[A][C]], f[join[A][C]]
    b = f[B]
    if r < p + q - b:
        return 'LEMMA1'
    if R * b > P * Q:
        return 'LEMMA2'
    return None


def worst_triple(n, join, meet, f):
    worst, arg = None, None
    for A in range(n):
        for B in range(n):
            for C in range(n):
                v = cleared_margin(join, meet, f, A, B, C)
                if worst is None or v > worst:
                    worst, arg = v, (A, B, C)
    return worst, arg


GRID = (1, 2, 3, 4, 6, 9)


def sweep(n, L):
    """(admissible count, max cleared margin, violations, lemma failures)."""
    leq, join, meet = L
    bot = min(range(n), key=lambda k: bin(leq[k]).count('1'))
    free = [k for k in range(n) if k != bot]
    nadm = 0
    worst = None
    viol = 0
    lem = 0
    for vals in itertools.product(GRID, repeat=len(free)):
        f = [1] * n
        for k, v in zip(free, vals):
            f[k] = v
        if hypothesis_failure(n, leq, join, meet, f) is not None:
            continue
        nadm += 1
        for A in range(n):
            for B in range(n):
                for C in range(n):
                    m = cleared_margin(join, meet, f, A, B, C)
                    if worst is None or m > worst:
                        worst = m
                    if m > 0:
                        viol += 1
                    if lemma_failure(n, join, meet, f, A, B, C) is not None:
                        lem += 1
    return nadm, worst, viol, lem


# =====================================================================
# 4. EXACT SPARSE POLYNOMIALS IN Z[p,P,q,Q,b] -- for Steps 4 and 5
# =====================================================================
VARS = ('p', 'P', 'q', 'Q', 'b')


class Poly(object):
    __slots__ = ('t',)

    def __init__(self, t=None):
        self.t = dict(t or {})

    @staticmethod
    def var(name):
        e = [0] * len(VARS)
        e[VARS.index(name)] = 1
        return Poly({tuple(e): 1})

    @staticmethod
    def const(c):
        return Poly({(0,) * len(VARS): c}) if c else Poly()

    def __add__(self, o):
        out = dict(self.t)
        for m, c in o.t.items():
            out[m] = out.get(m, 0) + c
            if out[m] == 0:
                del out[m]
        return Poly(out)

    def __neg__(self):
        return Poly(dict((m, -c) for m, c in self.t.items()))

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        out = {}
        for m1, c1 in self.t.items():
            for m2, c2 in o.t.items():
                m = tuple(a + b for a, b in zip(m1, m2))
                out[m] = out.get(m, 0) + c1 * c2
                if out[m] == 0:
                    del out[m]
        return Poly(out)

    def __eq__(self, o):
        return self.t == o.t

    def is_zero(self):
        return not self.t

    def __str__(self):
        if not self.t:
            return '0'
        parts = []
        for m in sorted(self.t, reverse=True):
            c = self.t[m]
            body = ''.join(VARS[i] * m[i] for i in range(len(VARS)))
            parts.append(('%+d' % c) + (body or ''))
        return ''.join(parts)


# =====================================================================
# 5. EXACT ARITHMETIC IN Q(sqrt 5) -- for the near-miss valuation
# =====================================================================
class Q5(object):
    """a + b*sqrt(5) with a, b rational.  Exact, ordered, no floats."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(self.a + o.a, self.b + o.b)

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Q5) else Q5(o)))

    def __mul__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__
    __radd__ = __add__

    def sign(self):
        """Exact sign of a + b*sqrt5.  sqrt5 is irrational, so a + b*sqrt5 = 0
        only when a = b = 0; otherwise compare a^2 with 5b^2."""
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if a > 0 and b > 0:
            return 1
        if a < 0 and b < 0:
            return -1
        d = a * a - 5 * b * b               # sign of (a - b sqrt5)(a + b sqrt5)
        if d == 0:
            return 0
        if a > 0:                           # b < 0: a + b sqrt5 > 0 iff a^2 > 5b^2
            return 1 if d > 0 else -1
        return -1 if d > 0 else 1           # a < 0, b > 0

    def __gt__(self, o):
        return (self - (o if isinstance(o, Q5) else Q5(o))).sign() > 0

    def __lt__(self, o):
        return (self - (o if isinstance(o, Q5) else Q5(o))).sign() < 0

    def __eq__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return self.a == o.a and self.b == o.b

    def __str__(self):
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return '%s*sqrt5' % self.b
        return '%s %s %s*sqrt5' % (self.a, '+' if self.b > 0 else '-', abs(self.b))


# =====================================================================
# THE OBJECTS PRINTED IN THE PAPER
# =====================================================================
# Section 3, cover relations.  0 = bottom, 4 (resp. 3) = top.
M3_COVERS = [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (3, 4)]
N5_COVERS = [(0, 1), (0, 2), (2, 3), (1, 4), (3, 4)]
B2_COVERS = [(0, 1), (0, 2), (1, 3), (2, 3)]
M3_NAMES = ('0', 'a', 'b', 'c', '1')
N5_NAMES = ('0', 'a', 'b', 'c', '1')
B2_NAMES = ('0', 'x', 'y', '1')

# Section 3, the four valuations on B_2 as (f(0), f(x), f(y), f(1)).
F1 = (2, 3, 3, 20)
F2 = (1, 3, 3, 4)
F3 = (1, 2, 3, 4)
F4 = (1, 1, 1, 1)
CLAIM_MARGIN_1 = Fraction(7, 30)
CLAIM_MARGIN_2 = Fraction(1, 4)
CLAIM_TRIPLE_1 = (1, 0, 2)          # (x, 0, y)
CLAIM_TRIPLE_2 = (1, 3, 2)          # (x, 1, y)

# Section 4, the near-miss valuation g on B_2.
G = (Q5(Fraction(5, 2), Fraction(-1, 2)), Q5(2), Q5(2), Q5(Fraction(5, 2), Fraction(1, 2)))
CLAIM_G_DIST = Q5(Fraction(-1, 2), Fraction(1, 2))     # (sqrt5 - 1)/2

# Section 3 ("A finite control"), the claimed sweep totals and enumerator control.
CLAIM_SWEEP_ADMISSIBLE = 1814
CLAIM_ISO_COUNTS = {2: 1, 3: 1, 4: 2, 5: 5, 6: 15}


def main():
    print('verify.py -- Supermodular Log-Submodular Valuations Obey the Jaccard')
    print('             Triangle Inequality on Every Lattice')
    print('             (answers Badica-Badica arXiv:2608.18194v1, Sec. 7, open problem 2)')
    print('stdlib only; exact integer / Fraction / Z[sqrt5] arithmetic; python 3.9+')
    print('')

    # ---------------- A. the printed objects are lattices ----------------
    print('=== A. the three lattices printed in Section 3 ===')
    built = {}
    for name, n, covers, names in (('M3', 5, M3_COVERS, M3_NAMES),
                                   ('N5', 5, N5_COVERS, N5_NAMES),
                                   ('B2', 4, B2_COVERS, B2_NAMES)):
        L = build(n, covers)
        built[name] = L
        if L is None:
            bad('%s_covers_define_a_lattice' % name.lower(), 'the closure is not a lattice')
            continue
        leq, join, meet = L
        nrel = sum(bin(m).count('1') for m in leq)
        bots = [k for k in range(n) if bin(leq[k]).count('1') == 1]
        tops = [k for k in range(n) if all(le(leq, i, k) for i in range(n))]
        ok('%s_covers_define_a_lattice' % name.lower(),
           '|L|=%d, %d order relations, unique bottom %s, unique top %s, all %d pairs have a '
           'unique join and a unique meet'
           % (n, nrel, names[bots[0]] if len(bots) == 1 else '?',
              names[tops[0]] if len(tops) == 1 else '?', n * n))
    M3, N5, B2 = built['M3'], built['N5'], built['B2']

    for name, L, names in (('m3', M3, M3_NAMES), ('n5', N5, N5_NAMES)):
        n = 5
        leq, join, meet = L
        d, w = is_distributive(n, join, meet)
        claim('%s_is_not_distributive' % name, not d,
              'witness (x,y,z)=(%s,%s,%s): x^(yvz)=%s but (x^y)v(x^z)=%s'
              % (names[w[0]], names[w[1]], names[w[2]],
                 names[meet[w[0]][join[w[1]][w[2]]]],
                 names[join[meet[w[0]][w[1]]][meet[w[0]][w[2]]]]) if w else 'distributive!')
    leq, join, meet = M3
    rc, wit = is_relatively_complemented(5, leq, join, meet)
    claim('m3_is_relatively_complemented', rc,
          'every element of every interval has a complement in that interval, so M_3 lacks '
          'exactly ONE of the two traits the open problem names (distributivity)')
    leq, join, meet = N5
    rc5, wit5 = is_relatively_complemented(5, leq, join, meet)
    claim('n5_is_not_relatively_complemented', not rc5,
          'no relative complement for %s in the interval [%s,%s], so N_5 lacks BOTH traits'
          % (N5_NAMES[wit5[1]], N5_NAMES[wit5[0]], N5_NAMES[wit5[2]]) if wit5 else 'it is!')

    # ---------------- B. the enumerator, and why |L|=5 is the first open size ----
    print('')
    print('=== B. the enumerator control, and the location of the first open size ===')
    reps = {}
    labelled = {}
    for n in (1, 2, 3, 4, 5, 6):
        labelled[n], reps[n] = all_lattices(n)
    got = dict((n, len(reps[n])) for n in reps if n >= 2)
    claim('lattice_iso_class_counts_match_oeis_a006966', got == CLAIM_ISO_COUNTS,
          'n=2..6 -> %s (OEIS A006966 a(2..6) = 1,1,2,5,15); naturally-labelled counts %s'
          % (', '.join(str(got[n]) for n in sorted(got)),
             ', '.join(str(labelled[n]) for n in sorted(got))))

    settled = []
    for n in (1, 2, 3, 4):
        for L in reps[n]:
            leq, join, meet = L
            ch = is_chain(n, leq)
            dist, _ = is_distributive(n, join, meet)
            rcp, _ = is_relatively_complemented(n, leq, join, meet)
            settled.append(ch or (dist and rcp))
    claim('every_lattice_of_order_at_most_four_is_settled_by_the_sources_own_theorems',
          all(settled),
          '%d of %d isomorphism classes of order <= 4 are a chain (every valuation is then '
          'modular, so th:mod applies) or are distributive AND relatively complemented (so '
          'th:super-log-sub applies)' % (sum(1 for s in settled if s), len(settled)))

    nd5 = [L for L in reps[5] if not is_distributive(5, L[1], L[2])[0]]
    hit_m3 = any(isomorphic(5, L, M3) for L in nd5)
    hit_n5 = any(isomorphic(5, L, N5) for L in nd5)
    claim('m3_and_n5_are_exactly_the_nondistributive_lattices_of_order_five',
          len(nd5) == 2 and hit_m3 and hit_n5,
          '%d of the 5 isomorphism classes of order 5 are non-distributive, and they are '
          'isomorphic to the printed M_3 (%s) and N_5 (%s); so |L|=5 is the first order the '
          'source\'s two theorems leave open, and M_3, N_5 are that whole instance'
          % (len(nd5), hit_m3, hit_n5))

    # ---------------- C. Steps 4 and 5 as exact polynomial identities ------
    print('')
    print('=== C. Steps 4 and 5 of the proof, as exact identities in Z[p,P,q,Q,b] ===')
    p, P, q, Q, b = (Poly.var(v) for v in VARS)
    one = Poly.const(1)
    # Step 4.  The sufficient inequality is  p/P + q/Q <= 1 + b(p+q-b)/(PQ).
    # Multiplying by PQ > 0 turns it into  pQ + qP <= PQ + b(p+q-b).
    cleared = (p * Q + q * P) - (P * Q + b * (p + q - b))
    target = -((P - p) * (Q - q) - (b - p) * (b - q))
    claim('step4_clearing_denominators_is_exactly_the_product_inequality',
          (cleared - target).is_zero(),
          'pQ+qP-PQ-b(p+q-b) = -[(P-p)(Q-q)-(b-p)(b-q)] identically; both sides expand to %s, '
          'so p/P+q/Q <= 1+b(p+q-b)/(PQ) IS (P-p)(Q-q) >= (b-p)(b-q)' % cleared)
    # Step 5.  The nonnegativity certificate: a sum of products of nonnegatives.
    diff = (P - p) * (Q - q) - (b - p) * (b - q)
    cert = (P - b) * (Q - q) + (b - p) * (Q - b)
    claim('step5_nonnegativity_certificate_is_an_identity',
          (diff - cert).is_zero(),
          '(P-p)(Q-q)-(b-p)(b-q) = (P-b)(Q-q)+(b-p)(Q-b) identically; monotonicity gives '
          'P>=b, Q>=b, b>=p, b>=q, so every one of the four factors is >= 0 and the sum is '
          '>= 0 termwise -- no case analysis, no ordering assumption beyond p,q <= b <= P,Q')
    claim('step4_and_step5_compose_to_the_triangle_inequality',
          (cleared + cert).is_zero(),
          'adding the two identities gives pQ+qP-PQ-b(p+q-b) = -[(P-b)(Q-q)+(b-p)(Q-b)] <= 0, '
          'which is the cleared form of S(A,B)+S(B,C) <= 1 + b(p+q-b)/(PQ) <= 1 + S(A,C)')

    # ---------------- D. the order-theoretic side conditions ---------------
    print('')
    print('=== D. the four order facts Lemma 1 and Lemma 2 use, on every lattice |L| <= 6 ===')
    print('    (Lemma 1 = the paper\'s Step 1, its equation (6) r >= p+q-b; '
          'Lemma 2 = the paper\'s Step 2, its equation (7) R <= PQ/b. '
          'The paper numbers these as STEPS, not as lemmas.)')
    counts = {'d1': 0, 'd2': 0, 'd3': 0, 'd4': 0}
    fails = {'d1': 0, 'd2': 0, 'd3': 0, 'd4': 0}
    triples = 0
    nclasses = sum(len(reps[n]) for n in (1, 2, 3, 4, 5, 6))
    for n in (1, 2, 3, 4, 5, 6):
        for (leq, join, meet) in reps[n]:
            for A in range(n):
                for B in range(n):
                    for C in range(n):
                        triples += 1
                        ab, bc = meet[A][B], meet[B][C]
                        AB, BC = join[A][B], join[B][C]
                        for key, cond in (
                                ('d1', le(leq, join[ab][bc], B)),
                                ('d2', le(leq, meet[meet[A][B]][C], meet[A][C])),
                                ('d3', le(leq, B, meet[AB][BC])),
                                ('d4', le(leq, join[A][C], join[join[A][B]][C]))):
                            if cond:
                                counts[key] += 1
                            else:
                                fails[key] += 1
    claim('lemma1_order_side_conditions_hold_at_every_triple',
          fails['d1'] == 0 and fails['d2'] == 0,
          '(A^B)v(B^C) <= B at %d/%d triples and A^B^C <= A^C at %d/%d triples, over all %d '
          'isomorphism classes of order <= 6'
          % (counts['d1'], triples, counts['d2'], triples, nclasses))
    claim('lemma2_order_side_conditions_hold_at_every_triple',
          fails['d3'] == 0 and fails['d4'] == 0,
          '(AvB)^(BvC) >= B at %d/%d triples and AvC <= AvBvC at %d/%d triples'
          % (counts['d3'], triples, counts['d4'], triples))

    # ---------------- E. both hypotheses are needed ------------------------
    print('')
    print('=== E. the two forced positives: drop either hypothesis and the TI breaks ===')
    n = 4
    leq, join, meet = B2
    for cname, f, want_fail, want_triple, want_margin, want_lemma in (
            ('control_a1_without_log_submodularity', F1, 'LOGSUBMODULAR',
             CLAIM_TRIPLE_1, CLAIM_MARGIN_1, 'LEMMA2'),
            ('control_a2_without_supermodularity', F2, 'SUPERMODULAR',
             CLAIM_TRIPLE_2, CLAIM_MARGIN_2, 'LEMMA1')):
        why = hypothesis_failure(n, leq, join, meet, list(f))
        w, arg = worst_triple(n, join, meet, list(f))
        m = true_margin(join, meet, list(f), *want_triple)
        lf = lemma_failure(n, join, meet, list(f), *want_triple)
        good = (why == want_fail and w > 0 and m == want_margin and lf == want_lemma)
        claim(cname, good,
              'f=(%s) on B_2 fails exactly %s; margin at (%s,%s,%s) is %s (paper: %s) > 0, so '
              'the triangle inequality is VIOLATED, and the step that breaks is %s (paper: %s)'
              % (','.join(str(x) for x in f), why,
                 B2_NAMES[want_triple[0]], B2_NAMES[want_triple[1]], B2_NAMES[want_triple[2]],
                 m, want_margin, lf, want_lemma))
    for cname, f in (('control_b1_modular_valuation_is_admissible', F3),
                     ('control_b2_constant_valuation_is_admissible', F4)):
        why = hypothesis_failure(n, leq, join, meet, list(f))
        w, arg = worst_triple(n, join, meet, list(f))
        claim(cname, why is None and w == 0,
              'f=(%s) on B_2 is admissible and the largest cleared margin over all 64 triples '
              'is exactly %d, so the triangle inequality holds everywhere with equality '
              'somewhere' % (','.join(str(x) for x in f), w))

    # ---------------- F. the sweep -----------------------------------------
    print('')
    print('=== F. the finite sweep of Section 3 (a control, not the proof) ===')
    for label, L, nn in (('m3', M3, 5), ('n5', N5, 5)):
        nadm, worst, viol, lem = sweep(nn, L)
        claim('sweep_%s_no_violation' % label, viol == 0 and lem == 0 and worst == 0,
              '%d admissible valuations on the printed %s over the grid %s with f(0)=1: '
              'largest cleared margin %d, triangle-inequality violations %d, lemma failures %d'
              % (nadm, label.upper(), ','.join(str(g) for g in GRID), worst, viol, lem))
    total = 0
    worst_all = None
    viol_all = 0
    lem_all = 0
    nlat = 0
    for n in (4, 5, 6):
        for L in reps[n]:
            nlat += 1
            nadm, worst, viol, lem = sweep(n, L)
            total += nadm
            viol_all += viol
            lem_all += lem
            if worst_all is None or worst > worst_all:
                worst_all = worst
    claim('sweep_orders_four_to_six_reproduces_the_papers_three_totals',
          total == CLAIM_SWEEP_ADMISSIBLE and worst_all == 0 and viol_all == 0 and lem_all == 0,
          '%d isomorphism classes, %d admissible valuations (paper: %d), largest cleared '
          'margin %d (paper: 0), violations %d, lemma failures %d'
          % (nlat, total, CLAIM_SWEEP_ADMISSIBLE, worst_all, viol_all, lem_all))

    # ---------------- G. the near-miss valuation ---------------------------
    print('')
    print('=== G. the near-miss valuation g of Section 4, exactly in Z[sqrt5] ===')
    leq, join, meet = B2
    g = list(G)
    mono = all((not le(leq, i, j)) or (not (g[i] > g[j])) for i in range(4) for j in range(4))
    claim('near_miss_g_is_monotone', mono,
          'g(0)=%s <= g(x)=g(y)=2 <= g(1)=%s' % (G[0], G[3]))
    sup = True
    worst_sup = None
    for i in range(4):
        for j in range(4):
            d = (g[join[i][j]] + g[meet[i][j]]) - (g[i] + g[j])
            if d.sign() < 0:
                sup = False
            if worst_sup is None or d < worst_sup:
                worst_sup = d
    claim('near_miss_g_is_supermodular', sup,
          'g(1)+g(0) = %s >= 4 = g(x)+g(y) on the only incomparable pair; smallest '
          'supermodularity slack over all 16 pairs is %s' % (G[0] + G[3], worst_sup))
    # 1/g supermodular:  1/g(i)+1/g(j) <= 1/g(ivj)+1/g(i^j).  Cleared by the
    # positive product g(i)g(j)g(ivj)g(i^j).
    recip = True
    for i in range(4):
        for j in range(4):
            J, M = g[join[i][j]], g[meet[i][j]]
            lhs = (g[i] + g[j]) * J * M
            rhs = (J + M) * g[i] * g[j]
            if (rhs - lhs).sign() < 0:
                recip = False
    claim('near_miss_g_has_supermodular_reciprocal', recip,
          '1/g(1)+1/g(0) = (g(1)+g(0))/(g(1)g(0)) = 5/5 = 1 >= 1/g(x)+1/g(y) = 1, with '
          'EQUALITY, and no pair fails; so g satisfies the source\'s necessary condition '
          'th:super-nec-cond in full')
    ls = g[3] * g[0]
    claim('near_miss_g_is_not_log_submodular', (ls - g[1] * g[2]).sign() > 0,
          'g(1)*g(0) = %s > 4 = g(x)*g(y), so g is OUTSIDE the hypothesis class of the '
          'theorem proved here' % ls)
    worst_m = None
    eq_at_claim = None
    for A in range(4):
        for B in range(4):
            for C in range(4):
                m = cleared_margin(join, meet, g, A, B, C)
                if worst_m is None or m > worst_m:
                    worst_m = m
                if (A, B, C) == CLAIM_TRIPLE_1:
                    eq_at_claim = m
    # d_g(x,y) = 1 - g(x^y)/g(xvy) = 1 - g(0)/g(1);  g(0)/g(1) = g(0)^2/(g(0)g(1)) = g(0)^2/5
    ratio = g[0] * g[0] * Q5(Fraction(1, 5))
    dxy = Q5(1) - ratio
    dx0 = Q5(1) - g[0] * Q5(Fraction(1, 2))  # 1 - g(0)/g(x) = 1 - g(0)/2
    claim('near_miss_g_satisfies_the_triangle_inequality_with_equality',
          worst_m.sign() <= 0 and eq_at_claim.sign() == 0 and dxy == CLAIM_G_DIST
          and (dx0 + dx0 - dxy) == Q5(0),
          'largest cleared margin over all 64 triples is %s (sign %d <= 0); at (x,0,y) it is '
          'exactly 0, and d_g(x,y) = %s = (sqrt5-1)/2 = d_g(x,0)+d_g(0,y) = 2*(%s). So '
          'log-submodularity is NOT necessary, the sufficient class here is STRICTLY smaller '
          'than the necessary class of th:super-nec-cond, and the source\'s open problem 4 '
          '(line 635) stays open in both directions'
          % (worst_m, worst_m.sign(), dxy, dx0))

    # ---------------- verdict ---------------------------------------------
    print('')
    note('NOT RE-RUN: the theorem of this paper quantifies over ALL lattices and all strictly '
         'positive monotone supermodular log-submodular valuations, and no program can range '
         'over that class. What is machine-checked above is (i) the two algebraic identities '
         'that carry Steps 4 and 5, as exact polynomial identities, and (ii) the order facts '
         'of Lemmas 1 and 2, on lattices of order <= 6 only. The finite sweep of section F is '
         'a CONTROL over orders 4,5,6 on a six-value integer grid; it is NOT an exhaustion of '
         'the admissible cone, and "largest margin 0" means "nothing beat the equality case", '
         'not "the cone was searched".')
    note('NOT RE-RUN: the source paper\'s own theorems (th:mod, th:super-log-sub, '
         'th:super-nec-cond) are quoted, not re-proved; and Step 3 of the proof (the sign '
         'split on p+q-b) plus the statement that d_{f,J} is in general only a PSEUDOmetric '
         'are hand arguments with no computational content. Nothing above touches the source\'s '
         'open problems 1, 3 or 4, or any literature claim.')
    print('')
    if FAILED:
        print('VERDICT: %d CHECK(S) FAILED of %d' % (FAILED, CHECKS + FAILED))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % CHECKS)
    return 0


if __name__ == '__main__':
    sys.exit(main())

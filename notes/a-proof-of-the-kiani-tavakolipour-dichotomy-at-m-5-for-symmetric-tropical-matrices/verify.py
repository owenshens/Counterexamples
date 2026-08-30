#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- re-derivation of every computational claim of

    "A Proof of the Kiani-Tavakolipour Dichotomy at m = 5 for Symmetric
     Tropical Matrices"  (paper.tex / paper.pdf in this folder)

Python 3.9+, STANDARD LIBRARY ONLY (fractions, itertools, sys).  No third-party
package, no external data file, no network, no float in any decisive step: every
comparison below is between exact integers or exact Fractions, and -infinity is
carried as the sentinel None rather than as a float.

The program reads the objects PRINTED IN THE PAPER -- the five matrices of
Section 7 and Table 1, the three sample certificates of Section 5, and the
census specification of Sections 4-5 -- and re-derives every number the paper
states about them.  It re-runs the census itself; nothing is taken on trust from
the paper or from any earlier run.

One line `PASS <name> [detail]` is printed per check, and the run closes with

    VERDICT: ALL <n> CHECKS PASS

exiting 0 if and only if every check passed.  What is NOT covered is printed by
the program itself, in the `NOT RE-RUN:` lines just before the verdict, and is
repeated in REVIEW_NOTE.md's `## Scope`.

Runtime: about 6 minutes single-threaded on a laptop (the m=5 census is ~2 min
of it and the n=5 unpooled control census ~3 min).
"""

import sys
from fractions import Fraction
from itertools import combinations, permutations

# ---------------------------------------------------------------------------
# check bookkeeping
# ---------------------------------------------------------------------------
_N_PASS = 0
_FAILED = []


def out(*a):
    print(*a)
    sys.stdout.flush()


def check(name, ok, detail=''):
    global _N_PASS
    if ok:
        _N_PASS += 1
        out('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        _FAILED.append(name)
        out('FAIL %s%s' % (name, (' ' + detail) if detail else ''))


# ---------------------------------------------------------------------------
# 0.  THE OBJECTS OF THE PAPER.  Every object below is printed in the paper
#     EXCEPT the source's own 4x4 and 9x9 examples, whose entries are
#     transcribed here from the e-print of [KT] and are not reprinted there.
#     None denotes -infinity.  Every other entry is an exact integer.
# ---------------------------------------------------------------------------
INF = None                                    # -infinity in Rmax

# Section 8, the source's own 4x4 example (nonsymmetric); entries from [KT].
PAPER_4x4 = [[0, 2, INF, 1],
             [1, 0, 3, INF],
             [INF, 2, 0, 4],
             [1, INF, 2, 0]]
PAPER_4x4_DELTAS = [0, 0, 6, 6, 10]           # as published

# Section 8, the source's own 9x9 example (nonsymmetric, all entries finite);
# entries from [KT].
PAPER_9x9 = [[-4, 2, 3, 4, 10, -1, -1, 3, -3],
             [-3, 3, 4, 3, 8, -1, 1, 6, 0],
             [-4, 1, 1, -6, 7, -1, -6, 10, -6],
             [-3, -4, -3, 9, 2, 1, -3, -2, 2],
             [-6, 2, 10, 6, -3, 5, 5, 1, -1],
             [4, 8, -5, 6, 0, -6, -1, 1, -4],
             [-2, 8, -5, -5, -4, 8, 9, 6, -3],
             [3, -2, -4, 8, -6, 3, -4, 7, 9],
             [5, -3, -4, 9, 9, 8, 10, -5, 5]]
PAPER_9x9_DELTAS = [0, 9, 18, 26, 38, 47, 56, 64, 72, 73]   # as published

# Section 8, the by-product: a NONSYMMETRIC 4x4 failing both m=4 inequalities.
BONUS_4x4 = [[1, INF, INF, INF],
             [INF, 0, 1, INF],
             [INF, INF, INF, 0],
             [INF, 0, INF, INF]]


def bonus_finite(c):
    """The same 4x4 with every -infinity replaced by the finite value c."""
    return [[c if v is INF else v for v in row] for row in BONUS_4x4]


# Section 9, the two symmetric 5x5 single-branch controls, printed there as
# C_I and C_II.  C_I fails branch (I) only; C_II fails branch (II) only.
def sym5(pairs, default):
    A = [[default] * 5 for _ in range(5)]
    for (i, j), v in pairs.items():
        A[i][j] = v
        A[j][i] = v
    return A


CTRL_I_ONLY = sym5({(0, 1): 1, (0, 2): 1, (1, 2): 1, (3, 4): 0}, -100)
CTRL_II_ONLY = sym5({(0, 1): 1, (1, 2): 1, (2, 3): 1, (3, 4): 1, (0, 4): 1,
                     (0, 0): 0, (1, 1): 0, (2, 2): 0, (3, 3): 0, (4, 4): 0}, -10)

# Section 8, "A tempting false lemma": the matrix that kills the cycle-type
# finiteness heuristic.
# 5-cycle edges 0, ALL loops and ALL chords -infinity.
C5_ALL_INF = sym5({(0, 1): 0, (1, 2): 0, (2, 3): 0, (3, 4): 0, (0, 4): 0}, INF)
for _i in range(5):
    C5_ALL_INF[_i][_i] = INF

# Section 5, the three sample certificates, transcribed from the paper.
# A cover is given as (support, image) exactly as the paper prints it.
CERTS = [
    # (m, P, Q, S, [(branch, R, lambda), ...])
    (5, ((0, 1, 2, 3, 4), (0, 1, 2, 3, 4)), ((0, 1, 2), (0, 1, 2)), ((0, 1), (0, 1)),
     [('I', ((0, 1, 2, 3), (0, 1, 2, 3)), Fraction(1, 2)),
      ('I', ((0, 1, 2, 4), (0, 1, 2, 4)), Fraction(1, 2))]),
    (5, ((0, 1, 2, 3, 4), (0, 1, 2, 3, 4)), ((0, 1, 2), (1, 0, 2)), ((0, 1), (0, 1)),
     [('I', ((0, 1, 2, 3), (1, 0, 2, 3)), Fraction(1, 2)),
      ('I', ((0, 1, 2, 4), (0, 1, 2, 4)), Fraction(1, 2))]),
    (4, ((0, 1, 2, 3), (0, 1, 2, 3)), ((0, 1), (0, 1)), ((0,), (0,)),
     [('I', ((0, 1, 2), (0, 1, 2)), Fraction(1, 2)),
      ('I', ((0, 1, 3), (0, 1, 3)), Fraction(1, 2))]),
]


# ---------------------------------------------------------------------------
# 1.  COVERS AND THEIR POSITION MULTISETS
# ---------------------------------------------------------------------------
def pos_key(i, j, sym):
    return (min(i, j), max(i, j)) if sym else (i, j)


def emult(support, image, sym):
    """E(.) as a dict position -> multiplicity."""
    d = {}
    for t in range(len(support)):
        p = pos_key(support[t], image[t], sym)
        d[p] = d.get(p, 0) + 1
    return d


def mkey(d):
    return tuple(sorted(d.items()))


def covers(verts, k, sym):
    """Every DISTINCT position multiset of a k-cover of `verts`, with one
    (support, image) description each.  A k-cover is a partial permutation whose
    support equals its image set and has size exactly k."""
    if k == 0:
        return [((), ())]
    seen = {}
    for S in combinations(verts, k):
        for pm in permutations(S):
            key = mkey(emult(S, pm, sym))
            if key not in seen:
                seen[key] = (S, pm)
    return list(seen.values())


def type_reps(m, sym):
    """One representative per relabelling class of m-covers of {0,...,m-1}."""
    seen = {}
    base = tuple(range(m))
    for pm in permutations(range(m)):
        d = emult(base, pm, sym)
        canon = None
        for s in permutations(range(m)):
            e = {}
            for (i, j), c in d.items():
                p = pos_key(s[i], s[j], sym)
                e[p] = e.get(p, 0) + c
            t = mkey(e)
            if canon is None or t < canon:
                canon = t
        if canon not in seen:
            seen[canon] = (base, pm)
    return list(seen.values())


_POOL_CACHE = {}


def pool(allowed, sym, ks):
    """For each k in `ks`, every distinct position multiset of a k-cover whose
    positions all lie in the SET `allowed`.  Multiplicities are unconstrained: a
    2-cycle on {i,j} uses the single position (i,j) twice.

    Enumerated over the vertices spanned by `allowed`, which is complete: a
    cover with positions inside `allowed` has its support inside those
    vertices."""
    ck = (tuple(sorted(allowed)), sym, tuple(ks))
    if ck in _POOL_CACHE:
        return _POOL_CACHE[ck]
    A = set(allowed)
    verts = sorted({v for p in A for v in p})
    res = {}
    for k in ks:
        acc = {}
        for S in combinations(verts, k):
            for pm in permutations(S):
                d = emult(S, pm, sym)
                if all(p in A for p in d):
                    acc[mkey(d)] = True
        res[k] = list(acc.keys())
    _POOL_CACHE[ck] = res
    return res


# ---------------------------------------------------------------------------
# 2.  DELTA_k BY EXACT BRUTE FORCE OVER PRINCIPAL SUBMATRICES
# ---------------------------------------------------------------------------
def deltas(A, upto):
    """delta_0..delta_upto with delta_0 = 0.  None means -infinity.  Exact
    integer arithmetic: no float anywhere."""
    n = len(A)
    d = [0]
    for k in range(1, upto + 1):
        best = None
        for S in combinations(range(n), k):
            for pm in permutations(S):
                s = 0
                ok = True
                for t in range(k):
                    a = A[S[t]][pm[t]]
                    if a is None:
                        ok = False
                        break
                    s += a
                if ok and (best is None or s > best):
                    best = s
        d.append(best)
    return d


def gt(x, y):
    """x > y in Rmax, with None = -infinity."""
    if x is None:
        return False
    if y is None:
        return True
    return x > y


def add(*xs):
    s = 0
    for x in xs:
        if x is None:
            return None
        s += x
    return s


def smul(c, x):
    return None if x is None else c * x


def fails_I(d, m):
    """delta_m + delta_{m-2} > 2 delta_{m-1}  (branch (I) of the paper, failing)."""
    return gt(add(d[m], d[m - 2]), smul(2, d[m - 1]))


def fails_II(d, m):
    """delta_m + 2 delta_{m-3} > 3 delta_{m-2}  (branch (II), failing)."""
    return gt(add(d[m], smul(2, d[m - 3])), smul(3, d[m - 2]))


# ---------------------------------------------------------------------------
# 3.  EXACT PHASE-1 SIMPLEX (Bland's rule; terminates, no cycling)
# ---------------------------------------------------------------------------
class LPBudget(Exception):
    pass


ITER_CAP = 200000


def solve_eq_nonneg(rows, rhs):
    """Find x >= 0 with rows . x = rhs, exactly, or return None if none exists.
    rows is a list of coefficient rows (integers or Fractions)."""
    m = len(rows)
    if m == 0:
        return None
    N = len(rows[0])
    if N == 0:
        return None
    T = [[Fraction(v) for v in rows[i]] for i in range(m)]
    b = [Fraction(v) for v in rhs]
    for i in range(m):
        if b[i] < 0:
            T[i] = [-v for v in T[i]]
            b[i] = -b[i]
    tot = N + m
    for i in range(m):
        T[i] = T[i] + [Fraction(1 if t == i else 0) for t in range(m)] + [b[i]]
    basis = [N + i for i in range(m)]
    inb = [False] * tot
    for v in basis:
        inb[v] = True
    obj = [Fraction(0)] * (tot + 1)
    for i in range(m):
        for j in range(tot + 1):
            obj[j] += T[i][j]
    for i in range(m):
        obj[N + i] -= 1
    it = 0
    while True:
        it += 1
        if it > ITER_CAP:
            raise LPBudget('simplex exceeded %d iterations' % ITER_CAP)
        piv = -1
        for j in range(tot):
            if obj[j] > 0 and not inb[j]:
                piv = j
                break
        if piv < 0:
            break
        br = -1
        best = None
        for i in range(m):
            if T[i][piv] > 0:
                r = T[i][tot] / T[i][piv]
                if best is None or r < best or (r == best and basis[i] < basis[br]):
                    best = r
                    br = i
        if br < 0:
            return None
        pv = T[br][piv]
        T[br] = [v / pv for v in T[br]]
        Tb = T[br]
        for i in range(m):
            if i != br and T[i][piv] != 0:
                f = T[i][piv]
                Ti = T[i]
                T[i] = [Ti[j] - f * Tb[j] for j in range(tot + 1)]
        if obj[piv] != 0:
            f = obj[piv]
            obj = [obj[j] - f * Tb[j] for j in range(tot + 1)]
        inb[basis[br]] = False
        inb[piv] = True
        basis[br] = piv
    if obj[tot] != 0:
        return None
    x = [Fraction(0)] * N
    for i in range(m):
        if basis[i] < N:
            x[basis[i]] = T[i][tot]
    return x


def gordan_dead(forms, D):
    """A certificate that {f(x) > 0 for all f in forms} has NO real solution:
    lambda >= 0 with sum(lambda) = 1 and sum lambda_j f_j = 0.  Returns the
    verified lambda, or None if no such lambda exists.

    Only this direction is used in the proof and it is elementary: if such a
    lambda exists then 0 = sum_j lambda_j f_j(x) > 0 for any strict x, absurd.
    Gordan's theorem supplies the converse, which is not on the proof path."""
    if not forms:
        return None
    N = len(forms)
    rows = [[forms[j][t] for j in range(N)] for t in range(D)]
    rows.append([1] * N)
    rhs = [0] * D + [1]
    lam = solve_eq_nonneg(rows, rhs)
    if lam is None:
        return None
    # re-verify, exactly, without reference to the solver's bookkeeping
    if any(v < 0 for v in lam):
        return None
    if sum(lam) != 1:
        return None
    for t in range(D):
        if sum(lam[j] * forms[j][t] for j in range(N)) != 0:
            return None
    return lam


# ---------------------------------------------------------------------------
# 4.  THE CENSUS
# ---------------------------------------------------------------------------
def build_system(m, mP, mQ, mS, sym, use_I=True, use_II=True, allowed=None):
    """The pooled constraint forms for one triple.  Returns (forms, tags, D)."""
    if allowed is None:
        allowed = sorted(set(mP) | set(mQ) | set(mS))
    idx = {p: t for t, p in enumerate(allowed)}
    D = len(allowed)

    def vec(d):
        v = [0] * D
        for p, c in d.items():
            v[idx[p]] += c
        return v

    vP, vQ, vS = vec(mP), vec(mQ), vec(mS)
    pl = pool(allowed, sym, (m - 1, m - 2))
    forms = []
    tags = []
    if use_I:
        for key in pl[m - 1]:
            vr = vec(dict(key))
            forms.append(tuple(vP[t] + vQ[t] - 2 * vr[t] for t in range(D)))
            tags.append(('I', key))
    if use_II:
        for key in pl[m - 2]:
            vq = vec(dict(key))
            forms.append(tuple(vP[t] + 2 * vS[t] - 3 * vq[t] for t in range(D)))
            tags.append(('II', key))
    return forms, tags, D, allowed


def simple_certificate(m, mP, mQ, mS, sym, allowed):
    """Look first for a certificate a reader can check by counting symbols:
      (a)  E(P) + E(Q)   = E(R1) + E(R2)          (branch (I), lambda = 1/2, 1/2)
      (b)  E(P) + 2 E(S) = E(Q1)+E(Q2)+E(Q3)      (branch (II), lambda = 1/3 x3)
    Returns ('I2'|'II3', lambda list, tag list) or None."""
    idx = {p: t for t, p in enumerate(allowed)}
    D = len(allowed)

    def vec(d):
        v = [0] * D
        for p, c in d.items():
            v[idx[p]] += c
        return tuple(v)

    vP, vQ, vS = vec(mP), vec(mQ), vec(mS)
    pl = pool(allowed, sym, (m - 1, m - 2))
    R = {}
    for key in pl[m - 1]:
        R[vec(dict(key))] = key
    tgt = tuple(vP[t] + vQ[t] for t in range(D))
    for a, ka in R.items():
        need = tuple(tgt[t] - a[t] for t in range(D))
        if min(need) >= 0 and need in R:
            return 'I2', [Fraction(1, 2), Fraction(1, 2)], [('I', ka), ('I', R[need])]
    Q3 = {}
    for key in pl[m - 2]:
        Q3[vec(dict(key))] = key
    tgt2 = tuple(vP[t] + 2 * vS[t] for t in range(D))
    for a, ka in Q3.items():
        r1 = tuple(tgt2[t] - a[t] for t in range(D))
        if min(r1) < 0:
            continue
        for b, kb in Q3.items():
            r2 = tuple(r1[t] - b[t] for t in range(D))
            if min(r2) >= 0 and r2 in Q3:
                return ('II3', [Fraction(1, 3)] * 3,
                        [('II', ka), ('II', kb), ('II', Q3[r2])])
    return None


def verify_multipliers(m, mP, mQ, mS, sym, allowed, lam, tags):
    """Re-substitute an arbitrary multiplier vector into the forms and confirm,
    exactly, that the nonnegative combination is identically zero."""
    idx = {p: t for t, p in enumerate(allowed)}
    D = len(allowed)

    def vec(d):
        v = [0] * D
        for p, c in d.items():
            v[idx[p]] += c
        return v

    vP, vQ, vS = vec(mP), vec(mQ), vec(mS)
    if any(v < 0 for v in lam) or sum(lam) != 1:
        return False
    acc = [Fraction(0)] * D
    for w, (branch, key) in zip(lam, tags):
        vr = vec(dict(key))
        if branch == 'I':
            f = [vP[t] + vQ[t] - 2 * vr[t] for t in range(D)]
        else:
            f = [vP[t] + 2 * vS[t] - 3 * vr[t] for t in range(D)]
        for t in range(D):
            acc[t] += w * f[t]
    return all(v == 0 for v in acc)


def norm_key(mP, mQ, mS, sym):
    """Relabel the index set order-preservingly.  Two triples with the same key
    have IDENTICAL constraint systems, so one may be solved for both.  Sound
    because relabelling the index set carries the pool bijectively onto the
    relabelled pool and permutes the coordinates of every form."""
    verts = sorted({v for d in (mP, mQ, mS) for p in d for v in p})
    r = {v: i for i, v in enumerate(verts)}

    def rk(p):
        return pos_key(r[p[0]], r[p[1]], sym)

    return tuple(tuple(sorted((rk(p), c) for p, c in d.items()))
                 for d in (mP, mQ, mS))


def census(m, sym, dedup=True):
    """The all-n census in the window V0 = [3m-5]."""
    V = list(range(m + (m - 2) + (m - 3)))
    Preps = type_reps(m, sym)
    Qs = covers(V, m - 2, sym)
    Ss = covers(V, m - 3, sym)
    total = len(Preps) * len(Qs) * len(Ss)
    classes = {}
    order = []
    for dP in Preps:
        mP = emult(dP[0], dP[1], sym)
        for dQ in Qs:
            mQ = emult(dQ[0], dQ[1], sym)
            for dS in Ss:
                mS = emult(dS[0], dS[1], sym)
                k = norm_key(mP, mQ, mS, sym) if dedup else (len(order),)
                if k not in classes:
                    classes[k] = 1
                    order.append((mP, mQ, mS))
                else:
                    classes[k] += 1
    dead = 0
    live = []
    shape = {'I2': 0, 'II3': 0, 'LP': 0}
    maxden = 1
    reverified = 0
    for (mP, mQ, mS) in order:
        allowed = sorted(set(mP) | set(mQ) | set(mS))
        sc = simple_certificate(m, mP, mQ, mS, sym, allowed)
        if sc is not None:
            kind, lam, tags = sc
        else:
            forms, tags, D, allowed = build_system(m, mP, mQ, mS, sym, allowed=allowed)
            lam = gordan_dead(forms, D)
            if lam is None:
                live.append((mP, mQ, mS))
                continue
            kind = 'LP'
            tags = [tags[j] for j in range(len(lam)) if lam[j] != 0]
            lam = [v for v in lam if v != 0]
        if not verify_multipliers(m, mP, mQ, mS, sym, allowed, lam, tags):
            live.append((mP, mQ, mS))
            continue
        reverified += 1
        shape[kind] += 1
        for v in lam:
            maxden = max(maxden, v.denominator)
        dead += 1
    return dict(V=len(V), nP=len(Preps), nQ=len(Qs), nS=len(Ss), total=total,
                nclass=len(order), dead=dead, live=live, shape=shape,
                maxden=maxden, reverified=reverified,
                covered=sum(classes[k] for k in classes))


def census_unpooled_n5():
    """The independent control: at n = 5, m = 5, use the FULL constraint set --
    every 4-cover and every 3-cover of [5], not just those inside Allowed --
    over the 15 symmetric positions."""
    n = 5
    allpos = [(i, j) for i in range(n) for j in range(i, n)]
    idx = {p: t for t, p in enumerate(allpos)}
    D = len(allpos)

    def vec(d):
        v = [0] * D
        for p, c in d.items():
            v[idx[p]] += c
        return v

    V = list(range(n))
    C5 = covers(V, 5, True)
    C4 = covers(V, 4, True)
    C3 = covers(V, 3, True)
    C2 = covers(V, 2, True)
    Preps = type_reps(5, True)
    v4 = [vec(emult(d[0], d[1], True)) for d in C4]
    v3 = [vec(emult(d[0], d[1], True)) for d in C3]
    dead = 0
    live = 0
    for dP in Preps:
        vP = vec(emult(dP[0], dP[1], True))
        for dQ in C3:
            vQ = vec(emult(dQ[0], dQ[1], True))
            for dS in C2:
                vS = vec(emult(dS[0], dS[1], True))
                forms = [tuple(vP[t] + vQ[t] - 2 * r[t] for t in range(D)) for r in v4]
                forms += [tuple(vP[t] + 2 * vS[t] - 3 * q[t] for t in range(D)) for q in v3]
                lam = gordan_dead(forms, D)
                if lam is None:
                    live += 1
                else:
                    dead += 1
    return dict(D=D, nC5=len(C5), nC4=len(C4), nC3=len(C3), nC2=len(C2),
                nP=len(Preps), forms=len(v4) + len(v3),
                total=len(Preps) * len(C3) * len(C2), dead=dead, live=live)


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    out('verify.py -- "The Kiani-Tavakolipour Dichotomy at m = 5 for Symmetric')
    out('              Tropical Matrices, in Additive Form"')
    out('exact integer / Fraction arithmetic only; -infinity carried as None')
    out('')

    # -------------------------------------------------------------------
    out('--- 1.  the census specification of Sections 4 and 5 ---')
    # -------------------------------------------------------------------
    check('window-size-3m-5', 3 * 5 - 5 == 10 and 5 + 3 + 2 == 10,
          '|Allowed| <= m + (m-2) + (m-3) = 3m-5 = 10 at m = 5')

    P5 = type_reps(5, True)
    check('P-classes-is-7', len(P5) == 7,
          'relabelling classes of a 5-cover = %d = p(5), the partitions of 5' % len(P5))
    cycle_types = sorted(tuple(sorted(len(c) for c in _cycles(d[1])))
                         for d in P5)
    parts5 = sorted(_partitions(5))
    check('P-classes-are-the-partitions-of-5', cycle_types == parts5,
          '%s' % (['+'.join(str(x) for x in t) for t in cycle_types],))

    V10 = list(range(10))
    C3_10 = covers(V10, 3, True)
    C2_10 = covers(V10, 2, True)
    check('C3-of-10-is-600', len(C3_10) == 600,
          '|C_3([10])| = %d  (120 triples of loops + 360 two-cycle-plus-loop + 120 triangles)'
          % len(C3_10))
    check('C2-of-10-is-90', len(C2_10) == 90,
          '|C_2([10])| = %d  (45 pairs of loops + 45 two-cycles)' % len(C2_10))
    check('census-size-378000', 7 * len(C3_10) * len(C2_10) == 378000,
          '7 x %d x %d = %d triples' % (len(C3_10), len(C2_10), 7 * len(C3_10) * len(C2_10)))

    V5 = list(range(5))
    n5 = {k: len(covers(V5, k, True)) for k in (2, 3, 4, 5)}
    check('C5-of-5-is-73', n5[5] == 73, '|C_5([5])| = %d' % n5[5])
    check('C4-of-5-is-85', n5[4] == 85, '|C_4([5])| = %d = 5 x 17' % n5[4])
    check('C3-of-5-is-50', n5[3] == 50, '|C_3([5])| = %d = 10 x 5' % n5[3])
    check('C2-of-5-is-20', n5[2] == 20, '|C_2([5])| = %d = 10 x 2' % n5[2])
    check('forms-per-LP-at-n5-is-135', n5[4] + n5[3] == 135,
          '%d + %d = 135 DISTINCT forms, not 180: for symmetric A permutations '
          'sharing a position multiset give the same form' % (n5[4], n5[3]))
    check('unpooled-n5-census-is-7000', 7 * n5[3] * n5[2] == 7000,
          '7 x %d x %d = %d distinct LPs, not 8400' % (n5[3], n5[2], 7 * n5[3] * n5[2]))

    # delta_4 is finite whenever both m=5 inequalities fail: every 5-cover
    # contains a 4-cover on its OWN positions.  Checked on all 7 classes.
    ok4 = True
    wit = []
    for d in P5:
        mP = emult(d[0], d[1], True)
        pl = pool(sorted(set(mP)), True, (4,))
        if not pl[4]:
            ok4 = False
        else:
            wit.append(len(pl[4]))
    check('every-5-cover-contains-a-4-cover-on-its-own-positions', ok4,
          'all 7 classes; counts of such 4-covers = %s -- this is why delta_4 '
          'is finite once delta_5 is' % (wit,))

    check('partial-perms-of-principal-4x4-is-65',
          sum(len(list(permutations(range(k)))) * len(list(combinations(range(4), k)))
              for k in range(5)) == 65,
          'sum_k C(4,k) k! = 1+4+12+24+24 = 65 (NOT 209 = sum_k C(4,k)^2 k!, '
          'which chooses row and column sets independently)')

    # -------------------------------------------------------------------
    out('')
    out('--- 2.  the matrices printed in the paper, delta_k by exact brute force ---')
    # -------------------------------------------------------------------
    d = deltas(PAPER_4x4, 4)
    check('source-4x4-delta-sequence', d == PAPER_4x4_DELTAS,
          'delta_0..delta_4 = %s, as published' % (d,))
    check('source-4x4-is-not-a-both-fail-object',
          fails_I(d, 4) and not fails_II(d, 4),
          '(I) 10+6 = 16 > 12 = 2x6 FAILS; (II) 10+2x0 = 10 <= 18 = 3x6 HOLDS')

    d9 = deltas(PAPER_9x9, 9)
    check('source-9x9-delta-sequence', d9 == PAPER_9x9_DELTAS,
          'delta_0..delta_9 = %s, all ten as published' % (d9,))
    check('source-9x9-both-fail-at-m4', fails_I(d9, 4) and fails_II(d9, 4),
          '(I) 38+18 = 56 > 52 = 2x26 FAILS; (II) 38+2x9 = 56 > 54 = 3x18 FAILS')
    check('source-9x9-branch-I-holds-at-m5', not fails_I(d9, 5),
          'delta_5+delta_3 = 47+26 = 73 <= 76 = 2x38 = 2 delta_4 (the right side '
          'is 2 delta_4, not 2 delta_5)')

    db = deltas(BONUS_4x4, 4)
    check('byproduct-4x4-delta-sequence', db == [0, 1, 1, 1, 2],
          'delta_0..delta_4 = %s' % (db,))
    check('byproduct-4x4-both-fail-at-m4', fails_I(db, 4) and fails_II(db, 4),
          'delta = %s, so at m = 4 branch (I) delta_4+delta_2 <= 2 delta_3 %s '
          'and branch (II) delta_4+2 delta_1 <= 3 delta_2 %s'
          % (db, 'FAILS' if fails_I(db, 4) else 'HOLDS',
             'FAILS' if fails_II(db, 4) else 'HOLDS'))
    for c in (-20, -50, -100):
        dc = deltas(bonus_finite(c), 4)
        check('byproduct-4x4-all-finite-at-%d' % c,
              dc == [0, 1, 1, 1, 2] and fails_I(dc, 4) and fails_II(dc, 4),
              'every -infinity replaced by %d: delta = %s, both branches still fail' % (c, dc))
    check('byproduct-4x4-is-nonsymmetric',
          any(BONUS_4x4[i][j] != BONUS_4x4[j][i] for i in range(4) for j in range(4)),
          'a_{1,3} = -infinity while a_{3,1} = 0, so the object does not '
          'contradict the theorem of this paper')
    check('m4-violation-forces-n-at-least-4', len(BONUS_4x4) == 4,
          'the conjecture ranges over 3 <= m <= n, so an m = 4 violation needs '
          'n >= 4; this object attains n = 4')

    dA = deltas(CTRL_I_ONLY, 5)
    check('control-fails-branch-I-only',
          dA == [0, -100, 2, 3, 2, 3] and fails_I(dA, 5) and not fails_II(dA, 5),
          'symmetric 5x5 C_I, delta_0..delta_5 = %s; branch (I) %s, branch (II) %s'
          % (dA, 'FAILS' if fails_I(dA, 5) else 'HOLDS',
             'FAILS' if fails_II(dA, 5) else 'HOLDS'))
    dB = deltas(CTRL_II_ONLY, 5)
    check('control-fails-branch-II-only',
          dB == [0, 0, 2, 2, 4, 5] and not fails_I(dB, 5) and fails_II(dB, 5),
          'symmetric 5x5 C_II, delta_0..delta_5 = %s; branch (I) %s, branch (II) %s'
          % (dB, 'FAILS' if fails_I(dB, 5) else 'HOLDS',
             'FAILS' if fails_II(dB, 5) else 'HOLDS'))

    dC = deltas(C5_ALL_INF, 5)
    check('minus-infinity-does-not-walk-down-the-cycle-types',
          dC == [0, None, 0, None, 0, 0],
          'symmetric 5-cycle with all loops and all chords -infinity: '
          'delta_2 = delta_4 = delta_5 = 0 are FINITE while delta_3 = -infinity, '
          'so "delta_5 finite forces delta_3 finite" is false')
    subs = []
    for pm in permutations((0, 1, 2)):
        w = 0
        okp = True
        for t, v in enumerate((0, 1, 2)):
            a = C5_ALL_INF[v][pm[t]]
            if a is None:
                okp = False
                break
            w += a
        subs.append(w if okp else None)
    check('minus-infinity-witness-3x3-block-is-all-minus-infinity',
          all(v is None for v in subs),
          'all 6 permutations of the principal 3x3 on {0,1,2} are -infinity')
    check('minus-infinity-witness-still-satisfies-the-dichotomy',
          not fails_I(dC, 5),
          '(I) reads -infinity <= 0 and holds trivially, so the theorem is untouched')

    # -------------------------------------------------------------------
    out('')
    out('--- 3.  the three sample certificates printed in Section 5 ---')
    # -------------------------------------------------------------------
    for ci, (m, dP, dQ, dS, terms) in enumerate(CERTS, 1):
        mP = emult(dP[0], dP[1], True)
        mQ = emult(dQ[0], dQ[1], True)
        mS = emult(dS[0], dS[1], True)
        allowed = sorted(set(mP) | set(mQ) | set(mS))
        tags = []
        lam = []
        legal = True
        for branch, dR, w in terms:
            mR = emult(dR[0], dR[1], True)
            k = m - 1 if branch == 'I' else m - 2
            if len(dR[0]) != k or not all(p in set(allowed) for p in mR):
                legal = False
            tags.append((branch, mkey(mR)))
            lam.append(w)
        okc = legal and verify_multipliers(m, mP, mQ, mS, True, allowed, lam, tags)
        lhs = dict(mP)
        for p, c in mQ.items():
            lhs[p] = lhs.get(p, 0) + c
        rhs = {}
        for branch, dR, w in terms:
            for p, c in emult(dR[0], dR[1], True).items():
                rhs[p] = rhs.get(p, 0) + c
        check('sample-certificate-%d' % ci, okc and lhs == rhs,
              'm=%d: E(P)+E(Q) = %s = %s = sum of the two 4-covers, multipliers '
              '%s (nonnegative, summing to 1)'
              % (m, _ms(lhs), _ms(rhs), [str(v) for v in lam]))

    # -------------------------------------------------------------------
    out('')
    out('--- 4.  anti-controls: the engine must NOT certify a live system ---')
    # -------------------------------------------------------------------
    for lbl, A, uI, uII in (('branch-I-only', CTRL_I_ONLY, True, False),
                            ('branch-II-only', CTRL_II_ONLY, False, True)):
        dd = deltas(A, 5)
        dP = _argmax_cover(A, 5)
        dQ = _argmax_cover(A, 3)
        dS = _argmax_cover(A, 2)
        mP = emult(dP[0], dP[1], True)
        mQ = emult(dQ[0], dQ[1], True)
        mS = emult(dS[0], dS[1], True)
        forms, tags, D, allowed = build_system(5, mP, mQ, mS, True, uI, uII)
        lam = gordan_dead(forms, D)
        check('anti-control-%s-is-live' % lbl, lam is None,
              '%d pooled constraints; the matrix itself is a strict solution, so '
              'no nonnegative certificate may exist -- engine returned %s'
              % (len(forms), 'none (correct)' if lam is None else 'ONE (UNSOUND)'))
    triv = gordan_dead([(1, 0), (0, 1)], 2)
    check('anti-control-trivially-live-system', triv is None,
          'forms x>0, y>0 admit a strict solution and the engine finds no certificate')
    triv2 = gordan_dead([(1, 0), (-1, 0)], 2)
    check('positive-control-trivially-dead-system',
          triv2 is not None and sum(triv2) == 1,
          'forms x>0, -x>0 are contradictory and the engine returns lambda = %s'
          % ([str(v) for v in triv2] if triv2 else None))

    # -------------------------------------------------------------------
    out('')
    out('--- 5.  the censuses (all n at once, in the window V0 = [3m-5]) ---')
    # -------------------------------------------------------------------
    r3 = census(3, True)
    check('census-m3-symmetric-all-n-silent',
          r3['live'] == [] and r3['dead'] == r3['nclass'] and r3['total'] == 12,
          'V0=[%d]  %d x %d x %d = %d triples in %d classes, dead %d, live 0 -- '
          'the source theorem for m=3 re-proved for every n by this route'
          % (r3['V'], r3['nP'], r3['nQ'], r3['nS'], r3['total'], r3['nclass'], r3['dead']))

    r4 = census(4, True)
    check('census-m4-symmetric-all-n-silent',
          r4['live'] == [] and r4['dead'] == r4['nclass'] and r4['total'] == 1470,
          'V0=[%d]  %d x %d x %d = %d triples in %d classes, dead %d, live 0 -- '
          'the source theorem for m=4 re-proved for every n'
          % (r4['V'], r4['nP'], r4['nQ'], r4['nS'], r4['total'], r4['nclass'], r4['dead']))

    r4n = census(4, False)
    check('forced-positive-census-m4-nonsymmetric-is-not-silent',
          len(r4n['live']) > 0 and r4n['total'] == 1470,
          'V0=[%d]  %d triples in %d classes, dead %d, LIVE %d > 0 -- symmetry is '
          'load-bearing and the engine detects it rather than certifying '
          'everything put to it; the by-product 4x4 checked above is an explicit '
          'both-fail object of this kind'
          % (r4n['V'], r4n['total'], r4n['nclass'], r4n['dead'], len(r4n['live'])))

    r5 = census(5, True)
    check('census-m5-symmetric-all-n-exhaustive',
          r5['total'] == 378000 and r5['covered'] == 378000,
          '%d triples, every one assigned to one of %d relabelling classes '
          '(%d + 0 = %d)' % (r5['total'], r5['nclass'], r5['covered'], r5['total']))
    check('census-m5-symmetric-all-n-silent',
          r5['live'] == [] and r5['dead'] == r5['nclass'],
          'dead %d of %d classes, LIVE 0 -- so no symmetric matrix over Rmax, at '
          'any n >= 5, fails both m=5 inequalities'
          % (r5['dead'], r5['nclass']))
    check('census-m5-multipliers-re-substituted',
          r5['reverified'] == r5['nclass'],
          '%d of %d certificates re-substituted into their own forms and found to '
          'sum to the zero form exactly, with lambda >= 0 and sum lambda = 1'
          % (r5['reverified'], r5['nclass']))
    check('census-m5-certificate-shapes',
          (r5['shape']['I2'], r5['shape']['II3'], r5['shape']['LP'])
          == (38524, 1016, 14360)
          and sum(r5['shape'].values()) == r5['nclass'],
          '%d classes are killed by the two-term identity E(P)+E(Q) = E(R1)+E(R2), '
          '%d by the three-term identity E(P)+2E(S) = E(Q1)+E(Q2)+E(Q3), and %d '
          'need a general nonnegative rational combination (38524 + 1016 + 14360 '
          '= 53900, as stated in Section 5 of the paper)'
          % (r5['shape']['I2'], r5['shape']['II3'], r5['shape']['LP']))
    check('census-m5-denominators-are-small', r5['maxden'] == 48,
          'largest denominator in any multiplier = %d, as stated in Section 5, so '
          'every certificate is a small rational vector' % r5['maxden'])

    ru = census_unpooled_n5()
    check('control-census-n5-unpooled-full-constraint-set',
          ru['live'] == 0 and ru['dead'] == ru['total'] == 7000,
          'n=5, D=%d symmetric positions, %d constraints per LP (%d 4-covers + %d '
          '3-covers, no pooling), %d x %d x %d = %d LPs, dead %d, live 0'
          % (ru['D'], ru['forms'], ru['nC4'], ru['nC3'], ru['nP'], ru['nC3'],
             ru['nC2'], ru['total'], ru['dead']))

    # -------------------------------------------------------------------
    out('')
    out('NOT RE-RUN: m >= 6.  This program checks m = 5 (and re-proves the '
        'published m = 3 and m = 4 cases); the conjecture ranges over '
        '3 <= m <= n and m >= 6 is untouched here.  For m = 6 the same reduction '
        'gives the window V0 = [13] and a raw triple count of about 1.9 x 10^8, '
        'which this program makes no attempt at.')
    out('NOT RE-RUN: the 378,000 individual multiplier vectors are not shipped '
        'as data.  They are regenerated and re-verified by this run, but no file '
        'in this folder lists them, so a referee wanting the full list must run '
        'this program.')
    out('NOT RE-RUN: the converse of the reduction.  Only the forward direction '
        'is checked and only it is used -- a pooled certificate is a certificate '
        'at every n.  The claim that a LIVE pooled triple always lifts to a real '
        'counterexample is not verified here and no conclusion of the paper rests '
        'on it.')
    out('NOT RE-RUN: no prior-art or literature search is performed by this '
        'program; it settles mathematics, not novelty.')
    out('')

    if _FAILED:
        out('VERDICT: %d CHECKS FAILED: %s' % (len(_FAILED), ', '.join(_FAILED)))
        return 1
    out('VERDICT: ALL %d CHECKS PASS' % _N_PASS)
    return 0


# ---------------------------------------------------------------------------
# small helpers used above
# ---------------------------------------------------------------------------
def _cycles(pm):
    n = len(pm)
    seen = [False] * n
    cyc = []
    for i in range(n):
        if seen[i]:
            continue
        c = []
        j = i
        while not seen[j]:
            seen[j] = True
            c.append(j)
            j = pm[j]
        cyc.append(tuple(c))
    return cyc


def _partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in _partitions(n - k, k):
            yield tuple(sorted((k,) + rest))


def _ms(d):
    parts = []
    for p in sorted(d):
        parts += ['%d%d' % p] * d[p]
    return '{' + ','.join(parts) + '}'


def _argmax_cover(A, k):
    n = len(A)
    best = None
    arg = None
    for S in combinations(range(n), k):
        for pm in permutations(S):
            s = 0
            ok = True
            for t in range(k):
                a = A[S[t]][pm[t]]
                if a is None:
                    ok = False
                    break
                s += a
            if ok and (best is None or s > best):
                best = s
                arg = (S, pm)
    return arg


if __name__ == '__main__':
    sys.exit(main())

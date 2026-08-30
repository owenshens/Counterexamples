#!/usr/bin/env python3
"""verify.py -- re-derives every quantity claimed in

    "A proof of the Chauhan-Shukla-Vinayak conjecture on the 3-cut complex of K_m box P_n"

from the objects PRINTED IN THE PAPER and from nothing else.

The paper prints two things, and they are the program's only input:

  (i)  the host graph, Section 2: G = K_m box P_n has vertex set {0,...,m-1} x {0,...,n-1}, and
       (i,j) ~ (i',j') iff (j = j' and i /= i') or (i = i' and |j - j'| = 1);
  (ii) the description T0 of the Alexander dual, Theorem 3: the maximal faces of dimension >= 2
       of D* = Delta_3(G)* are exactly the n column simplices, the (n-1)C(m,2) squares and the
       m(n-2) row triangles, every other maximal face being an edge.

Everything else -- the two closed forms A and B, the reduced homology of the two-column block, of
Y and of D*, the Alexander-duality bookkeeping, the ten cells of the source's table, the two
boundary lines already in print, and the source's own k = 2 theorem -- is COMPUTED here.

Contract: Python 3.9+, standard library only, no external data file, no network, no randomness, no
argument, no seed. All arithmetic is exact integer arithmetic; homology is computed by an integer
Smith normal form, so torsion is DETECTED and not assumed away. No floating-point value decides
any check. One `PASS <name> [detail]` line per check, then a single VERDICT line; exit 0 iff every
check passed.
"""

import sys
from itertools import combinations

# ---------------------------------------------------------------------------
# 0. the check harness
# ---------------------------------------------------------------------------
_PASS = 0
_FAIL = 0


def check(name, ok, detail=''):
    global _PASS, _FAIL
    if ok:
        _PASS += 1
        print('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        _FAIL += 1
        print('FAIL %s%s' % (name, (' ' + detail) if detail else ''))
    sys.stdout.flush()
    return bool(ok)


def head(title):
    print()
    print('=== %s ===' % title)
    sys.stdout.flush()


# ---------------------------------------------------------------------------
# 1. exact integer linear algebra: Smith normal form
# ---------------------------------------------------------------------------
def _gcd(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def _normalise_diagonal(diag):
    """The Smith divisors of a DIAGONAL integer matrix: repeatedly replace a non-dividing pair
    (a, b) by (gcd(a, b), lcm(a, b)). This terminates: the product of the entries is invariant and
    every replacement strictly decreases the smaller member of the pair, which cannot go on for
    ever in the positive integers."""
    d = sorted(abs(x) for x in diag if x)
    changed = True
    while changed:
        changed = False
        for i in range(len(d)):
            for j in range(i + 1, len(d)):
                a, b = d[i], d[j]
                if b % a:
                    g = _gcd(a, b)
                    d[i], d[j] = g, a // g * b
                    changed = True
        d.sort()
    return d


def _dense_snf(rows):
    """(rank, [Smith divisors]) of a small dense integer matrix given as a list of row lists."""
    M = [list(r) for r in rows]
    nr = len(M)
    nc = len(M[0]) if nr else 0
    diag = []
    t = 0
    while t < nr and t < nc:
        # pivot: nonzero entry of least absolute value in the remaining submatrix
        best = None
        for i in range(t, nr):
            Mi = M[i]
            for j in range(t, nc):
                v = Mi[j]
                if v and (best is None or abs(v) < best[0]):
                    best = (abs(v), i, j)
        if best is None:
            break
        _, pi, pj = best
        M[t], M[pi] = M[pi], M[t]
        for i in range(nr):
            M[i][t], M[i][pj] = M[i][pj], M[i][t]
        # clear the pivot row and column; repeat until both are clear
        while True:
            p = M[t][t]
            dirty = False
            for i in range(t + 1, nr):
                if M[i][t]:
                    q = M[i][t] // p
                    if q:
                        for j in range(t, nc):
                            M[i][j] -= q * M[t][j]
                    if M[i][t]:                       # remainder: swap in a smaller pivot
                        M[t], M[i] = M[i], M[t]
                        dirty = True
            for j in range(t + 1, nc):
                if M[t][j]:
                    q = M[t][j] // p
                    if q:
                        for i in range(t, nr):
                            M[i][j] -= q * M[i][t]
                    if M[t][j]:
                        for i in range(nr):
                            M[i][t], M[i][j] = M[i][j], M[i][t]
                        dirty = True
            if not dirty and not any(M[i][t] for i in range(t + 1, nr)) \
                         and not any(M[t][j] for j in range(t + 1, nc)):
                break
        diag.append(M[t][t])
        t += 1
    return len(diag), _normalise_diagonal(diag)


def smith(columns, nrows):
    """(rank, [divisors > 1]) of the integer matrix with `nrows` rows whose columns are the sparse
    dicts {row_index: value}.  Strategy: a boundary matrix is full of +-1 entries, so peel unit
    pivots off first -- each contributes an elementary divisor 1 and cannot create torsion -- then
    run a dense Smith normal form on whatever remains.  Exact integers throughout."""
    cols = {c: dict(d) for c, d in enumerate(columns) if d}
    rows = {}
    for c, d in cols.items():
        for r in d:
            assert 0 <= r < nrows, 'row index %d outside [0, %d)' % (r, nrows)
            rows.setdefault(r, set()).add(c)
    rank = 0
    while True:
        pivot = None
        for c, d in cols.items():
            for r, v in d.items():
                if v in (1, -1):
                    w = (len(d), len(rows[r]))
                    if pivot is None or w < pivot[0]:
                        pivot = (w, c, r)
            if pivot is not None and pivot[0][0] == 1:
                break
        if pivot is None:
            break
        _, pc, pr = pivot
        a = cols[pc][pr]
        for c2 in list(rows[pr]):
            if c2 == pc:
                continue
            f = cols[c2][pr] // a
            d2 = cols[c2]
            for rr, vv in cols[pc].items():
                nv = d2.get(rr, 0) - f * vv
                if nv:
                    if rr not in d2:
                        rows.setdefault(rr, set()).add(c2)
                    d2[rr] = nv
                elif rr in d2:
                    del d2[rr]
                    rows[rr].discard(c2)
        for rr in cols[pc]:
            rows[rr].discard(pc)
        del cols[pc]
        for c2 in list(rows[pr]):
            cols[c2].pop(pr, None)
        rows.pop(pr, None)
        cols = {c: d for c, d in cols.items() if d}
        rows = {r: s for r, s in rows.items() if s}
        rank += 1
    if cols:
        rlist = sorted({r for d in cols.values() for r in d})
        rpos = {r: i for i, r in enumerate(rlist)}
        clist = sorted(cols)
        dense = [[0] * len(clist) for _ in rlist]
        for j, c in enumerate(clist):
            for r, v in cols[c].items():
                dense[rpos[r]][j] = v
        r2, divs = _dense_snf(dense)
        rank += r2
        return rank, [d for d in divs if d > 1]
    return rank, []


# ---------------------------------------------------------------------------
# 2. simplicial complexes and reduced homology over Z
# ---------------------------------------------------------------------------
def closure(facets):
    """Every subset of every listed facet, as a set of sorted tuples, including the empty face."""
    out = set()
    for F in facets:
        F = tuple(sorted(F))
        for k in range(len(F) + 1):
            out.update(combinations(F, k))
    out.add(())
    return out


def reduced_homology(faces):
    """{k: (betti_k, [torsion divisors])} for the AUGMENTED chain complex, so the answer is reduced
    homology: C_{-1} = Z is the empty face and d_0 is the augmentation.  Uses
        b_k = |C_k| - rank d_k - rank d_{k+1},     Tors H_k = divisors > 1 of d_{k+1}.
    """
    byd = {}
    for f in faces:
        if f == ():
            continue
        byd.setdefault(len(f) - 1, []).append(f)
    if not byd:
        return {}
    top = max(byd)
    pos = {}
    for k in byd:
        byd[k].sort()
        pos[k] = {f: i for i, f in enumerate(byd[k])}
    ranks = {}
    divs = {}
    for k in range(0, top + 1):
        cols = []
        for f in byd[k]:
            d = {}
            if k == 0:
                d[0] = 1                                  # augmentation to the empty face
            else:
                for t, v in enumerate(f):
                    g = f[:t] + f[t + 1:]
                    d[pos[k - 1][g]] = (-1) ** t
            cols.append(d)
        nrows = 1 if k == 0 else len(byd[k - 1])
        ranks[k], divs[k] = smith(cols, nrows)
    out = {}
    for k in range(0, top + 1):
        rk1 = ranks.get(k + 1, 0)
        out[k] = (len(byd[k]) - ranks[k] - rk1, divs.get(k + 1, []))
    return out


def fmt_hom(H):
    parts = []
    for k in sorted(H):
        b, t = H[k]
        if b or t:
            s = 'Z^%d' % b
            if t:
                s += '+' + '+'.join('Z/%d' % d for d in t)
            parts.append('H~_%d=%s' % (k, s))
    return ', '.join(parts) if parts else 'all reduced homology 0'


def nonzero_degrees(H):
    return {k: v for k, v in H.items() if v[0] or v[1]}


# ---------------------------------------------------------------------------
# 3. THE OBJECT PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
def vid(m, i, j):
    """The vertex (i, j) of K_m box P_n as an integer, exactly as the paper indexes it."""
    return j * m + i


def adjacency(m, n):
    """The host graph, straight from the paper's Section 2: (i,j) ~ (i',j') iff same column, or
    same row and consecutive columns."""
    N = m * n
    adj = [set() for _ in range(N)]
    for j in range(n):
        for i in range(m):
            for i2 in range(m):
                if i2 != i:
                    adj[vid(m, i, j)].add(vid(m, i2, j))
            for j2 in (j - 1, j + 1):
                if 0 <= j2 < n:
                    adj[vid(m, i, j)].add(vid(m, i, j2))
    return adj


def edges(m, n):
    adj = adjacency(m, n)
    return {(u, v) for u in range(m * n) for v in adj[u] if u < v}


def T0_facets(m, n, with_rows=True, columns=None):
    """The three families the paper's Theorem 3 prints, and nothing else."""
    cols = range(n) if columns is None else list(columns)
    F = []
    for j in cols:                                              # n column simplices
        F.append(tuple(sorted(vid(m, i, j) for i in range(m))))
    cl = list(cols)
    for t in range(len(cl) - 1):                                # (n-1) C(m,2) squares
        j, jj = cl[t], cl[t + 1]
        if jj != j + 1:
            continue
        for a in range(m):
            for b in range(a + 1, m):
                F.append(tuple(sorted([vid(m, a, j), vid(m, b, j),
                                       vid(m, a, jj), vid(m, b, jj)])))
    if with_rows:                                               # m(n-2) row triangles
        for i in range(m):
            for j in range(n - 2):
                F.append(tuple(sorted([vid(m, i, j), vid(m, i, j + 1), vid(m, i, j + 2)])))
    return F


def Y_complex(m, n, with_rows=True, columns=None):
    return closure(T0_facets(m, n, with_rows, columns))


def all_pairs_skeleton(N):
    out = {()}
    out.update((v,) for v in range(N))
    out.update(combinations(range(N), 2))
    return out


# ---------------------------------------------------------------------------
# 4. the Alexander dual FROM THE RULE, independently of T0
# ---------------------------------------------------------------------------
def _disconnected(adj, T):
    """Is the subgraph of G induced on the tuple T disconnected?"""
    T = set(T)
    if not T:
        return False
    start = next(iter(T))
    seen = {start}
    stack = [start]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w in T and w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) != len(T)


def dual_rule_ok(adj, S, k):
    """S is a face of Delta_k(G)* iff no k-subset of S induces a disconnected subgraph.
    (V \\ S is in Delta_k(G) iff V \\ S sits inside some facet V \\ T, T a disconnected k-set,
    i.e. iff some k-subset T of S is disconnected.)"""
    if len(S) < k:
        return True
    return all(not _disconnected(adj, T) for T in combinations(S, k))


def dual_exhaustive(m, n, k=3):
    """Every subset of V tested against the rule.  Only for small mn: 2^(mn) subsets."""
    N = m * n
    adj = adjacency(m, n)
    out = set()
    for mask in range(1 << N):
        S = tuple(v for v in range(N) if mask >> v & 1)
        if dual_rule_ok(adj, S, k):
            out.add(S)
    return out


def dual_grown(m, n, k=3):
    """The same set, grown vertex by vertex: a set is a face iff it is a face with its largest
    vertex removed AND every k-subset containing that vertex is connected.  Independent of T0."""
    N = m * n
    adj = adjacency(m, n)
    out = {()}
    frontier = [()]
    while frontier:
        nxt = []
        for S in frontier:
            start = (S[-1] + 1) if S else 0
            for v in range(start, N):
                T = S + (v,)
                if len(T) < k:
                    ok = True
                else:
                    ok = all(not _disconnected(adj, W + (v,))
                             for W in combinations(S, k - 1))
                if ok:
                    out.add(T)
                    nxt.append(T)
        frontier = nxt
    return out


def primal_cut_complex(m, n, k=3):
    """Delta_k(G) itself, with no duality anywhere: facets are the complements of the k-subsets
    that induce a disconnected subgraph."""
    N = m * n
    adj = adjacency(m, n)
    allv = set(range(N))
    facets = [tuple(sorted(allv - set(T)))
              for T in combinations(range(N), k) if _disconnected(adj, T)]
    return facets, closure(facets)


# ---------------------------------------------------------------------------
# 5. the two closed forms, and Alexander duality as bookkeeping
# ---------------------------------------------------------------------------
def A_of(m, n):
    v = (m * m - 3 * m + 2) * (n - 1)
    assert v % 2 == 0
    return v // 2


def B_of(m, n):
    v = m * (m * n - m - 2) * (n - 2)
    assert v % 2 == 0
    return v // 2


def conjectured(m, n):
    """Conjecture 4.4, as a dictionary degree -> (betti, torsion)."""
    return {m * n - 5: (A_of(m, n), []), m * n - 4: (B_of(m, n), [])}


def dual_to_primal(m, n, Hdual):
    """H~_i(D) from H~_j(D*) by combinatorial Alexander duality and universal coefficients:
        H~_i(D) = H~^{N-i-3}(D*) = Z^{b_j} (+) Tors H~_{j-1}(D*),   j = N-i-3.
    """
    N = m * n
    out = {}
    top = max(Hdual) if Hdual else -1
    for j in range(0, top + 2):
        b = Hdual.get(j, (0, []))[0]
        prev = Hdual.get(j - 1, (0, []))[1]
        if b or prev:
            out[N - j - 3] = (b, list(prev))
    return out


# ---------------------------------------------------------------------------
# 6. a two-variable integer polynomial, for the identities the proof needs
# ---------------------------------------------------------------------------
class Poly:
    """Sparse Z[m, n]; exponents are (a, b) pairs.  Equality is coefficient equality, so an
    identity checked here is an identity, not a sample of values."""

    def __init__(self, terms=None):
        self.t = {}
        for e, c in (terms or {}).items():
            if c:
                self.t[e] = self.t.get(e, 0) + c
        self.t = {e: c for e, c in self.t.items() if c}

    @staticmethod
    def const(c):
        return Poly({(0, 0): c})

    @staticmethod
    def var(which):
        return Poly({((1, 0) if which == 'm' else (0, 1)): 1})

    def __add__(self, o):
        o = o if isinstance(o, Poly) else Poly.const(o)
        r = dict(self.t)
        for e, c in o.t.items():
            r[e] = r.get(e, 0) + c
        return Poly(r)

    def __neg__(self):
        return Poly({e: -c for e, c in self.t.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Poly) else Poly.const(o)))

    def __mul__(self, o):
        o = o if isinstance(o, Poly) else Poly.const(o)
        r = {}
        for e1, c1 in self.t.items():
            for e2, c2 in o.t.items():
                e = (e1[0] + e2[0], e1[1] + e2[1])
                r[e] = r.get(e, 0) + c1 * c2
        return Poly(r)

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, o):
        return Poly.const(o) - self

    def __eq__(self, o):
        return self.t == (o if isinstance(o, Poly) else Poly.const(o)).t

    def subs(self, m=None, n=None):
        r = {}
        for (a, b), c in self.t.items():
            k = c
            ea, eb = a, b
            if m is not None:
                k *= m ** a
                ea = 0
            if n is not None:
                k *= n ** b
                eb = 0
            r[(ea, eb)] = r.get((ea, eb), 0) + k
        return Poly(r)

    def __str__(self):
        if not self.t:
            return '0'
        out = []
        for (a, b), c in sorted(self.t.items(), key=lambda x: (-x[0][0] - x[0][1], x[0])):
            s = ('%+d' % c)
            if a:
                s += 'm' + ('^%d' % a if a > 1 else '')
            if b:
                s += 'n' + ('^%d' % b if b > 1 else '')
            out.append(s)
        return ''.join(out)


# ---------------------------------------------------------------------------
# 7. the checks
# ---------------------------------------------------------------------------
def group_engine():
    head('A. ENGINE CONTROLS -- the homology engine, on objects with known answers')
    H = reduced_homology(closure([(0, 1, 2, 3)]))
    check('engine-full-tetrahedron-acyclic', not nonzero_degrees(H),
          '[Delta^3: %s]' % fmt_hom(H))

    H = reduced_homology(closure([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]))
    check('engine-S2-boundary-of-tetrahedron', nonzero_degrees(H) == {2: (1, [])},
          '[%s]' % fmt_hom(H))

    H = reduced_homology(closure([(0, 1), (1, 2), (0, 2)]))
    check('engine-S1-triangle', nonzero_degrees(H) == {1: (1, [])}, '[%s]' % fmt_hom(H))

    H = reduced_homology(closure([(0,), (1,), (2,)]))
    check('engine-three-points-b0-is-2', nonzero_degrees(H) == {0: (2, [])}, '[%s]' % fmt_hom(H))

    # the standard 6-vertex triangulation of the real projective plane: 10 triangles
    rp2 = [(0, 1, 2), (0, 1, 3), (0, 2, 4), (0, 3, 5), (0, 4, 5),
           (1, 2, 5), (1, 3, 4), (1, 4, 5), (2, 3, 4), (2, 3, 5)]
    H = reduced_homology(closure(rp2))
    check('engine-RP2-torsion-visible', nonzero_degrees(H) == {1: (0, [2])},
          '[%s -- an engine that could not see torsion would print H~_1=Z^0 here]' % fmt_hom(H))

    # a 7-vertex torus (Moebius-Kantor / Csaszar): b_1 = 2, b_2 = 1, no torsion
    torus = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2),
             (0, 1, 5), (1, 2, 6), (2, 3, 0), (3, 4, 1), (4, 5, 2), (5, 6, 3), (6, 0, 4)]
    H = reduced_homology(closure(torus))
    check('engine-7-vertex-torus', nonzero_degrees(H) == {1: (2, []), 2: (1, [])},
          '[%s]' % fmt_hom(H))


def group_graph():
    head('B. THE HOST GRAPH, from the paper\'s definition')
    bad = []
    for m in range(2, 9):
        for n in range(2, 9):
            e = len(edges(m, n))
            want = n * m * (m - 1) // 2 + m * (n - 1)
            if e != want:
                bad.append((m, n, e, want))
    check('graph-edge-count-n*C(m,2)+m*(n-1)', not bad,
          '[49 cells 2<=m,n<=8, e.g. (5,5): %d edges]' % len(edges(5, 5)))

    bad = []
    for m in range(3, 8):
        for n in range(3, 8):
            adj = adjacency(m, n)
            degs = sorted(len(adj[v]) for v in range(m * n))
            interior = sum(1 for j in range(n) if 0 < j < n - 1) * m
            want = sorted([m] * (2 * m) + [m + 1] * interior)
            if degs != want:
                bad.append((m, n))
    check('graph-degree-profile-m-or-m+1', not bad,
          '[end columns have degree m, interior columns m+1]')

    bad = []
    for m in range(3, 8):
        for n in range(3, 8):
            adj = adjacency(m, n)
            seen, st = {0}, [0]
            while st:
                u = st.pop()
                for w in adj[u]:
                    if w not in seen:
                        seen.add(w)
                        st.append(w)
            if len(seen) != m * n:
                bad.append((m, n))
    check('graph-connected', not bad, '[25 cells 3<=m,n<=7]')


def group_T0(cells_small, cells_T0, cells_census):
    head('C. THEOREM 3 (T0): the dual is the three printed families over the full 1-skeleton')

    for (m, n) in cells_small:
        ex = dual_exhaustive(m, n)
        gr = dual_grown(m, n)
        check('dual-rule-two-enumerations-agree-m%dn%d' % (m, n), ex == gr,
              '[all 2^%d subsets tested against the rule = %d faces = the grown enumeration]'
              % (m * n, len(ex)))

    for (m, n) in cells_T0:
        N = m * n
        D = dual_grown(m, n)
        Y = Y_complex(m, n)
        P = all_pairs_skeleton(N)
        same = (D == (Y | P))
        check('T0-identity-m%dn%d' % (m, n), same,
              '[D* == Y u (all pairs): %s, %d faces]' % (same, len(D)))

    bad = []
    for (m, n) in cells_T0:
        N = m * n
        D = dual_grown(m, n)
        maximal = [f for f in D if f and not any(set(f) < set(g) for g in D)]
        big = [f for f in maximal if len(f) >= 3]
        colsimp = [f for f in big if len(f) == m and len({v // m for v in f}) == 1]
        squares = [f for f in big if len(f) == 4 and len({v // m for v in f}) == 2 and
                   len({v % m for v in f}) == 2]
        rowtri = [f for f in big if len(f) == 3 and len({v % m for v in f}) == 1]
        want = (n, (n - 1) * m * (m - 1) // 2, m * (n - 2))
        got = (len(colsimp), len(squares), len(rowtri))
        if got != want or len(colsimp) + len(squares) + len(rowtri) != len(big) \
           or not all(len(f) == 2 for f in maximal if len(f) < 3):
            bad.append((m, n, got, want))
    check('T0-maximal-face-census', not bad,
          '[%d cells: exactly n column simplices, (n-1)C(m,2) squares, m(n-2) row triangles, '
          'every other maximal face an edge]' % len(cells_T0))

    # the f-vector, and the reduced Euler characteristic, on cells far outside brute-force range
    for (m, n), want in cells_census:
        D = Y_complex(m, n) | all_pairs_skeleton(m * n)
        fv = {}
        for f in D:
            if f:
                fv[len(f) - 1] = fv.get(len(f) - 1, 0) + 1
        vec = [fv.get(k, 0) for k in range(max(fv) + 1)]
        chi = -1 + sum((-1) ** k * fv.get(k, 0) for k in range(max(fv) + 1))
        ok = (want is None or vec == want) and chi == A_of(m, n) - B_of(m, n)
        check('dual-f-vector-and-euler-m%dn%d' % (m, n), ok,
              '[f = %s ; chi~ = %d = A - B = %d - %d]'
              % (vec, chi, A_of(m, n), B_of(m, n)))


def group_polynomials():
    head('D. THE CLOSED FORMS, as identities in Z[m, n] -- coefficients, not samples')
    m = Poly.var('m')
    n = Poly.var('n')

    twoA = (m * m - 3 * m + 2) * (n - 1)
    check('poly-2A-two-forms-agree', twoA == (n - 1) * (m - 1) * (m - 2),
          '[(m^2-3m+2)(n-1) == 2(n-1)C(m-1,2) : %s]' % twoA)

    twoB = m * (m * n - m - 2) * (n - 2)
    twoE = ((m * n) * (m * n - 1) - n * m * (m - 1)) - (m * m * (n - 1) + m * (n - 2)) * 2
    check('poly-2E-equals-2B', twoE == twoB,
          '[2(C(mn,2) - nC(m,2) - m^2(n-1) - m(n-2)) == m(mn-m-2)(n-2) == %s]' % twoB)

    check('poly-block-euler-gives-C(m-1,2)', (m * m - 3 * m) + 2 == (m - 1) * (m - 2),
          '[2(1 + (m^2-3m)/2) == 2C(m-1,2)]')

    N = m * n
    check('poly-duality-degree-mn-5-goes-to-2', (N - (N - 5) - 3) == Poly.const(2),
          '[N - i - 3 = 2 at i = mn-5]')
    check('poly-duality-degree-mn-4-goes-to-1', (N - (N - 4) - 3) == Poly.const(1),
          '[N - i - 3 = 1 at i = mn-4]')

    check('poly-boundary-line-n2', twoA.subs(n=2) == (m - 1) * (m - 2) and
          twoB.subs(n=2) == Poly.const(0),
          '[A(m,2) = C(m-1,2), B(m,2) = 0 -- the prism over a clique, Bayer et al. 2024]')
    check('poly-boundary-line-m2', twoA.subs(m=2) == Poly.const(0) and
          twoB.subs(m=2) == 4 * (n - 2) * (n - 2),
          '[A(2,n) = 0, B(2,n) = 2(n-2)^2 -- the grid G(2,n), Cut Complexes II]')

    bad = [(mm, nn) for mm in range(3, 41) for nn in range(3, 41)
           if (mm * mm - 3 * mm + 2) * (nn - 1) % 2 or mm * (mm * nn - mm - 2) * (nn - 2) % 2
           or A_of(mm, nn) < 0 or B_of(mm, nn) < 0]
    check('AB-integral-and-nonnegative', not bad,
          '[1444 cells 3<=m,n<=40: both halved products are even and both values >= 0]')


def group_table():
    head('E. THE SOURCE\'S OWN TABLE, and the four values this project retracted')
    table = {(3, 3): (2, 6), (3, 4): (3, 21), (3, 5): (4, 45), (3, 6): (5, 78),
             (4, 3): (6, 12), (4, 4): (9, 40), (4, 5): (12, 84),
             (5, 3): (12, 20), (5, 4): (18, 65), (6, 3): (20, 30)}
    bad = [(c, v, (A_of(*c), B_of(*c))) for c, v in table.items()
           if (A_of(*c), B_of(*c)) != v]
    check('source-table-ten-numeric-cells', not bad,
          '[all 10 SageMath cells of Table "3-cut Km[]Pn" reproduced as (A,B) in the printed '
          'degrees mn-5, mn-4]')

    starred = [(4, 6), (6, 4), (5, 5), (5, 6), (6, 5), (6, 6)]
    check('starred-cells-have-values', all(A_of(*c) > 0 and B_of(*c) > 0 for c in starred),
          '[the 6 cells the source marks "not known": %s]'
          % ', '.join('(%d,%d)->%d/%d' % (c[0], c[1], A_of(*c), B_of(*c)) for c in starred))

    retracted = {(6, 4): 90, (5, 5): 180, (5, 6): 280, (6, 5): 240}
    corrected = {(6, 4): 96, (5, 5): 135, (5, 6): 230, (6, 5): 198}
    check('retracted-B-values-rejected',
          all(B_of(*c) != v for c, v in retracted.items()) and
          all(B_of(*c) == v for c, v in corrected.items()),
          '[the four mis-transcribed B values 90/180/280/240 are NOT the closed form; the '
          'corrected 96/135/230/198 are]')


def group_block(mrange):
    head('F. STEP 3: the two-column block B_k')
    bad = []
    for m in mrange:
        Y = Y_complex(m, 2, with_rows=False)
        verts = sorted(f[0] for f in Y if len(f) == 1)
        pairs = {f for f in Y if len(f) == 2}
        if pairs != set(combinations(verts, 2)) or len(verts) != 2 * m:
            bad.append(m)
    check('block-1-skeleton-is-complete-on-2m-vertices', not bad,
          '[m = %d..%d: every one of the C(2m,2) pairs is an edge of B_k, which is what makes '
          'H_1(B_k) generated by the triangles of K_2m]' % (mrange[0], mrange[-1]))

    for m in mrange:
        Y = Y_complex(m, 2, with_rows=False)
        H = reduced_homology(Y)
        want = {2: ((m - 1) * (m - 2) // 2, [])}
        check('block-homology-m%d' % m, nonzero_degrees(H) == want,
              '[%s ; expected H~_2 = Z^C(m-1,2) = Z^%d and nothing else, no torsion]'
              % (fmt_hom(H), (m - 1) * (m - 2) // 2))

    # THE INTEGRAL ARGUMENT the review stage supplied, checked in the explicit form the paper
    # gives it: every triangle of the complete 1-skeleton K_{2m} either is already a 2-face of
    # B_k or is the boundary of a named 2-chain of exactly three 2-faces.  A rank count over a
    # field could not settle this; a chain identity with +-1 coefficients does.
    def bdy(t):
        t = tuple(sorted(t))
        return {(t[1], t[2]): 1, (t[0], t[2]): -1, (t[0], t[1]): 1}

    def combo(chain, signs):
        acc = {}
        for s, t in zip(signs, chain):
            for e, v in bdy(t).items():
                acc[e] = acc.get(e, 0) + s * v
        return {e: v for e, v in acc.items() if v}

    for m in mrange:
        Y = Y_complex(m, 2, with_rows=False)
        twofaces = {f for f in Y if len(f) == 3}
        verts = sorted(f[0] for f in Y if len(f) == 1)
        v = lambda i, j: vid(m, i, j)
        already, coned, broken = 0, 0, []
        for T in combinations(verts, 3):
            if T in twofaces:
                already += 1
                continue
            c0 = [w for w in T if w // m == 0]
            c1 = [w for w in T if w // m == 1]
            if len(c0) == 2:                       # {(i,k),(i',k),(x,k+1)}, x not in {i,i'}
                i, ip = (w % m for w in sorted(c0))
                x = c1[0] % m
                chain = [(v(i, 0), v(ip, 0), v(x, 0)),
                         (v(ip, 0), v(x, 0), v(x, 1)),
                         (v(i, 0), v(x, 0), v(x, 1))]
            elif len(c1) == 2:                     # {(i,k),(x,k+1),(y,k+1)}, i not in {x,y}
                x, y = (w % m for w in sorted(c1))
                i = c0[0] % m
                chain = [(v(i, 1), v(x, 1), v(y, 1)),
                         (v(i, 0), v(i, 1), v(x, 1)),
                         (v(i, 0), v(i, 1), v(y, 1))]
            else:
                broken.append(('unclassified', T))
                continue
            if any(tuple(sorted(s)) not in twofaces for s in chain):
                broken.append(('chain not made of 2-faces', T))
                continue
            target = bdy(T)
            if not any(combo(chain, sg) == target
                       for sg in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1),
                                  (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)]):
                broken.append(('no +-1 combination', T))
                continue
            coned += 1
        check('block-every-K2m-triangle-bounds-m%d' % m, not broken,
              '[all C(2m,3) = %d triples: %d are already 2-faces of B_k and the other %d are '
              'boundaries of the paper\'s three-term chain with coefficients +-1, so H_1(B_k) = 0 '
              'INTEGRALLY and Z^A carries no torsion]'
              % (len(verts) * (len(verts) - 1) * (len(verts) - 2) // 6, already, coned))

    bad = []
    for m in mrange:
        for n in (4, 5):
            Bk = Y_complex(m, n, with_rows=False, columns=[0, 1])
            Bk1 = Y_complex(m, n, with_rows=False, columns=[1, 2])
            col = closure([tuple(sorted(vid(m, i, 1) for i in range(m)))])
            if (Bk & Bk1) != col:
                bad.append((m, n))
    check('block-intersection-is-the-full-shared-column-simplex', not bad,
          '[B_k n B_{k+1} = the closure of column k+1, a non-empty full simplex, so reduced '
          'Mayer-Vietoris splits with no connecting map]')


def group_Y(cells):
    head('G. STEP 4: Y, the subcomplex generated by the faces of dimension >= 2')
    bad = []
    for (m, n) in cells:
        Y0 = Y_complex(m, n, with_rows=False)
        cur = set(Y0)
        for i in range(m):
            for j in range(n - 2):
                T = tuple(sorted([vid(m, i, j), vid(m, i, j + 1), vid(m, i, j + 2)]))
                new = closure([T]) - cur
                newe = [f for f in new if len(f) == 2]
                if len(new) != 2 or len(newe) != 1 \
                   or newe[0] != tuple(sorted([vid(m, i, j), vid(m, i, j + 2)])):
                    bad.append((m, n, i, j))
                cur |= closure([T])
        if cur != Y_complex(m, n):
            bad.append((m, n, 'closure'))
    check('row-triangle-is-an-elementary-expansion', not bad,
          '[%d cells: each row triangle adds exactly its own skip edge and itself, nothing '
          'else -- the expansion is elementary, coefficient +-1]' % len(cells))

    for (m, n) in cells:
        H = reduced_homology(Y_complex(m, n))
        want = {2: (A_of(m, n), [])}
        check('Y-homology-m%dn%d' % (m, n), nonzero_degrees(H) == want,
              '[%s ; expected H~_1 = 0 and H~_2 = Z^A = Z^%d only, no torsion]'
              % (fmt_hom(H), A_of(m, n)))


def group_dual_homology(cells):
    head('H. STEP 1 + STEP 5: the dual end to end, and the conjecture')
    for (m, n) in cells:
        N = m * n
        D = Y_complex(m, n) | all_pairs_skeleton(N)
        HD = reduced_homology(D)
        HY = reduced_homology(Y_complex(m, n))
        E = len(set(combinations(range(N), 2)) - Y_complex(m, n))
        Eform = N * (N - 1) // 2 - n * m * (m - 1) // 2 - m * m * (n - 1) - m * (n - 2)
        wedge_ok = (HD.get(1, (0, []))[0] == HY.get(1, (0, []))[0] + E and
                    all(HD.get(k, (0, []))[0] == HY.get(k, (0, []))[0]
                        for k in range(2, max(list(HD) + list(HY)) + 1)))
        check('step1-wedge-and-edge-count-m%dn%d' % (m, n),
              E == Eform == B_of(m, n) and wedge_ok,
              '[E = %d = closed form = B ; b_1(D*) = b_1(Y) + E and b_k(D*) = b_k(Y) for k >= 2]'
              % E)

        want = {1: (B_of(m, n), []), 2: (A_of(m, n), [])}
        check('dual-homology-m%dn%d' % (m, n), nonzero_degrees(HD) == want,
              '[%s ; expected H~_1 = Z^B = Z^%d, H~_2 = Z^A = Z^%d, nothing else, no torsion]'
              % (fmt_hom(HD), B_of(m, n), A_of(m, n)))

        got = dual_to_primal(m, n, HD)
        check('conjecture-4.4-m%dn%d' % (m, n), got == conjectured(m, n),
              '[H~_%d = Z^%d, H~_%d = Z^%d, all else 0 -- exactly Conjecture 4.4]'
              % (m * n - 5, A_of(m, n), m * n - 4, B_of(m, n)))


def group_primal(cells):
    head('I. NO DUALITY AT ALL: Delta_3 computed directly, against the published table')
    published = {(3, 3): {4: (2, []), 5: (6, [])},
                 (3, 4): {7: (3, []), 8: (21, [])},
                 (4, 3): {7: (6, []), 8: (12, [])}}
    for (m, n) in cells:
        facets, faces = primal_cut_complex(m, n, 3)
        H = reduced_homology(faces)
        nz = nonzero_degrees(H)
        ok = nz == published[(m, n)] and nz == conjectured(m, n)
        check('primal-delta3-m%dn%d' % (m, n), ok,
              '[%d facets, %d faces; %s -- matches the source\'s table and the closed forms, '
              'computed with no Alexander duality anywhere]' % (len(facets), len(faces),
                                                               fmt_hom(H)))
        D = dual_grown(m, n)
        HD = reduced_homology(D)
        check('duality-reproduces-the-primal-answer-m%dn%d' % (m, n),
              dual_to_primal(m, n, HD) == nz,
              '[the dual route gives the same groups in the same degrees]')


def group_antifit(cells):
    head('J. ANTI-FIT: the source\'s own k = 2 theorem, which this work did not produce')
    for (m, n) in cells:
        adj = adjacency(m, n)
        D2 = dual_grown(m, n, k=2)
        clique = closure([tuple(sorted(S)) for S in
                          _maximal_cliques(adj, m * n)])
        check('delta2-dual-is-the-clique-complex-m%dn%d' % (m, n), D2 == clique,
              '[Delta_2(G)* = {S : every pair in S is an edge}, %d faces]' % len(D2))
        H = reduced_homology(D2)
        got = dual_to_primal(m, n, H)
        want = {m * n - 4: ((m - 1) * (n - 1), [])}
        check('delta2-is-a-wedge-of-(m-1)(n-1)-spheres-m%dn%d' % (m, n), got == want,
              '[%s -> H~_%d(Delta_2) = Z^%d: Chauhan-Shukla-Vinayak\'s own Theorem for k = 2, '
              'reproduced by the same machinery]' % (fmt_hom(H), m * n - 4,
                                                    (m - 1) * (n - 1)))

    facets, faces = primal_cut_complex(3, 3, 2)
    H = reduced_homology(faces)
    check('primal-delta2-m3n3', nonzero_degrees(H) == {5: (4, [])},
          '[%d facets, %d faces, %s = wedge of (m-1)(n-1) = 4 spheres S^5, computed directly]'
          % (len(facets), len(faces), fmt_hom(H)))


def _maximal_cliques(adj, N):
    """Every maximal clique, by plain Bron-Kerbosch.  Used only to build the clique complex a
    second way in the k = 2 anti-fit check."""
    out = []

    def bk(R, P, X):
        if not P and not X:
            out.append(tuple(sorted(R)))
            return
        for v in list(P):
            bk(R | {v}, P & adj[v], X & adj[v])
            P = P - {v}
            X = X | {v}

    bk(set(), set(range(N)), set())
    return out


def group_controls():
    head('K. FORCED NEGATIVES: the comparator must reject wrong predictions')
    m, n = 4, 4
    HD = reduced_homology(Y_complex(m, n) | all_pairs_skeleton(m * n))
    got = dual_to_primal(m, n, HD)
    truth = conjectured(m, n)

    wrong = {m * n - 5: (A_of(m, n) + 1, []), m * n - 4: (B_of(m, n), [])}
    check('control-rejects-A-plus-one', got != wrong and got == truth,
          '[at (4,4) the comparator says NO to Z^%d at degree %d]'
          % (A_of(m, n) + 1, m * n - 5))

    wrong = {m * n - 5: (A_of(m, n), []), m * n - 4: (B_of(m, n) + 1, [])}
    check('control-rejects-B-plus-one', got != wrong,
          '[and NO to Z^%d at degree %d]' % (B_of(m, n) + 1, m * n - 4))

    wrong = {m * n - 5: (B_of(m, n), []), m * n - 4: (A_of(m, n), [])}
    check('control-rejects-swapped-degrees', got != wrong,
          '[and NO to the two groups exchanged]')

    wrong = {m * n - 5: (A_of(m, n), [2]), m * n - 4: (B_of(m, n), [])}
    check('control-rejects-spurious-torsion', got != wrong,
          '[and NO to Z^%d + Z/2 at degree %d, so the freeness claim is a decided answer and '
          'not an unchecked default]' % (A_of(m, n), m * n - 5))


# ---------------------------------------------------------------------------
# 8. scope, and the driver
# ---------------------------------------------------------------------------
SCOPE = """
NOT RE-RUN: the theorem itself.  The paper's claim is a statement about ALL m, n >= 3 and it is
settled by the hand proof of Section 4, not by any census.  Every cell this program computes is a
finite check of that proof's intermediate claims and of its two closed forms.  A program cannot
verify an infinite family, and nothing below should be read as if it had.

NOT RE-RUN: the large census.  The discovery run computed the reduced homology of D* on 66
distinct cells, out to (10,10) at 100 vertices, on AWS instances.  This program re-derives the
same quantities on the much smaller list printed above, chosen so that a referee can run it in
under a minute on a laptop with nothing installed.  The cells (10,9) and (10,10) appear here only
in the f-vector and Euler-characteristic check, not in a homology computation.

NOT RE-RUN: Delta_3 computed directly, with no Alexander duality, at any cell beyond (3,3),
(3,4) and (4,3).  Those three are the three smallest cells of the source's table and they are
where the duality step is pinned to an independently computed answer; larger cells are checked
through the dual only.

NOT RE-RUN: the exhaustive-subset enumeration of D* beyond mn <= 16.  For larger cells the dual
is enumerated by the incremental rule instead, and the two enumerations are shown to agree
wherever both are affordable.

NOT RE-RUN: anything bibliographic.  The wording of Conjecture 4.4, its position at lines 871-880
of the source's 80,146-byte main .tex file, the numbering that follows from the shared theorem
counter, the contents of the source's table, the two boundary lines quoted from Bayer et al.
(2024) and Cut Complexes II (2025), and the prior-art search, were established by fetching and
reading those sources.  This program checks mathematics, not provenance, and fetches nothing.

NOT RE-RUN: the total k-cut complex Delta_3^t, and with it the source's sibling Conjectures 4.3,
4.5 and 4.6.  Nothing in this program touches the total operator.
"""


def main():
    print('verify.py -- the 3-cut complex of K_m box P_n')
    print('Python %s' % sys.version.split()[0])
    print('exact integer arithmetic only; homology by integer Smith normal form')

    group_engine()
    group_graph()
    group_T0(cells_small=[(3, 3), (3, 4), (4, 3), (3, 5), (5, 3), (4, 4)],
             cells_T0=[(3, 3), (3, 4), (4, 3), (4, 4), (5, 3), (3, 5), (5, 5), (6, 4), (4, 6),
                       (3, 7), (7, 3), (6, 6), (7, 5), (5, 7), (8, 4)],
             cells_census=[((3, 8), [24, 276, 110, 21]),
                           ((8, 3), [24, 276, 400, 266, 168, 84, 24, 3]),
                           ((10, 9), None), ((10, 10), None)])
    group_polynomials()
    group_table()
    group_block(mrange=[3, 4, 5, 6, 7])
    group_Y(cells=[(3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 12), (4, 3), (4, 4), (4, 5),
                   (4, 6), (5, 3), (5, 4), (5, 5), (6, 3), (6, 4), (7, 3), (8, 3), (9, 3),
                   (3, 20)])
    # ⭐ this list contains ALL SIX cells the source's table marks as not known -- (4,6), (6,4),
    # (5,5), (5,6), (6,5), (6,6) -- as well as (3,7) and (7,3), which lie outside that table.
    group_dual_homology(cells=[(3, 3), (3, 4), (4, 3), (4, 4), (3, 5), (5, 3), (3, 6), (6, 3),
                               (4, 5), (5, 4), (5, 5), (4, 6), (6, 4), (5, 6), (6, 5), (6, 6),
                               (3, 7), (7, 3)])
    group_primal(cells=[(3, 3), (3, 4), (4, 3)])
    group_antifit(cells=[(3, 3), (3, 4), (4, 3), (4, 4), (5, 3), (3, 5)])
    group_controls()

    print()
    print('=== SCOPE: what this program does NOT cover ===')
    for para in SCOPE.strip().split('\n\n'):
        print(para)
        print()
    print('VERDICT: ALL %d CHECKS PASS' % _PASS if not _FAIL
          else 'VERDICT: %d CHECKS FAILED of %d' % (_FAIL, _PASS + _FAIL))
    return 0 if not _FAIL else 1


if __name__ == '__main__':
    sys.exit(main())

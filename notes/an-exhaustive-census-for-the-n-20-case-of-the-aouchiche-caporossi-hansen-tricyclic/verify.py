#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- referee verification program for

    "An Exhaustive Census for the n=20 Case of the Aouchiche-Caporossi-Hansen
     Tricyclic-Energy Conjecture"

Standard library only (no numpy/sympy/nauty, no data files).  Every decision is
taken in exact integer or rational arithmetic; floating point is used only to
format numbers for printing.  Run:  python3 verify.py

VALUES TAKEN FROM THE PAPER (inputs, never checked against themselves)
---------------------------------------------------------------------
  * the paper's explicit model of P = P_20^{6,6,6} (three 6-cycles on
    {0..5},{6..11},{12..17} plus the edges {0,18},{6,18},{18,19},{12,19})
    and Q = P-{6,18}+{5,10};
  * the two factored characteristic polynomials, eq. (1) and eq. (2);
  * the decimal strings 27.2903786151832973536..., 27.2773174852287231601...,
    0.0130611299545741935..., the rounded values in Result 1 and the abstract's
    separation 1.3061e-2 / 1.31e-2;
  * the census total 3 288 208 176 and the claimed articulation-vertex counts
    5 (for P) and 6 (for Q).

INPUTS NOT TAKEN FROM THE PAPER (so that no reader mistakes them for quotes)
---------------------------------------------------------------------------
The paper prints no graph6 string anywhere; its only use of the word is the
phrase "for each graph6 record".  The three graph6 strings in the input block
below are supplied with this program and are NOT quotations:
  * PAPER_G6_P and PAPER_G6_Q are this program's own encodings of P and of Q
    (the PAPER_ prefix marks the paper's *graphs*, not a printed string).  They
    are not taken on trust: A3 pins the decoded PAPER_G6_P edge for edge to the
    paper's explicit model, and B2 pins the decoded PAPER_G6_Q to
    P-{6,18}+{5,10}, so a wrong string fails a check rather than sliding past.
  * G6_P_ALT_LABELLING is a differently labelled copy of P, carried as a
    relabelling test vector for the isomorphism machinery.  It is believed to be
    the output of nauty's labelg, but nauty is not available here, so this
    program neither verifies nor relies on that provenance: it does not check
    that the string is any canonical form.  Checks B8, B9 and B10 assert only
    (B8) that the record is a well-formed connected 20-vertex 22-edge graph,
    (B9) that it is isomorphic to P, by a cheap invariant and independently by
    an explicitly exhibited permutation, and (B10) that graph6 decoding and
    re-encoding round-trip on it and that no one-character variant of it names
    the same graph.  Nothing else in the program depends on this string, and no
    claim of the paper rests on those three checks.

WHAT THIS PROGRAM DERIVES (the checks)
--------------------------------------
  * decodes the graph6 strings and rebuilds P and Q from the paper's words:
    order, size, connectivity, tricyclicity (m=n+2, cycle rank 3), degree
    sequence, bipartiteness, the three vertex-disjoint 6-cycles and the
    attachment pattern of the ACH family, articulation vertices;
  * proves the alternatively labelled copy of P (see above) is isomorphic to P,
    and that Q is not, TWICE: once with a cheap invariant for the graph class,
    and once with a generic degree-preserving backtracking isomorphism search
    that exhibits an explicit permutation in the positive case and exhausts the
    search space in the negative case.  That the cheap invariant really is a
    *complete* invariant is not assumed: it is validated exhaustively against
    the backtracking search on T(10) (both directions, all five base
    multigraphs) and, in the merging direction that could hide a maximiser, on
    all 58461 configurations of T(20);
  * recomputes both characteristic polynomials from the adjacency matrices by
    Faddeev-LeVerrier in exact integer arithmetic, confirms each of them with a
    second, algorithmically independent engine (det(tI-A) by fraction-free
    Bareiss elimination at n+1 integer points), and compares them with the
    paper's factored forms; the low-order coefficients are additionally checked
    against independently enumerated edge, triangle, 4-cycle and disjoint-edge-
    pair counts (Sachs' coefficient theorem);
  * recomputes E(P), E(Q) and E(P)-E(Q) as *rigorous rational enclosures*:
    square-free decomposition over Z, root isolation by exact derivative
    interlacing, bisection with exact integer sign tests; the enclosures are
    then compared with the paper's digit strings;
  * re-derives the census total 3 288 208 176 independently of geng, by
    Burnside/Polya counting on the pair group S_20^(2) plus the multiset
    (Euler) logarithm that passes from all graphs to connected graphs; the
    counting machinery is itself validated against brute-force isomorph
    counting for n<=6 and against the classical totals A000088/A001349;
  * re-runs the energy maximisation exhaustively over two sub-families of the
    census that fit in the time budget (see the NOT REPRODUCED line printed at
    the end): the class T(20) of all connected 20-vertex 22-edge graphs with
    2<=deg<=3 (8157 isomorphism classes, about 2.5e-6 of the census, containing
    both P and Q), and the full one-edge-swap neighbourhood of P (3363
    connected graphs, most of them outside T(20)).  In both slices P is found
    to be the strict maximiser up to isomorphism and Q the strict runner-up.
    The energy scan over the whole census is NOT re-run; the program says so.

Total runtime about 6 minutes single-process (python3.9, 2021 laptop core).
Running with --quick omits sections E and F and therefore exits 2, not 0: a
zero exit status is only ever produced by the complete verification.
"""
import sys, time
from fractions import Fraction
from math import gcd
from itertools import combinations, combinations_with_replacement, permutations

sys.setrecursionlimit(10000)
_T0 = time.time()
_RESULTS = []          # (name, ok) in order
_SKIPPED = []


# ----------------------------------------------------------------------------
# INPUTS: everything quoted from the paper, in one block.
# ----------------------------------------------------------------------------
PAPER_N = 20
PAPER_M = 22
PAPER_G6_P = "ShEG?C@?G?_P????_?G?@??C?AKC???CC"
# not from the paper: an alternative labelling of P (believed to be labelg's
# canonical form, but that provenance is NOT checked here -- see the header).
G6_P_ALT_LABELLING = "SC`A@?_?O@A?@??_@?????EC?C_?E???k"
PAPER_G6_Q = "ShEG?C@?GG_P????_?G?@??C?AK????CC"
PAPER_CYCLES = [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11], [12, 13, 14, 15, 16, 17]]
PAPER_EXTRA_EDGES = [(0, 18), (6, 18), (18, 19), (12, 19)]
PAPER_Q_DEL = (6, 18)
PAPER_Q_ADD = (5, 10)
# eq. (1):  (x^2-4)(x^2-1)^4 (x^4-7x^2+8)(x^6-7x^4+12x^2-2)
PAPER_PHI_P = [([1, 0, -4], 1), ([1, 0, -1], 4), ([1, 0, -7, 0, 8], 1),
               ([1, 0, -7, 0, 12, 0, -2], 1)]
# eq. (2):  (x^2-1)^2 (x^16-20x^14+162x^12-686x^10+1638x^8-2219x^6+1621x^4-557x^2+64)
PAPER_PHI_Q = [([1, 0, -1], 2),
               ([1, 0, -20, 0, 162, 0, -686, 0, 1638, 0, -2219, 0, 1621, 0,
                 -557, 0, 64], 1)]
PAPER_E_P = "27.2903786151832973536"       # ... (truncated in the paper)
PAPER_E_Q = "27.2773174852287231601"       # ...
PAPER_GAP = "0.0130611299545741935"        # ...
PAPER_E_P_ROUND = "27.2903786151833"       # Result 1 / abstract
PAPER_E_Q_ROUND = "27.2773174852287"       # Result 1
PAPER_SEP_5SIG = Fraction(13061, 10 ** 6)  # abstract: approx 1.3061e-2
PAPER_SEP_3SIG = Fraction(131, 10 ** 4)    # status section: approx 1.31e-2
PAPER_CENSUS_CLASSES = 3288208176
PAPER_SHARDS = 997
PAPER_ART_P = 5
PAPER_ART_Q = 6


def check(name, ok, detail=""):
    """Record one check.  ok must be a real computed boolean."""
    ok = bool(ok)
    _RESULTS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line, flush=True)
    return ok


def skip(name, why):
    _SKIPPED.append(name)
    print("SKIP " + name + " [" + why + "]", flush=True)


def finish():
    n = len(_RESULTS)
    bad = [nm for nm, ok in _RESULTS if not ok]
    if _SKIPPED:
        print("WARNING: %d checks were SKIPPED (%s); this run is NOT the "
              "published verification" % (len(_SKIPPED), ",".join(_SKIPPED)))
    print("elapsed: %.1f s" % (time.time() - _T0))
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    if _SKIPPED:
        # a green exit status must be unobtainable without the load-bearing
        # sections (the T(20) maximisation and the invariant validation).
        print("VERDICT: INCOMPLETE RUN (%d skipped) -- exiting 2" % len(_SKIPPED))
        sys.exit(2)
    sys.exit(0)


# ----------------------------------------------------------------------------
# graph6 and elementary graph routines
# ----------------------------------------------------------------------------
def g6_decode(s):
    """graph6 -> (n, sorted edge list).  Raises on malformed input."""
    b = [ord(c) - 63 for c in s]
    if not b or min(b) < 0 or max(b) > 63:
        raise ValueError("bad graph6 characters")
    n = b[0]
    if not 0 <= n <= 62:
        raise ValueError("only n<=62 supported")
    need = (n * (n - 1) // 2 + 5) // 6
    if len(b) - 1 != need:
        raise ValueError("graph6 payload length %d, expected %d" % (len(b) - 1, need))
    bits = []
    for v in b[1:]:
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    edges = []
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                edges.append((i, j))
            idx += 1
    if any(bits[idx:]):                       # padding must be zero
        raise ValueError("nonzero graph6 padding bits")
    return n, sorted(edges)


def g6_encode(n, edges):
    es = set((min(a, b), max(a, b)) for a, b in edges)
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if (i, j) in es else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = 2 * v + b
        out.append(chr(v + 63))
    return "".join(out)


def adjacency(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        if a == b or not (0 <= a < n and 0 <= b < n):
            raise ValueError("bad edge")
        adj[a].append(b)
        adj[b].append(a)
    return adj


def n_components(n, adj, skip_v=None):
    seen = [False] * n
    comp = 0
    for s in range(n):
        if seen[s] or s == skip_v:
            continue
        comp += 1
        seen[s] = True
        st = [s]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w != skip_v and not seen[w]:
                    seen[w] = True
                    st.append(w)
    return comp


def is_connected(n, adj):
    return n > 0 and n_components(n, adj) == 1


def articulation_vertices(n, adj):
    """cut vertices, computed by brute-force removal (n is small here): v is a
    cut vertex iff deleting it increases the number of components."""
    base = n_components(n, adj)
    return [v for v in range(n) if n_components(n, adj, skip_v=v) > base]


def cycle_rank_by_spanning_forest(n, edges):
    """The cycle rank (first Betti number), computed by actually building a
    spanning forest and counting the edges left over -- not by restating
    m = n + 2.  Returns (rank, number of components)."""
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    extra = 0
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra == rb:
            extra += 1                        # closes a cycle: not a forest edge
        else:
            par[ra] = rb
    return extra, len({f(v) for v in range(n)})


def is_bipartite(n, adj):
    col = [-1] * n
    for s in range(n):
        if col[s] >= 0:
            continue
        col[s] = 0
        st = [s]
        while st:
            v = st.pop()
            for w in adj[v]:
                if col[w] < 0:
                    col[w] = 1 - col[v]
                    st.append(w)
                elif col[w] == col[v]:
                    return False
    return True


def induced_edges(edges, vs):
    s = set(vs)
    return sorted((a, b) for a, b in edges if a in s and b in s)


def is_cycle_on(vs, edges):
    """True iff the given edge set is exactly a single cycle through all of vs."""
    k = len(vs)
    if k < 3 or len(edges) != k:
        return False
    adj = {v: [] for v in vs}
    for a, b in edges:
        if a not in adj or b not in adj:
            return False
        adj[a].append(b)
        adj[b].append(a)
    if any(len(adj[v]) != 2 for v in vs):
        return False
    start = vs[0]
    prev, cur, cnt = start, adj[start][0], 1
    while cur != start:
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        prev, cur = cur, nxt
        cnt += 1
        if cnt > k:
            return False
    return cnt == k


# ----------------------------------------------------------------------------
# exact polynomial arithmetic over Z (coefficients high degree first)
# ----------------------------------------------------------------------------
def pnorm(p):
    i = 0
    while i < len(p) - 1 and p[i] == 0:
        i += 1
    p = list(p[i:])
    return p if p else [0]


def pzero(p):
    return len(p) == 1 and p[0] == 0


def pprim(p):
    p = pnorm(p)
    if pzero(p):
        return [0]
    g = 0
    for a in p:
        g = gcd(g, abs(a))
    if g > 1:
        p = [a // g for a in p]
    if p[0] < 0:
        p = [-a for a in p]
    return p


def pderiv(p):
    d = len(p) - 1
    if d <= 0:
        return [0]
    return pnorm([p[i] * (d - i) for i in range(d)])


def pmul(a, b):
    if pzero(a) or pzero(b):
        return [0]
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def ppow(a, k):
    r = [1]
    for _ in range(k):
        r = pmul(r, a)
    return r


def pdivmod(a, b):
    """polynomial division over Q; returns (quotient, remainder) as Fractions."""
    a = [Fraction(x) for x in pnorm(a)]
    b = [Fraction(x) for x in pnorm(b)]
    if len(b) == 1 and b[0] == 0:
        raise ZeroDivisionError("division by the zero polynomial")
    da, db = len(a) - 1, len(b) - 1
    if da < db:
        return [Fraction(0)], a
    q = [Fraction(0)] * (da - db + 1)
    r = a[:]
    for i in range(da - db + 1):
        c = r[i] / b[0]
        q[i] = c
        if c:
            for j in range(db + 1):
                r[i + j] -= c * b[j]
    rem = r[da - db + 1:]
    while len(rem) > 1 and rem[0] == 0:
        rem = rem[1:]
    return q, (rem if rem else [Fraction(0)])


def toint(p):
    den = 1
    for x in p:
        den = den * x.denominator // gcd(den, x.denominator)
    return pnorm([int(x * den) for x in p])


def pgcd(a, b):
    a, b = pprim(a), pprim(b)
    while not pzero(b):
        _, r = pdivmod(a, b)
        a, b = b, pprim(toint(r))
    return pprim(a)


def pquo_exact(a, b):
    q, r = pdivmod(a, b)
    if not all(x == 0 for x in r):
        raise ArithmeticError("inexact polynomial division")
    return toint(q)


def squarefree(p):
    """square-free decomposition over Z: [(f_i, i)] with prim(p) = prod f_i^i.

    Uses g_k = gcd(g_{k-1}, g_{k-1}') (characteristic 0), so g_k is the product
    of the irreducible factors of p with multiplicity reduced by k.  The caller
    must verify the result with sf_check (the program does)."""
    p = pprim(p)
    h, g = [], p
    while len(g) > 1:
        gn = pgcd(g, pderiv(g))
        h.append(pquo_exact(g, gn))
        g = gn
        if len(h) > 40:
            raise ArithmeticError("square-free decomposition did not terminate")
    out = []
    for i in range(len(h)):
        nxt = h[i + 1] if i + 1 < len(h) else [1]
        f = pprim(pquo_exact(h[i], nxt))
        if len(f) > 1:
            out.append((f, i + 1))
    return out


def sf_check(p, facs):
    r = [1]
    for f, m in facs:
        r = pmul(r, ppow(f, m))
    return pprim(r) == pprim(p)


def sign_at(p, num, den):
    """exact sign of p(num/den) for den>0: evaluates den^deg(p) * p(num/den)."""
    v, qp = 0, 1
    for a in p:
        v = v * num + a * qp if a else v * num
        qp *= den
    return (v > 0) - (v < 0)


def _mid(lo, hi):
    n = lo[0] * hi[1] + hi[0] * lo[1]
    d = 2 * lo[1] * hi[1]
    g = gcd(abs(n), d)
    if g > 1:
        n //= g
        d //= g
    return (n, d)


ROOT_BOUND = 32          # every adjacency eigenvalue of a graph with n<=32 is inside


def isolate(p, maxref=400):
    """Isolating intervals for the roots of a square-free integer polynomial.

    Returns a list of deg(p) disjoint intervals ((num,den),(num,den)), sorted,
    each containing exactly one root -- which simultaneously *proves* that all
    roots are real and simple.  Returns None if that proof is not reached.

    Method (exact, no floating point): the roots of p' separate the roots of p
    when p is real-rooted, so isolate p' recursively and use points strictly
    inside those intervals as separators; strict sign alternation of p over the
    resulting subdivision of [-B,B] certifies one root per part."""
    d = len(p) - 1
    if d <= 0:
        return []
    if d == 1:
        num, den = -p[1], p[0]
        if den < 0:
            num, den = -num, -den
        g = gcd(abs(num), den)
        if g > 1:
            num //= g
            den //= g
        return [((num, den), (num, den))]
    dv = pprim(pderiv(p))
    div = isolate(dv, maxref)
    if div is None:
        return None
    ivs = [[a, b] for a, b in div]
    for _ in range(maxref):
        pts = [(-ROOT_BOUND, 1)] + [_mid(lo, hi) for lo, hi in ivs] + [(ROOT_BOUND, 1)]
        sg = [sign_at(p, a, b) for a, b in pts]
        if 0 not in sg and all(sg[i] != sg[i + 1] for i in range(len(sg) - 1)):
            return [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
        for k, (lo, hi) in enumerate(ivs):        # one bisection everywhere
            if lo == hi:
                continue
            m = _mid(lo, hi)
            s = sign_at(dv, m[0], m[1])
            if s == 0:
                ivs[k] = [m, m]
            elif s == sign_at(dv, lo[0], lo[1]):
                ivs[k] = [m, hi]
            else:
                ivs[k] = [lo, m]
    return None


def refine(p, iv, eps, maxit=5000):
    """bisect an isolating interval until its width is <= eps (exact signs)."""
    lo, hi = iv
    if lo == hi:
        return lo, hi
    lo, hi = list(lo), list(hi)
    slo = sign_at(p, lo[0], lo[1])
    for _ in range(maxit):
        if Fraction(hi[0], hi[1]) - Fraction(lo[0], lo[1]) <= eps:
            break
        m = _mid(tuple(lo), tuple(hi))
        s = sign_at(p, m[0], m[1])
        if s == 0:
            return m, m
        if s == slo:
            lo = list(m)
        else:
            hi = list(m)
    return tuple(lo), tuple(hi)


def energy_interval(p, eps_total=Fraction(1, 10 ** 5)):
    """Rigorous rational enclosure [lo,hi] of sum |root| over all roots of the
    integer polynomial p, counted with multiplicity.  Returns None if the
    routine cannot certify that p has deg(p) real roots."""
    facs = squarefree(p)
    if not sf_check(p, facs):
        return None
    deg = len(pprim(p)) - 1
    if sum(m * (len(f) - 1) for f, m in facs) != deg:
        return None
    eps = Fraction(eps_total, 2 * max(deg, 1))
    lo = hi = Fraction(0)
    for f, m in facs:
        ivs = isolate(f)
        if ivs is None or len(ivs) != len(f) - 1:
            return None
        for iv in ivs:
            a, b = refine(f, iv, eps)
            A, B = Fraction(a[0], a[1]), Fraction(b[0], b[1])
            if A >= 0:
                l, h = A, B
            elif B <= 0:
                l, h = -B, -A
            else:
                l, h = Fraction(0), max(-A, B)
            lo += m * l
            hi += m * h
    return lo, hi


def charpoly(n, adj):
    """char. polynomial of the adjacency matrix, exact integers, high degree
    first, by Faddeev-LeVerrier with a sparse left multiplication by A."""
    c = [1] + [0] * n
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in adj[i]:
            M[i][j] += 1
    c[1] = -sum(M[i][i] for i in range(n))
    for k in range(1, n):
        N = [row[:] for row in M]
        for i in range(n):
            N[i][i] += c[k]
        Mn = [[0] * n for _ in range(n)]
        for i in range(n):
            r = Mn[i]
            for j in adj[i]:
                Nj = N[j]
                for t in range(n):
                    r[t] += Nj[t]
        tr = sum(Mn[i][i] for i in range(n))
        if tr % (k + 1):
            raise ArithmeticError("Faddeev-LeVerrier lost integrality")
        c[k + 1] = -tr // (k + 1)
        M = Mn
    return c


def trace_power(n, adj, k):
    """tr(A^k), exact."""
    M = [[1 if j == i else 0 for j in range(n)] for i in range(n)]
    for _ in range(k):
        Mn = [[0] * n for _ in range(n)]
        for i in range(n):
            r = Mn[i]
            for j in adj[i]:
                Mj = M[j]
                for t in range(n):
                    r[t] += Mj[t]
        M = Mn
    return sum(M[i][i] for i in range(n))


# ----------------------------------------------------------------------------
# SECOND, INDEPENDENT characteristic-polynomial engine.
#
# Faddeev-LeVerrier builds the coefficients *out of* the traces tr(A^k) via
# Newton's identities, so comparing the polynomial's Newton power sums with the
# traces is an identity of that algorithm, not a confirmation of it (and it
# constrains only the first few coefficients in any case).  The engine below
# shares no code and no algebra with it: it evaluates det(tI - A) at n+1 integer
# points by fraction-free (Bareiss) Gaussian elimination over Z.  Two degree-n
# polynomials that agree at n+1 points are equal, so agreement here pins down
# every coefficient of the characteristic polynomial.
# ----------------------------------------------------------------------------
def det_bareiss(M):
    """exact determinant of an integer matrix; every division is verified."""
    n = len(M)
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            piv = -1
            for i in range(k + 1, n):
                if M[i][k]:
                    piv = i
                    break
            if piv < 0:
                return 0                      # a zero column: det = 0
            M[k], M[piv] = M[piv], M[k]
            sign = -sign
        akk = M[k][k]
        Mk = M[k]
        for i in range(k + 1, n):
            Mi = M[i]
            aik = Mi[k]
            for j in range(k + 1, n):
                q, r = divmod(Mi[j] * akk - aik * Mk[j], prev)
                if r:
                    raise ArithmeticError("Bareiss elimination lost integrality")
                Mi[j] = q
            Mi[k] = 0
        prev = akk
    return sign * M[n - 1][n - 1]


def peval(p, t):
    """p(t), exact, coefficients high degree first."""
    v = 0
    for a in p:
        v = v * t + a
    return v


def charpoly_confirm(n, adj, cp):
    """True iff cp really is det(xI - A): checked at the n+1 points -1..n-1 with
    the independent Bareiss determinant engine."""
    if len(cp) != n + 1 or cp[0] != 1:
        return False
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in adj[i]:
            A[i][j] += 1
    for t in range(-1, n):
        M = [[(t if i == j else 0) - A[i][j] for j in range(n)] for i in range(n)]
        if det_bareiss(M) != peval(cp, t):
            return False
    return True


def count_c4_direct(n, adj):
    """number of 4-cycles, counted combinatorially (no matrix powers): every
    4-cycle u-a-v-b-u is seen once from the pair {u,v} and once from {a,b}."""
    nb = [set(a) for a in adj]
    tot = 0
    for u in range(n):
        for v in range(u + 1, n):
            c = len(nb[u] & nb[v])
            tot += c * (c - 1) // 2
    if tot % 2:
        raise ArithmeticError("4-cycle double count is odd")
    return tot // 2


def count_disjoint_edge_pairs(edges):
    """number of unordered pairs of edges with no common endpoint."""
    tot = 0
    for i in range(len(edges)):
        a, b = edges[i]
        for j in range(i + 1, len(edges)):
            c, d = edges[j]
            if a != c and a != d and b != c and b != d:
                tot += 1
    return tot


# ----------------------------------------------------------------------------
# SECOND, INDEPENDENT isomorphism engine: a generic degree-preserving
# backtracking search.  It uses nothing but "an isomorphism preserves degrees
# and adjacency", so it is independent of the degree-2 suppression invariant
# used elsewhere; it is complete, i.e. it returns None only when no isomorphism
# exists.  Positive answers come with an explicit permutation, which the caller
# re-verifies edge by edge.
# ----------------------------------------------------------------------------
def iso_map(n, E1, E2):
    """a permutation p with p(E1) = E2, or None if the graphs are not
    isomorphic (the search is exhaustive)."""
    A1 = [[0] * n for _ in range(n)]
    A2 = [[0] * n for _ in range(n)]
    for a, b in E1:
        A1[a][b] = A1[b][a] = 1
    for a, b in E2:
        A2[a][b] = A2[b][a] = 1
    d1 = [sum(r) for r in A1]
    d2 = [sum(r) for r in A2]
    if sorted(d1) != sorted(d2):
        return None
    order, placed = [], [False] * n
    for _ in range(n):                       # prefer vertices already attached
        best, key = -1, None
        for v in range(n):
            if placed[v]:
                continue
            k = (-sum(1 for w in order if A1[v][w]), -d1[v], v)
            if key is None or k < key:
                key, best = k, v
        order.append(best)
        placed[best] = True
    mp = [-1] * n
    used = [False] * n

    def rec(i):
        if i == n:
            return True
        v = order[i]
        for w in range(n):
            if used[w] or d2[w] != d1[v]:
                continue
            ok = True
            for t in range(i):
                u = order[t]
                if A1[v][u] != A2[w][mp[u]]:
                    ok = False
                    break
            if ok:
                mp[v] = w
                used[w] = True
                if rec(i + 1):
                    return True
                used[w] = False
                mp[v] = -1
        return False
    return list(mp) if rec(0) else None


def iso_certified(n, E1, E2):
    """True iff E1 and E2 are isomorphic, with the witnessing permutation
    re-verified by mapping the edge set explicitly."""
    mp = iso_map(n, E1, E2)
    if mp is None:
        return False
    if sorted(mp) != list(range(n)):
        return False
    img = sorted((min(mp[a], mp[b]), max(mp[a], mp[b])) for a, b in E1)
    return img == sorted((min(a, b), max(a, b)) for a, b in E2)


def g6_encoding_is_injective(s):
    """graph6 decoding is injective on well-formed strings: no single-character
    change to s can leave the decoded graph unaltered.  This is what makes the
    re-encoding round trip meaningful -- it says the paper's string cannot be
    perturbed without perturbing the graph it names."""
    n, E = g6_decode(s)
    for i in range(len(s)):
        for c in range(63, 127):
            if chr(c) == s[i]:
                continue
            t = s[:i] + chr(c) + s[i + 1:]
            try:
                n2, E2 = g6_decode(t)
            except ValueError:
                continue                      # rejected outright: fine
            if n2 == n and E2 == E:
                return False                  # two strings, one graph
    return True


# ----------------------------------------------------------------------------
# the class T(n) = connected graphs with n vertices, n+2 edges, 2 <= deg <= 3.
# Suppressing the degree-2 vertices of such a graph gives a connected multigraph
# on exactly 4 vertices in which every degree is 3 (loops counting twice); the
# graph is recovered by subdividing the 6 multi-edges.  P and Q both lie in T(20).
# ----------------------------------------------------------------------------
SLOTS = [(i, j) for i in range(4) for j in range(i, 4)]      # 10 slots, loops first


def enum_base_multigraphs():
    """exhaustive: every multiplicity vector over the 10 slots with all four
    degrees equal to 3 (a loop contributing 2) and connected support."""
    out = []

    def conn_ok(mult):
        par = list(range(4))

        def f(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for idx, (i, j) in enumerate(SLOTS):
            if mult[idx] > 0 and i != j:
                a, b = f(i), f(j)
                if a != b:
                    par[a] = b
        return len({f(i) for i in range(4)}) == 1

    def rec(idx, mult, deg):
        if idx == len(SLOTS):
            if all(d == 3 for d in deg) and conn_ok(mult):
                out.append(tuple(mult))
            return
        i, j = SLOTS[idx]
        for mm in range(0, 4):
            nd = deg[:]
            if i == j:
                nd[i] += 2 * mm
            else:
                nd[i] += mm
                nd[j] += mm
            if max(nd) > 3:
                break
            rec(idx + 1, mult + [mm], nd)
    rec(0, [], [0, 0, 0, 0])
    return out


def base_representatives(all_mults):
    """one labelled representative per isomorphism class of base multigraph."""
    seen, reps = set(), []
    for m in all_mults:
        c = canon_mult(m)
        if c not in seen:
            seen.add(c)
            reps.append(m)
    return reps


def base_edge_list(mult):
    E = []
    for idx, (i, j) in enumerate(SLOTS):
        E.extend([(i, j)] * mult[idx])
    return E


def subdivide(base_edges, tvec, n_target):
    """subdivide the k-th base edge with tvec[k] internal vertices; returns
    (n, sorted edges) or None if the result is not a simple graph of the
    required order/size."""
    nxt = 4
    edges = []
    for (i, j), t in zip(base_edges, tvec):
        if t == 0:
            if i == j:
                return None                  # an unsubdivided loop
            edges.append((i, j))
        else:
            path = [i] + list(range(nxt, nxt + t)) + [j]
            nxt += t
            edges.extend((path[a], path[a + 1]) for a in range(len(path) - 1))
    S = set()
    for a, b in edges:
        if a == b:
            return None
        k = (min(a, b), max(a, b))
        if k in S:
            return None                      # parallel edges left unsubdivided
        S.add(k)
    if nxt != n_target or len(S) != n_target + 2:
        return None
    return nxt, sorted(S)


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for f in range(total + 1):
        for rest in compositions(total - f, parts - 1):
            yield (f,) + rest


_PERM4 = list(permutations(range(4)))
_SLOTMAP = []                                  # slot index -> image slot index
for _p in _PERM4:
    _SLOTMAP.append([SLOTS.index((min(_p[i], _p[j]), max(_p[i], _p[j])))
                     for (i, j) in SLOTS])


def suppress(n, edges):
    """Suppress the degree-2 vertices.  Returns (mult_vector, slot->sorted tuple
    of subdivision counts) with the four degree-3 vertices taken in increasing
    order, or None if the graph is not in the class T(n)."""
    adj = adjacency(n, edges)
    deg = [len(a) for a in adj]
    if any(d not in (2, 3) for d in deg):
        return None
    br = [v for v in range(n) if deg[v] == 3]
    if len(br) != 4 or not is_connected(n, adj):
        return None
    idx = {v: i for i, v in enumerate(br)}
    seen, res = set(), dict((s, []) for s in SLOTS)
    for v in br:
        for w in adj[v]:
            prev, cur, t = v, w, 0
            walk = [(min(v, w), max(v, w))]
            while deg[cur] == 2:
                t += 1
                nx = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
                walk.append((min(cur, nx), max(cur, nx)))
                prev, cur = cur, nx
                if t > n:
                    return None
            key = tuple(sorted(walk))
            if key in seen:
                continue
            seen.add(key)
            a, b = idx[v], idx[cur]
            res[(min(a, b), max(a, b))].append(t)
    if sum(len(res[s]) for s in SLOTS) != 6:
        return None
    for s in SLOTS:
        res[s] = tuple(sorted(res[s]))
    return tuple(len(res[s]) for s in SLOTS), res


def signature(res):
    """Canonical form of a suppressed graph: the lexicographic minimum over the
    24 relabellings of the four branch vertices of the tuple of sorted
    subdivision-count lists per slot.  Two graphs of T(n) are isomorphic if and
    only if their signatures agree, i.e. this is claimed to be a complete
    invariant on T(n).  That claim is NOT taken on trust: check E7 validates it
    exhaustively at n = 10 in both directions against a generic backtracking
    isomorphism search, and check F2b validates the direction that matters for
    the maximisation (no merging of non-isomorphic graphs) on all 58461
    configurations of T(20) itself."""
    keys = [res[s] for s in SLOTS]
    best = None
    for mp in _SLOTMAP:
        cand = [None] * 10
        for k in range(10):
            cand[mp[k]] = keys[k]
        cand = tuple(cand)
        if best is None or cand < best:
            best = cand
    return best


def canon_mult(mult):
    best = None
    for mp in _SLOTMAP:
        cand = [0] * 10
        for k in range(10):
            cand[mp[k]] = mult[k]
        cand = tuple(cand)
        if best is None or cand < best:
            best = cand
    return best


# ----------------------------------------------------------------------------
# counting isomorphism classes exactly: Burnside/Polya on the pair group, then
# the multiset logarithm 1+G(x,y) = exp(sum_j C(x^j,y^j)/j).
# ----------------------------------------------------------------------------
def integer_partitions(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxp), 0, -1):
        for rest in integer_partitions(n - p, p):
            yield (p,) + rest


def pair_orbit_lengths(part):
    """orbit lengths of the induced action on unordered vertex pairs of a
    permutation with the given cycle type."""
    L = []
    for c in part:
        if c % 2:
            L += [c] * ((c - 1) // 2)
        else:
            L += [c] * ((c - 2) // 2) + [c // 2]
    for i in range(len(part)):
        for j in range(i + 1, len(part)):
            a, b = part[i], part[j]
            g = gcd(a, b)
            L += [a * b // g] * g
    return L


def all_graph_counts(N, K):
    """g[n][k] = number of graphs (up to isomorphism) with n vertices, k edges,
    for 1<=n<=N, 0<=k<=K.  Exact, by Burnside."""
    g = {}
    for n in range(1, N + 1):
        acc = [Fraction(0)] * (K + 1)
        for part in integer_partitions(n):
            cnt = {}
            for c in part:
                cnt[c] = cnt.get(c, 0) + 1
            z = 1
            for c, a in cnt.items():
                z *= (c ** a) * _fact(a)
            poly = [0] * (K + 1)
            poly[0] = 1
            for L in pair_orbit_lengths(part):
                if L <= K:
                    for t in range(K, L - 1, -1):
                        if poly[t - L]:
                            poly[t] += poly[t - L]
            for k in range(K + 1):
                if poly[k]:
                    acc[k] += Fraction(poly[k], z)
        row = []
        for k in range(K + 1):
            if acc[k].denominator != 1:
                raise ArithmeticError("Burnside average is not an integer")
            row.append(int(acc[k]))
        g[n] = row
    return g


def _fact(k):
    r = 1
    for i in range(2, k + 1):
        r *= i
    return r


def connected_graph_counts(N, K):
    """c[n][k] = number of *connected* graphs up to isomorphism with n vertices
    and k edges.  A graph is a multiset of connected graphs, so
        1 + G(x,y) = exp( sum_{j>=1} C(x^j,y^j)/j );
    take the logarithm of the Burnside table and undo the j-fold substitutions."""
    g = all_graph_counts(N, K)

    def zeros():
        return [[Fraction(0)] * (K + 1) for _ in range(N + 1)]

    G = zeros()
    for n in range(1, N + 1):
        for k in range(K + 1):
            G[n][k] = Fraction(g[n][k])

    def mul(A, B):
        R = zeros()
        for n1 in range(1, N + 1):
            for k1 in range(K + 1):
                a = A[n1][k1]
                if not a:
                    continue
                for n2 in range(1, N + 1 - n1):
                    for k2 in range(K + 1 - k1):
                        b = B[n2][k2]
                        if b:
                            R[n1 + n2][k1 + k2] += a * b
        return R

    Lg, Pw = zeros(), G                      # log(1+G) = sum (-1)^{m+1} G^m/m
    for m in range(1, N + 1):
        s = Fraction((-1) ** (m + 1), m)
        for n in range(1, N + 1):
            for k in range(K + 1):
                if Pw[n][k]:
                    Lg[n][k] += s * Pw[n][k]
        if m < N:
            Pw = mul(Pw, G)
    C = zeros()
    for n in range(1, N + 1):
        for k in range(K + 1):
            v = Lg[n][k]
            for j in range(2, n + 1):
                if n % j == 0 and k % j == 0:
                    v -= Fraction(1, j) * C[n // j][k // j]
            C[n][k] = v
    out = {}
    for n in range(1, N + 1):
        row = []
        for k in range(K + 1):
            if C[n][k].denominator != 1:
                raise ArithmeticError("connected count is not an integer")
            row.append(int(C[n][k]))
        out[n] = row
    return g, out


def brute_iso_classes(n, k=None, connected_only=False):
    """Brute-force count of isomorphism classes of graphs on n labelled
    vertices (all k, or one fixed k) by minimising over all n! relabellings.
    Independent of the Burnside/logarithm machinery.  Only for tiny n."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    P = len(pairs)
    pidx = dict((p, i) for i, p in enumerate(pairs))
    maps = []
    for p in permutations(range(n)):
        maps.append([pidx[(min(p[a], p[b]), max(p[a], p[b]))] for (a, b) in pairs])
    counts = {}
    seen = set()
    if k is None:
        it = range(1 << P)
    else:
        it = (sum(1 << i for i in c) for c in combinations(range(P), k))
    for mask in it:
        best = None
        for mp in maps:
            v = 0
            for i in range(P):
                if mask >> i & 1:
                    v |= 1 << mp[i]
            if best is None or v < best:
                best = v
        if best in seen:
            continue
        seen.add(best)
        edges = [pairs[i] for i in range(P) if mask >> i & 1]
        if connected_only and not is_connected(n, adjacency(n, edges)):
            continue
        counts[len(edges)] = counts.get(len(edges), 0) + 1
    return counts


# ----------------------------------------------------------------------------
# comparing rigorous enclosures with the paper's decimal strings
# ----------------------------------------------------------------------------
def _dec_parts(s):
    neg = s.startswith("-")
    if neg:
        s = s[1:]
    ip, _, fp = s.partition(".")
    v = Fraction(int(ip + fp), 10 ** len(fp))
    return (-v if neg else v), len(fp)


def brackets_truncated(iv, s):
    """the paper prints s followed by '...': every value in our enclosure must
    truncate to exactly the digits of s."""
    v0, d = _dec_parts(s)
    ulp = Fraction(1, 10 ** d)
    lo, hi = iv
    return v0 <= lo and hi < v0 + ulp


def brackets_rounded(iv, s):
    """s is a correctly rounded value: our enclosure must round to it."""
    v0, d = _dec_parts(s)
    half = Fraction(1, 2 * 10 ** d)
    lo, hi = iv
    return v0 - half < lo and hi < v0 + half


def approx_to(iv, v0, tol):
    lo, hi = iv
    return (v0 - tol) < lo and hi < (v0 + tol)


def dstr(x, k=22):
    """truncate a nonnegative Fraction to k decimals, as a string."""
    neg = x < 0
    if neg:
        x = -x
    m = x.numerator * 10 ** k // x.denominator
    s = str(m).rjust(k + 1, "0")
    return ("-" if neg else "") + s[:-k] + "." + s[-k:]


def newton_power_sums(c, K):
    """power sums of the roots of x^n + c[1]x^{n-1} + ... + c[n]."""
    n = len(c) - 1
    s = [n] + [0] * K
    for k in range(1, K + 1):
        acc = k * c[k] if k <= n else 0
        for i in range(1, k):
            if i <= n:
                acc += c[i] * s[k - i]
        s[k] = -acc
    return s


# ============================================================================
# SECTION A -- the exhibited object
# ============================================================================
def section_A():
    print("--- A. the exhibited graphs ---")
    words = []                                     # P rebuilt from the paper's words
    for cyc in PAPER_CYCLES:
        k = len(cyc)
        for i in range(k):
            a, b = cyc[i], cyc[(i + 1) % k]
            words.append((min(a, b), max(a, b)))
    words += [(min(a, b), max(a, b)) for a, b in PAPER_EXTRA_EDGES]
    words = sorted(set(words))

    n, E = g6_decode(PAPER_G6_P)
    check("A1.g6_P_wellformed_20_vertices_22_edges",
          n == PAPER_N and len(E) == PAPER_M, "n=%d m=%d" % (n, len(E)))
    # the round trip alone cannot fail (g6_decode rejects nonzero padding, so
    # re-encoding is forced); the content is the injectivity statement, which
    # says no single-character change to the input string names this graph.
    # (The string is this program's own, not a quotation: see the header.  What
    # ties P to the paper is A3, which matches the decoded edge set to the
    # paper's explicit model.)
    check("A2.g6_P_reencodes_to_its_own_string_and_encoding_is_injective",
          g6_encode(n, E) == PAPER_G6_P and g6_encoding_is_injective(PAPER_G6_P),
          "no 1-character variant of the string decodes to the same graph")
    check("A3.g6_P_is_the_paper_explicit_model", E == words,
          "3 six-cycles + %s" % (PAPER_EXTRA_EDGES,))
    adjP = adjacency(n, E)
    degP = sorted(len(a) for a in adjP)
    check("A4.P_connected", is_connected(n, adjP))
    rankP, compP = cycle_rank_by_spanning_forest(n, E)
    check("A5.P_tricyclic_m_eq_n_plus_2_cycle_rank_3",
          len(E) == n + 2 and compP == 1 and rankP == 3,
          "spanning-forest rank=%d components=%d" % (rankP, compP))
    check("A6.P_degree_sequence_16x2_4x3", degP == [2] * 16 + [3] * 4,
          "min=%d max=%d" % (degP[0], degP[-1]))
    check("A7.P_bipartite", is_bipartite(n, adjP))
    # the ACH construction: three disjoint 6-cycles and a path P_{n-18}=P_2,
    # two cycles joined to one end of the path, the third to the other end
    ok_cycles = all(len(c) == 6 and is_cycle_on(c, induced_edges(E, c))
                    for c in PAPER_CYCLES)
    cyc_vs = [v for c in PAPER_CYCLES for v in c]
    disjoint = len(set(cyc_vs)) == 18
    rest = sorted(set(range(n)) - set(cyc_vs))
    path_edge = (min(rest), max(rest)) in set(E) if len(rest) == 2 else False
    attach = {}
    for a, b in PAPER_EXTRA_EDGES:
        for u, v in ((a, b), (b, a)):
            if u in rest and v in cyc_vs:
                attach.setdefault(u, []).append(v)
    pattern = sorted(len(x) for x in attach.values()) == [1, 2]
    ends = len(rest) == 2 and set(attach) == set(rest)
    which = sorted(sorted(v) for v in attach.values())
    # "two cycles joined to one endpoint and the third cycle to the other" also
    # requires the three attachment edges to reach three DIFFERENT cycles;
    # counting attachments per endpoint does not say that on its own.
    cyc_of = {}
    for ci, c in enumerate(PAPER_CYCLES):
        for v in c:
            cyc_of[v] = ci
    hit = [sorted(cyc_of[v] for v in vs) for vs in attach.values()]
    three_cycles = (sorted(i for h in hit for i in h)
                    == list(range(len(PAPER_CYCLES)))
                    and all(len(set(h)) == len(h) for h in hit))
    check("A8.P_is_the_ACH_family_member_three_C6_plus_path_P2",
          ok_cycles and disjoint and path_edge and pattern and ends
          and three_cycles,
          "path=%s attachments=%s to distinct cycles %s" % (rest, which, hit))
    artP = articulation_vertices(n, adjP)
    check("A9.P_articulation_vertices_eq_%d" % PAPER_ART_P,
          len(artP) == PAPER_ART_P, "cut vertices %s" % (artP,))
    return n, E, adjP


def in_class_T(n, edges):
    """suppression data + signature, or (None,None) if the graph is not in T(n)."""
    r = suppress(n, edges)
    if r is None:
        return None, None
    return r[0], signature(r[1])


def same_class(s1, s2):
    """isomorphism test via the invariant validated by checks E7/F2b.  A missing
    invariant (a graph outside T) is never equal to anything -- in particular
    None != None."""
    return s1 is not None and s2 is not None and s1 == s2


def section_B(nP, EP):
    print("--- B. the runner-up Q and the isomorphism claims ---")
    nQ, EQ = g6_decode(PAPER_G6_Q)
    d = (min(PAPER_Q_DEL), max(PAPER_Q_DEL))
    a = (min(PAPER_Q_ADD), max(PAPER_Q_ADD))
    built = sorted((set(EP) - set([d])) | set([a]))
    check("B1.g6_Q_wellformed_20_vertices_22_edges",
          nQ == PAPER_N and len(EQ) == PAPER_M, "n=%d m=%d" % (nQ, len(EQ)))
    check("B2.g6_Q_equals_P_minus_%s_plus_%s" % (d, a), EQ == built)
    adjQ = adjacency(nQ, EQ)
    degQ = sorted(len(x) for x in adjQ)
    check("B3.Q_connected_and_tricyclic",
          is_connected(nQ, adjQ) and len(EQ) == nQ + 2)
    check("B4.Q_degree_sequence_16x2_4x3", degQ == [2] * 16 + [3] * 4)
    artQ = articulation_vertices(nQ, adjQ)
    check("B5.Q_articulation_vertices_eq_%d" % PAPER_ART_Q,
          len(artQ) == PAPER_ART_Q, "cut vertices %s" % (artQ,))
    mP, sP = in_class_T(nP, EP)
    mQ, sQ = in_class_T(nQ, EQ)
    check("B6.P_and_Q_lie_in_the_class_T20_deg_2_or_3",
          sP is not None and sQ is not None)
    # B7 is the paper's isomorphism claim (P and Q are distinct); B9 is a
    # relabelling test vector of this program's own.  Each is settled twice: by
    # the cheap invariant AND by the generic backtracking search, which shares
    # no reasoning with it.  Sections E/F validate the invariant.
    iso_PQ = iso_certified(PAPER_N, EP, EQ)
    check("B7.Q_is_not_isomorphic_to_P", sP is not None and sQ is not None
          and not same_class(sP, sQ) and not iso_PQ,
          "distinct invariants AND exhaustive backtracking search found no "
          "isomorphism")
    # B8/B9/B10 concern G6_P_ALT_LABELLING, which the paper does not print and
    # whose claimed labelg provenance is not verifiable without nauty.  They
    # therefore assert nothing about canonicity: only well-formedness,
    # isomorphism to P, and graph6 round-tripping/injectivity on the string
    # itself.  See the header block.
    nC, EC = g6_decode(G6_P_ALT_LABELLING)
    adjC = adjacency(nC, EC)
    mC, sC = in_class_T(nC, EC)
    check("B8.alt_labelling_g6_wellformed_connected_20_vertices_22_edges",
          nC == PAPER_N and len(EC) == PAPER_M and is_connected(nC, adjC))
    iso_PC = iso_certified(PAPER_N, EP, EC)
    check("B9.alt_labelling_g6_is_isomorphic_to_P",
          same_class(sC, sP) and iso_PC,
          "same signature AND explicit verified permutation; distinct "
          "labelling: %s" % (EC != EP,))
    check("B10.alt_labelling_g6_reencodes_to_its_own_string_and_encoding_"
          "is_injective",
          g6_encode(nC, EC) == G6_P_ALT_LABELLING
          and g6_encoding_is_injective(G6_P_ALT_LABELLING),
          "canonicity NOT asserted: no nauty here")
    return nQ, EQ, adjQ, sP, sQ


def expand_factors(facs):
    r = [1]
    for f, m in facs:
        r = pmul(r, ppow(f, m))
    return r


def spectral_moment_ok(n, adj, edges, cp, m):
    """Low-order coefficients of cp against *independently enumerated*
    subgraph counts (Sachs' coefficient theorem), plus the traces.

    Warning about what is and is not information here.  Faddeev-LeVerrier
    derives cp from the traces by Newton's identities, so "Newton power sums of
    cp == tr(A^k)" holds for any matrix whatever, correct or not; likewise
    cp[2] == -m and cp[3] == -2*(triangles) are identities once the power sums
    agree, and tr(A^4) - 2m - 4*sum C(d,2) is a nonnegative multiple of 8 for
    every graph.  Those clauses only exercise the implementations.  The parts
    that can genuinely fail are the comparisons against counts obtained without
    any matrix algebra: |E|, the 4-cycles counted from common neighbourhoods,
    and the pairs of disjoint edges, via
        cp[2] = -|E|,  cp[3] = -2*(triangles),
        cp[4] = #(pairs of disjoint edges) - 2*#(4-cycles).
    The whole polynomial is confirmed separately by charpoly_confirm."""
    ps = newton_power_sums(cp, 4)
    tr = [n] + [trace_power(n, adj, k) for k in range(1, 5)]
    c2 = sum(len(a) * (len(a) - 1) // 2 for a in adj)     # paths with two edges
    ok = all(ps[k] == tr[k] for k in range(1, 5))         # implementation only
    ok = ok and tr[1] == 0 and tr[2] == 2 * m and tr[3] % 6 == 0
    tri = tr[3] // 6                                     # triangles
    r4 = tr[4] - 2 * m - 4 * c2                          # = 8 * (4-cycles)
    if r4 % 8 or r4 < 0:                                 # identity; assert it
        raise ArithmeticError("tr(A^4) decomposition failed")
    c4_direct = count_c4_direct(n, adj)                  # no matrix powers
    dpairs = count_disjoint_edge_pairs(edges)
    ok = ok and len(edges) == m and r4 // 8 == c4_direct
    ok = ok and cp[2] == -m and cp[3] == -2 * tri
    ok = ok and cp[4] == dpairs - 2 * c4_direct          # Sachs, x^{n-4}
    return ok, (tri, c4_direct)


def section_C(nP, EP, adjP, nQ, EQ, adjQ):
    print("--- C. characteristic polynomials and rigorous energies ---")
    cpP = charpoly(nP, adjP)
    cpQ = charpoly(nQ, adjQ)
    check("C1.charpoly_P_from_matrix_equals_paper_eq1",
          cpP == expand_factors(PAPER_PHI_P), "deg=%d" % (len(cpP) - 1))
    check("C1b.charpoly_P_confirmed_by_independent_determinant_engine",
          charpoly_confirm(nP, adjP, cpP),
          "det(tI-A) by Bareiss at t=-1..%d agrees at all %d points"
          % (nP - 1, nP + 1))
    check("C2.charpoly_Q_from_matrix_equals_paper_eq2",
          cpQ == expand_factors(PAPER_PHI_Q))
    check("C2b.charpoly_Q_confirmed_by_independent_determinant_engine",
          charpoly_confirm(nQ, adjQ, cpQ))
    okP, dP = spectral_moment_ok(nP, adjP, EP, cpP, PAPER_M)
    okQ, dQ = spectral_moment_ok(nQ, adjQ, EQ, cpQ, PAPER_M)
    check("C3.charpoly_P_consistent_with_enumerated_subgraph_counts", okP,
          "triangles=%d 4-cycles=%d" % dP)
    check("C4.charpoly_Q_consistent_with_enumerated_subgraph_counts", okQ,
          "triangles=%d 4-cycles=%d" % dQ)
    check("C5.charpolys_of_P_and_Q_differ", cpP != cpQ,
          "so P and Q are not cospectral")
    eps = Fraction(1, 10 ** 28)
    ivP = energy_interval(cpP, eps)
    ivQ = energy_interval(cpQ, eps)
    check("C6.energy_of_P_enclosed_rigorously", ivP is not None,
          "" if ivP is None else "[%s, %s]" % (dstr(ivP[0]), dstr(ivP[1])))
    check("C7.energy_of_Q_enclosed_rigorously", ivQ is not None,
          "" if ivQ is None else "[%s, %s]" % (dstr(ivQ[0]), dstr(ivQ[1])))
    if ivP is None or ivQ is None:
        return None, None
    check("C8.energy_of_P_equals_paper_digits_%s" % PAPER_E_P,
          brackets_truncated(ivP, PAPER_E_P), dstr(ivP[0], 25))
    check("C9.energy_of_Q_equals_paper_digits_%s" % PAPER_E_Q,
          brackets_truncated(ivQ, PAPER_E_Q), dstr(ivQ[0], 25))
    check("C10.energy_of_P_rounds_to_result1_value_%s" % PAPER_E_P_ROUND,
          brackets_rounded(ivP, PAPER_E_P_ROUND))
    check("C11.energy_of_Q_rounds_to_result1_value_%s" % PAPER_E_Q_ROUND,
          brackets_rounded(ivQ, PAPER_E_Q_ROUND))
    gap = (ivP[0] - ivQ[1], ivP[1] - ivQ[0])
    check("C12.energy_of_P_strictly_exceeds_energy_of_Q", ivP[0] > ivQ[1],
          "gap >= %s" % dstr(gap[0], 22))
    check("C13.gap_equals_paper_digits_%s" % PAPER_GAP,
          brackets_truncated(gap, PAPER_GAP), dstr(gap[0], 25))
    check("C14.gap_matches_abstract_1.3061e-2",
          approx_to(gap, PAPER_SEP_5SIG, Fraction(5, 10 ** 7)))
    check("C15.gap_matches_status_section_1.31e-2",
          approx_to(gap, PAPER_SEP_3SIG, Fraction(5, 10 ** 5)))
    check("C16.gap_dwarfs_the_reported_solver_discrepancy_5.15e-13",
          gap[0] > Fraction(515, 10 ** 15) * 10 ** 6,
          "gap > 1e6 x 5.15e-13")
    return ivP, ivQ


# classical published totals, used only to validate the counting machinery:
# A000088 (graphs on n nodes) and A001349 (connected graphs on n nodes).
A000088 = {1: 1, 2: 2, 3: 4, 4: 11, 5: 34, 6: 156, 7: 1044}
A001349 = {1: 1, 2: 1, 3: 2, 4: 6, 5: 21, 6: 112, 7: 853}


def section_D():
    print("--- D. the census total, re-derived without geng ---")
    t = time.time()
    g, c = connected_graph_counts(PAPER_N, PAPER_M)
    check("D1.number_of_connected_20_vertex_22_edge_classes_eq_%d"
          % PAPER_CENSUS_CLASSES,
          c[PAPER_N][PAPER_M] == PAPER_CENSUS_CLASSES,
          "Burnside+multiset-log gives %d in %.1f s"
          % (c[PAPER_N][PAPER_M], time.time() - t))
    tri = [c[k][k + 2] for k in range(4, PAPER_N + 1) if k + 2 <= PAPER_M]
    # monotonicity alone is nearly unfalsifiable, so anchor the first three
    # entries of the tricyclic sequence to brute-force isomorph counting.
    anchor = [brute_iso_classes(4, 6, connected_only=True).get(6, 0),
              brute_iso_classes(5, 7, connected_only=True).get(7, 0),
              brute_iso_classes(6, 8, connected_only=True).get(8, 0)]
    check("D2.tricyclic_class_counts_anchored_at_n456_and_increasing",
          tri[:3] == anchor
          and all(x > 0 for x in tri) and all(tri[i] < tri[i + 1]
                                              for i in range(len(tri) - 1)),
          "brute force n=4,5,6 gives %s; n=4..20: %s"
          % (anchor, tri[:6] + ["..."] + tri[-1:]))
    ok88 = all(sum(g[n][k] for k in range(PAPER_M + 1)) == A000088[n]
               for n in range(1, 8))
    ok49 = all(sum(c[n][k] for k in range(PAPER_M + 1)) == A001349[n]
               for n in range(1, 8))
    check("D3.all_graph_totals_match_A000088_for_n_le_7", ok88,
          "sums %s" % ([sum(g[n][k] for k in range(PAPER_M + 1))
                        for n in range(1, 8)],))
    check("D4.connected_totals_match_A001349_for_n_le_7", ok49,
          "sums %s" % ([sum(c[n][k] for k in range(PAPER_M + 1))
                        for n in range(1, 8)],))
    b5 = brute_iso_classes(5)
    b5c = brute_iso_classes(5, connected_only=True)
    check("D5.brute_force_reproduces_all_graph_counts_at_n5",
          all(b5.get(k, 0) == g[5][k] for k in range(11)),
          "g(5,.)=%s" % (g[5][:11],))
    check("D6.brute_force_reproduces_connected_counts_at_n5",
          all(b5c.get(k, 0) == c[5][k] for k in range(11)),
          "c(5,.)=%s" % (c[5][:11],))
    b68 = brute_iso_classes(6, 8, connected_only=True)
    check("D7.brute_force_reproduces_connected_count_at_n6_m8",
          b68.get(8, 0) == c[6][8], "brute=%d table=%d" % (b68.get(8, 0), c[6][8]))
    return c


def enum_T(n_target, reps):
    """Every graph of T(n_target) arises by subdividing one of the base
    multigraphs; all compositions of n_target-4 into the 6 base edges are tried,
    so this enumeration is exhaustive by construction (validated at n=10 below).
    Yields (edges, mult, tvec)."""
    for mult in reps:
        be = base_edge_list(mult)
        for tv in compositions(n_target - 4, 6):
            r = subdivide(be, tv, n_target)
            if r is not None:
                yield r[1], mult, tv


def enum_degseq(n, target):
    """all labelled graphs on n vertices with the exact degree sequence given."""
    deg = [0] * n
    edges = []
    out = []

    def rec(v):
        if v == n:
            out.append(tuple(edges))
            return
        need = target[v] - deg[v]
        if need < 0:
            return
        cands = [w for w in range(v + 1, n) if deg[w] < target[w]]
        if need > len(cands):
            return
        for combo in combinations(cands, need):
            for w in combo:
                deg[w] += 1
                edges.append((v, w))
            deg[v] += need
            rec(v + 1)
            deg[v] -= need
            for w in combo:
                deg[w] -= 1
            if need:
                del edges[len(edges) - need:]
    rec(0)
    return out


def section_E(reps):
    """Completeness of the T-enumeration, cross-checked at n=10 -- the smallest
    order at which all five base multigraphs occur -- against an independent
    exhaustive enumeration.  Every graph of T(10) is isomorphic to one whose
    four degree-3 vertices are labelled 0,1,2,3, so enumerating all labelled
    graphs with degree sequence (3,3,3,3,2,2,2,2,2,2) meets every isomorphism
    class of T(10); each one is then required to suppress into the base family,
    and the two isomorphism-class sets are required to coincide."""
    print("--- E. the subclass enumerator, validated at n=10 ---")
    allm = enum_base_multigraphs()
    # enum_base_multigraphs only emits vectors it has already forced to be cubic
    # and connected, so re-testing its output proves nothing.  Re-derive the
    # family by a different route instead: every multiset of 6 slots, filtered
    # afterwards, with connectivity decided by the generic component counter.
    indep = set()
    for combo in combinations_with_replacement(range(len(SLOTS)), 6):
        d = [0, 0, 0, 0]
        mult = [0] * len(SLOTS)
        for idx in combo:
            mult[idx] += 1
            i, j = SLOTS[idx]
            if i == j:
                d[i] += 2
            else:
                d[i] += 1
                d[j] += 1
        if d != [3, 3, 3, 3]:
            continue
        e = [(i, j) for idx, (i, j) in enumerate(SLOTS)
             if mult[idx] and i != j]
        if n_components(4, adjacency(4, e)) != 1:
            continue
        indep.add(tuple(mult))
    check("E1.base_multigraph_family_reproduced_by_an_independent_enumeration",
          indep == set(allm) and len(allm) > 0,
          "%d labelled from each of the two enumerations, %d isomorphism "
          "classes" % (len(allm), len(reps)))
    # canon_mult was used to build reps, so re-applying it cannot fail; test the
    # 24 relabellings directly instead.
    pair_ok = True
    for x in range(len(reps)):
        for y in range(x + 1, len(reps)):
            for mp in _SLOTMAP:
                img = [0] * 10
                for k in range(10):
                    img[mp[k]] = reps[x][k]
                if tuple(img) == tuple(reps[y]):
                    pair_ok = False
    check("E2.base_representatives_are_pairwise_non_isomorphic", pair_ok,
          "no relabelling of the 4 branch vertices carries one onto another")
    nc = len(list(compositions(PAPER_N - 4, 6)))
    binom = _fact(PAPER_N - 4 + 5) // (_fact(5) * _fact(PAPER_N - 4))
    check("E3.all_compositions_of_16_into_6_parts_are_tried", nc == binom,
          "%d = C(21,5)" % nc)
    t = time.time()
    sig_sub = set()
    for edges, mult, tv in enum_T(10, reps):
        _, s = in_class_T(10, edges)
        if s is None:
            sig_sub = None
            break
        sig_sub.add(s)
    check("E4.subdivided_graphs_at_n10_all_lie_in_T10", sig_sub is not None,
          "%d isomorphism classes" % (0 if sig_sub is None else len(sig_sub)))
    brute = enum_degseq(10, [3, 3, 3, 3] + [2] * 6)
    basefam = {canon_mult(m) for m in reps}
    sig_bru, bad, conn = set(), 0, 0
    for E in brute:
        adj = adjacency(10, E)
        if not is_connected(10, adj):
            continue
        conn += 1
        r = suppress(10, E)
        if r is None or canon_mult(r[0]) not in basefam \
           or sum(sum(v) for v in r[1].values()) != 6:
            bad += 1
            continue
        sig_bru.add(signature(r[1]))
    check("E5.every_T10_graph_suppresses_into_the_base_family", bad == 0,
          "%d labelled graphs, %d connected, %d not accounted for"
          % (len(brute), conn, bad))
    check("E6.T10_enumerator_is_complete_and_sound",
          sig_sub is not None and sig_sub == sig_bru,
          "%d classes both ways, %.1f s"
          % (len(sig_bru), time.time() - t))
    # E6 applies the SAME invariant to both sides, so it would still pass if the
    # invariant were too coarse (both sets would collapse alike).  Completeness
    # of the invariant is load-bearing -- everything downstream keeps one
    # representative per signature -- so validate it here against the generic
    # backtracking isomorphism search, in both directions and exhaustively.
    # n = 10 is the smallest order at which all five base multigraphs occur.
    t = time.time()
    members, first = {}, {}
    for edges, mult, tv in enum_T(10, reps):
        r = suppress(10, edges)
        s = None if r is None else signature(r[1])
        members.setdefault(s, []).append(edges)
        first.setdefault(s, edges)
    merged = 0                     # same signature but NOT isomorphic
    for s, mem in members.items():
        rep = first[s]
        for E in mem:
            if E is not rep and not iso_certified(10, rep, E):
                merged += 1
    ks = list(first)
    split = 0                      # different signatures but isomorphic
    for x in range(len(ks)):
        for y in range(x + 1, len(ks)):
            if iso_certified(10, first[ks[x]], first[ks[y]]):
                split += 1
    check("E7.signature_is_a_complete_isomorphism_invariant_on_T10",
          merged == 0 and split == 0 and None not in first and len(ks) > 1,
          "%d classes: %d same-signature pairs non-isomorphic, %d "
          "distinct-signature pairs isomorphic, %d pairs tested, %.1f s"
          % (len(ks), merged, split, len(ks) * (len(ks) - 1) // 2,
             time.time() - t))


def expected_suppression(mult, tvec):
    """the (slot -> sorted counts) data that subdividing (mult,tvec) must yield."""
    res = dict((s, []) for s in SLOTS)
    for (i, j), t in zip(base_edge_list(mult), tvec):
        res[(min(i, j), max(i, j))].append(t)
    for s in SLOTS:
        res[s] = tuple(sorted(res[s]))
    return res


def section_F(reps, sigP, sigQ, ivP, ivQ, total_classes):
    print("--- F. exhaustive energy scan over T(20) ---")
    t = time.time()
    classes = {}
    nconf, bad_struct, bad_round = 0, 0, 0
    merged = 0                     # same signature, NOT isomorphic
    ncomp = len(list(compositions(PAPER_N - 4, 6)))
    grid = len(reps) * ncomp       # every (base, composition) pair enum_T sees
    nskipped = 0                   # rejected by subdivide as not simple / wrong n
    for mult in reps:
        be = base_edge_list(mult)
        for tv in compositions(PAPER_N - 4, 6):
            r0 = subdivide(be, tv, PAPER_N)
            if r0 is None:
                nskipped += 1
                continue
            edges = r0[1]
            nconf += 1
            adj = adjacency(PAPER_N, edges)
            deg = [len(a) for a in adj]
            if (len(edges) != PAPER_M or min(deg) < 2 or max(deg) > 3
                    or not is_connected(PAPER_N, adj)):
                bad_struct += 1
                continue
            r = suppress(PAPER_N, edges)
            if r is None or signature(r[1]) != signature(
                    expected_suppression(mult, tv)):
                bad_round += 1
                continue
            sg = signature(r[1])
            if sg in classes:
                # only one representative per signature is ever scanned, so a
                # signature collision between non-isomorphic graphs would drop a
                # candidate maximiser silently.  Refute that, here, at n = 20,
                # with the generic backtracking search.
                if not iso_certified(PAPER_N, classes[sg], edges):
                    merged += 1
            else:
                classes[sg] = edges
    # the enumeration is only exhaustive if nothing fell off the grid
    check("F1.every_enumerated_graph_is_connected_20_22_with_deg_2_or_3",
          bad_struct == 0 and nconf > 0 and nconf + nskipped == grid,
          "%d subdivision configurations + %d non-simple = %d = %d bases x "
          "C(21,5); %d rejected"
          % (nconf, nskipped, nconf + nskipped, len(reps), bad_struct))
    check("F2.suppression_inverts_subdivision_on_every_configuration",
          bad_round == 0, "%d isomorphism classes in T(20)" % len(classes))
    check("F2b.no_two_non_isomorphic_T20_configurations_share_a_signature",
          merged == 0,
          "%d of the %d configurations collided with a stored representative; "
          "every collision verified isomorphic by explicit permutation"
          % (nconf - len(classes), nconf))
    check("F3.P_and_Q_are_both_found_by_the_T20_enumeration",
          sigP in classes and sigQ in classes)
    print("    scanning %d classes (%.1f s so far) ..." % (len(classes), time.time() - t),
          flush=True)
    eps = Fraction(1, 10 ** 6)
    energies, failed = {}, 0
    for sig, edges in classes.items():
        cp = charpoly(PAPER_N, adjacency(PAPER_N, edges))
        iv = energy_interval(cp, eps)
        if iv is None:
            failed += 1
            continue
        energies[sig] = (iv, cp)
    check("F4.energy_enclosed_rigorously_for_every_class_of_T20", failed == 0,
          "%d classes, %.1f s" % (len(energies), time.time() - t))
    # the scan's charpolys all come from one engine; confirm a sample of them,
    # and the two that carry the claim, with the independent determinant engine.
    sample = [s for i, s in enumerate(classes) if i % 400 == 0]
    for s in (sigP, sigQ):
        if s in classes and s not in sample:
            sample.append(s)
    cpbad = [s for s in sample
             if s not in energies
             or not charpoly_confirm(PAPER_N, adjacency(PAPER_N, classes[s]),
                                     energies[s][1])]
    check("F4b.sampled_T20_charpolys_confirmed_by_the_determinant_engine",
          not cpbad and len(sample) > 2,
          "%d of %d classes re-derived as det(tI-A), including P and Q"
          % (len(sample), len(classes)))
    atop = [s for s, (iv, cp) in energies.items() if iv[1] >= ivP[0]]
    pdat, qdat = energies.get(sigP), energies.get(sigQ)
    check("F5.unique_maximizer_of_energy_over_T20_is_P",
          sigP is not None and atop == [sigP] and pdat is not None
          and pdat[1] == charpoly(PAPER_N, adjacency(*g6_decode(PAPER_G6_P))),
          "%d of %d classes reach E(P); strict runner-up margin"
          % (len(atop), len(energies)))
    rest = [(iv, s) for s, (iv, cp) in energies.items() if not same_class(s, sigP)]
    second = [s for iv, s in rest if iv[1] >= ivQ[0]]
    top2 = max([iv[1] for iv, s in rest] or [None])
    third = sorted((iv[1] for iv, s in rest if s not in second), reverse=True)
    check("F6.runner_up_over_T20_is_Q",
          sigQ is not None and second == [sigQ] and qdat is not None
          and top2 == qdat[0][1],
          "next distinct energy <= %s" % (dstr(third[0], 9) if third else "n/a"))
    vP, _ = _dec_parts(PAPER_E_P)
    vQ, _ = _dec_parts(PAPER_E_Q)
    okF7 = pdat is not None and qdat is not None
    if okF7:
        okF7 = (pdat[0][0] <= vP <= pdat[0][1]) and (qdat[0][0] <= vQ <= qdat[0][1])
    check("F7.T20_scan_encloses_the_paper_top_two_values", okF7,
          "scan (1e-6 precision) top two: %s, %s"
          % (dstr(pdat[0][0], 9) if pdat else "n/a",
             dstr(qdat[0][0], 9) if qdat else "n/a"))
    frac = Fraction(len(energies), total_classes)
    print("    T(20) covers %d of %d census classes (%.3g of the census);"
          % (len(energies), total_classes, float(frac)))
    print("    completeness of this slice rests on the degree-2 suppression "
          "argument, which section E validates exhaustively at n=10.")
    return len(energies)


def section_G(EP, sigP, sigQ, ivP, ivQ):
    """the full one-edge-swap neighbourhood of P: P-e+f for every edge e and
    every non-edge f.  This leaves the class T(20) (degree 4 and degree 1
    vertices occur), so it tests P against a differently shaped slice."""
    print("--- G. exhaustive one-edge-swap neighbourhood of P ---")
    t = time.time()
    Ps = set(EP)
    pairs = [(i, j) for i in range(PAPER_N) for j in range(i + 1, PAPER_N)]
    groups, nconn, ndis, outside = {}, 0, 0, 0
    seen_edge_sets = set()
    shape_bad = 0                  # the count identity below is automatic, so
    for e in sorted(Ps):           # also test what each neighbour actually is
        for f in pairs:
            if f in Ps:
                continue
            edges = sorted((Ps - {e}) | {f})
            seen_edge_sets.add(tuple(edges))
            if (len(edges) != PAPER_M
                    or len(set(edges) ^ Ps) != 2
                    or e in set(edges) or f not in set(edges)):
                shape_bad += 1
            adj = adjacency(PAPER_N, edges)
            if not is_connected(PAPER_N, adj):
                ndis += 1
                continue
            nconn += 1
            _, s = in_class_T(PAPER_N, edges)
            if s is None:
                outside += 1
            cp = tuple(charpoly(PAPER_N, adj))
            groups.setdefault(cp, []).append(s)
    nswap = len(Ps) * (len(pairs) - len(Ps))
    check("G1.swap_neighbourhood_enumerated",
          nconn + ndis == nswap and nconn > 0 and shape_bad == 0
          and len(seen_edge_sets) == nswap,
          "%d distinct one-swap graphs (= 22 x 168), each with 22 edges and "
          "symmetric difference 2 from P; %d connected, %d disconnected, %d "
          "outside T(20), %d spectra"
          % (len(seen_edge_sets), nconn, ndis, outside, len(groups)))
    eps = Fraction(1, 10 ** 6)
    energies, failed = {}, 0
    for cp in groups:
        iv = energy_interval(list(cp), eps)
        if iv is None:
            failed += 1
        else:
            energies[cp] = iv
    check("G2.energy_enclosed_rigorously_for_every_neighbour", failed == 0,
          "%d distinct spectra, %.1f s" % (len(energies), time.time() - t))
    viol, copiesP = [], 0
    for cp, sigs in groups.items():
        if energies.get(cp) is None or energies[cp][1] >= ivP[0]:
            for s in sigs:
                if same_class(s, sigP):
                    copiesP += 1
                else:
                    viol.append(cp)
    check("G3.no_swap_neighbour_beats_P_except_copies_of_P", not viol and failed == 0,
          "%d of %d neighbours are isomorphic to P, none exceeds E(P)"
          % (copiesP, nconn))
    others = [(iv, cp) for cp, iv in energies.items()
              if not all(same_class(s, sigP) for s in groups[cp])]
    best = max([iv[1] for iv, cp in others] or [None])
    at_best = [cp for iv, cp in others if iv[1] >= ivQ[0]]
    okQ = (len(at_best) == 1
           and all(same_class(s, sigQ) for s in groups[at_best[0]])
           and energies[at_best[0]][0] <= _dec_parts(PAPER_E_Q)[0]
           <= energies[at_best[0]][1])
    check("G4.best_swap_neighbour_not_isomorphic_to_P_is_Q", okQ,
          "max non-P energy %s" % (dstr(best, 9) if best is not None else "n/a"))
    return nconn


NOT_REPRODUCED = """
NOT REPRODUCED (stated explicitly):
  * The paper's energy scan over all %d isomorphism classes is NOT re-run here.
    Regenerating those classes and computing about 3.3e9 twenty-by-twenty
    symmetric eigensolves is thousands of CPU-hours; this program instead
    (i) re-derives the census TOTAL %d exactly, by Burnside/Polya counting
    plus the multiset logarithm, with no use of geng, and (ii) re-runs the
    energy maximisation exhaustively over two slices that do fit the budget:
    the class T(20) of connected 20-vertex 22-edge graphs with all degrees in
    {2,3} (%s isomorphism classes, containing both P and Q), and the complete
    one-edge-swap neighbourhood of P (%s connected graphs, most of them outside
    T(20)).  Classes of the census outside those two slices are not rescanned.
  * The paper's shard bookkeeping (%d geng shards, and the confirmation run
    with 512 shards) is not reproducible from a single-process stdlib program,
    and no per-shard record counts are published to compare against.
  * The paper's binary64 figures - the counts 0 and 1 of records within 1e-7 of
    the maximum, and the 5.15e-13 agreement between two eigensolvers on 20000
    sampled graphs - concern that floating-point pipeline itself; this program
    replaces floating point by exact rational enclosures for the graphs it does
    examine, so those two numbers are corroborated in spirit, not re-measured.
  * The paper prints no graph6 string, so none of the three graph6 strings used
    here is a quotation from it.  P's and Q's are pinned to the paper's own
    words by checks A3 and B2 respectively; the third is an alternative
    labelling of P whose claimed nauty/labelg canonicity is NOT verified (nauty
    is not used), so checks B8, B9 and B10 exercise this program's graph6 and
    isomorphism code on a relabelled copy of P and corroborate nothing printed
    in the paper.
  * The conjecture for n >= 22 is outside the paper's scope and is not tested.
"""


def guarded(name, fn, *a):
    """run a section; an exception is reported as a failed check, never as a
    traceback, so the verdict line is always emitted."""
    try:
        return fn(*a)
    except Exception as exc:                                # noqa: BLE001
        check(name + ".section_raised_no_exception",
              False, "%s: %s" % (type(exc).__name__, exc))
        return None


def main():
    quick = "--quick" in sys.argv
    print("verify.py: verification of the n=20 tricyclic-energy census")
    print("python %s | exact integer/rational arithmetic only"
          % sys.version.split()[0])
    a = guarded("A", section_A)
    if a is None:
        finish()
    nP, EP, adjP = a
    b = guarded("B", section_B, nP, EP)
    if b is None:
        finish()
    nQ, EQ, adjQ, sigP, sigQ = b
    c3 = guarded("C", section_C, nP, EP, adjP, nQ, EQ, adjQ)
    if c3 is None or c3[0] is None or c3[1] is None:
        finish()
    ivP, ivQ = c3
    cc = guarded("D", section_D)
    if cc is None:
        finish()
    total = cc[PAPER_N][PAPER_M]
    reps = base_representatives(enum_base_multigraphs())
    nT, ng = "?", "?"
    if quick:
        skip("E.T10_completeness_cross_check", "--quick")
        skip("F.T20_exhaustive_energy_scan", "--quick")
    else:
        guarded("E", section_E, reps)
        nT = guarded("F", section_F, reps, sigP, sigQ, ivP, ivQ, total)
    ng = guarded("G", section_G, EP, sigP, sigQ, ivP, ivQ)
    print(NOT_REPRODUCED % (PAPER_CENSUS_CLASSES, total, nT, ng, PAPER_SHARDS))
    finish()


if __name__ == "__main__":
    main()

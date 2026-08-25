#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Verification program for "A Third-Derivative Counterexample for Chromatic
Polynomials".

Standard library only, exact integer / Fraction arithmetic throughout.  No
floating point is used for any decision.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs -- these are NOT checks by themselves)
--------------------------------------------------------------------------
  G6_H          the graph6 string of the exhibited graph H, printed in the
                paper's reproduction section; terminals are vertices 0 and 1.
  K_WITNESS=3   the derivative order of the witness.
  X_WITNESS=-1/2 the witness point.
  PAPER_P, PAPER_DP, PAPER_DDP, PAPER_DDDP
                the four rationals P^(j)(H,-1/2), j=0..3, displayed in the
                proof of Theorem 1.
  PAPER_L3      the claimed value of L_H^{(3)}(-1/2).
  PAPER_N=32, PAPER_M=34, PAPER_DELTA=4, PAPER_LENGTHS=(8,8,9,9)
  NY_CONST=3.01 the constant in Ning--Yang Theorem 1.5 (x <= -3.01*Delta*k).
G6_H, PAPER_N, PAPER_M, PAPER_DELTA, PAPER_LENGTHS, PAPER_TERMINALS, PAPER_P,
PAPER_DP, PAPER_DDP, PAPER_DDDP and PAPER_L3 are each compared against an
independently recomputed value.  K_WITNESS and X_WITNESS are not recomputed
and are not claims: they are the choice of derivative order and evaluation
point at which everything else is computed.  NY_CONST is NOT recomputed
either -- it is transcribed from the statement of Ning--Yang Theorem 1.5,
which this program cannot see; the Step 8 check
"witness_point_is_not_inside_the_ning_yang_proved_range" is exact GIVEN that
transcribed constant and is disclosed as such in Step 8 and in the closing
scope note.

--------------------------------------------------------------------------
WHAT IS DERIVED HERE (the checks)
--------------------------------------------------------------------------
  * graph6 is decoded from scratch, re-encoded, and the decoded graph is
    counted, degree-profiled, connectivity-tested and *structurally
    decomposed*; an explicit bijection onto an independently constructed
    Theta(8,8,9,9) is exhibited and its edge sets compared.  The target of
    that bijection is built from the length multiset the PAPER claims, not
    from the multiset recovered from the object, so the comparison can
    actually fail (it does, for instance, for Theta(7,9,9,9)).
  * planarity is proved twice: by a Kuratowski branch-vertex count and by
    tracing the faces of an explicit rotation system and checking Euler.
  * P(H,q) is computed by TWO independent algorithms: (i) the paper's
    Lemma 2 product formula, with the path lengths read off the decoded
    graph, and (ii) a generic bucket-elimination proper-colouring count for
    q = 0..32 followed by exact Newton interpolation.  Both the lemma and
    the generic counter are themselves validated against literal brute-force
    enumeration on small graphs.
  * P^(j)(H,-1/2) for j=0..3 are computed from the interpolated integer
    coefficients AND, independently, by rational jet arithmetic on the
    unexpanded product formula.
  * L_H^{(3)}(-1/2) is computed by the paper's identity AND independently as
    6 * [t^3] log(f(-1/2+t)/f(-1/2)); its positivity (the refutation) is
    computed, never asserted.
  * Conjecture 3's hypotheses and its upper inequality are recomputed,
    including L_{K_32}^{(3)}(-1/2) by two routes and an exhaustive check that
    ALL 2448 spanning trees of H are chordal proper spanning subgraphs.
    The bare criterion "not (L3(H) < L3(K_32))" is the correct negation of the
    upper inequality but is DELIBERATELY NOT presented as evidence about H:
    since L_{K_n}^{(k)}(x) < 0 always, that criterion is met by many graphs
    that are not counterexamples to Conjecture 2 at all.  The program
    therefore also computes the corollary's actual route (L3(H) > 0 >
    L3(K_32)), which is sensitive to the exhibited object, and computes and
    reports the least order at which a noncomplete generalized theta graph
    already meets the bare criterion (order 4, C_4).
  * the graph6 codec itself is validated against the reference example of the
    format specification and against K_2,...,K_7.
  * an EXHAUSTIVE census of every generalized theta graph of order <= 34
    (87606 graphs) recomputes the sign of L^{(3)}(-1/2) for each.
  * negative controls: the whole load-bearing pipeline is re-run on C_32 and
    on a one-bit mutation of the graph6 string.
"""
import sys
import math
import itertools
import collections
from fractions import Fraction as Fr

# ==========================================================================
# INPUTS COPIED FROM THE PAPER
# ==========================================================================
G6_H = ("_PCGGD@_??_@?@??_?I?@_????G??G??C??@???G???g??@_???????_???G???@"
        "????C????G????I????C")
PAPER_N = 32
PAPER_M = 34
PAPER_DELTA = 4
PAPER_LENGTHS = (8, 8, 9, 9)
PAPER_TERMINALS = (0, 1)
K_WITNESS = 3
X_WITNESS = Fr(-1, 2)
PAPER_P = Fr(16528743341798025, 4294967296)
PAPER_DP = Fr(-38758016767164565, 536870912)
PAPER_DDP = Fr(5549131014505347, 4194304)
PAPER_DDDP = Fr(-1586794036505230809, 67108864)
PAPER_L3 = Fr(5255579329573016681520193060821596608,
              18016243181285229121994454910138997625)
NY_CONST = Fr(301, 100)          # Ning--Yang: (10) holds for x <= -3.01*Delta*k
CENSUS_MAX_ORDER = 34            # exhaustive census bound (see census section)

RESULTS = []


def check(name, ok, detail=""):
    """Record and print one check.  ok must be a genuine computed boolean."""
    ok = bool(ok)
    RESULTS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    sys.stdout.flush()
    return ok


def note(text):
    print("NOTE " + text)
    sys.stdout.flush()


# ==========================================================================
# 1. exact polynomial arithmetic (ascending coefficient lists)
# ==========================================================================
def ptrim(a):
    while len(a) > 1 and a[-1] == 0:
        a = a[:-1]
    return a


def padd(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i, x in enumerate(a):
        r[i] += x
    for i, x in enumerate(b):
        r[i] += x
    return ptrim(r)


def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return ptrim(r)


def ppow(a, k):
    r = [1]
    for _ in range(k):
        r = pmul(r, a)
    return r


def pdiv_by_q(a):
    """Exact division of a by the polynomial q.  Raises if not divisible."""
    if a[0] != 0:
        raise ValueError("not divisible by q: constant term %r" % (a[0],))
    return ptrim(a[1:]) if len(a) > 1 else [0]


def peval(a, x):
    r = Fr(0)
    for c in reversed(a):
        r = r * x + c
    return r


def pder(a):
    if len(a) <= 1:
        return [0]
    return [i * c for i, c in enumerate(a)][1:]


# ==========================================================================
# 2. graph6 codec (written from the format spec, not imported)
# ==========================================================================
def g6_decode(s):
    """Return (n, adj) with adj a list of sets.  Raises on malformed input."""
    if not s:
        raise ValueError("empty graph6 string")
    b = [ord(c) - 63 for c in s]
    for v in b:
        if v < 0 or v > 63:
            raise ValueError("graph6 byte out of range 63..126")
    n = b[0]
    if n == 63:
        raise ValueError("multi-byte order form not supported")
    need = n * (n - 1) // 2
    nbytes = (need + 5) // 6
    if len(b) - 1 != nbytes:
        raise ValueError("graph6 payload is %d bytes, expected %d"
                         % (len(b) - 1, nbytes))
    bits = []
    for v in b[1:]:
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    if any(bits[need:]):
        raise ValueError("graph6 padding bits are not zero")
    adj = [set() for _ in range(n)]
    i = 0
    for j in range(1, n):
        for k in range(j):
            if bits[i]:
                adj[k].add(j)
                adj[j].add(k)
            i += 1
    return n, adj


def g6_encode(n, adj):
    bits = []
    for j in range(1, n):
        for k in range(j):
            bits.append(1 if j in adj[k] else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for i in range(0, len(bits), 6):
        v = 0
        for bit in bits[i:i + 6]:
            v = v * 2 + bit
        out.append(chr(v + 63))
    return "".join(out)


# ==========================================================================
# 3. elementary graph utilities
# ==========================================================================
def edge_list(n, adj):
    return sorted((u, v) for u in range(n) for v in adj[u] if u < v)


def degrees(n, adj):
    return [len(adj[u]) for u in range(n)]


def is_symmetric_loopless(n, adj):
    for u in range(n):
        if u in adj[u]:
            return False
        for v in adj[u]:
            if not (0 <= v < n) or u not in adj[v]:
                return False
    return True


def components(n, adj, allowed=None):
    if allowed is None:
        allowed = set(range(n))
    seen = set()
    out = []
    for s in sorted(allowed):
        if s in seen:
            continue
        seen.add(s)
        stack = [s]
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w in allowed and w not in seen:
                    seen.add(w)
                    stack.append(w)
        out.append(sorted(comp))
    return out


def is_connected(n, adj):
    return n > 0 and len(components(n, adj)) == 1


def bipartite_2colouring(n, adj):
    """Return a proper 2-colouring, or None if an odd cycle exists."""
    col = [-1] * n
    for s in range(n):
        if col[s] != -1:
            continue
        col[s] = 0
        dq = collections.deque([s])
        while dq:
            u = dq.popleft()
            for w in adj[u]:
                if col[w] == -1:
                    col[w] = 1 - col[u]
                    dq.append(w)
                elif col[w] == col[u]:
                    return None
    return col


# ==========================================================================
# 4. generalized theta graphs: construction and structural recognition
# ==========================================================================
def build_theta(lengths):
    """Canonical Theta(s_1,...,s_p): terminals 0 and 1, then interiors.

    A length-1 path is the edge {0,1} itself.  Returns (n, adj)."""
    lengths = list(lengths)
    if len(lengths) < 2 or any(s < 1 for s in lengths):
        raise ValueError("need p>=2 paths of positive length")
    if sum(1 for s in lengths if s == 1) > 1:
        raise ValueError("two length-1 paths would give a multi-edge")
    n = 2 + sum(s - 1 for s in lengths)
    adj = [set() for _ in range(n)]

    def link(a, b):
        adj[a].add(b)
        adj[b].add(a)

    nxt = 2
    for s in lengths:
        prev = 0
        for _ in range(s - 1):
            link(prev, nxt)
            prev = nxt
            nxt += 1
        link(prev, 1)
    return n, adj


def theta_decompose(n, adj):
    """Recognise a generalized theta graph from scratch.

    Returns (terminals, paths) where each path is the vertex sequence from
    the first terminal to the second.  Returns None if the graph is not a
    generalized theta graph."""
    if not is_connected(n, adj):
        return None
    deg = degrees(n, adj)
    branch = [v for v in range(n) if deg[v] != 2]
    if len(branch) != 2:
        return None
    u, w = branch
    if deg[u] != deg[w] or deg[u] < 2:
        return None
    paths = []
    for start in sorted(adj[u]):
        if start == w:
            paths.append([u, w])
            continue
        seq = [u, start]
        prev, cur = u, start
        guard = 0
        while cur != w:
            nxts = [z for z in adj[cur] if z != prev]
            if len(nxts) != 1:
                return None
            prev, cur = cur, nxts[0]
            seq.append(cur)
            guard += 1
            if guard > n + 2:
                return None
        paths.append(seq)
    covered = set()
    for seq in paths:
        covered.update(seq)
    if covered != set(range(n)):
        return None
    if sum(len(s) - 2 for s in paths) != n - 2:
        return None          # interiors must be pairwise disjoint
    return (u, w), paths


def theta_isomorphism(n, adj, decomp, target_lengths):
    """Build an explicit bijection from the decoded graph onto the
    canonically constructed Theta(target_lengths) and return
    (mapping, canonical_n, canonical_adj, edge_sets_equal).

    target_lengths is the length multiset the PAPER claims, not the multiset
    just recovered from the decoded graph.  Building the canonical target from
    the recovered lengths instead would make the bijection fit by
    construction: it could then never fail, and in particular it would still
    report success for Theta(7,9,9,9).  A mismatch between the recovered and
    the claimed multiset is reported here as "no isomorphism"."""
    (u, w), paths = decomp
    order = sorted(range(len(paths)), key=lambda i: len(paths[i]))
    lengths = [len(paths[i]) - 1 for i in order]
    want = sorted(target_lengths)
    if lengths != want:
        cn, cadj = build_theta(want)
        return None, cn, cadj, False
    cn, cadj = build_theta(want)
    cdec = theta_decompose(cn, cadj)
    if cdec is None:
        return None, cn, cadj, False
    (cu, cw), cpaths = cdec
    corder = sorted(range(len(cpaths)), key=lambda i: len(cpaths[i]))
    phi = {u: cu, w: cw}
    for a, b in zip(order, corder):
        pa, pb = paths[a], cpaths[b]
        if len(pa) != len(pb):
            return None, cn, cadj, False
        for x, y in zip(pa[1:-1], pb[1:-1]):
            if x in phi:
                return None, cn, cadj, False
            phi[x] = y
    if sorted(phi.keys()) != list(range(n)) or sorted(phi.values()) != list(range(cn)):
        return None, cn, cadj, False
    mapped = sorted(tuple(sorted((phi[a], phi[b]))) for a, b in edge_list(n, adj))
    same = mapped == edge_list(cn, cadj)
    return phi, cn, cadj, same


# ==========================================================================
# 5. planarity, proved two independent ways
# ==========================================================================
def kuratowski_branch_bound(n, adj):
    """A subdivision of K_5 needs 5 vertices of degree >= 4 and a
    subdivision of K_{3,3} needs 6 vertices of degree >= 3.  Returns
    (n_deg_ge_3, n_deg_ge_4)."""
    deg = degrees(n, adj)
    return (sum(1 for d in deg if d >= 3), sum(1 for d in deg if d >= 4))


def trace_faces(n, adj, rot):
    """Count the faces of the embedding given by rotation system rot
    (rot[v] = cyclic list of neighbours of v).  Returns (n_faces, ok)."""
    for v in range(n):
        if sorted(rot[v]) != sorted(adj[v]) or len(set(rot[v])) != len(rot[v]):
            return 0, False
    pos = {v: {x: i for i, x in enumerate(rot[v])} for v in range(n)}
    darts = set((a, b) for a in range(n) for b in adj[a])
    seen = set()
    faces = 0
    for d in sorted(darts):
        if d in seen:
            continue
        faces += 1
        cur = d
        while cur not in seen:
            seen.add(cur)
            a, b = cur
            r = rot[b]
            cur = (b, r[(pos[b][a] - 1) % len(r)])
    return faces, len(seen) == len(darts)


def theta_rotation_system(n, adj, decomp):
    """The standard planar rotation system of a generalized theta graph:
    the p paths leave one terminal in some order and reach the other in the
    reverse order."""
    (u, w), paths = decomp
    rot = {}
    rot[u] = [seq[1] for seq in paths]
    rot[w] = [seq[-2] for seq in paths][::-1]
    for seq in paths:
        for i in range(1, len(seq) - 1):
            rot[seq[i]] = [seq[i - 1], seq[i + 1]]
    return rot


# ==========================================================================
# 6. route A to P(H,q): the paper's Lemma 2 product formula
# ==========================================================================
QM1 = [-1, 1]          # the polynomial q-1
QQ = [0, 1]            # the polynomial q


def A_poly(s):
    """A_s(q) = ((q-1)^s + (q-1)(-1)^s)/q, exact division."""
    num = padd(ppow(QM1, s), [((-1) ** s) * c for c in QM1])
    return pdiv_by_q(num)


def D_poly(s):
    """D_s(q) = ((q-1)^s - (-1)^s)/q, exact division."""
    num = padd(ppow(QM1, s), [-((-1) ** s)])
    return pdiv_by_q(num)


def theta_chromatic_lemma(lengths):
    """q * ( prod A_{s_i} + (q-1) prod D_{s_i} )."""
    pa, pd = [1], [1]
    for s in lengths:
        pa = pmul(pa, A_poly(s))
        pd = pmul(pd, D_poly(s))
    return pmul(QQ, padd(pa, pmul(QM1, pd)))


# ==========================================================================
# 7. literal brute force: enumerate every colouring (small graphs only)
# ==========================================================================
def brute_force_colourings(n, adj, q):
    """Count proper q-colourings by enumerating all q^n maps."""
    if q == 0:
        return 1 if n == 0 else 0
    total = 0
    nbr_before = [[w for w in adj[v] if w < v] for v in range(n)]
    col = [0] * n

    def rec(v):
        nonlocal total
        if v == n:
            total += 1
            return
        for c in range(q):
            ok = True
            for w in nbr_before[v]:
                if col[w] == c:
                    ok = False
                    break
            if ok:
                col[v] = c
                rec(v + 1)

    rec(0)
    return total


def brute_force_path_counts(s, q):
    """Count colourings of the interior of a length-s path with the two
    endpoint colours equal (A) and different (D)."""
    if s < 1:
        raise ValueError("s>=1")
    if q < 2:
        raise ValueError("the combinatorial reading of A_s,D_s needs q>=2")
    same = 0
    diff = 0
    for interior in itertools.product(range(q), repeat=s - 1):
        for c1, bucket in ((0, "A"), (1, "D")):
            seq = (0,) + interior + (c1,)
            if all(seq[i] != seq[i + 1] for i in range(len(seq) - 1)):
                if bucket == "A":
                    same += 1
                else:
                    diff += 1
    return same, diff


# ==========================================================================
# 8. route B, part 1: generic bucket elimination.  Counts proper
#    q-colourings of ANY graph; the cost is exponential in the induced width
#    only, and every vertex of H other than the two terminals has degree 2.
# ==========================================================================
def count_colourings_dp(n, adj, q):
    """Exact number of proper q-colourings, by variable elimination."""
    if n == 0:
        return 1
    if q == 0:
        return 0
    factors = []
    for u in range(n):
        for v in adj[u]:
            if u < v:
                tab = {}
                for a in range(q):
                    for b in range(q):
                        if a != b:
                            tab[(a, b)] = 1
                factors.append(((u, v), tab))
    alive = set(range(n))
    scalar = 1
    while alive:
        nbr = {v: set() for v in alive}
        for vs, _tab in factors:
            for a in vs:
                for b in vs:
                    if a != b:
                        nbr[a].add(b)
        v = min(alive, key=lambda z: (len(nbr[z]), z))
        involved = [(vs, t) for vs, t in factors if v in vs]
        rest = [(vs, t) for vs, t in factors if v not in vs]
        S = sorted(nbr[v])
        if len(S) > 4:
            raise RuntimeError("induced width too large for this verifier")
        newtab = {}
        for assign in itertools.product(range(q), repeat=len(S)):
            d = dict(zip(S, assign))
            tot = 0
            for cv in range(q):
                d[v] = cv
                pr = 1
                for vs, t in involved:
                    pr *= t.get(tuple(d[x] for x in vs), 0)
                    if pr == 0:
                        break
                tot += pr
            if tot:
                newtab[assign] = tot
        alive.discard(v)
        if S:
            factors = rest + [(tuple(S), newtab)]
        else:
            factors = rest
            scalar *= newtab.get((), 0)
    return scalar


# ==========================================================================
# 9. route B, part 2: exact Newton interpolation over the rationals
# ==========================================================================
def pmul_fr(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return r


def interpolate(xs, ys):
    """Newton form; returns the ascending coefficient list over Fraction."""
    if len(xs) != len(set(xs)):
        raise ValueError("interpolation nodes must be distinct")
    cs = [Fr(y) for y in ys]
    for k in range(1, len(xs)):
        for i in range(len(xs) - 1, k - 1, -1):
            cs[i] = (cs[i] - cs[i - 1]) / Fr(xs[i] - xs[i - k])
    poly = [Fr(0)]
    basis = [Fr(1)]
    for k in range(len(xs)):
        if len(poly) < len(basis):
            poly = poly + [Fr(0)] * (len(basis) - len(poly))
        for i, b in enumerate(basis):
            poly[i] += cs[k] * b
        basis = pmul_fr(basis, [Fr(-xs[k]), Fr(1)])
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def chromatic_polynomial_generic(n, adj):
    """P(G,q) for an arbitrary graph, from n+1 colouring counts."""
    xs = list(range(n + 1))
    ys = [count_colourings_dp(n, adj, q) for q in xs]
    coeffs = interpolate(xs, ys)
    if any(c.denominator != 1 for c in coeffs):
        raise ArithmeticError("interpolated chromatic polynomial is not integral")
    return [int(c) for c in coeffs], ys


# ==========================================================================
# 10. rational jet (truncated Taylor) arithmetic at a point
# ==========================================================================
class Jet(object):
    """Taylor coefficients c_i of f(x0 + t) = sum_i c_i t^i, i < K."""
    __slots__ = ("c", "K")

    def __init__(self, c, K):
        self.c = list(c) + [Fr(0)] * (K - len(c))
        self.K = K

    def __add__(self, o):
        return Jet([a + b for a, b in zip(self.c, o.c)], self.K)

    def __mul__(self, o):
        r = [Fr(0)] * self.K
        for i, a in enumerate(self.c):
            if a:
                for j in range(self.K - i):
                    if o.c[j]:
                        r[i + j] += a * o.c[j]
        return Jet(r, self.K)

    def scaled(self, k):
        return Jet([k * a for a in self.c], self.K)

    def inv(self):
        if self.c[0] == 0:
            raise ZeroDivisionError("jet is not invertible at the base point")
        r = [Fr(0)] * self.K
        r[0] = Fr(1) / self.c[0]
        for i in range(1, self.K):
            s = Fr(0)
            for j in range(1, i + 1):
                s += self.c[j] * r[i - j]
            r[i] = -s / self.c[0]
        return Jet(r, self.K)

    def power(self, k):
        r = Jet([Fr(1)], self.K)
        for _ in range(k):
            r = r * self
        return r


def theta_jet(lengths, x0, K):
    """The Lemma-2 expression evaluated in jet arithmetic at x0, i.e. the
    first K Taylor coefficients of P(Theta(lengths), x0+t) computed WITHOUT
    ever expanding the polynomial in the monomial basis."""
    jq = Jet([x0, Fr(1)], K)
    one = Jet([Fr(1)], K)
    jqm1 = jq + one.scaled(Fr(-1))
    jqinv = jq.inv()
    pa, pd = one, one
    for s in lengths:
        sg = Fr((-1) ** s)
        pa = pa * ((jqm1.power(s) + jqm1.scaled(sg)) * jqinv)
        pd = pd * ((jqm1.power(s) + one.scaled(-sg)) * jqinv)
    return jq * (pa + jqm1 * pd)


def jet_from_poly(coeffs, x0, K):
    """Taylor coefficients of a polynomial at x0, by repeated differentiation
    of the monomial expansion (independent of theta_jet)."""
    out = []
    cur = [Fr(c) for c in coeffs]
    for i in range(K):
        out.append(peval(cur, x0) / Fr(math.factorial(i)))
        cur = pder(cur)
        if not cur:
            cur = [Fr(0)]
    return Jet(out, K)


# ==========================================================================
# 11. logarithmic derivatives
# ==========================================================================
def log_derivatives_from_jet(jet, sign):
    """Given Taylor coefficients of f at x0 and sign = (-1)^n, return
    [*, L'(x0), ..., L^(K-1)(x0)] for L = log(sign*f), computed as the
    formal logarithm of the power series -- no derivative identity used.

    Entry 0 is NOT L(x0): the series is expanded about log(sign*f(x0)) and the
    additive constant is dropped, so entry 0 is always 0.  Only entries with
    index >= 1 are meaningful, and only those are ever read."""
    K = jet.K
    f0 = sign * jet.c[0]
    if f0 <= 0:
        raise ValueError("sign*f(x0) must be positive for log to be defined")
    v = [Fr(0)] + [sign * c / f0 for c in jet.c[1:]]
    res = [Fr(0)] * K
    vp = Jet([Fr(1)], K)
    vj = Jet(v, K)
    for k in range(1, K):
        vp = vp * vj
        s = Fr((-1) ** (k - 1), k)
        for i in range(K):
            res[i] += s * vp.c[i]
    return [res[i] * Fr(math.factorial(i)) for i in range(K)]


def L3_by_paper_identity(f, f1, f2, f3):
    """(log f)''' = f'''/f - 3 f' f''/f^2 + 2 (f'/f)^3, exactly as displayed
    in the paper's proof."""
    return f3 / f - 3 * f1 * f2 / f ** 2 + 2 * (f1 / f) ** 3


def L2_by_identity(f, f1, f2):
    return f2 / f - (f1 / f) ** 2


# ==========================================================================
# 12. chordality (LexBFS + verification of the perfect elimination ordering)
# ==========================================================================
def is_chordal(n, adj):
    label = {v: () for v in range(n)}
    un = set(range(n))
    order = []
    for i in range(n - 1, -1, -1):
        v = max(un, key=lambda z: (label[z], -z))
        un.discard(v)
        order.append(v)
        for w in adj[v]:
            if w in un:
                label[w] = label[w] + (i,)
    order.reverse()
    pos = {v: i for i, v in enumerate(order)}
    for v in order:
        later = [w for w in adj[v] if pos[w] > pos[v]]
        if not later:
            continue
        p = min(later, key=lambda z: pos[z])
        for w in later:
            if w != p and w not in adj[p]:
                return False
    return True


def bfs_spanning_tree(n, adj):
    parent = {0: None}
    dq = collections.deque([0])
    tedges = []
    while dq:
        v = dq.popleft()
        for w in sorted(adj[v]):
            if w not in parent:
                parent[w] = v
                tedges.append((min(v, w), max(v, w)))
                dq.append(w)
    tadj = [set() for _ in range(n)]
    for a, b in tedges:
        tadj[a].add(b)
        tadj[b].add(a)
    return tedges, tadj, len(parent) == n


def is_forest(n, adj):
    """A graph is a forest iff every component C has |E(C)| = |C| - 1."""
    for comp in components(n, adj):
        cs = set(comp)
        ec = sum(1 for v in comp for w in adj[v] if w in cs and v < w)
        if ec != len(comp) - 1:
            return False
    return True


def complete_graph(n):
    return n, [set(range(n)) - {v} for v in range(n)]


def cycle_graph(n):
    adj = [set() for _ in range(n)]
    for v in range(n):
        adj[v].add((v + 1) % n)
        adj[(v + 1) % n].add(v)
    return n, adj


# ==========================================================================
# 13. K_n and its closed-form logarithmic derivative
# ==========================================================================
def complete_chromatic(n):
    p = [1]
    for j in range(n):
        p = pmul(p, [-j, 1])
    return p


def Kn_log_derivative_closed(n, k, x):
    """L_{K_n}^{(k)}(x) = -(k-1)! sum_{j=0}^{n-1} (j-x)^{-k}."""
    return -Fr(math.factorial(k - 1)) * sum(Fr(1) / (Fr(j) - x) ** k
                                            for j in range(n))


# ==========================================================================
# 14. exhaustive census of generalized theta graphs of bounded order
# ==========================================================================
def enumerate_theta_parameter_sets(max_order):
    """Every multiset {s_1,...,s_p}, p>=2, s_i>=1, at most one s_i=1
    (otherwise the graph has a multi-edge), with
    order = 2 + sum(s_i - 1) <= max_order.  Yields sorted tuples, no
    repetitions."""
    budget = max_order - 2
    seen = set()

    def rec(maxs, rem, cur):
        if len(cur) >= 1:
            for extra in ((), (1,)):
                ss = tuple(sorted(cur + list(extra)))
                if len(ss) < 2:
                    continue
                if sum(s - 1 for s in ss) > budget:
                    continue
                if ss not in seen:
                    seen.add(ss)
                    yield ss
        top = min(maxs, rem + 1)
        for s in range(top, 1, -1):
            if s - 1 <= rem:
                for out in rec(s, rem - (s - 1), cur + [s]):
                    yield out

    for out in rec(budget + 1, budget, []):
        yield out


def census(max_order, x, k):
    """For every generalized theta graph of order <= max_order, compute the
    sign of L^{(k)}(x) exactly.  Returns (count, violators, nonpositive)."""
    count = 0
    violators = []
    nonpositive = []
    for ss in enumerate_theta_parameter_sets(max_order):
        order = 2 + sum(s - 1 for s in ss)
        if order < 3:
            continue
        count += 1
        P = theta_chromatic_lemma(list(ss))
        jet = jet_from_poly(P, x, k + 1)
        if (-1) ** order * jet.c[0] <= 0:
            nonpositive.append(ss)
            continue
        Ls = log_derivatives_from_jet(jet, (-1) ** order)
        if Ls[k] > 0:
            violators.append((order, ss, Ls[k]))
    violators.sort(key=lambda t: (t[0], t[1]))
    return count, violators, nonpositive


# ==========================================================================
# STEP 1: decode the exhibited object and print it back
# ==========================================================================
def step_decode(ctx):
    section("Step 1: the exhibited object H, decoded from its graph6 string")
    try:
        n, adj = g6_decode(G6_H)
        ok = True
        detail = ""
    except Exception as exc:                       # noqa: BLE001
        n, adj = 0, []
        ok = False
        detail = repr(exc)
    check("g6_string_is_wellformed", ok, detail)
    if not ok:
        raise SystemExit(verdict())
    ctx["n"], ctx["adj"] = n, adj
    deg = degrees(n, adj)
    E = edge_list(n, adj)
    profile = sorted(collections.Counter(deg).items())
    note("graph6 = " + G6_H)
    note("order n = %d, size m = %d, degree profile (deg:count) = %s"
         % (n, len(E), profile))
    note("edge list = " + " ".join("%d-%d" % e for e in E))
    check("g6_reencodes_to_the_same_string", g6_encode(n, adj) == G6_H,
          g6_encode(n, adj))
    # validate the decoder against the reference example in the graph6 format
    # specification (n=5 with edges 0-2, 0-4, 1-3, 3-4 encodes as "DQc")
    rn, radj = g6_decode("DQc")
    ok_dec = (rn == 5 and edge_list(rn, radj) == [(0, 2), (0, 4), (1, 3), (3, 4)])
    for t in range(2, 8):
        kn, kadj = complete_graph(t)
        bitsneeded = t * (t - 1) // 2
        allones = chr(t + 63) + "".join(
            chr(63 + sum(1 << (5 - i) for i in range(6)
                         if 6 * blk + i < bitsneeded))
            for blk in range((bitsneeded + 5) // 6))
        if g6_encode(kn, kadj) != allones:
            ok_dec = False
        dn, dadj = g6_decode(allones)
        if dn != t or dadj != kadj:
            ok_dec = False
    check("g6_codec_reproduces_the_specification_example_DQc_and_all_K_n",
          ok_dec and is_symmetric_loopless(n, adj),
          "reference edges %s" % (edge_list(rn, radj),))
    check("order_is_%d_as_claimed" % PAPER_N, n == PAPER_N, "n=%d" % n)
    check("size_is_%d" % PAPER_M, len(E) == PAPER_M, "m=%d" % len(E))
    check("max_degree_is_%d_as_claimed" % PAPER_DELTA,
          max(deg) == PAPER_DELTA, "Delta=%d" % max(deg))
    check("degree_profile_is_two_of_degree_4_and_thirty_of_degree_2",
          profile == [(2, 30), (4, 2)], str(profile))
    check("H_is_connected", is_connected(n, adj))
    ctx["deg"] = deg
    ctx["E"] = E


# ==========================================================================
# STEP 2: H really is the generalized theta graph Theta(8,8,9,9)
# ==========================================================================
def step_identify(ctx):
    section("Step 2: H is the generalized theta graph Theta(8,8,9,9)")
    n, adj = ctx["n"], ctx["adj"]
    decomp = theta_decompose(n, adj)
    check("H_is_recognised_as_a_generalized_theta_graph", decomp is not None)
    ctx["decomp"] = decomp
    if decomp is None:
        for nm in ("terminals_are_vertices_0_and_1_as_the_paper_states",
                   "path_count_is_4", "path_edge_lengths_are_8_8_9_9",
                   "terminals_are_nonadjacent_so_all_lengths_exceed_one",
                   "order_equals_2_plus_sum_of_lengths_minus_one",
                   "size_equals_sum_of_path_lengths",
                   "explicit_bijection_onto_independently_built_"
                   "Theta_8_8_9_9_is_an_isomorphism"):
            check(nm, False, "H is not a generalized theta graph")
        check("H_is_not_bipartite_because_8+9_is_an_odd_cycle_length",
              bipartite_2colouring(n, adj) is None)
        check("H_is_noncomplete", len(ctx["E"]) < n * (n - 1) // 2)
        return
    (u, w), paths = decomp
    lengths = tuple(sorted(len(p) - 1 for p in paths))
    note("terminals found = (%d,%d); internally disjoint paths of edge-lengths %s"
         % (u, w, str(lengths)))
    for p in paths:
        note("  path %s" % ("-".join(str(v) for v in p),))
    check("terminals_are_vertices_0_and_1_as_the_paper_states",
          tuple(sorted((u, w))) == PAPER_TERMINALS, "(%d,%d)" % (u, w))
    check("path_count_is_4", len(paths) == 4, "p=%d" % len(paths))
    check("path_edge_lengths_are_8_8_9_9", lengths == PAPER_LENGTHS,
          str(lengths))
    check("terminals_are_nonadjacent_so_all_lengths_exceed_one",
          w not in adj[u] and min(lengths) >= 2)
    check("order_equals_2_plus_sum_of_lengths_minus_one",
          n == 2 + sum(s - 1 for s in lengths),
          "2+%s=%d" % ("+".join(str(s - 1) for s in lengths),
                       2 + sum(s - 1 for s in lengths)))
    check("size_equals_sum_of_path_lengths", len(ctx["E"]) == sum(lengths),
          "%d" % sum(lengths))
    phi, cn, cadj, same = theta_isomorphism(n, adj, decomp, PAPER_LENGTHS)
    check("explicit_bijection_onto_independently_built_Theta_8_8_9_9_is_an_isomorphism",
          same and phi is not None,
          ("phi=" + ",".join("%d->%d" % (k, phi[k]) for k in sorted(phi))
           if phi else "no bijection"))
    ctx["decomp"] = decomp
    ctx["lengths"] = list(lengths)
    ctx["canon"] = (cn, cadj)
    cols = bipartite_2colouring(n, adj)
    check("H_is_not_bipartite_because_8+9_is_an_odd_cycle_length", cols is None)
    K = complete_graph(n)[1]
    check("H_is_noncomplete", len(ctx["E"]) < n * (n - 1) // 2 and adj != K,
          "m=%d < C(32,2)=%d" % (len(ctx["E"]), n * (n - 1) // 2))


# ==========================================================================
# STEP 3: planarity (two independent proofs)
# ==========================================================================
def step_planar(ctx):
    section("Step 3: H is planar (Kuratowski branch count and Euler formula)")
    n, adj, E = ctx["n"], ctx["adj"], ctx["E"]
    b3, b4 = kuratowski_branch_bound(n, adj)
    note("vertices of degree>=3: %d ; of degree>=4: %d" % (b3, b4))
    check("no_K5_subdivision_since_fewer_than_5_vertices_have_degree_at_least_4",
          b4 < 5, "%d < 5" % b4)
    check("no_K33_subdivision_since_fewer_than_6_vertices_have_degree_at_least_3",
          b3 < 6, "%d < 6" % b3)
    check("planar_edge_bound_m_le_3n_minus_6",
          len(E) <= 3 * n - 6, "%d <= %d" % (len(E), 3 * n - 6))
    if ctx.get("decomp") is None:
        check("face_trace_visits_every_dart_exactly_once", False,
              "no theta decomposition, no rotation system to trace")
        check("rotation_system_satisfies_Euler_n_minus_m_plus_f_equals_2",
              False, "no theta decomposition")
    else:
        rot = theta_rotation_system(n, adj, ctx["decomp"])
        faces, traced = trace_faces(n, adj, rot)
        note("explicit rotation system traced %d faces over %d darts"
             % (faces, 2 * len(E)))
        check("face_trace_visits_every_dart_exactly_once", traced)
        check("rotation_system_satisfies_Euler_n_minus_m_plus_f_equals_2",
              traced and n - len(E) + faces == 2,
              "%d-%d+%d=%d" % (n, len(E), faces, n - len(E) + faces))
    # the face trace must be able to report a non-planar answer: K_5 with the
    # ascending rotation system is exhibited as a control.
    k5n, k5adj = complete_graph(5)
    rot5 = {v: sorted(k5adj[v]) for v in range(5)}
    f5, t5 = trace_faces(k5n, k5adj, rot5)
    b3_5, b4_5 = kuratowski_branch_bound(k5n, k5adj)
    check("control_K5_fails_the_branch_bound_and_Euler_equality",
          b4_5 >= 5 and not (t5 and k5n - 10 + f5 == 2),
          "K5: deg>=4 count %d, Euler value %d" % (b4_5, k5n - 10 + f5))


# ==========================================================================
# STEP 4: validate both chromatic-polynomial engines against brute force
# ==========================================================================
SMALL_THETAS = [(1, 2), (2, 2), (1, 2, 2), (2, 3), (2, 2, 2), (1, 3, 3),
                (3, 3), (2, 2, 3), (2, 2, 2, 2), (1, 2, 2, 2)]


def step_validate_engines(ctx):
    section("Step 4: both chromatic-polynomial engines vs literal brute force")
    ok_ad = True
    bad = None
    for s in range(1, 7):
        for q in range(2, 6):
            a_bf, d_bf = brute_force_path_counts(s, q)
            if peval(A_poly(s), Fr(q)) != a_bf or peval(D_poly(s), Fr(q)) != d_bf:
                ok_ad = False
                bad = (s, q)
    check("A_s_and_D_s_match_brute_force_path_colour_counts", ok_ad,
          "first mismatch %s" % (bad,) if bad else "s=1..6, q=2..5")
    ok_int = True
    for s in range(1, 13):
        numA = padd(ppow(QM1, s), [((-1) ** s) * c for c in QM1])
        numD = padd(ppow(QM1, s), [-((-1) ** s)])
        if peval(numA, Fr(0)) != 0 or peval(numD, Fr(0)) != 0:
            ok_int = False
        if pmul(QQ, A_poly(s)) != numA or pmul(QQ, D_poly(s)) != numD:
            ok_int = False
    check("numerators_of_A_s_D_s_vanish_at_q_0_so_both_lie_in_Z_q", ok_int,
          "s=1..12, q*A_s and q*D_s reconstruct the numerators")
    ok_lem = True
    bad = None
    for ss in SMALL_THETAS:
        tn, tadj = build_theta(ss)
        Pl = theta_chromatic_lemma(list(ss))
        for q in range(0, tn + 1):
            if peval(Pl, Fr(q)) != brute_force_colourings(tn, tadj, q):
                ok_lem = False
                bad = (ss, q)
    check("lemma_2_product_formula_matches_brute_force_on_%d_small_thetas"
          % len(SMALL_THETAS), ok_lem, "first mismatch %s" % (bad,)
          if bad else "orders 3..6")
    ok_dp = True
    bad = None
    for ss in SMALL_THETAS:
        tn, tadj = build_theta(ss)
        for q in range(0, tn + 1):
            if count_colourings_dp(tn, tadj, q) != brute_force_colourings(tn, tadj, q):
                ok_dp = False
                bad = (ss, q)
    for extra in [complete_graph(4), cycle_graph(5), build_theta((2, 2, 2, 2, 2))]:
        tn, tadj = extra
        for q in range(0, tn + 1):
            if count_colourings_dp(tn, tadj, q) != brute_force_colourings(tn, tadj, q):
                ok_dp = False
                bad = ("extra", tn, q)
    check("generic_bucket_elimination_matches_brute_force_on_small_graphs",
          ok_dp, "first mismatch %s" % (bad,) if bad else "13 graphs, all q<=n")


# ==========================================================================
# STEP 5: P(H,q) by two independent algorithms
# ==========================================================================
def step_chromatic(ctx):
    section("Step 5: P(H,q) computed twice, independently")
    n, adj = ctx["n"], ctx["adj"]
    lengths = ctx.get("lengths")
    P_gen, values = chromatic_polynomial_generic(n, adj)
    ctx["P"] = P_gen
    ctx["values"] = values
    note("P(H,q) = " + " ".join("%+d*q^%d" % (c, i)
                               for i, c in enumerate(P_gen) if c))
    note("colouring counts P(H,q) for q=0..5: %s" % (values[:6],))
    if lengths is None:
        check("generic_colour_count_plus_interpolation_equals_lemma_2_formula",
              False, "no theta decomposition, lemma route unavailable")
    else:
        P_lemma = theta_chromatic_lemma(lengths)
        check("generic_colour_count_plus_interpolation_equals_lemma_2_formula",
              P_lemma == P_gen,
              "degree %d, %d coefficients agree"
              % (len(P_gen) - 1, len(P_gen)))
    check("P_H_has_degree_32_and_is_monic",
          len(P_gen) - 1 == PAPER_N and P_gen[-1] == 1,
          "deg=%d lead=%d" % (len(P_gen) - 1, P_gen[-1]))
    check("P_H_second_coefficient_is_minus_m",
          P_gen[-2] == -len(ctx["E"]), "a_31=%d, m=%d" % (P_gen[-2], len(ctx["E"])))
    # Integrality must be tested on the RATIONAL interpolant, before the int()
    # cast inside chromatic_polynomial_generic: "isinstance(c, int)" on the
    # already-cast list is true by construction and can never fail.
    raw = interpolate(list(range(n + 1)), values)
    deg_here = len(P_gen) - 1
    check("P_H_coefficients_are_integers_and_alternate_in_sign",
          all(c.denominator == 1 for c in raw)
          and [int(c) for c in raw] == P_gen
          and all((-1) ** (deg_here - i) * c >= 0 for i, c in enumerate(P_gen))
          and all(c != 0 for c in P_gen[1:]),
          "rational interpolant has %d coefficients, all with denominator 1"
          % len(raw))
    note("because the coefficients alternate and none of a_1..a_32 vanishes, "
         "(-1)^n P(H,x) = sum_i |a_i| |x|^i > 0 for every x<0, so L_H is "
         "defined on all of (-infinity,0)")
    ok_pos = all((-1) ** PAPER_N * peval(P_gen, Fr(-j, 8)) > 0
                 for j in range(1, 401))
    check("minus_one_to_the_n_times_P_is_positive_at_400_negative_rationals",
          ok_pos, "x = -1/8, -2/8, ..., -400/8 (400 points)")
    check("P_H_vanishes_at_q_0_1_2_and_is_positive_at_q_3",
          values[0] == 0 and values[1] == 0 and values[2] == 0 and values[3] > 0,
          "P(H,3)=%d" % values[3])
    check("chromatic_number_of_H_is_3", values[2] == 0 and values[3] > 0)
    # third, structurally different route on a few points: transfer matrices
    if lengths is None:
        check("transfer_matrix_count_agrees_at_q_3_4_5_6", False,
              "no theta decomposition")
    else:
        ok_tm = True
        for q in (3, 4, 5, 6):
            if transfer_matrix_count(lengths, q) != values[q]:
                ok_tm = False
        check("transfer_matrix_count_agrees_at_q_3_4_5_6", ok_tm)


def transfer_matrix_count(lengths, q):
    """Count proper q-colourings of Theta(lengths) with (J-I)^s matrices: a
    third algorithm, independent of both the lemma and the elimination DP."""
    M = [[0 if a == b else 1 for b in range(q)] for a in range(q)]

    def mm(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(q)) for j in range(q)]
                for i in range(q)]

    def mp(X, k):
        R = [[1 if i == j else 0 for j in range(q)] for i in range(q)]
        while k:
            if k & 1:
                R = mm(R, X)
            X = mm(X, X)
            k >>= 1
        return R

    mats = [mp(M, s) for s in lengths]
    tot = 0
    for c in range(q):
        for d in range(q):
            pr = 1
            for A in mats:
                pr *= A[c][d]
            tot += pr
    return tot


# ==========================================================================
# STEP 6: the four derivative values displayed in the paper
# ==========================================================================
def step_derivatives(ctx):
    section("Step 6: P^(j)(H,-1/2) for j=0,1,2,3")
    P = ctx.get("P")
    if P is None:
        check("derivatives_available", False, "no chromatic polynomial")
        return
    x = X_WITNESS
    d0 = P
    d1 = pder(d0)
    d2 = pder(d1)
    d3 = pder(d2)
    got = [peval(d0, x), peval(d1, x), peval(d2, x), peval(d3, x)]
    want = [PAPER_P, PAPER_DP, PAPER_DDP, PAPER_DDDP]
    names = ["P", "P_prime", "P_double_prime", "P_triple_prime"]
    for i in range(4):
        note("%s(-1/2) = %s" % (names[i], got[i]))
        check("paper_value_of_%s_at_minus_one_half" % names[i],
              got[i] == want[i],
              "computed %s vs paper %s" % (got[i], want[i]))
    jet_coeff = jet_from_poly(P, x, 4)
    if ctx.get("lengths") is None:
        check("jet_arithmetic_on_the_unexpanded_product_formula_reproduces_all_four",
              False, "no theta decomposition")
    else:
        jet_formula = theta_jet(ctx["lengths"], x, 4)
        check("jet_arithmetic_on_the_unexpanded_product_formula_reproduces_all_four",
              jet_formula.c == jet_coeff.c
              and [jet_formula.c[i] * math.factorial(i)
                   for i in range(4)] == got)
    ctx["fjet"] = jet_coeff
    ctx["derivs"] = got


# ==========================================================================
# STEP 7: the refutation of Conjecture 2 -- the load-bearing computation
# ==========================================================================
def step_refutation(ctx):
    section("Step 7: L_H^(3)(-1/2) and the refutation of Conjecture 2")
    P = ctx.get("P")
    if P is None:
        check("refutation_computable", False, "no chromatic polynomial")
        return
    n, x = ctx["n"], X_WITNESS
    f, f1, f2, f3 = ctx["derivs"]
    sgn = (-1) ** n
    check("n_is_even_so_L_H_equals_log_P_H", sgn == 1 and n % 2 == 0,
          "(-1)^%d = %d" % (n, sgn))
    check("minus_one_to_the_n_times_P_at_the_witness_point_is_positive",
          sgn * f > 0, "value %s" % (sgn * f,))
    l3_identity = L3_by_paper_identity(sgn * f, sgn * f1, sgn * f2, sgn * f3)
    l_series = log_derivatives_from_jet(ctx["fjet"], sgn)
    note("L_H^(3)(-1/2) = %s" % l3_identity)
    note("            = %s (decimal, display only)"
         % repr(float(l3_identity)))
    check("L3_from_the_paper_identity_equals_the_paper_value",
          l3_identity == PAPER_L3,
          "computed %s" % l3_identity)
    check("L3_from_formal_log_of_the_power_series_agrees_independently",
          l_series[3] == l3_identity)
    check("witness_satisfies_the_hypotheses_k_at_least_2_and_x_negative",
          K_WITNESS >= 2 and x < 0, "k=%d, x=%s" % (K_WITNESS, x))
    check("REFUTATION_L3_is_strictly_positive_contradicting_conjecture_2",
          l3_identity > 0,
          "sign of L_H^(3)(-1/2) is %s; value %s"
          % ("positive" if l3_identity > 0 else
             ("zero" if l3_identity == 0 else "negative"), l3_identity))
    l2 = L2_by_identity(sgn * f, sgn * f1, sgn * f2)
    check("L2_agrees_with_the_series_route", l2 == l_series[2])
    check("L2_at_the_witness_point_is_negative_so_the_failure_is_third_order",
          l2 < 0, "L_H''(-1/2) = %s" % l2)
    ctx["L3H"] = l3_identity


# ==========================================================================
# STEP 8: consistency with the proved range of Ning--Yang Theorem 1.5
# ==========================================================================
def L3_of_poly_at(P, n, x):
    jet = jet_from_poly(P, x, 4)
    if (-1) ** n * jet.c[0] <= 0:
        raise ValueError("(-1)^n P(x) is not positive at x=%s" % (x,))
    return log_derivatives_from_jet(jet, (-1) ** n)[3]


def step_ning_yang(ctx):
    section("Step 8: the witness point lies outside the proved range")
    P, n = ctx.get("P"), ctx["n"]
    if P is None:
        check("ning_yang_range_check_possible", False, "no polynomial")
        return
    delta = max(ctx["deg"])
    bound = -NY_CONST * Fr(delta) * Fr(K_WITNESS)
    note("Ning--Yang prove the inequality for x <= -3.01*Delta*k = %s" % bound)
    note("PROVENANCE of that bound: the constant 3.01 is TRANSCRIBED from "
         "Ning--Yang Theorem 1.5 and is NOT recomputed here (this program "
         "cannot see that paper); Delta=%d is read off the decoded graph and "
         "k=%d is the witness order, and the comparison below is exact given "
         "the transcribed constant" % (delta, K_WITNESS))
    check("witness_point_is_not_inside_the_ning_yang_proved_range",
          X_WITNESS > bound, "%s > %s" % (X_WITNESS, bound))
    deep = [bound, Fr(-37), Fr(-50), Fr(-100), Fr(-1000)]
    ok = True
    for x in deep:
        if L3_of_poly_at(P, n, x) >= 0:
            ok = False
    check("L3_of_H_is_negative_at_every_sampled_point_in_the_proved_range",
          ok, "sampled x in %s" % ([str(v) for v in deep],))
    coarse = [Fr(-j, 4) for j in range(1, 121)]
    pos = [x for x in coarse if L3_of_poly_at(P, n, x) > 0]
    note("on the grid x=-1/4,-2/4,...,-30 the set where L_H^(3)(x)>0 is "
         + (("%s ... %s (%d of %d points)"
             % (pos[0], pos[-1], len(pos), len(coarse))) if pos else "empty"))
    check("the_positivity_region_is_a_nonempty_proper_subset_of_the_coarse_grid",
          0 < len(pos) < len(coarse))
    check("the_witness_point_lies_in_the_computed_positivity_region",
          X_WITNESS in pos)
    fine = [Fr(-j, 1000) for j in range(300, 601)]
    signs = [L3_of_poly_at(P, n, x) > 0 for x in fine]
    idx = [i for i, s in enumerate(signs) if s]
    contiguous = bool(idx) and idx == list(range(idx[0], idx[-1] + 1))
    if idx:
        note("on the grid of step 1/1000 the positivity window is "
             "[%s, %s] (%d points), and the sign is negative on both sides"
             % (fine[idx[-1]], fine[idx[0]], len(idx)))
    check("positivity_window_is_a_single_interval_bracketed_by_negative_values",
          contiguous and idx[0] > 0 and idx[-1] < len(fine) - 1)


# ==========================================================================
# STEP 9: Corollary 3 -- the upper inequality of Conjecture 3 fails
# ==========================================================================
def all_theta_spanning_trees(n, decomp):
    """Every spanning tree of a generalized theta graph: leave one path whole
    and delete exactly one edge from each of the others."""
    (_u, _w), paths = decomp
    pedges = []
    for seq in paths:
        pedges.append([(min(seq[i], seq[i + 1]), max(seq[i], seq[i + 1]))
                       for i in range(len(seq) - 1)])
    allE = set(e for pe in pedges for e in pe)
    for keep in range(len(paths)):
        others = [i for i in range(len(paths)) if i != keep]
        for choice in itertools.product(*[pedges[i] for i in others]):
            yield allE - set(choice)


def step_conjecture3(ctx):
    section("Step 9: Conjecture 3, upper inequality, and its hypotheses")
    n, adj = ctx["n"], ctx["adj"]
    if "L3H" not in ctx:
        check("conjecture_3_check_possible", False, "no L3 for H")
        return
    x, k = X_WITNESS, K_WITNESS
    PK = complete_chromatic(n)
    l3k_poly = L3_of_poly_at(PK, n, x)
    l3k_closed = Kn_log_derivative_closed(n, k, x)
    note("L_{K_32}^(3)(-1/2) = %s" % l3k_closed)
    check("Kn_closed_form_matches_the_polynomial_route",
          l3k_poly == l3k_closed)
    check("Kn_chromatic_polynomial_matches_falling_factorial_values",
          all(peval(PK, Fr(q)) == 0 for q in range(n))
          and peval(PK, Fr(n)) == math.factorial(n))
    ok_prod = True
    for j in range(1, 40):
        xx = Fr(-j, 4)
        lhs = (-1) ** n * peval(PK, xx)
        rhs = Fr(1)
        for t in range(n):
            rhs *= (Fr(t) - xx)
        if lhs != rhs or lhs <= 0:
            ok_prod = False
    check("displayed_identity_minus_one_to_the_n_P_Kn_equals_prod_j_minus_x",
          ok_prod, "x = -1/4, ..., -39/4")
    check("L3_of_K32_is_negative_as_the_corollary_asserts", l3k_closed < 0)
    check("UPPER_INEQUALITY_OF_CONJECTURE_3_FAILS_for_H_k3_x_minus_half",
          not (ctx["L3H"] < l3k_closed),
          "L3(H)=%s is not < L3(K_32)=%s" % (float(ctx["L3H"]),
                                             float(l3k_closed)))
    # The criterion just applied -- "not (L3(G) < L3(K_n))" -- is the correct
    # negation of the upper inequality, but on its own it has almost no
    # discriminating power: L_{K_n}^{(k)}(x) < 0 always, so ANY graph whose
    # L^{(k)}(x) is merely greater than that bounded negative number fails it.
    # Two things are therefore computed rather than assumed.  First, the route
    # the paper's Corollary 3 actually takes: the failure is a consequence of
    # Theorem 1, i.e. L3(H) is positive while L3(K_32) is negative.  Unlike the
    # bare criterion, this DOES fail if the exhibited object is corrupted.
    check("upper_inequality_fails_via_the_corollarys_route_L3_H_positive_"
          "and_L3_K32_negative",
          ctx["L3H"] > 0 > l3k_closed,
          "L3(H)=%s > 0 > L3(K_32)=%s" % (float(ctx["L3H"]),
                                          float(l3k_closed)))
    # Second, the weakness of the bare criterion is measured, not hidden: the
    # least order at which a NONCOMPLETE generalized theta graph already fails
    # the upper inequality at k=3, x=-1/2 is computed by exhaustive search.
    small_fail = []
    for ss in enumerate_theta_parameter_sets(8):
        order = 2 + sum(s - 1 for s in ss)
        if order < 4 or sum(ss) == order * (order - 1) // 2:
            continue                      # skip complete graphs: Conj 3 needs
                                          # G noncomplete
        Ps = theta_chromatic_lemma(list(ss))
        if not (L3_of_poly_at(Ps, order, x)
                < Kn_log_derivative_closed(order, k, x)):
            small_fail.append((order, ss))
    small_fail.sort()
    least = small_fail[0][0] if small_fail else None
    note("SCOPE of Corollary 3: the bare criterion 'L^(3)(-1/2) >= "
         "L_{K_n}^(3)(-1/2)' is ALSO met by much smaller graphs -- the least "
         "order among noncomplete generalized theta graphs is %s, attained by "
         "%s (Theta(2,2) and Theta(1,3) are both C_4). So the upper "
         "inequality of Conjecture 3 is not specific to H; what is specific "
         "to H is that its L^(3) is POSITIVE, which is what refutes "
         "Conjecture 2." % (least, ", ".join("Theta%s" % (s,)
                                             for _o, s in small_fail
                                             if _o == least)))
    check("least_order_noncomplete_theta_graph_also_failing_the_upper_"
          "inequality_is_C4_of_order_4",
          least == 4 and (2, 2) in [s for _o, s in small_fail if _o == least],
          "least order = %s" % least)
    # chordality routine, validated on graphs whose status is known
    c4 = cycle_graph(4)
    k4 = complete_graph(4)
    check("chordality_routine_is_correct_on_C4_K4_C5_and_H",
          (not is_chordal(*c4)) and is_chordal(*k4)
          and (not is_chordal(*cycle_graph(5))) and (not is_chordal(n, adj)),
          "C4 no, K4 yes, C5 no, H no")
    tedges, tadj, spans = bfs_spanning_tree(n, adj)
    Eset = set(ctx["E"])
    check("Q_is_a_spanning_subgraph_of_H_with_n_minus_1_edges",
          spans and len(tedges) == n - 1 and set(tedges) <= Eset,
          "|E(Q)|=%d, spans=%s" % (len(tedges), spans))
    check("Q_is_a_tree_connected_and_acyclic",
          is_connected(n, tadj) and is_forest(n, tadj))
    check("Q_is_a_proper_subgraph_of_H", len(tedges) < len(ctx["E"]),
          "%d < %d" % (len(tedges), len(ctx["E"])))
    check("Q_is_chordal", is_chordal(n, tadj))
    PQ, _vals = chromatic_polynomial_generic(n, tadj)
    check("P_Q_equals_q_times_q_minus_1_to_the_n_minus_1",
          PQ == pmul(QQ, ppow(QM1, n - 1)))
    l3q = L3_of_poly_at(PQ, n, x)
    note("L_Q^(3)(-1/2) = %s for the BFS spanning tree Q" % l3q)
    check("lower_inequality_of_conjecture_3_still_holds_for_this_Q",
          l3q < ctx["L3H"], "%s < %s" % (float(l3q), float(ctx["L3H"])))
    check("H_is_connected_noncomplete_so_conjecture_3_applies",
          is_connected(n, adj) and len(ctx["E"]) < n * (n - 1) // 2)
    # the corollary claims ANY spanning tree works: enumerate all of them.
    if ctx.get("decomp") is None:
        check("EVERY_spanning_tree_of_H_is_a_chordal_proper_spanning_subgraph",
              False, "no theta decomposition")
        return
    total = 0
    ok_all = True
    seen = set()
    for tset in all_theta_spanning_trees(n, ctx["decomp"]):
        total += 1
        key = tuple(sorted(tset))
        if key in seen:
            ok_all = False
        seen.add(key)
        sadj = [set() for _ in range(n)]
        for a, b in tset:
            sadj[a].add(b)
            sadj[b].add(a)
        if not (len(tset) == n - 1 and tset <= Eset and len(tset) < len(Eset)
                and is_connected(n, sadj) and is_forest(n, sadj)
                and is_chordal(n, sadj)):
            ok_all = False
    paths = ctx["decomp"][1]
    expected = 0
    for i in range(len(paths)):
        prod = 1
        for j, p in enumerate(paths):
            if j != i:
                prod *= len(p) - 1
        expected += prod
    note("enumerated %d spanning trees of H (formula prod s_i * sum 1/s_i "
         "predicts %d)" % (total, expected))
    check("EVERY_spanning_tree_of_H_is_a_chordal_proper_spanning_subgraph",
          ok_all and total == expected == len(seen) and total > 0,
          "%d distinct spanning trees, all chordal proper spanning subgraphs"
          % total)


# ==========================================================================
# STEP 10: negative controls -- the pipeline must be able to say NO
# ==========================================================================
def flip_g6_bit(s, bitpos):
    b = [ord(c) - 63 for c in s]
    n = b[0]
    need = n * (n - 1) // 2
    if not 0 <= bitpos < need:
        raise ValueError("bit position outside the adjacency bits")
    byte = 1 + bitpos // 6
    off = bitpos % 6
    b[byte] ^= (1 << (5 - off))
    return "".join(chr(v + 63) for v in b)


def step_controls(ctx):
    section("Step 10: negative controls (non-vacuity of the pipeline)")
    n = ctx["n"]
    x = X_WITNESS
    cn, cadj = cycle_graph(n)
    PC, _v = chromatic_polynomial_generic(cn, cadj)
    check("control_pipeline_on_C32_reproduces_the_known_cycle_polynomial",
          PC == padd(ppow(QM1, n), QM1),
          "degree %d" % (len(PC) - 1))
    l3c = L3_of_poly_at(PC, cn, x)
    check("control_C32_satisfies_conjecture_2_at_the_same_point_L3_negative",
          l3c < 0, "L3(C_32)(-1/2) = %s" % l3c)
    ok_mut = True
    detail = []
    for bitpos in (0, 100, 400):
        g6m = flip_g6_bit(G6_H, bitpos)
        mn, madj = g6_decode(g6m)
        if madj == ctx["adj"]:
            ok_mut = False
            detail.append("bit %d did not change the graph" % bitpos)
            continue
        Pm, _vm = chromatic_polynomial_generic(mn, madj)
        l3m = L3_of_poly_at(Pm, mn, x)
        if l3m == PAPER_L3:
            ok_mut = False
            detail.append("bit %d still gave the paper value" % bitpos)
        else:
            detail.append("bit %d -> L3=%s" % (bitpos, repr(float(l3m))))
    check("control_one_bit_mutations_of_the_graph6_string_all_change_L3",
          ok_mut, "; ".join(detail))
    ok_near = True
    for ss in [(8, 8, 9, 10), (7, 9, 9, 9), (8, 8, 8, 9)]:
        order = 2 + sum(s - 1 for s in ss)
        if L3_of_poly_at(theta_chromatic_lemma(list(ss)), order, x) == PAPER_L3:
            ok_near = False
    check("control_neighbouring_theta_graphs_give_values_other_than_the_paper_value",
          ok_near, "Theta(8,8,9,10), Theta(7,9,9,9), Theta(8,8,8,9)")


# ==========================================================================
# STEP 11: exhaustive census of generalized theta graphs of small order
# ==========================================================================
def brute_parameter_sets(max_order):
    out = set()
    for p in range(2, max_order + 1):
        for ss in itertools.combinations_with_replacement(
                range(1, max_order + 1), p):
            if sum(1 for s in ss if s == 1) > 1:
                continue
            if 2 + sum(s - 1 for s in ss) > max_order:
                continue
            out.add(tuple(sorted(ss)))
    return out


def step_census(ctx):
    section("Step 11: exhaustive census of generalized theta graphs, order <= %d"
            % CENSUS_MAX_ORDER)
    ok_enum = True
    for mo in (8, 12, 14):
        if set(enumerate_theta_parameter_sets(mo)) != brute_parameter_sets(mo):
            ok_enum = False
    check("census_enumerator_is_complete_vs_brute_force_multiset_enumeration",
          ok_enum, "cross-checked at order bounds 8, 12, 14")
    ok_build = True
    tested = 0
    for ss in list(enumerate_theta_parameter_sets(10)):
        bn, badj = build_theta(ss)
        dec = theta_decompose(bn, badj)
        if len(ss) == 2:
            # a 2-path theta graph is a cycle: every degree is 2, so the two
            # terminals are not recoverable and the recogniser must decline.
            if dec is not None:
                ok_build = False
            continue
        tested += 1
        if dec is None or tuple(sorted(len(p) - 1 for p in dec[1])) != ss:
            ok_build = False
    check("every_enumerated_parameter_set_builds_a_graph_that_decomposes_back",
          ok_build, "%d graphs with p>=3 of order <= 10, plus cycles declined"
          % tested)
    count, viol, nonpos = census(CENSUS_MAX_ORDER, X_WITNESS, K_WITNESS)
    note("census examined %d generalized theta graphs (all p, all path "
         "lengths, order 3..%d)" % (count, CENSUS_MAX_ORDER))
    check("in_the_census_minus_one_to_the_n_times_P_is_positive_at_minus_half_"
          "for_every_member", not nonpos,
          "%d exceptions" % len(nonpos))
    check("census_found_at_least_one_violator_of_conjecture_2", len(viol) > 0,
          "%d violators" % len(viol))
    minorder = min(o for o, _s, _v in viol) if viol else None
    smallest = [s for o, s, _v in viol if o == minorder]
    note("the census finds %d graphs with L^(3)(-1/2)>0; the least order at "
         "which this happens is %s, attained exactly by %s"
         % (len(viol), minorder, ", ".join("Theta%s" % (s,) for s in smallest)))
    check("census_least_order_with_L3_positive_at_minus_one_half_is_32",
          minorder == PAPER_N,
          "least order with L^(3)(-1/2)>0 = %s; k=%d and x=%s only"
          % (minorder, K_WITNESS, X_WITNESS))
    check("the_paper_witness_is_one_of_the_least_order_violators_in_the_census",
          PAPER_LENGTHS in [tuple(s) for s in smallest],
          "least-order violators: %s" % (smallest,))
    note("SCOPE: the paper claims a single witness and that claim is "
         "re-derived above in full. The census in this step is an ADDITION, "
         "not a paper claim; it is exhaustive over generalized theta graphs "
         "of order <= %d only, and for each member it evaluates the sign of "
         "L^(k)(x) at the single pair k=%d, x=%s. NOT RE-RUN: orders above %d "
         "(the number of parameter multisets grows like the partition "
         "function); any graph outside the generalized theta family; and any "
         "other (k,x) with k>=2, x<0 -- so a census member of order below %d "
         "could still violate Conjecture 2 at some other (k,x) without this "
         "census seeing it, and the least-order statement above is a statement "
         "about the sign of L^(3)(-1/2), not about Conjecture 2 in full. ALSO "
         "NOT RE-RUN: the constant 3.01 of Ning--Yang Theorem 1.5 used in "
         "Step 8, which is transcribed from that paper and not recomputed."
         % (CENSUS_MAX_ORDER, K_WITNESS, X_WITNESS, CENSUS_MAX_ORDER,
            PAPER_N))


STEPS = [step_decode, step_identify, step_planar, step_validate_engines,
         step_chromatic, step_derivatives, step_refutation, step_ning_yang,
         step_conjecture3, step_controls, step_census]


def section(title):
    print("")
    print("=== " + title)
    sys.stdout.flush()


def main():
    ctx = {}
    for step in STEPS:
        try:
            step(ctx)
        except SystemExit:
            raise
        except Exception:                          # noqa: BLE001
            import traceback
            traceback.print_exc()
            check("step_%s_completed_without_an_exception" % step.__name__,
                  False, "see traceback above")
    return verdict()


def verdict():
    n = len(RESULTS)
    bad = [nm for nm, ok in RESULTS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    print("verification of the note: a third-derivative counterexample for")
    print("chromatic polynomials -- H = Theta(8,8,9,9), k=3, x=-1/2")
    print("python %s, exact arithmetic only" % sys.version.split()[0])
    try:
        code = main()
    except SystemExit:
        raise
    except Exception:                              # noqa: BLE001
        import traceback
        traceback.print_exc()
        check("verifier_ran_to_completion_without_an_exception", False)
        code = verdict()
    sys.exit(code)

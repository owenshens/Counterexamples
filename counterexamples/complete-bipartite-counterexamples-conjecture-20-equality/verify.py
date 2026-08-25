#!/usr/bin/env python3
"""Verification program for the note

    "Counterexamples to the Equality Clause in Conjecture 20 of
     Chen--Guo--Li--Wang".

Conjecture 20 (as quoted in the note): for a simple graph G on n vertices and
every integer k with 1 <= k <= 3n/4,

    mu_k(G) * mu_k(complement G) <= n(n-k),

with equality if and only if G or complement G is isomorphic to K_k v H, where
H is a disconnected graph on n-k vertices with at least k+1 components.

The decisive content verified here is the counterexample family: for t >= 2,
G = K_{2t,2t} with n = 4t, k = 3t attains equality, while the proposed
equality family is empty for those parameters.  All Laplacian spectra are
established with exact integer / rational arithmetic:

  * a claimed spectrum is proved by exact Gaussian elimination -- for each
    claimed eigenvalue lambda the nullity of L - lambda*I is computed over
    Fraction and matched against the claimed multiplicity, and the
    multiplicities are required to sum to n (which forces the multiset to be
    the whole spectrum);
  * independently, eigenvalue counts are obtained from the integer
    characteristic polynomial by an exact Taylor shift plus Descartes' rule of
    signs, which is exact (not merely an upper bound) because a Laplacian
    characteristic polynomial is real-rooted.

Standard library only.  No floating point is used anywhere.
"""

from fractions import Fraction
from itertools import combinations, permutations

# ----------------------------------------------------------------------
# verdict harness
# ----------------------------------------------------------------------

_RESULTS = []


def ck(name, ok, detail=""):
    ok = bool(ok)
    _RESULTS.append((name, ok))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         ("  --  " + detail) if detail else ""))
    return ok


def verdict():
    total = len(_RESULTS)
    failed = [nm for nm, ok in _RESULTS if not ok]
    if failed:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


# ----------------------------------------------------------------------
# graphs (adjacency as a list of sets of neighbours)
# ----------------------------------------------------------------------

def complete_bipartite(a, b):
    n = a + b
    adj = [set() for _ in range(n)]
    for u in range(a):
        for v in range(a, n):
            adj[u].add(v)
            adj[v].add(u)
    return adj


def disjoint_cliques(sizes):
    n = sum(sizes)
    adj = [set() for _ in range(n)]
    base = 0
    for s in sizes:
        for u, v in combinations(range(base, base + s), 2):
            adj[u].add(v)
            adj[v].add(u)
        base += s
    return adj


def complement(adj):
    n = len(adj)
    return [set(v for v in range(n) if v != u and v not in adj[u])
            for u in range(n)]


def edge_set(adj):
    """The edges of `adj` as one frozenset({u, v}) per edge.

    Both orientations (u, v) and (v, u) collapse to the same frozenset, so
    len(edge_set(adj)) is the number of edges -- it must NOT be halved.
    """
    return set(frozenset((u, v)) for u in range(len(adj)) for v in adj[u])


def degrees(adj):
    return [len(s) for s in adj]


def laplacian(adj):
    n = len(adj)
    L = [[0] * n for _ in range(n)]
    for u in range(n):
        L[u][u] = len(adj[u])
        for v in adj[u]:
            L[u][v] = -1
    return L


def components(adj, vertices=None):
    """Number of connected components of the subgraph induced on `vertices`."""
    if vertices is None:
        vertices = set(range(len(adj)))
    else:
        vertices = set(vertices)
    seen = set()
    count = 0
    for s in sorted(vertices):
        if s in seen:
            continue
        count += 1
        stack = [s]
        seen.add(s)
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w in vertices and w not in seen:
                    seen.add(w)
                    stack.append(w)
    return count


def isomorphic(adj1, adj2):
    """Brute-force isomorphism test (used only on 8-vertex graphs)."""
    n = len(adj1)
    if n != len(adj2):
        return False
    if sorted(degrees(adj1)) != sorted(degrees(adj2)):
        return False
    e2 = edge_set(adj2)
    if len(edge_set(adj1)) != len(e2):
        return False
    for p in permutations(range(n)):
        ok = True
        for u in range(n):
            for v in adj1[u]:
                if u < v and frozenset((p[u], p[v])) not in e2:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False


# ----------------------------------------------------------------------
# exact linear algebra
# ----------------------------------------------------------------------

def nullity(M):
    """Exact nullity of a square matrix over the rationals."""
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    rank = 0
    row = 0
    for col in range(n):
        piv = None
        for r in range(row, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        A[row], A[piv] = A[piv], A[row]
        pv = A[row][col]
        for r in range(row + 1, n):
            if A[r][col] != 0:
                f = A[r][col] / pv
                for c in range(col, n):
                    A[r][c] -= f * A[row][c]
        rank += 1
        row += 1
        if row == n:
            break
    return n - rank


def shifted_identity(M, lam):
    n = len(M)
    return [[Fraction(M[i][j]) - (Fraction(lam) if i == j else 0)
             for j in range(n)] for i in range(n)]


def spectrum_is(M, claimed):
    """`claimed` is a list of (eigenvalue, multiplicity) pairs.

    Returns (ok, detail).  Exact: each geometric multiplicity is computed by
    elimination, and their sum must be the matrix order, which leaves no room
    for any further eigenvalue and forces the algebraic multiplicities to
    agree.
    """
    n = len(M)
    got = []
    for lam, mult in claimed:
        d = nullity(shifted_identity(M, lam))
        got.append((lam, d, mult))
    if sum(d for _, d, _ in got) != n:
        return False, "multiplicities sum to %d, not %d" % (
            sum(d for _, d, _ in got), n)
    for lam, d, mult in got:
        if d != mult:
            return False, "eigenvalue %s has multiplicity %d, claimed %d" % (
                lam, d, mult)
    return True, "verified %s" % (
        ", ".join("%s^%d" % (lam, m) for lam, m in claimed),)


def charpoly(M):
    """Coefficients of det(x*I - M), descending powers, via Faddeev-LeVerrier.

    Exact over Fraction; for an integer matrix the result is integral.
    """
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    coeffs = [Fraction(1)]
    Mk = [[Fraction(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        # Mk <- A * Mk + coeffs[-1] * I
        NM = [[sum(A[i][t] * Mk[t][j] for t in range(n)) for j in range(n)]
              for i in range(n)]
        for i in range(n):
            NM[i][i] += coeffs[-1]
        Mk = NM
        tr = sum(sum(A[i][t] * Mk[t][i] for t in range(n)) for i in range(n))
        coeffs.append(-tr / k)
    return coeffs


def taylor_shift(coeffs, x):
    """Coefficients (descending) of p(y + x) for p given descending."""
    x = Fraction(x)
    cur = [Fraction(c) for c in coeffs]
    asc = []
    while cur:
        acc = Fraction(0)
        q = []
        for c in cur:
            acc = acc * x + c
            q.append(acc)
        asc.append(q.pop())
        cur = q
    return asc[::-1]


def sign_changes(coeffs):
    s = [c for c in coeffs if c != 0]
    return sum(1 for i in range(len(s) - 1) if (s[i] > 0) != (s[i + 1] > 0))


def count_gt(coeffs, x):
    """Number of eigenvalues strictly greater than x (with multiplicity).

    Exact for real-rooted p: Descartes' rule of signs is attained.
    """
    return sign_changes(taylor_shift(coeffs, x))


def mult_at(coeffs, x):
    q = taylor_shift(coeffs, x)
    m = 0
    for c in reversed(q):
        if c != 0:
            break
        m += 1
    return m


def mu_bracket(coeffs, n, k):
    """Locate the k-th largest root of a real-rooted p (all roots in [0, n]).

    Returns ('exact', Fraction) when the root is rational -- for a monic
    integer characteristic polynomial that means an integer -- and otherwise
    ('open', lo, hi) with lo < mu_k < hi and hi - lo = 1.
    """
    lo = Fraction(0)
    for c in range(0, n + 1):
        g = count_gt(coeffs, Fraction(c))
        if g < k <= g + mult_at(coeffs, Fraction(c)):
            return ("exact", Fraction(c))
        if g >= k:
            lo = Fraction(c)
    return ("open", lo, lo + 1)


def refine(coeffs, k, br):
    """One exact bisection step on an ('open', lo, hi) bracket."""
    if br[0] == "exact":
        return br
    lo, hi = br[1], br[2]
    mid = (lo + hi) / 2
    g = count_gt(coeffs, mid)
    if g < k <= g + mult_at(coeffs, mid):
        return ("exact", mid)
    if g >= k:
        return ("open", mid, hi)
    return ("open", lo, mid)


def bounds(br):
    if br[0] == "exact":
        return br[1], br[1]
    return br[1], br[2]


def mu_exact(coeffs, n, k):
    """The k-th largest root when it is rational, else None."""
    br = mu_bracket(coeffs, n, k)
    return br[1] if br[0] == "exact" else None


def product_vs(pG, pC, n, k, target, budget=60):
    """Compare mu_k(G)*mu_k(complement G) with `target`.

    Returns "lt", "eq", "gt", or "undecided".  Eigenvalues are non-negative,
    so the product of the interval endpoints brackets the product; the
    brackets are bisected exactly until the comparison is settled.
    """
    a = mu_bracket(pG, n, k)
    b = mu_bracket(pC, n, k)
    for _ in range(budget + 1):
        la, ha = bounds(a)
        lb, hb = bounds(b)
        if a[0] == "exact" and b[0] == "exact":
            p = la * lb
            return "lt" if p < target else ("eq" if p == target else "gt")
        if la * lb > target:
            return "gt"
        if ha * hb < target:
            return "lt"
        if a[0] == "open" and (b[0] == "exact" or (ha - la) >= (hb - lb)):
            a = refine(pG, k, a)
        else:
            b = refine(pC, k, b)
    return "undecided"


# ----------------------------------------------------------------------
# the conjecture's equality family
# ----------------------------------------------------------------------

def dominating_vertices(adj):
    n = len(adj)
    return [v for v in range(n) if len(adj[v]) == n - 1]


def in_equality_family(adj, k, min_components=None):
    """Is adj isomorphic to K_k v H with H on n-k vertices having at least
    `min_components` components?  Default min_components = k+1, the clause as
    quoted in the note.

    A vertex of the K_k part must have degree n-1, so the K_k part is a
    k-subset of the dominating set D.  Any two vertices of D are twins
    (their closed neighbourhoods are all of V), so the choice inside D is
    immaterial and taking the first k of D is without loss of generality.
    """
    n = len(adj)
    if min_components is None:
        min_components = k + 1
    if k < 1 or k > n - 1:
        return False
    D = dominating_vertices(adj)
    if len(D) < k:
        return False
    rest = set(range(n)) - set(D[:k])
    return components(adj, rest) >= min_components


def characterisation_holds(adj, k, min_components=None):
    return (in_equality_family(adj, k, min_components)
            or in_equality_family(complement(adj), k, min_components))


# ----------------------------------------------------------------------
# graph6
# ----------------------------------------------------------------------

def graph6_decode(s):
    data = [ord(c) - 63 for c in s]
    n = data[0]
    if not 0 <= n <= 62:
        raise ValueError("only graph6 orders up to 62 are handled")
    bits = []
    for b in data[1:]:
        for i in range(5, -1, -1):
            bits.append((b >> i) & 1)
    adj = [set() for _ in range(n)]
    pos = 0
    for j in range(1, n):
        for i in range(j):
            if pos < len(bits) and bits[pos]:
                adj[i].add(j)
                adj[j].add(i)
            pos += 1
    return adj


def graph6_encode(adj):
    n = len(adj)
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if j in adj[i] else 0)
    while len(bits) % 6:
        bits.append(0)
    out = chr(n + 63)
    for i in range(0, len(bits), 6):
        v = 0
        for b in bits[i:i + 6]:
            v = 2 * v + b
        out += chr(v + 63)
    return out


# ----------------------------------------------------------------------
# isomorph-free generation for small n (pure Python, orbit marking)
# ----------------------------------------------------------------------

def iso_classes(n):
    """One representative per isomorphism class of simple graphs on n vertices.

    Every labelled graph is encoded as a bitmask over the C(n,2) vertex pairs;
    the orbit of an unmarked mask under the full symmetric group is generated
    and marked, so each class is emitted exactly once.
    """
    pairs = list(combinations(range(n), 2))
    E = len(pairs)
    index = dict((p, i) for i, p in enumerate(pairs))
    maps = []
    for p in permutations(range(n)):
        m = [0] * E
        for i, (u, v) in enumerate(pairs):
            a, b = p[u], p[v]
            m[i] = index[(a, b) if a < b else (b, a)]
        maps.append(m)
    seen = bytearray(1 << E)
    reps = []
    for mask in range(1 << E):
        if seen[mask]:
            continue
        reps.append(mask)
        for m in maps:
            img = 0
            mm = mask
            i = 0
            while mm:
                if mm & 1:
                    img |= 1 << m[i]
                mm >>= 1
                i += 1
            seen[img] = 1
    out = []
    for mask in reps:
        adj = [set() for _ in range(n)]
        for i, (u, v) in enumerate(pairs):
            if mask >> i & 1:
                adj[u].add(v)
                adj[v].add(u)
        out.append(adj)
    return out


def join(a1, a2):
    n1, n2 = len(a1), len(a2)
    adj = [set(a1[u]) for u in range(n1)]
    adj += [set(w + n1 for w in a2[u]) for u in range(n2)]
    for u in range(n1):
        for v in range(n1, n1 + n2):
            adj[u].add(v)
            adj[v].add(u)
    return adj


def complete_graph(m):
    return disjoint_cliques([m]) if m else []


def empty_graph(m):
    return disjoint_cliques([1] * m)


def spec_list(claimed):
    out = []
    for lam, mult in claimed:
        out += [Fraction(lam)] * mult
    out.sort(reverse=True)
    return out


# ----------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------

T_RANGE = [2, 3, 4, 5]          # the family K_{2t,2t}, t = 2..5
T_CHARPOLY = [2, 3]             # cross-validated against characteristic polys
CENSUS_MAX_N = 6                # isomorph-free generation done for 2..6
PUBLISHED_G6 = ["G?zvf_", "GQhTUg", "GCXnf_"]


def main():
    print("Verifying the counterexample family K_{2t,2t} at n = 4t, k = 3t,")
    print("for t = %s, plus the five named eight-vertex graphs and an"
          % ", ".join(str(t) for t in T_RANGE))
    print("isomorph-free census for 2 <= n <= %d." % CENSUS_MAX_N)
    print("")

    fam = {}
    for t in T_RANGE:
        n, k, m = 4 * t, 3 * t, 2 * t
        G = complete_bipartite(m, m)
        Gc = complement(G)
        fam[t] = (n, k, m, G, Gc)

    # ---- 1. complement structure -------------------------------------
    bad = []
    for t in T_RANGE:
        n, k, m, G, Gc = fam[t]
        if edge_set(Gc) != edge_set(disjoint_cliques([m, m])):
            bad.append(t)
    ck("complement of K_{2t,2t} is the disjoint union of two copies of K_{2t}",
       not bad,
       "checked as labelled edge sets for t = %s"
       % ", ".join(str(t) for t in T_RANGE)
       if not bad else "mismatch at t = %s" % bad)

    # ---- 2. Laplacian spectrum of G ----------------------------------
    bad = []
    for t in T_RANGE:
        n, k, m, G, Gc = fam[t]
        ok, det = spectrum_is(laplacian(G),
                              [(4 * t, 1), (2 * t, 4 * t - 2), (0, 1)])
        if not ok:
            bad.append("t=%d: %s" % (t, det))
    ck("Spec_L(K_{2t,2t}) = {4t, (2t)^(4t-2), 0} by exact nullity of L - lambda I",
       not bad,
       "n = %s all confirmed" % ", ".join(str(4 * t) for t in T_RANGE)
       if not bad else "; ".join(bad))

    # ---- 3. Laplacian spectrum of the complement ---------------------
    bad = []
    for t in T_RANGE:
        n, k, m, G, Gc = fam[t]
        ok, det = spectrum_is(laplacian(Gc), [(2 * t, 4 * t - 2), (0, 2)])
        if not ok:
            bad.append("t=%d: %s" % (t, det))
    ck("Spec_L(2K_{2t}) = {(2t)^(4t-2), 0^2} by exact nullity of L - lambda I",
       not bad,
       "n = %s all confirmed" % ", ".join(str(4 * t) for t in T_RANGE)
       if not bad else "; ".join(bad))

    # ---- 4. the k-th eigenvalues -------------------------------------
    bad = []
    for t in T_RANGE:
        n, k, m, G, Gc = fam[t]
        eg = spec_list([(4 * t, 1), (2 * t, 4 * t - 2), (0, 1)])
        ec = spec_list([(2 * t, 4 * t - 2), (0, 2)])
        if not (len(eg) == len(ec) == n and 3 * t <= 4 * t - 2
                and eg[k - 1] == 2 * t and ec[k - 1] == 2 * t):
            bad.append(t)
    ck("mu_{3t}(G) = mu_{3t}(complement G) = 2t (index 3t lies in the block of "
       "multiplicity 4t-2)", not bad,
       "verified for t = %s" % ", ".join(str(t) for t in T_RANGE)
       if not bad else "fails at t = %s" % bad)

    # ---- 5. equality in the product ----------------------------------
    bad = []
    for t in T_RANGE:
        n, k, m, G, Gc = fam[t]
        eg = spec_list([(4 * t, 1), (2 * t, 4 * t - 2), (0, 1)])
        ec = spec_list([(2 * t, 4 * t - 2), (0, 2)])
        if eg[k - 1] * ec[k - 1] != n * (n - k) or 4 * t * t != n * (n - k):
            bad.append(t)
    ck("mu_k(G) mu_k(complement G) = 4t^2 = n(n-k) exactly at n = 4t, k = 3t",
       not bad,
       "products %s" % ", ".join("t=%d: %d" % (t, 4 * t * t) for t in T_RANGE)
       if not bad else "fails at t = %s" % bad)

    # ---- 6. k lies in the conjecture's range -------------------------
    bad = []
    for t in T_RANGE:
        n, k, m, G, Gc = fam[t]
        if not (1 <= k and Fraction(4 * k, 3) == n and 2 * k >= n):
            bad.append(t)
    ck("k = 3t satisfies 1 <= k <= 3n/4 (the endpoint k = 3n/4) and n/2 <= k, "
       "so it lies in the range where the inequality is already proved",
       not bad,
       "(n,k) = %s" % ", ".join("(%d,%d)" % (4 * t, 3 * t) for t in T_RANGE)
       if not bad else "fails at t = %s" % bad)

    # ---- 7. the proposed equality family is empty ---------------------
    bad = []
    for t in T_RANGE:
        n, k, m, G, Gc = fam[t]
        # a graph on n-k = t vertices has at most t components, and the clause
        # demands at least k+1 = 3t+1 of them
        if not (n - k == t and k + 1 == 3 * t + 1 and 3 * t + 1 > t):
            bad.append("t=%d arithmetic" % t)
        if in_equality_family(G, k) or in_equality_family(Gc, k):
            bad.append("t=%d predicate" % t)
    ck("the conjectured equality family is empty at (n,k) = (4t,3t): H would "
       "need 3t+1 components on t vertices", not bad,
       "n-k = t vertices vs k+1 = 3t+1 required components, t = %s"
       % ", ".join(str(t) for t in T_RANGE) if not bad else "; ".join(bad))

    # ---- 8. robustness of the refutation to the wording of the clause -
    bad = []
    for t in T_RANGE:
        n, k, m, G, Gc = fam[t]
        if max(degrees(G)) != 2 * t or max(degrees(Gc)) != 2 * t - 1:
            bad.append("t=%d degrees" % t)
        if max(degrees(G)) >= n - 1 or max(degrees(Gc)) >= n - 1:
            bad.append("t=%d dominating vertex" % t)
        for j in range(1, n):
            if in_equality_family(G, j, 2) or in_equality_family(Gc, j, 2):
                bad.append("t=%d j=%d" % (t, j))
    ck("neither K_{2t,2t} nor its complement is a join K_j v H for ANY j >= 1, "
       "since neither has a vertex of degree n-1, so the refutation survives "
       "any weaker component requirement (>= 2 instead of >= k+1)", not bad,
       "max degrees 2t and 2t-1, both < n-1 = 4t-1, for t = %s"
       % ", ".join(str(t) for t in T_RANGE) if not bad else "; ".join(bad))

    # ---- 9. independent route: characteristic polynomial + Descartes ---
    bad = []
    for t in T_CHARPOLY:
        n, k, m, G, Gc = fam[t]
        for adj, claimed in ((G, [(4 * t, 1), (2 * t, 4 * t - 2), (0, 1)]),
                             (Gc, [(2 * t, 4 * t - 2), (0, 2)])):
            p = charpoly(laplacian(adj))
            if any(c.denominator != 1 for c in p):
                bad.append("t=%d non-integral characteristic polynomial" % t)
            for lam, mult in claimed:
                if mult_at(p, lam) != mult:
                    bad.append("t=%d root multiplicity at %d" % (t, lam))
            tot = sum(mult for _, mult in claimed)
            if count_gt(p, Fraction(0)) != tot - dict(claimed)[0]:
                bad.append("t=%d positive-eigenvalue count" % t)
            if mu_exact(p, n, k) != 2 * t:
                bad.append("t=%d mu_k from the polynomial" % t)
    ck("root multiplicities of the integer characteristic polynomial "
       "independently reproduce the two spectra and give mu_{3t} = 2t",
       not bad,
       "exact Taylor shift + Descartes' rule (exact for real-rooted "
       "polynomials), t = %s" % ", ".join(str(t) for t in T_CHARPOLY)
       if not bad else "; ".join(sorted(set(bad))))

    # ---- 10. the displayed eight-vertex instance ----------------------
    K44 = complete_bipartite(4, 4)
    K44c = complement(K44)
    p44 = charpoly(laplacian(K44))
    p44c = charpoly(laplacian(K44c))
    a = mu_exact(p44, 8, 6)
    b = mu_exact(p44c, 8, 6)
    ck("mu_6(K_{4,4}) mu_6(2K_4) = 4 * 4 = 16 = 8(8-6), the displayed t = 2 "
       "instance",
       a == 4 and b == 4 and a * b == 16 == 8 * (8 - 6)
       and edge_set(K44c) == edge_set(disjoint_cliques([4, 4])),
       "mu_6 = %s and %s" % (a, b))

    # ---- 11. graph6 decoding of the three published codes -------------
    decoded = {}
    bad = []
    for code in PUBLISHED_G6:
        g = graph6_decode(code)
        decoded[code] = g
        if len(g) != 8:
            bad.append("%s has order %d" % (code, len(g)))
        if graph6_encode(g) != code:
            bad.append("%s does not re-encode to itself" % code)
    ck("the three graph6 codes quoted in the closing paragraph decode to "
       "eight-vertex graphs and re-encode to the same strings", not bad,
       "%s" % ", ".join("%s: %d edges" % (c, len(edge_set(decoded[c])))
                        for c in PUBLISHED_G6)
       if not bad else "; ".join(bad))

    # ---- 12. all five named eight-vertex graphs attain equality --------
    named = [("K_{4,4}", K44), ("2K_4", K44c)]
    named += [(c, decoded[c]) for c in PUBLISHED_G6]
    bad = []
    for label, g in named:
        pg = charpoly(laplacian(g))
        pc = charpoly(laplacian(complement(g)))
        r = product_vs(pg, pc, 8, 6, 16)
        aa, bb = mu_exact(pg, 8, 6), mu_exact(pc, 8, 6)
        if r != "eq":
            bad.append("%s: %s (mu_6 = %s, %s)" % (label, r, aa, bb))
    ck("each of the five graphs named for n = 8, k = 6 satisfies "
       "mu_6(G) mu_6(complement G) = 16 = n(n-k) exactly", not bad,
       "verified for %s" % ", ".join(lbl for lbl, _ in named)
       if not bad else "; ".join(bad))

    # ---- 13. the equality characterisation fails on all five ----------
    bad = []
    for label, g in named:
        if characterisation_holds(g, 6):
            bad.append(label)
    ck("for each of those five graphs neither G nor complement G is K_6 v H "
       "with H having at least k+1 = 7 components on n-k = 2 vertices, so the "
       "'only if' direction of Conjecture 20 fails", not bad,
       "characterisation fails for all five" if not bad
       else "unexpectedly in the family: %s" % ", ".join(bad))

    # ---- 13b. the stated structure of those three graphs ---------------
    g1, g2, g3 = (decoded[c] for c in PUBLISHED_G6)
    ok13b = (isomorphic(g3, complement(g3))
             and isomorphic(g1, complement(g2))
             and not isomorphic(g1, g3))
    ck("the closing paragraph's structural claims about the three codes hold: "
       "GCXnf_ is self-complementary and the other two form a complementary "
       "pair", ok13b,
       "orders %s and edge counts %s"
       % (", ".join("%d" % len(decoded[c]) for c in PUBLISHED_G6),
          ", ".join("%d" % len(edge_set(decoded[c])) for c in PUBLISHED_G6)))

    # ---- 14. the failure survives the weakest reading of the clause ----
    K8e = join(complete_graph(6), empty_graph(2))
    p8 = charpoly(laplacian(K8e))
    p8c = charpoly(laplacian(complement(K8e)))
    r1 = product_vs(p8, p8c, 8, 6, 16)
    # K_6 v 2K_1 is K_8 with the single edge {6, 7} removed: the two vertices
    # of the 2K_1 part are joined to all of the K_6 part but not to each other.
    # That identity is asserted edge by edge, and the count 28 - 1 = 27 follows;
    # len(edge_set(.)) is already the number of edges, so it is not halved.
    e14 = len(edge_set(K8e))
    ok14 = (edge_set(K8e) == edge_set(complete_graph(8)) - {frozenset((6, 7))}
            and e14 == 27
            and in_equality_family(K8e, 6, 2)
            and mu_exact(p8c, 8, 6) == 0
            and r1 == "lt")
    ck("weakest reading at n = 8, k = 6: H on 2 vertices with at least 2 "
       "components forces H = 2K_1, so the only candidate is K_6 v 2K_1 = "
       "K_8 minus an edge, and its product is 0, not 16", ok14,
       "K_8 - e has %d edges, mu_6(complement) = %s, comparison %s"
       % (e14, mu_exact(p8c, 8, 6), r1))

    # ---- 15. positive controls for the two predicates ------------------
    P = join(complete_graph(2), empty_graph(3))          # K_2 v 3K_1, n = 5
    pP = charpoly(laplacian(P))
    pPc = charpoly(laplacian(complement(P)))
    ok15 = (in_equality_family(P, 2)                     # 3 >= k+1 = 3
            and product_vs(pP, pPc, 5, 2, 5 * 3) == "eq"
            and mu_exact(pP, 5, 2) == 5 and mu_exact(pPc, 5, 2) == 3
            and not in_equality_family(P, 3)
            and not in_equality_family(K8e, 6))
    ck("positive control: the family predicate accepts K_2 v 3K_1 at k = 2 "
       "(and that graph does attain equality, mu_2 = 5 and 3, product 15 = "
       "5(5-2)), while rejecting cases whose component count is too small",
       ok15,
       "the predicate and the equality test are not vacuously false")

    # ---- 16..18. census for 2 <= n <= CENSUS_MAX_N --------------------
    known_counts = {2: 2, 3: 4, 4: 11, 5: 34, 6: 156}
    counts = {}
    pairs = 0
    violations = []
    undecided = []
    equalities = []
    unexplained = []
    for n in range(2, CENSUS_MAX_N + 1):
        reps = iso_classes(n)
        counts[n] = len(reps)
        for g in reps:
            pg = charpoly(laplacian(g))
            pc = charpoly(laplacian(complement(g)))
            for k in range(1, (3 * n) // 4 + 1):
                pairs += 1
                r = product_vs(pg, pc, n, k, n * (n - k))
                if r == "gt":
                    violations.append((n, k, graph6_encode(g)))
                elif r == "undecided":
                    undecided.append((n, k, graph6_encode(g)))
                elif r == "eq":
                    equalities.append((n, k, graph6_encode(g)))
                    if not characterisation_holds(g, k):
                        unexplained.append((n, k, graph6_encode(g)))
    print("")
    print("census: %s isomorphism classes for n = 2..%d, %d (graph, k) pairs "
          "evaluated," % (", ".join("%d:%d" % (n, counts[n])
                                    for n in sorted(counts)),
                          CENSUS_MAX_N, pairs))
    print("        %d equality cases found, %d of them outside the "
          "conjectured family." % (len(equalities), len(unexplained)))
    print("")

    ck("isomorph-free generation reproduces the known numbers of graphs up to "
       "isomorphism for n = 2..%d (2, 4, 11, 34, 156)" % CENSUS_MAX_N,
       all(counts.get(n) == known_counts[n]
           for n in range(2, CENSUS_MAX_N + 1)),
       "counts %s" % ", ".join("%d:%d" % (n, counts[n])
                               for n in sorted(counts)))

    ck("no graph with 2 <= n <= %d violates the product inequality "
       "mu_k(G) mu_k(complement G) <= n(n-k) for 1 <= k <= floor(3n/4)"
       % CENSUS_MAX_N,
       not violations and not undecided,
       "%d pairs evaluated, %d equality cases, no violation"
       % (pairs, len(equalities)) if not violations and not undecided
       else "violations %s; undecided %s" % (violations[:5], undecided[:5]))

    ck("every equality case with 2 <= n <= %d does satisfy the conjectured "
       "characterisation, so there is no counterexample of order at most %d "
       "and none of the five n = 8 graphs is beaten by a smaller one in this "
       "range" % (CENSUS_MAX_N, CENSUS_MAX_N),
       not unexplained,
       "%d equality cases, all explained by the K_k v H description"
       % len(equalities) if not unexplained
       else "unexplained: %s" % unexplained[:5])

    ck("the census actually evaluated a non-trivial number of (graph, k) "
       "pairs and every one of them was decided exactly",
       pairs >= 400 and not undecided,
       "%d pairs, 0 undecided" % pairs)

    print("")
    print("COVERAGE OF THE PUBLISHED RANGE 2 <= n <= 9: covered exhaustively "
          "here for n = 2, 3, 4, 5, 6 only.")
    print("NOT COVERED HERE: the orders n = 7, n = 8 and n = 9 were NOT "
          "enumerated (nauty was unavailable and pure-Python isomorph-free")
    print("generation was kept to n <= 6); the five named eight-vertex graphs "
          "were each checked individually, but their EXHAUSTIVENESS")
    print("('exactly five cases') and the minimality of K_{4,4} were not.")
    print("")

    rc = verdict()
    print("NOT RE-RUN HERE: the exhaustive census over n = 7, 8, 9 (no "
          "isomorph-free generation was performed above n = 6 and nauty/geng "
          "was not used), hence the published tally 'exactly five cases, all "
          "at n = 8, k = 6' is NOT confirmed as exhaustive, and the five named "
          "graphs are not checked to be pairwise non-isomorphic either (only "
          "that GCXnf_ is self-complementary, that G?zvf_ and GQhTUg are "
          "complementary, and that G?zvf_ and GCXnf_ differ), so the count "
          "five is not confirmed from below here either; the minimality "
          "claim that K_{4,4} has smallest order among counterexamples is "
          "confirmed only against orders n <= 6; the product inequality "
          "mu_k(G) mu_k(complement G) <= n(n-k) itself is nowhere proved here, "
          "only evaluated on the graphs enumerated or named above (and for "
          "1 <= k < n/2 it is untouched by the family, which lives at "
          "k = 3n/4); the family K_{2t,2t} is verified for t = 2..5 only, the "
          "general t being closed-form; the transcription of Conjecture 20, of "
          "Theorem 18 and of Remark 1 from the published article was not "
          "checked against the published text (no external source was "
          "consulted, so the equality clause is taken exactly as quoted in the "
          "note -- the refutation is however also checked to survive the "
          "weakest reading, at least 2 components instead of at least k+1); "
          "and the remark that Chen, Guo, Li and Wang had verified the "
          "conjecture for all graphs on at most nine vertices is reported, not "
          "examined, here.")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

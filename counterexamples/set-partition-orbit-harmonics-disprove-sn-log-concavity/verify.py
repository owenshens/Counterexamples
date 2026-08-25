#!/usr/bin/env python3
"""Verification of the counterexamples to the S_n-log-concavity conjecture
for set-partition orbit harmonics.

Setting.  For 1 <= m <= n let Pi_{n,m} be the set partitions of [n] with all
blocks of size <= m.  The degree-d piece of the orbit-harmonics quotient is
isomorphic, as an S_n-module, to the permutation module C[X_{n,m,d}] where
X_{n,m,d} = {P in Pi_{n,m} : |P| = n-d}.  Log-concavity in degree d would need
an S_n-injection C[X_{d-1}] tensor C[X_{d+1}] --> C[X_d] tensor C[X_d], hence
[sgn : source] <= [sgn : target] for the sign representation.

TAKEN FROM THE PAPER (inputs):
  * the family  k odd, k >= 5,  n = k^2 - 5,  d = n - k,  k+1 <= m <= n;
  * the asserted sign multiplicities: 4 in the source; 1 in the target for
    k = 5 and 0 in the target for k >= 7;
  * the four exhibited four-edge complement graphs
        C_4,  P_5^{(3,2)},  P_5^{(2,3)},  P_4 + K_2
    and the single exhibited five-edge complement graph
        K_{1,2} + K_{2,1} + K_2,
    entered below as literal edge lists;
  * the orbit criterion (Lemma 2 of the paper: [sgn : C[Y]] counts the orbits
    whose point stabiliser lies in A_n) and the grid-sign identity
        sgn_{E(G)} * sgn_{E(Z)} = sgn(sigma)^c * sgn(pi)^r
    for Z the complement of G inside K_{r,c}.

DERIVED HERE (checks):
  * the arithmetic of the family: |P|, |Q|, and the complement edge counts
    4 and 5, for many odd k;
  * that the exhibited edge lists really are the named graphs (edge count,
    no isolated vertex, degree sequences, connectivity, max degree 2) and,
    since those invariants do NOT determine the isomorphism type -- both
    P_4 + K_2 and K_{1,2} + K_{2,1} have 4 edges, sides (3,3), degrees
    (2,1,1) on each side and 2 components, and both K_{1,2}+K_{2,1}+K_2 and
    P_4 + 2K_2 have 5 edges, sides (4,4), degrees (2,1,1,1) on each side and
    3 components -- a per-name canonical-form match against each named graph
    rebuilt from its name out of paths, complete bipartite pieces and
    disjoint unions;
  * an independent exhaustive census, up to side-preserving isomorphism, of
    all four-edge and all five-edge bipartite complements, with each
    candidate's full side-preserving automorphism group recomputed and the
    sign condition evaluated on it.  The census returns exactly the paper's
    lists and exactly the multiplicities 4 / 1 / 0;
  * a from-scratch character-theoretic computation of the sign multiplicity
    in C[X_a] tensor C[X_b] for small n (summing sgn(g)|X_a^g||X_b^g| over
    S_n by cycle type, in exact integer arithmetic) which is compared with
    the graph-orbit census on every instance where both are available.  This
    validates the machinery used at n = 20 on cases small enough to brute
    force directly from the definition of Pi_{n,m}.
"""

import sys
from itertools import combinations, permutations

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return bool(ok)


# ----------------------------------------------------------------------
# Objects taken from the paper, as literal edge lists (left, right).
# ----------------------------------------------------------------------

SOURCE_COMPLEMENTS = {
    "C_4": [(0, 0), (0, 1), (1, 0), (1, 1)],
    "P_5^(3,2)": [(0, 0), (1, 0), (1, 1), (2, 1)],
    "P_5^(2,3)": [(0, 0), (0, 1), (1, 1), (1, 2)],
    "P_4+K_2": [(0, 0), (1, 0), (1, 1), (2, 2)],
}

TARGET_COMPLEMENT = {
    "K_{1,2}+K_{2,1}+K_2": [(0, 0), (0, 1), (1, 2), (2, 2), (3, 3)],
}

PAPER_SOURCE_MULT = 4
PAPER_TARGET_MULT = {5: 1}     # k = 5
PAPER_TARGET_MULT_LARGE = 0    # k >= 7

# The census enumerates complements exhaustively; the cases relevant to the
# theorem need only 4 and 5 edges.  Larger complements are refused outright
# rather than silently attempted.
MAX_COMPLEMENT_EDGES = 6

# Total number of checks main() must record, including the inventory check
# itself.  Kept explicit so a skipped check fails instead of vanishing.
EXPECTED_CHECKS = 32


# ----------------------------------------------------------------------
# The same graphs rebuilt from their NAMES, independently of the literal
# edge lists above, out of paths / complete bipartite pieces / disjoint
# unions.  The invariants checked below (edge count, bipartition sizes,
# degree sequences, component count) do not pin an isomorphism class:
# P_4 + K_2 and K_{1,2} + K_{2,1} share all of them, and so do
# K_{1,2} + K_{2,1} + K_2 and P_4 + 2K_2.  Comparing canonical forms
# name by name against these reconstructions does pin them.
# ----------------------------------------------------------------------

def complete_bipartite(a, b):
    """K_{a,b} as a sorted edge tuple."""
    return tuple(sorted((i, j) for i in range(a) for j in range(b)))


def bipartite_path(j, start_left):
    """The path on j vertices, sides alternating, first vertex on the left
    iff start_left.  Vertices are numbered along the path within each side."""
    verts, li, ri = [], 0, 0
    for t in range(j):
        if (t % 2 == 0) == bool(start_left):
            verts.append(("L", li))
            li += 1
        else:
            verts.append(("R", ri))
            ri += 1
    edges = []
    for t in range(j - 1):
        u, v = verts[t], verts[t + 1]
        edges.append((u[1], v[1]) if u[0] == "L" else (v[1], u[1]))
    return tuple(sorted(edges))


def disjoint_union(*parts):
    """Vertex-disjoint union, shifting each side's labels past the previous
    parts.  Every part must be nonempty and have no isolated vertex."""
    edges, dl, dc = [], 0, 0
    for p in parts:
        edges.extend((i + dl, j + dc) for (i, j) in p)
        dl += 1 + max(i for i, _ in p)
        dc += 1 + max(j for _, j in p)
    return tuple(sorted(edges))


NAMED_SOURCE = {
    "C_4": complete_bipartite(2, 2),
    "P_5^(3,2)": bipartite_path(5, True),
    "P_5^(2,3)": bipartite_path(5, False),
    "P_4+K_2": disjoint_union(bipartite_path(4, True),
                              complete_bipartite(1, 1)),
}

NAMED_TARGET = {
    "K_{1,2}+K_{2,1}+K_2": disjoint_union(complete_bipartite(1, 2),
                                          complete_bipartite(2, 1),
                                          complete_bipartite(1, 1)),
}


def used_sides(edges):
    return (len(set(i for i, _ in edges)), len(set(j for _, j in edges)))


def degrees(edges):
    left, right = {}, {}
    for i, j in edges:
        left[i] = left.get(i, 0) + 1
        right[j] = right.get(j, 0) + 1
    return sorted(left.values(), reverse=True), sorted(right.values(), reverse=True)


def components(edges):
    """Number of connected components of the graph on its non-isolated vertices."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in edges:
        for v in (("L", i), ("R", j)):
            parent.setdefault(v, v)
    for i, j in edges:
        a, b = find(("L", i)), find(("R", j))
        if a != b:
            parent[a] = b
    return len(set(find(v) for v in parent))


def check_exhibited_objects():
    """The literal edge lists are the graphs the paper names."""
    ok = True
    for name, edges in sorted(SOURCE_COMPLEMENTS.items()):
        if len(set(edges)) != 4 or len(edges) != 4:
            ok = False
    ck("source-complements-have-four-edges", ok,
       "4 graphs, 4 distinct edges each")

    sides = {name: used_sides(e) for name, e in SOURCE_COMPLEMENTS.items()}
    ck("source-complement-bipartition-sizes",
       sides == {"C_4": (2, 2), "P_5^(3,2)": (3, 2),
                 "P_5^(2,3)": (2, 3), "P_4+K_2": (3, 3)},
       str(sorted(sides.items())))

    degs = {name: degrees(e) for name, e in SOURCE_COMPLEMENTS.items()}
    ck("source-complement-degree-sequences",
       degs == {"C_4": ([2, 2], [2, 2]),
                "P_5^(3,2)": ([2, 1, 1], [2, 2]),
                "P_5^(2,3)": ([2, 2], [2, 1, 1]),
                "P_4+K_2": ([2, 1, 1], [2, 1, 1])},
       "max degree 2 throughout")

    comps = {name: components(e) for name, e in SOURCE_COMPLEMENTS.items()}
    ck("source-complement-component-counts",
       comps == {"C_4": 1, "P_5^(3,2)": 1, "P_5^(2,3)": 1, "P_4+K_2": 2},
       "three connected, one with two components")

    (name, edges), = TARGET_COMPLEMENT.items()
    ck("target-complement-shape",
       len(edges) == 5 and len(set(edges)) == 5
       and used_sides(edges) == (4, 4)
       and degrees(edges) == ([2, 1, 1, 1], [2, 1, 1, 1])
       and components(edges) == 3,
       name + ": 5 edges on (4,4) vertices, degrees (2,1,1,1) both sides, 3 parts")

    # The invariants above are not a complete isomorphism invariant, and the
    # ambiguity is not hypothetical: K_{1,2} + K_{2,1} matches every source
    # test applied to P_4 + K_2, and P_4 + 2K_2 matches every target test,
    # while both of those graphs FAIL the paper's sign condition.  So compare
    # canonical forms, name by name, with the graphs rebuilt from their names.
    got_src = {nm: canon(e) for nm, e in SOURCE_COMPLEMENTS.items()}
    want_src = {nm: canon(e) for nm, e in NAMED_SOURCE.items()}
    ck("source-complements-are-the-named-graphs", got_src == want_src,
       "each literal edge list is side-preserving isomorphic to the graph its "
       "name denotes, rebuilt from paths / complete bipartite pieces"
       if got_src == want_src
       else "MISMATCH at " + str(sorted(nm for nm in set(got_src) | set(want_src)
                                        if got_src.get(nm) != want_src.get(nm))))

    got_tgt = {nm: canon(e) for nm, e in TARGET_COMPLEMENT.items()}
    want_tgt = {nm: canon(e) for nm, e in NAMED_TARGET.items()}
    ck("target-complement-is-the-named-graph", got_tgt == want_tgt,
       "the literal edge list is side-preserving isomorphic to "
       "K_{1,2} + K_{2,1} + K_2 rebuilt from its name"
       if got_tgt == want_tgt else "MISMATCH")


def check_family_arithmetic(ks):
    """n = k^2-5, d = n-k gives block counts k+1, k-1, k and complements 4, 5."""
    bad = []
    for k in ks:
        n = k * k - 5
        d = n - k
        r_src, c_src = n - (d - 1), n - (d + 1)
        r_tgt = c_tgt = n - d
        if (r_src, c_src) != (k + 1, k - 1):
            bad.append((k, "src sides", r_src, c_src))
        if r_src * c_src - n != 4:
            bad.append((k, "src complement", r_src * c_src - n))
        if (r_tgt, c_tgt) != (k, k):
            bad.append((k, "tgt sides", r_tgt, c_tgt))
        if r_tgt * c_tgt - n != 5:
            bad.append((k, "tgt complement", r_tgt * c_tgt - n))
        if not (1 <= d - 1 and d + 1 <= n):
            bad.append((k, "degree range", d))
    ck("family-arithmetic", not bad,
       "k in %s: sides (k+1,k-1) and (k,k); |E(Z)| = 4 and 5%s"
       % (list(ks), "" if not bad else " VIOLATIONS " + str(bad)))

def count_partitions(n, blocks, m):
    """Set partitions of [n] into exactly `blocks` blocks, all of size <= m."""
    total = 0
    for lam in integer_partitions(n, m):
        if len(lam) != blocks:
            continue
        term = factorial(n)
        counts = {}
        for part in lam:
            term //= factorial(part)
            counts[part] = counts.get(part, 0) + 1
        for mult in counts.values():
            term //= factorial(mult)
        total += term
    return total


def check_first_member_dimensions():
    """At k=5 the three graded pieces in play are genuinely nonzero."""
    k, n, d, m = 5, 20, 15, 6
    # Recorded as a check rather than an early return: an early return here
    # would silently drop the three checks below and shrink the reported
    # check total instead of failing.
    ck("first-member-is-n-20", (k * k - 5, (k * k - 5) - k) == (n, d),
       "k = 5 gives n = k^2-5 = 20 and d = n-k = 15")
    dims = [count_partitions(n, n - deg, m) for deg in (d - 1, d, d + 1)]
    ck("first-member-graded-pieces-nonzero",
       all(x > 0 for x in dims) and len(dims) == 3,
       "n=20, m=6: dim R_14 = %d (6 blocks), dim R_15 = %d (5 blocks), "
       "dim R_16 = %d (4 blocks)" % tuple(dims))
    closed_form = factorial(20) // (factorial(5) ** 4 * factorial(4))
    ck("block-count-cap-is-respected",
       count_partitions(n, 4, 4) == 0
       and count_partitions(n, 3, 6) == 0
       and count_partitions(n, 4, 5) == closed_form,
       "4 blocks of size <= 4 and 3 of size <= 6 cannot cover 20 elements; "
       "4 blocks of size <= 5 gives %d, matching 20!/(5!^4 4!)" % closed_form)
    bell = [len(set_partitions(t)) for t in range(0, 8)]
    recomputed = [sum(count_partitions(t, b, t) for b in range(0, t + 1))
                  for t in range(0, 8)]
    ck("partition-counter-reproduces-bell-numbers",
       bell == recomputed and bell[7] == 877,
       "counts %s agree with direct enumeration of set partitions" % (bell,))


# ----------------------------------------------------------------------
# Exhaustive census of complements Z, up to side-preserving isomorphism.
# ----------------------------------------------------------------------

def parity(perm):
    """Sign of a permutation given as a tuple perm[i] = image of i."""
    seen = [False] * len(perm)
    swaps = 0
    for i in range(len(perm)):
        if not seen[i]:
            j = i
            length = 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                length += 1
            swaps += length - 1
    return -1 if swaps % 2 else 1


def edge_sign(edges_sorted, index, sigma, pi):
    """Sign of the permutation of E(Z) induced by (sigma, pi)."""
    image = [index[(sigma[i], pi[j])] for (i, j) in edges_sorted]
    return parity(tuple(image))


def core_automorphisms(edges, r0, c0):
    """All side-preserving automorphisms of Z restricted to its used vertices."""
    eset = frozenset(edges)
    out = []
    for sigma in permutations(range(r0)):
        moved = frozenset((sigma[i], j) for (i, j) in edges)
        for pi in permutations(range(c0)):
            if frozenset((i, pi[j]) for (i, j) in moved) == eset:
                out.append((sigma, pi))
    return out


def orbit_of(edges, r0, c0):
    """Every image of the edge set under side-preserving relabellings."""
    out = set()
    for sigma in permutations(range(r0)):
        moved = [(sigma[i], j) for (i, j) in edges]
        for pi in permutations(range(c0)):
            out.add(tuple(sorted((i, pi[j]) for (i, j) in moved)))
    return out


_CORE_CLASS_CACHE = {}


def core_classes(e):
    """Every e-edge bipartite graph with no isolated vertex, one per iso class.

    Returned as (canonical edge tuple, r0, c0); r0, c0 <= e always.  Each
    orbit is expanded once, so the canonical representative is provably the
    lexicographic minimum of the orbit and distinct classes never collide.

    Memoised: the answer depends only on e, so the wide-k sweeps below cost
    nothing beyond the first enumeration.  The cache is keyed on e and the
    returned list is never mutated by callers.
    """
    if e in _CORE_CLASS_CACHE:
        return _CORE_CLASS_CACHE[e]
    classes = []
    for r0 in range(1, e + 1):
        for c0 in range(1, e + 1):
            if r0 * c0 < e:
                continue
            cells = [(i, j) for i in range(r0) for j in range(c0)]
            seen = set()
            for sub in combinations(cells, e):
                if sub in seen:
                    continue
                if len(set(i for i, _ in sub)) != r0:
                    continue
                if len(set(j for _, j in sub)) != c0:
                    continue
                orb = orbit_of(sub, r0, c0)
                seen |= orb
                classes.append((min(orb), r0, c0))
    _CORE_CLASS_CACHE[e] = classes
    return classes


def contributes(edges, r0, c0, r, c):
    """Does the orbit of G = K_{r,c} \\ Z have its stabiliser inside A_n?

    Condition: sgn_{E(Z)}(sigma,pi) = sgn(sigma)^c sgn(pi)^r for every
    side-preserving automorphism of Z inside K_{r,c}.  Aut(Z) is
    Aut(Z-core) times the symmetric groups on the isolated vertices, and
    both sides of the condition are homomorphisms to {+1,-1}, so it is
    enough to test Aut(Z-core) and one transposition of isolated vertices
    on each side.
    """
    if r0 > r or c0 > c:
        return False
    es = sorted(edges)
    index = {cell: t for t, cell in enumerate(es)}
    for sigma, pi in core_automorphisms(es, r0, c0):
        lhs = edge_sign(es, index, sigma, pi)
        rhs = (parity(sigma) ** c) * (parity(pi) ** r)
        if lhs != rhs:
            return False
    if r - r0 >= 2 and c % 2 == 1:
        return False          # swap two isolated left vertices: 1 != (-1)^c
    if c - c0 >= 2 and r % 2 == 1:
        return False          # swap two isolated right vertices
    return True


def admissible(edges, r0, c0, r, c, m):
    """G = K_{r,c} \\ Z must have no isolated vertex and max degree <= m."""
    ldeg = [0] * r
    rdeg = [0] * c
    for i, j in edges:
        ldeg[i] += 1
        rdeg[j] += 1
    for i in range(r):
        g = c - ldeg[i]
        if g < 1 or g > m:
            return False
    for j in range(c):
        g = r - rdeg[j]
        if g < 1 or g > m:
            return False
    return True


def orbit_census(n, m, a, b):
    """[sgn : C[X_{n,m,a}] tensor C[X_{n,m,b}]] via the paper's orbit method.

    Returns (multiplicity, list of contributing complement classes) or None
    when the complement would need more edges than the census enumerates.
    """
    r, c = n - a, n - b
    if r < 1 or c < 1 or r > n or c > n:
        return (0, [])
    e = r * c - n
    if e < 0:
        return (0, [])
    if e > MAX_COMPLEMENT_EDGES:
        raise ValueError("complement too large to enumerate: %d edges" % e)
    if e == 0:
        good = contributes((), 0, 0, r, c) and admissible((), 0, 0, r, c, m)
        return (1 if good else 0, [((), 0, 0)] if good else [])
    hits = []
    for edges, r0, c0 in core_classes(e):
        if contributes(edges, r0, c0, r, c) and admissible(edges, r0, c0, r, c, m):
            hits.append((edges, r0, c0))
    return (len(hits), hits)


def cell_parity(r, c, sigma, pi):
    """Sign of the permutation of all r*c cells of K_{r,c} induced by (sigma, pi)."""
    cells = [(i, j) for i in range(r) for j in range(c)]
    index = {cell: t for t, cell in enumerate(cells)}
    return parity(tuple(index[(sigma[i], pi[j])] for (i, j) in cells))


def _transposition(size, i):
    p = list(range(size))
    p[i], p[i + 1] = p[i + 1], p[i]
    return tuple(p)


def sigma_pi_probes(r, c, budget=20000):
    """(sigma, pi) probes for K_{r,c}: all of S_r x S_c when that is cheap.

    Otherwise the adjacent transpositions of each side -- which GENERATE
    S_r x S_c -- plus a deterministic spread of rotations.  Both sides of the
    identity below are homomorphisms S_r x S_c -> {+1,-1}, so agreeing on a
    generating set is agreement everywhere; the rotations are redundant
    confirmation rather than a sample.
    """
    if factorial(r) * factorial(c) <= budget:
        return [(s, p) for s in permutations(range(r))
                for p in permutations(range(c))], True
    idr, idc = tuple(range(r)), tuple(range(c))
    out = [(_transposition(r, i), idc) for i in range(r - 1)]
    out += [(idr, _transposition(c, j)) for j in range(c - 1)]
    for s in (1, 2, 3, 5, 7):
        for t in (1, 2, 3, 5, 7):
            out.append((tuple((i + s) % r for i in range(r)),
                        tuple((j + t) % c for j in range(c))))
    return out, False


def check_grid_sign_identity():
    """Equation (3) of the argument, verified exhaustively on small grids.

    For every subgraph G of K_{r,c} and every side-preserving automorphism
    (sigma, pi) of G, the induced permutations of E(G) and of E(Z) with
    Z = K_{r,c} \\ G must satisfy
        sgn_{E(G)} * sgn_{E(Z)} = sgn(sigma)^c * sgn(pi)^r.

    The exhaustive sweep below reaches only 4x3, whereas `contributes`
    evaluates the right-hand side at the grids the theorem really uses --
    (k+1, k-1) and (k, k) for odd k up to 61, i.e. up to 62x60.  So the
    underlying cell-action law
        sgn(action of (sigma, pi) on all r*c cells) = sgn(sigma)^c sgn(pi)^r
    is verified at those grids as well; equation (3) then follows because
    E(G) and E(Z) are complementary (sigma, pi)-invariant subsets of the
    cells and the sign of a permutation is multiplicative over invariant
    blocks -- which is itself what the small-grid sweep confirms directly.
    """
    tested, bad = 0, 0
    for r, c in ((2, 2), (3, 2), (2, 3), (3, 3), (4, 2), (2, 4), (4, 3), (3, 4)):
        cells = [(i, j) for i in range(r) for j in range(c)]
        for size in range(len(cells) + 1):
            for g in combinations(cells, size):
                gset = frozenset(g)
                z = sorted(set(cells) - gset)
                gs = sorted(g)
                ig = {cell: t for t, cell in enumerate(gs)}
                iz = {cell: t for t, cell in enumerate(z)}
                for sigma in permutations(range(r)):
                    for pi in permutations(range(c)):
                        if frozenset((sigma[i], pi[j]) for (i, j) in g) != gset:
                            continue
                        lhs = (edge_sign(gs, ig, sigma, pi)
                               * edge_sign(z, iz, sigma, pi))
                        rhs = (parity(sigma) ** c) * (parity(pi) ** r)
                        tested += 1
                        if lhs != rhs:
                            bad += 1
    # The same law at the grids the theorem actually uses, k = 5, 7, 13, 61.
    real_tested, real_bad, exhaustive = 0, [], []
    for r, c in ((6, 4), (5, 5), (8, 6), (7, 7), (14, 12), (13, 13),
                 (62, 60), (61, 61)):
        probes, full = sigma_pi_probes(r, c)
        exhaustive.append(((r, c), full))
        for sigma, pi in probes:
            real_tested += 1
            if cell_parity(r, c, sigma, pi) != (parity(sigma) ** c) * (parity(pi) ** r):
                real_bad.append((r, c, sigma, pi))
    ck("grid-sign-identity", bad == 0 and tested > 1000
       and not real_bad and real_tested > 1000,
       "%d (subgraph, automorphism) pairs on grids up to 4x3, %d violations; "
       "the cell-action law also holds on %d (sigma,pi) probes at the real "
       "grids 6x4, 5x5, 8x6, 7x7, 14x12, 13x13, 62x60, 61x61 (%d of the 8 "
       "swept over all of S_r x S_c, the rest over a generating set), "
       "%d violations"
       % (tested, bad, real_tested,
          sum(1 for _, f in exhaustive if f), len(real_bad)))


def canon(edges):
    r0, c0 = used_sides(edges)
    relab_r = {v: t for t, v in enumerate(sorted(set(i for i, _ in edges)))}
    relab_c = {v: t for t, v in enumerate(sorted(set(j for _, j in edges)))}
    norm = tuple(sorted((relab_r[i], relab_c[j]) for (i, j) in edges))
    return min(orbit_of(norm, r0, c0)) if edges else ()


def check_census_source():
    """Independent census of the source multiplicity, and its identification."""
    k = 5
    n, d = k * k - 5, (k * k - 5) - k
    mult, hits = orbit_census(n, k + 1, d - 1, d + 1)
    ck("source-multiplicity-k5-is-4", mult == PAPER_SOURCE_MULT,
       "n=%d, |P|=%d, |Q|=%d: census gives %d contributing orbits"
       % (n, n - (d - 1), n - (d + 1), mult))

    found = set(canon(e) for e, _, _ in hits)
    claimed = set(canon(e) for e in SOURCE_COMPLEMENTS.values())
    ck("source-orbits-are-the-exhibited-graphs", found == claimed,
       "census classes coincide with C_4, P_5^(3,2), P_5^(2,3), P_4+K_2"
       if found == claimed else "MISMATCH")

    bad = [k2 for k2 in (7, 9, 11, 13)
           if orbit_census(k2 * k2 - 5, k2 + 1,
                           (k2 * k2 - 5) - k2 - 1,
                           (k2 * k2 - 5) - k2 + 1)[0] != PAPER_SOURCE_MULT]
    ck("source-multiplicity-is-4-for-larger-k", not bad,
       "k = 7, 9, 11, 13 all give 4" if not bad else "failed at k in " + str(bad))


def check_census_target():
    """Independent census of the target multiplicity: 1 at k=5, 0 for k>=7."""
    k = 5
    n, d = k * k - 5, (k * k - 5) - k
    mult, hits = orbit_census(n, k + 1, d, d)
    ck("target-multiplicity-k5-is-1", mult == PAPER_TARGET_MULT[5],
       "n=%d, |P|=|Q|=%d: census gives %d contributing orbit(s)"
       % (n, n - d, mult))

    found = set(canon(e) for e, _, _ in hits)
    claimed = set(canon(e) for e in TARGET_COMPLEMENT.values())
    ck("target-orbit-is-the-exhibited-graph", found == claimed,
       "the single class is K_{1,2}+K_{2,1}+K_2"
       if found == claimed else "MISMATCH")

    bad = []
    for k2 in (7, 9, 11, 13):
        n2, d2 = k2 * k2 - 5, (k2 * k2 - 5) - k2
        if orbit_census(n2, k2 + 1, d2, d2)[0] != PAPER_TARGET_MULT_LARGE:
            bad.append(k2)
    ck("target-multiplicity-is-0-for-k-at-least-7", not bad,
       "k = 7, 9, 11, 13 all give 0" if not bad else "failed at k in " + str(bad))


def check_m_range_and_conclusion():
    """The two multiplicities are constant over k+1 <= m <= n, and 4 > target."""
    k = 5
    n, d = 20, 15
    rows = [(m, orbit_census(n, m, d - 1, d + 1)[0], orbit_census(n, m, d, d)[0])
            for m in range(k + 1, n + 1)]
    ck("multiplicities-constant-on-claimed-m-range",
       all(s == 4 and t == 1 for _, s, t in rows) and len(rows) == 15,
       "every m from 6 to 20 gives source 4, target 1")

    strict = []
    for k2 in (5, 7, 9, 11, 13):
        n2, d2 = k2 * k2 - 5, (k2 * k2 - 5) - k2
        s = orbit_census(n2, k2 + 1, d2 - 1, d2 + 1)[0]
        t = orbit_census(n2, k2 + 1, d2, d2)[0]
        strict.append((k2, s, t, s > t))
    ck("sign-multiplicity-strictly-drops", all(x[3] for x in strict),
       "; ".join("k=%d: %d > %d" % (a, b, c) for a, b, c, _ in strict))
    ck("m-hypothesis-is-used",
       orbit_census(n, k, d - 1, d + 1)[0] == 0
       and orbit_census(n, k + 1, d - 1, d + 1)[0] == 4,
       "at m = k = 5 the source multiplicity is 0, so the hypothesis "
       "m >= k+1 is not vacuous")
    ck("source-nonzero-so-degree-is-interior",
       all(x[1] > 0 for x in strict),
       "both neighbouring graded pieces are nonzero")


def check_lemma_holds_without_admissibility():
    """The paper's Lemma 3 is a statement about the sign condition alone.

    `orbit_census` intersects the sign condition with admissibility (G has no
    isolated vertex and max degree <= m).  If admissibility were doing the
    discriminating work, the census could return the paper's counts while
    Lemma 3 as stated was false.  Here the sign condition is evaluated on its
    own, so the identification is not carried by the extra filter.
    """
    src = set(canon(e) for e, r0, c0 in core_classes(4)
              if contributes(e, r0, c0, 6, 4))
    tgt = set(canon(e) for e, r0, c0 in core_classes(5)
              if contributes(e, r0, c0, 5, 5))
    ck("lemma-holds-for-the-sign-condition-alone",
       src == set(canon(e) for e in SOURCE_COMPLEMENTS.values())
       and tgt == set(canon(e) for e in TARGET_COMPLEMENT.values()),
       "the sign condition by itself isolates %d source and %d target "
       "class(es); admissibility removes nothing here" % (len(src), len(tgt)))


def check_enumeration_completeness():
    """core_classes(e) is complete, against a from-scratch orbit enumeration.

    Every multiplicity reported here is a count of classes returned by
    core_classes, so a class MISSING from that enumeration would lower a
    multiplicity quietly rather than fail anything -- 4 could read as 3 and
    the target's 1 could read as 0, which is the direction that would flatter
    the paper.  Re-derive the classes by sweeping orbits inside a single fixed
    ambient K_{e,e} (a different code path from core_classes' (r0, c0) sweep;
    legitimate because an e-edge graph with no isolated vertex uses at most e
    vertices per side, which check_target_zero... verifies by brute force) and
    compare canonical forms and counts.
    """
    report = []
    ok = True
    for e in (4, 5):
        cells = [(i, j) for i in range(e) for j in range(e)]
        seen, scratch = set(), set()
        for sub in combinations(cells, e):
            if sub in seen:
                continue
            seen |= orbit_of(sub, e, e)
            scratch.add(canon(sub))
        want = set(canon(x) for x, _, _ in core_classes(e))
        report.append((e, len(scratch), len(want)))
        if scratch != want:
            ok = False
    # Counts pinned as well as compared, so that both code paths shrinking
    # together would still fail.
    ck("core-class-enumeration-is-complete",
       ok and report == [(4, 16, 16), (5, 34, 34)],
       "; ".join("e=%d: %d classes from a K_{%d,%d} orbit sweep, %d from the "
                 "(r0,c0) sweep" % (e, a, e, e, b) for e, a, b in report))


def check_source_four_for_every_odd_k():
    """Source multiplicity 4 for every odd k >= 5, not just a sample.

    Two ingredients.  (i) For the source grid r = k+1, c = k-1 are both even,
    so the right-hand side sgn(sigma)^c sgn(pi)^r of the grid-sign identity is
    identically 1 and the isolated-vertex vetoes (which need an odd side) can
    never fire; the surviving set of 4-edge classes is therefore the same in
    every even-by-even ambient grid, which is checked directly.  (ii) The full
    census, admissibility included, is then run over a wide range of odd k.
    """
    classes = core_classes(4)
    base = set(canon(e) for e in SOURCE_COMPLEMENTS.values())
    ambient = []
    for r in range(4, 19, 2):
        for c in range(4, 19, 2):
            got = set(canon(e) for e, r0, c0 in classes
                      if contributes(e, r0, c0, r, c))
            ambient.append(((r, c), got == base))
    wide = []
    for k2 in range(5, 62, 2):
        n2 = k2 * k2 - 5
        d2 = n2 - k2
        wide.append((k2, orbit_census(n2, k2 + 1, d2 - 1, d2 + 1)[0]))
    ck("source-multiplicity-four-for-every-odd-k-to-61",
       all(ok for _, ok in ambient)
       and all(v == PAPER_SOURCE_MULT for _, v in wide),
       "sign condition picks the same 4 classes in all %d even ambient grids "
       "up to 18x18 (so the answer cannot depend on k), and the full census "
       "gives 4 for every odd k from 5 to 61" % len(ambient))


def check_target_zero_for_every_odd_k_at_least_7():
    """Target multiplicity 0 for EVERY odd k >= 7, from a finite certificate.

    The complement Z has exactly 5 edges, so its core uses at most 5 vertices
    on each side.  For odd k >= 7 the target grid is K_{k,k}, hence Z has at
    least k - 5 >= 2 isolated left vertices while c = k is odd: transposing
    two of them fixes E(Z) but has sgn(sigma)^c = -1, so the orbit is vetoed.
    Nothing else about k enters, so verifying max r0 <= 5 settles all odd
    k >= 7 rather than a sample of them.
    """
    classes = core_classes(5)
    # The bound "at most five used vertices per side" must NOT be read off
    # core_classes' own loop range (r0 < e + 1); doing so would make the
    # certificate a restatement of the enumerator's construction instead of a
    # verification of the paper's premise.  So recompute the used sides from
    # the edge lists, and verify the underlying fact by brute force inside an
    # ambient K_{6,6}, which is large enough for a six-used-vertex violation
    # to be representable if one existed.
    labels_agree = all(used_sides(e) == (r0, c0) for e, r0, c0 in classes)
    max_r0 = max(used_sides(e)[0] for e, _, _ in classes)
    max_c0 = max(used_sides(e)[1] for e, _, _ in classes)
    ambient = [(i, j) for i in range(6) for j in range(6)]
    brute_max = max(max(used_sides(s)) for s in combinations(ambient, 5))
    certificate = (labels_agree and max_r0 <= 5 and max_c0 <= 5
                   and brute_max <= 5 and (7 - max_r0) >= 2)
    survivors = [(k2, e) for k2 in range(7, 62, 2)
                 for e, r0, c0 in classes
                 if contributes(e, r0, c0, k2, k2)]
    ck("target-multiplicity-zero-for-every-odd-k-at-least-7",
       certificate and not survivors,
       "all %d five-edge classes use at most %d left vertices (and every "
       "5-subset of K_{6,6} uses at most %d per side, so the bound is a fact "
       "about five-edge graphs, not an artefact of the enumerator), so every "
       "odd k >= 7 leaves >= 2 isolated left vertices against an odd c = k; "
       "0 of %d classes survive at k = 7,9,...,61"
       % (len(classes), max_r0, brute_max, len(classes)))


def check_m_independence_over_wide_k():
    """The two multiplicities do not move as m ranges over k+1 <= m <= n."""
    rows = []
    for k2 in (5, 7, 9, 11, 13):
        n2 = k2 * k2 - 5
        d2 = n2 - k2
        src = set(orbit_census(n2, m, d2 - 1, d2 + 1)[0]
                  for m in range(k2 + 1, n2 + 1))
        tgt = set(orbit_census(n2, m, d2, d2)[0]
                  for m in range(k2 + 1, n2 + 1))
        rows.append((k2, src, tgt))
    want_tgt = {5: {PAPER_TARGET_MULT[5]}}
    ok = all(src == {PAPER_SOURCE_MULT}
             and tgt == want_tgt.get(k2, {PAPER_TARGET_MULT_LARGE})
             for k2, src, tgt in rows)
    ck("multiplicities-independent-of-m-over-wide-k", ok,
       "; ".join("k=%d: source %s, target %s over m in [%d,%d]"
                 % (k2, sorted(s), sorted(t), k2 + 1, k2 * k2 - 5)
                 for k2, s, t in rows))


# ----------------------------------------------------------------------
# Independent computation of the same multiplicity for small n, straight
# from the definition, by exact character theory over S_n.
# ----------------------------------------------------------------------

def factorial(n):
    f = 1
    for i in range(2, n + 1):
        f *= i
    return f


def set_partitions(n):
    """All set partitions of {0,...,n-1} as tuples of sorted tuples."""
    if n == 0:
        return [()]
    out = []
    for smaller in set_partitions(n - 1):
        for t in range(len(smaller)):
            out.append(smaller[:t] + (smaller[t] + (n - 1,),) + smaller[t + 1:])
        out.append(smaller + ((n - 1,),))
    return out


def integer_partitions(n, cap=None):
    if cap is None:
        cap = n
    if n == 0:
        return [()]
    out = []
    for first in range(min(n, cap), 0, -1):
        for rest in integer_partitions(n - first, first):
            out.append((first,) + rest)
    return out


def cycle_rep(lam, n):
    """A permutation of {0,...,n-1} with cycle type lam, as a tuple."""
    g = [0] * n
    at = 0
    for length in lam:
        for t in range(length):
            g[at + t] = at + (t + 1) % length
        at += length
    return tuple(g)


def class_size(lam, n):
    counts = {}
    for j in lam:
        counts[j] = counts.get(j, 0) + 1
    denom = 1
    for j, mult in counts.items():
        denom *= (j ** mult) * factorial(mult)
    return factorial(n) // denom


def fixed_by(part, g):
    blocks = set(part)
    for block in part:
        if tuple(sorted(g[x] for x in block)) not in blocks:
            return False
    return True


def direct_sign_mult(n, m, a, b, all_parts):
    """[sgn : C[X_{n,m,a}] tensor C[X_{n,m,b}]] by Burnside/character sum."""
    ra, rb = n - a, n - b
    def family(r):
        return [P for P in all_parts
                if len(P) == r and max(len(B) for B in P) <= m] if 1 <= r <= n else []
    Xa, Xb = family(ra), family(rb)
    total = 0
    for lam in integer_partitions(n):
        g = cycle_rep(lam, n)
        fa = sum(1 for P in Xa if fixed_by(P, g))
        if fa == 0:
            continue
        fb = fa if (ra == rb) else sum(1 for P in Xb if fixed_by(P, g))
        if fb == 0:
            continue
        sign = -1 if (n - len(lam)) % 2 else 1
        total += class_size(lam, n) * sign * fa * fb
    if total % factorial(n) != 0:
        raise ArithmeticError("non-integral multiplicity")
    return total // factorial(n)


def check_against_character_theory():
    """Cross-validate the orbit census against the definition, small n."""
    agree, disagree, positives = 0, [], 0
    for n in range(4, 9):
        all_parts = set_partitions(n)
        for m in range(2, n + 1):
            for a in range(0, n):
                for b in range(a, n):
                    r, c = n - a, n - b
                    e = r * c - n
                    if e < 0 or e > 5:
                        continue
                    got = orbit_census(n, m, a, b)[0]
                    want = direct_sign_mult(n, m, a, b, all_parts)
                    if got == want:
                        agree += 1
                        if want > 0:
                            positives += 1
                    else:
                        disagree.append((n, m, a, b, got, want))
    # `positives` is asserted nonzero on purpose: agreement on instances that
    # are all zero would be agreement about nothing.
    ck("orbit-census-matches-character-theory",
       not disagree and agree > 50 and positives >= 20,
       "%d instances with 4<=n<=8 agree (%d of them with nonzero "
       "multiplicity)%s" % (agree, positives,
                            "" if not disagree else "; MISMATCHES " + str(disagree[:5])))
    return agree


def check_named_small_instances():
    """A few named small instances, both methods, including odd-sided ones."""
    cases = [(6, 6, 3, 3), (7, 7, 4, 3), (8, 8, 4, 5),
             (6, 6, 4, 2), (7, 3, 4, 3), (8, 8, 5, 5)]
    rows = []
    for (n, m, a, b) in cases:
        all_parts = set_partitions(n)
        rows.append(((n, m, a, b), orbit_census(n, m, a, b)[0],
                     direct_sign_mult(n, m, a, b, all_parts)))
    ck("named-small-instances-agree", all(x[1] == x[2] for x in rows),
       "; ".join("n=%d,m=%d,(a,b)=(%d,%d): %d" % (c[0], c[1], c[2], c[3], g)
                 for c, g, _ in rows))
    ck("small-instances-are-informative",
       sum(1 for _, g, _ in rows if g > 0) >= 2,
       "at least two of the named instances have nonzero sign multiplicity")


def main():
    check_exhibited_objects()
    check_family_arithmetic([5, 7, 9, 11, 13, 15, 17, 19, 21])
    check_first_member_dimensions()
    check_grid_sign_identity()
    check_census_source()
    check_census_target()
    check_m_range_and_conclusion()
    check_lemma_holds_without_admissibility()
    check_enumeration_completeness()
    check_source_four_for_every_odd_k()
    check_target_zero_for_every_odd_k_at_least_7()
    check_m_independence_over_wide_k()
    check_named_small_instances()
    check_against_character_theory()
    # Guards against a check being skipped rather than failing: a raised
    # exception or an early return would otherwise shrink the total silently.
    ck("check-inventory-complete", len(CHECKS) == EXPECTED_CHECKS - 1,
       "%d of the expected %d checks ran before this one"
       % (len(CHECKS), EXPECTED_CHECKS - 1))
    print("NOTE not re-run here: the isomorphism R(Pi_{n,m})_d = C[X_{n,m,d}] "
          "(Zhu, Prop. 3.29) is taken as given, so the orbit-harmonics "
          "quotient is never built; the multiplicities at n=20 and beyond come "
          "from the complement census, which is cross-validated against "
          "brute-force character theory only for 4 <= n <= 8.  The dependence "
          "on k is settled for all odd k (source: the surviving 4-edge classes "
          "are the same in every even ambient grid; target k >= 7: every "
          "5-edge core uses at most 5 vertices per side, so the isolated-vertex "
          "veto always fires), with censuses run explicitly up to k = 61.")
    print("NOTE scope: this program does not re-derive the paper's full claim.  "
          "What it re-proves by exhaustion is Lemma 3 of the paper (the "
          "four-edge / five-edge classification) together with the sign "
          "multiplicities and the family arithmetic that follow from it; the "
          "paper has four further load-bearing steps this program does not "
          "verify on its own terms.  (a) Lemma 2, the orbit criterion via Frobenius "
          "reciprocity, is assumed by the census rather than proved.  (b) So is "
          "the vanishing of NON-TRANSVERSE orbits (a block pair with "
          "|B and C| >= 2 admits an odd transposition in the stabiliser): the "
          "census enumerates transverse pairs only.  Steps (a) and (b) are the "
          "exact content of the census-vs-character-theory agreement, since "
          "direct_sign_mult counts every orbit straight from the definition -- "
          "so they are corroborated, but only for 4 <= n <= 8, not at n = 20.  "
          "(c) The closing inference -- complex S_n-representations are "
          "semisimple, so an injection would force source multiplicities to be "
          "dominated by target ones -- is pure representation theory and is "
          "not checked numerically; the checks establish only the strict drop "
          "4 > 1 and 4 > 0.  (d) Over the m-range k+1 <= m <= n, every m is "
          "swept for k = 5,...,13, but for k up to 61 only m = k+1 is run; "
          "that is the binding case because the admissibility cap 'every "
          "degree of G is at most m' is monotone in m, so larger m can only "
          "admit more.  No claim is made for even k or for k = 3, neither of "
          "which the paper asserts.")
    n_fail = sum(1 for _, ok in CHECKS if not ok)
    n = len(CHECKS)
    if n_fail == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
        return 0
    print("VERDICT: %d OF %d CHECKS FAILED" % (n_fail, n))
    return 1


if __name__ == "__main__":
    sys.exit(main())

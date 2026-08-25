#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- independent verification of

    "A Negative Answer to Mohammadi's Question on Powers of Cover Ideals"

    G is the graph on [8] with E(G) = {14,16,17,23,26,28,37,48,56,57,58}.
    Claim: G is Cohen-Macaulay over every field, J(G) has a 5-linear
    resolution, and beta_{1,a}(J(G)^2) = 1 for a = (1,1,2,2,2,2,1,1);
    hence J(G)^2 is generated in degree 10 but has a first syzygy in
    degree 12 and so has no linear resolution over any field.

Python 3.9, standard library only.  All arithmetic is exact: integers,
fractions.Fraction over Q, and integers mod p over F_p.  There is no
floating point anywhere in this file and no check is numeric.

--------------------------------------------------------------------------
TAKEN FROM THE PAPER (data, not verified here)
--------------------------------------------------------------------------
  * the edge list E(G) = {14,16,17,23,26,28,37,48,56,57,58};
  * display (1): the ten facets of Ind(G) in the paper's shelling order
    345, 135, 245, 125, 346, 247, 467, 138, 368, 678;
  * the paper's split of the maximal independent sets by intersection
    with {3,4} (used only as a cross-check of display (1));
  * the printed restriction lists
    35; 45; 15,25; 34; 24; 46,47; 13; 36,38; 67,68;
  * the printed minimal vertex covers C_1,...,C_10;
  * the multidegree a = (1,1,2,2,2,2,1,1) and m = x^a;
  * table (3): for each row, the singleton F, the cross-edge F, the list of
    eligible covers and the common vertex t;
  * the asserted numerals -- those the paper prints, together with the
    further comparison targets fixed here -- used ONLY as comparison
    targets on the right-hand side of a check:
      10 facets, each of size 3;  10 minimal vertex covers, each of size 5;
      53 minimal generators of J(G)^2, each of degree 10;
      K^a = <34,56>, two connected components, dim H~_0 = 1;
      beta_{1,a}(J(G)^2) = 1;
      Betti table of J(G): beta_{0,5}=10, beta_{1,6}=13, beta_{2,7}=4;
      total beta_1(J(G)^2) = 122 in degree 11 and exactly 1 in degree 12.

--------------------------------------------------------------------------
DERIVED HERE (this is what the checks actually decide)
--------------------------------------------------------------------------
  * all maximal independent sets of G, by brute force over the 256 subsets
    of [8]; their number, their sizes, and equality as a SET with display (1);
  * all minimal vertex covers of G, again by brute force over the 256
    subsets -- derived independently of complementation -- their number and
    sizes and equality as a set with the printed C_1..C_10; plus the test
    "S is a vertex cover  <=>  some C_i is contained in S" for all 256 S,
    which certifies that the ten printed monomials generate exactly the
    squarefree monomials of J(G) = intersection over edges of (x_i,x_j);
  * for each j = 2..10, the maximal elements of {F_i cap F_j : i < j};
    each list is compared with the paper's, and checked nonempty and pure of
    dimension 1, which is the shelling condition (hence Cohen-Macaulay over
    every field);
  * the 55 products x_{C_i}x_{C_j} (1<=i<=j<=10): the number of distinct
    monomials, whether each is minimal under divisibility, and their degrees;
  * the two identities of display (2), by multiplying exponent vectors;
  * the upper Koszul complex K^a = {F subset [8] : m/x_F in J(G)^2}, by
    testing all 256 F against divisibility by the 53 derived generators;
    its downward closedness, its face list, its vertex set, its number of
    connected components, and dim H~_0 over Q and over F_2, F_3, F_5, F_7,
    F_101 from boundary-matrix ranks;
  * table (3) row by row: the eligible-cover list re-derived from the
    singleton condition AND from the cross-edge condition (the two must
    agree with each other and with the printed list), the exponent-2
    coordinates of m/x_F, and the fact that every eligible cover contains t
    while m/x_F has x_t-exponent 1;
  * beta_{1,b}(J(G)^2) for b = a and for every pairwise lcm multidegree, by
    minimal-first-syzygy linear algebra
        beta_{1,b} = dim ker((S^53)_b -> S_b) - rank(sum_t x_t ker_{b-e_t}),
    with Gaussian elimination over Q and over F_2, F_3, F_5, F_7, F_101.
    This route is independent of \\cite[Thm 1.34]{MillerSturmfels}, the one
    citation that cannot be checked from inside this file.  The Taylor
    complex confines the support of beta_1 to pairwise lcms, so the scan
    is exhaustive; it yields
    the total beta_1 split by total degree and the set of degree-12
    multidegrees carrying a first syzygy;
  * the whole multigraded Betti table of J(G), via upper Koszul complexes
    over all 256 squarefree multidegrees and reduced simplicial homology, so
    5-linearity is decided here and not merely inherited from Eagon-Reiner;
  * the remark checks: the induced subgraph on {1,6,5,7} is a 4-cycle (so G
    is not chordal), x_1...x_8 lies in the symbolic square
    intersection over edges of (x_i,x_j)^2 but not in J(G)^2, and the ten
    generators of J(G) have gcd 1 (the Remark's separation from
    Ficarra-Moradi);
  * beta_{1,b}(J(G)^2) a THIRD time, with no field at all
    (beta1_combinatorial): the spanning set of (mZ)_b consists only of
    differences e_i - e_j inside a common index set, so its rank is |U| - c
    over every field.  This is what turns "over the six fields tested" into
    "over every field", which is what the paper asserts;
  * a cross-validation of the upper-Koszul engine against the
    Miller-Sturmfels-free syzygy engine on J(G) itself (section 7b): both
    must give the same beta_1 by degree.

--------------------------------------------------------------------------
HONEST ACCOUNTING OF WHAT IS ASSUMED
--------------------------------------------------------------------------
  * \\cite[Theorem~1.34]{MillerSturmfels} (beta_{i,b} = dim H~_{i-1}(K^b)) is
    the DEFINITION of the upper-Koszul engine used in sections 4 and 7.  It is
    not assumed for the headline -- section 6 decides beta_{1,a}(J(G)^2) by
    syzygy linear algebra and field-free combinatorics, and section 7b shows
    the two engines agree on beta_1(J(G)) -- but beta_{2,7}(J(G)) = 4 does
    rest on it.
  * "shellable => Cohen-Macaulay" and Eagon-Reiner are cited theorems.  What
    is verified here is the shelling itself (field-free) and, separately, the
    Betti table.
  * Checks whose label carries "[STRUCTURAL: ...]" cannot fail for any input;
    they document a step whose load-bearing content is decided elsewhere, and
    they should not be counted as evidence.
  * The wording of Mohammadi's Question 4.1 and of the Ha-Van Tuyl
    restatement are textual and outside the reach of any program; main()
    prints this and the other gaps before the verdict.
"""

import sys
from fractions import Fraction

N = 8
VERTS = tuple(range(1, N + 1))

# ---------------------------------------------------------------------------
# DATA FROM THE PAPER
# ---------------------------------------------------------------------------

# E(G), Theorem 1.
PAPER_EDGES = [(1, 4), (1, 6), (1, 7), (2, 3), (2, 6), (2, 8),
               (3, 7), (4, 8), (5, 6), (5, 7), (5, 8)]

# Display (1): the ten facets of Ind(G), in the paper's shelling order.
PAPER_FACETS = ["345", "135", "245", "125", "346",
                "247", "467", "138", "368", "678"]

# The paper's split of the maximal independent sets by intersection with {3,4}.
PAPER_SPLIT = {
    "both 3,4":        ["345", "346"],
    "3 but not 4":     ["135", "138", "368"],
    "4 but not 3":     ["245", "247", "467"],
    "neither 3 nor 4": ["125", "678"],
}

# The printed restriction lists, for j = 2,...,10.
PAPER_RESTRICTIONS = [["35"], ["45"], ["15", "25"], ["34"], ["24"],
                      ["46", "47"], ["13"], ["36", "38"], ["67", "68"]]

# The printed minimal vertex covers C_1,...,C_10.
PAPER_COVERS = ["12678", "24678", "13678", "34678", "12578",
                "13568", "12358", "24567", "12457", "12345"]

# The multidegree a of the claimed first syzygy; m = x^a.
PAPER_A = (1, 1, 2, 2, 2, 2, 1, 1)

# Table (3): (singleton F, cross-edge F, 1-based eligible cover indices, t).
PAPER_TABLE3 = [
    ("1", "35", [2, 4, 8], 7),
    ("2", "45", [3, 4, 6], 8),
    ("7", "46", [6, 7, 10], 1),
    ("8", "36", [8, 9, 10], 2),
]

# K^a as claimed in the paper: generated by the two edges 34 and 56.
PAPER_KOSZUL_GENS = ["34", "56"]

# Asserted numerals: those the paper prints, plus the further comparison
# targets fixed here.
EXP_N_FACETS = 10
EXP_FACET_SIZE = 3
EXP_N_COVERS = 10
EXP_COVER_SIZE = 5
EXP_N_PRODUCTS = 55
EXP_N_GENERATORS = 53
EXP_GEN_DEGREE = 10
EXP_KOSZUL_FACES = ["", "3", "4", "5", "6", "34", "56"]
EXP_KOSZUL_COMPONENTS = 2
EXP_BETA1_A = 1
EXP_BETTI_J = {(0, 5): 10, (1, 6): 13, (2, 7): 4}
EXP_TOTAL_BETA1_BY_DEGREE = {11: 122, 12: 1}

# Fields over which every homological quantity is recomputed.
# 0 denotes Q (exact rational arithmetic with fractions.Fraction).
FIELDS = [0, 2, 3, 5, 7, 101]


def field_name(p):
    return "Q" if p == 0 else "F_%d" % p


# ---------------------------------------------------------------------------
# sets of vertices and monomials
# ---------------------------------------------------------------------------
# A vertex set is a frozenset of elements of 1..8; it is printed in the
# paper's concatenated shorthand, e.g. frozenset({3,4,5}) -> "345".
# A monomial is an 8-tuple of nonnegative integer exponents.

def parse_set(s):
    """'345' -> frozenset({3,4,5}); '' -> the empty set."""
    return frozenset(int(ch) for ch in s)


def show_set(S):
    """frozenset({3,4,5}) -> '345'; the empty set -> '{}'."""
    return "".join(str(v) for v in sorted(S)) if S else "{}"


def show_sets(collection):
    """A canonical, order-independent rendering of a family of vertex sets."""
    return ",".join(sorted(show_set(S) for S in collection))


def x_set(S):
    """The squarefree monomial x_S as an exponent tuple."""
    return tuple(1 if v in S else 0 for v in VERTS)


def mon_mul(u, v):
    return tuple(a + b for a, b in zip(u, v))


def mon_div(u, v):
    """u/v as an exponent tuple; raises if v does not divide u."""
    w = tuple(a - b for a, b in zip(u, v))
    if min(w) < 0:
        raise ValueError("not divisible")
    return w


def divides(u, v):
    """Does the monomial u divide the monomial v?"""
    return all(a <= b for a, b in zip(u, v))


def mon_lcm(u, v):
    return tuple(max(a, b) for a, b in zip(u, v))


def deg(u):
    return sum(u)


def support(u):
    return frozenset(v for v, e in zip(VERTS, u) if e > 0)


def show_mon(u):
    if deg(u) == 0:
        return "1"
    out = []
    for v, e in zip(VERTS, u):
        if e == 1:
            out.append("x%d" % v)
        elif e > 1:
            out.append("x%d^%d" % (v, e))
    return "*".join(out)


# ---------------------------------------------------------------------------
# exact linear algebra
# ---------------------------------------------------------------------------

def matrix_rank(rows, ncols, p):
    """Rank of the matrix with the given rows (lists of ints) over the field
    Q (p == 0, arithmetic in fractions.Fraction) or F_p (p prime).

    Plain Gaussian elimination.  No floating point: over Q every entry is a
    Fraction, over F_p every entry is an int in [0,p).
    """
    if ncols == 0:
        return 0
    if p == 0:
        work = [[Fraction(int(a)) for a in r] for r in rows]
        zero = Fraction(0)
    else:
        work = [[int(a) % p for a in r] for r in rows]
        zero = 0
    nrows = len(work)
    rank = 0
    row = 0
    for col in range(ncols):
        piv = None
        for r in range(row, nrows):
            if work[r][col] != zero:
                piv = r
                break
        if piv is None:
            continue
        work[row], work[piv] = work[piv], work[row]
        pv = work[row][col]
        if p == 0:
            inv = Fraction(1) / pv
            work[row] = [e * inv for e in work[row]]
        else:
            inv = pow(pv, p - 2, p)
            work[row] = [(e * inv) % p for e in work[row]]
        for r in range(row + 1, nrows):
            f = work[r][col]
            if f != zero:
                if p == 0:
                    work[r] = [e - f * g for e, g in zip(work[r], work[row])]
                else:
                    work[r] = [(e - f * g) % p
                               for e, g in zip(work[r], work[row])]
        rank += 1
        row += 1
        if row == nrows:
            break
    return rank


# ---------------------------------------------------------------------------
# reduced simplicial homology of a complex given as a set of faces
# ---------------------------------------------------------------------------

def reduced_homology(faces, p):
    """dim_K H~_k(Delta; K) for every k, from ranks of the boundary maps of
    the augmented chain complex.

    'faces' is a set of frozensets.  The empty face is supplied here if the
    complex is nonvoid, so H~_{-1} is 1 exactly for Delta = {emptyset}.
    Returns a dict k -> dim, for k from -1 to dim(Delta).
    """
    faces = set(faces)
    if not faces:
        return {}                      # the void complex: no chain groups
    faces.add(frozenset())
    bysize = {}
    for F in faces:
        bysize.setdefault(len(F), []).append(F)
    maxsize = max(bysize)
    index = {}
    for s in bysize:
        bysize[s] = sorted(bysize[s], key=lambda F: sorted(F))
        index[s] = {F: i for i, F in enumerate(bysize[s])}

    # rank of boundary_k : C_k -> C_{k-1}, i.e. size k+1 -> size k
    ranks = {}
    for s in range(1, maxsize + 1):
        src = bysize.get(s, [])
        tgt = bysize.get(s - 1, [])
        if not src or not tgt:
            ranks[s - 1] = 0
            continue
        rows = []
        for F in src:
            vs = sorted(F)
            row = [0] * len(tgt)
            for i, v in enumerate(vs):
                G = F - {v}
                row[index[s - 1][G]] = -1 if i % 2 else 1
            rows.append(row)
        ranks[s - 1] = matrix_rank(rows, len(tgt), p)

    out = {}
    for k in range(-1, maxsize):
        dim_ck = len(bysize.get(k + 1, []))
        out[k] = dim_ck - ranks.get(k, 0) - ranks.get(k + 1, 0)
    return out


# ---------------------------------------------------------------------------
# derivations straight from E(G)
# ---------------------------------------------------------------------------

_ALL_SUBSETS = None


def all_subsets():
    """The 2^8 = 256 subsets of [8], as frozensets.  Built once; the list is
    never mutated by any caller."""
    global _ALL_SUBSETS
    if _ALL_SUBSETS is None:
        _ALL_SUBSETS = [frozenset(v for v in VERTS if mask >> (v - 1) & 1)
                        for mask in range(1 << N)]
    return _ALL_SUBSETS


def is_independent(S, edges):
    return all(not (i in S and j in S) for i, j in edges)


def is_cover(S, edges):
    return all(i in S or j in S for i, j in edges)


def maximal_independent_sets(edges):
    """Every maximal independent set of the graph, by brute force over all
    256 subsets: independent, and not contained in a larger independent set.
    """
    indep = [S for S in all_subsets() if is_independent(S, edges)]
    indep_set = set(indep)
    out = []
    for S in indep:
        if all((S | {v}) not in indep_set for v in VERTS if v not in S):
            out.append(S)
    return set(out)


def minimal_vertex_covers(edges):
    """Every minimal vertex cover, derived from E(G) directly (not by
    complementing the independent sets), by brute force over all 256 subsets.
    """
    covers = [S for S in all_subsets() if is_cover(S, edges)]
    cover_set = set(covers)
    out = []
    for S in covers:
        if all((S - {v}) not in cover_set for v in S):
            out.append(S)
    return set(out)


def maximal_elements(family):
    """The elements of 'family' maximal under inclusion."""
    fam = set(family)
    return set(S for S in fam
               if not any(S < T for T in fam))


def shelling_restrictions(facets):
    """For j = 2..len(facets), the facets of  F_j cap <F_1,...,F_{j-1}>,
    i.e. the maximal elements of {F_i cap F_j : i < j}.

    'facets' is the ordered list from display (1).  The subcomplex generated
    by F_1..F_{j-1} intersected with the simplex on F_j is generated by the
    pairwise intersections, so its facets are exactly those maximal elements.
    """
    out = []
    for j in range(1, len(facets)):
        inter = [facets[i] & facets[j] for i in range(j)]
        out.append(maximal_elements(inter))
    return out


# ---------------------------------------------------------------------------
# J(G)^2 : products of the cover monomials, and their minimalisation
# ---------------------------------------------------------------------------

def square_products(covers):
    """The 55 monomials x_{C_i} x_{C_j}, 1 <= i <= j <= 10, in order,
    together with the index pair that produced each one."""
    out = []
    for i in range(len(covers)):
        for j in range(i, len(covers)):
            out.append(((i + 1, j + 1),
                        mon_mul(x_set(covers[i]), x_set(covers[j]))))
    return out


def minimalise(monomials):
    """The minimal elements of a set of monomials under divisibility: the
    minimal monomial generating set of the ideal they generate."""
    ms = sorted(set(monomials))
    return sorted(u for u in ms
                  if not any(v != u and divides(v, u) for v in ms))


def in_ideal(mon, gens):
    """Is the monomial 'mon' in the monomial ideal generated by 'gens'?"""
    return any(divides(g, mon) for g in gens)


def upper_koszul_complex(b, gens):
    """K^b(I) = { F subset supp(b) : x^b / x_F  in  I },  Miller-Sturmfels'
    upper Koszul simplicial complex.  Returned as a set of frozensets.

    Every F is tested by explicit division and explicit divisibility against
    the generators; nothing is assumed about which F can occur.
    """
    supp = support(b)
    out = set()
    for F in all_subsets():
        if not F <= supp:
            continue
        q = mon_div(b, x_set(F))
        if in_ideal(q, gens):
            out.add(F)
    return out


# ---------------------------------------------------------------------------
# beta_1 in a multidegree by minimal-first-syzygy linear algebra
# ---------------------------------------------------------------------------

def beta1_multidegree(b, gens, p, report=None):
    """dim_K Tor_1^S(I, K)_b for the ideal I minimally generated by 'gens'.

    With F = S^gens -> I the minimal presentation and Z = ker F, one has
    beta_{1,b}(I) = dim Z_b - dim (m Z)_b,  (m Z)_b = sum_t x_t Z_{b-e_t}.

    The K-basis of (S^gens)_b is { e_i x^{b-a_i} : a_i <= b }, and F sends
    every such basis vector to x^b.  Hence dim Z_b = n_b - 1 when n_b >= 1,
    where n_b = #{ i : a_i <= b }, and Z_b is the sum-zero hyperplane in the
    coordinates indexed by those i.  Multiplication by x_t carries the basis
    vector e_i x^{(b-e_t)-a_i} to e_i x^{b-a_i}, so it is the coordinate
    inclusion of the index set I_t = { i : a_i <= b - e_t }; therefore
    x_t Z_{b-e_t} is the sum-zero subspace supported on I_t, spanned by the
    differences e_{i_0} - e_{i_k}, k >= 1.

    Rank of the span of all of those, over Q or F_p, is computed by Gaussian
    elimination.  No Miller-Sturmfels input is used anywhere here.
    """
    idx = [i for i, g in enumerate(gens) if divides(g, b)]
    n = len(idx)
    if n == 0:
        if report is not None:
            report.update(n_b=0, dim_Z=0, rank_mZ=0)
        return 0
    pos = {i: k for k, i in enumerate(idx)}
    rows = []
    it_sizes = []
    for t in range(N):
        if b[t] == 0:
            it_sizes.append(0)
            continue
        bt = list(b)
        bt[t] -= 1
        bt = tuple(bt)
        sub = [i for i in idx if divides(gens[i], bt)]
        it_sizes.append(len(sub))
        for k in range(1, len(sub)):
            row = [0] * n
            row[pos[sub[0]]] = 1
            row[pos[sub[k]]] = -1
            rows.append(row)
    r = matrix_rank(rows, n, p)
    if report is not None:
        report.update(n_b=n, dim_Z=n - 1, rank_mZ=r, it_sizes=it_sizes,
                      gen_indices=[i + 1 for i in idx])
    return (n - 1) - r


def beta1_combinatorial(b, gens):
    """The SAME beta_{1,b} as beta1_multidegree, computed with no field at all.

    The spanning set of (mZ)_b built in beta1_multidegree consists only of
    differences e_i - e_j with i,j in one common index set
    I_t = { i : a_i <= b - e_t }.  The span of all differences inside a family
    of index sets is, over ANY field, the sum-zero subspace of each connected
    component of the graph whose edges join the members of a common I_t.  Its
    dimension is therefore |U| - c, where U = union of the I_t and c is the
    number of those components -- a purely combinatorial number, the same over
    every field, prime or characteristic zero.  Hence

        beta_{1,b} = (n_b - 1) - (|U| - c)

    is INDEPENDENT OF THE FIELD.  This is what upgrades "over the six fields
    tested" to "over every field", for beta_1 and hence for the paper's
    conclusion.  (The same remark applies to the paper's own route, since
    dim H~_0 = #components - 1 over any field.)
    """
    idx = [i for i, g in enumerate(gens) if divides(g, b)]
    if not idx:
        return 0
    parent = {i: i for i in idx}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    used = set()
    for t in range(N):
        if b[t] == 0:
            continue
        bt = list(b)
        bt[t] -= 1
        sub = [i for i in idx if divides(gens[i], tuple(bt))]
        used.update(sub)
        for i in sub[1:]:
            ra, rb = find(sub[0]), find(i)
            if ra != rb:
                parent[ra] = rb
    comps = set(find(i) for i in used)
    rank = len(used) - len(comps)
    return (len(idx) - 1) - rank


# ---------------------------------------------------------------------------
# combinatorics of a simplicial complex
# ---------------------------------------------------------------------------

def is_downward_closed(faces):
    """Is the given family of faces closed under passing to subsets?"""
    fam = set(faces)
    for F in fam:
        for v in F:
            if (F - {v}) not in fam:
                return False
    return True


def complex_vertices(faces):
    out = set()
    for F in faces:
        out |= set(F)
    return frozenset(out)


def connected_components(faces):
    """Connected components of the 1-skeleton, as a sorted list of vertex
    sets.  Isolated vertices count as their own component."""
    verts = sorted(complex_vertices(faces))
    parent = {v: v for v in verts}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for F in faces:
        vs = sorted(F)
        for w in vs[1:]:
            a, b = find(vs[0]), find(w)
            if a != b:
                parent[a] = b
    buckets = {}
    for v in verts:
        buckets.setdefault(find(v), set()).add(v)
    return sorted((frozenset(s) for s in buckets.values()),
                  key=lambda s: sorted(s))


# ---------------------------------------------------------------------------
# table (3) of the paper, re-derived
# ---------------------------------------------------------------------------

def derive_table3_row(row, a, covers, gens):
    """Re-derive one row of table (3).

    The row names a singleton F = {r} and a cross-edge F = {p,q}, a list of
    'eligible covers' and a vertex t.  Both eligibility conditions are
    re-derived here from the paper's own reasoning:

      * F = {r}: the x_r-exponent of m/x_F is a_r - 1 = 0, so each factor
        cover must omit r.  Eligible = covers omitting r.
      * F = {p,q}: deg(m/x_F) = |a| - 2 = 10 = deg(x_{C_i} x_{C_j}), so
        divisibility forces equality, and every coordinate of exponent 2 in
        m/x_F must lie in both covers.  Eligible = covers containing all such
        coordinates.

    Returns a dict of derived facts.
    """
    rstr, pqstr, printed_idx, t = row
    r = int(rstr)
    pq = parse_set(pqstr)
    m = tuple(a)

    q1 = mon_div(m, x_set(frozenset([r])))
    elig1 = frozenset(i + 1 for i, C in enumerate(covers) if r not in C)

    q2 = mon_div(m, x_set(pq))
    exp2 = frozenset(v for v, e in zip(VERTS, q2) if e == 2)
    elig2 = frozenset(i + 1 for i, C in enumerate(covers) if exp2 <= C)

    return {
        "r": r, "pq": pq, "t": t,
        "printed": frozenset(printed_idx),
        "q_singleton": q1, "exp_r_in_q1": q1[r - 1],
        "q_pair": q2, "deg_q_pair": deg(q2), "exp2_coords": exp2,
        "elig_singleton": elig1, "elig_pair": elig2,
        "t_in_all": all(t in covers[i - 1] for i in sorted(elig1 | elig2)),
        "exp_t_q1": q1[t - 1], "exp_t_q2": q2[t - 1],
        "q1_in_square": in_ideal(q1, gens),
        "q2_in_square": in_ideal(q2, gens),
    }


# ---------------------------------------------------------------------------
# the full multigraded Betti table of a squarefree monomial ideal
# ---------------------------------------------------------------------------

_KOSZUL_CACHE = {}


def squarefree_koszul_complexes(gens):
    """[(b, K^b(I)) for all 256 squarefree multidegrees b].  Cached, because
    the complexes do not depend on the field; only their homology does."""
    if gens not in _KOSZUL_CACHE:
        _KOSZUL_CACHE[gens] = [(x_set(S), upper_koszul_complex(x_set(S), gens))
                               for S in all_subsets()]
    return _KOSZUL_CACHE[gens]


def betti_table_squarefree(gens, p):
    """beta_{i,j}(I) for a squarefree monomial ideal I, summed over all
    squarefree multidegrees b with |b| = j.

    beta_{i,b}(I) = dim H~_{i-1}(K^b(I); K).  The minimal free resolution of a
    squarefree monomial ideal is supported on the lcm lattice, whose elements
    are squarefree, so the 256 squarefree b exhaust the support -- nothing is
    missed by restricting to them.

    Returns dict (i, j) -> beta, omitting zeros.
    """
    out = {}
    for b, K in squarefree_koszul_complexes(tuple(gens)):
        h = reduced_homology(K, p)
        for k, d in h.items():
            if d:
                key = (k + 1, deg(b))
                out[key] = out.get(key, 0) + d
    return {k: v for k, v in out.items() if v}


def beta1_full_scan(gens, p):
    """beta_{1,b}(I) at every multidegree that can carry a first syzygy.

    The Taylor complex on the generators g_1,...,g_n is a (generally
    non-minimal) free resolution of I whose first free module has basis in the
    multidegrees lcm(g_i,g_j), i < j.  The minimal resolution is a multigraded
    direct summand in each multidegree, so beta_{1,b} = 0 unless b is such an
    lcm.  Scanning those lcms is therefore exhaustive.

    Returns (dict b -> beta_1 for the nonzero ones, number of lcms scanned).
    """
    lcms = set()
    for i in range(len(gens)):
        for j in range(i, len(gens)):
            lcms.add(mon_lcm(gens[i], gens[j]))
    nonzero = {}
    for b in sorted(lcms):
        v = beta1_multidegree(b, gens, p)
        if v:
            nonzero[b] = v
    return nonzero, len(lcms)


def beta1_full_scan_combinatorial(gens):
    """The same exhaustive pairwise-lcm scan as beta1_full_scan, but with the
    field-free route of beta1_combinatorial.  Returns (dict b -> beta_1 for the
    nonzero ones, number of lcms scanned).  Because no field enters, the
    resulting figures hold over EVERY field simultaneously."""
    lcms = set()
    for i in range(len(gens)):
        for j in range(i, len(gens)):
            lcms.add(mon_lcm(gens[i], gens[j]))
    nonzero = {}
    for b in sorted(lcms):
        v = beta1_combinatorial(b, gens)
        if v:
            nonzero[b] = v
    return nonzero, len(lcms)


# ---------------------------------------------------------------------------
# check harness
# ---------------------------------------------------------------------------

RESULTS = []

# Named derived facts, recorded by the sections and combined at the end.  The
# closing summary is a conjunction of THESE, not of "did every earlier check
# pass", so it decides something.
FACTS = {}


def check(ok, label):
    RESULTS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + label)
    return bool(ok)


def info(text):
    print("      " + text)


def head(text):
    print("")
    print("--- " + text)


# ---------------------------------------------------------------------------
# section 1: Ind(G), the minimal vertex covers, and display (1)
# ---------------------------------------------------------------------------

def section_graph():
    head("1. maximal independent sets and minimal vertex covers of G")
    edges = [tuple(sorted(e)) for e in PAPER_EDGES]
    info("E(G) taken from the paper: " +
         ",".join("%d%d" % e for e in edges) + "  (%d edges)" % len(edges))

    mis = maximal_independent_sets(edges)
    info("derived maximal independent sets (%d): %s"
         % (len(mis), show_sets(mis)))
    check(len(mis) == EXP_N_FACETS,
          "Ind(G) has %d facets (derived %d)" % (EXP_N_FACETS, len(mis)))
    sizes = sorted(set(len(S) for S in mis))
    check(sizes == [EXP_FACET_SIZE],
          "every facet of Ind(G) has size %d, so Ind(G) is pure of dimension "
          "2 (derived sizes %s)" % (EXP_FACET_SIZE, sizes))

    facets = [parse_set(s) for s in PAPER_FACETS]
    check(len(set(facets)) == len(facets) and set(facets) == mis,
          "display (1) is exactly the derived set of maximal independent sets")

    split = set()
    for grp in PAPER_SPLIT.values():
        split |= set(parse_set(s) for s in grp)
    check(split == mis,
          "the paper's split by intersection with {3,4} lists exactly the "
          "derived maximal independent sets")
    # Comparing only the UNION would leave the four LABELS of the split
    # unverified: a set listed in the wrong row, and the paper's count
    # 2+3+3+2=10, would both go unnoticed.  Test the classification itself.
    want = {"both 3,4": frozenset([3, 4]), "3 but not 4": frozenset([3]),
            "4 but not 3": frozenset([4]), "neither 3 nor 4": frozenset()}
    labels_ok = True
    for lab, grp in PAPER_SPLIT.items():
        for s in grp:
            if (parse_set(s) & frozenset([3, 4])) != want[lab]:
                labels_ok = False
                info("MISFILED: %s is listed under '%s' but meets {3,4} in %s"
                     % (s, lab, show_set(parse_set(s) & frozenset([3, 4]))))
    check(labels_ok and sum(len(g) for g in PAPER_SPLIT.values()) == len(mis),
          "each of the paper's four rows really contains the maximal "
          "independent sets with that intersection with {3,4}, and the row "
          "sizes 2+3+3+2 total %d" % len(mis))

    mvc = minimal_vertex_covers(edges)
    info("derived minimal vertex covers (%d): %s" % (len(mvc), show_sets(mvc)))
    check(len(mvc) == EXP_N_COVERS,
          "G has %d minimal vertex covers (derived %d)"
          % (EXP_N_COVERS, len(mvc)))
    csizes = sorted(set(len(S) for S in mvc))
    check(csizes == [EXP_COVER_SIZE],
          "every minimal vertex cover has size %d (derived sizes %s)"
          % (EXP_COVER_SIZE, csizes))

    covers = [parse_set(s) for s in PAPER_COVERS]
    check(len(set(covers)) == len(covers) and set(covers) == mvc,
          "the printed C_1..C_10 are exactly the derived minimal vertex covers")
    check(all(covers[i] == frozenset(VERTS) - facets[i]
              for i in range(len(facets))),
          "C_i is the complement of the i-th facet of display (1), for all i")

    bad = [S for S in all_subsets()
           if is_cover(S, edges) != any(C <= S for C in covers)]
    check(not bad,
          "for all 256 subsets S: x_S lies in the ideal generated by "
          "x_{C_1},...,x_{C_10}  iff  S is a vertex cover, so those ten "
          "monomials generate J(G) = intersection of (x_i,x_j) "
          "(%d discrepancies)  [STRUCTURAL: once the C_i are certified to be "
          "exactly the MINIMAL covers, 'every cover contains a minimal cover' "
          "and 'every superset of a cover is a cover' make this line "
          "unfailable; it records the bridge from covers to ideal generators, "
          "which itself needs the standard reduction that membership in a "
          "monomial ideal generated by squarefree monomials depends only on "
          "the support]" % len(bad))
    return edges, facets, covers


# ---------------------------------------------------------------------------
# section 2: display (1) is a shelling, hence G is Cohen-Macaulay
# ---------------------------------------------------------------------------

def section_shelling(facets):
    head("2. the order of display (1) is a shelling of Ind(G)")
    derived = shelling_restrictions(facets)
    printed = [set(parse_set(s) for s in row) for row in PAPER_RESTRICTIONS]
    check(len(derived) == len(printed),
          "one restriction list per facet j = 2..10 (derived %d, printed %d)"
          % (len(derived), len(printed)))

    all_match = True
    all_nonempty = True
    all_pure = True
    for j, (d, q) in enumerate(zip(derived, printed), start=2):
        info("j=%2d  F_j = %-4s  derived restriction facets: %s"
             % (j, show_set(facets[j - 1]), show_sets(d)))
        if d != q:
            all_match = False
            info("        PRINTED INSTEAD: %s" % show_sets(q))
        if not d or d == {frozenset()}:
            all_nonempty = False
        if sorted(set(len(S) for S in d)) != [2]:
            all_pure = False
    check(all_match,
          "all nine derived restriction lists agree with the paper "
          "(35; 45; 15,25; 34; 24; 46,47; 13; 36,38; 67,68)")
    check(all_nonempty,
          "every restriction complex is nonempty")
    check(all_pure,
          "every restriction complex is pure of dimension 1, i.e. all its "
          "facets have size 2 = |F_j| - 1")
    FACTS["shellable"] = bool(all_match and all_nonempty and all_pure)
    check(all_match and all_nonempty and all_pure,
          "display (1) is a shelling of the pure 2-dimensional complex "
          "Ind(G); shellable implies Cohen-Macaulay over EVERY field, so "
          "S/I(G) is Cohen-Macaulay over every field")


# ---------------------------------------------------------------------------
# section 3: the minimal generators of J(G)^2
# ---------------------------------------------------------------------------

def section_generators(covers):
    head("3. J(G)^2: the 55 products x_{C_i} x_{C_j} and their minimalisation")
    prods = square_products(covers)
    check(len(prods) == EXP_N_PRODUCTS,
          "there are %d products x_{C_i}x_{C_j} with 1<=i<=j<=10 (derived %d)"
          "  [STRUCTURAL: 10*11/2 is forced once #covers=10 is checked; this "
          "line cannot fail and is not evidence]"
          % (EXP_N_PRODUCTS, len(prods)))
    distinct = sorted(set(m for _, m in prods))
    info("derived number of distinct products: %d" % len(distinct))
    coincidences = {}
    for pair, m in prods:
        coincidences.setdefault(m, []).append(pair)
    repeats = {m: ps for m, ps in coincidences.items() if len(ps) > 1}
    for m, ps in sorted(repeats.items()):
        info("coincidence: %s arises from %s"
             % (show_mon(m), " and ".join("C%dC%d" % p for p in ps)))
    check(len(distinct) == EXP_N_GENERATORS,
          "the 55 products give %d distinct monomials (derived %d)"
          % (EXP_N_GENERATORS, len(distinct)))

    gens = minimalise(distinct)
    info("derived minimal generators of J(G)^2: %d" % len(gens))
    check(len(gens) == len(distinct),
          "no distinct product divides another, so all %d are minimal "
          "generators (derived %d minimal)  [STRUCTURAL: all products of two "
          "size-5 squarefree monomials have degree exactly 10, and distinct "
          "monomials of equal degree cannot divide one another, so this line "
          "cannot fail; the fact it records is still used below]"
          % (len(distinct), len(gens)))
    degs = sorted(set(deg(g) for g in gens))
    FACTS["square_generation_degree"] = degs[0] if len(degs) == 1 else None
    check(degs == [EXP_GEN_DEGREE],
          "every minimal generator of J(G)^2 has degree %d (derived degrees "
          "%s), so J(G)^2 is generated in degree %d"
          % (EXP_GEN_DEGREE, degs, EXP_GEN_DEGREE))
    info("first five generators: " + "; ".join(show_mon(g) for g in gens[:5]))
    info("last  five generators: " + "; ".join(show_mon(g) for g in gens[-5:]))
    return gens


# ---------------------------------------------------------------------------
# section 4: display (2) and the upper Koszul complex K^a
# ---------------------------------------------------------------------------

def section_koszul(covers, gens):
    head("4. display (2) and the upper Koszul complex K^a")
    m = tuple(PAPER_A)
    info("a = %s, m = x^a = %s, |a| = %d"
         % (str(PAPER_A), show_mon(m), deg(m)))

    # The two identities of display (2), as printed there.
    identities = [(("3", "4"), 6, 8), (("5", "6"), 4, 10)]
    for Fs, i, j in identities:
        F = parse_set("".join(Fs))
        lhs = mon_div(m, x_set(F))
        rhs = mon_mul(x_set(covers[i - 1]), x_set(covers[j - 1]))
        info("m/(x_%s x_%s) = %s ;  x_{C_%d} x_{C_%d} = %s"
             % (Fs[0], Fs[1], show_mon(lhs), i, j, show_mon(rhs)))
        check(lhs == rhs,
              "display (2): m/(x_%s x_%s) = x_{C_%d} x_{C_%d}"
              % (Fs[0], Fs[1], i, j))

    K = upper_koszul_complex(m, gens)
    info("derived K^a by testing all 2^8 = 256 subsets F of [8]: %d faces"
         % len(K))
    info("faces of K^a: " + show_sets(K))
    check(sorted(show_set(F) for F in K)
          == sorted(("{}" if s == "" else s) for s in EXP_KOSZUL_FACES),
          "K^a has exactly the faces {}, 3, 4, 5, 6, 34, 56")
    check(is_downward_closed(K),
          "K^a is closed under passing to subsets, hence a simplicial complex"
          "  [STRUCTURAL: true for every monomial ideal and every b, since "
          "m/x_{F'} = (m/x_F) x_{F setminus F'}; the paper proves it as a "
          "lemma, so this line documents rather than tests]")
    gen_faces = maximal_elements(K)
    check(show_sets(gen_faces) == show_sets(parse_set(s)
                                           for s in PAPER_KOSZUL_GENS),
          "K^a = <34,56> (derived facets %s)" % show_sets(gen_faces))
    vs = complex_vertices(K)
    check(vs == parse_set("3456"),
          "the vertex set of K^a is {3,4,5,6}; 1,2,7,8 are not vertices "
          "(derived %s)" % show_set(vs))
    comps = connected_components(K)
    info("derived connected components of K^a: " + show_sets(comps))
    check(len(comps) == EXP_KOSZUL_COMPONENTS,
          "K^a has %d connected components (derived %d), a disjoint union of "
          "two edges" % (EXP_KOSZUL_COMPONENTS, len(comps)))
    for p in FIELDS:
        h = reduced_homology(K, p)
        info("over %-5s reduced homology of K^a: %s"
             % (field_name(p),
                ", ".join("H~_%d=%d" % (k, h[k]) for k in sorted(h))))
        check(h.get(0, 0) == 1 and all(h[k] == 0 for k in h if k != 0),
              "dim H~_0(K^a) = 1 and all other reduced homology vanishes, "
              "over %s" % field_name(p))
    return m, K


# ---------------------------------------------------------------------------
# section 5: table (3) of the paper, re-derived row by row
# ---------------------------------------------------------------------------

def section_table3(covers, gens):
    head("5. table (3): the eight excluded faces of K^a")
    for row in PAPER_TABLE3:
        d = derive_table3_row(row, PAPER_A, covers, gens)
        r, t = d["r"], d["t"]
        info("row F={%d} or %s, printed eligible covers %s, printed t=%d"
             % (r, show_set(d["pq"]),
                ",".join("C%d" % i for i in sorted(d["printed"])), t))
        info("   m/x_%d = %s (x_%d-exponent %d);  m/x_%s = %s (degree %d, "
             "exponent-2 coordinates %s)"
             % (r, show_mon(d["q_singleton"]), r, d["exp_r_in_q1"],
                show_set(d["pq"]), show_mon(d["q_pair"]), d["deg_q_pair"],
                show_set(d["exp2_coords"])))
        info("   eligible from 'omit %d': %s ;  eligible from 'contain %s': %s"
             % (r, ",".join("C%d" % i for i in sorted(d["elig_singleton"])),
                show_set(d["exp2_coords"]),
                ",".join("C%d" % i for i in sorted(d["elig_pair"]))))
        check(d["exp_r_in_q1"] == 0,
              "row {%d}: the x_%d-exponent of m/x_%d is 0, so each factor "
              "cover must omit %d" % (r, r, r, r))
        check(d["deg_q_pair"] == 2 * EXP_COVER_SIZE,
              "row %s: deg(m/x_F) = %d = deg(x_{C_i}x_{C_j}), so divisibility "
              "would be equality  [STRUCTURAL: |a|=12 and |F|=2 are both given, "
              "so 12-2=10 cannot fail; the load-bearing line of this row is the "
              "eligible-cover comparison below]"
              % (show_set(d["pq"]), d["deg_q_pair"]))
        check(d["elig_singleton"] == d["printed"],
              "row {%d}: the covers omitting %d are exactly the printed %s"
              % (r, r, ",".join("C%d" % i for i in sorted(d["printed"]))))
        check(d["elig_pair"] == d["printed"],
              "row %s: the covers containing every exponent-2 coordinate %s "
              "are exactly the printed %s"
              % (show_set(d["pq"]), show_set(d["exp2_coords"]),
                 ",".join("C%d" % i for i in sorted(d["printed"]))))
        check(d["t_in_all"],
              "row {%d}/%s: every eligible cover contains t=%d, so a product "
              "of two of them has x_%d-exponent 2"
              % (r, show_set(d["pq"]), t, t))
        check(d["exp_t_q1"] == 1 and d["exp_t_q2"] == 1,
              "row {%d}/%s: m/x_F has x_%d-exponent 1 for both choices of F "
              "(derived %d and %d)"
              % (r, show_set(d["pq"]), t, d["exp_t_q1"], d["exp_t_q2"]))
        check(not d["q1_in_square"] and not d["q2_in_square"],
              "row {%d}/%s: independently of the table's argument, neither "
              "m/x_{%d} nor m/x_%s is divisible by any of the 53 generators, "
              "so neither F is a face of K^a"
              % (r, show_set(d["pq"]), r, show_set(d["pq"])))


# ---------------------------------------------------------------------------
# section 6: beta_1(J(G)^2) by linear algebra, independent of Thm 1.34
# ---------------------------------------------------------------------------

def section_syzygies(gens):
    head("6. beta_1(J(G)^2) by minimal-first-syzygy linear algebra")
    info("this route uses only exact linear algebra over the stated field; it "
         "does NOT invoke the upper-Koszul formula of Miller-Sturmfels")
    a = tuple(PAPER_A)
    for p in FIELDS:
        rep = {}
        v = beta1_multidegree(a, gens, p, report=rep)
        FACTS.setdefault("beta1_a", {})[field_name(p)] = v
        info("over %-5s at b = a: generators dividing x^a: %s ; dim Z_a = %d ; "
             "rank(sum_t x_t Z_{a-e_t}) = %d ; beta_1 = %d"
             % (field_name(p), rep.get("gen_indices"), rep["dim_Z"],
                rep["rank_mZ"], v))
        check(v == EXP_BETA1_A,
              "beta_{1,a}(J(G)^2) = %d over %s (derived %d)"
              % (EXP_BETA1_A, field_name(p), v))

    # ---- the six fields above are SIX SAMPLES; the paper says "any field".
    # beta1_combinatorial computes the same number with no field at all, so it
    # closes the extrapolation instead of hiding it.
    comb = beta1_combinatorial(a, gens)
    FACTS["beta1_a_combinatorial"] = comb
    info("field-free route (union-find on the index sets I_t, no linear "
         "algebra, no characteristic): beta_{1,a} = %d" % comb)
    check(comb == EXP_BETA1_A,
          "beta_{1,a}(J(G)^2) = %d by the FIELD-FREE route: the spanning set of "
          "(mZ)_a consists only of differences e_i-e_j inside a common I_t, so "
          "its rank is |U|-c over EVERY field; hence this single number holds "
          "in every characteristic, which is what the paper's 'over any field' "
          "needs and what six sampled fields could never give"
          % EXP_BETA1_A)
    check(all(FACTS["beta1_a"][field_name(p)] == comb for p in FIELDS),
          "the six field computations agree with the field-free value, so the "
          "Gaussian elimination is consistent with the combinatorial rank in "
          "every characteristic tested")

    head("6b. exhaustive scan of every multidegree that can carry a first "
         "syzygy")
    for p in FIELDS:
        nz, nlcm = beta1_full_scan(gens, p)
        by_deg = {}
        for b, v in nz.items():
            by_deg[deg(b)] = by_deg.get(deg(b), 0) + v
        total = sum(nz.values())
        FACTS.setdefault("beta1_by_degree", {})[field_name(p)] = dict(by_deg)
        info("over %-5s scanned %d distinct pairwise lcm multidegrees; "
             "%d carry a first syzygy; total beta_1 = %d; by total degree: %s"
             % (field_name(p), nlcm, len(nz), total,
                ", ".join("deg %d -> %d" % (d, by_deg[d])
                          for d in sorted(by_deg))))
        check(by_deg == EXP_TOTAL_BETA1_BY_DEGREE,
              "over %s the first Betti numbers of J(G)^2 by total degree match "
              "the comparison targets %s (derived %s).  NOTE: these "
              "totals are not printed in the paper; "
              "the theorem needs only beta_1 in degree 12 to be nonzero"
              % (field_name(p),
                 str(sorted(EXP_TOTAL_BETA1_BY_DEGREE.items())),
                 str(sorted(by_deg.items()))))
        deg12 = sorted(b for b in nz if deg(b) == 12)
        info("over %-5s degree-12 multidegrees with a first syzygy: %s"
             % (field_name(p), [str(b) for b in deg12]))
        check(deg12 == [a] and nz.get(a) == EXP_BETA1_A,
              "over %s the multidegree a = %s is the UNIQUE degree-12 "
              "multidegree carrying a first syzygy, with beta_1 = %d"
              % (field_name(p), str(a), EXP_BETA1_A))
        check(11 in by_deg and by_deg.get(12, 0) > 0,
              "over %s J(G)^2 is generated in degree 10 yet has a first "
              "syzygy in degree 12 as well as in degree 11, so its resolution "
              "is NOT 10-linear -- a 10-linear resolution would confine all "
              "first syzygies to degree 11" % field_name(p))

    # The same exhaustive scan with no field anywhere, so the sharpness claim
    # ("a is the unique degree-12 multidegree") is also a statement about every
    # field at once and not about six samples.
    nzc, nlcmc = beta1_full_scan_combinatorial(gens)
    by_deg_c = {}
    for b, v in nzc.items():
        by_deg_c[deg(b)] = by_deg_c.get(deg(b), 0) + v
    FACTS["beta1_by_degree_combinatorial"] = dict(by_deg_c)
    info("field-free scan: %d lcms, %d carry a first syzygy, by total degree %s"
         % (nlcmc, len(nzc),
            ", ".join("deg %d -> %d" % (d, by_deg_c[d])
                      for d in sorted(by_deg_c))))
    check(all(FACTS["beta1_by_degree"][field_name(p)] == by_deg_c
              for p in FIELDS),
          "the field-free scan reproduces the by-degree first Betti numbers of "
          "every one of the six fields, confirming that beta_1(J(G)^2) is "
          "field-independent (it is a count of connected components), so the "
          "six-field sample is not an extrapolation")
    deg12c = sorted(b for b in nzc if deg(b) == 12)
    check(deg12c == [a] and nzc.get(a) == EXP_BETA1_A and by_deg_c.get(12, 0) > 0,
          "FIELD-FREE: a = %s is the unique degree-12 multidegree carrying a "
          "first syzygy of J(G)^2, with beta_1 = 1, over EVERY field -- this is "
          "the line that supports the paper's 'no linear resolution over any "
          "field'" % str(a))


# ---------------------------------------------------------------------------
# section 7: the Betti table of J(G) itself, so 5-linearity is decided here
# ---------------------------------------------------------------------------

def section_betti_J(covers):
    head("7. the graded Betti table of J(G), computed directly")
    info("the paper obtains 5-linearity from Eagon-Reiner; here it is decided "
         "from the resolution itself, via upper Koszul complexes over all 256 "
         "squarefree multidegrees")
    jgens = [x_set(C) for C in covers]
    for p in FIELDS:
        tab = betti_table_squarefree(jgens, p)
        info("over %-5s nonzero graded Betti numbers of J(G): %s"
             % (field_name(p),
                ", ".join("beta_{%d,%d}=%d" % (k[0], k[1], tab[k])
                          for k in sorted(tab))))
        check(tab == EXP_BETTI_J,
              "over %s the Betti table of J(G) is beta_{0,5}=10, "
              "beta_{1,6}=13, beta_{2,7}=4 (derived %s).  NOTE: 13 and 4 are "
              "comparison targets fixed here and are not printed in the "
              "paper, which asserts only 5-linearity, "
              "which the next check decides on its own"
              % (field_name(p),
                 ", ".join("beta_{%d,%d}=%d" % (k[0], k[1], tab[k])
                           for k in sorted(tab))))
        FACTS.setdefault("J_linear_shift", {})[field_name(p)] = sorted(
            set(j - i for (i, j) in tab))
        check(bool(tab) and all(j == i + 5 for (i, j) in tab),
              "over %s every nonzero beta_{i,j}(J(G)) has j = i+5, i.e. the "
              "resolution of J(G) is 5-linear" % field_name(p))
        alt = sum((-1) ** i * v for (i, _), v in tab.items())
        check(alt == 1,
              "over %s the alternating sum of the free ranks is %d = rank of "
              "J(G), a consistency check on the resolution"
              % (field_name(p), alt))
        check(sum(v for (i, _), v in tab.items() if i == 0) == len(covers),
              "over %s beta_0(J(G)) equals the number of minimal vertex "
              "covers, %d" % (field_name(p), len(covers)))

    # ---- the upper-Koszul engine above IS \cite[Thm 1.34]{MillerSturmfels},
    # the paper's one uncheckable citation.  Cross-validate its first column
    # against the Miller-Sturmfels-FREE syzygy route on the same ten
    # generators: an index shift or a homology bug would show up here.
    head("7b. cross-validation of the upper-Koszul engine against the "
         "Miller-Sturmfels-free syzygy route, for J(G) itself")
    nzj, nlcmj = beta1_full_scan_combinatorial(jgens)
    by_deg_j = {}
    for b, v in nzj.items():
        by_deg_j[deg(b)] = by_deg_j.get(deg(b), 0) + v
    info("field-free syzygy route on the 10 cover monomials: %d lcms, total "
         "beta_1 = %d, by total degree %s"
         % (nlcmj, sum(nzj.values()),
            ", ".join("deg %d -> %d" % (d, by_deg_j[d])
                      for d in sorted(by_deg_j))))
    koszul_beta1 = {j: v for (i, j), v in
                    betti_table_squarefree(jgens, 0).items() if i == 1}
    check(by_deg_j == koszul_beta1 and bool(by_deg_j),
          "the first Betti numbers of J(G) by degree agree between the "
          "upper-Koszul engine (%s) and the syzygy engine (%s); the two routes "
          "are independent, so this validates the homological index convention "
          "of Theorem 1.34 as implemented rather than assuming it"
          % (str(sorted(koszul_beta1.items())), str(sorted(by_deg_j.items()))))
    check(sorted(by_deg_j) == [6],
          "every first syzygy of J(G) sits in degree 6 = 5+1 (field-free), so "
          "the first step of the resolution is 5-linear over EVERY field, not "
          "only over the six tested")


# ---------------------------------------------------------------------------
# section 8: the remark
# ---------------------------------------------------------------------------

def section_remark(edges, gens):
    head("8. remark: non-chordality, and the symbolic square")
    eset = set(edges)
    quad = [1, 6, 5, 7]
    present = sorted((i, j) for i, j in eset if i in quad and j in quad)
    info("edges induced on {1,5,6,7}: " +
         ",".join("%d%d" % e for e in present))
    cyc = [tuple(sorted((quad[k], quad[(k + 1) % 4]))) for k in range(4)]
    diag = [tuple(sorted((quad[0], quad[2]))), tuple(sorted((quad[1], quad[3])))]
    info("cycle edges 1-6-5-7-1: %s ; diagonals: %s"
         % (",".join("%d%d" % e for e in cyc),
            ",".join("%d%d" % e for e in diag)))
    check(all(e in eset for e in cyc) and all(e not in eset for e in diag)
          and len(present) == 4,
          "the induced subgraph on {1,5,6,7} is exactly the 4-cycle "
          "1-6-5-7-1 (the two diagonals 15 and 67 are absent), so G has a "
          "chordless 4-cycle and is not chordal")

    full = x_set(frozenset(VERTS))
    info("x_1...x_8 = %s, degree %d" % (show_mon(full), deg(full)))
    # J(G)^{(2)} = intersection over edges of (x_i,x_j)^2, the standard
    # description of the symbolic power of a cover ideal.  For a monomial u,
    # u in (x_i,x_j)^2 iff x_i^2 | u or x_i x_j | u or x_j^2 | u.
    ok_sym = True
    for i, j in edges:
        gens_ij = [mon_mul(x_set(frozenset([i])), x_set(frozenset([i]))),
                   mon_mul(x_set(frozenset([i])), x_set(frozenset([j]))),
                   mon_mul(x_set(frozenset([j])), x_set(frozenset([j])))]
        if not in_ideal(full, gens_ij):
            ok_sym = False
    check(ok_sym,
          "x_1...x_8 lies in (x_i,x_j)^2 for every one of the %d edges, hence "
          "in the symbolic square J(G)^{(2)} = intersection of (x_i,x_j)^2"
          % len(edges))
    check(not in_ideal(full, gens),
          "x_1...x_8 is divisible by none of the 53 generators of J(G)^2, so "
          "it is not in the ordinary square; therefore "
          "J(G)^2 != J(G)^{(2)}")
    check(deg(full) < EXP_GEN_DEGREE,
          "indeed deg(x_1...x_8) = %d < %d = the generation degree of J(G)^2"
          "  [STRUCTURAL: 8 < 10, two constants; cannot fail]"
          % (deg(full), EXP_GEN_DEGREE))

    # Remark, second half: the paper distinguishes its example from
    # Ficarra-Moradi by asserting that a graph cover ideal has gcd 1.  That is
    # checkable here for this G.
    gcd_exp = tuple(min(x_set(C)[k] for C in [parse_set(s) for s in
                                              PAPER_COVERS])
                    for k in range(N))
    missed = [v for v in VERTS if gcd_exp[v - 1] == 0]
    info("gcd of the ten cover monomials: %s ; vertices omitted by some cover: "
         "%s" % (show_mon(gcd_exp), "".join(str(v) for v in missed)))
    check(deg(gcd_exp) == 0,
          "the ten generators of J(G) have greatest common divisor 1 (every "
          "vertex is omitted by some minimal vertex cover, i.e. extends to a "
          "maximal independent set) -- the property the Remark uses to "
          "separate this example from Ficarra-Moradi's")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print("verification of: A Negative Answer to Mohammadi's Question on "
          "Powers of Cover Ideals")
    print("all arithmetic is exact (integers, fractions.Fraction over Q, "
          "integers mod p over F_p); no check is numeric, so no error bound "
          "is needed")
    print("fields used: " + ", ".join(field_name(p) for p in FIELDS))

    edges, facets, covers = section_graph()
    section_shelling(facets)
    gens = section_generators(covers)
    section_koszul(covers, gens)
    section_table3(covers, gens)
    section_syzygies(gens)
    section_betti_J(covers)
    section_remark(edges, gens)

    head("summary: the theorem's three assertions, from the recorded facts")
    shellable = FACTS.get("shellable") is True
    gd = FACTS.get("square_generation_degree")
    lin = FACTS.get("J_linear_shift", {})
    b1a = FACTS.get("beta1_a", {})
    bybd = FACTS.get("beta1_by_degree", {})
    fields = [field_name(p) for p in FIELDS]
    info("shellable = %s" % shellable)
    info("generation degree of J(G)^2 = %s" % gd)
    info("linear shifts of the resolution of J(G), by field: %s" % lin)
    info("beta_{1,a}(J(G)^2), by field: %s" % b1a)
    info("beta_1(J(G)^2) by total degree, by field: %s" % bybd)

    check(shellable,
          "assertion 1: Ind(G) is shellable, so G is Cohen-Macaulay over "
          "every field")
    check(bool(lin) and all(lin.get(f) == [5] for f in fields),
          "assertion 2: over every field tested, the resolution of J(G) is "
          "5-linear")
    check(bool(b1a) and all(b1a.get(f) == 1 for f in fields)
          and FACTS.get("beta1_a_combinatorial") == 1,
          "assertion 3: beta_{1,a}(J(G)^2) = 1 at a = (1,1,2,2,2,2,1,1), "
          "|a| = %d -- over the six fields tested AND, by the field-free "
          "union-find route, over every field" % sum(PAPER_A))
    check(gd == 10 and sum(PAPER_A) == 12 and bool(bybd)
          and all(bybd.get(f, {}).get(12, 0) > 0 for f in fields)
          and FACTS.get("beta1_by_degree_combinatorial", {}).get(12, 0) > 0,
          "conclusion: J(G)^2 is generated in degree %s but carries a first "
          "syzygy in degree 12, not 11, over every field tested, so J(G)^2 "
          "has no linear resolution over any field; Mohammadi's Question 4.1, "
          "restated as Ha-Van Tuyl Question 6.2, has a negative answer" % gd)

    head("what this program does NOT establish (stated, not hidden)")
    info("1. The wording of Mohammadi's Question 4.1 (Collect. Math. 65 (2014), "
         "paywalled, no preprint) and of the Ha-Van Tuyl Question 6.2 "
         "restatement quoted in the paper are TEXTUAL claims; no program can "
         "check them.  What is checked here is that the mathematical object "
         "satisfies the hypothesis (CM, 5-linear cover ideal) and violates the "
         "conclusion (square not linear).")
    info("2. \\cite[Theorem~1.34]{MillerSturmfels} is used as the definition of "
         "the upper-Koszul engine in sections 4 and 7.  It is NOT assumed for "
         "the headline: beta_{1,a}(J(G)^2)=1 is decided again in section 6 by "
         "syzygy linear algebra and once more field-free, and section 7b shows "
         "the two engines agree on beta_1(J(G)).  Higher Betti numbers of J(G) "
         "(beta_{2,7}=4) still rest on Theorem 1.34 as implemented.")
    info("3. 'Cohen-Macaulay over every field' is carried by SHELLABILITY, "
         "which is field-free and is checked in section 2; the step "
         "shellable => Cohen-Macaulay, and Eagon-Reiner (CM <=> linear "
         "resolution of the dual), are cited theorems, not verified here.  The "
         "six-field Betti table of section 7 corroborates but cannot by itself "
         "reach 'every field'; section 7b closes the beta_1 column field-free.")
    info("4. The Remark's comparison with Ficarra-Moradi (their Theorem D, the "
         "hitting set {a,c,f,y_1,...,y_{n-6}}, non-quadratic Alexander dual) is "
         "about THEIR construction and is not reproduced here.  Only the two "
         "claims about G itself are checked: the induced 4-cycle 1-6-5-7-1 and "
         "gcd of the generators = 1.")
    info("5. Only beta_0 and beta_1 are computed for J(G)^2.  The paper needs "
         "no more: generation in degree 10 plus one first syzygy in degree 12 "
         "already contradicts linearity.  Nothing here computes the full "
         "resolution of J(G)^2, and nothing needs to.")

    n = len(RESULTS)
    failed = n - sum(1 for r in RESULTS if r)
    print("")
    if failed == 0:
        print("VERDICT: ALL %d CHECKS PASS" % n)
        return 0
    print("VERDICT: %d OF %d CHECKS FAILED" % (failed, n))
    return 1


if __name__ == "__main__":
    sys.exit(main())

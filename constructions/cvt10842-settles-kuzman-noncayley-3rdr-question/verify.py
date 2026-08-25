#!/usr/bin/env python3
"""Verification of a non-Cayley vertex-transitive 3-rainbow-domination-regular
cubic graph.  Standard library only; every decision is an exact integer or set
comparison, no floating point anywhere.

TAKEN FROM THE PAPER (inputs, transcribed verbatim, nothing else assumed):
  * the sparse6 string of the exhibited graph X (four appendix lines joined);
  * the three colour classes C1, C2, C3 of the exhibited function, in the
    paper's labelling 1..108;
  * the numbers the paper asserts, kept in CLAIM_* only so that they can be
    contradicted: |V(X)|=108, degree 3, gamma_{r3}(X)=54, |Aut(X)|=216,
    |[A,A]|=54, three subgroups of index 2, |C_i|=18.

DERIVED HERE (recomputed from the witness alone):
  * decoding: order, edge count, simplicity of the raw record stream, degrees,
    connectivity, and that every vertex has the same finite eccentricity;
  * Aut(X) by exhaustive backtracking over every image of a BFS vertex order,
    then verified elementwise and verified to be closed under product and
    inverse; its order, its orbits, the vertex stabiliser order;
  * that f meets the definition of a 3-rainbow dominating function, its
    weight, its minimality, and the structural claim that each uncoloured
    vertex has exactly one neighbour in each class;
  * a LOWER bound: X box K_3 is verified here to be 5-regular on 324
    vertices, so any dominating set of it has size >= ceil(324/6) = 54.  The
    step that turns this into a bound on gamma_{r3}, namely that a 3RDF of
    weight w gives a dominating set of size w in X box K_3 (equivalently
    gamma_{r3}(G) = gamma(G box K_3)), is a LEMMA of the literature: it is
    cited, not machine-proved here.  It is exercised on a whole family of
    3-rainbow dominating functions of X (every automorphic image and every
    one-colour extension of the exhibited f; the family size is derived and
    printed) and confirmed as an exact equality on the small cubic fixtures.
    The exhibited function attains 54 (indeed as a
    perfect dominating set), so gamma_{r3}(X) = 54 = |V(X)|/2 and X is 3-RDR;
  * an independent randomised search for a lighter function, which reaches 54
    and never less;
  * the two ingredients of that bound tested independently on small cubic
    graphs: no 3RDF of weight |V|/2 - 1 exists on any of them, and their
    gamma_{r3} equals the domination number of their product with K_3, both
    by exhaustive enumeration;
  * the derived subgroup, the quotient A/[A,A], and a COMPLETE census of the
    subgroups of order 108: such a subgroup has index 2, so it contains every
    square and every commutator and is a union of cosets of R = <squares,
    commutators>, which makes a coset-level enumeration exhaustive; each
    survivor is materialised and re-verified closed, and independently
    confirmed to be a union of conjugacy classes.  All three are intransitive,
    so Aut(X) has no regular subgroup
    and X is not a Cayley graph.  The same census is run on a circulant that
    is a Cayley graph, where it does find regular subgroups, so the negative
    result is not an artefact of the method.

NOT RE-RUN: the identification of X with a specific entry of the published
cubic vertex-transitive census -- asserted in the paper's abstract and in the
first sentence of its proof -- which needs that external catalogue; and the
general lemma named above.

RUNTIME: a few minutes under CPython.  Nothing is printed while the randomised
descent (check 15) and the exhaustive small-cubic searches (checks 17 and 18)
run; those three account for nearly all of the time.
"""

import sys
from collections import deque
from itertools import combinations

# ---------------------------------------------------------------- paper inputs
# Sparse6 encoding of the exhibited graph, transcribed from the appendix
# (four strings concatenated without whitespace).
S6 = (":~?@k_GAA_WAC`WM@_gWGaWiF_wQDb?yNcHEL`g]GaWgPcpOUd`]Wc`e"
      "YexUJbGuMGw}OEhw_GiGchY]ggXMSdX]WeQ]jJAymJzApkjMsjYE\\fh}"
      "_grIcOIS}hbMvMRqwMbu{Njs~oKFAoz`CpZlNjYijjI{okTHRqcpOqst"
      "PsTJLTKlTmzXHlzeyVz|@VlhZVLp\\rK|bqDfIre@cstPfvd|dxEXf")

# The colour classes of the exhibited 3-rainbow dominating function, in the
# paper's labelling 1..108 (decoded sparse6 labels 0..107 shifted by +1).
C1 = [4, 5, 18, 19, 23, 33, 52, 54, 60, 61, 67, 74, 78, 79, 83, 86, 104, 107]
C2 = [8, 10, 11, 12, 21, 35, 40, 55, 58, 63, 65, 72, 80, 82, 84, 98, 99, 101]
C3 = [1, 20, 24, 26, 31, 34, 37, 38, 39, 41, 50, 56, 76, 85, 88, 97, 100, 103]

CLAIM_ORDER = 108          # |V(X)|
CLAIM_DEGREE = 3           # cubic
CLAIM_GAMMA_R3 = 54        # gamma_{r3}(X)
CLAIM_AUT = 216            # |Aut(X)|
CLAIM_DERIVED = 54         # |[A,A]|
CLAIM_INDEX2 = 3           # number of subgroups of index 2
CLAIM_CLASS_SIZE = 18      # |C_i|

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))


def finish():
    n = len(CHECKS)
    bad = [nm for nm, ok in CHECKS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


# --------------------------------------------------------------- sparse6 decode
def decode_sparse6(text):
    """Decode a sparse6 string into (n, sorted distinct edges, raw records).

    The raw record list is returned separately and NOT deduplicated, so that
    a loop or a repeated edge in the encoded stream stays observable instead
    of being silently absorbed by a set.
    """
    if text[0] != ":":
        raise ValueError("not sparse6")
    vals = []
    for ch in text[1:]:
        v = ord(ch) - 63
        if v < 0 or v > 63:
            raise ValueError("byte out of range")
        vals.append(v)
    if vals[0] == 63:
        if vals[1] == 63:
            n = 0
            for j in range(2, 8):
                n = (n << 6) | vals[j]
            pos = 8
        else:
            n = (vals[1] << 12) | (vals[2] << 6) | vals[3]
            pos = 4
    else:
        n = vals[0]
        pos = 1
    k = 1
    while (1 << k) < n:
        k += 1
    bits = []
    for v in vals[pos:]:
        for j in range(5, -1, -1):
            bits.append((v >> j) & 1)
    raw = []
    v = 0
    p = 0
    while p + k + 1 <= len(bits):
        b = bits[p]
        p += 1
        x = 0
        for _ in range(k):
            x = (x << 1) | bits[p]
            p += 1
        if b == 1:
            v += 1
        if x > v:
            v = x
        else:
            raw.append((min(x, v), max(x, v)))
        if v >= n:
            break
    return n, sorted(set(raw)), raw


# ------------------------------------------------------------ graph utilities
def adjacency(n, edges):
    """Adjacency lists; raises on a loop or a repeated edge."""
    adj = [[] for _ in range(n)]
    seen = set()
    for a, b in edges:
        if a == b:
            raise ValueError("loop")
        if (a, b) in seen:
            raise ValueError("repeated edge")
        seen.add((a, b))
        adj[a].append(b)
        adj[b].append(a)
    return [sorted(s) for s in adj]


def component_size(adj, start=0):
    seen = {start}
    stack = [start]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen)


def distance_rows(adj):
    n = len(adj)
    rows = []
    for s in range(n):
        d = [-1] * n
        d[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for w in adj[u]:
                if d[w] < 0:
                    d[w] = d[u] + 1
                    q.append(w)
        rows.append(d)
    return rows


def two_colouring(adj):
    """Return (is_bipartite, part sizes) computed by BFS 2-colouring."""
    n = len(adj)
    col = [-1] * n
    ok = True
    for s in range(n):
        if col[s] >= 0:
            continue
        col[s] = 0
        stack = [s]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if col[w] < 0:
                    col[w] = 1 - col[u]
                    stack.append(w)
                elif col[w] == col[u]:
                    ok = False
    return ok, (col.count(0), col.count(1))


# ------------------------------------------------- 3-rainbow domination basics
FULL = 7  # bitmask of {1,2,3}


def popcount(x):
    return bin(x).count("1")


def weight(f):
    return sum(popcount(x) for x in f)


def is_3rdf(adj, f):
    """Definition check: f(v) empty  =>  union of f over N(v) is {1,2,3}."""
    for v in range(len(adj)):
        if f[v] == 0:
            u = 0
            for w in adj[v]:
                u |= f[w]
            if u != FULL:
                return False
    return True


def product_with_k3(adj):
    """Cartesian product X box K_3; vertex (v,i) is indexed 3*v+i."""
    n = len(adj)
    h = [[] for _ in range(3 * n)]
    for v in range(n):
        for i in range(3):
            a = 3 * v + i
            for j in range(3):
                if j != i:
                    h[a].append(3 * v + j)
            for u in adj[v]:
                h[a].append(3 * u + i)
    return [sorted(s) for s in h]


def domination_multiplicities(h, dset):
    """For each vertex of h, how many members of dset lie in its closed nbhd."""
    m = [0] * len(h)
    ds = set(dset)
    for a in range(len(h)):
        c = 1 if a in ds else 0
        for b in h[a]:
            if b in ds:
                c += 1
        m[a] = c
    return m


def rdf_to_dominating_set(f):
    """(v,i) with colour i in f(v); a 3RDF of weight w gives w vertices."""
    out = []
    for v in range(len(f)):
        for i in range(3):
            if f[v] >> i & 1:
                out.append(3 * v + i)
    return out


# -------------------------------------------------- exhaustive / heuristic search
def exists_3rdf_of_weight(adj, m):
    """True iff some 3RDF of weight exactly m exists (exhaustive).

    Weight-m functions are exactly the m-subsets of V x {1,2,3}; adding a
    colour never destroys the defining implication, so testing weight exactly
    m also settles all weights below m.
    """
    n = len(adj)
    slots = [(v, i) for v in range(n) for i in range(3)]
    for pick in combinations(range(len(slots)), m):
        f = [0] * n
        for s in pick:
            v, i = slots[s]
            f[v] |= 1 << i
        if is_3rdf(adj, f):
            return True
    return False


def descent_search(adj, rounds, kicks, seed=999):
    """Deterministic randomised search for the lightest 3RDF (kick descent).

    Independent of the exhibited function: starts from noise, repairs the
    definition, strips every removable colour, then repeatedly deletes three
    random colours and repairs again.  Returns the best weight seen.
    """
    n = len(adj)
    state = seed
    best = 3 * n

    def rnd(k):
        nonlocal state
        state = (1103515245 * state + 12345) % (1 << 31)
        return state % k

    def seen_by(f, v):
        u = 0
        for w in adj[v]:
            u |= f[w]
        return u

    def repair(f):
        for _ in range(40 * n):
            bad = [v for v in range(n)
                   if f[v] == 0 and seen_by(f, v) != FULL]
            if not bad:
                return True
            v = bad[rnd(len(bad))]
            have = seen_by(f, v)
            miss = [i for i in range(3) if not (have >> i & 1)]
            tgt = v if rnd(4) == 0 or not adj[v] else adj[v][rnd(len(adj[v]))]
            f[tgt] |= 1 << miss[rnd(len(miss))]
        return is_3rdf(adj, f)

    def ok_around(f, v):
        """The 3RDF condition at v and at its neighbours only.

        Deleting a colour at v changes f only at v, so the defining
        implication can only break at v itself (if f(v) has just become
        empty) or at a neighbour of v: for any other vertex x the union of f
        over N(x) is untouched.  Given that f satisfied the definition before
        the deletion, this local test is therefore EXACTLY equivalent to a
        full is_3rdf sweep, at a fraction of the cost.
        """
        for x in (v,) + tuple(adj[v]):
            if f[x] == 0:
                u = 0
                for w in adj[x]:
                    u |= f[w]
                if u != FULL:
                    return False
        return True

    def strip(f):
        # f satisfies the definition on entry (repair() has just returned
        # True), and every deletion below is either kept -- leaving f a 3RDF --
        # or immediately undone, so the precondition of ok_around holds
        # throughout.
        for v in range(n):
            for i in range(3):
                if f[v] >> i & 1:
                    f[v] &= ~(1 << i)
                    if not ok_around(f, v):
                        f[v] |= 1 << i

    for _ in range(rounds):
        f = [0] * n
        for v in range(n):
            if rnd(4) == 0:
                f[v] = 1 << rnd(3)
        if not repair(f):
            continue
        strip(f)
        cur = weight(f)
        for _ in range(kicks):
            g = list(f)
            for _ in range(3):
                live = [(v, i) for v in range(n) for i in range(3)
                        if g[v] >> i & 1]
                v, i = live[rnd(len(live))]
                g[v] &= ~(1 << i)
            if not repair(g):
                continue
            strip(g)
            w = weight(g)
            if w <= cur:
                f, cur = g, w
        if cur < best:
            best = cur
    return best


# ------------------------------------------------------- automorphism group
def all_automorphisms(adj, cap=100000):
    """Exhaustive backtracking over every image of a BFS-ordered vertex list.

    Every vertex after the first has an already-mapped neighbour, so its image
    ranges over the (at most 3) neighbours of that image; partial maps are
    pruned by requiring all pairwise distances to agree.  Nothing outside the
    search space can be an automorphism, so the census is complete.
    """
    n = len(adj)
    dist = distance_rows(adj)
    order = [0]
    parent = [None] * n
    seen = [False] * n
    seen[0] = True
    q = deque([0])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True
                parent[w] = u
                order.append(w)
                q.append(w)
    if len(order) != n:
        raise ValueError("graph is disconnected")
    img = [-1] * n
    used = [False] * n
    found = []

    def rec(i):
        if len(found) > cap:
            raise ValueError("automorphism cap exceeded")
        if i == n:
            found.append(tuple(img))
            return
        u = order[i]
        du = dist[u]
        for c in adj[img[parent[u]]]:
            if used[c]:
                continue
            dc = dist[c]
            ok = True
            for j in range(i):
                w = order[j]
                if du[w] != dc[img[w]]:
                    ok = False
                    break
            if ok:
                img[u] = c
                used[c] = True
                rec(i + 1)
                used[c] = False
                img[u] = -1

    for w0 in range(n):
        img[0] = w0
        used[w0] = True
        rec(1)
        used[w0] = False
        img[0] = -1
    return found


# --------------------------------------------------- permutation group toolkit
def pmul(p, q):
    return tuple(p[x] for x in q)


def pinv(p):
    r = [0] * len(p)
    for i, x in enumerate(p):
        r[x] = i
    return tuple(r)


def preserves_edges(adj, p):
    for v in range(len(adj)):
        if sorted(p[x] for x in adj[v]) != adj[p[v]]:
            return False
    return True


def generated(gens, ident):
    grp = {ident}
    frontier = [ident]
    gens = list(gens)
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = pmul(g, a)
                if b not in grp:
                    grp.add(b)
                    nxt.append(b)
        frontier = nxt
    return grp


def is_closed(elems):
    s = set(elems)
    for a in s:
        if pinv(a) not in s:
            return False
        for b in s:
            if pmul(a, b) not in s:
                return False
    return True


def orbit_sizes(elems, n):
    lab = [-1] * n
    sizes = []
    for v in range(n):
        if lab[v] < 0:
            orb = set()
            for p in elems:
                orb.add(p[v])
            for x in orb:
                lab[x] = len(sizes)
            sizes.append(len(orb))
    return sorted(sizes)


def left_cosets(group, sub):
    cs = []
    rest = set(group)
    while rest:
        g = min(rest)
        c = frozenset(pmul(g, h) for h in sub)
        cs.append(c)
        rest -= c
    return cs


def conjugacy_classes(group):
    elems = sorted(group)
    classes = []
    left = set(elems)
    while left:
        g = min(left)
        c = frozenset(pmul(pmul(x, g), pinv(x)) for x in elems)
        classes.append(c)
        left -= c
    return classes


def residual(group, ident):
    """R = <squares, commutators>, i.e. the smallest subgroup with elementary
    abelian 2-quotient.  Returned together with the commutator set."""
    elems = sorted(group)
    comms = set()
    for a in elems:
        ia = pinv(a)
        for b in elems:
            comms.add(pmul(pmul(ia, pinv(b)), pmul(a, b)))
    return comms, generated(comms | {pmul(g, g) for g in elems}, ident)


def subgroups_of_order(group, m, resid=None):
    """Complete census of the subgroups of order m when |group|/m equals 2.

    A subgroup H of index 2 is the kernel of a homomorphism onto C_2, so H
    contains every square and every commutator.  Hence R = <squares,
    commutators> lies in H, H is a union of cosets of R, and H/R is an
    index-2 subgroup of the (elementary abelian, hence small) quotient
    group/R.  Enumerating the index-2 subgroups of that quotient at COSET
    level is therefore a complete census, and it costs O(q^2) in the index
    q = |group/R| instead of being exponential in the number of conjugacy
    classes.  Every survivor is still materialised and re-checked with
    is_closed, so nothing is taken on trust.
    """
    if 2 * m != len(group):
        raise ValueError("this census only applies to index 2")
    ident = tuple(range(len(next(iter(group)))))
    if resid is None:
        resid = residual(group, ident)[1]
    cos = left_cosets(group, resid)
    q = len(cos)
    if q % 2 or q > 16:
        raise ValueError("quotient by <squares, commutators> has index %d" % q)
    where = {}
    for i, c in enumerate(cos):
        for g in c:
            where[g] = i
    reps = [min(c) for c in cos]
    tab = [[where[pmul(reps[i], reps[j])] for j in range(q)] for i in range(q)]
    e = where[ident]
    out = []
    for pick in combinations([i for i in range(q) if i != e], q // 2 - 1):
        s = set(pick) | {e}
        if not all(tab[i][j] in s for i in s for j in s):
            continue
        h = set()
        for i in s:
            h |= set(cos[i])
        if len(h) == m and is_closed(h):
            out.append(frozenset(h))
    return out


def min_3rdf_weight(adj):
    """Exact gamma_{r3} by exhaustive search; only for tiny graphs."""
    for m in range(3 * len(adj) + 1):
        if exists_3rdf_of_weight(adj, m):
            return m
    return None


def min_dominating_size(h):
    """Exact domination number by exhaustive search; only for tiny graphs.

    The search space is unchanged -- every subset of every size, smallest size
    first -- but a subset is tested by OR-ing precomputed closed-neighbourhood
    bitmasks, which is the same predicate as
    min(domination_multiplicities(h, sub)) >= 1 written with exact integers.
    """
    n = len(h)
    masks = []
    for a in range(n):
        msk = 1 << a
        for b in h[a]:
            msk |= 1 << b
        masks.append(msk)
    full = (1 << n) - 1
    for m in range(n + 1):
        for sub in combinations(range(n), m):
            u = 0
            for a in sub:
                u |= masks[a]
            if u == full:
                return m
    return None


def small_cubic_fixtures():
    """A few connected cubic graphs, built here, for the n <= 2w lemma."""
    k4 = [sorted(set(range(4)) - {v}) for v in range(4)]
    k33 = [[3, 4, 5], [3, 4, 5], [3, 4, 5], [0, 1, 2], [0, 1, 2], [0, 1, 2]]
    prism = [[1, 2, 3], [0, 2, 4], [0, 1, 5],
             [0, 4, 5], [1, 3, 5], [2, 3, 4]]
    cube = [[] for _ in range(8)]
    for v in range(8):
        cube[v] = sorted(v ^ (1 << b) for b in range(3))
    pet = [[] for _ in range(10)]
    for i in range(5):
        pet[i] += [(i + 1) % 5, (i + 4) % 5, 5 + i]
        pet[5 + i] += [5 + (i + 2) % 5, 5 + (i + 3) % 5, i]
    pet = [sorted(s) for s in pet]
    return [("K4", k4), ("K33", k33), ("prism", prism),
            ("Q3", cube), ("Petersen", pet)]


# Built once: the two small-graph checks below share a single exhaustive pass
# over these fixtures instead of each rebuilding and re-searching them.
FIXTURES = small_cubic_fixtures()


def main():
    # ---- 1. the exhibited object: decode it, count it, print it back
    n, edges, raw = decode_sparse6(S6)
    adj = adjacency(n, edges)
    degs = sorted(len(a) for a in adj)
    ck("witness_decodes_to_cubic_graph_of_order_108",
       n == CLAIM_ORDER and degs[0] == degs[-1] == CLAIM_DEGREE
       and len(edges) == 3 * CLAIM_ORDER // 2,
       "n=%d edges=%d degrees from %d to %d"
       % (n, len(edges), degs[0], degs[-1]))
    ck("witness_is_connected",
       component_size(adj) == n,
       "component covers %d of %d vertices" % (component_size(adj), n))
    rows = distance_rows(adj)
    ecc = sorted(set(max(r) for r in rows))
    reachable = all(d >= 0 for r in rows for d in r)
    simple = len(set(raw)) == len(raw) and not any(a == b for a, b in raw)
    ck("witness_is_simple_and_has_one_common_finite_eccentricity",
       simple and reachable and len(ecc) == 1 and ecc[0] > 0,
       "raw records %d all distinct and loop-free: %s; all distances finite: "
       "%s; eccentricities occurring: %s"
       % (len(raw), simple, reachable, ecc))

    # ---- 2. hypotheses of the statement: vertex-transitive cubic graph
    auts = all_automorphisms(adj)
    aut_set = set(auts)
    valid = all(preserves_edges(adj, p) for p in auts)
    closed = is_closed(aut_set)
    ck("automorphism_census_is_a_group_of_order_216",
       len(aut_set) == len(auts) == CLAIM_AUT and valid and closed,
       "|Aut(X)|=%d, every member adjacency preserving: %s, closed under "
       "product and inverse: %s" % (len(aut_set), valid, closed))
    ck("X_is_vertex_transitive",
       orbit_sizes(auts, n) == [n],
       "orbit sizes %s" % (orbit_sizes(auts, n),))
    stab = [p for p in auts if p[0] == 0]
    ck("vertex_stabiliser_has_order_2",
       len(stab) == len(auts) // n == 2,
       "|Aut(X)|/|V|=%d, computed stabiliser order %d"
       % (len(auts) // n, len(stab)))

    # ---- 3. the colour classes and the exhibited function
    sets = [C1, C2, C3]
    flat = C1 + C2 + C3
    sizes = [len(set(c)) for c in sets]
    ck("colour_classes_are_disjoint_18_sets_of_valid_labels",
       sizes == [CLAIM_CLASS_SIZE] * 3 and len(set(flat)) == 54
       and min(flat) >= 1 and max(flat) <= n,
       "sizes %s, union %d, labels in [%d,%d]"
       % (sizes, len(set(flat)), min(flat), max(flat)))
    f = [0] * n
    for i, c in enumerate(sets):
        for lab in c:
            f[lab - 1] |= 1 << i
    empty = [v for v in range(n) if f[v] == 0]
    offenders = [v for v in empty
                 if any(len([u for u in adj[v] if f[u] >> i & 1]) != 1
                        for i in range(3))]
    ck("every_uncoloured_vertex_has_one_neighbour_in_each_class",
       not offenders and len(empty) == 54,
       "%d vertices have f(v) empty, %d of them fail the one-per-class count"
       % (len(empty), len(offenders)))
    legal = is_3rdf(adj, f)
    ck("f_is_a_3_rainbow_dominating_function_of_weight_54",
       legal and weight(f) == CLAIM_GAMMA_R3,
       "definition satisfied everywhere: %s, w(f)=%d" % (legal, weight(f)))
    bip, parts = two_colouring(adj)
    oneside = all(f[u] == 0 for v in range(n) if f[v] for u in adj[v])
    ck("coloured_set_and_its_complement_are_the_two_bipartition_classes",
       bip and parts == (54, 54) and oneside,
       "bipartite: %s with parts %s; coloured vertices confined to one side: "
       "%s" % (bip, parts, oneside))
    return n, adj, f, auts, aut_set


def check_gamma(n, adj, f, auts):
    """Lower bound via domination in X box K_3, then gamma_{r3}(X)=54."""
    h = product_with_k3(adj)
    hdeg = sorted(len(a) for a in h)
    lower = -(-len(h) // (hdeg[-1] + 1))     # ceiling division, exact integers
    ck("product_with_K3_is_5_regular_on_324_vertices",
       len(h) == 3 * n and hdeg[0] == hdeg[-1] == 5,
       "|V|=%d, every closed neighbourhood has %d vertices"
       % (len(h), hdeg[-1] + 1))
    dset = rdf_to_dominating_set(f)
    mult = domination_multiplicities(h, dset)
    ck("image_of_f_is_a_perfect_dominating_set_of_the_product",
       len(dset) == weight(f) and min(mult) == max(mult) == 1,
       "|D|=%d, every one of %d vertices dominated exactly %d time(s)"
       % (len(dset), len(h), max(mult)))
    # The step from the degree bound to gamma_r3 >= 54 -- that a 3RDF of weight
    # w yields a dominating set of size w in X box K_3 (equivalently
    # gamma_r3(G) = gamma(G box K_3)) -- is a lemma of the literature.  It is
    # NOT machine-proved here for all functions, and the detail line below says
    # so.  What IS machine-checked is the lemma's conclusion on a whole family
    # of 3RDFs of X, generated without reference to the paper's arithmetic:
    # the image of f under every automorphism of X, and every extension of f by
    # one further colour.  Each member must itself satisfy the definition and
    # must map to a dominating set of the product of exactly its own weight, so
    # a wrong product, a wrong image map or a broken witness fails here.
    family = []
    for p in auts:
        g = [0] * n
        for v in range(n):
            g[p[v]] = f[v]
        family.append(g)
    for v in range(n):
        for i in range(3):
            if not (f[v] >> i & 1):
                g = list(f)
                g[v] |= 1 << i
                family.append(g)
    fam_ok = 0
    for g in family:
        d = rdf_to_dominating_set(g)
        if (is_3rdf(adj, g) and len(d) == weight(g)
                and min(domination_multiplicities(h, d)) >= 1):
            fam_ok += 1
    fam_expected = len(auts) + 3 * n - weight(f)
    ck("degree_bound_forces_gamma_r3_at_least_54",
       lower == CLAIM_GAMMA_R3
       and fam_ok == len(family) == fam_expected,
       "ceil(%d/%d)=%d; the step 'a 3RDF of weight w yields a dominating set "
       "of size w in the product' is a cited lemma, not proved here, and its "
       "conclusion holds for %d of %d test functions (every automorphic image "
       "and every one-colour extension of f)"
       % (len(h), hdeg[-1] + 1, lower, fam_ok, len(family)))
    upper = weight(f) if is_3rdf(adj, f) else 3 * n
    ck("gamma_r3_equals_54_and_X_is_3_RDR",
       lower == upper == CLAIM_GAMMA_R3 == n // 2,
       "bounds meet: %d <= gamma_r3 <= %d, and |V(X)|/2 = %d"
       % (lower, upper, n // 2))
    best = descent_search(adj, 20, 600)
    ck("independent_search_attains_54_and_never_less",
       best == CLAIM_GAMMA_R3,
       "kick descent from random starts, ignorant of the paper's function, "
       "reaches weight %d" % best)
    survivors = []
    for v in range(n):
        for i in range(3):
            if f[v] >> i & 1:
                g = list(f)
                g[v] &= ~(1 << i)
                if is_3rdf(adj, g):
                    survivors.append((v + 1, i + 1))
    # is_3rdf(f) is part of the condition on purpose: for an f that is not a
    # 3RDF at all, no single deletion can make it one, so "no survivors" would
    # otherwise hold vacuously and this check would pass on a broken witness.
    ck("f_is_minimal_no_single_colour_can_be_dropped",
       is_3rdf(adj, f) and not survivors and weight(f) == CLAIM_GAMMA_R3,
       "f is a 3RDF: %s; all %d colours are essential; removable ones: %d"
       % (is_3rdf(adj, f), weight(f), len(survivors)))
    # A single exhaustive pass per fixture feeds both of the next two checks.
    # min_3rdf_weight scans the weights upwards from 0, so the lightest weight
    # it returns already settles whether a 3RDF of weight |V|/2 - 1 exists:
    # such a function exists iff that lightest weight is at most |V|/2 - 1
    # (adding a colour never destroys the defining implication).  The old
    # separate exists_3rdf_of_weight(g, |V|/2 - 1) pass repeated exactly this
    # enumeration, so it is folded in rather than dropped.
    pairs = []
    for name, g in FIXTURES:
        pairs.append((name, len(g), min_3rdf_weight(g),
                      min_dominating_size(product_with_k3(g))))
    bad = [nm for nm, nv, a, _ in pairs if a is None or a <= nv // 2 - 1]
    ck("lemma_n_le_2w_survives_exhaustive_small_cubic_search",
       bool(pairs) and not bad,
       "no 3RDF of weight |V|/2-1 on %s; exhaustive lightest weights %s "
       "against the |V|/2 thresholds %s"
       % (", ".join(nm for nm, _, _, _ in pairs),
          [a for _, _, a, _ in pairs],
          [nv // 2 for _, nv, _, _ in pairs]))
    ck("reduction_behind_the_lower_bound_holds_exactly_on_small_graphs",
       bool(pairs) and all(a == b and a is not None for _, _, a, b in pairs),
       "gamma_r3(G) vs gamma(G box K3): "
       + ", ".join("%s %d/%d" % (nm, a, b) for nm, _, a, b in pairs))


def check_noncayley(n, aut_set):
    """Derived subgroup, the quotient, and a complete census of order-108
    subgroups.  A Cayley graph on 108 vertices needs a regular -- hence
    transitive -- subgroup of order exactly 108 in its automorphism group."""
    ident = tuple(range(n))
    elems = sorted(aut_set)
    comms, resid = residual(aut_set, ident)
    der = generated(comms, ident)
    ck("derived_subgroup_has_order_54",
       len(der) == CLAIM_DERIVED and is_closed(der)
       and len(aut_set) % len(der) == 0,
       "%d distinct commutators generate a closed set of size %d, index %d"
       % (len(comms), len(der), len(aut_set) // len(der)))
    cos = left_cosets(aut_set, der)
    squares_in = all(pmul(g, g) in der for g in elems)
    ck("quotient_by_derived_subgroup_is_C2_x_C2",
       len(cos) == 4 and squares_in
       and all(len(c) == len(der) for c in cos),
       "cosets: %d, sizes %s, every square inside the derived subgroup: %s"
       % (len(cos), sorted(set(len(c) for c in cos)), squares_in))
    try:
        subs = subgroups_of_order(aut_set, len(aut_set) // 2, resid)
    except ValueError:
        subs = []
    # Independent cross-check: an index-2 subgroup is normal, so it must be a
    # union of conjugacy classes.  Recomputed from scratch, not reused above.
    classes = conjugacy_classes(aut_set)
    unions = all(c <= h or not (c & h) for h in subs for c in classes)
    ck("exactly_three_subgroups_of_index_2",
       len(subs) == CLAIM_INDEX2 and len(aut_set) // len(resid) == 4
       and unions and sum(len(c) for c in classes) == len(aut_set),
       "coset census over A/<squares, commutators> finds %d subgroups of "
       "order %d; that quotient has order %d; %d conjugacy classes and every "
       "census subgroup is a union of them: %s"
       % (len(subs), len(aut_set) // 2, len(aut_set) // len(resid),
          len(classes), unions))
    obs = [orbit_sizes(h, n) for h in subs]
    ck("each_index_2_subgroup_has_two_orbits_of_size_54",
       obs and all(o == [54, 54] for o in obs),
       "orbit size lists %s" % (obs,))
    regular = [h for h in subs if orbit_sizes(h, n) == [n]]
    ck("no_regular_subgroup_so_X_is_not_a_Cayley_graph",
       len(subs) == CLAIM_INDEX2 and not regular and len(aut_set) == 2 * n,
       "subgroups of order %d (=|V|): %d, transitive among them: %d"
       % (n, len(subs), len(regular)))
    # control: the same census, run on a graph that IS a Cayley graph, must
    # find a regular subgroup -- otherwise the test above proves nothing.
    m = 12
    ladder = [sorted({(v + 1) % m, (v - 1) % m, (v + m // 2) % m})
              for v in range(m)]
    la = all_automorphisms(ladder)
    lsubs = subgroups_of_order(set(la), len(la) // 2) if len(la) == 2 * m else []
    lreg = [h for h in lsubs if orbit_sizes(h, m) == [m]]
    ck("census_does_detect_regular_subgroups_on_a_known_Cayley_graph",
       len(la) == 2 * m and len(lreg) >= 1,
       "circulant on %d vertices with connection set {+-1, %d}: |Aut|=%d, "
       "regular subgroups found: %d" % (m, m // 2, len(la), len(lreg)))
    ck("X_answers_the_question_non_Cayley_vertex_transitive_3_RDR",
       orbit_sizes(sorted(aut_set), n) == [n] and not regular
       and len(subs) == CLAIM_INDEX2,
       "vertex-transitive: %s, regular subgroups: %d"
       % (orbit_sizes(sorted(aut_set), n) == [n], len(regular)))


if __name__ == "__main__":
    complete = True
    try:
        n, adj, f, auts, aut_set = main()
        check_gamma(n, adj, f, auts)
        check_noncayley(n, aut_set)
    except Exception as exc:          # a corrupted witness must FAIL, not crash
        complete = False
        print("exception: %r" % (exc,))
    ck("verification_ran_to_completion_on_the_stated_witness", complete)
    print("NOT RE-RUN: (i) the identification of X with the entry CVT[108,42] "
          "of the cubic vertex-transitive census -- asserted in the paper's "
          "abstract and in the first sentence of its proof -- and any "
          "minimality over that census, both of which would need the external "
          "census catalogue; (ii) the general lemma that a 3-rainbow "
          "dominating function of weight w yields a dominating set of size w "
          "in G box K_3, equivalently gamma_r3(G) = gamma(G box K_3), which is "
          "cited rather than machine-proved: its conclusion is checked above on "
          "every automorphic image and every one-colour extension of the "
          "exhibited function, and the equality is checked exhaustively on "
          "%d small cubic graphs, but it is not proved for all cubic graphs.  "
          "This program verifies the exhibited graph itself." % len(FIXTURES))
    finish()

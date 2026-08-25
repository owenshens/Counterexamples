#!/usr/bin/env python3
"""
verify.py -- referee verification for
  "The Ivan-Wang Conjecture on Induced N-Saturation for 1 <= n <= 5"
  (settles Conjecture 16 of Ivan-Wang affirmatively for n <= 5).

Standard library only.  Exact integer / bit-mask arithmetic throughout; no
floating point, no randomness, no external data, no tolerance parameter.

=====================================================================
VALUES TAKEN FROM THE PAPER (inputs; every one is compared against an
independently recomputed value, never against itself)
=====================================================================
  * the poset N itself: 4 elements a,b,c,d with a<b, c<b, c<d (a definition);
  * PAPER_SATSTAR      : sat*(n,N) = 2n for 1 <= n <= 5                (Thm 1);
  * PAPER_TABLE        : Table 1 columns |V_n|, |E_n|, tested,
                         independent, maximal for n = 3,4,5;
  * PAPER_LABELED      : 9, 60, 450 labeled minimum families           (Thm 1);
  * PAPER_ORBIT_COUNTS : 2, 4, 7 classes under Gamma_n = S_n x C_2     (Thm 1);
  * PAPER_ORBITS       : all 13 rows of Table 2 -- orbit size and the
                         canonical representative X, verbatim, as words
                         ("125" = {1,2,5}); the exhibited object;
  * PAPER_EDGE_FORMULA : |E_n| = 8^n - 3*7^n + 3*6^n - 5^n           (Lemma 4);
  * PAPER_PATTERNS,
    PAPER_REQUIRED     : the eight coordinate patterns and the three required
                         ones used in the proof of Lemma 4.

=====================================================================
WHAT THIS PROGRAM DERIVES (the checks)
=====================================================================
  1. N is a poset, its comparability graph is P4, every element has an
     incomparable partner, and N is self-dual (the facts Lemma 2 and the
     Gamma_n action rest on).
  2. The N-inducing 4-subsets of B_n are generated from the DEFINITION (all
     24 bijections onto N, induced sub-poset), for n = 1..6; the closed form
     of Lemma 4 is checked against that count; the P4 criterion of Lemma 2 is
     checked to give the identical set; no copy contains {} or [n] (Lemma 2,
     the computational content of the reduction Lemma 3); and each copy is
     shown to admit exactly ONE role assignment (A,B,C,D), which is what lets
     Lemma 4 count 4-subsets by counting role-assigned quadruples.
  2b. The proof mechanism of Lemma 4 itself: the eight listed coordinate
     patterns are exactly those compatible with A<=B, C<=B, C<=D, and over
     every one of the 8^n such quadruples (n <= 5) "induced copy in the roles
     (A,B,C,D)" is equivalent to "all three of 0001, 1100, 0111 occur".
  3. sat*(n,N) for n = 1,2,3,4 by exhaustive brute force over ALL 2^(2^n)
     subfamilies of B_n, using only the DEFINITION of induced-N-saturation --
     no lemma, no reduction, no hypergraph.  This yields the minimum size and
     the number of minimum families independently of the paper's method, and
     the reduction Lemma 3 is then verified as a set equality for n = 3,4.
  4. The complete census of Table 1 for n = 3,4,5, twice and independently:
     (A) an unpruned exhaustive scan of every one of the 57 / 6,476 /
         8,656,937 candidates X with |X| <= 2n-2, each tested by explicit
         4-subset membership; and
     (B) a pruned depth-first enumeration with incremental closure masks.
     The two must agree on the candidate count, the independent count by
     size, and the identical SET of maximal independent sets.
  5. Each of the 9 / 60 / 450 minimum families F = X u {{},[n]} is verified
     induced-N-saturated directly from the definition, and |F| = 2n.
  6. Gamma_n is built as 2*n! explicit permutations of V_n and VERIFIED to be
     a permutation group (every element a bijection of V_n, identity present,
     closed under composition, inverses present) -- without which the orbit
     language and the orbit-stabilizer identity would not be earned.  The
     census set is shown closed under it, orbits and canonical forms
     (lexicographically least image of the mu-sorted family) are computed, and
     the canonical-form classes are confronted with the orbits recomputed
     independently by closing the census under the action, so that "these
     classes are the Gamma_n-orbits" is computed and not true by construction.
     The resulting orbit count, orbit sizes, orbit-stabilizer identity and the
     13 canonical representatives are then compared with Table 2.  The paper's
     claim that Gamma_n preserves induced-N-saturation is tested on the
     complete brute-forced set of saturated families of every size (9 at
     n = 3, 118 at n = 4), not only on the minimum ones.

SCOPE.  The census of Table 1 is reproduced IN FULL: all 57 / 6,476 /
8,656,937 candidates for n = 3,4,5, nothing sampled or truncated.  The only
computation the paper's argument uses that is NOT redone here from the raw
definition is the n = 5 minimality: scanning all 2^(2^5) = 2^32 subfamilies of
B_5 is out of reach, so for n = 5 the step from "smallest maximal independent
set has 2n-2 elements" to "sat*(5,N) = 10" rests on the reduction of Lemma 3,
whose computational content is checked at n = 5 and whose full set equality is
checked at n = 3,4.  NOT RE-RUN: the paper's own first implementation, which is
not part of this bundle and is not named here; this file is an independent
reimplementation written from the paper's specification alone.  At n = 6 only
the 4-subset facts are checked (the Lemma 4 edge count, the Lemma 2 P4
criterion, role uniqueness, no extremes); no census, minimality or orbit claim
is checked at n = 6 and nothing at all for n >= 7.  The paper itself claims
nothing for n >= 6.

KNOWN-VACUOUS BOOKKEEPING (recorded so that no reader mistakes it for
evidence).  In orbits_partition_the_labeled_census the conjuncts
"len(union) == tot" and "union == census" hold by construction, because orb is
built by bucketing the very list that census is built from; the live content of
that check is the labeled total, the absence of repeats in the census list, and
check_orbits_are_orbits.  In table1_tested_candidates the equality
"tested == sum_k C(2^n-2,k)" is an arithmetic identity about
itertools.combinations rather than two independent routes to the number; its
real content is that the enumeration loop was not truncated, and both sides are
separately compared with the paper's value.  In comp_graph_is_P4 the "3 edges"
test is implied by the degree sequence [1,1,2,2], and in is_copy_in_roles the
"x != y" test is implied by the three required incomparabilities; both are kept
as redundant guards.

Runtime: NOT MEASURED -- no timing figure is claimed here, and any figure
quoted elsewhere for this file should be treated as unverified.  The run is
single process and finishes without intervention, but it is long: the dominant
costs are the n=5 census, computed twice over all 8,656,937 candidates with
1,690,127 independent-set nodes, and the n=6 pass over all 635,376 4-subsets of
B_6 with up to 24 role bijections each.  Under a bare CPython interpreter
expect minutes to tens of minutes rather than seconds; the program prints each
check as it completes and flushes, so progress is visible throughout and a
long silence at a stage is expected, not a hang.
Exit status is 0 if and only if every check passes.
"""

import sys
from itertools import combinations, permutations

_RESULTS = []


def check(name, ok, detail=""):
    _RESULTS.append(bool(ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    sys.stdout.flush()
    return bool(ok)


def info(msg):
    print("# " + str(msg))
    sys.stdout.flush()


def nCk(a, b):
    if b < 0 or b > a:
        return 0
    r = 1
    for i in range(b):
        r = r * (a - i) // (i + 1)
    return r


# =====================  VALUES TAKEN FROM THE PAPER  =====================
# Theorem 1: sat*(n,N) = 2n for 1 <= n <= 5.
PAPER_SATSTAR = {n: 2 * n for n in range(1, 6)}
# Theorem 1 / Table 1 last column: labeled minimum families.
PAPER_LABELED = {3: 9, 4: 60, 5: 450}
# Theorem 1 / Proposition 5: number of Gamma_n classes.
PAPER_ORBIT_COUNTS = {3: 2, 4: 4, 5: 7}
# Table 1: n -> (|V_n|, |E_n|, tested, independent, maximal)
PAPER_TABLE = {
    3: (6, 6, 57, 51, 9),
    4: (14, 156, 6476, 3110, 60),
    5: (30, 2550, 8656937, 1690127, 450),
}
# Lemma 4 closed form for |E_n|.
def PAPER_EDGE_FORMULA(n):
    return 8 ** n - 3 * 7 ** n + 3 * 6 ** n - 5 ** n


# Table 2, verbatim: (n, orbit size, canonical representative X as words).
PAPER_ORBITS = [
    (3, 6, "1,2,12,3"),
    (3, 3, "1,2,13,23"),
    (4, 24, "1,2,12,3,123,4"),
    (4, 6, "1,2,12,3,4,34"),
    (4, 24, "1,2,12,3,124,34"),
    (4, 6, "1,2,12,123,124,34"),
    (5, 120, "1,2,12,3,123,4,1234,5"),
    (5, 60, "1,2,12,3,123,4,5,45"),
    (5, 120, "1,2,12,3,123,4,1235,45"),
    (5, 60, "1,2,12,3,123,1234,1235,45"),
    (5, 30, "1,2,12,3,4,34,1234,5"),
    (5, 30, "1,2,12,3,4,34,125,345"),
    (5, 30, "1,2,12,123,1234,1235,45,345"),
]

# The poset N as a definition: elements 0,1,2,3 = a,b,c,d with
# a < b, c < b, c < d and no other relations.
N_LESS = [(0, 1), (2, 1), (2, 3)]
N_INCOMP = [(0, 2), (0, 3), (1, 3)]

def n_relation_matrix():
    """Strict order matrix of N built from N_LESS, closed under nothing."""
    R = [[False] * 4 for _ in range(4)]
    for a, b in N_LESS:
        R[a][b] = True
    return R


def check_N_structure():
    R = n_relation_matrix()
    # antisymmetric, irreflexive, transitive
    ok = True
    for i in range(4):
        if R[i][i]:
            ok = False
    for i in range(4):
        for j in range(4):
            if R[i][j] and R[j][i]:
                ok = False
            for k in range(4):
                if R[i][j] and R[j][k] and not R[i][k]:
                    ok = False
    check("N_is_a_strict_poset", ok)
    # the declared incomparable pairs really are incomparable, and
    # N_LESS u N_INCOMP accounts for all 6 unordered pairs exactly once
    pairs = set()
    for a, b in N_LESS:
        pairs.add((min(a, b), max(a, b)))
    for a, b in N_INCOMP:
        pairs.add((min(a, b), max(a, b)))
    inc_ok = all(not R[a][b] and not R[b][a] for a, b in N_INCOMP)
    check("N_pair_partition_complete", inc_ok and len(pairs) == 6,
          "pairs=%d" % len(pairs))
    # comparability graph of N is P4 (Lemma 2)
    adj = [[R[i][j] or R[j][i] for j in range(4)] for i in range(4)]
    deg = sorted(sum(1 for j in range(4) if adj[i][j]) for i in range(4))
    nedge = sum(1 for i in range(4) for j in range(i + 1, 4) if adj[i][j])
    check("N_comparability_graph_is_P4",
          deg == [1, 1, 2, 2] and nedge == 3 and graph_connected(adj, 4),
          "deg=%s edges=%d" % (deg, nedge))
    # every element of N is incomparable with some other element (Lemma 2 pf)
    lonely = [i for i in range(4)
              if all(adj[i][j] for j in range(4) if j != i)]
    check("N_every_element_has_an_incomparable_partner", not lonely,
          "violators=%s" % lonely)
    # N is self-dual (justifies complementation preserving saturation)
    dual_iso = False
    for p in permutations(range(4)):
        if all(R[j][i] == R[p[i]][p[j]] for i in range(4) for j in range(4)):
            dual_iso = True
            break
    check("N_is_self_dual", dual_iso)


def graph_connected(adj, k):
    seen = [False] * k
    seen[0] = True
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(k):
            if adj[u][v] and not seen[v]:
                seen[v] = True
                stack.append(v)
    return all(seen)


_P4PERMS = list(permutations(range(4)))


def is_copy_in_roles(r):
    """True iff the subsets r[0..3], taken in the roles (a,b,c,d) of N, form an
    INDUCED copy: exactly the relations a<b, c<b, c<d hold and the other three
    pairs are incomparable.  (Equal members are rejected automatically: every
    one of the six pairs is constrained.)"""
    for a, b in N_LESS:
        x, y = r[a], r[b]
        if x & y != x or x == y:      # need x a proper subset of y
            return False
    for a, b in N_INCOMP:
        x, y = r[a], r[b]
        if x & y == x or x & y == y:  # comparable -> not induced
            return False
    return True


def role_count(s):
    """How many of the 24 bijections onto (a,b,c,d) exhibit s as a copy of N."""
    return sum(1 for p in _P4PERMS
               if is_copy_in_roles([s[p[0]], s[p[1]], s[p[2]], s[p[3]]]))


def induces_N(s):
    """True iff the four subsets s[0..3] span an induced copy of N."""
    for p in _P4PERMS:
        if is_copy_in_roles([s[p[0]], s[p[1]], s[p[2]], s[p[3]]]):
            return True
    return False


def comp_graph_is_P4(s):
    """True iff the comparability graph of the four subsets is a path P4."""
    adj = [[False] * 4 for _ in range(4)]
    ne = 0
    for i in range(4):
        for j in range(i + 1, 4):
            x, y = s[i], s[j]
            if x & y == x or x & y == y:
                adj[i][j] = adj[j][i] = True
                ne += 1
    if ne != 3:
        return False
    deg = sorted(sum(1 for j in range(4) if adj[i][j]) for i in range(4))
    return deg == [1, 1, 2, 2] and graph_connected(adj, 4)


def scan_quads(n, both=True):
    """Generate the induced copies of N among the 4-subsets of B_n = 2^[n],
    straight from the definition.  If both, also collect the 4-subsets whose
    comparability graph is P4, to test the criterion of Lemma 2."""
    top = 1 << n
    by_def = []
    by_p4 = [] if both else None
    for q in combinations(range(top), 4):
        s = [q[0], q[1], q[2], q[3]]
        if induces_N(s):
            by_def.append(q)
        if both and comp_graph_is_P4(s):
            by_p4.append(q)
    return by_def, by_p4


def check_copies(n, both=True):
    """Checks on the generated copies; returns the definitional copy list.
    A check is only reported where it is not vacuous: the P4 criterion needs
    at least one 4-subset (n >= 2) and the extremes claim needs at least one
    copy (n >= 3)."""
    by_def, by_p4 = scan_quads(n, both)
    if both and nCk(1 << n, 4) > 0:
        check("lemma2_induced_N_iff_comparability_P4_n%d" % n,
              set(by_def) == set(by_p4),
              "byN=%d byP4=%d symdiff=%d"
              % (len(by_def), len(by_p4),
                 len(set(by_def) ^ set(by_p4))))
    if by_def:
        # Lemma 4 counts role-assigned quadruples (A,B,C,D); that equals the
        # number of 4-subsets only because the roles are unique.
        multi = [q for q in by_def
                 if role_count([q[0], q[1], q[2], q[3]]) != 1]
        check("each_copy_has_a_unique_role_assignment_n%d" % n, not multi,
              "copies with a role count != 1: %d of %d"
              % (len(multi), len(by_def)))
    if both and by_def:
        top = (1 << n) - 1
        bad = [q for q in by_def if 0 in q or top in q]
        check("lemma2_no_copy_contains_empty_or_full_n%d" % n, not bad,
              "violators=%d of %d copies" % (len(bad), len(by_def)))
    check("lemma4_edge_count_closed_form_n%d" % n,
          len(by_def) == PAPER_EDGE_FORMULA(n),
          "counted=%d formula=%d" % (len(by_def), PAPER_EDGE_FORMULA(n)))
    return by_def


PAPER_PATTERNS = [(0, 0, 0, 0), (0, 0, 0, 1), (0, 1, 0, 0), (0, 1, 0, 1),
                  (1, 1, 0, 0), (1, 1, 0, 1), (0, 1, 1, 1), (1, 1, 1, 1)]
PAPER_REQUIRED = [(0, 0, 0, 1), (1, 1, 0, 0), (0, 1, 1, 1)]


def check_lemma4_proof(nmax):
    """The mechanism of Lemma 4's proof, verified rather than assumed:
    (i) the eight listed coordinate patterns are exactly those compatible with
    the three prescribed inclusions A<=B, C<=B, C<=D; and (ii) for a quadruple
    with all patterns among the eight, being an induced copy of N in the roles
    (A,B,C,D) is equivalent to the three patterns 0001, 1100, 0111 all
    occurring -- checked over every one of the 8^n such quadruples."""
    allowed = [p for p in _BITS4
               if p[0] <= p[1] and p[2] <= p[1] and p[2] <= p[3]]
    check("lemma4_eight_patterns_are_exactly_the_compatible_ones",
          set(allowed) == set(PAPER_PATTERNS) and len(PAPER_PATTERNS) == 8,
          "computed=%d listed=%d" % (len(allowed), len(PAPER_PATTERNS)))
    req = set(PAPER_REQUIRED)
    for n in range(1, nmax + 1):
        bad = 0
        for word in _words(PAPER_PATTERNS, n):
            r = [0, 0, 0, 0]
            for i, p in enumerate(word):
                for j in range(4):
                    if p[j]:
                        r[j] |= 1 << i
            if is_copy_in_roles(r) != req.issubset(set(word)):
                bad += 1
        check("lemma4_induced_iff_the_three_patterns_occur_n%d" % n, bad == 0,
              "quadruples=%d mismatches=%d" % (8 ** n, bad))


_BITS4 = [(a, b, c, d) for a in (0, 1) for b in (0, 1)
          for c in (0, 1) for d in (0, 1)]


def _words(alphabet, n):
    if n == 0:
        yield ()
        return
    for w in _words(alphabet, n - 1):
        for a in alphabet:
            yield w + (a,)


def hypergraph(n, copies):
    """Recast the copies as a 4-uniform hypergraph on V_n = B_n - {{},[n]},
    vertex i  <->  the subset with mu = i+1."""
    top = (1 << n) - 1
    V = [s for s in range(1 << n) if s != 0 and s != top]
    m = len(V)
    if V != list(range(1, top)):
        raise AssertionError("V_n indexing broken")
    edges = []
    for q in copies:
        e = tuple(sorted(x - 1 for x in q))
        if e[0] < 0 or e[3] >= m:
            raise AssertionError("copy meets an extreme set")
        edges.append(e)
    return m, edges


def copy_masks(n, copies):
    """copies as bit masks over the 2^n subsets, plus, for each subset v, the
    masks of the copies through v with v deleted."""
    top = 1 << n
    masks = [sum(1 << x for x in q) for q in copies]
    red = [[] for _ in range(top)]
    for q, mk in zip(copies, masks):
        for v in q:
            red[v].append(mk ^ (1 << v))
    return masks, red


def family_is_N_free(F, masks):
    for mk in masks:
        if mk & F == mk:
            return False
    return True


def family_is_saturated(F, top, masks, red):
    """The definition: F contains no induced copy of N, and F u {S} contains
    one for every S in B_n - F."""
    if not family_is_N_free(F, masks):
        return False
    for v in range(top):
        if F >> v & 1:
            continue
        hit = False
        for r in red[v]:
            if r & F == r:
                hit = True
                break
        if not hit:
            return False
    return True


def family_bruteforce(n, copies):
    """EVERY one of the 2^(2^n) subfamilies of B_n, tested from the definition
    alone.  Returns (all saturated family masks, min size, count at min)."""
    top = 1 << n
    masks, red = copy_masks(n, copies)
    sat = []
    for F in range(1 << top):
        if family_is_saturated(F, top, masks, red):
            sat.append(F)
    sizes = [bin(F).count("1") for F in sat]
    best = min(sizes)
    return sat, best, sizes.count(best)


def triple_masks(edges):
    """M_n(T) of the paper: for a triple T, the mask of v with T u {v} an edge."""
    TM = {}
    for (a, b, c, d) in edges:
        for t, v in (((a, b, c), d), ((a, b, d), c),
                     ((a, c, d), b), ((b, c, d), a)):
            TM[t] = TM.get(t, 0) | (1 << v)
    return TM


def census_scan(m, edges, kmax):
    """(A) UNPRUNED exhaustive scan: literally every candidate X subset V_n
    with |X| <= kmax is generated and tested.  Independence is decided by
    explicit membership of the C(|X|,4) 4-subsets in the edge set; maximality
    by C_n(X) = union of M_n(T) over triples T of X, as in the paper."""
    ES = set(edges)
    TM = triple_masks(edges)
    FULL = (1 << m) - 1
    tested = 0
    ind = [0] * (kmax + 1)
    mxc = [0] * (kmax + 1)
    maximal = []
    for k in range(kmax + 1):
        for c in combinations(range(m), k):
            tested += 1
            indep = True
            for q in combinations(c, 4):
                if q in ES:
                    indep = False
                    break
            if not indep:
                continue
            ind[k] += 1
            cov = 0
            for t in combinations(c, 3):
                cov |= TM.get(t, 0)
            if (cov | sum(1 << v for v in c)) == FULL:
                mxc[k] += 1
                maximal.append(c)
    return tested, ind, mxc, maximal


def census_dfs(m, edges, kmax):
    """(B) PRUNED depth-first enumeration of the independent sets of size
    <= kmax (independence is downward closed, so every independent set is
    reached through independent prefixes).  The blocked mask is carried
    incrementally: C(X u {v}) = C(X) u union over pairs P of X of M(P u {v}).
    Completely different bookkeeping from census_scan."""
    TM = triple_masks(edges)
    FULL = (1 << m) - 1
    PT = {}
    for a in range(m):
        for b in range(a + 1, m):
            PT[a * m + b] = [TM.get((a, b, v), 0) for v in range(m)]
    ind = [0] * (kmax + 1)
    mxc = [0] * (kmax + 1)
    maximal = []
    sys.setrecursionlimit(10000)

    def rec(X, xs, blocked, lo):
        k = len(xs)
        ind[k] += 1
        if (blocked | X) == FULL:
            mxc[k] += 1
            maximal.append(tuple(xs))
        if k == kmax:
            return
        free = ~(blocked | X) & FULL & ~((1 << lo) - 1)
        rows = [PT[xs[i] * m + xs[j]]
                for i in range(k) for j in range(i + 1, k)]
        while free:
            low = free & -free
            v = low.bit_length() - 1
            free ^= low
            nb = blocked
            for row in rows:
                nb |= row[v]
            xs.append(v)
            rec(X | low, xs, nb, v + 1)
            xs.pop()

    rec(0, [], 0, 0)
    return ind, mxc, maximal


def group_elements(n, m):
    """Gamma_n = S_n x C_2 acting on V_n, as explicit vertex permutations.
    Vertex i is the subset with mu(S) = i+1; the nontrivial element of C_2
    complements."""
    full = (1 << n) - 1
    G = []
    for p in permutations(range(n)):
        for comp in (0, 1):
            img = []
            for i in range(m):
                s = i + 1
                t = 0
                for b in range(n):
                    if s >> b & 1:
                        t |= 1 << p[b]
                if comp:
                    t = full ^ t
                if not 1 <= t <= full - 1:
                    raise AssertionError("group left V_n")
                img.append(t - 1)
            G.append(tuple(img))
    return G


def canonical(X, G):
    """Lexicographically least image of the mu-sorted family under Gamma_n."""
    best = None
    for g in G:
        t = tuple(sorted(g[v] + 1 for v in X))
        if best is None or t < best:
            best = t
    return best


def word_of_mu(mu, n):
    return "".join(str(i + 1) for i in range(n) if mu >> i & 1)


def family_words(mus, n):
    return ",".join(word_of_mu(mu, n) for mu in mus)


def parse_family(text, n):
    """Decode a Table 2 representative.  Returns (mu tuple as printed, notes)
    where notes lists every well-formedness violation found."""
    notes = []
    mus = []
    for w in text.split(","):
        if not w or not all(ch.isdigit() for ch in w):
            notes.append("bad word %r" % w)
            continue
        ds = [int(ch) for ch in w]
        if any(d < 1 or d > n for d in ds):
            notes.append("element out of [n] in %r" % w)
        if ds != sorted(set(ds)):
            notes.append("word %r not a strictly increasing set" % w)
        mu = 0
        for d in ds:
            mu |= 1 << (d - 1)
        if mu == 0 or mu == (1 << n) - 1:
            notes.append("word %r is an extreme set, not in V_n" % w)
        mus.append(mu)
    if len(set(mus)) != len(mus):
        notes.append("repeated member")
    if list(mus) != sorted(mus):
        notes.append("members not printed in increasing mu order")
    return tuple(mus), notes


def check_minimum_families(n, maximal, copies):
    """Every minimum family F = X u {{},[n]} is re-tested against the raw
    DEFINITION of induced-N-saturation (no hypergraph, no lemma), and its size
    against the claim sat*(n,N) = 2n."""
    top = 1 << n
    masks, red = copy_masks(n, copies)
    base = 1 | (1 << (top - 1))
    bad_sat = 0
    bad_size = 0
    for X in maximal:
        F = base
        for v in X:
            F |= 1 << (v + 1)
        if bin(F).count("1") != PAPER_SATSTAR[n]:
            bad_size += 1
        if not family_is_saturated(F, top, masks, red):
            bad_sat += 1
    check("minimum_families_saturated_by_definition_n%d" % n,
          bad_sat == 0 and len(maximal) == PAPER_LABELED[n],
          "families=%d expected=%d failures=%d"
          % (len(maximal), PAPER_LABELED[n], bad_sat))
    check("minimum_family_size_is_2n_n%d" % n,
          bad_size == 0 and len(maximal) > 0,
          "wrong size=%d of %d" % (bad_size, len(maximal)))


def compute_orbits(maximal, G):
    """Group the labeled minimum families by canonical form."""
    orb = {}
    for X in maximal:
        orb.setdefault(canonical(X, G), set()).add(tuple(sorted(X)))
    return orb


def check_group_action(n, m, maximal, G):
    check("group_order_is_2_times_n_factorial_n%d" % n,
          len(G) == 2 * fact(n) and len(set(G)) == 2 * fact(n),
          "elements=%d distinct=%d expected=%d"
          % (len(G), len(set(G)), 2 * fact(n)))
    S = set(tuple(sorted(X)) for X in maximal)
    bad = 0
    for X in S:
        for g in G:
            if tuple(sorted(g[v] for v in X)) not in S:
                bad += 1
                break
    check("census_closed_under_Gamma_n%d" % n, bad == 0 and len(S) > 0,
          "families=%d with an image outside the census=%d" % (len(S), bad))


def fact(k):
    r = 1
    for i in range(2, k + 1):
        r *= i
    return r


def check_table2(n, maximal, G, orb):
    """The exhibited object: decode Table 2, print it back, and confront it
    with the orbits computed from the census."""
    rows = [(sz, w) for (nn, sz, w) in PAPER_ORBITS if nn == n]
    census = set(tuple(sorted(X)) for X in maximal)
    notes = []
    parsed = []
    for sz, w in rows:
        mus, nt = parse_family(w, n)
        notes += ["%s: %s" % (w, x) for x in nt]
        if len(mus) != 2 * n - 2:
            notes.append("%s: has %d members, expected %d"
                         % (w, len(mus), 2 * n - 2))
        parsed.append((sz, mus))
        info("  Table 2 row n=%d size=%-3d X = {%s}  mu = %s"
             % (n, sz, w, list(mus)))
    check("table2_representatives_wellformed_n%d" % n,
          not notes and len(rows) > 0,
          "; ".join(notes[:4]) if notes else "%d rows" % len(rows))
    outside = [family_words(mus, n) for _, mus in parsed
               if tuple(sorted(mu - 1 for mu in mus)) not in census]
    check("table2_representatives_are_maximal_independent_n%d" % n,
          not outside and len(parsed) > 0,
          "%d rows, not in census: %s" % (len(parsed), outside[:3]))
    lim = (1 << n) - 2
    usable = [mus for _, mus in parsed
              if mus and all(1 <= mu <= lim for mu in mus)]
    noncanon = [family_words(mus, n) for mus in usable
                if canonical([mu - 1 for mu in mus], G) != mus]
    noncanon += ["<not a family in V_%d>" % n] * (len(parsed) - len(usable))
    check("table2_representatives_are_canonical_n%d" % n,
          not noncanon and len(parsed) > 0,
          "%d rows, non-canonical: %s" % (len(parsed), noncanon[:3]))
    got = sorted((len(v), k) for k, v in orb.items())
    want = sorted((sz, mus) for sz, mus in parsed)
    check("table2_orbit_rows_match_computed_orbits_n%d" % n, got == want,
          "computed=%s" % ([(s, family_words(k, n)) for s, k in got],))
    check("orbit_count_n%d" % n, len(orb) == PAPER_ORBIT_COUNTS[n],
          "computed=%d paper=%d" % (len(orb), PAPER_ORBIT_COUNTS[n]))
    tot = sum(len(v) for v in orb.values())
    union = set()
    for v in orb.values():
        union |= v
    # NOTE on vacuity: orb was built by bucketing THIS list under canonical(),
    # and census is the set of the same list, so "len(union) == tot" and
    # "union == census" hold by construction and cannot fail -- they are
    # bookkeeping, not evidence.  The live conjuncts are the labeled total and
    # "len(maximal) == tot", which fails if the census list has a repeat.  The
    # substantive claim (these buckets really are the Gamma_n-orbits) is
    # check_orbits_are_orbits below, which recomputes them by group closure.
    check("orbits_partition_the_labeled_census_n%d" % n,
          tot == PAPER_LABELED[n] and len(union) == tot and union == census
          and len(maximal) == tot,
          "sum_of_orbit_sizes=%d distinct=%d census=%d listed=%d paper=%d"
          % (tot, len(union), len(census), len(maximal), PAPER_LABELED[n]))
    bad = []
    for k, v in orb.items():
        X = [mu - 1 for mu in k]
        stab = sum(1 for g in G if tuple(sorted(g[x] for x in X))
                   == tuple(sorted(X)))
        if len(v) * stab != len(G):
            bad.append((family_words(k, n), len(v), stab))
    check("orbit_stabilizer_identity_n%d" % n, not bad and len(orb) > 0,
          "%d orbits, violations=%s" % (len(orb), bad))


def check_census(n, m, edges, scanres, dfsres):
    tested, ind, mxc, maximal = scanres
    ind2, mxc2, maximal2 = dfsres
    pV, pE, pT, pI, pM = PAPER_TABLE[n]
    kmax = 2 * n - 2
    check("table1_V_size_n%d" % n, m == pV and m == (1 << n) - 2,
          "computed=%d paper=%d" % (m, pV))
    check("table1_edge_count_n%d" % n, len(edges) == pE,
          "computed=%d paper=%d" % (len(edges), pE))
    binsum = sum(nCk(m, k) for k in range(kmax + 1))
    check("table1_tested_candidates_n%d" % n,
          tested == pT and binsum == pT,
          "enumerated=%d binomial_sum=%d paper=%d" % (tested, binsum, pT))
    check("table1_independent_count_n%d" % n, sum(ind) == pI,
          "computed=%d paper=%d by_size=%s" % (sum(ind), pI, ind))
    check("table1_maximal_count_n%d" % n, sum(mxc) == pM,
          "computed=%d paper=%d" % (sum(mxc), pM))
    check("all_maximal_have_size_exactly_2n_minus_2_n%d" % n,
          sum(mxc[:kmax]) == 0 and mxc[kmax] == pM,
          "by_size=%s" % mxc)
    check("two_independent_census_algorithms_agree_n%d" % n,
          ind == ind2 and mxc == mxc2
          and set(tuple(sorted(x)) for x in maximal)
          == set(tuple(sorted(x)) for x in maximal2),
          "scan_ind=%d dfs_ind=%d scan_max=%d dfs_max=%d"
          % (sum(ind), sum(ind2), sum(mxc), sum(mxc2)))
    # sat*(n,N): smallest maximal independent set, plus the two extreme sets
    # that Lemma 2 forces into every saturated family.
    sizes = [k for k in range(kmax + 1) if mxc[k]]
    derived = (min(sizes) + 2) if sizes else None
    check("satstar_equals_2n_from_census_n%d" % n,
          derived == PAPER_SATSTAR[n],
          "derived=%s paper=%d" % (derived, PAPER_SATSTAR[n]))
    return maximal


def check_reduction(n, m, edges, sat_bf):
    """Lemma 3 as a set equality, for the n where brute force is possible:
    the saturated families are exactly {{},[n]} u X for X maximal independent
    (over ALL sizes, not merely the minimum ones)."""
    _, _, allmax = census_dfs(m, edges, m)
    top = 1 << n
    base = 1 | (1 << (top - 1))
    fromhyp = set()
    for X in allmax:
        F = base
        for v in X:
            F |= 1 << (v + 1)
        fromhyp.add(F)
    check("reduction_lemma3_set_equality_n%d" % n, fromhyp == set(sat_bf),
          "bruteforce=%d from_hypergraph=%d symdiff=%d"
          % (len(sat_bf), len(fromhyp), len(fromhyp ^ set(sat_bf))))


def check_gamma_is_a_group(n, m, G):
    """Gamma_n is claimed to be the GROUP S_n x C_2 acting on V_n.  Everything
    downstream (canonical form as a complete orbit invariant, the word
    "orbit", the orbit-stabilizer identity) is valid only for a group of
    BIJECTIONS, and group_order_is_2_times_n_factorial checks merely that the
    2*n! image tuples are distinct -- not that any one of them is a
    permutation, nor that the set is closed.  Verified here rather than
    assumed: every element is a bijection of V_n, the identity is present, the
    set is closed under composition, and every element's inverse is in it."""
    idl = list(range(m))
    ident = tuple(idl)
    nonbij = [g for g in G if sorted(g) != idl]
    Gs = set(G)
    has_id = ident in Gs
    closed = True
    for g in G:
        for h in G:
            if tuple(map(g.__getitem__, h)) not in Gs:
                closed = False
                break
        if not closed:
            break
    inv = True
    for g in G:
        ig = [0] * m
        for i, x in enumerate(g):
            ig[x] = i
        if tuple(ig) not in Gs:
            inv = False
            break
    check("gamma_is_a_permutation_group_n%d" % n,
          not nonbij and has_id and closed and inv and len(G) == 2 * fact(n),
          "elements=%d non_bijections=%d identity=%s closed=%s inverses=%s"
          % (len(G), len(nonbij), has_id, closed, inv))


def gamma_orbits_by_closure(maximal, G):
    """The orbits of the census under the ACTION, computed with no reference to
    canonical forms: repeatedly apply every group element until closure."""
    S = set(tuple(sorted(X)) for X in maximal)
    seen = set()
    orbits = []
    for X in sorted(S):
        if X in seen:
            continue
        comp = {X}
        stack = [X]
        while stack:
            Y = stack.pop()
            for g in G:
                Z = tuple(sorted(g[v] for v in Y))
                if Z not in comp:
                    comp.add(Z)
                    stack.append(Z)
        seen |= comp
        orbits.append(comp)
    return orbits


def check_orbits_are_orbits(n, maximal, G, orb):
    """compute_orbits DEFINES its classes as the fibres of canonical(), so
    "these classes are the Gamma_n-orbits" -- the content of the paper's
    orbit proposition -- would otherwise be true by construction rather than
    computed.  Here the orbits are recomputed by closing the census under the
    action and the two partitions are confronted."""
    orbits = gamma_orbits_by_closure(maximal, G)
    a = sorted(sorted(o) for o in orbits)
    b = sorted(sorted(v) for v in orb.values())
    check("canonical_classes_are_the_true_Gamma_orbits_n%d" % n,
          a == b and len(orbits) > 0,
          "closure_orbits=%d sizes=%s canonical_classes=%d sizes=%s"
          % (len(orbits), sorted(len(o) for o in orbits),
             len(orb), sorted(len(v) for v in orb.values())))


def gamma_full_elements(n):
    """Gamma_n as explicit permutations of ALL of B_n (V_n plus the two extreme
    sets), needed to test the paper's claim about saturated families rather
    than only about the vertex set of the hypergraph."""
    top = 1 << n
    full = top - 1
    G = []
    for p in permutations(range(n)):
        for comp in (0, 1):
            img = []
            for s in range(top):
                t = 0
                for b in range(n):
                    if s >> b & 1:
                        t |= 1 << p[b]
                if comp:
                    t = full ^ t
                img.append(t)
            if sorted(img) != list(range(top)):
                raise AssertionError("Gamma_n element is not a permutation"
                                     " of B_n")
            G.append(tuple(img))
    if len(set(G)) != 2 * fact(n):
        raise AssertionError("Gamma_n action on B_n is not faithful")
    return G


def check_saturation_closed_under_gamma(n, sat_bf):
    """The paper's claim that Gamma_n -- in particular complementation, N being
    self-dual -- PRESERVES induced-N-saturation.  census_closed_under_Gamma
    tests this only on the minimum families; here it is tested on the complete
    brute-forced set of saturated families of ALL sizes (9 at n=3, 118 at
    n=4), which is the actual statement."""
    top = 1 << n
    GF = gamma_full_elements(n)
    S = set(sat_bf)
    bad = 0
    for F in S:
        for g in GF:
            H = 0
            for s in range(top):
                if F >> s & 1:
                    H |= 1 << g[s]
            if H not in S:
                bad += 1
                break
    check("all_saturated_families_closed_under_Gamma_n%d" % n,
          bad == 0 and len(S) > 0,
          "saturated_families=%d group_order=%d with an image outside=%d"
          % (len(S), len(GF), bad))


def main():
    info("verify.py -- sat*(n,N) = 2n for 1 <= n <= 5 (Ivan-Wang Conj. 16)")
    info("stdlib only, exact integer arithmetic, no randomness")
    info("")
    info("--- the poset N ---")
    check_N_structure()

    info("--- the coordinate-pattern mechanism of Lemma 4 ---")
    check_lemma4_proof(5)

    info("--- induced copies of N in B_n, generated from the definition ---")
    copies = {}
    for n in range(1, 7):
        copies[n] = check_copies(n, both=True)
        info("  n=%d: %d induced copies of N among the %d 4-subsets of B_n"
             % (n, len(copies[n]), nCk(1 << n, 4)))

    info("--- sat*(n,N) by brute force over ALL 2^(2^n) subfamilies of B_n,"
         " from the definition alone ---")
    satbf = {}
    for n in (1, 2, 3, 4):
        sat, best, cnt = family_bruteforce(n, copies[n])
        satbf[n] = sat
        info("  n=%d: %d saturated families in all, minimum size %d, %d of"
             " minimum size" % (n, len(sat), best, cnt))
        check("bruteforce_satstar_equals_2n_n%d" % n,
              best == PAPER_SATSTAR[n],
              "computed=%d paper=%d" % (best, PAPER_SATSTAR[n]))
        if n in PAPER_LABELED:
            check("bruteforce_labeled_minimum_count_n%d" % n,
                  cnt == PAPER_LABELED[n],
                  "computed=%d paper=%d" % (cnt, PAPER_LABELED[n]))
        else:
            check("bruteforce_only_saturated_family_is_B_n_n%d" % n,
                  sat == [(1 << (1 << n)) - 1],
                  "families=%d" % len(sat))

    info("--- the census of Table 1, computed twice, independently ---")
    for n in (3, 4, 5):
        m, edges = hypergraph(n, copies[n])
        scanres = census_scan(m, edges, 2 * n - 2)
        dfsres = census_dfs(m, edges, 2 * n - 2)
        info("  n=%d: tested %d, independent %d, maximal %d"
             % (n, scanres[0], sum(scanres[1]), sum(scanres[2])))
        maximal = check_census(n, m, edges, scanres, dfsres)
        if n in satbf:
            check_reduction(n, m, edges, satbf[n])
        check_minimum_families(n, maximal, copies[n])
        G = group_elements(n, m)
        check_group_action(n, m, maximal, G)
        check_gamma_is_a_group(n, m, G)
        orb = compute_orbits(maximal, G)
        check_orbits_are_orbits(n, maximal, G, orb)
        check_table2(n, maximal, G, orb)
        if n in satbf:
            check_saturation_closed_under_gamma(n, satbf[n])

    info("")
    info("SCOPE: the census of Table 1 is reproduced in full for n=3,4,5 --"
         " every one of the 57 / 6,476 / 8,656,937 candidates was tested,"
         " nothing was sampled or truncated.")
    info("SCOPE: for n=1,2,3,4 sat*(n,N) is also obtained with no lemma at"
         " all, by scanning every subfamily of B_n. For n=5 that scan"
         " (2^32 subfamilies) is out of reach, so the n=5 minimum rests on"
         " the reduction of Lemma 3, whose computational content (no copy of"
         " N meets {} or [n]) is checked at n=5 and whose full set equality"
         " is checked at n=3,4.")
    info("SCOPE: at n=6 only the 4-subset facts are checked -- the Lemma 4"
         " edge count, the Lemma 2 P4 criterion, role uniqueness and the"
         " no-extremes fact. No census, no minimality and no orbit claim is"
         " verified at n=6, and nothing whatever for n >= 7. The paper claims"
         " nothing for n >= 6 either.")
    info("NOT RE-RUN: the paper's own first implementation, which this bundle"
         " does not ship; this file is an independent reimplementation written"
         " from the paper's specification alone. Also not re-run: the n=5"
         " brute-force scan of all 2^32 subfamilies of B_5, so sat*(5,N) here"
         " rests on the reduction of Lemma 3 as set out in the SCOPE line"
         " above; and no census, minimality or orbit computation at n >= 6.")
    return 0


def verdict():
    failed = _RESULTS.count(False)
    total = len(_RESULTS)
    if failed or total == 0:
        print("VERDICT: %d OF %d CHECKS FAILED" % (max(failed, 1), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:                       # never exit without a verdict
        import traceback
        traceback.print_exc(file=sys.stderr)
        check("no_internal_error", False,
              "%s: %s" % (type(exc).__name__, exc))
    sys.exit(verdict())

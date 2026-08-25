#!/usr/bin/env python3
"""Verification of nine mutually orthogonal cyclic 9-cycle systems of K_19.

Definitions used, as stated in the paper.  An l-cycle system of K_n is a
partition of E(K_n) into l-cycles.  Two such systems are ORTHOGONAL if every
cycle of one shares at most one edge with every cycle of the other.  A system on
vertex set Z_n is CYCLIC if it is invariant under x -> x+1.  mu'(l,n) is the
largest size of a pairwise orthogonal family of cyclic l-cycle systems of K_n.

TAKEN FROM THE PAPER (inputs, transcribed literally; nothing else is assumed):
  * n = 19, cycle length 9, vertex set Z_19.
  * The single base cycle       C = (0,1,3,18,5,13,4,9,12).
  * The multiplier              4, with the nine systems F_r = O(4^r C), r=0..8.
  * The printed signed difference sequence of C, (1,2,-4,6,8,-9,5,3,7).
  * The printed tail vector     a^(0) = (0,1,9,18,4,18,12,5,4).
  * The printed 4-by-9 table of residues a^(0)_d - a^(r)_d, r = 1..4.
  * The claim mu'(9,19) >= 9, and the quoted literature bound mu'(l,n) <= n-3
    whose value at n = 19 the paper reports as 16.
  * The stated prior record of eight such systems.

DERIVED HERE (computed from the inputs above, never copied from the paper):
  * C is decoded and printed back: 9 distinct vertices, 9 edges, 2-regular;
    its signed differences and its difference-d tails are recomputed.
  * 4 has multiplicative order 9 mod 19, so the nine base cycles 4^r C are
    formed and printed; each is developed into O(4^r C) by translation.
  * HYPOTHESES: each of the nine developments consists of 19 nine-cycles whose
    edges partition the 171 edges of K_19 exactly, and each is invariant under
    x -> x+1; the nine systems are pairwise distinct.
  * CONCLUSION (load-bearing, computed): all 36 unordered pairs are orthogonal,
    tested by intersecting the edge sets of all 19x19 pairs of cycles directly
    from the definition -- no lemma, no algebraic shortcut.  The certified
    family size is then recounted and mu'(9,19) >= 9 is read off it.
  * The paper's reduction of 36 pairs to 4 (multiply by 4^{-i}) is verified as a
    permutation of the nine systems, and its four printed table rows and their
    distinctness are recomputed.
  * STRUCTURAL LEMMA: the statement that makes the census exhaustive -- every
    cyclic 9-cycle system of K_19 is the translation orbit O(C) of one base
    cycle carrying exactly one edge of each difference class 1..9, and
    conversely -- is printed in full, its arithmetic steps are verified
    exhaustively on Z_19, and its conclusion is cross-checked end to end
    against a brute-force enumeration of ALL cycle systems in the smallest
    analogous case (n, l) = (7, 3), where no structural assumption is needed.
  * CENSUS: every cyclic 9-cycle system of K_19 is enumerated exhaustively
    (twice, under two different anchorings, and the results compared), the nine
    exhibited systems are located inside it, and the census is scanned for a
    tenth system orthogonal to all nine -- each rejection certified by
    exhibiting two cycles that share two edges.

Standard library only, exact integer arithmetic on Z_19, no floats, no data
files, no network.  Runs in a few seconds.
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


def finish():
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


# ----------------------------------------------------------------------
# Values TAKEN FROM THE PAPER (inputs).
# ----------------------------------------------------------------------
N = 19                     # order of the complete graph, vertex set Z_19
M = 9                      # cycle length
BASE = (0, 1, 3, 18, 5, 13, 4, 9, 12)          # the base cycle C
MULTIPLIER = 4                                  # F_r = 4^r F_0, r = 0..8
FAMILY_SIZE = 9                                 # claimed number of systems
PAPER_SIGNED_DIFFS = (1, 2, -4, 6, 8, -9, 5, 3, 7)
PAPER_A0 = (0, 1, 9, 18, 4, 18, 12, 5, 4)       # tails a^{(0)}_d, d = 1..9
PAPER_TABLE = {                                 # rows (a^{(0)}_d - a^{(r)}_d)
    1: (3, 6, 16, 18, 8, 4, 2, 1, 13),
    2: (16, 1, 12, 11, 5, 8, 9, 6, 2),
    3: (4, 14, 3, 6, 2, 9, 12, 5, 11),
    4: (11, 12, 2, 15, 6, 5, 14, 0, 4),
}
PAPER_LOWER = 9            # claimed  mu'(9,19) >= 9
PAPER_UPPER = 16           # cited    mu'(l,n) <= n-3, evaluated at n = 19


# ----------------------------------------------------------------------
# Basic arithmetic on Z_N.  Pure integers throughout; no floats anywhere.
# ----------------------------------------------------------------------
def edges_of(cyc):
    """Edge list of a closed vertex sequence, as sorted 2-tuples."""
    L = len(cyc)
    return [tuple(sorted((cyc[i], cyc[(i + 1) % L]))) for i in range(L)]


def diff_class(x, y):
    """The difference d in {1..N//2} of the edge {x,y}."""
    d = (y - x) % N
    return d if d <= N // 2 else N - d


def tails_of(cyc):
    """a_d for d = 1..9: tail of the difference-d edge oriented positively."""
    t = {}
    L = len(cyc)
    for i in range(L):
        x, y = cyc[i], cyc[(i + 1) % L]
        d = diff_class(x, y)
        tail = x if (y - x) % N == d else y
        if d in t:
            return None            # repeated difference: no tail vector
        t[d] = tail
    if sorted(t) != list(range(1, N // 2 + 1)):
        return None
    return tuple(t[d] for d in range(1, N // 2 + 1))


def scale(u, cyc):
    return tuple((u * v) % N for v in cyc)


def shift(s, cyc):
    return tuple((v + s) % N for v in cyc)


def eset(cyc):
    return frozenset(edges_of(cyc))


def develop(cyc):
    """The translation orbit O(C) as a list of edge sets."""
    return [eset(shift(t, cyc)) for t in range(N)]


def check_base_cycle_wellformed():
    """(1) Decode the exhibited base cycle, count it, print it back."""
    ok = (len(BASE) == M
          and len(set(BASE)) == M
          and all(isinstance(v, int) and 0 <= v < N for v in BASE))
    ed = edges_of(BASE)
    ok = ok and len(set(ed)) == M and all(a != b for a, b in ed)
    deg = {}
    for a, b in ed:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    ok = ok and set(deg.values()) == {2} and len(deg) == M
    print("  base cycle C = (" + ",".join(str(v) for v in BASE) + ")")
    print("  edges of C   = " + " ".join("%d-%d" % e for e in ed))
    return ck("base_cycle_is_a_9_cycle_on_Z19", ok,
              "%d distinct vertices, %d edges, all degrees 2" % (len(set(BASE)), len(ed)))


def check_signed_differences():
    """(2) Recompute the signed consecutive differences and compare to paper."""
    sd = []
    for i in range(M):
        r = (BASE[(i + 1) % M] - BASE[i]) % N
        sd.append(r if r <= N // 2 else r - N)
    ok = tuple(sd) == PAPER_SIGNED_DIFFS
    absd = sorted(abs(x) for x in sd)
    ok2 = absd == list(range(1, N // 2 + 1))
    return ck("signed_differences_match_and_hit_every_class", ok and ok2,
              "computed %s; |d| multiset %s" % (tuple(sd), tuple(absd)))


def check_tail_vector():
    """(3) Recompute the tail vector a^{(0)} and compare to the paper."""
    a0 = tails_of(BASE)
    ok = a0 is not None and a0 == PAPER_A0
    return ck("tail_vector_a0_matches_paper", ok, "computed %s" % (a0,))


def check_multiplier_order():
    """(4) The paper's claim that 4 has multiplicative order 9 modulo 19."""
    pows = [pow(MULTIPLIER, k, N) for k in range(1, 10)]
    order = None
    for k in range(1, N):
        if pow(MULTIPLIER, k, N) == 1:
            order = k
            break
    ok = (order == M and pows[8] == 1 and pows[0] != 1 and pows[2] == 7)
    return ck("multiplier_4_has_order_9_mod_19", ok,
              "order=%s, 4^3=%d, powers=%s" % (order, pows[2], tuple(pows)))


def family_bases():
    """The nine base cycles 4^r C, r = 0..8, derived from the paper's inputs."""
    return [scale(pow(MULTIPLIER, r, N), BASE) for r in range(FAMILY_SIZE)]


def check_family_is_cycle_systems(bases, systems):
    """(2) Every hypothesis: each F_r is a 9-cycle system of K_19."""
    ok = True
    detail = []
    allpairs = set()
    for x in range(N):
        for y in range(x + 1, N):
            allpairs.add((x, y))
    why = []
    for r, (b, S) in enumerate(zip(bases, systems)):
        if len(b) != M or len(set(b)) != M:
            ok = False
            why.append("r=%d base has %d distinct of %d vertices"
                       % (r, len(set(b)), len(b)))
        if len(S) != N:
            ok = False
            why.append("r=%d has %d parts, not %d" % (r, len(S), N))
        multi = []
        for c in S:
            if len(c) != M:
                ok = False
                why.append("r=%d has a part with %d edges, not %d"
                           % (r, len(c), M))
            multi.extend(c)
        if len(multi) != len(allpairs) or set(multi) != allpairs:
            ok = False
            why.append("r=%d covers %d edges (%d distinct) of %d"
                       % (r, len(multi), len(set(multi)), len(allpairs)))
        if len(set(multi)) != len(multi):
            ok = False           # some edge covered twice
            why.append("r=%d covers some edge twice" % r)
        detail.append(len(S))
    return ck("each_of_9_families_partitions_E_K19", ok,
              ("%d cycles of length %d covering all %d edges, per system; "
               "cycle counts %s" % (N, M, len(allpairs), tuple(detail))) if ok
              else "; ".join(sorted(set(why))[:4]))


def is_single_9_cycle(E):
    """True iff the edge set E is one connected 9-cycle: 9 loopless edges on 9
    vertices, every degree 2, and connected.  Read off the edges alone, so it
    does not trust the vertex sequence the edges were generated from."""
    adj = {}
    for a, b in E:
        if a == b:
            return False
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    if len(E) != M or len(adj) != M:
        return False
    if any(len(nb) != 2 for nb in adj.values()):
        return False
    seen = {next(iter(adj))}
    stack = list(seen)
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == M


def check_cyclicity(systems):
    """(2) Each system is a set of 19 DISTINCT cycles, each of them verified
    from its edge set alone to be one connected 9-cycle, and the set is
    invariant under the translation x -> x+1."""
    ok = True
    orbits = []
    why = []
    for r, S in enumerate(systems):
        SS = set(S)
        orbits.append(len(SS))
        if len(SS) != N:
            ok = False                      # orbit collapsed: fewer than 19
            why.append("r=%d has %d distinct cycles, not %d" % (r, len(SS), N))
        nbad = sum(1 for c in SS if not is_single_9_cycle(c))
        if nbad:
            ok = False                      # 9 edges, but not one 9-cycle
            why.append("r=%d: %d parts are not one connected 9-cycle"
                       % (r, nbad))
        shifted = set()
        for c in SS:
            shifted.add(frozenset(tuple(sorted(((a + 1) % N, (b + 1) % N)))
                                  for (a, b) in c))
        if shifted != SS:
            ok = False
            why.append("r=%d not fixed by x->x+1" % r)
    return ck("every_system_is_cyclic_under_translation", ok,
              ("all %d systems: %s distinct connected %d-cycles each, fixed by "
               "x->x+1" % (len(systems), tuple(orbits), M)) if ok
              else "; ".join(why[:4]))


def check_systems_distinct(systems):
    """(1) The nine systems are pairwise distinct as sets of cycles."""
    keys = [frozenset(S) for S in systems]
    ok = len(set(keys)) == FAMILY_SIZE
    return ck("nine_systems_pairwise_distinct", ok,
              "%d distinct among %d" % (len(set(keys)), len(keys)))


def orthogonal_bruteforce(S1, S2):
    """True iff no cycle of S1 shares two or more edges with a cycle of S2.

    This is the definition applied directly to the 19 x 19 pairs of cycles;
    it uses no lemma and no algebraic shortcut.
    """
    worst = 0
    for A in S1:
        for B in S2:
            k = len(A & B)
            if k > worst:
                worst = k
    return worst <= 1, worst


def check_pairwise_orthogonality(systems):
    """(3) LOAD-BEARING: all 36 unordered pairs are orthogonal, by brute force."""
    pairs = list(combinations(range(FAMILY_SIZE), 2))
    bad = []
    worst_overall = 0
    for i, j in pairs:
        good, worst = orthogonal_bruteforce(systems[i], systems[j])
        worst_overall = max(worst_overall, worst)
        if not good:
            bad.append((i, j, worst))
    ok = (len(pairs) == 36 and not bad)
    return ck("all_36_pairs_orthogonal_by_bruteforce", ok,
              "%d pairs tested over %d cycle-pairs each; max shared edges = %d%s"
              % (len(pairs), N * N, worst_overall,
                 "" if not bad else "; offenders " + str(bad[:4])))


def check_self_intersection_is_detectable(systems):
    """A control: a system is NOT orthogonal to itself, so the brute-force
    orthogonality test above is capable of returning False on this data."""
    good, worst = orthogonal_bruteforce(systems[0], systems[0])
    ok = (not good) and worst == M
    return ck("orthogonality_test_rejects_a_system_against_itself", ok,
              "max shared edges with itself = %d (expected %d)" % (worst, M))


def check_conclusion_lower_bound(bases, systems):
    """(3) The settlement: recount, from the verified data only, a pairwise
    orthogonal family of distinct cyclic 9-cycle systems, and read off the
    bound mu'(9,19) >= (its size)."""
    good_systems = []
    allpairs = set((x, y) for x in range(N) for y in range(x + 1, N))
    for S in systems:
        multi = []
        for c in S:
            multi.extend(c)
        if len(S) == N and len(multi) == len(allpairs) and set(multi) == allpairs \
                and all(len(c) == M for c in S):
            good_systems.append(S)
    bykey = {}
    for S in good_systems:
        bykey[frozenset(S)] = S
    reps = list(bykey.values())
    # a family certifies the bound only if it is a clique in the
    # orthogonality graph; otherwise its certified size is 0.
    clique = all(orthogonal_bruteforce(reps[i], reps[j])[0]
                 for i, j in combinations(range(len(reps)), 2))
    size = len(reps) if clique else 0
    ok = size >= PAPER_LOWER
    return ck("computed_family_size_certifies_mu_ge_9", ok,
              "certified mutually orthogonal cyclic 9-cycle systems: %d, so "
              "mu'(9,19) >= %d" % (size, size))


def check_paper_table(bases):
    """(4) Recompute the four printed rows (a^{(0)}_d - a^{(r)}_d) and the
    paper's assertion that each row holds nine distinct residues."""
    a0 = tails_of(bases[0])
    if a0 is None:
        return ck("printed_table_rows_reproduce_and_are_distinct", False,
                  "base cycle has no tail vector: it misses a difference class")
    ok = True
    detail = []
    rows_ok = []                    # rows that reproduced AND were distinct
    for r in sorted(PAPER_TABLE):
        ar = tails_of(bases[r])
        if ar is None:
            ok = False
            detail.append("r=%d base cycle has no tail vector" % r)
            continue
        row = tuple((a0[d] - ar[d]) % N for d in range(M))
        good = True
        if row != PAPER_TABLE[r]:
            ok = False
            good = False
            detail.append("r=%d computed %s != printed %s" % (r, row, PAPER_TABLE[r]))
        if len(set(row)) != M:
            ok = False
            good = False
            detail.append("r=%d row not distinct" % r)
        if good:
            rows_ok.append(r)
    # the count of verified rows is derived from rows_ok, not asserted
    ok = ok and len(PAPER_TABLE) > 0 and len(rows_ok) == len(PAPER_TABLE)
    return ck("printed_table_rows_reproduce_and_are_distinct", ok,
              ("rows r=%s reproduced -- %d of the %d printed rows -- each with "
               "%d distinct residues"
               % (",".join(str(r) for r in rows_ok), len(rows_ok),
                  len(PAPER_TABLE), M))
              if not detail else "; ".join(detail))


def orth_criterion(a, b):
    """The paper's Lemma: O(C) and O(D) are orthogonal iff the nine residues
    a_d - b_d are pairwise distinct.  Translation-invariant in both arguments."""
    return len({(a[d] - b[d]) % N for d in range(M)}) == M


def check_criterion_agrees_with_definition(census):
    """(4/5) The Lemma is the paper's engine for the census below, so test it
    against the raw definition on a deterministic sample containing BOTH
    orthogonal and non-orthogonal pairs."""
    step = max(1, len(census) // 12)
    samp = [census[i] for i in range(0, len(census), step)][:12]
    devs = [develop(c) for c in samp]
    tls = [tails_of(c) for c in samp]
    agree = 0
    disagree = []
    n_true = n_false = 0
    for i, j in combinations(range(len(samp)), 2):
        by_lemma = orth_criterion(tls[i], tls[j])
        by_def = orthogonal_bruteforce(devs[i], devs[j])[0]
        if by_lemma == by_def:
            agree += 1
        else:
            disagree.append((i, j, by_lemma, by_def))
        if by_def:
            n_true += 1
        else:
            n_false += 1
    ok = (not disagree) and n_true > 0 and n_false > 0
    return ck("lemma_criterion_matches_raw_definition_on_sample", ok,
              "%d/%d sample pairs agree (%d orthogonal, %d not)%s"
              % (agree, agree + len(disagree), n_true, n_false,
                 "" if not disagree else "; disagreements " + str(disagree[:3])))


# ----------------------------------------------------------------------
# THE STRUCTURAL LEMMA that makes the census exhaustive.
#
# LEMMA.  Let p be an odd prime and l = (p-1)/2.  Then every cyclic l-cycle
# system of K_p is the translation orbit O(C) of a single base cycle C that
# carries exactly one edge of each difference class 1..l; and conversely, the
# orbit of any such C is a cyclic l-cycle system.
#
#   (L1) |E(K_p)| = p(p-1)/2 = l*p, so an l-cycle system has exactly p cycles,
#        and a cyclic system is a set of p cycles permuted by Z_p.
#   (L2) p is prime, so a nonzero translation generates all of Z_p; every edge
#        of K_p has trivial stabiliser, so the edges split into l translation
#        orbits of size p -- exactly the difference classes.  A nonempty
#        translation-invariant set of edges therefore has at least p > l edges,
#        so no l-cycle can be fixed by a nonzero translation.
#   (L3) hence every orbit of cycles has size p (its size divides the prime p
#        and is not 1), and a system holding exactly p cycles is one orbit.
#   (L4) O(C) covers difference class d exactly p*m_d times, where m_d is the
#        number of difference-d edges of C; class d holds p edges, so an exact
#        partition forces m_d = 1 for every d -- and conversely, m_d = 1 for
#        every d makes O(C) an exact partition, by the orbit fact in (L2).
#
# (L1), (L2) and (L4) are verified exhaustively on Z_19 by the two checks below;
# (L3) is the group-theoretic step they feed.  The whole conclusion is then
# cross-checked against brute force in the smallest analogous case, p = 7,
# l = 3, where every partition of E(K_7) into 3-cycles can be enumerated with no
# structural assumption whatsoever.
# ----------------------------------------------------------------------
def binom(a, b):
    """Exact binomial coefficient, integer arithmetic only."""
    if b < 0 or b > a:
        return 0
    num = den = 1
    for k in range(b):
        num *= a - k
        den *= k + 1
    return num // den


def check_edge_orbits_under_translation():
    """LEMMA STEPS L1 and L2, verified exhaustively on Z_19: 19 is prime, every
    one of the 171 edges has trivial translation stabiliser, and the edge orbits
    are precisely the nine difference classes, each of size 19."""
    divisors = [k for k in range(1, N + 1) if N % k == 0]
    prime = (divisors == [1, N])
    alledges = [tuple(sorted((x, y)))
                for x in range(N) for y in range(N) if x < y]
    nedges = len(alledges)
    nstab = 0                      # count the tests; do not assert their number
    trivial_stab = True
    for (a, b) in alledges:
        for t in range(1, N):
            nstab += 1
            if tuple(sorted(((a + t) % N, (b + t) % N))) == (a, b):
                trivial_stab = False
    orbits = set()
    for (a, b) in alledges:
        orbits.add(frozenset(tuple(sorted(((a + t) % N, (b + t) % N)))
                             for t in range(N)))
    sizes = sorted(set(len(o) for o in orbits))
    classes = {}
    for e in alledges:
        classes.setdefault(diff_class(*e), set()).add(e)
    orbits_are_classes = (orbits ==
                          set(frozenset(v) for v in classes.values()))
    ok = (prime and nedges == N * (N - 1) // 2 and nedges == M * N
          and nstab == nedges * (N - 1)
          and trivial_stab and len(orbits) == M and sizes == [N]
          and orbits_are_classes
          and sorted(classes) == list(range(1, N // 2 + 1))
          and all(len(v) == N for v in classes.values()))
    return ck("edge_orbits_are_the_nine_difference_classes_of_size_19", ok,
              "%d prime (divisors %s); %d = %d*%d edges, so a %d-cycle system has "
              "exactly %d cycles; all %d = %dx%d nonzero-translation stabiliser "
              "tests trivial; %d edge orbits of sizes %s, equal to the difference "
              "classes -- so an invariant edge set has >= %d edges and no %d-cycle "
              "is translation-fixed"
              % (N, tuple(divisors), nedges, M, N, M, N, nstab, nedges, N - 1,
                 len(orbits), tuple(sizes), N, M))


def check_one_edge_per_difference_class():
    """LEMMA STEP L4, both directions, verified exhaustively.

    FORCED: enumerate every multiplicity vector (m_1..m_9) of nonnegative
    integers summing to 9 -- the possible difference multisets of a 9-edge base
    cycle -- and confirm that exactly one of them meets every class, the
    all-ones vector.  SUFFICIENT: for every edge e and its class d, the 19
    translates of e are the 19 edges of class d, each exactly once, so a base
    cycle with one edge per class develops to an exact partition of E(K_19).
    """
    half = N // 2
    counts = [0, 0, 0]          # total vectors, vectors meeting every class,
    #                             all-ones vectors

    def gen(i, left, vec):
        if i == half:
            if left == 0:
                counts[0] += 1
                if all(x >= 1 for x in vec):
                    counts[1] += 1
                    if all(x == 1 for x in vec):
                        counts[2] += 1
            return
        for k in range(left + 1):
            gen(i + 1, left - k, vec + (k,))

    gen(0, M, ())
    total, surj, ones = counts
    expected_total = binom(M + half - 1, half - 1)
    classes = {}
    for x in range(N):
        for y in range(x + 1, N):
            classes.setdefault(diff_class(x, y), set()).add((x, y))
    suff = True
    for d in range(1, half + 1):
        cls = classes[d]
        for (a, b) in sorted(cls):
            tr = [tuple(sorted(((a + t) % N, (b + t) % N))) for t in range(N)]
            if len(set(tr)) != N or set(tr) != cls:
                suff = False
    ok = (total == expected_total and surj == 1 and ones == 1 and suff)
    return ck("one_edge_per_difference_class_forced_and_sufficient", ok,
              "%d multiplicity vectors of %d edges over %d classes (= C(%d,%d) "
              "= %d); exactly %d meets every class and it is the all-ones vector "
              "(%d); and the %d translates of every edge exhaust its class once "
              "each = %s"
              % (total, M, half, M + half - 1, half - 1, expected_total,
                 surj, ones, N, suff))


def edges_of_gen(n, cyc):
    """Edge set of a closed vertex sequence in K_n, as a frozenset."""
    L = len(cyc)
    return frozenset(tuple(sorted((cyc[i], cyc[(i + 1) % L])))
                     for i in range(L))


def develop_gen(n, cyc):
    """The translation orbit of `cyc` in K_n, as a set of edge sets."""
    L = len(cyc)
    return frozenset(
        frozenset(tuple(sorted((((cyc[i] + t) % n), ((cyc[(i + 1) % L] + t) % n))))
                  for i in range(L))
        for t in range(n))


def translate_system(n, S):
    """Image of a set of cycles under x -> x+1 in K_n."""
    return frozenset(
        frozenset(tuple(sorted(((a + 1) % n, (b + 1) % n))) for (a, b) in c)
        for c in S)


def all_cycle_systems_bruteforce(n, m):
    """EVERY partition of E(K_n) into m-cycles, assuming nothing.

    A plain exact cover over the complete list of m-cycles of K_n: no
    cyclicity, no difference classes, no lemma.  Only tractable for small n,
    which is precisely why it can test the structural lemma itself rather than
    the lemma's own consequences.
    """
    alledges = frozenset(tuple(sorted((x, y)))
                         for x in range(n) for y in range(n) if x < y)
    cycles = set()
    for combo in combinations(range(n), m):
        for perm in permutations(combo[1:]):
            E = edges_of_gen(n, (combo[0],) + perm)
            if len(E) == m:
                cycles.add(E)
    byedge = {}
    for E in cycles:
        for e in E:
            byedge.setdefault(e, []).append(E)
    systems = []

    def rec(uncovered, chosen):
        if not uncovered:
            systems.append(frozenset(chosen))
            return
        e = min(uncovered)                     # each partition reached once
        for E in byedge.get(e, ()):
            if E <= uncovered:
                rec(uncovered - E, chosen + (E,))

    rec(alledges, ())
    return systems, len(cycles)


def check_structural_lemma_against_bruteforce():
    """The structural lemma has the same statement for every odd prime p with
    cycle length (p-1)/2.  At (p, l) = (7, 3) -- the smallest such case -- ALL
    partitions of E(K_7) into 3-cycles are enumerated with no structural
    assumption, the translation-invariant ones are selected, and the result is
    compared with the very routine used at n = 19.  This is the one test that
    can catch an error in the lemma itself rather than in its application."""
    n = 7
    m = (n - 1) // 2                             # the same shape as (19, 9)
    anchor_hi, anchor_lo = m, 1
    allsys, ncyc = all_cycle_systems_bruteforce(n, m)
    distinct = set(allsys)
    cyclic = set(S for S in distinct if translate_system(n, S) == S)
    lem = enumerate_census_gen(n, m, anchor_hi)
    lemsys = set(develop_gen(n, c) for c in lem)
    lem1 = enumerate_census_gen(n, m, anchor_lo)
    lemsys1 = set(develop_gen(n, c) for c in lem1)
    ok = (2 * m + 1 == n and 2 * M + 1 == N      # both are the (p,(p-1)/2) shape
          and len(distinct) == len(allsys)
          and len(cyclic) > 0                    # something to compare against
          and len(cyclic) < len(distinct)        # the invariance filter bites
          and len(lemsys) == len(lem)            # one output per system
          and len(lemsys1) == len(lem1)
          and lemsys == cyclic and lemsys1 == cyclic)
    return ck("structural_lemma_matches_bruteforce_in_the_K7_analogue", ok,
              "brute force over all %d %d-cycles of K_%d finds %d partitions of "
              "E(K_%d) into %d-cycles, %d of them invariant under x->x+1; the "
              "lemma-based enumerator used at n=%d returns %d base cycles "
              "(anchor d=%d) and %d (anchor d=%d); developments equal the "
              "brute-force cyclic systems = %s"
              % (ncyc, m, n, len(distinct), n, m, len(cyclic), N, len(lem),
                 anchor_hi, len(lem1), anchor_lo,
                 lemsys == cyclic and lemsys1 == cyclic))


def enumerate_census_gen(n, m, anchor):
    """Every cyclic m-cycle system of K_n, once each, for m = (n-1)/2.

    By the structural lemma above (verified, not assumed: see the two checks
    on Z_19 and the K_7 brute-force cross-check), such a system is a single
    translation orbit O(C) whose base cycle carries one edge of each difference
    1..m.  Anchoring the walk at the tail of the difference-`anchor` edge,
    traversed positively, from vertex 0, gives exactly one vertex sequence per
    system.
    """
    out = []
    rest = set(d for d in range(1, n // 2 + 1) if d != anchor)

    def rec(v, used, path, rem):
        if len(rem) == 1:
            d = next(iter(rem))
            if v == d or v == (n - d) % n:
                out.append(tuple(path))
            return
        for d in tuple(rem):
            rem.discard(d)
            for s in (d, n - d):
                w = (v + s) % n
                if w in used:
                    continue
                used.add(w)
                path.append(w)
                rec(w, used, path, rem)
                path.pop()
                used.discard(w)
            rem.add(d)

    rec(anchor, {0, anchor}, [0, anchor], rest)
    return out


def enumerate_census(anchor):
    """Every cyclic 9-cycle system of K_19, once each (n = 19, m = 9)."""
    return enumerate_census_gen(N, M, anchor)


def canonical(cyc, anchor=M):
    """Canonical representative of O(cyc): rotate/reflect/translate so the
    difference-`anchor` edge runs from 0 to `anchor`."""
    L = len(cyc)
    for i in range(L):
        x, y = cyc[i], cyc[(i + 1) % L]
        if (y - x) % N == anchor:
            return tuple((cyc[(i + k) % L] - x) % N for k in range(L))
        if (x - y) % N == anchor:
            return tuple((cyc[(i + 1 - k) % L] - y) % N for k in range(L))
    return None


def check_census_consistent(census9, census1, bases):
    """(5) Two independently anchored exhaustive enumerations of ALL cyclic
    9-cycle systems of K_19 must yield the same set of systems, and the nine
    exhibited systems must be members of it."""
    set9 = set(canonical(c) for c in census9)
    set1 = set(canonical(c) for c in census1)
    fam = set(canonical(b) for b in bases)
    ok = (len(set9) == len(census9) and len(set1) == len(census1)
          and set9 == set1 and len(fam) == FAMILY_SIZE and fam <= set9)
    every = all(tails_of(c) is not None for c in census9)
    ok = ok and every
    return ck("exhaustive_census_agrees_and_contains_the_nine", ok,
              "anchor d=9 gives %d systems, anchor d=1 gives %d, sets equal=%s, "
              "nine exhibited systems present=%s"
              % (len(set9), len(set1), set9 == set1, fam <= set9))


def check_family_is_maximal(census, bases):
    """(5) Census-wide: is there a tenth cyclic 9-cycle system of K_19
    orthogonal to all nine?  Each rejection is certified constructively by
    exhibiting two cycles that share two edges, so the answer does not rest on
    the Lemma."""
    nfam = len(bases)
    famT = [tails_of(b) for b in bases]
    famE = [eset(b) for b in bases]
    if any(t is None for t in famT):
        return ck("no_tenth_cyclic_system_extends_the_family", False,
                  "a family member is not a perfect difference cycle")
    famcanon = set(canonical(b) for b in bases)
    extenders = []
    unwitnessed = 0
    best_member = 0        # over the nine family members themselves
    best_outside = 0       # over every OTHER cyclic system: the real ceiling
    for c in census:
        a = tails_of(c)
        good = 0
        firstwit = None
        for r in range(nfam):
            b = famT[r]
            seen = {}
            wit = None
            for d in range(M):
                s = (a[d] - b[d]) % N
                if s in seen:
                    wit = s
                    break
                seen[s] = d
            if wit is None:
                good += 1                  # no early exit: count ALL nine
            elif firstwit is None:
                firstwit = (r, wit)
        if firstwit is None:
            extenders.append(c)
        else:
            r, wit = firstwit
            # c - wit shares its difference-d and difference-d' edges with 4^r C
            if len(eset(shift(-wit, c)) & famE[r]) < 2:
                unwitnessed += 1
        if canonical(c) in famcanon:
            best_member = max(best_member, good)
        else:
            best_outside = max(best_outside, good)
    ok = (unwitnessed == 0 and not extenders
          and best_member == nfam - 1 and best_outside < nfam)
    return ck("no_tenth_cyclic_system_extends_the_family", ok,
              "scanned all %d cyclic 9-cycle systems: %d extend the family; "
              "every rejection carries a 2-shared-edge witness (%d without); "
              "each family member is orthogonal to %d of the other %d; the best "
              "any system outside the family achieves is %d of %d"
              % (len(census), len(extenders), unwitnessed, best_member,
                 nfam - 1, best_outside, nfam))


def check_multiplier_conjugation(systems):
    """(4) The paper reduces the 36 pairs to four by multiplying by 4^{-i}.
    Verify the reduction as stated: 4^{-i} carries F_i to F_0 and F_j to
    F_{(j-i) mod 9}, as sets of cycles."""
    keys = [frozenset(S) for S in systems]
    inv = pow(MULTIPLIER, N - 2, N)        # N is prime, so this is 4^{-1}
    ok = True
    detail = []
    nhit = 0                       # images that landed correctly, counted
    for i in range(FAMILY_SIZE):
        u = pow(inv, i, N)
        for j in range(FAMILY_SIZE):
            mapped = frozenset(
                frozenset(tuple(sorted(((u * a) % N, (u * b) % N)))
                          for (a, b) in c) for c in systems[j])
            if mapped != keys[(j - i) % FAMILY_SIZE]:
                ok = False
                detail.append((i, j))
            else:
                nhit += 1
    ok = ok and nhit == FAMILY_SIZE * FAMILY_SIZE
    return ck("multiplication_by_4_inverse_i_shifts_the_index", ok,
              "all %d = %dx%d (i,j) images land on F_{(j-i) mod %d}"
              % (nhit, FAMILY_SIZE, FAMILY_SIZE, FAMILY_SIZE)
              if ok else "misses at " + str(detail[:4]))


def check_numeric_claims(census, systems):
    """(4) Every count a reader would rely on, recomputed.  The family size is
    taken from the derived data (distinct verified systems), not asserted."""
    derived_size = len({frozenset(S) for S in systems})
    facts = []
    facts.append(("n choose 2 = 171", (N * (N - 1)) // 2 == 171))
    facts.append(("9 * 19 = C(19,2)", M * N == (N * (N - 1)) // 2))
    facts.append(("C(9,2) = 36 pairs", len(list(combinations(range(derived_size), 2))) == 36))
    facts.append(("cited bound n-3 = 16", N - 3 == PAPER_UPPER))
    facts.append(("lower <= upper", PAPER_LOWER <= PAPER_UPPER))
    facts.append(("derived family size equals claimed 9", derived_size == PAPER_LOWER))
    facts.append(("derived family size exceeds previous record 8", derived_size > 8))
    facts.append(("census is nonempty and finite", 0 < len(census) < 10 ** 7))
    bad = [name for name, good in facts if not good]
    return ck("numeric_claims_recomputed", not bad,
              "%d/%d facts hold%s" % (len(facts) - len(bad), len(facts),
                                      "" if not bad else "; failed " + str(bad)))


def run(fn, *args):
    """Run one check; an exception counts as a failure, never as a pass."""
    try:
        return fn(*args)
    except Exception as exc:                                  # noqa: BLE001
        return ck(fn.__name__, False, "raised %s" % type(exc).__name__)


def main():
    print("Nine mutually orthogonal cyclic 9-cycle systems of K_19: verification")
    print("n = %d, cycle length = %d, multiplier = %d" % (N, M, MULTIPLIER))
    run(check_base_cycle_wellformed)
    run(check_signed_differences)
    run(check_tail_vector)
    run(check_multiplier_order)

    bases = family_bases()
    print("  the nine base cycles 4^r C, r = 0..8:")
    for r, b in enumerate(bases):
        print("    r=%d: (%s)" % (r, ",".join(str(v) for v in b)))
    systems = [develop(b) for b in bases]

    run(check_family_is_cycle_systems, bases, systems)
    run(check_cyclicity, systems)
    run(check_systems_distinct, systems)
    run(check_pairwise_orthogonality, systems)
    run(check_self_intersection_is_detectable, systems)
    run(check_conclusion_lower_bound, bases, systems)
    run(check_paper_table, bases)
    run(check_multiplier_conjugation, systems)

    print("  STRUCTURAL LEMMA, the statement that makes the census below")
    print("  exhaustive (verified here, not assumed).  For an odd prime p and")
    print("  cycle length l = (p-1)/2 -- so (p,l) = (19,9) -- every cyclic")
    print("  l-cycle system of K_p is the translation orbit O(C) of a single base")
    print("  cycle C carrying exactly one edge of each difference class 1..l, and")
    print("  conversely the orbit of any such C is a cyclic l-cycle system.")
    print("    (L1) |E(K_p)| = l*p, so an l-cycle system has exactly p cycles and")
    print("         Z_p permutes them.")
    print("    (L2) p prime, so a nonzero translation generates Z_p; every edge")
    print("         orbit has size p (the l difference classes), hence a nonempty")
    print("         invariant edge set has >= p > l edges and no l-cycle is fixed")
    print("         by a nonzero translation.")
    print("    (L3) so every orbit of cycles has size p (it divides the prime p")
    print("         and is not 1), and a system of exactly p cycles is one orbit.")
    print("    (L4) O(C) covers class d exactly p*m_d times with m_d the number of")
    print("         difference-d edges of C, and class d holds p edges; an exact")
    print("         partition forces m_d = 1 for all d, and conversely.")
    print("  (L1),(L2),(L4) are verified exhaustively on Z_19 by the next two")
    print("  checks; (L3) is the group step they feed; the whole conclusion is")
    print("  then cross-checked against brute force at (p,l) = (7,3).")
    run(check_edge_orbits_under_translation)
    run(check_one_edge_per_difference_class)
    run(check_structural_lemma_against_bruteforce)

    census9 = enumerate_census(9)
    census1 = enumerate_census(1)
    run(check_census_consistent, census9, census1, bases)
    run(check_criterion_agrees_with_definition, census9)
    run(check_family_is_maximal, census9, bases)
    run(check_numeric_claims, census9, systems)

    # Every number quoted in the scope paragraph is derived from the data just
    # computed, never a literal: the census size and the pair count especially.
    ncensus = len(set(canonical(c) for c in census9))
    npairs = len(list(combinations(range(FAMILY_SIZE), 2)))
    print("NOT RE-RUN HERE: (i) the report of eight mutually orthogonal cyclic "
          "%d-cycle systems attributed to earlier work, which lives in an "
          "external table; (ii) the upper bound mu'(l,n) <= n-3, quoted from "
          "the literature -- only its evaluation %d-3 = %d is checked; "
          "(iii) whether some entirely different family of ten mutually "
          "orthogonal cyclic %d-cycle systems of K_%d exists -- that is a "
          "maximum-clique question on the %d-vertex orthogonality graph "
          "and is out of budget; what is settled here is that THIS family of "
          "nine admits no tenth member; (iv) the census is exhaustive only "
          "modulo the structural lemma printed above: its arithmetic steps are "
          "verified exhaustively on Z_%d and its conclusion is confirmed "
          "against brute force only in the small analogue (p,l) = (7,3). No "
          "assumption-free enumeration of all %d-cycle systems of K_%d is "
          "attempted, and none of (i)-(iv) touches Theorem 1: the family of "
          "nine and its %d orthogonal pairs are established directly from the "
          "definition, with no lemma and no census."
          % (M, N, N - 3, M, N, ncensus, N, M, N, npairs))
    finish()


if __name__ == "__main__":
    main()

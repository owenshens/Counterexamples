#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- referee's verification program for

    "The p=3 Case of a Conjecture of Chauhan, Shukla, and Vinayak", proving Conjecture 4.2(ii) of Chauhan-Shukla-Vinayak at p=3.

Standard library only.  Exact integer / set / bitmask arithmetic throughout;
no floating point is used in any decision.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs, never re-derived -- they are the object
under test, and every one of them is a possible point of failure below):
  * the two edited blocks of L_n, before and after the repair, and their
    1-based position windows: n=13, positions 112-117; n=14, positions
    146-155 (equations (3) and (4) of the paper);
  * the claimed spanning sets S_13 (1 omitted triple) and S_14 (8 omitted
    triples) of Proposition 3;
  * the obstruction table of Section 2 (n, j, T_j, earlier witnesses T_i);
  * the audit table numbers: 169 / 238 facets, 14,196 / 28,203 pairs,
    0 exchange failures, 0 purity failures, 1 / 8 spanning facets;
  * the control numbers 16 and 25 at n=15,16 (section "Exact verification");
  * the abstract's homotopy types S^9 and \\/^8 S^10;
  * this program's own published in-block audit -- 15 of 15 and 44 of 45 pairs
    discharged, one residual pair (147,154) with four surviving witnesses --
    which checks D6 and D8 now GATE rather than merely print.

WHAT THIS PROGRAM DERIVES (the checks):
  * the facet set of Delta_3(C_n^3) from the definition, twice, by two
    independent routes (edge count vs. explicit graph connectivity);
  * that every triple the paper prints is a genuine facet, that the repaired
    block is a permutation of the original, and that the windows fit;
  * that the four listed obstructions are structurally real (the vertex 2 has
    no facet exchange candidate at all, in any order) and that the repaired
    block order -- and not the original block order -- discharges each of them
    from strictly inside the block;
  * an INDEPENDENT shelling order of Delta_3(C_n^3), built here from scratch,
    verified against Lemma 2's exchange condition for every one of the
    C(m,2) pairs, verified again by the separate Bjoerner-Wachs purity test,
    and its spanning facets counted -- giving shellability and the counts
    1, 8, 16, 25 without using the paper's order;
  * the full f-vector and reduced Euler characteristic of Delta_3(C_n^3) by
    enumeration over all 2^n vertex subsets, which forces the spanning count;
  * the reduced Betti numbers over GF(2) of Delta_3(C_n^3), which pin the
    rational homology exactly and exclude 2-torsion, independently of
    shellability;
  * the arithmetic identities the reader relies on: 4p+1=13, 6p-3=15,
    (n^2-4np-n+2)/2 = (n^2-13n+2)/2 at p=3, C(n-6,2)-20 = (n^2-13n+2)/2;
  * negative controls proving the verifier is not vacuous.

WHAT IS **NOT** RE-RUN -- see the "NOT REPRODUCED" block printed at the end.
"""

import itertools
import sys
import time

CHECKS = []


def check(name, ok, detail=""):
    """Record and print one check.  Every call must be able to print FAIL."""
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    sys.stdout.flush()
    return bool(ok)


def info(msg):
    print("# " + msg)
    sys.stdout.flush()


#############################################################################
# SECTION 0.  VALUES TAKEN FROM THE PAPER.  Inputs only -- nothing here is
# derived, and every entry is a place where a corruption must be detected.
#############################################################################

P = 3                       # the case of Conjecture 4.2(ii) settled

# (3): n = 13, positions 112--117 of L_13, before -> after
POS_13 = (112, 117)
BLOCK_13_OLD = [(3, 5, 12), (4, 5, 11), (4, 5, 12),
                (5, 9, 12), (5, 10, 12), (5, 11, 12)]
BLOCK_13_NEW = [(4, 5, 11), (5, 11, 12), (5, 10, 12),
                (5, 9, 12), (3, 5, 12), (4, 5, 12)]

# (4): n = 14, positions 146--155 of L_14, before -> after
POS_14 = (146, 155)
BLOCK_14_OLD = [(3, 5, 13), (4, 5, 11), (4, 5, 12), (4, 5, 13), (5, 9, 12),
                (5, 9, 13), (5, 10, 12), (5, 10, 13), (5, 11, 12), (5, 11, 13)]
BLOCK_14_NEW = [(4, 5, 11), (4, 5, 12), (5, 9, 12), (5, 10, 12), (5, 11, 12),
                (5, 11, 13), (5, 10, 13), (5, 9, 13), (3, 5, 13), (4, 5, 13)]

# Proposition 3: the claimed spanning facets, as omitted triples
SPAN_13 = [(4, 10, 12)]
SPAN_14 = [(0, 7, 13), (6, 12, 13), (1, 8, 13), (5, 12, 13),
           (2, 9, 13), (4, 10, 13), (4, 11, 13), (3, 10, 13)]

# Section 2 obstruction table: (n, j, T_j, [earlier witnesses T_i])
OBSTRUCTIONS = [
    (13, 112, (3, 5, 12), [(2, 5, 9)]),
    (13, 114, (4, 5, 12), [(2, 5, 9), (2, 5, 10)]),
    (14, 146, (3, 5, 13), [(2, 5, 9)]),
    (14, 149, (4, 5, 13), [(2, 5, 9), (2, 5, 10)]),
]

# audit table of the proof of Proposition 3, plus the n=15,16 controls
PAPER_FACETS = {13: 169, 14: 238, 15: 320, 16: 416}
PAPER_PAIRS = {13: 14196, 14: 28203}
PAPER_SPANNING = {13: 1, 14: 8, 15: 16, 16: 25}
PAPER_SPHERE_DIM = {13: 9, 14: 10}          # abstract: S^9 and \/^8 S^10
PAPER_BLOCK_SIZES = (6, 10)                 # abstract: "six and ten facets"
PAPER_CS_THRESHOLD = 15                     # 6p-3 at p=3
PAPER_CONJ_START = 13                       # 4p+1 at p=3

# The in-block pair audit that this program PUBLISHES in its own header:
# (pairs internal to the window, pairs discharged from in-block predecessors
# only, 1-based positions of the residual pairs), and the number of genuine
# facet witnesses available to each residual pair.  These are gated by D6/D8
# below.  Before the repair, D6 compared the pair count with C(k,2), which the
# counting loop that produces it makes true for EVERY possible input, and the
# discharge count "15 of 15 / 44 of 45" was printed but never gated -- so a
# reordering of the printed block that discharged only 13 of 15 pairs passed
# the whole program (verified: it did).
PAPER_INBLOCK = {13: (15, 15, []),
                 14: (45, 44, [(147, 154)])}
PAPER_RESIDUAL_WITNESSES = {13: {}, 14: {(147, 154): 4}}

CENSUS_MAX = 32       # independent shelling census range: n = 13..CENSUS_MAX
EULER_MAX = 22        # f-vector / Euler characteristic range
HOMOLOGY_MAX = 20     # GF(2) reduced-homology range


def verdict():
    n = len(CHECKS)
    bad = [c for c, ok in CHECKS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


#############################################################################
# SECTION 1.  The graph C_n^p and the facets of Delta_p'(G) with p' = 3.
#############################################################################

def circ_dist(n, a, b):
    """d_n(a,b) = min(|a-b|, n-|a-b|), exact integers."""
    d = abs(a - b)
    return d if d <= n - d else n - d


def adjacent(n, p, a, b):
    """a ~ b in C_n^p exactly when 0 < d_n(a,b) <= p."""
    d = circ_dist(n, a, b)
    return 0 < d <= p


def triple_is_disconnected_by_edgecount(n, p, T):
    """Route 1 (the paper's criterion): a 3-vertex graph is disconnected
    exactly when it has at most one edge."""
    e = 0
    for a, b in itertools.combinations(T, 2):
        if adjacent(n, p, a, b):
            e += 1
    return e <= 1


def triple_is_disconnected_by_search(n, p, T):
    """Route 2 (the definition): flood-fill the induced subgraph G[T] and ask
    whether it has more than one component.  Independent of route 1."""
    seen = {T[0]}
    stack = [T[0]]
    while stack:
        v = stack.pop()
        for w in T:
            if w not in seen and adjacent(n, p, v, w):
                seen.add(w)
                stack.append(w)
    return len(seen) < len(set(T))


def facet_triples(n, p=3):
    """The omitted triples of Delta_3(C_n^p), i.e. the disconnected 3-sets,
    as sorted tuples in lexicographic order."""
    out = []
    for T in itertools.combinations(range(n), 3):
        if triple_is_disconnected_by_edgecount(n, p, T):
            out.append(T)
    return out


def exchange_candidates(n, T, x):
    """The three triples (T \\ {y}) union {x}, y in T, for x not in T.
    F_r cap F_j = F_j \\ {x} holds exactly for T_r among these."""
    return [tuple(sorted([v for v in T if v != y] + [x])) for y in T]


class Complex(object):
    """Precomputed tables for Delta_3(C_n^p): facet triples, membership,
    bitmasks, and the exchange-candidate lists used by Lemma 2."""

    def __init__(self, n, p=3):
        self.n = n
        self.p = p
        self.facets = facet_triples(n, p)
        self.fset = set(self.facets)
        self.full = (1 << n) - 1
        self.mask = dict((T, (1 << T[0]) | (1 << T[1]) | (1 << T[2]))
                         for T in self.facets)
        self.exch = {}
        for T in self.facets:
            lst = []
            for x in range(n):
                if x in T:
                    continue
                for c in exchange_candidates(n, T, x):
                    if c in self.fset:
                        lst.append((x, c))
            self.exch[T] = lst

    def r_mask(self, T, placed):
        """R_j of Lemma 2 as a bitmask, given the set of earlier triples."""
        r = 0
        for x, c in self.exch[T]:
            if c in placed:
                r |= 1 << x
        return r


def check_exchange(cx, order):
    """Lemma 2 condition (2), evaluated for EVERY pair i<j.

    (T_i \\ T_j) cap R_j != empty, in bitmask form
        (mask(T_i) & ~mask(T_j) & R_j) != 0.
    Returns (n_failures, n_pairs, spanning_triples, first_failures)."""
    placed = set()
    prev = []
    fails = 0
    pairs = 0
    span = []
    examples = []
    for j, T in enumerate(order):
        rj = cx.r_mask(T, placed)
        mj = cx.mask[T]
        notj = cx.full & ~mj
        for i, mi in enumerate(prev):
            pairs += 1
            if (mi & notj & rj) == 0:
                fails += 1
                if len(examples) < 3:
                    examples.append((i + 1, order[i], j + 1, T))
        if j > 0 and rj == notj:
            span.append(T)
        placed.add(T)
        prev.append(mj)
    return fails, pairs, span, examples


def check_purity(cx, order):
    """Check (3) of the proof of Proposition 3, computed independently of
    check_exchange: for each j, every inclusion-maximal member of
    {F_i cap F_j : i<j} must have cardinality |F_j| - 1 = n - 4.

    F_i cap F_j = V \\ (T_i union T_j), so inclusion-maximal intersections
    correspond to inclusion-MINIMAL unions.  Returns (n_failures, examples)."""
    n = cx.n
    want = n - 4
    masks = [cx.mask[T] for T in order]
    fails = 0
    examples = []
    for j in range(1, len(order)):
        mj = masks[j]
        fam = sorted(set(masks[i] | mj for i in range(j)))
        for a in fam:
            minimal = True
            for b in fam:
                if b != a and (b & a) == b:
                    minimal = False
                    break
            if minimal:
                card = n - bin(a).count("1")
                if card != want:
                    fails += 1
                    if len(examples) < 3:
                        examples.append((j + 1, a, card))
    return fails, examples


#############################################################################
# SECTION 2.  An INDEPENDENT shelling order, built here from scratch.
# Nothing from the paper's order L_n enters this construction: the paper never
# prints L_n, so it cannot be rebuilt from the text (see NOT REPRODUCED).
#############################################################################

def build_shelling(cx, defer=()):
    """Greedily extend a partial order by any facet that still satisfies
    Lemma 2 condition (2) against ALL facets already placed.  Facets listed in
    `defer` are considered only after every other facet has been placed, which
    is what lets a prescribed set be driven into the spanning position.

    Returns the order, or None if the greedy run gets stuck."""
    dset = set(tuple(sorted(t)) for t in defer)
    pool = [T for T in cx.facets if T not in dset]
    tail = [T for T in cx.facets if T in dset]
    order = []
    placed = set()
    prev = []
    for stage in (pool, tail):
        rem = list(stage)
        while rem:
            chosen = None
            for T in rem:
                rj = cx.r_mask(T, placed)
                notj = cx.full & ~cx.mask[T]
                ok = True
                for mi in prev:
                    if (mi & notj & rj) == 0:
                        ok = False
                        break
                if ok:
                    chosen = T
                    break
            if chosen is None:
                return None
            rem.remove(chosen)
            order.append(chosen)
            placed.add(chosen)
            prev.append(cx.mask[chosen])
    return order


def spanning_necessary_candidates(cx):
    """Order-INDEPENDENT necessary condition for a facet to be spanning: for
    every x outside T some (T \\ {y}) union {x} must itself be a facet, else
    x can never enter R_j in any order.  Returns the set of triples passing."""
    out = set()
    for T in cx.facets:
        good = True
        for x in range(cx.n):
            if x in T:
                continue
            if not any(c in cx.fset for c in exchange_candidates(cx.n, T, x)):
                good = False
                break
        if good:
            out.add(T)
    return out


#############################################################################
# SECTION 3.  The whole complex: f-vector, reduced Euler characteristic, and
# reduced homology over GF(2).  These use only the DEFINITION of the complex,
# so they are independent of any shelling and of the paper's certificate.
#############################################################################

def face_flags(n, p=3):
    """isface[S] = 1 iff the vertex set with bitmask S is a face of
    Delta_3(C_n^p), i.e. S is contained in some facet V \\ T.  Equivalently
    the complement of S contains a disconnected triple, which we compute by
    an upward sweep (U contains one iff U is one, or U minus a vertex does)."""
    full = (1 << n) - 1
    has = bytearray(1 << n)
    for T in facet_triples(n, p):
        has[(1 << T[0]) | (1 << T[1]) | (1 << T[2])] = 1
    for u in range(1 << n):
        if has[u]:
            continue
        v = u
        while v:
            b = v & -v
            v ^= b
            if has[u ^ b]:
                has[u] = 1
                break
    isface = bytearray(1 << n)
    for s in range(1 << n):
        if has[full ^ s]:
            isface[s] = 1
    return isface


def f_vector(n, p=3):
    """f[k] = number of faces with exactly k vertices (f[0] = 1, the empty
    face).  Also returns the face-flag table for reuse."""
    isface = face_flags(n, p)
    f = [0] * (n + 1)
    for s in range(1 << n):
        if isface[s]:
            f[bin(s).count("1")] += 1
    return f, isface


def maximal_face_sizes(n, isface):
    """Cardinality census {k: count} of the INCLUSION-MAXIMAL faces, read off
    the face table by brute force over all 2^n vertex subsets.

    This is a third route to the facet set, independent of both
    facet_triples() and of any order: a complex that is pure of dimension n-4
    with m facets must return exactly {n-3: m}.  It replaces the earlier form
    of check C8, which computed n - len(set(T)) over triples emitted by
    itertools.combinations(range(n), 3) -- where len(set(T)) is 3 by
    construction, so that check could not fail for any input whatsoever."""
    full = (1 << n) - 1
    out = {}
    for s in range(1 << n):
        if not isface[s]:
            continue
        u = full & ~s
        ismax = True
        while u:
            b = u & -u
            u ^= b
            if isface[s | b]:
                ismax = False
                break
        if ismax:
            k = bin(s).count("1")
            out[k] = out.get(k, 0) + 1
    return out


def binom(a, b):
    """Exact integer binomial coefficient, no floating point."""
    if b < 0 or b > a:
        return 0
    r = 1
    for i in range(b):
        r = r * (a - i) // (i + 1)
    return r


def h_vector(f):
    """h_k = sum_{i<=k} (-1)^(k-i) C(d+1-i, k-i) f_i, where f_i counts the
    faces with exactly i vertices (so f_0 = 1).  For a shellable complex the
    h-vector is the distribution of the restriction sizes |R_j|, and h_{d+1}
    is the number of spanning facets."""
    d = max(k for k in range(len(f)) if f[k] > 0) - 1
    return [sum((-1) ** (k - i) * binom(d + 1 - i, k - i) * f[i]
                for i in range(k + 1)) for k in range(d + 2)]


def reduced_euler(f):
    """chi-tilde = -(sum_{faces S, including the empty face} (-1)^{|S|}).
    Exact integer arithmetic.  For a d-sphere this is (-1)^d."""
    return -sum((-1) ** k * f[k] for k in range(len(f)))


def gf2_reduced_betti(n, p=3):
    """Reduced Betti numbers of Delta_3(C_n^p) over GF(2), from the augmented
    simplicial chain complex (the empty face is the unique (-1)-cell).  Ranks
    are computed exactly by GF(2) elimination on rows held as Python integers.

    Returns (dim, f_vector, ranks, betti) with
        betti[k] = f[k+1] - rank(d_k) - rank(d_{k+1}).
    Over GF(2) rank(d) <= rank_Q(d), hence dim H_k(GF(2)) >= dim H_k(Q);
    a vanishing GF(2) group therefore forces the rational group to vanish."""
    f, isface = f_vector(n, p)
    d = max(k for k in range(n + 1) if f[k] > 0) - 1
    bysize = [[] for _ in range(n + 2)]
    for s in range(1 << n):
        if isface[s]:
            bysize[bin(s).count("1")].append(s)
    index = [dict((s, i) for i, s in enumerate(bysize[k]))
             for k in range(n + 2)]
    ranks = {}
    for k in range(0, d + 2):
        rows = bysize[k + 1]
        if not rows:
            ranks[k] = 0
            continue
        low = index[k]
        piv = {}
        for s in rows:
            r = 0
            u = s
            while u:
                b = u & -u
                u ^= b
                r |= 1 << low[s ^ b]
            while r:
                hb = r.bit_length() - 1
                if hb in piv:
                    r ^= piv[hb]
                else:
                    piv[hb] = r
                    break
        ranks[k] = len(piv)
    betti = {}
    for k in range(0, d + 1):
        betti[k] = f[k + 1] - ranks[k] - ranks.get(k + 1, 0)
    return d, f, ranks, betti


def fmt(T):
    return "{%s}" % ",".join(str(v) for v in T)


def fmt_list(L):
    return "(" + ", ".join(fmt(T) for T in L) + ")"


#############################################################################
# SECTION A.  The exhibited object: decode it, count it, print it back, and
# check it is well formed and is what the paper says it is.
#############################################################################

def section_A():
    info("=== SECTION A: the exhibited object, decoded and printed back ===")
    for n, pos, old, new, span in ((13, POS_13, BLOCK_13_OLD, BLOCK_13_NEW,
                                    SPAN_13),
                                   (14, POS_14, BLOCK_14_OLD, BLOCK_14_NEW,
                                    SPAN_14)):
        info("n=%d  L_%d positions %d-%d" % (n, n, pos[0], pos[1]))
        info("   before: %s" % fmt_list(old))
        info("   after : %s" % fmt_list(new))
        info("   S_%d  : {%s}  (%d facet%s)"
             % (n, ", ".join(fmt(T) for T in span), len(span),
                "" if len(span) == 1 else "s"))

        wf = all(len(T) == 3 and T[0] < T[1] < T[2] and
                 all(0 <= v < n for v in T) for T in old + new + span)
        check("A1_wellformed_triples_n%d" % n, wf,
              "all entries are strictly increasing 3-subsets of Z/%dZ" % n)

        dist = (len(set(old)) == len(old) and len(set(new)) == len(new)
                and len(set(span)) == len(span))
        check("A2_no_repetition_n%d" % n, dist,
              "old block %d/%d distinct, new block %d/%d distinct, S_%d %d/%d "
              "distinct" % (len(set(old)), len(old), len(set(new)), len(new),
                            n, len(set(span)), len(span)))

        # the abstract's "six and ten": the block sizes, plus how many
        # positions actually move (5 of 6 at n=13, 10 of 10 at n=14).
        want = PAPER_BLOCK_SIZES[0] if n == 13 else PAPER_BLOCK_SIZES[1]
        moved = sum(1 for a, b in zip(old, new) if a != b)
        a3 = (len(old) == want and len(new) == want
              and pos[1] - pos[0] + 1 == want)
        check("A3_block_size_n%d" % n, a3,
              ("block size %d = window %d-%d of length %d; %d of %d positions "
               "change" % (len(old), pos[0], pos[1], pos[1] - pos[0] + 1,
                           moved, len(old))) if a3 else
              ("abstract's count %d, listed block sizes %d and %d, window "
               "%d-%d has length %d -- these disagree"
               % (want, len(old), len(new), pos[0], pos[1],
                  pos[1] - pos[0] + 1)))

        a4 = sorted(old) == sorted(new)
        check("A4_new_block_is_permutation_of_old_n%d" % n, a4,
              ("multiset equality: the repair reorders, it does not change the "
               "facet set") if a4 else
              ("the replacement is NOT a permutation: only in old %s, only in "
               "new %s" % (fmt_list(sorted(set(old) - set(new))),
                           fmt_list(sorted(set(new) - set(old))))))

        # The length of L_n is DERIVED from the definition here, not read from
        # the paper's audit table; comparing the window against the paper's own
        # "169" would have been the paper checked against itself.
        m = len(facet_triples(n, P))
        a5 = 1 <= pos[0] <= pos[1] <= m and m == PAPER_FACETS[n]
        check("A5_window_inside_order_n%d" % n, a5,
              ("positions %d-%d lie in 1..%d, where %d is the facet count "
               "derived from the definition (the audit table also says %d)"
               % (pos[0], pos[1], m, m, PAPER_FACETS[n])) if a5 else
              ("positions %d-%d against a derived order length of %d (audit "
               "table says %d)" % (pos[0], pos[1], m, PAPER_FACETS[n])))

        check("A6_spanning_disjoint_from_block_n%d" % n,
              not (set(span) & set(new)),
              "no member of S_%d lies in the edited window, so its spanning "
              "status is untouched by the reordering" % n)


#############################################################################
# SECTION B.  Every hypothesis of the statement being acted on.
#############################################################################

def section_B():
    info("=== SECTION B: hypotheses of Conjecture 4.2(ii) at p=3 ===")

    check("B1_conjecture_range_start", 4 * P + 1 == PAPER_CONJ_START,
          "4p+1 = %d at p=%d, and the two exhibited cases are n=13,14 >= %d"
          % (4 * P + 1, P, 4 * P + 1))
    check("B2_cited_theorem_threshold", 6 * P - 3 == PAPER_CS_THRESHOLD,
          "6p-3 = %d at p=%d (paper says %d), so the range left open by "
          "[CS, Thm 1.4] is exactly {13,14}"
          % (6 * P - 3, P, PAPER_CS_THRESHOLD))
    check("B3_gap_is_two_cases",
          sorted(set(range(4 * P + 1, 6 * P - 3))) == [13, 14],
          "{n : 4p+1 <= n < 6p-3} = {13,14}")

    bad = [n for n in range(13, 201)
           if (n * n - 4 * n * P - n + 2) != (n * n - 13 * n + 2)]
    check("B4_exponent_reduces_at_p3", not bad,
          "(n^2-4np-n+2)/2 = (n^2-13n+2)/2 for n=13..200, exact integers")
    bad = [n for n in range(13, 201) if (n * n - 13 * n + 2) % 2 != 0]
    check("B5_exponent_is_an_integer", not bad,
          "n^2-13n+2 is even for every n=13..200")
    bad = [n for n in range(13, 201)
           if (n - 6) * (n - 7) // 2 - 20 != (n * n - 13 * n + 2) // 2]
    check("B6_binom_identity", not bad,
          "C(n-6,2)-20 = (n^2-13n+2)/2 for n=13..200; values 1,8,16,25 at "
          "n=13,14,15,16")

    # the two independent readings of "G[T] is disconnected" must agree
    mism = []
    for n in range(13, 19):
        for T in itertools.combinations(range(n), 3):
            a = triple_is_disconnected_by_edgecount(n, P, T)
            b = triple_is_disconnected_by_search(n, P, T)
            if a != b:
                mism.append((n, T))
    check("B7_disconnected_triple_criterion", not mism,
          "'at most one edge' agrees with explicit connectivity of G[T] on "
          "all %d triples for n=13..18"
          % sum(binom(n, 3) for n in range(13, 19))
          if not mism else "disagreements: %s" % mism[:3])

    # a non-facet must be rejected: {0,1,4} at n=13 induces two edges
    check("B8_facet_test_discriminates",
          not triple_is_disconnected_by_edgecount(13, P, (0, 1, 4))
          and triple_is_disconnected_by_edgecount(13, P, (0, 1, 5)),
          "{0,1,4} (2 edges) rejected, {0,1,5} (1 edge) accepted at n=13")

    for n, old, new, span in ((13, BLOCK_13_OLD, BLOCK_13_NEW, SPAN_13),
                              (14, BLOCK_14_OLD, BLOCK_14_NEW, SPAN_14)):
        cx = Complex(n, P)
        bad = [T for T in old + new if T not in cx.fset]
        check("B9_block_triples_are_facets_n%d" % n, not bad,
              "all %d distinct block triples induce at most one edge of C_%d^3"
              % (len(set(old)), n) if not bad else "not facets: %s"
              % fmt_list(bad))
        bad = [T for T in span if T not in cx.fset]
        check("B10_spanning_triples_are_facets_n%d" % n, not bad,
              "all %d member(s) of S_%d are genuine omitted triples"
              % (len(span), n) if not bad else "not facets: %s" % fmt_list(bad))


#############################################################################
# SECTION C.  The conclusion, COMPUTED: Delta_3(C_n^3) is shellable, with the
# spanning counts the paper claims.  The shelling used here is built by this
# program, because the paper's L_n is not printed anywhere (NOT REPRODUCED).
#############################################################################

def section_C(orders):
    info("=== SECTION C: shellability and spanning counts, computed ===")
    for n in (13, 14, 15, 16):
        cx = Complex(n, P)
        m = len(cx.facets)
        check("C1_facet_count_n%d" % n, m == PAPER_FACETS[n],
              "%d disconnected triples; audit table says %d"
              % (m, PAPER_FACETS[n]))

        order = build_shelling(cx)
        check("C2_shelling_constructed_n%d" % n,
              order is not None and len(order) == m,
              "independent greedy order of length %s"
              % (len(order) if order else "None"))
        if order is None:
            continue
        orders[n] = order

        # check (1) of the proof of Proposition 3: the list is exactly the set
        # of disconnected triples, without omission or repetition
        check("C3_list_is_exactly_the_facet_set_n%d" % n,
              len(set(order)) == m and set(order) == cx.fset,
              "%d listed, %d distinct, set equality with the enumerated "
              "facet set" % (len(order), len(set(order))))

        fails, pairs, span, ex = check_exchange(cx, order)
        if n in PAPER_PAIRS:
            check("C4_pair_count_n%d" % n, pairs == PAPER_PAIRS[n],
                  "C(%d,2) = %d pairs examined; paper says %d"
                  % (m, pairs, PAPER_PAIRS[n]))
        check("C5_exchange_condition_all_pairs_n%d" % n, fails == 0,
              "condition (2) holds for all %d pairs i<j (0 failures)" % pairs
              if fails == 0 else "%d failures, e.g. %s" % (fails, ex))

        pfails, pex = check_purity(cx, order)
        check("C6_purity_condition_n%d" % n, pfails == 0,
              "every inclusion-maximal F_i cap F_j has cardinality "
              "|F_j|-1 = %d (0 failures)" % (n - 4)
              if pfails == 0 else "%d failures, e.g. %s" % (pfails, pex))

        check("C7_spanning_count_n%d" % n, len(span) == PAPER_SPANNING[n],
              "%d spanning facets; paper says %d"
              % (len(span), PAPER_SPANNING[n]))

        # Purity read off the FACE table over all 2^n subsets.  The previous
        # form of this check computed n - len(set(T)) over triples emitted by
        # combinations(range(n), 3), where len(set(T)) is 3 by construction, so
        # it could not fail for any input; this form is a third, order-free
        # route to the facet set and reproduces its cardinality.
        f, isface = f_vector(n, P)
        maxfaces = maximal_face_sizes(n, isface)
        sizes = set(maxfaces)
        dim = (max(sizes) - 1) if sizes else -1
        nmax = sum(maxfaces.values())
        c8 = sizes == set([n - 3]) and dim == n - 4 and nmax == m
        check("C8_pure_of_dimension_n_minus_4_n%d" % n, c8,
              ("all %d inclusion-maximal faces found by sweeping the 2^%d "
               "vertex subsets have %d = n-3 vertices, so the complex is pure "
               "of dim %d = n-4, and their count reproduces the %d facets by a "
               "third independent route" % (nmax, n, n - 3, dim, m)) if c8 else
              ("inclusion-maximal faces have cardinalities %s (want {%d}) and "
               "number %d (facets %d): the complex is not pure of dimension "
               "n-4" % (sorted(maxfaces.items()), n - 3, nmax, m)))
        if n in PAPER_SPHERE_DIM:
            check("C9_sphere_dimension_n%d" % n, dim == PAPER_SPHERE_DIM[n],
                  "abstract asserts S^%d; computed dim = %d"
                  % (PAPER_SPHERE_DIM[n], dim))

        # the restriction sizes |R_j| of the order must reproduce the
        # h-vector computed from the f-vector of the complex, an identity
        # linking the order to the complex that a wrong order violates
        # (f was computed above, from the face table)
        hv = h_vector(f)
        dist = [0] * len(hv)
        placed = set()
        for T in order:
            s = bin(cx.r_mask(T, placed)).count("1")
            if s < len(dist):
                dist[s] += 1
            placed.add(T)
        c12 = dist == hv and hv[-1] == PAPER_SPANNING[n]
        check("C12_restriction_sizes_match_h_vector_n%d" % n, c12,
              ("|R_j| distribution %s equals the h-vector of the complex, "
               "whose last entry h_%d = %d is the spanning count"
               % (hv, dim + 1, hv[-1])) if c12 else
              ("distribution %s vs h-vector %s; h_%d = %d, paper claims %d"
               % (dist, hv, dim + 1, hv[-1], PAPER_SPANNING[n])))


def section_C2():
    """The paper's spanning SETS, not just their sizes."""
    info("=== SECTION C2: the claimed spanning sets S_13, S_14 ===")
    for n, span in ((13, SPAN_13), (14, SPAN_14)):
        cx = Complex(n, P)
        cands = spanning_necessary_candidates(cx)
        tgt = set(tuple(sorted(t)) for t in span)

        f, _ = f_vector(n, P)
        hv = h_vector(f)
        check("C13_spanning_set_size_n%d" % n,
              len(span) == PAPER_SPANNING[n] == hv[-1],
              "|S_%d| = %d, the paper's count is %d, and h_top of the complex "
              "is %d -- the number is forced by the complex"
              % (n, len(span), PAPER_SPANNING[n], hv[-1]))

        c10 = tgt <= cands
        check("C10_spanning_necessary_condition_n%d" % n, c10,
              ("every member of S_%d passes the order-independent test 'for "
               "each x outside T some (T\\{y})u{x} is a facet'; %d of %d "
               "facets pass it" % (n, len(cands), len(cx.facets))) if c10 else
              ("these members of S_%d can never be spanning in any order: %s"
               % (n, fmt_list(sorted(tgt - cands)))))

        order = build_shelling(cx, defer=span)
        ok = order is not None
        got = set()
        fails = pairs = -1
        if ok:
            fails, pairs, sp, _ = check_exchange(cx, order)
            got = set(sp)
        check("C11_S_n_realized_by_a_verified_shelling_n%d" % n,
              ok and fails == 0 and got == tgt,
              "a shelling order verified on all %d pairs has spanning set "
              "exactly S_%d" % (pairs, n)
              if ok and fails == 0 and got == tgt
              else "constructed=%s fails=%s got=%s"
                   % (ok, fails, sorted(got)))
    info("C11 shows S_13 and S_14 ARE the spanning set of some verified "
         "shelling; it does not certify that they are the spanning set of the "
         "paper's own order L-hat_n, which the paper never prints.")


#############################################################################
# SECTION D.  The repair mechanism: the obstruction table, and what the two
# edited blocks do about it.  Everything here is computed from the block
# listings alone, so it needs no knowledge of the unprinted prefix of L_n.
#############################################################################

def inblock_witnesses(cx, block, j, x):
    """Members of block[:j] of the form (block[j] \\ {y}) union {x}."""
    Tj = block[j]
    early = set(block[:j])
    return [c for c in exchange_candidates(cx.n, Tj, x)
            if c in cx.fset and c in early]


def section_D():
    info("=== SECTION D: the four obstructions and their repair ===")
    for n, pos, Tj, wits in OBSTRUCTIONS:
        cx = Complex(n, P)
        old, new, start = ((BLOCK_13_OLD, BLOCK_13_NEW, POS_13[0]) if n == 13
                           else (BLOCK_14_OLD, BLOCK_14_NEW, POS_14[0]))
        k = pos - start
        tag = "n%d_j%d" % (n, pos)

        d1 = 0 <= k < len(old) and old[k] == Tj
        check("D1_obstruction_position_%s" % tag, d1,
              ("position %d of L_%d is entry %d of the original block, and it "
               "holds %s as the table says" % (pos, n, k + 1, fmt(Tj))) if d1
              else ("the table puts %s at position %d, but entry %d of the "
                    "original block is %s"
                    % (fmt(Tj), pos, k + 1,
                       fmt(old[k]) if 0 <= k < len(old) else "out of range")))

        badw = [w for w in wits if w not in cx.fset]
        check("D2_obstruction_witnesses_are_facets_%s" % tag, not badw,
              ("witnesses %s are genuine omitted triples"
               % ", ".join(fmt(w) for w in wits)) if not badw else
              ("listed witnesses that are NOT facets: %s" % fmt_list(badw)))

        jn = new.index(Tj) if Tj in new else -1
        jo = old.index(Tj) if Tj in old else -1
        for w in wits:
            diff = sorted(set(w) - set(Tj))
            dead = [x for x in diff
                    if not any(c in cx.fset
                               for c in exchange_candidates(n, Tj, x))]
            live = [x for x in diff if x not in dead]
            wtag = "%s_Ti%s" % (tag, "".join(str(v) for v in w))

            d3 = len(dead) == 1 and len(live) == 1
            check("D3_dead_exchange_direction_%s" % wtag, d3,
                  ("T_i\\T_j = %s; vertex %s has NO facet of the form "
                   "(T_j\\{y})u{x}, so it can never enter R_j in any order; "
                   "the pair can only be discharged through %s"
                   % (diff, dead, live)) if d3 else
                  ("T_i\\T_j = %s does not split as one dead and one live "
                   "vertex: dead=%s live=%s" % (diff, dead, live)))

            wnew = [c for x in live for c in inblock_witnesses(cx, new, jn, x)]
            wold = [c for x in live for c in inblock_witnesses(cx, old, jo, x)]
            check("D4_repaired_block_discharges_%s" % wtag, bool(wnew),
                  ("in the repaired order %s precedes %s inside the block and "
                   "supplies %s" % (", ".join(fmt(c) for c in wnew), fmt(Tj),
                                    live)) if wnew else
                  ("NO in-block predecessor of %s in the repaired order "
                   "supplies %s, so the listed obstruction is not repaired "
                   "locally" % (fmt(Tj), live)))
            check("D5_original_block_does_not_discharge_%s" % wtag, not wold,
                  ("in the original order no in-block predecessor of %s has "
                   "the required form, which is exactly the obstruction the "
                   "table reports" % fmt(Tj)) if not wold else
                  ("the ORIGINAL block already supplies %s from %s, so the "
                   "table's obstruction is not what it says"
                   % (live, ", ".join(fmt(c) for c in wold))))

            # D5 is SILENT whenever T_j sits at the first position of the
            # original block: old[:0] is empty, so `not wold` holds for every
            # possible input, and two of the four listed obstructions (n=13
            # j=112 and n=14 j=146) are exactly of that kind.  The mechanism
            # the paper actually claims is checkable and cannot hold vacuously:
            # the discharging facet lies INSIDE the block, after T_j before the
            # repair and before T_j after it.
            mech = [(c, old.index(c), new.index(c))
                    for x in live
                    for c in exchange_candidates(n, Tj, x)
                    if c in cx.fset and c in old and c in new
                    and old.index(c) > jo and new.index(c) < jn]
            check("D5b_repair_moves_the_witness_across_%s" % wtag, bool(mech),
                  ("%s is in the block at position %d, after %s (position %d) "
                   "before the repair and before %s (position %d) after it: "
                   "the reordering, and nothing else, discharges the pair"
                   % (fmt(mech[0][0]), start + mech[0][1], fmt(Tj), start + jo,
                      fmt(Tj), start + jn)) if mech else
                  ("no facet (T_j\\{y})u{x} with x in %s moves from after %s "
                   "to before it inside the block, so the reordering does not "
                   "explain the repair of this pair" % (live, fmt(Tj))))


def inblock_audit(cx, block):
    """(pairs, pairs discharged from in-block predecessors, residual pairs) for
    a block, using in-block data only.  Factored out so that the negative
    control F9 can re-run it on deliberately reordered blocks."""
    n = cx.n
    total = 0
    done = 0
    residual = []
    for j in range(len(block)):
        Tj = block[j]
        rset = set()
        for x in range(n):
            if x in Tj:
                continue
            if inblock_witnesses(cx, block, j, x):
                rset.add(x)
        for i in range(j):
            Ti = block[i]
            total += 1
            if (set(Ti) - set(Tj)) & rset:
                done += 1
            else:
                residual.append((i, Ti, j, Tj))
    return total, done, residual


def section_D2():
    """Every pair internal to a repaired block, audited with in-block data
    only.  A pair not discharged inside the block must at least have a genuine
    facet available as an exchange witness, else the certificate would be
    refuted outright no matter what the unprinted prefix of L_n contains."""
    info("=== SECTION D2: in-block pair audit of the repaired blocks ===")
    for n, new, start in ((13, BLOCK_13_NEW, POS_13[0]),
                          (14, BLOCK_14_NEW, POS_14[0])):
        cx = Complex(n, P)
        total, done, residual = inblock_audit(cx, new)

        # The DISCHARGE count, not just the pair count, is gated here.  The
        # pair count alone is produced by the very loop that counts it, so
        # "total == C(k,2)" is true for every conceivable block and gates
        # nothing; the number that carries information is `done`.
        want_total, want_done, want_res = PAPER_INBLOCK[n]
        d6 = (total == len(new) * (len(new) - 1) // 2
              and total == want_total and done == want_done)
        check("D6_inblock_pairs_discharged_n%d" % n, d6,
              ("%d of %d pairs internal to the %d-position window are "
               "discharged from in-block predecessors alone, the count this "
               "program publishes for the printed block"
               % (done, total, len(new))) if d6 else
              ("%d of %d pairs discharged in-block; the printed block of "
               "equation (%d) gives %d of %d -- the content or the ORDER of "
               "the window has changed"
               % (done, total, 3 if n == 13 else 4, want_done, want_total)))

        blocked = []
        witcount = {}
        for i, Ti, j, Tj in residual:
            diff = sorted(set(Ti) - set(Tj))
            later = set(new[j:])          # known to sit at or after position j
            cands = [c for x in diff
                     for c in exchange_candidates(n, Tj, x)
                     if c in cx.fset and c not in later]
            info("   residual pair: position %d %s vs position %d %s, "
                 "T_i\\T_j = %s, %d facet exchange candidate(s) available "
                 "from the unprinted prefix: %s"
                 % (start + i, fmt(Ti), start + j, fmt(Tj), diff, len(cands),
                    ", ".join(fmt(c) for c in cands)))
            witcount[(start + i, start + j)] = len(cands)
            if not cands:
                blocked.append((start + i, Ti, start + j, Tj))

        check("D7_no_structurally_impossible_inblock_pair_n%d" % n,
              not blocked,
              "each of the %d residual in-block pair(s) still has at least one "
              "genuine facet that could discharge it from the prefix"
              % len(residual)
              if not blocked else "IMPOSSIBLE pairs: %s" % blocked)

        # The IDENTITY of the residual pairs and the NUMBER of witnesses each
        # retains are published in this program's header ("the single residual
        # pair, position 147 against position 154 at n=14, ... four genuine
        # facet witnesses"); they were printed by info() and gated by nothing.
        want_wit = PAPER_RESIDUAL_WITNESSES[n]
        d8 = (sorted(witcount) == sorted(want_res)
              and all(witcount[k] == want_wit[k] for k in witcount))
        check("D8_residual_pairs_are_the_published_ones_n%d" % n, d8,
              ("the residual in-block pairs are exactly %s, each with the "
               "published number of surviving facet witnesses %s"
               % (sorted(witcount) or "none",
                  sorted(want_wit.items()) or "-")) if d8 else
              ("residual pairs %s with witness counts %s; published: pairs %s "
               "with counts %s"
               % (sorted(witcount), sorted(witcount.items()),
                  sorted(want_res), sorted(want_wit.items()))))


#############################################################################
# SECTION E.  The homotopy claim, cross-checked from the DEFINITION of the
# complex alone: f-vector, reduced Euler characteristic, and reduced homology
# over GF(2).  No shelling and no certificate is used here.
#############################################################################

def section_E():
    info("=== SECTION E: Euler characteristic and GF(2) homology ===")
    for n in range(13, EULER_MAX + 1):
        f, _ = f_vector(n, P)
        chi = reduced_euler(f)
        beta = (n * n - 13 * n + 2) // 2
        dim = max(k for k in range(n + 1) if f[k] > 0) - 1
        sign = -1 if dim % 2 else 1
        ok = abs(chi) == beta and chi == sign * beta
        det = ("chi-tilde = %+d, dim = %d, so |chi-tilde| = %d = "
               "(n^2-13n+2)/2 with the sign (-1)^(n-4)" % (chi, dim, beta))
        if n in PAPER_SPANNING:
            if abs(chi) != PAPER_SPANNING[n]:
                ok = False
                det = ("|chi-tilde| = %d but the paper's spanning count is %d"
                       % (abs(chi), PAPER_SPANNING[n]))
            else:
                det += "; paper's spanning count %d" % PAPER_SPANNING[n]
        if abs(chi) != beta:
            ok = False
            det = ("|chi-tilde| = %d but (n^2-13n+2)/2 = %d"
                   % (abs(chi), beta))
        check("E1_reduced_euler_forces_multiplicity_n%d" % n, ok, det)
        if n == 13:
            info("   f-vector n=13: %s" % f)
        if n == 14:
            info("   f-vector n=14: %s" % f)

    for n in range(13, HOMOLOGY_MAX + 1):
        dim, f, ranks, betti = gf2_reduced_betti(n, P)
        beta = (n * n - 13 * n + 2) // 2
        top = betti.get(dim, None)
        lower = dict((k, v) for k, v in betti.items() if k != dim and v != 0)
        ok = (top == beta) and not lower
        det = ("reduced H_k(.;GF(2)) = 0 for k != %d and = GF(2)^%d for k = %d"
               % (dim, top, dim))
        if n in PAPER_SPANNING:
            det += ("; matches %s^%d with multiplicity %d"
                    % ("S", dim, PAPER_SPANNING[n]))
            ok = ok and top == PAPER_SPANNING[n]
        if not ok:
            det = "betti = %s, expected %d in degree %d" % (betti, beta, dim)
        check("E2_gf2_reduced_homology_n%d" % n, ok, det)
    info("Over GF(2) rank(d) <= rank_Q(d), so a vanishing GF(2) homology "
         "group forces the rational one to vanish; with the Euler "
         "characteristic this pins the rational homology of Delta_3(C_n^3) "
         "exactly, and by universal coefficients excludes all 2-torsion.")


#############################################################################
# SECTION F.  Negative controls.  Each of these deliberately corrupts an input
# and requires the corresponding test above to REJECT it, so that no check in
# this program can be passing vacuously.
#############################################################################

def section_F(orders):
    info("=== SECTION F: negative controls (every guard's failing path) ===")
    missing = [n for n in (13, 14) if n not in orders]
    if missing:
        # never silently skip: a missing order is itself a failure, and the
        # VERDICT line must still be printed rather than lost to a KeyError
        check("F0_orders_available_for_the_controls", False,
              "section C built no shelling for n=%s, so the negative controls "
              "cannot run" % missing)
        return
    for n in (13, 14):
        cx = Complex(n, P)
        order = orders[n]

        rev = list(reversed(order))
        f1, p1, s1, _ = check_exchange(cx, rev)
        check("F1_reversed_order_rejected_n%d" % n, f1 > 0,
              "reversing the verified shelling produces %d exchange failures "
              "out of %d pairs" % (f1, p1))
        pf1, _ = check_purity(cx, rev)
        check("F2_reversed_order_fails_purity_n%d" % n, pf1 > 0,
              "the same reversed order produces %d purity failures" % pf1)

        rot = [order[-1]] + order[:-1]
        f2, p2, _, _ = check_exchange(cx, rot)
        check("F3_rotated_order_rejected_n%d" % n, f2 > 0,
              "moving the last facet to position 1 produces %d exchange "
              "failures" % f2)

        _, _, span, _ = check_exchange(cx, order)
        check("F4_spanning_test_is_not_vacuous_n%d" % n,
              0 < len(span) < len(order) and order[0] not in span,
              "%d of %d facets are spanning, and the first facet (R_1 empty) "
              "is not among them" % (len(span), len(order)))

    # corrupting the exhibited object must be caught
    cx14 = Complex(14, P)
    bad_block = list(BLOCK_14_NEW)
    bad_block[7] = (0, 1, 2)                      # a connected triple
    check("F5_non_facet_in_block_rejected",
          not all(T in cx14.fset for T in bad_block)
          and triple_is_disconnected_by_edgecount(14, P, BLOCK_14_NEW[7]),
          "replacing %s by the connected triple {0,1,2} is rejected by the "
          "facet test of B9" % fmt(BLOCK_14_NEW[7]))

    bad_block = list(BLOCK_14_NEW)
    bad_block[0] = bad_block[1]                   # no longer a permutation
    f6a = sorted(BLOCK_14_OLD) != sorted(bad_block)
    f6b = sorted(BLOCK_14_OLD) == sorted(BLOCK_14_NEW)
    check("F6_non_permutation_rejected", f6a and f6b,
          "duplicating a block entry is rejected by the multiset test of A4"
          if f6a and f6b else
          "corrupted block %srejected; unmodified block %sa permutation"
          % ("" if f6a else "NOT ", "is " if f6b else "is NOT "))

    cx13 = Complex(13, P)
    moved = [(4, 5, 11), (5, 11, 12), (5, 10, 12),
             (3, 5, 12), (4, 5, 12), (5, 9, 12)]  # witness pushed past T_j
    j = moved.index((3, 5, 12))
    f7a = not inblock_witnesses(cx13, moved, j, 9)
    f7b = (3, 5, 12) in BLOCK_13_NEW and bool(
        inblock_witnesses(cx13, BLOCK_13_NEW,
                          BLOCK_13_NEW.index((3, 5, 12)), 9))
    check("F7_broken_repair_rejected", f7a and f7b,
          "moving {5,9,12} after {3,5,12} inside the n=13 block destroys the "
          "witness, and the D4 test detects it" if f7a and f7b else
          "corrupted block %srejected; the paper's block %sthe witness"
          % ("" if f7a else "NOT ", "supplies " if f7b else "does NOT supply "))

    cands = spanning_necessary_candidates(cx13)
    outside = [T for T in cx13.facets if T not in cands]
    f8 = len(outside) > 0 and SPAN_13[0] in cands
    check("F8_spanning_necessary_condition_discriminates", f8,
          ("%d of %d facets FAIL the order-independent spanning test (e.g. "
           "%s), while S_13 = %s passes it"
           % (len(outside), len(cx13.facets),
              fmt(outside[0]) if outside else "-", fmt(SPAN_13[0]))) if f8 else
          ("%d of %d facets fail the test and S_13 = %s is %s among them"
           % (len(outside), len(cx13.facets), fmt(SPAN_13[0]),
              "" if SPAN_13[0] in cands else "itself")))

    # F9: the count D6 now gates must be sensitive to the ORDER of the printed
    # block, not only to its content.  A transposition of two entries leaves the
    # block a permutation of itself (A4 still passes) and in most cases leaves
    # the four listed obstructions intact (D3-D5b still pass), so if the
    # discharge count never moved, D6 would still be gating nothing.
    for n, blk in ((13, BLOCK_13_NEW), (14, BLOCK_14_NEW)):
        cx = Complex(n, P)
        base_done = inblock_audit(cx, blk)[1]
        worse = []
        for a, b in itertools.combinations(range(len(blk)), 2):
            mm = list(blk)
            mm[a], mm[b] = mm[b], mm[a]
            d = inblock_audit(cx, mm)[1]
            if d < base_done:
                worse.append(((a + 1, b + 1), d))
        f9 = bool(worse) and base_done == PAPER_INBLOCK[n][1]
        check("F9_inblock_discharge_count_is_order_sensitive_n%d" % n, f9,
              ("%d of the %d transpositions of the repaired block strictly "
               "lower the in-block discharge count (swapping entries %d and %d "
               "gives %d instead of %d), so D6 rejects a reordered window"
               % (len(worse), len(blk) * (len(blk) - 1) // 2,
                  worse[0][0][0], worse[0][0][1], worse[0][1], base_done))
              if f9 else
              ("no transposition of the block lowers the discharge count "
               "(base %d, published %d), so D6 cannot see a reordering"
               % (base_done, PAPER_INBLOCK[n][1])))

    # F10: the maximal-face sweep that C8 now uses must read the object rather
    # than the loop bounds -- the same code on Delta_3(C_13^4) must disagree.
    mf3 = maximal_face_sizes(13, face_flags(13, P))
    mf4 = maximal_face_sizes(13, face_flags(13, 4))
    f10 = mf3 == {10: PAPER_FACETS[13]} and mf4 != mf3
    check("F10_maximal_face_sweep_discriminates", f10,
          "the sweep returns %s for Delta_3(C_13^3) and %s for "
          "Delta_3(C_13^4): the purity route depends on the graph"
          % (mf3, mf4))


#############################################################################
# SECTION G.  Census.  Theorem 1 asserts shellability and the wedge for EVERY
# n >= 13.  We re-derive it here, independently of both the paper's order and
# of [CS, Theorem 1.4], for n = 13..CENSUS_MAX -- one shelling per n, verified
# against condition (2) on all C(m,2) pairs, against the purity condition, and
# against the multiplicity (n^2-13n+2)/2.
#############################################################################

def section_G():
    info("=== SECTION G: independent census, n = 13..%d ===" % CENSUS_MAX)
    info("   n  facets      pairs  exch.fail  purity.fail  spanning  formula")
    ok_all = True
    rows = 0
    for n in range(13, CENSUS_MAX + 1):
        cx = Complex(n, P)
        m = len(cx.facets)
        order = build_shelling(cx)
        if order is None:
            check("G1_census_n%d" % n, False, "no shelling constructed")
            ok_all = False
            continue
        fails, pairs, span, ex = check_exchange(cx, order)
        pfails, _ = check_purity(cx, order)
        beta = (n * n - 13 * n + 2) // 2
        good = (len(order) == m and set(order) == cx.fset and fails == 0
                and pfails == 0 and len(span) == beta)
        info("  %2d  %6d  %9d  %9d  %11d  %8d  %7d"
             % (n, m, pairs, fails, pfails, len(span), beta))
        rows += 1
        if not good:
            ok_all = False
            check("G1_census_n%d" % n, False,
                  "m=%d fails=%d purity=%d spanning=%d expected=%d"
                  % (m, fails, pfails, len(span), beta))
    check("G1_census_shellable_with_right_multiplicity",
          ok_all and rows == CENSUS_MAX - 12,
          "for every n = 13..%d an explicit shelling of Delta_3(C_n^3) was "
          "built and verified on all C(m,2) pairs, and its spanning count is "
          "exactly (n^2-13n+2)/2 (%d values of n)" % (CENSUS_MAX, rows))


def not_reproduced():
    info("")
    info("=== NOT REPRODUCED ===")
    info("1. The paper's OWN two orders, L-hat_13 and L-hat_14, are not")
    info("   rebuilt here, and cannot be rebuilt from the manuscript: the")
    info("   base order L_n is specified only as 'the order of Chauhan and")
    info("   Shukla, Section 3 / Definition 3.10, specialized to p=3', the")
    info("   169 and 238 triples are never printed, and the paper's 'Exact")
    info("   verification' section states that no data file and no copy of")
    info("   the program that produced the certificates accompanies it --")
    info("   that program and its transcript are held in an archive")
    info("   available on request, and are not inputs to this file.")
    info("   Consequently the audit-table entries")
    info("   'exch. fail. 0' and 'purity fail. 0' FOR THE PAPER'S ORDERS, and")
    info("   the claim that S_13, S_14 are the spanning facets OF THOSE")
    info("   ORDERS, are not re-run.  What is re-run instead: the whole of")
    info("   the repair mechanism from the printed blocks (Section D), and")
    info("   the paper's CONCLUSION -- shellability of Delta_3(C_13^3) and")
    info("   Delta_3(C_14^3) with 1 and 8 spanning facets -- via orders this")
    info("   program builds and verifies itself (Sections C, E, G).")
    info("1a. Three further identities are NOT pinned, each verified by")
    info("   deliberate substitution to survive every check in this program:")
    info("   (i) the ORDER of the printed block is pinned only up to")
    info("   permutations with the same in-block discharge structure -- e.g.")
    info("   exchanging the two obstruction targets {3,5,12} and {4,5,12} in")
    info("   equation (3) still gives 15 of 15 and still discharges all four")
    info("   listed obstructions, so no check on the printed data can see it;")
    info("   (ii) the IDENTITY of the earlier witness triples T_i in the")
    info("   Section 2 table -- replacing {2,5,9} by {2,5,11} at n=13, j=112")
    info("   yields another facet with the same dead vertex 2 and a live")
    info("   vertex that the same reordering discharges, so D2-D5b cannot")
    info("   distinguish them; the table's claim that these are the ONLY")
    info("   obstruction pairs of L_n is likewise unverifiable without L_n;")
    info("   (iii) the members of S_13 and S_14 -- see 1b.")
    info("   All three are consequences of one omission: L_n is not printed.")
    info("1b. In particular the IDENTITY of the members of S_13 and S_14 is")
    info("   not pinned by this program.  Their NUMBER is forced (Sections C")
    info("   and E: h_{d+1} = |chi-tilde| = 1 and 8), every member is a")
    info("   genuine facet passing the order-independent spanning test (C10),")
    info("   and both sets are realised exactly by a shelling verified here")
    info("   (C11) -- but a different realisable set of the same size would")
    info("   also pass C10 and C11, because which facets are spanning depends")
    info("   on the order and the paper's order is not printed.")
    info("2. The census stops at n = %d.  Theorem 1 claims every n >= 13;" %
         CENSUS_MAX)
    info("   n > %d is not machine-checked here (the paper obtains it from" %
         CENSUS_MAX)
    info("   [CS, Theorem 1.4], not from computation).  The Euler-")
    info("   characteristic check stops at n = %d and the GF(2) homology at"
         % EULER_MAX)
    info("   n = %d, both for cost: they enumerate all 2^n vertex subsets."
         % HOMOLOGY_MAX)
    info("3. Integral ODD torsion in degrees below the top is not excluded by")
    info("   the GF(2) computation of Section E; it is excluded instead by")
    info("   shellability plus [Kozlov, Thm 12.3], which Sections C and G")
    info("   establish computationally for the two boundary cases.")


def main():
    t0 = time.time()
    info("verify.py -- Delta_3(C_n^3), p=3 case of CSV Conjecture 4.2(ii)")
    info("python %s" % sys.version.split()[0])
    info("configured ranges: shelling census n <= %d, f-vector / Euler "
         "characteristic n <= %d, GF(2) homology n <= %d"
         % (CENSUS_MAX, EULER_MAX, HOMOLOGY_MAX))
    info("COST: this is a long single-threaded run, not a quick smoke test.")
    info("   Section E enumerates all 2^n vertex subsets for n up to %d; its"
         % EULER_MAX)
    info("   GF(2) elimination at n = %d reduces boundary matrices of order"
         % HOMOLOGY_MAX)
    info("   10^5 rows by 10^5 columns whose rows are single Python integers")
    info("   of about 20 KB each; and Section G evaluates condition (2) on")
    info("   all C(m,2) pairs for n up to %d, more than 5*10^7 pairs in all."
         % CENSUS_MAX)
    info("   Budget hours of CPU time and several GB of memory.  The measured")
    info("   wall clock is printed on the 'elapsed' line at the very end, and")
    info("   that line -- not any figure quoted about this program elsewhere")
    info("   -- is this program's own statement of what it cost to run.")
    orders = {}
    section_A()
    section_B()
    section_C(orders)
    section_C2()
    section_D()
    section_D2()
    section_E()
    section_F(orders)
    section_G()
    not_reproduced()
    info("elapsed %.1f s" % (time.time() - t0))
    return verdict()


if __name__ == "__main__":
    sys.exit(main())

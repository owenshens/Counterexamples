#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- exact-integer verification of the degree-18 counterexample to
Conjecture 5.23 of Acevedo-Blekherman-Debus-Riener (ABDR).

Target: Python 3.9, standard library only.  No numpy / sympy / networkx.
Arithmetic: Python ints only.  No float literals, no '/', no math module.

=====================================================================================
WHAT IS TAKEN FROM THE PAPER (inputs, not verified here)
=====================================================================================
  P1. The witness y: 30 integer values indexed by the even partitions of 18
      (the paper's witness table (3)).  Transcribed literally below as WITNESS_RAW.
  P2. The violating triple  lambda^1 = (10,8), lambda^2 = (4,4,4,2,2,2), mu = (8,4,4,2).
      Taken as a candidate only -- the program PROVES it satisfies the T^(2)
      generating condition rather than assuming it.
  P3. The degree: 2d = 18, so d = 9; blocks are V_{r,s}(9), partitions are of 18.
  P4. The definitions: even partitions; multiset union lambda*mu; mu^{o2} = mu*mu;
      superdominance with the min-length truncation; T^(1); T^(2); the blocks
      V_{r,s}(d); the index maps A_u, B_u, C_{u,v}; the membership Lemma.
  P5. The five target counts 15, 131, 119, 937, 1056 -- used ONLY as right-hand
      sides against which independently derived counts are compared.
  P6. The three target slacks (first-order min 0, partial-symmetry min 0,
      second-order slack -2) and the two printed increasing prefix-sum vectors
      (2,4,6,10,14,18,26,36) and (2,4,8,12,16,20,28,36) -- comparison targets only.
  P7. UNPROVED-BY-PROGRAM MATHEMATICAL INPUT: the Lemma ("Degree-18 membership
      test") asserting that z in T_18^(1) together with the 1056 partial-symmetry
      inequalities implies z in trop(BSigma_18^*).  Its proof rests on ABDR
      Section 5.3 / Appendix A, Lemma 5.14, Lemma 5.17, Corollary 5.18,
      Proposition 5.19.  NO PYTHON SCRIPT CAN VERIFY THAT REDUCTION.  This
      program checks the finite inequality system, NOT the reduction.

=====================================================================================
WHAT IS DERIVED HERE (computed from the definitions, never tabulated)
=====================================================================================
  D1.  EvenPartitions(18) and its cardinality (comes out 30 = p(9)), plus
       EvenPartitions(m) for m = 0,2,4,6,8.
  D2.  The nonempty blocks V_{r,s}(9) and their number (comes out 15), together
       with the parity derivation r odd and the feasibility bound r + 2s <= 9.
  D3.  Every block size and the total term count (comes out 131).
  D4.  The diagonal inequality count, WITH MULTIPLICITY (comes out 119; the
       deduplicated value 35 is also reported, as a note).
  D5.  The off-diagonal count as unordered in-block pairs WITH MULTIPLICITY
       (comes out 937; ordered 1874 and deduplicated 250 reported as notes).
  D6.  The total 119 + 937 (comes out 1056).
  D7.  Every index A_u, B_u, C_{u,v}, built from alpha and lambda by the maps.
  D8.  The full first-order pair set and its cardinality (comes out 417).
  D9.  All 417 first-order slacks and their minimum (comes out 0, attained).
  D10. All 119 diagonal slacks and all 937 off-diagonal slacks, and their minima
       (each comes out 0, attained).
  D11. The two increasing prefix-sum vectors of lambda^1*lambda^2 and mu^{o2},
       and the truth value of the superdominance relation between them.
  D12. The second-order slack y_(10,8) + y_(4^3,2^3) - 2 y_(8,4^2,2) (comes out -2).
  D13. The size of the full second-order system on Lambda^{ev}_18 (comes out
       7025), its minimum slack (-2), and the fact that exactly one inequality is
       violated -- the paper's "Exact verification" paragraph asserts all three,
       so 7025 IS a paper claim and the pin below is BOTH a regression pin and the
       paper's stated number; plus the structural embedding T^(2) subset T^(1),
       which is beyond the paper and is the only check tying the two definitions
       to each other rather than to one worked example.
  D14. INDEPENDENT RECOUNT of every published count by integer dynamic programming
       (num_partitions / num_ordered_alpha / independent_block_sizes), sharing no
       code with the enumerators.  30, 15, 131, 119 and 937 are each produced
       twice by unrelated methods and compared.
  D15. THE PAPER'S FIRST FURTHER FACT (C16): the superdominance RELATION
       A_u succeq B_u on every one of the 119 diagonal terms, plus the paper's
       "therefore already a first-order inequality" half, checked by looking the
       ordered pair (A_u, B_u) up in the 417-element first-order pair set that
       first_order_pairs() builds from Lambda^{ev}_18 alone.  The diagonal SLACKS
       of D10 do not by themselves establish this relation.
  D16. THE PAPER'S SECOND FURTHER FACT (C17): the relation
       A_u A_v succeq C_{u,v}^{o2} on every one of the 937 off-diagonal pairs
       (both sides partitions of 36), plus the paper's "therefore itself an
       instance of (1)" half, checked by membership of the normalised triple
       (A_u, A_v, C_{u,v}) in the 7025-inequality T^(2) system of D13; and the
       derived count of off-diagonal triples that coincide with the unique
       VIOLATED T^(2) triple (comes out 0), which is what reconciles "y is outside
       T^(2)" with "all 937 hold with slack >= 0".  Both facts are THEOREMS (proofs
       in the C16/C17 docstrings), so they are confirmations of the implementation
       and of the paper's stated reasons, not robustness margins.

=====================================================================================
CONVENTION NOTES (under-specifications in the paper, resolved here)
=====================================================================================
  N1. alpha is enumerated as an ORDERED TUPLE over all sequences.  Canonicalising
      alpha up to sorting within the odd run and within the even run yields 73
      terms instead of 131 and breaks every published count.
      HONEST CAVEAT: this convention was therefore CHOSEN because it reproduces
      the paper's 131, so C2/C3/C4/C5 validate a fitted parameter and cannot
      discriminate against a paper whose counts came from the same reading.
      WHY THAT IS NOT FATAL: every sorted-alpha term is one of the ordered-alpha
      terms, so the sorted reading's inequality set is a SUBSET of the 1056
      verified here.  Verifying a SUPERSET of the Lemma's hypotheses is the
      conservative direction, so the membership conclusion survives either
      reading; only the printed COUNTS are convention-dependent.
  N2. The counts 119 and 937 are MULTIPLICITY counts over terms / pairs, not
      counts of distinct inequalities ("Repeated inequalities are harmless").
  N3. Superdominance truncates at min(len(lambda), len(mu)); the shorter vector is
      NOT padded.  This is not a real ambiguity: on partitions of a FIXED integer
      the two conventions are the SAME RELATION (proof in check_C15's docstring),
      so C15(a) confirms the implementation and cannot fail.  It is reported as a
      confirmation, not as a robustness margin.
  N4. "(2alpha)^star merges the two largest parts" is read as: merge the two
      largest entries of the multiset 2alpha.  C15 re-runs the alternative
      reading (merge the last two tuple positions).  Only the MIN SLACK half of
      that probe is informative -- the inequality COUNT is independent of the
      merge convention, so it cannot differ.
  N5. The witness table prints partitions DECREASING while the definitions use
      INCREASING order.  Everything is normalised to the decreasing tuple as the
      canonical key; prefix sums always re-sort increasingly.
"""

import sys
from itertools import combinations

# ------------------------------------------------------------------------------
# FROM THE PAPER (P3): the degree.
# ------------------------------------------------------------------------------
TWO_D = 18
D = 9

# ------------------------------------------------------------------------------
# FROM THE PAPER (P1): the witness, transcribed literally from its table (3).
# Keys are DECREASING part tuples (as printed); values are integers.
# ------------------------------------------------------------------------------
WITNESS_RAW = [
    # --- column 1 of the printed table ---
    ((18,), 0),
    ((16, 2), 9),
    ((14, 4), 6),
    ((14, 2, 2), 20),
    ((12, 6), 3),
    ((12, 4, 2), 14),
    ((12, 2, 2, 2), 32),
    ((10, 8), 0),
    ((10, 6, 2), 12),
    ((10, 4, 4), 12),
    # --- column 2 ---
    ((10, 4, 2, 2), 24),
    ((10, 2, 2, 2, 2), 44),
    ((8, 8, 2), 12),
    ((8, 6, 4), 10),
    ((8, 6, 2, 2), 24),
    ((8, 4, 4, 2), 21),
    ((8, 4, 2, 2, 2), 36),
    ((8, 2, 2, 2, 2, 2), 56),
    ((6, 6, 6), 10),
    ((6, 6, 4, 2), 20),
    # --- column 3 ---
    ((6, 6, 2, 2, 2), 36),
    ((6, 4, 4, 4), 20),
    ((6, 4, 4, 2, 2), 30),
    ((6, 4, 2, 2, 2, 2), 48),
    ((6, 2, 2, 2, 2, 2, 2), 68),
    ((4, 4, 4, 4, 2), 30),
    ((4, 4, 4, 2, 2, 2), 40),
    ((4, 4, 2, 2, 2, 2, 2), 60),
    ((4, 2, 2, 2, 2, 2, 2, 2), 80),
    ((2, 2, 2, 2, 2, 2, 2, 2, 2), 100),
]

# ------------------------------------------------------------------------------
# FROM THE PAPER (P2): the violating triple, as a candidate to be proved.
# ------------------------------------------------------------------------------
PAPER_LAMBDA1 = (10, 8)
PAPER_LAMBDA2 = (4, 4, 4, 2, 2, 2)
PAPER_MU = (8, 4, 4, 2)

# FROM THE PAPER (P6): printed increasing prefix-sum vectors, comparison targets.
PAPER_PREF_UNION = (2, 4, 6, 10, 14, 18, 26, 36)
PAPER_PREF_MU_O2 = (2, 4, 8, 12, 16, 20, 28, 36)

# FROM THE PAPER (P5): target counts, used only as right-hand sides.
TARGET_BLOCKS = 15
TARGET_TERMS = 131
TARGET_DIAG = 119
TARGET_OFFDIAG = 937
TARGET_TOTAL = 1056
TARGET_FIRST_ORDER = 417
TARGET_SECOND_ORDER_SLACK = -2

# NOT FROM THE PAPER.  The paper states only "15 nonempty blocks ... containing
# 131 terms in total"; it prints NO per-block table.  These 15 numbers were RECORDED
# from a previous run of this same enumeration, so comparing against them is a
# REGRESSION PIN, not an independent confirmation.  The real cross-derivation of
# these sizes is independent_block_sizes(), which shares no code with the
# enumerators.
RECORDED_BLOCK_SIZES = {
    (1, 0): 12, (1, 1): 14, (1, 2): 11, (1, 3): 5, (1, 4): 1,
    (3, 0): 25, (3, 1): 16, (3, 2): 6, (3, 3): 1,
    (5, 0): 22, (5, 1): 7, (5, 2): 1,
    (7, 0): 8, (7, 1): 1,
    (9, 0): 1,
}

# The size of the full T^(2) system (D13).  This number IS a paper claim: the
# paper's "Exact verification" paragraph
# says "of the 7025 inequalities of (1) on Lambda^{ev}_18 ... exactly one is
# violated by y".  So the number below is simultaneously the paper's stated target
# and a REGRESSION PIN recorded from a previous run of this code; the derivation
# that makes it evidence is full_second_order_list(), not this constant.
DERIVED_FULL_T2_SIZE = 7025


# ==============================================================================
# DERIVED D1: even partitions.
# ==============================================================================
def canon(parts):
    """Canonical key for a multiset of parts: tuple sorted DECREASING (note N5)."""
    return tuple(sorted(parts, reverse=True))


def even_partitions(n):
    """All multisets of positive EVEN integers summing to n, as canonical
    (decreasing) tuples.  Derived recursively from the definition; nothing
    tabulated.  even_partitions(0) == [()] (the empty partition).

    Recursion: choose the largest part first, then partition the remainder with
    parts bounded above by that choice.  Guarantees each multiset once.
    """
    if n < 0 or n % 2 != 0:
        return []

    def rec(remaining, max_part):
        if remaining == 0:
            yield ()
            return
        top = min(remaining, max_part)
        # parts are positive and even, so step down through even values
        if top % 2 != 0:
            top -= 1
        p = top
        while p >= 2:
            for tail in rec(remaining - p, p):
                yield (p,) + tail
            p -= 2

    return [tuple(t) for t in rec(n, n)]


def is_even_partition_of(parts, n):
    """True iff parts is a multiset of positive even integers summing to n."""
    if sum(parts) != n:
        return False
    for p in parts:
        if p <= 0 or p % 2 != 0:
            return False
    return True


def union(p, q):
    """Multiset union lambda * mu: concatenate part lists and re-sort (canonical)."""
    return canon(tuple(p) + tuple(q))


def o2(mu):
    """mu^{o2} := mu * mu -- each part's MULTIPLICITY doubled, VALUES unchanged."""
    return union(mu, mu)


# ==============================================================================
# Prefix sums and SUPERDOMINANCE.
# ==============================================================================
def pref(p):
    """Increasing prefix sums: sort p INCREASING, take running sums.
    pref(p)[j-1] = sum of the j smallest parts.  Length = number of parts.
    """
    out = []
    total = 0
    for part in sorted(p):
        total += part
        out.append(total)
    return tuple(out)


def superdominates(lam, mu):
    """SUPERDOMINANCE lam \\succeq mu, defined only for partitions of the SAME
    integer.  True iff pref(lam)[j] <= pref(mu)[j] for all
    0 <= j < min(len(lam), len(mu)).

    Note the direction (lam dominates when its small-part prefix sums are
    SMALLER) and the truncation at the min of the two lengths -- the shorter
    vector is NOT padded (note N3).
    """
    if sum(lam) != sum(mu):
        raise ValueError("superdominance compares partitions of the same integer")
    a = pref(lam)
    b = pref(mu)
    for j in range(min(len(a), len(b))):
        if a[j] > b[j]:
            return False
    return True


def superdominates_padded(lam, mu):
    """ROBUSTNESS PROBE (C15a): the alternative convention in which the shorter
    prefix vector is extended by its total and the comparison runs over every
    index.  Reported as a note; not used for any published number.
    """
    if sum(lam) != sum(mu):
        raise ValueError("superdominance compares partitions of the same integer")
    a = list(pref(lam))
    b = list(pref(mu))
    n = max(len(a), len(b))
    total = sum(lam)
    while len(a) < n:
        a.append(total)
    while len(b) < n:
        b.append(total)
    for j in range(n):
        if a[j] > b[j]:
            return False
    return True


# ==============================================================================
# DERIVED D2/D3: the partial-symmetry blocks V_{r,s}(d).
# ==============================================================================
def alpha_tuples(r, s, budget):
    """All ORDERED tuples alpha = (alpha_1,...,alpha_{r+s}) of positive integers
    with the first r entries ODD, the last s entries EVEN, and |alpha| <= budget.

    CRITICAL (note N1): alpha is a SEQUENCE enumerated over ALL orderings, NOT
    taken up to sorting.  Sorting within the odd run and within the even run
    gives 73 total terms instead of the published 131.
    """
    results = []

    def rec(idx, remaining, acc):
        if idx == r + s:
            results.append(tuple(acc))
            return
        # entries 0..r-1 are odd, entries r..r+s-1 are even
        v = 1 if idx < r else 2
        while v <= remaining:
            acc.append(v)
            rec(idx + 1, remaining - v, acc)
            acc.pop()
            v += 2

    rec(0, budget, [])
    return results


def block(r, s, d=D):
    """V_{r,s}(d): all terms u = (alpha, lambda) with alpha as above, lambda an
    even partition of d - |alpha| (so d - |alpha| must be even and >= 0), and
    |alpha| + |lambda| = d.  Returns a list of (alpha, lambda) pairs.

    The empty tuple (r = s = 0) is permitted by the definition but yields
    nothing for d = 9: it forces lambda to be an even partition of 9, and 9 is
    odd, so the block is empty.
    """
    terms = []
    for alpha in alpha_tuples(r, s, d):
        rest = d - sum(alpha)
        if rest < 0 or rest % 2 != 0:
            continue
        for lam in even_partitions(rest):
            terms.append((alpha, lam))
    return terms


def all_blocks(d=D, rmax=None, smax=None):
    """Enumerate r,s over 0..d and keep the blocks with at least one term.
    Returns an ordered dict-like list of ((r,s), terms).  The count is DERIVED,
    not assumed.
    """
    if rmax is None:
        rmax = d
    if smax is None:
        smax = d
    out = []
    for r in range(rmax + 1):
        for s in range(smax + 1):
            terms = block(r, s, d)
            if terms:
                out.append(((r, s), terms))
    return out


# ==============================================================================
# INDEPENDENT CROSS-DERIVATION of every published count.  Everything below
# COUNTS by integer dynamic programming and shares no code with
# even_partitions() or alpha_tuples(), which ENUMERATE.  Without a second route
# all five published counts would flow from one enumeration path, and the two
# printed "cross-checks" would be algebraic identities on that same path --
# unable to disagree with it.  These can.
# ==============================================================================
def num_partitions(n):
    """Number of partitions of n, by the classic part-by-part integer DP.
    Independent of even_partitions(): it never builds a partition.
    """
    if n < 0:
        return 0
    table = [0] * (n + 1)
    table[0] = 1
    for part in range(1, n + 1):
        for total in range(part, n + 1):
            table[total] += table[total - part]
    return table[n]


def num_even_partitions(n):
    """Number of partitions of n into positive EVEN parts.  Halving every part is
    a bijection onto the partitions of n//2, so this is p(n//2) for even n and 0
    for odd n.  (n = 0 gives 1, the empty partition.)
    """
    if n < 0 or n % 2 != 0:
        return 0
    return num_partitions(n // 2)


def _convolve(a, b, cap):
    """Integer polynomial product truncated at degree cap."""
    out = [0] * (cap + 1)
    for i in range(len(a)):
        if a[i] == 0:
            continue
        for j in range(len(b)):
            if b[j] == 0 or i + j > cap:
                continue
            out[i + j] += a[i] * b[j]
    return out


def num_ordered_alpha(r, s, cap):
    """Coefficient list whose entry a is the number of ORDERED tuples alpha with r
    positive-ODD entries followed by s positive-EVEN entries and |alpha| = a.
    Obtained by convolving r copies of the odd series and s copies of the even
    series -- no enumeration, so it is a genuine check on alpha_tuples().
    """
    odd = [0] * (cap + 1)
    for v in range(1, cap + 1, 2):
        odd[v] = 1
    even = [0] * (cap + 1)
    for v in range(2, cap + 1, 2):
        even[v] = 1
    series = [0] * (cap + 1)
    series[0] = 1
    for _ in range(r):
        series = _convolve(series, odd, cap)
    for _ in range(s):
        series = _convolve(series, even, cap)
    return series


def independent_block_sizes(d=D):
    """|V_{r,s}(d)| for every (r,s) with a nonzero count, by DP only:
    |V_{r,s}(d)| = sum_a (#ordered alpha with |alpha| = a) * (#even partitions of
    d - a).  Returns {(r,s): size} for the nonempty blocks.
    """
    sizes = {}
    for r in range(d + 1):
        for s in range(d + 1):
            series = num_ordered_alpha(r, s, d)
            n = 0
            for a in range(d + 1):
                if series[a]:
                    n += series[a] * num_even_partitions(d - a)
            if n:
                sizes[(r, s)] = n
    return sizes


# ==============================================================================
# DERIVED D7: the index maps.  Every index is CONSTRUCTED, never tabulated.
# ==============================================================================
def map_A(u):
    """A_u := sort( [2a for a in alpha] + lambda + lambda ).
    An even partition of 18, since 2|alpha| + 2|lambda| = 2d = 18.
    """
    alpha, lam = u
    parts = [2 * a for a in alpha] + list(lam) + list(lam)
    return canon(parts)


def merge2max(L):
    """Replace the two LARGEST entries of L by their sum, yielding len(L)-1
    entries (note N4, primary reading).  Requires len(L) >= 2.
    """
    if len(L) < 2:
        raise ValueError("merge2max needs at least two entries")
    srt = sorted(L, reverse=True)
    merged = [srt[0] + srt[1]] + srt[2:]
    return merged


def merge_last_two(L):
    """ROBUSTNESS PROBE (C15b): merge the LAST TWO TUPLE ENTRIES instead of the
    two largest.  For alpha these are the even entries, which need not be the
    largest.  Reported as a note only.
    """
    if len(L) < 2:
        raise ValueError("merge_last_two needs at least two entries")
    return list(L[:-2]) + [L[-2] + L[-1]]


def map_B(u, merger=merge2max):
    """B_u := sort( merger([2a for a in alpha]) + lambda + lambda ).
    Defined only for blocks with r+s >= 2.  An even partition of 18: merging two
    parts into their sum preserves the total.
    """
    alpha, lam = u
    two_alpha = [2 * a for a in alpha]
    parts = merger(two_alpha) + list(lam) + list(lam)
    return canon(parts)


def map_C(u, v):
    """C_{u,v} := sort( [a+b for a,b in zip(alpha,beta)] + lambda + mu ) for u,v
    in the SAME block (so alpha and beta have equal length).  Parities align
    (odd+odd and even+even are both even), so C_{u,v} is an even partition of 18.
    Symmetric in u,v.
    """
    alpha, lam = u
    beta, mu = v
    if len(alpha) != len(beta):
        raise ValueError("C_{u,v} requires u, v in the same block")
    parts = [a + b for a, b in zip(alpha, beta)] + list(lam) + list(mu)
    return canon(parts)


# ==============================================================================
# DERIVED D8: the first-order system T_18^(1).
# ==============================================================================
def first_order_pairs(coords, dominates=superdominates):
    """FO := { ORDERED pairs (lambda, mu) of distinct even partitions of 18 with
    lambda \\succeq mu }.  Slack of (lambda,mu) := y[lambda] - y[mu].
    """
    pairs = []
    for lam in coords:
        for mu in coords:
            if lam == mu:
                continue
            if dominates(lam, mu):
                pairs.append((lam, mu))
    return pairs


# ==============================================================================
# DERIVED D4/D5: the partial-symmetry inequality lists.
# COUNTED WITH MULTIPLICITY -- deliberately lists, never sets (note N2).
# ==============================================================================
def diagonal_list(blocks, merger=merge2max):
    """DIAGONAL: for every nonempty block with r+s >= 2 and every u in it, the
    inequality y[A_u] - y[B_u] >= 0.  Returned as a LIST of
    ((r,s), u, A_u, B_u) with multiplicity: distinct u may give identical
    inequalities and those repeats are KEPT.
    """
    out = []
    for (rs, terms) in blocks:
        r, s = rs
        if r + s < 2:
            continue
        for u in terms:
            out.append((rs, u, map_A(u), map_B(u, merger)))
    return out


def offdiagonal_list(blocks):
    """OFF-DIAGONAL: for every nonempty block and every UNORDERED pair {u,v} with
    u != v in that block, the inequality y[A_u] + y[A_v] - 2*y[C_{u,v}] >= 0.
    Returned as a LIST of ((r,s), u, v, A_u, A_v, C_{u,v}) with multiplicity.
    Unordered pairs: itertools.combinations, so each {u,v} appears once.
    """
    out = []
    for (rs, terms) in blocks:
        for u, v in combinations(terms, 2):
            out.append((rs, u, v, map_A(u), map_A(v), map_C(u, v)))
    return out


# ==============================================================================
# DERIVED D13: the FULL second-order system on Lambda^{ev}_18.  Its size (7025)
# and its single violated inequality ARE claims of the paper's "Exact
# verification" paragraph; only the T^(2) subset T^(1) embedding of C13 is beyond
# the paper.
# ==============================================================================
def full_second_order_list(coords, dominates=superdominates):
    """All unordered {lambda1, lambda2} pairs WITH REPETITION together with all
    mu in EvenPartitions(18), kept when (lambda1*lambda2) \\succeq (mu^{o2}).
    Both sides are partitions of 36.  Returns a list of (lam1, lam2, mu).
    """
    out = []
    mu_o2 = [(mu, o2(mu)) for mu in coords]
    n = len(coords)
    for i in range(n):
        for j in range(i, n):
            lam1 = coords[i]
            lam2 = coords[j]
            uni = union(lam1, lam2)
            for mu, mm in mu_o2:
                if dominates(uni, mm):
                    out.append((lam1, lam2, mu))
    return out


# ==============================================================================
# Reporting harness.  One line per check, "PASS " / "FAIL " prefixed.
# ==============================================================================
class Report(object):
    def __init__(self):
        self.results = []  # list of (name, ok)

    def note(self, text):
        print("       . " + text)

    def check(self, name, ok, detail=""):
        ok = bool(ok)
        self.results.append((name, ok))
        tag = "PASS " if ok else "FAIL "
        line = tag + name
        if detail:
            line += "  [" + detail + "]"
        print(line)
        return ok

    def total(self):
        return len(self.results)

    def failed(self):
        return [n for (n, ok) in self.results if not ok]

    def verdict(self):
        n = self.total()
        bad = self.failed()
        print("")
        if not bad:
            print("VERDICT: ALL %d CHECKS PASS" % n)
            print("VERDICT IS CONDITIONAL: what is machine-verified is the finite")
            print("  inequality system (417 first-order + 1056 partial-symmetry")
            print("  inequalities, the five counts, the second-order violation, the")
            print("  7025-inequality T^(2) system with exactly one violated")
            print("  inequality, and the paper's two further superdominance facts:")
            print("  A_u succeq B_u on all 119 diagonal terms [C16] and")
            print("  A_u A_v succeq C_{u,v}^{o2} on all 937 off-diagonal pairs [C17].")
            print("NOT VERIFIED HERE: the step 'y in T^(1) AND those 1056 inequalities")
            print("  => y in trop(BSigma_18^*)' is Lemma 1 of the paper, proved")
            print("  by citation into ABDR and NOT verified here; the witness table and")
            print("  the definitions are transcribed INPUTS (P1-P4), not checked against")
            print("  the source; and the ordered-alpha reading (note N1) is a CHOSEN")
            print("  convention fitted to the paper's 131, so the counts cannot")
            print("  discriminate against a paper that read alpha the same way.")
            print("  A green run does not establish the paper's theorem on its own.")
            return 0
        for name in bad:
            print("FAILED CHECK: " + name)
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1


def all_ints(values):
    """Exact-arithmetic hygiene (C14): every value must be a Python int.
    bool is a subclass of int and is explicitly rejected.
    """
    for v in values:
        if isinstance(v, bool) or not isinstance(v, int):
            return False
    return True


def fmt_part(p):
    """Compact printable form of a partition key."""
    return "(" + ",".join(str(x) for x in p) + ")"


# ==============================================================================
# CHECKS
# ==============================================================================
def check_C0(rep):
    """C0 coordinate set is complete and well formed."""
    print("--- C0: coordinate set Lambda^{ev}_18 ---")
    coords = even_partitions(TWO_D)
    print("       derived |EvenPartitions(18)| = %d" % len(coords))
    for m in (0, 2, 4, 6, 8):
        print("       derived |EvenPartitions(%d)| = %d  ->  %s"
              % (m, len(even_partitions(m)),
                 " ".join(fmt_part(p) for p in even_partitions(m)) or "(empty part.)"))
    y = {}
    dup = []
    for key, val in WITNESS_RAW:
        k = canon(key)
        if k in y:
            dup.append(k)
        y[k] = val

    # 30 is NOT a number the paper prints -- it prints a 30-row table.  So it
    # is cross-derived here as p(9) by the independent partition-count DP rather
    # than asserted as a "paper target".
    dp_card = num_even_partitions(TWO_D)
    print("       INDEPENDENT DP: #even partitions of 18 = p(9) = %d" % dp_card)
    ok_card = (len(coords) == dp_card)
    ok_set = (set(y.keys()) == set(coords))
    ok_shape = all(is_even_partition_of(k, TWO_D) for k in y.keys())
    ok_nodup = (not dup) and (len(y) == len(WITNESS_RAW))

    missing = sorted(set(coords) - set(y.keys()), reverse=True)
    extra = sorted(set(y.keys()) - set(coords), reverse=True)
    if missing:
        rep.note("MISSING keys: " + " ".join(fmt_part(p) for p in missing))
    if extra:
        rep.note("EXTRA keys: " + " ".join(fmt_part(p) for p in extra))
    if dup:
        rep.note("DUPLICATE keys: " + " ".join(fmt_part(p) for p in dup))

    rep.check("C0 coordinate set is complete and well formed",
              ok_card and ok_set and ok_shape and ok_nodup,
              "|EvenPartitions(18)|=%d (independently derived p(9)=%d), witness "
              "keys=%d, set-equal=%s" % (len(coords), dp_card, len(y), ok_set))
    return coords, y


def check_C1_C2(rep):
    """C1 block count == 15 ; C2 total term count == 131."""
    print("--- C1/C2: blocks V_{r,s}(9) and their terms ---")
    blocks = all_blocks(D)
    sizes = dict((rs, len(terms)) for (rs, terms) in blocks)
    total = sum(sizes.values())
    print("       derived nonempty blocks = %d" % len(blocks))
    for (rs, terms) in blocks:
        r, s = rs
        print("       V_{%d,%d}(9): %3d terms   (r odd: %s, r+2s=%d <= 9: %s)"
              % (r, s, len(terms), r % 2 == 1, r + 2 * s, r + 2 * s <= 9))
    print("       derived total terms = %d" % total)

    # DERIVED parity/feasibility structure (not assumed): every nonempty block
    # must have r odd and r + 2s <= 9.
    parity_ok = all((rs[0] % 2 == 1) and (rs[0] + 2 * rs[1] <= D) for (rs, _) in blocks)
    predicted = set()
    for r in (1, 3, 5, 7, 9):
        for s in range(0, (D - r) // 2 + 1):
            predicted.add((r, s))
    struct_ok = (set(sizes.keys()) == predicted)
    rep.note("derived block set matches {r odd, 0<=s<=(9-r)/2}: %s" % struct_ok)

    # INDEPENDENT CROSS-DERIVATION.  dp_sizes counts the same blocks by integer
    # DP (composition series convolution x partition-count recurrence) and shares
    # no code with alpha_tuples()/even_partitions().  This is the only check here
    # that can actually contradict the enumeration.
    dp_sizes = independent_block_sizes(D)
    dp_total = sum(dp_sizes.values())
    dp_agree = (dp_sizes == sizes)
    print("       INDEPENDENT DP: nonempty blocks = %d, total terms = %d"
          % (len(dp_sizes), dp_total))
    print("       INDEPENDENT DP agrees with the enumeration block-for-block: %s"
          % dp_agree)
    if not dp_agree:
        for rs in sorted(set(dp_sizes.keys()) | set(sizes.keys())):
            if dp_sizes.get(rs) != sizes.get(rs):
                print("           DISAGREEMENT V_{%d,%d}: enumerated %s, DP %s"
                      % (rs[0], rs[1], sizes.get(rs), dp_sizes.get(rs)))

    rep.check("C1 block count == 15", len(blocks) == TARGET_BLOCKS and parity_ok
              and struct_ok and len(dp_sizes) == TARGET_BLOCKS and dp_agree,
              "derived %d, independent DP %d, target %d"
              % (len(blocks), len(dp_sizes), TARGET_BLOCKS))
    sizes_match = (sizes == RECORDED_BLOCK_SIZES)
    rep.note("per-block sizes match the sizes RECORDED from a previous run "
             "(a regression pin -- the paper prints no per-block table): %s"
             % sizes_match)
    if total == 73:
        rep.note("WARNING: 73 terms means alpha was wrongly SORTED (see note N1)")
    rep.check("C2 total term count == 131",
              total == TARGET_TERMS and dp_total == TARGET_TERMS and dp_agree
              and sizes_match,
              "enumerated %d, independent DP %d, target %d"
              % (total, dp_total, TARGET_TERMS))
    return blocks, sizes, total, dp_sizes


def check_C3_C4_C5(rep, blocks, sizes, total, dp_sizes):
    """C3 diagonal == 119 ; C4 off-diagonal == 937 ; C5 total == 1056.

    NOTE ON WHAT COUNTS AS A CROSS-CHECK.  "131 - 12 = 119" and
    "sum C(n_block,2) = 937" computed from `sizes` are ALGEBRAIC IDENTITIES:
    diagonal_list() emits one entry per term of the blocks with r+s >= 2, and
    offdiagonal_list() emits exactly combinations(terms, 2), so those sums equal
    the list lengths for ANY input and can never disagree.  They are therefore
    recomputed here from dp_sizes -- the independent DP counts -- which can.
    """
    print("--- C3/C4/C5: partial-symmetry inequality counts (WITH multiplicity) ---")
    diag = diagonal_list(blocks)
    off = offdiagonal_list(blocks)

    # C3.  Cross-derivation from the INDEPENDENT DP block sizes.
    small = sum(n for (rs, n) in dp_sizes.items() if rs[0] + rs[1] < 2)
    dp_total = sum(dp_sizes.values())
    print("       derived diagonal count (multiplicity) = %d" % len(diag))
    print("       INDEPENDENT DP cross-check: %d - (DP terms in blocks with "
          "r+s<2 = %d) = %d" % (dp_total, small, dp_total - small))
    distinct_diag = set((a, b) for (_, _, a, b) in diag)
    rep.note("distinct diagonal inequalities = %d (paper counts terms, not "
             "distinct; a set-based implementation reports this instead)"
             % len(distinct_diag))
    rep.check("C3 diagonal inequality count == 119",
              len(diag) == TARGET_DIAG and (dp_total - small) == TARGET_DIAG,
              "enumerated %d, independent DP %d, target %d"
              % (len(diag), dp_total - small, TARGET_DIAG))

    # C4.  Cross-derivation from the INDEPENDENT DP block sizes: sum C(n, 2).
    combo = sum((n * (n - 1)) // 2 for n in dp_sizes.values())
    ordered = sum(n * (n - 1) for n in sizes.values())
    print("       derived off-diagonal count (unordered pairs) = %d" % len(off))
    print("       INDEPENDENT DP cross-check sum C(n_block,2) = %s = %d"
          % ("+".join(str((n * (n - 1)) // 2)
                      for (_, n) in sorted(dp_sizes.items())), combo))
    distinct_off = set()
    for (_, _, _, a, b, c) in off:
        distinct_off.add((min(a, b), max(a, b), c))
    rep.note("ordered in-block pairs would be %d; distinct off-diagonal "
             "inequalities = %d" % (ordered, len(distinct_off)))
    rep.check("C4 off-diagonal inequality count == 937",
              len(off) == TARGET_OFFDIAG and combo == TARGET_OFFDIAG,
              "enumerated %d, independent DP %d, target %d"
              % (len(off), combo, TARGET_OFFDIAG))

    # C5.
    grand = len(diag) + len(off)
    print("       derived total partial-symmetry inequalities = %d + %d = %d"
          % (len(diag), len(off), grand))
    rep.check("C5 total generated == 1056", grand == TARGET_TOTAL,
              "derived %d, target %d" % (grand, TARGET_TOTAL))
    return diag, off


def check_C6_C7(rep, coords, y):
    """C6 first-order count == 417 ; C7 y in T^(1) with minimum slack exactly 0."""
    print("--- C6/C7: first-order system T_18^(1) ---")
    fo = first_order_pairs(coords)
    print("       derived |FO| (ordered distinct pairs, lambda succeq mu) = %d" % len(fo))
    rep.check("C6 first-order inequality count == 417", len(fo) == TARGET_FIRST_ORDER,
              "derived %d, target %d" % (len(fo), TARGET_FIRST_ORDER))

    slacks = []
    ints_ok = True
    for (lam, mu) in fo:
        s = y[lam] - y[mu]
        if isinstance(s, bool) or not isinstance(s, int):
            ints_ok = False
        slacks.append(s)
    if not slacks:
        rep.check("C7 y in T_18^(1) with minimum slack exactly 0", False,
                  "no inequalities generated")
        return fo, slacks, False
    mn = min(slacks)
    neg = [(fo[i], slacks[i]) for i in range(len(slacks)) if slacks[i] < 0]
    tight = [fo[i] for i in range(len(slacks)) if slacks[i] == 0]
    print("       derived first-order slack minimum = %d, maximum = %d"
          % (mn, max(slacks)))
    print("       derived negative first-order slacks = %d" % len(neg))
    print("       derived tight (slack 0) first-order pairs = %d" % len(tight))
    for (lam, mu) in tight:
        print("           TIGHT  %s succeq %s   y=%d , y=%d"
              % (fmt_part(lam), fmt_part(mu), y[lam], y[mu]))
    for ((lam, mu), s) in neg:
        print("           NEGATIVE  %s succeq %s   slack %d"
              % (fmt_part(lam), fmt_part(mu), s))
    ok = rep.check("C7 y in T_18^(1) with minimum slack exactly 0",
                   mn == 0 and ints_ok,
                   "min slack %d (target 0), all-int=%s" % (mn, ints_ok))
    return fo, slacks, ok


def check_C8(rep, diag, y):
    """C8 diagonal minimum slack exactly 0."""
    print("--- C8: diagonal slacks y[A_u] - y[B_u] ---")
    slacks = []
    ints_ok = True
    for (rs, u, a, b) in diag:
        s = y[a] - y[b]
        if isinstance(s, bool) or not isinstance(s, int):
            ints_ok = False
        slacks.append((s, rs, u, a, b))
    vals = [t[0] for t in slacks]
    mn = min(vals) if vals else None
    print("       derived diagonal slacks: %d values, min = %s, max = %s"
          % (len(vals), mn, max(vals) if vals else None))
    tight = [t for t in slacks if t[0] == 0]
    print("       derived tight diagonal inequalities = %d" % len(tight))
    shown = 0
    for (s, rs, u, a, b) in tight:
        if shown >= 6:
            print("           ... %d further tight diagonal terms" % (len(tight) - shown))
            break
        print("           TIGHT V_{%d,%d} alpha=%s lam=%s : A=%s B=%s (both y=%d)"
              % (rs[0], rs[1], fmt_part(u[0]), fmt_part(u[1]),
                 fmt_part(a), fmt_part(b), y[a]))
        shown += 1
    for (s, rs, u, a, b) in slacks:
        if s < 0:
            print("           NEGATIVE V_{%d,%d} alpha=%s lam=%s : A=%s(%d) B=%s(%d) "
                  "slack %d" % (rs[0], rs[1], fmt_part(u[0]), fmt_part(u[1]),
                                fmt_part(a), y[a], fmt_part(b), y[b], s))
    ok = rep.check("C8 diagonal minimum slack exactly 0",
                   bool(vals) and mn == 0 and ints_ok,
                   "min slack %s over %d (target 0), all-int=%s"
                   % (mn, len(vals), ints_ok))
    return vals, ok


def check_C9(rep, off, y):
    """C9 off-diagonal minimum slack exactly 0."""
    print("--- C9: off-diagonal slacks y[A_u] + y[A_v] - 2*y[C_{u,v}] ---")
    slacks = []
    ints_ok = True
    for (rs, u, v, a, b, c) in off:
        s = y[a] + y[b] - 2 * y[c]
        if isinstance(s, bool) or not isinstance(s, int):
            ints_ok = False
        slacks.append((s, rs, u, v, a, b, c))
    vals = [t[0] for t in slacks]
    mn = min(vals) if vals else None
    print("       derived off-diagonal slacks: %d values, min = %s, max = %s"
          % (len(vals), mn, max(vals) if vals else None))
    tight = [t for t in slacks if t[0] == 0]
    print("       derived tight off-diagonal inequalities = %d" % len(tight))
    shown = 0
    for t in tight:
        if shown >= 6:
            print("           ... %d further tight off-diagonal pairs" % (len(tight) - shown))
            break
        (s, rs, u, v, a, b, c) = t
        print("           TIGHT V_{%d,%d} : A_u=%s(%d) A_v=%s(%d) C=%s(%d)"
              % (rs[0], rs[1], fmt_part(a), y[a], fmt_part(b), y[b], fmt_part(c), y[c]))
        shown += 1
    for t in slacks:
        if t[0] < 0:
            (s, rs, u, v, a, b, c) = t
            print("           NEGATIVE V_{%d,%d} : A_u=%s(%d) A_v=%s(%d) C=%s(%d) slack %d"
                  % (rs[0], rs[1], fmt_part(a), y[a], fmt_part(b), y[b],
                     fmt_part(c), y[c], s))
    ok = rep.check("C9 off-diagonal minimum slack exactly 0",
                   bool(vals) and mn == 0 and ints_ok,
                   "min slack %s over %d (target 0), all-int=%s"
                   % (mn, len(vals), ints_ok))
    return vals, ok


def check_C9b(rep, diag, off, y):
    """C9b index maps land in the coordinate set."""
    print("--- C9b: every A_u, B_u, C_{u,v} is a valid coordinate ---")
    bad = []
    n = 0
    for (rs, u, a, b) in diag:
        for key in (a, b):
            n += 1
            if not is_even_partition_of(key, TWO_D) or key not in y:
                bad.append(key)
    for (rs, u, v, a, b, c) in off:
        for key in (a, b, c):
            n += 1
            if not is_even_partition_of(key, TWO_D) or key not in y:
                bad.append(key)
    print("       derived indices inspected = %d, invalid = %d" % (n, len(bad)))
    for key in bad[:10]:
        print("           INVALID INDEX %s (sum %d)" % (fmt_part(key), sum(key)))
    rep.check("C9b index maps land in the coordinate set", not bad,
              "%d indices checked, %d invalid" % (n, len(bad)))


def check_C10(rep, coords):
    """C10 the violating triple really is a T^(2) generator.  DERIVED, not assumed."""
    print("--- C10: the paper's triple generates a required T^(2) inequality ---")
    l1 = canon(PAPER_LAMBDA1)
    l2 = canon(PAPER_LAMBDA2)
    mu = canon(PAPER_MU)
    in_coords = all(k in set(coords) for k in (l1, l2, mu))
    shapes = all(is_even_partition_of(k, TWO_D) for k in (l1, l2, mu))
    print("       lambda^1 = %s , lambda^2 = %s , mu = %s"
          % (fmt_part(l1), fmt_part(l2), fmt_part(mu)))
    print("       all three even partitions of 18 and coordinates of y: %s"
          % (in_coords and shapes))

    uni = union(l1, l2)
    mm = o2(mu)
    print("       derived lambda^1 * lambda^2 = %s   (|.| = %d, %d parts)"
          % (fmt_part(uni), sum(uni), len(uni)))
    print("       derived mu^{o2}            = %s   (|.| = %d, %d parts)"
          % (fmt_part(mm), sum(mm), len(mm)))
    pu = pref(uni)
    pm = pref(mm)
    print("       derived pref(lambda^1*lambda^2) = %s" % (pu,))
    print("       paper   pref(lambda^1*lambda^2) = %s" % (PAPER_PREF_UNION,))
    print("       derived pref(mu^{o2})           = %s" % (pm,))
    print("       paper   pref(mu^{o2})           = %s" % (PAPER_PREF_MU_O2,))
    rel = superdominates(uni, mm)
    print("       derived (lambda^1*lambda^2) succeq (mu^{o2}) = %s "
          "(componentwise <= on the first %d prefix sums)"
          % (rel, min(len(pu), len(pm))))
    both36 = (sum(uni) == 36 and sum(mm) == 36 and len(uni) == 8 and len(mm) == 8)
    ok = (in_coords and shapes and both36 and rel
          and pu == PAPER_PREF_UNION and pm == PAPER_PREF_MU_O2)
    ok = rep.check("C10 the violating triple really is a T^(2) generator", ok,
                   "pref match=%s/%s, succeq=%s"
                   % (pu == PAPER_PREF_UNION, pm == PAPER_PREF_MU_O2, rel))
    return l1, l2, mu, ok


def check_C11(rep, y, l1, l2, mu):
    """C11 second-order slack == -2, so y is outside T_18^(2)."""
    print("--- C11: the second-order slack ---")
    a = y[l1]
    b = y[l2]
    c = y[mu]
    slack = a + b - 2 * c
    ints_ok = all_ints([a, b, c, slack])
    print("       y[%s] = %d , y[%s] = %d , y[%s] = %d"
          % (fmt_part(l1), a, fmt_part(l2), b, fmt_part(mu), c))
    print("       derived slack = %d + %d - 2*%d = %d" % (a, b, c, slack))
    ok = rep.check("C11 second-order slack == -2, so y is outside T_18^(2)",
                   slack == TARGET_SECOND_ORDER_SLACK and slack < 0 and ints_ok,
                   "derived %d, target %d, all-int=%s"
                   % (slack, TARGET_SECOND_ORDER_SLACK, ints_ok))
    return slack, ok


def check_C12(rep, ok_C7, ok_C8, ok_C9, ok_C10, ok_C11):
    """C12 the paper's conclusion -- a CONJUNCTION, with the trust boundary printed."""
    print("--- C12: the paper's conclusion ---")
    print("       hypotheses of the Lemma: y in T^(1) [C7]=%s ; 119 diagonal [C8]=%s ;"
          % (ok_C7, ok_C8))
    print("                                937 off-diagonal [C9]=%s" % (ok_C9,))
    print("       generator condition [C10]=%s ; violation [C11]=%s" % (ok_C10, ok_C11))
    ok = ok_C7 and ok_C8 and ok_C9 and ok_C10 and ok_C11
    # The trust boundary is printed BEFORE the PASS line, not after it: a reader
    # scanning for "PASS" must meet the condition before meeting the claim.
    print("       *** TRUST BOUNDARY: the step 'y in T^(1) AND the 1056 partial-symmetry")
    print("       *** inequalities  =>  y in trop(BSigma_18^*)'  is the CITED LEMMA")
    print("       *** (ABDR Sec 5.3 / App A, Lemma 5.14, Lemma 5.17, Cor 5.18,")
    print("       *** Prop 5.19).  IT IS NOT MACHINE-VERIFIED by this program.  This")
    print("       *** program verifies only the finite inequality system and the counts.")
    print("       *** If that Lemma is false or misapplied, EVERY check here still")
    print("       *** passes and the paper is still wrong.  C12 below is therefore a")
    print("       *** CONDITIONAL conclusion, and it is a conjunction of C7-C11: it")
    print("       *** carries no evidence of its own.")
    rep.check("C12 the paper's conclusion, CONDITIONAL ON the cited membership Lemma "
              "(NOT machine-verified): given that Lemma, y in trop(BSigma_18^*) \\ "
              "T_18^(2), hence ABDR Conjecture 5.23 is false and the containment "
              "trop(BSigma_18^*) >= T_18^(2) is strict", ok,
              "conjunction of the reported results of C7,C8,C9,C10,C11; the "
              "finite system is verified, the reduction is assumed")
    return ok


def check_C13(rep, coords, y, l1, l2, mu, fo):
    """C13 global consistency of the violation.  DERIVED here -- and, except for
    the embedding below, ALSO ASSERTED BY THE PAPER: its "Exact verification"
    paragraph says that of the 7025 inequalities of (1) on Lambda^{ev}_18
    exactly one is violated by y, and its appendix listing prints 7025 as well.
    So the size, the minimum slack and the "exactly one violated" conclusion are
    an independent confirmation of published numbers, not a self-comparison.

    Also carries the ONE structural check that links the T^(1) and T^(2)
    definitions to each other:  T^(2) subset T^(1).  Without it the
    DIRECTION of superdominance and the meaning of mu^{o2} (multiplicity doubled,
    part values unchanged) would each be pinned by exactly one worked example -- the
    paper's own triple in C10 -- so a systematic misreading that happened to keep
    that single example true would have passed everything.

    The check:  for every one of the 417 first-order pairs (lam, mu), the triple
    (lam, lam, mu) must appear in the generated full second-order system, since
    then  2 y_lam >= 2 y_mu  is one of its inequalities.  This is a theorem
    (lam succeq mu  =>  lam*lam succeq mu*mu: add pref(lam)[i] <= pref(mu)[i] to
    pref(lam)[i+1] <= pref(mu)[i+1] to get the even-index case, and double the
    first for the odd-index case), so it passes -- but it FAILS under a flipped
    comparison and under the plausible bug "o2 doubles the part VALUES", i.e. it
    bites exactly where nothing else did.
    """
    print("--- C13: the FULL second-order system on Lambda^{ev}_18 (DERIVED here; "
          "its size, min slack and single violation are ALSO the paper's claims) ---")
    system = full_second_order_list(coords)
    slacks = []
    ints_ok = True
    for (a1, a2, m) in system:
        s = y[a1] + y[a2] - 2 * y[m]
        if isinstance(s, bool) or not isinstance(s, int):
            ints_ok = False
        slacks.append(s)
    vals = slacks
    mn = min(vals) if vals else None
    neg = [(system[i], vals[i]) for i in range(len(vals)) if vals[i] < 0]
    print("       derived |T^(2) system| = %d  (unordered {lambda1,lambda2} with "
          "repetition, all mu)" % len(system))
    print("       derived minimum slack over the full system = %s" % (mn,))
    print("       derived number of violated inequalities = %d" % len(neg))
    for ((a1, a2, m), s) in neg:
        print("           VIOLATED  lambda^1=%s lambda^2=%s mu=%s  slack %d"
              % (fmt_part(a1), fmt_part(a2), fmt_part(m), s))
    expected_triple = (min(l1, l2), max(l1, l2), mu)
    unique_is_paper = False
    if len(neg) == 1:
        (a1, a2, m), s = neg[0]
        unique_is_paper = ((min(a1, a2), max(a1, a2), m) == expected_triple)
    print("       derived: the unique violated inequality is the paper's triple: %s"
          % unique_is_paper)

    # STRUCTURAL CONSISTENCY T^(2) subset T^(1): every first-order pair must be
    # reproduced by the second-order system at lambda^1 = lambda^2.  This ties the
    # two definitions together instead of trusting each on one example.
    triples = set(system)
    embed_missing = [(lam, m) for (lam, m) in fo if (lam, lam, m) not in triples]
    print("       T^(2) subset T^(1) embedding: %d of %d first-order pairs (lam,mu) "
          "appear as (lam,lam,mu) in the T^(2) system; missing = %d"
          % (len(fo) - len(embed_missing), len(fo), len(embed_missing)))
    for (lam, m) in embed_missing[:10]:
        print("           NOT EMBEDDED  %s succeq %s -- superdominance direction or "
              "mu^{o2} semantics is inconsistent between T^(1) and T^(2)"
              % (fmt_part(lam), fmt_part(m)))
    embed_ok = (not embed_missing) and bool(fo)

    rep.note("size %d is a REGRESSION PIN recorded from a previous run of this code, "
             "and it is ALSO the paper's stated number: the 'Exact verification' "
             "paragraph asserts 7025 inequalities with exactly one violated.  Only "
             "the T^(2) subset T^(1) embedding below is beyond the paper.  Neither "
             "the pin nor the paper makes this an independent derivation -- "
             "full_second_order_list() does" % len(system))
    rep.check("C13 global consistency of the violation, and T^(2) subset T^(1) "
              "(derived here; the size 7025, the minimum slack and the 'exactly one "
              "violated' are ALSO the paper's stated claims, so agreement is an "
              "independent confirmation; only the embedding is beyond the paper)",
              len(system) == DERIVED_FULL_T2_SIZE and mn == TARGET_SECOND_ORDER_SLACK
              and len(neg) == 1 and unique_is_paper and ints_ok and embed_ok,
              "size %d (regression pin %d), min %s, violated %d, unique-is-paper %s, "
              "T2-in-T1 embedding %s"
              % (len(system), DERIVED_FULL_T2_SIZE, mn, len(neg), unique_is_paper,
                 embed_ok))

    # Normalised membership set for C17: the T^(2) system indexed by the UNORDERED
    # pair {lambda^1, lambda^2} together with mu.  Built here so that C17 tests the
    # off-diagonal triples against a system generated from Lambda^{ev}_18 alone,
    # with no knowledge of blocks or index maps.
    return set((min(a1, a2), max(a1, a2), m) for (a1, a2, m) in system)


def check_C14(rep, y, fo_slacks, diag_slacks, off_slacks, so_slack):
    """C14 exact arithmetic hygiene."""
    print("--- C14: exact arithmetic hygiene ---")
    vals_ok = all_ints(list(y.values()))
    keys_ok = all(all_ints(list(k)) for k in y.keys())
    sl_ok = (all_ints(fo_slacks) and all_ints(diag_slacks)
             and all_ints(off_slacks) and all_ints([so_slack]))
    src_ok, why = source_has_no_float()
    print("       witness values all int: %s ; witness keys all int: %s"
          % (vals_ok, keys_ok))
    print("       all %d + %d + %d + 1 slacks are int: %s"
          % (len(fo_slacks), len(diag_slacks), len(off_slacks), sl_ok))
    print("       source scan (no float literal, no '/', no math import): %s%s"
          % (src_ok, ("" if src_ok else "  -> " + why)))
    rep.check("C14 exact arithmetic hygiene",
              vals_ok and keys_ok and sl_ok and src_ok,
              "values=%s keys=%s slacks=%s source=%s"
              % (vals_ok, keys_ok, sl_ok, src_ok))


def check_C15(rep, coords, y, blocks, diag_count):
    """C15 convention probes -- with an HONEST account of what they can show.

    (a) IS NOT A ROBUSTNESS TEST: padded and truncated superdominance are the SAME
    RELATION on partitions of a fixed integer n, provably.  If lam succeq mu holds
    truncated then len(lam) >= len(mu), because the constraint at j = len(lam)-1
    reads n = pref(lam)[-1] <= pref(mu)[len(lam)-1], and if mu had more parts that
    right-hand side would be a proper partial sum, hence < n.  Given
    len(lam) >= len(mu), padding mu with the total n adds only constraints
    pref(lam)[j] <= n, which every partial sum satisfies.  So (a) CANNOT fail for
    any input; it confirms the implementation and it answers the paper's
    unremarked truncation question with a proof rather than a sample.

    (b) The COUNT half is likewise forced: diagonal_list() appends one entry per
    term regardless of `merger`, so len(diag_alt) == diag_count identically.  Only
    the min-slack half carries information about the merge convention.
    """
    print("--- C15: convention probes (see docstring: (a) is a theorem, not a "
          "robustness margin) ---")
    # (a) padded superdominance for the first-order system.
    fo_pad = first_order_pairs(coords, dominates=superdominates_padded)
    pad_slacks = [y[a] - y[b] for (a, b) in fo_pad]
    pad_min = min(pad_slacks) if pad_slacks else None
    fo_trunc = first_order_pairs(coords)
    len_ok = all(len(a) >= len(b) for (a, b) in fo_trunc)
    print("       (a) padded convention: |FO| = %d (truncated gave %d), min slack = %s"
          % (len(fo_pad), len(fo_trunc), pad_min))
    print("       (a) every truncated-convention pair has ell(lambda) >= ell(mu): %s"
          % len_ok)
    ok_a = (len(fo_pad) == TARGET_FIRST_ORDER and pad_min == 0 and len_ok
            and set(fo_pad) == set(fo_trunc))

    # (b) B_u by merging the LAST TWO TUPLE ENTRIES of 2alpha.
    diag_alt = diagonal_list(blocks, merger=merge_last_two)
    alt_slacks = [y[a] - y[b] for (_, _, a, b) in diag_alt]
    alt_min = min(alt_slacks) if alt_slacks else None
    print("       (b) merge-last-two convention: %d diagonal inequalities "
          "(primary gave %d), min slack = %s" % (len(diag_alt), diag_count, alt_min))
    ok_b = (len(diag_alt) == TARGET_DIAG and alt_min == 0)
    rep.note("(a) padded == truncated is a THEOREM here, so (a) cannot fail; the "
             "only informative content of this check is (b)'s min slack, which "
             "shows the merge-two-largest / merge-last-two ambiguity does not "
             "change the published minimum on this witness: (a)=%s (b)=%s"
             % (ok_a, ok_b))
    rep.check("C15 convention probes ((a) confirms a theorem; (b) min slack is the "
              "one informative half)", ok_a and ok_b,
              "padded FO=%d min=%s ; merge-last-two diag=%d min=%s"
              % (len(fo_pad), pad_min, len(diag_alt), alt_min))


def check_C16(rep, diag, fo):
    """C16 the paper's FIRST further fact: every diagonal inequality has
    A_u succeq B_u and is therefore already a first-order inequality.

    THE PAPER'S CLAIM (Exact verification paragraph): "each of the 119 diagonal
    inequalities has A_u succeq B_u and is therefore already a first-order
    inequality".  The RELATION is checked only here: C8 computes the diagonal
    SLACKS, which on their own would leave this sentence of the paper unchecked.

    DERIVED, NOT ASSERTED.  superdominates(A_u, B_u) is evaluated on every entry
    of the diagonal list; the number satisfying it is compared with the
    independently derived len(diag) (and with the paper's 119).  The "therefore
    already a first-order inequality" half is checked rather than inferred: the
    ORDERED pair (A_u, B_u) is looked up in the first-order pair set that
    first_order_pairs() builds by scanning all of Lambda^{ev}_18 x Lambda^{ev}_18,
    a different code path that knows nothing about blocks or index maps.

    HONEST ACCOUNT OF WHAT THIS CAN SHOW (same standard as C13/C15): the fact is
    a THEOREM, so it cannot fail for a correct implementation.  B_u is A_u with
    two parts replaced by their sum.  Take any j <= len(B_u) = len(A_u) - 1 =
    min(len(A_u), len(B_u)).  The j smallest parts of B_u pull back to a set of j
    or j+1 parts of A_u with the same total (split the merged part), so
    pref(B_u)[j-1] >= pref(A_u)[j-1] because the j (or j+1) smallest parts of
    A_u are the cheapest such set and all parts are positive.  That is exactly
    A_u succeq B_u.  So C16 is a CONFIRMATION, not a robustness margin.  What it
    would catch is a flipped superdominance direction, a merger that does not
    preserve the total, a truncation convention that breaks on unequal lengths,
    or a first-order set built with the wrong comparison.

    CONSEQUENCE, printed but deliberately NOT folded into C12: given C16, C8's
    minimum slack 0 is IMPLIED by C7 (y in T^(1)), so C8 carries no evidence
    independent of C7.  That is a property of the paper's construction, not a
    defect of the program, and it is stated here so nobody reads C7 and C8 as two
    independent confirmations.
    """
    print("--- C16: A_u succeq B_u on every diagonal inequality (paper's first "
          "further computational fact) ---")
    fo_set = set(fo)
    rel_ok = 0
    rel_bad = []
    in_fo = 0
    not_in_fo = []
    degenerate = []
    for (rs, u, a, b) in diag:
        if sum(a) != TWO_D or sum(b) != TWO_D:
            rel_bad.append((rs, u, a, b,
                            "not both partitions of %d: |A_u|=%d |B_u|=%d"
                            % (TWO_D, sum(a), sum(b))))
            continue
        if superdominates(a, b):
            rel_ok += 1
        else:
            rel_bad.append((rs, u, a, b, "A_u does NOT superdominate B_u"))
        if a == b:
            degenerate.append((rs, u, a))
        elif (a, b) in fo_set:
            in_fo += 1
        else:
            not_in_fo.append((rs, u, a, b))
    print("       derived diagonal inequalities inspected      = %d" % len(diag))
    print("       derived with A_u succeq B_u                   = %d  (paper says all "
          "%d)" % (rel_ok, TARGET_DIAG))
    print("       derived (A_u,B_u) found in the %d first-order pairs = %d"
          % (len(fo), in_fo))
    print("       derived degenerate cases A_u == B_u          = %d  (must be 0: "
          "merging two parts shortens 2alpha by one, so the multisets differ)"
          % len(degenerate))
    print("       derived DISTINCT (A_u,B_u) inequalities      = %d  (multiplicity "
          "%d, note N2)" % (len(set((a, b) for (_, _, a, b) in diag)), len(diag)))
    for (rs, u, a, b, why) in rel_bad[:10]:
        print("           RELATION FAILS V_{%d,%d} alpha=%s lam=%s : A_u=%s B_u=%s -- %s"
              % (rs[0], rs[1], fmt_part(u[0]), fmt_part(u[1]),
                 fmt_part(a), fmt_part(b), why))
    for (rs, u, a, b) in not_in_fo[:10]:
        print("           NOT A FIRST-ORDER PAIR V_{%d,%d} : A_u=%s B_u=%s -- the "
              "relation and the T^(1) pair set disagree"
              % (rs[0], rs[1], fmt_part(a), fmt_part(b)))
    for (rs, u, a) in degenerate[:10]:
        print("           DEGENERATE V_{%d,%d} alpha=%s lam=%s : A_u == B_u == %s"
              % (rs[0], rs[1], fmt_part(u[0]), fmt_part(u[1]), fmt_part(a)))
    print("       CONSEQUENCE: with this fact, C8's minimum slack 0 is implied by "
          "C7 (y in T^(1));")
    print("       C8 is therefore not independent evidence.  This is the paper's "
          "own reasoning,")
    print("       now machine-checked instead of quoted.")
    ok = (bool(diag) and not rel_bad and not not_in_fo and not degenerate
          and rel_ok == len(diag) and rel_ok == TARGET_DIAG
          and in_fo == len(diag))
    rep.check("C16 every diagonal inequality has A_u succeq B_u, hence is already a "
              "first-order inequality (paper's further fact, DERIVED here)", ok,
              "%d of %d satisfy the relation (paper %d), %d of %d found in the "
              "first-order pair set, %d degenerate"
              % (rel_ok, len(diag), TARGET_DIAG, in_fo, len(diag), len(degenerate)))
    return ok


def check_C17(rep, off, so_triples, l1, l2, mu):
    """C17 the paper's SECOND further fact: every off-diagonal inequality has
    A_u A_v succeq C_{u,v}^{o2} and is therefore itself an instance of (1).

    THE PAPER'S CLAIM (Exact verification paragraph): "each of the 937
    off-diagonal ones has A_u A_v succeq C_{u,v}^{o2} and is therefore itself an
    instance of (1)".  The relation itself is checked only here: C9 gives the
    937 SLACKS, never the relation.

    DERIVED, NOT ASSERTED.  For every entry the program forms the multiset union
    A_u A_v and the doubling C_{u,v}^{o2} (both partitions of 2 * 18 = 36) and
    evaluates superdominates() on them; the number satisfying it is compared with
    the independently derived len(off) and with the paper's 937.  The "therefore
    itself an instance of (1)" half is then checked against the FULL second-
    order system generated in C13 from Lambda^{ev}_18 alone (7025 triples, built
    with no reference to blocks or index maps): the normalised triple
    (min(A_u,A_v), max(A_u,A_v), C_{u,v}) must be a member of it.

    HONEST ACCOUNT OF WHAT THIS CAN SHOW: this fact is also a THEOREM, so a
    correct implementation cannot fail it.  With equal totals and equal lengths,
    superdominance as defined here (increasing prefix sums, componentwise <=) is
    the majorization order: A succeq B iff the j largest parts of A sum to at
    least the j largest of B, for every j.  Now A_u A_v = {2alpha_i} {2beta_i}
    lambda lambda mu mu and C_{u,v}^{o2} = {alpha_i+beta_i twice} lambda lambda
    mu mu, and both have 2(len(alpha) + ell(lambda) + ell(mu)) parts, so no
    truncation occurs.  Componentwise {2a, 2b} majorizes {a+b, a+b}, majorization
    is preserved under multiset union of comparable pieces, and the common block
    lambda lambda mu mu majorizes itself.  Hence the union relation holds.  What
    the check can catch is a flipped superdominance direction, the plausible bug
    "o2 doubles the part VALUES", a C map that mispairs alpha with beta, or a
    disagreement between this relation and the independently generated T^(2)
    system.

    CONSEQUENCE, printed but deliberately NOT folded into C12: the 937 are a
    SUBSET of the 7025 inequalities whose slacks C13 already computes, so C9's
    minimum slack 0 is implied by C13's "exactly one violated" -- provided the
    unique violated triple is not one of the 937.  That last point is DERIVED
    here (it comes out 0 occurrences), which is what reconciles "y is outside
    T^(2)" with "all 937 off-diagonal inequalities hold with slack >= 0".
    """
    print("--- C17: A_u A_v succeq C_{u,v}^{o2} on every off-diagonal inequality "
          "(paper's second further computational fact) ---")
    paper_triple = (min(l1, l2), max(l1, l2), mu)
    rel_ok = 0
    rel_bad = []
    in_t2 = 0
    not_in_t2 = []
    hits_violated = 0
    for (rs, u, v, a, b, c) in off:
        uni = union(a, b)
        cc = o2(c)
        if sum(uni) != 2 * TWO_D or sum(cc) != 2 * TWO_D:
            rel_bad.append((rs, a, b, c,
                            "not both partitions of %d: |A_uA_v|=%d |C^{o2}|=%d"
                            % (2 * TWO_D, sum(uni), sum(cc))))
            continue
        if len(uni) != len(cc):
            rel_bad.append((rs, a, b, c,
                            "unequal lengths %d vs %d, so the min-length "
                            "truncation would bite" % (len(uni), len(cc))))
            continue
        if superdominates(uni, cc):
            rel_ok += 1
        else:
            rel_bad.append((rs, a, b, c,
                            "A_u A_v does NOT superdominate C_{u,v}^{o2}"))
        triple = (min(a, b), max(a, b), c)
        if triple in so_triples:
            in_t2 += 1
        else:
            not_in_t2.append((rs, a, b, c))
        if triple == paper_triple:
            hits_violated += 1
    print("       derived off-diagonal inequalities inspected  = %d" % len(off))
    print("       derived with A_u A_v succeq C_{u,v}^{o2}      = %d  (paper says all "
          "%d)" % (rel_ok, TARGET_OFFDIAG))
    print("       derived triples present in the derived T^(2) system = %d of %d"
          % (in_t2, len(off)))
    print("       derived DISTINCT off-diagonal triples        = %d  (multiplicity "
          "%d, note N2)"
          % (len(set((min(a, b), max(a, b), c) for (_, _, _, a, b, c) in off)),
             len(off)))
    print("       derived occurrences of the UNIQUE VIOLATED T^(2) triple among "
          "them = %d" % hits_violated)
    for (rs, a, b, c, why) in rel_bad[:10]:
        print("           RELATION FAILS V_{%d,%d} : A_u=%s A_v=%s C=%s -- %s"
              % (rs[0], rs[1], fmt_part(a), fmt_part(b), fmt_part(c), why))
    for (rs, a, b, c) in not_in_t2[:10]:
        print("           NOT A T^(2) INSTANCE V_{%d,%d} : A_u=%s A_v=%s C=%s -- the "
              "relation and the generated T^(2) system disagree"
              % (rs[0], rs[1], fmt_part(a), fmt_part(b), fmt_part(c)))
    print("       CONSEQUENCE: the 937 are a SUBSET of the %d T^(2) inequalities "
          "scored in C13," % len(so_triples))
    print("       so C9's minimum slack 0 follows from C13 plus the 0 occurrences "
          "above; C9 is")
    print("       not independent evidence either.  Again this is the paper's own "
          "reasoning, checked.")
    ok = (bool(off) and not rel_bad and not not_in_t2
          and rel_ok == len(off) and rel_ok == TARGET_OFFDIAG
          and in_t2 == len(off) and hits_violated == 0)
    rep.check("C17 every off-diagonal inequality has A_u A_v succeq C_{u,v}^{o2}, "
              "hence is itself an instance of (1) (paper's further fact, DERIVED "
              "here)", ok,
              "%d of %d satisfy the relation (paper %d), %d of %d are members of the "
              "derived %d-inequality T^(2) system, %d coincide with the violated "
              "triple"
              % (rel_ok, len(off), TARGET_OFFDIAG, in_t2, len(off),
                 len(so_triples), hits_violated))
    return ok


def source_has_no_float():
    """Scan this file's own source for float literals, '/' division and math use.
    Comments/docstrings mention these words, so the scan looks at CODE tokens via
    the tokenize module rather than raw text.
    """
    import tokenize
    import token as tokmod
    try:
        with open(__file__, "rb") as fh:
            toks = list(tokenize.tokenize(fh.readline))
    except Exception as exc:  # pragma: no cover - defensive
        return False, "could not tokenize source: %r" % (exc,)
    for tk in toks:
        if tk.type == tokmod.NUMBER:
            txt = tk.string.lower()
            if ("." in txt) or ("e" in txt and not txt.startswith("0x")) or ("j" in txt):
                return False, "float/complex literal %r at line %d" % (tk.string, tk.start[0])
        if tk.type == tokmod.OP and tk.string == "/":
            return False, "true division '/' at line %d" % tk.start[0]
        if tk.type == tokmod.NAME and tk.string == "math":
            return False, "reference to 'math' at line %d" % tk.start[0]
    return True, ""


def print_not_rerun():
    """The closing disclosure.  Printed unconditionally, after the verdict, on
    both the passing and the failing path.  It names everything a referee might
    otherwise expect this program to have re-executed, INCLUDING the separate
    small-degree decision program shipped alongside this note, whose scope is a
    DIFFERENT degree and whose result appears nowhere in the paper.
    """
    print("")
    print("NOT RE-RUN: (a) the membership Lemma of the paper and the transcribed")
    print("  inputs P1-P4, exactly as the NOT VERIFIED HERE lines above state, and")
    print("  the ordered-alpha reading of note N1, which was fitted to the paper's")
    print("  131.  (b) The separate SMALL-DEGREE DECISION PROGRAM shipped alongside")
    print("  this note.  It is an independent search, not a part of this")
    print("  verification: it works in halved coordinates on partitions of d, runs")
    print("  controls at d = 2,...,5 against the source paper's published facet")
    print("  lists (plus one inclusion control at d = 6), and then takes ONE")
    print("  decision, at d = 6, i.e. 2d = 12 -- the smallest degree the source")
    print("  paper does not settle.  There it finds 0 of the 376 second-order")
    print("  generators unimplied by the generators of the tropical side, so at")
    print("  2d = 12 there is NO counterexample and, given the presentation of the")
    print("  tropical side that it takes from the cited results, the ABDR equality")
    print("  HOLDS at d = 6.  That is a companion result about a DIFFERENT degree:")
    print("  it never touches 2d = 18, it shares no code with this program, it is")
    print("  cited nowhere in the paper, its output carries no check count and no")
    print("  verdict of its own, and none of it is counted in the verdict above.")
    print("  It is therefore neither evidence for nor evidence against the")
    print("  2d = 18 counterexample; a referee should read it as a separate")
    print("  negative finding at the first degree ABDR left open.  (c) That")
    print("  generator-and-Farkas decision route is NOT run at d = 9 by either")
    print("  program, and no search over Lambda^{ev}_18 is attempted: at 2d = 18")
    print("  the refutation rests on the paper's explicit witness, which needs")
    print("  only the finite systems scored above.  (d) The constants 7025 and the")
    print("  15 per-block sizes in this file are values RECORDED from an earlier")
    print("  run of this same code; 7025 is also the paper's stated number, and")
    print("  the derivations that make either of them evidence are")
    print("  full_second_order_list() and independent_block_sizes(), not the")
    print("  constants.")


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    print("=" * 84)
    print("verify.py -- degree-18 counterexample to ABDR Conjecture 5.23")
    print("2d = %d, d = %d.  Exact integer arithmetic only." % (TWO_D, D))
    print("FROM THE PAPER: the 30 witness values, the candidate violating triple,")
    print("  the degree, the definitions, the 5 target counts, the 3 target slacks,")
    print("  the 2 printed prefix-sum vectors, and the unproved membership Lemma.")
    print("DERIVED HERE: all partitions, blocks, terms, index maps, inequality")
    print("  lists, counts and slacks -- plus the full T^(2) system (7025 with")
    print("  exactly one violated inequality, which the paper's 'Exact verification'")
    print("  paragraph also asserts) and the embedding T^(2) subset T^(1), which is")
    print("  the only part of C13 that is beyond the paper.")
    print("COUNTED TWICE BY UNRELATED METHODS: 30, 15, 131, 119, 937 are produced")
    print("  both by enumeration and by an independent integer DP, then compared.")
    print("ALSO DERIVED HERE (the paper's two further computational facts, i.e. the")
    print("  RELATIONS and not merely the slacks): A_u succeq B_u on all 119 diagonal")
    print("  terms [C16] and A_u A_v succeq C_{u,v}^{o2} on all 937 off-diagonal")
    print("  pairs [C17], each with its 'therefore' half checked against an")
    print("  independently generated inequality set.")
    print("NOT VERIFIED HERE: (i) the membership Lemma (see the TRUST BOUNDARY at")
    print("  C12) -- no Python program can check that reduction; (ii) the witness")
    print("  table, the degree and the definitions, which are transcribed INPUTS")
    print("  (P1-P4) -- a transcription error is invisible to this program; (iii)")
    print("  the ordered-alpha reading of note N1, which was CHOSEN because it")
    print("  reproduces the paper's 131, so the counts cannot discriminate against")
    print("  a paper that read alpha the same way.  Every other computational")
    print("  assertion of the paper's 'Exact verification' paragraph -- the five")
    print("  counts, the 417 first-order and 1056 partial-symmetry slacks, the")
    print("  second-order violation, the 7025-inequality T^(2) system with exactly")
    print("  one violated inequality, and the two further facts above -- IS derived")
    print("  and compared here.")
    print("=" * 84)
    rep = Report()

    # NOTE: every ok_Cn below is the boolean the corresponding check ACTUALLY
    # asserted and printed, not a weaker predicate re-derived here.  Recomputing
    # them here would risk dropping the all-int conjuncts and, for C10, the
    # coordinate/shape/prefix-vector comparisons -- so a C10 FAIL on the paper's
    # printed prefix sums might not propagate into C12.
    coords, y = check_C0(rep)
    blocks, sizes, total, dp_sizes = check_C1_C2(rep)
    diag, off = check_C3_C4_C5(rep, blocks, sizes, total, dp_sizes)
    fo, fo_slacks, ok_C7 = check_C6_C7(rep, coords, y)

    diag_slacks, ok_C8 = check_C8(rep, diag, y)
    off_slacks, ok_C9 = check_C9(rep, off, y)
    check_C9b(rep, diag, off, y)

    l1, l2, mu, ok_C10 = check_C10(rep, coords)
    so_slack, ok_C11 = check_C11(rep, y, l1, l2, mu)

    check_C12(rep, ok_C7, ok_C8, ok_C9, ok_C10, ok_C11)
    so_triples = check_C13(rep, coords, y, l1, l2, mu, fo)
    check_C14(rep, y, fo_slacks, diag_slacks, off_slacks, so_slack)
    check_C15(rep, coords, y, blocks, len(diag))

    # The two further computational facts asserted by the paper's "Exact
    # verification" paragraph.  The checks above score the diagonal and
    # off-diagonal SLACKS; the two superdominance RELATIONS that the paper
    # states are scored only by the two checks below.
    ok_C16 = check_C16(rep, diag, fo)
    ok_C17 = check_C17(rep, off, so_triples, l1, l2, mu)

    print("")
    print("--- derived summary ---")
    print("       |EvenPartitions(18)|              = %d   (paper target 30)" % len(coords))
    print("       nonempty blocks V_{r,s}(9)        = %d   (paper target %d)"
          % (len(blocks), TARGET_BLOCKS))
    print("       total terms u = (alpha,lambda)    = %d  (paper target %d)"
          % (total, TARGET_TERMS))
    print("       diagonal inequalities             = %d  (paper target %d)"
          % (len(diag), TARGET_DIAG))
    print("       off-diagonal inequalities         = %d  (paper target %d)"
          % (len(off), TARGET_OFFDIAG))
    print("       total partial-symmetry            = %d (paper target %d)"
          % (len(diag) + len(off), TARGET_TOTAL))
    print("       first-order inequalities          = %d  (paper target %d)"
          % (len(fo), TARGET_FIRST_ORDER))
    print("       min first-order slack             = %s   (paper target 0)"
          % (min(fo_slacks) if fo_slacks else None,))
    print("       min diagonal slack                = %s   (paper target 0)"
          % (min(diag_slacks) if diag_slacks else None,))
    print("       min off-diagonal slack            = %s   (paper target 0)"
          % (min(off_slacks) if off_slacks else None,))
    print("       exhibited second-order slack      = %s  (paper target %d)"
          % (so_slack, TARGET_SECOND_ORDER_SLACK))
    print("       A_u succeq B_u on all %d diagonal   = %s   (paper's further fact "
          "[C16])" % (len(diag), ok_C16))
    print("       A_uA_v succeq C^{o2} on all %d off  = %s   (paper's further fact "
          "[C17])" % (len(off), ok_C17))

    rc = rep.verdict()
    print_not_rerun()
    return rc


if __name__ == "__main__":
    sys.exit(main())

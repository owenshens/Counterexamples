#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- independent verifier for

    "A Counterexample for [65] to an Exhaustiveness Conjecture of
     Espinosa-Garcia and Pellicer"

Claim under test.  The partition P = {A_1,...,A_7} displayed in the proof is a
balanced Sidon-Ramsey partition of [65] into seven Sidon sets; neither P nor its
reflection rho(P), rho(x) = 66 - x, is one of the five partitions displayed on
p.6 of arXiv:2309.08553v1 or a reflection of one.  Hence the (unnumbered)
exhaustiveness conjecture quoted from that page is false.

Python 3.9+, standard library only.  All arithmetic is exact integer arithmetic
on ints and sets; there is no floating point anywhere in this file and no
numeric tolerance is used, so every check below is an exact combinatorial
decision, not an approximation.

-------------------------------------------------------------------------------
TAKEN FROM THE PAPER (data; transcribed verbatim, never used as a check result)
-------------------------------------------------------------------------------
  * W = [A_1,...,A_7], the exhibited partition, from the align* display in the
    proof (source .tex lines 65-71).
  * FIVE = the five balanced Sidon-Ramsey partitions of [65] displayed on p.6 of
    arXiv:2309.08553v1 (e-print tarball, file Further.tex, lines 34-77).  The
    counterexample paper cites these and reprints only their 10-element parts;
    the full seven-part lists are transcribed here from the EGP source so that
    membership in the conjectured family can be decided directly rather than
    through the paper's shortcut.
  * PAPER_SIGMA_ROW  = the second data row of the paper's table, |Sigma(A_i)|.
  * PAPER_SIZE_ROW   = the first data row of the paper's table, |A_i|.
  * PAPER_EIGHT      = the eight distinct 10-element sets the paper prints
    (source .tex lines 100-107), in the paper's row-major reading order.
  * PAPER_DIFF_TOTAL = 270, the paper's count of pairwise differences checked.
  * HAND_RHO_A2      = {1,6,8,16,22,33,42,45,46,64}, the set {66-a : a in A_2}
    written out by hand from the paper's A_2 display.  It is not a further datum
    of the paper: the paper states only that rho(A_2) contains {1,64}.  It is
    kept as an independent cross-check on the code's implementation of rho, and
    a referee can audit it by reflecting the ten entries of A_2 by hand.
  * rho(x) = 66 - x, the reflection defined in the quoted conjecture.
  * QUOTE_FRAGMENTS = the two sentences the paper quotes from p.6 of the
    e-print, and KNOWN_EGP_MD5 = the md5 of the Further.tex they and the five
    partitions were transcribed from.  These are used only by the optional
    section (8), which runs when the path to that file is given on the command
    line; without an argument section (8) reports itself SKIPPED and every
    mathematical check still runs.

-------------------------------------------------------------------------------
DERIVED HERE (everything the checks actually decide)
-------------------------------------------------------------------------------
  * Sigma(A) = {a+b : a,b in A, a <= b} for every part of every object below,
    its cardinality, and the comparison with binomial(|A|+1, 2).
  * The full multiset of representations a+b (a <= b) per part, and its maximum
    multiplicity, as an independent re-derivation of the Sidon property from the
    definition ("every integer has at most one representation").
  * The multiset of positive differences b - a per part, its size, and whether
    the differences are pairwise distinct; and the total difference count.
  * Disjointness, union = [65], part sizes, and the balance condition, for W,
    for rho(W), and for each of the five EGP partitions.
  * The complete list of balanced size profiles of 7 parts summing to 65,
    enumerated by brute force (no arithmetic shortcut).
  * F, the conjectured family: canonical (frozenset of frozensets) forms of the
    five displayed partitions together with their rho-images; |F|; whether F is
    rho-invariant; whether every member of F is a balanced Sidon-Ramsey
    partition of [65].
  * canon(W) in F?  canon(rho(W)) in F?  canon(rho(W)) == canon(W)?
  * The 10-element part slots of the five partitions, how many are distinct, the
    distinct sets themselves (compared against PAPER_EIGHT), which repeat, and
    whether any contains {2,65} or {1,64}.
  * rho(A_2), and the shortcut's conclusion re-derived over all ten members
    of F.
  * Optionally (section 8, when a path is supplied): the md5/sha256 of the
    supplied Further.tex, and whether all 35 transcribed parts and the two
    quoted sentences occur in it verbatim.

NOT decided here, and printed at the end as a numbered "gaps" block rather than
skipped silently: (i) that the quoted sentence means exhaustiveness rather than
merely that each displayed partition is balanced -- a reading, not a computation;
(ii) without a source file on the command line, that the five transcribed
partitions are ALL the balanced ones displayed on p.6 (a sixth would shrink F);
(iii) whether the quoted sentence survives into the version of record, Discrete
Appl. Math. 378 (2026), 120-124, doi:10.1016/j.dam.2025.07.002.  That article is
paywalled; the paper scopes its claim to arXiv v1, which is what is checked.
Gap (ii) is closed by section (8) when a Further.tex path is supplied: it parses
the five displays out of the source and requires equality with the transcription
in both directions.  Two checks below (rho-invariance of F, and rho being a
Sidon-preserving involution of [65]) are identities that no input can falsify;
they are labelled as such in their own output lines so they are not mistaken for
evidence about the paper.

Output: one "PASS "/"FAIL " line per check, indented lines for derived
intermediates, then "VERDICT: ...".  Exit status 0 iff every check passes.
"""

import re
import sys
from collections import Counter
from itertools import combinations
from math import comb

N = 65          # ground set is [N] = {1,...,N}
K = 7           # number of parts
REFL = 66       # rho(x) = REFL - x

# --- TAKEN FROM THE PAPER: the exhibited partition (proof display) -----------
W = [
    {1, 5, 13, 23, 26, 28, 37, 56, 57, 63},          # A_1
    {2, 20, 21, 24, 33, 44, 50, 58, 60, 65},         # A_2
    {3, 4, 27, 30, 32, 39, 43, 49, 64},              # A_3
    {6, 15, 16, 22, 36, 40, 48, 51, 53},             # A_4
    {7, 14, 19, 25, 29, 38, 52, 54, 55},             # A_5
    {8, 11, 17, 34, 41, 45, 46, 59, 61},             # A_6
    {9, 10, 12, 18, 31, 35, 42, 47, 62},             # A_7
]

# --- TAKEN FROM THE PAPER's table (compared against, never assigned) --------
PAPER_SIZE_ROW = [10, 10, 9, 9, 9, 9, 9]
PAPER_SIGMA_ROW = [55, 55, 45, 45, 45, 45, 45]
PAPER_DIFF_TOTAL = 270

# --- TAKEN FROM THE PAPER: eight distinct 10-element parts, row-major -------
PAPER_EIGHT = [
    {1, 3, 15, 22, 30, 33, 46, 50, 55, 56},
    {2, 4, 16, 23, 31, 34, 47, 51, 56, 57},
    {2, 9, 17, 21, 26, 27, 47, 49, 60, 63},
    {2, 6, 14, 24, 27, 29, 38, 57, 58, 64},
    {7, 10, 15, 19, 33, 43, 44, 50, 63, 65},
    {3, 9, 10, 29, 38, 40, 43, 53, 61, 65},
    {3, 5, 17, 24, 32, 35, 48, 52, 57, 58},
    {2, 6, 14, 16, 19, 37, 38, 44, 53, 64},
]

# --- INDEPENDENT HAND COMPUTATION: {66-a : a in A_2}, reflected by hand from
# --- the paper's A_2 display, as a cross-check on rho() below ---------------
HAND_RHO_A2 = {1, 6, 8, 16, 22, 33, 42, 45, 46, 64}

# --- Provenance of the EGP data, for the optional section (8) ---------------
# md5 of the Further.tex used for the transcription below, taken from the
# arXiv:2309.08553v1 e-print tarball.
KNOWN_EGP_MD5 = "b91fe6cef3bd7d0177a2e362a088d14e"
# sha256 of the same file.  md5 is collision-broken, so the stronger digest is
# pinned as well and both are required to match.
KNOWN_EGP_SHA256 = ("cfa3504997f84cfcadaa7170c3f79fa4245"
                    "dbad75c4d910ac1a75f65b8feecfd")
# The two sentences of the EGP source that delimit the region holding the five
# displayed partitions.  Used by parse_five_from_source() so that section (8)
# can check the transcription for COMPLETENESS (the source region contains no
# sixth partition) and not merely for inclusion.  Both must occur exactly once.
SRC_REGION_START = ("So far we have found $5$ balanced Sidon-Ramsey partitions "
                    "of $[65]$ with $7$ parts each")
SRC_REGION_END = "These partitions have two parts of $10$ elements each"
# The two sentences the counterexample paper quotes from p.6 of that e-print.
QUOTE_FRAGMENTS = [
    ("These partitions have two parts of $10$ elements each, and five parts "
     "with $9$ elements each.", "size sentence"),
    ("We conjecture that these and their reflected partitions (constructed by "
     "including the numbers $66-x$ instead of $x$ in each part) are all "
     "balanced Sidon-Ramsey partitions with those parameters.",
     "exhaustiveness conjecture"),
]

# --- TAKEN FROM arXiv:2309.08553v1, Further.tex lines 34-77 -----------------
P1 = [{1, 3, 15, 22, 30, 33, 46, 50, 55, 56}, {2, 9, 17, 21, 26, 27, 47, 49, 60, 63},
      {4, 7, 13, 23, 24, 28, 36, 54, 61}, {5, 12, 14, 20, 34, 44, 45, 57, 62},
      {6, 10, 16, 29, 38, 41, 43, 58, 59}, {8, 25, 31, 32, 35, 40, 51, 53, 65},
      {11, 18, 19, 37, 39, 42, 48, 52, 64}]
P2 = [{1, 3, 15, 22, 30, 33, 46, 50, 55, 56}, {2, 6, 14, 24, 27, 29, 38, 57, 58, 64},
      {5, 8, 23, 28, 39, 40, 47, 49, 53}, {10, 12, 13, 25, 34, 41, 45, 51, 59},
      {7, 16, 20, 21, 36, 42, 44, 54, 61}, {9, 11, 17, 32, 43, 48, 52, 62, 65},
      {4, 18, 19, 26, 31, 35, 37, 60, 63}]
P3 = [{2, 4, 16, 23, 31, 34, 47, 51, 56, 57}, {7, 10, 15, 19, 33, 43, 44, 50, 63, 65},
      {5, 12, 18, 21, 29, 39, 54, 58, 59}, {9, 11, 17, 28, 37, 40, 41, 55, 62},
      {8, 13, 14, 22, 24, 42, 45, 49, 64}, {3, 20, 26, 27, 30, 35, 46, 48, 60},
      {1, 6, 25, 32, 36, 38, 52, 53, 61}]
P4 = [{2, 4, 16, 23, 31, 34, 47, 51, 56, 57}, {3, 9, 10, 29, 38, 40, 43, 53, 61, 65},
      {5, 8, 22, 24, 37, 44, 45, 49, 55}, {14, 21, 27, 32, 35, 36, 52, 62, 64},
      {7, 11, 12, 26, 39, 42, 48, 50, 60}, {1, 13, 19, 20, 28, 30, 33, 54, 58},
      {6, 15, 17, 18, 25, 41, 46, 59, 63}]
P5 = [{3, 5, 17, 24, 32, 35, 48, 52, 57, 58}, {2, 6, 14, 16, 19, 37, 38, 44, 53, 64},
      {7, 10, 25, 30, 41, 42, 49, 51, 55}, {1, 13, 15, 18, 31, 39, 40, 46, 50},
      {4, 21, 22, 26, 28, 36, 47, 56, 59}, {8, 11, 23, 27, 29, 34, 54, 62, 63},
      {9, 12, 20, 33, 43, 45, 60, 61, 65}]
FIVE = [P1, P2, P3, P4, P5]

RESULTS = []    # list of (bool, str) filled by check()


def check(ok, label):
    """Record one check and print its PASS/FAIL line immediately."""
    RESULTS.append((bool(ok), label))
    print(("PASS " if ok else "FAIL ") + label)
    return bool(ok)


def say(*parts):
    """Print a derived intermediate quantity (indented, never a verdict line)."""
    print("    " + " ".join(str(p) for p in parts))


def sigma(A):
    """Sigma(A) = {a+b : a,b in A, a <= b}, derived from the definition."""
    s = sorted(A)
    out = set()
    for i in range(len(s)):
        for j in range(i, len(s)):
            out.add(s[i] + s[j])
    return out


def rep_counts(A):
    """Multiset of sums a+b with a <= b: sum -> number of representations.

    The paper's definition of Sidon is "every integer has at most one
    representation a+b with a,b in A, a <= b", so the maximum value of this
    counter is 1 exactly when A is Sidon.  This is derived independently of the
    |Sigma(A)| = C(|A|+1,2) criterion.
    """
    s = sorted(A)
    c = Counter()
    for i in range(len(s)):
        for j in range(i, len(s)):
            c[s[i] + s[j]] += 1
    return c


def diffs(A):
    """List of the C(|A|,2) positive differences b - a with a < b in A."""
    s = sorted(A)
    return [b - a for a, b in combinations(s, 2)]


def sidon_by_sigma(A):
    """Sidon test via |Sigma(A)| == C(|A|+1, 2)  (the paper's criterion)."""
    return len(sigma(A)) == comb(len(A) + 1, 2)


def sidon_by_reps(A):
    """Sidon test straight from the definition: no sum has two representations."""
    c = rep_counts(A)
    return (max(c.values()) if c else 0) <= 1


def sidon_by_diffs(A):
    """Sidon test via distinctness of the C(|A|,2) positive differences."""
    d = diffs(A)
    return len(d) == len(set(d))


def rho(A):
    """Reflection x -> 66 - x applied elementwise to a set."""
    return {REFL - x for x in A}


def rho_partition(P):
    """Reflection applied to every part of a partition."""
    return [rho(p) for p in P]


def canon(P):
    """Canonical form of an unordered family of parts."""
    return frozenset(frozenset(p) for p in P)


def show(P):
    """Readable canonical rendering of a partition: sorted list of sorted parts."""
    return sorted(sorted(p) for p in P)


def analyse(P):
    """Derive every property of a candidate balanced Sidon-Ramsey partition.

    Nothing here is assumed: sizes, disjointness, the union, the balance
    condition and the Sidon property of each part are all computed from P.
    """
    flat = [x for p in P for x in p]
    sizes = sorted(len(p) for p in P)
    d = {
        "nparts": len(P),
        "listed": len(flat),                       # counts repeats too
        "distinct": len(set(flat)),
        "sizes": sizes,
        "union_ok": set(flat) == set(range(1, N + 1)),
        "no_repeat": len(flat) == len(set(flat)),
        "balanced": (max(sizes) - min(sizes) <= 1) if sizes else False,
        "sigma_sizes": [len(sigma(p)) for p in P],
        "sigma_target": [comb(len(p) + 1, 2) for p in P],
        "sidon_sigma": [sidon_by_sigma(p) for p in P],
        "sidon_reps": [sidon_by_reps(p) for p in P],
        "sidon_diffs": [sidon_by_diffs(p) for p in P],
        "ndiffs": [len(diffs(p)) for p in P],
        "min_val": min(flat) if flat else None,
        "max_val": max(flat) if flat else None,
    }
    d["all_sidon"] = all(d["sidon_sigma"])
    d["criteria_agree"] = (d["sidon_sigma"] == d["sidon_reps"] == d["sidon_diffs"])
    d["is_bsrp"] = (d["nparts"] == K and d["listed"] == N and d["no_repeat"]
                    and d["union_ok"] and d["balanced"] and d["all_sidon"])
    return d


def all_profiles(nparts, total):
    """Every non-decreasing tuple of `nparts` positive integers summing to
    `total`, i.e. every partition of `total` into exactly `nparts` parts.

    Exhaustive: the recursion emits each such tuple exactly once.  The only
    pruning is the loop bound `v * slots <= rest`, which is exactly the
    condition that the remaining `slots` entries can all be >= v and still sum
    to `rest`; so no admissible tuple is skipped and every leaf is admissible.
    The balance condition is NOT used to prune, so filtering the output
    afterwards is a genuine search rather than a restatement of 65 = 7*9 + 2.
    """
    out = []

    def rec(prefix, lo, slots, rest):
        if slots == 1:
            if rest >= lo:
                out.append(tuple(prefix) + (rest,))
            return
        v = lo
        while v * slots <= rest:
            prefix.append(v)
            rec(prefix, v, slots - 1, rest - v)
            prefix.pop()
            v += 1

    rec([], 1, nparts, total)
    return out


def balanced_profiles(nparts, total):
    """The subset of all_profiles(nparts, total) that is balanced."""
    return [p for p in all_profiles(nparts, total) if max(p) - min(p) <= 1]


def report_partition(name, P, verbose=False):
    """Derive and print the properties of P; emit the checks for it.

    Returns the dict from analyse(P) so callers can reuse the derived values.
    """
    d = analyse(P)
    say(name + ": parts =", d["nparts"],
        "| integers listed =", d["listed"],
        "| distinct =", d["distinct"],
        "| min =", d["min_val"], "| max =", d["max_val"])
    say(name + ": derived size profile =", d["sizes"])
    say(name + ": derived |Sigma(A_i)| (parts in given order) =", d["sigma_sizes"])
    say(name + ": required C(|A_i|+1,2)                       =", d["sigma_target"])
    say(name + ": derived difference counts C(|A_i|,2) =", d["ndiffs"],
        "| total =", sum(d["ndiffs"]))
    if verbose:
        for i, p in enumerate(P, 1):
            c = rep_counts(p)
            say("  part %d (|A|=%d): %s" % (i, len(p), sorted(p)))
            say("    |Sigma| = %d, target %d, max sum-multiplicity = %d, "
                "differences %d all distinct = %s"
                % (len(sigma(p)), comb(len(p) + 1, 2), max(c.values()),
                   len(diffs(p)), sidon_by_diffs(p)))
    check(d["no_repeat"] and d["listed"] == N,
          "%s: the %d listed integers are pairwise distinct (derived %d listed, "
          "%d distinct)" % (name, N, d["listed"], d["distinct"]))
    check(d["union_ok"],
          "%s: union of the parts is exactly [1..%d]" % (name, N))
    check(d["nparts"] == K, "%s: exactly %d parts" % (name, K))
    check(d["balanced"],
          "%s: balanced -- max size - min size = %d <= 1"
          % (name, max(d["sizes"]) - min(d["sizes"])))
    check(d["all_sidon"],
          "%s: every part is Sidon by |Sigma(A)| = C(|A|+1,2)" % name)
    check(d["criteria_agree"] and all(d["sidon_reps"]),
          "%s: Sidon reconfirmed from the definition (no sum has two "
          "representations) and by distinctness of differences; all three "
          "criteria agree part by part" % name)
    return d


def section1_W():
    """(1),(2): W is a balanced Sidon-Ramsey partition of [65]; the paper's
    printed table row and difference count are reproduced from scratch."""
    print("=== (1)-(2)  the exhibited partition W = {A_1,...,A_7} ===")
    d = report_partition("W", W, verbose=True)

    say("paper's printed |A_i| row      =", PAPER_SIZE_ROW)
    say("derived |A_i| row (same order) =", [len(p) for p in W])
    check([len(p) for p in W] == PAPER_SIZE_ROW,
          "W: derived part sizes equal the paper's table row 10,10,9,9,9,9,9")
    check(sorted(d["sizes"]) == sorted(PAPER_SIZE_ROW) == [9, 9, 9, 9, 9, 10, 10],
          "W: derived size multiset is {10,10,9,9,9,9,9}")

    say("paper's printed |Sigma(A_i)| row =", PAPER_SIGMA_ROW)
    say("derived |Sigma(A_i)| row         =", d["sigma_sizes"])
    check(d["sigma_sizes"] == PAPER_SIGMA_ROW,
          "W: derived |Sigma(A_i)| row equals the paper's table row "
          "55,55,45,45,45,45,45")
    check(d["sigma_sizes"] == d["sigma_target"],
          "W: |Sigma(A_i)| = C(|A_i|+1,2) for every i, so every A_i is Sidon")

    total = sum(d["ndiffs"])
    say("derived per-part difference counts =", d["ndiffs"])
    say("derived total difference count =", total,
        "| paper's stated total =", PAPER_DIFF_TOTAL)
    check(total == PAPER_DIFF_TOTAL,
          "W: the derived total number of positive pairwise differences is %d, "
          "matching the paper's stated %d" % (total, PAPER_DIFF_TOTAL))
    check(d["ndiffs"] == [45, 45, 36, 36, 36, 36, 36],
          "W: per-part difference counts are 45,45 and 36 five times, as the "
          "paper states")
    check(all(d["sidon_diffs"]),
          "W: in every part the positive pairwise differences are distinct")
    check(d["is_bsrp"],
          "W is a balanced Sidon-Ramsey partition of [65] into 7 Sidon sets")
    return d


def section2_profiles():
    """(3): the size profile (9,9,9,9,9,10,10) is the only balanced one, so the
    paper's "65 = 7*9+2 forces two 10s and five 9s" is exhaustively confirmed."""
    print("=== (3)  balanced size profiles of 7 parts summing to 65 ===")
    allp = all_profiles(K, N)
    bal = balanced_profiles(K, N)
    say("partitions of 65 into exactly 7 positive parts (enumerated) =", len(allp))
    say("of those, balanced (max - min <= 1) =", len(bal), "->", bal)
    check(len(bal) == 1,
          "exactly one balanced size profile exists for 7 parts summing to 65 "
          "(derived count = %d)" % len(bal))
    check(bal == [(9, 9, 9, 9, 9, 10, 10)],
          "that unique balanced profile is (9,9,9,9,9,10,10)")
    check(len(allp) > 1 and any(max(p) - min(p) > 1 for p in allp),
          "the enumeration was not vacuous: %d profiles in total, of which %d "
          "are unbalanced" % (len(allp), len(allp) - len(bal)))
    return bal


def section3_five():
    """(4): each of the five displayed EGP partitions is itself a balanced
    Sidon-Ramsey partition of [65] with sizes (10,10,9,9,9,9,9)."""
    print("=== (4)  the five partitions displayed on p.6 of arXiv:2309.08553v1 ===")
    ds = []
    for i, P in enumerate(FIVE, 1):
        d = report_partition("P%d" % i, P)
        check(d["is_bsrp"],
              "P%d is a balanced Sidon-Ramsey partition of [65] into 7 Sidon sets"
              % i)
        check(d["sizes"] == [9, 9, 9, 9, 9, 10, 10],
              "P%d has size multiset {10,10,9,9,9,9,9} (derived %s)"
              % (i, d["sizes"]))
        ds.append(d)
    check(len({canon(P) for P in FIVE}) == 5,
          "the five displayed partitions are pairwise distinct as unordered "
          "families (derived %d distinct)" % len({canon(P) for P in FIVE}))
    return ds


def section4_family():
    """(5): build the conjectured family F = five displayed partitions together
    with their reflections, and establish |F| = 10 and rho(F) = F."""
    print("=== (5)  the conjectured family F = FIVE union rho(FIVE) ===")
    direct = {canon(P) for P in FIVE}
    reflected = {canon(rho_partition(P)) for P in FIVE}
    F = direct | reflected
    say("distinct canonical forms among the five displayed  =", len(direct))
    say("distinct canonical forms among their reflections   =", len(reflected))
    say("overlap (a displayed partition equal to a reflected one) =",
        len(direct & reflected))
    say("|F| =", len(F))
    check(len(F) == 10,
          "|F| = %d, i.e. the ten members of the conjectured family are "
          "pairwise distinct (expected 10)" % len(F))
    # The falsifiable content behind |F| = 10, checked rather than printed:
    # the five reflections are pairwise distinct AND none of them coincides with
    # a displayed partition.  (A self-reflective display would give |F| < 10.)
    check(len(reflected) == 5,
          "the five reflected partitions are pairwise distinct (derived %d)"
          % len(reflected))
    check(not (direct & reflected),
          "no displayed partition equals a reflected one, so the reflections "
          "are five genuinely new members (derived overlap %d)"
          % len(direct & reflected))
    F_reflected = {canon(rho_partition([set(p) for p in P])) for P in F}
    check(F_reflected == F,
          "F is rho-invariant: reflecting every member of F reproduces F "
          "exactly.  NOTE: this holds by construction -- F = D union rho(D) and "
          "rho is an involution, so rho(F) = rho(D) union D = F for ANY input; "
          "the check confirms the code, not the paper.  The falsifiable content "
          "is the two checks above")
    bad = []
    for P in sorted(F, key=lambda Q: show(Q)):
        d = analyse([set(p) for p in P])
        if not (d["is_bsrp"] and d["sizes"] == [9, 9, 9, 9, 9, 10, 10]):
            bad.append(show(P))
    check(len(F) == 10 and not bad,
          "all 10 members of F are balanced Sidon-Ramsey partitions of [65] "
          "with sizes (10,10,9,9,9,9,9) (derived over %d members)" % len(F))
    if bad:
        for b in bad:
            say("  offending member:", b)
    # rho is an involution of [65], a property of the reflection itself.
    ground = set(range(1, N + 1))
    say("rho([1..65]) has", len(rho(ground)), "elements; equals [1..65]:",
        rho(ground) == ground)
    check(rho(ground) == ground and all(rho(rho({x})) == {x} for x in ground),
          "rho(x) = 66 - x maps [1..65] onto itself and is an involution "
          "(identity in the integers: no input could make this fail)")
    parts_here = [p for P in FIVE + [W] for p in P]
    check(len(parts_here) == 6 * K
          and all(sidon_by_sigma(rho(p)) == sidon_by_sigma(p)
                  for p in parts_here)
          and all(sidon_by_sigma(p) for p in parts_here),
          "rho preserves the Sidon property on each of the %d parts occurring "
          "here, all %d of which are themselves Sidon (theorem, so this cannot "
          "fail either; the part count is asserted so the quantifier is not "
          "vacuous)" % (len(parts_here), len(parts_here)))
    return F


def section5_membership(F):
    """(6): the direct, shortcut-free refutation -- neither W nor rho(W) is a
    member of F, and the two are different partitions."""
    print("=== (6)  direct membership test of W and rho(W) against F ===")
    RW = rho_partition(W)
    cW, cRW = canon(W), canon(RW)
    say("W        =", show(W))
    say("rho(W)   =", show(RW))
    say("canon(W) in F      ->", cW in F)
    say("canon(rho(W)) in F ->", cRW in F)
    check(cW not in F,
          "W is NOT a member of the conjectured family F (direct set "
          "comparison against all 10 members)")
    check(cRW not in F,
          "rho(W) is NOT a member of the conjectured family F")
    check(cRW != cW,
          "rho(W) != W, so the counterexample paper exhibits two distinct "
          "omitted partitions")
    dR = report_partition("rho(W)", RW)
    check(dR["is_bsrp"] and dR["sizes"] == [9, 9, 9, 9, 9, 10, 10],
          "rho(W) is itself a balanced Sidon-Ramsey partition of [65] with "
          "sizes (10,10,9,9,9,9,9)")
    # How close does W come to any family member?  Purely informational, but it
    # shows the non-membership is not a near miss hidden by canonicalisation.
    best, arg = -1, None
    for P in F:
        ov = len(cW & P)
        if ov > best:
            best, arg = ov, P
    say("maximum number of parts W shares with any member of F =", best,
        "out of", K)
    say("closest member of F to W =", show(arg))
    say("parts of W absent from that member =",
        sorted(sorted(p) for p in (cW - arg)))
    check(best < K,
          "W shares at most %d of its %d parts with any single member of F"
          % (best, K))
    return cW, cRW


def section6_shortcut(F):
    """(7): the paper's own argument -- no 10-element part of the five displayed
    partitions contains {2,65} or {1,64}, while {2,65} lies in A_2."""
    print("=== (7)  the paper's shortcut argument, re-derived ===")
    slots = [p for P in FIVE for p in P if len(p) == 10]
    distinct = {frozenset(p) for p in slots}
    say("10-element part slots across the five partitions =", len(slots))
    say("distinct 10-element sets among them             =", len(distinct))
    check(len(slots) == 10,
          "the five displayed partitions contribute exactly 10 ten-element "
          "part slots (derived %d)" % len(slots))
    check(len(distinct) == 8,
          "only %d of those 10 slots are distinct sets, as the paper says"
          % len(distinct))
    paper_eight = {frozenset(s) for s in PAPER_EIGHT}
    check(len(paper_eight) == 8,
          "the eight sets printed in the paper are themselves pairwise distinct")
    check(distinct == paper_eight,
          "the distinct 10-element parts derived from the EGP source are "
          "exactly the eight sets printed in the paper")
    if distinct != paper_eight:
        for s in sorted(distinct - paper_eight, key=sorted):
            say("  derived but not printed:", sorted(s))
        for s in sorted(paper_eight - distinct, key=sorted):
            say("  printed but not derived:", sorted(s))
    counts = Counter(frozenset(p) for p in slots)
    repeats = sorted((sorted(s), c) for s, c in counts.items() if c > 1)
    for s, c in repeats:
        say("repeated 10-element part (occurs %d times):" % c, s)
    check([s for s, _ in repeats] == sorted(sorted(x) for x in PAPER_EIGHT[:2]),
          "the repeated 10-element parts are exactly the first two of the "
          "paper's eight, read row-major, each occurring twice")
    check(all(c == 2 for _, c in repeats) and len(repeats) == 2,
          "each of the two repeated parts occurs exactly twice (10 slots = "
          "8 distinct + 2 repeats)")
    # The paper says more than "they repeat": "the first is a part of the first
    # and of the second displayed partition, the second a part of the third and
    # of the fourth".  Derive the owning display indices and compare.
    owners = [[i for i, P in enumerate(FIVE, 1) if any(sorted(p) == s for p in P)]
              for s, _ in repeats]
    say("displayed partitions owning each repeated 10-element part =", owners)
    check(owners == [[1, 2], [3, 4]],
          "the first repeated part belongs to displays 1 and 2 and the second "
          "to displays 3 and 4, exactly as the paper states (derived %s)"
          % owners)
    return slots, distinct


def section7_pairs(F, slots):
    """(7, continued): the endpoint-pair obstruction that the paper actually
    uses, and its consequence re-derived over the whole family F."""
    print("=== (7b)  the {2,65} / {1,64} obstruction ===")
    pair, rpair = {2, 65}, {1, 64}
    say("rho({2,65}) =", sorted(rho(pair)), "so {1,64} is the reflected pair")
    check(rho(pair) == rpair,
          "rho({2,65}) = {1,64}, so the two pair conditions are reflections of "
          "each other")
    with_pair = [sorted(p) for p in slots if pair <= p]
    with_rpair = [sorted(p) for p in slots if rpair <= p]
    say("displayed 10-element slots containing {2,65} =", with_pair)
    say("displayed 10-element slots containing {1,64} =", with_rpair)
    check(not with_pair,
          "no 10-element part of the five displayed partitions contains {2,65}")
    check(not with_rpair,
          "no 10-element part of the five displayed partitions contains {1,64}")
    # Consequence, derived over all ten family members rather than argued.
    fam_pair = [show(P) for P in F
                if any(len(p) == 10 and pair <= set(p) for p in P)]
    fam_rpair = [show(P) for P in F
                 if any(len(p) == 10 and rpair <= set(p) for p in P)]
    say("members of F having a 10-element part containing {2,65} =", len(fam_pair))
    say("members of F having a 10-element part containing {1,64} =", len(fam_rpair))
    check(not fam_pair,
          "no member of F (displayed or reflected) has a 10-element part "
          "containing {2,65}")
    check(not fam_rpair,
          "no member of F has a 10-element part containing {1,64}")
    A2 = W[1]
    say("A_2 =", sorted(A2), "| |A_2| =", len(A2))
    check(len(A2) == 10 and pair <= A2,
          "{2,65} is contained in the 10-element part A_2 of W, hence W cannot "
          "be a member of F")
    say("derived rho(A_2)      =", sorted(rho(A2)))
    say("hand-computed value   =", sorted(HAND_RHO_A2))
    check(rho(A2) == HAND_RHO_A2,
          "derived rho(A_2) equals {1,6,8,16,22,33,42,45,46,64}, the set "
          "{66-a : a in A_2} reflected by hand from the paper's A_2 display; "
          "this cross-checks the code's rho, and a referee can redo it by hand")
    check(rpair <= rho(A2),
          "{1,64} is contained in rho(A_2), hence rho(W) cannot be a member "
          "of F either")
    tens_W = [p for p in W if len(p) == 10]
    say("10-element parts of W =", [sorted(p) for p in tens_W],
        "| any containing {1,64}:", any(rpair <= p for p in tens_W))
    check(len(tens_W) == 2 and not any(rpair <= p for p in tens_W),
          "neither 10-element part of W contains {1,64}, which is the paper's "
          "reason for rho(W) != W")


def parse_five_from_source(text):
    """Recover the five displayed partitions FROM the EGP source text.

    Returns (families, note).  `families` is a list, one entry per align*
    display found strictly between the two delimiting sentences, each entry a
    list of sets.  Nothing is looked up from FIVE, so comparing the result with
    FIVE tests the transcription in BOTH directions: every transcribed part is
    in the source (inclusion) and the source region holds nothing else -- no
    sixth partition and no extra part (completeness).  The region matters: the
    file contains a further align* display after the conjecture, holding the
    NON-balanced partition, which must not be swept in.
    """
    if text.count(SRC_REGION_START) != 1:
        return None, "opening sentence occurs %d times, expected 1" % \
            text.count(SRC_REGION_START)
    if text.count(SRC_REGION_END) != 1:
        return None, "closing sentence occurs %d times, expected 1" % \
            text.count(SRC_REGION_END)
    i = text.index(SRC_REGION_START)
    j = text.index(SRC_REGION_END)
    if not i < j:
        return None, "the two delimiting sentences occur in the wrong order"
    region = text[i:j]
    blocks = re.findall(r"\\begin\{align\*\}(.*?)\\end\{align\*\}", region,
                        re.S)
    families = []
    for blk in blocks:
        flat = blk.replace(" ", "").replace("\t", "").replace("\n", "")
        lits = re.findall(r"\\\{(\d+(?:,\d+)*)\\\}", flat)
        families.append([set(int(x) for x in lit.split(",")) for lit in lits])
    return families, "%d align* displays in the region" % len(blocks)


def optional_source_check(path):
    """Optional (8): re-verify the transcription against the EGP source itself.

    Run as `python3 verify.py /path/to/Further.tex`, where Further.tex is the
    file from the arXiv:2309.08553v1 e-print tarball
    (https://arxiv.org/e-print/2309.08553v1).  Without an argument the program
    still performs every mathematical check; only this provenance section is
    skipped, and it is announced as skipped rather than passed.
    """
    import hashlib
    print("=== (8)  optional: transcription re-checked against the EGP source ===")
    try:
        with open(path, "rb") as fh:
            raw = fh.read()
    except OSError as exc:
        # A path was asked for and could not be read: that is a failure, not a
        # skip, but it must not abort the run before the verdict line.
        check(False, "the supplied source path %r could be opened and read "
                     "(OSError: %s)" % (path, exc))
        return
    say("file =", path, "| bytes =", len(raw))
    say("md5    =", hashlib.md5(raw).hexdigest())
    say("sha256 =", hashlib.sha256(raw).hexdigest())
    check(hashlib.md5(raw).hexdigest() == KNOWN_EGP_MD5,
          "supplied Further.tex is byte-identical to the copy these five "
          "partitions were transcribed from (md5 %s)" % KNOWN_EGP_MD5)
    check(hashlib.sha256(raw).hexdigest() == KNOWN_EGP_SHA256,
          "the same file matches the pinned sha256 %s (md5 alone is "
          "collision-broken)" % KNOWN_EGP_SHA256)
    text = raw.decode("utf-8", "replace")
    setform = text.replace("\\", "").replace(" ", "").replace("\n", "")
    missing = []
    for i, P in enumerate(FIVE, 1):
        for p in P:
            lit = "{" + ",".join(str(x) for x in sorted(p)) + "}"
            if lit not in setform:
                missing.append((i, lit))
    say("part literals sought in the source =", sum(len(P) for P in FIVE),
        "| not found =", len(missing))
    for i, lit in missing:
        say("  P%d: literal absent from source:" % i, lit)
    check(not missing,
          "all 35 transcribed parts of the five displayed partitions occur "
          "verbatim (as increasing lists) in the EGP source")
    # The check above is one-directional: it cannot notice a SIXTH displayed
    # partition, or a part dropped from the transcription.  Either would shrink
    # F and could let a genuine member of the conjectured family pass as a
    # counterexample.  So parse the source region and require equality.
    parsed, note = parse_five_from_source(text)
    say("source region parse:", note)
    check(parsed is not None,
          "the EGP source region holding the five displays was located "
          "unambiguously (both delimiting sentences occur exactly once, in "
          "order)")
    if parsed is None:
        return
    say("displays found in the region =", len(parsed),
        "| parts per display =", [len(P) for P in parsed],
        "| parts in total =", sum(len(P) for P in parsed))
    check(len(parsed) == 5,
          "the source region displays exactly 5 partitions -- there is no "
          "sixth balanced partition on that page that the transcription "
          "omitted (derived %d)" % len(parsed))
    check([len(P) for P in parsed] == [7, 7, 7, 7, 7]
          and sum(len(P) for P in parsed) == 35,
          "each of those 5 displays holds exactly 7 parts, 35 in all (derived "
          "%s)" % [len(P) for P in parsed])
    src_form = [[sorted(p) for p in P] for P in parsed]
    our_form = [[sorted(p) for p in P] for P in FIVE]
    if src_form != our_form:
        for k in range(max(len(src_form), len(our_form))):
            a = src_form[k] if k < len(src_form) else None
            b = our_form[k] if k < len(our_form) else None
            if a != b:
                say("  display %d: source =" % (k + 1), a)
                say("  display %d: program =" % (k + 1), b)
    check(src_form == our_form,
          "FIVE as transcribed equals the five partitions parsed out of the "
          "source, display by display and part by part, in order -- so the "
          "transcription is complete as well as correct, and F is the whole "
          "conjectured family")
    alnum = "".join(ch for ch in text if ch.isalnum())
    for frag, desc in QUOTE_FRAGMENTS:
        key = "".join(ch for ch in frag if ch.isalnum())
        check(key in alnum,
              "the quoted %s occurs in the EGP source (compared after "
              "discarding non-alphanumeric characters)" % desc)


def main():
    print("verify.py -- counterexample to the Espinosa-Garcia / Pellicer")
    print("exhaustiveness conjecture for balanced Sidon-Ramsey partitions of [65]")
    print("Exact integer arithmetic only; no floating point, no tolerances.")
    print("")
    section1_W()
    print("")
    section2_profiles()
    print("")
    section3_five()
    print("")
    F = section4_family()
    print("")
    section5_membership(F)
    print("")
    slots, _ = section6_shortcut(F)
    print("")
    section7_pairs(F, slots)
    print("")
    if len(sys.argv) > 1:
        optional_source_check(sys.argv[1])
    else:
        print("=== (8)  optional source re-check: SKIPPED (no argument) ===")
        print("    Pass the path to Further.tex from the arXiv:2309.08553v1")
        print("    e-print tarball (https://arxiv.org/e-print/2309.08553v1) to")
        print("    re-verify the transcription of the five displayed partitions")
        print("    and of the quoted conjecture against the source itself:")
        print("        python3 verify.py /path/to/Further.tex")
        print("    The mathematical checks above do not depend on it.")
    print("")
    print("=== theorem ===")
    dW = analyse(W)
    dRW = analyse(rho_partition(W))
    check(dW["is_bsrp"] and dW["sizes"] == [9, 9, 9, 9, 9, 10, 10]
          and canon(W) not in F
          and dRW["is_bsrp"] and dRW["sizes"] == [9, 9, 9, 9, 9, 10, 10]
          and canon(rho_partition(W)) not in F
          and canon(rho_partition(W)) != canon(W),
          "THEOREM: W and rho(W) are two distinct balanced Sidon-Ramsey "
          "partitions of [65] with two parts of size 10 and five of size 9, "
          "and neither is among the five partitions displayed in "
          "arXiv:2309.08553v1 or their reflections -- so the exhaustiveness "
          "conjecture quoted on p.6 of that e-print is FALSE, subject to the "
          "three non-machine-checkable steps listed in the gaps block below")
    print("")
    print("=== gaps: steps between the checked facts and the paper's claim ===")
    print("    The checks above establish, exactly and without tolerance, that W")
    print("    and rho(W) are two distinct balanced Sidon-Ramsey partitions of")
    print("    [65] with profile (10,10,9,9,9,9,9) lying outside the ten-member")
    print("    family F.  Three steps from there to the paper's sentence are NOT")
    print("    covered by any check above and are not counted as passing:")
    print("")
    print("    (i)  READING OF THE QUOTED SENTENCE.  'We conjecture that these")
    print("         and their reflected partitions ... are all balanced Sidon-")
    print("         Ramsey partitions with those parameters' is read as an")
    print("         exhaustiveness claim (these are all of them).  Under the")
    print("         weak reading (each of these is one) the sentence is true and")
    print("         the refutation does not apply.  No program can settle this.")
    print("         Two textual facts favour the strong reading and a reader")
    print("         should weigh them: the preceding sentence already gives the")
    print("         parameters, so the weak reading is redundant; and the")
    print("         paragraph opens 'So far we have found $5$ ...', i.e. it is")
    print("         about how many such partitions exist.")
    if len(sys.argv) > 1:
        print("    (ii) TRANSCRIPTION OF THE FIVE DISPLAYED PARTITIONS: closed by")
        print("         section (8) above, which parsed them out of the supplied")
        print("         source and required equality with the transcription in")
        print("         both directions (no sixth display, no dropped part).")
    else:
        print("    (ii) TRANSCRIPTION OF THE FIVE DISPLAYED PARTITIONS.  F was")
        print("         built from partitions transcribed by hand from")
        print("         arXiv:2309.08553v1, Further.tex lines 34-77.  With no")
        print("         source file supplied, NOTHING here rules out a sixth")
        print("         displayed partition, or a dropped part, which would make")
        print("         F too small and could admit a genuine family member as a")
        print("         counterexample.  Run")
        print("             python3 verify.py /path/to/Further.tex")
        print("         to close this gap; section (8) then checks it in both")
        print("         directions, not merely that each transcribed part occurs.")
    print("    (iii) VERSION OF RECORD.  Whether the quoted sentence survives")
    print("         into Discrete Appl. Math. 378 (2026), 120-124,")
    print("         doi:10.1016/j.dam.2025.07.002 cannot be decided here: the")
    print("         article is paywalled.  The paper scopes its claim to")
    print("         arXiv:2309.08553v1, which is what is checked.  A human with")
    print("         journal access should read around p.123.")
    print("")
    failed = [lab for ok, lab in RESULTS if not ok]
    n = len(RESULTS)
    if failed:
        print("failing checks:")
        for lab in failed:
            print("    FAIL " + lab)
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())

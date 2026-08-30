#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent verification of the 5^1 2^30-CURD on 65 varieties printed in paper.tex.

Python 3.9+, STANDARD LIBRARY ONLY (hashlib, itertools, sys, collections).  No third-party package and no
external data file: the only input is the block of four base parallel classes reproduced verbatim below
from Section 2 of the paper, so a referee can run this with nothing installed and compare the string
against the printed page (the program prints its SHA-256).

All arithmetic is over the integers.  There is no floating-point decision anywhere in this file: the two
places where a ratio appears (the source's w and t formulas) are checked as exact divisions with
divmod().

What is checked
--------------
  * the source's own necessary conditions at (n, m) = (65, 5): w = 52 parallel classes, t = 4 five-blocks
    through a variety, mw = tn, the two bounds of the source's Theorem 2.5, and the block census;
  * the four printed base classes: each is a partition of the 65 varieties of type 5^1 2^30, and each
    five-block is the transversal predicted by the difference-matrix description in the paper;
  * the sigma-orbit argument that decides lambda = 1: sigma is an order-13 permutation, every pair orbit
    has size exactly 13, there are exactly 160 orbits, the base classes supply exactly 160 pairs, and
    their 160 orbit labels are pairwise distinct -- plus the hand-sized half of that argument, that the
    four five-blocks already use four distinct cross-labels on every column pair;
  * the whole development, checked DIRECTLY and independently of the orbit argument: 52 parallel classes,
    1612 blocks, every class a partition of type 5^1 2^30, all 2080 pairs covered with multiplicity
    multiset {1: 2080}, t = 4 at every variety, 52 distinct five-blocks pairwise meeting in at most one
    variety, and each family of thirteen five-blocks itself a parallel class;
  * two anti-controls, so that a PASS is not vacuous: a corrupted copy of the object and a copy with a
    repeated base class are both REJECTED by the same lambda test, with the damage named;
  * the arithmetic that distinguishes the 2004 near miss (Danziger-Stevens, EJC 11 #R23, Theorem 3.5 at
    (q, k, m) = (13, 5, 3)) from this object, and the published two-block-size bound v <= p_2^2;
  * the bookkeeping of the claim: which cell is closed, and that the source's other open cell has
    different parameters.

What is NOT re-run is printed by the program itself at the end,.
"""

import hashlib
import itertools
import sys
from collections import Counter, defaultdict

Q = 13          # the order of the cyclic group acting inside each column
M = 5           # the single large block size
COLS = 5        # number of columns; n = COLS * Q
N = COLS * Q    # 65 varieties
T = 4           # number of base classes = number of five-blocks through a variety

# ---------------------------------------------------------------------------
# THE OBJECT, copied verbatim from Section 2 of paper.tex.
# A variety is written c_x and means the integer u = 13*c + x, with c in 0..4 and x in Z_13.
# `Pk 5-BLOCK ...` gives the five-block of the base class P_k; the `Pk PAIRS` lines list its thirty
# 2-blocks, written a|b.  The design is the 52 classes sigma^s(P_k), k = 1..4, s = 0..12, where sigma
# adds 1 to x in every column.
# ---------------------------------------------------------------------------
OBJECT = """
P1 5-BLOCK 0_0 1_0 2_0 3_0 4_0
P1 PAIRS   0_1|0_8 0_2|0_12 0_3|4_12 0_4|1_10 0_5|3_2 0_6|3_4
P1 PAIRS   0_7|0_9 0_10|3_12 0_11|2_9 1_1|1_4 1_2|1_9 1_3|1_8
P1 PAIRS   1_5|3_8 1_6|2_1 1_7|3_3 1_11|2_7 1_12|3_11 2_2|2_12
P1 PAIRS   2_3|2_4 2_5|4_10 2_6|4_4 2_8|2_10 2_11|4_8 3_1|3_9
P1 PAIRS   3_5|4_3 3_6|3_7 3_10|4_7 4_1|4_2 4_5|4_11 4_6|4_9
P2 5-BLOCK 0_0 1_1 2_2 3_3 4_4
P2 PAIRS   0_1|1_6 0_2|4_9 0_3|2_6 0_4|3_8 0_5|1_12 0_6|2_7
P2 PAIRS   0_7|3_2 0_8|4_0 0_9|0_10 0_11|4_1 0_12|3_4 1_0|2_4
P2 PAIRS   1_2|3_10 1_3|4_2 1_4|2_3 1_5|1_9 1_7|2_0 1_8|2_5
P2 PAIRS   1_10|3_7 1_11|4_12 2_1|2_9 2_8|4_7 2_10|3_1 2_11|3_5
P2 PAIRS   2_12|3_11 3_0|3_9 3_6|4_11 3_12|4_5 4_3|4_8 4_6|4_10
P3 5-BLOCK 0_0 1_2 2_4 3_6 4_8
P3 PAIRS   0_1|0_10 0_2|1_6 0_3|4_1 0_4|0_9 0_5|2_10 0_6|1_4
P3 PAIRS   0_7|3_1 0_8|4_10 0_11|2_6 0_12|1_9 1_0|2_7 1_1|3_12
P3 PAIRS   1_3|4_0 1_5|2_3 1_7|4_11 1_8|3_0 1_10|3_4 1_11|4_3
P3 PAIRS   1_12|4_7 2_0|3_8 2_1|3_11 2_2|3_7 2_5|2_12 2_8|4_9
P3 PAIRS   2_9|3_2 2_11|4_6 3_3|3_10 3_5|4_12 3_9|4_5 4_2|4_4
P4 5-BLOCK 0_0 1_3 2_6 3_9 4_12
P4 PAIRS   0_1|1_10 0_2|2_1 0_3|2_10 0_4|4_1 0_5|2_2 0_6|4_7
P4 PAIRS   0_7|3_8 0_8|3_7 0_9|2_5 0_10|1_5 0_11|4_4 0_12|1_11
P4 PAIRS   1_0|1_1 1_2|2_7 1_4|4_11 1_6|1_8 1_7|4_9 1_9|3_10
P4 PAIRS   1_12|4_10 2_0|2_9 2_3|4_6 2_4|3_0 2_8|3_6 2_11|4_5
P4 PAIRS   2_12|4_8 3_1|3_11 3_2|3_4 3_3|4_2 3_5|4_0 3_12|4_3
"""

NOT_RE_RUN = [
    "NOT RE-RUN: the SEARCH. The object above was found by Algorithm X on a 360-element exact cover",
    "  (240 (class, variety) elements plus 120 leave orbits, 5,280 candidate rows), at seed 1 in 130",
    "  nodes. This program does not repeat that search and could not reproduce its choices; nothing here",
    "  depends on it, because a positive certificate of this kind is self-checking.",
    "NOT RE-VERIFIED: the 'only if' half of the source's Theorem 2.9, that a 5^1 2^((n-5)/2)-CURD forces",
    "  n in {5, 25, 65}. It is quoted, not re-proved. What is re-derived here is only its arithmetic at",
    "  n = 65 and the two bounds of the source's Theorem 2.5 at t = 4.",
    "NOT REPRODUCED: the source's own objects at n = 5 (its Example 1.2) and n = 25 (its Theorem 4.7 /",
    "  Corollary 4.11 at q = 5). Only their parameter arithmetic is re-derived below; no object at either",
    "  order is built, parsed or checked here, so the completed spectrum {5, 25, 65} rests on the source",
    "  for two of its three members.",
    "NOT COVERED: uniqueness, enumeration, or the number of solutions at n = 65. This is a construction,",
    "  not a classification; the exhausted fraction of the unrestricted space is 0, and no non-existence",
    "  claim is made or checked anywhere in this file.",
    "NOT COVERED: the source's OTHER open cell, a 4^1 2^12-CURD on n = 28 (its line 361). This program",
    "  re-derives only that its parameters differ from ours (m = 4, w = 21, t = 3) and asserts nothing",
    "  about its existence in either direction.",
    "NOT TRANSCRIBED: Inequality (6) of Danziger-Stevens (2004), their general two-block-size extremal",
    "  bound. Its text extraction is layout-mangled -- the fraction bar and the scope of the leading",
    "  lambda are unrecoverable -- so it was not transcribed with certainty and is not evaluated here.",
    "  Only the published bound p_2 > 1 => v <= p_2^2 is checked.",
    "NOT SEARCHED: prior art beyond the channels listed in the 'Bibliographic channels were incomplete'",
    "  paragraph of the accompanying note ",
    "  which full texts were unreachable. OpenAlex returned HTTP 429, MathSciNet and two Wiley full texts",
    "  were inaccessible, and the Handbook CURD section is print only. Nothing in this program bears on",
    "  those gaps.",
]

_passes = 0
_fails = 0


def check(name, ok, detail=""):
    global _passes, _fails
    if ok:
        _passes += 1
        print("PASS %s%s" % (name, (" " + detail) if detail else ""))
    else:
        _fails += 1
        print("FAIL %s%s" % (name, (" " + detail) if detail else ""))


# ---------------------------------------------------------------------------
# 0. helpers -- the coordinate system and the pair-orbit labels
# ---------------------------------------------------------------------------
def uid(tok):
    """'c_x' -> the integer 13*c + x, with bounds enforced."""
    c, x = tok.split("_")
    c, x = int(c), int(x)
    assert 0 <= c < COLS and 0 <= x < Q, tok
    return Q * c + x


def coords(u):
    return divmod(u, Q)          # (column, x)


def sigma(u):
    c, x = coords(u)
    return Q * c + (x + 1) % Q


def orbit_label(u, v):
    """The sigma-orbit of the pair {u, v}, as a hashable label.

    Same column c: ('S', c, d) with d = min(|x-y|, 13-|x-y|) in 1..6.
    Columns c < c':  ('C', c, c', (y-x) mod 13).
    """
    (c, x), (d, y) = coords(u), coords(v)
    if c == d:
        e = abs(x - y)
        return ("S", c, min(e, Q - e))
    if c > d:
        c, x, d, y = d, y, c, x
    return ("C", c, d, (y - x) % Q)


def all_orbit_labels():
    lab = set()
    for u, v in itertools.combinations(range(N), 2):
        lab.add(orbit_label(u, v))
    return lab


def parse(text):
    """-> [(five_block, [pair, ...]), ...] in k order, from the printed block only."""
    got = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        assert parts[0][0] == "P", line
        k = int(parts[0][1:])
        entry = got.setdefault(k, {"five": None, "pairs": []})
        if parts[1] == "5-BLOCK":
            assert entry["five"] is None, "two five-blocks for P%d" % k
            entry["five"] = tuple(sorted(uid(t) for t in parts[2:]))
        elif parts[1] == "PAIRS":
            for tok in parts[2:]:
                a, b = tok.split("|")
                entry["pairs"].append(tuple(sorted((uid(a), uid(b)))))
        else:
            raise AssertionError(line)
    return [(got[k]["five"], got[k]["pairs"]) for k in sorted(got)]


def develop(base):
    """base = [(five, pairs), ...] -> the list of parallel classes, each a list of blocks."""
    classes = []
    for five, pairs in base:
        for s in range(Q):
            sh = [u for u in range(N)]
            for _ in range(s):
                sh = [sigma(u) for u in sh]
            cl = [tuple(sorted(sh[u] for u in five))]
            for (u, v) in pairs:
                cl.append(tuple(sorted((sh[u], sh[v]))))
            classes.append(cl)
    return classes


def lambda_profile(classes):
    """The pair-multiplicity Counter over ALL C(65,2) pairs, zeros included."""
    mult = Counter()
    for cl in classes:
        for blk in cl:
            for p in itertools.combinations(blk, 2):
                mult[p] += 1
    prof = Counter()
    for p in itertools.combinations(range(N), 2):
        prof[mult[p]] += 1
    return mult, prof


def is_type_5_1_2_30(cl):
    flat = [u for blk in cl for u in blk]
    sizes = sorted(len(blk) for blk in cl)
    return (sorted(flat) == list(range(N))
            and sizes == [2] * 30 + [5])


# ---------------------------------------------------------------------------
# 1. the source's necessary conditions at (n, m) = (65, 5)
# ---------------------------------------------------------------------------
print("=== the cell: the source's own necessary conditions at n = 65, m = 5 ===")
n_pairs = N * (N - 1) // 2
check("n-m-and-the-pair-count", N == 65 and M == 5 and n_pairs == 2080,
      "n = %d varieties, m = %d, C(n,2) = %d pairs to cover exactly once" % (N, M, n_pairs))

num, den = N * (N - 1), M * M - 2 * M + N
w, rem = divmod(num, den)
check("w-from-the-first-necessary-condition", rem == 0 and w == 52,
      "w = n(n-1)/(m^2-2m+n) = %d/%d = %d parallel classes, an exact division" % (num, den, w))

num2, den2 = N - w - 1, M - 2
t, rem2 = divmod(num2, den2)
check("t-from-the-second-necessary-condition", rem2 == 0 and t == T == 4,
      "t = (n-w-1)/(m-2) = %d/%d = %d five-blocks through each variety, an exact division"
      % (num2, den2, t))

check("third-necessary-condition-mw-equals-tn", M * w == t * N == 260,
      "mw = %d and tn = %d" % (M * w, t * N))

bound_n = t ** 4 - t ** 3 + t ** 2 - t + 1
bound_m = t * (t * t - 2 * t + 2)
check("theorem-2-5-bounds-at-t-4", N <= bound_n and M <= bound_m,
      "n = %d <= t^4-t^3+t^2-t+1 = %d and m = %d <= t(t^2-2t+2) = %d" % (N, bound_n, M, bound_m))

per_class = M * (M - 1) // 2 + (N - M) // 2
check("block-and-pair-census-closes-exactly",
      per_class == 40 and w * per_class == n_pairs and w * 31 == 1612 and w * 30 == 1560,
      "each class covers %d pairs (10 inside the five-block, 30 pairs); %d x %d = %d = C(65,2); "
      "%d five-blocks + %d 2-blocks = %d blocks" % (per_class, w, per_class, w * per_class, w, w * 30,
                                                    w * 31))
check("partition-exponent-is-30", (N - M) // 2 == 30 and 2 * 30 + M == N,
      "the partition 5^1 2^((n-5)/2) reads 5^1 2^30 at n = 65")
check("n-65-is-one-of-the-three-admissible-orders", N in (5, 25, 65),
      "the source's Theorem 2.9 admits only n = 5, 25, 65; this is the third")

# the two orders the source itself settles -- arithmetic only, no object
for nn, ww, tt in ((5, 1, 1), (25, 15, 3)):
    a, ra = divmod(nn * (nn - 1), M * M - 2 * M + nn)
    b, rb = divmod(nn - a - 1, M - 2)
    ok = (ra == 0 and rb == 0 and a == ww and b == tt
          and a * (M * (M - 1) // 2 + (nn - M) // 2) == nn * (nn - 1) // 2)
    check("arithmetic-of-the-source-s-n-%d-cell" % nn, ok,
          "n = %d: w = %d, t = %d, and %d x %d = %d = C(%d,2) -- the source supplies the object, not us"
          % (nn, a, b, a, M * (M - 1) // 2 + (nn - M) // 2, nn * (nn - 1) // 2, nn))

# ---------------------------------------------------------------------------
# 2. the printed object
# ---------------------------------------------------------------------------
print("")
print("=== the four base parallel classes, read from the block printed in Section 2 ===")
base = parse(OBJECT)
check("four-base-classes-parsed", len(base) == T
      and all(len(f) == M and len(p) == 30 for f, p in base),
      "%d base classes, each one five-block of size 5 and 30 pairs" % len(base))

for k, (five, pairs) in enumerate(base, start=1):
    cl = [five] + pairs
    check("P%d-is-a-partition-of-type-5-1-2-30" % k, is_type_5_1_2_30(cl),
          "the 5 + 60 varieties listed are 0..64 each exactly once, block sizes 5, 2^30")

for k, (five, pairs) in enumerate(base, start=1):
    want = tuple(sorted(Q * c + ((k - 1) * c) % Q for c in range(COLS)))
    check("P%d-five-block-is-the-predicted-transversal" % k, five == want,
          "{13c + (%d)c mod 13 : c = 0..4} = %s, one variety per column"
          % (k - 1, "{" + ",".join("%d_%d" % coords(u) for u in five) + "}"))

# ---------------------------------------------------------------------------
# 3. sigma and the orbit argument that decides lambda = 1
# ---------------------------------------------------------------------------
print("")
print("=== sigma, the pair orbits, and the distinctness that certifies lambda = 1 ===")
img = [sigma(u) for u in range(N)]
ordr = 1
cur = img
while cur != list(range(N)):
    cur = [img[u] for u in cur]
    ordr += 1
check("sigma-is-a-permutation-of-order-13", sorted(img) == list(range(N)) and ordr == Q,
      "sigma: c_x -> c_(x+1 mod 13) permutes the 65 varieties and has order exactly %d" % ordr)

orbits = defaultdict(list)
for u, v in itertools.combinations(range(N), 2):
    orbits[orbit_label(u, v)].append((u, v))
check("every-pair-orbit-has-size-exactly-13",
      all(len(o) == Q for o in orbits.values()),
      "all %d orbits have size 13; no pair is fixed by any nontrivial power of sigma" % len(orbits))

same = sorted(l for l in orbits if l[0] == "S")
cross = sorted(l for l in orbits if l[0] == "C")
check("there-are-160-orbits-30-same-column-and-130-cross-column",
      len(orbits) == 160 and len(same) == 30 and len(cross) == 130,
      "5 columns x 6 distances = 30 same-column labels; C(5,2) x 13 = 130 cross-column labels")
check("160-orbits-times-13-is-the-whole-pair-set", 160 * Q == n_pairs,
      "160 x 13 = %d = C(65,2) = %d" % (160 * Q, n_pairs))

supplied = []
for five, pairs in base:
    for p in itertools.combinations(five, 2):
        supplied.append(orbit_label(*p))
    for (u, v) in pairs:
        supplied.append(orbit_label(u, v))
check("the-base-classes-supply-exactly-160-pairs", len(supplied) == T * per_class == 160,
      "4 classes x 40 pairs = %d pairs, exactly one per orbit if and only if the labels are distinct"
      % len(supplied))
check("the-160-orbit-labels-are-pairwise-distinct", len(set(supplied)) == 160,
      "%d distinct labels among %d supplied" % (len(set(supplied)), len(supplied)))
sup_same = [l for l in supplied if l[0] == "S"]
sup_cross = [l for l in supplied if l[0] == "C"]
check("same-column-labels-30-of-30", len(set(sup_same)) == len(sup_same) == 30,
      "every one of the 30 same-column orbits is used exactly once")
check("cross-column-labels-130-of-130", len(set(sup_cross)) == len(sup_cross) == 130,
      "every one of the 130 cross-column orbits is used exactly once")
check("lambda-1-certified-by-distinctness", set(supplied) == set(orbits),
      "the 160 supplied labels are exactly the 160 orbits, so developing by sigma covers each of the "
      "2080 pairs exactly once")

# the hand-sized half: the five-blocks alone already use 4 distinct cross-labels per column pair
hand = True
detail = []
for c, d in itertools.combinations(range(COLS), 2):
    labs = [((k - 1) * (d - c)) % Q for k in range(1, T + 1)]
    hand = hand and len(set(labs)) == T
    detail.append("(%d,%d):%s" % (c, d, "".join(str(x) for x in labs)))
check("five-blocks-use-four-distinct-cross-labels-on-every-column-pair", hand,
      "on columns (c,c') they use 0, delta, 2delta, 3delta with delta = c'-c != 0 mod 13, distinct "
      "because 13 is prime -- " + " ".join(detail))

leave = set(orbits) - set(l for five, _ in base for l in
                          (orbit_label(*p) for p in itertools.combinations(five, 2)))
check("the-2-blocks-must-cover-120-leave-orbits",
      len(leave) == 120 and len([l for l in leave if l[0] == "S"]) == 30
      and len([l for l in leave if l[0] == "C"]) == 90 and 120 * Q == w * 30,
      "160 - 40 = 120 orbits remain after the five-blocks (all 30 same-column plus 10 x 9 = 90 "
      "cross-column), and 120 x 13 = %d = 52 x 30 two-blocks" % (120 * Q))

# ---------------------------------------------------------------------------
# 4. the whole development, checked DIRECTLY
# ---------------------------------------------------------------------------
print("")
print("=== the 52 developed parallel classes, checked directly and not via the orbit argument ===")
classes = develop(base)
check("development-gives-52-parallel-classes", len(classes) == w == 52,
      "4 base classes x 13 powers of sigma = %d classes" % len(classes))
check("development-gives-1612-blocks", sum(len(cl) for cl in classes) == 1612,
      "52 x 31 = %d blocks: 52 of size 5 and 1560 of size 2" % sum(len(cl) for cl in classes))
check("every-developed-class-is-a-partition-of-type-5-1-2-30",
      all(is_type_5_1_2_30(cl) for cl in classes),
      "all 52 classes are resolutions of the 65 varieties with the same block-size multiset")

mult, prof = lambda_profile(classes)
covered = sum(1 for p in itertools.combinations(range(N), 2) if mult[p] >= 1)
check("all-2080-pairs-are-covered", covered == n_pairs,
      "%d of %d distinct pairs appear in some block" % (covered, n_pairs))
check("the-multiplicity-multiset-is-exactly-one-for-every-pair", dict(prof) == {1: n_pairs},
      "pair-multiplicity multiset = {1: %d}, so the design is pairwise balanced with lambda = 1"
      % n_pairs)

deg = Counter()
for cl in classes:
    for blk in cl:
        if len(blk) == M:
            for u in blk:
                deg[u] += 1
check("t-equals-4-at-every-variety", set(deg.values()) == {T} and len(deg) == N,
      "every one of the 65 varieties lies in exactly %d of the 52 five-blocks" % T)

fives = [blk for cl in classes for blk in cl if len(blk) == M]
check("the-52-five-blocks-are-distinct", len(set(fives)) == 52,
      "52 distinct five-blocks, one per parallel class")
check("five-blocks-pairwise-meet-in-at-most-one-variety",
      all(len(set(a) & set(b)) <= 1 for a, b in itertools.combinations(fives, 2)),
      "no two of the 52 five-blocks share two varieties, as lambda = 1 requires")

for k in range(T):
    fam = [fives[k * Q + s] for s in range(Q)]
    flat = sorted(u for blk in fam for u in blk)
    check("family-%d-of-thirteen-five-blocks-is-itself-a-parallel-class" % (k + 1),
          flat == list(range(N)),
          "variety 13c+x lies in sigma^s(P_%d) iff s = x - %d c mod 13, one solution in Z_13, so the "
          "13 five-blocks partition the 65 varieties" % (k + 1, k))

# ---------------------------------------------------------------------------
# 5. anti-controls: the same lambda test must REJECT damaged copies
# ---------------------------------------------------------------------------
print("")
print("=== anti-controls, so that the PASSes above are not vacuous ===")
bad = [(f, list(p)) for f, p in base]
(a1, b1), (a2, b2) = bad[0][1][0], bad[0][1][1]
bad[0][1][0], bad[0][1][1] = tuple(sorted((a1, b2))), tuple(sorted((a2, b1)))
_, prof_bad = lambda_profile(develop(bad))
check("anti-control-a-swap-of-two-varieties-is-rejected", dict(prof_bad) != {1: n_pairs},
      "exchanging %s with %s in P1 gives multiplicity multiset %s -- the test says INVALID"
      % ("%d_%d" % coords(b1), "%d_%d" % coords(b2), dict(sorted(prof_bad.items()))))

dup = [(f, list(p)) for f, p in base]
dup[1] = dup[0]
_, prof_dup = lambda_profile(develop(dup))
check("anti-control-a-repeated-base-class-is-rejected", dict(prof_dup) != {1: n_pairs},
      "setting P2 := P1 gives multiplicity multiset %s: only 120 of the 160 orbits are supplied, so "
      "520 pairs are doubled and 520 are missed" % dict(sorted(prof_dup.items())))

# ---------------------------------------------------------------------------
# 6. the 2004 near miss, distinguished by integer arithmetic alone
# ---------------------------------------------------------------------------
print("")
print("=== the near miss: Danziger-Stevens, EJC 11 (2004) #R23, Theorem 3.5 at (q,k,m) = (13,5,3) ===")
q, kk, mm = 13, 5, 3
check("r23-theorem-3-5-is-legal-at-13-5-3",
      2 <= kk <= q and 0 <= mm and mm * (kk - 1) <= q,
      "2 <= k = %d <= q = %d and 0 <= m = %d <= q/(k-1) = 13/4, so the theorem applies"
      % (kk, q, mm))
check("r23-lands-on-our-exact-order-and-partition",
      kk * q == N and mm * kk * (kk - 1) // 2 == 30 and q - mm * (kk - 1) == 1,
      "type q^k = 13^5 is %d points with partition 2^%d %d^%d -- our order and our partition"
      % (kk * q, mm * kk * (kk - 1) // 2, kk, q - mm * (kk - 1)))
within = kk * (q * (q - 1) // 2)
across = (kk * (kk - 1) // 2) * q * q
check("r23-object-is-group-divisible-so-390-pairs-are-covered-zero-times",
      within == 390 and across == 1690 and within + across == n_pairs,
      "%d within-group pairs are covered 0 times and only %d of the %d pairs are covered at all; a "
      "CURD must cover all %d" % (within, across, n_pairs, n_pairs))
lam = q - mm * (kk - 2)
res, rres = divmod(lam * across, per_class)
check("r23-index-is-4-not-1", lam == 4 and rres == 0 and res == 169,
      "lambda = q - m(k-2) = %d, giving %d resolution classes, not lambda = 1 with 52" % (lam, res))
witnesses = [(qq, mmm, kkk)
             for qq in range(1, 201) for kkk in range(2, min(qq, 200) + 1)
             for mmm in range(0, 201)
             if qq - mmm * (kkk - 2) == 1 and qq - mmm * (kkk - 1) == 1
             and not (mmm == 0 and qq == 1)]
check("r23-cannot-be-specialised-to-lambda-1-with-one-k-block-per-class", witnesses == [],
      "lambda = 1 needs q-m(k-2) = 1 while one k-block per class needs q-m(k-1) = 1; subtracting gives "
      "m = 0 hence q = 1, and an exhaustive sweep of 1 <= q <= 200, 0 <= m <= 200, 2 <= k <= q finds no "
      "other solution")
check("65-is-not-a-projective-plane-order-so-corollary-3-8-misses-it",
      all(x * x + x + 1 != N for x in range(1, 101)),
      "R23's only lambda = 1 CURD family lands on q^2+q+1 points; 7^2+7+1 = 57 and 8^2+8+1 = 73 "
      "straddle 65, so 65 is not of that form for any integer q")
p2, p3 = 30, 0
check("the-published-two-block-size-bound-does-not-exclude-65",
      p3 == 0 and p2 > 1 and N <= p2 * p2,
      "the reviewed exclusion bounds are p_2 > 1 => v <= p_2^2 and p_3 > 0 => v <= (3p_3)^2; here "
      "p_2 = %d and p_3 = %d, so the only applicable one reads v <= %d and 65 <= %d"
      % (p2, p3, p2 * p2, p2 * p2))

# ---------------------------------------------------------------------------
# 7. the bookkeeping of the claim
# ---------------------------------------------------------------------------
print("")
print("=== the bookkeeping of the claim ===")
n28 = 28
w28, r28 = divmod(n28 * (n28 - 1), 4 * 4 - 2 * 4 + n28)
t28, s28 = divmod(n28 - w28 - 1, 4 - 2)
check("the-source-s-other-open-cell-has-different-parameters",
      r28 == 0 and s28 == 0 and (w28, t28) == (21, 3) and (w28, t28) != (w, t),
      "the 4^1 2^12 cell at n = 28 has m = 4, w = %d, t = %d against our m = 5, w = %d, t = %d; this "
      "paper settles the n = 65 cell alone" % (w28, t28, w, t))
observed = (len(set(u for cl in classes for blk in cl for u in blk)),
            max(len(blk) for cl in classes for blk in cl),
            len(classes),
            max(deg.values()),
            max(prof))
check("the-parameters-read-off-the-object-are-the-claimed-ones",
      observed == (65, 5, 52, 4, 1),
      "measured from the developed blocks rather than from any formula: (n, m, w, t, lambda) = %s, which "
      "is the cell the paper claims and no other" % (observed,))

print("")
print("object sha256 (of the OBJECT string in this file, for comparison with the printed page): %s"
      % hashlib.sha256(OBJECT.encode("utf-8")).hexdigest())
print("")
for line in NOT_RE_RUN:
    print(line)
print("")
if _fails:
    print("VERDICT: %d CHECKS FAILED" % _fails)
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _passes)
sys.exit(0)

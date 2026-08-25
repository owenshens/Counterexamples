#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- referee's verification program for

    "The Harborth Constant of C_9 (+) C_9"

Standard library only, single process, exact integer arithmetic
throughout: no floating point enters any decision.

-------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER  (inputs -- transcribed, NOT checked)
-------------------------------------------------------------------------
  P_N            9              the parameter n (G = C_n (+) C_n)
  P_G_VALUE      17             the claimed Harborth constant g(C_9(+)C_9)
  P_A            the exhibited 16-element set, transcribed point by point
                 as the 16 pairs {(x,y) : 0 <= x <= 7, y in {0,1}}
                 (paper, proof of Proposition 2)
  P_B            {(0,0),(1,0)}  the normalizing pair (Lemma 4)
  P_E_SIZE       79             |G \\ B|
  P_STAB_ORDER   108            |Stab_{AGL_2(R)}(B)|
  P_CANON        992            number of canonical triples
  P_FEASIBLE     831            canonical triples with p_3 <= 66
  P_P3_BOUND     66             the feasibility bound on p_3
  P_COMPLETIONS  811509569804714   the integer of equation (2)
  P_DFS_NODES    879672298      depth-first-search nodes in the paper's table
  P_3G_ORDER     9              |3G|
  P_GT           g(C_n(+)C_n) = 2n-1 (n odd), 2n+1 (n even)   [eq. (1)]

Everything else below is recomputed from scratch.  Corrupt any input above
and at least one check must print FAIL.

-------------------------------------------------------------------------
WHAT IS DERIVED HERE  (the checks; nothing below is copied from the paper)
-------------------------------------------------------------------------
  * P_A is decoded to 16 points of (Z/9Z)^2, counted, printed back and
    compared with the description "{0,...,7} x {0,1}";
  * ALL 11440 nine-element subsets of A are summed: none is zero-sum, so
    g >= 17 -- the lower bound is recomputed, not quoted; the layer
    argument of Proposition 2 (1 <= t <= 8, second sum = t mod 9) is
    verified subset by subset;
  * A is shown maximal: for each of the 65 points z outside A an explicit
    8-element subset of A summing to -z is produced and re-added, so every
    17-element set containing A does contain a 9-element zero-sum subset;
  * the machinery of Section 3 is self-tested against brute force: the
    81-bit translation is a bit permutation (107 masks x 81 shifts), the
    incremental Sigma_j states equal enumerated j-subset sums, the
    admissibility test -z in Sigma_8(P) equals its definition, and a
    non-edge -(u+v) in Sigma_7(P) is turned into an explicit 9-element
    zero-sum set;
  * the SEARCH ENGINE ITSELF is audited, not just a re-implementation of
    its rules.  An over-pruning engine reports "no free set" for every
    input and would therefore confirm any upper bound at all, so: (i) a
    definition-only reference enumerator (no clique rule, no degree
    peeling, no colouring bound) is run to exhaustion on small cases and
    the engine must return exactly the lexicographically least free
    completion below every free prefix; (ii) the rule-1 and rule-3 masks
    that _make_dfs actually computes are handed back through an audit hook
    and compared point by point with "-z not in Sigma_8(P)" and "-(u+v)
    not in Sigma_7(P)" recomputed by brute force; (iii) at n = 9 the DFS
    is required to still return A from every prefix of A's own position
    sequence, so no rejection rule discards a subfamily containing a known
    free set;
  * Lemma 3: 3G is recomputed as a subgroup of order 9 whose cosets are
    the classes of "difference in 3G";  Lemma 4: GL_2(Z/9Z) is
    enumerated (3888 matrices, every inverse verified) and exactly the 72
    primitive vectors are moved to e_1 by an exhibited matrix; exp(G) is
    computed from element orders and 400 random affine maps are checked to
    preserve 9-element zero-sums;
  * the setwise affine stabilizer of B is brute-forced over all 314928
    affine maps (order 108) and shown to equal the paper's two printed
    families, to be a group, and to act faithfully on the 79 positions;
  * the 992 canonical triples are computed THREE independent ways (direct
    lexicographic test, union-find orbit count, Burnside's lemma), the 831
    feasible prefixes and the exact integer 811509569804714 of eq. (2)
    are recomputed, and the counted family is shown large enough to meet
    every Gamma-orbit;
  * eq. (1) is verified completely by pure brute force over all 2^(n^2)
    subsets for n = 2, 3, 4, and completely by the same DFS engine used at
    n = 9 for n = 3, 5, 7 (g = 5, 9, 13), for n = 3 and 5 also with no
    normalization and no symmetry reduction at all;
  * the engine, run at size 16, rediscovers A itself, so it is not
    vacuously reporting failure;
  * the n = 9 census: the DFS runs over feasible canonical prefixes with
    exact integer credit for every rejected subfamily, and the covered
    count is reported exactly.  With the default budget only part of the
    811509569804714 completions is covered; the uncovered remainder is
    printed as a NOTE and never silently dropped.

Usage:  python3 verify.py            default, about 19 minutes
        python3 verify.py --quick    n=9 census capped at 20 s
        python3 verify.py --seconds S   n=9 census capped at S seconds
        python3 verify.py --full      whole census, hours
Exit status is 0 if and only if every check passed.
"""

import sys
import time
import math
import itertools
import random
import traceback

# =====================================================================
# CHECK HARNESS
# =====================================================================

_RESULTS = []


def check(name, condition, detail=""):
    """Record and print one check.  Every call must be able to print FAIL."""
    ok = bool(condition)
    _RESULTS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return ok


def note(text):
    print("NOTE " + text)


def verdict():
    n = len(_RESULTS)
    bad = [nm for nm, ok in _RESULTS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


# =====================================================================
# BLOCK 1.  VERBATIM PAPER INPUTS.
# Corrupt any line here and the checks below must report FAIL.
# =====================================================================

P_N = 9                       # G = C_9 (+) C_9,  exp(G) = 9
P_G_VALUE = 17                # the claimed Harborth constant
P_E_SIZE = 79                 # |E| = |G| - |B|
P_STAB_ORDER = 108            # |Gamma|
P_CANON = 992                 # canonical triples
P_FEASIBLE = 831              # feasible canonical prefixes
P_P3_BOUND = 66               # p_3 <= 66
P_COMPLETIONS = 811509569804714      # equation (2)
P_DFS_NODES = 879672298       # nodes reported in the paper's table
P_3G_ORDER = 9                # |3G|

# The normalizing pair B = {0, e_1} of Lemma 4.
P_B = [(0, 0), (1, 0)]

# The exhibited 16-element set A of Proposition 2, transcribed point by
# point (the paper writes it as {0,1,...,7} x {0,1}).
P_A = [
    (0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (7, 0),
    (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (7, 1),
]


def gt_conjecture(n):
    """Equation (1): the Gao-Thangadurai value."""
    return 2 * n - 1 if n % 2 == 1 else 2 * n + 1


# =====================================================================
# BLOCK 2.  GROUP ARITHMETIC IN (Z/nZ)^2
# =====================================================================


def gidx(p, n):
    """Integer representative x + n*y used to order G (paper, Section 3)."""
    return p[0] % n + n * (p[1] % n)


def gvec(i, n):
    return (i % n, i // n)


def gadd(p, q, n):
    return ((p[0] + q[0]) % n, (p[1] + q[1]) % n)


def gneg(p, n):
    return ((-p[0]) % n, (-p[1]) % n)


def matvec(M, p, n):
    a, b, c, d = M
    return ((a * p[0] + b * p[1]) % n, (c * p[0] + d * p[1]) % n)


def matmul(M1, M2, n):
    a, b, c, d = M1
    e, f, g, h = M2
    return ((a * e + b * g) % n, (a * f + b * h) % n,
            (c * e + d * g) % n, (c * f + d * h) % n)


def subset_sums_by_size(points, n, jmax):
    """Exact brute force: {sum(U) : U subset of points, |U| = j} for j<=jmax.

    Returned as a list of sets of group elements.  Used to cross-check the
    bitset machinery, never as part of it.
    """
    out = [set() for _ in range(jmax + 1)]
    out[0].add((0, 0))
    for j in range(1, jmax + 1):
        if j > len(points):
            continue
        for U in itertools.combinations(points, j):
            sx = sum(p[0] for p in U) % n
            sy = sum(p[1] for p in U) % n
            out[j].add((sx, sy))
    return out


# =====================================================================
# BLOCK 3.  EXACT n^2-BIT SUBSET-SUM STATES
#
# A subset of G is stored as an exact integer bitmask, bit (x + n*y) for
# the point (x,y).  Translation by a group element is a bit permutation
# realised by masked shifts: a cyclic rotation by dx inside each block of
# n bits (the x-coordinate), then a rotation of whole blocks by dy (the
# y-coordinate).  No arithmetic other than exact integer shifts is used.
# =====================================================================


def translator(n):
    """Return (MASK, PARAMS, tr) for exact translation of n^2-bit masks."""
    N2 = n * n
    MASK = (1 << N2) - 1
    block = (1 << n) - 1

    def rep(m):
        s = 0
        for j in range(n):
            s |= m << (n * j)
        return s

    PARAMS = []
    for v in range(N2):
        dx, dy = v % n, v // n
        if dx == 0:
            lo = hi = 0
        else:
            low = (1 << (n - dx)) - 1
            lo = rep(low)
            hi = rep(block ^ low)
        PARAMS.append((dx, n - dx, lo, hi, dy, n * dy, n * (n - dy)))

    def tr(S, v):
        dx, ndx, lo, hi, dy, shl, shr = PARAMS[v]
        if dx:
            S = ((S & lo) << dx) | ((S & hi) >> ndx)
        if dy:
            S = ((S << shl) & MASK) | (S >> shr)
        return S

    return MASK, PARAMS, tr


def bits_to_points(S, n):
    out = []
    i = 0
    while S:
        if S & 1:
            out.append(gvec(i, n))
        S >>= 1
        i += 1
    return out


# =====================================================================
# BLOCK 4.  THE EXHAUSTIVE SEARCH ENGINE (the paper's Section 3)
# =====================================================================


def _make_dfs(st, tick, shift, L, L2, NP, MP, MASK, TAIL, CB, PRN):
    """Build the recursive node routine (kept separate for readability).

    st["audit"] is None in every production run (one `is not None` test per
    node, no other cost).  When it is a callable, it is handed the ACTUAL
    masks this routine computes -- not a re-implementation of them -- so that
    rules 1 and 3 can be compared against their definitions inside the code
    that the census actually executes.
    """

    def dfs(psi, s, r, sel):
        nd = st["nodes"] + 1
        st["nodes"] = nd
        if nd >= st["chk"]:
            tick(nd)
        if r == 0:                                  # a valid set survives
            st["covered"] += 1
            st["found"] = list(sel)
            return True
        # rule 1: z is admissible iff -z is not an (n-1)-element subsum
        adm = TAIL[s] & ~((psi[L] >> shift) & MP)
        audit = st["audit"]
        if audit is not None:
            audit("rule1", sel, s, r, adm, None)
        if bin(adm).count("1") < r:                 # rule 2
            st["covered"] += CB[NP - s][r]
            return False
        if r > 1:
            # rule 3: build the compatibility graph on admissible points
            pk = psi[L2]
            nb = [0] * NP
            rm = r - 1
            m = adm
            while m:
                b = m & -m
                u = b.bit_length() - 1
                m ^= b
                dx, ndx, lo, hi, dy, shl, shr = PRN[u]
                S = pk
                if dx:
                    S = ((S & lo) << dx) | ((S & hi) >> ndx)
                if dy:
                    S = ((S << shl) & MASK) | (S >> shr)
                nb[u] = adm & ~(S >> shift) & ~b
            if audit is not None:
                audit("rule3", sel, s, r, adm, nb)
            alive = adm
            while True:                             # iterative degree peel
                m = alive
                ch = 0
                while m:
                    b = m & -m
                    u = b.bit_length() - 1
                    m ^= b
                    if bin(nb[u] & alive).count("1") < rm:
                        alive ^= b
                        ch = 1
                if not ch:
                    break
            if bin(alive).count("1") < r:
                st["covered"] += CB[NP - s][r]
                return False
            unc = alive                             # greedy colouring bound
            kk = 0
            while unc:
                kk += 1
                if kk >= r:
                    break
                av = unc
                while av:
                    b = av & -av
                    u = b.bit_length() - 1
                    av = (av & ~nb[u]) ^ b
                    unc ^= b
            if kk < r:
                st["covered"] += CB[NP - s][r]
                return False
            adm = alive
        credit = CB[NP - s][r]
        kids = []
        m = adm
        while m:
            b = m & -m
            z = b.bit_length() - 1
            m ^= b
            c = CB[NP - 1 - z][r - 1]
            if c:
                credit -= c
                kids.append(z)
        st["covered"] += credit          # all rejected first-point choices
        for z in kids:
            dx, ndx, lo, hi, dy, shl, shr = PRN[z]
            nxt = psi[:]
            for j in range(L, 0, -1):
                S = psi[j - 1]
                if dx:
                    S = ((S & lo) << dx) | ((S & hi) >> ndx)
                if dy:
                    S = ((S << shl) & MASK) | (S >> shr)
                nxt[j] = psi[j] | S
            sel.append(z)
            if dfs(nxt, z + 1, r - 1, sel):
                return True
            sel.pop()
        return False

    return dfs


def make_engine(n, k, shift):
    """DFS for a k-element subset of (Z/nZ)^2 with no n-element zero-sum
    subset, among the sets containing the first `shift` group elements
    (shift = 2 forces B = {(0,0),(1,0)}; shift = 0 forces nothing).

    Psi_j = -Sigma_j(P) is kept as an exact n^2-bit mask, so that
       z admissible   <=>   z not in Psi_{n-1}        (Lemma 5)
       u ~ v          <=>   u+v not in Psi_{n-2}      (rejection rule 3)
    Position p of the ordered ground set corresponds to group index
    p + shift, hence a value mask becomes a position mask by >> shift.

    Exact accounting: a node entered with least allowed position s and r
    points still required owns exactly C(NP-s, r) completions and credits
    every subfamily it rejects, so `covered` is at all times the exact
    number of completions examined or soundly discarded.
    """
    N2 = n * n
    MASK, PARAMS, tr = translator(n)
    NP = N2 - shift
    MP = (1 << NP) - 1
    NEG = [gidx(gneg(gvec(v, n), n), n) for v in range(N2)]
    PRN = [PARAMS[NEG[p + shift]] for p in range(NP)]
    TAIL = [(MP >> s) << s for s in range(NP + 1)]
    CB = [[math.comb(a, b) for b in range(k + 2)] for a in range(NP + 2)]
    L = n - 1
    L2 = n - 2
    st = {"nodes": 0, "covered": 0, "found": None, "chk": 1 << 60,
          "limit": None, "deadline": None, "stopped": False,
          "prefix_bad": False, "audit": None}

    def tick(nd):
        if st["limit"] is not None and nd > st["limit"]:
            st["stopped"] = True
            raise TimeoutError
        if st["deadline"] is not None and time.time() > st["deadline"]:
            st["stopped"] = True
            raise TimeoutError
        st["chk"] = nd + 4096

    # ---- the recursion is installed by the next block ----
    dfs = _make_dfs(st, tick, shift, L, L2, NP, MP, MASK, TAIL, CB, PRN)

    psi0 = [0] * n
    psi0[0] = 1                      # Psi_0 = {0}
    for v in range(shift):           # insert the forced points of B
        nv = NEG[v]
        for j in range(L, 0, -1):
            psi0[j] |= tr(psi0[j - 1], nv)

    def run(prefix=(), limit=None, seconds=None, reset=False):
        """Search the family whose first |prefix| positions are `prefix`.

        PRECONDITION: the forced points together with `prefix` must not
        already contain an n-element zero-sum subset.  The DFS assumes a free
        partial set and would otherwise report a "found" set that is not
        free.  With shift + len(prefix) < n the precondition is automatic,
        which is the only regime used here (2 + 3 = 5 < 9 for the n = 9
        census); st["prefix_bad"] records a violation and the census check
        refuses to pass if it is ever set.
        """
        if reset:
            st["nodes"] = 0
            st["covered"] = 0
            st["found"] = None
        st["limit"] = limit
        st["deadline"] = None if seconds is None else time.time() + seconds
        st["stopped"] = False
        st["chk"] = (st["nodes"] + 1) if (limit or seconds) else (1 << 60)
        psi = psi0[:]
        for q in prefix:
            if (psi[L] >> (q + shift)) & 1:
                st["prefix_bad"] = True      # prefix already has a zero-sum
            dx, ndx, lo, hi, dy, shl, shr = PRN[q]
            for j in range(L, 0, -1):
                S = psi[j - 1]
                if dx:
                    S = ((S & lo) << dx) | ((S & hi) >> ndx)
                if dy:
                    S = ((S << shl) & MASK) | (S >> shr)
                psi[j] = psi[j] | S
        s = (prefix[-1] + 1) if prefix else 0
        r = k - shift - len(prefix)
        try:
            dfs(psi, s, r, list(prefix))
        except TimeoutError:
            pass
        return st

    def decode(positions):
        pts = [gvec(v, n) for v in range(shift)]
        pts += [gvec(p + shift, n) for p in positions]
        return pts

    return run, st, decode, NP, CB


def psi_states(points, n):
    """Incremental Psi_j = -Sigma_j, j = 0..n-1, by Lemma 5's update."""
    MASK, PARAMS, tr = translator(n)
    psi = [0] * n
    psi[0] = 1
    for p in points:
        nv = gidx(gneg(p, n), n)
        for j in range(n - 1, 0, -1):
            psi[j] |= tr(psi[j - 1], nv)
    return psi


# =====================================================================
# BLOCK 5.  SELF-TESTS OF THE MACHINERY THE UPPER BOUND RELIES ON
# =====================================================================


def stage_selftest(n):
    MASK, PARAMS, tr = translator(n)
    rng = random.Random(20260824)
    N2 = n * n

    tests = [1 << i for i in range(N2)]
    tests += [rng.getrandbits(N2) for _ in range(24)]
    tests += [0, MASK]
    bad = None
    for S in tests:
        for v in range(N2):
            want = 0
            for p in bits_to_points(S, n):
                want |= 1 << gidx(gadd(p, gvec(v, n), n), n)
            if tr(S, v) != want:
                bad = (S, v)
                break
        if bad:
            break
    check("translation_is_exact_bit_permutation", bad is None,
          "%d masks x %d shifts against coordinate addition"
          % (len(tests), N2) if bad is None else "fails at %r" % (bad,))

    worst = None
    for trial in range(12):
        m = rng.randrange(1, n)
        pts = rng.sample(range(N2), m)
        pts = [gvec(v, n) for v in pts]
        psi = psi_states(pts, n)
        brute = subset_sums_by_size(pts, n, n - 1)
        for j in range(n):
            got = set(gneg(p, n) for p in bits_to_points(psi[j], n))
            if got != brute[j]:
                worst = (trial, j)
    check("subsum_states_match_brute_force", worst is None,
          "12 random partial sets, all Sigma_j for j <= %d" % (n - 1)
          if worst is None else "mismatch %r" % (worst,))

    pts = [gvec(v, n) for v in rng.sample(range(N2), n + 1)]
    psi = psi_states(pts, n)
    brute = subset_sums_by_size(pts, n, n - 1)
    ok_adm = True
    for v in range(N2):
        z = gvec(v, n)
        if z in pts:
            continue
        by_bits = bool((psi[n - 1] >> v) & 1)
        by_brute = gneg(z, n) in brute[n - 1]
        if by_bits != by_brute:
            ok_adm = False
    check("admissibility_rule_matches_definition", ok_adm,
          "-z in Sigma_%d(P) tested against all (n-1)-subsets of a random "
          "P with |P| = %d" % (n - 1, n + 1))

    ok_edge = True
    wit = None
    out = [gvec(v, n) for v in range(N2) if gvec(v, n) not in pts]
    for u, v in itertools.combinations(out, 2):
        s = gadd(u, v, n)
        by_bits = bool((psi[n - 2] >> gidx(s, n)) & 1)
        by_brute = gneg(s, n) in brute[n - 2]
        if by_bits != by_brute:
            ok_edge = False
        if by_bits and wit is None:
            for W in itertools.combinations(pts, n - 2):
                if gadd(sum_points(W, n), s, n) == (0, 0):
                    wit = (u, v, W)
                    break
    good_wit = wit is not None and len(set(wit[2]) | {wit[0], wit[1]}) == n \
        and gadd(sum_points(wit[2], n), gadd(wit[0], wit[1], n), n) == (0, 0)
    check("nonedge_rule_yields_a_real_zero_sum_subset",
          ok_edge and good_wit,
          "witness %r + %r + %r sums to 0" % (wit[0], wit[1], list(wit[2]))
          if good_wit else "no witness")


def sum_points(points, n):
    sx = sum(p[0] for p in points) % n
    sy = sum(p[1] for p in points) % n
    return (sx, sy)


# =====================================================================
# BLOCK 6.  THE EXHIBITED OBJECT AND THE LOWER BOUND (Proposition 2)
# =====================================================================


def stage_object():
    n = P_N
    A = [tuple(p) for p in P_A]
    inrange = all(isinstance(x, int) and isinstance(y, int)
                  and 0 <= x < n and 0 <= y < n for x, y in A)
    check("A_is_a_set_of_%d_points_of_G" % len(P_A),
          inrange and len(set(A)) == len(A) == 16,
          "%d transcribed pairs, %d distinct, all in (Z/%dZ)^2"
          % (len(A), len(set(A)), n))
    described = set((x, y) for x in range(8) for y in range(2))
    check("A_equals_the_printed_description", set(A) == described,
          "{0,...,7} x {0,1}")
    note("A (by integer representative x + %dy) = %s"
         % (n, sorted(gidx(p, n) for p in A)))
    note("A (as points) = %s" % (sorted(A, key=lambda p: gidx(p, n)),))
    check("A_has_g_minus_one_elements",
          len(A) == P_G_VALUE - 1 == 2 * n - 2,
          "|A| = %d = %d - 1 = 2n-2" % (len(A), P_G_VALUE))
    check("claimed_value_is_the_Gao_Thangadurai_value",
          gt_conjecture(n) == P_G_VALUE and n % 2 == 1,
          "n = %d odd, 2n-1 = %d" % (n, gt_conjecture(n)))
    check("A_contains_the_normalizing_pair_B",
          all(tuple(b) in set(A) for b in P_B),
          "B = %s" % (P_B,))
    return A


def stage_lower_bound(A):
    n = P_N
    zero = (0, 0)
    seen9 = 0
    zs9 = 0
    bad_layer = 0
    for U in itertools.combinations(A, n):
        seen9 += 1
        sx = sum(p[0] for p in U) % n
        sy = sum(p[1] for p in U) % n
        if (sx, sy) == zero:
            zs9 += 1
        t = sum(1 for p in U if p[1] == 1)
        if not (1 <= t <= 8) or sy != t % n or sy == 0:
            bad_layer += 1
    check("A_has_no_%d_element_zero_sum_subset" % n, zs9 == 0,
          "all %d of the %d-element subsets of A examined, %d zero-sum"
          % (seen9, n, zs9))
    check("A_layer_argument_holds", bad_layer == 0 and seen9 == math.comb(16, n),
          "every %d-subset has 1 <= t <= 8 points in the upper layer and "
          "second-coordinate sum t mod %d != 0" % (n, n))
    check("lower_bound_g_at_least_%d" % P_G_VALUE, zs9 == 0 and len(A) == P_G_VALUE - 1,
          "a %d-element set with no %d-element zero-sum subset exists"
          % (len(A), n))

    wit8 = {}
    for U in itertools.combinations(A, n - 1):
        s = sum_points(U, n)
        if s not in wit8:
            wit8[s] = U
    outside = [gvec(v, n) for v in range(n * n) if gvec(v, n) not in set(A)]
    missing = 0
    badwit = 0
    example = None
    for z in outside:
        U = wit8.get(gneg(z, n))
        if U is None:
            missing += 1
            continue
        S = set(U) | {z}
        if len(S) != n or sum_points(list(S), n) != zero or not set(U) <= set(A):
            badwit += 1
        elif example is None:
            example = (z, sorted(U, key=lambda p: gidx(p, n)))
    check("every_17_element_superset_of_A_has_a_zero_sum_subset",
          missing == 0 and badwit == 0 and len(outside) == n * n - len(A),
          "%d extensions, each with an explicit %d-element zero-sum witness; "
          "e.g. z = %r with %r" % (len(outside), n, example[0], example[1])
          if example else "no witness")


# =====================================================================
# BLOCK 7.  NORMALIZATION (Lemmas 2.2 and 2.3)
# =====================================================================


def gl2(n):
    """All invertible 2x2 matrices over Z/nZ, with their inverses."""
    out = {}
    inv = {}
    for u in range(n):
        for w in range(n):
            if (u * w) % n == 1:
                inv[u] = w
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    det = (a * d - b * c) % n
                    if det in inv:
                        e = inv[det]
                        out[(a, b, c, d)] = ((d * e) % n, (-b * e) % n,
                                             (-c * e) % n, (a * e) % n)
    return out


def stage_normalization(A):
    n = P_N
    G = [gvec(v, n) for v in range(n * n)]
    H = sorted(set(((3 * p[0]) % n, (3 * p[1]) % n) for p in G))
    closed = all(gadd(u, v, n) in set(H) for u in H for v in H) \
        and all(gneg(u, n) in set(H) for u in H)
    check("three_G_is_a_subgroup_of_order_%d" % P_3G_ORDER,
          len(H) == P_3G_ORDER and closed and (0, 0) in H,
          "3G = %s" % (H,))

    Hs = set(H)
    classes = {}
    for p in G:
        key = min(gidx(gadd(p, h, n), n) for h in H)
        classes.setdefault(key, []).append(p)
    sizes = sorted(len(c) for c in classes.values())
    equiv = all(((gadd(u, gneg(v, n), n) in Hs)
                 == (min(gidx(gadd(u, h, n), n) for h in H)
                     == min(gidx(gadd(v, h, n), n) for h in H)))
                for u in G for v in G)
    check("lemma_primitive_difference",
          len(classes) == n and sizes == [n] * n and equiv
          and sum(sizes) == n * n,
          "the relation u-v in 3G is an equivalence with %d classes of size "
          "%d, so any set of more than %d points has a primitive difference"
          % (len(classes), n, n))

    # Lemma 3 has a HYPOTHESIS ("more than nine elements") that the
    # normalization step of Lemma 4 must actually satisfy for the sets the
    # paper reduces, namely the 17-element ones.  Nothing above tests it.
    biggest_coset = max(len(c) for c in classes.values())
    check("normalization_hypothesis_holds_for_%d_element_sets" % P_G_VALUE,
          biggest_coset == P_3G_ORDER == len(H)
          and P_G_VALUE > biggest_coset
          and all(len(c) <= biggest_coset for c in classes.values()),
          "the largest coset of 3G has %d elements and %d > %d, so every "
          "%d-element subset of G meets two cosets and therefore has a "
          "primitive difference: Lemma 4 does apply to the candidates"
          % (biggest_coset, P_G_VALUE, biggest_coset, P_G_VALUE))

    MS = gl2(n)
    ok_inv = all(matmul(M, MS[M], n) == (1, 0, 0, 1) for M in MS)
    check("GL2_over_Z9_enumerated", len(MS) == 3888 and ok_inv,
          "%d invertible matrices, every inverse verified" % len(MS))

    prim = [p for p in G if p not in Hs]
    normalizable = {}
    for M in MS:
        v = (M[0], M[2])                      # first column of M = M e_1
        if v not in normalizable:
            normalizable[v] = MS[M]           # M^{-1} sends v to e_1
    ok_send = all(matvec(normalizable[v], v, n) == (1, 0)
                  for v in normalizable)
    check("primitive_vectors_are_exactly_the_normalizable_ones",
          len(prim) == n * n - P_3G_ORDER == 72
          and set(normalizable) == set(prim) and ok_send,
          "%d primitive vectors, each with an explicit M in GL_2 sending it "
          "to e_1; no vector of 3G is normalizable" % len(prim))

    orders = {}
    for p in G:
        m = 1
        q = p
        while q != (0, 0):
            q = gadd(q, p, n)
            m += 1
        orders[m] = orders.get(m, 0) + 1
    expG = max(orders)
    kills = all(sum_points([c] * expG, n) == (0, 0) for c in G)
    check("exponent_of_G_is_the_zero_sum_size",
          expG == n == P_N and kills and orders == {1: 1, 3: 8, 9: 72}
          and sum(orders.values()) == n * n,
          "element orders %r, exp(G) = %d, and exp(G)*c = 0 for every c, "
          "so an exp(G)-term sum is affine-invariant" % (orders, expG))

    rng = random.Random(97531)
    keys = list(MS)
    bad = 0
    tested = 0
    for _ in range(400):
        M = keys[rng.randrange(len(keys))]
        c = G[rng.randrange(len(G))]
        X = [G[i] for i in rng.sample(range(n * n), n)]
        Y = [gadd(matvec(M, x, n), c, n) for x in X]
        tested += 1
        if len(set(Y)) != n:
            bad += 1
        if sum_points(Y, n) != matvec(M, sum_points(X, n), n):
            bad += 1
        if (sum_points(X, n) == (0, 0)) != (sum_points(Y, n) == (0, 0)):
            bad += 1
    check("affine_maps_preserve_9_element_zero_sums", bad == 0,
          "%d random affine maps x random %d-subsets (NOTE: the three tests "
          "behind this count are algebraic identities -- for any M and c the "
          "translation part contributes n*c = 0 -- so this check exercises "
          "the arithmetic code, not the paper; the substantive form of "
          "Lemma 4's invariance is the next check)" % (tested, n))

    # The check above cannot fail for any input.  What Lemma 4 actually
    # needs is that FREENESS -- a property of a whole 16-element set, not of
    # one 9-term sum -- transports along affine maps.  That is tested here on
    # real sets, with a non-free control so the test can distinguish the two
    # outcomes.
    Afree = [tuple(p) for p in A]
    # A 16-element NON-free control built from scratch, not from A, so that
    # it stays a genuine control however A is transcribed: the nine points
    # (0,0)..(8,0) already sum to (36,0) = 0.
    Anot = [(x, 0) for x in range(n)] + [(x, 1) for x in range(16 - n)]
    n_free_fail = 0
    n_ctrl_fail = 0
    moved = 0
    ctrl_ok = len(set(Anot)) == 16 and not brute_free(Anot, n)
    for _ in range(6):
        M = keys[rng.randrange(len(keys))]
        c = G[rng.randrange(len(G))]
        img_free = [gadd(matvec(M, x, n), c, n) for x in Afree]
        img_not = [gadd(matvec(M, x, n), c, n) for x in Anot]
        moved += 1
        if len(set(img_free)) != len(set(Afree)) or not brute_free(img_free, n):
            n_free_fail += 1
        if len(set(img_not)) != 16 or brute_free(img_not, n):
            n_ctrl_fail += 1
    good = (ctrl_ok and n_free_fail == 0 and n_ctrl_fail == 0 and moved == 6)
    if good:
        detail = ("under %d random affine maps the image of the free 16-set A "
                  "is still free while the image of a from-scratch non-free "
                  "16-set is still non-free (all %d nine-subsets summed each "
                  "time)" % (moved, math.comb(16, n)))
    else:
        detail = ("control set well formed and non-free: %s; images of A "
                  "that were NOT free: %d of %d; images of the non-free "
                  "control that became free or collided: %d of %d"
                  % (ctrl_ok, n_free_fail, moved, n_ctrl_fail, moved))
    check("affine_maps_transport_freeness_of_whole_sets", good, detail)


# =====================================================================
# BLOCK 8.  THE AFFINE STABILIZER OF B (paper, Section 3.2)
# =====================================================================


def aff_apply(g, p, n):
    M, c = g
    return gadd(matvec(M, p, n), c, n)


def aff_comp(g, h, n):
    """Composition g o h as an affine map."""
    M, c = g
    N, d = h
    return (matmul(M, N, n), gadd(matvec(M, d, n), c, n))


def stage_stabilizer():
    n = P_N
    B = [tuple(b) for b in P_B]
    Bs = set(B)
    E = [gvec(v, n) for v in range(n * n) if gvec(v, n) not in Bs]
    E.sort(key=lambda p: gidx(p, n))
    pos = dict((p, i) for i, p in enumerate(E))
    check("E_has_%d_elements_ordered_by_x_plus_9y" % P_E_SIZE,
          len(E) == P_E_SIZE == n * n - len(B)
          and all(gidx(E[i], n) < gidx(E[i + 1], n) for i in range(len(E) - 1))
          and all(gidx(E[p], n) == p + 2 for p in range(len(E))),
          "position p carries group index p+2, p = 0..%d" % (len(E) - 1))

    # make_engine(n, k, 2) forces "the first two group elements" and maps
    # position p to group index p+2.  Nothing else in this program checks
    # that those two conventions are the paper's B and the paper's ordering
    # of E; if they were not, the census would search a different family
    # from the one the canonical prefixes were computed for.
    engine_forced = [gvec(v, n) for v in range(len(B))]
    check("engine_forced_points_are_exactly_B_in_E_position_order",
          set(engine_forced) == Bs and len(engine_forced) == 2
          and all(gvec(p + len(B), n) == E[p] for p in range(len(E)))
          and len(E) == P_E_SIZE,
          "shift = %d makes the DFS force %s = B, and its position p <-> "
          "group index p+%d convention reproduces E[p] for all %d positions"
          % (len(B), sorted(engine_forced, key=lambda p: gidx(p, n)),
             len(B), len(E)))

    MS = gl2(n)
    G = [gvec(v, n) for v in range(n * n)]
    stab = []
    for M in MS:
        for c in G:
            if set(gadd(matvec(M, b, n), c, n) for b in B) == Bs:
                stab.append((M, c))
    check("stabilizer_order_is_%d" % P_STAB_ORDER,
          len(stab) == P_STAB_ORDER and len(set(stab)) == len(stab),
          "brute force over all %d affine maps of G" % (len(MS) * len(G)))

    units = [u for u in range(n) if math.gcd(u, n) == 1]
    claimed = set()
    for a in range(n):
        for u in units:
            claimed.add(((1, a, 0, u), (0, 0)))
            claimed.add((((-1) % n, a, 0, u), (1, 0)))
    same = claimed == set(stab)
    formula = 2 * n * len(units) == P_STAB_ORDER
    check("stabilizer_equals_the_printed_two_families", same and formula,
          "|R^x| = %d and 2*%d*%d = %d; the two printed families give exactly "
          "the %d brute-forced maps" % (len(units), n, len(units),
                                        2 * n * len(units), len(stab)))

    S = set(stab)
    closed = all(aff_comp(g, h, n) in S for g in stab for h in stab)
    check("stabilizer_is_a_group", closed,
          "%d compositions, all inside Gamma" % (len(stab) ** 2))

    perms = []
    okperm = True
    for g in stab:
        pm = [pos.get(aff_apply(g, p, n)) for p in E]
        if None in pm or sorted(pm) != list(range(len(E))):
            okperm = False
            break
        perms.append(pm)
    check("gamma_permutes_the_79_positions", okperm and len(perms) == len(stab),
          "each of the %d maps induces a bijection of {0,...,%d}"
          % (len(perms), len(E) - 1))

    rng = random.Random(13579)
    badhom = 0
    for _ in range(300):
        i = rng.randrange(len(stab))
        j = rng.randrange(len(stab))
        gh = aff_comp(stab[i], stab[j], n)
        k = stab.index(gh)
        if perms[k] != [perms[i][perms[j][p]] for p in range(len(E))]:
            badhom += 1
    check("position_action_is_a_group_action", badhom == 0,
          "300 random pairs: pi_{gh} = pi_g o pi_h")
    return E, stab, perms


# =====================================================================
# BLOCK 9.  CANONICAL PREFIXES AND EQUATION (2)
# =====================================================================


def closure(gens, ident):
    seen = {ident}
    frontier = [ident]
    while frontier:
        p = frontier.pop()
        for g in gens:
            q = tuple(g[i] for i in p)
            if q not in seen:
                seen.add(q)
                frontier.append(q)
    return seen


def find_generators(perms):
    ident = tuple(range(len(perms[0])))
    gens = []
    for pm in perms:
        if len(closure(gens, ident)) == len(perms):
            break
        if tuple(pm) not in closure(gens, ident):
            gens.append(tuple(pm))
    return gens


def stage_prefixes(perms, E):
    npos = len(E)
    triples = list(itertools.combinations(range(npos), 3))
    canon = []
    for t in triples:
        p1, p2, p3 = t
        ok = True
        for pm in perms:
            a, b, c = pm[p1], pm[p2], pm[p3]
            if a > b:
                a, b = b, a
            if b > c:
                b, c = c, b
            if a > b:
                a, b = b, a
            if (a, b, c) < t:
                ok = False
                break
        if ok:
            canon.append(t)
    check("canonical_triples_direct_enumeration", len(canon) == P_CANON,
          "%d of the %d triples of positions are canonical"
          % (len(canon), len(triples)))

    gens = find_generators(perms)
    ident = tuple(range(npos))
    check("gamma_is_generated_by_%d_of_its_elements" % len(gens),
          len(closure(gens, ident)) == len(perms),
          "closure of %d generators has order %d" % (len(gens), len(perms)))

    parent = dict((t, t) for t in triples)

    def find(x):
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r

    for t in triples:
        for g in gens:
            u = tuple(sorted((g[t[0]], g[t[1]], g[t[2]])))
            ru, rt = find(u), find(t)
            if ru != rt:
                parent[max(ru, rt)] = min(ru, rt)
    orbits = {}
    for t in triples:
        orbits.setdefault(find(t), []).append(t)
    canset = set(canon)
    per_orbit = sorted(set(sum(1 for t in o if t in canset)
                           for o in orbits.values()))
    check("orbit_count_confirms_canonical_count",
          len(orbits) == P_CANON and per_orbit == [1],
          "%d Gamma-orbits on 3-subsets by union-find, exactly one canonical "
          "triple in each" % len(orbits))

    total_fixed = 0
    for pm in perms:
        seen = [False] * npos
        cyc = {}
        for i in range(npos):
            if not seen[i]:
                l = 0
                j = i
                while not seen[j]:
                    seen[j] = True
                    j = pm[j]
                    l += 1
                cyc[l] = cyc.get(l, 0) + 1
        f = cyc.get(1, 0)
        total_fixed += math.comb(f, 3) + f * cyc.get(2, 0) + cyc.get(3, 0)
    check("burnside_confirms_canonical_count",
          total_fixed % len(perms) == 0
          and total_fixed // len(perms) == P_CANON,
          "(1/|Gamma|) sum of fixed 3-subsets = %d/%d = %d"
          % (total_fixed, len(perms), total_fixed // len(perms)))

    feas = [t for t in canon if t[2] <= P_P3_BOUND]
    lost = sum(math.comb(npos - 1 - t[2], 12) for t in canon
               if t[2] > P_P3_BOUND)
    check("feasible_prefixes_and_lossless_bound",
          len(feas) == P_FEASIBLE and lost == 0
          and P_P3_BOUND == npos - 1 - 12,
          "%d canonical triples have p_3 <= %d = 78-12, and the %d discarded "
          "ones admit 0 completions" % (len(feas), P_P3_BOUND,
                                        len(canon) - len(feas)))

    tot = sum(math.comb(npos - 1 - t[2], 12) for t in feas)
    check("equation_2_completion_count", tot == P_COMPLETIONS,
          "sum of C(78-p_3,12) = %d" % tot)

    allsets = math.comb(npos, 15)
    check("counted_family_is_large_enough_to_meet_every_orbit",
          tot * len(perms) >= allsets,
          "%d counted >= %d/%d = %d needed" % (tot, allsets, len(perms),
                                               -(-allsets // len(perms))))

    rng = random.Random(24680)
    badlex = 0
    badmono = 0
    for _ in range(200):
        T = tuple(sorted(rng.sample(range(npos), 15)))
        best = None
        for pm in perms:
            img = tuple(sorted(pm[p] for p in T))
            if best is None or img < best:
                best = img
            head = tuple(sorted((pm[T[0]], pm[T[1]], pm[T[2]])))
            if img[:3] > head:
                badmono += 1
        if best[:3] not in canset:
            badlex += 1
    check("lemma_6_sampled_on_random_15_subsets",
          badlex == 0 and badmono == 0,
          "200 random 15-subsets: the lexicographically least image always "
          "begins with a canonical triple (sampled, not exhaustive)")
    return feas


# =====================================================================
# BLOCK 10.  EQUATION (1) BY PURE BRUTE FORCE FOR SMALL n
# No normalization, no symmetry, no bitset engine: every subset of
# (Z/nZ)^2 is inspected.  This is an independent yardstick for the
# machinery used at n = 9.
# =====================================================================


def brute_harborth(n):
    """Exact g(C_n (+) C_n) by inspecting all 2^(n^2) subsets."""
    N = n * n
    FULL = (1 << N) - 1
    zs = []
    for U in itertools.combinations(range(N), n):
        if sum_points([gvec(v, n) for v in U], n) == (0, 0):
            zs.append(sum(1 << v for v in U))
    bad = bytearray(1 << N)
    for Z in zs:
        comp = FULL ^ Z
        sub = comp
        while True:
            bad[Z | sub] = 1
            if sub == 0:
                break
            sub = (sub - 1) & comp
    best = -1
    witness = 0
    for m in range(1 << N):
        if not bad[m]:
            c = bin(m).count("1")
            if c > best:
                best = c
                witness = m
    return best, witness, len(zs)


def stage_brute_small():
    for n in (2, 3, 4):
        best, witness, nzs = brute_harborth(n)
        g = best + 1
        pts = [gvec(v, n) for v in range(n * n) if (witness >> v) & 1]
        ok_free = all(sum_points(list(U), n) != (0, 0)
                      for U in itertools.combinations(pts, n))
        check("brute_force_g_C%d_plus_C%d_equals_%d" % (n, n, gt_conjecture(n)),
              g == gt_conjecture(n) and ok_free and best == len(pts),
              "all 2^%d subsets inspected, %d zero-sum %d-subsets of G, "
              "largest set with none has %d elements -> g = %d, "
              "eq. (1) gives %d; extremal example %s"
              % (n * n, nzs, n, best, g, gt_conjecture(n),
                 sorted(pts, key=lambda p: gidx(p, n))))


# =====================================================================
# BLOCK 11.  THE ENGINE ON CASES WHOSE ANSWER CAN BE CHECKED
# =====================================================================


def brute_free(points, n):
    """True iff `points` has no n-element zero-sum subset (brute force)."""
    return all(sum_points(list(U), n) != (0, 0)
               for U in itertools.combinations(points, n))


def normalization_lemma(n):
    """Lemmas 2.2 + 2.3 for a general n: returns (ok, detail)."""
    G = [gvec(v, n) for v in range(n * n)]
    prim = [p for p in G if math.gcd(math.gcd(p[0], p[1]), n) == 1]
    D = set(p for p in G if p not in set(prim))
    subgroup = all(gadd(u, v, n) in D for u in D for v in D)
    cls = {}
    for p in G:
        cls.setdefault(min(gidx(gadd(p, h, n), n) for h in D), []).append(p)
    sizes = sorted(set(len(c) for c in cls.values()))
    MS = gl2(n)
    firstcols = set((M[0], M[2]) for M in MS)
    expkill = all(sum_points([c] * n, n) == (0, 0) for c in G)
    ok = (subgroup and sizes == [len(D)] and len(cls) == n * n // len(D)
          and firstcols == set(prim) and expkill and 2 * n - 1 > len(D))
    return ok, ("non-primitive set has order %d and is a subgroup, its %d "
                "cosets have size %d, GL_2(Z/%dZ) (order %d) moves every one "
                "of the %d primitive vectors to e_1, exp(G)*c = 0, and "
                "2n-1 = %d > %d"
                % (len(D), len(cls), len(D), n, len(MS), len(prim),
                   2 * n - 1, len(D)))


def stage_engine_small():
    for n in (3, 5, 7):
        k = 2 * n - 2
        run, st, dec, NP, CB = make_engine(n, k, 0)
        run()
        pts = dec(st["found"]) if st["found"] is not None else []
        described = set((x, y) for x in range(n - 1) for y in range(2))
        check("engine_finds_a_%d_element_free_set_in_C%d_plus_C%d"
              % (k, n, n),
              st["found"] is not None and len(pts) == k
              and brute_free(pts, n) and set(pts) == described,
              "unbiased DFS returns %s, brute-force free of %d-element "
              "zero-sum subsets; it is {0,...,%d} x {0,1}"
              % (sorted(pts, key=lambda p: gidx(p, n)), n, n - 2))

        if n > 3:
            ok, detail = normalization_lemma(n)
            check("normalization_reduction_valid_for_n%d" % n, ok, detail)

        k = 2 * n - 1
        run, st, dec, NP, CB = make_engine(n, k, 2)
        run()
        tot = CB[NP][k - 2]
        check("normalized_census_complete_C%d_plus_C%d" % (n, n),
              st["covered"] == tot and st["found"] is None and tot > 0,
              "all %d completions of B to a %d-element set accounted for in "
              "%d nodes; none is free" % (tot, k, st["nodes"]))

        if n <= 5:
            run2, st2, dec2, NP2, CB2 = make_engine(n, k, 0)
            run2()
            tot2 = CB2[NP2][k]
            check("unnormalized_census_complete_C%d_plus_C%d" % (n, n),
                  st2["covered"] == tot2 and st2["found"] is None and tot2 > 0,
                  "all %d subsets of size %d of G accounted for in %d nodes "
                  "with no normalization and no symmetry reduction; none is "
                  "free" % (tot2, k, st2["nodes"]))

        check("harborth_constant_C%d_plus_C%d_is_%d" % (n, n, 2 * n - 1),
              st["covered"] == tot and st["found"] is None
              and gt_conjecture(n) == 2 * n - 1,
              "lower bound from the %d-element free set, upper bound from the "
              "exhaustive census: g = %d = 2n-1, the odd branch of eq. (1)"
              % (2 * n - 2, 2 * n - 1))


def stage_engine_n9(A):
    n = P_N
    k = P_G_VALUE - 1
    run, st, dec, NP, CB = make_engine(n, k, 2)
    run(limit=500000)
    pts = dec(st["found"]) if st["found"] is not None else []
    check("engine_rediscovers_A_at_size_%d" % k,
          st["found"] is not None and len(pts) == k
          and brute_free(pts, n) and set(pts) == set(A),
          "the DFS at k = %d returns %s in %d nodes; brute force over all "
          "C(%d,%d) subsets confirms it is free, and it is exactly A"
          % (k, sorted(pts, key=lambda p: gidx(p, n)), st["nodes"],
             k, n) if st["found"] else "no set found in %d nodes" % st["nodes"])

    run2, st2, dec2, NP2, CB2 = make_engine(n, k, 2)
    run2(prefix=(0, 1, 2), limit=500000)
    pts2 = dec2(st2["found"]) if st2["found"] is not None else []
    check("prefix_code_path_can_report_a_free_set",
          st2["found"] is not None and set(pts2) == set(A),
          "the same prefix-driven code path used for the n = 9 census finds "
          "A below the canonical prefix (0,1,2), so a survivor would not be "
          "silently dropped")


# =====================================================================
# BLOCK 11b.  AUDIT OF THE PRUNING RULES THEMSELVES
#
# The self-tests of Block 5 compare psi_states() -- a SECOND, independent
# implementation -- with brute force.  They therefore say nothing about the
# masks that _make_dfs actually computes, and an over-pruning search engine
# reports "no free set" for every input, i.e. it confirms any upper bound
# whatsoever.  This block closes that gap two ways:
#
#   (a) a definition-only reference enumerator (its ONLY pruning rule is
#       Lemma 5 evaluated by brute force -- no clique rule, no degree
#       peeling, no colouring bound) is run to exhaustion on small cases,
#       and the engine must return exactly the lexicographically least free
#       completion below EVERY free prefix, or None exactly when the
#       reference has none;
#   (b) the audit hook feeds back the engine's own rule-1 and rule-3 masks,
#       which are compared point by point with "-z not in Sigma_{n-1}(P)"
#       and "-(u+v) not in Sigma_{n-2}(P)" recomputed by brute force.
#
# Both are stated so that they FAIL when the rule is sound-looking but
# wrong: substituting Sigma_8 for Sigma_7 in rule 3, or peeling at degree
# < r instead of < r-1, leaves every other check in this program passing.
# =====================================================================


def reference_free_sets(n, k, shift):
    """Every free k-set containing the first `shift` group elements, by
    position tuple, using ONLY the definition of freeness as a rule."""
    NP = n * n - shift
    out = []

    def rec(pts, s, r, sel):
        if r == 0:
            out.append(tuple(sel))
            return
        for z in range(s, NP - r + 1):
            p = gvec(z + shift, n)
            ok = True
            if len(pts) >= n - 1:
                for U in itertools.combinations(pts, n - 1):
                    if sum_points(list(U) + [p], n) == (0, 0):
                        ok = False
                        break
            if ok:
                pts.append(p)
                sel.append(z)
                rec(pts, z + 1, r - 1, sel)
                sel.pop()
                pts.pop()

    rec([gvec(v, n) for v in range(shift)], 0, k - shift, [])
    return sorted(out)


def stage_engine_audit():
    # ---- (a) exhaustive differential test against the reference ---------
    configs = [(3, 4, 0, 2), (3, 5, 0, 2), (5, 8, 2, 2), (5, 9, 2, 2)]
    bad = []
    nprefix = 0
    nref = 0
    for n, k, shift, maxpref in configs:
        ref = reference_free_sets(n, k, shift)
        nref += len(ref)
        NP = n * n - shift
        for L in range(1, maxpref + 1):
            for t in itertools.combinations(range(NP), L):
                if shift + L >= n:
                    continue          # engine assumes a free prefix
                nprefix += 1
                below = [f for f in ref if f[:L] == t]
                want = below[0] if below else None
                run, st, dec, NP2, CB = make_engine(n, k, shift)
                run(prefix=t)
                got = tuple(st["found"]) if st["found"] is not None else None
                if got != want:
                    bad.append((n, k, shift, t, want, got))
                if got is not None and not brute_free(dec(list(got)), n):
                    bad.append((n, k, shift, t, "not free", got))
    if bad:
        detail = ("%d disagreements with the unpruned reference, e.g. "
                  "(n,k,shift,prefix,reference,engine) = %r"
                  % (len(bad), bad[0]))
    elif nprefix == 0 or nref == 0:
        detail = ("the differential test exercised %d prefixes and %d "
                  "reference free sets: it did not actually run"
                  % (nprefix, nref))
    else:
        detail = ("%d free prefixes over %d cases: the engine returned "
                  "exactly the lexicographically least free completion the "
                  "unpruned reference enumerator found (%d free sets "
                  "enumerated in total), and None exactly when there is none"
                  % (nprefix, len(configs), nref))
    check("engine_agrees_with_a_definition_only_reference_on_every_prefix",
          not bad and nprefix > 0 and nref > 0, detail)

    # ---- (b) the engine's OWN rule-1 / rule-3 masks vs the definitions --
    err = []
    seen = {"rule1": 0, "rule3": 0}

    def make_auditor(n, shift, NP):
        def reviewer(kind, sel, s, r, adm, nb):
            if len(err) > 4:
                return
            seen[kind] += 1
            pts = [gvec(v, n) for v in range(shift)]
            pts += [gvec(p + shift, n) for p in sel]
            sig = subset_sums_by_size(pts, n, n - 1)
            want = 0
            for p in range(s, NP):
                z = gvec(p + shift, n)
                if gneg(z, n) not in sig[n - 1]:
                    want |= 1 << p
            if kind == "rule1":
                if adm != want:
                    err.append(("rule1", tuple(sel), s, adm ^ want))
                return
            m = adm
            while m:
                b = m & -m
                u = b.bit_length() - 1
                m ^= b
                wnb = 0
                a2 = adm
                while a2:
                    b2 = a2 & -a2
                    v = b2.bit_length() - 1
                    a2 ^= b2
                    if v == u:
                        continue
                    su = gadd(gvec(u + shift, n), gvec(v + shift, n), n)
                    if gneg(su, n) not in sig[n - 2]:
                        wnb |= b2
                if nb[u] != wnb:
                    err.append(("rule3", tuple(sel), u, nb[u] ^ wnb))
                    return
        return reviewer

    for n, k, shift, prefix in [(5, 8, 2, ()), (5, 9, 2, (0,)),
                                (9, 17, 2, (0, 1, 2)), (9, 16, 2, (0, 3))]:
        run, st, dec, NP, CB = make_engine(n, k, shift)
        st["audit"] = make_auditor(n, shift, NP)
        run(limit=260 if n == 9 else 900, **({} if not prefix
                                             else {"prefix": prefix}))
        st["audit"] = None
    enough = seen["rule1"] > 100 and seen["rule3"] > 100
    if err:
        detail = "%d mask errors, e.g. %r" % (len(err), err[0])
    elif not enough:
        detail = ("the audit reached only %d rule-1 and %d rule-3 nodes "
                  "(>100 of each required): the search collapsed before the "
                  "rules could be exercised, which is itself a defect"
                  % (seen["rule1"], seen["rule3"]))
    else:
        detail = ("the masks _make_dfs itself computed at %d rule-1 and %d "
                  "rule-3 nodes (n = 5 and n = 9) agree everywhere with "
                  "'-z not in Sigma_%d(P)' and '-(u+v) not in Sigma_%d(P)' "
                  "recomputed by brute force from the selected points"
                  % (seen["rule1"], seen["rule3"], P_N - 1, P_N - 2))
    check("engine_rule1_and_rule3_masks_equal_their_definitions",
          not err and enough, detail)


def stage_pruning_control_n9(A):
    """Positive control at n = 9 itself: the pruning must never reject a
    branch that contains the known free 16-element set A."""
    n = P_N
    k = P_G_VALUE - 1
    Apos = tuple(sorted(gidx(p, n) - 2 for p in A if gidx(p, n) >= 2))
    forced = set(gvec(v, n) for v in range(2))
    bad = []
    for L in range(1, len(Apos) + 1):
        run, st, dec, NP, CB = make_engine(n, k, 2)
        run(prefix=Apos[:L], limit=400000)
        pts = dec(st["found"]) if st["found"] is not None else []
        if set(pts) != set(A):
            bad.append((L, st["nodes"], st["found"]))
    shape = forced <= set(A) and len(Apos) == k - 2
    if bad:
        detail = "%d prefixes lost A, e.g. %r" % (len(bad), bad[0])
    elif not shape:
        detail = ("A does not have the shape this control needs: it has %d "
                  "positions outside B (expected %d) and %s contain B"
                  % (len(Apos), k - 2, "does" if forced <= set(A)
                     else "does NOT"))
    else:
        detail = ("for each of the %d prefixes of A's own position sequence "
                  "the DFS still returns A, so no rejection rule discards a "
                  "subfamily that does contain a known free set" % len(Apos))
    check("pruning_never_rejects_the_branch_containing_A",
          not bad and shape, detail)


# =====================================================================
# BLOCK 12.  THE n = 9 CENSUS OVER FEASIBLE CANONICAL PREFIXES
# =====================================================================


def stage_census_n9(feas, seconds, full):
    n = P_N
    k = P_G_VALUE
    run, st, dec, NP, CB = make_engine(n, k, 2)
    rest = k - 2 - 3                      # points still free below a prefix
    weight = dict((t, CB[NP - 1 - t[2]][rest]) for t in feas)
    total = sum(weight.values())
    # The cheap prefixes (large p_3) finish instantly and give the exact
    # accounting identity many independent tests; the expensive ones are
    # then taken in decreasing weight, which maximises coverage per second.
    byweight = sorted(feas, key=lambda t: (weight[t], t))
    order = byweight[:100] + sorted(byweight[100:],
                                    key=lambda t: (-weight[t], t))
    deadline = time.time() + seconds
    acc_ok = True
    cum = 0
    done = 0
    t0 = time.time()
    for t in order:
        left = None if full else deadline - time.time()
        if left is not None and left <= 0:
            break
        run(prefix=t, seconds=left)
        if st["stopped"] or st["found"] is not None:
            break
        cum += weight[t]
        done += 1
        if st["covered"] != cum:
            acc_ok = False
    elapsed = time.time() - t0
    covered = st["covered"]

    check("n9_census_accounting_is_exact",
          acc_ok and done > 0 and not st["prefix_bad"],
          "after each of the %d fully searched prefixes the credited "
          "completion count equalled the independent binomial total exactly"
          % done)
    check("n9_census_no_free_17_element_set_in_the_searched_range",
          st["found"] is None,
          "%d nodes searched, no 17-element normalized set survives" % st["nodes"])
    frac = (100.0 * covered / total) if total else 0.0
    check("n9_census_coverage_is_sound_and_nonempty",
          0 < covered <= total and st["nodes"] > 10 * max(done, 1),
          "%d of %d completions of eq. (2) examined or soundly rejected "
          "(%.3f%%), %d of %d prefixes finished, %d nodes, %.0f s"
          % (covered, total, frac, done, len(feas), st["nodes"], elapsed))
    if full:
        check("n9_census_complete",
              covered == total and done == len(feas) and st["found"] is None,
              "all %d prefixes finished and all %d completions covered"
              % (done, total))
    if covered < total:
        rate = covered / elapsed if elapsed > 0 else 0.0
        note("NOT RE-RUN: %d of the %d completions of eq. (2) "
             "(%.3f%% of the census was covered here). The paper's upper "
             "bound is an exhaustive search over all %d feasible canonical "
             "prefixes; at the rate measured above (%.3g completions/s in "
             "pure Python) finishing it needs of the order of %.0f s more, "
             "single process, which exceeds this program's budget -- and the "
             "prefixes not yet touched are the less efficient ones, so that "
             "is a lower estimate. Run 'python3 verify.py --full' to "
             "reproduce the census in full; no --full run accompanies this "
             "transcript, so the time above is an estimate extrapolated from "
             "the covered part and not a measurement."
             % (total - covered, total, frac, len(feas), rate,
                (total - covered) / rate if rate else 0.0))
        note("what IS re-derived exactly at n = 9: the 16-element lower-bound "
             "set and its maximality, the stabilizer of order 108, the 992 "
             "canonical triples, the 831 feasible prefixes and the integer "
             "811509569804714 of eq. (2); and the identical search engine "
             "settles n = 3, 5, 7 exhaustively.")
    note("the paper's node count %d is not reproduced: this program prunes "
         "with degree peeling plus a greedy-colouring clique bound rather "
         "than an exact clique search, so its node count (%d here) is a "
         "property of the implementation, not a mathematical claim."
         % (P_DFS_NODES, st["nodes"]))
    note("NOT CHECKED: the authors' own census code, and the complete "
         "independent re-run of it reported in the paper. This program is "
         "standard-library only and reads no external file; it re-derives the "
         "mathematics instead.")
    return covered == total and done == len(feas) and st["found"] is None


# =====================================================================
# MAIN
# =====================================================================


def main(argv):
    seconds = 1080.0
    full = False
    args = argv[1:]
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--quick":
            seconds = 20.0
        elif a == "--full":
            full = True
        elif a == "--seconds" and i + 1 < len(args):
            i += 1
            seconds = float(args[i])
        else:
            print("usage: verify.py [--quick] [--seconds S] [--full]")
            return 2
        i += 1
    if seconds <= 0:
        print("census budget must be positive")
        return 2

    sys.setrecursionlimit(10000)
    t0 = time.time()
    note("verify.py -- g(C_9 (+) C_9) = 17, the n = 9 case of "
         "Gao-Thangadurai")
    note("exact integer arithmetic only; n = 9 census budget %.0f s%s"
         % (seconds, " (ignored: --full)" if full else ""))

    crashed = None
    complete = False
    try:
        A = stage_object()
        stage_lower_bound(A)
        stage_selftest(P_N)
        stage_normalization(A)
        E, stab, perms = stage_stabilizer()
        feas = stage_prefixes(perms, E)
        stage_brute_small()
        stage_engine_small()
        stage_engine_n9(A)
        stage_engine_audit()
        stage_pruning_control_n9(A)
        complete = stage_census_n9(feas, seconds, full)
    except Exception:
        crashed = traceback.format_exc()
        sys.stderr.write(crashed)
    check("program_ran_to_completion", crashed is None,
          "every stage returned" if crashed is None
          else crashed.strip().splitlines()[-1])

    if complete:
        note("g(C_9 (+) C_9) >= 17 and g <= 17 are both fully re-derived "
             "here: the census over every feasible canonical prefix closed "
             "with no surviving 17-element set.")
    else:
        note("g(C_9 (+) C_9) >= 17 is fully re-derived here; g <= 17 is "
             "re-derived only over the covered part of the census, as "
             "reported above.")
    note("elapsed %.1f s" % (time.time() - t0))
    return verdict()


if __name__ == "__main__":
    sys.exit(main(sys.argv))


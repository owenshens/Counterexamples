#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_vr_q6_r4.py -- ancillary verification for the Vietoris-Rips complex
K = X^{6,4} = VR(Q_6; 4) and the equivariant Euler-characteristic computation
Psi(g) = I(F_g; -1) over the 65 conjugacy classes of H_6 = F_2^6 : S_6.

Python 3.9, standard library ONLY.  Exact integer arithmetic throughout;
no floating point is used to decide any assertion.

---------------------------------------------------------------------------
TAKEN FROM THE PAPER (inputs, not re-derived here)
---------------------------------------------------------------------------
  P1. Parameters n = 6, r = 4, and the definition
        VR(Q_n; r) = { sigma subseteq Q_n : d(x,y) <= r for all x,y in sigma }
      (purely PAIRWISE; empty set and singletons are faces).
  P2. The group G = H_6 = Aut(Q_6) = F_2^6 semidirect S_6 acting by
      x |-> a + pi x, with |G| = 46080 (also re-derived as 2^6 * 6! in C2).
  P3. LEMMA lem:fixed-orbit, ASSUMED and NOT reproved:
        Psi(g) = I(F_g; -1) = J(F_g)  for every g in G,
      together with the definition of O_g (the g-orbits on Q_6 of internal
      Hamming diameter <= 4; larger orbits are DISCARDED) and the edge rule
      {O,O'} in E(F_g) iff diam(O union O') > 4.
      This program computes only the right-hand side; it never touches
      homology.
  P4. The recurrence (eq:recurrence) J(F) = J(F - v) - J(F - N[v]) with
      J(empty graph) = 1, and multiplicativity of J over connected
      components.  This is the prescribed algorithm; C6 tests our
      implementation of it against independent brute force anyway.
  P5. The class-size formula (eq:class-size)
        |C_{alpha,beta}| = 2^6 * 6! / ( 2^{l(alpha)+l(beta)} z_alpha z_beta )
      with z_lambda = prod_j j^{m_j} m_j! and l = number of parts.
  P6. The representative recipe: lay the parts of alpha then the parts of
      beta on consecutive coordinate blocks of {0..5}; each block of length L
      becomes an L-cycle of pi; a = 0 on every alpha block; exactly one
      coordinate of a set to 1 on every beta block.
  P7. (eq:homology), cited by the paper from GMW Appendix C.3:
        reduced H_i(X^{6,4}; Q) = Q^239 (i = 7), Q^14 (i = 15), 0 otherwise.
      Used ONLY in C9 (as the sum 239 + 14 = 253) and to license the
      character-averaging step.  NOT verified here -- see the LIMITATIONS
      block printed at the end of the run.
  P8. Character averaging: (1/|G|) sum_{g in G} tr(g|V) = dim V^G for a
      finite-dimensional Q[G]-representation V.
  P9. THE WITNESS: the full 65-row Appendix A table (type, |C|, Psi),
      transcribed verbatim below as WITNESS_TABLE_TEXT.  This is the object
      under test -- it is a CLAIM, never an input to any computation.

---------------------------------------------------------------------------
DERIVED FROM SCRATCH BY THIS PROGRAM (nothing below is read off the paper)
---------------------------------------------------------------------------
  D1. All bipartitions (alpha; beta) of 6, and the count 65 itself.
  D2. Each representative (a, pi) explicitly, plus recovery of the
      (alpha,beta) invariant from it by cycle decomposition and a-parity.
  D3. Every class size from (eq:class-size), the total 46080, and the
      identity 46080 = 2^6 * 6!.
  D4. The g-orbits on Q_6 for all 65 classes, their internal Hamming
      diameters, the retained vertex set O_g and |O_g|.
  D5. Every adjacency of F_g, from actual Hamming diameters of orbit unions.
  D6. All 65 values Psi = J(F_g) by the recurrence, in exact integers.
  D7. An independent brute-force I(F_g;-1) for every class with |O_g| <= 20,
      plus small-graph unit tests of the recurrence.
  D8. The sums S1 = sum |C| and S2 = sum |C| * Psi, and S2 == 3 * S1.
  D9. Psi(identity) = 253, compared against 239 + 14 (P7).
  D10. A FULL CENSUS of all 6! * 2^6 = 46080 elements of G (C12): each
      element classified by its own recovered (alpha;beta) invariant, so
      (eq:class-size) is checked CLASS BY CLASS rather than only in total,
      sum_{g in G} Psi(g) = 138240 is recomputed over ELEMENTS with no appeal
      to the class-size formula, and Psi is shown constant on each type block
      by re-deriving F_g and J on extra non-representative elements.
  D11. A second exact evaluation of every J(F_g) by an independent recurrence
      driver (minimum-degree branch vertex, no component multiplicativity),
      which is what covers the classes too large for enumerative brute force
      (C6b).

---------------------------------------------------------------------------
PROVENANCE OF EVERY HARDCODED PAPER NUMBER  (each is only ever the RIGHT-hand
side of a comparison against a derived value; none is fed into a computation)
---------------------------------------------------------------------------
  N = 6, R = 4                      the paper's parameter point (P1); C11
                                    checks (6,4) is the minimal admissible pair
  PAPER_H7_DIM = 239                (eq:homology); C9 only
  PAPER_H15_DIM = 14                (eq:homology); C9 only
  PAPER_HOMOLOGY_DEGREES = (7, 15)  (eq:homology); C9 (signs) and C10 (prose)
  PAPER_GROUP_ORDER = 46080         |G|; C2 vs the derived class-size total
  PAPER_S1 = 46080                  Proposition; C2, C8 vs derived sums
  PAPER_S2 = 138240                 Proposition; C8 and C12 vs derived sums
  PAPER_AVERAGE = 3                 Theorem; C8 vs divmod(S2, S1)
  bare 65                           C1, C7, C12 vs the generated class count
  bare 253                          C9 vs the derived Psi at the identity
  bare 3 in "S2 == 3 * S1"          C8 vs derived S1, S2
  "111111|-"                        C9/C12 label of the identity class
  WITNESS_TABLE_TEXT (130 numbers)  Appendix A, verbatim; C7 only
  the 10 hand J values in
      unit_tests_recurrence()       textbook values of I(F;-1) on tiny graphs
That is the complete list.

NOT DERIVABLE HERE (stated as limitations, not verified):
  the individual numbers 239 and 14; the claim that reduced homology of
  X^{6,4} vanishes outside degrees {7, 15}; Lemma lem:fixed-orbit itself;
  and WHICH of degree 7 or 15 carries the repeated trivial constituent
  (the paper states this is undetermined -- C10 checks that the inference
  stops exactly there).
"""

import sys
import math

# ---------------------------------------------------------------------------
# Fixed parameters of the paper's computation (P1).
# ---------------------------------------------------------------------------
N = 6                    # cube dimension
R = 4                    # Vietoris-Rips parameter
NPTS = 1 << N            # 64 points of Q_6

# ---------------------------------------------------------------------------
# P9 -- THE WITNESS.  Appendix A, Table 1: "Exact values of Psi(g)=I(F_g;-1)
# on all conjugacy classes."  65 rows, in the paper's reading order
# (left column-group, middle, right; row by row).  Format: type  |C|  Psi.
# This is the ONLY large block of paper numbers the program holds, and it is
# used solely as the claim to be checked in C7.
# ---------------------------------------------------------------------------
WITNESS_TABLE_TEXT = """
-|111111      1     1
-|21111      30     1
-|2211      180    17
-|222       120     1
-|3111      160     1
-|321       960     1
-|33        640     1
-|411       720     1
-|42       1440     1
-|51       2304     1
-|6        3840     1
1|11111       6     1
1|2111      120     5
1|221       360    17
1|311       480     1
1|32        960     5
1|41       1440     5
1|5        2304     1
11|1111      15    53
11|211      180     5
11|22       180     5
11|31       480     5
11|4        720     5
111|111      20    53
111|21      120     9
111|3       160     5
1111|11      15     9
1111|2       30     9
11111|1       6     9
111111|-      1   253
2|1111       30     3
2|211       360    -1
2|22        360     3
2|31        960     3
2|4        1440    -1
21|111      120     3
21|21       720     3
21|3        960     3
211|11      180     7
211|2       360     3
2111|1      120     7
21111|-      30    11
22|11       180     9
22|2        360     1
221|1       360     9
2211|-      180    29
222|-       120     3
3|111       160    -1
3|21        960     3
3|3        1280     5
31|11       480     3
31|2        960     3
311|1       480     3
3111|-      160     7
32|1        960     1
321|-       960     5
33|-        640    13
4|11        720     3
4|2        1440     3
41|1       1440     3
411|-       720     7
42|-       1440     5
5|1        2304    -1
51|-       2304     3
6|-        3840     3
"""

# P7 -- (eq:homology) from GMW Appendix C.3, cited by the paper.
PAPER_H7_DIM = 239
PAPER_H15_DIM = 14
PAPER_HOMOLOGY_DEGREES = (7, 15)

# P2 -- the paper's stated group order (independently re-derived in C2).
PAPER_GROUP_ORDER = 46080

# The paper's two claimed sums (checked in C8, computed from derived data).
PAPER_S1 = 46080
PAPER_S2 = 138240
PAPER_AVERAGE = 3


# ===========================================================================
# Reporting harness
# ===========================================================================
RESULTS = []          # list of (name, ok, detail)


def check(name, ok, detail=""):
    """Record one named check.  Returns ok so callers can branch."""
    RESULTS.append((name, bool(ok), detail))
    print(("PASS " if ok else "FAIL ") + name + ((" -- " + detail) if detail else ""))
    return bool(ok)


def info(msg):
    """Intermediate quantity, printed so a human can follow the derivation."""
    print("      | " + msg)


def section(title):
    print("")
    print("=== " + title + " " + "=" * max(0, 68 - len(title)))


# ===========================================================================
# D1 -- partitions and bipartitions, generated from scratch
# ===========================================================================
def partitions(n, maxpart=None):
    """All partitions of n as tuples of parts in weakly DECREASING order."""
    if maxpart is None:
        maxpart = n
    if n == 0:
        return [()]
    out = []
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, first):
            out.append((first,) + rest)
    return out


def partition_count_pentagonal(n):
    """p(0..n) by Euler's pentagonal-number recurrence
        p(m) = sum_{k>=1} (-1)^{k-1} [ p(m - k(3k-1)/2) + p(m - k(3k+1)/2) ].

    This shares NO code with the `partitions` generator above, so C1 can use it
    to TEST that generator.  (The check it replaces compared
    len(all_bipartitions(6)) with sum_k len(partitions(k))*len(partitions(6-k)),
    which is an identity of all_bipartitions' own double loop and therefore
    could not fail for any implementation of `partitions`.)"""
    p = [0] * (n + 1)
    p[0] = 1
    for m in range(1, n + 1):
        total = 0
        k = 1
        while True:
            g1 = k * (3 * k - 1) // 2
            g2 = k * (3 * k + 1) // 2
            if g1 > m and g2 > m:
                break
            sign = 1 if (k % 2) else -1
            if g1 <= m:
                total += sign * p[m - g1]
            if g2 <= m:
                total += sign * p[m - g2]
            k += 1
        p[m] = total
    return p


def type_string(alpha, beta):
    """Appendix notation: parts concatenated weakly decreasing, '-' for empty,
    the two halves joined by '|'.  e.g. ((2,2);(1,1)) -> '22|11'."""
    left = "".join(str(p) for p in alpha) if alpha else "-"
    right = "".join(str(p) for p in beta) if beta else "-"
    return left + "|" + right


def all_bipartitions(n):
    """Every (alpha; beta) with |alpha| + |beta| = n, alpha and beta
    partitions.  Generated -- the count is NOT read from the paper."""
    out = []
    for k in range(0, n + 1):
        for alpha in partitions(k):
            for beta in partitions(n - k):
                out.append((alpha, beta))
    return out


def z_lambda(lam):
    """z_lambda = prod_j j^{m_j} * m_j!  (P5)."""
    mult = {}
    for p in lam:
        mult[p] = mult.get(p, 0) + 1
    z = 1
    for j, m in mult.items():
        z *= (j ** m) * math.factorial(m)
    return z


def class_size(alpha, beta):
    """(eq:class-size), with an EXACT-division assertion (P5)."""
    num = (2 ** N) * math.factorial(N)
    den = (2 ** (len(alpha) + len(beta))) * z_lambda(alpha) * z_lambda(beta)
    if num % den != 0:
        raise AssertionError(
            "class-size division not exact for %s: %d / %d"
            % (type_string(alpha, beta), num, den))
    return num // den


# ===========================================================================
# D2 -- representatives (a, pi) by the paper's block recipe (P6),
#       and INDEPENDENT recovery of (alpha, beta) out of them (C3)
# ===========================================================================
def build_representative(alpha, beta):
    """Return (a, perm) where perm[i] = pi(i) on coordinates {0..N-1} and a is
    a bit vector packed as an int (bit j = a_j).

    Recipe (P6): parts of alpha then parts of beta on CONSECUTIVE coordinate
    blocks; each block of length L becomes an L-cycle; a = 0 on alpha blocks;
    exactly one coordinate (the first) of a set to 1 on each beta block.
    A length-1 beta block therefore puts a 1 in that single coordinate."""
    perm = [None] * N
    a = 0
    c = 0
    for L in alpha:
        for t in range(L):
            perm[c + t] = c + ((t + 1) % L)
        c += L
    for L in beta:
        for t in range(L):
            perm[c + t] = c + ((t + 1) % L)
        a |= 1 << c            # single 1 on the first coordinate of the block
        c += L
    if c != N:
        raise AssertionError("block layout covered %d of %d coordinates" % (c, N))
    if sorted(perm) != list(range(N)):
        raise AssertionError("block layout is not a permutation: %r" % (perm,))
    return a, tuple(perm)


def cycles_of(perm):
    """Cycle decomposition of perm as a list of lists of coordinates."""
    seen = [False] * len(perm)
    out = []
    for i in range(len(perm)):
        if seen[i]:
            continue
        cyc = []
        j = i
        while not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = perm[j]
        out.append(cyc)
    return out


def recover_type(a, perm):
    """Recompute the bipartition invariant from (a, pi): cycles of pi with
    EVEN a-parity go to alpha, ODD a-parity to beta; both sorted decreasing.
    Also returns the total of the cycle lengths for a covering check."""
    ev, od = [], []
    total = 0
    for cyc in cycles_of(perm):
        parity = 0
        for j in cyc:
            parity ^= (a >> j) & 1
        (od if parity else ev).append(len(cyc))
        total += len(cyc)
    return tuple(sorted(ev, reverse=True)), tuple(sorted(od, reverse=True)), total


# ===========================================================================
# Q_6 geometry: action, Hamming distance, orbits, diameters
# ===========================================================================
def popcount(x):
    return bin(x).count("1")


def hamming(x, y):
    return popcount(x ^ y)


def apply_g(a, perm, x):
    """x |-> a + pi x with (pi x)_{pi(j)} = x_j, i.e. (pi x)_i = x_{pi^{-1}(i)}.
    The spec notes the other convention gives a conjugate of g^{-1}, which has
    the same orbits up to Hamming isometry, so every derived quantity is
    unchanged."""
    y = 0
    for j in range(N):
        if (x >> j) & 1:
            y |= 1 << perm[j]
    return y ^ a


def orbits_of(a, perm):
    """D4: the g-orbits on Q_6, each as a tuple listing the cycle in the order
    visited starting from its least element."""
    seen = [False] * NPTS
    out = []
    for x0 in range(NPTS):
        if seen[x0]:
            continue
        orb = []
        x = x0
        while not seen[x]:
            seen[x] = True
            orb.append(x)
            x = apply_g(a, perm, x)
        out.append(tuple(orb))
    return out


def diameter(points):
    """Internal Hamming diameter; diam(singleton) = 0, diam(empty) = 0."""
    best = 0
    pts = tuple(points)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            d = hamming(pts[i], pts[j])
            if d > best:
                best = d
    return best


def cross_max(o1, o2):
    """max over x in o1, y in o2 of d(x,y)."""
    best = 0
    for x in o1:
        for y in o2:
            d = hamming(x, y)
            if d > best:
                best = d
    return best


# ===========================================================================
# D5 -- the graph F_g
# ===========================================================================
def build_Fg(a, perm):
    """Vertices = the g-orbits of internal diameter <= R (O_g); orbits of
    diameter >= R+1 are DISCARDED and are NOT vertices.  Edge {O,O'} for
    O != O' iff diam(O union O') > R.

    The O_g filter is load-bearing: applying the union rule to unfiltered
    orbits would create loops (diam(O union O) = diam(O) > R) and give wrong
    answers.  We therefore build the vertex set FIRST and only ever join
    distinct retained orbits.

    Returns dict: kept (list of orbit tuples), dropped (list of orbit tuples),
    kept_diam, dropped_diam, adj (list of int bitmasks)."""
    orbs = orbits_of(a, perm)
    kept, dropped, kdiam, ddiam = [], [], [], []
    for o in orbs:
        d = diameter(o)
        if d <= R:
            kept.append(o)
            kdiam.append(d)
        else:
            dropped.append(o)
            ddiam.append(d)
    m = len(kept)
    adj = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if diameter(kept[i] + kept[j]) > R:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return {"kept": kept, "dropped": dropped, "kept_diam": kdiam,
            "dropped_diam": ddiam, "adj": adj, "n": m}


# ===========================================================================
# D6 -- J(F) = I(F;-1) via (eq:recurrence), memoised, component-multiplicative
# ===========================================================================
def components(mask, adj):
    """Connected components of the induced subgraph on `mask`, as bitmasks."""
    comps = []
    rem = mask
    while rem:
        v = (rem & -rem).bit_length() - 1
        comp = 1 << v
        frontier = 1 << v
        while frontier:
            newf = 0
            f = frontier
            while f:
                u = (f & -f).bit_length() - 1
                f &= f - 1
                newf |= adj[u] & mask & ~comp
            comp |= newf
            frontier = newf
        comps.append(comp)
        rem &= ~comp
    return comps


def J_value(mask, adj, memo):
    """J(F[mask]) = I(F[mask]; -1) by (eq:recurrence) (P4).

      J(empty graph) = 1;
      J multiplies over connected components;
      inside a component, branch on a MAXIMUM-DEGREE vertex v:
          J = J(F - v) - J(F - N[v])   with N[v] the CLOSED neighbourhood.

    A one-vertex component falls out of the same branch:
    J = J(empty) - J(empty) = 0.  Exact integers only."""
    if mask == 0:
        return 1
    if mask in memo:
        return memo[mask]
    comps = components(mask, adj)
    if len(comps) > 1:
        r = 1
        for c in comps:
            r *= J_value(c, adj, memo)
    else:
        best_v, best_d = -1, -1
        m = mask
        while m:
            v = (m & -m).bit_length() - 1
            m &= m - 1
            d = popcount(adj[v] & mask)
            if d > best_d:
                best_d, best_v = d, v
        minus_v = mask & ~(1 << best_v)
        minus_closed = mask & ~((1 << best_v) | adj[best_v])
        r = J_value(minus_v, adj, memo) - J_value(minus_closed, adj, memo)
    memo[mask] = r
    return r


def J_of_graph(adj):
    """J of the whole graph given by adjacency bitmask list."""
    n = len(adj)
    return J_value((1 << n) - 1 if n else 0, adj, {})


class _BudgetExceeded(Exception):
    """Raised by J_value_alt when its memo exceeds the allowed size."""


def J_value_alt(mask, adj, memo, budget):
    """A SECOND exact evaluator of J = I(F;-1), deliberately different from
    J_value so that the two can be cross-checked on the LARGE classes where
    subset/independent-set enumeration is infeasible:

      * it branches on a MINIMUM-degree vertex, not a maximum-degree one;
      * it does NOT use multiplicativity over connected components at all,
        so a bug in `components` cannot be shared between the two;
      * it terminates only via the recurrence itself and the mask == 0 base.

    Same recurrence (eq:recurrence), same exact integer arithmetic.  Raises
    _BudgetExceeded if the memo grows past `budget`, so a class that is too
    expensive is REPORTED as un-cross-checked instead of silently skipped."""
    if mask == 0:
        return 1
    got = memo.get(mask)
    if got is not None:
        return got
    if len(memo) > budget:
        raise _BudgetExceeded()
    best_v, best_d = -1, -1
    m = mask
    while m:
        v = (m & -m).bit_length() - 1
        m &= m - 1
        d = popcount(adj[v] & mask)
        if best_v < 0 or d < best_d:
            best_d, best_v = d, v
    minus_v = mask & ~(1 << best_v)
    minus_closed = mask & ~((1 << best_v) | adj[best_v])
    r = (J_value_alt(minus_v, adj, memo, budget)
         - J_value_alt(minus_closed, adj, memo, budget))
    memo[mask] = r
    return r


def J_of_graph_alt(adj, budget):
    """J via the alternate evaluator.  Returns (value_or_None, ok, nodes) where
    nodes is the number of memo entries consumed, so a caller can enforce a
    TOTAL budget across all 65 classes and keep the whole run bounded."""
    n = len(adj)
    memo = {}
    try:
        v = J_value_alt((1 << n) - 1 if n else 0, adj, memo, budget)
        return v, True, len(memo)
    except _BudgetExceeded:
        return None, False, len(memo)


# ===========================================================================
# D7 -- two INDEPENDENT brute forces for I(F;-1) (no recurrence, no memo)
# ===========================================================================
def brute_I_minus1_by_subsets(adj):
    """Enumerate EVERY subset of V(F), test independence, sum (-1)^{|S|}.
    Includes the empty set.  Only used for small vertex counts."""
    n = len(adj)
    total = 0
    for s in range(1 << n):
        ok = True
        m = s
        while m:
            v = (m & -m).bit_length() - 1
            m &= m - 1
            if adj[v] & s:
                ok = False
                break
        if ok:
            total += -1 if (popcount(s) & 1) else 1
    return total


def brute_I_minus1_by_extension(adj):
    """Enumerate only the independent sets, by recursive extension over the
    vertex order, and accumulate (-1)^{|S|}.  Structurally different from both
    the recurrence and the subset scan."""
    n = len(adj)

    def rec(idx, forbidden):
        if idx == n:
            return 1
        # branch 1: idx not in S
        total = rec(idx + 1, forbidden)
        # branch 2: idx in S (allowed iff not forbidden by an earlier choice)
        if not ((forbidden >> idx) & 1):
            total -= rec(idx + 1, forbidden | adj[idx])
        return total

    return rec(0, 0)


def independent_set_polynomial(adj):
    """I(F;t) coefficient list [c_0, c_1, ...], c_k = # independent k-sets,
    by explicit enumeration.  Used to print I(F;-1) as an alternating sum for
    the identity class, so a reader can audit the sign convention."""
    n = len(adj)
    coeffs = [0] * (n + 1)

    def rec(idx, forbidden, size):
        coeffs[size] += 1
        for v in range(idx, n):
            if not ((forbidden >> v) & 1):
                rec(v + 1, forbidden | adj[v] | (1 << v), size + 1)

    rec(0, 0, 0)
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


# ===========================================================================
# The witness, parsed.  Used ONLY as the claim in C7 (and for the structural
# audit of the table itself).  Never fed into any computation.
# ===========================================================================
def parse_witness():
    """Return (rows, dup_types) where rows maps type string -> (|C|, Psi)."""
    rows = {}
    dups = []
    order = []
    for line in WITNESS_TABLE_TEXT.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 3:
            raise AssertionError("malformed witness row: %r" % (line,))
        t, c, psi = parts[0], int(parts[1]), int(parts[2])
        if t in rows:
            dups.append(t)
        rows[t] = (c, psi)
        order.append(t)
    return rows, dups, order


def parse_type_string(t):
    """Inverse of type_string: '22|11' -> ((2,2),(1,1)); '-' -> ().
    Returns None if the string is not a well-formed bipartition of N."""
    if t.count("|") != 1:
        return None
    left, right = t.split("|")

    def half(s):
        if s == "-":
            return ()
        if not s.isdigit() or "0" in s:
            return None
        return tuple(int(ch) for ch in s)

    al, be = half(left), half(right)
    if al is None or be is None:
        return None
    if sum(al) + sum(be) != N:
        return None
    if list(al) != sorted(al, reverse=True) or list(be) != sorted(be, reverse=True):
        return None
    return (al, be)


# ===========================================================================
# Per-class derivation, carrying the C3 / C4 / C5 validations
# ===========================================================================
def derive_one(alpha, beta):
    """Build everything for one class.  Returns (record, problems) where
    problems is a list of (check_name, message) for validations that failed."""
    probs = []
    t = type_string(alpha, beta)
    a, perm = build_representative(alpha, beta)

    # --- C3: the built (a,pi) really has type (alpha;beta) ---
    ral, rbe, tot = recover_type(a, perm)
    if (ral, rbe) != (alpha, beta):
        probs.append(("C3", "%s: recovered %s" % (t, type_string(ral, rbe))))
    if tot != N:
        probs.append(("C3", "%s: cycle lengths sum to %d, not %d" % (t, tot, N)))

    # --- C4: orbits form a genuine partition of Q_6, each a true g-cycle ---
    orbs = orbits_of(a, perm)
    seen = set()
    total = 0
    for o in orbs:
        total += len(o)
        for x in o:
            if x in seen:
                probs.append(("C4", "%s: point %d in two orbits" % (t, x)))
            seen.add(x)
        if apply_g(a, perm, o[-1]) != o[0]:
            probs.append(("C4", "%s: orbit %r does not close up" % (t, o)))
        for k in range(len(o) - 1):
            if apply_g(a, perm, o[k]) != o[k + 1]:
                probs.append(("C4", "%s: orbit %r not a g-trajectory" % (t, o)))
    if total != NPTS or len(seen) != NPTS:
        probs.append(("C4", "%s: orbit sizes total %d, %d distinct points"
                      % (t, total, len(seen))))

    # C4, independent recomputation: materialise g as a permutation P of the 64
    # points, check it IS a bijection (so g really is a group element), and
    # recover its cycles with the generic cycle routine.  This does not reuse
    # orbits_of, so an error in that routine cannot hide here.
    P = [apply_g(a, perm, x) for x in range(NPTS)]
    if sorted(P) != list(range(NPTS)):
        probs.append(("C4", "%s: g is not a bijection of Q_6" % t))
    else:
        cyc_sets = set(frozenset(c) for c in cycles_of(P))
        orb_sets = set(frozenset(o) for o in orbs)
        if cyc_sets != orb_sets:
            probs.append(("C4", "%s: orbits_of disagrees with the cycle "
                          "decomposition of the point permutation" % t))
        if sum(len(c) for c in cyc_sets) != NPTS:
            probs.append(("C4", "%s: point-permutation cycles cover %d points"
                          % (t, sum(len(c) for c in cyc_sets))))

    F = build_Fg(a, perm)
    probs.extend(validate_Fg(t, F))
    psi = J_of_graph(F["adj"])
    rec = {"type": t, "alpha": alpha, "beta": beta, "size": class_size(alpha, beta),
           "a": a, "perm": perm, "orbits": orbs, "F": F, "psi": psi,
           "n_kept": F["n"], "n_dropped": len(F["dropped"])}
    return rec, probs


def validate_Fg(t, F):
    """C5: vertex set is exactly the diameter-<=R orbits; adjacency is
    symmetric, loop-free, and reproducible from diam(O union O') > R.  The
    cross-pair reformulation is checked too (it is only valid BECAUSE the
    O_g filter has already been applied)."""
    probs = []
    kept, dropped, adj = F["kept"], F["dropped"], F["adj"]
    for i, o in enumerate(kept):
        d = diameter(o)
        if d > R:
            probs.append(("C5", "%s: kept orbit %r has diameter %d > %d"
                          % (t, o, d, R)))
        if F["kept_diam"][i] != d:
            probs.append(("C5", "%s: cached diameter mismatch on %r" % (t, o)))
    for o in dropped:
        d = diameter(o)
        if d <= R:
            probs.append(("C5", "%s: dropped orbit %r has diameter %d <= %d"
                          % (t, o, d, R)))
    m = len(kept)
    for i in range(m):
        if (adj[i] >> i) & 1:
            probs.append(("C5", "%s: loop at vertex %d" % (t, i)))
        if adj[i] >> m:
            probs.append(("C5", "%s: adjacency bit out of range at %d" % (t, i)))
        for j in range(m):
            if i == j:
                continue
            ij = bool((adj[i] >> j) & 1)
            ji = bool((adj[j] >> i) & 1)
            if ij != ji:
                probs.append(("C5", "%s: asymmetric edge %d-%d" % (t, i, j)))
            want = diameter(kept[i] + kept[j]) > R
            if ij != want:
                probs.append(("C5", "%s: edge %d-%d is %s, union rule says %s"
                              % (t, i, j, ij, want)))
            if want != (cross_max(kept[i], kept[j]) > R):
                probs.append(("C5", "%s: union rule and cross-pair rule "
                              "disagree on %d-%d" % (t, i, j)))
    return probs


# ===========================================================================
# C6 part 1 -- unit tests of the recurrence on hand-computed small graphs,
#              plus a deterministic randomized cross-check of the whole
#              memo/component/branch machinery against both brute forces.
# ===========================================================================
def edge_list_to_adj(n, edges):
    adj = [0] * n
    for (u, v) in edges:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def unit_tests_recurrence():
    """Return (problems, log_lines).  Hand values: J(empty)=1, J(K_1)=0,
    J(K_2)=-1, J(k isolated)=0, J(K_3)=-2, J(P_3)=-1, J(C_4)=-1,
    J(K_{1,3})=-1."""
    probs, log = [], []
    cases = [("empty graph", 0, [], 1),
             ("K_1 (single vertex)", 1, [], 0),
             ("K_2 (single edge)", 2, [(0, 1)], -1),
             ("2 isolated vertices", 2, [], 0),
             ("5 isolated vertices", 5, [], 0),
             ("K_3 (triangle)", 3, [(0, 1), (1, 2), (0, 2)], -2),
             ("P_3 (path)", 3, [(0, 1), (1, 2)], -1),
             ("C_4 (4-cycle)", 4, [(0, 1), (1, 2), (2, 3), (3, 0)], -1),
             ("K_{1,3} (star)", 4, [(0, 1), (0, 2), (0, 3)], -1),
             ("K_2 + K_2 (disjoint)", 4, [(0, 1), (2, 3)], 1)]
    for (name, n, edges, want) in cases:
        adj = edge_list_to_adj(n, edges)
        got = J_of_graph(adj)
        b1 = brute_I_minus1_by_subsets(adj)
        b2 = brute_I_minus1_by_extension(adj)
        log.append("%-22s n=%d  J=%-4d hand=%-4d subsets=%-4d extension=%d"
                   % (name, n, got, want, b1, b2))
        if not (got == want == b1 == b2):
            probs.append(("C6", "%s: J=%d hand=%d subsets=%d extension=%d"
                          % (name, got, want, b1, b2)))
    # deterministic pseudo-random graphs (LCG, no floats, no random module)
    state = 20240819
    for trial in range(40):
        n = 3 + (trial % 8)
        adj = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                state = (1103515245 * state + 12345) % (1 << 31)
                if (state >> 16) & 1:
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
        got, b1, b2 = (J_of_graph(adj), brute_I_minus1_by_subsets(adj),
                       brute_I_minus1_by_extension(adj))
        if not (got == b1 == b2):
            probs.append(("C6", "random graph #%d (n=%d): J=%d subsets=%d "
                          "extension=%d" % (trial, n, got, b1, b2)))
    log.append("40 deterministic pseudo-random graphs (n = 3..10): "
               "recurrence == subset scan == independent-set enumeration")
    return probs, log


# ===========================================================================
# C1, C2
# ===========================================================================
def run_C1():
    """C1: exactly 65 conjugacy classes of H_6, generated from scratch, all
    distinct.  The count 65 is NOT read from the paper."""
    section("C1  class enumeration count")
    classes = all_bipartitions(N)
    pcounts = [len(partitions(k)) for k in range(N + 1)]
    info("p(k) for k=0..6 (derived): %s" % (pcounts,))
    info("sum_k p(k)*p(6-k) = %s = %d"
         % (" + ".join("%d*%d" % (pcounts[k], pcounts[N - k])
                       for k in range(N + 1)),
            sum(pcounts[k] * pcounts[N - k] for k in range(N + 1))))
    types = [type_string(al, be) for (al, be) in classes]
    info("generated %d (alpha;beta) pairs, %d distinct type strings"
         % (len(classes), len(set(types))))

    # Independent test of the partition generator itself (see the docstring of
    # partition_count_pentagonal for why the old "count agrees with
    # sum_k p(k)p(6-k)" check could not fail).
    euler = partition_count_pentagonal(N)
    info("p(k) from Euler's pentagonal recurrence (shares no code): %s"
         % (euler,))
    shape_bad = []
    for k in range(N + 1):
        ps = partitions(k)
        if len(set(ps)) != len(ps):
            shape_bad.append("p(%d): duplicate partitions generated" % k)
        for lam in ps:
            if sum(lam) != k:
                shape_bad.append("p(%d): %r sums to %d" % (k, lam, sum(lam)))
            if list(lam) != sorted(lam, reverse=True):
                shape_bad.append("p(%d): %r not weakly decreasing" % (k, lam))
            if any(part <= 0 for part in lam):
                shape_bad.append("p(%d): %r has a non-positive part" % (k, lam))
    info("generated partitions checked for shape/duplication: %d partitions, "
         "%d violations" % (sum(len(partitions(k)) for k in range(N + 1)),
                            len(shape_bad)))

    ok = True
    ok &= check("C1 generated bipartition count == 65", len(classes) == 65,
                "derived %d" % len(classes))
    ok &= check("C1 all 65 type strings pairwise distinct",
                len(set(types)) == len(types),
                "%d distinct of %d" % (len(set(types)), len(types)))
    ok &= check("C1 every generated pair is a bipartition of 6",
                all(sum(al) + sum(be) == N for (al, be) in classes))
    ok &= check("C1 partition generator agrees with Euler's pentagonal "
                "recurrence for p(0..6), and every generated partition is a "
                "distinct weakly-decreasing partition of its k",
                pcounts == euler and not shape_bad,
                "generator %s, Euler %s, %d shape violations%s"
                % (pcounts, euler, len(shape_bad),
                   ("; first: " + shape_bad[0]) if shape_bad else ""))
    return classes, ok


def run_C2(classes):
    """C2: class sizes positive integers from (eq:class-size), summing to |G|."""
    section("C2  class sizes integral and summing to |G|")
    sizes = {}
    bad = []
    for (al, be) in classes:
        num = (2 ** N) * math.factorial(N)
        den = (2 ** (len(al) + len(be))) * z_lambda(al) * z_lambda(be)
        if num % den != 0:
            bad.append("%s: %d %% %d != 0" % (type_string(al, be), num, den))
            continue
        s = num // den
        if s <= 0:
            bad.append("%s: size %d not positive" % (type_string(al, be), s))
        sizes[type_string(al, be)] = s
    total = sum(sizes.values())
    info("all 65 divisions exact: %s" % (not bad))
    info("sum of derived class sizes = %d" % total)
    info("2^6 * 6! = %d * %d = %d" % (2 ** N, math.factorial(N),
                                      2 ** N * math.factorial(N)))
    info("largest class %d, smallest class %d" % (max(sizes.values()),
                                                 min(sizes.values())))
    ok = True
    ok &= check("C2 every class size is an exact positive integer",
                not bad, "; ".join(bad[:3]))
    ok &= check("C2 sum of class sizes == 46080", total == PAPER_S1,
                "derived %d" % total)
    ok &= check("C2 46080 == 2^6 * 6! (independent of the paper)",
                46080 == 2 ** N * math.factorial(N))
    ok &= check("C2 derived total equals the paper's stated |G|",
                total == PAPER_GROUP_ORDER, "derived %d, paper %d"
                % (total, PAPER_GROUP_ORDER))
    return sizes, ok


# ===========================================================================
# The main derivation: representatives, orbits, F_g, Psi -- with C3/C4/C5
# ===========================================================================
def run_derivation(classes):
    section("D2-D6  representatives, orbits, F_g and Psi = J(F_g) "
            "(with C3, C4, C5)")
    records = []
    probs = []
    for (al, be) in classes:
        rec, p = derive_one(al, be)
        records.append(rec)
        probs.extend(p)
    print("      | %-12s %-6s %-7s %-5s %-5s %-6s %s"
          % ("type", "|C|", "orbits", "|Og|", "drop", "edges", "Psi"))
    for rec in records:
        adj = rec["F"]["adj"]
        nedges = sum(popcount(m) for m in adj) // 2
        print("      | %-12s %-6d %-7d %-5d %-5d %-6d %d"
              % (rec["type"], rec["size"], len(rec["orbits"]), rec["n_kept"],
                 rec["n_dropped"], nedges, rec["psi"]))
    empties = [r["type"] for r in records if r["n_kept"] == 0]
    info("classes with EMPTY F_g (|O_g| = 0): %d -- %s"
         % (len(empties), " ".join(empties)))
    info("their Psi values: %s"
         % sorted(set(r["psi"] for r in records if r["n_kept"] == 0)))
    info("|O_g| ranges from %d to %d over the 65 classes"
         % (min(r["n_kept"] for r in records),
            max(r["n_kept"] for r in records)))
    for cname in ("C3", "C4", "C5"):
        msgs = [m for (c, m) in probs if c == cname]
        label = {"C3": "C3 representative has the claimed type",
                 "C4": "C4 orbit partition valid",
                 "C5": "C5 F_g well formed"}[cname]
        check(label, not msgs,
              ("%d violations, first: %s" % (len(msgs), msgs[0])) if msgs
              else "all 65 classes")
    ok = not probs
    return records, ok


# ===========================================================================
# C6 -- recurrence vs brute force
# ===========================================================================
def run_C6(records, brute_limit=20, subset_limit=14):
    section("C6  J agrees with brute force")
    uprobs, ulog = unit_tests_recurrence()
    for line in ulog:
        info(line)
    check("C6 recurrence unit tests on hand-computed small graphs",
          not uprobs, "; ".join(m for (_, m) in uprobs[:3]))

    mismatches = []
    n_ext = n_sub = 0
    for rec in records:
        if rec["n_kept"] > brute_limit:
            continue
        adj = rec["F"]["adj"]
        b2 = brute_I_minus1_by_extension(adj)
        n_ext += 1
        if b2 != rec["psi"]:
            mismatches.append("%s: recurrence %d, extension %d"
                              % (rec["type"], rec["psi"], b2))
        if rec["n_kept"] <= subset_limit:
            b1 = brute_I_minus1_by_subsets(adj)
            n_sub += 1
            if b1 != rec["psi"]:
                mismatches.append("%s: recurrence %d, subset scan %d"
                                  % (rec["type"], rec["psi"], b1))
    info("classes with |O_g| <= %d re-done by independent-set enumeration: %d"
         % (brute_limit, n_ext))
    info("classes with |O_g| <= %d re-done by full subset scan: %d"
         % (subset_limit, n_sub))
    ok = check("C6 brute force reproduces J on every class with |O_g| <= %d"
               % brute_limit, not mismatches, "; ".join(mismatches[:3]))

    # ---- C6b: a SECOND exact evaluator on the classes brute force cannot
    # reach.  Enumeration is infeasible above ~20 vertices (F_g for 2211|- is
    # 40 vertices / 124 edges), which left 16 classes -- every large Psi in the
    # paper, including Psi(id) = 253 -- resting on a single algorithm.  The
    # alternate evaluator uses the same recurrence but a MINIMUM-degree branch
    # vertex and NO component multiplicativity, so a bug in `components` or in
    # the branch selection cannot be shared.
    # The budget is PER CLASS.  An induced-subgraph memo on v vertices can never
    # hold more than 2^v masks, so a budget of 2.5e6 GUARANTEES completion for
    # every class with |O_g| <= 21 -- i.e. for at least the same 49 classes that
    # enumerative brute force reaches, which is what the hard check below
    # requires.  Classes above that are best-effort and are NAMED either way, so
    # the coverage of this check can never be silently zero.
    alt_budget = 2500000
    alt_bad, alt_done, alt_skipped = [], [], []
    for rec in sorted(records, key=lambda r: r["n_kept"]):
        val, done, _nodes = J_of_graph_alt(rec["F"]["adj"], alt_budget)
        if not done:
            alt_skipped.append(rec["type"])
            continue
        alt_done.append(rec["type"])
        if val != rec["psi"]:
            alt_bad.append("%s: primary %d, min-degree/no-component %d"
                           % (rec["type"], rec["psi"], val))
    info("classes re-evaluated by the min-degree / no-component-split "
         "recurrence: %d of %d (per-class node budget %d)"
         % (len(alt_done), len(records), alt_budget))
    if alt_skipped:
        info("classes NOT cross-checked by the second evaluator (per-class "
             "node budget exhausted): %s" % " ".join(alt_skipped))
    big = [r["type"] for r in records if r["n_kept"] > brute_limit]
    big_done = [t for t in big if t in set(alt_done)]
    info("of the %d classes with |O_g| > %d (which enumerative brute force "
         "cannot reach), %d were cross-checked by the second evaluator: %s"
         % (len(big), brute_limit, len(big_done),
            " ".join(big_done) if big_done else "NONE"))
    small_n = sum(1 for r in records if r["n_kept"] <= brute_limit)
    ok &= check("C6b second exact evaluator (min-degree branch, no component "
                "multiplicativity) agrees with J on every class it reached, "
                "and reached at least the %d classes brute force covers"
                % small_n,
                not alt_bad and len(alt_done) >= small_n,
                ("%d disagreements, first: %s" % (len(alt_bad), alt_bad[0]))
                if alt_bad else "%d of %d classes agreed"
                % (len(alt_done), len(records)))

    # A third, fully explicit view for a few small classes: the independent-set
    # polynomial I(F_g;t) by coefficient, whose alternating sum must be Psi.
    poly_bad = []
    shown = 0
    for rec in sorted(records, key=lambda r: r["n_kept"]):
        if rec["n_kept"] == 0 or rec["n_kept"] > 12:
            continue
        coeffs = independent_set_polynomial(rec["F"]["adj"])
        alt = sum((-1) ** k * c for k, c in enumerate(coeffs))
        if alt != rec["psi"]:
            poly_bad.append("%s: alt sum %d vs Psi %d"
                            % (rec["type"], alt, rec["psi"]))
        if shown < 6:
            info("I(F;t) for %-12s |O_g|=%-3d coeffs %s -> I(F;-1) = %d"
                 % (rec["type"], rec["n_kept"], coeffs, alt))
            shown += 1
    ok &= check("C6 independent-set polynomial's alternating sum equals Psi "
                "on every class with 1 <= |O_g| <= 12", not poly_bad,
                "; ".join(poly_bad[:3]))
    return ok and not uprobs


# ===========================================================================
# C7 -- the appendix table is the CLAIM, joined to derived values by type
# ===========================================================================
def run_C7(records, witness, dups, order):
    section("C7  appendix table agreement (table = claim, not input)")
    derived = dict((r["type"], (r["size"], r["psi"])) for r in records)
    info("witness rows parsed: %d (duplicates: %s)"
         % (len(order), dups if dups else "none"))
    malformed = [t for t in witness if parse_type_string(t) is None]
    info("witness rows that are well-formed bipartitions of 6: %d of %d"
         % (len(witness) - len(malformed), len(witness)))
    missing = sorted(set(derived) - set(witness))
    extra = sorted(set(witness) - set(derived))
    csize_bad = [(t, witness[t][0], derived[t][0]) for t in sorted(derived)
                 if t in witness and witness[t][0] != derived[t][0]]
    psi_bad = [(t, witness[t][1], derived[t][1]) for t in sorted(derived)
               if t in witness and witness[t][1] != derived[t][1]]
    matched = [t for t in derived if t in witness]
    info("types matched on both sides: %d ; |C| agreements: %d ; "
         "Psi agreements: %d" % (len(matched), len(matched) - len(csize_bad),
                                 len(matched) - len(psi_bad)))
    for (t, w, d) in csize_bad:
        info("  |C| MISMATCH %-12s paper %d, derived %d" % (t, w, d))
    for (t, w, d) in psi_bad:
        info("  Psi MISMATCH %-12s paper %d, derived %d" % (t, w, d))
    ok = True
    ok &= check("C7 witness table has 65 rows, no duplicates, all well formed",
                len(order) == 65 and not dups and not malformed,
                "rows=%d dups=%s malformed=%s" % (len(order), dups, malformed))
    ok &= check("C7 no derived class missing from the table", not missing,
                " ".join(missing[:5]))
    ok &= check("C7 no table row without a derived class", not extra,
                " ".join(extra[:5]))
    ok &= check("C7 all 65 |C| values agree", not csize_bad and
                len(matched) == 65, "%d mismatches" % len(csize_bad))
    ok &= check("C7 all 65 Psi values agree", not psi_bad and
                len(matched) == 65, "%d mismatches" % len(psi_bad))
    return ok


# ===========================================================================
# C8 -- the two sums, from DERIVED data only
# ===========================================================================
def run_C8(records):
    section("C8  the two sums")
    S1 = sum(r["size"] for r in records)
    S2 = sum(r["size"] * r["psi"] for r in records)
    info("S1 = sum |C|            = %d   (derived)" % S1)
    info("S2 = sum |C| * Psi      = %d   (derived)" % S2)
    pos = sum(r["size"] * r["psi"] for r in records if r["psi"] > 0)
    neg = sum(r["size"] * r["psi"] for r in records if r["psi"] < 0)
    info("positive contribution %d, negative contribution %d, net %d"
         % (pos, neg, pos + neg))
    info("classes with Psi < 0: %s"
         % " ".join("%s(%d)" % (r["type"], r["psi"])
                    for r in records if r["psi"] < 0))
    idc = [r for r in records if r["type"] == "111111|-"]
    if idc:
        info("identity class contributes |C|*Psi = %d * %d = %d of S2"
             % (idc[0]["size"], idc[0]["psi"], idc[0]["size"] * idc[0]["psi"]))
    info("distinct Psi values over the 65 classes: %s"
         % sorted(set(r["psi"] for r in records)))
    q, rem = divmod(S2, S1)
    info("S2 = %d * S1 + %d  (exact integer division, no floats)" % (q, rem))
    ok = True
    ok &= check("C8 S1 == 46080", S1 == PAPER_S1, "derived %d" % S1)
    ok &= check("C8 S2 == 138240", S2 == PAPER_S2, "derived %d" % S2)
    ok &= check("C8 S2 == 3 * S1 exactly (integer decision)",
                S2 == 3 * S1 and rem == 0 and q == PAPER_AVERAGE,
                "quotient %d remainder %d" % (q, rem))
    return S1, S2, ok


# ===========================================================================
# C9 -- identity class cross-check against (eq:homology)
# ===========================================================================
def run_C9(records):
    section("C9  identity-class cross-check against (eq:homology)")
    idt = "111111|-"
    rec = [r for r in records if r["type"] == idt]
    ok = check("C9 the identity class 111111|- was generated", len(rec) == 1)
    if not ok:
        return False
    rec = rec[0]
    a, perm = rec["a"], rec["perm"]
    info("identity representative: a = %s (bits %s), pi = %s"
         % (a, format(a, "06b"), perm))
    info("orbits: %d, all singletons: %s"
         % (len(rec["orbits"]), all(len(o) == 1 for o in rec["orbits"])))
    adj = rec["F"]["adj"]
    nedges = sum(popcount(m) for m in adj) // 2
    degs = set(popcount(m) for m in adj)
    far = sum(1 for x in range(NPTS) for y in range(x + 1, NPTS)
              if hamming(x, y) >= R + 1)
    info("|O_g| = %d vertices, %d edges, degree set %s"
         % (rec["n_kept"], nedges, sorted(degs)))
    info("pairs {x,y} in Q_6 with d(x,y) >= 5 counted directly: %d" % far)
    ok &= check("C9 identity F_g is the 64-vertex d>=5 graph, 7-regular",
                rec["n_kept"] == NPTS and nedges == far and degs == {7},
                "n=%d edges=%d far=%d degs=%s"
                % (rec["n_kept"], nedges, far, sorted(degs)))
    ok &= check("C9 derived Psi(111111|-) == +253 (sign convention "
                "(-1)^{i+1}, so Psi(id) = -chi~(K) = +253)",
                rec["psi"] == 253, "derived %d" % rec["psi"])
    info("paper's (eq:homology), cited from GMW C.3: dim H_7 = %d, "
         "dim H_15 = %d, sum = %d"
         % (PAPER_H7_DIM, PAPER_H15_DIM, PAPER_H7_DIM + PAPER_H15_DIM))

    # Three former checks lived here:
    #     "253 == 239 + 14"                 -- a literal against two constants
    #                                          in this file; could not fail;
    #     "derived Psi(id) == 239 + 14"      -- the SAME condition as
    #                                          rec["psi"] == 253 above;
    #     "the sign is + not -"              -- rec["psi"] > 0 and != -253,
    #                                          implied by rec["psi"] == 253,
    #                                          with the != -253 half dead.
    # They are replaced by ONE check that actually has content: reconstruct the
    # reduced Euler characteristic BOTH ways and compare.  Derived side:
    # Psi = sum_i (-1)^{i+1} tr, so at g = id, chi~(K) = -Psi(id).  Paper side:
    # chi~(K) = sum_i (-1)^i dim H~_i, evaluated on (eq:homology) using the
    # PAPER'S DEGREES with their own signs -- so this bites on the PARITY of
    # the two degrees as well as on the two dimensions.  With degrees (7,15)
    # both odd it must equal -(239 + 14); had GMW reported a nonzero group in
    # an even degree, the two sides would differ and this check would fail.
    dims = (PAPER_H7_DIM, PAPER_H15_DIM)
    chi_paper = sum((-1) ** d * m
                    for d, m in zip(PAPER_HOMOLOGY_DEGREES, dims))
    chi_derived = -rec["psi"]
    info("chi~(K) from the derived Psi(id) (= -Psi(id)) = %d" % chi_derived)
    info("chi~(K) from (eq:homology) as sum_i (-1)^i dim H~_i over degrees "
         "%s = %s = %d"
         % (PAPER_HOMOLOGY_DEGREES,
            " + ".join("(-1)^%d*%d" % (d, m)
                       for d, m in zip(PAPER_HOMOLOGY_DEGREES, dims)),
            chi_paper))
    ok &= check("C9 chi~(K) from the derived Psi equals sum_i (-1)^i dim H~_i "
                "from (eq:homology), and both cited degrees are odd so the "
                "signs are +1 in Psi",
                chi_derived == chi_paper
                and chi_paper == -(PAPER_H7_DIM + PAPER_H15_DIM),
                "derived %d, paper %d, -(239+14) = %d"
                % (chi_derived, chi_paper, -(PAPER_H7_DIM + PAPER_H15_DIM)))
    return ok


# ===========================================================================
# C10 -- what the average 3 does and does not license
# ===========================================================================
def run_C10(S1, S2):
    section("C10  the conclusion is a valid inference")
    avg = S2 // S1
    info("average (1/|G|) sum_g Psi(g) = %d, decided by %d == %d * %d"
         % (avg, S2, avg, S1))
    info("character averaging (P8) plus reduced homology concentrated in "
         "degrees 7 and 15 -- BOTH ODD, so (-1)^{i+1} = +1 for both -- gives")
    info("    d7 + d15 = %d   where d_i = dim H_i(K;Q)^G" % avg)
    splits = [(d7, avg - d7) for d7 in range(avg + 1)]
    info("all (d7,d15) with d7 + d15 = %d and both >= 0: %s" % (avg, splits))
    all_have_big = all(max(s) >= 2 for s in splits)
    info("every split has a coordinate >= 2: %s  ==>  at least one of H_7, "
         "H_15 is NOT multiplicity-free" % all_have_big)
    info("number of surviving splits = %d > 1  ==>  WHICH degree carries the "
         "repeat is UNDETERMINED by this computation" % len(splits))
    info("the average is %s, which is why no split can have both parts <= 1"
         % ("odd" if avg % 2 else "even"))
    ok = True
    # Two former checks here were the same predicate written twice
    # ("no split has both parts <= 1" and "every split has a part >= 2"), and
    # the first one additionally demanded avg % 2 == 1, which the CONCLUSION
    # does not need -- an average of 4 would still force a part >= 2 yet would
    # have FAILED that check.  Merged into the one condition the inference
    # actually uses.  The "both degrees are odd" check that also lived here
    # tested the literal tuple (7,15) against itself; it now has real content
    # in C9, where it is compared against the derived Euler characteristic.
    ok &= check("C10 every nonnegative split of the derived average has a part "
                ">= 2, so at least one of H_7, H_15 is not multiplicity-free",
                all_have_big and avg > 0,
                "avg = %d, splits %s" % (avg, splits))
    ok &= check("C10 more than one split survives, so the degree is not "
                "determined (this program does NOT report which)",
                len(splits) > 1, "%d splits" % len(splits))
    return ok


# ===========================================================================
# C11 -- minimality of the parameter point (6,4)
# ===========================================================================
def run_C11():
    section("C11  minimality of the parameter point")
    pairs = sorted((n, r) for n in range(1, 9) for r in range(0, n + 1)
                   if 4 <= r <= n - 2)
    info("all (n,r) with n <= 8 and 4 <= r <= n-2: %s" % (pairs,))
    info("r >= 4 and r <= n-2 force n >= 6; at n = 6 the only r is 4")
    # "no admissible pair has n < 6" used to be a separate check, but `pairs` is
    # sorted and the check below already pins pairs[0] == (6,4), so it could not
    # fail independently.  Folded in as an explicit conjunct instead.
    ok = check("C11 lexicographic minimum of {(n,r) : 4 <= r <= n-2} is (6,4) "
               "and no admissible pair has n < 6",
               bool(pairs) and pairs[0] == (N, R)
               and not any(n < 6 for (n, r) in pairs),
               "min %s, %d admissible pairs with n <= 8"
               % (pairs[0] if pairs else None, len(pairs)))
    return ok


# ===========================================================================
# C12 -- the bridge from "65 representatives" to the paper's "sum over g in G"
# ===========================================================================
def run_C12(records, sample_per_type=3):
    """Everything above evaluates Psi on ONE element per type.  The paper's
    Proposition is a sum over all 46080 ELEMENTS of G.  Two links were
    previously unasserted anywhere in this file:

      (a) that the (eq:class-size) VALUE equals the ACTUAL number of elements
          of that type.  C2 only checks that the 65 formula values SUM to
          46080, and a sum can match while individual terms are wrong in
          compensating directions -- so a paper error in (eq:class-size) of
          that shape would have survived every other check here;
      (b) that Psi = J(F_g) is CONSTANT on a type block, i.e. that "cycle
          lengths plus cycle parities" really is a complete conjugacy
          invariant (the paper cites GMW 5.2 for this).  If a block were a
          union of two genuine classes with different J, the weighted sum
          would be wrong and nothing else here would notice.

    Both are cheap to settle by brute force, so we do: enumerate every one of
    the 6! * 2^6 = 46080 pairs (a, pi), classify each by its OWN recovered
    (alpha;beta) invariant, and compare the census with the derived data class
    by class.  Then re-derive F_g and J from scratch on extra, non-representative
    elements of each block."""
    import itertools

    section("C12  full census of all 46080 elements of G "
            "(bridge: 65 classes -> sum over g in G)")
    derived_size = dict((r["type"], r["size"]) for r in records)
    derived_psi = dict((r["type"], r["psi"]) for r in records)
    rep_of = dict((r["type"], (r["a"], r["perm"])) for r in records)

    census = {}
    samples = {}
    total_elts = 0
    bad_total = []
    for perm in itertools.permutations(range(N)):
        for a in range(NPTS):
            ral, rbe, tot = recover_type(a, perm)
            if tot != N:
                bad_total.append((a, perm, tot))
            t = type_string(ral, rbe)
            census[t] = census.get(t, 0) + 1
            total_elts += 1
            if rep_of.get(t) != (a, perm):
                lst = samples.setdefault(t, [])
                if len(lst) < sample_per_type:
                    lst.append((a, perm))

    info("elements enumerated: %d ; distinct types found by census: %d"
         % (total_elts, len(census)))
    info("census total vs 2^6 * 6! : %d vs %d"
         % (total_elts, 2 ** N * math.factorial(N)))
    size_bad = sorted(t for t in set(census) | set(derived_size)
                      if census.get(t) != derived_size.get(t))
    for t in size_bad[:5]:
        info("  |C| CENSUS MISMATCH %-12s census %s, (eq:class-size) %s"
             % (t, census.get(t), derived_size.get(t)))
    info("types whose census count equals the (eq:class-size) value: %d of %d"
         % (len(census) - len(size_bad), len(census)))

    # The Proposition as literally written: a sum over g in G, with no appeal
    # to the class-size formula at all.
    S2_census = sum(census[t] * derived_psi[t] for t in census
                    if t in derived_psi)
    unmatched = sorted(t for t in census if t not in derived_psi)
    info("sum_{g in G} Psi(g) recomputed element-by-element from the census "
         "= %d" % S2_census)
    if unmatched:
        info("census types with NO derived Psi (would silently drop from the "
             "sum): %s" % " ".join(unmatched))

    # (b): Psi must not depend on which element of the block we picked.
    inv_bad = []
    n_probe = 0
    probed_types = []
    for t in sorted(samples):
        for (a, perm) in samples[t]:
            F = build_Fg(a, perm)
            val = J_of_graph(F["adj"])
            n_probe += 1
            if val != derived_psi[t]:
                inv_bad.append("%s: representative Psi %d, element "
                               "(a=%d, pi=%r) gives %d"
                               % (t, derived_psi[t], a, perm, val))
            if F["n"] != next(r["n_kept"] for r in records if r["type"] == t):
                inv_bad.append("%s: |O_g| differs on a second element of the "
                               "same block" % t)
        probed_types.append(t)
    no_sample = sorted(t for t in derived_psi if t not in samples)
    info("blocks probed with up to %d EXTRA (non-representative) elements: "
         "%d of 65, %d elements re-derived from scratch"
         % (sample_per_type, len(probed_types), n_probe))
    info("blocks with no second element to probe (block size 1): %s"
         % (" ".join(no_sample) if no_sample else "none"))

    ok = True
    ok &= check("C12 census enumerated exactly 2^6*6! = 46080 elements and "
                "every one has cycle lengths summing to 6",
                total_elts == 2 ** N * math.factorial(N) and not bad_total,
                "%d elements, %d bad" % (total_elts, len(bad_total)))
    ok &= check("C12 the census finds exactly 65 types, and they are the 65 "
                "generated classes", len(census) == 65
                and set(census) == set(derived_size),
                "census %d types, derived %d types"
                % (len(census), len(derived_size)))
    ok &= check("C12 (eq:class-size) is correct CLASS BY CLASS, not merely in "
                "total: census count == |C| for all 65", not size_bad,
                "%d mismatches: %s" % (len(size_bad), " ".join(size_bad[:5])))
    ok &= check("C12 sum_{g in G} Psi(g) == 138240 computed over ELEMENTS "
                "(the Proposition as written, no class-size formula used)",
                S2_census == PAPER_S2 and not unmatched,
                "census sum %d, %d unmatched types" % (S2_census, len(unmatched)))
    ok &= check("C12 Psi is constant on each type block: J(F_g) re-derived on "
                "%d extra elements agrees with the representative" % n_probe,
                not inv_bad and n_probe > 0,
                ("%d disagreements, first: %s" % (len(inv_bad), inv_bad[0]))
                if inv_bad else "%d elements over %d blocks"
                % (n_probe, len(probed_types)))
    return ok


# ===========================================================================
# Limitations -- printed every run, pass or fail.  This file does NOT verify
# the theorem end to end and must not claim to.
# ===========================================================================
def print_limitations():
    section("LIMITATIONS  (what this file does NOT establish)")
    for line in [
        "1. Lemma lem:fixed-orbit, Psi(g) = I(F_g;-1), is ASSUMED. This program",
        "   computes only the right-hand side; it never builds a chain complex",
        "   and never computes homology.",
        "2. The individual dimensions 239 (degree 7) and 14 (degree 15) are",
        "   NOT derived here. They are cited by the paper from GMW Appendix",
        "   C.3. C9 tests only their SUM, 239 + 14 = 253, against the derived",
        "   Psi at the identity class.",
        "3. The claim that reduced H_i(X^{6,4};Q) VANISHES outside {7,15} is",
        "   likewise external and unverified here. The averaging step",
        "   (1/|G|) sum_g Psi(g) = d7 + d15 needs BOTH that vanishing AND the",
        "   fact that 7 and 15 are odd (so both signs are +1). A nonzero",
        "   EVEN-degree group in GMW's computation would enter with the",
        "   opposite sign and break the identification -- this file cannot",
        "   close that gap. Doing so would mean building the rational chain",
        "   complex of VR(Q_6;4) with faces up to dimension 15 on 64 vertices.",
        "4. Reduced vs unreduced homology: the theorem is stated with",
        "   unreduced H_7, H_15 while the computation is with reduced homology.",
        "   They agree in positive degrees, and 7, 15 > 0, so the step is",
        "   sound -- but it is a step, not an identity of definitions.",
        "5. WHICH of degree 7 or 15 carries the repeated trivial constituent",
        "   is NOT determined (C10). Only 'at least one of H_7, H_15 is not",
        "   multiplicity-free' follows from d7 + d15 = 3.",
        "6. Convention: the action is taken as (pi x)_i = x_{pi^{-1}(i)}. The",
        "   other convention yields a conjugate of g^{-1}; orbits of g and of",
        "   g^{-1} coincide and conjugation by a coordinate permutation is a",
        "   Hamming isometry, so every quantity above is unchanged.",
        "7. The theorem is stated over EVERY field of characteristic zero. The",
        "   base-change step -- H_i(K;k) = H_i(K;Q) tensor k, and invariants",
        "   commuting with it via the Reynolds idempotent |G|^{-1} sum_g g --",
        "   is prose in the paper and is NOT modelled here. Everything this",
        "   file computes is over Z and concerns Q-coefficients only. The step",
        "   is standard and correct, but it is assumed, not verified.",
        "8. C6's enumerative brute force cannot reach the classes with",
        "   |O_g| > 20. Those are covered instead by C6b, a second exact",
        "   evaluator using a different branch rule and no component",
        "   multiplicativity, under a node budget; any class the budget did",
        "   not reach is named in the C6b output rather than silently passing.",
    ]:
        print("      | " + line)


def main():
    # J_value recurses at most once per deleted vertex (<= 64 deep) and the
    # brute forces at most |O_g|+1 deep; this is headroom, not a requirement.
    sys.setrecursionlimit(20000)
    print("verify_vr_q6_r4.py -- K = X^{6,4} = VR(Q_6;4), G = H_6 = "
          "F_2^6 : S_6")
    print("Checking: Psi(g) = I(F_g;-1) on all conjugacy classes; "
          "sum |C| = 46080; sum |C| Psi = 138240; average 3.")
    print("FROM THE PAPER: n=6, r=4; the group and its action; Lemma "
          "lem:fixed-orbit (assumed);")
    print("                (eq:recurrence); (eq:class-size); the "
          "representative recipe;")
    print("                (eq:homology) 239 + 14 from GMW C.3; character "
          "averaging; and the")
    print("                65-row Appendix A table, held ONLY as the claim "
          "under test.")
    print("DERIVED HERE:   the 65 bipartitions, representatives, class sizes, "
          "orbits, diameters,")
    print("                O_g, every edge of F_g, all 65 Psi values, both "
          "sums, and Psi(id).")

    classes, ok1 = run_C1()
    sizes, ok2 = run_C2(classes)
    records, ok345 = run_derivation(classes)
    ok6 = run_C6(records)
    witness, dups, order = parse_witness()
    ok7 = run_C7(records, witness, dups, order)
    S1, S2, ok8 = run_C8(records)
    ok9 = run_C9(records)
    ok10 = run_C10(S1, S2)
    ok11 = run_C11()
    ok12 = run_C12(records)

    # consistency of the two independent class-size computations
    section("cross-check: class sizes from run_C2 vs from derive_one")
    bad = [r["type"] for r in records if sizes.get(r["type"]) != r["size"]]
    check("class sizes agree between the two computation sites", not bad,
          " ".join(bad[:5]))

    # The verdict is computed from the RESULTS log; this confirms the log did
    # not lose a section's outcome.
    section("bookkeeping: section return values vs the recorded check log")
    section_ok = all([ok1, ok2, ok345, ok6, ok7, ok8, ok9, ok10, ok11, ok12])
    log_ok = all(okk for (_, okk, _) in RESULTS)
    info("AND of section return values: %s ; AND of the %d logged checks: %s"
         % (section_ok, len(RESULTS), log_ok))
    # NOTE: this is BOOKKEEPING, not mathematics.  It compares two summaries of
    # the same run and passes whenever they agree -- including when both are
    # False.  It is named accordingly so it cannot be mistaken for a
    # mathematical check in the final count.
    check("bookkeeping: section return values agree with the recorded check log",
          section_ok == log_ok)

    print_limitations()
    return verdict()


def verdict():
    section("VERDICT")
    total = len(RESULTS)
    failed = [nm for (nm, okk, _) in RESULTS if not okk]
    for nm in failed:
        print("      | FAILED: " + nm)
    if failed:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:                      # noqa: BLE001
        import traceback
        traceback.print_exc()
        # An exception is a FAILED verification, not an absent one: record it
        # so the transcript still ends in the contracted verdict line.
        check("no unhandled exception during verification", False,
              "%s: %s" % (type(exc).__name__, exc))
        sys.exit(verdict())

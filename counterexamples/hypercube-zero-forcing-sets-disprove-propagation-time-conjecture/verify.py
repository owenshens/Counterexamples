#!/usr/bin/env python3
"""Independent verification of two explicit minimum zero forcing sets of
hypercubes whose propagation times exceed a conjectured maximum.

Setting.  Vertices of Q_d are the integers 0..2^d-1; x and y are adjacent
exactly when x xor y is a power of two.  For a filled set C, a filled vertex
u forces v when N(u) \\ C = {v}.  With C_0 = S, D_t = Phi(C_{t-1}) and
C_t = C_{t-1} union D_t; if C_T = V and C_{T-1} != V then pt(Q_d,S) = T.  PT(G)
is the maximum of pt(G,S) over all MINIMUM zero forcing sets S.

TAKEN FROM THE PAPER (inputs, transcribed verbatim, nothing else):
  * the two exhibited vertex sets S_6 (subset of Q_6) and S_7 (subset of Q_7);
  * the conjectured formula PT(Q_d) = 2^(d-2) that is being refuted;
  * the asserted propagation times 18 and 43, and the asserted layer
    sequences (D_1,...,D_18) and (D_1,...,D_43), used only as data to be
    compared against independently recomputed layers;
  * the cited identity Z(Q_d) = 2^(d-1) for the zero forcing number, used
    only where explicitly named.

DERIVED HERE (all decisions computed from scratch with exact integers, no
floating point):
  * well-formedness and cardinality of the two exhibited sets;
  * that each set is a zero forcing set of the correct hypercube, by running
    the forcing process to closure;
  * the propagation times, computed layer by layer, and the strict
    inequalities against the conjectured formula (the refutation);
  * that every recomputed layer equals the layer printed in the paper, that
    the layers partition the complement of S, and that each forced vertex has
    a private forcing witness;
  * inclusion-wise minimality of both sets (no single-vertex deletion forces);
  * Z(Q_d) = 2^(d-1) exhaustively for d <= 4, the upper bound Z(Q_d) <=
    2^(d-1) constructively for 2 <= d <= 7, and PT(Q_d) = 2^(d-2)
    exhaustively over all minimum zero forcing sets for d <= 4;
  * a deterministic (randomness-free) ascent over Q_5 confirming that 2^(5-2)
    is attained there and is not exceeded by anything the ascent reaches.

NOT RE-RUN (printed at the end of the run as well): the lower bound
Z(Q_d) >= 2^(d-1) for d = 5,6,7, which the paper takes from the literature;
and the exhaustive census of all minimum zero forcing sets of Q_5, Q_6, Q_7,
which is far beyond a single-process budget.
"""
import sys
from itertools import combinations

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
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


# ---------------------------------------------------------------- paper data
S6_PAPER = [2, 5, 6, 7, 10, 13, 16, 17, 18, 20, 21, 22, 25, 28, 29, 30, 31,
            33, 35, 36, 38, 43, 44, 49, 50, 51, 56, 57, 58, 60, 62, 63]

S7_PAPER = [0, 1, 2, 3, 4, 10, 11, 13, 14, 15, 21, 22, 24, 27, 30, 31,
            34, 36, 37, 38, 39, 41, 45, 46, 48, 49, 50, 51, 55, 56, 57, 58,
            63, 64, 65, 66, 67, 68, 69, 71, 74, 77, 79, 85, 88, 90, 93, 98,
            100, 101, 102, 103, 108, 109, 112, 113, 115, 119, 120, 122, 123,
            124, 126, 127]

# The layer sequence (D_1,...,D_18) as printed in the paper.
LAYERS6_PAPER = [[61], [52], [4], [14], [26], [19],
                 [1], [37], [32], [41], [59],
                 [27, 42, 55], [11, 24, 47], [8, 9, 12], [0],
                 [3, 40, 48], [15, 23, 34, 45, 53], [39, 46, 54]]

# The layer sequence (D_1,...,D_43) as printed in the paper.
LAYERS7_PAPER = [[5], [7], [47], [35], [19], [59], [43], [111], [99], [33],
                 [53], [17], [9], [75], [83], [114], [106], [96], [116], [80],
                 [72], [78], [95], [125], [105], [40], [60], [92], [89], [29],
                 [12, 23], [6, 18, 20, 54, 87], [70, 86], [82], [118],
                 [110, 117],
                 [97, 107], [32, 42], [104], [121],
                 [61, 73, 81, 91], [25, 44, 52, 62, 76, 84, 94],
                 [8, 16, 26, 28]]

PT6_PAPER = 18          # asserted pt(Q_6, S_6)
PT7_PAPER = 43          # asserted pt(Q_7, S_7)


def conjectured_PT(d):
    """The conjectured formula PT(Q_d) = 2^(d-2), taken from the paper."""
    return 2 ** (d - 2)


def cited_Z(d):
    """The cited identity Z(Q_d) = 2^(d-1)."""
    return 2 ** (d - 1)


# ------------------------------------------------------- graph and mechanics
def hypercube(d):
    """Neighbourhood bitmasks of Q_d: N(u) = {u xor 2^i : 0 <= i < d}."""
    return [sum(1 << (u ^ (1 << i)) for i in range(d)) for u in range(1 << d)]


def to_mask(vertices):
    m = 0
    for v in vertices:
        m |= 1 << v
    return m


def to_list(mask):
    out = []
    while mask:
        b = mask & -mask
        out.append(b.bit_length() - 1)
        mask ^= b
    return out


def popcount(mask):
    return bin(mask).count("1")


def forcing_layers(nbr, full, start):
    """Run the synchronous forcing process.  Returns (layers, closure) where
    layers is the list of masks D_1, D_2, ... and closure is the final filled
    set (equal to full exactly when start is a zero forcing set)."""
    filled = start
    layers = []
    while filled != full:
        newly = 0
        rest = filled
        while rest:
            b = rest & -rest
            u = b.bit_length() - 1
            rest ^= b
            free = nbr[u] & ~filled
            if free and (free & (free - 1)) == 0:
                newly |= free
        if newly == 0:
            break
        layers.append(newly)
        filled |= newly
    return layers, filled


def is_forcing(nbr, full, start):
    return forcing_layers(nbr, full, start)[1] == full


def nbr_from_edges(n, edges):
    adj = [0] * n
    for a, b in edges:
        adj[a] |= 1 << b
        adj[b] |= 1 << a
    return adj


# ------------------------------------------------------------------ check 1/2
def check_wellformed(tag, d, listed):
    """Decode the exhibited set, count it, print it back."""
    n = 1 << d
    ok_range = all(isinstance(v, int) and 0 <= v < n for v in listed)
    ok_distinct = len(set(listed)) == len(listed)
    card = len(set(listed))
    ok_card = card == cited_Z(d)
    dec = sorted(set(listed))
    print("info: S_%d decoded (%d vertices of Q_%d): %s"
          % (d, card, d, ",".join(str(v) for v in dec)))
    print("info: S_%d as binary words: %s"
          % (d, " ".join(format(v, "0" + str(d) + "b") for v in dec[:8])
             + (" ..." if card > 8 else "")))
    return ck("%s_wellformed" % tag,
              ok_range and ok_distinct and ok_card,
              "%d distinct vertices in [0,%d), 2^(d-1)=%d"
              % (card, n, cited_Z(d)))


# -------------------------------------------------------------------- check 3
def check_hypercube_structure():
    """The neighbourhood table really is Q_d: d-regular, d*2^(d-1) edges,
    matches the independent xor-popcount definition, connected, bipartite."""
    bad = []
    for d in range(2, 8):
        n = 1 << d
        nbr = hypercube(d)
        if any(popcount(nbr[u]) != d for u in range(n)):
            bad.append("degree d=%d" % d)
        edges = sum(popcount(nbr[u]) for u in range(n)) // 2
        if edges != d * (1 << (d - 1)):
            bad.append("edges d=%d" % d)
        for x in range(n):
            ref = to_mask(y for y in range(n)
                          if popcount(x ^ y) == 1)
            if ref != nbr[x]:
                bad.append("adjacency d=%d" % d)
                break
        seen, stack = 1, [0]
        while stack:
            u = stack.pop()
            for v in to_list(nbr[u] & ~seen):
                seen |= 1 << v
                stack.append(v)
        if popcount(seen) != n:
            bad.append("connected d=%d" % d)
        if any(popcount(x ^ y) % 2 == 0
               for x in range(n) for y in to_list(nbr[x])):
            bad.append("bipartite d=%d" % d)
    return ck("hypercube_structure", not bad,
              "d=2..7 regular, connected, bipartite"
              if not bad else "; ".join(sorted(set(bad))))


# -------------------------------------------------------------------- check 4
def check_engine_sanity():
    """The forcing engine reproduces textbook values on graphs whose zero
    forcing sets and propagation times are known by hand."""
    cases = []
    # path on 5 vertices: one endpoint forces in 4 rounds, the centre never
    p5 = nbr_from_edges(5, [(0, 1), (1, 2), (2, 3), (3, 4)])
    full5 = (1 << 5) - 1
    cases.append(("P5 endpoint pt=4",
                  forcing_layers(p5, full5, to_mask([0])) == (
                      [to_mask([1]), to_mask([2]), to_mask([3]),
                       to_mask([4])], full5)))
    cases.append(("P5 centre not forcing",
                  not is_forcing(p5, full5, to_mask([2]))))
    # 5-cycle: an adjacent pair forces in 2 rounds, a non-adjacent pair fails
    c5 = nbr_from_edges(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
    lay, clo = forcing_layers(c5, full5, to_mask([0, 1]))
    cases.append(("C5 adjacent pair pt=2", clo == full5 and len(lay) == 2))
    cases.append(("C5 spread pair not forcing",
                  not is_forcing(c5, full5, to_mask([0, 2]))))
    # complete graph on 4 vertices: Z = 3 and pt = 1
    k4 = nbr_from_edges(4, [(a, b) for a in range(4) for b in range(a + 1, 4)])
    full4 = (1 << 4) - 1
    lay, clo = forcing_layers(k4, full4, to_mask([0, 1, 2]))
    cases.append(("K4 pt=1", clo == full4 and len(lay) == 1))
    cases.append(("K4 two vertices not forcing",
                  not is_forcing(k4, full4, to_mask([0, 1]))))
    # empty start forces nothing
    cases.append(("empty set forces nothing",
                  forcing_layers(p5, full5, 0) == ([], 0)))
    bad = [name for name, ok in cases if not ok]
    return ck("forcing_engine_sanity", not bad,
              "%d hand-checked cases" % len(cases)
              if not bad else "; ".join(bad))


# -------------------------------------------------------------------- check 5
def check_is_zero_forcing(tag, d, listed):
    """Hypothesis of the theorem: the exhibited set forces all of Q_d."""
    n = 1 << d
    full = (1 << n) - 1
    nbr = hypercube(d)
    lay, clo = forcing_layers(nbr, full, to_mask(listed))
    missed = popcount(full) - popcount(clo)
    return ck("%s_is_zero_forcing" % tag, clo == full,
              "closure = all %d vertices in %d rounds" % (n, len(lay))
              if clo == full else "%d vertices never forced" % missed)


# ------------------------------------------------------------------ check 6/7
def check_propagation_time(tag, d, listed, claimed):
    """Compute pt(Q_d, S) from the recurrence and compare with the paper."""
    n = 1 << d
    full = (1 << n) - 1
    nbr = hypercube(d)
    lay, clo = forcing_layers(nbr, full, to_mask(listed))
    got = len(lay) if clo == full else None
    counted = popcount(to_mask(listed)) + sum(popcount(m) for m in lay)
    # C_{T-1} != V is automatic: the loop stops the first time C_t = V, and
    # every layer is nonempty, so no round is wasted.
    ok = (clo == full and got == claimed and counted == n
          and all(m != 0 for m in lay))
    return ck("%s_pt_equals_%d" % (tag, claimed), ok,
              "computed pt(Q_%d,S)=%s, paper says %d" % (d, got, claimed))


# -------------------------------------------------------------------- check 8
def check_refutation(tag, d, listed):
    """THE LOAD-BEARING CHECK.  pt is recomputed and compared, as integers,
    with the conjectured maximum 2^(d-2).  Since PT(Q_d) is a maximum over
    minimum zero forcing sets and S is one of them, pt(Q_d,S) > 2^(d-2)
    contradicts PT(Q_d) = 2^(d-2)."""
    n = 1 << d
    full = (1 << n) - 1
    nbr = hypercube(d)
    lay, clo = forcing_layers(nbr, full, to_mask(listed))
    if clo != full:
        return ck("%s_refutes_conjectured_formula" % tag, False,
                  "not a zero forcing set, nothing is refuted")
    pt = len(lay)
    conj = conjectured_PT(d)
    minimum_card = len(set(listed)) == cited_Z(d)
    strict = pt > conj
    return ck("%s_refutes_conjectured_formula" % tag,
              strict and minimum_card,
              "PT(Q_%d) >= pt = %d > 2^(%d-2) = %d, excess %d, |S| = Z(Q_%d)"
              % (d, pt, d, conj, pt - conj, d)
              if strict and minimum_card
              else "pt=%d vs conjectured %d, |S|=%d vs Z=%d"
              % (pt, conj, len(set(listed)), cited_Z(d)))


# ----------------------------------------------------------------- check 9/10
def check_layers_match_paper(tag, d, listed, paper_layers):
    """Every layer D_t printed in the paper equals the recomputed layer."""
    n = 1 << d
    full = (1 << n) - 1
    nbr = hypercube(d)
    lay, clo = forcing_layers(nbr, full, to_mask(listed))
    got = [to_list(m) for m in lay]
    want = [sorted(x) for x in paper_layers]
    if len(got) != len(want):
        return ck("%s_layers_match_paper" % tag, False,
                  "%d recomputed layers, %d printed" % (len(got), len(want)))
    diffs = [t + 1 for t in range(len(got)) if got[t] != want[t]]
    if diffs:
        t = diffs[0] - 1
        return ck("%s_layers_match_paper" % tag, False,
                  "first mismatch at D_%d: computed %s, printed %s"
                  % (diffs[0], got[t], want[t]))
    print("info: all %d layers of Q_%d agree with the printed table; "
          "layer sizes %s" % (len(got), d, [len(x) for x in got]))
    return ck("%s_layers_match_paper" % tag, True,
              "D_1..D_%d identical to the paper" % len(got))


# ------------------------------------------------------------------- check 11
def check_layer_partition(tag, d, listed, paper_layers):
    """The layers printed in the paper are pairwise disjoint, disjoint from S,
    and their union is exactly V_d minus S -- checked on the paper's own data,
    independently of the recomputation."""
    n = 1 << d
    s = to_mask(listed)
    acc = 0
    overlap = False
    for lay in paper_layers:
        m = to_mask(lay)
        if len(set(lay)) != len(lay) or (m & acc) or (m & s) or m == 0:
            overlap = True
        acc |= m
    total = popcount(acc) + popcount(s)
    ok = (not overlap) and acc == ((1 << n) - 1) & ~s and total == n
    return ck("%s_paper_layer_partition" % tag, ok,
              "%d + %d = %d = 2^%d, layers pairwise disjoint"
              % (popcount(s), popcount(acc), total, d)
              if ok else "union/disjointness failure (overlap=%s, covered=%d)"
              % (overlap, popcount(acc)))


# ------------------------------------------------------------------- check 12
def check_witnesses(tag, d, listed):
    """Each vertex of each layer has a witness u in C_{t-1} with
    N(u) \\ C_{t-1} = {v}, and witnesses of distinct targets in the same layer
    are distinct -- the independence of the simultaneous forces."""
    n = 1 << d
    full = (1 << n) - 1
    nbr = hypercube(d)
    filled = to_mask(listed)
    lay, clo = forcing_layers(nbr, full, filled)
    total, bad = 0, []
    for t, m in enumerate(lay):
        wit = {}
        for v in to_list(m):
            ws = [u for u in to_list(filled)
                  if nbr[u] & ~filled == (1 << v)]
            if not ws:
                bad.append("D_%d target %d has no witness" % (t + 1, v))
                continue
            wit[v] = ws
            total += 1
        seen = {}
        for v, ws in wit.items():
            for u in ws:
                if u in seen:
                    bad.append("D_%d: vertex %d witnesses both %d and %d"
                               % (t + 1, u, seen[u], v))
                seen[u] = v
        filled |= m
    ok = (not bad) and total == n - len(set(listed)) and filled == full
    return ck("%s_witness_certificates" % tag, ok,
              "%d forces, each with a private witness" % total
              if ok else "; ".join(bad[:3]) or "count mismatch")


# ------------------------------------------------------------------- check 13
def check_inclusionwise_minimal(tag, d, listed):
    """No proper subset obtained by deleting one vertex is zero forcing, so S
    is inclusion-wise minimal.  (Were some deletion forcing, Z(Q_d) would be
    below 2^(d-1) and the cited identity would be contradicted.)"""
    n = 1 << d
    full = (1 << n) - 1
    nbr = hypercube(d)
    s = to_mask(listed)
    survivors = []
    for v in to_list(s):
        if is_forcing(nbr, full, s & ~(1 << v)):
            survivors.append(v)
    return ck("%s_inclusionwise_minimal" % tag, not survivors,
              "all %d single-vertex deletions fail to force" % popcount(s)
              if not survivors
              else "deleting %s still forces Q_%d" % (survivors[:5], d))


# ------------------------------------------------------------------- check 14
def check_halfcube_upper_bound():
    """Z(Q_d) <= 2^(d-1) for 2 <= d <= 7, constructively: the half cube
    {v : v < 2^(d-1)} has 2^(d-1) vertices and every one of them has its
    partner as unique unfilled neighbour, so it forces in a single round."""
    rows, bad = [], []
    for d in range(2, 8):
        n = 1 << d
        full = (1 << n) - 1
        nbr = hypercube(d)
        half = to_mask(range(n // 2))
        lay, clo = forcing_layers(nbr, full, half)
        if clo != full or len(lay) != 1 or popcount(half) != cited_Z(d):
            bad.append("d=%d" % d)
        rows.append("d=%d:|S|=%d,pt=%d" % (d, popcount(half), len(lay)))
    return ck("halfcube_upper_bound", not bad,
              "Z(Q_d) <= 2^(d-1) for d=2..7 (" + " ".join(rows) + ")"
              if not bad else "failed at " + ",".join(bad))


# ------------------------------------------------------------------- check 15
def check_zfn_exhaustive_small_d():
    """Z(Q_d) = 2^(d-1) by brute force for d = 2,3,4: no subset of size below
    2^(d-1) is zero forcing, and some subset of that size is."""
    rows, bad = [], []
    for d in (2, 3, 4):
        n = 1 << d
        full = (1 << n) - 1
        nbr = hypercube(d)
        z = None
        for k in range(1, n + 1):
            hit = False
            for comb in combinations(range(n), k):
                if is_forcing(nbr, full, to_mask(comb)):
                    hit = True
                    break
            if hit:
                z = k
                break
        if z != cited_Z(d):
            bad.append("Z(Q_%d)=%s not %d" % (d, z, cited_Z(d)))
        rows.append("Z(Q_%d)=%d" % (d, z))
    return ck("zfn_exhaustive_small_d", not bad,
              "exhaustive: " + " ".join(rows)
              if not bad else "; ".join(bad))


# ------------------------------------------------------------------- check 16
def check_PT_exhaustive_small_d():
    """PT(Q_d) = 2^(d-2) for d = 2,3,4 by a complete census of every minimum
    zero forcing set.  This is the regime where the formula is reported to
    hold, so a counterexample here would indict the formula earlier than the
    paper claims, and its absence supports the minimality of dimension 6."""
    rows, bad = [], []
    for d in (2, 3, 4):
        n = 1 << d
        full = (1 << n) - 1
        nbr = hypercube(d)
        best, count = 0, 0
        for comb in combinations(range(n), cited_Z(d)):
            lay, clo = forcing_layers(nbr, full, to_mask(comb))
            if clo == full:
                count += 1
                if len(lay) > best:
                    best = len(lay)
        if best != conjectured_PT(d):
            bad.append("PT(Q_%d)=%d not %d" % (d, best, conjectured_PT(d)))
        rows.append("PT(Q_%d)=%d over %d minimum sets" % (d, best, count))
    return ck("PT_exhaustive_small_d", not bad,
              "; ".join(rows) if not bad else "; ".join(bad))


# ------------------------------------------------------------------- check 17
def check_PT_q5_search():
    """Q_5 is beyond a complete census (C(32,16) = 601080390 candidates), so a
    fully deterministic ascent is run instead: from each of the 32 xor
    translates of the half cube, repeatedly apply the first swap of one chosen
    and one unchosen vertex that increases (vertices forced, rounds used),
    scanning vertices in increasing order.  No randomness, so the result is
    identical on every interpreter.  The best value found must equal 2^(5-2),
    the value the formula predicts for Q_5: less would mean the formula is not
    even attained, more would mean the formula already fails at d = 5."""
    d, n = 5, 32
    full = (1 << n) - 1
    nbr = hypercube(d)

    def score(mask):
        lay, clo = forcing_layers(nbr, full, mask)
        return (popcount(clo), len(lay))

    best, seeds = 0, 0
    for shift in range(n):
        cur = to_mask((v ^ shift) for v in range(n // 2))
        val = score(cur)
        seeds += 1
        moved = True
        while moved:
            moved = False
            for a in to_list(cur):
                for b in to_list(full & ~cur):
                    cand = cur ^ (1 << a) ^ (1 << b)
                    sc = score(cand)
                    if sc > val:
                        cur, val, moved = cand, sc, True
                        break
                if moved:
                    break
        if val[0] == n and val[1] > best:
            best = val[1]
    target = conjectured_PT(d)
    return ck("PT_q5_deterministic_search", best == target,
              "best pt over %d deterministic ascents = %d = 2^(5-2)"
              % (seeds, best) if best == target
              else "search found %d, formula predicts %d" % (best, target))


def guarded(label, fn, *args):
    """A check that raises is a failed check, never a silent exit."""
    try:
        return fn(*args)
    except Exception as exc:
        return ck(label + "_raised", False,
                  type(exc).__name__ + ": " + str(exc)[:90])


def main():
    print("info: hypercube maximum propagation time; two exhibited minimum "
          "zero forcing sets, in Q_6 and Q_7")
    plan = [
        ("s6_wellformed", check_wellformed, ("s6", 6, S6_PAPER)),
        ("s7_wellformed", check_wellformed, ("s7", 7, S7_PAPER)),
        ("hypercube_structure", check_hypercube_structure, ()),
        ("forcing_engine_sanity", check_engine_sanity, ()),
        ("s6_is_zero_forcing", check_is_zero_forcing, ("s6", 6, S6_PAPER)),
        ("s7_is_zero_forcing", check_is_zero_forcing, ("s7", 7, S7_PAPER)),
        ("s6_pt", check_propagation_time, ("s6", 6, S6_PAPER, PT6_PAPER)),
        ("s7_pt", check_propagation_time, ("s7", 7, S7_PAPER, PT7_PAPER)),
        ("s6_refutes", check_refutation, ("s6", 6, S6_PAPER)),
        ("s7_refutes", check_refutation, ("s7", 7, S7_PAPER)),
        ("s6_layers", check_layers_match_paper,
         ("s6", 6, S6_PAPER, LAYERS6_PAPER)),
        ("s7_layers", check_layers_match_paper,
         ("s7", 7, S7_PAPER, LAYERS7_PAPER)),
        ("s6_partition", check_layer_partition,
         ("s6", 6, S6_PAPER, LAYERS6_PAPER)),
        ("s7_partition", check_layer_partition,
         ("s7", 7, S7_PAPER, LAYERS7_PAPER)),
        ("s6_witnesses", check_witnesses, ("s6", 6, S6_PAPER)),
        ("s7_witnesses", check_witnesses, ("s7", 7, S7_PAPER)),
        ("s6_minimal", check_inclusionwise_minimal, ("s6", 6, S6_PAPER)),
        ("s7_minimal", check_inclusionwise_minimal, ("s7", 7, S7_PAPER)),
        ("halfcube_upper_bound", check_halfcube_upper_bound, ()),
        ("zfn_exhaustive_small_d", check_zfn_exhaustive_small_d, ()),
        ("PT_exhaustive_small_d", check_PT_exhaustive_small_d, ()),
        ("PT_q5_search", check_PT_q5_search, ()),
    ]
    for label, fn, args in plan:
        guarded(label, fn, *args)
    print("info: NOT RE-RUN -- the lower bound Z(Q_d) >= 2^(d-1) for "
          "d = 5,6,7 is taken from the literature, not recomputed here; "
          "only the upper bound is reproved above, and only d <= 4 is "
          "settled exhaustively.")
    print("info: NOT RE-RUN -- no complete census of the minimum zero "
          "forcing sets of Q_5, Q_6 or Q_7 is attempted (C(32,16), C(64,32) "
          "and C(128,64) candidates), so the exact values of PT(Q_6) and "
          "PT(Q_7) are not determined; only the lower bounds 18 and 43 are "
          "established, which is exactly what the paper claims.")
    print("info: NOT RE-RUN -- the separate interval assertion of the "
          "original conjecture is outside the paper's claim and is not "
          "examined.")
    finish()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verification program for a settlement of the Coprime Buratti-Horak-Rosa
Conjecture at the support {1,6,18}.  Standard library only, exact integer
arithmetic, no data files, no network.

Setting: the vertices of K_v are 0,...,v-1, the cyclic length of an edge is
l_v(x,y) = min(|x-y|, v-|x-y|), and a Hamiltonian path P = (p_0,...,p_{v-1})
realizes the multiset Delta_v(P) of the v-1 cyclic lengths of its edges.
A multiset L of size v-1 is realizable if some Hamiltonian path realizes it.

TAKEN FROM THE PAPER (inputs, and nothing else is assumed):
  * the support {1,6,18} and the three orders left open, v = 47, 49, 59;
  * the census table |T_v| and |T_v^+| and the totals 4123 and 3667;
  * the split of the census by reduction that the paper reports, 3792 targets
    settled by s = 6^-1, 43 by s = 18^-1 and 288 by s = 1, used only as a
    number to compare against when the run is asked, with --paper-protocol,
    to follow the paper's own single-budget protocol;
  * the reduction table of images of (1,6,18) under multiplication by the
    inverses of 6 and of 18, modulo each v;
  * the one certificate the paper prints in full, a 59-term linear
    realization at v = 59 for (a,b,c) = (13,6,39), and the claim that
    multiplying it by 6 modulo 59 realizes {1^13, 6^6, 18^39};
  * the search rules of the paper's Section on certificate verification,
    reimplemented here from that description.

DERIVED HERE (every line below is computed, never asserted):
  * the census, by enumeration and by a second independent enumeration,
    compared with the paper's table, its binomial formulas and its totals;
  * that all hypotheses of the theorem hold at the three orders, and that
    they are equivalent to the conjecture's own conditions for this support;
  * that the printed certificate is a permutation of 0,...,58, is linear,
    realizes the stated reduced multiset, and that its image under
    multiplication by 6 is a Hamiltonian path of K_59 whose cyclic-length
    multiset is exactly {1^13, 6^6, 18^39}, the conclusion for that target;
  * the reduction table, recomputed from modular inverses, and the scaling
    lemma, checked on all 58 units of Z_59;
  * the closed-form certificates at the corners of each T_v;
  * an independent re-run of the search: by default over the whole census,
    every target of every T_v, a Hamiltonian path is produced and then
    verified exactly against the target multiset.  (--stride=N re-runs a
    deterministic part instead, namely all 456 targets with a zero exponent
    plus every N-th other target; the number of targets actually attempted
    is printed, and what was left out is named on the NOT RE-RUN line.)  The
    node cap climbs a ladder whose last rung is, unless --budget= says
    otherwise, the paper's own budget of 10^5 nodes per starting vertex, so
    that no reduction here is allotted less than the paper allots itself;
    the ladder actually in force, and whether it meets, falls short of or
    exceeds the paper's budget, is printed by the run itself, and if the
    re-run still fails to settle a target it is named on a RESIDUAL line,
    with the budget it was given.  --paper-protocol drops the cheap rungs
    and searches at a single rung equal to that budget, which is exactly the
    protocol the paper describes, so that with it the split of the census by
    reduction printed by the run is directly comparable with the paper's
    reported 3792/43/288, and the run prints that comparison;
  * that the search engine itself is sound and complete enough to be
    evidence: it agrees with unpruned exhaustive enumeration on a small
    analogue, it fails, as it must, on an unrealizable multiset, and the
    cap ladder demonstrably escalates from a rung too small to settle a
    target to a rung that settles it.

A line beginning NOT RE-RUN reports what this program does not reproduce.
"""

import sys
from math import gcd

# ---------------------------------------------------------------- paper inputs
ORDERS = (47, 49, 59)              # the three orders left open, from the paper
SUPPORT = (1, 6, 18)               # the support the paper treats
CENSUS_TABLE = {47: (1128, 990), 49: (1225, 1081), 59: (1770, 1596)}
CENSUS_TOTAL = (4123, 3667)
# the paper's reported split of the census by reduction: s -> number of targets
# it settles.  Compared against only when the run is following the paper's own
# protocol (a single node cap, reductions in the paper's order, whole census).
PAPER_SPLIT = {6: 3792, 18: 43, 1: 288}
# rows of the paper's reduction table: v -> (s=6^-1 image, s=18^-1 image),
# each image being the lengths that 1, 6, 18 are sent to, in that order
IMAGE_TABLE = {47: ((8, 1, 3), (13, 16, 1)),
               49: ((8, 1, 3), (19, 16, 1)),
               59: ((10, 1, 3), (23, 20, 1))}
# the one certificate the paper prints, at v = 59, (a, b, c) = (13, 6, 39)
EXHIBITED_ORDER = 59
EXHIBITED_ABC = (13, 6, 39)
EXHIBITED_PATH = (
    0, 1, 2, 3, 4, 5, 6, 9, 12, 15, 18, 8, 11, 14, 17, 7, 10, 13, 16, 19,
    22, 25, 28, 31, 21, 24, 27, 30, 20, 23, 26, 29, 32, 35, 38, 41, 51, 48,
    58, 55, 45, 42, 52, 49, 46, 56, 53, 50, 47, 57, 54, 44, 34, 37, 40, 43,
    33, 36, 39)
EXHIBITED_LINEAR_MULTISET = {1: 6, 3: 39, 10: 13}
EXHIBITED_SCALE = 6                # p_i = 6 x_i mod 59

# how much of the census the re-run covers: stride 1 is the whole of it, which
# is the default, so that no target of the census is left on the paper's word.
# --stride=N (N > 1) asks for a cheaper partial pass instead.
DEFAULT_STRIDE = 1
# the paper's own search budget: "each starting vertex is allotted 10^5 nodes
# of the search tree".  Taken from the paper, and used as the last rung below.
PAPER_BUDGET_PER_START = 10 ** 5
# the re-run climbs this ladder of node caps: cheap rungs first, so that easy
# targets cost little, and a last rung equal to the paper's budget, so that on
# each reduction this search is allotted no less than the paper allots itself.
# Overridable with --budget=N, which replaces the last rung (and is reported).
# --paper-protocol drops the cheap rungs and uses the single rung N alone, which
# is what the paper does; the ladder actually in force is always printed.
CAP_LADDER = (500, 5000, PAPER_BUDGET_PER_START)
# how many residual (unsettled) targets to name individually
RESIDUAL_PRINT_LIMIT = 100

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    sys.stdout.flush()
    return bool(ok)


def hat(g, v):
    """Cyclic length of a difference g modulo v."""
    g %= v
    return min(g, v - g)


def counts(items):
    """Multiset of a list, as a dict, using exact integer counting."""
    out = {}
    for x in items:
        out[x] = out.get(x, 0) + 1
    return out


def cyclic_multiset(path, v):
    """Multiset Delta_v(P) of the v-1 cyclic edge lengths of a path."""
    return counts([hat(path[i] - path[i + 1], v) for i in range(len(path) - 1)])


def linear_multiset(path):
    """Multiset of absolute differences along a path."""
    return counts([abs(path[i] - path[i + 1]) for i in range(len(path) - 1)])


def is_permutation(path, v):
    return len(path) == v and sorted(path) == list(range(v))


def target_multiset(a, b, c):
    """{1^a, 6^b, 18^c} as a dict, zero exponents omitted."""
    out = {}
    for length, mult in zip(SUPPORT, (a, b, c)):
        if mult:
            out[length] = out.get(length, 0) + mult
    return out


def inverse_mod(s, v):
    """Inverse of s modulo v by exhaustive exact search; None if none exists."""
    s %= v
    for t in range(v):
        if (s * t) % v == 1:
            return t
    return None


def triples(v, positive_only):
    """T_v, or T_v^+, enumerated: (a,b,c) >= 0 with a+b+c = v-1."""
    lo = 1 if positive_only else 0
    out = []
    for a in range(lo, v):
        for b in range(lo, v - a):
            c = v - 1 - a - b
            if c >= lo:
                out.append((a, b, c))
    return out


def binom2(n):
    return n * (n - 1) // 2


def check_census():
    """Enumerate T_v and T_v^+ and compare with the paper's table."""
    tot_all = tot_pos = 0
    ok_rows, ok_binom, detail = True, True, []
    for v in ORDERS:
        n_all = len(triples(v, False))
        n_pos = len(triples(v, True))
        tot_all += n_all
        tot_pos += n_pos
        if (n_all, n_pos) != CENSUS_TABLE[v]:
            ok_rows = False
        if n_all != binom2(v + 1) or n_pos != binom2(v - 2):
            ok_binom = False
        detail.append("v=%d:%d/%d" % (v, n_all, n_pos))
    ck("census_rows_match_table", ok_rows, " ".join(detail))
    ck("census_binomial_formulas", ok_binom,
       "|T_v|=C(v+1,2), |T_v^+|=C(v-2,2)")
    ck("census_totals_4123_3667", (tot_all, tot_pos) == CENSUS_TOTAL,
       "enumerated %d / %d" % (tot_all, tot_pos))
    # triples with a zero exponent, i.e. support strictly inside {1,6,18}
    degenerate = tot_all - tot_pos
    expected = sum(3 * (v - 2) + 3 for v in ORDERS)
    ck("census_degenerate_count", degenerate == expected == 456,
       "%d triples have a zero exponent" % degenerate)
    return tot_all, tot_pos


def check_hypotheses():
    """Every hypothesis of the theorem, at each of the three orders."""
    ok, detail = True, []
    for v in ORDERS:
        half = v // 2
        h1 = 18 <= half                      # supp(L) inside {1,...,floor(v/2)}
        h2 = gcd(v, 18) == 1                 # coprimality of the whole support
        h3 = all(gcd(v, x) == 1 for x in SUPPORT)
        h4 = all(x <= half for x in SUPPORT)
        if not (h1 and h2 and h3 and h4):
            ok = False
        detail.append("v=%d half=%d" % (v, half))
    ck("orders_satisfy_theorem_hypotheses", ok, " ".join(detail))
    # the theorem's two hypotheses are exactly the conjecture's, for this
    # support: check the equivalence order by order over a wide range
    same = True
    for v in range(3, 601):
        stated = (18 <= v // 2) and gcd(v, 18) == 1
        conjecture = all(x <= v // 2 and gcd(v, x) == 1 for x in SUPPORT)
        if stated != conjecture:
            same = False
    ck("hypotheses_equal_conjecture_form", same, "checked v=3..600")
    # the census covers every case the conjecture asserts at these orders:
    # enumerate, by a second and independent method, all multisets of size
    # v-1 with support inside {1,6,18} and compare with T_v
    complete, detail2 = True, []
    for v in ORDERS:
        seen = set()
        for a in range(v + 1):
            for b in range(v + 1):
                for c in range(v + 1):
                    if a + b + c == v - 1:
                        seen.add((a, b, c))
        if seen != set(triples(v, False)):
            complete = False
        detail2.append("v=%d:%d" % (v, len(seen)))
    ck("census_covers_every_asserted_case", complete, " ".join(detail2))
    return True


def check_exhibited():
    """Decode, count and re-print the one certificate the paper displays."""
    v = EXHIBITED_ORDER
    x = list(EXHIBITED_PATH)
    a, b, c = EXHIBITED_ABC
    print("exhibited linear realization at v=%d, (a,b,c)=(%d,%d,%d):"
          % (v, a, b, c))
    print("  " + ", ".join(str(t) for t in x))
    ck("exhibited_is_permutation", is_permutation(x, v),
       "%d entries, sorted = 0..%d" % (len(x), v - 1))
    ck("exhibited_length_is_v", len(x) == v and a + b + c == v - 1,
       "path length %d, a+b+c = %d" % (len(x), a + b + c))
    # linear: every step stays within floor(v/2), so cyclic length = |difference|
    steps = [abs(x[i] - x[i + 1]) for i in range(v - 1)]
    ck("exhibited_is_linear", all(d <= v // 2 for d in steps),
       "max step %d <= %d" % (max(steps), v // 2))
    lin = linear_multiset(x)
    ck("exhibited_linear_multiset", lin == EXHIBITED_LINEAR_MULTISET,
       " ".join("%d^%d" % kv for kv in sorted(lin.items())))
    # the reduction: p_i = 6 x_i mod 59 must be a Hamiltonian path of K_59
    p = [(EXHIBITED_SCALE * t) % v for t in x]
    ck("exhibited_image_is_permutation", is_permutation(p, v),
       "6x mod 59 is a bijection")
    got = cyclic_multiset(p, v)
    want = target_multiset(a, b, c)
    ck("exhibited_image_realizes_target", got == want,
       "Delta = " + " ".join("%d^%d" % kv for kv in sorted(got.items()))
       + " ; target " + " ".join("%d^%d" % kv for kv in sorted(want.items())))
    # the multiset it realizes obeys the conjecture's own two conditions
    ck("exhibited_target_is_admissible",
       all(g <= v // 2 and gcd(v, g) == 1 for g in got) and
       sum(got.values()) == v - 1,
       "all lengths <= %d, coprime to %d, %d edges"
       % (v // 2, v, sum(got.values())))
    return True


def check_image_table():
    """Recompute the paper's reduction table from scratch."""
    ok_tab, ok_bound, detail = True, True, []
    for v in ORDERS:
        row = []
        for s_base in (6, 18):
            s = inverse_mod(s_base, v)
            if s is None:
                ok_tab = False
                continue
            image = tuple(hat(s * g, v) for g in SUPPORT)
            row.append(image)
        if tuple(row) != IMAGE_TABLE[v]:
            ok_tab = False
        detail.append("v=%d %s" % (v, row))
        # s = 1 row is the identity
        if tuple(hat(g, v) for g in SUPPORT) != SUPPORT:
            ok_tab = False
        # the paper's claim about the lengths it displays: every entry of the
        # printed table, and of the s=1 row, lies in {1,...,floor(v/2)}.  The
        # numbers tested here are the paper's own, not this program's output,
        # so a table entry outside the range makes this check fail.
        displayed = [g for img in IMAGE_TABLE[v] for g in img] + list(SUPPORT)
        if not all(1 <= g <= v // 2 for g in displayed):
            ok_bound = False
    ck("reduction_table_recomputed", ok_tab, " ".join(detail))
    ck("reduction_images_within_half_v", ok_bound,
       "every length printed in the table lies in 1..floor(v/2)")
    # multiplication by s is an automorphism of Z_v when gcd(s,v)=1: verify
    # that it permutes Z_v and sends lengths as the table says, at each order
    ok_auto = True
    for v in ORDERS:
        for s_base in (1, 6, 18):
            s = inverse_mod(s_base, v)
            if sorted((s * t) % v for t in range(v)) != list(range(v)):
                ok_auto = False
            if hat(s * s_base, v) != 1:
                ok_auto = False
        # the coprimality hypothesis is load-bearing, not decoration: every
        # residue sharing a factor with v must fail to permute Z_v
        for s in range(v):
            if gcd(s, v) == 1:
                continue
            if sorted((s * t) % v for t in range(v)) == list(range(v)):
                ok_auto = False
    ck("scaling_is_automorphism", ok_auto,
       "units permute Z_v, non-units do not")
    return True


def check_corners():
    """p_i = g i mod v settles the corners of T_v, for g = 1, 6, 18."""
    ok, detail = True, []
    for v in ORDERS:
        for j, g in enumerate(SUPPORT):
            p = [(g * i) % v for i in range(v)]
            abc = [0, 0, 0]
            abc[j] = v - 1
            good = (is_permutation(p, v) and
                    cyclic_multiset(p, v) == target_multiset(*abc))
            if not good:
                ok = False
            detail.append("%d:%d^%d%s" % (v, g, v - 1, "" if good else "!"))
    ck("corner_closed_forms", ok, " ".join(detail))
    return True


def build_neighbours(v, lengths):
    """Bitmask neighbour lists of the graph on 0..v-1 joining u,w iff
    |u-w| is one of the given lengths."""
    nbr = [0] * v
    for u in range(v):
        m = 0
        for d in lengths:
            if u - d >= 0:
                m |= 1 << (u - d)
            if u + d < v:
                m |= 1 << (u + d)
        nbr[u] = m
    return nbr


def connected_and_degrees_ok(mask, start, nbr):
    """True if the vertices in mask induce a connected subgraph and at most
    one vertex other than start has fewer than two neighbours in it."""
    seen = 1 << start
    frontier = seen
    while frontier:
        nxt = 0
        f = frontier
        while f:
            low = f & -f
            nxt |= nbr[low.bit_length() - 1]
            f ^= low
        nxt &= mask & ~seen
        seen |= nxt
        frontier = nxt
    if seen != mask:
        return False
    thin = 0
    m = mask & ~(1 << start)
    while m:
        low = m & -m
        i = low.bit_length() - 1
        m ^= low
        deg = nbr[i] & mask
        if deg & (deg - 1) == 0:      # zero or one neighbour
            thin += 1
            if thin > 1:
                return False
    return True


class _Budget(Exception):
    pass


def find_linear(v, mult, budget_per_start):
    """Search for a linear realization of the multiset mult on 0..v-1, by the
    rules of the paper: start at t, extend to p +- d for an unused length d,
    candidates ordered by how many extensions they admit then by label, with
    connectivity and thin-vertex pruning, and a node cap per start.  The
    second component returned is the total number of search-tree nodes visited
    over all starting vertices, the quantity the paper counts."""
    if sum(mult.values()) != v - 1 or any(d > v // 2 for d in mult):
        return None, 0
    keys = tuple(sorted(mult))
    all_mask = (1 << v) - 1
    cache = {}

    def nbr_for(rem):
        key = tuple(d for d in keys if rem[d] > 0)
        got = cache.get(key)
        if got is None:
            got = build_neighbours(v, key)
            cache[key] = got
        return key, got

    spent = [0]                        # nodes spent by the current start
    total = [0]                        # nodes spent over all starts

    def dfs(p, used, rem, depth, path):
        if depth == v:
            return True
        spent[0] += 1
        total[0] += 1
        if spent[0] > budget_per_start:
            raise _Budget
        key, nbr = nbr_for(rem)
        scored = []
        for d in key:
            for q in (p - d, p + d):
                if 0 <= q < v and not (used >> q) & 1:
                    rem[d] -= 1
                    nbr2 = nbr_for(rem)[1]
                    free = all_mask & ~used & ~(1 << q)
                    deg = bin(nbr2[q] & free).count("1")
                    keep = (depth + 1 == v) or (
                        deg > 0 and
                        connected_and_degrees_ok(free | (1 << q), q, nbr2))
                    rem[d] += 1
                    if keep:
                        scored.append((deg, q, d))
        scored.sort()
        for _, q, d in scored:
            rem[d] -= 1
            path.append(q)
            if dfs(q, used | (1 << q), rem, depth + 1, path):
                return True
            path.pop()
            rem[d] += 1
        return False

    for t in range(v):
        rem = dict(mult)
        spent[0] = 0
        path = [t]
        try:
            if dfs(t, 1 << t, rem, 1, path):
                return path, total[0]
        except _Budget:
            pass
    return None, total[0]


def reduced_multiset(v, a, b, c, s_base):
    """The multiset {1^a,6^b,18^c} pushed through multiplication by s_base^-1."""
    s = inverse_mod(s_base, v)
    mult = {}
    for g, m in zip(SUPPORT, (a, b, c)):
        if m:
            h = hat(s * g, v)
            mult[h] = mult.get(h, 0) + m
    return mult


def settle(v, a, b, c, caps=CAP_LADDER):
    """Try the paper's three reductions, in its order, at increasing node caps.
    Returns (linear realization, scale s_base, cap that produced it, total
    nodes spent).  The path is returned unverified: the caller verifies it.
    (None, None, last cap, nodes) if no reduction produced a path even at the
    last rung of the ladder."""
    spent = 0
    for cap in caps:
        for s_base in (6, 18, 1):
            mult = reduced_multiset(v, a, b, c, s_base)
            x, n = find_linear(v, mult, cap)
            spent += n
            if x is not None:
                return x, s_base, cap, spent
    return None, None, caps[-1], spent


def verify_target(v, a, b, c, x, s_base):
    """Exact check that x is a linear realization of the reduced multiset and
    that its image under multiplication by s_base is a Hamiltonian path of K_v
    realizing {1^a,6^b,18^c}.  Returns a list of failed conditions."""
    bad = []
    if not is_permutation(x, v):
        bad.append("linear-not-permutation")
    if any(abs(x[i] - x[i + 1]) > v // 2 for i in range(len(x) - 1)):
        bad.append("not-linear")
    if linear_multiset(x) != reduced_multiset(v, a, b, c, s_base):
        bad.append("wrong-reduced-multiset")
    p = [(s_base * t) % v for t in x]
    if not is_permutation(p, v):
        bad.append("image-not-permutation")
    if cyclic_multiset(p, v) != target_multiset(a, b, c):
        bad.append("wrong-cyclic-multiset")
    return bad


def sample_targets(stride):
    """A deterministic sample of the census: every triple with a zero exponent
    (support strictly inside {1,6,18}), plus every stride-th other triple."""
    out = []
    for v in ORDERS:
        chosen = []
        for i, t in enumerate(triples(v, False)):
            if min(t) == 0 or i % stride == 0:
                chosen.append(t)
        out.append((v, chosen))
    return out


def check_ladder():
    """The cap ladder must actually escalate.  A rung of one node cannot settle
    anything at v=47, so a target settled by a later rung exercises the
    escalation path, and the path that rung returns must verify exactly."""
    v, (a, b, c) = 47, (46, 0, 0)
    x0, _s0, _cap0, n0 = settle(v, a, b, c, caps=(1,))
    x1, s1, cap1, n1 = settle(v, a, b, c, caps=(1, 500))
    ok = (x0 is None and x1 is not None and cap1 == 500 and
          not verify_target(v, a, b, c, x1, s1))
    ck("cap_ladder_escalates_to_later_rung", ok,
       "v=47 (46,0,0): a rung of 1 node returns no path after %d nodes, the "
       "next rung returns a verified one after %d, credited to s=%s; default "
       "ladder %s, whose last rung is the paper's %d nodes per start"
       % (n0, n1, s1, list(CAP_LADDER), PAPER_BUDGET_PER_START))
    return True


def check_search(stride, caps, single_rung, census_size):
    """Re-run the paper's search over the census (the whole of it at stride 1)
    and verify every path it returns in exact integer arithmetic."""
    attempted = settled = 0
    unsettled, defects, by_scale, by_rung = [], [], {}, {}
    degen_total = degen_settled = 0
    worst_nodes, worst_target, node_total = 0, None, 0
    for v, chosen in sample_targets(stride):
        for (a, b, c) in chosen:
            attempted += 1
            deg = min((a, b, c)) == 0
            degen_total += 1 if deg else 0
            x, s_base, cap, nodes = settle(v, a, b, c, caps=caps)
            node_total += nodes
            if nodes > worst_nodes:
                worst_nodes, worst_target = nodes, (v, a, b, c)
            if x is None:
                unsettled.append((v, a, b, c))
                continue
            bad = verify_target(v, a, b, c, x, s_base)
            if bad:
                defects.append((v, a, b, c, ",".join(bad)))
                continue
            settled += 1
            degen_settled += 1 if deg else 0
            by_scale[s_base] = by_scale.get(s_base, 0) + 1
            by_rung[cap] = by_rung.get(cap, 0) + 1
            if settled % 250 == 0:
                print("  ... %d targets settled and verified so far" % settled)
                sys.stdout.flush()
    ck("search_returns_verified_paths", not defects,
       "%d verified, %d failed verification" % (settled, len(defects)))
    # every attempted target must have been settled AND verified: a target for
    # which a path was produced but failed verification counts against this
    ck("search_settles_every_attempted_target",
       settled == attempted and attempted > 0 and not unsettled,
       "%d of %d attempted targets settled at up to %d nodes per start; %s"
       % (settled, attempted, caps[-1],
          "none unsettled" if not unsettled else
          "%d unsettled, named on the RESIDUAL lines below" % len(unsettled)))
    ck("degenerate_support_targets_all_settled",
       degen_total == 456 and degen_settled == 456,
       "%d of %d zero-exponent triples settled" % (degen_settled, degen_total))
    print("scales used by the re-run: " +
          " ".join("s=%d:%d" % kv for kv in sorted(by_scale.items())))
    print("cap rung that settled each target: " +
          " ".join("<=%d:%d" % kv for kv in sorted(by_rung.items())))
    # the paper's split of the census by reduction is reproducible only under the
    # paper's own protocol: one node cap, its order of reductions, whole census.
    split = tuple(by_scale.get(s, 0) for s in (6, 18, 1))
    paper = tuple(PAPER_SPLIT[s] for s in (6, 18, 1))
    compared = (single_rung and tuple(caps) == (PAPER_BUDGET_PER_START,) and
                attempted == census_size and not unsettled and not defects)
    if compared:
        print("PAPER-PROTOCOL SPLIT: this run used a single rung of %d nodes per "
              "starting vertex and tried the reductions in the paper's order over "
              "all %d targets, which is the paper's own protocol, so its split by "
              "reduction is directly comparable with the paper's %d/%d/%d; this "
              "run got %d/%d/%d, which %s"
              % (caps[-1], attempted, paper[0], paper[1], paper[2],
                 split[0], split[1], split[2],
                 "agrees with the paper" if split == paper else
                 "does NOT agree, so the candidate ordering reimplemented here "
                 "differs from the paper's in some detail the paper does not "
                 "state; the certificates verified above are unaffected, since "
                 "each was checked exactly against its own target"))
    else:
        print("SPLIT NOT COMPARABLE: the paper's %d/%d/%d split by reduction is "
              "not tested by this run (%s); rerun with --paper-protocol and no "
              "--stride= to compare against it"
              % (paper[0], paper[1], paper[2],
                 "; ".join(
                     ([] if single_rung and
                      tuple(caps) == (PAPER_BUDGET_PER_START,) else
                      ["the cap schedule in force here is %s, not the paper's "
                       "single budget of %d nodes per start, and a target is "
                       "credited to whichever reduction settles it first at the "
                       "rung that settles it"
                       % (list(caps), PAPER_BUDGET_PER_START)]) +
                     ([] if attempted == census_size else
                      ["only %d of the %d targets were attempted"
                       % (attempted, census_size)]) +
                     ([] if not unsettled and not defects else
                      ["%d target(s) were not settled and verified"
                       % (len(unsettled) + len(defects))]))))
    # a rung above the paper's budget settles a target without showing that the
    # paper's own allotment would have sufficed for it: say so, by count
    over_paper = sum(k for cap, k in by_rung.items()
                     if cap > PAPER_BUDGET_PER_START)
    if over_paper:
        print("ABOVE PAPER BUDGET: %d of the %d settled target(s) needed a "
              "rung larger than the paper's %d nodes per starting vertex, so "
              "for those targets this transcript shows a certificate exists "
              "but does not show the paper's own allotment finds it"
              % (over_paper, settled, PAPER_BUDGET_PER_START))
    print("nodes visited by the re-run: %d in total over %d attempted target(s), "
          "most for one target %d%s (the paper reports 4.11e9 in total and "
          "1.51e7 for its largest target, under its single budget of %d nodes "
          "per start%s)"
          % (node_total, attempted, worst_nodes,
             "" if worst_target is None else " at v=%d (a,b,c)=(%d,%d,%d)"
             % worst_target, PAPER_BUDGET_PER_START,
             " and over the whole census, which is the schedule in force here "
             "too, so the two totals are comparable" if compared else
             " rather than the cap schedule in force here"))
    for (v, a, b, c, why) in defects[:RESIDUAL_PRINT_LIMIT]:
        print("  DEFECT v=%d (a,b,c)=(%d,%d,%d) a path was produced but "
              "failed exact verification: %s" % (v, a, b, c, why))
    if unsettled:
        if caps[-1] >= PAPER_BUDGET_PER_START:
            why = ("That budget is at least the paper's own %d nodes per "
                   "start, so either the candidate ordering reimplemented "
                   "here differs from the paper's in some detail it does not "
                   "state, or the paper's claim that every one of the %d "
                   "targets is settled within that budget does not hold at "
                   "these triples."
                   % (PAPER_BUDGET_PER_START, CENSUS_TOTAL[0]))
        else:
            why = ("That budget is below the paper's own %d nodes per start, "
                   "so this says nothing about the paper's claim; drop "
                   "--budget= to test it." % PAPER_BUDGET_PER_START)
        print("RESIDUAL: %d attempted target(s) were not settled by any of "
              "the three reductions within %d search-tree nodes per starting "
              "vertex, and are named below as v (a,b,c).  %s  Nothing checked "
              "elsewhere in this transcript is affected, and every path this "
              "re-run did produce was verified exactly; a referee should read "
              "the named list, not the bare word FAIL, as the finding."
              % (len(unsettled), caps[-1], why))
        for (v, a, b, c) in unsettled[:RESIDUAL_PRINT_LIMIT]:
            print("  RESIDUAL v=%d (a,b,c)=(%d,%d,%d)" % (v, a, b, c))
        if len(unsettled) > RESIDUAL_PRINT_LIMIT:
            print("  RESIDUAL ... and %d more not named individually"
                  % (len(unsettled) - RESIDUAL_PRINT_LIMIT))
    return attempted, settled, unsettled, compared


def brute_linear_exists(v, mult):
    """Exhaustive, unpruned enumeration: does a linear realization exist?"""
    def dfs(p, used, rem, depth):
        if depth == v:
            return True
        for d in list(rem):
            if rem[d] == 0:
                continue
            for q in (p - d, p + d):
                if 0 <= q < v and not (used >> q) & 1:
                    rem[d] -= 1
                    if dfs(q, used | (1 << q), rem, depth + 1):
                        rem[d] += 1
                        return True
                    rem[d] += 1
        return False
    for t in range(v):
        if dfs(t, 1 << t, dict(mult), 1):
            return True
    return False


def check_search_engine():
    """The heuristic search must agree with exhaustive enumeration, and must
    fail on a multiset that no Hamiltonian path can realize."""
    v, sup = 9, (1, 2, 4)               # small analogue, every length <= 4
    agree, realizable, unrealizable = True, 0, 0
    for a in range(v):
        for b in range(v - a):
            c = v - 1 - a - b
            mult = {}
            for g, m in zip(sup, (a, b, c)):
                if m:
                    mult[g] = mult.get(g, 0) + m
            brute = brute_linear_exists(v, mult)
            found, _ = find_linear(v, mult, 20000)
            if brute:
                realizable += 1
                if found is None or linear_multiset(found) != mult:
                    agree = False
            else:
                unrealizable += 1
                if found is not None:
                    agree = False
    ck("search_agrees_with_exhaustive_enumeration", agree,
       "v=9 support {1,2,4}: %d realizable, %d not" % (realizable, unrealizable))
    ck("realizability_is_not_automatic", unrealizable > 0,
       "%d of %d small targets have no linear realization"
       % (unrealizable, realizable + unrealizable))
    # a target no path can realize: steps of 7 in Z_49 leave 7 components
    stuck, _ = find_linear(49, {7: 48}, 200000)
    ck("search_rejects_disconnected_target", stuck is None,
       "no linear realization of {7^48} at v=49")
    return True


def check_automorphism_instance():
    """The reduction lemma, on the certificate the paper prints: for every s
    coprime to v, sP realizes the scaled multiset."""
    v = EXHIBITED_ORDER
    p = [(EXHIBITED_SCALE * t) % v for t in EXHIBITED_PATH]
    base = cyclic_multiset(p, v)
    ok, tried = True, 0
    for s in range(1, v):
        if gcd(s, v) != 1:
            continue
        tried += 1
        want = {}
        for g, m in base.items():
            h = hat(s * g, v)
            want[h] = want.get(h, 0) + m
        q = [(s * t) % v for t in p]
        if not is_permutation(q, v) or cyclic_multiset(q, v) != want:
            ok = False
    ck("scaling_preserves_realization", ok and tried == v - 1,
       "checked all %d units of Z_%d" % (tried, v))
    return True


def main(argv):
    stride = DEFAULT_STRIDE             # 1: the whole census, and the default
    if "--full" in argv:                # kept as an explicit name for stride 1
        stride = 1
    budget = PAPER_BUDGET_PER_START
    # --paper-protocol: one cap per starting vertex and no cheaper rungs, which
    # is the protocol the paper states for itself; with the default budget it
    # reproduces the paper's search exactly, so the split by reduction and the
    # node counts it produces are directly comparable with the paper's.
    single_rung = "--paper-protocol" in argv
    for arg in argv:
        if arg.startswith("--stride="):
            stride = max(1, int(arg.split("=", 1)[1]))
        if arg.startswith("--budget="):
            budget = max(1, int(arg.split("=", 1)[1]))
    if single_rung:
        caps = (budget,)
    else:
        caps = tuple(sorted(set(c for c in CAP_LADDER if c < budget) | {budget}))
    sys.setrecursionlimit(10000)
    print("orders %s, support %s" % (list(ORDERS), list(SUPPORT)))
    if single_rung:
        print("search at the single cap %d nodes per starting vertex, with no "
              "cheaper rungs, and the reductions tried in the paper's order: "
              "this is the paper's own protocol, and the cap is %s the paper's "
              "own budget of %d"
              % (caps[0],
                 "exactly" if caps[0] == PAPER_BUDGET_PER_START else "NOT",
                 PAPER_BUDGET_PER_START))
    else:
        print("search cap ladder %s nodes per starting vertex; the last rung is "
              "the paper's own budget of %d unless --budget= says otherwise "
              "(--paper-protocol searches at that single budget with no cheaper "
              "rungs, which is the paper's own protocol)"
              % (list(caps), PAPER_BUDGET_PER_START))
    total_all, _ = check_census()
    check_hypotheses()
    check_exhibited()
    check_image_table()
    check_automorphism_instance()
    check_corners()
    check_search_engine()
    check_ladder()
    attempted, settled, unsettled, compared = check_search(
        stride, caps, single_rung, total_all)
    if attempted == total_all:
        print("re-ran the search on all %d targets of the census, the whole of "
              "it (last rung %d nodes per start); the paper's own pass over the "
              "census spent 4.11e9 nodes" % (attempted, caps[-1]))
    else:
        print("re-ran the search on %d of the %d targets of the census "
              "(stride %d, last rung %d nodes per start); drop --stride= for "
              "all %d, the default, which attempts the hard targets too and "
              "therefore costs more: the paper's own pass over the census spent "
              "4.11e9 nodes"
              % (attempted, total_all, stride, caps[-1], total_all))
    if caps[-1] == PAPER_BUDGET_PER_START:
        rung_note = ("exactly the paper's own budget, so on each reduction "
                     "this search is allotted neither more nor less than the "
                     "paper allots itself")
    elif caps[-1] > PAPER_BUDGET_PER_START:
        rung_note = ("ABOVE the paper's own budget of %d, so a target settled "
                     "here need not be settled within the paper's allotment "
                     "and any ABOVE PAPER BUDGET line printed above names how "
                     "many such targets there were"
                     % PAPER_BUDGET_PER_START)
    else:
        rung_note = ("BELOW the paper's own budget of %d, so a target left "
                     "unsettled here need not be unsettled in the paper"
                     % PAPER_BUDGET_PER_START)
    items = []
    if attempted < total_all:
        items.append("the remaining %d targets of the census; dropping "
                     "--stride= attempts all %d at this same cap schedule, and "
                     "is the default, but this transcript covers only the %d "
                     "attempted here, so the settled count above is a lower "
                     "bound on the paper's %d and nothing here rules out a hard "
                     "target among those not attempted"
                     % (total_all - attempted, total_all, attempted, total_all))
    if compared:
        items.append("the paper's node counts (1.51e7 for its largest target, "
                     "4.11e9 in total), which are counts of the paper's own run: "
                     "this run searched at the paper's single budget of %d nodes "
                     "per start, in the paper's order of reductions, over the "
                     "whole census, and printed its own totals above, but a "
                     "difference between the two totals would show only that the "
                     "candidate ordering reimplemented here differs from the "
                     "paper's in some detail the paper does not state.  The "
                     "paper's 3792/43/288 split of targets by reduction is not "
                     "on this list: it is compared with this run's split on the "
                     "PAPER-PROTOCOL SPLIT line above" % caps[-1])
    else:
        items.append("the paper's node counts (1.51e7 for the largest target, "
                     "4.11e9 in total) and its 3792/43/288 split of targets by "
                     "reduction: the cap schedule in force here is %s nodes per "
                     "starting vertex, whose last rung is %s%s, and a target is "
                     "credited to whichever reduction settles it first at the "
                     "rung that settles it, so both the node counts and the "
                     "split printed above are this program's and not the "
                     "paper's; the SPLIT NOT COMPARABLE line above says which "
                     "of those obstacles applies to this run, and running with "
                     "--paper-protocol and no --stride= removes them all"
                     % (list(caps), rung_note,
                        ", the cheaper rungs being tried first"
                        if len(caps) > 1 else ""))
    items.append("the cited prior work, which covers the orders outside "
                 "{47,49,59} and which is what makes these three the only open "
                 "ones")
    items.append("the claim that (13,6,39) is one of ten triples that prior "
                 "work leaves open at v=59, a fact about that reference")
    if unsettled:
        items.append("certificates for the %d target(s) on the RESIDUAL lines "
                     "above, which this re-run did not produce" % len(unsettled))
    numerals = ("(i)", "(ii)", "(iii)", "(iv)", "(v)", "(vi)")
    print("NOT RE-RUN: " + "; ".join("%s %s" % (numerals[i], text)
                                     for i, text in enumerate(items)) + ".")
    n = len(CHECKS)
    k = sum(1 for _, o in CHECKS if not o)
    if k:
        # name the failing checks, so that the verdict is never a bare count
        print("FAILED CHECKS: " +
              ", ".join(name for name, o in CHECKS if not o))
        print("VERDICT: %d OF %d CHECKS FAILED" % (k, n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
